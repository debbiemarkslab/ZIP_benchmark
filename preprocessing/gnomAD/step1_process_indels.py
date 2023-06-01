import argparse 
from tqdm import tqdm 
import time 
import os 
import pandas as  pd 
import requests 
from collections import namedtuple 
import gzip 
import re 
from Bio import SeqIO
from io import StringIO
from multiprocessing import Pool 
import sys 
sys.path.append("../utils")
from gnomad_utils import add_mutation_type_annotations, add_frequency_annotations, drop_invalid_hgvsp_and_add_variants, add_depth_annotation

def drop_columns(
    df,
    to_keep=[
        "SYMBOL",
        "Gene",
        "HGVSp",
        "HGVSc",
        "AC",
        "AC_popmax",
        "AC_raw",
        "AF",
        "AF_popmax",
        "AF_raw",
        "protein_variant",
        "BIOTYPE",
        "UNIPARC"
    ]
):
    drop_columns = []
    for column in df.columns:
        if column in to_keep:
            continue 
        else:
            drop_columns.append(column)
    df = df.drop(drop_columns,axis=1)
    return df 

def get_uniprot_sequence_from_web(query_id):
    baseurl = "http://www.uniprot.org/uniparc/"+query_id+".fasta"
    response = requests.post(baseurl)
    cData=''.join(response.text)
    Seq = StringIO(cData)
    pSeq = list(SeqIO.parse(Seq, 'fasta'))
    return pSeq[0].seq 

# Constants and functions for parsing gnomAD file and VEP column 
# From PyVCF:
infos = {}
# Metadata parsers/constants
RESERVED_INFO = {
    'AA': 'String', 'AC': 'Integer', 'AF': 'Float', 'AN': 'Integer',
    'BQ': 'Float', 'CIGAR': 'String', 'DB': 'Flag', 'DP': 'Integer',
    'END': 'Integer', 'H2': 'Flag', 'H3': 'Flag', 'MQ': 'Float',
    'MQ0': 'Integer', 'NS': 'Integer', 'SB': 'String', 'SOMATIC': 'Flag',
    'VALIDATED': 'Flag', '1000G': 'Flag',

    # Keys used for structural variants
    'IMPRECISE': 'Flag', 'NOVEL': 'Flag', 'SVTYPE': 'String',
    'SVLEN': 'Integer', 'CIPOS': 'Integer', 'CIEND': 'Integer',
    'HOMLEN': 'Integer', 'HOMSEQ': 'String', 'BKPTID': 'String',
    'MEINFO': 'String', 'METRANS': 'String', 'DGVID': 'String',
    'DBVARID': 'String', 'DBRIPID': 'String', 'MATEID': 'String',
    'PARID': 'String', 'EVENT': 'String', 'CILEN': 'Integer',
    'DPADJ': 'Integer', 'CN': 'Integer', 'CNADJ': 'Integer',
    'CICN': 'Integer', 'CICNADJ': 'Integer'
}

# Conversion between value in file and Python value
field_counts = {
    '.': None,  # Unknown number of values
    'A': -1,  # Equal to the number of alternate alleles in a given record
    'G': -2,  # Equal to the number of genotypes in a given record
    'R': -3,  # Equal to the number of alleles including reference in a given record
}

def vcf_field_count(num_str):
        """Cast vcf header numbers to integer or None"""
        if num_str is None:
            return None
        elif num_str not in field_counts:
            # Fixed, specified number
            return int(num_str)
        else:
            return field_counts[num_str]
_Info = namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])

info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),\s*
            Number=(?P<number>-?\d+|\.|[AGR])?,\s*
            Type=(?P<type>Integer|Float|Flag|Character|String),\s*
            Description="(?P<desc>[^"]*)"
            (?:,\s*Source="(?P<source>[^"]*)")?
            (?:,\s*Version="?(?P<version>[^"]*)"?)?
            >''', re.VERBOSE)

def read_info(info_string):
        '''Read a meta-information INFO line.'''
        match = info_pattern.match(info_string)
        if not match:
            raise SyntaxError(
                "One of the INFO lines is malformed: %s" % info_string)

        num = vcf_field_count(match.group('number'))

        info = _Info(match.group('id'), num,
                     match.group('type'), match.group('desc'),
                     match.group('source'), match.group('version'))

        return (match.group('id'), info)

def _map(func, iterable, bad=['.', '']):
        '''``map``, but make bad values None.'''
        return [func(x) if x not in bad else None
                for x in iterable]
def _parse_info(info_str):
        '''Parse the INFO field of a VCF entry into a dictionary of Python
        types.

        '''
        if info_str == '.':
            return {}

        entries = info_str.split(';')
        retdict = {}
        # TODO: Loop over only certain columns 
        for entry in entries:
            entry = entry.split('=', 1)
            ID = entry[0]
            try:
                entry_type = infos[ID].type
            except KeyError:
                try:
                    entry_type = RESERVED_INFO[ID]
                except KeyError:
                    if entry[1:]:
                        entry_type = 'String'
                    else:
                        entry_type = 'Flag'

            if entry_type == 'Integer':
                vals = entry[1].split(',')
                try:
                    val = _map(int, vals)
                # Allow specified integers to be flexibly parsed as floats.
                # Handles cases with incorrectly specified header types.
                except ValueError:
                    val = _map(float, vals)
            elif entry_type == 'Float':
                vals = entry[1].split(',')
                val = _map(float, vals)
            elif entry_type == 'Flag':
                val = True
            elif entry_type in ('String', 'Character'):
                try:
                    vals = entry[1].split(',') # commas are reserved characters indicating multiple values
                    val = _map(str, vals)
                except IndexError:
                    entry_type = 'Flag'
                    val = True

            try:
                if infos[ID].num == 1 and entry_type not in ( 'Flag', ):
                    val = val[0]
            except KeyError:
                pass

            retdict[ID] = val

        return retdict

def get_info_dfs(df, multiprocess=False):
    if multiprocess:
        # Apply to all columns
        tqdm.pandas()
        n_cpus = len(os.sched_getaffinity(0))
        pool = Pool(n_cpus)
        info_dicts = tqdm(pool.imap(_parse_info, df.iloc[:len(df)]["INFO"].values))
        df_info = pd.DataFrame(list(info_dicts))
        return df_info 
    else:
        info_dicts = []
        for info_str in df["INFO"].values:
            info_dicts.append(_parse_info(info_str))
        return pd.DataFrame(info_dicts) 

def get_vep_dfs(df, df_info, infos):
    vep_string = infos['vep'].desc.split("Format: ")[-1]
    vep_fields = vep_string.split("|")
    vep_len = len(vep_fields)
    info_example = df.iloc[0]["INFO"]
    def split_and_check(s):
        s_split = s.split("|")
        assert len(s_split) == vep_len
        return s_split
    vep_split = df_info['vep'].explode()
    # Keep note of original index so that we can merge back later
    vep_split = vep_split.reset_index().rename(columns={'index': 'vep_index'})
    # Split the vep annotations into columns: https://stackoverflow.com/questions/35491274/split-a-pandas-column-of-lists-into-multiple-columns
    split_result = vep_split['vep'].apply(split_and_check)
    df_vep = pd.DataFrame(split_result.to_list(), index=vep_split["vep_index"], columns=vep_fields)
    df_vep = df_vep.replace("", pd.NA)
    return df_vep

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gnomad_database_file", help="Path to gnomAD database file. The file extension should be .vcf.bgz") 
    parser.add_argument("--gnomad_output_folder", help="Path to output folder for processed gnomad files.", default="../../processed_data/gnomAD")
    parser.add_argument("--exclude_insertions", action="store_true", help="Exclude insertions from the filtered files")
    parser.add_argument("--exclude_deletions", action="store_true", help="Exclude deletions from the filtered files")
    parser.add_argument("--exclude_duplications", action="store_true", help="Exclude duplications from the filtered files")
    parser.add_argument("--exclude_delins", action="store_true", help="Exclude delins from the filtered files")
    args = parser.parse_args()
    if not os.path.exists(args.gnomad_output_folder):
        os.mkdir(args.gnomad_output_folder)

    print(f"Filtering data from {args.gnomad_database_file} and converting it into a dataframe. File will be saved to {args.gnomad_output_folder}/gnomad_filtered_indels.csv.gz")
    start_time = time.time()
    # Find header:
    for line in gzip.open(args.gnomad_database_file, mode='rt'):
        if line.startswith("##"):
            continue
        # Break after finding header
        header = line.strip().split("\t")
        # print(f"Header={header}")
        break

    indicator_columns = ["inframe_dup","inframe_ins","inframe_del","inframe_delins","inframe_synon","inframe_single_sub", "terminating","unknown", "is_frameshift"]
    conflict_columns = ["inframe_dup","inframe_ins","inframe_del","inframe_delins","inframe_synon","inframe_single_sub"]
    infos = {}
    for line in gzip.open(args.gnomad_database_file, mode='rt'):
        # Stop after file header
        if not line.startswith("##"):
            break
        if line.startswith("##INFO"):
            # Get ID, Number, Type, Description
            # Note: Can also use PyVCF
            info_id, info = read_info(line.strip())
            infos[info_id] = info


    def get_filtered_dataframe(pool_args):
        i, df = pool_args
        df_info = get_info_dfs(df)
        df_vep = get_vep_dfs(df, df_info, infos)
        filtered_df_info = drop_columns(df_info)
        filtered_df_info[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']] = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']].reset_index(drop=True)
        filtered_df_vep = drop_columns(df_vep)
        df_merged_vep = filtered_df_info.merge(filtered_df_vep, left_index=True, right_on='vep_index')
        filtered_df = drop_invalid_hgvsp_and_add_variants(df_merged_vep)
        if "AC" in filtered_df.columns:
            filtered_df["AC"] = filtered_df["AC"].apply(lambda x: x[0] if type(x) != float else x)
            filtered_df["AC_raw"] = filtered_df["AC_raw"].apply(lambda x: x[0] if type(x) != float else x)
            filtered_df["AC_popmax"] = filtered_df["AC_popmax"].apply(lambda x: x[0] if type(x) != float else x)
        filtered_df["AF"] = filtered_df["AF"].apply(lambda x: x[0] if type(x) != float else x)
        filtered_df["AF_raw"] = filtered_df["AF_raw"].apply(lambda x: x[0] if type(x) != float else x)
        filtered_df["AF_popmax"] = filtered_df["AF_popmax"].apply(lambda x: x[0] if type(x) != float else x)
        # Dropping alleles with zero population max frequency (have a non-zero AF_raw value, but are dropped because of low-confidence). We could also switch and use AF_raw instead to include those lower confidence variants
        filtered_df = filtered_df[filtered_df["AF_popmax"].notnull() & (filtered_df["AF_popmax"] != 0.0)]
        if len(filtered_df) == 0:
            print(f"chunk {i} was empty after filtering")
            return None
        add_mutation_type_annotations(filtered_df)
        add_frequency_annotations(filtered_df)
        filtered_df["mutation_depth"] = filtered_df.apply(add_depth_annotation, axis=1)
        rare_df = filtered_df[filtered_df["frequency"] == "rare"]
        common_df = filtered_df[filtered_df["frequency"] == "common"]
        indiv_variant_dfs = filtered_df.groupby("vep_index", group_keys=False)
        rare_variant_dfs = rare_df.groupby("vep_index", group_keys=False)
        common_variant_dfs = common_df.groupby("vep_index", group_keys=False)
        # Getting genetic variant level classifications. Conflicts within protein annotations for one variant are labeled True in the 'conflicts' column. 
        compressed_df = indiv_variant_dfs[indicator_columns].any()
        rare_variant_compressed = rare_variant_dfs[indicator_columns].any()
        common_variant_compressed = common_variant_dfs[indicator_columns].any()
        compressed_df["mutation_depth"] = indiv_variant_dfs["mutation_depth"].unique()
        rare_variant_compressed["mutation_depth"] = rare_variant_dfs["mutation_depth"].unique()
        common_variant_compressed["mutation_depth"] = common_variant_dfs["mutation_depth"].unique()
        def check_conflict(x):
            if x.sum() > 1:
                return True 
            else:
                return False 
        # Checking for conflicts across types of mutations, removing conflicts, frameshift variants, terminating variants and variants with unknown effects
        compressed_df["conflicts"] = compressed_df[conflict_columns].apply(check_conflict, axis=1)
        rare_variant_compressed["conflicts"] = rare_variant_compressed[conflict_columns].apply(check_conflict, axis=1)
        common_variant_compressed["conflicts"] = common_variant_compressed[conflict_columns].apply(check_conflict, axis=1)
        filtered_compressed_df = compressed_df[(compressed_df["terminating"] == False) & (compressed_df["unknown"] == False) & (compressed_df["conflicts"] == False) & (compressed_df["is_frameshift"] == False)].copy()
        filtered_compressed_df["mutation_depth"] = filtered_compressed_df["mutation_depth"].apply(lambda x:x[0])
        full_filtered_df = filtered_df.loc[filtered_compressed_df.index] 
        # Applies filters to leave out certain types of mutations if arguments are specified 
        if args.exclude_insertions:
            full_filtered_df = full_filtered_df[full_filtered_df["inframe_ins"] == False]
        if args.exclude_deletions:
            full_filtered_df = full_filtered_df[full_filtered_df["inframe_del"] == False]
        if args.exclude_delins:
            full_filtered_df = full_filtered_df[full_filtered_df["inframe_delins"] == False]
        if args.exclude_duplications:
            full_filtered_df = full_filtered_df[full_filtered_df["duplication"] == False]
        return full_filtered_df

    chunk_dfs = pd.read_csv(args.gnomad_database_file, delimiter="\t", comment="#", compression='gzip', names=header, chunksize=5000)
    counter = 0 
    n_cpus = len(os.sched_getaffinity(0))
    pool = Pool(n_cpus)
    stats = pool.imap(get_filtered_dataframe, enumerate(chunk_dfs), n_cpus*2)
    # print(stats)


    total_records = 0

    df = next(stats)
    while df is None:
        df = next(stats)
    df.to_csv(f"{args.gnomad_output_folder}/gnomad_filtered_indels.csv.gz",compression='gzip', index=False, header=True)
    total_records += len(df)
    for df in tqdm(stats):
        if df is None:
            continue
        df.to_csv(f"{args.gnomad_output_folder}/gnomad_filtered_indels.csv.gz", mode='a', compression='gzip', index=False, header=False)
        total_records += len(df)

    print(f"Total entries that passed filtering: {total_records}")
    print(f"Filtering and Processing completed in {time.time() - start_time} seconds, gnomAD filtered indels saved to {args.gnomad_output_folder}/gnomad_filtered_indels.csv.gz")
