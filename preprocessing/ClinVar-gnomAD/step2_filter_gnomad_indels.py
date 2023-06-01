import argparse 
import os 
from glob import glob 
import pandas as pd
import matplotlib.pyplot as plt 
from Bio.SeqUtils import IUPACData
from tqdm.auto import tqdm 
from collections import namedtuple
from requests_futures.sessions import FuturesSession
from multiprocessing import Pool
import time 
import sys 
sys.path.append("../utils")
from mutated_sequence_utils import hgvs_regex, refseq_regex, insert, delete, delins, duplicate 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--filtered_gnomad_file", help="Path to filtered gnomAD file output from process_indels_gnomad_data.py. The file extension should be .csv.gz", default=None) 
    parser.add_argument("--gnomad_output_folder", help="Path to output folder for processed gnomad files. Default is ./gnomad_output", default="./gnomad_output")
    parser.add_argument("--allele_frequency_threshold", help="Threshold for determining common or rare variants. Variants with AF values above threshold will be considered common. Default is 0.05", default=0.05, type=float)
    args = parser.parse_args()
    if args.filtered_gnomad_file == None:
        if not os.path.exists(f"{args.gnomad_output_folder}/gnomad_filtered_indels.csv.gz"):
            raise ValueError(f"Please specify a valid path to a filtered gnomAD file. The file extension should be .csv.gz. This file is ouput from step1_process_indels.py.")
        else:
            filtered_gnomad_file = f"{args.gnomad_output_folder}/gnomad_filtered_indels.csv.gz"
    else:
        filtered_gnomad_file = args.filtered_gnomad_file

    start_time = time.time()
    print(f"Now saving file of common gnomAD variants with allele frequency > {args.allele_frequency_threshold} to {args.gnomad_output_folder}/gnomad_filtered_common_indels.csv")

    dfs = []

    for df in tqdm(pd.read_csv(
        f"{filtered_gnomad_file}",
        compression="gzip",
        chunksize=10000,
    )):
        dfs.append(df.loc[(df["FILTER"] == "PASS") & (df['inframe_delins'] | df['inframe_ins'] | df['inframe_del'] | df['inframe_dup']) & (df['AF_popmax'] > args.allele_frequency_threshold)])
    df = pd.concat(dfs)
    all_common_indels_df = df
    df.to_csv(f"{args.gnomad_output_folder}/gnomad_filtered_common_indels.csv", index=False)

    print(f"Common gnomad variants saved to {args.gnomad_output_folder}/gnomad_filtered_common_indels.csv. Now saving deduped and aa_filtered files to {args.gnomad_output_folder}/gnomad_filtered_common_indels_dedup.csv and {args.gnomad_output_folder}/gnomad_filtered_common_indels_dedup_lt3aa.csv")

    # For each uniparc id, we fetch the canonical sequence, or the longest sequence if the canonical sequence can't be found  
    def response_hook(resp, *args, **kwargs):
        # parse the json storing the result on the response object
        resp.data = resp.json()

    session = FuturesSession()
    session.hooks['response'] = response_hook
    uniparc_info = {} 
    futures = {}
    for uniparc_id in tqdm(all_common_indels_df["UNIPARC"]):
        if uniparc_id in uniparc_info:
            continue
        futures[uniparc_id] = session.get(f"https://rest.uniprot.org/uniparc/{uniparc_id}")

    for uniparc_id, future in tqdm(futures.items()):
        response = future.result()
        if response.status_code != 200:
            print('response status {0}'.format(response.status_code))
            continue
        uniparc_info[uniparc_id] = response.data

    # Function to take in a group of variants and determine the best uniparc id to use for each variant 
    def get_uniparc_seq(uniparc_id):
        return uniparc_info[uniparc_id]["sequence"]["value"]

    def choose_best_uniparc_group(group):
        best_uniparc = choose_best_uniparc(group["UNIPARC"])
        return group[group["UNIPARC"] == best_uniparc]

    def choose_best_uniparc(uniparc_ids):
        if len(uniparc_ids) == 1:
            return uniparc_ids.iloc[0]
        # choose by the reviewed UniProt ID activity, then the longest sequence
        uniprot_ids = list(
            (
                (uniparc_id, cr["id"]), 
                (cr["active"], uniparc_info[uniparc_id]["sequence"]["length"])
            ) for uniparc_id in uniparc_ids for cr in uniparc_info[uniparc_id]["uniParcCrossReferences"] if cr["database"] == "UniProtKB/Swiss-Prot")
        if len(uniprot_ids) > 0:
            uniprot_ids.sort(key=lambda x: x[1], reverse=True)
            if len(uniprot_ids) > 1:
                print("IDs:\t", uniprot_ids, sep='\t')
            return uniprot_ids[0][0][0]

        # choose by the reviewed UniProt ID isoform activity, then the longest sequence
        uniprot_ids = list(
            (
                (uniparc_id, cr["id"]),
                (cr["active"], uniparc_info[uniparc_id]["sequence"]["length"])
            ) for uniparc_id in uniparc_ids for cr in uniparc_info[uniparc_id]["uniParcCrossReferences"] if cr["database"] == "UniProtKB/Swiss-Prot protein isoforms")
        if len(uniprot_ids) > 0:
            uniprot_ids.sort(key=lambda x: x[1], reverse=True)
            if len(uniprot_ids) > 1:
                print("Isoforms:", uniprot_ids, sep='\t')
            return uniprot_ids[0][0][0]

        # choose by X in sequence, named gene, named protein, unreviewed UniProt ID (TrEMBL) activity, non-fragment, non-cDNA, the shortest gene name, the then the longest sequence
        uniprot_ids = list(
            (
                i, 
                ("X" not in uniparc_info[uniparc_id]["sequence"]["value"], cr["active"], "geneName" in cr, "proteinName" in cr, "Fragment" not in cr["proteinName"] if "proteinName" in cr else False, "cDNA" not in cr["proteinName"] if "proteinName" in cr else False, -len(cr["geneName"]) if "geneName" in cr else -99, uniparc_info[uniparc_id]["sequence"]["length"])
            ) for i, uniparc_id in enumerate(uniparc_ids) for cr in uniparc_info[uniparc_id]["uniParcCrossReferences"] if cr["database"] == "UniProtKB/TrEMBL")
        if len(uniprot_ids) > 0:
            uniprot_ids.sort(key=lambda x: x[1], reverse=True)
            return uniparc_ids.iloc[uniprot_ids[0][0]]

    chosen_uniparc_df = all_common_indels_df.drop_duplicates(subset=["UNIPARC", "protein_variant"]).groupby(["#CHROM", "POS", "REF", "ALT"]).apply(choose_best_uniparc_group)
    chosen_uniparc_df.index = chosen_uniparc_df.index.droplevel([0, 1, 2, 3])
    chosen_uniparc_df["protein_sequence"] = chosen_uniparc_df["UNIPARC"].apply(get_uniparc_seq)

    chosen_uniparc_df['protein_variant_1letter'] = chosen_uniparc_df['protein_variant']
    for aa3, aa1 in tqdm(IUPACData.protein_letters_3to1.items()):
        chosen_uniparc_df['protein_variant_1letter'] = chosen_uniparc_df['protein_variant_1letter'].str.replace(aa3, aa1, regex=False)
    assert not chosen_uniparc_df['protein_variant_1letter'].isna().any()

    # Check later that all sequences have been mutated
    chosen_uniparc_df['protein_sequence_mutated'] = chosen_uniparc_df['protein_sequence'] 

    # Insertions
    # insert() will check all tagged insertions in correct format
    chosen_uniparc_df.loc[chosen_uniparc_df['inframe_ins'], 'protein_sequence_mutated'] = \
        chosen_uniparc_df[chosen_uniparc_df['inframe_ins']].apply(lambda row: insert(row['protein_variant_1letter'], SEQ=row['protein_sequence']), axis=1)
    # Deletions on their own
    chosen_uniparc_df.loc[chosen_uniparc_df['inframe_del'], 'protein_sequence_mutated'] = \
    chosen_uniparc_df[chosen_uniparc_df['inframe_del']].apply(lambda row: delete(row['protein_variant_1letter'], SEQ=row['protein_sequence']), axis=1)
    # Duplications on their own
    chosen_uniparc_df.loc[chosen_uniparc_df['inframe_dup'], 'protein_sequence_mutated'] = \
    chosen_uniparc_df[chosen_uniparc_df['inframe_dup']].apply(lambda row: duplicate(row['protein_variant_1letter'], SEQ=row['protein_sequence']), axis=1)
    # Sanity check on remaining rows 
    rows_not_mutated = chosen_uniparc_df[chosen_uniparc_df['protein_sequence'] == chosen_uniparc_df['protein_sequence_mutated']]
    assert rows_not_mutated['inframe_delins'].all()
    rows_not_mutated['protein_variant_1letter']
    # Delins on their own (copied to clinvar_utils)
    # Check delins() function in clinvar_utils to see how it works
    chosen_uniparc_df.loc[chosen_uniparc_df['inframe_delins'], 'protein_sequence_mutated'] = \
        chosen_uniparc_df[chosen_uniparc_df['inframe_delins']].apply(lambda row: delins(row['protein_variant_1letter'], SEQ=row['protein_sequence']), axis=1)

    # We've mutated every sequence
    rows_not_mutated = chosen_uniparc_df[chosen_uniparc_df['protein_sequence'] == chosen_uniparc_df['protein_sequence_mutated']]
    assert len(rows_not_mutated) == 0
    # Save the deduped and deduped less than three amino acid variants to csvs 
    chosen_uniparc_df.to_csv(f"{args.gnomad_output_folder}/gnomad_filtered_common_indels_dedup.csv",index=False)
    chosen_uniparc_df[chosen_uniparc_df["mutation_depth"] <= 3].to_csv(f"{args.gnomad_output_folder}/gnomad_filtered_common_indels_dedup_lt3aa.csv", index=False)
    print(f"Deduplication complete. Total time taken is {time.time() - start_time}")