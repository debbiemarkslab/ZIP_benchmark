import argparse 
import pandas as pd 
from Bio.SeqUtils import IUPACData
from Bio import SeqIO
from tqdm import tqdm 
from collections import defaultdict
import sys 
import gzip 
import os 
import pickle 
sys.path.append("..")
from utils.mutated_sequence_utils import hgvs_regex, refseq_regex, insert, delete, delins, duplicate

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="process clinvar data")
    parser.add_argument("--input_file", type=str, help="Path to input clinvar query file")
    parser.add_argument("--refseq_sequence_file", type=str, help="Path to refseq sequence file")
    parser.add_argument("--grch38_protein_file", type=str, help="Path to GRCh38 protein fasta file")
    parser.add_argument("--max_context_len", type=int, default=1022, help="Adds a column flagging sequences with length greater than the given value. Useful for marking sequences below the maximum context window of a language model")
    parser.add_argument("--max_len", type=int, default=6000, help="Sequences greater than this length are filtered out.")
    parser.add_argument("--max_depth", type=int, default=3, help="Mutations with depth greater than this value are filtered out")
    parser.add_argument("--save_benign", action="store_true", help="Save benign variants as well. Defalt is just pathogenic")
    parser.add_argument("--uniprot_refseq_mapping_file", type=str, help="Path to uniprot mapping file")
    parser.add_argument("--output_folder", type=str, help="Path to output folder")
    args = parser.parse_args()
    
    pathogenic_categories = ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]
    benign_categories = ["Benign", "Likely benign", "Benign/Likely benign"]

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    df_clinvar_query = pd.read_table(args.input_file)
    # Removing extra empty column 
    df_clinvar_query = df_clinvar_query.drop(columns=["Unnamed: 15"]).copy()

    review_stars_mapping = {
    "practice guideline": 4,
    "reviewed by expert panel": 3, 
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion provided": 0,
    "no interpretation for the single variant": 0,
    "no assertion criteria provided": 0,
    }

    df_clinvar_query["Stars"] = df_clinvar_query["Review status"].apply(review_stars_mapping.get)
    split_significance_regex = r"^(.+)\((.+)\)$"
    df_clinvar_query[['Clinical significance', 'Last reviewed']] = df_clinvar_query["Clinical significance (Last reviewed)"].str.extract(split_significance_regex)
    df_clinvar_query['hgvs_valid'] = df_clinvar_query['Name'].str.contains(hgvs_regex)
    df_clinvar_query['protein_variant'] = df_clinvar_query['Name'].str.extract(hgvs_regex)

    # Adding flag columns for kinds of mutations 
    filter_cols = []
    synon_stop_regex = "=|\*|Ter"
    df_clinvar_query['inframe_dup'] = df_clinvar_query['protein_variant'].str.contains("dup")
    df_clinvar_query['inframe_synon_stop'] = df_clinvar_query['protein_variant'].str.contains(synon_stop_regex)
    filter_cols.append('inframe_synon_stop')
    df_clinvar_query['question_mark'] = df_clinvar_query['protein_variant'].str.contains("\?")
    filter_cols.append('question_mark')
    df_clinvar_query['square_bracket'] = df_clinvar_query['protein_variant'].str.contains("\[\d+\]")
    filter_cols.append('square_bracket')
    df_clinvar_query['hgvs_frameshift'] = df_clinvar_query['protein_variant'].str.contains("fs")
    filter_cols.append('hgvs_frameshift')

    df_clinvar_query['invalid_aa'] = df_clinvar_query['protein_variant'].str.contains("X")
    filter_cols.append('invalid_aa')

    # Filtering out variants we don't want to consider 
    df_clinvar_clean = df_clinvar_query[df_clinvar_query['hgvs_valid']].copy()
    for col in filter_cols:
        df_clinvar_clean = df_clinvar_clean[~df_clinvar_clean[col].astype(bool)]

    # Adding more flags for insertion, deletion and delins mutations 
    df_clinvar_clean['inframe_delins'] = df_clinvar_clean['protein_variant'].str.contains("delins")
    df_clinvar_clean['inframe_ins'] = df_clinvar_clean['protein_variant'].str.contains("ins") & ~df_clinvar_clean['inframe_delins']
    df_clinvar_clean['inframe_del'] = df_clinvar_clean['protein_variant'].str.contains("del") & ~df_clinvar_clean['inframe_delins']

    # Cast all to bool
    df_clinvar_clean['inframe_delins'] = df_clinvar_clean['inframe_delins'].astype(bool)
    df_clinvar_clean['inframe_ins'] = df_clinvar_clean['inframe_ins'].astype(bool)
    df_clinvar_clean['inframe_del'] = df_clinvar_clean['inframe_del'].astype(bool)
    df_clinvar_clean['inframe_dup'] = df_clinvar_clean['inframe_dup'].astype(bool)

    # Getting wt positions and amino acids 
    pos_regex = "([A-Z][a-z]{2})(\d+)"
    df_clinvar_clean[['aa_wt_start_3letter', 'pos_start']] = df_clinvar_clean['protein_variant'].str.extract(pos_regex)
    df_clinvar_clean['pos_start'] = df_clinvar_clean['pos_start'].astype('int')
    df_clinvar_clean['aa_wt_start'] = df_clinvar_clean['aa_wt_start_3letter'].apply(lambda s: IUPACData.protein_letters_3to1[s])  # Using [s] instead of .get since .get silently returns NaNs
    df_clinvar_clean['protein_variant_1letter'] = df_clinvar_clean['protein_variant']
    for aa3, aa1 in IUPACData.protein_letters_3to1.items():
        df_clinvar_clean['protein_variant_1letter'] = df_clinvar_clean['protein_variant_1letter'].str.replace(aa3, aa1, regex=False)
    
    # Adding in refseq ids and sequences 
    matches = df_clinvar_clean["Name"].str.contains(refseq_regex, regex=True)
    df_clinvar_clean['refseq_id'] = df_clinvar_clean["Name"].str.extract(refseq_regex)
    df_refseq_mapping = pd.read_csv(args.refseq_sequence_file, delimiter='\t', low_memory=False)
    column_mapper = {'protein': 'refseq_protein_id', 'mrna': 'refseq_mrna_id'} # Can also rename 'gene'
    df_merged = df_clinvar_clean.merge(df_refseq_mapping.rename(columns=column_mapper), how='left', left_on='refseq_id', right_on='refseq_mrna_id')
    df_merged['has_mrna_mapping'] = ~df_merged['refseq_mrna_id'].isna()
    df_clinvar_clean_mapped = df_merged[df_merged['has_mrna_mapping']].copy()
    fasta_38_dict = SeqIO.to_dict(SeqIO.parse(args.grch38_protein_file, "fasta"))
    df_clinvar_clean_mapped['protein_sequence'] = df_clinvar_clean_mapped['refseq_protein_id'].apply(lambda s: str(fasta_38_dict[s].seq))
    
    # Adding in uniprot IDs 
    # There may be multiple uniprot ids (values) for each refseq id (key), which is why we store them as a set and .add() each one
    refseq_uniprot_dict = defaultdict(set)
    with gzip.open(args.uniprot_refseq_mapping_file, mode='rt') as f:
        for i, line in enumerate(tqdm(f)):
            refseq, uniprot_id = line.split("\t")
            if refseq.startswith("NP_"):
                refseq_uniprot_dict[refseq].add(uniprot_id.strip())  
    df_clinvar_clean_mapped['uniprot_id'] = df_clinvar_clean_mapped['refseq_protein_id'].apply(refseq_uniprot_dict.get)

    # Adding mutated sequences 
    df_clinvar_clean_mapped['protein_sequence_mutated'] = df_clinvar_clean_mapped['protein_sequence'] 
    # Insertions
    df_clinvar_clean_mapped.loc[df_clinvar_clean_mapped['inframe_ins'], 'protein_sequence_mutated'] = \
    df_clinvar_clean_mapped[df_clinvar_clean_mapped['inframe_ins']].apply(lambda row: insert(row['protein_variant_1letter'], SEQ=row['protein_sequence']), axis=1)
    # Deletions
    df_clinvar_clean_mapped.loc[df_clinvar_clean_mapped['inframe_del'], 'protein_sequence_mutated'] = \
    df_clinvar_clean_mapped[df_clinvar_clean_mapped['inframe_del']].apply(lambda row: delete(row['protein_variant_1letter'], SEQ=row['protein_sequence']), axis=1)
    # Duplications
    df_clinvar_clean_mapped.loc[df_clinvar_clean_mapped['inframe_dup'], 'protein_sequence_mutated'] = \
    df_clinvar_clean_mapped[df_clinvar_clean_mapped['inframe_dup']].apply(lambda row: duplicate(row['protein_variant_1letter'], SEQ=row['protein_sequence']), axis=1)
    # Delins
    df_clinvar_clean_mapped.loc[df_clinvar_clean_mapped['inframe_delins'], 'protein_sequence_mutated'] = \
    df_clinvar_clean_mapped[df_clinvar_clean_mapped['inframe_delins']].apply(lambda row: delins(row['protein_variant_1letter'], SEQ=row['protein_sequence']), axis=1)
    
    # Adding mutation depth
    df_clinvar_clean_mapped['mutation_depth'] = df_clinvar_clean_mapped['protein_sequence_mutated'].str.len() - df_clinvar_clean_mapped['protein_sequence'].str.len()

    df_clinvar_clean_mapped[f'filter_lt_{args.max_depth}_aas'] = (df_clinvar_clean_mapped['mutation_depth'].abs() <= args.max_depth)
    df_clinvar_clean_mapped[f'filter_lt_{args.max_context_len}_seqlen'] = (df_clinvar_clean_mapped['protein_sequence_mutated'].str.len() <= args.max_context_len) & (df_clinvar_clean_mapped['protein_sequence'].str.len() <= args.max_context_len)
    df_clinvar_clean_mapped[f'filter_lt_{args.max_len}_seqlen'] = (df_clinvar_clean_mapped['protein_sequence_mutated'].str.len() <= args.max_len) & (df_clinvar_clean_mapped['protein_sequence'].str.len() <= args.max_len)
    df_subset_scoring = df_clinvar_clean_mapped[df_clinvar_clean_mapped[f'filter_lt_{args.max_len}_seqlen'] & df_clinvar_clean_mapped['filter_lt_3_aas']]

    # Filtering out genes with no pathogenic variants 
    def filter_significance(group):
        return (group['Clinical significance'].isin(pathogenic_categories).sum() >= 1)
    rows_filtered_genes = df_subset_scoring.groupby("Gene(s)").filter(filter_significance) 
    significant_rows = rows_filtered_genes[rows_filtered_genes["Clinical significance"].isin(pathogenic_categories+benign_categories)].copy()
    mapper = {"Likely benign": "Benign", "Likely pathogenic": "Pathogenic", "Benign/Likely benign": "Benign", "Pathogenic/Likely pathogenic": "Pathogenic"}
    significant_rows['ClinSigSimple'] = significant_rows['Clinical significance'].apply(lambda s: mapper.get(s, s))

    significant_rows["label"] = significant_rows["ClinSigSimple"] 
    significant_rows["HGVSp"] = significant_rows["protein_variant"]
    if args.save_benign:
        significant_rows.to_csv(args.output_folder + os.sep + "clinvar_filtered_significant_variants.csv", index=False)
    else:
        significant_rows = significant_rows[significant_rows["ClinSigSimple"] != "Benign"]
        significant_rows.to_csv(args.output_folder + os.sep + "clinvar_filtered_significant_pathogenic_variants.csv", index=False)
    # rows_filtered_genes.to_csv(args.output_folder + os.sep + "clinvar_filtered_variants.csv", index=False)
    # df_clinvar_clean_mapped.to_csv(args.output_folder + os.sep + "clinvar_unfiltered_variants.csv", index=False)
