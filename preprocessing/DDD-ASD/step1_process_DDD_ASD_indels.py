import argparse 
import pandas as pd
import os 
from glob import glob 
from tqdm import tqdm 
from pyfaidx import Fasta 
from Bio import Entrez, SeqIO 
from Bio.SeqUtils import IUPACData
import sys 
sys.path.append("..")
from utils.mutated_sequence_utils import delete, insert, duplicate, delins

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description='Process DDD indels') 
    parser.add_argument("--input_file", type=str, help="file with combined DDD/ASD variants")
    parser.add_argument("--output_folder", type=str, help="Folder where processed data will be saved", default="../../processed_data/DDD")
    parser.add_argument("--refseq_reference_folder", type=str, help="mapping file for refseq variants")
    parser.add_argument("--entrez_email", type=str, help="email to use for entrez queries")
    parser.add_argument("--canonical_isoform_file", type=str, help="file containing canonical isoforms for refseq proteins")
    parser.add_argument("--length_threshold", type=int, default=6000)
    df = pd.read_csv(parser.parse_args().input_file, sep="\t") 
    args = parser.parse_args()
    df = pd.read_csv(args.input_file, sep="\t", low_memory=False)
    df = df.dropna(subset="Consequence")
    df.replace("-", pd.NA, inplace=True)
    df_inframe = df[df["Consequence"].str.contains("inframe")]
    print(df_inframe)
    exit()
    records = dict()
    proteins_query = set(df_inframe["protein"].unique())

    for filepath in tqdm(glob(f'{args.refseq_reference_folder}/*.faa')):
        fasta = Fasta(filepath)
        headers = fasta.keys()
        for protein in set(headers) & proteins_query:
            records[protein] = str(fasta[protein])
    df_inframe_mapped = df_inframe.copy()
    df_inframe_mapped['protein_sequence'] = df_inframe_mapped['protein'].apply(lambda s: str(records[s]) if s in records else pd.NA)
    missing = df_inframe_mapped[df_inframe_mapped['protein_sequence'].isna()]["protein"].unique()
    # Write out missing proteins to file
    Entrez.email = args.entrez_email
    with open("./temp_missing_seqs.faa", "w") as f:
        handle = Entrez.efetch(db="protein", id=missing, rettype="fasta", retmode="fasta")
        f.write(handle.read())
        handle.close()
    for filepath in glob('./temp_missing_seqs.faa'):
        fasta = Fasta(filepath)
        headers = fasta.keys()
        for protein in set(headers) & proteins_query:
            records[protein] = str(fasta[protein])
    os.remove("./temp_missing_seqs.faa") 
    os.remove("./temp_missing_seqs.faa.fai")
    canonical_list = args.canonical_isoform_file
    with open(canonical_list, "r") as f:
        canonical_refseq_proteins = set([line.strip() for line in f.readlines()])
    df_inframe_mapped["is_RefSeq_Select_canonical"] = df_inframe_mapped["protein"].isin(canonical_refseq_proteins)
    df_inframe_mapped_dedup = df_inframe_mapped.drop_duplicates(subset=["id", "pos", "ref", "alt", "mutant", "protein"])  # Note: Which one to keep?
    # Per-gene filter out non-canonical isoforms
    # For each gene, for each person exclude duplicates by keeping the canonical isoform.
    # If there are no canonical isoforms, keep all of them
    def filter_indiv(indiv):
        if indiv["is_RefSeq_Select_canonical"].sum() > 1:
            print(f"More than one canonical isoform per person per gene")
            assert False
        elif indiv["is_RefSeq_Select_canonical"].sum() == 0:
            print(f"No isoforms for gene {indiv['gene'].iloc[0]}, person id {indiv['id'].iloc[0]}, keeping all {len(indiv)} variants")
            return indiv
        return indiv[indiv["is_RefSeq_Select_canonical"]] if indiv["is_RefSeq_Select_canonical"].sum() == 1 else indiv

    filter_groups = []
    for group in df_inframe_mapped_dedup.groupby("gene"):
        id_group = group[1].groupby("id")
        for group_id, group_content in id_group:
            filter_group = filter_indiv(group_content)
            filter_groups.append(filter_group)
    df_inframe_mapped_isoforms = pd.concat(filter_groups)

    # Removing variants with NaN HGVS annotations 
    diffs = df_inframe_mapped_isoforms["HGVSp"].str.split(":p.").str[0] != df_inframe_mapped_isoforms["protein"]
    print(f"Number of differences between HGVSp and protein: {diffs.sum()}")
    df_inframe_mapped_isoforms_valid = df_inframe_mapped_isoforms[~diffs].copy()

    # Dropping Terminating mutations 
    df_inframe_mapped_clean = df_inframe_mapped_isoforms_valid[~df_inframe_mapped_isoforms_valid["HGVSp"].str.contains("Ter")].copy()

    # Adding protein variant column 
    df_inframe_mapped_clean["protein_variant"] = df_inframe_mapped_clean["HGVSp"].str.split(":p.", regex=False).str[1]
    df_inframe_mapped_clean["protein_variant_1letter"] = df_inframe_mapped_clean["protein_variant"]
    for aa3, aa1 in IUPACData.protein_letters_3to1.items():
        df_inframe_mapped_clean["protein_variant_1letter"] = df_inframe_mapped_clean["protein_variant_1letter"].str.replace(aa3, aa1)
        
    # Adding mutation type annotations 
    df_inframe_mapped_clean["inframe_delins"] = df_inframe_mapped_clean["HGVSp"].str.contains("delins")
    df_inframe_mapped_clean["inframe_deletion"] = df_inframe_mapped_clean["HGVSp"].str.contains("del") & ~df_inframe_mapped_clean["inframe_delins"]
    df_inframe_mapped_clean["inframe_insertion"] = df_inframe_mapped_clean["HGVSp"].str.contains("ins") & ~df_inframe_mapped_clean["inframe_delins"]
    df_inframe_mapped_clean["inframe_dup"] = df_inframe_mapped_clean["HGVSp"].str.contains("dup")

    # Finding mismatches between variant annotations and sequences 
    df_inframe_mapped_clean[["aa_start", "start_pos"]] = df_inframe_mapped_clean["protein_variant_1letter"].str.extract(r"([A-Z])(\d+).*")
    matching = df_inframe_mapped_clean.apply(lambda row: row["aa_start"] == row["protein_sequence"][int(row["start_pos"])-1], axis=1)
    df_inframe_mapped_clean_mutated = df_inframe_mapped_clean[matching].copy()
    df_inframe_mapped_clean_mutated["protein_sequence_mutated"] = pd.NA

    # Deletions
    df_inframe_mapped_clean_mutated.loc[df_inframe_mapped_clean_mutated["inframe_deletion"], "protein_sequence_mutated"] = \
        df_inframe_mapped_clean_mutated.loc[df_inframe_mapped_clean_mutated["inframe_deletion"]].apply(lambda row: delete(SEQ=row["protein_sequence"], mut=row["protein_variant_1letter"]), axis=1)
    # Insertions
    df_inframe_mapped_clean_mutated.loc[df_inframe_mapped_clean_mutated["inframe_insertion"], "protein_sequence_mutated"] = \
        df_inframe_mapped_clean_mutated.loc[df_inframe_mapped_clean_mutated["inframe_insertion"]].apply(lambda row: insert(SEQ=row["protein_sequence"], mut=row["protein_variant_1letter"]), axis=1)
    # Duplication
    df_inframe_mapped_clean_mutated.loc[df_inframe_mapped_clean_mutated["inframe_dup"], "protein_sequence_mutated"] = \
        df_inframe_mapped_clean_mutated.loc[df_inframe_mapped_clean_mutated["inframe_dup"]].apply(lambda row: duplicate(SEQ=row["protein_sequence"], mut=row["protein_variant_1letter"]), axis=1)
    # Delins
    df_inframe_mapped_clean_mutated.loc[df_inframe_mapped_clean_mutated["inframe_delins"], "protein_sequence_mutated"] = \
        df_inframe_mapped_clean_mutated.loc[df_inframe_mapped_clean_mutated["inframe_delins"]].apply(lambda row: delins(SEQ=row["protein_sequence"], mut=row["protein_variant_1letter"]), axis=1)
    df_inframe_mapped_clean_mutated.to_csv(f"{args.output_folder}/ddd_inframe_indels.tsv", sep="\t", index=False)
    df_inframe_mapped_clean_mutated["mutation_depth"] = df_inframe_mapped_clean_mutated["protein_sequence"].str.len() - df_inframe_mapped_clean_mutated["protein_sequence_mutated"].str.len()
    df_inframe_lt3aa = df_inframe_mapped_clean_mutated[df_inframe_mapped_clean_mutated["mutation_depth"].abs() < 3]
    # filtering out very long proteins 
    df_inframe_lt3aa = df_inframe_lt3aa[df_inframe_lt3aa["protein_sequence"].str.len() < args.length_threshold]
    df_inframe_lt3aa.to_csv(f"{args.output_folder}/ddd_inframe_indels_lt3aa.tsv", sep="\t", index=False)
