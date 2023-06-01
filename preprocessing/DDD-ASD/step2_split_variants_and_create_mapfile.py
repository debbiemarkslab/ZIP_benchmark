import argparse 
import pandas as pd 
import sys 
import os 
sys.path.append("..")
from utils.mapfile_utils import create_mapfile

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--indel_variant_file", type=str, help="Path to csv file with indel variants", default="../../processed_data/DDD/ddd_inframe_indels_lt3aa.tsv")
    parser.add_argument("--output_folder", type=str, help="Folder to write individual variant files to", default="../../processed_data/DDD/by_protein")
    parser.add_argument("--mapfile", type=str, help="Path to write mapfile to", default="../../processed_data/DDD/mapfiles/ddd_inframe_indels_lt3aa_mapfile.csv")
    parser.add_argument("--alignment_folder", default=None, type=str, help="Path to folder containing alignments for each sequence")
    args = parser.parse_args()
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    df_inframe_mutated = pd.read_csv(args.indel_variant_file, sep="\t") 
    df_unique_protein = df_inframe_mutated.drop_duplicates(subset=["protein_sequence"])
    df_inframe_mutated["mutation_hgvs"] = df_inframe_mutated["protein_variant"]
    df_inframe_dedup_id = df_inframe_mutated.drop_duplicates(subset=["protein_variant","protein_sequence","protein_sequence_mutated","mutation_hgvs"])
    for protein in df_inframe_dedup_id["protein"].unique():
        subset = df_inframe_dedup_id[df_inframe_dedup_id["protein"] == protein]
        subset.to_csv(f"{args.output_folder}/{protein}.csv")
    create_mapfile(args.output_folder, args.mapfile, args.alignment_folder)    
