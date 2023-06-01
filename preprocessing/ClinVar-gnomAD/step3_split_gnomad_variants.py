import pandas as pd 
import os 
import argparse 
from tqdm import tqdm
import sys
sys.path.append("../utils")
from mapfile_utils import create_mapfile
f"""
This script takes one of the gnomAD csvs processed by step1_process_indels.py and step2_filter_indels.py and splits it into multiple csvs, one for each uniparc id.
If you have alignments for each uniparc sequence (to use e.g. with Tranception), you can pass the path to the folder containing the alignments with the --alignment_folder argument.
The alignment folder must have the format <folder>/<uniparc_id>.fa for each uniparc id.
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split gnomAD data by uniparc id")
    parser.add_argument("-input_file", type=str, help="Path to processed gnomAD csv file")
    parser.add_argument("-output_folder", type=str, help="Path to output folder where individual uniparc csv files will be saved")
    parser.add_argument("--alignment_folder", type=str, help="Path to folder containing alignments for each uniparc id sequence")
    parser.add_argument("-mapping_file", type=str, help="Path to save the mapping file to")
    args = parser.parse_args()
    if not os.path.isdir(args.output_folder):
        os.mkdir(args.output_folder)

    df = pd.read_csv(args.input_file)
    uniq_col = "UNIPARC"
    print("Splitting data by uniparc id and saving to individual csv files")
    for uniparc_id in tqdm(df[uniq_col].unique()):
        subset = df[df[uniq_col] == uniparc_id]
        protein_sequences = subset["protein_sequence"].unique()
        # Check that there is a single unique target sequence for this uniparc id
        assert len(protein_sequences) == 1 
        # Adding labels here (all gnomAD variants common variants are assumed benign)
        subset["label"] = "Benign"
        subset.to_csv(os.path.join(args.output_folder, f"{uniparc_id}.csv"), index=False) 