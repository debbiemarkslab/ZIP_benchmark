import pandas as pd 
import os 
import argparse 
from tqdm import tqdm
import sys
sys.path.append("../utils")
from mapfile_utils import create_mapfile
f"""
This script creates the mapping file for the individual mutational scans. 
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make mapfile for deep mutational scans")
    parser.add_argument("-input_folder", type=str, help="Path to folder containing mutational scan csvs")
    parser.add_argument("--alignment_folder", type=str, help="Path to folder containing alignments for each mutational scan")
    parser.add_argument("-mapping_file", type=str, help="Path to save the mapping file to")
    args = parser.parse_args()
    create_mapfile(args.input_folder, args.mapping_file, args.alignment_folder)
    