import argparse 
import os 
import pandas as pd 
import numpy as np 
from tqdm import tqdm 
"""
Because the CAPSD mutational scan is very large, we also provide this script to split a mutational scan into chunks for parallelization. 
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split CAPSD mutational scan into chunks')
    parser.add_argument('--input_file', type=str, help='Input DMS file')
    parser.add_argument('--output_dir', type=str, help='Output directory')
    parser.add_argument('--num_chunks', type=int, default=100, help='Size of each chunk')
    parser.add_argument("--target_sequence_id", type=str, help="DMS id to chunk")
    parser.add_argument("--original_mapping_file", type=str, help="Path to original mapping file for mutational scans")
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    DMS_data = pd.read_csv(args.input_file)
    split_dfs = np.array_split(DMS_data, args.num_chunks)
    for i in tqdm(range(args.num_chunks)):
        split_dfs[i].to_csv(os.path.join(args.output_dir, f'chunk_{i}.csv'), index=False)
    
    # Create mapping file for each chunk
    orig_map_df = pd.read_csv(args.original_mapping_file)
    orig_row = orig_map_df[orig_map_df["target_sequence_id"] == args.target_sequence_id]
    df_dict = {key:[] for key in orig_map_df.columns}
    df_dict["target_sequence_id"] = [f"{args.target_sequence_id}_chunk_{i}" for i in range(args.num_chunks)]
    df_dict["target_sequence_filename"] = [f"chunk_{i}.csv" for i in range(args.num_chunks)]
    for key in orig_map_df.columns:
        if key == "target_sequence_id" or key == "target_sequence_filename":
            continue
        df_dict[key] = [orig_row[key].values[0] for i in range(args.num_chunks)]
    chunk_map_df = pd.DataFrame(df_dict)

    chunk_map_df.to_csv(os.path.join(args.output_dir, f'mapfile_{args.target_sequence_id}_chunked.csv'), index=False)