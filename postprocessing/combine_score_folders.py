import os 
import pandas as pd 
import argparse 

"""
this script takes in a list of folders with individual per-protein scores for a model 
and combines them into a single pandas dataframe with the same columns as the individual files.
This assumes that the input score files have columns 'protein_sequence_mutated', 'model_score', and 
'label'.
For example, we can pass in the score folder for both ClinVar_Tranception_Medium and gnomAD_Tranception_Medium 
to combine them into one file with all the scores.
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Combine per-protein scores into a single dataframe')
    parser.add_argument('--score_folder', type=str, help='Folder containing per-protein variant scores')
    parser.add_argument('--output_file', type=str, help='Path to output file')
    args = parser.parse_args()
    combined_dfs = [] 
    for folder in args.score_folders:
        assert os.path.isdir(folder), f"{folder} is not a valid folder"
        for file in os.listdir(folder):
            assert os.path.isfile(os.path.join(folder,file)), f"{file} is not a valid file"
            df = pd.read_csv(os.path.join(folder,file))
            assert "protein_sequence_mutated" in df.columns, f"protein_sequence_mutated not in {file}"
            assert "model_score" in df.columns, f"model_score not in {file}"
            assert "label" in df.columns, f"label not in {file}"
            combined_dfs.append(df[["protein_sequence_mutated","model_score","label"]])
    combined_df = pd.concat(combined_dfs)
    combined_df.to_csv(args.output_file,index=False)

    
