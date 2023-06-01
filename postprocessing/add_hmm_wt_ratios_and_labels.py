import os 
import pandas as pd 
import argparse 

"""
this script takes in a list of folders with individual per-protein scores for the HMM model and adds
the wild type ratio column (difference in predicted log-prob between wild type and mutant sequences) and the label column. 
This is done after the scoring here (rather than simply including the ratio and labels in original outputs) because the HMM
scoring is done using jackhmmer's internal forward-backward algorithm, which does not output the wild type ratio or allow for the output of 
extraneous columns. 
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add wild type ratio and label columns to per-protein scores')
    parser.add_argument('--hmm_score_folder', type=str, help='folder containing per-protein scores for HMM model')
    parser.add_argument('--output_folder', type=str, help='Path to output post-processed per-protein scores')
    parser.add_argument("--dataset_folder", type=str, help="Path to folder containing per-protein variants that were scored")
    parser.add_argument("--label_column", default="label", type=str, help="Name of column in dataset_folder csvs that contains labels")
    parser.add_argument("--mutated_sequence_column", default="protein_sequence_mutated", type=str, help="Name of column in dataset_folder csvs that contains mutated sequences")
    args = parser.parse_args()

    assert os.path.isdir(args.hmm_score_folder), f"{args.hmm_score_folder} is not a valid folder"
    if not os.path.isdir(args.output_folder):
        os.mkdir(args.output_folder)

    
    for file in os.listdir(args.hmm_score_folder):
        if not os.path.exists(os.path.join(args.dataset_folder,file)):
            print(f"Warning: {file} not found in {args.dataset_folder}. Skipping.") 
            continue
        df = pd.read_csv(os.path.join(args.hmm_score_folder,file))
        # removing extra whitespace from seq_name column
        df["seq_name"] = df["seq_name"].apply(lambda x: x.replace(" ",""))
        wt_logprob = df[df["seq_name"] == "WT"]["logprob"].values[0]
        df["model_score"] = df["logprob"].astype(float) - float(wt_logprob)
        # Dropping WT score here to match other output formats 
        df = df[df["seq_name"] != "WT"]
        df = df.rename(columns={"seq":"mutated_sequence"})
        variant_df = pd.read_csv(os.path.join(args.dataset_folder,file)).rename(columns={args.mutated_sequence_column:"mutated_sequence"})
        df = df.merge(variant_df[["mutated_sequence",args.label_column]], on="mutated_sequence", how="left")
        if args.label_column != "label":
            df = df.rename(columns={args.label_column:"label"})
        df.to_csv(os.path.join(args.output_folder,file),index=False)