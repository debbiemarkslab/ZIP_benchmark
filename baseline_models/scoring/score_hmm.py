import argparse 
import pandas as pd 
import os 
import numpy as np 
import glob 
from Bio.SeqIO.FastaIO import SimpleFastaParser 
"""
This script scores a folder of variants using the HMM model.
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Score a folder of variants using the HMM model')   
    #We may pass in all required information about the dataset via the provided reference files, or specify all relevant fields manually
    parser.add_argument('--dataset_reference_file', type=str, help='Path to reference file with list of target sequences (and filepaths to their associated variants) to score')
    parser.add_argument('--target_sequence_index', type=int, help='Index of sequence and variants to score in reference file')
    parser.add_argument('--target_sequence_index_range_start', type=int, help="Start of range of experiment indexes to run (used in place of a single experiment index)")
    parser.add_argument('--target_sequence_index_range_end', type=int, help="End of range of experiment indexes to run (used in place of a single experiment index)")
    parser.add_argument("--hmmer_path", type=str, help="Path to hmmer installation")
    parser.add_argument("--reformat_script_path", type=str, help="Path to folder containing reformat.pl script")
    parser.add_argument('--dataset_folder', type=str, help='Path to folder that contains the protein variants for each target sequence')
    parser.add_argument('--output_scores_folder', default='./', type=str, help='Name of folder to write model scores to')
    parser.add_argument("--intermediate_outputs_folder", type=str, default="./intermediate_outputs", help="Path to folder to write intermediate outputs to")
    parser.add_argument('--MSA_folder', default='.', type=str, help='Path to MSA for neighborhood scoring')
    parser.add_argument('--MSA_weights_folder', default=None, type=str, help='Path to MSA weights for neighborhood scoring')
    parser.add_argument("--label_column", type=str, default="label", help="Name of column in dataset file that contains the label (optional)")

    #Fields to be passed manually if reference file is not used
    parser.add_argument('--target_seq', default=None, type=str, help='Full wild type sequence that is mutated in the experiment')
    parser.add_argument('--target_sequence_file_name', default=None, type=str, help='Name of experiment file')
    parser.add_argument('--MSA_filename', default=None, type=str, help='Name of MSA (eg., a2m) file constructed on the wild type sequence')
    parser.add_argument('--MSA_weight_file_name', default=None, type=str, help='Weight of sequences in the MSA (optional)')
    parser.add_argument('--MSA_start', default=None, type=int, help='Sequence position that the MSA starts at (1-indexing)')
    parser.add_argument('--MSA_end', default=None, type=int, help='Sequence position that the MSA ends at (1-indexing)')
    parser.add_argument("--mutant_column",default="mutant",type=str)
    parser.add_argument("--mutated_sequence_column",default="mutated_sequence",type=str)
    args = parser.parse_args()

    if not os.path.exists(args.output_scores_folder):
        os.makedirs(args.output_scores_folder)
    if not os.path.exists(args.intermediate_outputs_folder):    
        os.makedirs(args.intermediate_outputs_folder)
    if not os.path.exists(args.intermediate_outputs_folder + os.sep + "hmm"):
        os.mkdir(args.intermediate_outputs_folder + os.sep + "hmm")
    if not os.path.exists(args.intermediate_outputs_folder + os.sep + "fa"):
        os.mkdir(args.intermediate_outputs_folder + os.sep + "fa")
    if not os.path.exists(args.intermediate_outputs_folder + os.sep + "a2m"):
        os.mkdir(args.intermediate_outputs_folder + os.sep + "a2m")
    if not os.path.exists(args.intermediate_outputs_folder + os.sep + "raw_scores"):
        os.mkdir(args.intermediate_outputs_folder + os.sep + "raw_scores")
    hmm_path = args.intermediate_outputs_folder + os.sep + "hmm"
    fa_path = args.intermediate_outputs_folder + os.sep + "fa"
    raw_score_path = args.intermediate_outputs_folder + os.sep + "raw_scores"
    
    if args.dataset_reference_file:
        target_seq_mapfile = pd.read_csv(args.dataset_reference_file)
        list_target_sequences = target_seq_mapfile["target_sequence_id"].tolist()
        if args.target_sequence_index != None:
            target_sequence_ids=[list_target_sequences[args.target_sequence_index]]
            print(f"Computing HMM scores for {len(target_sequence_ids)} target sequences)")
            target_seqs = [target_seq_mapfile["target_seq"][target_seq_mapfile["target_sequence_id"]==target_sequence_ids[0]].values[0].upper()]
            target_sequence_file_names = [target_seq_mapfile["target_sequence_filename"][target_seq_mapfile["target_sequence_id"]==target_sequence_ids[0]].values[0]]
            MSA_data_files = [target_seq_mapfile["MSA_filename"].tolist()[args.target_sequence_index]]
            MSA_data_files = [args.MSA_folder + os.sep + MSA_data_file if type(MSA_data_file) != float else None for MSA_data_file in MSA_data_files]
            MSA_weight_file_names = [target_seq_mapfile["weight_file_name"].tolist()[args.target_sequence_index]]
            MSA_weight_file_names = [args.MSA_weights_folder + os.sep + MSA_weight_file_name if type(MSA_weight_file_name) != float else None for MSA_weight_file_name in MSA_weight_file_names]
            MSA_starts = [target_seq_mapfile["MSA_start"].tolist()[args.target_sequence_index]]
            MSA_starts = [int(MSA_start) - 1 if not np.isnan(MSA_start) else None for MSA_start in MSA_starts] # MSA_start typically based on 1-indexing 
            MSA_ends = [target_seq_mapfile["MSA_end"].tolist()[args.target_sequence_index]]
            MSA_ends = [int(MSA_end) if not np.isnan(MSA_end) else None for MSA_end in MSA_ends]
            if "mutant_column" in target_seq_mapfile.columns:
                mutant_columns = [target_seq_mapfile["mutant_column"][args.target_sequence_index]]
            else:
                mutant_columns = [args.mutant_column]
            if "mutated_sequence_column" in target_seq_mapfile.columns:
                mutated_sequence_columns = [target_seq_mapfile["mutated_sequence_column"][args.target_sequence_index]]
            else:
                mutated_sequence_columns = [args.mutated_sequence_column]
        elif args.target_sequence_index_range_start != None and args.target_sequence_index_range_end != None:
            target_sequence_ids = list_target_sequences[args.target_sequence_index_range_start:args.target_sequence_index_range_end]
            print(f"Computing HMM scores for {len(target_sequence_ids)} target sequences")
            target_seqs = [target_seq_mapfile["target_seq"][target_seq_mapfile["target_sequence_id"] == target_sequence_id].values[0].upper() for target_sequence_id in target_sequence_ids]
            target_sequence_file_names = [target_seq_mapfile["target_sequence_filename"][target_seq_mapfile["target_sequence_id"] == target_sequence_id].values[0] for target_sequence_id in target_sequence_ids]
            MSA_data_files = target_seq_mapfile["MSA_filename"].tolist()[args.target_sequence_index_range_start:args.target_sequence_index_range_end]
            MSA_data_files = [args.MSA_folder + os.sep + MSA_data_file if type(MSA_data_file) != float else None for MSA_data_file in MSA_data_files]
            MSA_weight_file_names = target_seq_mapfile["weight_file_name"].tolist()[args.target_sequence_index_range_start:args.target_sequence_index_range_end]
            MSA_weight_file_names = [args.MSA_weights_folder + os.sep + MSA_weight_file_name if type(MSA_weight_file_name) != float else None for MSA_weight_file_name in MSA_weight_file_names]
            MSA_starts = target_seq_mapfile["MSA_start"].tolist()[args.target_sequence_index_range_start:args.target_sequence_index_range_end]
            MSA_starts = [int(MSA_start) - 1 if not np.isnan(MSA_start) else None for MSA_start in MSA_starts] # MSA_start typically based on 1-indexing 
            MSA_ends = target_seq_mapfile["MSA_end"].tolist()[args.target_sequence_index_range_start:args.target_sequence_index_range_end]
            MSA_ends = [int(MSA_end) if not np.isnan(MSA_end) else None for MSA_end in MSA_ends]
            if "mutant_column" in target_seq_mapfile.columns:
                mutant_columns = [target_seq_mapfile["mutant_column"][i] for i in range(args.target_sequence_index_range_start, args.target_sequence_index_range_end)]
            else:
                mutant_columns = [args.mutant_column for i in range(args.target_sequence_index_range_start, args.target_sequence_index_range_end)]
            if "mutated_sequence_column" in target_seq_mapfile.columns:
                mutated_sequence_columns = [target_seq_mapfile["mutated_sequence_column"][i] for i in range(args.target_sequence_index_range_start, args.target_sequence_index_range_end)]
            else:
                mutated_sequence_columns = [args.mutated_sequence_column for i in range(args.target_sequence_index_range_start, args.target_sequence_index_range_end)]
        else:
            raise ValueError("Must pass in either a target_sequence_index value or values for target_sequence_index_range_start and target_sequence_index_range_end")
    else:
        target_seqs=[args.target_seq]
        target_sequence_file_names=[args.target_sequence_file_name]
        target_sequence_ids = [target_sequence_file_names.split(".")[0]]
        MSA_data_files = [args.MSA_folder + os.sep + args.MSA_filename if args.MSA_folder is not None else None]
        MSA_weight_file_names = [args.MSA_weights_folder + os.sep + args.MSA_weight_file_name if args.MSA_weights_folder is not None else None]
        MSA_starts = [args.MSA_start - 1] # MSA_start based on 1-indexing
        MSA_ends = [args.MSA_end]
        mutant_columns = [args.mutant_column]
        mutated_sequence_columns = [args.mutated_sequence_column]

    # checking all a3m files in alignment and calling reformat if an a2m or hmm file are not already generated 
    print("Checking for alignment files and reformatting a3m files if needed")
    a2m_MSA_files = []
    for i,MSA_data_file in enumerate(MSA_data_files):
        if MSA_data_file is None:
            print(f"No alignment found for {target_sequence_ids[i]}. Skipping.")
            continue 
        basename = os.path.splitext(MSA_data_file)[0]
        if os.path.exists(hmm_path + os.sep + basename + ".hmm"):
            continue
        if MSA_data_file.endswith(".a2m"):
            if os.path.exists(args.MSA_folder + os.sep + MSA_data_file):
                a2m_MSA_files.append(basename + ".a2m")
                continue
            else:
                if os.path.exists(args.MSA_folder + os.sep + basename + ".a3m"):
                    print(f"Warning: {MSA_data_file} not found. Reformatting {basename + '.a3m'} instead.")
                    os.system(f"{args.reformat_script_path}/reformat.pl a3m a2m {args.MSA_folder + os.sep + basename + '.a3m'} {args.MSA_folder + os.sep + basename + '.a2m'}")
                    a2m_MSA_files.append(basename + ".a2m")
    # Building HMMs from a2ms if HMMs do not already exist
    print("Building HMMs if they don't already exist")
    for a2m_MSA_file in a2m_MSA_files:
        print("Building HMM for " + a2m_MSA_file)
        basename = os.path.splitext(a2m_MSA_file)[0]
        if os.path.exists(hmm_path + os.sep + basename + ".hmm"):
            raise ValueError(f"{hmm_path + os.sep + basename + '.hmm'} already exists. This should not happen")
        os.system(f"{args.hmmer_path}" + os.sep + "bin" + os.sep + 'hmmbuild --amino ' + f"{hmm_path + os.sep + basename + '.hmm'} {args.MSA_folder + os.sep + a2m_MSA_file}")
    print("Writing out variants to fasta file for hmm scoring")
    for i, target_sequence_id in enumerate(target_sequence_ids):
        if not os.path.exists(args.dataset_folder + os.sep + target_sequence_id + ".csv"):
            print(f"Warning: {target_sequence_id}.csv not found in {args.dataset_folder}. Skipping.")
            continue
        variant_df = pd.read_csv(args.dataset_folder + os.sep + target_sequence_id + ".csv")
        with open(fa_path + os.sep + target_sequence_id + ".fa", "w+") as f:
            f.write(f">WT\n{target_seqs[i]}\n")
            for j, row in variant_df.iterrows():
                f.write(f">{row['HGVSp']}\n{row[mutated_sequence_columns[i]]}\n")
    print("Scoring variants with HMM forward-backward algorithm") 
    for i,target_sequence_id in enumerate(target_sequence_ids):
        hmm_file = hmm_path + os.sep + target_sequence_id + ".hmm"
        fa_file = fa_path + os.sep + target_sequence_id + ".fa"
        output_file = raw_score_path + os.sep + target_sequence_id + ".csv"
        os.system(f"{args.hmmer_path}" + os.sep + "src" + os.sep + f"generic_fwdback_example {hmm_file} {fa_file} > {output_file}")

    # Postprocessing scores to have label columns, wt ratio scores 
    print("Postprocessing scores to match output format")
    for i, target_sequence_id in enumerate(target_sequence_ids):
        if not os.path.exists(os.path.join(args.dataset_folder,target_sequence_id + ".csv")):
            print(f"Warning: {target_sequence_id + '.csv'} not found in {args.dataset_folder}. Skipping.") 
            continue
        df = pd.read_csv(os.path.join(raw_score_path,target_sequence_id + ".csv"))
        # removing extra whitespace from seq_name column
        df["seq_name"] = df["seq_name"].apply(lambda x: x.replace(" ",""))
        wt_logprob = df[df["seq_name"] == "WT"]["logprob"].values[0]
        df["model_score"] = df["logprob"].astype(float) - float(wt_logprob)
        # Occasionaly there are characters in sequences outside the alphabet of the hmm model, result in logprobs of -inf for both sequences and nan model scores. 
        # We set those NaNs equal to 0 here, assuming no difference between the mutants as it's outside the hmm alphabet. 
        df["model_score"] = df["model_score"].fillna(0)
        # Dropping WT score here to match other output formats 
        df = df[df["seq_name"] != "WT"]
        df = df.rename(columns={"seq":"mutated_sequence"})
        variant_df = pd.read_csv(os.path.join(args.dataset_folder,target_sequence_id + ".csv"))
        if args.mutated_sequence_column != "mutated_sequence":
            variant_df = variant_df.rename(columns={args.mutated_sequence_column:"mutated_sequence"})
        df = df.merge(variant_df[["mutated_sequence",args.label_column]], on="mutated_sequence", how="left")
        if args.label_column != "label":
            df = df.rename(columns={args.label_column:"label"})
        df.to_csv(os.path.join(args.output_scores_folder,target_sequence_id + ".csv"), index=False)