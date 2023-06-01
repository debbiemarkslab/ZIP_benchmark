import os
import argparse
import json
import pandas as pd
import pickle 
import torch
import sys 
import numpy as np 
sys.path.append(os.path.dirname(os.path.dirname(__file__)) + "/models")

from transformers import PreTrainedTokenizerFast
import tranception
from tranception import config, model_pytorch


def main():
    """
    Main script to score sets of mutated protein sequences (substitutions or indels) with Tranception.
    """
    parser = argparse.ArgumentParser(description='Tranception scoring')
    parser.add_argument('--checkpoint', type=str, help='Path of Tranception model checkpoint')
    parser.add_argument('--model_framework', default='pytorch', type=str, help='Underlying framework [pytorch|JAX]')
    parser.add_argument('--batch_size_inference', default=20, type=int, help='Batch size for inference')

    #We may pass in all required information about the dataset via the provided reference files, or specify all relevant fields manually
    parser.add_argument('--dataset_reference_file', default=None, type=str, help='Path to reference file with list of target sequences (and filepaths to their associated variants) to score')
    parser.add_argument('--target_sequence_index', default=None, type=int, help='Index of sequence and variants to score in reference file')
    parser.add_argument('--target_sequence_index_range_start', default=None, type=int, help="Start of range of experiment indexes to run (used in place of a single experiment index)")
    parser.add_argument('--target_sequence_index_range_end', default=None, type=int, help="End of range of experiment indexes to run (used in place of a single experiment index)")
    #Fields to be passed manually if reference file is not used
    parser.add_argument('--target_seq', default=None, type=str, help='Full wild type sequence that is mutated in the experiment')
    parser.add_argument('--target_sequence_file_name', default=None, type=str, help='Name of experiment file')
    parser.add_argument('--MSA_filename', default=None, type=str, help='Name of MSA (eg., a2m) file constructed on the wild type sequence')
    parser.add_argument('--MSA_weight_file_name', default=None, type=str, help='Weight of sequences in the MSA (optional)')
    parser.add_argument('--MSA_start', default=None, type=int, help='Sequence position that the MSA starts at (1-indexing)')
    parser.add_argument('--MSA_end', default=None, type=int, help='Sequence position that the MSA ends at (1-indexing)')
    parser.add_argument("--mutant_column",default="mutant",type=str)
    parser.add_argument("--mutated_sequence_column",default="mutated_sequence",type=str)

    parser.add_argument('--dataset_folder', type=str, help='Path to folder that contains the protein variants for each target sequence')
    parser.add_argument('--output_scores_folder', default='./', type=str, help='Name of folder to write model scores to')
    parser.add_argument('--deactivate_scoring_mirror', action='store_true', help='Whether to deactivate sequence scoring from both directions (Left->Right and Right->Left)')
    parser.add_argument('--indel_mode', action='store_true', help='Flag to be used when scoring insertions and deletions. Otherwise assumes substitutions')
    parser.add_argument('--scoring_window', default=None, type=str, help='Sequence window selection mode (when sequence length longer than model context size)')
    parser.add_argument('--num_workers', default=10, type=int, help='Number of workers for model scoring data loader')
    parser.add_argument('--inference_time_retrieval', action='store_true', help='Whether to perform inference-time retrieval')
    parser.add_argument('--retrieval_inference_weight', default=0.6, type=float, help='Coefficient (alpha) used when aggregating autoregressive transformer and retrieval')
    parser.add_argument('--MSA_folder', default='.', type=str, help='Path to MSA for neighborhood scoring')
    parser.add_argument('--MSA_weights_folder', default=None, type=str, help='Path to MSA weights for neighborhood scoring')
    parser.add_argument('--clustal_omega_location', default=None, type=str, help='Path to Clustal Omega (only needed with scoring indels with retrieval)')
    parser.add_argument("--clustal_suffix", default="", type=str, help="A suffix added when generating clustal files. Only needed if you are using the same alignment for multiple concurrent jobs, where overwriting may be a concern.")
    parser.add_argument("--score_wt_only", action="store_true", help="Will only score wildtype sequences if passed")
    parser.add_argument("--tokenizer_filepath", type=str, default="./tranception/utils/tokenizers/Basic_tokenizer", help="Path to tokenizer file")
    parser.add_argument("--label_column", type=str, default="label", help="Name of column in dataset file that contains the label (optional)")
    args = parser.parse_args()
    
    model_name = args.checkpoint.split("/")[-1]
    
    tokenizer = PreTrainedTokenizerFast(tokenizer_file=args.tokenizer_filepath,
                                                unk_token="[UNK]",
                                                sep_token="[SEP]",
                                                pad_token="[PAD]",
                                                cls_token="[CLS]",
                                                mask_token="[MASK]"
                                            )

    if args.dataset_reference_file:
        target_seq_mapfile = pd.read_csv(args.dataset_reference_file)
        list_target_sequences = target_seq_mapfile["target_sequence_id"].tolist()
        if args.target_sequence_index != None:
            target_sequence_ids=[list_target_sequences[args.target_sequence_index]]
            print("Compute scores for experiment: "+str(target_sequence_ids))
            target_seqs = [target_seq_mapfile["target_seq"][target_seq_mapfile["target_sequence_id"]==target_sequence_ids[0]].values[0].upper()]
            target_sequence_file_names = [target_seq_mapfile["target_sequence_filename"][target_seq_mapfile["target_sequence_id"]==target_sequence_ids[0]].values[0]]
            if args.inference_time_retrieval:
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
            print("Computing Scores for: " + str(target_sequence_ids))
            target_seqs = [target_seq_mapfile["target_seq"][target_seq_mapfile["target_sequence_id"] == target_sequence_id].values[0].upper() for target_sequence_id in target_sequence_ids]
            target_sequence_file_names = [target_seq_mapfile["target_sequence_filename"][target_seq_mapfile["target_sequence_id"] == target_sequence_id].values[0] for target_sequence_id in target_sequence_ids]
            if args.inference_time_retrieval:
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
        if args.inference_time_retrieval:
            MSA_data_files = [args.MSA_folder + os.sep + args.MSA_filename if args.MSA_folder is not None else None]
            MSA_weight_file_names = [args.MSA_weights_folder + os.sep + args.MSA_weight_file_name if args.MSA_weights_folder is not None else None]
            MSA_starts = [args.MSA_start - 1] # MSA_start based on 1-indexing
            MSA_ends = [args.MSA_end]
        mutant_columns = [args.mutant_column]
        mutated_sequence_columns = [args.mutated_sequence_column]


    config = json.load(open(args.checkpoint+os.sep+'config.json'))
    config = tranception.config.TranceptionConfig(**config)
    config.attention_mode="tranception"
    config.position_embedding="grouped_alibi"
    config.tokenizer = tokenizer
    config.clustal_suffix = args.clustal_suffix
    # Iterating through each target protein sequence
    for i in range(len(target_sequence_ids)):
        print(f"Scoring: {target_sequence_ids[i]}")
        if args.indel_mode and len(target_seqs[i]) > 1022:
            if args.scoring_window == "optimal":
                print(f"Warning: Target sequence for {target_sequence_ids[i]} is greater than 1022 and indel mode is enabled, but scoring window is set to 'optimal', which will miss part of the sequence. Use scoring window 'sliding' to account for extra sequence length.")
                config.scoring_window = args.scoring_window
            else:
                print(f"Target sequence for {target_sequence_ids[i]} is greater than 1022 and indel mode is enabled. Setting scoring window to 'sliding' to account for extra sequence length.")
                config.scoring_window = "sliding" 
        else:
            config.scoring_window = 'optimal'
        if args.inference_time_retrieval and MSA_data_files[i] != None:
            config.retrieval_aggregation_mode = "aggregate_indel" if args.indel_mode else "aggregate_substitution"
            config.MSA_filename=MSA_data_files[i]
            config.full_protein_length=len(target_seqs[i])
            config.MSA_weight_file_name=MSA_weight_file_names[i]
            config.retrieval_inference_weight=args.retrieval_inference_weight
            config.MSA_start = MSA_starts[i]
            config.MSA_end = MSA_ends[i]
            if args.indel_mode:
                config.clustal_omega_location = args.clustal_omega_location
        else:
            if args.inference_time_retrieval:
                print(f"Warning: Retrieval set but alignment file not present for {MSA_data_files[i]}. Running without retrieval")
            config.retrieval_aggregation_mode = None
        
        device = torch.device(f"cuda:0" if torch.cuda.is_available() else "cpu")
            
        if args.model_framework=="pytorch":
            model = tranception.model_pytorch.TranceptionLMHeadModel.from_pretrained(pretrained_model_name_or_path=args.checkpoint,config=config)
            if torch.cuda.is_available():
                model.cuda()
        model.eval()
        
        if not os.path.isdir(args.output_scores_folder):
            os.makedirs(args.output_scores_folder,exist_ok=True)
        scoring_filename = args.output_scores_folder
        scoring_filename += os.sep + target_sequence_ids[i] + '.csv'
        
        if not args.score_wt_only:
            variant_data = pd.read_csv(args.dataset_folder + os.sep + target_sequence_file_names[i], low_memory=False)
            if mutant_columns[i] not in variant_data.columns:
                if mutated_sequence_columns[i] not in variant_data.columns:
                    print(f"Error: Neither {mutant_columns[i]} nor {mutated_sequence_columns[i]} are in the protein variant dataframe for {target_sequence_ids[i]}")
            if mutant_columns[i] in variant_data.columns:
                variant_data = variant_data.rename(columns={mutant_columns[i]:"mutant"})
            if mutated_sequence_columns[i] in variant_data.columns:
                variant_data = variant_data.rename(columns={mutated_sequence_columns[i]:"mutated_sequence"})

            # If protein variant file was saved with an index, it gets set to Unnamed: 0, so we remove that here (causes errors with too much memory allocation otherwise)
            if args.label_column in variant_data.columns:
                labels = variant_data[args.label_column]
            else:
                print("Warning: No label column found in protein variant file. Labels will not be output with scores")
            if "Unnamed: 0" in variant_data.columns:
                variant_data = variant_data.drop(["Unnamed: 0"], axis=1)
            if "mutant" in variant_data.columns:
                if "mutated_sequence" in variant_data.columns:
                    variant_data_subset = variant_data[["mutant", "mutated_sequence"]]
                else:
                    variant_data_subset = variant_data[["mutant"]]
            else:
                if "mutated_sequence" in variant_data.columns:
                    variant_data_subset = variant_data[["mutated_sequence"]]
                else:
                    print("Must pass in either mutant or mutated_sequence columns with protein variant file")
                    exit() 
            all_scores = model.score_mutants(
                            DMS_data=variant_data_subset, 
                            target_seq=target_seqs[i], 
                            scoring_mirror=not args.deactivate_scoring_mirror, 
                            batch_size_inference=args.batch_size_inference,  
                            num_workers=args.num_workers, 
                            indel_mode=args.indel_mode,
                            )
            all_scores["model_score"] = all_scores["avg_score"]
            if args.label_column in variant_data.columns:
                all_scores["label"] = labels
            all_scores.to_csv(scoring_filename, index=False)
        else:
            variant_data_subset = pd.DataFrame({"mutated_sequence":[target_seqs[i]]})
            all_scores = model.score_mutants(
                DMS_data=variant_data_subset,
                target_seq=None,
                scoring_mirror=not args.deactivate_scoring_mirror,
                batch_size_inference=args.batch_size_inference,
                num_workers=args.num_workers,
                indel_mode=args.indel_mode,
                )
            all_scores["model_score"] = all_scores["avg_score"]
            all_scores.to_csv(scoring_filename, index=False)

if __name__ == '__main__':
    main()
