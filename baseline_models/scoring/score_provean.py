import argparse 
import os 
import pandas as pd 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re 
import sys 
sys.path.append("../../preprocessing")
from utils.mutated_sequence_utils import *
amino_acid_map = {
    "Ala":"A",
    "Asx":"B",
    "Cys":"C",
    "Asp":"D",
    "Glu":"E",
    "Phe":"F",
    "Gly":"G",
    "His":"H",
    "Ile":"I",
    "Lys":"K",
    "Leu":"L",
    "Met":"M",
    "Asn":"N",
    "Pro":"P",
    "Gln":"Q",
    "Arg":"R",
    "Ser":"S",
    "Thr":"T",
    "Val":"V",
    "Trp":"W",
    "Tyr":"Y",
    "Glx":"Z",
}

def get_mutated_sequence(focus_seq, mutant, start_idx=1, AA_vocab="ACDEFGHIKLMNPQRSTVWY"):
    """
    Helper function that mutates an input sequence (focus_seq) via an input mutation triplet (substitutions only).
    Mutation triplet are typically based on 1-indexing: start_idx is used for switching to 0-indexing.
    """
    mutated_seq = list(focus_seq)
    for mutation in mutant.split(":"):
        try:
            from_AA, position, to_AA = mutation[0], int(mutation[1:-1]), mutation[-1]
        except:
            print("Issue with mutant: "+str(mutation))
        relative_position = position - start_idx
        assert (from_AA==focus_seq[relative_position]), "Invalid from_AA or mutant position: "+ str(mutation)+" from_AA: "+ str(from_AA) + " relative pos: "+ str(relative_position) + " focus_seq: "+ str(focus_seq)
        assert (to_AA in AA_vocab) , "Mutant to_AA is invalid: "+ str(mutation)
        mutated_seq[relative_position] = to_AA
    return "".join(mutated_seq)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Provean scoring')
    parser.add_argument('--provean_script_path', type=str, help='Path to provean scoring bash script')
    parser.add_argument('--dataset_reference_file', type=str, help='Path to map file for input dataset')
    parser.add_argument('--dataset_folder', type=str, help='Path of folder containing containing the variants for each target sequence')
    parser.add_argument('--fasta_folder', default="./fasta_files", type=str, help="Path to write fasta formatted proteins to so they can be passed in to Provean")
    parser.add_argument('--keep_fasta', action='store_true', help="If passed, fasta files generated will not be deleted at end of run")
    parser.add_argument('--target_sequence_index', default=None, type=int, help='Index of target sequence in reference dataframe')
    parser.add_argument('--target_sequence_range_start', default=None, type=int, help="Start of range of target sequence indexes to run (used in place of a single target sequence index)")
    parser.add_argument('--target_sequence_range_end', default=None, type=int, help="End of range of target sequence indexes to run (used in place of a single target sequence index)")
    parser.add_argument('--output_raw_scores_folder', default="./output_scores", type=str, help='Name of folder to write model scores to')
    parser.add_argument("--output_scores_folder", type=str, help="Directory to save scores to after processing Provean's output")
    parser.add_argument('--supporting_set_folder',type=str, default="./output_supporting_sets")
    parser.add_argument('--indel_mode', action='store_true', help='Whether to score sequences with insertions and deletions')
    parser.add_argument('--name_column', default='Name', type=str, help='Column containing hgvs descriptions of variants')
    parser.add_argument('--convert_three_letter_to_single', default=True, type=bool, help="Whether to convert hgvs annotations from three letter abbreviations to single letter ones")
    parser.add_argument('--num_threads', type=int, default=1, help="Number of threads to run provean with")
    parser.add_argument("--label_column", type=str, default="label", help="Column in csv score files containing the label")
    args = parser.parse_args()

    if not os.path.isdir(args.fasta_folder):
        os.makedirs(args.fasta_folder, exist_ok=True)
    if not os.path.isdir(args.output_scores_folder):
        os.makedirs(args.output_scores_folder, exist_ok=True)
    if not os.path.isdir(args.supporting_set_folder):
        os.makedirs(args.supporting_set_folder, exist_ok=True)
    if not os.path.isdir(args.output_raw_scores_folder):
        os.makedirs(args.output_raw_scores_folder, exist_ok=True)

    target_seq_mapfile = pd.read_csv(args.dataset_reference_file)
    list_target_sequences = target_seq_mapfile["target_sequence_id"]
    if args.target_sequence_index != None:
        target_sequence_ids = pd.Series([list_target_sequences[args.target_sequence_index]])
    elif args.target_sequence_range_start != None and args.target_sequence_range_end != None:
        target_sequence_ids = list_target_sequences[args.target_sequence_range_start:args.target_sequence_range_end]
    else:
        print("Error: Must specify either target_sequence_index or target_sequence_range_start and target_sequence_range_end in command line arguments")
        exit()
    print("Computing scores for: {} with Provean".format(target_sequence_ids))
    target_sequence_file_names = [target_seq_mapfile["target_sequence_filename"][target_seq_mapfile["target_sequence_id"] == target_sequence_id].values[0] for target_sequence_id in target_sequence_ids]
    target_sequence_file_names_partial = [os.path.splitext(target_sequence_file_name)[0] for target_sequence_file_name in target_sequence_file_names]
    target_seqs = [target_seq_mapfile["target_seq"][target_seq_mapfile["target_sequence_id"] == target_sequence_id].values[0].upper() for target_sequence_id in target_sequence_ids]
    target_sequence_ids = target_sequence_ids.reset_index(drop=True)
    for i in range(len(target_sequence_ids)):
        variant_data = pd.read_csv(args.variant_data_folder + os.sep + target_sequence_file_names[i], low_memory=False)
        if variant_data[args.name_column].values[0][:3] == "(p." and variant_data[args.name_column].values[0][-1] == ")":
            variant_data["variations"] = variant_data[args.name_column].str.extract("\(p\.(.*)\)")
        elif variant_data[args.name_column].values[0][:3] == "p.(" and variant_data[args.name_column].values[0][-1] == ")":
            variant_data["variations"] = variant_data[args.name_column].str.extract("p\.\((.*)\)")
        else:
            variant_data["variations"] = variant_data[args.name_column]
        if args.convert_three_letter_to_single:
            for key in amino_acid_map:
                variant_data["variations"] = variant_data["variations"].str.replace(key, amino_acid_map[key])
        wild_type_seq = SeqRecord(Seq(target_seqs[i]),id=target_sequence_ids[i])
        variant_filename = f"{args.fasta_folder}/{target_sequence_file_names_partial[i]}.var"
        wild_type_filename = f"{args.fasta_folder}/{target_sequence_file_names_partial[i]}.fasta"
        with open(variant_filename,"w+") as file:
            for j in range(len(variant_data)):
                file.write(variant_data.iloc[j]["variations"] + "\n")
        SeqIO.write(wild_type_seq, wild_type_filename, "fasta")
        if os.path.exists(f"{args.supporting_set_folder}/{target_sequence_file_names_partial[i]}.sss"):
            print(f"Found saved supporting set for {target_sequence_file_names_partial[i]}")
            os.system(f"bash {args.provean_script_path} -q {wild_type_filename} -v {variant_filename} --supporting_set {args.supporting_set_folder}/{target_sequence_file_names_partial[i]}.sss --num_threads {args.num_threads} > {args.output_raw_scores_folder}/{target_sequence_file_names_partial[i]}.txt")
        else:
            os.system(f"bash {args.provean_script_path} -q {wild_type_filename} -v {variant_filename} --save_supporting_set {args.supporting_set_folder}/{target_sequence_file_names_partial[i]}.sss --num_threads {args.num_threads} > {args.output_raw_scores_folder}/{target_sequence_file_names_partial[i]}.txt")
        if not args.keep_fasta:
            os.remove(wild_type_filename)
            os.remove(variant_filename)
        
        target_sequence_name = target_sequence_file_names_partial[i]
        target_seq = target_seq_mapfile[target_seq_mapfile['target_sequence_id'] == target_sequence_name]["target_seq"].values[0]
        result_dict = {"mutant":[],"provean_score":[],"mutated_sequence":[]}
        if args.label_column in variant_data.columns:
            result_dict["label"] = variant_data[args.label_column]
        with open(f"{args.output_raw_scores_folder}/{target_sequence_file_names_partial[i]}.txt","r") as score_file:
            lines = score_file.readlines() 
            score_start = lines.index('# VARIATION\tSCORE\n')
            scores = lines[score_start+1:]
            for score in scores:
                variant, val = score.strip().split("\t")
                result_dict["provean_score"].append(val)
                result_dict["mutant"].append(variant)
                if "delins" in variant:
                    mutated_sequence = delins(variant, target_seq)
                elif "ins" in variant:
                    mutated_sequence = insert(variant, target_seq)
                elif "del" in variant:
                    mutated_sequence = delete(variant, target_seq)
                elif "dup" in variant:
                    mutated_sequence = duplicate(variant, target_seq)
                elif re.match("\w\d*\w",variant):
                    mutated_sequence = substitution(variant, target_seq)
                else:
                    print(variant)
                    print("ERROR")
                    exit()
                result_dict["mutated_sequence"].append(mutated_sequence)
            df = pd.DataFrame(result_dict)
            df["model_score"] = df["provean_score"]
            assert len(df) != 0
            df.to_csv(args.output_scores_folder + os.sep + target_sequence_name + ".csv", index=False)
    