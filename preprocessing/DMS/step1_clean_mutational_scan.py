import argparse 
import sys 
sys.path.append("..")
from utils.DMS_utils import DMS_file_cleanup 
from mutalyzer.protein import protein_description

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Clean up DMS data')
    parser.add_argument('--input_file', type=str, help='Input DMS file')
    parser.add_argument('--output_file', type=str, help='Output DMS file')
    parser.add_argument('--target_seq', type=str, help='Target sequence')
    parser.add_argument('--start_idx', type=int, default=1, help='Start index of target sequence')
    parser.add_argument('--end_idx', type=int, default=None, help='End index of target sequence')
    parser.add_argument('--DMS_mutant_column', type=str, default='mutant', help='Name of column containing mutant')
    parser.add_argument('--DMS_phenotype_name', type=str, default='score', help='Name of column containing phenotype')
    parser.add_argument('--DMS_directionality', type=int, default=1, help='Directionality of phenotype')
    parser.add_argument("--seq_suffix", default=None, type=str, help="If mutants are only a portion of the sequence, can pass in this suffix to add on to end of mutated sequences. (See Russ DMS for an example)")
    parser.add_argument('--AA_vocab', type=str, default="ACDEFGHIKLMNPQRSTVWY", help='Amino acid vocabulary')
    args = parser.parse_args()

    if args.end_idx is None:
        end_idx = args.start_idx + len(args.target_seq) - 1
    else:
        end_idx = args.end_idx 
    DMS_data = DMS_file_cleanup(args.input_file, args.target_seq, args.start_idx, end_idx, args.DMS_mutant_column, args.DMS_phenotype_name, args.DMS_directionality, args.AA_vocab)
    if args.seq_suffix != None:
        DMS_data["mutated_sequence"] = DMS_data["mutated_sequence"] + args.seq_suffix
    DMS_data["protein_sequence"] = args.target_seq
    # The first argument is the position of the stop codon, which is just used to determine if the mutation is inframe or not. All of 
    # ours are inframe, so we just pass in 3 which returns an inframe description. The indexing removes p.( and ) from the description.
    DMS_data["HGVSp"] = DMS_data["mutated_sequence"].apply(lambda x: protein_description(3, args.target_seq, x)[0][3:-1]) 
    DMS_data.to_csv(args.output_file, index=False)
