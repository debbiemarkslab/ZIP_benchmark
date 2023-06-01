import argparse 
import os 
from tqdm import tqdm 

def get_focus_columns_only(alignment_file, output_file):
    file_handle = open(alignment_file)
    out_handle = open(output_file,"a+")
    for line in file_handle.readlines():
        if line.startswith(">"):
            out_handle.write(line)
        else:
            out = "".join([char for char in line if (not char.islower() and not char == ".")])
            out_handle.write(out)
    file_handle.close()
    out_handle.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get focus columns only from alignment file')
    parser.add_argument('--alignment_folder', type=str, help='Alignment file')
    parser.add_argument('--output_folder', type=str, help='Output file')
    args = parser.parse_args()
    for file in tqdm(os.listdir(args.alignment_folder)):
        basename = os.path.splitext(file)[0]
        alignment_file = f"{args.alignment_folder}/{file}"
        output_file = f"{args.output_folder}/{basename}.a2m"
        get_focus_columns_only(alignment_file, output_file)