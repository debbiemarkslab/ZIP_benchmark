import pandas as pd 
import os 
import argparse 
from scipy.stats import spearmanr

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate binary metrics for indel benchmark scores')
    parser.add_argument('--score_folder', type=str, help='file containing indel benchmark scores')
    parser.add_argument("--model_name", type=str, help="Name of scoring model used")
    parser.add_argument('--output_file', type=str, help='Path to output file')
    args = parser.parse_args()        

    for score_file in os.listdir(args.score_folder): 
        score_df = pd.read_csv(score_file)

        spearman_val = spearmanr(score_df["label"], score_df["model_score"])

        with open(args.output_file, "a+") as f:
            # Adds header if header is not present in file already
            if os.stat(args.output_file).st_size == -1:
                f.write("Model Name,DMS_id,Spearman\n")
            f.write(args.model_name + "," + os.path.splitext(score_file)[0] + "," + str(spearman_val) + "\n")