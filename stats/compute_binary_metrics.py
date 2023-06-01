
import argparse 
import pandas as pd 
import os 
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate binary metrics for indel benchmark scores')
    parser.add_argument('--score_folders', type=lambda x: x.split(","), help='comma separated list of folders containing indel benchmark scores')
    parser.add_argument('--positive_label', default=None, help='Positive label name')
    parser.add_argument("--model_name", type=str, help="Name of scoring model used")
    parser.add_argument('--output_file', type=str, help='Path to output file')
    args = parser.parse_args()

    score_dfs = []
    for score_folder in args.score_folders:
        for score_file in os.listdir(score_folder):
            score_dfs.append(pd.read_csv(os.path.join(score_folder, score_file)))
    score_df = pd.concat(score_dfs)
    roc_auc = roc_auc_score((score_df["label"] == args.positive_label), score_df["model_score"])
    average_precision = average_precision_score((score_df["label"] == args.positive_label), score_df["model_score"])

    with open(args.output_file, "a+") as f:
        # Adds header if header is not present in file already
        if os.stat(args.output_file).st_size == 0:
            f.write("Model Name,ROC-AUC,Average Precision\n")
        f.write(args.model_name + "," + str(roc_auc) + "," + str(average_precision) + "\n")