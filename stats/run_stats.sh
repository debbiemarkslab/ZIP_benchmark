#! /bin/bash
# Fail on first line failure
set -e

## Tranception Large retrieval 

# ClinVar/gnomAD
python compute_binary_metrics.py --score_folders="../outputs/output_scores/ClinVar_Tranception_Large_0.6_retrieval/,../outputs/output_scores/gnomAD_Tranception_Large_0.6_retrieval/" \
--positive_label="Benign" --model_name="Tranception_Large_0.6_retrieval" --output_file="./binary_stats.csv"
# DMS 
python compute_DMS_regression_metrics.py --score_folder="../outputs/output_scores/DMS_Tranception_Large_0.6_retrieval" \
--model_name="Tranception_Large_0.6_retrieval" --output_file="./regression_stats.csv"

## HMMs 

# ClinVar/gnomAD
python compute_binary_metrics.py --score_folders="../outputs/output_scores/ClinVar_HMM_postprocessed,../outputs/output_scores/gnomAD_HMM_postprocessed" \
--positive_label="Benign" --model_name="HMM" --output_file="./binary_stats.csv"

# DMS
python compute_DMS_regression_metrics.py --score_folder="../outputs/output_scores/DMS_HMM_postprocessed" \
--model_name="HMM" --output_file="./regression_stats.csv"
