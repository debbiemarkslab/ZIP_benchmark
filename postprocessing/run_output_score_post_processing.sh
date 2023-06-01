#!/bin/bash

# MODEL_NAMES=( 'Provean' 'Progen2_XLarge' 'Progen2_Medium' 'HMM' \
# 'Tranception_Medium' 'Tranception_Medium_0.6_retrieval' 'Tranception_Large' 'Tranception_Large_0.6_retrieval' )

# for MODEL_NAME in "${MODEL_NAMES[@]}"
# do
#     echo "Combining ClinVar and gnomAD score folders for $MODEL_NAME" 
#     python3 combine_score_folders.py --score_folders="../outputs/output_scores/ClinVar_$MODEL_NAME,../outputs/output_scores/gnomAD_$MODEL_NAME" \
#     --output_folder="../outputs/combined_output_scores/ClinVar_gnomAD_$MODEL_NAME"
    
#     # Note: haven't added in DMS scores yet but should follow this name convention 
#     # echo "Combining DMS scores for $MODEL_NAME" 
#     # python3 combine_score_folders.py --score_folders="../outputs/output_scores/DMS_$MODEL_NAME" \
#     # --output_folder="../outputs/combined_output_scores/DMS_$MODEL_NAME"
    
#     # Note: may want DDD here too if we split (or can keep it combined in preprocessing)
# done

# HMM_SCORE_FOLDER="../outputs/output_scores/gnomAD_HMM"
# OUTPUT_FOLDER="../outputs/output_scores/gnomAD_HMM_postprocessed"
# DATASET_FOLDER="../processed_data/gnomAD/by_uniparc" 
HMM_SCORE_FOLDER="../outputs/output_scores/ClinVar_HMM"
OUTPUT_FOLDER="../outputs/output_scores/ClinVar_HMM_postprocessed"
DATASET_FOLDER="../processed_data/ClinVar/by_refseq"
LABEL_COLUMN="ClinSigSimple"
python add_hmm_wt_ratios_and_labels.py --hmm_score_folder=$HMM_SCORE_FOLDER --output_folder=$OUTPUT_FOLDER \
--dataset_folder=$DATASET_FOLDER --label_column=$LABEL_COLUMN