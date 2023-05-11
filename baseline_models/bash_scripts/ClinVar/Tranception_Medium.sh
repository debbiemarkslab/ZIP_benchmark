#!/bin/bash
# NOTE: The model checkpoints are not kept in the repository due to their large size. Download them separately from the link in the README and place
# them in the model_checkpoints folder, or place them somewhere else and change the MODEL_CHECKPOINT path to point there. 
set -e # fail fully on first line failure (from Joost slurm_for_ml)
# These determine the range of indexes to run. The indexes correspond to the rows of the mapping file for the ClinVar data. 
INDEX_RANGE_START=0 
INDEX_RANGE_END=847


NAME="ClinVar_Tranception_Large"
REPO_ROOT="../../../"
OUTPUT_FOLDER="${REPO_ROOT}/outputs/output_scores/${NAME}"
PROTEIN_FOLDER="${REPO_ROOT}/processed_data/ClinVar/by_refseq/"
MAP_FILE="${REPO_ROOT}/processed_data/ClinVar/2023-05-05-tranception_mapping_by_refseq_id.csv"
MODEL_CHECKPOINT="../model_checkpoints/Tranception_Medium"
TOKENIZER_PATH="../model_checkpoints/tokenizers/Tranception/Basic_tokenizer"
python ../../scoring/score_tranception.py \
--dataset_reference_file=$MAP_FILE --dataset_folder=$PROTEIN_FOLDER --target_sequence_reference_range_start=$INDEX_RANGE_START \
--target_sequence_reference_range_end=$INDEX_RANGE_END \
--output_scores_folder=$OUTPUT_FOLDER --indel_mode \
--checkpoint=$MODEL_CHECKPOINT --tokenizer_filepath=$TOKENIZER_PATH
