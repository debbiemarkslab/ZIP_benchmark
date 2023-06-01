#!/bin/bash
# NOTE: The model checkpoints are not kept in the repository due to their large size. Download them separately from the link in the README and place
# them in the model_checkpoints folder, or place them somewhere else and change the MODEL_CHECKPOINT path to point there. 
set -e # fail fully on first line failure (from Joost slurm_for_ml)
# These determine the range of indexes to run. The indexes correspond to the rows of the mapping file for the ClinVar data. 
INDEX_RANGE_START=$2 
INDEX_RANGE_END=$3


NAME=$1
REPO_ROOT="../../../"
OUTPUT_FOLDER="${REPO_ROOT}/outputs/output_scores/${NAME}"
MODEL_CHECKPOINT=$4
PROTEIN_FOLDER=$5
MAP_FILE=$6
LABEL_COLUMN=$7
TOKENIZER_PATH="../model_checkpoints/tokenizers/Tranception/Basic_tokenizer"
python ../../scoring/score_tranception.py \
--dataset_reference_file=$MAP_FILE --dataset_folder=$PROTEIN_FOLDER --target_sequence_index_range_start=$INDEX_RANGE_START \
--target_sequence_index_range_end=$INDEX_RANGE_END \
--output_scores_folder=$OUTPUT_FOLDER --indel_mode \
--checkpoint=$MODEL_CHECKPOINT --tokenizer_filepath=$TOKENIZER_PATH --label_column=$LABEL_COLUMN
