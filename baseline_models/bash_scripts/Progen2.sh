#!/bin/bash
INDEX_RANGE_START=$2 
INDEX_RANGE_END=$3

NAME=$1
MODEL_PATH=$4
PROTEIN_FOLDER=$5
MAP_FILE=$6
LABEL_COLUMN=$7
OUTPUT_SCORES_FOLDER="../../../outputs/output_scores/${NAME}"
TOKENIZER_PATH="../model_checkpoints/tokenizers/Progen2/tokenizer.json"
python ../../score_progen.py \
--Progen2_model_name_or_path=$MODEL_PATH --dataset_reference_file=$MAP_FILE --dataset_folder=$PROTEIN_FOLDER \
--output_scores_folder=$OUTPUT_SCORES_FOLDER --target_sequence_index_range_start=$INDEX_RANGE_START --target_sequence_index_range_end=$INDEX_RANGE_END \
--indel_mode \
--tokenizer_file=$TOKENIZER_PATH --label_column=$LABEL_COLUMN
