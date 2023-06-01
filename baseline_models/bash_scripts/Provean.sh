#!/bin/bash
# These determine the range of indexes to run. The indexes correspond to the rows of the mapping file for the gnomAD data. 
# Note: Provean does a complete BLAST search for each protein to generate an alignment, which it then uses to compute 
# it's final score. The BLAST search is very slow (often several hours per search), so this script is quite slow and 
# can be sped up by splitting the range of indexes to run into multiple scripts and running them in parallel. 
INDEX_RANGE_START=$2 
INDEX_RANGE_END=$3

NAME=$1
OUTPUT_SCORES_FOLDER="../../../outputs/output_scores/${NAME}"
OUTPUT_RAW_SCORES_FOLDER="../../../outputs/intermediate_outputs/${NAME}/raw_output_scores" 
FASTA_FOLDER="../../../outputs/intermediate_outputs/${NAME}/fasta_files"
SUPPORTING_SET_FOLDER="../../../outputs/intermediate_outputs/${NAME}/output_supporting_sets"
PROTEIN_FOLDER=$4
MAP_FILE=$5
NUM_THREADS=4 
HGVS_COLUMN=$6 
LABEL_COLUMN=$7
PROVEAN_SCRIPT_PATH="../../models/provean/provean.sh" 

python ../../scoring/score_provean.py \
--provean_script_path=$PROVEAN_SCRIPT_PATH --dataset_reference_file=$MAP_FILE --dataset_folder=$PROTEIN_FOLDER \
--fasta_folder=$FASTA_FOLDER --target_sequence_index_range_start=$INDEX_RANGE_START --target_sequence_index_range_end=$INDEX_RANGE_END --num_threads=$NUM_THREADS --keep_fasta --indel_mode \
--name_column=$HGVS_COLUMN --output_scores_folder=$OUTPUT_SCORES_FOLDER --output_raw_scores_folder=$OUTPUT_RAW_SCORES_FOLDER --supporting_set_folder=$SUPPORTING_SET_FOLDER --label_column=$LABEL_COLUMN