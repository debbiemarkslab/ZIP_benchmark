#!/bin/bash 

INDEX_RANGE_START=0 
INDEX_RANGE_END=697 
REPO_ROOT="../../.." 
PROTEIN_FOLDER="${REPO_ROOT}/processed_data/gnomAD/by_uniparc"
MAP_FILE="${REPO_ROOT}/processed_data/gnomAD/mapfiles/gnomad_filtered_common_indels_dedup_lt3aa_mapfile.csv"
MSA_FOLDER="${REPO_ROOT}/alignments/focus_column_only/gnomAD"
FULL_ALIGNMENT_FOLDER="${REPO_ROOT}/alignments/full/gnomAD" 
VARIANT_CSV="${REPO_ROOT}/processed_data/gnomAD/gnomad_filtered_common_indels_dedup_lt3aa.csv"
# Provean 
# bash ../Provean.sh "gnomAD_Provean" $INDEX_RANGE_START $INDEX_RANGE_END $PROTEIN_FOLDER $MAP_FILE "protein_variant" 

# Progen2 XLarge 
# bash ../Progen2.sh "gnomAD_Progen2_XLarge" $INDEX_RANGE_START $INDEX_RANGE_END "../../model_checkpoints/Progen2_XLarge" $PROTEIN_FOLDER $MAP_FILE 

# Progen2 Medium
# bash ../Progen2.sh "gnomAD_Progen2_Medium" $INDEX_RANGE_START $INDEX_RANGE_END "../../model_checkpoints/Progen2_Medium" $PROTEIN_FOLDER $MAP_FILE

# HMM 
bash ../HMM.sh "gnomAD_HMM" $INDEX_RANGE_START $INDEX_RANGE_END $PROTEIN_FOLDER $MAP_FILE $FULL_ALIGNMENT_FOLDER "label" 

# Tranception Medium with 0.6 retrieval 
# bash ../Tranception_0.6_retrieval.sh "gnomAD_Tranception_Medium_0.6_retrieval" $INDEX_RANGE_START $INDEX_RANGE_END \
# "../model_checkpoints/Tranception_Medium" $PROTEIN_FOLDER $MAP_FILE $MSA_FOLDER

# Tranception Large with 0.6 retrieval 
# bash ../Tranception_0.6_retrieval.sh "gnomAD_Tranception_Large_0.6_retrieval" $INDEX_RANGE_START $INDEX_RANGE_END \
# "../model_checkpoints/Tranception_Large" $PROTEIN_FOLDER $MAP_FILE $MSA_FOLDER

# Tranception Medium 
# bash ../Tranception.sh "gnomAD_Tranception_Medium" $INDEX_RANGE_START $INDEX_RANGE_END \
# "../model_checkpoints/Tranception_Medium" $PROTEIN_FOLDER $MAP_FILE

# Tranception Large 
# bash ../Tranception.sh "gnomAD_Tranception_Large" $INDEX_RANGE_START $INDEX_RANGE_END \
# "../model_checkpoints/Tranception_Large" $PROTEIN_FOLDER $MAP_FILE