#!/bin/bash 

# Change ranges to run portions of the clinvar variants 
# and comment out methods below to only run a subset of the baselines 
INDEX_RANGE_START=0 
INDEX_RANGE_END=7 
REPO_ROOT="../../.." 
PROTEIN_FOLDER="${REPO_ROOT}/processed_data/DMS"
MAP_FILE="${REPO_ROOT}/processed_data/DMS/mapfiles/DMS_mapping.csv"
MSA_FOLDER="${REPO_ROOT}/alignments/focus_column_only/DMS"
FULL_ALIGNMENT_FOLDER="${REPO_ROOT}/alignments/full/DMS" 
VARIANT_CSV="${REPO_ROOT}/processed_data/"

# Provean 
# bash ../Provean.sh "DMS_Provean" $INDEX_RANGE_START $INDEX_RANGE_END $PROTEIN_FOLDER $MAP_FILE "protein_variant" 

# Progen2 XLarge 
# bash ../Progen2.sh "DMS_Progen2_XLarge" $INDEX_RANGE_START $INDEX_RANGE_END "../../model_checkpoints/Progen2_XLarge" $PROTEIN_FOLDER $MAP_FILE 

# Progen2 Medium
# bash ../Progen2.sh "DMS_Progen2_Medium" $INDEX_RANGE_START $INDEX_RANGE_END "../../model_checkpoints/Progen2_Medium" $PROTEIN_FOLDER $MAP_FILE

# HMM 
bash ../HMM.sh "DMS_HMM" $INDEX_RANGE_START $INDEX_RANGE_END $PROTEIN_FOLDER $MAP_FILE $FULL_ALIGNMENT_FOLDER "label" 

# Tranception Medium with 0.6 retrieval 
# bash ../Tranception_0.6_retrieval.sh "DMS_Tranception_Medium_0.6_retrieval" $INDEX_RANGE_START $INDEX_RANGE_END \
# "../model_checkpoints/Tranception_Medium" $PROTEIN_FOLDER $MAP_FILE $MSA_FOLDER

# Tranception Large with 0.6 retrieval 
# bash ../Tranception_0.6_retrieval.sh "DMS_Tranception_Large_0.6_retrieval" $INDEX_RANGE_START $INDEX_RANGE_END \
# "../model_checkpoints/Tranception_Large" $PROTEIN_FOLDER $MAP_FILE $MSA_FOLDER

# Tranception Medium 
# bash ../Tranception.sh "DMS_Tranception_Medium" $INDEX_RANGE_START $INDEX_RANGE_END \
# "../model_checkpoints/Tranception_Medium" $PROTEIN_FOLDER $MAP_FILE

# Tranception Large 
# bash ../Tranception.sh "DMS_Tranception_Large" $INDEX_RANGE_START $INDEX_RANGE_END \
# "../model_checkpoints/Tranception_Large" $PROTEIN_FOLDER $MAP_FILE