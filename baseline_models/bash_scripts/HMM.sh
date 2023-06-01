#!/bin/bash

INDEX_RANGE_START=$2 
INDEX_RANGE_END=$3

NAME=$1
REPO_ROOT="../../.."
OUTPUT_FOLDER="${REPO_ROOT}/outputs/output_scores/${NAME}"
INTERMEDIATE_OUTPUTS_FOLDER="${REPO_ROOT}/outputs/intermediate_outputs/${NAME}"
PROTEIN_FOLDER=$4
MAP_FILE=$5
MSA_FOLDER=$6
LABEL_COLUMN=$7
HMMER_PATH="/n/groups/marks/software/jackhmmer/hmmer-3.3.1"
REFORMAT_SCRIPT_PATH="/n/groups/marks/software/hhsuite/hh-suite3/scripts"
python ../../scoring/score_hmm.py \
--dataset_reference_file=$MAP_FILE --dataset_folder=$PROTEIN_FOLDER --target_sequence_index_range_start=$INDEX_RANGE_START \
--target_sequence_index_range_end=$INDEX_RANGE_END \
--hmmer_path=$HMMER_PATH --reformat_script_path=$REFORMAT_SCRIPT_PATH \
--output_scores_folder=$OUTPUT_FOLDER --intermediate_outputs_folder=$INTERMEDIATE_OUTPUTS_FOLDER \
--label_column=$LABEL_COLUMN

# reformat_script_path="/n/groups/marks/software/hhsuite/hh-suite3/scripts"
# hmmer_path="/n/groups/marks/software/jackhmmer/hmmer-3.3.1" 
# name=$1 
# map_file=$2
# ali_path=$3
# variant_folder=$4

# hmm_path=../../../outputs/intermediate_outputs/$name/hmm
# fa_path=../../../outputs/intermediate_outputs/$name/fa
# csv_path=../../../outputs/output_scores/$name

# mkdir -p $hmm_path
# mkdir -p $fa_path
# mkdir -p $csv_path

# # convert a3ms to a2ms if needed, including insertions in a2m output 
# # skip a2ms that are already present and have the correct number of sequences
# # also skips entirely if there are not a2ms present in the file 
# while read -r a3mfile; do
#     a2mfile="${a3mfile/.a3m/.a2m}"
#     hmmfile="${hmm_path}/${a3mfile/.a3m/.hmm}"
#     if [[ -f $a2mfile && -s $a2mfile ]] && [[ ! -f $hmmfile && ! -s $hmmfile ]] && (( "$(grep -c '>' $a3mfile)" == "$(grep -c '>' $a2mfile)" )); then
#         continue
#     fi
#     $reformat_script_path/reformat.pl a3m a2m $a3mfile $a2mfile
# done <<< "$(ls $ali_path/*.a3m)"

# # build hmms if needed
# a2m_count=ls $ali_path/*.a2m 2>/dev/null | wc -l
# if [[ $a2m_count -ne 0 ]]; then
#     while read -r alifile; do
#         ali_name="$(basename $alifile)"
#         hmmfile="$hmm_path/$ali_name"
#         hmmfile="${hmmfile/.a2m/.hmm}"
#         if [[ -f $hmmfile && -s $hmmfile ]]; then
#             continue
#         fi
#         echo "hmmbuild $hmmfile $alifile"
#         $hmmer_path/bin/hmmbuild --amino "$hmmfile" "$alifile"
#     done <<< "$(ls $ali_path/*.a2m)"
# fi

# fapath=$fa_path alikey=$ali_key variant_folder=$variant_folder mapfile=$map_file python -c '
# import pandas as pd
# import os
# from bio.seqio.fastaio import simplefastaparser
# map_df = pd.read_csv(os.environ["mapfile"])
# for seq in map_df["target_seq"].unique():
#     seq_name = map_df[map_df["target_seq"] == seq]["target_sequence_id"].values[0]
#     if not os.path.exists(os.environ["variant_folder"] + os.sep + f"{seq_name}.csv"):
#         print(f"no variant file present for {seq_name}")
#         continue
#     with open(os.environ["fapath"] + os.sep + f"{seq_name}.fa", "w+") as f:
#         with open(os.environ["variant_folder"] + os.sep + f"{seq_name}.csv", "r") as v:
#             subset = pd.read_csv(v)
#         f.write(f">wt\n{seq}\n")
#         for row in subset.itertuples():
#             f.write(f">{row.protein_variant}\n{row.protein_sequence_mutated}\n")
# '

# # run fwdback algorithm
# while read -r hmmfile; do
#     hmm_name="$(basename $hmmfile)"
#     seqfile="$fa_path/${hmm_name/.hmm/.fa}"
#     csvfile="$csv_path/${hmm_name/.hmm/.csv}"
#     echo "$hmmer_path/src/generic_fwdback_example $hmmfile $seqfile > $csvfile"
#     $hmmer_path/src/generic_fwdback_example $hmmfile $seqfile > $csvfile || echo "error in generic_fwdback_example $hmmfile $seqfile > $csvfile"
# done <<< "$(ls $hmm_path/*.hmm)"
