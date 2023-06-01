#! /bin/bash

# This script runs the preprocessing of gnomAD data for the indels benchmark.

# Replace this path with the path to your copy of the gnomAD database file
# Note that this version has additional headers stripped out to make parsing much faster. The 
# function that does this is strip_gnomad_headers in utils/gnomad_utils.py. Run that function
# on your original gnomAD database first to significantly speed up the preprocessing.
GNOMAD_DATABASE_FILE="/n/groups/marks/databases/gnomad/v2_exomes_GRCh38/latest_gnomad.exomes.r2.1.1.sites.liftover_grch38_AF_vep_only.vcf.bgz"

# python step1_process_indels.py --gnomad_database_file=$GNOMAD_DATABASE_FILE
# python step2_filter_indels.py 
# python step3_split_variants_and_create_mapfile.py --input_file="../../processed_data/gnomAD/gnomad_filtered_common_indels_dedup_lt3aa.csv" \
# --output_folder="../../processed_data/gnomAD/by_uniparc" --mapping_file="../../processed_data/gnomAD/mapfiles/gnomad_filtered_common_indels_dedup_lt3aa_mapfile.csv" \
# --alignment_folder="../../alignments/focus_column_only/gnomAD"

python step3_split_variants_and_create_mapfile.py --input_file="../../processed_data/gnomAD/gnomad_filtered_singleton_indels_dedup_lt3aa_random10.csv" \
--output_folder="../../processed_data/gnomAD/by_uniparc_singletons/" --mapping_file="../../processed_data/gnomAD/mapfiles/gnomad_filtered_singleton_indels_dedup_lt3aa_random10_mapfile.csv" \
--alignment_folder="../../alignments/focus_column_only/gnomAD"
