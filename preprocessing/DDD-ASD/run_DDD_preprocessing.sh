#! /bin/bash

# This script runs the preprocessing of DDD data for the indels benchmark.


# DDD_RAW_FILE="/n/groups/marks/projects/Poly/data/31k_trios/41586_2020_2832_MOESM3_ESM.txt"
DDD_RAW_FILE="/n/groups/marks/projects/Poly/data/31k_trios/20221009_ddd_annotated_all.tsv"
ASD_RAW_FILE="/n/groups/marks/databases/sfari/data/denovo/zhou2022.cleaned.vep.tsv"
ASD_REPLICATE_RAW_FILE="/n/groups/marks/databases/sfari/data/denovo/zhou2022_spark_replicate.tsv"
python convert_SPARK_to_vcf.py --input_tsv_file=$ASD_RAW_FILE --output_file="ASD_annotated.vcf"
python convert_SPARK_to_vcf.py --input_tsv_file=$ASD_REPLICATE_RAW_FILE --output_file="ASD_replicate_annotated.vcf"
# INDEX_FILE="../../data/mappings/ddd_asd_annotation_file.tsv"
# COMBINED_FILE="./ASD-DDD_annotated.tsv"
# REFSEQ_REFERENCE_FOLDER="/n/groups/marks/databases/ukbiobank/users/rose/data/mapping/refseq_mapping_gff_updated"
# CANONICAL_ISOFORM_FILE="/n/groups/marks/projects/indels_human/data/ClinVar/raw/RefSeq_canonical_isoforms_list.txt"
# ENTREZ_EMAIL="Daniel_Ritter@hms.harvard.edu"
# MSA_FOLDER="../../alignments/focus_column_only/DDD" 
# COMBINED_FILE="/n/groups/marks/projects/Poly/data/31k_trios/20221108_ddd_asd_annotated_all.tsv" 
# python annotate_raw_DDD_and_ASD_files.py --raw_DDD_file=$DDD_RAW_FILE --raw_ASD_file=$ASD_RAW_FILE --raw_ASD_replicate_file=$ASD_REPLICATE_RAW_FILE --index_file=$INDEX_FILE --output_file=$COMBINED_FILE
# python step1_process_DDD_ASD_indels.py --input_file=$COMBINED_FILE --refseq_reference_folder=$REFSEQ_REFERENCE_FOLDER --entrez_email=$ENTREZ_EMAIL --canonical_isoform_file=$CANONICAL_ISOFORM_FILE
# python step2_split_variants_and_create_mapfile.py --alignment_folder=$MSA_FOLDER 