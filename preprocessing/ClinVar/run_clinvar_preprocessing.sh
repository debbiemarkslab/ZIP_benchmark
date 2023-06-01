#!/bin/bash 

ROOT_DIR="../.."

CLINVAR_QUERY_FILE="${ROOT_DIR}/data/ClinVar/raw/clinvar_result_inframe_indels_20230206.txt"
GRCH38_PROTEIN_FILE="${ROOT_DIR}/data/ClinVar/raw/GRCh38_latest_protein.faa" 
REFSEQ_SEQUENCE_MAPPING="${ROOT_DIR}/data/mappings/refseq_mapping_gff.tsv" 
REFSEQ_UNIPROT_MAPPING="${ROOT_DIR}/data/reference_files/gene_refseq_uniprotkb_collab.gz"
OUTPUT_FOLDER="${ROOT_DIR}/processed_data/ClinVar"
ALIGNMENT_FOLDER="${ROOT_DIR}/alignments/focus_column_only/ClinVar"
# python step1_process_clinvar_indels.py --input_file=$CLINVAR_QUERY_FILE \
# --refseq_sequence_file=$REFSEQ_SEQUENCE_MAPPING --uniprot_refseq_mapping_file=$REFSEQ_UNIPROT_MAPPING \
# --grch38_protein_file=$GRCH38_PROTEIN_FILE --output_folder=$OUTPUT_FOLDER

PER_PROTEIN_FOLDER="${ROOT_DIR}/processed_data/ClinVar/by_refseq"
MAPFILE="${ROOT_DIR}/processed_data/ClinVar/mapfiles/ClinVar_mapping.csv"
python step2_split_variants_create_mapfile.py --input_file="$OUTPUT_FOLDER/clinvar_filtered_significant_pathogenic_variants.csv" \
--alignment_folder=$ALIGNMENT_FOLDER --output_folder=$PER_PROTEIN_FOLDER --mapping_file=$MAPFILE