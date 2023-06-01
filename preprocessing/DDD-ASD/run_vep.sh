#!/bin/bash

REFDIR=/n/groups/marks/projects/indels_human/data/vep_cache
# ANNOTATION_FILE="/n/groups/marks/projects/indels_human/data/iceland_denovo_data/merged_vcfs.vcf"
ANNOTATION_FILE="./ASD_annotated.vcf"
# OUTPUT_FILE="/n/groups/marks/projects/indels_human/data/iceland_denovo_data/merged_vcfs_vep_annotated.tsv"
OUTPUT_FILE="./ASD_vep_annotated.tsv"
CACHE_DIR="/n/groups/marks/projects/indels_human/data/vep_cache"

/n/groups/marks/projects/indels_human/software/ensembl-vep/vep -i $ANNOTATION_FILE -o $OUTPUT_FILE \
    --tab -v --offline --refseq --fork 8 --force_overwrite --buffer_size 20000 --format vcf \
    --no_stats \
    --assembly GRCh38 --use_transcript_ref \
    --fasta $REFDIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --hgvs  \
    --dir_cache=$CACHE_DIR --protein --biotype

ANNOTATION_FILE="./ASD_replicate_annotated.vcf" 
OUTPUT_FILE="./ASD_replicate_vep_annotated.tsv"

/n/groups/marks/projects/indels_human/software/ensembl-vep/vep -i $ANNOTATION_FILE -o $OUTPUT_FILE \
    --tab -v --offline --refseq --fork 8 --force_overwrite --buffer_size 20000 --format vcf \
    --no_stats \
    --assembly GRCh38 --use_transcript_ref \
    --fasta $REFDIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --hgvs  \
    --dir_cache=$CACHE_DIR --protein --biotype
# ANNOTATION_FILE="./DDD.vcf" 
# OUTPUT_FILE="./DDD_vep_annotated.tsv"

# /n/groups/marks/projects/indels_human/software/ensembl-vep/vep -i $ANNOTATION_FILE -o $OUTPUT_FILE \
    # --tab -v --offline --refseq --fork 8 --force_overwrite --buffer_size 20000 --format vcf \
    # --no_stats \
    # --assembly GRCh37 --use_transcript_ref \
    # --fasta $REFDIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    # --hgvs \
    # --dir_cache=$CACHE_DIR --protein --biotype
