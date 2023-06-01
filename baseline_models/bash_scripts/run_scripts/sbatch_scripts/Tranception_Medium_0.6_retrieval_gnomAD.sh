#!/bin/bash
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 1-00:00                        # Runtime in D-HH:MM format
#SBATCH -p gpu_quad   #,gpu_marks,gpu,gpu_requeue        # Partition to run in
# If on gpu_quad, use teslaV100s
# If on gpu_requeue, use teslaM40 or a100?
# If on gpu, any of them are fine (teslaV100, teslaM40, teslaK80) although K80 sometimes is too slow
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1,vram:48G
#SBATCH --qos=gpu_quad
#SBATCH --mem=32G                          # Memory total in MB (for all cores)

#SBATCH --mail-type=TIME_LIMIT_80,TIME_LIMIT,FAIL,ARRAY_TASKS
#SBATCH --mail-user="Daniel_Ritter@hms.harvard.edu"

#SBATCH -o slurm_files/slurm-%j.out                 # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH --job-name="tranception_large_scoring_2023_04_17_provean_clinvar_60R"

# Job array-specific
#SBATCH --output=slurm_files/slurm-lvn-%A_%3a-%x.out   # Nice tip: using %3a to pad to 3 characters (23 -> 023)
#SBATCH --error=slurm_files/slurm-lvn-%A_%3a-%x.err   # Optional: Redirect STDERR to its own file
#SBATCH --array=0-696:250
################################################################################

set -e # fail fully on first line failure (from Joost slurm_for_ml)

echo "hostname: $(hostname)"
echo "Running from: $(pwd)"
echo "GPU available: $(nvidia-smi)"
module load gcc/6.2.0
module load cuda/10.2 
module load miniconda3/4.10.3

source activate /n/groups/marks/software/anaconda_o2/envs/indels_env
# Ideally could get this just from slurm job numbers, but apparently, when using step sizes, it just gives back the greatest multiple of the step size, not the highest number of the range you enter
MAX_DMS_ID=696
DMS_RANGE_START=$SLURM_ARRAY_TASK_ID 
if (($SLURM_ARRAY_TASK_ID+$SLURM_ARRAY_TASK_STEP < $MAX_DMS_ID+1)); then
    DMS_RANGE_END="$(($SLURM_ARRAY_TASK_ID+$SLURM_ARRAY_TASK_STEP))"
else 
    DMS_RANGE_END="$(($MAX_DMS_ID+1))" 
fi 
REPO_ROOT="$(git rev-parse --show-toplevel)" 
PROTEIN_FOLDER="${REPO_ROOT}/processed_data/gnomAD/by_uniparc"
MAP_FILE="${REPO_ROOT}/processed_data/gnomAD/mapfiles/gnomad_filtered_common_indels_dedup_lt3aa_mapfile.csv"
MSA_FOLDER="${REPO_ROOT}/alignments/focus_column_only/gnomAD"
NAME="gnomAD_Tranception_Medium_0.6_retrieval"
OUTPUT_FOLDER="${REPO_ROOT}/outputs/output_scores/${NAME}"
MODEL_CHECKPOINT="../../model_checkpoints/Tranception_Medium"
TOKENIZER_PATH="../../model_checkpoints/tokenizers/Tranception/Basic_tokenizer"
BATCH_SIZE=1
CLUSTAL_LOCATION="${REPO_ROOT}/clustal/clustalo-1.2.4-Ubuntu"

/n/groups/marks/software/anaconda_o2/envs/indels_env/bin/python ../../scoring/score_tranception.py \
--dataset_reference_file=$MAP_FILE --dataset_folder=$PROTEIN_FOLDER --target_sequence_index_range_start=$DMS_RANGE_START \
--target_sequence_index_range_end=$DMS_RANGE_END \
--output_scores_folder=$OUTPUT_FOLDER --indel_mode \
--inference_time_retrieval --checkpoint=$MODEL_CHECKPOINT \
--clustal_omega_location=$CLUSTAL_LOCATION --MSA_folder=$MSA_FOLDER --batch_size_inference=$BATCH_SIZE \
--tokenizer_filepath=$TOKENIZER_PATH --num_workers=4