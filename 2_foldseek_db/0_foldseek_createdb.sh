#!/bin/bash -l
#SBATCH --partition=page
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=foldseek_createDB
#SBATCH --output=/groups/scienceai/tronkko/temp/logs/%j_out.txt
#SBATCH --error=/groups/scienceai/tronkko/temp/logs/%j_err.txt
#SBATCH --time=1-00:00:00
#SBATCH --gres=gpu:1

module purge
module load astro
module load cuda/12.8

# Export cuda environment variables
export CUDA_HOME=/software/astro/cuda/12.8/
export PATH=/software/astro/cuda/12.8/bin/${PATH:+:${PATH}}$
export CPATH=/software/astro/cuda/12.8/include:$CPATH
export LIBRARY_PATH=/software/astro/cuda/12.8/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=/software/astro/cuda/12.8/lib64/${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}$
export CUDA_LAUNCH_BLOCKING=1
export GGML_CUDA_LOG=1

# Detect node name
NODE="${SLURMD_NODENAME:-$HOSTNAME}"

echo "Running on node: $NODE"
if [[ "$NODE" == "node263" ]]; then
    echo "Activating env for node263"
    conda activate foldseek_env
    WDIR="/groups/scienceai/tronkko/foldseek/foldseek"
elif [[ "$NODE" == "node264" ]]; then
    echo "Activating env for node264"
    conda activate /lustre/hpc/scienceai/tronkko/foldseek_env_V100
    WDIR="/groups/scienceai/tronkko/foldseek/foldseek_V100/foldseek"
fi

cd $WDIR
export PATH="$WDIR/build/bin":$PATH

# Path to ProstT5 model weights
MODEL_WEIGHTS="/groups/scienceai/tronkko/foldseek/weights"
FOLDSEEK_DIR="/groups/scienceai/tronkko/WP2_SweetProteinDiscovery/2_foldseek_db"
TMP_DIR="/groups/scienceai/tronkko/foldseek/tmp_gpu"
mkdir -p $TMP_DIR

INPUT_DB=$1                      # Path to input database (full fasta file)

# check which GPU is available
echo "Checking available GPU..."
if command -v nvidia-smi &> /dev/null; then
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -n 1)
    echo "NVIDIA GPU detected: $GPU_NAME"
    # report gpu index/number (0, 1, etc)
    GPU_INDEX=$(nvidia-smi --query-gpu=index --format=csv,noheader | head -n 1)
    echo "Using GPU index: $GPU_INDEX"
    echo "nvidia-smi output:"
    nvidia-smi
elif command -v rocm-smi &> /dev/null; then
    GPU_NAME=$(rocm-smi --showproductname | grep 'Product Name' | awk -F ': ' '{print $2}' | head -n 1)
    echo "AMD GPU detected: $GPU_NAME"
    # report gpu index/number (0, 1, etc)
    GPU_INDEX=$(rocm-smi --showid | grep 'GPU ID' | awk -F ': ' '{print $2}' | head -n 1)
    echo "Using GPU index: $GPU_INDEX"
else
    echo "No compatible GPU found. Exiting."
    exit 1
fi

echo "CUDA VISIBLE DEVICES"
echo $CUDA_VISIBLE_DEVICES

# Build input paths
BATCH="${SLURM_ARRAY_TASK_ID}"   # Batch is based on array task id
DB_NAME=$(basename "$INPUT_DB")  # Parse database name without extension
DB_NAME="${DB_NAME%.*}"          
DB_SOURCE="${DB_NAME##*_}"       # Parse database source - OneKP or UniProt

# Build batch folder path
BATCH_FOLDER="/groups/scienceai/tronkko/WP2_SweetProteinDiscovery/1_viridiplantae_prot_db/3_database_batches/${DB_SOURCE}_batches"

# Build input and output paths for the batch 
DB_BATCH="${BATCH_FOLDER}/${DB_NAME}_batch_${BATCH}.fasta"
DB_OUTPUT_DIR="${FOLDSEEK_DIR}/${DB_NAME}/${DB_NAME}_batch_${BATCH}"
DB_OUTPUT="${DB_OUTPUT_DIR}/${DB_NAME}_batch_${BATCH}"

# Print input and outputs
echo -e "Input database:\t$INPUT_DB"
echo -e "Batch number:\t$BATCH"
echo -e "Database name:\t$DB_NAME"
echo -e "Database batch file:\t$DB_BATCH"
echo -e "Output folder:\t$DB_OUTPUT_DIR"
echo -e "Output path:\t$DB_OUTPUT"
mkdir -p $DB_OUTPUT_DIR

echo "Launching Foldseek DB creation using ProstT5 model."
start_time=$(date +%s)
foldseek createdb $DB_BATCH $DB_OUTPUT --prostt5-model $MODEL_WEIGHTS --gpu 1 --threads $SLURM_CPUS_PER_TASK --write-lookup 1
end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "Database created at $DB_OUTPUT"
echo "Runtime: ${runtime} seconds"

# Move logs
LOG_DIR="${FOLDSEEK_DIR}/${DB_NAME}/logs"
mkdir -p $LOG_DIR
mv /groups/scienceai/tronkko/temp/logs/${SLURM_JOB_ID}_out.txt ${LOG_DIR}/${DB_SOURCE}_createdb_batch_${BATCH}_${SLURM_JOB_ID}_out.txt
mv /groups/scienceai/tronkko/temp/logs/${SLURM_JOB_ID}_err.txt ${LOG_DIR}/${DB_SOURCE}_createdb_batch_${BATCH}_${SLURM_JOB_ID}_err.txt
