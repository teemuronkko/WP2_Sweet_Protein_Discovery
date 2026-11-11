#!/bin/bash -l
#SBATCH --partition=page
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=foldseek_createDB
#SBATCH --output=~/temp/logs/%j_out.txt
#SBATCH --error=~/temp/logs/%j_err.txt
#SBATCH --time=1-00:00:00
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
    WDIR="~/foldseek/foldseek"
elif [[ "$NODE" == "node264" ]]; then
    echo "Activating env for node264"
    conda activate ~/foldseek_env_V100
    WDIR="~/foldseek/foldseek_V100/foldseek"
fi

cd $WDIR
export PATH="$WDIR/build/bin":$PATH

# Path to ProstT5 model weights
MODEL_WEIGHTS="~/foldseek/weights"
FOLDSEEK_DIR="~/WP2_Sweet_Protein_Discovery/foldseek_db"
TMP_DIR="~/foldseek/tmp"
mkdir -p $TMP_DIR

INPUT_DB=$1                      # Path to input database
DB_NAME=$(basename "$INPUT_DB")  # Parse database name without extension
DB_NAME="${DB_NAME%.*}"          
DB_OUTPUT="${FOLDSEEK_DIR}/${DB_NAME}/${DB_NAME}"

# Print input and outputs
echo -e "Input database:\t$INPUT_DB"
echo -e "Database name:\t$DB_NAME"
echo -e "Output path:\t$DB_OUTPUT"
mkdir -p $DB_OUTPUT_DIR

echo "Launching Foldseek DB creation."
start_time=$(date +%s)
foldseek createdb $INPUT_DB $DB_OUTPUT
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Database created at $DB_OUTPUT"
echo "Runtime: ${runtime} seconds"

# Move logs
LOG_DIR="${FOLDSEEK_DIR}/${DB_NAME}/logs"
mkdir -p $LOG_DIR
mv ~/temp/logs/${SLURM_JOB_ID}_out.txt ${LOG_DIR}/${DB_SOURCE}_createdb_query_${SLURM_JOB_ID}_out.txt
mv ~/temp/logs/${SLURM_JOB_ID}_err.txt ${LOG_DIR}/${DB_SOURCE}_createdb_query_${SLURM_JOB_ID}_err.txt
