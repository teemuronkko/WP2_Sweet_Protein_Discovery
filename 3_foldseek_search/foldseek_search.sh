#!/bin/bash -l
#SBATCH --partition=page
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --job-name=foldseek_search
#SBATCH --output=/groups/scienceai/tronkko/temp/logs/%j_out.txt
#SBATCH --error=/groups/scienceai/tronkko/temp/logs/%j_err.txt
#SBATCH --time=1-00:00:00
export CUDA_VISIBLE_DEVICES=""
USER=$(whoami)
if [[ "$USER" == "pqh443" ]]; then
    REPO_DIR="/Users/pqh443/Documents/Git_projects/WP2_SweetProteinDiscovery"
elif [[ "$USER" == "tronkko" ]]; then
    REPO_DIR="/groups/scienceai/tronkko/WP2_SweetProteinDiscovery"
    conda activate foldseek_env
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
else
    exit 0
fi

# Set foldseek paths and activate conda env
FOLDSEEK_DIR="/groups/scienceai/tronkko/foldseek"
MODEL_WEIGHTS="$FOLDSEEK_DIR/weights"
TMP_DIR="$FOLDSEEK_DIR/tmp"
mkdir -p $MODEL_WEIGHTS
mkdir -p $TMP_DIR

# Input and output paths
INPUT_DB=$1
TARGET_DB=$2
OUTPUT=$3
DB_NAME=$(basename "$TARGET_DB")
DB_NAME="${DB_NAME%.*}"
OUTPUT_DIR=$(dirname "$OUTPUT")
echo "Input DB: $INPUT_DB"
echo "TARGET_DB: $TARGET_DB"
echo "OUTPUT: $OUTPUT"
echo "TMP_DIR: $TMP_DIR"
mkdir -p ${OUTPUT_DIR}/logs/

# Launch foldseek easy-search
foldseek easy-search "$INPUT_DB" "$TARGET_DB" $OUTPUT $TMP_DIR \
	--exhaustive-search 1 \
	--min-ungapped-score 30 \
	--threads $SLURM_CPUS_PER_TASK \
	-e 0.05 \
	--alignment-type 2 \
	-c 0.50 \
	--cov-mode 2 \
	--format-output "query,target,evalue,bits,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,qseq,tseq,q3di,t3di,qheader,theader,qaln,taln,q3dialn,t3dialn,qcov,tcov" \
	--format-mode 4

# Move log files
LOG_PATH="/groups/scienceai/tronkko/temp/logs"
mv "${LOG_PATH}/${SLURM_JOB_ID}_out.txt" ${OUTPUT_DIR}/logs/foldseek_search_${DB_NAME}_${SLURM_JOB_ID}_out.txt
mv "${LOG_PATH}/${SLURM_JOB_ID}_err.txt" ${OUTPUT_DIR}/logs/foldseek_search_${DB_NAME}_${SLURM_JOB_ID}_err.txt
