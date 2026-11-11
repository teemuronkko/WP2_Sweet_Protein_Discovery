#!/bin/bash
#SBATCH --account=cu_0027
#SBATCH --job-name=OF3_INFER_JSON
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --output=/dev/null

WORK_DIR="$1"
BATCH_ID=$SLURM_ARRAY_TASK_ID

LOG_FILE="${WORK_DIR}/logs/infer_batch_${BATCH_ID}.log"
exec > "$LOG_FILE" 2>&1

echo -e $WORK_DIR
echo "above is workdir"

# Directories
JSON_DIR="${WORK_DIR}/jsons/batch_${BATCH_ID}"
PRED_DIR="${WORK_DIR}/predictions/batch_${BATCH_ID}"
BENCH_FILE="${WORK_DIR}/benchmarks/batch_${BATCH_ID}_infer.csv"

# Write Header
echo "sample_name,step,duration_sec" > "$BENCH_FILE"

# --- CONFIGURATION ---
OF3_ROOT="/dcai/software/structure-prediction/openfold3"
CONTAINER="${OF3_ROOT}/containers/openfoldconsortium+openfold3+stable.sqsh"
CONV_SCRIPT="${OF3_ROOT}/scripts/safe_convert.py"

# Mounts (Mirror Mounts)
MOUNTS="${OF3_ROOT}:${OF3_ROOT},${WORK_DIR}:${WORK_DIR}"

module load Miniconda3

echo "Batch: $BATCH_ID"

for json in "$JSON_DIR"/*.json; do
    stem=$(basename "$json" .json)
    
    # Skip if final PDB exists
    if [[ -f "${PRED_DIR}/${stem}/${stem}_model_1.pdb" ]] || \
       [[ -f "${PRED_DIR}/${stem}/output.pdb" ]]; then 
        echo "Skipping $stem (PDB exists)"
        continue
    fi

    echo "Processing $stem..."
    SAMPLE_START=$(date +%s)

    JSON_PATH="${JSON_DIR}/${stem}.json"
    OUTPUT_PATH="${PRED_DIR}/${stem}"
    CKPT_PATH="${OF3_ROOT}/params/of3-p2-155k.pt"

    # --- 1. INFERENCE (Generates .cif) ---
    INFER_START=$(date +%s)
    
    srun --overlap --container-image="$CONTAINER" --container-mounts="$MOUNTS" \
        run_openfold predict \
        --query_json "$JSON_PATH" \
        --output_dir "$OUTPUT_PATH" \
        --inference_ckpt_path "$CKPT_PATH" \
        --use_msa_server False \
        --use_templates False 
    #   --runner_yaml "$WORK_DIR/runner.yml"

    INFER_END=$(date +%s)
    echo "${stem},inference,$((INFER_END - INFER_START))" >> "$BENCH_FILE"

    # --- 2. CONVERSION (CIF -> PDB) ---
    # We run this inside the container to ensure we have the right python libraries
    CONV_START=$(date +%s)
    
    echo "Converting $stem to PDB..."
    srun --overlap --container-image="$CONTAINER" --container-mounts="$MOUNTS" \
         python3 "$CONV_SCRIPT" "$OUTPUT_PATH"

    CONV_END=$(date +%s)
    echo "${stem},conversion,$((CONV_END - CONV_START))" >> "$BENCH_FILE"

    # --- TOTAL TIME ---
    echo "${stem},total_loop,$((CONV_END - SAMPLE_START))" >> "$BENCH_FILE"

done
