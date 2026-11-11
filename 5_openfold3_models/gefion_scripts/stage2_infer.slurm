#!/bin/bash
#SBATCH --account=cu_0027
#SBATCH --job-name=OF3_INFER
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

# Directories
MSA_DIR="${WORK_DIR}/msas/batch_${BATCH_ID}"
JSON_DIR="${WORK_DIR}/jsons/batch_${BATCH_ID}"
PRED_DIR="${WORK_DIR}/predictions/batch_${BATCH_ID}"
BENCH_FILE="${WORK_DIR}/benchmarks/batch_${BATCH_ID}_infer.csv"

# Write Header
echo "sample_name,step,duration_sec" > "$BENCH_FILE"

# --- CONFIGURATION ---
OF3_ROOT="/dcai/software/structure-prediction/openfold3"
CONTAINER="${OF3_ROOT}/containers/openfoldconsortium+openfold3+stable.sqsh"
JSON_SCRIPT="${OF3_ROOT}/scripts/create_json.py"
CONV_SCRIPT="${OF3_ROOT}/scripts/safe_convert.py"  # <--- NEW: Converter Script

# Mounts (Mirror Mounts)
MOUNTS="${OF3_ROOT}:${OF3_ROOT},${WORK_DIR}:${WORK_DIR}"

module load Miniconda3

echo "Batch: $BATCH_ID"

for a3m in "$MSA_DIR"/*.a3m; do
    stem=$(basename "$a3m" .a3m)
    
    # Skip if final PDB exists
    if [[ -f "${PRED_DIR}/${stem}/${stem}_model_1.pdb" ]] || \
       [[ -f "${PRED_DIR}/${stem}/output.pdb" ]]; then 
        echo "Skipping $stem (PDB exists)"
        continue
    fi

    echo "Processing $stem..."
    SAMPLE_START=$(date +%s)

    # --- 1. JSON GEN ---
    TEMP_FASTA="${JSON_DIR}/${stem}.temp.fasta"
    head -n 2 "$a3m" > "$TEMP_FASTA"
    
    MSA_PATH="${MSA_DIR}/${stem}.a3m"
    JSON_PATH="${JSON_DIR}/${stem}.json"
    OUTPUT_PATH="${PRED_DIR}/${stem}"
    CKPT_PATH="${OF3_ROOT}/params/of3-p2-155k.pt"

    python3 "$JSON_SCRIPT" \
        "$TEMP_FASTA" "$MSA_PATH" "$JSON_PATH" \
        "$OF3_ROOT" "$OF3_ROOT"

    rm "$TEMP_FASTA"

    # --- 2. INFERENCE (Generates .cif) ---
    INFER_START=$(date +%s)
    
    srun --overlap --container-image="$CONTAINER" --container-mounts="$MOUNTS" \
        run_openfold predict \
        --query_json "$JSON_PATH" \
        --output_dir "$OUTPUT_PATH" \
        --inference_ckpt_path "$CKPT_PATH" \
        --use_msa_server False \
        --use_templates False

    INFER_END=$(date +%s)
    echo "${stem},inference,$((INFER_END - INFER_START))" >> "$BENCH_FILE"

    # --- 3. CONVERSION (CIF -> PDB) ---
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
