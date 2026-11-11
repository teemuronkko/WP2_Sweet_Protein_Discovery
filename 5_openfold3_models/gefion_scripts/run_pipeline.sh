#!/bin/bash
# Usage: ./run_pipeline.sh -i <inputs> -o <output> [-t fasta|a3m] [-b 200] [-m 16] [-n 8]

# --- DEFAULTS ---
BATCH_SIZE=20
TYPE="a3m"
THROTTLE_MSA=16
THROTTLE_INF=320

# --- PARSE ARGS ---
while getopts "i:o:t:b:m:n:h" opt; do
  case $opt in
    i) INPUT_DIR=$(realpath "$OPTARG") ;;
    o) WORK_DIR=$(realpath "$OPTARG") ;;
    t) TYPE="$OPTARG" ;;
    b) BATCH_SIZE="$OPTARG" ;;
    m) THROTTLE_MSA="$OPTARG" ;;
    n) THROTTLE_INF="$OPTARG" ;;
    h) echo "Usage: $0 ..."; exit 0 ;;
    *) echo "Invalid option"; exit 1 ;;
  esac
done

if [[ -z "$INPUT_DIR" ]] || [[ -z "$WORK_DIR" ]]; then
    echo "Error: Input and Output required."
    exit 1
fi

SCRIPT_DIR="/dcai/projects/cu_0027/scripts"

echo "=================================="
echo "      OPENFOLD3 PIPELINE          "
echo "=================================="
echo "Mode:    $TYPE"
echo "Limits:  MSA=%${THROTTLE_MSA} | INF=%${THROTTLE_INF}"
echo "Dependency: Streaming (1-to-1)"
echo "=================================="

# --- STEP 1: SHARDING ---
echo "[1/3] Sharding inputs..."
mkdir -p "$WORK_DIR"/{inputs,msas,predictions,jsons,logs,benchmarks}

batch_num=1
count=0
mkdir -p "${WORK_DIR}/inputs/batch_${batch_num}"
mkdir -p "${WORK_DIR}/msas/batch_${batch_num}"

if [[ "$TYPE" == "fasta" ]]; then PATTERN="*.fasta"; else PATTERN="*.a3m"; fi

while IFS= read -r -d '' f; do
    if [[ "$TYPE" == "fasta" ]]; then
        ln -sf "$f" "${WORK_DIR}/inputs/batch_${batch_num}/"
    else
        ln -sf "$f" "${WORK_DIR}/msas/batch_${batch_num}/"
    fi

    ((count++))
    if (( count >= BATCH_SIZE )); then
        ((batch_num++))
        mkdir -p "${WORK_DIR}/inputs/batch_${batch_num}"
        mkdir -p "${WORK_DIR}/msas/batch_${batch_num}"
        count=0
    fi
done < <(find "$INPUT_DIR" -maxdepth 1 -name "$PATTERN" -print0)

# Cleanup empty last batch
if [ -z "$(ls -A ${WORK_DIR}/inputs/batch_${batch_num} 2>/dev/null)" ] && \
   [ -z "$(ls -A ${WORK_DIR}/msas/batch_${batch_num} 2>/dev/null)" ]; then
   rmdir "${WORK_DIR}/inputs/batch_${batch_num}" 2>/dev/null
   rmdir "${WORK_DIR}/msas/batch_${batch_num}" 2>/dev/null
   ((batch_num--))
fi

echo "      Created $batch_num batches."

# --- STEP 2: SUBMIT JOBS ---
echo "[2/3] Submitting Slurm Arrays..."

if [[ "$TYPE" == "fasta" ]]; then
    
    # 1. Submit MSA
    JID_MSA=$(sbatch --parsable --array=1-${batch_num}%${THROTTLE_MSA} "${SCRIPT_DIR}/stage1_msa.slurm" "$WORK_DIR")
    echo "      Submitted MSA Job ID: $JID_MSA"
    
    # 2. Submit Inference (Use AFTERCORR for 1-to-1 dependency)
    JID_INF=$(sbatch --parsable --dependency=aftercorr:${JID_MSA} --array=1-${batch_num}%${THROTTLE_INF} "${SCRIPT_DIR}/stage2_infer.slurm" "$WORK_DIR")
    echo "      Submitted INF Job ID: $JID_INF (Starts as each MSA batch finishes)"

elif [[ "$TYPE" == "a3m" ]]; then
    JID_INF=$(sbatch --parsable --array=1-${batch_num}%${THROTTLE_INF} "${SCRIPT_DIR}/stage2_infer.slurm" "$WORK_DIR")
    echo "      Submitted INF Job ID: $JID_INF"
fi

echo "Done."
