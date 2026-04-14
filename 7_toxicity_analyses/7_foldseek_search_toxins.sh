TMP_DIR="./tmp"
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
    --min-seq_id 0.5 \
	--threads 8 \
	-e 1e-7 \
	--alignment-type 2 \
	-c 0.80 \
	--cov-mode 1 \
	--format-output "query,target,evalue,bits,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,qseq,tseq,q3di,t3di,qheader,theader,qaln,taln,q3dialn,t3dialn,qcov,tcov" \
	--format-mode 4

# Move log files
LOG_PATH="/groups/scienceai/tronkko/temp/logs"
mv "${LOG_PATH}/${SLURM_JOB_ID}_out.txt" ${OUTPUT_DIR}/logs/foldseek_search_${DB_NAME}_${SLURM_JOB_ID}_out.txt
mv "${LOG_PATH}/${SLURM_JOB_ID}_err.txt" ${OUTPUT_DIR}/logs/foldseek_search_${DB_NAME}_${SLURM_JOB_ID}_err.txt
