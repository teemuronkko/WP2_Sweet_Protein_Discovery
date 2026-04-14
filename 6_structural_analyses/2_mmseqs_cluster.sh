# Input and output paths
SWEET_PROTEINS=("Thaumatin-I_II")
for SWEET_PROTEIN in "${SWEET_PROTEINS[@]}"; do
	INPUT_DB=./1_tm_scores/merged_hits/${SWEET_PROTEIN}_merged_hits.fasta
	OUTPUT_DIR=./2_mmseqs/clustered_hits/${SWEET_PROTEIN}
	OUTPUT=${OUTPUT_DIR}/${SWEET_PROTEIN}
	TMP_DIR="./tmp"
	mkdir -p $OUTPUT_DIR
	mkdir -p $TMP_DIR
	echo "INPUT_DB: $INPUT_DB"
	echo "TMP_DIR: $TMP_DIR"

	# Launch mmseqs2 easy-cluster
	mmseqs easy-cluster "$INPUT_DB" $OUTPUT $TMP_DIR \
		-s 7.5 \
		-c 0.90 \
		--min-seq-id 0.90 \
		--threads 8 \
		-e 0.001 \
		--cov-mode 0 \
		--cluster-mode 2
done