# Input and output paths
SWEET_PROTEINS=("Honey_truffle_sweetener" "MNEI_1IV7")
for SWEET_PROTEIN in "${SWEET_PROTEINS[@]}"; do
	INPUT_DB=./0_${SWEET_PROTEIN}_hits.fasta

	# Launch mmseqs2 easy-cluster
	for id in 0.3 0.4 0.5 0.6 0.7; do
		OUTPUT_DIR=./2_mmseqs_clustered_hits/${SWEET_PROTEIN}
		OUTPUT=${OUTPUT_DIR}/${SWEET_PROTEIN}_${id}
		TMP_DIR="${OUTPUT_DIR}/tmp"
		mkdir -p $TMP_DIR
		mkdir -p $OUTPUT_DIR
		echo "INPUT_DB: $INPUT_DB"
		echo "TMP_DIR: $TMP_DIR"
    	mmseqs easy-cluster "$INPUT_DB" $OUTPUT $TMP_DIR -s 7.5 -c 0.9 --cov-mode 0 --cluster-mode 2 --min-seq-id ${id} --threads 8 -e 0.001
	done
done