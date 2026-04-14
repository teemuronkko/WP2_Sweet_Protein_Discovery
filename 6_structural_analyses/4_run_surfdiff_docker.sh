#! /bin/bash
export INPUT_DATA_PATH="$1"
export OUTPUT_PATH="$2"
QUERY=$3
shift 3
SUBJECTS=("$@")   

mkdir -p $OUTPUT_PATH
docker run --cpus=10 --rm -v $OUTPUT_PATH:/workspace -v $INPUT_DATA_PATH:/input surfdiff compare --query-id $QUERY --subject-ids "${SUBJECTS[@]}" -p /input/pdb_dir/ -j 1