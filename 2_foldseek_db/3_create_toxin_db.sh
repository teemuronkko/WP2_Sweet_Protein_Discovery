#!/bin/bash -l
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
foldseek createdb $INPUT_DB $DB_OUTPUT