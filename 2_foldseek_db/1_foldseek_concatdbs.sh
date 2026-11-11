#!/bin/bash -l
#SBATCH --partition=page
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
#SBATCH --job-name=foldseek_concatDB
#SBATCH --output=~/temp/logs/%j_out.txt
#SBATCH --error=~/temp/logs/%j_err.txt
#SBATCH --time=1:00:00
USER=$(whoami)
if [[ "$USER" == "pqh443" ]]; then
    repo_path="~/Documents/Git_projects/WP2_Sweet_Protein_Discovery"
elif [[ "$USER" == "tronkko" ]]; then
    repo_path="~/WP2_Sweet_Protein_Discovery"
    conda activate foldseek_env
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
else
    exit 0
fi

DB=$1
DB_OUTPUT_PATH="${repo_path}/foldseek_db/foldseek_${DB}"
DB_OUTPUT_TEMP_PATH="${repo_path}/foldseek_db/temp_foldseek_${DB}"
INPUT_FOLDER="${repo_path}/foldseek_db/viridiplantae_prot_db_${DB}"
echo "Output path: $DB_OUTPUT_PATH"
echo "Input folder: $INPUT_FOLDER"
mkdir -p $DB_OUTPUT_PATH
mkdir -p $DB_OUTPUT_TEMP_PATH
DB_OUTPUT="${DB_OUTPUT_PATH}/foldseek_${DB}"
DB_OUTPUT_TEMP="${DB_OUTPUT_TEMP_PATH}/foldseek_${DB}"

DB_PATHS=(${INPUT_FOLDER}/viridiplantae_prot_db_${DB}_batch_*)
DBS=()
if [ ${#DB_PATHS[@]} -eq 0 ]; then
    echo "Error: No matching database folders found for pattern: ${INPUT_FOLDER}/viridiplantae_prot_db_${DB}_batch_*"
    exit 1
fi
echo "Found ${#DB_PATHS[@]} database folders to process."

for path in "${DB_PATHS[@]}"; do
    DB_NAME="${path##*/}"
    if [[ "$path" == *"_1" ]]; then
        DB_1="${path}/${DB_NAME}"
    fi

    if [[ "$path" == *"_2" ]]; then
        DB_2="${path}/${DB_NAME}"
    fi

    if [[ -d "$path" && ! "$path" =~ _[12]$ ]]; then
        DBS+=("${path}/${DB_NAME}")
    fi
done

echo -e "Concatenate ${DB_1} with ${DB_2}"
foldseek concatdbs $DB_1 $DB_2 $DB_OUTPUT_TEMP
foldseek concatdbs ${DB_1}_ss ${DB_2}_ss ${DB_OUTPUT_TEMP}_ss 
foldseek concatdbs ${DB_1}_h ${DB_2}_h ${DB_OUTPUT_TEMP}_h

for DB in "${DBS[@]}"; do
    foldseek concatdbs $DB $DB_OUTPUT_TEMP $DB_OUTPUT
    echo -e "First merge done\n\n\n"
    foldseek concatdbs ${DB}_ss ${DB_OUTPUT_TEMP}_ss ${DB_OUTPUT}_ss
    echo -e "Second merge _ss done\n\n\n"
    foldseek concatdbs ${DB}_h ${DB_OUTPUT_TEMP}_h ${DB_OUTPUT}_h
    echo -e "Third merge _h done\n\n\n"
    rm -rf $DB_OUTPUT_TEMP_PATH
    mv $DB_OUTPUT_PATH $DB_OUTPUT_TEMP_PATH
    mkdir -p $DB_OUTPUT_PATH
done

mv $DB_OUTPUT_TEMP_PATH $DB_OUTPUT_PATH

if [[ "$USER" == "tronkko" ]]; then
    mv ~/temp/logs/${SLURM_JOB_ID}_out.txt ${DB_OUTPUT_PATH}/merge_${DB}_${SLURM_JOB_ID}_out.txt
    mv ~/temp/logs/${SLURM_JOB_ID}_err.txt ${DB_OUTPUT_PATH}/merge_${DB}_${SLURM_JOB_ID}_err.txt
fi

