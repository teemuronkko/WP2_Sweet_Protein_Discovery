#!/bin/bash -l
#SBATCH --partition=standard
#SBATCH --qos=normal
#SBATCH --exclude=ilfgridgpun02fl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/projects/ilfgrid/people/pqh443/AF3/temp/logs/%j_out.txt
#SBATCH --error=/projects/ilfgrid/people/pqh443/AF3/temp/logs/%j_err.txt
#SBATCH --job-name=of3_model
#SBATCH --gres=gpu:1

Help()
{
  echo ""
  echo "Usage: sbatch [OPTIONS(0)...] [ : [OPTIONS(N)...]] $0 [input_json_txt]"
  echo ""
  echo "  SLURM script to submit an OpenFold 3 modelling job."
  echo "  NOTE: -s and -e parameters must be given before positional arguments!"
  echo ""
  echo "  OPTIONS         : Arguments to pass for sbatch command (e.g., --job-name=test_json, --qos=prio, --partition=priority)"
  echo "  -s              : Start index for processing lines in the input_json_txt file (default: 1)"
  echo "  -e              : End index for processing lines in the input_json_txt file (default: number of lines in the file)"
  echo "  -h              : Display this help message"
  echo "  input_fasta_txt : Path to a .txt file containing paths to one or more json files (one path per line) "
}

module purge
module load miniconda/24.5.0
module load cuda/12.8
module load gcc
conda activate /projects/ilfgrid/apps/openfold3_preview2/envs/of3_env
source activate /projects/ilfgrid/apps/openfold3_preview2/envs/of3_env

export CUDA_HOME=/opt/software/cuda/12.8
export LIBRARY_PATH=/opt/software/cuda/12.8/lib64:/opt/software/cuda/12.8/nvvm/lib64:/opt/software/cuda/12.8/extras/CUPTI/lib64:/opt/software/cuda/12.8/extras/Debugger/lib64
export LD_LIBRARY_PATH=/opt/software/cuda/12.8/lib64:/opt/software/cuda/12.8/extras/Debugger/lib64:/opt/software/cuda/12.8/extras/CUPTI/lib64
export LDFLAGS="-L/projects/ilfgrid/apps/openfold3_preview2/envs/of3_env/lib $LDFLAGS"
export CFLAGS="-I/projects/ilfgrid/apps/openfold3_preview2/envs/of3_env/include $CFLAGS"

export PYTHONNOUSERSITE=1
export WORKSPACE="/projects/ilfgrid/apps/openfold3_preview2"
export PIP_CACHE_DIR="$WORKSPACE/cache/pip"
export CONDA_PKGS_DIRS="$WORKSPACE/cache/conda"
export TORCH_HOME="$WORKSPACE/cache/torch"
export HF_HOME="$WORKSPACE/cache/huggingface"
export CONDARC="$WORKSPACE/.condarc"
export TRITON_CACHE_DIR="$WORKSPACE/cache/triton"
export CUTLASS_PATH="/projects/ilfgrid/apps/openfold3_preview2/envs/of3_env/lib/python3.10/site-packages/cutlass_library/source"
export OPENFOLD_CACHE="/projects/ilfgrid/apps/openfold3_preview2/openfold3_cache"
export TORCH_EXTENSIONS_DIR="$CONDA_PREFIX/torch_extensions"

# Set default start and end line numbers
start_index=1
end_index=1

# Parse command-line options
while getopts ":s:e:h" opt; do
    case $opt in
      s)
         start_index="$OPTARG"
         ;;
      e)
         end_index="$OPTARG"
         ;;
      h)
         Help
         exit;;
  esac
done

# Validate start and end values
if [ "$start_index" -gt "$end_index" ]; then
    echo -e "Error: Start value cannot be greater than end value."
    exit 1
fi

# Shift to process positional arguments
shift $((OPTIND - 1))

# Get variables from this script
JSON_PATHS=$1

# Path to temporary log file folder - has to be the same folder as "--output" and "--error" path
LOG_PATH="/projects/ilfgrid/people/pqh443/OpenFold3/temp/logs"
mkdir -p $LOG_PATH

# Check if the json path file exists
if [ ! -e "$JSON_PATHS" ]; then
    echo "json path file does not exist: $JSON_PATHS"
    exit 1
fi

# Output directory for OpenFold predictions
OF3_OUTPUT_DIR="/projects/ilfgrid/people/pqh443/OpenFold3/Models"
OF3_DIR="/projects/ilfgrid/apps/openfold3_preview2"
CKPT_PATH="/projects/ilfgrid/apps/openfold3_preview2/openfold3_cache/of3-p2-155k.pt"

# Counter for current line number
current_line=0

# Open the file and read line by line
while IFS= read -r line; do

    # Increase line counter
    ((current_line++))

    # Check if the current line is within the specified range and move to the next line if not
    if [ "$current_line" -lt "$start_index" ] || [ "$current_line" -gt "$end_index" ]; then
        continue
    fi
    
    # Read in a json file from the list
    JSON_FILE=$line

    # Get the basename of the json file
    BASENAME=`basename $JSON_FILE .json`

    # Path to the log file
    LOG_FILE="${LOG_PATH}/${BASENAME}_log.txt"

    # Redirect standard output and standard error to the log file
    exec > "${LOG_FILE}" 2>&1

    # Start time 
    start=`date +%s`

    # Print output for SLURM log file
    echo -e ">> Running the following command:\n${0} ${@}\n"
    echo -e ">> Start time:  $(date)"
    echo -e "   Machine:     $(hostname)"
    echo -e "   Directory:   $(pwd)"
    echo -e "   GPU:         ${SLURM_JOB_GPUS}"

    # Echo output
    echo -e "\n>> OF3 input:"
    echo -e "   .json:      $(readlink -f $JSON_FILE)"
    
    # Make output directory
    mkdir -p "${OF3_OUTPUT_DIR}/${BASENAME}"
    
    # Check if prediction already exists and skip file if so
    if [ -f "${OF3_OUTPUT_DIR}/${BASENAME}/${BASENAME}/seed_42/${BASENAME}_seed_42_sample_5_model.cif" ]; then    
	    # Close the log file and stop redirecting
        exec >&-

        # Remove log file
        rm $LOG_FILE 
        
        # Continue to the next line in the file
        continue
    fi

    # Check if .json file exists and start prediction if so
    if [ -e $JSON_FILE ]; then
        echo -e ">> Starting prediction for $BASENAME \n"

        # Launch OF3 modelling 
        run_openfold predict \
            --query_json "$JSON_FILE" \
            --output_dir "${OF3_OUTPUT_DIR}/${BASENAME}" \
            --inference_ckpt_path "$CKPT_PATH" \
            --use_msa_server False \
            --use_templates False 
    else
        echo -e ">> Prediction is skipped because $JSON_FILE does not exist.\n"
    fi

    if [ -f "${OF3_OUTPUT_DIR}/${BASENAME}/${BASENAME}/seed_42/${BASENAME}_seed_42_sample_5_model.cif" ]; then
        pred_end=`date +%s`
        pred_runtime=$((pred_end-pred_start))
        echo -e ">> Prediction for $BASENAME finished in $pred_runtime seconds.\n"
    else     
        echo -e ">> Prediction failed, or ${OF3_OUTPUT_DIR}/${BASENAME}/${BASENAME}/seed_42/${BASENAME}_seed_42_sample_5_model.cif doesn't exist.\n"
    fi

    end=`date +%s`
    runtime=$((end-start))
    echo -e "\n>> Job finished at: $(date)"
    echo -e "\n>> Runtime: $runtime s"

    # Move log files
    mv "${LOG_FILE}" "${OF3_OUTPUT_DIR}/${BASENAME}/OF3_pred_${SLURM_JOB_ID}_log.txt"

    # Make sure the output files have the correct rights
    chmod 776 --recursive "${OF3_OUTPUT_DIR}/${BASENAME}/"

    # Close the log file and stop redirecting
    exec >&-

done < "$JSON_PATHS"

# Move log file for batch job
mv "${LOG_PATH}/${SLURM_JOB_ID}_err.txt" "${AF_OUTPUT_DIR}/OF3_pred_${SLURM_JOB_ID}_err.txt"
mv "${LOG_PATH}/${SLURM_JOB_ID}_out.txt" "${AF_OUTPUT_DIR}/OF3_pred_${SLURM_JOB_ID}_out.txt"
