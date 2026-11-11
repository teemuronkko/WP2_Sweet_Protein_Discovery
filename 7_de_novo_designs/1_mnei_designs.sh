#!/bin/bash
source activate mlfold
folder_with_pdbs="./1_input_pdb/"
output_dir="./1_output"
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

# Input paths
path_for_parsed_chains=$output_dir"/parsed_pdb.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdb.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdb.jsonl"
chains_to_design="A"
fixed_positions="2 23 28 41 65 67 68"
num_seqs=20
seed=10
python ./1_ProteinMPNN_helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains
python ./1_ProteinMPNN_helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
python ./1_ProteinMPNN_helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

# Run ProteinMPNN over multiple sampling temperatures
for temp in 0.1 0.15 0.2 0.25 0.3; do
    folder_suffix="temp_${temp//./}"
    
    echo "Running ProteinMPNN at sampling temp: $temp"
    
    python ./protein_mpnn_run.py \
        --jsonl_path "$path_for_parsed_chains" \
        --chain_id_jsonl "$path_for_assigned_chains" \
        --fixed_positions_jsonl "$path_for_fixed_positions" \
        --out_folder "$output_dir/fixed_$folder_suffix" \
        --num_seq_per_target "$num_seqs" \
        --sampling_temp "$temp" \
        --seed "$seed" \
        --batch_size 1 \
        --use_soluble_model

    python ./protein_mpnn_run.py \
        --jsonl_path "$path_for_parsed_chains" \
        --chain_id_jsonl "$path_for_assigned_chains" \
        --out_folder "$output_dir/free_$folder_suffix" \
        --num_seq_per_target "$num_seqs" \
        --sampling_temp "$temp" \
        --seed "$seed" \
        --batch_size 1 \
        --use_soluble_model
done