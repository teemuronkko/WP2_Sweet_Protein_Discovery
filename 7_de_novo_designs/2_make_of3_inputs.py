import glob 
import os
import sys

# Custom imports from utils.py
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_path = script_dir.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"
sys.path.append(repo_path)
from utils import fasta_to_dict, make_of3_json

if __name__ == "__main__":
    # Make fastas for all sequences in the de novo design 
    temps = ["01", "015", "02", "025", "03"]
    modes = ["fixed", "free"]
    for mode in modes:
        for temp in temps:
            fasta_path = f"{repo_path}/7_de_novo_designs/output/{mode}_temp_{temp}/seqs/MNEI_1IV7.fa"
            fasta_dict = fasta_to_dict(fasta_path)
            count = 0
            for seq_name, seq in fasta_dict.items():
                if "MNEI_1IV7" in seq_name:
                    continue
                else:
                    count += 1
                    if count > 20:
                        continue
                    seq_name = f"{mode}_{temp}_MNEI_1IV7_{count}"
                fa_output = f"{repo_path}/7_de_novo_designs/fastas/{seq_name}.fasta"
                os.makedirs(os.path.dirname(fa_output), exist_ok=True)
                if os.path.exists(fa_output):
                    continue
                else:
                    with open(fa_output, 'w') as f:
                        f.write(f">{seq_name}\n{seq}\n")

    # Make OF3 input jsons for all sequences in the de novo design
    fasta_dir = f"{repo_path}/7_de_novo_designs/fastas"
    output_dir = f"{repo_path}/7_de_novo_designs/of3_jsons"
    os.makedirs(output_dir, exist_ok=True)
    for fasta in glob.glob(f"{fasta_dir}/*.fasta"):
        fasta_dict = fasta_to_dict(fasta)
        for seq_name, sweet_protein_seq in fasta_dict.items():
            chain_dict = {
                "A": {"seq": sweet_protein_seq}    
            }
            if not os.path.exists(f"{output_dir}/{seq_name}.json"):
                make_of3_json(chain_dict, model_name=seq_name, num_recycles=3, json_path=f"{output_dir}/{seq_name}.json", use_paired_msas=True, use_msas=False)