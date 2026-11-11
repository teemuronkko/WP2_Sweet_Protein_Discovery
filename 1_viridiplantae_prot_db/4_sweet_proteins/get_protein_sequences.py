import os
import sys 

# Custom imports from utils.py
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_path = script_dir.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"
sys.path.append(repo_path)
from utils import read_seq_from_pdb


if __name__ == "__main__":
    sweet_protein_structure_dir = f"{repo_path}/1_viridiplantae_prot_db/4_sweet_proteins/cleaned_pdbs/"
    sweet_protein_structures = os.listdir(sweet_protein_structure_dir)
    sweet_protein_structures = [f for f in sweet_protein_structures if f.endswith(".pdb")]
    sweet_protein_fastas = f"{repo_path}/1_viridiplantae_prot_db/4_sweet_proteins/sweet_protein_fastas/"
    os.makedirs(sweet_protein_fastas, exist_ok=True)

    for pdb_file in sweet_protein_structures:
        pdb_path = os.path.join(sweet_protein_structure_dir, pdb_file)
        if not pdb_path.endswith(".pdb"):
            continue
        sequence = read_seq_from_pdb(pdb_path)
        fasta_path = os.path.join(sweet_protein_fastas, pdb_file.replace(".pdb", ".fasta"))
        with open(fasta_path, "w") as f:
            f.write(f">{pdb_file.replace('.pdb', '')}\n")
            f.write(sequence + "\n")