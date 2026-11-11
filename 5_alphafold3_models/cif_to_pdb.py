import os
import sys

# Custom imports from utils.py
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_path = script_dir.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"
sys.path.append(repo_path)
from utils import cif_to_pdb

if __name__ == "__main__":
    # Apply above function on all CIF files in a directory
    dbs = ["OneKP", "UniProt"]
    success = 0
    cif_total = 0

    # List all CIF files in the directories and convert them
    for db in dbs: 
        model_path = f"{repo_path}/5_alphafold3_models/{db}/models"
        folders = [f for f in os.listdir(model_path) if os.path.isdir(os.path.join(model_path, f))]
        cif_files = []
        for folder in folders:
            folder_path = os.path.join(model_path, folder)
            cif_files += [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.cif')]
    
        # Print number of CIF files found
        print(f"Found {len(cif_files)} CIF files in {model_path}")
        cif_total += len(cif_files)
        for cif_path in cif_files:
            pdb_file = os.path.basename(cif_path).replace('.cif', '.pdb')
            pdb_path = os.path.join(os.path.dirname(cif_path), pdb_file)
            cif_to_pdb(cif_path, pdb_path)
            if os.path.exists(pdb_path):
                success += 1

    print(f"Successfully converted {success} CIF files to PDB format.")
    print(f"Total CIF files processed: {cif_total}")