
from urllib import request
import os
import sys

# Custom imports from utils.py
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_path = script_dir.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"
sys.path.append(repo_path)
from utils import crop_pdb, cif_to_pdb, get_afdb_structure

if __name__ == "__main__":
    # Make a directory for pdb files
    pdb_dir = f"{repo_path}/1_viridiplantae_prot_db/4_sweet_proteins/sweet_protein_pdbs"
    cif_dir = f"{repo_path}/1_viridiplantae_prot_db/4_sweet_proteins/cifs"
    os.makedirs(pdb_dir, exist_ok=True)

    # Dictionary mapping sweet protein names to their PDB IDs or file paths
    sweet_proteins = {
        "MNEI": "1IV7", 
        "Brazzein": "4HE7",
        "Thaumatin-I": "1LR2", 
        "Thaumatin-II": "3AOK",
        "Miraculin": f"{pdb_dir}/Miraculin.pdb",
        "Honey_truffle_sweetener": f"{pdb_dir}/Honey_truffle_sweetener.pdb", 
    }

    ### Get protein structures for Miraculin and Honey truffle sweetener  
    # For Honey truffle sweetener, convert from CIF to PDB and save
    cif_to_pdb(f"{cif_dir}/jf4c04368_si_003.cif", f"{pdb_dir}/Honey_truffle_sweetener.pdb")
    # Crop Miraculin structure from AlphaFold prediction from AlphaFold DB
    miraculin_pdb = get_afdb_structure("P13087", f"{pdb_dir}/AF-P13087-F1-model_v6.pdb")
    crop_pdb(miraculin_pdb, f"{pdb_dir}/Miraculin.pdb", "A", 30, 220)
    os.remove(miraculin_pdb)

    ### Download PDB files for other sweet proteins from RCSB PDB
    for protein, pdb_id in sweet_proteins.items():
        # Skip already processed local files
        if pdb_id.startswith(f"{pdb_dir}/"):
            continue  

        # Download PDB file from RCSB PDB
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        outfile = f"{pdb_dir}/{protein}_{pdb_id}.pdb"
        request.urlretrieve(url, outfile)
        print("Downloaded:", outfile)

        # For MNEI, use chain A only
        if protein == "MNEI":
            crop_pdb(f"{pdb_dir}/{protein}_{pdb_id}.pdb", outfile, "A", None, None)
