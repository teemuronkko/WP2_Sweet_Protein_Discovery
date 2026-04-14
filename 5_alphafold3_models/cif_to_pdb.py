# Helper script to convert CIF files to PDB format
import os
from Bio.PDB import MMCIFParser, PDBIO

def cif_to_pdb(cif_file, pdb_file):
    """
    Converts a CIF file to PDB format.

    Args:
        cif_file (str): Path to the input CIF file.
        pdb_file (str): Path to the output PDB file.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_file)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

if __name__ == "__main__":

    # Get directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_name = "WP2_SweetProteinDiscovery"
    repo_dir = script_dir.split(repo_name)[0] + repo_name

    # Apply above function on all CIF files in a directory
    dbs = ["OneKP", "UniProt"]
    success = 0
    cif_total = 0

    # List all CIF files in the directories and convert them
    for db in dbs:
        
        model_path = f"{repo_dir}/5_alphafold3_models/{db}/models"
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