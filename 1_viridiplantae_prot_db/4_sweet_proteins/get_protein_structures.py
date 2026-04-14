
from Bio.PDB import MMCIFParser, PDBParser, PDBIO
from urllib import request
import os

def crop_pdb(pdb: str, output_file: str, chain_id: str, start_idx: int, end_idx: int) -> str:
    """
    Crops a PDB file to include a specific range and chain, given the start and end residue indices.
    Args:
        pdb (str): Path to the input PDB file.
        output_file (str): Path to save the cropped PDB file.
        chain_id (str): Chain identifier to retain.
        start_idx (int): Starting residue index (inclusive).
        end_idx (int): Ending residue index (inclusive).

    Returns:
        str: Path to the cropped PDB file.
    """
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("structure", pdb)
    for model in structure:
        for chain in list(model):
            if chain.id != chain_id:
                model.detach_child(chain.id)
                continue

            for residue in list(chain):
                if start_idx is None and end_idx is None:
                    continue
                elif not residue.id[1] in range(start_idx, end_idx + 1):
                    chain.detach_child(residue.id)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

    if os.path.exists(output_file):
        print(f"Cropped PDB saved to {output_file}")
        return output_file
    else:
        print("Error: Cropped PDB file was not created.")
        return None
    

def cif_to_pdb(cif_file: str, pdb_file: str) -> str:
    """
    Converts a CIF file to PDB format.
    Args:
        cif_file (str): Path to the input CIF file.
        pdb_file (str): Path to save the converted PDB file.

    Returns:
        str: Path to the converted PDB file.
    """
    parser = MMCIFParser(QUIET=True)   # QUIET avoids warnings
    structure = parser.get_structure("structure", cif_file)

    # Write PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

    if os.path.exists(pdb_file):
        print(f"Converted PDB saved to {pdb_file}")
        return pdb_file
    else:
        print("Error: Converted PDB file was not created.")
        return None
    
def get_afdb_structure(structure_id: str, output_file: str) -> str:
    """
    Downloads a protein structure from the AlphaFold Database.
    Args:
        structure_id (str): AlphaFold structure ID.
        output_file (str): Path to save the downloaded PDB file.
    Returns:
        str: Path to the downloaded PDB file.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{structure_id}-F1-model_v6.pdb"
    request.urlretrieve(url, output_file)
    return output_file


if __name__ == "__main__":

    # Get directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_name = "WP2_SweetProteinDiscovery"
    repo_dir = script_dir.split(repo_name)[0] + repo_name

    # Make a directory for pdb files
    pdb_dir = f"{repo_dir}/1_viridiplantae_prot_db/4_sweet_proteins/sweet_protein_pdbs"
    cif_dir = f"{repo_dir}/1_viridiplantae_prot_db/4_sweet_proteins/cifs"
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
