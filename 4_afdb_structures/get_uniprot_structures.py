import requests
import os
import pandas as pd
from tqdm import tqdm

def get_afdb_structure(uid: str, output_dir: str, version: int = 6) -> str:
    """
    Downloads a protein structure from the AlphaFold Database.
    Args:
        uid (str): UniProt ID of the protein.
        output_dir (str): Directory to save the downloaded PDB file.
        version (int): Version of the AlphaFold model to download. Default is 6.
    Returns:
        str: Path to the downloaded PDB file.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v{version}.pdb"
    if output_dir is None:
        output_dir = "."
    os.makedirs(output_dir, exist_ok=True)
    
    output_file = os.path.join(output_dir, f"AF-{uid}-F1-model_v{version}.pdb")   
    if os.path.exists(output_file):
        print(f"AlphaFold structure for {uid} already exists at {output_file}")
        return output_file

    response = requests.get(url)
    response.raise_for_status()
    with open(output_file, "wb") as f:
        f.write(response.content)

    if os.path.exists(output_file):
        print(f"Downloaded AlphaFold structure for {uid} to {output_file}")
        return output_file
    else:
        print("Error: PDB file was not created.")
        return 
    
def fasta_to_dict(fasta_path):
    """
    Reads a fasta file to a dictionary with sequence names as keys and sequences as items.

    Args:
        fasta_path (str): Path to input fasta file

    Returns:
        dict: Dictionary of sequences
    """
    fasta_dict = {}
    with open(fasta_path, 'r') as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                if seq:
                    fasta_dict[key] = seq
                    seq = ''
                key = line.replace(">", "").strip()
            else:
                seq += line.strip()
        fasta_dict[key] = seq
    return fasta_dict

if __name__ == "__main__":

    # Get directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_name = "WP2_SweetProteinDiscovery"
    repo_dir = script_dir.split(repo_name)[0] + repo_name

    # Make output folders for missing structures
    missing_structures_dir = f"{repo_dir}/5_openfold3_models/fastas/uniprot_missing_structures"
    os.makedirs(missing_structures_dir, exist_ok=True)
    os.makedirs(os.path.join(missing_structures_dir, "toxins"), exist_ok=True)
    os.makedirs(os.path.join(missing_structures_dir, "allergens"), exist_ok=True)
    os.makedirs(os.path.join(missing_structures_dir, "uniprot_hits"), exist_ok=True)

    # Read in fasta files to dictionaries
    toxin_dict = fasta_to_dict(f"{repo_dir}/1_viridiplantae_prot_db/2_uniprot_data/uniprot_toxins.fasta")
    allergen_dict = fasta_to_dict(f"{repo_dir}/1_viridiplantae_prot_db/2_uniprot_data/uniprot_allergens.fasta")

    # Read UniProt hits from the Foldseek search results
    uniprot_tsv_path = f"{repo_dir}/3_foldseek_search/foldseek_UniProt_search.tsv"
    uniprot_hits_df = pd.read_csv(uniprot_tsv_path, sep="\t")
    uniprot_hits_df = uniprot_hits_df[["target", "tseq"]].drop_duplicates()
    output_dir = f"{repo_dir}/4_afdb_structures/uniprot_hits_afdb_structures/pdb_dir"
    os.makedirs(output_dir, exist_ok=True)

    process_uniprot_hits = True
    process_toxins = True
    process_allergens = True

    # Loop through each UniProt ID and download the corresponding AFDB structure
    if process_uniprot_hits:
        for _, row in tqdm(uniprot_hits_df.iterrows(), total=len(uniprot_hits_df), desc="Downloading AFDB structures"):
            uniprot_id = row["target"].split("|")[1]
            try:
                if not os.path.exists(os.path.join(output_dir, f"AF-{uniprot_id}-F1-model_v6.pdb")):
                    get_afdb_structure(uniprot_id, output_dir)
                else:
                    print(f"AlphaFold structure for {uniprot_id} already exists. Skipping download.")
            except requests.HTTPError as e:
                print(f"Failed to download structure for {uniprot_id}: {e}")
                if "X" not in row["tseq"]:
                    with open(os.path.join(missing_structures_dir, "uniprot_hits", f"{uniprot_id}.fasta"), "w") as f:
                        f.write(f">{uniprot_id}\n{row['tseq']}\n")

    # UniProt toxins 
    if process_toxins:
        toxin_fa = f"{repo_dir}/1_viridiplantae_prot_db/2_uniprot_data/uniprot_toxins.fasta"
        output_dir = f"{repo_dir}/4_afdb_structures/toxindb_afdb_structures/pdb_dir"
        os.makedirs(output_dir, exist_ok=True)
        with open(toxin_fa, "r") as input_f:
            for line in tqdm(input_f, desc="Downloading toxin AFDB structures"):
                if line.startswith(">"):
                    uniprot_id = line.split("|")[1]
                    try:
                        get_afdb_structure(uniprot_id, output_dir)
                    except requests.HTTPError as e:
                        print(f"Failed to download structure for {uniprot_id}: {e}")
                        seq_name = line.replace(">", "").strip()
                        seq = toxin_dict[seq_name]
                        if "X" not in seq:
                            with open(os.path.join(missing_structures_dir, "toxins", f"{uniprot_id}.fasta"), "w") as f:
                                f.write(f">{uniprot_id}\n{toxin_dict[seq_name]}\n")
    
    # UniProt allergens 
    if process_allergens:
        allergen_fa = f"{repo_dir}/1_viridiplantae_prot_db/2_uniprot_data/uniprot_allergens.fasta"
        output_dir = f"{repo_dir}/4_afdb_structures/allergendb_afdb_structures/pdb_dir"
        os.makedirs(output_dir, exist_ok=True)
        with open(allergen_fa, "r") as input_f:
            for line in tqdm(input_f, desc="Downloading allergen AFDB structures"):
                if line.startswith(">"):
                    uniprot_id = line.split("|")[1]
                    try:
                        get_afdb_structure(uniprot_id, output_dir)
                    except requests.HTTPError as e:
                        print(f"Failed to download structure for {uniprot_id}: {e}")
                        seq_name = line.replace(">", "").strip()
                        seq = allergen_dict[seq_name]
                        if "X" not in seq:
                            with open(os.path.join(missing_structures_dir, "allergens", f"{uniprot_id}.fasta"), "w") as f:
                                f.write(f">{uniprot_id}\n{seq}\n")
