# Script to remove redundant FASTA entries from UniProtKB hits or sequences with non-canonical amino acids.
import os

def check_canonical_aa(sequence):
    """
    Check if a sequence contains only canonical amino acids.

    Args:
        sequence (str): Amino acid sequence.
    Returns:
        bool: True if sequence contains only canonical amino acids, False otherwise.
    """
    AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
    return all(residue in AMINO_ACIDS for residue in sequence)

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

    # Read in list of redundant UniProt IDs to skip
    redundant_uniprot_ids = set()
    redundant_ids_path = f"{repo_dir}/4_afdb_structures/proteins_to_remove_from_UniProtKB.txt"
    if os.path.exists(redundant_ids_path):
        with open(redundant_ids_path, "r") as f:
            for line in f:
                redundant_uniprot_ids.add(line.strip())

    uniprot_dir = f"{repo_dir}/5_openfold3_models/fastas/uniprot_missing_structures/"
    redundant_fasta_dir = f"{repo_dir}/5_openfold3_models/fastas/redundant_uniprot_sequences/"
    os.makedirs(redundant_fasta_dir, exist_ok=True)

    # Set maximum sequence length to 1000
    max_sequence_length = 1000

    # List fasta files in uniprot_dir
    uniprot_fasta_files = [os.path.join(root, f) for root, dirs, files in os.walk(uniprot_dir) for f in files if f.endswith(".fasta")]
    print("Total FASTA files before filtering:", len(uniprot_fasta_files))
    uniprot_ids_to_remove = set()

    for fasta_file in uniprot_fasta_files:
        fasta_dict = fasta_to_dict(fasta_file)
        for seq_id, sequence in fasta_dict.items():
            uniprot_id = seq_id.split("|")[1] if "|" in seq_id else seq_id
            if uniprot_id in redundant_uniprot_ids or len(sequence) > max_sequence_length:
                print(f"Removing redundant UniProt ID: {uniprot_id}, sequence length: {len(sequence)}")
                os.rename(fasta_file, os.path.join(redundant_fasta_dir, os.path.basename(fasta_file)))
                uniprot_ids_to_remove.add(uniprot_id)
                continue

            if not check_canonical_aa(sequence):
                print(f"Removing non-canonical UniProt ID: {uniprot_id}, sequence length: {len(sequence)}")
                os.rename(fasta_file, os.path.join(redundant_fasta_dir, os.path.basename(fasta_file)))
                uniprot_ids_to_remove.add(uniprot_id)
                continue

    # Same filtering for onekp_hits
    onekp_dir = f"{repo_dir}/5_openfold3_models/fastas/onekp_hits/"
    redundant_onekp_fasta_dir = f"{repo_dir}/5_openfold3_models/fastas/redundant_onekp_sequences/"
    os.makedirs(redundant_onekp_fasta_dir, exist_ok=True)
    onekp_fasta_files = [os.path.join(root, f) for root, dirs, files in os.walk(onekp_dir) for f in files if f.endswith(".fasta")]
    print("\n\n\n\nTotal OneKP FASTA files before filtering:", len(onekp_fasta_files))
    onekp_ids_to_remove = set()

    for fasta_file in onekp_fasta_files:
        fasta_dict = fasta_to_dict(fasta_file)
        for seq_id, sequence in fasta_dict.items():
            onekp_id = seq_id
            if "+" in onekp_id:
                # Replace + with _
                onekp_fix_id = onekp_id.replace("+", "_")
                # Rename fasta file accordingly
                new_fasta_file = os.path.join(os.path.dirname(fasta_file), onekp_fix_id + ".fasta")
                os.rename(fasta_file, new_fasta_file)
                print(f"Renaming {fasta_file} to {new_fasta_file}")
                fasta_file = new_fasta_file

            if not check_canonical_aa(sequence) or len(sequence) > max_sequence_length:
                print(f"Removing non-canonical OneKP ID: {onekp_id}, sequence length: {len(sequence)}")
                os.rename(fasta_file, os.path.join(redundant_onekp_fasta_dir, os.path.basename(fasta_file)))
                onekp_ids_to_remove.add(onekp_id)
                continue

    print("Total UniProt FASTA files removed:", len(uniprot_ids_to_remove))
    print("Total OneKP FASTA files removed:", len(onekp_ids_to_remove))