
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

input_file = "/Users/pqh443/Documents/Git_projects/WP2_SweetProteinDiscovery/5_alphafold3_models/missing_onekp_models.txt"

# Read lines 
with open(input_file, "r") as infile:
    lines = infile.readlines()
lines = [line.strip() for line in lines if line.strip()]
lines = [line.split("skipping ")[1] for line in lines]

fasta_p = "/Users/pqh443/Documents/Git_projects/WP2_SweetProteinDiscovery/1_viridiplantae_prot_db/viridiplantae_prot_db_OneKP.fasta"

# Read fasta 
seqs = fasta_to_dict(fasta_p)

for seq in lines:
    sequence = seqs[seq]
    print(len(sequence))