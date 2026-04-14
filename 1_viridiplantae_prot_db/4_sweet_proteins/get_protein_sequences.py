# For file in cleaned_pdbs, read sequence from PDB and save to FASTA in sweet_protein_fastas
import os

def read_seq_from_pdb(pdb_path, pdb_id = "structure", chain_id = "A"):
    """
    Reads the sequence from a PDB file.

    Args:
        pdb_path (str): Path to the PDB file.

    Returns:
        str: Amino acid sequence extracted from the PDB file.
    """

    from Bio.PDB import PDBParser
    from Bio.SeqUtils import seq1

    pdbparser = PDBParser()
    structure = pdbparser.get_structure(pdb_id, pdb_path)
    chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}

    return chains[chain_id]


# Get directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_name = "WP2_SweetProteinDiscovery"
repo_dir = script_dir.split(repo_name)[0] + repo_name

sweet_protein_structure_dir = f"{repo_dir}/1_viridiplantae_prot_db/4_sweet_proteins/cleaned_pdbs/"
sweet_protein_structures = os.listdir(sweet_protein_structure_dir)
sweet_protein_structures = [f for f in sweet_protein_structures if f.endswith(".pdb")]
sweet_protein_fastas = f"{repo_dir}/1_viridiplantae_prot_db/4_sweet_proteins/sweet_protein_fastas/"
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