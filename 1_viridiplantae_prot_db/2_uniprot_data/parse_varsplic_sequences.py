import pandas as pd
import os
import sys

# Custom imports from utils.py
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_path = script_dir.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"
sys.path.append(repo_path)
from utils import read_sequences

def read_uniprot_ids(file_path):
    """
    Reads UniProt IDs from a TSV file. 

    Args:
        file_path (str): Path to the TSV file.
    Returns:
        set: A set of UniProt IDs.
    """
    df = pd.read_csv(file_path, sep="\t", usecols=["Accession"], dtype=str)
    return set(df["Accession"])

def read_uniprot_ids_from_fasta(file_path):
    """
    Reads UniProt IDs from a FASTA file. 

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        set: A set of UniProt IDs.
    """
    ids = set()
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                id = line.split('|')[1].split('-')[0]
                ids.add(id)
    return ids


if __name__ == "__main__":

    # Get script directory and build output file path
    output_file = f"{script_dir}/uniprot_viridiplantae.tsv"

    # Read all UniProt fasta IDs and sequences from varsplic fasta file
    varsplic_ids = read_uniprot_ids_from_fasta(f"{script_dir}/uniprot_sprot_varsplic.fasta")
    varsplic_sequences = read_sequences(f"{script_dir}/uniprot_sprot_varsplic.fasta")
    uniprot_ids = read_uniprot_ids(f"{script_dir}/uniprot_viridiplantae.tsv")
    toxin_ids = read_uniprot_ids_from_fasta(f"{script_dir}/uniprot_toxins.fasta")
    print(f"Total varsplic IDs: {len(varsplic_ids)}")
    print(f"Total uniprot plant IDs: {len(uniprot_ids)}")
    print(f"Total toxin IDs: {len(toxin_ids)}")
    print(f"Total non-toxin uniprot plant IDs: {len(uniprot_ids - toxin_ids)}")

    # Remove toxin sequences from uniprot_ids
    uniprot_ids = uniprot_ids - toxin_ids

    # Read uniprot plant sequneces and keep only non-toxin ones
    uniprot_plant_sequences = read_sequences(f"{script_dir}/uniprot_viridiplantae.fasta")
    uniprot_plant_sequences = uniprot_plant_sequences[uniprot_plant_sequences['ID'].apply(lambda x: x.split('|')[1].split('-')[0] in uniprot_ids)]

    # Get intersection of both sets
    common_ids = uniprot_ids.intersection(varsplic_ids)
    print(f"Total common IDs between uniprot plant and varsplic: {len(common_ids)}")

    ids_to_remove = []
    for index, row in varsplic_sequences.iterrows():
        if row['ID'].split('|')[1].split('-')[0] not in common_ids:
            ids_to_remove.append(index)

    # Drop rows with IDs not in common_ids
    varsplic_sequences = varsplic_sequences.drop(ids_to_remove).reset_index(drop=True)
    print(f"Total varsplic sequences after filtering: {len(varsplic_sequences)}")
    print(f"Final number of sequences to be written: {len(uniprot_plant_sequences) + len(varsplic_sequences)}")

    # Save sequences to a new fasta file
    with open(f"{script_dir}/uniprot_viridiplantae_varsplic.fasta", 'w') as f:
        for index, row in varsplic_sequences.iterrows():
            f.write(f">{row['ID']}\n{row['Sequence']}\n")

    # Make a merged fasta file of both original plant fasta and varsplic fasta
    with open(f"{script_dir}/uniprot_viridiplantae_all.fasta", 'w') as outfile:
        # Write original plant fasta sequences
        for index, row in uniprot_plant_sequences.iterrows():
            outfile.write(f">{row['ID']}\n{row['Sequence']}\n")

        # Write filtered varsplic sequences
        with open(f"{script_dir}/uniprot_viridiplantae_varsplic.fasta", 'r') as infile:
            for line in infile:
                outfile.write(line)