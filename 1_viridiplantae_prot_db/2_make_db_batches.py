# Script to split sequence databases to smaller batches for ProstT5 processing
import os
from Bio import SeqIO   

def split_fasta_to_batches(input_fasta, output_dir, batch_size = 500000):
    """
    Splits a fasta file into smaller batches of specified size.
    Args:
        input_fasta (str): Path to the input fasta file.
        output_dir (str): Directory to save the output batches.
        batch_size (int): Number of sequences per batch.

    Returns:
        None
    """
    os.makedirs(output_dir, exist_ok=True)
    batch = []
    batch_number = 1
    fname = os.path.basename(input_fasta).split(".")[0]

    for i, record in enumerate(SeqIO.parse(input_fasta, "fasta")):
        batch.append(record)
        if (i + 1) % batch_size == 0:
            batch_file = os.path.join(output_dir, f"{fname}_batch_{batch_number}.fasta")
            SeqIO.write(batch, batch_file, "fasta")
            print(f"Written {batch_file} with {len(batch)} sequences.")
            batch = []
            batch_number += 1

    # Write remaining sequences
    if batch:
        batch_file = os.path.join(output_dir, f"{fname}_batch_{batch_number}.fasta")
        SeqIO.write(batch, batch_file, "fasta")
        print(f"Written {batch_file} with {len(batch)} sequences.")

if __name__ == "__main__":

    # Get directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_name = "WP2_SweetProteinDiscovery"
    repo_dir = script_dir.split(repo_name)[0] + repo_name
    batch_size = 500000

    # Paths to OneKP and UniProt Viridiplantae fasta files
    uniprot_viridiplantae_path = f"{repo_dir}/1_viridiplantae_prot_db/viridiplantae_prot_db_UniProt.fasta"
    oneKP_all_transcripts_path = f"{repo_dir}/1_viridiplantae_prot_db/viridiplantae_prot_db_OneKP.fasta"

    # Output directories for batches
    uniprot_output_dir = f"{repo_dir}/1_viridiplantae_prot_db/3_database_batches/UniProt_batches/"
    oneKP_output_dir = f"{repo_dir}/1_viridiplantae_prot_db/3_database_batches/OneKP_batches/"

    # Split the fasta files into batches
    print("Splitting UniProt Viridiplantae sequences into batches...")
    split_fasta_to_batches(uniprot_viridiplantae_path, uniprot_output_dir, batch_size=batch_size)
    print("Splitting OneKP Viridiplantae sequences into batches...")
    split_fasta_to_batches(oneKP_all_transcripts_path, oneKP_output_dir, batch_size=batch_size)