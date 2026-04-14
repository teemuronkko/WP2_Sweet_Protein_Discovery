import requests
from requests.adapters import HTTPAdapter, Retry
import re
import os
import gzip
import shutil
import pandas as pd
import subprocess

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def gunzip(path, output_path=None, remove_original=False):
    """
    Unzip a .gz compressed file.
    
    Args:
        path (str): Path to the .gz file.
        output_path (str): Path to save the uncompressed file.
        remove_original (bool): Whether to remove the original .gz file after unzipping.
        
    Returns:
        path (str): Path to the uncompressed file.
    """
    if output_path is None:
        output_path = path[:-3]  # Remove .gz extension

    # Unzip the file
    with gzip.open(path, 'rb') as input_file:
        with open(output_path, 'wb') as output_file:
            shutil.copyfileobj(input_file, output_file)

    # Remove original file if specified
    if remove_original and os.path.exists(output_path):
        os.remove(path)

    return output_path

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)

def fetch_uniprot_fasta(url, output_file):
    # If output_file doesn't exist, create it
    if not os.path.exists(output_file):
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

    for batch, total in get_batch(url):
        # Write entries to the output_file
        with open(output_file, "a") as f:
            f.write(batch.text)

def get_variants(fasta_url: str, output_path: str, remove_original: bool = False):
    """
    Fetches fasta sequences from the given URL using wget.

    Args:
        fasta_url (str): URL to fetch fasta sequences from.
        output_path (str): Path to save the uncompressed fasta file.
        remove_original (bool): Whether to remove the original compressed file after unzipping.

    Returns:
        pd.DataFrame: DataFrame with columns 'ID', 'Sequence'.
    """

    # Build paths
    output_dir = os.path.dirname(output_path)
    fasta_path = os.path.join(output_dir, "variants.fasta.gz")

    # Download the fasta file using wget
    subprocess.run(["wget", "-O", fasta_path, fasta_url], check=True)

    # Unzip the downloaded file
    gunzip(fasta_path, output_path, remove_original=remove_original)

    return pd.read_csv(output_path, sep="\t")


if __name__ == "__main__":

    # Get script directory and build output file path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_file = f"{script_dir}/uniprot_viridiplantae.fasta"
    url = "https://rest.uniprot.org/uniprotkb/search?format=fasta&includeIsoform=true&query=%28%28taxonomy_id%3A33090%29%29&size=500"

    # Fetch UniProt plant proteins and save to file
    fetch_uniprot_fasta(url, output_file)

    # Fetch toxin proteins
    toxin_url = "https://rest.uniprot.org/uniprotkb/search?format=fasta&query=%28%28keyword%3AKW-0800%29%29&size=500"
    toxin_output_file = f"{script_dir}/uniprot_toxins.fasta"
    fetch_uniprot_fasta(toxin_url, toxin_output_file)

    # Fetch allergen proteins
    allergen_url = "https://rest.uniprot.org/uniprotkb/search?format=fasta&query=%28%28keyword%3AKW-0020%29%29&size=500"
    allergen_output_file = f"{script_dir}/uniprot_allergens.fasta"
    fetch_uniprot_fasta(allergen_url, allergen_output_file)

    # Fetch variants and print the number of sequences
    variants_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz"
    variants_path = f"{script_dir}/uniprot_sprot_varsplic.fasta"
    variant_path = get_variants(variants_url, output_path=variants_path, remove_original=True)
