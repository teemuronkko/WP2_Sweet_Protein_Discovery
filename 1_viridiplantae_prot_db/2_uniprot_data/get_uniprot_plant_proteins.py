
import requests
from requests.adapters import HTTPAdapter, Retry
import re
import os

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

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

def fetch_uniprot_plant_proteins(url, output_file, columns):
    # Open file if it doesnt exist and write header
    if not os.path.exists(output_file):
        with open(output_file, "w") as f:
            f.write("\t".join(columns) + "\n")

    for batch, total in get_batch(url):
        results = []
        for line in batch.text.splitlines()[1:]:
            results.append(line.split("\t"))

        # Write results in the file after each batch
        with open(output_file, "a") as f:
            for entry in results:
                f.write("\t".join(entry) + "\n")

        print(len(results), "/", total, "rows processed.")

if __name__ == "__main__":
    # Get script directory and build output file path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_file = f"{script_dir}/uniprot_viridiplantae.tsv"

    # Define UniProt query URL and columns to fetch
    url = "https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Csequence%2Cft_var_seq&format=tsv&query=%28%28taxonomy_id%3A33090%29%29&size=500"
    cols = ["Accession","Reviewed", "ID","Protein names","Gene names","Organism","Length","Sequence","Alternative Sequence"]

    # Fetch UniProt plant proteins and save to file
    fetch_uniprot_plant_proteins(url, output_file, cols)
