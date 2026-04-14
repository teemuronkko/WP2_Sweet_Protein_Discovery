# Helper script to save FASTA sequences for OneKP proteins with hits found using foldseek-search.
import os
import pandas as pd

if __name__ == "__main__":

    # Get directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_name = "WP2_SweetProteinDiscovery"
    repo_dir = script_dir.split(repo_name)[0] + repo_name
    output_dir = f"{repo_dir}/5_openfold3_models/fastas/onekp_hits/"
    os.makedirs(output_dir, exist_ok = True)

    # Read OneKP results
    onekp_results_path = f"{repo_dir}/3_foldseek_search/foldseek_OneKP_search.tsv"
    onekp_df = pd.read_csv(onekp_results_path, sep = "\t", header = 0)
    onekp_df = onekp_df[["target", "tseq"]].drop_duplicates()

    for index, row in onekp_df.iterrows():
        id = row["target"]
        seq = row["tseq"]
        if "X" in seq:
            print("WARNING: ", id, "contains an X.")
        with open(os.path.join(output_dir, f"{id}.fasta"), "w") as f:
            f.write(f">{id}\n{seq}\n")