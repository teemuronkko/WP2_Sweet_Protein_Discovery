

uniprot_redundant = "/Users/pqh443/Documents/Git_projects/WP2_SweetProteinDiscovery/5_openfold3_models/fastas/redundant_uniprot_sequences"

# Read file names of .fasta files in the directory
import os
fasta_files = [f.split('.')[0] for f in os.listdir(uniprot_redundant) if f.endswith('.fasta')]
fasta_files = [f.lower() for f in fasta_files]
fasta_files = [f for f in fasta_files if f]

af3_path = "/Users/pqh443/Documents/Git_projects/WP2_SweetProteinDiscovery/5_alphafold3_models/AF3_out"

# list folders in af3_path
af3_folders = [f for f in os.listdir(af3_path) if os.path.isdir(os.path.join(af3_path, f))]

# Find redudnant folders    
for folder in af3_folders:
    folder_lower = folder.lower()
    if folder_lower in fasta_files:
        print(f"Removing redundant folder: {folder}")
        folder_path = os.path.join(af3_path, folder)
        os.system(f"rm -rf '{folder_path}'")

    # If folder starts with "onekp", move to OneKP folder
    if folder_lower.startswith("onekp"):
        outfolder = "/Users/pqh443/Documents/Git_projects/WP2_SweetProteinDiscovery/5_alphafold3_models/OneKP/models"
        outpath = os.path.join(outfolder, folder)
        print(f"Moving OneKP folder: {folder} to {outpath}")
        folder_path = os.path.join(af3_path, folder)
        os.rename(folder_path, outpath)


# Remove jsons for redundant sequences
af3_in_path = "/Users/pqh443/Documents/Git_projects/WP2_SweetProteinDiscovery/5_alphafold3_models/AF3_in"
af3_jsons = [f for f in os.listdir(af3_in_path) if f.endswith('.json')]

for json_file in af3_jsons:
    json_name = json_file.split('.')[0].lower()
    if json_name in fasta_files:
        print(f"Removing redundant json file: {json_file}")
        json_path = os.path.join(af3_in_path, json_file)
        os.remove(json_path)

    if json_name.startswith("onekp"):
        # move to json folder
        outfolder = "/Users/pqh443/Documents/Git_projects/WP2_SweetProteinDiscovery/5_alphafold3_models/OneKP/jsons"
        outpath = os.path.join(outfolder, json_file)
        print(f"Moving OneKP json file: {json_file} to {outpath}")
        json_path = os.path.join(af3_in_path, json_file)
        os.rename(json_path, outpath)