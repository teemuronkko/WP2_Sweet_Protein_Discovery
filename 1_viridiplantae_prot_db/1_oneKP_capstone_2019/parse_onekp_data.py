import pandas as pd 
import os
from Bio import SeqIO
import tarfile
import shutil

def read_id(tsv_path):
    """
    Reads a TSV file containing ID, Gene, and Description columns for each protein. 
    
    Args:
        tsv_path (str): Path to the TSV file.

    Returns:
        pd.DataFrame: DataFrame with columns 'ID', 'Gene', and 'Description'.
    """
    data = []
    with open(tsv_path, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) >= 3:
                data.append({"ID": cols[0], "Gene": cols[1], "Description": cols[2]})
    return pd.DataFrame(data)

def read_sequences(fasta_path):
    """
    Reads a fasta file containing protein sequences.
    
    Args:
        fasta_path (str): Path to the fasta file.

    Returns:
        pd.DataFrame: DataFrame with columns 'ID', 'Sequence'.
    """
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences.append({"ID": record.id, "Sequence": str(record.seq)})

    return pd.DataFrame(sequences)

def unzip_tarbz2_file(path, output_path):
    """
    Unzips a .tar.bz2 compressed file.

    Args:
        path (str): Path to the .tar.bz2 file.
        output_path (str): Path to save the uncompressed file.

    Returns:
        path (str): Path to the uncompressed file/folder.
    """

    if output_path is None:
        output_path = path.replace(".tar.bz2", "")
    
    # Remove existing folder if it exists
    if os.path.exists(output_path):
        shutil.rmtree(output_path)

    with tarfile.open(path, "r:bz2") as tar:
        extract_folder = output_path
        os.makedirs(extract_folder, exist_ok=True)
        tar.extractall(path=extract_folder)

    return output_path

def read_index(index_path):
    """
    Reads in index file mapping IDs to unique identifiers.
    
    Args:
        index_path (str): Path to the index file.

    Returns:
        pd.DataFrame: DataFrame mapping original IDs to unique identifiers.

    """
    index_df = pd.read_csv(index_path, sep="\t", header=None, names=["ID", "OneKP_ID"])
    return index_df

def gunzip(path, output_path=None):
    """
    Unzips a .gz compressed file.

    Args:
        path (str): Path to the .gz file.
        output_path (str): Path to save the uncompressed file.

    Returns:
        path (str): Path to the uncompressed file.
    """
    import gzip
    import shutil

    if output_path is None:
        output_path = path.replace(".gz", "")

    with gzip.open(path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return output_path

def process_folder(folder):
    output_path = None
    for file in os.listdir(folder):
        if file.endswith("translated.tar.bz2"):
            bz_path = os.path.join(folder, file)
            output_path = os.path.join(folder, file.replace(".tar.bz2", ""))

            # Remove existing folder if it exists
            if os.path.exists(output_path):
                shutil.rmtree(output_path)
            
            # Unzip the tar.bz2 file
            output_path = unzip_tarbz2_file(bz_path, output_path)
            break

    index_df = None    
    for file in os.listdir(folder):
        if file.endswith(".index"):
            index_df = read_index(os.path.join(folder, file))
            break

    seq_df, id_df = None, None
    for file in os.listdir(output_path):
        if file.endswith(".prots.out"):
            tsv_path = os.path.join(output_path, file)
            seq_df = read_sequences(tsv_path)
        elif file.endswith(".names.out"):
            id_df = read_id(os.path.join(output_path, file))  

    if seq_df is None or id_df is None:
        for file in os.listdir(folder):
            if file.endswith("-translated-reference-names.tsv.gz") and id_df is None:
                tsv_path = os.path.join(folder, file)
                tsv_path = gunzip(tsv_path)
                id_df = read_id(tsv_path)
                if os.path.exists(tsv_path):
                    os.remove(tsv_path)
                continue
            elif file.endswith("-translated-protein.fa.gz") and seq_df is None:
                fasta_path = os.path.join(folder, file)
                fasta_path = gunzip(fasta_path)
                seq_df = read_sequences(fasta_path)
                # remove unzipped file
                if os.path.exists(fasta_path):
                    os.remove(fasta_path)

    # Make sure index_df, seq_df, and id_df have all ID column as str
    if index_df is not None:
        index_df["ID"] = index_df["ID"].astype(str)
    if seq_df is not None:
        seq_df["ID"] = seq_df["ID"].astype(str)
    if id_df is not None:
        id_df["ID"] = id_df["ID"].astype(str)   

    if seq_df is not None and id_df is not None:
        merged_df = pd.merge(seq_df, id_df, on="ID", how="left")  
        merged_df = pd.merge(merged_df, index_df, on="ID", how="left") if index_df is not None else merged_df
        if "OneKP_ID" not in merged_df.columns:
            merged_df["OneKP_ID"] = None
    else:
        raise ValueError("Required files not found in the folder.")
    
    merged_df["Organism"] = folder.split("/")[-1]

    # Loop through the rows and if OneKP_ID is NaN, manually parse Id
    for idx, row in merged_df.iterrows():
        if pd.isna(row["OneKP_ID"]):
            original_id = row["ID"]
            organism_name_split = row["Organism"].split("-")
            organism_id = organism_name_split[0].strip()
            organism_name = organism_name_split[1].strip()
            merged_df.at[idx, "OneKP_ID"] = f"scaffold-{organism_id}-{2000001 + int(original_id)}-{organism_name}"
            if idx == 0:
                print("New name", merged_df.at[idx, "OneKP_ID"])

    # Remove the unzipped folder 
    if output_path is not None and os.path.exists(output_path):
        shutil.rmtree(output_path)

    return merged_df

def parse_transcripts(input_folder, output_dir, output_file):
    """
    Parses the protein transcipts from the OneKP dataset. This code was made for the 
    transcript data downloaded from the Google Drive link provided by OneKP project.
    The original Google Drive link is https://drive.google.com/drive/u/0/folders/18AOvneP_1l5uzE7tVWPKVR9MAkhrtA2N.
    This shared folder was first copied to a private Google Drive folder using Colab notebook and this
    code uses the copied folder as input to parse the transcripts.

    Args:
        input_folder (str): Path to the input folder containing OneKP transcript data.
        output_dir (str): Path to the output directory to save individual TSV files.
        output_file (str): Path to the combined output TSV file.

    Returns:
        None
    """

    # Make output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Process each folder in the input directory
    for folder in os.listdir(input_folder):
        folder_path = os.path.join(input_folder, folder)
        if os.path.isdir(folder_path):
            try:
                if os.path.exists(os.path.join(output_dir, f"{folder}.tsv")):
                    continue
                temp_df = process_folder(folder_path)
                temp_df.to_csv(os.path.join(output_dir, f"{folder}.tsv"), sep="\t", index=False)
            except Exception as e:
                print(f"Error processing folder {folder_path.split('/')[-1]}. Error message: {e}")

    # Combine all individual TSV files into one and save it as oneKP_all_transcripts.tsv
    if not os.path.exists(output_file):
        all_dfs = []
        for file in os.listdir(output_dir):
            if file.endswith(".tsv"):
                file_path = os.path.join(output_dir, file)
                temp_df = pd.read_csv(file_path, sep="\t")
                all_dfs.append(temp_df)

        combined_df = pd.concat(all_dfs, ignore_index=True)
        combined_df = combined_df.drop(columns=["ID"])
        combined_df = combined_df.rename(columns={"OneKP_ID": "Accession"})
        combined_df.to_csv(output_file, sep="\t", index=False)
    else:
        print(f"Combined output file {output_file} already exists. Skipping combining step.")

if __name__ == "__main__":

    # Get directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Define input and output paths
    output_dir = f"{script_dir}/parsed_transcripts"
    output_file = f"{script_dir}/oneKP_all_transcripts.tsv"
    google_folder = "/Users/pqh443/Library/CloudStorage/GoogleDrive-ronkko.teemu@gmail.com/My Drive/onekp-data_copied/"

    # Run the parsing function
    parse_transcripts(google_folder, output_dir, output_file)