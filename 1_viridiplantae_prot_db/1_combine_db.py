from Bio import SeqIO
import json
import os
import gzip
import hashlib

def gzip_file(input_path: str, gzip_path: str = None) -> str:
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")

    if gzip_path is None:
        gzip_path = input_path + ".gz"

    try:
        with open(input_path, 'rb') as f_in:
            with gzip.open(gzip_path, 'wb') as f_out:
                f_out.writelines(f_in)
        return gzip_path
    except Exception as e:
        if os.path.exists(gzip_path):
            os.remove(gzip_path)
        raise

def md5(file_path: str, block_size=2**20) -> str:
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    hasher = hashlib.md5()
    try:
        with open(file_path, 'rb') as f:
            while True:
                chunk = f.read(block_size)
                if not chunk:
                    break
                hasher.update(chunk)
        md5 = hasher.hexdigest()
    
    except Exception as e:
        raise e

    md5_file_path = file_path + ".md5"
    try:
        with open(md5_file_path, 'w') as md5_file:
            md5_file.write(md5)
    except Exception as e:
        raise e

    return md5_file_path

def read_sequences(fasta_path, db_id: str = None, max_seq_length: int = None):
    """
    Reads a fasta file containing protein sequences.
    
    Args:
        fasta_path (str): Path to the fasta file.
        db_id (str): Optional database identifier to prepend to sequence IDs.
        max_seq_length (int): Optional maximum sequence length to filter sequences.

    Returns:
        tuple: Tuple of 'ID' and 'Sequence
    """
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        if len(record.id) == 0 or len(record.seq) == 0:
            print(f"Skipping entry with empty ID or sequence: ID='{record.id}', Sequence length={len(record.seq)}")
            continue
        elif (max_seq_length is not None) and (len(record.seq) > max_seq_length):
            continue
        if db_id is not None:
            seq_id = f"{db_id}_{record.id}"
        else:
            seq_id = record.id

        sequences.append((seq_id, str(record.seq)))
    return sequences

def split_seq(seq, accession, window, step):
    """
    Splits an input sequence into overlapping fragments given a specified window size and step size.

    Args:
        seq (str): The input sequence to be split.
        accession (str): The accession identifier for the sequence.
        window (int): The size of each fragment window.
        step (int): The step size for overlapping fragments.

    Returns:
        dict: A dictionary where keys are fragment identifiers (accession_start-end) and values are the corresponding sequence fragments.
    """

    # Empty dictionary to hold fragments
    fragments = {}

    # If the sequence is shorter than or equal to the window size, return it as a single fragment
    if len(seq) <= window:
        return {accession: seq}
    else:
        # Loop through the sequence to create overlapping fragments
        for i in range(0, len(seq), step):

            # Determine the end index of the fragment
            if i + window > len(seq):
                end = len(seq)
            else:
                end = i + window

            # Extract the fragment and add it to the dictionary if it is longer than window - step
            frag = seq[i:end]
            if len(frag) > window-step:
                accession_frag = f"{accession}_frag_{i+1}-{i+len(frag)}"
                fragments[accession_frag] = frag
    return fragments

def combine_dbs(
        uniprot_viridiplantae_path: str,
        oneKP_all_transcripts_path: str,
        toxin_uniprot_path: str,
        out_db_path: str,
        json_file: str,
        max_seq_length: int,
        window: int,
        step: int,
        deduplicate: bool = False,
        compress: bool = False):
    """
    Combines UniProt Viridiplantae and OneKP transcript sequences into a single fasta file,
    while tracking duplicate sequences and their associated IDs.

    Args:
        uniprot_viridiplantae_path (str): Path to the UniProt Viridiplantae fasta file.
        oneKP_all_transcripts_path (str): Path to the OneKP transcripts tsv file.
        toxin_uniprot_path (str): Path to the UniProt toxins fasta file.
        out_db_path (str): Path to output combined fasta file.
        json_file (str): Path to output JSON file mapping sequences to their IDs.
        max_seq_length (int): Maximum sequence length.
        window (int): Window size for overlapping fragments when splitting long sequences.
        step (int): Step size for overlapping fragments when splitting long sequences.
        deduplicate (bool): Whether to deduplicate sequences.
        compress (bool): Whether the fasta file should be compressed.

    Returns:
        None; Saves the combined fasta file and JSON mapping file.
    """

    # Read oneKP_all_transcripts.tsv line by line and save "Sequence" and "Accession" to a list of tuples
    print("Reading OneKP transcript sequences...")
    seq_tuples = []
    with open(oneKP_all_transcripts_path, 'r') as file:
        next(file)  # Skip header line
        for line in file:
            line = line.strip().split('\t')
            accession = line[3]
            accession = "OneKP_" + accession
            sequence = line[0]
            if len(sequence) == 0 or len(accession) == 0 or not isinstance(sequence, str):
                print(f"Skipping entry with empty sequence or accession: Accession='{accession}', Sequence length={len(sequence)}")
                continue
            elif (max_seq_length is not None) and (len(sequence) > max_seq_length):
                continue
            
            seq_tuples.append((accession, sequence))
    print(f"Total OneKP transcript sequences read: {len(seq_tuples)}\n")

    # Read UniProt Viridiplantae fasta file into tuplesgrep
    print("Reading UniProt Viridiplantae sequences...")
    uniprot_tuples = read_sequences(uniprot_viridiplantae_path, db_id="UniProt", max_seq_length=max_seq_length)
    print(f"Total UniProt Viridiplantae sequences read: {len(uniprot_tuples)}\n")

    # Read toxin UniProt fasta file
    print("Reading UniProt toxin sequences...")
    toxin_tuples = read_sequences(toxin_uniprot_path, db_id="Toxin", max_seq_length=max_seq_length)
    print(f"Total UniProt toxin sequences read: {len(toxin_tuples)}\n")

    # Combine the datasets
    seq_tuples = toxin_tuples + uniprot_tuples + seq_tuples
    print(f"Total combined sequences before processing: {len(seq_tuples)}\n")

    # Remove uniprot_tuples and toxin_tuples to free up memory
    del uniprot_tuples
    del toxin_tuples

    if deduplicate:
        # Add "_dedup" suffix to output file name
        out_db_path = out_db_path.replace(".fasta", "_dedup.fasta")

        # Create a mapping of sequences to their accessions as a dictionary
        print("Mapping sequences to accessions...")
        seq_map = {}
        for accession, sequence in seq_tuples:
            cleaned_sequence = sequence.strip()
            if cleaned_sequence in seq_map:
                seq_map[cleaned_sequence].append(accession)
            else:
                seq_map[cleaned_sequence] = [accession]
        print(f"Total unique sequences after mapping: {len(seq_map)}")

        # Open fasta file for writing sequences
        print("Writing combined fasta database and JSON mapping file...")
        alt_keys = {}
        with open(out_db_path, 'w') as fasta_handle:
            # Loop through the sequence map to write data
            for sequence, id_list in seq_map.items():
                
                # Determine the accession and Alternative_IDs
                toxin_id = next((id for id in id_list if id.startswith("Toxin_")), None)
                if toxin_id:
                    id_list.remove(toxin_id)
                    accession = toxin_id
                    if len(id_list) > 0 and not all(id.startswith("Toxin_") for id in id_list):
                        print(f"Warning: Toxin ID {toxin_id} has alternative IDs: {id_list}")
                else:
                    accession = id_list.pop(0) if id_list else None

                # Continue if sequence is longer than max_seq_length
                if len(sequence) > max_seq_length:
                    continue

                # If there are alternative IDs, add them to alt_keys
                if len(id_list) > 0:
                    alt_keys[accession] = id_list

                # Split long sequences into overlapping fragments
                if window is not None and step is not None:
                    split_sequences = split_seq(sequence, accession, window=window, step=step)
                    for frag_accession, frag_sequence in split_sequences.items():
                        fasta_header = f">{frag_accession}\n"
                        fasta_handle.write(fasta_header)
                        fasta_handle.write(f"{frag_sequence}\n")
                else:
                    # Write fasta entry
                    fasta_header = f">{accession}\n"
                    fasta_handle.write(fasta_header)
                    fasta_handle.write(f"{sequence}\n")

        # Free up memory and remove seq_map
        del seq_map

        # Write seq_map into json file
        with open(json_file, 'w') as json_handle:
            json.dump(alt_keys, json_handle)

    else:
        # Write all sequences directly to fasta file without deduplication
        print("Writing combined fasta database without deduplication...")
        with open(out_db_path, 'w') as fasta_handle:
            for accession, sequence in seq_tuples:
                if len(sequence) > max_seq_length:
                    continue
                fasta_header = f">{accession}\n"
                fasta_handle.write(fasta_header)
                fasta_handle.write(f"{sequence}\n")

    # Make separate fasta files for OneKP, UniProt and Toxin
    with open(f"{out_db_path.replace('.fasta', '_OneKP.fasta')}", 'w') as oneKP_handle, \
         open(f"{out_db_path.replace('.fasta', '_UniProt.fasta')}", 'w') as uniprot_handle, \
         open(f"{out_db_path.replace('.fasta', '_Toxin.fasta')}", 'w') as toxin_handle: 
        with open(out_db_path, 'r') as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                if record.id.startswith("OneKP_"):
                    SeqIO.write(record, oneKP_handle, "fasta")
                elif record.id.startswith("UniProt_"):
                    SeqIO.write(record, uniprot_handle, "fasta")
                elif record.id.startswith("Toxin_"):
                    SeqIO.write(record, toxin_handle, "fasta")
 
    if compress:
        # Generate md5 checksum for the fasta file and compress the fasta db file
        print("Generating md5 checksum and compressing the fasta database...")
        md5_path = md5(out_db_path)
        gzip_path = gzip_file(out_db_path)

        # Unzip gzip file to check file integrity and generate md5 checksum
        print("Verifying gzip integrity by unzipping and generating md5 checksum...")
        fasta_out = out_db_path.replace(".fasta", "_unzipped.fasta")
        with gzip.open(gzip_path, 'rt') as gz_handle, open(fasta_out, 'w') as out_handle:
            for line in gz_handle:
                out_handle.write(line)
        md5_unzipped_path = md5(fasta_out)

        print("MD5 Checksum Results:")
        print(f"MD5 checksum for {out_db_path}: {md5_path}")
        with open(md5_path, 'r') as f:
            md5_checksum = f.read().strip()
        print(f"MD5 checksum value: {md5_checksum}")
        print(f"MD5 checksum for {fasta_out}: {md5_unzipped_path}")
        with open(md5_unzipped_path, 'r') as f:
            md5_checksum_unzipped = f.read().strip()
        print(f"MD5 checksum value: {md5_checksum_unzipped}")

        if md5_checksum == md5_checksum_unzipped:
            print("MD5 checksums match between gzipped and unzipped files.")
        else:
            print("MD5 checksums do NOT match between gzipped and unzipped files.")
            os.remove(gzip_path)

        os.remove(fasta_out)
        os.remove(md5_unzipped_path)
    return 


if __name__ == "__main__":

    # Get directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_name = "WP2_SweetProteinDiscovery"
    repo_dir = script_dir.split(repo_name)[0] + repo_name

    # Paths to OneKP and UniProt Viridiplantae fasta files
    oneKP_all_transcripts_path = f"{repo_dir}/1_viridiplantae_prot_db/1_oneKP_capstone_2019/oneKP_all_transcripts.tsv"
    uniprot_viridiplantae_path = f"{repo_dir}/1_viridiplantae_prot_db/2_uniprot_data/uniprot_viridiplantae_all.fasta"
    uniprot_toxins_path = f"{repo_dir}/1_viridiplantae_prot_db/2_uniprot_data/uniprot_toxins.fasta"

    # Do not deduplicate sequences
    deduplicate = False

    # Define output filenames
    os.makedirs(f"{repo_dir}/1_viridiplantae_prot_db/", exist_ok=True)
    out_db_path = f"{repo_dir}/1_viridiplantae_prot_db/viridiplantae_prot_db.fasta"
    json_file = f"{repo_dir}/1_viridiplantae_prot_db/duplicate_sequences.json"

    # Set maximum sequence length
    max_sequence_length = 1000
    window = None
    step = None

    print("Combining UniProt Viridiplantae and OneKP transcript sequences into a single fasta database...")
    print(f"UniProt Viridiplantae path: {uniprot_viridiplantae_path}")
    print(f"OneKP all transcripts path: {oneKP_all_transcripts_path}")
    print(f"UniProt toxins path: {uniprot_toxins_path}")
    print(f"Output database path: {out_db_path}")
    print(f"Output JSON mapping file: {json_file}")
    combine_dbs(
        uniprot_viridiplantae_path, 
        oneKP_all_transcripts_path,
        uniprot_toxins_path,
        out_db_path,
        json_file,
        max_seq_length=max_sequence_length,
        window=window,
        step=step,
        deduplicate=deduplicate,
        compress = False
    )