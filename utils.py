import os
import gemmi
import pandas as pd
import sys
import numpy as np

# Custom imports from utils.py
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_path = script_dir.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"
sys.path.append(repo_path)

# These residues are retrived from Table 3 of the following paper:
# Vo, Phillip, Daniel E. Connors, Brendan P. Sharkey, Stephen A. Gravina, Guy Servant, Chase T. McFarland, Tyler T. Pitkanen et al. 
# "Sweet Protein Allosteric Binding and Activation of the Human T1R2/R3 Sweet Receptor: A Simulation Model Validated by in Vitro Receptor Activation Assay."
# Biochemistry (2026).
# For Miraculin, the important histidine residues were chosen based on the following paper:
# Koizumi, Ayako, Asami Tsuchiya, Ken-ichiro Nakajima, Keisuke Ito, Tohru Terada, Akiko Shimizu-Ibuka, Loïc Briand, Tomiko Asakura, Takumi Misaka, and Keiko Abe. 
# "Human sweet taste receptor mediates acid-induced sweetness of miraculin."
# Proceedings of the National Academy of Sciences 108, no. 40 (2011): 16819-16824.
residues_to_check = {
    "Thaumatin-I_1LR2": [("K", 66), ("K", 77), ("K", 96), ("K", 105), ("K", 136), ("Y", 168), ("K", 186), ("D", 20)],
    "Thaumatin-II_3AOK": [("R", 66), ("K", 77), ("R", 81), ("K", 96), ("K", 105), ("K", 136), ("Y", 168), ("K", 186), ("D", 20)],
    "MNEI_1IV7": [('E', 22), ('Q', 27), ("C", 40), ("Y", 64), ("S", 66), ("D", 67), ("E", 1)],
    "Brazzein_4HE7": [("D", 1), ("K", 4), ("Y", 38), ("D", 39), ("E", 40), ("K", 41), ("R", 42), ("D", 49), ("Y", 53), ("D", 28), ("E", 35), ("R", 32)],
    "Honey_truffle_sweetener": [("R", 19), ("T", 20), ("K", 25), ("Q", 30), ("R", 56), ("Y", 61), ("Q", 63), ("R", 65), ("R", 74), ("S", 78), ("E", 88)],
    "Miraculin": [("H", 29), ("H", 59)]
}

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
    import shutil
    import gzip

    if output_path is None:
        output_path = path[:-3] 

    with gzip.open(path, 'rb') as input_file:
        with open(output_path, 'wb') as output_file:
            shutil.copyfileobj(input_file, output_file)

    if remove_original and os.path.exists(output_path):
        os.remove(path)

    return output_path

def gzip_file(input_path: str, gzip_path: str = None) -> str:
    """
    Compresses a file using gzip.

    Args:
        input_path (str): Path to the input file to be compressed.
        gzip_path (str, optional): Path to save the compressed file. If None, it will be saved as input_path + ".gz".

    Returns:
        str: Path to the compressed gzip file.
    """

    import gzip
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
        raise e

def md5(file_path: str, block_size=2**20) -> str:
    """
    Calculates the MD5 checksum of a file and saves it to a .md5 file.

    Args:
        file_path (str): Path to the input file.
        block_size (int): Size of the blocks to read from the file. Default is 1 MB.

    Returns:
        str: Path to the generated .md5 file containing the checksum.
    """
    import hashlib
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

def read_sequences(fasta_path):
    """
    Reads a fasta file containing protein sequences.
    
    Args:
        fasta_path (str): Path to the fasta file.

    Returns:
        pd.DataFrame: DataFrame with columns 'ID', 'Sequence'.
    """
    from Bio import SeqIO
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences.append({"ID": record.id, "Sequence": str(record.seq)})

    return pd.DataFrame(sequences)

def read_sequence_tuples(fasta_path, db_id: str = None, max_seq_length: int = None):
    """
    Reads a fasta file containing protein sequences.
    
    Args:
        fasta_path (str): Path to the fasta file.
        db_id (str): Optional database identifier to prepend to sequence IDs.
        max_seq_length (int): Optional maximum sequence length to filter sequences.

    Returns:
        tuple: Tuple of 'ID' and 'Sequence
    """
    from Bio import SeqIO
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

def crop_pdb(pdb: str, output_file: str, chain_id: str, start_idx: int, end_idx: int) -> str:
    """
    Crops a PDB file to include a specific range and chain, given the start and end residue indices.
    Args:
        pdb (str): Path to the input PDB file.
        output_file (str): Path to save the cropped PDB file.
        chain_id (str): Chain identifier to retain.
        start_idx (int): Starting residue index (inclusive).
        end_idx (int): Ending residue index (inclusive).

    Returns:
        str: Path to the cropped PDB file.
    """
    from Bio.PDB import PDBParser, PDBIO
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("structure", pdb)
    for model in structure:
        for chain in list(model):
            if chain.id != chain_id:
                model.detach_child(chain.id)
                continue

            for residue in list(chain):
                if start_idx is None and end_idx is None:
                    continue
                elif not residue.id[1] in range(start_idx, end_idx + 1):
                    chain.detach_child(residue.id)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

    if os.path.exists(output_file):
        print(f"Cropped PDB saved to {output_file}")
        return output_file
    else:
        print("Error: Cropped PDB file was not created.")
        return None
    

def cif_to_pdb(cif_file: str, pdb_file: str) -> str:
    """
    Converts a CIF file to PDB format.
    Args:
        cif_file (str): Path to the input CIF file.
        pdb_file (str): Path to save the converted PDB file.

    Returns:
        str: Path to the converted PDB file.
    """
    from Bio.PDB import MMCIFParser, PDBIO
    parser = MMCIFParser(QUIET=True)   # QUIET avoids warnings
    structure = parser.get_structure("structure", cif_file)

    # Write PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

    if os.path.exists(pdb_file):
        print(f"Converted PDB saved to {pdb_file}")
        return pdb_file
    else:
        print("Error: Converted PDB file was not created.")
        return None
    
def cif_to_pdb_gemmi(cif_file, pdb_file):
    """
    Converts an OpenFold3 CIF file to PDB format using Gemmi.
    Handles missing occupancy and multi-character chain IDs natively.

    Args:
        cif_file (str): Path to the input CIF file.
        pdb_file (str): Path to save the converted PDB file.

    Returns:
        str: Path to the converted PDB file.
    """
    try:
        structure = gemmi.read_structure(cif_file)
        structure.write_pdb(pdb_file)
    except Exception as e:
        print(f"Error converting {cif_file}: {e}")
        return None
    if os.path.exists(pdb_file):
        print(f"Converted PDB saved to {pdb_file}")
        return pdb_file
    
def get_afdb_structure(uid: str, output_dir: str, version: int = 6) -> str:
    """
    Downloads a protein structure from the AlphaFold Database.
    Args:
        uid (str): UniProt ID of the protein.
        output_dir (str): Directory to save the downloaded PDB file.
        version (int): Version of the AlphaFold model to download. Default is 6.
    Returns:
        str: Path to the downloaded PDB file.
    """
    import requests
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v{version}.pdb"
    if output_dir is None:
        output_dir = "."
    os.makedirs(output_dir, exist_ok=True)
    
    output_file = os.path.join(output_dir, f"AF-{uid}-F1-model_v{version}.pdb")   
    if os.path.exists(output_file):
        print(f"AlphaFold structure for {uid} already exists at {output_file}")
        return output_file

    response = requests.get(url)
    response.raise_for_status()
    with open(output_file, "wb") as f:
        f.write(response.content)

    if os.path.exists(output_file):
        print(f"Downloaded AlphaFold structure for {uid} to {output_file}")
        return output_file
    else:
        print("Error: PDB file was not created.")
        return None
    
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

def check_canonical_aa(sequence):
    """
    Check if a sequence contains only canonical amino acids.

    Args:
        sequence (str): Amino acid sequence.
    Returns:
        bool: True if sequence contains only canonical amino acids, False otherwise.
    """
    AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
    return all(residue in AMINO_ACIDS for residue in sequence)

def get_best_of3_model(json_paths):
    """
    Gets the best model from a list of JSON files based on the highest ranking score.
    Args:
        json_paths (list): List of paths to JSON files containing model rankings.
    Returns:
        str: Name of the best model.
    """
    import json
    best_model = None
    best_score = float('-inf')

    for json_path in json_paths:
        with open(json_path, 'r') as f:
            data = json.load(f)
        score = data.get('sample_ranking_score', float('-inf'))
        if score > best_score:
            best_score = score
            best_model = os.path.basename(json_path).replace('_confidences_aggregated.json', '.pdb')
    
    return best_model

def pdb_len(pdb_file, chain_id="A"):
    """
    Calculates the length of a protein from its PDB file.

    Args:
        pdb_file (str): Path to PDB file
    Returns:
        int: Length of the protein (number of residues)
    """
    from Bio.PDB import PDBParser, is_aa
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    chain = structure[0][chain_id]
    return sum(1 for res in chain if is_aa(res, standard=True))

def avg_plddt(pdb_file, chain_id="A", mode = "CA"):
    """
    Calculates the average pLDDT score for a given chain in a PDB file.

    Args:
        pdb_file (str): Path to PDB file
        chain_id (str): Chain identifier to calculate pLDDT for. Default is "A".
        mode (str): Mode of calculation. Default is "CA" (use CA atoms' B-factors).
                    Otherwise use all atoms.

    Returns:
        float: Average pLDDT score for the specified chain.
    """
    from Bio.PDB import PDBParser

    # Check if file exists
    if pdb_file is None or not os.path.exists(pdb_file):
        print(f"Error: PDB file {pdb_file} does not exist.")
        return None

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    chain = structure[0][chain_id]

    plddt_scores = []
    for residue in chain:
        for atom in residue:
            if atom.get_name() == "CA" and mode == "CA":
                b_factor = atom.get_bfactor()
                plddt_scores.append(b_factor)
                break
            elif mode == "all_atoms":
                b_factor = atom.get_bfactor()
                plddt_scores.append(b_factor)

    if plddt_scores:
        return sum(plddt_scores) / len(plddt_scores)
    else:
        return 0.0
    
def calculate_tm_score(pdb_file_1, pdb_file_2, outfile, norm_length = 0):
    """
    Outputs the TM-score between two proteins given their PDB files.

    Args:
        pdb_file_1 (str): Path to first PDB file
        pdb_file_2 (str): Path to second PDB file
        outfile (str): Path to output file to save TMalign results
        norm_length (int, optional): Normalization length for TM-score calculation. 
                                     If 0, no normalization is applied,
                                     If 1, pdb_file_1 length is used, 
                                     and if 2, pdb_file_2 length is used. Defaults to 0.
    Returns:
        float: TM-score between the two proteins (normalized as specified).
    """
    import subprocess

    command = ["TMalign", pdb_file_1, pdb_file_2]

    if norm_length == 1:
        len_1 = pdb_len(pdb_file_1)
        command.extend(["-u", str(len_1)])
    elif norm_length == 2:
        len_2 = pdb_len(pdb_file_2)
        command.extend(["-u", str(len_2)])

    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running TMalign:", result.stderr)
        return None, None
    
    # Parse TM-score from the output
    tm_score = None
    rmsd = None
    line_c = 0
    start = -1
    end = -1
    for line in result.stdout.splitlines():
        line_c += 1
        if norm_length != 0 and line.startswith("TM-score=") and "user-specified LN" in line:
            tm_score = float(line.split()[1])

        if line.startswith("Aligned length="):
            rmsd = line.split(",")[1].split("=")[1].strip()
            rmsd = float(rmsd)

        if line_c == 21:
            # Find position of first "." or ":"
            start_a = line.find(":")
            start_b = line.find(".")
            if start_a == -1 or (start_b != -1 and start_b < start_a):
                start = start_b
            else:
                start = start_a
            
            # Find position of last "." or ":"
            end_a = line.rfind(":")
            end_b = line.rfind(".")
            if end_a == -1 or (end_b != -1 and end_b > end_a):
                end = end_b
            else:
                end = end_a

        elif line_c == 22 and start != -1 and end != -1:
            start -= line[:start].count("-")
            end -= line[:end+1].count("-")

    # Save output to file
    with open(outfile, "w") as f:
        f.write(result.stdout)

    return tm_score, rmsd, start, end

def calculate_grantham_score(res1, res2, matrix_csv_path):
    """
    Calculates the Grantham score between two lists of amino acid residues.

    Args:
        res1 (list): The first list of amino acid residues (single-letter codes).
        res2 (list): The second list of amino acid residues (single-letter codes).
        matrix_csv_path (str): Path to the CSV file containing the Grantham matrix.
    
    Returns:
        float: The average Grantham score between the two set residues.
    """
    grantham_matrix = pd.read_csv(matrix_csv_path, index_col=0, sep = ";")
    grantham_dict = {}
    for r1 in grantham_matrix.index:
        for r2 in grantham_matrix.columns:
            value = grantham_matrix.loc[r1, r2]
            if pd.isna(value):
                value = grantham_matrix.loc[r2, r1]
            grantham_dict[(r1, r2)] = value
            grantham_dict[(r2, r1)] = value  
            grantham_dict[(r2, r2)] = 0  
        grantham_dict[(r1, r1)] = 0  
            
    scores = []
    for r1, r2 in zip(res1, res2):
        if (r1, r2) in grantham_dict:
            scores.append(grantham_dict[(r1, r2)])
        else:
            raise ValueError(f"Residues {r1} and/or {r2} not found in Grantham matrix.")

    return sum(scores) / len(scores)

def read_aligned_residues(tm_align_res_file):
    """
    Reads aligned residues from a TM-align result file.

    Args:
        tm_align_res_file (str): Path to the TM-align result file.

    Returns:
        dict: A dictionary where keys are positions in the query sequence and values are dictionaries containing the aligned residues from the query and target sequences.
    """
    aligned_residues = {}
    with open(tm_align_res_file, "r") as f:
        lines = f.readlines()
        query_seq = lines[19]
        align_info = lines[20]
        target_seq = lines[21]

    position = -1
    for index, residue in enumerate(query_seq):
        if query_seq[index] != "-":
             position += 1
        if align_info[index] == ":":
            aligned_residues[position] = {"query": residue, "target": target_seq[index]}
    
    return aligned_residues


def make_of3_json(chain_dict, model_name, num_recycles=3, json_path="of3_input.json", use_paired_msas=False, use_msas=True):
    """
    Makes an input json file for OF3 given a chain_dict.

    Args:
        chain_dict (dict): A dictionary where keys are chain IDs and values are dictionaries with keys
            "seq" (str): The sequence of the chain.
            "msa_path" (str): The path to the MSA file for the chain.
        model_name (str): The name of the model to use in OF3.
        num_recycles (int, optional): The number of recycles to use in OF3. Defaults to 3.
        json_path (str, optional): The path to save the output json file. Defaults to "of3_input.json".
        use_paired_msas (bool, optional): Whether to use paired MSAs in OF3. Defaults to False.
        use_msas (bool, optional): Whether to use MSAs in OF3. Defaults to True.


    Returns:
        None: The function saves the input json file to the specified path.
    """
    import json
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    json_dict = {"queries": {model_name: {}}}
    for chain in chain_dict.keys():
        if "chains" not in json_dict["queries"][model_name]:
            json_dict["queries"][model_name]["chains"] = []

        if "msa_path" not in chain_dict[chain]:
            json_dict["queries"][model_name]["chains"].append({
                "molecule_type": "protein",
                "chain_ids": [chain],
                "sequence": chain_dict[chain]["seq"],
            })
        else:
            json_dict["queries"][model_name]["chains"].append({
                "molecule_type": "protein",
                "chain_ids": [chain],
                "sequence": chain_dict[chain]["seq"],
                "main_msa_file_paths": [chain_dict[chain]["msa_path"]],
            })
    
    json_dict["queries"][model_name]["model"] = {
        "num_recycles": num_recycles
    }
    if not use_paired_msas:
        json_dict["queries"][model_name]["use_paired_msas"] = use_paired_msas

    if not use_msas:
        json_dict["queries"][model_name]["use_msas"] = use_msas

    with open(json_path, 'w') as f:
        json.dump(json_dict, f, indent=4)

def get_contacts(pdb, chain1_id, chain2_id, cutoff=5.0):
    """
    Get the number of contacts between two chains using Biopython's NeighborSearch.

    Args:
        pdb (str): Path to the PDB file.
        chain1_id (str): Identifier for the first chain (e.g., "A").
        chain2_id (str): Identifier for the second chain (e.g., "B").
        cutoff (float): Distance cutoff in angstroms to consider atoms as being in contact.
    
    Returns:
        list: List of unique residue pairs in contact, where each pair is represented as a tuple (res1_id, res1_name, res2_id, res2_name).
    """
    from Bio.PDB import PDBParser, NeighborSearch

    # Parser to read the PDB file and extract the specified chains and their atoms
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb)
    model = structure[0]
    
    # Read the specified chains and extract their atoms into flat lists
    chain1 = model[chain1_id]
    chain2 = model[chain2_id]
    atoms_chain1 = list(chain1.get_atoms())
    atoms_chain2 = list(chain2.get_atoms())
    
    # Loop through atoms in chain 1 and find close atoms in chain 2 using NeighborSearch
    ns = NeighborSearch(atoms_chain2)
    contact_set = set()    
    for atom1 in atoms_chain1:
        close_atoms = ns.search(atom1.coord, cutoff)
        for atom2 in close_atoms:
            res1 = atom1.get_parent()
            res2 = atom2.get_parent()
            contact_set.add((
                res1.get_id()[1], res1.get_resname(),
                res2.get_id()[1], res2.get_resname()
            ))
    return list(contact_set)

def calculate_distance(atom1, atom2):
    """
    Helper function to calculate the Euclidean distance between two atoms.
    atom1: numpy array with the coordinates of the first atom
    atom2: numpy array with the coordinates of the second atom

    Returns the Euclidean distance between the two atoms.
    """
    return np.linalg.norm(atom1 - atom2)

def find_disulfide_pairs(pdb_file, min_distance=0.0, max_distance=3.0, chain = None, atom = "SG", min_gap = 3, use_actual_numbering = False):
    """
    Find disulfide bridges in a protein structure given a PDB file.
    pdb_file: path to the PDB file
    min_distance: minimum distance between sulfur atoms to form a disulfide bond
    max_distance: maximum distance between sulfur atoms to form a disulfide bond
    chain: chain identifier of the protein
    atom: atom type for cysteine sulfur atoms
    min_gap: minimum gap between cysteines to form a disulfide bond

    return: a dictionary with chain identifiers as keys 
    and a list of tuples with the residue numbers of cysteines forming disulfide bonds.
    """

    # Dictionary to store positions of cysteine residues
    cysteine_positions = {}
    
    # Open PDB file and extract sulfur atoms of cysteines
    with open(pdb_file, 'r') as f:
        for line in f:
            # Get sulfur atoms of cysteines
            if line.startswith('ATOM') and line[17:20] == 'CYS' and line[12:16].strip() == atom:
                # Extract chain, residue number, and atom coordinates
                chain_id = line[21]
                if chain and line[21] != chain:
                    continue

                # Get coordinates of the sulfur atom
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                # Get residue number
                residue_number = int(line[22:26].strip())

                # Add coordinates to dictionary
                if chain_id not in cysteine_positions:
                    cysteine_positions[chain_id] = []
                cysteine_positions[chain_id].append((residue_number, np.array([x, y, z])))

    # If dictionary is empty, return an empty dict
    if chain and chain not in cysteine_positions:
        return {chain: []}
    if not cysteine_positions:
        return {}

    # Get numbering of cysteines residues, starting from 1 for each chain
    cysteine_order = {}
    for chain, positions in cysteine_positions.items():
        c = 1
        for position in positions:
            if chain not in cysteine_order:
                cysteine_order[chain] = {}
            if use_actual_numbering:
                cysteine_order[chain][position[0]] = position[0]
            else:
                cysteine_order[chain][position[0]] = c 
            c += 1

    # Identify pairs of closest cysteines forming disulfide bridges
    disulfide_pairs = {}
    for chain, positions in cysteine_positions.items():

        # If there are more than one cysteine in the chain, find the closest pairs
        if len(positions) > 1:
            # List to store the distances between all cysteine pairs
            distances = []

            # Calculate the distance between all cysteine pairs
            for i in range(len(positions)):
                for j in range(i + 1, len(positions)):
                    # Make sure that neighboring cysteines are not paired
                    if abs(positions[i][0] - positions[j][0]) < min_gap:
                        continue
                    distance = calculate_distance(positions[i][1], positions[j][1])
                    if min_distance <= distance and distance <= max_distance:
                        distances.append((distance, positions[i][0], positions[j][0]))

            # Sort the distances to find the closest pairs
            distances.sort()
            used_residues = set()
            disulfide_pairs[chain] = []

            # Greedily select the closest pairs without reusing cysteines
            for dist, res1, res2 in distances:
                if res1 not in used_residues and res2 not in used_residues:
                    disulfide_pairs[chain].append((res1, res2))
                    used_residues.add(res1)
                    used_residues.add(res2)

    # Sort the disulfide pairs by the first residue number
    for chain, pairs in disulfide_pairs.items():
        disulfide_pairs[chain] = sorted(pairs, key=lambda x: x[0])
        disulfide_pairs[chain] = [(cysteine_order[chain][res1], cysteine_order[chain][res2]) for res1, res2 in disulfide_pairs[chain]]
                   
    return disulfide_pairs

def sequence_identity(seq1, seq2):
    """
    Simple function to calculate sequence identity between two sequences.
    Assumes sequences are of the same length and are aligned. If sequences are of different lengths, it calculates identity based on the length of the longer sequence.
    
    Args:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
    
    Returns:
        float: Sequence identity between the two sequences.
    """
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    length = max(len(seq1), len(seq2))
    identity = matches / length if length > 0 else 0
    return identity