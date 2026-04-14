import os
import pandas as pd
from tqdm import tqdm
import json

def get_best_model(json_paths):
    """
    Gets the best model from a list of JSON files based on the highest ranking score.
    Args:
        json_paths (list): List of paths to JSON files containing model rankings.
    Returns:
        str: Name of the best model.
    """
    best_model = None
    best_score = float('-inf')

    for json_path in json_paths:
        with open(json_path, 'r') as f:
            data = json.load(f)
        score = data.get('sample_ranking_score', float('-inf'))
        if score > best_score:
            best_score = score
            best_model = os.path.basename(json_path).replace('_confidences_aggregated.json', '_model.pdb')
    
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

def crop_pdb(pdb_file, start_residue, end_residue, output_file, chain_id="A"):
    """
    Crops a PDB file to include only residues from start_residue to end_residue.

    Args:
        pdb_file (str): Path to input PDB file
        start_residue (int): Starting residue number (1-based)
        end_residue (int): Ending residue number (1-based)
        output_file (str): Path to output cropped PDB file
        chain_id (str): Chain identifier to crop. Default is "A".

    Returns:
        str: Path to cropped PDB file
    """

    from Bio.PDB import PDBParser, PDBIO, Select

    # Check if file exists
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file {pdb_file} does not exist.")
        return None

    class ResidueSelect(Select):
        def accept_residue(self, residue):
            res_id = residue.get_id()[1]
            return start_residue <= res_id <= end_residue

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    if chain_id:
        structure = structure[0][chain_id]

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, ResidueSelect())
    return output_file

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
    

if __name__ == "__main__":

    # Get directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_name = "WP2_SweetProteinDiscovery"
    repo_dir = script_dir.split(repo_name)[0] + repo_name
    outdir = f"{repo_dir}/6_structural_analyses/1_tm_scores"
    structure_outdir = f"{repo_dir}/6_structural_analyses/1_tm_scores/cropped_structures"
    uniprot_outdir = f"{outdir}/UniProt"
    uniprot_cropped = f"{structure_outdir}/UniProt/"
    onekp_outdir = f"{outdir}/OneKP"
    onekp_cropped = f"{structure_outdir}/OneKP/"
    google_drive = "/Users/pqh443/Library/CloudStorage/GoogleDrive-ronkko.teemu@gmail.com/My Drive/WP2_SweetProteinDiscovery/"
    os.makedirs(uniprot_outdir, exist_ok=True)
    os.makedirs(onekp_outdir, exist_ok=True)
    os.makedirs(uniprot_cropped, exist_ok=True)
    os.makedirs(onekp_cropped, exist_ok=True)

    if not os.path.exists(f"{outdir}/uniprot_hits_afdb_tm_scores.tsv"):
        # Read in UniProt hits from Foldseek search results
        fs_results_df = pd.read_csv(f"{repo_dir}/3_foldseek_search/foldseek_UniProt_search.tsv", sep="\t", header = 0)
        fs_results_df["tm_score"] = None
        fs_results_df["RMSD"] = None
        fs_results_df["tm_align_start"] = None
        fs_results_df["tm_align_end"] = None
        fs_results_df["avg_pLDDT"] = None

        # Loop over the structures
        for index, row in tqdm(fs_results_df.iterrows(), total=len(fs_results_df)):
            pdb1 = f"{repo_dir}/1_viridiplantae_prot_db/4_sweet_proteins/cleaned_pdbs/{row['query']}.pdb"
            target = row["target"].split("|")[1]
            pdb2 = f"{repo_dir}/4_afdb_structures/uniprot_hits_afdb_structures/AF-{target}-F1-model_v6.pdb"
            pdb1_name = pdb1.split("/")[-1].replace(".pdb", "")
            pdb2_name = pdb2.split("/")[-1].replace(".pdb", "")
            output_file = f"{uniprot_outdir}/" + pdb1_name + "____" + target + ".txt"
            tm_score, rmsd, start, end = calculate_tm_score(pdb1, pdb2, output_file, norm_length=1)
            fs_results_df.at[index, "tm_score"] = tm_score
            fs_results_df.at[index, "RMSD"] = rmsd
            fs_results_df.at[index, "tm_align_start"] = start
            fs_results_df.at[index, "tm_align_end"] = end
            cropped = crop_pdb(pdb2, start+1, end+1, f"{uniprot_cropped}/cropped_afdb_{pdb1_name}___{target}.pdb")
            fs_results_df.at[index, "avg_pLDDT"] = avg_plddt(cropped, mode = "CA")

        # Save the results
        fs_results_df.to_csv(f"{outdir}/uniprot_hits_afdb_tm_scores.tsv", sep="\t", index=False)

    if not os.path.exists(f"{outdir}/uniprot_hits_missing_tm_scores.tsv"):
        # AF3: Read in UniProt hits from Foldseek search results
        fs_results_df = pd.read_csv(f"{repo_dir}/3_foldseek_search/foldseek_UniProt_search.tsv", sep="\t", header = 0)
        fs_results_df["af3_tm_score"] = None
        fs_results_df["of3_tm_score"] = None
        fs_results_df["af3_RMSD"] = None
        fs_results_df["of3_RMSD"] = None
        fs_results_df["af3_tm_align_start"] = None
        fs_results_df["af3_tm_align_end"] = None
        fs_results_df["of3_tm_align_start"] = None
        fs_results_df["of3_tm_align_end"] = None
        fs_results_df["af3_avg_pLDDT"] = None
        fs_results_df["of3_avg_pLDDT"] = None

        # Read in the UniProt hits that already have AFDB TM-scores and drop them from fs_results_df
        uniprot_hits_w_structures = pd.read_csv(f"{outdir}/uniprot_hits_afdb_tm_scores.tsv", sep="\t", header=0)
        uniprot_hits_w_structures = uniprot_hits_w_structures[~uniprot_hits_w_structures["tm_score"].isna()]

        # Drop these rows from fs_results_df
        fs_results_df = fs_results_df[~fs_results_df["target"].isin(uniprot_hits_w_structures["target"])]
        fs_results_df = fs_results_df.reset_index(drop=True)

        # Loop over the structures
        for index, row in tqdm(fs_results_df.iterrows(), total=len(fs_results_df)):

            # AF3 TM-score calculation
            pdb1 = f"{repo_dir}/1_viridiplantae_prot_db/4_sweet_proteins/cleaned_pdbs/{row['query']}.pdb"
            target = row["target"].split("|")[1]
            pdb2 = f"{google_drive}/5_alphafold3_models/UniProt/models/{target.lower()}/{target.lower()}_model.pdb"
            if os.path.exists(pdb2):
                pdb1_name = pdb1.split("/")[-1].replace(".pdb", "")
                pdb2_name = pdb2.split("/")[-1].replace(".pdb", "")
                output_file = f"{uniprot_outdir}/AF3_" + pdb1_name + "____" + target + ".txt"
                tm_score, rmsd, start, end = calculate_tm_score(pdb1, pdb2, output_file, norm_length=1)
                fs_results_df.at[index, "af3_tm_score"] = tm_score
                fs_results_df.at[index, "af3_RMSD"] = rmsd
                fs_results_df.at[index, "af3_tm_align_start"] = start
                fs_results_df.at[index, "af3_tm_align_end"] = end
                cropped = crop_pdb(pdb2, start+1, end+1, f"{uniprot_cropped}/cropped_af3_{pdb1_name}___{target}.pdb")
                fs_results_df.at[index, "af3_avg_pLDDT"] = avg_plddt(cropped, mode = "CA")
            else:
                print(pdb2, "does not exist, skipping", target)

            # OF3 TM-score calculation
            json_folder = f"{google_drive}/5_openfold3_models/UniProt/{target}/{target}.temp/seed_42"
            if not os.path.exists(json_folder):
                print(f"JSON folder {json_folder} does not exist, skipping {target}")
                continue
            json_paths = [os.path.join(json_folder, f) for f in os.listdir(json_folder) if f.endswith("_confidences_aggregated.json")]
            best_model = get_best_model(json_paths)
            pdb2 = f"{json_folder}/{best_model}"
            if os.path.exists(pdb2):
                output_file = f"{uniprot_outdir}/OF3_" + pdb1_name + "____" + target + ".txt"
                tm_score, rmsd, start, end = calculate_tm_score(pdb1, pdb2, output_file, norm_length=1)
                fs_results_df.at[index, "of3_tm_score"] = tm_score
                fs_results_df.at[index, "of3_RMSD"] = rmsd
                fs_results_df.at[index, "of3_tm_align_start"] = start
                fs_results_df.at[index, "of3_tm_align_end"] = end
                cropped = crop_pdb(pdb2, start+1, end+1, f"{uniprot_cropped}/cropped_of3_{pdb1_name}___{target}.pdb")
                fs_results_df.at[index, "of3_avg_pLDDT"] = avg_plddt(cropped, mode = "CA") 
            else:
                print(pdb2, " does not exist, skipping ", target)

        fs_results_df.to_csv(f"{outdir}/uniprot_hits_missing_tm_scores.tsv", sep="\t", index=False)


    if not os.path.exists(f"{outdir}/onekp_hits_tm_scores.tsv"):
        # AF3: Read in OneKP hits from Foldseek search results
        fs_results_df = pd.read_csv(f"{repo_dir}/3_foldseek_search/foldseek_OneKP_search.tsv", sep="\t", header = 0)
        fs_results_df["af3_tm_score"] = None
        fs_results_df["of3_tm_score"] = None
        fs_results_df["af3_RMSD"] = None
        fs_results_df["of3_RMSD"] = None
        fs_results_df["af3_tm_align_start"] = None
        fs_results_df["af3_tm_align_end"] = None
        fs_results_df["of3_tm_align_start"] = None
        fs_results_df["of3_tm_align_end"] = None
        fs_results_df["af3_avg_pLDDT"] = None
        fs_results_df["of3_avg_pLDDT"] = None

        # Loop over the structures
        for index, row in tqdm(fs_results_df.iterrows(), total=len(fs_results_df)):

            # AF3 TM-score calculation
            pdb1 = f"{repo_dir}/1_viridiplantae_prot_db/4_sweet_proteins/cleaned_pdbs/{row['query']}.pdb"
            target = row["target"]
            target_cleaned = target.replace("+", "_").lower()
            pdb2 = f"{google_drive}/5_alphafold3_models/OneKP/models/{target_cleaned}/{target_cleaned}_model.pdb"
            if os.path.exists(pdb2):
                pdb1_name = pdb1.split("/")[-1].replace(".pdb", "")
                pdb2_name = pdb2.split("/")[-1].replace(".pdb", "")
                output_file = f"{onekp_outdir}/AF3_" + pdb1_name + "____" + target + ".txt"
                tm_score, rmsd, start, end = calculate_tm_score(pdb1, pdb2, output_file, norm_length=1)
                fs_results_df.at[index, "af3_tm_score"] = tm_score
                fs_results_df.at[index, "af3_RMSD"] = rmsd
                fs_results_df.at[index, "af3_tm_align_start"] = start
                fs_results_df.at[index, "af3_tm_align_end"] = end
                cropped = crop_pdb(pdb2, start+1, end+1, f"{onekp_cropped}/cropped_af3_{pdb1_name}___{target}.pdb")
                fs_results_df.at[index, "af3_avg_pLDDT"] = avg_plddt(cropped, mode = "CA")
            else:
                print(pdb2, "does not exist, skipping", target)

            # OF3 TM-score calculation
            target_cleaned = target.replace("+", "_")
            json_folder = f"{google_drive}/5_openfold3_models/OneKP/{target_cleaned}/{target_cleaned}.temp/seed_42"
            if not os.path.exists(json_folder):
                print(f"JSON folder {json_folder} does not exist, skipping {target}")
                continue
            json_paths = [os.path.join(json_folder, f) for f in os.listdir(json_folder) if f.endswith("_confidences_aggregated.json")]
            best_model = get_best_model(json_paths)
            pdb2 = f"{json_folder}/{best_model}"
            if os.path.exists(pdb2):
                output_file = f"{onekp_outdir}/OF3_" + pdb1_name + "____" + target + ".txt"
                tm_score, rmsd, start, end = calculate_tm_score(pdb1, pdb2, output_file, norm_length=1)
                fs_results_df.at[index, "of3_tm_score"] = tm_score
                fs_results_df.at[index, "of3_RMSD"] = rmsd
                fs_results_df.at[index, "of3_tm_align_start"] = start
                fs_results_df.at[index, "of3_tm_align_end"] = end
                cropped = crop_pdb(pdb2, start+1, end+1, f"{onekp_cropped}/cropped_of3_{pdb1_name}___{target}.pdb")
                fs_results_df.at[index, "of3_avg_pLDDT"] = avg_plddt(cropped, mode = "CA")  
            else:
                print(pdb2, " does not exist, skipping ", target)

        fs_results_df.to_csv(f"{outdir}/onekp_hits_tm_scores.tsv", sep="\t", index=False)