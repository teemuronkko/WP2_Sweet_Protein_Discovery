import os
import pandas as pd
from tqdm import tqdm
import json
import sys

# Custom imports from utils.py
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_path = script_dir.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"
sys.path.append(repo_path)
from utils import get_best_of3_model, crop_pdb, avg_plddt, calculate_tm_score, gzip_file, gunzip, pdb_len

if __name__ == "__main__":
    outdir = f"{repo_path}/6_structural_analyses/1_tm_scores"
    structure_outdir = f"{repo_path}/6_structural_analyses/1_tm_scores/cropped_structures"
    uniprot_outdir = f"{outdir}/UniProt"
    uniprot_cropped = f"{structure_outdir}/UniProt/"
    onekp_outdir = f"{outdir}/OneKP"
    onekp_cropped = f"{structure_outdir}/OneKP/"
    google_drive = "~/Library/CloudStorage/GoogleDrive-ronkko.teemu@gmail.com/My Drive/WP2_Sweet_Protein_Discovery/"
    os.makedirs(uniprot_outdir, exist_ok=True)
    os.makedirs(onekp_outdir, exist_ok=True)
    os.makedirs(uniprot_cropped, exist_ok=True)
    os.makedirs(onekp_cropped, exist_ok=True)

    if not os.path.exists(f"{outdir}/uniprot_hits_afdb_tm_scores.tsv"):
        # Read in UniProt hits from Foldseek search results
        fs_results_df = pd.read_csv(f"{repo_path}/3_foldseek_search/foldseek_UniProt_search.tsv", sep="\t", header = 0)
        fs_results_df["tm_score"] = None
        fs_results_df["RMSD"] = None
        fs_results_df["tm_align_start"] = None
        fs_results_df["tm_align_end"] = None
        fs_results_df["avg_pLDDT"] = None

        # Loop over the structures
        for index, row in tqdm(fs_results_df.iterrows(), total=len(fs_results_df)):
            pdb1 = f"{repo_path}/1_viridiplantae_prot_db/4_sweet_proteins/cleaned_pdbs/{row['query']}.pdb"
            target = row["target"].split("|")[1]
            pdb2 = f"{repo_path}/4_afdb_structures/uniprot_hits_afdb_structures/AF-{target}-F1-model_v6.pdb"
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
        fs_results_df = pd.read_csv(f"{repo_path}/3_foldseek_search/foldseek_UniProt_search.tsv", sep="\t", header = 0)
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
            pdb1 = f"{repo_path}/1_viridiplantae_prot_db/4_sweet_proteins/cleaned_pdbs/{row['query']}.pdb"
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
            best_model = get_best_of3_model(json_paths)
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
        fs_results_df = pd.read_csv(f"{repo_path}/3_foldseek_search/foldseek_OneKP_search.tsv", sep="\t", header = 0)
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
            pdb1 = f"{repo_path}/1_viridiplantae_prot_db/4_sweet_proteins/cleaned_pdbs/{row['query']}.pdb"
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
            best_model = get_best_of3_model(json_paths)
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