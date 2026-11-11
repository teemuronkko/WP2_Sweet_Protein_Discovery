import pandas as pd 
import os

if __name__ == "__main__":
    
    # Load 0_de_novo_candidates.tsv and 0_foldseek_candidates.tsv
    repo_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    repo_path = repo_path.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"
    de_novo_candidates = pd.read_csv(os.path.join(repo_path, "8_protein_selection", "0_de_novo_candidates.tsv"), sep = "\t")
    foldseek_candidates = pd.read_csv(os.path.join(repo_path, "8_protein_selection", "2_round", "0_foldseek_candidates.tsv"), sep = "\t")
    foldseek_candidates = foldseek_candidates.drop(columns = ["length_diff", "interface_residues", "grantham_score", "model", "target", "seq_length"])
    foldseek_candidates = foldseek_candidates.rename(columns = {"cropped_tseq": "sequence"})

    # Make fasta file for all MNEI and Honey_truffle_sweetener hits
    with open(os.path.join(repo_path, "8_protein_selection", "2_round", "0_MNEI_1IV7_hits.fasta"), "w") as fasta_file:
        for index, row in foldseek_candidates.iterrows():
            if row["query"] == "MNEI_1IV7":
                fasta_file.write(f">{row['target_id']}\n{row['sequence']}\n")
        for index, row in de_novo_candidates.iterrows():
            fasta_file.write(f">{row['id']}\n{row['sequence']}\n")
        fasta_file.write(">MNEI_1IV7\nGEWEIIDIGPFTQNLGKFAVDEENKIGQYGRLTFNKVIRPCMKKTIYENEGFREIKGYEYQLYVYASDKLFRADISEDYKTRGRKLLRFNGPVPPP\n")

    with open(os.path.join(repo_path, "8_protein_selection", "2_round", "0_Honey_truffle_sweetener_hits.fasta"), "w") as fasta_file:
        for index, row in foldseek_candidates.iterrows():
            if row["query"] == "Honey_truffle_sweetener":
                fasta_file.write(f">{row['target_id']}\n{row['sequence']}\n")
        fasta_file.write(">Honey_truffle_sweetener\nMPDLSSFITIKNNSNHVFTRTAIYSKYAAVQWSPEPQLSISPGKWDLFILKDILSIRGTSGYVQYRVGDGPGWVRVTFSSLVGADEVAEWSSGDLPDGFVLQKPVRTGSRPLQATFEATKQ\n")

    de_novo_candidates["query"] = "MNEI_1IV7"
    de_novo_candidates["source"] = "ProteinMPNN"
    de_novo_candidates["pdb"] = de_novo_candidates["id"]  
    de_novo_candidates["best_model"] = de_novo_candidates["id"] + "_seed_42_sample_" + de_novo_candidates["model"].astype(str) + "_model"
    de_novo_candidates = de_novo_candidates.drop(columns = ["model"])
    de_novo_candidates["disulfide_bridges"] = [[] for i in range(len(de_novo_candidates))]
    de_novo_candidates["all_cysteines_paired"] = True
    de_novo_candidates = de_novo_candidates.rename(columns = {
        'id': 'target_id',
        'residue_similarity_scores': 'similarity_scores'
    })
    de_novo_candidates = de_novo_candidates[foldseek_candidates.columns] 
    all_candidates = pd.concat([foldseek_candidates, de_novo_candidates], ignore_index = True)
    all_candidates = all_candidates[all_candidates["query"].isin(["MNEI_1IV7", "Honey_truffle_sweetener"])]

    mnei_row = {
        "target_id": "MNEI_1IV7",
        "query": "MNEI_1IV7_native",
        "sequence": "GEWEIIDIGPFTQNLGKFAVDEENKIGQYGRLTFNKVIRPCMKKTIYENEGFREIKGYEYQLYVYASDKLFRADISEDYKTRGRKLLRFNGPVPPP",
        "source": "Native_MNEI_1IV7",
        "pdb": "MNEI_1IV7",
        "best_model": "MNEI_1IV7",
    }

    hts_row = {
        "target_id": "Honey_truffle_sweetener",
        "query": "Honey_truffle_sweetener_native",
        "sequence": "MPDLSSFITIKNNSNHVFTRTAIYSKYAAVQWSPEPQLSISPGKWDLFILKDILSIRGTSGYVQYRVGDGPGWVRVTFSSLVGADEVAEWSSGDLPDGFVLQKPVRTGSRPLQATFEATKQ",
        "source": "Native_Honey_truffle_sweetener",
        "pdb": "Honey_truffle_sweetener",
        "best_model": "Honey_truffle_sweetener",
       
    }

    all_candidates = pd.concat([all_candidates, pd.DataFrame([mnei_row, hts_row])], ignore_index=True)
    all_candidates.to_csv(os.path.join(repo_path, "8_protein_selection", "2_round", "1_merged_candidates.tsv"), sep = "\t", index = False)