import pandas as pd 
import os

if __name__ == "__main__":
    
    # Load 0_de_novo_candidates.tsv and 0_foldseek_candidates.tsv
    repo_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    repo_path = repo_path.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"
    de_novo_candidates = pd.read_csv(os.path.join(repo_path, "8_protein_selection", "0_de_novo_candidates.tsv"), sep = "\t")
    foldseek_candidates = pd.read_csv(os.path.join(repo_path, "8_protein_selection", "0_foldseek_candidates.tsv"), sep = "\t")
    foldseek_candidates = foldseek_candidates.drop(columns = ["length_diff", "interface_residues", "grantham_score", "model", "target", "seq_length"])
    foldseek_candidates = foldseek_candidates.rename(columns = {"cropped_tseq": "sequence"})

    # First, make fasta file for all MNEI hits
    with open(os.path.join(repo_path, "8_protein_selection", "0_MNEI_1IV7_hits.fasta"), "w") as fasta_file:
        for index, row in foldseek_candidates.iterrows():
            if row["query"] == "MNEI_1IV7":
                fasta_file.write(f">{row['target_id']}\n{row['sequence']}\n")
        for index, row in de_novo_candidates.iterrows():
            fasta_file.write(f">{row['id']}\n{row['sequence']}\n")
        fasta_file.write(">MNEI_1IV7\nGEWEIIDIGPFTQNLGKFAVDEENKIGQYGRLTFNKVIRPCMKKTIYENEGFREIKGYEYQLYVYASDKLFRADISEDYKTRGRKLLRFNGPVPPP\n")

    # Then, make fasta file for all Thaumatin hits
    with open(os.path.join(repo_path, "8_protein_selection", "0_Thaumatin_hits.fasta"), "w") as fasta_file:
        for index, row in foldseek_candidates.iterrows():
            if row["query"] in ["Thaumatin-I_1LR2", "Thaumatin-II_3AOK"]:
                fasta_file.write(f">{row['target_id']}\n{row['sequence']}\n")
        fasta_file.write(">Thaumatin-I_1LR2\nATFEIVNRCSYTVWAAASKGDAALDAGGRQLNSGESWTINVEPGTNGGKIWARTDCYFDDSGSGICKTGDCGGLLRCKRFGRPPTTLAEFSLNQYGKDYIDISNIKGFNVPMNFSPTTRGCRGVRCAADIVGQCPAKLKAPGGGCNDACTVFQTSEYCCTTGKCGPTEYSRFFKRLCPDAFSYVLDKPTTVTCPGSSNYRVTFCPTA\n")
        fasta_file.write(">Thaumatin-II_3AOK\nATFEIVNRCSYTVWAAASKGDAALDAGGRQLNSGESWTINVEPGTKGGKIWARTDCYFDDSGRGICRTGDCGGLLQCKRFGRPPTTLAEFSLNQYGKDYIDISNIKGFNVPMDFSPTTRGCRGVRCAADIVGQCPAKLKAPGGGCNDACTVFQTSEYCCTTGKCGPTEYSRFFKRLCPDAFSYVLDKPTTVTCPGSSNYRVTFCPTA\n")

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
    all_candidates = all_candidates[all_candidates["query"].isin(["MNEI_1IV7", "Thaumatin-I_1LR2", "Thaumatin-II_3AOK"])]

    # Add native MNEI sequence to the dataframe
    mnei_row = {
        "target_id": "MNEI_1IV7",
        "query": "MNEI_1IV7_native",
        "sequence": "GEWEIIDIGPFTQNLGKFAVDEENKIGQYGRLTFNKVIRPCMKKTIYENEGFREIKGYEYQLYVYASDKLFRADISEDYKTRGRKLLRFNGPVPPP",
        "source": "Native_MNEI_1IV7",
        "pdb": "MNEI_1IV7",
        "best_model": "MNEI_1IV7",
    }

    thaumatin_i_row = {
        "target_id": "Thaumatin-I_1LR2",
        "query": "Thaumatin-I_1LR2_native",
        "sequence": "ATFEIVNRCSYTVWAAASKGDAALDAGGRQLNSGESWTINVEPGTNGGKIWARTDCYFDDSGSGICKTGDCGGLLRCKRFGRPPTTLAEFSLNQYGKDYIDISNIKGFNVPMNFSPTTRGCRGVRCAADIVGQCPAKLKAPGGGCNDACTVFQTSEYCCTTGKCGPTEYSRFFKRLCPDAFSYVLDKPTTVTCPGSSNYRVTFCPTA",
        "source": "Native_Thaumatin-I",
        "pdb": "Thaumatin-I_1LR2",
        "best_model": "Thaumatin-I_1LR2",
    }

    thaumatin_ii_row = {
        "target_id": "Thaumatin-II_3AOK",
        "query": "Thaumatin-II_3AOK_native",
        "sequence": "ATFEIVNRCSYTVWAAASKGDAALDAGGRQLNSGESWTINVEPGTKGGKIWARTDCYFDDSGRGICRTGDCGGLLQCKRFGRPPTTLAEFSLNQYGKDYIDISNIKGFNVPMDFSPTTRGCRGVRCAADIVGQCPAKLKAPGGGCNDACTVFQTSEYCCTTGKCGPTEYSRFFKRLCPDAFSYVLDKPTTVTCPGSSNYRVTFCPTA",
        "source": "Native_Thaumatin-II",
        "pdb": "Thaumatin-II_3AOK",
        "best_model": "Thaumatin-II_3AOK",
    }

    all_candidates = pd.concat([all_candidates, pd.DataFrame([mnei_row, thaumatin_i_row, thaumatin_ii_row])], ignore_index=True)
    all_candidates.to_csv(os.path.join(repo_path, "8_protein_selection", "1_merged_candidates.tsv"), sep = "\t", index = False)