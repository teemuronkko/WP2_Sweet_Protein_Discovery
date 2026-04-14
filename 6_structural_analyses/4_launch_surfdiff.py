import os
import subprocess

if __name__ == "__main__":

    # Get directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_name = "WP2_SweetProteinDiscovery"
    repo_dir = script_dir.split(repo_name)[0] + repo_name
    bash_script_path = f"{repo_dir}/6_structural_analyses/4_run_surfdiff_docker.sh"
    mode = "batch"

    # Path to sweet protein structure folder
    sweet_protein_dir = f"{repo_dir}/1_viridiplantae_prot_db/4_sweet_proteins/cleaned_pdbs"

    # Define input directories
    input_dir = f"{repo_dir}/6_structural_analyses/3_representative_structures/"
    output_dir = f"{repo_dir}/6_structural_analyses/tmp/"
    output_final = f"{repo_dir}/6_structural_analyses/4_surfdiff_results/"
    os.makedirs(output_final, exist_ok=True)

    # Loop through each sweet protein
    for sweet_protein in os.listdir(sweet_protein_dir):
        if "Miraculin" not in sweet_protein:
            continue
        print(f"Processing sweet protein: {sweet_protein}")
        sweet_protein_name = sweet_protein.replace(".pdb", "")
        sweet_protein_input_dir = f"{input_dir}/{sweet_protein_name}/"
        pdb_dir = f"{input_dir}/{sweet_protein_name}/pdb_dir/"

        # Copy cleaned sweet protein PDB file to pdb_dir
        query_pdb = os.path.join(sweet_protein_dir, sweet_protein)
        query_pdb_copy = os.path.join(pdb_dir, sweet_protein)
        if not os.path.exists(query_pdb_copy):
            os.system(f"cp {query_pdb} {query_pdb_copy}")

        # Make output directory for the sweet protein
        sweet_protein_output_dir = os.path.join(output_final, sweet_protein_name)
        os.makedirs(sweet_protein_output_dir, exist_ok=True)

        query_pdb = query_pdb_copy
        subject_pdbs = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if f != sweet_protein and f.startswith("cropped_") and f.endswith(".pdb") and sweet_protein_name in f]
        print(f"Found {len(subject_pdbs)} subject PDBs for {sweet_protein_name}.")

        # Make sure output_dir exists
        os.makedirs(output_dir, exist_ok=True)

        subject_pdbs = [s.replace(".pdb", "").split("/")[-1] for s in subject_pdbs]
        query_pdb = query_pdb.replace(".pdb", "").split("/")[-1]
        if mode == "single":
            for pdb_1 in subject_pdbs:
                pdb_1 = [pdb_1]
                cmd = [
                    "bash",
                    bash_script_path,
                    sweet_protein_input_dir,
                    output_dir,
                    query_pdb,
                    *pdb_1,
                ]
                subprocess.run(cmd)

                sweet_protein_ind_dir = os.path.join(sweet_protein_output_dir, pdb_1[0])
                os.makedirs(sweet_protein_ind_dir, exist_ok=True)
                os.system(f"mv {output_dir} {sweet_protein_ind_dir}")

        elif mode == "batch":
        
            subject_pdbs = [s.replace(".pdb", "").split("/")[-1] for s in subject_pdbs]
            query_pdb = query_pdb.replace(".pdb", "").split("/")[-1]
            cmd = [
                "bash",
                bash_script_path,
                sweet_protein_input_dir,
                output_dir,
                query_pdb,
                *subject_pdbs,
            ]
            subprocess.run(cmd)
            
            # Move results to output directory
            os.system(f"mv {output_dir} {sweet_protein_output_dir}")