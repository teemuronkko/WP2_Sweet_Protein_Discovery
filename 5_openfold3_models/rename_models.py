import os 
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_path = script_dir.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"

if __name__ == "__main__":
    model_dir = f"{repo_path}/5_openfold3_models"
    dbs = ["UniProt", "OneKP"]
    for db in dbs:
        model_path = f"{model_dir}/{db}/predictions/"
        if not os.path.exists(model_path):
            continue

        output_dir = f"{repo_path}/5_openfold3_models/{db}/models"
        os.makedirs(output_dir, exist_ok=True)
        batches = [f for f in os.listdir(model_path) if os.path.isdir(os.path.join(model_path, f))]
        for batch in batches:
            batch_path = os.path.join(model_path, batch)
            models = [f for f in os.listdir(batch_path) if os.path.isdir(os.path.join(batch_path, f))]
            for model in models:
                model_path = os.path.join(batch_path, model)
                model_path_full = f"{model_path}/{model}.temp/seed_42/"

                new_path = os.path.join(output_dir, model)
                os.makedirs(new_path, exist_ok=True)
                all_files = os.listdir(model_path_full)
                for file in all_files:
                    new_file = file.replace(".temp_seed_42_sample_", "_")
                    src = os.path.join(model_path_full, file)
                    dst = os.path.join(new_path, new_file)
                    print("Moving", src, "to", dst)
                    os.rename(src, dst)
                    
                print("Removing folder", model_path)
                has_files = any(files for _, _, files in os.walk(model_path))
                if not has_files:
                    os.rmdir(model_path)