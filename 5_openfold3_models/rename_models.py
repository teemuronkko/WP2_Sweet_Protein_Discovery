import os 
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
            best_model = os.path.basename(json_path).replace('_confidences_aggregated.json', '.pdb')
    
    return best_model

if __name__ == "__main__":

    # Get directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_name = "WP2_SweetProteinDiscovery"
    repo_dir = script_dir.split(repo_name)[0] + repo_name

    model_dir = f"{repo_dir}/5_openfold3_models"
    dbs = ["UniProt", "OneKP"]

    for db in dbs:
        model_path = f"{model_dir}/{db}/predictions/"
        if not os.path.exists(model_path):
            continue

        output_dir = f"{repo_dir}/5_openfold3_models/{db}/models"
        os.makedirs(output_dir, exist_ok=True)

        # List batches in model_path
        batches = [f for f in os.listdir(model_path) if os.path.isdir(os.path.join(model_path, f))]
        for batch in batches:
            batch_path = os.path.join(model_path, batch)
            models = [f for f in os.listdir(batch_path) if os.path.isdir(os.path.join(batch_path, f))]
            for model in models:
                model_path = os.path.join(batch_path, model)
                model_path_full = f"{model_path}/{model}.temp/seed_42/"

                new_path = os.path.join(output_dir, model)
                os.makedirs(new_path, exist_ok=True)

                # List all files in model_path_full
                all_files = os.listdir(model_path_full)
                for file in all_files:
                    new_file = file.replace(".temp_seed_42_sample_", "_")
                    src = os.path.join(model_path_full, file)
                    dst = os.path.join(new_path, new_file)
                    print("Moving", src, "to", dst)
                    os.rename(src, dst)
                    
                # Remove empty folder
                print("Removing folder", model_path)
                # List all files in model_path to check if empty
                has_files = any(files for _, _, files in os.walk(model_path))
                if not has_files:
                    print("Folder is empty, removing:", model_path)
                    os.rmdir(model_path)
                else:
                    print("Folder is not empty, not removing:", model_path)



    # f"{repo_dir}/5_openfold3_models/{db}/predictions/{batch}/{model)/{model}.temp/seed_42/{model}.temp_seed_42_sample_5_model.cif

if False:
    json_files = [f for f in os.listdir(model_path_full) if f.endswith('_confidences_aggregated.json')]
    json_paths = [os.path.join(model_path_full, f) for f in json_files]
    best_model = get_best_model(json_paths)
    for json_path in json_paths:
        print(f"Processed {json_path}")

    if best_model is None:
        print(f"No valid models found for {model} in {db} {batch}.")
       # continue
    else:
        print(f"Best model for {model} in {db} {batch} is {best_model}.")