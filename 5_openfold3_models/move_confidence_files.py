import os
import shutil

def move_confidence_files(source_dir, destination_dir, pattern="_confidences.json"):
    """
    Move confidence files from source directory to destination directory.

    Args:
        source_dir (str): The directory to search for confidence files.
        destination_dir (str): The directory to move confidence files to.
        pattern (str): The file pattern to match (default is "_confidences.json").
    """
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)
        print(f"Created folder: {destination_dir}")

    count = 0
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            if file.endswith(pattern):
                source_path = os.path.join(root, file)
                dest_path = os.path.join(destination_dir, file)

                if os.path.exists(dest_path):
                    name, ext = os.path.splitext(file)
                    dest_path = os.path.join(destination_dir, f"{name}_{count}{ext}")

                shutil.move(source_path, dest_path)
                count += 1

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_path = script_dir.split("WP2_Sweet_Protein_Discovery")[0] + "/WP2_Sweet_Protein_Discovery"
    search_in = f"{repo_path}/5_openfold3_models/onekp_hits"
    move_to = f"{repo_path}/5_openfold3_models/confidence_files/onekp_hits"
    move_confidence_files(search_in, move_to, pattern="_confidences.json")