import os
import shutil

def move_confidence_files(source_dir, destination_dir):
    # Create the destination directory if it doesn't exist
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)
        print(f"Created folder: {destination_dir}")

    count = 0
    # Walk through the directory and subdirectories
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            if file.endswith("_confidences.json") and "summary_confidences.json" not in file:
                source_path = os.path.join(root, file)
                dest_path = os.path.join(destination_dir, file)

                # Handle potential filename collisions
                if os.path.exists(dest_path):
                    name, ext = os.path.splitext(file)
                    dest_path = os.path.join(destination_dir, f"{name}_{count}{ext}")

                shutil.move(source_path, dest_path)
                #print(f"Moved: {file} to {dest_path}")
                count += 1

    print(f"\nTask complete. Total files moved: {count}")

# Define your paths here
search_in = "/Users/pqh443/Documents/Git_projects/WP2_SweetProteinDiscovery/5_alphafold3_models/UniProt"
move_to = "/Users/pqh443/Documents/Git_projects/WP2_SweetProteinDiscovery/5_alphafold3_models/confidence_files/UniProt"

move_confidence_files(search_in, move_to)