#!/bin/bash

# Get current directory
current_dir=$(pwd)



# Check if 'files' directory exists
files_dir="$current_dir/files"
if [ ! -d "$files_dir" ]; then
  echo "Error: 'files' directory not found in the current directory."
  exit 1
fi

# Define arrays for folders
noise_folders=("Noise_0.05" )
angle_folders=("Angle_45" "Angle_180")

# Create all folders (angle and noise directories)
for angle_folder in "${angle_folders[@]}"; do
  for noise_folder in "${noise_folders[@]}"; do
    folder_path="$current_dir/$angle_folder/$noise_folder"
    
    # Ensure the directories exist (create if missing)
    mkdir -p "$folder_path"
  done
done

# Copy files from the 'files' directory to each noise folder
for angle_folder in "${angle_folders[@]}"; do
  for noise_folder in "${noise_folders[@]}"; do
    folder_path="$current_dir/$angle_folder/$noise_folder"
    
    # Copy all contents from 'files' to the current folder
    cp -r "$files_dir/"* "$folder_path/"
  done
done
