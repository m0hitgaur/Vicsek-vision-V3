#!/bin/bash

# Get current directory
current_dir=$(pwd)

# Define arrays
noise_folders=("Noise_1" "Noise_0.05" "Noise_0.5")
angle_folders=("Angle_45" "Angle_90" "Angle_120" "Angle_180")

# Loop through arrays and run commands in parallel
for angle_folder in "${angle_folders[@]}"; do
  for noise_folder in "${noise_folders[@]}"; do
    folder_path="$current_dir/$angle_folder/$noise_folder"

    if [ -f "$folder_path/localdensity.cpp" ]; then
      echo "Processing folder: $folder_path"
      # Run commands in background with '&'
      (
        cd "$folder_path"
        g++ localdensity.cpp -o density.out
        if [ $? -ne 0 ]; then
          echo "Error: Compilation failed in $folder_path"
        else
          if [ -f density.out ]; then
            ./density.out
            if [ $? -ne 0 ]; then
              echo "Error: Program execution failed in $folder_path"
            fi
          else
            echo "Error: density.out not found after compilation in $folder_path"
          fi
        fi
        cd "$current_dir"
      ) &
    else
      echo "localdensity.cpp not found in $angle_folder/$noise_folder"
    fi
  done
done

# Wait for all background processes to finish
wait