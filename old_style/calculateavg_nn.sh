#!/bin/bash

# Get current directory
current_dir=$(pwd)

# Define arrays
noise_folders=("Noise_0.05")
angle_folders=("Angle_45"  "Angle_180")
# Loop through arrays and run commands in parallel
for angle_folder in "${angle_folders[@]}"; do
  for noise_folder in "${noise_folders[@]}"; do
    folder_path="$current_dir/$angle_folder/$noise_folder"

    if [ -f "$folder_path/calculateavg_nneighbour.cpp" ]; then
      echo "Processing folder: $folder_path"
      # Run commands in background with '&'
      (
        cd "$folder_path"
        g++ calculateavg_nneighbour.cpp -o nn.out
        if [ $? -ne 0 ]; then
          echo "Error: Compilation failed in $folder_path"
        else
          if [ -f nn.out ]; then
            ./nn.out
            if [ $? -ne 0 ]; then
              echo "Error: Program execution failed in $folder_path"
            fi
          else
            echo "Error: nn.out not found after compilation in $folder_path"
          fi
        fi
        cd "$current_dir"
      ) &
    else
      echo "calculateavg_nneighbour.cpp not found in $angle_folder/$noise_folder"
    fi
  done
done

# Wait for all background processes to finish
wait