#!/bin/bash

# Author: Maryam KafiKang
# Created: 2024
# Description: This script runs TNet model on FAVITES output data using SLURM job scheduler

#SBATCH --job-name=TNet
#SBATCH --output=logs/tnet_%A_%a.out
#SBATCH --error=logs/tnet_%A_%a.out
#SBATCH --nodes=4
#SBATCH --array=1-500%100
#SBATCH --partition=lo-core
#SBATCH --mem=32G

# Create logs directory if not exists
mkdir -p logs

i=$SLURM_ARRAY_TASK_ID
input_dir="Favites_log"

# Define full path to the test directory
full_path="$input_dir/FAVITES_output_logSI_contemp_T200_N100_E1_$i"
# Construct the input file path
fav_input="${full_path}/error_free_files/phylogenetic_trees/Tnet_Tree.0"
# Define the output file path
output_file_path="tnet_output_updated/FAVITES_output_logSI_contemp_T200_N100_E1_${i}_TNet_Network.txt"

# Extract the directory path from the output file path
output_dir=$(dirname "$output_file_path")

# Create the directory if it doesn't exist
mkdir -p "$output_dir"

# Check if output file already exists
if [[ -f "$output_file_path" ]]; then
    echo "Skipping: Output file already exists for node $i -> $output_file_path"
    exit 0
fi

# Check if input file exists before running the python script
if [[ -f "$fav_input" ]]; then
    start_time=$(date +%s.%N)

    # Run TNet model
    python tnet.py "$fav_input" "$output_file_path" -t 1000 -rs -mx

    # End time
    end_time=$(date +%s.%N)

    # Compute elapsed time in seconds
    runtime=$(echo "$end_time - $start_time" | bc)

    # Append runtime to the output file
    echo -e "\nRunning Time (seconds): $runtime" >> "$output_file_path"

else
    echo "Input file not found: $fav_input"
    exit 1
fi
