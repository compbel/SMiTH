#!/bin/bash
#SBATCH --job-name=startus.txt
#SBATCH --output=logs/st_%A_%a.out
#SBATCH --error=logs/st_%A_%a.out
#SBATCH --nodes=1
#SBATCH --array=1-500%50  # SLURM requires increasing arrays
#SBATCH --partition=lo-core
#SBATCH --mem=16G

mkdir -p logs
module load r

# Reverse SLURM task ID to run from 500 to 1
i=$SLURM_ARRAY_TASK_ID

count_hosts() {
    local tree_file="$1"
    # Extract host IDs from leaf nodes (format: number_Nnumber:branchlength)
    # First get all leaf nodes, then extract just the host numbers (everything before _N)
    grep -o '[0-9]*_N[0-9]*:[0-9.]*' "$tree_file" | grep -o '^[0-9]*' | sort -u | wc -l
}
# Set paths to required files
TREE_FILE="log_50/FAVITES_output_logSI_contemp_T200_N100_E1_${i}/error_free_files/phylogenetic_trees/Tnet_Tree.0"
OUTPUT_PREFIX="stratus_output/FAVITES_output_logSI_contemp_T200_N100_E1_${i}/"
MULTIPLE_SAMPLING_FILE="STraTUS/host_data/FAVITES_output_logSI_contemp_T200_N100_E1_${i}.csv"

# Check if the output already exists
OUTPUT_FILE="${OUTPUT_PREFIX}_annotations.csv"
if [[ -f "$OUTPUT_FILE" ]]; then
    echo "_annotations in Output file for job $i already exists. Skipping..."
    exit 0
fi

if [[ -f "$TREE_FILE" ]]; then
    # Count unique hosts
    host_count=$(count_hosts "$TREE_FILE")
    echo "Number of unique hosts in tree: $host_count"
    
    # Check if host count is between 5 and 30
    if [ "$host_count" -ge 5 ] && [ "$host_count" -le 30 ]; then
        # Create output directory if it does not exist
        mkdir -p "$OUTPUT_PREFIX"

        # Define log file for execution time
        TIME_LOG="${OUTPUT_PREFIX}Runtime_$i.txt"

        # Construct the command
        COMMAND="Rscript STraTUSCommandLine.R \
            -m $MULTIPLE_SAMPLING_FILE \
            -s 10000 \
            $TREE_FILE $OUTPUT_PREFIX"

        # Start timer
        start_time=$(date +%s)

        # Run STraTUS and suppress standard output
        eval $COMMAND

        # End timer
        end_time=$(date +%s)

        # Compute elapsed time
        elapsed_time=$((end_time - start_time))

        # Save execution time to log file
        echo "Run $i: $elapsed_time seconds" >> $TIME_LOG

        echo "STraTUS execution completed for job $i in $elapsed_time seconds."
    else
        echo "Skipping job $i: Number of hosts ($host_count) is not between 5 and 30"
    fi
else
    echo "Tree file not found for job $i"
fi
