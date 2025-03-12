#!/bin/bash
#SBATCH --job-name=clustering_job      # Job name
#SBATCH --output=clustering_job.out    # Output file
#SBATCH --error=clustering_job.err     # Error file
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=1              # Number of CPU cores per task
#SBATCH --time=03:00:00                # Time limit (hh:mm:ss)
#SBATCH --mem=4G                       # Memory per node (4GB)
#SBATCH --partition=general            # Partition to submit to (adjust as needed)

# Run your Python script
python3 clustring.py