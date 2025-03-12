#!/bin/bash
#SBATCH --job-name=cas.txt
#SBATCH --output=logs/cas_%A_%a.out
#SBATCH --error=logs/cas_%A_%a.out
#SBATCH --nodes=1
#SBATCH --array=1-500%20   
#SBATCH --partition=lo-core
#SBATCH --mem=32G

mkdir -p logs

TEST_NUM=$SLURM_ARRAY_TASK_ID

# Run R script with arguments
python run_cassiopeia.py $TEST_NUM