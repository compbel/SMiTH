#!/bin/bash
#SBATCH --job-name=phyloscanner.txt
#SBATCH --output=logs/ph_%A_%a.out
#SBATCH --error=logs/ph_%A_%a.out
#SBATCH --nodes=1
#SBATCH --array=1-500%20   
#SBATCH --partition=lo-core
#SBATCH --mem=32G

mkdir -p logs
module load r
conda activate phyloscanner-env

TEST_NUM=$SLURM_ARRAY_TASK_ID

# Run R script with arguments
python phyloscanner_run_favlog.py $TEST_NUM