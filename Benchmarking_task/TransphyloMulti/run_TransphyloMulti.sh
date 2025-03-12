#!/bin/bash
#SBATCH --job-name=transM.txt
#SBATCH --output=logs/transM_%A_%a.out
#SBATCH --error=logs/transM_%A_%a.out
#SBATCH --array=1-500%20   
#SBATCH --partition=lo-core
#SBATCH --mem=32G

mkdir -p logs

SHAPE=0.1
TEST_NUM=$SLURM_ARRAY_TASK_ID

# Run R script with arguments
Rscript run_TransphyloMulti.R $SHAPE $TEST_NUM