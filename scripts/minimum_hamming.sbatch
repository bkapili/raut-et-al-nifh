#! /bin/bash
#
#SBATCH --job-name=min_ham
#
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-10

Rscript minimum_hamming.R ${SLURM_ARRAY_TASK_ID}