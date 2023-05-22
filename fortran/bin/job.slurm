#!/bin/bash
#SBATCH --job-name=KL
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:59:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=youremail@yourhost.edu
#SBATCH --output=KL_%A_%a.out
#SBATCH --error=KL_%A_%a.err
#SBATCH --mem=64GB
module load intel openmpi

ulimit -s unlimited

export MKL_NUM_THREADS=1
export NAG_KUSARI_FILE=lic.txt

./main.exe $SLURM_ARRAY_TASK_ID



