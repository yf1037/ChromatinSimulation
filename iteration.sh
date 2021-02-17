#!/bin/bash
#SBATCH --job-name=iteration
#SBATCH --partition=cpu_long
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=20-00:00:00
#SBATCH --mem=4GB
#SBATCH --mail-type=END
#SBATCH --mail-user=Yi.Fu2@nyumc.org
#SBATCH --mail-user=yf1037@nyu.edu
#SBATCH --output=MYC_%j.out

echo "Running iteration.sh"

# Load all the packages you need
module load slurm/current
module load anaconda3
source activate py36

# run the python driver script
python iteration.py
