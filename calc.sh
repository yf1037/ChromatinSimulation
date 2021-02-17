#!/bin/bash
#SBATCH --job-name=GenerateHiCMYC
#SBATCH --partition=cpu_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:10:00
#SBATCH --mem=30GB
#SBATCH --mail-type=END
#SBATCH --mail-user=Yi.Fu2@nyumc.org
#SBATCH --mail-user=yf1037@nyu.edu
#SBATCH --output=MYC_%j.out

echo "Running calc.sh"

module purge
module load miniconda3
source activate py36

ID=$1

python ContactFrequency3.py "${ID}"

