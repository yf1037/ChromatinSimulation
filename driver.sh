#!/bin/bash
#SBATCH --job-name=optimization_driver  # Job name
#SBATCH --ntasks=1 # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb # Job memory request
#SBATCH --time=20-00:00:00 # Time limit hrs:min:sec
#SBATCH -p cpu_long


# Load all the packages you need
module load slurm/current
module load anaconda3
source activate py36

# run the python driver script
python driver.py
