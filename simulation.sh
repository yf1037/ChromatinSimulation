#!/bin/bash
#SBATCH --array 1-10
#SBATCH --ntasks=1 # Run on a single CPU
#SBATCH --time=1:00:00 # Time limit hrs:min:sec
#SBATCH --job-name=SimulationHiCMYC
#SBATCH --partition=gpu8_short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=50GB
#SBATCH --mail-type=END
#SBATCH --mail-user=Yi.Fu2@nyumc.org
#SBATCH --mail-user=yf1037@nyu.edu
#SBATCH --output=MYC_%j.out

echo "Running simulation.sh"

module purge
module load openmm/7.4.0
module load cuda/9.0

PARAM_ID=$1

echo "PARAM_ID ${PARAM_ID}"
echo "ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID}"

python SimulatorCode.py "${PARAM_ID}" "${SLURM_ARRAY_TASK_ID}"
#for testing
#module load anaconda3
#source activate py36

#python test.py "${PARAM_ID}" "${SLURM_ARRAY_TASK_ID}"
