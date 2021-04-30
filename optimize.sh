#!/bin/bash

echo "Running optimize.sh"

PARAM_JID=$1

#submit array jobs 20 dat per job * 10 jobs
OPTIZE_STRING=$(sbatch "simulation.sh" "${PARAM_JID}")
OPTIZE_JID=${OPTIZE_STRING##* }
echo "submitted simulation job: dependant on ${OPTIZE_JID}"

#calculate correlation from 200 dat
CALC_JID=$(sbatch --wait --dependency "afterok:$OPTIZE_JID" "calc.sh" "${PARAM_JID}")

echo "Submited calculation job"
