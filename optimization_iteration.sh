#!/bin/bash
# Run simulations for the first time, N=num of params+1
OPTIZE_JIDS=""

echo "Running optimization_iteration.sh"

#for N in {0..$(eval echo $1)}
for ((N=0;N<$1;N++))
do
	echo "Submitting initial simulation job ${N}"
	/bin/bash optimize.sh ${N}
	#OPTIZE_JIDS="$OPTIZE_JIDS:$(sbatch "optimize.sh" "${N}")"
	#each simulation is separated into 10 array jobs, saving 20 blocks per job
	#save an array of OPTIZE_JIDS, so next job will wait until all jobs finish to start
done

echo "Submitted initial simulation jobs"

PARAM_JID=$(sbatch --wait "iteration.sh")

echo "Submitted optimization job"

# Run the script to generate parameters
#PARAM_RESULT=$(sbatch "generate_params.sh")

# This extracts the Job ID for the previous command
#PARAM_JID=${PARAM_RESULT##* }

# Wait until the previous batch script completes before running this
# I think that the --wait flag should pause this script until the job completes
#OPTIIZE_JID=$(sbatch --wait --dependency "afterok:$PARAM_JID" "optimize.sh" "${PARAM_JID}")


