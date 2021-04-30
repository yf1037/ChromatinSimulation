#!/bin/bash
# Run simulations for the first time, N=num of params

echo "Running optimization_iteration.sh"

SIM=$2

#for N in {0..$(eval echo $1)}
for ((N=0;N<$1;N++))
do
	echo "Submitting initial simulation job ${N}"
	/bin/bash optimize.sh ${N} ${SIM} 0
	#OPTIZE_JIDS="$OPTIZE_JIDS:$(sbatch "optimize.sh" "${N}")"
	#each simulation is separated into 10 array jobs, saving 20 blocks per job
	#save an array of OPTIZE_JIDS, so next job will wait until all jobs finish to start
done

/bin/bash optimize.sh ${N} ${SIM} 1

echo "Submitted initial simulation jobs"

sbatch iteration.sh ${SIM}

echo "Submitted optimization job"
