#!/bin/sh

# measure execution time of MAGI as a function of Ntot repeatedly
NITER=10

PROBLEM=80

for NTOT in 131072 262144 524288 1048576 2097152 4194304 8388608 16777216 33554432 67108864 134217728 268435456 536870912
do
    for (( ii = 0; ii < $NITER; ii++ ))
    do
	sbatch --dependency=singleton --export=PROBLEM=${PROBLEM},NTOT=${NTOT} sh/slurm/magi_exec.sh
    done
done

exit 0
