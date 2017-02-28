#!/bin/bash
#SBATCH -J jobname            # name of job
#SBATCH -t 00:10:00           # upper limit of elapsed time
#SBATCH -p normal             # partition name
#SBATCH -N 1                  # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH -n 2                  # number of total MPI processes, set to SLURM_NTASKS
#SBATCH --ntasks-per-socket=1 # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET
#SBATCH --get-user-env        # retrieve the login environment variables

# set execution file and input argument(s)
EXEC=build/hello
OPTION=

# set socket ID for numactl
PROCS_PER_NODE=`expr $SLURM_NTASKS / $SLURM_JOB_NUM_NODES`
TEMPID=`expr $SLURM_PROCID % $PROCS_PER_NODE` # MPI rank of the current process is $SLURM_PROCID ($*_COMM_WORLD_LOCAL_RANK: MV2 for MVAPICH2, OMPI for OpenMPI)
SOCKET=`expr $TEMPID / $SLURM_NTASKS_PER_SOCKET`

# set number of OpenMP threads per MPI process
export OMP_NUM_THREADS=`expr $SLURM_CPUS_ON_NODE / $PROCS_PER_NODE`

# start logging
cd $SLURM_SUBMIT_DIR
TIME=`date`
echo "start: $TIME"

# set stdout and stderr
STDOUT=log/$SLURM_JOB_NAME.$SLURM_JOB_ID.out
STDERR=log/$SLURM_JOB_NAME.$SLURM_JOB_ID.err

# execute the job
if [ `which numactl` ]; then
    # mpirun with numactl
    echo "mpirun -np $SLURM_NTASKS numactl --cpunodebind=$SOCKET --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    mpirun -np $SLURM_NTASKS numactl --cpunodebind=$SOCKET --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    # mpirun without numactl
    echo "mpirun -np $SLURM_NTASKS $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    mpirun -np $SLURM_NTASKS $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi

# finish logging
TIME=`date`
echo "finish: $TIME"
