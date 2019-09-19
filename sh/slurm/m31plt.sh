#!/bin/bash
#SBATCH -J m31plt      # name of job
#SBATCH -t 08:00:00    # upper limit of elapsed time
#SBATCH -p normal      # partition name
#SBATCH --nodes=1      # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=16    # number of total MPI processes, set to SLURM_NTASKS
#SBATCH --get-user-env # retrieve the login environment variables

EXEC=py/shownws.py
# EXEC=py/showgap.py

SPECIFY="--specify"

TARGET="--files=nws-test-m95-orbit4"
# TARGET="--files=nws-test-m95-orbit5"
# TARGET="--files=nws-continue"

PDFFIGS="--pdf"
CONTINUE="--continued"

# use numactl if available
if [ `which numactl` ]; then
    NUMACTL="numactl --localalloc"
fi

# set number of MPI processes per node
PROCS_PER_NODE=`expr $SLURM_NTASKS / $SLURM_JOB_NUM_NODES`

# set number of MPI processes per socket
SOCKETS_PER_NODE=`lscpu | grep NUMA | grep CPU | wc -l`
PROCS_PER_SOCKET=`expr $PROCS_PER_NODE / $SOCKETS_PER_NODE`
if [ $PROCS_PER_SOCKET -lt 1 ]; then
    PROCS_PER_SOCKET=1
fi

# visualize outputs
MPL_CFG_DIR=/tmp/matplotlib
mpiexec -n $SLURM_NTASKS sh/mplwrapper.sh $PROCS_PER_NODE $PROCS_PER_SOCKET $MPL_CFG_DIR $EXEC $TARGET $PDFFIGS $CONTINUE $SPECIFY

exit 0
