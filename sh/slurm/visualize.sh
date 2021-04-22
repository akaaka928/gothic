#!/bin/bash
#SBATCH -J visualize          # name of job
#SBATCH -t 02:00:00           # upper limit of elapsed time
#SBATCH -p regular            # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --get-user-env  # retrieve the login environment variables


# job execution via SLURM
# set stdout and stderr
STDOUT=log/$SLURM_JOB_NAME.o${SLURM_JOB_ID}
STDERR=log/$SLURM_JOB_NAME.e${SLURM_JOB_ID}

# set number of MPI processes per node
if [ -z "$SLURM_NTASKS" ]; then
    Ncore=`grep processor /proc/cpuinfo | wc -l` # flat MPI
    SLURM_NTASKS=`expr $Ncore / 2` # use 50% of the system
fi
PROCS_PER_NODE=`expr $SLURM_NTASKS / $SLURM_JOB_NUM_NODES`

# set number of MPI processes per socket
SOCKETS_PER_NODE=`lscpu | grep NUMA | grep CPU | wc -l`
PROCS_PER_SOCKET=`expr $PROCS_PER_NODE / $SOCKETS_PER_NODE`
if [ $PROCS_PER_SOCKET -lt 1 ]; then
    PROCS_PER_SOCKET=1
fi

# set tmp for matplotlib
MPL_CFG_DIR=/tmp/matplotlib

# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_NUM_NODES nodes"
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs per node"
TIME=`date`
echo "start: $TIME"

# module purge
# module load intelmpi
# module load texlive

# job execution by using SLURM
echo $LOADEDMODULES | grep -e openmpi -e intelmpi -e mvapich -e mvapich2
if [ $? -eq 0 ]; then
    # mpiexec -n $SLURM_NTASKS sh/mplwrapper.sh $PROCS_PER_NODE $PROCS_PER_SOCKET $MPL_CFG_DIR py/plot.py
    # mpiexec -n $SLURM_NTASKS sh/mplwrapper.sh $PROCS_PER_NODE $PROCS_PER_SOCKET $MPL_CFG_DIR py/draw.py
    # srun sh/mplwrapper.sh $PROCS_PER_NODE $PROCS_PER_SOCKET $MPL_CFG_DIR py/plot.py
    # srun sh/mplwrapper.sh $PROCS_PER_NODE $PROCS_PER_SOCKET $MPL_CFG_DIR py/draw.py
    mpiexec -n $SLURM_NTASKS --oversubscribe sh/mplwrapper.sh $PROCS_PER_NODE $PROCS_PER_SOCKET $MPL_CFG_DIR py/plot.py
    mpiexec -n $SLURM_NTASKS --oversubscribe sh/mplwrapper.sh $PROCS_PER_NODE $PROCS_PER_SOCKET $MPL_CFG_DIR py/draw.py
else
    echo "/usr/bin/time -f 'wall clock time: %e s\nmemory usage: %M KiB' python py/plot.py 1>>$STDOUT 2>>$STDERR"
    /usr/bin/time -f "wall clock time: %e s\nmemory usage: %M KiB" python py/plot.py 1>>$STDOUT 2>>$STDERR
    echo "/usr/bin/time -f 'wall clock time: %e s\nmemory usage: %M KiB' python py/draw.py 1>>$STDOUT 2>>$STDERR"
    /usr/bin/time -f "wall clock time: %e s\nmemory usage: %M KiB" python py/draw.py 1>>$STDOUT 2>>$STDERR
fi

# finish logging
TIME=`date`
echo "finish: $TIME"

exit 0

