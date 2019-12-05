#!/bin/bash
###############################################################
#SBATCH -J subaru              # name of job
#SBATCH -t 24:00:00            # upper limit of elapsed time
#SBATCH -p normal              # partition name
#SBATCH --nodes=1              # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=16            # number of total MPI processes, set to SLURM_NTASKS
#SBATCH --ntasks-per-socket=8  # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET
#SBATCH --get-user-env         # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/subaru
###############################################################
# set number of grid points for density maps
if [ -z "$NX" ]; then
    NX=256
    # NX=128
    # NX=32
fi
if [ -z "$NY" ]; then
    NY=256
    # NY=128
    # NY=32
fi
if [ -z "$NV" ]; then
    NV=256
    # NV=128
    # NV=32
fi
###############################################################
if [ -z "$FILE" ]; then
    FILE=nws-continue
fi
###############################################################
if [ -z "$START" ]; then
    # START=0
    # START=50
    START=32
fi
###############################################################
if [ -z "$END" ]; then
    # END=560
    # END=66
    # END=100
    END=112
fi
###############################################################
if [ -z "$INTERVAL" ]; then
    INTERVAL=1
fi
###############################################################
# set input arguments
OPTION="-file=$FILE -start=$START -end=$END -interval=$INTERVAL -nx=$NX -ny=$NY -nv=$NV"
###############################################################


###############################################################
# job execution via SLURM
###############################################################
# set number of MPI processes per node
PROCS_PER_NODE=`expr $SLURM_NTASKS / $SLURM_JOB_NUM_NODES`
###############################################################
# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_NUM_NODES nodes"
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs per node"
TIME=`date`
echo "start: $TIME"
###############################################################
# execute the job
echo "mpiexec -n $SLURM_NTASKS sh/wrapper.sh $EXEC log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTION"
mpiexec -n $SLURM_NTASKS sh/wrapper.sh $EXEC log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTION
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
