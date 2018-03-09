#!/bin/bash
###############################################################
#SBATCH -J analysis           # name of job
#SBATCH -t 02:00:00           # upper limit of elapsed time
#SBATCH -p normal             # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=16           # number of total MPI processes, set to SLURM_NTASKS
##SBATCH --ntasks-per-socket=8 # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET
#SBATCH --ntasks-per-socket=16 # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET
#SBATCH --get-user-env        # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/extract
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=124
fi
###############################################################
# set number of N-body particles per bin to estimate density
if [ -z "$NCRIT" ]; then
    NCRIT=2048
fi
###############################################################
# set number of grid points for density maps
if [ -z "$NX" ]; then
    NX=256
    # NX=128
fi
if [ -z "$NY" ]; then
    NY=256
    # NY=128
fi
if [ -z "$NZ" ]; then
    NZ=256
    # NZ=128
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 124 ]; then
    FILE=halocusp_run
    # FINISH=14000.0
    FINISH=4400.0
    INTERVAL=100.0
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# count up number of snapshot files
if [ -z "$START" ]; then
    START=0
fi
if [ -z "$END" ]; then
    END=`echo "scale=0; $FINISH / $INTERVAL" | bc`
fi
if [ -z "$INCREMENT" ]; then
    INCREMENT=1
fi
###############################################################
# set input arguments
OPTION="-file=$FILE -start=$START -end=$END -interval=$INCREMENT -ncrit=$NCRIT -nx=$NX -xmin=$XMIN -xmax=$XMAX -ny=$NY -ymin=$YMIN -ymax=$YMAX -nz=$NZ -zmin=$ZMIN -zmax=$ZMAX"
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
