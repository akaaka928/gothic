#!/bin/bash
###############################################################
#SBATCH -J m31obs              # name of job
#SBATCH -t 02:00:00            # upper limit of elapsed time
#SBATCH -p normal              # partition name
#SBATCH --nodes=1              # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=16            # number of total MPI processes, set to SLURM_NTASKS
#SBATCH --ntasks-per-socket=16 # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET
# #SBATCH --ntasks=1            # number of total MPI processes, set to SLURM_NTASKS
# #SBATCH --ntasks-per-socket=1 # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET
#SBATCH --get-user-env         # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/m31obs
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=112
fi
###############################################################
# set number of N-body particles per bin to estimate density
if [ -z "$NCRIT" ]; then
    NCRIT=2048
fi
###############################################################
# set number of grid points for density maps
if [ -z "$NX" ]; then
    # NX=512
    NX=256
fi
if [ -z "$NY" ]; then
    # NY=512
    NY=256
fi
if [ -z "$NZ" ]; then
    # NZ=512
    NZ=256
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# reproduction of Kirihara et al. (2017) in the disk coordinate system
if [ $PROBLEM -eq 112 ]; then
    FILE=k17disk
    FINISH=1246.875
    INTERVAL=3.125
    XIMIN=-5.0
    XIMAX=7.0
    ETAMIN=-8.5
    ETAMAX=3.5
    DMIN=-80.0
    DMAX=130.0
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
OPTION="-file=$FILE -start=$START -end=$END -interval=$INCREMENT -ncrit=$NCRIT -nxi=$NX -ximin=$XIMIN -ximax=$XIMAX -neta=$NY -etamin=$ETAMIN -etamax=$ETAMAX -nD=$NZ -Dmin=$DMIN -Dmax=$DMAX"
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
