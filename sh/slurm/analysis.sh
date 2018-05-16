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
    # PROBLEM=2
    # PROBLEM=14
    # PROBLEM=22
    # PROBLEM=26
    # PROBLEM=112
    # PROBLEM=130
    # PROBLEM=131
    # PROBLEM=132
    PROBLEM=133
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
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
    FILE=hernquist
    FINISH=23.5
    INTERVAL=0.5
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere
if [ $PROBLEM -eq 14 ]; then
    FILE=hd
    FINISH=104.0
    INTERVAL=8.0
    XMAX=5.0
    YMAX=5.0
    ZMAX=5.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere
if [ $PROBLEM -eq 22 ]; then
    FILE=disk
    FINISH=104.0
    INTERVAL=8.0
    XMAX=5.0
    YMAX=5.0
    ZMAX=5.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# dynamical stability of multi components galaxy model (only spherical components)
if [ $PROBLEM -eq 26 ]; then
    FILE=etg
    FINISH=1175.0
    INTERVAL=25.0
    XMAX=5.0
    YMAX=5.0
    ZMAX=5.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# reproduction of Kirihara et al. (2017) in the disk coordinate system
if [ $PROBLEM -eq 112 ]; then
    FILE=k17disk
    FINISH=1246.875
    INTERVAL=3.125
    XMAX=128.0
    YMAX=128.0
    ZMAX=128.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 130 ]; then
    FILE=halocusp_run
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 131 ]; then
    FILE=halocore1_run
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 132 ]; then
    FILE=halocore2_run
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 133 ]; then
    FILE=halocore3_run
    FINISH=14000.0
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
