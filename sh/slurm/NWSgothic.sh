#!/bin/bash
###############################################################
#SBATCH -J NWSgothic          # name of job
#SBATCH -t 24:00:00           # upper limit of elapsed time
#SBATCH -p normal             # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=2            # number of total MPI processes, set to SLURM_NTASKS (must be equal to number of GPUs)
#SBATCH --ntasks-per-socket=2 # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET (must be equal to number of GPUs per socket)
#SBATCH --get-user-env        # retrieve the login environment variables
###############################################################
PROCS=$SLURM_NTASKS
###############################################################


###############################################################
# generation of the simulation
if [ -z "$GEN" ]; then
    GEN=1
    # GEN=2
fi
###############################################################
# name of the series
if [ -z "$SERIES" ]; then
    SERIES=cusp
fi
###############################################################
# number of runs in this generation
if [ -z "$NRUN" ]; then
    NRUN=64
    # NRUN=4
fi
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/gothic
###############################################################
# topology of MPI processes
if [ -z "$NX" ]; then
    NX=2
fi
if [ -z "$NY" ]; then
    NY=1
fi
if [ -z "$NZ" ]; then
    NZ=1
fi
if [ $SLURM_NTASKS -eq 1 ]; then
    NX=1
    NY=1
    NZ=1
fi
PROCS=`expr $NX \* $NY \* $NZ`
if [ $PROCS -ne $SLURM_NTASKS ]; then
    echo "product of $NX, $NY, and $NZ must be equal to the number of total MPI processes ($SLURM_NTASKS)"
    exit 1
fi
###############################################################
# value of accuracy controling parameter: GADGET MAC by Springel (2005)
if [ -z "$ABSERR" ]; then
    # ABSERR=1.250000000e-1
    # ABSERR=6.250000000e-2
    # ABSERR=3.125000000e-2
    # ABSERR=1.562500000e-2
    # ABSERR=7.812500000e-3
    # ABSERR=3.906250000e-3
    ABSERR=1.953125000e-3
    # ABSERR=9.765625000e-4
    # ABSERR=4.882812500e-4
    # ABSERR=2.441406250e-4
    # ABSERR=1.220703125e-4
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
if [ -z "$SPHEPOT" ]; then
    SPHEPOT=m31
fi
if [ -z "$DISKPOT" ]; then
    DISKPOT=m31
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -jobID=$SLURM_JOB_ID"
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
FULL=${SERIES}-gen${GEN}.lst
LIST=${SERIES}-gen${GEN}-split.lst
if [ $PROCS -lt 10 ]; then
    split -a 1 -d -n l/$PROCS $FULL $LIST
else
    split -d -n l/$PROCS $FULL $LIST
fi
###############################################################
# execute the job
mpiexec -n $SLURM_NTASKS sh/split.sh $EXEC log/${SERIES}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $LIST $OPTION
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
