#!/bin/bash
###############################################################
#SBATCH -J NWSrun-cmaes # name of job
#SBATCH -t 02:00:00     # upper limit of elapsed time
#SBATCH -p normal       # partition name
#SBATCH --nodes=1       # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --get-user-env  # retrieve the login environment variables
###############################################################


###############################################################
# number of runs in this generation
if [ -z "$NRUN" ]; then
    NRUN=64
fi
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/cmaes
###############################################################


###############################################################
# problem specific configurations
###############################################################
# number of parameters
if [ -z "$NDIM" ]; then
    NDIM=9
fi
###############################################################
initial prediction
###############################################################
sigma??
###############################################################
tolerance values
###############################################################
# set input arguments
OPTION="-file=$FILE -ndim=$NDIM -lambda=$NRUN"
###############################################################


###############################################################
# job execution via SLURM
###############################################################
# set stdout and stderr
STDOUT=log/${FILE}_$SLURM_JOB_NAME.o${SLURM_JOB_ID}
STDERR=log/${FILE}_$SLURM_JOB_NAME.e${SLURM_JOB_ID}
###############################################################
# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs"
TIME=`date`
echo "start: $TIME"
###############################################################
