#!/bin/bash
###############################################################
#SBATCH -J NWScmaes     # name of job
#SBATCH -t 02:00:00     # upper limit of elapsed time
#SBATCH -p normal       # partition name
#SBATCH --nodes=1       # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --get-user-env  # retrieve the login environment variables
###############################################################


###############################################################
# generation of the sequences
if [ -z "$GEN" ]; then
    GEN=1
fi
###############################################################
# name of the series
if [ -z "$FILE" ]; then
    FILE=cusp
fi
###############################################################
# initial prediction of the mean value
if [ -z "$INIT" ]; then
    INIT=cfg/$FILE/mean.ini
    SIGMA=2.0
fi
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
# tolerance values
###############################################################
# set input arguments
OPTION="-file=$FILE -ndim=$NDIM -lambda=$NRUN"
###############################################################
if [ $GEN -gt 1 ]; then
    PREV=`expr $GEN - 1`
fi
###############################################################
# specify the previous generation
if [ -n "$PREV" ]; then
    OPTION="$OPTION -restart=$PREV"
fi
###############################################################
# specify mean distribution of the first generation
if [ -n "$INIT" ]; then
    OPTION="$OPTION -init=$INIT"
fi
###############################################################
# specify the value of sigma
if [ -n "$SIGMA" ]; then
    OPTION="$OPTION -sigma=$SIGMA"
fi
###############################################################
# specify tolerance values
if [ -n "$FTOL" ]; then
    OPTION="$OPTION -ftol=$FTOL"
fi
if [ -n "$HTOL" ]; then
    OPTION="$OPTION -htol=$HTOL"
fi
if [ -n "$XTOL" ]; then
    OPTION="$OPTION -xtol=$XTOL"
fi
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
mkdir cfg/${FILE}/gen${GEN}
###############################################################
# execute the job
if [ `which numactl` ]; then
    echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
