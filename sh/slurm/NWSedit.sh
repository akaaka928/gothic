#!/bin/bash
###############################################################
#SBATCH -J NWSedit      # name of job
#SBATCH -t 02:00:00     # upper limit of elapsed time
#SBATCH -p normal       # partition name
#SBATCH --nodes=1       # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --get-user-env  # retrieve the login environment variables
###############################################################


###############################################################
# generation of the simulation
if [ -z "$GENERATION" ]; then
    GENERATION=1
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
fi
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/editor
###############################################################
# dump file generation interval (in units of minute)
if [ -z "$SAVE" ]; then
    SAVE=55.0
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# gravitational softening length
if [ -z "$EPS" ]; then
    EPS=1.5625e-2
fi
###############################################################
# safety parameter for time steps
if [ -z "$ETA" ]; then
    ETA=0.5
fi
###############################################################
# interval of snapshot files (in units of astrophysical unit)
if [ -z "$INTERVAL" ]; then
    INTERVAL=25.0
fi
###############################################################
# final time of the simulation (in units of astrophysical unit)
if [ -z "$FINISH" ]; then
    FINISH=14000.0
fi
###############################################################


###############################################################
# job execution via SLURM
###############################################################
# set stdout and stderr
STDOUT=log/${SERIES}_$SLURM_JOB_NAME.o${SLURM_JOB_ID}
STDERR=log/${SERIES}_$SLURM_JOB_NAME.e${SLURM_JOB_ID}
###############################################################
# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs"
TIME=`date`
echo "start: $TIME"
###############################################################


# edit $NRUN models in generation $GENERATION
###############################################################
for ii in `seq 1 $NRUN`
do
    # set model ID
    ID=`expr $ii - 1`
    FILE=${SERIES}-gen${GENERATION}-run${ID}
    CFG=$SERIES/gen${GENERATION}/run${ID}-edit.cfg

    # set input arguments
    OPTION="-file=$FILE -list=$CFG -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE"

    # execute the job
    if [ `which numactl` ]; then
	echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
	numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
    else
	echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
	$EXEC $OPTION 1>>$STDOUT 2>>$STDERR
    fi
done
###############################################################


###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
