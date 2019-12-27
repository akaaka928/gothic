#!/bin/bash
###############################################################
#SBATCH -J NWSinit-magi # name of job
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
# number of runs in this generation
if [ -z "$NRUN" ]; then
    NRUN=64
fi
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/magi
###############################################################
# set environmental variables for OpenMP
OMP_OPT_ENV="env OMP_DISPLAY_ENV=verbose OMP_PLACES=cores"
# OMP_OPT_ENV="env KMP_AFFINITY=verbose,granularity=core,scatter"
###############################################################
# number of N-body particles
if [ -z "$NTOT" ]; then
    NTOT=1048576
fi
###############################################################
# dump file generation interval (in units of minute)
if [ -z "$SAVE" ]; then
    SAVE=55.0
fi
###############################################################
# denoise distribution function
if [ -z "$DENOISE_DF" ]; then
    DENOISE_DF=1
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
STDOUT=log/${FILE}_$SLURM_JOB_NAME.o${SLURM_JOB_ID}
STDERR=log/${FILE}_$SLURM_JOB_NAME.e${SLURM_JOB_ID}
###############################################################
# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs"
TIME=`date`
echo "start: $TIME"
###############################################################


# generate $NRUN models in generation $GENERATION
###############################################################
for ii in `seq 1 $NRUN`
do
    # set model ID
    ID=`expr $ii - 1`
    FILE=gen${GENERATION}-run${ID}-magi
    CONFIG=nws-find/gen${GENERATION}/run${ID}-magi.cfg

    # set input arguments
    OPTION="-file=$FILE -config=$CONFIG -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE -denoisingDistributionFunction=$DENOISE_DF"

    # execute the job
    if [ `which numactl` ]; then
	echo "$OMP_OPT_ENV numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
	$OMP_OPT_ENV numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
    else
	echo "$OMP_OPT_ENV $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
	$OMP_OPT_ENV $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
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
