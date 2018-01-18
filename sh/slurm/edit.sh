#!/bin/bash
#SBATCH -J editor             # name of job
#SBATCH -t 00:30:00           # upper limit of elapsed time
#SBATCH -p normal             # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --get-user-env        # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/editor
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=0
fi
###############################################################
# dump file generation interval (in units of minute)
if [ -z "$SAVE" ]; then
    # SAVE=2.0
    # SAVE=60.0
    SAVE=140.0
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# dynamical stability of a disk in a spherical potential field
if [ $PROBLEM -eq 0 ]; then
    FILE=disk
    CFG=external/split.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=11.75
    INTERVAL=0.25
fi
###############################################################
# set input arguments
OPTION="-file=$FILE -list=$CFG -config=$CONFIG -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE"
###############################################################


###############################################################
# job execution via SLURM
###############################################################
# set stdout and stderr
STDOUT=log/$SLURM_JOB_NAME.$SLURM_JOB_ID.out
STDERR=log/$SLURM_JOB_NAME.$SLURM_JOB_ID.err
###############################################################
# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs"
TIME=`date`
echo "start: $TIME"
###############################################################
# execute the job
if [ `which numactl` ]; then
    # run with numactl
    echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    # run without numactl
    echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
