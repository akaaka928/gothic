#!/bin/bash
#SBATCH -J observation        # name of job
#SBATCH -t 24:00:00           # upper limit of elapsed time
#SBATCH -p normal             # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --get-user-env        # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/observation
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=0
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# cold collapse of a uniform sphere
if [ $PROBLEM -eq 0 ]; then
    FILE=nws
    FINISH=14000.0
    INI=0
    END=560
    INTERVAL=1

    LOGM_DM=7.69897
    LOGRS_DM=0.0
    LOGM_STAR=6.0
    LOGR0_STAR=-0.522878745
    LOGRT_STAR=0.149162743

    THETA=0.0
    VRAD=100.0
    VTAN=100.0
    VANGLE=0.0
    EPS=1.562500e-2
fi
###############################################################
# set input arguments
OPTION="-file=$FILE -start=$INI -end=$END -interval=$INTERVAL -logM_dm=$LOGM_DM -logrs_dm=$LOGRS_DM -logM_star=$LOGM_STAR -logr0_star=$LOGR0_STAR -logrt_star=$LOGRT_STAR -theta=$THETA -vrad=$VRAD -vtan=$VTAN -vangle=$VANGLE -eps=$EPS -ft=$FINISH"
###############################################################
# set environmental variables for OpenMP
OMP_OPT_ENV="env OMP_DISPLAY_ENV=verbose OMP_PLACES=cores"
# OMP_OPT_ENV="env KMP_AFFINITY=verbose,granularity=core,scatter"
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
# execute the job
if [ `which numactl` ]; then
    # run with numactl
    echo "$OMP_OPT_ENV numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    $OMP_OPT_ENV numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    # run without numactl
    echo "$OMP_OPT_ENV $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    $OMP_OPT_ENV $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
