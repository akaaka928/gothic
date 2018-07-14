#!/bin/bash
###############################################################
#SBATCH -J m31plt             # name of job
#SBATCH -t 02:00:00           # upper limit of elapsed time
#SBATCH -p normal             # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --get-user-env        # retrieve the login environment variables
###############################################################


###############################################################
# job execution via SLURM
###############################################################
# set stdout and stderr
STDOUT=log/$SLURM_JOB_NAME.o${SLURM_JOB_ID}
STDERR=log/$SLURM_JOB_NAME.e${SLURM_JOB_ID}
###############################################################
# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_NUM_NODES nodes"
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs per node"
TIME=`date`
echo "start: $TIME"
###############################################################
# execute the job
echo "/usr/bin/time -f 'wall clock time: %e s\nmemory usage: %M KiB' python py/shownws.py 1>>$STDOUT 2>>$STDERR"
/usr/bin/time -f "wall clock time: %e s\nmemory usage: %M KiB" python py/shownws.py 1>>$STDOUT 2>>$STDERR
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
