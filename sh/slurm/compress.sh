#!/bin/bash
###############################################################
#SBATCH -J compress            # name of job
#SBATCH -t 24:00:00            # upper limit of elapsed time
#SBATCH -p normal              # partition name
#SBATCH --nodes=1              # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=16            # number of total MPI processes, set to SLURM_NTASKS
#SBATCH --ntasks-per-socket=8  # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET
#SBATCH --get-user-env         # retrieve the login environment variables
###############################################################


SRCDIR=dat
DSTDIR=cdat
if [ ! -e $DSTDIR ]; then
    mkdir $DSTDIR
fi

LOG=compress.log
LIST=compress_list.txt
if [ -e $LIST ]; then
    rm -f $LIST
fi
# ls -1 $SRCDIR/*.snp???.h5 > $LIST
ls -1 $SRCDIR/*.h5 > $LIST


SUBLIST=compress_sub
split -d -n l/$SLURM_NTASKS $LIST $SUBLIST
rm -f $LIST

mpiexec -n $SLURM_NTASKS sh/slurm/compress_sub.sh $SUBLIST $SRCDIR $DSTDIR


exit 0
