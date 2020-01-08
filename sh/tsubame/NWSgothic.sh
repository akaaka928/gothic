#!/bin/sh
#$ -cwd
#$ -l f_node=16
#$ -l h_rt=24:00:00
#$ -N NWSgothic
###############################################################
NGPUS_PER_NODE=4
NGPUS_PER_SOCKET=2
###############################################################
PROCS=`expr $NQUEUES \* $NGPUS_PER_NODE`
###############################################################


###############################################################
# generation of the simulation
if [ -z "$GEN" ]; then
    GEN=1
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
if [ -z "$EXEC" ]; then
    EXEC=bin/gothic
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
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ -z "$SPHEPOT" ]; then
    SPHEPOT=m31
fi
if [ -z "$DISKPOT" ]; then
    DISKPOT=m31
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -jobID=$JOB_ID"
###############################################################


###############################################################
# job execution via UNIVA Grid Engine
###############################################################
# set stdout and stderr
STDOUT=log/${SERIES}_$REQUEST.o$JOB_ID
STDERR=log/${SERIES}_$REQUEST.e$JOB_ID
###############################################################
# load modules
. /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:/gs/hs1/jh180045/share/opt/Modules
module load intel/19.0.0.117 cuda/9.2.148 openmpi/2.1.2-opa10.9
module load cub phdf5
module list 1>>$STDOUT 2>>$STDERR
export OMP_NUM_THREADS=7
###############################################################
cat $PE_HOSTFILE 1>>$STDOUT 2>>$STDERR
TIME=`date`
echo "start: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################


###############################################################
FULL=${SERIES}-gen${GEN}.lst
LIST=${SERIES}-gen${GEN}-split.lst
split -d -n l/$PROCS $FULL $LIST
# mpiexec -n $SLURM_NTASKS sh/slurm/compress_sub.sh $SUBLIST $SRCDIR $DSTDIR
###############################################################


###############################################################
# execute the job
# if [ `which numactl` ]; then
#     # run with numactl
#     echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
#     numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
# else
#     # run without numactl
#     echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
#     $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
# fi
mpirun -np $PROCS -npernode $NGPUS_PER_NODE -x LD_LIBRARY_PATH -x PATH -x OMP_NUM_THREADS sh/split.sh $EXEC log/${SERIES}_${REQUEST} $JOB_ID $NGPUS_PER_NODE $NGPUS_PER_SOCKET $LIST $OPTION
# mpirun -tag-output -mca pml cm -mca mtl psm2 -np $PROCS -npernode $NGPUS_PER_NODE --bind-to core -cpu-set 0,7,14,21 -x CUDA_VISIBLE_DEVICES=0,1,2,3 -x PSM2_CUDA=1 -x PSM2_GPUDIRECT=1 -x LD_LIBRARY_PATH -x PATH $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
###############################################################



###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################
exit 0
###############################################################
