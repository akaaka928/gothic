#!/bin/bash
#PBS -A GALAXY
#PBS -q gpu
#PBS -N gothic
#PBS -b 16
#PBS -l elapstim_req=24:00:00
#PBS -T mvapich
#PBS -v NQSV_MPI_VER=2.3.1/intel-cuda10.1
#PBS -M ymiki@cc.u-tokyo.ac.jp
#PBS -m be
###############################################################
NGPUS_PER_NODE=4
NGPUS_PER_SOCKET=4
###############################################################
NQUEUES=16
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
# set minimum dt
if [ -z "$DTMIN" ]; then
    DTMIN=1.220703125e-2
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -jobID=$PBS_JOBID -dtmin=$DTMIN"
###############################################################


###############################################################
# job execution via PBS
###############################################################
# set stdout and stderr
STDOUT=log/${SERIES}_$PBS_JOBNAME.o${PBS_JOBID}
STDERR=log/${SERIES}_$PBS_JOBNAME.e${PBS_JOBID}
###############################################################
# start logging
cd $PBS_O_WORKDIR
TIME=`date`
echo "start: $TIME"
###############################################################
module purge
export MODULEPATH=$MODULEPATH:/work/CSPP/ymiki/opt/module
module load mvapich/2.3.1/intel-cuda10.1
module load phdf5
module load cub
###############################################################


###############################################################
FULL=${SERIES}-gen${GEN}.lst
LIST=${SERIES}-gen${GEN}-split.lst
# split -d -n l/$PROCS $FULL $LIST
###############################################################


###############################################################
# execute the job
mpiexec ${NQSII_MPIOPTS} -np ${PROCS} -ppn ${NGPUS_PER_NODE} -genv MV2_NUM_HCAS 4 -genv MV2_ENABLE_AFFINITY 0 sh/cygnus/NWSsplit.sh $EXEC log/${SERIES}_${REQUEST} $PBS_JOBID $NGPUS_PER_NODE $NGPUS_PER_SOCKET $LIST $OPTION
###############################################################


###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################
exit 0
###############################################################
