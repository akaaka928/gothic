#!/bin/bash
#PBS -A GALAXY
#PBS -q gpu
#PBS -N gothic
#PBS -b 1
#PBS -l elapstim_req=24:00:00
#PBS -T mvapich
#PBS -v NQSV_MPI_VER=2.3.1/intel-cuda10.1
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/gothic
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=0
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
    # ABSERR=6.103515625e-5
fi
###############################################################
if [ -z "$REBUILD" ]; then
    REBUILD=16
fi
###############################################################
if [ -z "$BRENT" ]; then
    BRENT=1.0
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# cold collapse of a uniform sphere
if [ $PROBLEM -eq 0 ]; then
    FILE=nws-live
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$PBS_JOBID"
# -file=$FILE
###############################################################


###############################################################
# job execution via NQSV
###############################################################
# set stdout and stderr
STDOUT=log/${FILE}_$PBS_JOBNAME.o${PBS_JOBID}
STDERR=log/${FILE}_$PBS_JOBNAME.e${PBS_JOBID}
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
# execute the job
PROCS=4
PROCS_PER_NODE=4
# PROCS=2
# PROCS_PER_NODE=2
mpiexec ${NQSII_MPIOPTS} -np ${PROCS} -genv MV2_NUM_HCAS 4 sh/cygnus/split.sh $EXEC ${FILE}_$PBS_JOBNAME $PBS_JOBID $PROCS_PER_NODE $PROCS_PER_NODE $OPTION
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
