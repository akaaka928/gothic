#!/bin/bash
#PBS -A GALAXY
#PBS -q gpu
#PBS -N mass90-prada-snap160_vel140
#PBS -b 1
#PBS -l elapstim_req=24:00:00
#PBS -T openmpi
#PBS -v NQSV_MPI_VER=gdr/4.0.3/intel19.0.5-cuda10.2
###############################################################
# MASS=mass80
# MASS=mass85
# MASS=mass90
MASS=mass95
###############################################################
# HALO=test
HALO=prada
# HALO=ishiyama
# HALO=gilman
###############################################################
EPOCH=snap160_vel140
# EPOCH=snap180_vel140
# EPOCH=snap200_vel140
# EPOCH=snap220_vel140
# EPOCH=snap240_vel140
# EPOCH=snap260_vel140
# EPOCH=snap280_vel140
# EPOCH=snap300_vel140
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
if [ -z "$REBUILD" ]; then
    REBUILD=16
fi
###############################################################
if [ -z "$BRENT" ]; then
    BRENT=1.0
fi
###############################################################
if [ -z "$SPHEPOT" ]; then
    SPHEPOT=m31
fi
if [ -z "$DISKPOT" ]; then
    DISKPOT=m31
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# FILE=${MASS}-${HALO}-${EPOCH}_${ORBIT}
FILE=${MASS}-${HALO}-${EPOCH}
###############################################################
# set input arguments
# OPTION="-absErr=$ABSERR -file=$FILE -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$PBS_JOBID"
OPTION="-absErr=$ABSERR -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$PBS_JOBID"
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
export MODULEPATH=/work/GALAXY/share/opt/modules:/work/GALAXY/ymiki/opt/module:$MODULEPATH
module load openmpi/$NQSV_MPI_VER
module load phdf5
module load cub
###############################################################
# execute the job
PROCS=4
PROCS_PER_NODE=4
# PROCS=2
# PROCS_PER_NODE=2
mpiexec ${NQSII_MPIOPTS} -np ${PROCS} -x UCX_MAX_RNDV_RAILS=4 sh/cygnus/nws/split.sh $EXEC ${FILE} ${FILE}_$PBS_JOBNAME $PBS_JOBID $PROCS_PER_NODE $PROCS_PER_NODE $OPTION
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
