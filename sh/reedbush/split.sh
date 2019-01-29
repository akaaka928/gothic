#!/bin/bash
###############################################################
#PBS -q l-regular
#PBS -l select=1:mpiprocs=4:ompthreads=1
#PBS -W group_list=jh180045l
#PBS -l walltime=24:00:00
#PBS -N gothic
###############################################################
PROCS=4
PROCS_PER_NODE=4
PROCS_PER_SOCKET=2
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/gothic
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=200
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


###############################################################
# problem specific configurations
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 200 ]; then
    FILE=w3
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 210 ]; then
    FILE=w5
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 220 ]; then
    FILE=w7
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 230 ]; then
    FILE=m12
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 240 ]; then
    FILE=m09
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 250 ]; then
    FILE=a00
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 260 ]; then
    FILE=a05
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 270 ]; then
    FILE=a10
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 280 ]; then
    FILE=a15
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 290 ]; then
    FILE=a20
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -Nx=$NX -Ny=$NY -Nz=$NZ -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$PBS_JOBID"
###############################################################


###############################################################
# job execution via PBS
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
export MODULEPATH=$MODULEPATH:/lustre/jh180045l/share/opt/Modules
. /etc/profile.d/modules.sh
module purge
module load pbsutils
module load intel
# module load cuda9/9.0.176 mvapich2/gdr/2.3a/gnu phdf5/mv2
module load cuda9/9.1.85 openmpi/gdr/2.1.2/intel phdf5/ompi
# module load cuda9 intel-mpi phdf5/impi
module load cub
###############################################################
# execute the job
mpirun -np $PROCS ./wrapper.sh $EXEC $FILE $PBS_JOBID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
