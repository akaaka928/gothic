#!/bin/sh
###############################################################
#PBS -q u-short
#PBS -l select=1:ncpus=36:mpiprocs=1:ompthreads=36
#PBS -W group_list=jh180045u
#PBS -l walltime=00:10:00
#PBS -N editor
###############################################################
# l-small for jh180045l
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/editor
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=0
    # PROBLEM=1
    PROBLEM=2
    # PROBLEM=10
    # PROBLEM=11
    # PROBLEM=12
    # PROBLEM=21
    # PROBLEM=22
    # PROBLEM=23
    # PROBLEM=24
fi
###############################################################
# dump file generation interval (in units of minute)
if [ -z "$SAVE" ]; then
    SAVE=55.0
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
# dynamical stability of two disks in a spherical potential field
if [ $PROBLEM -eq 1 ]; then
    FILE=ltg_disk
    CFG=external/ltg_split.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1175.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of an exponential disk in a spherical potential field (NFW sphere)
if [ $PROBLEM -eq 2 ]; then
    FILE=disk
    CFG=external/hd_split.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=104.0
    INTERVAL=8.0
fi
###############################################################
# dynamical stability of M31 in the observed coordinate
if [ $PROBLEM -eq 10 ]; then
    FILE=m31_obs
    CFG=gss/obs.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=3175.0
    INTERVAL=25.0
fi
###############################################################
# reproduction of Miki et al. (2016) in the disk coordinate system
if [ $PROBLEM -eq 11 ]; then
    FILE=m16king
    CFG=gss/satellite.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=12.5
    # FINISH=100.0
    # FINISH=1050.0
    # INTERVAL=1.5625
    FINISH=1246.875
    INTERVAL=3.125
fi
###############################################################
# reproduction of Kirihara et al. (2017) in the disk coordinate system
if [ $PROBLEM -eq 12 ]; then
    FILE=k17disk
    CFG=gss/satellite.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=12.5
    # FINISH=100.0
    # FINISH=1050.0
    # INTERVAL=1.5625
    FINISH=1246.875
    INTERVAL=3.125
fi
###############################################################
# GSS simulation in observed frame with live M31
if [ $PROBLEM -eq 20 ]; then
    FILE=gss
    CFG=gss/obs_run.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1246.875
    INTERVAL=3.125
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 21 ]; then
    FILE=halocore1_run
    CFG=pbh/halocore1_run.cfg
    EPS=9.765625000e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 22 ]; then
    FILE=halocore2_run
    CFG=pbh/halocore2_run.cfg
    EPS=9.765625000e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 23 ]; then
    FILE=halocore3_run
    CFG=pbh/halocore3_run.cfg
    EPS=9.765625000e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 24 ]; then
    FILE=halocusp_run
    CFG=pbh/halocusp_run.cfg
    EPS=9.765625000e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# set input arguments
OPTION="-file=$FILE -list=$CFG -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE"
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
# module load mvapich2/gdr/2.3a/gnu phdf5/mv2
module load openmpi/gdr/2.1.2/intel phdf5/ompi
# module load intel-mpi phdf5/impi
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
exit 0
###############################################################
