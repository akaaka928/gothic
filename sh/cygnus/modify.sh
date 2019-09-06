#!/bin/bash
#PBS -A GALAXY
#PBS -q gpu
#PBS -N editor
#PBS -b 1
#PBS -l elapstim_req=00:5:00
#PBS -v OMP_NUM_THREADS=24
#PBS -T mvapich
#PBS -v NQSV_MPI_VER=2.3.1/intel-cuda10.1
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/editor
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=68
    # PROBLEM=69
    # PROBLEM=70
    # PROBLEM=71

    # PROBLEM=108
    # PROBLEM=109
    # PROBLEM=110
    # PROBLEM=111

    # collision with live halo (orbit 4 and 5)
    PROBLEM=148
    # PROBLEM=149
    # PROBLEM=150
    # PROBLEM=151

    # collision with live halo (orbit 6 and 7)
    # PROBLEM=188
    # PROBLEM=189
    # PROBLEM=190
    # PROBLEM=191
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
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 13 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
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
# preparation for throwing DM subhalo into NW stream
if [ $PROBLEM -eq 30 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=7400.0
    INTERVAL=25.0
fi
###############################################################
# continue N-body simulation to reproduce NW stream
if [ $PROBLEM -eq 31 ]; then
    FILE=nws-continue
    CFG=nws/continue.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 40 ]; then
    FILE=nws-test-m7_0-orbit0
    CFG=nws/test/m7_0-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 41 ]; then
    FILE=nws-test-m7_0-orbit1
    CFG=nws/test/m7_0-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 42 ]; then
    FILE=nws-test-m7_5-orbit0
    CFG=nws/test/m7_5-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 43 ]; then
    FILE=nws-test-m7_5-orbit1
    CFG=nws/test/m7_5-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 44 ]; then
    FILE=nws-test-m8_0-orbit0
    CFG=nws/test/m8_0-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 45 ]; then
    FILE=nws-test-m8_0-orbit1
    CFG=nws/test/m8_0-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 46 ]; then
    FILE=nws-test-m8_5-orbit0
    CFG=nws/test/m8_5-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 47 ]; then
    FILE=nws-test-m8_5-orbit1
    CFG=nws/test/m8_5-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 48 ]; then
    FILE=nws-test-m9_0-orbit0
    CFG=nws/test/m9_0-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 49 ]; then
    FILE=nws-test-m9_0-orbit1
    CFG=nws/test/m9_0-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 50 ]; then
    FILE=nws-test-m9_5-orbit0
    CFG=nws/test/m9_5-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 51 ]; then
    FILE=nws-test-m9_5-orbit1
    CFG=nws/test/m9_5-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 60 ]; then
    FILE=nws-live-m7_0-orbit0
    CFG=nws/live/m7_0-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 61 ]; then
    FILE=nws-live-m7_0-orbit1
    CFG=nws/live/m7_0-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 62 ]; then
    FILE=nws-live-m7_5-orbit0
    CFG=nws/live/m7_5-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 63 ]; then
    FILE=nws-live-m7_5-orbit1
    CFG=nws/live/m7_5-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 64 ]; then
    FILE=nws-live-m8_0-orbit0
    CFG=nws/live/m8_0-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 65 ]; then
    FILE=nws-live-m8_0-orbit1
    CFG=nws/live/m8_0-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 66 ]; then
    FILE=nws-live-m8_5-orbit0
    CFG=nws/live/m8_5-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 67 ]; then
    FILE=nws-live-m8_5-orbit1
    CFG=nws/live/m8_5-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 68 ]; then
    FILE=nws-live-m9_0-orbit0
    CFG=nws/live/m9_0-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 69 ]; then
    FILE=nws-live-m9_0-orbit1
    CFG=nws/live/m9_0-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 70 ]; then
    FILE=nws-live-m9_5-orbit0
    CFG=nws/live/m9_5-orbit0.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 71 ]; then
    FILE=nws-live-m9_5-orbit1
    CFG=nws/live/m9_5-orbit1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=6600.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 108 ]; then
    FILE=nws-live-m9_0-orbit2
    CFG=nws/live/m9_0-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 109 ]; then
    FILE=nws-live-m9_0-orbit3
    CFG=nws/live/m9_0-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 110 ]; then
    FILE=nws-live-m9_5-orbit2
    CFG=nws/live/m9_5-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 111 ]; then
    FILE=nws-live-m9_5-orbit3
    CFG=nws/live/m9_5-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 148 ]; then
    FILE=nws-live-m9_0-orbit4
    CFG=nws/live/m9_0-orbit4.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 149 ]; then
    FILE=nws-live-m9_0-orbit5
    CFG=nws/live/m9_0-orbit5.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 150 ]; then
    FILE=nws-live-m9_5-orbit4
    CFG=nws/live/m9_5-orbit4.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 151 ]; then
    FILE=nws-live-m9_5-orbit5
    CFG=nws/live/m9_5-orbit5.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 188 ]; then
    FILE=nws-live-m9_0-orbit6
    CFG=nws/live/m9_0-orbit6.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 189 ]; then
    FILE=nws-live-m9_0-orbit7
    CFG=nws/live/m9_0-orbit7.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 190 ]; then
    FILE=nws-live-m9_5-orbit6
    CFG=nws/live/m9_5-orbit6.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 191 ]; then
    FILE=nws-live-m9_5-orbit7
    CFG=nws/live/m9_5-orbit7.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# set input arguments
OPTION="-file=$FILE -list=$CFG -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE"
###############################################################
# set environmental variables for OpenMP
# OMP_OPT_ENV="env OMP_DISPLAY_ENV=verbose OMP_PLACES=cores"
OMP_OPT_ENV="env KMP_AFFINITY=verbose,granularity=core,scatter"
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
