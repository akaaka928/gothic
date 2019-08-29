#!/bin/sh
#$ -cwd
#$ -l q_core=1
#$ -l h_rt=0:05:00
#$ -N editor
#$ -hold_jid magi
###############################################################
# NCORES_PER_NODE=28 # for f_node
# NCORES_PER_NODE=14 # for h_node
# NCORES_PER_NODE=7 # for q_node
NCORES_PER_NODE=4 # for q_core
###############################################################
export OMP_NUM_THREADS=$NCORES_PER_NODE
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/editor
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=30
    # PROBLEM=31
    PROBLEM=32
    # PROBLEM=33
    # PROBLEM=40
    # PROBLEM=41
    # PROBLEM=42
    # PROBLEM=43
    # PROBLEM=44
    # PROBLEM=45
    # PROBLEM=46
    # PROBLEM=47
    # PROBLEM=48
    # PROBLEM=49
    # PROBLEM=50
    # PROBLEM=51
    # PROBLEM=60
    # PROBLEM=61
    # PROBLEM=62
    # PROBLEM=63
    # PROBLEM=64
    # PROBLEM=65
    # PROBLEM=66
    # PROBLEM=67
    # PROBLEM=68
    # PROBLEM=69
    # PROBLEM=70
    # PROBLEM=71

    # PROBLEM=80
    # PROBLEM=81
    # PROBLEM=82
    # PROBLEM=83
    # PROBLEM=84
    # PROBLEM=85
    # PROBLEM=86
    # PROBLEM=87
    # PROBLEM=88
    # PROBLEM=89
    # PROBLEM=90
    # PROBLEM=91


    # PROBLEM=100
    # PROBLEM=101
    # PROBLEM=102
    # PROBLEM=103
    # PROBLEM=104
    # PROBLEM=105
    # PROBLEM=106
    # PROBLEM=107
    # PROBLEM=108
    # PROBLEM=109
    # PROBLEM=110
    # PROBLEM=111
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
# preparation for throwing DM subhalo into NW stream
if [ $PROBLEM -eq 32 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=3300.0
    INTERVAL=25.0
fi
###############################################################
# continue N-body simulation to reproduce NW stream
if [ $PROBLEM -eq 33 ]; then
    FILE=nws-continue
    CFG=nws/continue.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
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
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 80 ]; then
    FILE=nws-test-m7_0-orbit2
    CFG=nws/test/m7_0-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 81 ]; then
    FILE=nws-test-m7_0-orbit3
    CFG=nws/test/m7_0-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 82 ]; then
    FILE=nws-test-m7_5-orbit2
    CFG=nws/test/m7_5-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 83 ]; then
    FILE=nws-test-m7_5-orbit3
    CFG=nws/test/m7_5-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 84 ]; then
    FILE=nws-test-m8_0-orbit2
    CFG=nws/test/m8_0-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 85 ]; then
    FILE=nws-test-m8_0-orbit3
    CFG=nws/test/m8_0-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 86 ]; then
    FILE=nws-test-m8_5-orbit2
    CFG=nws/test/m8_5-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 87 ]; then
    FILE=nws-test-m8_5-orbit3
    CFG=nws/test/m8_5-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 88 ]; then
    FILE=nws-test-m9_0-orbit2
    CFG=nws/test/m9_0-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 89 ]; then
    FILE=nws-test-m9_0-orbit3
    CFG=nws/test/m9_0-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 90 ]; then
    FILE=nws-test-m9_5-orbit2
    CFG=nws/test/m9_5-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 91 ]; then
    FILE=nws-test-m9_5-orbit3
    CFG=nws/test/m9_5-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 100 ]; then
    FILE=nws-live-m7_0-orbit2
    CFG=nws/live/m7_0-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 101 ]; then
    FILE=nws-live-m7_0-orbit3
    CFG=nws/live/m7_0-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 102 ]; then
    FILE=nws-live-m7_5-orbit2
    CFG=nws/live/m7_5-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 103 ]; then
    FILE=nws-live-m7_5-orbit3
    CFG=nws/live/m7_5-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 104 ]; then
    FILE=nws-live-m8_0-orbit2
    CFG=nws/live/m8_0-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 105 ]; then
    FILE=nws-live-m8_0-orbit3
    CFG=nws/live/m8_0-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 106 ]; then
    FILE=nws-live-m8_5-orbit2
    CFG=nws/live/m8_5-orbit2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
    INTERVAL=25.0
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 107 ]; then
    FILE=nws-live-m8_5-orbit3
    CFG=nws/live/m8_5-orbit3.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=10700.0
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
# set input arguments
OPTION="-file=$FILE -list=$CFG -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE"
###############################################################


###############################################################
# job execution via UNIVA Grid Engine
###############################################################
# set stdout and stderr
STDOUT=log/${FILE}_$REQUEST.o$JOB_ID
STDERR=log/${FILE}_$REQUEST.e$JOB_ID
###############################################################
# load modules
. /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:/gs/hs1/jh180045/share/opt/Modules
module load intel/19.0.0.117 cuda/9.2.148 openmpi/2.1.2-opa10.9
module load phdf5
module list 1>>$STDOUT 2>>$STDERR
###############################################################
cat $PE_HOSTFILE 1>>$STDOUT 2>>$STDERR
TIME=`date`
echo "start: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################
# execute the job
if [ `which numactl` ]; then
    # run with numactl
    echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    # run without numactl
    echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################
exit 0
###############################################################
