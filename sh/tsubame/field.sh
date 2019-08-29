#!/bin/sh
#$ -cwd
#$ -l s_gpu=1
#$ -l h_rt=12:00:00
#$ -N gothic_ext
#$ -hold_jid magi,modify,gothic_ext
###############################################################


###############################################################
# global configurations
###############################################################
if [ -z "$EXEC" ]; then
    EXEC=bin/gothic
fi
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=30
    # PROBLEM=31
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
# dynamical stability of a disk in a spherical potential field
if [ $PROBLEM -eq 0 ]; then
    FILE=disk
    SPHEPOT=va15
fi
###############################################################
# reproduction of Miki et al. (2016) in the disk coordinate system
if [ $PROBLEM -eq 11 ]; then
    FILE=m16king
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# reproduction of Kirihara et al. (2017) in the disk coordinate system
if [ $PROBLEM -eq 12 ]; then
    FILE=k17disk
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 13 ]; then
    FILE=k18nws
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# dynamical stability of an exponential disk in a spherical potential field (NFW sphere)
if [ $PROBLEM -eq 22 ]; then
    FILE=disk
    SPHEPOT=hd
fi
###############################################################
# preparation for throwing DM subhalo into NW stream
if [ $PROBLEM -eq 30 ]; then
    FILE=k18nws
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# continue N-body simulation to reproduce NW stream
if [ $PROBLEM -eq 31 ]; then
    FILE=nws-continue
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 40 ]; then
    FILE=nws-test-m7_0-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 41 ]; then
    FILE=nws-test-m7_0-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 42 ]; then
    FILE=nws-test-m7_5-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 43 ]; then
    FILE=nws-test-m7_5-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 44 ]; then
    FILE=nws-test-m8_0-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 45 ]; then
    FILE=nws-test-m8_0-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 46 ]; then
    FILE=nws-test-m8_5-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 47 ]; then
    FILE=nws-test-m8_5-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 48 ]; then
    FILE=nws-test-m9_0-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 49 ]; then
    FILE=nws-test-m9_0-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 50 ]; then
    FILE=nws-test-m9_5-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 51 ]; then
    FILE=nws-test-m9_5-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 60 ]; then
    FILE=nws-live-m7_0-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 61 ]; then
    FILE=nws-live-m7_0-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 62 ]; then
    FILE=nws-live-m7_5-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 63 ]; then
    FILE=nws-live-m7_5-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 64 ]; then
    FILE=nws-live-m8_0-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 65 ]; then
    FILE=nws-live-m8_0-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 66 ]; then
    FILE=nws-live-m8_5-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 67 ]; then
    FILE=nws-live-m8_5-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 68 ]; then
    FILE=nws-live-m9_0-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 69 ]; then
    FILE=nws-live-m9_0-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 70 ]; then
    FILE=nws-live-m9_5-orbit0
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 71 ]; then
    FILE=nws-live-m9_5-orbit1
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -jobID=$JOB_ID"
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
module load cub phdf5
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
