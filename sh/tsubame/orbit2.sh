#!/bin/sh
#$ -cwd
#$ -l s_gpu=1
#$ -l h_rt=24:00:00
#$ -N nws-test-m7^0-orbit2
#$ -hold_jid nws-test-m7^0-orbit2
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
    PROBLEM=80
    # PROBLEM=82
    # PROBLEM=84
    # PROBLEM=86
    # PROBLEM=88
    # PROBLEM=90

    # PROBLEM=100
    # PROBLEM=102
    # PROBLEM=104
    # PROBLEM=106
    # PROBLEM=108
    # PROBLEM=110
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
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 80 ]; then
    FILE=nws-test-m7_0-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 82 ]; then
    FILE=nws-test-m7_5-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 84 ]; then
    FILE=nws-test-m8_0-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 86 ]; then
    FILE=nws-test-m8_5-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 88 ]; then
    FILE=nws-test-m9_0-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (test-particle) and NW stream
if [ $PROBLEM -eq 90 ]; then
    FILE=nws-test-m9_5-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 100 ]; then
    FILE=nws-live-m7_0-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 102 ]; then
    FILE=nws-live-m7_5-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 104 ]; then
    FILE=nws-live-m8_0-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 106 ]; then
    FILE=nws-live-m8_5-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 108 ]; then
    FILE=nws-live-m9_0-orbit2
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# collision between DM subhalo (N-body) and NW stream
if [ $PROBLEM -eq 110 ]; then
    FILE=nws-live-m9_5-orbit2
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
