#!/bin/sh
#$ -cwd
#$ -l s_gpu=1
#$ -l h_rt=24:00:00
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
    # reproduction of Komiyama et al. (2018)
    PROBLEM=13

    # continue N-body simulation to compare results with DM sub-halo free simulation
    # PROBLEM=40
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
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 13 ]; then
    FILE=k18nws
fi
###############################################################
# continue N-body simulation to reproduce NW stream
if [ $PROBLEM -eq 40 ]; then
    FILE=nws-continue
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$JOB_ID"
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
