#!/bin/sh
###############################################################
if [ $# -lt 1 ]; then
    echo "$# input(s) is/are detected while at least 1 input is required to specify <test problem>" 1>&2
    exit 1
fi
JOB_ID=$$
PROBLEM=$1
###############################################################
NX=1
NY=1
NZ=1
if [ $# -eq 4 ]; then
    NX=$2
    NY=$3
    NZ=$4
fi
PROCS=`expr $NX \* $NY \* $NZ`
###############################################################
DO_MEMCHK=0
if [ $# -eq 2 ]; then
    if [ $PROCS -eq 1 ]; then
	DO_MEMCHK=$2
    fi
fi
if [ $# -eq 5 ]; then
    DO_MEMCHK=$5
fi
###############################################################
PROCS_PER_NODE=4
PROCS_PER_SOCKET=2
###############################################################
#
#
###############################################################
# global configurations
###############################################################
MEMCHK=cuda-memcheck
EXE=bin/gothic
###############################################################
# accuracy contoling parameters
###############################################################
# GADGET MAC by Springel (2005)
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
#
#
###############################################################
# specify the target problem
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 100 ]; then
    FILE=m12iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 101 ]; then
    FILE=m12ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 102 ]; then
    FILE=m12ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 103 ]; then
    FILE=m12ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 104 ]; then
    FILE=m12ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 105 ]; then
    FILE=m12ra1_8_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 110 ]; then
    FILE=m13iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 111 ]; then
    FILE=m13ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 112 ]; then
    FILE=m13ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 113 ]; then
    FILE=m13ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 114 ]; then
    FILE=m13ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 120 ]; then
    FILE=m14iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 121 ]; then
    FILE=m14ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 122 ]; then
    FILE=m14ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 123 ]; then
    FILE=m14ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 124 ]; then
    FILE=m14ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 130 ]; then
    FILE=m15iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 131 ]; then
    FILE=m15ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 132 ]; then
    FILE=m15ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 133 ]; then
    FILE=m15ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 134 ]; then
    FILE=m15ra1_4_run
    SINK=1
    COM_GROUP_HEAD=0
    COM_GROUP_NUM=16777216
    COM_RESET_INTERVAL=8
    REPLACE_CENTRAL_BH=1
    CENTRAL_BH_ID=16777216
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 140 ]; then
    FILE=m11iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 141 ]; then
    FILE=m11ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 142 ]; then
    FILE=m11ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 143 ]; then
    FILE=m11ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 144 ]; then
    FILE=m11ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 145 ]; then
    FILE=m11ra1_8_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 150 ]; then
    FILE=m10iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 151 ]; then
    FILE=m10ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 152 ]; then
    FILE=m10ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 153 ]; then
    FILE=m10ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 154 ]; then
    FILE=m10ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 155 ]; then
    FILE=m10ra1_8_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 160 ]; then
    FILE=m09iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 161 ]; then
    FILE=m09ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 162 ]; then
    FILE=m09ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 163 ]; then
    FILE=m09ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 164 ]; then
    FILE=m09ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 165 ]; then
    FILE=m09ra1_8_run
    SINK=1
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -Nsink=$SINK -Nx=$NX -Ny=$NY -Nz=$NZ -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$JOB_ID -com_group_head=$COM_GROUP_HEAD -com_group_num=$COM_GROUP_NUM -com_reset_interval=$COM_RESET_INTERVAL -replace_central_bh=$REPLACE_CENTRAL_BH -central_bh_id=$CENTRAL_BH_ID"
###############################################################
#
#
###############################################################
# execute numerical simulation
###############################################################
# LOG=$JOB_NAME.o$JOB_ID
LOG=log/$FILE.l
STDOUT=log/$FILE.o$JOB_ID
STDERR=log/$FILE.e$JOB_ID
TIME=`date`
echo "$TIME: $EXE start" >> $LOG
if [ $PROCS -eq 1 ]; then
    if [ $DO_MEMCHK -eq 0 ]; then
	if [ `which numactl` ]; then
	    echo "numactl --localalloc $EXE $OPTION 1>>$STDOUT 2>>$STDERR" >> $LOG
	    numactl --localalloc $EXE $OPTION 1>>$STDOUT 2>>$STDERR
	else
	    # run without numactl
	    echo "$EXE $OPTION 1>>$STDOUT 2>>$STDERR" >>$LOG
	    $EXE $OPTION 1>>$STDOUT 2>>$STDERR
	fi
    else
	echo "$MEMCHK $EXE $OPTION 1>>$STDOUT 2>>$STDERR" >> $LOG
	$MEMCHK $EXE $OPTION 1>>$STDOUT 2>>$STDERR
    fi
else
    if [ $DO_MEMCHK -eq 0 ]; then
	echo "mpiexec -n $PROCS sh/wrapper.sh $EXE log/$FILE $JOB_ID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION" >> $LOG
	mpiexec -n $PROCS sh/wrapper.sh $EXE log/$FILE $JOB_ID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION
    else
	echo "$MEMCHK $MPIEXE $EXE $PROCS $LOG $OPTION 1>>$STDOUT 2>>$STDERR" >> $LOG
	$MEMCHK $MPIEXE $EXE $PROCS $LOG $OPTION 1>>$STDOUT 2>>$STDERR
    fi
fi
TIME=`date`
echo "$TIME: $EXE finish" >> $LOG
###############################################################
