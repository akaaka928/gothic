#!/bin/sh
if [ $# -lt 5 ]; then
    echo "$# input(s) is/are detected while at least 5 inputs are required to specify <EXEC>, <LOGNAME>, <JOB_ID>, <PROCS_PER_NODE>, and <PROCS_PER_SOCKET>" 1>&2
    exit 1
fi

# obtain rank of MPI process
RANK=${MV2_COMM_WORLD_RANK:=${PMI_RANK:=${OMPI_COMM_WORLD_RANK:=0}}}

# set stdout and stderr for each MPI process
LOGNAME=$2
JOB_ID=$3
STDOUT=${LOGNAME}.o${JOB_ID}_${RANK}
STDERR=${LOGNAME}.e${JOB_ID}_${RANK}

# configuration on NUMA node
PROCS_PER_NODE=$4
PROCS_PER_SOCKET=$5

# execute job with numactl --localalloc
EXEC=$1
OPTION="`echo $@ | sed -e "s|$EXEC||" -e "s|$LOGNAME||" -e "s/$JOB_ID//" -e "s/$PROCS_PER_NODE//" -e "s/$PROCS_PER_SOCKET//"`"
if [ `which numactl` ]; then
    TEMPID=`expr $RANK % $PROCS_PER_NODE`
    SOCKET=`expr $TEMPID / $PROCS_PER_SOCKET`

    if [ $RANK -eq 0 ]; then
	FILE=nws-live-m9_0-orbit4
	DEVID=0
    fi
    if [ $RANK -eq 1 ]; then
	FILE=nws-live-m9_0-orbit5
	DEVID=2
    fi
    if [ $RANK -eq 2 ]; then
	FILE=nws-live-m9_0-orbit6
	DEVID=1
    fi
    if [ $RANK -eq 3 ]; then
	FILE=nws-live-m9_0-orbit7
	DEVID=3
    fi

    # echo "numactl --cpunodebind=$SOCKET --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    numactl --cpunodebind=$SOCKET --localalloc $EXEC -deviceID=$DEVID -file=$FILE $OPTION 1>>$STDOUT 2>>$STDERR
else
    # echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi

exit 0
