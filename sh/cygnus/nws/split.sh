#!/bin/sh
if [ $# -lt 6 ]; then
    echo "$# input(s) is/are detected while at least 6 inputs are required to specify <EXEC>, <SERIES>, <LOGNAME>, <JOB_ID>, <PROCS_PER_NODE>, and <PROCS_PER_SOCKET>" 1>&2
    exit 1
fi

# obtain rank of MPI process
RANK=${MV2_COMM_WORLD_RANK:=${PMI_RANK:=${OMPI_COMM_WORLD_RANK:=0}}}

# set stdout and stderr for each MPI process
SERIES=$2
LOGNAME=$3
JOB_ID=$4
STDOUT=${LOGNAME}.o${JOB_ID}_${RANK}
STDERR=${LOGNAME}.e${JOB_ID}_${RANK}

# configuration on NUMA node
PROCS_PER_NODE=$5
PROCS_PER_SOCKET=$6

# execute job with numactl --localalloc
EXEC=$1
OPTION="`echo $@ | sed -e "s|$EXEC||" -e "s|$SERIES||" -e "s|$LOGNAME||" -e "s/$JOB_ID//" -e "s/$PROCS_PER_NODE//" -e "s/$PROCS_PER_SOCKET//"`"


TEMPID=`expr $RANK % $PROCS_PER_NODE`
# SOCKET=`expr $TEMPID / $PROCS_PER_SOCKET`

FILE=${SERIES}_orbit${RANK}
DEVID=${TEMPID}
numactl --cpunodebind=0 --localalloc $EXEC -deviceID=$DEVID -file=$FILE $OPTION 1>>$STDOUT 2>>$STDERR

exit 0
