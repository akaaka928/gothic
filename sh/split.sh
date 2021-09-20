#!/bin/sh
if [ $# -lt 6 ]; then
    echo "$# input(s) is/are detected while at least 6 inputs are required to specify <EXEC>, <LOGNAME>, <JOB_ID>, <PROCS_PER_NODE>, <PROCS_PER_SOCKET>, and <LIST>" 1>&2
    exit 1
fi


RANK=${MV2_COMM_WORLD_RANK:=${PMI_RANK:=${OMPI_COMM_WORLD_RANK:=0}}}
SIZE=${MV2_COMM_WORLD_SIZE:=${PMI_SIZE:=${OMPI_COMM_WORLD_SIZE:=0}}}

# set stdout and stderr for each MPI process
EXEC=$1
LOGNAME=$2
JOB_ID=$3
STDOUT=${LOGNAME}.o${JOB_ID}_${RANK}
STDERR=${LOGNAME}.e${JOB_ID}_${RANK}
# STDOUT=${LOGNAME}.o${JOB_ID}
# STDERR=${LOGNAME}.e${JOB_ID}

# configuration on NUMA node
PROCS_PER_NODE=$4
PROCS_PER_SOCKET=$5


MAX=`expr $SIZE - 1`
digit=${#MAX}
int="`echo ${#RANK}`"
if [ "$int" -le "$digit" ]; then
    rem=`expr $digit - $int`
    zeros=""
    count=0
    while [ "$count" -lt "$rem" ]; do
	zeros="`echo ${zeros}0`"
	count=`expr $count + 1`
    done
fi
NAME=${6}
LIST=${NAME}${zeros}${RANK}


OPTION="`echo $@ | sed -e "s|$EXEC||" -e "s|$LOGNAME||" -e "s/$JOB_ID//" -e "s/$PROCS_PER_NODE//" -e "s/$PROCS_PER_SOCKET//" -e "s|$NAME||"`"
if [ `which numactl` ]; then
    TEMPID=`expr $RANK % $PROCS_PER_NODE`
    SOCKET=`expr $TEMPID / $PROCS_PER_SOCKET`
    while read LINE; do
	FILE=${LINE}
	numactl --cpunodebind=$SOCKET --localalloc $EXEC -file=$FILE $OPTION -deviceID=$TEMPID 1>>$STDOUT 2>>$STDERR
    done < $LIST
    rm $LIST
else
    while read LINE; do
	FILE=${LINE}
	$EXEC -file=$FILE $OPTION -deviceID=$TEMPID 1>>$STDOUT 2>>$STDERR
    done < $LIST
    rm $LIST
fi

exit 0
