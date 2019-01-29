#!/bin/sh
if [ $# -lt 5 ]; then
    echo "$# input(s) is/are detected while at least 5 inputs are required to specify <EXEC>, <FILENAME>, <JOB_ID>, <PROCS_PER_NODE>, and <PROCS_PER_SOCKET>" 1>&2
    exit 1
fi

# obtain rank of MPI process
RANK=${MV2_COMM_WORLD_RANK:=${PMI_RANK:=${OMPI_COMM_WORLD_RANK:=0}}}

# assign jobs
FILENAME=$2
if [ $RANK -eq 0 ]; then
    FILE=${FILENAME}iso
fi
if [ $RANK -eq 1 ]; then
    FILE=${FILENAME}ra2
fi
if [ $RANK -eq 2 ]; then
    FILE=${FILENAME}ra1
fi
if [ $RANK -eq 3 ]; then
    FILE=${FILENAME}ra1_2
fi
if [ $RANK -eq 4 ]; then
    FILE=${FILENAME}ra1_4
fi
if [ $RANK -eq 5 ]; then
    FILE=${FILENAME}ra1_8
fi

# set stdout and stderr for each MPI process
JOB_ID=$3
LOGNAME=log/${FILE}_gothic
STDOUT=${LOGNAME}.o${JOB_ID}
STDERR=${LOGNAME}.e${JOB_ID}

# configuration on NUMA node
PROCS_PER_NODE=$4
PROCS_PER_SOCKET=$5

# execute job with numactl --localalloc
EXEC=$1
OPTION="`echo $@ | sed -e "s|$EXEC||" -e "s|$LOGNAME||" -e "s/$JOB_ID//" -e "s/$PROCS_PER_NODE//" -e "s/$PROCS_PER_SOCKET//"`"
if [ `which numactl` ]; then
    TEMPID=`expr $RANK % $PROCS_PER_NODE`
    SOCKET=`expr $TEMPID / $PROCS_PER_SOCKET`
    numactl --cpunodebind=$SOCKET --localalloc $EXEC -file=$FILE -deviceID=$RANK $OPTION 1>>$STDOUT 2>>$STDERR
else
    $EXEC -file=$FILE -deviceID=$RANK $OPTION 1>>$STDOUT 2>>$STDERR
fi

exit 0
