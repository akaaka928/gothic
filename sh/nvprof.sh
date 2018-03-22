#!/bin/sh
if [ $# -lt 6 ]; then
    echo "$# input(s) is/are detected while at least 5 inputs are required to specify <EXEC>, <LOGNAME>, <JOB_ID>, <PROCS_PER_NODE>, <PROCS_PER_SOCKET>, and <NVPROF_TIMEOUT>" 1>&2
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

NVPROF_TIMEOUT=$6

# execute job with numactl --localalloc
EXEC=$1
OPTION="`echo $@ | sed -e "s|$EXEC||" -e "s|$LOGNAME||" -e "s/$JOB_ID//" -e "s/$PROCS_PER_NODE//" -e "s/$PROCS_PER_SOCKET//" -e "s/$NVPROF_TIMEOUT//"`"
if [ `which numactl` ]; then
    TEMPID=`expr $RANK % $PROCS_PER_NODE`
    SOCKET=`expr $TEMPID / $PROCS_PER_SOCKET`
    echo "numactl --cpunodebind=$SOCKET --localalloc nvprof --timeout $NVPROF_TIMEOUT --force-overwrite -o ${LOGNAME}.p${JOB_ID}_r${RANK}.nvprof $EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    numactl --cpunodebind=$SOCKET --localalloc nvprof --timeout $NVPROF_TIMEOUT --force-overwrite -o ${LOGNAME}.p${JOB_ID}_r${RANK}.nvprof $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    echo "nvprof --timeout $NVPROF_TIMEOUT --force-overwrite -o ${LOGNAME}.p${JOB_ID}_r${RANK}.nvprof $EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    nvprof --timeout $NVPROF_TIMEOUT --force-overwrite -o ${LOGNAME}.p${JOB_ID}_r${RANK}.nvprof $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi

exit 0
