#!/bin/sh
##################################################################
if [ $# -lt 3 ]; then
    echo "$# input(s) is/are detected while at least 3 inputs are required to specify <# of MPI processes per socket>, <execution binary> and <log file>" 1>&2
    exit 1
fi
##################################################################
PROCS_PER_SOCKET=$1
LOG=$2
EXEC=$3
OPTION="`echo $@ | sed -e "s/$PROCS_PER_SOCKET//" -e "s|$LOG||" -e "s|$EXEC||"`"
##################################################################
# Machine specific values
SOCKETS_PER_NODE=2
HOSTNAME=`hostname`
if [ $HOSTNAME = augustus ]; then
SOCKETS_PER_NODE=1
fi
##################################################################
PROCS_PER_NODE=`expr $PROCS_PER_SOCKET \* $SOCKETS_PER_NODE`
##################################################################
# implementation for mvapich2
MYRANK=0
TEMPID=0
SOCKET=0
if [ $PROCS_PER_SOCKET -gt 1 ]; then
    MYRANK=$MV2_COMM_WORLD_LOCAL_RANK
    TEMPID=`expr $MYRANK % $PROCS_PER_NODE`
    SOCKET=`expr $TEMPID / $PROCS_PER_SOCKET`
fi
# echo "rank $MYRANK: numactl --cpunodebind=$SOCKET --interleave=$SOCKET $EXEC $OPTION"
# exit 1
##################################################################
# echo "rank $MYRANK: numactl --cpunodebind=$SOCKET --interleave=$SOCKET $EXEC $OPTION" >> $LOG
# numactl --cpunodebind=$SOCKET --interleave=$SOCKET $EXEC $OPTION
echo "rank $MYRANK: numactl --cpunodebind=$SOCKET --localalloc $EXEC $OPTION" >> $LOG
numactl --cpunodebind=$SOCKET --localalloc $EXEC $OPTION
##################################################################
