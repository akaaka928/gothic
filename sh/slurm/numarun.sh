#!/bin/sh
###############################################################
if [ $# -lt 3 ]; then
    echo "$# input(s) is/are detected while at least 3 inputs are required to specify <# of MPI processes per node>, <# of MPI processes per socket>, and <execution binary>" 1>&2
    exit 1
fi
###############################################################

###############################################################
PROCS_PER_NODE=$1
PROCS_PER_SOCKET=$2
EXEC=$3
###############################################################
OPTION="`echo $@ | sed -e "s/$PROCS_PER_NODE//" -e "s/$PROCS_PER_SOCKET//" -e "s|$EXEC||"`"
###############################################################

###############################################################
# set socket ID for numactl
TEMPID=`expr $MV2_COMM_WORLD_LOCAL_RANK % $PROCS_PER_NODE` # OMPI_COMM_WORLD_LOCAL_RANK for OpenMPI
SOCKET=`expr $TEMPID / $PROCS_PER_SOCKET`
###############################################################

###############################################################
numactl --cpunodebind=$SOCKET --localalloc $EXEC $OPTION
###############################################################
