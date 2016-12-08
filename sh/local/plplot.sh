#!/bin/sh
##################################################################
# http://program.station.ez-net.jp/special/handbook/sh/variables/name-ext.asp
##################################################################
if [ $# -lt 3 ]; then
    echo "$# inputs are detected while at least 3 inputs are required to execute this script" 1>&2
    exit 1
fi
##################################################################
EXEC=$1
FILENAME=$2
EXTENSION=$3
##################################################################
OPTION="`echo $@ | sed -e "s|$EXEC||" -e "s/$FILENAME//" -e "s/$EXTENSION//"`"
##################################################################
HOSTNAME=`hostname`
HOSTFILE=host/$HOSTNAME
NODES=2
NCORE_PER_NODE=8
if [ $HOSTNAME = augustus ]; then
NODES=1
NCORE_PER_NODE=4
fi
##################################################################
NCORE=`expr $NODES \* $NCORE_PER_NODE`
##################################################################
mpirun -np $NCORE -hostfile $HOSTFILE $EXEC -file=$FILENAME $OPTION -dev ${EXTENSION}cairo
##################################################################
