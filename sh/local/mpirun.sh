#!/bin/sh
##################################################################
if [ $# -lt 3 ]; then
    echo "$# input(s) is/are detected while at least 3 inputs are required to specify <binary>, <# of MPI processes> and <log file>" 1>&2
    exit 1
fi
##################################################################
EXEC=$1
PROCS=$2
LOG=$3
##################################################################
OPTION="`echo $@ | sed -e "s|$EXEC||" -e "s/$PROCS//" -e "s|$LOG||"`"
##################################################################
HOSTNAME=`hostname`
HOSTFILE=host/$HOSTNAME
# PROCS=`expr $NODES \* $PROCS_PER_NODE`
##################################################################
PROCS_PER_SOCKET=2
if [ $HOSTNAME = augustus ]; then
PROCS_PER_SOCKET=1
fi
# if [ $HOSTNAME = mw ]; then
# PROCS_PER_SOCKET=1
# fi
##################################################################
if [ `which numactl` ]; then
    # run with numactl
    echo "mpirun_rsh -np $PROCS -hostfile $HOSTFILE sh/local/numarun.sh $PROCS_PER_SOCKET $LOG $EXEC $OPTION" >> $LOG
    mpirun_rsh -np $PROCS -hostfile $HOSTFILE sh/local/numarun.sh $PROCS_PER_SOCKET $LOG $EXEC $OPTION
else
    # run without numactl
    echo "mpirun_rsh -np $PROCS -hostfile $HOSTFILE $EXEC $OPTION" >> $LOG
    mpirun_rsh -np $PROCS -hostfile $HOSTFILE $EXEC $OPTION
fi
##################################################################
# # collecting data with nvprof
# echo "mpirun_rsh -np $PROCS -hostfile $HOSTFILE nvprof --timeout 4 --force-overwrite -o $LOG.$$.%q{MV2_COMM_WORLD_RANK}.nvprof $EXEC $OPTION" >> $LOG
# mpirun_rsh -np $PROCS -hostfile $HOSTFILE nvprof --timeout 4 --force-overwrite -o $LOG.$$.%q{MV2_COMM_WORLD_RANK}.nvprof $EXEC $OPTION
# # echo "mpirun_rsh -np $PROCS -hostfile $HOSTFILE nvprof --timeout 3 --device-buffer-size 128 --device-cdp-buffer-size 128 --system-profiling on --force-overwrite -o $LOG.$$.%q{MV2_COMM_WORLD_RANK}.nvprof $EXEC $OPTION" >> $LOG
# # mpirun_rsh -np $PROCS -hostfile $HOSTFILE nvprof --timeout 3 --device-buffer-size 128 --device-cdp-buffer-size 128 --system-profiling on --force-overwrite -o $LOG.$$.%q{MV2_COMM_WORLD_RANK}.nvprof $EXEC $OPTION
#################################################################
