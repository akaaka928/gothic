#!/bin/sh
if [ $# -lt 4 ]; then
    echo "$# input(s) is/are detected while at least 4 inputs are required to specify <PROCS_PER_NODE>, <PROCS_PER_SOCKET>, <MATPLOTLIB_CONFIG_DIR>, and EXEC with options" 1>&2
    exit 1
fi

# obtain rank of MPI process
RANK=${MV2_COMM_WORLD_RANK:=${PMI_RANK:=${OMPI_COMM_WORLD_RANK:=0}}}

# configuration on NUMA node
PROCS_PER_NODE=$1
PROCS_PER_SOCKET=$2
if [ `which numactl` ]; then
    TEMPID=`expr $RANK % $PROCS_PER_NODE`
    SOCKET=`expr $TEMPID / $PROCS_PER_SOCKET`
    NUMACTL="numactl --cpunodebind=$SOCKET --localalloc"
fi

# set tex cache directory for matplotlib
MPL_CFG_DIR=${3}_p${$}_r${RANK}
mkdir -p $MPL_CFG_DIR

# execute job with numactl --localalloc
EXEC="`echo $@ | sed -e "s|$3||" -e "s/$PROCS_PER_NODE//" -e "s/$PROCS_PER_SOCKET//"`"
echo "MPLCONFIGDIR=$MPL_CFG_DIR $NUMACTL python $EXEC"
MPLCONFIGDIR=$MPL_CFG_DIR $NUMACTL python $EXEC

rm -rf $MPL_CFG_DIR

exit 0
