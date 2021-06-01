#!/bin/sh

# read input key-value pairs and extract options
OPTION=""
for arg in "$@"
do
	IFS='=' read -r key val <<< "$arg"
	case $key in
		--wrapper-Nprocs_node) PROCS_PER_NODE=${val};;
		--wrapper-Nprocs_socket) PROCS_PER_SOCKET=${val};;
		--wrapper-EXEC-list) EXEC_LIST=${val};;
		--wrapper-packetID) PACKET_ID=${val};;
		--wrapper-digits) digit=${val};;
		--wrapper-series) SERIES=${val};;
		--wrapper-logdir) LOGDIR=${val};;
		--wrapper-omp_env) OMP_ENV=${val};;
		--wrapper-index) INDEX=${val};;
		*) OPTION="$OPTION $arg";;
	esac
done

# obtain rank of MPI process
MPI_RANK=${MV2_COMM_WORLD_RANK:=${PMI_RANK:=${OMPI_COMM_WORLD_RANK:=${PMIX_RANK:=0}}}}

# set number of MPI processes per node (if necessary)
if [ -z "$PROCS_PER_NODE" ]; then
	CORES_PER_NODE=`grep processor /proc/cpuinfo | wc -l`
	MPI_SIZE=${MV2_COMM_WORLD_SIZE:=${PMI_SIZE:=${OMPI_COMM_WORLD_SIZE:=${OMPI_UNIVERSE_SIZE:=0}}}}
	PROCS_PER_NODE=`expr $MPI_SIZE / $CORES_PER_NODE`
	if [ $PROCS_PER_NODE -lt 1 ]; then
		# when $MPI_SIZE < $CORES_PER_NODE; i.e. a single node can cover all MPI processes
		PROCS_PER_NODE=$MPI_SIZE
	fi
fi
DEVICE_ID=`expr $MPI_RANK % $PROCS_PER_NODE`

# set number of MPI processes per socket (if necessary)
if [ -z "$PROCS_PER_SOCKET" ]; then
	SOCKETS_PER_NODE=`lscpu | grep NUMA | grep CPU | wc -l`
	PROCS_PER_SOCKET=`expr $PROCS_PER_NODE / $SOCKETS_PER_NODE`
	if [ $PROCS_PER_SOCKET -lt 1 ]; then
		# when $PROCS_PER_NODE < $SOCKETS_PER_NODE; i.e. a single socket can cover all MPI processes in this node
		PROCS_PER_SOCKET=$PROCS_PER_NODE
	fi
fi

# configuration on NUMA node
if [ `which numactl` ]; then
    TEMPID=`expr $MPI_RANK % $PROCS_PER_NODE`
    SOCKET=`expr $TEMPID / $PROCS_PER_SOCKET`
    NUMACTL="numactl --cpunodebind=$SOCKET --localalloc"
fi


# pick up the EXEC
int="`echo ${#PACKET_ID}`"
if [ "$int" -le "$digit" ]; then
    rem=`expr $digit - $int`
    zeros=""
    count=0
    while [ "$count" -lt "$rem" ]; do
	zeros="`echo ${zeros}0`"
	count=`expr $count + 1`
    done
fi
LIST=${EXEC_LIST}${zeros}${PACKET_ID}
LINE=`expr $MPI_RANK + 1`
EXEC=`sed -n ${LINE}P ${LIST}`

# set filename
FILE=${SERIES}_${MPI_RANK}

# set stdout and stderr
STDOUT=${LOGDIR}/${FILE}_${INDEX}.log
STDERR=${LOGDIR}/${FILE}_${INDEX}.log

# execute job
echo "rank ${MPI_RANK} on ${HOSTNAME}"
echo "$OMP_ENV $NUMACTL $EXEC -file=${FILE} -jobID=${INDEX} -deviceID=${DEVICE_ID} $OPTION 1>>$STDOUT 2>>$STDERR"
$OMP_ENV $NUMACTL $EXEC -file=${FILE} -jobID=${INDEX} -deviceID=${DEVICE_ID} $OPTION 1>>$STDOUT 2>>$STDERR

exit 0
