#!/bin/sh

# read input key-value pairs and extract options with EXEC
EXEC=""
for arg in "$@"
do
	IFS='=' read -r key val <<< "$arg"
	case $key in
		--wrapper-Nprocs_node) PROCS_PER_NODE=${val};;
		--wrapper-Nprocs_socket) PROCS_PER_SOCKET=${val};;
		--wrapper-series) SERIES=${val};;
		--wrapper-logdir) LOGDIR=${val};;
		--wrapper-packetID) PACKET_ID=${val};;
		# --wrapper-omp_env) OMP_ENV=${val};;
		*) EXEC="$EXEC $arg";;
	esac
done

# obtain rank of MPI process
MPI_RANK=${MV2_COMM_WORLD_RANK:=${PMI_RANK:=${OMPI_COMM_WORLD_RANK:=${PMIX_RANK:=0}}}}

# set number of MPI processes per node (if necessary)
if [ -z "$PROCS_PER_NODE" ]; then
	CORES_PER_NODE=`LANG=C /usr/bin/lscpu | /usr/bin/sed -n 's/^CPU(s): *//p'`
	MPI_SIZE=${MV2_COMM_WORLD_SIZE:=${PMI_SIZE:=${OMPI_COMM_WORLD_SIZE:=${OMPI_UNIVERSE_SIZE:=0}}}}
	PROCS_PER_NODE=`expr $MPI_SIZE / $CORES_PER_NODE`
	if [ $PROCS_PER_NODE -lt 1 ]; then
		# when $MPI_SIZE < $CORES_PER_NODE; i.e. a single node can cover all MPI processes
		PROCS_PER_NODE=$MPI_SIZE
	fi
fi
DEVICE_ID=`expr $MPI_RANK % $PROCS_PER_NODE`

# configuration on NUMA node
if [ `which numactl` ]; then
	DOMAINS_PER_NODE=`LANG=C /usr/bin/lscpu | /usr/bin/sed -n 's/^NUMA node(s): *//p'`
	SOCKETS_PER_NODE=`LANG=C /usr/bin/lscpu | /usr/bin/sed -n 's/^Socket(s): *//p'`
	DOMAINS_PER_SOCKET=`expr $DOMAINS_PER_NODE / $SOCKETS_PER_NODE`

	# set number of MPI processes per socket (if necessary)
	if [ -z "$PROCS_PER_SOCKET" ]; then
		PROCS_PER_SOCKET=`expr $PROCS_PER_NODE / $SOCKETS_PER_NODE`
		if [ $PROCS_PER_SOCKET -lt 1 ]; then
			# when $PROCS_PER_NODE < $SOCKETS_PER_NODE; i.e. a single socket can cover all MPI processes in this node
			PROCS_PER_SOCKET=$PROCS_PER_NODE
		fi
	fi

    TEMPID=`expr $MPI_RANK % $PROCS_PER_NODE`
    SOCKET=`expr $TEMPID / $PROCS_PER_SOCKET`
    NUMAID=`expr $SOCKET \* $DOMAINS_PER_SOCKET`
    NUMACTL="numactl --cpunodebind=$NUMAID --localalloc"
fi


# set filename
FILE=${SERIES}_${MPI_RANK}

# set stdout and stderr
STDOUT=${LOGDIR}/${FILE}_${PACKET_ID}.log
STDERR=${LOGDIR}/${FILE}_${PACKET_ID}.log

# OMP_ENV="env KMP_AFFINITY=verbose,granularity=core,balanced" # for Intel Compilers
OMP_ENV="env OMP_DISPLAY_ENV=verbose OMP_PLACES=cores" # for GCC

# execute job
echo "rank ${MPI_RANK} on ${HOSTNAME}"
echo "$OMP_ENV $NUMACTL $EXEC -file=${FILE} -jobID=${PACKET_ID} -deviceID=${DEVICE_ID} 1>>$STDOUT 2>>$STDERR"
$OMP_ENV $NUMACTL $EXEC -file=${FILE} -jobID=${PACKET_ID} -deviceID=${DEVICE_ID} 1>>$STDOUT 2>>$STDERR

exit 0
