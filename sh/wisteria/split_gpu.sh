#!/bin/sh

# read input key-value pairs and extract EXEC (with options)
EXEC=""
for arg in "$@"
do
	IFS='=' read -r key val <<< "$arg"
	case $key in
		--wrapper-Nprocs_node) PROCS_PER_NODE=${val};;
		--wrapper-Nprocs_socket) PROCS_PER_SOCKET=${val};;
		--wrapper-stdout) STDOUT=${val};;
		--wrapper-stderr) STDERR=${val};;
		--wrapper-omp_env) OMP_ENV=${val};;
		*) EXEC="$EXEC $arg";;
	esac
done

# set stdout and stderr (if necessary)
if [ -z "$STDOUT" ]; then
	STDOUT=stdout.o${$}
fi
if [ -z "$STDERR" ]; then
	STDERR=stdout.e${$}
fi

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

# # set stdout and stderr for each MPI process (if desired)
# STDOUT=${STDOUT}_${MPI_RANK}
# STDERR=${STDERR}_${MPI_RANK}

# execute job
echo "$OMP_ENV $NUMACTL $EXEC 1>>$STDOUT 2>>$STDERR"
$OMP_ENV $NUMACTL $EXEC 1>>$STDOUT 2>>$STDERR

exit 0
