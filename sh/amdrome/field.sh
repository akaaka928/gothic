#!/bin/bash
#SBATCH -J gothic             # name of job
#SBATCH -t 24:00:00           # upper limit of elapsed time
#SBATCH -p amdrome            # partition name
#SBATCH -w amd2     # use compute node equips NVIDIA A100 PCIe
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=1
##SBATCH --gpu-freq=1410 # requested GPU frequency: 210 MHz--1410 MHz (bin = 15 MHz)
#SBATCH --mem=2048 # in units of MB


if [ -z "$SLURM_NTASKS" ]; then
    SLURM_NTASKS=1
fi


# global configurations
EXEC=bin/gothic

# problem ID
if [ -z "$PROBLEM" ]; then
    # reproduction of Komiyama et al. (2018)
    PROBLEM=13

    # continue N-body simulation to compare results with DM sub-halo free simulation
    # PROBLEM=40
fi

# topology of MPI processes
if [ -z "$NX" ]; then
    NX=2
fi
if [ -z "$NY" ]; then
    NY=1
fi
if [ -z "$NZ" ]; then
    NZ=1
fi
if [ $SLURM_NTASKS -eq 1 ]; then
    NX=1
    NY=1
    NZ=1
fi
PROCS=`expr $NX \* $NY \* $NZ`
if [ $PROCS -ne $SLURM_NTASKS ]; then
    echo "product of $NX, $NY, and $NZ must be equal to the number of total MPI processes ($SLURM_NTASKS)"
    exit 1
fi

# value of accuracy controling parameter: GADGET MAC by Springel (2005)
if [ -z "$ABSERR" ]; then
    # ABSERR=1.250000000e-1
    # ABSERR=6.250000000e-2
    # ABSERR=3.125000000e-2
    # ABSERR=1.562500000e-2
    # ABSERR=7.812500000e-3
    # ABSERR=3.906250000e-3
    ABSERR=1.953125000e-3
    # ABSERR=9.765625000e-4
    # ABSERR=4.882812500e-4
    # ABSERR=2.441406250e-4
    # ABSERR=1.220703125e-4
    # ABSERR=6.103515625e-5
fi

if [ -z "$REBUILD" ]; then
    REBUILD=16
fi

if [ -z "$BRENT" ]; then
    BRENT=1.0
fi


if [ -z "$SPHEPOT" ]; then
    SPHEPOT=m31
fi

if [ -z "$DISKPOT" ]; then
    DISKPOT=m31
fi


## problem specific configurations
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 13 ]; then
    FILE=k18nws
fi
# continue N-body simulation to reproduce NW stream
if [ $PROBLEM -eq 40 ]; then
    FILE=nws-continue
fi

# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -Nx=$NX -Ny=$NY -Nz=$NZ -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$SLURM_JOB_ID"


# job execution via SLURM
# set number of MPI processes per node
PROCS_PER_NODE=`expr $SLURM_NTASKS / $SLURM_JOB_NUM_NODES`

export MODULEPATH=$HOME/opt/modules:$MODULEPATH
module purge
module load cuda
module load openmpi
module load phdf5

# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_NUM_NODES nodes"
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs per node"
TIME=`date`
echo "start: $TIME"

# execute the job
if [ $PROCS -gt 1 ]; then
    echo "mpiexec -n $SLURM_NTASKS sh/wrapper.sh $EXEC log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTION"
    mpiexec -n $SLURM_NTASKS sh/wrapper.sh $EXEC log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTION
else
    # set stdout and stderr
    STDOUT=log/${FILE}_$SLURM_JOB_NAME.o${SLURM_JOB_ID}
    STDERR=log/${FILE}_$SLURM_JOB_NAME.e${SLURM_JOB_ID}
    if [ `which numactl` ]; then
	# run with numactl
	echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
	numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
    else
	# run without numactl
	echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
	$EXEC $OPTION 1>>$STDOUT 2>>$STDERR
    fi
fi

# finish logging
TIME=`date`
echo "finish: $TIME"

exit 0
