#!/bin/bash
#PJM -g gr31
#PJM -N gothic
#PJM -L rscgrp=share
#PJM -L gpu=1
#PJM --mpi proc=1
#PJM -L elapse=3:00:00
#PJM -s


# global configurations
EXEC=bin/gothic

# problem ID
if [ -z "$PROBLEM" ]; then
    # reproduction of Komiyama et al. (2018)
    # PROBLEM=13

    # test-run for on-the-fly analysis (initial-condition)
    # PROBLEM=20

    # test-run for on-the-fly analysis (present-day)
    # continue N-body simulation to compare results with DM sub-halo free simulation
    PROBLEM=40
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
if [ $PJM_MPI_PROC -eq 1 ]; then
    NX=1
    NY=1
    NZ=1
fi
PROCS=`expr $NX \* $NY \* $NZ`
if [ $PROCS -ne $PJM_MPI_PROC ]; then
    echo "product of $NX, $NY, and $NZ must be equal to the number of total MPI processes ($PJM_MPI_PROC)"
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
# test-run for on-the-fly analysis
if [ $PROBLEM -eq 20 ]; then
    FILE=nws
fi
# continue N-body simulation to reproduce NW stream
if [ $PROBLEM -eq 40 ]; then
    FILE=nws-continue
fi

# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -Nx=$NX -Ny=$NY -Nz=$NZ -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$PJM_JOBID"


# set maximum number of MPI processes per node
PROCS_PER_NODE=8
PROCS_PER_SOCKET=4


# default module settings for Aquarius
module purge
module load cuda/11.1
module load gcc/8.3.1
module load ompi/4.1.1

# load my modules
module use /work/gr31/share/ymiki/opt/modules/lib
module load cub/1.14.0
module load phdf5/1.12.1

# start logging
cd $PJM_O_WORKDIR
TIME=`date`
echo "start: $TIME"

# execute the job
if [ $PROCS -gt 1 ]; then
    echo "mpiexec sh/wrapper.sh $EXEC log/${FILE}_${PJM_JOBNAME} $PJM_JOBID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION"
    mpiexec sh/wrapper.sh $EXEC log/${FILE}_${PJM_JOBNAME} $PJM_JOBID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION
else
    # set stdout and stderr
    STDOUT=log/${FILE}_${PJM_JOBNAME}.o${PJM_JOBID}
    STDERR=log/${FILE}_${PJM_JOBNAME}.e${PJM_JOBID}
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
