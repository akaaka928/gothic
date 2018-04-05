#!/bin/bash
###############################################################
#PBS -q u-short
#PBS -l select=1:ncpus=36:mpiprocs=36:ompthreads=1
#PBS -W group_list=jh180045u
#PBS -l walltime=01:00:00
#PBS -N analysis
###############################################################
# l-small for jh180045l
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/extract
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=2
    # PROBLEM=14
    # PROBLEM=22
    # PROBLEM=112
    # PROBLEM=124
fi
###############################################################
# set number of N-body particles per bin to estimate density
if [ -z "$NCRIT" ]; then
    NCRIT=2048
fi
###############################################################
# set number of grid points for density maps
if [ -z "$NX" ]; then
    # NX=512
    NX=256
fi
if [ -z "$NY" ]; then
    # NY=512
    NY=256
fi
if [ -z "$NZ" ]; then
    # NZ=512
    NZ=256
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
    FILE=hernquist
    FINISH=23.5
    INTERVAL=0.5
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere
if [ $PROBLEM -eq 14 ]; then
    FILE=hd
    FINISH=104.0
    INTERVAL=8.0
    XMAX=5.0
    YMAX=5.0
    ZMAX=5.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere
if [ $PROBLEM -eq 22 ]; then
    FILE=disk
    FINISH=104.0
    INTERVAL=8.0
    XMAX=5.0
    YMAX=5.0
    ZMAX=5.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# reproduction of Kirihara et al. (2017) in the disk coordinate system
if [ $PROBLEM -eq 112 ]; then
    FILE=k17disk
    FINISH=1246.875
    INTERVAL=3.125
    XMAX=128.0
    YMAX=128.0
    ZMAX=128.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 124 ]; then
    FILE=halocusp_run
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
fi
###############################################################
# count up number of snapshot files
if [ -z "$START" ]; then
    START=0
fi
if [ -z "$END" ]; then
    END=`echo "scale=0; $FINISH / $INTERVAL" | bc`
fi
if [ -z "$INCREMENT" ]; then
    INCREMENT=1
fi
###############################################################
# set input arguments
OPTION="-file=$FILE -start=$START -end=$END -interval=$INCREMENT -ncrit=$NCRIT -nx=$NX -xmin=$XMIN -xmax=$XMAX -ny=$NY -ymin=$YMIN -ymax=$YMAX -nz=$NZ -zmin=$ZMIN -zmax=$ZMAX"
###############################################################


###############################################################
# job execution via PBS
###############################################################
# set number of MPI processes per node
PROCS_PER_NODE=36
PROCS_PER_SOCKET=18
###############################################################
# start logging
cd $PBS_O_WORKDIR
TIME=`date`
echo "start: $TIME"
###############################################################
export MODULEPATH=$MODULEPATH:/lustre/jh180045l/share/opt/Modules
. /etc/profile.d/modules.sh
module load intel intel-mpi
module load phdf5/impi
# module load mvapich2/gdr/2.2/gnu
# module load phdf5/mv2
###############################################################
# execute the job
if [ $PROCS -gt 1 ]; then
    echo "mpirun -n $PROCS -f ${PBS_NODEFILE} sh/wrapper.sh $EXEC log/${FILE}_${PBS_JOBNAME} $PBS_JOBID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION"
    mpirun -n $PROCS -f ${PBS_NODEFILE} sh/wrapper.sh $EXEC log/${FILE}_${PBS_JOBNAME} $PBS_JOBID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION
else
    # set stdout and stderr
    STDOUT=log/$PBS_JOBNAME.$PBS_JOBID.out
    STDERR=log/$PBS_JOBNAME.$PBS_JOBID.err
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
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
