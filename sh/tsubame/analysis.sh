#!/bin/sh
#$ -cwd
#$ -l q_core=1
#$ -l h_rt=1:00:00
#$ -N analysis
#$ -hold_jid gothic
###############################################################
# NCORES_PER_NODE=28 # for f_node
# NCORES_PER_NODE=14 # for h_node
# NCORES_PER_NODE=7 # for q_node
NCORES_PER_NODE=4 # for q_core
###############################################################
PROCS=`expr $NQUEUES \* $NCORES_PER_NODE`
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
# job execution via UNIVA Grid Engine
###############################################################
# set number of MPI processes per node
PROCS_PER_NODE=$NCORES_PER_NODE
PROCS_PER_SOCKET=14
###############################################################
# load modules
. /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:/gs/hs1/jh180045/share/opt/Modules
module load intel/19.0.0.117 cuda/9.2.148 openmpi/2.1.2-opa10.9
module load phdf5
###############################################################
# start logging
cat $PE_HOSTFILE
TIME=`date`
echo "start: $TIME"
###############################################################
# execute the job
if [ $PROCS -gt 1 ]; then
    echo "mpirun -n $PROCS sh/wrapper.sh $EXEC log/${FILE}_${REQUEST} $JOB_ID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION"
    mpirun -n $PROCS sh/wrapper.sh $EXEC log/${FILE}_${REQUEST} $JOB_ID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION
else
    # set stdout and stderr
    STDOUT=log/${FILE}_${REQUEST}.o${JOB_ID}
    STDERR=log/${FILE}_${REQUEST}.e${JOB_ID}
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
