#!/bin/bash
#$ -cwd
#$ -l f_node=1
#$ -l h_rt=00:10:00
#$ -N m31obs
###############################################################
NCPUS_PER_NODE=28
# NGPUS_PER_NODE=14
# NGPUS_PER_NODE=7
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/m31obs
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=112
    PROBLEM=113
fi
###############################################################
# set number of N-body particles per bin to estimate density
if [ -z "$NCRIT" ]; then
    NCRIT=2048
fi
###############################################################
# set number of grid points for density maps
if [ -z "$NX3D" ]; then
    NX3D=256
fi
if [ -z "$NY3D" ]; then
    NY3D=256
fi
if [ -z "$NZ3D" ]; then
    NZ3D=256
fi
if [ -z "$NX" ]; then
    NX=1024
fi
if [ -z "$NY" ]; then
    NY=1024
fi
if [ -z "$NZ" ]; then
    NZ=1024
fi
if [ -z "$NV" ]; then
    NV=1024
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# reproduction of Kirihara et al. (2017) in the disk coordinate system
if [ $PROBLEM -eq 112 ]; then
    FILE=k17disk
    FINISH=1246.875
    INTERVAL=3.125
    XIMIN=-5.0
    XIMAX=7.0
    ETAMIN=-8.5
    ETAMAX=3.5
    DMIN=-80.0
    DMAX=130.0
    VMIN=-700.0
    VMAX=-100.0
fi
###############################################################
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 113 ]; then
    FILE=k18nws
    FINISH=14000.0
    INTERVAL=25.0
    XIMIN=-9.0
    XIMAX=7.0
    ETAMIN=-3.0
    ETAMAX=13.0
    DMIN=-50.0
    DMAX=200.0
    VMIN=-540.0
    VMAX=-200.0
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
OPTION="-file=$FILE -start=$START -end=$END -interval=$INCREMENT -ncrit=$NCRIT -nxi=$NX -ximin=$XIMIN -ximax=$XIMAX -neta=$NY -etamin=$ETAMIN -etamax=$ETAMAX -nD=$NZ -Dmin=$DMIN -Dmax=$DMAX -nvlos=$NV -vlosmin=$VMIN -vlosmax=$VMAX -nx3D=$NX3D -ny3D=$NY3D -nz3D=$NZ3D"
###############################################################


###############################################################
# job execution via UNIVA Grid Engine
###############################################################
# set number of MPI processes
PROCS=`expr $NQUEUES \* $NCPUS_PER_NODE`
###############################################################
# set stdout and stderr
STDOUT=log/${FILE}_$REQUEST.o$JOB_ID
STDERR=log/${FILE}_$REQUEST.e$JOB_ID
###############################################################
TIME=`date`
echo "start: $TIME"
###############################################################
# load modules
. /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:/gs/hs1/jh180045/share/opt/Modules
module load intel/19.0.0.117 cuda/9.2.148 openmpi/2.1.2-opa10.9
module load cub phdf5
module list
###############################################################
cat $PE_HOSTFILE
###############################################################
# execute the job
mpirun -n $PROCS -npernode $NCPUS_PER_NODE -x LD_LIBRARY_PATH $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
