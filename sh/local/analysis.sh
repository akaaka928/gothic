#!/bin/sh
###############################################################
if [ $# -lt 1 ]; then
    echo "$# input(s) is/are detected while at least 1 input is required to specify <test problem>" 1>&2
    exit 1
fi
JOB_ID=$$
PROBLEM=$1
###############################################################
# set number of MPI processes per node
PROCS=16
if [ $# -gt 2 ]; then
    PROCS=$2
fi
PROCS_PER_NODE=32
PROCS_PER_SOCKET=16
###############################################################
#
#
###############################################################
# global configurations
###############################################################
EXEC=bin/extract
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=2
    PROBLEM=10
    # PROBLEM=14
    # PROBLEM=22
    # PROBLEM=26
    # PROBLEM=27
    # PROBLEM=62
    # PROBLEM=130
    # PROBLEM=131
    # PROBLEM=132
    # PROBLEM=133
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
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
    FILE=hernquist
    FINISH=23.5
    INTERVAL=0.5
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    VMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
fi
###############################################################
# dynamical stability of a King sphere within an NFW halo
if [ $PROBLEM -eq 10 ]; then
    FILE=hb
    FINISH=23.5
    INTERVAL=0.5
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    VMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
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
    VMAX=5.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
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
    VMAX=5.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
fi
###############################################################
# dynamical stability of multi components galaxy model (only spherical components)
if [ $PROBLEM -eq 26 ]; then
    FILE=etg
    # FINISH=75.0
    FINISH=375.0
    # FINISH=1175.0
    # FINISH=3175.0
    INTERVAL=25.0
    XMAX=5.0
    YMAX=5.0
    ZMAX=5.0
    VMAX=5.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
fi
###############################################################
# dynamical stability of an M31 model (NFW halo, Hernquist bulge, and exponential disk)
# basically, this is Fardal et al. (2007) model
# stellar halo: Gilbert et al. (2012): \Sigma \propto R^-2.2; Rmin = 9kpc, Rmax = 176kpc; Ibata et al. (2014, ApJ, 780, 128): total stellar mass of the smooth halo is ~8e+9 Msun
# disk: Toomre's Q-value is set to reproduce Tenjes et al. (2017): Q_min = 1.8 @ 12-13 kpc
if [ $PROBLEM -eq 27 ]; then
    FILE=m31
    # FINISH=75.0
    FINISH=375.0
    # FINISH=1175.0
    # FINISH=3175.0
    INTERVAL=25.0
    XMAX=15.0
    YMAX=15.0
    ZMAX=15.0
    VMAX=400.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
fi
###############################################################
# dynamical stability of a progenitor model for NW stream determined by Komiyama et al. (2018)
if [ $PROBLEM -eq 62 ]; then
    FILE=satellite
    INTERVAL=100.0
    FINISH=14000.0
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    VMAX=5.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 130 ]; then
    FILE=halocusp_run
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=3.0
x    YMAX=3.0
    ZMAX=3.0
    VMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 131 ]; then
    FILE=halocore1_run
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    VMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 132 ]; then
    FILE=halocore2_run
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    VMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 133 ]; then
    FILE=halocore3_run
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=3.0
    YMAX=3.0
    ZMAX=3.0
    VMAX=3.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
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
OPTION="-file=$FILE -start=$START -end=$END -interval=$INCREMENT -ncrit=$NCRIT -nx=$NX -xmin=$XMIN -xmax=$XMAX -ny=$NY -ymin=$YMIN -ymax=$YMAX -nz=$NZ -zmin=$ZMIN -zmax=$ZMAX -nv=$NV -vmin=$VMIN -vmax=$VMAX -nx3D=$NX3D -ny3D=$NY3D -nz3D=$NZ3D"
###############################################################
#
#
###############################################################
# start logging
TIME=`date`
echo "start: $TIME"
###############################################################
# execute the job
echo "mpiexec -n $PROCS sh/wrapper.sh $EXEC log/${FILE}_analysis $JOB_ID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION"
mpiexec -n $PROCS sh/wrapper.sh $EXEC log/${FILE}_analysis $JOB_ID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
