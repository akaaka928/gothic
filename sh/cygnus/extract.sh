#!/bin/bash
#PBS -A GALAXY
#PBS -q gpu
#PBS -N analysis
#PBS -b 1
#PBS -l elapstim_req=01:00:00
#PBS -T mvapich
#PBS -v NQSV_MPI_VER=2.3.1/intel-cuda10.1
###############################################################
NUM_NODES=1
PROCS_PER_SOCKET=12
PROCS_PER_NODE=`expr $PROCS_PER_SOCKET \* 2`
# PROCS_PER_NODE=$PROCS_PER_SOCKET
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/extract
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=2
    # PROBLEM=10
    # PROBLEM=14
    # PROBLEM=22
    # PROBLEM=26
    PROBLEM=27
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
if [ -z "$NE" ]; then
    NE=1024
fi
if [ -z "$NJ" ]; then
    NJ=1024
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
fi
###############################################################
# dynamical stability of an M31 model (NFW halo, Hernquist bulge, and exponential disk)
# basically, this is Fardal et al. (2007) model
# stellar halo: Gilbert et al. (2012): \Sigma \propto R^-2.2; Rmin = 9kpc, Rmax = 176kpc; Ibata et al. (2014, ApJ, 780, 128): total stellar mass of the smooth halo is ~8e+9 Msun
# disk: Toomre's Q-value is set to reproduce Tenjes et al. (2017): Q_min = 1.8 @ 12-13 kpc
if [ $PROBLEM -eq 27 ]; then
    FILE=m31
    # FINISH=75.0
    # FINISH=375.0
    FINISH=1175.0
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 130 ]; then
    FILE=halocusp_run
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
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
    EMIN=-5.0e+5 # minmax.txt x100
    EMAX=-0.0e+3 # minmax.txt x100
    JMIN=3.1e-3 # minmax.txt x10
    JMAX=3.1e+4 # minmax.txt x10
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 134 ]; then
    FILE=m15ra1_4_run
    FINISH=2000.0
    INTERVAL=25.0
    XMAX=1000.0
    YMAX=1000.0
    ZMAX=1000.0
    VMAX=1000.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.0e+7
    EMAX=-3.1e+5
    JMIN=8.0e-1
    JMAX=7.0e+5
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 200 ]; then
    FILE=w3iso
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.15
    YMAX=0.15
    ZMAX=0.15
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+2 # minmax.txt x100
    EMAX=-2.0e+1 # minmax.txt x100
    JMIN=1.0e-4 # minmax.txt x10
    JMAX=7.0e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 201 ]; then
    FILE=w3ra2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.15
    YMAX=0.15
    ZMAX=0.15
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+2 # minmax.txt x100
    EMAX=-2.5e+1 # minmax.txt x100
    JMIN=1.0e-4 # minmax.txt x10
    JMAX=5.0e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 202 ]; then
    FILE=w3ra1
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.15
    YMAX=0.15
    ZMAX=0.15
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+2 # minmax.txt x100
    EMAX=-2.5e+1 # minmax.txt x100
    JMIN=3.1e-5 # minmax.txt x10
    JMAX=3.5e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 203 ]; then
    FILE=w3ra1_2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.15
    YMAX=0.15
    ZMAX=0.15
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.3e+2 # minmax.txt x100
    EMAX=-2.8e+1 # minmax.txt x100
    JMIN=1.0e-4 # minmax.txt x10
    JMAX=2.5e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 204 ]; then
    FILE=w3ra1_4
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.15
    YMAX=0.15
    ZMAX=0.15
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.3e+2 # minmax.txt x100
    EMAX=-2.9e+1 # minmax.txt x100
    JMIN=3.1e-5 # minmax.txt x10
    JMAX=1.4e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 205 ]; then
    FILE=w3ra1_8
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.15
    YMAX=0.15
    ZMAX=0.15
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.4e+2 # minmax.txt x100
    EMAX=-2.9e+1 # minmax.txt x100
    JMIN=3.1e-5 # minmax.txt x10
    JMAX=8.0e-2 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 210 ]; then
    FILE=w5iso
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.22
    YMAX=0.22
    ZMAX=0.22
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+2 # minmax.txt x100
    EMAX=-2.0e+1 # minmax.txt x100
    JMIN=1.0e-4 # minmax.txt x10
    JMAX=7.0e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 211 ]; then
    FILE=w5ra2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.22
    YMAX=0.22
    ZMAX=0.22
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+2 # minmax.txt x100
    EMAX=-1.8e+1 # minmax.txt x100
    JMIN=1.0e-4 # minmax.txt x10
    JMAX=5.0e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 212 ]; then
    FILE=w5ra1
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.22
    YMAX=0.22
    ZMAX=0.22
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+2 # minmax.txt x100
    EMAX=-1.8e+1 # minmax.txt x100
    JMIN=8.0e-5 # minmax.txt x10
    JMAX=4.0e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 213 ]; then
    FILE=w5ra1_2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.22
    YMAX=0.22
    ZMAX=0.22
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+2 # minmax.txt x100
    EMAX=-1.8e+1 # minmax.txt x100
    JMIN=8.0e-5 # minmax.txt x10
    JMAX=3.0e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 214 ]; then
    FILE=w5ra1_4
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.22
    YMAX=0.22
    ZMAX=0.22
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.3e+2 # minmax.txt x100
    EMAX=-1.8e+1 # minmax.txt x100
    JMIN=5.0e-5 # minmax.txt x10
    JMAX=1.5e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 215 ]; then
    FILE=w5ra1_8
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.22
    YMAX=0.22
    ZMAX=0.22
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.3e+2 # minmax.txt x100
    EMAX=-1.8e+1 # minmax.txt x100
    JMIN=1.0e-5 # minmax.txt x10
    JMAX=8.0e-2 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 220 ]; then
    FILE=w7iso
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.34
    YMAX=0.34
    ZMAX=0.34
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.9e+2 # minmax.txt x100
    EMAX=-1.0e+1 # minmax.txt x100
    JMIN=8.0e-5 # minmax.txt x10
    JMAX=9.0e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 221 ]; then
    FILE=w7ra2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.34
    YMAX=0.34
    ZMAX=0.34
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.9e+2 # minmax.txt x100
    EMAX=-1.0e+1 # minmax.txt x100
    JMIN=8.0e-5 # minmax.txt x10
    JMAX=5.0e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 222 ]; then
    FILE=w7ra1
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.34
    YMAX=0.34
    ZMAX=0.34
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.9e+2 # minmax.txt x100
    EMAX=-1.0e+1 # minmax.txt x100
    JMIN=4.0e-5 # minmax.txt x10
    JMAX=4.0e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 223 ]; then
    FILE=w7ra1_2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.34
    YMAX=0.34
    ZMAX=0.34
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.9e+2 # minmax.txt x100
    EMAX=-1.0e+1 # minmax.txt x100
    JMIN=4.0e-5 # minmax.txt x10
    JMAX=3.0e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 224 ]; then
    FILE=w7ra1_4
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.34
    YMAX=0.34
    ZMAX=0.34
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.9e+2 # minmax.txt x100
    EMAX=-1.0e+1 # minmax.txt x100
    JMIN=4.0e-5 # minmax.txt x10
    JMAX=1.5e-1 # minmax.txt x10
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 225 ]; then
    FILE=w7ra1_8
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=0.34
    YMAX=0.34
    ZMAX=0.34
    VMAX=20.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.8e+2 # minmax.txt x100
    EMAX=-1.0e+1 # minmax.txt x100
    JMIN=1.0e-5 # minmax.txt x10
    JMAX=8.0e-2 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 230 ]; then
    FILE=m12iso
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=200.0
    YMAX=200.0
    ZMAX=200.0
    VMAX=300.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.2e+5 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-2 # minmax.txt x10
    JMAX=3.4e+4 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 231 ]; then
    FILE=m12ra2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=200.0
    YMAX=200.0
    ZMAX=200.0
    VMAX=300.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.2e+5 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-2 # minmax.txt x10
    JMAX=2.2e+4 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 232 ]; then
    FILE=m12ra1
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=200.0
    YMAX=200.0
    ZMAX=200.0
    VMAX=300.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.2e+5 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-2 # minmax.txt x10
    JMAX=1.6e+4 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 233 ]; then
    FILE=m12ra1_2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=200.0
    YMAX=200.0
    ZMAX=200.0
    VMAX=300.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.2e+5 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-2 # minmax.txt x10
    JMAX=1.1e+4 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 234 ]; then
    FILE=m12ra1_4
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=200.0
    YMAX=200.0
    ZMAX=200.0
    VMAX=300.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.2e+5 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-2 # minmax.txt x10
    JMAX=2.1e+4 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 235 ]; then
    FILE=m12ra1_8
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=200.0
    YMAX=200.0
    ZMAX=200.0
    VMAX=300.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.2e+5 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=1.0e-2 # minmax.txt x10
    JMAX=2.8e+4 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 240 ]; then
    FILE=m09iso
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=7.0
    YMAX=7.0
    ZMAX=7.0
    VMAX=40.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-4 # minmax.txt x10
    JMAX=3.3e+2 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 241 ]; then
    FILE=m09ra2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=7.0
    YMAX=7.0
    ZMAX=7.0
    VMAX=40.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-5 # minmax.txt x10
    JMAX=2.1e+2 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 242 ]; then
    FILE=m09ra1
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=7.0
    YMAX=7.0
    ZMAX=7.0
    VMAX=40.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-5 # minmax.txt x10
    JMAX=1.5e+2 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 243 ]; then
    FILE=m09ra1_2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=7.0
    YMAX=7.0
    ZMAX=7.0
    VMAX=40.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-5 # minmax.txt x10
    JMAX=1.0e+2 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 244 ]; then
    FILE=m09ra1_4
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=15.0
    YMAX=15.0
    ZMAX=15.0
    VMAX=40.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-5 # minmax.txt x10
    JMAX=2.7e+2 # minmax.txt x10
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 245 ]; then
    FILE=m09ra1_8
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=15.0
    YMAX=15.0
    ZMAX=15.0
    VMAX=40.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-1.6e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-5 # minmax.txt x10
    JMAX=2.4e+2 # minmax.txt x10
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 250 ]; then
    FILE=a00iso
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=5.0
    YMAX=5.0
    ZMAX=5.0
    VMAX=30.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-2.4e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=2.0e-3 # minmax.txt x10
    JMAX=4.9e+2 # minmax.txt x10
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 251 ]; then
    FILE=a00ra2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=5.0
    YMAX=5.0
    ZMAX=5.0
    VMAX=30.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-2.4e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-4 # minmax.txt x10
    JMAX=1.6e+2 # minmax.txt x10
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 252 ]; then
    FILE=a00ra1
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=5.0
    YMAX=5.0
    ZMAX=5.0
    VMAX=30.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-2.4e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-4 # minmax.txt x10
    JMAX=1.1e+2 # minmax.txt x10
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 253 ]; then
    FILE=a00ra1_2
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=5.0
    YMAX=5.0
    ZMAX=5.0
    VMAX=30.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-2.4e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-4 # minmax.txt x10
    JMAX=7.0e+1 # minmax.txt x10
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 254 ]; then
    FILE=a00ra1_4
    FINISH=14000.0
    INTERVAL=100.0
    XMAX=10.0
    YMAX=10.0
    ZMAX=10.0
    VMAX=30.0
    XMIN=-$XMAX
    YMIN=-$YMAX
    ZMIN=-$ZMAX
    VMIN=-$VMAX
    EMIN=-2.4e+3 # minmax.txt x100
    EMAX=-0.0e+1 # minmax.txt x100
    JMIN=3.1e-4 # minmax.txt x10
    JMAX=1.7e+2 # minmax.txt x10
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
# INCREMENT=$END
###############################################################
# set input arguments
OPTION="-file=$FILE -start=$START -end=$END -interval=$INCREMENT -ncrit=$NCRIT -nx=$NX -xmin=$XMIN -xmax=$XMAX -ny=$NY -ymin=$YMIN -ymax=$YMAX -nz=$NZ -zmin=$ZMIN -zmax=$ZMAX -nv=$NV -vmin=$VMIN -vmax=$VMAX -nx3D=$NX3D -ny3D=$NY3D -nz3D=$NZ3D -nE=$NE -Emin=$EMIN -Emax=$EMAX -nJ=$NJ -Jmin=$JMIN -Jmax=$JMAX"
###############################################################


###############################################################
# job execution via NQSV
###############################################################
# set stdout and stderr
STDOUT=log/${FILE}_$PBS_JOBNAME.o${PBS_JOBID}
STDERR=log/${FILE}_$PBS_JOBNAME.e${PBS_JOBID}
###############################################################
# start logging
cd $PBS_O_WORKDIR
TIME=`date`
echo "start: $TIME"
###############################################################
module purge
export MODULEPATH=$MODULEPATH:/work/CSPP/ymiki/opt/module
module load mvapich/2.3.1/intel-cuda10.1
module load gsl lis
module load phdf5
###############################################################
# set number of MPI processes per node
PROCS=`expr $NUM_NODES \* $PROCS_PER_NODE`
###############################################################
# execute the job
echo "mpiexec ${NQSII_MPIOPTS} -genv MV2_NUM_HCAS 4 -np $PROCS sh/wrapper.sh $EXEC log/${FILE}_${PBS_JOBNAME} $PBS_JOBID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION"
mpiexec ${NQSII_MPIOPTS} -genv MV2_NUM_HCAS 4 -np $PROCS sh/wrapper.sh $EXEC log/${FILE}_${PBS_JOBNAME} $PBS_JOBID $PROCS_PER_NODE $PROCS_PER_SOCKET $OPTION
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
