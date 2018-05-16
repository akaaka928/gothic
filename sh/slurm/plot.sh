#!/bin/bash
###############################################################
#SBATCH -J plplot             # name of job
#SBATCH -t 00:30:00           # upper limit of elapsed time
#SBATCH -p normal             # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=16           # number of total MPI processes, set to SLURM_NTASKS
##SBATCH --ntasks-per-socket=8 # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET
#SBATCH --ntasks-per-socket=16 # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET
#SBATCH --get-user-env        # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
PLTENE=bin/energy
PLTMAP=bin/distribution
PLTCDF=bin/cdf
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=2
    # PROBLEM=27
    # PROBLEM=111
    # PROBLEM=112
    # PROBLEM=124
fi
###############################################################
# file extension of figures
if [ -z "$EXTENSION" ]; then
    EXTENSION=png
    # EXTENSION=svg
fi
###############################################################
# set number of N-body particles per bin to estimate density
if [ -z "$NCRIT" ]; then
    # NCRIT=8
    # NCRIT=32
    # NCRIT=64
    # NCRIT=128
    # NCRIT=256
    # NCRIT=512
    # NCRIT=1024
    NCRIT=2048
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# cold collapse of a uniform sphere
if [ $PROBLEM -eq 0 ]; then
    FILE=ccuni
    FINISH=11.75
    INTERVAL=0.25
fi
###############################################################
# dynamical stability of a King sphere
if [ $PROBLEM -eq 1 ]; then
    FILE=king
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
    FILE=hernquist
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of an NFW sphere with small truncation radius
if [ $PROBLEM -eq 3 ]; then
    FILE=nfw
    FINISH=23.5
    INTERVAL=0.5
    # FINISH=0.1
    # INTERVAL=0.1
fi
###############################################################
# dynamical stability of an Einasto profile
if [ $PROBLEM -eq 4 ]; then
    FILE=einasto
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Plummer profile
if [ $PROBLEM -eq 5 ]; then
    FILE=plummer
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Burkert profile
if [ $PROBLEM -eq 6 ]; then
    FILE=burkert
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Moore profile
if [ $PROBLEM -eq 7 ]; then
    FILE=moore
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Two-power sphere
if [ $PROBLEM -eq 8 ]; then
    FILE=twopower
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a King sphere within an Einasto profile
if [ $PROBLEM -eq 10 ]; then
    FILE=hb
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere and a King sphere
if [ $PROBLEM -eq 11 ]; then
    FILE=hbd
    FINISH=47.5
    INTERVAL=0.5
    # FINISH=1648.0
    # INTERVAL=16.0
fi
###############################################################
# dynamical stability of exponential disks (thick/thin disk) in an NFW sphere and a King sphere
if [ $PROBLEM -eq 12 ]; then
    FILE=hbdd
    FINISH=47.5
    INTERVAL=0.5
    # FINISH=1400.0
    # INTERVAL=8.0
fi
###############################################################
# dynamical stability of thick exponential disk and thin Sersic disk in an Einasto sphere and a King sphere
if [ $PROBLEM -eq 13 ]; then
    FILE=ekes
    FINISH=47.5
    INTERVAL=0.5
    # FINISH=1400.0
    # INTERVAL=8.0
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007)
if [ $PROBLEM -eq 20 ]; then
    FILE=m31
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=1175.0
    INTERVAL=25.0
    # FINISH=5175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model
if [ $PROBLEM -eq 22 ]; then
    FILE=galaxy
    FINISH=75.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008)
if [ $PROBLEM -eq 23 ]; then
    FILE=mw
    FINISH=75.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of M31 model (Sofue 2015; Gilbert et al. 2012)
if [ $PROBLEM -eq 24 ]; then
    FILE=s15
    FINISH=75.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model with fixed number of particles
if [ $PROBLEM -eq 25 ]; then
    FILE=compare
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=3175.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (only spherical components)
if [ $PROBLEM -eq 26 ]; then
    FILE=etg
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=1175.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of an M31 model (NFW halo, de Vaucouleurs bulge, and exponential disk)
if [ $PROBLEM -eq 27 ]; then
    FILE=m31
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=1175.0
    # FINISH=3175.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk)
if [ $PROBLEM -eq 28 ]; then
    FILE=ltg
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=1175.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (Vasiliev & Athanassoula 2015) with fixed number of particles
if [ $PROBLEM -eq 29 ]; then
    FILE=va15
    FINISH=100.0
    INTERVAL=4.0
fi
###############################################################
# time evolution of MW/A defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 30 ]; then
    FILE=kd95a
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of MW/B defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 31 ]; then
    FILE=kd95b
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of MW/C defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 32 ]; then
    FILE=kd95c
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of MW/D defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 33 ]; then
    FILE=kd95d
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of M31/A defined in Widrow et al. (2003)
if [ $PROBLEM -eq 34 ]; then
    FILE=w03a
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of M31/D defined in Widrow et al. (2003)
if [ $PROBLEM -eq 35 ]; then
    FILE=w03d
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of MWa defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 36 ]; then
    FILE=mwa
    FINISH=124.0
    INTERVAL=4.0
    # FINISH=1016.0
    # INTERVAL=8.0
fi
###############################################################
# time evolution of MWb defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 37 ]; then
    FILE=mwb
    FINISH=124.0
    INTERVAL=4.0
    # FINISH=1016.0
    # INTERVAL=8.0
fi
###############################################################
# time evolution of MW like galaxy (based on Widrow et al. 2008, Bedorf et al. 2014)
if [ $PROBLEM -eq 38 ]; then
    FILE=bonsai
    FINISH=124.0
    INTERVAL=4.0
    # FINISH=1176.0
    # INTERVAL=8.0
fi
###############################################################
# dynamical stability of a Plummer profile in table form
if [ $PROBLEM -eq 40 ]; then
    FILE=tplummer
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a Two-power slope profile in table form
if [ $PROBLEM -eq 41 ]; then
    FILE=dblpower
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a de Vaucouleurs sphere in table form
if [ $PROBLEM -eq 42 ]; then
    FILE=bulgetbl
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a spherical Sersic profile
if [ $PROBLEM -eq 50 ]; then
    FILE=deVaucouleurs
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a projected Two-power model
if [ $PROBLEM -eq 51 ]; then
    FILE=prjTwoPow
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (based on Cole & Binney 2017)
if [ $PROBLEM -eq 80 ]; then
    FILE=cb17
    FINISH=1175.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (based on Cole & Binney 2017)
if [ $PROBLEM -eq 81 ]; then
    FILE=cb17_core
    FINISH=1175.0
    INTERVAL=25.0
    PROBLEM=80
fi
###############################################################
# reproduction of Miki et al. (2016) in the disk coordinate system
if [ $PROBLEM -eq 111 ]; then
    FILE=m16king
    FINISH=1246.875
    INTERVAL=3.125
fi
###############################################################
# reproduction of Kirihara et al. (2017) in the disk coordinate system
if [ $PROBLEM -eq 112 ]; then
    FILE=k17disk
    FINISH=1246.875
    INTERVAL=3.125
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 121 ]; then
    FILE=halocore1_run
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 122 ]; then
    FILE=halocore2_run
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 123 ]; then
    FILE=halocore3_run
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 124 ]; then
    FILE=halocusp_run
    # FINISH=14000.0
    FINISH=4400.0
    INTERVAL=100.0
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
OPTENE="-file=$FILE -problem=$PROBLEM -start=$START -end=$END -interval=$INCREMENT -dev ${EXTENSION}cairo"
OPTCDF="-file=$FILE -dev ${EXTENSION}cairo"
OPTMAP="-file=$FILE -problem=$PROBLEM -ncrit=$NCRIT -start=$START -end=$END -interval=$INCREMENT -dev ${EXTENSION}cairo"
###############################################################


###############################################################
# job execution via SLURM
###############################################################
# set number of MPI processes per node
PROCS_PER_NODE=`expr $SLURM_NTASKS / $SLURM_JOB_NUM_NODES`
###############################################################
# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_NUM_NODES nodes"
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs per node"
TIME=`date`
echo "start: $TIME"
###############################################################
# execute the job
echo "mpiexec -n $SLURM_NTASKS sh/wrapper.sh $PLTENE log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTENE"
mpiexec -n $SLURM_NTASKS sh/wrapper.sh $PLTENE log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTENE

echo "mpiexec -n $SLURM_NTASKS sh/wrapper.sh $PLTMAP log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTMAP"
mpiexec -n $SLURM_NTASKS sh/wrapper.sh $PLTMAP log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTMAP

if [ -e dat/$FILE.direct000.dat ]; then
    echo "mpiexec -n $SLURM_NTASKS sh/wrapper.sh $PLTCDF log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTCDF"
    mpiexec -n $SLURM_NTASKS sh/wrapper.sh $PLTCDF log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTCDF
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
