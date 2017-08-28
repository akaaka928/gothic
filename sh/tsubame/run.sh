#!/bin/sh
#$ -cwd
#$ -l s_gpu=1
#$ -l h_rt=1:00:00
##$ -N gothic
##$ -hold_jid magi
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/gothic
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=20
    # PROBLEM=26
    # PROBLEM=28
    # PROBLEM=80
    # PROBLEM=81
fi
###############################################################
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
###############################################################
# value of accuracy controling parameter: GADGET MAC by Springel (2005)
if [ -z "$ABSERR" ]; then
    # ABSERR=1.250000000e-1
    # ABSERR=6.250000000e-2
    # ABSERR=3.125000000e-2
    # ABSERR=1.562500000e-2
    # ABSERR=7.812500000e-3
    ABSERR=3.906250000e-3
    # ABSERR=1.953125000e-3
    # ABSERR=9.765625000e-4
    # ABSERR=4.882812500e-4
    # ABSERR=2.441406250e-4
    # ABSERR=1.220703125e-4
fi
###############################################################
# value of accuracy controling parameter: opening criterion by Barnes & Hut (1986)
if [ -z "$THETA" ]; then
    # THETA=0.9
    # THETA=0.8
    # THETA=0.7
    # THETA=0.6
    # THETA=0.5
    THETA=0.4
    # THETA=0.3
fi
###############################################################
# value of accuracy controling parameter: multipole moment MAC by Warren & Salmon (1993)
if [ -z "$ACCERR" ]; then
    # ACCERR=1.250000e-1
    # ACCERR=6.250000e-2
    # ACCERR=3.125000e-2
    ACCERR=1.562500e-2
    # ACCERR=7.812500e-3
    # ACCERR=3.906250e-3
    # ACCERR=1.953125e-3
    # ACCERR=9.765625e-4
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# cold collapse of a uniform sphere
if [ $PROBLEM -eq 0 ]; then
    FILE=ccuni
fi
###############################################################
# dynamical stability of a King sphere
if [ $PROBLEM -eq 1 ]; then
    FILE=king
fi
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
    FILE=hernquist
fi
###############################################################
# dynamical stability of an NFW sphere small truncation radius
if [ $PROBLEM -eq 3 ]; then
    FILE=nfw
fi
###############################################################
# dynamical stability of an Einasto profile
if [ $PROBLEM -eq 4 ]; then
    FILE=einasto
fi
###############################################################
# dynamical stability of a Plummer profile
if [ $PROBLEM -eq 5 ]; then
    FILE=plummer
fi
###############################################################
# dynamical stability of a Burkert profile
if [ $PROBLEM -eq 6 ]; then
    FILE=burkert
fi
###############################################################
# dynamical stability of a Moore profile
if [ $PROBLEM -eq 7 ]; then
    FILE=moore
fi
###############################################################
# dynamical stability of a Two-power sphere
if [ $PROBLEM -eq 8 ]; then
    FILE=twopower
fi
###############################################################
# dynamical stability of a King sphere within an Einasto profile
if [ $PROBLEM -eq 10 ]; then
    FILE=hb
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere and a King sphere
if [ $PROBLEM -eq 11 ]; then
    FILE=hbd
fi
###############################################################
# dynamical stability of exponential disks (thick/thin disk) in an NFW sphere and a King sphere
if [ $PROBLEM -eq 12 ]; then
    FILE=hbdd
fi
###############################################################
# dynamical stability of thick exponential disk and thin Sersic disk in an Einasto sphere and a King sphere
if [ $PROBLEM -eq 13 ]; then
    FILE=ekes
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007)
if [ $PROBLEM -eq 20 ]; then
    FILE=m31
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007) in unit of GalactICS system
if [ $PROBLEM -eq 21 ]; then
    FILE=m31_gics
fi
###############################################################
# dynamical stability of multi components galaxy model
if [ $PROBLEM -eq 22 ]; then
    FILE=galaxy
fi
###############################################################
# dynamical stability of MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008)
if [ $PROBLEM -eq 23 ]; then
    FILE=mw
fi
###############################################################
# dynamical stability of M31 model (Sofue 2015; Gilbert et al. 2012)
if [ $PROBLEM -eq 24 ]; then
    FILE=s15
fi
###############################################################
# dynamical stability of multi components galaxy model with fixed number of particles
if [ $PROBLEM -eq 25 ]; then
    FILE=compare
fi
###############################################################
# dynamical stability of multi components galaxy model (only spherical components)
if [ $PROBLEM -eq 26 ]; then
    FILE=spherical
fi
###############################################################
# dynamical stability of an M31 model (NFW halo, de Vaucouleurs bulge, and exponential disk)
if [ $PROBLEM -eq 27 ]; then
    FILE=m31_mod
fi
###############################################################
# dynamical stability of multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk)
if [ $PROBLEM -eq 28 ]; then
    FILE=ltg
fi
###############################################################
# dynamical stability of multi components galaxy model (Vasiliev & Athanassoula 2015) with fixed number of particles
if [ $PROBLEM -eq 29 ]; then
    FILE=va15
fi
###############################################################
# time evolution of MW/A defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 30 ]; then
    FILE=kd95a
fi
###############################################################
# time evolution of MW/B defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 31 ]; then
    FILE=kd95b
fi
###############################################################
# time evolution of MW/C defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 32 ]; then
    FILE=kd95c
fi
###############################################################
# time evolution of MW/D defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 33 ]; then
    FILE=kd95d
fi
###############################################################
# time evolution of M31/A defined in Widrow et al. (2003)
if [ $PROBLEM -eq 34 ]; then
    FILE=w03a
fi
###############################################################
# time evolution of M31/D defined in Widrow et al. (2003)
if [ $PROBLEM -eq 35 ]; then
    FILE=w03d
fi
###############################################################
# time evolution of MWa defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 36 ]; then
    FILE=mwa
fi
###############################################################
# time evolution of MWb defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 37 ]; then
    FILE=mwb
fi
###############################################################
# time evolution of MW like galaxy (based on Widrow et al. 2008, Bedorf et al. 2014)
if [ $PROBLEM -eq 38 ]; then
    FILE=bonsai
fi
###############################################################
# dynamical stability of a Plummer profile in table form
if [ $PROBLEM -eq 40 ]; then
    FILE=tplummer
fi
###############################################################
# dynamical stability of a Two-power slope profile in table form
if [ $PROBLEM -eq 41 ]; then
    FILE=dblpower
fi
###############################################################
# dynamical stability of a de Vaucouleurs sphere in table form
if [ $PROBLEM -eq 42 ]; then
    FILE=bulgetbl
fi
###############################################################
# dynamical stability of a spherical Sersic profile
if [ $PROBLEM -eq 50 ]; then
    FILE=deVaucouleurs
fi
###############################################################
# dynamical stability of a projected Two-power model
if [ $PROBLEM -eq 51 ]; then
    FILE=prjTwoPow
fi
###############################################################
# dynamical stability of multi components galaxy model (based on Cole & Binney 2017)
if [ $PROBLEM -eq 80 ]; then
    FILE=cb17
fi
###############################################################
# dynamical stability of multi components galaxy model (based on Cole & Binney 2017)
if [ $PROBLEM -eq 81 ]; then
    FILE=cb17_core
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -accErr=$ACCERR -theta=$THETA -file=$FILE -Nx=$NX -Ny=$NY -Nz=$NZ -jobID=$SLURM_JOB_ID"
###############################################################


###############################################################
# job execution via UNIVA Grid Engine
###############################################################
# set stdout and stderr
STDOUT=log/$REQUEST.o$JOB_ID
STDERR=log/$REQUEST.e$JOB_ID
###############################################################
# load modules
. /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:$HOME/opt/Modules
module load intel intel-mpi phdf5
module load gsl
module load cuda cub
module list 1>>$STDOUT 2>>$STDERR
###############################################################
cat $PE_HOSTFILE 1>>$STDOUT 2>>$STDERR
TIME=`date`
echo "start: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################
# execute the job
if [ `which numactl` ]; then
    # run with numactl
    echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    # run without numactl
    echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################
