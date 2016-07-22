#!/bin/bash
###############################################################
#PBS -S /bin/bash
#PBS -N GOTHIC
#PBS -A GALAXY
#PBS -q tcag
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1
#PBS -l place=scatter
#PBS -l walltime=24:00:00
#PBS -j oe
###############################################################
NODES=1
PROCS_PER_NODE=1
PROCS=`expr $NODES \* $PROCS_PER_NODE`
###############################################################
EXEC=bin/gothic
###############################################################
PROBLEM=10
###############################################################
#
#
###############################################################
# accuracy contoling parameters
###############################################################
# opening criterion by Barnes & Hut (1986)
# THETA=0.9
# THETA=0.8
THETA=0.7
# THETA=0.6
# THETA=0.5
# THETA=0.4
# THETA=0.3
###############################################################
# multipole moment MAC by Warren & Salmon (1993)
# ACCERR=1.250000e-1
# ACCERR=6.250000e-2
# ACCERR=3.125000e-2
# ACCERR=1.562500e-2
ACCERR=7.812500e-3
# ACCERR=3.906250e-3
# ACCERR=1.953125e-3
# ACCERR=9.765625e-4
###############################################################
# GADGET MAC by Springel (2005)
# ABSERR=1.250000000e-1
# ABSERR=6.250000000e-2
# ABSERR=3.125000000e-2
# ABSERR=1.562500000e-2
ABSERR=7.812500000e-3
# ABSERR=3.906250000e-3
# ABSERR=1.953125000e-3
# ABSERR=9.765625000e-4
# ABSERR=4.882812500e-4
# ABSERR=2.441406250e-4
# ABSERR=1.220703125e-4
###############################################################
#
#
###############################################################
# specify the target problem
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
#
#
###############################################################
# execute numerical simulation
###############################################################
. /opt/Modules/default/init/bash
module load intel/15.0.5 intelmpi/5.1.1
module load cuda/7.5.18 cuda/samples_7.5.18
# module load ddt/5.1 reports/5.1
###############################################################
cd $PBS_O_WORKDIR
STDOUT=log/$PBS_JOBNAME.$PBS_JOBID.out
STDERR=log/$PBS_JOBNAME.$PBS_JOBID.err
echo "used hosts:" >> $STDOUT
cat $PBS_NODEFILE >> $STDOUT
TIME=`date`
echo "start: $TIME" >> $STDOUT
###############################################################
# export ALLINEA_NO_TIMEOUT
# perf-report mpirun -np $PROCS -hostfile $PBS_NODEFILE -perhost $PROCS_PER_NODE $EXEC -absErr=$ABSERR -accErr=$ACCERR -theta=$THETA -file=$FILE -jobID=$PBS_JOBID 1>>$STDOUT 2>>$STDERR
mpirun -np $PROCS -hostfile $PBS_NODEFILE -perhost $PROCS_PER_NODE $EXEC -absErr=$ABSERR -accErr=$ACCERR -theta=$THETA -file=$FILE -jobID=$PBS_JOBID 1>>$STDOUT 2>>$STDERR
###############################################################
TIME=`date`
echo "finish: $TIME" >> $STDOUT
###############################################################
