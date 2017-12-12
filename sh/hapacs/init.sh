#!/bin/bash
###############################################################
#PBS -S /bin/bash
#PBS -N MAGI
#PBS -A GALAXY
#PBS -q tcag
#PBS -l select=1:ncpus=16:mpiprocs=1:ompthreads=16
#PBS -l place=scatter
#PBS -l walltime=00:30:00
#PBS -j oe
###############################################################
NTOT=8388608
SAVE=55.0
###############################################################
INI=bin/magi
###############################################################
PROBLEM=10
###############################################################
#
#
###############################################################
# problem specific configurations
###############################################################
# dynamical stability of a King sphere
if [ $PROBLEM -eq 1 ]; then
FILE=king
CONFIG=single/king.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
FILE=hernquist
CONFIG=single/hernquist.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of an NFW sphere with small truncation radius
if [ $PROBLEM -eq 3 ]; then
FILE=nfw
CONFIG=single/nfw.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of an Einasto profile
if [ $PROBLEM -eq 4 ]; then
FILE=einasto
CONFIG=single/einasto.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a Plummer profile
if [ $PROBLEM -eq 5 ]; then
FILE=plummer
CONFIG=single/plummer.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a Burkert profile
if [ $PROBLEM -eq 6 ]; then
FILE=burkert
CONFIG=single/burkert.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a Moore profile
if [ $PROBLEM -eq 7 ]; then
FILE=moore
CONFIG=single/moore.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a Two-power sphere
if [ $PROBLEM -eq 8 ]; then
FILE=twopower
CONFIG=single/twopower.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a King sphere within an Einasto profile
if [ $PROBLEM -eq 10 ]; then
FILE=hb
CONFIG=multi/hb.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere and a King sphere
if [ $PROBLEM -eq 11 ]; then
FILE=hbd
CONFIG=multi/hbd.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
# FINISH=1648.0
# INTERVAL=16.0
fi
###############################################################
# dynamical stability of exponential disks (thick/thin disk) in an NFW sphere and a King sphere
if [ $PROBLEM -eq 12 ]; then
FILE=hbdd
CONFIG=multi/hbdd.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
# FINISH=1400.0
# INTERVAL=8.0
fi
###############################################################
# dynamical stability of thick exponential disk and thin Sersic disk in an Einasto sphere and a King sphere
if [ $PROBLEM -eq 13 ]; then
FILE=ekes
CONFIG=multi/ekes.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=94.0
INTERVAL=2.0
# FINISH=1400.0
# INTERVAL=8.0
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007)
if [ $PROBLEM -eq 20 ]; then
FILE=m31
CONFIG=galaxy/f07.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=75.0
INTERVAL=25.0
# FINISH=5175.0
# INTERVAL=25.0
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007) in unit of GalactICS system
if [ $PROBLEM -eq 21 ]; then
FILE=m31_gics
CONFIG=galaxy/f07_gics.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=5175.0
INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model
if [ $PROBLEM -eq 22 ]; then
FILE=galaxy
CONFIG=galaxy/galaxy.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=75.0
INTERVAL=25.0
# FINISH=3175.0
# INTERVAL=25.0
fi
###############################################################
# dynamical stability of MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008)
if [ $PROBLEM -eq 23 ]; then
FILE=mw
CONFIG=galaxy/mw.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=75.0
INTERVAL=25.0
# FINISH=3175.0
# INTERVAL=25.0
fi
###############################################################
# dynamical stability of M31 model (Sofue 2015; Gilbert et al. 2012)
if [ $PROBLEM -eq 24 ]; then
FILE=s15
CONFIG=galaxy/s15.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=75.0
INTERVAL=25.0
# FINISH=3175.0
# INTERVAL=25.0
fi
###############################################################
# time evolution of MW/A defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 30 ]; then
FILE=kd95a
CONFIG=gics/kd95/MWa.cfg
EPS=1.5625e-2
ETA=0.5
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MW/B defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 31 ]; then
FILE=kd95b
CONFIG=gics/kd95/MWb.cfg
EPS=1.5625e-2
ETA=0.5
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MW/C defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 32 ]; then
FILE=kd95c
CONFIG=gics/kd95/MWc.cfg
EPS=1.5625e-2
ETA=0.5
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MW/D defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 33 ]; then
FILE=kd95d
CONFIG=gics/kd95/MWd.cfg
EPS=1.5625e-2
ETA=0.5
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of M31/A defined in Widrow et al. (2003)
if [ $PROBLEM -eq 34 ]; then
FILE=w03a
CONFIG=gics/w03/M31a.cfg
EPS=1.5625e-2
ETA=0.5
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of M31/D defined in Widrow et al. (2003)
if [ $PROBLEM -eq 35 ]; then
FILE=w03d
CONFIG=gics/w03/M31d.cfg
EPS=1.5625e-2
ETA=0.5
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MWa defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 36 ]; then
FILE=mwa
CONFIG=gics/wd05/MWa.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=124.0
INTERVAL=4.0
# FINISH=1016.0
# INTERVAL=8.0
fi
###############################################################
# time evolution of MWb defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 37 ]; then
FILE=mwb
CONFIG=gics/wd05/MWb.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=124.0
INTERVAL=4.0
# FINISH=1016.0
# INTERVAL=8.0
fi
###############################################################
# time evolution of MW like galaxy (based on Widrow et al. 2008; Bedorf et al. 2014)
if [ $PROBLEM -eq 38 ]; then
FILE=bonsai
CONFIG=gics/bonsai/b14.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=124.0
INTERVAL=4.0
# FINISH=1176.0
# INTERVAL=8.0
fi
###############################################################
# dynamical stability of a Plummer profile in table form
if [ $PROBLEM -eq 40 ]; then
FILE=tplummer
CONFIG=table/tplummer.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=1575.0
INTERVAL=25.0
fi
###############################################################
# dynamical stability of a Two-power slope profile in table form
if [ $PROBLEM -eq 41 ]; then
FILE=dblpower
CONFIG=table/dblpower.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=1575.0
INTERVAL=25.0
fi
###############################################################
# dynamical stability of a de Vaucouleurs sphere in table form
if [ $PROBLEM -eq 42 ]; then
FILE=bulgetbl
CONFIG=table/bulgetbl.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=1575.0
INTERVAL=25.0
fi
###############################################################
# dynamical stability of a spherical Sersic profile
if [ $PROBLEM -eq 50 ]; then
FILE=deVaucouleurs
CONFIG=abel/deVaucouleurs.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=1575.0
INTERVAL=25.0
fi
###############################################################
# dynamical stability of a projected Two-power model
if [ $PROBLEM -eq 51 ]; then
FILE=prjTwoPow
CONFIG=abel/prjTwoPow.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=1575.0
INTERVAL=25.0
fi
###############################################################
#
#
###############################################################
# generate initial condition
###############################################################
NODES=1
PROCS_PER_NODE=1
EXEC=$INI
###############################################################
PROCS=`expr $NODES \* $PROCS_PER_NODE`
###############################################################
. /opt/Modules/default/init/bash
export MODULEPATH=$MODULEPATH:/work/GALAXY/ymiki/opt/Modules
module load intel/16.0.4 mvapich2/2.2_intel_cuda-8.0.44 cuda/8.0.44 cuda-samples/8.0.44
module load gsl hdf5
###############################################################
cd $PBS_O_WORKDIR
STDOUT=log/$PBS_JOBNAME.$PBS_JOBID.out
STDERR=log/$PBS_JOBNAME.$PBS_JOBID.err
echo "used hosts:" >> $STDOUT
cat $PBS_NODEFILE >> $STDOUT
TIME=`date`
echo "start: $TIME" >> $STDOUT
###############################################################
# # for IntelMPI
# mpiexec -np $PROCS -hostfile $PBS_NODEFILE -perhost $PROCS_PER_NODE -l -exitinfo $EXEC -file=$FILE -config=$CONFIG -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>$STDOUT 2>>$STDERR
###############################################################
# for mvapich2
mpiexec -np $PROCS -hostfile $PBS_NODEFILE -l -exitinfo $EXEC -file=$FILE -config=$CONFIG -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>$STDOUT 2>>$STDERR
###############################################################
TIME=`date`
echo "finish: $TIME" >> $STDOUT
###############################################################
