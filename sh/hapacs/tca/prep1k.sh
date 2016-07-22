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
SAVE=1380.0
###############################################################
INI=bin/magi
###############################################################
PROBLEM=1
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
FINISH=4.0
INTERVAL=8.0
fi
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
FILE=hernquist
CONFIG=single/hernquist.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=4.0
INTERVAL=8.0
fi
###############################################################
# dynamical stability of an NFW sphere with small truncation radius
if [ $PROBLEM -eq 3 ]; then
FILE=nfw
CONFIG=single/nfw.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=4.0
INTERVAL=8.0
fi
###############################################################
# dynamical stability of an Einasto profile
if [ $PROBLEM -eq 4 ]; then
FILE=einasto
CONFIG=single/einasto.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=4.0
INTERVAL=8.0
fi
###############################################################
# dynamical stability of a Plummer profile
if [ $PROBLEM -eq 5 ]; then
FILE=plummer
CONFIG=single/plummer.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=4.0
INTERVAL=8.0
fi
###############################################################
# dynamical stability of a Burkert profile
if [ $PROBLEM -eq 6 ]; then
FILE=burkert
CONFIG=single/burkert.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=8.0
INTERVAL=16.0
fi
###############################################################
# dynamical stability of a Moore profile
if [ $PROBLEM -eq 7 ]; then
FILE=moore
CONFIG=single/moore.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=4.0
INTERVAL=8.0
fi
###############################################################
# dynamical stability of a Two-power sphere
if [ $PROBLEM -eq 8 ]; then
FILE=twopower
CONFIG=single/twopower.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=2.0
INTERVAL=4.0
fi
###############################################################
# dynamical stability of a King sphere within an Einasto profile
if [ $PROBLEM -eq 10 ]; then
FILE=hb
CONFIG=multi/hb.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=4.0
INTERVAL=8.0
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere and a King sphere
if [ $PROBLEM -eq 11 ]; then
FILE=hbd
CONFIG=multi/hbd.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=8.0
INTERVAL=16.0
fi
###############################################################
# dynamical stability of exponential disks (thick/thin disk) in an NFW sphere and a King sphere
if [ $PROBLEM -eq 12 ]; then
FILE=hbdd
CONFIG=multi/hbdd.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=8.0
INTERVAL=16.0
fi
###############################################################
# dynamical stability of thick exponential disk and thin Sersic disk in an Einasto sphere and a King sphere
if [ $PROBLEM -eq 13 ]; then
FILE=ekes
CONFIG=multi/ekes.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=4.0
INTERVAL=8.0
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007)
if [ $PROBLEM -eq 20 ]; then
FILE=m31
CONFIG=galaxy/f07.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=6.25
INTERVAL=25.0
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007) in unit of GalactICS system
if [ $PROBLEM -eq 21 ]; then
FILE=m31_gics
CONFIG=galaxy/f07_gics.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=4.88896113
INTERVAL=2.5
fi
###############################################################
# dynamical stability of multi components galaxy model
if [ $PROBLEM -eq 22 ]; then
FILE=galaxy
CONFIG=galaxy/galaxy.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=6.25
INTERVAL=25.0
fi
###############################################################
# dynamical stability of MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008)
if [ $PROBLEM -eq 23 ]; then
FILE=mw
CONFIG=galaxy/mw.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=6.25
INTERVAL=25.0
fi
###############################################################
# dynamical stability of M31 model (Sofue 2015; Gilbert et al. 2012)
if [ $PROBLEM -eq 24 ]; then
FILE=s15
CONFIG=galaxy/s15.cfg
EPS=1.5625e-2
ETA=0.5
FINISH=6.25
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
module load intel/15.0.5 intelmpi/5.1.1
module load cuda/7.5.18 cuda/samples_7.5.18
###############################################################
cd $PBS_O_WORKDIR
STDOUT=log/$PBS_JOBNAME.$PBS_JOBID.out
STDERR=log/$PBS_JOBNAME.$PBS_JOBID.err
echo "used hosts:" >> $STDOUT
cat $PBS_NODEFILE >> $STDOUT
TIME=`date`
echo "start: $TIME" >> $STDOUT
###############################################################
mpirun -np $PROCS -hostfile $PBS_NODEFILE -perhost $PROCS_PER_NODE $EXEC -file=$FILE -config=$CONFIG -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>$STDOUT 2>>$STDERR
###############################################################
TIME=`date`
echo "finish: $TIME" >> $STDOUT
###############################################################
