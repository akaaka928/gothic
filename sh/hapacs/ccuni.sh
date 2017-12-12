#!/bin/bash
###############################################################
#PBS -S /bin/bash
#PBS -N ccuni
#PBS -A GALAXY
#PBS -q tcag
#PBS -l select=1:ncpus=16:mpiprocs=1:ompthreads=16
#PBS -l place=scatter
#PBS -l walltime=00:30:00
#PBS -j oe
###############################################################
NTOT=131072
SAVE=140.0
###############################################################
INI=bin/uniformsphere
###############################################################
#
#
###############################################################
# problem specific configurations
###############################################################
FILE=ccuni
UNIT=-1
MTOT=1.0
VIRIAL=0.2
RAD=1.0
EPS=1.5625e-2
ETA=0.125
FINISH=11.75
INTERVAL=0.25
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
module load intel/16.0.4 mvapich2/2.2_intel_cuda-7.5.18
module load cuda/7.5.18 cuda-samples/7.5.18
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
# mpiexec -np $PROCS -hostfile $PBS_NODEFILE -perhost $PROCS_PER_NODE -l -exitinfo $EXEC -file=$FILE -unit=$UNIT -Ntot=$NTOT -Mtot=$MTOT -virial=$VIRIAL -rad=$RAD -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>$STDOUT 2>>$STDERR
###############################################################
# for mvapich2
mpiexec -np $PROCS -hostfile $PBS_NODEFILE -l -exitinfo $EXEC -file=$FILE -unit=$UNIT -Ntot=$NTOT -Mtot=$MTOT -virial=$VIRIAL -rad=$RAD -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>$STDOUT 2>>$STDERR
###############################################################
TIME=`date`
echo "finish: $TIME" >> $STDOUT
###############################################################
