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
# THETA=0.6
# ACCERR=7.812500e-3
# ABSERR=7.812500e-3
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
echo "0" > doc/$FILE.cfg.dat
cp ../../dat/$FILE.tmp0.h5 dat/
###############################################################
mpirun -np $PROCS -hostfile $PBS_NODEFILE -perhost $PROCS_PER_NODE $EXEC -absErr=$VALUE -accErr=$VALUE -theta=$VALUE -file=$FILE -jobID=$PBS_JOBID 1>>$STDOUT 2>>$STDERR
###############################################################
TIME=`date`
echo "finish: $TIME" >> $STDOUT
###############################################################
