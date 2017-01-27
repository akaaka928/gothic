#!/bin/bash
###############################################################
#PBS -S /bin/bash
#PBS -N PLplot
#PBS -A GALAXY
#PBS -q tcag
#PBS -l select=1:ncpus=16:mpiprocs=16:ompthreads=1
#PBS -l place=scatter
#PBS -l walltime=00:30:00
#PBS -j oe
###############################################################
# set # of N-body particles per bin to estimate density
# NCRIT=128
# NCRIT=256
# NCRIT=512
# NCRIT=1024
NCRIT=2048
###############################################################
EXTENSION=svg
###############################################################
PLTENE=bin/plot.energy
PLTMAP=bin/plot.distribution
PLTCDF=bin/plot.cdf
###############################################################
# PROBLEM=3
# TARGET=0
###############################################################
#
#
###############################################################
# specific configurations
###############################################################
# dynamical stability of a King sphere
if [ $PROBLEM -eq 1 ]; then
FILE=king
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
FILE=hernquist
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of an NFW sphere with small truncation radius
if [ $PROBLEM -eq 3 ]; then
FILE=nfw
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of an Einasto profile
if [ $PROBLEM -eq 4 ]; then
FILE=einasto
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a Plummer profile
if [ $PROBLEM -eq 5 ]; then
FILE=plummer
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a Burkert profile
if [ $PROBLEM -eq 6 ]; then
FILE=burkert
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a Moore profile
if [ $PROBLEM -eq 7 ]; then
FILE=moore
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of a King sphere within an Einasto profile
if [ $PROBLEM -eq 10 ]; then
FILE=hb
FINISH=94.0
INTERVAL=2.0
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere and a King sphere
if [ $PROBLEM -eq 11 ]; then
FILE=hbd
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# dynamical stability of exponential disks (thick/thin disk) in an NFW sphere and a King sphere
if [ $PROBLEM -eq 12 ]; then
FILE=hbdd
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# dynamical stability of thick exponential disk and thin Sersic disk in an Einasto sphere and a King sphere
if [ $PROBLEM -eq 13 ]; then
FILE=ekes
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007)
if [ $PROBLEM -eq 20 ]; then
FILE=m31
FINISH=5175.0
INTERVAL=25.0
fi
###############################################################
# time evolution of MW/A defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 30 ]; then
FILE=kd95a
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MW/B defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 31 ]; then
FILE=kd95b
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MW/C defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 32 ]; then
FILE=kd95c
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MW/D defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 33 ]; then
FILE=kd95d
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of M31/A defined in Widrow et al. (2003)
if [ $PROBLEM -eq 34 ]; then
FILE=w03a
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of M31/D defined in Widrow et al. (2003)
if [ $PROBLEM -eq 35 ]; then
FILE=w03d
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MWa defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 36 ]; then
FILE=mwa
FINISH=1016.0
INTERVAL=8.0
fi
###############################################################
# time evolution of MWb defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 37 ]; then
FILE=mwb
FINISH=1016.0
INTERVAL=8.0
fi
###############################################################
# time evolution of MW like galaxy (based on Widrow et al. 2008, Bedorf et al. 2014)
if [ $PROBLEM -eq 38 ]; then
FILE=bonsai
FINISH=1176.0
INTERVAL=8.0
fi
###############################################################
#
#
###############################################################
# plot figures
###############################################################
START=0
END=`echo "scale=0; $FINISH / $INTERVAL" | bc`
INCREMENT=1
###############################################################
NODES=1
PROCS_PER_NODE=16
###############################################################
PROCS=`expr $NODES \* $PROCS_PER_NODE`
###############################################################
. /opt/Modules/default/init/bash
export MODULEPATH=$MODULEPATH:/work/GALAXY/ymiki/opt/Modules
module load intel/15.0.5 intelmpi/5.1.1
module load cuda/7.5.18 cuda/samples_7.5.18
module load gsl hdf5 plplot
###############################################################
cd $PBS_O_WORKDIR
STDOUT=log/$PBS_JOBNAME.$PBS_JOBID.out
STDERR=log/$PBS_JOBNAME.$PBS_JOBID.err
echo "used hosts:" >> $STDOUT
cat $PBS_NODEFILE >> $STDOUT
TIME=`date`
echo "start: $TIME" >> $STDOUT
###############################################################
if [ $TARGET -eq 0 ]; then
mpirun -np $PROCS -hostfile $PBS_NODEFILE -perhost $PROCS_PER_NODE $PLTENE -file=$FILE -start=$START -end=$END -interval=$INCREMENT -problem=$PROBLEM               -dev ${EXTENSION}cairo 1>>$STDOUT 2>>$STDERR
fi
if [ $TARGET -eq 1 ]; then
mpirun -np $PROCS -hostfile $PBS_NODEFILE -perhost $PROCS_PER_NODE $PLTMAP -file=$FILE -start=$START -end=$END -interval=$INCREMENT -problem=$PROBLEM -ncrit=$NCRIT -dev ${EXTENSION}cairo 1>>$STDOUT 2>>$STDERR
fi
if [ $TARGET -eq 2 ]; then
mpirun -np $PROCS -hostfile $PBS_NODEFILE -perhost $PROCS_PER_NODE $PLTCDF -file=$FILE                                                                              -dev ${EXTENSION}cairo 1>>$STDOUT 2>>$STDERR
fi
###############################################################
TIME=`date`
echo "finish: $TIME" >> $STDOUT
###############################################################
