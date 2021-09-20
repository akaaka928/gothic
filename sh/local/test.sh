#!/bin/sh
###############################################################
if [ $# -lt 2 ]; then
    echo "$# inputs are detected while at least 2 inputs are required to specify <test problem> and <extension>" 1>&2
    exit 1
fi
PROBLEM=$1
JOB_ID=$$
if [ $# -ge 3 ]; then
    PROCS=$3
else
    PROCS=1
fi
###############################################################
MPIEXEC=sh/local/mpiexec.sh
EXE=bin/gothic
PLTENE=bin/plot.energy
PLTMAP=bin/plot.distribution
PLTMUL=bin/plot.multipole
PLTCDF=bin/plot.cdf
###############################################################
PLTSH=sh/local/plplot.sh
EXTENSION=$2
###############################################################
# NTOT=40
# NTOT=256
# NTOT=512
# NTOT=1024
# NTOT=2048
# NTOT=4096
# NTOT=8192
# NTOT=16384
# NTOT=32768
# NTOT=65536
NTOT=131072
# NTOT=262144
# NTOT=524288
# NTOT=1048576
# NTOT=2097152
# NTOT=4194304
# NTOT=8388608
# NTOT=16777216
# NTOT=33554432
# NTOT=67108864
###############################################################
# SAVE=30.0
SAVE=59.0
###############################################################
# THETA=0.9
# THETA=0.8
# THETA=0.7
# THETA=0.6
# THETA=0.5
# THETA=0.4
THETA=0.3
###############################################################
# ACCERR=1.250000e-1
# ACCERR=6.250000e-2
# ACCERR=3.125000e-2
# ACCERR=1.562500e-2
# ACCERR=7.812500e-3
ACCERR=3.906250e-3
# ACCERR=1.953125e-3
# ACCERR=9.765625e-4
###############################################################
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
###############################################################
# NCRIT=128
# NCRIT=256
NCRIT=512
# NCRIT=1024
# NCRIT=2048
###############################################################
# cold collapse of a uniform sphere
if [ $PROBLEM -eq 0 ]; then
INI=bin/uniformsphere
FILE=ccuni
MTOT=1.0
SIGMA=0.25
RAD=1.0
EPS=1.5625e-2
ETA=0.125
# FINISH=0.0
# FINISH=1.75
# FINISH=3.75
FINISH=11.75
INTERVAL=0.25
fi
###############################################################
# dynamical stability of a king sphere
if [ $PROBLEM -eq 1 ]; then
INI=bin/magi
FILE=king
CONFIG=bulge.cfg
EPS=1.5625e-2
ETA=0.125
# FINISH=3.5
FINISH=23.5
# FINISH=95.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
INI=bin/magi
FILE=hernquist
CONFIG=halo.cfg
EPS=1.5625e-2
ETA=0.125
# FINISH=0.0
# FINISH=23.5
FINISH=47.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of an NFW sphere with small truncation radius
if [ $PROBLEM -eq 3 ]; then
INI=bin/magi
FILE=nfw
CONFIG=cusp.bench.cfg
EPS=1.5625e-2
ETA=0.125
# FINISH=0.0
FINISH=3.5
# FINISH=23.5
# FINISH=47.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere and a King sphere
if [ $PROBLEM -eq 4 ]; then
INI=bin/magi
FILE=galaxy
CONFIG=galaxy.cfg
EPS=1.5625e-2
ETA=0.125
# FINISH=0.0
# FINISH=3.5
# FINISH=23.5
# FINISH=47.5
# FINISH=95.5
# INTERVAL=0.5
# FINISH=490.0
FINISH=1584.0
INTERVAL=16.0
fi
###############################################################
# dynamical stability of a Plummer profile
if [ $PROBLEM -eq 5 ]; then
INI=bin/magi
FILE=plummer
CONFIG=cluster.cfg
EPS=1.5625e-2
ETA=0.125
# FINISH=0.0
# FINISH=3.5
FINISH=23.5
INTERVAL=0.5
# INTERVAL=0.125
fi
###############################################################
# dynamical stability of a Burkert profile
if [ $PROBLEM -eq 6 ]; then
INI=bin/magi
FILE=burkert
CONFIG=dwarf.cfg
EPS=1.5625e-2
ETA=0.125
# FINISH=0.0
# FINISH=3.5
FINISH=23.5
INTERVAL=0.5
# INTERVAL=0.125
fi
###############################################################
# dynamical stability of a Moore profile
if [ $PROBLEM -eq 7 ]; then
INI=bin/magi
FILE=moore
CONFIG=fmm.cfg
EPS=1.5625e-2
ETA=0.125
# FINISH=0.0
# FINISH=3.5
FINISH=23.5
INTERVAL=0.5
# INTERVAL=0.125
fi
###############################################################
# dynamical stability of an Einasto profile
if [ $PROBLEM -eq 8 ]; then
INI=bin/magi
FILE=einasto
CONFIG=cosmo.cfg
EPS=1.5625e-2
ETA=0.125
# FINISH=0.0
# FINISH=3.5
FINISH=23.5
INTERVAL=0.5
# INTERVAL=0.125
fi
###############################################################
# dynamical stability of exponential disks (thick/thin disk) in an NFW sphere and a King sphere
if [ $PROBLEM -eq 9 ]; then
INI=bin/magi
FILE=mw
CONFIG=mw.cfg
EPS=1.5625e-2
ETA=0.125
# FINISH=23.5
# INTERVAL=0.5
FINISH=1648.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MW like galaxy (based on Widrow et al. 2008, Bedorf et al. 2014)
if [ $PROBLEM -eq 10 ]; then
INI=bin/magi
FILE=bonsai
CONFIG=bonsai.cfg
EPS=1.5625e-2
ETA=0.125
# FINISH=124.0
# INTERVAL=4.0
FINISH=1176.0
INTERVAL=8.0
fi
###############################################################
START=0
END=`echo "scale=0; $FINISH / $INTERVAL" | bc`
INCREMENT=1
###############################################################
# LOG=$JOB_NAME.o$JOB_ID
LOG=log/$FILE.l
STDOUT=log/$FILE.o$JOB_ID
STDERR=log/$FILE.e$JOB_ID
###############################################################
TIME=`date`
echo "$TIME: $INI start" >> $LOG
if [ $PROBLEM -eq 0 ]; then
echo "$INI -file=$FILE -Ntot=$NTOT -Mtot=$MTOT -sigma=$SIGMA -rad=$RAD -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>$STDOUT 2>>$STDERR" >> $LOG
$INI -file=$FILE -Ntot=$NTOT -Mtot=$MTOT -sigma=$SIGMA -rad=$RAD -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>$STDOUT 2>>$STDERR
fi
if [ $PROBLEM -ge 1 ]; then
echo "$INI -file=$FILE -config=$CONFIG -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>$STDOUT 2>>$STDERR" >> $LOG
$INI -file=$FILE -config=$CONFIG -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>$STDOUT 2>>$STDERR
fi
TIME=`date`
echo "$TIME: $INI finish" >> $LOG
###############################################################
TIME=`date`
echo "$TIME: $EXE start" >> $LOG
if [ $PROCS -eq 1 ]; then
echo "$EXE -absErr=$ABSERR -accErr=$ACCERR -theta=$THETA -file=$FILE -jobID=$JOB_ID 1>>$STDOUT 2>>$STDERR" >> $LOG
$EXE -absErr=$ABSERR -accErr=$ACCERR -theta=$THETA -file=$FILE -jobID=$JOB_ID 1>>$STDOUT 2>>$STDERR
else
echo "$MPIEXEC $EXE $PROCS $LOG -absErr=$ABSERR -accErr=$ACCERR -theta=$THETA -file=$FILE -jobID=$JOB_ID 1>>$STDOUT 2>>$STDERR" >> $LOG
$MPIEXEC $EXE $PROCS $LOG -absErr=$ABSERR -accErr=$ACCERR -theta=$THETA -file=$FILE -jobID=$JOB_ID 1>>$STDOUT 2>>$STDERR
fi
TIME=`date`
echo "$TIME: $EXE finish" >> $LOG
###############################################################
TIME=`date`
echo "$TIME: $PLTENE start" >> $LOG
echo "$PLTSH $PLTENE $FILE $EXTENSION -problem=$PROBLEM -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR" >> $LOG
$PLTSH $PLTENE $FILE $EXTENSION -problem=$PROBLEM -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR
TIME=`date`
echo "$TIME: $PLTENE finish" >> $LOG
###############################################################
if [ -e dat/$FILE.app000.dat ]; then
    echo "$TIME: $PLTCDF start" >> $LOG
    echo "$PLTSH $PLTCDF $FILE $EXTENSION -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR" >> $LOG
    $PLTSH $PLTCDF $FILE $EXTENSION -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR
    TIME=`date`
    echo "$TIME: $PLTCDF finish" >> $LOG
fi
###############################################################
echo "$TIME: $PLTMAP start" >> $LOG
echo "$PLTSH $PLTMAP $FILE $EXTENSION -problem=$PROBLEM -ncrit=$NCRIT -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR" >> $LOG
$PLTSH $PLTMAP $FILE $EXTENSION -problem=$PROBLEM -ncrit=$NCRIT -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR
TIME=`date`
echo "$TIME: $PLTMAP finish" >> $LOG
###############################################################
# echo "$TIME: $PLTMUL start" >> $LOG
# echo "$PLTSH $PLTMUL $FILE $EXTENSION -ncrit=$NCRIT -problem=$PROBLEM -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR" >> $LOG
# $PLTSH $PLTMUL $FILE $EXTENSION -problem=$PROBLEM -ncrit=$NCRIT -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR
# TIME=`date`
# echo "$TIME: $PLTMUL finish" >> $LOG
###############################################################
###############################################################
