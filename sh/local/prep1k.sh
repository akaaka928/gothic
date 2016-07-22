#!/bin/sh
###############################################################
if [ $# -ne 1 ]; then
    echo "$# inputs are detected while 1 input is required to specify the test problem" 1>&2
    exit 1
fi
PROBLEM=$1
JOB_ID=$$
###############################################################
# global configurations
###############################################################
INI=bin/magi
###############################################################
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
# NTOT=134217728
# NTOT=268435456
# NTOT=536870912
# NTOT=1073741824
# NTOT=2147483648
# NTOT=4294967296
# NTOT=8589934592
SAVE=60.0
###############################################################
#
#
###############################################################
# specific configurations
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
# 6.25 Myr
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
# LOG=$JOB_NAME.o$JOB_ID
LOG=log/$FILE.l
STDOUT=log/$FILE.o$JOB_ID
STDERR=log/$FILE.e$JOB_ID
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
