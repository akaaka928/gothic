#!/bin/sh
###############################################################
if [ $# -ne 2 ]; then
    echo "$# inputs are detected while 2 inputs are required to specify the test problem" 1>&2
    exit 1
fi
PROBLEM=$1
JOB_ID=$$
###############################################################
# global configurations
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
# NTOT=131072
# NTOT=262144
# NTOT=524288
# NTOT=1048576
# NTOT=2097152
# NTOT=4194304
NTOT=8388608
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
###############################################################
PLTEXE=bin/elapsed
BRKDWN=bin/breakdown
###############################################################
PLTSH=sh/local/plplot.sh
EXTENSION=$2
###############################################################
#
#
###############################################################
# specific configurations
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
# dynamical stability of an NFW sphere with small truncation radius
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
#
#
###############################################################
# plot figures
###############################################################
# LOG=$JOB_NAME.o$JOB_ID
LOG=log/$FILE.l
STDOUT=log/$FILE.o$JOB_ID
STDERR=log/$FILE.e$JOB_ID
###############################################################
TIME=`date`
echo "$TIME: $PLTEXE start" >> $LOG
echo "$PLTSH $PLTEXE $FILE $EXTENSION -Ntot=$NTOT 1>>$STDOUT 2>>$STDERR" >> $LOG
$PLTSH $PLTEXE $FILE $EXTENSION -Ntot=$NTOT 1>>$STDOUT 2>>$STDERR
TIME=`date`
echo "$TIME: $PLTEXE finish" >> $LOG
###############################################################
TIME=`date`
echo "$TIME: $PLTEXE start" >> $LOG
echo "$BRKDWN $PLTEXE $FILE $EXTENSION 1>>$STDOUT 2>>$STDERR" >> $LOG
$BRKDWN $PLTEXE $FILE $EXTENSION 1>>$STDOUT 2>>$STDERR
TIME=`date`
echo "$TIME: $PLTEXE finish" >> $LOG
###############################################################
