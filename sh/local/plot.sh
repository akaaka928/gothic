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
PLTENE=bin/energy
PLTMAP=bin/distribution
PLTMUL=bin/multipole
PLTCDF=bin/cdf
###############################################################
PLTSH=sh/local/plplot.sh
EXTENSION=$2
###############################################################
# set # of N-body particles per bin to estimate density
# NCRIT=8
# NCRIT=32
# NCRIT=64
# NCRIT=128
# NCRIT=256
# NCRIT=512
# NCRIT=1024
NCRIT=2048
###############################################################
#
#
###############################################################
# specific configurations
###############################################################
# cold collapse of a uniform sphere
if [ $PROBLEM -eq 0 ]; then
FILE=ccuni
FINISH=11.75
INTERVAL=0.25
fi
###############################################################
# dynamical stability of a King sphere
if [ $PROBLEM -eq 1 ]; then
FILE=king
FINISH=23.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
FILE=hernquist
FINISH=23.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of an NFW sphere with small truncation radius
if [ $PROBLEM -eq 3 ]; then
FILE=nfw
FINISH=23.5
INTERVAL=0.5
# FINISH=0.1
# INTERVAL=0.1
fi
###############################################################
# dynamical stability of an Einasto profile
if [ $PROBLEM -eq 4 ]; then
FILE=einasto
FINISH=23.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Plummer profile
if [ $PROBLEM -eq 5 ]; then
FILE=plummer
FINISH=23.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Burkert profile
if [ $PROBLEM -eq 6 ]; then
FILE=burkert
FINISH=23.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Moore profile
if [ $PROBLEM -eq 7 ]; then
FILE=moore
FINISH=23.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Two-power sphere
if [ $PROBLEM -eq 8 ]; then
FILE=twopower
FINISH=23.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of a King sphere within an Einasto profile
if [ $PROBLEM -eq 10 ]; then
FILE=hb
FINISH=23.5
INTERVAL=0.5
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere and a King sphere
if [ $PROBLEM -eq 11 ]; then
FILE=hbd
FINISH=47.5
INTERVAL=0.5
# FINISH=1648.0
# INTERVAL=16.0
fi
###############################################################
# dynamical stability of exponential disks (thick/thin disk) in an NFW sphere and a King sphere
if [ $PROBLEM -eq 12 ]; then
FILE=hbdd
FINISH=47.5
INTERVAL=0.5
# FINISH=1400.0
# INTERVAL=8.0
fi
###############################################################
# dynamical stability of thick exponential disk and thin Sersic disk in an Einasto sphere and a King sphere
if [ $PROBLEM -eq 13 ]; then
FILE=ekes
FINISH=47.5
INTERVAL=0.5
# FINISH=1400.0
# INTERVAL=8.0
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007)
if [ $PROBLEM -eq 20 ]; then
FILE=m31
# FINISH=75.0
# INTERVAL=25.0
FINISH=1175.0
INTERVAL=25.0
# FINISH=5175.0
# INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model
if [ $PROBLEM -eq 22 ]; then
FILE=galaxy
FINISH=75.0
INTERVAL=25.0
# FINISH=3175.0
# INTERVAL=25.0
fi
###############################################################
# dynamical stability of MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008)
if [ $PROBLEM -eq 23 ]; then
FILE=mw
FINISH=75.0
INTERVAL=25.0
# FINISH=3175.0
# INTERVAL=25.0
fi
###############################################################
# dynamical stability of M31 model (Sofue 2015; Gilbert et al. 2012)
if [ $PROBLEM -eq 24 ]; then
FILE=s15
FINISH=75.0
INTERVAL=25.0
# FINISH=3175.0
# INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model with fixed number of particles
if [ $PROBLEM -eq 25 ]; then
FILE=compare
# FINISH=75.0
# INTERVAL=25.0
FINISH=3175.0
INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (only spherical components)
if [ $PROBLEM -eq 26 ]; then
FILE=etg
# FINISH=75.0
# INTERVAL=25.0
FINISH=1175.0
INTERVAL=25.0
# FINISH=3175.0
# INTERVAL=25.0
fi
###############################################################
# dynamical stability of an M31 model (NFW halo, de Vaucouleurs bulge, and exponential disk)
if [ $PROBLEM -eq 27 ]; then
FILE=m31_mod
FINISH=75.0
INTERVAL=25.0
# FINISH=3175.0
# INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk)
if [ $PROBLEM -eq 28 ]; then
FILE=ltg
# FINISH=75.0
# INTERVAL=25.0
FINISH=1175.0
INTERVAL=25.0
# FINISH=3175.0
# INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (Vasiliev & Athanassoula 2015) with fixed number of particles
if [ $PROBLEM -eq 29 ]; then
if [ $NTOT -ne 1000000 ]; then
    NTOT=1000000
fi
FILE=va15
FINISH=100.0
INTERVAL=4.0
fi
###############################################################
# time evolution of MW/A defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 30 ]; then
FILE=kd95a
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MW/B defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 31 ]; then
FILE=kd95b
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MW/C defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 32 ]; then
FILE=kd95c
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MW/D defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 33 ]; then
FILE=kd95d
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of M31/A defined in Widrow et al. (2003)
if [ $PROBLEM -eq 34 ]; then
FILE=w03a
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of M31/D defined in Widrow et al. (2003)
if [ $PROBLEM -eq 35 ]; then
FILE=w03d
# FINISH=124.0
# INTERVAL=4.0
FINISH=2032.0
INTERVAL=16.0
fi
###############################################################
# time evolution of MWa defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 36 ]; then
FILE=mwa
FINISH=124.0
INTERVAL=4.0
# FINISH=1016.0
# INTERVAL=8.0
fi
###############################################################
# time evolution of MWb defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 37 ]; then
FILE=mwb
FINISH=124.0
INTERVAL=4.0
# FINISH=1016.0
# INTERVAL=8.0
fi
###############################################################
# time evolution of MW like galaxy (based on Widrow et al. 2008, Bedorf et al. 2014)
if [ $PROBLEM -eq 38 ]; then
FILE=bonsai
FINISH=124.0
INTERVAL=4.0
# FINISH=1176.0
# INTERVAL=8.0
fi
###############################################################
# dynamical stability of a Plummer profile in table form
if [ $PROBLEM -eq 40 ]; then
FILE=tplummer
FINISH=1575.0
INTERVAL=25.0
fi
###############################################################
# dynamical stability of a Two-power slope profile in table form
if [ $PROBLEM -eq 41 ]; then
FILE=dblpower
FINISH=1575.0
INTERVAL=25.0
fi
###############################################################
# dynamical stability of a de Vaucouleurs sphere in table form
if [ $PROBLEM -eq 42 ]; then
FILE=bulgetbl
FINISH=1575.0
INTERVAL=25.0
fi
###############################################################
# dynamical stability of a spherical Sersic profile
if [ $PROBLEM -eq 50 ]; then
FILE=deVaucouleurs
FINISH=1575.0
INTERVAL=25.0
fi
###############################################################
# dynamical stability of a projected Two-power model
if [ $PROBLEM -eq 51 ]; then
FILE=prjTwoPow
FINISH=1575.0
INTERVAL=25.0
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
# LOG=$JOB_NAME.o$JOB_ID
LOG=log/$FILE.l
STDOUT=log/$FILE.o$JOB_ID
STDERR=log/$FILE.e$JOB_ID
###############################################################
TIME=`date`
echo "$TIME: $PLTENE start" >> $LOG
echo "$PLTSH $PLTENE $FILE $EXTENSION -problem=$PROBLEM -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR" >> $LOG
$PLTSH $PLTENE $FILE $EXTENSION -problem=$PROBLEM -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR
TIME=`date`
echo "$TIME: $PLTENE finish" >> $LOG
###############################################################
if [ -e dat/$FILE.direct000.dat ]; then
    echo "$TIME: $PLTCDF start" >> $LOG
    echo "$PLTSH $PLTCDF $FILE $EXTENSION 1>>$STDOUT 2>>$STDERR" >> $LOG
    $PLTSH $PLTCDF $FILE $EXTENSION 1>>$STDOUT 2>>$STDERR
    TIME=`date`
    echo "$TIME: $PLTCDF finish" >> $LOG
fi
###############################################################
echo "$TIME: $PLTMAP start" >> $LOG
echo "$PLTSH $PLTMAP $FILE $EXTENSION -ncrit=$NCRIT -problem=$PROBLEM -start=$START -end=$END -interval=$INCREMENT 1>>$STDOUT 2>>$STDERR" >> $LOG
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
