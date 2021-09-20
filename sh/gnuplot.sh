#!/bin/sh
###############################################################
LOGDIR=bench
# LOGFILE=${LOGDIR}/ccuni.o5498
# LOGFILE=${LOGDIR}/ccuni.err000.dat
# LOGFILE=${LOGDIR}/nfw.time00597651.bare.dat
# LOGFILE=${LOGDIR}/m31.time00597654.bare.dat
# LOGFILE=${LOGDIR}/galaxy.time00597653.bare.dat
# LOGFILE=${LOGDIR}/plummer.time00597655.bare.dat
# LOGFILE=${LOGDIR}/016k/m31.time00608636.bare.dat
LOGFILE=${LOGDIR}/128k/m31.time00608633.bare.dat
###############################################################
# NUM=1048576
# TIME=1.667870e-01
###############################################################
BLOCK_TIME_STEP=1
###############################################################
PLTDIR=plt
# PLTFILE=time.gp
# PLTFILE=err.gp
# PLTFILE=acc.gp
# PLTFILE=compare.gp
PLTFILE=breakdown.gp
###############################################################
# gnuplot -e "filename='$LOGFILE'" -e "num=$NUM" -e "t_ful=$TIME" $PLTDIR/$PLTFILE
gnuplot -e "filename='$LOGFILE'" -e "block_time_step='$BLOCK_TIME_STEP'" $PLTDIR/$PLTFILE
###############################################################
