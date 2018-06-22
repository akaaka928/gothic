#!/bin/bash
###############################################################
GROUP=0
TAG=block
INDEX=0
###############################################################
EXEC=bin/gothic
PREV=magi
###############################################################
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=5.0000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.5000000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2500000000000e-1 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.2500000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.1250000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5625000000000e-2 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.8125000000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.9062500000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.7656250000000e-4 -l h_rt=0:10:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=4.8828125000000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=4.8828125000000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=4.8828125000000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=4.8828125000000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=4.8828125000000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=4.8828125000000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=4.8828125000000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=4.8828125000000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=4.8828125000000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=4.8828125000000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.4414062500000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.4414062500000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.4414062500000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.4414062500000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.4414062500000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.4414062500000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.4414062500000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.4414062500000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.4414062500000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=2.4414062500000e-4 -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2207031250000e-4 -l h_rt=0:20:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2207031250000e-4 -l h_rt=0:20:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2207031250000e-4 -l h_rt=0:20:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2207031250000e-4 -l h_rt=0:20:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.2207031250000e-4 -l h_rt=0:20:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.1035156250000e-5 -l h_rt=0:25:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.1035156250000e-5 -l h_rt=0:25:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.1035156250000e-5 -l h_rt=0:25:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.1035156250000e-5 -l h_rt=0:25:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=6.1035156250000e-5 -l h_rt=0:25:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.0517578125000e-5 -l h_rt=0:30:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.0517578125000e-5 -l h_rt=0:30:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.0517578125000e-5 -l h_rt=0:30:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.0517578125000e-5 -l h_rt=0:30:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.0517578125000e-5 -l h_rt=0:30:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5258789062500e-5 -l h_rt=1:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5258789062500e-5 -l h_rt=1:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5258789062500e-5 -l h_rt=1:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5258789062500e-5 -l h_rt=1:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.5258789062500e-5 -l h_rt=1:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.6293945312500e-6 -l h_rt=2:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=7.6293945312500e-6 -l h_rt=2:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.8146972656250e-6 -l h_rt=3:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=3.8146972656250e-6 -l h_rt=3:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9073486328125e-6 -l h_rt=4:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=1.9073486328125e-6 -l h_rt=4:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
##########
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.5367431640625e-7 -l h_rt=5:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N ${TAG}${GROUP}_${INDEX} -v ABSERR=9.5367431640625e-7 -l h_rt=5:00:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=${TAG}${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
###############################################################
exit 0
###############################################################
