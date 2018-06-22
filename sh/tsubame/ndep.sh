#!/bin/bash
###############################################################
GROUP=0
INDEX=0
###############################################################
EXEC=bin/gothic
PREV=magi
###############################################################
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_001k -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_002k -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_004k -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_008k -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_016k -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_032k -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_064k -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_128k -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_256k -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_512k -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_001M -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_002M -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_004M -l h_rt=0:05:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_008M -l h_rt=0:06:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_016M -l h_rt=0:15:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
qsub -g jh180045 -N job${GROUP}_${INDEX} -v ABSERR=1.9531250000000e-3,FILE=m31_032M -l h_rt=0:30:00 -hold_jid ${PREV} sh/tsubame/gothic.sh;PREV=job${GROUP}_${INDEX};INDEX=`expr $INDEX + 1`
###############################################################
exit 0
###############################################################
