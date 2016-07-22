#!/bin/sh
###############################################################
set -e
###############################################################
jobid0=`echo sleep 10 | qsub -N test0                            -v FILE=nfw,EXEC=bin/tot0256sub16loop01ijp02 sh/comq/exec.sh`
jobid1=`echo sleep 10 | qsub -N test1 -W depend=afterany:$jobid0 -v FILE=nfw,EXEC=bin/tot0256sub16loop01ijp04 sh/comq/exec.sh`
jobid2=`echo sleep 10 | qsub -N test2 -W depend=afterany:$jobid1 -v FILE=nfw,EXEC=bin/tot0256sub16loop01ijp08 sh/comq/exec.sh`
###############################################################
