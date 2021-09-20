#!/bin/bash
set -e
jobid0=`echo sleep 10 | qsub                           sh/reedbush/gothic.sh`
jobid1=`echo sleep 10 | qsub -W depend=afterok:$jobid0 sh/reedbush/gothic.sh`
jobid2=`echo sleep 10 | qsub -W depend=afterok:$jobid1 sh/reedbush/gothic.sh`
jobid3=`echo sleep 10 | qsub -W depend=afterok:$jobid2 sh/reedbush/gothic.sh`
jobid4=`echo sleep 10 | qsub -W depend=afterok:$jobid3 sh/reedbush/gothic.sh`
jobid5=`echo sleep 10 | qsub -W depend=afterok:$jobid4 sh/reedbush/gothic.sh`
jobid6=`echo sleep 10 | qsub -W depend=afterok:$jobid5 sh/reedbush/gothic.sh`
jobid7=`echo sleep 10 | qsub -W depend=afterok:$jobid6 sh/reedbush/gothic.sh`
jobid8=`echo sleep 10 | qsub -W depend=afterok:$jobid7 sh/reedbush/gothic.sh`
jobid9=`echo sleep 10 | qsub -W depend=afterok:$jobid8 sh/reedbush/gothic.sh`
exit 0
