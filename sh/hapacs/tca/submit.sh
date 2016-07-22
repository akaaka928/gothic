#!/bin/sh
###############################################################
set -e
###############################################################
PROBLEM=1
###############################################################
jobid0=`echo sleep 10 | qsub -N plot.energy                             -v PROBLEM=${PROBLEM},TARGET=0 sh/hapacs/tca/plot.sh`
jobid1=`echo sleep 10 | qsub -N plot.profile -W depend=afterany:$jobid0 -v PROBLEM=${PROBLEM},TARGET=1 sh/hapacs/tca/plot.sh`
# jobid2=`echo sleep 10 | qsub -N plot.cdf     -W depend=afterany:$jobid1 -v PROBLEM=${PROBLEM},TARGET=2 sh/hapacs/tca/plot.sh`
###############################################################
