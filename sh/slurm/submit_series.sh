#!/bin/sh
###############################################################
set -e
###############################################################
PROBLEM=20
# PROBLEM=28
# PROBLEM=81
###############################################################
jobid0=`echo sleep 10 | sbatch                              --export=PROBLEM=${PROBLEM} sh/slurm/init.sh | awk '{print $4}'`
jobid1=`echo sleep 10 | sbatch --dependency=afterok:$jobid0 --export=PROBLEM=${PROBLEM} sh/slurm/exec.sh | awk '{print $4}'`
jobid2=`echo sleep 10 | sbatch --dependency=afterok:$jobid1 --export=PROBLEM=${PROBLEM} sh/slurm/plot.sh | awk '{print $4}'`
###############################################################
