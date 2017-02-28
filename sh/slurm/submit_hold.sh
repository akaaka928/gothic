#!/bin/sh
###############################################################
set -e
###############################################################
jobid0=`echo sleep 10 | sbatch                               sh/slurm/exec.sh | awk '{print $4}'`
jobid1=`echo sleep 10 | sbatch --dependency=afterany:$jobid0 sh/slurm/exec.sh | awk '{print $4}'`
jobid2=`echo sleep 10 | sbatch --dependency=afterany:$jobid1 sh/slurm/exec.sh | awk '{print $4}'`
jobid3=`echo sleep 10 | sbatch --dependency=afterany:$jobid2 sh/slurm/exec.sh | awk '{print $4}'`
jobid4=`echo sleep 10 | sbatch --dependency=afterany:$jobid3 sh/slurm/exec.sh | awk '{print $4}'`
###############################################################
