#!/bin/sh
###############################################################
set -e
###############################################################
sbatch                        sh/slurm/exec.sh
sbatch --dependency=singleton sh/slurm/exec.sh
sbatch --dependency=singleton sh/slurm/exec.sh
sbatch --dependency=singleton sh/slurm/exec.sh
sbatch --dependency=singleton sh/slurm/exec.sh
###############################################################
