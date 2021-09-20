#!/bin/sh
#SBATCH -J visualize          # name of job
#SBATCH -t 02:00:00           # upper limit of elapsed time
#SBATCH -p amdrome            # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=1

source $HOME/.bash_profile
conda activate intel
python py/plot.py
python py/draw.py
conda deactivate

exit 0
