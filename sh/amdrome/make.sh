#!/bin/bash
#SBATCH -J make_gothic # name of job
#SBATCH -t 00:10:00 # upper limit of elapsed time
#SBATCH -p amdrome  # partition name
#SBATCH -w amd2     # use compute node equips NVIDIA A100 PCIe
#SBATCH --nodes=1   # number of nodes, set to SLURM_JOB_NUM_NODES



export MODULEPATH=$HOME/opt/modules:$MODULEPATH
module purge

module load openmpi
module load phdf5
module load gsl lis
module load cuda cub

cd $SLURM_SUBMIT_DIR
module list
make
# make gothic


exit 0
