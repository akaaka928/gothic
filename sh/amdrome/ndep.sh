#!/bin/bash
#SBATCH -J ndep             # name of job
#SBATCH -p amdrome            # partition name
#SBATCH -w amd2     # use compute node equips NVIDIA A100 PCIe
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=1
#SBATCH --gpu-freq=1410 # requested GPU frequency: 210 MHz--1410 MHz (bin = 15 MHz)
#SBATCH --mem=8192 # in units of MB


export MODULEPATH=$HOME/opt/modules:$MODULEPATH
module purge
module load openmpi
module load phdf5
module load gsl lis
module load cuda cub
cd $SLURM_SUBMIT_DIR
module list



EXEC=bin/gothic
ACC=1.953125000e-3
MODEL=m31

INDEX=0
for NUM in 001k 002k 004k 008k 016k 032k 064k 128k 256k 512k 001M 002M 004M 008M 016M 032M
do
    for ii in `seq 1 100`
    do
	TARGET=${MODEL}_$NUM
	srun numactl --cpunodebind=0 --localalloc $EXEC -absErr=$ACC -file=$TARGET -jobID=$INDEX 1>>log/$TARGET.o$INDEX 2>>log/$TARGET.e$INDEX
	INDEX=`expr $INDEX + 1`
    done
done
