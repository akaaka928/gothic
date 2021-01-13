#!/bin/bash
#SBATCH -J adep             # name of job
#SBATCH -p amdrome            # partition name
#SBATCH -w amd2     # use compute node equips NVIDIA A100 PCIe
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=1
#SBATCH --gpu-freq=1410 # requested GPU frequency: 210 MHz--1410 MHz (bin = 15 MHz)
#SBATCH --mem=2048 # in units of MB


export MODULEPATH=$HOME/opt/modules:$MODULEPATH
module purge
module load openmpi
module load phdf5
module load gsl lis
module load cuda cub
cd $SLURM_SUBMIT_DIR
module list



EXEC=bin/gothic
MODEL=m31

INDEX=0
for ACC in 5.0000000000000e-01 2.5000000000000e-01 1.2500000000000e-01 6.2500000000000e-02 3.1250000000000e-02 1.5625000000000e-02 7.8125000000000e-03 3.9062500000000e-03 1.9531250000000e-03 9.7656250000000e-04 4.8828125000000e-04 2.4414062500000e-04 1.2207031250000e-04 6.1035156250000e-05 3.0517578125000e-05 1.5258789062500e-05 7.6293945312500e-06 3.8146972656250e-06 1.9073486328125e-06 9.5367431640625e-07
do
    for ii in `seq 1 100`
    do
	srun numactl --cpunodebind=0 --localalloc $EXEC -absErr=$ACC -file=$MODEL -jobID=$INDEX 1>>log/$MODEL.o$INDEX 2>>log/$MODEL.e$INDEX
	INDEX=`expr $INDEX + 1`
    done
done
