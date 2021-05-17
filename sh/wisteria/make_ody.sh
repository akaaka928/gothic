#!/bin/bash
#PJM -g gr16
#PJM -N make4ody
#PJM -L rscgrp=debug-o
#PJM -L node=1
#PJM --mpi proc=48
#PJM -L elapse=0:10:00

module purge
module load fj fjmpi
module load gsl
module load phdf5

# load myown modules
export MODULEPATH="/work/gr16/share/modules/lib:$MODULEPATH"
module load lis

make dir
make clean
make -j init extend anal m31
