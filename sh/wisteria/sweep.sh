#!/bin/bash
#PJM -g gr16
#PJM -N sweep
#PJM -L rscgrp=share
#PJM -L gpu=1
#PJM --mpi proc=1
#PJM -L elapse=1:00:00
#PJM -s



# to get missing list:
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 5.0000000000000e-01 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 2.5000000000000e-01 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 1.2500000000000e-01 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 6.2500000000000e-02 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 3.1250000000000e-02 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 1.5625000000000e-02 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 7.8125000000000e-03 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 3.9062500000000e-03 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 9.7656250000000e-04 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 4.8828125000000e-04 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 2.4414062500000e-04 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 1.2207031250000e-05 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 6.1035156250000e-05 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 3.0517578125000e-05 | wc -l  <-- should be 8
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 7.6293945312500e-06 | wc -l  <-- should be 1
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 3.8146972656250e-06 | wc -l  <-- should be 1
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 1.9073486328125e-06 | wc -l  <-- should be 1
# $ cat log/m31_0.acceleration.block.cc80.8388608.time.log | grep 9.5367431640625e-07 | wc -l  <-- should be 1


# ABSERR=3.1250000000000e-02
# FILE=m31_2


FILE=m31_5
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_11
# FILE=m31_12


if [ -z "$ABSERR" ]; then
    ABSERR=1.953125000e-3
fi


EXEC=bin/gothic


module purge
module load cuda/11.2
module load gcc/8.3.1
module load ompi/4.1.1
module use /work/gr16/share/modules/lib
module load phdf5/1.12.0
module load cub/1.12.1
module list

STDOUT=log/${FILE}_${PJM_JOBNAME}.o${PJM_JOBID}
STDERR=log/${FILE}_${PJM_JOBNAME}.e${PJM_JOBID}


OPTION="-absErr=$ABSERR -file=$FILE -jobID=$PJM_JOBID"
numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
