#!/bin/bash
###############################################################
INDEX=0
###############################################################
DEVICE=0
# DEVICE=1
# DEVICE=2
# DEVICE=3
###############################################################
if [ "$DEVICE" -le "1" ]; then
    CPUNODE=0
else
    CPUNODE=8
fi
###############################################################
EXEC=bin/gothic
FILE=m31
###############################################################
TARGET=log/${FILE}_cpu${CPUNODE}_dev${DEVICE}
###############################################################
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=5.0000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.5000000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2500000000000e-1 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.2500000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.1250000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5625000000000e-2 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.8125000000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.9062500000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9531250000000e-3 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.7656250000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=4.8828125000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=4.8828125000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=4.8828125000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=4.8828125000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=4.8828125000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=4.8828125000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=4.8828125000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=4.8828125000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=4.8828125000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=4.8828125000000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.4414062500000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.4414062500000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.4414062500000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.4414062500000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.4414062500000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.4414062500000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.4414062500000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.4414062500000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.4414062500000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=2.4414062500000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2207031250000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2207031250000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2207031250000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2207031250000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.2207031250000e-4 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.1035156250000e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.1035156250000e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.1035156250000e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.1035156250000e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=6.1035156250000e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.0517578125000e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.0517578125000e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.0517578125000e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.0517578125000e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.0517578125000e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5258789062500e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5258789062500e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5258789062500e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5258789062500e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.5258789062500e-5 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.6293945312500e-6 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=7.6293945312500e-6 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.8146972656250e-6 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=3.8146972656250e-6 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9073486328125e-6 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=1.9073486328125e-6 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
##########
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.5367431640625e-7 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=9.5367431640625e-7 -file=${FILE} -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
###############################################################
exit 0
###############################################################
