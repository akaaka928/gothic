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
ABSERR=1.953125000e-3
###############################################################
TARGET=log/${FILE}_cpu${CPUNODE}_dev${DEVICE}
###############################################################
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_001k -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_002k -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_004k -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_008k -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_016k -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_032k -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_064k -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_128k -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_256k -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_512k -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_001M -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_002M -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_004M -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_008M -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_016M -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
numactl --cpunodebind=${CPUNODE} --localalloc $EXEC -absErr=${ABSERR} -file=${FILE}_032M -deviceID=${DEVICE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err;INDEX=`expr $INDEX + 1`
###############################################################
exit 0
###############################################################
