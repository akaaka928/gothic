#!/bin/sh
###############################################################
EXE=bin/gothic
###############################################################
ABSERR=1.953125e-3
###############################################################
FILE=a1b1g35
JOB_ID=0
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g34
JOB_ID=1
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g33
JOB_ID=2
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g32
JOB_ID=3
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g31
JOB_ID=4
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g30
JOB_ID=5
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g29
JOB_ID=6
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g28
JOB_ID=7
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g27
JOB_ID=8
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g26
JOB_ID=9
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g25
JOB_ID=10
numactl --cpunodebind=0 --localalloc $EXE -file=$FILE -absErr=$ABSERR -jobID=$JOB_ID 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
###############################################################
exit 0
###############################################################
