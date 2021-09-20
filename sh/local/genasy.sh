#!/bin/sh
###############################################################
INI=bin/magi
###############################################################
NTOT=16777216
SAVE=55.0
###############################################################
EPS=1.5625e-2
ETA=0.5
FINISH=23.5
INTERVAL=0.5
###############################################################
FILE=a1b1g35
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g34
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g33
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g32
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g31
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g30
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g29
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g28
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g27
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g26
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
FILE=a1b1g25
numactl --localalloc $INI -file=$FILE -config=anisotropy/$FILE.cfg -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE 1>>log/$FILE.o$JOB_ID 2>>log/$FILE.e$JOB_ID
###############################################################
exit 0
###############################################################
