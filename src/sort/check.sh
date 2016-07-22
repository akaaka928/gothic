#!/bin/sh
#########################################
LOGDIR=log
#########################################
SRC=$LOGDIR/input.dat
CHK=$LOGDIR/cmd.dat
DST=$LOGDIR/output.dat
#########################################
sort -n -k 1 $SRC -o $CHK
diff $DST $CHK
#########################################
