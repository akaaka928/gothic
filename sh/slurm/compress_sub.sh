#!/bin/sh
if [ $# -lt 3 ]; then
    echo "$# input(s) is/are detected while at least 3 inputs are required to specify <LIST>, <SRCDIR>, and <DSTDIR>" 1>&2
    exit 1
fi


RANK=${MV2_COMM_WORLD_RANK:=${PMI_RANK:=${OMPI_COMM_WORLD_RANK:=0}}}
SIZE=${MV2_COMM_WORLD_SIZE:=${PMI_SIZE:=${OMPI_COMM_WORLD_SIZE:=0}}}


MAX=`expr $SIZE - 1`
digit=${#MAX}
int="`echo ${#RANK}`"
if [ "$int" -le "$digit" ]; then
    rem=`expr $digit - $int`
    zeros=""
    count=0
    while [ "$count" -lt "$rem" ]; do
	zeros="`echo ${zeros}0`"
	count=`expr $count + 1`
    done
fi
LIST=${1}${zeros}${RANK}

SRCDIR=$2
DSTDIR=$3

# use numactl if available
if [ `which numactl` ]; then
    NUMACTL="numactl --localalloc"
fi

while read LINE; do
    FILE=${LINE##*/}
    SRC=$SRCDIR/$FILE
    DST=$DSTDIR/$FILE

    $NUMACTL h5repack -f SHUF -f GZIP=9 -i $SRC -o $DST
    $NUMACTL h5diff -n 1 $SRC $DST
    if [ $? -eq 0 ]; then
    	# $? is 0 if h5diff returns no error
    	$NUMACTL rm -f $SRC
    	echo $SRC is compressed: $DST | tee -a $LOG
    fi
done < $LIST
rm -f $LIST

exit 0
