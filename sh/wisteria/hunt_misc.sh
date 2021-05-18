#!/bin/sh

# modules settings for Aquarius
module purge
module load cuda
module load gcc ompi-cuda
export MODULEPATH="/work/gr16/share/modules/lib:$MODULEPATH"
module load phdf5
module load cub
module list


LOGDIR=bench
LOG=${LOGDIR}/make_misc.log
ERR=${LOGDIR}/make_misc.err
make clean


LIST0=${LOGDIR}/bench_misc0.sh
LIST1=${LOGDIR}/bench_misc1.sh
LIST2=${LOGDIR}/bench_misc2.sh
LIST3=${LOGDIR}/bench_misc3.sh
LIST4=${LOGDIR}/bench_misc4.sh
LIST5=${LOGDIR}/bench_misc5.sh
LIST6=${LOGDIR}/bench_misc6.sh
LIST7=${LOGDIR}/bench_misc7.sh


for TARGET in $LIST0 $LIST1 $LIST2 $LIST3 $LIST4 $LIST5 $LIST6 $LIST7
do
    echo "#!/bin/bash" > $TARGET
    echo "#PJM -g gr16" >> $TARGET
    echo "#PJM -N hunt_misc" >> $TARGET
    echo "#PJM -L rscgrp=share" >> $TARGET
    echo "#PJM -L gpu=1" >> $TARGET
    echo "#PJM -L elapse=0:03:00" >> $TARGET
    echo "" >> $TARGET
    echo "module purge" >> $TARGET
    echo "module load cuda gcc ompi-cuda" >> $TARGET
    echo "export MODULEPATH=/work/gr16/share/modules/lib:\$MODULEPATH" >> $TARGET
    echo "module load phdf5 cub" >> $TARGET
    echo "module list" >> $TARGET
    echo "" >> $TARGET
    echo "cd \$PJM_O_WORKDIR" >> $TARGET
done



FILE=m31
ABSERR=1.953125000e-3


INDEX=0
for NTHREADS in 512 256 128 1024
do
    # pad 0s for NTHREADS
    digit=4
    input="`echo ${#NTHREADS}`"
    if [ "$input" -le "$digit" ]; then
	rem=`expr $digit - $input`
	zeros=""
	count=0
	while [ "$count" -lt "$rem" ]; do
	    zeros="`echo ${zeros}0`"
	    count=`expr $count + 1`
	done
    fi
    TOTZEROS=$zeros

    for TSUB in 32 16 8 4 2 1
    do
	# pad 0s for TSUB
	digit=2
	input="`echo ${#TSUB}`"
	if [ "$input" -le "$digit" ]; then
	    rem=`expr $digit - $input`
	    zeros=""
	    count=0
	    while [ "$count" -lt "$rem" ]; do
		zeros="`echo ${zeros}0`"
		count=`expr $count + 1`
	    done
	fi
	SUBZEROS=$zeros

	for USE_WS in 1 0
	do

	    for USE_WR in 1 0
	    do

		for SM_PREF in 1 0 # only for HUNT_OPTIMAL_NEIGHBOUR=1
		do

		    # logging
		    EXEC=bin/tot${TOTZEROS}${NTHREADS}sub${SUBZEROS}${TSUB}ws${USE_WS}wr${USE_WR}sm${SM_PREF}
		    echo "## generate $EXEC" >> $LOG

		    make gothic MEASURE_ELAPSED_TIME=1 HUNT_OPTIMAL_WALK_TREE=0 HUNT_OPTIMAL_INTEGRATE=0 HUNT_OPTIMAL_MAKE_TREE=1 HUNT_OPTIMAL_MAKE_NODE=1 HUNT_OPTIMAL_NEIGHBOUR=1 HUNT_OPTIMAL_SEPARATION=0 NUM_NTHREADS=$NTHREADS NUM_TSUB=$TSUB USE_WARPSHUFFLE=$USE_WS USE_WARPREDUCE=$USE_WR PREF_SHARED_MEM=$SM_PREF ADOPT_GADGET_TYPE_MAC=1 1>>$LOG 2>>$ERR

		    # rename the executable
		    mv bin/gothic $EXEC
		    make clean

		    # generate job lists instead of running the execution file
		    if [ -e $EXEC ]; then
			TARGET=log/${EXEC##*/}
			kind=`expr $INDEX % 8`
			if [ "$kind" -eq "0" ]; then
			    echo "numactl --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err" >> $LIST0
			fi
			if [ "$kind" -eq "1" ]; then
			    echo "numactl --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err" >> $LIST1
			fi
			if [ "$kind" -eq "2" ]; then
			    echo "numactl --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err" >> $LIST2
			fi
			if [ "$kind" -eq "3" ]; then
			    echo "numactl --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err" >> $LIST3
			fi
			if [ "$kind" -eq "4" ]; then
			    echo "numactl --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err" >> $LIST4
			fi
			if [ "$kind" -eq "5" ]; then
			    echo "numactl --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err" >> $LIST5
			fi
			if [ "$kind" -eq "6" ]; then
			    echo "numactl --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err" >> $LIST6
			fi
			if [ "$kind" -eq "7" ]; then
			    echo "numactl --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -jobID=$INDEX 1>>${TARGET}_${INDEX}.log 2>>${TARGET}_${INDEX}.err" >> $LIST7
			fi
			INDEX=`expr $INDEX + 1`
		    fi
		done
	    done
	done
    done
done


for TARGET in $LIST0 $LIST1 $LIST2 $LIST3 $LIST4 $LIST5 $LIST6 $LIST7
do
    chmod +x $TARGET
done

exit 0
