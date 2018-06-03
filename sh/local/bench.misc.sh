#!/bin/bash
###############################################################
LOGDIR=bench
LOG=${LOGDIR}/bench.misc.log
ERR=${LOGDIR}/bench.misc.err
###############################################################
echo "benchmark on $HOSTNAME" > $LOG
echo "benchmark on $HOSTNAME" > $ERR
TIME=`date`
echo "$TIME: compilation start" >> $LOG
echo "$TIME: compilation start" >> $ERR
###############################################################
LIST0=${LOGDIR}/bench.misc0.sh
LIST1=${LOGDIR}/bench.misc1.sh
LIST2=${LOGDIR}/bench.misc2.sh
LIST3=${LOGDIR}/bench.misc3.sh
LIST4=${LOGDIR}/bench.misc4.sh
LIST5=${LOGDIR}/bench.misc5.sh
LIST6=${LOGDIR}/bench.misc6.sh
LIST7=${LOGDIR}/bench.misc7.sh
###############################################################
echo "#!/bin/sh" > $LIST0
echo "#!/bin/sh" > $LIST1
echo "#!/bin/sh" > $LIST2
echo "#!/bin/sh" > $LIST3
echo "#!/bin/sh" > $LIST4
echo "#!/bin/sh" > $LIST5
echo "#!/bin/sh" > $LIST6
echo "#!/bin/sh" > $LIST7
###############################################################
echo "###############################################################" >> $LIST0
echo "###############################################################" >> $LIST1
echo "###############################################################" >> $LIST2
echo "###############################################################" >> $LIST3
echo "###############################################################" >> $LIST4
echo "###############################################################" >> $LIST5
echo "###############################################################" >> $LIST6
echo "###############################################################" >> $LIST7
###############################################################
echo "ulimit -t 3600" >> $LIST0
echo "ulimit -t 3600" >> $LIST1
echo "ulimit -t 3600" >> $LIST2
echo "ulimit -t 3600" >> $LIST3
echo "ulimit -t 3600" >> $LIST4
echo "ulimit -t 3600" >> $LIST5
echo "ulimit -t 3600" >> $LIST6
echo "ulimit -t 3600" >> $LIST7
###############################################################
echo "###############################################################" >> $LIST0
echo "###############################################################" >> $LIST1
echo "###############################################################" >> $LIST2
echo "###############################################################" >> $LIST3
echo "###############################################################" >> $LIST4
echo "###############################################################" >> $LIST5
echo "###############################################################" >> $LIST6
echo "###############################################################" >> $LIST7
###############################################################
make clean
###############################################################
FILE=m31
ABSERR=1.953125000e-3
###############################################################
INDEX=0
###############################################################
for NTHREADS in 1024 512 256 128
do
    ###########################################################
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
    ###########################################################
    for TSUB in 32 16 8 4 2 1
    do
	#######################################################
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
	#######################################################
	for USE_WS in 1 0
	do
	    ###################################################
	    for SM_PREF in 1 0
	    do
		###############################################
		# logging
		EXEC=bin/tot${TOTZEROS}${NTHREADS}sub${SUBZEROS}${TSUB}ws${USE_WS}sm${SM_PREF}
		echo "## generate $EXEC" >> $LOG
		echo "## generate $EXEC" >> $ERR
		###############################################
		# compile the N-body code w/ neighbor searching
		make -j gothic MEASURE_ELAPSED_TIME=1 HUNT_OPTIMAL_WALK_TREE=0 HUNT_OPTIMAL_INTEGRATE=0 HUNT_OPTIMAL_MAKE_TREE=1 HUNT_OPTIMAL_MAKE_NODE=1 HUNT_OPTIMAL_NEIGHBOUR=1 HUNT_OPTIMAL_SEPARATION=0 NUM_NTHREADS=$NTHREADS NUM_TSUB=$TSUB USE_WARPSHUFFLE=$USE_WS PREF_SHARED_MEM=$SM_PREF ADOPT_GADGET_TYPE_MAC=1 1>>$LOG 2>>$ERR
		###############################################
		# rename the executable
		mv bin/gothic $EXEC
		###############################################
		# clean up related object files
		rm -f bin/obj/*.o
		###############################################
		# generate job lists instead of running the execution file
		if [ -e $EXEC ]; then
		    # echo $LINE       # ファイルに書いてある内容を１行ずつ取得
		    # echo ${LINE##*/} # $LINE のファイル名を取得
		    # echo ${LINE%.*}  # $LINE の拡張子を除いたファイル名（フォルダ込）を取得
		    # echo ${LINE##*.} # $LINE の拡張子を取得
		    # echo ${LINE%/*}  # $LINE のフォルダ名を取得
		    # echo ${LINE##*/} # $LINE のファイル名を取得
		    TARGET=log/${EXEC##*/}
		    kind=`expr $INDEX % 8`
		    if [ "$kind" -eq "0" ]; then
			echo "numactl --cpunodebind=0 --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -deviceID=0 -jobID=$INDEX 1>>${TARGET}.${INDEX}.log 2>>${TARGET}.${INDEX}.err" >> $LIST0
		    fi
		    if [ "$kind" -eq "1" ]; then
			echo "numactl --cpunodebind=0 --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -deviceID=1 -jobID=$INDEX 1>>${TARGET}.${INDEX}.log 2>>${TARGET}.${INDEX}.err" >> $LIST1
		    fi
		    if [ "$kind" -eq "2" ]; then
			echo "numactl --cpunodebind=8 --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -deviceID=2 -jobID=$INDEX 1>>${TARGET}.${INDEX}.log 2>>${TARGET}.${INDEX}.err" >> $LIST2
		    fi
		    if [ "$kind" -eq "3" ]; then
			echo "numactl --cpunodebind=8 --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -deviceID=3 -jobID=$INDEX 1>>${TARGET}.${INDEX}.log 2>>${TARGET}.${INDEX}.err" >> $LIST3
		    fi
		    if [ "$kind" -eq "4" ]; then
			echo "numactl --cpunodebind=0 --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -deviceID=0 -jobID=$INDEX 1>>${TARGET}.${INDEX}.log 2>>${TARGET}.${INDEX}.err" >> $LIST4
		    fi
		    if [ "$kind" -eq "5" ]; then
			echo "numactl --cpunodebind=0 --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -deviceID=1 -jobID=$INDEX 1>>${TARGET}.${INDEX}.log 2>>${TARGET}.${INDEX}.err" >> $LIST5
		    fi
		    if [ "$kind" -eq "6" ]; then
			echo "numactl --cpunodebind=8 --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -deviceID=2 -jobID=$INDEX 1>>${TARGET}.${INDEX}.log 2>>${TARGET}.${INDEX}.err" >> $LIST6
		    fi
		    if [ "$kind" -eq "7" ]; then
			echo "numactl --cpunodebind=8 --localalloc $EXEC -absErr=$ABSERR -file=${FILE} -deviceID=3 -jobID=$INDEX 1>>${TARGET}.${INDEX}.log 2>>${TARGET}.${INDEX}.err" >> $LIST7
		    fi
		    INDEX=`expr $INDEX + 1`
		fi
		###############################################
	    done
	    ###################################################
	done
	#######################################################
    done
    ###########################################################
done
###############################################################
TIME=`date`
echo "$TIME: compilation finish" >> $LOG
echo "$TIME: compilation finish" >> $ERR
###############################################################
echo "###############################################################" >> $LIST0
echo "###############################################################" >> $LIST1
echo "###############################################################" >> $LIST2
echo "###############################################################" >> $LIST3
echo "###############################################################" >> $LIST4
echo "###############################################################" >> $LIST5
echo "###############################################################" >> $LIST6
echo "###############################################################" >> $LIST7
###############################################################
echo "exit 0" >> $LIST0
echo "exit 0" >> $LIST1
echo "exit 0" >> $LIST2
echo "exit 0" >> $LIST3
echo "exit 0" >> $LIST4
echo "exit 0" >> $LIST5
echo "exit 0" >> $LIST6
echo "exit 0" >> $LIST7
###############################################################
echo "###############################################################" >> $LIST0
echo "###############################################################" >> $LIST1
echo "###############################################################" >> $LIST2
echo "###############################################################" >> $LIST3
echo "###############################################################" >> $LIST4
echo "###############################################################" >> $LIST5
echo "###############################################################" >> $LIST6
echo "###############################################################" >> $LIST7
###############################################################
chmod +x $LIST0
chmod +x $LIST1
chmod +x $LIST2
chmod +x $LIST3
chmod +x $LIST4
chmod +x $LIST5
chmod +x $LIST6
chmod +x $LIST7
###############################################################
exit 0
###############################################################
