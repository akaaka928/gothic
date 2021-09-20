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
###############################################################
echo "#!/bin/sh" > $LIST0
echo "###############################################################" >> $LIST0
###############################################################
make clean
###############################################################
# FILE=m31
PROBLEM=27
ABSERR=1.953125000e-3
PREV=magi
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
		make gothic MEASURE_ELAPSED_TIME=1 HUNT_OPTIMAL_WALK_TREE=0 HUNT_OPTIMAL_INTEGRATE=0 HUNT_OPTIMAL_MAKE_TREE=1 HUNT_OPTIMAL_MAKE_NODE=1 HUNT_OPTIMAL_NEIGHBOUR=1 HUNT_OPTIMAL_SEPARATION=0 NUM_NTHREADS=$NTHREADS NUM_TSUB=$TSUB USE_WARPSHUFFLE=$USE_WS PREF_SHARED_MEM=$SM_PREF ADOPT_GADGET_TYPE_MAC=1 1>>$LOG 2>>$ERR
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
		    echo "qsub -g jh180045 -N job_${INDEX} -v EXEC=${EXEC},PROBLEM=${PROBLEM},ABSERR=${ABSERR} -hold_jid ${PREV} sh/tsubame/gothic.sh" >> $LIST0
		    PREV=job_${INDEX}
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
echo "exit 0" >> $LIST0
echo "###############################################################" >> $LIST0
###############################################################
chmod +x $LIST0
###############################################################
exit 0
###############################################################
