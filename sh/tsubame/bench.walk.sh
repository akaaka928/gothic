#!/bin/bash
###############################################################
LOGDIR=bench
LOG=${LOGDIR}/bench.walk.log
ERR=${LOGDIR}/bench.walk.err
###############################################################
echo "benchmark on $HOSTNAME" > $LOG
echo "benchmark on $HOSTNAME" > $ERR
TIME=`date`
echo "$TIME: compilation start" >> $LOG
echo "$TIME: compilation start" >> $ERR
###############################################################
LIST0=${LOGDIR}/bench.walk0.sh
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
for NTHREADS in 512 256 1024 128
do
    ###########################################################
    # for N_block = 2 with P100
    if [ $NTHREADS -eq  128 ]; then
	# MAX_NLOOP=14
	MAX_NLOOP=6
    fi
    if [ $NTHREADS -eq  256 ]; then
	MAX_NLOOP=6
    fi
    if [ $NTHREADS -eq  512 ]; then
	MAX_NLOOP=2
    fi
    if [ $NTHREADS -eq 1024 ]; then
	# same with NTHREADS = 512
	MAX_NLOOP=1
    fi
    ###########################################################
    # # for N_block = 3 with P100
    # if [ $NTHREADS -eq  128 ]; then
    # 	# MAX_NLOOP=8
    # 	MAX_NLOOP=6
    # fi
    # if [ $NTHREADS -eq  256 ]; then
    # 	MAX_NLOOP=3
    # fi
    # if [ $NTHREADS -eq  512 ]; then
    # 	# same with NTHREADS = 256
    # 	MAX_NLOOP=1
    # fi
    # if [ $NTHREADS -eq 1024 ]; then
    # 	# same with NTHREADS = 256
    # 	MAX_NLOOP=1
    # fi
    ###########################################################
    # # for N_block = 4 with P100
    # if [ $NTHREADS -eq  128 ]; then
    # 	MAX_NLOOP=6
    # fi
    # if [ $NTHREADS -eq  256 ]; then
    # 	MAX_NLOOP=2
    # fi
    # if [ $NTHREADS -eq  512 ]; then
    # 	# same with NTHREADS = 256
    # 	MAX_NLOOP=1
    # fi
    # if [ $NTHREADS -eq 1024 ]; then
    # 	# same with NTHREADS = 256
    # 	MAX_NLOOP=1
    # fi
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
    for (( NLOOP = 1 ; NLOOP <= $MAX_NLOOP ; NLOOP += 1 ))
    do
	#######################################################
	# pad 0s for NLOOP
	digit=2
	input="`echo ${#NLOOP}`"
	if [ "$input" -le "$digit" ]; then
	    rem=`expr $digit - $input`
	    zeros=""
	    count=0
	    while [ "$count" -lt "$rem" ]; do
		zeros="`echo ${zeros}0`"
		count=`expr $count + 1`
	    done
	fi
	LOOPZEROS=$zeros
	#######################################################
	# for NWARP in 2 4 1 8 16 32
	for NWARP in 2 4 1 8
	do
	    ###################################################
	    # pad 0s for NWARP
	    digit=2
	    input="`echo ${#NWARP}`"
	    if [ "$input" -le "$digit" ]; then
		rem=`expr $digit - $input`
		zeros=""
		count=0
		while [ "$count" -lt "$rem" ]; do
		    zeros="`echo ${zeros}0`"
		    count=`expr $count + 1`
		done
	    fi
	    WARPZEROS=$zeros
	    ###################################################
	    # for TSUB in 32 16 8 4 2 1
	    for TSUB in 32 16 8
	    do
		###############################################
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
		###############################################
		# logging
		EXEC=bin/tot${TOTZEROS}${NTHREADS}sub${SUBZEROS}${TSUB}loop${LOOPZEROS}${NLOOP}ijp${WARPZEROS}${NWARP}
		echo "## generate $EXEC" >> $LOG
		echo "## generate $EXEC" >> $ERR
		###############################################
		# compile the N-body code w/ neighbor searching
		make gothic MEASURE_ELAPSED_TIME=1 HUNT_OPTIMAL_WALK_TREE=1 HUNT_OPTIMAL_INTEGRATE=1 HUNT_OPTIMAL_MAKE_TREE=0 HUNT_OPTIMAL_MAKE_NODE=0 HUNT_OPTIMAL_NEIGHBOUR=0 HUNT_OPTIMAL_SEPARATION=1 NUM_NTHREADS=$NTHREADS NUM_TSUB=$TSUB NUM_NLOOP=${NLOOP} NUM_NWARP=${NWARP} ADOPT_GADGET_TYPE_MAC=1 1>>$LOG 2>>$ERR
		###############################################
		# rename the executable
		mv bin/gothic $EXEC
		###############################################
		# clean up related object files
		rm -f bin/obj/*.o
		###############################################
		# generate job lists instead of running the execution file
		if [ -e $EXEC ]; then
		    ###########################################
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
		    ###########################################
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
