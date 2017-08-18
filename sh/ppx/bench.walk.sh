#!/bin/bash
###############################################################
LOGDIR=bench
###############################################################
LIST=${LOGDIR}/bench.walk.sh
###############################################################
echo "#!/bin/sh" > $LIST
echo "###############################################################" >> $LIST
echo "set -e" >> $LIST
echo "###############################################################" >> $LIST
###############################################################
make clean
###############################################################
INDEX=0
PREV=0
###############################################################
# for NTHREADS in 512 256 1024 128
for NTHREADS in 256 512 128
# for NTHREADS in 512 1024
do
    ###########################################################
    if [ $NTHREADS -eq  128 ]; then
	MAX_NLOOP=8
	# MAX_NLOOP=10
	# MAX_NLOOP=14
    fi
    if [ $NTHREADS -eq  256 ]; then
	MAX_NLOOP=6
    fi
    if [ $NTHREADS -eq  512 ]; then
	MAX_NLOOP=2
    fi
    if [ $NTHREADS -eq 1024 ]; then
	MAX_NLOOP=1
    fi
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
	# for NWARP in 2
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
	    # for TSUB in 32 16
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
		# pad 0s for INDEX
		digit=3
		input="`echo ${#INDEX}`"
		if [ "$input" -le "$digit" ]; then
		    rem=`expr $digit - $input`
		    zeros=""
		    count=0
		    while [ "$count" -lt "$rem" ]; do
			zeros="`echo ${zeros}0`"
			count=`expr $count + 1`
		    done
		fi
		IDXZEROS=$zeros
		###############################################
		# logging
		NAME=walk${IDXZEROS}${INDEX}_tot${TOTZEROS}${NTHREADS}sub${SUBZEROS}${TSUB}loop${LOOPZEROS}${NLOOP}ijp${WARPZEROS}${NWARP}
		EXEC=bin/${NAME}
		###############################################
		# compile the N-body code w/ neighbor searching
		make -j gothic MEASURE_ELAPSED_TIME=1 HUNT_OPTIMAL_WALK_TREE=1 HUNT_OPTIMAL_INTEGRATE=1 HUNT_OPTIMAL_MAKE_TREE=0 HUNT_OPTIMAL_MAKE_NODE=0 HUNT_OPTIMAL_NEIGHBOUR=0 HUNT_OPTIMAL_SEPARATION=1 NUM_NTHREADS=$NTHREADS NUM_TSUB=$TSUB NUM_NLOOP=${NLOOP} NUM_NWARP=${NWARP} ADOPT_GADGET_TYPE_MAC=1 ADOPT_WS93_TYPE_MAC=1 1>${LOGDIR}/${NAME}.log 2>${LOGDIR}/${NAME}.err
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
		    # if [ $INDEX != $PREV ]; then
		    if [ `expr $INDEX % 5` = 0 ]; then
			echo "#" >> $LIST
			echo "# jobid${INDEX}=\`echo sleep 10 | sbatch                               --job-name=${NAME} --export=PROBLEM=20,EXEC=${EXEC} sh/ppx/exec.sh | awk '{print \$4}'\`" >> $LIST
		    else
			echo "# jobid${INDEX}=\`echo sleep 10 | sbatch --dependency=afterany:\$jobid${PREV} --job-name=${NAME} --export=PROBLEM=20,EXEC=${EXEC} sh/ppx/exec.sh | awk '{print \$4}'\`" >> $LIST
		    fi
		    PREV=$INDEX
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
echo "###############################################################" >> $LIST
echo "exit 0" >> $LIST
echo "###############################################################" >> $LIST
###############################################################
chmod +x $LIST
###############################################################
exit 0
###############################################################
