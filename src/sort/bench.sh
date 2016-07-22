#!/bin/bash
###############################################################
LIST=bench/exec.txt
LOG=bench/bench.log
ERR=bench/bench.err
###############################################################
echo "benchmark on $HOSTNAME" > $LOG
TIME=`date`
echo "$TIME: benchmark start" >> $LOG
###############################################################
# clean up
make clean
make dir
rm -f $LIST
touch $LIST
###############################################################
echo "#!/bin/sh" > $LIST
###############################################################
for NTHREADS in 1024 512 256 128 64 32
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
    for MODE in 0 1 2 3 4
    do
	#######################################################
	for NKEYS in 4 2 1
	do
	    ###################################################
	    for MANY_BITS in 0 1
	    do
		###############################################
		for SMPREF in 0 1
		do
		    ###########################################
		    for WARP_SHUFFLE in 0 1
		    do
			#######################################
			# compile the radix sort code
			make HUNT_OPTIMAL_PARAMETER_SETS=1 NUM_NTHREADS=$NTHREADS NUM_ELEMENTS=$NKEYS CHECK_MORE_BITS_PER_STEP=$MANY_BITS PREFERR_SHARED_MEMORY_GLOBAL=$SMPREF USE_WARP_SHUFFLE_INSTRUCTION=$WARP_SHUFFLE ALGORITHM=$MODE
			#######################################
			# rename the executable
			if [ -e bin/radix ]; then
			    EXEC=bin/thrd${TOTZEROS}${NTHREADS}key${NKEYS}bit${MANY_BITS}SM${SMPREF}WS${WARP_SHUFFLE}mode${MODE}
			    mv bin/radix $EXEC
			    echo "$EXEC 1>>$LOG 2>>$ERR" >> $LIST
			fi
			#######################################
			# clean up related object files
			rm -f bin/obj/*.o
			#######################################
			# # run execution file
			# $EXEC 1>>$LOG 2>>$ERR
			#######################################
		    done
		    ###########################################
		done
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
echo "$TIME: benchmark finish" >> $LOG
###############################################################
exit 0
###############################################################
