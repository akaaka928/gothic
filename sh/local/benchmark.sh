#!/bin/bash
###############################################################
LOG=bench/bench.log
###############################################################
echo "benchmark on $HOSTNAME" > $LOG
TIME=`date`
echo "$TIME: benchmark start" >> $LOG
###############################################################
# compile initial condition generators
make clean
make init
# clean up related object files
rm -f bin/obj/*.o
# generate initial conditions
# sh/local/init.sh 0
# sh/local/init.sh 2
sh/local/init.sh 1
sh/local/init.sh 3
###############################################################
# LOG0=bench/bench.ccuni.log
# LOG2=bench/bench.hernquist.log
LOG1=bench/bench.king.log
LOG3=bench/bench.nfw.log
###############################################################
# sample copied from sh/local/exec.sh
# ACCERR=1.562500e-2
# FILE=ccuni
# FILE=king
# FILE=hernquist
# FILE=nfw
# LOG=log/$FILE.l
# STDOUT=log/$FILE.o$JOB_ID
# STDERR=log/$FILE.e$JOB_ID
# $EXE -accErr=$ACCERR -theta=$THETA -file=$FILE -jobID=$JOB_ID 1>>$STDOUT 2>>$STDERR
###############################################################
THETA=0.6
ACCERR=3.906250e-3
ABSERR=3.906250000e-3
###############################################################
for NTHREADS in 1024 512 256 128
do
    ###########################################################
    if [ $NTHREADS -eq  128 ]; then
	MAX_NLOOP=10
    fi
    if [ $NTHREADS -eq  256 ]; then
	MAX_NLOOP=4
    fi
    if [ $NTHREADS -eq  512 ]; then
	MAX_NLOOP=1
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
	for NWARP in 1 2 4 8 16 32
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
	    # for TSUB in 32 16 8 4
	    for TSUB in 32 16
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
		# compile the N-body code w/ localized swith is offed
		make calc MEASURE_ELAPSED_TIME=1 LOCALIZE_I_PARTICLES=0 NUM_NTHREADS=$NTHREADS NUM_TSUB=$TSUB NUM_NLOOP=$NLOOP NUM_NWARP=$NWARP
		###############################################
		# rename the executable
		EXEC=bin/tot${TOTZEROS}${NTHREADS}sub${SUBZEROS}${TSUB}loop${LOOPZEROS}${NLOOP}ijp${WARPZEROS}${NWARP}
		mv bin/gothic $EXEC
		###############################################
		# clean up related object files
		rm -f bin/obj/*.o
		###############################################
		# run execution file
		$EXEC -theta=$THETA -accErr=$ACCERR -absErr=$ABSERR -file=nfw       -jobID=3 1>>$LOG3 2>>${EXEC}.err.3
		$EXEC -theta=$THETA -accErr=$ACCERR -absErr=$ABSERR -file=king      -jobID=1 1>>$LOG1 2>>${EXEC}.err.1
		# $EXEC -theta=$THETA -accErr=$ACCERR -absErr=$ABSERR -file=hernquist -jobID=2 1>>$LOG2 2>>${EXEC}.err.2
		# $EXEC -theta=$THETA -accErr=$ACCERR -absErr=$ABSERR -file=ccuni     -jobID=0 1>>$LOG0 2>>${EXEC}.err.0
		###############################################
		for (( NEIGHBOR_PHKEY_LEVEL = 1 ; NEIGHBOR_PHKEY_LEVEL <= 4 ; NEIGHBOR_PHKEY_LEVEL += 1 ))
		do
		    ###########################################
		    # pad 0s for NEIGHBOR_PHKEY_LEVEL
		    digit=2
		    input="`echo ${#NEIGHBOR_PHKEY_LEVEL}`"
		    if [ "$input" -le "$digit" ]; then
			rem=`expr $digit - $input`
			zeros=""
			count=0
			while [ "$count" -lt "$rem" ]; do
			    zeros="`echo ${zeros}0`"
			    count=`expr $count + 1`
			done
		    fi
		    JUMPZEROS=$zeros
		    ###########################################
		    # compile the N-body code w/ localized swith is offed
		    make calc MEASURE_ELAPSED_TIME=1 LOCALIZE_I_PARTICLES=1 NUM_NTHREADS=$NTHREADS NUM_TSUB=$TSUB NUM_NLOOP=$NLOOP NUM_NWARP=$NWARP LEV_NEIGHBOR=$NEIGHBOR_PHKEY_LEVEL
		    ###########################################
		    # rename the executable
		    EXEC=bin/tot${TOTZEROS}${NTHREADS}sub${SUBZEROS}${TSUB}loop${LOOPZEROS}${NLOOP}ijp${WARPZEROS}${NWARP}jump${JUMPZEROS}${NEIGHBOR_PHKEY_LEVEL}
		    mv bin/gothic $EXEC
		    ###########################################
		    # clean up related object files
		    rm -f bin/obj/*.o
		    ###########################################
		    # run execution file
		    $EXEC -theta=$THETA -accErr=$ACCERR -absErr=$ABSERR -file=nfw       -jobID=3 1>>$LOG3 2>>${EXEC}.err.3
		    $EXEC -theta=$THETA -accErr=$ACCERR -absErr=$ABSERR -file=king      -jobID=1 1>>$LOG1 2>>${EXEC}.err.1
		    # $EXEC -theta=$THETA -accErr=$ACCERR -absErr=$ABSERR -file=hernquist -jobID=2 1>>$LOG2 2>>${EXEC}.err.2
		    # $EXEC -theta=$THETA -accErr=$ACCERR -absErr=$ABSERR -file=ccuni     -jobID=0 1>>$LOG0 2>>${EXEC}.err.0
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
