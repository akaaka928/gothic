#!/bin/sh

# modules settings for Aquarius
module purge
module load cuda
module load gcc ompi-cuda

# load my modules
module use /work/gr16/share/modules/lib
module load phdf5
module load cub
module list


LOGDIR=bench
LOG=${LOGDIR}/make_walk.log
ERR=${LOGDIR}/make_walk.err
make clean


# specification of NVIDIA A100
SHARED_MEM_PER_SM=167936 # in units of byte, for NVIDIA A100
MAX_REGISTERS_PER_SM=65536
MAX_BLOCKS_PER_SM=32

# typical value
REGISTERS_PER_BLOCK=64 # typical value, not always correct


# compilation
LIST="EXEC_ALL"
# INDEX=0
# for NTHREADS in 1024
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

    # maximum number of blocks per SM (based on GPU spec)
    MAX_NBLOCKS_SM=$MAX_BLOCKS_PER_SM
    MIN_NBLOCKS_SM=2

    # maximum number of blocks per SM (based on register usage)
    MAX_BLOCKS_PER_SM_REGISTER=`echo "scale=0; (2 * $MAX_REGISTERS_PER_SM) / ($REGISTERS_PER_BLOCK * $NTHREADS)" | bc`# 2 is safety factor (REGISTERS_PER_BLOCK is just a typical value)
    if [ $MAX_NBLOCKS_SM -gt $MAX_BLOCKS_PER_SM_REGISTER ]; then
		MAX_NBLOCKS_SM=$MAX_BLOCKS_PER_SM_REGISTER
    fi

    for USE_WS in 1 0
    do
		# maximum number of blocks per SM (based on shared memory usage)
		MAX_BLOCKS_PER_SM_SHARED_MEM=`echo "scale=0; $SHARED_MEM_PER_SM / (4 * $NTHREADS * (10 - $USE_WS))" | bc`
		if [ $MAX_NBLOCKS_SM -gt $MAX_BLOCKS_PER_SM_SHARED_MEM ]; then
			MAX_NBLOCKS_SM=$MAX_BLOCKS_PER_SM_SHARED_MEM
		fi

		if [ $MAX_NBLOCKS_SM -gt 4 ]; then
			MAX_NBLOCKS_SM=4
		fi
		if [ $MIN_NBLOCKS_SM -gt $MAX_NBLOCKS_SM ]; then
			MIN_NBLOCKS_SM=1
		fi

		for (( NBLOCKS_SM = $MIN_NBLOCKS_SM ; NBLOCKS_SM <= $MAX_NBLOCKS_SM ; NBLOCKS_SM += 1 ))
		do
			MAX_NLOOP=`echo "scale=0; ($SHARED_MEM_PER_SM / (4 * $NTHREADS * $NBLOCKS_SM)) - (6 - $USE_WS)" | bc`
			if [ $MAX_NLOOP -gt 4 ]; then
				MAX_NLOOP=4
			fi

			for (( NLOOP = 1 ; NLOOP <= $MAX_NLOOP ; NLOOP += 1 ))
			do
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

				# for NWARP in 2 4 1 8 16 32
				# for NWARP in 2 4 1 8
				for NWARP in 2 1 4
				do
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

					# for TSUB in 32 16 8 4 2 1
					for TSUB in 32 16
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

						for USE_WR in 1 0
						do

							for USE_L2_ASIDE in 1 0
							do
								MIN_L2_TREELEV=5
								MAX_L2_TREELEV=7
								if [ $USE_L2_ASIDE -eq 0 ]; then
									MIN_L2_TREELEV=0
									MAX_L2_TREELEV=0
								fi

								for (( L2_TREELEV = $MIN_L2_TREELEV ; L2_TREELEV <= $MAX_L2_TREELEV ; L2_TREELEV += 1 ))
								do
									# logging
									EXEC=bin/tot${TOTZEROS}${NTHREADS}sub${SUBZEROS}${TSUB}blk${NBLOCKS_SM}loop${LOOPZEROS}${NLOOP}ijp${WARPZEROS}${NWARP}ws${USE_WS}wr${USE_WR}l2a${USE_L2_ASIDE}lev${L2_TREELEV}
									echo "## generate $EXEC" >> $LOG

									# compile the N-body code w/ neighbor searching
									make gothic MEASURE_ELAPSED_TIME=1 HUNT_OPTIMAL_WALK_TREE=1 HUNT_OPTIMAL_INTEGRATE=1 HUNT_OPTIMAL_MAKE_TREE=0 HUNT_OPTIMAL_MAKE_NODE=0 HUNT_OPTIMAL_NEIGHBOUR=0 HUNT_OPTIMAL_SEPARATION=1 HUNT_OPTIMAL_L2LEVEL=1 NUM_NTHREADS=$NTHREADS NUM_TSUB=$TSUB NUM_NLOOP=${NLOOP} NUM_NWARP=${NWARP} NUM_BLOCKS_SM=${NBLOCKS_SM} USE_WARPSHUFFLE=${USE_WS} USE_WARPREDUCE=${USE_WR} USE_L2SETASIDE=${USE_L2_ASIDE} NUM_TREELEV_L2=${L2_TREELEV} ADOPT_GADGET_TYPE_MAC=1 1>>$LOG 2>>$ERR

									if [ -e bin/gothic ]; then
										# rename the executable
										mv bin/gothic $EXEC

										# generate job lists instead of running the execution file
										echo "${EXEC}" >> $LIST
									fi
									make clean
									# if [ -e $EXEC ]; then
									# 	TARGET=log/${EXEC##*/}
									# 	echo "srun -t ${TIME} numactl --cpunodebind=0 --localalloc ${EXEC} -absErr=${ABSERR} -file=${FILE} -jobID=$INDEX 1>>${TARGET}.o 2>>${TARGET}.e" >> $LIST
									# 	INDEX=`expr $INDEX + 1`
									# fi
								done
							done
						done
					done
				done
			done
		done
    done
done


exit 0
