#!/bin/bash
###############################################################
LOGDIR=bench
LOG=${LOGDIR}/check.misc.log
ERR=${LOGDIR}/check.misc.err
###############################################################
echo "check on $HOSTNAME" > $LOG
echo "check on $HOSTNAME" > $ERR
TIME=`date`
echo "$TIME: compilation start" >> $LOG
echo "$TIME: compilation start" >> $ERR
###############################################################
LIST0=${LOGDIR}/check.misc0.sh
LIST1=${LOGDIR}/check.misc1.sh
LIST2=${LOGDIR}/check.misc2.sh
LIST3=${LOGDIR}/check.misc3.sh
LIST4=${LOGDIR}/check.misc4.sh
LIST5=${LOGDIR}/check.misc5.sh
LIST6=${LOGDIR}/check.misc6.sh
LIST7=${LOGDIR}/check.misc7.sh
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
echo "set -e" >> $LIST0
echo "set -e" >> $LIST1
echo "set -e" >> $LIST2
echo "set -e" >> $LIST3
echo "set -e" >> $LIST4
echo "set -e" >> $LIST5
echo "set -e" >> $LIST6
echo "set -e" >> $LIST7
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
echo "TARGET=nfw"       >> $LIST0
echo "TARGET=king"      >> $LIST1
echo "TARGET=hernquist" >> $LIST2
echo "TARGET=plummer"   >> $LIST3
echo "TARGET=hb"        >> $LIST4
echo "TARGET=hbd"       >> $LIST5
echo "TARGET=hbdd"      >> $LIST6
echo "TARGET=ekes"      >> $LIST7
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
# compile initial condition generators
make clean
# make init
# # clean up related object files
# rm -f bin/obj/*.o
# # generate initial conditions
# # sh/local/init.sh 0
# sh/local/init.sh 1
# sh/local/init.sh 2
# sh/local/init.sh 3
###############################################################
INDEX=0
PREID=0
PREZEROS=00
###############################################################
for NTHREADS in 128 256 512 1024
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
    # for TSUB in 32 16 8 4 2 1
    for TSUB in 32
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
	for SM_PREF in 1 0
	do
	    ###################################################
	    # logging
	    EXEC=bin/chk.tot${TOTZEROS}${NTHREADS}sub${SUBZEROS}${TSUB}sm${SM_PREF}
	    echo "## generate $EXEC" >> $LOG
	    echo "## generate $EXEC" >> $ERR
	    ###################################################
	    # compile the N-body code w/ neighbor searching
	    make gothic CHECK_FUNCTION_CALLS=1 LOCALIZE_I_PARTICLES=1 BRUTE_FORCE_LOCALIZER=1 FACILE_NEIGHBOR_SEARCH=1 CHECK_CALL_MAKE_TREE=1 CHECK_CALL_MAKE_NODE=1 CHECK_CALL_NEIGHBOUR=1 NUM_NTHREADS=$NTHREADS NUM_TSUB=$TSUB USE_WARPSHUFFLE=0 PREF_SHARED_MEM=$SM_PREF ADOPT_GADGET_TYPE_MAC=0 ADOPT_WS93_TYPE_MAC=1 1>>$LOG 2>>$ERR
	    ###################################################
	    # rename the executable
	    mv bin/gothic $EXEC
	    ###################################################
	    # clean up related object files
	    rm -f bin/obj/*.o
	    ###################################################
	    # generate job lists instead of running the execution file
	    if [ -e $EXEC ]; then
		###############################################
		# pad 0s for INDEX
		digit=2
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
		if [ $INDEX != $PREID ]; then
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST0
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST1
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST2
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST3
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST4
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST5
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST6
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST7
		else
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST0
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST1
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST2
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST3
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST4
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST5
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST6
		    echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},EXEC=${EXEC} sh/hapacs/comq/bench.exec.sh\`" >> $LIST7
		fi
		###############################################
		PREZEROS=$IDXZEROS
		PREID=$INDEX
		INDEX=`expr $INDEX + 1`
		###############################################
	    fi
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
