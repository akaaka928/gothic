#!/bin/bash
###############################################################
LOGDIR=bench
###############################################################
LIST0=${LOGDIR}/prec0.sh
LIST1=${LOGDIR}/prec1.sh
LIST2=${LOGDIR}/prec2.sh
LIST3=${LOGDIR}/prec3.sh
LIST4=${LOGDIR}/prec4.sh
LIST5=${LOGDIR}/prec5.sh
LIST6=${LOGDIR}/prec6.sh
LIST7=${LOGDIR}/prec7.sh
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
INDEX=0
PREID=0
PREZEROS=00
###############################################################
# for VALUE in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
for VALUE in 2.500000e-1 1.250000e-1 6.250000e-1 3.125000e-2 1.562500e-2 7.812500e-3 3.906250e-3 1.953125e-3 9.765625e-4
do
    ###########################################################
    # generate job lists instead of running the execution file
    ###########################################################
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
    ###########################################################
    if [ $INDEX != $PREID ]; then
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST0
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST1
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST2
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST3
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST4
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST5
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST6
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST7
    else
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST0
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST1
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST2
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST3
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST4
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST5
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST6
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} sh/hapacs/comq/measure.exec.sh\`" >> $LIST7
    fi
    ###########################################################
    PREZEROS=$IDXZEROS
    PREID=$INDEX
    INDEX=`expr $INDEX + 1`
    ###########################################################
done
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
