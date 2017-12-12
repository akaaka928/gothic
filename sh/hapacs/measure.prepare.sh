#!/bin/bash
###############################################################
LOGDIR=bench
###############################################################
LIST00=${LOGDIR}/time00.sh
LIST01=${LOGDIR}/time01.sh
LIST02=${LOGDIR}/time02.sh
LIST03=${LOGDIR}/time03.sh
LIST04=${LOGDIR}/time04.sh
LIST05=${LOGDIR}/time05.sh
LIST06=${LOGDIR}/time06.sh
LIST07=${LOGDIR}/time07.sh
LIST08=${LOGDIR}/time08.sh
LIST09=${LOGDIR}/time09.sh
LIST10=${LOGDIR}/time10.sh
LIST11=${LOGDIR}/time11.sh
LIST12=${LOGDIR}/time12.sh
LIST13=${LOGDIR}/time13.sh
LIST14=${LOGDIR}/time14.sh
LIST15=${LOGDIR}/time15.sh
###############################################################
echo "#!/bin/sh" > $LIST00
echo "#!/bin/sh" > $LIST01
echo "#!/bin/sh" > $LIST02
echo "#!/bin/sh" > $LIST03
echo "#!/bin/sh" > $LIST04
echo "#!/bin/sh" > $LIST05
echo "#!/bin/sh" > $LIST06
echo "#!/bin/sh" > $LIST07
echo "#!/bin/sh" > $LIST08
echo "#!/bin/sh" > $LIST09
echo "#!/bin/sh" > $LIST10
echo "#!/bin/sh" > $LIST11
echo "#!/bin/sh" > $LIST12
echo "#!/bin/sh" > $LIST13
echo "#!/bin/sh" > $LIST14
echo "#!/bin/sh" > $LIST15
###############################################################
echo "###############################################################" >> $LIST00
echo "###############################################################" >> $LIST01
echo "###############################################################" >> $LIST02
echo "###############################################################" >> $LIST03
echo "###############################################################" >> $LIST04
echo "###############################################################" >> $LIST05
echo "###############################################################" >> $LIST06
echo "###############################################################" >> $LIST07
echo "###############################################################" >> $LIST08
echo "###############################################################" >> $LIST09
echo "###############################################################" >> $LIST10
echo "###############################################################" >> $LIST11
echo "###############################################################" >> $LIST12
echo "###############################################################" >> $LIST13
echo "###############################################################" >> $LIST14
echo "###############################################################" >> $LIST15
###############################################################
echo "set -e" >> $LIST00
echo "set -e" >> $LIST01
echo "set -e" >> $LIST02
echo "set -e" >> $LIST03
echo "set -e" >> $LIST04
echo "set -e" >> $LIST05
echo "set -e" >> $LIST06
echo "set -e" >> $LIST07
echo "set -e" >> $LIST08
echo "set -e" >> $LIST09
echo "set -e" >> $LIST10
echo "set -e" >> $LIST11
echo "set -e" >> $LIST12
echo "set -e" >> $LIST13
echo "set -e" >> $LIST14
echo "set -e" >> $LIST15
###############################################################
echo "###############################################################" >> $LIST00
echo "###############################################################" >> $LIST01
echo "###############################################################" >> $LIST02
echo "###############################################################" >> $LIST03
echo "###############################################################" >> $LIST04
echo "###############################################################" >> $LIST05
echo "###############################################################" >> $LIST06
echo "###############################################################" >> $LIST07
echo "###############################################################" >> $LIST08
echo "###############################################################" >> $LIST09
echo "###############################################################" >> $LIST10
echo "###############################################################" >> $LIST11
echo "###############################################################" >> $LIST12
echo "###############################################################" >> $LIST13
echo "###############################################################" >> $LIST14
echo "###############################################################" >> $LIST15
###############################################################
echo "TARGET=king"      >> $LIST00
echo "TARGET=hernquist" >> $LIST01
echo "TARGET=nfw"       >> $LIST02
echo "TARGET=einasto"   >> $LIST03
echo "TARGET=plummer"   >> $LIST04
echo "TARGET=burkert"   >> $LIST05
echo "TARGET=moore"     >> $LIST06
echo "TARGET=twopower"  >> $LIST07
echo "TARGET=hb"        >> $LIST08
echo "TARGET=hbd"       >> $LIST09
echo "TARGET=hbdd"      >> $LIST10
echo "TARGET=ekes"      >> $LIST11
echo "TARGET=m31"       >> $LIST12
echo "TARGET=galaxy"    >> $LIST13
echo "TARGET=mw"        >> $LIST14
echo "TARGET=s15"       >> $LIST15
###############################################################
echo "###############################################################" >> $LIST00
echo "###############################################################" >> $LIST01
echo "###############################################################" >> $LIST02
echo "###############################################################" >> $LIST03
echo "###############################################################" >> $LIST04
echo "###############################################################" >> $LIST05
echo "###############################################################" >> $LIST06
echo "###############################################################" >> $LIST07
echo "###############################################################" >> $LIST08
echo "###############################################################" >> $LIST09
echo "###############################################################" >> $LIST10
echo "###############################################################" >> $LIST11
echo "###############################################################" >> $LIST12
echo "###############################################################" >> $LIST13
echo "###############################################################" >> $LIST14
echo "###############################################################" >> $LIST15
###############################################################
INDEX=0
PREID=0
PREZEROS=00
###############################################################
# for VALUE in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
for VALUE in 2.500000e-1 1.250000e-1 6.250000e-2 3.125000e-2 1.562500e-2 7.812500e-3 3.906250e-3 1.953125e-3 9.765625e-4 4.8828125e-4 2.44140625e-4 1.220703125e-4 6.103515625e-5 3.0517578125e-5 1.52587890625e-5 7.62939453125e-6 3.814697265625e-6 1.9073486328125e-6
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
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST00
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST01
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST02
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST03
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST04
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST05
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST06
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST07
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST08
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST09
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST10
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST11
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST12
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST13
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST14
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -W depend=afterany:\$jobid${PREZEROS}${PREID} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST15
    else
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST00
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST01
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST02
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST03
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST04
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST05
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST06
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST07
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST08
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST09
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST10
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST11
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST12
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST13
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST14
	echo "jobid${IDXZEROS}${INDEX}=\`echo sleep 10 | qsub -N \${TARGET}${IDXZEROS}${INDEX} -v FILE=\${TARGET},VALUE=${VALUE} ./measure.exec.sh\`" >> $LIST15
    fi
    ###########################################################
    PREZEROS=$IDXZEROS
    PREID=$INDEX
    INDEX=`expr $INDEX + 1`
    ###########################################################
done
###############################################################
echo "###############################################################" >> $LIST00
echo "###############################################################" >> $LIST01
echo "###############################################################" >> $LIST02
echo "###############################################################" >> $LIST03
echo "###############################################################" >> $LIST04
echo "###############################################################" >> $LIST05
echo "###############################################################" >> $LIST06
echo "###############################################################" >> $LIST07
echo "###############################################################" >> $LIST08
echo "###############################################################" >> $LIST09
echo "###############################################################" >> $LIST10
echo "###############################################################" >> $LIST11
echo "###############################################################" >> $LIST12
echo "###############################################################" >> $LIST13
echo "###############################################################" >> $LIST14
echo "###############################################################" >> $LIST15
###############################################################
echo "exit 0" >> $LIST00
echo "exit 0" >> $LIST01
echo "exit 0" >> $LIST02
echo "exit 0" >> $LIST03
echo "exit 0" >> $LIST04
echo "exit 0" >> $LIST05
echo "exit 0" >> $LIST06
echo "exit 0" >> $LIST07
echo "exit 0" >> $LIST08
echo "exit 0" >> $LIST09
echo "exit 0" >> $LIST10
echo "exit 0" >> $LIST11
echo "exit 0" >> $LIST12
echo "exit 0" >> $LIST13
echo "exit 0" >> $LIST14
echo "exit 0" >> $LIST15
###############################################################
echo "###############################################################" >> $LIST00
echo "###############################################################" >> $LIST01
echo "###############################################################" >> $LIST02
echo "###############################################################" >> $LIST03
echo "###############################################################" >> $LIST04
echo "###############################################################" >> $LIST05
echo "###############################################################" >> $LIST06
echo "###############################################################" >> $LIST07
echo "###############################################################" >> $LIST08
echo "###############################################################" >> $LIST09
echo "###############################################################" >> $LIST10
echo "###############################################################" >> $LIST11
echo "###############################################################" >> $LIST12
echo "###############################################################" >> $LIST13
echo "###############################################################" >> $LIST14
echo "###############################################################" >> $LIST15
###############################################################
chmod +x $LIST00
chmod +x $LIST01
chmod +x $LIST02
chmod +x $LIST03
chmod +x $LIST04
chmod +x $LIST05
chmod +x $LIST06
chmod +x $LIST07
chmod +x $LIST08
chmod +x $LIST09
chmod +x $LIST10
chmod +x $LIST11
chmod +x $LIST12
chmod +x $LIST13
chmod +x $LIST14
chmod +x $LIST15
###############################################################
exit 0
###############################################################
