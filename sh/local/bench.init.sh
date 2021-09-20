#!/bin/sh

# measure execution time of MAGI as a function of Ntot repeatedly
NITER=10

# for NTOT in 131072 262144 524288 1048576 2097152 4194304 8388608 16777216 33554432 67108864 134217728 268435456 536870912 1073741824
for NTOT in 131072 262144 524288 1048576 2097152 4194304 8388608 16777216 33554432 67108864 134217728 268435456 536870912
do
    for (( ii = 0; ii < $NITER; ii++ ))
    do
	/usr/bin/rm -rf dat doc
	make dir
	sh/local/init.bench.sh 28 $NTOT
	if [ $ii -eq 0 ]; then
	    ls -lSr dat/*.h5 > log/size.N${NTOT}.txt
	fi
    done
done

exit 0
