#!/bin/sh

for MASS in 70 75 80 85 90 95
do
    for HALO in p g i
    do
	for ORBIT in 00 01 02 03 04 05 06 07 10 11 12 13 14 15 16 17
	do
	    CFG=o${ORBIT}-m${MASS}${HALO}.cfg

	    echo "1" > $CFG
	    echo "2" >> $CFG
	    echo "k18nws nws/continue.param" >> $CFG
	    echo "m${MASS}${HALO} nws/subhalo${ORBIT}.param" >> $CFG
	done
    done
done

exit 0
