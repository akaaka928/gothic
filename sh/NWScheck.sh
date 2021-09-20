#!/bin/sh

###############################################################
# generation of the sequences
if [ -z "$GEN" ]; then
    GEN=1
    # GEN=2
    # GEN=3
fi
###############################################################
# name of the series
if [ -z "$FILE" ]; then
    FILE=cusp
fi
###############################################################


h5dump -a info/time  dat/${FILE}-gen${GEN}-run*_score.h5 | grep -v DATA | grep -v ATTR | grep -v } > ${FILE}-gen${GEN}-time.log
h5dump -a info/score dat/${FILE}-gen${GEN}-run*_score.h5 | grep -v DATA | grep -v ATTR | grep -v } > ${FILE}-gen${GEN}-score.log

exit 0
