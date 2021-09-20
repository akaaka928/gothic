#!/bin/sh
###############################################################
if [ $# -ne 2 ]; then
    echo "$# inputs are detected while 2 inputs are required to specify <input value> and <digits>" 1>&2
    exit 1
fi
###############################################################
VALUE=$1
DIGIT=$2
###############################################################
INPUT="`echo ${#VALUE}`"
if [ "$INPUT" -le "$DIGIT" ]; then
    REM=`expr $DIGIT - $INPUT`
    ZEROS=""
    COUNT=0
    while [ "$COUNT" -lt "$REM" ]; do
	ZEROS="`echo ${ZEROS}0`"
	COUNT=`expr $COUNT + 1`
    done
fi
###############################################################
echo "${ZEROS}${VALUE}"
###############################################################
