#!/bin/sh
###############################################################
PROCS=16
###############################################################
EXEC=bin/extract
###############################################################
NCRIT=2048
NX3D=256
NY3D=256
NZ3D=256
NX=1024
NY=1024
NZ=1024
NV=1024
###############################################################
FINISH=23.5
INTERVAL=0.5
XMAX=3.0
YMAX=3.0
ZMAX=3.0
VMAX=3.0
XMIN=-$XMAX
YMIN=-$YMAX
ZMIN=-$ZMAX
VMIN=-$VMAX
###############################################################
START=0
END=`echo "scale=0; $FINISH / $INTERVAL" | bc`
INCREMENT=1
###############################################################
OPTION="-start=$START -end=$END -interval=$INCREMENT -ncrit=$NCRIT -nx=$NX -xmin=$XMIN -xmax=$XMAX -ny=$NY -ymin=$YMIN -ymax=$YMAX -nz=$NZ -zmin=$ZMIN -zmax=$ZMAX -nv=$NV -vmin=$VMIN -vmax=$VMAX -nx3D=$NX3D -ny3D=$NY3D -nz3D=$NZ3D"
###############################################################
FILE=a1b1g35
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
FILE=a1b1g34
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
FILE=a1b1g33
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
FILE=a1b1g32
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
FILE=a1b1g31
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
FILE=a1b1g30
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
FILE=a1b1g29
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
FILE=a1b1g28
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
FILE=a1b1g27
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
FILE=a1b1g26
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
FILE=a1b1g25
mpiexec -n $PROCS $EXEC -file=$FILE $OPTION 1>>log/${FILE}_anal.o 2>>log/${FILE}_anal.e
###############################################################
exit 0
###############################################################
