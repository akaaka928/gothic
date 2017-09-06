#!/bin/sh

NUM_NODE=1
JOB_NAME=gothic
JOB_TIME=00:30:00

# global configurations
EXEC=bin/gothic

# PROBLEM=20
PROBLEM=2

NX=2
NY=1
NZ=1
PROCS=`expr $NX \* $NY \* $NZ`

# ABSERR=1.250000000e-1
# ABSERR=6.250000000e-2
# ABSERR=3.125000000e-2
# ABSERR=1.562500000e-2
# ABSERR=7.812500000e-3
ABSERR=3.906250000e-3
# ABSERR=1.953125000e-3
# ABSERR=9.765625000e-4
# ABSERR=4.882812500e-4
# ABSERR=2.441406250e-4
# ABSERR=1.220703125e-4


# problem specific configurations
if [ $PROBLEM -eq 0 ]; then
    FILE=ccuni
fi
if [ $PROBLEM -eq 1 ]; then
    FILE=king
fi
if [ $PROBLEM -eq 2 ]; then
    FILE=hernquist
fi
if [ $PROBLEM -eq 3 ]; then
    FILE=nfw
fi
if [ $PROBLEM -eq 4 ]; then
    FILE=einasto
fi
if [ $PROBLEM -eq 5 ]; then
    FILE=plummer
fi
if [ $PROBLEM -eq 6 ]; then
    FILE=burkert
fi
if [ $PROBLEM -eq 7 ]; then
    FILE=moore
fi
if [ $PROBLEM -eq 8 ]; then
    FILE=twopower
fi
if [ $PROBLEM -eq 10 ]; then
    FILE=hb
fi
if [ $PROBLEM -eq 11 ]; then
    FILE=hbd
fi
if [ $PROBLEM -eq 12 ]; then
    FILE=hbdd
fi
if [ $PROBLEM -eq 13 ]; then
    FILE=ekes
fi
if [ $PROBLEM -eq 20 ]; then
    FILE=m31
fi
if [ $PROBLEM -eq 21 ]; then
    FILE=m31_gics
fi
if [ $PROBLEM -eq 22 ]; then
    FILE=galaxy
fi
if [ $PROBLEM -eq 23 ]; then
    FILE=mw
fi
if [ $PROBLEM -eq 24 ]; then
    FILE=s15
fi
if [ $PROBLEM -eq 25 ]; then
    FILE=compare
fi
if [ $PROBLEM -eq 26 ]; then
    FILE=spherical
fi
if [ $PROBLEM -eq 27 ]; then
    FILE=m31_mod
fi
if [ $PROBLEM -eq 28 ]; then
    FILE=ltg
fi
if [ $PROBLEM -eq 29 ]; then
    FILE=va15
fi
if [ $PROBLEM -eq 30 ]; then
    FILE=kd95a
fi
if [ $PROBLEM -eq 31 ]; then
    FILE=kd95b
fi
if [ $PROBLEM -eq 32 ]; then
    FILE=kd95c
fi
if [ $PROBLEM -eq 33 ]; then
    FILE=kd95d
fi
if [ $PROBLEM -eq 34 ]; then
    FILE=w03a
fi
if [ $PROBLEM -eq 35 ]; then
    FILE=w03d
fi
if [ $PROBLEM -eq 36 ]; then
    FILE=mwa
fi
if [ $PROBLEM -eq 37 ]; then
    FILE=mwb
fi
if [ $PROBLEM -eq 38 ]; then
    FILE=bonsai
fi
if [ $PROBLEM -eq 40 ]; then
    FILE=tplummer
fi
if [ $PROBLEM -eq 41 ]; then
    FILE=dblpower
fi
if [ $PROBLEM -eq 42 ]; then
    FILE=bulgetbl
fi
if [ $PROBLEM -eq 50 ]; then
    FILE=deVaucouleurs
fi
if [ $PROBLEM -eq 51 ]; then
    FILE=prjTwoPow
fi
if [ $PROBLEM -eq 80 ]; then
    FILE=cb17
fi
if [ $PROBLEM -eq 81 ]; then
    FILE=cb17_core
fi

# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -Nx=$NX -Ny=$NY -Nz=$NZ -jobID=$SLURM_JOB_ID"

export MV2_USE_CUDA=1

module purge
module load cuda
module load gcc mvapich2-gdr

export MODULEPATH=$MODULEPATH:/home/ymiki/opt/Modules
module load phdf5
module load cub

STDOUT=log/$JOB_NAME.$SLURM_JOB_ID.out
STDERR=log/$JOB_NAME.$SLURM_JOB_ID.err

echo "srun --job-name=$JOB_NAME --time=$JOB_TIME --partition=comq --nodes=$NUM_NODE --ntasks=$PROCS --ntasks-per-socket=1 $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
srun --job-name=$JOB_NAME --time=$JOB_TIME --partition=comq --nodes=$NUM_NODE --ntasks=$PROCS --ntasks-per-socket=1 $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
