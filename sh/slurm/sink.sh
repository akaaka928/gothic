#!/bin/bash
###############################################################
#SBATCH -J gothic             # name of job
#SBATCH -t 24:00:00           # upper limit of elapsed time
# #SBATCH -t 00:10:00           # upper limit of elapsed time
#SBATCH -p normal             # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=1            # number of total MPI processes, set to SLURM_NTASKS (must be equal to number of GPUs)
#SBATCH --ntasks-per-socket=1 # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET (must be equal to number of GPUs per socket)
# #SBATCH --ntasks=2            # number of total MPI processes, set to SLURM_NTASKS (must be equal to number of GPUs)
# #SBATCH --ntasks-per-socket=2 # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET (must be equal to number of GPUs per socket)
#SBATCH --get-user-env        # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/gothic
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=100
fi
###############################################################
# topology of MPI processes
if [ -z "$NX" ]; then
    NX=2
fi
if [ -z "$NY" ]; then
    NY=1
fi
if [ -z "$NZ" ]; then
    NZ=1
fi
if [ $SLURM_NTASKS -eq 1 ]; then
    NX=1
    NY=1
    NZ=1
fi
PROCS=`expr $NX \* $NY \* $NZ`
if [ $PROCS -ne $SLURM_NTASKS ]; then
    echo "product of $NX, $NY, and $NZ must be equal to the number of total MPI processes ($SLURM_NTASKS)"
    exit 1
fi
###############################################################
# value of accuracy controling parameter: GADGET MAC by Springel (2005)
if [ -z "$ABSERR" ]; then
    # ABSERR=1.250000000e-1
    # ABSERR=6.250000000e-2
    # ABSERR=3.125000000e-2
    # ABSERR=1.562500000e-2
    # ABSERR=7.812500000e-3
    # ABSERR=3.906250000e-3
    ABSERR=1.953125000e-3
    # ABSERR=9.765625000e-4
    # ABSERR=4.882812500e-4
    # ABSERR=2.441406250e-4
    # ABSERR=1.220703125e-4
    # ABSERR=6.103515625e-5
fi
###############################################################
if [ -z "$REBUILD" ]; then
    REBUILD=16
fi
###############################################################
if [ -z "$BRENT" ]; then
    BRENT=1.0
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 100 ]; then
    FILE=m12iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 101 ]; then
    FILE=m12ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 102 ]; then
    FILE=m12ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 103 ]; then
    FILE=m12ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 104 ]; then
    FILE=m12ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 105 ]; then
    FILE=m12ra1_8_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 110 ]; then
    FILE=m13iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 111 ]; then
    FILE=m13ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 112 ]; then
    FILE=m13ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 113 ]; then
    FILE=m13ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 114 ]; then
    FILE=m13ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 120 ]; then
    FILE=m14iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 121 ]; then
    FILE=m14ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 122 ]; then
    FILE=m14ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 123 ]; then
    FILE=m14ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 124 ]; then
    FILE=m14ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 130 ]; then
    FILE=m15iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 131 ]; then
    FILE=m15ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 132 ]; then
    FILE=m15ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 133 ]; then
    FILE=m15ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 134 ]; then
    FILE=m15ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 140 ]; then
    FILE=m11iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 141 ]; then
    FILE=m11ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 142 ]; then
    FILE=m11ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 143 ]; then
    FILE=m11ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 144 ]; then
    FILE=m11ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 145 ]; then
    FILE=m11ra1_8_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 150 ]; then
    FILE=m10iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 151 ]; then
    FILE=m10ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 152 ]; then
    FILE=m10ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 153 ]; then
    FILE=m10ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 154 ]; then
    FILE=m10ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 155 ]; then
    FILE=m10ra1_8_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 160 ]; then
    FILE=m09iso_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 161 ]; then
    FILE=m09ra1_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 162 ]; then
    FILE=m09ra2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 163 ]; then
    FILE=m09ra1_2_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 164 ]; then
    FILE=m09ra1_4_run
    SINK=1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 165 ]; then
    FILE=m09ra1_8_run
    SINK=1
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -Nsink=$SINK -Nx=$NX -Ny=$NY -Nz=$NZ -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$SLURM_JOB_ID"
###############################################################


###############################################################
# job execution via SLURM
###############################################################
# set number of MPI processes per node
PROCS_PER_NODE=`expr $SLURM_NTASKS / $SLURM_JOB_NUM_NODES`
###############################################################
# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_NUM_NODES nodes"
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs per node"
TIME=`date`
echo "start: $TIME"
###############################################################
# execute the job
if [ $PROCS -gt 1 ]; then
    echo "mpiexec -n $SLURM_NTASKS sh/wrapper.sh $EXEC log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTION"
    mpiexec -n $SLURM_NTASKS sh/wrapper.sh $EXEC log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTION
else
    # set stdout and stderr
    STDOUT=log/${FILE}_$SLURM_JOB_NAME.o${SLURM_JOB_ID}
    STDERR=log/${FILE}_$SLURM_JOB_NAME.e${SLURM_JOB_ID}
    if [ `which numactl` ]; then
	# run with numactl
	echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
	numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
    else
	# run without numactl
	echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
	$EXEC $OPTION 1>>$STDOUT 2>>$STDERR
    fi
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
