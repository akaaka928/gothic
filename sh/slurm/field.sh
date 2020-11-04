#!/bin/bash
###############################################################
#SBATCH -J gothic_ext         # name of job
#SBATCH -t 24:00:00           # upper limit of elapsed time
##SBATCH -t 00:10:00           # upper limit of elapsed time
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
    # PROBLEM=0
    # PROBLEM=11
    # PROBLEM=12
    # PROBLEM=13
    # PROBLEM=22
    # PROBLEM=30
    # PROBLEM=31
    # PROBLEM=32
    # PROBLEM=33
    # PROBLEM=34
    # PROBLEM=35
    # PROBLEM=36
    # PROBLEM=37
    PROBLEM=40
    # PROBLEM=41
    # PROBLEM=42
    # PROBLEM=43
    # PROBLEM=44
    # PROBLEM=45
    # PROBLEM=46
    # PROBLEM=47
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
# dynamical stability of a disk in a spherical potential field
if [ $PROBLEM -eq 0 ]; then
    FILE=disk
    SPHEPOT=va15
fi
###############################################################
# reproduction of Miki et al. (2016) in the disk coordinate system
if [ $PROBLEM -eq 11 ]; then
    FILE=m16king
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# reproduction of Kirihara et al. (2017) in the disk coordinate system
if [ $PROBLEM -eq 12 ]; then
    FILE=k17disk
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 13 ]; then
    FILE=k18nws
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# dynamical stability of an exponential disk in a spherical potential field (NFW sphere)
if [ $PROBLEM -eq 22 ]; then
    FILE=disk
    SPHEPOT=hd
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 30 ]; then
    FILE=subhalo0-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 31 ]; then
    FILE=subhalo1-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 32 ]; then
    FILE=subhalo2-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 33 ]; then
    FILE=subhalo3-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 34 ]; then
    FILE=subhalo4-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 35 ]; then
    FILE=subhalo5-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 36 ]; then
    FILE=subhalo6-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 37 ]; then
    FILE=subhalo7-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 40 ]; then
    FILE=subhalo10-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 41 ]; then
    FILE=subhalo11-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 42 ]; then
    FILE=subhalo12-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 43 ]; then
    FILE=subhalo13-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 44 ]; then
    FILE=subhalo14-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 45 ]; then
    FILE=subhalo15-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 46 ]; then
    FILE=subhalo16-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# throwing DM subhalo toward NW stream
if [ $PROBLEM -eq 47 ]; then
    FILE=subhalo17-backward
    SPHEPOT=m31
    DISKPOT=m31
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -pot_file_sphe=$SPHEPOT -pot_file_disk=$DISKPOT -Nx=$NX -Ny=$NY -Nz=$NZ -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$SLURM_JOB_ID"
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
