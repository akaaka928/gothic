#!/bin/bash
###############################################################
#SBATCH -J m31obs              # name of job
#SBATCH -t 02:00:00            # upper limit of elapsed time
#SBATCH -p normal              # partition name
#SBATCH --nodes=1              # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=16             # number of total MPI processes, set to SLURM_NTASKS
#SBATCH --ntasks-per-socket=8  # number of MPI processes per socket, set to SLURM_NTASKS_PER_SOCKET
#SBATCH --get-user-env         # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/m31obs
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=111
    # PROBLEM=112
    # PROBLEM=113
    # PROBLEM=114
    PROBLEM=115
    # PROBLEM=116
fi
###############################################################
# set number of N-body particles per bin to estimate density
if [ -z "$NCRIT" ]; then
    NCRIT=2048
fi
###############################################################
# set number of grid points for density maps
if [ -z "$NX3D" ]; then
    NX3D=256
fi
if [ -z "$NY3D" ]; then
    NY3D=256
fi
if [ -z "$NZ3D" ]; then
    NZ3D=256
fi
if [ -z "$NX" ]; then
    NX=1024
fi
if [ -z "$NY" ]; then
    NY=1024
fi
if [ -z "$NZ" ]; then
    NZ=1024
fi
if [ -z "$NV" ]; then
    NV=1024
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# reproduction of Miki et al. (2016) in the disk coordinate system
if [ $PROBLEM -eq 111 ]; then
    FILE=m16king
    FINISH=1246.875
    INTERVAL=3.125
    XIMIN=-5.0
    XIMAX=7.0
    ETAMIN=-8.5
    ETAMAX=3.5
    DMIN=-80.0
    DMAX=130.0
    VMIN=-700.0
    VMAX=-100.0
fi
###############################################################
# reproduction of Kirihara et al. (2017) in the disk coordinate system
if [ $PROBLEM -eq 112 ]; then
    FILE=k17disk
    FINISH=1246.875
    INTERVAL=3.125
    XIMIN=-5.0
    XIMAX=7.0
    ETAMIN=-8.5
    ETAMAX=3.5
    DMIN=-80.0
    DMAX=130.0
    VMIN=-700.0
    VMAX=-100.0
fi
###############################################################
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 113 ]; then
    # FILE=nws-continue
    FILE=nws-test-m9_5-orbit0
    FINISH=6600.0
    # FILE=k18nws
    # FINISH=14000.0
    INTERVAL=25.0
    XIMIN=-9.0
    XIMAX=7.0
    ETAMIN=-3.0
    ETAMAX=13.0
    DMIN=-50.0
    DMAX=200.0
    VMIN=-540.0
    VMAX=-200.0
fi
###############################################################
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 114 ]; then
    FILE=nws-continue
    # FILE=nws-test-m9_5-orbit0
    # FILE=nws-test-m9_5-orbit1
    FINISH=6600.0
    INTERVAL=25.0
    # XIMIN=-11.0
    # XIMAX=21.0
    XIMIN=-16.0
    XIMAX=26.0
    # ETAMIN=-14.0
    # ETAMAX=18.0
    ETAMIN=-20.0
    ETAMAX=22.0
    DMIN=-150.0
    DMAX=350.0
    VMIN=-650.0
    VMAX=-50.0

    NX3D=256
    NY3D=256
    NZ3D=256
    NX=2048
    NY=2048
    NZ=2048
    NV=2048
fi
###############################################################
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 115 ]; then
    FILE=nws-continue
    # FILE=nws-test-m9_5-orbit2
    # FILE=nws-test-m9_5-orbit3
    FINISH=10700.0
    INTERVAL=25.0
    XIMIN=-11.0
    XIMAX=21.0
    ETAMIN=-14.0
    ETAMAX=18.0
    DMIN=-150.0
    DMAX=350.0
    VMIN=-650.0
    VMAX=-50.0

    NX3D=256
    NY3D=256
    NZ3D=256
    NX=2048
    NY=2048
    NZ=2048
    NV=2048
fi
###############################################################
# reproduction of Komiyama et al. (2018) in the disk coordinate system (show whole domain)
if [ $PROBLEM -eq 116 ]; then
    FILE=k18nws
    FINISH=14000.0
    INTERVAL=25.0
    XIMIN=-22.0
    XIMAX=22.0
    ETAMIN=-22.0
    ETAMAX=22.0
    DMIN=-150.0
    DMAX=350.0
    VMIN=-650.0
    VMAX=-50.0

    NX3D=256
    NY3D=256
    NZ3D=256
    NX=2048
    NY=2048
    NZ=2048
    NV=2048
fi
###############################################################
# count up number of snapshot files
if [ -z "$START" ]; then
    START=0
fi
if [ -z "$END" ]; then
    END=`echo "scale=0; $FINISH / $INTERVAL" | bc`
fi
if [ -z "$INCREMENT" ]; then
    INCREMENT=1
fi
# START=325
# END=325
###############################################################
# set input arguments
OPTION="-file=$FILE -start=$START -end=$END -interval=$INCREMENT -ncrit=$NCRIT -nxi=$NX -ximin=$XIMIN -ximax=$XIMAX -neta=$NY -etamin=$ETAMIN -etamax=$ETAMAX -nD=$NZ -Dmin=$DMIN -Dmax=$DMAX -nvlos=$NV -vlosmin=$VMIN -vlosmax=$VMAX -nx3D=$NX3D -ny3D=$NY3D -nz3D=$NZ3D"
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
echo "mpiexec -n $SLURM_NTASKS sh/wrapper.sh $EXEC log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTION"
mpiexec -n $SLURM_NTASKS sh/wrapper.sh $EXEC log/${FILE}_${SLURM_JOB_NAME} $SLURM_JOB_ID $PROCS_PER_NODE $SLURM_NTASKS_PER_SOCKET $OPTION
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
exit 0
###############################################################
