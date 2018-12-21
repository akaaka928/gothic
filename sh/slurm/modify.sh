#!/bin/bash
###############################################################
#SBATCH -J editor             # name of job
#SBATCH -t 00:30:00           # upper limit of elapsed time
#SBATCH -p normal             # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --get-user-env        # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/editor
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=0
    # PROBLEM=1
    # PROBLEM=2
    # PROBLEM=10
    PROBLEM=11
    # PROBLEM=12
    # PROBLEM=13
    # PROBLEM=21
    # PROBLEM=22
    # PROBLEM=23
    # PROBLEM=24
fi
###############################################################
# dump file generation interval (in units of minute)
if [ -z "$SAVE" ]; then
    SAVE=55.0
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# dynamical stability of a disk in a spherical potential field
if [ $PROBLEM -eq 0 ]; then
    FILE=disk
    CFG=external/split.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=11.75
    INTERVAL=0.25
fi
###############################################################
# dynamical stability of two disks in a spherical potential field
if [ $PROBLEM -eq 1 ]; then
    FILE=ltg_disk
    CFG=external/ltg_split.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1175.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of an exponential disk in a spherical potential field (NFW sphere)
if [ $PROBLEM -eq 2 ]; then
    FILE=disk
    CFG=external/hd_split.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=104.0
    INTERVAL=8.0
fi
###############################################################
# dynamical stability of M31 in the observed coordinate
if [ $PROBLEM -eq 10 ]; then
    FILE=m31_obs
    CFG=gss/obs.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=3175.0
    INTERVAL=25.0
fi
###############################################################
# reproduction of Miki et al. (2016) in the disk coordinate system
if [ $PROBLEM -eq 11 ]; then
    FILE=m16king
    CFG=gss/satellite.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=12.5
    # FINISH=100.0
    # FINISH=1050.0
    # INTERVAL=1.5625
    FINISH=1246.875
    INTERVAL=3.125
fi
###############################################################
# reproduction of Kirihara et al. (2017) in the disk coordinate system
if [ $PROBLEM -eq 12 ]; then
    FILE=k17disk
    CFG=gss/satellite.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=12.5
    # FINISH=100.0
    # FINISH=1050.0
    # INTERVAL=1.5625
    FINISH=1246.875
    INTERVAL=3.125
fi
###############################################################
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 13 ]; then
    FILE=nws
    CFG=nws/satellite.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# GSS simulation in observed frame with live M31
if [ $PROBLEM -eq 20 ]; then
    FILE=gss
    CFG=gss/obs_run.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1246.875
    INTERVAL=3.125
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 30 ]; then
    FILE=halocusp_run
    CFG=fornax/halocusp_run.cfg
    EPS=9.765625000e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 31 ]; then
    FILE=halocore1_run
    CFG=fornax/halocore1_run.cfg
    EPS=9.765625000e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 32 ]; then
    FILE=halocore2_run
    CFG=fornax/halocore2_run.cfg
    EPS=9.765625000e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 33 ]; then
    FILE=halocore3_run
    CFG=fornax/halocore3_run.cfg
    EPS=9.765625000e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 100 ]; then
    FILE=m12iso_run
    CFG=losscone/m12iso_run.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 101 ]; then
    FILE=m12ra1_run
    CFG=losscone/m12ra1_run.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 102 ]; then
    FILE=m12ra2_run
    CFG=losscone/m12ra2_run.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 103 ]; then
    FILE=m12ra1_2_run
    CFG=losscone/m12ra1_2_run.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 104 ]; then
    FILE=m12ra1_4_run
    CFG=losscone/m12ra1_4_run.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 105 ]; then
    FILE=m12ra1_8_run
    CFG=losscone/m12ra1_8_run.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 110 ]; then
    FILE=m13iso_run
    CFG=losscone/m13iso_run.cfg
    EPS=2.5e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 111 ]; then
    FILE=m13ra1_run
    CFG=losscone/m13ra1_run.cfg
    EPS=2.5e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 112 ]; then
    FILE=m13ra2_run
    CFG=losscone/m13ra2_run.cfg
    EPS=2.5e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 113 ]; then
    FILE=m13ra1_2_run
    CFG=losscone/m13ra1_2_run.cfg
    EPS=2.5e-1
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 114 ]; then
    FILE=m13ra1_4_run
    CFG=losscone/m13ra1_4_run.cfg
    EPS=2.5e-1
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 120 ]; then
    FILE=m14iso_run
    CFG=losscone/m14iso_run.cfg
    EPS=5.0e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 121 ]; then
    FILE=m14ra1_run
    CFG=losscone/m14ra1_run.cfg
    EPS=5.0e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 122 ]; then
    FILE=m14ra2_run
    CFG=losscone/m14ra2_run.cfg
    EPS=5.0e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 123 ]; then
    FILE=m14ra1_2_run
    CFG=losscone/m14ra1_2_run.cfg
    EPS=5.0e-1
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 124 ]; then
    FILE=m14ra1_4_run
    CFG=losscone/m14ra1_4_run.cfg
    EPS=5.0e-1
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 130 ]; then
    FILE=m15iso_run
    CFG=losscone/m15iso_run.cfg
    EPS=1.0
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 131 ]; then
    FILE=m15ra1_run
    CFG=losscone/m15ra1_run.cfg
    EPS=1.0
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 132 ]; then
    FILE=m15ra2_run
    CFG=losscone/m15ra2_run.cfg
    EPS=1.0
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 133 ]; then
    FILE=m15ra1_2_run
    CFG=losscone/m15ra1_2_run.cfg
    EPS=1.0
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 134 ]; then
    FILE=m15ra1_4_run
    CFG=losscone/m15ra1_4_run.cfg
    EPS=1.0
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 140 ]; then
    FILE=m11iso_run
    CFG=losscone/m11iso_run.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 141 ]; then
    FILE=m11ra1_run
    CFG=losscone/m11ra1_run.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 142 ]; then
    FILE=m11ra2_run
    CFG=losscone/m11ra2_run.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 143 ]; then
    FILE=m11ra1_2_run
    CFG=losscone/m11ra1_2_run.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 144 ]; then
    FILE=m11ra1_4_run
    CFG=losscone/m11ra1_4_run.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 145 ]; then
    FILE=m11ra1_8_run
    CFG=losscone/m11ra1_8_run.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 150 ]; then
    FILE=m10iso_run
    CFG=losscone/m10iso_run.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 151 ]; then
    FILE=m10ra1_run
    CFG=losscone/m10ra1_run.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 152 ]; then
    FILE=m10ra2_run
    CFG=losscone/m10ra2_run.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 153 ]; then
    FILE=m10ra1_2_run
    CFG=losscone/m10ra1_2_run.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 154 ]; then
    FILE=m10ra1_4_run
    CFG=losscone/m10ra1_4_run.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 155 ]; then
    FILE=m10ra1_8_run
    CFG=losscone/m10ra1_8_run.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 160 ]; then
    FILE=m09iso_run
    CFG=losscone/m09iso_run.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 161 ]; then
    FILE=m09ra1_run
    CFG=losscone/m09ra1_run.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 162 ]; then
    FILE=m09ra2_run
    CFG=losscone/m09ra2_run.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 163 ]; then
    FILE=m09ra1_2_run
    CFG=losscone/m09ra1_2_run.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 164 ]; then
    FILE=m09ra1_4_run
    CFG=losscone/m09ra1_4_run.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 165 ]; then
    FILE=m09ra1_8_run
    CFG=losscone/m09ra1_8_run.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# set input arguments
OPTION="-file=$FILE -list=$CFG -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE"
###############################################################


###############################################################
# job execution via SLURM
###############################################################
# set stdout and stderr
STDOUT=log/${FILE}_$SLURM_JOB_NAME.o${SLURM_JOB_ID}
STDERR=log/${FILE}_$SLURM_JOB_NAME.e${SLURM_JOB_ID}
###############################################################
# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs"
TIME=`date`
echo "start: $TIME"
###############################################################
# execute the job
if [ `which numactl` ]; then
    # run with numactl
    echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    # run without numactl
    echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
