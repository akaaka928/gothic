#!/bin/bash
###############################################################
#SBATCH -J accretion           # name of job
#SBATCH -t 24:00:00            # upper limit of elapsed time
#SBATCH -p normal              # partition name
#SBATCH --nodes=1              # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --get-user-env         # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/accretion
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    PROBLEM=100
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
    FILE=hernquist
    # CONFIG=single/hernquist.cfg
    CONFIG=anisotropy/hernquist_beta.cfg
    BH_MASS=0.002
fi
###############################################################
# dynamical stability of an NFW sphere with small truncation radius
if [ $PROBLEM -eq 3 ]; then
    FILE=nfw
    CONFIG=single/nfw.cfg
    # CONFIG=anisotropy/nfw_beta.cfg
    BH_MASS=0.002
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 100 ]; then
    FILE=m12iso
    CONFIG=losscone/m12iso.cfg
    BH_MASS=2.926550e+8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 101 ]; then
    FILE=m12ra1
    CONFIG=losscone/m12ra1.cfg
    BH_MASS=2.926550e+8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 102 ]; then
    FILE=m12ra2
    CONFIG=losscone/m12ra2.cfg
    BH_MASS=2.926550e+8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 103 ]; then
    FILE=m12ra3
    CONFIG=losscone/m12ra3.cfg
    BH_MASS=2.926550e+8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 104 ]; then
    FILE=m12ra4
    CONFIG=losscone/m12ra4.cfg
    BH_MASS=2.926550e+8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 105 ]; then
    FILE=m12ra5
    CONFIG=losscone/m12ra5.cfg
    BH_MASS=2.926550e+8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 106 ]; then
    FILE=m12ra3_4
    CONFIG=losscone/m12ra3_4.cfg
    BH_MASS=2.926550e+8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 107 ]; then
    FILE=m12ra1_2
    CONFIG=losscone/m12ra1_2.cfg
    BH_MASS=2.926550e+8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 110 ]; then
    FILE=m13iso
    CONFIG=losscone/m13iso.cfg
    BH_MASS=1.316911e+9
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 120 ]; then
    FILE=m14iso
    CONFIG=losscone/m14iso.cfg
    BH_MASS=3.703017e+9
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 130 ]; then
    FILE=m15iso
    CONFIG=losscone/m15iso.cfg
    BH_MASS=1.030934e+10
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 140 ]; then
    FILE=m11iso
    CONFIG=losscone/m11iso.cfg
    BH_MASS=1.525018e+6
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 150 ]; then
    FILE=m10iso
    CONFIG=losscone/m10iso.cfg
    BH_MASS=2.685220e+4
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 160 ]; then
    FILE=m09iso
    CONFIG=losscone/m09iso.cfg
    BH_MASS=4.712758e+2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 170 ]; then
    FILE=m08iso
    CONFIG=losscone/m08iso.cfg
    BH_MASS=8.298123
fi
###############################################################
# set input arguments
OPTION="-file=$FILE -config=$CONFIG -BH_mass=$BH_MASS"
###############################################################
# set environmental variables for OpenMP
OMP_OPT_ENV="env OMP_DISPLAY_ENV=verbose OMP_PLACES=cores"
# OMP_OPT_ENV="env KMP_AFFINITY=verbose,granularity=core,scatter"
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
    echo "$OMP_OPT_ENV numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    $OMP_OPT_ENV numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    # run without numactl
    echo "$OMP_OPT_ENV $EXEC $OPTION 1>>$STDOUT 2>>$STDERR"
    $OMP_OPT_ENV $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME"
###############################################################
