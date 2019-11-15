#!/bin/bash
#PBS -A GALAXY
#PBS -q gpu
#PBS -N magi
#PBS -b 1
#PBS -l elapstim_req=00:10:00
#PBS -v OMP_NUM_THREADS=24
#PBS -T mvapich
#PBS -v NQSV_MPI_VER=2.3.1/intel-cuda10.1
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/magi
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # M31 model
    PROBLEM=27

    # progenitor model by Komiyama et al. (2018)
    # PROBLEM=62

    # sub-halo using c-M relation by Prada et al. (2012)
    # PROBLEM=100
    # PROBLEM=101
    # PROBLEM=102
    # PROBLEM=103
    # PROBLEM=104
    # PROBLEM=105

    # sub-halo using c-M relation by Ishiyama & Ando (2019)
    # PROBLEM=110
    # PROBLEM=111
    # PROBLEM=112
    # PROBLEM=113
    # PROBLEM=114
    # PROBLEM=115

    # sub-halo using c-M relation by Gilman et al. (2019)
    # PROBLEM=120
    # PROBLEM=121
    # PROBLEM=122
    # PROBLEM=123
    # PROBLEM=124
    # PROBLEM=125
fi
###############################################################
# number of N-body particles
if [ -z "$NTOT" ]; then
    NTOT=1048576
fi
###############################################################
# dump file generation interval (in units of minute)
if [ -z "$SAVE" ]; then
    SAVE=55.0
fi
###############################################################
# denoise distribution function
if [ -z "$DENOISE_DF" ]; then
    DENOISE_DF=1
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# dynamical stability of an M31 model (NFW halo, Hernquist bulge, and exponential disk)
# basically, this is Fardal et al. (2007) model
# stellar halo: Gilbert et al. (2012): \Sigma \propto R^-2.2; Rmin = 9kpc, Rmax = 176kpc; Ibata et al. (2014, ApJ, 780, 128): total stellar mass of the smooth halo is ~8e+9 Msun
# disk: Toomre's Q-value is set to reproduce Tenjes et al. (2017): Q_min = 1.8 @ 12-13 kpc
if [ $PROBLEM -eq 27 ]; then
    FILE=m31
    CONFIG=galaxy/m31.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1175.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a progenitor model for NW stream determined by Komiyama et al. (2018)
if [ $PROBLEM -eq 62 ]; then
    FILE=sat
    CONFIG=nws/k18.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a sub-halo model (10^7 Msun, c-M relation by Prada et al. 2012)
if [ $PROBLEM -eq 100 ]; then
    NTOT=209715
    FILE=prada_mass70
    CONFIG=nws/subhalo/mass70-prada.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a sub-halo model (10^7.5 Msun, c-M relation by Prada et al. 2012)
if [ $PROBLEM -eq 101 ]; then
    NTOT=663178
    FILE=prada_mass75
    CONFIG=nws/subhalo/mass75-prada.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^8 Msun, c-M relation by Prada et al. 2012)
if [ $PROBLEM -eq 102 ]; then
    NTOT=2097152
    FILE=prada_mass80
    CONFIG=nws/subhalo/mass80-prada.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^8.5 Msun, c-M relation by Prada et al. 2012)
if [ $PROBLEM -eq 103 ]; then
    NTOT=6631777
    FILE=prada_mass85
    CONFIG=nws/subhalo/mass85-prada.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^9 Msun, c-M relation by Prada et al. 2012)
if [ $PROBLEM -eq 104 ]; then
    NTOT=20971520
    FILE=prada_mass90
    CONFIG=nws/subhalo/mass90-prada.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^9.5 Msun, c-M relation by Prada et al. 2012)
if [ $PROBLEM -eq 105 ]; then
    NTOT=66317769
    FILE=prada_mass95
    CONFIG=nws/subhalo/mass95-prada.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a sub-halo model (10^7 Msun, c-M relation by Ishiyama & Ando 2019)
if [ $PROBLEM -eq 110 ]; then
    NTOT=209715
    FILE=ishiyama_mass70
    CONFIG=nws/subhalo/mass70-ishiyama.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a sub-halo model (10^7.5 Msun, c-M relation by Ishiyama & Ando 2019)
if [ $PROBLEM -eq 111 ]; then
    NTOT=663178
    FILE=ishiyama_mass75
    CONFIG=nws/subhalo/mass75-ishiyama.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^8 Msun, c-M relation by Ishiyama & Ando 2019)
if [ $PROBLEM -eq 112 ]; then
    NTOT=2097152
    FILE=ishiyama_mass80
    CONFIG=nws/subhalo/mass80-ishiyama.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^8.5 Msun, c-M relation by Ishiyama & Ando 2019)
if [ $PROBLEM -eq 113 ]; then
    NTOT=6631777
    FILE=ishiyama_mass85
    CONFIG=nws/subhalo/mass85-ishiyama.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^9 Msun, c-M relation by Ishiyama & Ando 2019)
if [ $PROBLEM -eq 114 ]; then
    NTOT=20971520
    FILE=ishiyama_mass90
    CONFIG=nws/subhalo/mass90-ishiyama.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^9.5 Msun, c-M relation by Ishiyama & Ando 2019)
if [ $PROBLEM -eq 115 ]; then
    NTOT=66317769
    FILE=ishiyama_mass95
    CONFIG=nws/subhalo/mass95-ishiyama.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a sub-halo model (10^7 Msun, c-M relation by Gilman et al. 2019)
if [ $PROBLEM -eq 120 ]; then
    NTOT=209715
    FILE=gilman_mass70
    CONFIG=nws/subhalo/mass70-gilman.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a sub-halo model (10^7.5 Msun, c-M relation by Gilman et al. 2019)
if [ $PROBLEM -eq 121 ]; then
    NTOT=663178
    FILE=gilman_mass75
    CONFIG=nws/subhalo/mass75-gilman.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^8 Msun, c-M relation by Gilman et al. 2019)
if [ $PROBLEM -eq 122 ]; then
    NTOT=2097152
    FILE=gilman_mass80
    CONFIG=nws/subhalo/mass80-gilman.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^8.5 Msun, c-M relation by Gilman et al. 2019)
if [ $PROBLEM -eq 123 ]; then
    NTOT=6631777
    FILE=gilman_mass85
    CONFIG=nws/subhalo/mass85-gilman.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^9 Msun, c-M relation by Gilman et al. 2019)
if [ $PROBLEM -eq 124 ]; then
    NTOT=20971520
    FILE=gilman_mass90
    CONFIG=nws/subhalo/mass90-gilman.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a subhalo model (10^9.5 Msun, c-M relation by Gilman et al. 2019)
if [ $PROBLEM -eq 125 ]; then
    NTOT=66317769
    FILE=gilman_mass95
    CONFIG=nws/subhalo/mass95-gilman.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=25.0
fi
###############################################################
# set input arguments
OPTION="-file=$FILE -config=$CONFIG -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE -denoisingDistributionFunction=$DENOISE_DF"
###############################################################
# set environmental variables for OpenMP
# OMP_OPT_ENV="env OMP_DISPLAY_ENV=verbose OMP_PLACES=cores"
OMP_OPT_ENV="env KMP_AFFINITY=verbose,granularity=core,scatter"
###############################################################


###############################################################
# job execution via NQSV
###############################################################
# set stdout and stderr
STDOUT=log/${FILE}_$PBS_JOBNAME.o${PBS_JOBID}
STDERR=log/${FILE}_$PBS_JOBNAME.e${PBS_JOBID}
###############################################################
# start logging
cd $PBS_O_WORKDIR
TIME=`date`
echo "start: $TIME"
###############################################################
module purge
export MODULEPATH=$MODULEPATH:/work/CSPP/ymiki/opt/module
module load mvapich/2.3.1/intel-cuda10.1
module load gsl lis
module load phdf5
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
exit 0
###############################################################
