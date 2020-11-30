#!/bin/bash
#SBATCH -J gothic             # name of job
#SBATCH -t 24:00:00           # upper limit of elapsed time
#SBATCH -p amdrome            # partition name
#SBATCH -w amd2     # use compute node equips NVIDIA A100 PCIe
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --ntasks=1            # number of total MPI processes, set to SLURM_NTASKS (must be equal to number of GPUs)
#SBATCH --gpus-per-node=1
##SBATCH --gpu-freq=210 # requested GPU frequency: 210 MHz--1410 MHz (bin = 15 MHz)




# global configurations
EXEC=bin/gothic

# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=2
    # PROBLEM=10
    # PROBLEM=14
    # PROBLEM=26
    PROBLEM=27
    # PROBLEM=62
fi

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

if [ -z "$REBUILD" ]; then
    REBUILD=16
fi

if [ -z "$BRENT" ]; then
    BRENT=1.0
fi


###############################################################
# problem specific configurations
###############################################################
# cold collapse of a uniform sphere
if [ $PROBLEM -eq 0 ]; then
    FILE=ccuni
fi
###############################################################
# dynamical stability of a King sphere
if [ $PROBLEM -eq 1 ]; then
    FILE=king
fi
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
    FILE=hernquist
    REBUILD=8
    BRENT=1.0
fi
###############################################################
# dynamical stability of an NFW sphere small truncation radius
if [ $PROBLEM -eq 3 ]; then
    FILE=nfw
fi
###############################################################
# dynamical stability of an Einasto profile
if [ $PROBLEM -eq 4 ]; then
    FILE=einasto
fi
###############################################################
# dynamical stability of a Plummer profile
if [ $PROBLEM -eq 5 ]; then
    FILE=plummer
fi
###############################################################
# dynamical stability of a Burkert profile
if [ $PROBLEM -eq 6 ]; then
    FILE=burkert
fi
###############################################################
# dynamical stability of a Moore profile
if [ $PROBLEM -eq 7 ]; then
    FILE=moore
fi
###############################################################
# dynamical stability of a Two-power sphere
if [ $PROBLEM -eq 8 ]; then
    FILE=twopower
fi
###############################################################
# dynamical stability of a King sphere within an Einasto profile
if [ $PROBLEM -eq 10 ]; then
    FILE=hb
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere and a King sphere
if [ $PROBLEM -eq 11 ]; then
    FILE=hbd
fi
###############################################################
# dynamical stability of exponential disks (thick/thin disk) in an NFW sphere and a King sphere
if [ $PROBLEM -eq 12 ]; then
    FILE=hbdd
fi
###############################################################
# dynamical stability of thick exponential disk and thin Sersic disk in an Einasto sphere and a King sphere
if [ $PROBLEM -eq 13 ]; then
    FILE=ekes
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere
if [ $PROBLEM -eq 14 ]; then
    FILE=hd
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007)
if [ $PROBLEM -eq 20 ]; then
    FILE=m31
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007) in unit of GalactICS system
if [ $PROBLEM -eq 21 ]; then
    FILE=m31_gics
fi
###############################################################
# dynamical stability of multi components galaxy model
if [ $PROBLEM -eq 22 ]; then
    FILE=galaxy
fi
###############################################################
# dynamical stability of MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008)
if [ $PROBLEM -eq 23 ]; then
    FILE=mw
fi
###############################################################
# dynamical stability of M31 model (Sofue 2015; Gilbert et al. 2012)
if [ $PROBLEM -eq 24 ]; then
    FILE=s15
fi
###############################################################
# dynamical stability of multi components galaxy model with fixed number of particles
if [ $PROBLEM -eq 25 ]; then
    FILE=compare
fi
###############################################################
# dynamical stability of multi components galaxy model (only spherical components)
if [ $PROBLEM -eq 26 ]; then
    FILE=etg
fi
###############################################################
# dynamical stability of an M31 model (NFW halo, Hernquist bulge, and exponential disk)
# basically, this is Fardal et al. (2007) model
# stellar halo: Gilbert et al. (2012): \Sigma \propto R^-2.2; Rmin = 9kpc, Rmax = 176kpc; Ibata et al. (2014, ApJ, 780, 128): total stellar mass of the smooth halo is ~8e+9 Msun
# disk: Toomre's Q-value is set to reproduce Tenjes et al. (2017): Q_min = 1.8 @ 12-13 kpc
if [ $PROBLEM -eq 27 ]; then
    FILE=m31
fi
###############################################################
# dynamical stability of multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk)
if [ $PROBLEM -eq 28 ]; then
    FILE=ltg
fi
###############################################################
# dynamical stability of multi components galaxy model (Vasiliev & Athanassoula 2015) with fixed number of particles
if [ $PROBLEM -eq 29 ]; then
    FILE=va15
fi
###############################################################
# time evolution of MW/A defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 30 ]; then
    FILE=kd95a
fi
###############################################################
# time evolution of MW/B defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 31 ]; then
    FILE=kd95b
fi
###############################################################
# time evolution of MW/C defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 32 ]; then
    FILE=kd95c
fi
###############################################################
# time evolution of MW/D defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 33 ]; then
    FILE=kd95d
fi
###############################################################
# time evolution of M31/A defined in Widrow et al. (2003)
if [ $PROBLEM -eq 34 ]; then
    FILE=w03a
fi
###############################################################
# time evolution of M31/D defined in Widrow et al. (2003)
if [ $PROBLEM -eq 35 ]; then
    FILE=w03d
fi
###############################################################
# time evolution of MWa defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 36 ]; then
    FILE=mwa
fi
###############################################################
# time evolution of MWb defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 37 ]; then
    FILE=mwb
fi
###############################################################
# time evolution of MW like galaxy (based on Widrow et al. 2008, Bedorf et al. 2014)
if [ $PROBLEM -eq 38 ]; then
    FILE=bonsai
fi
###############################################################
# dynamical stability of a Plummer profile in table form
if [ $PROBLEM -eq 40 ]; then
    FILE=tplummer
fi
###############################################################
# dynamical stability of a Two-power slope profile in table form
if [ $PROBLEM -eq 41 ]; then
    FILE=dblpower
fi
###############################################################
# dynamical stability of a de Vaucouleurs sphere in table form
if [ $PROBLEM -eq 42 ]; then
    FILE=bulgetbl
fi
###############################################################
# dynamical stability of a spherical Sersic profile
if [ $PROBLEM -eq 50 ]; then
    FILE=deVaucouleurs
fi
###############################################################
# dynamical stability of a projected Two-power model
if [ $PROBLEM -eq 51 ]; then
    FILE=prjTwoPow
fi
###############################################################
# dynamical stability of a progenitor model for NW stream determined by Komiyama et al. (2018)
if [ $PROBLEM -eq 62 ]; then
    FILE=satellite
fi
###############################################################
# dynamical stability of multi components galaxy model (based on Cole & Binney 2017)
if [ $PROBLEM -eq 80 ]; then
    FILE=cb17
fi
###############################################################
# dynamical stability of multi components galaxy model (based on Cole & Binney 2017)
if [ $PROBLEM -eq 81 ]; then
    FILE=cb17_core
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 100 ]; then
    FILE=m12iso
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 101 ]; then
    FILE=m12ra1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 102 ]; then
    FILE=m12ra2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 103 ]; then
    FILE=m12ra1_2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 104 ]; then
    FILE=m12ra1_4
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 105 ]; then
    FILE=m12ra1_8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 110 ]; then
    FILE=m13iso
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 111 ]; then
    FILE=m13ra1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 112 ]; then
    FILE=m13ra2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 113 ]; then
    FILE=m13ra1_2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 114 ]; then
    FILE=m13ra1_4
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 120 ]; then
    FILE=m14iso
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 121 ]; then
    FILE=m14ra1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 122 ]; then
    FILE=m14ra2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 123 ]; then
    FILE=m14ra1_2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 124 ]; then
    FILE=m14ra1_4
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 130 ]; then
    FILE=m15iso
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 131 ]; then
    FILE=m15ra1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 132 ]; then
    FILE=m15ra2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 133 ]; then
    FILE=m15ra1_2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 134 ]; then
    FILE=m15ra1_4
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 140 ]; then
    FILE=m11iso
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 141 ]; then
    FILE=m11ra1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 142 ]; then
    FILE=m11ra2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 143 ]; then
    FILE=m11ra1_2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 144 ]; then
    FILE=m11ra1_4
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 145 ]; then
    FILE=m11ra1_8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 150 ]; then
    FILE=m10iso
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 151 ]; then
    FILE=m10ra1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 152 ]; then
    FILE=m10ra2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 153 ]; then
    FILE=m10ra1_2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 154 ]; then
    FILE=m10ra1_4
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 155 ]; then
    FILE=m10ra1_8
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 160 ]; then
    FILE=m09iso
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 161 ]; then
    FILE=m09ra1
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 162 ]; then
    FILE=m09ra2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 163 ]; then
    FILE=m09ra1_2
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 164 ]; then
    FILE=m09ra1_4
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 165 ]; then
    FILE=m09ra1_8
fi

# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -Nx=$NX -Ny=$NY -Nz=$NZ -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$SLURM_JOB_ID"


# job execution via SLURM
# set number of MPI processes per node
PROCS_PER_NODE=`expr $SLURM_NTASKS / $SLURM_JOB_NUM_NODES`

export MODULEPATH=$HOME/opt/modules:$MODULEPATH
module purge
module load cuda cub
module load openmpi
module load phdf5

# start logging
cd $SLURM_SUBMIT_DIR
echo "use $SLURM_JOB_NUM_NODES nodes"
echo "use $SLURM_JOB_CPUS_PER_NODE CPUs per node"
TIME=`date`
echo "start: $TIME"

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

# finish logging
TIME=`date`
echo "finish: $TIME"

exit 0
