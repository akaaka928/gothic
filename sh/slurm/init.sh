#!/bin/bash
#SBATCH -J magi               # name of job
#SBATCH -t 02:00:00           # upper limit of elapsed time
#SBATCH -p normal             # partition name
#SBATCH --nodes=1             # number of nodes, set to SLURM_JOB_NUM_NODES
#SBATCH --get-user-env        # retrieve the login environment variables
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/magi
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=2
    # PROBLEM=3
    # PROBLEM=10
    # PROBLEM=14
    # PROBLEM=15
    # PROBLEM=21
    # PROBLEM=23
    # PROBLEM=26
    # PROBLEM=27
    # PROBLEM=33
    # PROBLEM=39
    # PROBLEM=46
    # PROBLEM=60
    # PROBLEM=61
    # PROBLEM=62
    # PROBLEM=63
    # PROBLEM=64
    # PROBLEM=90
    # PROBLEM=91
    # PROBLEM=92
    # PROBLEM=93
    # PROBLEM=94
    # PROBLEM=95
    # PROBLEM=96
    # PROBLEM=97
    # PROBLEM=98
    # PROBLEM=99
    # PROBLEM=100
    PROBLEM=1000
fi
###############################################################
# number of N-body particles
if [ -z "$NTOT" ]; then
    # NTOT=128
    # NTOT=256
    # NTOT=512
    # NTOT=1024
    # NTOT=2048
    # NTOT=4096
    # NTOT=8192
    # NTOT=16384
    # NTOT=32768
    # NTOT=65536
    # NTOT=131072
    # NTOT=262144
    # NTOT=524288
    # NTOT=1048576
    # NTOT=1900000
    # NTOT=2097152
    # NTOT=4194304
    # NTOT=8388608
    NTOT=16777216
    # NTOT=33554432
    # NTOT=67108864
    # NTOT=134217728
    # NTOT=268435456
    # NTOT=536870912
    # NTOT=1073741824
    # NTOT=2147483648
    # NTOT=4294967296
    # NTOT=8589934592
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
# cold collapse of a uniform sphere
if [ $PROBLEM -eq 0 ]; then
    EXEC=bin/uniformsphere
    FILE=ccuni
    UNIT=-1
    MTOT=1.0
    VIRIAL=0.2
    RAD=1.0
    EPS=1.5625e-2
    ETA=0.5
    FINISH=11.75
    INTERVAL=0.25
fi
###############################################################
# dynamical stability of a King sphere
if [ $PROBLEM -eq 1 ]; then
    FILE=king
    CONFIG=single/king.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Hernquist sphere
if [ $PROBLEM -eq 2 ]; then
    FILE=hernquist
    # CONFIG=single/hernquist.cfg
    CONFIG=anisotropy/hernquist_beta.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of an NFW sphere with small truncation radius
if [ $PROBLEM -eq 3 ]; then
    FILE=nfw
    CONFIG=single/nfw.cfg
    # CONFIG=anisotropy/nfw_beta.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of an Einasto profile
if [ $PROBLEM -eq 4 ]; then
    FILE=einasto
    CONFIG=single/einasto.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Plummer profile
if [ $PROBLEM -eq 5 ]; then
    FILE=plummer
    CONFIG=single/plummer.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Burkert profile
if [ $PROBLEM -eq 6 ]; then
    FILE=burkert
    CONFIG=single/burkert.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Moore profile
if [ $PROBLEM -eq 7 ]; then
    FILE=moore
    CONFIG=single/moore.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a Two-power sphere
if [ $PROBLEM -eq 8 ]; then
    FILE=twopower
    CONFIG=single/twopower.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a King sphere within an Einasto profile
if [ $PROBLEM -eq 10 ]; then
    FILE=hb
    # CONFIG=multi/hb.cfg
    CONFIG=anisotropy/hb_beta.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere and a King sphere
if [ $PROBLEM -eq 11 ]; then
    FILE=hbd
    CONFIG=multi/hbd.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=47.5
    INTERVAL=0.5
    # FINISH=1648.0
    # INTERVAL=16.0
fi
###############################################################
# dynamical stability of exponential disks (thick/thin disk) in an NFW sphere and a King sphere
if [ $PROBLEM -eq 12 ]; then
    FILE=hbdd
    CONFIG=multi/hbdd.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=47.5
    INTERVAL=0.5
    # FINISH=1400.0
    # INTERVAL=8.0
fi
###############################################################
# dynamical stability of thick exponential disk and thin Sersic disk in an Einasto sphere and a King sphere
if [ $PROBLEM -eq 13 ]; then
    FILE=ekes
    CONFIG=multi/ekes.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=47.5
    INTERVAL=0.5
    # FINISH=1400.0
    # INTERVAL=8.0
fi
###############################################################
# dynamical stability of an exponential disk in an NFW sphere
if [ $PROBLEM -eq 14 ]; then
    FILE=hd
    CONFIG=multi/hd.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=104.0
    INTERVAL=8.0
    NTOT=600034
fi
###############################################################
# dynamical stability of a Hernquist bulge with a central MBH
if [ $PROBLEM -eq 15 ]; then
    FILE=hernquist_bh
    CONFIG=anisotropy/hernquist_beta_bh.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007)
if [ $PROBLEM -eq 20 ]; then
    FILE=m31
    CONFIG=galaxy/f07.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=1175.0
    INTERVAL=25.0
    # FINISH=5175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of an M31 model determined by Fardal et al. (2007) in unit of GalactICS system
if [ $PROBLEM -eq 21 ]; then
    FILE=m31_gics
    CONFIG=galaxy/f07_gics.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=5175.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model
if [ $PROBLEM -eq 22 ]; then
    FILE=galaxy
    CONFIG=galaxy/galaxy.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=75.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008)
if [ $PROBLEM -eq 23 ]; then
    FILE=mw
    CONFIG=galaxy/mw.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=75.0
    INTERVAL=25.0
    FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of M31 model (Sofue 2015; Gilbert et al. 2012)
if [ $PROBLEM -eq 24 ]; then
    FILE=s15
    CONFIG=galaxy/s15.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=75.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model with fixed number of particles
if [ $PROBLEM -eq 25 ]; then
    if [ $NTOT -lt 4194304 ]; then
	NTOT=4194304
    fi
    FILE=compare
    CONFIG=galaxy/compare.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=3175.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (only spherical components)
if [ $PROBLEM -eq 26 ]; then
    FILE=etg
    CONFIG=galaxy/etg.cfg
    EPS=1.5625e-2
    ETA=0.5
    INTERVAL=25.0
    # FINISH=75.0
    FINISH=375.0
    # FINISH=1175.0
    # FINISH=3175.0
fi
###############################################################
# dynamical stability of an M31 model (NFW halo, Hernquist bulge, and exponential disk)
# basically, this is Fardal et al. (2007) model
# stellar halo: Gilbert et al. (2012): \Sigma \propto R^-2.2; Rmin = 9kpc, Rmax = 176kpc; Ibata et al. (2014, ApJ, 780, 128): total stellar mass of the smooth halo is ~8e+9 Msun
# disk: Toomre's Q-value is set to reproduce Tenjes et al. (2017): Q_min = 1.8 @ 12-13 kpc
if [ $PROBLEM -eq 27 ]; then
    FILE=m31
    CONFIG=galaxy/m31.cfg
    EPS=1.5625e-2
    # EPS=7.8125e-3
    ETA=0.5
    # FINISH=75.0
    # FINISH=375.0
    FINISH=1175.0
    # FINISH=3175.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk)
if [ $PROBLEM -eq 28 ]; then
    FILE=ltg
    CONFIG=galaxy/ltg.cfg
    EPS=1.5625e-2
    ETA=0.5
    # ETA=0.25
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=1175.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (Vasiliev & Athanassoula 2015) with fixed number of particles
if [ $PROBLEM -eq 29 ]; then
    if [ $NTOT -ne 1000000 ]; then
	NTOT=1000000
    fi
    FILE=va15
    CONFIG=galaxy/va15.cfg
    EPS=0.01
    ETA=1.0
    FINISH=100.0
    INTERVAL=4.0
fi
###############################################################
# time evolution of MW/A defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 30 ]; then
    FILE=kd95a
    CONFIG=gics/kd95/MWa.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of MW/B defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 31 ]; then
    FILE=kd95b
    CONFIG=gics/kd95/MWb.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of MW/C defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 32 ]; then
    FILE=kd95c
    CONFIG=gics/kd95/MWc.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of MW/D defined in Kuijken & Dubinski (1995)
if [ $PROBLEM -eq 33 ]; then
    FILE=kd95d
    CONFIG=gics/kd95/MWd.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of M31/A defined in Widrow et al. (2003)
if [ $PROBLEM -eq 34 ]; then
    FILE=w03a
    CONFIG=gics/w03/M31a.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of M31/D defined in Widrow et al. (2003)
if [ $PROBLEM -eq 35 ]; then
    FILE=w03d
    CONFIG=gics/w03/M31d.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=124.0
    # INTERVAL=4.0
    FINISH=2032.0
    INTERVAL=16.0
fi
###############################################################
# time evolution of MWa defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 36 ]; then
    FILE=mwa
    CONFIG=gics/wd05/MWa.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=124.0
    INTERVAL=4.0
    # FINISH=1016.0
    # INTERVAL=8.0
fi
###############################################################
# time evolution of MWb defined in Widrow & Dubinski (2005)
if [ $PROBLEM -eq 37 ]; then
    FILE=mwb
    CONFIG=gics/wd05/MWb.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=124.0
    INTERVAL=4.0
    # FINISH=1016.0
    # INTERVAL=8.0
fi
###############################################################
# time evolution of MW like galaxy (based on Widrow et al. 2008; Bedorf et al. 2014)
if [ $PROBLEM -eq 38 ]; then
    FILE=bonsai
    CONFIG=gics/bonsai/b14.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=124.0
    INTERVAL=4.0
    # FINISH=1176.0
    # INTERVAL=8.0
fi
###############################################################
# dynamical stability of multi components galaxy model
if [ $PROBLEM -eq 39 ]; then
    if [ $NTOT -ne 200000 ]; then
	NTOT=200000
    fi
    FILE=diskhalo
    CONFIG=debug/galaxy.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a Plummer profile in table form
if [ $PROBLEM -eq 40 ]; then
    FILE=tplummer
    CONFIG=table/tplummer.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a Two-power slope profile in table form
if [ $PROBLEM -eq 41 ]; then
    FILE=dblpower
    CONFIG=table/dblpower.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a de Vaucouleurs sphere in table form
if [ $PROBLEM -eq 42 ]; then
    FILE=bulgetbl
    CONFIG=table/bulgetbl.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a Hernquist profile in table form
if [ $PROBLEM -eq 43 ]; then
    FILE=thernquist
    CONFIG=table/thernquist.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a Plummer profile in analytic form
if [ $PROBLEM -eq 44 ]; then
    FILE=aplummer
    CONFIG=table/aplummer.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a Hernquist profile in analytic form
if [ $PROBLEM -eq 45 ]; then
    FILE=ahernquist
    CONFIG=table/ahernquist.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a disk in table form within NFW halo
if [ $PROBLEM -eq 46 ]; then
    FILE=thd
    CONFIG=table/thd.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of a spherical Sersic profile
if [ $PROBLEM -eq 50 ]; then
    FILE=deVaucouleurs
    CONFIG=abel/deVaucouleurs.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a projected Two-power model
if [ $PROBLEM -eq 51 ]; then
    FILE=prjTwoPow
    CONFIG=abel/prjTwoPow.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1575.0
    INTERVAL=25.0
fi
###############################################################
# dynamical stability of a progenitor model for GSS determined by Miki et al. (2016)
if [ $PROBLEM -eq 60 ]; then
    FILE=sat
    CONFIG=gss/m16.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=1575.0
    INTERVAL=25.0
    # FINISH=5175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of a progenitor model for GSS determined by Kirihara et al. (2017)
if [ $PROBLEM -eq 61 ]; then
    FILE=sat
    CONFIG=gss/k17.cfg
    EPS=1.5625e-2
    ETA=0.5
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=1575.0
    INTERVAL=25.0
    # FINISH=5175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of a progenitor model for NW stream determined by Komiyama et al. (2018)
if [ $PROBLEM -eq 62 ]; then
    FILE=sat
    CONFIG=nws/k18.cfg
    EPS=1.5625e-2
    ETA=0.5
    # INTERVAL=25.0
    # FINISH=75.0
    # FINISH=1575.0
    # FINISH=5175.0
    # FINISH=5175.0
    INTERVAL=100.0
    FINISH=14000.0
fi
###############################################################
# dynamical stability of a progenitor model for NW stream
if [ $PROBLEM -eq 63 ]; then
    FILE=sat
    CONFIG=nws/core.cfg
    EPS=1.5625e-2
    ETA=0.5
    INTERVAL=100.0
    FINISH=14000.0
fi
###############################################################
# dynamical stability of a progenitor model for NW stream
if [ $PROBLEM -eq 64 ]; then
    FILE=sat
    CONFIG=nws/cusp.cfg
    EPS=1.5625e-2
    ETA=0.5
    INTERVAL=100.0
    FINISH=14000.0
fi
###############################################################
# dynamical stability of Eridanus II (Crnojevic et al. 2016; Li et al. 2017)
if [ $PROBLEM -eq 70 ]; then
    FILE=eri2
    CONFIG=dwarf/EriII.cfg
    EPS=2.0
    ETA=0.5
    FINISH=1175.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of a King sphere with a central cusp
if [ $PROBLEM -eq 71 ]; then
    FILE=king
    CONFIG=single/king_slope.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=23.5
    INTERVAL=0.5
fi
###############################################################
# dynamical stability of multi components galaxy model (based on Cole & Binney 2017)
if [ $PROBLEM -eq 80 ]; then
    FILE=cb17
    CONFIG=galaxy/cb17.cfg
    EPS=1.5625e-2
    ETA=0.5
    # ETA=0.25
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=1175.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# dynamical stability of multi components galaxy model (based on Cole & Binney 2017)
if [ $PROBLEM -eq 81 ]; then
    FILE=cb17_core
    CONFIG=galaxy/cb17_core.cfg
    EPS=1.5625e-2
    ETA=0.5
    # ETA=0.25
    # FINISH=75.0
    # INTERVAL=25.0
    FINISH=1175.0
    INTERVAL=25.0
    # FINISH=3175.0
    # INTERVAL=25.0
fi
###############################################################
# Fornax simulation (globular cluster model 0)
if [ $PROBLEM -eq 90 ]; then
    FILE=gc0
    CONFIG=fornax/gc0.cfg
    NTOT=100000
    EPS=9.765625000e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax 1 (globular cluster)
if [ $PROBLEM -eq 91 ]; then
    FILE=for1
    CONFIG=fornax/gc1.cfg
    NTOT=100000
    EPS=9.765625000e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax 2 (globular cluster)
if [ $PROBLEM -eq 92 ]; then
    FILE=for2
    CONFIG=fornax/gc2.cfg
    NTOT=100000
    EPS=4.882812500e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax 3 (globular cluster)
if [ $PROBLEM -eq 93 ]; then
    FILE=for3
    CONFIG=fornax/gc3.cfg
    NTOT=100000
    EPS=1.220703125e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax 4 (globular cluster)
if [ $PROBLEM -eq 94 ]; then
    FILE=for4
    CONFIG=fornax/gc4.cfg
    NTOT=100000
    EPS=6.103515625e-5
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax 5 (globular cluster)
if [ $PROBLEM -eq 95 ]; then
    FILE=for5
    CONFIG=fornax/gc5.cfg
    NTOT=100000
    EPS=2.441406250e-4
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 96 ]; then
    FILE=halocusp
    CONFIG=fornax/halocusp.cfg
    NTOT=11000000
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 97 ]; then
    FILE=halocore1
    CONFIG=fornax/halocore1.cfg
    NTOT=11000000
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 98 ]; then
    FILE=halocore2
    CONFIG=fornax/halocore2.cfg
    NTOT=11000000
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Fornax simulation
if [ $PROBLEM -eq 99 ]; then
    FILE=halocore3
    CONFIG=fornax/halocore3.cfg
    NTOT=11000000
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 100 ]; then
    FILE=m12iso
    CONFIG=losscone/m12iso.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 101 ]; then
    FILE=m12ra1
    CONFIG=losscone/m12ra1.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 102 ]; then
    FILE=m12ra2
    CONFIG=losscone/m12ra2.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 103 ]; then
    FILE=m12ra1_2
    CONFIG=losscone/m12ra1_2.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 104 ]; then
    FILE=m12ra1_4
    CONFIG=losscone/m12ra1_4.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 105 ]; then
    FILE=m12ra1_8
    CONFIG=losscone/m12ra1_8.cfg
    EPS=6.25e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 110 ]; then
    FILE=m13iso
    CONFIG=losscone/m13iso.cfg
    EPS=2.5e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 111 ]; then
    FILE=m13ra1
    CONFIG=losscone/m13ra1.cfg
    EPS=2.5e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 112 ]; then
    FILE=m13ra2
    CONFIG=losscone/m13ra2.cfg
    EPS=2.5e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 113 ]; then
    FILE=m13ra1_2
    CONFIG=losscone/m13ra1_2.cfg
    EPS=2.5e-1
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 114 ]; then
    FILE=m13ra1_4
    CONFIG=losscone/m13ra1_4.cfg
    EPS=2.5e-1
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# # DM accretion onto central MBH
# if [ $PROBLEM -eq 115 ]; then
#     FILE=m13ra1_8
#     CONFIG=losscone/m13ra1_8.cfg
#     EPS=7.8125e-3
#     ETA=0.5
#     FINISH=1000.0
#     INTERVAL=25.0
# fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 120 ]; then
    FILE=m14iso
    CONFIG=losscone/m14iso.cfg
    EPS=5.0e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 121 ]; then
    FILE=m14ra1
    CONFIG=losscone/m14ra1.cfg
    EPS=5.0e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 122 ]; then
    FILE=m14ra2
    CONFIG=losscone/m14ra2.cfg
    EPS=5.0e-1
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 123 ]; then
    FILE=m14ra1_2
    CONFIG=losscone/m14ra1_2.cfg
    EPS=5.0e-1
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 124 ]; then
    FILE=m14ra1_4
    CONFIG=losscone/m14ra1_4.cfg
    EPS=5.0e-1
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# # DM accretion onto central MBH
# if [ $PROBLEM -eq 125 ]; then
#     FILE=m14ra1_8
#     CONFIG=losscone/m14ra1_8.cfg
#     EPS=7.8125e-3
#     ETA=0.5
#     FINISH=1000.0
#     INTERVAL=25.0
# fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 130 ]; then
    FILE=m15iso
    CONFIG=losscone/m15iso.cfg
    EPS=1.0
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 131 ]; then
    FILE=m15ra1
    CONFIG=losscone/m15ra1.cfg
    EPS=1.0
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 132 ]; then
    FILE=m15ra2
    CONFIG=losscone/m15ra2.cfg
    EPS=1.0
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 133 ]; then
    FILE=m15ra1_2
    CONFIG=losscone/m15ra1_2.cfg
    EPS=1.0
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 134 ]; then
    FILE=m15ra1_4
    CONFIG=losscone/m15ra1_4.cfg
    EPS=1.0
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# # DM accretion onto central MBH
# if [ $PROBLEM -eq 135 ]; then
#     FILE=m15ra1_8
#     CONFIG=losscone/m15ra1_8.cfg
#     EPS=7.8125e-3
#     ETA=0.5
#     FINISH=1000.0
#     INTERVAL=25.0
# fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 140 ]; then
    FILE=m11iso
    CONFIG=losscone/m11iso.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 141 ]; then
    FILE=m11ra1
    CONFIG=losscone/m11ra1.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 142 ]; then
    FILE=m11ra2
    CONFIG=losscone/m11ra2.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 143 ]; then
    FILE=m11ra1_2
    CONFIG=losscone/m11ra1_2.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 144 ]; then
    FILE=m11ra1_4
    CONFIG=losscone/m11ra1_4.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 145 ]; then
    FILE=m11ra1_8
    CONFIG=losscone/m11ra1_8.cfg
    EPS=3.125e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 150 ]; then
    FILE=m10iso
    CONFIG=losscone/m10iso.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 151 ]; then
    FILE=m10ra1
    CONFIG=losscone/m10ra1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 152 ]; then
    FILE=m10ra2
    CONFIG=losscone/m10ra2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 153 ]; then
    FILE=m10ra1_2
    CONFIG=losscone/m10ra1_2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 154 ]; then
    FILE=m10ra1_4
    CONFIG=losscone/m10ra1_4.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 155 ]; then
    FILE=m10ra1_8
    CONFIG=losscone/m10ra1_8.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 160 ]; then
    FILE=m09iso
    CONFIG=losscone/m09iso.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 161 ]; then
    FILE=m09ra1
    CONFIG=losscone/m09ra1.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 162 ]; then
    FILE=m09ra2
    CONFIG=losscone/m09ra2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=1000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 163 ]; then
    FILE=m09ra1_2
    CONFIG=losscone/m09ra1_2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 164 ]; then
    FILE=m09ra1_4
    CONFIG=losscone/m09ra1_4.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# DM accretion onto central MBH
if [ $PROBLEM -eq 165 ]; then
    FILE=m09ra1_8
    CONFIG=losscone/m09ra1_8.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=2000.0
    INTERVAL=25.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 200 ]; then
    NTOT=1000000
    FILE=w3iso
    CONFIG=roi/king/w3iso.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 201 ]; then
    NTOT=1000000
    FILE=w3ra2
    CONFIG=roi/king/w3ra2.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 202 ]; then
    NTOT=1000000
    FILE=w3ra1
    CONFIG=roi/king/w3ra1.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 203 ]; then
    NTOT=1000000
    FILE=w3ra1_2
    CONFIG=roi/king/w3ra1_2.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 204 ]; then
    NTOT=1000000
    FILE=w3ra1_4
    CONFIG=roi/king/w3ra1_4.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 205 ]; then
    NTOT=1000000
    FILE=w3ra1_8
    CONFIG=roi/king/w3ra1_8.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 210 ]; then
    NTOT=1000000
    FILE=w5iso
    CONFIG=roi/king/w5iso.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 211 ]; then
    NTOT=1000000
    FILE=w5ra2
    CONFIG=roi/king/w5ra2.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 212 ]; then
    NTOT=1000000
    FILE=w5ra1
    CONFIG=roi/king/w5ra1.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 213 ]; then
    NTOT=1000000
    FILE=w5ra1_2
    CONFIG=roi/king/w5ra1_2.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 214 ]; then
    NTOT=1000000
    FILE=w5ra1_4
    CONFIG=roi/king/w5ra1_4.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 215 ]; then
    NTOT=1000000
    FILE=w5ra1_8
    CONFIG=roi/king/w5ra1_8.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 220 ]; then
    NTOT=1000000
    FILE=w7iso
    CONFIG=roi/king/w7iso.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 221 ]; then
    NTOT=1000000
    FILE=w7ra2
    CONFIG=roi/king/w7ra2.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 222 ]; then
    NTOT=1000000
    FILE=w7ra1
    CONFIG=roi/king/w7ra1.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 223 ]; then
    NTOT=1000000
    FILE=w7ra1_2
    CONFIG=roi/king/w7ra1_2.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 224 ]; then
    NTOT=1000000
    FILE=w7ra1_4
    CONFIG=roi/king/w7ra1_4.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 225 ]; then
    NTOT=1000000
    FILE=w7ra1_8
    CONFIG=roi/king/w7ra1_8.cfg
    EPS=4.8828125e-4
    ETA=0.25
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 230 ]; then
    FILE=m12iso
    CONFIG=roi/nfw/m12iso.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 231 ]; then
    FILE=m12ra2
    CONFIG=roi/nfw/m12ra2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 232 ]; then
    FILE=m12ra1
    CONFIG=roi/nfw/m12ra1.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 233 ]; then
    FILE=m12ra1_2
    CONFIG=roi/nfw/m12ra1_2.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 234 ]; then
    FILE=m12ra1_4
    CONFIG=roi/nfw/m12ra1_4.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 235 ]; then
    FILE=m12ra1_8
    CONFIG=roi/nfw/m12ra1_8.cfg
    EPS=1.5625e-2
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 240 ]; then
    FILE=m09iso
    CONFIG=roi/nfw/m09iso.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 241 ]; then
    FILE=m09ra2
    CONFIG=roi/nfw/m09ra2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 242 ]; then
    FILE=m09ra1
    CONFIG=roi/nfw/m09ra1.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 243 ]; then
    FILE=m09ra1_2
    CONFIG=roi/nfw/m09ra1_2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 244 ]; then
    FILE=m09ra1_4
    CONFIG=roi/nfw/m09ra1_4.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 245 ]; then
    FILE=m09ra1_8
    CONFIG=roi/nfw/m09ra1_8.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 250 ]; then
    FILE=a00iso
    CONFIG=roi/dehnen/a00iso.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 251 ]; then
    FILE=a00ra2
    CONFIG=roi/dehnen/a00ra2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 252 ]; then
    FILE=a00ra1
    CONFIG=roi/dehnen/a00ra1.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 253 ]; then
    FILE=a00ra1_2
    CONFIG=roi/dehnen/a00ra1_2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 254 ]; then
    FILE=a00ra1_4
    CONFIG=roi/dehnen/a00ra1_4.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 260 ]; then
    FILE=a05iso
    CONFIG=roi/dehnen/a05iso.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 261 ]; then
    FILE=a05ra2
    CONFIG=roi/dehnen/a05ra2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 262 ]; then
    FILE=a05ra1
    CONFIG=roi/dehnen/a05ra1.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 263 ]; then
    FILE=a05ra1_2
    CONFIG=roi/dehnen/a05ra1_2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 264 ]; then
    FILE=a05ra1_4
    CONFIG=roi/dehnen/a05ra1_4.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 265 ]; then
    FILE=a05ra1_8
    CONFIG=roi/dehnen/a05ra1_8.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 270 ]; then
    FILE=a10iso
    CONFIG=roi/dehnen/a10iso.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 271 ]; then
    FILE=a10ra2
    CONFIG=roi/dehnen/a10ra2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 272 ]; then
    FILE=a10ra1
    CONFIG=roi/dehnen/a10ra1.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 273 ]; then
    FILE=a10ra1_2
    CONFIG=roi/dehnen/a10ra1_2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 274 ]; then
    FILE=a10ra1_4
    CONFIG=roi/dehnen/a10ra1_4.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 275 ]; then
    FILE=a10ra1_8
    CONFIG=roi/dehnen/a10ra1_8.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 280 ]; then
    FILE=a15iso
    CONFIG=roi/dehnen/a15iso.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 281 ]; then
    FILE=a15ra2
    CONFIG=roi/dehnen/a15ra2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 282 ]; then
    FILE=a15ra1
    CONFIG=roi/dehnen/a15ra1.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 283 ]; then
    FILE=a15ra1_2
    CONFIG=roi/dehnen/a15ra1_2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 284 ]; then
    FILE=a15ra1_4
    CONFIG=roi/dehnen/a15ra1_4.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 285 ]; then
    FILE=a15ra1_8
    CONFIG=roi/dehnen/a15ra1_8.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 290 ]; then
    FILE=a20iso
    CONFIG=roi/dehnen/a20iso.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 291 ]; then
    FILE=a20ra2
    CONFIG=roi/dehnen/a20ra2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 292 ]; then
    FILE=a20ra1
    CONFIG=roi/dehnen/a20ra1.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 293 ]; then
    FILE=a20ra1_2
    CONFIG=roi/dehnen/a20ra1_2.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 294 ]; then
    FILE=a20ra1_4
    CONFIG=roi/dehnen/a20ra1_4.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 295 ]; then
    FILE=a20ra1_8
    CONFIG=roi/dehnen/a20ra1_8.cfg
    EPS=7.8125e-3
    ETA=0.5
    FINISH=14000.0
    INTERVAL=100.0
fi
###############################################################
# dynamical stability of multi components galaxy model (only spherical components) with gaseous component
if [ $PROBLEM -eq 1000 ]; then
    FILE=etg_gas
    CONFIG=gas/etg.cfg
    GAS_CONFIG=gas/etg_gas.cfg
    EPS=1.5625e-2
    ETA=0.5
    INTERVAL=25.0
    FINISH=1175.0
fi
###############################################################
# set input arguments
if [ $PROBLEM -ge 1 ]; then
    OPTION="-file=$FILE -config=$CONFIG -Ntot=$NTOT -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE -denoisingDistributionFunction=$DENOISE_DF"
else
    OPTION="-file=$FILE -unit=$UNIT -Ntot=$NTOT -Mtot=$MTOT -virial=$VIRIAL -rad=$RAD -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE"
fi
###############################################################
if [ $PROBLEM -ge 1000 ]; then
    OPTION="$OPTION -gas_config=$GAS_CONFIG"
fi
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
exit 0
###############################################################
