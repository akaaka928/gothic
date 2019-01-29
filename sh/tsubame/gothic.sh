#!/bin/sh
#$ -cwd
#$ -l s_gpu=1
#$ -l h_rt=0:20:00
#$ -N gothic
#$ -hold_jid magi,gothic
###############################################################


###############################################################
# global configurations
###############################################################
if [ -z "$EXEC" ]; then
    EXEC=bin/gothic
fi
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=2
    # PROBLEM=20
    # PROBLEM=26
    PROBLEM=27
    # PROBLEM=80
    # PROBLEM=81
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
    if [ -z "$FILE" ]; then
	FILE=m31
    fi
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
# Long-time stability of globular cluster
if [ $PROBLEM -eq 200 ]; then
    FILE=w3iso
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 201 ]; then
    FILE=w3ra2
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 202 ]; then
    FILE=w3ra1
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 203 ]; then
    FILE=w3ra1_2
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 204 ]; then
    FILE=w3ra1_4
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 205 ]; then
    FILE=w3ra1_8
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 210 ]; then
    FILE=w5iso
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 211 ]; then
    FILE=w5ra2
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 212 ]; then
    FILE=w5ra1
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 213 ]; then
    FILE=w5ra1_2
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 214 ]; then
    FILE=w5ra1_4
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 215 ]; then
    FILE=w5ra1_8
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 220 ]; then
    FILE=w7iso
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 221 ]; then
    FILE=w7ra2
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 222 ]; then
    FILE=w7ra1
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 223 ]; then
    FILE=w7ra1_2
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 224 ]; then
    FILE=w7ra1_4
fi
###############################################################
# Long-time stability of globular cluster
if [ $PROBLEM -eq 225 ]; then
    FILE=w7ra1_8
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 230 ]; then
    FILE=m12iso
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 231 ]; then
    FILE=m12ra2
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 232 ]; then
    FILE=m12ra1
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 233 ]; then
    FILE=m12ra1_2
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 234 ]; then
    FILE=m12ra1_4
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 235 ]; then
    FILE=m12ra1_8
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 240 ]; then
    FILE=m09iso
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 241 ]; then
    FILE=m09ra2
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 242 ]; then
    FILE=m09ra1
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 243 ]; then
    FILE=m09ra1_2
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 244 ]; then
    FILE=m09ra1_4
fi
###############################################################
# Long-time stability of DM halo
if [ $PROBLEM -eq 245 ]; then
    FILE=m09ra1_8
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 250 ]; then
    FILE=a00iso
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 251 ]; then
    FILE=a00ra2
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 252 ]; then
    FILE=a00ra1
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 253 ]; then
    FILE=a00ra1_2
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 254 ]; then
    FILE=a00ra1_4
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 260 ]; then
    FILE=a05iso
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 261 ]; then
    FILE=a05ra2
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 262 ]; then
    FILE=a05ra1
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 263 ]; then
    FILE=a05ra1_2
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 264 ]; then
    FILE=a05ra1_4
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 265 ]; then
    FILE=a05ra1_8
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 270 ]; then
    FILE=a10iso
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 271 ]; then
    FILE=a10ra2
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 272 ]; then
    FILE=a10ra1
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 273 ]; then
    FILE=a10ra1_2
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 274 ]; then
    FILE=a10ra1_4
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 275 ]; then
    FILE=a10ra1_8
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 280 ]; then
    FILE=a15iso
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 281 ]; then
    FILE=a15ra2
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 282 ]; then
    FILE=a15ra1
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 283 ]; then
    FILE=a15ra1_2
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 284 ]; then
    FILE=a15ra1_4
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 285 ]; then
    FILE=a15ra1_8
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 290 ]; then
    FILE=a20iso
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 291 ]; then
    FILE=a20ra2
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 292 ]; then
    FILE=a20ra1
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 293 ]; then
    FILE=a20ra1_2
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 294 ]; then
    FILE=a20ra1_4
fi
###############################################################
# Long-time stability of Dehnen model
if [ $PROBLEM -eq 295 ]; then
    FILE=a20ra1_8
fi
###############################################################
# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -Nx=$NX -Ny=$NY -Nz=$NZ -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$JOB_ID"
###############################################################


###############################################################
# job execution via UNIVA Grid Engine
###############################################################
# set stdout and stderr
STDOUT=log/${FILE}_$REQUEST.o$JOB_ID
STDERR=log/${FILE}_$REQUEST.e$JOB_ID
###############################################################
# load modules
. /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:/gs/hs1/jh180045/share/opt/Modules
module load intel cuda openmpi
module load cub phdf5/ompi
module list 1>>$STDOUT 2>>$STDERR
###############################################################
cat $PE_HOSTFILE 1>>$STDOUT 2>>$STDERR
TIME=`date`
echo "start: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################
# execute the job
if [ `which numactl` ]; then
    # run with numactl
    echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    # run without numactl
    echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################
