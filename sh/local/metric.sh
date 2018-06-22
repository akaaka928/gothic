#!/bin/sh
###############################################################
if [ $# -lt 1 ]; then
    echo "$# input(s) is/are detected while at least 1 input is required to specify <test problem>" 1>&2
    exit 1
fi
JOB_ID=$$
PROBLEM=$1
###############################################################
#
#
###############################################################
# global configurations
###############################################################
EXEC=bin/gothic
###############################################################
#
#
###############################################################
# upper limit of execution time for nvprof --metric
if [ -z "$NVPROF_METRIC_TIMEOUT" ]; then
    NVPROF_METRIC_TIMEOUT=3600
fi
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # PROBLEM=2
    PROBLEM=27
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
REBUILD=16
BRENT=1.0
###############################################################
#
#
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
# dynamical stability of an M31 model (NFW halo, de Vaucouleurs bulge, and exponential disk)
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
# set input arguments
OPTION="-absErr=$ABSERR -file=$FILE -rebuild_interval=$REBUILD -brent_frac=$BRENT -jobID=$JOB_ID"
###############################################################
# nvprof --query-metrics
NVPROF_METRICS="--metrics flop_sp_efficiency,flop_dp_efficiency,gld_efficiency,gst_efficiency,shared_efficiency,sm_efficiency,branch_efficiency,warp_execution_efficiency,inst_fp_32,inst_fp_64,inst_integer,inst_bit_convert,inst_control,inst_compute_ld_st,inst_misc,inst_inter_thread_communication,flop_count_sp_fma,flop_count_sp_add,flop_count_sp_mul,flop_count_sp_special,flop_count_dp_fma,flop_count_dp_add,flop_count_dp_mul,single_precision_fu_utilization,double_precision_fu_utilization,special_fu_utilization,cf_fu_utilization,cf_issued,cf_executed,ldst_fu_utilization,ldst_issued,ldst_executed,tex_fu_utilization,shared_load_throughput,shared_store_throughput,shared_utilization,global_hit_rate,tex_cache_hit_rate,l2_tex_read_hit_rate,l2_tex_write_hit_rate,local_hit_rate,gld_requested_throughput,gld_throughput,gst_requested_throughput,gst_throughput,tex_cache_throughput,tex_utilization,l2_read_throughput,l2_tex_read_throughput,l2_write_throughput,l2_tex_write_throughput,l2_atomic_throughput,l2_utilization,local_memory_overhead,local_load_throughput,local_store_throughput,dram_read_throughput,dram_write_throughput,dram_utilization,ecc_throughput,sysmem_read_throughput,sysmem_read_utilization,sysmem_write_throughput,sysmem_write_utilization,sysmem_utilization,issue_slot_utilization,issue_slots,issued_ipc,ipc,stall_inst_fetch,stall_exec_dependency,stall_memory_dependency,stall_texture,stall_sync,stall_other,stall_constant_memory_dependency,stall_pipe_busy,stall_memory_throttle,stall_not_selected"
###############################################################
#
#
###############################################################
# execute numerical simulation
###############################################################
# start logging
LOG=log/${FILE}.l
PRFOUT=log/${FILE}_gothic.m${JOB_ID}
TIME=`date`
echo "start: $TIME" >> $LOG
###############################################################
# execute the job
if [ `which numactl` ]; then
    # run with numactl
    echo "numactl --cpunodebind=0 --localalloc nvprof --timeout $NVPROF_METRIC_TIMEOUT --kernels \"::calcAcc_kernel:\" $NVPROF_METRICS --kernels \"::calcMultipole_kernel:\" $NVPROF_METRICS --kernels \"::makeTree_kernel:\" $NVPROF_METRICS --kernels \"::calcPHkey_kernel:\" $NVPROF_METRICS --kernels \"::adjustTimeStep_kernel:\" $NVPROF_METRICS --kernels \"::prediction_kernel:\" $NVPROF_METRICS $EXEC $OPTION 1>>$PRFOUT 2>>&1" >> $LOG
    numactl --cpunodebind=0 --localalloc nvprof --timeout $NVPROF_METRIC_TIMEOUT --kernels "::calcAcc_kernel:" $NVPROF_METRICS --kernels "::calcMultipole_kernel:" $NVPROF_METRICS --kernels "::makeTree_kernel:" $NVPROF_METRICS --kernels "::calcPHkey_kernel:" $NVPROF_METRICS --kernels "::adjustTimeStep_kernel:" $NVPROF_METRICS --kernels "::prediction_kernel:" $NVPROF_METRICS $EXEC $OPTION 1>>$PRFOUT 2>>&1
else
    # run without numactl
    echo "nvprof --timeout $NVPROF_METRIC_TIMEOUT --kernels \"::calcAcc_kernel:\" $NVPROF_METRICS --kernels \"::calcMultipole_kernel:\" $NVPROF_METRICS --kernels \"::makeTree_kernel:\" $NVPROF_METRICS --kernels \"::calcPHkey_kernel:\" $NVPROF_METRICS --kernels \"::adjustTimeStep_kernel:\" $NVPROF_METRICS --kernels \"::prediction_kernel:\" $NVPROF_METRICS $EXEC $OPTION 1>>$PRFOUT 2>>&1" >> $LOG
    nvprof --timeout $NVPROF_METRIC_TIMEOUT --kernels "::calcAcc_kernel:" $NVPROF_METRICS --kernels "::calcMultipole_kernel:" $NVPROF_METRICS --kernels "::makeTree_kernel:" $NVPROF_METRICS --kernels "::calcPHkey_kernel:" $NVPROF_METRICS --kernels "::adjustTimeStep_kernel:" $NVPROF_METRICS --kernels "::prediction_kernel:" $NVPROF_METRICS $EXEC $OPTION 1>>$PRFOUT 2>>&1
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME" >> $LOG
###############################################################
