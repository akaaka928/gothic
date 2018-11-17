#################################################################################################
# last updated on 2018/09/10 (Mon) 14:31:34
# Makefile for C Programming
# Gravitational octree code for collisionless N-body simulations on GPUs
#################################################################################################


#################################################################################################
## Overall options
USEDBG	:= 0
USEDP	:= 0
MKOPREP	:= 0
# USEPAPI	:= 1
VERBOSE	:=
# VERBOSE	:= @
#################################################################################################
# Macros for Pre-Processor
DEBUG	:= -DNDEBUG
# PROFILE	:= -pg
#################################################################################################
# Execution options
FORCE_SINGLE_GPU_RUN	:= 1
DISABLE_AUTO_TUNING	:= 0
COMBINE_WITH_J_PARALLEL	:= 1
MPI_AUTO_TUNE_FOR_RMAX	:= 0
SKIP_UNUSED_LET_BUILD	:= 1
RECTANGULAR_BOX_FOR_LET	:= 1
ENCLOSING_BALL_FOR_LET	:= 1
COMMUNICATION_VIA_HOST	:= 0
USE_MPI_GET_FOR_LET	:= 1
USE_MPI_GET_FOR_EXCG	:= 0
CARE_EXTERNAL_PARTICLES	:= 0
ACC_ACCUMULATION_IN_DP	:= 0
KAHAN_SUM_CORRECTION	:= 0
MONITOR_ENERGY_ERROR	:= 1
MONITOR_LETGEN_TIME	:= 1
DIVERT_GEOMETRIC_CENTER	:= 1
ADOPT_BLOCK_TIME_STEP	:= 1
ADOPT_VECTOR_ACC_MAC	:= 0
ADOPT_GADGET_TYPE_MAC	:= 1
ADOPT_WS93_TYPE_MAC	:= 0
IJ_PARALLELIZED_WALK	:= 1
CLOCK_BASED_AUTO_TUNING	:= 1
OMIT_VELOCITY_BASED_DT	:= 0
REPORT_COMPUTE_RATE	:= 1
USE_COOPERATIVE_GROUPS	:= 0
GET_NBLOCKS_PER_SM_AUTO	:= 1
DISABLE_NVML_FOR_CLOCK	:= 1
DATAFILE_FORMAT_HDF5	:= 1
HDF5_FOR_ZINDAIJI	:= 1
PREPARE_XDMF_FILES	:= 1
DUMPFILE_IN_TIPSY	:= 0
DUMPFILE_AS_GALACTICS	:= 0
DUMPFILE_AS_GADGET	:= 0
USE_OFFICIAL_SFMT	:= 1
USE_OFFICIAL_SFMT_JUMP	:= 1
USE_LIS_FOR_MAGI	:= 1
SET_EXTERNAL_FIELD	:= 0
SET_EXTERNAL_FIELD_DISK	:= 1
ADAPTIVE_EXTERNAL_FIELD	:= 0
USE_OSIPKOV_MERRITT	:= 1
USE_ZH78_RETROGRADING	:= 0
#################################################################################################
# Debugging options
EVALUATE_FORCE_ERROR	:= 0
#################################################################################################
# Benchmark options
REPORT_ELAPSED_TIME	:= 1
REPORT_CLOCK_FREQUENCY	:= 0
MEASURE_EXEC_METRICS	:= 0
MEASURE_ELAPSED_TIME	:= 0
HUNT_OPTIMAL_WALK_TREE	:= 0
HUNT_OPTIMAL_INTEGRATE	:= 0
HUNT_OPTIMAL_MAKE_TREE	:= 0
HUNT_OPTIMAL_MAKE_NODE	:= 0
HUNT_OPTIMAL_NEIGHBOUR	:= 0
HUNT_OPTIMAL_SEPARATION	:= 0
MEASURE_NI_DEPENDENCE	:= 0
#################################################################################################
# Checking options
CHECK_FUNCTION_CALLS	:= 0
CHECK_CALL_MAKE_TREE	:= 0
CHECK_CALL_MAKE_NODE	:= 0
CHECK_CALL_NEIGHBOUR	:= 0
#################################################################################################


#################################################################################################
## Executing environment
HOSTNAME	:= $(shell hostname)
MYDIR	:= $(HOME)
MYINC	:= $(MYDIR)/inc
MYLIB	:= $(MYDIR)/lib
#################################################################################################
# Reedbush
ifeq ($(findstring reedbush, $(HOSTNAME)), reedbush)
MYDIR	:= /lustre/jh180045l/$(USER)
MYINC	:= $(MYDIR)/inc
MYLIB	:= $(MYDIR)/lib
endif
#################################################################################################
# TSUBAME 3.0
ifeq ($(findstring login, $(HOSTNAME)), login)
MYDIR	:= /gs/hs1/jh180045/$(USER)
MYINC	:= $(MYDIR)/inc
MYLIB	:= $(MYDIR)/lib
ifeq ($(USE_MPI_GET_FOR_LET), 1)
# tentative treatment to switch off one-sided communication on TSUBAME 3.0
USE_MPI_GET_FOR_LET	:= 0
endif
endif
#################################################################################################
# Server at ITC
ifeq ($(findstring pn, $(HOSTNAME)), pn)
MYDIR	:= $(HOME)/pn
MYINC	:= $(MYDIR)/inc
MYLIB	:= $(MYDIR)/lib
endif
#################################################################################################
include	$(MYINC)/common.mk
#################################################################################################
ifeq ($(USEHDF5), 0)
DATAFILE_FORMAT_HDF5	:= 0
endif
#################################################################################################
ifeq ($(GDR_AVAILABLE), 0)
COMMUNICATION_VIA_HOST	:= 0
endif
#################################################################################################


#################################################################################################
ifeq ($(MEASURE_EXEC_METRICS), 1)
MEASURE_ELAPSED_TIME	:= 1
endif
#################################################################################################
ifeq ($(MEASURE_ELAPSED_TIME), 1)
CCARG	+= -DEXEC_BENCHMARK
CUARG	+= -DEXEC_BENCHMARK
# CCARG	+= -DTFIN_MIN_BENCHMARK="(60)"
CCARG	+= -DTFIN_MIN_BENCHMARK="(1430)"
DEBUG	:= -DNDEBUG
PROFILE	:=
USEDBG	:= 0
USEPAPI	:= 0
EVALUATE_FORCE_ERROR	:= 0
# REPORT_ELAPSED_TIME	:= 0
REPORT_CLOCK_FREQUENCY	:= 0
endif
#################################################################################################
ifeq ($(EVALUATE_FORCE_ERROR), 1)
REPORT_ELAPSED_TIME	:= 0
REPORT_CLOCK_FREQUENCY	:= 0
endif
#################################################################################################
ifeq ($(REPORT_ELAPSED_TIME), 1)
CCARG	+= -DREPORT_TOTAL_ELAPSED_TIME
endif
#################################################################################################
ifneq ($(USENVML), 1)
DISABLE_NVML_FOR_CLOCK	:= 1
endif
ifeq ($(DISABLE_NVML_FOR_CLOCK), 1)
REPORT_CLOCK_FREQUENCY	:= 0
CCARG	+= -DDISABLE_NVML_FOR_CLOCK_FREQ
CUARG	+= -DDISABLE_NVML_FOR_CLOCK_FREQ
endif
#################################################################################################
ifeq ($(REPORT_CLOCK_FREQUENCY), 1)
CCARG	+= -DREPORT_GPU_CLOCK_FREQUENCY
CUARG	+= -DREPORT_GPU_CLOCK_FREQUENCY
endif
#################################################################################################
ifeq ($(CHECK_FUNCTION_CALLS), 1)
DEBUG	:=
PROFILE	:=
USEDBG	:= 0
USEPAPI	:= 0
EVALUATE_FORCE_ERROR	:= 0
REPORT_ELAPSED_TIME	:= 0
endif
#################################################################################################
ifeq ($(MEASURE_EXEC_METRICS), 1)
CCARG	+= -DCOUNT_INTERACTIONS
CUARG	+= -DCOUNT_INTERACTIONS
endif
#################################################################################################
ifeq ($(MONITOR_ENERGY_ERROR), 1)
CCARG	+= -DMONITOR_ENERGY_ERROR
endif
#################################################################################################
ifeq ($(CARE_EXTERNAL_PARTICLES), 1)
CCARG	+= -DCARE_EXTERNAL_PARTICLES
CUARG	+= -DCARE_EXTERNAL_PARTICLES
endif
#################################################################################################
ifeq ($(ACC_ACCUMULATION_IN_DP), 1)
ifeq ($(USEDP), 0)
KAHAN_SUM_CORRECTION	:= 0
CCARG	+= -DDPADD_FOR_ACC
CUARG	+= -DDPADD_FOR_ACC
endif
endif
#################################################################################################
ifeq ($(KAHAN_SUM_CORRECTION), 1)
CCARG	+= -DKAHAN_SUM_CORRECTION
CUARG	+= -DKAHAN_SUM_CORRECTION
endif
#################################################################################################
ifeq ($(FORCE_SINGLE_GPU_RUN), 1)
CCARG	+= -DSERIALIZED_EXECUTION
CUARG	+= -DSERIALIZED_EXECUTION
COMBINE_WITH_J_PARALLEL	:= 0
MPI_AUTO_TUNE_FOR_RMAX	:= 0
SKIP_UNUSED_LET_BUILD	:= 0
RECTANGULAR_BOX_FOR_LET	:= 0
ENCLOSING_BALL_FOR_LET	:= 0
COMMUNICATION_VIA_HOST	:= 0
USE_MPI_GET_FOR_LET	:= 0
USE_MPI_GET_FOR_EXCG	:= 0
CARE_EXTERNAL_PARTICLES	:= 0
else
DISABLE_AUTO_TUNING	:= 1
# direct solver is not parallelized
EVALUATE_FORCE_ERROR	:= 0
ifeq ($(COMBINE_WITH_J_PARALLEL), 1)
CCARG	+= -DSWITCH_WITH_J_PARALLELIZATION
CUARG	+= -DSWITCH_WITH_J_PARALLELIZATION
endif
ifeq ($(MPI_AUTO_TUNE_FOR_RMAX), 1)
CCARG	+= -DMPI_MAX_FOR_RMAX_IN_AUTO_TUNING
CUARG	+= -DMPI_MAX_FOR_RMAX_IN_AUTO_TUNING
endif
ifeq ($(MONITOR_LETGEN_TIME), 1)
CCARG	+= -DMONITOR_LETGEN_TIME
CUARG	+= -DMONITOR_LETGEN_TIME
endif
ifeq ($(COMMUNICATION_VIA_HOST), 1)
CCARG	+= -DMPI_VIA_HOST
CUARG	+= -DMPI_VIA_HOST
endif
ifeq ($(USE_MPI_GET_FOR_LET), 1)
CCARG	+= -DMPI_ONE_SIDED_FOR_LET
CUARG	+= -DMPI_ONE_SIDED_FOR_LET
endif
ifeq ($(USE_MPI_GET_FOR_EXCG), 1)
CCARG	+= -DMPI_ONE_SIDED_FOR_EXCG
CUARG	+= -DMPI_ONE_SIDED_FOR_EXCG
endif
ifeq ($(SKIP_UNUSED_LET_BUILD), 1)
CCARG	+= -DSKIP_UNUSED_LET_GENERATION
CUARG	+= -DSKIP_UNUSED_LET_GENERATION
endif
ifeq ($(ENCLOSING_BALL_FOR_LET), 1)
CCARG	+= -DUSE_ENCLOSING_BALL_FOR_LET
CUARG	+= -DUSE_ENCLOSING_BALL_FOR_LET
DIVERT_GEOMETRIC_CENTER	:= 1
else
RECTANGULAR_BOX_FOR_LET	:= 1
DIVERT_GEOMETRIC_CENTER	:= 0
endif
ifeq ($(DIVERT_GEOMETRIC_CENTER), 1)
CCARG	+= -DRETURN_CENTER_BY_PHKEY_GENERATOR
CUARG	+= -DRETURN_CENTER_BY_PHKEY_GENERATOR
endif
ifeq ($(RECTANGULAR_BOX_FOR_LET), 1)
CCARG	+= -DUSE_RECTANGULAR_BOX_FOR_LET
CUARG	+= -DUSE_RECTANGULAR_BOX_FOR_LET
# ENCLOSING_BALL_FOR_LET	:= 0
# DIVERT_GEOMETRIC_CENTER	:= 0
endif
endif
#################################################################################################
ifeq ($(DISABLE_AUTO_TUNING), 1)
CCARG	+= -DDISABLE_AUTO_TUNING
CUARG	+= -DDISABLE_AUTO_TUNING
endif
#################################################################################################
ifeq ($(ADOPT_BLOCK_TIME_STEP), 1)
CCARG	+= -DBLOCK_TIME_STEP
CUARG	+= -DBLOCK_TIME_STEP
endif
#################################################################################################
ifeq ($(ADOPT_VECTOR_ACC_MAC), 1)
CCARG	+= -DYMIKI_MAC
CUARG	+= -DYMIKI_MAC
# YMIKI_MAC is very similar to GADGET MAC
ADOPT_GADGET_TYPE_MAC	:= 1
# LET generator for YMIKI_MAC is not yet implementated (also, would not been implemented in future)
FORCE_SINGLE_GPU_RUN	:= 1
endif
#################################################################################################
ifeq ($(ADOPT_GADGET_TYPE_MAC), 1)
CCARG	+= -DGADGET_MAC
CUARG	+= -DGADGET_MAC
ADOPT_WS93_TYPE_MAC	:= 0
endif
#################################################################################################
ifeq ($(ADOPT_WS93_TYPE_MAC), 1)
CCARG	+= -DWS93_MAC
CUARG	+= -DWS93_MAC
endif
#################################################################################################
ifeq ($(EVALUATE_FORCE_ERROR), 1)
CCARG	+= -DCOMPARE_WITH_DIRECT_SOLVER
CUARG	+= -DCOMPARE_WITH_DIRECT_SOLVER
endif
#################################################################################################
ifeq ($(IJ_PARALLELIZED_WALK), 1)
CCARG	+= -DIJ_PARALLELIZATION
CUARG	+= -DIJ_PARALLELIZATION
endif
#################################################################################################
ifeq ($(CLOCK_BASED_AUTO_TUNING), 1)
CCARG	+= -DUSE_CLOCK_CYCLES_FOR_BRENT_METHOD
CUARG	+= -DUSE_CLOCK_CYCLES_FOR_BRENT_METHOD
endif
#################################################################################################
ifeq ($(OMIT_VELOCITY_BASED_DT), 1)
CCARG	+= -DOMIT_VELOCITY_FOR_TIME_STEP
CUARG	+= -DOMIT_VELOCITY_FOR_TIME_STEP
endif
#################################################################################################
ifeq ($(REPORT_COMPUTE_RATE), 1)
CCARG	+= -DREPORT_COMPUTE_RATE
endif
#################################################################################################
ifeq ($(COOPERATIVE_GROUPS_AVAILABLE), 1)
ifeq ($(USE_COOPERATIVE_GROUPS), 1)
CCARG	+= -DUSE_COOPERATIVE_GROUPS
CUARG	+= -DUSE_COOPERATIVE_GROUPS
endif
else
USE_COOPERATIVE_GROUPS	:= 0
endif
#################################################################################################
ifeq ($(GET_NBLOCKS_PER_SM_AUTO), 1)
CCARG	+= -DUSE_OCCUPANCY_CALCULATOR
CUARG	+= -DUSE_OCCUPANCY_CALCULATOR
endif
#################################################################################################
ifeq ($(USEPAPI), 0)
PAPIOPT	:=
endif
#################################################################################################
ifeq ($(USELIS), 0)
USE_LIS_FOR_MAGI	:= 0
endif
ifeq ($(USE_LIS_FOR_MAGI), 1)
CCARG	+= -DUSE_LIS -DDISABLE_MPI_LIS
endif
#################################################################################################
ifeq ($(USESFMT), 0)
USE_OFFICIAL_SFMT	:= 0
USE_OFFICIAL_SFMT_JUMP	:= 0
endif
ifeq ($(USE_OFFICIAL_SFMT), 1)
CCARG	+= -DUSE_SFMT
endif
#################################################################################################
ifeq ($(USESMTJ), 0)
USE_OFFICIAL_SFMT_JUMP	:= 0
endif
ifeq ($(USE_OFFICIAL_SFMT_JUMP), 1)
CCARG	+= -DUSE_SFMTJUMP
endif
#################################################################################################
ifeq ($(SET_EXTERNAL_FIELD), 1)
CCARG	+= -DSET_EXTERNAL_POTENTIAL_FIELD
CUARG	+= -DSET_EXTERNAL_POTENTIAL_FIELD
else
SET_EXTERNAL_FIELD_DISK	:= 0
ADAPTIVE_EXTERNAL_FIELD	:= 0
endif
#################################################################################################
ifeq ($(SET_EXTERNAL_FIELD_DISK), 1)
CCARG	+= -DSET_EXTERNAL_POTENTIAL_FIELD_DISK
CUARG	+= -DSET_EXTERNAL_POTENTIAL_FIELD_DISK
endif
#################################################################################################
ifeq ($(ADAPTIVE_EXTERNAL_FIELD), 1)
CCARG	+= -DADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
CUARG	+= -DADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
endif
#################################################################################################
ifeq ($(USE_OSIPKOV_MERRITT), 1)
CCARG	+= -DUSE_OSIPKOV_MERRITT_METHOD
endif
#################################################################################################
ifeq ($(USE_ZH78_RETROGRADING), 1)
CCARG	+= -DUSE_ZANG_HOHL_1978_EQ5
endif
#################################################################################################
NUM_NTHREADS	:= 512
NUM_TSUB	:= 32
NUM_NWARP	:= 4
NUM_NLOOP	:= 1
LEV_NEIGHBOR	:= 1
USE_WARPSHUFFLE	:= 1
PREF_SHARED_MEM	:= 1
#################################################################################################
ifeq ($(MEASURE_ELAPSED_TIME), 1)
#################################################################################################
ifeq ($(HUNT_OPTIMAL_SEPARATION), 1)
HUNT_OPTIMAL_WALK_TREE	:= 1
HUNT_OPTIMAL_INTEGRATE	:= 1
HUNT_OPTIMAL_MAKE_TREE	:= 0
HUNT_OPTIMAL_MAKE_NODE	:= 0
HUNT_OPTIMAL_NEIGHBOUR	:= 0
endif
#################################################################################################
ifeq ($(HUNT_OPTIMAL_WALK_TREE), 1)
CCARG	+= -DHUNT_WALK_PARAMETER
CUARG	+= -DHUNT_WALK_PARAMETER
CCARG	+= -DNTHREADS="($(NUM_NTHREADS))" -DTSUB="($(NUM_TSUB))" -DNLOOP="($(NUM_NLOOP))"
CUARG	+= -DNTHREADS="($(NUM_NTHREADS))" -DTSUB="($(NUM_TSUB))" -DNLOOP="($(NUM_NLOOP))"
ifeq ($(IJ_PARALLELIZED_WALK), 1)
CCARG	+= -DNWARP="($(NUM_NWARP))"
CUARG	+= -DNWARP="($(NUM_NWARP))"
endif
ifeq ($(USE_WARPSHUFFLE), 1)
CCARG	+= -DUSE_WARP_SHUFFLE_FUNC
CUARG	+= -DUSE_WARP_SHUFFLE_FUNC
endif
endif
#################################################################################################
ifeq ($(HUNT_OPTIMAL_MAKE_TREE), 1)
CCARG	+= -DHUNT_MAKE_PARAMETER
CUARG	+= -DHUNT_MAKE_PARAMETER
CCARG	+= -DNTHREADS_MAKE_TREE="($(NUM_NTHREADS))" -DNTHREADS_LINK_TREE="($(NUM_NTHREADS))" -DNTHREADS_TRIM_TREE="($(NUM_NTHREADS))"
CUARG	+= -DNTHREADS_MAKE_TREE="($(NUM_NTHREADS))" -DNTHREADS_LINK_TREE="($(NUM_NTHREADS))" -DNTHREADS_TRIM_TREE="($(NUM_NTHREADS))"
CCARG	+= -DNTHREADS_INIT_LINK="($(NUM_NTHREADS))" -DNTHREADS_INIT_CELL="($(NUM_NTHREADS))" -DNTHREADS_INIT_NODE="($(NUM_NTHREADS))"
CUARG	+= -DNTHREADS_INIT_LINK="($(NUM_NTHREADS))" -DNTHREADS_INIT_CELL="($(NUM_NTHREADS))" -DNTHREADS_INIT_NODE="($(NUM_NTHREADS))"
CCARG	+= -DNTHREADS_INIT_BODY="($(NUM_NTHREADS))" -DNTHREADS_COPY_BODY="($(NUM_NTHREADS))"
CUARG	+= -DNTHREADS_INIT_BODY="($(NUM_NTHREADS))" -DNTHREADS_COPY_BODY="($(NUM_NTHREADS))"
CCARG	+= -DNTHREADS_PH="($(NUM_NTHREADS))" -DNTHREADS_PHSORT="($(NUM_NTHREADS))"
CUARG	+= -DNTHREADS_PH="($(NUM_NTHREADS))" -DNTHREADS_PHSORT="($(NUM_NTHREADS))"
ifeq ($(USE_WARPSHUFFLE), 1)
CCARG	+= -DUSE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
CUARG	+= -DUSE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
endif
endif
#################################################################################################
ifeq ($(HUNT_OPTIMAL_MAKE_NODE), 1)
CCARG	+= -DHUNT_NODE_PARAMETER
CUARG	+= -DHUNT_NODE_PARAMETER
CCARG	+= -DNTHREADS_MAC="($(NUM_NTHREADS))" -DTSUB_MAC="($(NUM_TSUB))"
CUARG	+= -DNTHREADS_MAC="($(NUM_NTHREADS))" -DTSUB_MAC="($(NUM_TSUB))"
ifeq ($(USE_WARPSHUFFLE), 1)
CCARG	+= -DUSE_WARP_SHUFFLE_FUNC_MAC
CUARG	+= -DUSE_WARP_SHUFFLE_FUNC_MAC
endif
endif
#################################################################################################
ifeq ($(HUNT_OPTIMAL_NEIGHBOUR), 1)
CCARG	+= -DHUNT_FIND_PARAMETER
CUARG	+= -DHUNT_FIND_PARAMETER
CCARG	+= -DNTHREADS_SHRINK="($(NUM_NTHREADS))" -DNTHREADS_FACILE_NS="($(NUM_NTHREADS))" -DNTHREADS_NEIGHBOR="($(NUM_NTHREADS))" -DTSUB_NEIGHBOR="($(NUM_TSUB))"
CUARG	+= -DNTHREADS_SHRINK="($(NUM_NTHREADS))" -DNTHREADS_FACILE_NS="($(NUM_NTHREADS))" -DNTHREADS_NEIGHBOR="($(NUM_NTHREADS))" -DTSUB_NEIGHBOR="($(NUM_TSUB))"
ifeq ($(USE_WARPSHUFFLE), 1)
CCARG	+= -DUSE_WARP_SHUFFLE_FUNC_NEIGHBOR
CUARG	+= -DUSE_WARP_SHUFFLE_FUNC_NEIGHBOR
endif
ifeq ($(PREF_SHARED_MEM), 1)
CCARG	+= -DSMEM_PREF_FOR_NEIGHBOR_SEARCH
CUARG	+= -DSMEM_PREF_FOR_NEIGHBOR_SEARCH
endif
endif
#################################################################################################
ifeq ($(HUNT_OPTIMAL_INTEGRATE), 1)
CCARG	+= -DHUNT_TIME_PARAMETER
CUARG	+= -DHUNT_TIME_PARAMETER
CCARG	+= -DNTHREADS_TIME="($(NUM_NTHREADS))"
CUARG	+= -DNTHREADS_TIME="($(NUM_NTHREADS))"
ifeq ($(USE_WARPSHUFFLE), 1)
CCARG	+= -DUSE_WARP_SHUFFLE_FUNC_TIME
CUARG	+= -DUSE_WARP_SHUFFLE_FUNC_TIME
endif
endif
#################################################################################################
ifeq ($(ADOPT_BLOCK_TIME_STEP), 1)
ifeq ($(MEASURE_NI_DEPENDENCE), 1)
CCARG	+= -DSHOW_NI_DEPENDENCE
CUARG	+= -DSHOW_NI_DEPENDENCE
endif
endif
#################################################################################################
endif
#################################################################################################
ifeq ($(CHECK_FUNCTION_CALLS), 1)
#################################################################################################
ifeq ($(CHECK_CALL_MAKE_TREE), 1)
CCARG	+= -DNTHREADS_MAKE_TREE="($(NUM_NTHREADS))" -DNTHREADS_LINK_TREE="($(NUM_NTHREADS))" -DNTHREADS_TRIM_TREE="($(NUM_NTHREADS))"
CUARG	+= -DNTHREADS_MAKE_TREE="($(NUM_NTHREADS))" -DNTHREADS_LINK_TREE="($(NUM_NTHREADS))" -DNTHREADS_TRIM_TREE="($(NUM_NTHREADS))"
CCARG	+= -DNTHREADS_PH="($(NUM_NTHREADS))"
CUARG	+= -DNTHREADS_PH="($(NUM_NTHREADS))"
ifeq ($(USE_WARPSHUFFLE), 1)
CCARG	+= -DUSE_WARP_SHUFFLE_FUNC_MAKE_TREE -DUSE_WARP_SHUFFLE_FUNC_LINK_TREE
CUARG	+= -DUSE_WARP_SHUFFLE_FUNC_MAKE_TREE -DUSE_WARP_SHUFFLE_FUNC_LINK_TREE
endif
endif
#################################################################################################
ifeq ($(CHECK_CALL_MAKE_NODE), 1)
CCARG	+= -DNTHREADS_MAC="($(NUM_NTHREADS))" -DTSUB_MAC="($(NUM_TSUB))"
CUARG	+= -DNTHREADS_MAC="($(NUM_NTHREADS))" -DTSUB_MAC="($(NUM_TSUB))"
ifeq ($(USE_WARPSHUFFLE), 1)
CCARG	+= -DUSE_WARP_SHUFFLE_FUNC_MAC
CUARG	+= -DUSE_WARP_SHUFFLE_FUNC_MAC
endif
endif
#################################################################################################
ifeq ($(CHECK_CALL_NEIGHBOUR), 1)
CCARG	+= -DNTHREADS_SHRINK="($(NUM_NTHREADS))" -DNTHREADS_FACILE_NS="($(NUM_NTHREADS))" -DNTHREADS_NEIGHBOR="($(NUM_NTHREADS))" -DTSUB_NEIGHBOR="($(NUM_TSUB))"
CUARG	+= -DNTHREADS_SHRINK="($(NUM_NTHREADS))" -DNTHREADS_FACILE_NS="($(NUM_NTHREADS))" -DNTHREADS_NEIGHBOR="($(NUM_NTHREADS))" -DTSUB_NEIGHBOR="($(NUM_TSUB))"
ifeq ($(USE_WARPSHUFFLE), 1)
CCARG	+= -DUSE_WARP_SHUFFLE_FUNC_NEIGHBOR
CUARG	+= -DUSE_WARP_SHUFFLE_FUNC_NEIGHBOR
endif
ifeq ($(PREF_SHARED_MEM), 1)
CCARG	+= -DSMEM_PREF_FOR_NEIGHBOR_SEARCH
CUARG	+= -DSMEM_PREF_FOR_NEIGHBOR_SEARCH
endif
endif
#################################################################################################
endif
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
CCARG	+= -DUSE_HDF5_FORMAT
CUARG	+= -DUSE_HDF5_FORMAT
else
HDF5_FOR_ZINDAIJI	:= 0
endif
#################################################################################################
ifeq ($(DUMPFILE_IN_TIPSY), 1)
CCARG	+= -DWRITE_IN_TIPSY_FORMAT
endif
#################################################################################################
ifeq ($(DUMPFILE_AS_GALACTICS), 1)
CCARG	+= -DWRITE_IN_GALACTICS_FORMAT
endif
#################################################################################################
ifeq ($(DUMPFILE_AS_GADGET), 1)
CCARG	+= -DWRITE_IN_GADGET_FORMAT
endif
#################################################################################################
ifeq ($(HDF5_FOR_ZINDAIJI), 1)
CCARG	+= -DHDF5_FOR_ZINDAIJI
endif
#################################################################################################
ifeq ($(PREPARE_XDMF_FILES), 1)
CCARG	+= -DPREPARE_XDMF_FILES
endif
#################################################################################################


#################################################################################################
## Source Directories
SRCDIR	:= src
MAINDIR	:= $(SRCDIR)/main
MISCDIR	:= $(SRCDIR)/misc
UTILDIR	:= $(SRCDIR)/util
FILEDIR	:= $(SRCDIR)/file
TIMEDIR	:= $(SRCDIR)/time
SORTDIR	:= $(SRCDIR)/sort
TREEDIR	:= $(SRCDIR)/tree
INITDIR	:= $(SRCDIR)/init
PLOTDIR	:= $(SRCDIR)/plot
PARADIR	:= $(SRCDIR)/para
ANALDIR	:= $(SRCDIR)/anal
VPATH	:= $(MAINDIR) $(MISCDIR) $(UTILDIR) $(FILEDIR) $(TIMEDIR) $(SORTDIR) $(TREEDIR) $(INITDIR) $(PLOTDIR) $(PARADIR) $(ANALDIR)
#################################################################################################
## Executables for collisionless N-body code
GOTHIC	:= $(BINDIR)/gothic
MKCOLD	:= $(BINDIR)/uniformsphere
MAGI	:= $(BINDIR)/magi
EDITOR	:= $(BINDIR)/editor
ifneq ($(PLPLOT), 0)
PLTENE	:= $(BINDIR)/energy
PLTDST	:= $(BINDIR)/distribution
PLTCDF	:= $(BINDIR)/cdf
PLTELP	:= $(BINDIR)/elapsed
PLTDEP	:= $(BINDIR)/ndep
PLTBRK	:= $(BINDIR)/breakdown
PLTFLP	:= $(BINDIR)/performance
PLTRAD	:= $(BINDIR)/ball
PLTDF	:= $(BINDIR)/df
PLTJET	:= $(BINDIR)/needle
PLTDISK	:= $(BINDIR)/disk
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
PLTCMP	:= $(BINDIR)/comparison
endif
endif
ANALACT	:= $(BINDIR)/action
ANALERR	:= $(BINDIR)/error
ANALPRF	:= $(BINDIR)/extract
M31OBS	:= $(BINDIR)/m31obs
OPTCFG	:= $(BINDIR)/showOptConfig
SAMPLE	:= $(BINDIR)/sample
#################################################################################################
## Plugins for VisIt to read HDF5 files
PLGDIR	:= plugins
TAGBODY	:= GOTHIC_split
TAGSNAP	:= GOTHIC_snp
TAGPLOT	:= GOTHIC_plt
TAGPM31	:= GOTHIC_m31
TAGAERR	:= GOTHIC_err
TAGDUMP	:= GOTHIC
TAGDISK	:= MAGI_disk
TAGDIST	:= MAGI_df
TAGPROF	:= MAGI_profile
DIRBODY	:= $(PLGDIR)/$(TAGBODY)
DIRSNAP	:= $(PLGDIR)/$(TAGSNAP)
DIRPLOT	:= $(PLGDIR)/$(TAGPLOT)
DIRPM31	:= $(PLGDIR)/$(TAGPM31)
DIRAERR	:= $(PLGDIR)/$(TAGAERR)
DIRDUMP	:= $(PLGDIR)/$(TAGDUMP)
DIRDISK	:= $(PLGDIR)/$(TAGDISK)
DIRDIST	:= $(PLGDIR)/$(TAGDIST)
DIRPROF	:= $(PLGDIR)/$(TAGPROF)
XMLBODY	:= $(DIRBODY)/$(TAGBODY).xml
XMLSNAP	:= $(DIRSNAP)/$(TAGSNAP).xml
XMLPLOT	:= $(DIRPLOT)/$(TAGPLOT).xml
XMLPM31	:= $(DIRPM31)/$(TAGPM31).xml
XMLAERR	:= $(DIRAERR)/$(TAGAERR).xml
XMLDUMP	:= $(DIRDUMP)/$(TAGDUMP).xml
XMLDISK	:= $(DIRDISK)/$(TAGDISK).xml
XMLDIST	:= $(DIRDIST)/$(TAGDIST).xml
XMLPROF	:= $(DIRPROF)/$(TAGPROF).xml
SRCBODY	:= $(DIRBODY)/$(TAGBODY)CommonPluginInfo.C $(DIRBODY)/$(TAGBODY)EnginePluginInfo.C $(DIRBODY)/$(TAGBODY)MDServerPluginInfo.C
SRCSNAP	:= $(DIRSNAP)/$(TAGSNAP)CommonPluginInfo.C $(DIRSNAP)/$(TAGSNAP)EnginePluginInfo.C $(DIRSNAP)/$(TAGSNAP)MDServerPluginInfo.C
SRCPLOT	:= $(DIRPLOT)/$(TAGPLOT)CommonPluginInfo.C $(DIRPLOT)/$(TAGPLOT)EnginePluginInfo.C $(DIRPLOT)/$(TAGPLOT)MDServerPluginInfo.C
SRCPM31	:= $(DIRPM31)/$(TAGPM31)CommonPluginInfo.C $(DIRPM31)/$(TAGPM31)EnginePluginInfo.C $(DIRPM31)/$(TAGPM31)MDServerPluginInfo.C
SRCAERR	:= $(DIRAERR)/$(TAGAERR)CommonPluginInfo.C $(DIRAERR)/$(TAGAERR)EnginePluginInfo.C $(DIRAERR)/$(TAGAERR)MDServerPluginInfo.C
SRCDUMP	:= $(DIRDUMP)/$(TAGDUMP)CommonPluginInfo.C $(DIRDUMP)/$(TAGDUMP)EnginePluginInfo.C $(DIRDUMP)/$(TAGDUMP)MDServerPluginInfo.C
SRCDISK	:= $(DIRDISK)/$(TAGDISK)CommonPluginInfo.C $(DIRDISK)/$(TAGDISK)EnginePluginInfo.C $(DIRDISK)/$(TAGDISK)MDServerPluginInfo.C
SRCDIST	:= $(DIRDIST)/$(TAGDIST)CommonPluginInfo.C $(DIRDIST)/$(TAGDIST)EnginePluginInfo.C $(DIRDIST)/$(TAGDIST)MDServerPluginInfo.C
SRCPROF	:= $(DIRPROF)/$(TAGPROF)CommonPluginInfo.C $(DIRPROF)/$(TAGPROF)EnginePluginInfo.C $(DIRPROF)/$(TAGPROF)MDServerPluginInfo.C
SRCBODY	+= $(DIRBODY)/$(TAGBODY)PluginInfo.C $(DIRBODY)/$(TAGBODY)PluginInfo.h $(DIRBODY)/avt$(TAGBODY)FileFormat.C $(DIRBODY)/avt$(TAGBODY)FileFormat.h
SRCSNAP	+= $(DIRSNAP)/$(TAGSNAP)PluginInfo.C $(DIRSNAP)/$(TAGSNAP)PluginInfo.h $(DIRSNAP)/avt$(TAGSNAP)FileFormat.C $(DIRSNAP)/avt$(TAGSNAP)FileFormat.h
SRCPLOT	+= $(DIRPLOT)/$(TAGPLOT)PluginInfo.C $(DIRPLOT)/$(TAGPLOT)PluginInfo.h $(DIRPLOT)/avt$(TAGPLOT)FileFormat.C $(DIRPLOT)/avt$(TAGPLOT)FileFormat.h
SRCPM31	+= $(DIRPM31)/$(TAGPM31)PluginInfo.C $(DIRPM31)/$(TAGPM31)PluginInfo.h $(DIRPM31)/avt$(TAGPM31)FileFormat.C $(DIRPM31)/avt$(TAGPM31)FileFormat.h
SRCAERR	+= $(DIRAERR)/$(TAGAERR)PluginInfo.C $(DIRAERR)/$(TAGAERR)PluginInfo.h $(DIRAERR)/avt$(TAGAERR)FileFormat.C $(DIRAERR)/avt$(TAGAERR)FileFormat.h
SRCDUMP	+= $(DIRDUMP)/$(TAGDUMP)PluginInfo.C $(DIRDUMP)/$(TAGDUMP)PluginInfo.h $(DIRDUMP)/avt$(TAGDUMP)FileFormat.C $(DIRDUMP)/avt$(TAGDUMP)FileFormat.h
SRCDISK	+= $(DIRDISK)/$(TAGDISK)PluginInfo.C $(DIRDISK)/$(TAGDISK)PluginInfo.h $(DIRDISK)/avt$(TAGDISK)FileFormat.C $(DIRDISK)/avt$(TAGDISK)FileFormat.h
SRCDIST	+= $(DIRDIST)/$(TAGDIST)PluginInfo.C $(DIRDIST)/$(TAGDIST)PluginInfo.h $(DIRDIST)/avt$(TAGDIST)FileFormat.C $(DIRDIST)/avt$(TAGDIST)FileFormat.h
SRCPROF	+= $(DIRPROF)/$(TAGPROF)PluginInfo.C $(DIRPROF)/$(TAGPROF)PluginInfo.h $(DIRPROF)/avt$(TAGPROF)FileFormat.C $(DIRPROF)/avt$(TAGPROF)FileFormat.h
#################################################################################################
## source files for collisionless N-body code relying on Barnes--Hut tree
NBDYSRC	:= gothic.c
FINDSRC	:= showOptConfig.c
#################################################################################################
ALLCLIB	:= allocate.c
CNVTLIB	:= convert.c
TUNELIB	:= tune.c
BRNTLIB	:= brent.c
UTILGPU	:= allocate_dev.cu
#################################################################################################
FILELIB	:= io.c
#################################################################################################
# TREESRC	:= peano.c
UTILGPU	+= peano_dev.cu
#################################################################################################
# TREESRC	:= make.c
TREEGPU	:= neighbor_dev.cu
ifeq ($(FORCE_SINGLE_GPU_RUN), 0)
LETHOST	:= let.c
LET_GPU	:= let_dev.cu
LET_GPU	+= adv_dev.cu make_dev.cu walk_dev.cu
GEO_GPU	:= geo_dev.cu
ifeq ($(ENCLOSING_BALL_FOR_LET), 1)
GEO_GPU	+= icom_dev.cu
endif
else
TREEGPU	+= adv_dev.cu make_dev.cu walk_dev.cu
endif
ifeq ($(MPI_AUTO_TUNE_FOR_RMAX), 1)
LET_GPU	+= shrink_dev.cu
else
UTILGPU	+= shrink_dev.cu
endif
ifeq ($(SET_EXTERNAL_FIELD), 1)
POT_EXT	:= potential_dev.cu
endif
#################################################################################################
COLDSRC	:= uniformsphere.c
MAGISRC	:= magi.c
SMPLSRC	:= sample.c
DISKLIB	:= potdens.c diskDF.c
MAGILIB	:= profile.c eddington.c king.c abel.c blas.c spline.c table.c
ifeq ($(SET_EXTERNAL_FIELD), 1)
MAGILIB	+= external.c
endif
EDITSRC	:= editor.c
#################################################################################################
PENESRC	:= energy.c
DISTSRC	:= distribution.c
PCDFSRC	:= cdf.c
PCDFLIB	:= cdflib.c
PELPSRC	:= elapsed.c
PDEPSRC	:= ndep.c
PBRKSRC	:= breakdown.c
PCMPSRC	:= comparison.c
PFLPSRC	:= performance.c
PRADSRC	:= ball.c
PDFSRC	:= df.c
PJETSRC	:= needle.c
PDSKSRC	:= disk.c
#################################################################################################
LETHOST	+= exchange.c mpicfg.c
EXCGGPU	:= exchange_dev.cu
#################################################################################################
AACTSRC	:= action.c
AERRSRC	:= error.c
APRFSRC	:= extract.c
AM31SRC	:= m31obs.c
M31LIB	:= m31coord.c
#################################################################################################


#################################################################################################
## Object files
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJMPGT	:= $(patsubst %.c,  $(OBJDIR)/%.mpi.hdf5.o, $(notdir $(NBDYSRC) $(FILELIB) $(ALLCLIB) $(CNVTLIB)))
ifeq ($(SET_EXTERNAL_FIELD), 1)
OBJMPGT	+= $(patsubst %.cu, $(OBJDIR)/%.mpi.hdf5.o, $(notdir $(POT_EXT)))
endif
else
OBJMPGT	:= $(patsubst %.c,  $(OBJDIR)/%.mpi.o,      $(notdir $(NBDYSRC) $(FILELIB)))
OBJMPGT	+= $(patsubst %.c,  $(OBJDIR)/%.o,          $(notdir $(ALLCLIB) $(CNVTLIB)))
ifeq ($(SET_EXTERNAL_FIELD), 1)
OBJMPGT	+= $(patsubst %.cu, $(OBJDIR)/%.mpi.o,      $(notdir $(POT_EXT)))
endif
endif
OBJMPGT	+= $(patsubst %.c,  $(OBJDIR)/%.o,          $(notdir $(TREESRC) $(TUNELIB)))
OBJMPGT	+= $(patsubst %.cu, $(OBJDIR)/%.o,          $(notdir $(UTILGPU) $(TREEGPU)))
ifeq ($(FORCE_SINGLE_GPU_RUN), 0)
OBJMPGT	+= $(patsubst %.cu, $(OBJDIR)/%.mpi.o,      $(notdir $(LET_GPU)))
OBJMPGT	+= $(patsubst %.cu, $(OBJDIR)/%.o,          $(notdir $(GEO_GPU)))
OBJMPGT	+= $(patsubst %.c,  $(OBJDIR)/%.mpi.o,      $(notdir $(LETHOST)))
OBJMPGT	+= $(patsubst %.cu, $(OBJDIR)/%.mpi.o,      $(notdir $(EXCGGPU)))
endif
OBJMPGT	+= $(patsubst %.c,  $(OBJDIR)/%.o,          $(notdir $(BRNTLIB)))
#################################################################################################
OBJFIND	:= $(patsubst %.c, $(OBJDIR)/%.o, $(notdir $(FINDSRC)))
OBJSMPL	:= $(patsubst %.c, $(OBJDIR)/%.o, $(notdir $(SMPLSRC)))
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
ifeq ($(USE_OFFICIAL_SFMT), 1)
ifeq ($(USE_OFFICIAL_SFMT_JUMP), 1)
OBJCOLD	:= $(patsubst %.c, $(OBJDIR)/%.mpi.smtj.hdf5.o, $(notdir $(COLDSRC)))
else
OBJCOLD	:= $(patsubst %.c, $(OBJDIR)/%.mpi.sfmt.hdf5.o, $(notdir $(COLDSRC)))
endif
else
OBJCOLD	:= $(patsubst %.c, $(OBJDIR)/%.mpi.gsl.hdf5.o,  $(notdir $(COLDSRC)))
endif
OBJCOLD	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,      $(notdir $(FILELIB) $(ALLCLIB)))
else
ifeq ($(USE_OFFICIAL_SFMT), 1)
OBJCOLD	:= $(patsubst %.c, $(OBJDIR)/%.sfmt.o,          $(notdir $(COLDSRC)))
else
OBJCOLD	:= $(patsubst %.c, $(OBJDIR)/%.gsl.o,           $(notdir $(COLDSRC)))
endif
OBJCOLD	+= $(patsubst %.c, $(OBJDIR)/%.o,               $(notdir $(ALLCLIB)))
OBJCOLD	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,           $(notdir $(FILELIB)))
endif
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
ifeq ($(USE_OFFICIAL_SFMT), 1)
ifeq ($(USE_OFFICIAL_SFMT_JUMP), 1)
OBJMAGI	:= $(patsubst %.c, $(OBJDIR)/%.ompmpi.gsl.smtj.hdf5.o, $(notdir $(MAGISRC)))
else
OBJMAGI	:= $(patsubst %.c, $(OBJDIR)/%.ompmpi.gsl.sfmt.hdf5.o, $(notdir $(MAGISRC)))
endif
else
OBJMAGI	:= $(patsubst %.c, $(OBJDIR)/%.ompmpi.gsl.hdf5.o,      $(notdir $(MAGISRC)))
endif
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,             $(notdir $(FILELIB) $(ALLCLIB)))
else
ifeq ($(USE_OFFICIAL_SFMT), 1)
ifeq ($(USE_OFFICIAL_SFMT_JUMP), 1)
OBJMAGI	:= $(patsubst %.c, $(OBJDIR)/%.ompmpi.gsl.smtj.o,      $(notdir $(MAGISRC)))
else
OBJMAGI	:= $(patsubst %.c, $(OBJDIR)/%.ompmpi.gsl.sfmt.o,      $(notdir $(MAGISRC)))
endif
else
OBJMAGI	:= $(patsubst %.c, $(OBJDIR)/%.ompmpi.gsl.o,           $(notdir $(MAGISRC)))
endif
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.o,                      $(notdir $(ALLCLIB)))
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,                  $(notdir $(FILELIB)))
endif
ifeq ($(USE_OFFICIAL_SFMT), 1)
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.omp.gsl.sfmt.o,         $(notdir $(DISKLIB)))
else
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.omp.gsl.o,              $(notdir $(DISKLIB)))
endif
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.omp.o,                  $(notdir $(MAGILIB)))
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJEDIT	:= $(patsubst %.c, $(OBJDIR)/%.ompmpi.hdf5.o, $(notdir $(EDITSRC)))
OBJEDIT	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,    $(notdir $(FILELIB) $(ALLCLIB)))
else
OBJEDIT	:= $(patsubst %.c, $(OBJDIR)/%.ompmpi.o,      $(notdir $(EDITSRC)))
OBJEDIT	+= $(patsubst %.c, $(OBJDIR)/%.o,             $(notdir $(ALLCLIB)))
OBJEDIT	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,         $(notdir $(FILELIB)))
endif
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJPENE	:= $(patsubst %.c, $(OBJDIR)/%.mpi.pl.hdf5.o, $(notdir $(PENESRC)))
OBJPENE	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,    $(notdir $(FILELIB) $(ALLCLIB)))
else
OBJPENE	:= $(patsubst %.c, $(OBJDIR)/%.mpi.pl.o,      $(notdir $(PENESRC)))
OBJPENE	+= $(patsubst %.c, $(OBJDIR)/%.o,             $(notdir $(ALLCLIB)))
OBJPENE	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,         $(notdir $(FILELIB)))
endif
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJDIST	:= $(patsubst %.c, $(OBJDIR)/%.mpi.gsl.pl.hdf5.o, $(notdir $(DISTSRC)))
OBJDIST	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,        $(notdir $(FILELIB) $(ALLCLIB)))
else
OBJDIST	:= $(patsubst %.c, $(OBJDIR)/%.mpi.gsl.pl.o,      $(notdir $(DISTSRC)))
OBJDIST	+= $(patsubst %.c, $(OBJDIR)/%.o,                 $(notdir $(ALLCLIB)))
OBJDIST	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,             $(notdir $(FILELIB)))
endif
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJPCDF	:= $(patsubst %.c, $(OBJDIR)/%.mpi.pl.hdf5.o, $(notdir $(PCDFSRC)))
OBJPCDF	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,    $(notdir $(FILELIB) $(ALLCLIB)))
else
OBJPCDF	:= $(patsubst %.c, $(OBJDIR)/%.mpi.pl.o,      $(notdir $(PCDFSRC)))
OBJPCDF	+= $(patsubst %.c, $(OBJDIR)/%.o,             $(notdir $(ALLCLIB)))
OBJPCDF	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,         $(notdir $(FILELIB)))
endif
OBJPCDF	+= $(patsubst %.c, $(OBJDIR)/%.o,             $(notdir $(PCDFLIB)))
#################################################################################################
OBJPELP	:= $(patsubst %.c, $(OBJDIR)/%.mpi.pl.o, $(notdir $(PELPSRC)))
OBJPELP	+= $(patsubst %.c, $(OBJDIR)/%.o,        $(notdir $(PCDFLIB)))
#################################################################################################
OBJPDEP	:= $(patsubst %.c, $(OBJDIR)/%.pl.o, $(notdir $(PDEPSRC)))
#################################################################################################
OBJPBRK	:= $(patsubst %.c, $(OBJDIR)/%.pl.o, $(notdir $(PBRKSRC)))
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJPCMP	:= $(patsubst %.c, $(OBJDIR)/%.mpi.pl.hdf5.o, $(notdir $(PCMPSRC)))
endif
#################################################################################################
OBJPFLP	:= $(patsubst %.c, $(OBJDIR)/%.pl.o, $(notdir $(PFLPSRC)))
#################################################################################################
OBJPRAD	:= $(patsubst %.c, $(OBJDIR)/%.pl.o, $(notdir $(PRADSRC)))
#################################################################################################
OBJPDF	:= $(patsubst %.c, $(OBJDIR)/%.mpi.pl.hdf5.o, $(notdir $(PDFSRC)))
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJPJET	:= $(patsubst %.c, $(OBJDIR)/%.mpi.pl.hdf5.o, $(notdir $(PJETSRC)))
OBJPJET	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,    $(notdir $(FILELIB)))
else
OBJPJET	:= $(patsubst %.c, $(OBJDIR)/%.pl.o,          $(notdir $(PJETSRC)))
OBJPJET	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,         $(notdir $(FILELIB)))
endif
OBJPJET	+= $(patsubst %.c, $(OBJDIR)/%.o,             $(notdir $(ALLCLIB)))
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJPDSK	:= $(patsubst %.c, $(OBJDIR)/%.mpi.pl.hdf5.o, $(notdir $(PDSKSRC)))
else
OBJPDSK	:= $(patsubst %.c, $(OBJDIR)/%.pl.o,          $(notdir $(PDSKSRC)))
endif
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJAACT	:= $(patsubst %.c, $(OBJDIR)/%.mpi.gsl.pl.hdf5.o, $(notdir $(AACTSRC)))
OBJAACT	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,        $(notdir $(FILELIB) $(ALLCLIB)))
else
OBJAACT	:= $(patsubst %.c, $(OBJDIR)/%.mpi.gsl.pl.o,      $(notdir $(PACTSRC)))
OBJAACT	+= $(patsubst %.c, $(OBJDIR)/%.o,		  $(notdir $(ALLCLIB)))
OBJAACT	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,		  $(notdir $(FILELIB)))
endif
OBJAACT	+= $(patsubst %.c, $(OBJDIR)/%.omp.o,		  $(notdir spline.c))
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJAERR	:= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o, $(notdir $(AERRSRC) $(FILELIB) $(ALLCLIB)))
else
OBJAERR	:= $(patsubst %.c, $(OBJDIR)/%.mpi.o,      $(notdir $(AERRSRC) $(FILELIB)))
OBJAERR	+= $(patsubst %.c, $(OBJDIR)/%.o,          $(notdir $(ALLCLIB)))
endif
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJAPRF	:= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o, $(notdir $(APRFSRC) $(FILELIB) $(ALLCLIB)))
else
OBJAPRF	:= $(patsubst %.c, $(OBJDIR)/%.mpi.o,      $(notdir $(APRFSRC) $(FILELIB)))
OBJAPRF	+= $(patsubst %.c, $(OBJDIR)/%.o,          $(notdir $(ALLCLIB)))
endif
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJAM31	:= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o, $(notdir $(AM31SRC) $(FILELIB) $(ALLCLIB)))
else
OBJAM31	:= $(patsubst %.c, $(OBJDIR)/%.mpi.o,      $(notdir $(AM31SRC) $(FILELIB)))
OBJAM31	+= $(patsubst %.c, $(OBJDIR)/%.o,          $(notdir $(ALLCLIB)))
endif
OBJAM31	+= $(patsubst %.c, $(OBJDIR)/%.o,          $(notdir $(M31LIB)))
#################################################################################################


#################################################################################################
## Rules
#################################################################################################
all:	TAGS $(GOTHIC) $(MAGI) $(EDITOR) $(PLTENE) $(PLTDST) $(ANALPRF)
#################################################################################################
.PHONY:	gothic init magi cold editor plot bench sample disk anal
gothic:	TAGS $(GOTHIC)
init:	TAGS $(MAGI) $(MKCOLD) $(EDITOR)
magi:	TAGS $(MAGI)
cold:	TAGS $(MKCOLD)
editor:	TAGS $(EDITOR)
plot:	TAGS $(PLTENE) $(PLTDST) $(PLTCDF)
bench:	TAGS $(OPTCFG) $(PLTELP) $(PLTDEP) $(PLTBRK) $(PLTCMP) $(PLTFLP) $(PLTRAD)
sample:	TAGS $(SAMPLE) $(PLTDF)
disk:	TAGS $(PLTJET) $(PLTDISK)
anal:	TAGS $(ANALACT) $(ANALERR) $(ANALPRF)
m31:	TAGS $(M31OBS)
sass:	TAGS $(GOTHIC).sass
#################################################################################################
## Making TAGS file for Emacs
TAGS:
	$(VERBOSE)find $(SRCDIR) -name "*.[ch]" -print | etags -
#################################################################################################
ifeq ($(USE_COOPERATIVE_GROUPS), 1)
$(OBJDIR)/make_dev.o:	make_dev.cu
	$(VERBOSE)$(CU) $(CUFLAG) $(CUARG) $(DEBUG) $(CUDBG) $(CUOPT)          $(CUINC)           $(PREC) $(UNIT) $(PROFILE) --device-c -o $@ -c $<
$(OBJDIR)/make_dev.mpi.o:	make_dev.cu
	$(VERBOSE)$(CU) $(CUFLAG) $(CUARG) $(DEBUG) $(CUDBG) $(CUOPT) $(CUMPI) $(CUINC) $(MPIINC) $(PREC) $(UNIT) $(PROFILE) --device-c -o $@ -c $<
endif
#################################################################################################
ifeq ($(USE_LIS_FOR_MAGI), 1)
$(OBJDIR)/magi.ompmpi.gsl.smtj.hdf5.o:	magi.c
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCARG) $(DEBUG) $(CCDBG) $(CCOPT) $(CCINC) $(OMPINC) $(GSLINC) $(LISINC) $(SMTJINC) $(HDF5INC) $(PREC) $(UNIT) $(PROFILE) -o $@ -c $<
$(OBJDIR)/magi.ompmpi.gsl.sfmt.hdf5.o:	magi.c
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCARG) $(DEBUG) $(CCDBG) $(CCOPT) $(CCINC) $(OMPINC) $(GSLINC) $(LISINC) $(SFMTINC) $(HDF5INC) $(PREC) $(UNIT) $(PROFILE) -o $@ -c $<
$(OBJDIR)/magi.ompmpi.gsl.hdf5.o:	magi.c
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCARG) $(DEBUG) $(CCDBG) $(CCOPT) $(CCINC) $(OMPINC) $(GSLINC) $(LISINC) $(HDF5INC) $(PREC) $(UNIT) $(PROFILE) -o $@ -c $<
$(OBJDIR)/magi.ompmpi.gsl.smtj.o:	magi.c
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCARG) $(DEBUG) $(CCDBG) $(CCOPT) $(CCINC) $(OMPINC) $(GSLINC) $(LISINC) $(SMTJINC) $(PREC) $(UNIT) $(PROFILE) -o $@ -c $<
$(OBJDIR)/magi.ompmpi.gsl.sfmt.o:	magi.c
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCARG) $(DEBUG) $(CCDBG) $(CCOPT) $(CCINC) $(OMPINC) $(GSLINC) $(LISINC) $(SFMTINC) $(PREC) $(UNIT) $(PROFILE) -o $@ -c $<
$(OBJDIR)/magi.ompmpi.gsl.o:	magi.c
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCARG) $(DEBUG) $(CCDBG) $(CCOPT) $(CCINC) $(OMPINC) $(GSLINC) $(LISINC) $(PREC) $(UNIT) $(PROFILE) -o $@ -c $<
$(OBJDIR)/potdens.omp.gsl.sfmt.o:	potdens.c
	$(VERBOSE)$(CC)    $(CCFLAG) $(CCARG) $(DEBUG) $(CCDBG) $(CCOPT) $(CCINC) $(OMPINC) $(GSLINC) $(LISINC) $(SFMTINC) $(PREC) $(UNIT) $(PROFILE) -o $@ -c $<
$(OBJDIR)/potdens.omp.gsl.o:	potdens.c
	$(VERBOSE)$(CC)    $(CCFLAG) $(CCARG) $(DEBUG) $(CCDBG) $(CCOPT) $(CCINC) $(OMPINC) $(GSLINC) $(LISINC) $(PREC) $(UNIT) $(PROFILE) -o $@ -c $<
endif
#################################################################################################
## CUDA/C code
# multiple GPUs code
ifeq ($(USE_COOPERATIVE_GROUPS), 1)
ifeq ($(FORCE_SINGLE_GPU_RUN), 1)
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(GOTHIC):	$(OBJMPGT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)cudalib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(CU) $(CUFLAG) $(CUARG) $(DEBUG) $(CUDBG) $(CUOPT) $(CUINC) $(PREC) $(UNIT) $(PROFILE) --device-link $(OBJMPGT) --output-file $(OBJDIR)/gothic_link.o
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMPGT) $(OBJDIR)/gothic_link.o -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)timer -l$(LIBPREC)cudalib -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(CULIB) -lcudadevrt $(CCLIB)
else
$(GOTHIC):	$(OBJMPGT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)cudalib.a
	$(VERBOSE)$(CU) $(CUFLAG) $(CUARG) $(DEBUG) $(CUDBG) $(CUOPT) $(CUINC) $(PREC) $(UNIT) $(PROFILE) --device-link $(OBJMPGT) --output-file $(OBJDIR)/gothic_link.o
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMPGT) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)timer -l$(LIBPREC)cudalib -l$(LIBPREC)mpilib $(CULIB) -lcudadevrt $(CCLIB)
endif
else
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(GOTHIC):	$(OBJMPGT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)cudalib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(CU) $(CUFLAG) $(CUARG) $(DEBUG) $(CUDBG) $(CUOPT) $(CUMPI) $(CUINC) $(MPIINC) $(PREC) $(UNIT) $(PROFILE) --device-link $(OBJMPGT) --output-file $(OBJDIR)/gothic_link.o
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMPGT) $(OBJDIR)/gothic_link.o -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)timer -l$(LIBPREC)cudalib -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(CULIB) -lcudadevrt $(CCLIB)
else
$(GOTHIC):	$(OBJMPGT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)cudalib.a
	$(VERBOSE)$(CU) $(CUFLAG) $(CUARG) $(DEBUG) $(CUDBG) $(CUOPT) $(CUMPI) $(CUINC) $(MPIINC) $(PREC) $(UNIT) $(PROFILE) --device-link $(OBJMPGT) --output-file $(OBJDIR)/gothic_link.o
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMPGT) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)timer -l$(LIBPREC)cudalib -l$(LIBPREC)mpilib $(CULIB) -lcudadevrt $(CCLIB)
endif
endif
else
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(GOTHIC):	$(OBJMPGT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)cudalib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMPGT) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)timer -l$(LIBPREC)cudalib -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(CULIB) $(CCLIB)
else
$(GOTHIC):	$(OBJMPGT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)cudalib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMPGT) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)timer -l$(LIBPREC)cudalib -l$(LIBPREC)mpilib $(CULIB) $(CCLIB)
endif
endif
#################################################################################################
## C code
$(OPTCFG):	$(OBJFIND)	$(MYLIB)/lib$(LIBPREC)myutil.a
	$(VERBOSE)$(CC)    $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJFIND) -L$(MYLIB) -l$(LIBPREC)myutil $(CCLIB)
$(SAMPLE):	$(OBJSMPL)
	$(VERBOSE)$(CC)    $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJSMPL) $(CCLIB)
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
ifeq ($(USE_OFFICIAL_SFMT), 1)
ifeq ($(USE_OFFICIAL_SFMT_JUMP), 1)
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a	$(MYLIB)/libsfmt$(SFMTPER).a	$(MYLIB)/libsmtj$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_sfmt$(SFMTPER) -l$(LIBPREC)timer -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(SMTJLIB) $(GSLLIB) $(CCLIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a	$(MYLIB)/libsfmt$(SFMTPER).a	$(MYLIB)/libsmtj$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_sfmt$(SFMTPER) -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)hdf5lib -l$(LIBPREC)ompmpilib $(HDF5LIB) $(SMTJLIB) $(GSLLIB) $(LISLIB) $(OMPLIB) $(CCLIB)
else
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a	$(MYLIB)/libsfmt$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_sfmt$(SFMTPER) -l$(LIBPREC)timer -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(SFMTLIB) $(GSLLIB) $(CCLIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a	$(MYLIB)/libsfmt$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_sfmt$(SFMTPER) -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)hdf5lib -l$(LIBPREC)ompmpilib $(HDF5LIB) $(SFMTLIB) $(GSLLIB) $(LISLIB) $(OMPLIB) $(CCLIB)
endif
else
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_gsl.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_gsl -l$(LIBPREC)timer -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(GSLLIB) $(CCLIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_gsl.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_gsl -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)hdf5lib -l$(LIBPREC)ompmpilib $(HDF5LIB) $(GSLLIB) $(LISLIB) $(OMPLIB) $(CCLIB)
endif
$(EDITOR):	$(OBJEDIT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJEDIT) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)hdf5lib -l$(LIBPREC)ompmpilib $(HDF5LIB) $(OMPLIB) $(CCLIB)
$(PLTENE):	$(OBJPENE)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPENE) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(CCLIB)
$(PLTDST):	$(OBJDIST)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJDIST) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(GSLLIB) $(CCLIB)
$(PLTCDF):	$(OBJPCDF)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPCDF) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(CCLIB)
$(ANALACT):	$(OBJAACT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJAACT) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(GSLLIB) $(CCLIB)
$(ANALERR):	$(OBJAERR)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJAERR) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(CCLIB)
$(ANALPRF):	$(OBJAPRF)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJAPRF) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(CCLIB)
$(M31OBS):	$(OBJAM31)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJAM31) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rotate -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(CCLIB)
else
ifeq ($(USE_OFFICIAL_SFMT), 1)
ifeq ($(USE_OFFICIAL_SFMT_JUMP), 1)
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/libsfmt$(SFMTPER).a	$(MYLIB)/libsmtj$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_sfmt$(SFMTPER) -l$(LIBPREC)timer -l$(LIBPREC)mpilib $(SMTJLIB) $(GSLLIB) $(CCLIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/libsfmt$(SFMTPER).a	$(MYLIB)/libsmtj$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_sfmt$(SFMTPER) -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)ompmpilib $(SMTJLIB) $(GSLLIB) $(LISLIB) $(OMPLIB) $(CCLIB)
else
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/libsfmt$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_sfmt$(SFMTPER) -l$(LIBPREC)timer -l$(LIBPREC)mpilib $(SFMTLIB) $(GSLLIB) $(CCLIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/libsfmt$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_sfmt$(SFMTPER) -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)ompmpilib $(SFMTLIB) $(GSLLIB) $(LISLIB) $(OMPLIB) $(CCLIB)
endif
else
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_gsl.a	$(MYLIB)/lib$(LIBPREC)timer.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)rand_gsl -l$(LIBPREC)timer $(GSLLIB) $(CCLIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_gsl.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rand_gsl -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)ompmpilib $(GSLLIB) $(LISLIB) $(OMPLIB) $(CCLIB)
endif
$(EDITOR):	$(OBJEDIT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJEDIT) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)ompmpilib $(OMPLIB) $(CCLIB)
$(PLTENE):	$(OBJPENE)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPENE) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib -l$(LIBPREC)mpilib $(CCLIB)
$(PLTDST):	$(OBJDIST)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJDIST) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib -l$(LIBPREC)mpilib $(GSLLIB) $(CCLIB)
$(PLTCDF):	$(OBJPCDF)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPCDF) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib -l$(LIBPREC)mpilib $(CCLIB)
$(ANALACT):	$(OBJAACT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJAACT) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib $(GSLLIB) $(OMPLIB) $(CCLIB)
$(ANALERR):	$(OBJAERR)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJAERR) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib $(CCLIB)
$(ANALPRF):	$(OBJAPRF)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJAPRF) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib $(CCLIB)
$(M31OBS):	$(OBJAM31)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rotate.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJAM31) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)rotate -l$(LIBPREC)mpilib $(CCLIB)
endif
$(PLTELP):	$(OBJPELP)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)mpilib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPELP) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)plplotlib -l$(LIBPREC)mpilib $(CCLIB)
$(PLTDEP):	$(OBJPDEP)	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPDEP) -L$(MYLIB) -l$(LIBPREC)plplotlib $(CCLIB)
$(PLTBRK):	$(OBJPBRK)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPBRK) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)plplotlib $(CCLIB)
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(PLTCMP):	$(OBJPCMP)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPCMP) -L$(MYLIB) -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib $(HDF5LIB) $(CCLIB)
endif
$(PLTFLP):	$(OBJPFLP)	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPFLP) -L$(MYLIB) -l$(LIBPREC)plplotlib $(CCLIB)
$(PLTRAD):	$(OBJPRAD)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPRAD) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib $(CCLIB)
$(PLTDF):	$(OBJPDF)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a	$(MYLIB)/lib$(LIBPREC)constants.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPDF) -L$(MYLIB) -l$(LIBPREC)constants -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib $(HDF5LIB) $(CCLIB)
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(PLTJET):	$(OBJPJET)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPJET) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib -l$(LIBPREC)mpilib $(HDF5LIB) $(CCLIB)
else
$(PLTJET):	$(OBJPJET)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPJET) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib $(CCLIB)
endif
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(PLTDISK):	$(OBJPDSK)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPDSK) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib $(HDF5LIB) $(CCLIB)
else
$(PLTDISK):	$(OBJPDSK)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPDSK) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)plplotlib $(CCLIB)
endif
#################################################################################################
# sass file
$(GOTHIC).sass:	$(GOTHIC)
	$(VERBOSE)$(CUASM) $< > $@
#################################################################################################
# original version
check-syntax:
	$(CC) $(DEBUG) $(CCFLAG) $(CCARG) -fsyntax-only $(CCINC) $(OMPINC) $(MPIINC) $(GSLINC) $(SFMTINC) $(SMTJINC) $(PLINC) -I/$(HOME)/tau/include $(HDF5INC) $(LISINC) $(PREC) $(UNIT) $(CHK_SOURCES)
check-syntax-cu:
	$(CU) $(DEBUG) $(CUARG) -Xcompiler -fsyntax-only $(CCINC) $(CUINC) $(OMPINC) $(MPIINC) $(GSLINC) $(SFMTINC) $(SMTJINC) $(PLINC) -I/$(HOME)/tau/include $(HDF5INC) $(PREC) $(UNIT) -o /dev/null -c $(CHK_SOURCES)
#################################################################################################
# depend:
# 	makedepend -- $(CCINC) $(CUINC) $(OMPINC) $(MPIINC) $(GSLINC) $(PLINC) -I/$(HOME)/tau/include $(HDF5INC) -- $(CC_FILE) $(CU_FILE)
#################################################################################################
dir:
	$(VERBOSE)mkdir -p $(BINDIR)
	$(VERBOSE)mkdir -p $(OBJDIR)
	$(VERBOSE)mkdir -p $(ASMDIR)
	$(VERBOSE)mkdir -p $(PTXDIR)
	$(VERBOSE)mkdir -p dat doc fig log bench
#################################################################################################
zip:
	if [ ! -d pub ]; then \
	$(VERBOSE)mkdir pub; \
	fi
	$(VERBOSE)tar cvJf $(DATE)tree.tar.xz \
	README.md LICENSE.txt Makefile $(SRCDIR) sh cfg plt py plugins/README \
	$(XMLBODY) $(SRCBODY) $(XMLSNAP) $(SRCSNAP) $(XMLAERR) $(SRCAERR) $(XMLDUMP) $(SRCDUMP) $(XMLDISK) $(SRCDISK) $(XMLDIST) $(SRCDIST) $(XMLPROF) $(SRCPROF)
	$(VERBOSE)mv       $(DATE)tree.tar.xz pub/
#################################################################################################
clean:
	$(VERBOSE)rm -f TAGS gmon.out
	$(VERBOSE)rm -f $(ASMDIR)/*.s
	$(VERBOSE)rm -f $(BINDIR)/*.log
	$(VERBOSE)rm -f $(OBJDIR)/*.log
	$(VERBOSE)rm -f $(OBJDIR)/*.o
	$(VERBOSE)rm -f $(GOTHIC) $(OBJMPGT) $(GOTHIC).sass
	$(VERBOSE)rm -f $(MKCOLD) $(MAGI) $(EDITOR)
ifneq ($(PLPLOT), 0)
	$(VERBOSE)rm -f $(PLTENE) $(PLTDST) $(PLTCDF)
	$(VERBOSE)rm -f $(PLTELP) $(PLTDEP) $(PLTBRK) $(PLTCMP) $(PLTFLP)
	$(VERBOSE)rm -f $(PLTRAD) $(PLTJET) $(PLTDF) $(PLTDISK)
endif
	$(VERBOSE)rm -f $(OPTCFG)
	$(VERBOSE)rm -f $(ANALACT) $(ANALERR) $(ANALPRF)
	$(VERBOSE)rm -f $(M31OBS)
	$(VERBOSE)rm -f $(SAMPLE)
#################################################################################################
visit:	$(DIRBODY)/Makefile $(DIRSNAP)/Makefile $(DIRPLOT)/Makefile $(DIRPM31)/Makefile $(DIRAERR)/Makefile $(DIRDUMP)/Makefile $(DIRDISK)/Makefile $(DIRDIST)/Makefile $(DIRPROF)/Makefile
#################################################################################################
$(DIRBODY)/Makefile:	$(XMLBODY) $(SRCBODY)
	$(VERBOSE)cd $(DIRBODY);\
	$(VERBOSE)xml2info $(TAGBODY).xml;\
	$(VERBOSE)xml2cmake $(TAGBODY).xml;\
	$(VERBOSE)cmake -DCMAKE_BUILD_TYPE:STRING=Debug;\
	$(VERBOSE)make
	$(VERBOSE)touch $@
$(DIRSNAP)/Makefile:	$(XMLSNAP) $(SRCSNAP)
	$(VERBOSE)cd $(DIRSNAP);\
	$(VERBOSE)xml2info $(TAGSNAP).xml;\
	$(VERBOSE)xml2cmake $(TAGSNAP).xml;\
	$(VERBOSE)cmake -DCMAKE_BUILD_TYPE:STRING=Debug;\
	$(VERBOSE)make
	$(VERBOSE)touch $@
$(DIRPLOT)/Makefile:	$(XMLPLOT) $(SRCPLOT)
	$(VERBOSE)cd $(DIRPLOT);\
	$(VERBOSE)xml2info $(TAGPLOT).xml;\
	$(VERBOSE)xml2cmake $(TAGPLOT).xml;\
	$(VERBOSE)cmake -DCMAKE_BUILD_TYPE:STRING=Debug;\
	$(VERBOSE)make
	$(VERBOSE)touch $@
$(DIRPM31)/Makefile:	$(XMLPM31) $(SRCPM31)
	$(VERBOSE)cd $(DIRPM31);\
	$(VERBOSE)xml2info $(TAGPM31).xml;\
	$(VERBOSE)xml2cmake $(TAGPM31).xml;\
	$(VERBOSE)cmake -DCMAKE_BUILD_TYPE:STRING=Debug;\
	$(VERBOSE)make
	$(VERBOSE)touch $@
$(DIRAERR)/Makefile:	$(XMLAERR) $(SRCAERR)
	$(VERBOSE)cd $(DIRAERR);\
	$(VERBOSE)xml2info $(TAGAERR).xml;\
	$(VERBOSE)xml2cmake $(TAGAERR).xml;\
	$(VERBOSE)cmake -DCMAKE_BUILD_TYPE:STRING=Debug;\
	$(VERBOSE)make
	$(VERBOSE)touch $@
$(DIRDUMP)/Makefile:	$(XMLDUMP) $(SRCDUMP)
	$(VERBOSE)cd $(DIRDUMP);\
	$(VERBOSE)xml2info $(TAGDUMP).xml;\
	$(VERBOSE)xml2cmake $(TAGDUMP).xml;\
	$(VERBOSE)cmake -DCMAKE_BUILD_TYPE:STRING=Debug;\
	$(VERBOSE)make
	$(VERBOSE)touch $@
$(DIRDISK)/Makefile:	$(XMLDISK) $(SRCDISK)
	$(VERBOSE)cd $(DIRDISK);\
	$(VERBOSE)xml2info $(TAGDISK).xml;\
	$(VERBOSE)xml2cmake $(TAGDISK).xml;\
	$(VERBOSE)cmake -DCMAKE_BUILD_TYPE:STRING=Debug;\
	$(VERBOSE)make
	$(VERBOSE)touch $@
$(DIRDIST)/Makefile:	$(XMLDIST) $(SRCDIST)
	$(VERBOSE)cd $(DIRDIST);\
	$(VERBOSE)xml2info $(TAGDIST).xml;\
	$(VERBOSE)xml2cmake $(TAGDIST).xml;\
	$(VERBOSE)cmake -DCMAKE_BUILD_TYPE:STRING=Debug;\
	$(VERBOSE)make
	$(VERBOSE)touch $@
$(DIRPROF)/Makefile:	$(XMLPROF) $(SRCPROF)
	$(VERBOSE)cd $(DIRPROF);\
	$(VERBOSE)xml2info $(TAGPROF).xml;\
	$(VERBOSE)xml2cmake $(TAGPROF).xml;\
	$(VERBOSE)cmake -DCMAKE_BUILD_TYPE:STRING=Debug;\
	$(VERBOSE)make
	$(VERBOSE)touch $@
#################################################################################################
rmvisit:
	$(VERBOSE)rm  -f $(DIRBODY)/CMakeCache.txt $(DIRBODY)/CMakeLists.txt $(DIRBODY)/Makefile $(DIRBODY)/cmake_install.cmake $(DIRBODY)/make.log
	$(VERBOSE)rm -rf $(DIRBODY)/CMakeFiles
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libE$(TAGBODY)Database_???.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libI$(TAGBODY)Database.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libM$(TAGBODY)Database.so
	$(VERBOSE)rm  -f $(DIRSNAP)/CMakeCache.txt $(DIRSNAP)/CMakeLists.txt $(DIRSNAP)/Makefile $(DIRSNAP)/cmake_install.cmake $(DIRSNAP)/make.log
	$(VERBOSE)rm -rf $(DIRSNAP)/CMakeFiles
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libE$(TAGSNAP)Database_???.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libI$(TAGSNAP)Database.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libM$(TAGSNAP)Database.so
	$(VERBOSE)rm  -f $(DIRPLOT)/CMakeCache.txt $(DIRPLOT)/CMakeLists.txt $(DIRPLOT)/Makefile $(DIRPLOT)/cmake_install.cmake $(DIRPLOT)/make.log
	$(VERBOSE)rm -rf $(DIRPLOT)/CMakeFiles
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libE$(TAGPLOT)Database_???.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libI$(TAGPLOT)Database.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libM$(TAGPLOT)Database.so
	$(VERBOSE)rm  -f $(DIRPM31)/CMakeCache.txt $(DIRPM31)/CMakeLists.txt $(DIRPM31)/Makefile $(DIRPM31)/cmake_install.cmake $(DIRPM31)/make.log
	$(VERBOSE)rm -rf $(DIRPM31)/CMakeFiles
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libE$(TAGPM31)Database_???.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libI$(TAGPM31)Database.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libM$(TAGPM31)Database.so
	$(VERBOSE)rm  -f $(DIRAERR)/CMakeCache.txt $(DIRAERR)/CMakeLists.txt $(DIRAERR)/Makefile $(DIRAERR)/cmake_install.cmake $(DIRAERR)/make.log
	$(VERBOSE)rm -rf $(DIRAERR)/CMakeFiles
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libE$(TAGAERR)Database_???.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libI$(TAGAERR)Database.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libM$(TAGAERR)Database.so
	$(VERBOSE)rm  -f $(DIRDUMP)/CMakeCache.txt $(DIRDUMP)/CMakeLists.txt $(DIRDUMP)/Makefile $(DIRDUMP)/cmake_install.cmake $(DIRDUMP)/make.log
	$(VERBOSE)rm -rf $(DIRDUMP)/CMakeFiles
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libE$(TAGDUMP)Database_???.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libI$(TAGDUMP)Database.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libM$(TAGDUMP)Database.so
	$(VERBOSE)rm  -f $(DIRDISK)/CMakeCache.txt $(DIRDISK)/CMakeLists.txt $(DIRDISK)/Makefile $(DIRDISK)/cmake_install.cmake $(DIRDISK)/make.log
	$(VERBOSE)rm -rf $(DIRDISK)/CMakeFiles
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libE$(TAGDISK)Database_???.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libI$(TAGDISK)Database.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libM$(TAGDISK)Database.so
	$(VERBOSE)rm  -f $(DIRDIST)/CMakeCache.txt $(DIRDIST)/CMakeLists.txt $(DIRDIST)/Makefile $(DIRDIST)/cmake_install.cmake $(DIRDIST)/make.log
	$(VERBOSE)rm -rf $(DIRDIST)/CMakeFiles
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libE$(TAGDIST)Database_???.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libI$(TAGDIST)Database.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libM$(TAGDIST)Database.so
	$(VERBOSE)rm  -f $(DIRPROF)/CMakeCache.txt $(DIRPROF)/CMakeLists.txt $(DIRPROF)/Makefile $(DIRPROF)/cmake_install.cmake $(DIRPROF)/make.log
	$(VERBOSE)rm -rf $(DIRPROF)/CMakeFiles
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libE$(TAGPROF)Database_???.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libI$(TAGPROF)Database.so
	$(VERBOSE)rm  -f $(HOME)/.visit/2.*/linux-x86_64/plugins/databases/libM$(TAGPROF)Database.so
#################################################################################################


#################################################################################################
## Dependency
#################################################################################################
COMMON_DEP	:=	Makefile	$(MYINC)/common.mk	$(MYINC)/macro.h
#################################################################################################
## $(MAINDIR)/*
GOTHIC_DEP	:=	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/name.h	$(MYINC)/constants.h	$(MYINC)/timer.h	$(MYINC)/cudalib.h
GOTHIC_DEP	+=	$(MISCDIR)/device.h	$(MISCDIR)/benchmark.h	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(MISCDIR)/allocate_dev.h	$(MISCDIR)/convert.h	$(MISCDIR)/tune.h
GOTHIC_DEP	+=	$(FILEDIR)/io.h
GOTHIC_DEP	+=	$(SORTDIR)/peano.h
GOTHIC_DEP	+=	$(TREEDIR)/make.h	$(TREEDIR)/let.h	$(TREEDIR)/buf_inc.h	$(TREEDIR)/make_dev.h	$(TREEDIR)/walk_dev.h
GOTHIC_DEP	+=	$(TIMEDIR)/adv_dev.h
ifeq ($(FORCE_SINGLE_GPU_RUN), 0)
GOTHIC_DEP	+=	$(MYINC)/mpilib.h	$(PARADIR)/mpicfg.h	$(PARADIR)/exchange.h
GOTHIC_DEP	+=	$(TREEDIR)/geo_dev.h	$(TREEDIR)/let_dev.h	$(PARADIR)/exchange_dev.h
ifeq ($(ENCLOSING_BALL_FOR_LET), 1)
GOTHIC_DEP	+=	$(TREEDIR)/icom_dev.h
endif
endif
ifeq ($(SET_EXTERNAL_FIELD), 1)
GOTHIC_DEP	+=	$(TREEDIR)/potential_dev.h
endif
GOTHIC_DEP	+=	$(MISCDIR)/brent.h
GOTHIC_DEP	+=	$(SORTDIR)/peano_dev.h
GOTHIC_DEP	+=	$(TREEDIR)/neighbor_dev.h
GOTHIC_DEP	+=	$(TREEDIR)/shrink_dev.h
$(OBJDIR)/gothic.mpi.hdf5.o:	$(GOTHIC_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/gothic.mpi.o:		$(GOTHIC_DEP)
$(OBJDIR)/showOptConfig.o:	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/name.h
#################################################################################################
## $(MISCDIR)/*
$(OBJDIR)/allocate.mpi.hdf5.o:	$(COMMON_DEP)	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h
$(OBJDIR)/allocate.mpi.o:	$(COMMON_DEP)	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h
$(OBJDIR)/allocate_dev.o:	$(COMMON_DEP)	$(MISCDIR)/structure.h	$(MISCDIR)/allocate_dev.h	$(MYINC)/cudalib.h	$(TREEDIR)/walk_dev.h
$(OBJDIR)/convert.mpi.hdf5.o:	$(COMMON_DEP)	$(MISCDIR)/structure.h	$(MISCDIR)/convert.h
$(OBJDIR)/convert.o:		$(COMMON_DEP)	$(MISCDIR)/structure.h	$(MISCDIR)/convert.h
$(OBJDIR)/brent.o:		$(COMMON_DEP)	$(MISCDIR)/brent.h
$(OBJDIR)/tune.o:		$(COMMON_DEP)	$(MISCDIR)/tune.h
#################################################################################################
## $(FILEDIR)/*
IO_DEP	:=	$(COMMON_DEP)	$(MYINC)/constants.h	$(MYINC)/name.h	$(MYINC)/mpilib.h
IO_DEP	+=	$(MISCDIR)/benchmark.h	$(MISCDIR)/structure.h	$(MISCDIR)/tune.h	$(MISCDIR)/brent.h	$(FILEDIR)/io.h
ifeq ($(HUNT_OPTIMAL_WALK_TREE), 1)
IO_DEP	+=	$(TREEDIR)/walk_dev.h
IO_DEP	+=	$(TREEDIR)/neighbor_dev.h
IO_DEP	+=	$(TREEDIR)/shrink_dev.h
endif
ifeq ($(HUNT_OPTIMAL_MAKE_TREE), 1)
IO_DEP	+=	$(SORTDIR)/peano_dev.h	$(TREEDIR)/make_dev.h
endif
ifeq ($(HUNT_OPTIMAL_MAKE_NODE), 1)
IO_DEP	+=	$(TREEDIR)/make_dev.h
endif
ifeq ($(HUNT_OPTIMAL_INTEGRATE), 1)
IO_DEP	+=	$(TIMEDIR)/adv_dev.h
endif
ifeq ($(HUNT_OPTIMAL_NEIGHBOUR), 1)
IO_DEP	+=	$(TREEDIR)/neighbor_dev.h
endif
ifeq ($(SET_EXTERNAL_FIELD), 1)
IO_DEP	+=	$(INITDIR)/external.h
endif
$(OBJDIR)/io.mpi.hdf5.o:	$(IO_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/io.mpi.o:		$(IO_DEP)
#################################################################################################
## $(SORTDIR)/*
PEANO_DEP	:=	$(COMMON_DEP)	$(MISCDIR)/benchmark.h	$(MISCDIR)/structure.h	$(SORTDIR)/peano.h	$(TREEDIR)/make.h
$(OBJDIR)/peano.o:	$(PEANO_DEP)
$(OBJDIR)/peano_dev.o:	$(PEANO_DEP)	$(MYINC)/cudalib.h	$(UTILDIR)/gsync_dev.cu	$(SORTDIR)/peano_dev.h	$(PARADIR)/exchange_dev.h	$(UTILDIR)/compare_vec4_inc.cu	$(UTILDIR)/compare_vec4_inc.cuh	$(UTILDIR)/compare_vec4_del.cuh
#################################################################################################
## $(TIMEDIR)/*
ADV_DEV_DEP	:=	$(COMMON_DEP)	$(MYINC)/timer.h	$(MYINC)/cudalib.h
ADV_DEV_DEP	+=	$(MISCDIR)/structure.h	$(MISCDIR)/benchmark.h	$(MISCDIR)/device.h	$(TIMEDIR)/adv_dev.h	$(TREEDIR)/walk_dev.h
$(OBJDIR)/adv_dev.mpi.o:	$(ADV_DEV_DEP)	$(MYINC)/mpilib.h	$(PARADIR)/mpicfg.h
$(OBJDIR)/adv_dev.o:		$(ADV_DEV_DEP)
#################################################################################################
## $(TREEDIR)/*
TREE_DEP	:=	$(COMMON_DEP)	$(MISCDIR)/benchmark.h	$(MISCDIR)/structure.h	$(TREEDIR)/make.h
$(OBJDIR)/geo_dev.o:	$(TREE_DEP)	$(MYINC)/cudalib.h	$(TREEDIR)/walk_dev.h	$(TREEDIR)/geo_dev.h
$(OBJDIR)/let.mpi.o:	$(TREE_DEP)	$(MYINC)/mpilib.h	$(TREEDIR)/let.h	$(PARADIR)/mpicfg.h
TREE_DEV_DEP	:=	$(TREE_DEP)	$(MYINC)/cudalib.h	$(MISCDIR)/device.h
TREE_BUF_DEP	:=	$(TREEDIR)/buf_inc.h	$(TREEDIR)/buf_inc.cu
TREE_LET_DEP	:=	$(MYINC)/mpilib.h	$(TREEDIR)/let.h	$(TREEDIR)/let_dev.h
$(OBJDIR)/let_dev.mpi.o:	$(TREE_DEV_DEP)	$(TREE_BUF_DEP)	$(TREE_LET_DEP)	$(MYINC)/timer.h	$(TREEDIR)/walk_dev.h
MAKE_DEV_DEP	:=	$(UTILDIR)/gsync_dev.cu	$(TREEDIR)/make_dev.h	$(TREEDIR)/let.h
MAKE_DEV_DEP	+=	$(UTILDIR)/scan_inc.cu	$(UTILDIR)/scan_inc.cuh	$(UTILDIR)/scan_del.cuh
MAKE_DEV_DEP	+=	$(UTILDIR)/scan_tsub_inc.cu	$(UTILDIR)/scan_tsub_inc.cuh	$(UTILDIR)/scan_tsub_del.cuh
MAKE_DEV_DEP	+=	$(UTILDIR)/compare_tsub_inc.cu	$(UTILDIR)/compare_tsub_inc.cuh	$(UTILDIR)/compare_tsub_del.cuh
$(OBJDIR)/make_dev.mpi.o:	$(TREE_DEV_DEP)	$(SORTDIR)/peano.h	$(SORTDIR)/peano_dev.h	$(MAKE_DEV_DEP)	$(MYINC)/timer.h	$(MISCDIR)/device.h
$(OBJDIR)/make_dev.o:		$(TREE_DEV_DEP)	$(SORTDIR)/peano.h	$(SORTDIR)/peano_dev.h	$(MAKE_DEV_DEP)
$(OBJDIR)/make.o:		$(TREE_DEP)	$(SORTDIR)/peano.h
TREE_RSORT_DEP	:=	$(UTILDIR)/gsync_dev.cu
NEIGHBOR_DEV_DEP	:=	$(TREE_RSORT_DEP)	$(MYINC)/cudalib.h	$(TREEDIR)/make_dev.h	$(TREEDIR)/neighbor_dev.h
$(OBJDIR)/neighbor_dev.mpi.o:	$(TREE_DEP)	$(NEIGHBOR_DEV_DEP)
$(OBJDIR)/neighbor_dev.o:	$(TREE_DEP)	$(NEIGHBOR_DEV_DEP)
SHRINK_DEV_DEP	:=	$(MYINC)/cudalib.h	$(MISCDIR)/device.h	$(TREEDIR)/walk_dev.h	$(TREEDIR)/make_dev.h	$(TREEDIR)/shrink_dev.h
SHRINK_DEV_DEP	+=	$(MISCDIR)/brent.h
SHRINK_DEV_DEP	+=	$(TREEDIR)/neighbor_dev.h
$(OBJDIR)/shrink_dev.mpi.o:	$(TREE_DEP)	$(SHRINK_DEV_DEP)
$(OBJDIR)/shrink_dev.o:		$(TREE_DEP)	$(SHRINK_DEV_DEP)
WALK_DEV_DEP	:=	$(TREE_DEV_DEP)	$(TREE_BUF_DEP)	$(MYINC)/timer.h	$(TREEDIR)/walk_dev.h	$(TREEDIR)/seb_dev.cu
ifeq ($(SET_EXTERNAL_FIELD), 1)
WALK_DEV_DEP	+=	$(TREEDIR)/potential_dev.h
endif
$(OBJDIR)/walk_dev.mpi.o:	$(WALK_DEV_DEP)	$(TREE_LET_DEP)	$(PARADIR)/mpicfg.h	$(MISCDIR)/tune.h
$(OBJDIR)/walk_dev.o:		$(WALK_DEV_DEP)
ICOM_DEV_DEP	:=	$(TREE_DEV_DEP)	$(TREEDIR)/make_dev.h	$(SORTDIR)/peano_dev.h	$(TREEDIR)/icom_dev.h
ICOM_DEV_DEP	+=	$(UTILDIR)/compare_inc.cu	$(UTILDIR)/compare_inc.cuh	$(UTILDIR)/compare_del.cuh
ICOM_DEV_DEP	+=	$(UTILDIR)/scan_inc.cu	$(UTILDIR)/scan_inc.cuh	$(UTILDIR)/scan_del.cuh
$(OBJDIR)/icom_dev.o:	$(ICOM_DEV_DEP)
$(OBJDIR)/potential_dev.mpi.hdf5.o:	$(TREE_DEP)	$(MYINC)/name.h	$(MYINC)/cudalib.h	$(TREEDIR)/walk_dev.h	$(TREEDIR)/potential_dev.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/potential_dev.mpi.o:		$(TREE_DEP)	$(MYINC)/name.h	$(MYINC)/cudalib.h	$(TREEDIR)/walk_dev.h	$(TREEDIR)/potential_dev.h
#################################################################################################
## $(INITDIR)/*
$(OBJDIR)/sample.o:	$(COMMON_DEP)	$(MYINC)/name.h
$(OBJDIR)/blas.omp.o:	$(COMMON_DEP)	$(INITDIR)/blas.h
$(OBJDIR)/spline.omp.o:	$(COMMON_DEP)	$(INITDIR)/spline.h
SPLINE_DEP	:=	$(COMMON_DEP)	$(INITDIR)/profile.h	$(INITDIR)/table.h	$(INITDIR)/spline.h
$(OBJDIR)/abel.omp.o:	$(SPLINE_DEP)	$(MYINC)/name.h	$(INITDIR)/magi.h	$(INITDIR)/abel.h
$(OBJDIR)/table.omp.o:	$(SPLINE_DEP)	$(MYINC)/name.h
ifeq ($(SET_EXTERNAL_FIELD_DISK), 1)
$(OBJDIR)/external.omp.o:	$(SPLINE_DEP)	$(INITDIR)/external.h	$(INITDIR)/potdens.h
else
$(OBJDIR)/external.omp.o:	$(SPLINE_DEP)	$(INITDIR)/external.h
endif
PROFILE_DEP	:=	$(COMMON_DEP)	$(MYINC)/constants.h	$(INITDIR)/profile.h	$(INITDIR)/magi.h
$(OBJDIR)/eddington.omp.o:	$(PROFILE_DEP)	$(INITDIR)/eddington.h
$(OBJDIR)/king.omp.o:		$(PROFILE_DEP)	$(INITDIR)/king.h
$(OBJDIR)/profile.omp.o:	$(PROFILE_DEP)	$(MYINC)/name.h
DISK_DEP	:=	$(PROFILE_DEP)	$(INITDIR)/magi.h	$(INITDIR)/spline.h	$(INITDIR)/blas.h	$(INITDIR)/potdens.h
$(OBJDIR)/potdens.omp.gsl.o:	$(DISK_DEP)	$(INITDIR)/abel.h
$(OBJDIR)/potdens.omp.gsl.sfmt.o:	$(DISK_DEP)	$(INITDIR)/abel.h
$(OBJDIR)/diskDF.omp.gsl.sfmt.o:	$(DISK_DEP)	$(MYINC)/rand.h	$(MISCDIR)/structure.h	$(INITDIR)/diskDF.h
$(OBJDIR)/diskDF.omp.gsl.o:	$(DISK_DEP)	$(MYINC)/rand.h	$(MISCDIR)/structure.h	$(INITDIR)/diskDF.h
IOFILE_DEP	:=	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/constants.h	$(MYINC)/timer.h	$(MYINC)/name.h
IOFILE_DEP	+=	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(MISCDIR)/tune.h	$(MISCDIR)/brent.h	$(FILEDIR)/io.h
$(OBJDIR)/uniformsphere.mpi.smtj.gsl.hdf5.o:	$(IOFILE_DEP)	$(MYINC)/rand.h	$(MYINC)/hdf5lib.h	$(MYINC)/sfmtjump_polynomial.h
$(OBJDIR)/uniformsphere.mpi.sfmt.gsl.hdf5.o:	$(IOFILE_DEP)	$(MYINC)/rand.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/uniformsphere.mpi.gsl.hdf5.o:		$(IOFILE_DEP)	$(MYINC)/rand.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/uniformsphere.sfmt.o:			$(IOFILE_DEP)	$(MYINC)/rand.h
$(OBJDIR)/uniformsphere.gsl.o:			$(IOFILE_DEP)	$(MYINC)/rand.h
MAGI_DEP	:=	$(MYINC)/rotate.h	$(MYINC)/rand.h	$(INITDIR)/magi.h	$(INITDIR)/king.h	$(INITDIR)/profile.h	$(INITDIR)/eddington.h
MAGI_DEP	+=	$(INITDIR)/table.h	$(INITDIR)/abel.h	$(INITDIR)/potdens.h	$(INITDIR)/diskDF.h
ifeq ($(SET_EXTERNAL_FIELD), 1)
MAGI_DEP	+=	$(INITDIR)/external.h
endif
$(OBJDIR)/magi.ompmpi.gsl.smtj.hdf5.o:	$(IOFILE_DEP)	$(MAGI_DEP)	$(MYINC)/hdf5lib.h	$(MYINC)/sfmtjump_polynomial.h
$(OBJDIR)/magi.ompmpi.gsl.sfmt.hdf5.o:	$(IOFILE_DEP)	$(MAGI_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/magi.ompmpi.gsl.hdf5.o:	$(IOFILE_DEP)	$(MAGI_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/magi.ompmpi.gsl.smtj.o:	$(IOFILE_DEP)	$(MAGI_DEP)	$(MYINC)/sfmtjump_polynomial.h
$(OBJDIR)/magi.ompmpi.gsl.sfmt.o:	$(IOFILE_DEP)	$(MAGI_DEP)
$(OBJDIR)/magi.ompmpi.gsl.o:		$(IOFILE_DEP)	$(MAGI_DEP)
$(OBJDIR)/editor.ompmpi.hdf5.o:	$(IOFILE_DEP)	$(MYINC)/rotate.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/editor.ompmpi.o:	$(IOFILE_DEP)	$(MYINC)/rotate.h
#################################################################################################
## $(PLOTDIR)/*
PLOT_DEP	:=	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/constants.h	\
			$(MYINC)/plplotlib.h	$(MYINC)/mpilib.h	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(FILEDIR)/io.h
$(OBJDIR)/cdf.mpi.pl.hdf5.o:	$(PLOT_DEP)	$(MYINC)/name.h	$(PLOTDIR)/cdflib.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/cdf.mpi.pl.o:		$(PLOT_DEP)	$(MYINC)/name.h	$(PLOTDIR)/cdflib.h
$(OBJDIR)/distribution.mpi.gsl.pl.hdf5.o:	$(PLOT_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/distribution.mpi.gsl.pl.o:		$(PLOT_DEP)
$(OBJDIR)/energy.mpi.pl.hdf5.o:	$(PLOT_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/energy.mpi.pl.o:	$(PLOT_DEP)
$(OBJDIR)/ndep.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h
$(OBJDIR)/ball.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/myutil.h	$(MYINC)/name.h	$(MYINC)/constants.h
$(OBJDIR)/breakdown.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/myutil.h	$(MYINC)/name.h	$(MISCDIR)/benchmark.h
$(OBJDIR)/elapsed.mpi.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/myutil.h	$(MYINC)/name.h	$(MYINC)/mpilib.h	$(PLOTDIR)/cdflib.h
$(OBJDIR)/comparison.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/df.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h	$(MYINC)/constants.h	$(INITDIR)/magi.h	$(INITDIR)/eddington.h
$(OBJDIR)/performance.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h
$(OBJDIR)/needle.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h	$(MYINC)/constants.h	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(FILEDIR)/io.h
$(OBJDIR)/disk.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h	$(MYINC)/constants.h
$(OBJDIR)/cdflib.o:	$(COMMON_DEP)	$(PLOTDIR)/cdflib.h
#################################################################################################
## $(PARADIR)/*
PARA_DEP	:=	$(COMMON_DEP)	$(MYINC)/name.h	$(MYINC)/mpilib.h
PARA_DEP	+=	$(MISCDIR)/structure.h	$(PARADIR)/mpicfg.h
EXCG_DEP	:=	$(PARA_DEP)	$(TIMEDIR)/adv_dev.h	$(PARADIR)/exchange.h	$(MISCDIR)/tune.h	$(MISCDIR)/brent.h
$(OBJDIR)/exchange.mpi.o:	$(PARA_DEP)	$(PARADIR)/exchange.h
$(OBJDIR)/exchange_dev.mpi.o:	$(EXCG_DEP)	$(MYINC)/cudalib.h	$(MYINC)/timer.h	$(MISCDIR)/benchmark.h	$(UTILDIR)/gsync_dev.cu	$(SORTDIR)/peano_dev.h	$(PARADIR)/exchange_dev.h
$(OBJDIR)/mpicfg.mpi.o:		$(PARA_DEP)	$(TREEDIR)/make.h
#################################################################################################
# ## $(ANALDIR)/*
ANAL_DEP	:=	$(COMMON_DEP)	$(MYINC)/constants.h	$(MYINC)/name.h	$(MYINC)/mpilib.h	$(MYINC)/myutil.h
ANAL_DEP	+=	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(FILEDIR)/io.h
M31_DEP	:=	$(ANALDIR)/m31coord.h
$(OBJDIR)/action.mpi.gsl.pl.hdf5.o:	$(ANAL_DEP)	$(INITDIR)/spline.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/action.mpi.gsl.pl.o:		$(ANAL_DEP)	$(INITDIR)/spline.h
$(OBJDIR)/error.mpi.hdf5.o:	$(ANAL_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/error.mpi.o:		$(ANAL_DEP)
$(OBJDIR)/extract.mpi.hdf5.o:	$(ANAL_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/extract.mpi.o:		$(ANAL_DEP)
$(OBJDIR)/m31obs.mpi.hdf5.o:	$(ANAL_DEP)	$(M31_DEP)	$(MYINC)/rotate.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/m31obs.mpi.o:		$(ANAL_DEP)	$(M31_DEP)	$(MYINC)/rotate.h
$(OBJDIR)/m31coord.o:		$(COMMON_DEP)	$(MYINC)/constants.h	$(MYINC)/rotate.h	$(MISCDIR)/structure.h	$(M31_DEP)
#################################################################################################
#################################################################################################
