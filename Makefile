#################################################################################################
# last updated on 2016/09/06(Tue) 18:04:52
# Makefile for C Programming
# Calculation Code for OcTree Collisionless N-body Simulation on GPUs
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
USETCA	:= 1
#################################################################################################
# Macros for Pre-Processor
DEBUG	:= -DNDEBUG
# PROFILE	:= -pg
#################################################################################################
# Execution options
MONITOR_ENERGY_ERROR	:= 1
FORCE_SINGLE_GPU_RUN	:= 0
COMMUNICATION_VIA_HOST	:= 1
CONSTRUCT_LET_ON_GPU	:= 1
DIVERT_GEOMETRIC_CENTER	:= 1
CONSTRUCT_TREE_ON_GPU	:= 1
ADOPT_BLOCK_TIME_STEP	:= 1
ADOPT_GADGET_TYPE_MAC	:= 1
ADOPT_WS93_TYPE_MAC	:= 1
IJ_PARALLELIZED_WALK	:= 1
COMPARE_REBUILD_TIMING	:= 1
GENERATE_PHKEY_ON_GPU	:= 1
CALC_MULTIPOLE_ON_GPU	:= 1
LOCALIZE_I_PARTICLES	:= 1
USE_BRENT_METHOD	:= 1
BRUTE_FORCE_LOCALIZER	:= 1
FACILE_NEIGHBOR_SEARCH	:= 1
ADAPTIVE_PHYEY_JUMP	:= 0
DATAFILE_FORMAT_HDF5	:= 1
HDF5_FOR_ZINDAIJI	:= 1
APPEND_ASCII_ICDATA	:= 0
DUMPFILE_FOR_BONSAI	:= 0
#################################################################################################
# Debugging options
EVALUATE_FORCE_ERROR	:= 0
DEBUG_PARALLEL_HDF5	:= 0
DEBUG_LETGEN_ON_GPU	:= 0
DEBUG_TREE_TRAVERSAL	:= 0
DEBUG_MULTIPOLE_GPU	:= 0
ALLOCATE_LETBUFFER	:= 0
#################################################################################################
# Benchmark options
REPORT_ELAPSED_TIME	:= 1
MEASURE_ELAPSED_TIME	:= 0
MEASURE_EXEC_METRICS	:= 0
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
ifeq ($(findstring hapacs, $(HOSTNAME)), hapacs)
MYDIR	:= /work/GALAXY/$(USER)
ifeq ($(USETCA), 1)
MYINC	:= $(MYDIR)/inc.tca
MYLIB	:= $(MYDIR)/lib.tca
else
MYINC	:= $(MYDIR)/inc.comq
MYLIB	:= $(MYDIR)/lib.comq
endif
endif
#################################################################################################
include	$(MYINC)/common.mk
#################################################################################################
ifeq ($(USEHDF5), 0)
DATAFILE_FORMAT_HDF5	:= 0
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
CCARG	+= -DTFIN_MIN_BENCHMARK="(60)"
DEBUG	:= -DNDEBUG
PROFILE	:=
USEDBG	:= 0
USEPAPI	:= 0
EVALUATE_FORCE_ERROR	:= 0
DEBUG_PARALLEL_HDF5	:= 0
DEBUG_LETGEN_ON_GPU	:= 0
DEBUG_TREE_TRAVERSAL	:= 0
DEBUG_MULTIPOLE_GPU	:= 0
REPORT_ELAPSED_TIME	:= 0
endif
#################################################################################################
ifeq ($(EVALUATE_FORCE_ERROR), 1)
REPORT_ELAPSED_TIME	:= 0
endif
#################################################################################################
ifeq ($(REPORT_ELAPSED_TIME), 1)
CCARG	+= -DREPORT_TOTAL_ELAPSED_TIME
endif
#################################################################################################
ifeq ($(CHECK_FUNCTION_CALLS), 1)
DEBUG	:=
PROFILE	:=
USEDBG	:= 0
USEPAPI	:= 0
EVALUATE_FORCE_ERROR	:= 0
DEBUG_PARALLEL_HDF5	:= 0
DEBUG_LETGEN_ON_GPU	:= 0
DEBUG_TREE_TRAVERSAL	:= 0
DEBUG_MULTIPOLE_GPU	:= 0
REPORT_ELAPSED_TIME	:= 0
endif
#################################################################################################
ifeq ($(DEBUG_PARALLEL_HDF5), 1)
CCARG	+= -DDBG_PARALLEL_HDF5
CUARG	+= -DDBG_PARALLEL_HDF5
# DEBUG	:=
REPORT_ELAPSED_TIME	:= 0
endif
#################################################################################################
ifeq ($(DEBUG_LETGEN_ON_GPU), 1)
CCARG	+= -DDBG_LETGEN_ON_GPU
CUARG	+= -DDBG_LETGEN_ON_GPU
DEBUG	:=
REPORT_ELAPSED_TIME	:= 0
endif
#################################################################################################
ifeq ($(DEBUG_TREE_TRAVERSAL), 1)
CCARG	+= -DDBG_TREE_WALK
CUARG	+= -DDBG_TREE_WALK
DEBUG	:=
REPORT_ELAPSED_TIME	:= 0
endif
#################################################################################################
ifeq ($(DEBUG_MULTIPOLE_GPU), 1)
CCARG	+= -DDBG_CALC_MULTIPOLE
CUARG	+= -DDBG_CALC_MULTIPOLE
DEBUG	:=
REPORT_ELAPSED_TIME	:= 0
endif
#################################################################################################
ifeq ($(ALLOCATE_LETBUFFER), 1)
CCARG	+= -DALLOCATE_LETBUFFER
CUARG	+= -DALLOCATE_LETBUFFER
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
ifeq ($(FORCE_SINGLE_GPU_RUN), 1)
CCARG	+= -DSERIALIZED_EXECUTION
CUARG	+= -DSERIALIZED_EXECUTION
COMMUNICATION_VIA_HOST	:= 0
CONSTRUCT_LET_ON_GPU	:= 0
else
# direct solver is not parallelized
EVALUATE_FORCE_ERROR	:= 0
ifeq ($(CONSTRUCT_LET_ON_GPU), 1)
CCARG	+= -DBUILD_LET_ON_DEVICE
CUARG	+= -DBUILD_LET_ON_DEVICE
ifeq ($(COMMUNICATION_VIA_HOST), 1)
CCARG	+= -DLET_COMMUNICATION_VIA_HOST
CUARG	+= -DLET_COMMUNICATION_VIA_HOST
endif
ifeq ($(DIVERT_GEOMETRIC_CENTER), 1)
CCARG	+= -DRETURN_CENTER_BY_PHKEY_GENERATOR
CUARG	+= -DRETURN_CENTER_BY_PHKEY_GENERATOR
endif
endif
endif
#################################################################################################
ifeq ($(COMPARE_REBUILD_TIMING), 1)
CCARG	+= -DWALK_TREE_COMBINED_MODEL
CUARG	+= -DWALK_TREE_COMBINED_MODEL
endif
#################################################################################################
ifeq ($(CONSTRUCT_TREE_ON_GPU), 1)
CCARG	+= -DMAKE_TREE_ON_DEVICE
CUARG	+= -DMAKE_TREE_ON_DEVICE
endif
#################################################################################################
ifeq ($(ADOPT_BLOCK_TIME_STEP), 1)
CCARG	+= -DBLOCK_TIME_STEP
CUARG	+= -DBLOCK_TIME_STEP
else
ADAPTIVE_PHYEY_JUMP	:= 0
endif
#################################################################################################
ifeq ($(ADOPT_GADGET_TYPE_MAC), 1)
CCARG	+= -DGADGET_MAC
CUARG	+= -DGADGET_MAC
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
ifeq ($(GENERATE_PHKEY_ON_GPU), 1)
CCARG	+= -DGENERATE_PHKEY_ON_DEVICE
CUARG	+= -DGENERATE_PHKEY_ON_DEVICE
endif
#################################################################################################
ifeq ($(CALC_MULTIPOLE_ON_GPU), 1)
CCARG	+= -DCALC_MULTIPOLE_ON_DEVICE
CUARG	+= -DCALC_MULTIPOLE_ON_DEVICE
endif
#################################################################################################
ifeq ($(LOCALIZE_I_PARTICLES), 1)
CCARG	+= -DLOCALIZE_I_PARTICLES
CUARG	+= -DLOCALIZE_I_PARTICLES
else
USE_BRENT_METHOD	:= 0
BRUTE_FORCE_LOCALIZER	:= 0
FACILE_NEIGHBOR_SEARCH	:= 0
ADAPTIVE_PHYEY_JUMP	:= 0
endif
#################################################################################################
ifeq ($(BRUTE_FORCE_LOCALIZER), 1)
CCARG	+= -DBRUTE_FORCE_LOCALIZATION
CUARG	+= -DBRUTE_FORCE_LOCALIZATION
ifeq ($(FACILE_NEIGHBOR_SEARCH), 1)
CCARG	+= -DFACILE_NEIGHBOR_SEARCH
CUARG	+= -DFACILE_NEIGHBOR_SEARCH
ifeq ($(USE_BRENT_METHOD), 1)
CCARG	+= -DUSE_BRENT_METHOD
CUARG	+= -DUSE_BRENT_METHOD
endif
endif
else
ADAPTIVE_PHYEY_JUMP	:= 0
FACILE_NEIGHBOR_SEARCH	:= 0
endif
#################################################################################################
ifeq ($(ADAPTIVE_PHYEY_JUMP), 1)
CCARG	+= -DUSE_VARIABLE_NEIGHBOR_LEVEL
CUARG	+= -DUSE_VARIABLE_NEIGHBOR_LEVEL
endif
#################################################################################################
ifeq ($(IJ_PARALLELIZED_WALK), 1)
CCARG	+= -DIJ_PARALLELIZATION
CUARG	+= -DIJ_PARALLELIZATION
endif
#################################################################################################
ifeq ($(USEPAPI), 0)
PAPIOPT	:=
endif
#################################################################################################
NUM_NTHREADS	:= 512
NUM_TSUB	:= 32
NUM_NWARP	:= 4
NUM_NLOOP	:= 1
LEV_NEIGHBOR	:= 1
USE_WARPSHUFFLE	:= 1
PREF_SHARED_MEM	:= 1
LENGTH_RMAX	:= 0.8
LENGTH_FACTOR	:= 0.1
#################################################################################################
ifeq ($(MEASURE_ELAPSED_TIME), 1)
#################################################################################################
ifeq ($(HUNT_OPTIMAL_SEPARATION), 1)
LOCALIZE_I_PARTICLES	:= 1
BRUTE_FORCE_LOCALIZER	:= 1
FACILE_NEIGHBOR_SEARCH	:= 1
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
ifeq ($(LOCALIZE_I_PARTICLES), 1)
CCARG	+= -DNEIGHBOR_PHKEY_LEVEL="($(LEV_NEIGHBOR))"
CUARG	+= -DNEIGHBOR_PHKEY_LEVEL="($(LEV_NEIGHBOR))"
endif
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
CCARG	+= -DUSE_WARP_SHUFFLE_FUNC_MAKE_TREE -DUSE_WARP_SHUFFLE_FUNC_LINK_TREE
CUARG	+= -DUSE_WARP_SHUFFLE_FUNC_MAKE_TREE -DUSE_WARP_SHUFFLE_FUNC_LINK_TREE
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
else
HDF5_FOR_ZINDAIJI	:= 0
endif
#################################################################################################
ifeq ($(APPEND_ASCII_ICDATA), 1)
CCARG	+= -DAPPEND_ASCII_ICDATA
endif
#################################################################################################
ifeq ($(DUMPFILE_FOR_BONSAI), 1)
CCARG	+= -DDUMPFILE_FOR_BONSAI
endif
#################################################################################################
ifeq ($(HDF5_FOR_ZINDAIJI), 1)
CCARG	+= -DHDF5_FOR_ZINDAIJI
endif
#################################################################################################


#################################################################################################
## Source Directories
SRCDIR	:= src
MAINDIR	:= $(SRCDIR)/main
MISCDIR	:= $(SRCDIR)/misc
FILEDIR	:= $(SRCDIR)/file
TIMEDIR	:= $(SRCDIR)/time
SORTDIR	:= $(SRCDIR)/sort
TREEDIR	:= $(SRCDIR)/tree
INITDIR	:= $(SRCDIR)/init
PLOTDIR	:= $(SRCDIR)/plot
PARADIR	:= $(SRCDIR)/para
VPATH	:= $(MAINDIR) $(MISCDIR) $(FILEDIR) $(TIMEDIR) $(SORTDIR) $(TREEDIR) $(INITDIR) $(PLOTDIR) $(PARADIR)
#################################################################################################
## Executables for collisionless N-body code
GOTHIC	:= $(BINDIR)/gothic
MKCOLD	:= $(BINDIR)/uniformsphere
MAGI	:= $(BINDIR)/magi
PLTENE	:= $(BINDIR)/plot.energy
PLTDST	:= $(BINDIR)/plot.distribution
PLTCDF	:= $(BINDIR)/plot.cdf
PLTMUL	:= $(BINDIR)/plot.multipole
PLTELP	:= $(BINDIR)/plot.time
PLTDEP	:= $(BINDIR)/plot.ndep
PLTBRK	:= $(BINDIR)/plot.breakdown
PLTFLP	:= $(BINDIR)/plot.performance
PLTRAD	:= $(BINDIR)/plot.ball
PLTDF	:= $(BINDIR)/plot.df
PLTJET	:= $(BINDIR)/plot.needle
PLTDISK	:= $(BINDIR)/plot.disk
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
PLTCMP	:= $(BINDIR)/plot.comparison
endif
OPTCFG	:= $(BINDIR)/showOptConfig
SAMPLE	:= $(BINDIR)/sample
#################################################################################################
## Plugins for VisIt to read HDF5 files
PLGDIR	:= plugins
TAGBODY	:= GOTHIC_split
TAGSNAP	:= GOTHIC_snp
TAGDUMP	:= GOTHIC
TAGDISK	:= MAGI_disk
TAGDIST	:= MAGI_df
TAGPROF	:= MAGI_profile
DIRBODY	:= $(PLGDIR)/$(TAGBODY)
DIRSNAP	:= $(PLGDIR)/$(TAGSNAP)
DIRDUMP	:= $(PLGDIR)/$(TAGDUMP)
DIRDISK	:= $(PLGDIR)/$(TAGDISK)
DIRDIST	:= $(PLGDIR)/$(TAGDIST)
DIRPROF	:= $(PLGDIR)/$(TAGPROF)
XMLBODY	:= $(DIRBODY)/$(TAGBODY).xml
XMLSNAP	:= $(DIRSNAP)/$(TAGSNAP).xml
XMLDUMP	:= $(DIRDUMP)/$(TAGDUMP).xml
XMLDISK	:= $(DIRDISK)/$(TAGDISK).xml
XMLDIST	:= $(DIRDIST)/$(TAGDIST).xml
XMLPROF	:= $(DIRPROF)/$(TAGPROF).xml
SRCBODY	:= $(DIRBODY)/$(TAGBODY)CommonPluginInfo.C $(DIRBODY)/$(TAGBODY)EnginePluginInfo.C $(DIRBODY)/$(TAGBODY)MDServerPluginInfo.C
SRCSNAP	:= $(DIRSNAP)/$(TAGSNAP)CommonPluginInfo.C $(DIRSNAP)/$(TAGSNAP)EnginePluginInfo.C $(DIRSNAP)/$(TAGSNAP)MDServerPluginInfo.C
SRCDUMP	:= $(DIRDUMP)/$(TAGDUMP)CommonPluginInfo.C $(DIRDUMP)/$(TAGDUMP)EnginePluginInfo.C $(DIRDUMP)/$(TAGDUMP)MDServerPluginInfo.C
SRCDISK	:= $(DIRDISK)/$(TAGDISK)CommonPluginInfo.C $(DIRDISK)/$(TAGDISK)EnginePluginInfo.C $(DIRDISK)/$(TAGDISK)MDServerPluginInfo.C
SRCDIST	:= $(DIRDIST)/$(TAGDIST)CommonPluginInfo.C $(DIRDIST)/$(TAGDIST)EnginePluginInfo.C $(DIRDIST)/$(TAGDIST)MDServerPluginInfo.C
SRCPROF	:= $(DIRPROF)/$(TAGPROF)CommonPluginInfo.C $(DIRPROF)/$(TAGPROF)EnginePluginInfo.C $(DIRPROF)/$(TAGPROF)MDServerPluginInfo.C
SRCBODY	+= $(DIRBODY)/$(TAGBODY)PluginInfo.C $(DIRBODY)/$(TAGBODY)PluginInfo.h $(DIRBODY)/avt$(TAGBODY)FileFormat.C $(DIRBODY)/avt$(TAGBODY)FileFormat.h
SRCSNAP	+= $(DIRSNAP)/$(TAGSNAP)PluginInfo.C $(DIRSNAP)/$(TAGSNAP)PluginInfo.h $(DIRSNAP)/avt$(TAGSNAP)FileFormat.C $(DIRSNAP)/avt$(TAGSNAP)FileFormat.h
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
BRNTLIB	:= brent.c
UTILGPU	:= allocate_dev.cu
#################################################################################################
FILELIB	:= io.c
#################################################################################################
TREEGPU	:= adv_dev.cu
#################################################################################################
TREESRC	:= peano.c
ifeq ($(GENERATE_PHKEY_ON_GPU), 1)
UTILGPU	+= peano_dev.cu
endif
#################################################################################################
TREESRC	+= make.c macutil.c
TREEGPU	+= make_dev.cu walk_dev.cu neighbor_dev.cu
ifeq ($(BRUTE_FORCE_LOCALIZER), 1)
UTILGPU	+= shrink_dev.cu
else
UTILGPU	+= stat_dev.cu
endif
ifeq ($(FORCE_SINGLE_GPU_RUN), 0)
LETHOST	:= let.c
LET_GPU	:= let_dev.cu geo_dev.cu
endif
#################################################################################################
COLDSRC	:= uniformsphere.c
MAGISRC	:= magi.c
SMPLSRC	:= sample.c
# DISKLIB	:= disk.c
DISKLIB	:= potdens.c diskDF.c
MAGILIB	:= profile.c eddington.c king.c abel.c blas.c spline.c table.c
#################################################################################################
PENESRC	:= plot.energy.c
DISTSRC	:= plot.distribution.c
PCDFSRC	:= plot.cdf.c
PMULSRC	:= plot.multipole.c
PCDFLIB	:= cdflib.c
PELPSRC	:= plot.time.c
PDEPSRC	:= plot.ndep.c
PBRKSRC	:= plot.breakdown.c
PCMPSRC	:= plot.comparison.c
PFLPSRC	:= plot.performance.c
PRADSRC	:= plot.ball.c
PDFSRC	:= plot.df.c
PJETSRC	:= plot.needle.c
PDSKSRC	:= plot.disk.c
#################################################################################################
LETHOST	+= exchange.c mpicfg.c
#################################################################################################


#################################################################################################
## Object files
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJMPGT	:= $(patsubst %.c,  $(OBJDIR)/%.mpi.hdf5.o, $(notdir $(NBDYSRC) $(FILELIB) $(ALLCLIB) $(CNVTLIB)))
else
OBJMPGT	:= $(patsubst %.c,  $(OBJDIR)/%.mpi.o,      $(notdir $(NBDYSRC) $(FILELIB)))
OBJMPGT	+= $(patsubst %.c,  $(OBJDIR)/%.o,          $(notdir $(ALLCLIB) $(CNVTLIB)))
endif
OBJMPGT	+= $(patsubst %.c,  $(OBJDIR)/%.o,          $(notdir $(TREESRC)))
OBJMPGT	+= $(patsubst %.cu, $(OBJDIR)/%.o,          $(notdir $(UTILGPU)))
ifeq ($(FORCE_SINGLE_GPU_RUN), 1)
OBJMPGT	+= $(patsubst %.cu, $(OBJDIR)/%.o,          $(notdir $(TREEGPU)))
else
OBJMPGT	+= $(patsubst %.cu, $(OBJDIR)/%.mpi.o,      $(notdir $(TREEGPU) $(LET_GPU)))
OBJMPGT	+= $(patsubst %.c,  $(OBJDIR)/%.mpi.o,      $(notdir $(LETHOST)))
endif
ifeq ($(USE_BRENT_METHOD), 1)
OBJMPGT	+= $(patsubst %.c,  $(OBJDIR)/%.o,          $(notdir $(BRNTLIB)))
endif
#################################################################################################
OBJFIND	:= $(patsubst %.c, $(OBJDIR)/%.o, $(notdir $(FINDSRC)))
OBJSMPL	:= $(patsubst %.c, $(OBJDIR)/%.o, $(notdir $(SMPLSRC)))
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJCOLD	:= $(patsubst %.c, $(OBJDIR)/%.mpi.gsl.hdf5.o,      $(notdir $(COLDSRC)))
OBJCOLD	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,          $(notdir $(FILELIB) $(ALLCLIB)))
else
OBJCOLD	:= $(patsubst %.c, $(OBJDIR)/%.gsl.o,               $(notdir $(COLDSRC)))
OBJCOLD	+= $(patsubst %.c, $(OBJDIR)/%.o,                   $(notdir $(ALLCLIB)))
OBJCOLD	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,               $(notdir $(FILELIB)))
endif
#################################################################################################
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJMAGI	:= $(patsubst %.c, $(OBJDIR)/%.ompmpi.gsl.hdf5.o, $(notdir $(MAGISRC)))
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,        $(notdir $(FILELIB) $(ALLCLIB)))
else
OBJMAGI	:= $(patsubst %.c, $(OBJDIR)/%.ompmpi.gsl.o,      $(notdir $(MAGISRC)))
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.o,                 $(notdir $(ALLCLIB)))
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,             $(notdir $(FILELIB)))
endif
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.omp.gsl.o,         $(notdir $(DISKLIB)))
OBJMAGI	+= $(patsubst %.c, $(OBJDIR)/%.omp.o,             $(notdir $(MAGILIB)))
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
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJPMUL	:= $(patsubst %.c, $(OBJDIR)/%.mpi.gsl.pl.hdf5.o, $(notdir $(PMULSRC)))
OBJPMUL	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,        $(notdir $(FILELIB) $(ALLCLIB)))
else
OBJPMUL	:= $(patsubst %.c, $(OBJDIR)/%.mpi.gsl.pl.o,      $(notdir $(PMULSRC)))
OBJPMUL	+= $(patsubst %.c, $(OBJDIR)/%.o,                 $(notdir $(ALLCLIB)))
OBJPMUL	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,             $(notdir $(FILELIB)))
endif
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


#################################################################################################
## Rules
#################################################################################################
all:	TAGS $(GOTHIC) $(MAGI) $(MKCOLD) $(PLTENE) $(PLTDST) $(PLTCDF)# $(PLTMUL)# $(OBJASM)
#################################################################################################
.PHONY:	gothic init magi cold edf plot bench sample disk
gothic:	TAGS $(GOTHIC)
init:	TAGS $(MAGI) $(MKCOLD)
magi:	TAGS $(MAGI)
cold:	TAGS $(MKCOLD)
edf:	TAGS $(MAGI)
plot:	TAGS $(PLTENE) $(PLTDST) $(PLTCDF)# $(PLTMUL)
bench:	TAGS $(OPTCFG) $(PLTELP) $(PLTDEP) $(PLTBRK) $(PLTCMP) $(PLTFLP) $(PLTRAD)
sample:	TAGS $(SAMPLE) $(PLTDF)
disk:	TAGS $(PLTJET) $(PLTDISK)
sass:	TAGS $(GOTHIC).sass
#################################################################################################
## Making TAGS file for Emacs
TAGS:
	$(VERBOSE)find $(SRCDIR) -name "*.[ch]" -print | etags -
#################################################################################################
## CUDA/C code
# multiple GPUs code
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(GOTHIC):	$(OBJMPGT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)cudalib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMPGT) $(CCLIB) $(CULIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)timer -l$(LIBPREC)mpilib -l$(LIBPREC)cudalib -l$(LIBPREC)hdf5lib $(HDF5LIB)
else
$(GOTHIC):	$(OBJMPGT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)cudalib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMPGT) $(CCLIB) $(CULIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)timer -l$(LIBPREC)mpilib -l$(LIBPREC)cudalib
endif
#################################################################################################
## C code
$(OPTCFG):	$(OBJFIND)	$(MYLIB)/lib$(LIBPREC)myutil.a
	$(VERBOSE)$(CC)    $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJFIND) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil
$(SAMPLE):	$(OBJSMPL)
	$(VERBOSE)$(CC)    $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJSMPL) $(CCLIB)
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)ompmpilib $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(PLTENE):	$(OBJPENE)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPENE) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(PLTDST):	$(OBJDIST)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJDIST) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib $(GSLLIB) -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(PLTCDF):	$(OBJPCDF)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPCDF) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(PLTMUL):	$(OBJPMUL)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPMUL) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib $(GSLLIB) -l$(LIBPREC)hdf5lib $(HDF5LIB)
else
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)timer.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib $(GSLLIB) -l$(LIBPREC)timer
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)ompmpilib $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)rotate
$(PLTENE):	$(OBJPENE)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPENE) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib
$(PLTDST):	$(OBJDIST)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJDIST) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib $(GSLLIB)
$(PLTCDF):	$(OBJPCDF)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPCDF) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib
$(PLTMUL):	$(OBJPMUL)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPMUL) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib $(GSLLIB)
endif
$(PLTELP):	$(OBJPELP)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)mpilib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPELP) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib
$(PLTDEP):	$(OBJPDEP)	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPDEP) $(CCLIB) -L$(MYLIB) -l$(LIBPREC)plplotlib
$(PLTBRK):	$(OBJPBRK)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPBRK) $(CCLIB) -L$(MYLIB) -l$(LIBPREC)plplotlib -l$(LIBPREC)myutil
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(PLTCMP):	$(OBJPCMP)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPCMP) $(CCLIB) -L$(MYLIB) -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib $(HDF5LIB)
endif
$(PLTFLP):	$(OBJPFLP)	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPFLP) $(CCLIB) -L$(MYLIB) -l$(LIBPREC)plplotlib
$(PLTRAD):	$(OBJPRAD)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPRAD) $(CCLIB) -L$(MYLIB) -l$(LIBPREC)plplotlib -l$(LIBPREC)myutil -l$(LIBPREC)constants
$(PLTDF):	$(OBJPDF)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a	$(MYLIB)/lib$(LIBPREC)constants.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPDF) $(CCLIB) -L$(MYLIB) -l$(LIBPREC)constants -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib $(HDF5LIB)
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(PLTJET):	$(OBJPJET)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPJET) $(CCLIB) -L$(MYLIB) -l$(LIBPREC)plplotlib -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)hdf5lib $(HDF5LIB)
else
$(PLTJET):	$(OBJPJET)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPJET) $(CCLIB) -L$(MYLIB) -l$(LIBPREC)plplotlib -l$(LIBPREC)myutil -l$(LIBPREC)constants
endif
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
$(PLTDISK):	$(OBJPDSK)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPDSK) $(CCLIB) -L$(MYLIB) -l$(LIBPREC)plplotlib -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)hdf5lib $(HDF5LIB)
else
$(PLTDISK):	$(OBJPDSK)	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a
	$(VERBOSE)$(CC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPDSK) $(CCLIB) -L$(MYLIB) -l$(LIBPREC)plplotlib -l$(LIBPREC)myutil -l$(LIBPREC)constants
endif
#################################################################################################
# sass file
$(GOTHIC).sass:	$(GOTHIC)
	$(VERBOSE)$(CUASM) $< > $@
#################################################################################################
# original version
check-syntax:
	$(CC) $(DEBUG) $(CCFLAG) $(CCARG) -fsyntax-only $(CCINC) $(OMPINC) $(MPIINC) $(GSLINC) $(PLINC) -I/$(HOME)/tau/include $(HDF5INC) $(PREC) $(UNIT) $(CHK_SOURCES)
check-syntax-cu:
	$(CU) $(DEBUG) $(CUARG) -Xcompiler -fsyntax-only $(CCINC) $(CUINC) $(OMPINC) $(MPIINC) $(GSLINC) $(PLINC) -I/$(HOME)/tau/include $(HDF5INC) $(PREC) $(UNIT) -o /dev/null -c $(CHK_SOURCES)
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
	README Makefile $(SRCDIR) sh cfg host plt plugins/README \
	$(XMLBODY) $(SRCBODY) $(XMLSNAP) $(SRCSNAP) $(XMLDUMP) $(SRCDUMP) $(XMLDISK) $(SRCDISK) $(XMLDIST) $(SRCDIST) $(XMLPROF) $(SRCPROF)
	$(VERBOSE)mv       $(DATE)tree.tar.xz pub/
#################################################################################################
clean:
	$(VERBOSE)rm -f TAGS gmon.out
	$(VERBOSE)rm -f $(ASMDIR)/*.s
	$(VERBOSE)rm -f $(BINDIR)/*.log
	$(VERBOSE)rm -f $(OBJDIR)/*.log
	$(VERBOSE)rm -f $(GOTHIC) $(OBJMPGT) $(GOTHIC).sass
	$(VERBOSE)rm -f $(MKCOLD) $(OBJCOLD)
	$(VERBOSE)rm -f $(MAGI)   $(OBJMAGI)
	$(VERBOSE)rm -f $(PLTENE) $(OBJPENE)
	$(VERBOSE)rm -f $(PLTDST) $(OBJDIST)
	$(VERBOSE)rm -f $(PLTCDF) $(OBJPCDF)
	$(VERBOSE)rm -f $(PLTMUL) $(OBJPMUL)
	$(VERBOSE)rm -f $(PLTELP) $(OBJPELP)
	$(VERBOSE)rm -f $(PLTDEP) $(OBJPDEP)
	$(VERBOSE)rm -f $(PLTBRK) $(OBJPBRK)
	$(VERBOSE)rm -f $(PLTCMP) $(OBJPCMP)
	$(VERBOSE)rm -f $(PLTFLP) $(OBJPFLP)
	$(VERBOSE)rm -f $(PLTRAD) $(OBJPRAD)
	$(VERBOSE)rm -f $(PLTJET) $(OBJPJET)
	$(VERBOSE)rm -f $(OPTCFG) $(OBJFIND)
	$(VERBOSE)rm -f $(SAMPLE) $(OBJSMPL)
	$(VERBOSE)rm -f $(PLTDF)  $(OBJPDF)
	$(VERBOSE)rm -f $(PLTDISK) $(OBJPDSK)
	$(VERBOSE)rm -f $(OBJDIR)/*.o
#################################################################################################
visit:	$(DIRBODY)/Makefile $(DIRSNAP)/Makefile $(DIRDUMP)/Makefile $(DIRDISK)/Makefile $(DIRDIST)/Makefile $(DIRPROF)/Makefile
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
## Dependence
#################################################################################################
COMMON_DEP	:=	Makefile	$(MYINC)/common.mk	$(MYINC)/macro.h
#################################################################################################
## $(MAINDIR)/*
GOTHIC_DEP	:=	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/name.h	$(MYINC)/constants.h	$(MYINC)/timer.h	$(MYINC)/cudalib.h
GOTHIC_DEP	+=	$(MISCDIR)/device.h	$(MISCDIR)/benchmark.h	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(MISCDIR)/allocate_dev.h	$(MISCDIR)/convert.h
GOTHIC_DEP	+=	$(FILEDIR)/io.h
GOTHIC_DEP	+=	$(SORTDIR)/peano.h
GOTHIC_DEP	+=	$(TREEDIR)/macutil.h	$(TREEDIR)/make.h	$(TREEDIR)/let.h	$(TREEDIR)/buf_inc.h	$(TREEDIR)/make_dev.h	$(TREEDIR)/walk_dev.h
GOTHIC_DEP	+=	$(TIMEDIR)/adv_dev.h
ifeq ($(FORCE_SINGLE_GPU_RUN), 0)
GOTHIC_DEP	+=	$(MYINC)/mpilib.h	$(PARADIR)/mpicfg.h	$(PARADIR)/exchange.h
GOTHIC_DEP	+=	$(TREEDIR)/geo_dev.h	$(TREEDIR)/let_dev.h
endif
ifeq ($(USE_BRENT_METHOD), 1)
GOTHIC_DEP	+=	$(MISCDIR)/brent.h
endif
ifeq ($(GENERATE_PHKEY_ON_GPU), 1)
GOTHIC_DEP	+=	$(SORTDIR)/peano_dev.h
endif
ifeq ($(BRUTE_FORCE_LOCALIZER), 1)
ifeq ($(LOCALIZE_I_PARTICLES), 1)
GOTHIC_DEP	+=	$(TREEDIR)/neighbor_dev.h
endif
GOTHIC_DEP	+=	$(TREEDIR)/shrink_dev.h
else
GOTHIC_DEP	+=	$(TREEDIR)/stat_dev.h
endif
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
#################################################################################################
## $(FILEDIR)/*
IO_DEP	:=	$(COMMON_DEP)	$(MYINC)/constants.h	$(MYINC)/name.h	$(MYINC)/mpilib.h
IO_DEP	+=	$(MISCDIR)/benchmark.h	$(MISCDIR)/structure.h	$(FILEDIR)/io.h
ifeq ($(HUNT_OPTIMAL_WALK_TREE), 1)
IO_DEP	+=	$(TREEDIR)/walk_dev.h
ifeq ($(BRUTE_FORCE_LOCALIZER), 1)
ifeq ($(LOCALIZE_I_PARTICLES), 1)
IO_DEP	+=	$(TREEDIR)/neighbor_dev.h
endif
IO_DEP	+=	$(TREEDIR)/shrink_dev.h
else
IO_DEP	+=	$(TREEDIR)/stat_dev.h
endif
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
ifeq ($(BRUTE_FORCE_LOCALIZER), 1)
ifeq ($(LOCALIZE_I_PARTICLES), 1)
IO_DEP	+=	$(TREEDIR)/neighbor_dev.h
endif
endif
endif
$(OBJDIR)/io.mpi.hdf5.o:	$(IO_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/io.mpi.o:		$(IO_DEP)
#################################################################################################
## $(SORTDIR)/*
PEANO_DEP	:=	$(COMMON_DEP)	$(MISCDIR)/benchmark.h	$(MISCDIR)/structure.h	$(SORTDIR)/peano.h	$(TREEDIR)/make.h
$(OBJDIR)/peano.o:	$(PEANO_DEP)
$(OBJDIR)/peano_dev.o:	$(PEANO_DEP)	$(TREEDIR)/macutil.h	$(MYINC)/cudalib.h	$(MISCDIR)/gsync_dev.cu	$(SORTDIR)/peano_dev.h
#################################################################################################
## $(TIMEDIR)/*
ADV_DEV_DEP	:=	$(COMMON_DEP)	$(MYINC)/timer.h	$(MYINC)/cudalib.h
ADV_DEV_DEP	+=	$(MISCDIR)/structure.h	$(MISCDIR)/benchmark.h	$(MISCDIR)/device.h	$(TIMEDIR)/adv_dev.h	$(TREEDIR)/walk_dev.h
$(OBJDIR)/adv_dev.mpi.o:	$(ADV_DEV_DEP)	$(MYINC)/mpilib.h	$(PARADIR)/mpicfg.h
$(OBJDIR)/adv_dev.o:		$(ADV_DEV_DEP)
#################################################################################################
## $(TREEDIR)/*
$(OBJDIR)/macutil.o:	$(COMMON_DEP)	$(SORTDIR)/peano.h	$(TREEDIR)/macutil.h	$(TREEDIR)/make.h
TREE_DEP	:=	$(COMMON_DEP)	$(MISCDIR)/benchmark.h	$(MISCDIR)/structure.h	$(TREEDIR)/make.h
$(OBJDIR)/geo_dev.o:	$(TREE_DEP)	$(MYINC)/cudalib.h	$(TREEDIR)/walk_dev.h	$(TREEDIR)/geo_dev.h
$(OBJDIR)/let.mpi.o:	$(TREE_DEP)	$(MYINC)/mpilib.h	$(TREEDIR)/let.h	$(PARADIR)/mpicfg.h
TREE_DEV_DEP	:=	$(TREE_DEP)	$(MYINC)/cudalib.h	$(MISCDIR)/device.h	$(TREEDIR)/macutil.h
TREE_BUF_DEP	:=	$(TREEDIR)/buf_inc.h	$(TREEDIR)/buf_inc.cu
TREE_LET_DEP	:=	$(MYINC)/mpilib.h	$(TREEDIR)/let.h	$(TREEDIR)/let_dev.h	$(PARADIR)/mpicfg.h
$(OBJDIR)/let_dev.mpi.o:	$(TREE_DEV_DEP)	$(TREE_BUF_DEP)	$(TREE_LET_DEP)	$(MYINC)/timer.h	$(TREEDIR)/walk_dev.h
MAKE_DEV_DEP	:=	$(MISCDIR)/gsync_dev.cu	$(TREEDIR)/make_dev.h	$(TREEDIR)/let.h	$(TREEDIR)/make_inc.cu	$(TREEDIR)/make_inc.h	$(TREEDIR)/make_del.h
$(OBJDIR)/make_dev.mpi.o:	$(TREE_DEV_DEP)	$(SORTDIR)/peano.h	$(SORTDIR)/peano_dev.h	$(MAKE_DEV_DEP)
$(OBJDIR)/make_dev.o:		$(TREE_DEV_DEP)	$(SORTDIR)/peano.h	$(SORTDIR)/peano_dev.h	$(MAKE_DEV_DEP)
$(OBJDIR)/make.o:		$(TREE_DEP)	$(SORTDIR)/peano.h	$(TREEDIR)/macutil.h
TREE_RSORT_DEP	:=	$(MISCDIR)/gsync_dev.cu	$(SORTDIR)/radix_dev.h	$(SORTDIR)/radix_inc.cu	$(SORTDIR)/radix_inc.h	$(SORTDIR)/radix_del.h
NEIGHBOR_DEV_DEP	:=	$(TREE_RSORT_DEP)	$(MYINC)/cudalib.h	$(TREEDIR)/make_dev.h	$(TREEDIR)/neighbor_dev.h
$(OBJDIR)/neighbor_dev.mpi.o:	$(TREE_DEP)	$(NEIGHBOR_DEV_DEP)
$(OBJDIR)/neighbor_dev.o:	$(TREE_DEP)	$(NEIGHBOR_DEV_DEP)
SHRINK_DEV_DEP	:=	$(MYINC)/cudalib.h	$(MISCDIR)/device.h	$(TREEDIR)/walk_dev.h	$(TREEDIR)/make_dev.h	$(TREEDIR)/shrink_dev.h
ifeq ($(LOCALIZE_I_PARTICLES), 1)
ifeq ($(USE_BRENT_METHOD), 1)
SHRINK_DEV_DEP	+=	$(MISCDIR)/brent.h
endif
SHRINK_DEV_DEP	+=	$(TREEDIR)/neighbor_dev.h
endif
$(OBJDIR)/shrink_dev.mpi.o:	$(TREE_DEP)	$(SHRINK_DEV_DEP)
$(OBJDIR)/shrink_dev.o:		$(TREE_DEP)	$(SHRINK_DEV_DEP)
ifeq ($(LOCALIZE_I_PARTICLES), 1)
$(OBJDIR)/stat_dev.o:	$(TREE_DEP)	$(MYINC)/cudalib.h	$(MISCDIR)/device.h	$(TREEDIR)/walk_dev.h	$(TREEDIR)/stat_dev.h	$(SORTDIR)/peano.h
else
$(OBJDIR)/stat_dev.o:	$(TREE_DEP)	$(MYINC)/cudalib.h	$(MISCDIR)/device.h	$(TREEDIR)/walk_dev.h	$(TREEDIR)/stat_dev.h
endif
WALK_DEV_DEP	:=	$(TREE_DEV_DEP)	$(TREE_BUF_DEP)	$(MYINC)/timer.h	$(TREEDIR)/walk_dev.h	$(TREEDIR)/seb_dev.cu
$(OBJDIR)/walk_dev.mpi.o:	$(WALK_DEV_DEP)	$(TREE_LET_DEP)
$(OBJDIR)/walk_dev.o:		$(WALK_DEV_DEP)
#################################################################################################
## $(INITDIR)/*
$(OBJDIR)/sample.o:	$(COMMON_DEP)
$(OBJDIR)/blas.omp.o:	$(COMMON_DEP)	$(INITDIR)/blas.h
$(OBJDIR)/spline.omp.o:	$(COMMON_DEP)	$(INITDIR)/spline.h
SPLINE_DEP	:=	$(COMMON_DEP)	$(INITDIR)/profile.h	$(INITDIR)/table.h	$(INITDIR)/spline.h
$(OBJDIR)/abel.omp.o:	$(SPLINE_DEP)	$(INITDIR)/magi.h	$(INITDIR)/abel.h
$(OBJDIR)/table.omp.o:	$(SPLINE_DEP)
PROFILE_DEP	:=	$(COMMON_DEP)	$(MYINC)/constants.h	$(INITDIR)/profile.h
$(OBJDIR)/eddington.omp.o:	$(PROFILE_DEP)	$(INITDIR)/magi.h	$(INITDIR)/eddington.h
$(OBJDIR)/king.omp.o:		$(PROFILE_DEP)	$(INITDIR)/magi.h	$(INITDIR)/king.h
$(OBJDIR)/profile.omp.o:	$(PROFILE_DEP)
DISK_DEP	:=	$(PROFILE_DEP)	$(INITDIR)/magi.h	$(INITDIR)/spline.h	$(INITDIR)/blas.h	$(INITDIR)/potdens.h
$(OBJDIR)/potdens.omp.gsl.o:	$(DISK_DEP)
$(OBJDIR)/diskDF.omp.gsl.o:	$(DISK_DEP)	$(MISCDIR)/structure.h	$(INITDIR)/diskDF.h
$(OBJDIR)/disk.omp.gsl.o:	$(DISK_DEP)	$(MISCDIR)/structure.h	$(INITDIR)/disk.h
IOFILE_DEP	:=	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/constants.h	$(MYINC)/timer.h	$(MYINC)/name.h
IOFILE_DEP	+=	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(FILEDIR)/io.h
$(OBJDIR)/uniformsphere.mpi.gsl.hdf5.o:	$(IOFILE_DEP)			$(MYINC)/hdf5lib.h
$(OBJDIR)/uniformsphere.gsl.o:		$(IOFILE_DEP)
MAGI_DEP	:=	$(MYINC)/rotate.h	$(INITDIR)/magi.h	$(INITDIR)/king.h	$(INITDIR)/profile.h	$(INITDIR)/eddington.h
MAGI_DEP	+=	$(INITDIR)/table.h	$(INITDIR)/abel.h	$(INITDIR)/potdens.h	$(INITDIR)/diskDF.h
$(OBJDIR)/magi.ompmpi.gsl.hdf5.o:	$(IOFILE_DEP)	$(MAGI_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/magi.omp.gsl.o:		$(IOFILE_DEP)	$(MAGI_DEP)
#################################################################################################
## $(PLOTDIR)/*
PLOT_DEP	:=	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/constants.h	\
			$(MYINC)/plplotlib.h	$(MYINC)/mpilib.h	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(FILEDIR)/io.h
$(OBJDIR)/plot.cdf.mpi.pl.hdf5.o:		$(PLOT_DEP)	$(MYINC)/name.h	$(PLOTDIR)/cdflib.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/plot.cdf.mpi.pl.o:			$(PLOT_DEP)	$(MYINC)/name.h	$(PLOTDIR)/cdflib.h
$(OBJDIR)/plot.distribution.mpi.gsl.pl.hdf5.o:	$(PLOT_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/plot.distribution.mpi.gsl.pl.o:	$(PLOT_DEP)
$(OBJDIR)/plot.energy.mpi.pl.hdf5.o:		$(PLOT_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/plot.energy.mpi.pl.o:			$(PLOT_DEP)
$(OBJDIR)/plot.multipole.mpi.gsl.pl.hdf5.o:	$(PLOT_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/plot.multipole.mpi.gsl.pl.o:		$(PLOT_DEP)
$(OBJDIR)/plot.ndep.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h
$(OBJDIR)/plot.ball.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/myutil.h	$(MYINC)/name.h	$(MYINC)/constants.h
$(OBJDIR)/plot.breakdown.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/myutil.h	$(MYINC)/name.h	$(MISCDIR)/benchmark.h
$(OBJDIR)/plot.time.mpi.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/myutil.h	$(MYINC)/name.h	$(MYINC)/mpilib.h	$(PLOTDIR)/cdflib.h
$(OBJDIR)/plot.comparison.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/plot.df.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h	$(MYINC)/constants.h
$(OBJDIR)/plot.performance.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h
$(OBJDIR)/plot.needle.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h	$(MYINC)/constants.h	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(FILEDIR)/io.h
$(OBJDIR)/plot.disk.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h	$(MYINC)/constants.h
$(OBJDIR)/cdflib.o:	$(COMMON_DEP)	$(PLOTDIR)/cdflib.h
#################################################################################################
## $(PARADIR)/*
PARA_DEP	:=	$(COMMON_DEP)	$(MYINC)/name.h	$(MYINC)/mpilib.h
PARA_DEP	+=	$(MISCDIR)/structure.h	$(TREEDIR)/make.h	$(PARADIR)/mpicfg.h
$(OBJDIR)/exchange.mpi.o:	$(PARA_DEP)	$(TIMEDIR)/adv_dev.h	$(PARADIR)/exchange.h
$(OBJDIR)/mpicfg.mpi.o:		$(PARA_DEP)
#################################################################################################
#################################################################################################
