#################################################################################################
# last updated on 2018/01/17 (Wed) 20:16:54
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
# Macros for Pre-Processor
DEBUG	:= -DNDEBUG
# PROFILE	:= -pg
#################################################################################################
# Execution options
FORCE_SINGLE_GPU_RUN	:= 1
ENCLOSING_BALL_FOR_LET	:= 1
COMMUNICATION_VIA_HOST	:= 1
USE_MPI_PUT_FOR_LET	:= 0
USE_MPI_PUT_FOR_EXCG	:= 0
CARE_EXTERNAL_PARTICLES	:= 0
ACC_ACCUMULATION_IN_DP	:= 0
KAHAN_SUM_CORRECTION	:= 0
MONITOR_ENERGY_ERROR	:= 1
MONITOR_LETGEN_TIME	:= 1
MEASURE_BY_CUDA_EVENT	:= 0
DIVERT_GEOMETRIC_CENTER	:= 1
ADOPT_BLOCK_TIME_STEP	:= 1
ADOPT_VECTOR_ACC_MAC	:= 0
ADOPT_GADGET_TYPE_MAC	:= 1
ADOPT_WS93_TYPE_MAC	:= 1
IJ_PARALLELIZED_WALK	:= 1
DATAFILE_FORMAT_HDF5	:= 1
HDF5_FOR_ZINDAIJI	:= 0
DUMPFILE_IN_TIPSY	:= 0
DUMPFILE_AS_GALACTICS	:= 0
USE_OFFICIAL_SFMT	:= 1
USE_OFFICIAL_SFMT_JUMP	:= 1
SET_EXTERNAL_FIELD	:= 1
SET_EXTERNAL_FIELD_SPHE	:= 1
SET_EXTERNAL_FIELD_DISK	:= 0
#################################################################################################
# Debugging options
EVALUATE_FORCE_ERROR	:= 0
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
MYINC	:= $(MYDIR)/inc.tca
MYLIB	:= $(MYDIR)/lib.tca
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
# CCARG	+= -DTFIN_MIN_BENCHMARK="(60)"
CCARG	+= -DTFIN_MIN_BENCHMARK="(1430)"
DEBUG	:= -DNDEBUG
PROFILE	:=
USEDBG	:= 0
USEPAPI	:= 0
EVALUATE_FORCE_ERROR	:= 0
# REPORT_ELAPSED_TIME	:= 0
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
ifeq ($(MEASURE_BY_CUDA_EVENT), 1)
CCARG	+= -DUSE_CUDA_EVENT
CUARG	+= -DUSE_CUDA_EVENT
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
COMMUNICATION_VIA_HOST	:= 0
else
# direct solver is not parallelized
EVALUATE_FORCE_ERROR	:= 0
ifeq ($(MONITOR_LETGEN_TIME), 1)
CCARG	+= -DMONITOR_LETGEN_TIME
CUARG	+= -DMONITOR_LETGEN_TIME
endif
ifeq ($(COMMUNICATION_VIA_HOST), 1)
CCARG	+= -DMPI_VIA_HOST
CUARG	+= -DMPI_VIA_HOST
endif
ifeq ($(USE_MPI_PUT_FOR_LET), 1)
CCARG	+= -DMPI_ONE_SIDED_FOR_LET
CUARG	+= -DMPI_ONE_SIDED_FOR_LET
endif
ifeq ($(USE_MPI_PUT_FOR_EXCG), 1)
CCARG	+= -DMPI_ONE_SIDED_FOR_EXCG
CUARG	+= -DMPI_ONE_SIDED_FOR_EXCG
endif
ifeq ($(ENCLOSING_BALL_FOR_LET), 1)
CCARG	+= -DUSE_ENCLOSING_BALL_FOR_LET
CUARG	+= -DUSE_ENCLOSING_BALL_FOR_LET
DIVERT_GEOMETRIC_CENTER	:= 1
endif
ifeq ($(DIVERT_GEOMETRIC_CENTER), 1)
CCARG	+= -DRETURN_CENTER_BY_PHKEY_GENERATOR
CUARG	+= -DRETURN_CENTER_BY_PHKEY_GENERATOR
endif
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
ifeq ($(USEPAPI), 0)
PAPIOPT	:=
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
SET_EXTERNAL_FIELD_SPHE	:= 0
SET_EXTERNAL_FIELD_DISK	:= 0
endif
#################################################################################################
ifeq ($(SET_EXTERNAL_FIELD_SPHE), 1)
CCARG	+= -DSET_EXTERNAL_POTENTIAL_FIELD_SPHERICAL
CUARG	+= -DSET_EXTERNAL_POTENTIAL_FIELD_SPHERICAL
endif
#################################################################################################
ifeq ($(SET_EXTERNAL_FIELD_DISK), 1)
CCARG	+= -DSET_EXTERNAL_POTENTIAL_FIELD_DISK
CUARG	+= -DSET_EXTERNAL_POTENTIAL_FIELD_DISK
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
ifeq ($(HDF5_FOR_ZINDAIJI), 1)
CCARG	+= -DHDF5_FOR_ZINDAIJI
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
PLTENE	:= $(BINDIR)/plot.energy
PLTACT	:= $(BINDIR)/plot.action
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
ANALERR	:= $(BINDIR)/anal.error
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
TAGAERR	:= GOTHIC_err
TAGDUMP	:= GOTHIC
TAGDISK	:= MAGI_disk
TAGDIST	:= MAGI_df
TAGPROF	:= MAGI_profile
DIRBODY	:= $(PLGDIR)/$(TAGBODY)
DIRSNAP	:= $(PLGDIR)/$(TAGSNAP)
DIRAERR	:= $(PLGDIR)/$(TAGAERR)
DIRDUMP	:= $(PLGDIR)/$(TAGDUMP)
DIRDISK	:= $(PLGDIR)/$(TAGDISK)
DIRDIST	:= $(PLGDIR)/$(TAGDIST)
DIRPROF	:= $(PLGDIR)/$(TAGPROF)
XMLBODY	:= $(DIRBODY)/$(TAGBODY).xml
XMLSNAP	:= $(DIRSNAP)/$(TAGSNAP).xml
XMLAERR	:= $(DIRAERR)/$(TAGAERR).xml
XMLDUMP	:= $(DIRDUMP)/$(TAGDUMP).xml
XMLDISK	:= $(DIRDISK)/$(TAGDISK).xml
XMLDIST	:= $(DIRDIST)/$(TAGDIST).xml
XMLPROF	:= $(DIRPROF)/$(TAGPROF).xml
SRCBODY	:= $(DIRBODY)/$(TAGBODY)CommonPluginInfo.C $(DIRBODY)/$(TAGBODY)EnginePluginInfo.C $(DIRBODY)/$(TAGBODY)MDServerPluginInfo.C
SRCSNAP	:= $(DIRSNAP)/$(TAGSNAP)CommonPluginInfo.C $(DIRSNAP)/$(TAGSNAP)EnginePluginInfo.C $(DIRSNAP)/$(TAGSNAP)MDServerPluginInfo.C
SRCAERR	:= $(DIRAERR)/$(TAGAERR)CommonPluginInfo.C $(DIRAERR)/$(TAGAERR)EnginePluginInfo.C $(DIRAERR)/$(TAGAERR)MDServerPluginInfo.C
SRCDUMP	:= $(DIRDUMP)/$(TAGDUMP)CommonPluginInfo.C $(DIRDUMP)/$(TAGDUMP)EnginePluginInfo.C $(DIRDUMP)/$(TAGDUMP)MDServerPluginInfo.C
SRCDISK	:= $(DIRDISK)/$(TAGDISK)CommonPluginInfo.C $(DIRDISK)/$(TAGDISK)EnginePluginInfo.C $(DIRDISK)/$(TAGDISK)MDServerPluginInfo.C
SRCDIST	:= $(DIRDIST)/$(TAGDIST)CommonPluginInfo.C $(DIRDIST)/$(TAGDIST)EnginePluginInfo.C $(DIRDIST)/$(TAGDIST)MDServerPluginInfo.C
SRCPROF	:= $(DIRPROF)/$(TAGPROF)CommonPluginInfo.C $(DIRPROF)/$(TAGPROF)EnginePluginInfo.C $(DIRPROF)/$(TAGPROF)MDServerPluginInfo.C
SRCBODY	+= $(DIRBODY)/$(TAGBODY)PluginInfo.C $(DIRBODY)/$(TAGBODY)PluginInfo.h $(DIRBODY)/avt$(TAGBODY)FileFormat.C $(DIRBODY)/avt$(TAGBODY)FileFormat.h
SRCSNAP	+= $(DIRSNAP)/$(TAGSNAP)PluginInfo.C $(DIRSNAP)/$(TAGSNAP)PluginInfo.h $(DIRSNAP)/avt$(TAGSNAP)FileFormat.C $(DIRSNAP)/avt$(TAGSNAP)FileFormat.h
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
# TREESRC	+= make.c macutil.c
# TREESRC	:= make.c
TREEGPU	:= neighbor_dev.cu
UTILGPU	+= shrink_dev.cu
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
PENESRC	:= plot.energy.c
PACTSRC	:= plot.action.c
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
EXCGGPU	:= exchange_dev.cu
#################################################################################################
AERRSRC	:= anal.error.c
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
OBJPACT	:= $(patsubst %.c, $(OBJDIR)/%.mpi.gsl.pl.hdf5.o, $(notdir $(PACTSRC)))
OBJPACT	+= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o,        $(notdir $(FILELIB) $(ALLCLIB)))
else
OBJPACT	:= $(patsubst %.c, $(OBJDIR)/%.mpi.gsl.pl.o,      $(notdir $(PACTSRC)))
OBJPACT	+= $(patsubst %.c, $(OBJDIR)/%.o,		  $(notdir $(ALLCLIB)))
OBJPACT	+= $(patsubst %.c, $(OBJDIR)/%.mpi.o,		  $(notdir $(FILELIB)))
endif
OBJPACT	+= $(patsubst %.c, $(OBJDIR)/%.omp.o,		  $(notdir spline.c))
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
ifeq ($(DATAFILE_FORMAT_HDF5), 1)
OBJAERR	:= $(patsubst %.c, $(OBJDIR)/%.mpi.hdf5.o, $(notdir $(AERRSRC) $(FILELIB) $(ALLCLIB)))
else
OBJPCDF	:= $(patsubst %.c, $(OBJDIR)/%.mpi.o,      $(notdir $(AERRSRC) $(FILELIB)))
OBJPCDF	+= $(patsubst %.c, $(OBJDIR)/%.o,          $(notdir $(ALLCLIB)))
endif
#################################################################################################


#################################################################################################
## Rules
#################################################################################################
all:	TAGS $(GOTHIC) $(MAGI) $(EDITOR) $(PLTENE) $(PLTACT) $(PLTDST)# $(MKCOLD) $(PLTCDF)# $(PLTMUL)# $(OBJASM)
#################################################################################################
.PHONY:	gothic init magi cold editor plot bench sample disk anal
gothic:	TAGS $(GOTHIC)
init:	TAGS $(MAGI) $(MKCOLD) $(EDITOR)
magi:	TAGS $(MAGI)
cold:	TAGS $(MKCOLD)
editor:	TAGS $(EDITOR)
plot:	TAGS $(PLTENE) $(PLTACT) $(PLTDST) $(PLTCDF)# $(PLTMUL)
bench:	TAGS $(OPTCFG) $(PLTELP) $(PLTDEP) $(PLTBRK) $(PLTCMP) $(PLTFLP) $(PLTRAD)
sample:	TAGS $(SAMPLE) $(PLTDF)
disk:	TAGS $(PLTJET) $(PLTDISK)
anal:	TAGS $(ANALERR)
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
ifeq ($(USE_OFFICIAL_SFMT), 1)
ifeq ($(USE_OFFICIAL_SFMT_JUMP), 1)
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a	$(MYLIB)/libsfmt$(SFMTPER).a	$(MYLIB)/libsmtj$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)rand_sfmt$(SFMTPER) $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)hdf5lib $(HDF5LIB) $(SMTJLIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a	$(MYLIB)/libsfmt$(SFMTPER).a	$(MYLIB)/libsmtj$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)ompmpilib -l$(LIBPREC)rand_sfmt$(SFMTPER) $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)hdf5lib $(HDF5LIB) $(SMTJLIB)
else
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a	$(MYLIB)/libsfmt$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)rand_sfmt$(SFMTPER) $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)hdf5lib $(HDF5LIB) $(SFMTLIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a	$(MYLIB)/libsfmt$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)ompmpilib -l$(LIBPREC)rand_sfmt$(SFMTPER) $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)hdf5lib $(HDF5LIB) $(SFMTLIB)
endif
else
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_gsl.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)rand_gsl $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_gsl.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)ompmpilib -l$(LIBPREC)rand_gsl $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)hdf5lib $(HDF5LIB)
endif
$(EDITOR):	$(OBJEDIT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJEDIT) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)ompmpilib -l$(LIBPREC)timer -l$(LIBPREC)rotate -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(PLTENE):	$(OBJPENE)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPENE) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(PLTACT):	$(OBJPACT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPACT) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib $(GSLLIB) -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(PLTDST):	$(OBJDIST)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJDIST) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib $(GSLLIB) -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(PLTCDF):	$(OBJPCDF)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPCDF) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(PLTMUL):	$(OBJPMUL)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPMUL) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib $(GSLLIB) -l$(LIBPREC)hdf5lib $(HDF5LIB)
$(ANALERR):	$(OBJAERR)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)hdf5lib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJAERR) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)hdf5lib $(HDF5LIB)
else
ifeq ($(USE_OFFICIAL_SFMT), 1)
ifeq ($(USE_OFFICIAL_SFMT_JUMP), 1)
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/libsfmt$(SFMTPER).a	$(MYLIB)/libsmtj$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)rand_sfmt$(SFMTPER) $(GSLLIB) -l$(LIBPREC)timer $(SMTJLIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/libsfmt$(SFMTPER).a	$(MYLIB)/libsmtj$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)ompmpilib -l$(LIBPREC)rand_sfmt$(SFMTPER) $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)rotate $(SMTJLIB)
else
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/libsfmt$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)rand_sfmt$(SFMTPER) $(GSLLIB) -l$(LIBPREC)timer $(SFMTLIB)
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_sfmt$(SFMTPER).a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a	$(MYLIB)/libsfmt$(SFMTPER).a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)ompmpilib -l$(LIBPREC)rand_sfmt$(SFMTPER) $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)rotate $(SFMTLIB)
endif
else
$(MKCOLD):	$(OBJCOLD)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)rand_gsl.a	$(MYLIB)/lib$(LIBPREC)timer.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJCOLD) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)rand_gsl $(GSLLIB) -l$(LIBPREC)timer
$(MAGI):	$(OBJMAGI)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)rand_gsl.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJMAGI) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)ompmpilib -l$(LIBPREC)rand_gsl $(GSLLIB) -l$(LIBPREC)timer -l$(LIBPREC)rotate
endif
$(EDITOR):	$(OBJEDIT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)ompmpilib.a	$(MYLIB)/lib$(LIBPREC)timer.a	$(MYLIB)/lib$(LIBPREC)rotate.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJEDIT) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)ompmpilib -l$(LIBPREC)timer -l$(LIBPREC)rotate
$(PLTENE):	$(OBJPENE)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPENE) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib
$(PLTACT):	$(OBJPACT)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPACT) $(CCLIB) $(OMPLIB) -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib $(GSLLIB)
$(PLTDST):	$(OBJDIST)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJDIST) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib $(GSLLIB)
$(PLTCDF):	$(OBJPCDF)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPCDF) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib
$(PLTMUL):	$(OBJPMUL)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a	$(MYLIB)/lib$(LIBPREC)plplotlib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJPMUL) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib -l$(LIBPREC)plplotlib $(GSLLIB)
$(ANALERR):	$(OBJAERR)	$(MYLIB)/lib$(LIBPREC)myutil.a	$(MYLIB)/lib$(LIBPREC)constants.a	$(MYLIB)/lib$(LIBPREC)mpilib.a
	$(VERBOSE)$(MPICC) $(CCFLAG) $(CCDBG) $(PROFILE) -o $@ $(OBJAERR) $(CCLIB)           -L$(MYLIB) -l$(LIBPREC)myutil -l$(LIBPREC)constants -l$(LIBPREC)mpilib
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
	$(CC) $(DEBUG) $(CCFLAG) $(CCARG) -fsyntax-only $(CCINC) $(OMPINC) $(MPIINC) $(GSLINC) $(SFMTINC) $(SMTJINC) $(PLINC) -I/$(HOME)/tau/include $(HDF5INC) $(PREC) $(UNIT) $(CHK_SOURCES)
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
	README.txt LICENSE.txt Makefile $(SRCDIR) sh cfg host plt py plugins/README \
	$(XMLBODY) $(SRCBODY) $(XMLSNAP) $(SRCSNAP) $(XMLAERR) $(SRCAERR) $(XMLDUMP) $(SRCDUMP) $(XMLDISK) $(SRCDISK) $(XMLDIST) $(SRCDIST) $(XMLPROF) $(SRCPROF)
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
	$(VERBOSE)rm -f $(EDITOR) $(OBJEDIT)
	$(VERBOSE)rm -f $(PLTENE) $(OBJPENE)
	$(VERBOSE)rm -f $(PLTACT) $(OBJPACT)
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
	$(VERBOSE)rm -f $(ANALERR) $(OBJAERR)
	$(VERBOSE)rm -f $(OBJDIR)/*.o
#################################################################################################
visit:	$(DIRBODY)/Makefile $(DIRSNAP)/Makefile $(DIRAERR)/Makefile $(DIRDUMP)/Makefile $(DIRDISK)/Makefile $(DIRDIST)/Makefile $(DIRPROF)/Makefile
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
GOTHIC_DEP	+=	$(TREEDIR)/macutil.h	$(TREEDIR)/make.h	$(TREEDIR)/let.h	$(TREEDIR)/buf_inc.h	$(TREEDIR)/make_dev.h	$(TREEDIR)/walk_dev.h
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
$(OBJDIR)/peano_dev.o:	$(PEANO_DEP)	$(TREEDIR)/macutil.h	$(MYINC)/cudalib.h	$(UTILDIR)/gsync_dev.cu	$(SORTDIR)/peano_dev.h	$(PARADIR)/exchange_dev.h	$(UTILDIR)/compare_vec4_inc.cu	$(UTILDIR)/compare_vec4_inc.cuh	$(UTILDIR)/compare_vec4_del.cuh
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
TREE_LET_DEP	:=	$(MYINC)/mpilib.h	$(TREEDIR)/let.h	$(TREEDIR)/let_dev.h
$(OBJDIR)/let_dev.mpi.o:	$(TREE_DEV_DEP)	$(TREE_BUF_DEP)	$(TREE_LET_DEP)	$(MYINC)/timer.h	$(TREEDIR)/walk_dev.h
MAKE_DEV_DEP	:=	$(UTILDIR)/gsync_dev.cu	$(TREEDIR)/make_dev.h	$(TREEDIR)/let.h
MAKE_DEV_DEP	+=	$(UTILDIR)/scan_inc.cu	$(UTILDIR)/scan_inc.cuh	$(UTILDIR)/scan_del.cuh
MAKE_DEV_DEP	+=	$(UTILDIR)/scan_tsub_inc.cu	$(UTILDIR)/scan_tsub_inc.cuh	$(UTILDIR)/scan_tsub_del.cuh
MAKE_DEV_DEP	+=	$(UTILDIR)/compare_tsub_inc.cu	$(UTILDIR)/compare_tsub_inc.cuh	$(UTILDIR)/compare_tsub_del.cuh
$(OBJDIR)/make_dev.mpi.o:	$(TREE_DEV_DEP)	$(SORTDIR)/peano.h	$(SORTDIR)/peano_dev.h	$(MAKE_DEV_DEP)	$(MYINC)/timer.h	$(MISCDIR)/device.h
$(OBJDIR)/make_dev.o:		$(TREE_DEV_DEP)	$(SORTDIR)/peano.h	$(SORTDIR)/peano_dev.h	$(MAKE_DEV_DEP)
$(OBJDIR)/make.o:		$(TREE_DEP)	$(SORTDIR)/peano.h	$(TREEDIR)/macutil.h
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
$(OBJDIR)/external.omp.o:	$(SPLINE_DEP)	$(INITDIR)/external.h
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
$(OBJDIR)/plot.cdf.mpi.pl.hdf5.o:		$(PLOT_DEP)	$(MYINC)/name.h	$(PLOTDIR)/cdflib.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/plot.cdf.mpi.pl.o:			$(PLOT_DEP)	$(MYINC)/name.h	$(PLOTDIR)/cdflib.h
$(OBJDIR)/plot.action.mpi.gsl.pl.hdf5.o:	$(PLOT_DEP)	$(INITDIR)/spline.h	$(MYINC)/hdf5lib.h
$(OBJDIR)/plot.action.mpi.gsl.pl.o:		$(PLOT_DEP)	$(INITDIR)/spline.h
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
$(OBJDIR)/plot.df.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h	$(MYINC)/constants.h	$(INITDIR)/magi.h	$(INITDIR)/eddington.h
$(OBJDIR)/plot.performance.pl.o:	$(COMMON_DEP)	$(MYINC)/plplotlib.h
$(OBJDIR)/plot.needle.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h	$(MYINC)/constants.h	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(FILEDIR)/io.h
$(OBJDIR)/plot.disk.mpi.pl.hdf5.o:	$(COMMON_DEP)	$(MYINC)/myutil.h	$(MYINC)/plplotlib.h	$(MYINC)/name.h	$(MYINC)/hdf5lib.h	$(MYINC)/constants.h
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
# PARA_DEP	:=	$(COMMON_DEP)	$(MYINC)/name.h	$(MYINC)/mpilib.h
# PARA_DEP	+=	$(MISCDIR)/structure.h	$(PARADIR)/mpicfg.h
# EXCG_DEP	:=	$(PARA_DEP)	$(TIMEDIR)/adv_dev.h	$(PARADIR)/exchange.h	$(MISCDIR)/tune.h	$(MISCDIR)/brent.h
AERR_DEP	:=	$(COMMON_DEP)	$(MYINC)/constants.h	$(MYINC)/name.h	$(MYINC)/mpilib.h
AERR_DEP	+=	$(MISCDIR)/structure.h	$(MISCDIR)/allocate.h	$(FILEDIR)/io.h
$(OBJDIR)/anal.error.mpi.hdf5.o:	$(AERR_DEP)	$(MYINC)/hdf5lib.h
$(OBJDIR)/anal.error.mpi.o:		$(AERR_DEP)
#################################################################################################
#################################################################################################
