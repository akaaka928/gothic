/************************************************************************* \
 *                                                                       *
                  last updated on 2017/01/18(Wed) 11:05:29
 *                                                                       *
 *    N-body code based on Barnes--Hut tree                              *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>/* for on-the-fly monitoring of execution time of multiple functions */
#include <unistd.h>/* to check the existence of files */
#include <mpi.h>
//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
#include "macro.h"
#include "myutil.h"
#include "name.h"
#include "constants.h"
#include "timer.h"
#include "cudalib.h"
//-------------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
#include "mpilib.h"
#include "../para/mpicfg.h"
#endif//SERIALIZED_EXECUTION
//-------------------------------------------------------------------------
#include "../misc/device.h"
#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../misc/allocate_dev.h"
#include "../misc/convert.h"
#include "../misc/tune.h"
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#include "../misc/brent.h"
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
//-------------------------------------------------------------------------
#include "../file/io.h"
//-------------------------------------------------------------------------
#include "../sort/peano.h"
#ifdef  GENERATE_PHKEY_ON_DEVICE
#include "../sort/peano_dev.h"
#endif//GENERATE_PHKEY_ON_DEVICE
//-------------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
#include "../para/exchange.h"
#ifdef  EXCHANGE_USING_GPUS
#include "../para/exchange_dev.h"
#endif//EXCHANGE_USING_GPUS
#endif//SERIALIZED_EXECUTION
//-------------------------------------------------------------------------
#include "../tree/macutil.h"
#include "../tree/make.h"
#include "../tree/let.h"
#include "../tree/buf_inc.h"
#include "../tree/make_dev.h"
#include "../tree/walk_dev.h"
#ifdef  BRUTE_FORCE_LOCALIZATION
#ifdef  LOCALIZE_I_PARTICLES
#include "../tree/neighbor_dev.h"
#endif//LOCALIZE_I_PARTICLES
#include "../tree/shrink_dev.h"
#else///BRUTE_FORCE_LOCALIZATION
#include "../tree/stat_dev.h"
#endif//BRUTE_FORCE_LOCALIZATION
#ifndef SERIALIZED_EXECUTION
#include "../tree/geo_dev.h"
#include "../tree/let_dev.h"
#endif//SERIALIZED_EXECUTION
//-------------------------------------------------------------------------
#include "../time/adv_dev.h"
//-------------------------------------------------------------------------
/* #define LEAP_FROG_INTEGRATOR */
//-------------------------------------------------------------------------
#   if  defined(LEAP_FROG_INTEGRATOR) && defined(BLOCK_TIME_STEP)
#undef          LEAP_FROG_INTEGRATOR
#endif//defined(LEAP_FROG_INTEGRATOR) && defined(BLOCK_TIME_STEP)
//-------------------------------------------------------------------------
#ifndef WS93_MAC
real theta2;
#endif//WS93_MAC
//-------------------------------------------------------------------------
#   if  defined(WALK_TREE_COMBINED_MODEL) || defined(WALK_TREE_TOTAL_SUM_MODEL) || defined(WALK_TREE_GEOMETRIC_PROGRESSION_MODEL) || defined(COUNT_INTERACTIONS)
#include <math.h>
#endif//defined(WALK_TREE_COMBINED_MODEL) || defined(WALK_TREE_TOTAL_SUM_MODEL) || defined(WALK_TREE_GEOMETRIC_PROGRESSION_MODEL) || defined(COUNT_INTERACTIONS)
//-------------------------------------------------------------------------
#   if  defined(HUNT_MAKE_PARAMETER) || (defined(HUNT_FIND_PARAMETER) && defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES))
int treeBuildCalls = 0;
#endif//defined(HUNT_MAKE_PARAMETER) || (defined(HUNT_FIND_PARAMETER) && defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES))
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void setMultipoleMoment
(const int bottomLev, const soaTreeCell cell_dev, const int numNode, const soaTreeNode node_dev
 , const int num, const iparticle ibody, const soaMakeTreeBuf buf
#ifdef  CALC_MULTIPOLE_ON_DEVICE
 , const deviceProp devProp
#   if  !defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
 , const soaTreeNode node_hst, real *bmax_root
#endif//!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
#ifdef  COUNT_INTERACTIONS
 , const soaTreeCell cell_hst
#endif//COUNT_INTERACTIONS
#else///CALC_MULTIPOLE_ON_DEVICE
 , const soaTreeNode node_hst, real *bmax_root, const soaTreeCell cell_hst
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
 , const real eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifndef SERIALIZED_EXECUTION
#ifdef  CARE_EXTERNAL_PARTICLES
 , domainLocation *location
#endif//CARE_EXTERNAL_PARTICLES
 , measuredTime *elapsed
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
 , tree_stats *level
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate multipole moment of tree nodes */
  //-----------------------------------------------------------------------
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
  /* calculate multipole moment on device */
#ifndef SERIALIZED_EXECUTION
  double tmac = 0.0;
#endif//SERIALIZED_EXECUTION
  calcMultipole_dev(bottomLev, cell_dev,
		    num, ibody, numNode, node_dev, buf, devProp
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
		    , eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifndef SERIALIZED_EXECUTION
#ifndef BUILD_LET_ON_DEVICE
		    , node_hst, bmax_root
#endif//BUILD_LET_ON_DEVICE
#ifdef  CARE_EXTERNAL_PARTICLES
		    , location
#endif//CARE_EXTERNAL_PARTICLES
		    , &tmac
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
		    , cell_hst, level
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
		    , execTime
#endif//EXEC_BENCHMARK
		    );
#ifndef SERIALIZED_EXECUTION
  elapsed->sum_excg    += tmac;
  elapsed->sum_rebuild += tmac;
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
#else///CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
  /* calculate multipole moment on host */
  setRealParticles
    (num, node_hst.jtag,
#ifdef  BLOCK_TIME_STEP
     ibody.jpos,
#else///BLOCK_TIME_STEP
     ibody.pos,
#endif//BLOCK_TIME_STEP
     numNode, node_hst.jpos, node_hst.mj, node_hst.bmax
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
     , eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  WS93_MAC
     , node_hst.mr2
#endif//WS93_MAC
#ifdef  EXEC_BENCHMARK
     , execTime
#endif//EXEC_BENCHMARK
     );
  int fail_hst = 0;
  calcMultipole
    (list_hst, cell_hst, leaf_hst,
#ifdef  BLOCK_TIME_STEP
     ibody.jpos,
#else///BLOCK_TIME_STEP
     ibody.pos,
#endif//BLOCK_TIME_STEP
     cell_hst.ptag, node_hst.more, node_hst.node2cell, node_hst.jpos, node_hst.mj, node_hst.bmax,
     buf.more0, buf.more1, buf.rjmax, &buf.fail
#ifdef  WS93_MAC
     , node_hst.mr2
#endif//WS93_MAC
#ifdef  COUNT_INTERACTIONS
     , level
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
     , execTime
#endif//EXEC_BENCHMARK
     );
  *bmax_root = node_hst.bmax[0];
  if( buf_hst.fail != 0 ){
    __KILL__(stderr, "ERROR: Nc_buf exceeds bufSize of %d at least %d times.\nPLEASE re-simulate after increasing NUM_ALLOC_MACBUF defined in src/tree/make.h.\n", NUM_ALLOC_MACBUF, buf_hst.fail);
  }
  //-----------------------------------------------------------------------
  /* copy tree nodes from host to device */
  setTreeNode_dev((size_t)numNode, node_dev.more, node_hst.more,
		  node_dev.jpos, node_dev.mj,
		  node_hst.jpos, node_hst.mj
#ifdef  EXEC_BENCHMARK
		  , execTime
#endif//EXEC_BENCHMARK
		  );
  //-----------------------------------------------------------------------
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
static inline void updateDomain
(int *num, const int num_max, int *Ni,
#ifdef  EXCHANGE_USING_GPUS
 iparticle *ibody0_dev, iparticle *ibody1_dev, iparticle *ibody0_hst, iparticle *ibody1_hst, const ulong Ntot,
 particlePos pos_hst, particlePos pos_dev, domainCfg domain, domainDecomposeKey key,
 sampling sample, samplePos loc, samplePos ful, soaPHsort soa, const deviceProp devProp, const deviceInfo devInfo,
 autoTuningParam *exchangeParam, double *exchangeInterval, MPIinfo orm[restrict], MPIinfo rep[restrict],
#else///EXCHANGE_USING_GPUS
 iparticle ibody0_dev, iparticle ibody0_hst, iparticle ibody1_hst,
 const int ndim, MPIinfo *orbCfg, domainCfg *domCfg, domainDecomposeKey *domDecKey,
 const real samplingRate, const int sampleNumMax, real *sampleLoc, real *sampleFul, int *sampleRecvNum, int *sampleRecvDsp, real *domainMin, real *domainMax,
#endif//EXCHANGE_USING_GPUS
#ifdef  CARE_EXTERNAL_PARTICLES
 domainLocation *location,
#endif//CARE_EXTERNAL_PARTICLES
 const MPIcfg_tree letcfg, sendCfg *iparticleSendBuf, recvCfg *iparticleRecvBuf, measuredTime *measured
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
 , brentStatus *status, brentMemory *memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#ifdef  EXCHANGE_USING_GPUS
  /* clear counters */
  *exchangeInterval = 0.0;
  initStatVal(&((*exchangeParam).linearStats));  initGuessTime(&((*exchangeParam).linearGuess));
  initStatVal(&((*exchangeParam). powerStats));  initGuessTime(&((*exchangeParam). powerGuess));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  initStatVal  (&((*exchangeParam).parabolicStats));
  initGuessTime(&((*exchangeParam).parabolicGuess));
#endif//USE_PARABOLIC_GROWTH_MODEL
  /* execute particle exchanging */
  exchangeParticles_dev
    (*num, Ntot, num_max, num, ibody0_dev, ibody1_dev, ibody0_hst, ibody1_hst,
     pos_hst, pos_dev, key, iparticleSendBuf, iparticleRecvBuf,
     orm, rep, domain, letcfg,
     measured->sum_excg,
/* #ifdef  MONITOR_LETGEN_TIME */
/*      measured->genTree + measured->calcAcc + measured->calcMAC + measured->makeLET, */
/* #else///MONITOR_LETGEN_TIME */
/*      measured->genTree + measured->calcAcc + measured->calcMAC, */
/* #endif//MONITOR_LETGEN_TIME */
     sample, loc, ful, soa, devProp, devInfo, measured
#ifdef  CARE_EXTERNAL_PARTICLES
     , location
#endif//CARE_EXTERNAL_PARTICLES
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
     , status, memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  EXEC_BENCHMARK
     , execTime
#endif//EXEC_BENCHMARK
     );
#else///EXCHANGE_USING_GPUS
  exchangeParticles(*num, num_max, num, ibody0_dev, ibody0_hst, ibody1_hst,
		    domDecKey, iparticleSendBuf, iparticleRecvBuf,
		    samplingRate, sampleNumMax,
#ifdef  MONITOR_LETGEN_TIME
		    measured->genTree + measured->calcAcc + measured->calcMAC + measured->makeLET,
#else///MONITOR_LETGEN_TIME
		    measured->genTree + measured->calcAcc + measured->calcMAC,
#endif//MONITOR_LETGEN_TIME
		    sampleLoc, sampleFul, sampleRecvNum, sampleRecvDsp, domainMin, domainMax,
		    ndim, orbCfg, domCfg, letcfg, measured
#ifdef  CARE_EXTERNAL_PARTICLES
		    , location
#endif//CARE_EXTERNAL_PARTICLES
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
		    , status, memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  EXEC_BENCHMARK
		    , execTime
#endif//EXEC_BENCHMARK
		    );
#endif//EXCHANGE_USING_GPUS
  *Ni = *num;
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
#endif//SERIALIZED_EXECUTION
//-------------------------------------------------------------------------
static inline void buildTreeStructure
(int num,
#   if  !defined(MAKE_TREE_ON_DEVICE) || !defined(GENERATE_PHKEY_ON_DEVICE)
 const soaPHsort soaPH_hst,
#endif//!defined(MAKE_TREE_ON_DEVICE) || !defined(GENERATE_PHKEY_ON_DEVICE)
#ifdef  GENERATE_PHKEY_ON_DEVICE
 iparticle *ibody0_dev, iparticle *ibody1_dev, const soaPHsort soaPH_dev, const deviceProp devProp,
#ifdef  CUB_AVAILABLE
 const soaPHsort soaPH_pre,
#endif//CUB_AVAILABLE
#else///GENERATE_PHKEY_ON_DEVICE
 iparticle *ibody0_hst, iparticle *ibody1_hst, const iparticle ibody_dev, identity *tag,
#endif//GENERATE_PHKEY_ON_DEVICE
 int *leafLev, int *numCell, int *numNode,
#ifdef  MAKE_TREE_ON_DEVICE
 int *leafLev_dev, int *scanNum_dev, int *numCell_dev, int *numNode_dev, const soaMakeTreeBuf buf
#else///MAKE_TREE_ON_DEVICE
 int *remCell, const soaTreeCell cell_hst, const soaTreeNode node_hst
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
 , const soaTreeCell cell_dev, const soaTreeNode node_dev
#endif//defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
   , domainLocation *location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
 , struct timeval *start
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#   if  defined(HUNT_MAKE_PARAMETER) || (defined(HUNT_FIND_PARAMETER) && defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES))
  treeBuildCalls++;
#endif//defined(HUNT_MAKE_PARAMETER) || (defined(HUNT_FIND_PARAMETER) && defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES))
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* sort N-body particles in Peano-Hilbert order */
#ifdef  GENERATE_PHKEY_ON_DEVICE
  sortParticlesPHcurve_dev(num, ibody0_dev, ibody1_dev, soaPH_dev, devProp
#ifdef  CUB_AVAILABLE
			   , soaPH_pre
#endif//CUB_AVAILABLE
#ifndef MAKE_TREE_ON_DEVICE
			   , soaPH_hst
#endif//MAKE_TREE_ON_DEVICE
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
			   , location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
			   , start
#ifdef  EXEC_BENCHMARK
			   , execTime
#endif//EXEC_BENCHMARK
			   );
#else///GENERATE_PHKEY_ON_DEVICE
  sortParticlesPHcurve(num, ibody0_hst, ibody1_hst, tag, peano
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
		       , location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
		       , start
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );
#endif//GENERATE_PHKEY_ON_DEVICE
  //-----------------------------------------------------------------------
  /* build a breadth-first octree structure */
#ifdef  MAKE_TREE_ON_DEVICE
  makeTreeStructure_dev(num, soaPH_dev.key, leafLev, leafLev_dev, scanNum_dev, numCell, numCell_dev, cell_dev, numNode, numNode_dev, node_dev, buf, devProp
#ifdef  EXEC_BENCHMARK
			, execTime
#endif//EXEC_BENCHMARK
			);
#else///MAKE_TREE_ON_DEVICE
  makeTree(cell_hst.level, numCell, remCell, cell_hst.cell, cell_hst.hkey, cell_hst.parent, cell_hst.children, cell_hst.leaf,
	   num, soaPH_hst.key, numNode, cell_hst.ptag, node_hst.more, node_hst.jtag, node_hst.node2cell
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	   , node_hst.niSub
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#ifdef  EXEC_BENCHMARK
	   , execTime
#endif//EXEC_BENCHMARK
	   );
  for(int levelIdx = (NUM_PHKEY_LEVEL - 1); levelIdx >= 0; levelIdx--)
    if( cell_hst.level[levelIdx].num != 0 ){
      *leafLev = levelIdx;
      break;
    }
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
  /* copy tree cells from host to device if necessary */
#   if  defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(MAKE_TREE_ON_DEVICE)
  setTreeNode_dev((size_t)(*numNode), node_dev, node_hst, num
#ifdef  EXEC_BENCHMARK
		  , execTime
#endif//EXEC_BENCHMARK
		  );
  setTreeCell_dev(*numCell, cell_dev, cell_hst
#ifdef  EXEC_BENCHMARK
		  , execTime
#endif//EXEC_BENCHMARK
		  );
#endif//defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(MAKE_TREE_ON_DEVICE)
  //-----------------------------------------------------------------------
  /* set N-body particles on device */
#ifndef GENERATE_PHKEY_ON_DEVICE
  copyParticle_hst2dev(num, *ibody0_hst, ibody_dev
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );
#endif//GENERATE_PHKEY_ON_DEVICE
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void configDistribution
(const int num, const int inumPerLane, int *Ngrp, const int maxNgrp, laneinfo *laneInfo_hst, laneinfo *laneInfo_dev
#   if  defined(BLOCK_TIME_STEP) || (defined(LOCALIZE_I_PARTICLES) && defined(BRUTE_FORCE_LOCALIZATION))
 , iparticle ibody_dev
#endif//defined(BLOCK_TIME_STEP) || (defined(LOCALIZE_I_PARTICLES) && defined(BRUTE_FORCE_LOCALIZATION))
#ifdef  LOCALIZE_I_PARTICLES
#ifdef  USE_BRENT_METHOD
 , brentStatus *brent
#else///USE_BRENT_METHOD
#ifdef  CUB_AVAILABLE
 , soaCUBreal util
#endif//CUB_AVAILABLE
#endif//USE_BRENT_METHOD
#       ifdef  BRUTE_FORCE_LOCALIZATION
 , int *inum_dev, int *inum_hst
#ifndef FACILE_NEIGHBOR_SEARCH
 , const soaTreeCell cell, const soaTreeNode node, const soaMakeTreeBuf makeBuf, const soaNeighborSearchBuf searchBuf, deviceProp devProp
#endif//FACILE_NEIGHBOR_SEARCH
#       else///BRUTE_FORCE_LOCALIZATION
 , PHint *peano
#       endif//BRUTE_FORCE_LOCALIZATION
#endif//LOCALIZE_I_PARTICLES
#ifdef  BLOCK_TIME_STEP
 /* , real *laneTime_dev */
 , double *laneTime_dev
#       ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
 , real *dt_dev, real *dt_hst, int *dtInfo_num, histogram_dt **dtInfo
#       endif//USE_VARIABLE_NEIGHBOR_LEVEL
#endif//BLOCK_TIME_STEP
 , const struct timeval start, measuredTime *time
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
  static bool initialized = false;
  if( initialized )    brentCalc1st(brent, 1.0e-13);
  else                 initialized = true;
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
  //-----------------------------------------------------------------------
  updateParticleGroups(num, laneInfo_hst, inumPerLane, maxNgrp, Ngrp
#ifdef  LOCALIZE_I_PARTICLES
#       ifdef  BRUTE_FORCE_LOCALIZATION
		       , ibody_dev, inum_dev, inum_hst
#ifdef  USE_BRENT_METHOD
		       , (real)brent->u.pos
#else///USE_BRENT_METHOD
#ifdef  CUB_AVAILABLE
		       , util
#endif//CUB_AVAILABLE
#endif//USE_BRENT_METHOD
#ifndef FACILE_NEIGHBOR_SEARCH
		       , cell, node, makeBuf, searchBuf, devProp
#endif//FACILE_NEIGHBOR_SEARCH
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
#       else///BRUTE_FORCE_LOCALIZATION
		       , peano
#       endif//BRUTE_FORCE_LOCALIZATION
#endif//LOCALIZE_I_PARTICLES
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
		       , ibody_dev, dt_dev, dt_hst, dtInfo_num, dtInfo
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
		       );
  //-----------------------------------------------------------------------
  commitParticleGroups(*Ngrp, laneInfo_hst, laneInfo_dev);
  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
  setLaneTime_dev(*Ngrp, laneInfo_dev, laneTime_dev, ibody_dev
#ifdef  EXEC_BENCHMARK
		  , execTime
#endif//EXEC_BENCHMARK
		  );
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
  static struct timeval finish;
  gettimeofday(&finish, NULL);
  time->makeTree = calcElapsedTimeInSec(start, finish);
#ifndef SERIALIZED_EXECUTION
  time->sum_excg    += time->makeTree;
  time->sum_rebuild += time->makeTree;
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef EXEC_BENCHMARK
//-------------------------------------------------------------------------
static inline void dumpSnapshot
(const int unit, const int num, iparticle ibody_dev, iparticle ibody_hst,
#ifdef  USE_HDF5_FORMAT
 nbody_hdf5 *body, hdf5struct hdf5type,
#ifdef  MONITOR_ENERGY_ERROR
 energyError *relEneErr,
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
 const double time, const ulong steps, char *file, const uint present, uint *previous
#ifdef  LEAP_FROG_INTEGRATOR
 , const double dt
#endif//LEAP_FROG_INTEGRATOR
#ifndef SERIALIZED_EXECUTION
 , MPIcfg_dataio *iocfg, const ulong Ntot
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
#ifdef  COMPARE_WITH_DIRECT_SOLVER
 , const int Ni, const iparticle ibody_direct_dev, const deviceProp devProp, acceleration *direct_hst, const int Ngrp, laneinfo *laneInfo_dev
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
 , const real eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
 , char *accfile
#endif//COMPARE_WITH_DIRECT_SOLVER
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set particle data on host */
  copyParticle_dev2hst(num, ibody_dev, ibody_hst
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  copyAry4Snapshot(0, num, ibody_hst, *body);
#       ifdef  LEAP_FROG_INTEGRATOR
  backVel4Snapshot(0, num,            *body, (real)dt);
#       endif//LEAP_FROG_INTEGRATOR
#else///USE_HDF5_FORMAT
#       ifdef  LEAP_FROG_INTEGRATOR
  backVel(0, num, ibody_hst, (real)dt);
#       endif//LEAP_FROG_INTEGRATOR
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  /* output the snapshot file */
#ifdef  SERIALIZED_EXECUTION
  writeSnapshot        (unit, time, steps, num,
#ifdef  USE_HDF5_FORMAT
			body,
#else///USE_HDF5_FORMAT
			ibody_hst,
#endif//USE_HDF5_FORMAT
			file, present
#ifdef  USE_HDF5_FORMAT
			, hdf5type
#ifdef  MONITOR_ENERGY_ERROR
			, relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
			);
#else///SERIALIZED_EXECUTION
  writeSnapshotParallel(unit, time, steps, num,
#ifdef  USE_HDF5_FORMAT
			body,
#else///USE_HDF5_FORMAT
			ibody_hst,
#endif//USE_HDF5_FORMAT
			file, present, iocfg, Ntot
#ifdef  USE_HDF5_FORMAT
			, hdf5type
#ifdef  MONITOR_ENERGY_ERROR
			, relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
			);
#endif//SERIALIZED_EXECUTION
  *previous = present;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  COMPARE_WITH_DIRECT_SOLVER
  //-----------------------------------------------------------------------
  /* output exact and approximated acceleration and potential */
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.dat", DATAFOLDER, file, BRUTE, present);
  if( 0 != access(filename, F_OK) ){
    //---------------------------------------------------------------------
    static soaTreeNode dummy_node;
    static soaTreeWalkBuf dummy_buf;
    double null;
    //---------------------------------------------------------------------
    calcGravity_dev(
#ifdef  BLOCK_TIME_STEP
		    Ngrp, NULL,
#endif//BLOCK_TIME_STEP
		    Ngrp, laneInfo_dev, Ni, ibody_direct_dev, dummy_node, dummy_buf, NULL, devProp
#ifndef SERIALIZED_EXECUTION
		    , NULL, 0, 1, NULL, letcfg, NULL, NULL, 0, NULL, NULL, NULL
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
		    , NULL
#endif//COUNT_INTERACTIONS
		    , &null
#ifdef  EXEC_BENCHMARK
		    , NULL
#endif//EXEC_BENCHMARK
#ifdef  COMPARE_WITH_DIRECT_SOLVER
		    , false
#endif//COMPARE_WITH_DIRECT_SOLVER
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
		    , eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
		    );
    //---------------------------------------------------------------------
    copyAccel_dev2hst(Ni, ibody_direct_dev.acc, direct_hst);
    writeApproxAccel(time, steps, num,
#ifdef  USE_HDF5_FORMAT
		     (*body).idx,
#endif//USE_HDF5_FORMAT
		     direct_hst, filename);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.%s.%.3u.dat", DATAFOLDER, file, accfile, present);
  writeApproxAccel(time, steps, num,
#ifdef  USE_HDF5_FORMAT
		   (*body).idx,
#endif//USE_HDF5_FORMAT
		   ibody_hst.acc, filename);
  //-----------------------------------------------------------------------
#endif//COMPARE_WITH_DIRECT_SOLVER
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifndef COMPARE_WITH_DIRECT_SOLVER
static inline void dumpRestarter
(const int num, iparticle ibody_dev, iparticle ibody_hst, const double time, const double dt, const ulong steps, char *file, int *last, double *formerTime
#ifdef  USE_HDF5_FORMAT
 , hdf5struct hdf5type
 , rebuildTree rebuild, measuredTime measured
#ifdef  WALK_TREE_COMBINED_MODEL
 , autoTuningParam rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
 , brentStatus status, brentMemory memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
 , energyError relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
#ifndef SERIALIZED_EXECUTION
 , MPIcfg_dataio *iocfg, const MPIinfo mpi, const ulong Ntot
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output tentative results of the current running simulation */
  //-----------------------------------------------------------------------
  /* set particle data on host */
  copyParticle_dev2hst(num, ibody_dev, ibody_hst
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );
#ifdef  LEAP_FROG_INTEGRATOR
  backVel(0, num, body, (real)dt);
#endif//LEAP_FROG_INTEGRATOR
  //-----------------------------------------------------------------------
  /* output the dump file */
#ifdef  SERIALIZED_EXECUTION
  writeTentativeData        (time, dt, steps, num, ibody_hst, file, last
#ifdef  USE_HDF5_FORMAT
			     , hdf5type
			     , rebuild, measured
#ifdef  WALK_TREE_COMBINED_MODEL
			     , rebuildParam
/* 			     , linearStats, linearGuess, powerStats, powerGuess */
/* #ifdef  USE_PARABOLIC_GROWTH_MODEL */
/* 			     , parabolicStats, parabolicGuess */
/* #endif//USE_PARABOLIC_GROWTH_MODEL */
#endif//WALK_TREE_COMBINED_MODEL
#ifdef  USE_BRENT_METHOD
			     , status, memory
#endif//USE_BRENT_METHOD
#ifdef  MONITOR_ENERGY_ERROR
			     , relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
			     );
#else///SERIALIZED_EXECUTION
  writeTentativeDataParallel(time, dt, steps, num, ibody_hst, file, last, iocfg, Ntot
#ifdef  USE_HDF5_FORMAT
			     , hdf5type
			     , rebuild, measured
#ifdef  WALK_TREE_COMBINED_MODEL
			     , rebuildParam
/* 			     , linearStats, linearGuess, powerStats, powerGuess */
/* #ifdef  USE_PARABOLIC_GROWTH_MODEL */
/* 			     , parabolicStats, parabolicGuess */
/* #endif//USE_PARABOLIC_GROWTH_MODEL */
#endif//WALK_TREE_COMBINED_MODEL
#ifdef  USE_BRENT_METHOD
			     , status, memory
#endif//USE_BRENT_METHOD
#ifdef  MONITOR_ENERGY_ERROR
			     , relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
			     );
  if( mpi.rank == 0 )
#endif//SERIALIZED_EXECUTION
    {
      updateConfigFile(*last, file);
      *formerTime = getPresentTimeInMin();
    }
#ifndef SERIALIZED_EXECUTION
  chkMPIerr(MPI_Bcast(formerTime, 1, MPI_DOUBLE, 0, mpi.comm));
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
#endif//COMPARE_WITH_DIRECT_SOLVER
//-------------------------------------------------------------------------
#endif//EXEC_BENCHMARK
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef COUNT_INTERACTIONS
void analyzeWalkStatistics(int Ni, const int Ngrp, int * restrict results, walk_stats *stats);
void analyzeTreeMetrics(tree_metrics *metric);
#endif//COUNT_INTERACTIONS
//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //-----------------------------------------------------------------------
  /* parallelized region employing MPI start */
#ifndef SERIALIZED_EXECUTION
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);
  static MPIcfg_dataio iocfg;
  createMPIcfg_dataio(&iocfg, mpi);
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
#ifdef  REPORT_TOTAL_ELAPSED_TIME
  static struct timeval timeInit;
  gettimeofday(&timeInit, NULL);
#endif//REPORT_TOTAL_ELAPSED_TIME
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* configure the details of the numerical simulation */
  //-----------------------------------------------------------------------
  /* read command line arguments */
  if( argc < 4 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 4);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -jobID=<int>\n");
#ifdef  GADGET_MAC
    __FPRINTF__(stderr, "          -absErr=<real>\n");
#else///GADGET_MAC
#ifdef  WS93_MAC
    __FPRINTF__(stderr, "          -accErr=<real>\n");
#else///WS93_MAC
    __FPRINTF__(stderr, "          -theta=<real>\n");
#endif//WS93_MAC
#endif//GADGET_MAC
    __FPRINTF__(stderr, "          -dropPrevTune=<int> (optional)\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }
  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv,  "file", &file));
  double tmp;
#ifdef  GADGET_MAC
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "absErr", &tmp));  const real absErr = (real)tmp;
#else///GADGET_MAC
#ifdef  WS93_MAC
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "accErr", &tmp));  const real accErr = (real)tmp;
#else///WS93_MAC
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv,  "theta", &tmp));  const real  theta = (real)tmp;
  theta2 = theta * theta;
#endif//WS93_MAC
#endif//GADGET_MAC
#ifdef  EXEC_BENCHMARK
  static int jobID;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "jobID", &jobID));
#endif//EXEC_BENCHMARK
#if 1
  static int dropPrevTune;
  if( optionalCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "dropPrevTune", &dropPrevTune)) != myUtilAvail )
    dropPrevTune = 0;
#endif
  //-----------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
  static int Nx;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "Nx", &Nx));
  static int Ny;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "Ny", &Ny));
  static int Nz;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "Nz", &Nz));
  if( (Nx * Ny * Nz) != mpi.size ){
    __KILL__(stderr, "ERROR: Invalid input arguments: Nx = %d, Ny = %d, Nz = %d while mpi.size = %d\n", Nx, Ny, Nz, mpi.size);
  }/* if( (Nx * Ny * Nz) != mpi.size ){ */
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
  /* read settings about the simulation */
  static ulong Ntot;
  static real eps, eta;
  static double ft, snapshotInterval, saveInterval;
  static int unit;
#ifdef  SERIALIZED_EXECUTION
  readSettings        (&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file);
#else///SERIALIZED_EXECUTION
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file, mpi);
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
  /* read setting to dump tentative results of the simulation */
  static int last;
#ifdef  SERIALIZED_EXECUTION
  readConfigFile        (&last, file);
#else///SERIALIZED_EXECUTION
  readConfigFileParallel(&last, file, mpi);
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* activate the attatched accelerator device(s) */
  //-----------------------------------------------------------------------
  int devIdx;
  deviceInfo devInfo;
  deviceProp devProp;
#ifdef  SERIALIZED_EXECUTION
  openAcceleratorGPUs(&devIdx, &devInfo, &devProp, 0, 1, 0, 1);
#else///SERIALIZED_EXECUTION
  openAcceleratorGPUs(&devIdx, &devInfo, &devProp, 0, 1, mpi.rank, mpi.size);
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
  /* set CUDA streams */
  cudaStream_t *stream;
  kernelStream sinfo;
  setCUDAstreams_dev(&stream, &sinfo, &devInfo, &devProp);
  //-----------------------------------------------------------------------
  /* set CUDA events */
  /* maximum # of CUDA events = # of GPUs for tree traversal, # of GPUs - 1 for LET construction */
  /* this code assumes that # of GPUs = # of MPI processes (1 GPU / 1 MPI process) */
#   if  defined(USE_CUDA_EVENT) && !defined(SERIALIZED_EXECUTION)
  cudaEvent_t *iniCalcAcc, *finCalcAcc;
#ifdef  MONITOR_LETGEN_TIME
  cudaEvent_t *iniMakeLET, *finMakeLET;
#endif//MONITOR_LETGEN_TIME
  allocateCUDAevents_dev(&iniCalcAcc, &finCalcAcc
#ifdef  MONITOR_LETGEN_TIME
			, &iniMakeLET, &finMakeLET
#endif//MONITOR_LETGEN_TIME
			 , mpi.size);
#else///defined(USE_CUDA_EVENT) && !defined(SERIALIZED_EXECUTION)
#   if  defined(USE_CUDA_EVENT) && defined(PRINT_PSEUDO_PARTICLE_INFO)
  cudaEvent_t *iniCalcAcc, *finCalcAcc;
  allocateCUDAevents_dev(&iniCalcAcc, &finCalcAcc, 1);
#endif//defined(USE_CUDA_EVENT) && defined(PRINT_PSEUDO_PARTICLE_INFO)
#endif//defined(USE_CUDA_EVENT) && !defined(SERIALIZED_EXECUTION)
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set the number of N-body particles */
  /* this procedure is obvious for this serial process execution, but useful for the future upgrade (parallelization) */
  //-----------------------------------------------------------------------
  /* assume num >= Ni */
  static int num, Ni;
  //-----------------------------------------------------------------------
#ifdef  SERIALIZED_EXECUTION
  num = Ni = (int)Ntot;
  const int num_max = num;
#else///SERIALIZED_EXECUTION
  static MPIcfg_tree letcfg;
  setNodeConfig(Ntot, &num, &Ni, mpi, &letcfg, devIdx);
  /* const int num_max = num * MAX_FACTOR_FROM_EQUIPARTITION; */
  const int num_max = (int)ceilf((float)num * MAX_FACTOR_FROM_EQUIPARTITION);
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set global constants */
  //-----------------------------------------------------------------------
  /* set gravitational softening and opening criterion */
  extern real newton;
  setGlobalConstants_walk_dev_cu
    (newton, eps * eps
#ifndef WS93_MAC
     , theta2
#endif//WS93_MAC
     );
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  setGlobalConstants_make_dev_cu(
#ifdef  GADGET_MAC
				 absErr, newton
#else///GADGET_MAC
#     ifdef  WS93_MAC
                                 accErr
#     endif//WS93_MAC
#endif//GADGET_MAC
				 );
#else///CALC_MULTIPOLE_ON_DEVICE
#     ifdef  WS93_MAC
  setGlobalConstants_make_c     (UNITY / (real)(1.0e-30 + (double)accErr));
#     endif//WS93_MAC
#endif//CALC_MULTIPOLE_ON_DEVICE
#   if  !defined(WS93_MAC) && !defined(GADGET_MAC)
  setGlobalConstants_macutil_c(theta);
#endif//!defined(WS93_MAC) && !defined(GADGET_MAC)
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && defined(SMEM_PREF_FOR_NEIGHBOR_SEARCH)
  setGlobalConstants_neighbor_dev_cu();
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && defined(SMEM_PREF_FOR_NEIGHBOR_SEARCH)
#ifndef SERIALIZED_EXECUTION
  setGlobalConstants_let_dev_cu(
#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
				theta2
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
				);
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation */
  //-----------------------------------------------------------------------
  /* declaration of array to contain whole information of whole N-body particles */
#ifdef  USE_HDF5_FORMAT
  nbody_hdf5 body_snapshot;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
  ulong *hdf5_idx;
  const muse alloc_snap = allocSnapshotArray(&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx, num_max, &body_snapshot);
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  /* declarations of arrays to contain Peano--Hilbert key of N-body particles */
#ifndef MAKE_TREE_ON_DEVICE
  PHinfo *list;
#endif//MAKE_TREE_ON_DEVICE
  PHint *peano;
#ifdef  GENERATE_PHKEY_ON_DEVICE
  PHint *peano_dev;
  int   *  tag_dev;
  float4 *min_dev, *max_dev;
#ifndef SERIALIZED_EXECUTION
  float4 *box_min, *box_max, *min_hst, *max_hst;
#endif//SERIALIZED_EXECUTION
  int *gsync_ph0, *gsync_ph1;
  soaPHsort soaPH_dev, soaPH_hst;
#ifdef  CUB_AVAILABLE
  soaPHsort soaPH_pre;
  PHint *peano_pre;
  int   *  tag_pre;
  void *phsort_temp_storage;
#endif//CUB_AVAILABLE
  const muse alloc_phkey = allocPeanoHilbertKey_dev
    (num_max, &tag_dev, &peano_dev, &peano, &min_dev, &max_dev, &gsync_ph0, &gsync_ph1,
#ifndef SERIALIZED_EXECUTION
    &box_min, &box_max, &min_hst, &max_hst,
#endif//SERIALIZED_EXECUTION
#ifndef CALC_MULTIPOLE_ON_DEVICE
     &list,
#endif//CALC_MULTIPOLE_ON_DEVICE
     &soaPH_dev, &soaPH_hst,
#ifdef  CUB_AVAILABLE
     &soaPH_pre, &phsort_temp_storage, &tag_pre, &peano_pre,
#endif//CUB_AVAILABLE
     devProp
     );
#ifndef CALC_MULTIPOLE_ON_DEVICE
  soaPH_hst.level = list;
#endif//CALC_MULTIPOLE_ON_DEVICE
#else///GENERATE_PHKEY_ON_DEVICE
  identity *tag;
  const muse alloc_phkey = allocPeanoHilbertKey
    (num_max, &peano, &tag
#ifndef CALC_MULTIPOLE_ON_DEVICE
     , &list
#endif//CALC_MULTIPOLE_ON_DEVICE
     );
#endif//GENERATE_PHKEY_ON_DEVICE
  //-----------------------------------------------------------------------
  /* declarations of arrays to contain information of i-particles */
  iparticle     ibody0,  ibody1,  ibody0_dev;
  ulong        *  idx0, *  idx1, *  idx0_dev;
  position     *  pos0, *  pos1, *  pos0_dev;
  acceleration *  acc0, *  acc1, *  acc0_dev;
#ifdef  GENERATE_PHKEY_ON_DEVICE
  iparticle                       ibody1_dev;
  ulong                          *  idx1_dev;
  position                       *  pos1_dev;
  acceleration                   *  acc1_dev;
#endif//GENERATE_PHKEY_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
  real *neighbor_dev;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  position *encBall_dev, *encBall_hst;
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  DPADD_FOR_ACC
  DPacc *tmp_dev;
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
  acceleration *res_dev;
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
#ifdef  BLOCK_TIME_STEP
  velocity     *  vel0, *  vel1, *  vel0_dev;
  ibody_time   * time0, * time1, * time0_dev;
#       ifdef  GENERATE_PHKEY_ON_DEVICE
  velocity                       *  vel1_dev;
  ibody_time                     * time1_dev;
#       else///GENERATE_PHKEY_ON_DEVICE
  position *jpos_dev;
  velocity *jvel_dev;
#       endif//GENERATE_PHKEY_ON_DEVICE
  const muse alloc_ibody0    = allocParticleDataSoA_hst(num_max, &ibody0    , &idx0    , &pos0    , &acc0    , &vel0    , &time0    );
  const muse alloc_ibody1    = allocParticleDataSoA_hst(num_max, &ibody1    , &idx1    , &pos1    , &acc1    , &vel1    , &time1    );
  const muse alloc_ibody_dev = allocParticleDataSoA_dev(num_max
							, &ibody0_dev, &idx0_dev, &pos0_dev, &acc0_dev, &vel0_dev, &time0_dev
#       ifdef  GENERATE_PHKEY_ON_DEVICE
							, &ibody1_dev, &idx1_dev, &pos1_dev, &acc1_dev, &vel1_dev, &time1_dev
#       else///GENERATE_PHKEY_ON_DEVICE
							, &jpos_dev, &jvel_dev
#       endif//GENERATE_PHKEY_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
							, &neighbor_dev
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
							, &encBall_dev, &encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  DPADD_FOR_ACC
							, &tmp_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
							, &res_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
							);
#   if  !defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(GENERATE_PHKEY_ON_DEVICE)
  ibody0.jpos = pos1;  ibody0.jvel = vel1;
  ibody1.jpos = pos0;  ibody1.jvel = vel0;
#endif//!defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(GENERATE_PHKEY_ON_DEVICE)
#else///BLOCK_TIME_STEP
  real *vx0, *vx1, *vx0_dev;
  real *vy0, *vy1, *vy0_dev;
  real *vz0, *vz1, *vz0_dev;
#ifdef  GENERATE_PHKEY_ON_DEVICE
  real             *vx1_dev;
  real             *vy1_dev;
  real             *vz1_dev;
#endif//GENERATE_PHKEY_ON_DEVICE
  const muse alloc_ibody0    = allocParticleDataSoA_hst(num_max, &ibody0    , &idx0    , &pos0    , &acc0    , &vx0    , &vy0    , &vz0   );
  const muse alloc_ibody1    = allocParticleDataSoA_hst(num_max, &ibody1    , &idx1    , &pos1    , &acc1    , &vx1    , &vy1    , &vz1   );
  const muse alloc_ibody_dev = allocParticleDataSoA_dev(num_max
							, &ibody0_dev, &idx0_dev, &pos0_dev, &acc0_dev, &vx0_dev, &vy0_dev, &vz0_dev
#       ifdef  GENERATE_PHKEY_ON_DEVICE
							, &ibody1_dev, &idx1_dev, &pos1_dev, &acc1_dev, &vx1_dev, &vy1_dev, &vz1_dev
#       endif//GENERATE_PHKEY_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
							, &neighbor_dev
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
							, &encBall_dev, &encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  DPADD_FOR_ACC
							, &tmp_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
							, &res_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
							);
#endif//BLOCK_TIME_STEP
  iparticle_treeinfo treeinfo, treeinfo_dev;
#ifndef MAKE_TREE_ON_DEVICE
  int *jtag;
#endif//MAKE_TREE_ON_DEVICE
  int *jtag_dev;
#ifdef  COUNT_INTERACTIONS
  int *Nj, *Nj_dev, *Nbuf, *Nbuf_dev;
#endif//COUNT_INTERACTIONS
  const muse alloc_jtag = allocParticleInfoSoA_hst
    (num_max, &treeinfo
#ifndef MAKE_TREE_ON_DEVICE
     , &jtag
#endif//MAKE_TREE_ON_DEVICE
#ifdef  COUNT_INTERACTIONS
     , &Nj, &Nbuf
#endif//COUNT_INTERACTIONS
     );
  const muse alloc_jtag_dev = allocParticleInfoSoA_dev
    (num_max, &treeinfo_dev, &jtag_dev
#ifdef  COUNT_INTERACTIONS
     , &Nj_dev, &Nbuf_dev
#endif//COUNT_INTERACTIONS
     );
  //-----------------------------------------------------------------------
  /* declarations of variables to contain information on groups of i-particles */
  laneinfo *laneInfo_hst, *laneInfo_dev;
  /* real *laneTime_dev; */
  double *laneTime_dev;
  int inumPerLane, maxNgrp;
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
  int *inum_hst, *inum_dev;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
  real *ibody_dt_dev, *ibody_dt_hst;
  histogram_dt *dtInfo;
  int dtInfo_num;
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
#   if  defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
  soaCUBreal soaCUBneighbor;
  void *neighbor_temp_storage;
  real *neighbor_cubout;
#endif//defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
  const muse alloc_lane_dev = allocParticleGroups(&laneInfo_hst, &laneInfo_dev, &laneTime_dev,
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
						  &inum_hst, &inum_dev,
#   if  defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
						  &soaCUBneighbor, &neighbor_temp_storage, &neighbor_cubout, ibody0_dev,
#endif//defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
						  &ibody_dt_dev, &ibody_dt_hst, &dtInfo, &dtInfo_num,
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
						  &inumPerLane, &maxNgrp, num_max, devProp);
  //-----------------------------------------------------------------------
  /* declarations of arrays to contain information of tree-cells */
#ifdef  MAKE_TREE_ON_DEVICE
  PHint    *hkey_dev;
  uint     *parent_dev, *children_dev;
#else///MAKE_TREE_ON_DEVICE
  PHint    *hkey;
  uint     *parent, *children;
#endif//MAKE_TREE_ON_DEVICE
#   if  !defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(MAKE_TREE_ON_DEVICE)
  treecell *cell;
  bool     *leaf;
  uint     *node;
#endif//!defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(MAKE_TREE_ON_DEVICE)
  soaTreeCell soaCell_dev, soaCell_hst;
#ifndef MAKE_TREE_ON_DEVICE
  const muse alloc_cell = allocTreeCell(&hkey, &parent, &children
#ifndef CALC_MULTIPOLE_ON_DEVICE
					, &cell, &leaf, &node
#endif//CALC_MULTIPOLE_ON_DEVICE
					, &soaCell_hst);
#else///MAKE_TREE_ON_DEVICE
  const muse alloc_cell = {0, 0};
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
  treecell *cell_dev;
  bool     *leaf_dev;
  uint     *node_dev;
  PHinfo   *list_dev;
#ifndef MAKE_TREE_ON_DEVICE
  treecell *cell;
  bool     *leaf;
  uint     *node;
#else///MAKE_TREE_ON_DEVICE
#ifdef  COUNT_INTERACTIONS
  treecell *cell;
  bool     *leaf;
  uint     *node;
  PHinfo   *list;
#endif//COUNT_INTERACTIONS
#endif//MAKE_TREE_ON_DEVICE
#ifdef  MAKE_TREE_ON_DEVICE
  int *bottomLev_dev, *numCell_dev, *numNode_dev, *scanNum_dev;
#endif//MAKE_TREE_ON_DEVICE
  const muse alloc_cell_dev =
    allocTreeCell_dev(&cell_dev, &leaf_dev, &node_dev, &list_dev,
#ifdef  MAKE_TREE_ON_DEVICE
		      &hkey_dev, &parent_dev, &children_dev, &bottomLev_dev, &numCell_dev, &numNode_dev, &scanNum_dev,
#else///MAKE_TREE_ON_DEVICE
		      &cell    , &leaf    , &node    , &list    ,
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
		      &cell    , &leaf    , &node    , &list    ,
#endif//defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
		      devProp,
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
		      &soaCell_dev, &soaCell_hst);
#endif//defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
  //-----------------------------------------------------------------------
  /* declarations of arrays to contain information of tree-nodes */
  uint    *more_dev;
  jparticle *pj_dev;
  jmass     *mj_dev;
#ifdef CALC_MULTIPOLE_ON_DEVICE
  real *bmax_dev;
  int *node2cell_dev;
  int *gsync0, *gsync1;
#       ifdef  WS93_MAC
  real *mr2_dev;
#       endif//WS93_MAC
#else///CALC_MULTIPOLE_ON_DEVICE
  real *bmax;
#       ifdef  WS93_MAC
  real *mr2;
#       endif//WS93_MAC
#endif//CALC_MULTIPOLE_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  int *niSub_dev;
#ifndef MAKE_TREE_ON_DEVICE
  int *niSub_hst;
#endif//MAKE_TREE_ON_DEVICE
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
  uint *more;
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
  jparticle *pj;
  jmass     *mj;
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
#ifdef  MAKE_TREE_ON_DEVICE
  int *gmem_make_tree_dev, *gsync0_make_tree_dev, *gsync1_make_tree_dev, *gsync2_make_tree_dev, *gsync3_make_tree_dev;
  int *gmem_link_tree_dev, *gsync0_link_tree_dev, *gsync1_link_tree_dev;
#else///MAKE_TREE_ON_DEVICE
  int *node2cell;
#endif//MAKE_TREE_ON_DEVICE
#ifdef  GADGET_MAC
  real *mac_dev;
#endif//GADGET_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  int *gmem_external_dev, *gsync0_external_dev, *gsync1_external_dev;
  domainLocation location;
  float *diameter_dev, *diameter_hst;
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  int  *more0Buf, *more1Buf, *makeFail;
  real *rjmaxBuf;
  soaTreeNode soaNode_dev, soaNode_hst;
  soaMakeTreeBuf soaMakeBuf;
  const muse alloc_node_dev = allocTreeNode_dev(&more_dev, &pj_dev, &mj_dev,
#ifdef  CALC_MULTIPOLE_ON_DEVICE
						&bmax_dev, &node2cell_dev, &gsync0, &gsync1, devProp,
#       ifdef  WS93_MAC
						&mr2_dev,
#       endif//WS93_MAC
#else///CALC_MULTIPOLE_ON_DEVICE
						&bmax,
#       ifdef  WS93_MAC
						&mr2    ,
#       endif//WS93_MAC
#endif//CALC_MULTIPOLE_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
						&niSub_dev,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
						&more,
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
						&pj, &mj,
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
#ifdef  MAKE_TREE_ON_DEVICE
						&gmem_make_tree_dev, &gsync0_make_tree_dev, &gsync1_make_tree_dev, &gsync2_make_tree_dev, &gsync3_make_tree_dev,
						&gmem_link_tree_dev, &gsync0_link_tree_dev, &gsync1_link_tree_dev,
#else///MAKE_TREE_ON_DEVICE
						&node2cell,
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
						&niSub_hst,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#endif//MAKE_TREE_ON_DEVICE
#ifdef  GADGET_MAC
						&mac_dev,
#endif//GADGET_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
						&gmem_external_dev, &gsync0_external_dev, &gsync1_external_dev, &diameter_dev, &diameter_hst, &location, CAST_R2F(eps), CAST_R2F(eta),
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
						&more0Buf, &more1Buf, &rjmaxBuf, &makeFail, &soaNode_dev, &soaNode_hst, &soaMakeBuf);
  //-----------------------------------------------------------------------
  soaNode_dev.jtag = jtag_dev;
#ifndef MAKE_TREE_ON_DEVICE
  soaNode_hst.jtag = jtag;
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  int *gsync_ns0, *gsync_ns1;
  uint *freeLst_ns;
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
  uint *freeNum_ns;
  int *active_ns;
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
  soaNeighborSearchBuf soaNeighbor_dev;
  const muse alloc_neighbor_dev =
    allocNeighborSearch_dev(&gsync_ns0, &gsync_ns1, &freeLst_ns,
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
			    &freeNum_ns, &active_ns,
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
			    &soaNeighbor_dev, devProp);
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  //-----------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
  sendCfg *iparticleSendBuf;
  recvCfg *iparticleRecvBuf;
  int *sampleRecvNum, *sampleRecvDsp;
#ifdef  EXCHANGE_USING_GPUS
  static MPIinfo ormCfg[3], repCfg[3];
  float *dxmin, *dxmax, *dymin, *dymax, *dzmin, *dzmax;
  float *sxmin, *sxmax, *symin, *symax, *szmin, *szmax;
  MPI_Request *dmreq;
  domainCfg domCfg;
  sampling sample;
  const muse alloc_dd = allocateORMtopology(&dxmin, &dxmax, &dymin, &dymax, &dzmin, &dzmax, &dmreq,
					    &sxmin, &sxmax, &symin, &symax, &szmin, &szmax,
					    &iparticleSendBuf, &iparticleRecvBuf, &sampleRecvNum, &sampleRecvDsp,
					    ormCfg, repCfg, Nx, Ny, Nz, &letcfg, &domCfg, &sample, Ntot);
  /* float4 *pmin_hst, *pmax_hst, *pmin_dev, *pmax_dev; */
  /* int *gsync_box0, *gsync_box1; */
  /* soaBoxSize soaBox; */
  /* const muse alloc_box = allocateBoxSize_dev(&pmin_hst, &pmax_hst, &pmin_dev, &pmax_dev, &gsync_box0, &gsync_box1, &soaBox, devProp); */
  checkBoxSize_dev(devProp);
  float *x0hst, *x1hst, *y0hst, *y1hst, *z0hst, *z1hst;  int *idhst;
  float *x0dev, *x1dev, *y0dev, *y1dev, *z0dev, *z1dev;  int *iddev;
  samplePos samplePos0, samplePos1;
  const muse alloc_spl = allocateSamplePos(&x0hst, &x1hst, &y0hst, &y1hst, &z0hst, &z1hst, &idhst,
					   &x0dev, &x1dev, &y0dev, &y1dev, &z0dev, &z1dev, &iddev, &samplePos0, &samplePos1, sample);
  float *xhst, *yhst, *zhst;  particlePos particlePos_hst;
  float *xdev, *ydev, *zdev;  particlePos particlePos_dev;
  int *rank_hst, *rank_dev, *idx_dev;
  domainDecomposeKey domDecKey;
  const muse alloc_pos = allocateParticlePosition(&xhst, &yhst, &zhst, &particlePos_hst, &xdev, &ydev, &zdev, &particlePos_dev, &rank_hst, &rank_dev, &idx_dev, &domDecKey, Ntot);
#else///EXCHANGE_USING_GPUS
  domainDecomposeKey *domDecKey;
  domDecKey = (domainDecomposeKey *)malloc(num_max * sizeof(domainDecomposeKey));
  const muse alloc_dd = {num_max * sizeof(domainDecomposeKey), 0};
  if( domDecKey == NULL ){    __KILL__(stderr, "ERROR: failure to allocate domDecKey\n");  }
  int ndim;
  MPIinfo   *orbCfg;
  domainCfg *domCfg;
  real samplingRate;
  int sampleNumMax;
  real *sampleLoc, *sampleFul;
  real *domainMin, *domainMax;
  configORBtopology(&domCfg, &iparticleSendBuf, &iparticleRecvBuf, &letcfg, &ndim, &orbCfg,
		    Ntot, &samplingRate, &sampleNumMax, &sampleLoc, &sampleFul,
		    &sampleRecvNum, &sampleRecvDsp, &domainMin, &domainMax);
#endif//EXCHANGE_USING_GPUS
#ifdef  BUILD_LET_ON_DEVICE
  soaGEO soaGEO_dev;
  real *r2geo_dev;
/* #ifdef  CUB_AVAILABLE */
/*   void *r2geo_temp_storage; */
/* #endif//CUB_AVAILABLE */
  const muse alloc_r2geo = allocGeometricEnclosingBall_dev(&r2geo_dev,
/* #ifdef  CUB_AVAILABLE */
/* 							   &r2geo_temp_storage, */
/* #endif//CUB_AVAILABLE */
							   &soaGEO_dev, num_max);
  int *numSend_hst, *numSend_dev;
#endif//BUILD_LET_ON_DEVICE
  static domainInfo *nodeInfo;
  static position   *ipos;
#ifdef  GADGET_MAC
  static real *amin;
#endif//GADGET_MAC
  int Nstream_let;
  cudaStream_t *stream_let;
#ifdef  ALLOCATE_LETBUFFER
  uint *buf4LET;
  soaTreeWalkBuf soaWalk_dev;
#endif//ALLOCATE_LETBUFFER
  const muse alloc_LETtopology = configLETtopology(&nodeInfo, &ipos,
#ifdef  GADGET_MAC
						   &amin,
#endif//GADGET_MAC
#ifdef  BUILD_LET_ON_DEVICE
						   &numSend_hst, &numSend_dev,
#endif//BUILD_LET_ON_DEVICE
#ifdef  ALLOCATE_LETBUFFER
						   &buf4LET, &soaWalk_dev,
#endif//ALLOCATE_LETBUFFER
/* #ifdef  USE_CUDA_EVENT */
/* 						   &iniMakeLET, &finMakeLET, */
/* #endif//USE_CUDA_EVENT */
						   &stream_let, &Nstream_let, devProp, letcfg);
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
#ifdef COMPARE_WITH_DIRECT_SOLVER
  static char accfile[16];
#ifdef  GADGET_MAC
  sprintf(accfile, "%s.%.2d", APPG2, -(int)FLOOR(LOG2(absErr)));
#else///GADGET_MAC
#ifdef  WS93_MAC
  sprintf(accfile, "%s.%.2d", APPWS, -(int)FLOOR(LOG2(accErr)));
#else///WS93_MAC
  sprintf(accfile, "%s.%.2d", APPBH,  (int)FLOOR( theta * TEN));
#endif//WS93_MAC
#endif//GADGET_MAC
  acceleration *direct, *direct_dev;
#ifdef  GADGET_MAC
  acceleration *direct_old_dev;
#endif//GADGET_MAC
  const muse alloc_acc_dev = allocAccel_dev(num, &direct_dev, &direct
#ifdef  GADGET_MAC
					    , &direct_old_dev
#endif//GADGET_MAC
					    );
  iparticle ibody_direct_dev;
#      ifndef GENERATE_PHKEY_ON_DEVICE
  ibody_direct_dev     = ibody0_dev;
  ibody_direct_dev.acc = direct_dev;
#ifdef  GADGET_MAC
  ibody_direct_dev.acc_old = direct_old_dev;
#endif//GADGET_MAC
#      endif//GENERATE_PHKEY_ON_DEVICE
#endif//COMPARE_WITH_DIRECT_SOLVER
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* initialize the simulation run */
  //-----------------------------------------------------------------------
  /* variables for automatic optimization */
  static struct timeval start;
  static measuredTime elapsed;
#ifndef SERIALIZED_EXECUTION
  elapsed.sum_excg = elapsed.sum_rebuild = elapsed.excg = 0.0;
#endif//SERIALIZED_EXECUTION
#ifdef  WALK_TREE_COMBINED_MODEL
  static autoTuningParam rebuildParam;
  initStatVal(&rebuildParam.linearStats);  initGuessTime(&rebuildParam.linearGuess);
  initStatVal(&rebuildParam. powerStats);  initGuessTime(&rebuildParam. powerGuess);
#       ifdef  USE_PARABOLIC_GROWTH_MODEL
  initStatVal  (&rebuildParam.parabolicStats);
  initGuessTime(&rebuildParam.parabolicGuess);
#ifdef  USE_ADDITIONAL_SWITCH
  static int useParabolicGuess = 0;
#endif//USE_ADDITIONAL_SWITCH
#       endif//USE_PARABOLIC_GROWTH_MODEL
#endif//WALK_TREE_COMBINED_MODEL
  //-----------------------------------------------------------------------
  /* variables to control tree rebuild interval */
  static rebuildTree rebuild;
  rebuild.reuse = 1;
  rebuild.interval = 0.0;
#ifdef  BLOCK_TIME_STEP
  rebuild.adjust = false;
#endif//BLOCK_TIME_STEP
#   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
  rebuild.avg = 0.0;
  rebuild.var = 0.0;
#endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
  //-----------------------------------------------------------------------
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
  static brentStatus brentDistance;  brentDistance.initialized = false;
  static brentMemory brentHistory = {0.0, 0.0, 0, 0};
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
  //-----------------------------------------------------------------------
  /* read initial condition */
#   if  defined(MONITOR_ENERGY_ERROR) && defined(USE_HDF5_FORMAT)
  static energyError relEneErr;
#endif//defined(MONITOR_ENERGY_ERROR) && defined(USE_HDF5_FORMAT)
  /* static real time, dt; */
  static double time, dt;
  static ulong steps = 0;
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
#endif//USE_HDF5_FORMAT
#ifdef  SERIALIZED_EXECUTION
  readTentativeData        (&time, &dt, &steps, num, ibody0, file, last
#ifdef  USE_HDF5_FORMAT
			    , hdf5type
			    , &dropPrevTune, &rebuild, &elapsed
#ifdef  WALK_TREE_COMBINED_MODEL
			    , &rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
			    , &brentDistance, &brentHistory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
			    , &relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
			    );
#else///SERIALIZED_EXECUTION
  readTentativeDataParallel(&time, &dt, &steps, &num, ibody0, file, last, &iocfg, Ntot
#ifdef  USE_HDF5_FORMAT
			    , hdf5type
			    , &dropPrevTune, &rebuild, &elapsed
#ifdef  WALK_TREE_COMBINED_MODEL
			    , &rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
			    , &brentDistance, &brentHistory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
			    , &relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
			    );
#endif//SERIALIZED_EXECUTION
#ifndef  BLOCK_TIME_STEP
  real *dt_dev;
  const muse alloc_dt_dev = allocTimeStep_dev(&dt_dev);
#endif//BLOCK_TIME_STEP
#ifdef  REPORT_TOTAL_ELAPSED_TIME
  static ulong stepsInit;
  stepsInit = steps;
#endif//REPORT_TOTAL_ELAPSED_TIME
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < num; ii++)
    if( ibody0.pos[ii].x > 2.499790e-1f )
      ibody0.pos[ii].m = ZERO;
#endif
  //-----------------------------------------------------------------------
  /* set up for output files */
  static double formerTime;
#ifdef  SERIALIZED_EXECUTION
  formerTime = getPresentTimeInMin();
#else///SERIALIZED_EXECUTION
  if( mpi.rank == 0 )
    formerTime = getPresentTimeInMin();
  chkMPIerr(MPI_Bcast(&formerTime, 1, MPI_DOUBLE, 0, mpi.comm));
#endif//SERIALIZED_EXECUTION
#   if  defined(BLOCK_TIME_STEP) || !defined(EXEC_BENCHMARK)
  /* static real invSnapshotInterval; */
  /* invSnapshotInterval = UNITY / snapshotInterval; */
  static double invSnapshotInterval;
  invSnapshotInterval = 1.0 / snapshotInterval;
  static uint previous;
  previous = (uint)(time * invSnapshotInterval);
#endif//defined(BLOCK_TIME_STEP) || !defined(EXEC_BENCHMARK)
  //-----------------------------------------------------------------------
  /* check size of memory in use */
  muse       body_mem = {alloc_phkey.host   + alloc_ibody0.host   + alloc_ibody1.host   + alloc_ibody_dev.host,
			 alloc_phkey.device + alloc_ibody0.device + alloc_ibody1.device + alloc_ibody_dev.device};
  const muse tree_mem = {alloc_jtag.host   + alloc_jtag_dev.host   + alloc_cell.host   + alloc_cell_dev.host   + alloc_node_dev.host,
			 alloc_jtag.device + alloc_jtag_dev.device + alloc_cell.device + alloc_cell_dev.device + alloc_node_dev.device};
  muse       misc_mem = {alloc_lane_dev.host, alloc_lane_dev.device};
#ifdef  USE_HDF5_FORMAT
  body_mem.host   += alloc_snap.host;
  body_mem.device += alloc_snap.device;
#endif//USE_HDF5_FORMAT
  static muse used_mem;
  used_mem.host   = body_mem.host   + tree_mem.host   + misc_mem.host  ;
  used_mem.device = body_mem.device + tree_mem.device + misc_mem.device;
#ifndef SERIALIZED_EXECUTION
  const muse  let_mem = {alloc_dd.host	 + alloc_LETtopology.host  ,
			 alloc_dd.device + alloc_LETtopology.device};
  used_mem.host   += let_mem.host;
  used_mem.device += let_mem.device;
#ifdef  EXCHANGE_USING_GPUS
  const muse box_mem = {alloc_spl.host   + alloc_pos.host  ,
			alloc_spl.device + alloc_pos.device};
  used_mem.host   += box_mem.host;
  used_mem.device += box_mem.device;
#endif//EXCHANGE_USING_GPUS
#ifdef  BUILD_LET_ON_DEVICE
  const muse ball_mem = {alloc_r2geo.host,
			 alloc_r2geo.device};
  used_mem.host   += ball_mem.host;
  used_mem.device += ball_mem.device;
#endif//BUILD_LET_ON_DEVICE
#endif//SERIALIZED_EXECUTION
#ifdef  COMPARE_WITH_DIRECT_SOLVER
  used_mem.host   += alloc_acc_dev.host;
  used_mem.device += alloc_acc_dev.device;
#endif//COMPARE_WITH_DIRECT_SOLVER
#ifndef BLOCK_TIME_STEP
  used_mem.host   += alloc_dt_dev.host  ;
  used_mem.device += alloc_dt_dev.device;
#endif//BLOCK_TIME_STEP
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  misc_mem.host   += alloc_neighbor_dev.host;
  misc_mem.device += alloc_neighbor_dev.device;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  //-----------------------------------------------------------------------
  /* declarations of arrays to prevent for stack overflow */
  int *fail_dev;
  uint *buffer, *freeLst;
#ifndef ALLOCATE_LETBUFFER
  soaTreeWalkBuf soaWalk_dev;
#endif//ALLOCATE_LETBUFFER
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
  uint *freeNum;
  int *active;
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#ifndef USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  unsigned long long int *cycles_hst, *cycles_dev;
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  unsigned long long int *cycles_let_hst, *cycles_let_dev;
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#endif//USE_CUDA_EVENT
  const muse alloc_buf_dev = allocTreeBuffer_dev
    (&fail_dev, &buffer, &freeLst,
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
     &freeNum, &active,
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#ifndef USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
     &cycles_hst, &cycles_dev,
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
     &cycles_let_hst, &cycles_let_dev,
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#endif//USE_CUDA_EVENT
     &soaWalk_dev, num_max, used_mem, devProp);
  used_mem.host   += alloc_buf_dev.host;
  used_mem.device += alloc_buf_dev.device;
  //-----------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && !defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) && !defined(TIME_BASED_MODIFICATION)
  soaMakeBuf.ubuf_external = buffer;
  soaMakeBuf.Nbuf_external = soaWalk_dev.bufSize * (size_t)(NGROUPS * NBLOCKS_PER_SM * devProp.numSM);
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && !defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) && !defined(TIME_BASED_MODIFICATION)
  //-----------------------------------------------------------------------
  /* output memory usage */
#ifndef SERIALIZED_EXECUTION
  if( mpi.rank == 0 )
#endif//SERIALIZED_EXECUTION
    {
      //-------------------------------------------------------------------
      static FILE *fp;
      static char memlog[256];
      sprintf(memlog, "%s/%s.%s.txt", LOGFOLDER, file, "memory");
      fp = fopen(memlog, "w");
      if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", memlog);  }
      //-------------------------------------------------------------------
      size_t free, used, total;
      queryFreeDeviceMemory(&free, &total);      used = total - free;
      //-------------------------------------------------------------------
      fprintf(fp, "Allocated memory on device per process: %zu B (%zu GiB, %zu MiB)\n", used_mem.device, used_mem.device >> 30, used_mem.device >> 20);
      fprintf(fp, "    total memory on device per process: %zu B (%zu GiB, %zu MiB)\n"       , total, total >> 30, total >> 20);
      fprintf(fp, "     used memory on device per process: %zu B (%zu GiB, %zu MiB; %lf%%)\n",  used,  used >> 30,  used >> 20, 100.0 * (double)used / (double)total);
      fprintf(fp, "     free memory on device per process: %zu B (%zu GiB, %zu MiB; %lf%%)\n",  free,  free >> 30,  free >> 20, 100.0 * (double)free / (double)total);
      //-------------------------------------------------------------------
      fprintf(fp, "Allocated memory on   host per process: %zu B (%zu GiB, %zu MiB)\n", used_mem.host  , used_mem.host   >> 30, used_mem.host   >> 20);
      //-------------------------------------------------------------------
      fprintf(fp, "\nBreakdown of allocated memory on device:\n");
      fprintf(fp, "\t%zu MiB (%zu GiB) for N-body particles\n"  , body_mem.device >> 20, body_mem.device >> 30);
      fprintf(fp, "\t%zu MiB (%zu GiB) for tree structure\n"    , tree_mem.device >> 20, tree_mem.device >> 30);
      fprintf(fp, "\t%zu MiB for miscellaneous data\n", misc_mem.device >> 20);
      fprintf(fp, "\t%zu MiB (%zu GiB) for tree walk buffer\n"  , alloc_buf_dev.device >> 20, alloc_buf_dev.device >> 30);
#ifndef SERIALIZED_EXECUTION
      fprintf(fp, "\t%zu MiB (%zu B) for LET related data\n"  ,  let_mem.device >> 20,  let_mem.device);
#ifdef  EXCHANGE_USING_GPUS
      fprintf(fp, "\t%zu MiB (%zu KiB) for domain decomposition related data\n", box_mem.device >> 20, box_mem.device >> 10);
#endif//EXCHANGE_USING_GPUS
#ifdef  BUILD_LET_ON_DEVICE
      fprintf(fp, "\t%zu MiB (%zu KiB) for enclosing ball related data\n", ball_mem.device >> 20, ball_mem.device >> 10);
#endif//BUILD_LET_ON_DEVICE
#endif//SERIALIZED_EXECUTION
#ifdef  COMPARE_WITH_DIRECT_SOLVER
      fprintf(fp, "\t%zu MiB (%zu GiB) for accuracy test data\n", alloc_acc_dev.device >> 20, alloc_acc_dev.device >> 30);
#endif//COMPARE_WITH_DIRECT_SOLVER
#ifndef BLOCK_TIME_STEP
      fprintf(fp, "\t%zu B (%zu KiB) for time step value\n", alloc_dt_dev.device >> 20, alloc_dt_dev.device >> 10);
#endif//BLOCK_TIME_STEP
      //-------------------------------------------------------------------
      fprintf(fp, "\nBreakdown of allocated memory on host:\n");
      fprintf(fp, "\t%zu MiB (%zu GiB) for N-body particles\n"  , body_mem.host >> 20, body_mem.host >> 30);
      fprintf(fp, "\t%zu MiB (%zu GiB) for tree structure\n"    , tree_mem.host >> 20, tree_mem.host >> 30);
      fprintf(fp, "\t%zu MiB for miscellaneous data\n", misc_mem.host >> 20);
      fprintf(fp, "\t%zu MiB (%zu B) for tree walk buffer\n", alloc_buf_dev.host >> 20, alloc_buf_dev.host);
#ifndef SERIALIZED_EXECUTION
      fprintf(fp, "\t%zu MiB (%zu B) for LET related data\n", let_mem.host >> 20, let_mem.host);
#ifdef  EXCHANGE_USING_GPUS
      fprintf(fp, "\t%zu MiB (%zu KiB) for domain decomposition related data\n", box_mem.host >> 20, box_mem.host >> 10);
#endif//EXCHANGE_USING_GPUS
#ifdef  BUILD_LET_ON_DEVICE
      fprintf(fp, "\t%zu MiB (%zu B) for enclosing ball related data\n", ball_mem.host >> 20, ball_mem.host);
#endif//BUILD_LET_ON_DEVICE
#endif//SERIALIZED_EXECUTION
#ifdef  COMPARE_WITH_DIRECT_SOLVER
      fprintf(fp, "\t%zu MiB (%zu GiB) for accuracy test data\n", alloc_acc_dev.host >> 20, alloc_acc_dev.host >> 30);
#endif//COMPARE_WITH_DIRECT_SOLVER
#ifndef BLOCK_TIME_STEP
      fprintf(fp, "\t%zu B (%zu KiB) for time step value\n", alloc_dt_dev.host >> 20, alloc_dt_dev.host >> 10);
#endif//BLOCK_TIME_STEP
      //-------------------------------------------------------------------
      fclose(fp);
      //-------------------------------------------------------------------
    }/* if( mpi.rank == 0 ) */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* preparation of the benchmark */
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  /* declaration of counters */
  static wall_clock_time execTime[BENCHMARK_STEPS];
#ifdef  COUNT_INTERACTIONS
  static tree_metrics    treeProp[BENCHMARK_STEPS];
#endif//COUNT_INTERACTIONS
  //-----------------------------------------------------------------------
  {
    //---------------------------------------------------------------------
    /* zero-clear of all counters */
    //---------------------------------------------------------------------
    const wall_clock_time zero_time =
      {0.0, 0.0, 0.0,/* calcGravity_dev, calcMultipole, makeTree */
#ifdef  BLOCK_TIME_STEP
       0.0, 0.0, 0.0, 0.0,/* prediction_dev, correction_dev, setLaneTime_dev, adjustParticleTime_dev */
#else///BLOCK_TIME_STEP
       0.0, 0.0,/* advPos_dev, advVel_dev */
#endif//BLOCK_TIME_STEP
       0.0, 0.0,/* setTimeStep_dev, sortParticlesPHcurve */
       0.0, 0.0,/* copyParticle_dev2hst, copyParticle_hst2dev */
       0.0, 0.0/* setTreeNode_dev, setTreeCell_dev */
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
       , 0.0, 0.0/* examineNeighbor_dev, searchNeighbor_dev */
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  HUNT_MAKE_PARAMETER
       , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/* genPHkey_kernel, rsortKey_library, sortBody_kernel, makeTree_kernel, linkTree_kernel, trimTree_kernel */
       , 0.0, 0.0, 0.0, 0.0, 0.0/* initTreeLink_kernel, initTreeCell_kernel, initTreeNode_kernel, initTreeBody_kernel, copyRealBody_kernel */
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
       , 0.0, 0.0, 0.0, 0.0/* searchNeighbor_kernel, sortNeighbors, countNeighbors_kernel, commitNeighbors */
#endif//HUNT_FIND_PARAMETER
      };
#ifdef  COUNT_INTERACTIONS
    const tree_stats tree_zero = {ZERO, ZERO, ZERO, ZERO, 0, 0};
    tree_metrics zero_tree;
    for(uint ll = 0; ll < MAXIMUM_PHKEY_LEVEL; ll++)
      zero_tree.level[ll] = tree_zero;
    zero_tree.total = tree_zero;
    zero_tree.Ninteractions = 0;
#endif//COUNT_INTERACTIONS
    //---------------------------------------------------------------------
    for(uint ll = 0; ll < BENCHMARK_STEPS; ll++){
      execTime[ll] = zero_time;
#ifdef  COUNT_INTERACTIONS
      treeProp[ll] = zero_tree;
#endif//COUNT_INTERACTIONS
    }/* for(uint ll = 0; ll < BENCHMARK_STEPS; ll++){ */
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  static uint bench_begin;
  bench_begin = steps;
  //-----------------------------------------------------------------------
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  /* set N-body particles on device */
  copyParticle_hst2dev(num, ibody0, ibody0_dev
#ifdef  EXEC_BENCHMARK
		       , &execTime[0]
#endif//EXEC_BENCHMARK
		       );
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
#ifdef  EXCHANGE_USING_GPUS
  static loadImbalanceDetector balancer;
  balancer.enable  = false;
  balancer.execute = false;
  static double exchangeInterval = 0.0;/* corresponding to rebuildTree.interval */
  static autoTuningParam exchangeParam;
  initStatVal(&exchangeParam.linearStats);  initGuessTime(&exchangeParam.linearGuess);
  initStatVal(&exchangeParam. powerStats);  initGuessTime(&exchangeParam. powerGuess);
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  initStatVal  (&exchangeParam.parabolicStats);
  initGuessTime(&exchangeParam.parabolicGuess);
#endif//USE_PARABOLIC_GROWTH_MODEL
#endif//EXCHANGE_USING_GPUS
  //-----------------------------------------------------------------------
  updateDomain(&num, num_max, &Ni,
#ifdef  EXCHANGE_USING_GPUS
	       &ibody0_dev, &ibody1_dev, &ibody0, &ibody1, Ntot,
	       particlePos_hst, particlePos_dev, domCfg, domDecKey,
	       sample, samplePos0, samplePos1, soaPH_dev, devProp, devInfo,
	       &exchangeParam, &exchangeInterval, ormCfg, repCfg,
#else///EXCHANGE_USING_GPUS
	       ibody0_dev, ibody0, ibody1,
	       ndim, orbCfg, domCfg, domDecKey,
 samplingRate, sampleNumMax, sampleLoc, sampleFul, sampleRecvNum, sampleRecvDsp, domainMin, domainMax,
#endif//EXCHANGE_USING_GPUS
#ifdef  CARE_EXTERNAL_PARTICLES
	       &location,
#endif//CARE_EXTERNAL_PARTICLES
	       letcfg, iparticleSendBuf, iparticleRecvBuf, &elapsed
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
	       , &brentDistance, &brentHistory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  EXEC_BENCHMARK
	       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
      );
  //-----------------------------------------------------------------------
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* start the numerical simulation */
  //-----------------------------------------------------------------------
#ifndef MAKE_TREE_ON_DEVICE
  int remCell = 0;
#endif//MAKE_TREE_ON_DEVICE
  int numCell = NUM_ALLOC_TREE_CELL;
  int numNode = NUM_ALLOC_TREE_NODE;
  //-----------------------------------------------------------------------
  int bottomLev;
  buildTreeStructure(num,
#   if  !defined(MAKE_TREE_ON_DEVICE) || !defined(GENERATE_PHKEY_ON_DEVICE)
		     soaPH_hst,
#endif//!defined(MAKE_TREE_ON_DEVICE) || !defined(GENERATE_PHKEY_ON_DEVICE)
#ifdef  GENERATE_PHKEY_ON_DEVICE
		     &ibody0_dev, &ibody1_dev, soaPH_dev, devProp,
#ifdef  CUB_AVAILABLE
		     soaPH_pre,
#endif//CUB_AVAILABLE
#else///GENERATE_PHKEY_ON_DEVICE
		     &ibody0, &ibody1, ibody0_dev, tag,
#endif//GENERATE_PHKEY_ON_DEVICE
		     &bottomLev, &numCell, &numNode,
#ifdef  MAKE_TREE_ON_DEVICE
		     bottomLev_dev, scanNum_dev, numCell_dev, numNode_dev, soaMakeBuf
#else///MAKE_TREE_ON_DEVICE
		     &remCell, soaCell_hst, soaNode_hst
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
		     , soaCell_dev, soaNode_dev
#endif//defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
		     , &location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
		     , &start
#ifdef  EXEC_BENCHMARK
		     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		     );
#if 0
  printf("numNode = %d @ rank %d\n", numNode, letcfg.rank);
  MPI_Finalize();
  exit(0);
#endif
  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
  prediction_dev(num, time, ibody0_dev
#ifndef CALC_MULTIPOLE_ON_DEVICE
		 , ibody0
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef  EXEC_BENCHMARK
		 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		 );
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
#   if  !defined(CALC_MULTIPOLE_ON_DEVICE) || (!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE))
  real bmax_root;
#endif//!defined(CALC_MULTIPOLE_ON_DEVICE) || (!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE))
#if 0
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && defined(TIME_BASED_MODIFICATION)
  /* since dt is initialized as zero by MAGI, the below statement cause division by zero */
  /* since elapsed = 0.0f in here, wrong dtmin never cause problem */
  if( location.dtmin > dt ){
    location.dtmin = dt;
    location.dtinv = 1.0f / location.dtmin;
  }/* if( location.dtmin < dt ){ */
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && defined(TIME_BASED_MODIFICATION)
#endif
  setMultipoleMoment
    (NUM_PHKEY_LEVEL - 1, soaCell_dev, numNode, soaNode_dev, num
#ifdef  CALC_MULTIPOLE_ON_DEVICE
     , ibody0_dev, soaMakeBuf, devProp
#   if  !defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
     , soaNode_hst, &bmax_root
#endif//!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
#ifdef  COUNT_INTERACTIONS
     , soaCell_hst
#endif//COUNT_INTERACTIONS
#else///CALC_MULTIPOLE_ON_DEVICE
     , ibody, soaMakeBuf_hst, soaNode_hst, &bmax_root, soaCell_hst
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
     , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifndef SERIALIZED_EXECUTION
#ifdef  CARE_EXTERNAL_PARTICLES
     , &location
#endif//CARE_EXTERNAL_PARTICLES
     , &elapsed
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
     , treeProp[steps - bench_begin].level
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
#if 0
  printParticleInfo_dev(num, ibody0_dev, mpi);
#endif
#if 0
  printJtag_dev(num, soaNode_dev, mpi);
#endif
#if 0
  printRealBody_dev(num, soaNode_dev, mpi);
#endif
#if 0
  printTreeNode_dev(numNode, soaNode_dev, mpi);
#endif
#   if  !defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(SERIALIZED_EXECUTION)
  bmax_root = bmax[0];
#endif//!defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(SERIALIZED_EXECUTION)
                                                                           
/* #ifndef SERIALIZED_EXECUTION */
/*   elapsedCalc += ; */
/* #endif//SERIALIZED_EXECUTION */
                                                                           
  //-----------------------------------------------------------------------
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
  if( dropPrevTune == 1 ){
    /* first call of examineParticleSeparation() */
    /* therefore, brentDistance.a, brentDistance.b, brentDistance.x.pos are set by calling brentInit1st() */
    examineParticleSeparation
      (num, ibody0_dev
#ifdef  USE_BRENT_METHOD
       , &brentDistance
#else///USE_BRENT_METHOD
#ifdef  CUB_AVAILABLE
       , soaCUBneighbor
#endif//CUB_AVAILABLE
#endif//USE_BRENT_METHOD
#ifndef FACILE_NEIGHBOR_SEARCH
       , soaCell_dev, soaNode_dev, soaMakeBuf, soaNeighbor_dev, devProp
#endif//FACILE_NEIGHBOR_SEARCH
#ifdef  EXEC_BENCHMARK
       , execTime
#endif//EXEC_BENCHMARK
       );
    brentDistance.u = brentDistance.x;
 }/* if( dropPrevTune == 1 ){ */
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
  //-----------------------------------------------------------------------
  /* commit i-particle groups */
  int Ngrp;
  /* first call of configDistribution() */
  /* do not call brentCalc1st() */
  configDistribution
    (num, inumPerLane, &Ngrp, maxNgrp, laneInfo_hst, laneInfo_dev
#   if  defined(BLOCK_TIME_STEP) || (defined(LOCALIZE_I_PARTICLES) && defined(BRUTE_FORCE_LOCALIZATION))
     , ibody0_dev
#endif//defined(BLOCK_TIME_STEP) || (defined(LOCALIZE_I_PARTICLES) && defined(BRUTE_FORCE_LOCALIZATION))
#ifdef  LOCALIZE_I_PARTICLES
#ifdef  USE_BRENT_METHOD
     , &brentDistance
#else///USE_BRENT_METHOD
#ifdef  CUB_AVAILABLE
     , soaCUBneighbor
#endif//CUB_AVAILABLE
#endif//USE_BRENT_METHOD
#       ifdef  BRUTE_FORCE_LOCALIZATION
     , inum_dev, inum_hst
#ifndef FACILE_NEIGHBOR_SEARCH
     , soaCell_dev, soaNode_dev, soaMakeBuf, soaNeighbor_dev, devProp
#endif//FACILE_NEIGHBOR_SEARCH
#       else///BRUTE_FORCE_LOCALIZATION
     , peano
#       endif//BRUTE_FORCE_LOCALIZATION
#endif//LOCALIZE_I_PARTICLES
#ifdef  BLOCK_TIME_STEP
     , laneTime_dev
#       ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
     , ibody_dt_dev, ibody_dt_hst, &dtInfo_num, &dtInfo
#       endif//USE_VARIABLE_NEIGHBOR_LEVEL
#endif//BLOCK_TIME_STEP
     , start, &elapsed
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
  static double reduce = 1.0;
  //-----------------------------------------------------------------------
  /* preparation to construct LET */
#ifndef SERIALIZED_EXECUTION
#ifdef  BUILD_LET_ON_DEVICE
  calc_r2max_dev(Ngrp, laneInfo_dev, &ibody0_dev, soaGEO_dev
#ifdef  EXEC_BENCHMARK
		 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		 );
  shareNodePosition(letcfg.size, nodeInfo,    ipos, *(ibody0_dev.encBall_hst),
#ifdef  GADGET_MAC
		    amin, ibody0_dev.amin,
#endif//GADGET_MAC
		    letcfg);
  guessLETpartition(letcfg.size, nodeInfo, numNode, *(ibody0_dev.encBall_hst), letcfg);
#else///BUILD_LET_ON_DEVICE
  static position encBall;
  encBall.x = pj[0].x;
  encBall.y = pj[0].y;
  encBall.z = pj[0].z;
  encBall.m = bmax_root * bmax_root;
  shareNodePosition(letcfg.size, nodeInfo,    ipos, encBall,
#ifdef  GADGET_MAC
		    amin, ibody0_dev.amin,
#endif//GADGET_MAC
		    letcfg);
  guessLETpartition(letcfg.size, nodeInfo, numNode, encBall, letcfg);
#endif//BUILD_LET_ON_DEVICE
#if 0
  printf("r2max = %e, amin = %e for rank %d\n", ibody0_dev.encBall_hst->m, ibody0_dev.amin, letcfg.rank);
  MPI_Finalize();
  exit(0);
#endif
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  if( steps == 0 ){
    //---------------------------------------------------------------------
#ifdef  GADGET_MAC
    /* set Barnes-Hut MAC */
    enforceBarnesHutMAC_dev(Ni, ibody0_dev, numNode, soaNode_dev);
    /* force calculation based on opening criterion (theta = 1.0) */
    calcGravity_dev(
#ifdef  BLOCK_TIME_STEP
		    Ngrp, &reduce,
#endif//BLOCK_TIME_STEP
		    Ngrp, laneInfo_dev, Ni, ibody0_dev, soaNode_dev, soaWalk_dev, &sinfo, devProp, &elapsed.walkTree[0]
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
		    , file
#endif//PRINT_PSEUDO_PARTICLE_INFO
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  USE_CUDA_EVENT
		    , iniCalcAcc, finCalcAcc
#else///USE_CUDA_EVENT
		    , cycles_hst, cycles_dev
#endif//USE_CUDA_EVENT
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifndef SERIALIZED_EXECUTION
		    , &elapsed, numNode
#ifdef  LET_COMMUNICATION_VIA_HOST
		    , soaNode_hst
#endif//LET_COMMUNICATION_VIA_HOST
		    , letcfg.size, nodeInfo, Nstream_let, stream_let, letcfg
#ifdef  MONITOR_LETGEN_TIME
#ifdef  USE_CUDA_EVENT
		    , iniMakeLET, finMakeLET
#else///USE_CUDA_EVENT
		    , cycles_let_hst, cycles_let_dev
#endif//USE_CUDA_EVENT
#endif//MONITOR_LETGEN_TIME
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
		    , treeinfo_dev
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
		    , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
#ifdef  COMPARE_WITH_DIRECT_SOLVER
		    , true
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
		    , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#endif//COMPARE_WITH_DIRECT_SOLVER
		    );
    /* set GADGET-MAC */
    recoverGADGET_MAC_dev(numNode, soaNode_dev);
    /* reset stop watch */
    elapsed.walkTree[0] = 0.0;
#ifdef  EXEC_BENCHMARK
    execTime[steps - bench_begin].calcGravity_dev = 0.0;
#endif//EXEC_BENCHMARK
#ifndef SERIALIZED_EXECUTION
#ifdef  BUILD_LET_ON_DEVICE
    calc_r2max_dev(Ngrp, laneInfo_dev, &ibody0_dev, soaGEO_dev
#ifdef  EXEC_BENCHMARK
		   , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		   );
    shareNodePosition(letcfg.size, nodeInfo,    ipos, *(ibody0_dev.encBall_hst),
#ifdef  GADGET_MAC
		      amin, ibody0_dev.amin,
#endif//GADGET_MAC
		      letcfg);
#else///BUILD_LET_ON_DEVICE
    encBall.x = pj[0].x;
    encBall.y = pj[0].y;
    encBall.z = pj[0].z;
    encBall.m = bmax_root * bmax_root;
    shareNodePosition(letcfg.size, nodeInfo, ipos, encBall,
#ifdef  GADGET_MAC
		      amin, ibody0_dev.amin,
#endif//GADGET_MAC
		      letcfg);
#endif//BUILD_LET_ON_DEVICE
#endif//SERIALIZED_EXECUTION

#endif//GADGET_MAC
    //---------------------------------------------------------------------
    /* calculate gravitational acceleration and potential */
    calcGravity_dev(
#ifdef  BLOCK_TIME_STEP
		    Ngrp, &reduce,
#endif//BLOCK_TIME_STEP
		    Ngrp, laneInfo_dev, Ni, ibody0_dev, soaNode_dev, soaWalk_dev, &sinfo, devProp, &elapsed.walkTree[0]
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
		    , file
#endif//PRINT_PSEUDO_PARTICLE_INFO
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  USE_CUDA_EVENT
		    , iniCalcAcc, finCalcAcc
#else///USE_CUDA_EVENT
		    , cycles_hst, cycles_dev
#endif//USE_CUDA_EVENT
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifndef SERIALIZED_EXECUTION
		    , &elapsed, numNode
#ifdef  LET_COMMUNICATION_VIA_HOST
		    , soaNode_hst
#endif//LET_COMMUNICATION_VIA_HOST
		    , letcfg.size, nodeInfo, Nstream_let, stream_let, letcfg
#ifdef  MONITOR_LETGEN_TIME
#ifdef  USE_CUDA_EVENT
		    , iniMakeLET, finMakeLET
#else///USE_CUDA_EVENT
		    , cycles_let_hst, cycles_let_dev
#endif//USE_CUDA_EVENT
#endif//MONITOR_LETGEN_TIME
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
		    , treeinfo_dev
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
		    , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
#ifdef  COMPARE_WITH_DIRECT_SOLVER
		    , true
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
		    , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#endif//COMPARE_WITH_DIRECT_SOLVER
		    );
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
    brentDistance.u.val += elapsed.walkTree[0];
    brentHistory.totNum += Ngrp;
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
    rebuild.interval += 1.0;
#ifdef  WALK_TREE_COMBINED_MODEL
    linearModel(&(rebuildParam.linearGuess), &(rebuildParam.linearStats), rebuild.interval, elapsed.walkTree[0], 1.0);
    powerModel (&(rebuildParam. powerGuess), &(rebuildParam. powerStats), rebuild.interval, elapsed.walkTree[0], 1.0);
#       ifdef  USE_PARABOLIC_GROWTH_MODEL
    parabolicModel(&(rebuildParam.parabolicGuess), &(rebuildParam.parabolicStats), rebuild.interval, elapsed.walkTree[0], 1.0);
#ifdef  USE_ADDITIONAL_SWITCH
    useParabolicGuess |= (elapsed.walkTree[0] < elapsed.makeTree);
#endif//USE_ADDITIONAL_SWITCH
#       endif//USE_PARABOLIC_GROWTH_MODEL
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
    rebuild.avg += reduce;
    rebuild.var += reduce * reduce;
#endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
    //---------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
    copyCounters_dev2hst(Ni, treeinfo_dev, treeinfo);
    analyzeTreeMetrics(&treeProp[steps - bench_begin]);
    analyzeWalkStatistics(Ni, Ngrp, Nj  , &(treeProp[steps - bench_begin].Nj  ));
    analyzeWalkStatistics(Ni, Ngrp, Nbuf, &(treeProp[steps - bench_begin].Nbuf));
    for(int ii = 0; ii < Ni; ii++)
      treeProp[steps - bench_begin].Ninteractions += (ulong)Nj[ii];
    treeProp[steps - bench_begin].bufSize = soaWalk_dev.bufSize;
#endif//COUNT_INTERACTIONS
    //---------------------------------------------------------------------
#ifndef EXEC_BENCHMARK
#         if  defined(GENERATE_PHKEY_ON_DEVICE) && defined(COMPARE_WITH_DIRECT_SOLVER)
    ibody_direct_dev     = ibody0_dev;
    ibody_direct_dev.acc = direct_dev;
#ifdef  GADGET_MAC
    ibody_direct_dev.acc_old = direct_old_dev;
#endif//GADGET_MAC
#      endif//defined(GENERATE_PHKEY_ON_DEVICE) && defined(COMPARE_WITH_DIRECT_SOLVER)
    dumpSnapshot
      (unit, num, ibody0_dev, ibody0,
#ifdef  USE_HDF5_FORMAT
       &body_snapshot, hdf5type,
#ifdef  MONITOR_ENERGY_ERROR
       &relEneErr,
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
       time, steps, file, 0, &previous
#ifdef  LEAP_FROG_INTEGRATOR
       , ZERO
#endif//LEAP_FROG_INTEGRATOR
#ifndef SERIALIZED_EXECUTION
       , &iocfg, Ntot
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
#ifdef  COMPARE_WITH_DIRECT_SOLVER
       , Ni, ibody_direct_dev, devProp, direct, Ngrp, laneInfo_dev
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
       , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
       , accfile
#endif//COMPARE_WITH_DIRECT_SOLVER
       );
#endif//EXEC_BENCHMARK
    //---------------------------------------------------------------------
    /* initialize time step */
#ifndef BLOCK_TIME_STEP
    setTimeStep_dev(num, ibody0, eta, eps, dt_dev, &dt
#ifndef SERIALIZED_EXECUTION
		    , letcfg
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
		    , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		    );
#endif//BLOCK_TIME_STEP
#ifdef  LEAP_FROG_INTEGRATOR
    if( dt > eps * eta )
      dt = eps * eta;
#endif//LEAP_FROG_INTEGRATOR
    //---------------------------------------------------------------------
#if 0
    MPI_Finalize();
    /* fprintf(stderr, "# for normal run:  enable 2108-th line of tree/make_dev.cu\n"); */
    /* fprintf(stderr, "# for normal run: disable 1548-th line of main/gothic.c\n"); */
    /* fprintf(stderr, "# for normal run: disable  441-th line of sort/peano_dev.cu\n"); */
    /* fprintf(stderr, "# for normal run: disable  903-th line of sort/peano_dev.cu\n"); */
    /* fprintf(stderr, "# for normal run: disable 1344-th line of para/exchange_dev.cu\n"); */
    fprintf(stderr, "# for normal run: disable 2362-th line of main/gothic.c\n");
    fprintf(stderr, "# exit for debug @ %s(%d)\n", __FILE__, __LINE__);
    exit(0);
#endif
    //---------------------------------------------------------------------
  }/* if( steps == 0 ){ */
  //-----------------------------------------------------------------------
  /* initialize time step */
#ifdef  BLOCK_TIME_STEP
  adjustParticleTime_dev(Ngrp, laneInfo_dev, laneTime_dev, eps, eta, ibody0_dev
#ifdef  EXEC_BENCHMARK
			 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
			 );
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
#ifdef  LEAP_FROG_INTEGRATOR
  /* time integration for velocity to implement leap-frog method */
  advVel_dev(Ni, ibody0_dev, HALF * (real)dt
#ifdef  EXEC_BENCHMARK
	     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	     );
#endif//LEAP_FROG_INTEGRATOR
  //-----------------------------------------------------------------------
#if defined(DBG_TREE_WALK)
  __KILL__(stdout, "kill to dump results quickly\n");
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate time evolution */
  //-----------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
  bool existNewTree = false;
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
  while( time < ft ){
    //---------------------------------------------------------------------
    __NOTE__("t = %e(%zu step(s)), ft = %e\n", time, steps, ft);
    //---------------------------------------------------------------------
#if 0
    __FPRINTF__(stdout, "%zu-th step: excg = %e, sum_excg = %e, sum_rebuild = %e: interval = %d\n", steps, elapsed.excg, elapsed.sum_excg, elapsed.sum_rebuild, (int)nearbyint(exchangeInterval));
#endif
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
    //---------------------------------------------------------------------
    static double currentTime;
#ifdef  SERIALIZED_EXECUTION
    currentTime = getPresentTimeInMin();
#else///SERIALIZED_EXECUTION
    if( mpi.rank == 0 )
      currentTime = getPresentTimeInMin();
    chkMPIerr(MPI_Bcast(&currentTime, 1, MPI_DOUBLE, 0, mpi.comm));
#endif//SERIALIZED_EXECUTION
    if( (steps >= (bench_begin + BENCHMARK_STEPS - 1)) || (currentTime >= (formerTime + TFIN_MIN_BENCHMARK)) ){
      //-------------------------------------------------------------------
      /* output results of benchmark (wall clock time) */
      dumpBenchmark (jobID, file, 1 + steps - bench_begin, execTime);
      /* output results of benchmark (statistics of tree) */
#ifdef  COUNT_INTERACTIONS
      dumpStatistics(jobID, file, 1 + steps - bench_begin, MAXIMUM_PHKEY_LEVEL, treeProp);
#endif//COUNT_INTERACTIONS
      //-------------------------------------------------------------------
      /* exit loop about time evolution */
      break;
      //-------------------------------------------------------------------
    }/* if( (steps >= (bench_begin + BENCHMARK_STEPS - 1)) || (currentTime >= (formerTime + TFIN_MIN_BENCHMARK)) ){ */
    //---------------------------------------------------------------------
#endif//EXEC_BENCHMARK
    //---------------------------------------------------------------------
#   if  defined(BLOCK_TIME_STEP) || !defined(EXEC_BENCHMARK)
    static uint present;
#endif//defined(BLOCK_TIME_STEP) || !defined(EXEC_BENCHMARK)
    steps++;
    //---------------------------------------------------------------------
/* #   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP) */
/*     bool resetAllTimeStep = false; */
/* #endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP) */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* update octree structure if necessary */
    //---------------------------------------------------------------------
    /* copy position of N-body particles from device to host if necessary */
#ifndef CALC_MULTIPOLE_ON_DEVICE
    copyParticle_dev2hst(num, ibody0_dev, ibody0
#ifdef  EXEC_BENCHMARK
			 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
			 );
#endif//CALC_MULTIPOLE_ON_DEVICE
    //---------------------------------------------------------------------
    /* rebuild tree structure if required */
#ifndef SERIALIZED_EXECUTION
#ifdef  EXCHANGE_USING_GPUS
    /* detect slow down */
    if( !rebuild.reuse ){
      exchangeInterval += 1.0;
      linearModel(&(exchangeParam.linearGuess), &(exchangeParam.linearStats), exchangeInterval, elapsed.sum_rebuild, 1.0);
      powerModel (&(exchangeParam. powerGuess), &(exchangeParam. powerStats), exchangeInterval, elapsed.sum_rebuild, 1.0);
#ifdef  USE_PARABOLIC_GROWTH_MODEL
      parabolicModel(&(exchangeParam.parabolicGuess), &(exchangeParam.parabolicStats), exchangeInterval, elapsed.sum_rebuild, 1.0);
#endif//USE_PARABOLIC_GROWTH_MODEL
      balancer.execute = false;
      /* detect slow down of gravity calculation */
      if( balancer.enable ){
	double rchi2 = exchangeParam.linearGuess.rchisq;
	double guess = exchangeParam.linearGuess.time;
	if( rchi2 > exchangeParam.powerGuess.rchisq ){
	  rchi2 = exchangeParam.powerGuess.rchisq;
	  guess = exchangeParam.powerGuess.time;
	}/* if( rchi2 > exchangeParam.powerGuess.rchisq ){ */
#ifdef  USE_PARABOLIC_GROWTH_MODEL
	if( rchi2 > exchangeParam.parabolicGuess.rchisq ){
	  rchi2 = exchangeParam.parabolicGuess.rchisq;
	  guess = exchangeParam.parabolicGuess.time;
	}/* if( rchi2 > exchangeParam.parabolicGuess.rchisq ){ */
#endif//USE_PARABOLIC_GROWTH_MODEL
	if( guess > elapsed.excg )
	  balancer.execute = true;
#if 0
	__FPRINTF__(stdout, "%zu-th step: rchi2 = %e, guess = %e, exchange = %e\n", steps, rchi2, guess, elapsed.excg);
#endif
      }/* if( balancer.enable ){ */
    }/* if( !rebuild.reuse ){ */
    chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &(balancer.execute), 1, MPI_BOOL, MPI_LOR, letcfg.comm));

    /* detect load imbalance */
    if( balancer.enable & !balancer.execute ){
      balancer.tmin = balancer.tmax = elapsed.sum_excg;
      chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &(balancer.tmin), 1, MPI_DOUBLE, MPI_MIN, letcfg.comm));
      chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &(balancer.tmax), 1, MPI_DOUBLE, MPI_MAX, letcfg.comm));
#if 0
      if( mpi.rank == 0 )
	__FPRINTF__(stdout, "%zu-th step: tmin / tmax = %e, crit = %e\n", steps, balancer.tmin / balancer.tmax, loadImbalanceCrit);
#endif
      if( balancer.tmin < (loadImbalanceCrit * balancer.tmax) )
	balancer.execute = true;
    }/* if( balancer.enable & !balancer.execute ){ */
    __NOTE__("balancer.execute = %d @ rank %d\n", balancer.execute, letcfg.rank);
    //---------------------------------------------------------------------
    if( balancer.execute ){
      //-------------------------------------------------------------------
#if 1
      __FPRINTF__(stdout, "particle exchanging @ %zu-th step\n", steps);
#endif
      updateDomain(&num, num_max, &Ni,
		   &ibody0_dev, &ibody1_dev, &ibody0, &ibody1, Ntot,
		   particlePos_hst, particlePos_dev, domCfg, domDecKey,
		   sample, samplePos0, samplePos1, soaPH_dev, devProp, devInfo,
		   &exchangeParam, &exchangeInterval, ormCfg, repCfg,
#ifdef  CARE_EXTERNAL_PARTICLES
		   &location,
#endif//CARE_EXTERNAL_PARTICLES
		   letcfg, iparticleSendBuf, iparticleRecvBuf, &elapsed
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
		   , &brentDistance, &brentHistory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  EXEC_BENCHMARK
		   , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		   );
      //-------------------------------------------------------------------
      rebuild.reuse = 0;
      balancer.enable = false;
      //-------------------------------------------------------------------
    }/* if( balancer.execute ){ */
    //---------------------------------------------------------------------
#else///EXCHANGE_USING_GPUS
    //---------------------------------------------------------------------
    /* TENTATIVE IMPLEMENTATION: if the tree structure will be rebuild, then always exchange N-body particles among all MPI processes */
    chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &rebuild.reuse, 1, MPI_INT, MPI_LAND, letcfg.comm));
    if( !rebuild.reuse )
      updateDomain(&num, num_max, &Ni,
		   ibody0_dev, ibody0, ibody1,
		   ndim, orbCfg, domCfg, domDecKey,
		   samplingRate, sampleNumMax, sampleLoc, sampleFul, sampleRecvNum, sampleRecvDsp, domainMin, domainMax,
#ifdef  CARE_EXTERNAL_PARTICLES
		   &location,
#endif//CARE_EXTERNAL_PARTICLES
		   letcfg, iparticleSendBuf, iparticleRecvBuf, &elapsed
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
		   , &brentDistance, &brentHistory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  EXEC_BENCHMARK
		   , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		   );
#endif//EXCHANGE_USING_GPUS
#endif//SERIALIZED_EXECUTION
    if( !rebuild.reuse ){
      //-------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) && defined(EXCHANGE_USING_GPUS)
      elapsed.sum_rebuild = 0.0;
#endif//!defined(SERIALIZED_EXECUTION) && defined(EXCHANGE_USING_GPUS)
      //-------------------------------------------------------------------
#   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
      /* if Var(reduce) < (Avg(reduce))^2, set time of all particles same */
      /* if( ((rebuildInterval * reduceVar) / (reduceAvg * reduceAvg) - 1.0) < (1.0 - DBL_EPSILON) ) */
      if( DBL_EPSILON + rebuild.interval * rebuild.var < (2.0 * rebuild.avg * rebuild.avg) )
	rebuild.adjust = true;
	/* resetAllTimeStep = true; */
      rebuild.avg = 0.0;
      rebuild.var = 0.0;
#endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
      rebuild.interval = 0.0;
#ifdef  WALK_TREE_COMBINED_MODEL
      initStatVal(&(rebuildParam.   linearStats));      initGuessTime(&(rebuildParam.   linearGuess));
      initStatVal(&(rebuildParam.    powerStats));      initGuessTime(&(rebuildParam.    powerGuess));
#       ifdef  USE_PARABOLIC_GROWTH_MODEL
      initStatVal(&(rebuildParam.parabolicStats));      initGuessTime(&(rebuildParam.parabolicGuess));
#ifdef  USE_ADDITIONAL_SWITCH
      useParabolicGuess = 0;
#endif//USE_ADDITIONAL_SWITCH
#       endif//USE_PARABOLIC_GROWTH_MODEL
#endif//WALK_TREE_COMBINED_MODEL
#ifdef  WALK_TREE_TOTAL_SUM_MODEL
      elapsed.incSum = 0.0;
#endif//WALK_TREE_TOTAL_SUM_MODEL
#ifdef  EXEC_BENCHMARK
      printf("#rebuild @ %zu-th step\n", steps);
      fflush(stdout);
#endif//EXEC_BENCHMARK
      //-------------------------------------------------------------------
      __NOTE__("rebuild @ %zu-th step, t = %e, sum_rebuild = %e @ rank %d\n", steps, time, elapsed.sum_rebuild, mpi.rank);
      //-------------------------------------------------------------------
      /* copy position of N-body particles from device to host if necessary */
#ifndef SERIALIZED_EXECUTION
#   if  !defined(GENERATE_PHKEY_ON_DEVICE) && defined(CALC_MULTIPOLE_ON_DEVICE)
      copyParticle_dev2hst(num, ibody0_dev, ibody0
#ifdef  EXEC_BENCHMARK
			   , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
			   );
#endif//!defined(GENERATE_PHKEY_ON_DEVICE) && defined(CALC_MULTIPOLE_ON_DEVICE)
#endif//SERIALIZED_EXECUTION
      //-------------------------------------------------------------------
      buildTreeStructure(num,
#   if  !defined(MAKE_TREE_ON_DEVICE) || !defined(GENERATE_PHKEY_ON_DEVICE)
			 soaPH_hst,
#endif//!defined(MAKE_TREE_ON_DEVICE) || !defined(GENERATE_PHKEY_ON_DEVICE)
#ifdef  GENERATE_PHKEY_ON_DEVICE
			 &ibody0_dev, &ibody1_dev, soaPH_dev, devProp,
#ifdef  CUB_AVAILABLE
			 soaPH_pre,
#endif//CUB_AVAILABLE
#else///GENERATE_PHKEY_ON_DEVICE
			 &ibody0, &ibody1, ibody0_dev, tag,
#endif//GENERATE_PHKEY_ON_DEVICE
			 &bottomLev, &numCell, &numNode,
#ifdef  MAKE_TREE_ON_DEVICE
			 bottomLev_dev, scanNum_dev, numCell_dev, numNode_dev, soaMakeBuf
#else///MAKE_TREE_ON_DEVICE
			 &remCell, soaCell_hst, soaNode_hst
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
			 , soaCell_dev, soaNode_dev
#endif//defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
			 , &location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
			 , &start
#ifdef  EXEC_BENCHMARK
			 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
			 );
      //-------------------------------------------------------------------
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#ifdef  BLOCK_TIME_STEP
      copySortedParticles_dev(num, ibody0_dev);
#endif//BLOCK_TIME_STEP
      setMultipoleMoment
	(bottomLev, soaCell_dev, numNode, soaNode_dev, num
#ifdef  CALC_MULTIPOLE_ON_DEVICE
	 , ibody0_dev, soaMakeBuf, devProp
#   if  !defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
	 , soaNode_hst, &bmax_root
#endif//!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
#ifdef  COUNT_INTERACTIONS
	 , soaCell_hst
#endif//COUNT_INTERACTIONS
#else///CALC_MULTIPOLE_ON_DEVICE
	 , ibody, soaMakeBuf_hst, soaNode_hst, &bmax_root, soaCell_hst
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	 , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifndef SERIALIZED_EXECUTION
#ifdef  CARE_EXTERNAL_PARTICLES
	 , &location
#endif//CARE_EXTERNAL_PARTICLES
	 , &elapsed
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
	 , treeProp[steps - bench_begin].level
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
	 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	 );
#   if  !defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(SERIALIZED_EXECUTION)
      bmax_root = bmax[0];
#endif//!defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(SERIALIZED_EXECUTION)

#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
      //-------------------------------------------------------------------
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
      brentDistance.u.val /= (double)brentHistory.totNum;
      brentHistory.totNum = 0;
      brentHistory.interval++;
      static bool initialized = false;
      if( initialized ){
	//-----------------------------------------------------------------
	/* from the second time, managing brentStatus is necessary */
	//-----------------------------------------------------------------
	/* check # of sequences executed and the evolution of elapsed time */
	bool perturbBrent = false;
	if( brentHistory.interval > BRENT_METHOD_LAUNCH ){
	  //---------------------------------------------------------------
	  /* if increase of the elapsed time exceeds the allowable range, then perturb brentStatus */
	  if( brentDistance.u.val > BRENT_METHOD_ALLOW * brentHistory.previous )	    perturbBrent = true;
	  //---------------------------------------------------------------
	  /* when the elapsed time grows continuously, parturb the brentStatus */
	  if( brentDistance.u.val > brentHistory.previous ){
	    brentHistory.degraded++;
	    if( brentHistory.degraded > BRENT_METHOD_MODIFY )	      perturbBrent = true;
	  }/* if( brentDistance.u.val > brentHistory.previous ){ */
	  else	                                                    brentHistory.degraded = 0;
	  //---------------------------------------------------------------
	}/* if( brentHistory.interval > BRENT_METHOD_LAUNCH ){ */
	brentHistory.previous = brentDistance.u.val;
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	if( perturbBrent ){
	  //---------------------------------------------------------------
	  examineParticleSeparation
	    (num, ibody0_dev
#ifdef  USE_BRENT_METHOD
	     , &brentDistance
#else///USE_BRENT_METHOD
#ifdef  CUB_AVAILABLE
	     , soaCUBneighbor
#endif//CUB_AVAILABLE
#endif//USE_BRENT_METHOD
#ifndef FACILE_NEIGHBOR_SEARCH
	     , soaCell_dev, soaNode_dev, soaMakeBuf, soaNeighbor_dev, devProp
#endif//FACILE_NEIGHBOR_SEARCH
#ifdef  EXEC_BENCHMARK
	     , execTime
#endif//EXEC_BENCHMARK
	     );
	  //---------------------------------------------------------------
	  brentHistory.interval = 0;
	  brentHistory.degraded = 0;
	  //---------------------------------------------------------------
	}/* if( brentHistory.degraded > BRENT_METHOD_MODIFY ){ */
	else{
	  //---------------------------------------------------------------
	  /* when the brentStatus is conserved */
	  //---------------------------------------------------------------
	  brentCalc2nd(&brentDistance);
	  //---------------------------------------------------------------
	}/* else{ */
      }/* if( initialized ){ */
      else{
	//-----------------------------------------------------------------
	/* the first execution of tree rebuild */
	//-----------------------------------------------------------------
	brentDistance.x = brentDistance.u;
	brentInit2nd(&brentDistance);
	initialized = true;
	//-----------------------------------------------------------------
      }/* else{ */
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
      //-------------------------------------------------------------------
      /* from the second time, brentCalc1st() is called internally */
      configDistribution
	(num, inumPerLane, &Ngrp, maxNgrp, laneInfo_hst, laneInfo_dev
#   if  defined(BLOCK_TIME_STEP) || (defined(LOCALIZE_I_PARTICLES) && defined(BRUTE_FORCE_LOCALIZATION))
	 , ibody0_dev
#endif//defined(BLOCK_TIME_STEP) || (defined(LOCALIZE_I_PARTICLES) && defined(BRUTE_FORCE_LOCALIZATION))
#ifdef  LOCALIZE_I_PARTICLES
#ifdef  USE_BRENT_METHOD
	 , &brentDistance
#else///USE_BRENT_METHOD
#ifdef  CUB_AVAILABLE
	 , soaCUBneighbor
#endif//CUB_AVAILABLE
#endif//USE_BRENT_METHOD
#       ifdef  BRUTE_FORCE_LOCALIZATION
	 , inum_dev, inum_hst
#ifndef FACILE_NEIGHBOR_SEARCH
	 , soaCell_dev, soaNode_dev, soaMakeBuf, soaNeighbor_dev, devProp
#endif//FACILE_NEIGHBOR_SEARCH
#       else///BRUTE_FORCE_LOCALIZATION
	 , peano
#       endif//BRUTE_FORCE_LOCALIZATION
#endif//LOCALIZE_I_PARTICLES
#ifdef  BLOCK_TIME_STEP
	 , laneTime_dev
#       ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
	 , ibody_dt_dev, ibody_dt_hst, &dtInfo_num, &dtInfo
#       endif//USE_VARIABLE_NEIGHBOR_LEVEL
#endif//BLOCK_TIME_STEP
     , start, &elapsed
#ifdef  EXEC_BENCHMARK
	 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	 );
      //-------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
      existNewTree = true;
#ifdef  EXCHANGE_USING_GPUS
      if( !balancer.execute )      	balancer.enable  =  true;
      else	                        balancer.execute = false;
#endif//EXCHANGE_USING_GPUS
#endif//SERIALIZED_EXECUTION
      //-------------------------------------------------------------------
    }/* if( !rebuild.reuse ){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* share information on domain decomposition in the next time step */
    //---------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) && defined(EXCHANGE_USING_GPUS)
    chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &(balancer.enable), 1, MPI_BOOL, MPI_LOR, letcfg.comm));
    __NOTE__("balancer.enable = %d @ rank %d\n", balancer.enable, letcfg.rank);
    chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &existNewTree, 1, MPI_BOOL, MPI_LOR, letcfg.comm));
#endif//!defined(SERIALIZED_EXECUTION) && defined(EXCHANGE_USING_GPUS)
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* time integration */
    //---------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
    //---------------------------------------------------------------------
    setLaneTime_dev(Ngrp, laneInfo_dev, laneTime_dev, ibody0_dev
#ifdef  EXEC_BENCHMARK
		    , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		    );
    int grpNum = 0;
    setTimeStep_dev(Ngrp, laneInfo_dev, laneTime_dev, &grpNum, ibody0_dev,
		    time, &time, &dt, rebuild.adjust, invSnapshotInterval, previous, &present
/* #ifdef  CUB_AVAILABLE */
/* 		    , soaCUBtimeBuf */
/* #endif//CUB_AVAILABLE */
#ifndef SERIALIZED_EXECUTION
		    , letcfg
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
		    , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		    );
/* #   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP) */
/*     adjustAllTimeStep = resetAllTimeStep; */
/* #endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP) */
#ifdef  BLOCK_TIME_STEP
    rebuild.adjust = false;
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
#else///BLOCK_TIME_STEP
    //---------------------------------------------------------------------
#ifndef LEAP_FROG_INTEGRATOR
    /* initialize time step */
    setTimeStep_dev(Ni, ibody0, eta, eps, dt_dev, &dt
#ifndef SERIALIZED_EXECUTION
		    , letcfg
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
		    , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		    );
#endif//LEAP_FROG_INTEGRATOR
    //---------------------------------------------------------------------
    time += dt;
    //---------------------------------------------------------------------
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && defined(TIME_BASED_MODIFICATION)
    location.elapsed = (float)((double)location.elapsed + dt);
    if( location.dtmin > dt ){
      location.dtmin = dt;
      location.dtinv = 1.0f / location.dtmin;
    }/* if( location.dtmin < dt ){ */
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && defined(TIME_BASED_MODIFICATION)
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* orbit integration (predict step) */
    //---------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
    //---------------------------------------------------------------------
    prediction_dev(num, time, ibody0_dev
#ifndef CALC_MULTIPOLE_ON_DEVICE
		   , ibody0
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef  EXEC_BENCHMARK
		   , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		   );
    //---------------------------------------------------------------------
#else///BLOCK_TIME_STEP
    //---------------------------------------------------------------------
#ifndef LEAP_FROG_INTEGRATOR
    //---------------------------------------------------------------------
    /* predict velocity and update position of i-particles */
    advVel_dev(Ni, ibody0_dev, HALF * (real)dt
#ifdef  EXEC_BENCHMARK
	       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	       );
    advPos_dev(Ni, ibody0_dev,        (real)dt
#ifdef  EXEC_BENCHMARK
	       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	       );
    //---------------------------------------------------------------------
#else///LEAP_FROG_INTEGRATOR
    //---------------------------------------------------------------------
    /* update position of i-particles */
    advPos_dev(Ni, ibody0_dev, (real)dt
#ifdef  EXEC_BENCHMARK
	       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	       );
    //---------------------------------------------------------------------
#endif//LEAP_FROG_INTEGRATOR
    //---------------------------------------------------------------------
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* calculate multipole moment of tree nodes */
    //---------------------------------------------------------------------
    setMultipoleMoment
      (bottomLev, soaCell_dev, numNode, soaNode_dev, num
#ifdef  CALC_MULTIPOLE_ON_DEVICE
       , ibody0_dev, soaMakeBuf, devProp
#   if  !defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
       , soaNode_hst, &bmax_root
#endif//!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
#ifdef  COUNT_INTERACTIONS
       , soaCell_hst
#endif//COUNT_INTERACTIONS
#else///CALC_MULTIPOLE_ON_DEVICE
       , ibody, soaMakeBuf_hst, soaNode_hst, &bmax_root, soaCell_hst
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
       , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifndef SERIALIZED_EXECUTION
#ifdef  CARE_EXTERNAL_PARTICLES
       , &location
#endif//CARE_EXTERNAL_PARTICLES
       , &elapsed
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
       , treeProp[steps - bench_begin].level
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
#   if  !defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(SERIALIZED_EXECUTION)
    bmax_root = bmax[0];
#endif//!defined(CALC_MULTIPOLE_ON_DEVICE) && !defined(SERIALIZED_EXECUTION)
                                                                           
/* #ifndef SERIALIZED_EXECUTION */
/*     elapsedCalc += ; */
/* #endif//SERIALIZED_EXECUTION */
                                                                           
    //---------------------------------------------------------------------
    /* preparation to construct LET */
#ifndef SERIALIZED_EXECUTION
#ifdef  BUILD_LET_ON_DEVICE
    calc_r2max_dev(Ngrp, laneInfo_dev, &ibody0_dev, soaGEO_dev
#ifdef  EXEC_BENCHMARK
		   , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		   );
    shareNodePosition(letcfg.size, nodeInfo,    ipos, *(ibody0_dev.encBall_hst),
#ifdef  GADGET_MAC
		      amin, ibody0_dev.amin,
#endif//GADGET_MAC
		      letcfg);
    /* this function must be called when tree structure is rebuild */
    if( existNewTree )
      guessLETpartition(letcfg.size, nodeInfo, numNode, *(ibody0_dev.encBall_hst), letcfg);
#else///BUILD_LET_ON_DEVICE
    encBall.x = pj[0].x;
    encBall.y = pj[0].y;
    encBall.z = pj[0].z;
    encBall.m = bmax_root * bmax_root;
    shareNodePosition(letcfg.size, nodeInfo, ipos, encBall,
#ifdef  GADGET_MAC
		      amin, ibody0_dev.amin,
#endif//GADGET_MAC
		      letcfg);
    /* this function must be called when tree structure is rebuild */
    if( existNewTree )
      guessLETpartition(letcfg.size, nodeInfo, numNode, encBall, letcfg);
#endif//BUILD_LET_ON_DEVICE
    existNewTree = false;
#endif//SERIALIZED_EXECUTION
    //---------------------------------------------------------------------
    /* calculate acceleration of i-particles by j-cells */
#ifdef  SHOW_NI_DEPENDENCE
    printf("#file: %s\n", file);
    printf("#NTHREADS = %d, TSUB = %d, NWARP = %d\n", NTHREADS, TSUB, NWARP);
    printf("#grpNum\tNi_active\ttime(sec)\n");
    for(grpNum = Ngrp; grpNum >= NGROUPS; grpNum = (int)((double)grpNum * M_SQRT1_2)){
#endif//SHOW_NI_DEPENDENCE
      //---------------------------------------------------------------------
      calcGravity_dev(
#ifdef  BLOCK_TIME_STEP
		      grpNum, &reduce,
#endif//BLOCK_TIME_STEP
		      Ngrp, laneInfo_dev, Ni, ibody0_dev, soaNode_dev, soaWalk_dev, &sinfo, devProp, &elapsed.walkTree[rebuild.reuse]
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
		      , file
#endif//PRINT_PSEUDO_PARTICLE_INFO
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  USE_CUDA_EVENT
		      , iniCalcAcc, finCalcAcc
#else///USE_CUDA_EVENT
		      , cycles_hst, cycles_dev
#endif//USE_CUDA_EVENT
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifndef SERIALIZED_EXECUTION
		      , &elapsed, numNode
#ifdef  LET_COMMUNICATION_VIA_HOST
		      , soaNode_hst
#endif//LET_COMMUNICATION_VIA_HOST
		      , letcfg.size, nodeInfo, Nstream_let, stream_let, letcfg
#ifdef  MONITOR_LETGEN_TIME
#ifdef  USE_CUDA_EVENT
		      , iniMakeLET, finMakeLET
#else///USE_CUDA_EVENT
		      , cycles_let_hst, cycles_let_dev
#endif//USE_CUDA_EVENT
#endif//MONITOR_LETGEN_TIME
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
		      , treeinfo_dev
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
		      , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
#ifdef  COMPARE_WITH_DIRECT_SOLVER
		      , true
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
		      , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#endif//COMPARE_WITH_DIRECT_SOLVER
		      );
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
      brentDistance.u.val += elapsed.walkTree[rebuild.reuse];
#ifdef  BLOCK_TIME_STEP
      brentHistory.totNum += grpNum;
#else///BLOCK_TIME_STEP
      brentHistory.totNum += Ngrp;
#endif//BLOCK_TIME_STEP
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
/* #ifndef SERIALIZED_EXECUTION */
/*       elapsedCalc += elapsed.walkTree[reuse]; */
/* #endif//SERIALIZED_EXECUTION */
      //-------------------------------------------------------------------
#ifdef  SHOW_NI_DEPENDENCE
      //-------------------------------------------------------------------
      int Ni_active = 0;
      for(int ii = 0; ii < grpNum; ii++)
	Ni_active += laneInfo_hst[ii].num;
      printf("%d\t%d\t%le\n", grpNum, Ni_active, execTime[steps - bench_begin].calcGravity_dev);
      fflush(stdout);
      execTime[steps - bench_begin].calcGravity_dev = 0.0;
    }/* for(grpNum = Ngrp; grpNum >= NGROUPS; grpNum = (int)((double)grpNum * M_SQRT1_2)){ */
    steps = bench_begin + BENCHMARK_STEPS - 1;/* kill other benchmark */
    //---------------------------------------------------------------------
#endif//SHOW_NI_DEPENDENCE
    //---------------------------------------------------------------------
    rebuild.interval += 1.0;
#   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
    rebuild.avg += reduce;
    rebuild.var += reduce * reduce;
#endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
    //---------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
    linearModel(&(rebuildParam.linearGuess), &(rebuildParam.linearStats), rebuild.interval, elapsed.walkTree[rebuild.reuse] * reduce, reduce);
    powerModel (&(rebuildParam. powerGuess), &(rebuildParam. powerStats), rebuild.interval, elapsed.walkTree[rebuild.reuse] * reduce, reduce);
#       ifdef  USE_PARABOLIC_GROWTH_MODEL
    parabolicModel(&(rebuildParam.parabolicGuess), &(rebuildParam.parabolicStats), rebuild.interval, elapsed.walkTree[rebuild.reuse] * reduce, reduce);
#ifdef  USE_ADDITIONAL_SWITCH
     useParabolicGuess |= (elapsed.walkTree[rebuild.reuse] < elapsed.makeTree);
#endif//USE_ADDITIONAL_SWITCH
#       endif//USE_PARABOLIC_GROWTH_MODEL
    rebuild.reuse = 1;
#       ifdef  USE_PARABOLIC_GROWTH_MODEL
    if( rebuild.interval > 3.5 ){
      double rchi2 = rebuildParam.linearGuess.rchisq;
      double guess = rebuildParam.linearGuess.time;
      if( rchi2 >     rebuildParam.powerGuess.rchisq ){	rchi2 =     rebuildParam.powerGuess.rchisq;	guess =     rebuildParam.powerGuess.time;      }
#ifdef  USE_ADDITIONAL_SWITCH
      if( useParabolicGuess )
#endif//USE_ADDITIONAL_SWITCH
      if( rchi2 > rebuildParam.parabolicGuess.rchisq ){	rchi2 = rebuildParam.parabolicGuess.rchisq;	guess = rebuildParam.parabolicGuess.time;      }
      if( guess > (elapsed.makeTree * reduce) )
	rebuild.reuse = 0;
    }
#       else///USE_PARABOLIC_GROWTH_MODEL
    if( (rebuild.interval > 3.5) &&
	((rebuildParam.linearGuess.rchisq < 1.0e-30 + rebuildParam.powerGuess.rchisq) ? (rebuildParam.linearGuess.time) : (rebuildParam.powerGuess.time)) > (elapsed.makeTree * reduce) )
	rebuild.reuse = 0;
#       endif//USE_PARABOLIC_GROWTH_MODEL
#endif//WALK_TREE_COMBINED_MODEL
    //---------------------------------------------------------------------
#ifdef  WALK_TREE_TOTAL_SUM_MODEL
#ifdef  BLOCK_TIME_STEP
    if( !rebuild.reuse )
      elapsed.walkTree[0] /= reduce;
#endif//BLOCK_TIME_STEP
    const double diff = elapsed.walkTree[rebuild.reuse] - reduce * elapsed.walkTree[0];
    if( diff > 0.0 )
      elapsed.incSum += diff;
    reuse = ((elapsed.incSum + 0.5 * rebuildInterval * diff) < elapsed.makeTree) ? 1 : 0;
#endif//WALK_TREE_TOTAL_SUM_MODEL
    //---------------------------------------------------------------------
#ifdef  WALK_TREE_ARITHMETIC_PROGRESSION_MODEL
    /* based on model assuming arithmetic progression */
    /* if BLOCK_TIME_STEP is not defined, then reduce is unity */
#ifdef  BLOCK_TIME_STEP
    if( !rebuild.reuse )
      elapsed.walkTree[0] /= reduce;
#endif//BLOCK_TIME_STEP
    rebuild.reuse = 1;
    if( (rebuild.interval > 1.5) &&
    	( ((rebuild.interval * rebuild.interval) * (elapsed.walkTree[rebuild.reuse] - reduce * elapsed.walkTree[0])) > (2.0 * (rebuild.interval - 1.0) * elapsed.makeTree)) )
      rebuild.reuse = 0;
#endif//WALK_TREE_ARITHMETIC_PROGRESSION_MODEL
    //---------------------------------------------------------------------
#ifdef  WALK_TREE_GEOMETRIC_PROGRESSION_MODEL
    /* based on model assuming geometric progression */
    /* if BLOCK_TIME_STEP is not defined, then reduce is unity */
#ifdef  BLOCK_TIME_STEP
    if( !rebuild.reuse )
      elapsed.walkTree[0] /= reduce;
#endif//BLOCK_TIME_STEP
    rebuild.reuse = 1;
    if( rebuild.interval > 1.5 ){
      const double tinit = reduce * elapsed.walkTree[0];
      const double ratio = elapsed.walkTree[rebuild.reuse] / tinit;/* = ratio^(n-1) */
      if( ratio > DBL_EPSILON + 1.0 ){
	const double growth = pow(ratio, 1.0 / (rebuild.interval - 1.0));
	if( (rebuildInterval * rebuildInterval * elapsed.walkTree[rebuild.reuse]) > ((growth - 1.0) * elapsed.makeTree + (ratio * growth - 1.0) * tinit) )
	  rebuild.reuse = 0;
      }/* if( ratio > DBL_EPSILON + 1.0 ){ */
    }/* if( rebuild.interval > 1.5 ){ */
#endif//WALK_TREE_GEOMETRIC_PROGRESSION_MODEL
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* orbit integration (correct step) */
    //---------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
    //---------------------------------------------------------------------
    correction_dev(grpNum, laneInfo_dev, laneTime_dev, eps, eta, ibody0_dev, rebuild.reuse
#ifdef  EXEC_BENCHMARK
		   , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
		   );
    //---------------------------------------------------------------------
#else///BLOCK_TIME_STEP
    //---------------------------------------------------------------------
#ifndef LEAP_FROG_INTEGRATOR
    //---------------------------------------------------------------------
    /* correct velocity of i-particles */
    advVel_dev(Ni, ibody0_dev, HALF * (real)dt
#ifdef  EXEC_BENCHMARK
	       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	       );
    //---------------------------------------------------------------------
#else///LEAP_FROG_INTEGRATOR
    //---------------------------------------------------------------------
    /* update velocity of i-particles */
    advVel_dev(Ni, ibody0_dev, (real)dt
#ifdef  EXEC_BENCHMARK
	       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	       );
    //---------------------------------------------------------------------
#endif//LEAP_FROG_INTEGRATOR
    //---------------------------------------------------------------------
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
    copyCounters_dev2hst(Ni, treeinfo_dev, treeinfo);
    analyzeTreeMetrics(&treeProp[steps - bench_begin]);
    analyzeWalkStatistics(Ni, grpNum, Nj  , &(treeProp[steps - bench_begin].Nj  ));
    analyzeWalkStatistics(Ni, grpNum, Nbuf, &(treeProp[steps - bench_begin].Nbuf));
    for(int ii = 0; ii < Ni; ii++)
      treeProp[steps - bench_begin].Ninteractions += (ulong)Nj[ii];
    treeProp[steps - bench_begin].bufSize = soaWalk_dev.bufSize;
#endif//COUNT_INTERACTIONS
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
#ifndef EXEC_BENCHMARK
    //---------------------------------------------------------------------
#ifndef BLOCK_TIME_STEP
    /* output time evolution of numerical results */
    present = (uint)(time * invSnapshotInterval);
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
    /* save tentative results to restart the simulation */
    static double currentTime;
#ifdef  SERIALIZED_EXECUTION
    currentTime = getPresentTimeInMin();
#else///SERIALIZED_EXECUTION
    if( mpi.rank == 0 )
      currentTime = getPresentTimeInMin();
    chkMPIerr(MPI_Bcast(&currentTime, 1, MPI_DOUBLE, 0, mpi.comm));
#endif//SERIALIZED_EXECUTION
    if( currentTime > formerTime + saveInterval ){
      //-------------------------------------------------------------------
#ifndef COMPARE_WITH_DIRECT_SOLVER
      dumpRestarter
	(num, ibody0_dev, ibody0, time, dt, steps, file, &last, &formerTime
#ifdef  USE_HDF5_FORMAT
	 , hdf5type
	 , rebuild, elapsed
#ifdef  WALK_TREE_COMBINED_MODEL
	 , rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
	 , brentDistance, brentHistory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
	 , relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
#ifndef SERIALIZED_EXECUTION
	 , &iocfg, mpi, Ntot
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
	 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	 );
#endif//COMPARE_WITH_DIRECT_SOLVER
      //-------------------------------------------------------------------
      rebuild.reuse = 0;
#ifdef  BLOCK_TIME_STEP
      rebuild.adjust = false;
#endif//BLOCK_TIME_STEP
      //-------------------------------------------------------------------
    }/* if( currentTime > formerTime + saveInterval ){ */
    //---------------------------------------------------------------------
    /* output snapshot file to analyze and visualize the results */
    if( present != previous ){
      //-------------------------------------------------------------------
#   if  defined(GENERATE_PHKEY_ON_DEVICE) && defined(COMPARE_WITH_DIRECT_SOLVER)
      ibody_direct_dev     = ibody0_dev;
      ibody_direct_dev.acc = direct_dev;
#ifdef  GADGET_MAC
      ibody_direct_dev.acc_old = direct_old_dev;
#endif//GADGET_MAC
#endif//defined(GENERATE_PHKEY_ON_DEVICE) && defined(COMPARE_WITH_DIRECT_SOLVER)
      dumpSnapshot
	(unit, num, ibody0_dev, ibody0,
#ifdef  USE_HDF5_FORMAT
	 &body_snapshot, hdf5type,
#ifdef  MONITOR_ENERGY_ERROR
	 &relEneErr,
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
	 time, steps, file, present, &previous
#ifdef  LEAP_FROG_INTEGRATOR
	 , dt
#endif//LEAP_FROG_INTEGRATOR
#ifndef SERIALIZED_EXECUTION
	 , &iocfg, Ntot
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
	 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
#ifdef  COMPARE_WITH_DIRECT_SOLVER
	 , Ni, ibody_direct_dev, devProp, direct, Ngrp, laneInfo_dev
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	 , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
	 , accfile
#endif//COMPARE_WITH_DIRECT_SOLVER
	 );
      //-------------------------------------------------------------------
      rebuild.reuse = 0;
#ifdef  BLOCK_TIME_STEP
      rebuild.adjust = false;
#endif//BLOCK_TIME_STEP
      //-------------------------------------------------------------------
    }
    //---------------------------------------------------------------------
#else///EXEC_BENCHMARK
#       ifdef  BLOCK_TIME_STEP
    previous = present;
#       endif//BLOCK_TIME_STEP
#endif//EXEC_BENCHMARK
    //---------------------------------------------------------------------
  }/* while( time < ft ){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  COMPARE_WITH_DIRECT_SOLVER
  //-----------------------------------------------------------------------
  {
    //---------------------------------------------------------------------
    int currentNum = 0;
    static char filename[128];
    FILE *fp;
    sprintf(filename, "%s/%s.%s.txt", DOCUMENTFOLDER, file, ERRFILE_NUM);
    if( 0 == access(filename, F_OK) ){
      fp = fopen(filename, "r");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
      int checker = 1;
      checker &= (1 == fscanf(fp, "%d", &currentNum));
      fclose(fp);
      if( !checker ){	__KILL__(stderr, "ERROR: read failure from \"%s\"\n", filename);      }
    }/* if( 0 == access(filename, F_OK) ){ */
    //---------------------------------------------------------------------
    currentNum++;
    fp = fopen(filename, "w");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    fprintf(fp, "%d\n", currentNum);
    fclose(fp);
    //---------------------------------------------------------------------
    sprintf(filename, "%s/%s.%s.txt", DOCUMENTFOLDER, file, ERRFILE_LST);
    fp = fopen(filename, "a");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    fprintf(fp, "%s\t%u\t%s\t%.13e\n", accfile, 1 + previous,
#ifdef  GADGET_MAC
#ifdef  YMIKI_MAC
	    "VecAcc_MAC", absErr
#else///YMIKI_MAC
	    "GADGET2_MAC", absErr
#endif//YMIKI_MAC
#else///GADGET_MAC
#ifdef  WS93_MAC
	    "Multipole_MAC", accErr
#else///WS93_MAC
	    "Opening_angle", theta
#endif//WS93_MAC
#endif//GADGET_MAC
	    );
    fclose(fp);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#endif//COMPARE_WITH_DIRECT_SOLVER
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#   if  !defined(EXEC_BENCHMARK) && !defined(COMPARE_WITH_DIRECT_SOLVER)
  //-----------------------------------------------------------------------
  /* output final stage of numerical results */
  //-----------------------------------------------------------------------
  dumpRestarter
    (num, ibody0_dev, ibody0, time, dt, steps, file, &last, &formerTime
#ifdef  USE_HDF5_FORMAT
     , hdf5type
     , rebuild, elapsed
#ifdef  WALK_TREE_COMBINED_MODEL
     , rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
     , brentDistance, brentHistory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
     , relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
#ifndef SERIALIZED_EXECUTION
     , &iocfg, mpi, Ntot
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
  //-----------------------------------------------------------------------
#endif//!defined(EXEC_BENCHMARK) && !defined(COMPARE_WITH_DIRECT_SOLVER)
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  freeSnapshotArray(hdf5_pos, hdf5_vel, hdf5_acc, hdf5_m, hdf5_pot, hdf5_idx);
#endif//USE_HDF5_FORMAT
#ifdef  GENERATE_PHKEY_ON_DEVICE
  freePeanoHilbertKey_dev
    (tag_dev, peano_dev, peano, min_dev, max_dev, gsync_ph0, gsync_ph1
#ifndef SERIALIZED_EXECUTION
     , box_min, box_max, min_hst, max_hst
#endif//SERIALIZED_EXECUTION
#ifndef CALC_MULTIPOLE_ON_DEVICE
     , list
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef  CUB_AVAILABLE
     , phsort_temp_storage, tag_pre, peano_pre
#endif//CUB_AVAILABLE
     );
#else///GENERATE_PHKEY_ON_DEVICE
  freePeanoHilbertKey
    (peano, tag
#ifndef CALC_MULTIPOLE_ON_DEVICE
     , list
#endif//CALC_MULTIPOLE_ON_DEVICE
     );
#endif//GENERATE_PHKEY_ON_DEVICE
#ifdef  BLOCK_TIME_STEP
  freeParticleDataSoA_hst(idx0    , pos0    , acc0    , vel0    , time0);
  freeParticleDataSoA_hst(idx1    , pos1    , acc1    , vel1    , time1);
  freeParticleDataSoA_dev(idx0_dev, pos0_dev, acc0_dev, vel0_dev, time0_dev,
#       ifdef  GENERATE_PHKEY_ON_DEVICE
			  idx1_dev, pos1_dev, acc1_dev, vel1_dev, time1_dev
#       else///GENERATE_PHKEY_ON_DEVICE
			  jpos_dev, jvel_dev
#       endif//GENERATE_PHKEY_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
			  , neighbor_dev
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
			  , encBall_dev, encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  DPADD_FOR_ACC
			  , tmp_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
			  , res_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
			  );
#else///BLOCK_TIME_STEP
  freeParticleDataSoA_hst(idx0    , pos0    , acc0    , vx0    , vy0    , vz0    );
  freeParticleDataSoA_hst(idx1    , pos1    , acc1    , vx1    , vy1    , vz1    );
  freeParticleDataSoA_dev(idx0_dev, pos0_dev, acc0_dev, vx0_dev, vy0_dev, vz0_dev
#       ifdef  GENERATE_PHKEY_ON_DEVICE
			  , idx1_dev, pos1_dev, acc1_dev, vx1_dev, vy1_dev, vz1_dev
#       endif//GENERATE_PHKEY_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
			  , neighbor_dev
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
			  , encBall_dev, encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  DPADD_FOR_ACC
			  , tmp_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
			  , res_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
			  );
#endif//BLOCK_TIME_STEP
  freeParticleInfoSoA_hst(
#ifndef MAKE_TREE_ON_DEVICE
   jtag
#ifdef  COUNT_INTERACTIONS
   , Nj, Nbuf
#endif//COUNT_INTERACTIONS
#else///MAKE_TREE_ON_DEVICE
#ifdef  COUNT_INTERACTIONS
   Nj, Nbuf
#endif//COUNT_INTERACTIONS
#endif//MAKE_TREE_ON_DEVICE
/* jtag */
/* #ifdef  COUNT_INTERACTIONS */
/* 			  , Nj, Nbuf */
/* #endif//COUNT_INTERACTIONS */
			  );
  freeParticleInfoSoA_dev(jtag_dev
#ifdef  COUNT_INTERACTIONS
			  , Nj_dev, Nbuf_dev
#endif//COUNT_INTERACTIONS
			  );
  freeParticleGroups(laneInfo_hst, laneInfo_dev, laneTime_dev
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
		     , inum_hst, inum_dev
#   if  defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
		     , neighbor_temp_storage, neighbor_cubout
#endif//defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
		     , ibody_dt_dev, ibody_dt_hst, dtInfo
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
		     );
#ifndef MAKE_TREE_ON_DEVICE
  freeTreeCell(hkey, parent, children
#ifndef CALC_MULTIPOLE_ON_DEVICE
	       , cell, leaf, node
#endif//CALC_MULTIPOLE_ON_DEVICE
	       );
#endif//MAKE_TREE_ON_DEVICE
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  freeTreeCell_dev(cell_dev, leaf_dev, node_dev, list_dev,
#ifdef  MAKE_TREE_ON_DEVICE
		   hkey_dev, parent_dev, children_dev, bottomLev_dev, numCell_dev, numNode_dev, scanNum_dev
#else///MAKE_TREE_ON_DEVICE
                   cell    , leaf    , node, list
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
                   , cell, leaf, node, list
#endif//defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
		   );
#endif//CALC_MULTIPOLE_ON_DEVICE
  freeTreeNode_dev(more_dev, pj_dev, mj_dev,
#ifdef  CALC_MULTIPOLE_ON_DEVICE
		   bmax_dev, node2cell_dev, gsync0, gsync1,
#       ifdef  WS93_MAC
		   mr2_dev,
#       endif//WS93_MAC
#else///CALC_MULTIPOLE_ON_DEVICE
		   bmax,
#       ifdef  WS93_MAC
		   mr2    ,
#       endif//WS93_MAC
#endif//CALC_MULTIPOLE_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
		   niSub_dev,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
		   more,
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
		   pj, mj,
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
#ifdef  MAKE_TREE_ON_DEVICE
		   gmem_make_tree_dev, gsync0_make_tree_dev, gsync1_make_tree_dev, gsync2_make_tree_dev, gsync3_make_tree_dev,
		   gmem_link_tree_dev, gsync0_link_tree_dev, gsync1_link_tree_dev,
#else///MAKE_TREE_ON_DEVICE
		   node2cell,
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
		   niSub_hst,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#endif//MAKE_TREE_ON_DEVICE
#ifdef  GADGET_MAC
		   mac_dev,
#endif//GADGET_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
		   gmem_external_dev, gsync0_external_dev, gsync1_external_dev, diameter_dev, diameter_hst,
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
		   more0Buf, more1Buf, rjmaxBuf, makeFail);
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  freeNeighborSearch_dev(gsync_ns0, gsync_ns1, freeLst_ns
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
			 , freeNum_ns, active_ns
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
			 );
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  freeTreeBuffer_dev(fail_dev, buffer, freeLst
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
		     , freeNum, active
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#ifndef USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
		     , cycles_hst, cycles_dev
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
		     , cycles_let_hst, cycles_let_dev
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#endif//USE_CUDA_EVENT
		     );
#ifdef  COMPARE_WITH_DIRECT_SOLVER
  freeAccel_dev(direct_dev, direct
#ifdef  GADGET_MAC
		, direct_old_dev
#endif//GADGET_MAC
		);
#endif//COMPARE_WITH_DIRECT_SOLVER
#ifndef BLOCK_TIME_STEP
  freeTimeStep_dev(dt_dev);
#endif//BLOCK_TIME_STEP
/* #   if  defined(CUB_AVAILABLE) && defined(BLOCK_TIME_STEP) */
/*   freeTimeStep_dev(time_temp_storage, laneInfo_cub, laneTime_cub); */
/* #endif//defined(CUB_AVAILABLE) && defined(BLOCK_TIME_STEP) */
#ifndef SERIALIZED_EXECUTION
  releaseLETtopology(nodeInfo, ipos,
#ifdef  GADGET_MAC
		     amin,
#endif//GADGET_MAC
#ifdef  BUILD_LET_ON_DEVICE
		     numSend_hst, numSend_dev,
#endif//BUILD_LET_ON_DEVICE
#ifdef  ALLOCATE_LETBUFFER
		     buf4LET,
#endif//ALLOCATE_LETBUFFER
		     stream_let, Nstream_let);
#ifdef  EXCHANGE_USING_GPUS
  releaseORMtopology(dxmin, dxmax, dymin, dymax, dzmin, dzmax, dmreq,
		     sxmin, sxmax, symin, symax, szmin, szmax,
		     iparticleSendBuf, iparticleRecvBuf, sampleRecvNum, sampleRecvDsp,
		     ormCfg, repCfg, letcfg.rank);
  /* releaseBoxSize_dev(pmin_hst, pmax_hst, pmin_dev, pmax_dev, gsync_box0, gsync_box1); */
  releaseSamplePos(x0hst, x1hst, y0hst, y1hst, z0hst, z1hst, idhst,
		   x0dev, x1dev, y0dev, y1dev, z0dev, z1dev, iddev);
  releaseParticlePosition(xhst, yhst, zhst, xdev, ydev, zdev, rank_hst, rank_dev, idx_dev);
#else///EXCHANGE_USING_GPUS
  free(domDecKey);
  removeORBtopology(domCfg, iparticleSendBuf, iparticleRecvBuf, ndim, orbCfg,
		    sampleLoc, sampleFul, sampleRecvNum, sampleRecvDsp, domainMin, domainMax);
#endif//EXCHANGE_USING_GPUS
#ifdef  BUILD_LET_ON_DEVICE
  freeGeometricEnclosingBall_dev(r2geo_dev
/* #ifdef  CUB_AVAILABLE */
/* 				 , r2geo_temp_storage */
/* #endif//CUB_AVAILABLE */
				 );
#endif//BUILD_LET_ON_DEVICE
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
  /* destroy CUDA events */
#   if  defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
  releaseCUDAevents_dev(iniCalcAcc, finCalcAcc
#ifdef  MONITOR_LETGEN_TIME
			, iniMakeLET, finMakeLET
#endif//MONITOR_LETGEN_TIME
#ifndef SERIALIZED_EXECUTION
			, mpi.size
#else///SERIALIZED_EXECUTION
			, 1
#endif//SERIALIZED_EXECUTION
			);
#endif//defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
  //-----------------------------------------------------------------------
  /* destroy CUDA streams */
  for(int ii = 0; ii < sinfo.num; ii++){
    mycudaStreamDestroy(stream[ii]);
  }/* for(int ii = 0; ii < sinfo.num; ii++){ */
  free(stream);
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  REPORT_TOTAL_ELAPSED_TIME
  steps -= stepsInit;
#ifndef SERIALIZED_EXECUTION
  MPI_Barrier(mpi.comm);
#endif//SERIALIZED_EXECUTION
  static struct timeval timeExit;
  gettimeofday(&timeExit, NULL);
  double totalElapsed = calcElapsedTimeInSec(timeInit, timeExit);
#ifndef SERIALIZED_EXECUTION
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &totalElapsed, 1, MPI_DOUBLE, MPI_MAX, mpi.comm));
  if( mpi.rank == 0 )
#endif//SERIALIZED_EXECUTION
    {
      //-------------------------------------------------------------------
      fprintf(stdout, "# %s is used with accuracy controlling parameter of %e.\n",
#ifdef  GADGET_MAC
	      "Acceleration MAC", absErr
#else///GADGET_MAC
#ifdef  WS93_MAC
	      "Multipole MAC", accErr
#else///WS93_MAC
	      "Opening criterion", theta
#endif//WS93_MAC
#endif//GADGET_MAC
	      );
      fprintf(stdout, "# Total elapsed time is %e sec., Ntot = %zu, %zu steps in total.\n", totalElapsed, Ntot, steps);
      fprintf(stdout, "# Elapsed time per step is %e sec.\n", totalElapsed / (double)steps);
      //-------------------------------------------------------------------
      FILE *fp;
      static char filename[128];
      sprintf(filename, "%s/%s.%s.%s.cc%d.%zu.time.log", LOGFOLDER, file,
#ifdef  GADGET_MAC
#ifdef  YMIKI_MAC
	      "vector_acc",
#else///YMIKI_MAC
	      "acceleration",
#endif//YMIKI_MAC
#else///GADGET_MAC
#ifdef  WS93_MAC
	      "multipole",
#else///WS93_MAC
	      "theta",
#endif//WS93_MAC
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
	      "block",
#else///BLOCK_TIME_STEP
	      "share",
#endif//BLOCK_TIME_STEP
	      GPUVER, Ntot);
      fp = fopen(filename, "a");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
      fprintf(fp, "%.13e\t%e\t%zu\t%e\n",
#ifdef  GADGET_MAC
	      absErr,
#else///GADGET_MAC
#ifdef  WS93_MAC
	      accErr,
#else///WS93_MAC
	      theta,
#endif//WS93_MAC
#endif//GADGET_MAC
	      totalElapsed, steps, totalElapsed / (double)steps);
      fclose(fp);
      //-------------------------------------------------------------------
    }
#endif//REPORT_TOTAL_ELAPSED_TIME
  //-----------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
  exitMPI();
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef COUNT_INTERACTIONS
//-------------------------------------------------------------------------
void analyzeWalkStatistics(int Ni, const int Ngrp, int * restrict results, walk_stats *stats)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  stats->Ngroup = Ngrp;
  stats->min    = INT_MAX;
  stats->max    = 0;
  stats->mean   = 0.0f;
  stats->sdev   = 0.0f;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ni; ii++){
    //---------------------------------------------------------------------
    const int dat = results[ii];
    //---------------------------------------------------------------------
    if( dat != 0 ){
      if( stats->min > dat )	stats->min = dat;
      if( stats->max < dat )	stats->max = dat;
    }/* if( dat != 0 ){ */
    //---------------------------------------------------------------------
    stats->mean += dat;
    stats->sdev += dat * dat;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ni; ii++){ */
  //-----------------------------------------------------------------------
  if( (Ngrp * (TSUB / NWARP)) < Ni )
    Ni = Ngrp * (TSUB / NWARP);
  const float inv = 1.0f / (float)Ni;
  stats->mean *= inv;
  stats->sdev *= inv;
  stats->sdev -= stats->mean * stats->mean;
  stats->sdev  = sqrtf(stats->sdev);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void analyzeTreeStatistics(tree_stats *stats)
{
  //-----------------------------------------------------------------------
  const real inv = UNITY / (real)stats->nodeNum;
  //-----------------------------------------------------------------------
  stats->mjMean *= inv;  stats->r2Mean *= inv;  stats->mjSdev -= stats->mjMean * stats->mjMean;  stats->mjSdev = SQRT(stats->mjSdev);
  stats->mjSdev *= inv;  stats->r2Sdev *= inv;  stats->r2Sdev -= stats->r2Mean * stats->r2Mean;  stats->r2Sdev = SQRT(stats->r2Sdev);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void analyzeTreeMetrics(tree_metrics *metric)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(uint ii = 0; ii < MAXIMUM_PHKEY_LEVEL; ii++){
    //---------------------------------------------------------------------
    if( metric->level[ii].nodeNum != 0 ){
      //-------------------------------------------------------------------
      metric->total.nodeNum += metric->level[ii].nodeNum;
      metric->total.cellNum += metric->level[ii].cellNum;
      metric->total.mjMean  += metric->level[ii].mjMean;
      metric->total.mjSdev  += metric->level[ii].mjSdev;
      metric->total.r2Mean  += metric->level[ii].r2Mean;
      metric->total.r2Sdev  += metric->level[ii].r2Sdev;
      //-------------------------------------------------------------------
      analyzeTreeStatistics(&metric->level[ii]);
      //-------------------------------------------------------------------
    }
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  analyzeTreeStatistics(&metric->total);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//COUNT_INTERACTIONS
//-------------------------------------------------------------------------
