/**
 * @file gothic.c
 *
 * @brief Source code for GOTHIC (Gravitational Oct-Tree code accelerated by HIerarchical time step Controlling)
 * @detail See Miki & Umemura (2017), New Astronomy, 52, 65--81 for implementation details
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/11/21 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

/**
 * @def OUTPUT_EXCHANGE_SIGNAL
 *
 * @brief a switch to output signals of particle exchanging
 */
/* #define OUTPUT_EXCHANGE_SIGNAL */

#ifdef  SERIALIZED_EXECUTION
#define OUTPUT_MEMORY_USAGE
#else///SERIALIZED_EXECUTION
/* #define MONITOR_SIMULATION_STATUS */
/* #define SHARED_AUTO_TUNER */

#ifndef DISABLE_AUTO_TUNING
/**
 * @def DISABLE_EXCG_BODIES_BEFORE_SATURATION
 *
 * @brief disable particle exchanging before auto-tuning reaches saturation
 */
#define DISABLE_EXCG_BODIES_BEFORE_SATURATION
#endif//DISABLE_AUTO_TUNING
#endif//SERIALIZED_EXECUTION

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>/**< for on-the-fly monitoring of execution time of multiple functions */
#include <unistd.h>/**< to check the existence of files */
#include <math.h>/**< for auto-tuning */
#include <mpi.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "myutil.h"
#include "name.h"
#include "constants.h"
#include "timer.h"
#include "cudalib.h"

#ifndef SERIALIZED_EXECUTION
#include "mpilib.h"
#include "../para/mpicfg.h"
#endif//SERIALIZED_EXECUTION

#include "../misc/device.h"
#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../misc/allocate_dev.h"
#include "../misc/convert.h"
#include "../misc/tune.h"
#include "../misc/brent.h"

#include "../file/io.h"

#include "../sort/peano.h"
#include "../sort/peano_dev.h"

#ifndef SERIALIZED_EXECUTION
#include "../para/exchange.h"
#include "../para/exchange_dev.h"
#endif//SERIALIZED_EXECUTION

#include "../tree/make.h"
#include "../tree/let.h"
#include "../tree/buf_inc.h"
#include "../tree/make_dev.h"
#include "../tree/walk_dev.h"
#include "../tree/neighbor_dev.h"
#include "../tree/shrink_dev.h"
#ifndef SERIALIZED_EXECUTION
#include "../tree/geo_dev.h"
#include "../tree/let_dev.h"
#ifdef  USE_ENCLOSING_BALL_FOR_LET
#include "../tree/icom_dev.h"
#endif//USE_ENCLOSING_BALL_FOR_LET
#endif//SERIALIZED_EXECUTION

#include "../time/adv_dev.h"

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
#include "../tree/potential_dev.h"
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#define LEAP_FROG_INTEGRATOR
#   if  defined(LEAP_FROG_INTEGRATOR) && defined(BLOCK_TIME_STEP)
#undef          LEAP_FROG_INTEGRATOR
#endif//defined(LEAP_FROG_INTEGRATOR) && defined(BLOCK_TIME_STEP)


#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
real theta2;
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)


#   if  defined(HUNT_MAKE_PARAMETER) || defined(HUNT_FIND_PARAMETER)
int treeBuildCalls = 0;
#endif//defined(HUNT_MAKE_PARAMETER) || defined(HUNT_FIND_PARAMETER)


#ifdef  SWITCH_WITH_J_PARALLELIZATION
static inline void selectCommunicationMode
(const int grpNum, const int totNum, laneinfo *laneInfo_hst, const MPIcfg_tree mpi,
 bool *transferMode, int * restrict Ni_local, int * restrict Ni_total, int * restrict Ni_list, int * restrict head_list, int * restrict grpNum_list)
{
  __NOTE__("%s\n", "start");

  int transfer = ((float)grpNum < (FCRIT_J_PARALLELIZATION * (float)totNum)) ? 1 : 0;
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &transfer, 1, MPI_INT, MPI_SUM, mpi.comm));

  *transferMode = ((float)transfer > (FSIZE_J_PARALLELIZATION * (float)mpi.size));
  if( *transferMode ){
    __NOTE__("%s\n", "i-particle transfer mode");

    *Ni_local = (grpNum > 0) ? (laneInfo_hst[grpNum - 1].head + laneInfo_hst[grpNum - 1].num) : 0;
    chkMPIerr(MPI_Allgather(Ni_local, 1, MPI_INT, Ni_list, 1, MPI_INT, mpi.comm));
    chkMPIerr(MPI_Allgather(&grpNum, 1, MPI_INT, grpNum_list, 1, MPI_INT, mpi.comm));
    head_list[0] = 0;
    for(int ii = 1; ii < mpi.size; ii++)
      head_list[ii] = head_list[ii - 1] + Ni_list[ii - 1];
    *Ni_total = head_list[mpi.size - 1] + Ni_list[mpi.size - 1];
    __NOTE__("Ni_local = %d, Ni_total = %d\n", *Ni_local, *Ni_total);

    if( *Ni_total > NMAX_J_PARALLELIZATION )
      *transferMode = false;/**< #this is a tentative treatment before developing double buffer mode for i-particle transfer mode in src/tree/walk_dev.cu */
  }/* if( *transferMode ){ */

  __NOTE__("%s\n", "end");
}
#endif//SWITCH_WITH_J_PARALLELIZATION


/**
 * @fn setMultipoleMoment
 *
 * @brief Call function to calculate multipole moment on GPU.
 *
 * @sa calcMultipole_dev
 */
static inline void setMultipoleMoment
(const int bottomLev, const soaTreeCell cell_dev, const int numNode, const soaTreeNode node_dev
 , const int num, const iparticle ibody, const soaMakeTreeBuf buf, const deviceProp devProp
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
 , const soaTreeCell cell_hst, tree_stats *level
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
 )
{
  __NOTE__("%s\n", "start");


  /** calculate multipole moment of tree nodes */
#ifndef SERIALIZED_EXECUTION
  double tmac = 0.0;
#endif//SERIALIZED_EXECUTION

  calcMultipole_dev(bottomLev, cell_dev,
		    num, ibody, numNode, node_dev, buf, devProp
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
		    , eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifndef SERIALIZED_EXECUTION
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


  __NOTE__("%s\n", "end");
}


#ifndef SERIALIZED_EXECUTION
/**
 * @fn updateDomain
 *
 * @brief Call function to exchanging N-body particles among MPI processes.
 *
 * @sa initStatVal
 * @sa initGuessTime
 * @sa exchangeParticles_dev
 * @sa exchangeParticles
 */
static inline void updateDomain
(int *num, const int num_max, int *Ni,
 iparticle *ibody0_dev, iparticle *ibody1_dev, iparticle *ibody0_hst, iparticle *ibody1_hst, const ulong Ntot,
 sendDom domBoundary, particlePos pos_hst, particlePos pos_dev, domainCfg domain, domainDecomposeKey key,
 sampling sample, samplePos loc, samplePos ful, soaPHsort soa, const deviceProp devProp, const deviceInfo devInfo,
#ifndef DISABLE_AUTO_TUNING
 autoTuningParam *exchangeParam, double *exchangeInterval,
#endif//DISABLE_AUTO_TUNING
 MPIinfo orm[restrict], MPIinfo rep[restrict],
#ifdef  CARE_EXTERNAL_PARTICLES
 domainLocation *location,
#endif//CARE_EXTERNAL_PARTICLES
 const MPIcfg_tree letcfg, sendCfg *iparticleSendBuf, recvCfg *iparticleRecvBuf, measuredTime *measured
 , brentStatus *status, brentMemory *memory
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
)
{
  __NOTE__("%s\n", "start");


  /** clear counters */
#ifndef DISABLE_AUTO_TUNING
  *exchangeInterval = 0.0;
  initStatVal(&((*exchangeParam).linearStats));  initGuessTime(&((*exchangeParam).linearGuess));
  initStatVal(&((*exchangeParam). powerStats));  initGuessTime(&((*exchangeParam). powerGuess));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  initStatVal  (&((*exchangeParam).parabolicStats));
  initGuessTime(&((*exchangeParam).parabolicGuess));
#endif//USE_PARABOLIC_GROWTH_MODEL
#endif//DISABLE_AUTO_TUNING

  /** execute particle exchanging */
  exchangeParticles_dev
    (*num, Ntot, num_max, num, ibody0_dev, ibody1_dev, ibody0_hst, ibody1_hst,
     domBoundary, pos_hst, pos_dev, key, iparticleSendBuf, iparticleRecvBuf,
     orm, rep, domain, letcfg,
     measured->sum_excg,
     sample, loc, ful, soa, devProp, devInfo, measured
#ifdef  CARE_EXTERNAL_PARTICLES
     , location
#endif//CARE_EXTERNAL_PARTICLES
     , status, memory
#ifdef  EXEC_BENCHMARK
     , execTime
#endif//EXEC_BENCHMARK
     );

  *Ni = *num;


  __NOTE__("%s\n", "end");
}
#endif//SERIALIZED_EXECUTION


/**
 * @fn buildTreeStructure
 *
 * @brief Build tree structure on GPU.
 *
 * @sa sortParticlesPHcurve_dev
 * @sa makeTreeStructure_dev
 */
static inline void buildTreeStructure
(int num, iparticle *ibody0_dev, iparticle *ibody1_dev, const soaPHsort soaPH_dev, const deviceProp devProp,
#ifdef  CUB_AVAILABLE
 const soaPHsort soaPH_pre,
#endif//CUB_AVAILABLE
#ifndef SERIALIZED_EXECUTION
 int *numOld,
#endif//SERIALIZED_EXECUTION
 int *leafLev, int *numCell, int *numNode, int *leafLev_dev, int *scanNum_dev, int *numCell_dev, int *numNode_dev, const soaMakeTreeBuf buf, const soaTreeCell cell_dev, const soaTreeNode node_dev
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
   , domainLocation *location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
 , struct timespec *start
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
 )
{
  __NOTE__("%s\n", "start");


#   if  defined(HUNT_MAKE_PARAMETER) || defined(HUNT_FIND_PARAMETER)
  treeBuildCalls++;
#endif//defined(HUNT_MAKE_PARAMETER) || defined(HUNT_FIND_PARAMETER)

  /** sort N-body particles in Peano-Hilbert order */
  sortParticlesPHcurve_dev(num, ibody0_dev, ibody1_dev, soaPH_dev, devProp
#ifdef  CUB_AVAILABLE
			   , soaPH_pre
#endif//CUB_AVAILABLE
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
			   , location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
			   , start
#ifdef  EXEC_BENCHMARK
			   , execTime
#endif//EXEC_BENCHMARK
			   );

  /** build a breadth-first octree structure */
  makeTreeStructure_dev(num, soaPH_dev.key, leafLev, leafLev_dev, scanNum_dev, numCell, numCell_dev, cell_dev, numNode, numNode_dev, node_dev, buf, devProp
#ifndef SERIALIZED_EXECUTION
			, *numOld
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
			, execTime
#endif//EXEC_BENCHMARK
			);
#ifndef SERIALIZED_EXECUTION
  __NOTE__("num = %d, numOld = %d\n", num, *numOld);
  *numOld = num;
#endif//SERIALIZED_EXECUTION

#ifndef NDEBUG
  if( *leafLev == MAXIMUM_PHKEY_LEVEL )
    printPHkey_location(num, soaPH_dev.key, *ibody0_dev);
#endif//NDEBUG

  __NOTE__("%s\n", "end");
}


/**
 * @fn configDistribution
 *
 * @brief Configure particle distribution.
 *
 * @sa brentCalc1st
 * @sa updateParticleGroups
 * @sa commitParticleGroups
 * @sa setLaneTime_dev
 */
static inline void configDistribution
(const int num, const int inumPerLane, int *Ngrp, const int maxNgrp, laneinfo *laneInfo_hst, laneinfo *laneInfo_dev, iparticle ibody_dev, brentStatus *brent, int *inum_dev, int *inum_hst
#ifdef  BLOCK_TIME_STEP
 , double *laneTime_dev
#endif//BLOCK_TIME_STEP
 , const struct timespec start, measuredTime *time
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
 )
{
  __NOTE__("%s\n", "start");


#ifndef DISABLE_AUTO_TUNING
  static bool initialized = false;
  if( initialized )    brentCalc1st(brent, 1.0e-13);
  else                 initialized = true;
#endif//DISABLE_AUTO_TUNING

  updateParticleGroups(num, laneInfo_hst, inumPerLane, maxNgrp, Ngrp, ibody_dev, inum_dev, inum_hst, CAST_D2R(brent->u.pos)
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );

  commitParticleGroups(*Ngrp, laneInfo_hst, laneInfo_dev);

#ifdef  BLOCK_TIME_STEP
  setLaneTime_dev(*Ngrp, laneInfo_dev, laneTime_dev, ibody_dev
#ifdef  EXEC_BENCHMARK
		  , execTime
#endif//EXEC_BENCHMARK
		  );
#endif//BLOCK_TIME_STEP

  static struct timespec finish;
  clock_gettime(CLOCK_MONOTONIC_RAW, &finish);
  time->makeTree = calcElapsedTimeInSec(start, finish);
#ifndef SERIALIZED_EXECUTION
  time->sum_excg    += time->makeTree;
  time->sum_rebuild += time->makeTree;
#endif//SERIALIZED_EXECUTION


  __NOTE__("%s\n", "end");
}


#ifndef EXEC_BENCHMARK
/**
 * @fn dumpSnapshot
 *
 * @brief Write snapshot.
 *
 * @sa copyParticle_dev2hst
 * @sa copyAry4Snapshot
 * @sa backVel4Snapshot
 * @sa backVel
 * @sa writeSnapshot
 * @sa writeSnapshotParallel
 * @sa calcGravity_dev
 * @sa copyAccel_dev2hst
 * @sa writeApproxAccel
 */
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
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
 , gpu_clock *deviceMonitors, int *monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  REPORT_COMPUTE_RATE
 , struct timespec *clock_prev, double *time_prev, ulong *step_prev, const double ft, double *brent_rate, ulong *brent_prev, ulong *rebuild_freq
#endif//REPORT_COMPUTE_RATE
#ifdef  COMPARE_WITH_DIRECT_SOLVER
 , const int Ni, const iparticle ibody_direct_dev, const deviceProp devProp, acceleration *direct_hst, const int Ngrp, laneinfo *laneInfo_dev
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
 , const real eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 , const potential_field sphe
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
 , const disk_potential disk
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
 , char *accfile
#endif//COMPARE_WITH_DIRECT_SOLVER
 )
{
  __NOTE__("%s\n", "start");


  /** set particle data on host */
  copyParticle_dev2hst(num, ibody_dev, ibody_hst
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );

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

  /** report speed of N-body simulation */
#ifdef  REPORT_COMPUTE_RATE
  static struct timespec clock;
  clock_gettime(CLOCK_MONOTONIC_RAW, &clock);
  const double elapsed = calcElapsedTimeInSec(*clock_prev, clock);
  const double proceed = time - (*time_prev);
  const double speed = elapsed / ((present > 0) ? ((double)(steps - (*step_prev))) : (1.0));
  const double speed_run = proceed / elapsed;
  const double guess = (ft - time) / speed_run;
  const double complete = (time / ft) * 100.0;
  const double brent_avg = *brent_rate / ((present > 0) ? ((double)(steps - (*step_prev))) : (1.0));
  const double rebuild_interval = ((present > 0) ? ((double)(steps - (*step_prev))) : (1.0)) / (double)(*rebuild_freq);
  *clock_prev = clock;
  *time_prev = time;
  *step_prev = steps;
  *brent_rate = 0.0;
  *brent_prev = steps;
  *rebuild_freq = 0.0;
#endif//REPORT_COMPUTE_RATE

  /** output the snapshot file */
#ifdef  SERIALIZED_EXECUTION
  writeSnapshot
    (unit, time, steps, num, file, present
#ifdef  USE_HDF5_FORMAT
     , body, hdf5type
#ifdef  MONITOR_ENERGY_ERROR
     , relEneErr
#endif//MONITOR_ENERGY_ERROR
#else///USE_HDF5_FORMAT
     , ibody_hst
#endif//USE_HDF5_FORMAT
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
     , deviceMonitors, *monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  REPORT_COMPUTE_RATE
     , speed, speed_run, complete, guess, brent_avg, rebuild_interval
#endif//REPORT_COMPUTE_RATE
     );
#else///SERIALIZED_EXECUTION
  writeSnapshotParallel
    (unit, time, steps, num, file, present, iocfg, Ntot
#ifdef  USE_HDF5_FORMAT
     , body, hdf5type
#ifdef  MONITOR_ENERGY_ERROR
     , relEneErr
#endif//MONITOR_ENERGY_ERROR
#else///USE_HDF5_FORMAT
     , ibody_hst
#endif//USE_HDF5_FORMAT
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
     , deviceMonitors, *monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  REPORT_COMPUTE_RATE
     , speed, speed_run, complete, guess, brent_avg, rebuild_interval
#endif//REPORT_COMPUTE_RATE
     );
#endif//SERIALIZED_EXECUTION
  *previous = present;
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
  *monitor_step = 0;
#endif//REPORT_GPU_CLOCK_FREQUENCY


#ifdef  COMPARE_WITH_DIRECT_SOLVER
 /** output exact and approximated acceleration and potential */
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.dat", DATAFOLDER, file, BRUTE, present);
  if( 0 != access(filename, F_OK) ){
    static soaTreeNode dummy_node;
    static soaTreeWalkBuf dummy_buf;
    double null;

    calcGravity_dev(
#ifdef  BLOCK_TIME_STEP
		    Ngrp, NULL,
#endif//BLOCK_TIME_STEP
		    Ngrp, laneInfo_dev, ibody_direct_dev, dummy_node, dummy_buf, NULL, devProp, &null
#   if  !defined(BLOCK_TIME_STEP) || defined(COMPARE_WITH_DIRECT_SOLVER) || defined(COUNT_INTERACTIONS) || defined(PRINT_PSEUDO_PARTICLE_INFO)
		    , Ni
#endif//!defined(BLOCK_TIME_STEP) || defined(COMPARE_WITH_DIRECT_SOLVER) || defined(COUNT_INTERACTIONS) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
		    , pot_tbl_sphe
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
		    , pot_tbl_disk
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifndef SERIALIZED_EXECUTION
		    , NULL, 0, 1, NULL, letcfg, NULL, NULL, 0, NULL, NULL, NULL
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
		    , NULL
#endif//COUNT_INTERACTIONS
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

    copyAccel_dev2hst(Ni, ibody_direct_dev.acc, direct_hst);
    writeApproxAccel(time, steps, num,
#ifdef  USE_HDF5_FORMAT
		     (*body).idx,
#endif//USE_HDF5_FORMAT
		     direct_hst, filename);
  }/* if( 0 != access(filename, F_OK) ){ */

  sprintf(filename, "%s/%s.%s.%.3u.dat", DATAFOLDER, file, accfile, present);
  writeApproxAccel(time, steps, num,
#ifdef  USE_HDF5_FORMAT
		   (*body).idx,
#endif//USE_HDF5_FORMAT
		   ibody_hst.acc, filename);
#endif//COMPARE_WITH_DIRECT_SOLVER


  __NOTE__("%s\n", "end");
}


#ifndef COMPARE_WITH_DIRECT_SOLVER
/**
 * @fn dumpRestarter
 *
 * @brief Write files for restarting the simulation.
 *
 * @sa copyParticle_dev2hst
 * @sa backVel
 * @sa writeTentativeData
 * @sa writeTentativeDataParallelo
 * @sa updateConfigFile
 */
static inline void dumpRestarter
(const int num, iparticle ibody_dev, iparticle ibody_hst, const double time, const double dt, const ulong steps, const double elapsed, char *file, int *last, double *formerTime
#ifdef  USE_HDF5_FORMAT
 , hdf5struct hdf5type, rebuildTree rebuild, measuredTime measured, autoTuningParam rebuildParam, brentStatus status, brentMemory memory
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
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
 , gpu_clock *deviceMonitors, int *monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
 )
{
  __NOTE__("%s\n", "start");


  /** output tentative results of the current running simulation */
  /** set particle data on host */
  copyParticle_dev2hst(num, ibody_dev, ibody_hst
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );
#ifdef  LEAP_FROG_INTEGRATOR
  backVel(0, num, body, (real)dt);
#endif//LEAP_FROG_INTEGRATOR

  /** output the dump file */
#ifdef  SERIALIZED_EXECUTION
  writeTentativeData        (time, dt, steps, elapsed, num, ibody_hst, file, last
#ifdef  USE_HDF5_FORMAT
			     , hdf5type, rebuild, measured, rebuildParam, status, memory
#ifdef  MONITOR_ENERGY_ERROR
			     , relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
			     , true, deviceMonitors, *monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
			     );
#else///SERIALIZED_EXECUTION
  writeTentativeDataParallel(time, dt, steps, elapsed, num, ibody_hst, file, last, iocfg, Ntot
#ifdef  USE_HDF5_FORMAT
			     , hdf5type, rebuild, measured, rebuildParam, status, memory
#ifdef  MONITOR_ENERGY_ERROR
			     , relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
			     , true, deviceMonitors, *monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
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


  __NOTE__("%s\n", "end");
}
#endif//COMPARE_WITH_DIRECT_SOLVER
#endif//EXEC_BENCHMARK


#ifdef  COUNT_INTERACTIONS
void analyzeWalkStatistics(int Ni, const int Ngrp, int * restrict results, walk_stats *stats);
void analyzeTreeMetrics(tree_metrics *metric);
#endif//COUNT_INTERACTIONS


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
#ifndef SERIALIZED_EXECUTION
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);
  static MPIcfg_dataio iocfg;
  createMPIcfg_dataio(&iocfg, mpi);
#endif//SERIALIZED_EXECUTION

#ifdef  REPORT_TOTAL_ELAPSED_TIME
  static struct timespec timeInit;
  clock_gettime(CLOCK_MONOTONIC_RAW, &timeInit);
#endif//REPORT_TOTAL_ELAPSED_TIME


  /** configure the details of the numerical simulation */
  /** read command line arguments */
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
#ifdef  DISABLE_AUTO_TUNING
    __FPRINTF__(stderr, "          -rebuild_interval=<real>\n");
    __FPRINTF__(stderr, "          -brent_frac=<real>\n");
#endif//DISABLE_AUTO_TUNING
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
    __FPRINTF__(stderr, "          -pot_file_sphe=<char *>\n");
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
    __FPRINTF__(stderr, "          -pot_file_disk=<char *>\n");
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    __FPRINTF__(stderr, "          -dropPrevTune=<int> (optional)\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }
  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv,  "file", &file));
  double tmp;
#ifdef  GADGET_MAC
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "absErr", &tmp));  const real absErr = CAST_D2R(tmp);
#else///GADGET_MAC
#ifdef  WS93_MAC
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "accErr", &tmp));  const real accErr = CAST_D2R(tmp);
#else///WS93_MAC
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv,  "theta", &tmp));  const real  theta = CAST_D2R(tmp);
  theta2 = theta * theta;
#endif//WS93_MAC
#endif//GADGET_MAC
#ifdef  DISABLE_AUTO_TUNING
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "rebuild_interval", &tmp));  const double rebuild_interval = tmp;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv,       "brent_frac", &tmp));  const double       brent_frac = tmp;
#endif//DISABLE_AUTO_TUNING
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  char *pot_file_sphe;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv,  "pot_file_sphe", &pot_file_sphe));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
  char *pot_file_disk;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv,  "pot_file_disk", &pot_file_disk));
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  EXEC_BENCHMARK
  static int jobID;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "jobID", &jobID));
#endif//EXEC_BENCHMARK
  static int dropPrevTune;
  if( optionalCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "dropPrevTune", &dropPrevTune)) != myUtilAvail )
    dropPrevTune = 0;

#ifndef SERIALIZED_EXECUTION
  static int Nx;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "Nx", &Nx));
  static int Ny;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "Ny", &Ny));
  static int Nz;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "Nz", &Nz));
  if( (Nx * Ny * Nz) != mpi.size ){
    __KILL__(stderr, "ERROR: Invalid input arguments: Nx = %d, Ny = %d, Nz = %d while mpi.size = %d\n", Nx, Ny, Nz, mpi.size);
  }/* if( (Nx * Ny * Nz) != mpi.size ){ */
#endif//SERIALIZED_EXECUTION


  /** read settings about the simulation */
  static ulong Ntot;
  static real eps, eta;
  static double ft, snapshotInterval, saveInterval;
  static int unit;
#ifdef  SERIALIZED_EXECUTION
  readSettings        (&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file);
#else///SERIALIZED_EXECUTION
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file, mpi);
#endif//SERIALIZED_EXECUTION

#if 0
  saveInterval = 10.0;
#endif

  /** read setting to dump tentative results of the simulation */
  static int last;
#ifdef  SERIALIZED_EXECUTION
  readConfigFile        (&last, file);
#else///SERIALIZED_EXECUTION
  readConfigFileParallel(&last, file, mpi);
#endif//SERIALIZED_EXECUTION


  /** activate the attatched accelerator device(s) */
  int devIdx;
  deviceInfo devInfo;
  deviceProp devProp;
  int headIdx;
  if( optionalCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "deviceID", &headIdx)) != myUtilAvail )
    headIdx = 0;
#ifdef  SERIALIZED_EXECUTION
  openAcceleratorGPUs(&devIdx, &devInfo, &devProp, headIdx, 1, 0, 1);
#else///SERIALIZED_EXECUTION
  openAcceleratorGPUs(&devIdx, &devInfo, &devProp, headIdx, 1, mpi.rank, mpi.size);
#endif//SERIALIZED_EXECUTION
  __NOTE__("devIdx = %d\n", devIdx);

  /** set CUDA streams */
  cudaStream_t *stream;
  kernelStream sinfo;
  setCUDAstreams_dev(&stream, &sinfo, &devInfo);

  /** set the number of N-body particles */
  /** this procedure is obvious for this serial process execution, but useful for the future upgrade (parallelization) */
  /** assume num >= Ni */
  static int num, Ni;

#ifdef  SERIALIZED_EXECUTION
  num = Ni = (int)Ntot;
  const int num_max = num;
#else///SERIALIZED_EXECUTION
  static MPIcfg_tree letcfg;
  setNodeConfig(Ntot, &num, &Ni, mpi, &letcfg
#   if  !defined(NDEBUG) || defined(REPORT_GPU_MPI_BINDING)
		, devIdx
#endif//!defined(NDEBUG) || defined(REPORT_GPU_MPI_BINDING)
		);
  const int num_max = (int)ceilf((float)num * MAX_FACTOR_FROM_EQUIPARTITION);
#endif//SERIALIZED_EXECUTION


  /** set global constants */
  /** set gravitational softening and opening criterion */
  extern real newton;
  setGlobalConstants_walk_dev_cu
    (newton, eps * eps
#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
     , theta2
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
     );
  setGlobalConstants_make_dev_cu(
#ifdef  GADGET_MAC
				 absErr, newton
#else///GADGET_MAC
#     ifdef  WS93_MAC
                                 accErr
#     endif//WS93_MAC
#endif//GADGET_MAC
				 );
#ifdef  SMEM_PREF_FOR_NEIGHBOR_SEARCH
  setGlobalConstants_neighbor_dev_cu();
#endif//SMEM_PREF_FOR_NEIGHBOR_SEARCH
#ifndef SERIALIZED_EXECUTION
  setGlobalConstants_let_dev_cu(
#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
				theta2
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
				);
#endif//SERIALIZED_EXECUTION

#   if  GPUVER >= 70
  setCacheConfig_adv_dev_cu();
#ifndef SERIALIZED_EXECUTION
  setCacheConfig_geo_dev_cu();
#endif//SERIALIZED_EXECUTION
#endif//GPUVER >= 70


  /** memory allocation */
  /** declaration of array to contain whole information of whole N-body particles */
#ifdef  USE_HDF5_FORMAT
  nbody_hdf5 body_snapshot;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
  ulong *hdf5_idx;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  real *hdf5_acc_ext, *hdf5_pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  const muse alloc_snap = allocSnapshotArray
    (&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     &hdf5_acc_ext, &hdf5_pot_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     num_max, &body_snapshot);
#endif//USE_HDF5_FORMAT

  /** declarations of arrays to contain Peano--Hilbert key of N-body particles */
  PHint *peano, *peano_dev;
  int *tag_dev;
  float4 *min_dev, *max_dev;
#ifndef SERIALIZED_EXECUTION
  float4 *box_min, *box_max, *min_hst, *max_hst;
#endif//SERIALIZED_EXECUTION
  int *gsync_ph0, *gsync_ph1;
  soaPHsort soaPH_dev, soaPH_hst;
#ifdef  CUB_AVAILABLE
  soaPHsort soaPH_pre;
  PHint *peano_pre;
  int *tag_pre;
  void *phsort_temp_storage;
#endif//CUB_AVAILABLE
#   if  !defined(SERIALIZED_EXECUTION) && defined(USE_OCCUPANCY_CALCULATOR)
  checkBoxSize_dev(&soaPH_dev);
#endif//!defined(SERIALIZED_EXECUTION) && defined(USE_OCCUPANCY_CALCULATOR)
  const muse alloc_phkey = allocPeanoHilbertKey_dev
    (num_max, &tag_dev, &peano_dev, &peano, &min_dev, &max_dev, &gsync_ph0, &gsync_ph1,
#ifndef SERIALIZED_EXECUTION
    &box_min, &box_max, &min_hst, &max_hst,
#endif//SERIALIZED_EXECUTION
     &soaPH_dev, &soaPH_hst,
#ifdef  CUB_AVAILABLE
     &soaPH_pre, &phsort_temp_storage, &tag_pre, &peano_pre,
#endif//CUB_AVAILABLE
     devProp);

  /** declarations of arrays to contain information of i-particles */
  iparticle     ibody0,  ibody1,  ibody0_dev,  ibody1_dev;
  ulong        *  idx0, *  idx1, *  idx0_dev, *  idx1_dev;
  position     *  pos0, *  pos1, *  pos0_dev, *  pos1_dev;
  acceleration *  acc0, *  acc1, *  acc0_dev, *  acc1_dev;
  real *neighbor_dev;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  acceleration *acc_ext0, *acc_ext1, *acc_ext_dev;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  position *encBall_dev, *encBall_hst;
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
  position *box_min_hst, *box_max_hst, *icom_hst;
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
  DPacc *tmp_dev;
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
  acceleration *res_dev;
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
#ifdef  SWITCH_WITH_J_PARALLELIZATION
#ifdef  MPI_VIA_HOST
  iparticle     ibody_dist0,  ibody_dist1;
  ulong        *  idx_dist0, *  idx_dist1;
  position     *  pos_dist0, *  pos_dist1;
  acceleration *  acc_dist0, *  acc_dist1;
#endif//MPI_VIA_HOST
  iparticle     ibody_dist0_dev,  ibody_dist1_dev;
  ulong        *  idx_dist0_dev, *  idx_dist1_dev;
  position     *  pos_dist0_dev, *  pos_dist1_dev;
  acceleration *  acc_dist0_dev, *  acc_dist1_dev;
  real *neighbor_dist_dev;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  acceleration *acc_ext_dist0, *acc_ext_dist1, *acc_ext_dist_dev;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  position *encBall_dist_dev, *encBall_dist_hst;
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
  position *box_min_dist_hst, *box_max_dist_hst, *icom_dist_hst;
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
  DPacc *tmp_dist_dev;
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
  acceleration *res_dist_dev;
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
#endif//SWITCH_WITH_J_PARALLELIZATION
#ifdef  BLOCK_TIME_STEP
  velocity     * vel0, * vel1, * vel0_dev, * vel1_dev;
  ibody_time   *time0, *time1, *time0_dev, *time1_dev;
  const muse alloc_ibody0    = allocParticleDataSoA_hst
    (num_max, &ibody0, &idx0, &pos0, &acc0, &vel0, &time0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext0
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  const muse alloc_ibody1    = allocParticleDataSoA_hst
    (num_max, &ibody1, &idx1, &pos1, &acc1, &vel1, &time1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext1
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  const muse alloc_ibody_dev = allocParticleDataSoA_dev
    (num_max
     , &ibody0_dev, &idx0_dev, &pos0_dev, &acc0_dev, &vel0_dev, &time0_dev
     , &ibody1_dev, &idx1_dev, &pos1_dev, &acc1_dev, &vel1_dev, &time1_dev
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext_dev
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     , &neighbor_dev
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
     , &encBall_dev, &encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     , &box_min_hst, &box_max_hst, &icom_hst
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
     , &tmp_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     , &res_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     );
#ifdef  SWITCH_WITH_J_PARALLELIZATION
  velocity   * vel_dist0_dev, * vel_dist1_dev;
  ibody_time *time_dist0_dev, *time_dist1_dev;
#ifdef  MPI_VIA_HOST
  velocity     * vel_dist0, * vel_dist1;
  ibody_time   *time_dist0, *time_dist1;
  const muse alloc_ibody_dist0    = allocParticleDataSoA_hst
    (NMAX_J_PARALLELIZATION, &ibody_dist0, &idx_dist0, &pos_dist0, &acc_dist0, &vel_dist0, &time_dist0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext_dist0
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  const muse alloc_ibody_dist1    = allocParticleDataSoA_hst
    (NMAX_J_PARALLELIZATION, &ibody_dist1, &idx_dist1, &pos_dist1, &acc_dist1, &vel_dist1, &time_dist1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext_dist1
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
#endif//MPI_VIA_HOST
  const muse alloc_ibody_dist_dev = allocParticleDataSoA_dev
    (NMAX_J_PARALLELIZATION
     , &ibody_dist0_dev, &idx_dist0_dev, &pos_dist0_dev, &acc_dist0_dev, &vel_dist0_dev, &time_dist0_dev
     , &ibody_dist1_dev, &idx_dist1_dev, &pos_dist1_dev, &acc_dist1_dev, &vel_dist1_dev, &time_dist1_dev
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext_dist_dev
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     , &neighbor_dist_dev
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
     , &encBall_dist_dev, &encBall_dist_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     , &box_min_dist_hst, &box_max_dist_hst, &icom_dist_hst
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
     , &tmp_dist_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     , &res_dist_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     );
#endif//SWITCH_WITH_J_PARALLELIZATION
#else///BLOCK_TIME_STEP
  real *vx0, *vx1, *vx0_dev, *vx1_dev;
  real *vy0, *vy1, *vy0_dev, *vy1_dev;
  real *vz0, *vz1, *vz0_dev, *vz1_dev;
  const muse alloc_ibody0    = allocParticleDataSoA_hst
    (num_max, &ibody0, &idx0, &pos0, &acc0, &vx0, &vy0, &vz0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext0
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  const muse alloc_ibody1    = allocParticleDataSoA_hst
    (num_max, &ibody1, &idx1, &pos1, &acc1, &vx1, &vy1, &vz1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext1
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  const muse alloc_ibody_dev = allocParticleDataSoA_dev
    (num_max
     , &ibody0_dev, &idx0_dev, &pos0_dev, &acc0_dev, &vx0_dev, &vy0_dev, &vz0_dev
     , &ibody1_dev, &idx1_dev, &pos1_dev, &acc1_dev, &vx1_dev, &vy1_dev, &vz1_dev
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext_dev
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     , &neighbor_dev
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
     , &encBall_dev, &encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     , &box_min_hst, &box_max_hst, &icom_hst
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
     , &tmp_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     , &res_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     );
#ifdef  SWITCH_WITH_J_PARALLELIZATION
  real *vx_dist0_dev, *vx_dist1_dev;
  real *vy_dist0_dev, *vy_dist1_dev;
  real *vz_dist0_dev, *vz_dist1_dev;
#ifdef  MPI_VIA_HOST
  real *vx_dist0, *vx_dist1;
  real *vy_dist0, *vy_dist1;
  real *vz_dist0, *vz_dist1;
  const muse alloc_ibody_dist0    = allocParticleDataSoA_hst
    (NMAX_J_PARALLELIZATION, &ibody_dist0, &idx_dist0, &pos_dist0, &acc_dist0, &vx_dist0, &vy_dist0, &vz_dist0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext_dist0
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  const muse alloc_ibody_dist1    = allocParticleDataSoA_hst
    (NMAX_J_PARALLELIZATION, &ibody_dist1, &idx_dist1, &pos_dist1, &acc_dist1, &vx_dist1, &vy_dist1, &vz_dist1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext_dist1
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
#endif//MPI_VIA_HOST
  const muse alloc_ibody_dist_dev = allocParticleDataSoA_dev
    (NMAX_J_PARALLELIZATION
     , &ibody_dist0_dev, &idx_dist0_dev, &pos_dist0_dev, &acc_dist0_dev, &vx_dist0_dev, &vy_dist0_dev, &vz_dist0_dev
     , &ibody_dist1_dev, &idx_dist1_dev, &pos_dist1_dev, &acc_dist1_dev, &vx_dist1_dev, &vy_dist1_dev, &vz_dist1_dev
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext_dist_dev
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     , &neighbor_dist_dev
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
     , &encBall_dist_dev, &encBall_dist_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     , &box_min_dist_hst, &box_max_dist_hst, &icom_dist_hst
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
     , &tmp_dist_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     , &res_dist_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     );
#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//BLOCK_TIME_STEP

  iparticle_treeinfo treeinfo_dev;
  int *jtag_dev;
#ifdef  COUNT_INTERACTIONS
  iparticle_treeinfo treeinfo;
  int *Nj, *Nj_dev, *Nbuf, *Nbuf_dev;
#endif//COUNT_INTERACTIONS
  const muse alloc_jtag = allocParticleInfoSoA_hst
    (num_max
#ifdef  COUNT_INTERACTIONS
     , &treeinfo, &Nj, &Nbuf
#endif//COUNT_INTERACTIONS
     );
  const muse alloc_jtag_dev = allocParticleInfoSoA_dev
    (num_max, &treeinfo_dev, &jtag_dev
#ifdef  COUNT_INTERACTIONS
     , &Nj_dev, &Nbuf_dev
#endif//COUNT_INTERACTIONS
     );

  /** declarations of variables to contain information on groups of i-particles */
  laneinfo *laneInfo_hst, *laneInfo_dev;
  double *laneTime_dev;
  int inumPerLane, maxNgrp;
  int *inum_hst, *inum_dev;
  const muse alloc_lane_dev =
    allocParticleGroups(&laneInfo_hst, &laneInfo_dev, &laneTime_dev, &inum_hst, &inum_dev,
#ifdef  SWITCH_WITH_J_PARALLELIZATION
			true,
#endif//SWITCH_WITH_J_PARALLELIZATION
			&inumPerLane, &maxNgrp, num_max);
#ifdef  SWITCH_WITH_J_PARALLELIZATION
  int inumPerLane_ext, maxNgrp_ext;
  laneinfo *laneInfo_ext_hst, *laneInfo_ext_dev;
  const muse alloc_lane_ext_dev =
    allocParticleGroups(&laneInfo_ext_hst, &laneInfo_ext_dev, NULL, NULL, NULL, false, &inumPerLane_ext, &maxNgrp_ext, NMAX_J_PARALLELIZATION);
#endif//SWITCH_WITH_J_PARALLELIZATION

  /** declarations of arrays to contain information of tree-cells */
  PHint    *hkey_dev;
  uint     *parent_dev, *children_dev;
  soaTreeCell soaCell_dev;
  treecell *cell_dev;
  bool     *leaf_dev;
  uint     *node_dev;
  PHinfo   *list_dev;
#ifdef  COUNT_INTERACTIONS
  soaTreeCell soaCell_hst;
  treecell *cell;
  bool     *leaf;
  uint     *node;
  PHinfo   *list;
#endif//COUNT_INTERACTIONS
  int *bottomLev_dev, *numCell_dev, *numNode_dev, *scanNum_dev;
  const muse alloc_cell_dev =
    allocTreeCell_dev(&soaCell_dev, &cell_dev, &leaf_dev, &node_dev, &list_dev,
		      &hkey_dev, &parent_dev, &children_dev, &bottomLev_dev, &numCell_dev, &numNode_dev, &scanNum_dev
#ifdef  COUNT_INTERACTIONS
		      , &soaCell_hst, &cell    , &leaf    , &node    , &list
#endif//COUNT_INTERACTIONS
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
		      , devProp
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
		      );

  /** declarations of arrays to contain information of tree-nodes */
  soaTreeNode soaNode_dev;
  uint    *more_dev;
  jparticle *pj_dev;
  jmass     *mj_dev;
  real *bmax_dev;
  int *node2cell_dev;
  int *gsync0, *gsync1;
#ifdef  WS93_MAC
  real *mr2_dev;
#endif//WS93_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
  soaTreeNode soaNode_hst;
  uint *more;
  jparticle *pj;
  jmass     *mj;
#endif//!defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
  int *gmem_make_tree_dev, *gsync0_make_tree_dev, *gsync1_make_tree_dev, *gsync2_make_tree_dev, *gsync3_make_tree_dev;
  int *gmem_link_tree_dev, *gsync0_link_tree_dev, *gsync1_link_tree_dev;
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
  soaMakeTreeBuf soaMakeBuf;

/* #ifdef  WS93_MAC */
/*   __NOTE__("WS93_MAC is defined\n"); */
/* #endif//WS93_MAC */
/* #ifdef  GADGET_MAC */
/*   __NOTE__("GADGET_MAC is defined\n"); */
/* #endif//GADGET_MAC */
/* #ifdef  SERIALIZED_EXECUTION */
/*   __NOTE__("SERIALIZED_EXECUTION is defined\n"); */
/* #endif//SERIALIZED_EXECUTION */
/* #ifdef  MPI_VIA_HOST */
/*   __NOTE__("MPI_VIA_HOST is defined\n"); */
/* #endif//MPI_VIA_HOST */
/* #ifdef  CARE_EXTERNAL_PARTICLES */
/*   __NOTE__("CARE_EXTERNAL_PARTICLES is defined\n"); */
/* #endif//CARE_EXTERNAL_PARTICLES */
  /* __NOTE__("address of dev is %p\n", &soaNode_dev); */
  /* __NOTE__("address of more_dev is %p\n", &more_dev); */
  /* __NOTE__("address of pj_dev is %p\n", &pj_dev); */
  /* __NOTE__("address of mj_dev is %p\n", &mj_dev); */
  /* __NOTE__("address of bmax_dev is %p\n", &bmax_dev); */
  /* __NOTE__("address of n2c_dev is %p\n", &node2cell_dev); */
  /* __NOTE__("address of gsync0 is %p\n", &gsync0); */
  /* __NOTE__("address of gsync1 is %p\n", &gsync1); */
  /* __NOTE__("address of gmem_make_tree is %p\n", &gmem_make_tree_dev); */
  /* __NOTE__("address of gsync0_make_tree is %p\n", &gsync0_make_tree_dev); */
  /* __NOTE__("address of gsync1_make_tree is %p\n", &gsync1_make_tree_dev); */
  /* __NOTE__("address of gsync2_make_tree is %p\n", &gsync2_make_tree_dev); */
  /* __NOTE__("address of gsync3_make_tree is %p\n", &gsync3_make_tree_dev); */
  /* __NOTE__("address of gmem_link_tree is %p\n", &gmem_link_tree_dev); */
  /* __NOTE__("address of gsync0_link_tree is %p\n", &gsync0_link_tree_dev); */
  /* __NOTE__("address of gsync1_link_tree is %p\n", &gsync1_link_tree_dev); */
  /* __NOTE__("address of mac_dev is %p\n", &mac_dev); */
  /* __NOTE__("address of more0Buf is %p\n", &more0Buf); */
  /* __NOTE__("address of more1Buf is %p\n", &more1Buf); */
  /* __NOTE__("address of rjmaxBuf is %p\n", &rjmaxBuf); */
  /* __NOTE__("address of fail_dev is %p\n", &makeFail); */
  /* __NOTE__("address of buf is %p\n", &soaMakeBuf); */
  const muse alloc_node_dev = allocTreeNode_dev
    (&soaNode_dev, &more_dev, &pj_dev, &mj_dev, &bmax_dev, &node2cell_dev, &gsync0, &gsync1,
#ifdef  WS93_MAC
     &mr2_dev,
#endif//WS93_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
     &soaNode_hst, &more, &pj, &mj,
#endif//!defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
     &gmem_make_tree_dev, &gsync0_make_tree_dev, &gsync1_make_tree_dev, &gsync2_make_tree_dev, &gsync3_make_tree_dev,
     &gmem_link_tree_dev, &gsync0_link_tree_dev, &gsync1_link_tree_dev,
#ifdef  GADGET_MAC
     &mac_dev,
#endif//GADGET_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
     &gmem_external_dev, &gsync0_external_dev, &gsync1_external_dev, &diameter_dev, &diameter_hst, &location, CAST_R2F(eps), CAST_R2F(eta),
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
     &more0Buf, &more1Buf, &rjmaxBuf, &makeFail, &soaMakeBuf, devProp);
  soaNode_dev.jtag = jtag_dev;


#ifndef SERIALIZED_EXECUTION
  sendCfg *iparticleSendBuf;
  recvCfg *iparticleRecvBuf;
  int *sampleRecvNum, *sampleRecvDsp;

  static MPIinfo ormCfg[3], repCfg[3];
  float *dxmin, *dxmax, *dymin, *dymax, *dzmin, *dzmax;
  float *sxmin, *sxmax, *symin, *symax, *szmin, *szmax;
  MPI_Request *dmreq;
  domainCfg domCfg;
  sampling sample;
  const muse alloc_dd = allocateORMtopology
    (&dxmin, &dxmax, &dymin, &dymax, &dzmin, &dzmax, &dmreq,
     &sxmin, &sxmax, &symin, &symax, &szmin, &szmax,
     &iparticleSendBuf, &iparticleRecvBuf, &sampleRecvNum, &sampleRecvDsp,
     ormCfg, repCfg, Nx, Ny, Nz, &letcfg, &domCfg, &sample, Ntot);

  /* float4 *pmin_hst, *pmax_hst, *pmin_dev, *pmax_dev; */
  /* int *gsync_box0, *gsync_box1; */
  /* soaBoxSize soaBox; */
  /* const muse alloc_box = allocateBoxSize_dev(&pmin_hst, &pmax_hst, &pmin_dev, &pmax_dev, &gsync_box0, &gsync_box1, &soaBox, devProp); */
#ifndef USE_OCCUPANCY_CALCULATOR
  checkBoxSize_dev(devProp);
#endif//USE_OCCUPANCY_CALCULATOR
  float *x0hst, *x1hst, *y0hst, *y1hst, *z0hst, *z1hst;  int *idhst;
  float *x0dev, *x1dev, *y0dev, *y1dev, *z0dev, *z1dev;  int *iddev;
  samplePos samplePos0, samplePos1;
  const muse alloc_spl = allocateSamplePos
    (&x0hst, &x1hst, &y0hst, &y1hst, &z0hst, &z1hst, &idhst,
     &x0dev, &x1dev, &y0dev, &y1dev, &z0dev, &z1dev, &iddev, &samplePos0, &samplePos1, sample);
  float *xhst, *yhst, *zhst;  particlePos particlePos_hst;
  float *xdev, *ydev, *zdev;  particlePos particlePos_dev;
  int *rank_hst, *rank_dev, *idx_dev;
  domainDecomposeKey domDecKey;
  const muse alloc_pos = allocateParticlePosition
    (&xhst, &yhst, &zhst, &particlePos_hst, &xdev, &ydev, &zdev, &particlePos_dev, &rank_hst, &rank_dev, &idx_dev, &domDecKey, num_max);

  float *xmin_dev, *xmax_dev, *ymin_dev, *ymax_dev, *zmin_dev, *zmax_dev;
  float *xmin_hst, *xmax_hst, *ymin_hst, *ymax_hst, *zmin_hst, *zmax_hst;
  int *domrank_dev, *domrank_hst;
  int *numNew_dev, *numNew_hst;
  int4 *gmem_dom;
  int *gsync0_dom, *gsync1_dom;
  sendDom domBoundary;
  const muse alloc_dom = allocateDomainPos
    (&xmin_dev, &xmax_dev, &ymin_dev, &ymax_dev, &zmin_dev, &zmax_dev,
     &xmin_hst, &xmax_hst, &ymin_hst, &ymax_hst, &zmin_hst, &zmax_hst,
     &domrank_dev, &domrank_hst, &numNew_dev, &numNew_hst, &gmem_dom, &gsync0_dom, &gsync1_dom,
     &domBoundary, letcfg.size, devProp);

  soaGEO soaGEO_dev;
  real *r2geo_dev;
  const muse alloc_r2geo = allocGeometricEnclosingBall_dev(&r2geo_dev, &soaGEO_dev, num_max);
  int *numSend_hst, *numSend_dev;
  static domainInfo *nodeInfo;
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
  static position *min, *max;
#endif//USE_RECTANGULAR_BOX_FOR_LET
  static position   *ipos;
#ifdef  GADGET_MAC
  static real *amin;
#endif//GADGET_MAC
  int Nstream_let;
  cudaStream_t *stream_let;
  const muse alloc_LETtopology = configLETtopology
    (&nodeInfo,
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     &min, &max,
#endif//USE_RECTANGULAR_BOX_FOR_LET
     &ipos,
#ifdef  GADGET_MAC
     &amin,
#endif//GADGET_MAC
     &numSend_hst, &numSend_dev, &stream_let, &Nstream_let, devProp, letcfg);

#   if  defined(USE_ENCLOSING_BALL_FOR_LET) && !defined(OCTREE_BASED_SEARCH)
  soaEncBall soaEB;
  void *appEncBall_dev;
  int *gsync0_EB, *gsync1_EB;
  const muse alloc_appEnc = allocApproxEnclosingBall_dev(&appEncBall_dev, &gsync0_EB, &gsync1_EB, &soaEB, devProp);
#endif//defined(USE_ENCLOSING_BALL_FOR_LET) && !defined(OCTREE_BASED_SEARCH)

#ifdef  SWITCH_WITH_J_PARALLELIZATION
  int *Ni_list, *head_list, *grpNum_list, *grpNum_disp;
  const muse alloc_extBodyInfo = allocateExternalBodyInfo(&Ni_list, &head_list, &grpNum_list, &grpNum_disp, letcfg);
  static bool transferMode = false;
  static int Ni_local, Ni_total;
#endif//SWITCH_WITH_J_PARALLELIZATION

#ifdef  MPI_ONE_SIDED_FOR_LET
  /** create MPI Window objects for one-sided communication */
  const size_t let_win_size = (size_t)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE);/**< corresponding to # of elements of tree nodes allocated in allocTreeNode_dev() in src/tree/make_dev.cu */
#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_create(soaNode_dev.more, sizeof(uint)      * let_win_size, sizeof(uint)	    , letcfg.info, letcfg.comm, &(letcfg.win_more)));
  chkMPIerr(MPI_Win_create(soaNode_dev.jpos, sizeof(jparticle) * let_win_size, sizeof(jparticle), letcfg.info, letcfg.comm, &(letcfg.win_jpos)));
  chkMPIerr(MPI_Win_create(soaNode_dev.mj  , sizeof(jmass)     * let_win_size, sizeof(jmass)    , letcfg.info, letcfg.comm, &(letcfg.win_mass)));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_create(soaNode_hst.more, sizeof(uint)      * let_win_size, sizeof(uint)	    , letcfg.info, letcfg.comm, &(letcfg.win_more)));
  chkMPIerr(MPI_Win_create(soaNode_hst.jpos, sizeof(jparticle) * let_win_size, sizeof(jparticle), letcfg.info, letcfg.comm, &(letcfg.win_jpos)));
  chkMPIerr(MPI_Win_create(soaNode_hst.mj  , sizeof(jmass)     * let_win_size, sizeof(jmass)    , letcfg.info, letcfg.comm, &(letcfg.win_mass)));
#endif//MPI_VIA_HOST
#endif//MPI_ONE_SIDED_FOR_LET

#ifdef  MPI_ONE_SIDED_FOR_EXCG

#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_create(samplePos0.x_dev, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos0.win_x)));
  chkMPIerr(MPI_Win_create(samplePos0.y_dev, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos0.win_y)));
  chkMPIerr(MPI_Win_create(samplePos0.z_dev, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos0.win_z)));
  chkMPIerr(MPI_Win_create(samplePos1.x_dev, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos1.win_x)));
  chkMPIerr(MPI_Win_create(samplePos1.y_dev, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos1.win_y)));
  chkMPIerr(MPI_Win_create(samplePos1.z_dev, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos1.win_z)));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_create(samplePos0.x_hst, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos0.win_x)));
  chkMPIerr(MPI_Win_create(samplePos0.y_hst, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos0.win_y)));
  chkMPIerr(MPI_Win_create(samplePos0.z_hst, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos0.win_z)));
  chkMPIerr(MPI_Win_create(samplePos1.x_hst, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos1.win_x)));
  chkMPIerr(MPI_Win_create(samplePos1.y_hst, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos1.win_y)));
  chkMPIerr(MPI_Win_create(samplePos1.z_hst, sizeof(float) * sample.Nmax, sizeof(float), letcfg.info, letcfg.comm, &(samplePos1.win_z)));
#endif//MPI_VIA_HOST

  const size_t body_win_size = (size_t)num_max + (size_t)(((num_max % NTHREADS) != 0) ? (NTHREADS - (num_max & NTHREADS)) : (0));/**< corresponding to # of particles per node allocated in allocParticleDataSoA_dev() in src/misc/allocate_dev.cu */

#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_create(ibody0_dev.pos, sizeof(position) * body_win_size, sizeof(position), letcfg.info, letcfg.comm, &(ibody0_dev.win_ipos)));
  chkMPIerr(MPI_Win_create(ibody1_dev.pos, sizeof(position) * body_win_size, sizeof(position), letcfg.info, letcfg.comm, &(ibody1_dev.win_ipos)));
#ifdef  GADGET_MAC
  chkMPIerr(MPI_Win_create(ibody0_dev.acc, sizeof(acceleration) * body_win_size, sizeof(acceleration), letcfg.info, letcfg.comm, &(ibody0_dev.win_iacc)));
  chkMPIerr(MPI_Win_create(ibody1_dev.acc, sizeof(acceleration) * body_win_size, sizeof(acceleration), letcfg.info, letcfg.comm, &(ibody1_dev.win_iacc)));
#endif//GADGET_MAC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  chkMPIerr(MPI_Win_create(ibody0_dev.ext, sizeof(acceleration) * body_win_size, sizeof(acceleration), letcfg.info, letcfg.comm, &(ibody0_dev.win_iext)));
  chkMPIerr(MPI_Win_create(ibody1_dev.ext, sizeof(acceleration) * body_win_size, sizeof(acceleration), letcfg.info, letcfg.comm, &(ibody1_dev.win_iext)));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_create(ibody0_dev.vel, sizeof(velocity) * body_win_size, sizeof(position), letcfg.info, letcfg.comm, &(ibody0_dev.win_ivel)));
  chkMPIerr(MPI_Win_create(ibody1_dev.vel, sizeof(velocity) * body_win_size, sizeof(position), letcfg.info, letcfg.comm, &(ibody1_dev.win_ivel)));
  chkMPIerr(MPI_Win_create(ibody0_dev.time, sizeof(ibody_time) * body_win_size, sizeof(ibody_time), letcfg.info, letcfg.comm, &(ibody0_dev.win_time)));
  chkMPIerr(MPI_Win_create(ibody1_dev.time, sizeof(ibody_time) * body_win_size, sizeof(ibody_time), letcfg.info, letcfg.comm, &(ibody1_dev.win_time)));
#else///BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_create(ibody0_dev.vx, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody0_dev.win_vx)));
  chkMPIerr(MPI_Win_create(ibody1_dev.vx, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody1_dev.win_vx)));
  chkMPIerr(MPI_Win_create(ibody0_dev.vy, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody0_dev.win_vy)));
  chkMPIerr(MPI_Win_create(ibody1_dev.vy, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody1_dev.win_vy)));
  chkMPIerr(MPI_Win_create(ibody0_dev.vz, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody0_dev.win_vz)));
  chkMPIerr(MPI_Win_create(ibody1_dev.vz, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody1_dev.win_vz)));
#endif//BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_create(ibody0_dev.idx, sizeof(ulong) * body_win_size, sizeof(ulong), letcfg.info, letcfg.comm, &(ibody0_dev.win_idx)));
  chkMPIerr(MPI_Win_create(ibody1_dev.idx, sizeof(ulong) * body_win_size, sizeof(ulong), letcfg.info, letcfg.comm, &(ibody1_dev.win_idx)));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_create(ibody0.pos, sizeof(position) * body_win_size, sizeof(position), letcfg.info, letcfg.comm, &(ibody0.win_ipos)));
  chkMPIerr(MPI_Win_create(ibody1.pos, sizeof(position) * body_win_size, sizeof(position), letcfg.info, letcfg.comm, &(ibody1.win_ipos)));
#ifdef  GADGET_MAC
  chkMPIerr(MPI_Win_create(ibody0.acc, sizeof(acceleration) * body_win_size, sizeof(acceleration), letcfg.info, letcfg.comm, &(ibody0.win_iacc)));
  chkMPIerr(MPI_Win_create(ibody1.acc, sizeof(acceleration) * body_win_size, sizeof(acceleration), letcfg.info, letcfg.comm, &(ibody1.win_iacc)));
#endif//GADGET_MAC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  chkMPIerr(MPI_Win_create(ibody0.ext, sizeof(acceleration) * body_win_size, sizeof(acceleration), letcfg.info, letcfg.comm, &(ibody0.win_iext)));
  chkMPIerr(MPI_Win_create(ibody1.ext, sizeof(acceleration) * body_win_size, sizeof(acceleration), letcfg.info, letcfg.comm, &(ibody1.win_iext)));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_create(ibody0.vel, sizeof(velocity) * body_win_size, sizeof(position), letcfg.info, letcfg.comm, &(ibody0.win_ivel)));
  chkMPIerr(MPI_Win_create(ibody1.vel, sizeof(velocity) * body_win_size, sizeof(position), letcfg.info, letcfg.comm, &(ibody1.win_ivel)));
  chkMPIerr(MPI_Win_create(ibody0.time, sizeof(ibody_time) * body_win_size, sizeof(ibody_time), letcfg.info, letcfg.comm, &(ibody0.win_time)));
  chkMPIerr(MPI_Win_create(ibody1.time, sizeof(ibody_time) * body_win_size, sizeof(ibody_time), letcfg.info, letcfg.comm, &(ibody1.win_time)));
#else///BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_create(ibody0.vx, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody0.win_vx)));
  chkMPIerr(MPI_Win_create(ibody1.vx, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody1.win_vx)));
  chkMPIerr(MPI_Win_create(ibody0.vy, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody0.win_vy)));
  chkMPIerr(MPI_Win_create(ibody1.vy, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody1.win_vy)));
  chkMPIerr(MPI_Win_create(ibody0.vz, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody0.win_vz)));
  chkMPIerr(MPI_Win_create(ibody1.vz, sizeof(real) * body_win_size, sizeof(real), letcfg.info, letcfg.comm, &(ibody1.win_vz)));
#endif//BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_create(ibody0.idx, sizeof(ulong) * body_win_size, sizeof(ulong), letcfg.info, letcfg.comm, &(ibody0.win_idx)));
  chkMPIerr(MPI_Win_create(ibody1.idx, sizeof(ulong) * body_win_size, sizeof(ulong), letcfg.info, letcfg.comm, &(ibody1.win_idx)));
#endif//MPI_VIA_HOST
#endif//MPI_ONE_SIDED_FOR_EXCG

#endif//SERIALIZED_EXECUTION


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
#endif//COMPARE_WITH_DIRECT_SOLVER


  /** initialize the simulation run */
  /** variables for automatic optimization */
  static struct timespec start;
  static measuredTime elapsed;
  static autoTuningParam rebuildParam;
  initStatVal(&rebuildParam.linearStats);  initGuessTime(&rebuildParam.linearGuess);
  initStatVal(&rebuildParam. powerStats);  initGuessTime(&rebuildParam. powerGuess);
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  initStatVal  (&rebuildParam.parabolicStats);
  initGuessTime(&rebuildParam.parabolicGuess);
#ifdef  USE_ADDITIONAL_SWITCH
#ifndef DISABLE_AUTO_TUNING
  static int useParabolicGuess = 0;
#endif//DISABLE_AUTO_TUNING
#endif//USE_ADDITIONAL_SWITCH
#endif//USE_PARABOLIC_GROWTH_MODEL

  /** variables to control tree rebuild interval */
  static rebuildTree rebuild;
  rebuild.reuse = 1;
  rebuild.interval = 0.0;
#ifdef  BLOCK_TIME_STEP
  rebuild.adjust = false;
#ifdef  FORCE_ADJUSTING_PARTICLE_TIME_STEPS
  rebuild.avg = 0.0;
  rebuild.var = 0.0;
#endif//FORCE_ADJUSTING_PARTICLE_TIME_STEPS
#endif//BLOCK_TIME_STEP
  static brentStatus brentDistance;  brentDistance.initialized = false;
  static brentMemory brentHistory = {0.0, 0.0, 0, 0};
#ifdef  SHARED_AUTO_TUNER
  static double shared_brent_u;
  static int shared_brent_num;
#endif//SHARED_AUTO_TUNER

  /** read initial condition */
#   if  defined(MONITOR_ENERGY_ERROR) && defined(USE_HDF5_FORMAT)
  static energyError relEneErr;
#endif//defined(MONITOR_ENERGY_ERROR) && defined(USE_HDF5_FORMAT)
  static double time, dt;
  static ulong steps = 0;
  static double prevElapsed = 0.0;
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
#endif//USE_HDF5_FORMAT
#ifdef  SERIALIZED_EXECUTION
  readTentativeData
    (&time, &dt, &steps, &prevElapsed, num, ibody0, file, last
#ifdef  USE_HDF5_FORMAT
     , hdf5type, &dropPrevTune, &rebuild, &elapsed, &rebuildParam, &brentDistance, &brentHistory
#ifdef  MONITOR_ENERGY_ERROR
     , &relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
     );
#else///SERIALIZED_EXECUTION
  readTentativeDataParallel
    (&time, &dt, &steps, &prevElapsed, &num, ibody0, file, last, &iocfg, Ntot
#ifdef  USE_HDF5_FORMAT
     , hdf5type, &dropPrevTune, &rebuild, &elapsed, &rebuildParam, &brentDistance, &brentHistory
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
#ifdef  REPORT_COMPUTE_RATE
  static struct timespec clock_prev;
  clock_gettime(CLOCK_MONOTONIC_RAW, &clock_prev);
#ifndef EXEC_BENCHMARK
  static double time_prev;  time_prev = time;
  static ulong  step_prev;  step_prev = steps;
#endif//EXEC_BENCHMARK
  static double brent_rate = 0.0;/**< brent->u.pos / brent->b */
  static ulong brent_prev;/**< step of the previous update of brentStatus */
  brent_prev = steps;
  static ulong rebuild_freq = 0;/**< number of tree rebuild */
#endif//REPORT_COMPUTE_RATE

  /** set up for output files */
  static double formerTime;
#ifdef  SERIALIZED_EXECUTION
  formerTime = getPresentTimeInMin();
#else///SERIALIZED_EXECUTION
  if( mpi.rank == 0 )
    formerTime = getPresentTimeInMin();
  chkMPIerr(MPI_Bcast(&formerTime, 1, MPI_DOUBLE, 0, mpi.comm));
#endif//SERIALIZED_EXECUTION
#   if  defined(BLOCK_TIME_STEP) || !defined(EXEC_BENCHMARK)
  static double invSnapshotInterval;
  invSnapshotInterval = 1.0 / snapshotInterval;
  static uint previous;
  previous = (uint)(time * invSnapshotInterval);
#endif//defined(BLOCK_TIME_STEP) || !defined(EXEC_BENCHMARK)


  /** read numeric table for fixed potential field */
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  potential_field pot_tbl_sphe;
  pot2 *pot_tbl_sphe_Phi;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  real *pot_tbl_sphe_rad;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  const muse alloc_ext_pot_sphe = readFixedPotentialTableSpherical
    (unit, pot_file_sphe, &pot_tbl_sphe, &pot_tbl_sphe_Phi
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     , &pot_tbl_sphe_rad
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#ifdef  USE_HDF5_FORMAT
     , hdf5type
#endif//USE_HDF5_FORMAT
     );
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
  disk_potential pot_tbl_disk;
  real *pot_tbl_disk_Phi;
  pot2 *pot_tbl_disk_sphe_Phi;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  real *pot_tbl_disk_RR, *pot_tbl_disk_zz;
  real *pot_tbl_disk_sphe_rad;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  disk_grav *pot_tbl_disk_FRz;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  const muse alloc_ext_pot_disk = readFixedPotentialTableDisk
    (unit, pot_file_disk, &pot_tbl_disk_Phi, &pot_tbl_disk_sphe_Phi,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     &pot_tbl_disk_RR, &pot_tbl_disk_zz, &pot_tbl_disk_sphe_rad,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     &pot_tbl_disk_FRz,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     &pot_tbl_disk
#ifdef  USE_HDF5_FORMAT
   , hdf5type
#endif//USE_HDF5_FORMAT
     );
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD


  /** check size of memory in use */
  muse       body_mem = {alloc_phkey.host   + alloc_ibody0.host   + alloc_ibody1.host   + alloc_ibody_dev.host,
			 alloc_phkey.device + alloc_ibody0.device + alloc_ibody1.device + alloc_ibody_dev.device};
  const muse tree_mem = {alloc_jtag.host   + alloc_jtag_dev.host   + alloc_cell_dev.host   + alloc_node_dev.host,
			 alloc_jtag.device + alloc_jtag_dev.device + alloc_cell_dev.device + alloc_node_dev.device};
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
  const muse box_mem = {alloc_spl.host   + alloc_pos.host   + alloc_dom.host  ,
			alloc_spl.device + alloc_pos.device + alloc_dom.device};
  used_mem.host   += box_mem.host;
  used_mem.device += box_mem.device;
#   if  defined(USE_ENCLOSING_BALL_FOR_LET) && !defined(OCTREE_BASED_SEARCH)
  const muse ball_mem = {alloc_r2geo.host + alloc_appEnc.host, alloc_r2geo.device + alloc_appEnc.device};
#else///defined(USE_ENCLOSING_BALL_FOR_LET) && !defined(OCTREE_BASED_SEARCH)
  const muse ball_mem = {alloc_r2geo.host, alloc_r2geo.device};
#endif//defined(USE_ENCLOSING_BALL_FOR_LET) && !defined(OCTREE_BASED_SEARCH)
  used_mem.host   += ball_mem.host;
  used_mem.device += ball_mem.device;
#ifdef  SWITCH_WITH_J_PARALLELIZATION
#ifndef MPI_VIA_HOST
  const muse  ext_mem = {alloc_extBodyInfo.host   + alloc_ibody_dist_dev.host   + alloc_lane_ext_dev.host  ,
			 alloc_extBodyInfo.device + alloc_ibody_dist_dev.device + alloc_lane_ext_dev.device};
#else///MPI_VIA_HOST
  const muse  ext_mem = {alloc_extBodyInfo.host   + alloc_ibody_dist_dev.host   + alloc_lane_ext_dev.host   + alloc_ibody_dist0.host   + alloc_ibody_dist1.host  ,
			 alloc_extBodyInfo.device + alloc_ibody_dist_dev.device + alloc_lane_ext_dev.device + alloc_ibody_dist0.device + alloc_ibody_dist1.device};
#endif//MPI_VIA_HOST
  used_mem.host   += ext_mem.host;
  used_mem.device += ext_mem.device;
#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//SERIALIZED_EXECUTION
#ifdef  COMPARE_WITH_DIRECT_SOLVER
  used_mem.host   += alloc_acc_dev.host;
  used_mem.device += alloc_acc_dev.device;
#endif//COMPARE_WITH_DIRECT_SOLVER
#ifndef BLOCK_TIME_STEP
  used_mem.host   += alloc_dt_dev.host  ;
  used_mem.device += alloc_dt_dev.device;
#endif//BLOCK_TIME_STEP
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  used_mem.host   += alloc_ext_pot_sphe.host  ;
  used_mem.device += alloc_ext_pot_sphe.device;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
  used_mem.host   += alloc_ext_pot_disk.host  ;
  used_mem.device += alloc_ext_pot_disk.device;
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  /** declarations of arrays to prevent for stack overflow */
  int *fail_dev;
  uint *buffer, *freeLst;
  soaTreeWalkBuf soaWalk_dev;
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
  uint *freeNum;
  int *active;
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  unsigned long long int *cycles_hst, *cycles_dev;
#ifndef SERIALIZED_EXECUTION
  unsigned long long int *cycles_dist_hst, *cycles_dist_dev;
#endif//SERIALIZED_EXECUTION
#endif//defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  unsigned long long int *cycles_let_hst, *cycles_let_dev;
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  const muse alloc_buf_dev = allocTreeBuffer_dev
    (&fail_dev, &buffer, &freeLst,
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
     &freeNum, &active,
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
     &cycles_hst, &cycles_dev,
#ifndef SERIALIZED_EXECUTION
     &cycles_dist_hst, &cycles_dist_dev,
#endif//SERIALIZED_EXECUTION
#endif//defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
     &cycles_let_hst, &cycles_let_dev,
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
     &soaWalk_dev, num_max, devProp);
  used_mem.host   += alloc_buf_dev.host;
  used_mem.device += alloc_buf_dev.device;

#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && !defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) && !defined(TIME_BASED_MODIFICATION)
  soaMakeBuf.ubuf_external = buffer;
  soaMakeBuf.Nbuf_external = soaWalk_dev.bufSize * (size_t)(NGROUPS * NBLOCKS_PER_SM * devProp.numSM);
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && !defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) && !defined(TIME_BASED_MODIFICATION)

#ifdef  OUTPUT_MEMORY_USAGE
  /** output memory usage */
  {
    size_t free_mem, used, total;
    queryFreeDeviceMemory(&free_mem, &total);      used = total - free_mem;

    fprintf(stdout, "device ID: %d (%s on %s)\n\n", devIdx, devProp.name, devProp.host);

    fprintf(stdout, "Allocated memory on device per process: %zu B (%zu GiB, %zu MiB)\n", used_mem.device, used_mem.device >> 30, used_mem.device >> 20);
    fprintf(stdout, "    total memory on device per process: %zu B (%zu GiB, %zu MiB)\n"       , total, total >> 30, total >> 20);
    fprintf(stdout, "     used memory on device per process: %zu B (%zu GiB, %zu MiB; %lf%%)\n",  used,  used >> 30,  used >> 20, 100.0 * (double)used / (double)total);
    fprintf(stdout, "     free memory on device per process: %zu B (%zu GiB, %zu MiB; %lf%%)\n",  free_mem,  free_mem >> 30,  free_mem >> 20, 100.0 * (double)free_mem / (double)total);

    fprintf(stdout, "Allocated memory on   host per process: %zu B (%zu GiB, %zu MiB)\n", used_mem.host  , used_mem.host   >> 30, used_mem.host   >> 20);

    fprintf(stdout, "\nBreakdown of allocated memory on device:\n");
    fprintf(stdout, "\t%zu MiB (%zu GiB) for N-body particles\n"  , body_mem.device >> 20, body_mem.device >> 30);
    fprintf(stdout, "\t%zu MiB (%zu GiB) for tree structure\n"    , tree_mem.device >> 20, tree_mem.device >> 30);
    fprintf(stdout, "\t%zu MiB for miscellaneous data\n", misc_mem.device >> 20);
    fprintf(stdout, "\t%zu MiB (%zu GiB) for tree walk buffer\n"  , alloc_buf_dev.device >> 20, alloc_buf_dev.device >> 30);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
    fprintf(stdout, "\t%zu MiB (%zu KiB) for external spherical potential-field\n", alloc_ext_pot_sphe.device >> 20, alloc_ext_pot_sphe.device >> 10);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
    fprintf(stdout, "\t%zu MiB (%zu KiB) for external disk potential-field\n", alloc_ext_pot_disk.device >> 20, alloc_ext_pot_disk.device >> 10);
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifndef SERIALIZED_EXECUTION
    fprintf(stdout, "\t%zu MiB (%zu B) for LET related data\n"  ,  let_mem.device >> 20,  let_mem.device);
    fprintf(stdout, "\t%zu MiB (%zu KiB) for domain decomposition related data\n", box_mem.device >> 20, box_mem.device >> 10);
    fprintf(stdout, "\t%zu MiB (%zu KiB) for enclosing ball related data\n", ball_mem.device >> 20, ball_mem.device >> 10);
#ifdef  SWITCH_WITH_J_PARALLELIZATION
    fprintf(stdout, "\t%zu MiB (%zu B) for external N-body particles\n"  ,  ext_mem.device >> 20,  ext_mem.device);
#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//SERIALIZED_EXECUTION
#ifdef  COMPARE_WITH_DIRECT_SOLVER
    fprintf(stdout, "\t%zu MiB (%zu GiB) for accuracy test data\n", alloc_acc_dev.device >> 20, alloc_acc_dev.device >> 30);
#endif//COMPARE_WITH_DIRECT_SOLVER
#ifndef BLOCK_TIME_STEP
    fprintf(stdout, "\t%zu B (%zu KiB) for time step value\n", alloc_dt_dev.device >> 20, alloc_dt_dev.device >> 10);
#endif//BLOCK_TIME_STEP

    fprintf(stdout, "\nBreakdown of allocated memory on host:\n");
    fprintf(stdout, "\t%zu MiB (%zu GiB) for N-body particles\n"  , body_mem.host >> 20, body_mem.host >> 30);
    fprintf(stdout, "\t%zu MiB (%zu GiB) for tree structure\n"    , tree_mem.host >> 20, tree_mem.host >> 30);
    fprintf(stdout, "\t%zu MiB for miscellaneous data\n", misc_mem.host >> 20);
    fprintf(stdout, "\t%zu MiB (%zu B) for tree walk buffer\n", alloc_buf_dev.host >> 20, alloc_buf_dev.host);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
    fprintf(stdout, "\t%zu MiB (%zu KiB) for external spherical potential-field\n", alloc_ext_pot_sphe.host >> 20, alloc_ext_pot_sphe.host >> 10);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
    fprintf(stdout, "\t%zu MiB (%zu KiB) for external disk potential-field\n", alloc_ext_pot_disk.host >> 20, alloc_ext_pot_disk.host >> 10);
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifndef SERIALIZED_EXECUTION
    fprintf(stdout, "\t%zu MiB (%zu B) for LET related data\n", let_mem.host >> 20, let_mem.host);
    fprintf(stdout, "\t%zu MiB (%zu KiB) for domain decomposition related data\n", box_mem.host >> 20, box_mem.host >> 10);
    fprintf(stdout, "\t%zu MiB (%zu B) for enclosing ball related data\n", ball_mem.host >> 20, ball_mem.host);
#ifdef  SWITCH_WITH_J_PARALLELIZATION
    fprintf(stdout, "\t%zu MiB (%zu B) for external N-body particles\n"  ,  ext_mem.host >> 20,  ext_mem.host);
#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//SERIALIZED_EXECUTION
#ifdef  COMPARE_WITH_DIRECT_SOLVER
    fprintf(stdout, "\t%zu MiB (%zu GiB) for accuracy test data\n", alloc_acc_dev.host >> 20, alloc_acc_dev.host >> 30);
#endif//COMPARE_WITH_DIRECT_SOLVER
#ifndef BLOCK_TIME_STEP
    fprintf(stdout, "\t%zu B (%zu KiB) for time step value\n", alloc_dt_dev.host >> 20, alloc_dt_dev.host >> 10);
#endif//BLOCK_TIME_STEP

    fprintf(stdout, "\n");
    fflush(NULL);
  }
#endif//OUTPUT_MEMORY_USAGE

#ifdef  EXEC_BENCHMARK
  /** preparation of the benchmark */
  /** declaration of counters */
  static wall_clock_time execTime[BENCHMARK_STEPS];/**< zero-clear by static modifier */
#ifdef  COUNT_INTERACTIONS
  static tree_metrics    treeProp[BENCHMARK_STEPS];/**< zero-clear by static modifier */
#endif//COUNT_INTERACTIONS

  static uint bench_begin;
  bench_begin = steps;
#endif//EXEC_BENCHMARK

#ifdef  REPORT_GPU_CLOCK_FREQUENCY
  /** preparation of the GPU monitoring */
  /** declaration of counters */
  static gpu_clock deviceMonitors[CLOCK_RECORD_STEPS];/**< zero-clear by static modifier */
  static int monitor_step;
  monitor_step = 0;
#endif//REPORT_GPU_CLOCK_FREQUENCY


  /** set N-body particles on device */
  copyParticle_hst2dev(num, ibody0, ibody0_dev
#ifdef  EXEC_BENCHMARK
		       , &execTime[0]
#endif//EXEC_BENCHMARK
		       );


#ifndef SERIALIZED_EXECUTION
#ifndef DISABLE_AUTO_TUNING
  static loadImbalanceDetector balancer;
  balancer.enable  = false;
  balancer.execute = false;
  static double exchangeInterval = 0.0;/**< corresponding to rebuildTree.interval */
  static autoTuningParam exchangeParam;
  initStatVal(&exchangeParam.linearStats);  initGuessTime(&exchangeParam.linearGuess);
  initStatVal(&exchangeParam. powerStats);  initGuessTime(&exchangeParam. powerGuess);
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  initStatVal  (&exchangeParam.parabolicStats);
  initGuessTime(&exchangeParam.parabolicGuess);
#endif//USE_PARABOLIC_GROWTH_MODEL
#endif//DISABLE_AUTO_TUNING

#if 1
  if( steps == 0 )
    elapsed.sum_excg = 1.0f;
#endif
  updateDomain
    (&num, num_max, &Ni,
     &ibody0_dev, &ibody1_dev, &ibody0, &ibody1, Ntot,
     domBoundary, particlePos_hst, particlePos_dev, domCfg, domDecKey,
     sample, samplePos0, samplePos1, soaPH_dev, devProp, devInfo,
#ifndef DISABLE_AUTO_TUNING
     &exchangeParam, &exchangeInterval,
#endif//DISABLE_AUTO_TUNING
     ormCfg, repCfg,
#ifdef  CARE_EXTERNAL_PARTICLES
     &location,
#endif//CARE_EXTERNAL_PARTICLES
     letcfg, iparticleSendBuf, iparticleRecvBuf, &elapsed, &brentDistance, &brentHistory
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
#endif//SERIALIZED_EXECUTION


  /** start the numerical simulation */
  int numCell = NUM_ALLOC_TREE_CELL;
  int numNode = NUM_ALLOC_TREE_NODE;

  int bottomLev;
#ifndef SERIALIZED_EXECUTION
  int numOld = num;
#endif//SERIALIZED_EXECUTION
  buildTreeStructure
    (num, &ibody0_dev, &ibody1_dev, soaPH_dev, devProp,
#ifdef  CUB_AVAILABLE
     soaPH_pre,
#endif//CUB_AVAILABLE
#ifndef SERIALIZED_EXECUTION
     &numOld,
#endif//SERIALIZED_EXECUTION
     &bottomLev, &numCell, &numNode, bottomLev_dev, scanNum_dev, numCell_dev, numNode_dev, soaMakeBuf, soaCell_dev, soaNode_dev
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
     , &location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
     , &start
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
#ifdef  REPORT_COMPUTE_RATE
  rebuild_freq++;
#endif//REPORT_COMPUTE_RATE

#ifdef  BLOCK_TIME_STEP
  prediction_dev
    (num, time, ibody0_dev
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
#endif//BLOCK_TIME_STEP

  setMultipoleMoment
    (NUM_PHKEY_LEVEL - 1, soaCell_dev, numNode, soaNode_dev, num, ibody0_dev, soaMakeBuf, devProp
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
     , soaCell_hst, treeProp[steps - bench_begin].level
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );

  if( dropPrevTune == 1 ){
    /** first call of examineParticleSeparation() */
    /** therefore, brentDistance.a, brentDistance.b, brentDistance.x.pos are set by calling brentInit1st() */
    examineParticleSeparation
      (num, ibody0_dev, &brentDistance
#ifdef  MPI_MAX_FOR_RMAX_IN_AUTO_TUNING
       , letcfg.comm
#endif//MPI_MAX_FOR_RMAX_IN_AUTO_TUNING
#ifdef  EXEC_BENCHMARK
       , execTime
#endif//EXEC_BENCHMARK
       );
    brentDistance.u = brentDistance.x;
  }/* if( dropPrevTune == 1 ){ */

  /** commit i-particle groups */
  int Ngrp;
  /** first call of configDistribution() */
  /** do not call brentCalc1st() */
#ifdef  DISABLE_AUTO_TUNING
  brentDistance.u.pos = brent_frac * brentDistance.b;
#endif//DISABLE_AUTO_TUNING
  configDistribution
    (num, inumPerLane, &Ngrp, maxNgrp, laneInfo_hst, laneInfo_dev,
     ibody0_dev, &brentDistance, inum_dev, inum_hst
#ifdef  BLOCK_TIME_STEP
     , laneTime_dev
#endif//BLOCK_TIME_STEP
     , start, &elapsed
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
  static double reduce = 1.0;


#ifndef SERIALIZED_EXECUTION
  /** find center of enclosing ball */
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
  getEnclosingBox_dev(Ngrp, Ngrp, laneInfo_hst, ibody0_dev, soaPH_dev, devProp);
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  USE_ENCLOSING_BALL_FOR_LET
  getApproxEnclosingBall_dev
    (num, ibody0_dev
#ifdef  OCTREE_BASED_SEARCH
     , soaCell_dev, soaNode_dev, soaMakeBuf
#else///OCTREE_BASED_SEARCH
     , soaEB, soaPH_dev, devProp
#endif//OCTREE_BASED_SEARCH
     , stream_let[Nstream_let - 1]
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
#endif//USE_ENCLOSING_BALL_FOR_LET

  /** preparation to construct LET */
  calc_r2max_dev
    (Ngrp, laneInfo_dev, &ibody0_dev, soaGEO_dev
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
  shareNodePosition
    (letcfg.size, nodeInfo,
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     min, *(ibody0_dev.min_hst),
     max, *(ibody0_dev.max_hst),
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  USE_ENCLOSING_BALL_FOR_LET
     ipos, *(ibody0_dev.encBall_hst),
#endif//USE_ENCLOSING_BALL_FOR_LET
#ifdef  GADGET_MAC
     amin, ibody0_dev.amin,
#endif//GADGET_MAC
     letcfg);
  guessLETpartition(letcfg.size, nodeInfo, numNode,
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
		    *(ibody0_dev.icom_hst),
#else///USE_RECTANGULAR_BOX_FOR_LET
		    *(ibody0_dev.encBall_hst),
#endif//USE_RECTANGULAR_BOX_FOR_LET
		    letcfg);
#endif//SERIALIZED_EXECUTION


  if( steps == 0 ){
#ifdef  GADGET_MAC
    /** set Barnes-Hut MAC */
    enforceBarnesHutMAC_dev(Ni, ibody0_dev, numNode, soaNode_dev);
    /** force calculation based on opening criterion (theta = 1.0) */
    calcGravity_dev
      (
#ifdef  BLOCK_TIME_STEP
       Ngrp, &reduce,
#endif//BLOCK_TIME_STEP
       Ngrp, laneInfo_dev, ibody0_dev, soaNode_dev, soaWalk_dev, &sinfo, devProp, &elapsed.walkTree[0]
#   if  !defined(BLOCK_TIME_STEP) || defined(COMPARE_WITH_DIRECT_SOLVER) || defined(COUNT_INTERACTIONS) || defined(PRINT_PSEUDO_PARTICLE_INFO)
       , Ni
#endif//!defined(BLOCK_TIME_STEP) || defined(COMPARE_WITH_DIRECT_SOLVER) || defined(COUNT_INTERACTIONS) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
       , pot_tbl_sphe
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
       , pot_tbl_disk
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
       , file
#endif//PRINT_PSEUDO_PARTICLE_INFO
#   if  defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
       , cycles_hst, cycles_dev
#ifndef SERIALIZED_EXECUTION
       , cycles_dist_hst, cycles_dist_dev
#endif//SERIALIZED_EXECUTION
#endif//defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifndef SERIALIZED_EXECUTION
       , &elapsed, numNode
#ifdef  MPI_VIA_HOST
       , soaNode_hst
#endif//MPI_VIA_HOST
       , letcfg.size, nodeInfo, Nstream_let, stream_let, letcfg
#ifdef  MONITOR_LETGEN_TIME
       , cycles_let_hst, cycles_let_dev
#endif//MONITOR_LETGEN_TIME
#ifdef  SWITCH_WITH_J_PARALLELIZATION
       , false, Ni_local, Ni_total
       , Ni_list, head_list, grpNum_list, grpNum_disp
       , maxNgrp_ext, laneInfo_ext_dev, laneInfo_ext_hst
       , laneInfo_hst, ibody_dist0_dev
#ifdef  MPI_VIA_HOST
       , ibody_dist0, ibody_dist1
#endif//MPI_VIA_HOST
#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
       , treeinfo_dev
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
       , deviceMonitors, &monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  COMPARE_WITH_DIRECT_SOLVER
       , true
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
       , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#endif//COMPARE_WITH_DIRECT_SOLVER
       );
    /** set GADGET-MAC */
    recoverGADGET_MAC_dev(numNode, soaNode_dev);
    /** reset stop watch */
    elapsed.walkTree[0] = 0.0;
#ifndef SERIALIZED_EXECUTION
    /** find center of enclosing ball */
#ifdef  USE_ENCLOSING_BALL_FOR_LET
    getApproxEnclosingBall_dev
      (num, ibody0_dev
#ifdef  OCTREE_BASED_SEARCH
       , soaCell_dev, soaNode_dev, soaMakeBuf
#else///OCTREE_BASED_SEARCH
       , soaEB, soaPH_dev, devProp
#endif//OCTREE_BASED_SEARCH
       , stream_let[Nstream_let - 1]
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
#endif//USE_ENCLOSING_BALL_FOR_LET
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
    getEnclosingBox_dev(Ngrp, Ngrp, laneInfo_hst, ibody0_dev, soaPH_dev, devProp);
#endif//USE_RECTANGULAR_BOX_FOR_LET

    calc_r2max_dev
      (Ngrp, laneInfo_dev, &ibody0_dev, soaGEO_dev
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
    shareNodePosition
      (letcfg.size, nodeInfo,
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
       min, *(ibody0_dev.min_hst),
       max, *(ibody0_dev.max_hst),
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  USE_ENCLOSING_BALL_FOR_LET
       ipos, *(ibody0_dev.encBall_hst),
#endif//USE_ENCLOSING_BALL_FOR_LET
#ifdef  GADGET_MAC
       amin, ibody0_dev.amin,
#endif//GADGET_MAC
       letcfg);
#endif//SERIALIZED_EXECUTION
#endif//GADGET_MAC


    /** calculate gravitational acceleration and potential */
    calcGravity_dev
      (
#ifdef  BLOCK_TIME_STEP
       Ngrp, &reduce,
#endif//BLOCK_TIME_STEP
       Ngrp, laneInfo_dev, ibody0_dev, soaNode_dev, soaWalk_dev, &sinfo, devProp, &elapsed.walkTree[0]
#   if  !defined(BLOCK_TIME_STEP) || defined(COMPARE_WITH_DIRECT_SOLVER) || defined(COUNT_INTERACTIONS) || defined(PRINT_PSEUDO_PARTICLE_INFO)
       , Ni
#endif//!defined(BLOCK_TIME_STEP) || defined(COMPARE_WITH_DIRECT_SOLVER) || defined(COUNT_INTERACTIONS) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
       , pot_tbl_sphe
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
       , pot_tbl_disk
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
       , file
#endif//PRINT_PSEUDO_PARTICLE_INFO
#   if  defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
       , cycles_hst, cycles_dev
#ifndef SERIALIZED_EXECUTION
       , cycles_dist_hst, cycles_dist_dev
#endif//SERIALIZED_EXECUTION
#endif//defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifndef SERIALIZED_EXECUTION
       , &elapsed, numNode
#ifdef  MPI_VIA_HOST
       , soaNode_hst
#endif//MPI_VIA_HOST
       , letcfg.size, nodeInfo, Nstream_let, stream_let, letcfg
#ifdef  MONITOR_LETGEN_TIME
       , cycles_let_hst, cycles_let_dev
#endif//MONITOR_LETGEN_TIME
#ifdef  SWITCH_WITH_J_PARALLELIZATION
       , false, Ni_local, Ni_total
       , Ni_list, head_list, grpNum_list, grpNum_disp
       , maxNgrp_ext, laneInfo_ext_dev, laneInfo_ext_hst
       , laneInfo_hst, ibody_dist0_dev
#ifdef  MPI_VIA_HOST
       , ibody_dist0, ibody_dist1
#endif//MPI_VIA_HOST
#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
       , treeinfo_dev
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
       , deviceMonitors, &monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  COMPARE_WITH_DIRECT_SOLVER
       , true
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
       , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#endif//COMPARE_WITH_DIRECT_SOLVER
       );

#ifndef DISABLE_AUTO_TUNING
#ifndef SHARED_AUTO_TUNER
#ifdef  USE_CLOCK_CYCLES_FOR_BRENT_METHOD
    brentDistance.u.val += (double)(*cycles_hst);
#else///USE_CLOCK_CYCLES_FOR_BRENT_METHOD
    brentDistance.u.val += elapsed.walkTree[0];
#endif//USE_CLOCK_CYCLES_FOR_BRENT_METHOD
    brentHistory.totNum += Ngrp;
#else///SHARED_AUTO_TUNER
#ifdef  USE_CLOCK_CYCLES_FOR_BRENT_METHOD
    shared_brent_u = (double)(*cycles_hst);
#else///USE_CLOCK_CYCLES_FOR_BRENT_METHOD
    shared_brent_u = elapsed.walkTree[0];
#endif//USE_CLOCK_CYCLES_FOR_BRENT_METHOD
    shared_brent_num = Ngrp;
    chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &shared_brent_u, 1, MPI_DOUBLE, MPI_SUM, letcfg.comm));
    chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &shared_brent_num, 1, MPI_INT, MPI_SUM, letcfg.comm));
    brentDistance.u.val += shared_brent_u;
    brentHistory.totNum += shared_brent_num;
#endif//SHARED_AUTO_TUNER
#endif//DISABLE_AUTO_TUNING

    rebuild.interval += 1.0;
#ifndef DISABLE_AUTO_TUNING
    linearModel(&(rebuildParam.linearGuess), &(rebuildParam.linearStats), rebuild.interval, elapsed.walkTree[0], 1.0);
    powerModel (&(rebuildParam. powerGuess), &(rebuildParam. powerStats), rebuild.interval, elapsed.walkTree[0], 1.0);
#ifdef  USE_PARABOLIC_GROWTH_MODEL
    parabolicModel(&(rebuildParam.parabolicGuess), &(rebuildParam.parabolicStats), rebuild.interval, elapsed.walkTree[0], 1.0);
#ifdef  USE_ADDITIONAL_SWITCH
    useParabolicGuess |= (elapsed.walkTree[0] < elapsed.makeTree);
#endif//USE_ADDITIONAL_SWITCH
#endif//USE_PARABOLIC_GROWTH_MODEL

#   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
    rebuild.avg += reduce;
    rebuild.var += reduce * reduce;
#endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
#endif//DISABLE_AUTO_TUNING

#ifdef  COUNT_INTERACTIONS
    copyCounters_dev2hst(Ni, treeinfo_dev, treeinfo);
    analyzeTreeMetrics(&treeProp[steps - bench_begin]);
    analyzeWalkStatistics(Ni, Ngrp, Nj  , &(treeProp[steps - bench_begin].Nj  ));
    analyzeWalkStatistics(Ni, Ngrp, Nbuf, &(treeProp[steps - bench_begin].Nbuf));
    for(int ii = 0; ii < Ni; ii++)
      treeProp[steps - bench_begin].Ninteractions += (ulong)Nj[ii];
    treeProp[steps - bench_begin].bufSize = soaWalk_dev.bufSize;
#endif//COUNT_INTERACTIONS

#ifndef EXEC_BENCHMARK
#ifdef  COMPARE_WITH_DIRECT_SOLVER
    ibody_direct_dev     = ibody0_dev;
    ibody_direct_dev.acc = direct_dev;
#ifdef  GADGET_MAC
    ibody_direct_dev.acc_old = direct_old_dev;
#endif//GADGET_MAC
#endif//COMPARE_WITH_DIRECT_SOLVER

#ifdef  REPORT_COMPUTE_RATE
    brent_rate = brentDistance.u.pos / brentDistance.b;
#endif//REPORT_COMPUTE_RATE
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
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
       , deviceMonitors, &monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  REPORT_COMPUTE_RATE
       , &clock_prev, &time_prev, &step_prev, ft, &brent_rate, &brent_prev, &rebuild_freq
#endif//REPORT_COMPUTE_RATE
#ifdef  COMPARE_WITH_DIRECT_SOLVER
       , Ni, ibody_direct_dev, devProp, direct, Ngrp, laneInfo_dev
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
       , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
       , pot_tbl_sphe
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
       , pot_tbl_disk
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
       , accfile
#endif//COMPARE_WITH_DIRECT_SOLVER
       );
#endif//EXEC_BENCHMARK


    /** initialize time step */
#ifndef BLOCK_TIME_STEP
    setTimeStep_dev
      (num, ibody0, eta, eps, dt_dev, &dt
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
  }/* if( steps == 0 ){ */


  /** initialize time step */
#ifdef  BLOCK_TIME_STEP
  adjustParticleTime_dev
    (Ngrp, laneInfo_dev, laneTime_dev, eps, eta, ibody0_dev
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
#endif//BLOCK_TIME_STEP

#ifdef  LEAP_FROG_INTEGRATOR
  /** time integration for velocity to implement leap-frog method */
  advVel_dev
    (Ni, ibody0_dev, HALF * (real)dt
#ifdef  EXEC_BENCHMARK
     , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
     );
#endif//LEAP_FROG_INTEGRATOR


  /** calculate time evolution */
#ifndef SERIALIZED_EXECUTION
  bool existNewTree = false;
#endif//SERIALIZED_EXECUTION
#ifndef DISABLE_AUTO_TUNING
#ifdef  SERIALIZED_EXECUTION
  static const double brent_method_allow = BRENT_METHOD_ALLOW;
#else///SERIALIZED_EXECUTION
#ifdef  DISABLE_EXCG_BODIES_BEFORE_SATURATION
  bool allowBodyExcg = false;
  static double brent_method_allow = BRENT_METHOD_ALLOW_LAUNCH;
#else///DISABLE_EXCG_BODIES_BEFORE_SATURATION
  static const double brent_method_allow = BRENT_METHOD_ALLOW;
#endif//DISABLE_EXCG_BODIES_BEFORE_SATURATION
#endif//SERIALIZED_EXECUTION
#   if  NMAX_FOR_PERTURBATION_ABOUT_BRENT > 0
  static int brent_perturb_count = 0;
#endif//NMAX_FOR_PERTURBATION_ABOUT_BRENT > 0
#endif//DISABLE_AUTO_TUNING
  while( time < ft ){
#ifdef  MONITOR_SIMULATION_STATUS
    if( mpi.rank == 0 ){
      static char date[64];
      getPresentDateInStrings(date);
      __FPRINTF__(stdout, "t = %e(%zu step(s)) on %s", time, steps, date);
    }/* if( mpi.rank == 0 ){ */
#endif//MONITOR_SIMULATION_STATUS
    __NOTE__("t = %e(%zu step(s)), ft = %e\n", time, steps, ft);


#ifdef  EXEC_BENCHMARK
    static double currentTime;
#ifdef  SERIALIZED_EXECUTION
    currentTime = getPresentTimeInMin();
#else///SERIALIZED_EXECUTION
    if( mpi.rank == 0 )
      currentTime = getPresentTimeInMin();
    chkMPIerr(MPI_Bcast(&currentTime, 1, MPI_DOUBLE, 0, mpi.comm));
#endif//SERIALIZED_EXECUTION

    if( (steps >= (bench_begin + BENCHMARK_STEPS - 1)) || (currentTime >= (formerTime + TFIN_MIN_BENCHMARK)) ){
      /** output results of benchmark (wall clock time) */
      dumpBenchmark (jobID, file, 1 + steps - bench_begin, execTime);
#ifdef  COUNT_INTERACTIONS
      /** output results of benchmark (statistics of tree) */
      dumpStatistics(jobID, file, 1 + steps - bench_begin, MAXIMUM_PHKEY_LEVEL, treeProp);
#endif//COUNT_INTERACTIONS

      /** exit loop about time evolution */
      break;
    }/* if( (steps >= (bench_begin + BENCHMARK_STEPS - 1)) || (currentTime >= (formerTime + TFIN_MIN_BENCHMARK)) ){ */
#endif//EXEC_BENCHMARK


#   if  defined(BLOCK_TIME_STEP) || !defined(EXEC_BENCHMARK)
    static uint present;
#endif//defined(BLOCK_TIME_STEP) || !defined(EXEC_BENCHMARK)
    steps++;

    /** update octree structure if necessary */
    /** rebuild tree structure if required */
#ifndef SERIALIZED_EXECUTION
#ifndef DISABLE_AUTO_TUNING
    /** detect slow down */
    if( !rebuild.reuse ){
      exchangeInterval += 1.0;
      linearModel(&(exchangeParam.linearGuess), &(exchangeParam.linearStats), exchangeInterval, elapsed.sum_rebuild, 1.0);
      powerModel (&(exchangeParam. powerGuess), &(exchangeParam. powerStats), exchangeInterval, elapsed.sum_rebuild, 1.0);
#ifdef  USE_PARABOLIC_GROWTH_MODEL
      parabolicModel(&(exchangeParam.parabolicGuess), &(exchangeParam.parabolicStats), exchangeInterval, elapsed.sum_rebuild, 1.0);
#endif//USE_PARABOLIC_GROWTH_MODEL
      balancer.execute = false;

      /** detect slow down of gravity calculation */
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
	if( guess > elapsed.excg ){
	  balancer.execute = true;
#ifdef  OUTPUT_EXCHANGE_SIGNAL
	  __FPRINTF__(stdout, "exchange by detecting slow down @ %zu-th step\n", steps);
#endif//OUTPUT_EXCHANGE_SIGNAL
	}/* if( guess > elapsed.excg ){ */
      }/* if( balancer.enable ){ */
    }/* if( !rebuild.reuse ){ */

    chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &(balancer.execute), 1, MPI_BOOL, MPI_LOR, letcfg.comm));
    /** detect load imbalance */
    if( (balancer.enable) && (!balancer.execute) ){
      balancer.tmin = balancer.tmax = elapsed.sum_excg;
      chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &(balancer.tmax), 1, MPI_DOUBLE, MPI_MAX, letcfg.comm));
      if( balancer.tmax > minimumExecutionTime ){
	chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &(balancer.tmin), 1, MPI_DOUBLE, MPI_MIN, letcfg.comm));
	if( balancer.tmin < (loadImbalanceCrit * balancer.tmax) ){
	  balancer.execute = true;
#ifdef  OUTPUT_EXCHANGE_SIGNAL
	  __FPRINTF__(stdout, "exchange by detecting load imbalance (tmin = %e, tmax = %e) @ %zu-th step\n", balancer.tmin, balancer.tmax, steps);
#endif//OUTPUT_EXCHANGE_SIGNAL
	}/* if( balancer.tmin < (loadImbalanceCrit * balancer.tmax) ){ */
      }/* if( balancer.tmax > minimumExecutionTime ){ */
    }/* if( (balancer.enable) && (!balancer.execute) ){ */
    __NOTE__("balancer.execute = %d @ rank %d\n", balancer.execute, letcfg.rank);

    if( balancer.execute ){
#ifdef  MONITOR_SIMULATION_STATUS
      if( mpi.rank == 0 ){
	static char date[64];
	getPresentDateInStrings(date);
	__FPRINTF__(stdout, "domain re-decomposition on %s", date);
      }/* if( mpi.rank == 0 ){ */
#endif//MONITOR_SIMULATION_STATUS
      updateDomain
	(&num, num_max, &Ni,
	 &ibody0_dev, &ibody1_dev, &ibody0, &ibody1, Ntot,
	 domBoundary, particlePos_hst, particlePos_dev, domCfg, domDecKey,
	 sample, samplePos0, samplePos1, soaPH_dev, devProp, devInfo,
#ifndef DISABLE_AUTO_TUNING
	 &exchangeParam, &exchangeInterval,
#endif//DISABLE_AUTO_TUNING
	 ormCfg, repCfg,
#ifdef  CARE_EXTERNAL_PARTICLES
	 &location,
#endif//CARE_EXTERNAL_PARTICLES
	 letcfg, iparticleSendBuf, iparticleRecvBuf, &elapsed, &brentDistance, &brentHistory
#ifdef  EXEC_BENCHMARK
	 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	 );

      rebuild.reuse = 0;
      balancer.enable = false;

#ifdef  DISABLE_EXCG_BODIES_BEFORE_SATURATION
      allowBodyExcg = false;
      brent_method_allow = BRENT_METHOD_ALLOW_LAUNCH;
#endif//DISABLE_EXCG_BODIES_BEFORE_SATURATION
    }/* if( balancer.execute ){ */
#else///DISABLE_AUTO_TUNING
    if( !rebuild.reuse ){
#ifdef  MONITOR_SIMULATION_STATUS
      if( mpi.rank == 0 ){
	static char date[64];
	getPresentDateInStrings(date);
	__FPRINTF__(stdout, "domain re-decomposition on %s", date);
      }/* if( mpi.rank == 0 ){ */
#endif//MONITOR_SIMULATION_STATUS
      updateDomain
	(&num, num_max, &Ni,
	 &ibody0_dev, &ibody1_dev, &ibody0, &ibody1, Ntot,
	 domBoundary, particlePos_hst, particlePos_dev, domCfg, domDecKey,
	 sample, samplePos0, samplePos1, soaPH_dev, devProp, devInfo,
#ifndef DISABLE_AUTO_TUNING
	 &exchangeParam, &exchangeInterval,
#endif//DISABLE_AUTO_TUNING
	 ormCfg, repCfg,
#ifdef  CARE_EXTERNAL_PARTICLES
	 &location,
#endif//CARE_EXTERNAL_PARTICLES
	 letcfg, iparticleSendBuf, iparticleRecvBuf, &elapsed, &brentDistance, &brentHistory
#ifdef  EXEC_BENCHMARK
	 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	 );
    }/* if( !rebuild.reuse ){ */
#endif//DISABLE_AUTO_TUNING
#endif//SERIALIZED_EXECUTION


    /** rebuild the tree structure */
    if( !rebuild.reuse ){
#ifndef DISABLE_AUTO_TUNING
#ifndef SERIALIZED_EXECUTION
      elapsed.sum_rebuild = 0.0;
#endif//SERIALIZED_EXECUTION

#   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
      /* if Var(reduce) < (Avg(reduce))^2, set time of all particles same */
      /* if( ((rebuildInterval * reduceVar) / (reduceAvg * reduceAvg) - 1.0) < (1.0 - DBL_EPSILON) ) */
      if( DBL_EPSILON + rebuild.interval * rebuild.var < (2.0 * rebuild.avg * rebuild.avg) )
	rebuild.adjust = true;

      rebuild.avg = 0.0;
      rebuild.var = 0.0;
#endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
#endif//DISABLE_AUTO_TUNING
      rebuild.interval = 0.0;

#ifndef DISABLE_AUTO_TUNING
      initStatVal(&(rebuildParam.   linearStats));      initGuessTime(&(rebuildParam.   linearGuess));
      initStatVal(&(rebuildParam.    powerStats));      initGuessTime(&(rebuildParam.    powerGuess));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
      initStatVal(&(rebuildParam.parabolicStats));      initGuessTime(&(rebuildParam.parabolicGuess));
#ifdef  USE_ADDITIONAL_SWITCH
      useParabolicGuess = 0;
#endif//USE_ADDITIONAL_SWITCH
#endif//USE_PARABOLIC_GROWTH_MODEL
#endif//DISABLE_AUTO_TUNING

#ifdef  EXEC_BENCHMARK
      printf("#rebuild @ %zu-th step\n", steps);
      fflush(stdout);
#endif//EXEC_BENCHMARK
#ifndef SERIALIZED_EXECUTION
      __NOTE__("rebuild @ %zu-th step, t = %e, sum_rebuild = %e\n", steps, time, elapsed.sum_rebuild);
#endif//SERIALIZED_EXECUTION

      buildTreeStructure
	(num, &ibody0_dev, &ibody1_dev, soaPH_dev, devProp,
#ifdef  CUB_AVAILABLE
	 soaPH_pre,
#endif//CUB_AVAILABLE
#ifndef SERIALIZED_EXECUTION
	 &numOld,
#endif//SERIALIZED_EXECUTION
	 &bottomLev, &numCell, &numNode, bottomLev_dev, scanNum_dev, numCell_dev, numNode_dev, soaMakeBuf, soaCell_dev, soaNode_dev
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
	 , &location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
	 , &start
#ifdef  EXEC_BENCHMARK
	 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	 );
#ifdef  REPORT_COMPUTE_RATE
      rebuild_freq++;
#endif//REPORT_COMPUTE_RATE

#ifndef DISABLE_AUTO_TUNING
      brentDistance.u.val /= (double)brentHistory.totNum;
      brentHistory.totNum = 0;
      brentHistory.interval++;
#endif//DISABLE_AUTO_TUNING
      static bool initialized = false;
      if( initialized ){
	/** from the second time, managing brentStatus is necessary */
	/** check # of sequences executed and the evolution of elapsed time */
#ifndef DISABLE_AUTO_TUNING
	bool perturbBrent = false;
	if( brentHistory.interval > BRENT_METHOD_LAUNCH ){
	  /** when the elapsed time grows continuously, parturb the brentStatus */
	  if( brentDistance.u.val > brentHistory.previous ){
	    brentHistory.degraded++;
	    perturbBrent = (brentHistory.degraded > BRENT_METHOD_MODIFY) || (brentDistance.u.val > (brent_method_allow * brentHistory.previous));
	  }/* if( brentDistance.u.val > brentHistory.previous ){ */
	  else
	    brentHistory.degraded = 0;
	}/* if( brentHistory.interval > BRENT_METHOD_LAUNCH ){ */
	brentHistory.previous = brentDistance.u.val;

	if( perturbBrent ){
#ifdef  DISABLE_EXCG_BODIES_BEFORE_SATURATION
	  if( brentHistory.degraded > BRENT_METHOD_MODIFY ){
	    allowBodyExcg = true;
	    brent_method_allow = BRENT_METHOD_ALLOW;
	  }/* if( brentHistory.degraded > BRENT_METHOD_MODIFY ){ */
#endif//DISABLE_EXCG_BODIES_BEFORE_SATURATION
#if 0
	  __FPRINTF__(stdout, "perturb: interval = %d, degraded = %d @ %zu-step.\n", brentHistory.interval, brentHistory.degraded, steps);
#endif
#endif//DISABLE_AUTO_TUNING

#ifdef  REPORT_COMPUTE_RATE
	  brent_rate += (brentDistance.u.pos / brentDistance.b) * (double)(steps - brent_prev);
	  /* if( fpclassify(brent_rate) == FP_NAN ){ */
	  /*   __FPRINTF__(stderr, "brent_rate = %e: pos = %e, b = %e, step = %zu, prev = %zu\n", brent_rate, brentDistance.u.pos, brentDistance.b, steps, brent_prev); */
	  /* }/\* if( fpclassify(brent_rate) == FP_NAN ){ *\/ */
	  brent_prev = steps;
#endif//REPORT_COMPUTE_RATE

	  examineParticleSeparation
	    (num, ibody0_dev, &brentDistance
#ifdef  MPI_MAX_FOR_RMAX_IN_AUTO_TUNING
	     , letcfg.comm
#endif//MPI_MAX_FOR_RMAX_IN_AUTO_TUNING
#ifdef  EXEC_BENCHMARK
	     , execTime
#endif//EXEC_BENCHMARK
	     );

#ifndef DISABLE_AUTO_TUNING
#   if  NMAX_FOR_PERTURBATION_ABOUT_BRENT > 0
	  brent_perturb_count++;
	  if( brent_perturb_count > NMAX_FOR_PERTURBATION_ABOUT_BRENT ){
	    brentInit1st(&brentDistance, brentDistance.a, brentDistance.b);
	    brentInit2nd(&brentDistance);
	    brentDistance.u = brentDistance.x;
	    brent_perturb_count = 0;
	  }/* if( brent_perturb_count > NMAX_FOR_PERTURBATION_ABOUT_BRENT ){ */
#endif//NMAX_FOR_PERTURBATION_ABOUT_BRENT > 0

	  brentHistory.interval = 0;
	  brentHistory.degraded = 0;
	}/* if( perturbBrent ){ */
	else{
	  /** when the brentStatus is conserved */
	  brentCalc2nd(&brentDistance);
	}/* else{ */
#endif//DISABLE_AUTO_TUNING
      }/* if( initialized ){ */
      else{
	/** the first execution of tree rebuild */
#ifndef DISABLE_AUTO_TUNING
	brentDistance.x = brentDistance.u;
	brentInit2nd(&brentDistance);
#endif//DISABLE_AUTO_TUNING
	initialized = true;
      }/* else{ */

      /** from the second time, brentCalc1st() is called internally */
#ifdef  DISABLE_AUTO_TUNING
      brentDistance.u.pos = brent_frac * brentDistance.b;
#endif//DISABLE_AUTO_TUNING
      configDistribution
	(num, inumPerLane, &Ngrp, maxNgrp, laneInfo_hst, laneInfo_dev,
	 ibody0_dev, &brentDistance, inum_dev, inum_hst
#ifdef  BLOCK_TIME_STEP
	 , laneTime_dev
#endif//BLOCK_TIME_STEP
	 , start, &elapsed
#ifdef  EXEC_BENCHMARK
	 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	 );

#ifndef SERIALIZED_EXECUTION
      existNewTree = true;
#ifndef DISABLE_AUTO_TUNING
      if( !balancer.execute )      	balancer.enable  =  true;
      else	                        balancer.execute = false;
#endif//DISABLE_AUTO_TUNING
#endif//SERIALIZED_EXECUTION
    }/* if( !rebuild.reuse ){ */


    /** share information on domain decomposition in the next time step */
#ifndef SERIALIZED_EXECUTION
    chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &existNewTree, 1, MPI_BOOL, MPI_LOR, letcfg.comm));
#ifndef DISABLE_AUTO_TUNING
#ifdef  DISABLE_EXCG_BODIES_BEFORE_SATURATION
    balancer.enable &= allowBodyExcg;
    chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &(balancer.enable), 1, MPI_BOOL, MPI_LAND, letcfg.comm));
#else///DISABLE_EXCG_BODIES_BEFORE_SATURATION
    chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &(balancer.enable), 1, MPI_BOOL, MPI_LOR, letcfg.comm));
#endif//DISABLE_EXCG_BODIES_BEFORE_SATURATION
#endif//DISABLE_AUTO_TUNING
#endif//SERIALIZED_EXECUTION


    /** time integration */
#ifdef  BLOCK_TIME_STEP
    setLaneTime_dev
      (Ngrp, laneInfo_dev, laneTime_dev, ibody0_dev
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
    int grpNum = 0;
    setTimeStep_dev
      (Ngrp, laneInfo_dev, laneTime_dev, &grpNum, ibody0_dev,
       time, &time, &dt, rebuild.adjust, invSnapshotInterval, previous, &present
#ifndef SERIALIZED_EXECUTION
       , letcfg
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
    rebuild.adjust = false;
#else///BLOCK_TIME_STEP
#ifndef LEAP_FROG_INTEGRATOR
    /* initialize time step */
    setTimeStep_dev
      (Ni, ibody0, eta, eps, dt_dev, &dt
#ifndef SERIALIZED_EXECUTION
       , letcfg
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
#endif//LEAP_FROG_INTEGRATOR
    time += dt;
#endif//BLOCK_TIME_STEP

#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && defined(TIME_BASED_MODIFICATION)
    location.elapsed = (float)((double)location.elapsed + dt);
    if( location.dtmin > dt ){
      location.dtmin = dt;
      location.dtinv = 1.0f / location.dtmin;
    }/* if( location.dtmin < dt ){ */
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && defined(TIME_BASED_MODIFICATION)


    /** orbit integration (predict step) */
#ifdef  BLOCK_TIME_STEP
    prediction_dev
      (num, time, ibody0_dev
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
#else///BLOCK_TIME_STEP
#ifndef LEAP_FROG_INTEGRATOR
    /** predict velocity and update position of i-particles */
    advVel_dev
      (Ni, ibody0_dev, HALF * (real)dt
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
    advPos_dev
      (Ni, ibody0_dev,        (real)dt
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
#else///LEAP_FROG_INTEGRATOR
    /** update position of i-particles */
    advPos_dev
      (Ni, ibody0_dev, (real)dt
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
#endif//LEAP_FROG_INTEGRATOR
#endif//BLOCK_TIME_STEP


    /** calculate multipole moment of tree nodes */
    setMultipoleMoment
      (bottomLev, soaCell_dev, numNode, soaNode_dev, num, ibody0_dev, soaMakeBuf, devProp
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
       , soaCell_hst, treeProp[steps - bench_begin].level
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );


#ifndef SERIALIZED_EXECUTION

#ifdef  SWITCH_WITH_J_PARALLELIZATION
    selectCommunicationMode(grpNum, Ngrp, laneInfo_hst, letcfg, &transferMode, &Ni_local, &Ni_total, Ni_list, head_list, grpNum_list);

    if( !transferMode )
#endif//SWITCH_WITH_J_PARALLELIZATION
      {
	/** preparation to construct LET */
	/** find center of enclosing ball */
#ifdef  USE_ENCLOSING_BALL_FOR_LET
	getApproxEnclosingBall_dev
	  (num, ibody0_dev
#ifdef  OCTREE_BASED_SEARCH
	   , soaCell_dev, soaNode_dev, soaMakeBuf
#else///OCTREE_BASED_SEARCH
	   , soaEB, soaPH_dev, devProp
#endif//OCTREE_BASED_SEARCH
	   , stream_let[Nstream_let - 1]
#ifdef  EXEC_BENCHMARK
	   , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	   );
#endif//USE_ENCLOSING_BALL_FOR_LET
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
	getEnclosingBox_dev(existNewTree ? Ngrp : grpNum, Ngrp, laneInfo_hst, ibody0_dev, soaPH_dev, devProp);
#endif//USE_RECTANGULAR_BOX_FOR_LET

	/*     calc_r2max_dev */
	/*       (Ngrp, laneInfo_dev, &ibody0_dev, soaGEO_dev */
	/* #ifdef  EXEC_BENCHMARK */
	/*        , &execTime[steps - bench_begin] */
	/* #endif//EXEC_BENCHMARK */
	/*        ); */
	calc_r2max_dev
	  (existNewTree ? Ngrp : grpNum, laneInfo_dev, &ibody0_dev, soaGEO_dev
#ifdef  EXEC_BENCHMARK
	   , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
	   );
	shareNodePosition
	  (letcfg.size, nodeInfo,
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
	   min, *(ibody0_dev.min_hst),
	   max, *(ibody0_dev.max_hst),
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  USE_ENCLOSING_BALL_FOR_LET
	   ipos, *(ibody0_dev.encBall_hst),
#endif//USE_ENCLOSING_BALL_FOR_LET
#ifdef  GADGET_MAC
	   amin, ibody0_dev.amin,
#endif//GADGET_MAC
	   letcfg);
	/** this function must be called when tree structure is rebuild */
	if( existNewTree )
	  guessLETpartition(letcfg.size, nodeInfo, numNode,
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
			    *(ibody0_dev.icom_hst),
#else///USE_RECTANGULAR_BOX_FOR_LET
			    *(ibody0_dev.encBall_hst),
#endif//USE_RECTANGULAR_BOX_FOR_LET
			    letcfg);
	existNewTree = false;
      }
#endif//SERIALIZED_EXECUTION


    /** calculate acceleration of i-particles by j-cells */
#ifdef  SHOW_NI_DEPENDENCE
    printf("#file: %s\n", file);
    printf("#NTHREADS = %d, TSUB = %d, NWARP = %d\n", NTHREADS, TSUB, NWARP);
    printf("#grpNum\tNi_active\ttime(sec)\n");
    for(grpNum = Ngrp; grpNum >= NGROUPS; grpNum = (int)((double)grpNum * M_SQRT1_2)){
#endif//SHOW_NI_DEPENDENCE
      calcGravity_dev
	(
#ifdef  BLOCK_TIME_STEP
	 grpNum, &reduce,
#endif//BLOCK_TIME_STEP
	 Ngrp, laneInfo_dev, ibody0_dev, soaNode_dev, soaWalk_dev, &sinfo, devProp, &elapsed.walkTree[rebuild.reuse]
#   if  !defined(BLOCK_TIME_STEP) || defined(COMPARE_WITH_DIRECT_SOLVER) || defined(COUNT_INTERACTIONS) || defined(PRINT_PSEUDO_PARTICLE_INFO)
	 , Ni
#endif//!defined(BLOCK_TIME_STEP) || defined(COMPARE_WITH_DIRECT_SOLVER) || defined(COUNT_INTERACTIONS) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
	 , pot_tbl_sphe
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
	 , pot_tbl_disk
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
	 , file
#endif//PRINT_PSEUDO_PARTICLE_INFO
#   if  defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
	 , cycles_hst, cycles_dev
#ifndef SERIALIZED_EXECUTION
	 , cycles_dist_hst, cycles_dist_dev
#endif//SERIALIZED_EXECUTION
#endif//defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifndef SERIALIZED_EXECUTION
	 , &elapsed, numNode
#ifdef  MPI_VIA_HOST
	 , soaNode_hst
#endif//MPI_VIA_HOST
	 , letcfg.size, nodeInfo, Nstream_let, stream_let, letcfg
#ifdef  MONITOR_LETGEN_TIME
	 , cycles_let_hst, cycles_let_dev
#endif//MONITOR_LETGEN_TIME
#ifdef  SWITCH_WITH_J_PARALLELIZATION
	 , transferMode, Ni_local, Ni_total
	 , Ni_list, head_list, grpNum_list, grpNum_disp
	 , maxNgrp_ext, laneInfo_ext_dev, laneInfo_ext_hst
	 , laneInfo_hst, ibody_dist0_dev
#ifdef  MPI_VIA_HOST
	 , ibody_dist0, ibody_dist1
#endif//MPI_VIA_HOST
#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
	 , treeinfo_dev
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
	 , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
	 , deviceMonitors, &monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  COMPARE_WITH_DIRECT_SOLVER
	 , true
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	 , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#endif//COMPARE_WITH_DIRECT_SOLVER
	 );

#ifndef DISABLE_AUTO_TUNING
#ifndef SHARED_AUTO_TUNER
/* #ifdef  SWITCH_WITH_J_PARALLELIZATION */
/*       if( !transferMode ) */
/* #endif//SWITCH_WITH_J_PARALLELIZATION */
	{
#ifdef  USE_CLOCK_CYCLES_FOR_BRENT_METHOD
	  brentDistance.u.val += (double)(*cycles_hst);
#else///USE_CLOCK_CYCLES_FOR_BRENT_METHOD
	  brentDistance.u.val += elapsed.walkTree[rebuild.reuse];
#endif//USE_CLOCK_CYCLES_FOR_BRENT_METHOD
#ifdef  BLOCK_TIME_STEP
	  brentHistory.totNum += grpNum;
#else///BLOCK_TIME_STEP
	  brentHistory.totNum += Ngrp;
#endif//BLOCK_TIME_STEP
	}
#else///SHARED_AUTO_TUNER
#ifdef  USE_CLOCK_CYCLES_FOR_BRENT_METHOD
      shared_brent_u = (double)(*cycles_hst);
#else///USE_CLOCK_CYCLES_FOR_BRENT_METHOD
      shared_brent_u = elapsed.walkTree[rebuild.reuse];
#endif//USE_CLOCK_CYCLES_FOR_BRENT_METHOD
      chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &shared_brent_u, 1, MPI_DOUBLE, MPI_SUM, letcfg.comm));
#ifdef  SWITCH_WITH_J_PARALLELIZATION
      shared_brent_num = grpNum_disp[letcfg.size - 1] + grpNum_list[letcfg.size - 1];
#else///SWITCH_WITH_J_PARALLELIZATION
#ifdef  BLOCK_TIME_STEP
      shared_brent_num = grpNum;
#else///BLOCK_TIME_STEP
      shared_brent_num = Ngrp;
#endif//BLOCK_TIME_STEP
      chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &shared_brent_num, 1, MPI_INT, MPI_SUM, letcfg.comm));
#endif//SWITCH_WITH_J_PARALLELIZATION
      brentDistance.u.val += shared_brent_u;
      brentHistory.totNum += shared_brent_num;
#endif//SHARED_AUTO_TUNER
#endif//DISABLE_AUTO_TUNING


#ifdef  SHOW_NI_DEPENDENCE
      int Ni_active = 0;
      for(int ii = 0; ii < grpNum; ii++)
	Ni_active += laneInfo_hst[ii].num;
      printf("%d\t%d\t%le\n", grpNum, Ni_active, execTime[steps - bench_begin].calcGravity_dev);
      fflush(stdout);
    }/* for(grpNum = Ngrp; grpNum >= NGROUPS; grpNum = (int)((double)grpNum * M_SQRT1_2)){ */
    steps = bench_begin + BENCHMARK_STEPS - 1;/**< kill other benchmark */
#endif//SHOW_NI_DEPENDENCE


    rebuild.interval += 1.0;
#ifndef DISABLE_AUTO_TUNING
#   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
    rebuild.avg += reduce;
    rebuild.var += reduce * reduce;
#endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)

    linearModel(&(rebuildParam.linearGuess), &(rebuildParam.linearStats), rebuild.interval, elapsed.walkTree[rebuild.reuse] * reduce, reduce);
    powerModel (&(rebuildParam. powerGuess), &(rebuildParam. powerStats), rebuild.interval, elapsed.walkTree[rebuild.reuse] * reduce, reduce);
#ifdef  USE_PARABOLIC_GROWTH_MODEL
    parabolicModel(&(rebuildParam.parabolicGuess), &(rebuildParam.parabolicStats), rebuild.interval, elapsed.walkTree[rebuild.reuse] * reduce, reduce);
#ifdef  USE_ADDITIONAL_SWITCH
     useParabolicGuess |= (elapsed.walkTree[rebuild.reuse] < elapsed.makeTree);
#endif//USE_ADDITIONAL_SWITCH
#endif//USE_PARABOLIC_GROWTH_MODEL
#endif//DISABLE_AUTO_TUNING

#ifdef  DISABLE_AUTO_TUNING
    rebuild.reuse = rebuild.interval < rebuild_interval;
#else///DISABLE_AUTO_TUNING
    rebuild.reuse = 1;
#ifdef  USE_PARABOLIC_GROWTH_MODEL
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
    }/* if( rebuild.interval > 3.5 ){ */
#else///USE_PARABOLIC_GROWTH_MODEL
    if( (rebuild.interval > 3.5) &&
	((rebuildParam.linearGuess.rchisq < 1.0e-30 + rebuildParam.powerGuess.rchisq) ? (rebuildParam.linearGuess.time) : (rebuildParam.powerGuess.time)) > (elapsed.makeTree * reduce) )
	rebuild.reuse = 0;
#endif//USE_PARABOLIC_GROWTH_MODEL
#endif//DISABLE_AUTO_TUNING


    /** orbit integration (correct step) */
#ifdef  BLOCK_TIME_STEP
    correction_dev
      (grpNum, laneInfo_dev, laneTime_dev, eps, eta, ibody0_dev, rebuild.reuse
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
#else///BLOCK_TIME_STEP
#ifndef LEAP_FROG_INTEGRATOR
    /** correct velocity of i-particles */
    advVel_dev
      (Ni, ibody0_dev, HALF * (real)dt
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
#else///LEAP_FROG_INTEGRATOR
    /** update velocity of i-particles */
    advVel_dev
      (Ni, ibody0_dev, (real)dt
#ifdef  EXEC_BENCHMARK
       , &execTime[steps - bench_begin]
#endif//EXEC_BENCHMARK
       );
#endif//LEAP_FROG_INTEGRATOR
#endif//BLOCK_TIME_STEP


#ifdef  COUNT_INTERACTIONS
    copyCounters_dev2hst(Ni, treeinfo_dev, treeinfo);
    analyzeTreeMetrics(&treeProp[steps - bench_begin]);
    analyzeWalkStatistics(Ni, grpNum, Nj  , &(treeProp[steps - bench_begin].Nj  ));
    analyzeWalkStatistics(Ni, grpNum, Nbuf, &(treeProp[steps - bench_begin].Nbuf));
    for(int ii = 0; ii < Ni; ii++)
      treeProp[steps - bench_begin].Ninteractions += (ulong)Nj[ii];
    treeProp[steps - bench_begin].bufSize = soaWalk_dev.bufSize;
#endif//COUNT_INTERACTIONS


#ifndef EXEC_BENCHMARK
    /** output time evolution of numerical results */
#ifndef BLOCK_TIME_STEP
    present = (uint)(time * invSnapshotInterval);
#endif//BLOCK_TIME_STEP

    /** save tentative results to restart the simulation */
    static double currentTime;
#ifdef  SERIALIZED_EXECUTION
    currentTime = getPresentTimeInMin();
#else///SERIALIZED_EXECUTION
    if( mpi.rank == 0 )
      currentTime = getPresentTimeInMin();
    chkMPIerr(MPI_Bcast(&currentTime, 1, MPI_DOUBLE, 0, mpi.comm));
#endif//SERIALIZED_EXECUTION

    if( currentTime > formerTime + saveInterval ){
#ifndef COMPARE_WITH_DIRECT_SOLVER
#ifndef SERIALIZED_EXECUTION
      MPI_Barrier(mpi.comm);
#endif//SERIALIZED_EXECUTION
      static struct timespec timeDump;
      clock_gettime(CLOCK_MONOTONIC_RAW, &timeDump);
      double dumpElapsed = calcElapsedTimeInSec(timeInit, timeDump);
#ifndef SERIALIZED_EXECUTION
      chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &dumpElapsed, 1, MPI_DOUBLE, MPI_MAX, mpi.comm));
#endif//SERIALIZED_EXECUTION
      dumpElapsed += prevElapsed;
      dumpRestarter
	(num, ibody0_dev, ibody0, time, dt, steps, dumpElapsed, file, &last, &formerTime
#ifdef  USE_HDF5_FORMAT
	 , hdf5type, rebuild, elapsed, rebuildParam, brentDistance, brentHistory
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
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
	 , deviceMonitors, &monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
	 );
#endif//COMPARE_WITH_DIRECT_SOLVER
      rebuild.reuse = 0;
#ifdef  BLOCK_TIME_STEP
      rebuild.adjust = false;
#endif//BLOCK_TIME_STEP
    }/* if( currentTime > formerTime + saveInterval ){ */

    /** output snapshot file to analyze and visualize the results */
    if( present != previous ){
#ifdef  COMPARE_WITH_DIRECT_SOLVER
      ibody_direct_dev     = ibody0_dev;
      ibody_direct_dev.acc = direct_dev;
#ifdef  GADGET_MAC
      ibody_direct_dev.acc_old = direct_old_dev;
#endif//GADGET_MAC
#endif//COMPARE_WITH_DIRECT_SOLVER

#ifdef  REPORT_COMPUTE_RATE
      brent_rate += (brentDistance.u.pos / brentDistance.b) * (double)(steps - brent_prev);
      brent_prev = steps;
#endif//REPORT_COMPUTE_RATE
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
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
	 , deviceMonitors, &monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  REPORT_COMPUTE_RATE
	 , &clock_prev, &time_prev, &step_prev, ft, &brent_rate, &brent_prev, &rebuild_freq
#endif//REPORT_COMPUTE_RATE
#ifdef  COMPARE_WITH_DIRECT_SOLVER
	 , Ni, ibody_direct_dev, devProp, direct, Ngrp, laneInfo_dev
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	 , eps * eps
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
	 , pot_tbl_sphe
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
	 , pot_tbl_disk
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
	 , accfile
#endif//COMPARE_WITH_DIRECT_SOLVER
	 );

      rebuild.reuse = 0;
#ifdef  BLOCK_TIME_STEP
      rebuild.adjust = false;
#endif//BLOCK_TIME_STEP
    }/* if( present != previous ){ */

#else///EXEC_BENCHMARK
#       ifdef  BLOCK_TIME_STEP
    previous = present;
#       endif//BLOCK_TIME_STEP
#endif//EXEC_BENCHMARK
  }/* while( time < ft ){ */


#ifdef  COMPARE_WITH_DIRECT_SOLVER
  {
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

    currentNum++;
    fp = fopen(filename, "w");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    fprintf(fp, "%d\n", currentNum);
    fclose(fp);

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
  }
#endif//COMPARE_WITH_DIRECT_SOLVER


#   if  !defined(EXEC_BENCHMARK) && !defined(COMPARE_WITH_DIRECT_SOLVER)
  /** output final stage of numerical results */
#ifndef SERIALIZED_EXECUTION
      MPI_Barrier(mpi.comm);
#endif//SERIALIZED_EXECUTION
      static struct timespec timeDump;
      clock_gettime(CLOCK_MONOTONIC_RAW, &timeDump);
      double dumpElapsed = calcElapsedTimeInSec(timeInit, timeDump);
#ifndef SERIALIZED_EXECUTION
      chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &dumpElapsed, 1, MPI_DOUBLE, MPI_MAX, mpi.comm));
#endif//SERIALIZED_EXECUTION
      dumpElapsed += prevElapsed;
  dumpRestarter
    (num, ibody0_dev, ibody0, time, dt, steps, dumpElapsed, file, &last, &formerTime
#ifdef  USE_HDF5_FORMAT
     , hdf5type, rebuild, elapsed, rebuildParam, brentDistance, brentHistory
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
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
     , deviceMonitors, &monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
     );
#endif//!defined(EXEC_BENCHMARK) && !defined(COMPARE_WITH_DIRECT_SOLVER)


#ifndef SERIALIZED_EXECUTION

#ifdef  MPI_ONE_SIDED_FOR_LET
  chkMPIerr(MPI_Win_free(&(letcfg.win_more)));
  chkMPIerr(MPI_Win_free(&(letcfg.win_jpos)));
  chkMPIerr(MPI_Win_free(&(letcfg.win_mass)));
#endif//MPI_ONE_SIDED_FOR_LET

#ifdef  MPI_ONE_SIDED_FOR_EXCG

  chkMPIerr(MPI_Win_free(&(samplePos0.win_x)));  chkMPIerr(MPI_Win_free(&(samplePos1.win_x)));
  chkMPIerr(MPI_Win_free(&(samplePos0.win_y)));  chkMPIerr(MPI_Win_free(&(samplePos1.win_y)));
  chkMPIerr(MPI_Win_free(&(samplePos0.win_z)));  chkMPIerr(MPI_Win_free(&(samplePos1.win_z)));

#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_free(&(ibody0_dev.win_ipos)));
  chkMPIerr(MPI_Win_free(&(ibody1_dev.win_ipos)));
#ifdef  GADGET_MAC
  chkMPIerr(MPI_Win_free(&(ibody0_dev.win_iacc)));
  chkMPIerr(MPI_Win_free(&(ibody1_dev.win_iacc)));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_free(&(ibody0_dev.win_ivel)));
  chkMPIerr(MPI_Win_free(&(ibody1_dev.win_ivel)));
  chkMPIerr(MPI_Win_free(&(ibody0_dev.win_time)));
  chkMPIerr(MPI_Win_free(&(ibody1_dev.win_time)));
#else///BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_free(&(ibody0_dev.win_vx)));
  chkMPIerr(MPI_Win_free(&(ibody1_dev.win_vx)));
  chkMPIerr(MPI_Win_free(&(ibody0_dev.win_vy)));
  chkMPIerr(MPI_Win_free(&(ibody1_dev.win_vy)));
  chkMPIerr(MPI_Win_free(&(ibody0_dev.win_vz)));
  chkMPIerr(MPI_Win_free(&(ibody1_dev.win_vz)));
#endif//BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_free(&(ibody0_dev.win_idx)));
  chkMPIerr(MPI_Win_free(&(ibody1_dev.win_idx)));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_free(&(ibody0.win_ipos)));
  chkMPIerr(MPI_Win_free(&(ibody1.win_ipos)));
#ifdef  GADGET_MAC
  chkMPIerr(MPI_Win_free(&(ibody0.win_iacc)));
  chkMPIerr(MPI_Win_free(&(ibody1.win_iacc)));
#endif//GADGET_MAC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  chkMPIerr(MPI_Win_free(&(ibody0.win_iext)));
  chkMPIerr(MPI_Win_free(&(ibody1.win_iext)));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_free(&(ibody0.win_ivel)));
  chkMPIerr(MPI_Win_free(&(ibody1.win_ivel)));
  chkMPIerr(MPI_Win_free(&(ibody0.win_time)));
  chkMPIerr(MPI_Win_free(&(ibody1.win_time)));
#else///BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_free(&(ibody0.win_vx)));
  chkMPIerr(MPI_Win_free(&(ibody1.win_vx)));
  chkMPIerr(MPI_Win_free(&(ibody0.win_vy)));
  chkMPIerr(MPI_Win_free(&(ibody1.win_vy)));
  chkMPIerr(MPI_Win_free(&(ibody0.win_vz)));
  chkMPIerr(MPI_Win_free(&(ibody1.win_vz)));
#endif//BLOCK_TIME_STEP
  chkMPIerr(MPI_Win_free(&(ibody0.win_idx)));
  chkMPIerr(MPI_Win_free(&(ibody1.win_idx)));
#endif//MPI_VIA_HOST
#endif//MPI_ONE_SIDED_FOR_EXCG

#endif//SERIALIZED_EXECUTION


  /** memory deallocation */
#ifdef  USE_HDF5_FORMAT
  freeSnapshotArray
    (hdf5_pos, hdf5_vel, hdf5_acc, hdf5_m, hdf5_pot, hdf5_idx
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , hdf5_acc_ext, hdf5_pot_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
#endif//USE_HDF5_FORMAT

  freePeanoHilbertKey_dev
    (tag_dev, peano_dev, peano, min_dev, max_dev, gsync_ph0, gsync_ph1
#ifndef SERIALIZED_EXECUTION
     , box_min, box_max, min_hst, max_hst
#endif//SERIALIZED_EXECUTION
#ifdef  CUB_AVAILABLE
     , phsort_temp_storage, tag_pre, peano_pre
#endif//CUB_AVAILABLE
     );

#ifdef  BLOCK_TIME_STEP
  freeParticleDataSoA_hst
    (idx0, pos0, acc0, vel0, time0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext0
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  freeParticleDataSoA_hst
    (idx1, pos1, acc1, vel1, time1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext1
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  freeParticleDataSoA_dev
    (idx0_dev, pos0_dev, acc0_dev, vel0_dev, time0_dev,
     idx1_dev, pos1_dev, acc1_dev, vel1_dev, time1_dev
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext_dev
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     , neighbor_dev
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
     , encBall_dev, encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     , box_min_hst, box_max_hst, icom_hst
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
     , tmp_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     , res_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     );
#ifdef  SWITCH_WITH_J_PARALLELIZATION
#ifdef  MPI_VIA_HOST
  freeParticleDataSoA_hst
    (idx_dist0, pos_dist0, acc_dist0, vel_dist0, time_dist0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext_dist0
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  freeParticleDataSoA_hst
    (idx_dist1, pos_dist1, acc_dist1, vel_dist1, time_dist1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext_dist1
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
#endif//MPI_VIA_HOST
  freeParticleDataSoA_dev
    (idx_dist0_dev, pos_dist0_dev, acc_dist0_dev, vel_dist0_dev, time_dist0_dev,
     idx_dist1_dev, pos_dist1_dev, acc_dist1_dev, vel_dist1_dev, time_dist1_dev
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext_dist_dev
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     , neighbor_dist_dev
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
     , encBall_dist_dev, encBall_dist_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     , box_min_dist_hst, box_max_dist_hst, icom_dist_hst
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
     , tmp_dist_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     , res_dist_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     );
#endif//SWITCH_WITH_J_PARALLELIZATION
#else///BLOCK_TIME_STEP
  freeParticleDataSoA_hst
    (idx0, pos0, acc0, vx0, vy0, vz0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext0
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  freeParticleDataSoA_hst
    (idx1, pos1, acc1, vx1, vy1, vz1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext1
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  freeParticleDataSoA_dev
    (idx0_dev, pos0_dev, acc0_dev, vx0_dev, vy0_dev, vz0_dev
     , idx1_dev, pos1_dev, acc1_dev, vx1_dev, vy1_dev, vz1_dev
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext_dev
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     , neighbor_dev
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
     , encBall_dev, encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     , box_min_hst, box_max_hst, icom_hst
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
     , tmp_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     , res_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     );
#ifdef  SWITCH_WITH_J_PARALLELIZATION
#ifdef  MPI_VIA_HOST
  freeParticleDataSoA_hst
    (idx_dist0, pos_dist0, acc_dist0, vx_dist0, vy_dist0, vz_dist0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext_dist0
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  freeParticleDataSoA_hst
    (idx_dist1, pos_dist1, acc_dist1, vx_dist1, vy_dist1, vz_dist1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext_dist1
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
#endif//MPI_VIA_HOST
  freeParticleDataSoA_dev
    (idx_dist0_dev, pos_dist0_dev, acc_dist0_dev, vx_dist0_dev, vy_dist0_dev, vz_dist0_dev
     , idx1_dev, pos1_dev, acc1_dev, vx1_dev, vy1_dev, vz1_dev
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , acc_ext_dist_dev
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     , neighbor_dist_dev
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
     , encBall_dist_dev, encBall_dist_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     , box_min_dist_hst, box_max_dist_hst, icom_dist_hst
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
     , tmp_dist_dev
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     , res_dist_dev
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
     );
#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//BLOCK_TIME_STEP

#ifdef  COUNT_INTERACTIONS
  freeParticleInfoSoA_hst(Nj, Nbuf);
#endif//COUNT_INTERACTIONS
  freeParticleInfoSoA_dev
    (jtag_dev
#ifdef  COUNT_INTERACTIONS
     , Nj_dev, Nbuf_dev
#endif//COUNT_INTERACTIONS
     );

  freeParticleGroups(laneInfo_hst, laneInfo_dev, laneTime_dev, inum_hst, inum_dev
#ifdef  SWITCH_WITH_J_PARALLELIZATION
		     , true
#endif//SWITCH_WITH_J_PARALLELIZATION
		     );
#ifdef  SWITCH_WITH_J_PARALLELIZATION
  freeParticleGroups(laneInfo_ext_hst, laneInfo_ext_dev, NULL, NULL, NULL, false);
#endif//SWITCH_WITH_J_PARALLELIZATION

  freeTreeCell_dev
    (cell_dev, leaf_dev, node_dev, list_dev,
     hkey_dev, parent_dev, children_dev, bottomLev_dev, numCell_dev, numNode_dev, scanNum_dev
#ifdef  COUNT_INTERACTIONS
     , cell, leaf, node, list
#endif//COUNT_INTERACTIONS
     );

  freeTreeNode_dev
    (more_dev, pj_dev, mj_dev, bmax_dev, node2cell_dev, gsync0, gsync1,
#       ifdef  WS93_MAC
     mr2_dev,
#       endif//WS93_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
     more, pj, mj,
#endif//!defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
     gmem_make_tree_dev, gsync0_make_tree_dev, gsync1_make_tree_dev, gsync2_make_tree_dev, gsync3_make_tree_dev,
     gmem_link_tree_dev, gsync0_link_tree_dev, gsync1_link_tree_dev,
#ifdef  GADGET_MAC
     mac_dev,
#endif//GADGET_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
     gmem_external_dev, gsync0_external_dev, gsync1_external_dev, diameter_dev, diameter_hst,
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
     more0Buf, more1Buf, rjmaxBuf, makeFail);

  freeTreeBuffer_dev
    (fail_dev, buffer, freeLst
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
     , freeNum, active
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
     , cycles_hst, cycles_dev
#ifndef SERIALIZED_EXECUTION
     , cycles_dist_hst, cycles_dist_dev
#endif//SERIALIZED_EXECUTION
#endif//defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
     , cycles_let_hst, cycles_let_dev
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
     );

#ifdef  COMPARE_WITH_DIRECT_SOLVER
  freeAccel_dev
    (direct_dev, direct
#ifdef  GADGET_MAC
     , direct_old_dev
#endif//GADGET_MAC
     );
#endif//COMPARE_WITH_DIRECT_SOLVER

#ifndef BLOCK_TIME_STEP
  freeTimeStep_dev(dt_dev);
#endif//BLOCK_TIME_STEP

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  freeSphericalPotentialTable_dev
    (pot_tbl_sphe_Phi
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     , pot_tbl_sphe_rad
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     );
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
  freeDiskPotentialTable_dev
    (pot_tbl_disk_Phi,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     pot_tbl_disk_RR, pot_tbl_disk_zz
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     pot_tbl_disk_FRz
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     );
  freeSphericalPotentialTable_dev
    (pot_tbl_disk_sphe_Phi
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     , pot_tbl_disk_sphe_rad
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     );
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#ifndef SERIALIZED_EXECUTION
  releaseLETtopology
    (nodeInfo,
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
     min, max,
#endif//USE_RECTANGULAR_BOX_FOR_LET
     ipos,
#ifdef  GADGET_MAC
     amin,
#endif//GADGET_MAC
     numSend_hst, numSend_dev, stream_let, Nstream_let);

#   if  defined(USE_ENCLOSING_BALL_FOR_LET) && !defined(OCTREE_BASED_SEARCH)
  freeApproxEnclosingBall_dev(appEncBall_dev, gsync0_EB, gsync1_EB);
#endif//defined(USE_ENCLOSING_BALL_FOR_LET) && !defined(OCTREE_BASED_SEARCH)

  releaseORMtopology
    (dxmin, dxmax, dymin, dymax, dzmin, dzmax, dmreq,
     sxmin, sxmax, symin, symax, szmin, szmax,
     iparticleSendBuf, iparticleRecvBuf, sampleRecvNum, sampleRecvDsp,
     ormCfg, repCfg, letcfg.rank);
  releaseSamplePos
    (x0hst, x1hst, y0hst, y1hst, z0hst, z1hst, idhst,
     x0dev, x1dev, y0dev, y1dev, z0dev, z1dev, iddev);
  releaseParticlePosition
    (xhst, yhst, zhst, xdev, ydev, zdev, rank_hst, rank_dev, idx_dev);
  releaseDomainPos
    (xmin_dev, xmax_dev, ymin_dev, ymax_dev, zmin_dev, zmax_dev,
     xmin_hst, xmax_hst, ymin_hst, ymax_hst, zmin_hst, zmax_hst,
     domrank_dev, domrank_hst, numNew_dev, numNew_hst, gmem_dom, gsync0_dom, gsync1_dom);

  freeGeometricEnclosingBall_dev(r2geo_dev);

#ifdef  SWITCH_WITH_J_PARALLELIZATION
  releaseExternalBodyInfo(Ni_list, head_list, grpNum_list, grpNum_disp);
#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//SERIALIZED_EXECUTION

  /** destroy CUDA streams */
  for(int ii = 0; ii < sinfo.num; ii++)
    mycudaStreamDestroy(stream[ii]);
  free(stream);

#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
#endif//USE_HDF5_FORMAT


#ifdef  REPORT_TOTAL_ELAPSED_TIME
#ifndef SERIALIZED_EXECUTION
  MPI_Barrier(mpi.comm);
#endif//SERIALIZED_EXECUTION
  static struct timespec timeExit;
  clock_gettime(CLOCK_MONOTONIC_RAW, &timeExit);
  double totalElapsed = calcElapsedTimeInSec(timeInit, timeExit);
#ifndef SERIALIZED_EXECUTION
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &totalElapsed, 1, MPI_DOUBLE, MPI_MAX, mpi.comm));
  if( mpi.rank == 0 )
#endif//SERIALIZED_EXECUTION
    {
      totalElapsed += prevElapsed;
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
    }
#endif//REPORT_TOTAL_ELAPSED_TIME


#ifndef SERIALIZED_EXECUTION
  exitMPI();
#endif//SERIALIZED_EXECUTION


  return (0);
}


#ifdef  COUNT_INTERACTIONS
/**
 * @fn analyzeWalkStatistics
 *
 * @brief Analyze statistics about tree traversal.
 *
 * @param (Ni) number of i-particles
 * @param (Ngrp) number of particle groups
 * @param (results) measured results
 * @return (stats) analyzed results
 */
void analyzeWalkStatistics(int Ni, const int Ngrp, int * restrict results, walk_stats *stats)
{
  __NOTE__("%s\n", "start");


  stats->Ngroup = Ngrp;
  stats->min    = INT_MAX;
  stats->max    = 0;
  stats->mean   = 0.0f;
  stats->sdev   = 0.0f;

  for(int ii = 0; ii < Ni; ii++){
    const int dat = results[ii];

    if( dat != 0 ){
      if( stats->min > dat )	stats->min = dat;
      if( stats->max < dat )	stats->max = dat;
    }/* if( dat != 0 ){ */

    stats->mean += dat;
    stats->sdev += dat * dat;
  }/* for(int ii = 0; ii < Ni; ii++){ */

  if( (Ngrp * (TSUB / NWARP)) < Ni )
    Ni = Ngrp * (TSUB / NWARP);
  const float inv = 1.0f / (float)Ni;
  stats->mean *= inv;
  stats->sdev *= inv;
  stats->sdev -= stats->mean * stats->mean;
  stats->sdev  = sqrtf(stats->sdev);


  __NOTE__("%s\n", "end");
}


/**
 * @fn analyzeTreeStatistics
 *
 * @brief Analyze statistics about tree structure.
 *
 * @return (stats) analyzed results
 */
static inline void analyzeTreeStatistics(tree_stats *stats)
{
  const real inv = UNITY / (real)stats->nodeNum;

  stats->mjMean *= inv;  stats->r2Mean *= inv;  stats->mjSdev -= stats->mjMean * stats->mjMean;  stats->mjSdev = SQRT(stats->mjSdev);
  stats->mjSdev *= inv;  stats->r2Sdev *= inv;  stats->r2Sdev -= stats->r2Mean * stats->r2Mean;  stats->r2Sdev = SQRT(stats->r2Sdev);
}


/**
 * @fn analyzeTreeMetrics
 *
 * @brief Analyze metrics about tree structure.
 *
 * @return (metrics) analyzed results
 *
 * @sa analyzeTreeStatistics
 */
void analyzeTreeMetrics(tree_metrics *metric)
{
  __NOTE__("%s\n", "start");


  for(uint ii = 0; ii < MAXIMUM_PHKEY_LEVEL; ii++)
    if( metric->level[ii].nodeNum != 0 ){
      metric->total.nodeNum += metric->level[ii].nodeNum;
      metric->total.cellNum += metric->level[ii].cellNum;
      metric->total.mjMean  += metric->level[ii].mjMean;
      metric->total.mjSdev  += metric->level[ii].mjSdev;
      metric->total.r2Mean  += metric->level[ii].r2Mean;
      metric->total.r2Sdev  += metric->level[ii].r2Sdev;

      analyzeTreeStatistics(&metric->level[ii]);
    }/* if( metric->level[ii].nodeNum != 0 ){ */

  analyzeTreeStatistics(&metric->total);


  __NOTE__("%s\n", "end");
}
#endif//COUNT_INTERACTIONS
