/**
 * @file adv_dev.h
 *
 * @brief Header file for orbit integration of N-body particles
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/10/26 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef ADV_DEV_H
#define ADV_DEV_H


#include <sys/time.h>

#include "macro.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../tree/walk_dev.h"

#ifndef SERIALIZED_EXECUTION
#include <mpi.h>
#include "mpilib.h"
#include "../para/mpicfg.h"
#endif//SERIALIZED_EXECUTION


#ifndef HUNT_TIME_PARAMETER
/* the below macro is enabled in the default option; switched off in the parameter survey mode to use -D and -U from Makefile */
#define USE_WARP_SHUFFLE_FUNC_TIME
#endif//HUNT_TIME_PARAMETER

#   if  defined(USE_WARP_SHUFFLE_FUNC_TIME) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_TIME
#endif//defined(USE_WARP_SHUFFLE_FUNC_TIME) && (GPUGEN < 30)


/**
 * @def NTHREADS_TIME
 *
 * @brief number of threads per block for setTimeStep_kernel, prediction_kernel, advPos_kernel, and advVel_kernel
 */
#ifndef NTHREADS_TIME
#   if  GPUGEN >= 60
#define NTHREADS_TIME (128)
#else///GPUGEN >= 60
#   if  GPUGEN >= 52
#define NTHREADS_TIME (512)
#else///GPUGEN >= 52
#   if  GPUGEN >= 30
#define NTHREADS_TIME (128)
#else///GPUGEN >= 30
#define NTHREADS_TIME (512)
#endif//GPUGEN >= 30
#endif//GPUGEN >= 52
#endif//GPUGEN >= 60
#endif//NTHREADS_TIME

/** NTHREADS_TIME must be equal or greater than NTHREADS, to set massless particles */
#   if  NTHREADS_TIME < NTHREADS
#undef  NTHREADS_TIME
#define NTHREADS_TIME   NTHREADS
#endif//NTHREADS_TIME


/* list of functions appeared in ``adv_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
#ifndef BLOCK_TIME_STEP
  muse allocTimeStep_dev(real **dt_dev);
  void  freeTimeStep_dev(real  *dt_dev);
#endif//BLOCK_TIME_STEP

#ifdef  BLOCK_TIME_STEP
  void setTimeStep_dev
  (const int Ngrp, laneinfo * RESTRICT laneInfo_dev, double * RESTRICT laneTime_dev, int *grpNum, const iparticle pi,
   const double told, double *tnew, double *dt, bool adjustAllTimeStep, const double invSnapshotInterval, const uint previous, uint *present
#ifndef SERIALIZED_EXECUTION
   , MPIcfg_tree mpi
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
   );
#else///BLOCK_TIME_STEP
  void setTimeStep_dev(const int Ni, iparticle ibody, const real eta, const real eps, real *dt_dev, double *dt_hst
#ifndef SERIALIZED_EXECUTION
		       , MPIcfg_tree mpi
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
		       , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		       );
#endif//BLOCK_TIME_STEP

#ifdef  BLOCK_TIME_STEP
  void prediction_dev(const int Nj, const double tnew, const iparticle pi
#ifdef  EXEC_BENCHMARK
		      , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		      );
  void correction_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, const real eps, const real eta, const iparticle pi, const int reuseTree
#ifdef  EXEC_BENCHMARK
		      , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		      );
  void adjustParticleTime_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, const real eps, const real eta, const iparticle pi
#ifdef  EXEC_BENCHMARK
			      , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
			      );
  void setLaneTime_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, const iparticle pi
#ifdef  EXEC_BENCHMARK
		       , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		       );
#else///BLOCK_TIME_STEP
  void advPos_dev(const int Ni, iparticle ibody, const real dt
#ifdef  EXEC_BENCHMARK
		  , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		  );
  void advVel_dev(const int Ni, iparticle ibody, const real dt
#ifdef  EXEC_BENCHMARK
		  , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		  );
#endif//BLOCK_TIME_STEP

  void copyParticle_hst2dev(const int Ni, iparticle hst, iparticle dev
#ifdef  EXEC_BENCHMARK
			    , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
			    );
  void copyParticle_dev2hst(const int Ni, iparticle dev, iparticle hst
#ifdef  EXEC_BENCHMARK
			    , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
			    );

  void copyParticleAsync_hst2dev(const int Ni, iparticle hst, iparticle dev, cudaStream_t stream);
  void copyParticleAsync_dev2hst(const int Ni, iparticle dev, iparticle hst, cudaStream_t stream);

#ifdef  COMPARE_WITH_DIRECT_SOLVER
  void copyAccel_dev2hst(const int Ni, acceleration * RESTRICT dev, acceleration * RESTRICT hst);
#endif//COMPARE_WITH_DIRECT_SOLVER

#ifdef  COUNT_INTERACTIONS
  void copyCounters_dev2hst(const int Ni, iparticle_treeinfo dev, iparticle_treeinfo hst);
#endif//COUNT_INTERACTIONS

#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//ADV_DEV_H
