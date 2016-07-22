/*************************************************************************\
 *                                                                       *
                  last updated on 2016/03/11(Fri) 14:19:18
 *                                                                       *
 *    Header File for orbit integration of N-body simulation             *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef ADV_DEV_H
#define ADV_DEV_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef _SYS_TIME_H
#      include <sys/time.h>
#endif//_SYS_TIME_H
//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef BENCHMARK_H
#       include "../misc/benchmark.h"
#endif//BENCHMARK_H
//-------------------------------------------------------------------------
#ifndef STRUCTURE_H
#       include "../misc/structure.h"
#endif//STRUCTURE_H
//-------------------------------------------------------------------------
#ifndef WALK_DEV_H
#       include "../tree/walk_dev.h"
#endif//WALK_DEV_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef HUNT_TIME_PARAMETER
/* the below macro is enabled in the default option; switched off in the parameter survey mode to use -D and -U from Makefile */
#define USE_WARP_SHUFFLE_FUNC_TIME
#endif//HUNT_TIME_PARAMETER
//-------------------------------------------------------------------------
#   if  defined(USE_WARP_SHUFFLE_FUNC_TIME) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_TIME
#endif//defined(USE_WARP_SHUFFLE_FUNC_TIME) && (GPUGEN < 30)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef NTHREADS_TIME
#          if  GPUGEN >= 52
#define NTHREADS_TIME (512)
#       else///GPUGEN >= 52
#          if  GPUGEN >= 30
#define NTHREADS_TIME (128)
#       else///GPUGEN >= 30
#define NTHREADS_TIME (512)
#       endif//GPUGEN >= 30
#       endif//GPUGEN >= 52
#endif//NTHREADS_TIME
//-------------------------------------------------------------------------
/* NTHREADS_TIME must be equal or greater than NTHREADS, to set massless particles */
#   if  NTHREADS_TIME < NTHREADS
#undef  NTHREADS_TIME
#define NTHREADS_TIME   NTHREADS
#endif//NTHREADS_TIME
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
//-------------------------------------------------------------------------
#if  !defined(MPI_INCLUDED) && !defined(OMPI_MPI_H)
#       include <mpi.h>
#endif
//-------------------------------------------------------------------------
#ifndef MPILIB_H
#       include <mpilib.h>
#endif//MPILIB_H
//-------------------------------------------------------------------------
#ifndef MPICFG_H
#       include "../para/mpicfg.h"
#endif//MPICFG_H
//-------------------------------------------------------------------------
#endif//SERIALIZED_EXECUTION
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #ifdef  CUB_AVAILABLE */
/* typedef struct */
/* { */
/*   laneinfo *info; */
/*   double *time; */
/*   void *temp_storage; */
/*   size_t temp_storage_size; */
/* } soaCUBtime; */
/* #endif//CUB_AVAILABLE */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "adv_dev.cu"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
#ifndef BLOCK_TIME_STEP
  muse allocTimeStep_dev(real **dt_dev);
  void  freeTimeStep_dev(real  *dt_dev);
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
/* #   if  defined(CUB_AVAILABLE) && defined(BLOCK_TIME_STEP) */
/*   muse allocTimeStep_dev(void **temp_storage, laneinfo **info, double **time, soaCUBtime *buf, laneinfo *laneInfo_dev, double *laneTime_dev, const int Ngrp); */
/*   void  freeTimeStep_dev(void  *temp_storage, laneinfo  *info, double  *time); */
/* #endif//defined(CUB_AVAILABLE) && defined(BLOCK_TIME_STEP) */
  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
  void setTimeStep_dev
  (const int Ngrp, laneinfo * RESTRICT laneInfo_dev, double * RESTRICT laneTime_dev, int *grpNum, const iparticle pi,
   const double told, double *tnew, double *dt, bool adjustAllTimeStep, const double invSnapshotInterval, const uint previous, uint *present
/* #ifdef  CUB_AVAILABLE */
/*    , soaCUBtime buf */
/* #endif//CUB_AVAILABLE */
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
  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
/*   void prediction_dev(const int Nj, const real tnew, const iparticle pi */
/* #ifndef CALC_MULTIPOLE_ON_DEVICE */
/* 		      , const iparticle pi_hst */
/* #endif//CALC_MULTIPOLE_ON_DEVICE */
/* #ifdef  EXEC_BENCHMARK */
/* 		      , wall_clock_time *elapsed */
/* #endif//EXEC_BENCHMARK */
/* 		      ); */
  void prediction_dev(const int Nj, const double tnew, const iparticle pi
#ifndef CALC_MULTIPOLE_ON_DEVICE
		      , const iparticle pi_hst
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef  EXEC_BENCHMARK
		      , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		      );
/*   void correction_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, real * RESTRICT laneTime, const real eps, const real eta, const iparticle pi, const int reuseTree */
/* #ifdef  EXEC_BENCHMARK */
/* 		      , wall_clock_time *elapsed */
/* #endif//EXEC_BENCHMARK */
/* 		      ); */
  void correction_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, const real eps, const real eta, const iparticle pi, const int reuseTree
#ifdef  EXEC_BENCHMARK
		      , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		      );
/*   void adjustParticleTime_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, real * RESTRICT laneTime, const real eps, const real eta, const iparticle pi */
/* #ifdef  EXEC_BENCHMARK */
/* 			      , wall_clock_time *elapsed */
/* #endif//EXEC_BENCHMARK */
/* 			      ); */
  void adjustParticleTime_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, const real eps, const real eta, const iparticle pi
#ifdef  EXEC_BENCHMARK
			      , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
			      );
/*   void setLaneTime_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, real * RESTRICT laneTime, const iparticle pi */
/* #ifdef  EXEC_BENCHMARK */
/* 		       , wall_clock_time *elapsed */
/* #endif//EXEC_BENCHMARK */
/* 		       ); */
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
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
#ifdef  COMPARE_WITH_DIRECT_SOLVER
  void copyAccel_dev2hst(const int Ni, acceleration * RESTRICT dev, acceleration * RESTRICT hst);
#endif//COMPARE_WITH_DIRECT_SOLVER
  //-----------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
  void copyCounters_dev2hst(const int Ni, iparticle_treeinfo dev, iparticle_treeinfo hst);
#endif//COUNT_INTERACTIONS
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//ADV_DEV_H
//-------------------------------------------------------------------------
