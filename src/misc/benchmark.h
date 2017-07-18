/**
 * @file benchmark.h
 *
 * @brief Header file for benchmark of GOTHIC
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/06/27 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef BENCHMARK_H
#define BENCHMARK_H
#ifdef EXEC_BENCHMARK


#include "macro.h"


/**
 * @def BENCHMARK_STEPS
 *
 * @brief Number of time steps in benchmark
 */
/* #define BENCHMARK_STEPS (8) */
/* #define BENCHMARK_STEPS (16) */
/* #define BENCHMARK_STEPS (32) */
/* #define BENCHMARK_STEPS (64) */
/* #define BENCHMARK_STEPS (128) */
/* #define BENCHMARK_STEPS (256) */
/* #define BENCHMARK_STEPS (512) */
/* #define BENCHMARK_STEPS (1024) */
#define BENCHMARK_STEPS (2048)
/* #define BENCHMARK_STEPS (4096) */
/* #define BENCHMARK_STEPS (8192) */


/* macro to specify the file name */
#define WALLCLOCK "time"


/**
 * @struct wall_clock_time
 *
 * @brief structure to store results of the measurements
 */
typedef struct
{
  double calcGravity_dev;
  double calcMultipole;
  double makeTree;
#ifdef  BLOCK_TIME_STEP
  double prediction_dev;
  double correction_dev;
  double setLaneTime_dev;
  double adjustParticleTime_dev;
#else///BLOCK_TIME_STEP
  double advPos_dev;
  double advVel_dev;
#endif//BLOCK_TIME_STEP
  double setTimeStep_dev;
  double sortParticlesPHcurve;
  double copyParticle_dev2hst;
  double copyParticle_hst2dev;
  double setTreeNode_dev;
  double setTreeCell_dev;
  double examineNeighbor_dev;
  double  searchNeighbor_dev;
#ifdef  HUNT_MAKE_PARAMETER
  double genPHkey_kernel;
  double rsortKey_library;
  double sortBody_kernel;
  double makeTree_kernel;
  double linkTree_kernel;
  double trimTree_kernel;
  double initTreeLink_kernel;
  double initTreeCell_kernel;
  double initTreeNode_kernel;
  double initTreeBody_kernel;
  double copyRealBody_kernel;
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
  double searchNeighbor_kernel;
  double sortNeighbors;
  double countNeighbors_kernel;
  double commitNeighbors;
#endif//HUNT_FIND_PARAMETER
} wall_clock_time;


/**
 * @struct tree_stats
 *
 * @brief structure to store statistics of the tree cells and nodes
 */
typedef struct
{
  real mjMean, mjSdev;
  real r2Mean, r2Sdev;
  int cellNum, nodeNum;
} tree_stats;
/**
 * @struct walk_stats
 *
 * @brief structure to store statistics of the tree traversal
 */
typedef struct
{
  int Ngroup, min, max;
  float mean, sdev;
} walk_stats;
/**
 * @struct tree_stats
 *
 * @brief structure to store statistics of the tree structure
 */
typedef struct
{
  tree_stats level[21];/* 21 is the possible maximum of MAXIMUM_PHKEY_LEVEL */
  tree_stats total;
  walk_stats Nbuf;
  walk_stats Nj;
  ulong Ninteractions;/* ULONG_MAX is ~10^20, would be enough size */
  size_t bufSize;
} tree_metrics;


#include <time.h>
#include "timer.h"
static struct timespec _benchIni, _benchFin;

/**
 * @fn initStopwatch
 *
 * @brief Initialize benchmarking.
 */
static inline void initStopwatch(void)
{
#ifdef  __CUDACC__
  cudaDeviceSynchronize();
#endif//__CUDACC__
  clock_gettime(CLOCK_MONOTONIC_RAW, &_benchIni);
}
/**
 * @fn stopStopwatch
 *
 * @brief Finalize benchmarking.
 *
 * @return (result) measured execution time in units of second
 */
static inline void stopStopwatch(double *result)
{
#ifdef  __CUDACC__
  cudaDeviceSynchronize();
#endif//__CUDACC__
  clock_gettime(CLOCK_MONOTONIC_RAW, &_benchFin);
  *result += calcElapsedTimeInSec(_benchIni, _benchFin);
}


#endif//EXEC_BENCHMARK
#endif//BENCHMARK_H
