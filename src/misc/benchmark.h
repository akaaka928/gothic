/**
 * @file benchmark.h
 *
 * @brief Header file for benchmark of GOTHIC
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
/* #define BENCHMARK_STEPS (2048) */
/* #define BENCHMARK_STEPS (4096) */
/* #define BENCHMARK_STEPS (8192) */
#define BENCHMARK_STEPS (16384)


/* macro to specify the file name */
#define WALLCLOCK "time"


/**
 * @struct wall_clock_time
 *
 * @brief structure to store results of the measurements
 */
typedef struct
{
  double calcGravity_dev;/**< execution time of calcGravity_dev() in src/tree/walk_dev.cu */
  double calcMultipole_dev;/**< execution time of calcMultipole_dev() in src/tree/make_dev.cu */
  double makeTree;/**< execution time of makeTreeStructure_dev() in src/tree/make_dev.cu */
#ifdef  BLOCK_TIME_STEP
  double prediction_dev;/**< execution time of prediction_dev() in src/time/adv_dev.cu */
  double correction_dev;/**< execution time of correction_dev() in src/time/adv_dev.cu */
  double setLaneTime_dev;/**< execution time of setLaneTime_dev() in src/time/adv_dev.cu */
  double adjustParticleTime_dev;/**< execution time of adjustParticleTime_dev() in src/time/adv_dev.cu */
#else///BLOCK_TIME_STEP
  double advPos_dev;/**< execution time of advPos_dev() in src/time/adv_dev.cu */
  double advVel_dev;/**< execution time of advVel_dev() in src/time/adv_dev.cu */
#endif//BLOCK_TIME_STEP
  double setTimeStep_dev;/**< execution time of setTimeStep_dev() in src/time/adv_dev.cu */
  double sortParticlesPHcurve;/**< execution time of sortParticlesPHcurve_dev() in src/sort/peano_dev.cu */
  double copyParticle_dev2hst;/**< execution time of copyParticle_dev2hst() in src/time/adv_dev.cu */
  double copyParticle_hst2dev;/**< execution time of copyParticle_hst2dev() in src/time/adv_dev.cu */
  double examineNeighbor_dev;/**< execution time of examineParticleSeparation() in src/tree/shrink_dev.cu */
  double  searchNeighbor_dev;/**< execution time of updateParticleGroups() in src/tree/shrink_dev.cu */
#ifdef  HUNT_MAKE_PARAMETER
  double genPHkey_kernel;/**< execution time of calcPHkey_kernel() in src/sort/peano_dev.cu */
  double rsortKey_library;/**< execution time of Peano--Hilbert key sorting in src/sort/peano_dev.cu */
  double sortBody_kernel;/**< execution time of sortParticlesPHcurve_kernel() in src/sort/peano_dev.cu */
  double makeTree_kernel;/**< execution time of makeTree_kernel() in src/tree/make_dev.cu */
  double linkTree_kernel;/**< execution time of linkTree_kernel() in src/tree/make_dev.cu */
  double trimTree_kernel;/**< execution time of trimTree_kernel() in src/tree/make_dev.cu */
  double initTreeLink_kernel;/**< execution time of initTreeLink_kernel() in src/tree/make_dev.cu */
  double initTreeCell_kernel;/**< execution time of initTreeCell_kernel() in src/tree/make_dev.cu */
  double initTreeNode_kernel;/**< execution time of initTreeNode_kernel() in src/tree/make_dev.cu */
  double initTreeBody_kernel;/**< execution time of initTreeBody_kernel() in src/tree/make_dev.cu */
  double copyRealBody_kernel;/**< execution time of copyRealBody_kernel() in src/tree/make_dev.cu */
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
  double searchNeighbor_kernel;/**< execution time of facileNeighborSearching_dev() in src/tree/neighbor_dev.cu */
  double sortNeighbors;/**< execution time of reduction about rmax in src/tree/shrink_dev.cu */
  double countNeighbors_kernel;/**< execution time of countContinuousNeighbor_kernel() in src/tree/shrink_dev.cu */
  double commitNeighbors;/**< execution time of while() loop in updateParticleGroups in src/tree/shrink_dev.cu */
#endif//HUNT_FIND_PARAMETER
#ifndef SERIALIZED_EXECUTION
  double calcAcc_kernel;/**< execution time of calcAcc_kernel() in src/tree/walk_dev.cu; measured by using clock64() function and clock frequency */
  double exchangeLET;/**< elapsed time for MPI communications in calcGravity_dev() in src/tree/walk_dev.cu */
#ifdef  MONITOR_LETGEN_TIME
  double makeLET_kernel;/**< execution time of makeLET_kernel() in src/tree/let_dev.cu; measured by using clock64() function and clock frequency */
#endif//MONITOR_LETGEN_TIME
#ifdef  USE_ENCLOSING_BALL_FOR_LET
  double encBall_dev;/**< execution time of getApproxEnclosingBall_dev() in src/tree/icom_dev.cu */
#endif//USE_ENCLOSING_BALL_FOR_LET

  double calc_r2max_dev;/**< execution time of calc_r2max_dev() in src/tree/geo_dev.cu */




  double excgBody_dev;/**< execution time of exchangeParticles_dev() in src/para/exchange_dev.cu */


#if 1
  なるべく process ごとの値を生で保存しておいて，後から load balancing も取り出せるようにしたい;

  HUNT_????_PARAMETER を設定して，各種パラメータを見つけるためのジョブも準備したい;

  各関数の処理時間を正確に測定しようとすると，auto-tuning によって実行パラメータが変化していくので，breakdown を正確に測定する際には必要最低限の同期処理をつっこむだけにしておきたい;
#endif

#endif//SERIALIZED_EXECUTION
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
