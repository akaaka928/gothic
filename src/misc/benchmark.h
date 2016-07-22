/*************************************************************************\
 *                                                                       *
                  last updated on 2016/03/09(Wed) 13:51:11
 *                                                                       *
 *    Header File for Benchmark of N-body Simulation                     *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef BENCHMARK_H
#define BENCHMARK_H
//-------------------------------------------------------------------------
#ifdef EXEC_BENCHMARK
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define BENCHMARK_STEPS (8) */
/* #define BENCHMARK_STEPS (16) */
/* #define BENCHMARK_STEPS (32) */
/* #define BENCHMARK_STEPS (64) */
/* #define BENCHMARK_STEPS (128) */
/* #define BENCHMARK_STEPS (256) */
/* #define BENCHMARK_STEPS (512) */
#define BENCHMARK_STEPS (1024)
/* #define BENCHMARK_STEPS (2048) */
/* #define BENCHMARK_STEPS (4096) */
/* #define BENCHMARK_STEPS (8192) */
//-------------------------------------------------------------------------
#define WALLCLOCK "time"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
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
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
  double examineNeighbor_dev;
  double  searchNeighbor_dev;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
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
//-------------------------------------------------------------------------
typedef struct
{
  real mjMean, mjSdev;
  real r2Mean, r2Sdev;
  int cellNum, nodeNum;
} tree_stats;
typedef struct
{
  int Ngroup, min, max;
  float mean, sdev;
} walk_stats;
typedef struct
{
  tree_stats level[21];/* 21 is the possible maximum of MAXIMUM_PHKEY_LEVEL */
  tree_stats total;
  walk_stats Nbuf;
  walk_stats Nj;
  ulong Ninteractions;/* ULONG_MAX is ~10^20, would be enough size */
  int bufSize;
} tree_metrics;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef _SYS_TIME_H
#      include <sys/time.h>
#endif//_SYS_TIME_H
//-------------------------------------------------------------------------
#ifndef TIMER_H
#      include <timer.h>
#endif//TIMER_H
//-------------------------------------------------------------------------
static struct timeval _benchIni, _benchFin;
//-------------------------------------------------------------------------
static inline void initStopwatch(void)
{
#ifdef  __CUDACC__
  cudaDeviceSynchronize();
#endif//__CUDACC__
  gettimeofday(&_benchIni, NULL);
}
//-------------------------------------------------------------------------
static inline void stopStopwatch(double *result)
{
#ifdef  __CUDACC__
  cudaDeviceSynchronize();
#endif//__CUDACC__
  gettimeofday(&_benchFin, NULL);
  *result += calcElapsedTimeInSec(_benchIni, _benchFin);
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//EXEC_BENCHMARK
//-------------------------------------------------------------------------
#endif//BENCHMARK_H
//-------------------------------------------------------------------------
