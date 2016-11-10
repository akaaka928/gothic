/*************************************************************************\
 *                                                                       *
                  last updated on 2016/11/01(Tue) 10:46:01
 *                                                                       *
 *    Header File for tree traversal based on octree structure           *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef WALK_DEV_H
#define WALK_DEV_H
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/device.h"
//-------------------------------------------------------------------------
#include "../tree/macutil.h"
#include "../tree/make.h"
#include "../tree/buf_inc.h"
//-------------------------------------------------------------------------
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#ifndef SERIALIZED_EXECUTION
#include "../para/mpicfg.h"
#include "../tree/let.h"
#endif//SERIALIZED_EXECUTION
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
//-------------------------------------------------------------------------
#ifdef  COMPARE_WITH_DIRECT_SOLVER
#       include <stdbool.h>
#endif//COMPARE_WITH_DIRECT_SOLVER
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define PRINT_PSEUDO_PARTICLE_INFO */
//-------------------------------------------------------------------------
#   if  GPUGEN >= 52
/* #define ADOPT_SMALLEST_ENCLOSING_BALL */
#define ADOPT_APPROXIMATED_ENCLOSING_BALL
/* #define COMPARE_ENCLOSING_BALLS */
#endif//GPUGEN >= 52
//-------------------------------------------------------------------------
#define ADOPT_ENCLOSING_BALL
//-------------------------------------------------------------------------
#   if  defined(ADOPT_APPROXIMATED_ENCLOSING_BALL) && defined(ADOPT_SMALLEST_ENCLOSING_BALL)
#       undef   ADOPT_APPROXIMATED_ENCLOSING_BALL
#endif//defined(ADOPT_APPROXIMATED_ENCLOSING_BALL) && defined(ADOPT_SMALLEST_ENCLOSING_BALL)
//-------------------------------------------------------------------------
#   if  defined(COMPARE_ENCLOSING_BALLS) && (defined(ADOPT_SMALLEST_ENCLOSING_BALL) || defined(ADOPT_APPROXIMATED_ENCLOSING_BALL))
#       undef   COMPARE_ENCLOSING_BALLS
#endif//defined(COMPARE_ENCLOSING_BALLS) && (defined(ADOPT_SMALLEST_ENCLOSING_BALL) || defined(ADOPT_APPROXIMATED_ENCLOSING_BALL))
//-------------------------------------------------------------------------
#   if  !defined(ADOPT_ENCLOSING_BALL) && (defined(ADOPT_SMALLEST_ENCLOSING_BALL) || defined(ADOPT_APPROXIMATED_ENCLOSING_BALL) || defined(COMPARE_ENCLOSING_BALLS))
#        define  ADOPT_ENCLOSING_BALL
#endif//!defined(ADOPT_ENCLOSING_BALL) && (defined(ADOPT_SMALLEST_ENCLOSING_BALL) || defined(ADOPT_APPROXIMATED_ENCLOSING_BALL) || defined(COMPARE_ENCLOSING_BALLS))
//-------------------------------------------------------------------------
#ifndef HUNT_WALK_PARAMETER
/* the below macro is enabled in the default option; switched off in the parameter survey mode to use -D and -U from Makefile */
#define USE_WARP_SHUFFLE_FUNC
#endif//HUNT_WALK_PARAMETER
//-------------------------------------------------------------------------
#   if  defined(USE_WARP_SHUFFLE_FUNC) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC
#endif//defined(USE_WARP_SHUFFLE_FUNC) && (GPUGEN < 30)
//-------------------------------------------------------------------------
#ifdef  IJ_PARALLELIZATION
/* NWARP must be 1, 2, 4, 8, 16, or 32 */
#ifndef NWARP
#          if  GPUGEN >= 52
#define NWARP (4)
#       else///GPUGEN >= 52
/* #          if  GPUGEN >= 30 */
/* #define NWARP (2) */
/* #       else///GPUGEN >= 30 */
#define NWARP (1)
/* #       endif//GPUGEN >= 30 */
#       endif//GPUGEN >= 52
#endif//NWARP
#else///IJ_PARALLELIZATION
/* NWARP must be unity */
#ifndef NWARP
#define NWARP (1)
#endif//NWARP
#endif//IJ_PARALLELIZATION
//-------------------------------------------------------------------------
#ifndef NTHREADS
#          if  GPUGEN >= 30
#define NTHREADS (512)
#       else///GPUGEN >= 30
#define NTHREADS (256)
#       endif//GPUGEN >= 30
#endif//NTHREADS
//-------------------------------------------------------------------------
/* NTHREADS must be equal or smaller than 1024 (limitation comes from reduction defined in ../tree/geo_dev.cu) */
/* NTHREADS must be equal or smaller than 512 (limitation comes from NLOOP, NQUEUE defined in ../tree/walk_dev.h) */
#   if  NTHREADS > 512
#undef  NTHREADS
#define NTHREADS  (512)
#endif//NTHREADS > 512
//-------------------------------------------------------------------------
#ifdef  DIV_NTHREADS
#undef  DIV_NTHREADS
#endif//DIV_NTHREADS
#   if  NTHREADS == 1024
#define DIV_NTHREADS(a) ((a) >> 10)
#endif//NTHREADS == 1024
#   if  NTHREADS ==  512
#define DIV_NTHREADS(a) ((a) >>  9)
#endif//NTHREADS ==  512
#   if  NTHREADS ==  256
#define DIV_NTHREADS(a) ((a) >>  8)
#endif//NTHREADS ==  256
#   if  NTHREADS ==  128
#define DIV_NTHREADS(a) ((a) >>  7)
#endif//NTHREADS ==  128
//-------------------------------------------------------------------------
/* TSUB must be 1, 2, 4, 8, 16, 32 */
#ifndef TSUB
#define TSUB (32)
/* #define TSUB (16) */
#endif//TSUB
/* TSUB should be equal or greater than 16 for the performance issue */
#   if  TSUB < 16
#undef  TSUB
#define TSUB  (16)
#endif//TSUB < 16
/* TSUB must be equal or greater than NWARP */
#   if  TSUB < NWARP
#undef  TSUB
#define TSUB   NWARP
#endif//TSUB > NWARP
/* TSUB must be equal or smaller than NTHREADS */
#   if  TSUB > NTHREADS
#undef  TSUB
#define TSUB   NTHREADS
#endif//TSUB > NTHREADS
//-------------------------------------------------------------------------
#ifdef  DIV_TSUB
#undef  DIV_TSUB
#endif//DIV_TSUB
#   if  TSUB == 32
#define DIV_TSUB(a) ((a) >> 5)
#endif//TSUB == 32
#   if  TSUB == 16
#define DIV_TSUB(a) ((a) >> 4)
#endif//TSUB == 16
#   if  TSUB ==  8
#define DIV_TSUB(a) ((a) >> 3)
#endif//TSUB ==  8
#   if  TSUB ==  4
#define DIV_TSUB(a) ((a) >> 2)
#endif//TSUB ==  4
#   if  TSUB ==  2
#define DIV_TSUB(a) ((a) >> 1)
#endif//TSUB ==  2
#   if  TSUB ==  1
#define DIV_TSUB(a) ((a))
#endif//TSUB ==  1
//-------------------------------------------------------------------------
#ifdef  DIV_NWARP
#undef  DIV_NWARP
#endif//DIV_NWARP
#   if  NWARP == 32
#define DIV_NWARP(a) ((a) >> 5)
#endif//NWARP == 32
#   if  NWARP == 16
#define DIV_NWARP(a) ((a) >> 4)
#endif//NWARP == 16
#   if  NWARP ==  8
#define DIV_NWARP(a) ((a) >> 3)
#endif//NWARP ==  8
#   if  NWARP ==  4
#define DIV_NWARP(a) ((a) >> 2)
#endif//NWARP ==  4
#   if  NWARP ==  2
#define DIV_NWARP(a) ((a) >> 1)
#endif//NWARP ==  2
#   if  NWARP ==  1
#define DIV_NWARP(a) ((a))
#endif//NWARP ==  1
//-------------------------------------------------------------------------
#define NGROUPS (DIV_TSUB(NTHREADS))
//-------------------------------------------------------------------------
#ifndef NLOOP
/* #define NLOOP (6) */
/* #define NLOOP (5) */
/* #define NLOOP (4) */
/* #define NLOOP (3) */
/* #define NLOOP (2) */
#define NLOOP (1)
#endif//NLOOP
//-------------------------------------------------------------------------
/* determined by capacity of the shared memory (48KB per SM --> 24KB per block) */
#ifndef NBLOCKS_PER_SM
#define NBLOCKS_PER_SM (2)
#endif//NBLOCKS_PER_SM
//-------------------------------------------------------------------------
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
//-------------------------------------------------------------------------
#   if  NBLOCKS_PER_SM == 2
/* 6144 = 6 * 1024 = 12 * 1024 / 2 = (48KiB / sizeof(float)) / NBLOCKS_PER_SM */
#   if  NLOOP > ((6144 / (5 * NTHREADS)) - 2)
#undef  NLOOP
#define NLOOP   ((6144 / (5 * NTHREADS)) - 2)
#endif//NLOOP > ((6144 / (5 * NTHREADS)) - 2)
#       ifdef  USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(6144) - 5 * (1 + NLOOP))
#       else///USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(6144) - 5 * (1 + NLOOP) - 1)
#       endif//USE_WARP_SHUFFLE_FUNC
#else///NBLOCKS_PER_SM == 2
/* 12288 = 12 * 1024 = 48KiB / sizeof(float) */
#   if  NLOOP > ((6144 / (5 * NTHREADS)) - 2)
#undef  NLOOP
#define NLOOP   ((6144 / (5 * NTHREADS)) - 2)
#endif//NLOOP > ((6144 / (5 * NTHREADS)) - 2)
#       ifdef  USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(6144) - 5 * (1 + NLOOP))
#       else///USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(6144) - 5 * (1 + NLOOP) - 1)
#       endif//USE_WARP_SHUFFLE_FUNC
#endif//NBLOCKS_PER_SM == 2
//-------------------------------------------------------------------------
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
//-------------------------------------------------------------------------
#   if  NBLOCKS_PER_SM == 2
/* 6144 = 6 * 1024 = 12 * 1024 / 2 = (48KiB / sizeof(float)) / NBLOCKS_PER_SM */
#   if  NLOOP > ((6144 / (4 * NTHREADS)) - 2)
#undef  NLOOP
#define NLOOP   ((6144 / (4 * NTHREADS)) - 2)
#endif//NLOOP > ((6144 / (4 * NTHREADS)) - 2)
#       ifdef  USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(6144) - 4 * (1 + NLOOP))
#       else///USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(6144) - 4 * (1 + NLOOP) - 1)
#       endif//USE_WARP_SHUFFLE_FUNC
#else///NBLOCKS_PER_SM == 2
/* 12288 = 12 * 1024 = 48KiB / sizeof(float) */
#   if  NLOOP > ((12288 / (4 * NTHREADS * NBLOCKS_PER_SM)) - 2)
#undef  NLOOP
#define NLOOP   ((12288 / (4 * NTHREADS * NBLOCKS_PER_SM)) - 2)
#endif//NLOOP > ((12288 / (4 * NTHREADS * NBLOCKS_PER_SM)) - 2)
#       ifdef  USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(12288 / NBLOCKS_PER_SM) - 4 * (1 + NLOOP))
#       else///USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(12288 / NBLOCKS_PER_SM) - 4 * (1 + NLOOP) - 1)
#       endif//USE_WARP_SHUFFLE_FUNC
#endif//NBLOCKS_PER_SM == 2
//-------------------------------------------------------------------------
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
//-------------------------------------------------------------------------
/* NSTOCK must be a power of two */
/* #define NSTOCK (1) */
/* #define NSTOCK (2) */
#define NSTOCK (4)
//-------------------------------------------------------------------------
#define IDX_SHIFT_BITS (8)
#define IDX_SHIFT_MASK (0xff)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef CUDALIB_H
typedef struct __align__(16){  uint x, y, z, w;} uint4;
#endif//CUDALIB_H
//-------------------------------------------------------------------------
typedef struct
{
  cudaStream_t *stream;
  int idx, num;
} kernelStream;
//-------------------------------------------------------------------------
#ifdef  ADOPT_SMALLEST_ENCLOSING_BALL
#define  NDIM_SEB (3)
#define NDIM2_SEB (9)
#   if  (NDIM_SEB * (NDIM_SEB + 1)) > TSUB
#undef  NDIM_SEB
#undef NDIM2_SEB
#endif//(NDIM_SEB * (NDIM_SEB + 1)) > TSUB
typedef struct __align__(16)
{
  real x, y, z;
  int  support;
} pos4seb;
typedef union __align__(16)
{
  pos4seb pos;
  real4   r4;
  real    val[NDIM_SEB + 1];
  int     idx[NDIM_SEB + 1];
} dat4seb;
#endif//ADOPT_SMALLEST_ENCLOSING_BALL
typedef union __align__(16)
{
  jparticle pos;
  position  pi;
  acceleration ai;
  uint idx[NSTOCK];
  real val[4];
#ifdef  ADOPT_SMALLEST_ENCLOSING_BALL
  dat4seb seb;
#endif//ADOPT_SMALLEST_ENCLOSING_BALL
} jnode;
//-------------------------------------------------------------------------
typedef union
{
  uint i;
  real r;
} uint_real;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "walk_dev.cu"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  muse setCUDAstreams_dev(cudaStream_t **stream, kernelStream *sinfo, deviceInfo *info, deviceProp *prop
/* #   if  defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)) */
/* 			  , cudaEvent_t **iniEvent, cudaEvent_t **finEvent */
/* #endif//defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)) */
			  );
  //-----------------------------------------------------------------------
#   if  defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
  muse allocateCUDAevents_dev
  (cudaEvent_t **iniWalk, cudaEvent_t **finWalk
#ifdef  MONITOR_LETGEN_TIME
   , cudaEvent_t **iniMake, cudaEvent_t **finMake
#endif//MONITOR_LETGEN_TIME
   , const int Ngpu);
  void  releaseCUDAevents_dev
  (cudaEvent_t  *iniWalk, cudaEvent_t  *finWalk
#ifdef  MONITOR_LETGEN_TIME
   , cudaEvent_t  *iniMake, cudaEvent_t  *finMake
#endif//MONITOR_LETGEN_TIME
   , const int Ngpu);
#endif//defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
  //-----------------------------------------------------------------------
  void  freeTreeBuffer_dev
  (int  *failure, uint  *buffer, uint  *freeLst
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
   , uint  *freeNum, int  *active
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#ifndef USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
   , unsigned long long int  *cycles_hst, unsigned long long int  *cycles_dev
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
   , unsigned long long int  *cycles_let_hst, unsigned long long int  *cycles_let_dev
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#endif//USE_CUDA_EVENT
   );
  muse allocTreeBuffer_dev
  (int **failure, uint **buffer, uint **freeLst,
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
   uint **freeNum, int **active,
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#ifndef USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
   unsigned long long int **cycles_hst, unsigned long long int **cycles_dev,
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
   unsigned long long int **cycles_let_hst, unsigned long long int **cycles_let_dev,
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#endif//USE_CUDA_EVENT
   soaTreeWalkBuf *buf, const int num_max, const muse used, const deviceProp gpu);
  //-----------------------------------------------------------------------
  void calcGravity_dev
  (const int grpNum
#ifdef  BLOCK_TIME_STEP
   , double *reduce, const int totNum
#endif//BLOCK_TIME_STEP
   , laneinfo * RESTRICT laneInfo, const int Ni, const iparticle pi, const soaTreeNode tree, const soaTreeWalkBuf buf
   , kernelStream *sinfo, deviceProp devProp, double *time
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
   , char *file
#endif//PRINT_PSEUDO_PARTICLE_INFO
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  USE_CUDA_EVENT
   , cudaEvent_t *iniCalcAcc, cudaEvent_t *finCalcAcc
#else///USE_CUDA_EVENT
   , unsigned long long int *cycles_hst, unsigned long long int *cycles_dev
#endif//USE_CUDA_EVENT
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#ifndef SERIALIZED_EXECUTION
   , double *twalk, const int pjNum
#ifdef  LET_COMMUNICATION_VIA_HOST
   , const soaTreeNode tree_hst
#endif//LET_COMMUNICATION_VIA_HOST
   , const int Nlet, domainInfo *let, const int Nstream_let, cudaStream_t stream_let[], MPIcfg_tree mpi
#ifdef  MONITOR_LETGEN_TIME
   , double *tlet
#ifdef  USE_CUDA_EVENT
   , cudaEvent_t *iniMakeLET, cudaEvent_t *finMakeLET
#else///USE_CUDA_EVENT
   , unsigned long long int *cycles_let_hst, unsigned long long int *cycles_let_dev
#endif//USE_CUDA_EVENT
#endif//MONITOR_LETGEN_TIME
#endif//SERIALIZED_EXECUTION
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#ifdef  COUNT_INTERACTIONS
   , iparticle_treeinfo treeInfo
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
#ifdef  COMPARE_WITH_DIRECT_SOLVER
   , const bool approxGravity
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
   , const real eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#endif//COMPARE_WITH_DIRECT_SOLVER
   );
  //-----------------------------------------------------------------------
  void setGlobalConstants_walk_dev_cu(const real newton_hst, const real eps2_hst
#ifndef WS93_MAC
				      , const real theta2_hst
#endif//WS93_MAC
				      );
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//WALK_DEV_H
//-------------------------------------------------------------------------
