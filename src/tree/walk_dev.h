/**
 * @file walk_dev.h
 *
 * @brief Header file for tree traversal based on octree structure on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2021/04/23 (Fri)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef WALK_DEV_H
#define WALK_DEV_H


#ifdef  COMPARE_WITH_DIRECT_SOLVER
#include <stdbool.h>
#endif//COMPARE_WITH_DIRECT_SOLVER

#include "macro.h"
#include "cudalib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/device.h"

#include "../tree/make.h"
#include "../tree/buf_inc.h"

#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#ifndef SERIALIZED_EXECUTION
#include "../misc/tune.h"
#include "../para/mpicfg.h"
#include "../tree/let.h"
#ifdef  SWITCH_WITH_J_PARALLELIZATION

/**
 * @def NMAX_J_PARALLELIZATION
 *
 * @brief maximum number of activated i-particles to employ j-parallelization
 */
#define NMAX_J_PARALLELIZATION NUM_BODY_MAX

/**
 * @def FCRIT_J_PARALLELIZATION
 *
 * @brief consider to exploit j-parallelization when grpNum < FCRIT_J_PARALLELIZATION * totNum
 */
/* #define FCRIT_J_PARALLELIZATION (2.5e-1f) */
/* #define FCRIT_J_PARALLELIZATION (1.25e-1f) */
/* #define FCRIT_J_PARALLELIZATION (6.25e-2f) */
/* #define FCRIT_J_PARALLELIZATION (3.125e-2f) */
/* #define FCRIT_J_PARALLELIZATION (1.5625e-2f) */
/* #define FCRIT_J_PARALLELIZATION (7.8125e-3f) */
#define FCRIT_J_PARALLELIZATION (3.90625e-3f)
/* #define FCRIT_J_PARALLELIZATION (1.953125e-3f) */
/* #define FCRIT_J_PARALLELIZATION (9.765625e-4f) */

/**
 * @def FSIZE_J_PARALLELIZATION
 *
 * @brief consider to exploit j-parallelization when more than FSIZE_J_PARALLELIZATION * mpi.size processes prefer to exploit j-parallelization
 */
#define FSIZE_J_PARALLELIZATION (0.5f)

#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//SERIALIZED_EXECUTION
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)


/**
 * @def ADOPT_MIKI2013_PERFORMANCE_MODEL
 *
 * @brief activates performance model based on Miki et al. (2013, Computer Physics Communications, 184, 2159-2168)
 */
/* #define ADOPT_MIKI2013_PERFORMANCE_MODEL */


/**
 * @def PARTIAL_SUM_ACCELERATION
 *
 * @brief activates accumulation with partial sum to reduce floating-point number operation errors
 */
#define PARTIAL_SUM_ACCELERATION

/**
 * @def ACCURATE_ACCUMULATION
 *
 * @brief activates Kahan summation to accumulate partial sums to reduce floating-point number operation errors
 */
/**
 * @def ACCURATE_PARTIAL_SUM
 *
 * @brief activates Kahan summation to calculate partial sums to reduce floating-point number operation errors; this options is much slower than the original implementation
 */
#ifdef  PARTIAL_SUM_ACCELERATION
#define ACCURATE_ACCUMULATION
/* #define ACCURATE_PARTIAL_SUM */
#endif//PARTIAL_SUM_ACCELERATION

#   if  defined(ACCURATE_PARTIAL_SUM) && defined(ACCURATE_ACCUMULATION)
#undef          ACCURATE_PARTIAL_SUM
#endif//defined(ACCURATE_PARTIAL_SUM) && defined(ACCURATE_ACCUMULATION)


/**
 * @def MERGE_QUEUED_TREE_NODES
 *
 * @brief activates merging to queued tree nodes to maximize capacity of buffers; this options is much slower than the original implementation
 */
/* #define MERGE_QUEUED_TREE_NODES */


/* #define PRINT_PSEUDO_PARTICLE_INFO */


/**
 * @def ADOPT_SMALLEST_ENCLOSING_BALL
 *
 * @brief uses the smallest enclosing ball (Fischer et al. 2003)
 */
/**
 * @def ADOPT_APPROXIMATED_ENCLOSING_BALL
 *
 * @brief uses the efficient bounding sphere (Ritter 1990)
 */
/**
 * @def ADOPT_ENCLOSING_BALL
 *
 * @brief uses the sphere centered on the geometric center of the enclosing rectangular cuboid
 */
/**
 * @def COMPARE_ENCLOSING_BALLS
 *
 * @brief uses the smaller ball of GEO and CMP
 * @detail see Appendix B in Miki & Umemura (2017), New Astronomy, 52, 65--81
 */


#   if  GPUVER >= 52
/* #define ADOPT_SMALLEST_ENCLOSING_BALL */
#define ADOPT_APPROXIMATED_ENCLOSING_BALL
/* #define COMPARE_ENCLOSING_BALLS */
#endif//GPUVER >= 52


#define ADOPT_ENCLOSING_BALL


#   if  defined(ADOPT_APPROXIMATED_ENCLOSING_BALL) && defined(ADOPT_SMALLEST_ENCLOSING_BALL)
#       undef   ADOPT_APPROXIMATED_ENCLOSING_BALL
#endif//defined(ADOPT_APPROXIMATED_ENCLOSING_BALL) && defined(ADOPT_SMALLEST_ENCLOSING_BALL)

#   if  defined(COMPARE_ENCLOSING_BALLS) && (defined(ADOPT_SMALLEST_ENCLOSING_BALL) || defined(ADOPT_APPROXIMATED_ENCLOSING_BALL))
#       undef   COMPARE_ENCLOSING_BALLS
#endif//defined(COMPARE_ENCLOSING_BALLS) && (defined(ADOPT_SMALLEST_ENCLOSING_BALL) || defined(ADOPT_APPROXIMATED_ENCLOSING_BALL))

#   if  !defined(ADOPT_ENCLOSING_BALL) && (defined(ADOPT_SMALLEST_ENCLOSING_BALL) || defined(ADOPT_APPROXIMATED_ENCLOSING_BALL) || defined(COMPARE_ENCLOSING_BALLS))
#        define  ADOPT_ENCLOSING_BALL
#endif//!defined(ADOPT_ENCLOSING_BALL) && (defined(ADOPT_SMALLEST_ENCLOSING_BALL) || defined(ADOPT_APPROXIMATED_ENCLOSING_BALL) || defined(COMPARE_ENCLOSING_BALLS))


#ifndef HUNT_WALK_PARAMETER
/* the below macro is enabled in the default option; switched off in the parameter survey mode to use -D and -U from Makefile */
#define USE_WARP_SHUFFLE_FUNC
#define USE_WARP_REDUCE_FUNC
#define USE_L2_SET_ASIDE_POLICY
#endif//HUNT_WALK_PARAMETER

#   if  defined(USE_WARP_SHUFFLE_FUNC) && (GPUVER < 30)
#undef          USE_WARP_SHUFFLE_FUNC
#endif//defined(USE_WARP_SHUFFLE_FUNC) && (GPUVER < 30)

#   if  defined(USE_WARP_REDUCE_FUNC) && !defined(ENABLE_WARP_REDUCE_FUNC)
#undef          USE_WARP_REDUCE_FUNC
#endif//defined(USE_WARP_REDUCE_FUNC) && !defined(ENABLE_WARP_REDUCE_FUNC)

#   if  defined(USE_L2_SET_ASIDE_POLICY) && !defined(ENABLE_L2_SET_ASIDE)
#undef          USE_L2_SET_ASIDE_POLICY
#endif//defined(USE_L2_SET_ASIDE_POLICY) && !defined(ENABLE_L2_SET_ASIDE)


#ifdef  USE_L2_SET_ASIDE_POLICY
/**
 * @def NLEVEL_TREE_NODE_L2_PERSISTING
 *
 * @brief number of levels for tree nodes persisting on L2 cache (40 MB)
 */
#ifndef NLEVEL_TREE_NODE_L2_PERSISTING
#define NLEVEL_TREE_NODE_L2_PERSISTING (6)
#endif//NLEVEL_TREE_NODE_L2_PERSISTING
/* more(4 byte) + jpos(16 byte) + mj(4 byte or 8 byte for w/o or w/ INDIVIDUAL_GRAVITATIONAL_SOFTENING) per tree node */
/* Lev = 0: Nnode = 1 -> 24 byte or 28 byte */
/* Lev = 1: Nnode = 1 + 8 = 9 -> 216 byte or 252 byte */
/* Lev = 2: Nnode = 9 + 64 = 73 -> 1752 byte or 2044 byte */
/* Lev = 3: Nnode = 73 + 512 = 585 -> ~14 KiB or ~16 KiB */
/* Lev = 4: Nnode = 585 + 4096 = 4681 -> ~110 KiB or ~128 KiB */
/* Lev = 5: Nnode = 4681 + 32768 = 37449 -> ~878 KiB or ~1 MiB */
/* Lev = 6: Nnode = 37449 + 262144 = 299593 -> ~6.9 MiB or ~8 MiB */
/** jpos(16 byte) per tree node */
/** Lev = 0: Nnode = 1 -> 16 byte */
/** Lev = 1: Nnode = 1 + 8 = 9 -> 144 byte */
/** Lev = 2: Nnode = 9 + 64 = 73 -> 1168 byte */
/** Lev = 3: Nnode = 73 + 512 = 585 -> 9360 byte ~9 KiB */
/** Lev = 4: Nnode = 585 + 4096 = 4681 -> 74896 byte ~73 KiB */
/** Lev = 5: Nnode = 4681 + 32768 = 37449 -> 599184 byte ~585 KiB */
/** Lev = 6: Nnode = 37449 + 262144 = 299593 -> 4793488 byte ~4.6 MiB */
/** Lev = 7: Nnode = 299593 + 2097152 = 2396745 -> 38347920 byte ~37 MiB */
#endif//USE_L2_SET_ASIDE_POLICY






/**
 * @def NWARP
 *
 * @brief number of threads share an i-particle
 * @detail must be 1, 2, 4, 8, 16, or 32 to exploit implicit synchronization within a warp
 */
#ifdef  IJ_PARALLELIZATION
#ifndef NWARP
#   if  GPUVER >= 80
#define NWARP (1)
#else///GPUVER >= 80
#   if  GPUVER >= 70
#define NWARP (2)
#else///GPUVER >= 70
#   if  GPUVER >= 60
#define NWARP (1)
#else///GPUVER >= 60
#   if  GPUVER >= 52
#define NWARP (4)
#else///GPUVER >= 52
#define NWARP (1)
#endif//GPUVER >= 52
#endif//GPUVER >= 60
#endif//GPUVER >= 70
#endif//GPUVER >= 80
#endif//NWARP
#else///IJ_PARALLELIZATION
/** NWARP must be unity */
#ifndef NWARP
#define NWARP (1)
#endif//NWARP
#endif//IJ_PARALLELIZATION


#ifdef  DPADD_FOR_ACC
#undef  NWARP
#define NWARP (1)
#endif//DPADD_FOR_ACC


/**
 * @def NBLOCKS_PER_SM
 *
 * @brief number of blocks per SM
 * @detail determined by capacity of the shared memory
 */
#ifndef NBLOCKS_PER_SM
#   if  GPUVER >= 80
#define NBLOCKS_PER_SM (3)
#else///GPUVER >= 80
#define NBLOCKS_PER_SM (2)
#endif//GPUVER >= 80
#endif//NBLOCKS_PER_SM


/**
 * @def NTHREADS
 *
 * @brief number of threads per block for calcAcc_kernel
 */
#ifndef NTHREADS
#   if  GPUVER >= 80
#define NTHREADS (256)
#else///GPUVER >= 80
#   if  GPUVER >= 30
#define NTHREADS (512)
#else///GPUVER >= 30
#define NTHREADS (256)
#endif//GPUVER >= 30
#endif//GPUVER >= 80
#endif//NTHREADS

/** NTHREADS must be equal or smaller than 1024 (limitation comes from reduction defined in ../tree/geo_dev.cu) */
#   if  GPUVER >= 70
/** NTHREADS * NBLOCKS_PER_SM <= 2048 for V100 */
#   if  NBLOCKS_PER_SM == 2
#   if  NTHREADS > 1024
#undef  NTHREADS
#define NTHREADS  (1024)
#endif//NTHREADS > 512
#else///NBLOCKS_PER_SM == 2
/** for NBLOCKS_PER_SM is 3 or 4 */
#   if  NTHREADS > 512
#undef  NTHREADS
#define NTHREADS  (512)
#endif//NTHREADS > 512
#endif//NBLOCKS_PER_SM == 2
#else///GPUVER >= 70
/** NTHREADS must be equal or smaller than 512 (limitation comes from NLOOP defined in ../tree/walk_dev.h) */
#   if  NBLOCKS_PER_SM == 2
#   if  NTHREADS > 512
#undef  NTHREADS
#define NTHREADS  (512)
#endif//NTHREADS > 512
#else///NBLOCKS_PER_SM == 2
/** for NBLOCKS_PER_SM is 3 or 4 */
#   if  NTHREADS > 256
#undef  NTHREADS
#define NTHREADS  (256)
#endif//NTHREADS > 256
#endif//NBLOCKS_PER_SM == 2
#endif//GPUVER >= 70

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


/**
 * @def TSUB
 *
 * @brief number of threads share operations
 * @detail must be 1, 2, 4, 8, 16, or 32 to exploit implicit synchronization within a warp
 */
#ifndef TSUB
#define TSUB (32)
#endif//TSUB

/** TSUB should be equal or greater than 16 for the performance issue */
#ifndef HUNT_WALK_PARAMETER
#   if  TSUB < 16
#undef  TSUB
#define TSUB  (16)
#endif//TSUB < 16
#endif//HUNT_WALK_PARAMETER
/** TSUB must be equal or greater than NWARP */
#   if  TSUB < NWARP
#undef  TSUB
#define TSUB   NWARP
#endif//TSUB > NWARP
/** TSUB must be equal or smaller than NTHREADS */
#   if  TSUB > NTHREADS
#undef  TSUB
#define TSUB   NTHREADS
#endif//TSUB > NTHREADS

#ifdef  DIV_TSUB
#undef  DIV_TSUB
#endif//DIV_TSUB
#   if  TSUB == 32
#define DIV_TSUB(a) ((a) >> 5)
#define SHFL_MASK_TSUB SHFL_MASK_32
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

#define NGROUPS (DIV_TSUB(NTHREADS))

/**
 * @def NLOOP
 *
 * @brief a parameter to increase arithmetic intensity
 */
#ifndef NLOOP
#   if  GPUVER >= 80
#define NLOOP (1)
#else///GPUVER >= 80
#   if  GPUVER >= 70
#define NLOOP (3)
#else///GPUVER >= 70
#define NLOOP (1)
#endif//GPUVER >= 70
#endif//GPUVER >= 80
#endif//NLOOP


#define SMEM_SIZE_CALC_ACC (SMEM_SIZE_SM_PREF / NBLOCKS_PER_SM)
#   if  SMEM_SIZE_CALC_ACC > MAX_SMEM_SIZE_BLOCK
#undef  SMEM_SIZE_CALC_ACC
#define SMEM_SIZE_CALC_ACC (MAX_SMEM_SIZE_BLOCK)
#endif//SMEM_SIZE_CALC_ACC > MAX_SMEM_SIZE_BLOCK


#   if  SMEM_SIZE_CALC_ACC > MAX_STATIC_SMEM_SIZE_BLOCK
/* use dynamic allocation for shared memory */
#define DYNAMIC_SMEM_SIZE_CALC_ACC SMEM_SIZE_CALC_ACC
#else///SMEM_SIZE_CALC_ACC > MAX_STATIC_SMEM_SIZE_BLOCK
/* use static allocation for shared memory */
#define DYNAMIC_SMEM_SIZE_CALC_ACC SMEM_SIZE
#endif//SMEM_SIZE_CALC_ACC > MAX_STATIC_SMEM_SIZE_BLOCK



/** capacity of shared memory per block in units of sizeof(float) */
#define NSM4TRAVERSAL (SMEM_SIZE_CALC_ACC >> 2)

#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING

#   if  NLOOP > ((NSM4TRAVERSAL / (5 * NTHREADS)) - 2)
#undef  NLOOP
#define NLOOP   ((NSM4TRAVERSAL / (5 * NTHREADS)) - 2)
#endif//NLOOP > ((NSM4TRAVERSAL / (5 * NTHREADS)) - 2)
#ifdef  USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(NSM4TRAVERSAL) - 5 * (1 + NLOOP))
#else///USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(NSM4TRAVERSAL) - 5 * (1 + NLOOP) - 1)
#endif//USE_WARP_SHUFFLE_FUNC

#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING

#   if  NLOOP > ((NSM4TRAVERSAL / (4 * NTHREADS)) - 2)
#undef  NLOOP
#define NLOOP   ((NSM4TRAVERSAL / (4 * NTHREADS)) - 2)
#endif//NLOOP > ((NSM4TRAVERSAL / (4 * NTHREADS)) - 2)
#ifdef  USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(NSM4TRAVERSAL) - 4 * (1 + NLOOP))
#else///USE_WARP_SHUFFLE_FUNC
#define NQUEUE  (DIV_NTHREADS(NSM4TRAVERSAL) - 4 * (1 + NLOOP) - 1)
#endif//USE_WARP_SHUFFLE_FUNC

#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING


/**
 * @def NSTOCK
 *
 * @brief a parameter to increase arithmetic intensity
 * @detail must be a power of two
 */
#define NSTOCK (4)
/* #define NSTOCK (2) */
/* #define NSTOCK (1) */


/**
 * @def IDX_SHIFT_BITS
 *
 * @brief a parameter to execute multiple prefix sum in a single operation
 * @detail IDX_SHIFT_BITS bits is used to memory a value
 */
#define IDX_SHIFT_BITS (8)

/**
 * @def IDX_SHIFT_MASK
 *
 * @brief a parameter to execute multiple prefix sum in a single operation
 * @detail determined by IDX_SHIFT_BITS
 */
#define IDX_SHIFT_MASK (0xff)


#ifndef CUDALIB_H
typedef struct __align__(16){  uint x, y, z, w;} uint4;
#endif//CUDALIB_H


/**
 * @struct kernelStream
 *
 * @brief structure for CUDA stream
 */
typedef struct
{
  cudaStream_t *stream;
  int idx, num;
} kernelStream;


#ifdef  ADOPT_SMALLEST_ENCLOSING_BALL

#define  NDIM_SEB (3)
#define NDIM2_SEB (9)
#   if  (NDIM_SEB * (NDIM_SEB + 1)) > TSUB
#undef   NDIM_SEB
#undef  NDIM2_SEB
#endif//(NDIM_SEB * (NDIM_SEB + 1)) > TSUB

/**
 * @struct pos4seb
 *
 * @brief structure for the smallest enclosing ball
 */
typedef struct __align__(16)
{
  real x, y, z;
  int  support;
} pos4seb;

/**
 * @struct dat4seb
 *
 * @brief structure for the smallest enclosing ball
 */
typedef union __align__(16)
{
  pos4seb pos;
  real4   r4;
  real    val[NDIM_SEB + 1];
  int     idx[NDIM_SEB + 1];
} dat4seb;
#endif//ADOPT_SMALLEST_ENCLOSING_BALL


/**
 * @union jnode
 *
 * @brief union for recycling registers in tree traversal
 */
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


/**
 * @union uint_real
 *
 * @brief union for recycling registers in tree traversal
 */
typedef union
{
  uint i;
  real r;
} uint_real;


/* list of functions appeared in ``walk_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  muse setCUDAstreams_dev(cudaStream_t **stream, kernelStream *sinfo, deviceInfo *info);

#ifdef  USE_L2_SET_ASIDE_POLICY
  void setL2persistent_walk_dev(kernelStream sinfo, const soaTreeNode tree);
#endif//USE_L2_SET_ASIDE_POLICY

  muse allocParticleDataSoA_dev
  (const int num
#ifdef  BLOCK_TIME_STEP
   , iparticle *body0, ulong **idx0, position **pos0, acceleration **acc0, velocity **vel0, ibody_time **ti0
   , iparticle *body1, ulong **idx1, position **pos1, acceleration **acc1, velocity **vel1, ibody_time **ti1
#else///BLOCK_TIME_STEP
   , iparticle *body0, ulong **idx0, position **pos0, acceleration **acc0, real **vx0, real **vy0, real **vz0
   , iparticle *body1, ulong **idx1, position **pos1, acceleration **acc1, real **vx1, real **vy1, real **vz1
#endif//BLOCK_TIME_STEP
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
   , acceleration **acc_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
   , real **neighbor
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
   , position **encBall, position **encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
   , position **box_min_hst, position **box_max_hst, position **icom_hst
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
   , DPacc **tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
   , acceleration **res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
   );
  void  freeParticleDataSoA_dev
  (ulong  *idx0, position  *pos0, acceleration  *acc0
#ifdef  BLOCK_TIME_STEP
   , velocity  *vel0, ibody_time  *ti0
   , ulong  *idx1, position  *pos1, acceleration  *acc1, velocity  *vel1, ibody_time  *ti1
#else///BLOCK_TIME_STEP
   , real  *vx0, real  *vy0, real  *vz0
   , ulong  *idx1, position  *pos1, acceleration  *acc1, real  *vx1, real  *vy1, real  *vz1
#endif//BLOCK_TIME_STEP
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
   , acceleration  *acc_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
   , real  *neighbor
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
   , position  *encBall, position  *encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
   , position  *box_min_hst, position  *box_max_hst, position  *icom_hst
#endif//USE_RECTANGULAR_BOX_FOR_LET
#ifdef  DPADD_FOR_ACC
   , DPacc  *tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
   , acceleration  *res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
   );

  void  freeTreeBuffer_dev
  (int  *failure, uint  *buffer, uint  *freeLst
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
   , uint  *freeNum, int  *active
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
   , unsigned long long int  *cycles_hst, unsigned long long int  *cycles_dev
#ifndef SERIALIZED_EXECUTION
   , unsigned long long int  *cycles_dist_hst, unsigned long long int  *cycles_dist_dev
#endif//SERIALIZED_EXECUTION
#endif//defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
   , unsigned long long int  *cycles_let_hst, unsigned long long int  *cycles_let_dev
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
   );
  muse allocTreeBuffer_dev
  (int **failure, uint **buffer, uint **freeLst,
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
   uint **freeNum, int **active,
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
   unsigned long long int **cycles_hst, unsigned long long int **cycles_dev,
#ifndef SERIALIZED_EXECUTION
   unsigned long long int **cycles_dist_hst, unsigned long long int **cycles_dist_dev,
#endif//SERIALIZED_EXECUTION
#endif//defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
   unsigned long long int **cycles_let_hst, unsigned long long int **cycles_let_dev,
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
   soaTreeWalkBuf *buf, const int num_max, const deviceProp gpu);

  void calcGravity_dev
  (const int grpNum
#ifdef  BLOCK_TIME_STEP
   , double *reduce, const int totNum
#endif//BLOCK_TIME_STEP
   , laneinfo * RESTRICT laneInfo, const iparticle pi, const soaTreeNode tree, const soaTreeWalkBuf buf
   , kernelStream *sinfo, deviceProp devProp, double *time
#   if  !defined(BLOCK_TIME_STEP) || defined(COMPARE_WITH_DIRECT_SOLVER) || defined(COUNT_INTERACTIONS) || defined(PRINT_PSEUDO_PARTICLE_INFO)
   , const int Ni
#endif//!defined(BLOCK_TIME_STEP) || defined(COMPARE_WITH_DIRECT_SOLVER) || defined(COUNT_INTERACTIONS) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
   , const potential_field sphe
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
   , const disk_potential disk
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
   , char *file
#endif//PRINT_PSEUDO_PARTICLE_INFO
#   if  defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
   , unsigned long long int *cycles_hst, unsigned long long int *cycles_dev
#ifndef SERIALIZED_EXECUTION
   , unsigned long long int *cycles_dist_hst, unsigned long long int *cycles_dist_dev
#endif//SERIALIZED_EXECUTION
#endif//defined(USE_CLOCK_CYCLES_FOR_BRENT_METHOD) || !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#ifndef SERIALIZED_EXECUTION
   , measuredTime *measured, const int pjNum
#ifdef  MPI_VIA_HOST
   , const soaTreeNode tree_hst
#endif//MPI_VIA_HOST
   , const int Nlet, domainInfo *let, const int Nstream_let, cudaStream_t stream_let[], MPIcfg_tree mpi
#ifdef  MONITOR_LETGEN_TIME
   , unsigned long long int *cycles_let_hst, unsigned long long int *cycles_let_dev
#endif//MONITOR_LETGEN_TIME
#ifdef  SWITCH_WITH_J_PARALLELIZATION
   , const bool transferMode, const int Ni_local, const int Ni_total
   , int * RESTRICT Ni_list, int * RESTRICT head_list, int * RESTRICT grpNum_list, int * RESTRICT displs
   , const int maxNgrp_ext, laneinfo * RESTRICT laneInfo_ext, laneinfo * RESTRICT laneInfo_ext_hst
   , laneinfo * RESTRICT laneInfo_hst, const iparticle pi_ext
#ifdef  MPI_VIA_HOST
   , const iparticle pi_ext_hst_loc, const iparticle pi_ext_hst_ful
#endif//MPI_VIA_HOST
#endif//SWITCH_WITH_J_PARALLELIZATION
#endif//SERIALIZED_EXECUTION
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#ifdef  COUNT_INTERACTIONS
   , iparticle_treeinfo treeInfo
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
   , gpu_clock *clockInfo, int *recordStep
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  COMPARE_WITH_DIRECT_SOLVER
   , const bool approxGravity
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
   , const real eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#endif//COMPARE_WITH_DIRECT_SOLVER
   );

  void setGlobalConstants_walk_dev_cu(const real newton_hst, const real eps2_hst
#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
				      , const real theta2_hst
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
				      );

#ifdef  SET_SINK_PARTICLES
  void set_loss_cone_dev(const int Nbh, const sinkparticle bh, const iparticle pi, const real lightspeed);
#endif//SET_SINK_PARTICLES

#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//WALK_DEV_H
