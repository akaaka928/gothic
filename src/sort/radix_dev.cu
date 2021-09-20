/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:44:17
 *                                                                       *
 *    Utility tool for stable sort based on radix sort                   *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef RADIX_DEV_CU
#define RADIX_DEV_CU
//-------------------------------------------------------------------------
#define USE_LSD_GRID_SORT_FUNC
//-------------------------------------------------------------------------
#define USE_SPLITTED_GRID_SORT_FUNC
//-------------------------------------------------------------------------
/* #define USE_SLICED_GRID_SORT_FUNC */
//-------------------------------------------------------------------------
/* #define USE_MODIFIED_GRID_SORT_FUNC */
//-------------------------------------------------------------------------
#ifndef SORT_ELEMENTS_PER_THREAD
#define SORT_ELEMENTS_PER_THREAD (4)
#endif//SORT_ELEMENTS_PER_THREAD
//-------------------------------------------------------------------------
#define PERFORMANCE_COMPARISON_MODE
//-------------------------------------------------------------------------
#define PERFORMANCE_COMPARISON_WITH_THRUST
//-------------------------------------------------------------------------
#define IMPLEMENTATION_CHECK_MODE
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <helper_cuda.h>
//-------------------------------------------------------------------------
#ifdef  IMPLEMENTATION_CHECK_MODE
#include <math.h>/* <-- for isfinite in test code */
#define TEST_SORT_ONLY_MODE
#ifdef  TEST_SORT_ONLY_MODE
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#endif//TEST_SORT_ONLY_MODE
#endif//IMPLEMENTATION_CHECK_MODE
//-------------------------------------------------------------------------
#ifdef  PERFORMANCE_COMPARISON_MODE
#define MEASURE_ITERATION_NUM (16)
#endif//PERFORMANCE_COMPARISON_MODE
//-------------------------------------------------------------------------
#include "macro.h"
#include "cudalib.h"
//-------------------------------------------------------------------------
#include "../misc/gsync_dev.cu"
#include "../sort/radix_dev.h"
//-------------------------------------------------------------------------
/* limitation from capacity of shared memory */
/* 48 KiB = 48 * 1024 B = 49152 B */
/* 16 KiB = 16 * 1024 B = 16384 B */
#ifdef  SMEM_PREF_FOR_GLOBAL_SORT
#       define SHARED_MEM_CAPACITY (49152)
#else///SMEM_PREF_FOR_GLOBAL_SORT
#       define SHARED_MEM_CAPACITY (16384)
#endif//SMEM_PREF_FOR_GLOBAL_SORT
//-------------------------------------------------------------------------
/* usage of shared memory per block (in units of B) */
#define USAGE_KEY32 (4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD)
#define USAGE_KEY64 (8 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD)
#define USAGE_IDX   (4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD)
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
#   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#define USAGE_SMEM (4 * NTHREADS_SORT)
#else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#define USAGE_SMEM (4 * 7 * RADIX_SORT_NUM_BUCKET)
#endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#else///USE_WARP_SHUFFLE_FUNC_SORT
#define USAGE_SMEM (4 * RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD)
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
#define USAGE_SBUF (4 * 4 * NTHREADS_SORT)
#else///RADIX_SORT_CHECK_BITS > 2
#define USAGE_SBUF (0)
#endif//RADIX_SORT_CHECK_BITS > 2
//-------------------------------------------------------------------------
#define NBLOCKS_PER_SM_SORT32        (SHARED_MEM_CAPACITY / (USAGE_KEY32             + USAGE_SMEM + USAGE_SBUF))
#define NBLOCKS_PER_SM_SORT64        (SHARED_MEM_CAPACITY / (USAGE_KEY64             + USAGE_SMEM + USAGE_SBUF))
#define NBLOCKS_PER_SM_SORT32_BY_KEY (SHARED_MEM_CAPACITY / (USAGE_KEY32 + USAGE_IDX + USAGE_SMEM + USAGE_SBUF))
#define NBLOCKS_PER_SM_SORT64_BY_KEY (SHARED_MEM_CAPACITY / (USAGE_KEY64 + USAGE_IDX + USAGE_SMEM + USAGE_SBUF))
//-------------------------------------------------------------------------
/* 64k registers per SM; a remedy to prevent register spill */
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  NBLOCKS_PER_SM_SORT32 > (65536 / (128 * NTHREADS_SORT))
#undef  NBLOCKS_PER_SM_SORT32
#define NBLOCKS_PER_SM_SORT32   (65536 / (128 * NTHREADS_SORT))
#endif//NBLOCKS_PER_SM_SORT32 > (65536 / (128 * NTHREADS_SORT))
#else///USE_MASK_INSTEAD_OF_IF
#   if  NBLOCKS_PER_SM_SORT32 > (65536 / ( 96 * NTHREADS_SORT))
#undef  NBLOCKS_PER_SM_SORT32
#define NBLOCKS_PER_SM_SORT32   (65536 / ( 96 * NTHREADS_SORT))
#endif//NBLOCKS_PER_SM_SORT32 > (65536 / ( 96 * NTHREADS_SORT))
#endif//USE_MASK_INSTEAD_OF_IF
//-------------------------------------------------------------------------
#   if  NBLOCKS_PER_SM_SORT32 < 1
#undef  NBLOCKS_PER_SM_SORT32
#define NBLOCKS_PER_SM_SORT32 (1)
#endif//NBLOCKS_PER_SM_SORT32
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define KEY_BITS (32)
#include "../sort/radix_inc.cu"
#define SORT_ONLY
#include "../sort/radix_inc.cu"
//-------------------------------------------------------------------------
#undef KEY_BITS
#undef SORT_ONLY
//-------------------------------------------------------------------------
#define KEY_BITS (64)
#include "../sort/radix_inc.cu"
#define SORT_ONLY
#include "../sort/radix_inc.cu"
//-------------------------------------------------------------------------
#undef KEY_BITS
#undef SORT_ONLY
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* radix sort for 32 bits unsigned integer (uint) */
//-------------------------------------------------------------------------
extern "C"
__global__ void radixSortTsub_uint(const int num, uint *key)
{
  //-----------------------------------------------------------------------
  __shared__ uint key_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  __shared__ uint    smem[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
  __shared__ uint4_array sbuf[NTHREADS_SORT];
#endif//RADIX_SORT_CHECK_BITS > 2
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------
  const int lane = tidx & (TSUB_SORT - 1);
  const int head = tidx - lane;
  const int hp = head * SORT_ELEMENTS_PER_THREAD;
  //-----------------------------------------------------------------------
#pragma unroll
  for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){
    const int idx = hp + lane + ii * TSUB_SORT;
    key_sm[idx] = (idx < num) ? (key[idx]) : (0xffffffff);
  }/* for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
  //-----------------------------------------------------------------------
  __radixSortTsub32(32 / RADIX_SORT_CHECK_BITS, lane, &key_sm[hp]
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
		    , tidx, &smem[hp]
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
		    , &sbuf[head]
#endif//RADIX_SORT_CHECK_BITS > 2
		    );
  //-----------------------------------------------------------------------
#pragma unroll
  for(int ii = lane; ii < num; ii += TSUB_SORT)
    key[hp + ii] = key_sm[hp + ii];
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* extern "C" */
/* __global__ void radixSortTsub_uint_by_key(const int num, READ_ONLY uint * RESTRICT key, uint * RESTRICT dat) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __shared__  int idx_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   __shared__ uint key_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   __shared__ uint smem[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* #   if  RADIX_SORT_CHECK_BITS > 2 */
/*   __shared__ uint4_array sbuf[NTHREADS_SORT]; */
/* #endif//RADIX_SORT_CHECK_BITS > 2 */
/*   //----------------------------------------------------------------------- */
/*   const int tidx = THREADIDX_X1D; */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (TSUB_SORT - 1); */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/*   //----------------------------------------------------------------------- */
/* #pragma unroll */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = hp + lane + ii * TSUB_SORT; */
/*     idx_sm[idx] = idx; */
/*     key_sm[idx] = (idx < num) ? (key[idx]) : (0xffffffff); */
/*   }/\* for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   __radixSortTsub32idx(32 / RADIX_SORT_CHECK_BITS, lane, &key_sm[hp] */
/* 		       , &idx_sm[hp] */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 		       , tidx, &smem[hp] */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* #   if  RADIX_SORT_CHECK_BITS > 2 */
/* 		       , &sbuf[head] */
/* #endif//RADIX_SORT_CHECK_BITS > 2 */
/* 		       ); */
/*   //----------------------------------------------------------------------- */
/* #pragma unroll */
/*   for(int ii = lane; ii < num; ii += TSUB_SORT) */
/*     key_sm[hp + ii] = dat[idx_sm[hp + ii]]; */
/*   //----------------------------------------------------------------------- */
/* #pragma unroll */
/*   for(int ii = lane; ii < num; ii += TSUB_SORT) */
/*     dat[hp + ii] = key_sm[hp + ii]; */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
extern "C"
__global__ void radixSortBlck_uint(const int num, uint *key)
{
  //-----------------------------------------------------------------------
  __shared__ uint key_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  __shared__ uint_int smem[       RADIX_SORT_NUM_OF_I32       *  NTHREADS_SORT       * SORT_ELEMENTS_PER_THREAD];
#else///USE_WARP_SHUFFLE_FUNC_SORT
  __shared__ uint_int smem[(1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * SORT_ELEMENTS_PER_THREAD];
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
  __shared__ uint4_array sbuf[NTHREADS_SORT];
#endif//RADIX_SORT_CHECK_BITS > 2
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int lane = tidx & (warpSize - 1);
#  if  SORT_ELEMENTS_PER_THREAD != 1
  const int head = tidx - lane;
  const int hp = head * SORT_ELEMENTS_PER_THREAD;
#endif//SORT_ELEMENTS_PER_THREAD != 1
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------
#pragma unroll
  for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){
    const int idx = tidx + ii * NTHREADS_SORT;
    key_sm[idx] = (idx < num) ? (key[idx]) : (0xffffffff);
  }/* for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
  //-----------------------------------------------------------------------
  __radixSortBlck32(32 / RADIX_SORT_CHECK_BITS, tidx,
#   if  SORT_ELEMENTS_PER_THREAD != 1
		    hp,
#endif//SORT_ELEMENTS_PER_THREAD != 1
		    key_sm, smem
#   if  RADIX_SORT_CHECK_BITS > 2
		    , sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
		    , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
		    );
  //-----------------------------------------------------------------------
#pragma unroll
  for(int ii = tidx; ii < num; ii += NTHREADS_SORT)
    key[ii] = key_sm[ii];
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* extern "C" */
/* __global__ void radixSortBlck_uint_by_key(const int num, READ_ONLY uint * RESTRICT key, uint * RESTRICT dat) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __shared__  int idx_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   __shared__ uint key_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   __shared__ uint_int smem[       RADIX_SORT_NUM_OF_I32       *  NTHREADS_SORT       * SORT_ELEMENTS_PER_THREAD]; */
/* #else///USE_WARP_SHUFFLE_FUNC_SORT */
/*   __shared__ uint_int smem[(1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* #   if  RADIX_SORT_CHECK_BITS > 2 */
/*   __shared__ uint4_array sbuf[NTHREADS_SORT]; */
/* #endif//RADIX_SORT_CHECK_BITS > 2 */
/*   //----------------------------------------------------------------------- */
/*   const int tidx = THREADIDX_X1D; */
/* #   if  SORT_ELEMENTS_PER_THREAD != 1 */
/*   const int lane = tidx & (warpSize - 1); */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #endif//SORT_ELEMENTS_PER_THREAD != 1 */
/*   //----------------------------------------------------------------------- */
/* #pragma unroll */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = tidx + ii * NTHREADS_SORT; */
/*     idx_sm[idx] = idx; */
/*     key_sm[idx] = (idx < num) ? (key[idx]) : (0xffffffff); */
/*   }/\* for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   __radixSortBlck32idx(32 / RADIX_SORT_CHECK_BITS, tidx, */
/* #   if  SORT_ELEMENTS_PER_THREAD != 1 */
/* 		       hp, */
/* #endif//SORT_ELEMENTS_PER_THREAD != 1 */
/* 		       key_sm, idx_sm, smem */
/* #   if  RADIX_SORT_CHECK_BITS > 2 */
/* 		       , sbuf */
/* #endif//RADIX_SORT_CHECK_BITS > 2 */
/* 		       ); */
/*   //----------------------------------------------------------------------- */
/* #pragma unroll */
/*   for(int ii = tidx; ii < num; ii += NTHREADS_SORT) */
/*     key_sm[ii] = dat[idx_sm[ii]]; */
/*   __syncthreads(); */
/*   //----------------------------------------------------------------------- */
/* #pragma unroll */
/*   for(int ii = tidx; ii < num; ii += NTHREADS_SORT) */
/*     dat[ii] = key_sm[ii]; */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
extern "C"
__global__ void __launch_bounds__(NTHREADS_SORT, NBLOCKS_PER_SM_SORT32) radixSortGrid_uint
     (const int num, uint * RESTRICT key,
      int * RESTRICT numIterFul, int2 * RESTRICT infoFul, int * RESTRICT numIterLoc, int2 * RESTRICT infoLoc, const int remLocMax,
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
      int * RESTRICT checkedCounts,
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
      int * RESTRICT fail_dev,
      uint * RESTRICT key_tmp_gm,
      int * RESTRICT scan_gm, int * RESTRICT bucket_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul,
      int * RESTRICT gsync0, int * RESTRICT gsync1
#ifdef  USE_MODIFIED_GRID_SORT_FUNC
      , int * RESTRICT gsync0Loc, int * RESTRICT gsync1Loc, int * RESTRICT remFulAdd, int * RESTRICT remLocAdd
#endif//USE_MODIFIED_GRID_SORT_FUNC
      )
{
  //-----------------------------------------------------------------------
  /* 4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes */
  __shared__ uint key_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  /* 4 * RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes */
  __shared__ uint_int smem[RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#else///USE_WARP_SHUFFLE_FUNC_SORT
  /* max of */
  /* (1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * SORT_ELEMENTS_PER_THREAD, */
  /* 7 * RADIX_SORT_NUM_BUCKET, */
  /* and NTHREADS_SORT */
  /* RADIX_SORT_CHECK_BITS == 4 --> 8 * SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT / 32 <= NTHREADS_SORT */
  /* RADIX_SORT_CHECK_BITS == 2 --> 2 * SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT / 32 <= NTHREADS_SORT / 4 */
  /* therefore, maximum is NTHREADS_SORT (because, typically this value is >= 128 = 8 * 16 and RADIX_SORT_NUM_BUCKET <= 16) */
#   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[            NTHREADS_SORT];
#else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[7 * RADIX_SORT_NUM_BUCKET];
#endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
  __shared__ uint4_array sbuf[NTHREADS_SORT];
#endif//RADIX_SORT_CHECK_BITS > 2
  //-----------------------------------------------------------------------
  const int bnum =   GRIDDIM_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  const int lane = tidx & (warpSize - 1);
#   if  SORT_ELEMENTS_PER_THREAD != 1
  const int head = tidx - lane;
  const int hp = head * SORT_ELEMENTS_PER_THREAD;
#endif//SORT_ELEMENTS_PER_THREAD != 1
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* data preparation */
  //-----------------------------------------------------------------------
  if( gidx == 0 ){
    //---------------------------------------------------------------------
    int2 info = {num, 0};
    infoFul[0] = info;
    numIterFul[0] = 32 / RADIX_SORT_CHECK_BITS;
    /* numIterFul[0] = 1; */
    //---------------------------------------------------------------------
  }/* if( gidx == 0 ){ */
  //-----------------------------------------------------------------------
  /* data conversion is necessary for int or float */
  globalSync(tidx, bidx, bnum, gsync0, gsync1);
  //-----------------------------------------------------------------------
#ifdef  USE_MODIFIED_GRID_SORT_FUNC
  __radixSortGrid32_mod
    (numIterFul, infoFul, 1, numIterLoc, infoLoc, remLocMax,
     remFulAdd, remLocAdd,
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
     checkedCounts,
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
     fail_dev, tidx, lane,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
     hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
     key, key_tmp_gm, key_sm,
#   if  RADIX_SORT_CHECK_BITS > 2
     sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
     smem, scan_gm, bucket_gm, tail_gm, scanFul,
     gidx, bidx, bnum, gsync0, gsync1, gsync0Loc, gsync1Loc
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
#else///USE_MODIFIED_GRID_SORT_FUNC
  __radixSortGrid32
    (numIterFul, infoFul, 1, numIterLoc, infoLoc, remLocMax,
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
     checkedCounts,
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
     fail_dev, tidx, lane,
#   if  SORT_ELEMENTS_PER_THREAD != 1
     hp,
#endif//SORT_ELEMENTS_PER_THREAD != 1
     key, key_tmp_gm, key_sm,
#   if  RADIX_SORT_CHECK_BITS > 2
     sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
     smem, scan_gm, bucket_gm, tail_gm, scanFul,
     gidx, bidx, bnum, gsync0, gsync1
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
#endif//USE_MODIFIED_GRID_SORT_FUNC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifdef  USE_LSD_GRID_SORT_FUNC
//-------------------------------------------------------------------------
extern "C"
__global__ void __launch_bounds__(NTHREADS_SORT, NBLOCKS_PER_SM_SORT32) radixSortGridLSD_uint
     (const int num, uint * RESTRICT key, uint * RESTRICT key_tmp_gm,
      int * RESTRICT scan_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  //-----------------------------------------------------------------------
  /* 4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes */
  __shared__ uint key_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  /* 4 * RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes */
  __shared__ uint_int smem[RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#else///USE_WARP_SHUFFLE_FUNC_SORT
#   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[            NTHREADS_SORT];
#else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[7 * RADIX_SORT_NUM_BUCKET];
#endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
  __shared__ uint4_array sbuf[NTHREADS_SORT];
#endif//RADIX_SORT_CHECK_BITS > 2
  //-----------------------------------------------------------------------
  const int bnum =   GRIDDIM_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  const int lane = tidx & (warpSize - 1);
#   if  SORT_ELEMENTS_PER_THREAD != 1
  const int head = tidx - lane;
  const int hp = head * SORT_ELEMENTS_PER_THREAD;
#endif//SORT_ELEMENTS_PER_THREAD != 1
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------

  /* //----------------------------------------------------------------------- */
  /* /\* data preparation *\/ */
  /* //----------------------------------------------------------------------- */
  /* /\* data conversion is necessary for int or float *\/ */
  /* globalSync(tidx, bidx, bnum, gsync0, gsync1); */
  //-----------------------------------------------------------------------
  __radixSortLsdGrid32
    (32 / RADIX_SORT_CHECK_BITS, num,
     tidx, lane,
#   if  SORT_ELEMENTS_PER_THREAD != 1
     hp,
#endif//SORT_ELEMENTS_PER_THREAD != 1
     key, key_tmp_gm, key_sm,
#   if  RADIX_SORT_CHECK_BITS > 2
     sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
     smem, scan_gm, tail_gm, scanFul,
     gidx, bidx, bnum, gsync0, gsync1
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_LSD_GRID_SORT_FUNC
//-------------------------------------------------------------------------
#   if  defined(USE_LSD_GRID_SORT_FUNC) && defined(USE_SPLITTED_GRID_SORT_FUNC)
//-------------------------------------------------------------------------
extern "C"
__global__ void radixSortGridLSD_uint_1st
(const int iter, const int num, uint * RESTRICT key, uint * RESTRICT key_tmp_gm, int * RESTRICT scan_gm, int * RESTRICT scanFul)
{
  //-----------------------------------------------------------------------
  /* 4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes */
  __shared__ uint key_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  /* 4 * RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes */
  __shared__ uint_int smem[RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#else///USE_WARP_SHUFFLE_FUNC_SORT
#   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[            NTHREADS_SORT];
#else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[7 * RADIX_SORT_NUM_BUCKET];
#endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
  __shared__ uint4_array sbuf[NTHREADS_SORT];
#endif//RADIX_SORT_CHECK_BITS > 2
  //-----------------------------------------------------------------------
  const int bnum =   GRIDDIM_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------
  const int lane = tidx & (warpSize - 1);
#   if  SORT_ELEMENTS_PER_THREAD != 1
  const int head = tidx - lane;
  const int hp = head * SORT_ELEMENTS_PER_THREAD;
#endif//SORT_ELEMENTS_PER_THREAD != 1
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __radixSortLsdGrid32_1st
    (iter, num,
     tidx, lane,
#   if  SORT_ELEMENTS_PER_THREAD != 1
     hp,
#endif//SORT_ELEMENTS_PER_THREAD != 1
     key, key_tmp_gm, key_sm,
#   if  RADIX_SORT_CHECK_BITS > 2
     sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
     smem, scan_gm, scanFul, bidx, bnum
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
__global__ void radixSortGridLSD_uint_2nd(const int bnum, int * RESTRICT tail_gm, int * RESTRICT scanFul)
{
  //-----------------------------------------------------------------------
  __shared__ uint_int smem[NTHREADS_SORT_ACCUMULATION];
  //-----------------------------------------------------------------------
  const int bidx =  BLOCKIDX_X1D;
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int lane = tidx & (warpSize - 1);
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __radixSortLsdGrid32_2nd
    (tidx, lane, smem, tail_gm, scanFul, bidx, bnum
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
__global__ void radixSortGridLSD_uint_3rd
(const int iter, const int num, uint * RESTRICT key, uint * RESTRICT key_tmp_gm,
 READ_ONLY int * RESTRICT scan_gm, READ_ONLY int * RESTRICT tail_gm, READ_ONLY int * RESTRICT scanFul)
{
  //-----------------------------------------------------------------------
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  /* 4 * RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes */
  __shared__ uint_int smem[RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#else///USE_WARP_SHUFFLE_FUNC_SORT
#   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[            NTHREADS_SORT];
#else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[7 * RADIX_SORT_NUM_BUCKET];
#endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#endif//USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
  const int bnum =   GRIDDIM_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------
  const int lane = tidx & (warpSize - 1);
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __radixSortLsdGrid32_3rd
    (iter, num, tidx, lane,
     key, key_tmp_gm,
     smem, scan_gm, tail_gm, scanFul, bidx, bnum
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//defined(USE_LSD_GRID_SORT_FUNC) && defined(USE_SPLITTED_GRID_SORT_FUNC)
//-------------------------------------------------------------------------
/* extern "C" */
/* __global__ void __launch_bounds__(NTHREADS_SORT, NBLOCKS_PER_SM_SORT32) radixSortGrid_uint_by_key */
/*      (const int num, uint * RESTRICT key, int * RESTRICT idx, uint * RESTRICT dat, */
/*       int * RESTRICT numIterFul, int2 * RESTRICT infoFul, int * RESTRICT numIterLoc, int2 * RESTRICT infoLoc, const int remLocMax, */
/* #ifdef  ATOMIC_BASED_JOB_ASSIGNMENT */
/*       int * RESTRICT checkedCounts, */
/* #endif//ATOMIC_BASED_JOB_ASSIGNMENT */
/*       int * RESTRICT fail_dev, uint * RESTRICT key_tmp_gm, int * RESTRICT idx_tmp_gm, */
/*       int * RESTRICT scan_gm, int * RESTRICT bucket_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul, */
/*       int * RESTRICT gsync0, int * RESTRICT gsync1 */
/*       ) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   /\* 4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes *\/ */
/*   /\* 4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes *\/ */
/*   __shared__  int idx_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   __shared__ uint key_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   /\* 4 * 4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes *\/ */
/*   __shared__ uint_int smem[RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #else///USE_WARP_SHUFFLE_FUNC_SORT */
/* #   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET) */
/*   __shared__ uint_int smem[            NTHREADS_SORT]; */
/* #else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET) */
/*   __shared__ uint_int smem[7 * RADIX_SORT_NUM_BUCKET]; */
/* #endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET) */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* #   if  RADIX_SORT_CHECK_BITS > 2 */
/*   __shared__ uint4_array sbuf[NTHREADS_SORT]; */
/* #endif//RADIX_SORT_CHECK_BITS > 2 */
/*   //----------------------------------------------------------------------- */
/*   const int bnum =   GRIDDIM_X1D; */
/*   const int bidx =  BLOCKIDX_X1D; */
/*   const int tidx = THREADIDX_X1D; */
/*   const int gidx = GLOBALIDX_X1D; */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (warpSize - 1); */
/* #   if  SORT_ELEMENTS_PER_THREAD != 1 */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #endif//SORT_ELEMENTS_PER_THREAD != 1 */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* data preparation *\/ */
/*   //----------------------------------------------------------------------- */
/*   if( gidx == 0 ){ */
/*     //--------------------------------------------------------------------- */
/*     int2 info = {num, 0}; */
/*     infoFul[0] = info; */
/*     numIterFul[0] = 32 / RADIX_SORT_CHECK_BITS; */
/*     //--------------------------------------------------------------------- */
/*   }/\* if( gidx == 0 ){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   /\* data conversion is necessary for int or float *\/ */
/*   globalSync(tidx, bidx, bnum, gsync0, gsync1); */
/*   //----------------------------------------------------------------------- */
/*   __radixSortGrid32idx */
/*     (numIterFul, infoFul, 1, numIterLoc, infoLoc, remLocMax, */
/* #ifdef  ATOMIC_BASED_JOB_ASSIGNMENT */
/*      checkedCounts, */
/* #endif//ATOMIC_BASED_JOB_ASSIGNMENT */
/*      fail_dev, tidx, lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD != 1 */
/*      hp, */
/* #endif//SORT_ELEMENTS_PER_THREAD != 1 */
/*      key, key_tmp_gm, key_sm, */
/*      idx, idx_tmp_gm, idx_sm, */
/* #   if  RADIX_SORT_CHECK_BITS > 2 */
/*      sbuf, */
/* #endif//RADIX_SORT_CHECK_BITS > 2 */
/*      smem, scan_gm, bucket_gm, tail_gm, scanFul, */
/*      gidx, bidx, bnum, gsync0, gsync1); */
/*   //----------------------------------------------------------------------- */
/*   globalSync(tidx, bidx, bnum, gsync0, gsync1); */
/*   for(int ii = gidx; ii < num; ii += bnum * NTHREADS_SORT) */
/*     key_tmp_gm[ii] = dat[idx[ii]]; */
/*   //----------------------------------------------------------------------- */
/*   globalSync(tidx, bidx, bnum, gsync0, gsync1); */
/*   for(int ii = gidx; ii < num; ii += bnum * NTHREADS_SORT) */
/*     dat[ii] = key_tmp_gm[ii]; */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
//-------------------------------------------------------------------------
extern "C"
__global__ void __launch_bounds__(NTHREADS_SORT, NBLOCKS_PER_SM_SORT32) radixSortGrid_uint_msd
     (const int num, uint * RESTRICT key,
      int * RESTRICT numIterFul, int2 * RESTRICT infoFul, int * RESTRICT numIterLoc, int2 * RESTRICT infoLoc, const int remLocMax,
      int * RESTRICT remLoc_dev, int * RESTRICT fail_dev, uint * RESTRICT key_tmp_gm,
      int * RESTRICT scan_gm, int * RESTRICT bucket_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul,
      int * RESTRICT gsync0, int * RESTRICT gsync1
      )
{
  //-----------------------------------------------------------------------
  __shared__ uint   key_sm[                        NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  __shared__ uint_int smem[RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#else///USE_WARP_SHUFFLE_FUNC_SORT
#   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[            NTHREADS_SORT];
#else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[7 * RADIX_SORT_NUM_BUCKET];
#endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
  __shared__ uint4_array sbuf[NTHREADS_SORT];
#endif//RADIX_SORT_CHECK_BITS > 2
  //-----------------------------------------------------------------------
  const int bnum =   GRIDDIM_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  const int lane = tidx & (warpSize - 1);
#   if  SORT_ELEMENTS_PER_THREAD != 1
  const int head = tidx - lane;
  const int hp = head * SORT_ELEMENTS_PER_THREAD;
#endif//SORT_ELEMENTS_PER_THREAD != 1
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* data preparation */
  //-----------------------------------------------------------------------
  if( gidx == 0 ){
    //---------------------------------------------------------------------
    int2 info = {num, 0};
    infoFul[0] = info;
    numIterFul[0] = 32 / RADIX_SORT_CHECK_BITS;
    //---------------------------------------------------------------------
  }/* if( gidx == 0 ){ */
  //-----------------------------------------------------------------------
  /* data conversion is necessary for int or float */
  globalSync(tidx, bidx, bnum, gsync0, gsync1);
  //-----------------------------------------------------------------------
  __radixSortGrid32_msd
    (numIterFul, infoFul, 1, numIterLoc, infoLoc, remLocMax,
     remLoc_dev, fail_dev, tidx, lane,
#   if  SORT_ELEMENTS_PER_THREAD != 1
     hp,
#endif//SORT_ELEMENTS_PER_THREAD != 1
     key, key_tmp_gm, key_sm,
#   if  RADIX_SORT_CHECK_BITS > 2
     sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
     smem, scan_gm, bucket_gm, tail_gm, scanFul
     , gidx, bidx, bnum, gsync0, gsync1
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_SPLITTED_GRID_SORT_FUNC
//-------------------------------------------------------------------------
#ifdef  USE_SLICED_GRID_SORT_FUNC
//-------------------------------------------------------------------------
extern "C"
__global__ void radixSortGrid_mod_scatter_uint
(int * RESTRICT numIterFul, int2 * RESTRICT infoFul, int remFul, uint * RESTRICT key, uint * RESTRICT key_tmp_gm,
 int * RESTRICT snum, int2 * RESTRICT disp, int * RESTRICT scan_gm, int * RESTRICT scanFul)
{
  //-----------------------------------------------------------------------
  __shared__ uint   key_sm[                        NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  __shared__ uint_int smem[RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#else///USE_WARP_SHUFFLE_FUNC_SORT
#   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[            NTHREADS_SORT];
#else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[7 * RADIX_SORT_NUM_BUCKET];
#endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
  __shared__ uint4_array sbuf[NTHREADS_SORT];
#endif//RADIX_SORT_CHECK_BITS > 2
  //-----------------------------------------------------------------------
  const int bnum =   GRIDDIM_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  const int lane = tidx & (warpSize - 1);
#   if  SORT_ELEMENTS_PER_THREAD != 1
  const int head = tidx - lane;
  const int hp = head * SORT_ELEMENTS_PER_THREAD;
#endif//SORT_ELEMENTS_PER_THREAD != 1
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __radixSortGrid32_mod_scatter
    (numIterFul, infoFul, remFul, snum, disp,
     gidx, bidx, bnum, tidx, lane,
#   if  SORT_ELEMENTS_PER_THREAD != 1
     hp,
#endif//SORT_ELEMENTS_PER_THREAD != 1
     key, key_tmp_gm, key_sm,
#   if  RADIX_SORT_CHECK_BITS > 2
     sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
     smem, scan_gm, scanFul
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
  //-----------------------------------------------------------------------
#if 0
  if( gidx == 0 )
    printf("snum = %d\n", *snum);
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* __global__ void radixSortGrid_mod_scatter_uint */
/* (int * RESTRICT numIterFul, int2 * RESTRICT infoFul, int remFul, uint * RESTRICT key, uint * RESTRICT key_tmp_gm, */
/*  int * RESTRICT snum, int2 * RESTRICT disp, int * RESTRICT scan_gm, int * RESTRICT scanFul) */
extern "C"
__global__ void radixSortGrid_mod_gather_uint
(int * RESTRICT numIterSrc, int2 * RESTRICT infoSrc, int            remSrc,
 int * RESTRICT numIterDst, int2 * RESTRICT infoDst, int * RESTRICT remDstAdd,
 int * RESTRICT numIterLoc, int2 * RESTRICT infoLoc, int * RESTRICT remLocAdd,
 int * RESTRICT snum, const int scanned,
 uint * RESTRICT key, uint * RESTRICT key_tmp_gm,
 int * RESTRICT scan_gm, int * RESTRICT bucket_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul)
{
  //-----------------------------------------------------------------------
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  __shared__ uint_int smem[RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#else///USE_WARP_SHUFFLE_FUNC_SORT
#   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[            NTHREADS_SORT];
#else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[7 * RADIX_SORT_NUM_BUCKET];
#endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
  __shared__ uint4_array sbuf[NTHREADS_SORT];
#endif//RADIX_SORT_CHECK_BITS > 2
  //-----------------------------------------------------------------------
  const int bnum =   GRIDDIM_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  const int lane = tidx & (warpSize - 1);
#   if  SORT_ELEMENTS_PER_THREAD != 1
  const int head = tidx - lane;
  const int hp = head * SORT_ELEMENTS_PER_THREAD;
#endif//SORT_ELEMENTS_PER_THREAD != 1
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __radixSortGrid32_mod_gather
    (numIterSrc, infoSrc, remSrc, numIterDst, infoDst, remDstAdd, numIterLoc, infoLoc, remLocAdd,
     snum, scanned,
     gidx, bidx, bnum, tidx, lane,
#   if  SORT_ELEMENTS_PER_THREAD != 1
     hp,
#endif//SORT_ELEMENTS_PER_THREAD != 1
     key, key_tmp_gm,
#   if  RADIX_SORT_CHECK_BITS > 2
     sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
     smem, scan_gm, bucket_gm, tail_gm, scanFul
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
  //-----------------------------------------------------------------------
#if 0
  if( gidx == 0 )
    printf("snum = %d\n", *snum);
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_SLICED_GRID_SORT_FUNC
//-------------------------------------------------------------------------
/* extern "C" */
/* __global__ void __launch_bounds__(NTHREADS_SORT, NBLOCKS_PER_SM_SORT32) radixSortGrid_uint_by_key_msd */
/*      (const int num, uint * RESTRICT key, int * RESTRICT idx, */
/*       int * RESTRICT numIterFul, int2 * RESTRICT infoFul, int * RESTRICT numIterLoc, int2 * RESTRICT infoLoc, const int remLocMax, */
/*       int * RESTRICT remLoc_dev, int * RESTRICT fail_dev, uint * RESTRICT key_tmp_gm, int * RESTRICT idx_tmp_gm, */
/*       int * RESTRICT scan_gm, int * RESTRICT bucket_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul, */
/*       int * RESTRICT gsync0, int * RESTRICT gsync1 */
/*       ) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   /\* 4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes *\/ */
/*   /\* 4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes *\/ */
/*   __shared__  int idx_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   __shared__ uint key_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   /\* 4 * 4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes *\/ */
/*   __shared__ uint_int smem[RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #else///USE_WARP_SHUFFLE_FUNC_SORT */
/* #   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET) */
/*   __shared__ uint_int smem[            NTHREADS_SORT]; */
/* #else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET) */
/*   __shared__ uint_int smem[7 * RADIX_SORT_NUM_BUCKET]; */
/* #endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET) */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* #   if  RADIX_SORT_CHECK_BITS > 2 */
/*   __shared__ uint4_array sbuf[NTHREADS_SORT]; */
/* #endif//RADIX_SORT_CHECK_BITS > 2 */
/*   //----------------------------------------------------------------------- */
/*   const int bnum =   GRIDDIM_X1D; */
/*   const int bidx =  BLOCKIDX_X1D; */
/*   const int tidx = THREADIDX_X1D; */
/*   const int gidx = GLOBALIDX_X1D; */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (warpSize - 1); */
/* #   if  SORT_ELEMENTS_PER_THREAD != 1 */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #endif//SORT_ELEMENTS_PER_THREAD != 1 */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* data preparation *\/ */
/*   //----------------------------------------------------------------------- */
/*   if( gidx == 0 ){ */
/*     //--------------------------------------------------------------------- */
/*     int2 info = {num, 0}; */
/*     infoFul[0] = info; */
/*     numIterFul[0] = 32 / RADIX_SORT_CHECK_BITS; */
/*     //--------------------------------------------------------------------- */
/*   }/\* if( gidx == 0 ){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   for(int ii = gidx; ii < num; ii += bnum * NTHREADS_SORT) */
/*     idx[ii] = ii; */
/*   //----------------------------------------------------------------------- */
/*   /\* data conversion is necessary for int or float *\/ */
/*   globalSync(tidx, bidx, bnum, gsync0, gsync1); */
/*   //----------------------------------------------------------------------- */
/*   __radixSortGrid32idx_msd */
/*     (numIterFul, infoFul, 1, numIterLoc, infoLoc, remLocMax, */
/*      remLoc_dev, fail_dev, tidx, lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD != 1 */
/*      hp, */
/* #endif//SORT_ELEMENTS_PER_THREAD != 1 */
/*      key, key_tmp_gm, key_sm, */
/*      idx, idx_tmp_gm, idx_sm, */
/* #   if  RADIX_SORT_CHECK_BITS > 2 */
/*      sbuf, */
/* #endif//RADIX_SORT_CHECK_BITS > 2 */
/*      smem, scan_gm, bucket_gm, tail_gm, scanFul, */
/*      gidx, bidx, bnum, gsync0, gsync1); */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
#   if  defined(USE_SLICED_GRID_SORT_FUNC) || defined(USE_SPLITTED_GRID_SORT_FUNC)
//-------------------------------------------------------------------------
extern "C"
__global__ void radixSortGrid_uint_lsd
(uint * RESTRICT key, READ_ONLY int * RESTRICT numIterLoc, READ_ONLY int2 * RESTRICT infoLoc)
{
  //-----------------------------------------------------------------------
  __shared__ uint   key_sm[                        NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  __shared__ uint_int smem[RADIX_SORT_NUM_OF_I32 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD];
#else///USE_WARP_SHUFFLE_FUNC_SORT
#   if  NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[            NTHREADS_SORT];
#else///NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
  __shared__ uint_int smem[7 * RADIX_SORT_NUM_BUCKET];
#endif//NTHREADS_SORT > (7 * RADIX_SORT_NUM_BUCKET)
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
  __shared__ uint4_array sbuf[NTHREADS_SORT];
#endif//RADIX_SORT_CHECK_BITS > 2
  //-----------------------------------------------------------------------
  const int bidx =  BLOCKIDX_X1D;
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------
  const int lane = tidx & (warpSize - 1);
#   if  SORT_ELEMENTS_PER_THREAD != 1
  const int hp = (tidx - lane) * SORT_ELEMENTS_PER_THREAD;
#endif//SORT_ELEMENTS_PER_THREAD != 1
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    int2 info =    infoLoc[bidx];
    smem[0].i = numIterLoc[bidx];
    smem[1].i = info.x;
    smem[2].i = info.y;
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  __syncthreads();
  const int numIter = smem[0].i;
  const int rem     = smem[1].i;
  const int head    = smem[2].i;
  //-----------------------------------------------------------------------
  __radixSortGrid32_lsd
    (numIter, rem, head, tidx, lane,
#   if  SORT_ELEMENTS_PER_THREAD != 1
     hp,
#endif//SORT_ELEMENTS_PER_THREAD != 1
     key, key_sm,
#   if  RADIX_SORT_CHECK_BITS > 2
     sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
     smem
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//defined(USE_SLICED_GRID_SORT_FUNC) || defined(USE_SPLITTED_GRID_SORT_FUNC)
//-------------------------------------------------------------------------






//-------------------------------------------------------------------------
/* #define RADIX_SORT_WITHIN_WARP */
/* __device__ __forceinline__ void radixSort32bitsWarp */
/* (const int tidx, const int bidx, const int bnum, volatile int * gsync0, volatile int * gsync1) */
/* { */
/*   //----------------------------------------------------------------------- */
/* #include "../sort/radix_src.cu" */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
/* #define RADIX_SORT_WITHIN_BLOCK */
/* __device__ __forceinline__ void radixSort32bitsBlock */
/* (const int tidx, const int bidx, const int bnum, volatile int * gsync0, volatile int * gsync1) */
/* #include "../sort/radix_src.cu" */
//-------------------------------------------------------------------------
/* #define RADIX_SORT_WITHIN_GRID */
/* __device__ __forceinline__ void radixSort32bitsGrid */
/* (const int tidx, const int bidx, const int bnum, volatile int * gsync0, volatile int * gsync1) */
/* #include "../sort/radix_src.cu" */
//-------------------------------------------------------------------------
/* #define RADIX_SORT_LARGE_ARRAY */
/* __device__ __forceinline__ void radixSort32bitsArray */
/* (const int tidx, const int bidx, const int bnum, volatile int * gsync0, volatile int * gsync1) */
/* #include "../sort/radix_src.cu" */
//-------------------------------------------------------------------------
/* #undef RADIX_SORT_WITHIN_WARP */
/* #undef RADIX_SORT_WITHIN_BLOCK */
/* #undef RADIX_SORT_WITHIN_GRID */
/* #undef RADIX_SORT_LARGE_ARRAY */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define RADIX_SORT_WITHIN_WARP */
/* __device__ __forceinline__ void radixSort64bitsWarp */
/* (const int tidx, const int bidx, const int bnum, volatile int * gsync0, volatile int * gsync1) */
/* #include "../sort/radix_src.cu" */
//-------------------------------------------------------------------------
/* #define RADIX_SORT_WITHIN_BLOCK */
/* __device__ __forceinline__ void radixSort64bitsBlock */
/* (const int tidx, const int bidx, const int bnum, volatile int * gsync0, volatile int * gsync1) */
/* #include "../sort/radix_src.cu" */
//-------------------------------------------------------------------------
/* #define RADIX_SORT_WITHIN_GRID */
/* __device__ __forceinline__ void radixSort64bitsGrid */
/* (const int tidx, const int bidx, const int bnum, volatile int * gsync0, volatile int * gsync1) */
/* #include "../sort/radix_src.cu" */
//-------------------------------------------------------------------------
/* #define RADIX_SORT_LARGE_ARRAY */
/* __device__ __forceinline__ void radixSort64bitsArray */
/* (const int tidx, const int bidx, const int bnum, volatile int * gsync0, volatile int * gsync1) */
/* #include "../sort/radix_src.cu" */
//-------------------------------------------------------------------------
/* #undef RADIX_SORT_WITHIN_WARP */
/* #undef RADIX_SORT_WITHIN_BLOCK */
/* #undef RADIX_SORT_WITHIN_GRID */
/* #undef RADIX_SORT_LARGE_ARRAY */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  IMPLEMENTATION_CHECK_MODE
//-------------------------------------------------------------------------
/* extern "C" */
/* __global__ void radixSortUnsigned32bits(const int num, uint *dat) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   const int tidx = THREADIDX_X1D; */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const int gidx = GLOBALIDX_X1D; */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (TSUB_SORT - 1); */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   __shared__ uint smem[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const uint key_src = (gidx < num) ? (dat[gidx]) : (0xffffffff); */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__  int idx_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ uint key_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = hp + lane + ii * TSUB_SORT; */
/* #ifndef TEST_SORT_ONLY_MODE */
/*     idx_src[idx] = idx; */
/* #endif//TEST_SORT_ONLY_MODE */
/*     key_src[idx] = (idx < num) ? (dat[idx]) : (0xffffffff); */
/*   } */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__  int idx_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ uint key_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*   __radixSortTsub32(32 / RADIX_SORT_CHECK_BITS, lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 		    key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 		    &key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 		    , tidx, &smem[hp] */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 		    ); */
/* #else///TEST_SORT_ONLY_MODE */
/*   __radixSortTsub32idx(32 / RADIX_SORT_CHECK_BITS, lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 		       gidx, &idx_dst[hp], key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 		       &idx_src[hp], &idx_dst[hp], &key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 		       , tidx, &smem[hp] */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 		       ); */
/* #endif//TEST_SORT_ONLY_MODE */
/*   //----------------------------------------------------------------------- */
/* #pragma unroll */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = hp + lane + ii * TSUB_SORT; */
/* /\* #ifndef TEST_SORT_ONLY_MODE *\/ */
/* /\*     idx_src[idx] = idx; *\/ */
/* /\* #endif//TEST_SORT_ONLY_MODE *\/ */
/*     dat[idx] = key_dst[idx]; */
/*   } */
/*   //----------------------------------------------------------------------- */
/* /\* #   if  SORT_ELEMENTS_PER_THREAD == 1 *\/ */
/* /\*   if( gidx < num ) *\/ */
/* /\* #ifdef  TEST_SORT_ONLY_MODE *\/ */
/* /\*     printf("%d\t%u\t%u\n", gidx, key_src, key_dst[tidx]); *\/ */
/* /\* #else///TEST_SORT_ONLY_MODE *\/ */
/* /\*     printf("%d\t%u\t%d\t%u\t%u\n", gidx, key_src, idx_dst[tidx], key_dst[tidx], dat[idx_dst[tidx]]); *\/ */
/* /\* #endif//TEST_SORT_ONLY_MODE *\/ */
/* /\* #else///SORT_ELEMENTS_PER_THREAD == 1 *\/ */
/* /\*   if( tidx == 0 ) *\/ */
/* /\*     for(int ii = 0; ii < num; ii++) *\/ */
/* /\* #ifdef  TEST_SORT_ONLY_MODE *\/ */
/* /\*       printf("%d\t%u\t%u\n", ii, dat[ii], key_dst[ii]); *\/ */
/* /\* #else///TEST_SORT_ONLY_MODE *\/ */
/* /\*       printf("%d\t%u\t%d\t%u\t%u\n", ii, dat[ii], idx_dst[ii], key_dst[ii], dat[idx_dst[ii]]); *\/ */
/* /\* #endif//TEST_SORT_ONLY_MODE *\/ */
/* /\* #endif//SORT_ELEMENTS_PER_THREAD == 1 *\/ */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
/* __global__ void radixSortSigned32bits(const int num, int *dat) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   const int tidx = THREADIDX_X1D; */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const int gidx = GLOBALIDX_X1D; */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (TSUB_SORT - 1); */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   __shared__ uint smem[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const int src = (gidx < num) ? (dat[gidx]) : INT_MAX; */
/*   uint_int tmp;  tmp.i = src; */
/*   const uint key_src = flip32int(tmp.u); */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__  int idx_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ uint key_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = hp + lane + ii * TSUB_SORT; */
/*     const int src = (idx < num) ? (dat[idx]) : INT_MAX; */
/*     uint_int tmp;  tmp.i = src; */
/* #ifndef TEST_SORT_ONLY_MODE */
/*     idx_src[idx] = idx; */
/* #endif//TEST_SORT_ONLY_MODE */
/*     key_src[idx] = flip32int(tmp.u); */
/*   } */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__  int idx_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ uint key_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*   __radixSort32bitsWarp(lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 			key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 			&key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 			, smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 			); */
/* #else///TEST_SORT_ONLY_MODE */
/*   __radixSort32bits_withidx_Warp(lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 				 gidx, &idx_dst[hp], key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 				 &idx_src[hp], &idx_dst[hp], &key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 				 , smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 				 ); */
/* #endif//TEST_SORT_ONLY_MODE */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   if( gidx < num ) */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*     printf("%d\t%d\t%u\t%d\n", gidx, src, key_dst[tidx], flip32int(key_dst[tidx])); */
/* #else///TEST_SORT_ONLY_MODE */
/*     printf("%d\t%d\t%d\t%u\t%d\t%d\n", gidx, src, idx_dst[tidx], key_dst[tidx], flip32int(key_dst[tidx]), __shfl(src, idx_dst[tidx], TSUB_SORT)); */
/* #endif//TEST_SORT_ONLY_MODE */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/*   if( lane == 0 ) */
/*     for(int ii = 0; ii < num; ii++) */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*       printf("%d\t%d\t%u\t%d\n", ii, dat[ii], key_dst[ii], flip32int(key_dst[ii])); */
/* #else///TEST_SORT_ONLY_MODE */
/*       printf("%d\t%d\t%d\t%u\t%d\t%d\n", ii, dat[ii], idx_dst[ii], key_dst[ii], flip32int(key_dst[ii]), dat[idx_dst[ii]]); */
/* #endif//TEST_SORT_ONLY_MODE */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
/* __global__ void radixSortFloating32bits(const int num, float *dat) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   const int tidx = THREADIDX_X1D; */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const int gidx = GLOBALIDX_X1D; */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (TSUB_SORT - 1); */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   __shared__ uint smem[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const float src = (gidx < num) ? (dat[gidx]) : FLT_MAX; */
/*   uint_float tmp;  tmp.f = src; */
/*   const uint key_src = flip32flt(tmp.u); */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__  int idx_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ uint key_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = hp + lane + ii * TSUB_SORT; */
/*     const float src = (idx < num) ? (dat[idx]) : FLT_MAX; */
/*     uint_float tmp;  tmp.f = src; */
/* #ifndef TEST_SORT_ONLY_MODE */
/*     idx_src[idx] = idx; */
/* #endif//TEST_SORT_ONLY_MODE */
/*     key_src[idx] = flip32flt(tmp.u); */
/*   } */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__  int idx_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ uint key_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*   __radixSort32bitsWarp(lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 			key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 			&key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 			, smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 			); */
/* #else///TEST_SORT_ONLY_MODE */
/*   __radixSort32bits_withidx_Warp(lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 				 gidx, &idx_dst[hp], key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 				 &idx_src[hp], &idx_dst[hp], &key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 				 , smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 				 ); */
/* #endif//TEST_SORT_ONLY_MODE */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   if( gidx < num ){ */
/*     tmp.u = undo32flt(key_dst[tidx]); */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*     printf("%d\t%e\t%u\t%e\n", gidx, src, key_dst[tidx], tmp.f); */
/* #else///TEST_SORT_ONLY_MODE */
/*     printf("%d\t%e\t%d\t%u\t%e\t%e\n", gidx, src, idx_dst[tidx], key_dst[tidx], tmp.f, __shfl(src, idx_dst[tidx], TSUB_SORT)); */
/* #endif//TEST_SORT_ONLY_MODE */
/*   } */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/*   if( lane == 0 ) */
/*     for(int ii = 0; ii < num; ii++){ */
/*       uint_float tmp; */
/*       tmp.u = undo32flt(key_dst[ii]); */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*       printf("%d\t%e\t%u\t%e\n", ii, dat[ii], key_dst[ii], tmp.f); */
/* #else///TEST_SORT_ONLY_MODE */
/*       printf("%d\t%e\t%d\t%u\t%e\t%e\n", ii, dat[ii], idx_dst[ii], key_dst[ii], tmp.f, dat[idx_dst[ii]]); */
/* #endif//TEST_SORT_ONLY_MODE */
/*     } */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
/* __global__ void radixSortUnsigned64bits(const int num, ulong *dat) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   const int tidx = THREADIDX_X1D; */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const int gidx = GLOBALIDX_X1D; */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (TSUB_SORT - 1); */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   __shared__ uint smem[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const ulong key_src = (gidx < num) ? (dat[gidx]) : (0xffffffffffffffff); */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__ int   idx_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ ulong key_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = hp + lane + ii * TSUB_SORT; */
/* #ifndef TEST_SORT_ONLY_MODE */
/*     idx_src[idx] = idx; */
/* #endif//TEST_SORT_ONLY_MODE */
/*     key_src[idx] = (idx < num) ? (dat[idx]) : (0xffffffffffffffff); */
/*   } */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__   int idx_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ ulong key_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*   __radixSort64bitsWarp(lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 			key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 			&key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 			, smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 			); */
/* #else///TEST_SORT_ONLY_MODE */
/*   __radixSort64bits_withidx_Warp(lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 				 gidx, &idx_dst[hp], key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 				 &idx_src[hp], &idx_dst[hp], &key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 				 , smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 				 ); */
/* #endif//TEST_SORT_ONLY_MODE */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   if( gidx < num ) */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*     printf("%d\t%lu\t%lu\n", gidx, key_src, key_dst[tidx]); */
/* #else///TEST_SORT_ONLY_MODE */
/*     printf("%d\t%lu\t%d\t%lu\t%lu\n", gidx, key_src, idx_dst[tidx], key_dst[tidx], dat[idx_dst[tidx]]); */
/* #endif//TEST_SORT_ONLY_MODE */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/*   if( lane == 0 ) */
/*     for(int ii = 0; ii < num; ii++) */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*       printf("%d\t%lu\t%lu\n", ii, dat[ii], key_dst[ii]); */
/* #else///TEST_SORT_ONLY_MODE */
/*       printf("%d\t%lu\t%d\t%lu\t%lu\n", ii, dat[ii], idx_dst[ii], key_dst[ii], dat[idx_dst[ii]]); */
/* #endif//TEST_SORT_ONLY_MODE */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
/* __global__ void radixSortSigned64bits(const int num, long *dat) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   const int tidx = THREADIDX_X1D; */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const int gidx = GLOBALIDX_X1D; */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (TSUB_SORT - 1); */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   __shared__ uint smem[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const long src = (gidx < num) ? (dat[gidx]) : LONG_MAX; */
/*   ulong_long tmp;  tmp.i = src; */
/*   const ulong key_src = flip64int(tmp.u); */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__   int idx_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ ulong key_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = hp + lane + ii * TSUB_SORT; */
/*     const long src = (idx < num) ? (dat[idx]) : LONG_MAX; */
/*     ulong_long tmp;  tmp.i = src; */
/* #ifndef TEST_SORT_ONLY_MODE */
/*     idx_src[idx] = idx; */
/* #endif//TEST_SORT_ONLY_MODE */
/*     key_src[idx] = flip64int(tmp.u); */
/*   } */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__   int idx_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ ulong key_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*   __radixSort64bitsWarp(lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 			key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 			&key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 			, smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 			); */
/* #else///TEST_SORT_ONLY_MODE */
/*   __radixSort64bits_withidx_Warp(lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 				 gidx, &idx_dst[hp], key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 				 &idx_src[hp], &idx_dst[hp], &key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 				 , smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 				 ); */
/* #endif//TEST_SORT_ONLY_MODE */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   if( gidx < num ) */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*     printf("%d\t%ld\t%lu\t%ld\n", gidx, src, key_dst[tidx], flip64int(key_dst[tidx])); */
/* #else///TEST_SORT_ONLY_MODE */
/*     printf("%d\t%ld\t%d\t%lu\t%ld\t%ld\n", gidx, src, idx_dst[tidx], key_dst[tidx], flip64int(key_dst[tidx]), __shfl(src, idx_dst[tidx], TSUB_SORT)); */
/* #endif//TEST_SORT_ONLY_MODE */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/*   if( lane == 0 ) */
/*     for(int ii = 0; ii < num; ii++) */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*       printf("%d\t%ld\t%lu\t%ld\n", ii, dat[ii], key_dst[ii], flip64int(key_dst[ii])); */
/* #else///TEST_SORT_ONLY_MODE */
/*       printf("%d\t%ld\t%d\t%lu\t%ld\t%ld\n", ii, dat[ii], idx_dst[ii], key_dst[ii], flip64int(key_dst[ii]), dat[idx_dst[ii]]); */
/* #endif//TEST_SORT_ONLY_MODE */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
/* __global__ void radixSortFloating64bits(const int num, double *dat) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   const int tidx = THREADIDX_X1D; */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const int gidx = GLOBALIDX_X1D; */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (TSUB_SORT - 1); */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   __shared__ uint smem[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   const double src = (gidx < num) ? (dat[gidx]) : DBL_MAX; */
/*   ulong_double tmp;  tmp.f = src; */
/*   const ulong key_src = flip64flt(tmp.u); */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__   int idx_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ ulong key_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = hp + lane + ii * TSUB_SORT; */
/*     const double src = (idx < num) ? (dat[idx]) : DBL_MAX; */
/*     ulong_double tmp;  tmp.f = src; */
/* #ifndef TEST_SORT_ONLY_MODE */
/*     idx_src[idx] = idx; */
/* #endif//TEST_SORT_ONLY_MODE */
/*     key_src[idx] = flip64flt(tmp.u); */
/*   } */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__   int idx_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__ ulong key_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*   __radixSort64bitsWarp(lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 			key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 			&key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 			, smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 			); */
/* #else///TEST_SORT_ONLY_MODE */
/*   __radixSort64bits_withidx_Warp(lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/* 				 gidx, &idx_dst[hp], key_src, &key_dst[hp] */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/* 				 &idx_src[hp], &idx_dst[hp], &key_src[hp], &key_dst[hp] */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 				 , smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 				 ); */
/* #endif//TEST_SORT_ONLY_MODE */
/*   //----------------------------------------------------------------------- */
/* #   if  SORT_ELEMENTS_PER_THREAD == 1 */
/*   if( gidx < num ){ */
/*     tmp.u = undo64flt(key_dst[tidx]); */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*     printf("%d\t%e\t%lu\t%e\n", gidx, src, key_dst[tidx], tmp.f); */
/* #else///TEST_SORT_ONLY_MODE */
/*     printf("%d\t%e\t%d\t%lu\t%e\t%e\n", gidx, src, idx_dst[tidx], key_dst[tidx], tmp.f, __shfl(src, idx_dst[tidx], TSUB_SORT)); */
/* #endif//TEST_SORT_ONLY_MODE */
/*   } */
/* #else///SORT_ELEMENTS_PER_THREAD == 1 */
/*   if( lane == 0 ) */
/*     for(int ii = 0; ii < num; ii++){ */
/*       ulong_double tmp; */
/*       tmp.u = undo64flt(key_dst[ii]); */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*       printf("%d\t%e\t%lu\t%e\n", ii, dat[ii], key_dst[ii], tmp.f); */
/* #else///TEST_SORT_ONLY_MODE */
/*       printf("%d\t%e\t%d\t%lu\t%e\t%e\n", ii, dat[ii], idx_dst[ii], key_dst[ii], tmp.f, dat[idx_dst[ii]]); */
/* #endif//TEST_SORT_ONLY_MODE */
/*     } */
/* #endif//SORT_ELEMENTS_PER_THREAD == 1 */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
/* extern "C" */
/* __global__ void radixSortUnsigned32bits_block(const int num, uint *dat) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   const int tidx = THREADIDX_X1D; */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (warpSize - 1); */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   __shared__ uint smem[4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/*   //----------------------------------------------------------------------- */
/*   __shared__ uint key_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   __shared__ uint key_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #ifndef TEST_SORT_ONLY_MODE */
/*   __shared__  int idx_src[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//TEST_SORT_ONLY_MODE */
/*   __shared__  int idx_dst[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = hp + lane + ii * warpSize; */
/* #ifndef TEST_SORT_ONLY_MODE */
/*     idx_src[idx] = idx; */
/* #endif//TEST_SORT_ONLY_MODE */
/*     key_src[idx] = (idx < num) ? (dat[idx]) : (0xffffffff); */
/*   } */
/*   //----------------------------------------------------------------------- */
/* #ifdef  TEST_SORT_ONLY_MODE */
/*   __radixSortBlck32(32 / RADIX_SORT_CHECK_BITS, tidx, */
/* #   if  SORT_ELEMENTS_PER_THREAD != 1 */
/* 		    hp, */
/* #endif//SORT_ELEMENTS_PER_THREAD != 1 */
/* 		    idx_dst, key_src, key_dst */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 		    , smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 		    ); */
/* #else///TEST_SORT_ONLY_MODE */
/*   __radixSortBlck32idx(32 / RADIX_SORT_CHECK_BITS, tidx, */
/* #   if  SORT_ELEMENTS_PER_THREAD != 1 */
/* 		       hp, */
/* #endif//SORT_ELEMENTS_PER_THREAD != 1 */
/* 		       idx_src, */
/* 		       idx_dst, key_src, key_dst */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/* 		       , smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/* 		       ); */
/* #endif//TEST_SORT_ONLY_MODE */
/*   //----------------------------------------------------------------------- */
/* #pragma unroll */
/*   for(int ii = 0; ii < SORT_ELEMENTS_PER_THREAD; ii++){ */
/*     const int idx = hp + lane + ii * warpSize; */
/* /\* #ifndef TEST_SORT_ONLY_MODE *\/ */
/* /\*     idx_src[idx] = idx; *\/ */
/* /\* #endif//TEST_SORT_ONLY_MODE *\/ */
/*     dat[idx] = key_dst[idx]; */
/*   } */
/*   //----------------------------------------------------------------------- */
/* /\* #   if  SORT_ELEMENTS_PER_THREAD == 1 *\/ */
/* /\*   const int gidx = GLOBALIDX_X1D; *\/ */
/* /\*   if( gidx < num ) *\/ */
/* /\* #ifdef  TEST_SORT_ONLY_MODE *\/ */
/* /\*     /\\* printf("%d\t%u\t%u\n", gidx, key_src, key_dst[tidx]); *\\/ *\/ */
/* /\*     printf("%d\t%u\t%u\n", gidx, dat[gidx], key_dst[gidx]); *\/ */
/* /\* #else///TEST_SORT_ONLY_MODE *\/ */
/* /\*     printf("%d\t%u\t%d\t%u\t%u\n", gidx, key_src, idx_dst[tidx], key_dst[tidx], dat[idx_dst[tidx]]); *\/ */
/* /\* #endif//TEST_SORT_ONLY_MODE *\/ */
/* /\* #else///SORT_ELEMENTS_PER_THREAD == 1 *\/ */
/* /\*   if( tidx == 0 ) *\/ */
/* /\*     for(int ii = 0; ii < num; ii++) *\/ */
/* /\* #ifdef  TEST_SORT_ONLY_MODE *\/ */
/* /\*       printf("%d\t%u\t%u\n", ii, dat[ii], key_dst[ii]); *\/ */
/* /\* #else///TEST_SORT_ONLY_MODE *\/ */
/* /\*       printf("%d\t%u\t%d\t%u\t%u\n", ii, dat[ii], idx_dst[ii], key_dst[ii], dat[idx_dst[ii]]); *\/ */
/* /\* #endif//TEST_SORT_ONLY_MODE *\/ */
/* /\* #endif//SORT_ELEMENTS_PER_THREAD == 1 *\/ */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
/* extern "C" */
/* __global__ void __launch_bounds__(NTHREADS_SORT, NBLOCKS_PER_SM_SORT32) radixSortUnsigned32bits_grid */
/*      (const int num, uint * RESTRICT dat, */
/*       int * RESTRICT numIterFul, int2 * RESTRICT infoFul, */
/*       int * RESTRICT numIterLoc, int2 * RESTRICT infoLoc, const int remLocMax, int * RESTRICT fail_dev, */
/*       uint * RESTRICT key_tmp_gm, */
/*       int * RESTRICT scan_gm, int * RESTRICT bucket_gm, int * RESTRICT tail_gm, */
/*       int * RESTRICT scanFul, */
/*       int * RESTRICT gsync0, int * RESTRICT gsync1 */
/*       ) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   const int bnum =   GRIDDIM_X1D; */
/*   const int bidx =  BLOCKIDX_X1D; */
/*   const int tidx = THREADIDX_X1D; */
/*   const int gidx = GLOBALIDX_X1D; */
/*   //----------------------------------------------------------------------- */
/*   const int lane = tidx & (warpSize - 1); */
/* #   if  SORT_ELEMENTS_PER_THREAD != 1 */
/*   const int head = tidx - lane; */
/*   const int hp = head * SORT_ELEMENTS_PER_THREAD; */
/* #endif//SORT_ELEMENTS_PER_THREAD != 1 */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*   /\* 4 * 4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes *\/ */
/*   __shared__ uint smem[4 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/*   //----------------------------------------------------------------------- */
/*   /\* 4 * 2 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes *\/ */
/*   __shared__  int idx_src_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   __shared__  int idx_dst_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   /\* 4 * 2 * NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD bytes *\/ */
/*   __shared__ uint key_src_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   __shared__ uint key_dst_sm[NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD]; */
/*   //----------------------------------------------------------------------- */
/* #if 0 */
/*   if( tidx == 0 ) */
/*     printf("bidx = %d, bnum = %d, tidx = %d, gidx = %d\n", bidx, bnum, tidx, gidx); */
/* #endif */
/*   //----------------------------------------------------------------------- */


/*   //----------------------------------------------------------------------- */
/*   /\* data preparation *\/ */
/*   //----------------------------------------------------------------------- */
/*   if( gidx == 0 ){ */
/*     //--------------------------------------------------------------------- */
/*     int2 info = {num, 0}; */
/*     infoFul[0] = info; */
/*     numIterFul[0] = 32 / RADIX_SORT_CHECK_BITS; */
/*     //--------------------------------------------------------------------- */
/*   }/\* if( gidx == 0 ){ *\/ */
/*   //----------------------------------------------------------------------- */
/* #if 0 */
/*   if( tidx == 0 ) */
/*     printf("bidx = %d, bnum = %d, tidx = %d, gidx = %d\n", bidx, bnum, tidx, gidx); */
/* #endif */
/*   //----------------------------------------------------------------------- */
/*   /\* data conversion is necessary for int or float *\/ */
/*   globalSync(tidx, bidx, bnum, gsync0, gsync1); */
/*   //----------------------------------------------------------------------- */
/*   __radixSortGrid32 */
/*     (numIterFul, infoFul, 1, */
/*      numIterLoc, infoLoc, remLocMax, fail_dev, */
/*      tidx, lane, */
/* #   if  SORT_ELEMENTS_PER_THREAD != 1 */
/*      hp, */
/* #endif//SORT_ELEMENTS_PER_THREAD != 1 */
/*      dat, key_tmp_gm */
/* #ifndef USE_WARP_SHUFFLE_FUNC_SORT */
/*      , smem */
/* #endif//USE_WARP_SHUFFLE_FUNC_SORT */
/*      , idx_src_sm , idx_dst_sm, key_src_sm, key_dst_sm */
/*      , scan_gm, bucket_gm, tail_gm */
/*      , scanFul */
/*      , gidx, bidx, bnum, gsync0, gsync1); */
/*   //----------------------------------------------------------------------- */
/* #if 0 */
/*   if( tidx == 0 ) */
/*     printf("bidx = %d, bnum = %d, tidx = %d, gidx = %d\n", bidx, bnum, tidx, gidx); */
/* #endif */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
int main(void)
{
  //-----------------------------------------------------------------------
  int devIdx;
  deviceInfo devInfo;
  deviceProp devProp;
  openAcceleratorGPUs(&devIdx, &devInfo, &devProp, 0, 1, 0, 1);
  //-----------------------------------------------------------------------
  const int blockNum = devProp.numSM * NBLOCKS_PER_SM_SORT32;
  //-----------------------------------------------------------------------
#ifdef  PERFORMANCE_COMPARISON_MODE
  fprintf(stderr, "blockNum = %d, Ttot = %d, Nvec = %d, Nbit = %d, SMpref = %d, useWS = %d, mode = %d, mask = %d\n", blockNum, NTHREADS_SORT, SORT_ELEMENTS_PER_THREAD, RADIX_SORT_CHECK_BITS,
#ifdef  SMEM_PREF_FOR_GLOBAL_SORT
	  1,
#else///SMEM_PREF_FOR_GLOBAL_SORT
	  0,
#endif//SMEM_PREF_FOR_GLOBAL_SORT
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	  1,
#else///USE_WARP_SHUFFLE_FUNC_SORT
	  0,
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#ifdef  USE_LSD_GRID_SORT_FUNC
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
	  4,
	  /* radixSortGridLSD_uint_1st<<<BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT), NTHREADS_SORT>>> */
	  /* radixSortGridLSD_uint_2nd<<<RADIX_SORT_NUM_BUCKET, NTHREADS_SORT_ACCUMULATION>>> */
	  /* radixSortGridLSD_uint_3rd<<<BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT), NTHREADS_SORT>>> */
#else///USE_SPLITTED_GRID_SORT_FUNC
	  3,
	  /* radixSortGridLSD_uint<<<blockNum, NTHREADS_SORT>>> */
#endif//USE_SPLITTED_GRID_SORT_FUNC
#else///USE_LSD_GRID_SORT_FUNC
#ifdef  USE_SLICED_GRID_SORT_FUNC
	  2,
	  /* radixSortGrid_mod_scatter_uint<<<bnum, NTHREADS_SORT>>> */
	  /* radix_sort_grid_mod_collect<<<RADIX_SORT_NUM_BUCKET * (*snum_hst), NTHREADS_SORT>>> */
	  /* radixSortGrid_mod_gather_uint<<<bnum, NTHREADS_SORT>>> */
#else///USE_SLICED_GRID_SORT_FUNC
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
	  1,
	  /* radixSortGrid_uint_msd<<<blockNum, NTHREADS_SORT>>> */
	  /* radixSortGrid_uint_lsd<<<*remLoc_hst, NTHREADS_SORT>>>(uint_dev, iterLoc, infoLoc); */
#else///USE_SPLITTED_GRID_SORT_FUNC
	  0,
	  /* radixSortGrid_uint       <<<blockNum, NTHREADS_SORT>>> */
#endif//USE_SPLITTED_GRID_SORT_FUNC
#endif//USE_SLICED_GRID_SORT_FUNC
#endif//USE_LSD_GRID_SORT_FUNC
#ifdef  USE_MASK_INSTEAD_OF_IF
	  1
#else///USE_MASK_INSTEAD_OF_IF
	  0
#endif//USE_MASK_INSTEAD_OF_IF
	  );
  fflush(stderr);
#endif//PERFORMANCE_COMPARISON_MODE
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  FILE *fp;
  char file[32];
#ifdef  PERFORMANCE_COMPARISON_WITH_THRUST
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
  sprintf(file, "log/bench.split.dat");
#else///USE_SPLITTED_GRID_SORT_FUNC
  sprintf(file, "log/bench.merge.dat");
#endif//USE_SPLITTED_GRID_SORT_FUNC
  fp = fopen(file, "w");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", file);
  }
  fprintf(fp, "# num\ttime(s)\t# of pairs per sec.\n");
  fprintf(fp, "# thrust_sort\tthrust_stable_sort\toriginal radix sort\n");
  fprintf(fp, "#\tNTHREADS_SORT = %d, RADIX_SORT_CHECK_BITS = %d, SORT_ELEMENTS_PER_THREAD = %d\n", NTHREADS_SORT, RADIX_SORT_CHECK_BITS, SORT_ELEMENTS_PER_THREAD);
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
  fprintf(fp, "#\tUSE_SPLITTED_GRID_SORT_FUNC is enabled\n");
#endif//USE_SPLITTED_GRID_SORT_FUNC
#ifdef  PERFORMANCE_COMPARISON_MODE
  fprintf(fp, "#\tMEASURE_ITERATION_NUM = %d\n", MEASURE_ITERATION_NUM);
#endif//PERFORMANCE_COMPARISON_MODE
  fclose(fp);
#endif//PERFORMANCE_COMPARISON_WITH_THRUST
  //-----------------------------------------------------------------------
  cudaEvent_t eventIni;  checkCudaErrors(cudaEventCreate(&eventIni));
  cudaEvent_t eventFin;  checkCudaErrors(cudaEventCreate(&eventFin));
  //-----------------------------------------------------------------------
#ifndef PERFORMANCE_COMPARISON_MODE
  const int num = SORT_ELEMENTS_PER_THREAD *     TSUB_SORT;
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * 2; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * 3; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * 4; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum + 15; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum *     2; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum *     4; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum *     8; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum *    16; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum *    32; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum *    64; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum *   128; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum *  1024; */
  /* const int num = SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum * 16384; */
  /* const int num = 1048576; */
#else///PERFORMANCE_COMPARISON_MODE
  /* const int num = SORT_ELEMENTS_PER_THREAD *     TSUB_SORT; */
  for(int num = SORT_ELEMENTS_PER_THREAD *     TSUB_SORT; num <= SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * blockNum * 1024; num <<= 1)
#endif//PERFORMANCE_COMPARISON_MODE
    {
      //-------------------------------------------------------------------
      uint   * uint_dev;  mycudaMalloc    ((void **)& uint_dev, num * sizeof(  uint));
      uint   * uint_hst;  mycudaMallocHost((void **)& uint_hst, num * sizeof(  uint));
      /* int    *  int_dev;  mycudaMalloc    ((void **)&  int_dev, num * sizeof(   int)); */
      /* int    *  int_hst;  mycudaMallocHost((void **)&  int_hst, num * sizeof(   int)); */
      /* float  *  flt_dev;  mycudaMalloc    ((void **)&  flt_dev, num * sizeof( float)); */
      /* float  *  flt_hst;  mycudaMallocHost((void **)&  flt_hst, num * sizeof( float)); */
      /* ulong  *ulong_dev;  mycudaMalloc    ((void **)&ulong_dev, num * sizeof( ulong)); */
      /* ulong  *ulong_hst;  mycudaMallocHost((void **)&ulong_hst, num * sizeof( ulong)); */
      /* long   * long_dev;  mycudaMalloc    ((void **)& long_dev, num * sizeof(  long)); */
      /* long   * long_hst;  mycudaMallocHost((void **)& long_hst, num * sizeof(  long)); */
      /* double *  dbl_dev;  mycudaMalloc    ((void **)&  dbl_dev, num * sizeof(double)); */
      /* double *  dbl_hst;  mycudaMallocHost((void **)&  dbl_hst, num * sizeof(double)); */
      //-------------------------------------------------------------------
#ifndef PERFORMANCE_COMPARISON_MODE
      sprintf(file, "log/input.dat");
      fp = fopen(file, "w");
      if( fp == NULL ){
	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", file);
      }
#endif//PERFORMANCE_COMPARISON_MODE
      //-------------------------------------------------------------------
      for(int ii = 0; ii < num; ii++){
	uint_hst[ii] = ((rand() & 0xffff) << 16) | (rand() & 0xffff);
	/* int_hst [ii] = ((rand() & 0xffff) << 16) | (rand() & 0xffff); */
	/* while( true ){ */
	/*   const uint p = ((rand() & 0xffff) << 16) | (rand() & 0xffff); */
	/*   flt_hst [ii] = *(float *)&p; */
	/*   if( isfinite(flt_hst[ii]) != 0 ) */
	/* 	break; */
	/* } */
	/* ulong_hst[ii] = ((ulong)(rand() & 0xffff) << 48) | ((ulong)(rand() & 0xffff) << 32) | ((ulong)(rand() & 0xffff) << 16) | (ulong)(rand() & 0xffff); */
	/* long_hst [ii] = (( long)(rand() & 0xffff) << 48) | (( long)(rand() & 0xffff) << 32) | (( long)(rand() & 0xffff) << 16) | ( long)(rand() & 0xffff); */
	/* while( true ){ */
	/*   const ulong p = ((ulong)(rand() & 0xffff) << 48) | ((ulong)(rand() & 0xffff) << 32) | ((ulong)(rand() & 0xffff) << 16) | (ulong)(rand() & 0xffff); */
	/*   dbl_hst [ii] = *(double *)&p; */
	/*   if( isfinite(dbl_hst[ii]) != 0 ) */
	/* 	break; */
	/* } */
#ifndef PERFORMANCE_COMPARISON_MODE
	fprintf(fp, "%u\t%x\n", uint_hst[ii], uint_hst[ii]);
#endif//PERFORMANCE_COMPARISON_MODE
      }/* for(int ii = 0; ii < num; ii++){ */
#ifndef PERFORMANCE_COMPARISON_MODE
      fclose(fp);
#endif//PERFORMANCE_COMPARISON_MODE
      //-------------------------------------------------------------------
      checkCudaErrors(cudaMemcpy( uint_dev,  uint_hst, num * sizeof(  uint), cudaMemcpyHostToDevice));
      /* checkCudaErrors(cudaMemcpy(  int_dev,   int_hst, num * sizeof(   int), cudaMemcpyHostToDevice)); */
      /* checkCudaErrors(cudaMemcpy(  flt_dev,   flt_hst, num * sizeof( float), cudaMemcpyHostToDevice)); */
      /* checkCudaErrors(cudaMemcpy(ulong_dev, ulong_hst, num * sizeof( ulong), cudaMemcpyHostToDevice)); */
      /* checkCudaErrors(cudaMemcpy( long_dev,  long_hst, num * sizeof(  long), cudaMemcpyHostToDevice)); */
      /* checkCudaErrors(cudaMemcpy(  dbl_dev,   dbl_hst, num * sizeof(double), cudaMemcpyHostToDevice)); */
      //-------------------------------------------------------------------
      /* preparation for performance measurements */
#ifdef  PERFORMANCE_COMPARISON_WITH_THRUST
      float time_thrust_sort, time_thrust_stable_sort;
#endif//PERFORMANCE_COMPARISON_WITH_THRUST
      float time_original;
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
#ifdef  PERFORMANCE_COMPARISON_WITH_THRUST
      //-------------------------------------------------------------------
      /* performance measurements for sort in CUDA thrust */
      time_thrust_sort = 0.0f;
#ifdef  PERFORMANCE_COMPARISON_MODE
      for(int ii = 0; ii < MEASURE_ITERATION_NUM; ii++)
#endif//PERFORMANCE_COMPARISON_MODE
	{
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaMemcpy( uint_dev,  uint_hst, num * sizeof(  uint), cudaMemcpyHostToDevice));
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaEventRecord(eventIni, 0));
	  thrust::sort((thrust::device_ptr<uint>)uint_dev, (thrust::device_ptr<uint>)(uint_dev + num));
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaEventRecord(eventFin, 0));
	  checkCudaErrors(cudaEventSynchronize(eventFin));
	  float time = 0.0f;
	  checkCudaErrors(cudaEventElapsedTime(&time, eventIni, eventFin));
	  //---------------------------------------------------------------
	  time_thrust_sort += time;
	  //---------------------------------------------------------------
	}
      //-------------------------------------------------------------------
      /* performance measurements for stable_sort in CUDA thrust */
      time_thrust_stable_sort = 0.0f;
#ifdef  PERFORMANCE_COMPARISON_MODE
      for(int ii = 0; ii < MEASURE_ITERATION_NUM; ii++)
#endif//PERFORMANCE_COMPARISON_MODE
	{
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaMemcpy( uint_dev,  uint_hst, num * sizeof(  uint), cudaMemcpyHostToDevice));
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaEventRecord(eventIni, 0));
	  thrust::stable_sort((thrust::device_ptr<uint>)uint_dev, (thrust::device_ptr<uint>)(uint_dev + num));
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaEventRecord(eventFin, 0));
	  checkCudaErrors(cudaEventSynchronize(eventFin));
	  float time = 0.0f;
	  checkCudaErrors(cudaEventElapsedTime(&time, eventIni, eventFin));
	  //---------------------------------------------------------------
	  time_thrust_stable_sort += time;
	  //---------------------------------------------------------------
	}
      //-------------------------------------------------------------------
#endif//PERFORMANCE_COMPARISON_WITH_THRUST
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
      /* shared memory / L1 data cache configuration for radix sort within a device */
      //-------------------------------------------------------------------
      if( num > SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT ){
	//-----------------------------------------------------------------
#ifdef  SMEM_PREF_FOR_GLOBAL_SORT
#ifdef  TEST_SORT_ONLY_MODE
#ifdef  USE_LSD_GRID_SORT_FUNC
	checkCudaErrors(cudaFuncSetCacheConfig(radixSortGridLSD_uint    , cudaFuncCachePreferShared));
#endif//USE_LSD_GRID_SORT_FUNC
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
	checkCudaErrors(cudaFuncSetCacheConfig(radixSortGrid_uint_msd, cudaFuncCachePreferShared));
	checkCudaErrors(cudaFuncSetCacheConfig(radixSortGrid_uint_lsd, cudaFuncCachePreferShared));
#else///USE_SPLITTED_GRID_SORT_FUNC
	checkCudaErrors(cudaFuncSetCacheConfig(radixSortGrid_uint    , cudaFuncCachePreferShared));
#endif//USE_SPLITTED_GRID_SORT_FUNC
#else///TEST_SORT_ONLY_MODE
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
	/* checkCudaErrors(cudaFuncSetCacheConfig(radixSortGrid_uint_by_key_msd, cudaFuncCachePreferShared)); */
	/* checkCudaErrors(cudaFuncSetCacheConfig(radixSortGrid_uint_by_key_lsd, cudaFuncCachePreferShared)); */
#else///USE_SPLITTED_GRID_SORT_FUNC
	/* checkCudaErrors(cudaFuncSetCacheConfig(radixSortGrid_uint_by_key    , cudaFuncCachePreferShared)); */
#endif//USE_SPLITTED_GRID_SORT_FUNC
#endif//TEST_SORT_ONLY_MODE
#endif//SMEM_PREF_FOR_GLOBAL_SORT
	//-----------------------------------------------------------------
      }/* if( num > SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT ){ */
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
      /* start performance measurements for original library */
      time_original = 0.0f;
      //-------------------------------------------------------------------
      if( num <= SORT_ELEMENTS_PER_THREAD *     TSUB_SORT ){
	//-----------------------------------------------------------------
	/* radix sort within a warp */
	//-----------------------------------------------------------------
#ifdef  PERFORMANCE_COMPARISON_MODE
      for(int ii = 0; ii < MEASURE_ITERATION_NUM; ii++)
#endif//PERFORMANCE_COMPARISON_MODE
	{
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaMemcpy( uint_dev,  uint_hst, num * sizeof(  uint), cudaMemcpyHostToDevice));
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaEventRecord(eventIni, 0));
#ifdef  TEST_SORT_ONLY_MODE
	  radixSortTsub_uint       <<<1, TSUB_SORT>>>(num, uint_dev);
#else///TEST_SORT_ONLY_MODE
	  radixSortTsub_uint_by_key<<<1, TSUB_SORT>>>(num, uint_dev, uint * RESTRICT dat);
#endif//TEST_SORT_ONLY_MODE
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaEventRecord(eventFin, 0));
	  checkCudaErrors(cudaEventSynchronize(eventFin));
	  float time = 0.0f;
	  checkCudaErrors(cudaEventElapsedTime(&time, eventIni, eventFin));
	  //---------------------------------------------------------------
	  time_original += time;
	  //---------------------------------------------------------------
	}
	//-----------------------------------------------------------------
      }/* if( num <= SORT_ELEMENTS_PER_THREAD *     TSUB_SORT ){ */
      //-------------------------------------------------------------------
      else{
	//-----------------------------------------------------------------
	if( num <= SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT ){
	  //---------------------------------------------------------------
	  /* radix sort within a block */
	  //---------------------------------------------------------------
#ifdef  PERFORMANCE_COMPARISON_MODE
	  for(int ii = 0; ii < MEASURE_ITERATION_NUM; ii++)
#endif//PERFORMANCE_COMPARISON_MODE
	    {
	      //-----------------------------------------------------------
	      checkCudaErrors(cudaMemcpy( uint_dev,  uint_hst, num * sizeof(  uint), cudaMemcpyHostToDevice));
	      //-----------------------------------------------------------
	      checkCudaErrors(cudaEventRecord(eventIni, 0));
#ifdef  TEST_SORT_ONLY_MODE
	      radixSortBlck_uint       <<<1, NTHREADS_SORT>>>(num, uint_dev);
#else///TEST_SORT_ONLY_MODE
	      radixSortBlck_uint_by_key<<<1, NTHREADS_SORT>>>(num, uint_dev, uint * RESTRICT dat);
#endif//TEST_SORT_ONLY_MODE
	      //-----------------------------------------------------------
	      checkCudaErrors(cudaEventRecord(eventFin, 0));
	      checkCudaErrors(cudaEventSynchronize(eventFin));
	      float time = 0.0f;
	      checkCudaErrors(cudaEventElapsedTime(&time, eventIni, eventFin));
	      //-----------------------------------------------------------
	      time_original += time;
	      //-----------------------------------------------------------
	    }
	  //---------------------------------------------------------------
	}/* if( num <= SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT ){ */
	//-----------------------------------------------------------------
	else{
	  //---------------------------------------------------------------
	  /* radix sort within a grid */
	  //---------------------------------------------------------------
	  /* preparation for global synchronization */
	  int *gsync0;      mycudaMalloc((void **)&gsync0, blockNum * sizeof(int));
	  int *gsync1;      mycudaMalloc((void **)&gsync1, blockNum * sizeof(int));
	  initGsync_kernel<<<1, blockNum>>>(blockNum, gsync0, gsync1);
	  getLastCudaError("initGsync_kernel");
#ifdef  USE_MODIFIED_GRID_SORT_FUNC
	  int *gsync0Loc;      mycudaMalloc((void **)&gsync0Loc, blockNum * sizeof(int));
	  int *gsync1Loc;      mycudaMalloc((void **)&gsync1Loc, blockNum * sizeof(int));
	  initGsync_kernel<<<1, blockNum>>>(blockNum, gsync0Loc, gsync1Loc);
#endif//USE_MODIFIED_GRID_SORT_FUNC
#   if  defined(USE_MODIFIED_GRID_SORT_FUNC) || defined(USE_SLICED_GRID_SORT_FUNC)
	  int *remFulAdd;	  mycudaMalloc    ((void **)&remFulAdd, sizeof(int));
	  int *remLocAdd;	  mycudaMalloc    ((void **)&remLocAdd, sizeof(int));
	  int *remLocHst;	  mycudaMallocHost((void **)&remLocHst, sizeof(int));
#endif//defined(USE_MODIFIED_GRID_SORT_FUNC) || defined(USE_SLICED_GRID_SORT_FUNC)
	  //---------------------------------------------------------------
	  /* memory allocation for working stack */
	  uint *dat_tmp;
	  int *dst_gm, *idx_gm, *sum_gm, *sumFul;
	  mycudaMalloc((void **)&dat_tmp, num * sizeof(uint));
#   if  !defined(USE_MODIFIED_GRID_SORT_FUNC) && !defined(USE_SLICED_GRID_SORT_FUNC)
	  mycudaMalloc((void **)&dst_gm, RADIX_SORT_NUM_BUCKET * sizeof(int));
	  mycudaMalloc((void **)&idx_gm, RADIX_SORT_NUM_BUCKET * sizeof(int));
#else///!defined(USE_MODIFIED_GRID_SORT_FUNC) && !defined(USE_SLICED_GRID_SORT_FUNC)
#ifdef  USE_SLICED_GRID_SORT_FUNC
	  mycudaMalloc((void **)&dst_gm, RADIX_SORT_NUM_BUCKET * RADIX_SORT_BOX_NUM_ALLOCATED * sizeof(int));
	  mycudaMalloc((void **)&idx_gm, RADIX_SORT_NUM_BUCKET * RADIX_SORT_BOX_NUM_ALLOCATED * sizeof(int));
#else///USE_SLICED_GRID_SORT_FUNC
	  mycudaMalloc((void **)&dst_gm, RADIX_SORT_NUM_BUCKET * blockNum * sizeof(int));
	  mycudaMalloc((void **)&idx_gm, RADIX_SORT_NUM_BUCKET * blockNum * sizeof(int));
#endif//USE_SLICED_GRID_SORT_FUNC
#endif//!defined(USE_MODIFIED_GRID_SORT_FUNC) && !defined(USE_SLICED_GRID_SORT_FUNC)
#ifdef  USE_SLICED_GRID_SORT_FUNC
	  mycudaMalloc((void **)&sum_gm, RADIX_SORT_NUM_BUCKET * RADIX_SORT_BOX_NUM_ALLOCATED * BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) * sizeof(int));
	  mycudaMalloc((void **)&sumFul, RADIX_SORT_NUM_BUCKET * RADIX_SORT_BOX_NUM_ALLOCATED * BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) * sizeof(int));
#else///USE_SLICED_GRID_SORT_FUNC
	  mycudaMalloc((void **)&sum_gm, RADIX_SORT_NUM_BUCKET * blockNum * BLOCKSIZE(num, blockNum * SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) * sizeof(int));
	  mycudaMalloc((void **)&sumFul, RADIX_SORT_NUM_BUCKET * blockNum * BLOCKSIZE(num, blockNum * SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) * sizeof(int));
	  /* mycudaMalloc((void **)&sum_gm, RADIX_SORT_NUM_BUCKET * BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) * sizeof(int)); */
	  /* mycudaMalloc((void **)&sumFul, RADIX_SORT_NUM_BUCKET * BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) * sizeof(int)); */
#endif//USE_SLICED_GRID_SORT_FUNC
	  //---------------------------------------------------------------
	  /* memory allocation for buckets used in MSD radix sort phase */
	  const int guess = RADIX_SORT_REM_LOC_MAX_GUESS << ((1 + (int)ilog2((uint)BLOCKSIZE(num, NTHREADS_SORT * SORT_ELEMENTS_PER_THREAD)) / RADIX_SORT_CHECK_BITS) * RADIX_SORT_CHECK_BITS);

#ifdef  USE_SLICED_GRID_SORT_FUNC
	  int  *snum_dev;	  mycudaMalloc    ((void **)&snum_dev, sizeof(int));
	  int  *snum_hst;	  mycudaMallocHost((void **)&snum_hst, sizeof(int));
	  /* iterFul0 and iterFul1 should be implemented; then, global synchronization function would be removed */
	  int  *iterFul0;	  mycudaMalloc((void **)&iterFul0, guess * sizeof(int));
	  int  *iterFul1;	  mycudaMalloc((void **)&iterFul1, guess * sizeof(int));
	  int2 *infoFul0;	  mycudaMalloc((void **)&infoFul0, guess * sizeof(int2));
	  int2 *infoFul1;	  mycudaMalloc((void **)&infoFul1, guess * sizeof(int2));

	  int2 *disp_dev;	  mycudaMalloc((void **)&disp_dev, RADIX_SORT_BOX_NUM_ALLOCATED * sizeof(int2));
#else///USE_SLICED_GRID_SORT_FUNC
	  int *iterFul;      mycudaMalloc((void **)&iterFul, guess * sizeof(int));      int2 *infoFul;      mycudaMalloc((void **)&infoFul, guess * sizeof(int2));
#endif//USE_SLICED_GRID_SORT_FUNC
	  int *iterLoc;      mycudaMalloc((void **)&iterLoc, guess * sizeof(int));      int2 *infoLoc;      mycudaMalloc((void **)&infoLoc, guess * sizeof(int2));
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
	  int *remLoc_dev;      mycudaMalloc    ((void **)&remLoc_dev, sizeof(int));
	  int *remLoc_hst;      mycudaMallocHost((void **)&remLoc_hst, sizeof(int));
#else///USE_SPLITTED_GRID_SORT_FUNC
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
	  int *checkedCounts;      mycudaMalloc((void **)&checkedCounts, sizeof(int));
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
#endif//USE_SPLITTED_GRID_SORT_FUNC
          //---------------------------------------------------------------
	  /* memory allocation for error detection */
	  int *fail_dev;      mycudaMalloc((void **)&fail_dev, sizeof(int));
	  int fail_hst = 0;
	  checkCudaErrors(cudaMemcpy(fail_dev, &fail_hst, sizeof(int), cudaMemcpyHostToDevice));
	  //---------------------------------------------------------------
#ifdef  PERFORMANCE_COMPARISON_MODE
	  for(int ii = 0; ii < MEASURE_ITERATION_NUM; ii++)
#endif//PERFORMANCE_COMPARISON_MODE
	    {
	      //-----------------------------------------------------------
	      checkCudaErrors(cudaMemcpy( uint_dev,  uint_hst, num * sizeof(  uint), cudaMemcpyHostToDevice));
	      //-----------------------------------------------------------
	      checkCudaErrors(cudaEventRecord(eventIni, 0));
	      //-----------------------------------------------------------
#ifdef  TEST_SORT_ONLY_MODE
#ifdef  USE_LSD_GRID_SORT_FUNC
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
#pragma unroll
	      for(int iter = 0; iter < 32 / RADIX_SORT_CHECK_BITS; iter++){
		radixSortGridLSD_uint_1st<<<BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT), NTHREADS_SORT>>>
		  (iter, num, uint_dev, dat_tmp, sum_gm, sumFul);
		radixSortGridLSD_uint_2nd<<<RADIX_SORT_NUM_BUCKET, NTHREADS_SORT_ACCUMULATION>>>
		  (BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT), idx_gm, sumFul);
		radixSortGridLSD_uint_3rd<<<BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT), NTHREADS_SORT>>>
		  (iter, num, uint_dev, dat_tmp, sum_gm, idx_gm, sumFul);
	      }/* for(int iter = 0; iter < 32 / RADIX_SORT_CHECK_BITS; iter++){ */

#else///USE_SPLITTED_GRID_SORT_FUNC
	      radixSortGridLSD_uint<<<blockNum, NTHREADS_SORT>>>
		(num, uint_dev, dat_tmp,
		 sum_gm, idx_gm, sumFul, gsync0, gsync1);
#endif//USE_SPLITTED_GRID_SORT_FUNC
#else///USE_LSD_GRID_SORT_FUNC
#ifdef  USE_SLICED_GRID_SORT_FUNC
	      //-----------------------------------------------------------
	      /* initialization */
	      radixSortGrid_kickoff<<<1, 1>>>(num, 32, iterFul0, infoFul0, iterFul1, infoFul1, remFulAdd, remLocAdd, snum_dev);
	      int remFul0 = 1;
	      int remFul1 = 0;
	      //-----------------------------------------------------------
	      /* # of blocks must satisfy: <= NTHREADS_SORT / 4, <= RADIX_SORT_BOX_NUM_ALLOCATED */
	      int bnum = BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
	      if( bnum > RADIX_SORT_BOX_NUM_ALLOCATED )		bnum = RADIX_SORT_BOX_NUM_ALLOCATED;
	      if( bnum > (NTHREADS_SORT >> 2) )		bnum = (NTHREADS_SORT >> 2);
	      //-----------------------------------------------------------
	      /* MSD radix sort phase */
	      while( remFul0 > 0 || remFul1 > 0 ){
		//---------------------------------------------------------
		/* 1st phase */
		//---------------------------------------------------------
		/* sweep iterFul0 and infoFul0 */
		while( remFul0 > 0 ){
		  radixSortGrid_mod_scatter_uint<<<bnum, NTHREADS_SORT>>>(iterFul0, infoFul0, remFul0, uint_dev, dat_tmp, snum_dev, disp_dev, sum_gm, sumFul);
		  checkCudaErrors(cudaMemcpy(snum_hst, snum_dev, sizeof(int), cudaMemcpyDeviceToHost));
		  radix_sort_grid_mod_collect<<<RADIX_SORT_NUM_BUCKET * (*snum_hst), NTHREADS_SORT>>>(dst_gm, idx_gm, sumFul, disp_dev);
		  radixSortGrid_mod_gather_uint<<<bnum, NTHREADS_SORT>>>(iterFul0, infoFul0, remFul0, iterFul1, infoFul1, remFulAdd, iterLoc, infoLoc, remLocAdd, snum_dev, *snum_hst, uint_dev, dat_tmp, sum_gm, dst_gm, idx_gm, sumFul);
		  remFul0 -= RADIX_SORT_BOX_NUM_ALLOCATED;
		}/* while( remFul0 > 0 ){ */
		radixSortGrid_finalize<<<1, 1>>>(guess, remFulAdd, remLocAdd, fail_dev);
		checkCudaErrors(cudaMemcpy(&remFul1, remFulAdd, sizeof(int), cudaMemcpyDeviceToHost));
		radixSortGrid_switch<<<1, 1>>>(remFulAdd);
		//---------------------------------------------------------

		//---------------------------------------------------------
		/* 2nd phase */
		//---------------------------------------------------------
		/* sweep iterFul1 and infoFul1 */
		while( remFul1 > 0 ){
		  radixSortGrid_mod_scatter_uint<<<bnum, NTHREADS_SORT>>>(iterFul1, infoFul1, remFul1, uint_dev, dat_tmp, snum_dev, disp_dev, sum_gm, sumFul);
		  checkCudaErrors(cudaMemcpy(snum_hst, snum_dev, sizeof(int), cudaMemcpyDeviceToHost));
		  radix_sort_grid_mod_collect<<<RADIX_SORT_NUM_BUCKET * (*snum_hst), NTHREADS_SORT>>>(dst_gm, idx_gm, sumFul, disp_dev);
		  radixSortGrid_mod_gather_uint<<<bnum, NTHREADS_SORT>>>(iterFul1, infoFul1, remFul1, iterFul0, infoFul0, remFulAdd, iterLoc, infoLoc, remLocAdd, snum_dev, *snum_hst, uint_dev, dat_tmp, sum_gm, dst_gm, idx_gm, sumFul);
		  remFul1 -= RADIX_SORT_BOX_NUM_ALLOCATED;
		}/* while( remFul1 > 0 ){ */
		radixSortGrid_finalize<<<1, 1>>>(guess, remFulAdd, remLocAdd, fail_dev);
		checkCudaErrors(cudaMemcpy(&remFul0, remFulAdd, sizeof(int), cudaMemcpyDeviceToHost));
		radixSortGrid_switch<<<1, 1>>>(remFulAdd);
		//---------------------------------------------------------
	      }/* while( remFul0 > 0 || remFul1 > 0 ){ */
	      //-----------------------------------------------------------
	      /* LSD radix sort phase */
	      checkCudaErrors(cudaMemcpy(remLocHst, remLocAdd, sizeof(int), cudaMemcpyDeviceToHost));
	      radixSortGrid_uint_lsd<<<*remLocHst, NTHREADS_SORT>>>(uint_dev, iterLoc, infoLoc);
	      //-----------------------------------------------------------
#else///USE_SLICED_GRID_SORT_FUNC
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
	      radixSortGrid_uint_msd<<<blockNum, NTHREADS_SORT>>>
		(num, uint_dev, iterFul, infoFul, iterLoc, infoLoc, guess, remLoc_dev, fail_dev,
		 dat_tmp, sum_gm, dst_gm, idx_gm, sumFul, gsync0, gsync1);
	      checkCudaErrors(cudaMemcpy(remLoc_hst, remLoc_dev, sizeof(int), cudaMemcpyDeviceToHost));
	      radixSortGrid_uint_lsd<<<*remLoc_hst, NTHREADS_SORT>>>(uint_dev, iterLoc, infoLoc);
#else///USE_SPLITTED_GRID_SORT_FUNC
	      radixSortGrid_uint       <<<blockNum, NTHREADS_SORT>>>
		(num, uint_dev,
		 iterFul, infoFul, iterLoc, infoLoc, guess,
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
		 checkedCounts,
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
		 fail_dev,
		 dat_tmp,
		 sum_gm, dst_gm, idx_gm, sumFul, gsync0, gsync1
#ifdef  USE_MODIFIED_GRID_SORT_FUNC
		 , gsync0Loc, gsync1Loc, remFulAdd, remLocAdd
#endif//USE_MODIFIED_GRID_SORT_FUNC
		 );
#endif//USE_SPLITTED_GRID_SORT_FUNC
#endif//USE_SLICED_GRID_SORT_FUNC
#endif//USE_LSD_GRID_SORT_FUNC
#else///TEST_SORT_ONLY_MODE
	      radixSortGrid_uint_by_key<<<blockNum, NTHREADS_SORT>>>
		(num, uint_dev, int * RESTRICT idx, uint * RESTRICT dat,
		 iterFul, infoFul, iterLoc, infoLoc, guess,
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
		 checkedCounts,
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
		 fail_dev,
		 dat_tmp, int * RESTRICT idx_tmp_gm,
		 sum_gm, dst_gm, idx_gm, sumFul,
		 gsync0, gsync1);
#endif//TEST_SORT_ONLY_MODE
	      //-----------------------------------------------------------
	      checkCudaErrors(cudaEventRecord(eventFin, 0));
	      checkCudaErrors(cudaEventSynchronize(eventFin));
	      float time = 0.0f;
	      checkCudaErrors(cudaEventElapsedTime(&time, eventIni, eventFin));
	      //-----------------------------------------------------------
	      time_original += time;
	      //-----------------------------------------------------------
	    }
	  //---------------------------------------------------------------
	  /* error detection */
	  checkCudaErrors(cudaMemcpy(&fail_hst, fail_dev, sizeof(int), cudaMemcpyDeviceToHost));
	  if( fail_hst != 0 ){
	    __KILL__(stderr, "ERROR: remLoc(%d) exceeds the allocated size of %d.\nPLEASE re-simulate after increasing the value ``RADIX_SORT_REM_LOC_MAX_GUESS''(%d) defined in src/sort/radix_dev.h\n", fail_hst, guess, RADIX_SORT_REM_LOC_MAX_GUESS);
	  }/* if( fail_hst != 0 ){ */
	  //---------------------------------------------------------------
	  /* memory deallocation */
	  mycudaFree(gsync0);      mycudaFree(gsync1);
#ifdef  USE_MODIFIED_GRID_SORT_FUNC
	  mycudaFree(gsync0Loc);	  mycudaFree(gsync1Loc);
#endif//USE_MODIFIED_GRID_SORT_FUNC
#   if  defined(USE_MODIFIED_GRID_SORT_FUNC) || defined(USE_SLICED_GRID_SORT_FUNC)
	  mycudaFree    (remFulAdd);
	  mycudaFree    (remLocAdd);
	  mycudaFreeHost(remLocHst);
#endif//defined(USE_MODIFIED_GRID_SORT_FUNC) || defined(USE_SLICED_GRID_SORT_FUNC)
	  mycudaFree(dat_tmp);
	  mycudaFree(dst_gm);      mycudaFree(idx_gm);
	  mycudaFree(sum_gm);      mycudaFree(sumFul);
#ifdef  USE_SLICED_GRID_SORT_FUNC
	  mycudaFree(snum_dev);	  mycudaFreeHost(snum_hst);
	  mycudaFree(iterFul0);	  mycudaFree(infoFul0);
	  mycudaFree(iterFul1);	  mycudaFree(infoFul1);
	  mycudaFree(disp_dev);
#else///USE_SLICED_GRID_SORT_FUNC
	  mycudaFree(iterFul);      mycudaFree(infoFul);
#endif//USE_SLICED_GRID_SORT_FUNC
	  mycudaFree(iterLoc);      mycudaFree(infoLoc);
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
	  mycudaFree(remLoc_dev);      mycudaFreeHost(remLoc_hst);
#else///USE_SPLITTED_GRID_SORT_FUNC
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
	  mycudaFree(checkedCounts);
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
#endif//USE_SPLITTED_GRID_SORT_FUNC
	  mycudaFree(fail_dev);
	  //---------------------------------------------------------------
	}/* else{ */
	//-----------------------------------------------------------------
      }/* else{ */
      //-------------------------------------------------------------------
      /* check sorted results using CUDA thrust */
      bool success = thrust::is_sorted((thrust::device_ptr<uint>)uint_dev, (thrust::device_ptr<uint>)(uint_dev + num));
      if( !success )
	fprintf(stderr, "ERROR: sorted result is not match with CUDA thrust for num = %d\n", num);
      //-------------------------------------------------------------------
      /* unit of the measured time is milli-second */
#ifdef  PERFORMANCE_COMPARISON_WITH_THRUST
      time_thrust_sort        *= 1.0e-3f;
      time_thrust_stable_sort *= 1.0e-3f;
#endif//PERFORMANCE_COMPARISON_WITH_THRUST
      time_original           *= 1.0e-3f;
#ifdef  PERFORMANCE_COMPARISON_MODE
#ifdef  PERFORMANCE_COMPARISON_WITH_THRUST
      time_thrust_sort        /= (float)MEASURE_ITERATION_NUM;
      time_thrust_stable_sort /= (float)MEASURE_ITERATION_NUM;
#endif//PERFORMANCE_COMPARISON_WITH_THRUST
      time_original           /= (float)MEASURE_ITERATION_NUM;
#endif//PERFORMANCE_COMPARISON_MODE
      //-------------------------------------------------------------------
#ifdef  PERFORMANCE_COMPARISON_WITH_THRUST
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
      sprintf(file, "log/bench.split.dat");
#else///USE_SPLITTED_GRID_SORT_FUNC
      sprintf(file, "log/bench.merge.dat");
#endif//USE_SPLITTED_GRID_SORT_FUNC
#else///PERFORMANCE_COMPARISON_WITH_THRUST
      sprintf(file, "bench/time.dat");
#endif//PERFORMANCE_COMPARISON_WITH_THRUST
      fp = fopen(file, "a");
      if( fp == NULL ){
	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", file);
      }
#ifdef  PERFORMANCE_COMPARISON_WITH_THRUST
      fprintf(fp, "%d\t%e\t%e\t%e\t%e\t%e\t%e\n", num,
	      time_thrust_sort       , (float)num / time_thrust_sort,
	      time_thrust_stable_sort, (float)num / time_thrust_stable_sort,
	      time_original          , (float)num / time_original);
#else///PERFORMANCE_COMPARISON_WITH_THRUST
      /*  1. # of elements */
      /*  2. elapsed time (sec) */
      /*  3. sorting rate (pairs per time) */
      /*  4. # of elements per thread (SORT_ELEMENTS_PER_THREAD) */
      /*  5. # of threads per block (NTHREADS_SORT) */
      /*  6. # of bits per step (RADIX_SORT_CHECK_BITS) */
      /*  7. SMEM_PREF_FOR_GLOBAL_SORT (1 is SM pref.) */
      /*  8. USE_WARP_SHUFFLE_FUNC_SORT (1 uses warp shuffle instructions) */
      /*  9. mode (USE_LSD_GRID_SORT_FUNC etc.) */
      /* 10. USE_MASK_INSTEAD_OF_IF (1 uses mask instead of if().) */
      fprintf(fp, "%d\t%e\t%e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
	      num, time_original, (float)num / time_original,
	      SORT_ELEMENTS_PER_THREAD, NTHREADS_SORT, RADIX_SORT_CHECK_BITS,
#ifdef  SMEM_PREF_FOR_GLOBAL_SORT
	      1,
#else///SMEM_PREF_FOR_GLOBAL_SORT
	      0,
#endif//SMEM_PREF_FOR_GLOBAL_SORT
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	      1,
#else///USE_WARP_SHUFFLE_FUNC_SORT
	      0,
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#ifdef  USE_LSD_GRID_SORT_FUNC
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
	      4,
	      /* radixSortGridLSD_uint_1st<<<BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT), NTHREADS_SORT>>> */
	      /* radixSortGridLSD_uint_2nd<<<RADIX_SORT_NUM_BUCKET, NTHREADS_SORT_ACCUMULATION>>> */
	      /* radixSortGridLSD_uint_3rd<<<BLOCKSIZE(num, SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT), NTHREADS_SORT>>> */
#else///USE_SPLITTED_GRID_SORT_FUNC
	      3,
	      /* radixSortGridLSD_uint<<<blockNum, NTHREADS_SORT>>> */
#endif//USE_SPLITTED_GRID_SORT_FUNC
#else///USE_LSD_GRID_SORT_FUNC
#ifdef  USE_SLICED_GRID_SORT_FUNC
	      2,
	      /* radixSortGrid_mod_scatter_uint<<<bnum, NTHREADS_SORT>>> */
	      /* radix_sort_grid_mod_collect<<<RADIX_SORT_NUM_BUCKET * (*snum_hst), NTHREADS_SORT>>> */
	      /* radixSortGrid_mod_gather_uint<<<bnum, NTHREADS_SORT>>> */
#else///USE_SLICED_GRID_SORT_FUNC
#ifdef  USE_SPLITTED_GRID_SORT_FUNC
	      1,
	      /* radixSortGrid_uint_msd<<<blockNum, NTHREADS_SORT>>> */
	      /* radixSortGrid_uint_lsd<<<*remLoc_hst, NTHREADS_SORT>>>(uint_dev, iterLoc, infoLoc); */
#else///USE_SPLITTED_GRID_SORT_FUNC
	      0,
	      /* radixSortGrid_uint       <<<blockNum, NTHREADS_SORT>>> */
#endif//USE_SPLITTED_GRID_SORT_FUNC
#endif//USE_SLICED_GRID_SORT_FUNC
#endif//USE_LSD_GRID_SORT_FUNC
#ifdef  USE_MASK_INSTEAD_OF_IF
	      1
#else///USE_MASK_INSTEAD_OF_IF
	      0
#endif//USE_MASK_INSTEAD_OF_IF
	      );
#endif//PERFORMANCE_COMPARISON_WITH_THRUST
      fclose(fp);
      //-------------------------------------------------------------------
      checkCudaErrors(cudaMemcpy( uint_hst,  uint_dev, num * sizeof(  uint), cudaMemcpyDeviceToHost));
#ifndef PERFORMANCE_COMPARISON_MODE
      sprintf(file, "log/output.dat");
      fp = fopen(file, "w");
      if( fp == NULL ){
	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", file);
      }
      for(int ii = 0; ii < num; ii++)
	fprintf(fp, "%u\t%x\n", uint_hst[ii], uint_hst[ii]);
      fclose(fp);
#endif//PERFORMANCE_COMPARISON_MODE
      //-------------------------------------------------------------------
      mycudaFree( uint_dev);  mycudaFreeHost( uint_hst);
      /* mycudaFree(  int_dev);  mycudaFreeHost(  int_hst); */
      /* mycudaFree(  flt_dev);  mycudaFreeHost(  flt_hst); */
      /* mycudaFree(ulong_dev);  mycudaFreeHost(ulong_hst); */
      /* mycudaFree( long_dev);  mycudaFreeHost( long_hst); */
      /* mycudaFree(  dbl_dev);  mycudaFreeHost(  dbl_hst); */
      //-------------------------------------------------------------------
    }
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaEventDestroy(eventIni));
  checkCudaErrors(cudaEventDestroy(eventFin));
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//IMPLEMENTATION_CHECK_MODE
//-------------------------------------------------------------------------





//-------------------------------------------------------------------------
#endif//RADIX_DEV_CU
//-------------------------------------------------------------------------
