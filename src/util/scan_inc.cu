/**
 * @file scan_inc.cu
 *
 * @brief Source code for parallel prefix sum library on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/11/30 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <helper_cuda.h>

#include "macro.h"
#include "cudalib.h"

#include "../util/gsync_dev.cu"
#include "../util/scan_inc.cuh"


#ifndef SCAN_INC_CU_MULTI_CALL
#define SCAN_INC_CU_MULTI_CALL

#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && !defined(_COOPERATIVE_GROUPS_H_)
#include <cooperative_groups.h>
using namespace cooperative_groups;
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && !defined(_COOPERATIVE_GROUPS_H_)

/**
 * @fn prefixSumWarp
 *
 * @brief Get parallel (inclusive) prefix sum within a warp.
 * @detail implicit synchronization within 32 threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type prefixSumWarp
(Type val, const int lane
#ifndef USE_WARP_SHUFFLE_FUNC_SCAN_INC
 , volatile Type * smem, const int tidx
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_INC

  Type tmp;
  tmp = __SHFL_UP(SHFL_MASK_32, val,  1, warpSize);  if( lane >=  1 )    val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  tmp = __SHFL_UP(SHFL_MASK_32, val,  2, warpSize);  if( lane >=  2 )    val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  tmp = __SHFL_UP(SHFL_MASK_32, val,  4, warpSize);  if( lane >=  4 )    val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  tmp = __SHFL_UP(SHFL_MASK_32, val,  8, warpSize);  if( lane >=  8 )    val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  tmp = __SHFL_UP(SHFL_MASK_32, val, 16, warpSize);  if( lane >= 16 )    val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

#else///USE_WARP_SHUFFLE_FUNC_SCAN_INC

  smem[tidx] = val;  if( lane >=  1 )    val += smem[tidx -  1];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;  if( lane >=  2 )    val += smem[tidx -  2];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;  if( lane >=  4 )    val += smem[tidx -  4];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;  if( lane >=  8 )    val += smem[tidx -  8];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;  if( lane >= 16 )    val += smem[tidx - 16];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;

#endif//USE_WARP_SHUFFLE_FUNC_SCAN_INC

  return (val);
}


/**
 * @fn totalSumWarp
 *
 * @brief Get total sum within a warp.
 * @detail implicit synchronization within 32 threads (a warp) is assumed.
 */
#ifdef  USE_WARP_REDUCE_FUNC_SCAN_INC
__device__ __forceinline__      int totalSumWarp(     int val){  return (__reduce_add_sync(SHFL_MASK_32, val));}
__device__ __forceinline__ unsigned totalSumWarp(unsigned val){  return (__reduce_add_sync(SHFL_MASK_32, val));}
#endif//USE_WARP_REDUCE_FUNC_SCAN_INC
template <typename Type>
__device__ __forceinline__ Type totalSumWarp
(Type val
#ifndef USE_WARP_SHUFFLE_FUNC_SCAN_INC
 , volatile Type * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_INC

  Type tmp;
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  1, warpSize);  val += tmp;
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  2, warpSize);  val += tmp;
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  4, warpSize);  val += tmp;
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  8, warpSize);  val += tmp;
  tmp = __SHFL_XOR(SHFL_MASK_32, val, 16, warpSize);  val += tmp;
  val = __SHFL(SHFL_MASK_32, val, 0, warpSize);

#else///USE_WARP_SHUFFLE_FUNC_SCAN_INC

  smem[tidx] = val;
  val += smem[tidx ^  1];  smem[tidx] = val;
  val += smem[tidx ^  2];  smem[tidx] = val;
  val += smem[tidx ^  4];  smem[tidx] = val;
  val += smem[tidx ^  8];  smem[tidx] = val;
  val += smem[tidx ^ 16];  smem[tidx] = val;
  val = smem[head];

#endif//USE_WARP_SHUFFLE_FUNC_SCAN_INC

  return (val);
}
#endif//SCAN_INC_CU_MULTI_CALL


/**
 * @fn PREFIX_SUM_BLCK
 *
 * @brief Get parallel (inclusive) prefix sum within a block.
 */
template <typename Type>
__device__ __forceinline__ Type PREFIX_SUM_BLCK(Type val, volatile Type * __restrict__ smem, const int lane, const int tidx)
{
  /** 1. prefix sum within warp */
  Type scan = prefixSumWarp(val, lane
#ifndef USE_WARP_SHUFFLE_FUNC_SCAN_INC
			 , smem, tidx
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_INC
			 );
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_INC
  smem[tidx] = scan;
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_INC


  /** 2. prefix sum about tail of each warp */
#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (NTHREADS_SCAN_INC > 32)
  thread_block_tile<(NTHREADS_SCAN_INC >> 5)> tile = tiled_partition<(NTHREADS_SCAN_INC >> 5)>(this_thread_block());
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (NTHREADS_SCAN_INC > 32)
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_SCAN_INC <= 1024 --> NTHREADS_SCAN_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_SCAN_INC >> 5) ){
    val = smem[tidx * warpSize + warpSize - 1];

#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_INC
    const int groupSize = NTHREADS_SCAN_INC >> 5;
    Type tmp;
#   if  (NTHREADS_SCAN_INC >> 5) >=  2
    tmp = __SHFL_UP(SHFL_MASK_SCAN_INC, val,  1, groupSize);    if( lane >=  1 )      val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  (NTHREADS_SCAN_INC >> 5) >=  4
    tmp = __SHFL_UP(SHFL_MASK_SCAN_INC, val,  2, groupSize);    if( lane >=  2 )      val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  (NTHREADS_SCAN_INC >> 5) >=  8
    tmp = __SHFL_UP(SHFL_MASK_SCAN_INC, val,  4, groupSize);    if( lane >=  4 )      val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  (NTHREADS_SCAN_INC >> 5) >= 16
    tmp = __SHFL_UP(SHFL_MASK_SCAN_INC, val,  8, groupSize);    if( lane >=  8 )      val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  (NTHREADS_SCAN_INC >> 5) == 32
    tmp = __SHFL_UP(SHFL_MASK_SCAN_INC, val, 16, groupSize);    if( lane >= 16 )      val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#endif//(NTHREADS_SCAN_INC >> 5) == 32
#endif//(NTHREADS_SCAN_INC >> 5) >= 16
#endif//(NTHREADS_SCAN_INC >> 5) >=  8
#endif//(NTHREADS_SCAN_INC >> 5) >=  4
#endif//(NTHREADS_SCAN_INC >> 5) >=  2
    smem[tidx] = val;
#else///USE_WARP_SHUFFLE_FUNC_SCAN_INC
    smem[tidx] = val;
#   if  (NTHREADS_SCAN_INC >> 5) >=  2
    if( lane >=  1 )      val += smem[tidx -  1];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    smem[tidx] = val;
#   if  (NTHREADS_SCAN_INC >> 5) >=  4
    if( lane >=  2 )      val += smem[tidx -  2];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    smem[tidx] = val;
#   if  (NTHREADS_SCAN_INC >> 5) >=  8
    if( lane >=  4 )      val += smem[tidx -  4];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    smem[tidx] = val;
#   if  (NTHREADS_SCAN_INC >> 5) >= 16
    if( lane >=  8 )      val += smem[tidx -  8];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    smem[tidx] = val;
#   if  (NTHREADS_SCAN_INC >> 5) == 32
    if( lane >= 16 )      val += smem[tidx - 16];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    smem[tidx] = val;
#endif//(NTHREADS_SCAN_INC >> 5) == 32
#endif//(NTHREADS_SCAN_INC >> 5) >= 16
#endif//(NTHREADS_SCAN_INC >> 5) >=  8
#endif//(NTHREADS_SCAN_INC >> 5) >=  4
#endif//(NTHREADS_SCAN_INC >> 5) >=  2
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_INC
  }/* if( tidx < (NTHREADS_SCAN_INC >> 5) ){ */
  __syncthreads();


  /** 3. prefix sum within a block */
  /** warpSize = 32 = 2^5 */
  if( tidx >= warpSize )
    scan += smem[(tidx >> 5) - 1];
  __syncthreads();


  /** 4. upload calculated prefix sum */
  smem[tidx] = scan;
  __syncthreads();

  return (scan);
}


/**
 * @fn PREFIX_SUM_GRID
 *
 * @brief Get parallel (inclusive) prefix sum within a grid.
 *
 * @return (scanNum) total sum of val
 */
template <typename Type>
__device__ __forceinline__ Type PREFIX_SUM_GRID
(Type val, volatile Type * __restrict__ smem, const int lane, const int tidx,
 Type * __restrict__ scanNum,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  Type scan = PREFIX_SUM_BLCK(val, smem, lane, tidx);
  *scanNum = smem[NTHREADS_SCAN_INC - 1];

  /** global prefix sum is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local prefix sum via global memory */
    /** store data on the global memory */
    if( tidx == (NTHREADS_SCAN_INC - 1) )
      gmem[bidx] = scan;

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** calculate global prefix sum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_SCAN_INC);
      Type head = 0;
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_SCAN_INC;

	/** load from the global memory */
	Type pidx = ((target < bnum) ? gmem[target] : 0);

	/** calculate local prefix sum */
	PREFIX_SUM_BLCK(pidx, smem, lane, tidx);

	/** store to the global memory */
	if( target < bnum )
	  gmem[target] = head + smem[tidx];

	head += smem[NTHREADS_SCAN_INC - 1];
	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 2 )
      smem[tidx] = (tidx == 0) ? ((bidx > 0) ? gmem[bidx - 1] : 0) : (gmem[bnum - 1]);
    __syncthreads();

    scan += smem[0];
    *scanNum = smem[1];
    __syncthreads();

    /** upload calculated prefix sum */
    smem[tidx] = scan;
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (scan);
}


/**
 * @fn PREFIX_SUM_GRID_WITH_PARTITION
 *
 * @brief Get parallel (inclusive) prefix sum within a grid.
 *
 * @return (headIdx) head index of the scan in the corresponding block
 * @return (scanNum) total sum of val
 */
template <typename Type>
__device__ __forceinline__ Type PREFIX_SUM_GRID_WITH_PARTITION
(Type val, volatile Type * __restrict__ smem, const int lane, const int tidx,
 Type * __restrict__ headIdx, Type * __restrict__ scanNum,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  Type scan = PREFIX_SUM_BLCK(val, smem, lane, tidx);
  *headIdx = 0;
  *scanNum = smem[NTHREADS_SCAN_INC - 1];

  /** global prefix sum is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local prefix sum via global memory */
    /** store data on the global memory */
    if( tidx == (NTHREADS_SCAN_INC - 1) )
      gmem[bidx] = scan;

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** calculate global prefix sum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_SCAN_INC);
      Type head = 0;
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_SCAN_INC;

	/** load from the global memory */
	Type pidx = ((target < bnum) ? gmem[target] : 0);

	/** calculate local prefix sum */
	PREFIX_SUM_BLCK(pidx, smem, lane, tidx);

	/** store to the global memory */
	if( target < bnum )
	  gmem[target] = head + smem[tidx];

	head += smem[NTHREADS_SCAN_INC - 1];
	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 2 )
      smem[tidx] = (tidx == 0) ? ((bidx > 0) ? gmem[bidx - 1] : 0) : (gmem[bnum - 1]);
    __syncthreads();

    *scanNum = smem[1];
    *headIdx = smem[0];
    scan += (*headIdx);
    __syncthreads();

    /** upload calculated prefix sum */
    smem[tidx] = scan;
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (scan);
}


/**
 * @fn TOTAL_SUM_BLCK
 *
 * @brief Get total sum within a block.
 */
#ifdef  USE_WARP_REDUCE_FUNC_SCAN_INC
__device__ __forceinline__ int TOTAL_SUM_BLCK(int val, volatile int * __restrict__ smem, const int tidx, const int head)
{
  /** 1. total sum within a warp */
  val = totalSumWarp(val);
  if( tidx == head )
    smem[tidx] = val;


  /** 2. reduction of partial sum */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_SCAN_INC <= 1024 --> NTHREADS_SCAN_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_SCAN_INC >> 5) ){
    val = smem[tidx * warpSize];

    const uint mask = __activemask();
    smem[tidx] = __reduce_add_sync(mask, val);
  }/* if( tidx < (NTHREADS_SCAN_INC >> 5) ){ */
  __syncthreads();
  val = smem[0];
  __syncthreads();

  return (val);
}
__device__ __forceinline__ unsigned TOTAL_SUM_BLCK(unsigned val, volatile unsigned * __restrict__ smem, const int tidx, const int head)
{
  /** 1. total sum within a warp */
  val = totalSumWarp(val);
  if( tidx == head )
    smem[tidx] = val;


  /** 2. reduction of partial sum */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_SCAN_INC <= 1024 --> NTHREADS_SCAN_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_SCAN_INC >> 5) ){
    val = smem[tidx * warpSize];

    const uint mask = __activemask();
    smem[tidx] = __reduce_add_sync(mask, val);
  }/* if( tidx < (NTHREADS_SCAN_INC >> 5) ){ */
  __syncthreads();
  val = smem[0];
  __syncthreads();

  return (val);
}
#endif//USE_WARP_REDUCE_FUNC_SCAN_INC
template <typename Type>
__device__ __forceinline__ Type TOTAL_SUM_BLCK(Type val, volatile Type * __restrict__ smem, const int tidx, const int head)
{
  /** 1. total sum within a warp */
  val = totalSumWarp(val
#ifndef USE_WARP_SHUFFLE_FUNC_SCAN_INC
			   , smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_INC
			 );
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_INC
  if( tidx == head )
    smem[tidx] = val;
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_INC


  /** 2. reduction of partial sum */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_SCAN_INC <= 1024 --> NTHREADS_SCAN_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_SCAN_INC >> 5) ){
    val = smem[tidx * warpSize];

#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_INC
    const int groupSize = NTHREADS_SCAN_INC >> 5;
    Type tmp;
#   if  (NTHREADS_SCAN_INC >> 5) >=  2
    tmp = __SHFL_XOR(SHFL_MASK_SCAN_INC, val,  1, groupSize);    val += tmp;
#   if  (NTHREADS_SCAN_INC >> 5) >=  4
    tmp = __SHFL_XOR(SHFL_MASK_SCAN_INC, val,  2, groupSize);    val += tmp;
#   if  (NTHREADS_SCAN_INC >> 5) >=  8
    tmp = __SHFL_XOR(SHFL_MASK_SCAN_INC, val,  4, groupSize);    val += tmp;
#   if  (NTHREADS_SCAN_INC >> 5) >= 16
    tmp = __SHFL_XOR(SHFL_MASK_SCAN_INC, val,  8, groupSize);    val += tmp;
#   if  (NTHREADS_SCAN_INC >> 5) == 32
    tmp = __SHFL_XOR(SHFL_MASK_SCAN_INC, val, 16, groupSize);    val += tmp;
#endif//(NTHREADS_SCAN_INC >> 5) == 32
#endif//(NTHREADS_SCAN_INC >> 5) >= 16
#endif//(NTHREADS_SCAN_INC >> 5) >=  8
#endif//(NTHREADS_SCAN_INC >> 5) >=  4
#endif//(NTHREADS_SCAN_INC >> 5) >=  2
    if( tidx == 0 )
      smem[0] = val;
#else///USE_WARP_SHUFFLE_FUNC_SCAN_INC
    smem[tidx] = val;
#   if  (NTHREADS_SCAN_INC >> 5) >=  2
    val += smem[tidx ^  1];    smem[tidx] = val;
#   if  (NTHREADS_SCAN_INC >> 5) >=  4
    val += smem[tidx ^  2];    smem[tidx] = val;
#   if  (NTHREADS_SCAN_INC >> 5) >=  8
    val += smem[tidx ^  4];    smem[tidx] = val;
#   if  (NTHREADS_SCAN_INC >> 5) >= 16
    val += smem[tidx ^  8];    smem[tidx] = val;
#   if  (NTHREADS_SCAN_INC >> 5) == 32
    val += smem[tidx ^ 16];    smem[tidx] = val;
#endif//(NTHREADS_SCAN_INC >> 5) == 32
#endif//(NTHREADS_SCAN_INC >> 5) >= 16
#endif//(NTHREADS_SCAN_INC >> 5) >=  8
#endif//(NTHREADS_SCAN_INC >> 5) >=  4
#endif//(NTHREADS_SCAN_INC >> 5) >=  2
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_INC
  }/* if( tidx < (NTHREADS_SCAN_INC >> 5) ){ */
  __syncthreads();
  val = smem[0];
  __syncthreads();

  return (val);
}


/**
 * @fn TOTAL_SUM_GRID
 *
 * @brief Get total sum within a grid.
 */
template <typename Type>
__device__ __forceinline__ Type TOTAL_SUM_GRID
(Type val, volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  val = TOTAL_SUM_BLCK(val, smem, tidx, head);

  /** global prefix sum is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local prefix sum via global memory */
    /** store data on the global memory */
    if( tidx == 0 )
      gmem[bidx] = val;

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** calculate global prefix sum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_SCAN_INC);
      Type sum = 0;
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_SCAN_INC;

	/** load from the global memory */
	Type subset = ((target < bnum) ? gmem[target] : 0);

	/** calculate partial sum */
	sum += TOTAL_SUM_BLCK(subset, smem, tidx, head);

	__syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */
      if( tidx == 0 )
	gmem[0] = sum;
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx == 0 )
      smem[tidx] = gmem[0];
    __syncthreads();

    /** upload calculated total sum */
    val = smem[0];
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (val);
}
