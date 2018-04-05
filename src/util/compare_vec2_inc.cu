/**
 * @file compare_vec2_inc.cu
 *
 * @brief Source code for comparing values in 2-components vector on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/04/03 (Tue)
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
#include "../util/vector_inc.cu"
#include "../util/comparison_inc.cu"
#include "../util/compare_vec2_inc.cuh"


#ifndef COMPARE_VEC2_INC_CU_MULTI_CALL
#define COMPARE_VEC2_INC_CU_MULTI_CALL
/**
 * @fn getVec2MinWarp
 *
 * @brief Get minimum value within a warp.
 * @detail implicit synchronization within 32 threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type getVec2MinWarp
(Type val
#ifndef USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC
 , volatile Type * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC

  Type tmp, loc;

  loc = val.x;
  tmp = __SHFL_XOR(loc,  1, warpSize);  loc = getMinVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  2, warpSize);  loc = getMinVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  4, warpSize);  loc = getMinVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  8, warpSize);  loc = getMinVal(loc, tmp);
  tmp = __SHFL_XOR(loc, 16, warpSize);  loc = getMinVal(loc, tmp);
  val.x = __SHFL(loc, 0, warpSize);

  loc = val.y;
  tmp = __SHFL_XOR(loc,  1, warpSize);  loc = getMinVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  2, warpSize);  loc = getMinVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  4, warpSize);  loc = getMinVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  8, warpSize);  loc = getMinVal(loc, tmp);
  tmp = __SHFL_XOR(loc, 16, warpSize);  loc = getMinVal(loc, tmp);
  val.y = __SHFL(loc, 0, warpSize);

#else///USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC

  Type tmp;

  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  1]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  2]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  4]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  8]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^ 16]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
  val = ldvec(smem[head]);

#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC

  return (val);
}


/**
 * @fn getVec2MaxWarp
 *
 * @brief Get maximum value within a warp.
 * @detail implicit synchronization within 32 threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type getVec2MaxWarp
(Type val
#ifndef USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC
 , volatile Type * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC

  Type tmp, loc;

  loc = val.x;
  tmp = __SHFL_XOR(loc,  1, warpSize);  loc = getMaxVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  2, warpSize);  loc = getMaxVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  4, warpSize);  loc = getMaxVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  8, warpSize);  loc = getMaxVal(loc, tmp);
  tmp = __SHFL_XOR(loc, 16, warpSize);  loc = getMaxVal(loc, tmp);
  val.x = __SHFL(loc, 0, warpSize);

  loc = val.y;
  tmp = __SHFL_XOR(loc,  1, warpSize);  loc = getMaxVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  2, warpSize);  loc = getMaxVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  4, warpSize);  loc = getMaxVal(loc, tmp);
  tmp = __SHFL_XOR(loc,  8, warpSize);  loc = getMaxVal(loc, tmp);
  tmp = __SHFL_XOR(loc, 16, warpSize);  loc = getMaxVal(loc, tmp);
  val.y = __SHFL(loc, 0, warpSize);

#else///USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC

  Type tmp;

  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  1]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  2]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  4]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  8]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^ 16]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
  val = ldvec(smem[head]);

#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC

  return (val);
}
#endif//COMPARE_VEC2_INC_CU_MULTI_CALL


/**
 * @fn GET_VEC2_MIN_BLCK
 *
 * @brief Get minimum value within a block.
 */
template <typename Type>
__device__ __forceinline__ Type GET_VEC2_MIN_BLCK(Type val, volatile Type * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  Type ret = getVec2MinWarp(val, smem, tidx, head);


  /** 2. reduction among warps */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_VEC2_INC <= 1024 --> NTHREADS_COMPARE_VEC2_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_COMPARE_VEC2_INC >> 5) ){
    val = ldvec(smem[tidx * warpSize]);
    Type tmp;

    stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  2
    tmp = ldvec(smem[tidx ^  1]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  4
    tmp = ldvec(smem[tidx ^  2]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  8
    tmp = ldvec(smem[tidx ^  4]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >= 16
    tmp = ldvec(smem[tidx ^  8]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) == 32
    tmp = ldvec(smem[tidx ^ 16]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  2

  }/* if( tidx < (NTHREADS_COMPARE_VEC2_INC >> 5) ){ */
  __syncthreads();

  ret = ldvec(smem[0]);
  __syncthreads();

  return (ret);
}


/**
 * @fn GET_VEC2_MAX_BLCK
 *
 * @brief Get maximum value within a block.
 */
template <typename Type>
__device__ __forceinline__ Type GET_VEC2_MAX_BLCK(Type val, volatile Type * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  Type ret = getVec2MaxWarp(val, smem, tidx, head);


  /** 2. reduction among warps */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_VEC2_INC <= 1024 --> NTHREADS_COMPARE_VEC2_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_COMPARE_VEC2_INC >> 5) ){
    val = ldvec(smem[tidx * warpSize]);
    Type tmp;

    stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  2
    tmp = ldvec(smem[tidx ^  1]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  4
    tmp = ldvec(smem[tidx ^  2]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  8
    tmp = ldvec(smem[tidx ^  4]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >= 16
    tmp = ldvec(smem[tidx ^  8]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) == 32
    tmp = ldvec(smem[tidx ^ 16]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  2
  }/* if( tidx < (NTHREADS_COMPARE_VEC2_INC >> 5) ){ */
  __syncthreads();

  ret = ldvec(smem[0]);
  __syncthreads();

  return (ret);
}


/**
 * @fn GET_VEC2_MIN_MAX_BLCK
 *
 * @brief Get minimum and maximum values within a block.
 */
template <typename Type>
  __device__ __forceinline__
  void GET_VEC2_MIN_MAX_BLCK
(Type * __restrict__ min, Type * __restrict__ max, volatile Type * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  *min = getVec2MinWarp(*min, smem, tidx, head);
  *max = getVec2MaxWarp(*max, smem, tidx, head);


  /** 2. reduction among warps */
  __syncthreads();
  if( tidx == head ){
    stvec((Type *)smem,                                     head >> 5 , *min);/**< := smem[                                   (tidx / warpSize)] = min; */
    stvec((Type *)smem, (NTHREADS_COMPARE_VEC2_INC >> 1) + (head >> 5), *max);/**< := smem[(NTHREADS_COMPARE_VEC2_INC >> 1) + (tidx / warpSize)] = max; */
  }/* if( tidx == head ){ */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_VEC2_INC <= 1024 --> NTHREADS_COMPARE_VEC2_INC >> 5 <= 32 = warpSize */
  /** get minimum */
  if( tidx < (NTHREADS_COMPARE_VEC2_INC >> 5) ){
    Type val = ldvec(smem[tidx]);
    Type tmp;

#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  2
    tmp = ldvec(smem[tidx ^  1]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  4
    tmp = ldvec(smem[tidx ^  2]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  8
    tmp = ldvec(smem[tidx ^  4]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >= 16
    tmp = ldvec(smem[tidx ^  8]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) == 32
    tmp = ldvec(smem[tidx ^ 16]);  val.x = getMinVal(val.x, tmp.x);  val.y = getMinVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  2
  }/* if( tidx < (NTHREADS_COMPARE_VEC2_INC >> 5) ){ */

  /** get maximum */
  if( (head >= (NTHREADS_COMPARE_VEC2_INC >> 1)) && ((tidx - head) < (NTHREADS_COMPARE_VEC2_INC >> 5)) ){
    Type val = ldvec(smem[tidx]);
    Type tmp;

#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  2
    tmp = ldvec(smem[tidx ^  1]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  4
    tmp = ldvec(smem[tidx ^  2]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >=  8
    tmp = ldvec(smem[tidx ^  4]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) >= 16
    tmp = ldvec(smem[tidx ^  8]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_VEC2_INC >> 5) == 32
    tmp = ldvec(smem[tidx ^ 16]);  val.x = getMaxVal(val.x, tmp.x);  val.y = getMaxVal(val.y, tmp.y);  stvec((Type *)smem, tidx, val);
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_VEC2_INC >> 5) >=  2
  }/* if( (head >= (NTHREADS_COMPARE_VEC2_INC >> 1)) && ((tidx - head) < (NTHREADS_COMPARE_VEC2_INC >> 5)) ){ */

  __syncthreads();

  *min = ldvec(smem[                             0]);
  *max = ldvec(smem[NTHREADS_COMPARE_VEC2_INC >> 1]);

  __syncthreads();
}


/**
 * @fn GET_VEC2_MIN_GRID
 *
 * @brief Get minimum value within a grid.
 */
template <typename Type>
__device__ __forceinline__ Type GET_VEC2_MIN_GRID
(Type val, volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  Type ret = GET_VEC2_MIN_BLCK(val, smem, tidx, head);


  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** stvec data on the global memory */
    if( tidx == 0 )
      stvec((Type *)gmem, bidx, ret);

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get minimum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_VEC2_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_VEC2_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldvec(gmem[target]) : ret);

	/** calculate local reduction */
	tmp = GET_VEC2_MIN_BLCK(tmp, smem, tidx, head);
	ret.x = getMinVal(ret.x, tmp.x);
	ret.y = getMinVal(ret.y, tmp.y);

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stvec((Type *)gmem, bidx, ret);
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx == 0 )
      stvec((Type *)smem, 0, ldvec(gmem[0]));
    __syncthreads();

    /** upload obtained result */
    ret = ldvec(smem[0]);
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (ret);
}


/**
 * @fn GET_VEC2_MAX_GRID
 *
 * @brief Get maximum value within a grid.
 */
template <typename Type>
__device__ __forceinline__ Type GET_VEC2_MAX_GRID
(Type val, volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  Type ret = GET_VEC2_MAX_BLCK(val, smem, tidx, head);


  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** stvec data on the global memory */
    if( tidx == 0 )
      stvec((Type *)gmem, bidx, ret);

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get maximum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_VEC2_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_VEC2_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldvec(gmem[target]) : ret);

	/** calculate local reduction */
	tmp = GET_VEC2_MAX_BLCK(tmp, smem, tidx, head);
	ret.x = getMaxVal(ret.x, tmp.x);
	ret.y = getMaxVal(ret.y, tmp.y);

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stvec((Type *)gmem, bidx, ret);
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx == 0 )
      stvec((Type *)smem, 0, ldvec(gmem[0]));
    __syncthreads();

    /** upload obtained result */
    ret = ldvec(smem[0]);
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (ret);
}


/**
 * @fn GET_VEC2_MIN_MAX_GRID
 *
 * @brief Get minimum and maximum values within a grid.
 */
template <typename Type>
__device__ __forceinline__
void GET_VEC2_MIN_MAX_GRID
(Type * __restrict__ min, Type * __restrict__ max, volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem_min, volatile Type * __restrict__ gmem_max, const int bidx, const int bnum,
 int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  GET_VEC2_MIN_MAX_BLCK(min, max, smem, tidx, head);

  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** stvec data on the global memory */
    if( tidx == 0 ){
      stvec((Type *)gmem_min, bidx, *min);
      stvec((Type *)gmem_max, bidx, *max);
    }/* if( tidx == 0 ){ */

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get minimum by a representative block */
    if( bidx == 0 ){
      Type ret = *min;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_VEC2_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_VEC2_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldvec(gmem_min[target]) : ret);

	/** calculate local reduction */
	tmp = GET_VEC2_MIN_BLCK(tmp, smem, tidx, head);
	ret.x = getMinVal(ret.x, tmp.x);
	ret.y = getMinVal(ret.y, tmp.y);

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stvec((Type *)gmem_min, 0, ret);
    }/* if( bidx == 0 ){ */


    /** get maximum by a representative block */
    if( bidx == (bnum - 1) ){
      Type ret = *max;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_VEC2_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_VEC2_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldvec(gmem_max[target]) : ret);

	/** calculate local reduction */
	tmp = GET_VEC2_MAX_BLCK(tmp, smem, tidx, head);
	ret.x = getMaxVal(ret.x, tmp.x);
	ret.y = getMaxVal(ret.y, tmp.y);

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stvec((Type *)gmem_max, 0, ret);
    }/* if( bidx == (bnum - 1) ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 2 ){
      stvec((Type *)smem, tidx, (tidx == 0) ? ldvec(gmem_min[0]) : ldvec(gmem_max[0]));
    }/* if( tidx < 2 ){ */
    __syncthreads();

    /** upload obtained result */
    *min = ldvec(smem[0]);
    *max = ldvec(smem[1]);
    __syncthreads();
  }/* if( bnum > 1 ){ */
}
