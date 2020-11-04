/**
 * @file compare_inc.cu
 *
 * @brief Source code for comparing values on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/11/04 (Wed)
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
#include "../util/compare_inc.cuh"
#include "../util/comparison_inc.cu"


#ifndef COMPARE_INC_CU_MULTI_CALL
#define COMPARE_INC_CU_MULTI_CALL

#   if  (GPUGEN >= 70) && !defined(_COOPERATIVE_GROUPS_H_)
#include <cooperative_groups.h>
using namespace cooperative_groups;
#endif//(GPUGEN >= 70) && !defined(_COOPERATIVE_GROUPS_H_)


#ifdef  USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
__device__ __forceinline__ uint flipFP32(const uint src){  uint mask = -int(src >> 31)   | 0x80000000;  return (src ^ mask);}
__device__ __forceinline__ uint undoFP32(const uint src){  uint mask = ((src >> 31) - 1) | 0x80000000;  return (src ^ mask);}
#endif//USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC


/**
 * @fn getMinWarp
 *
 * @brief Get minimum value within a warp.
 * @detail implicit synchronization within 32 threads (a warp) is assumed.
 */
#ifdef  USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
__device__ __forceinline__      int getMinWarp(     int val){  return (__reduce_min_sync(SHFL_MASK_32, val));}
__device__ __forceinline__ unsigned getMinWarp(unsigned val){  return (__reduce_min_sync(SHFL_MASK_32, val));}
__device__ __forceinline__    float getMinWarp(   float val){
  union {uint u; float f;} tmp;
  tmp.f = val;
  tmp.u = undoFP32(getMinWarp(flipFP32(tmp.u)));
  return (tmp.f);
}
#endif//USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
template <typename Type>
__device__ __forceinline__ Type getMinWarp
(Type val
#ifndef USE_WARP_SHUFFLE_FUNC_COMPARE_INC
 , volatile Type * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_INC

  Type tmp;
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  1, warpSize);  val = getMinVal(val, tmp);
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  2, warpSize);  val = getMinVal(val, tmp);
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  4, warpSize);  val = getMinVal(val, tmp);
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  8, warpSize);  val = getMinVal(val, tmp);
  tmp = __SHFL_XOR(SHFL_MASK_32, val, 16, warpSize);  val = getMinVal(val, tmp);
  val = __SHFL(SHFL_MASK_32, val, 0, warpSize);

#else///USE_WARP_SHUFFLE_FUNC_COMPARE_INC

  smem[tidx] = val;
  val = getMinVal(val, smem[tidx ^  1]);  smem[tidx] = val;
  val = getMinVal(val, smem[tidx ^  2]);  smem[tidx] = val;
  val = getMinVal(val, smem[tidx ^  4]);  smem[tidx] = val;
  val = getMinVal(val, smem[tidx ^  8]);  smem[tidx] = val;
  val = getMinVal(val, smem[tidx ^ 16]);  smem[tidx] = val;
  val = smem[head];

#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_INC

  return (val);
}


/**
 * @fn getMaxWarp
 *
 * @brief Get maximum value within a warp.
 * @detail implicit synchronization within 32 threads (a warp) is assumed.
 */
#ifdef  USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
__device__ __forceinline__      int getMaxWarp(     int val){  return (__reduce_max_sync(SHFL_MASK_32, val));}
__device__ __forceinline__ unsigned getMaxWarp(unsigned val){  return (__reduce_max_sync(SHFL_MASK_32, val));}
__device__ __forceinline__    float getMaxWarp(   float val){
  union {uint u; float f;} tmp;
  tmp.f = val;
  tmp.u = undoFP32(getMaxWarp(flipFP32(tmp.u)));
  return (tmp.f);
}
#endif//USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
template <typename Type>
__device__ __forceinline__ Type getMaxWarp
(Type val
#ifndef USE_WARP_SHUFFLE_FUNC_COMPARE_INC
 , volatile Type * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_INC

  Type tmp;
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  1, warpSize);  val = getMaxVal(val, tmp);
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  2, warpSize);  val = getMaxVal(val, tmp);
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  4, warpSize);  val = getMaxVal(val, tmp);
  tmp = __SHFL_XOR(SHFL_MASK_32, val,  8, warpSize);  val = getMaxVal(val, tmp);
  tmp = __SHFL_XOR(SHFL_MASK_32, val, 16, warpSize);  val = getMaxVal(val, tmp);
  val = __SHFL(SHFL_MASK_32, val, 0, warpSize);

#else///USE_WARP_SHUFFLE_FUNC_COMPARE_INC

  smem[tidx] = val;
  val = getMaxVal(val, smem[tidx ^  1]);  smem[tidx] = val;
  val = getMaxVal(val, smem[tidx ^  2]);  smem[tidx] = val;
  val = getMaxVal(val, smem[tidx ^  4]);  smem[tidx] = val;
  val = getMaxVal(val, smem[tidx ^  8]);  smem[tidx] = val;
  val = getMaxVal(val, smem[tidx ^ 16]);  smem[tidx] = val;
  val = smem[head];

#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_INC

  return (val);
}


/**
 * @fn getMinlocWarp
 *
 * @brief Get minimum value with location within a warp.
 * @detail implicit synchronization within 32 threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type getMinlocWarp
(Type val, volatile Type * smem, const int tidx, const int head)
{
  Type tmp;
  stloc((Type *)smem, tidx, val);  tmp = ldloc(smem[tidx ^  1]);  if( tmp.val < val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  stloc((Type *)smem, tidx, val);  tmp = ldloc(smem[tidx ^  2]);  if( tmp.val < val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  stloc((Type *)smem, tidx, val);  tmp = ldloc(smem[tidx ^  4]);  if( tmp.val < val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  stloc((Type *)smem, tidx, val);  tmp = ldloc(smem[tidx ^  8]);  if( tmp.val < val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  stloc((Type *)smem, tidx, val);  tmp = ldloc(smem[tidx ^ 16]);  if( tmp.val < val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  stloc((Type *)smem, tidx, val);

  val = ldloc(smem[head]);
  return (val);
}


/**
 * @fn getMaxlocWarp
 *
 * @brief Get maximum value with location within a warp.
 * @detail implicit synchronization within 32 threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type getMaxlocWarp
(Type val, volatile Type * smem, const int tidx, const int head)
{
  Type tmp;

  stloc((Type *)smem, tidx, val);  tmp = ldloc(smem[tidx ^  1]);  if( tmp.val > val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  stloc((Type *)smem, tidx, val);  tmp = ldloc(smem[tidx ^  2]);  if( tmp.val > val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  stloc((Type *)smem, tidx, val);  tmp = ldloc(smem[tidx ^  4]);  if( tmp.val > val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  stloc((Type *)smem, tidx, val);  tmp = ldloc(smem[tidx ^  8]);  if( tmp.val > val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  stloc((Type *)smem, tidx, val);  tmp = ldloc(smem[tidx ^ 16]);  if( tmp.val > val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  stloc((Type *)smem, tidx, val);

  val = ldloc(smem[head]);

  return (val);
}
#endif//COMPARE_INC_CU_MULTI_CALL


/**
 * @fn GET_MIN_BLCK
 *
 * @brief Get minimum value within a block.
 */
#ifdef  USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
__device__ __forceinline__ int GET_MIN_BLCK(int val, volatile int * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  int ret = getMinWarp(val);
  smem[tidx] = ret;


  /** 2. reduction among warps */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_INC <= 1024 --> NTHREADS_COMPARE_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_COMPARE_INC >> 5) ){
    val = smem[tidx * warpSize];

    const uint mask = __activemask();
    smem[tidx] = __reduce_min_sync(mask, val);
  }/* if( tidx < (NTHREADS_COMPARE_INC >> 5) ){ */
  __syncthreads();

  ret = smem[0];
  __syncthreads();

  return (ret);
}
__device__ __forceinline__ unsigned GET_MIN_BLCK(unsigned val, volatile unsigned * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  unsigned ret = getMinWarp(val);
  smem[tidx] = ret;


  /** 2. reduction among warps */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_INC <= 1024 --> NTHREADS_COMPARE_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_COMPARE_INC >> 5) ){
    val = smem[tidx * warpSize];

    const uint mask = __activemask();
    smem[tidx] = __reduce_min_sync(mask, val);
  }/* if( tidx < (NTHREADS_COMPARE_INC >> 5) ){ */
  __syncthreads();

  ret = smem[0];
  __syncthreads();

  return (ret);
}
__device__ __forceinline__ float GET_MIN_BLCK(float val, volatile float * __restrict__ smem, const int tidx, const int head)
{
  union {uint u; float f;} tmp;
  tmp.f = val;
  tmp.u = undoFP32(GET_MIN_BLCK(flipFP32(tmp.u), (unsigned *)smem, tidx, head));
  return (tmp.f);
}
#endif//USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
template <typename Type>
__device__ __forceinline__ Type GET_MIN_BLCK(Type val, volatile Type * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  Type ret = getMinWarp(val
#ifndef USE_WARP_SHUFFLE_FUNC_COMPARE_INC
			, smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_INC
			);
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_INC
  smem[tidx] = ret;
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_INC


  /** 2. reduction among warps */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_INC <= 1024 --> NTHREADS_COMPARE_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_COMPARE_INC >> 5) ){
    val = smem[tidx * warpSize];

#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_INC
    const int groupSize = NTHREADS_COMPARE_INC >> 5;
    Type tmp;
#   if  (NTHREADS_COMPARE_INC >> 5) >=  2
    tmp = __SHFL_XOR(SHFL_MASK_COMPARE_INC, val,  1, groupSize);    val = getMinVal(val, tmp);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  4
    tmp = __SHFL_XOR(SHFL_MASK_COMPARE_INC, val,  2, groupSize);    val = getMinVal(val, tmp);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  8
    tmp = __SHFL_XOR(SHFL_MASK_COMPARE_INC, val,  4, groupSize);    val = getMinVal(val, tmp);
#   if  (NTHREADS_COMPARE_INC >> 5) >= 16
    tmp = __SHFL_XOR(SHFL_MASK_COMPARE_INC, val,  8, groupSize);    val = getMinVal(val, tmp);
#   if  (NTHREADS_COMPARE_INC >> 5) == 32
    tmp = __SHFL_XOR(SHFL_MASK_COMPARE_INC, val, 16, groupSize);    val = getMinVal(val, tmp);
#endif//(NTHREADS_COMPARE_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_INC >> 5) >=  2
    smem[tidx] = val;
#else///USE_WARP_SHUFFLE_FUNC_COMPARE_INC
    smem[tidx] = val;
#   if  (NTHREADS_COMPARE_INC >> 5) >=  2
    val = getMinVal(val, smem[tidx ^  1]);    smem[tidx] = val;
#   if  (NTHREADS_COMPARE_INC >> 5) >=  4
    val = getMinVal(val, smem[tidx ^  2]);    smem[tidx] = val;
#   if  (NTHREADS_COMPARE_INC >> 5) >=  8
    val = getMinVal(val, smem[tidx ^  4]);    smem[tidx] = val;
#   if  (NTHREADS_COMPARE_INC >> 5) >= 16
    val = getMinVal(val, smem[tidx ^  8]);    smem[tidx] = val;
#   if  (NTHREADS_COMPARE_INC >> 5) == 32
    val = getMinVal(val, smem[tidx ^ 16]);    smem[tidx] = val;
#endif//(NTHREADS_COMPARE_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_INC >> 5) >=  2
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_INC
  }/* if( tidx < (NTHREADS_COMPARE_INC >> 5) ){ */
  __syncthreads();

  ret = smem[0];
  __syncthreads();

  return (ret);
}


/**
 * @fn GET_MAX_BLCK
 *
 * @brief Get maximum value within a block.
 */
#ifdef  USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
__device__ __forceinline__ int GET_MAX_BLCK(int val, volatile int * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  int ret = getMaxWarp(val);
  smem[tidx] = ret;


  /** 2. reduction among warps */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_INC <= 1024 --> NTHREADS_COMPARE_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_COMPARE_INC >> 5) ){
    val = smem[tidx * warpSize];

    const uint mask = __activemask();
    smem[tidx] = __reduce_max_sync(mask, val);
  }/* if( tidx < (NTHREADS_COMPARE_INC >> 5) ){ */
  __syncthreads();

  ret = smem[0];
  __syncthreads();

  return (ret);
}
__device__ __forceinline__ unsigned GET_MAX_BLCK(unsigned val, volatile unsigned * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  unsigned ret = getMaxWarp(val);
  smem[tidx] = ret;


  /** 2. reduction among warps */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_INC <= 1024 --> NTHREADS_COMPARE_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_COMPARE_INC >> 5) ){
    val = smem[tidx * warpSize];

    const uint mask = __activemask();
    smem[tidx] = __reduce_max_sync(mask, val);
  }/* if( tidx < (NTHREADS_COMPARE_INC >> 5) ){ */
  __syncthreads();

  ret = smem[0];
  __syncthreads();

  return (ret);
}
__device__ __forceinline__ float GET_MAX_BLCK(float val, volatile float * __restrict__ smem, const int tidx, const int head)
{
  union {uint u; float f;} tmp;
  tmp.f = val;
  tmp.u = undoFP32(GET_MAX_BLCK(flipFP32(tmp.u), (unsigned *)smem, tidx, head));
  return (tmp.f);
}
#endif//USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
template <typename Type>
__device__ __forceinline__ Type GET_MAX_BLCK(Type val, volatile Type * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  Type ret = getMaxWarp(val
#ifndef USE_WARP_SHUFFLE_FUNC_COMPARE_INC
			, smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_INC
			);
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_INC
  smem[tidx] = ret;
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_INC


  /** 2. reduction among warps */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_INC <= 1024 --> NTHREADS_COMPARE_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_COMPARE_INC >> 5) ){
    val = smem[tidx * warpSize];

#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_INC
    const int groupSize = NTHREADS_COMPARE_INC >> 5;
    Type tmp;
#   if  (NTHREADS_COMPARE_INC >> 5) >=  2
    tmp = __SHFL_XOR(SHFL_MASK_COMPARE_INC, val,  1, groupSize);    val = getMaxVal(val, tmp);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  4
    tmp = __SHFL_XOR(SHFL_MASK_COMPARE_INC, val,  2, groupSize);    val = getMaxVal(val, tmp);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  8
    tmp = __SHFL_XOR(SHFL_MASK_COMPARE_INC, val,  4, groupSize);    val = getMaxVal(val, tmp);
#   if  (NTHREADS_COMPARE_INC >> 5) >= 16
    tmp = __SHFL_XOR(SHFL_MASK_COMPARE_INC, val,  8, groupSize);    val = getMaxVal(val, tmp);
#   if  (NTHREADS_COMPARE_INC >> 5) == 32
    tmp = __SHFL_XOR(SHFL_MASK_COMPARE_INC, val, 16, groupSize);    val = getMaxVal(val, tmp);
#endif//(NTHREADS_COMPARE_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_INC >> 5) >=  2
    smem[tidx] = val;
#else///USE_WARP_SHUFFLE_FUNC_COMPARE_INC
    smem[tidx] = val;
#   if  (NTHREADS_COMPARE_INC >> 5) >=  2
    val = getMaxVal(val, smem[tidx ^  1]);    smem[tidx] = val;
#   if  (NTHREADS_COMPARE_INC >> 5) >=  4
    val = getMaxVal(val, smem[tidx ^  2]);    smem[tidx] = val;
#   if  (NTHREADS_COMPARE_INC >> 5) >=  8
    val = getMaxVal(val, smem[tidx ^  4]);    smem[tidx] = val;
#   if  (NTHREADS_COMPARE_INC >> 5) >= 16
    val = getMaxVal(val, smem[tidx ^  8]);    smem[tidx] = val;
#   if  (NTHREADS_COMPARE_INC >> 5) == 32
    val = getMaxVal(val, smem[tidx ^ 16]);    smem[tidx] = val;
#endif//(NTHREADS_COMPARE_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_INC >> 5) >=  2
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_INC
  }/* if( tidx < (NTHREADS_COMPARE_INC >> 5) ){ */
  __syncthreads();

  ret = smem[0];
  __syncthreads();

  return (ret);
}


/**
 * @fn GET_MINLOC_BLCK
 *
 * @brief Get minimum value with location within a block.
 */
template <typename Type>
__device__ __forceinline__ Type GET_MINLOC_BLCK(Type val, volatile Type * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  Type ret = getMinlocWarp(val, smem, tidx, head);

  /** 2. reduction among warps */
#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (NTHREADS_COMPARE_INC > 32)
  thread_block_tile<NTHREADS_COMPARE_INC >> 5> tile = tiled_partition<NTHREADS_COMPARE_INC >> 5>(this_thread_block());
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (NTHREADS_COMPARE_INC > 32)
  __syncthreads();
  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_INC <= 1024 --> NTHREADS_COMPARE_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_COMPARE_INC >> 5) ){
    val = ldloc(smem[tidx * warpSize]);
    stloc((Type *)smem, tidx, val);

#   if  (NTHREADS_COMPARE_INC >> 5) >=  2
    Type tmp;
    tmp = ldloc(smem[tidx ^  1]);    if( tmp.val < val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  4
    tmp = ldloc(smem[tidx ^  2]);    if( tmp.val < val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  8
    tmp = ldloc(smem[tidx ^  4]);    if( tmp.val < val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >= 16
    tmp = ldloc(smem[tidx ^  8]);    if( tmp.val < val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) == 32
    tmp = ldloc(smem[tidx ^ 16]);    if( tmp.val < val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#endif//(NTHREADS_COMPARE_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_INC >> 5) >=  2

  }/* if( tidx < (NTHREADS_COMPARE_INC >> 5) ){ */
  __syncthreads();

  ret = ldloc(smem[0]);
  __syncthreads();

  return (ret);
}


/**
 * @fn GET_MAXLOC_BLCK
 *
 * @brief Get maximum value with location within a block.
 */
template <typename Type>
__device__ __forceinline__ Type GET_MAXLOC_BLCK(Type val, volatile Type * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  Type ret = getMaxlocWarp(val, smem, tidx, head);

  /** 2. reduction among warps */
#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (NTHREADS_COMPARE_INC > 32)
  thread_block_tile<NTHREADS_COMPARE_INC >> 5> tile = tiled_partition<NTHREADS_COMPARE_INC >> 5>(this_thread_block());
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (NTHREADS_COMPARE_INC > 32)
  __syncthreads();
  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_INC <= 1024 --> NTHREADS_COMPARE_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_COMPARE_INC >> 5) ){
    val = ldloc(smem[tidx * warpSize]);
    stloc((Type *)smem, tidx, val);

#   if  (NTHREADS_COMPARE_INC >> 5) >=  2
    Type tmp;
    tmp = ldloc(smem[tidx ^  1]);    if( tmp.val > val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  4
    tmp = ldloc(smem[tidx ^  2]);    if( tmp.val > val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  8
    tmp = ldloc(smem[tidx ^  4]);    if( tmp.val > val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >= 16
    tmp = ldloc(smem[tidx ^  8]);    if( tmp.val > val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) == 32
    tmp = ldloc(smem[tidx ^ 16]);    if( tmp.val > val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#endif//(NTHREADS_COMPARE_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_INC >> 5) >=  2

  }/* if( tidx < (NTHREADS_COMPARE_INC >> 5) ){ */
  __syncthreads();

  ret = ldloc(smem[0]);
  __syncthreads();

  return (ret);
}


/**
 * @fn GET_MINLOC_MAXLOC_BLCK
 *
 * @brief Get minimum and maximum values with their locations within a block.
 */
template <typename Type>
__device__ __forceinline__ void GET_MINLOC_MAXLOC_BLCK(Type * __restrict__ minloc, Type * __restrict__ maxloc, volatile Type * __restrict__ smem, const int tidx, const int head)
{
  /** 1. reduction within warp */
  *minloc = getMinlocWarp(*minloc, smem, tidx, head);
  *maxloc = getMaxlocWarp(*maxloc, smem, tidx, head);

  /** 2. reduction among warps */
  /** warpSize = 32 = 2^5; NTHREADS_COMPARE_INC <= 1024 --> NTHREADS_COMPARE_INC >> 5 <= 32 = warpSize */
#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (NTHREADS_COMPARE_INC > 32)
  thread_block_tile<NTHREADS_COMPARE_INC >> 5> tile = tiled_partition<NTHREADS_COMPARE_INC >> 5>(this_thread_block());
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (NTHREADS_COMPARE_INC > 32)
  __syncthreads();
  if( tidx == head ){
    stloc((Type *)smem,                                head >> 5 , *minloc);/**< := smem[                              (tidx / warpSize)] = minloc; */
    stloc((Type *)smem, (NTHREADS_COMPARE_INC >> 1) + (head >> 5), *maxloc);/**< := smem[(NTHREADS_COMPARE_INC >> 1) + (tidx / warpSize)] = maxloc; */
  }/* if( tidx == head ){ */
  __syncthreads();

  /** get minimum */
  if( tidx < (NTHREADS_COMPARE_INC >> 5) ){
    Type val = ldloc(smem[tidx]);

#   if  (NTHREADS_COMPARE_INC >> 5) >=  2
    Type tmp;
    tmp = ldloc(smem[tidx ^  1]);    if( tmp.val < val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  4
    tmp = ldloc(smem[tidx ^  2]);    if( tmp.val < val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  8
    tmp = ldloc(smem[tidx ^  4]);    if( tmp.val < val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >= 16
    tmp = ldloc(smem[tidx ^  8]);    if( tmp.val < val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) == 32
    tmp = ldloc(smem[tidx ^ 16]);    if( tmp.val < val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#endif//(NTHREADS_COMPARE_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_INC >> 5) >=  2
  }/* if( tidx < (NTHREADS_COMPARE_INC >> 5) ){ */

  /** get maximum */
  if( (head >= (NTHREADS_COMPARE_INC >> 1)) && ((tidx - head) < (NTHREADS_COMPARE_INC >> 5)) ){
    Type val = ldloc(smem[tidx]);

#   if  (NTHREADS_COMPARE_INC >> 5) >=  2
    Type tmp;
    tmp = ldloc(smem[tidx ^  1]);    if( tmp.val > val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  4
    tmp = ldloc(smem[tidx ^  2]);    if( tmp.val > val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >=  8
    tmp = ldloc(smem[tidx ^  4]);    if( tmp.val > val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) >= 16
    tmp = ldloc(smem[tidx ^  8]);    if( tmp.val > val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#   if  (NTHREADS_COMPARE_INC >> 5) == 32
    tmp = ldloc(smem[tidx ^ 16]);    if( tmp.val > val.val )      val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    tile.sync();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    stloc((Type *)smem, tidx, val);
#endif//(NTHREADS_COMPARE_INC >> 5) == 32
#endif//(NTHREADS_COMPARE_INC >> 5) >= 16
#endif//(NTHREADS_COMPARE_INC >> 5) >=  8
#endif//(NTHREADS_COMPARE_INC >> 5) >=  4
#endif//(NTHREADS_COMPARE_INC >> 5) >=  2
  }/* if( (head >= (NTHREADS_COMPARE_INC >> 1)) && ((tidx - head) < (NTHREADS_COMPARE_INC >> 5)) ){ */

  __syncthreads();

  *minloc = ldloc(smem[                        0]);
  *maxloc = ldloc(smem[NTHREADS_COMPARE_INC >> 1]);

  __syncthreads();
}


/**
 * @fn GET_MIN_GRID
 *
 * @brief Get minimum value within a grid.
 */
#ifdef  USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
__device__ __forceinline__ float GET_MIN_GRID(float val, volatile float * __restrict__ smem, const int tidx, const int head,
					      volatile float * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  union {uint u; float f;} tmp;
  tmp.f = val;
  tmp.u = undoFP32(GET_MIN_GRID(flipFP32(tmp.u), (unsigned *)smem, tidx, head, (unsigned *)gmem, bidx, bnum, gsync0, gsync1));
  return (tmp.f);
}
#endif//USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
template <typename Type>
__device__ __forceinline__ Type GET_MIN_GRID
(Type val, volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  Type ret = GET_MIN_BLCK(val, smem, tidx, head);


  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 )
      gmem[bidx] = ret;

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** calculate global prefix sum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? gmem[target] : ret);

	/** calculate local reduction */
	tmp = GET_MIN_BLCK(tmp, smem, tidx, head);
	ret = getMinVal(ret, tmp);

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	gmem[0] = ret;
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx == 0 )
      smem[0] = gmem[0];
    __syncthreads();

    /** upload obtained result */
    ret = smem[0];
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (ret);
}


/**
 * @fn GET_MAX_GRID
 *
 * @brief Get maximum value within a grid.
 */
#ifdef  USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
__device__ __forceinline__ float GET_MAX_GRID(float val, volatile float * __restrict__ smem, const int tidx, const int head,
					      volatile float * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  union {uint u; float f;} tmp;
  tmp.f = val;
  tmp.u = undoFP32(GET_MAX_GRID(flipFP32(tmp.u), (unsigned *)smem, tidx, head, (unsigned *)gmem, bidx, bnum, gsync0, gsync1));
  return (tmp.f);
}
#endif//USE_WARP_REDUCE_FUNCTIONS_COMPARE_INC
template <typename Type>
__device__ __forceinline__ Type GET_MAX_GRID
(Type val, volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  Type ret = GET_MAX_BLCK(val, smem, tidx, head);


  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 )
      gmem[bidx] = ret;

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** calculate global prefix sum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? gmem[target] : ret);

	/** calculate local reduction */
	tmp = GET_MAX_BLCK(tmp, smem, tidx, head);
	ret = getMaxVal(ret, tmp);

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	gmem[0] = ret;
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx == 0 )
      smem[0] = gmem[0];
    __syncthreads();

    /** upload obtained result */
    ret = smem[0];
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (ret);
}


/**
 * @fn GET_MINLOC_GRID
 *
 * @brief Get minimum value with location within a grid.
 */
template <typename Type>
__device__ __forceinline__ Type GET_MINLOC_GRID
(Type val, volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  Type ret = GET_MINLOC_BLCK(val, smem, tidx, head);


  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 )
      stloc((Type *)gmem, bidx, ret);

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** calculate global prefix sum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem, 0, ret);
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx == 0 )
      stloc((Type *)smem, 0, ldloc(gmem[0]));
    __syncthreads();

    /** upload obtained result */
    ret = ldloc(smem[0]);
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (ret);
}


/**
 * @fn GET_MAXLOC_GRID
 *
 * @brief Get maximum value with location within a grid.
 */
template <typename Type>
__device__ __forceinline__ Type GET_MAXLOC_GRID
(Type val, volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  Type ret = GET_MAXLOC_BLCK(val, smem, tidx, head);


  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 )
      stloc((Type *)gmem, bidx, ret);

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** calculate global prefix sum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem, 0, ret);
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx == 0 )
      stloc((Type *)smem, 0, ldloc(gmem[0]));
    __syncthreads();

    /** upload obtained result */
    ret = ldloc(smem[0]);
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (ret);
}


/**
 * @fn GET_MINLOC_MAXLOC_GRID
 *
 * @brief Get minimum and maximum values with their locations within a grid.
 */
template <typename Type>
__device__ __forceinline__ void GET_MINLOC_MAXLOC_GRID
(Type * __restrict__ minloc, Type * __restrict__ maxloc, volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem_minloc, volatile Type * __restrict__ gmem_maxloc, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  GET_MINLOC_MAXLOC_BLCK(minloc, maxloc, smem, tidx, head);


  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 ){
      stloc((Type *)gmem_minloc, bidx, *minloc);
      stloc((Type *)gmem_maxloc, bidx, *maxloc);
    }/* if( tidx == 0 ){ */

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get minimum with its location by a representative block */
    if( bidx == 0 ){
      Type ret = *minloc;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_minloc[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_minloc, 0, ret);
    }/* if( bidx == 0 ){ */


    /** get maximum with its location by a representative block */
    if( bidx == (bnum - 1) ){
      Type ret = *maxloc;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_maxloc[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_maxloc, 0, ret);
    }/* if( bidx == (bnum - 1) ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 2 ){
      stloc((Type *)smem, tidx, (tidx == 0) ? ldloc(gmem_minloc[0]) : ldloc(gmem_maxloc[0]));
    }/* if( tidx < 2 ){ */
    __syncthreads();

    /** upload obtained result */
    *minloc = ldloc(smem[0]);
    *maxloc = ldloc(smem[1]);
    __syncthreads();
  }/* if( bnum > 1 ){ */
}


/**
 * @fn GET_MINLOC_2VALS_GRID
 *
 * @brief Get minimum values with their locations within a grid.
 */
template <typename Type>
__device__ __forceinline__ void GET_MINLOC_2VALS_GRID
(Type * __restrict__ xmin, Type * __restrict__ ymin,
 volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem_xmin, volatile Type * __restrict__ gmem_ymin,
 const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  *xmin = GET_MINLOC_BLCK(*xmin, smem, tidx, head);
  *ymin = GET_MINLOC_BLCK(*ymin, smem, tidx, head);

  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 ){
      stloc((Type *)gmem_xmin, bidx, *xmin);
      stloc((Type *)gmem_ymin, bidx, *ymin);
    }/* if( tidx == 0 ){ */

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get minimum with its location by a representative block */
    if( bidx == 0 ){
      Type ret = *xmin;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_xmin[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_xmin, 0, ret);
    }/* if( bidx == 0 ){ */

    /** get minimum with its location by a representative block */
    if( bidx == (1 % bnum) ){
      Type ret = *ymin;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_ymin[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_ymin, 0, ret);
    }/* if( bidx == (1 % bnum) ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 2 ){
      stloc((Type *)smem, tidx, (tidx == 0) ? ldloc(gmem_xmin[0]) : ldloc(gmem_ymin[0]));
    }/* if( tidx < 2 ){ */
    __syncthreads();

    /** upload obtained result */
    *xmin = ldloc(smem[0]);
    *ymin = ldloc(smem[1]);
    __syncthreads();
  }/* if( bnum > 1 ){ */
}


/**
 * @fn GET_MAXLOC_2VALS_GRID
 *
 * @brief Get maximum values with their locations within a grid.
 */
template <typename Type>
__device__ __forceinline__ void GET_MAXLOC_2VALS_GRID
(Type * __restrict__ xmax, Type * __restrict__ ymax,
 volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem_xmax, volatile Type * __restrict__ gmem_ymax,
 const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  *xmax = GET_MAXLOC_BLCK(*xmax, smem, tidx, head);
  *ymax = GET_MAXLOC_BLCK(*ymax, smem, tidx, head);

  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 ){
      stloc((Type *)gmem_xmax, bidx, *xmax);
      stloc((Type *)gmem_ymax, bidx, *ymax);
    }/* if( tidx == 0 ){ */

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get maximum with its location by a representative block */
    if( bidx == 0 ){
      Type ret = *xmax;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_xmax[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_xmax, 0, ret);
    }/* if( bidx == 0 ){ */

    /** get maximum with its location by a representative block */
    if( bidx == (1 % bnum) ){
      Type ret = *ymax;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_ymax[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_ymax, 0, ret);
    }/* if( bidx == (1 % bnum) ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 2 ){
      stloc((Type *)smem, tidx, (tidx == 0) ? ldloc(gmem_xmax[0]) : ldloc(gmem_ymax[0]));
    }/* if( tidx < 2 ){ */
    __syncthreads();

    /** upload obtained result */
    *xmax = ldloc(smem[0]);
    *ymax = ldloc(smem[1]);
    __syncthreads();
  }/* if( bnum > 1 ){ */
}


/**
 * @fn GET_MINLOC_MAXLOC_2VALS_GRID
 *
 * @brief Get minimum and maximum values with their locations within a grid.
 */
template <typename Type>
__device__ __forceinline__ void GET_MINLOC_MAXLOC_2VALS_GRID
(Type * __restrict__ xmin, Type * __restrict__ xmax, Type * __restrict__ ymin, Type * __restrict__ ymax,
 volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem_xmin, volatile Type * __restrict__ gmem_xmax, volatile Type * __restrict__ gmem_ymin, volatile Type * __restrict__ gmem_ymax,
 const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  GET_MINLOC_MAXLOC_BLCK(xmin, xmax, smem, tidx, head);
  GET_MINLOC_MAXLOC_BLCK(ymin, ymax, smem, tidx, head);

  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 ){
      stloc((Type *)gmem_xmin, bidx, *xmin);	  stloc((Type *)gmem_xmax, bidx, *xmax);
      stloc((Type *)gmem_ymin, bidx, *ymin);	  stloc((Type *)gmem_ymax, bidx, *ymax);
    }/* if( tidx == 0 ){ */

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get minimum with its location by a representative block */
    if( bidx == 0 ){
      Type ret = *xmin;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_xmin[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_xmin, 0, ret);
    }/* if( bidx == 0 ){ */

    /** get maximum with its location by a representative block */
    if( bidx == (1 % bnum) ){
      Type ret = *xmax;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_xmax[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_xmax, 0, ret);
    }/* if( bidx == (1 % bnum) ){ */

    /** get minimum with its location by a representative block */
    if( bidx == (2 % bnum) ){
      Type ret = *ymin;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_ymin[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_ymin, 0, ret);
    }/* if( bidx == (2 % bnum) ){ */

    /** get maximum with its location by a representative block */
    if( bidx == (3 % bnum) ){
      Type ret = *ymax;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_ymax[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_ymax, 0, ret);
    }/* if( bidx == (3 % bnum) ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 4 ){
      stloc((Type *)smem, tidx, (tidx < 2) ? ((tidx == 0) ? ldloc(gmem_xmin[0]) : ldloc(gmem_xmax[0])) : ((tidx == 2) ? ldloc(gmem_ymin[0]) : ldloc(gmem_ymax[0])));
    }/* if( tidx < 4 ){ */
    __syncthreads();

    /** upload obtained result */
    *xmin = ldloc(smem[0]);    *xmax = ldloc(smem[1]);
    *ymin = ldloc(smem[2]);    *ymax = ldloc(smem[3]);
    __syncthreads();
  }/* if( bnum > 1 ){ */
}


/**
 * @fn GET_MINLOC_3VALS_GRID
 *
 * @brief Get minimum values with their locations within a grid.
 */
template <typename Type>
__device__ __forceinline__ void GET_MINLOC_3VALS_GRID
(Type * __restrict__ xmin, Type * __restrict__ ymin, Type * __restrict__ zmin,
 volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem_xmin, volatile Type * __restrict__ gmem_ymin, volatile Type * __restrict__ gmem_zmin,
 const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  *xmin = GET_MINLOC_MAXLOC_BLCK(*xmin, smem, tidx, head);
  *ymin = GET_MINLOC_MAXLOC_BLCK(*ymin, smem, tidx, head);
  *zmin = GET_MINLOC_MAXLOC_BLCK(*zmin, smem, tidx, head);

  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 ){
      stloc((Type *)gmem_xmin, bidx, *xmin);
      stloc((Type *)gmem_ymin, bidx, *ymin);
      stloc((Type *)gmem_zmin, bidx, *zmin);
    }/* if( tidx == 0 ){ */

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get minimum with its location by a representative block */
    if( bidx == 0 ){
      Type ret = *xmin;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_xmin[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_xmin, 0, ret);
    }/* if( bidx == 0 ){ */

    /** get minimum with its location by a representative block */
    if( bidx == (1 % bnum) ){
      Type ret = *ymin;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_ymin[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_ymin, 0, ret);
    }/* if( bidx == (1 % bnum) ){ */

    /** get minimum with its location by a representative block */
    if( bidx == (2 % bnum) ){
      Type ret = *zmin;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_zmin[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_zmin, 0, ret);
    }/* if( bidx == (2 % bnum) ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 3 ){
      stloc((Type *)smem, tidx, (tidx < 2) ? ((tidx == 0) ? ldloc(gmem_xmin[0]) : ldloc(gmem_ymin[0])) : ldloc(gmem_zmin[0]));
    }/* if( tidx < 3 ){ */
    __syncthreads();

    /** upload obtained result */
    *xmin = ldloc(smem[0]);
    *ymin = ldloc(smem[1]);
    *zmin = ldloc(smem[2]);
    __syncthreads();
  }/* if( bnum > 1 ){ */
}


/**
 * @fn GET_MAXLOC_3VALS_GRID
 *
 * @brief Get maximum values with their locations within a grid.
 */
template <typename Type>
__device__ __forceinline__ void GET_MAXLOC_3VALS_GRID
(Type * __restrict__ xmax, Type * __restrict__ ymax, Type * __restrict__ zmax,
 volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem_xmax, volatile Type * __restrict__ gmem_ymax, volatile Type * __restrict__ gmem_zmax,
 const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  *xmax = GET_MAXLOC_BLCK(*xmax, smem, tidx, head);
  *ymax = GET_MAXLOC_BLCK(*ymax, smem, tidx, head);
  *zmax = GET_MAXLOC_BLCK(*zmax, smem, tidx, head);

  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 ){
      stloc((Type *)gmem_xmax, bidx, *xmax);
      stloc((Type *)gmem_ymax, bidx, *ymax);
      stloc((Type *)gmem_zmax, bidx, *zmax);
    }/* if( tidx == 0 ){ */

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get maximum with its location by a representative block */
    if( bidx == 0 ){
      Type ret = *xmax;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_xmax[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_xmax, 0, ret);
    }/* if( bidx == 0 ){ */

    /** get maximum with its location by a representative block */
    if( bidx == (1 % bnum) ){
      Type ret = *ymax;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_ymax[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_ymax, 0, ret);
    }/* if( bidx == (1 % bnum) ){ */

    /** get maximum with its location by a representative block */
    if( bidx == (2 % bnum) ){
      Type ret = *zmax;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_zmax[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_zmax, 0, ret);
    }/* if( bidx == (2 % bnum) ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 3 ){
      stloc((Type *)smem, tidx, (tidx < 2) ? ((tidx == 0) ? ldloc(gmem_xmax[0]) : ldloc(gmem_ymax[0])) : ldloc(gmem_zmax[0]));
    }/* if( tidx < 3 ){ */
    __syncthreads();

    /** upload obtained result */
    *xmax = ldloc(smem[0]);
    *ymax = ldloc(smem[1]);
    *zmax = ldloc(smem[2]);
    __syncthreads();
  }/* if( bnum > 1 ){ */
}


/**
 * @fn GET_MINLOC_MAXLOC_3VALS_GRID
 *
 * @brief Get minimum and maximum values with their locations within a grid.
 */
template <typename Type>
__device__ __forceinline__ void GET_MINLOC_MAXLOC_3VALS_GRID
(Type * __restrict__ xmin, Type * __restrict__ xmax, Type * __restrict__ ymin, Type * __restrict__ ymax, Type * __restrict__ zmin, Type * __restrict__ zmax,
 volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem_xmin, volatile Type * __restrict__ gmem_xmax, volatile Type * __restrict__ gmem_ymin, volatile Type * __restrict__ gmem_ymax, volatile Type * __restrict__ gmem_zmin, volatile Type * __restrict__ gmem_zmax,
 const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  GET_MINLOC_MAXLOC_BLCK(xmin, xmax, smem, tidx, head);
  GET_MINLOC_MAXLOC_BLCK(ymin, ymax, smem, tidx, head);
  GET_MINLOC_MAXLOC_BLCK(zmin, zmax, smem, tidx, head);

  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local reduction via global memory */
    /** store data on the global memory */
    if( tidx == 0 ){
      stloc((Type *)gmem_xmin, bidx, *xmin);	   stloc((Type *)gmem_xmax, bidx, *xmax);
      stloc((Type *)gmem_ymin, bidx, *ymin);	   stloc((Type *)gmem_ymax, bidx, *ymax);
      stloc((Type *)gmem_zmin, bidx, *zmin);	   stloc((Type *)gmem_zmax, bidx, *zmax);
    }/* if( tidx == 0 ){ */

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get minimum with its location by a representative block */
    if( bidx == 0 ){
      Type ret = *xmin;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_xmin[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_xmin, 0, ret);
    }/* if( bidx == 0 ){ */

    /** get maximum with its location by a representative block */
    if( bidx == (1 % bnum) ){
      Type ret = *xmax;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_xmax[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_xmax, 0, ret);
    }/* if( bidx == (1 % bnum) ){ */

    /** get minimum with its location by a representative block */
    if( bidx == (2 % bnum) ){
      Type ret = *ymin;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_ymin[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_ymin, 0, ret);
    }/* if( bidx == (2 % bnum) ){ */

    /** get maximum with its location by a representative block */
    if( bidx == (3 % bnum) ){
      Type ret = *ymax;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_ymax[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_ymax, 0, ret);
    }/* if( bidx == (3 % bnum) ){ */

    /** get minimum with its location by a representative block */
    if( bidx == (4 % bnum) ){
      Type ret = *zmin;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_zmin[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MINLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val < ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_zmin, 0, ret);
    }/* if( bidx == (4 % bnum) ){ */

    /** get maximum with its location by a representative block */
    if( bidx == (5 % bnum) ){
      Type ret = *zmax;
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_COMPARE_INC);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_COMPARE_INC;

	/** load from the global memory */
	Type tmp = ((target < bnum) ? ldloc(gmem_zmax[target]) : ret);

	/** calculate local reduction */
	tmp = GET_MAXLOC_BLCK(tmp, smem, tidx, head);
	if( tmp.val > ret.val )
	  ret = tmp;

	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      /** share global reduction via global memory */
      if( tidx == 0 )
	stloc((Type *)gmem_zmax, 0, ret);
    }/* if( bidx == (5 % bnum) ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 6 ){
      stloc((Type *)smem, tidx, (tidx < 2) ? ((tidx == 0) ? ldloc(gmem_xmin[0]) : ldloc(gmem_xmax[0])) : ((tidx < 4) ? ((tidx == 2) ? ldloc(gmem_ymin[0]) : ldloc(gmem_ymax[0])) : ((tidx == 4) ? ldloc(gmem_zmin[0]) : ldloc(gmem_zmax[0]))));
    }/* if( tidx < 6 ){ */
    __syncthreads();

    /** upload obtained result */
    *xmin = ldloc(smem[0]);    *xmax = ldloc(smem[1]);
    *ymin = ldloc(smem[2]);    *ymax = ldloc(smem[3]);
    *zmin = ldloc(smem[4]);    *zmax = ldloc(smem[5]);
    __syncthreads();
  }/* if( bnum > 1 ){ */
}
