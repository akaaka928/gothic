/**
 * @file compare_tsub_nwarp_inc.cu
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

#include "../util/compare_tsub_nwarp_inc.cuh"
#include "../util/comparison_inc.cu"

#   if  (GPUGEN >= 70) && !defined(_COOPERATIVE_GROUPS_H_)
#include <cooperative_groups.h>
using namespace cooperative_groups;
#endif//(GPUGEN >= 70) && !defined(_COOPERATIVE_GROUPS_H_)


#ifdef  USE_WARP_REDUCE_FUNCTIONS_COMPARE_TSUB_NWARP_INC
__device__ __forceinline__ uint flipFP32(const uint src){  uint mask = -int(src >> 31)   | 0x80000000;  return (src ^ mask);}
__device__ __forceinline__ uint undoFP32(const uint src){  uint mask = ((src >> 31) - 1) | 0x80000000;  return (src ^ mask);}
#endif//USE_WARP_REDUCE_FUNCTIONS_COMPARE_TSUB_NWARP_INC


/**
 * @fn GET_MIN_TSUB_NWARP
 *
 * @brief Get minimum value within a group of TSUB_TN_COMPARE_INC threads (NWARP_TN_COMPARE_INC continuous threads have the identical value).
 * @detail implicit synchronization within TSUB_TN_COMPARE_INC (<= 32) threads (a warp) is assumed.
 */
#ifdef  USE_WARP_REDUCE_FUNCTIONS_COMPARE_TSUB_NWARP_INC
__device__ __forceinline__      int GET_MIN_TSUB_NWARP(     int val, const uint mask){  return (__reduce_min_sync(mask, val));}
__device__ __forceinline__ unsigned GET_MIN_TSUB_NWARP(unsigned val, const uint mask){  return (__reduce_min_sync(mask, val));}
__device__ __forceinline__    float GET_MIN_TSUB_NWARP(   float val, const uint mask){
  union {uint u; float f;} tmp;
  tmp.f = val;
  tmp.u = undoFP32(GET_MIN_TSUB_NWARP(flipFP32(tmp.u), mask));
  return (tmp.f);
}
#endif//USE_WARP_REDUCE_FUNCTIONS_COMPARE_TSUB_NWARP_INC
template <typename Type>
__device__ __forceinline__ Type GET_MIN_TSUB_NWARP
(Type val
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
 , const uint mask
#else///USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
 , volatile Type * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC

  Type tmp;
#   if  TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  tmp = __SHFL_XOR(mask, val,      NWARP_TN_COMPARE_INC, TSUB_TN_COMPARE_INC);  val = getMinVal(val, tmp);
#   if  TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
  tmp = __SHFL_XOR(mask, val,  2 * NWARP_TN_COMPARE_INC, TSUB_TN_COMPARE_INC);  val = getMinVal(val, tmp);
#   if  TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
  tmp = __SHFL_XOR(mask, val,  4 * NWARP_TN_COMPARE_INC, TSUB_TN_COMPARE_INC);  val = getMinVal(val, tmp);
#   if  TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
  tmp = __SHFL_XOR(mask, val,  8 * NWARP_TN_COMPARE_INC, TSUB_TN_COMPARE_INC);  val = getMinVal(val, tmp);
#   if  TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
  tmp = __SHFL_XOR(mask, val, 16 * NWARP_TN_COMPARE_INC, TSUB_TN_COMPARE_INC);  val = getMinVal(val, tmp);
#endif//TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  val = __SHFL(mask, val, 0, TSUB_TN_COMPARE_INC);

#else///USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC

  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  val = getMinVal(val, smem[tidx ^ (     NWARP_TN_COMPARE_INC)]);  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
  val = getMinVal(val, smem[tidx ^ ( 2 * NWARP_TN_COMPARE_INC)]);  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
  val = getMinVal(val, smem[tidx ^ ( 4 * NWARP_TN_COMPARE_INC)]);  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
  val = getMinVal(val, smem[tidx ^ ( 8 * NWARP_TN_COMPARE_INC)]);  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
  val = getMinVal(val, smem[tidx ^ (16 * NWARP_TN_COMPARE_INC)]);  smem[tidx] = val;
#endif//TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  val = smem[head];

#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC

  return (val);
}


/**
 * @fn GET_MAX_TSUB_NWARP
 *
 * @brief Get maximum value within a group of TSUB_TN_COMPARE_INC threads (NWARP_TN_COMPARE_INC continuous threads have the identical value).
 * @detail implicit synchronization within TSUB_TN_COMPARE_INC (<= 32) threads (a warp) is assumed.
 */
#ifdef  USE_WARP_REDUCE_FUNCTIONS_COMPARE_TSUB_NWARP_INC
__device__ __forceinline__      int GET_MAX_TSUB_NWARP(     int val, const uint mask){  return (__reduce_max_sync(mask, val));}
__device__ __forceinline__ unsigned GET_MAX_TSUB_NWARP(unsigned val, const uint mask){  return (__reduce_max_sync(mask, val));}
__device__ __forceinline__    float GET_MAX_TSUB_NWARP(   float val, const uint mask){
  union {uint u; float f;} tmp;
  tmp.f = val;
  tmp.u = undoFP32(GET_MAX_TSUB_NWARP(flipFP32(tmp.u), mask));
  return (tmp.f);
}
#endif//USE_WARP_REDUCE_FUNCTIONS_COMPARE_TSUB_NWARP_INC
template <typename Type>
__device__ __forceinline__ Type GET_MAX_TSUB_NWARP
(Type val
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
 , const uint mask
#else///USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
 , volatile Type * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC

  Type tmp;
#   if  TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  tmp = __SHFL_XOR(mask, val,      NWARP_TN_COMPARE_INC, TSUB_TN_COMPARE_INC);  val = getMaxVal(val, tmp);
#   if  TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
  tmp = __SHFL_XOR(mask, val,  2 * NWARP_TN_COMPARE_INC, TSUB_TN_COMPARE_INC);  val = getMaxVal(val, tmp);
#   if  TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
  tmp = __SHFL_XOR(mask, val,  4 * NWARP_TN_COMPARE_INC, TSUB_TN_COMPARE_INC);  val = getMaxVal(val, tmp);
#   if  TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
  tmp = __SHFL_XOR(mask, val,  8 * NWARP_TN_COMPARE_INC, TSUB_TN_COMPARE_INC);  val = getMaxVal(val, tmp);
#   if  TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
  tmp = __SHFL_XOR(mask, val, 16 * NWARP_TN_COMPARE_INC, TSUB_TN_COMPARE_INC);  val = getMaxVal(val, tmp);
#endif//TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  val = __SHFL(mask, val, 0, TSUB_TN_COMPARE_INC);

#else///USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC

  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  val = getMaxVal(val, smem[tidx ^ (     NWARP_TN_COMPARE_INC)]);  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
  val = getMaxVal(val, smem[tidx ^ ( 2 * NWARP_TN_COMPARE_INC)]);  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
  val = getMaxVal(val, smem[tidx ^ ( 4 * NWARP_TN_COMPARE_INC)]);  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
  val = getMaxVal(val, smem[tidx ^ ( 8 * NWARP_TN_COMPARE_INC)]);  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
  val = getMaxVal(val, smem[tidx ^ (16 * NWARP_TN_COMPARE_INC)]);  smem[tidx] = val;
#endif//TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  val = smem[head];

#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC

  return (val);
}


/**
 * @fn GET_MINLOC_TSUB_NWARP
 *
 * @brief Get minimum value with location within a group of TSUB_TN_COMPARE_INC threads (NWARP_TN_COMPARE_INC continuous threads have the identical value).
 * @detail implicit synchronization within TSUB_TN_COMPARE_INC (<= 32) threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type GET_MINLOC_TSUB_NWARP
(Type val, volatile Type * smem, const int tidx, const int head)
{
  Type tmp;
  smem[tidx] = val;

#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_TN_COMPARE_INC < 32)
  thread_block_tile<TSUB_TN_COMPARE_INC> tile = tiled_partition<TSUB_TN_COMPARE_INC>(this_thread_block());
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_TN_COMPARE_INC < 32)

#   if  TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ (     NWARP_TN_COMPARE_INC)];  if( tmp.val < val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 2 * NWARP_TN_COMPARE_INC)];  if( tmp.val < val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 4 * NWARP_TN_COMPARE_INC)];  if( tmp.val < val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 8 * NWARP_TN_COMPARE_INC)];  if( tmp.val < val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ (16 * NWARP_TN_COMPARE_INC)];  if( tmp.val < val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#endif//TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)

  val = smem[head];
  return (val);
}


/**
 * @fn GET_MAXLOC_TSUB_NWARP
 *
 * @brief Get maximum value with location within a group of TSUB_TN_COMPARE_INC threads (NWARP_TN_COMPARE_INC continuous threads have the identical value).
 * @detail implicit synchronization within TSUB_TN_COMPARE_INC (<= 32) threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type GET_MAXLOC_TSUB_NWARP
(Type val, volatile Type * smem, const int tidx, const int head)
{
  Type tmp;
  smem[tidx] = val;

#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_TN_COMPARE_INC < 32)
  thread_block_tile<TSUB_TN_COMPARE_INC> tile = tiled_partition<TSUB_TN_COMPARE_INC>(this_thread_block());
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_TN_COMPARE_INC < 32)

#   if  TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ (     NWARP_TN_COMPARE_INC)];  if( tmp.val > val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 2 * NWARP_TN_COMPARE_INC)];  if( tmp.val > val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 4 * NWARP_TN_COMPARE_INC)];  if( tmp.val > val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 8 * NWARP_TN_COMPARE_INC)];  if( tmp.val > val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ (16 * NWARP_TN_COMPARE_INC)];  if( tmp.val > val.val )    val = tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#endif//TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)

  val = smem[head];
  return (val);
}
