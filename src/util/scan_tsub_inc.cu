/**
 * @file scan_tsub_inc.cu
 *
 * @brief Source code for parallel prefix sum library on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/11/16 (Mon)
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

#include "../util/scan_tsub_inc.cuh"

#   if  (GPUGEN >= 70) && !defined(_COOPERATIVE_GROUPS_H_)
#include <cooperative_groups.h>
using namespace cooperative_groups;
#endif//(GPUGEN >= 70) && !defined(_COOPERATIVE_GROUPS_H_)


/**
 * @fn PREFIX_SUM_TSUB
 *
 * @brief Get parallel (inclusive) prefix sum within a group of TSUB_SCAN_INC threads.
 * @detail implicit synchronization within TSUB_SCAN_INC (<= 32) threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type PREFIX_SUM_TSUB
(Type val, const int lane
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
 , const uint mask
#else///USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
 , volatile Type * smem, const int tidx
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
 )
{
#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_SCAN_INC < 32)
  thread_block_tile<TSUB_SCAN_INC> tile = tiled_partition<TSUB_SCAN_INC>(this_thread_block());
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_SCAN_INC < 32)

#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC

  Type tmp;
#   if  TSUB_SCAN_INC >=  2
  tmp = __SHFL_UP(mask, val,  1, TSUB_SCAN_INC);  if( lane >=  1 )    val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC == 32
  __syncwarp();
#else///TSUB_SCAN_INC == 32
  tile.sync();
#endif//TSUB_SCAN_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC >=  4
  tmp = __SHFL_UP(mask, val,  2, TSUB_SCAN_INC);  if( lane >=  2 )    val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC == 32
  __syncwarp();
#else///TSUB_SCAN_INC == 32
  tile.sync();
#endif//TSUB_SCAN_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC >=  8
  tmp = __SHFL_UP(mask, val,  4, TSUB_SCAN_INC);  if( lane >=  4 )    val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC == 32
  __syncwarp();
#else///TSUB_SCAN_INC == 32
  tile.sync();
#endif//TSUB_SCAN_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC >= 16
  tmp = __SHFL_UP(mask, val,  8, TSUB_SCAN_INC);  if( lane >=  8 )    val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC == 32
  __syncwarp();
#else///TSUB_SCAN_INC == 32
  tile.sync();
#endif//TSUB_SCAN_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC == 32
  tmp = __SHFL_UP(mask, val, 16, TSUB_SCAN_INC);  if( lane >= 16 )    val += tmp;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#endif//TSUB_SCAN_INC == 32
#endif//TSUB_SCAN_INC >= 16
#endif//TSUB_SCAN_INC >=  8
#endif//TSUB_SCAN_INC >=  4
#endif//TSUB_SCAN_INC >=  2

#else///USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC

  smem[tidx] = val;
#   if  TSUB_SCAN_INC >=  2
  if( lane >=  1 )    val += smem[tidx -  1];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC == 32
  __syncwarp();
#else///TSUB_SCAN_INC == 32
  tile.sync();
#endif//TSUB_SCAN_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_SCAN_INC >=  4
  if( lane >=  2 )    val += smem[tidx -  2];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC == 32
  __syncwarp();
#else///TSUB_SCAN_INC == 32
  tile.sync();
#endif//TSUB_SCAN_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_SCAN_INC >=  8
  if( lane >=  4 )    val += smem[tidx -  4];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC == 32
  __syncwarp();
#else///TSUB_SCAN_INC == 32
  tile.sync();
#endif//TSUB_SCAN_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_SCAN_INC >= 16
  if( lane >=  8 )    val += smem[tidx -  8];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_INC == 32
  __syncwarp();
#else///TSUB_SCAN_INC == 32
  tile.sync();
#endif//TSUB_SCAN_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#   if  TSUB_SCAN_INC == 32
  if( lane >= 16 )    val += smem[tidx - 16];
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  smem[tidx] = val;
#endif//TSUB_SCAN_INC == 32
#endif//TSUB_SCAN_INC >= 16
#endif//TSUB_SCAN_INC >=  8
#endif//TSUB_SCAN_INC >=  4
#endif//TSUB_SCAN_INC >=  2

#endif//USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC

  return (val);
}


/**
 * @fn TOTAL_SUM_TSUB
 *
 * @brief Get total sum within a group of TSUB_SCAN_INC threads.
 * @detail implicit synchronization within TSUB_SCAN_INC (<= 32) threads (a warp) is assumed.
 */
#ifdef  USE_WARP_REDUCE_FUNC_SCAN_TSUB_INC
__device__ __forceinline__      int TOTAL_SUM_TSUB(     int val, const uint mask){  return (__reduce_add_sync(mask, val));}
__device__ __forceinline__ unsigned TOTAL_SUM_TSUB(unsigned val, const uint mask){  return (__reduce_add_sync(mask, val));}
#endif//USE_WARP_REDUCE_FUNC_SCAN_TSUB_INC
template <typename Type>
__device__ __forceinline__ Type TOTAL_SUM_TSUB
(Type val
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
 , const uint mask
#else///USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
 , volatile Type * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC

  Type tmp;
#   if  TSUB_SCAN_INC >=  2
  tmp = __SHFL_XOR(mask, val,  1, TSUB_SCAN_INC);  val += tmp;
#   if  TSUB_SCAN_INC >=  4
  tmp = __SHFL_XOR(mask, val,  2, TSUB_SCAN_INC);  val += tmp;
#   if  TSUB_SCAN_INC >=  8
  tmp = __SHFL_XOR(mask, val,  4, TSUB_SCAN_INC);  val += tmp;
#   if  TSUB_SCAN_INC >= 16
  tmp = __SHFL_XOR(mask, val,  8, TSUB_SCAN_INC);  val += tmp;
#   if  TSUB_SCAN_INC == 32
  tmp = __SHFL_XOR(mask, val, 16, TSUB_SCAN_INC);  val += tmp;
#endif//TSUB_SCAN_INC == 32
#endif//TSUB_SCAN_INC >= 16
#endif//TSUB_SCAN_INC >=  8
#endif//TSUB_SCAN_INC >=  4
#endif//TSUB_SCAN_INC >=  2
  val = __SHFL(mask, val, 0, TSUB_SCAN_INC);

#else///USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC

  smem[tidx] = val;
#   if  TSUB_SCAN_INC >=  2
  val += smem[tidx ^  1];  smem[tidx] = val;
#   if  TSUB_SCAN_INC >=  4
  val += smem[tidx ^  2];  smem[tidx] = val;
#   if  TSUB_SCAN_INC >=  8
  val += smem[tidx ^  4];  smem[tidx] = val;
#   if  TSUB_SCAN_INC >= 16
  val += smem[tidx ^  8];  smem[tidx] = val;
#   if  TSUB_SCAN_INC == 32
  val += smem[tidx ^ 16];  smem[tidx] = val;
#endif//TSUB_SCAN_INC == 32
#endif//TSUB_SCAN_INC >= 16
#endif//TSUB_SCAN_INC >=  8
#endif//TSUB_SCAN_INC >=  4
#endif//TSUB_SCAN_INC >=  2
  val = smem[head];

#endif//USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC

  return (val);
}
