/**
 * @file compare_tsub_nwarp_inc.cu
 *
 * @brief Source code for comparing values on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/06/01 (Fri)
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


/**
 * @fn GET_MIN_TSUB_NWARP
 *
 * @brief Get minimum value within a group of TSUB_TN_COMPARE_INC threads (NWARP_TN_COMPARE_INC continuous threads have the identical value).
 * @detail implicit synchronization within TSUB_TN_COMPARE_INC (<= 32) threads (a warp) is assumed.
 */
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

#   if  (__CUDA_ARCH__ >= 700) && (TSUB_TN_COMPARE_INC < 32)
  thread_block_tile<TSUB_TN_COMPARE_INC> tile = tiled_partition<TSUB_TN_COMPARE_INC>(this_thread_block());
#endif//(__CUDA_ARCH__ >= 700) && (TSUB_TN_COMPARE_INC < 32)

#   if  TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ (     NWARP_TN_COMPARE_INC)];  if( tmp.val < val.val )    val = tmp;
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 2 * NWARP_TN_COMPARE_INC)];  if( tmp.val < val.val )    val = tmp;
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 4 * NWARP_TN_COMPARE_INC)];  if( tmp.val < val.val )    val = tmp;
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 8 * NWARP_TN_COMPARE_INC)];  if( tmp.val < val.val )    val = tmp;
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ (16 * NWARP_TN_COMPARE_INC)];  if( tmp.val < val.val )    val = tmp;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
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

#   if  (__CUDA_ARCH__ >= 700) && (TSUB_TN_COMPARE_INC < 32)
  thread_block_tile<TSUB_TN_COMPARE_INC> tile = tiled_partition<TSUB_TN_COMPARE_INC>(this_thread_block());
#endif//(__CUDA_ARCH__ >= 700) && (TSUB_TN_COMPARE_INC < 32)

#   if  TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ (     NWARP_TN_COMPARE_INC)];  if( tmp.val > val.val )    val = tmp;
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 2 * NWARP_TN_COMPARE_INC)];  if( tmp.val > val.val )    val = tmp;
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 4 * NWARP_TN_COMPARE_INC)];  if( tmp.val > val.val )    val = tmp;
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ ( 8 * NWARP_TN_COMPARE_INC)];  if( tmp.val > val.val )    val = tmp;
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_COMPARE_INC == 32
  __syncwarp();
#else///TSUB_TN_COMPARE_INC == 32
  tile.sync();
#endif//TSUB_TN_COMPARE_INC == 32
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
  tmp = smem[tidx ^ (16 * NWARP_TN_COMPARE_INC)];  if( tmp.val > val.val )    val = tmp;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#endif//TSUB_TN_COMPARE_INC == (32 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= (16 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 8 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 4 * NWARP_TN_COMPARE_INC)
#endif//TSUB_TN_COMPARE_INC >= ( 2 * NWARP_TN_COMPARE_INC)

  val = smem[head];
  return (val);
}
