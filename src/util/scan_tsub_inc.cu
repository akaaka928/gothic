/**
 * @file scan_tsub_inc.cu
 *
 * @brief Source code for parallel prefix sum library on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/04/04 (Wed)
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


/**
 * @fn PREFIX_SUM_TSUB
 *
 * @brief Get parallel (inclusive) prefix sum within a group of TSUB_SCAN_INC threads.
 * @detail implicit synchronization within TSUB_SCAN_INC (<= 32) threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type PREFIX_SUM_TSUB
(Type val, const int lane
#ifndef USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
 , volatile Type * smem, const int tidx
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC

  Type tmp;
#   if  TSUB_SCAN_INC >=  2
  tmp = __SHFL_UP(val,  1, TSUB_SCAN_INC);  if( lane >=  1 )    val += tmp;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB_SCAN_INC >=  4
  tmp = __SHFL_UP(val,  2, TSUB_SCAN_INC);  if( lane >=  2 )    val += tmp;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB_SCAN_INC >=  8
  tmp = __SHFL_UP(val,  4, TSUB_SCAN_INC);  if( lane >=  4 )    val += tmp;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB_SCAN_INC >= 16
  tmp = __SHFL_UP(val,  8, TSUB_SCAN_INC);  if( lane >=  8 )    val += tmp;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB_SCAN_INC == 32
  tmp = __SHFL_UP(val, 16, TSUB_SCAN_INC);  if( lane >= 16 )    val += tmp;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#endif//TSUB_SCAN_INC == 32
#endif//TSUB_SCAN_INC >= 16
#endif//TSUB_SCAN_INC >=  8
#endif//TSUB_SCAN_INC >=  4
#endif//TSUB_SCAN_INC >=  2

#else///USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC

  smem[tidx] = val;
#   if  TSUB_SCAN_INC >=  2
  if( lane >=  1 )    val += smem[tidx -  1];
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_SCAN_INC >=  4
  if( lane >=  2 )    val += smem[tidx -  2];
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_SCAN_INC >=  8
  if( lane >=  4 )    val += smem[tidx -  4];
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_SCAN_INC >= 16
  if( lane >=  8 )    val += smem[tidx -  8];
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  smem[tidx] = val;
#   if  TSUB_SCAN_INC == 32
  if( lane >= 16 )    val += smem[tidx - 16];
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
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
template <typename Type>
__device__ __forceinline__ Type TOTAL_SUM_TSUB
(Type val
#ifndef USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
 , volatile Type * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
 )
{
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC

  Type tmp;
#   if  TSUB_SCAN_INC >=  2
  tmp = __SHFL_XOR(val,  1, TSUB_SCAN_INC);  val += tmp;
#   if  TSUB_SCAN_INC >=  4
  tmp = __SHFL_XOR(val,  2, TSUB_SCAN_INC);  val += tmp;
#   if  TSUB_SCAN_INC >=  8
  tmp = __SHFL_XOR(val,  4, TSUB_SCAN_INC);  val += tmp;
#   if  TSUB_SCAN_INC >= 16
  tmp = __SHFL_XOR(val,  8, TSUB_SCAN_INC);  val += tmp;
#   if  TSUB_SCAN_INC == 32
  tmp = __SHFL_XOR(val, 16, TSUB_SCAN_INC);  val += tmp;
#endif//TSUB_SCAN_INC == 32
#endif//TSUB_SCAN_INC >= 16
#endif//TSUB_SCAN_INC >=  8
#endif//TSUB_SCAN_INC >=  4
#endif//TSUB_SCAN_INC >=  2
  val = __SHFL(val, 0, TSUB_SCAN_INC);

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
