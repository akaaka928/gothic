/**
 * @file scan_vec4_tsub_nwarp_inc.cu
 *
 * @brief Source code for parallel prefix sum library for 4-components vector on GPU
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

#include "../util/vector_inc.cu"
#include "../util/scan_vec4_tsub_nwarp_inc.cuh"

#   if  (GPUGEN >= 70) && !defined(_COOPERATIVE_GROUPS_H_)
#include <cooperative_groups.h>
using namespace cooperative_groups;
#endif//(GPUGEN >= 70) && !defined(_COOPERATIVE_GROUPS_H_)


/**
 * @fn PREFIX_SUM_VEC4_TSUB_NWARP
 *
 * @brief Get parallel (inclusive) prefix sum within a group of TSUB_TN_SCAN_VEC4_INC threads.
 * @detail implicit synchronization within TSUB_TN_SCAN_VEC4_INC (<= 32) threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type PREFIX_SUM_VEC4_TSUB_NWARP(Type val, const int lane, volatile Type * smem, const int tidx)
{
#   if  (__CUDA_ARCH__ >= 700) && (TSUB_TN_SCAN_VEC4_INC < 32)
  thread_block_tile<TSUB_TN_SCAN_VEC4_INC> tile = tiled_partition<TSUB_TN_SCAN_VEC4_INC>(this_thread_block());
#endif//(__CUDA_ARCH__ >= 700) && (TSUB_TN_SCAN_VEC4_INC < 32)
  Type tmp;

  stvec((Type *)smem, tidx, val);
#   if  TSUB_TN_SCAN_VEC4_INC >=  2
#   if  NWARP_TN_SCAN_VEC4_INC <  2
  if( lane >=  1 ){    tmp = ldvec(smem[tidx -  1]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);  }
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_SCAN_VEC4_INC == 32
  __syncwarp();
#else///TSUB_TN_SCAN_VEC4_INC == 32
  tile.sync();
#endif//TSUB_TN_SCAN_VEC4_INC == 32
#endif//__CUDA_ARCH__ >= 700
#endif//NWARP_TN_SCAN_VEC4_INC <  2
#   if  TSUB_TN_SCAN_VEC4_INC >=  4
#   if  NWARP_TN_SCAN_VEC4_INC <  4
  if( lane >=  2 ){    tmp = ldvec(smem[tidx -  2]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);  }
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_SCAN_VEC4_INC == 32
  __syncwarp();
#else///TSUB_TN_SCAN_VEC4_INC == 32
  tile.sync();
#endif//TSUB_TN_SCAN_VEC4_INC == 32
#endif//__CUDA_ARCH__ >= 700
#endif//NWARP_TN_SCAN_VEC4_INC <  4
#   if  TSUB_TN_SCAN_VEC4_INC >=  8
#   if  NWARP_TN_SCAN_VEC4_INC <  8
  if( lane >=  4 ){    tmp = ldvec(smem[tidx -  4]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);  }
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_SCAN_VEC4_INC == 32
  __syncwarp();
#else///TSUB_TN_SCAN_VEC4_INC == 32
  tile.sync();
#endif//TSUB_TN_SCAN_VEC4_INC == 32
#endif//__CUDA_ARCH__ >= 700
#endif//NWARP_TN_SCAN_VEC4_INC <  8
#   if  TSUB_TN_SCAN_VEC4_INC >= 16
#   if  NWARP_TN_SCAN_VEC4_INC < 16
  if( lane >=  8 ){    tmp = ldvec(smem[tidx -  8]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);  }
#   if  __CUDA_ARCH__ >= 700
#   if  TSUB_TN_SCAN_VEC4_INC == 32
  __syncwarp();
#else///TSUB_TN_SCAN_VEC4_INC == 32
  tile.sync();
#endif//TSUB_TN_SCAN_VEC4_INC == 32
#endif//__CUDA_ARCH__ >= 700
#endif//NWARP_TN_SCAN_VEC4_INC < 16
#   if  TSUB_TN_SCAN_VEC4_INC == 32
#   if  NWARP_TN_SCAN_VEC4_INC < 32
  if( lane >= 16 ){    tmp = ldvec(smem[tidx - 16]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#endif//NWARP_TN_SCAN_VEC4_INC < 32
#endif//TSUB_TN_SCAN_VEC4_INC == 32
#endif//TSUB_TN_SCAN_VEC4_INC >= 16
#endif//TSUB_TN_SCAN_VEC4_INC >=  8
#endif//TSUB_TN_SCAN_VEC4_INC >=  4
#endif//TSUB_TN_SCAN_VEC4_INC >=  2

  return (val);
}


/**
 * @fn TOTAL_SUM_VEC4_TSUB_NWARP
 *
 * @brief Get total sum within a group of TSUB_TN_SCAN_VEC4_INC threads.
 * @detail implicit synchronization within TSUB_TN_SCAN_VEC4_INC (<= 32) threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type TOTAL_SUM_VEC4_TSUB_NWARP(Type val, volatile Type * smem, const int tidx, const int head)
{
  Type tmp;

  stvec((Type *)smem, tidx, val);
#   if  TSUB_TN_SCAN_VEC4_INC >=  2
#   if  NWARP_TN_SCAN_VEC4_INC <  2
  tmp = ldvec(smem[tidx ^  1]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  val.w += tmp.w;  stvec((Type *)smem, tidx, val);
#endif//NWARP_TN_SCAN_VEC4_INC <  2
#   if  TSUB_TN_SCAN_VEC4_INC >=  4
#   if  NWARP_TN_SCAN_VEC4_INC <  4
  tmp = ldvec(smem[tidx ^  2]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  val.w += tmp.w;  stvec((Type *)smem, tidx, val);
#endif//NWARP_TN_SCAN_VEC4_INC <  4
#   if  TSUB_TN_SCAN_VEC4_INC >=  8
#   if  NWARP_TN_SCAN_VEC4_INC <  8
  tmp = ldvec(smem[tidx ^  4]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  val.w += tmp.w;  stvec((Type *)smem, tidx, val);
#endif//NWARP_TN_SCAN_VEC4_INC <  8
#   if  TSUB_TN_SCAN_VEC4_INC >= 16
#   if  NWARP_TN_SCAN_VEC4_INC < 16
  tmp = ldvec(smem[tidx ^  8]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  val.w += tmp.w;  stvec((Type *)smem, tidx, val);
#endif//NWARP_TN_SCAN_VEC4_INC < 16
#   if  TSUB_TN_SCAN_VEC4_INC == 32
#   if  NWARP_TN_SCAN_VEC4_INC < 32
  tmp = ldvec(smem[tidx ^ 16]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  val.w += tmp.w;  stvec((Type *)smem, tidx, val);
#endif//NWARP_TN_SCAN_VEC4_INC < 32
#endif//TSUB_TN_SCAN_VEC4_INC == 32
#endif//TSUB_TN_SCAN_VEC4_INC >= 16
#endif//TSUB_TN_SCAN_VEC4_INC >=  8
#endif//TSUB_TN_SCAN_VEC4_INC >=  4
#endif//TSUB_TN_SCAN_VEC4_INC >=  2

  val = ldvec(smem[head]);

  return (val);
}
