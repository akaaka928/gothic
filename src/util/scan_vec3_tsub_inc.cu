/**
 * @file scan_vec3_tsub_inc.cu
 *
 * @brief Source code for parallel prefix sum library for 3-components vector on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/09/14 (Mon)
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
#include "../util/scan_vec3_tsub_inc.cuh"

#   if  (GPUGEN >= 70) && !defined(_COOPERATIVE_GROUPS_H_)
#include <cooperative_groups.h>
using namespace cooperative_groups;
#endif//(GPUGEN >= 70) && !defined(_COOPERATIVE_GROUPS_H_)


/**
 * @fn PREFIX_SUM_VEC3_TSUB
 *
 * @brief Get parallel (inclusive) prefix sum within a group of TSUB_SCAN_VEC3_INC threads.
 * @detail implicit synchronization within TSUB_SCAN_VEC3_INC (<= 32) threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type PREFIX_SUM_VEC3_TSUB(Type val, const int lane, volatile Type * smem, const int tidx)
{
#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_SCAN_VEC3_INC < 32)
  thread_block_tile<TSUB_SCAN_VEC3_INC> tile = tiled_partition<TSUB_SCAN_VEC3_INC>(this_thread_block());
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_SCAN_VEC3_INC < 32)
  Type tmp;

  stvec((Type *)smem, tidx, val);
#   if  TSUB_SCAN_VEC3_INC >=  2
  if( lane >=  1 ){    tmp = ldvec(smem[tidx -  1]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    stvec((Type *)smem, tidx, val);  }
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_VEC3_INC == 32
  __syncwarp();
#else///TSUB_SCAN_VEC3_INC == 32
  tile.sync();
#endif//TSUB_SCAN_VEC3_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_VEC3_INC >=  4
  if( lane >=  2 ){    tmp = ldvec(smem[tidx -  2]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    stvec((Type *)smem, tidx, val);  }
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_VEC3_INC == 32
  __syncwarp();
#else///TSUB_SCAN_VEC3_INC == 32
  tile.sync();
#endif//TSUB_SCAN_VEC3_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_VEC3_INC >=  8
  if( lane >=  4 ){    tmp = ldvec(smem[tidx -  4]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    stvec((Type *)smem, tidx, val);  }
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_VEC3_INC == 32
  __syncwarp();
#else///TSUB_SCAN_VEC3_INC == 32
  tile.sync();
#endif//TSUB_SCAN_VEC3_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_VEC3_INC >= 16
  if( lane >=  8 ){    tmp = ldvec(smem[tidx -  8]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    stvec((Type *)smem, tidx, val);  }
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_VEC3_INC == 32
  __syncwarp();
#else///TSUB_SCAN_VEC3_INC == 32
  tile.sync();
#endif//TSUB_SCAN_VEC3_INC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#   if  TSUB_SCAN_VEC3_INC == 32
  if( lane >= 16 ){    tmp = ldvec(smem[tidx - 16]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    stvec((Type *)smem, tidx, val);  }
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#endif//TSUB_SCAN_VEC3_INC == 32
#endif//TSUB_SCAN_VEC3_INC >= 16
#endif//TSUB_SCAN_VEC3_INC >=  8
#endif//TSUB_SCAN_VEC3_INC >=  4
#endif//TSUB_SCAN_VEC3_INC >=  2

  return (val);
}


/**
 * @fn TOTAL_SUM_VEC3_TSUB
 *
 * @brief Get total sum within a group of TSUB_SCAN_VEC3_INC threads.
 * @detail implicit synchronization within TSUB_SCAN_VEC3_INC (<= 32) threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type TOTAL_SUM_VEC3_TSUB(Type val, volatile Type * smem, const int tidx, const int head)
{
  Type tmp;

  stvec((Type *)smem, tidx, val);
#   if  TSUB_SCAN_VEC3_INC >=  2
  tmp = ldvec(smem[tidx ^  1]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  stvec((Type *)smem, tidx, val);
#   if  TSUB_SCAN_VEC3_INC >=  4
  tmp = ldvec(smem[tidx ^  2]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  stvec((Type *)smem, tidx, val);
#   if  TSUB_SCAN_VEC3_INC >=  8
  tmp = ldvec(smem[tidx ^  4]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  stvec((Type *)smem, tidx, val);
#   if  TSUB_SCAN_VEC3_INC >= 16
  tmp = ldvec(smem[tidx ^  8]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  stvec((Type *)smem, tidx, val);
#   if  TSUB_SCAN_VEC3_INC == 32
  tmp = ldvec(smem[tidx ^ 16]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  stvec((Type *)smem, tidx, val);
#endif//TSUB_SCAN_VEC3_INC == 32
#endif//TSUB_SCAN_VEC3_INC >= 16
#endif//TSUB_SCAN_VEC3_INC >=  8
#endif//TSUB_SCAN_VEC3_INC >=  4
#endif//TSUB_SCAN_VEC3_INC >=  2

  val = ldvec(smem[head]);

  return (val);
}
