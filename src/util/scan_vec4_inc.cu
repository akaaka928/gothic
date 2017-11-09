/**
 * @file scan_vec4_inc.cu
 *
 * @brief Source code for parallel prefix sum library for 4-components vector on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/11/09 (Thu)
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
#include "../util/scan_vec4_inc.cuh"


#ifndef SCAN_VEC4_INC_CU_MULTI_CALL
#define SCAN_VEC4_INC_CU_MULTI_CALL
/**
 * @fn prefixSumVec4Warp
 *
 * @brief Get parallel (inclusive) prefix sum within a warp.
 * @detail implicit synchronization within 32 threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type prefixSumVec4Warp(Type val, const int lane, volatile Type * smem, const int tidx)
{
  Type tmp;
  stvec((Type *)smem, tidx, val);

  if( lane >=  1 ){    tmp = ldvec(smem[tidx -  1]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);  }
  if( lane >=  2 ){    tmp = ldvec(smem[tidx -  2]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);  }
  if( lane >=  4 ){    tmp = ldvec(smem[tidx -  4]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);  }
  if( lane >=  8 ){    tmp = ldvec(smem[tidx -  8]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);  }
  if( lane >= 16 ){    tmp = ldvec(smem[tidx - 16]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);  }

  return (val);
}


/**
 * @fn totalSumVec4Warp
 *
 * @brief Get total sum within a warp.
 * @detail implicit synchronization within 32 threads (a warp) is assumed.
 */
template <typename Type>
__device__ __forceinline__ Type totalSumVec4Warp(Type val, volatile Type * smem, const int tidx, const int head)
{
  Type tmp;

  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  1]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  val.w += tmp.w;  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  2]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  val.w += tmp.w;  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  4]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  val.w += tmp.w;  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^  8]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  val.w += tmp.w;  stvec((Type *)smem, tidx, val);
  tmp = ldvec(smem[tidx ^ 16]);  val.x += tmp.x;  val.y += tmp.y;  val.z += tmp.z;  val.w += tmp.w;  stvec((Type *)smem, tidx, val);

  val = ldvec(smem[head]);
  return (val);
}
#endif//SCAN_VEC4_INC_CU_MULTI_CALL


/**
 * @fn PREFIX_SUM_VEC4_BLCK
 *
 * @brief Get parallel (inclusive) prefix sum within a block.
 */
template <typename Type>
__device__ __forceinline__ Type PREFIX_SUM_VEC4_BLCK(Type val, volatile Type * __restrict__ smem, const int lane, const int tidx)
{
  /** 1. prefix sum within warp */
  Type scan = prefixSumVec4Warp(val, lane, smem, tidx);


  /** 2. prefix sum about tail of each warp */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_SCAN_VEC4_INC <= 1024 --> NTHREADS_SCAN_VEC4_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_SCAN_VEC4_INC >> 5) ){
    val = ldvec(smem[tidx * warpSize + warpSize - 1]);

    stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_SCAN_VEC4_INC >> 5) >=  2
    Type tmp;
    if( lane >=  1 ){      tmp = ldvec(smem[tidx -  1]);      val.x += tmp.x;      val.y += tmp.y;      val.z += tmp.z;      val.w += tmp.w;      stvec((Type *)smem, tidx, val);    }
#   if  (NTHREADS_SCAN_VEC4_INC >> 5) >=  4
    if( lane >=  2 ){      tmp = ldvec(smem[tidx -  2]);      val.x += tmp.x;      val.y += tmp.y;      val.z += tmp.z;      val.w += tmp.w;      stvec((Type *)smem, tidx, val);    }
#   if  (NTHREADS_SCAN_VEC4_INC >> 5) >=  8
    if( lane >=  4 ){      tmp = ldvec(smem[tidx -  4]);      val.x += tmp.x;      val.y += tmp.y;      val.z += tmp.z;      val.w += tmp.w;      stvec((Type *)smem, tidx, val);    }
#   if  (NTHREADS_SCAN_VEC4_INC >> 5) >= 16
    if( lane >=  8 ){      tmp = ldvec(smem[tidx -  8]);      val.x += tmp.x;      val.y += tmp.y;      val.z += tmp.z;      val.w += tmp.w;      stvec((Type *)smem, tidx, val);    }
#   if  (NTHREADS_SCAN_VEC4_INC >> 5) == 32
    if( lane >= 16 ){      tmp = ldvec(smem[tidx - 16]);      val.x += tmp.x;      val.y += tmp.y;      val.z += tmp.z;      val.w += tmp.w;      stvec((Type *)smem, tidx, val);    }
#endif//(NTHREADS_SCAN_VEC4_INC >> 5) == 32
#endif//(NTHREADS_SCAN_VEC4_INC >> 5) >= 16
#endif//(NTHREADS_SCAN_VEC4_INC >> 5) >=  8
#endif//(NTHREADS_SCAN_VEC4_INC >> 5) >=  4
#endif//(NTHREADS_SCAN_VEC4_INC >> 5) >=  2

  }/* if( tidx < (NTHREADS_SCAN_VEC4_INC >> 5) ){ */
  __syncthreads();


  /** 3. prefix sum within a block */
  /** warpSize = 32 = 2^5 */
  if( tidx >= warpSize ){
    Type tmp = ldvec(smem[(tidx >> 5) - 1]);
    scan.x += tmp.x;    scan.y += tmp.y;    scan.z += tmp.z;    scan.w += tmp.w;
  }/* if( tidx >= warpSize ){ */
  __syncthreads();


  /** 4. upload calculated prefix sum */
  stvec((Type *)smem, tidx, scan);
  __syncthreads();

  return (scan);
}


/**
 * @fn PREFIX_SUM_VEC4_GRID
 *
 * @brief Get parallel (inclusive) prefix sum within a grid.
 *
 * @return (scanNum) total sum of val
 */
template <typename Type>
__device__ __forceinline__ Type PREFIX_SUM_VEC4_GRID
(Type val, volatile Type * __restrict__ smem, const int lane, const int tidx,
 Type * __restrict__ scanNum,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  Type scan = PREFIX_SUM_VEC4_BLCK(val, smem, lane, tidx);
  *scanNum = ldvec(smem[NTHREADS_SCAN_VEC4_INC - 1]);

  /** global prefix sum is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    const Type zero = {0, 0, 0, 0};

    /** share local prefix sum via global memory */
    /** store data on the global memory */
    if( tidx == (NTHREADS_SCAN_VEC4_INC - 1) )
      stvec((Type *)gmem, bidx, scan);

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** calculate global prefix sum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_SCAN_VEC4_INC);
      Type head = zero;
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_SCAN_VEC4_INC;

	/** load from the global memory */
	Type pidx = ((target < bnum) ? ldvec(gmem[target]) : zero);

	/** calculate local prefix sum */
	PREFIX_SUM_VEC4_BLCK(pidx, smem, lane, tidx);

	/** store to the global memory */
	if( target < bnum ){
	  Type tmp = ldvec(smem[tidx]);
	  tmp.x += head.x;	  tmp.y += head.y;	  tmp.z += head.z;	  tmp.w += head.w;
	  stvec((Type *)gmem, target, tmp);
	}/* if( target < bnum ){ */

	head.x += smem[NTHREADS_SCAN_VEC4_INC - 1].x;
	head.y += smem[NTHREADS_SCAN_VEC4_INC - 1].y;
	head.z += smem[NTHREADS_SCAN_VEC4_INC - 1].z;
	head.w += smem[NTHREADS_SCAN_VEC4_INC - 1].w;
	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 2 )
      stvec((Type *)smem, tidx, (tidx == 0) ? ((bidx > 0) ? ldvec(gmem[bidx - 1]) : zero) : ldvec(gmem[bnum - 1]));
    __syncthreads();

    Type psub = ldvec(smem[0]);
    scan.x += psub.x;    scan.y += psub.y;    scan.z += psub.z;    scan.w += psub.w;
    psub = ldvec(smem[1]);
    scanNum->x = psub.x;    scanNum->y = psub.y;    scanNum->z = psub.z;    scanNum->w = psub.w;
    __syncthreads();

    /** upload calculated prefix sum */
    stvec((Type *)smem, tidx, scan);
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (scan);
}


/**
 * @fn PREFIX_SUM_VEC4_GRID_WITH_PARTITION
 *
 * @brief Get parallel (inclusive) prefix sum within a grid.
 *
 * @return (headIdx) head index of the scan in the corresponding block
 * @return (scanNum) total sum of val
 */
template <typename Type>
__device__ __forceinline__ Type PREFIX_SUM_VEC4_GRID_WITH_PARTITION
(Type val, volatile Type * __restrict__ smem, const int lane, const int tidx,
 Type * __restrict__ headIdx, Type * __restrict__ scanNum,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  const Type zero = {0, 0, 0, 0};

  Type scan = PREFIX_SUM_VEC4_BLCK(val, smem, lane, tidx);
  *headIdx = zero;
  *scanNum = ldvec(smem[NTHREADS_SCAN_VEC4_INC - 1]);

  /** global prefix sum is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local prefix sum via global memory */
    /** store data on the global memory */
    if( tidx == (NTHREADS_SCAN_VEC4_INC - 1) )
      stvec((Type *)gmem, bidx, scan);

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** calculate global prefix sum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_SCAN_VEC4_INC);
      Type head = zero;
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_SCAN_VEC4_INC;

	/** load from the global memory */
	Type pidx = ((target < bnum) ? ldvec(gmem[target]) : zero);

	/** calculate local prefix sum */
	PREFIX_SUM_VEC4_BLCK(pidx, smem, lane, tidx);

	/** store to the global memory */
	if( target < bnum ){
	  Type tmp = ldvec(smem[tidx]);
	  tmp.x += head.x;	  tmp.y += head.y;	  tmp.z += head.z;	  tmp.w += head.w;
	  stvec((Type *)gmem, target, tmp);
	}/* if( target < bnum ){ */

	head.x += smem[NTHREADS_SCAN_VEC4_INC - 1].x;
	head.y += smem[NTHREADS_SCAN_VEC4_INC - 1].y;
	head.z += smem[NTHREADS_SCAN_VEC4_INC - 1].z;
	head.w += smem[NTHREADS_SCAN_VEC4_INC - 1].w;
	if( loop != (Nloop - 1) )
	  __syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx < 2 )
      stvec((Type *)smem, tidx, (tidx == 0) ? ((bidx > 0) ? ldvec(gmem[bidx - 1]) : zero) : ldvec(gmem[bnum - 1]));
    __syncthreads();

    Type psub = ldvec(smem[0]);
    *scanNum  = ldvec(smem[1]);
    *headIdx  = psub;
    scan.x += psub.x;    scan.y += psub.y;    scan.z += psub.z;    scan.w += psub.w;
    __syncthreads();

    /** upload calculated prefix sum */
    stvec((Type *)smem, tidx, scan);
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (scan);
}


/**
 * @fn TOTAL_SUM_VEC4_BLCK
 *
 * @brief Get total sum within a block.
 */
template <typename Type>
__device__ __forceinline__ Type TOTAL_SUM_VEC4_BLCK(Type val, volatile Type * __restrict__ smem, const int tidx, const int head)
{
  /** 1. total sum within a warp */
  val = totalSumVec4Warp(val, smem, tidx, head);


  /** 2. reduction of partial sum */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_SCAN_VEC4_INC <= 1024 --> NTHREADS_SCAN_VEC4_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_SCAN_VEC4_INC >> 5) ){
    val = ldvec(smem[tidx * warpSize]);

    stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_SCAN_VEC4_INC >> 5) >=  2
    Type tmp;
    tmp = ldvec(smem[tidx ^  1]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_SCAN_VEC4_INC >> 5) >=  4
    tmp = ldvec(smem[tidx ^  2]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_SCAN_VEC4_INC >> 5) >=  8
    tmp = ldvec(smem[tidx ^  4]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_SCAN_VEC4_INC >> 5) >= 16
    tmp = ldvec(smem[tidx ^  8]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);
#   if  (NTHREADS_SCAN_VEC4_INC >> 5) == 32
    tmp = ldvec(smem[tidx ^ 16]);    val.x += tmp.x;    val.y += tmp.y;    val.z += tmp.z;    val.w += tmp.w;    stvec((Type *)smem, tidx, val);
#endif//(NTHREADS_SCAN_VEC4_INC >> 5) == 32
#endif//(NTHREADS_SCAN_VEC4_INC >> 5) >= 16
#endif//(NTHREADS_SCAN_VEC4_INC >> 5) >=  8
#endif//(NTHREADS_SCAN_VEC4_INC >> 5) >=  4
#endif//(NTHREADS_SCAN_VEC4_INC >> 5) >=  2

  }/* if( tidx < (NTHREADS_SCAN_VEC4_INC >> 5) ){ */
  __syncthreads();
  val = ldvec(smem[0]);
  __syncthreads();

  return (val);
}


/**
 * @fn TOTAL_SUM_VEC4_GRID
 *
 * @brief Get total sum within a grid.
 */
template <typename Type>
__device__ __forceinline__ Type TOTAL_SUM_VEC4_GRID
(Type val, volatile Type * __restrict__ smem, const int tidx, const int head,
 volatile Type * __restrict__ gmem, const int bidx, const int bnum, int * __restrict__ gsync0, int * __restrict__ gsync1)
{
  val = TOTAL_SUM_VEC4_BLCK(val, smem, tidx, head);

  /** global prefix sum is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    const Type zero = {0, 0, 0, 0};

    /** share local prefix sum via global memory */
    /** store data on the global memory */
    if( tidx == 0 )
      stvec((Type *)gmem, bidx, val);

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** calculate global prefix sum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_SCAN_VEC4_INC);
      Type sum = zero;
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_SCAN_VEC4_INC;

	/** load from the global memory */
	Type subset = ((target < bnum) ? ldvec(gmem[target]) : zero);

	/** calculate partial sum */
	subset = TOTAL_SUM_VEC4_BLCK(subset, smem, tidx, head);
	sum.x += subset.x;	sum.y += subset.y;	sum.z += subset.z;	sum.w += subset.w;

	__syncthreads();
      }/* for(int loop = 0; loop < Nloop; loop++){ */
      if( tidx == 0 )
	stvec((Type *)gmem, 0, sum);
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx == 0 )
      stvec((Type *)smem, tidx, ldvec(gmem[0]));
    __syncthreads();

    /** upload calculated total sum */
    val = ldvec(smem[0]);
    __syncthreads();
  }/* if( bnum > 1 ){ */

  return (val);
}
