/**
 * @file buf_inc.cu
 *
 * @brief Source code for booking buffer on global memory
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/01 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <helper_cuda.h>

#include "macro.h"
#include "cudalib.h"

#include "../misc/device.h"

#include "make.h"
#include "buf_inc.h"
#include "walk_dev.h"


#ifdef  USE_SMID_TO_GET_BUFID
/**
 * @fn getSMidx
 *
 * @brief Get index of the current SM (streaming multiprocessor).
 * @detail implementation is taken from https://gist.github.com/allanmac/4751080
 *
 * @return index of the current SM
 */
__device__ __forceinline__ uint getSMidx()
{
  uint rr;
  asm("mov.u32 %0, %%smid;" : "=r"(rr));
  return (rr);
}


/**
 * @fn occupyBuffer
 *
 * @brief Book a free buffer.
 *
 * @param (tidx) thread index within a block
 * @param (freeLst) list of index for the buffer on global memory
 * @return (bufIdx) head index of buffer allocated on shared memory
 * @return index of buffer in freeLst
 *
 * @sa getSMidx
 */
__device__ __forceinline__ int occupyBuffer(const int tidx, uint * RESTRICT freeLst, uint * RESTRICT bufIdx)
{
  int target = 0;

  if( tidx == 0 ){
#   if  NBLOCKS_PER_SM == 2
    target         = (getSMidx() * NBLOCKS_PER_SM) ^ 1;
#else///NBLOCKS_PER_SM == 2
    const int head =  getSMidx() * NBLOCKS_PER_SM;
#endif//NBLOCKS_PER_SM == 2

    uint tmp = UINT_MAX;
    while( tmp == UINT_MAX ){
#   if  NBLOCKS_PER_SM == 2
      target ^= 1;
      tmp = atomicExch(&freeLst[target], UINT_MAX);
#else///NBLOCKS_PER_SM == 2
      for(int ii = 0; ii < NBLOCKS_PER_SM; ii++){
	tmp = atomicExch(&freeLst[head + ii], UINT_MAX);

	if( tmp != UINT_MAX ){
	  target = head + ii;
	  break;
	}/* if( tmp != UINT_MAX ){ */
      }/* for(int ii = 0; ii < NBLOCKS_PER_SM; ii++){ */
#endif//NBLOCKS_PER_SM == 2
    }/* while( tmp == UINT_MAX ){ */

    bufIdx[0] = tmp;
  }/* if( tidx == 0 ){ */

  __syncthreads();

  return (target);
}

/**
 * @fn releaseBuffer
 *
 * @brief Release the current buffer.
 *
 * @param (tidx) thread index within a block
 * @return (freeLst) list of index for the buffer on global memory
 * @param (bufIdx) head index of buffer allocated on shared memory
 * @param (target) index of buffer in freeLst
 */
__device__ __forceinline__ void releaseBuffer(const int tidx, uint *freeLst, uint bufIdx, const int target)
{
  __syncthreads();

  if( tidx == 0 )
    atomicExch(&freeLst[target], bufIdx);
}
#else///USE_SMID_TO_GET_BUFID
#ifdef  TRY_MODE_ABOUT_BUFFER
/**
 * @fn occupyBuffer
 *
 * @brief Book a free buffer.
 *
 * @param (tidx) thread index within a block
 * @param (bidx) index of the thread-block
 * @param (bufNum) number of buffers
 * @param (freeLst) list of index for the buffer on global memory
 * @return (bufIdx) head index of buffer allocated on shared memory
 * @return index of buffer in freeLst
 */
__device__ __forceinline__ int occupyBuffer(const int tidx, const int bidx, const int bufNum, uint * RESTRICT freeLst, uint * RESTRICT bufIdx)
{
  int target = 0;

  if( tidx == 0 ){
    target = bidx % bufNum;

    while( true ){
      const uint tmp = atomicExch(&freeLst[target], UINT_MAX);

      if( tmp != UINT_MAX ){
	bufIdx[0] = tmp;
	break;
      }/* if( tmp != UINT_MAX ){ */

      target++;
      target %= bufNum;
    }/* while( true ){ */
  }/* if( tidx == 0 ){ */

  __syncthreads();

  return (target);
}

/**
 * @fn releaseBuffer
 *
 * @brief Release the current buffer.
 *
 * @param (tidx) thread index within a block
 * @return (freeLst) list of index for the buffer on global memory
 * @param (bufIdx) head index of buffer allocated on shared memory
 * @param (target) index of buffer in freeLst
 */
__device__ __forceinline__ void releaseBuffer(const int tidx, uint *freeLst, uint bufIdx, const int target)
{
  __syncthreads();

  if( tidx == 0 )
    atomicExch(&freeLst[target], bufIdx);
}
#else///TRY_MODE_ABOUT_BUFFER

/**
 * @fn occupyBuffer
 *
 * @brief Book a free buffer.
 *
 * @param (tidx) thread index within a block
 * @param (freeNum) number of free buffers on global memory
 * @param (freeLst) list of index for the buffer on global memory
 * @return (bufIdx) head index of buffer allocated on shared memory
 * @param (active) a shared pointer to lock the shared array
 */
__device__ __forceinline__ void  occupyBuffer(const int tidx, uint * RESTRICT freeNum, uint * RESTRICT freeLst, uint * RESTRICT bufIdx, int * RESTRICT active)
{
  if( tidx == 0 ){
    /** lock the shared array */
    while( true )
      if( atomicAnd(active, 0) )
	break;

    /** pick up a free buffer */
    const uint target = atomicDec(freeNum, UINT_MAX) - 1;
    bufIdx[0] = atomicExch(&freeLst[target], UINT_MAX);    /**< bufIdx[0] = freeLst[target]; freeLst[target] = UINT_MAX; */

    /** release the shared array */
    atomicOr(active, 1);
  }/* if( tidx == 0 ){ */

  __syncthreads();
}

/**
 * @fn releaseBuffer
 *
 * @brief Release the current buffer.
 *
 * @param (tidx) thread index within a block
 * @return (freeNum) number of free buffers on global memory
 * @return (freeLst) list of index for the buffer on global memory
 * @param (bufIdx) head index of buffer allocated on shared memory
 * @param (active) a shared pointer to lock the shared array
 */
__device__ __forceinline__ void releaseBuffer(const int tidx, uint * RESTRICT freeNum, uint * RESTRICT freeLst, const int bufIdx, int * RESTRICT active)
{
  __syncthreads();

  if( tidx == 0 ){
    /** lock the shared array */
    while( true )
      if( atomicAnd(active, 0) )
	break;

    /** release the buffer */
    const uint target = atomicInc(freeNum, UINT_MAX);
    freeLst[target] = (uint)bufIdx;

    /** release the shared array */
    atomicOr(active, 1);
  }/* if( tidx == 0 ){ */
}
#endif//TRY_MODE_ABOUT_BUFFER
#endif//USE_SMID_TO_GET_BUFID
