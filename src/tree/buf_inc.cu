/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:47:01
 *                                                                       *
 *    Octree N-body calculation for collisionless systems on NVIDIA GPUs *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <helper_cuda.h>
//-------------------------------------------------------------------------
#include "macro.h"
#include "cudalib.h"
//-------------------------------------------------------------------------
#include "../misc/device.h"
//-------------------------------------------------------------------------
#include "make.h"
#include "buf_inc.h"
#include "walk_dev.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/*    tidx: thread index within a block */
/* freeNum: number of free buffers on global memory */
/* freeLst: list of index for the buffer on global memory */
/*  bufIdx: head index of buffer allocated on shared memory */
//-------------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID
//-------------------------------------------------------------------------
__device__ __forceinline__ uint getSMidx()
{
  //-----------------------------------------------------------------------
  uint rr;
  asm("mov.u32 %0, %%smid;" : "=r"(rr));
  return (rr);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ int occupyBuffer(const int tidx, uint * RESTRICT freeLst, uint * RESTRICT bufIdx)
{
  //-----------------------------------------------------------------------
  int target = 0;
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
#   if  NBLOCKS_PER_SM == 2
    target         = (getSMidx() * NBLOCKS_PER_SM) ^ 1;
#else///NBLOCKS_PER_SM == 2
    const int head =  getSMidx() * NBLOCKS_PER_SM;
#endif//NBLOCKS_PER_SM == 2
    //---------------------------------------------------------------------
    uint tmp = UINT_MAX;
#if 0
    int counter = 0;
#endif
    while( tmp == UINT_MAX ){
      //-------------------------------------------------------------------
#   if  NBLOCKS_PER_SM == 2
      //-------------------------------------------------------------------
      target ^= 1;
      tmp = atomicExch(&freeLst[target], UINT_MAX);
#if 0
      if( tmp == UINT_MAX )
	printf("wait(%d) on SM %d\n", counter, target >> 1);
      counter++;
#endif
      //-------------------------------------------------------------------
#else///NBLOCKS_PER_SM == 2
      //-------------------------------------------------------------------
      for(int ii = 0; ii < NBLOCKS_PER_SM; ii++){
	//-----------------------------------------------------------------
	tmp = atomicExch(&freeLst[head + ii], UINT_MAX);
	//-----------------------------------------------------------------
	if( tmp != UINT_MAX ){
	  target = head + ii;
	  break;
	}/* if( tmp != UINT_MAX ){ */
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < NBLOCKS_PER_SM; ii++){ */
      //-------------------------------------------------------------------
#endif//NBLOCKS_PER_SM == 2
      //-------------------------------------------------------------------
    }/* while( tmp == UINT_MAX ){ */
    //---------------------------------------------------------------------
    bufIdx[0] = tmp;
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  return (target);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void releaseBuffer(const int tidx, uint *freeLst, uint bufIdx, const int target)
{
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    atomicExch(&freeLst[target], bufIdx);
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///USE_SMID_TO_GET_BUFID
//-------------------------------------------------------------------------
#ifdef  TRY_MODE_ABOUT_BUFFER
//-------------------------------------------------------------------------
__device__ __forceinline__ int occupyBuffer(const int tidx, const int bidx, const int bufNum, uint * RESTRICT freeLst, uint * RESTRICT bufIdx)
{
  //-----------------------------------------------------------------------
  int target = 0;
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    target = bidx % bufNum;
    //---------------------------------------------------------------------
    while( true ){
      //-------------------------------------------------------------------
      const uint tmp = atomicExch(&freeLst[target], UINT_MAX);
      //-------------------------------------------------------------------
      if( tmp != UINT_MAX ){
	bufIdx[0] = tmp;
	break;
      }/* if( tmp != UINT_MAX ){ */
      //-------------------------------------------------------------------
      target++;
      target %= bufNum;
      //-------------------------------------------------------------------
    }/* while( true ){ */
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  return (target);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void releaseBuffer(const int tidx, uint *freeLst, uint bufIdx, const int target)
{
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    atomicExch(&freeLst[target], bufIdx);
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///TRY_MODE_ABOUT_BUFFER
//-------------------------------------------------------------------------
__device__ __forceinline__ void  occupyBuffer(const int tidx, uint * RESTRICT freeNum, uint * RESTRICT freeLst, uint * RESTRICT bufIdx, int * RESTRICT active)
{
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    /* lock the shared array */
    while( true )
      if( atomicAnd(active, 0) )
	break;
    //---------------------------------------------------------------------
    /* pick up a free buffer */
    const uint target = atomicDec(freeNum, UINT_MAX) - 1;
    /* bufIdx [0]      = freeLst[target]; */
    /* freeLst[target] = UINT_MAX; */
    bufIdx[0] = atomicExch(&freeLst[target], UINT_MAX);
#ifdef  DBG_TREE_WALK
    printf("%u buffers are free, index of %u is in use\n", target + 1, bufIdx[0]);
#endif//DBG_TREE_WALK
    //---------------------------------------------------------------------
    /* release the shared array */
    atomicOr(active, 1);
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void releaseBuffer(const int tidx, uint * RESTRICT freeNum, uint * RESTRICT freeLst, const int bufIdx, int * RESTRICT active)
{
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    /* lock the shared array */
    while( true )
      if( atomicAnd(active, 0) )
	break;
    //---------------------------------------------------------------------
    /* release the buffer */
    const uint target = atomicInc(freeNum, UINT_MAX);
    freeLst[target] = (uint)bufIdx;
#ifdef  DBG_TREE_WALK
    printf("%u buffers are free, index of %u is released\n", target + 1, bufIdx);
#endif//DBG_TREE_WALK
    //---------------------------------------------------------------------
    /* release the shared array */
    atomicOr(active, 1);
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//TRY_MODE_ABOUT_BUFFER
//-------------------------------------------------------------------------
#endif//USE_SMID_TO_GET_BUFID
//-------------------------------------------------------------------------
