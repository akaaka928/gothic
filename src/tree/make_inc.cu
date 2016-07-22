/*************************************************************************\
 *                                                                       *
                  last updated on 2015/11/26(Thu) 18:01:41
 *                                                                       *
 *    Constructing octree structure for collisionless systems            *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <helper_cuda.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../misc/gsync_dev.cu"
//-------------------------------------------------------------------------
#include "../tree/make_inc.h"
//-------------------------------------------------------------------------
#ifndef MAKE_INC_CU_MULTI_CALL
#define MAKE_INC_CU_FIRST_CALL
#endif//MAKE_INC_CU_MULTI_CALL
//-------------------------------------------------------------------------





//-------------------------------------------------------------------------
/* parallel prefix sum within a warp (assuming implicit synchronization) */
/* type of prefix sum is inclusive */
//-------------------------------------------------------------------------
#ifdef  MAKE_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_INC
__device__ __forceinline__  int prefixSumMakeWarp(const int psum,                                      const int lane)
#else///USE_WARP_SHUFFLE_FUNC_MAKE_INC
__device__ __forceinline__ void prefixSumMakeWarp(const int psum, volatile int * smem, const int tidx, const int lane)
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_INC
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_INC
  //-----------------------------------------------------------------------
  int val = psum;
  int tmp;
  /* #   if  warpSize >=  2 */  tmp = __shfl_up(val,  1, warpSize);  if( lane >=  1 )    val += tmp;
  /* #   if  warpSize >=  4 */  tmp = __shfl_up(val,  2, warpSize);  if( lane >=  2 )    val += tmp;
  /* #   if  warpSize >=  8 */  tmp = __shfl_up(val,  4, warpSize);  if( lane >=  4 )    val += tmp;
  /* #   if  warpSize >= 16 */  tmp = __shfl_up(val,  8, warpSize);  if( lane >=  8 )    val += tmp;
  /* #   if  warpSize == 32 */  tmp = __shfl_up(val, 16, warpSize);  if( lane >= 16 )    val += tmp;
  return (val);
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_MAKE_INC
  //-----------------------------------------------------------------------
  smem[tidx] = psum;
  //-----------------------------------------------------------------------
  /* #   if  warpSize >=  2 */  if( lane >=  1 )    smem[tidx] += smem[tidx -  1];
  /* #   if  warpSize >=  4 */  if( lane >=  2 )    smem[tidx] += smem[tidx -  2];
  /* #   if  warpSize >=  8 */  if( lane >=  4 )    smem[tidx] += smem[tidx -  4];
  /* #   if  warpSize >= 16 */  if( lane >=  8 )    smem[tidx] += smem[tidx -  8];
  /* #   if  warpSize == 32 */  if( lane >= 16 )    smem[tidx] += smem[tidx - 16];
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_INC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//MAKE_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__device__ __forceinline__ void PREFIX_SUM_MAKE_INC_BLCK(int val, const int tidx, const int lane, volatile int * RESTRICT psum)
{
  //-----------------------------------------------------------------------
  /* 1. prefix sum within warp */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_INC
  int tmp;
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_INC
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_INC
  //-----------------------------------------------------------------------
  /* calculate inclusive prefix sum */
  /* #   if  warpSize >=  2 */  tmp = __shfl_up(val,  1, warpSize);  if( lane >=  1 )    val += tmp;
  /* #   if  warpSize >=  4 */  tmp = __shfl_up(val,  2, warpSize);  if( lane >=  2 )    val += tmp;
  /* #   if  warpSize >=  8 */  tmp = __shfl_up(val,  4, warpSize);  if( lane >=  4 )    val += tmp;
  /* #   if  warpSize >= 16 */  tmp = __shfl_up(val,  8, warpSize);  if( lane >=  8 )    val += tmp;
  /* #   if  warpSize >= 32 */  tmp = __shfl_up(val, 16, warpSize);  if( lane >= 16 )    val += tmp;
  /* return calculated inclusive prefix sum */
  psum[tidx] = val;
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_MAKE_INC
  //-----------------------------------------------------------------------
  /* load index */
  psum[tidx] = val;
  /* calculate inclusive prefix sum */
  /* #   if  warpSize >=  2 */  if( lane >=  1 )    psum[tidx] += psum[tidx -  1];
  /* #   if  warpSize >=  4 */  if( lane >=  2 )    psum[tidx] += psum[tidx -  2];
  /* #   if  warpSize >=  8 */  if( lane >=  4 )    psum[tidx] += psum[tidx -  4];
  /* #   if  warpSize >= 16 */  if( lane >=  8 )    psum[tidx] += psum[tidx -  8];
  /* #   if  warpSize >= 32 */  if( lane >= 16 )    psum[tidx] += psum[tidx - 16];
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC_MAKE_INC
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 2. prefix sum about tail of each warp */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_INC
  int scan = val;
#else///USE_WARP_SHUFFLE_FUNC_MAKE_INC
  int scan = psum[tidx];
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_INC
  __syncthreads();
  /* warpSize = 32 = 2^5; NTHREADS_MAKE_INC <= 1024 --> NTHREADS_MAKE_INC >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_MAKE_INC >> 5) ){
    //---------------------------------------------------------------------
    val = psum[tidx * warpSize + warpSize - 1];
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_INC
    const int groupSize = NTHREADS_MAKE_INC >> 5;
#   if  (NTHREADS_MAKE_INC >> 5) >=  2
    tmp = __shfl_up(val,  1, groupSize);    if( lane >=  1 )      val += tmp;
#   if  (NTHREADS_MAKE_INC >> 5) >=  4
    tmp = __shfl_up(val,  2, groupSize);    if( lane >=  2 )      val += tmp;
#   if  (NTHREADS_MAKE_INC >> 5) >=  8
    tmp = __shfl_up(val,  4, groupSize);    if( lane >=  4 )      val += tmp;
#   if  (NTHREADS_MAKE_INC >> 5) >= 16
    tmp = __shfl_up(val,  8, groupSize);    if( lane >=  8 )      val += tmp;
#   if  (NTHREADS_MAKE_INC >> 5) == 32
    tmp = __shfl_up(val, 16, groupSize);    if( lane >= 16 )      val += tmp;
#endif//(NTHREADS_MAKE_INC >> 5) == 32
#endif//(NTHREADS_MAKE_INC >> 5) >= 16
#endif//(NTHREADS_MAKE_INC >> 5) >=  8
#endif//(NTHREADS_MAKE_INC >> 5) >=  4
#endif//(NTHREADS_MAKE_INC >> 5) >=  2
    psum[tidx] = val;
#else///USE_WARP_SHUFFLE_FUNC_MAKE_INC
    psum[tidx] = val;
#   if  (NTHREADS_MAKE_INC >> 5) >=  2
    if( lane >=  1 )      psum[tidx] += psum[tidx -  1];
#   if  (NTHREADS_MAKE_INC >> 5) >=  4
    if( lane >=  2 )      psum[tidx] += psum[tidx -  2];
#   if  (NTHREADS_MAKE_INC >> 5) >=  8
    if( lane >=  4 )      psum[tidx] += psum[tidx -  4];
#   if  (NTHREADS_MAKE_INC >> 5) >= 16
    if( lane >=  8 )      psum[tidx] += psum[tidx -  8];
#   if  (NTHREADS_MAKE_INC >> 5) == 32
    if( lane >= 16 )      psum[tidx] += psum[tidx - 16];
#endif//(NTHREADS_MAKE_INC >> 5) == 32
#endif//(NTHREADS_MAKE_INC >> 5) >= 16
#endif//(NTHREADS_MAKE_INC >> 5) >=  8
#endif//(NTHREADS_MAKE_INC >> 5) >=  4
#endif//(NTHREADS_MAKE_INC >> 5) >=  2
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_INC
    //---------------------------------------------------------------------
  }/* if( tidx < (NTHREADS_MAKE_INC >> 5) ){ */
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 3. prefix sum within a block */
  //-----------------------------------------------------------------------
  /* warpSize = 32 = 2^5 */
  if( tidx >= warpSize )
    scan += psum[(tidx >> 5) - 1];
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 4. upload calculate prefix sum */
  //-----------------------------------------------------------------------
  psum[tidx] = scan;
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__device__ __forceinline__ void PREFIX_SUM_MAKE_INC_GRID
(int * RESTRICT headIdx, int * RESTRICT scanNum,
 const int bnum, const int bidx, const int tidx, const int lane,
 volatile int * RESTRICT smem, volatile int * RESTRICT gmem, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  //-----------------------------------------------------------------------
  /* global prefix sum is necessary if # of blocks is greater than unity */
  if( bnum > 1 ){
    //---------------------------------------------------------------------
    /* share local prefix sum via global memory */
    //---------------------------------------------------------------------
    /* store data on the global memory */
    if( tidx == (NTHREADS_MAKE_INC - 1) )
      gmem[bidx] = *scanNum;
      /* gmem[bidx] = smem[tidx]; */
    /* global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* calculate global prefix sum by a representative block */
    //---------------------------------------------------------------------
    if( bidx == 0 ){
      //-------------------------------------------------------------------
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_MAKE_INC);
      int head = 0;
      for(int loop = 0; loop < Nloop; loop++){
	//-----------------------------------------------------------------
#if 1
	//-----------------------------------------------------------------
	const int target = tidx + loop * NTHREADS_MAKE_INC;
	//-----------------------------------------------------------------
	/* load from the global memory */
	int pidx = ((target < bnum) ? gmem[target] : 0);
	/* calculate local prefix sum */
	PREFIX_SUM_MAKE_INC_BLCK(pidx, tidx, lane, smem);
	/* store to the global memory */
	if( target < bnum )
	  gmem[target] = head + smem[tidx];
	//-----------------------------------------------------------------
	head += smem[NTHREADS_MAKE_INC - 1];
	if( loop != (Nloop - 1) )
	  __syncthreads();
	//-----------------------------------------------------------------
#else
	//-----------------------------------------------------------------
	/* load from the global memory */
	int pidx = head + ((target < bnum) ? gmem[target] : 0);
	/* calculate local prefix sum */
	PREFIX_SUM_MAKE_INC_BLCK(pidx, tidx, lane, smem);
	/* store to the global memory */
	if( target < bnum )
	  gmem[target] = smem[tidx];
	//-----------------------------------------------------------------
	head = smem[NTHREADS_MAKE_INC - 1];
	if( loop != (Nloop - 1) )
	  __syncthreads();
	//-----------------------------------------------------------------
#endif
	//-----------------------------------------------------------------
      }/* for(int loop = 0; loop < Nloop; loop++){ */
      //-------------------------------------------------------------------
    }/* if( bidx == 0 ){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    /* load from the global memory */
    if( tidx < 2 )
      smem[tidx] = (tidx == 0) ? ((bidx > 0) ? gmem[bidx - 1] : 0) : (gmem[bnum - 1]);
    __syncthreads();
    *headIdx = smem[0];
    *scanNum = smem[1];
    /* __syncthreads(); */
    //---------------------------------------------------------------------
  }/* if( bnum > 1 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  MAKE_INC_CU_FIRST_CALL
#undef  MAKE_INC_CU_FIRST_CALL
#define MAKE_INC_CU_MULTI_CALL
#endif//MAKE_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------
