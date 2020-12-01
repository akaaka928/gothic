/**
 * @file gsync_dev.cu
 *
 * @brief Utility tool for inter-block GPU synchronization
 *        (based on GPU Lock-Free Synchronization by Xiao & Feng 2009)
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/12/01 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#ifndef GSYNC_DEV_CU
#define GSYNC_DEV_CU


#include <stdio.h>
#include <stdlib.h>
#include <helper_cuda.h>

#include "macro.h"
#include "cudalib.h"


#define USE_ATOMIC_OPS_FOR_GSYNC


/**
 * @fn globalSync
 *
 * @brief Function to execute global synchronization.
 *
 * @param (tidx) thread index within a block
 * @param (bidx) block index
 * @param (bnum) number of blocks
 * @param (gsync0) temporary array on the global memory
 * @param (gsync1) temporary array on the global memory
 */
#ifdef  USE_ATOMIC_OPS_FOR_GSYNC
__device__ __forceinline__ void globalSync(const int tidx, const int bidx, const int bnum, int * gsync0, int * gsync1)
#else///USE_ATOMIC_OPS_FOR_GSYNC
__device__ __forceinline__ void globalSync(const int tidx, const int bidx, const int bnum, volatile int * gsync0, volatile int * gsync1)
#endif//USE_ATOMIC_OPS_FOR_GSYNC
{
/* #ifndef NDEBUG */
/*   if(tidx == 0) */
/*     printf("l%d: bidx = %d: sync %d blocks..\n", __LINE__, bidx, bnum); */
/* #endif//NDEBUG */
  /** Phase 0. tell */
  __syncthreads();
#ifdef  USE_ATOMIC_OPS_FOR_GSYNC
  if(tidx == 0)
    atomicOr(&gsync0[bidx], 1);
#else///USE_ATOMIC_OPS_FOR_GSYNC
  if(tidx == 0)
    gsync0[bidx] = 1;
#endif//USE_ATOMIC_OPS_FOR_GSYNC


  /** Phase 1. watch */
  if(bidx == 0){
    const int tnum = BLOCKDIM_X1D;
#ifdef  USE_ATOMIC_OPS_FOR_GSYNC
    for(int ii = tidx; ii < bnum; ii += tnum)
      while(true)
	if(atomicAnd(&gsync0[ii], 0) == 1)
	  break;
#else///USE_ATOMIC_OPS_FOR_GSYNC
    for(int ii = tidx; ii < bnum; ii += tnum)
      while(true)
      	if(gsync0[ii] == 1){
      	  gsync0[ii] = 0;
      	  break;
      	}
#endif//USE_ATOMIC_OPS_FOR_GSYNC

    __syncthreads();

#ifdef  USE_ATOMIC_OPS_FOR_GSYNC
    for(int ii = tidx; ii < bnum; ii += tnum)
      atomicOr(&gsync1[ii], 1);
#else///USE_ATOMIC_OPS_FOR_GSYNC
    for(int ii = tidx; ii < bnum; ii += tnum)
      gsync1[ii] = 1;
#endif//USE_ATOMIC_OPS_FOR_GSYNC
  }/* if( bidx == 0 ){ */


  /** Phase 2. check */
#ifdef  USE_ATOMIC_OPS_FOR_GSYNC
  if(tidx == 0)
    while(true)
      if(atomicAnd(&gsync1[bidx], 0) == 1)
  	break;
#else///USE_ATOMIC_OPS_FOR_GSYNC
  if(tidx == 0)
    while(true)
      if(gsync1[bidx] == 1){
  	gsync1[bidx] = 0;
  	break;
      }/* if( gsync1[bidx] ){ */
#endif//USE_ATOMIC_OPS_FOR_GSYNC

/* #ifndef NDEBUG */
/*   if(tidx == 0) */
/*     printf("l%d: bidx = %d: %d blocks synchronized\n", __LINE__, bidx, bnum); */
/* #endif//NDEBUG */
  __syncthreads();
}


/**
 * @fn initGsync_kernel
 *
 * @brief Initialize arrays for inter-block GPU synchronization.
 *
 * @param (num) number of blocks
 * @return (gsync0) temporary array on the global memory
 * @return (gsync1) temporary array on the global memory
 */
__global__ static void initGsync_kernel(const int num, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  const int gidx = GLOBALIDX_X1D;

  if( gidx < num ){
    gsync0[gidx] = 0;
    gsync1[gidx] = 0;
  }/* if( gidx < num ){ */
}


#endif//GSYNC_DEV_CU
