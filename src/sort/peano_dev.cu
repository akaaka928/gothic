/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/14(Thu) 11:15:50
 *                                                                       *
 *    Key generation of Peano-Hilbert space filling curve                *
 *    sort N-body particles to ovey Peano-Hilbert space filling curve    *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <helper_cuda.h>
//-------------------------------------------------------------------------
/* library load for radix sort */
#ifdef  CUB_AVAILABLE
#include <cub/device/device_radix_sort.cuh>
#else///CUB_AVAILABLE
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#endif//CUB_AVAILABLE
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "peano.h"
#include "peano_dev.h"
//-------------------------------------------------------------------------
#include "../misc/gsync_dev.cu"
//-------------------------------------------------------------------------
#include "../tree/macutil.h"
//-------------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
#include "../tree/make.h"/* <-- required to read NLEAF */
#endif//MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------
/* in L1 cache preferred configuration, capacity of shared memory is 16KiB per SM */
/* real4 smem[NTHREADS_PH] corresponds 16 * NTHREADS_PH bytes */
#define NBLOCKS_PER_SM_PH (1024 / NTHREADS_PH)
//-------------------------------------------------------------------------
#define REGISTERS_PER_THREAD_PH (40)
/* calcPHkey_kernel uses 32 registers @ Tesla M2090, Ttot = 1024 (registers are spilled to local memory) */
/* calcPHkey_kernel uses 47 registers @ Tesla M2090, Ttot =  512 */
/* calcPHkey_kernel uses 36 registers @ Tesla M2090, Ttot =  128, 256 */
#   if  GPUVER == 20
#undef  REGISTERS_PER_THREAD_PH
#          if  NTHREADS_PH == 1024
#define REGISTERS_PER_THREAD_PH (32)
#       else///NTHREADS_PH == 1024
#          if  NTHREADS_PH ==  512
#define REGISTERS_PER_THREAD_PH (47)
#       else///NTHREADS_PH ==  512
#define REGISTERS_PER_THREAD_PH (36)
#       endif//NTHREADS_PH ==  512
#       endif//NTHREADS_PH == 1024
#endif//GPUVER == 20
/* calcPHkey_kernel uses 38 registers @ Tesla K20X */
/* #   if  GPUVER == 35 */
/* #undef  REGISTERS_PER_THREAD_PH */
/* #define REGISTERS_PER_THREAD_PH (38) */
/* #endif//GPUVER == 35 */
/* calcPHkey_kernel uses 40 registers @ GTX 750 Ti */
/* #   if  GPUVER == 50 */
/* #undef  REGISTERS_PER_THREAD_PH */
/* #define REGISTERS_PER_THREAD_PH (40) */
/* #endif//GPUVER == 50 */
/* calcPHkey_kernel uses 40 registers @ GTX 970 */
/* #   if  GPUVER == 52 */
/* #undef  REGISTERS_PER_THREAD_PH */
/* #define REGISTERS_PER_THREAD_PH (40) */
/* #endif//GPUVER == 52 */
//-------------------------------------------------------------------------
/* limitation from number of registers */
#   if  NBLOCKS_PER_SM_PH > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_PH * NTHREADS_PH))
#undef  NBLOCKS_PER_SM_PH
#define NBLOCKS_PER_SM_PH   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_PH * NTHREADS_PH))
#endif//NBLOCKS_PER_SM_PH > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_PH * NTHREADS_PH))
//-------------------------------------------------------------------------
/* maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_PH > (MAX_THREADS_PER_SM / NTHREADS_PH)
#undef  NBLOCKS_PER_SM_PH
#define NBLOCKS_PER_SM_PH   (MAX_THREADS_PER_SM / NTHREADS_PH)
#endif//NBLOCKS_PER_SM_PH > (MAX_THREADS_PER_SM / NTHREADS_PH)
//-------------------------------------------------------------------------
/* maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_PH > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_PH
#define NBLOCKS_PER_SM_PH   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_PH > MAX_BLOCKS_PER_SM
//-------------------------------------------------------------------------
/* maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_PH > ((MAX_WARPS_PER_SM * 32) / NTHREADS_PH)
#undef  NBLOCKS_PER_SM_PH
#define NBLOCKS_PER_SM_PH   ((MAX_WARPS_PER_SM * 32) / NTHREADS_PH)
#endif//NBLOCKS_PER_SM_PH > ((MAX_WARPS_PER_SM * 32) / NTHREADS_PH)
//-------------------------------------------------------------------------
/* # of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_PH < 1
#undef  NBLOCKS_PER_SM_PH
#define NBLOCKS_PER_SM_PH  (1)
#endif//NBLOCKS_PER_SM_PH < 1
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* works less equal than 63 (= 3 * 21) bits key */
/* based on Raman & Wise (2007), IEEE Trans. Comput., C99, 13 */
__device__ __forceinline__ PHint dilate3D(const PHint val)
{
  //-----------------------------------------------------------------------
  PHint ret = val;
  //-----------------------------------------------------------------------
#   if  MAXIMUM_PHKEY_LEVEL > 16
  ret = (ret * 0x100000001) & 0x7fff00000000ffff;/* ((x << 32) + x) = (2^32 + 1) * x = (16^8 + 1) * x */
#endif//MAXIMUM_PHKEY_LEVEL > 16
  ret = (ret * 0x000010001) & 0x00ff0000ff0000ff;/* ((x << 16) + x) = (2^16 + 1) * x = (16^4 + 1) * x */
  ret = (ret * 0x000000101) & 0x700f00f00f00f00f;/* ((x <<  8) + x) = (2^8  + 1) * x = (16^2 + 1) * x */
  ret = (ret * 0x000000011) & 0x30c30c30c30c30c3;/* ((x <<  4) + x) = (2^4  + 1) * x = (16^1 + 1) * x */
  ret = (ret * 0x000000005) & 0x1249249249249249;/* ((x <<  2) + x) = (2^2  + 1) * x = (       5) * x */
  //-----------------------------------------------------------------------
  return (ret);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* 1. estimate box size */
/* 2. calculate Peano--Hilbert key */
//-------------------------------------------------------------------------
__global__ void __launch_bounds__(NTHREADS_PH, NBLOCKS_PER_SM_PH) calcPHkey_kernel
     (const int num, READ_ONLY position * RESTRICT ipos, real4 * RESTRICT min_all, real4 * RESTRICT max_all,
      const int nlevel, int * RESTRICT idx, PHint * RESTRICT key, int * RESTRICT gsync0, int * RESTRICT gsync1
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
      , position * RESTRICT center_dev
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
      )
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int bnum =   GRIDDIM_X1D;
  //-----------------------------------------------------------------------
  const int uidx = tidx + (NTHREADS_PH >> 1);
  const int hidx = tidx - (tidx & (warpSize - 1));
  //-----------------------------------------------------------------------
  const int ihead = (num *      bidx ) / bnum;
  const int itail = (num * (1 + bidx)) / bnum;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate required box size to contain all N-body particles in the local domain */
  //-----------------------------------------------------------------------
  real4 min = { REAL_MAX,  REAL_MAX,  REAL_MAX,  REAL_MAX};
  real4 max = {-REAL_MAX, -REAL_MAX, -REAL_MAX, -REAL_MAX};
  //-----------------------------------------------------------------------
  for(int ih = ihead; ih < itail; ih += NTHREADS_PH){
    //---------------------------------------------------------------------
    const int ii = ih + tidx;
    //---------------------------------------------------------------------
    if( ii < itail ){
      //-------------------------------------------------------------------
      const position pi = ipos[ii];
      //-------------------------------------------------------------------
      if( pi.x < min.x )	min.x = pi.x;
      if( pi.y < min.y )	min.y = pi.y;
      if( pi.z < min.z )	min.z = pi.z;
      //-------------------------------------------------------------------
      if( pi.x > max.x )	max.x = pi.x;
      if( pi.y > max.y )	max.y = pi.y;
      if( pi.z > max.z )	max.z = pi.z;
      //-------------------------------------------------------------------
    }/* if( ii < itail ){ */
    //---------------------------------------------------------------------
  }/* for(int ih = ihead; ih < itail; ih += NTHREADS_PH){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* get the box size within a warp */
  //-----------------------------------------------------------------------
  __shared__  real4 smem[NTHREADS_PH];
  real4 tmp;
  //-----------------------------------------------------------------------
  /* get minimum */
  smem[tidx] = min;
  tmp = smem[tidx ^  1];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
  tmp = smem[tidx ^  2];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
  tmp = smem[tidx ^  4];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
  tmp = smem[tidx ^  8];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
  tmp = smem[tidx ^ 16];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
  /* min = smem[hidx]; */
  //-----------------------------------------------------------------------
  /* get maximum */
  smem[tidx] = max;
  tmp = smem[tidx ^  1];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
  tmp = smem[tidx ^  2];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
  tmp = smem[tidx ^  4];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
  tmp = smem[tidx ^  8];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
  tmp = smem[tidx ^ 16];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
  /* max = smem[hidx]; */
  //-----------------------------------------------------------------------
/* #ifndef NDEBUG */
/*   if( (tidx & (warpSize - 1)) == 0 ){ */
/*     printf("min(%d:%d): %f, %f, %f\n", bidx, tidx / warpSize, min.x, min.y, min.z); */
/*     printf("max(%d:%d): %f, %f, %f\n", bidx, tidx / warpSize, max.x, max.y, max.z); */
/*   } */
/* #endif//NDEBUG */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* get the box size within a block */
  //-----------------------------------------------------------------------
  /* 1024 = 32^2 = warpSize^2 is the maximum of the number of threads; therefore, Nwarp <= warpSize */
  const int Nwarp = NTHREADS_PH >> 5;/* := NTHREADS_PH / warpSize; */
  __syncthreads();
  //-----------------------------------------------------------------------
  if( hidx == tidx ){
    //---------------------------------------------------------------------
    smem[                     (tidx >> 5)] = min;/* := smem[                     (tidx / warpSize)] = min; */
    smem[(NTHREADS_PH >> 1) + (tidx >> 5)] = max;/* := smem[(NTHREADS_PH >> 1) + (tidx / warpSize)] = max; */
    //---------------------------------------------------------------------
  }/* if( hidx == tidx ){ */
  __syncthreads();
  //-----------------------------------------------------------------------
  if( tidx < Nwarp ){
    //---------------------------------------------------------------------
    /* get minimum */
    min = smem[tidx];
#   if  NTHREADS_PH >=   64
    tmp = smem[tidx ^  1];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
#   if  NTHREADS_PH >=  128
    tmp = smem[tidx ^  2];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
#   if  NTHREADS_PH >=  256
    tmp = smem[tidx ^  4];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
#   if  NTHREADS_PH >=  512
    tmp = smem[tidx ^  8];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
#   if  NTHREADS_PH == 1024
    tmp = smem[tidx ^ 16];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
#endif//NTHREADS_PH == 1024
#endif//NTHREADS_PH >=  512
#endif//NTHREADS_PH >=  256
#endif//NTHREADS_PH >=  128
#endif//NTHREADS_PH >=   64
    min = smem[0];
    //---------------------------------------------------------------------
    /* get maximum */
    max = smem[uidx];
    /* smem[tidx] = max; */
#   if  NTHREADS_PH >=   64
    tmp = smem[uidx ^  1];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[uidx] = max;
#   if  NTHREADS_PH >=  128
    tmp = smem[uidx ^  2];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[uidx] = max;
#   if  NTHREADS_PH >=  256
    tmp = smem[uidx ^  4];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[uidx] = max;
#   if  NTHREADS_PH >=  512
    tmp = smem[uidx ^  8];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[uidx] = max;
#   if  NTHREADS_PH == 1024
    tmp = smem[uidx ^ 16];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[uidx] = max;
#endif//NTHREADS_PH == 1024
#endif//NTHREADS_PH >=  512
#endif//NTHREADS_PH >=  256
#endif//NTHREADS_PH >=  128
#endif//NTHREADS_PH >=   64
    max = smem[NTHREADS_PH >> 1];
    //---------------------------------------------------------------------
  }/* if( tidx < Nwarp ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  /* min = smem[               0]; */
  /* max = smem[NTHREADS_PH >> 1]; */
  //-----------------------------------------------------------------------
/* #ifndef NDEBUG */
/*   if( tidx  == 0 ){ */
/*     printf("min(%d): %f, %f, %f\n", bidx, min.x, min.y, min.z); */
/*     printf("max(%d): %f, %f, %f\n", bidx, max.x, max.y, max.z); */
/*   } */
/* #endif//NDEBUG */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* get the box size within a grid */
  //-----------------------------------------------------------------------
  if( bnum > 1 ){
    //---------------------------------------------------------------------
    if( tidx == 0 ){
      //-------------------------------------------------------------------
      min_all[bidx] = min;
      max_all[bidx] = max;
      //-------------------------------------------------------------------
#if 0
      printf("%d-th block:\t[%e:%e], [%e:%e], [%e:%e]\n", bidx, min.x, max.x, min.y, max.y, min.z, max.z);
#endif
      //-------------------------------------------------------------------
    }/* if( tidx == 0 ){ */
    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* get minimum */
    //---------------------------------------------------------------------
    /* bnum is set to (# of SMs) * NBLOCKS_PER_SM_PH <= NTHREADS_PH */
    /* ---> # of working block is 1 (bidx = 0) */
    if( bidx == 0 ){
      //-------------------------------------------------------------------
      smem[tidx] = min;
      //-------------------------------------------------------------------
      if( tidx < bnum  ){
	//-----------------------------------------------------------------
	min = min_all[tidx];
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* reduction within a warp */
	//-----------------------------------------------------------------
	smem[tidx] = min;
	if( bnum >  1 ){
	  tmp = smem[tidx ^  1];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
	if( bnum >  2 ){
	  tmp = smem[tidx ^  2];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
	if( bnum >  4 ){
	  tmp = smem[tidx ^  4];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
	if( bnum >  8 ){
	  tmp = smem[tidx ^  8];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
	if( bnum > 16 ){
	  tmp = smem[tidx ^ 16];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
	}}}}}
	//-----------------------------------------------------------------
#if 0
	if( tidx == hidx )
	  printf("%d-th warp (min):\t% e, % e, % e\n", tidx >> 5, min.x, min.y, min.z);
#endif
	//-----------------------------------------------------------------
      }/* if( tidx < bnum  ){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* reduction within a block (final step) */
      //-------------------------------------------------------------------
      if( bnum > warpSize ){
	//-----------------------------------------------------------------
	const int Ndata = BLOCKSIZE(bnum, warpSize);
	__syncthreads();
	//-----------------------------------------------------------------
	if( (hidx == tidx) && (tidx < bnum) )
	  smem[(tidx >> 5)] = min;/* := (tidx / warpSize) */
	__syncthreads();
	//-----------------------------------------------------------------
	/* bnum <= 1024 --->>> Ndata <= 32 */
	if( tidx < Ndata ){
	  //---------------------------------------------------------------
	  min = smem[tidx];
	  if( Ndata >  1 ){
	    tmp = smem[tidx ^  1];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
	  if( Ndata >  2 ){
	    tmp = smem[tidx ^  2];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
	  if( Ndata >  4 ){
	    tmp = smem[tidx ^  4];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
	  if( Ndata >  8 ){
	    tmp = smem[tidx ^  8];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
	  if( Ndata > 16 ){
	    tmp = smem[tidx ^ 16];  if( tmp.x < min.x )    min.x = tmp.x;  if( tmp.y < min.y )    min.y = tmp.y;  if( tmp.z < min.z )    min.z = tmp.z;  smem[tidx] = min;
	  }}}}}
	  //---------------------------------------------------------------
	}/* if( tidx < Ndata ){ */
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------
      }/* if( bnum > warpSize ){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      if( tidx == 0 )
	min_all[0] = min;
      //-------------------------------------------------------------------
    }/* if( bidx == 0 ){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* get maximum */
    //---------------------------------------------------------------------
    /* bnum is set to (# of SMs) * NBLOCKS_PER_SM_PH <= NTHREADS_PH */
    /* ---> # of working block is 1 (bidx = bnum - 1) */
    if( bidx == (bnum - 1) ){
      //-------------------------------------------------------------------
      smem[tidx] = max;
      //-------------------------------------------------------------------
      if( tidx < bnum  ){
	//-----------------------------------------------------------------
	max = max_all[tidx];
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* reduction within a warp */
	//-----------------------------------------------------------------
	smem[tidx] = max;
	if( bnum >  1 ){
	  tmp = smem[tidx ^  1];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
	if( bnum >  2 ){
	  tmp = smem[tidx ^  2];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
	if( bnum >  4 ){
	  tmp = smem[tidx ^  4];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
	if( bnum >  8 ){
	  tmp = smem[tidx ^  8];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
	if( bnum > 16 ){
	  tmp = smem[tidx ^ 16];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
	}}}}}
	//-----------------------------------------------------------------
#if 0
	if( tidx == hidx )
	  printf("%d-th warp (max):\t% e, % e, % e\n", tidx >> 5, max.x, max.y, max.z);
#endif
	//-----------------------------------------------------------------
      }/* if( tidx < bnum  ){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* reduction within a block (final step) */
      //-------------------------------------------------------------------
      if( bnum > warpSize ){
	//-----------------------------------------------------------------
	const int Ndata = BLOCKSIZE(bnum, warpSize);
	__syncthreads();
	//-----------------------------------------------------------------
	if( (hidx == tidx) && (tidx < bnum) )
	  smem[(tidx >> 5)] = max;/* := (tidx / warpSize) */
	__syncthreads();
	//-----------------------------------------------------------------
	/* bnum <= 1024 --->>> Ndata <= 32 */
	if( tidx < Ndata ){
	  //---------------------------------------------------------------
	  max = smem[tidx];
	  if( Ndata >  1 ){
	    tmp = smem[tidx ^  1];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
	  if( Ndata >  2 ){
	    tmp = smem[tidx ^  2];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
	  if( Ndata >  4 ){
	    tmp = smem[tidx ^  4];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
	  if( Ndata >  8 ){
	    tmp = smem[tidx ^  8];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
	  if( Ndata > 16 ){
	    tmp = smem[tidx ^ 16];  if( tmp.x > max.x )    max.x = tmp.x;  if( tmp.y > max.y )    max.y = tmp.y;  if( tmp.z > max.z )    max.z = tmp.z;  smem[tidx] = max;
	  }}}}}
	  //---------------------------------------------------------------
	}/* if( tidx < Ndata ){ */
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------
      }/* if( bnum > warpSize ){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      if( tidx == 0 )
	max_all[0] = max;
      //-------------------------------------------------------------------
    }/* if( bidx == (bnum - 1) ){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------
    min = min_all[0];
    max = max_all[0];
    //---------------------------------------------------------------------
  }/* if( bnum > 1 ){ */
  else{
    //---------------------------------------------------------------------
    min = smem[               0];
    max = smem[NTHREADS_PH >> 1];
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  //-----------------------------------------------------------------------
  if( (tidx + bidx) == 0 ){
    //---------------------------------------------------------------------
    position center;
    center.x = HALF * (min.x + max.x);
    center.y = HALF * (min.y + max.y);
    center.z = HALF * (min.z + max.z);
    center.m = ZERO;
    //---------------------------------------------------------------------
    *center_dev = center;
    //---------------------------------------------------------------------
  }/* if( (tidx + bidx) == 0 ){ */
  //-----------------------------------------------------------------------
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set box */
  //-----------------------------------------------------------------------
  real diameter = ZERO;  real length;
  length = max.x - min.x;  if( diameter < length )    diameter = length;
  length = max.y - min.y;  if( diameter < length )    diameter = length;
  length = max.z - min.z;  if( diameter < length )    diameter = length;
  diameter = LDEXP(UNITY, (int)CEIL(LOG2(diameter)));
  const real dinv = UNITY / diameter;
  //-----------------------------------------------------------------------
#if 0
  if( (tidx + bidx) == 0 ){
    printf("xmin = %e, xmax = %e, x0 = %e\n", min.x, max.x, HALF * (min.x + max.x));
    printf("ymin = %e, ymax = %e, y0 = %e\n", min.y, max.y, HALF * (min.y + max.y));
    printf("zmin = %e, zmax = %e, z0 = %e\n", min.z, max.z, HALF * (min.z + max.z));
    printf("diameter = %e by %d particles\n", diameter, num);
  }
#endif
  //-----------------------------------------------------------------------
  min.x = HALF * (min.x + max.x - diameter);
  min.y = HALF * (min.y + max.y - diameter);
  min.z = HALF * (min.z + max.z - diameter);
  //-----------------------------------------------------------------------
#if 0
  if( (tidx + bidx) == 0 ){
    printf("xmin = %e, xmax = %e\n", min.x, min.x + diameter);
    printf("ymin = %e, ymax = %e\n", min.y, min.y + diameter);
    printf("zmin = %e, zmax = %e\n", min.z, min.z + diameter);
  }
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate Peano--Hilbert key */
  //-----------------------------------------------------------------------
  const PHint keymax = (PHint)1 << nlevel;
  const  real dscale = dinv * (real)keymax;
  //-----------------------------------------------------------------------
  for(int ih = ihead; ih < itail; ih += NTHREADS_PH){
    //---------------------------------------------------------------------
    const int ii = ih + tidx;
    //---------------------------------------------------------------------
    if( ii < itail ){
      //-------------------------------------------------------------------
      const position pi = ipos[ii];
      PHint px = (PHint)(dscale * (pi.x - min.x));
      PHint py = (PHint)(dscale * (pi.y - min.y));
      PHint pz = (PHint)(dscale * (pi.z - min.z));
#if 0
      PHint tkey = ((dilate3D(px) << 2) | (dilate3D(py) << 1) | (dilate3D(pz)));
#endif
      PHint tkey = 0;
      //-------------------------------------------------------------------
      for(int jj = nlevel - 1; jj >= 0; jj--){
	//-----------------------------------------------------------------
	/* get xi, yi, zi from given position */
	PHint xi = (px >> jj) & 1;
	PHint yi = (py >> jj) & 1;
	PHint zi = (pz >> jj) & 1;
	//-----------------------------------------------------------------
	/* turn px, py, and pz */
	px ^= -( xi & ((!yi) |   zi));
	py ^= -((xi & (  yi  |   zi)) | (yi & (!zi)));
	pz ^= -((xi &  (!yi) & (!zi)) | (yi & (!zi)));
	//-----------------------------------------------------------------
	/* append 3 bits to the key */
	tkey |= ((xi << 2) | ((xi ^ yi) << 1) | ((xi ^ zi) ^ yi)) << (3 * jj);
	//-----------------------------------------------------------------
	/* if zi == 1, then rotate uncyclic (x->z->y->x) */
	if( zi ){	    PHint pt = px;	    px = py;	    py = pz;	  pz = pt;	}
	else{
	  /* if yi == 0, then exchange x and z */
	  if( !yi ){	    PHint pt = px;	    px = pz;	                  pz = pt;	}
	}
	//-----------------------------------------------------------------
      }/* for(int jj = nlevel - 1; jj >= 0; jj--){ */
      //-------------------------------------------------------------------
      idx[ii] = ii;
      key[ii] = tkey;
      //-------------------------------------------------------------------
    }/* if( ii < itail ){ */
    //---------------------------------------------------------------------
  }/* for(int ih = ihead; ih < itail; ih += bnum * NTHREADS_PH){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
extern "C"
muse allocPeanoHilbertKey_dev
(const int num, int **idx_dev, PHint **key_dev, PHint **key_hst, real4 **minall, real4 **maxall, int **gsync0, int **gsync1,
#ifndef CALC_MULTIPOLE_ON_DEVICE
 PHinfo **info_hst,
#endif//CALC_MULTIPOLE_ON_DEVICE
 soaPHsort *dev, soaPHsort *hst,
#ifdef  CUB_AVAILABLE
 soaPHsort *pre, void **temp_storage, int **idx_pre, PHint **key_pre,
#endif//CUB_AVAILABLE
 const deviceProp devProp)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* Peano--Hilbert key of N-body particles */
  //-----------------------------------------------------------------------
  /* the size of the array is set to be a multiple of NTHREADS_PH */
  size_t size = (size_t)num;
  if( (num % NTHREADS_PH) != 0 )
    size += (size_t)(NTHREADS_PH - (num % NTHREADS_PH));
  //-----------------------------------------------------------------------
  /* memory allocation and simple confirmation */
  mycudaMalloc    ((void **)idx_dev, num * sizeof(  int));  alloc.device += num * sizeof(  int);
  mycudaMalloc    ((void **)key_dev, num * sizeof(PHint));  alloc.device += num * sizeof(PHint);
  mycudaMallocHost((void **)key_hst, num * sizeof(PHint));  alloc.host   += num * sizeof(PHint);
  //-----------------------------------------------------------------------
#ifdef  CUB_AVAILABLE
  mycudaMalloc    ((void **)idx_pre, num * sizeof(  int));  alloc.device += num * sizeof(  int);
  mycudaMalloc    ((void **)key_pre, num * sizeof(PHint));  alloc.device += num * sizeof(PHint);
#endif//CUB_AVAILABLE
  //-----------------------------------------------------------------------
  mycudaMalloc((void **)minall, devProp.numSM * NBLOCKS_PER_SM_PH * sizeof(real4));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_PH * sizeof(real4);
  mycudaMalloc((void **)maxall, devProp.numSM * NBLOCKS_PER_SM_PH * sizeof(real4));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_PH * sizeof(real4);
  mycudaMalloc((void **)gsync0, devProp.numSM * NBLOCKS_PER_SM_PH * sizeof(  int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_PH * sizeof(  int);
  mycudaMalloc((void **)gsync1, devProp.numSM * NBLOCKS_PER_SM_PH * sizeof(  int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_PH * sizeof(  int);
  //-----------------------------------------------------------------------
  initGsync_kernel<<<1, devProp.numSM * NBLOCKS_PER_SM_PH>>>(devProp.numSM * NBLOCKS_PER_SM_PH, *gsync0, *gsync1);
  getLastCudaError("initGsync_kernel");
  //-----------------------------------------------------------------------
  dev->idx = *idx_dev;
  dev->key = *key_dev;  hst->key = *key_hst;
  dev->min = *minall;
  dev->max = *maxall;
  dev->gsync0 = *gsync0;
  dev->gsync1 = *gsync1;
  //-----------------------------------------------------------------------
#ifdef  CUB_AVAILABLE
  size_t temp_storage_size = 0;
  pre->idx = *idx_pre;
  pre->key = *key_pre;
  pre->min = *minall;
  pre->max = *maxall;
  pre->gsync0 = *gsync0;
  pre->gsync1 = *gsync1;
  *temp_storage = NULL;
  cub::DeviceRadixSort::SortPairs(*temp_storage, temp_storage_size, pre->key, dev->key, pre->idx, dev->idx, num);
  mycudaMalloc(temp_storage, temp_storage_size);  alloc.device += temp_storage_size;
  dev->temp_storage = *temp_storage;  dev->temp_storage_size = temp_storage_size;
  pre->temp_storage = *temp_storage;  pre->temp_storage_size = temp_storage_size;
#endif//CUB_AVAILABLE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifndef CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
  /* Properties of Peano--Hilbert key of tree cells in the hierarchical structure */
  //-----------------------------------------------------------------------
  /* memory allocation and simple confirmation */
  *info = (PHinfo *)malloc(NUM_PHKEY_LEVEL * sizeof(PHinfo));  if( *info == NULL ){    __KILL__(stderr, "ERROR: failure to allocate info");  }
  alloc.host +=            NUM_PHKEY_LEVEL * sizeof(PHinfo);
  //-----------------------------------------------------------------------
  initPHinfo(*info);
  //-----------------------------------------------------------------------
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* error checking before running the kernel */
  //-----------------------------------------------------------------------
  struct cudaFuncAttributes funcAttr;
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, calcPHkey_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_PH ){
    //---------------------------------------------------------------------
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in calcPHkey_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_PH);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_PH = %d.\n", GPUGEN, GPUVER, NTHREADS_PH);
    fflush (stderr);
    //---------------------------------------------------------------------
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAC ){ */
  int regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS_PH);
  if( regLimit > (MAX_REGISTERS_PER_SM / NTHREADS_PH) )
    regLimit = (MAX_REGISTERS_PER_SM / NTHREADS_PH);
  int memLimit = (16 * 1024) / funcAttr.sharedSizeBytes;
  int Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / NTHREADS_PH) )    Nblck = MAX_THREADS_PER_SM / NTHREADS_PH;
  if( Nblck >   MAX_BLOCKS_PER_SM                      )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / NTHREADS_PH) )    Nblck = ((MAX_WARPS_PER_SM * 32) / NTHREADS_PH);
  if( Nblck != NBLOCKS_PER_SM_PH ){
    //---------------------------------------------------------------------
    __KILL__(stderr, "ERROR: # of blocks per SM for calcPHkey_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_PH defined in src/sort/peano_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_PH, funcAttr.numRegs, REGISTERS_PER_THREAD_PH, GPUGEN, GPUVER);
    //---------------------------------------------------------------------
  }/* if( Nblck != NBLOCKS_PER_SM ){ */
  //-----------------------------------------------------------------------
  if( (devProp.numSM * NBLOCKS_PER_SM_PH) > NTHREADS_PH ){
    __KILL__(stderr, "ERROR: product (%d) of devProp.numSM(%d) * NBLOCKS_PER_SM_PH(%d) must be smaller than NTHREADS_PH(%d) to use shared memory.\n", devProp.numSM * NBLOCKS_PER_SM_PH, devProp.numSM, NBLOCKS_PER_SM_PH, NTHREADS_PH);
  }/* if( (devProp.numSM * NBLOCKS_PER_SM_PH) > NTHREADS_PH ){ */
  //-----------------------------------------------------------------------



  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
void  freePeanoHilbertKey_dev
(int  *idx_dev, PHint  *key_dev, PHint  *key_hst, real4  *minall, real4  *maxall, int  *gsync0, int  *gsync1
#ifndef CALC_MULTIPOLE_ON_DEVICE
 , PHinfo  *info_hst
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef  CUB_AVAILABLE
 , void  *temp_storage, int  *idx_pre, PHint  *key_pre
#endif//CUB_AVAILABLE
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  mycudaFree    (idx_dev);
  mycudaFree    (key_dev);
  mycudaFreeHost(key_hst);
  //-----------------------------------------------------------------------
  mycudaFree(minall);
  mycudaFree(maxall);
  mycudaFree(gsync0);
  mycudaFree(gsync1);
  //-----------------------------------------------------------------------
#ifndef CALC_MULTIPOLE_ON_DEVICE
  free(info_hst);
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifdef  CUB_AVAILABLE
  mycudaFree(temp_storage);
  mycudaFree(idx_pre);
  mycudaFree(key_pre);
#endif//CUB_AVAILABLE
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__global__ void sortParticlesPHcurve_kernel
(const int num, READ_ONLY int * RESTRICT old,
 ulong        * RESTRICT didx , READ_ONLY ulong        * RESTRICT sidx ,
 position     * RESTRICT dpos , READ_ONLY position     * RESTRICT spos ,
 acceleration * RESTRICT dacc , READ_ONLY acceleration * RESTRICT sacc ,
#ifdef  BLOCK_TIME_STEP
 velocity     * RESTRICT dvel , READ_ONLY velocity     * RESTRICT svel ,
 ibody_time   * RESTRICT dtime, READ_ONLY ibody_time   * RESTRICT stime
#else///BLOCK_TIME_STEP
 real         * RESTRICT dvx  , READ_ONLY real         * RESTRICT svx  ,
 real         * RESTRICT dvy  , READ_ONLY real         * RESTRICT svy  ,
 real         * RESTRICT dvz  , READ_ONLY real         * RESTRICT svz
#endif//BLOCK_TIME_STEP
)
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int ii = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( ii < num ){
    //---------------------------------------------------------------------
    /* load old tag */
    const int jj = old[ii];
    //---------------------------------------------------------------------
    didx [ii] = sidx [jj];
    dpos [ii] = spos [jj];
    dacc [ii] = sacc [jj];
#ifdef  BLOCK_TIME_STEP
    dvel [ii] = svel [jj];
    dtime[ii] = stime[jj];
#else///BLOCK_TIME_STEP
    dvx  [ii] = svx  [jj];
    dvy  [ii] = svy  [jj];
    dvz  [ii] = svz  [jj];
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
  }/* if( ii < num ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__global__ void sortParticlesPHcurve_kernel_offset
(const int num, READ_ONLY int * RESTRICT old,
 ulong        * RESTRICT didx , READ_ONLY ulong        * RESTRICT sidx ,
 position     * RESTRICT dpos , READ_ONLY position     * RESTRICT spos ,
 acceleration * RESTRICT dacc , READ_ONLY acceleration * RESTRICT sacc ,
#ifdef  BLOCK_TIME_STEP
 velocity     * RESTRICT dvel , READ_ONLY velocity     * RESTRICT svel ,
 ibody_time   * RESTRICT dtime, READ_ONLY ibody_time   * RESTRICT stime
#else///BLOCK_TIME_STEP
 real         * RESTRICT dvx  , READ_ONLY real         * RESTRICT svx  ,
 real         * RESTRICT dvy  , READ_ONLY real         * RESTRICT svy  ,
 real         * RESTRICT dvz  , READ_ONLY real         * RESTRICT svz
#endif//BLOCK_TIME_STEP
 , const int offset
)
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int ii = offset + GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( ii < num ){
    //---------------------------------------------------------------------
    /* load old tag */
    const int jj = old[ii];
    //---------------------------------------------------------------------
    didx [ii] = sidx [jj];
    dpos [ii] = spos [jj];
    dacc [ii] = sacc [jj];
#ifdef  BLOCK_TIME_STEP
    dvel [ii] = svel [jj];
    dtime[ii] = stime[jj];
#else///BLOCK_TIME_STEP
    dvx  [ii] = svx  [jj];
    dvy  [ii] = svy  [jj];
    dvz  [ii] = svz  [jj];
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
  }/* if( ii < num ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
extern "C"
void sortParticlesPHcurve_dev(const int num, iparticle * RESTRICT src, iparticle * RESTRICT dst, soaPHsort dev, const deviceProp devProp
#ifdef  CUB_AVAILABLE
			      , soaPHsort pre
#endif//CUB_AVAILABLE
#ifndef MAKE_TREE_ON_DEVICE
			      , const soaPHsort hst
#endif//MAKE_TREE_ON_DEVICE
			      , struct timeval *start
#ifdef EXEC_BENCHMARK
			      , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#ifndef EXEC_BENCHMARK
  checkCudaErrors(cudaDeviceSynchronize());
  gettimeofday(start, NULL);
#else///EXEC_BENCHMARK
  initStopwatch();
  *start = _benchIni;
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* generate Peano--Hilbert key */
  //-----------------------------------------------------------------------
#ifdef  CUB_AVAILABLE
  calcPHkey_kernel<<<devProp.numSM * NBLOCKS_PER_SM_PH, NTHREADS_PH>>>(num, (*src).pos, dev.min, dev.max, MAXIMUM_PHKEY_LEVEL, pre.idx, pre.key, dev.gsync0, dev.gsync1
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
								       , (*dst).encBall
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
								       );
#else///CUB_AVAILABLE
  calcPHkey_kernel<<<devProp.numSM * NBLOCKS_PER_SM_PH, NTHREADS_PH>>>(num, (*src).pos, dev.min, dev.max, MAXIMUM_PHKEY_LEVEL, dev.idx, dev.key, dev.gsync0, dev.gsync1
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
								       , (*dst).encBall
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
								       );
#endif//CUB_AVAILABLE
  getLastCudaError("calcPHkey_kernel");
#if 0
  checkCudaErrors(cudaDeviceSynchronize());
  fprintf(stderr, "kill to debug calcPHkey_kernel\n");
  fflush(NULL);
  exit(0);
#endif
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  checkCudaErrors(cudaMemcpy((*dst).encBall_hst, (*dst).encBall, sizeof(position), cudaMemcpyDeviceToHost));
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->genPHkey_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* sort Peano--Hilbert key */
  //-----------------------------------------------------------------------
#ifdef  CUB_AVAILABLE
  cub::DeviceRadixSort::SortPairs(dev.temp_storage, dev.temp_storage_size, pre.key, dev.key, pre.idx, dev.idx, num);
#else///CUB_AVAILABLE
  thrust::stable_sort_by_key((thrust::device_ptr<PHint>)dev.key, (thrust::device_ptr<PHint>)(dev.key + num), (thrust::device_ptr<int>)dev.idx);
#endif//CUB_AVAILABLE
#ifndef MAKE_TREE_ON_DEVICE
  checkCudaErrors(cudaMemcpy(hst.key, dev.key, num * sizeof(PHint), cudaMemcpyDeviceToHost));
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->rsortKey_library));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* sort N-body particles using Peano--Hilbert key */
  //-----------------------------------------------------------------------
  int Nrem = BLOCKSIZE(num, NTHREADS_PHSORT);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    sortParticlesPHcurve_kernel<<<Nrem, NTHREADS_PHSORT>>>
      (num, dev.idx,
       (*dst).idx , (*src).idx,
       (*dst).pos , (*src).pos,
       (*dst).acc , (*src).acc,
#ifdef  BLOCK_TIME_STEP
       (*dst).vel , (*src).vel,
       (*dst).time, (*src).time
#else///BLOCK_TIME_STEP
       (*dst).vx  , (*src).vx,
       (*dst).vy  , (*src).vy,
       (*dst).vz  , (*src).vz
#endif//BLOCK_TIME_STEP
       );
  else{
    //---------------------------------------------------------------------
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;
    //---------------------------------------------------------------------
    for(int iter = 0; iter < Niter; iter++){
      //-------------------------------------------------------------------
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;
      //-------------------------------------------------------------------
      int Nsub = Nblck * NTHREADS_PHSORT;
      sortParticlesPHcurve_kernel_offset<<<Nblck, NTHREADS_PHSORT>>>
	(num, dev.idx,
	 (*dst).idx , (*src).idx,
	 (*dst).pos , (*src).pos,
	 (*dst).acc , (*src).acc,
#ifdef  BLOCK_TIME_STEP
	 (*dst).vel , (*src).vel,
	 (*dst).time, (*src).time
#else///BLOCK_TIME_STEP
	 (*dst).vx  , (*src).vx,
	 (*dst).vy  , (*src).vy,
	 (*dst).vz  , (*src).vz
#endif//BLOCK_TIME_STEP
	 , hidx
	 );
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
  getLastCudaError("sortParticlesPHcurve_kernel");
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->sortBody_kernel));
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* swap the list structure */
  //-----------------------------------------------------------------------
  iparticle _tmp;
  _tmp = *src;
  *src = *dst;
  *dst = _tmp;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
#ifndef HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->sortParticlesPHcurve));
#else///HUNT_MAKE_PARAMETER
  elapsed->sortParticlesPHcurve += elapsed->genPHkey_kernel + elapsed->rsortKey_library + elapsed->sortBody_kernel;
#endif//HUNT_MAKE_PARAMETER
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  defined(BLOCK_TIME_STEP) && defined(LOCALIZE_I_PARTICLES) && defined(BRUTE_FORCE_LOCALIZATION)
//-------------------------------------------------------------------------
__global__ void copySortedParticles_kernel(const int Ni, READ_ONLY position * RESTRICT src, position * RESTRICT dst)
{
  //-----------------------------------------------------------------------
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( gidx < Ni )
    dst[gidx] = src[gidx];
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
void copySortedParticles_dev(const int Ni, const iparticle pi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  int Nrem = BLOCKSIZE(Ni, 1024);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    copySortedParticles_kernel<<<Nrem, 1024>>>(Ni, pi.pos, pi.jpos);
  else{
    //---------------------------------------------------------------------
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;
    //---------------------------------------------------------------------
    for(int iter = 0; iter < Niter; iter++){
      //-------------------------------------------------------------------
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;
      //-------------------------------------------------------------------
      int Nsub = Nblck * 1024;
      copySortedParticles_kernel<<<Nblck, 1024>>>(Nsub, &pi.pos[hidx], &pi.jpos[hidx]);
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  /* copySortedParticles_kernel<<<BLOCKSIZE(Ni, 1024), 1024>>>(Ni, pi.pos, pi.jpos); */
  getLastCudaError("copySortedParticles_dev");
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//defined(BLOCK_TIME_STEP) && defined(LOCALIZE_I_PARTICLES) && defined(BRUTE_FORCE_LOCALIZATION)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------
/* Properties of Peano--Hilbert key of tree cells in the hierarchical structure */
/* NOTE: syncthreads() and if statements of ii < NUM_PHKEY_LEVEL can be removed; however, remained for safety */
//-------------------------------------------------------------------------
__global__ void initPHinfo_kernel(PHinfo *gm)
{
  //-----------------------------------------------------------------------
  const int ii = THREADIDX_X1D;
  __shared__ PHinfo sm[NUM_PHKEY_LEVEL];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* initialize fundamental properties */
  //-----------------------------------------------------------------------
  if( ii < NUM_PHKEY_LEVEL ){
    //---------------------------------------------------------------------
    sm[ii].level = MAXIMUM_PHKEY_LEVEL - ii;
    //---------------------------------------------------------------------
    ulong ntmp = 1;
    int jj = NUM_PHKEY_LEVEL - 1;
    if( ii == 0 )
      while( true ){
	//-----------------------------------------------------------------
	sm[jj].nmax = (int)ntmp;
	//-----------------------------------------------------------------
	ntmp *= NLEAF;    jj--;
	if( ntmp > INT_MAX )	ntmp = INT_MAX;
	//-----------------------------------------------------------------
	if( jj < 0 )	break;
	//-----------------------------------------------------------------
      }/* while( true ){ */
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  __syncthreads();
  if( ii < NUM_PHKEY_LEVEL )
    gm[ii] = sm[ii];
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
void initPHinfo_dev(PHinfo *info)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  initPHinfo_kernel<<<1, NUM_PHKEY_LEVEL>>>(info);
  getLastCudaError("initPHinfo_kernel");
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------
