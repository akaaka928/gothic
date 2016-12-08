/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:48:22
 *                                                                       *
 *    Generation of enclosing ball containing all N-body particles       *
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
#include "macro.h"
#include "cudalib.h"
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/gsync_dev.cu"
//-------------------------------------------------------------------------
/* in L1 cache preferred configuration, capacity of shared memory is 16KiB per SM */
/* real4 smem[NTHREADS_EB] corresponds 16 * NTHREADS_EB bytes */
#define NBLOCKS_PER_SM_EB (1024 / NTHREADS_EB)
//-------------------------------------------------------------------------
#define REGISTERS_PER_THREAD_EB (40)
/* calcPHkey_kernel uses 32 registers @ Tesla M2090, Ttot = 1024 (registers are spilled to local memory) */
/* calcPHkey_kernel uses 47 registers @ Tesla M2090, Ttot =  512 */
/* calcPHkey_kernel uses 36 registers @ Tesla M2090, Ttot =  128, 256 */
#   if  GPUVER == 20
#undef  REGISTERS_PER_THREAD_EB
#          if  NTHREADS_EB == 1024
#define REGISTERS_PER_THREAD_EB (32)
#       else///NTHREADS_EB == 1024
#          if  NTHREADS_EB ==  512
#define REGISTERS_PER_THREAD_EB (47)
#       else///NTHREADS_EB ==  512
#define REGISTERS_PER_THREAD_EB (36)
#       endif//NTHREADS_EB ==  512
#       endif//NTHREADS_EB == 1024
#endif//GPUVER == 20
/* calcPHkey_kernel uses 38 registers @ Tesla K20X */
/* #   if  GPUVER == 35 */
/* #undef  REGISTERS_PER_THREAD_EB */
/* #define REGISTERS_PER_THREAD_EB (38) */
/* #endif//GPUVER == 35 */
/* calcPHkey_kernel uses 40 registers @ GTX 750 Ti */
/* #   if  GPUVER == 50 */
/* #undef  REGISTERS_PER_THREAD_EB */
/* #define REGISTERS_PER_THREAD_EB (40) */
/* #endif//GPUVER == 50 */
/* calcPHkey_kernel uses 40 registers @ GTX 970 */
/* #   if  GPUVER == 52 */
/* #undef  REGISTERS_PER_THREAD_EB */
/* #define REGISTERS_PER_THREAD_EB (40) */
/* #endif//GPUVER == 52 */
//-------------------------------------------------------------------------
/* limitation from number of registers */
#   if  NBLOCKS_PER_SM_EB > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_EB * NTHREADS_EB))
#undef  NBLOCKS_PER_SM_EB
#define NBLOCKS_PER_SM_EB   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_EB * NTHREADS_EB))
#endif//NBLOCKS_PER_SM_EB > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_EB * NTHREADS_EB))
//-------------------------------------------------------------------------
/* maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_EB > (MAX_THREADS_PER_SM / NTHREADS_EB)
#undef  NBLOCKS_PER_SM_EB
#define NBLOCKS_PER_SM_EB   (MAX_THREADS_PER_SM / NTHREADS_EB)
#endif//NBLOCKS_PER_SM_EB > (MAX_THREADS_PER_SM / NTHREADS_EB)
//-------------------------------------------------------------------------
/* maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_EB > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_EB
#define NBLOCKS_PER_SM_EB   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_EB > MAX_BLOCKS_PER_SM
//-------------------------------------------------------------------------
/* maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_EB > ((MAX_WARPS_PER_SM * 32) / NTHREADS_EB)
#undef  NBLOCKS_PER_SM_EB
#define NBLOCKS_PER_SM_EB   ((MAX_WARPS_PER_SM * 32) / NTHREADS_EB)
#endif//NBLOCKS_PER_SM_EB > ((MAX_WARPS_PER_SM * 32) / NTHREADS_EB)
//-------------------------------------------------------------------------
/* # of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_EB < 1
#undef  NBLOCKS_PER_SM_EB
#define NBLOCKS_PER_SM_EB  (1)
#endif//NBLOCKS_PER_SM_EB < 1
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct __align__(8){
  real val;
  int idx;
} rloc;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* 1. estimate box size */
/* 2. calculate Peano--Hilbert key */
//-------------------------------------------------------------------------
__global__ void __launch_bounds__(NTHREADS_EB, NBLOCKS_PER_SM_EB) getMinLocMaxLoc
     (const int num, READ_ONLY position * RESTRICT ipos, rloc * RESTRICT min_all, rloc * RESTRICT max_all,
      const int nlevel, int * RESTRICT idx, PHint * RESTRICT key, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  //-----------------------------------------------------------------------
  /* static allocation of the shared memory */
  //-----------------------------------------------------------------------
  __shared__ rloc smem[NTHREADS_EB];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int bnum =   GRIDDIM_X1D;
  //-----------------------------------------------------------------------
  const int uidx = tidx + (NTHREADS_EB >> 1);
  const int hidx = tidx - (tidx & (warpSize - 1));
  //-----------------------------------------------------------------------
  const int ihead = (num *      bidx ) / bnum;
  const int itail = (num * (1 + bidx)) / bnum;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate required box size to contain all N-body particles in the local domain */
  //-----------------------------------------------------------------------
  rloc xmin = { REAL_MAX, INT_MAX};  rloc ymin = { REAL_MAX, INT_MAX};  rloc zmin = { REAL_MAX, INT_MAX};
  rloc xmax = {-REAL_MAX, INT_MAX};  rloc ymax = {-REAL_MAX, INT_MAX};  rloc zmax = {-REAL_MAX, INT_MAX};
  //-----------------------------------------------------------------------
  for(int ih = ihead; ih < itail; ih += NTHREADS_EB){
    //---------------------------------------------------------------------
    const int ii = ih + tidx;
    //---------------------------------------------------------------------
    if( ii < itail ){
      //-------------------------------------------------------------------
      const position pi = ipos[ii];
      //-------------------------------------------------------------------
      if( pi.x < xmin.val ){	xmin.val = pi.x;	xmin.idx = ii;      }
      if( pi.y < ymin.val ){	ymin.val = pi.y;	ymin.idx = ii;      }
      if( pi.z < zmin.val ){	zmin.val = pi.z;	zmin.idx = ii;      }
      //-------------------------------------------------------------------
      if( pi.x > xmax.val ){	xmax.val = pi.x;	xmax.idx = ii;      }
      if( pi.y > ymax.val ){	ymax.val = pi.y;	ymax.idx = ii;      }
      if( pi.z > zmax.val ){	zmax.val = pi.z;	zmax.idx = ii;      }
      //-------------------------------------------------------------------
    }/* if( ii < itail ){ */
    //---------------------------------------------------------------------
  }/* for(int ih = ihead; ih < itail; ih += NTHREADS_EB){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* get the box size within a warp */
  //-----------------------------------------------------------------------
  rloc tmp;
  //-----------------------------------------------------------------------
  /* get minimum */
  smem[tidx] = xmin;
  tmp = smem[tidx ^  1];  if( tmp.val < xmin.val )    xmin = tmp;  smem[tidx] = xmin;
  tmp = smem[tidx ^  2];  if( tmp.val < xmin.val )    xmin = tmp;  smem[tidx] = xmin;
  tmp = smem[tidx ^  4];  if( tmp.val < xmin.val )    xmin = tmp;  smem[tidx] = xmin;
  tmp = smem[tidx ^  8];  if( tmp.val < xmin.val )    xmin = tmp;  smem[tidx] = xmin;
  tmp = smem[tidx ^ 16];  if( tmp.val < xmin.val )    xmin = tmp;  smem[tidx] = xmin;
  smem[tidx] = ymin;
  tmp = smem[tidx ^  1];  if( tmp.val < ymin.val )    ymin = tmp;  smem[tidx] = ymin;
  tmp = smem[tidx ^  2];  if( tmp.val < ymin.val )    ymin = tmp;  smem[tidx] = ymin;
  tmp = smem[tidx ^  4];  if( tmp.val < ymin.val )    ymin = tmp;  smem[tidx] = ymin;
  tmp = smem[tidx ^  8];  if( tmp.val < ymin.val )    ymin = tmp;  smem[tidx] = ymin;
  tmp = smem[tidx ^ 16];  if( tmp.val < ymin.val )    ymin = tmp;  smem[tidx] = ymin;
  smem[tidx] = zmin;
  tmp = smem[tidx ^  1];  if( tmp.val < zmin.val )    zmin = tmp;  smem[tidx] = zmin;
  tmp = smem[tidx ^  2];  if( tmp.val < zmin.val )    zmin = tmp;  smem[tidx] = zmin;
  tmp = smem[tidx ^  4];  if( tmp.val < zmin.val )    zmin = tmp;  smem[tidx] = zmin;
  tmp = smem[tidx ^  8];  if( tmp.val < zmin.val )    zmin = tmp;  smem[tidx] = zmin;
  tmp = smem[tidx ^ 16];  if( tmp.val < zmin.val )    zmin = tmp;  smem[tidx] = zmin;
  //-----------------------------------------------------------------------
  /* get maximum */
  smem[tidx] = xmax;
  tmp = smem[tidx ^  1];  if( tmp.val > xmax.val )    xmax = tmp;  smem[tidx] = xmax;
  tmp = smem[tidx ^  2];  if( tmp.val > xmax.val )    xmax = tmp;  smem[tidx] = xmax;
  tmp = smem[tidx ^  4];  if( tmp.val > xmax.val )    xmax = tmp;  smem[tidx] = xmax;
  tmp = smem[tidx ^  8];  if( tmp.val > xmax.val )    xmax = tmp;  smem[tidx] = xmax;
  tmp = smem[tidx ^ 16];  if( tmp.val > xmax.val )    xmax = tmp;  smem[tidx] = xmax;
  smem[tidx] = ymax;
  tmp = smem[tidx ^  1];  if( tmp.val > ymax.val )    ymax = tmp;  smem[tidx] = ymax;
  tmp = smem[tidx ^  2];  if( tmp.val > ymax.val )    ymax = tmp;  smem[tidx] = ymax;
  tmp = smem[tidx ^  4];  if( tmp.val > ymax.val )    ymax = tmp;  smem[tidx] = ymax;
  tmp = smem[tidx ^  8];  if( tmp.val > ymax.val )    ymax = tmp;  smem[tidx] = ymax;
  tmp = smem[tidx ^ 16];  if( tmp.val > ymax.val )    ymax = tmp;  smem[tidx] = ymax;
  smem[tidx] = zmax;
  tmp = smem[tidx ^  1];  if( tmp.val > zmax.val )    zmax = tmp;  smem[tidx] = zmax;
  tmp = smem[tidx ^  2];  if( tmp.val > zmax.val )    zmax = tmp;  smem[tidx] = zmax;
  tmp = smem[tidx ^  4];  if( tmp.val > zmax.val )    zmax = tmp;  smem[tidx] = zmax;
  tmp = smem[tidx ^  8];  if( tmp.val > zmax.val )    zmax = tmp;  smem[tidx] = zmax;
  tmp = smem[tidx ^ 16];  if( tmp.val > zmax.val )    zmax = tmp;  smem[tidx] = zmax;
  //-----------------------------------------------------------------------
/* #ifndef NDEBUG */
/*   if( (tidx & (warpSize - 1)) == 0 ){ */
/*     printf("min(%d:%d): %f, %f, %f\n", bidx, tidx / warpSize, xmin.val, ymin.val, zmin.val); */
/*     printf("max(%d:%d): %f, %f, %f\n", bidx, tidx / warpSize, xmax.val, ymax.val, zmax.val); */
/*   } */
/* #endif//NDEBUG */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* get the box size within a block */
  //-----------------------------------------------------------------------
  /* 1024 = 32^2 = warpSize^2 is the maximum of the number of threads; therefore, Nwarp <= warpSize */
  const int Nwarp = NTHREADS_EB >> 5;/* := NTHREADS_EB / warpSize; */
  __syncthreads();
  //-----------------------------------------------------------------------
  if( hidx == tidx ){
    //---------------------------------------------------------------------
    smem[                     (tidx >> 5)] = xmin;/* := smem[                     (tidx / warpSize)] = xmin; */
    smem[(NTHREADS_EB >> 1) + (tidx >> 5)] = xmax;/* := smem[(NTHREADS_EB >> 1) + (tidx / warpSize)] = xmax; */
    //---------------------------------------------------------------------
  }/* if( hidx == tidx ){ */
  __syncthreads();
  //-----------------------------------------------------------------------
  if( tidx < Nwarp ){
    //---------------------------------------------------------------------
    /* get minimum */
    xmin = smem[tidx];
#   if  NTHREADS_EB >=   64
    tmp = smem[tidx ^  1];  if( tmp.val < xmin.val )	xmin = tmp;	smem[tidx] = xmin;
#   if	NTHREADS_EB >=	128
    tmp = smem[tidx ^  2];  if( tmp.val < xmin.val )	xmin = tmp;	smem[tidx] = xmin;
#   if	NTHREADS_EB >=	256
    tmp = smem[tidx ^  4];  if( tmp.val < xmin.val )	xmin = tmp;	smem[tidx] = xmin;
#   if	NTHREADS_EB >=	512
    tmp = smem[tidx ^  8];  if( tmp.val < xmin.val )	xmin = tmp;	smem[tidx] = xmin;
#   if	NTHREADS_EB == 1024
    tmp = smem[tidx ^ 16];  if( tmp.val < xmin.val )	xmin = tmp;	smem[tidx] = xmin;
#endif//NTHREADS_EB == 1024
#endif//NTHREADS_EB >=  512
#endif//NTHREADS_EB >=  256
#endif//NTHREADS_EB >=  128
#endif//NTHREADS_EB >=   64
    xmin = smem[0];
    //---------------------------------------------------------------------
    /* get maximum */
    xmax = smem[uidx];
#   if  NTHREADS_EB >=   64
    tmp = smem[uidx ^  1];  if( tmp.val > xmax.val )	xmax = tmp;	smem[uidx] = xmax;
#   if	NTHREADS_EB >=	128
    tmp = smem[uidx ^  2];  if( tmp.val > xmax.val )	xmax = tmp;	smem[uidx] = xmax;
#   if	NTHREADS_EB >=	256
    tmp = smem[uidx ^  4];  if( tmp.val > xmax.val )	xmax = tmp;	smem[uidx] = xmax;
#   if	NTHREADS_EB >=	512
    tmp = smem[uidx ^  8];  if( tmp.val > xmax.val )	xmax = tmp;	smem[uidx] = xmax;
#   if	NTHREADS_EB == 1024
    tmp = smem[uidx ^ 16];  if( tmp.val > xmax.val )	xmax = tmp;	smem[uidx] = xmax;
#endif//NTHREADS_EB == 1024
#endif//NTHREADS_EB >=  512
#endif//NTHREADS_EB >=  256
#endif//NTHREADS_EB >=  128
#endif//NTHREADS_EB >=   64
    xmax = smem[NTHREADS_PH >> 1];
    //---------------------------------------------------------------------
  }/* if( tidx < Nwarp ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  if( hidx == tidx ){
    //---------------------------------------------------------------------
    smem[                     (tidx >> 5)] = ymin;/* := smem[                     (tidx / warpSize)] = ymin; */
    smem[(NTHREADS_EB >> 1) + (tidx >> 5)] = ymax;/* := smem[(NTHREADS_EB >> 1) + (tidx / warpSize)] = ymax; */
    //---------------------------------------------------------------------
  }/* if( hidx == tidx ){ */
  __syncthreads();
  //-----------------------------------------------------------------------
  if( tidx < Nwarp ){
    //---------------------------------------------------------------------
    /* get minimum */
    ymin = smem[tidx];
#   if  NTHREADS_EB >=   64
    tmp = smem[tidx ^  1];  if( tmp.val < ymin.val )	ymin = tmp;	smem[tidx] = ymin;
#   if	NTHREADS_EB >=	128
    tmp = smem[tidx ^  2];  if( tmp.val < ymin.val )	ymin = tmp;	smem[tidx] = ymin;
#   if	NTHREADS_EB >=	256
    tmp = smem[tidx ^  4];  if( tmp.val < ymin.val )	ymin = tmp;	smem[tidx] = ymin;
#   if	NTHREADS_EB >=	512
    tmp = smem[tidx ^  8];  if( tmp.val < ymin.val )	ymin = tmp;	smem[tidx] = ymin;
#   if	NTHREADS_EB == 1024
    tmp = smem[tidx ^ 16];  if( tmp.val < ymin.val )	ymin = tmp;	smem[tidx] = ymin;
#endif//NTHREADS_EB == 1024
#endif//NTHREADS_EB >=  512
#endif//NTHREADS_EB >=  256
#endif//NTHREADS_EB >=  128
#endif//NTHREADS_EB >=   64
    ymin = smem[0];
    //---------------------------------------------------------------------
    /* get maximum */
    ymax = smem[uidx];
#   if  NTHREADS_EB >=   64
    tmp = smem[uidx ^  1];  if( tmp.val > ymax.val )	ymax = tmp;	smem[uidx] = ymax;
#   if	NTHREADS_EB >=	128
    tmp = smem[uidx ^  2];  if( tmp.val > ymax.val )	ymax = tmp;	smem[uidx] = ymax;
#   if	NTHREADS_EB >=	256
    tmp = smem[uidx ^  4];  if( tmp.val > ymax.val )	ymax = tmp;	smem[uidx] = ymax;
#   if	NTHREADS_EB >=	512
    tmp = smem[uidx ^  8];  if( tmp.val > ymax.val )	ymax = tmp;	smem[uidx] = ymax;
#   if	NTHREADS_EB == 1024
    tmp = smem[uidx ^ 16];  if( tmp.val > ymax.val )	ymax = tmp;	smem[uidx] = ymax;
#endif//NTHREADS_EB == 1024
#endif//NTHREADS_EB >=  512
#endif//NTHREADS_EB >=  256
#endif//NTHREADS_EB >=  128
#endif//NTHREADS_EB >=   64
    ymax = smem[NTHREADS_PH >> 1];
    //---------------------------------------------------------------------
  }/* if( tidx < Nwarp ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  if( hidx == tidx ){
    //---------------------------------------------------------------------
    smem[                     (tidx >> 5)] = zmin;/* := smem[                     (tidx / warpSize)] = zmin; */
    smem[(NTHREADS_EB >> 1) + (tidx >> 5)] = zmax;/* := smem[(NTHREADS_EB >> 1) + (tidx / warpSize)] = zmax; */
    //---------------------------------------------------------------------
  }/* if( hidx == tidx ){ */
  __syncthreads();
  //-----------------------------------------------------------------------
  if( tidx < Nwarp ){
    //---------------------------------------------------------------------
    /* get minimum */
    zmin = smem[tidx];
#   if  NTHREADS_EB >=   64
    tmp = smem[tidx ^  1];  if( tmp.val < zmin.val )	zmin = tmp;	smem[tidx] = zmin;
#   if	NTHREADS_EB >=	128
    tmp = smem[tidx ^  2];  if( tmp.val < zmin.val )	zmin = tmp;	smem[tidx] = zmin;
#   if	NTHREADS_EB >=	256
    tmp = smem[tidx ^  4];  if( tmp.val < zmin.val )	zmin = tmp;	smem[tidx] = zmin;
#   if	NTHREADS_EB >=	512
    tmp = smem[tidx ^  8];  if( tmp.val < zmin.val )	zmin = tmp;	smem[tidx] = zmin;
#   if	NTHREADS_EB == 1024
    tmp = smem[tidx ^ 16];  if( tmp.val < zmin.val )	zmin = tmp;	smem[tidx] = zmin;
#endif//NTHREADS_EB == 1024
#endif//NTHREADS_EB >=  512
#endif//NTHREADS_EB >=  256
#endif//NTHREADS_EB >=  128
#endif//NTHREADS_EB >=   64
    zmin = smem[0];
    //---------------------------------------------------------------------
    /* get maximum */
    zmax = smem[uidx];
#   if  NTHREADS_EB >=   64
    tmp = smem[uidx ^  1];  if( tmp.val > zmax.val )	zmax = tmp;	smem[uidx] = zmax;
#   if	NTHREADS_EB >=	128
    tmp = smem[uidx ^  2];  if( tmp.val > zmax.val )	zmax = tmp;	smem[uidx] = zmax;
#   if	NTHREADS_EB >=	256
    tmp = smem[uidx ^  4];  if( tmp.val > zmax.val )	zmax = tmp;	smem[uidx] = zmax;
#   if	NTHREADS_EB >=	512
    tmp = smem[uidx ^  8];  if( tmp.val > zmax.val )	zmax = tmp;	smem[uidx] = zmax;
#   if	NTHREADS_EB == 1024
    tmp = smem[uidx ^ 16];  if( tmp.val > zmax.val )	zmax = tmp;	smem[uidx] = zmax;
#endif//NTHREADS_EB == 1024
#endif//NTHREADS_EB >=  512
#endif//NTHREADS_EB >=  256
#endif//NTHREADS_EB >=  128
#endif//NTHREADS_EB >=   64
    zmax = smem[NTHREADS_PH >> 1];
    //---------------------------------------------------------------------
  }/* if( tidx < Nwarp ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
/* #ifndef NDEBUG */
/*   if( tidx  == 0 ){ */
/*     printf("min(%d): %f, %f, %f\n", bidx, xmin.val, ymin.val, zmin.val); */
/*     printf("max(%d): %f, %f, %f\n", bidx, xmax.val, ymax.val, zmax.val); */
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
      min_all[bidx] = xmin;      min_all[bidx + bnum] = ymin;      min_all[bidx + 2 * bnum] = zmin;
      max_all[bidx] = xmax;      max_all[bidx + bnum] = ymax;      max_all[bidx + 2 * bnum] = zmax;
      //-------------------------------------------------------------------
#if 0
      printf("%d-th block:\t[%e:%e], [%e:%e], [%e:%e]\n", bidx, xmin.val, xmax.val, ymin.val, ymax.val, zmin.val, zmax.val);
#endif
      //-------------------------------------------------------------------
    }/* if( tidx == 0 ){ */
    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* get minimum */
    //---------------------------------------------------------------------
    /* bnum is set to (# of SMs) * NBLOCKS_PER_SM_EB <= NTHREADS_EB */
    /* ---> # of working block is 1 (bidx = 0) */
    if( bidx == 0 ){
      //-------------------------------------------------------------------
      smem[tidx] = xmin;
      //-------------------------------------------------------------------
      if( tidx < bnum  ){
	//-----------------------------------------------------------------
	/* reduction within a warp */
	//-----------------------------------------------------------------
	xmin = min_all[tidx];
	smem[tidx] = xmin;
	if( bnum >  1 ){
	  tmp = smem[tidx ^  1];  if( tmp.val < xmin.val )	   xmin = tmp;	smem[tidx] = xmin;
	if( bnum >  2 ){
	  tmp = smem[tidx ^  2];  if( tmp.val < xmin.val )	   xmin = tmp;	smem[tidx] = xmin;
	if( bnum >  4 ){
	  tmp = smem[tidx ^  4];  if( tmp.val < xmin.val )	   xmin = tmp;	smem[tidx] = xmin;
	if( bnum >  8 ){
	  tmp = smem[tidx ^  8];  if( tmp.val < xmin.val )	   xmin = tmp;	smem[tidx] = xmin;
	if( bnum > 16 ){
	  tmp = smem[tidx ^ 16];  if( tmp.val < xmin.val )	   xmin = tmp;	smem[tidx] = xmin;
	}}}}}
	//-----------------------------------------------------------------
	ymin = min_all[tidx + bnum];
	smem[tidx] = ymin;
	if( bnum >  1 ){
	  tmp = smem[tidx ^  1];  if( tmp.val < ymin.val )	   ymin = tmp;	smem[tidx] = ymin;
	if( bnum >  2 ){
	  tmp = smem[tidx ^  2];  if( tmp.val < ymin.val )	   ymin = tmp;	smem[tidx] = ymin;
	if( bnum >  4 ){
	  tmp = smem[tidx ^  4];  if( tmp.val < ymin.val )	   ymin = tmp;	smem[tidx] = ymin;
	if( bnum >  8 ){
	  tmp = smem[tidx ^  8];  if( tmp.val < ymin.val )	   ymin = tmp;	smem[tidx] = ymin;
	if( bnum > 16 ){
	  tmp = smem[tidx ^ 16];  if( tmp.val < ymin.val )	   ymin = tmp;	smem[tidx] = ymin;
	}}}}}
	//-----------------------------------------------------------------
	zmin = min_all[tidx + 2 * bnum];
	smem[tidx] = zmin;
	if( bnum >  1 ){
	  tmp = smem[tidx ^  1];  if( tmp.val < zmin.val )	   zmin = tmp;	smem[tidx] = zmin;
	if( bnum >  2 ){
	  tmp = smem[tidx ^  2];  if( tmp.val < zmin.val )	   zmin = tmp;	smem[tidx] = zmin;
	if( bnum >  4 ){
	  tmp = smem[tidx ^  4];  if( tmp.val < zmin.val )	   zmin = tmp;	smem[tidx] = zmin;
	if( bnum >  8 ){
	  tmp = smem[tidx ^  8];  if( tmp.val < zmin.val )	   zmin = tmp;	smem[tidx] = zmin;
	if( bnum > 16 ){
	  tmp = smem[tidx ^ 16];  if( tmp.val < zmin.val )	   zmin = tmp;	smem[tidx] = zmin;
	}}}}}
	//-----------------------------------------------------------------
#if 0
	if( tidx == hidx )
	  printf("%d-th warp (min):\t% e, % e, % e\n", tidx >> 5, xmin.val, ymin.val, zmin.val);
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
	  smem[(tidx >> 5)] = xmin;/* := (tidx / warpSize) */
	__syncthreads();
	//-----------------------------------------------------------------
	/* bnum <= 1024 --->>> Ndata <= 32 */
	if( tidx < Ndata ){
	  //---------------------------------------------------------------
	  xmin = smem[tidx];
	  if( Ndata >  1 ){
	    tmp = smem[tidx ^  1];  if( tmp.val < xmin.val )	 xmin = tmp;	 smem[tidx] = xmin;
	  if( Ndata >  2 ){
	    tmp = smem[tidx ^  2];  if( tmp.val < xmin.val )	 xmin = tmp;	 smem[tidx] = xmin;
	  if( Ndata >  4 ){
	    tmp = smem[tidx ^  4];  if( tmp.val < xmin.val )	 xmin = tmp;	 smem[tidx] = xmin;
	  if( Ndata >  8 ){
	    tmp = smem[tidx ^  8];  if( tmp.val < xmin.val )	 xmin = tmp;	 smem[tidx] = xmin;
	  if( Ndata > 16 ){
	    tmp = smem[tidx ^ 16];  if( tmp.val < xmin.val )	 xmin = tmp;	 smem[tidx] = xmin;
	  }}}}}
	  //---------------------------------------------------------------
	}/* if( tidx < Ndata ){ */
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	if( (hidx == tidx) && (tidx < bnum) )
	  smem[(tidx >> 5)] = ymin;/* := (tidx / warpSize) */
	__syncthreads();
	//-----------------------------------------------------------------
	/* bnum <= 1024 --->>> Ndata <= 32 */
	if( tidx < Ndata ){
	  //---------------------------------------------------------------
	  ymin = smem[tidx];
	  if( Ndata >  1 ){
	    tmp = smem[tidx ^  1];  if( tmp.val < ymin.val )	 ymin = tmp;	 smem[tidx] = ymin;
	  if( Ndata >  2 ){
	    tmp = smem[tidx ^  2];  if( tmp.val < ymin.val )	 ymin = tmp;	 smem[tidx] = ymin;
	  if( Ndata >  4 ){
	    tmp = smem[tidx ^  4];  if( tmp.val < ymin.val )	 ymin = tmp;	 smem[tidx] = ymin;
	  if( Ndata >  8 ){
	    tmp = smem[tidx ^  8];  if( tmp.val < ymin.val )	 ymin = tmp;	 smem[tidx] = ymin;
	  if( Ndata > 16 ){
	    tmp = smem[tidx ^ 16];  if( tmp.val < ymin.val )	 ymin = tmp;	 smem[tidx] = ymin;
	  }}}}}
	  //---------------------------------------------------------------
	}/* if( tidx < Ndata ){ */
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	if( (hidx == tidx) && (tidx < bnum) )
	  smem[(tidx >> 5)] = zmin;/* := (tidx / warpSize) */
	__syncthreads();
	//-----------------------------------------------------------------
	/* bnum <= 1024 --->>> Ndata <= 32 */
	if( tidx < Ndata ){
	  //---------------------------------------------------------------
	  zmin = smem[tidx];
	  if( Ndata >  1 ){
	    tmp = smem[tidx ^  1];  if( tmp.val < zmin.val )	 zmin = tmp;	 smem[tidx] = zmin;
	  if( Ndata >  2 ){
	    tmp = smem[tidx ^  2];  if( tmp.val < zmin.val )	 zmin = tmp;	 smem[tidx] = zmin;
	  if( Ndata >  4 ){
	    tmp = smem[tidx ^  4];  if( tmp.val < zmin.val )	 zmin = tmp;	 smem[tidx] = zmin;
	  if( Ndata >  8 ){
	    tmp = smem[tidx ^  8];  if( tmp.val < zmin.val )	 zmin = tmp;	 smem[tidx] = zmin;
	  if( Ndata > 16 ){
	    tmp = smem[tidx ^ 16];  if( tmp.val < zmin.val )	 zmin = tmp;	 smem[tidx] = zmin;
	  }}}}}
	  //---------------------------------------------------------------
	}/* if( tidx < Ndata ){ */
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------
      }/* if( bnum > warpSize ){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      if( tidx == 0 ){
	min_all[0] = xmin;
	min_all[1] = ymin;
	min_all[2] = zmin;
      }/* if( tidx == 0 ){ */
      //-------------------------------------------------------------------
    }/* if( bidx == 0 ){ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* get maximum */
    //---------------------------------------------------------------------
    /* bnum is set to (# of SMs) * NBLOCKS_PER_SM_EB <= NTHREADS_EB */
    /* ---> # of working block is 1 (bidx = bnum - 1) */
    if( bidx == (bnum - 1) ){
      //-------------------------------------------------------------------
      smem[tidx] = xmax;
      //-------------------------------------------------------------------
      if( tidx < bnum  ){
	//-----------------------------------------------------------------
	/* reduction within a warp */
	//-----------------------------------------------------------------
	xmax = max_all[tidx];
	smem[tidx] = xmax;
	if( bnum >  1 ){
	  tmp = smem[tidx ^  1];  if( tmp.val > xmax.val )	   xmax = tmp;	smem[tidx] = xmax;
	if( bnum >  2 ){
	  tmp = smem[tidx ^  2];  if( tmp.val > xmax.val )	   xmax = tmp;	smem[tidx] = xmax;
	if( bnum >  4 ){
	  tmp = smem[tidx ^  4];  if( tmp.val > xmax.val )	   xmax = tmp;	smem[tidx] = xmax;
	if( bnum >  8 ){
	  tmp = smem[tidx ^  8];  if( tmp.val > xmax.val )	   xmax = tmp;	smem[tidx] = xmax;
	if( bnum > 16 ){
	  tmp = smem[tidx ^ 16];  if( tmp.val > xmax.val )	   xmax = tmp;	smem[tidx] = xmax;
	}}}}}
	//-----------------------------------------------------------------
	ymax = max_all[tidx + bnum];
	smem[tidx] = ymax;
	if( bnum >  1 ){
	  tmp = smem[tidx ^  1];  if( tmp.val > ymax.val )	   ymax = tmp;	smem[tidx] = ymax;
	if( bnum >  2 ){
	  tmp = smem[tidx ^  2];  if( tmp.val > ymax.val )	   ymax = tmp;	smem[tidx] = ymax;
	if( bnum >  4 ){
	  tmp = smem[tidx ^  4];  if( tmp.val > ymax.val )	   ymax = tmp;	smem[tidx] = ymax;
	if( bnum >  8 ){
	  tmp = smem[tidx ^  8];  if( tmp.val > ymax.val )	   ymax = tmp;	smem[tidx] = ymax;
	if( bnum > 16 ){
	  tmp = smem[tidx ^ 16];  if( tmp.val > ymax.val )	   ymax = tmp;	smem[tidx] = ymax;
	}}}}}
	//-----------------------------------------------------------------
	zmax = max_all[tidx + 2 * bnum];
	smem[tidx] = zmax;
	if( bnum >  1 ){
	  tmp = smem[tidx ^  1];  if( tmp.val > zmax.val )	   zmax = tmp;	smem[tidx] = zmax;
	if( bnum >  2 ){
	  tmp = smem[tidx ^  2];  if( tmp.val > zmax.val )	   zmax = tmp;	smem[tidx] = zmax;
	if( bnum >  4 ){
	  tmp = smem[tidx ^  4];  if( tmp.val > zmax.val )	   zmax = tmp;	smem[tidx] = zmax;
	if( bnum >  8 ){
	  tmp = smem[tidx ^  8];  if( tmp.val > zmax.val )	   zmax = tmp;	smem[tidx] = zmax;
	if( bnum > 16 ){
	  tmp = smem[tidx ^ 16];  if( tmp.val > zmax.val )	   zmax = tmp;	smem[tidx] = zmax;
	}}}}}
	//-----------------------------------------------------------------
#if 0
	if( tidx == hidx )
	  printf("%d-th warp (max):\t% e, % e, % e\n", tidx >> 5, xmax.val, ymax.val, zmax.val);
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
	  smem[(tidx >> 5)] = xmax;/* := (tidx / warpSize) */
	__syncthreads();
	//-----------------------------------------------------------------
	/* bnum <= 1024 --->>> Ndata <= 32 */
	if( tidx < Ndata ){
	  //---------------------------------------------------------------
	  xmax = smem[tidx];
	  if( Ndata >  1 ){
	    tmp = smem[tidx ^  1];  if( tmp.val > xmax.val )	xmax = tmp;	smem[tidx] = xmax;
	  if( Ndata >  2 ){
	    tmp = smem[tidx ^  2];  if( tmp.val > xmax.val )	xmax = tmp;	smem[tidx] = xmax;
	  if( Ndata >  4 ){
	    tmp = smem[tidx ^  4];  if( tmp.val > xmax.val )	xmax = tmp;	smem[tidx] = xmax;
	  if( Ndata >  8 ){
	    tmp = smem[tidx ^  8];  if( tmp.val > xmax.val )	xmax = tmp;	smem[tidx] = xmax;
	  if( Ndata > 16 ){
	    tmp = smem[tidx ^ 16];  if( tmp.val > xmax.val )	xmax = tmp;	smem[tidx] = xmax;
	  }}}}}
	  //---------------------------------------------------------------
	}/* if( tidx < Ndata ){ */
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	if( (hidx == tidx) && (tidx < bnum) )
	  smem[(tidx >> 5)] = ymax;/* := (tidx / warpSize) */
	__syncthreads();
	//-----------------------------------------------------------------
	/* bnum <= 1024 --->>> Ndata <= 32 */
	if( tidx < Ndata ){
	  //---------------------------------------------------------------
	  ymax = smem[tidx];
	  if( Ndata >  1 ){
	    tmp = smem[tidx ^  1];  if( tmp.val > ymax.val )	ymax = tmp;	smem[tidx] = ymax;
	  if( Ndata >  2 ){
	    tmp = smem[tidx ^  2];  if( tmp.val > ymax.val )	ymax = tmp;	smem[tidx] = ymax;
	  if( Ndata >  4 ){
	    tmp = smem[tidx ^  4];  if( tmp.val > ymax.val )	ymax = tmp;	smem[tidx] = ymax;
	  if( Ndata >  8 ){
	    tmp = smem[tidx ^  8];  if( tmp.val > ymax.val )	ymax = tmp;	smem[tidx] = ymax;
	  if( Ndata > 16 ){
	    tmp = smem[tidx ^ 16];  if( tmp.val > ymax.val )	ymax = tmp;	smem[tidx] = ymax;
	  }}}}}
	  //---------------------------------------------------------------
	}/* if( tidx < Ndata ){ */
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	if( (hidx == tidx) && (tidx < bnum) )
	  smem[(tidx >> 5)] = zmax;/* := (tidx / warpSize) */
	__syncthreads();
	//-----------------------------------------------------------------
	/* bnum <= 1024 --->>> Ndata <= 32 */
	if( tidx < Ndata ){
	  //---------------------------------------------------------------
	  zmax = smem[tidx];
	  if( Ndata >  1 ){
	    tmp = smem[tidx ^  1];  if( tmp.val > zmax.val )	zmax = tmp;	smem[tidx] = zmax;
	  if( Ndata >  2 ){
	    tmp = smem[tidx ^  2];  if( tmp.val > zmax.val )	zmax = tmp;	smem[tidx] = zmax;
	  if( Ndata >  4 ){
	    tmp = smem[tidx ^  4];  if( tmp.val > zmax.val )	zmax = tmp;	smem[tidx] = zmax;
	  if( Ndata >  8 ){
	    tmp = smem[tidx ^  8];  if( tmp.val > zmax.val )	zmax = tmp;	smem[tidx] = zmax;
	  if( Ndata > 16 ){
	    tmp = smem[tidx ^ 16];  if( tmp.val > zmax.val )	zmax = tmp;	smem[tidx] = zmax;
	  }}}}}
	  //---------------------------------------------------------------
	}/* if( tidx < Ndata ){ */
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------
      }/* if( bnum > warpSize ){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      if( tidx == 0 ){
	max_all[0] = xmax;
	max_all[1] = ymax;
	max_all[2] = zmax;
      }/* if( tidx == 0 ){ */
      //-------------------------------------------------------------------
    }/* if( bidx == (bnum - 1) ){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------
    xmin = min_all[0];    ymin = min_all[1];    zmin = min_all[2];
    xmax = max_all[0];    ymax = max_all[1];    zmax = max_all[2];
    //---------------------------------------------------------------------
  }/* if( bnum > 1 ){ */
  /* else{ */
  /*   //--------------------------------------------------------------------- */
  /*   min = smem[               0]; */
  /*   max = smem[NTHREADS_EB >> 1]; */
  /*   //--------------------------------------------------------------------- */
  /* }/\* else{ *\/ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------







//-------------------------------------------------------------------------
#ifdef  BUILD_LET_ON_DEVICE
//-------------------------------------------------------------------------
/* set pseudo i-particle based on geometrical center of all i-particles */
//-------------------------------------------------------------------------
extern "C"
void setPseudoIparticle4LET
  (const cudaStream_t stream, domainInfo *let, MPIcfg_tree mpi,
   const soaTreeNode tree, const int numSendGuess, const soaTreeWalkBuf buf
   )
{
  //-----------------------------------------------------------------------
  makeLET_kernel<<<1, NTHREADS_MAKE_LET, SMEM_SIZE, stream>>>
    (let->icom,
#ifdef  GADGET_MAC
     let->amin,
#endif//GADGET_MAC
     &let->numSend_dev,
     tree.more, tree.jpos, tree.mj,
     &(tree.more[let->headSend]), &(tree.jpos[let->headSend]), &(tree.mj[let->headSend]),
#ifndef USE_SMID_TO_GET_BUFID
#ifndef TRY_MODE_ABOUT_BUFFER
     buf.active,
#endif//TRY_MODE_ABOUT_BUFFER
     buf.freeNum,
#endif//USE_SMID_TO_GET_BUFID
     buf.freeLst, buf.buffer, NGROUPS * buf.bufSize, buf.fail
#ifdef  MONITOR_LETGEN_TIME
     , cycles
#endif//MONITOR_LETGEN_TIME
     );
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaMemcpyAsync(&let->numSend, &let->numSend_dev, sizeof(int), cudaMemcpyDeviceToHost, stream));
  if( let->numSend > numSendGuess ){
    __KILL__(stderr, "ERROR: predicted size of send buffer(%d) is not sufficient for true size of that(%d) @ rank %d for rand %d.\n\tsuggestion: consider increasing \"LETSIZE_OVERESTIMATION_FACTOR\" defined in src/tree/let.h (current value is %f).\n", numSendGuess, let->numSend, mpi.rank, let->rank, LETSIZE_OVERESTIMATION_FACTOR);
  }/* if( *numSend_hst > numSendGuess ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//BUILD_LET_ON_DEVICE
//-------------------------------------------------------------------------
