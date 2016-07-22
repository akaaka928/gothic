/*************************************************************************\
 *                                                                       *
                  last updated on 2015/11/16(Mon) 12:15:33
 *                                                                       *
 *    Neighbor searching using breadth-first tree                        *
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
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/gsync_dev.cu"
//-------------------------------------------------------------------------
#include "make.h"
#include "make_dev.h"/* to read NBLOCKS_PER_SM_MAC */
#include "neighbor_dev.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef FACILE_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------
/* radix sort library (must be called after neighbor_dev.h, because settings are written in neighbor_dev.h) */
#include "../sort/radix_dev.h"
#include "../sort/radix_inc.cu"
#define SORT_ONLY
#include "../sort/radix_inc.cu"
#undef  SORT_ONLY
//-------------------------------------------------------------------------
/* initialized by flipped FLT_MAX */
/* FLT_MAX = 3.402823e+38 , 0x7f7fffff         --> 0xff7fffff         (flipped) */
/* DBL_MAX = 1.797693e+308, 0x7fefffffffffffff --> 0xffefffffffffffff (flipped) */
#ifdef  DOUBLE_PRECISION
#define FLIPPED_REAL_MAX (0xffefffffffffffff)
#else///DOUBLE_PRECISION
#define FLIPPED_REAL_MAX (0xff7fffff)
#endif//DOUBLE_PRECISION
//-------------------------------------------------------------------------
#   if  NBUF_NEIGHBOR == 4
typedef union __align__(16)
{
  int4 i4;
  int  ia[4];
} alignedInt;
typedef union __align__(16)
{
  real4 r4;
  real  ra[4];
} alignedFlt;
#endif//NBUF_NEIGHBOR == 4
//-------------------------------------------------------------------------
#   if  NBUF_NEIGHBOR == 2
typedef union __align__(8)
{
  int2 i4;
  int  ia[2];
} alignedInt;
typedef union __align__(8)
{
  real2 r4;
  real  ra[2];
} alignedFlt;
#endif//NBUF_NEIGHBOR == 2
//-------------------------------------------------------------------------
typedef union
{
  int  i;
  real r;
} int_real;
//-------------------------------------------------------------------------
#endif//FACILE_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  FACILE_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------
__global__ void facileNeighborSearching_kernel(const int Ni, READ_ONLY position * RESTRICT ibody, real * RESTRICT neighbor_length)
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  __shared__ position ipos[NTHREADS_FACILE_NS + NEIGHBOR_NUM];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  if( gidx < Ni ){
    //---------------------------------------------------------------------
    /* load position of i-particles and neighbor candidates */
    //---------------------------------------------------------------------
    /* int idx = gidx - (NEIGHBOR_NUM >> 1);    if( idx < 0 )      idx = 0; */
    /* ipos[tidx] = ibody[idx]; */
    ipos[tidx] = ibody[gidx];
    if( tidx < NEIGHBOR_NUM ){
      //-------------------------------------------------------------------
      /* idx = gidx - (NEIGHBOR_NUM >> 1) + NTHREADS_FACILE_NS; */
      /* if( idx > (Ni - 1) ) */
      /* 	idx = Ni - 1; */
      /* ipos[NTHREADS_FACILE_NS + tidx] = ibody[idx]; */
      int idx = gidx + NTHREADS_FACILE_NS;      if( idx > (Ni - 1) )      	idx = Ni - 1;
      ipos[NTHREADS_FACILE_NS + tidx] = ibody[idx];
      //-------------------------------------------------------------------
    }/* if( tidx < NEIGHBOR_NUM ){ */
    __syncthreads();
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* calculate distance with NEIGHBOR_NUM particles and remember the maximum */
    //---------------------------------------------------------------------
    /* idx = (NEIGHBOR_NUM >> 1) + tidx; */
    /* const position pi = ipos[idx]; */
    const position pi = ipos[tidx];
    real r2max = ZERO;
    //---------------------------------------------------------------------
/*     /\* this part may be unnecessary *\/ */
/* #pragma unroll */
/*     for(int ii = 0; ii < (NEIGHBOR_NUM >> 1); ii++){ */
/*       //------------------------------------------------------------------- */
/*       const position pj = ipos[idx - ii]; */
/*       //------------------------------------------------------------------- */
/*       const real dx = pj.x - pi.x; */
/*       const real dy = pj.y - pi.y; */
/*       const real dz = pj.z - pi.z; */
/*       const real r2 = 1.0e-30f + dx * dx + dy * dy + dz * dz; */
/*       //------------------------------------------------------------------- */
/*       if( r2 > r2max ) */
/* 	r2max = r2; */
/*       //------------------------------------------------------------------- */
/*     }/\* for(int ii = 0; ii < (NEIGHBOR_NUM >> 1); ii++){ *\/ */
    //---------------------------------------------------------------------
/* #pragma unroll */
/*     for(int ii = 0; ii < (NEIGHBOR_NUM >> 1); ii++){ */
#pragma unroll
    for(int ii = 0; ii < NEIGHBOR_NUM; ii++){
      //-------------------------------------------------------------------
      const position pj = ipos[tidx + 1 + ii];
      //-------------------------------------------------------------------
      const real dx = pj.x - pi.x;
      const real dy = pj.y - pi.y;
      const real dz = pj.z - pi.z;
      const real r2 = 1.0e-30f + dx * dx + dy * dy + dz * dz;
      //-------------------------------------------------------------------
      if( r2 > r2max )
	r2max = r2;
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NEIGHBOR_NUM; ii++){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* store the derived guess about the length of neighbor arm */
    //---------------------------------------------------------------------
    neighbor_length[gidx] = r2max * RSQRT(r2max);
    //---------------------------------------------------------------------
  }/* if( gidx < Ni ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
void facileNeighborSearching_dev(const int Ni, const iparticle pi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  int Nrem = BLOCKSIZE(Ni, NTHREADS_FACILE_NS);
  //-----------------------------------------------------------------------
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    facileNeighborSearching_kernel<<<Nrem, NTHREADS_FACILE_NS>>>(Ni, pi.pos, pi.neighbor);
  //-----------------------------------------------------------------------
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
      int Nsub = Nblck * NTHREADS_FACILE_NS;
      facileNeighborSearching_kernel<<<Nblck, NTHREADS_FACILE_NS>>>(Nsub, &pi.pos[hidx], &pi.neighbor[hidx]);
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
  getLastCudaError("facileNeighborSearching_kernel");
  //-----------------------------------------------------------------------
  /* facileNeighborSearching_kernel<<<BLOCKSIZE(Ni, NTHREADS_FACILE_NS), NTHREADS_FACILE_NS>>>(Ni, pi.pos, pi.neighbor); */
  /* getLastCudaError("facileNeighborSearching_kernel"); */
  /* /\* cudaDeviceSynchronize(); *\/ */
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///FACILE_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------
/* parallel prefix sum within a group of TSUB threads (TSUB <= 32 to use implicit synchronization) */
/* type of prefix sum is inclusive */
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
__device__ __forceinline__  int prefixSumTsub(const int psum,                                           const int lane)
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
__device__ __forceinline__ void prefixSumTsub(const int psum, volatile int_real * smem, const int tidx, const int lane)
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  //-----------------------------------------------------------------------
  int val = psum;
  int tmp;
#   if  TSUB_NEIGHBOR >=  2
  tmp = __shfl_up(val,  1, TSUB_NEIGHBOR);  if( lane >=  1 )    val += tmp;
#   if  TSUB_NEIGHBOR >=  4
  tmp = __shfl_up(val,  2, TSUB_NEIGHBOR);  if( lane >=  2 )    val += tmp;
#   if  TSUB_NEIGHBOR >=  8
  tmp = __shfl_up(val,  4, TSUB_NEIGHBOR);  if( lane >=  4 )    val += tmp;
#   if  TSUB_NEIGHBOR >= 16
  tmp = __shfl_up(val,  8, TSUB_NEIGHBOR);  if( lane >=  8 )    val += tmp;
#   if  TSUB_NEIGHBOR == 32
  tmp = __shfl_up(val, 16, TSUB_NEIGHBOR);  if( lane >= 16 )    val += tmp;
#endif//TSUB_NEIGHBOR == 32
#endif//TSUB_NEIGHBOR >= 16
#endif//TSUB_NEIGHBOR >=  8
#endif//TSUB_NEIGHBOR >=  4
#endif//TSUB_NEIGHBOR >=  2
  return (val);
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  //-----------------------------------------------------------------------
  smem[tidx].i = psum;
  //-----------------------------------------------------------------------
# if  TSUB_NEIGHBOR >=  2
  if( lane >=  1 )    smem[tidx].i += smem[tidx -  1].i;
# if  TSUB_NEIGHBOR >=  4
  if( lane >=  2 )    smem[tidx].i += smem[tidx -  2].i;
# if  TSUB_NEIGHBOR >=  8
  if( lane >=  4 )    smem[tidx].i += smem[tidx -  4].i;
# if  TSUB_NEIGHBOR >= 16
  if( lane >=  8 )    smem[tidx].i += smem[tidx -  8].i;
# if  TSUB_NEIGHBOR == 32
  if( lane >= 16 )    smem[tidx].i += smem[tidx - 16].i;
#endif//TSUB_NEIGHBOR == 32
#endif//TSUB_NEIGHBOR >= 16
#endif//TSUB_NEIGHBOR >=  8
#endif//TSUB_NEIGHBOR >=  4
#endif//TSUB_NEIGHBOR >=  2
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
__device__ __forceinline__ real accumulateRealTsub(const real  sum)
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
__device__ __forceinline__ void accumulateRealTsub(      real *sum, volatile int_real * smem, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  //-----------------------------------------------------------------------
  real val = sum;
#   if  TSUB_NEIGHBOR >=  2
  real tmp;
  tmp = __shfl_xor(val,  1, TSUB_NEIGHBOR);  val += tmp;
#   if  TSUB_NEIGHBOR >=  4
  tmp = __shfl_xor(val,  2, TSUB_NEIGHBOR);  val += tmp;
#   if  TSUB_NEIGHBOR >=  8
  tmp = __shfl_xor(val,  4, TSUB_NEIGHBOR);  val += tmp;
#   if  TSUB_NEIGHBOR >= 16
  tmp = __shfl_xor(val,  8, TSUB_NEIGHBOR);  val += tmp;
#   if  TSUB_NEIGHBOR == 32
  tmp = __shfl_xor(val, 16, TSUB_NEIGHBOR);  val += tmp;
#endif//TSUB_NEIGHBOR == 32
#endif//TSUB_NEIGHBOR >= 16
#endif//TSUB_NEIGHBOR >=  8
#endif//TSUB_NEIGHBOR >=  4
#endif//TSUB_NEIGHBOR >=  2
  return (__shfl(val, 0, TSUB_NEIGHBOR));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  //-----------------------------------------------------------------------
  smem[tidx].r = *sum;
  //-----------------------------------------------------------------------
#   if  TSUB_NEIGHBOR >=  2
  real tmp;
  tmp = smem[tidx ^  1].r;  *sum += tmp;  smem[tidx].r = *sum;
#   if  TSUB_NEIGHBOR >=  4
  tmp = smem[tidx ^  2].r;  *sum += tmp;  smem[tidx].r = *sum;
#   if  TSUB_NEIGHBOR >=  8
  tmp = smem[tidx ^  4].r;  *sum += tmp;  smem[tidx].r = *sum;
#   if  TSUB_NEIGHBOR >= 16
  tmp = smem[tidx ^  8].r;  *sum += tmp;  smem[tidx].r = *sum;
#   if  TSUB_NEIGHBOR == 32
  tmp = smem[tidx ^ 16].r;  *sum += tmp;  smem[tidx].r = *sum;
#endif//TSUB_NEIGHBOR == 32
#endif//TSUB_NEIGHBOR >= 16
#endif//TSUB_NEIGHBOR >=  8
#endif//TSUB_NEIGHBOR >=  4
#endif//TSUB_NEIGHBOR >=  2
  //-----------------------------------------------------------------------
  *sum = smem[head].r;
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
__device__ __forceinline__  int accumulateIntTsub(const int  sum)
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
__device__ __forceinline__ void accumulateIntTsub(      int *sum, volatile int_real * smem, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  //-----------------------------------------------------------------------
  int val = sum;
#   if  TSUB_NEIGHBOR >=  2
  int tmp;
  tmp = __shfl_xor(val,  1, TSUB_NEIGHBOR);  val += tmp;
#   if  TSUB_NEIGHBOR >=  4
  tmp = __shfl_xor(val,  2, TSUB_NEIGHBOR);  val += tmp;
#   if  TSUB_NEIGHBOR >=  8
  tmp = __shfl_xor(val,  4, TSUB_NEIGHBOR);  val += tmp;
#   if  TSUB_NEIGHBOR >= 16
  tmp = __shfl_xor(val,  8, TSUB_NEIGHBOR);  val += tmp;
#   if  TSUB_NEIGHBOR == 32
  tmp = __shfl_xor(val, 16, TSUB_NEIGHBOR);  val += tmp;
#endif//TSUB_NEIGHBOR == 32
#endif//TSUB_NEIGHBOR >= 16
#endif//TSUB_NEIGHBOR >=  8
#endif//TSUB_NEIGHBOR >=  4
#endif//TSUB_NEIGHBOR >=  2
  return (__shfl(val, 0, TSUB_NEIGHBOR));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  //-----------------------------------------------------------------------
  smem[tidx].i = *sum;
  //-----------------------------------------------------------------------
#   if  TSUB_NEIGHBOR >=  2
  int tmp;
  tmp = smem[tidx ^  1].i;  *sum += tmp;  smem[tidx].i = *sum;
#   if  TSUB_NEIGHBOR >=  4
  tmp = smem[tidx ^  2].i;  *sum += tmp;  smem[tidx].i = *sum;
#   if  TSUB_NEIGHBOR >=  8
  tmp = smem[tidx ^  4].i;  *sum += tmp;  smem[tidx].i = *sum;
#   if  TSUB_NEIGHBOR >= 16
  tmp = smem[tidx ^  8].i;  *sum += tmp;  smem[tidx].i = *sum;
#   if  TSUB_NEIGHBOR == 32
  tmp = smem[tidx ^ 16].i;  *sum += tmp;  smem[tidx].i = *sum;
#endif//TSUB_NEIGHBOR == 32
#endif//TSUB_NEIGHBOR >= 16
#endif//TSUB_NEIGHBOR >=  8
#endif//TSUB_NEIGHBOR >=  4
#endif//TSUB_NEIGHBOR >=  2
  //-----------------------------------------------------------------------
  *sum = smem[head].i;
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  //-----------------------------------------------------------------------
}
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
__device__ __forceinline__ int occupyBuffer(const int tidx, uint *freeLst, uint *bufIdx)
{
  //-----------------------------------------------------------------------
  int target = 0;
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
#   if  NBLOCKS_PER_SM_NEIGHBOR == 2
    target         = getSMidx() * NBLOCKS_PER_SM_NEIGHBOR;
#else///NBLOCKS_PER_SM_NEIGHBOR == 2
    const int head = getSMidx() * NBLOCKS_PER_SM_NEIGHBOR;
#endif//NBLOCKS_PER_SM_NEIGHBOR == 2
    //---------------------------------------------------------------------
    uint tmp = UINT_MAX;
    while( tmp == UINT_MAX ){
      //-------------------------------------------------------------------
#   if  NBLOCKS_PER_SM_NEIGHBOR == 2
      //-------------------------------------------------------------------
      target ^= 1;
      tmp = atomicExch(&freeLst[target], UINT_MAX);
      //-------------------------------------------------------------------
#else///NBLOCKS_PER_SM_NEIGHBOR == 2
      //-------------------------------------------------------------------
      for(int ii = 0; ii < NBLOCKS_PER_SM_NEIGHBOR; ii++){
	//-----------------------------------------------------------------
	tmp = atomicExch(&freeLst[head + ii], UINT_MAX);
	//-----------------------------------------------------------------
	if( tmp != UINT_MAX ){
	  target = head + ii;
	  break;
	}/* if( tmp != UINT_MAX ){ */
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < NBLOCKS_PER_SM_NEIGHBOR; ii++){ */
      //-------------------------------------------------------------------
#endif//NBLOCKS_PER_SM_NEIGHBOR == 2
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
  if( tidx == 0 )
    atomicExch(&freeLst[target], bufIdx);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///USE_SMID_TO_GET_BUFID
//-------------------------------------------------------------------------
#ifdef  TRY_MODE_ABOUT_BUFFER
//-------------------------------------------------------------------------
__device__ __forceinline__ int occupyBuffer(const int tidx, const int bidx, const int bufNum, uint *freeLst, uint *bufIdx)
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
  if( tidx == 0 )
    atomicExch(&freeLst[target], bufIdx);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///TRY_MODE_ABOUT_BUFFER
//-------------------------------------------------------------------------
__device__ __forceinline__ void  occupyBuffer(const int tidx, uint *freeNum, uint *freeLst, uint *bufIdx, int *active)
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
__device__ __forceinline__ void releaseBuffer(const int tidx, uint *freeNum, uint *freeLst, const int bufIdx, int *active)
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


//-------------------------------------------------------------------------
/* calculate distance between an i-particle and the corresponding NEIGHBOR_NUM-th neighbor particle */
//-------------------------------------------------------------------------
/* Ni              :: input          :: Number of i-particles */
/* pi              :: input          :: position and mass of N-body particles */
/* neighbor_length ::         output :: distance between an i-particle and the corresponding NEIGHBOR_NUM-th neighbor particle */
/* cell            :: input          :: head index and number of N-body particles contained in the corresponding tree cell */
/* leaf            :: input          :: a flag to remember the corresponding tree cell is either leaf(true) of node(false) */
/* node            :: input          :: head index and number of pseudo particles of the corresponding tree cell (a leaf cell contains up to NCRIT particles) */
/* more            :: input          :: head index and number of child particles of the corresponding j-particle (a leaf cell contains up to NCRIT particles) */
/* node2cell       :: input          :: index of the tree cell corresponding a pseudo particle */
/* pj              :: input          :: position and squared radius of pseudo N-body particle as j-particles */
/* bmax            :: input          :: size of pseudo N-body particle as j-particles */
//-------------------------------------------------------------------------
/* NOTE for future upgrade: if this kernel is too slow; then customize to employ radix sort with Nvec = 4 (the optimal one) */
//-------------------------------------------------------------------------
__global__ void __launch_bounds__(NTHREADS_NEIGHBOR, NBLOCKS_PER_SM_NEIGHBOR) searchNeighbors_kernel
(const int Ni, READ_ONLY position * RESTRICT pi, real * RESTRICT neighbor_length,
 READ_ONLY treecell * RESTRICT cell, READ_ONLY bool * RESTRICT leaf, READ_ONLY uint * RESTRICT node, READ_ONLY uint * RESTRICT more, READ_ONLY int * RESTRICT node2cell, READ_ONLY int * RESTRICT niSub,
 READ_ONLY jparticle * RESTRICT pj, READ_ONLY real * RESTRICT bmax,
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
 int * RESTRICT active, uint * RESTRICT freeNum,
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) &&  defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
 const int freeNum,
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) &&  defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
 uint * RESTRICT freeLst,
 volatile int * RESTRICT more0Buf, volatile int * RESTRICT more1Buf, volatile real * RESTRICT rjminBuf,
 int * RESTRICT overflow, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  /* const int bidx =  BLOCKIDX_X1D; */
  /* const int bnum =   GRIDDIM_X1D; */
  //-----------------------------------------------------------------------
  const int lane = tidx & (TSUB_NEIGHBOR - 1);
  const int head = tidx - lane;
#ifndef USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  const int tail = head + TSUB_NEIGHBOR - 1;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  const int hbuf = (head / TSUB_NEIGHBOR) * TSUB_NEIGHBOR * NBUF_NEIGHBOR;/* head index of the shared array close and queue within a thread group */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  int smem;
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  __shared__  int_real  smem[NTHREADS_NEIGHBOR];
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
  /* __shared__ jparticle pj_sm[NTHREADS_NEIGHBOR]; */
  __shared__      real rjbuf[NTHREADS_NEIGHBOR * NBUF_NEIGHBOR];
  __shared__  int      list0[NTHREADS_NEIGHBOR * NBUF_NEIGHBOR];
  __shared__  int      list1[NTHREADS_NEIGHBOR * NBUF_NEIGHBOR];
#ifndef REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
  __shared__  int      pjidx[NTHREADS_NEIGHBOR * NBUF_NEIGHBOR];
#endif//REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
  //-----------------------------------------------------------------------
  /* arrays and variables for radix sorting */
  __shared__ uint sort_tmp[NTHREADS_NEIGHBOR * SORT_ELEMENTS_PER_THREAD];
  __shared__  int sort_sub[NTHREADS_NEIGHBOR * SORT_ELEMENTS_PER_THREAD];/* contains number of child particles */
  __shared__ uint sort_dat[NEIGHBOR_NUM_INC * NGROUPS_NEIGHBOR];
  __shared__  int sort_num[NEIGHBOR_NUM_INC * NGROUPS_NEIGHBOR];/* contains number of child particles */
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  __shared__ uint sort_smem[NTHREADS_NEIGHBOR * SORT_ELEMENTS_PER_THREAD];
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
  __shared__ uint4_array sort_sbuf[NTHREADS_SORT];
#endif//RADIX_SORT_CHECK_BITS > 2
  const int hp_sort = head * SORT_ELEMENTS_PER_THREAD;
  const int hp_rmin = (head / TSUB_NEIGHBOR) * NEIGHBOR_NUM_INC;
  union {float f; uint u;} flip;
  //-----------------------------------------------------------------------
  /* head index of remote buffers */
#ifdef  USE_SMID_TO_GET_BUFID_NEIGHBOR
  const int target = occupyBuffer(tidx, freeLst, sort_tmp);
#else///USE_SMID_TO_GET_BUFID_NEIGHBOR
#ifdef  TRY_MODE_ABOUT_BUFFER_NEIGHBOR
  const int target = occupyBuffer(tidx, BLOCKIDX_X1D, freeNum, freeLst, sort_tmp);
#else///TRY_MODE_ABOUT_BUFFER_NEIGHBOR
  occupyBuffer(tidx, freeNum, freeLst, sort_tmp, active);
#endif//TRY_MODE_ABOUT_BUFFER_NEIGHBOR
#endif//USE_SMID_TO_GET_BUFID_NEIGHBOR
  const int bufIdx = (int)sort_tmp[0];
  __syncthreads();
  const int bufHead = ((head / TSUB_NEIGHBOR) + bufIdx * NGROUPS_NEIGHBOR) * NUM_ALLOC_NEIGHBORBUF;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate neighbor length of all i-particles */
  //-----------------------------------------------------------------------
  /* /\* loop to set a maximum number for # of blocks *\/ */
  /* for(int ii = 0; ii < bnum * BLOCKSIZE(Ni, bnum * NGROUPS_NEIGHBOR); ii += bnum){ */
    //---------------------------------------------------------------------
    /* load i-particle */
    const      int iidx = gidx / TSUB_NEIGHBOR;
    const position ipos = (iidx < Ni) ? (pi[iidx]) : (pi[Ni - 1]);
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* estimate size of neighbor particle distribution */
    //---------------------------------------------------------------------
    /* load the root cell */
    const treecell root = cell[0];
    /* initialize list of examined tree nodes and related variables */
    int inum = root.num;
    int Ntry = 1;
    if( lane == 0 )	    list0[hbuf] = 0;
    //---------------------------------------------------------------------
    bool not_yet = true;
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* pick up NI_NEIGHBOR_ESTIMATE i-particles in maximum to estimate bmax */
    while( inum > NI_NEIGHBOR_ESTIMATE ){
      //-------------------------------------------------------------------
#if 1
      real rmax = REAL_MAX;
#else
      real rmax = 10.0f;
#endif
      //-------------------------------------------------------------------
      int Nmin = 0;
#pragma unroll
      for(int jj = lane; jj < NEIGHBOR_NUM_INC; jj += TSUB_SORT){
	sort_dat[hp_rmin + jj] = FLIPPED_REAL_MAX;
	sort_num[hp_rmin + jj] = 0;
      }
      //-------------------------------------------------------------------
      int Nloc = 0;
      int Nbuf = 0;
      //-------------------------------------------------------------------
      /* leaf node detector */
      int lowest = 1;
      //-------------------------------------------------------------------
      int Niter = BLOCKSIZE(Ntry, TSUB_NEIGHBOR * NBUF_NEIGHBOR);
      int chead;
      for(int iter = 0; iter < Niter; iter++){
	//-----------------------------------------------------------------
	const int Nsweep = (Ntry > (TSUB_NEIGHBOR * NBUF_NEIGHBOR)) ? (TSUB_NEIGHBOR * NBUF_NEIGHBOR) : Ntry;
	const int ibuf_loop = BLOCKSIZE(Nsweep, TSUB_NEIGHBOR);
	for(int ibuf = 0; ibuf < ibuf_loop; ibuf++){
	  //---------------------------------------------------------------
	  int cnum = 0;
	  if( (lane + ibuf * TSUB_NEIGHBOR) < Nsweep ){
	    //-------------------------------------------------------------
	    /* load a tree node corresponding the tree cell */
	    uint more_tmp = node[list0[hbuf + lane + ibuf * TSUB_NEIGHBOR]];
	    const int nodenum  = 1 + (more_tmp >> IDXBITS);
	    const int nodehead =      more_tmp  & IDXMASK;
	    //-------------------------------------------------------------
	    /* load all child nodes of the tree cell */
	    more_tmp = more[nodehead];
	    cnum  = 1 + (more_tmp >> IDXBITS);
	    chead =      more_tmp  & IDXMASK;
	    for(int jj = 1; jj < nodenum; jj++)
	      cnum += (1 + (more[nodehead + jj] >> IDXBITS));
	    //-------------------------------------------------------------
	    /* check whether the node is a leaf or not */
	    lowest &= (more_tmp == nodehead);
	    //-------------------------------------------------------------
	  }/* if( (lane + ibuf * TSUB_NEIGHBOR) < Nsweep ){ */
	  //---------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  smem = prefixSumTsub(cnum, lane);
	  const int lend = BLOCKSIZE(__shfl(smem, TSUB_NEIGHBOR - 1, TSUB_NEIGHBOR), NBUF_NEIGHBOR * TSUB_NEIGHBOR);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  prefixSumTsub(cnum, smem, tidx, lane);
	  const int lend = BLOCKSIZE(       smem[tail].i,                  NBUF_NEIGHBOR * TSUB_NEIGHBOR);
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  for(int ll = 0; ll < lend; ll++){
	    //-------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    const int unum =
	      (smem         <= (TSUB_NEIGHBOR * NBUF_NEIGHBOR)) ? cnum :
	      ((smem         >= (cnum + TSUB_NEIGHBOR * NBUF_NEIGHBOR)) ? (0) : (cnum + TSUB_NEIGHBOR * NBUF_NEIGHBOR - smem));
	    const int shead = hbuf + smem         - cnum;
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    const int unum =
	      (smem[tidx].i <= (TSUB_NEIGHBOR * NBUF_NEIGHBOR)) ? cnum :
	      ((smem[tidx].i >= (cnum + TSUB_NEIGHBOR * NBUF_NEIGHBOR)) ? (0) : (cnum + TSUB_NEIGHBOR * NBUF_NEIGHBOR - smem[tidx].i));
	    const int shead = hbuf + smem[tidx].i - cnum;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    for(int jj = 0; jj < unum; jj++){
#ifdef  REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
	      list1[shead + jj] = chead;
#else///REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
	      pjidx[shead + jj] = chead;
#endif//REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
	      chead++;
	    }/* for(int jj = 0; jj < unum; jj++){ */
	    cnum -= unum;
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    const int Ntmp = smem         - (NBUF_NEIGHBOR * TSUB_NEIGHBOR);/* Ntmp is a temporal buffer */
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    const int Ntmp = smem[tidx].i - (NBUF_NEIGHBOR * TSUB_NEIGHBOR);/* Ntmp is a temporal buffer */
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    /* pick up candidate tree nodes */
	    //-------------------------------------------------------------
#   if  NBUF_NEIGHBOR == 4
	    alignedFlt rjmax_loc = { REAL_MIN,  REAL_MIN,  REAL_MIN,  REAL_MIN};
	    alignedFlt rjmin_loc = { REAL_MAX,  REAL_MAX,  REAL_MAX,  REAL_MAX};
	    alignedInt pjidx_loc = {NULL_NODE, NULL_NODE, NULL_NODE, NULL_NODE};
	    alignedInt nisub_loc = {        0,         0,         0,         0};
#endif//NBUF_NEIGHBOR == 4
#   if  NBUF_NEIGHBOR == 2
	    alignedFlt rjmax_loc = { REAL_MIN,  REAL_MIN};
	    alignedFlt rjmin_loc = { REAL_MAX,  REAL_MAX};
	    alignedInt pjidx_loc = {NULL_NODE, NULL_NODE};
	    alignedInt nisub_loc = {        0,         0};
#endif//NBUF_NEIGHBOR == 2
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    const int stail = (__shfl(smem, TSUB_NEIGHBOR - 1, TSUB_NEIGHBOR) < (TSUB_NEIGHBOR * NBUF_NEIGHBOR)) ? (__shfl(smem, TSUB_NEIGHBOR - 1, TSUB_NEIGHBOR)) : (TSUB_NEIGHBOR * NBUF_NEIGHBOR);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    const int stail = (       smem[tail].i                  < (TSUB_NEIGHBOR * NBUF_NEIGHBOR)) ? (       smem[tail].i                 ) : (TSUB_NEIGHBOR * NBUF_NEIGHBOR);
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
#if 1
#pragma unroll
	    for(int kk = 0; kk < NBUF_NEIGHBOR; kk++){
	      //-----------------------------------------------------------
	      const int jj = lane + kk * TSUB_NEIGHBOR;
	      if( jj >= stail )		break;
	      /* if( jj > stail )		break; */
	      //-----------------------------------------------------------
#ifdef  REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
	      const int kidx = list1[hbuf + jj];
#else///REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
	      const int kidx = pjidx[hbuf + jj];
#endif//REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
	      const jparticle jpos = pj[kidx];
	      //-----------------------------------------------------------
	      const real dx = jpos.x - ipos.x;
	      const real dy = jpos.y - ipos.y;
	      const real dz = jpos.z - ipos.z;
	      const real d2 = 1.0e-30f + dx * dx + dy * dy + dz * dz;
	      const real dr = d2 * RSQRT(d2);
	      //-----------------------------------------------------------
	      const real rjmax = bmax[kidx] + dr;	      /* the possible maximum distance between the i-particle and particles within the j-node */
	      const real rjmin = -rjmax + (TWO * (UNITY - EPSILON)) * dr;
	      //-----------------------------------------------------------
#if 0
	      if( iidx == 0 )
		printf("kk = %d, lane = %d, jj = %d: kidx = %d (%e, %e, %e): rmax = %e, rjmax = %e, rjmax = %e, d2 = %e, dr = %e\n", kk, lane, jj, kidx, jpos.x, jpos.y, jpos.z, rmax, rjmax, rjmin, d2, dr);
#endif
	      //-----------------------------------------------------------
	      if( rjmin < rmax ){
		//---------------------------------------------------------
		pjidx_loc.ia[kk] = kidx;
		rjmin_loc.ra[kk] = rjmin;
		rjmax_loc.ra[kk] = rjmax;
		//---------------------------------------------------------
		nisub_loc.ia[kk] = niSub[kidx];
		//---------------------------------------------------------
	      }/* if( rjmin < rmax ){ */
	      //-----------------------------------------------------------
	    }/* for(int kk = 0; kk < NBUF_NEIGHBOR; kk++){ */
#else
	    for(int jj = lane; jj < stail; jj += TSUB_NEIGHBOR){
	      //-----------------------------------------------------------
#ifdef  REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
	      const int kidx = list1[hbuf + jj];
#else///REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
	      const int kidx = pjidx[hbuf + jj];
#endif//REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
	      const jparticle jpos = pj[kidx];
	      //-----------------------------------------------------------
	      const real dx = jpos.x - ipos.x;
	      const real dy = jpos.y - ipos.y;
	      const real dz = jpos.z - ipos.z;
	      const real d2 = 1.0e-30f + dx * dx + dy * dy + dz * dz;
	      const real dr = d2 * RSQRT(d2);
	      //-----------------------------------------------------------
	      const real rjmax = bmax[kidx] + dr;	      /* the possible maximum distance between the i-particle and particles within the j-node */
	      const real rjmin = -rjmax + (TWO * (UNITY - EPSILON)) * dr;
	      //-----------------------------------------------------------
	      if( rjmin < rmax ){
		//---------------------------------------------------------
		const int itmp = jj / TSUB_NEIGHBOR;
		pjidx_loc.ia[itmp] = kidx;
		rjmin_loc.ra[itmp] = rjmin;
		rjmax_loc.ra[itmp] = rjmax;
		//---------------------------------------------------------
		/* nisub_loc.ia[itmp] = cell[node2cell[kidx]].num; */
		nisub_loc.ia[itmp] = niSub[kidx];
		//---------------------------------------------------------
	      }/* if( rjmin < rmax ){ */
	      //-----------------------------------------------------------
	    }/* for(int jj = lane; jj < stail; jj += TSUB_NEIGHBOR){ */
#endif
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    /* update rmax by sorting an array storing rjmax */
	    //-------------------------------------------------------------
#pragma unroll
	    for(int jj = lane; jj < TSUB_SORT * SORT_ELEMENTS_PER_THREAD; jj += TSUB_SORT){
	      sort_tmp[hp_sort + jj] = FLIPPED_REAL_MAX;
	      sort_sub[hp_sort + jj] =                0;
	    }
	    //-------------------------------------------------------------
	    int Nsort = 0;
#pragma unroll
	    for(int jj = 0; jj < NBUF_NEIGHBOR; jj++){
	      //-----------------------------------------------------------
	      const int share = ( rjmin_loc.ra[jj] < rmax ) ? 1 : 0;
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      smem = prefixSumTsub(share, lane);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      prefixSumTsub(share, smem, tidx, lane);
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      //-----------------------------------------------------------
	      if( share ){
		//---------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
		const int dst = hp_sort + Nsort + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
		const int dst = hp_sort + Nsort + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
		//---------------------------------------------------------
		flip.f = rjmax_loc.ra[jj];
		sort_tmp[dst] = flip32flt(flip.u);
		sort_sub[dst] = nisub_loc.ia[jj];
		//---------------------------------------------------------
#if 0
		if( iidx == 0 )
		  printf("lane = %d, dst = %d, flip.f = %e\n", lane, dst, flip.f);
#endif
		//---------------------------------------------------------
	      }/* if( share ){ */
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      Nsort += __shfl(smem, TSUB_NEIGHBOR - 1, TSUB_NEIGHBOR);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      Nsort +=        smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      //-----------------------------------------------------------
	    }/* for(int jj = 0; jj < NBUF_NEIGHBOR; jj++){ */
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    /* get NEIGHBOR_NUM-th minimum of rjmax */
	    //-------------------------------------------------------------
	    if( (Nsort + Nmin) > TSUB_SORT * SORT_ELEMENTS_PER_THREAD ){
	      //-----------------------------------------------------------
	      __radixSortTsub32idx(32 / RADIX_SORT_CHECK_BITS, lane, &sort_tmp[hp_sort], &sort_sub[hp_sort]
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
				   , tidx, &sort_smem[hp_sort]
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
				   , &sbuf[head]
#endif//RADIX_SORT_CHECK_BITS > 2
				   );
	      //-----------------------------------------------------------
	      Nsort = NEIGHBOR_NUM_INC;
	      //-----------------------------------------------------------
	    }/* if( (Nsort + Nmin) > TSUB_SORT * SORT_ELEMENTS_PER_THREAD ){ */
	    //-------------------------------------------------------------
#pragma unroll
	    for(int jj = lane; jj < Nmin; jj += TSUB_NEIGHBOR){
	      sort_tmp[hp_sort + Nsort + jj] = sort_dat[hp_rmin + jj];
	      sort_sub[hp_sort + Nsort + jj] = sort_num[hp_rmin + jj];
	    }
	    //-------------------------------------------------------------
	    __radixSortTsub32idx(32 / RADIX_SORT_CHECK_BITS, lane, &sort_tmp[hp_sort], &sort_sub[hp_sort]
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
				 , tidx, &sort_smem[hp_sort]
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
				 , &sbuf[head]
#endif//RADIX_SORT_CHECK_BITS > 2
				 );
	    //-------------------------------------------------------------
	    /* derive Nmin := the minimum number of cells to contain at least NEIGHBOR_NUM + 1 particles */
	    Nmin = ((Nsort + Nmin) < NEIGHBOR_NUM_INC) ? (Nsort + Nmin) : NEIGHBOR_NUM_INC;
	    int Ninc = sort_sub[hp_sort];
	    for(int jj = 1; jj < Nmin; jj++){
	      //-----------------------------------------------------------
	      if( Ninc >= NEIGHBOR_NUM_INC ){
		Nmin = jj;
		break;
	      }/* if( Ninc >= NEIGHBOR_NUM_INC ){ */
	      //-----------------------------------------------------------
	      Ninc += sort_sub[hp_sort + jj];
	      //-----------------------------------------------------------
	    }/* for(int jj = 0; jj < Nmin; jj++){ */
	    //-------------------------------------------------------------
#pragma unroll
	    for(int jj = lane; jj < Nmin; jj += TSUB_NEIGHBOR){
	      sort_dat[hp_rmin + jj] = sort_tmp[hp_sort + jj];
	      sort_num[hp_rmin + jj] = sort_sub[hp_sort + jj];
	    }
	    //-------------------------------------------------------------
	    flip.u = undo32flt(sort_dat[hp_rmin + Nmin - 1]);
	    rmax = flip.f;
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    /* recheck local buffer (is really rjmin smaller than rmax ?) */
	    //-------------------------------------------------------------
#pragma unroll
	    for(int jj = 0; jj < NBUF_NEIGHBOR; jj++){
	      //-----------------------------------------------------------
	      const int share = ( rjmin_loc.ra[jj] < rmax ) ? 1 : 0;
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      smem = prefixSumTsub(share, lane);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      prefixSumTsub(share, smem, tidx, lane);
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      //-----------------------------------------------------------
	      if( share ){
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
		const int dst = hbuf + Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
		const int dst = hbuf + Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
		list1[dst] = pjidx_loc.ia[jj];
		rjbuf[dst] = rjmin_loc.ra[jj];
	      }/* if( share ){ */
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      Nloc += __shfl(smem, TSUB_NEIGHBOR - 1, TSUB_NEIGHBOR);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      Nloc +=        smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	      //-----------------------------------------------------------
#ifndef REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
	      if( Nloc > ((NBUF_NEIGHBOR - 1) * TSUB_NEIGHBOR) )
#endif//REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
		{
		  //-------------------------------------------------------
		  for(int kk = lane; kk < Nloc; kk += TSUB_NEIGHBOR){
		    more1Buf[bufHead + Nbuf + kk] = list1[hbuf + kk];
		    rjminBuf[bufHead + Nbuf + kk] = rjbuf[hbuf + kk];
		  }/* for(int kk = lane; kk < Nloc; kk += TSUB_NEIGHBOR){ */
		  //-------------------------------------------------------
		  Nbuf += Nloc;
		  Nloc = 0;
		  //-------------------------------------------------------
		}/* if( Nloc > ((NBUF_NEIGHBOR - 1) * TSUB_NEIGHBOR) ){ */
	      //-----------------------------------------------------------
	    }/* for(int jj = 0; jj < NBUF_NEIGHBOR; jj++){ */
	    //-------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    smem         = Ntmp;/* Ntmp is a temporal buffer */
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    smem[tidx].i = Ntmp;/* Ntmp is a temporal buffer */
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    //-------------------------------------------------------------
	  }/* for(int ll = 0; ll < lend; ll++){ */
	  //---------------------------------------------------------------
	}/* for(int ibuf = 0; ibuf < NBUF_NEIGHBOR; ibuf++){ */
	//-----------------------------------------------------------------
	Ntry -= Nsweep;
	//-----------------------------------------------------------------
	/* copy data from global memory to shared memory */
	const int Ncopy = (Ntry < (TSUB_NEIGHBOR * NBUF_NEIGHBOR)) ? (Ntry) : (TSUB_NEIGHBOR * NBUF_NEIGHBOR);
	for(int jj = lane; jj < Ncopy; jj += TSUB_NEIGHBOR)
	  list0[hbuf + jj] = more0Buf[(bufHead + (TSUB_NEIGHBOR * NBUF_NEIGHBOR) * (iter + 1)) + jj];
	//-----------------------------------------------------------------
      }/* for(int iter = 0; iter < Niter; iter++){ */
      //-------------------------------------------------------------------
      if( Nbuf != 0 )
	for(int ll = lane; ll < Nloc; ll += TSUB_NEIGHBOR){
	  more1Buf[bufHead + Nbuf + ll] = list1[hbuf + ll];
	  rjminBuf[bufHead + Nbuf + ll] = rjbuf[hbuf + ll];
	}/* for(int ll = lane; ll < Nloc; ll += TSUB_NEIGHBOR){ */
      if( Nbuf != 0 )
	for(int ll = lane; ll < TSUB_NEIGHBOR * NBUF_NEIGHBOR; ll += TSUB_NEIGHBOR){
	  list1[hbuf + ll] = more1Buf[bufHead + ll];
	  rjbuf[hbuf + ll] = rjminBuf[bufHead + ll];
	}/* for(int ll = lane; ll < TSUB_NEIGHBOR * NBUF_NEIGHBOR; ll += TSUB_NEIGHBOR){ */
      //-------------------------------------------------------------------
      Ntry = Nbuf + Nloc;
      if( (lane == 0) && (Ntry > NUM_ALLOC_NEIGHBORBUF) )
	atomicAdd(overflow, 1);
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
      /* if the all candidate pseudo-particles are real-particles (which means more[idx] = idx), then set neighbor length here */
      /* this escape sequence can evade infinite loop (if all nodes reach each leaf level and # of particles within the corresponding cells is greater than NI_NEIGHBOR_ESTIMATE, the while loop never finishes ) */
      //-------------------------------------------------------------------
      /* leaf node detection */
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
      smem = prefixSumTsub(lowest, lane);
      lowest = __shfl(smem, TSUB_NEIGHBOR - 1, TSUB_NEIGHBOR);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
      prefixSumTsub(lowest, smem, tidx, lane);
      lowest =        smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
      //-------------------------------------------------------------------
      /* when the all candidate particles are real particles... */
      if( lowest == TSUB_NEIGHBOR ){
	//-----------------------------------------------------------------
#pragma unroll
	for(int jj = lane; jj < TSUB_SORT * SORT_ELEMENTS_PER_THREAD; jj += TSUB_SORT)
	  sort_tmp[hp_sort + jj] = FLIPPED_REAL_MAX;
	//-----------------------------------------------------------------
	Niter = BLOCKSIZE(Ntry, TSUB_SORT * SORT_ELEMENTS_PER_THREAD - NEIGHBOR_NUM_INC);
	for(int iter = 0; iter < Niter; iter++){
	  //---------------------------------------------------------------
	  const int krem = (Ntry < (TSUB_SORT * SORT_ELEMENTS_PER_THREAD - NEIGHBOR_NUM_INC)) ? Ntry : (TSUB_SORT * SORT_ELEMENTS_PER_THREAD - NEIGHBOR_NUM_INC);
	  const int knum = BLOCKSIZE(krem, TSUB_NEIGHBOR);
	  //---------------------------------------------------------------
	  for(int ki = 0; ki < knum; ki++){
	    //-------------------------------------------------------------
	    const int pidx = lane + ki * TSUB_NEIGHBOR;
	    if( pidx < krem ){
	      //-----------------------------------------------------------
	      flip.f = (pidx < TSUB_NEIGHBOR * NBUF_NEIGHBOR) ? rjbuf[hbuf + pidx] : rjminBuf[bufHead + pidx];
	      sort_tmp[hp_sort + NEIGHBOR_NUM_INC + pidx] = flip32flt(flip.u);
	      //-----------------------------------------------------------
	    }/* if( pidx < krem ){ */
	    //-------------------------------------------------------------
	  }/* for(int ki = 0; ki < knum; ki++){ */
	  //---------------------------------------------------------------
	  /* sort distance */
	  __radixSortTsub32(32 / RADIX_SORT_CHECK_BITS, lane, &sort_tmp[hp_sort]
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
			    , tidx, &sort_smem[hp_sort]
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
			    , &sort_sbuf[head]
#endif//RADIX_SORT_CHECK_BITS > 2
			    );
	  //---------------------------------------------------------------
	}/* for(int iter = 0; iter < Niter; iter++){ */
	//-----------------------------------------------------------------
	if( (iidx < Ni) && (lane == 0) ){
	  //---------------------------------------------------------------
	  /* get NEIGHBOR_NUM-th distance */
	  /* flip.u  = undo32flt(sort_tmp[hp_sort + NEIGHBOR_NUM - 1]); */
	  /* i-particle itself (r2 = ZERO; that means the minimum distance --> the head component in the sorted array) must be eliminated */
	  flip.u  = undo32flt(sort_tmp[hp_sort + NEIGHBOR_NUM]);
	  flip.f *= RSQRT(flip.f);
	  //---------------------------------------------------------------
	  /* return the length of the neighbor arm */
	  neighbor_length[iidx] = flip.f;
	  //---------------------------------------------------------------
	}/* if( (iidx < Ni) && (lane == 0) ){ */
	//-----------------------------------------------------------------
	not_yet = false;
	break;
	//-----------------------------------------------------------------
      }/* if( lowest == TSUB_NEIGHBOR ){ */
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
      /* list up all child nodes that satisfy rjmin < rmax */
      //-------------------------------------------------------------------
      inum = 0;
      //-------------------------------------------------------------------
      Nloc = 0;
      Nbuf = 0;
      //-------------------------------------------------------------------
      Niter = BLOCKSIZE(Ntry, NBUF_NEIGHBOR * TSUB_NEIGHBOR);
      for(int iter = 0; iter < Niter; iter++){
	//-----------------------------------------------------------------
	const int krem = (Ntry < (NBUF_NEIGHBOR * TSUB_NEIGHBOR)) ? Ntry : (NBUF_NEIGHBOR * TSUB_NEIGHBOR);
	const int knum = BLOCKSIZE(krem, TSUB_NEIGHBOR);
	//-----------------------------------------------------------------
	for(int ki = 0; ki < knum; ki++){
	  //---------------------------------------------------------------
	  int cellIdx = lane + ki * TSUB_NEIGHBOR;
	  int  add = 0;
	  int iadd = 0;
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* select distant tree cells */
	  //---------------------------------------------------------------
	  if( cellIdx < krem ){
	    //-------------------------------------------------------------
	    /* when the current node must be taken into account */
	    cellIdx += hbuf;
	    if( rjbuf[cellIdx] < rmax ){
	      //-----------------------------------------------------------
	      /* count up total number of contained i-particles */
	      cellIdx = node2cell[list1[cellIdx]];
	      iadd = cell[cellIdx].num;
	      //-----------------------------------------------------------
	      add = 1;
	      //-----------------------------------------------------------
	    }/* if( rjbuf[cellIdx] < rmax ){ */
	    //-------------------------------------------------------------
	  }/* if( cellIdx < krem ){ */
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* remove duplicated tree cells */
	  //---------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  smem = prefixSumTsub(add, lane);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  prefixSumTsub(add, smem, tidx, lane);
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  if( add ){
	    //-------------------------------------------------------------
	    /* test uploading... */
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    const int smidx = Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    const int smidx = Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    list0[hbuf + smidx] = cellIdx;
	    //-------------------------------------------------------------
	    /* if detect duplication, upload flag is turned off */
	    if( ((smidx > 0) && (list0[hbuf + smidx - 1] == cellIdx)) || ((Nbuf > 0) && (smidx == 0) && (more0Buf[bufHead + Nbuf - 1] == cellIdx)) ){
	      add  = 0;
	      iadd = 0;
	    }/* if( ((smidx > 0) && (list0[hbuf + smidx - 1] == cellIdx)) || ((Nbuf > 0) && (smidx == 0) && (more0Buf[bufHead + Nbuf - 1] == cellIdx)) ){ */
	    //-------------------------------------------------------------
	  }/* if( add ){ */
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* save tree cells on the local buffer */
	  //---------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  smem = prefixSumTsub(add, lane);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  prefixSumTsub(add, smem, tidx, lane);
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  if( add ){
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    const int smidx = Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    const int smidx = Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	    list0[hbuf + smidx] = cellIdx;
	  }/* if( add ){ */
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  Nloc += __shfl(smem, TSUB_NEIGHBOR - 1, TSUB_NEIGHBOR);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  Nloc +=        smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* move data to the remote buffer if necessary */
	  if( Nloc > ((NBUF_NEIGHBOR - 1) * TSUB_NEIGHBOR) ){
	    for(int ll = lane; ll < Nloc; ll += TSUB_NEIGHBOR)
	      more0Buf[bufHead + Nbuf + ll] = list0[hbuf + ll];
	    Nbuf += Nloc;
	    Nloc  = 0;
	  }/* if( Nloc > ((NBUF_NEIGHBOR - 1) * TSUB_NEIGHBOR) ){ */
	  //---------------------------------------------------------------
	  /* sum up iadd within TSUB_NEIGHBOR threads */
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  iadd = accumulateIntTsub(iadd);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  accumulateIntTsub(&iadd, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  inum += iadd;
	  //---------------------------------------------------------------
	}/* for(int ki = 0; ki < knum; ki++){ */
	//-----------------------------------------------------------------
	Ntry -= krem;
	//-----------------------------------------------------------------
	/* copy data from remote buffer to local buffer */
	const int Ncopy = (Ntry < (TSUB_NEIGHBOR * NBUF_NEIGHBOR)) ? (Ntry) : (TSUB_NEIGHBOR * NBUF_NEIGHBOR);
	for(int jj = lane; jj < Ncopy; jj += TSUB_NEIGHBOR){
	  rjbuf[hbuf + jj] = rjminBuf[bufHead + NBUF_NEIGHBOR * TSUB_NEIGHBOR * (iter + 1) + jj];
	  list1[hbuf + jj] = more1Buf[bufHead + NBUF_NEIGHBOR * TSUB_NEIGHBOR * (iter + 1) + jj];
	}/* for(int jj = lane; jj < Ncopy; jj += TSUB_NEIGHBOR){ */
	//-----------------------------------------------------------------
      }/* for(int iter = 0; iter < Niter1; iter++){ */
      //-------------------------------------------------------------------
      if( Nbuf != 0 )
	for(int ll = lane; ll < Nloc; ll += TSUB_NEIGHBOR)
	  more0Buf[bufHead + Nbuf + ll] = list0[hbuf + ll];
      if( Nbuf != 0 )
	for(int ll = lane; ll < NBUF_NEIGHBOR * TSUB_NEIGHBOR; ll += TSUB_NEIGHBOR)
	  list0[hbuf + ll] = more0Buf[bufHead + ll];
      //-------------------------------------------------------------------
      Ntry = Nloc + Nbuf;
      if( (lane == 0) && (Ntry > NUM_ALLOC_NEIGHBORBUF) )
	atomicAdd(overflow, 1);
      //-------------------------------------------------------------------
    }/* while( inum > NI_NEIGHBOR_ESTIMATE ){ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* check positions of all the pick upped i-particles */
    //---------------------------------------------------------------------
    if( not_yet ){
      //-------------------------------------------------------------------
#pragma unroll
      for(int jj = lane; jj < TSUB_NEIGHBOR * SORT_ELEMENTS_PER_THREAD; jj += TSUB_NEIGHBOR)
	sort_tmp[hp_sort + jj] = FLIPPED_REAL_MAX;
      /* Since NI_NEIGHBOR_ESTIMATE <= TSUB_NEIGHBOR * NBUF_NEIGHBOR, Ntry1 is less than TSUB_NEIGHBOR * NBUF_NEIGHBOR */
      /* load index of the pick upped i-particles to list1 */
      const int Niter = BLOCKSIZE(Ntry, TSUB_NEIGHBOR);
      int Ncand = 0;
      for(int iter = 0; iter < Niter; iter++){
	//-----------------------------------------------------------------
	treecell cand;
	int pnum = 0;
	if( lane < Ntry ){
	  cand = cell[list0[hbuf + iter * TSUB_NEIGHBOR + lane]];
	  pnum = cand.num;
	}/* if( lane < Ntry ){ */
	//-----------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	smem = prefixSumTsub(pnum, lane);
	for(int jj = 0; jj < pnum; jj++)
	  list1[hbuf + Ncand + smem         - pnum + jj] = cand.head + jj;
	Ncand += __shfl(smem, TSUB_NEIGHBOR - 1, TSUB_NEIGHBOR);
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	prefixSumTsub(pnum, smem, tidx, lane);
	for(int jj = 0; jj < pnum; jj++)
	  list1[hbuf + Ncand + smem[tidx].i - pnum + jj] = cand.head + jj;
	Ncand +=        smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
        //-----------------------------------------------------------------
	Ntry -= TSUB_NEIGHBOR;
	//-----------------------------------------------------------------
      }/* for(int iter = 0; iter < Niter; iter++){ */
      //-------------------------------------------------------------------
      for(int jj = lane; jj < Ncand; jj += TSUB_NEIGHBOR){
	//-----------------------------------------------------------------
	const position ipos_can = pi[list1[hbuf + jj]];
	//-----------------------------------------------------------------
	const real dx = ipos_can.x - ipos.x;
	const real dy = ipos_can.y - ipos.y;
	const real dz = ipos_can.z - ipos.z;
	const real r2 = 1.0e-30f + dx * dx + dy * dy + dz * dz;
	//-----------------------------------------------------------------
	flip.f = r2;
	sort_tmp[hp_sort + jj] = flip32flt(flip.u);
	//-----------------------------------------------------------------
      }/* for(int jj = lane; jj < Ncand; jj += TSUB_NEIGHBOR){ */
      //-------------------------------------------------------------------
      /* sort distance */
      __radixSortTsub32(32 / RADIX_SORT_CHECK_BITS, lane, &sort_tmp[hp_sort]
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
			, tidx, &sort_smem[hp_sort]
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
			, &sort_sbuf[head]
#endif//RADIX_SORT_CHECK_BITS > 2
			);
      //-------------------------------------------------------------------
      if( (iidx < Ni) && (lane == 0) ){
	//-----------------------------------------------------------------
	/* get NEIGHBOR_NUM-th distance */
	/* flip.u  = undo32flt(sort_tmp[hp_sort + NEIGHBOR_NUM - 1]); */
	/* i-particle itself (r2 = ZERO; that means the minimum distance --> the head component in the sorted array) must be eliminated */
	flip.u  = undo32flt(sort_tmp[hp_sort + NEIGHBOR_NUM]);
	flip.f *= RSQRT(flip.f);
	//-----------------------------------------------------------------
	/* return the length of the neighbor arm */
#if 1
	neighbor_length[iidx] = flip.f;
#if 0
	if( iidx < 32 )
	  printf("iidx = %d: arm = %e from (%e, %e, %e)\n", iidx, flip.f, ipos.x, ipos.y, ipos.z);
#endif
#else
	neighbor_length[iidx] = flip.f * RSQRT(ipos.x * ipos.x + ipos.y * ipos.y + ipos.z * ipos.z);
#endif
	//-----------------------------------------------------------------
      }/* if( (iidx < Ni) && (lane == 0) ){ */
      //-------------------------------------------------------------------
    }/* if( not_yet ){ */
    //---------------------------------------------------------------------
  /* }/\* for(int ii = 0; ii < bnum * BLOCKSIZE(Ni, bnum * NGROUPS_NEIGHBOR); ii += bnum){ *\/ */
  //-----------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID_NEIGHBOR
  releaseBuffer(tidx, freeLst, (uint)bufIdx, target);
#else///USE_SMID_TO_GET_BUFID_NEIGHBOR
#ifdef  TRY_MODE_ABOUT_BUFFER_NEIGHBOR
  releaseBuffer(tidx, freeLst, (uint)bufIdx, target);
#else///TRY_MODE_ABOUT_BUFFER_NEIGHBOR
  releaseBuffer(tidx, freeNum, freeLst, bufIdx, active);
#endif//TRY_MODE_ABOUT_BUFFER_NEIGHBOR
#endif//USE_SMID_TO_GET_BUFID_NEIGHBOR
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* arrays to store properties of tree cells (allocated on the global memory) */
//-------------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID_NEIGHBOR
//-------------------------------------------------------------------------
/* complicated treatments is a remedy for ``not contiguous'' case of smid */
//-------------------------------------------------------------------------
__global__ void initFreeLst4NS(const int numLanes, uint * RESTRICT freeLst, const int numFul, READ_ONLY int * RESTRICT smid)
{
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------
  if( tidx < numFul )
    freeLst[tidx] = INT_MAX;
  //-----------------------------------------------------------------------
  if( tidx < numLanes ){
    //---------------------------------------------------------------------
    const int target = (tidx % NBLOCKS_PER_SM_NEIGHBOR) + smid[tidx / NBLOCKS_PER_SM_NEIGHBOR] * NBLOCKS_PER_SM_NEIGHBOR;
    //---------------------------------------------------------------------
    freeLst[target] = (uint)tidx;
    //---------------------------------------------------------------------
  }/* if( tidx < numLanes ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///USE_SMID_TO_GET_BUFID_NEIGHBOR
//-------------------------------------------------------------------------
__global__ void initFreeLst4NS
(const int numLanes, uint * RESTRICT freeLst
#ifndef TRY_MODE_ABOUT_BUFFER_NEIGHBOR
 , uint * RESTRICT freeNum, int * RESTRICT active
#endif//TRY_MODE_ABOUT_BUFFER_NEIGHBOR
 )
{
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------
  if( tidx < numLanes ){
    //---------------------------------------------------------------------
#ifdef  TRY_MODE_ABOUT_BUFFER_NEIGHBOR
    freeLst[tidx] = (uint)tidx;
#else///TRY_MODE_ABOUT_BUFFER_NEIGHBOR
    freeLst[tidx] = (uint)(numLanes - (tidx + 1));
#endif//TRY_MODE_ABOUT_BUFFER_NEIGHBOR
    //---------------------------------------------------------------------
#ifndef TRY_MODE_ABOUT_BUFFER_NEIGHBOR
    if( tidx == 0 ){
      *freeNum = (uint)numLanes;
      *active  = 1;
    }
#endif//TRY_MODE_ABOUT_BUFFER_NEIGHBOR
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_SMID_TO_GET_BUFID_NEIGHBOR
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* void  freeTreeBuffer_dev */
/* (int  *failure, uint  *buffer, uint  *freeLst */
/* #   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER) */
/*  , uint  *freeNum, int  *active */
/* #endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER) */
/*  ) */
//-------------------------------------------------------------------------
/* extern "C" */
/* muse allocNeighborSearch_dev(int **gsync0, int **gsync1, deviceProp devProp) */
extern "C"
muse allocNeighborSearch_dev
(int **gsync0, int **gsync1, uint **freeLst,
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
 uint **freeNum, int **active,
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
 soaNeighborSearchBuf *buf, deviceProp devProp)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const int nblocks = devProp.numSM * NBLOCKS_PER_SM_NEIGHBOR;
  //-----------------------------------------------------------------------
  mycudaMalloc((void **)gsync0, nblocks * sizeof(int));
  mycudaMalloc((void **)gsync1, nblocks * sizeof(int));
  alloc.device +=               nblocks * sizeof(int);
  alloc.device +=               nblocks * sizeof(int);
  //-----------------------------------------------------------------------
  initGsync_kernel<<<1, nblocks>>>(nblocks, *gsync0, *gsync1);
  getLastCudaError("initGsync_kernel");
  //-----------------------------------------------------------------------
  buf->gsync0 = *gsync0;
  buf->gsync1 = *gsync1;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID_NEIGHBOR
  int last = 0;
  int num = 0;
  int *smid_dev;  mycudaMalloc    ((void **)&smid_dev, sizeof(int) * devProp.numSM);
  int *smid_hst;  mycudaMallocHost((void **)&smid_hst, sizeof(int) * devProp.numSM);
  for(int ii = 0; ii < 64; ii++)
    if( devProp.smid[ii] != -1 ){
      smid_hst[num] = devProp.smid[ii];      num++;
      last = ii;
    }
  last++;
  mycudaMalloc((void **)freeLst, (NBLOCKS_PER_SM_NEIGHBOR * last) * sizeof(uint));
  alloc.device                += (NBLOCKS_PER_SM_NEIGHBOR * last) * sizeof(uint);
#else///USE_SMID_TO_GET_BUFID_NEIGHBOR
  mycudaMalloc((void **)freeLst, nblocks * sizeof(uint));
  alloc.device                += nblocks * sizeof(uint);
#endif//USE_SMID_TO_GET_BUFID_NEIGHBOR
  buf->freeLst = *freeLst;
  //-----------------------------------------------------------------------
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
  mycudaMalloc((void **)freeNum, sizeof(uint));  alloc.device += sizeof(uint);  buf->freeNum = *freeNum;
  mycudaMalloc((void **) active, sizeof( int));  alloc.device += sizeof( int);  buf-> active = * active;
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
  //-----------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID_NEIGHBOR
  checkCudaErrors(cudaMemcpy(smid_dev, smid_hst, sizeof(int) * devProp.numSM, cudaMemcpyHostToDevice));
  initFreeLst4NS<<<1, NBLOCKS_PER_SM_NEIGHBOR * last>>>(nblocks, *freeLst, NBLOCKS_PER_SM_NEIGHBOR * last, smid_dev);
  mycudaFree    (smid_dev);
  mycudaFreeHost(smid_hst);
#else///USE_SMID_TO_GET_BUFID_NEIGHBOR
  initFreeLst4NS<<<1, nblocks>>>(nblocks, *freeLst
#ifndef TRY_MODE_ABOUT_BUFFER_NEIGHBOR
			      , *freeNum, *active
#endif//TRY_MODE_ABOUT_BUFFER_NEIGHBOR
			      );
#endif//USE_SMID_TO_GET_BUFID_NEIGHBOR
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
void  freeNeighborSearch_dev
(int  *gsync0, int  *gsync1, uint  *freeLst
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
 , uint  *freeNum, int  *active
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  mycudaFree(gsync0);
  mycudaFree(gsync1);
  //-----------------------------------------------------------------------
  mycudaFree(freeLst);
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
  mycudaFree(freeNum);
  mycudaFree( active);
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* calculate distance between an i-particle and the corresponding NEIGHBOR_NUM-th neighbor particle */
//-------------------------------------------------------------------------
/* Ni              :: input          :: Number of i-particles */
/* pi              :: input          :: position and mass of N-body particles */
/* neighbor_length ::         output :: distance between an i-particle and the corresponding NEIGHBOR_NUM-th neighbor particle */
/* cell            :: input          :: head index and number of N-body particles contained in the corresponding tree cell */
/* leaf            :: input          :: a flag to remember the corresponding tree cell is either leaf(true) of node(false) */
/* node            :: input          :: head index and number of pseudo particles of the corresponding tree cell (a leaf cell contains up to NCRIT particles) */
/* more            :: input          :: head index and number of child particles of the corresponding j-particle (a leaf cell contains up to NCRIT particles) */
/* node2cell       :: input          :: index of the tree cell corresponding a pseudo particle */
/* pj              :: input          :: position and squared radius of pseudo N-body particle as j-particles */
/* bmax            :: input          :: size of pseudo N-body particle as j-particles */
//-------------------------------------------------------------------------
extern "C"
void searchNeighbors_dev
(const int Ni, const iparticle pi, const soaTreeCell cell, const soaTreeNode node,
 const soaMakeTreeBuf makeBuf, const soaNeighborSearchBuf searchBuf, deviceProp devProp)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set pseudo j-particles */
  //-----------------------------------------------------------------------
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
  /* searchNeighbors_kernel<<<devProp.numSM * NBLOCKS_PER_SM_NEIGHBOR, NTHREADS_NEIGHBOR>>> */
  searchNeighbors_kernel<<<BLOCKSIZE(Ni, NTHREADS_NEIGHBOR / TSUB_NEIGHBOR), NTHREADS_NEIGHBOR>>>
    (Ni, pi.pos, pi.neighbor,
     cell.cell, cell.leaf, cell.ptag, node.more, node.node2cell, node.niSub, node.jpos, node.bmax,
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
     searchBuf.active, searchBuf.freeNum,
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) &&  defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
     searchBuf.freeNum,
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) &&  defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
     searchBuf.freeLst, makeBuf.more0, makeBuf.more1, makeBuf.rjmax, makeBuf.fail, searchBuf.gsync0, searchBuf.gsync1);
  getLastCudaError("searchNeighbors_kernel");
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
  //-----------------------------------------------------------------------
  int fail_hst;
  checkCudaErrors(cudaMemcpy(&fail_hst, makeBuf.fail, sizeof(int), cudaMemcpyDeviceToHost));
  if( fail_hst != 0 ){
    __KILL__(stderr, "ERROR: buffer (%d elements per %d threads group) overflow at least %d times.\nPLEASE re-simulate after increasing NUM_ALLOC_MACBUF defined in src/tree/make.h.\n", NUM_ALLOC_NEIGHBORBUF, TSUB_NEIGHBOR, fail_hst);
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//FACILE_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  SMEM_PREF_FOR_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------
extern "C"
void setGlobalConstants_neighbor_dev_cu(void)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#ifdef  FACILE_NEIGHBOR_SEARCH
  checkCudaErrors(cudaFuncSetCacheConfig(facileNeighborSearching_kernel, cudaFuncCachePreferShared));
#else///FACILE_NEIGHBOR_SEARCH
  checkCudaErrors(cudaFuncSetCacheConfig(        searchNeighbors_kernel, cudaFuncCachePreferShared));
#endif//FACILE_NEIGHBOR_SEARCH
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//SMEM_PREF_FOR_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------
