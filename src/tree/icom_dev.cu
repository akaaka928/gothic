/**
 * @file icom_dev.cu
 *
 * @brief Source code to generate Enclosing Ball containing all N-body particles on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/08/28 (Mon)
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
#include <math.h>
#include <sys/time.h>
#include <helper_cuda.h>

#include "macro.h"
#include "cudalib.h"

#include "../misc/device.h"
#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "icom_dev.h"

#ifdef  OCTREE_BASED_SEARCH
#include "../tree/make.h"
#include "../tree/make_dev.h"
#else///OCTREE_BASED_SEARCH
#include "../sort/peano_dev.h"
#endif//OCTREE_BASED_SEARCH


#   if  defined(USE_WARP_SHUFFLE_FUNC_EB) && !defined(USE_WARP_SHUFFLE_FUNC_COMPARE_INC)
#define  USE_WARP_SHUFFLE_FUNC_COMPARE_INC
#endif//defined(USE_WARP_SHUFFLE_FUNC_EB) && !defined(USE_WARP_SHUFFLE_FUNC_COMPARE_INC)
#define NTHREADS_COMPARE_INC NTHREADS_EB
#include "../util/compare_inc.cu"


#ifdef  OCTREE_BASED_SEARCH

#   if  defined(USE_WARP_SHUFFLE_FUNC_EB) && !defined(USE_WARP_SHUFFLE_FUNC_SCAN_INC)
#define  USE_WARP_SHUFFLE_FUNC_SCAN_INC
#endif//defined(USE_WARP_SHUFFLE_FUNC_EB) && !defined(USE_WARP_SHUFFLE_FUNC_SCAN_INC)
#define NTHREADS_SCAN_INC NTHREADS_EB
#include "../util/scan_inc.cu"


/* #define NTHREADS_COMPARE_VEC3_INC NTHREADS_EB */
/* #include "../util/compare_vec3_inc.cu" */

/**
 * @union alignedInt
 *
 * @brief union for store multiple integers in registers
 */
/**
 * @union alignedFlt
 *
 * @brief union for store multiple floats in registers
 */
#   if  NBUF_EB == 4
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
#endif//NBUF_EB == 4
#   if  NBUF_EB == 2
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
#endif//NBUF_EB == 2

/**
 * @union int_real
 *
 * @brief union for switching integer and float
 */
typedef union
{
  int  i;
  float f;
} int_float;


/**
 * @fn findFarthestParticle
 *
 * @brief Find the farthest particle from the specified point.
 *
 * @param (cen) position of center of the enclosing sphere
 * @return (r2max) squared distance and index of the N-body particle which locates farthest from the center
 * @param (cell) head index and number of N-body particles contained in the corresponding tree cell
 * @param (pi) position and mass of N-body particles
 * @param (node) head index and number of pseudo particles of the corresponding tree cell (a leaf cell contains up to NCRIT particles)
 * @param (more) head index and number of child particles of the corresponding j-particle (a leaf cell contains up to NCRIT particles)
 * @param (node2cell) index of the tree cell corresponding a pseudo particle
 * @param (pj) position and squared radius of pseudo N-body particle as j-particles
 * @param (bmax) size of pseudo N-body particle as j-particles
 * @param (more0Buf) buffer for more
 * @param (more1Buf) buffer for more
 * @param (rjmaxBuf) buffer for rjmax
 * @return (overflow) detector of buffer overflow
 * @param (bufSize) size of buffers on global memory
 */
__device__ __forceinline__ float findFarthestParticle
(const position cen, floc * RESTRICT r2max,
 READ_ONLY treecell * RESTRICT cell, READ_ONLY position * RESTRICT pi,
 READ_ONLY uint * RESTRICT node, READ_ONLY uint * RESTRICT more, READ_ONLY int * RESTRICT node2cell,
 READ_ONLY jparticle * RESTRICT pj, READ_ONLY real * RESTRICT bmax,
 int * RESTRICT more0Buf, int * RESTRICT more1Buf, float * RESTRICT rjmaxBuf, int * RESTRICT overflow, const size_t bufSize)
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;

  const int lane = tidx & (warpSize - 1);
  const int head = tidx - lane;
  const int tail = NTHREADS_EB - 1;

  __shared__  int_float smem [NTHREADS_EB];
  __shared__      float rjbuf[NTHREADS_EB * NBUF_EB];
  __shared__  int       list0[NTHREADS_EB * NBUF_EB];
  __shared__  int       list1[NTHREADS_EB * NBUF_EB];
  __shared__  int       pjidx[NTHREADS_EB * NBUF_EB];


  /** focus only on the root cell */
  const      int cidx = 0;
  const treecell root = cell[cidx];


  uint more_tmp = more[node[cidx] & IDXMASK];
  int cnum  = (more_tmp >> IDXBITS) + 1;
  int chead = (more_tmp  & IDXMASK);


  /** estimate size of particle distribution */
  /** initialize list of examined tree nodes and related variables */
  int inum = root.num;
  int Ntry = 1;
  if( tidx == 0 )
    list0[0] = cidx;
  __syncthreads();


  /** pick up NI_R2MAX_ESTIMATE i-particles in maximum to estimate bmax by reduction within a block */
  while( inum > NI_R2MAX_ESTIMATE ){
    float rmin = 0.0f;
    int Nloc = 0;
    int Nbuf = 0;

    int Niter = BLOCKSIZE(Ntry, NTHREADS_EB * NBUF_EB);
    for(int iter = 0; iter < Niter; iter++){
      const int Nsweep = (Ntry > (NTHREADS_EB * NBUF_EB)) ? (NTHREADS_EB * NBUF_EB) : Ntry;
      const int ibuf_loop = BLOCKSIZE(Nsweep, NTHREADS_EB);
      for(int ibuf = 0; ibuf < ibuf_loop; ibuf++){
	cnum = 0;
	const int target = tidx + ibuf * NTHREADS_EB;
	if( target < Nsweep ){
	  /** load a tree node corresponding the tree cell */
	  more_tmp = node[list0[target]];
	  const int nodenum  = 1 + (more_tmp >> IDXBITS);
	  const int nodehead =      more_tmp  & IDXMASK;

	  /** load all child nodes of the tree cell */
	  more_tmp = more[nodehead];
	  cnum  = 1 + (more_tmp >> IDXBITS);
	  chead =      more_tmp  & IDXMASK;
	  for(int jj = 1; jj < nodenum; jj++)
	    cnum += (1 + (more[nodehead + jj] >> IDXBITS));
	}/* if( target < Nsweep ){ */


	PREFIX_SUM_BLCK(cnum, (int *)smem, lane, tidx);
	const int lend = BLOCKSIZE(smem[tail].i, NTHREADS_EB * NBUF_EB);


	for(int ll = 0; ll < lend; ll++){
	  const int unum =
	    (smem[tidx].i <= (NTHREADS_EB * NBUF_EB)) ? cnum :
	    ((smem[tidx].i >= (cnum + NTHREADS_EB * NBUF_EB)) ? (0) : (cnum + NTHREADS_EB * NBUF_EB - smem[tidx].i));
	  const int shead = smem[tidx].i - cnum;

	  for(int jj = 0; jj < unum; jj++){
	    pjidx[shead + jj] = chead;
	    chead++;
	  }/* for(int jj = 0; jj < unum; jj++){ */
	  cnum -= unum;
	  const int Ntmp = smem[tidx].i - (NTHREADS_EB * NBUF_EB);/**< Ntmp is a temporal buffer */


	  /** pick up candidate tree nodes */
#   if  NBUF_EB == 4
	  alignedFlt rjmax_loc = { FLT_MIN ,  FLT_MIN ,  FLT_MIN ,  FLT_MIN };
	  alignedInt pjidx_loc = {NULL_NODE, NULL_NODE, NULL_NODE, NULL_NODE};
#endif//NBUF_EB == 4
#   if  NBUF_EB == 2
	  alignedFlt rjmax_loc = { FLT_MIN ,  FLT_MIN };
	  alignedInt pjidx_loc = {NULL_NODE, NULL_NODE};
#endif//NBUF_EB == 2
	  const int stail = (smem[tail].i < (NTHREADS_EB * NBUF_EB)) ? (smem[tail].i) : (NTHREADS_EB * NBUF_EB);


#pragma unroll
	  for(int kk = 0; kk < NBUF_EB; kk++){
	    const int jj = tidx + kk * NTHREADS_EB;
	    if( jj >= stail )
	      break;

	    const int kidx = pjidx[jj];
	    const jparticle jpos = pj[kidx];

	    const float dx = CAST_R2F(jpos.x - cen.x);
	    const float dy = CAST_R2F(jpos.y - cen.y);
	    const float dz = CAST_R2F(jpos.z - cen.z);
	    const float d2 = FLT_MIN + dx * dx + dy * dy + dz * dz;
	    const float dr = d2 * rsqrtf(d2);

	    const float rjmax = CAST_R2F(bmax[kidx]) + dr;
	    const float rjmin = -rjmax + (2.0f * (1.0f - FLT_EPSILON)) * dr;
	    rmin = fmaxf(rjmin, rmin);

	    if( rjmax > rmin ){
	      pjidx_loc.ia[kk] = kidx;
	      rjmax_loc.ra[kk] = rjmax;
	    }/* if( rjmax > rmin ){ */
	  }/* for(int kk = 0; kk < NBUF_EB; kk++){ */

	  /** share rmin within NTHREADS_EB threads */
	  rmin = GET_MAX_BLCK(rmin, (float *)smem, tidx, head);


	  /** recheck local buffer (is really rjmax greater than rmin ?) */
#pragma unroll
	  for(int jj = 0; jj < NBUF_EB; jj++){
	    const int share = (rjmax_loc.ra[jj] > rmin) ? 1 : 0;
	    PREFIX_SUM_BLCK(share, (int *)smem, lane, tidx);

	    if( share ){
	      const int dst = Nloc + smem[tidx].i - 1;
	      list1[dst] = pjidx_loc.ia[jj];
	      rjbuf[dst] = rjmax_loc.ra[jj];
	    }/* if( share ){ */
	    Nloc += smem[tail].i;


	    if( Nloc > ((NBUF_EB - 1) * NTHREADS_EB) ){
	      for(int kk = tidx; kk < Nloc; kk += NTHREADS_EB){
		more1Buf[Nbuf + kk] = list1[kk];
		rjmaxBuf[Nbuf + kk] = rjbuf[kk];
	      }/* for(int kk = tidx; kk < Nloc; kk += NTHREADS_EB){ */

	      Nbuf += Nloc;
	      Nloc = 0;
	    }/* if( Nloc > ((NBUF_EB - 1) * TSUB_EB) ){ */
	  }/* for(int jj = 0; jj < NBUF_EB; jj++){ */

	  smem[tidx].i = Ntmp;/**< Ntmp is a temporal buffer */
	}/* for(int ll = 0; ll < lend; ll++){ */
      }/* for(int ibuf = 0; ibuf < NBUF_EB; ibuf++){ */


      Ntry -= Nsweep;


      /** copy data from global memory to shared memory */
      const int Ncopy = (Ntry < (NTHREADS_EB * NBUF_EB)) ? (Ntry) : (NTHREADS_EB * NBUF_EB);
      for(int jj = tidx; jj < Ncopy; jj += NTHREADS_EB)
	list0[jj] = more0Buf[jj + ((NTHREADS_EB * NBUF_EB) * (iter + 1))];
    }/* for(int iter = 0; iter < Niter; iter++){ */

    if( Nbuf != 0 ){
      for(int ll = tidx; ll < Nloc; ll += NTHREADS_EB){
	more1Buf[Nbuf + ll] = list1[ll];
	rjmaxBuf[Nbuf + ll] = rjbuf[ll];
      }/* for(int ll = tidx; ll < Nloc; ll += NTHREADS_EB){ */

      for(int ll = tidx; ll < NTHREADS_EB * NBUF_EB; ll += NTHREADS_EB){
	list1[ll] = more1Buf[ll];
	rjbuf[ll] = rjmaxBuf[ll];
      }/* for(int ll = lane; ll < TSUB_EB * NBUF_EB; ll += TSUB_EB){ */
    }/* if( Nbuf != 0 ){ */

    Ntry = Nbuf + Nloc;
    if( (tidx == 0) && (Ntry > bufSize) )
      atomicAdd(overflow, 1);


    /** list up all child nodes that satisfy rjmax > rmin */
    inum = 0;
    Nloc = 0;
    Nbuf = 0;

    Niter = BLOCKSIZE(Ntry, NTHREADS_EB * NBUF_EB);
    for(int iter = 0; iter < Niter; iter++){
      const int krem = (Ntry < (NBUF_EB * NTHREADS_EB)) ? Ntry : (NBUF_EB * NTHREADS_EB);
      const int knum = BLOCKSIZE(krem, NTHREADS_EB);

      for(int ki = 0; ki < knum; ki++){
	int cellIdx = tidx + ki * NTHREADS_EB;
	int  add = 0;
	int iadd = 0;

	/** select distant tree cells */
	if( cellIdx < krem ){
	  /** when the current node must be taken into account */
	  if( rjbuf[cellIdx] > rmin ){
	    /** count up total number of contained i-particles */
	    cellIdx = node2cell[list1[cellIdx]];
	    iadd = cell[cellIdx].num;

	    add = 1;
	  }/* if( rjbuf[cellIdx] > rmin ){ */
	}/* if( cellIdx < krem ){ */

	/** remove duplicated tree cells */
	PREFIX_SUM_BLCK(add, (int *)smem, lane, tidx);

	if( add ){
	  /** test uploading... */
	  const int smidx = Nloc + smem[tidx].i - 1;
	  list0[smidx] = cellIdx;

	  /** if detect duplication, upload flag is turned off */
	  if( ((smidx > 0) && (list0[smidx - 1] == cellIdx)) || ((Nbuf > 0) && (smidx == 0) && (more0Buf[Nbuf - 1] == cellIdx)) ){
	    add  = 0;
	    iadd = 0;
	  }/* if( ((smidx > 0) && (list0[smidx - 1] == cellIdx)) || ((Nbuf > 0) && (smidx == 0) && (more0Buf[Nbuf - 1] == cellIdx)) ){ */
	}/* if( add ){ */

	/** save tree cells on the local buffer */
	PREFIX_SUM_BLCK(add, (int *)smem, lane, tidx);

	if( add ){
	  const int smidx = Nloc + smem[tidx].i - 1;
	  list0[smidx] = cellIdx;
	}/* if( add ){ */

	Nloc += smem[tail].i;

	/** move data to the remote buffer if necessary */
	if( Nloc > ((NBUF_EB - 1) * NTHREADS_EB) ){
	  for(int ll = tidx; ll < Nloc; ll += NTHREADS_EB)
	    more0Buf[Nbuf + ll] = list0[ll];
	  Nbuf += Nloc;
	  Nloc  = 0;
	}/* if( Nloc > ((NBUF_EB - 1) * NTHREADS_EB) ){ */

	/** sum up iadd within NTHREADS_EB threads */

	iadd = TOTAL_SUM_BLCK(iadd, (int *)smem, tidx, head);

	inum += iadd;
      }/* for(int ki = 0; ki < knum; ki++){ */

      Ntry -= krem;

      /** copy data from remote buffer to local buffer */
      const int Ncopy = (Ntry < (NTHREADS_EB * NBUF_EB)) ? (Ntry) : (NTHREADS_EB * NBUF_EB);
      for(int jj = tidx; jj < Ncopy; jj += NTHREADS_EB){
	rjbuf[jj] = rjmaxBuf[jj + NBUF_EB * NTHREADS_EB * (iter + 1)];
	list1[jj] = more1Buf[jj + NBUF_EB * NTHREADS_EB * (iter + 1)];
      }/* for(int jj = tidx; jj < Ncopy; jj += NTHREADS_EB){ */
    }/* for(int iter = 0; iter < Niter1; iter++){ */

    if( Nbuf != 0 ){
      for(int ll = tidx; ll < Nloc; ll += NTHREADS_EB)
	more0Buf[Nbuf + ll] = list0[ll];

      for(int ll = tidx; ll < NBUF_EB * NTHREADS_EB; ll += NTHREADS_EB)
	list0[ll] = more0Buf[ll];
    }/* if( Nbuf != 0 ){ */

    Ntry = Nloc + Nbuf;
    if( (tidx == 0) && (Ntry > bufSize) )
      atomicAdd(overflow, 1);
  }/* while( inum > NI_R2MAX_ESTIMATE ){ */


  /** load index of the pick upped i-particles to list1 */
  const int Niter = BLOCKSIZE(Ntry, NTHREADS_EB);
  int Ncand = 0;
  for(int iter = 0; iter < Niter; iter++){
    treecell cand;
    int pnum = 0;
    if( tidx < Ntry ){
      cand = cell[list0[tidx + iter * NTHREADS_EB]];
      pnum = cand.num;
    }/* if( tidx < Ntry ){ */

    PREFIX_SUM_BLCK(pnum, (int *)smem, lane, tidx);

    for(int jj = 0; jj < pnum; jj++)
      list1[Ncand + smem[tidx].i - pnum + jj] = cand.head + jj;
    Ncand +=        smem[tail].i;

    Ntry -= NTHREADS_EB;
  }/* for(int iter = 0; iter < Niter; iter++){ */


  r2max->val = -FLT_MAX;
  r2max->idx =  INT_MAX;
  for(int jj = tidx; jj < Ncand; jj += NTHREADS_EB){
    const int idx = list1[jj];
    const position ipos = pi[idx];

    const float dx = CAST_R2F(ipos.x - cen.x);
    const float dy = CAST_R2F(ipos.y - cen.y);
    const float dz = CAST_R2F(ipos.z - cen.z);
    const float r2 = FLT_MIN + dx * dx + dy * dy + dz * dz;

    if( r2 > r2max->val ){
      r2max->val = r2;
      r2max->idx = idx;
    };/* if( r2 > r2max->val ){ */

  }/* for(int jj = tidx; jj < Ncand; jj += NTHREADS_EB){ */


  *r2max = GET_MAXLOC_BLCK(*r2max, (floc *)rjbuf, tidx, head);


  return (r2max->val);
}


/**
 * @fn getApproxEnclosingBall_kernel
 *
 * @brief Get approximated enclosing ball (efficient bounding sphere proposed by Ritter 1990).
 *
 * @return (ebs) position of center and squared radius of the enclosing sphere
 * @param (num) number of N-body particles
 * @param (ipos) position and mass of N-body particles
 * @param (cell) head index and number of N-body particles contained in the corresponding tree cell
 * @param (node) head index and number of pseudo particles of the corresponding tree cell (a leaf cell contains up to NCRIT particles)
 * @param (more) head index and number of child particles of the corresponding j-particle (a leaf cell contains up to NCRIT particles)
 * @param (node2cell) index of the tree cell corresponding a pseudo particle
 * @param (pj) position and squared radius of pseudo N-body particle as j-particles
 * @param (bmax) size of pseudo N-body particle as j-particles
 * @param (more0Buf) buffer for more
 * @param (more1Buf) buffer for more
 * @param (rjmaxBuf) buffer for rjmax
 * @return (overflow) detector of buffer overflow
 * @param (bufSize) size of buffers on global memory
 */
__global__ void getApproxEnclosingBall_kernel
(position * RESTRICT ebs, const int num, READ_ONLY position * RESTRICT ipos,
 READ_ONLY treecell * RESTRICT cell, READ_ONLY uint * RESTRICT node, READ_ONLY uint * RESTRICT more, READ_ONLY int * RESTRICT node2cell,
 READ_ONLY jparticle * RESTRICT pj, READ_ONLY real * RESTRICT bmax,
 int * RESTRICT more0Buf, int * RESTRICT more1Buf, float * RESTRICT rjmaxBuf, int * RESTRICT overflow, const size_t bufSize)
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  /* const int head = tidx - (tidx & (warpSize - 1)); */

  __shared__ position smem;
  if( tidx == 0 )
    smem = *ebs;

  __syncthreads();
  position cen = smem;
  cen.m = ZERO;

  floc r2max;
  /* float rold = CAST_R2F(cen.m) * rsqrtf(CAST_R2F(cen.m)); */
  float rold = 0.0f;
  while( findFarthestParticle(cen, &r2max, cell, ipos, node, more, node2cell, pj, bmax, more0Buf, more1Buf, rjmaxBuf, overflow, bufSize) > cen.m ){
    /** update the sphere */
    const float dinv = rsqrtf(r2max.val);
    const float dmax = dinv * r2max.val;

    /** calculate the displacement */
    const float rnew = 0.5f * (rold + dmax);
    const float disp = (dmax - rnew) * dinv;

#if 1
    /** work around to escape from the infinite loop */
    if ( disp < FLT_EPSILON )
      break;
#endif

    /** shift the center of the sphere */

    if( tidx == 0 )
      smem = ipos[r2max.idx];
    __syncthreads();
    const position far = smem;

    cen.x += (far.x - cen.x) * disp;
    cen.y += (far.y - cen.y) * disp;
    cen.z += (far.z - cen.z) * disp;

    /** enlarge the sphere */
    cen.m = rnew * rnew;
    rold = rnew;
  }

  if( tidx == 0 )
    *ebs = cen;
}

#else///OCTREE_BASED_SEARCH

/**# future plan: check only activated i-particles using laneInfo */
/* #include "walk_dev.h"/\**< to read NTHREADS, NWARP and TSUB *\/ */


#include "../util/gsync_dev.cu"



#   if  GPUGEN >= 60
/** capacity of shared memory is 64KiB per SM on newer GPUs */
/** real4 smem[NTHREADS_EB] corresponds 16 * NTHREADS_EB bytes */
#define NBLOCKS_PER_SM_EB (4096 / NTHREADS_EB)
#else///GPUGEN >= 60
/** in L1 cache preferred configuration, capacity of shared memory is 16KiB per SM on older GPUs */
/** real4 smem[NTHREADS_EB] corresponds 16 * NTHREADS_EB bytes */
#define NBLOCKS_PER_SM_EB (1024 / NTHREADS_EB)
#endif//GPUGEN >= 60

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

/** limitation from number of registers */
#   if  NBLOCKS_PER_SM_EB > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_EB * NTHREADS_EB))
#undef  NBLOCKS_PER_SM_EB
#define NBLOCKS_PER_SM_EB   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_EB * NTHREADS_EB))
#endif//NBLOCKS_PER_SM_EB > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_EB * NTHREADS_EB))

/** maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_EB > (MAX_THREADS_PER_SM / NTHREADS_EB)
#undef  NBLOCKS_PER_SM_EB
#define NBLOCKS_PER_SM_EB   (MAX_THREADS_PER_SM / NTHREADS_EB)
#endif//NBLOCKS_PER_SM_EB > (MAX_THREADS_PER_SM / NTHREADS_EB)

/** maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_EB > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_EB
#define NBLOCKS_PER_SM_EB   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_EB > MAX_BLOCKS_PER_SM

/** maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_EB > ((MAX_WARPS_PER_SM * 32) / NTHREADS_EB)
#undef  NBLOCKS_PER_SM_EB
#define NBLOCKS_PER_SM_EB   ((MAX_WARPS_PER_SM * 32) / NTHREADS_EB)
#endif//NBLOCKS_PER_SM_EB > ((MAX_WARPS_PER_SM * 32) / NTHREADS_EB)

/** # of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_EB < 1
#undef  NBLOCKS_PER_SM_EB
#define NBLOCKS_PER_SM_EB  (1)
#endif//NBLOCKS_PER_SM_EB < 1


/**
 * @union alignedInt
 *
 * @brief union for store multiple integers in registers
 */
/**
 * @union alignedFlt
 *
 * @brief union for store multiple floats in registers
 */
#   if  NBUF_EB == 4
typedef union __align__(16)
{
  int4 i4;
  int  ia[4];
} alignedInt;
typedef union __align__(16)
{
  float4 f4;
  float  fa[4];
} alignedFlt;
#endif//NBUF_EB == 4
#   if  NBUF_EB == 2
typedef union __align__(8)
{
  int2 i4;
  int  ia[2];
} alignedInt;
typedef union __align__(8)
{
  float2 f4;
  float  fa[2];
} alignedFlt;
#endif//NBUF_EB == 2


/**
 * @fn findFarthestParticle
 *
 * @brief Find the farthest particle from the specified point.
 *
 * @param (cen) position of center of the enclosing sphere
 * @return (r2max) squared distance and index of the N-body particle which locates farthest from the center
 * @param (ipos) position and mass of N-body particles
 */
__device__ __forceinline__ float findFarthestParticle
(const position cen, floc * RESTRICT r2max, READ_ONLY position * RESTRICT ipos, const int ihead, const int itail,
   volatile floc * RESTRICT smem, volatile floc * RESTRICT gmem, const int tidx, const int head, const int bidx, const int bnum, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  r2max->val = -FLT_MAX;
  r2max->idx =  INT_MAX;

  /** local reduction */
  for(int ih = ihead; ih < itail; ih += NTHREADS_EB){
    const int ii = ih + tidx;
    if( ii < itail ){
      const position pi = ipos[ii];
      const float px = CAST_R2F(pi.x - cen.x);
      const float py = CAST_R2F(pi.y - cen.y);
      const float pz = CAST_R2F(pi.z - cen.z);
      const float r2 = FLT_MIN + px * px + py * py + pz * pz;

      if( r2 > r2max->val ){
	r2max->val = r2;
	r2max->idx = ii;
      }/* if( r2 > r2max->val ){ */

    }/* if( ii < itail ){ */
  }/* for(int ih = ihead; ih < itail; ih += NTHREADS_EB){ */

  /** global reduction */
  *r2max = GET_MAXLOC_GRID(*r2max, smem, tidx, head, gmem, bidx, bnum, gsync0, gsync1);

  return (r2max->val);
}


/**
 * @fn getApproxEnclosingBall_kernel
 *
 * @brief Get approximated enclosing ball (efficient bounding sphere proposed by Ritter 1990).
 *
 * @return (ebs) position of center and squared radius of the enclosing sphere
 * @param (num) number of N-body particles
 * @param (ipos) position and mass of N-body particles
 */
__global__ void __launch_bounds__(NTHREADS_EB, NBLOCKS_PER_SM_EB) getApproxEnclosingBall_kernel
     (position * RESTRICT ebs, const int num, READ_ONLY position * RESTRICT ipos, float4 * RESTRICT min_all, float4 * RESTRICT max_all, floc * RESTRICT gmem,
      int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int bnum =   GRIDDIM_X1D;

  const int head = tidx - (tidx & (warpSize - 1));

  const int ihead = (num *      bidx ) / bnum;
  const int itail = (num * (1 + bidx)) / bnum;


  /** static allocation of the shared memory */
  __shared__ floc smem[NTHREADS_EB];

  __shared__ position ebs_sm;
  if( tidx == 0 )
    ebs_sm = *ebs;

  __syncthreads();
  position cen = ebs_sm;
  cen.m = ZERO;

  floc r2max;
  /* float rold = CAST_R2F(cen.m) * rsqrtf(CAST_R2F(cen.m)); */
  float rold = 0.0f;
  while( findFarthestParticle(cen, &r2max, ipos, ihead, itail, smem, gmem, tidx, head, bidx, bnum, gsync0, gsync1) > cen.m ){
    /** update the sphere */
    const float dinv = rsqrtf(r2max.val);
    const float dmax = dinv * r2max.val;

    /** calculate the displacement */
    const float rnew = 0.5f * (rold + dmax);
    const float disp = (dmax - rnew) * dinv;

#if 1
    /** work around to escape from the infinite loop */
    if ( disp < FLT_EPSILON )
      break;
#endif

    /** shift the center of the sphere */

    if( tidx == 0 )
      ebs_sm = ipos[r2max.idx];
    __syncthreads();
    const position far = ebs_sm;

    cen.x += (far.x - cen.x) * disp;
    cen.y += (far.y - cen.y) * disp;
    cen.z += (far.z - cen.z) * disp;

    /** enlarge the sphere */
    cen.m = rnew * rnew;
    rold = rnew;
  }

  if( tidx + bidx == 0 )
    *ebs = cen;
}


static inline void checkDeadLockCondition(void (*func)(position * RESTRICT ebs, const int num, READ_ONLY position * RESTRICT ipos, float4 * RESTRICT min_all, float4 * RESTRICT max_all, floc * RESTRICT gmem, int * RESTRICT gsync0, int * RESTRICT gsync1), const char *name, const char *file, const int line, const char *call, const int Nregs, const int Ttot, const bool warpShuffle, const int Bguess)
{
  struct cudaFuncAttributes funcAttr;
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, func));

  if( funcAttr.numRegs != REGISTERS_PER_THREAD_EB ){
    fprintf(stderr, "%s(%d): %s\n", file, line, call);
    fprintf(stderr, "warning: # of registers used (%d) in %s is not match with the predicted value (%d).\n", funcAttr.numRegs, name, Nregs);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, Ttot = %d.\n", GPUGEN, GPUVER, Ttot);
    fprintf(stderr, "note: warp shuffle instruction is %s\n", warpShuffle ? " enabled" : "disabled");
    fflush (stderr);
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_EB ){ */

  int regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * Ttot);
  if( regLimit > (MAX_REGISTERS_PER_SM / Ttot) )
    regLimit = (MAX_REGISTERS_PER_SM / Ttot);

  int memLimit = SMEM_SIZE_SM_PREF / funcAttr.sharedSizeBytes;

  int Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / Ttot) )    Nblck = MAX_THREADS_PER_SM / Ttot;
  if( Nblck >   MAX_BLOCKS_PER_SM               )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / Ttot) )    Nblck = ((MAX_WARPS_PER_SM * 32) / Ttot);

  if( Nblck != Bguess ){
    __KILL__(stderr, "ERROR: # of blocks per SM for %s is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of # of blocks in %s is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", name, Nblck, regLimit, memLimit, file, Bguess, funcAttr.numRegs, Nregs, GPUGEN, GPUVER);
  }/* if( Nblck != Bguess ){ */
}


/**
 * @fn allocApproxEnclosingBall_dev
 *
 * @brief Memory allocation for enclosing ball generator for LET.
 */
extern "C"
muse allocApproxEnclosingBall_dev(void **dev, const deviceProp devProp)
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  /** memory allocation and simple confirmation */
  const size_t num = devProp.numSM * NBLOCKS_PER_SM_EB;
  mycudaMalloc(dev, num * sizeof(floc));
  alloc.device +=   num * sizeof(floc);

  /** error checking before running the kernel */
  checkDeadLockCondition(getApproxEnclosingBall_kernel, "getApproxEnclosingBall_kernel", (const char *)__FILE__, __LINE__, (const char *)__func__, REGISTERS_PER_THREAD_EB, NTHREADS_EB,
#ifdef  USE_WARP_SHUFFLE_FUNC_EB
			   true,
#else///USE_WARP_SHUFFLE_FUNC_EB
			   false,
#endif//USE_WARP_SHUFFLE_FUNC_EB
			   NBLOCKS_PER_SM_EB);


  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freeApproxEnclosingBall_dev
 *
 * @brief Memory deallocation for enclosing ball generator for LET.
 */
extern "C"
void  freeApproxEnclosingBall_dev(void *dev)
{
  __NOTE__("%s\n", "start");

  mycudaFree(dev);

  __NOTE__("%s\n", "end");
}

#endif//OCTREE_BASED_SEARCH


/**
 * @fn getApproxEnclosingBall_dev
 *
 * @brief Get approximated enclosing ball (efficient bounding sphere proposed by Ritter 1990).
 *
 * @return (ebs) position of center and squared radius of the enclosing sphere
 * @param (num) number of N-body particles
 * @param (ipos) position and mass of N-body particles
 * @param (cell) head index and number of N-body particles contained in the corresponding tree cell
 * @param (node) head index and number of pseudo particles of the corresponding tree cell (a leaf cell contains up to NCRIT particles)
 * @param (more) head index and number of child particles of the corresponding j-particle (a leaf cell contains up to NCRIT particles)
 * @param (node2cell) index of the tree cell corresponding a pseudo particle
 * @param (pj) position and squared radius of pseudo N-body particle as j-particles
 * @param (bmax) size of pseudo N-body particle as j-particles
 * @param (more0Buf) buffer for more
 * @param (more1Buf) buffer for more
 * @param (rjmaxBuf) buffer for rjmax
 * @return (overflow) detector of buffer overflow
 * @param (bufSize) size of buffers on global memory
 */
extern "C"
void getApproxEnclosingBall_dev
(const int num, const iparticle body
#ifdef  OCTREE_BASED_SEARCH
 , const soaTreeCell cell, const soaTreeNode node, const soaMakeTreeBuf buf
#else///OCTREE_BASED_SEARCH
 , void *gmem, const soaPHsort soa, const deviceProp devProp
#endif//OCTREE_BASED_SEARCH
 , const cudaStream_t stream
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
 )
{

#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK


#ifdef  OCTREE_BASED_SEARCH

  getApproxEnclosingBall_kernel<<<1, NTHREADS_EB, SMEM_SIZE, stream>>>
    (body.encBall, num, body.pos,
     cell.cell, cell.ptag, node.more, node.node2cell, node.jpos, node.bmax,
     buf.more0, buf.more1, buf.rjmax, buf.fail, buf.Nbuf);
  getLastCudaError("getApproxEnclosingBall_kernel");

  int fail_hst;
  checkCudaErrors(cudaMemcpy(&fail_hst, buf.fail, sizeof(int), cudaMemcpyDeviceToHost));
  if( fail_hst != 0 ){
    __KILL__(stderr, "ERROR: buffer overflow at least %d times.\nPLEASE re-simulate after increasing NUM_ALLOC_MACBUF defined in src/tree/make.h.\n", fail_hst);
  }/* if( fail_hst != 0 ){ */

#else///OCTREE_BASED_SEARCH

  getApproxEnclosingBall_kernel<<<devProp.numSM * NBLOCKS_PER_SM_EB, NTHREADS_EB, SMEM_SIZE, stream>>>
    (body.encBall, num, body.pos, soa.min, soa.max, (floc *)gmem, soa.gsync0, soa.gsync1);
  getLastCudaError("getApproxEnclosingBall_kernel");

#endif//OCTREE_BASED_SEARCH


  checkCudaErrors(cudaMemcpy(body.encBall_hst, body.encBall, sizeof(position), cudaMemcpyDeviceToHost));


  /**# EBS の初期推定値は毎ステップ固定 or 前ステップで採用した値で良いと思う． */
  /**# つまり，getEBS() については初期の推定値を外から渡してあげる形式に書き換える． */
  /**# そうすると，後半部分の振る舞いが geometric の場合と書き変わるだけだから，define を使って関数の中身を ON/OFF するだけで切り換えられることになるので，その方が楽だと思う． */
  /**# EBS を与えるコードができていれば，そいつに center を渡して上げるだけで簡単にできるはず */
  /**# box size は src/para/exchange_dev.cu 中の getBoxSize_kernel() を動かせば結果が得られる． */
  /**# この box size は tree rebuild 毎に評価する程度にサボっても良いと思う． */
#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->getApproxEnclosingBall_kernel));
#endif//EXEC_BENCHMARK
}
