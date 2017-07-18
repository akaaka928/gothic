/**
 * @file icom_dev.cu
 *
 * @brief Source code to generate Enclosing Ball containing all N-body particles on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/06/26 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

/**# Smallest Enclosing Ball の実装は後回しにして，その他の Enclosing Balls の実装を先に済ませる． */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <helper_cuda.h>

#include "macro.h"
#include "cudalib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "../util/gsync_dev.cu"

/** in L1 cache preferred configuration, capacity of shared memory is 16KiB per SM */
/** real4 smem[NTHREADS_EB] corresponds 16 * NTHREADS_EB bytes */
#define NBLOCKS_PER_SM_EB (1024 / NTHREADS_EB)

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



#define NTHREADS_COMPARE_VEC3_INC NTHREADS_EB
#include "../util/compare_vec3_inc.cu"
/* getMinMax -->> GET_VEC3_MIN_MAX_GRID */


#define NTHREADS_COMPARE_INC NTHREADS_EB
#include "../util/compare_inc.cu"
/* getMaxLoc -->> GET_MAXLOC_3VALS_GRID */
/* getMinLocMaxLoc -->> GET_MINLOC_MAXLOC_3VALS_GRID */


/**
 * @fn findFurthestParticle
 *
 * @brief Find the most distant particle.
 * @detail algorithm is based on Ritter (1990), ``An efficient bounding sphere'', Graphics Gems
 */
__device__ __forceinline__ float findFurthestParticle
(floc * RESTRICT r2max, READ_ONLY position * RESTRICT ipos, const int ihead, const int itail,
   volatile floc * RESTRICT smem, volatile floc * RESTRICT gmem, const int tidx, const int head, const int bidx, const int bnum, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  floc r2max = {-FLT_MAX, INT_MAX};

  for(int ih = ihead; ih < itail; ih += NTHREADS_EB){
    const int ii = ih + tidx;
    if( ii < itail ){
      const position pi = ipos[ii];
      const float px = CAST_R2F(pi.x);
      const float py = CAST_R2F(pi.y);
      const float pz = CAST_R2F(pi.z);
      const float r2 = px * px + py * py + pz * pz;

      if( r2 > r2max.val ){
	r2max.val = r2;
	r2max.idx = ii;
      }/* if( r2 > r2max.val ){ */

    }/* if( ii < itail ){ */
  }/* for(int ih = ihead; ih < itail; ih += NTHREADS_EB){ */


  /* __shared__ floc smem[NTHREADS_EB]; */
  r2max = GET_MAXLOC_GRID(r2max, smem, tidx, head, gmem, bidx, bnum, gsync0, gsync1);



  return (r2max.val);
}


/**
 * @fn getEBS
 *
 * @brief Get efficient bounding sphere proposed by Ritter (1990)
 */
__global__ void __launch_bounds__(NTHREADS_EB, NBLOCKS_PER_SM_EB) getEBS
     (position * RESTRICT ebs, const int num, READ_ONLY position * RESTRICT ipos, float4 * RESTRICT min_all, float4 * RESTRICT max_all, floc * RESTRICT gmem,
      int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  /** static allocation of the shared memory */
  __shared__ float4 smem[NTHREADS_EB];

  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int bnum =   GRIDDIM_X1D;

  const int head = tidx - (tidx & (warpSize - 1));
  /* const int uidx = tidx + (NTHREADS_EB >> 1); */
  /* const int hidx = tidx - (tidx & (warpSize - 1)); */

  const int ihead = (num *      bidx ) / bnum;
  const int itail = (num * (1 + bidx)) / bnum;


  /** calculate required box size to contain all N-body particles in the local domain */
  floc xmin = { FLT_MAX, INT_MAX};  floc ymin = { FLT_MAX, INT_MAX};  floc zmin = { FLT_MAX, INT_MAX};
  floc xmax = {-FLT_MAX, INT_MAX};  floc ymax = {-FLT_MAX, INT_MAX};  floc zmax = {-FLT_MAX, INT_MAX};

  for(int ih = ihead; ih < itail; ih += NTHREADS_EB){
    const int ii = ih + tidx;
    if( ii < itail ){
      const position pi = ipos[ii];
      const float px = CAST_R2F(pi.x);
      const float py = CAST_R2F(pi.y);
      const float pz = CAST_R2F(pi.z);

      if( px < xmin.val ){	xmin.val = px;	xmin.idx = ii;	    }
      if( py < ymin.val ){	ymin.val = py;	ymin.idx = ii;	    }
      if( pz < zmin.val ){	zmin.val = pz;	zmin.idx = ii;	    }

      if( px > xmax.val ){	xmax.val = px;	xmax.idx = ii;	    }
      if( py > ymax.val ){	ymax.val = py;	ymax.idx = ii;	    }
      if( pz > zmax.val ){	zmax.val = pz;	zmax.idx = ii;	    }
    }/* if( ii < itail ){ */
  }/* for(int ih = ihead; ih < itail; ih += NTHREADS_EB){ */


  position cen;
  float rold;
  {
    float4 min = {xmin.val, ymin.val, zmin.val, 0.0f};
    float4 max = {xmax.val, ymax.val, zmax.val, 0.0f};

    GET_VEC3_MIN_MAX_GRID(&min, &max, smem, tidx, head, min_all, max_all, bidx, bnum, gsync0, gsync1);

    cen.x = HALF * (min.x + max.x);
    cen.y = HALF * (min.y + max.y);
    cen.z = HALF * (min.z + max.z);
    rold = SEB_TINY + 0.5f * fmaxf(max.x - min.x, fmaxf(max.y - min.y, max.z - min.z));
  }
  cen.m = rold * rold;

  floc r2max;
  while( findFurthestParticle(&r2max, ipos, ihead, itail, (floc *)smem, gmem, tidx, head, bidx, bnum, gsync0, gsync1) > cen.m ){
    /** update the sphere */
    const float dinv = rsqrtf(r2max.val);
    const float dmax = dinv * r2max.val;

    /** calculate the displacement */
    const float rnew = 0.5f * (rold + dmax);
    const float disp = (dmax - rnew) * dinv;

#if 1
    /** work around to escape from the infinite loop */
    if ( disp < EPSILON )
      break;
#endif

    /** shift the center of the sphere */

    if( tidx == 0 )
      smem[0] = ipos[r2max.idx].pi;
    __syncthreads();
    const position far = smem[0];

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

















/**
 * @fn getMaximumFloatWarp
 *
 * @brief Get maximum value within a warp.
 * @detail implicit synchronization within a warp (a group of 32 threads) is assumed
 */
__device__ __forceinline__ float getMaximumFloatWarp
(
#ifdef  USE_WARP_SHUFFLE_FUNC_EB
 const float max
#else///USE_WARP_SHUFFLE_FUNC_EB
 float max, volatile float * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC_EB
 )
{
  float tmp;
#ifdef  USE_WARP_SHUFFLE_FUNC_EB
  float val = max;
  tmp = __shfl_xor(val,  1, warpSize);  val = fmaxf(val, tmp);
  tmp = __shfl_xor(val,  2, warpSize);  val = fmaxf(val, tmp);
  tmp = __shfl_xor(val,  4, warpSize);  val = fmaxf(val, tmp);
  tmp = __shfl_xor(val,  8, warpSize);  val = fmaxf(val, tmp);
  tmp = __shfl_xor(val, 16, warpSize);  val = fmaxf(val, tmp);
  return (__shfl(val, 0, warpSize));
#else///USE_WARP_SHUFFLE_FUNC_EB
  smem[tidx] = max;
  tmp = smem[tidx ^  1];  max = fmaxf(max, tmp);  smem[tidx] = max;
  tmp = smem[tidx ^  2];  max = fmaxf(max, tmp);  smem[tidx] = max;
  tmp = smem[tidx ^  4];  max = fmaxf(max, tmp);  smem[tidx] = max;
  tmp = smem[tidx ^  8];  max = fmaxf(max, tmp);  smem[tidx] = max;
  tmp = smem[tidx ^ 16];  max = fmaxf(max, tmp);  smem[tidx] = max;
  return (smem[head]);
#endif///USE_WARP_SHUFFLE_FUNC_EB
}

/**
 * @fn getMaximumFloatBlock
 *
 * @brief Get maximum value within a block.
 */
__device__ __forceinline__ float getMaximumFloatBlock(float val, const int tidx, const int lane, volatile float * RESTRICT smem)
{
  /** 1. maximum within a warp */
  val = getMaximumFloatWarp(val
#ifndef USE_WARP_SHUFFLE_FUNC_EB
			    , smem, tidx, lane
#endif//USE_WARP_SHUFFLE_FUNC_EB
			    );
  smem[tidx] = val;

  /** 2. maximum among warps */
  __syncthreads();

  /** warpSize = 32 = 2^5; NTHREADS_EB <= 1024 --> NTHREADS_EB >> 5 <= 32 = warpSize */
  if( tidx < (NTHREADS_EB >> 5) ){
    val = smem[tidx * warpSize];
    float tmp;
#ifdef  USE_WARP_SHUFFLE_FUNC_EB
    const int groupSize = NTHREADS_EB >> 5;
#   if  (NTHREADS_EB >> 5) >=  2
    tmp = __shfl_xor(val,  1, groupSize);    val = fmaxf(val, tmp);
#   if  (NTHREADS_EB >> 5) >=  4
    tmp = __shfl_xor(val,  2, groupSize);    val = fmaxf(val, tmp);
#   if  (NTHREADS_EB >> 5) >=  8
    tmp = __shfl_xor(val,  4, groupSize);    val = fmaxf(val, tmp);
#   if  (NTHREADS_EB >> 5) >= 16
    tmp = __shfl_xor(val,  8, groupSize);    val = fmaxf(val, tmp);
#   if  (NTHREADS_EB >> 5) == 32
    tmp = __shfl_xor(val, 16, groupSize);    val = fmaxf(val, tmp);
#endif//(NTHREADS_EB >> 5) == 32
#endif//(NTHREADS_EB >> 5) >= 16
#endif//(NTHREADS_EB >> 5) >=  8
#endif//(NTHREADS_EB >> 5) >=  4
#endif//(NTHREADS_EB >> 5) >=  2
    smem[tidx] = val;
#else///USE_WARP_SHUFFLE_FUNC_EB
    smem[tidx] = val;
#   if  (NTHREADS_EB >> 5) >=  2
    tmp = smem[tidx ^  1];    max = fmaxf(max, tmp);    smem[tidx] = max;
#   if  (NTHREADS_EB >> 5) >=  4
    tmp = smem[tidx ^  2];    max = fmaxf(max, tmp);    smem[tidx] = max;
#   if  (NTHREADS_EB >> 5) >=  8
    tmp = smem[tidx ^  4];    max = fmaxf(max, tmp);    smem[tidx] = max;
#   if  (NTHREADS_EB >> 5) >= 16
    tmp = smem[tidx ^  8];    max = fmaxf(max, tmp);    smem[tidx] = max;
#   if  (NTHREADS_EB >> 5) == 32
    tmp = smem[tidx ^ 16];    max = fmaxf(max, tmp);    smem[tidx] = max;
#endif//(NTHREADS_EB >> 5) == 32
#endif//(NTHREADS_EB >> 5) >= 16
#endif//(NTHREADS_EB >> 5) >=  8
#endif//(NTHREADS_EB >> 5) >=  4
#endif//(NTHREADS_EB >> 5) >=  2
#endif//USE_WARP_SHUFFLE_FUNC_EB
  }/* if( tidx < (NTHREADS_EB >> 5) ){ */
  __syncthreads();


  /** 3. return the resultant maximum */
  val = smem[0];
  /* __syncthreads(); */

  return (val);
}


/**
 * @fn getMaximumFloatGrid
 *
 * @brief Get maximum value within a grid.
 */
__device__ __forceinline__ float getMaximumFloatGrid
  (float val,
   const int bnum, const int bidx, const int tidx, const int lane,
   volatile float * RESTRICT smem, volatile float * RESTRICT gmem, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  val = getMaximumFloatBlock(val, tidx, lane, smem);

  /** global reduction is necessary only if number of blocks is greater than unity */
  if( bnum > 1 ){
    /** share local maximum via global memory */
    /** store data on the global memory */
    if( tidx == 0 )
      gmem[bidx] = val;

    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);


    /** get maximum by a representative block */
    if( bidx == 0 ){
      const int Nloop = BLOCKSIZE(bnum, NTHREADS_EB);
      for(int loop = 0; loop < Nloop; loop++){
	const int target = tidx + loop * NTHREADS_EB;

	/** load from the global memory */
	float tmp = ((target < bnum) ? (gmem[target]) : (val));

	/** get maximum */
	val = getMaximumFloatBlock(tmp, tidx, lane, smem);
      }/* for(int loop = 0; loop < Nloop; loop++){ */

      if( tidx == 0 )
	gmem[0] = val;
    }/* if( bidx == 0 ){ */


    /** global synchronization within bnum blocks */
    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** load from the global memory */
    if( tidx == 0 )
      smem[0] = gmem[0];
    __syncthreads();

    val = smem[0];
    /* __syncthreads(); */
  }/* if( bnum > 1 ){ */

  return (val);
}



/**# この関数は，各 tree node の中で重心から一番遠い粒子を選び出すことができる． */
/**# 従って，中心点を1つ決めた時に，その位置から最も遠い粒子を選び出すアルゴリズムに書き換えることも容易なはずで，tree 構造を使うので全探索しなくてもすむようになる． */
/**# つまり，findFurthestParticle の置き換えとして使う． */
/**# また，TSUB_MAC を排除した方が高速な可能性があるので，そういうことにも気をつけながらコードを読み直してみる． */

/**# ある特定の particle を検出すれば良いだけなので，single block で探索するコードにすれば良いと思う．<-- この方針で実装する．つまり，全体同期関数は不要 */

/**
 * @fn calcMultipole_kernel
 *
 * @brief Calculate multipole moment(s) of tree cells based on the width-first search.
 *
 * @param (bottomLev) bottom level of the PH-hierarchy
 * @param (level) head index and number of tree cells contained in the corresponding hierarchy
 * @param (cell) head index and number of N-body particles contained in the corresponding tree cell
 * @param (leaf) a flag to remember the corresponding tree cell is either leaf(true) of node(false)
 * @param (pi) position and mass of N-body particles
 * @param (node) head index and number of pseudo particles of the corresponding tree cell (a leaf cell contains up to NCRIT particles)
 * @param (more) head index and number of child particles of the corresponding j-particle (a leaf cell contains up to NCRIT particles)
 * @param (node2cell) index of the tree cell corresponding a pseudo particle
 * @return (pj) position and squared radius of pseudo N-body particle as j-particles
 * @return (mj) mass of pseudo N-body particle as j-particles
 * @return (mr2) mass times r squared of pseudo N-body particle as j-particles
 * @return (bmax) size of pseudo N-body particle as j-particles
 */
#define NTHREADS_MAKE_INC NTHREADS_EB
#include "../tree/make_inc.cu"
__global__ void findFurthestParticle_kernel
     (const int bottomLev, READ_ONLY PHinfo * RESTRICT level,
      READ_ONLY treecell * RESTRICT cell, READ_ONLY bool * RESTRICT leaf, READ_ONLY position * RESTRICT pi,
      READ_ONLY uint * RESTRICT node, READ_ONLY uint * RESTRICT more, READ_ONLY int * RESTRICT node2cell,
      jparticle * RESTRICT pj, jmass * RESTRICT mj, real * RESTRICT bmax,
      int * RESTRICT more0Buf, int * RESTRICT more1Buf, real * RESTRICT rjmaxBuf, int * RESTRICT overflow
      )
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;

  const int lane = tidx & (warpSize - 1);
/*   const int head = tidx - lane; */
/* #ifndef USE_WARP_SHUFFLE_FUNC_EB */
  const int tail = NTHREADS_EB - 1;
/* #endif//USE_WARP_SHUFFLE_FUNC_EB */
/*   const int hbuf = DIV_TSUB_EB(head) * TSUB_EB * NBUF_EB;/\**< head index of the shared array close and queue within a thread group *\/ */

/* #ifdef  USE_WARP_SHUFFLE_FUNC_EB */
/*   int smem; */
/* #else///USE_WARP_SHUFFLE_FUNC_EB */
  __shared__  int smem[NTHREADS_EB];
/* #endif//USE_WARP_SHUFFLE_FUNC_EB */

  __shared__ jparticle pj_sm[NTHREADS_EB];
  __shared__      real rjbuf[NTHREADS_EB * NBUF_EB];
  __shared__  int      list0[NTHREADS_EB * NBUF_EB];
  __shared__  int      list1[NTHREADS_EB * NBUF_EB];
  __shared__  int      pjidx[NTHREADS_EB * NBUF_EB];  /**< pjidx would be able to be removed by integrating with list; however, reducing SM usage from 24KB to 20KB has no effect @ Ttot = 512 */

  /** head index of remote buffers */
  const int bufHead = (DIV_TSUB_EB(head) + bidx * NGROUPS_EB) * NUM_ALLOC_EBBUF;







  const      int cidx = 0;
  const treecell root = cell[cidx];


  /* jcom.x, jcom.y, jcom.z からの距離が最大の particle を pick up する実装になっているはずなので，その部分を抽出する． */
  /* TSUB_EB threads が 1 個の tree node を担当するような実装になっているので，その部分を改変してあげる． */
  float4 fcen = {0.0f, 0.0f, 0.0f, 0.0f};


  /** estimate size of particle distribution */
  /** initialize list of examined tree nodes and related variables */
  int inum = root.num;
  int Ntry = 1;

  /* list0, list1 は shared memory に置いておくとして，global memory 上に退避するための配列も作っておいて，適宜やり取りしてあげるべきだと思う． */
  if( tidx == 0 )
    list0_sm[0] = cidx;


  __syncthreads();


  /** pick up NTHREADS_EB i-particles in maximum to estimate bmax */
  /* ここを NTHREADS_EB 以下にすることで，block 内の reduction に持ち込めるようにする */
  while( inum > NTHREADS_EB ){
    float fmin = 0.0f;
    int Nloc = 0;
    int Nbuf = 0;


    int Niter = BLOCKSIZE(Ntry, NTHREADS_EB * NBUF_EB * bnum);
    /* head = NTHREADS_EB * NBUF_EB * bidx からの sweep を Niter 回繰り返した後で reduction すれば良い */
    for(int iter = 0; iter < Niter; iter++){
      const int Nsweep = (Ntry > (NTHREADS_EB * NBUF_EB * bnum)) ? (NTHREADS_EB * NBUF_EB * bsub) : Ntry;


      /** set number of working blocks by using Ntry */
      const int bsub = BLOCKSIZE(Nsweep, NTHREADS_EB);


      /** get number of elements for this block */
      const int Nini = (Nsweep *  bidx     ) / bsub;
      const int Nfin = (Nsweep * (bidx + 1)) / bsub;
      const int Nget = Nfin - Nini;


      /** load data on list0_gm to list0 (on shared memory) */
      const int ibuf_loop = BLOCKSIZE(Nget, NTHREADS_EB);
      for(int ibuf = 0; ibuf < ibuf_loop; ibuf++)
	list0[tidx + ibuf * NTHREADS_EB] = list0_gm[Nini + tidx + ibuf * NTHREADS_EB];


      for(int ibuf = 0; ibuf < ibuf_loop; ibuf++){
	cnum = 0;
	if( (tidx + ibuf * NTHREADS_EB) < Nget ){
	  /** load a tree node corresponding the tree cell */
	  more_tmp = node[list0[tidx + ibuf * NTHREADS_EB]];
	  const int nodenum  = 1 + (more_tmp >> IDXBITS);
	  const int nodehead =      more_tmp  & IDXMASK;

	  /** load all child nodes of the tree cell */
	  more_tmp = more[nodehead];
	  cnum  = 1 + (more_tmp >> IDXBITS);
	  chead =      more_tmp  & IDXMASK;
	  for(int jj = 1; jj < nodenum; jj++)
	    cnum += (1 + (more[nodehead + jj] >> IDXBITS));
	}/* if( (tidx + ibuf * NTHREADS_EB) < Nget ){ */


	/** scan within a block */
	PREFIX_SUM_MAKE_INC_BLCK(cnum, tidx, lane, smem);
	const int lend = BLOCKSIZE(smem[tail], NBUF_EB * NTHREADS_EB);


	/** upload loaded child nodes to pjidx */
	for(int ll = 0; ll < lend; ll++){
	  const int unum =
	    (smem[tidx] <= (NTHREADS_EB * NBUF_EB)) ? cnum :
	    ((smem[tidx] >= (cnum + NTHREADS_EB * NBUF_EB)) ? (0) : (cnum + NTHREADS_EB * NBUF_EB - smem[tidx]));
	  const int shead = smem[tidx] - cnum;

	  for(int jj = 0; jj < unum; jj++){
	    pjidx[shead + jj] = chead;
	    chead++;
	  }/* for(int jj = 0; jj < unum; jj++){ */


	  cnum -= unum;
	  const int Ntmp = smem[tidx] - (NBUF_EB * NTHREADS_EB);/**< Ntmp is a temporal buffer */


	  /** pick up candidate tree nodes */
#   if  NBUF_EB == 4
	  alignedFlt fjmax_loc = { FLT_MIN ,  FLT_MIN ,  FLT_MIN ,  FLT_MIN };
	  alignedInt pjidx_loc = {NULL_NODE, NULL_NODE, NULL_NODE, NULL_NODE};
#endif//NBUF_EB == 4
#   if  NBUF_EB == 2
	  alignedFlt fjmax_loc = { FLT_MIN ,  FLT_MIN };
	  alignedInt pjidx_loc = {NULL_NODE, NULL_NODE};
#endif//NBUF_EB == 2


	  const int stail = (smem[tail] < (NTHREADS_EB * NBUF_EB)) ? (smem[tail]) : (NTHREADS_EB * NBUF_EB);


#pragma unroll
	  for(int kk = 0; kk < NBUF_EB; kk++){
	    const int jj = tidx + kk * NTHREADS_EB;
	    if( jj >= stail )
	      break;

	    const int kidx = pjidx[jj];
	    const jparticle jpos = pj[kidx];

	    const float dx = CAST_R2F(jpos.x) - fcen.x;
	    const float dy = CAST_R2F(jpos.y) - fcen.y;
	    const float dz = CAST_R2F(jpos.z) - fcen.z;
	    const float d2 = FLT_MIN + dx * dx + dy * dy + dz * dz;
	    const float dr = d2 * rsqrtf(d2);

	    const float fjmax = CAST_R2F(bmax[kidx]) + dr;
	    const float fjmin = -rjmax + (2.0f * (1.0f - FLT_EPSILON)) * dr;
	    fmin = fmaxf(rjmin, rmin);

	    if( fjmax > fmin ){
	      pjidx_loc.ia[kk] = kidx;
	      fjmax_loc.fa[kk] = fjmax;
	    }/* if( fjmax > fmin ){ */
	  }/* for(int kk = 0; kk < NBUF_EB; kk++){ */


	  /** share fmin within the grid */
	  fmin = getMaximumFloatGrid(fmin, bsub, bidx, tidx, lane, (float *)smem, gmem, gsync_sub0, gsync_sub1);







	  /** recheck local buffer (is really rjmax greater than rmin ?) */
#pragma unroll
	  for(int jj = 0; jj < NBUF_EB; jj++){
	    const int share = ( fjmax_loc.ra[jj] > fmin ) ? 1 : 0;
	    PREFIX_SUM_MAKE_INC_BLCK(share, tidx, lane, smem);

	    if( share ){
	      const int dst = Nloc + smem[tidx] - 1;
	      list1[dst] = pjidx_loc.ia[jj];
	      fjbuf[dst] = rjmax_loc.ra[jj];
	    }/* if( share ){ */
	    Nloc += smem[tail];


	    /* あと，この -1 の出所がよくわからない． */
	    if( Nloc > ((NBUF_EB - 1) * NTHREADS_EB) ){
	      for(int kk = tidx; kk < Nloc; kk += NTHREADS_EB){
		/* bufHead をうまくずらしておくか，きちんと parallel prefix sum を計算して，仕切りを入手しておく． */
		/* どう考えても，parallel prefix sum among blocks の出番 */
		more1Buf[bufHead + Nbuf + kk] = list1[hbuf + kk];
		fjmaxBuf[bufHead + Nbuf + kk] = fjbuf[hbuf + kk];
	      }/* for(int kk = lane; kk < Nloc; kk += TSUB_EB){ */

	      Nbuf += Nloc;
	      Nloc = 0;
	    }/* if( Nloc > ((NBUF_EB - 1) * TSUB_EB) ){ */
	  }/* for(int jj = 0; jj < NBUF_EB; jj++){ */

#ifdef  USE_WARP_SHUFFLE_FUNC_EB
	  smem         = Ntmp;/**< Ntmp is a temporal buffer */
#else///USE_WARP_SHUFFLE_FUNC_EB
	  smem[tidx].i = Ntmp;/**< Ntmp is a temporal buffer */
#endif//USE_WARP_SHUFFLE_FUNC_EB
	}/* for(int ll = 0; ll < lend; ll++){ */
      }/* for(int ibuf = 0; ibuf < NBUF_EB; ibuf++){ */

      Ntry -= Nsweep;

      /** copy data from global memory to shared memory */
      const int Ncopy = (Ntry < (TSUB_EB * NBUF_EB)) ? (Ntry) : (TSUB_EB * NBUF_EB);
      for(int jj = lane; jj < Ncopy; jj += TSUB_EB)
	list0[hbuf + jj] = more0Buf[(bufHead + (TSUB_EB * NBUF_EB) * (iter + 1)) + jj];
    }/* for(int iter = 0; iter < Niter; iter++){ */

    if( Nbuf != 0 ){
      for(int ll = lane; ll < Nloc; ll += TSUB_EB){
	more1Buf[bufHead + Nbuf + ll] = list1[hbuf + ll];
	rjmaxBuf[bufHead + Nbuf + ll] = rjbuf[hbuf + ll];
      }/* for(int ll = lane; ll < Nloc; ll += TSUB_EB){ */

      for(int ll = lane; ll < TSUB_EB * NBUF_EB; ll += TSUB_EB){
	list1[hbuf + ll] = more1Buf[bufHead + ll];
	rjbuf[hbuf + ll] = rjmaxBuf[bufHead + ll];
      }/* for(int ll = lane; ll < TSUB_EB * NBUF_EB; ll += TSUB_EB){ */
    }/* if( Nbuf != 0 ){ */

    Ntry = Nbuf + Nloc;
    if( (lane == 0) && (Ntry > NUM_ALLOC_EBBUF) )
      atomicAdd(overflow, 1);


    /** list up all child nodes that satisfy rjmax > rmin */
    inum = 0;
    Nloc = 0;
    Nbuf = 0;

    Niter = BLOCKSIZE(Ntry, NBUF_EB * TSUB_EB);
    for(int iter = 0; iter < Niter; iter++){
      const int krem = (Ntry < (NBUF_EB * TSUB_EB)) ? Ntry : (NBUF_EB * TSUB_EB);
      const int knum = BLOCKSIZE(krem, TSUB_EB);

      for(int ki = 0; ki < knum; ki++){
	int cellIdx = lane + ki * TSUB_EB;
	int  add = 0;
	int iadd = 0;

	/** select distant tree cells */
	if( cellIdx < krem ){
	  /** when the current node must be taken into account */
	  cellIdx += hbuf;
	  if( rjbuf[cellIdx] > rmin ){
	    /** count up total number of contained i-particles */
	    cellIdx = node2cell[list1[cellIdx]];
	    iadd = cell[cellIdx].num;

	    add = 1;
	  }/* if( rjbuf[cellIdx] > rmin ){ */
	}/* if( cellIdx < krem ){ */

	/** remove duplicated tree cells */
#ifdef  USE_WARP_SHUFFLE_FUNC_EB
	smem = prefixSumTsub(add, lane);
#else///USE_WARP_SHUFFLE_FUNC_EB
	prefixSumTsub(add, smem, tidx, lane);
#endif//USE_WARP_SHUFFLE_FUNC_EB

	if( add ){
	  /** test uploading... */
#ifdef  USE_WARP_SHUFFLE_FUNC_EB
	  const int smidx = Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_EB
	  const int smidx = Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_EB
	  list0[hbuf + smidx] = cellIdx;

	  /** if detect duplication, upload flag is turned off */
	  if( ((smidx > 0) && (list0[hbuf + smidx - 1] == cellIdx)) || ((Nbuf > 0) && (smidx == 0) && (more0Buf[bufHead + Nbuf - 1] == cellIdx)) ){
	    add  = 0;
	    iadd = 0;
	  }/* if( ((smidx > 0) && (list0[hbuf + smidx - 1] == cellIdx)) || ((Nbuf > 0) && (smidx == 0) && (more0Buf[bufHead + Nbuf - 1] == cellIdx)) ){ */
	}/* if( add ){ */

	/** save tree cells on the local buffer */
#ifdef  USE_WARP_SHUFFLE_FUNC_EB
	smem = prefixSumTsub(add, lane);
#else///USE_WARP_SHUFFLE_FUNC_EB
	prefixSumTsub(add, smem, tidx, lane);
#endif//USE_WARP_SHUFFLE_FUNC_EB

	if( add ){
#ifdef  USE_WARP_SHUFFLE_FUNC_EB
	  const int smidx = Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_EB
	  const int smidx = Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_EB
	  list0[hbuf + smidx] = cellIdx;
	}/* if( add ){ */

#ifdef  USE_WARP_SHUFFLE_FUNC_EB
	Nloc += __shfl(smem, TSUB_EB - 1, TSUB_EB);
#else///USE_WARP_SHUFFLE_FUNC_EB
	Nloc +=        smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_EB

	/** move data to the remote buffer if necessary */
	if( Nloc > ((NBUF_EB - 1) * TSUB_EB) ){
	  for(int ll = lane; ll < Nloc; ll += TSUB_EB)
	    more0Buf[bufHead + Nbuf + ll] = list0[hbuf + ll];
	  Nbuf += Nloc;
	  Nloc  = 0;
	}/* if( Nloc > ((NBUF_EB - 1) * TSUB_EB) ){ */

	/** sum up iadd within TSUB_EB threads */
#ifdef  USE_WARP_SHUFFLE_FUNC_EB
	iadd = accumulateIntTsub(iadd);
#else///USE_WARP_SHUFFLE_FUNC_EB
	accumulateIntTsub(&iadd, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_EB
	inum += iadd;
      }/* for(int ki = 0; ki < knum; ki++){ */

      Ntry -= krem;

      /** copy data from remote buffer to local buffer */
      const int Ncopy = (Ntry < (TSUB_EB * NBUF_EB)) ? (Ntry) : (TSUB_EB * NBUF_EB);
      for(int jj = lane; jj < Ncopy; jj += TSUB_EB){
	rjbuf[hbuf + jj] = rjmaxBuf[bufHead + NBUF_EB * TSUB_EB * (iter + 1) + jj];
	list1[hbuf + jj] = more1Buf[bufHead + NBUF_EB * TSUB_EB * (iter + 1) + jj];
      }/* for(int jj = lane; jj < Ncopy; jj += TSUB_EB){ */
    }/* for(int iter = 0; iter < Niter1; iter++){ */

    if( Nbuf != 0 ){
      for(int ll = lane; ll < Nloc; ll += TSUB_EB)
	more0Buf[bufHead + Nbuf + ll] = list0[hbuf + ll];

      for(int ll = lane; ll < NBUF_EB * TSUB_EB; ll += TSUB_EB)
	list0[hbuf + ll] = more0Buf[bufHead + ll];
    }/* if( Nbuf != 0 ){ */

    Ntry = Nloc + Nbuf;
    if( (lane == 0) && (Ntry > NUM_ALLOC_EBBUF) )
      atomicAdd(overflow, 1);


    /* Ntry の値を変化させるので，この段階で全体同期が必要． */
    /* つまり，全員で inum の値を共有してあげる必要がある． */


    globalSync(tidx, bidx, bnum, gsync0, gsync1);


  }/* while( inum > NI_BMAX_ESTIMATE ){ */







  /** check positions of all the pick upped i-particles */
  real jbmax = ZERO;
  /** Since NI_BMAX_ESTIMATE <= TSUB_EB * NBUF_EB, Ntry1 is less than TSUB_EB * NBUF_EB */
  /** load index of the pick upped i-particles to list1 */
  const int Niter = BLOCKSIZE(Ntry, TSUB_EB);
  int Ncand = 0;
  for(int iter = 0; iter < Niter; iter++){
    treecell cand;
    int pnum = 0;
    if( lane < Ntry ){
      cand = cell[list0[hbuf + iter * TSUB_EB + lane]];
      pnum = cand.num;
    }/* if( lane < Ntry ){ */

#ifdef  USE_WARP_SHUFFLE_FUNC_EB
    smem = prefixSumTsub(pnum, lane);
    for(int jj = 0; jj < pnum; jj++)
      list1[hbuf + Ncand + smem         - pnum + jj] = cand.head + jj;
    Ncand += __shfl(smem, TSUB_EB - 1, TSUB_EB);
#else///USE_WARP_SHUFFLE_FUNC_EB
    prefixSumTsub(pnum, smem, tidx, lane);
    for(int jj = 0; jj < pnum; jj++)
      list1[hbuf + Ncand + smem[tidx].i - pnum + jj] = cand.head + jj;
    Ncand +=        smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_EB

    Ntry -= TSUB_EB;
  }/* for(int iter = 0; iter < Niter; iter++){ */

  for(int jj = lane; jj < Ncand; jj += TSUB_EB){
    const position ipos = pi[list1[hbuf + jj]];

    const real dx = ipos.x - jcom.x;
    const real dy = ipos.y - jcom.y;
    const real dz = ipos.z - jcom.z;
    const real r2 = FLT_MIN + dx * dx + dy * dy + dz * dz;

    jbmax = FMAX(jbmax, r2);
  }/* for(int jj = lane; jj < Ncand; jj += TSUB_EB){ */

#ifdef  USE_WARP_SHUFFLE_FUNC_EB
  jbmax = getMaximumRealTsub(jbmax);
#else///USE_WARP_SHUFFLE_FUNC_EB
  getMaximumRealTsub(&jbmax, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_EB
  jbmax *= RSQRT(jbmax);
  if( lane == 0 )
    bmax[jidx] = jbmax;

  jcom.w  = jbmax * jbmax;

  if( lane == 0 )
    pj[jidx] = jcom;



}
