/**
 * @file peano_dev.cu
 *
 * @brief Source code for calculating Peano--Hilbert space-filling curve
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/06/27 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <helper_cuda.h>

#ifdef  CUB_AVAILABLE
#include <cub/device/device_radix_sort.cuh>
#else///CUB_AVAILABLE
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#endif//CUB_AVAILABLE

#include "macro.h"
#include "cudalib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "../util/gsync_dev.cu"

#include "peano.h"
#include "peano_dev.h"

#include "../tree/macutil.h"

#include "../tree/make.h"/**< required to read NLEAF */

#ifndef SERIALIZED_EXECUTION
#include "../para/exchange_dev.h"/**< required to read NBLOCKS_PER_SM_BOX */
#endif//SERIALIZED_EXECUTION


/** in L1 cache preferred configuration, capacity of shared memory is 16KiB per SM */
/** real4 smem[NTHREADS_PH] corresponds 16 * NTHREADS_PH bytes */
#define NBLOCKS_PER_SM_PH (1024 / NTHREADS_PH)

#define REGISTERS_PER_THREAD_PH (37)
/** calcPHkey_kernel uses 32 registers @ Tesla M2090, Ttot = 1024 (registers are spilled to local memory) */
/** calcPHkey_kernel uses 47 registers @ Tesla M2090, Ttot =  512 */
/** calcPHkey_kernel uses 36 registers @ Tesla M2090, Ttot =  128, 256 */
#   if  GPUVER == 20
#undef  REGISTERS_PER_THREAD_PH
#   if  NTHREADS_PH == 1024
#define REGISTERS_PER_THREAD_PH (32)
#else///NTHREADS_PH == 1024
#   if  NTHREADS_PH ==  512
#define REGISTERS_PER_THREAD_PH (47)
#else///NTHREADS_PH ==  512
#define REGISTERS_PER_THREAD_PH (36)
#endif//NTHREADS_PH ==  512
#endif//NTHREADS_PH == 1024
#endif//GPUVER == 20

/** limitation from number of registers */
#   if  NBLOCKS_PER_SM_PH > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_PH * NTHREADS_PH))
#undef  NBLOCKS_PER_SM_PH
#define NBLOCKS_PER_SM_PH   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_PH * NTHREADS_PH))
#endif//NBLOCKS_PER_SM_PH > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_PH * NTHREADS_PH))

/** maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_PH > (MAX_THREADS_PER_SM / NTHREADS_PH)
#undef  NBLOCKS_PER_SM_PH
#define NBLOCKS_PER_SM_PH   (MAX_THREADS_PER_SM / NTHREADS_PH)
#endif//NBLOCKS_PER_SM_PH > (MAX_THREADS_PER_SM / NTHREADS_PH)

/** maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_PH > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_PH
#define NBLOCKS_PER_SM_PH   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_PH > MAX_BLOCKS_PER_SM

/** maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_PH > ((MAX_WARPS_PER_SM * 32) / NTHREADS_PH)
#undef  NBLOCKS_PER_SM_PH
#define NBLOCKS_PER_SM_PH   ((MAX_WARPS_PER_SM * 32) / NTHREADS_PH)
#endif//NBLOCKS_PER_SM_PH > ((MAX_WARPS_PER_SM * 32) / NTHREADS_PH)

/** number of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_PH < 1
#undef  NBLOCKS_PER_SM_PH
#define NBLOCKS_PER_SM_PH  (1)
#endif//NBLOCKS_PER_SM_PH < 1


/**
 * @fn dilate3D
 *
 * @brief Calculate Morton key in 3D.
 * @detail works less equal than 63 (= 3 * 21) bits key
 * based on Raman & Wise (2007), IEEE Trans. Comput., C99, 13
 */
__device__ __forceinline__ PHint dilate3D(const PHint val)
{
  PHint ret = val;

#   if  MAXIMUM_PHKEY_LEVEL > 16
  ret = (ret * 0x100000001) & 0x7fff00000000ffff;/**< ((x << 32) + x) = (2^32 + 1) * x = (16^8 + 1) * x */
#endif//MAXIMUM_PHKEY_LEVEL > 16
  ret = (ret * 0x000010001) & 0x00ff0000ff0000ff;/**< ((x << 16) + x) = (2^16 + 1) * x = (16^4 + 1) * x */
  ret = (ret * 0x000000101) & 0x700f00f00f00f00f;/**< ((x <<  8) + x) = (2^8  + 1) * x = (16^2 + 1) * x */
  ret = (ret * 0x000000011) & 0x30c30c30c30c30c3;/**< ((x <<  4) + x) = (2^4  + 1) * x = (16^1 + 1) * x */
  ret = (ret * 0x000000005) & 0x1249249249249249;/**< ((x <<  2) + x) = (2^2  + 1) * x = (       5) * x */

  return (ret);
}


#define NTHREADS_COMPARE_VEC3_INC NTHREADS_PH
#include "../util/compare_vec3_inc.cu"
/**
 * @fn calcPHkey_kernel
 *
 * @brief Calculate Peano--Hilbert space-filling curve on GPU.
 */
__global__ void __launch_bounds__(NTHREADS_PH, NBLOCKS_PER_SM_PH) calcPHkey_kernel
     (const int num, READ_ONLY position * RESTRICT ipos, float4 * RESTRICT min_all, float4 * RESTRICT max_all,
      const int nlevel, int * RESTRICT idx, PHint * RESTRICT key, int * RESTRICT gsync0, int * RESTRICT gsync1
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
      , position * RESTRICT center_dev
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  SHARE_PH_BOX_BOUNDARY
      , READ_ONLY float4 * RESTRICT box_min, READ_ONLY float4 * RESTRICT box_max
#endif//SHARE_PH_BOX_BOUNDARY
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
      , float * RESTRICT length
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
      )
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int bnum =   GRIDDIM_X1D;

  const int hidx = tidx - (tidx & (warpSize - 1));

  const int ihead = (num *      bidx ) / bnum;
  const int itail = (num * (1 + bidx)) / bnum;


  /** calculate required box size to contain all N-body particles in the local domain */
  float4 min = { FLT_MAX,  FLT_MAX,  FLT_MAX,  FLT_MAX};
  float4 max = {-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX};

  for(int ih = ihead; ih < itail; ih += NTHREADS_PH){
    const int ii = ih + tidx;
    if( ii < itail ){
      const position pi = ipos[ii];
      const float px = CAST_R2F(pi.x);
      const float py = CAST_R2F(pi.y);
      const float pz = CAST_R2F(pi.z);

      min.x = fminf(min.x, px);      max.x = fmaxf(max.x, px);
      min.y = fminf(min.y, py);      max.y = fmaxf(max.y, py);
      min.z = fminf(min.z, pz);      max.z = fmaxf(max.z, pz);
    }/* if( ii < itail ){ */
  }/* for(int ih = ihead; ih < itail; ih += NTHREADS_PH){ */


  /** get the box size within a warp */
  __shared__ float4 smem[NTHREADS_PH];
  GET_VEC3_MIN_MAX_GRID(&min, &max, smem, tidx, hidx, min_all, max_all, bidx, bnum, gsync0, gsync1);


#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  if( (tidx + bidx) == 0 ){
    position center;
    center.x = CAST_F2R(HALF * (min.x + max.x));
    center.y = CAST_F2R(HALF * (min.y + max.y));
    center.z = CAST_F2R(HALF * (min.z + max.z));
    center.m = ZERO;

    *center_dev = center;
  }/* if( (tidx + bidx) == 0 ){ */
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR


#ifdef  SHARE_PH_BOX_BOUNDARY
  float4 tmp;
  tmp = *box_min;  min.x = fminf(min.x, tmp.x);  min.y = fminf(min.y, tmp.y);  min.z = fminf(min.z, tmp.z);
  tmp = *box_max;  max.x = fmaxf(max.x, tmp.x);  max.y = fmaxf(max.y, tmp.y);  max.z = fmaxf(max.z, tmp.z);
#endif//SHARE_PH_BOX_BOUNDARY


  /** set box */
  float diameter = 0.0f;
  diameter = fmaxf(diameter, max.x - min.x);
  diameter = fmaxf(diameter, max.y - min.y);
  diameter = fmaxf(diameter, max.z - min.z);
  diameter = ldexpf(1.0f, (int)ceilf(log2f(diameter)));
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  /** store diameter to global memory (purpose: detect external particles) */
  if( (bidx + tidx) == 0 )
    *length = diameter;
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  const float dinv = 1.0f / diameter;

  min.x = 0.5f * (min.x + max.x - diameter);
  min.y = 0.5f * (min.y + max.y - diameter);
  min.z = 0.5f * (min.z + max.z - diameter);


  /** calculate Peano--Hilbert key */
  const PHint keymax = (PHint)1 << nlevel;
  const float dscale = dinv * (float)keymax;

  for(int ih = ihead; ih < itail; ih += NTHREADS_PH){
    const int ii = ih + tidx;
    if( ii < itail ){
      const position pi = ipos[ii];
      PHint px = (PHint)(dscale * (CAST_R2F(pi.x) - min.x));
      PHint py = (PHint)(dscale * (CAST_R2F(pi.y) - min.y));
      PHint pz = (PHint)(dscale * (CAST_R2F(pi.z) - min.z));
      PHint tkey = 0;

      for(int jj = nlevel - 1; jj >= 0; jj--){
	/** get xi, yi, zi from given position */
	PHint xi = (px >> jj) & 1;
	PHint yi = (py >> jj) & 1;
	PHint zi = (pz >> jj) & 1;

	/** turn px, py, and pz */
	px ^= -( xi & ((!yi) |   zi));
	py ^= -((xi & (  yi  |   zi)) | (yi & (!zi)));
	pz ^= -((xi &  (!yi) & (!zi)) | (yi & (!zi)));

	/** append 3 bits to the key */
	tkey |= ((xi << 2) | ((xi ^ yi) << 1) | ((xi ^ zi) ^ yi)) << (3 * jj);

	/** if zi == 1, then rotate uncyclic (x->z->y->x) */
	if( zi ){	    PHint pt = px;	    px = py;	    py = pz;	  pz = pt;	}
	else{
	  /** if yi == 0, then exchange x and z */
	  if( !yi ){	    PHint pt = px;	    px = pz;	                  pz = pt;	}
	}
      }/* for(int jj = nlevel - 1; jj >= 0; jj--){ */

      idx[ii] = ii;
      key[ii] = tkey;
    }/* if( ii < itail ){ */
  }/* for(int ih = ihead; ih < itail; ih += bnum * NTHREADS_PH){ */
}


/**
 * @fn allocPeanoHilbertKey_dev
 *
 * @brief Memory allocation for PH-key.
 */
extern "C"
muse allocPeanoHilbertKey_dev
(const int num, int **idx_dev, PHint **key_dev, PHint **key_hst, float4 **minall, float4 **maxall, int **gsync0, int **gsync1,
#ifndef SERIALIZED_EXECUTION
 float4 **box_min, float4 **box_max, float4 **min_hst, float4 **max_hst,
#endif//SERIALIZED_EXECUTION
 soaPHsort *dev, soaPHsort *hst,
#ifdef  CUB_AVAILABLE
 soaPHsort *pre, void **temp_storage, int **idx_pre, PHint **key_pre,
#endif//CUB_AVAILABLE
 const deviceProp devProp)
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  /** Peano--Hilbert key of N-body particles */
  /** the size of the array is set to be a multiple of NTHREADS_PH */
  size_t size = (size_t)num;
  if( (num % NTHREADS_PH) != 0 )
    size += (size_t)(NTHREADS_PH - (num % NTHREADS_PH));

  /** memory allocation and simple confirmation */
  mycudaMalloc    ((void **)idx_dev, num * sizeof(  int));  alloc.device += num * sizeof(  int);
  mycudaMalloc    ((void **)key_dev, num * sizeof(PHint));  alloc.device += num * sizeof(PHint);
  mycudaMallocHost((void **)key_hst, num * sizeof(PHint));  alloc.host   += num * sizeof(PHint);

#ifdef  CUB_AVAILABLE
  mycudaMalloc    ((void **)idx_pre, num * sizeof(  int));  alloc.device += num * sizeof(  int);
  mycudaMalloc    ((void **)key_pre, num * sizeof(PHint));  alloc.device += num * sizeof(PHint);
#endif//CUB_AVAILABLE

#ifdef  SERIALIZED_EXECUTION
  const size_t num_ph = devProp.numSM * NBLOCKS_PER_SM_PH;
#else///SERIALIZED_EXECUTION
  const size_t num_ph = devProp.numSM * ((NBLOCKS_PER_SM_PH > NBLOCKS_PER_SM_BOX) ? (NBLOCKS_PER_SM_PH) : (NBLOCKS_PER_SM_BOX));
  mycudaMalloc    ((void **)box_min, sizeof(float4));  alloc.device += sizeof(float4);  dev->box_min = *box_min;
  mycudaMalloc    ((void **)box_max, sizeof(float4));  alloc.device += sizeof(float4);  dev->box_max = *box_max;
  mycudaMallocHost((void **)min_hst, sizeof(float4));  alloc.host   += sizeof(float4);  dev->min_hst = *min_hst;
  mycudaMallocHost((void **)max_hst, sizeof(float4));  alloc.host   += sizeof(float4);  dev->max_hst = *max_hst;
#endif//SERIALIZED_EXECUTION
  mycudaMalloc((void **)minall, num_ph * sizeof(float4));  alloc.device += num_ph * sizeof(float4);
  mycudaMalloc((void **)maxall, num_ph * sizeof(float4));  alloc.device += num_ph * sizeof(float4);
  mycudaMalloc((void **)gsync0, num_ph * sizeof(int   ));  alloc.device += num_ph * sizeof(int);
  mycudaMalloc((void **)gsync1, num_ph * sizeof(int   ));  alloc.device += num_ph * sizeof(int);

  initGsync_kernel<<<1, num_ph>>>(num_ph, *gsync0, *gsync1);
  getLastCudaError("initGsync_kernel");

  dev->idx = *idx_dev;
  dev->key = *key_dev;  hst->key = *key_hst;
  dev->min = *minall;
  dev->max = *maxall;
  dev->gsync0 = *gsync0;
  dev->gsync1 = *gsync1;

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


  /** error checking before running the kernel */
  struct cudaFuncAttributes funcAttr;
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, calcPHkey_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_PH ){
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in calcPHkey_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_PH);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_PH = %d.\n", GPUGEN, GPUVER, NTHREADS_PH);
    fflush (stderr);
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
    __KILL__(stderr, "ERROR: # of blocks per SM for calcPHkey_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_PH defined in src/sort/peano_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_PH, funcAttr.numRegs, REGISTERS_PER_THREAD_PH, GPUGEN, GPUVER);
  }/* if( Nblck != NBLOCKS_PER_SM ){ */

  if( (devProp.numSM * NBLOCKS_PER_SM_PH) > NTHREADS_PH ){
    __KILL__(stderr, "ERROR: product (%d) of devProp.numSM(%d) * NBLOCKS_PER_SM_PH(%d) must be smaller than NTHREADS_PH(%d) to use shared memory.\n", devProp.numSM * NBLOCKS_PER_SM_PH, devProp.numSM, NBLOCKS_PER_SM_PH, NTHREADS_PH);
  }/* if( (devProp.numSM * NBLOCKS_PER_SM_PH) > NTHREADS_PH ){ */


  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freePeanoHilbertKey_dev
 *
 * @brief Memory deallocation for PH-key.
 */
extern "C"
void  freePeanoHilbertKey_dev
(int  *idx_dev, PHint  *key_dev, PHint  *key_hst, float4  *minall, float4  *maxall, int  *gsync0, int  *gsync1
#ifndef SERIALIZED_EXECUTION
 , float4  *box_min, float4  *box_max, float4  *min_hst, float4  *max_hst
#endif//SERIALIZED_EXECUTION
#ifdef  CUB_AVAILABLE
 , void  *temp_storage, int  *idx_pre, PHint  *key_pre
#endif//CUB_AVAILABLE
 )
{
  __NOTE__("%s\n", "start");


  mycudaFree    (idx_dev);
  mycudaFree    (key_dev);
  mycudaFreeHost(key_hst);

  mycudaFree(minall);
  mycudaFree(maxall);
  mycudaFree(gsync0);
  mycudaFree(gsync1);

#ifndef SERIALIZED_EXECUTION
  mycudaFree    (box_min);
  mycudaFree    (box_max);
  mycudaFreeHost(min_hst);
  mycudaFreeHost(max_hst);
#endif//SERIALIZED_EXECUTION

#ifdef  CUB_AVAILABLE
  mycudaFree(temp_storage);
  mycudaFree(idx_pre);
  mycudaFree(key_pre);
#endif//CUB_AVAILABLE


  __NOTE__("%s\n", "end");
}


/**
 * @fn sortParticlesPHcurve_kernel
 *
 * @brief Sort N-body particles by PH-key.
 */
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
  const int ii = GLOBALIDX_X1D;

  if( ii < num ){
    /** load old tag */
    const int jj = old[ii];

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
  }/* if( ii < num ){ */
}

/**
 * @fn sortParticlesPHcurve_kernel_offset
 *
 * @brief Sort N-body particles by PH-key.
 */
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
  const int ii = offset + GLOBALIDX_X1D;

  if( ii < num ){
    /** load old tag */
    const int jj = old[ii];

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
  }/* if( ii < num ){ */
}


/**
 * @fn sortParticlesPHcurve_dev
 *
 * @brief Sort N-body particles by PH-key.
 *
 * @sa calcPHkey_kernel
 * @sa sortParticlesPHcurve_kernel
 * @sa sortParticlesPHcurve_kernel_offset
 */
extern "C"
void sortParticlesPHcurve_dev(const int num, iparticle * RESTRICT src, iparticle * RESTRICT dst, soaPHsort dev, const deviceProp devProp
#ifdef  CUB_AVAILABLE
			      , soaPHsort pre
#endif//CUB_AVAILABLE
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
			      , domainLocation *location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
			      , struct timespec *start
#ifdef EXEC_BENCHMARK
			      , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
)
{
  __NOTE__("%s\n", "start");


#ifndef EXEC_BENCHMARK
  checkCudaErrors(cudaDeviceSynchronize());
  clock_gettime(CLOCK_MONOTONIC_RAW, start);
#else///EXEC_BENCHMARK
  initStopwatch();
  *start = _benchIni;
#endif//EXEC_BENCHMARK


  /** generate Peano--Hilbert key */
#ifdef  CUB_AVAILABLE
  calcPHkey_kernel<<<devProp.numSM * NBLOCKS_PER_SM_PH, NTHREADS_PH>>>(num, (*src).pos, dev.min, dev.max, MAXIMUM_PHKEY_LEVEL, pre.idx, pre.key, dev.gsync0, dev.gsync1
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
								       , (*dst).encBall
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  SHARE_PH_BOX_BOUNDARY
								       , dev.box_min, dev.box_max
#endif//SHARE_PH_BOX_BOUNDARY
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
								       , (*location).diameter_dev
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
								       );
#else///CUB_AVAILABLE
  calcPHkey_kernel<<<devProp.numSM * NBLOCKS_PER_SM_PH, NTHREADS_PH>>>(num, (*src).pos, dev.min, dev.max, MAXIMUM_PHKEY_LEVEL, dev.idx, dev.key, dev.gsync0, dev.gsync1
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
								       , (*dst).encBall
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  SHARE_PH_BOX_BOUNDARY
								       , dev.box_min, dev.box_max
#endif//SHARE_PH_BOX_BOUNDARY
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
								       , (*location).diameter_dev
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
								       );
#endif//CUB_AVAILABLE
  getLastCudaError("calcPHkey_kernel");

#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  checkCudaErrors(cudaMemcpy((*dst).encBall_hst, (*dst).encBall, sizeof(position), cudaMemcpyDeviceToHost));
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR

#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  checkCudaErrors(cudaMemcpy((*location).diameter_hst, (*location).diameter_dev, sizeof(float), cudaMemcpyDeviceToHost));
#ifdef  TIME_BASED_MODIFICATION
  location->linv = 2.0f / (*(*location).diameter_hst);
#else///TIME_BASED_MODIFICATION
  location->dL_L = location->eta * location->eps / (*(*location).diameter_hst);
#endif//TIME_BASED_MODIFICATION
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->genPHkey_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER


  /** sort Peano--Hilbert key */
#ifdef  CUB_AVAILABLE
  cub::DeviceRadixSort::SortPairs(dev.temp_storage, dev.temp_storage_size, pre.key, dev.key, pre.idx, dev.idx, num);
#else///CUB_AVAILABLE
  thrust::stable_sort_by_key((thrust::device_ptr<PHint>)dev.key, (thrust::device_ptr<PHint>)(dev.key + num), (thrust::device_ptr<int>)dev.idx);
#endif//CUB_AVAILABLE

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->rsortKey_library));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER


  /** sort N-body particles using Peano--Hilbert key */
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
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;
    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

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

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("sortParticlesPHcurve_kernel");

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->sortBody_kernel));
#endif//HUNT_MAKE_PARAMETER


  /** swap the list structure */
  iparticle _tmp;
  _tmp = *src;
  *src = *dst;
  *dst = _tmp;


#ifdef  EXEC_BENCHMARK
#ifndef HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->sortParticlesPHcurve));
#else///HUNT_MAKE_PARAMETER
  elapsed->sortParticlesPHcurve += elapsed->genPHkey_kernel + elapsed->rsortKey_library + elapsed->sortBody_kernel;
#endif//HUNT_MAKE_PARAMETER
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}


#ifdef  BLOCK_TIME_STEP
/**
 * @fn copySortedParticles_kernel
 *
 * @brief Copy sorted N-body particles.
 */
__global__ void copySortedParticles_kernel(const int Ni, READ_ONLY position * RESTRICT src, position * RESTRICT dst)
{
  const int gidx = GLOBALIDX_X1D;

  if( gidx < Ni )
    dst[gidx] = src[gidx];
}


/**
 * @fn copySortedParticles_dev
 *
 * @brief Copy sorted N-body particles.
 *
 * @sa copySortedParticles_kernel
 */
extern "C"
void copySortedParticles_dev(const int Ni, const iparticle pi)
{
  __NOTE__("%s\n", "start");

  int Nrem = BLOCKSIZE(Ni, 1024);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    copySortedParticles_kernel<<<Nrem, 1024>>>(Ni, pi.pos, pi.jpos);
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * 1024;
      copySortedParticles_kernel<<<Nblck, 1024>>>(Nsub, &pi.pos[hidx], &pi.jpos[hidx]);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */
  getLastCudaError("copySortedParticles_dev");

  __NOTE__("%s\n", "end");
}
#endif//BLOCK_TIME_STEP


/**
 * @fn initPHinfo_kernel
 *
 * @brief Initialize properties of PH-key of tree cells.
 * @detail NOTE: syncthreads and if statements of ii < NUM_PHKEY_LEVEL can be removed; however, remained for safety
 */
__global__ void initPHinfo_kernel(PHinfo *gm)
{
  const int ii = THREADIDX_X1D;
  __shared__ PHinfo sm[NUM_PHKEY_LEVEL];

  /** initialize fundamental properties */
  if( ii < NUM_PHKEY_LEVEL ){
    sm[ii].level = MAXIMUM_PHKEY_LEVEL - ii;

    ulong ntmp = 1;
    int jj = NUM_PHKEY_LEVEL - 1;
    if( ii == 0 )
      while( true ){
	sm[jj].nmax = (int)ntmp;

	ntmp *= NLEAF;    jj--;
	if( ntmp > INT_MAX )	ntmp = INT_MAX;

	if( jj < 0 )	break;
      }/* while( true ){ */
  }/* if( ii < NUM_PHKEY_LEVEL ){ */

  __syncthreads();

  if( ii < NUM_PHKEY_LEVEL )
    gm[ii] = sm[ii];
}


/**
 * @fn initPHinfo_dev
 *
 * @brief Initialize properties of PH-key of tree cells.
 */
extern "C"
void initPHinfo_dev(PHinfo *info)
{
  __NOTE__("%s\n", "start");

  initPHinfo_kernel<<<1, NUM_PHKEY_LEVEL>>>(info);
  getLastCudaError("initPHinfo_kernel");

  __NOTE__("%s\n", "end");
}
