/**
 * @file geo_dev.cu
 *
 * @brief Source code for generating enclosing ball
 * whose center is the geometric on of the enclosing rectangular cuboid
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/29 (Wed)
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
#include <sys/time.h>
#include <helper_cuda.h>

#include <thrust/reduce.h>
#include <thrust/device_ptr.h>

#include "macro.h"
#include "cudalib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "make.h"
#include "walk_dev.h"
#include "geo_dev.h"


/**
 * @fn init_r2max_kernel
 *
 * @brief Initialize array for r2.
 */
__global__ void init_r2max_kernel(const int num, real *r2)
{
  const int tidx = THREADIDX_X1D;

  if( tidx < num )
    r2[tidx] = ZERO;
}

/**
 * @fn init_amin_kernel
 *
 * @brief Initialize array for amin.
 */
__global__ void init_amin_kernel(const int num, real *amin)
{
  const int tidx = THREADIDX_X1D;

  if( tidx < num )
    amin[tidx] = REAL_MAX;
}


#ifdef  USE_WARP_SHUFFLE_FUNC
#define USE_WARP_SHUFFLE_FUNC_COMPARE_INC
#endif//USE_WARP_SHUFFLE_FUNC
#define NTHREADS_COMPARE_INC NTHREADS
#include "../util/compare_inc.cu"


__device__ __forceinline__ float atomicMin(float* addr, float val)
{
  int* addr_as_i = (int*)addr;
  int old = *addr_as_i, assumed;

  do{
    assumed = old;
    old = atomicCAS(addr_as_i, assumed, __float_as_int(fminf(val, __int_as_float(assumed))));
  } while( assumed != old );

  return (__int_as_float(old));
}
__device__ __forceinline__ float atomicMax(float* addr, float val)
{
  int* addr_as_i = (int*)addr;
  int old = *addr_as_i, assumed;

  do{
    assumed = old;
    old = atomicCAS(addr_as_i, assumed, __float_as_int(fmaxf(val, __int_as_float(assumed))));
  } while( assumed != old );

  return (__int_as_float(old));
}


/**
 * @fn calc_r2max_kernel
 *
 * @brief Get r2max in the array.
 */
__global__ void calc_r2max_kernel
(const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo, READ_ONLY position * RESTRICT cen_dev, position * RESTRICT jpos, real * RESTRICT r2_dev, const bool singleCall)
{
  const int tidx = THREADIDX_X1D;
  const int lane    = tidx          & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);

  __shared__ real smem[NTHREADS];

  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];

  real r2 = ZERO;
  position cen = *cen_dev;

  if( lane < info.num ){
    const int idx = info.head + lane;
    const position pi = jpos[idx];

    const real dx = pi.x - cen.x;
    const real dy = pi.y - cen.y;
    const real dz = pi.z - cen.z;
    r2 = FLT_MIN + dx * dx + dy * dy + dz * dz;
  }/* if( lane < info.num ){ */

  r2 = GET_MAX_BLCK(r2, smem, tidx, tidx - lane);
  if( tidx == 0 ){
    if( singleCall )      r2_dev[BLOCKIDX_X1D] = r2;
    else      atomicMax(&r2_dev[BLOCKIDX_X1D], r2);
  }/* if( tidx == 0 ){ */
}


#ifdef  GADGET_MAC
/**
 * @fn calc_amin_kernel
 *
 * @brief Get amin in the array.
 */
__global__ void calc_amin_kernel
(const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo, acceleration * RESTRICT iacc, real * RESTRICT amin_dev, const bool singleCall)
{
  const int tidx = THREADIDX_X1D;
  const int lane    = tidx          & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);

  __shared__ real smem[NTHREADS];

  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];

  real a2 = REAL_MAX;

  if( lane < info.num ){
    const int idx = info.head + lane;
    const acceleration ai = iacc[idx];

    a2 = FLT_MIN + ai.x * ai.x + ai.y * ai.y + ai.z * ai.z;
  }/* if( lane < info.num ){ */

  a2 = GET_MIN_BLCK(a2, smem, tidx, tidx - lane);
  if( tidx == 0 ){
    a2 *= RSQRT(a2);
    if( singleCall )      amin_dev[BLOCKIDX_X1D] = a2;
    else      atomicMin(&amin_dev[BLOCKIDX_X1D], a2);
  }/* if( tidx == 0 ){ */
}
#endif//GADGET_MAC


/**
 * @fn calc_r2max_dev
 *
 * @brief Get r2max in the array.
 */
extern "C"
void calc_r2max_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, iparticle *pi, soaGEO dev
#ifdef  EXEC_BENCHMARK
		    , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		    )
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  (*pi).encBall_hst->m = ZERO;
#ifdef  GADGET_MAC
  pi->amin = REAL_MAX;
#endif//GADGET_MAC

  /** thread-block structure must be identical to tree traversal */
#ifndef SERIALIZED_EXECUTION
  if( Ngrp != 0 )
#endif//SERIALIZED_EXECUTION
    {
      int Nrem = BLOCKSIZE(Ngrp, NWARP * NGROUPS);

      /** NOTE: pi.jpos contains position of i-particle on the time of the gravity calculation */
      if( Nrem <= dev.Nblock ){
	calc_r2max_kernel<<<Nrem, NTHREADS>>>(BLOCKSIZE(Ngrp, NGROUPS) * NGROUPS, laneInfo, (*pi).encBall,
#ifdef  BLOCK_TIME_STEP
					      (*pi).jpos,
#else///BLOCK_TIME_STEP
					      (*pi).pos,
#endif//BLOCK_TIME_STEP
					      dev.r2, true);

	/** reduction using library */
	(*pi).encBall_hst->m = thrust::reduce((thrust::device_ptr<real>)dev.r2, (thrust::device_ptr<real>)(dev.r2 + Nrem), ZERO, thrust::maximum<real>());

#ifdef  GADGET_MAC
	calc_amin_kernel<<<Nrem, NTHREADS>>>(BLOCKSIZE(Ngrp, NGROUPS) * NGROUPS, laneInfo, (*pi).acc, dev.r2, true);
	pi->amin = thrust::reduce((thrust::device_ptr<real>)dev.r2, (thrust::device_ptr<real>)(dev.r2 + Nrem), REAL_MAX, thrust::minimum<real>());
#endif//GADGET_MAC
      }/* if( Nrem <= dev.Nblock ){ */
      else{
	const int Niter = BLOCKSIZE(Nrem, dev.Nblock);
	int hidx = 0;
	init_r2max_kernel<<<BLOCKSIZE(dev.Nblock, 1024), 1024>>>(dev.Nblock, dev.r2);

	for(int iter = 0; iter < Niter; iter++){
	  int Nblck = dev.Nblock;
	  if( Nblck > Nrem )	    Nblck = Nrem;

	  int Nsub = Nblck * NWARP * NGROUPS;
	  calc_r2max_kernel<<<Nblck, NTHREADS>>>(BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx], (*pi).encBall,
#ifdef  BLOCK_TIME_STEP
						 (*pi).jpos,
#else///BLOCK_TIME_STEP
						 (*pi).pos,
#endif//BLOCK_TIME_STEP
						 dev.r2, false);

	  hidx += Nsub;
	  Nrem -= Nblck;
	}/* for(int iter = 0; iter < Niter; iter++){ */

	/** reduction using library */
	(*pi).encBall_hst->m = thrust::reduce((thrust::device_ptr<real>)dev.r2, (thrust::device_ptr<real>)(dev.r2 + dev.Nblock), ZERO, thrust::maximum<real>());

#ifdef  GADGET_MAC
	/** initialization */
	Nrem = BLOCKSIZE(Ngrp, NWARP * NGROUPS);
	hidx = 0;
	init_amin_kernel<<<BLOCKSIZE(dev.Nblock, 1024), 1024>>>(dev.Nblock, dev.r2);

	/** reduction within a block */
	for(int iter = 0; iter < Niter; iter++){
	  int Nblck = dev.Nblock;
	  if( Nblck > Nrem )	    Nblck = Nrem;

	  int Nsub = Nblck * NWARP * NGROUPS;
	  calc_amin_kernel<<<Nblck, NTHREADS>>>(BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx], (*pi).acc, dev.r2, false);

	  hidx += Nsub;
	  Nrem -= Nblck;
	}/* for(int iter = 0; iter < Niter; iter++){ */

	/** reduction using library */
	pi->amin = thrust::reduce((thrust::device_ptr<real>)dev.r2, (thrust::device_ptr<real>)(dev.r2 + dev.Nblock), REAL_MAX, thrust::minimum<real>());
#endif//GADGET_MAC
      }/* else{ */

      getLastCudaError("calc_r2max_kernel");

    }

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->calc_r2max_dev));
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}


/**
 * @fn allocGeometricEnclosingBall_dev
 *
 * @brief Memory allocation for geometric enclosing ball.
 */
extern "C"
muse allocGeometricEnclosingBall_dev(real **r2_dev, soaGEO *dev, const int num_max)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};

  int Nblock = BLOCKSIZE(NWARP * num_max, NTHREADS);
  if( Nblock > MAX_BLOCKS_PER_GRID )
    Nblock = MAX_BLOCKS_PER_GRID;

  mycudaMalloc((void **)r2_dev, Nblock * sizeof(real));  alloc.device += Nblock * sizeof(real);

  dev->r2 = *r2_dev;
  dev->Nblock = Nblock;

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freeGeometricEnclosingBall_dev
 *
 * @brief Memory deallocation for geometric enclosing ball.
 */
extern "C"
void  freeGeometricEnclosingBall_dev(real  *r2_dev)
{
  __NOTE__("%s\n", "start");

  mycudaFree(r2_dev);

  __NOTE__("%s\n", "end");
}
