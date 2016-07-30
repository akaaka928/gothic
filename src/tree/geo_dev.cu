/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/29(Fri) 11:18:40
 *                                                                       *
 *    Generation of enclosing ball containing all N-body particles       *
 *   the center is the geometric one of the enclosing rectangular cuboid *
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
/* library load to get the maximum value */
/* #ifdef  CUB_AVAILABLE */
/* #include <cub/device/device_reduce.cuh> */
/* #else///CUB_AVAILABLE */
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
/* #endif//CUB_AVAILABLE */
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "../tree/make.h"
#include "../tree/walk_dev.h"
#include "geo_dev.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__global__ void init_r2max_kernel(const int num, real *r2)
{
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------
  if( tidx < num )
    r2[tidx] = ZERO;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__global__ void init_amin_kernel(const int num, real *amin)
{
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------
  if( tidx < num )
    amin[tidx] = REAL_MAX;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* minimum value within a block */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
__device__ __forceinline__ real getMin_block(const real val, const int tidx, const int lane, volatile real * RESTRICT smem)
{
  //-----------------------------------------------------------------------
  /* 1. reduction within a warp */
  //-----------------------------------------------------------------------
  real min = val;
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  real tmp;
  tmp = __shfl_up(min,  1, warpSize);  if( tmp < min )    min = tmp;
  tmp = __shfl_up(min,  2, warpSize);  if( tmp < min )    min = tmp;
  tmp = __shfl_up(min,  4, warpSize);  if( tmp < min )    min = tmp;
  tmp = __shfl_up(min,  8, warpSize);  if( tmp < min )    min = tmp;
  tmp = __shfl_up(min, 16, warpSize);  if( tmp < min )    min = tmp;
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx] = min;
  if( smem[tidx -  1] < min ){    min = smem[tidx -  1];    smem[tidx] = min;  }
  if( smem[tidx -  2] < min ){    min = smem[tidx -  2];    smem[tidx] = min;  }
  if( smem[tidx -  4] < min ){    min = smem[tidx -  4];    smem[tidx] = min;  }
  if( smem[tidx -  8] < min ){    min = smem[tidx -  8];    smem[tidx] = min;  }
  if( smem[tidx - 16] < min ){    min = smem[tidx - 16];    smem[tidx] = min;  }
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  /* return calculated inclusive prefix sum */
  smem[tidx] = min;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#   if  NTHREADS >= 64
  //-----------------------------------------------------------------------
  /* 2. reduction among multiple warps */
  //-----------------------------------------------------------------------
  __syncthreads();
  /* warpSize = 32 = 2^5 */
  if( tidx < (NTHREADS >> 5) ){
    //---------------------------------------------------------------------
    min = smem[tidx * warpSize];
#ifdef  USE_WARP_SHUFFLE_FUNC
#   if  NTHREADS >=   64
    const int groupSize = NTHREADS >> 5;
    tmp = __shfl_up(min,  1, groupSize);    if( tmp < min )      min = tmp;
#   if  NTHREADS >=  128
    tmp = __shfl_up(min,  2, groupSize);    if( tmp < min )      min = tmp;
#   if  NTHREADS >=  256
    tmp = __shfl_up(min,  4, groupSize);    if( tmp < min )      min = tmp;
#   if  NTHREADS >=  512
    tmp = __shfl_up(min,  8, groupSize);    if( tmp < min )      min = tmp;
#   if  NTHREADS == 1024
    tmp = __shfl_up(min, 16, groupSize);    if( tmp < min )      min = tmp;
#endif//NTHREADS == 1024
#endif//NTHREADS >=  512
#endif//NTHREADS >=  256
#endif//NTHREADS >=  128
#endif//NTHREADS >=   64
#else///USE_WARP_SHUFFLE_FUNC
    smem[tidx] = min;
#   if  NTHREADS >=   64
    if( smem[tidx -  1] < min ){      min = smem[tidx -  1];      smem[tidx] = min;    }
#   if  NTHREADS >=  128
    if( smem[tidx -  2] < min ){      min = smem[tidx -  2];      smem[tidx] = min;    }
#   if  NTHREADS >=  256
    if( smem[tidx -  4] < min ){      min = smem[tidx -  4];      smem[tidx] = min;    }
#   if  NTHREADS >=  512
    if( smem[tidx -  8] < min ){      min = smem[tidx -  8];      smem[tidx] = min;    }
#   if  NTHREADS == 1024
    if( smem[tidx - 16] < min ){      min = smem[tidx - 16];      smem[tidx] = min;    }
#endif//NTHREADS == 1024
#endif//NTHREADS >=  512
#endif//NTHREADS >=  256
#endif//NTHREADS >=  128
#endif//NTHREADS >=   64
    min = smem[tidx];
#endif//USE_WARP_SHUFFLE_FUNC
    //---------------------------------------------------------------------
    smem[tidx] = min;
    //---------------------------------------------------------------------
  }/* if( tidx < (NTHREADS >> 5) ){ */
  __syncthreads();
  //-----------------------------------------------------------------------
#endif//NTHREADS >= 64
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  min = smem[0];
  return (min);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* maximum value within a block */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
__device__ __forceinline__ real getMax_block(const real val, const int tidx, const int lane, volatile real * RESTRICT smem)
{
  //-----------------------------------------------------------------------
  /* 1. reduction within a warp */
  //-----------------------------------------------------------------------
  real max = val;
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  real tmp;
  tmp = __shfl_up(max,  1, warpSize);  if( tmp > max )    max = tmp;
  tmp = __shfl_up(max,  2, warpSize);  if( tmp > max )    max = tmp;
  tmp = __shfl_up(max,  4, warpSize);  if( tmp > max )    max = tmp;
  tmp = __shfl_up(max,  8, warpSize);  if( tmp > max )    max = tmp;
  tmp = __shfl_up(max, 16, warpSize);  if( tmp > max )    max = tmp;
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx] = max;
  if( smem[tidx -  1] > max ){    max = smem[tidx -  1];    smem[tidx] = max;  }
  if( smem[tidx -  2] > max ){    max = smem[tidx -  2];    smem[tidx] = max;  }
  if( smem[tidx -  4] > max ){    max = smem[tidx -  4];    smem[tidx] = max;  }
  if( smem[tidx -  8] > max ){    max = smem[tidx -  8];    smem[tidx] = max;  }
  if( smem[tidx - 16] > max ){    max = smem[tidx - 16];    smem[tidx] = max;  }
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  /* return calculated inclusive prefix sum */
  smem[tidx] = max;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#   if  NTHREADS >= 64
  //-----------------------------------------------------------------------
  /* 2. reduction among multiple warps */
  //-----------------------------------------------------------------------
  __syncthreads();
  /* warpSize = 32 = 2^5 */
  if( tidx < (NTHREADS >> 5) ){
    //---------------------------------------------------------------------
    max = smem[tidx * warpSize];
#ifdef  USE_WARP_SHUFFLE_FUNC
#   if  NTHREADS >=   64
    const int groupSize = NTHREADS >> 5;
    tmp = __shfl_up(max,  1, groupSize);    if( tmp > max )      max = tmp;
#   if  NTHREADS >=  128
    tmp = __shfl_up(max,  2, groupSize);    if( tmp > max )      max = tmp;
#   if  NTHREADS >=  256
    tmp = __shfl_up(max,  4, groupSize);    if( tmp > max )      max = tmp;
#   if  NTHREADS >=  512
    tmp = __shfl_up(max,  8, groupSize);    if( tmp > max )      max = tmp;
#   if  NTHREADS == 1024
    tmp = __shfl_up(max, 16, groupSize);    if( tmp > max )      max = tmp;
#endif//NTHREADS == 1024
#endif//NTHREADS >=  512
#endif//NTHREADS >=  256
#endif//NTHREADS >=  128
#endif//NTHREADS >=   64
#else///USE_WARP_SHUFFLE_FUNC
    smem[tidx] = max;
#   if  NTHREADS >=   64
    if( smem[tidx -  1] > max ){      max = smem[tidx -  1];      smem[tidx] = max;    }
#   if  NTHREADS >=  128
    if( smem[tidx -  2] > max ){      max = smem[tidx -  2];      smem[tidx] = max;    }
#   if  NTHREADS >=  256
    if( smem[tidx -  4] > max ){      max = smem[tidx -  4];      smem[tidx] = max;    }
#   if  NTHREADS >=  512
    if( smem[tidx -  8] > max ){      max = smem[tidx -  8];      smem[tidx] = max;    }
#   if  NTHREADS == 1024
    if( smem[tidx - 16] > max ){      max = smem[tidx - 16];      smem[tidx] = max;    }
#endif//NTHREADS == 1024
#endif//NTHREADS >=  512
#endif//NTHREADS >=  256
#endif//NTHREADS >=  128
#endif//NTHREADS >=   64
    max = smem[tidx];
#endif//USE_WARP_SHUFFLE_FUNC
    //---------------------------------------------------------------------
    smem[tidx] = max;
    //---------------------------------------------------------------------
  }/* if( tidx < (NTHREADS >> 5) ){ */
  __syncthreads();
  //-----------------------------------------------------------------------
#endif//NTHREADS >= 64
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  max = smem[0];
  return (max);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__device__ __forceinline__ float atomicMin(float* addr, float val)
{
  //-----------------------------------------------------------------------
  int* addr_as_i = (int*)addr;
  int old = *addr_as_i, assumed;
  //-----------------------------------------------------------------------
  do{
    assumed = old;
    old = atomicCAS(addr_as_i, assumed, __float_as_int(fminf(val, __int_as_float(assumed))));
  } while( assumed != old );
  //-----------------------------------------------------------------------
  return (__int_as_float(old));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ float atomicMax(float* addr, float val)
{
  //-----------------------------------------------------------------------
  int* addr_as_i = (int*)addr;
  int old = *addr_as_i, assumed;
  //-----------------------------------------------------------------------
  do{
    assumed = old;
    old = atomicCAS(addr_as_i, assumed, __float_as_int(fmaxf(val, __int_as_float(assumed))));
  } while( assumed != old );
  //-----------------------------------------------------------------------
  return (__int_as_float(old));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__global__ void calc_r2max_kernel
(const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo, READ_ONLY position * RESTRICT cen_dev, position * RESTRICT jpos, real * RESTRICT r2_dev, const bool singleCall)
{
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int lane    = tidx          & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);
  //-----------------------------------------------------------------------
  /* const int head = tidx - lane; */
  __shared__ real smem[NTHREADS];
  //-----------------------------------------------------------------------
  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];
  //-----------------------------------------------------------------------
  real r2 = ZERO;
  position cen = *cen_dev;
  //-----------------------------------------------------------------------
  if( lane < info.num ){
    //---------------------------------------------------------------------
    const int idx = info.head + lane;
    const position pi = jpos[idx];
    //---------------------------------------------------------------------
    const real dx = pi.x - cen.x;
    const real dy = pi.y - cen.y;
    const real dz = pi.z - cen.z;
    r2 = dx * dx + dy * dy + dz * dz;
    //---------------------------------------------------------------------
  }/* if( lane < info.num ){ */
  //-----------------------------------------------------------------------
  r2 = getMax_block(r2, tidx, tidx & (warpSize - 1), smem);
  if( tidx == 0 ){
    if( singleCall )      r2_dev[BLOCKIDX_X1D] = r2;
    else      atomicMax(&r2_dev[BLOCKIDX_X1D], r2);
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifdef  GADGET_MAC
__global__ void calc_amin_kernel
(const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo, acceleration * RESTRICT iacc, real * RESTRICT amin_dev, const bool singleCall)
{
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int lane    = tidx          & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);
  //-----------------------------------------------------------------------
  /* const int head = tidx - lane; */
  __shared__ real smem[NTHREADS];
  //-----------------------------------------------------------------------
  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];
  //-----------------------------------------------------------------------
  real a2 = REAL_MAX;
  //-----------------------------------------------------------------------
  if( lane < info.num ){
    //---------------------------------------------------------------------
    const int idx = info.head + lane;
    const acceleration ai = iacc[idx];
    //---------------------------------------------------------------------
    a2 = 1.0e-30f + ai.x * ai.x + ai.y * ai.y + ai.z * ai.z;
    //---------------------------------------------------------------------
  }/* if( lane < info.num ){ */
  //-----------------------------------------------------------------------
  a2 = getMin_block(a2, tidx, tidx & (warpSize - 1), smem);
  if( tidx == 0 ){
    a2 *= RSQRT(a2);
    if( singleCall )      amin_dev[BLOCKIDX_X1D] = a2;
    else      atomicMin(&amin_dev[BLOCKIDX_X1D], a2);
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
}
#endif//GADGET_MAC
//-------------------------------------------------------------------------
extern "C"
void calc_r2max_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, iparticle *pi, soaGEO dev
#ifdef  EXEC_BENCHMARK
		    , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		    )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  (*pi).encBall_hst->m = ZERO;
#ifdef  GADGET_MAC
  pi->amin = REAL_MAX;
#endif//GADGET_MAC
  //-----------------------------------------------------------------------
  /* thread-block structure must be identical to tree traversal */
#ifndef SERIALIZED_EXECUTION
  if( Ngrp != 0 )
#endif//SERIALIZED_EXECUTION
    {
      //-------------------------------------------------------------------
      int Nrem = BLOCKSIZE(Ngrp, NWARP * NGROUPS);
      //-------------------------------------------------------------------
      /* NOTE: pi.jpos contains position of i-particle on the time of the gravity calculation */
      if( Nrem <= dev.Nblock ){
	//-----------------------------------------------------------------
	calc_r2max_kernel<<<Nrem, NTHREADS>>>(BLOCKSIZE(Ngrp, NGROUPS) * NGROUPS, laneInfo, (*pi).encBall, (*pi).jpos, dev.r2, true);
	//-----------------------------------------------------------------
	/* reduction using library */
/* #ifdef  CUB_AVAILABLE */
/* 	cub::DeviceReduce::Max(dev.temp_storage, dev.temp_storage_size, dev.r2, (*pi).encBall_hst->m, Nrem); */
/* #else///CUB_AVAILABLE */
	(*pi).encBall_hst->m = thrust::reduce((thrust::device_ptr<real>)dev.r2, (thrust::device_ptr<real>)(dev.r2 + Nrem), ZERO, thrust::maximum<real>());
/* #endif//CUB_AVAILABLE */
	//-----------------------------------------------------------------
#ifdef  GADGET_MAC
	calc_amin_kernel<<<Nrem, NTHREADS>>>(BLOCKSIZE(Ngrp, NGROUPS) * NGROUPS, laneInfo, (*pi).acc, dev.r2, true);
	pi->amin = thrust::reduce((thrust::device_ptr<real>)dev.r2, (thrust::device_ptr<real>)(dev.r2 + Nrem), REAL_MAX, thrust::minimum<real>());
#endif//GADGET_MAC
	//-----------------------------------------------------------------
      }/* if( Nrem <= dev.Nblock ){ */
      //-------------------------------------------------------------------
      else{
	//-----------------------------------------------------------------
	const int Niter = BLOCKSIZE(Nrem, dev.Nblock);
	int hidx = 0;
	init_r2max_kernel<<<BLOCKSIZE(dev.Nblock, 1024), 1024>>>(dev.Nblock, dev.r2);
	//-----------------------------------------------------------------
	for(int iter = 0; iter < Niter; iter++){
	  //---------------------------------------------------------------
	  int Nblck = dev.Nblock;
	  if( Nblck > Nrem )	    Nblck = Nrem;
	  //---------------------------------------------------------------
	  int Nsub = Nblck * NWARP * NGROUPS;
	  calc_r2max_kernel<<<Nblck, NTHREADS>>>(BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx], (*pi).encBall, (*pi).jpos, dev.r2, false);
	  //-------------------------------------------------------------------
	  hidx += Nsub;
	  Nrem -= Nblck;
	//-----------------------------------------------------------------
	}/* for(int iter = 0; iter < Niter; iter++){ */
	//-----------------------------------------------------------------
	/* reduction using library */
/* #ifdef  CUB_AVAILABLE */
/* 	cub::DeviceReduce::Max(dev.temp_storage, dev.temp_storage_size, dev.r2, (*pi).encBall_hst->m, dev.Nblock); */
/* #else///CUB_AVAILABLE */
	(*pi).encBall_hst->m = thrust::reduce((thrust::device_ptr<real>)dev.r2, (thrust::device_ptr<real>)(dev.r2 + dev.Nblock), ZERO, thrust::maximum<real>());
/* #endif//CUB_AVAILABLE */
	//-----------------------------------------------------------------
#ifdef  GADGET_MAC
	//-----------------------------------------------------------------
	/* initialization */
	Nrem = BLOCKSIZE(Ngrp, NWARP * NGROUPS);
	hidx = 0;
	init_amin_kernel<<<BLOCKSIZE(dev.Nblock, 1024), 1024>>>(dev.Nblock, dev.r2);
	/* reduction within a block */
	for(int iter = 0; iter < Niter; iter++){
	  //---------------------------------------------------------------
	  int Nblck = dev.Nblock;
	  if( Nblck > Nrem )	    Nblck = Nrem;
	  //---------------------------------------------------------------
	  int Nsub = Nblck * NWARP * NGROUPS;
	  calc_amin_kernel<<<Nblck, NTHREADS>>>(BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx], (*pi).acc, dev.r2, false);
	  //-------------------------------------------------------------------
	  hidx += Nsub;
	  Nrem -= Nblck;
	//-----------------------------------------------------------------
	}/* for(int iter = 0; iter < Niter; iter++){ */
	/* reduction using library */
	pi->amin = thrust::reduce((thrust::device_ptr<real>)dev.r2, (thrust::device_ptr<real>)(dev.r2 + dev.Nblock), REAL_MAX, thrust::minimum<real>());
	//-----------------------------------------------------------------
#endif//GADGET_MAC
	//-----------------------------------------------------------------
      }/* else{ */
      //-------------------------------------------------------------------
      getLastCudaError("calc_r2max_kernel");
      //-------------------------------------------------------------------
    }
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->calc_r2max_dev));
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
extern "C"
muse allocGeometricEnclosingBall_dev
(real **r2_dev
/* #ifdef  CUB_AVAILABLE */
/*  , void **temp_storage */
/* #endif//CUB_AVAILABLE */
 , soaGEO *dev, const int num_max)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  int Nblock = DIV_NTHREADS(NWARP * num_max);
  if( Nblock > MAX_BLOCKS_PER_GRID )
    Nblock = MAX_BLOCKS_PER_GRID;
  //-----------------------------------------------------------------------
  /* memory allocation and simple confirmation */
  mycudaMalloc((void **)r2_dev, Nblock * sizeof(real));  alloc.device += Nblock * sizeof(real);
  //-----------------------------------------------------------------------
  dev->r2 = *r2_dev;
  dev->Nblock = Nblock;
  //-----------------------------------------------------------------------
/* #ifdef  CUB_AVAILABLE */
/*   real r2max_tmp; */
/*   size_t temp_storage_size = 0; */
/*   *temp_storage = NULL; */
/*   cub::DeviceReduce::Max(*temp_storage, temp_storage_size, dev.r2, &r2max_tmp, Nblock); */
/*   mycudaMalloc(temp_storage, temp_storage_size);  alloc.device += temp_storage_size; */
/*   dev->temp_storage = *temp_storage;  dev->temp_storage_size = temp_storage_size; */
/* #endif//CUB_AVAILABLE */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
void  freeGeometricEnclosingBall_dev
(real  *r2_dev
/* #ifdef  CUB_AVAILABLE */
/*  , void  *temp_storage */
/* #endif//CUB_AVAILABLE */
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  mycudaFree(r2_dev);
/* #ifdef  CUB_AVAILABLE */
/*   mycudaFree(temp_storage); */
/* #endif//CUB_AVAILABLE */
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
