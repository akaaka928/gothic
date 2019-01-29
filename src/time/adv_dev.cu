/**
 * @file adv_dev.cu
 *
 * @brief Source code for orbit integration of N-body particles
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/12/28 (Fri)
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
#include <helper_cuda.h>

#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/find.h>
#include <thrust/sort.h>

#include "macro.h"
#include "cudalib.h"
#include "sys/time.h"
#include "timer.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/device.h"

#ifndef SERIALIZED_EXECUTION
#include <mpi.h>
#include "mpilib.h"
#include "../para/mpicfg.h"
#endif//SERIALIZED_EXECUTION

#include "../tree/walk_dev.h"

#include "adv_dev.h"


#ifndef BLOCK_TIME_STEP
/**
 * @fn allocTimeStep_dev
 *
 * @brief Memory allocation for time step.
 */
extern "C"
muse allocTimeStep_dev(real **dt_dev)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};
  mycudaMalloc((void **)dt_dev, 1 * sizeof(real));
  alloc.device +=               1 * sizeof(real);

  __NOTE__("%s\n", "end");
  return (alloc);
}

/**
 * @fn freeTimeStep_dev
 *
 * @brief Memory deallocation for time step.
 */
extern "C"
void  freeTimeStep_dev(real  *dt_dev)
{
  __NOTE__("%s\n", "start");

  mycudaFree(dt_dev);

  __NOTE__("%s\n", "end");
}
#endif//BLOCK_TIME_STEP


#ifdef  BLOCK_TIME_STEP
/**
 * @fn adjustTimeStep_kernel
 *
 * @brief Adjust time step for N-body particles.
 */
__global__ void adjustTimeStep_kernel
(const double tnew, const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo, velocity * RESTRICT ivel, ibody_time * RESTRICT time)
{
  const int lane    = THREADIDX_X1D & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);

  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//__CUDA_ARCH__ >= 700

  if( lane < info.num ){
    const int idx = info.head + lane;

    ibody_time ti = time[idx];
    ti.t1 = tnew;
    time[idx] = ti;

    velocity vi = ivel[idx];
    vi.dt = (real)(ti.t1 - ti.t0);
    ivel[idx] = vi;
  }/* if( lane < info.num ){ */
}


struct greater_than
{
  double val;
greater_than(double val) : val(val) {}
  __host__ __device__
  int operator()(const double &x) const {
    return (x > val);
  }
};


/**
 * @fn setTimeStep_dev
 *
 * @brief Set time step for N-body particles.
 *
 * @sa adjustTimeStep_kernel
 */
extern "C"
void setTimeStep_dev
(const int Ngrp, laneinfo * RESTRICT laneInfo_dev, double * RESTRICT laneTime_dev, int *grpNum, const iparticle pi,
 const double told, double *tnew, double *dt, bool adjustAllTimeStep, const double invSnapshotInterval, const uint previous, uint *present
#ifndef SERIALIZED_EXECUTION
 , MPIcfg_tree mpi
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
)
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  thrust::stable_sort_by_key((thrust::device_ptr<double>)laneTime_dev, (thrust::device_ptr<double>)(laneTime_dev + Ngrp), (thrust::device_ptr<laneinfo>)laneInfo_dev);
  checkCudaErrors(cudaMemcpy(tnew, laneTime_dev, sizeof(double), cudaMemcpyDeviceToHost));
#ifndef SERIALIZED_EXECUTION
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, tnew, 1, MPI_DOUBLE, MPI_MIN, mpi.comm));
#endif//SERIALIZED_EXECUTION

  *dt = (*tnew) - told;
  *present = (uint)((*tnew) * invSnapshotInterval);
  if( *present != previous )
    adjustAllTimeStep = true;

  if( !adjustAllTimeStep ){
    thrust::device_vector<double>::iterator iter1 =                 (thrust::device_ptr<double>)laneTime_dev;
    thrust::device_vector<double>::iterator iter2 = thrust::find_if((thrust::device_ptr<double>)laneTime_dev, (thrust::device_ptr<double>)(laneTime_dev + Ngrp), greater_than(*tnew));
    *grpNum = thrust::distance(iter1, iter2);
  }/* if( !adjustAllTimeStep ){ */
  else
    *grpNum = Ngrp;

#ifndef SERIALIZED_EXECUTION
  if( *grpNum > 0 )
#endif//SERIALIZED_EXECUTION
    {
      int Nrem = BLOCKSIZE(*grpNum, NWARP * NGROUPS);
      if( Nrem <= MAX_BLOCKS_PER_GRID )
	adjustTimeStep_kernel<<<Nrem, NTHREADS>>>(*tnew, BLOCKSIZE(*grpNum, NGROUPS) * NGROUPS, laneInfo_dev, pi.vel, pi.time);
      else{
	const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
	int hidx = 0;

	for(int iter = 0; iter < Niter; iter++){
	  int Nblck = MAX_BLOCKS_PER_GRID;
	  if( Nblck > Nrem )	    Nblck = Nrem;

	  int Nsub = Nblck * NWARP * NGROUPS;
	  adjustTimeStep_kernel<<<Nblck, NTHREADS>>>(*tnew, BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo_dev[hidx], pi.vel, pi.time);

	  hidx += Nsub;
	  Nrem -= Nblck;
	}/* for(int iter = 0; iter < Niter; iter++){ */
      }/* else{ */

      getLastCudaError("adjustTimeStep_kernel");
    }

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->setTimeStep_dev));
#endif//EXEC_BENCHMARK


#ifdef  SHOW_NI_DEPENDENCE
  laneinfo *laneInfo_hst;  mycudaMallocHost((void **)&laneInfo_hst, sizeof(laneinfo) * Ngrp);
  double   *laneTime_hst;  mycudaMallocHost((void **)&laneTime_hst, sizeof(double  ) * Ngrp);
  checkCudaErrors(cudaMemcpy(laneInfo_hst, laneInfo_dev, sizeof(laneinfo) * Ngrp, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(laneTime_hst, laneTime_dev, sizeof(double  ) * Ngrp, cudaMemcpyDeviceToHost));
  printf("#group ID\tnew time\n");
  int Ni_active = 0;
  for(int ii = 0; ii < Ngrp; ii++){
    Ni_active += laneInfo_hst[ii].num;
    printf("%d\t%d\t%e\n", ii, Ni_active, laneTime_hst[ii]);
  }/* for(int ii = 0; ii < Ngrp; ii++){ */
  printf("#*grpNum is %d\n\n", *grpNum);
  fflush(stdout);
  mycudaFreeHost(laneInfo_hst);
  mycudaFreeHost(laneTime_hst);
#endif//SHOW_NI_DEPENDENCE


  __NOTE__("%s\n", "end");
}
#else///BLOCK_TIME_STEP
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
#   if  NTHREADS_TIME ==   64
#define SHFL_MASK_WARPS SHFL_MASK_02
#endif//NTHREADS_TIME ==   64
#   if  NTHREADS_TIME ==  128
#define SHFL_MASK_WARPS SHFL_MASK_04
#endif//NTHREADS_TIME ==  128
#   if  NTHREADS_TIME ==  256
#define SHFL_MASK_WARPS SHFL_MASK_08
#endif//NTHREADS_TIME ==  256
#   if  NTHREADS_TIME ==  512
#define SHFL_MASK_WARPS SHFL_MASK_16
#endif//NTHREADS_TIME ==  512
#   if  NTHREADS_TIME == 1024
#define SHFL_MASK_WARPS SHFL_MASK_32
#endif//NTHREADS_TIME == 1024
#endif//USE_WARP_SHUFFLE_FUNC_TIME
/**
 * @fn setTimeStep_kernel
 *
 * @brief Set time step for N-body particles.
 */
__global__ void setTimeStep_kernel
(const int Ni,
#ifndef OMIT_VELOCITY_FOR_TIME_STEP
 READ_ONLY real * RESTRICT vix, READ_ONLY real * RESTRICT viy, READ_ONLY real * RESTRICT viz,
#endif//OMIT_VELOCITY_FOR_TIME_STEP
 READ_ONLY acceleration * RESTRICT iacc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 READ_ONLY acceleration * RESTRICT iacc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
 const real eta, const real eps,
 real * RESTRICT dt)
{
  const int tidx = THREADIDX_X1D;

  /** calculate time step of individual N-body particle */
  real dtloc = REAL_MAX;

  for(int ih = 0; ih < Ni; ih += NTHREADS_TIME){
    const int ii = ih + tidx;

    if( ii < Ni ){
#ifndef OMIT_VELOCITY_FOR_TIME_STEP
      const real         vx = vix [ii];
      const real         vy = viy [ii];
      const real         vz = viz [ii];
#endif//OMIT_VELOCITY_FOR_TIME_STEP
#ifndef SET_EXTERNAL_POTENTIAL_FIELD
      const acceleration ai = iacc[ii];
#else///SET_EXTERNAL_POTENTIAL_FIELD
      const acceleration slf = iacc    [ii];
      const acceleration ext = iacc_ext[ii];
      const acceleration ai = {slf.x + ext.x, slf.y + ext.y, slf.z + ext.z, slf.pot + ext.pot};
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#ifndef OMIT_VELOCITY_FOR_TIME_STEP
      const real v2 = FLT_MIN +	  vx *	 vx +	vy *   vy +   vz *   vz;
#endif//OMIT_VELOCITY_FOR_TIME_STEP
      const real a2 = FLT_MIN + ai.x * ai.x + ai.y * ai.y + ai.z * ai.z;

#ifndef OMIT_VELOCITY_FOR_TIME_STEP
      const real vdt =      eps * RSQRT(v2);
#endif//OMIT_VELOCITY_FOR_TIME_STEP
      const real adt = SQRT(eps * RSQRT(a2));

#ifndef OMIT_VELOCITY_FOR_TIME_STEP
      dtloc = FMIN(dtloc, FMIN(vdt, adt));
#else///OMIT_VELOCITY_FOR_TIME_STEP
      dtloc = FMIN(dtloc, adt);
#endif//OMIT_VELOCITY_FOR_TIME_STEP
    }/* if( ii < Ni ){ */
  }/* for(int ih = 0; ih < Ni; ih += NTHREADS_TIME){ */
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//__CUDA_ARCH__ >= 700


  /** find the minimum time step */
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
  __shared__ real dtmin[32];
#else///USE_WARP_SHUFFLE_FUNC_TIME
  __shared__ real dtmin[NTHREADS_TIME];
#endif//USE_WARP_SHUFFLE_FUNC_TIME

  /** find minimum time step within a warp */
  real dttmp;
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
  dttmp = __SHFL_XOR(SHFL_MASK_32, dtloc,  1, warpSize);  dtloc = FMIN(dtloc, dttmp);
  dttmp = __SHFL_XOR(SHFL_MASK_32, dtloc,  2, warpSize);  dtloc = FMIN(dtloc, dttmp);
  dttmp = __SHFL_XOR(SHFL_MASK_32, dtloc,  4, warpSize);  dtloc = FMIN(dtloc, dttmp);
  dttmp = __SHFL_XOR(SHFL_MASK_32, dtloc,  8, warpSize);  dtloc = FMIN(dtloc, dttmp);
  dttmp = __SHFL_XOR(SHFL_MASK_32, dtloc, 16, warpSize);  dtloc = FMIN(dtloc, dttmp);
  if( (tidx & (warpSize - 1)) == 0 )
    dtmin[tidx / warpSize] = dtloc;
#else///USE_WARP_SHUFFLE_FUNC_TIME
  dtmin[tidx] = dtloc;
  dttmp = dtmin[tidx ^  1];  dtloc = FMIN(dtloc, dttmp);  dtmin[tidx] = dtloc;/**< w/ \pm  1 */
  dttmp = dtmin[tidx ^  2];  dtloc = FMIN(dtloc, dttmp);  dtmin[tidx] = dtloc;/**< w/ \pm  2 */
  dttmp = dtmin[tidx ^  4];  dtloc = FMIN(dtloc, dttmp);  dtmin[tidx] = dtloc;/**< w/ \pm  4 */
  dttmp = dtmin[tidx ^  8];  dtloc = FMIN(dtloc, dttmp);  dtmin[tidx] = dtloc;/**< w/ \pm  8 */
  dttmp = dtmin[tidx ^ 16];  dtloc = FMIN(dtloc, dttmp);  dtmin[tidx] = dtloc;/**< w/ \pm 16 */
#endif//USE_WARP_SHUFFLE_FUNC_TIME

  /** warpSize^2 = 32^2 = 1024 is the maximum of the number of threads */
  __syncthreads();

  /** warpSize = 32 = 2^5 -->> NTHREADS_TIME / warpSize = NTHREADS_TIME >> 5 */
  if( tidx < (NTHREADS_TIME >> 5) ){
    /** share the minimum time step in each warp */
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
    dtloc = dtmin[tidx];
#else///USE_WARP_SHUFFLE_FUNC_TIME
    dttmp = dtmin[tidx * warpSize];
#endif//USE_WARP_SHUFFLE_FUNC_TIME

    /** find the minimum time step within the whole threads */
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
#   if  NTHREADS_TIME >=   64
    dttmp = __SHFL_XOR(SHFL_MASK_WARPS, dtloc,  1, NTHREADS_TIME >> 5);  dtloc = FMIN(dtloc, dttmp);
#   if  NTHREADS_TIME >=  128
    dttmp = __SHFL_XOR(SHFL_MASK_WARPS, dtloc,  2, NTHREADS_TIME >> 5);  dtloc = FMIN(dtloc, dttmp);
#   if  NTHREADS_TIME >=  256
    dttmp = __SHFL_XOR(SHFL_MASK_WARPS, dtloc,  4, NTHREADS_TIME >> 5);  dtloc = FMIN(dtloc, dttmp);
#   if  NTHREADS_TIME >=  512
    dttmp = __SHFL_XOR(SHFL_MASK_WARPS, dtloc,  8, NTHREADS_TIME >> 5);  dtloc = FMIN(dtloc, dttmp);
#   if  NTHREADS_TIME == 1024
    dttmp = __SHFL_XOR(SHFL_MASK_WARPS, dtloc, 16, NTHREADS_TIME >> 5);  dtloc = FMIN(dtloc, dttmp);
#endif//NTHREADS_TIME == 1024
#endif//NTHREADS_TIME >=  512
#endif//NTHREADS_TIME >=  256
#endif//NTHREADS_TIME >=  128
#endif//NTHREADS_TIME >=   64
#else///USE_WARP_SHUFFLE_FUNC_TIME
    dtmin[tidx] = dttmp;
#   if  NTHREADS_TIME >=   64
    dttmp = dtmin[tidx ^  1];    dtloc = FMIN(dtloc, dttmp);	dtmin[tidx] = dtloc;
#   if  NTHREADS_TIME >=  128
    dttmp = dtmin[tidx ^  2];	 dtloc = FMIN(dtloc, dttmp);	dtmin[tidx] = dtloc;
#   if  NTHREADS_TIME >=  256
    dttmp = dtmin[tidx ^  4];	 dtloc = FMIN(dtloc, dttmp);	dtmin[tidx] = dtloc;
#   if  NTHREADS_TIME >=  512
    dttmp = dtmin[tidx ^  8];	 dtloc = FMIN(dtloc, dttmp);	dtmin[tidx] = dtloc;
#   if  NTHREADS_TIME == 1024
    dttmp = dtmin[tidx ^ 16];	 dtloc = FMIN(dtloc, dttmp);	dtmin[tidx] = dtloc;
#endif//NTHREADS_TIME == 1024
#endif//NTHREADS_TIME >=  512
#endif//NTHREADS_TIME >=  256
#endif//NTHREADS_TIME >=  128
#endif//NTHREADS_TIME >=   64
    dtloc = dtmin[0];
#endif//USE_WARP_SHUFFLE_FUNC_TIME
  }/* if( tidx < (NTHREADS_TIME >> 5) ){ */


  if( tidx == 0 )
    *dt = LDEXP(UNITY, (int)FLOOR(LOG2(eta * dtloc)));
}


/**
 * @fn setTimeStep_kernel
 *
 * @brief Set time step for N-body particles.
 *
 * @sa setTimeStep_kernel
 */
extern "C"
void setTimeStep_dev(const int Ni, iparticle ibody, const real eta, const real eps, real *dt_dev, double *dt_hst
#ifndef SERIALIZED_EXECUTION
		     , MPIcfg_tree mpi
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
		     , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		     )
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  setTimeStep_kernel<<<1, NTHREADS_TIME>>>
    (Ni,
#ifndef OMIT_VELOCITY_FOR_TIME_STEP
     ibody.vx, ibody.vy, ibody.vz,
#endif//OMIT_VELOCITY_FOR_TIME_STEP
     ibody.acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     ibody.acc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     eta, eps, dt_dev);
  getLastCudaError("setTimeStep_kernel");

  real dt_tmp;
  checkCudaErrors(cudaMemcpy(&dt_tmp, dt_dev, sizeof(real), cudaMemcpyDeviceToHost));
  *dt_hst = (double)dt_tmp;

#ifndef SERIALIZED_EXECUTION
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, dt_hst, 1, MPI_DOUBLE, MPI_MIN, mpi.comm));
#endif//SERIALIZED_EXECUTION

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->setTimeStep_dev));
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}
#endif//BLOCK_TIME_STEP


#ifdef  BLOCK_TIME_STEP
/**
 * @fn prediction_kernel
 *
 * @brief Time integration (prediction phase: pi, vi, ai --> pj, vj).
 */
__global__ void prediction_kernel
(const int Nj, const double tnew,
 READ_ONLY position * RESTRICT ipos, READ_ONLY velocity * RESTRICT ivel, READ_ONLY ibody_time * RESTRICT time, READ_ONLY acceleration * RESTRICT iacc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 READ_ONLY acceleration * RESTRICT iacc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
 position * RESTRICT jpos, velocity * RESTRICT jvel)
{
  const int jj = GLOBALIDX_X1D;
  position pj = {ZERO, ZERO, ZERO, ZERO};
  velocity vj = {ZERO, ZERO, ZERO, ZERO};

  if( jj < Nj ){
    /** load information of all i-particles */
    pj = ipos[jj];
    vj = ivel[jj];
#ifndef SET_EXTERNAL_POTENTIAL_FIELD
    const acceleration aj = iacc[jj];
#else///SET_EXTERNAL_POTENTIAL_FIELD
    const acceleration slf = iacc    [jj];
    const acceleration ext = iacc_ext[jj];
    const acceleration aj = {slf.x + ext.x, slf.y + ext.y, slf.z + ext.z, slf.pot + ext.pot};
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    const real dt = (real)(tnew - time[jj].t0);
    const real dt_2 = HALF * dt;

    /** predict position and velocity of j-particle based on 2nd-order Runge-Kutta scheme */
    vj.x += dt_2 * aj.x;    pj.x += dt * vj.x;
    vj.y += dt_2 * aj.y;    pj.y += dt * vj.y;
    vj.z += dt_2 * aj.z;    pj.z += dt * vj.z;
  }/* if( jj < Nj ){ */
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//__CUDA_ARCH__ >= 700

  jpos[jj] = pj;
  jvel[jj] = vj;
}


#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
__device__ __forceinline__ double getMinimumDblTsub(const double  min)
#else///USE_WARP_SHUFFLE_FUNC_TIME
__device__ __forceinline__ void   getMinimumDblTsub(      double *min, volatile double * smem, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC_TIME
{
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
  /** 32 threads in a warp call this function at the same time */
  union {int2 i; double d;} val, tmp;  val.d = min;
#   if  DIV_NWARP(TSUB) >=  2
  tmp.i.x = __SHFL_XOR(SHFL_MASK_32, val.i.x,  1, DIV_NWARP(TSUB));  tmp.i.y = __SHFL_XOR(SHFL_MASK_32, val.i.y,  1, DIV_NWARP(TSUB));  val.d = fmin(val.d, tmp.d);
#   if  DIV_NWARP(TSUB) >=  4
  tmp.i.x = __SHFL_XOR(SHFL_MASK_32, val.i.x,  2, DIV_NWARP(TSUB));  tmp.i.y = __SHFL_XOR(SHFL_MASK_32, val.i.y,  2, DIV_NWARP(TSUB));  val.d = fmin(val.d, tmp.d);
#   if  DIV_NWARP(TSUB) >=  8
  tmp.i.x = __SHFL_XOR(SHFL_MASK_32, val.i.x,  4, DIV_NWARP(TSUB));  tmp.i.y = __SHFL_XOR(SHFL_MASK_32, val.i.y,  4, DIV_NWARP(TSUB));  val.d = fmin(val.d, tmp.d);
#   if  DIV_NWARP(TSUB) >= 16
  tmp.i.x = __SHFL_XOR(SHFL_MASK_32, val.i.x,  8, DIV_NWARP(TSUB));  tmp.i.y = __SHFL_XOR(SHFL_MASK_32, val.i.y,  8, DIV_NWARP(TSUB));  val.d = fmin(val.d, tmp.d);
#   if  DIV_NWARP(TSUB) == 32
  tmp.i.x = __SHFL_XOR(SHFL_MASK_32, val.i.x, 16, DIV_NWARP(TSUB));  tmp.i.y = __SHFL_XOR(SHFL_MASK_32, val.i.y, 16, DIV_NWARP(TSUB));  val.d = fmin(val.d, tmp.d);
#endif//DIV_NWARP(TSUB) == 32
#endif//DIV_NWARP(TSUB) >= 16
#endif//DIV_NWARP(TSUB) >=  8
#endif//DIV_NWARP(TSUB) >=  4
#endif//DIV_NWARP(TSUB) >=  2
  tmp.i.x = __SHFL(SHFL_MASK_32, val.i.x, 0, DIV_NWARP(TSUB));
  tmp.i.y = __SHFL(SHFL_MASK_32, val.i.y, 0, DIV_NWARP(TSUB));
  return (tmp.d);
#else///USE_WARP_SHUFFLE_FUNC_TIME
  smem[tidx] = *min;
#   if  DIV_NWARP(TSUB) >=  2
  double tmp;
  tmp = smem[tidx ^  1];  *min = fmin(*min, tmp);  smem[tidx] = *min;
#   if  DIV_NWARP(TSUB) >=  4
  tmp = smem[tidx ^  2];  *min = fmin(*min, tmp);  smem[tidx] = *min;
#   if  DIV_NWARP(TSUB) >=  8
  tmp = smem[tidx ^  4];  *min = fmin(*min, tmp);  smem[tidx] = *min;
#   if  DIV_NWARP(TSUB) >= 16
  tmp = smem[tidx ^  8];  *min = fmin(*min, tmp);  smem[tidx] = *min;
#   if  DIV_NWARP(TSUB) >= 32
  tmp = smem[tidx ^ 16];  *min = fmin(*min, tmp);  smem[tidx] = *min;
#endif//DIV_NWARP(TSUB) >= 32
#endif//DIV_NWARP(TSUB) >= 16
#endif//DIV_NWARP(TSUB) >=  8
#endif//DIV_NWARP(TSUB) >=  4
#endif//DIV_NWARP(TSUB) >=  2
  *min = smem[head];
#endif//USE_WARP_SHUFFLE_FUNC_TIME
}


/**
 * @fn setParticleTime
 *
 * @brief Calculate particle time step.
 */
__device__ __forceinline__ real setParticleTime
(const acceleration ai,
#ifndef OMIT_VELOCITY_FOR_TIME_STEP
 const velocity vi,
#endif//OMIT_VELOCITY_FOR_TIME_STEP
 const real eps, const real eta)
{
  /** estimate the required time step to resolve eps */
#ifndef OMIT_VELOCITY_FOR_TIME_STEP
  const real v2 = FLT_MIN + vi.x * vi.x + vi.y * vi.y + vi.z * vi.z;  const real vdt =      eps * RSQRT(v2);
#endif//OMIT_VELOCITY_FOR_TIME_STEP
  const real a2 = FLT_MIN + ai.x * ai.x + ai.y * ai.y + ai.z * ai.z;  const real adt = SQRT(eps * RSQRT(a2));

  /** set new time step */
#ifndef OMIT_VELOCITY_FOR_TIME_STEP
  return (LDEXP(UNITY, (int)FLOOR(LOG2(eta * FMIN(vdt, adt)))));
#else///OMIT_VELOCITY_FOR_TIME_STEP
  return (LDEXP(UNITY, (int)FLOOR(LOG2(eta * adt))));
#endif//OMIT_VELOCITY_FOR_TIME_STEP
}


#ifdef  SET_SINK_PARTICLES
extern "C"
void set_sink_particle_dev(const int Nsink, const sinkparticle sink, ulong *tag_hst)
{
  __NOTE__("%s\n", "start");
  checkCudaErrors(cudaMemcpy(sink.tag, tag_hst, sizeof(ulong) * Nsink, cudaMemcpyHostToDevice));
  __NOTE__("%s\n", "end");
}


__global__ void correct_MBH_kernel
(const int Nbh, READ_ONLY int * RESTRICT list, position * RESTRICT pos, velocity * RESTRICT vel, velocity * RESTRICT BHsum, READ_ONLY real * RESTRICT lmax2,
 READ_ONLY position * RESTRICT ipos, READ_ONLY velocity * RESTRICT ivel, READ_ONLY acceleration * RESTRICT iacc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 READ_ONLY acceleration * RESTRICT iacc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
 READ_ONLY position * RESTRICT jpos, READ_ONLY velocity * RESTRICT jvel
)
{
  const int tidx = THREADIDX_X1D;

  if( tidx < Nbh ){
    const int idx = list[tidx];
    pos[tidx] = jpos[idx];

    velocity vinew = jvel[idx];
#ifndef SET_EXTERNAL_POTENTIAL_FIELD
    const acceleration ai = iacc[idx];
#else///SET_EXTERNAL_POTENTIAL_FIELD
    const acceleration slf = iacc    [idx];
    const acceleration ext = iacc_ext[idx];
    const acceleration ai = {slf.x + ext.x, slf.y + ext.y, slf.z + ext.z, slf.pot + ext.pot};
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    vinew.dt *= HALF;
    vinew.x += vinew.dt * ai.x;
    vinew.y += vinew.dt * ai.y;
    vinew.z += vinew.dt * ai.z;

    vinew.dt = lmax2[tidx];
    vel[tidx] = vinew;

    const velocity zero_mom = {ZERO, ZERO, ZERO, ZERO};
    BHsum[tidx] = zero_mom;
  }/* if( tidx < Nbh ){ */
}
__global__ void update_MBH_kernel(const int Nbh, READ_ONLY int * RESTRICT list, READ_ONLY velocity * RESTRICT BHsum, position * RESTRICT ipos, velocity * RESTRICT ivel, const double time)
{
  const int tidx = THREADIDX_X1D;

  if( tidx < Nbh ){
    const velocity mom = BHsum[tidx];
    if( mom.dt > ZERO ){/**< mom.dt contains accreted mass in this time step */
      const int idx = list[tidx];
      position pi = ipos[idx];
      velocity vi = ivel[idx];

      const real minv = UNITY / (FLT_MIN + pi.m + mom.dt);
      vi.x = (mom.x + pi.m * vi.x) * minv;
      vi.y = (mom.y + pi.m * vi.y) * minv;
      vi.z = (mom.z + pi.m * vi.z) * minv;
      pi.m += mom.dt;

      ipos[idx] = pi;
      ivel[idx] = vi;

      printf("%e\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", time, idx, mom.dt, pi.m, pi.x, pi.y, pi.z, vi.x, vi.y, vi.z);
    }/* if( mom.m > ZERO ){ */
  }/* if( tidx < Nbh ){ */
}

#ifdef  RESET_CENTER_OF_MASS
__global__ void recorrect_central_MBH_kernel(const int Nbh, READ_ONLY int * RESTRICT list, READ_ONLY ulong * RESTRICT pidx, position * RESTRICT pos, velocity * RESTRICT vel, const ulong central_MBH)
{
  const int tidx = THREADIDX_X1D;

  if( tidx < Nbh ){
    const int idx = list[tidx];
    const ulong tag = pidx[idx];
    if( tag == central_MBH ){
      position pi = pos[idx];
      pi.x = ZERO;
      pi.y = ZERO;
      pi.z = ZERO;
      pos[tidx] = pi;
      velocity vi = vel[idx];
      vi.x = ZERO;
      vi.y = ZERO;
      vi.z = ZERO;
      vel[idx] = vi;
    }/* if( tag == central_MBH ){ */
  }/* if( tidx < Nbh ){ */
}
#endif//RESET_CENTER_OF_MASS
#endif//SET_SINK_PARTICLES


/**
 * @fn correction_kernel
 *
 * @brief Time integration (correction phase: pj, vj, ai --> pi, vi).
 */
__global__ void correction_kernel
(const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, const real eps, const real eta,
 position * RESTRICT ipos, velocity * RESTRICT ivel, ibody_time * RESTRICT time, READ_ONLY acceleration * RESTRICT iacc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 READ_ONLY acceleration * RESTRICT iacc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  SET_SINK_PARTICLES
 const int Nsink, READ_ONLY int * RESTRICT sinkidx, READ_ONLY position * RESTRICT sinkpos, READ_ONLY velocity * RESTRICT sinkvel, velocity * RESTRICT BHsum,
#endif//SET_SINK_PARTICLES
 READ_ONLY position * RESTRICT jpos, READ_ONLY velocity * RESTRICT jvel,
 const int reuseTree)
{
  const int tidx = THREADIDX_X1D;
  const int lane    = tidx          & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);

#ifndef USE_WARP_SHUFFLE_FUNC_TIME
  const int head = tidx - lane;
  __shared__ double smem[NTHREADS];
#endif//USE_WARP_SHUFFLE_FUNC_TIME

  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//__CUDA_ARCH__ >= 700

  ibody_time ti = {0.0, DBL_MAX};

  if( lane < info.num ){
    const int idx = info.head + lane;

    /** load pj, vj, ti, and ai */
#ifndef SET_EXTERNAL_POTENTIAL_FIELD
    const acceleration ai = iacc[idx];
#else///SET_EXTERNAL_POTENTIAL_FIELD
    const acceleration slf = iacc    [idx];
    const acceleration ext = iacc_ext[idx];
    const acceleration ai = {slf.x + ext.x, slf.y + ext.y, slf.z + ext.z, slf.pot + ext.pot};
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    velocity vi = jvel[idx];
    ti = time[idx];

    /** store pi */
#ifndef SET_SINK_PARTICLES
    ipos[idx] = jpos[idx];
#else///SET_SINK_PARTICLES
    position pi = jpos[idx];
    ipos[idx] = pi;
#endif//SET_SINK_PARTICLES

    /** update vi */
    vi.dt *= HALF;
    vi.x += vi.dt * ai.x;
    vi.y += vi.dt * ai.y;
    vi.z += vi.dt * ai.z;

    /** set new time step */
    vi.dt = setParticleTime
      (ai,
#ifndef OMIT_VELOCITY_FOR_TIME_STEP
    vi,
#endif//OMIT_VELOCITY_FOR_TIME_STEP
    eps, eta);
    /** store vi */
    ivel[idx] = vi;

    /** store ti */
    ti.t0 = ti.t1;
    ti.t1 += (double)vi.dt;
    time[idx] = ti;

#ifdef  SET_SINK_PARTICLES
#pragma unroll
    for(int ii = 0; ii < Nsink; ii++){
      if( sinkidx[ii] != idx ){
	const position ps = sinkpos[ii];
	const velocity vs = sinkvel[ii];/**< vs.dt contains l_max^2 */
	const real dx = pi.x - ps.x;
	const real dy = pi.y - ps.y;
	const real dz = pi.z - ps.z;
	const real vx = vi.x - vs.x;
	const real vy = vi.y - vs.y;
	const real vz = vi.z - vs.z;
	const real lx = dy * vz - dz * vy;
	const real ly = dz * vx - dx * vz;
	const real lz = dx * vy - dy * vx;

	/** pick up (1) separation is less than epsilon, (2) specific angular momentum is less than l_max, and (3) particle is approaching to MBH */
	/** how to eliminate BH particle?? (is (x, v) < 0 sufficient??; if no, add r^2 > ZERO) */
	const bool accrete =
	  /* ((dx * dx + dy * dy + dz * dz) <= (eta * eps * eta * eps)) && */
	  ((dx * dx + dy * dy + dz * dz) <= ((vx * vx + vy * vy + vz * vz) * (vi.dt * vi.dt))) &&
	  ((lx * lx + ly * ly + lz * lz) <= vs.dt) &&
	  ((dx * vx + dy * vy + dz * vz) < ZERO);

	if( accrete ){
	  const real mass = pi.m;
	  pi.m = ZERO;
	  ipos[idx] = pi;
	  atomicAdd(&(BHsum[ii].dt), mass       );
	  atomicAdd(&(BHsum[ii].x ), mass * vi.x);
	  atomicAdd(&(BHsum[ii].y ), mass * vi.y);
	  atomicAdd(&(BHsum[ii].z ), mass * vi.z);
	}/* if( accrete ){ */
      }/* if( sinkidx[ii] != idx ){ */
    }/* for(int ii = 0; ii < Nsink; ii++){ */

#endif//SET_SINK_PARTICLES
  }/* if( lane < info.num ){ */
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//__CUDA_ARCH__ >= 700


  /** get minimum ti.t1 of this group (TSUB threads) */
  if( reuseTree ){
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
    double tmin = getMinimumDblTsub(ti.t1);
#else///USE_WARP_SHUFFLE_FUNC_TIME
    double tmin = ti.t1;
    getMinimumDblTsub(&tmin, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_TIME
    if( lane == 0 )
      laneTime[laneIdx] = tmin;
  }/* if( reuseTree ){ */
}


/**
 * @fn adjustParticleTime_kernel
 *
 * @brief Adjust particle time step.
 */
__global__ void adjustParticleTime_kernel
(const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, const real eps, const real eta,
 velocity * RESTRICT ivel, ibody_time * RESTRICT time, READ_ONLY acceleration * RESTRICT iacc
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 , READ_ONLY acceleration * RESTRICT iacc_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
)
{
  const int tidx = THREADIDX_X1D;
  const int lane    = tidx          & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);

#ifndef USE_WARP_SHUFFLE_FUNC_TIME
  const int head = tidx - lane;
  __shared__ double smem[NTHREADS];
#endif//USE_WARP_SHUFFLE_FUNC_TIME

  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//__CUDA_ARCH__ >= 700

  ibody_time ti = {0.0, DBL_MAX};

  if( lane < info.num ){
    const int idx = info.head + lane;

    /** load pj, vj, ti, and ai */
#ifndef SET_EXTERNAL_POTENTIAL_FIELD
    const acceleration ai = iacc[idx];
#else///SET_EXTERNAL_POTENTIAL_FIELD
    const acceleration slf = iacc    [idx];
    const acceleration ext = iacc_ext[idx];
    const acceleration ai = {slf.x + ext.x, slf.y + ext.y, slf.z + ext.z, slf.pot + ext.pot};
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    velocity vi = ivel[idx];
    ti = time[idx];

    /** set new time step */
    vi.dt = setParticleTime
      (ai,
#ifndef OMIT_VELOCITY_FOR_TIME_STEP
    vi,
#endif//OMIT_VELOCITY_FOR_TIME_STEP
    eps, eta);
    /** store vi */
    ivel[idx] = vi;

    /** store ti */
    ti.t1 = ti.t0 + (double)vi.dt;
    time[idx] = ti;
  }/* if( lane < info.num ){ */
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//__CUDA_ARCH__ >= 700

  /** get minimum ti.t1 of this group (TSUB threads) */
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
  double tmin = getMinimumDblTsub(ti.t1);
#else//USE_WARP_SHUFFLE_FUNC_TIME
  double tmin = ti.t1;
  getMinimumDblTsub(&tmin, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_TIME
  if( lane == 0 )
    laneTime[laneIdx] = tmin;
}


/**
 * @fn setLaneTime_kernel
 *
 * @brief Set time step of group of N-body particles.
 */
__global__ void setLaneTime_kernel(const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, READ_ONLY ibody_time * RESTRICT time)
{
  const int tidx = THREADIDX_X1D;
  const int lane    = tidx          & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);

#ifndef USE_WARP_SHUFFLE_FUNC_TIME
  const int head = tidx - lane;
  __shared__ double smem[NTHREADS];
#endif//USE_WARP_SHUFFLE_FUNC_TIME

  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//__CUDA_ARCH__ >= 700

  ibody_time ti = {0.0, DBL_MAX};

  if( lane < info.num )
    ti = time[info.head + lane];
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//__CUDA_ARCH__ >= 700

  /** get minimum ti.t1 of this group (TSUB threads) */
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
  double tmin = getMinimumDblTsub(ti.t1);
#else///USE_WARP_SHUFFLE_FUNC_TIME
  double tmin = ti.t1;
  getMinimumDblTsub(&tmin, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_TIME
  if( lane == 0 )
    laneTime[laneIdx] = tmin;
}


/**
 * @fn prediction_dev
 *
 * @brief Time integration (prediction phase: pi, vi, ai --> pj, vj).
 *
 * @sa prediction_kernel
 */
extern "C"
void prediction_dev(const int Nj, const double tnew, const iparticle pi
#   if  defined(SET_SINK_PARTICLES) && defined(RESET_CENTER_OF_MASS)
		    , const int Nbh, const sinkparticle bh, const bool replace_central_bh, const ulong central_MBH_idx
#endif//defined(SET_SINK_PARTICLES) && defined(RESET_CENTER_OF_MASS)
#ifdef  EXEC_BENCHMARK
		    , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		    )
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  int Nrem = BLOCKSIZE(Nj, NTHREADS_TIME);

  if( Nrem <= MAX_BLOCKS_PER_GRID )
    prediction_kernel<<<Nrem, NTHREADS_TIME>>>
      (Nj, tnew, pi.pos, pi.vel, pi.time, pi.acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
       pi.acc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
       pi.jpos, pi.jvel);
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_TIME;
      prediction_kernel<<<Nblck, NTHREADS_TIME>>>
	(Nsub, tnew, &pi.pos[hidx], &pi.vel[hidx], &pi.time[hidx], &pi.acc[hidx],
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
	 &pi.acc_ext[hidx],
#endif//SET_EXTERNAL_POTENTIAL_FIELD
	 &pi.jpos[hidx], &pi.jvel[hidx]);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("prediction_kernel");

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->prediction_dev));
#endif//EXEC_BENCHMARK

#   if  defined(SET_SINK_PARTICLES) && defined(RESET_CENTER_OF_MASS)
  if( replace_central_bh ){
    /** place central BH at the center with null velocity */
    recorrect_central_MBH_kernel<<<1, Nbh>>>(Nbh, bh.list, pi.idx, pi.jpos, pi.jvel, central_MBH_idx);
    getLastCudaError("recorrect_central_MBH_kernel");
  }/* if( replace_central_bh ){ */
#endif//defined(SET_SINK_PARTICLES) && defined(RESET_CENTER_OF_MASS)


  __NOTE__("%s\n", "end");
}


/**
 * @fn correction_dev
 *
 * @brief Time integration (correction phase: pj, vj, ai --> pi, vi).
 *
 * @sa correction_kernel
 */
extern "C"
void correction_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, const real eps, const real eta, const iparticle pi, const int reuseTree
#ifdef  SET_SINK_PARTICLES
		    , const int Nbh, const sinkparticle bh, const double time
#ifdef  RESET_CENTER_OF_MASS
		    , const bool replace_central_bh, const ulong central_MBH_idx
#endif//RESET_CENTER_OF_MASS
#endif//SET_SINK_PARTICLES
#ifdef  EXEC_BENCHMARK
		    , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		    )
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

#ifdef  SET_SINK_PARTICLES
  /** remember BH position and velocity in the previous and current time step */
  correct_MBH_kernel<<<1, Nbh>>>(Nbh, bh.list, bh.pos, bh.vel, bh.mom, bh.lmax2,
				 pi.pos, pi.vel, pi.acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
				 pi.acc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
				 pi.jpos, pi.jvel);
  getLastCudaError("correct_MBH_kernel");

#ifdef  RESET_CENTER_OF_MASS
  if( replace_central_bh ){
    /** place central BH at the center with null velocity */
    recorrect_central_MBH_kernel<<<1, Nbh>>>(Nbh, bh.list, pi.idx, bh.pos, bh.vel, central_MBH_idx);
    getLastCudaError("recorrect_central_MBH_kernel");
  }/* if( replace_central_bh ){ */
#endif//RESET_CENTER_OF_MASS
#endif//SET_SINK_PARTICLES

  /** thread-block structure must be identical to tree traversal */
#ifndef SERIALIZED_EXECUTION
  if( Ngrp > 0 )
#endif//SERIALIZED_EXECUTION
    {
      int Nrem = BLOCKSIZE(Ngrp, NWARP * NGROUPS);

      if( Nrem <= MAX_BLOCKS_PER_GRID )
	correction_kernel<<<Nrem, NTHREADS>>>
	  (BLOCKSIZE(Ngrp, NGROUPS) * NGROUPS, laneInfo, laneTime, eps, eta, pi.pos, pi.vel, pi.time, pi.acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
	   pi.acc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  SET_SINK_PARTICLES
	   Nbh, bh.list, bh.pos, bh.vel, bh.mom,
#endif//SET_SINK_PARTICLES
	   pi.jpos, pi.jvel, reuseTree);
      else{
	const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
	int hidx = 0;

	for(int iter = 0; iter < Niter; iter++){
	  int Nblck = MAX_BLOCKS_PER_GRID;
	  if( Nblck > Nrem )	    Nblck = Nrem;

	  int Nsub = Nblck * NWARP * NGROUPS;
	  correction_kernel<<<Nblck, NTHREADS>>>
	    (BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx], &laneTime[hidx], eps, eta, pi.pos, pi.vel, pi.time, pi.acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
	     pi.acc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  SET_SINK_PARTICLES
	     Nbh, bh.list, bh.pos, bh.vel, bh.mom,
#endif//SET_SINK_PARTICLES
	     pi.jpos, pi.jvel, reuseTree);

	  hidx += Nsub;
	  Nrem -= Nblck;
	}/* for(int iter = 0; iter < Niter; iter++){ */
      }/* else{ */

      getLastCudaError("correction_kernel");
    }

#ifdef  SET_SINK_PARTICLES
  /** update BH mass and velocity */
  update_MBH_kernel<<<1, Nbh>>>(Nbh, bh.list, bh.mom, pi.pos, pi.vel, time);
  getLastCudaError("update_MBH_kernel");

#ifdef  RESET_CENTER_OF_MASS
  if( replace_central_bh ){
    /** place central BH at the center with null velocity */
    recorrect_central_MBH_kernel<<<1, Nbh>>>(Nbh, bh.list, pi.idx, pi.pos, pi.vel, central_MBH_idx);
    getLastCudaError("recorrect_central_MBH_kernel");
  }/* if( replace_central_bh ){ */
#endif//RESET_CENTER_OF_MASS
#endif//SET_SINK_PARTICLES

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->correction_dev));
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}


/**
 * @fn adjustParticleTime_dev
 *
 * @brief Adjust particle time step.
 *
 * @sa adjustParticleTime_kernel
 */
extern "C"
void adjustParticleTime_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, const real eps, const real eta, const iparticle pi
#ifdef  EXEC_BENCHMARK
			    , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
			    )
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  int Nrem = BLOCKSIZE(Ngrp, NWARP * NGROUPS);

  if( Nrem <= MAX_BLOCKS_PER_GRID )
    adjustParticleTime_kernel<<<Nrem, NTHREADS>>>
      (BLOCKSIZE(Ngrp, NGROUPS) * NGROUPS, laneInfo, laneTime, eps, eta, pi.vel, pi.time, pi.acc
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
       , pi.acc_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
       );
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NWARP * NGROUPS;
      adjustParticleTime_kernel<<<Nblck, NTHREADS>>>
	(BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx], &laneTime[hidx], eps, eta, pi.vel, pi.time, pi.acc
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
	 , pi.acc_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
	 );

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("adjustParticleTime_kernel");

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->adjustParticleTime_dev));
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}


/**
 * @fn setLaneTime_dev
 *
 * @brief Set time step of group of N-body particles.
 *
 * @sa setLaneTime_kernel
 */
extern "C"
void setLaneTime_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, double * RESTRICT laneTime, const iparticle pi
#ifdef  EXEC_BENCHMARK
		     , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		     )
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  int Nrem = BLOCKSIZE(Ngrp, NWARP * NGROUPS);

  if( Nrem <= MAX_BLOCKS_PER_GRID )
    setLaneTime_kernel<<<Nrem, NTHREADS>>>(BLOCKSIZE(Ngrp, NGROUPS) * NGROUPS, laneInfo, laneTime, pi.time);
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NWARP * NGROUPS;
      setLaneTime_kernel<<<Nblck, NTHREADS>>>(BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx], &laneTime[hidx], pi.time);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("setLaneTime_kernel");

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->setLaneTime_dev));
#else///EXEC_BENCHMARK
  cudaDeviceSynchronize();
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}
#else///BLOCK_TIME_STEP
/**
 * @fn advPos_kernel
 *
 * @brief Orbit integral for position.
 */
__global__ void advPos_kernel
(const int Ni, position * RESTRICT ipos, READ_ONLY real * RESTRICT vx, READ_ONLY real * RESTRICT vy, READ_ONLY real * RESTRICT vz, const real dt)
{
  const int ii = GLOBALIDX_X1D;
  position pi = {ZERO, ZERO, ZERO, ZERO};

  if( ii < Ni ){
    /** load an i-particle */
    pi = ipos[ii];

    pi.x += dt * vx[ii];
    pi.y += dt * vy[ii];
    pi.z += dt * vz[ii];
  }/* if( ii < Ni ){ */
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//__CUDA_ARCH__ >= 700

  ipos[ii] = pi;
}


/**
 * @fn advVel_kernel
 *
 * @brief Orbit integral for velocity.
 */
__global__ void advVel_kernel
(const int Ni, READ_ONLY acceleration * RESTRICT iacc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 READ_ONLY acceleration * RESTRICT iacc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
 real * RESTRICT vx, real * RESTRICT vy, real * RESTRICT vz, const real dt)
{
  const int ii = GLOBALIDX_X1D;

  if( ii < Ni ){
    /** load acceleration of i-particle */
#ifndef SET_EXTERNAL_POTENTIAL_FIELD
    const acceleration ai = iacc[ii];
#else///SET_EXTERNAL_POTENTIAL_FIELD
    const acceleration slf = iacc    [ii];
    const acceleration ext = iacc_ext[ii];
    const acceleration ai = {slf.x + ext.x, slf.y + ext.y, slf.z + ext.z, slf.pot + ext.pot};
#endif//SET_EXTERNAL_POTENTIAL_FIELD

    /** update velocity */
    vx[ii] += dt * ai.x;
    vy[ii] += dt * ai.y;
    vz[ii] += dt * ai.z;
  }/* if( ii < Ni ){ */
}


/**
 * @fn advPos_dev
 *
 * @brief Orbit integral for position.
 *
 * @sa advPos_kernel
 */
extern "C"
void advPos_dev(const int Ni, iparticle ibody, const real dt
#ifdef  EXEC_BENCHMARK
		, wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		)
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  int Nrem = BLOCKSIZE(Ni, NTHREADS_TIME);

  if( Nrem <= MAX_BLOCKS_PER_GRID )
    advPos_kernel<<<Nrem, NTHREADS_TIME>>>(Ni, ibody.pos, ibody.vx, ibody.vy, ibody.vz, dt);
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_TIME;
      advPos_kernel<<<Nblck, NTHREADS_TIME>>>(Nsub, &ibody.pos[hidx], &ibody.vx[hidx], &ibody.vy[hidx], &ibody.vz[hidx], dt);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("advPos_kernel");

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->advPos_dev));
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}


/**
 * @fn advVel_dev
 *
 * @brief Orbit integral for velocity.
 *
 * @sa advVel_kernel
 */
extern "C"
void advVel_dev(const int Ni, iparticle ibody, const real dt
#ifdef  EXEC_BENCHMARK
		, wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		)
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  int Nrem = BLOCKSIZE(Ni, NTHREADS_TIME);

  if( Nrem <= MAX_BLOCKS_PER_GRID )
    advVel_kernel<<<Nrem, NTHREADS_TIME>>>
      (Ni, ibody.acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
       ibody.acc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
       ibody.vx, ibody.vy, ibody.vz, dt);
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_TIME;
      advVel_kernel<<<Nblck, NTHREADS_TIME>>>
	(Nsub, &ibody.acc[hidx],
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
	 &ibody.acc_ext[hidx],
#endif//SET_EXTERNAL_POTENTIAL_FIELD
	 &ibody.vx[hidx], &ibody.vy[hidx], &ibody.vz[hidx], dt);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("advVel_kernel");

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->advVel_dev));
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}
#endif//BLOCK_TIME_STEP


/**
 * @fn copyParticle_hst2dev
 *
 * @brief Copy N-body particles from host to device.
 */
extern "C"
void copyParticle_hst2dev(const int Ni, iparticle hst, iparticle dev
#ifdef  EXEC_BENCHMARK
			  , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
			  )
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  /** send i-particles from the host to the device using the default CUDA stream */
  checkCudaErrors(cudaMemcpy(dev. idx, hst. idx, sizeof(       ulong) * Ni, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(dev. pos, hst. pos, sizeof(    position) * Ni, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(dev. acc, hst. acc, sizeof(acceleration) * Ni, cudaMemcpyHostToDevice));
#ifdef  BLOCK_TIME_STEP
  checkCudaErrors(cudaMemcpy(dev. vel, hst. vel, sizeof(    velocity) * Ni, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(dev.time, hst.time, sizeof(  ibody_time) * Ni, cudaMemcpyHostToDevice));
#else///BLOCK_TIME_STEP
  checkCudaErrors(cudaMemcpy(dev.  vx, hst.  vx, sizeof(        real) * Ni, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(dev.  vy, hst.  vy, sizeof(        real) * Ni, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(dev.  vz, hst.  vz, sizeof(        real) * Ni, cudaMemcpyHostToDevice));
#endif//BLOCK_TIME_STEP

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  checkCudaErrors(cudaMemcpy(dev.acc_ext, hst.acc_ext, sizeof(acceleration) * Ni, cudaMemcpyHostToDevice));
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->copyParticle_hst2dev));
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}


/**
 * @fn copyParticle_dev2hst
 *
 * @brief Copy N-body particles from device to host
 */
extern "C"
void copyParticle_dev2hst(const int Ni, iparticle dev, iparticle hst
#ifdef  EXEC_BENCHMARK
			  , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
			  )
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  /** send i-particles from the device to the host using the default CUDA stream */
  checkCudaErrors(cudaMemcpy(hst. idx, dev. idx, sizeof(       ulong) * Ni, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(hst. pos, dev. pos, sizeof(    position) * Ni, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(hst. acc, dev. acc, sizeof(acceleration) * Ni, cudaMemcpyDeviceToHost));
#ifdef  BLOCK_TIME_STEP
  checkCudaErrors(cudaMemcpy(hst. vel, dev. vel, sizeof(    velocity) * Ni, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(hst.time, dev.time, sizeof(  ibody_time) * Ni, cudaMemcpyDeviceToHost));
#else///BLOCK_TIME_STEP
  checkCudaErrors(cudaMemcpy(hst.  vx, dev.  vx, sizeof(        real) * Ni, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(hst.  vy, dev.  vy, sizeof(        real) * Ni, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(hst.  vz, dev.  vz, sizeof(        real) * Ni, cudaMemcpyDeviceToHost));
#endif//BLOCK_TIME_STEP

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  checkCudaErrors(cudaMemcpy(hst.acc_ext, dev.acc_ext, sizeof(acceleration) * Ni, cudaMemcpyDeviceToHost));
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->copyParticle_dev2hst));
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}


/**
 * @fn copyParticleAsync_hst2dev
 *
 * @brief Copy N-body particles from host to device.
 */
extern "C"
void copyParticleAsync_hst2dev(const int Ni, iparticle hst, iparticle dev, cudaStream_t stream)
{
  __NOTE__("%s\n", "start");

  /** send i-particles from the host to the device using a CUDA stream */
  checkCudaErrors(cudaMemcpyAsync(dev. idx, hst. idx, sizeof(       ulong) * Ni, cudaMemcpyHostToDevice, stream));
  checkCudaErrors(cudaMemcpyAsync(dev. pos, hst. pos, sizeof(    position) * Ni, cudaMemcpyHostToDevice, stream));
  checkCudaErrors(cudaMemcpyAsync(dev. acc, hst. acc, sizeof(acceleration) * Ni, cudaMemcpyHostToDevice, stream));
#ifdef  BLOCK_TIME_STEP
  checkCudaErrors(cudaMemcpyAsync(dev. vel, hst. vel, sizeof(    velocity) * Ni, cudaMemcpyHostToDevice, stream));
  checkCudaErrors(cudaMemcpyAsync(dev.time, hst.time, sizeof(  ibody_time) * Ni, cudaMemcpyHostToDevice, stream));
#else///BLOCK_TIME_STEP
  checkCudaErrors(cudaMemcpyAsync(dev.  vx, hst.  vx, sizeof(        real) * Ni, cudaMemcpyHostToDevice, stream));
  checkCudaErrors(cudaMemcpyAsync(dev.  vy, hst.  vy, sizeof(        real) * Ni, cudaMemcpyHostToDevice, stream));
  checkCudaErrors(cudaMemcpyAsync(dev.  vz, hst.  vz, sizeof(        real) * Ni, cudaMemcpyHostToDevice, stream));
#endif//BLOCK_TIME_STEP

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  checkCudaErrors(cudaMemcpyAsync(dev.acc_ext, hst.acc_ext, sizeof(acceleration) * Ni, cudaMemcpyHostToDevice, stream));
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
}


/**
 * @fn copyParticleAsync_dev2hst
 *
 * @brief Copy N-body particles from device to host.
 */
extern "C"
void copyParticleAsync_dev2hst(const int Ni, iparticle dev, iparticle hst, cudaStream_t stream)
{
  __NOTE__("%s\n", "start");

  /** send i-particles from the device to the host using the default CUDA stream */
  checkCudaErrors(cudaMemcpyAsync(hst. idx, dev. idx, sizeof(       ulong) * Ni, cudaMemcpyDeviceToHost, stream));
  checkCudaErrors(cudaMemcpyAsync(hst. pos, dev. pos, sizeof(    position) * Ni, cudaMemcpyDeviceToHost, stream));
  checkCudaErrors(cudaMemcpyAsync(hst. acc, dev. acc, sizeof(acceleration) * Ni, cudaMemcpyDeviceToHost, stream));
#ifdef  BLOCK_TIME_STEP
  checkCudaErrors(cudaMemcpyAsync(hst. vel, dev. vel, sizeof(    velocity) * Ni, cudaMemcpyDeviceToHost, stream));
  checkCudaErrors(cudaMemcpyAsync(hst.time, dev.time, sizeof(  ibody_time) * Ni, cudaMemcpyDeviceToHost, stream));
#else///BLOCK_TIME_STEP
  checkCudaErrors(cudaMemcpyAsync(hst.  vx, dev.  vx, sizeof(        real) * Ni, cudaMemcpyDeviceToHost, stream));
  checkCudaErrors(cudaMemcpyAsync(hst.  vy, dev.  vy, sizeof(        real) * Ni, cudaMemcpyDeviceToHost, stream));
  checkCudaErrors(cudaMemcpyAsync(hst.  vz, dev.  vz, sizeof(        real) * Ni, cudaMemcpyDeviceToHost, stream));
#endif//BLOCK_TIME_STEP

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  checkCudaErrors(cudaMemcpyAsync(hst.acc_ext, dev.acc_ext, sizeof(acceleration) * Ni, cudaMemcpyDeviceToHost, stream));
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
}


#ifdef  COMPARE_WITH_DIRECT_SOLVER
/**
 * @fn copyAccel_dev2hst
 *
 * @brief Copy particle acceleration from device to host.
 */
extern "C"
void copyAccel_dev2hst(const int Ni, acceleration * RESTRICT dev, acceleration * RESTRICT hst)
{
  __NOTE__("%s\n", "start");

  checkCudaErrors(cudaMemcpy(hst, dev, sizeof(acceleration) * Ni, cudaMemcpyDeviceToHost));

  __NOTE__("%s\n", "end");
}
#endif//COMPARE_WITH_DIRECT_SOLVER


#ifdef  COUNT_INTERACTIONS
/**
 * @fn copyCounters_dev2hst
 *
 * @brief Copy metric counters from device to host.
 */
extern "C"
void copyCounters_dev2hst(const int Ni, iparticle_treeinfo dev, iparticle_treeinfo hst)
{
  __NOTE__("%s\n", "start");

  checkCudaErrors(cudaMemcpy(hst.  Nj, dev.  Nj, sizeof(int) * Ni, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(hst.Nbuf, dev.Nbuf, sizeof(int) * Ni, cudaMemcpyDeviceToHost));

  __NOTE__("%s\n", "end");
}
#endif//COUNT_INTERACTIONS


#   if  GPUVER >= 70
extern "C"
void setCacheConfig_adv_dev_cu(void)
{
  __NOTE__("%s\n", "start");

  /* remove shared memory if __global__ function does not use */
#ifdef  BLOCK_TIME_STEP
  checkCudaErrors(cudaFuncSetAttribute(adjustTimeStep_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(prediction_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
  checkCudaErrors(cudaFuncSetAttribute(correction_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(adjustParticleTime_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(setLaneTime_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
#endif//USE_WARP_SHUFFLE_FUNC_TIME
#else///BLOCK_TIME_STEP
  checkCudaErrors(cudaFuncSetAttribute(advPos_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(advVel_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
#endif//BLOCK_TIME_STEP

  __NOTE__("%s\n", "end");
}
#endif//GPUVER >= 70
