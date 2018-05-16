/**
 * @file shrink_dev.cu
 *
 * @brief Source code to split i-particle groups on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/05/08 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

/**
 * @def OUTPUT_AUTO_TUNING_PARAMETER
 *
 * @brief a switch to output values in Brent's method
 */
/* #define OUTPUT_AUTO_TUNING_PARAMETER */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <helper_cuda.h>

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

#ifdef  MPI_MAX_FOR_RMAX_IN_AUTO_TUNING
#include <mpi.h>
#include "mpilib.h"
/* #define USE_FIXED_RMAX_VALUE */
#endif//MPI_MAX_FOR_RMAX_IN_AUTO_TUNING

#include "macro.h"
#include "cudalib.h"

#include "../misc/structure.h"
#include "../misc/device.h"

#include "make.h"
#include "walk_dev.h"
#include "make_dev.h"

#include "../misc/brent.h"
#include "neighbor_dev.h"

#include "shrink_dev.h"


static const laneinfo nullInfo = {NUM_BODY_MAX, 0};


/* #define NEIGHBOR_NUM_LANE (TSUB / NWARP) */
#define NEIGHBOR_NUM_LANE (DIV_NWARP(TSUB))


/**
 * @fn initLaneTime_kernel
 *
 * @brief Initialize time step of particle group.
 */
__global__ void initLaneTime_kernel(const int num, double *laneTime)
{
  const int idx = GLOBALIDX_X1D;

  if( idx < num )
    laneTime[idx] = DBL_MAX;
}


/**
 * @fn freeParticleGroups
 *
 * @brief Memory deallocation for particle groups.
 */
extern "C"
void freeParticleGroups
(laneinfo  *laneInfo_hst, laneinfo  *laneInfo_dev, double  *laneTime_dev, int  *inum_hst, int  *inum_dev
#ifdef  SWITCH_WITH_J_PARALLELIZATION
 , const bool forLocal
#endif//SWITCH_WITH_J_PARALLELIZATION
)
{
  __NOTE__("%s\n", "start");

  mycudaFreeHost(laneInfo_hst);
  mycudaFree    (laneInfo_dev);

#ifdef  SWITCH_WITH_J_PARALLELIZATION
  if( forLocal )
#endif//SWITCH_WITH_J_PARALLELIZATION
    {
      mycudaFree    (laneTime_dev);
      mycudaFreeHost(    inum_hst);
      mycudaFree    (    inum_dev);
    }

  __NOTE__("%s\n", "end");
}


/**
 * @fn allocParticleGroups
 *
 * @brief Memory allocation for particle groups.
 *
 * @sa initLaneTime_kernel
 */
extern "C"
muse allocParticleGroups
(laneinfo **laneInfo_hst, laneinfo **laneInfo_dev, double **laneTime_dev, int **inum_hst, int **inum_dev,
#ifdef  SWITCH_WITH_J_PARALLELIZATION
 const bool forLocal,
#endif//SWITCH_WITH_J_PARALLELIZATION
 int *inumPerLane, int *maxNgrp, const int num_max)
{
  __NOTE__("%s\n", "start");


#   if  GPUGEN >= 70
  /* remove shared memory if __global__ function does not use */
  checkCudaErrors(cudaFuncSetCacheConfig(initLaneTime_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, 0));
#endif//GPUGEN >= 70

  muse alloc = {0, 0};

  /** number of i-particles per lane (a group of TSUB threads) */
  *inumPerLane = NEIGHBOR_NUM_LANE;
  *maxNgrp = BLOCKSIZE(num_max, *inumPerLane) * NUM_IGROUP_SAFETY_FACTOR;

  /** memory allocation for laneInfo */
  mycudaMallocHost((void **)laneInfo_hst, (size_t)(*maxNgrp) * sizeof(laneinfo));  alloc.host   += (size_t)(*maxNgrp) * sizeof(laneinfo);
  mycudaMalloc    ((void **)laneInfo_dev, (size_t)(*maxNgrp) * sizeof(laneinfo));  alloc.device += (size_t)(*maxNgrp) * sizeof(laneinfo);

  /** initialize laneInfo */
  for(int ii = 0; ii < *maxNgrp; ii++)
    (*laneInfo_hst)[ii] = nullInfo;
  checkCudaErrors(cudaMemcpy(*laneInfo_dev, *laneInfo_hst, (*maxNgrp) * sizeof(laneinfo), cudaMemcpyHostToDevice));

#ifdef  SWITCH_WITH_J_PARALLELIZATION
  if( forLocal )
#endif//SWITCH_WITH_J_PARALLELIZATION
    {
      /** memory allocation to estimate number of i-particle within a group */
      mycudaMallocHost((void **)inum_hst, (size_t)num_max * sizeof(int));      alloc.host   += (size_t)num_max * sizeof(int);
      mycudaMalloc    ((void **)inum_dev, (size_t)num_max * sizeof(int));      alloc.device += (size_t)num_max * sizeof(int);

      /** memory allocation for laneTime */
      mycudaMalloc    ((void **)laneTime_dev, (size_t)(*maxNgrp) * sizeof(double));
      alloc.device +=                         (size_t)(*maxNgrp) * sizeof(double);

      /** initialize laneTime */
      int Nrem = BLOCKSIZE(*maxNgrp, 1024);
      if( Nrem <= MAX_BLOCKS_PER_GRID )
	initLaneTime_kernel<<<Nrem, 1024>>>(*maxNgrp, *laneTime_dev);
      else{
	const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
	int hidx = 0;

	for(int iter = 0; iter < Niter; iter++){
	  int Nblck = MAX_BLOCKS_PER_GRID;
	  if( Nblck > Nrem )	Nblck = Nrem;

	  int Nsub = Nblck * 1024;
	  initLaneTime_kernel<<<Nblck, 1024>>>(Nsub, &((*laneTime_dev)[hidx]));

	  hidx += Nsub;
	  Nrem -= Nblck;
	}/* for(int iter = 0; iter < Niter; iter++){ */
      }/* else{ */
    }

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn countContinuousNeighbor_kernel
 *
 * @brief Count up continuous particles.
 */
__global__ void countContinuousNeighbor_kernel(const int Ni, position * RESTRICT ibody, const int inumPerLane, const real r2max, int * RESTRICT inum)
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;

  __shared__ position ipos[NTHREADS_SHRINK + NEIGHBOR_NUM_LANE];

  if( gidx < Ni ){
    /** load position of i-particles and neighbor candidates */
    ipos[tidx] = ibody[gidx];

    if( tidx < NEIGHBOR_NUM_LANE ){
      int idx = gidx + NTHREADS_SHRINK;      if( idx > (Ni - 1) )	idx = Ni - 1;
      ipos[NTHREADS_SHRINK + tidx] = ibody[idx];
    }/* if( tidx < NEIGHBOR_NUM_LANE ){ */
    __syncthreads();

    /** calculate distance with NEIGHBOR_NUM_LANE particles and remember the maximum */
    const position pi = ipos[tidx];
    int nmax = inumPerLane;

#pragma unroll
    for(int ii = 0; ii < NEIGHBOR_NUM_LANE; ii++){
      const position pj = ipos[tidx + 1 + ii];

      const real dx = pj.x - pi.x;
      const real dy = pj.y - pi.y;
      const real dz = pj.z - pi.z;
      const real r2 = FLT_MIN + dx * dx + dy * dy + dz * dz;

      if( r2 > r2max ){
	nmax = 1 + ii;
	break;
      }/* if( r2 > r2max ){ */
    }/* for(int ii = 0; ii < NEIGHBOR_NUM_LANE; ii++){ */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700

    /** commit the derived maximum number as neighbor particles within rmax */
    inum[gidx] = nmax;
  }/* if( gidx < Ni ){ */
}


/**
 * @fn examineParticleSeparation
 *
 * @brief Examine optimal separation for particle groups.
 *
 * @sa facileNeighborSearching_dev
 * @sa brentPerturb
 * @sa brentInit1st
 */
extern "C"
void examineParticleSeparation(const int Ni, iparticle body_dev, brentStatus *brent
#ifdef  MPI_MAX_FOR_RMAX_IN_AUTO_TUNING
			       , MPI_Comm mpi_comm
#endif//MPI_MAX_FOR_RMAX_IN_AUTO_TUNING
#ifdef  EXEC_BENCHMARK
			       , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
			       )
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  /** set a criterion to judge whether unify i-particles or not */
  facileNeighborSearching_dev(Ni, body_dev);
#ifdef  HUNT_FIND_PARAMETER
  stopStopwatch(&(elapsed->searchNeighbor_kernel));
  initStopwatch();
#endif//HUNT_FIND_PARAMETER

#ifndef DISABLE_AUTO_TUNING
#ifndef MPI_MAX_FOR_RMAX_IN_AUTO_TUNING
  const double rmax = CAST_R2D(thrust::reduce((thrust::device_ptr<real>)body_dev.neighbor, (thrust::device_ptr<real>)(body_dev.neighbor + Ni), ZERO, thrust::maximum<real>()));
  const double rmin = rmax * NEIGHBOR_LENGTH_SHRINK_FACTOR;
#else///MPI_MAX_FOR_RMAX_IN_AUTO_TUNING
#ifndef USE_FIXED_RMAX_VALUE
  double rmax = CAST_R2D(thrust::reduce((thrust::device_ptr<real>)body_dev.neighbor, (thrust::device_ptr<real>)(body_dev.neighbor + Ni), ZERO, thrust::maximum<real>()));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &rmax, 1, MPI_DOUBLE, MPI_MAX, mpi_comm));
  const double rmin = rmax * NEIGHBOR_LENGTH_SHRINK_FACTOR;
#else///USE_FIXED_RMAX_VALUE
  const double rmax = CAST_R2D(thrust::reduce((thrust::device_ptr<real>)body_dev.neighbor, (thrust::device_ptr<real>)(body_dev.neighbor + Ni), ZERO, thrust::maximum<real>()));
  const double rmin = rmax * NEIGHBOR_LENGTH_SHRINK_FACTOR;
#endif//USE_FIXED_RMAX_VALUE
#endif//MPI_MAX_FOR_RMAX_IN_AUTO_TUNING

#ifndef USE_FIXED_RMAX_VALUE
  if( brent->initialized ){
    if( (brent->u.pos >= rmin) && (brent->u.pos <= rmax) )
      brentPerturb(brent, rmin, rmax);
    else{
      brentInit1st(brent, rmin, rmax);
      brent->x = brent->u;
      brentInit2nd(brent);
    }/* else{ */
  }
  else{
    brent->initialized = true;
    brentInit1st(brent, rmin, rmax);
  }/* else{ */
#else///USE_FIXED_RMAX_VALUE
  if( !brent->initialized ){
    brent->initialized = true;
    brentInit1st(brent, rmin, rmax);
  }/* if( !brent->initialized ){ */
  /* brent->u.pos = brent->x.pos = brent->w.pos = brent->v.pos = brent->a = brent->b = 0.5 * (rmin + rmax); */
  brent->u.pos = brent->x.pos = brent->w.pos = brent->v.pos = brent->a = brent->b = 4.6e+2;
#endif//USE_FIXED_RMAX_VALUE

#ifdef  OUTPUT_AUTO_TUNING_PARAMETER
  __FPRINTF__(stdout, "min, max = %e, %e\n", rmin, rmax);
#endif//OUTPUT_AUTO_TUNING_PARAMETER
#else///DISABLE_AUTO_TUNING
  brent->b = CAST_R2D(thrust::reduce((thrust::device_ptr<real>)body_dev.neighbor, (thrust::device_ptr<real>)(body_dev.neighbor + Ni), ZERO, thrust::maximum<real>()));
#endif//DISABLE_AUTO_TUNING

#ifdef  HUNT_FIND_PARAMETER
  stopStopwatch(&(elapsed->sortNeighbors));
#endif//HUNT_FIND_PARAMETER


#ifdef  EXEC_BENCHMARK
#ifndef HUNT_FIND_PARAMETER
  stopStopwatch(&(elapsed->examineNeighbor_dev));
#else///HUNT_FIND_PARAMETER
  elapsed->examineNeighbor_dev += elapsed->searchNeighbor_kernel + elapsed->sortNeighbors;
#endif//HUNT_FIND_PARAMETER
#endif//EXEC_BENCHMARK
  __NOTE__("%s\n", "end");
}


/**
 * @fn updateParticleGroups
 *
 * @brief Update particle groups.
 *
 * @sa countContinuousNeighbor_kernel
 */
extern "C"
void updateParticleGroups
(const int Ni, laneinfo *laneInfo, const int inumPerLane, const int maxNgrp, int *Ngrp,
 const iparticle body_dev, int *inum_dev, int *inum_hst, const real rmax
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
 )
{
  __NOTE__("%s\n", "start");


#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

#ifdef  OUTPUT_AUTO_TUNING_PARAMETER
  __FPRINTF__(stdout, "brent_rad = %e\n", rmax);
#endif//OUTPUT_AUTO_TUNING_PARAMETER

  int Nrem = BLOCKSIZE(Ni, NTHREADS_SHRINK);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    countContinuousNeighbor_kernel<<<Nrem, NTHREADS_SHRINK>>>(Ni, body_dev.pos, inumPerLane, rmax * rmax, inum_dev);
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_SHRINK;
      countContinuousNeighbor_kernel<<<Nblck, NTHREADS_SHRINK>>>(Nsub, &body_dev.pos[hidx], inumPerLane, rmax * rmax, &inum_dev[hidx]);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("countContinuousNeighbor_kernel");

#ifdef  HUNT_FIND_PARAMETER
  stopStopwatch(&(elapsed->countNeighbors_kernel));
  initStopwatch();
#endif//HUNT_FIND_PARAMETER
  checkCudaErrors(cudaMemcpy(inum_hst, inum_dev, sizeof(int) * Ni, cudaMemcpyDeviceToHost));


  int irem = Ni;
  int hidx = 0;
  *Ngrp = 0;
  while( irem > 0 ){
    if( *Ngrp >= maxNgrp ){
      __KILL__(stderr, "ERROR: buffer overflow is predicted, %d particles remain\nIncrease NUM_IGROUP_SAFETY_FACTOR(%d) in tree/shrink_dev.h or decrease NEIGHBOR_LENGTH_SHRINK_FACTOR(%e) defined in tree/shrink_dev.h\n",
	       irem, NUM_IGROUP_SAFETY_FACTOR, NEIGHBOR_LENGTH_SHRINK_FACTOR);
    }/* if( *Ngrp >= maxNgrp ){ */

    const int num = IMIN(irem, inum_hst[hidx]);

    laneInfo[*Ngrp].head = hidx;
    laneInfo[*Ngrp].num  =  num;
    *Ngrp += 1;

    irem -= num;
    hidx += num;
  }/* while( irem > 0 ){ */

#ifndef NDEBUG
  fprintf(stdout, "Ni = %d, inumPerLane = %d, maxNgrp = %d, *Ngrp = %d\n", Ni, inumPerLane, maxNgrp, *Ngrp);
  fflush(stdout);
#endif//NDEBUG


#ifdef  HUNT_FIND_PARAMETER
  stopStopwatch(&(elapsed->commitNeighbors));
#endif//HUNT_FIND_PARAMETER
#ifdef  EXEC_BENCHMARK
#ifndef HUNT_FIND_PARAMETER
  stopStopwatch(&(elapsed->searchNeighbor_dev));
#else///HUNT_FIND_PARAMETER
  elapsed->searchNeighbor_dev += elapsed->countNeighbors_kernel + elapsed->commitNeighbors;
#endif//HUNT_FIND_PARAMETER
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}


/**
 * @fn commitParticleGroups
 *
 * @brief Commit particle groups.
 */
extern "C"
void commitParticleGroups(const int Ngrp, laneinfo *laneInfo_hst, laneinfo *laneInfo_dev)
{
  __NOTE__("%s\n", "start");


  /** copy lane information from host to device */
  const int totNgrp = NGROUPS * BLOCKSIZE(Ngrp, NGROUPS);
  for(int ii = Ngrp; ii < totNgrp; ii++)
    laneInfo_hst[ii] = nullInfo;

  checkCudaErrors(cudaMemcpy(laneInfo_dev, laneInfo_hst, totNgrp * sizeof(laneinfo), cudaMemcpyHostToDevice));

#ifndef NDEBUG
  fprintf(stdout, "Ngrp = %d, totNgrp = %d\n", Ngrp, totNgrp);
  fflush(stdout);
#endif//NDEBUG

#ifndef BLOCK_TIME_STEP
  cudaDeviceSynchronize();
#endif//BLOCK_TIME_STEP


  __NOTE__("%s\n", "end");
}
