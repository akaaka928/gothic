/*************************************************************************\
 *                                                                       *
                  last updated on 2016/01/20(Wed) 14:17:09
 *                                                                       *
 *    Constructing octree structure for collisionless systems            *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
/* #define PRINT_OUT_NEIGHBOR_LENGTH */
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <helper_cuda.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#include "../misc/device.h"
//-------------------------------------------------------------------------
#include "make.h"
#include "walk_dev.h"
#include "make_dev.h"
//-------------------------------------------------------------------------
#ifdef  LOCALIZE_I_PARTICLES
#       ifdef  USE_BRENT_METHOD
#include "../misc/brent.h"
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#       else///USE_BRENT_METHOD
#ifdef  CUB_AVAILABLE
#include <cub/device/device_radix_sort.cuh>
#else///CUB_AVAILABLE
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
# endif//CUB_AVAILABLE
#       endif//USE_BRENT_METHOD
#include "neighbor_dev.h"
#endif//LOCALIZE_I_PARTICLES
//-------------------------------------------------------------------------
#include "shrink_dev.h"
//-------------------------------------------------------------------------
static const laneinfo nullInfo = {NUM_BODY_MAX, 0};
//-------------------------------------------------------------------------
#define NEIGHBOR_NUM_LANE (TSUB / NWARP)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__global__ void initLaneTime_kernel(const int num, double *laneTime)
{
  //-----------------------------------------------------------------------
  const int idx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( idx < num )
    laneTime[idx] = DBL_MAX;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
extern "C"
void freeParticleGroups
(laneinfo  *laneInfo_hst, laneinfo  *laneInfo_dev, double  *laneTime_dev
#ifdef  LOCALIZE_I_PARTICLES
 , int  *inum_hst, int  *inum_dev
#   if  defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
 , void  *temp_storage, real  *outCub
#endif//defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
#endif//LOCALIZE_I_PARTICLES
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  mycudaFreeHost(laneInfo_hst);
  mycudaFree    (laneInfo_dev);
  mycudaFree    (laneTime_dev);
#ifdef  LOCALIZE_I_PARTICLES
  mycudaFreeHost(    inum_hst);
  mycudaFree    (    inum_dev);
#   if  defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
  mycudaFree    (temp_storage);
  mycudaFree    (      outCub);
#endif//defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
#endif//LOCALIZE_I_PARTICLES
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
muse allocParticleGroups
(laneinfo **laneInfo_hst, laneinfo **laneInfo_dev, double **laneTime_dev
#ifdef  LOCALIZE_I_PARTICLES
 , int **inum_hst, int **inum_dev
#   if  defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
 , soaCUBreal *util, void **temp_storage, real **outCub, iparticle body_dev
#endif//defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
#endif//LOCALIZE_I_PARTICLES
 , int *inumPerLane, int *maxNgrp, const int num_max, deviceProp devProp)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* number of i-particles per lane (a group of TSUB threads) */
  *inumPerLane = NEIGHBOR_NUM_LANE;
  *maxNgrp = BLOCKSIZE(num_max, *inumPerLane) * NUM_IGROUP_SAFETY_FACTOR;
  //-----------------------------------------------------------------------
  /* memory allocation for laneInfo */
  mycudaMallocHost((void **)laneInfo_hst, (size_t)(*maxNgrp) * sizeof(laneinfo));  alloc.host   += (size_t)(*maxNgrp) * sizeof(laneinfo);
  mycudaMalloc    ((void **)laneInfo_dev, (size_t)(*maxNgrp) * sizeof(laneinfo));  alloc.device += (size_t)(*maxNgrp) * sizeof(laneinfo);
  mycudaMalloc    ((void **)laneTime_dev, (size_t)(*maxNgrp) * sizeof(  double));  alloc.device += (size_t)(*maxNgrp) * sizeof(  double);
  //-----------------------------------------------------------------------
  /* memory allocation to estimate number of i-particle within a group */
#ifdef  LOCALIZE_I_PARTICLES
  mycudaMallocHost((void **)inum_hst, (size_t)num_max * sizeof(int));  alloc.host   += (size_t)num_max * sizeof(int);
  mycudaMalloc    ((void **)inum_dev, (size_t)num_max * sizeof(int));  alloc.device += (size_t)num_max * sizeof(int);
#endif//LOCALIZE_I_PARTICLES
  //-----------------------------------------------------------------------
  /* memory allocation to use CUB */
#   if  defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
  mycudaMalloc((void **)outCub, (size_t)num_max * sizeof(real));  alloc.device += (size_t)num_max * sizeof(real);
  size_t temp_storage_size = 0;
  *temp_storage = NULL;
  cub::DeviceRadixSort::SortKeys(*temp_storage, temp_storage_size, body_dev.neighbor, *outCub, num_max);
  mycudaMalloc(temp_storage, temp_storage_size);  alloc.device += temp_storage_size;
  util->temp_storage = *temp_storage;
  util->temp_storage_size = temp_storage_size;
  util->out = *outCub;
#endif//defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
  //-----------------------------------------------------------------------
  /* initialize laneInfo */
  for(int ii = 0; ii < *maxNgrp; ii++)
    (*laneInfo_hst)[ii] = nullInfo;
  checkCudaErrors(cudaMemcpy(*laneInfo_dev, *laneInfo_hst, (*maxNgrp) * sizeof(laneinfo), cudaMemcpyHostToDevice));
  //-----------------------------------------------------------------------
  /* initialize laneTime */
  int Nrem = BLOCKSIZE(*maxNgrp, 1024);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    initLaneTime_kernel<<<Nrem, 1024>>>(*maxNgrp, *laneTime_dev);
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
      initLaneTime_kernel<<<Nblck, 1024>>>(Nsub, &((*laneTime_dev)[hidx]));
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  /* initLaneTime_kernel<<<BLOCKSIZE(*maxNgrp, 1024), 1024>>>(*maxNgrp, *laneTime_dev); */
  //-----------------------------------------------------------------------
  /* set SM/L1 configuraetion for neighbor searching */
#   if  defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  setGlobalConstants_neighbor_dev_cu();
#endif//defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#if 1
//-------------------------------------------------------------------------
__global__ void countContinuousNeighbor_kernel(const int Ni, position * RESTRICT ibody, const int inumPerLane, const real r2max, int * RESTRICT inum)
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  __shared__ position ipos[NTHREADS_SHRINK + NEIGHBOR_NUM_LANE];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  if( gidx < Ni ){
    //---------------------------------------------------------------------
    /* load position of i-particles and neighbor candidates */
    //---------------------------------------------------------------------
    /* int idx = gidx - NEIGHBOR_NUM_LANE;    if( idx < 0 )      idx = 0; */
    /* ipos[tidx] = ibody[idx]; */
    ipos[tidx] = ibody[gidx];
    if( tidx < NEIGHBOR_NUM_LANE ){
      //-------------------------------------------------------------------
      /* idx = gidx - NEIGHBOR_NUM_LANE + NTHREADS_SHRINK; */
      /* ipos[NTHREADS_SHRINK + tidx] = ibody[idx]; */
      int idx = gidx + NTHREADS_SHRINK;      if( idx > (Ni - 1) )	idx = Ni - 1;
      ipos[NTHREADS_SHRINK + tidx] = ibody[idx];
      //-------------------------------------------------------------------
    }/* if( tidx < NEIGHBOR_NUM_LANE ){ */
    __syncthreads();
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* calculate distance with NEIGHBOR_NUM_LANE particles and remember the maximum */
    //---------------------------------------------------------------------
    const position pi = ipos[tidx];
    int nmax = inumPerLane;
    //---------------------------------------------------------------------
#pragma unroll
    for(int ii = 0; ii < NEIGHBOR_NUM_LANE; ii++){
      //-------------------------------------------------------------------
      const position pj = ipos[tidx + 1 + ii];
      //-------------------------------------------------------------------
      const real dx = pj.x - pi.x;
      const real dy = pj.y - pi.y;
      const real dz = pj.z - pi.z;
      const real r2 = 1.0e-30f + dx * dx + dy * dy + dz * dz;
      //-------------------------------------------------------------------
      if( r2 > r2max ){
	nmax = 1 + ii;
	break;
      }
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NEIGHBOR_NUM_LANE; ii++){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* commit the derived maximum number as neighbor particles within rmax */
    //---------------------------------------------------------------------
    inum[gidx] = nmax;
    //---------------------------------------------------------------------
  }/* if( gidx < Ni ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else
//-------------------------------------------------------------------------
__global__ void countContinuousNeighbor_kernel(const int Ni, position * RESTRICT ipos, const int inumPerLane, const real r2max, int * RESTRICT inum)
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( gidx < Ni ){
    //---------------------------------------------------------------------
    const position pi = ipos[gidx];
    //---------------------------------------------------------------------
    const int num = (inumPerLane < (Ni - gidx)) ? inumPerLane : (Ni - gidx);
    int nmax = num;
    for(int jj = 1; jj < num; jj++){
      //-------------------------------------------------------------------
      const position pj = ipos[gidx + jj];
      const real dx = pj.x - pi.x;
      const real dy = pj.y - pi.y;
      const real dz = pj.z - pi.z;
      const real r2 = dx * dx + dy * dy + dz * dz;
      //-------------------------------------------------------------------
      if( r2 > r2max ){
	nmax = jj;
	break;
      }
      //-------------------------------------------------------------------
    }/* for(int jj = 1; jj < num; jj++){ */
    //---------------------------------------------------------------------
    inum[gidx] = nmax;
    //---------------------------------------------------------------------
  }/* if( tidx < Ni ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  PRINT_OUT_NEIGHBOR_LENGTH
__global__ void storeRadius(const int Ni, READ_ONLY position * RESTRICT pi, real * RESTRICT rad)
{
  //-----------------------------------------------------------------------
  const int idx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( idx < Ni ){
    //---------------------------------------------------------------------
    const position ipos = pi[idx];
    const real r2 = ipos.x * ipos.x + ipos.y * ipos.y + ipos.z * ipos.z;
    const real rr = r2 * RSQRT(r2);
    //---------------------------------------------------------------------
    rad[idx] = rr;
    //---------------------------------------------------------------------
  }/* if( idx < Ni ){ */
  //-----------------------------------------------------------------------
}
#endif//PRINT_OUT_NEIGHBOR_LENGTH
//-------------------------------------------------------------------------
extern "C"
void examineParticleSeparation
(const int Ni, iparticle body_dev
#ifdef  USE_BRENT_METHOD
 , brentStatus *brent
#else///USE_BRENT_METHOD
#ifdef  CUB_AVAILABLE
 , soaCUBreal util
#endif//CUB_AVAILABLE
 , real *rmax
#endif//USE_BRENT_METHOD
#ifndef FACILE_NEIGHBOR_SEARCH
 , const soaTreeCell cell, const soaTreeNode node, const soaMakeTreeBuf makeBuf, const soaNeighborSearchBuf searchBuf, deviceProp devProp
#endif//FACILE_NEIGHBOR_SEARCH
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set a criterion to judge whether unify i-particles or not */
  //-----------------------------------------------------------------------
#ifdef  FACILE_NEIGHBOR_SEARCH
  facileNeighborSearching_dev(Ni, body_dev);
#else///FACILE_NEIGHBOR_SEARCH
  searchNeighbors_dev(Ni, body_dev, cell, node, makeBuf, searchBuf, devProp);
#endif//FACILE_NEIGHBOR_SEARCH
#ifdef  HUNT_FIND_PARAMETER
  stopStopwatch(&(elapsed->searchNeighbor_kernel));
  initStopwatch();
#endif//HUNT_FIND_PARAMETER
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifndef PRINT_OUT_NEIGHBOR_LENGTH
  //-----------------------------------------------------------------------
#ifndef USE_BRENT_METHOD
#ifdef  CUB_AVAILABLE
  cub::DeviceRadixSort::SortKeys(util.temp_storage, util.temp_storage_size, body_dev.neighbor, util.out, Ni);
#else///CUB_AVAILABLE
  thrust::sort((thrust::device_ptr<real>)body_dev.neighbor, (thrust::device_ptr<real>)(body_dev.neighbor + Ni));
#endif//CUB_AVAILABLE
#endif//USE_BRENT_METHOD
  //-----------------------------------------------------------------------
#else///PRINT_OUT_NEIGHBOR_LENGTH
  //-----------------------------------------------------------------------
  real *neighbor_hst;
  mycudaMallocHost((void **)&neighbor_hst, sizeof(real) * Ni);
  checkCudaErrors(cudaMemcpy(neighbor_hst, body_dev.neighbor, sizeof(real) * Ni, cudaMemcpyDeviceToHost));
  real *rad_dev;  mycudaMalloc    ((void **)&rad_dev, sizeof(real) * Ni);
  real *rad_hst;  mycudaMallocHost((void **)&rad_hst, sizeof(real) * Ni);
  storeRadius<<<BLOCKSIZE(Ni, 1024), 1024>>>(Ni, body_dev.pos, rad_dev);
  checkCudaErrors(cudaMemcpy(rad_hst, rad_dev, sizeof(real) * Ni, cudaMemcpyDeviceToHost));
  for(int ii = 0; ii < Ni; ii++)
    fprintf(stderr, "%d\t%e\t%e\n", ii, rad_hst[ii], neighbor_hst[ii]);
  mycudaFreeHost(neighbor_hst);
  mycudaFree    (rad_dev);
  mycudaFreeHost(rad_hst);
  exit(0);
  //-----------------------------------------------------------------------
#endif//PRINT_OUT_NEIGHBOR_LENGTH
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  USE_BRENT_METHOD
  //-----------------------------------------------------------------------
  const real rmax = thrust::reduce((thrust::device_ptr<real>)body_dev.neighbor, (thrust::device_ptr<real>)(body_dev.neighbor + Ni), ZERO, thrust::maximum<real>());
  const real rmin = rmax * NEIGHBOR_LENGTH_SHRINK_FACTOR;
  //-----------------------------------------------------------------------
/*   real rmin, rmax; */
/* #ifdef  CUB_AVAILABLE */
/*   checkCudaErrors(cudaMemcpy(&rmax, &util.out[Ni - 1], sizeof(real), cudaMemcpyDeviceToHost)); */
/* #else///CUB_AVAILABLE */
/*   checkCudaErrors(cudaMemcpy(&rmax, &body_dev.neighbor[Ni - 1], sizeof(real), cudaMemcpyDeviceToHost)); */
/* #endif//CUB_AVAILABLE */
/*   rmin = rmax * NEIGHBOR_LENGTH_SHRINK_FACTOR; */
  static bool initialized = false;
  if( initialized )
    brentPerturb(brent, (double)rmin, (double)rmax);
  else{
    initialized = true;
    brentInit1st(brent, (double)rmin, (double)rmax);
  }/* else{ */
  //-----------------------------------------------------------------------
#else///USE_BRENT_METHOD
  //-----------------------------------------------------------------------
#ifdef  CUB_AVAILABLE
  checkCudaErrors(cudaMemcpy(rmax, &util.out[(int)ceilf((float)Ni * (UNITY - NEIGHBOR_LENGTH_UPPER_FRACTION)) - 1], sizeof(real), cudaMemcpyDeviceToHost));
#else///CUB_AVAILABLE
  checkCudaErrors(cudaMemcpy(rmax, &body_dev.neighbor[(int)ceilf((float)Ni * (UNITY - NEIGHBOR_LENGTH_UPPER_FRACTION)) - 1], sizeof(real), cudaMemcpyDeviceToHost));
#endif//CUB_AVAILABLE
  *rmax *= (UNITY + NEIGHBOR_LENGTH_EXTEND_FRACTION);
#if 0
  fprintf(stdout, "target = %d: rmax = %10.8e\n", (int)ceilf((float)Ni * (UNITY - NEIGHBOR_LENGTH_UPPER_FRACTION)) - 1, *rmax);
  fflush(stdout);
  exit(0);
#endif
  //-----------------------------------------------------------------------
#endif//USE_BRENT_METHOD
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  HUNT_FIND_PARAMETER
  stopStopwatch(&(elapsed->sortNeighbors));
#endif//HUNT_FIND_PARAMETER
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
#ifndef HUNT_FIND_PARAMETER
  stopStopwatch(&(elapsed->examineNeighbor_dev));
#else///HUNT_FIND_PARAMETER
  elapsed->examineNeighbor_dev += elapsed->searchNeighbor_kernel + elapsed->sortNeighbors;
#endif//HUNT_FIND_PARAMETER
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
void updateParticleGroups
(const int Ni, laneinfo *laneInfo, const int inumPerLane, const int maxNgrp, int *Ngrp
#ifdef  LOCALIZE_I_PARTICLES
 , const iparticle body_dev, int *inum_dev, int *inum_hst
#ifdef  USE_BRENT_METHOD
 , const real rmax
#endif//USE_BRENT_METHOD
#   if  !defined(FACILE_NEIGHBOR_SEARCH) && !defined(USE_BRENT_METHOD)
 , const soaTreeCell cell, const soaTreeNode node, const soaMakeTreeBuf makeBuf, const soaNeighborSearchBuf searchBuf, deviceProp devProp
#endif//!defined(FACILE_NEIGHBOR_SEARCH) && !defined(USE_BRENT_METHOD)
#endif//LOCALIZE_I_PARTICLES
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#   if  defined(LOCALIZE_I_PARTICLES) && !defined(USE_BRENT_METHOD)
  real rmax;
  examineParticleSeparation
    (Ni, body_dev, &rmax
#ifndef FACILE_NEIGHBOR_SEARCH
     , cell, node, makeBuf, searchBuf, devProp
#endif//FACILE_NEIGHBOR_SEARCH
#ifdef  EXEC_BENCHMARK
     , elapsed
#endif//EXEC_BENCHMARK
     );
#endif//defined(LOCALIZE_I_PARTICLES) && !defined(USE_BRENT_METHOD)
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  LOCALIZE_I_PARTICLES
  //-----------------------------------------------------------------------
  int Nrem = BLOCKSIZE(Ni, NTHREADS_SHRINK);
  //-----------------------------------------------------------------------
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    countContinuousNeighbor_kernel<<<Nrem, NTHREADS_SHRINK>>>(Ni, body_dev.pos, inumPerLane, rmax * rmax, inum_dev);
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
      int Nsub = Nblck * NTHREADS_SHRINK;
      countContinuousNeighbor_kernel<<<Nblck, NTHREADS_SHRINK>>>(Nsub, &body_dev.pos[hidx], inumPerLane, rmax * rmax, &inum_dev[hidx]);
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
  getLastCudaError("countContinuousNeighbor_kernel");
  //-----------------------------------------------------------------------
  /* countContinuousNeighbor_kernel<<<BLOCKSIZE(Ni, NTHREADS_SHRINK), NTHREADS_SHRINK>>>(Ni, body_dev.pos, inumPerLane, rmax * rmax, inum_dev); */
  /* getLastCudaError("countContinuousNeighbor_kernel"); */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  HUNT_FIND_PARAMETER
  stopStopwatch(&(elapsed->countNeighbors_kernel));
  initStopwatch();
#endif//HUNT_FIND_PARAMETER
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaMemcpy(inum_hst, inum_dev, sizeof(int) * Ni, cudaMemcpyDeviceToHost));
  //-----------------------------------------------------------------------
#endif//LOCALIZE_I_PARTICLES
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  int irem = Ni;
  int hidx = 0;
  *Ngrp = 0;
  //-----------------------------------------------------------------------
  while( irem > 0 ){
    //---------------------------------------------------------------------
    if( *Ngrp >= maxNgrp ){
#ifdef  LOCALIZE_I_PARTICLES
#       ifdef  USE_BRENT_METHOD
      __KILL__(stderr, "ERROR: buffer overflow is predicted, %d particles remain\nIncrease NUM_IGROUP_SAFETY_FACTOR(%d) in tree/shrink_dev.h or decrease NEIGHBOR_LENGTH_SHRINK_FACTOR(%e) defined in tree/shrink_dev.h\n",
	       irem, NUM_IGROUP_SAFETY_FACTOR, NEIGHBOR_LENGTH_SHRINK_FACTOR);
#       else///USE_BRENT_METHOD
      __KILL__(stderr, "ERROR: buffer overflow is predicted, %d particles remain\nIncrease NUM_IGROUP_SAFETY_FACTOR(%d) in tree/shrink_dev.h or decrease NEIGHBOR_LENGTH_EXTEND_FRACTION(%e) or NEIGHBOR_LENGTH_UPPER_FRACTION(%e) defined in tree/shrink_dev.h\n",
	       irem, NUM_IGROUP_SAFETY_FACTOR, NEIGHBOR_LENGTH_EXTEND_FRACTION, NEIGHBOR_LENGTH_UPPER_FRACTION);
#       endif//USE_BRENT_METHOD
#else///LOCALIZE_I_PARTICLES
      __KILL__(stderr, "ERROR: buffer overflow is predicted, %d particles remain\nIncrease NUM_IGROUP_SAFETY_FACTOR(%d) in tree/shrink_dev.h\n",
	       irem, NUM_IGROUP_SAFETY_FACTOR);
#endif//LOCALIZE_I_PARTICLES
    }
    //---------------------------------------------------------------------
    const int num = IMIN(irem,
#ifdef  LOCALIZE_I_PARTICLES
			 inum_hst[hidx]
#else///LOCALIZE_I_PARTICLES
			 inumPerLane
#endif//LOCALIZE_I_PARTICLES
			 );
    //---------------------------------------------------------------------
    laneInfo[*Ngrp].head = hidx;
    laneInfo[*Ngrp].num  =  num;
    *Ngrp += 1;
    //---------------------------------------------------------------------
    irem -= num;
    hidx += num;
    //---------------------------------------------------------------------
  }/* while( irem > 0 ){ */
  //-----------------------------------------------------------------------
#ifndef NDEBUG
  fprintf(stdout, "Ni = %d, inumPerLane = %d, maxNgrp = %d, *Ngrp = %d\n", Ni, inumPerLane, maxNgrp, *Ngrp);
  fflush(stdout);
#endif//NDEBUG
  //-----------------------------------------------------------------------
#if 0
  printf("%s(%d): %s: tentative kill for performance tuning\n", __FILE__, __LINE__, __func__);
  exit(0);
#endif
  //-----------------------------------------------------------------------
#if 0
  printf("# rmax = %e, *Ngrp = %d\n", rmax, *Ngrp);
  fflush(stdout);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
void commitParticleGroups(const int Ngrp, laneinfo *laneInfo_hst, laneinfo *laneInfo_dev)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* copy lane information from host to device */
  //-----------------------------------------------------------------------
  const int totNgrp = NGROUPS * BLOCKSIZE(Ngrp, NGROUPS);
  for(int ii = Ngrp; ii < totNgrp; ii++)
    laneInfo_hst[ii] = nullInfo;
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaMemcpy(laneInfo_dev, laneInfo_hst, totNgrp * sizeof(laneinfo), cudaMemcpyHostToDevice));
  //-----------------------------------------------------------------------
#ifndef NDEBUG
  fprintf(stdout, "Ngrp = %d, totNgrp = %d\n", Ngrp, totNgrp);
  /* for(int ii = 0; ii < totNgrp; ii++) */
  /*   printf("%d\t%d\t%d\n", ii, laneInfo_hst[ii].head, laneInfo_hst[ii].num); */
  fflush(stdout);
#endif//NDEBUG
  //-----------------------------------------------------------------------
#ifndef BLOCK_TIME_STEP
  cudaDeviceSynchronize();
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < totNgrp; ii++)
    fprintf(stderr, "%d\t%d\t%d\n", ii, laneInfo_hst[ii].head, laneInfo_hst[ii].num);
  fflush(NULL);
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
