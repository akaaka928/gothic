/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:58:09
 *                                                                       *
 *    Constructing octree structure for collisionless systems            *
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
#include <math.h>
#include <helper_cuda.h>
//-------------------------------------------------------------------------
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/find.h>
#include <thrust/sort.h>
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
//-------------------------------------------------------------------------
#include "macro.h"
#include "cudalib.h"
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#include "../misc/device.h"
//-------------------------------------------------------------------------
#ifdef  LOCALIZE_I_PARTICLES
#include "../sort/peano.h"
#endif//LOCALIZE_I_PARTICLES
//-------------------------------------------------------------------------
#include "make.h"
#include "walk_dev.h"
#include "stat_dev.h"
//-------------------------------------------------------------------------
static const laneinfo nullInfo = {NUM_BODY_MAX, 0};
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
void  freeParticleGroups(laneinfo  *laneInfo_hst, laneinfo  *laneInfo_dev, double  *laneTime_dev
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
			 , real  *dt_dev, real  *dt_hst, histogram_dt  *dtInfo
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  mycudaFreeHost(laneInfo_hst);
  mycudaFree    (laneInfo_dev);
  mycudaFree    (laneTime_dev);
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
  mycudaFree    (dt_dev);
  mycudaFreeHost(dt_hst);
  free(dtInfo);
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
muse allocParticleGroups
(laneinfo **laneInfo_hst, laneinfo **laneInfo_dev, double **laneTime_dev
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
 , real **dt_dev, real **dt_hst, histogram_dt **dtInfo, int *dtInfo_num
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
 , int *inumPerLane, int *maxNgrp, const int num_max, deviceProp devProp)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* number of i-particles per lane (a group of TSUB threads) */
  *inumPerLane = TSUB / NWARP;
  *maxNgrp = BLOCKSIZE(num_max, *inumPerLane) * NUM_IGROUP_SAFETY_FACTOR;
  //-----------------------------------------------------------------------
  /* memory allocation for laneInfo */
  mycudaMallocHost((void **)laneInfo_hst, (size_t)(*maxNgrp) * sizeof(laneinfo));  alloc.host   += (size_t)(*maxNgrp) * sizeof(laneinfo);
  mycudaMalloc    ((void **)laneInfo_dev, (size_t)(*maxNgrp) * sizeof(laneinfo));  alloc.device += (size_t)(*maxNgrp) * sizeof(laneinfo);
  mycudaMalloc    ((void **)laneTime_dev, (size_t)(*maxNgrp) * sizeof(  double));  alloc.device += (size_t)(*maxNgrp) * sizeof(  double);
  //-----------------------------------------------------------------------
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
  *dtInfo_num = MAXIMUM_PHKEY_LEVEL;
  *dtInfo = (histogram_dt *)malloc((*dtInfo_num) * sizeof(histogram_dt));
  alloc.host += (*dtInfo_num) * sizeof(histogram_dt);
  if( *dtInfo == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dtInfo\n");  }
  /* the size of the array is set to be a multiple of NTHREADS */
  size_t size = (size_t)num_max;
  if( (num_max % NTHREADS) != 0 )
    size += (size_t)(NTHREADS - (num_max % NTHREADS));
  mycudaMalloc    ((void **)dt_dev, size * sizeof(real));  alloc.device += size * sizeof(real);
  mycudaMallocHost((void **)dt_hst, size * sizeof(real));  alloc.host   += size * sizeof(real);
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
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
  }
  /* initLaneTime_kernel<<<BLOCKSIZE(*maxNgrp, 1024), 1024>>>(*maxNgrp, *laneTime_dev); */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
//-------------------------------------------------------------------------
struct greater_than
{
  real val;

greater_than(real val) : val(val) {}

  __host__ __device__
  int operator()(const real &x) const {
    return (x > val);
  }
};
//-------------------------------------------------------------------------
static inline int bisection(const real val, const int num, histogram_dt tab[])
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = num - 1;
  //-----------------------------------------------------------------------
  /* prohibit extraporation */
  if( val >= tab[rr].dt )    return (rr);
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)(ll + rr)) >> 1;
    //---------------------------------------------------------------------
    if( val <= tab[cc].dt )      rr = (int)cc;
    else                         ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr )
      return (ll);
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__global__ void cpDt_kernel
(const int Ni, READ_ONLY velocity * RESTRICT vel, real * RESTRICT dt)
{
  //-----------------------------------------------------------------------
  const int ii = GLOBALIDX_X1D;
  if( ii < Ni )
    dt[ii] = vel[ii].dt;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
//-------------------------------------------------------------------------
extern "C"
void updateParticleGroups
(const int Ni, laneinfo *laneInfo, const int inumPerLane, const int maxNgrp, int *Ngrp
#ifdef  LOCALIZE_I_PARTICLES
 , PHint *peano
#endif//LOCALIZE_I_PARTICLES
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
 , const iparticle body, real *dt_dev, real *dt_hst, int *gnum, histogram_dt **dtInfo
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* determine VARIABLE_NEIGHBOR_LEVEL */
  //-----------------------------------------------------------------------
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
  //-----------------------------------------------------------------------
  /* sort time step of N-body particles */
  //-----------------------------------------------------------------------
  /* body.velocity.dt --> dt_dev */
  int Nrem = BLOCKSIZE(Ni, 1024);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    cpDt_kernel<<<Nrem, 1024>>>(Ni, body.vel, dt_dev);
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
      cpDt_kernel<<<Nblck, 1024>>>(Nsub, &body.vel[hidx], &dt_dev[hidx]);
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }
  /* cpDt_kernel<<<BLOCKSIZE(Ni, 1024), 1024>>>(Ni, body.vel, dt_dev); */
  getLastCudaError("cpDt_kernel");
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaMemcpy(dt_hst, dt_dev, sizeof(real) * Ni, cudaMemcpyDeviceToHost));
  thrust::stable_sort((thrust::device_ptr<real>)dt_dev, (thrust::device_ptr<real>)(dt_dev + Ni));
  //-----------------------------------------------------------------------
  /* get histogram and mode about time step of N-body particles */
  int pidx = 0;
  int grem = *gnum;
  *gnum = 0;
  int nmax = -1;
  /* int mode = -1; */
  while( true ){
    //---------------------------------------------------------------------
    if( grem <= 0 ){
      const int nadd = 4;
      *dtInfo = (histogram_dt *)realloc(*dtInfo, sizeof(histogram_dt) * (nadd + (*gnum)));
      if( *dtInfo == NULL ){	__KILL__(stderr, "ERROR: failure to rescale dtInfo\n");      }
      grem = nadd;
    }
    //---------------------------------------------------------------------
    real dt_tmp;
    checkCudaErrors(cudaMemcpy(&dt_tmp, &dt_dev[pidx], sizeof(real), cudaMemcpyDeviceToHost));
    //---------------------------------------------------------------------
    /* count up # of N-body particles that share the time step */
    thrust::device_vector<real>::iterator iter1 =                 (thrust::device_ptr<real>)(dt_dev + pidx);
    thrust::device_vector<real>::iterator iter2 = thrust::find_if((thrust::device_ptr<real>)(dt_dev + pidx), (thrust::device_ptr<real>)(dt_dev + Ni), greater_than(dt_tmp));
    const int num = thrust::distance(iter1, iter2);
    //---------------------------------------------------------------------
    if( num > nmax ){
      nmax = num;
      /* mode = *gnum; */
    }
    //---------------------------------------------------------------------
    (*dtInfo)[*gnum].dt  = dt_tmp;
    (*dtInfo)[*gnum].num = num;
    grem--;    *gnum += 1;
    //---------------------------------------------------------------------
    pidx += num;
    if( pidx >= Ni )
      break;
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
  /* determine keyJump for each time step */
#if 1
  float Ntot = 0.0f;
  const float log2Ntot = log2f((float)nmax);
#endif
#if 1
  ulong costMax = 0;
  for(int ii = 0; ii < *gnum; ii++){
    ulong cost = (ulong)(*dtInfo)[ii].num * (ulong)1 << ((*gnum) - 1 - ii);
    if( cost > costMax )
      costMax = cost;
    (*dtInfo)[ii].jump = cost;
  }/* for(int ii = 0; ii < *gnum; ii++){ */
  costMax /= 1000;
#endif
  for(int ii = 0; ii < *gnum; ii++){
    //---------------------------------------------------------------------
#if 1
    Ntot += (float)(*dtInfo)[ii].num;
#if 1
    int splitLevel = MINIMUM_NEIGHBOR_PHKEY_LEVEL + (int)floorf(ONE_THIRD * (log2Ntot - log2f(Ntot) + (float)(2 * ((*gnum) - 1 - ii))));
#else
    int splitLevel = MINIMUM_NEIGHBOR_PHKEY_LEVEL + (int)floorf(ONE_THIRD * (log2Ntot - log2f(Ntot)));
#endif
#else
    /* a base to change PH-key split level set to be 10 (this value can be arbitrary) */
    int splitLevel = (int)floorf(log10f((float)nmax / (float)(*dtInfo)[ii].num));
#endif

#if 1
    const int tmpLevel = (int)ceilf(ONE_THIRD * (float)(2 * ((*gnum) - 1 - ii)));
    if( tmpLevel > splitLevel )
      splitLevel = tmpLevel;
#endif


    splitLevel = (int)ceilf((float)(2 * ii) * ONE_THIRD);


    if( splitLevel > MAXIMUM_NEIGHBOR_PHKEY_LEVEL )
      splitLevel = MAXIMUM_NEIGHBOR_PHKEY_LEVEL;
#if 1
    if( (*dtInfo)[ii].jump >= costMax )
      splitLevel = 0;
#endif
    //---------------------------------------------------------------------
    (*dtInfo)[ii].jump = (PHint)1 << ((MAXIMUM_PHKEY_LEVEL - splitLevel) * 3);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < *gnum; ii++){ */
  //-----------------------------------------------------------------------
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* for(int levelIdx = (NUM_PHKEY_LEVEL - 1); levelIdx >= 0; levelIdx--) */
  /*   if( level[levelIdx].num != 0 ){ */
  /*     *bottomLev = levelIdx; */
  /*     break; */
  /*   } */
#   if  defined(LOCALIZE_I_PARTICLES) && !defined(USE_VARIABLE_NEIGHBOR_LEVEL)
  const int splitLevel =  NEIGHBOR_PHKEY_LEVEL;
  const PHint keyJump = (PHint)1 << ((MAXIMUM_PHKEY_LEVEL - splitLevel) * 3);
#endif//defined(LOCALIZE_I_PARTICLES) && !defined(USE_VARIABLE_NEIGHBOR_LEVEL)
  //-----------------------------------------------------------------------
  int irem = Ni;
  int hidx = 0;
  *Ngrp = 0;
  //-----------------------------------------------------------------------
  while( irem > 0 ){
    //---------------------------------------------------------------------
    if( *Ngrp >= maxNgrp ){
#ifdef  LOCALIZE_I_PARTICLES
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
      __KILL__(stderr, "ERROR: buffer overflow is predicted, %d particles remain\nIncrease NUM_IGROUP_SAFETY_FACTOR(%d) in tree/stat_dev.h or decrease MAXIMUM_NEIGHBOR_PHKEY_LEVEL(%d) defined in tree/stat_dev.h\n",
	       irem, NUM_IGROUP_SAFETY_FACTOR, MAXIMUM_NEIGHBOR_PHKEY_LEVEL);
#else///USE_VARIABLE_NEIGHBOR_LEVEL
      __KILL__(stderr, "ERROR: buffer overflow is predicted, %d particles remain (splitLev = %d, allowed keyJump is %zu)\nIncrease NUM_IGROUP_SAFETY_FACTOR(%d) in tree/stat_dev.h or decrease NEIGHBOR_PHKEY_LEVEL(%d) defined in tree/stat_dev.h or Makefile\n",
	       irem, splitLevel, keyJump, NUM_IGROUP_SAFETY_FACTOR, NEIGHBOR_PHKEY_LEVEL);
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
#else///LOCALIZE_I_PARTICLES
      __KILL__(stderr, "ERROR: buffer overflow is predicted, %d particles remain\nIncrease NUM_IGROUP_SAFETY_FACTOR(%d) in tree/stat_dev.h\n",
	       irem, NUM_IGROUP_SAFETY_FACTOR);
#endif//LOCALIZE_I_PARTICLES
    }
    //---------------------------------------------------------------------
#ifndef LOCALIZE_I_PARTICLES
    //---------------------------------------------------------------------
    const int  num = IMIN(irem, inumPerLane);
    //---------------------------------------------------------------------
#else///LOCALIZE_I_PARTICLES
    //---------------------------------------------------------------------
    const int ntmp = IMIN(irem, inumPerLane);
    const PHint keyHead = peano[hidx];
    int num = ntmp;
    //---------------------------------------------------------------------
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
    /* set keyJump based on the particle time step */
    real dt_min = dt_hst[hidx];
    int gidx = bisection(dt_min, *gnum, *dtInfo);
    PHint keyJump = (*dtInfo)[gidx].jump;
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
    //---------------------------------------------------------------------
    for(int ii = 1; ii < ntmp; ii++){
      //-------------------------------------------------------------------
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
      /* reset keyJump based on the particle time step, if required */
      const real dt_tmp = dt_hst[hidx + ii];
      if( dt_tmp < dt_min ){
	dt_min = dt_tmp;
	gidx = bisection(dt_min, gidx, *dtInfo);
	keyJump = (*dtInfo)[gidx].jump;
      }
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
      //-------------------------------------------------------------------
      /* check PH-key difference */
      if( peano[hidx + ii] > (keyHead + keyJump) ){
	num = ii;
	break;
      }/* if( peano[hidx + ii] > (keyHead + keyJump) ){ */
      //-------------------------------------------------------------------
    }/* for(int ii = 1; ii < ntmp; ii++){ */
    //---------------------------------------------------------------------
#endif//LOCALIZE_I_PARTICLES
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
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
  *gnum += grem;
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
  //-----------------------------------------------------------------------


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

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
