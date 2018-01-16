/**
 * @file potential_dev.cu
 *
 * @brief Source code for calculating gravity from external fixed-potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/16 (Tue)
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
#include <helper_cuda.h>

#include "macro.h"
#include "cudalib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "make.h"
#include "walk_dev.h"
#include "potential_dev.h"


/**
 * @fn allocSphericalPotentialTable_dev
 *
 * @brief Memory allocation on the accelerator device for external fixed-potential field by spherical components.
 */
extern "C"
muse allocSphericalPotentialTable_dev(real **rad, pot2 **Phi, const int Nr)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};

  mycudaMalloc((void **)rad, Nr * sizeof(real));  alloc.device += Nr * sizeof(real);
  mycudaMalloc((void **)Phi, Nr * sizeof(pot2));  alloc.device += Nr * sizeof(pot2);

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn  freeSphericalPotentialTable_dev
 *
 * @brief Memory deallocation on the accelerator device for external fixed-potential field by spherical components.
 */
extern "C"
void  freeSphericalPotentialTable_dev(real  *rad, pot2  *Phi)
{
  __NOTE__("%s\n", "start");

  mycudaFree(rad);
  mycudaFree(Phi);

  __NOTE__("%s\n", "end");
}


/**
 * @fn allocSphericalPotentialTable_hst
 *
 * @brief Memory allocation on the host for external fixed-potential field by spherical components.
 */
extern "C"
muse allocSphericalPotentialTable_hst(real **rad, pot2 **Phi, const int Nr)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};

  mycudaMallocHost((void **)rad, Nr * sizeof(real));  alloc.host += Nr * sizeof(real);
  mycudaMallocHost((void **)Phi, Nr * sizeof(pot2));  alloc.host += Nr * sizeof(pot2);

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn  freeSphericalPotentialTable_hst
 *
 * @brief Memory deallocation on the host for external fixed-potential field by spherical components.
 */
extern "C"
void  freeSphericalPotentialTable_hst(real  *rad, pot2  *Phi)
{
  __NOTE__("%s\n", "start");

  mycudaFreeHost(rad);
  mycudaFreeHost(Phi);

  __NOTE__("%s\n", "end");
}


/**
 * @fn setSphericalPotentialTable_dev
 *
 * @brief Set external fixed-potential field by spherical components on the device.
 */
extern "C"
void setSphericalPotentialTable_dev(real *rad_hst, pot2 *Phi_hst, real *rad_dev, pot2 *Phi_dev, const int Nr)
{
  __NOTE__("%s\n", "start");

  checkCudaErrors(cudaMemcpy(rad_dev, rad_hst, sizeof(real) * Nr, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(Phi_dev, Phi_hst, sizeof(pot2) * Nr, cudaMemcpyHostToDevice));

  __NOTE__("%s\n", "end");
}














/**
 * @fn calcExternalForce_spherical
 *
 * @brief Get F_r = -\nabla \Phi(r) and \Phi(r) by an external potential field at the specified point based on cubic spline interpolation in 1D.
 *
 * @return (F_r) F_r = -\nabla \Phi(r)
 * @return (Phi) \Phi(r)
 * @param (rr) position of the specified point
 * @param (logrmin) minimum of r-table in logarithmic scale
 * @param (invlogrbin) inverse of logrbin
 * @param (num) number of data points
 * @param (xx) position of data points
 * @param (yy) value of data points and coefficients in cubic spline interpolation
 */
__device__ __forceinline__
void calcExternalForce_spherical(real * RESTRICT F_r, real * RESTRICT Phi, const real rr, const real logrmin, const real invlogrbin, const int num, real * RESTRICT xx, pot2 * RESTRICT yy)
{
  int ii = (int)FLOOR((LOG10(rr) - logrmin) * invlogrbin);
  ii = (ii >=       0 ) ? ii :       0;
  ii = (ii < (num - 1)) ? ii : num - 2;

  const real xl = xx[ii];
  const real xr = xx[ii + 1];

  const real dx = xr - xl;
  const real dx2 = dx * dx;
  const real dxinv = RSQRT(dx2);
  const real aa = (rr - xl) * dxinv;

  const pot2 yl = yy[ii];
  const pot2 yr = yy[ii + 1];

  *Phi =  ((UNITY - aa) * yl.val + aa * (yr.val + (aa - UNITY) * dx2 * ((aa + UNITY) * yr.dr2 - (aa - TWO) * yl.dr2) * ONE_SIXTH));
  *F_r = -((yr.val - yl.val) * dxinv + dx * ((THREE * aa * aa - UNITY) * yr.dr2 - (TWO + THREE * aa * (aa - TWO)) * yl.dr2) * ONE_SIXTH);
}


/**
 * @fn calcExternalGravity_kernel
 *
 * @brief Calculate gravitational force from external potential field.
 *
 * @param (logrmin) minimum of r-table in logarithmic scale
 * @param (invlogrbin) inverse of logrbin
 * @param (num) number of data points
 * @param (xx) position of data points
 * @param (yy) value of data points and coefficients in cubic spline interpolation
 */
__global__ void calcExternalGravity_kernel
(acceleration * RESTRICT acc, READ_ONLY position * RESTRICT pos,
#ifdef  BLOCK_TIME_STEP
 const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo,
#endif//BLOCK_TIME_STEP
 const real logrmin, const real invlogrbin, const int num, real * RESTRICT xx, pot2 * RESTRICT yy)
{
#ifdef  BLOCK_TIME_STEP
  const int lane    = THREADIDX_X1D & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);

  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];

  if( lane < info.num ){
    const int ii = info.head + lane;
#else///BLOCK_TIME_STEP
    const int ii = GLOBALIDX_X1D;
#endif//BLOCK_TIME_STEP

    /** load particle position and acceleration */
    const position pi = pos[ii];
    acceleration ai = acc[ii];

    /** set current location */
    const real r2 = 1.0e-30f + pi.x * pi.x + pi.y * pi.y + pi.z * pi.z;
    const real rinv = RSQRT(r2);
    const real rr = r2 * rinv;

    /** evaluate gravitational field from the external potential field */
    real F_r, Phi;
    calcExternalForce_spherical(&F_r, &Phi, rr, logrmin, invlogrbin, num, xx, yy);
    F_r *= rinv;

    /** accumulate the external force */
    ai.x   += F_r * ai.x;
    ai.y   += F_r * ai.y;
    ai.z   += F_r * ai.z;
    ai.pot += Phi;

    /** store acceleration */
    acc[ii] = ai;
#ifdef  BLOCK_TIME_STEP
  }/* if( lane < info.num ){ */
#endif//BLOCK_TIME_STEP
}


/**
 * @fn calcExternalGravity_dev
 *
 * @brief Calculate gravitational force from external potential field.
 */
extern "C"
void calcExternalGravity_dev
  (const iparticle pi, const potential_field sphe
#ifdef  BLOCK_TIME_STEP
   , const int thrd, const int grpNum, laneinfo * RESTRICT laneInfo
#else///BLOCK_TIME_STEP
   , const int Ni
#endif//BLOCK_TIME_STEP
   )
{
  __NOTE__("%s\n", "start");


#ifndef BLOCK_TIME_STEP
  const int thrd = NTHREADS;
#endif//BLOCK_TIME_STEP

  /** evaluate gravitational force by external potential field */
#ifdef  BLOCK_TIME_STEP
  int Nrem = BLOCKSIZE(grpNum, NWARP * NGROUPS);
#else///BLOCK_TIME_STEP
  int Nrem = BLOCKSIZE(Ni, NTHREADS);
#endif//BLOCK_TIME_STEP

  if( Nrem <= MAX_BLOCKS_PER_GRID ){
    /** when grid splitting is not required... */
#   if  defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
    if( grpNum != 0 )
#endif//defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
      calcExternalGravity_kernel<<<Nrem, thrd>>>
	(pi.acc, pi.pos,
#ifdef  BLOCK_TIME_STEP
	 BLOCKSIZE(grpNum, NGROUPS) * NGROUPS, laneInfo,
#endif//BLOCK_TIME_STEP
	 sphe.logrmin, sphe.invlogrbin, sphe.num, sphe.rad, sphe.Phi);
  }/* if( Nrem <= MAX_BLOCKS_PER_GRID ){ */
  else{
    /** when grid splitting is required... */
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

#ifdef  BLOCK_TIME_STEP
      int Nsub = Nblck * NWARP * NGROUPS;
#else///BLOCK_TIME_STEP
      int Nsub = Nblck * NTHREADS;
#endif//BLOCK_TIME_STEP

      calcExternalGravity_kernel<<<Nblck, thrd>>>
	(
#ifdef  BLOCK_TIME_STEP
	 pi.acc, pi.pos, BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx],
#else///BLOCK_TIME_STEP
	 &pi.acc[hidx], &pi.pos[hidx],
#endif//BLOCK_TIME_STEP
	 sphe.logrmin, sphe.invlogrbin, sphe.num, sphe.rad, sphe.Phi);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("calcExternalGravity_kernel");


  __NOTE__("%s\n", "end");
}
