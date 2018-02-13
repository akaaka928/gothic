/**
 * @file potential_dev.cu
 *
 * @brief Source code for calculating gravity from external fixed-potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/02/13 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <helper_cuda.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
/* The maximum number of elements in a chunk is 2^32-1 which is equal to 4,294,967,295 */
/* The maximum size for any chunk is 4GB */
#define MAXIMUM_CHUNK_SIZE      ((hsize_t)1 << 31)
#define MAXIMUM_CHUNK_SIZE_4BIT ((hsize_t)1 << 30)
#define MAXIMUM_CHUNK_SIZE_8BIT ((hsize_t)1 << 29)
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "cudalib.h"
#include "name.h"

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
static inline muse allocSphericalPotentialTable_dev
(pot2 **Phi,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real **rad,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const int Nr)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};

  mycudaMalloc((void **)Phi, Nr * sizeof(pot2));  alloc.device += Nr * sizeof(pot2);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaMalloc((void **)rad, Nr * sizeof(real));  alloc.device += Nr * sizeof(real);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn  freeSphericalPotentialTable_dev
 *
 * @brief Memory deallocation on the accelerator device for external fixed-potential field by spherical components.
 */
extern "C"
void  freeSphericalPotentialTable_dev
(pot2  *Phi
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
, real  *rad
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 )
{
  __NOTE__("%s\n", "start");

  mycudaFree(Phi);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaFree(rad);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
}


/**
 * @fn allocSphericalPotentialTable_hst
 *
 * @brief Memory allocation on the host for external fixed-potential field by spherical components.
 */
static inline muse allocSphericalPotentialTable_hst
(pot2 **Phi,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real **rad,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const int Nr)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};

  mycudaMallocHost((void **)Phi, Nr * sizeof(pot2));  alloc.host += Nr * sizeof(pot2);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaMallocHost((void **)rad, Nr * sizeof(real));  alloc.host += Nr * sizeof(real);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn  freeSphericalPotentialTable_hst
 *
 * @brief Memory deallocation on the host for external fixed-potential field by spherical components.
 */
static inline void  freeSphericalPotentialTable_hst
(pot2  *Phi
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 , real  *rad
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 )
{
  __NOTE__("%s\n", "start");

  mycudaFreeHost(Phi);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaFreeHost(rad);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
}


/**
 * @fn setSphericalPotentialTable_dev
 *
 * @brief Set external fixed-potential field by spherical components on the device.
 */
static inline void setSphericalPotentialTable_dev
(pot2 *Phi_hst, pot2 *Phi_dev,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real *rad_hst, real *rad_dev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const int Nr)
{
  __NOTE__("%s\n", "start");

  checkCudaErrors(cudaMemcpy(Phi_dev, Phi_hst, sizeof(pot2) * Nr, cudaMemcpyHostToDevice));
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  checkCudaErrors(cudaMemcpy(rad_dev, rad_hst, sizeof(real) * Nr, cudaMemcpyHostToDevice));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
}


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
/**
 * @fn allocDiskPotentialTable_dev
 *
 * @brief Memory allocation on the accelerator device for external fixed-potential field by disk components.
 */
static inline muse allocDiskPotentialTable_dev
(real **Phi,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real **RR, real **zz,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 disk_grav **FRz,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const disk_potential disk)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaMalloc((void **)RR , disk.maxLev *  disk.NR                      * sizeof(real));  alloc.device += disk.maxLev *  disk.NR                      * sizeof(real);
  mycudaMalloc((void **)zz , disk.maxLev *                  disk.Nz      * sizeof(real));  alloc.device += disk.maxLev *                  disk.Nz      * sizeof(real);
  mycudaMalloc((void **)Phi, disk.maxLev * (disk.NR + 1) * (disk.Nz + 1) * sizeof(real));  alloc.device += disk.maxLev * (disk.NR + 1) * (disk.Nz + 1) * sizeof(real);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaMalloc((void **)Phi, (disk.NR + 1) * (disk.Nz + 1) * sizeof(real));  alloc.device += (disk.NR + 1) * (disk.Nz + 1) * sizeof(real);
  mycudaMalloc((void **)FRz, (disk.NR + 1) * (disk.Nz + 1) * sizeof(disk_grav));  alloc.device += (disk.NR + 1) * (disk.Nz + 1) * sizeof(disk_grav);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freeDiskPotentialTable_dev
 *
 * @brief Memory deallocation on the accelerator device for external fixed-potential field by disk components.
 */
extern "C"
void  freeDiskPotentialTable_dev
(real  *Phi,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real  *RR, real  *zz
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 disk_grav *FRz
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 )
{
  __NOTE__("%s\n", "start");

  mycudaFree(Phi);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaFree(RR);
  mycudaFree(zz);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaFree(FRz);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
}


/**
 * @fn allocDiskPotentialTable_hst
 *
 * @brief Memory allocation on the host for external fixed-potential field by disk components.
 */
static inline muse allocDiskPotentialTable_hst
(real **Phi,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real **RR, real **zz,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 disk_grav **FRz,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const disk_potential disk)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaMallocHost((void **)RR , disk.maxLev *  disk.NR                      * sizeof(real));  alloc.device += disk.maxLev *  disk.NR                      * sizeof(real);
  mycudaMallocHost((void **)zz , disk.maxLev *                  disk.Nz      * sizeof(real));  alloc.device += disk.maxLev *                  disk.Nz      * sizeof(real);
  mycudaMallocHost((void **)Phi, disk.maxLev * (disk.NR + 1) * (disk.Nz + 1) * sizeof(real));  alloc.device += disk.maxLev * (disk.NR + 1) * (disk.Nz + 1) * sizeof(real);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaMallocHost((void **)Phi, (disk.NR + 1) * (disk.Nz + 1) * sizeof(real));  alloc.device += (disk.NR + 1) * (disk.Nz + 1) * sizeof(real);
  mycudaMallocHost((void **)FRz, (disk.NR + 1) * (disk.Nz + 1) * sizeof(disk_grav));  alloc.device += (disk.NR + 1) * (disk.Nz + 1) * sizeof(disk_grav);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freeDiskPotentialTable_hst
 *
 * @brief Memory deallocation on the host for external fixed-potential field by disk components.
 */
static inline void  freeDiskPotentialTable_hst
(real  *Phi,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real  *RR, real  *zz
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 disk_grav *FRz
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 )
{
  __NOTE__("%s\n", "start");

  mycudaFreeHost(Phi);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaFreeHost(RR);
  mycudaFreeHost(zz);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  mycudaFreeHost(FRz);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
}


/**
 * @fn setDiskPotentialTable_dev
 *
 * @brief Set external fixed-potential field by disk components on the device.
 */
static inline void setDiskPotentialTable_dev
(real *Phi_hst, real *Phi_dev,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real *RR_hst, real *RR_dev, real *zz_hst, real *zz_dev,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 disk_grav *FRz_hst, disk_grav *FRz_dev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const disk_potential disk)
{
  __NOTE__("%s\n", "start");

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  checkCudaErrors(cudaMemcpy(Phi_dev, Phi_hst, sizeof(real) * disk.maxLev * (disk.NR + 1) * (disk.Nz + 1), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy( RR_dev,  RR_hst, sizeof(real) * disk.maxLev *  disk.NR                     , cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy( zz_dev,  zz_hst, sizeof(real) * disk.maxLev *                  disk.Nz     , cudaMemcpyHostToDevice));
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  checkCudaErrors(cudaMemcpy(Phi_dev, Phi_hst, sizeof(real) * (disk.NR + 1) * (disk.Nz + 1), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(FRz_dev, FRz_hst, sizeof(disk_grav) * (disk.NR + 1) * (disk.Nz + 1), cudaMemcpyHostToDevice));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
}
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK


#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
/**
 * @fn bisec
 *
 * @brief Execute bisection.
 *
 * @param (val) the target value
 * @param (num) number of data points
 * @param (tab) array contains data points
 * @return (ratio) parameter for linear interpolation
 * @return lower index of the corresponding data point
 */
__device__ __forceinline__
int bisec(const real val, const int num, real * RESTRICT tab, real * RESTRICT ratio, real * RESTRICT dxinv)
{
  int ll =       0;
  int rr = num - 1;

  /** prohibit extraporation */
  if( val < tab[ll] + EPSILON ){    *ratio =  ZERO;    *dxinv = UNITY / (tab[ll + 1] - tab[ll]);    return (ll    );  }
  if( val > tab[rr] - EPSILON ){    *ratio = UNITY;    *dxinv = UNITY / (tab[rr + 1] - tab[rr]);    return (rr - 1);  }

  while( true ){
    const int cc = (ll + rr) >> 1;

    if( (tab[cc] - val) * (tab[ll] - val) <= ZERO)      rr = cc;
    else                                                ll = cc;

    if( (1 + ll) == rr ){
      *dxinv = UNITY / (tab[rr] - tab[ll]);
      *ratio = (val - tab[ll]) * (*dxinv);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
  }/* while( true ){ */
}
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD


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
void calcExternalForce_spherical
(real * RESTRICT F_r, real * RESTRICT Phi, const real rr, const int num,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real * RESTRICT xx,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const real logrmin, const real invlogrbin,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 pot2 * RESTRICT yy)
{
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  real aa, dxinv;
  int ii = bisec(rr, num, xx, &aa, &dxinv);

  const real xl = xx[ii];
  const real xr = xx[ii + 1];
  const pot2 yl = yy[ii];
  const pot2 yr = yy[ii + 1];

  const real dx = xr - xl;
  const real dx2 = dx * dx;

  *Phi =  ((UNITY - aa) * yl.val + aa * (yr.val + (aa - UNITY) * dx2 * ((aa + UNITY) * yr.dr2 - (aa - TWO) * yl.dr2) * ONE_SIXTH));
  *F_r = -((yr.val - yl.val) * dxinv + dx * ((THREE * aa * aa - UNITY) * yr.dr2 - (TWO + THREE * aa * (aa - TWO)) * yl.dr2) * ONE_SIXTH);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  const real logr = (LOG10(rr) - logrmin) * invlogrbin;
  int ii = (int)FLOOR(logr);
  ii = (ii >=       0 ) ? ii :       0;
  ii = (ii < (num - 1)) ? ii : num - 2;
  const real aa = logr - (real)ii;

  const pot2 yl = yy[ii];
  const pot2 yr = yy[ii + 1];

  *Phi = (UNITY - aa) * yl.Phi + aa * yr.Phi;
  *F_r = (UNITY - aa) * yl.Fr  + aa * yr.Fr;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
}


/**
 * @fn calcExternalGravity_kernel
 *
 * @brief Calculate gravitational force from external potential field by spherical component(s).
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
 const int num,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real * RESTRICT xx,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const real logrmin, const real invlogrbin,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 pot2 * RESTRICT yy)
{
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  extern __shared__ real xx_sm[];
  for(int ii = THREADIDX_X1D; ii < num; ii += NTHREADS)
    if( ii < num )
      xx_sm[ii] = xx[ii];
  __syncthreads();
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

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

    /** load particle position */
    const position pi = pos[ii];

    /** set current location */
    const real r2 = 1.0e-30f + pi.x * pi.x + pi.y * pi.y + pi.z * pi.z;
    const real rinv = RSQRT(r2);
    const real rr = r2 * rinv;

    /** evaluate gravitational field from the external potential field */
    real F_r, pot;
    calcExternalForce_spherical
      (&F_r, &pot, rr, num,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
       xx_sm,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
       logrmin, invlogrbin,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
       yy);
    F_r *= rinv;

    /** calculate the external force by spherical components */
    acceleration ai;
    ai.x   = F_r * pi.x;
    ai.y   = F_r * pi.y;
    ai.z   = F_r * pi.z;
    ai.pot = pot;

    /** store acceleration */
#ifndef SET_EXTERNAL_POTENTIAL_FIELD_DISK
    acc[ii] = ai;
#else///SET_EXTERNAL_POTENTIAL_FIELD_DISK
    atomicAdd(&(acc[ii].x  ), ai.x  );
    atomicAdd(&(acc[ii].y  ), ai.y  );
    atomicAdd(&(acc[ii].z  ), ai.z  );
    atomicAdd(&(acc[ii].pot), ai.pot);
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#ifdef  BLOCK_TIME_STEP
  }/* if( lane < info.num ){ */
#endif//BLOCK_TIME_STEP
}


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK

#ifdef  DOUBLE_PRECISION
#define THREE_HALVES (1.5)
#else///DOUBLE_PRECISION
#define THREE_HALVES (1.5f)
#endif//DOUBLE_PRECISION

/**
 * @fn calcExternalDiskGravity_kernel
 *
 * @brief Calculate gravitational force from external potential field by disk component(s).
 *
 * @param (logrmin) minimum of r-table in logarithmic scale
 * @param (invlogrbin) inverse of logrbin
 * @param (num) number of data points
 * @param (xx) position of data points
 * @param (yy) value of data points and coefficients in cubic spline interpolation
 */
__global__ void calcExternalDiskGravity_kernel
(acceleration * RESTRICT acc, READ_ONLY position * RESTRICT pos,
#ifdef  BLOCK_TIME_STEP
 const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo,
#endif//BLOCK_TIME_STEP
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const int maxLev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const int NR, const int Nz, const real hh, const real hinv, READ_ONLY real * RESTRICT Phi,
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 READ_ONLY disk_grav * RESTRICT FRz,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const int num,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real * RESTRICT xx,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 const real logrmin, const real invlogrbin,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 pot2 * RESTRICT yy)
{
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  extern __shared__ real xx_sm[];
  for(int ii = THREADIDX_X1D; ii < num; ii += NTHREADS)
    if( ii < num )
      xx_sm[ii] = xx[ii];
  __syncthreads();
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

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

    /** load particle position */
    const position pi = pos[ii];

    /** set current location */
    const real R2 = 1.0e-30f + pi.x * pi.x + pi.y * pi.y;
    const real Rinv = RSQRT(R2);
    const real RR = R2 * Rinv;
    const real z2 = 1.0e-30f + pi.z * pi.z;
    const real zinv = RSQRT(z2);
    const real zz = z2 * zinv;

    acceleration ai;

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    /** find appropriate grid level */
    const int levR = (int)FLOOR(LOG2(hh * ((real)NR - THREE_HALVES) * Rinv));
    const int levz = (int)FLOOR(LOG2(hh * ((real)Nz - THREE_HALVES) * zinv));
    int lev = (levz < levR) ? levz : levR;
    lev = (lev < (maxLev - 1)) ? lev : (maxLev - 1);

    if( lev >= 0 ){
      /** find the corresponding grid location */
      const int jj = 1 + (int)FMAX(FLOOR(LDEXP(RR * hinv, lev) - HALF), ZERO);
      const int kk = 1 + (int)FMAX(FLOOR(LDEXP(zz * hinv, lev) - HALF), ZERO);

      /** load disk potential at 5 points */
      const real mc = Phi[INDEX(maxLev, NR + 1, Nz + 1, lev, jj - 1, kk    )];
      const real cm = Phi[INDEX(maxLev, NR + 1, Nz + 1, lev, jj    , kk - 1)];
      const real cc = Phi[INDEX(maxLev, NR + 1, Nz + 1, lev, jj    , kk    )];
      const real cp = Phi[INDEX(maxLev, NR + 1, Nz + 1, lev, jj    , kk + 1)];
      const real pc = Phi[INDEX(maxLev, NR + 1, Nz + 1, lev, jj + 1, kk    )];

      /** prepare for interpolation */
      const real dlinv = LDEXP(hinv, lev - 1);/**< = 2^(lev - 1) / h = 1 / (2 h 2^-lev) = 1 / (2 dR) = 1 / (2 dz) */
      const real dPhidR = (pc - mc) * dlinv;
      const real dPhidz = (cp - cm) * dlinv;

      const real dl = LDEXP(hh, -lev);
      const real R0 = dl * (HALF + (real)(jj - 1));
      const real z0 = dl * (HALF + (real)(kk - 1));
      const real pot = cc + dPhidR * (RR - R0) + dPhidz * (zz - z0);

      /** evaluate gravitational field from the external potential field */
      const real FR = Rinv * dPhidR;
      const real Fz = zinv * dPhidz;

      /** calculate the external force by disk components */
      ai.x = FR * pi.x;
      ai.y = FR * pi.y;
      ai.z = Fz * pi.z;
      ai.pot = pot;
    }/* if( lev >= 0 ){ */
    else{
      /** set current location */
      const real r2 = R2 + pi.z * pi.z;
      const real rinv = RSQRT(r2);
      const real rr = r2 * rinv;

      /** evaluate gravitational field from the external potential field */
      real F_r, pot;
      calcExternalForce_spherical(&F_r, &pot, rr, num, xx_sm, yy);
      F_r *= rinv;

      /** calculate the external force by disk components */
      ai.x   = F_r * pi.x;
      ai.y   = F_r * pi.y;
      ai.z   = F_r * pi.z;
      ai.pot = pot;
    }/* else{ */
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    const real Rref = RR * hinv;
    const real zref = zz * hinv;
    const int jj = (int)FLOOR(Rref);
    const int kk = (int)FLOOR(zref);
    if( (jj < NR) && (kk < Nz) ){
      const real Phi_ll = Phi[INDEX2D(NR + 1, Nz + 1,     jj,     kk)];
      const real Phi_lu = Phi[INDEX2D(NR + 1, Nz + 1,     jj, 1 + kk)];
      const real Phi_ul = Phi[INDEX2D(NR + 1, Nz + 1, 1 + jj,     kk)];
      const real Phi_uu = Phi[INDEX2D(NR + 1, Nz + 1, 1 + jj, 1 + kk)];

      const disk_grav FRz_ll = FRz[INDEX2D(NR + 1, Nz + 1,     jj,     kk)];
      const disk_grav FRz_lu = FRz[INDEX2D(NR + 1, Nz + 1,     jj, 1 + kk)];
      const disk_grav FRz_ul = FRz[INDEX2D(NR + 1, Nz + 1, 1 + jj,     kk)];
      const disk_grav FRz_uu = FRz[INDEX2D(NR + 1, Nz + 1, 1 + jj, 1 + kk)];

      const real aa = Rref - (real)jj;
      const real bb = zref - (real)kk;

      const real pot =  ((UNITY - bb) * Phi_ll   + bb * Phi_lu  ) * (UNITY - aa) + ((UNITY - bb) * Phi_ul   + bb * Phi_uu  ) * aa;
      const real FR  = (((UNITY - bb) * FRz_ll.R + bb * FRz_lu.R) * (UNITY - aa) + ((UNITY - bb) * FRz_ul.R + bb * FRz_uu.R) * aa) * Rinv;
      const real Fz  = (((UNITY - bb) * FRz_ll.z + bb * FRz_lu.z) * (UNITY - aa) + ((UNITY - bb) * FRz_ul.z + bb * FRz_uu.z) * aa) * zinv;

      /** calculate the external force by disk components */
      ai.x = FR * pi.x;
      ai.y = FR * pi.y;
      ai.z = Fz * pi.z;
      ai.pot = pot;
    }/* if( (jj < NR) && (kk < Nz) ){ */
    else{
      /** set current location */
      const real r2 = R2 + pi.z * pi.z;
      const real rinv = RSQRT(r2);
      const real rr = r2 * rinv;

      /** evaluate gravitational field from the external potential field */
      real F_r, pot;
      calcExternalForce_spherical(&F_r, &pot, rr, num, logrmin, invlogrbin, yy);
      F_r *= rinv;

      /** calculate the external force by disk components */
      ai.x   = F_r * pi.x;
      ai.y   = F_r * pi.y;
      ai.z   = F_r * pi.z;
      ai.pot = pot;
    }/* else{ */
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

    /** store acceleration */
    atomicAdd(&(acc[ii].x  ), ai.x  );
    atomicAdd(&(acc[ii].y  ), ai.y  );
    atomicAdd(&(acc[ii].z  ), ai.z  );
    atomicAdd(&(acc[ii].pot), ai.pot);

#ifdef  BLOCK_TIME_STEP
  }/* if( lane < info.num ){ */
#endif//BLOCK_TIME_STEP
}
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK


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
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
   , const disk_potential disk
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
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
      calcExternalGravity_kernel
	<<<Nrem, thrd
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	, sphe.num * sizeof(real)
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	>>>
	(pi.acc_ext, pi.pos,
#ifdef  BLOCK_TIME_STEP
	 BLOCKSIZE(grpNum, NGROUPS) * NGROUPS, laneInfo,
#endif//BLOCK_TIME_STEP
	 sphe.num,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 sphe.rad,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 sphe.logrmin, sphe.invlogrbin,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 sphe.Phi);
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

      calcExternalGravity_kernel
	<<<Nblck, thrd
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	, sphe.num * sizeof(real)
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	>>>
	(
#ifdef  BLOCK_TIME_STEP
	 pi.acc_ext, pi.pos, BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx],
#else///BLOCK_TIME_STEP
	 &pi.acc_ext[hidx], &pi.pos[hidx],
#endif//BLOCK_TIME_STEP
	 sphe.num,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 sphe.rad,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 sphe.logrmin, sphe.invlogrbin,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 sphe.Phi);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("calcExternalGravity_kernel");


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
#ifdef  BLOCK_TIME_STEP
  Nrem = BLOCKSIZE(grpNum, NWARP * NGROUPS);
#else///BLOCK_TIME_STEP
  Nrem = BLOCKSIZE(Ni, NTHREADS);
#endif//BLOCK_TIME_STEP

  if( Nrem <= MAX_BLOCKS_PER_GRID ){
    /** when grid splitting is not required... */
#   if  defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
    if( grpNum != 0 )
#endif//defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
      calcExternalDiskGravity_kernel
	<<<Nrem, thrd
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	, disk.sphe.num * sizeof(real)
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	>>>
	(pi.acc_ext, pi.pos,
#ifdef  BLOCK_TIME_STEP
	 BLOCKSIZE(grpNum, NGROUPS) * NGROUPS, laneInfo,
#endif//BLOCK_TIME_STEP
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.maxLev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.NR, disk.Nz, disk.hh, disk.hinv, disk.Phi,
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.FRz,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.sphe.num,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.sphe.rad,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.sphe.logrmin, disk.sphe.invlogrbin,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.sphe.Phi);
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

      calcExternalDiskGravity_kernel
	<<<Nblck, thrd
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	, disk.sphe.num * sizeof(real)
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	>>>
	(
#ifdef  BLOCK_TIME_STEP
	 pi.acc_ext, pi.pos, BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx],
#else///BLOCK_TIME_STEP
	 &pi.acc_ext[hidx], &pi.pos[hidx],
#endif//BLOCK_TIME_STEP
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.maxLev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.NR, disk.Nz, disk.hh, disk.hinv, disk.Phi,
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.FRz,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.sphe.num,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.sphe.rad,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.sphe.logrmin, disk.sphe.invlogrbin,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	 disk.sphe.Phi);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

#if 0
  cudaDeviceSynchronize();
  exit(0);
#endif

  getLastCudaError("calcExternalDiskGravity_kernel");
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK


  __NOTE__("%s\n", "end");
}


/**
 * @fn  readFixedPotentialTableSpherical
 *
 * @brief Read fixed potential field (of spherical symmetric components) represented in cubic spline interpolation for GOTHIC.
 *
 * @param (unit) unit system of the potential field
 * @param (cfg) name of the input file
 * @return (pot_tbl) superposed potential field for cubic spline interpolation (only for spherical symmetric components)
 * @return (rad) radius for cubic spline interpolation
 * @return (Phi) potential field for cubic spline interpolation
 */
extern "C"
muse  readFixedPotentialTableSpherical
(const int unit, char file[], potential_field *pot_tbl, pot2 **Phi
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 , real **rad
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#ifdef  USE_HDF5_FORMAT
 , hdf5struct type
#endif//USE_HDF5_FORMAT
 )
{
  __NOTE__("%s\n", "start");


  char cfgfile[128];
  sprintf(cfgfile, "%s/%s.%s.cfg", DATAFOLDER, file, "ext_pot_sphe");
  FILE *fp_cfg;
  fp_cfg = fopen(cfgfile, "r");
  if( fp_cfg == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", cfgfile);  }

  int Nread;
  bool success_cfg = true;
  success_cfg &= (1 == fscanf(fp_cfg, "%d", &Nread));

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  if( Nread != 1 ){
    __KILL__(stderr, "ERROR: Nread = %d; however, reading and superposing multiple potential tables in GOTHIC are not yet supported.\n", Nread);
  }/* if( Nread != 1 ){ */
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  /* open an existing file with read only option */
  char filename[128];
#ifdef  USE_HDF5_FORMAT

  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "ext_pot");
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* read attribute data */
  hid_t attribute;
  /* read flag about DOUBLE_PRECISION */
  int useDP;
  attribute = H5Aopen(target, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* read flag about unit system */
  int unit_pot;
  hid_t group = H5Gopen(target, "unit_system", H5P_DEFAULT);
  attribute = H5Aopen(group, "unit", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &unit_pot));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));

  /* simple error checks */
#ifdef  DOUBLE_PRECISION
  if( useDP != 1 ){    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);  }
#else///DOUBLE_PRECISION
  if( useDP != 0 ){    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);  }
#endif//DOUBLE_PRECISION
  if( unit_pot != unit ){
    __KILL__(stderr, "ERROR: unit system of the potential field (%d) differs with that in the simulation run (%d)\n", unit_pot, unit);
  }/* if( unit_pot != unit ){ */

#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  group = H5Gopen(target, "spherical", H5P_DEFAULT);
  /* read # of data points */
  attribute = H5Aopen(group, "num", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &pot_tbl->num));
  chkHDF5err(H5Aclose(attribute));
  const int num = pot_tbl->num;
  /* read log_10(r_min) */
  attribute = H5Aopen(group, "logrmin", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, type.real, &pot_tbl->logrmin));
  chkHDF5err(H5Aclose(attribute));
  /* read logrbin */
  attribute = H5Aopen(group, "logrbin", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, type.real, &pot_tbl->logrbin));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  /* memory allocation on the accelerator device */
  muse alloc_tbl;
  pot2 *Phi_hst;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  real *rad_hst;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  alloc_tbl = allocSphericalPotentialTable_dev(Phi, num);
  allocSphericalPotentialTable_hst (&Phi_hst, num);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  /* memory allocation on the host as a temporary buffer */
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  pot2 *Phi_tmp;
  if( Nread == 1 )
    Phi_tmp = Phi_hst;
  else{
    Phi_tmp = (pot2 *)malloc(num * sizeof(pot2));    if( Phi_tmp == NULL ){      __KILL__(stderr, "ERROR: failure to allocate Phi_tmp\n");    }

    const pot2 zero = {ZERO, ZERO};
    for(int ii = 0; ii < num; ii++)
      Phi_hst[ii] = zero;
  }/* else{ */
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  for(int ii = 0; ii < Nread; ii++){
    hid_t dataset;

    /* open an existing group */
    char list[64];
    success_cfg &= (1 == fscanf(fp_cfg, "%s", list));
    group = H5Gopen(target, list, H5P_DEFAULT);

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    /* read # of data points */
    attribute = H5Aopen(group, "num", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &pot_tbl->num));
    chkHDF5err(H5Aclose(attribute));

    const int num = pot_tbl->num;

    alloc_tbl = allocSphericalPotentialTable_dev(Phi, rad, num);
    allocSphericalPotentialTable_hst(&Phi_hst, &rad_hst, num);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    /* read radius */
    dataset = H5Dopen(group, "r", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rad_hst));
    chkHDF5err(H5Dclose(dataset));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

    /* read potential */
    dataset = H5Dopen(group, "Phi(r)", H5P_DEFAULT);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    chkHDF5err(H5Dread(dataset, type.pot2, H5S_ALL, H5S_ALL, H5P_DEFAULT, Phi_hst));
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    chkHDF5err(H5Dread(dataset, type.pot2, H5S_ALL, H5S_ALL, H5P_DEFAULT, Phi_tmp));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    chkHDF5err(H5Dclose(dataset));

    chkHDF5err(H5Gclose(group));

#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    if( Nread > 1 )
      for(int jj = 0; jj < num; jj++){
	Phi_hst[jj].Phi += Phi_tmp[jj].Phi;
	Phi_hst[jj].Fr  += Phi_tmp[jj].Fr;
      }/* for(int jj = 0; jj < num; jj++){ */
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  }/* for(int ii = 0; ii < Nread; ii++){ */
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  if( Nread != 1 )
    free(Phi_tmp);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  /* close the file */
  chkHDF5err(H5Fclose(target));

#else///USE_HDF5_FORMAT

  /* read numeric table for superposed spherical components */
  sprintf(filename, "%s/%s.%s.%s", DATAFOLDER, file, "pot", "sphe");
  FILE *fp;
  fp = fopen(filename, "rb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

  int unit_pot, num;
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  real rmin, rbin;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  bool success = true;
  size_t tmp;

  tmp = 1;  if( tmp != fread(&unit_pot, sizeof(int), tmp, fp) )    success = false;
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  tmp = 1;  if( tmp != fread(&num, sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&rmin, sizeof(real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&rbin, sizeof(real), tmp, fp) )    success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  /* simple error check */
  if( unit_pot != unit ){
    __KILL__(stderr, "ERROR: unit system of the potential field (%d) differs with that in the simulation run (%d)\n", unit_pot, unit);
  }/* if( unit_pot != unit ){ */

#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  pot_tbl->num = num;
  pot_tbl->logrmin = rmin;
  pot_tbl->logrbin = rbin;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  /* memory allocation on the accelerator device */
  muse alloc_tbl;
  pot2 *Phi_hst;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  real *rad_hst;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  alloc_tbl = allocSphericalPotentialTable_dev(rad, Phi, num);
  allocSphericalPotentialTable_hst(&rad_hst, &Phi_hst, num);

  tmp = num;  if( tmp != fread(rad_hst, sizeof(real), tmp, fp) )    success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  int list;
  success_cfg &= (1 == fscanf(fp_cfg, "%d", &list));
  if( (Nread == 1) && (list == READ_SUPERPOSED_TABLE_SPHE) ){
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = 1;    if( tmp != fread(&num, sizeof(int), tmp, fp) )      success = false;

    pot_tbl->num = num;
    alloc_tbl = allocSphericalPotentialTable_dev(rad, Phi, num);
    allocSphericalPotentialTable_hst(&rad_hst, &Phi_hst, num);

    tmp = num;    if( tmp != fread(rad_hst, sizeof(real), tmp, fp) )      success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = num;    if( tmp != fread(Phi_hst, sizeof(pot2), tmp, fp) )      success = false;

    if( success != true ){      __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);    }
    fclose(fp);
  }/* if( (Nread == 1) && (list == READ_SUPERPOSED_TABLE_SPHE) ){ */
  else{
    /* close the superposed file */
    if( success != true ){      __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);    }
    fclose(fp);

#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    pot2 *Phi_tmp;
    Phi_tmp = (pot2 *)malloc(num * sizeof(pot2));    if( Phi_tmp == NULL ){      __KILL__(stderr, "ERROR: failure to allocate Phi_tmp\n");    }

    const pot2 zero = {ZERO, ZERO};
    for(int ii = 0; ii < num; ii++)
      Phi_hst[ii] = zero;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

    /* open another file */
    for(int ii = 0; ii < Nread; ii++){
      sprintf(filename, "%s/%s.%s.%d", DATAFOLDER, file, "pot", list);
      fp = fopen(filename, "rb");

      tmp = 1;      if( tmp != fread(&unit_pot, sizeof(int), tmp, fp) )	success = false;
      tmp = 1;      if( tmp != fread(&num, sizeof(int), tmp, fp) )	success = false;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
      pot_tbl->num = num;
      alloc_tbl = allocSphericalPotentialTable_dev(rad, Phi, num);
      allocSphericalPotentialTable_hst(&rad_hst, &Phi_hst, num);

      tmp = num;      if( tmp != fread(rad_hst, sizeof(real), tmp, fp) )	success = false;
      tmp = num;      if( tmp != fread(Phi_hst, sizeof(pot2), tmp, fp) )	success = false;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
      tmp = 1;      if( tmp != fread(&rmin, sizeof(real), tmp, fp) )	success = false;
      tmp = 1;      if( tmp != fread(&rbin, sizeof(real), tmp, fp) )	success = false;

      tmp = num;      if( tmp != fread(Phi_tmp, sizeof(pot2), tmp, fp) )	success = false;

      for(int jj = 0; jj < num; jj++){
	Phi_hst[jj].val += Phi_tmp[jj].val;
	Phi_hst[jj].dr2 += Phi_tmp[jj].dr2;
      }/* for(int jj = 0; jj < *Nr; jj++){ */
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

      if( success != true ){	__KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);      }
      fclose(fp);

      if( ii < (Nread - 1) )
	success_cfg &= (1 == fscanf(fp_cfg, "%d", &list));
    }/* for(int ii = 0; ii < Nread; ii++){ */
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    free(Phi_tmp);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  }/* else{ */

#endif//USE_HDF5_FORMAT

  setSphericalPotentialTable_dev
    (Phi_hst, *Phi,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     rad_hst, *rad,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     pot_tbl->num);
  freeSphericalPotentialTable_hst
    (Phi_hst
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     , rad_hst
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     );

  pot_tbl->Phi = *Phi;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  pot_tbl->rad = *rad;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  pot_tbl->invlogrbin = UNITY / pot_tbl->logrbin;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  if( success_cfg != true ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", cfgfile);  }
  fclose(fp_cfg);

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  if( (pot_tbl->num * sizeof(real)) > (SMEM_SIZE_L1_PREF >> 1) )
    checkCudaErrors(cudaFuncSetCacheConfig(calcExternalGravity_kernel, cudaFuncCachePreferShared));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD


  __NOTE__("%s\n", "end");

  return (alloc_tbl);
}

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK

/**
 * @fn  readFixedPotentialTableDisk
 *
 * @brief Read fixed potential field (of disk components) represented in cubic spline interpolation for GOTHIC.
 *
 */
extern "C"
muse  readFixedPotentialTableDisk
(const int unit, char file[], real **Phi_dev, pot2 **Phi_sphe_dev,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 real **RR_dev, real **zz_dev, real **rad_sphe_dev,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 disk_grav **FRz_dev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
 disk_potential *disk
#ifdef  USE_HDF5_FORMAT
 , hdf5struct type
#endif//USE_HDF5_FORMAT
 )
{
  __NOTE__("%s\n", "start");

  muse alloc_disk, alloc_sphe;
  real *Phi_hst;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  real *RR_hst, *zz_hst;
  real *rad_sphe_hst;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  disk_grav *FRz_hst;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  pot2 *Phi_sphe_hst;

  /* open an existing file with read only option */
  char filename[128];
#ifdef  USE_HDF5_FORMAT
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "ext_disk");
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* read attribute data */
  hid_t attribute;
  /* read flag about DOUBLE_PRECISION */
  int useDP;
  attribute = H5Aopen(target, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* read flag about unit system */
  int unit_pot;
  hid_t group = H5Gopen(target, "unit_system", H5P_DEFAULT);
  attribute = H5Aopen(group, "unit", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &unit_pot));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));

  /* simple error checks */
#ifdef  DOUBLE_PRECISION
  if( useDP != 1 ){    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);  }
#else///DOUBLE_PRECISION
  if( useDP != 0 ){    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);  }
#endif//DOUBLE_PRECISION
  if( unit_pot != unit ){
    __KILL__(stderr, "ERROR: unit system of the potential field (%d) differs with that in the simulation run (%d)\n", unit_pot, unit);
  }/* if( unit_pot != unit ){ */


  /* read potential table of superposed disk components */
  group = H5Gopen(target, "2D", H5P_DEFAULT);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* read # of nested levels */
  attribute = H5Aopen(group, "maxLev", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &disk->maxLev));
  chkHDF5err(H5Aclose(attribute));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* read # of R grids */
  attribute = H5Aopen(group, "NR", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &disk->NR));
  chkHDF5err(H5Aclose(attribute));
  /* read # of z grids */
  attribute = H5Aopen(group, "Nz", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &disk->Nz));
  chkHDF5err(H5Aclose(attribute));
  /* read hh */
  attribute = H5Aopen(group, "hh", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, type.real, &disk->hh));
  chkHDF5err(H5Aclose(attribute));

  /* memory allocation on the accelerator device */
  alloc_disk = allocDiskPotentialTable_dev
    (Phi_dev,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     RR_dev, zz_dev,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     FRz_dev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     *disk);
  allocDiskPotentialTable_hst
    (&Phi_hst,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     &RR_hst, &zz_hst,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     &FRz_hst,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     *disk);

  hid_t dataset;
  /* read \Phi(R, z) */
  dataset = H5Dopen(group, "Phi", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, Phi_hst));
  chkHDF5err(H5Dclose(dataset));
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* read R */
  dataset = H5Dopen(group, "R", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, RR_hst));
  chkHDF5err(H5Dclose(dataset));
  /* write z */
  dataset = H5Dopen(group, "z", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, zz_hst));
  chkHDF5err(H5Dclose(dataset));
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* read F(R, z) */
  dataset = H5Dopen(group, "FRz", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.disk_grav, H5S_ALL, H5S_ALL, H5P_DEFAULT, FRz_hst));
  chkHDF5err(H5Dclose(dataset));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  chkHDF5err(H5Gclose(group));


  /* read spherical averaged potential profile */
  group = H5Gopen(target, "spherical", H5P_DEFAULT);
  /* read # of data points */
  attribute = H5Aopen(group, "num", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &disk->sphe.num));
  chkHDF5err(H5Aclose(attribute));
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* read log_10(r_min) */
  attribute = H5Aopen(group, "logrmin", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, type.real, &disk->sphe.logrmin));
  chkHDF5err(H5Aclose(attribute));
  /* read logrbin */
  attribute = H5Aopen(group, "logrbin", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, type.real, &disk->sphe.logrbin));
  chkHDF5err(H5Aclose(attribute));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  alloc_sphe = allocSphericalPotentialTable_dev
    (Phi_sphe_dev,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     rad_sphe_dev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     disk->sphe.num);
  allocSphericalPotentialTable_hst
    (&Phi_sphe_hst,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     &rad_sphe_hst,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     disk->sphe.num);

  /* read radius */
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  dataset = H5Dopen(group, "r", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rad_sphe_hst));
  chkHDF5err(H5Dclose(dataset));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* read potential */
  dataset = H5Dopen(group, "Phi(r)", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.pot2, H5S_ALL, H5S_ALL, H5P_DEFAULT, Phi_sphe_hst));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Gclose(group));

  /* close the file */
  chkHDF5err(H5Fclose(target));

#else///USE_HDF5_FORMAT

  /* read numeric table for superposed disk components */
  sprintf(filename, "%s/%s.%s.dat", DATAFOLDER, file, "ext_disk");
  FILE *fp;
  fp = fopen(filename, "rb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

  int unit_pot;

  bool success = true;
  size_t tmp;

  tmp = 1;  if( tmp != fread(&unit_pot, sizeof(int), tmp, fp) )    success = false;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  tmp = 1;  if( tmp != fread(&disk->maxLev, sizeof(int), tmp, fp) )    success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  tmp = 1;  if( tmp != fread(&disk->NR, sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&disk->Nz, sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&disk->hh, sizeof(real), tmp, fp) )    success = false;

  /* simple error check */
  if( unit_pot != unit ){
    __KILL__(stderr, "ERROR: unit system of the potential field (%d) differs with that in the simulation run (%d)\n", unit_pot, unit);
  }/* if( unit_pot != unit ){ */

  /* memory allocation on the accelerator device */
  alloc_disk = allocDiskPotentialTable_dev
    (Phi_dev,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     RR_dev, zz_dev,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     FRz_dev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     *disk);
  allocDiskPotentialTable_hst
    (&Phi_hst,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     &RR_hst, &zz_hst,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     FRz_hst,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     *disk);

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  tmp = disk->maxLev *  disk->NR                      ;  if( tmp != fread( RR_hst, sizeof(real), tmp, fp) )    success = false;
  tmp = disk->maxLev *                   disk->Nz     ;  if( tmp != fread( zz_hst, sizeof(real), tmp, fp) )    success = false;
  tmp = disk->maxLev * (disk->NR + 1) * (disk->Nz + 1);  if( tmp != fread(Phi_hst, sizeof(real), tmp, fp) )    success = false;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  tmp = (disk->NR + 1) * (disk->Nz + 1);  if( tmp != fread(Phi_hst, sizeof(real), tmp, fp) )    success = false;
  tmp = (disk->NR + 1) * (disk->Nz + 1);  if( tmp != fread(FRz_hst, sizeof(disk_grav), tmp, fp) )    success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  tmp = 1;  if( tmp != fread(&(disk->sphe.num), sizeof(int), tmp, fp) )    success = false;
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  tmp = 1;  if( tmp != fread(&(disk->sphe.logrmin), sizeof(real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&(disk->sphe.logrbin), sizeof(real), tmp, fp) )    success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  alloc_sphe = allocSphericalPotentialTable_dev
    (Phi_sphe_dev,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     rad_sphe_dev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     disk->sphe.num);
  allocSphericalPotentialTable_hst
    (&Phi_sphe_hst,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     &rad_sphe_hst,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     disk->sphe.num);

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  tmp = disk->sphe.num;  if( tmp != fwrite(rad_sphe_hst, sizeof(real), tmp, fp) )    success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  tmp = disk->sphe.num;  if( tmp != fwrite(Phi_sphe_hst, sizeof(pot2), tmp, fp) )    success = false;

  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }
  fclose(fp);

#endif//USE_HDF5_FORMAT

  setDiskPotentialTable_dev
    (Phi_hst, *Phi_dev,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     RR_hst, *RR_dev, zz_hst, *zz_dev,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     FRz_hst, *FRz_dev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     *disk);

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  disk->RR  = * RR_dev;
  disk->zz  = * zz_dev;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  disk->FRz = *FRz_dev;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  disk->Phi = *Phi_dev;
  disk->hinv = UNITY / disk->hh;

  setSphericalPotentialTable_dev
    (Phi_sphe_hst, *Phi_sphe_dev,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     rad_sphe_hst, *rad_sphe_dev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     disk->sphe.num);

  disk->sphe.Phi = *Phi_sphe_dev;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  disk->sphe.rad = *rad_sphe_dev;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  disk->sphe.invlogrbin = UNITY / disk->sphe.logrbin;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  if( (disk->sphe.num * sizeof(real)) > (SMEM_SIZE_L1_PREF >> 1) )
    checkCudaErrors(cudaFuncSetCacheConfig(calcExternalDiskGravity_kernel, cudaFuncCachePreferShared));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  freeDiskPotentialTable_hst
    (Phi_hst,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     RR_hst, zz_hst
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     FRz_hst
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     );
  freeSphericalPotentialTable_hst
    (Phi_sphe_hst
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     , rad_sphe_hst
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     );


  __NOTE__("%s\n", "end");
  const muse alloc_tbl = {alloc_disk.host + alloc_sphe.host, alloc_disk.device + alloc_sphe.device};

  return (alloc_tbl);
}
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
