/**
 * @file seb_dev.cu
 *
 * @brief Source code to generate Smallest Enclosing Ball on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/04/04 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#ifndef SEB_DEV_CU
#define SEB_DEV_CU


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <helper_cuda.h>

#include "macro.h"
#include "cudalib.h"

#include "walk_dev.h"


#ifdef  DOUBLE_PRECISION
#define SEB_EPS  (1.0e-15)
#define SEB_EPS2 (1.0e-30)
#define SEB_TINY (DBL_MIN)
#else///DOUBLE_PRECISION
#define SEB_EPS  (1.0e-6f)
#define SEB_EPS2 (1.0e-12f)
#define SEB_TINY (FLT_MIN)
#endif//DOUBLE_PRECISION


/**
 * @fn getMinLocRealTsub
 *
 * @brief Get minimum value with location within a group of TSUB threads.
 * @detail implicit synchronization within TSUB (<= 32) threads is assumed
 * continuous NWARP threads have the same value as input
 */
#ifdef  USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ real getMinLocRealTsub(const real min, int *loc)
#else///USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ real getMinLocRealTsub(      real min, int *loc, volatile uint_real * smem, volatile int * sidx, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC
{
#ifdef  USE_WARP_SHUFFLE_FUNC
  real val =  min;
  int  idx = *loc;
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  int  buf;
  tmp = __SHFL_XOR(val,      NWARP, TSUB);  buf = __SHFL_XOR(idx,      NWARP, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= ( 4 * NWARP)
  tmp = __SHFL_XOR(val,  2 * NWARP, TSUB);  buf = __SHFL_XOR(idx,  2 * NWARP, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= ( 8 * NWARP)
  tmp = __SHFL_XOR(val,  4 * NWARP, TSUB);  buf = __SHFL_XOR(idx,  4 * NWARP, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= (16 * NWARP)
  tmp = __SHFL_XOR(val,  8 * NWARP, TSUB);  buf = __SHFL_XOR(idx,  8 * NWARP, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB == (32 * NWARP)
  tmp = __SHFL_XOR(val, 16 * NWARP, TSUB);  buf = __SHFL_XOR(idx, 16 * NWARP, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
  *loc =  __SHFL(idx, 0, TSUB);
  return (__SHFL(val, 0, TSUB));
#else///USE_WARP_SHUFFLE_FUNC
  smem[tidx].r =  min;
  sidx[tidx]   = *loc;
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = smem[tidx ^ (     NWARP)].r;  if( tmp < min ){    min = tmp;    sidx[tidx] = sidx[tidx ^ (     NWARP)];  }  smem[tidx].r = min;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= ( 4 * NWARP)
  tmp = smem[tidx ^ ( 2 * NWARP)].r;  if( tmp < min ){    min = tmp;    sidx[tidx] = sidx[tidx ^ ( 2 * NWARP)];  }  smem[tidx].r = min;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= ( 8 * NWARP)
  tmp = smem[tidx ^ ( 4 * NWARP)].r;  if( tmp < min ){    min = tmp;    sidx[tidx] = sidx[tidx ^ ( 4 * NWARP)];  }  smem[tidx].r = min;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= (16 * NWARP)
  tmp = smem[tidx ^ ( 8 * NWARP)].r;  if( tmp < min ){    min = tmp;    sidx[tidx] = sidx[tidx ^ ( 8 * NWARP)];  }  smem[tidx].r = min;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB == (32 * NWARP)
  tmp = smem[tidx ^ (16 * NWARP)].r;  if( tmp < min ){    min = tmp;    sidx[tidx] = sidx[tidx ^ (16 * NWARP)];  }  smem[tidx].r = min;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
  *loc = sidx[head];
  return (smem[head].r);
#endif///USE_WARP_SHUFFLE_FUNC
}


/**
 * @fn getMaxLocRealTsub
 *
 * @brief Get maximum value with location within a group of TSUB threads.
 * @detail implicit synchronization within TSUB (<= 32) threads is assumed
 * continuous NWARP threads have the same value as input
 */
#ifdef  USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ real getMaxLocRealTsub(const real max, int *loc)
#else///USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ real getMaxLocRealTsub(      real max, int *loc, volatile uint_real * smem, volatile int * sidx, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC
{
#ifdef  USE_WARP_SHUFFLE_FUNC
  real val =  max;
  int  idx = *loc;
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  int  buf;
  tmp = __SHFL_XOR(val,      NWARP, TSUB);  buf = __SHFL_XOR(idx,      NWARP, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= ( 4 * NWARP)
  tmp = __SHFL_XOR(val,  2 * NWARP, TSUB);  buf = __SHFL_XOR(idx,  2 * NWARP, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= ( 8 * NWARP)
  tmp = __SHFL_XOR(val,  4 * NWARP, TSUB);  buf = __SHFL_XOR(idx,  4 * NWARP, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= (16 * NWARP)
  tmp = __SHFL_XOR(val,  8 * NWARP, TSUB);  buf = __SHFL_XOR(idx,  8 * NWARP, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB == (32 * NWARP)
  tmp = __SHFL_XOR(val, 16 * NWARP, TSUB);  buf = __SHFL_XOR(idx, 16 * NWARP, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
  *loc =  __SHFL(idx, 0, TSUB);
  return (__SHFL(val, 0, TSUB));
#else///USE_WARP_SHUFFLE_FUNC
  smem[tidx].r =  max;
  sidx[tidx]   = *loc;
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = smem[tidx ^ (     NWARP)].r;  if( tmp > max ){    max = tmp;    sidx[tidx] = sidx[tidx ^ (     NWARP)];  }  smem[tidx].r = max;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= ( 4 * NWARP)
  tmp = smem[tidx ^ ( 2 * NWARP)].r;  if( tmp > max ){    max = tmp;    sidx[tidx] = sidx[tidx ^ ( 2 * NWARP)];  }  smem[tidx].r = max;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= ( 8 * NWARP)
  tmp = smem[tidx ^ ( 4 * NWARP)].r;  if( tmp > max ){    max = tmp;    sidx[tidx] = sidx[tidx ^ ( 4 * NWARP)];  }  smem[tidx].r = max;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB >= (16 * NWARP)
  tmp = smem[tidx ^ ( 8 * NWARP)].r;  if( tmp > max ){    max = tmp;    sidx[tidx] = sidx[tidx ^ ( 8 * NWARP)];  }  smem[tidx].r = max;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#   if  TSUB == (32 * NWARP)
  tmp = smem[tidx ^ (16 * NWARP)].r;  if( tmp > max ){    max = tmp;    sidx[tidx] = sidx[tidx ^ (16 * NWARP)];  }  smem[tidx].r = max;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
  *loc = sidx[head];
  return (smem[head].r);
#endif///USE_WARP_SHUFFLE_FUNC
}


#ifdef  ADOPT_APPROXIMATED_ENCLOSING_BALL
/**
 * @fn findFurthestParticle
 *
 * @brief Find the most distant particle.
 * @detail algorithm is based on Ritter (1990), ``An efficient bounding sphere'', Graphics Gems
 *
 * @sa getMaxLocRealTsub
 */
__device__ __forceinline__ real findFurthestParticle(const int lane, const position pi, const position pj, int *idx
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * RESTRICT smem, volatile int * RESTRICT sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
						     )
{
  *idx = lane;

  const real dx = pi.x - pj.x;
  const real dy = pi.y - pj.y;
  const real dz = pi.z - pj.z;

  const real r2max = getMaxLocRealTsub(SEB_TINY + dx * dx + dy * dy + dz * dz, idx
#ifndef USE_WARP_SHUFFLE_FUNC
				       , smem, sidx, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
				       );

  return (r2max);
}


/**
 * @fn approxSEB
 *
 * @brief Find the approximated smallest enclosing ball of a given set of data points.
 * @detail algorithm is based on Ritter (1990), ``An efficient bounding sphere'', Graphics Gems
 *
 * @sa findFurthestParticle
 */
__device__ __forceinline__ void approxSEB
(const int lane, jnode * RESTRICT ipos, const position pi, position * RESTRICT cen
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * RESTRICT smem, volatile int * RESTRICT sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
 )
{
  int idx = lane;
  real rold, d2max;

  {
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
    const real xmin = GET_MIN_TSUB_NWARP(pi.x                          );    const real xmax = GET_MAX_TSUB_NWARP(pi.x                          );
    const real ymin = GET_MIN_TSUB_NWARP(pi.y                          );    const real ymax = GET_MAX_TSUB_NWARP(pi.y                          );
    const real zmin = GET_MIN_TSUB_NWARP(pi.z                          );    const real zmax = GET_MAX_TSUB_NWARP(pi.z                          );
#else///USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
    const real xmin = GET_MIN_TSUB_NWARP(pi.x, (real *)smem, tidx, head);    const real xmax = GET_MAX_TSUB_NWARP(pi.x, (real *)smem, tidx, head);
    const real ymin = GET_MIN_TSUB_NWARP(pi.y, (real *)smem, tidx, head);    const real ymax = GET_MAX_TSUB_NWARP(pi.y, (real *)smem, tidx, head);
    const real zmin = GET_MIN_TSUB_NWARP(pi.z, (real *)smem, tidx, head);    const real zmax = GET_MAX_TSUB_NWARP(pi.z, (real *)smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC

    cen->x = HALF * (xmin + xmax);
    cen->y = HALF * (ymin + ymax);
    cen->z = HALF * (zmin + zmax);
    rold = SEB_TINY + HALF * FMAX(xmax - xmin, FMAX(ymax - ymin, zmax - zmin));
  }
  cen->m = rold * rold;

  while( (d2max = findFurthestParticle(lane, pi, *cen, &idx
#ifndef USE_WARP_SHUFFLE_FUNC
				      , smem, sidx, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
				       )) > cen->m ){
    /** update the sphere */
    const real dinv = RSQRT(d2max);
    const real dmax = dinv * d2max;

    /** calculate the displacement */
    const real rnew = HALF * (rold + dmax);
    const real disp = (dmax - rnew) * dinv;

#if 1
    /** work around to escape from the infinite loop */
    if ( disp < EPSILON )
      break;
#endif

    /** shift the center of the sphere */
    const position far = ipos[idx].pi;
    cen->x += (far.x - cen->x) * disp;
    cen->y += (far.y - cen->y) * disp;
    cen->z += (far.z - cen->z) * disp;

    /** enlarge the sphere */
    cen->m = rnew * rnew;
    rold = rnew;
    idx = lane;
  }
}
#endif//ADOPT_APPROXIMATED_ENCLOSING_BALL


#ifdef  ADOPT_SMALLEST_ENCLOSING_BALL
/** algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641 */


/** NOTE: translation from C code to CUDA code (to utilize real4 type) */
/** real qq[NDIM_SEB][NDIM_SEB] + real uu[NDIM_SEB] --> real4 qq[NDIM_SEB] */
/**      qq[INDEX2D(NDIM_SEB, NDIM_SEB, ii, jj)] = qq[ii][jj], qq[ii + NDIM2_SEB] = uu[ii]; */
/** real rr[NDIM_SEB][NDIM_SEB] + real ww[NDIM_SEB] --> real4 rr[NDIM_SEB] */
/**      rr[INDEX2D(NDIM_SEB, NDIM_SEB, ii, jj)] = rr[ii][jj], rr[ii + NDIM2_SEB] = ww[ii]; */
/** real cen[NDIM_SEB], r2 --> real4 cen */
/**      cen.x = cen[0], cen.y = cen[1], cen.z = cen[2], cen.w = r2; */
/** real cen2aff[NDIM_SEB], r2_to_aff --> real4 cen2aff */
/**      cen2aff.x = cen2aff[0], cen2aff.y = cen2aff[1], cen2aff.z = cen2aff[2], cen2aff.w = r2_to_aff; */
#define QQ(ii, jj) (qq[INDEX2D(NDIM_SEB, NDIM_SEB, ii, jj)])
#define RR(ii, jj) (rr[INDEX2D(NDIM_SEB, NDIM_SEB, ii, jj)])
#define UU(ii) (qq[NDIM2_SEB + (ii)])
#define WW(ii) (rr[NDIM2_SEB + (ii)])


/**
 * @fn initQR
 *
 * @brief Initialize support set and the corresponding matrices Q and R.
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 *
 * @param (lane) index of threads within a warp
 * @return (qq) on shared memory, orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim]
 * @return (rr) on shared memory, rectangular matrix R[Ndim][rank] and vector w[Ndim]
 * @param (idx) index of the specified data point as a support set
 * @return (seb) property of the smallest enclosing ball
 * @return (mem) indices of data points included in support set
 * @return (rank) rank of matrix A[rank][rank] = QR, rank = # of support points - 1
 */
__device__ __forceinline__ void initQR
(const int lane, volatile real * RESTRICT qq, volatile real * RESTRICT rr, const int idx, pos4seb * RESTRICT seb, volatile int * RESTRICT mem, int * RESTRICT rank)
{
  /** initialize Q to the identity matrix */
  /** initialize R to null */
  /** initialize uu and ww to null */
  if( lane < (NDIM2_SEB + NDIM_SEB) ){
    qq[lane] = ZERO;
    rr[lane] = ZERO;
  }/* if( lane < (NDIM2_SEB + NDIM_SEB) ){ */

  if( lane < NDIM_SEB )
    QQ(lane, lane) = UNITY;

  /** commit the inputted data point as a support set */
  if( lane == 0 )
    mem[0] = idx;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  *rank = 0;

  /** initialize the identifier of the support set */
  seb->support = ((lane != idx) ? false : true);
}


/**
 * @fn initSEB
 *
 * @brief Initialize the smallest enclosing ball (SEB).
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ void initSEB
(const int lane, jnode * RESTRICT pi, pos4seb * RESTRICT pos, real4 * RESTRICT cen, real * RESTRICT qq, real * RESTRICT rr, volatile int * RESTRICT mem, int * RESTRICT rank
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * RESTRICT smem, volatile int * RESTRICT sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
 )
{
  /** set a tentative center */
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
  const real xmin = GET_MIN_TSUB_NWARP(pos->x                          );  const real xmax = GET_MAX_TSUB_NWARP(pos->x                          );
  const real ymin = GET_MIN_TSUB_NWARP(pos->y                          );  const real ymax = GET_MAX_TSUB_NWARP(pos->y                          );
  const real zmin = GET_MIN_TSUB_NWARP(pos->z                          );  const real zmax = GET_MAX_TSUB_NWARP(pos->z                          );
#else///USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
  const real xmin = GET_MIN_TSUB_NWARP(pos->x, (real *)smem, tidx, head);  const real xmax = GET_MAX_TSUB_NWARP(pos->x, (real *)smem, tidx, head);
  const real ymin = GET_MIN_TSUB_NWARP(pos->y, (real *)smem, tidx, head);  const real ymax = GET_MAX_TSUB_NWARP(pos->y, (real *)smem, tidx, head);
  const real zmin = GET_MIN_TSUB_NWARP(pos->z, (real *)smem, tidx, head);  const real zmax = GET_MAX_TSUB_NWARP(pos->z, (real *)smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC

  cen->x = HALF * (xmin + xmax);
  cen->y = HALF * (ymin + ymax);
  cen->z = HALF * (zmin + zmax);

  /** find farthest point from the tentative center */
  const real dx = pos->x - cen->x;
  const real dy = pos->y - cen->y;
  const real dz = pos->z - cen->z;
  int idx = lane;
  cen->w = getMaxLocRealTsub(dx * dx + dy * dy + dz * dz, &idx
#ifndef USE_WARP_SHUFFLE_FUNC
			     , smem, sidx, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
			     );

  /** commit the farthest point as a support point */
  initQR(lane, qq, rr, idx, pos, mem, rank);
}


/**
 * @fn getShortestVector
 *
 * @brief Returns the square of the distance from the current center to aff(T).
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ real getShortestVector
(jnode * RESTRICT pi, const pos4seb pos, const real4 cen, real4 * RESTRICT cen2aff, volatile int * RESTRICT mem, const int rank, real * RESTRICT qq)
{
  /** compute vector from pos to origin */
  const position org = pi[mem[rank]].pi;
  cen2aff->x = org.x - cen.x;
  cen2aff->y = org.y - cen.y;
  cen2aff->z = org.z - cen.z;

  /** remove projections of ww onto the affine hull */
  for(int ii = 0; ii < rank; ii++){
    const real qx = QQ(ii, 0);
    const real qy = QQ(ii, 1);
    const real qz = QQ(ii, 2);

    const real scale = cen2aff->x * qx + cen2aff->y * qy + cen2aff->z * qz;

    cen2aff->x -= scale * qx;
    cen2aff->y -= scale * qy;
    cen2aff->z -= scale * qz;
  }/* for(int ii = 0; ii < rank; ii++){ */

  return (cen2aff->x * cen2aff->x + cen2aff->y * cen2aff->y + cen2aff->z * cen2aff->z);
}


/**
 * @fn findAffineCoefficients
 *
 * @brief Calculate affine coefficients of the current center with current support points.
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ void findAffineCoefficients
(const int lane, jnode * RESTRICT pi, const pos4seb pos, const real4 cen, volatile real * RESTRICT aff,
 const int rank, volatile int * RESTRICT mem, volatile real * RESTRICT qq, volatile real * RESTRICT rr)
{
  /** update uu and rr */
  if( lane < NDIM_SEB ){
    /** compute relative position of pos; i.e., uu = pos - origin */
    const position org = pi[mem[rank]].pi;
    if( lane == 0 ){
      UU(0) = cen.x - org.x;
      UU(1) = cen.y - org.y;
      UU(2) = cen.z - org.z;
    }/* if( lane == 0 ){ */

    /** calculate Q^{-1} uu into ww, Q is an orthogonal matrix */
    WW(lane) = QQ(lane, 0) * UU(0) + QQ(lane, 1) * UU(1) + QQ(lane, 2) * UU(2);
  }/* if( lane < NDIM_SEB ){ */
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700

  /** calculate lambda by backward substitution */
  /** TENTATIVE IMPLEMENTATION based on atomic function */
  real affRem = UNITY;
  for(int ii = rank - 1; ii >= 0; ii--){
    if( ((ii + 1) <= lane) && (lane < rank) )
      atomicAdd((real *)&WW(ii), -aff[lane] * RR(lane, ii));
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700

    const real tmp = WW(ii) / RR(ii, ii);
    affRem -= tmp;
    if( lane == 0 )
      aff[ii] = tmp;
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  }/* for(int ii = rank - 1; ii >= 0; ii--){ */

  if( lane == 0 )
    aff[rank] = affRem;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
}


/**
 * @fn givensRotation
 *
 * @brief Execute Givens rotation.
 * @detail algorithm is taken from Algorithm 4 in Anderson (2000)
 * COPYSIGN function is skipped in Fischer et al. (2003)
 */
__device__ __forceinline__ void givensRotation(const real ff, const real gg, real * RESTRICT cc, real * RESTRICT ss)
{
  if( gg == ZERO ){
    *cc = COPYSIGN(UNITY, ff);
    *ss = ZERO;
  }/* if( gg == ZERO ){ */
  else if( ff == ZERO ){
    *cc = ZERO;
    *ss = COPYSIGN(UNITY, gg);
  }/* else if( ff == ZERO ){ */
  else if( FABS(ff) >= FABS(gg) ){
    const real tt = gg / ff;
    *cc = COPYSIGN(RSQRT(UNITY + tt * tt), ff);
    *ss = (*cc) * tt;
  }/* else if( FABS(ff) >= FABS(gg) ){ */
  else{
    const real tt = ff / gg;
    *ss = COPYSIGN(RSQRT(UNITY + tt * tt), ff);
    *cc = (*ss) * tt;
  }/* else{ */
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700
}


/**
 * @fn clearHessenberg
 *
 * @brief Zero clear subdiagonal entries in R via Givens rotation
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ void clearHessenberg
(const int lane, const int rank, const int idx, volatile real * RESTRICT qq, volatile real * RESTRICT rr)
{
  for(int ii = idx; ii < rank; ii++){
    /** compute Givens coefficients */
    real cc, ss;
    givensRotation(RR(ii, ii), RR(ii, ii + 1), &cc, &ss);

    /** rotate R-rows */
    if( lane == ii )
      RR(ii, ii) = cc * RR(ii, ii) + ss * RR(ii, ii + 1);/**< the other one is an implicit zero */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700
    if( ((ii + 1) <= lane) && (lane < rank) ){
      const real ff = RR(lane, ii    );
      const real gg = RR(lane, ii + 1);

      RR(lane, ii    ) = cc * ff + ss * gg;
      RR(lane, ii + 1) = cc * gg - ss * ff;
    }/* if( ((ii + 1) <= lane) && (lane < rank) ){ */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700

    /** rotate Q-columns */
    if( lane < NDIM_SEB ){
      const real ff = QQ(ii    , lane);
      const real gg = QQ(ii + 1, lane);

      QQ(ii    , lane) = cc * ff + ss * gg;
      QQ(ii + 1, lane) = cc * gg - ss * ff;
    }/* if( lane < NDIM_SEB ){ */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  }/* for(int ii = idx; ii < qr->rank; ii++){ */
}


/**
 * @fn updateQR
 *
 * @brief Update matrices Q and R.
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ void updateQR(const int lane, const int rank, volatile real * RESTRICT qq, volatile real * RESTRICT rr)
{
  /** compute w = Q^{-1} u */
  if( lane < NDIM_SEB )
    WW(lane) = QQ(lane, 0) * UU(0) + QQ(lane, 1) * UU(1) + QQ(lane, 2) * UU(2);
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700

  /** rotate w */
  for(int jj = NDIM_SEB - 1; jj > 0; jj--){
    /** compute Givens coefficients */
    real cc, ss;
    givensRotation(WW(jj - 1), WW(jj), &cc, &ss);

    /** rotate w */
    if( lane == rank )
      WW(jj - 1) = cc * WW(jj - 1) + ss * WW(jj);
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700

    /** rotate two R-rows */
    if( lane == (jj - 1) ){
      RR(jj - 1, jj    )  = -ss * RR(jj - 1, jj - 1);
      RR(jj - 1, jj - 1) *=  cc;
    }/* if( lane == (jj - 1) ){ */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700
    if( (jj <= lane) && (lane < rank) ){
      const real ff = RR(lane, jj - 1);
      const real gg = RR(lane, jj    );

      RR(lane, jj - 1) = cc * ff + ss * gg;
      RR(lane, jj    ) = cc * gg - ss * ff;
    }/* if( (jj <= lane) && (lane < rank) ){ */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700

    /** rotate two Q-columns */
    if( lane < NDIM_SEB ){
      const real ff = QQ(jj - 1, lane);
      const real gg = QQ(jj    , lane);

      QQ(jj - 1, lane) = cc * ff + ss * gg;
      QQ(jj    , lane) = cc * gg - ss * ff;
    }/* if( lane < NDIM_SEB ){ */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  }/* for(int jj = NDIM_SEB - 1; jj > 0; jj--){ */

  /** update R */
  if( lane < rank )
    RR(lane, 0) += WW(0);
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700

  /** clear subdiagonal entries */
  clearHessenberg(lane, rank, 0, qq, rr);
}


/**
 * @fn removePoint
 *
 * @brief Remove the specified point from support set.
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ void removePoint
(const int lane, jnode * RESTRICT pi, const int idx, pos4seb * RESTRICT pos, volatile int * RESTRICT mem,
 int * RESTRICT rank, volatile real * RESTRICT qq, volatile real * RESTRICT rr)
{
  /** remove the point and update Q and R */
  if( lane == mem[idx] )
    pos->support = false;
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700

  if( idx == *rank ){
    /** remove origin if the remove point is the origin */
    if( lane == mem[(*rank) - 1] ){
      const position org = pi[mem[*rank]].pi;
      UU(0) = org.x - pos->x;
      UU(1) = org.y - pos->y;
      UU(2) = org.z - pos->z;
    }/* if( lane == mem.idx[*rank] ){ */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700

    (*rank)--;

    updateQR(lane, *rank, qq, rr);
  }/* if( idx == *rank ){ */
  else{
    /** general case: delete column vector from R */
    /** shift righter column of R one line to the left */
    const real tmp = ((lane < NDIM_SEB) ? (RR(idx, lane)) : (ZERO));
    for(int ii = idx + 1; ii < *rank; ii++){
      if( lane < NDIM_SEB )
	RR(ii - 1, lane) = RR(ii, lane);
      if( lane == NDIM_SEB )
	mem[ii - 1] = mem[ii];
#   if  __CUDA_ARCH__ >= 700
      __syncwarp();
#endif//__CUDA_ARCH__ >= 700
    }/* for(int ii = idx + 1; ii < *rank; ii++){ */

    if( lane == NDIM_SEB )
      mem[(*rank) - 1] = mem[*rank];
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700
    (*rank)--;

    if( lane < NDIM_SEB )
      RR(*rank, lane) = tmp;
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700

    clearHessenberg(lane, *rank, idx, qq, rr);
  }/* else{ */
}


/**
 * @fn dropSupport
 *
 * @brief Remove a data point from support set if it has non-positive affine coefficient.
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ bool dropSupport
(const int lane, jnode * RESTRICT pi, pos4seb * RESTRICT pos, const real4 cen, real * RESTRICT aff,
 int * RESTRICT rank, volatile int * RESTRICT mem, real * RESTRICT qq, real * RESTRICT rr)
{
  /** find affine coefficients of the center */
  findAffineCoefficients(lane, pi, *pos, cen, aff, *rank, mem, qq, rr);

  /** find a non-positive coefficient and drop it */
  int idx = 0;
  real min = UNITY;
  for(int ii = 0; ii < (*rank) + 1; ii++)
    if( aff[ii] < min ){
      min = aff[ii];
      idx = ii;
    }/* if( aff[ii] < min ){ */
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700

  /** drop a point with non-positive coefficient */
  if( min <= ZERO ){
    removePoint(lane, pi, idx, pos, mem, rank, qq, rr);
    return (true);
  }/* if( min <= ZERO ){ */

  return (false);
}


/**
 * @fn findStopper
 *
 * @brief Find a stopper against shrinking the SEB.
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ void findStopper
(const int lane, const pos4seb pos, const real4 cen, const real4 cen2aff, real * RESTRICT frac, int * RESTRICT idx
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * RESTRICT smem, volatile int * RESTRICT sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
 )
{
  const real tinyVal = SEB_EPS * cen.w * cen2aff.w * RSQRT(SEB_TINY + cen.w * cen2aff.w);
  *idx  = -1;
  real bound = UNITY;/**< this is the maximum case; i.e., there is no stopper */

  /** data points in support set are not candidate as a stopper */
  if( !pos.support ){
    const real dx = pos.x - cen.x;
    const real dy = pos.y - cen.y;
    const real dz = pos.z - cen.z;

    const real disp = dx * cen2aff.x + dy * cen2aff.y + dz * cen2aff.z;
    const real diff = cen2aff.w - disp;
    /** diff >= 0 is the necessary condition to be a stopper */
    if( diff >= tinyVal ){
      bound = HALF * (cen.w - dx * dx - dy * dy - dz * dz) * RSQRT(SEB_TINY + diff * diff);
      *idx  = lane;
    }/* if( diff >= tinyVal ){ */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  }/* if( !pos.support ){ */
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700

  *frac = getMinLocRealTsub(bound, idx
#ifndef USE_WARP_SHUFFLE_FUNC
			    , smem, sidx, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
			    );
}


/**
 * @fn appendColumn
 *
 * @brief Append a new column to R.
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ void appendColumn
(const int lane, const int rank, volatile real * RESTRICT qq, volatile real * RESTRICT rr)
{
  /** compute new column appendend to R */
  if( lane < NDIM_SEB )
    RR(rank, lane) = QQ(lane, 0) * UU(0) + QQ(lane, 1) * UU(1) + QQ(lane, 2) * UU(2);

  /** update QR-decomposition */
#pragma unroll
  for(int jj = NDIM_SEB - 1; jj > rank; jj--){
    real cc, ss;
    givensRotation(RR(rank, jj - 1), RR(rank, jj), &cc, &ss);

    /** rotate one R-entry (another entry is an implicit zero) */
    if( lane == 0 )
      RR(rank, jj - 1) = cc * RR(rank, jj - 1) + ss * RR(rank, jj);
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700

    /** rotate two Q-columns */
    if( lane < NDIM_SEB ){
      const real ff = QQ(jj - 1, lane);
      const real gg = QQ(jj    , lane);

      QQ(jj - 1, lane) = cc * ff + ss * gg;
      QQ(jj    , lane) = cc * gg - ss * ff;
    }/* if( lane < NDIM_SEB ){ */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  }/* for(int jj = Ndim - 1; jj > rank; jj--){ */
}


/**
 * @fn addSupport
 *
 * @brief Add the specified stopper to the support set.
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ void addSupport
(const int lane, jnode * RESTRICT pi, pos4seb * RESTRICT pos, const int idx,
 int * RESTRICT rank, volatile int * RESTRICT mem, volatile real * RESTRICT qq, real * RESTRICT rr
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * RESTRICT smem
#endif//USE_WARP_SHUFFLE_FUNC
 )
{
  /** compute u := pos[idx] - origin */
  if( lane == idx ){
    const position org = pi[mem[*rank]].pi;
    UU(0) = pos->x - org.x;
    UU(1) = pos->y - org.y;
    UU(2) = pos->z - org.z;
    pos->support = true;
    mem[(*rank) + 1] = mem[*rank];
    mem[ *rank     ] = idx;
  }/* if( lane == idx ){ */
#   if  __CUDA_ARCH__ >= 700
  __syncwarp();
#endif//__CUDA_ARCH__ >= 700

  /** append new column u to R and update QR-decomposition */
  appendColumn(lane, *rank, qq, rr);

  /** move origin index and insert new index */
  (*rank)++;
}


/**
 * @fn findSEB
 *
 * @brief Find the smallest enclosing ball of a given set of data points.
 * @detail algorithm is based on Fischer et al. (2003), ``Fast Smallest-Enclosing-Ball Computation in High Dimensions'', Proc. 11th European Symposium on Algorithms (ESA), 630-641
 */
__device__ __forceinline__ void findSEB
(const int lane, jnode * RESTRICT pi, pos4seb * RESTRICT pos, real4 * RESTRICT cen, real * RESTRICT qq, real * RESTRICT rr, volatile int * RESTRICT mem, real * RESTRICT aff
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * RESTRICT smem, volatile int * RESTRICT sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  /** define variables */
  real4 cen2aff;/**< vector from the current center to the circumcenter of aff(T) */
  int rank;

  initSEB(lane, pi, pos, cen, qq, rr, mem, &rank
#ifndef USE_WARP_SHUFFLE_FUNC
	  , smem, sidx, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
	  );

  while( true ){
    /** compute a walking direction */
    while( (cen2aff.w = getShortestVector(pi, *pos, *cen, &cen2aff, mem, rank, qq)) <= (SEB_EPS2 * cen->w) )
      if( !dropSupport(lane, pi, pos, *cen, aff, &rank, mem, qq, rr) )
	/** if dropSupport returns false, it means the center lies in the convex hull; i.e., SEB is already found */
	return;

    /** find the stopper in data points */
    real stopFrac;
    int  stopIdx;
    findStopper(lane, *pos, *cen, cen2aff, &stopFrac, &stopIdx
#ifndef USE_WARP_SHUFFLE_FUNC
		, smem, sidx, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
		);

    /** if a stopping point exists */
    if( (stopIdx >= 0) && ((1 + rank) <= NDIM_SEB) ){
      /** move the center */
      cen->x += stopFrac * cen2aff.x;
      cen->y += stopFrac * cen2aff.y;
      cen->z += stopFrac * cen2aff.z;

      /** update the radius */
      const position org = pi[mem[rank]].pi;
      const real dx = org.x - cen->x;
      const real dy = org.y - cen->y;
      const real dz = org.z - cen->z;
      cen->w = dx * dx + dy * dy + dz * dz;

      /** add the stopper to the support set */
      addSupport(lane, pi, pos, stopIdx, &rank, mem, qq, rr
#ifndef USE_WARP_SHUFFLE_FUNC
		 , smem
#endif//USE_WARP_SHUFFLE_FUNC
		 );
    }/* if( (stopIdx >= 0) && ((1 + rank) <= NDIM_SEB) ){ */

    /** if there is no stopping point */
    else{
      /** move the center */
      cen->x += cen2aff.x;
      cen->y += cen2aff.y;
      cen->z += cen2aff.z;

      /** update the radius */
      const position org = pi[mem[rank]].pi;
      const real dx = org.x - cen->x;
      const real dy = org.y - cen->y;
      const real dz = org.z - cen->z;
      cen->w = dx * dx + dy * dy + dz * dz;

      if( !dropSupport(lane, pi, pos, *cen, aff, &rank, mem, qq, rr) )
	return;
    }/* else{ */
#   if  __CUDA_ARCH__ >= 700
    __syncwarp();
#endif//__CUDA_ARCH__ >= 700
  }/* while( true ){ */
}
#endif//ADOPT_SMALLEST_ENCLOSING_BALL


#endif//SEB_DEV_CU
