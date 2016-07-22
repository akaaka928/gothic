/*************************************************************************\
 *                                                                       *
                  last updated on 2016/04/14(Thu) 16:27:54
 *                                                                       *
 *    Generate Smallest Enclosing Ball                                   *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef SEB_DEV_CU
#define SEB_DEV_CU
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
#include "walk_dev.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* NOTE: translation from C code to CUDA code (to utilize real4 type) */
//-------------------------------------------------------------------------
/* real qq[NDIM_SEB][NDIM_SEB] + real uu[NDIM_SEB] --> real4 qq[NDIM_SEB] */
/*      qq[INDEX2D(NDIM_SEB, NDIM_SEB, ii, jj)] = qq[ii][jj], qq[ii + NDIM2_SEB] = uu[ii]; */
/* real rr[NDIM_SEB][NDIM_SEB] + real ww[NDIM_SEB] --> real4 rr[NDIM_SEB] */
/*      rr[INDEX2D(NDIM_SEB, NDIM_SEB, ii, jj)] = rr[ii][jj], rr[ii + NDIM2_SEB] = ww[ii]; */
/* real cen[NDIM_SEB], r2 --> real4 cen */
/*      cen.x = cen[0], cen.y = cen[1], cen.z = cen[2], cen.w = r2; */
/* real cen2aff[NDIM_SEB], r2_to_aff --> real4 cen2aff */
/*      cen2aff.x = cen2aff[0], cen2aff.y = cen2aff[1], cen2aff.z = cen2aff[2], cen2aff.w = r2_to_aff; */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  DOUBLE_PRECISION
#       define SEB_EPS  (1.0e-15)
#       define SEB_EPS2 (1.0e-30)
#       define SEB_TINY (DBL_MIN)
#else///DOUBLE_PRECISION
#       define SEB_EPS  (1.0e-6f)
#       define SEB_EPS2 (1.0e-12f)
#       define SEB_TINY (FLT_MIN)
#endif//DOUBLE_PRECISION
//-------------------------------------------------------------------------
#define QQ(ii, jj) (qq[INDEX2D(NDIM_SEB, NDIM_SEB, ii, jj)])
#define RR(ii, jj) (rr[INDEX2D(NDIM_SEB, NDIM_SEB, ii, jj)])
#define UU(ii) (qq[NDIM2_SEB + (ii)])
#define WW(ii) (rr[NDIM2_SEB + (ii)])
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* minimum value within a group of TSUB threads (TSUB <= 32 to use implicit synchronization) */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
/* continuous NWARP threads have the same value as input */
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ real getMinLocRealTsub(const real min, int *loc)
#else///USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ real getMinLocRealTsub(const real min, int *loc, volatile uint_real * smem, volatile int * sidx, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  real val =  min;
  int  idx = *loc;
#if 1
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  int  buf;
  tmp = __shfl_xor(val,      NWARP, TSUB);  buf = __shfl_xor(idx,      NWARP, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  TSUB >= ( 4 * NWARP)
  tmp = __shfl_xor(val,  2 * NWARP, TSUB);  buf = __shfl_xor(idx,  2 * NWARP, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  TSUB >= ( 8 * NWARP)
  tmp = __shfl_xor(val,  4 * NWARP, TSUB);  buf = __shfl_xor(idx,  4 * NWARP, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  TSUB >= (16 * NWARP)
  tmp = __shfl_xor(val,  8 * NWARP, TSUB);  buf = __shfl_xor(idx,  8 * NWARP, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  TSUB == (32 * NWARP)
  tmp = __shfl_xor(val, 16 * NWARP, TSUB);  buf = __shfl_xor(idx, 16 * NWARP, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
#else
#   if  TSUB >=  2
  real tmp;
  int  buf;
  tmp = __shfl_xor(val,  1, TSUB);  buf = __shfl_xor(idx,  1, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  TSUB >=  4
  tmp = __shfl_xor(val,  2, TSUB);  buf = __shfl_xor(idx,  2, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  TSUB >=  8
  tmp = __shfl_xor(val,  4, TSUB);  buf = __shfl_xor(idx,  4, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  TSUB >= 16
  tmp = __shfl_xor(val,  8, TSUB);  buf = __shfl_xor(idx,  8, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#   if  TSUB == 32
  tmp = __shfl_xor(val, 16, TSUB);  buf = __shfl_xor(idx, 16, TSUB);  if( tmp < val ){    val = tmp;    idx = buf;  }
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
#endif
  *loc =  __shfl(idx, 0, TSUB);
  return (__shfl(val, 0, TSUB));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx].r =  min;
  sidx[tidx]   = *loc;
  //-----------------------------------------------------------------------
#if 1
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = smem[tidx ^ (     NWARP)].r;  if( tmp < *min ){    *min = tmp;    sidx[tidx] = sidx[tidx ^ (     NWARP)];  }  smem[tidx].r = *min;
#   if  TSUB >= ( 4 * NWARP)
  tmp = smem[tidx ^ ( 2 * NWARP)].r;  if( tmp < *min ){    *min = tmp;    sidx[tidx] = sidx[tidx ^ ( 2 * NWARP)];  }  smem[tidx].r = *min;
#   if  TSUB >= ( 8 * NWARP)
  tmp = smem[tidx ^ ( 4 * NWARP)].r;  if( tmp < *min ){    *min = tmp;    sidx[tidx] = sidx[tidx ^ ( 4 * NWARP)];  }  smem[tidx].r = *min;
#   if  TSUB >= (16 * NWARP)
  tmp = smem[tidx ^ ( 8 * NWARP)].r;  if( tmp < *min ){    *min = tmp;    sidx[tidx] = sidx[tidx ^ ( 8 * NWARP)];  }  smem[tidx].r = *min;
#   if  TSUB == (32 * NWARP)
  tmp = smem[tidx ^ (16 * NWARP)].r;  if( tmp < *min ){    *min = tmp;    sidx[tidx] = sidx[tidx ^ (16 * NWARP)];  }  smem[tidx].r = *min;
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
#else
#   if  TSUB >=  2
  real tmp;
  tmp = smem[tidx ^  1].r;  if( tmp < *min ){    *min = tmp;    sidx[tidx] = sidx[tidx ^  1];  }  smem[tidx].r = *min;
#   if  TSUB >=  4
  tmp = smem[tidx ^  2].r;  if( tmp < *min ){    *min = tmp;    sidx[tidx] = sidx[tidx ^  2];  }  smem[tidx].r = *min;
#   if  TSUB >=  8
  tmp = smem[tidx ^  4].r;  if( tmp < *min ){    *min = tmp;    sidx[tidx] = sidx[tidx ^  4];  }  smem[tidx].r = *min;
#   if  TSUB >= 16
  tmp = smem[tidx ^  8].r;  if( tmp < *min ){    *min = tmp;    sidx[tidx] = sidx[tidx ^  8];  }  smem[tidx].r = *min;
#   if  TSUB == 32
  tmp = smem[tidx ^ 16].r;  if( tmp < *min ){    *min = tmp;    sidx[tidx] = sidx[tidx ^ 16];  }  smem[tidx].r = *min;
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
#endif
  //-----------------------------------------------------------------------
  *loc = sidx[head];
  return (smem[head].r);
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* maximum value within a group of TSUB threads (TSUB <= 32 to use implicit synchronization) */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
/* continuous NWARP threads have the same value as input */
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ real getMaxLocRealTsub(const real max, int *loc)
#else///USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ real getMaxLocRealTsub(const real max, int *loc, volatile uint_real * smem, volatile int * sidx, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  real val =  max;
  int  idx = *loc;
#if 1
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  int  buf;
  tmp = __shfl_xor(val,      NWARP, TSUB);  buf = __shfl_xor(idx,      NWARP, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  TSUB >= ( 4 * NWARP)
  tmp = __shfl_xor(val,  2 * NWARP, TSUB);  buf = __shfl_xor(idx,  2 * NWARP, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  TSUB >= ( 8 * NWARP)
  tmp = __shfl_xor(val,  4 * NWARP, TSUB);  buf = __shfl_xor(idx,  4 * NWARP, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  TSUB >= (16 * NWARP)
  tmp = __shfl_xor(val,  8 * NWARP, TSUB);  buf = __shfl_xor(idx,  8 * NWARP, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  TSUB == (32 * NWARP)
  tmp = __shfl_xor(val, 16 * NWARP, TSUB);  buf = __shfl_xor(idx, 16 * NWARP, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
#else
#   if  TSUB >=  2
  real tmp;
  int  buf;
  tmp = __shfl_xor(val,  1, TSUB);  buf = __shfl_xor(idx,  1, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  TSUB >=  4
  tmp = __shfl_xor(val,  2, TSUB);  buf = __shfl_xor(idx,  2, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  TSUB >=  8
  tmp = __shfl_xor(val,  4, TSUB);  buf = __shfl_xor(idx,  4, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  TSUB >= 16
  tmp = __shfl_xor(val,  8, TSUB);  buf = __shfl_xor(idx,  8, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#   if  TSUB == 32
  tmp = __shfl_xor(val, 16, TSUB);  buf = __shfl_xor(idx, 16, TSUB);  if( tmp > val ){    val = tmp;    idx = buf;  }
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
#endif
  *loc =  __shfl(idx, 0, TSUB);
  return (__shfl(val, 0, TSUB));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx].r =  max;
  sidx[tidx]   = *loc;
  //-----------------------------------------------------------------------
#if 1
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = smem[tidx ^ (     NWARP)].r;  if( tmp > *max ){    *max = tmp;    sidx[tidx] = sidx[tidx ^ (     NWARP)];  }  smem[tidx].r = *max;
#   if  TSUB >= ( 4 * NWARP)
  tmp = smem[tidx ^ ( 2 * NWARP)].r;  if( tmp > *max ){    *max = tmp;    sidx[tidx] = sidx[tidx ^ ( 2 * NWARP)];  }  smem[tidx].r = *max;
#   if  TSUB >= ( 8 * NWARP)
  tmp = smem[tidx ^ ( 4 * NWARP)].r;  if( tmp > *max ){    *max = tmp;    sidx[tidx] = sidx[tidx ^ ( 4 * NWARP)];  }  smem[tidx].r = *max;
#   if  TSUB >= (16 * NWARP)
  tmp = smem[tidx ^ ( 8 * NWARP)].r;  if( tmp > *max ){    *max = tmp;    sidx[tidx] = sidx[tidx ^ ( 8 * NWARP)];  }  smem[tidx].r = *max;
#   if  TSUB == (32 * NWARP)
  tmp = smem[tidx ^ (16 * NWARP)].r;  if( tmp > *max ){    *max = tmp;    sidx[tidx] = sidx[tidx ^ (16 * NWARP)];  }  smem[tidx].r = *max;
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
#else
#   if  TSUB >=  2
  real tmp;
  tmp = smem[tidx ^  1].r;  if( tmp > *max ){    *max = tmp;    sidx[tidx] = sidx[tidx ^  1];  }  smem[tidx].r = *max;
#   if  TSUB >=  4
  tmp = smem[tidx ^  2].r;  if( tmp > *max ){    *max = tmp;    sidx[tidx] = sidx[tidx ^  2];  }  smem[tidx].r = *max;
#   if  TSUB >=  8
  tmp = smem[tidx ^  4].r;  if( tmp > *max ){    *max = tmp;    sidx[tidx] = sidx[tidx ^  4];  }  smem[tidx].r = *max;
#   if  TSUB >= 16
  tmp = smem[tidx ^  8].r;  if( tmp > *max ){    *max = tmp;    sidx[tidx] = sidx[tidx ^  8];  }  smem[tidx].r = *max;
#   if  TSUB == 32
  tmp = smem[tidx ^ 16].r;  if( tmp > *max ){    *max = tmp;    sidx[tidx] = sidx[tidx ^ 16];  }  smem[tidx].r = *max;
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
#endif
  //-----------------------------------------------------------------------
  *loc = sidx[head];
  return (smem[head].r);
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* initialize support set and the corresponding matrices Q and R */
//-------------------------------------------------------------------------
/* lane ::  input ::               :: index of threads within a warp */
/* qq   :: output :: shared memory :: orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim] */
/* rr   :: output :: shared memory :: rectangular matrix R[Ndim][rank] and vector w[Ndim] */
/* idx  ::  input ::               :: index of the specified data point as a support set */
/* mem  :: output ::               :: indices of data points included in support set */
/* rank :: output ::               :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
//-------------------------------------------------------------------------
__device__ __forceinline__ void initQR
(const int lane, volatile real * RESTRICT qq, volatile real * RESTRICT rr, const int idx, pos4seb * RESTRICT seb, volatile int * RESTRICT mem, int * RESTRICT rank)
{
  //-----------------------------------------------------------------------
  /* initialize Q to the identity matrix */
  /* initialize R to null */
  /* initialize uu and ww to null */
  //-----------------------------------------------------------------------
  if( lane < (NDIM2_SEB + NDIM_SEB) ){
    //---------------------------------------------------------------------
    qq[lane] = ZERO;
    rr[lane] = ZERO;
    //---------------------------------------------------------------------
  }/* if( lane < (NDIM2_SEB + NDIM_SEB) ){ */
  //-----------------------------------------------------------------------
  if( lane < NDIM_SEB )
    QQ(lane, lane) = UNITY;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit the inputted data point as a support set */
  if( lane == 0 )
    mem[0] = idx;
  *rank = 0;
  //-----------------------------------------------------------------------
  /* initialize the identifier of the support set */
  seb->support = ((lane != idx) ? false : true);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* initialize the smallest enclosing ball (SEB) */
//-------------------------------------------------------------------------
/* lane ::  input ::               :: index of threads within a warp */
/* qq   :: output :: shared memory :: orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim] */
/* rr   :: output :: shared memory :: rectangular matrix R[Ndim][rank] and vector w[Ndim] */
/* dat  ::  input ::               :: position of data points for (lane)-th data */
/* cen  :: output ::               :: position of the SEB's center and squared radius of that */
/* mem  :: output ::               :: indices of data points included in support set */
/* rank :: output ::               :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
//-------------------------------------------------------------------------
__device__ __forceinline__ void initSEB
(const int lane, jnode * RESTRICT pi, pos4seb * RESTRICT pos, real4 * RESTRICT cen, real * RESTRICT qq, real * RESTRICT rr, int * RESTRICT mem, int * RESTRICT rank
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * RESTRICT smem, volatile int * RESTRICT sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  //-----------------------------------------------------------------------
  /* set a tentative center to the first point in pos */
  //-----------------------------------------------------------------------
  const position ipos = pi[0].pi;
  cen->x = ipos.x;
  cen->y = ipos.y;
  cen->z = ipos.z;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* find farthest point from the tentative center */
  //-----------------------------------------------------------------------
  const real dx = pos->x - cen->x;
  const real dy = pos->y - cen->y;
  const real dz = pos->z - cen->z;
  int idx = lane;
  cen->w = getMaxLocRealTsub(dx * dx + dy * dy + dz * dz, &idx
#ifndef USE_WARP_SHUFFLE_FUNC
			     , smem, sidx, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
			     );
  //-----------------------------------------------------------------------
#if 0
  if( lane == THREADIDX_X1D )
    printf("lane = %d:\tidx = %d, cen->w = %e, loc = %e\n", lane, idx, cen->w, dx * dx + dy * dy + dz * dz);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit the farthest point as a support point */
  //-----------------------------------------------------------------------
  initQR(lane, qq, rr, idx, pos, mem, rank);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* returns the square of the distance from the current center to aff(T) */
//-------------------------------------------------------------------------
/* dat     ::  input :: position of data points for (lane)-th data */
/* cen     ::  input :: position of the SEB's center and squared radius of that */
/* cen2aff :: output :: vector from the current center to the circumcenter of aff(T) */
/* mem     ::  input :: indices of data points included in support set */
/* rank    ::  input :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* qq      ::  input :: orthogonal matrix Q[Ndim][Ndim] */
//-------------------------------------------------------------------------
__device__ __forceinline__ real getShortestVector
(jnode * RESTRICT pi, const pos4seb pos, const real4 cen, real4 * RESTRICT cen2aff, int * RESTRICT mem, const int rank, real * RESTRICT qq)
{
  //-----------------------------------------------------------------------
  /* compute vector from pos to origin */
  const position org = pi[mem[rank]].pi;
  cen2aff->x = org.x - cen.x;
  cen2aff->y = org.y - cen.y;
  cen2aff->z = org.z - cen.z;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* remove projections of ww onto the affine hull */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < rank; ii++){
    //---------------------------------------------------------------------
    const real qx = QQ(ii, 0);
    const real qy = QQ(ii, 1);
    const real qz = QQ(ii, 2);
    //---------------------------------------------------------------------
    const real scale = cen2aff->x * qx + cen2aff->y * qy + cen2aff->z * qz;
    //---------------------------------------------------------------------
    cen2aff->x -= scale * qx;
    cen2aff->y -= scale * qy;
    cen2aff->z -= scale * qz;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < rank; ii++){ */
  //-----------------------------------------------------------------------
  return (cen2aff->x * cen2aff->x + cen2aff->y * cen2aff->y + cen2aff->z * cen2aff->z);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* calculate affine coefficients of the current center with current support points */
//-------------------------------------------------------------------------
/* lane ::  input ::               :: index of threads within a warp */
/* dat  ::  input ::               :: position of data points for (lane)-th data */
/* cen  ::  input ::               :: position of the SEB's center and squared radius of that */
/* aff  :: output ::               :: affine coefficents lambda[rank + 1] */
/* rank ::  input ::               :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* mem  ::  input ::               :: indices of data points included in support set */
/* qq   :: output :: shared memory :: orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim] */
/* rr   :: output :: shared memory :: rectangular matrix R[Ndim][rank] and vector w[Ndim] */
//-------------------------------------------------------------------------
__device__ __forceinline__ void findAffineCoefficients
(const int lane, jnode * RESTRICT pi, const pos4seb pos, const real4 cen, dat4seb * RESTRICT aff,
 const int rank, int * RESTRICT mem, volatile real * RESTRICT qq, volatile real * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  /* update uu and rr */
  //-----------------------------------------------------------------------
  if( lane < NDIM_SEB ){
    //---------------------------------------------------------------------
    /* compute relative position of pos; i.e., uu = pos - origin */
    if( lane == 0 ){
      const position org = pi[mem[rank]].pi;
      UU(0) = cen.x - org.x;
      UU(1) = cen.y - org.y;
      UU(2) = cen.z - org.z;
    }/* if( lane == 0 ){ */
    //---------------------------------------------------------------------
    /* calculate Q^{-1} uu into ww, Q is an orthogonal matrix */
    WW(lane) = QQ(lane, 0) * UU(0) + QQ(lane, 1) * UU(1) + QQ(lane, 2) * UU(2);
    //---------------------------------------------------------------------
  }/* if( lane < NDIM_SEB ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate lambda by backward substitution */
  /* TENTATIVE IMPLEMENTATION based on atomic function */
  //-----------------------------------------------------------------------
  real affRem = UNITY;
  for(int ii = rank - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    if( ((ii + 1) <= lane) && (lane < rank) )
      atomicAdd((real *)&WW(ii), -aff->val[lane] * RR(lane, ii));
    //---------------------------------------------------------------------
    const real tmp = WW(ii) / RR(ii, ii);
    aff->val[ii]  = tmp;
    affRem       -= tmp;
    //---------------------------------------------------------------------
  }/* for(int ii = rank - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  aff->val[rank] = affRem;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* algorithm is taken from Algorithm 4 in Anderson (2000) */
/* COPYSIGN function is skipped in Fischer et al. (2003) */
//-------------------------------------------------------------------------
__device__ __forceinline__ void givensRotation(const real ff, const real gg, real * RESTRICT cc, real * RESTRICT ss)
{
  //-----------------------------------------------------------------------
  if( gg == ZERO ){
    //---------------------------------------------------------------------
    *cc = COPYSIGN(UNITY, ff);
    *ss = ZERO;
    //---------------------------------------------------------------------
  }/* if( gg == ZERO ){ */
  //-----------------------------------------------------------------------
  else if( ff == ZERO ){
    //---------------------------------------------------------------------
    *cc = ZERO;
    *ss = COPYSIGN(UNITY, gg);
    //---------------------------------------------------------------------
  }/* else if( ff == ZERO ){ */
  //-----------------------------------------------------------------------
  else if( FABS(ff) >= FABS(gg) ){
    //---------------------------------------------------------------------
    const real tt = gg / ff;
    *cc = COPYSIGN(RSQRT(UNITY + tt * tt), ff);
    *ss = (*cc) * tt;
    //---------------------------------------------------------------------
  }/* else if( FABS(ff) >= FABS(gg) ){ */
  //-----------------------------------------------------------------------
  else{
    //---------------------------------------------------------------------
    const real tt = ff / gg;
    *ss = COPYSIGN(RSQRT(UNITY + tt * tt), ff);
    *cc = (*ss) * tt;
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* zero clear subdiagonal entries in R via Givens rotation */
//-------------------------------------------------------------------------
/* lane ::  input ::               :: index of threads within a warp */
/* rank ::  input ::               :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* idx  ::  input ::               :: index of the head column which has non-zero subdiagonal element */
/* qq   :: output :: shared memory :: orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim] */
/* rr   :: output :: shared memory :: rectangular matrix R[Ndim][rank] and vector w[Ndim] */
//-------------------------------------------------------------------------
__device__ __forceinline__ void clearHessenberg
(const int lane, const int rank, const int idx, volatile real * RESTRICT qq, volatile real * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  for(int ii = idx; ii < rank; ii++){
    //---------------------------------------------------------------------
    /* compute Givens coefficients */
    real cc, ss;
    givensRotation(RR(ii, ii), RR(ii, ii + 1), &cc, &ss);
    //---------------------------------------------------------------------
    /* rotate R-rows */
    if( lane == ii )
      RR(ii, ii) = cc * RR(ii, ii) + ss * RR(ii, ii + 1);/* the other one is an implicit zero */
    if( ((ii + 1) <= lane) && (lane < rank) ){
      //-------------------------------------------------------------------
      const real ff = RR(lane, ii    );
      const real gg = RR(lane, ii + 1);
      //-------------------------------------------------------------------
      RR(lane, ii    ) = cc * ff + ss * gg;
      RR(lane, ii + 1) = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* if( ((ii + 1) <= lane) && (lane < rank) ){ */
    //---------------------------------------------------------------------
    /* rotate Q-columns */
    if( lane < NDIM_SEB ){
      //-------------------------------------------------------------------
      const real ff = QQ(ii    , lane);
      const real gg = QQ(ii + 1, lane);
      //-------------------------------------------------------------------
      QQ(ii    , lane) = cc * ff + ss * gg;
      QQ(ii + 1, lane) = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* if( lane < NDIM_SEB ){ */
    //---------------------------------------------------------------------
  }/* for(int ii = idx; ii < qr->rank; ii++){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* update matrices Q and R */
//-------------------------------------------------------------------------
/* lane ::  input ::               :: index of threads within a warp */
/* rank ::  input ::               :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* qq   :: output :: shared memory :: orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim] */
/* rr   :: output :: shared memory :: rectangular matrix R[Ndim][rank] and vector w[Ndim] */
//-------------------------------------------------------------------------
__device__ __forceinline__ void updateQR(const int lane, const int rank, volatile real * RESTRICT qq, volatile real * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  /* compute w = Q^{-1} u */
  //-----------------------------------------------------------------------
  if( lane < NDIM_SEB )
    WW(lane) = QQ(lane, 0) * UU(0) + QQ(lane, 1) * UU(1) + QQ(lane, 2) * UU(2);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* rotate w */
  //-----------------------------------------------------------------------
  for(int jj = NDIM_SEB - 1; jj > 0; jj--){
    //---------------------------------------------------------------------
    /* compute Givens coefficients */
    real cc, ss;
    givensRotation(WW(jj - 1), WW(jj), &cc, &ss);
    //---------------------------------------------------------------------
    /* rotate w */
    if( lane == rank )
      WW(jj - 1) = cc * WW(jj - 1) + ss * WW(jj);
    //---------------------------------------------------------------------
    /* rotate two R-rows */
    if( lane == (jj - 1) ){
      RR(jj - 1, jj    )  = -ss * RR(jj - 1, jj - 1);
      RR(jj - 1, jj - 1) *=  cc;
    }/* if( lane == (jj - 1) ){ */
    if( (jj <= lane) && (lane < rank) ){
      //-------------------------------------------------------------------
      const real ff = RR(lane, jj - 1);
      const real gg = RR(lane, jj    );
      //-------------------------------------------------------------------
      RR(lane, jj - 1) = cc * ff + ss * gg;
      RR(lane, jj    ) = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* if( (jj <= lane) && (lane < rank) ){ */
    //---------------------------------------------------------------------
    /* rotate two Q-columns */
    if( lane < NDIM_SEB ){
      //-------------------------------------------------------------------
      const real ff = QQ(jj - 1, lane);
      const real gg = QQ(jj    , lane);
      //-------------------------------------------------------------------
      QQ(jj - 1, lane) = cc * ff + ss * gg;
      QQ(jj    , lane) = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* if( lane < NDIM_SEB ){ */
    //---------------------------------------------------------------------
  }/* for(int jj = NDIM_SEB - 1; jj > 0; jj--){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* update R */
  //-----------------------------------------------------------------------
  if( lane < rank )
    RR(lane, 0) += WW(0);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* clear subdiagonal entries */
  //-----------------------------------------------------------------------
  clearHessenberg(lane, rank, 0, qq, rr);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* remove the specified point from support set */
//-------------------------------------------------------------------------
/* lane ::  input ::               :: index of threads within a warp */
/* idx  ::  input ::               :: index of the support point to be removed */
/* dat  :: output ::               :: position of data points for (lane)-th data */
/* mem  :: output ::               :: indices of data points included in support set */
/* rank :: output ::               :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* qq   :: output :: shared memory :: orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim] */
/* rr   :: output :: shared memory :: rectangular matrix R[Ndim][rank] and vector w[Ndim] */
//-------------------------------------------------------------------------
__device__ __forceinline__ void removePoint
(const int lane, jnode * RESTRICT pi, const int idx, pos4seb * RESTRICT pos, volatile int * RESTRICT mem,
 int * RESTRICT rank, volatile real * RESTRICT qq, volatile real * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  /* remove the point and update Q and R */
  //-----------------------------------------------------------------------
  if( lane == mem[idx] )
    pos->support = false;
  //-----------------------------------------------------------------------
  if( idx == *rank ){
    //---------------------------------------------------------------------
    /* remove origin if the remove point is the origin */
    //---------------------------------------------------------------------
    if( lane == mem[(*rank) - 1] ){
      const position org = pi[mem[*rank]].pi;
      UU(0) = org.x - pos->x;
      UU(1) = org.y - pos->y;
      UU(2) = org.z - pos->z;
    }/* if( lane == mem.idx[*rank] ){ */
    //---------------------------------------------------------------------
    (*rank)--;
    //---------------------------------------------------------------------
    updateQR(lane, *rank, qq, rr);
    //---------------------------------------------------------------------
  }/* if( idx == *rank ){ */
  //-----------------------------------------------------------------------
  else{
    //---------------------------------------------------------------------
    /* general case: delte column vector from R */
    //---------------------------------------------------------------------
    /* shift righter column of R one line to the left */
    const real tmp = ((lane < NDIM_SEB) ? (RR(idx, lane)) : (ZERO));
#pragma unroll
    for(int ii = idx + 1; ii < *rank; ii++){
      //-------------------------------------------------------------------
      if( lane < NDIM_SEB )
	RR(ii - 1, lane) = RR(ii, lane);
      if( lane == NDIM_SEB )
	mem[ii - 1] = mem[ii];
      //-------------------------------------------------------------------
    }/* for(int ii = idx + 1; ii < *rank; ii++){ */
    //---------------------------------------------------------------------
    if( lane == NDIM_SEB )
      mem[(*rank) - 1] = mem[*rank];
    (*rank)--;
    //---------------------------------------------------------------------
    if( lane < NDIM_SEB )
      RR(*rank, lane) = tmp;
    //---------------------------------------------------------------------
    clearHessenberg(lane, *rank, idx, qq, rr);
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* remove a data point from support set if it has non-positive affine coefficient */
//-------------------------------------------------------------------------
/* lane ::  input ::               :: index of threads within a warp */
/* dat  :: output ::               :: position of data points for (lane)-th data */
/* cen  ::  input ::               :: position of the SEB's center and squared radius of that */
/* aff  :: output ::               :: affine coefficents lambda[rank + 1] */
/* rank :: output ::               :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* mem  :: output ::               :: indices of data points included in support set */
/* qq   :: output :: shared memory :: orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim] */
/* rr   :: output :: shared memory :: rectangular matrix R[Ndim][rank] and vector w[Ndim] */
//-------------------------------------------------------------------------
__device__ __forceinline__ bool dropSupport
(const int lane, jnode * RESTRICT pi, pos4seb * RESTRICT pos, const real4 cen, dat4seb * RESTRICT aff,
 int * RESTRICT rank, int * RESTRICT mem, real * RESTRICT qq, real * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  /* find affine coefficients of the center */
  //-----------------------------------------------------------------------
  findAffineCoefficients(lane, pi, *pos, cen, aff, *rank, mem, qq, rr);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* find a non-positive coefficient and drop it */
  //-----------------------------------------------------------------------
  int idx = 0;
  real min = UNITY;
#pragma unroll
  for(int ii = 0; ii < (*rank) + 1; ii++)
    if( aff->idx[ii] < min ){
      min = aff->val[ii];
      idx = ii;
    }/* if( aff[ii] < min ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* drop a point with non-positive coefficient */
  //-----------------------------------------------------------------------
  if( min <= ZERO ){
    //---------------------------------------------------------------------
    removePoint(lane, pi, idx, pos, mem, rank, qq, rr);
    return (true);
    //---------------------------------------------------------------------
  }/* if( min <= ZERO ){ */
  //-----------------------------------------------------------------------
  return (false);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* find a stopper against shrinking the SEB */
//-------------------------------------------------------------------------
/* lane    ::  input :: index of threads within a warp */
/* dat     ::  input :: position of data points for (lane)-th data */
/* cen     ::  input :: position of the SEB's center and squared radius of that */
/* cen2aff ::  input :: vector from the current center to the circumcenter of aff(T) */
/* frac    :: output :: fraction of final radius of the SEB against the initial radius of that */
/* idx     :: output :: index of the stopper */
//-------------------------------------------------------------------------
__device__ __forceinline__ void findStopper
(const int lane, const pos4seb pos, const real4 cen, const real4 cen2aff, real * RESTRICT frac, int * RESTRICT idx
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * RESTRICT smem, volatile int * RESTRICT sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  //-----------------------------------------------------------------------
  const real tinyVal = SEB_EPS * cen.w * cen2aff.w * RSQRT(SEB_TINY + cen.w * cen2aff.w);
  *idx  = -1;
  //-----------------------------------------------------------------------
  real bound = UNITY;/* this is the maximum case; i.e., there is no stopper */
  /* data points in support set are not candidate as a stopper */
  if( !pos.support ){
    //---------------------------------------------------------------------
    const real dx = pos.x - cen.x;
    const real dy = pos.y - cen.y;
    const real dz = pos.z - cen.z;
    //---------------------------------------------------------------------
    const real disp = dx * cen2aff.x + dy * cen2aff.y + dz * cen2aff.z;
    const real diff = cen2aff.w - disp;
    /* diff >= 0 is the necessary condition to be a stopper */
    if( diff >= tinyVal ){
      bound = HALF * (cen.w - dx * dx - dy * dy - dz * dz) * RSQRT(SEB_TINY + diff * diff);
      *idx  = lane;
    }/* if( diff >= tinyVal ){ */
    //---------------------------------------------------------------------
  }/* if( !pos.support ){ */
  //-----------------------------------------------------------------------
  *frac = getMinLocRealTsub(bound, idx
#ifndef USE_WARP_SHUFFLE_FUNC
			    , smem, sidx, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
			    );
  //-----------------------------------------------------------------------
#if 0
  if( lane == GLOBALIDX_X1D )
    printf("lane = %d:\tidx = %d, frac = %e, bound = %e\n", lane, *idx, *frac, bound);
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* append a new column to R */
//-------------------------------------------------------------------------
/* lane ::  input ::               :: index of threads within a warp */
/* rank ::  input ::               :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* qq   :: output :: shared memory :: orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim] */
/* rr   :: output :: shared memory :: rectangular matrix R[Ndim][rank] and vector w[Ndim] */
//-------------------------------------------------------------------------
__device__ __forceinline__ void appendColumn
(const int lane, const int rank, volatile real * RESTRICT qq, volatile real * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  /* compute new column appendend to R */
  //-----------------------------------------------------------------------
  if( lane < NDIM_SEB )
    RR(rank, lane) = QQ(lane, 0) * UU(0) + QQ(lane, 1) * UU(1) + QQ(lane, 2) * UU(2);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* update QR-decomposition */
  //-----------------------------------------------------------------------
#pragma unroll
  for(int jj = NDIM_SEB - 1; jj > rank; jj--){
    //---------------------------------------------------------------------
    real cc, ss;
    givensRotation(RR(rank, jj - 1), RR(rank, jj), &cc, &ss);
    //---------------------------------------------------------------------
    /* rotate one R-entry (another entry is an implicit zero) */
    if( lane == 0 )
      RR(rank, jj - 1) = cc * RR(rank, jj - 1) + ss * RR(rank, jj);
    //---------------------------------------------------------------------
    /* rotate two Q-columns */
    if( lane < NDIM_SEB ){
      //-------------------------------------------------------------------
      const real ff = QQ(jj - 1, lane);
      const real gg = QQ(jj    , lane);
      //-------------------------------------------------------------------
      QQ(jj - 1, lane) = cc * ff + ss * gg;
      QQ(jj    , lane) = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* if( lane < NDIM_SEB ){ */
    //---------------------------------------------------------------------
  }/* for(int jj = Ndim - 1; jj > rank; jj--){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* add the specified stopper to the support set */
//-------------------------------------------------------------------------
/* lane ::  input ::               :: index of threads within a warp */
/* pos  :: output ::               :: position of data points for (lane)-th data */
/* idx  ::  input ::               :: index of the specified data point as a support set */
/* rank :: output ::               :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* mem  :: output ::               :: indices of data points included in support set */
/* qq   :: output :: shared memory :: orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim] */
/* rr   :: output :: shared memory :: rectangular matrix R[Ndim][rank] and vector w[Ndim] */
//-------------------------------------------------------------------------
__device__ __forceinline__ void addSupport
(const int lane, jnode * RESTRICT pi, pos4seb * RESTRICT pos, const int idx,
 int * RESTRICT rank, volatile int * RESTRICT mem, volatile real * RESTRICT qq, real * RESTRICT rr
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * RESTRICT smem
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  //-----------------------------------------------------------------------
  /* compute u := pos[idx] - origin */
  //-----------------------------------------------------------------------
  if( lane == idx ){
    const position org = pi[mem[*rank]].pi;
    UU(0) = pos->x - org.x;
    UU(1) = pos->y - org.y;
    UU(2) = pos->z - org.z;
    pos->support = true;
  }/* if( lane == idx ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* append new column u to R and update QR-decomposition */
  //-----------------------------------------------------------------------
  appendColumn(lane, *rank, qq, rr);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* move origin index and insert new index */
  //-----------------------------------------------------------------------
  if( lane == 0 ){
    mem[(*rank) + 1] = mem[*rank];
    mem[ *rank     ] = idx;
  }/* if( lane == 0 ){ */
  //-----------------------------------------------------------------------
  (*rank)++;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* find the smallest enclosing ball of a given set of data points */
//-------------------------------------------------------------------------
/* lane ::  input ::               :: index of threads within a warp */
/* pos  ::  input ::               :: position of data points for (lane)-th data */
/* cen  :: output ::               :: position of the SEB's center and squared radius of that */
/* qq   :: output :: shared memory :: orthogonal matrix Q[Ndim][Ndim] and vector u[Ndim] */
/* rr   :: output :: shared memory :: rectangular matrix R[Ndim][rank] and vector w[Ndim] */
//-------------------------------------------------------------------------
__device__ __forceinline__ void findSEB
(const int lane, jnode * RESTRICT pi, pos4seb * RESTRICT pos, real4 * RESTRICT cen, real * RESTRICT qq, real * RESTRICT rr, int * RESTRICT mem
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * RESTRICT smem, volatile int * RESTRICT sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  //-----------------------------------------------------------------------
  /* define variables */
  //-----------------------------------------------------------------------
  real4 cen2aff;/* vector from the current center to the circumcenter of aff(T) */
  dat4seb aff;/* affine coefficents, maximum # of support points is d + 1 */
  /* dat4seb mem;/\* indices of points in support set *\/ */
  //-----------------------------------------------------------------------
  int rank;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  initSEB(lane, pi, pos, cen, qq, rr, mem, &rank);
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
#if 0
    if( GLOBALIDX_X1D == 28224 ){
      printf("rank = %d:\tsupport = %d", rank, mem.idx[0]);
      for(int ii = 1; ii < rank + 1; ii++)
	printf(", %d", mem.idx[ii]);
      printf("\n");
    }
#endif
    //---------------------------------------------------------------------
    /* compute a walking direction */
    //---------------------------------------------------------------------
    while( (cen2aff.w = getShortestVector(pi, *pos, *cen, &cen2aff, mem, rank, qq)) <= (SEB_EPS2 * cen->w) )
      if( !dropSupport(lane, pi, pos, *cen, &aff, &rank, mem, qq, rr) )
	/* if dropSupport returns false, it means the center lies in the convex hull; i.e., SEB is already found */
	return;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* find the stopper in data points */
    //---------------------------------------------------------------------
    real stopFrac;
    int  stopIdx;
    findStopper(lane, *pos, *cen, cen2aff, &stopFrac, &stopIdx
#ifndef USE_WARP_SHUFFLE_FUNC
		, smem, sidx, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
		);
    //---------------------------------------------------------------------
    /* if a stopping point exists */
    if( (stopIdx >= 0) && ((1 + rank) <= NDIM_SEB) ){
      //-------------------------------------------------------------------
      /* move the center */
      cen->x += stopFrac * cen2aff.x;
      cen->y += stopFrac * cen2aff.y;
      cen->z += stopFrac * cen2aff.z;
      //-------------------------------------------------------------------
      /* update the radius */
      const position org = pi[mem[rank]].pi;
      const real dx = org.x - cen->x;
      const real dy = org.y - cen->y;
      const real dz = org.z - cen->z;
      cen->w = dx * dx + dy * dy + dz * dz;
      //-------------------------------------------------------------------
      /* add the stopper to the support set */
      addSupport(lane, pi, pos, stopIdx, &rank, mem, qq, rr);
      //-------------------------------------------------------------------
    }/* if( (stopIdx >= 0) && ((1 + rank) <= NDIM_SEB) ){ */
    //---------------------------------------------------------------------
    /* if there is no stopping point */
    else{
      //-------------------------------------------------------------------
      /* move the center */
      cen->x += cen2aff.x;
      cen->y += cen2aff.y;
      cen->z += cen2aff.z;
      //-------------------------------------------------------------------
      /* update the radius */
      const position org = pi[mem[rank]].pi;
      const real dx = org.x - cen->x;
      const real dy = org.y - cen->y;
      const real dz = org.z - cen->z;
      cen->w = dx * dx + dy * dy + dz * dz;
      //-------------------------------------------------------------------
      if( !dropSupport(lane, pi, pos, *cen, &aff, &rank, mem, qq, rr) )
	return;
      //-------------------------------------------------------------------
    }/* else{ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//SEB_DEV_CU
//-------------------------------------------------------------------------
