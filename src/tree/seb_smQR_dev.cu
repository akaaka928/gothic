/*************************************************************************\
 *                                                                       *
                  last updated on 2016/04/13(Wed) 11:16:38
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
/*      qq[0].x = qq[0][0], qq[0].y = qq[0][1], qq[0].z = qq[0][2], qq[0].w = uu[0]; */
/*      qq[1].x = qq[1][0], qq[1].y = qq[1][1], qq[1].z = qq[1][2], qq[1].w = uu[1]; */
/*      qq[2].x = qq[2][0], qq[2].y = qq[2][1], qq[2].z = qq[2][2], qq[2].w = uu[2]; */
/* real rr[NDIM_SEB][NDIM_SEB] + real ww[NDIM_SEB] --> real4 rr[NDIM_SEB] */
/*      rr[0].x = rr[0][0], rr[0].y = rr[0][1], rr[0].z = rr[0][2], rr[0].w = ww[0]; */
/*      rr[1].x = rr[1][0], rr[1].y = rr[1][1], rr[1].z = rr[1][2], rr[1].w = ww[1]; */
/*      rr[2].x = rr[2][0], rr[2].y = rr[2][1], rr[2].z = rr[2][2], rr[2].w = ww[2]; */
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


//-------------------------------------------------------------------------
/* maximum value within a group of TSUB threads (TSUB <= 32 to use implicit synchronization) */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ real getSumRealTsub(const real  sum)
#else///USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ void getSumRealTsub(      real *sum, volatile uint_real * smem, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  real val = sum;
#   if  TSUB >=  2
  real tmp;
  tmp = __shfl_xor(val,  1, TSUB);  val += tmp;
#   if  TSUB >=  4
  tmp = __shfl_xor(val,  2, TSUB);  val += tmp;
#   if  TSUB >=  8
  tmp = __shfl_xor(val,  4, TSUB);  val += tmp;
#   if  TSUB >= 16
  tmp = __shfl_xor(val,  8, TSUB);  val += tmp;
#   if  TSUB == 32
  tmp = __shfl_xor(val, 16, TSUB);  val += tmp;
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
  return (__shfl(val, 0, TSUB));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx].r = *sum;
  //-----------------------------------------------------------------------
#   if  TSUB >=  2
  real tmp;
  tmp = smem[tidx ^  1].r;  smem[tidx].r += tmp;
#   if  TSUB >=  4
  tmp = smem[tidx ^  2].r;  smem[tidx].r += tmp;
#   if  TSUB >=  8
  tmp = smem[tidx ^  4].r;  smem[tidx].r += tmp;
#   if  TSUB >= 16
  tmp = smem[tidx ^  8].r;  smem[tidx].r += tmp;
#   if  TSUB == 32
  tmp = smem[tidx ^ 16].r;  smem[tidx].r += tmp;
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
  //-----------------------------------------------------------------------
  *sum = smem[head].r;
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* minimum value within a group of TSUB threads (TSUB <= 32 to use implicit synchronization) */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
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
  *loc =  __shfl(idx, 0, TSUB);
  return (__shfl(val, 0, TSUB));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx].r =  min;
  sidx[tidx]   = *loc;
  //-----------------------------------------------------------------------
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
  *loc =  __shfl(idx, 0, TSUB);
  return (__shfl(val, 0, TSUB));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx].r =  max;
  sidx[tidx]   = *loc;
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  *loc = sidx[head];
  return (smem[head].r);
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__device__ __forceinline__ real rdot (real4 aa, real4 bb){  return (aa.x * bb.x + aa.y * bb.y + aa.z * bb.z);}
__device__ __forceinline__ real rnrm2(real4 aa          ){  return (aa.x * aa.x + aa.y * aa.y + aa.z * aa.z);}
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
(const int lane, dat4seb * RESTRICT qq, dat4seb * RESTRICT rr, const int idx, pos4seb * RESTRICT seb, dat4seb * RESTRICT mem, int * RESTRICT rank)
{
  //-----------------------------------------------------------------------
  /* initialize Q to the identity matrix */
  /* initialize R to null */
  /* initialize uu and ww to null */
  //-----------------------------------------------------------------------
  if( lane < NDIM_SEB ){
    //---------------------------------------------------------------------
    const real4 r4zero = {ZERO, ZERO, ZERO, ZERO};
    rr[lane].r4 = r4zero;
    qq[lane].r4 = r4zero;
    if( lane == 0 )      rr[0].r4.x = UNITY;
    if( lane == 1 )      rr[1].r4.y = UNITY;
    if( lane == 2 )      rr[2].r4.z = UNITY;
    //---------------------------------------------------------------------
  }/* if( lane < NDIM_SEB ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit the inputted data point as a support set */
  mem->idx[0] = idx;
  *rank = 0;
  //-----------------------------------------------------------------------
  /* initialize the identifier of the support set */
  seb->support = false;
  if( lane == idx )
    seb->support = true;
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
(const int lane, pos4seb * RESTRICT pos, real4 * RESTRICT cen, dat4seb * RESTRICT qq, dat4seb * RESTRICT rr, dat4seb * RESTRICT mem, int * RESTRICT rank
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * smem, volatile int * sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  //-----------------------------------------------------------------------
  /* set a tentative center to the first point in pos */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  cen->x = __shfl(pos->x, 0, TSUB);
  cen->y = __shfl(pos->y, 0, TSUB);
  cen->z = __shfl(pos->z, 0, TSUB);
#else///USE_WARP_SHUFFLE_FUNC
  if( lane == 0 ){
    smem[0].r = pos->x;
    smem[1].r = pos->y;
    smem[2].r = pos->z;
  }/* if( lane == 0 ){ */
  cen->x = smem[0].r;
  cen->y = smem[1].r;
  cen->z = smem[2].r;
#endif//USE_WARP_SHUFFLE_FUNC
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
(const pos4seb pos, const real4 cen, real4 * RESTRICT cen2aff, const dat4seb mem, const int rank, dat4seb * RESTRICT qq
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * smem, const int lane
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  //-----------------------------------------------------------------------
  /* compute vector from pos to origin */
  real4 org;
  const int src = mem.idx[rank];
#ifdef  USE_WARP_SHUFFLE_FUNC
  org.x = __shfl(pos.x, src, TSUB);
  org.y = __shfl(pos.y, src, TSUB);
  org.z = __shfl(pos.z, src, TSUB);
#else///USE_WARP_SHUFFLE_FUNC
  if( lane == src ){
    smem[0].r = pos.x;
    smem[1].r = pos.y;
    smem[2].r = pos.z;
  }/* if( lane == src ){ */
  org.x = smem[0].r;
  org.y = smem[1].r;
  org.z = smem[2].r;
#endif//USE_WARP_SHUFFLE_FUNC
  cen2aff->x = org.x - cen.x;
  cen2aff->y = org.y - cen.y;
  cen2aff->z = org.z - cen.z;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* remove projections of ww onto the affine hull */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < rank; ii++){
    //---------------------------------------------------------------------
    const real4 mat = {qq[ii].r4.x, qq[ii].r4.y, qq[ii].r4.z, ZERO};
    //---------------------------------------------------------------------
    const real scale = rdot(*cen2aff, mat);
    //---------------------------------------------------------------------
    cen2aff->x -= scale * mat.x;
    cen2aff->y -= scale * mat.y;
    cen2aff->z -= scale * mat.z;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < rank; ii++){ */
  //-----------------------------------------------------------------------
  return (rnrm2(*cen2aff));
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
(const int lane, const pos4seb pos, const real4 cen, dat4seb * RESTRICT aff,
 const int rank, const dat4seb mem, dat4seb * RESTRICT qq, dat4seb * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  /* compute relative position of pos; i.e., uu = pos - origin */
  //-----------------------------------------------------------------------
  const int src = mem.idx[rank];
#ifdef  USE_WARP_SHUFFLE_FUNC
  real4 org;
  org.x = __shfl(pos.x, src, TSUB);  if( lane == 0 )    qq[lane].r4.w = cen.x - org.x;
  org.y = __shfl(pos.y, src, TSUB);  if( lane == 1 )    qq[lane].r4.w = cen.y - org.y;
  org.z = __shfl(pos.z, src, TSUB);  if( lane == 2 )    qq[lane].r4.w = cen.z - org.z;
#else///USE_WARP_SHUFFLE_FUNC
  if( lane == src ){
    smem[0].r = pos.x;
    smem[1].r = pos.y;
    smem[2].r = pos.z;
  }/* if( lane == src ){ */
  if( lane == 0 )    qq[lane].r4.w = cen.x - smem[lane].r;
  if( lane == 1 )    qq[lane].r4.w = cen.y - smem[lane].r;
  if( lane == 2 )    qq[lane].r4.w = cen.z - smem[lane].r;
#endif//USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate Q^{-1} uu into ww, Q is an orthogonal matrix */
  //-----------------------------------------------------------------------
  if( lane < NDIM_SEB )
    rr[lane].r4.w = qq[lane].r4.x * qq[0].r4.w + qq[lane].r4.y * qq[1].r4.w + qq[lane].r4.z * qq[2].r4.w;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate lambda by backward substitution */
  /* TENTATIVE IMPLEMENTATION based on atomic function */
  //-----------------------------------------------------------------------
  real affRem = UNITY;
#pragma unroll
  for(int ii = rank - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    if( ((ii + 1) <= lane) && (lane < rank) )
      atomicAdd(&(rr[ii].r4.w), -aff->val[lane] * rr[lane].val[ii]);
    //---------------------------------------------------------------------
    const real tmp = rr[ii].r4.w / rr[ii].val[ii];
    aff->val[ii] = tmp;
    affRem -= tmp;
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
__device__ __forceinline__ void clearHessenberg(const int lane, const int rank, const int idx, dat4seb * RESTRICT qq, dat4seb * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  for(int ii = idx; ii < rank; ii++){
    //---------------------------------------------------------------------
    /* compute Givens coefficients */
    real cc, ss;
    givensRotation(rr[ii].val[ii], rr[ii].val[ii + 1], &cc, &ss);
    //---------------------------------------------------------------------
    /* rotate R-rows */
    if( lane == rank )
      rr[ii].val[ii] = cc * rr[ii].val[ii] + ss * rr[ii].val[ii + 1];/* the other one is an implicit zero */
    if( ((ii + 1) <= lane) && (lane < rank) ){
      //-------------------------------------------------------------------
      const real ff = rr[lane].val[ii    ];
      const real gg = rr[lane].val[ii + 1];
      //-------------------------------------------------------------------
      rr[lane].val[ii    ] = cc * ff + ss * gg;
      rr[lane].val[ii + 1] = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* if( ((ii + 1) <= lane) && (lane < rank) ){ */
    //---------------------------------------------------------------------
    /* rotate Q-columns */
    if( lane < NDIM_SEB ){
      //-------------------------------------------------------------------
      const real ff = qq[ii    ].val[lane];
      const real gg = qq[ii + 1].val[lane];
      //-------------------------------------------------------------------
      qq[ii    ].val[lane] = cc * ff + ss * gg;
      qq[ii + 1].val[lane] = cc * gg - ss * ff;
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
__device__ __forceinline__ void updateQR(const int lane, const int rank, dat4seb * RESTRICT qq, dat4seb * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  /* compute w = Q^{-1} u */
  //-----------------------------------------------------------------------
  if( lane < NDIM_SEB )
    rr[lane].r4.w = qq[lane].r4.x * qq[0].r4.w + qq[lane].r4.y * qq[1].r4.w + qq[lane].r4.z * qq[2].r4.w;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* rotate w */
  //-----------------------------------------------------------------------
  for(int jj = NDIM_SEB - 1; jj > 0; jj--){
    //---------------------------------------------------------------------
    /* compute Givens coefficients */
    real cc, ss;
    givensRotation(rr[jj - 1].r4.w, rr[jj].r4.w, &cc, &ss);
    //---------------------------------------------------------------------
    /* rotate w */
    if( lane == rank )
      rr[jj - 1].r4.w = cc * rr[jj - 1].r4.w + ss * rr[jj].r4.w;
    //---------------------------------------------------------------------
    /* rotate two R-rows */
    if( lane == (jj - 1) ){
      rr[lane].val[jj    ]  = -ss * rr[lane].val[jj - 1];
      rr[lane].val[jj - 1] *=  cc;
    }/* if( lane == (jj - 1) ){ */
    if( lane < rank ){
      //-------------------------------------------------------------------
      const real ff = rr[lane].val[jj - 1];
      const real gg = rr[lane].val[jj    ];
      //-------------------------------------------------------------------
      rr[lane].val[jj - 1] = cc * ff + ss * gg;
      rr[lane].val[jj    ] = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* if( lane < rank ){ */
    //---------------------------------------------------------------------
    /* rotate two Q-columns */
    if( lane < NDIM_SEB ){
      //-------------------------------------------------------------------
      const real ff = qq[jj - 1].val[lane];
      const real gg = qq[jj    ].val[lane];
      //-------------------------------------------------------------------
      qq[jj - 1].val[lane] = cc * ff + ss * gg;
      qq[jj    ].val[lane] = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* if( lane < NDIM_SEB ){ */
    //---------------------------------------------------------------------
  }/* for(int jj = NDIM_SEB - 1; jj > 0; jj--){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* update R */
  //-----------------------------------------------------------------------
  if( lane < rank )
    rr[lane].r4.x += rr[0].r4.w;
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
(const int lane, const int idx, pos4seb * RESTRICT pos, dat4seb * RESTRICT mem,
 int * RESTRICT rank, dat4seb * RESTRICT qq, dat4seb * RESTRICT rr
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * smem
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  //-----------------------------------------------------------------------
  /* remove the point and update Q and R */
  //-----------------------------------------------------------------------
  if( lane == mem->idx[idx] )
    pos->support = false;
  //-----------------------------------------------------------------------
  if( idx == *rank ){
    //---------------------------------------------------------------------
    /* remove origin if the remove point is the origin */
    //---------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
    const int src = mem->idx[*rank];
    real4 org;
    org.x = __shfl(pos->x, src, TSUB);
    org.y = __shfl(pos->y, src, TSUB);
    org.z = __shfl(pos->z, src, TSUB);
    if( lane == mem->idx[(*rank) - 1] ){
      qq[0].r4.w = org.x - pos->x;
      qq[1].r4.w = org.y - pos->y;
      qq[2].r4.w = org.z - pos->z;
    }/* if( lane == mem.idx[*rank] ){ */
#else///USE_WARP_SHUFFLE_FUNC
    if( lane == mem.idx[*rank] ){
      smem[0].r = pos->x;
      smem[1].r = pos->y;
      smem[2].r = pos->z;
    }/* if( lane == mem.idx[*rank] ){ */
    if( lane == mem.idx[(*rank) - 1] ){
      qq[0].r4.w = smem[0].r - pos->x;
      qq[1].r4.w = smem[1].r - pos->y;
      qq[2].r4.w = smem[2].r - pos->z;
    }/* if( lane == mem.idx[*rank] ){ */
#endif//USE_WARP_SHUFFLE_FUNC
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
    dat4seb tmp = rr[idx];
    if( lane == 0 )
      for(int ii = idx + 1; ii < *rank; ii++){
	//-----------------------------------------------------------------
	rr[ii - 1].r4.x = rr[ii].r4.x;
	rr[ii - 1].r4.y = rr[ii].r4.y;
	rr[ii - 1].r4.z = rr[ii].r4.z;
	//-----------------------------------------------------------------
      }/* for(int ii = idx + 1; ii < *rank; ii++){ */
    //---------------------------------------------------------------------
#pragma unroll
    for(int ii = idx + 1; ii < *rank; ii++)
      mem->idx[ii - 1] = mem->idx[ii];
    mem->idx[(*rank) - 1] = mem->idx[*rank];
    //---------------------------------------------------------------------
    (*rank)--;
    if( lane == *rank ){
      rr[lane].r4.x = tmp.r4.x;
      rr[lane].r4.y = tmp.r4.y;
      rr[lane].r4.z = tmp.r4.z;
    }/* if( lane == *rank ){ */
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
(const int lane, pos4seb * RESTRICT pos, const real4 cen, dat4seb * RESTRICT aff,
 int * RESTRICT rank, dat4seb * RESTRICT mem, dat4seb * RESTRICT qq, dat4seb * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  /* find affine coefficients of the center */
  //-----------------------------------------------------------------------
  findAffineCoefficients(lane, *pos, cen, aff, *rank, *mem, qq, rr);
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
    removePoint(lane, idx, pos, mem, rank, qq, rr);
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
 , volatile uint_real * smem, volatile int * sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  //-----------------------------------------------------------------------
  *frac = UNITY;
  *idx  = -1;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const real tinyVal = SEB_EPS * cen.w * cen2aff.w * RSQRT(SEB_TINY + cen.w * cen2aff.w);
  //-----------------------------------------------------------------------
  real bound = UNITY;/* this is the maximum case; i.e., there is no stopper */
  /* data points in support set are not candidate as a stopper */
  if( !pos.support ){
    //---------------------------------------------------------------------
    real4 cen2rr;
    cen2rr.x = pos.x - cen.x;
    cen2rr.y = pos.y - cen.y;
    cen2rr.z = pos.z - cen.z;
    //---------------------------------------------------------------------
    const real disp = rdot(cen2rr, cen2aff);
    const real diff = cen2aff.w - disp;
    /* diff >= 0 is the necessary condition to be a stopper */
    if( diff >= tinyVal ){
      bound = HALF * (cen.w - rnrm2(cen2rr)) * RSQRT(SEB_TINY + diff * diff);
      *idx  = lane;
    }
    //---------------------------------------------------------------------
  }/* if( !pos.support ){ */
  //-----------------------------------------------------------------------
  *frac = getMinLocRealTsub(bound, idx
#ifndef USE_WARP_SHUFFLE_FUNC
			    , smem, sidx, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
			    );
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
(const int lane, const int rank, dat4seb * RESTRICT qq, dat4seb * RESTRICT rr)
{
  //-----------------------------------------------------------------------
  /* compute new column appendend to R */
  //-----------------------------------------------------------------------
  if( lane < NDIM_SEB )
    rr[rank].val[lane] = qq[lane].r4.x * qq[0].r4.w + qq[lane].r4.y * qq[1].r4.w + qq[lane].r4.z * qq[2].r4.w;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* update QR-decomposition */
  //-----------------------------------------------------------------------
#pragma unroll
  for(int jj = NDIM_SEB - 1; jj > rank; jj--){
    //---------------------------------------------------------------------
    real cc, ss;
    givensRotation(rr[rank].val[jj - 1], rr[rank].val[jj], &cc, &ss);
    //---------------------------------------------------------------------
    /* rotate one R-entry (another entry is an implicit zero) */
    if( lane == rank )
      rr[lane].val[jj - 1] = cc * rr[lane].val[jj - 1] + ss * rr[lane].val[jj];
    //---------------------------------------------------------------------
    /* rotate two Q-columns */
    if( lane < NDIM_SEB ){
      //-------------------------------------------------------------------
      const real ff = qq[jj - 1].val[lane];
      const real gg = qq[jj    ].val[lane];
      //-------------------------------------------------------------------
      qq[jj - 1].val[lane] = cc * ff + ss * gg;
      qq[jj    ].val[lane] = cc * gg - ss * ff;
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
(const int lane, pos4seb * RESTRICT pos, const int idx,
 int * RESTRICT rank, dat4seb * RESTRICT mem, dat4seb * RESTRICT qq, dat4seb * RESTRICT rr
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * smem
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  //-----------------------------------------------------------------------
  /* compute u := pos[idx] - origin */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  const int src = mem->idx[*rank];
  real4 org;
  org.x = __shfl(pos->x, src, TSUB);
  org.y = __shfl(pos->y, src, TSUB);
  org.z = __shfl(pos->z, src, TSUB);
  if( lane == idx ){
    qq[0].r4.w = pos->x - org.x;
    qq[1].r4.w = pos->y - org.y;
    qq[2].r4.w = pos->z - org.z;
    pos->support = true;
  }/* if( lane == idx ){ */
#else///USE_WARP_SHUFFLE_FUNC
  if( lane == mem.idx[*rank] ){
    smem[0].r = pos->x;
    smem[1].r = pos->y;
    smem[2].r = pos->z;
  }/* if( lane == mem.idx[*rank] ){ */
  if( lane == idx ){
    qq[0].r4.w = pos->x - smem[0].r;
    qq[1].r4.w = pos->x - smem[1].r;
    qq[2].r4.w = pos->x - smem[2].r;
    pos->support = true;
  }/* if( lane == idx ){ */
#endif//USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* append new column u to R and update QR-decomposition */
  //-----------------------------------------------------------------------
  appendColumn(lane, *rank, qq, rr);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* move origin index and insert new index */
  //-----------------------------------------------------------------------
  mem->idx[(*rank) + 1] = mem->idx[*rank];
  mem->idx[ *rank     ] = idx;
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
(const int lane, pos4seb * RESTRICT pos, real4 * RESTRICT cen, dat4seb * RESTRICT qq, dat4seb * RESTRICT rr
#ifndef USE_WARP_SHUFFLE_FUNC
 , volatile uint_real * smem, volatile int * sidx, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
)
{
  //-----------------------------------------------------------------------
  /* define variables */
  //-----------------------------------------------------------------------
  real4 cen2aff;/* vector from the current center to the circumcenter of aff(T) */
  dat4seb aff;/* affine coefficents, maximum # of support points is d + 1 */
  dat4seb mem;/* indices of points in support set */
  //-----------------------------------------------------------------------
  int rank;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  initSEB(lane, pos, cen, qq, rr, &mem, &rank);
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    /* compute a walking direction */
    //---------------------------------------------------------------------
    while( (cen2aff.w = getShortestVector(*pos, *cen, &cen2aff, mem, rank, qq)) <= (SEB_EPS2 * cen->w) )
      if( !dropSupport(lane, pos, *cen, &aff, &rank, &mem, qq, rr) )
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
#ifdef  USE_WARP_SHUFFLE_FUNC
      const int src = mem.idx[rank];
      const real dx = __shfl(pos->x, src, TSUB) - cen->x;
      const real dy = __shfl(pos->y, src, TSUB) - cen->y;
      const real dz = __shfl(pos->z, src, TSUB) - cen->z;
#else///USE_WARP_SHUFFLE_FUNC
      if( lane == mem.idx[rank] ){
	smem[0].r = pos->x;
	smem[1].r = pos->y;
	smem[2].r = pos->z;
      }/* if( lane == mem.idx[rank] ){ */
      const real dx = smem[0].r - cen->x;
      const real dy = smem[0].r - cen->y;
      const real dz = smem[0].r - cen->z;
#endif//USE_WARP_SHUFFLE_FUNC
      cen->w = dx * dx + dy * dy + dz * dz;
      //-------------------------------------------------------------------
      /* add the stopper to the support set */
      addSupport(lane, pos, stopIdx, &rank, &mem, qq, rr);
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
#ifdef  USE_WARP_SHUFFLE_FUNC
      const int src = mem.idx[rank];
      const real dx = __shfl(pos->x, src, TSUB) - cen->x;
      const real dy = __shfl(pos->y, src, TSUB) - cen->y;
      const real dz = __shfl(pos->z, src, TSUB) - cen->z;
#else///USE_WARP_SHUFFLE_FUNC
      if( lane == mem.idx[rank] ){
	smem[0].r = pos->x;
	smem[1].r = pos->y;
	smem[2].r = pos->z;
      }/* if( lane == mem.idx[rank] ){ */
      const real dx = smem[0].r - cen->x;
      const real dy = smem[0].r - cen->y;
      const real dz = smem[0].r - cen->z;
#endif//USE_WARP_SHUFFLE_FUNC
      cen->w = dx * dx + dy * dy + dz * dz;
      //-------------------------------------------------------------------
      if( !dropSupport(lane, pos, *cen, &aff, &rank, &mem, qq, rr) )
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
