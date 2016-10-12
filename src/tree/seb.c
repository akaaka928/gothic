/*************************************************************************\
 *                                                                       *
                  last updated on 2016/10/11(Tue) 17:08:43
 *                                                                       *
 *    Test Code to generate Smallest Enclosing Ball                      *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#define NDIM_SEB (3)
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//-------------------------------------------------------------------------
#include <gsl/gsl_rng.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <myutil.h>
//-------------------------------------------------------------------------
gsl_rng *GSLRand;
#define UNIRAND ((real)gsl_rng_uniform(GSLRand))
#define RANDVAL (TWO * (UNIRAND) - UNITY)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  DOUBLE_PRECISION
static const double sebEps  = 1.0e-15;
static const double sebEps2 = 1.0e-30;
static const double sebTiny = DBL_MIN;
#else///DOUBLE_PRECISION
static const  float sebEps  = 1.0e-6f;
static const  float sebEps2 = 1.0e-12f;
static const  float sebTiny = FLT_MIN;
#endif//DOUBLE_PRECISION
//-------------------------------------------------------------------------
#define SEB_AFFINE_ORIGIN(Ndim, mem, rank, jj) (pos[INDEX2D(Ndat, Ndim, mem[rank], jj)])
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline real rdot(const int Ndim, real * restrict aa, real * restrict bb)
{
  //-----------------------------------------------------------------------
  real sum = ZERO;
  //-----------------------------------------------------------------------
  for(int jj = 0; jj < Ndim; jj++)
    sum += aa[jj] * bb[jj];
  //-----------------------------------------------------------------------
  return (sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline real rnrm2(const int Ndim, real * restrict vec)
{
  //-----------------------------------------------------------------------
  real sum = ZERO;
  //-----------------------------------------------------------------------
  for(int jj = 0; jj < Ndim; jj++)
    sum += vec[jj] * vec[jj];
  //-----------------------------------------------------------------------
  return (sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* initialize support set and the corresponding matrices Q and R */
//-------------------------------------------------------------------------
/* Ndim ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* Ndat ::  input :: Number of data points */
/* sup  :: output :: true for data points in support set; false otherwise */
/* qq   :: output :: orthogonal matrix Q[Ndim][Ndim] */
/* rr   :: output :: rectangular matrix R[Ndim][rank] */
/* uu   :: output :: vector u[Ndim] to store the relative position pos[Ndim] - origin[Ndim] */
/* ww   :: output :: vector w[Ndim] to store the matrix vector multiplication Qu */
/* idx  ::  input :: index of the specified data point as a support set */
/* mem  :: output :: indices of data points included in support set */
/* rank :: output :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
//-------------------------------------------------------------------------
void initQR
(const int Ndim, const int Ndat, bool * restrict sup,
 real qq[restrict][NDIM_SEB], real rr[restrict][NDIM_SEB], real uu[restrict], real ww[restrict],
 const int idx, int mem[restrict], int * restrict rank)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start ");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* initialize Q to the identity matrix */
  for(int ii = 0; ii < Ndim; ii++)
    for(int jj = 0; jj < Ndim; jj++)
      qq[ii][jj] = (jj != ii) ? (ZERO) : (UNITY);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ndim; ii++)
    for(int jj = 0; jj < Ndim; jj++)
      rr[ii][jj] = ZERO;
  for(int ii = 0; ii < Ndim; ii++)
    uu[ii] = ww[ii] = ZERO;
  //-----------------------------------------------------------------------
  /* commit the inputted data point as a support set */
  mem[0] = idx;
  *rank = 0;
  //-----------------------------------------------------------------------
  /* initialize the identifier of the support set */
  for(int ii = 0; ii < Ndat; ii++)
    sup[ii] = false;
  sup[idx] = true;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* initialize the smallest enclosing ball (SEB) */
//-------------------------------------------------------------------------
/* Ndim ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* Ndat ::  input :: Number of data points */
/* pos  ::  input :: position of data points p[Ndat * Ndim] */
/* cen  :: output :: position of the SEB's center */
/* r2   :: output :: square of the radius of the SEB */
/* sup  :: output :: true for data points in support set; false otherwise */
/* qq   :: output :: orthogonal matrix Q[Ndim][Ndim] */
/* rr   :: output :: rectangular matrix R[Ndim][rank] */
/* uu   :: output :: vector u[Ndim] to store the relative position pos[Ndim] - origin[Ndim] */
/* ww   :: output :: vector w[Ndim] to store the matrix vector multiplication Qu */
/* mem  :: output :: indices of data points included in support set */
/* rank :: output :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
//-------------------------------------------------------------------------
void initSEB
(const int Ndim, const int Ndat, real * restrict pos, real cen[restrict], real * restrict r2,
 bool * restrict sup, real qq[restrict][NDIM_SEB], real rr[restrict][NDIM_SEB], real uu[restrict], real ww[restrict], int mem[restrict], int * restrict rank)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set a tentative center to the first point in pos */
  //-----------------------------------------------------------------------
  for(int jj = 0; jj < Ndim; jj++)
    cen[jj] = pos[INDEX2D(Ndat, Ndim, 0, jj)];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* find farthest point from the tentative center */
  //-----------------------------------------------------------------------
  real r2max = ZERO;
  int idx = 0;
  for(int ii = 1; ii < Ndat; ii++){
    //---------------------------------------------------------------------
    real r2 = ZERO;
    for(int jj = 0; jj < Ndim; jj++){
      const real dr = pos[INDEX2D(Ndat, Ndim, ii, jj)] - cen[jj];
      r2 += dr * dr;
    }/* for(int jj = 0; jj < Ndim; jj++){ */
    //---------------------------------------------------------------------
    if( r2 >= r2max ){
      r2max = r2;
      idx = ii;
    }/* if( r2 >= r2max ){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 1; ii < Ndat; ii++){ */
  //-----------------------------------------------------------------------
  *r2 = r2max;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit the farthest point as a support point */
  //-----------------------------------------------------------------------
  initQR(Ndim, Ndat, sup, qq, rr, uu, ww, idx, mem, rank);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* returns the square of the distance from the current center to aff(T) */
//-------------------------------------------------------------------------
/* Ndim    ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* pos     ::  input :: position of data points p[Ndat * Ndim] */
/* cen     ::  input :: position of the SEB's center */
/* cen2aff :: output :: vector from the current center to the circumcenter of aff(T) */
/* mem     ::  input :: indices of data points included in support set */
/* rank    ::  input :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* qq      ::  input :: orthogonal matrix Q[Ndim][Ndim] */
//-------------------------------------------------------------------------
static inline real getShortestVector
(const int Ndim, real * restrict pos, real cen[restrict], real cen2aff[restrict], int mem[restrict], const int rank, real qq[restrict][NDIM_SEB])
{
  //-----------------------------------------------------------------------
  /* compute vector from pos to origin */
  for(int jj = 0; jj < Ndim; jj++)
    cen2aff[jj] = SEB_AFFINE_ORIGIN(Ndim, mem, rank, jj) - cen[jj];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* remove projections of ww onto the affine hull */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < rank; ii++){
    //---------------------------------------------------------------------
    const real scale = rdot(Ndim, cen2aff, qq[ii]);
    //---------------------------------------------------------------------
    for(int jj = 0; jj < Ndim; jj++)
      cen2aff[jj] -= scale * qq[ii][jj];
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < rank; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  return (rnrm2(Ndim, cen2aff));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* calculate affine coefficients of the current center with current support points */
//-------------------------------------------------------------------------
/* Ndim ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* pos  ::  input :: position of data points p[Ndat * Ndim] */
/* cen  ::  input :: position of the SEB's center */
/* aff  :: output :: affine coefficents lambda[rank + 1] */
/* rank ::  input :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* mem  ::  input :: indices of data points included in support set */
/* uu   :: output :: vector u[Ndim] to store the relative position pos[Ndim] - origin[Ndim] */
/* ww   :: output :: vector w[Ndim] to store the matrix vector multiplication Qu */
/* qq   ::  input :: orthogonal matrix Q[Ndim][Ndim] */
/* rr   ::  input :: rectangular matrix R[Ndim][rank] */
//-------------------------------------------------------------------------
static inline void findAffineCoefficients
(const int Ndim, real * restrict pos, real cen[restrict], real aff[restrict],
 const int rank, int mem[restrict], real uu[restrict], real ww[restrict], real qq[restrict][NDIM_SEB], real rr[restrict][NDIM_SEB])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* compute relative position of pos; i.e., uu = pos - origin */
  //-----------------------------------------------------------------------
  for(int jj = 0; jj < Ndim; jj++)
    uu[jj] = cen[jj] - SEB_AFFINE_ORIGIN(Ndim, mem, rank, jj);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate Q^{-1} uu into ww, Q is an orthogonal matrix */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ndim; ii++){
    ww[ii] = ZERO;
    for(int jj = 0; jj < Ndim; jj++)
      ww[ii] += qq[ii][jj] * uu[jj];
  }/* for(int ii = 0; ii < Ndim; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate lambda by backward substitution */
  //-----------------------------------------------------------------------
  real affRem = UNITY;
  for(int ii = rank - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    for(int jj = ii + 1; jj < rank; jj++)
      ww[ii] -= aff[jj] * rr[jj][ii];
    //---------------------------------------------------------------------
    const real tmp = ww[ii] / rr[ii][ii];
    aff[ii] = tmp;
    affRem -= tmp;
    //---------------------------------------------------------------------
  }/* for(int ii = rank - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  aff[rank] = affRem;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* algorithm is taken from Algorithm 4 in Anderson (2000) */
/* COPYSIGN function is skipped in Fischer et al. (2003) */
//-------------------------------------------------------------------------
static inline void givensRotation(const real ff, const real gg, real * restrict cc, real * restrict ss)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

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

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* zero clear subdiagonal entries in R via Givens rotation */
//-------------------------------------------------------------------------
/* Ndim ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* rank ::  input :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* idx  ::  input :: index of the head column which has non-zero subdiagonal element */
/* qq   :: output :: orthogonal matrix Q[Ndim][Ndim] */
/* rr   :: output :: rectangular matrix R[Ndim][rank] */
//-------------------------------------------------------------------------
static inline void clearHessenberg(const int Ndim, const int rank, const int idx, real qq[restrict][NDIM_SEB], real rr[restrict][NDIM_SEB])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int ii = idx; ii < rank; ii++){
    //---------------------------------------------------------------------
    /* compute Givens coefficients */
    real cc, ss;
    givensRotation(rr[ii][ii], rr[ii][ii + 1], &cc, &ss);
    //---------------------------------------------------------------------
    /* rotate R-rows */
    rr[ii][ii] = cc * rr[ii][ii] + ss * rr[ii][ii + 1];/* the other one is an implicit zero */
    for(int jj = ii + 1; jj < rank; jj++){
      //-------------------------------------------------------------------
      const real ff = rr[jj][ii    ];
      const real gg = rr[jj][ii + 1];
      //-------------------------------------------------------------------
      rr[jj][ii    ] = cc * ff + ss * gg;
      rr[jj][ii + 1] = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* for(int jj = ii + 1; jj < rank; jj++){ */
    //---------------------------------------------------------------------
    /* rotate Q-columns */
    for(int jj = 0; jj < Ndim; jj++){
      //-------------------------------------------------------------------
      const real ff = qq[ii    ][jj];
      const real gg = qq[ii + 1][jj];
      //-------------------------------------------------------------------
      qq[ii    ][jj] = cc * ff + ss * gg;
      qq[ii + 1][jj] = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < Ndim; jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = idx; ii < qr->rank; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* update matrices Q and R */
//-------------------------------------------------------------------------
/* Ndim ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* uu   ::  input :: vector u[Ndim] to store the relative position pos[Ndim] - origin[Ndim] */
/* ww   :: output :: vector w[Ndim] to store the matrix vector multiplication Qu */
/* rank ::  input :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* qq   :: output :: orthogonal matrix Q[Ndim][Ndim] */
/* rr   :: output :: rectangular matrix R[Ndim][rank] */
//-------------------------------------------------------------------------
void updateQR(const int Ndim, real uu[restrict], real ww[restrict], const int rank, real qq[restrict][NDIM_SEB], real rr[restrict][NDIM_SEB])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* compute w = Q^{-1} u */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ndim; ii++){
    //---------------------------------------------------------------------
    ww[ii] = ZERO;
    for(int jj = 0; jj < Ndim; jj++)
      ww[ii] += qq[ii][jj] * uu[jj];
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ndim; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* rotate w */
  //-----------------------------------------------------------------------
  for(int jj = Ndim - 1; jj > 0; jj--){
    //---------------------------------------------------------------------
    /* compute Givens coefficients */
    real cc, ss;
    givensRotation(ww[jj - 1], ww[jj], &cc, &ss);
    //---------------------------------------------------------------------
    /* rotate w */
    ww[jj - 1] = cc * ww[jj - 1] + ss * ww[jj];
    //---------------------------------------------------------------------
    /* rotate two R-rows */
    rr[jj - 1][jj    ]  = -ss * rr[jj - 1][jj - 1];
    rr[jj - 1][jj - 1] *=  cc;
    for(int ii = jj; ii < rank; ii++){
      //-------------------------------------------------------------------
      const real ff = rr[ii][jj - 1];
      const real gg = rr[ii][jj    ];
      //-------------------------------------------------------------------
      rr[ii][jj - 1] = cc * ff + ss * gg;
      rr[ii][jj    ] = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* for(int ii = jj; ii < Ndim; ii++){ */
    //---------------------------------------------------------------------
    /* rotate two Q-columns */
    for(int kk = 0; kk < Ndim; kk++){
      //-------------------------------------------------------------------
      const real ff = qq[jj - 1][kk];
      const real gg = qq[jj    ][kk];
      //-------------------------------------------------------------------
      qq[jj - 1][kk] = cc * ff + ss * gg;
      qq[jj    ][kk] = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* for(int kk = 0; kk < Ndim; kk++){ */
    //---------------------------------------------------------------------
  }/* for(int jj = Ndim - 1; jj > 0; jj--){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* update R */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < rank; ii++)
    rr[ii][0] += ww[0];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* clear subdiagonal entries */
  //-----------------------------------------------------------------------
  clearHessenberg(Ndim, rank, 0, qq, rr);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* remove the specified point from support set */
//-------------------------------------------------------------------------
/* Ndim ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* idx  ::  input :: index of the support point to be removed */
/* pos  ::  input :: position of data points p[Ndat * Ndim] */
/* mem  :: output :: indices of data points included in support set */
/* sup  :: output :: true for data points in support set; false otherwise */
/* rank :: output :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* uu   :: output :: vector u[Ndim] to store the relative position pos[Ndim] - origin[Ndim] */
/* ww   :: output :: vector w[Ndim] to store the matrix vector multiplication Qu */
/* qq   :: output :: orthogonal matrix Q[Ndim][Ndim] */
/* rr   :: output :: rectangular matrix R[Ndim][rank] */
//-------------------------------------------------------------------------
static inline void removePoint
(const int Ndim, const int idx, real * restrict pos, bool * restrict sup, int mem[restrict],
 int * restrict rank, real uu[restrict], real ww[restrict], real qq[restrict][NDIM_SEB], real rr[restrict][NDIM_SEB])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* remove the point and update Q and R */
  //-----------------------------------------------------------------------
  sup[mem[idx]] = false;
  //-----------------------------------------------------------------------
  if( idx == *rank ){
    //---------------------------------------------------------------------
    /* remove origin if the remove point is the origin */
    //---------------------------------------------------------------------
    for(int jj = 0; jj < Ndim; jj++)
      uu[jj] = SEB_AFFINE_ORIGIN(Ndim, mem, *rank, jj) - pos[INDEX2D(Ndat, Ndim, mem[(*rank) - 1], jj)];
    //---------------------------------------------------------------------
    (*rank)--;
    //---------------------------------------------------------------------
    updateQR(Ndim, uu, ww, *rank, qq, rr);
    //---------------------------------------------------------------------
  }/* if( idx == *rank ){ */
  //-----------------------------------------------------------------------
  else{
    //---------------------------------------------------------------------
    /* general case: delte column vector from R */
    //---------------------------------------------------------------------
    /* shift righter column of R one line to the left */
    real tmp[NDIM_SEB];
    for(int jj = 0; jj < Ndim; jj++)
      tmp[jj] = rr[idx][jj];
    for(int ii = idx + 1; ii < *rank; ii++){
      //-------------------------------------------------------------------
      for(int jj = 0; jj < Ndim; jj++)
	rr[ii - 1][jj] = rr[ii][jj];
      mem[ii - 1] = mem[ii];
      //-------------------------------------------------------------------
    }/* for(int ii = idx + 1; ii < *rank; ii++){ */
    //---------------------------------------------------------------------
    mem[(*rank) - 1] = mem[*rank];
    (*rank)--;
    for(int jj = 0; jj < Ndim; jj++)
      rr[*rank][jj] = tmp[jj];
    //---------------------------------------------------------------------
    clearHessenberg(Ndim, *rank, idx, qq, rr);
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* remove a data point from support set if it has non-positive affine coefficient */
//-------------------------------------------------------------------------
/* Ndim ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* pos  ::  input :: position of data points p[Ndat * Ndim] */
/* sup  :: output :: true for data points in support set; false otherwise */
/* cen  ::  input :: position of the SEB's center */
/* aff  :: output :: affine coefficents lambda[rank + 1] */
/* rank :: output :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* mem  :: output :: indices of data points included in support set */
/* uu   :: output :: vector u[Ndim] to store the relative position pos[Ndim] - origin[Ndim] */
/* ww   :: output :: vector w[Ndim] to store the matrix vector multiplication Qu */
/* qq   :: output :: orthogonal matrix Q[Ndim][Ndim] */
/* rr   :: output :: rectangular matrix R[Ndim][rank] */
//-------------------------------------------------------------------------
static inline bool dropSupport
(const int Ndim, real * restrict pos, bool * restrict sup, real cen[restrict], real aff[restrict],
 int * restrict rank, int mem[restrict], real uu[restrict], real ww[restrict], real qq[restrict][NDIM_SEB], real rr[restrict][NDIM_SEB])
{
  //-----------------------------------------------------------------------
  /* find affine coefficients of the center */
  //-----------------------------------------------------------------------
  findAffineCoefficients(Ndim, pos, cen, aff, *rank, mem, uu, ww, qq, rr);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* find a non-positive coefficient and drop it */
  //-----------------------------------------------------------------------
  int idx = 0;
  real min = UNITY;
  for(int ii = 0; ii < (*rank) + 1; ii++)
    if( aff[ii] < min ){
      min = aff[ii];
      idx = ii;
    }/* if( aff[ii] < min ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* drop a point with non-positive coefficient */
  //-----------------------------------------------------------------------
  if( min <= ZERO ){
    //---------------------------------------------------------------------
    removePoint(Ndim, idx, pos, sup, mem, rank, uu, ww, qq, rr);
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
/* Ndim      ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* Ndat      ::  input :: Number of data points */
/* pos       ::  input :: position of data points p[Ndat * Ndim] */
/* sup       ::  input :: true for data points in support set; false otherwise */
/* r2        ::  input :: square of the radius of the SEB */
/* r2_to_aff ::  input :: square of the distance from the current center to the circumcenter of aff(T) */
/* cen       ::  input :: position of the SEB's center */
/* cen2aff   ::  input :: vector from the current center to the circumcenter of aff(T) */
/* frac      :: output :: fraction of final radius of the SEB against the initial radius of that */
/* idx       :: output :: index of the stopper */
//-------------------------------------------------------------------------
static inline void findStopper
(const int Ndim, const int Ndat, real * restrict pos, bool * restrict sup,
 const real r2, const real r2_to_aff, real cen[restrict], real cen2aff[restrict],
 real * restrict frac, int * restrict idx)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  *frac = UNITY;/* this is the maximum case; i.e., there is no stopper */
  *idx  = -1;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const real tinyVal = sebEps * r2 * r2_to_aff * RSQRT(sebTiny + r2 * r2_to_aff);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ndat; ii++){
    //---------------------------------------------------------------------
    /* data points in support set are not candidate as a stopper */
    if( !sup[ii] ){
      //-------------------------------------------------------------------
      real cen2rr[NDIM_SEB];
      for(int jj = 0; jj < Ndim; jj++)
	cen2rr[jj] = pos[INDEX2D(Ndat, Ndim, ii, jj)] - cen[jj];
      //-------------------------------------------------------------------
      const real disp = rdot(Ndim, cen2rr, cen2aff);
      const real diff = r2_to_aff - disp;
      /* diff >= 0 is the necessary condition to be a stopper */
      if( diff >= tinyVal ){
	//-----------------------------------------------------------------
	const real bound = HALF * (r2 - rnrm2(Ndim, cen2rr)) * RSQRT(sebTiny + diff * diff);
	if( bound < *frac ){
	  //---------------------------------------------------------------
	  *frac = bound;
	  *idx = ii;
	  //---------------------------------------------------------------
	}/* if( bound < *frac ){ */
	//-----------------------------------------------------------------
      }/* if( diff > tinyVal ){ */
      //-------------------------------------------------------------------
    }/* if( !qr.sup[ii] ){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ndat; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* append a new column to R */
//-------------------------------------------------------------------------
/* Ndim ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* rank ::  input :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* uu   ::  input :: vector u[Ndim] to store the relative position pos[Ndim] - origin[Ndim] */
/* qq   :: output :: orthogonal matrix Q[Ndim][Ndim] */
/* rr   :: output :: rectangular matrix R[Ndim][rank] */
//-------------------------------------------------------------------------
static inline void appendColumn
(const int Ndim, const int rank, real uu[restrict], real qq[restrict][NDIM_SEB], real rr[restrict][NDIM_SEB])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* compute new column appendend to R */
  //-----------------------------------------------------------------------
  for(int jj = 0; jj < Ndim; jj++){
    //---------------------------------------------------------------------
    rr[rank][jj] = ZERO;
    for(int kk = 0; kk < Ndim; kk++)
      rr[rank][jj] += qq[jj][kk] * uu[kk];
    //---------------------------------------------------------------------
  }/* for(int jj = 0; jj < Ndim; jj++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* update QR-decomposition */
  //-----------------------------------------------------------------------
  for(int jj = Ndim - 1; jj > rank; jj--){
    //---------------------------------------------------------------------
    real cc, ss;
    givensRotation(rr[rank][jj - 1], rr[rank][jj], &cc, &ss);
    //---------------------------------------------------------------------
    /* rotate one R-entry (another entry is an implicit zero) */
    rr[rank][jj - 1] = cc * rr[rank][jj - 1] + ss * rr[rank][jj];
    //---------------------------------------------------------------------
    /* rotate two Q-columns */
    for(int kk = 0; kk < Ndim; kk++){
      //-------------------------------------------------------------------
      const real ff = qq[jj - 1][kk];
      const real gg = qq[jj    ][kk];
      //-------------------------------------------------------------------
      qq[jj - 1][kk] = cc * ff + ss * gg;
      qq[jj    ][kk] = cc * gg - ss * ff;
      //-------------------------------------------------------------------
    }/* for(int kk = 0; kk < Ndim; kk++){ */
    //---------------------------------------------------------------------
  }/* for(int jj = Ndim - 1; jj > rank; jj--){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* add the specified stopper to the support set */
//-------------------------------------------------------------------------
/* Ndim ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* pos  ::  input :: position of data points p[Ndat * Ndim] */
/* sup  :: output :: true for data points in support set; false otherwise */
/* idx  ::  input :: index of the specified data point as a support set */
/* rank :: output :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* mem  :: output :: indices of data points included in support set */
/* uu   :: output :: vector u[Ndim] to store the relative position pos[Ndim] - origin[Ndim] */
/* qq   :: output :: orthogonal matrix Q[Ndim][Ndim] */
/* rr   :: output :: rectangular matrix R[Ndim][rank] */
//-------------------------------------------------------------------------
static inline void addSupport
(const int Ndim, real * restrict pos, bool * restrict sup, const int idx,
 int * restrict rank, int mem[restrict], real uu[restrict], real qq[restrict][NDIM_SEB], real rr[restrict][NDIM_SEB])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* compute u := pos[idx] - origin */
  //-----------------------------------------------------------------------
  for(int jj = 0; jj < Ndim; jj++)
    uu[jj] = pos[INDEX2D(Ndat, Ndim, idx, jj)] - SEB_AFFINE_ORIGIN(Ndim, mem, *rank, jj);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* append new column u to R and update QR-decomposition */
  //-----------------------------------------------------------------------
  appendColumn(Ndim, *rank, uu, qq, rr);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* move origin index and insert new index */
  //-----------------------------------------------------------------------
  sup[idx] = true;
  mem[(*rank) + 1] = mem[*rank];
  mem[ *rank     ] = idx;
  //-----------------------------------------------------------------------
  (*rank)++;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* fix floating-point computation origin error on r2 */
//-------------------------------------------------------------------------
/* Ndim :: input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* pos  :: input :: position of data points p[Ndat * Ndim] */
/* cen  :: input :: position of the SEB's center */
/* r2   :: input :: square of the radius of the SEB */
/* rank :: input :: rank of matrix A[rank][rank] = QR, rank = # of support points - 1 */
/* mem  :: input :: indices of data points included in support set */
//-------------------------------------------------------------------------
static inline void finalizeSEBfinder(const int Ndim, real * restrict pos, real cen[restrict], real * restrict r2, const int rank, int mem[restrict])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#if 1
  //-----------------------------------------------------------------------
  /* output information on the support set */
  //-----------------------------------------------------------------------
  fprintf(stdout, "rank = %d for the SEB\n", rank);
  for(int ii = 0; ii < rank + 1; ii++)
    fprintf(stdout, "%d-th support: %d-th particle:\t(%e, %e, %e), %e from center\n", ii, mem[ii],
	    pos[INDEX2D(num, NDIM_SEB, mem[ii], 0)], pos[INDEX2D(num, NDIM_SEB, mem[ii], 1)], pos[INDEX2D(num, NDIM_SEB, mem[ii], 2)],
	    SQRT((pos[INDEX2D(num, NDIM_SEB, mem[ii], 0)] - cen[0]) * (pos[INDEX2D(num, NDIM_SEB, mem[ii], 0)] - cen[0]) +
		 (pos[INDEX2D(num, NDIM_SEB, mem[ii], 1)] - cen[1]) * (pos[INDEX2D(num, NDIM_SEB, mem[ii], 1)] - cen[1]) +
		 (pos[INDEX2D(num, NDIM_SEB, mem[ii], 2)] - cen[2]) * (pos[INDEX2D(num, NDIM_SEB, mem[ii], 2)] - cen[2])));
  //-----------------------------------------------------------------------
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* fix r2 by checking all data points contained in the support set */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < rank + 1; ii++){
    //---------------------------------------------------------------------
    real tmp = ZERO;
    for(int jj = 0; jj < Ndim; jj++){
      //-------------------------------------------------------------------
      const real dr = pos[INDEX2D(Ndat, Ndim, mem[ii], jj)] - cen[jj];
      tmp += dr * dr;
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < Ndim; jj++){ */
    //---------------------------------------------------------------------
    if( tmp > *r2 )
      *r2 = tmp;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < rank + 1; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* find the smallest enclosing ball of a given set of data points */
//-------------------------------------------------------------------------
/* Ndim ::  input :: Number of dimension, must be satisfy ``Ndim <= NDIM_SEB'' */
/* Ndat ::  input :: Number of data points */
/* pos  ::  input :: position of data points p[Ndat * Ndim] */
/* sup  ::        :: true for data points in support set; false otherwise */
/* cen  :: output :: position of the SEB's center */
/* r2   :: output :: square of the radius of the SEB */
//-------------------------------------------------------------------------
void findSEB
(const int Ndim, const int Ndat, real * restrict pos, bool * restrict sup, real cen[restrict], real * restrict r2)
{
  //-----------------------------------------------------------------------
  /* define variables */
  //-----------------------------------------------------------------------
  real qq[NDIM_SEB][NDIM_SEB], rr[NDIM_SEB][NDIM_SEB];/* Q and R */
  real uu[NDIM_SEB], ww[NDIM_SEB];
  real aff[NDIM_SEB + 1];/* affine coefficents, maximum # of support points is d + 1 */
  real cen2aff[NDIM_SEB];/* vector from the current center to the circumcenter of aff(T) */
  int mem[NDIM_SEB + 1];/* indices of points in support set */
  //-----------------------------------------------------------------------
  int rank;
  real r2_to_aff;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  initSEB(Ndim, Ndat, pos, cen, r2, sup, qq, rr, uu, ww, mem, &rank);
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
#if 1
    fprintf(stdout, "rank = %d:\tsupport = %d", rank, mem[0]);
    for(int ii = 1; ii < rank + 1; ii++)
      fprintf(stdout, ", %d", mem[ii]);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    //---------------------------------------------------------------------
    /* compute a walking direction */
    //---------------------------------------------------------------------
    while( (r2_to_aff = getShortestVector(Ndim, pos, cen, cen2aff, mem, rank, qq)) <= sebEps2 * (*r2) )
      if( !dropSupport(Ndim, pos, sup, cen, aff, &rank, mem, uu, ww, qq, rr) ){
	/* if dropSupport returns false, it means the center lies in the convex hull; i.e., SEB is already found */
	finalizeSEBfinder(Ndim, pos, cen, r2, rank, mem);
	return;
      }/* if( !dropSupport(Ndim, pos, sup, cen, aff, &rank, mem, uu, ww, qq, rr) ){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* find the stopper in data points */
    //---------------------------------------------------------------------
    real stopFrac;
    int  stopIdx;
    findStopper(Ndim, Ndat, pos, sup, *r2, r2_to_aff, cen, cen2aff, &stopFrac, &stopIdx);
    //---------------------------------------------------------------------
    /* if a stopping point exists */
    if( (stopIdx >= 0) && ((1 + rank) <= Ndim) ){
      //-------------------------------------------------------------------
      /* move the center */
      for(int ii = 0; ii < Ndim; ii++)
	cen[ii] += stopFrac * cen2aff[ii];
      //-------------------------------------------------------------------
      /* update the radius */
      *r2 = ZERO;
      for(int ii = 0; ii < Ndim; ii++){
	const real dr = pos[INDEX2D(Ndat, Ndim, mem[rank], ii)] - cen[ii];
	*r2 += dr * dr;
      }/* for(int ii = 0; ii < Ndim; ii++){ */
      //-------------------------------------------------------------------
      /* add the stopper to the support set */
      addSupport(Ndim, pos, sup, stopIdx, &rank, mem, uu, qq, rr);
      //-------------------------------------------------------------------
    }/* if( (stopIdx >= 0) && ((1 + rank) <= Ndim) ){ */
    //---------------------------------------------------------------------
    /* if there is no stopping point */
    else{
      //-------------------------------------------------------------------
      /* move the center */
      for(int ii = 0; ii < Ndim; ii++)
	cen[ii] += cen2aff[ii];
      //-------------------------------------------------------------------
      /* update the radius */
      *r2 = ZERO;
      for(int ii = 0; ii < Ndim; ii++){
	const real dr = pos[INDEX2D(Ndat, Ndim, mem[rank], ii)] - cen[ii];
	*r2 += dr * dr;
      }/* for(int ii = 0; ii < Ndim; ii++){ */
      //-------------------------------------------------------------------
      if( !dropSupport(Ndim, pos, sup, cen, aff, &rank, mem, uu, ww, qq, rr) ){
	finalizeSEBfinder(Ndim, pos, cen, r2, rank, mem);
	return;
      }/* if( !dropSupport(Ndim, pos, sup, cen, aff, &rank, mem, uu, ww, qq, rr) ){ */
      //-------------------------------------------------------------------
    }/* else{ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //-----------------------------------------------------------------------
  /* initialization */
  //-----------------------------------------------------------------------
  if( argc < 2 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 2);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -num=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 2 ){ */
  //-----------------------------------------------------------------------
  int num;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "num", &num));
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* GSL initialization */
  //-----------------------------------------------------------------------
  /* initialize random number provided by GSL */
  const gsl_rng_type *RandType;
  gsl_rng_env_setup();
  RandType = gsl_rng_mt19937;
  GSLRand  = gsl_rng_alloc(RandType);
  gsl_rng_set(GSLRand, 5489);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation */
  //-----------------------------------------------------------------------
  /* data points */
  real *pos;  pos = (real *)malloc(sizeof(real) * (size_t)num * (size_t)NDIM_SEB);
  bool *sup;  sup = (bool *)malloc(sizeof(bool) * (size_t)num);
  if( pos == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pos\n");  }
  if( sup == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sup\n");  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* initialization of data points */
  /* generate randomly distributed points in [-1/2, 1/2] */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++)
    for(int jj = 0; jj < NDIM_SEB; jj++)
      pos[INDEX2D(num, NDIM_SEB, ii, jj)] = HALF * RANDVAL;
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < num; ii++){
    fprintf(stdout, "%d", ii);
    for(int jj = 0; jj < NDIM_SEB; jj++)
      fprintf(stdout, "\t%e", pos[INDEX2D(num, NDIM_SEB, ii, jj)]);
    fprintf(stdout, "\n");
  }
#endif
  //-----------------------------------------------------------------------
#if 1
  for(int jj = 0; jj < NDIM_SEB; jj++)    pos[INDEX2D(num, NDIM_SEB, 0, jj)] =  HALF;
  for(int jj = 0; jj < NDIM_SEB; jj++)    pos[INDEX2D(num, NDIM_SEB, 1, jj)] = -HALF;
#endif
#if 1
  for(int ii = 2; ii < num >> 1; ii++)
    for(int jj = 0; jj < NDIM_SEB; jj++)
      pos[INDEX2D(num, NDIM_SEB, ii, jj)] = pos[INDEX2D(num, NDIM_SEB, 0, jj)];
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  real cen[NDIM_SEB];
  real r2;
  //-----------------------------------------------------------------------
  findSEB(NDIM_SEB, num, pos, sup, cen, &r2);
  //-----------------------------------------------------------------------
  fprintf(stdout, "rad = %e, center = (%e, %e, %e)\n", SQRT(r2), cen[0], cen[1], cen[2]);
  //-----------------------------------------------------------------------
  /* confirmation */
  real max = ZERO;
  int maxIdx = 0;
  for(int ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    real r2ful = ZERO;
    for(int jj = 0; jj < NDIM_SEB; jj++){
      const real dr = pos[INDEX2D(num, NDIM_SEB, ii, jj)] - cen[jj];
      r2ful += dr * dr;
    }/* for(int jj = 0; jj < NDIM_SEB; jj++){ */
    if( r2ful > max ){
      max = r2ful;
      maxIdx = ii;
    }/* if( r2ful > max ){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < num; ii++){ */
  fprintf(stdout, "maximum rad = %e with %d-th particle\n", SQRT(max), maxIdx);
  //-----------------------------------------------------------------------



  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(sup);
  free(pos);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
