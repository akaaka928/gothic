/*************************************************************************\
 *                                                                       *
                  last updated on 2016/04/25(Mon) 14:53:00
 *                                                                       *
 *    Making Initial Condition Code of N-body Simulation                 *
 *       King sphere                                                     *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
/* #define TEST_CUBIC_SPLINE */
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <constants.h>
//-------------------------------------------------------------------------
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "magi.h"
#include "king.h"
#ifdef  TEST_CUBIC_SPLINE
#include "spline.h"
#include "blas.h"
#endif//TEST_CUBIC_SPLINE
//-------------------------------------------------------------------------
extern const real newton;
//-------------------------------------------------------------------------
static const double inv3 = 1.0 / 3.0;
static const double inv6 = 1.0 / 6.0;
//-------------------------------------------------------------------------
static const double convergence = 1.0e-4;
/* static const double convergence = 1.0e-5; */
static const double extreme     = 0.25;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define NADD_KING (16384) */
/* #define NADD_KING (131072) */
#define NADD_KING (1048576)
/* #define NADD_KING (4 * 1048576) */
//-------------------------------------------------------------------------
static inline void allocateArray(const int num, double **rad, double **rho, double **psi, double **dr1, double **dr2)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  *rad = (double *)malloc(sizeof(double) * num);  if( *rad == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rad");  }
  *rho = (double *)malloc(sizeof(double) * num);  if( *rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rho");  }
  *psi = (double *)malloc(sizeof(double) * num);  if( *psi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate psi");  }
  *dr1 = (double *)malloc(sizeof(double) * num);  if( *dr1 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dr1");  }
  *dr2 = (double *)malloc(sizeof(double) * num);  if( *dr2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dr2");  }
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void enlargeArray(const int num, double **rad, double **rho, double **psi, double **dr1, double **dr2)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  *rad = realloc(*rad, sizeof(double) * num);  if( *rad == NULL ){    __KILL__(stderr, "ERROR: failure to rescale rad");  }
  *rho = realloc(*rho, sizeof(double) * num);  if( *rho == NULL ){    __KILL__(stderr, "ERROR: failure to rescale rho");  }
  *psi = realloc(*psi, sizeof(double) * num);  if( *psi == NULL ){    __KILL__(stderr, "ERROR: failure to rescale psi");  }
  *dr1 = realloc(*dr1, sizeof(double) * num);  if( *dr1 == NULL ){    __KILL__(stderr, "ERROR: failure to rescale dr1");  }
  *dr2 = realloc(*dr2, sizeof(double) * num);  if( *dr2 == NULL ){    __KILL__(stderr, "ERROR: failure to rescale dr2");  }
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void  releaseArray(               double  *rad, double  *rho, double  *psi, double  *dr1, double  *dr2)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  free(rad);
  free(rho);
  free(psi);
  free(dr1);
  free(dr2);
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline double kingErrFunc(const double psi, const double sqrt_psi)
{
  //-----------------------------------------------------------------------
  return (exp(psi) * erf(sqrt_psi));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double kingFunc0(const double sqrt_psi, const double kingExpErrFunc)
{
  //-----------------------------------------------------------------------
  return (kingExpErrFunc - M_2_SQRTPI * sqrt_psi);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double getKingDensity(const double psi, const double sqrt_psi, const double kingExpErrFunc, const double rho1)
{
  //-----------------------------------------------------------------------
  return (rho1 * (kingExpErrFunc - M_2_SQRTPI * sqrt_psi * (1.0 + 2.0 * psi * inv3)));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline double calcKingDensity(double psi, double rho1)
{
  //-----------------------------------------------------------------------
  const double sqrt_psi = sqrt(psi);
  return (rho1 * (kingErrFunc(psi, sqrt_psi) - M_2_SQRTPI * sqrt_psi * (1.0 + 2.0 * psi * inv3)));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
#define kingFunc1(rad, uu, yy) _kingFunc1(rad, yy)
static inline double _kingFunc1(const double rad, const double yy)
{
  //-----------------------------------------------------------------------
  return ((rad > DBL_EPSILON) ? (yy / (rad * rad)) : (0.0));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#define kingFunc2(rad, psi, dWdr, rho0inv, rho1) _kingFunc2(rad, psi, rho0inv, rho1)
static inline double _kingFunc2(const double rad, const double psi, const double rho0inv, const double rho1)
{
  //-----------------------------------------------------------------------
  return  (-9.0 * rad * rad * rho0inv * calcKingDensity(psi, rho1));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void rungeKutta4thForKing(const double rr, double *uu, double *yy, const double hh, const double rho0inv, const double rho1)
{
  //-----------------------------------------------------------------------
  const double k1u = kingFunc1(rr,            *uu,                  *yy                                );
  const double k1y = kingFunc2(rr,            *uu,                  *yy,                  rho0inv, rho1);
  const double k2u = kingFunc1(rr + 0.5 * hh, *uu + 0.5 * hh * k1u, *yy + 0.5 * hh * k1y               );
  const double k2y = kingFunc2(rr + 0.5 * hh, *uu + 0.5 * hh * k1u, *yy + 0.5 * hh * k1y, rho0inv, rho1);
  const double k3u = kingFunc1(rr + 0.5 * hh, *uu + 0.5 * hh * k2u, *yy + 0.5 * hh * k2y               );
  const double k3y = kingFunc2(rr + 0.5 * hh, *uu + 0.5 * hh * k2u, *yy + 0.5 * hh * k2y, rho0inv, rho1);
  const double k4u = kingFunc1(rr +       hh, *uu +       hh * k3u, *yy +       hh * k3y               );
  const double k4y = kingFunc2(rr +       hh, *uu +       hh * k3u, *yy +       hh * k3y, rho0inv, rho1);
  //-----------------------------------------------------------------------
  *uu += (hh * inv6) * (k1u + 2.0 * k2u + 2.0 * k3u + k4u);
  *yy += (hh * inv6) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
/* W is the dimensionless King parameter Psi / sigma^2 */
static inline void solvePoissonEqOfKingDF(const double W0, double **Rad, double **Psi, double **Rho, double **Dr1, double **Dr2, int *num, int *rem)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  double *rad, *psi, *rho, *dr1, *dr2;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set initial condition */
  //-----------------------------------------------------------------------
  double hh = MINRAD * 0.25;
  /* double hh = 1.0 / 1024.0; */
  double hnew = hh;
  //-----------------------------------------------------------------------
  /* assumption is rcal := r / r0, r0 = 1.0, sigma = 1.0 */
  //-----------------------------------------------------------------------
  /* see Eq (4.106) in Galactic Dynamics, second edtion, page 305 */
  /* rho0 = (real)2.25 * sigma * sigma / (pi * newton * r0 * r0); */
  const double rho0 = 2.25 * M_1_PI;
  const double rho0inv = 1.0 / rho0;
  //-----------------------------------------------------------------------
  /* uu stores the dimensionless King parameter W, yy stores r^2 * dW/dr */
  double uu = W0;
  /* boundary condition is dW/dr(r = 0) = 0: assume core profile at the center */
  double yy = 0.0;
  //-----------------------------------------------------------------------
  int ii = 0;
  rad = *Rad;  rad[ii] = 0.0;
  psi = *Psi;  psi[ii] = uu;
  rho = *Rho;  rho[ii] = calcKingDensity(psi[ii], 1.0);
  dr1 = *Dr1;  dr1[ii] = 0.0;
  dr2 = *Dr2;  dr2[ii] = 0.0;
  //-----------------------------------------------------------------------
  const double rho1 = rho0 / rho[ii];
  rho[ii] = calcKingDensity(psi[ii], rho1);
  //-----------------------------------------------------------------------
  double rold = rho[ii];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* solve Poisson Equation using 4th-order Runge-Kutta method */
  //-----------------------------------------------------------------------
#ifdef  KING_PROGRESS_REPORT_ON
  fprintf(stdout, "# Poisson solver for King model start\n");
  fflush(stdout);
#endif//KING_PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    hh = hnew;
    const double uold = uu;
    const double yold = yy;
    const double uoldinv = 1.0 / uold;
    const double yoldinv = 1.0 / (1.0e-30 + fabs(yold));
    const double roldinv = 1.0 / rold;
    //---------------------------------------------------------------------
    double dens, sqrtW, kingExpErrFunc;
    //---------------------------------------------------------------------
    while( true ){
      //-------------------------------------------------------------------
      uu = uold;
      yy = yold;
      rungeKutta4thForKing(rad[ii], &uu, &yy, hh, rho0inv, rho1);
      //-------------------------------------------------------------------
      sqrtW = sqrt(uu);
      kingExpErrFunc = kingErrFunc(uu, sqrtW);
      dens = getKingDensity(uu, sqrtW, kingExpErrFunc, rho1);
      //-------------------------------------------------------------------
      /* convergence tests */
#if 1
      const double udiff = fabs(uu - uold) * uoldinv;
      const double ydiff = fabs(yy - yold) * yoldinv;
      double  diff = (udiff > ydiff) ? udiff : ydiff;
      /* const double rdiff = (convergence * dens > DBL_EPSILON) ? (fabs(dens - rold) * roldinv) : (0.0); */
      /* const double rdiff = (dens > DBL_EPSILON * rho1) ? (fabs(dens - rold) * roldinv) : (0.0); */
      const double rdiff = (convergence * dens > DBL_EPSILON * rho[0]) ? (fabs(dens - rold) * roldinv) : (0.0);
      if( rdiff > diff )
	diff = rdiff;
#else
      const double diff = fabs(uu - uold) * uoldinv;
#endif
      if( diff < convergence ){
	if( diff < extreme * convergence )	  hnew = hh * 2.0;
	break;
      }
      else{
	hh *= 0.5;
	hnew = hh;
      }
      //-------------------------------------------------------------------
    }
    //---------------------------------------------------------------------
    if( *rem == 0 ){
      enlargeArray((*num) + NADD_KING, Rad, Rho, Psi, Dr1, Dr2);
      rad = *Rad;
      psi = *Psi;
      rho = *Rho;
      dr1 = *Dr1;
      dr2 = *Dr2;
      *rem += NADD_KING;
    }
    ii++;    *num += 1;    *rem -= 1;
    rad[ii] = rad[ii - 1] + hh;
    //---------------------------------------------------------------------
    /* if( (uu > DBL_EPSILON) && (dens > rho[0] * 1.0e-12) ){ */
    if( (uu > DBL_EPSILON) && (dens > DBL_EPSILON) ){
      //-------------------------------------------------------------------
      psi[ii] = uu;
      rho[ii] = dens;
      rold = dens;
      //-------------------------------------------------------------------
      const double tmp = kingFunc0(sqrtW, kingExpErrFunc);
      const double rinv = 1.0 / rad[ii];
      const double dW_dr = yy * rinv * rinv;
      dr1[ii] = rho1 *  tmp * dW_dr;
#if 1
      const double d2W_dr2 = -9.0 * dens * rho0inv - 2.0 * dW_dr * rinv;
      dr2[ii] = rho1 * (tmp * d2W_dr2 + kingExpErrFunc * dW_dr * dW_dr);
#else
      dr2[ii] = rho1 * ((kingExpErrFunc * dW_dr - 2.0 * tmp * rinv) * dW_dr - 9.0 * dens * rho0inv * tmp);
#endif
    }
    //---------------------------------------------------------------------
    else{
      psi[ii] = 0.0;
      rho[ii] = 0.0;
      dr1[ii] = 0.0;
      dr2[ii] = 0.0;
      break;
    }
    //---------------------------------------------------------------------
#ifdef  KING_PROGRESS_REPORT_ON
    if( (ii % KING_PROGRESS_REPORT_ON) == 0 ){
      fprintf(stdout, "# Poisson solver for King model: %d steps finished: rad = %e, rho = %e, W = %e\n", ii, rad[ii], rho[ii], psi[ii]);
      fflush(stdout);
    }
#endif//KING_PROGRESS_REPORT_ON
    //---------------------------------------------------------------------
#if 0
    if( (ii & 127) == 0 )
      fprintf(stderr, "%e\t%e\t%e\n", rad[ii], rho[ii], psi[ii]);
#endif
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#ifdef  KING_PROGRESS_REPORT_ON
  fprintf(stdout, "# Poisson solver for King model finish: %d elements used: rt / r0 = %e, c = %e\n#\n#\n", *num, rad[*num - 1], log10(rad[*num - 1]));
  fflush(stdout);
#endif//KING_PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------
#ifdef  KING_PROGRESS_REPORT_ON
  fprintf(stdout, "# King model: rho1 is %le w/o scaling\n", rho1);
  fflush(stdout);
#endif//KING_PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#if 0
  for(int jj = 1; jj < *num - 1; jj += 100)
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", rad[jj], dr1[jj], (rho[jj + 1] - rho[jj - 1]) / (rad[jj + 1] - rad[jj - 1]), dr2[jj], (dr1[jj + 1] - dr1[jj - 1]) / (rad[jj + 1] - rad[jj - 1]));
  exit(0);
#endif
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static void rescaleKingSphere(const double Mtot, const double r0, double *rt, const int num, double *rad, double *rho, double *dr1, double *dr2, double *enc)
{
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 1; ii < num - 1; ii += 100)
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", rad[ii], dr1[ii], (rho[ii + 1] - rho[ii - 1]) / (rad[ii + 1] - rad[ii - 1]), dr2[ii], (dr1[ii + 1] - dr1[ii - 1]) / (rad[ii + 1] - rad[ii - 1]));
  exit(0);
#endif
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  double MencCalc = 2.0 * M_PI;/* 4pi / 2 */
  enc[0] = MencCalc * rad[1] * rad[1] * rad[1] * (rho[0] + rho[1]) * inv3;
  for(int ii = 1; ii < num - 1; ii++)
    enc[ii] = enc[ii - 1] + MencCalc * rad[ii] * rad[ii + 1] * (rad[ii + 1] - rad[ii]) * (rho[ii] + rho[ii + 1]);
    /* enc[ii] = enc[ii - 1] + MencCalc * rad[ii + 1] * rad[ii + 1] * (rad[ii + 1] - rad[ii]) * (rho[ii] + rho[ii + 1]) + enc[ii - 1]; */
  MencCalc = enc[num - 1] = enc[num - 2];
  double r0Calc = 1.0;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* evaluate scaling factors */
  //-----------------------------------------------------------------------
  const double lengthUnit = r0   /   r0Calc;
  const double   massUnit = Mtot / MencCalc;
  //-----------------------------------------------------------------------
  const double    rhoUnit = massUnit / (lengthUnit * lengthUnit * lengthUnit);
  //-----------------------------------------------------------------------
  const double  drho_drUnit  = rhoUnit / lengthUnit;
  const double d2rho_dr2Unit = rhoUnit / (lengthUnit * lengthUnit);
  //-----------------------------------------------------------------------
  *rt = rad[num - 1] * lengthUnit;
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < num; ii++){
    rad[ii] *=    lengthUnit;
    rho[ii] *=       rhoUnit;
    enc[ii] *=      massUnit;
    dr1[ii] *=   drho_drUnit;
    dr2[ii] *= d2rho_dr2Unit;
  }/* for(int ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 1; ii < num - 1; ii += 100)
    fprintf(stderr, "%e\t%e\t%e\t%e\n", rad[ii], rho[ii], dr1[ii], dr2[ii]);
  exit(0);
#endif
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 1; ii < num - 1; ii += 100)
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", rad[ii], dr1[ii], (rho[ii + 1] - rho[ii - 1]) / (rad[ii + 1] - rad[ii - 1]), dr2[ii], (dr1[ii + 1] - dr1[ii - 1]) / (rad[ii + 1] - rad[ii - 1]));
  exit(0);
#endif
  //-----------------------------------------------------------------------
#ifdef  KING_PROGRESS_REPORT_ON
  fprintf(stdout, "# King model:   lengthUnit is %le\n", lengthUnit);
  fprintf(stdout, "# King model:     massUnit is %le\n",   massUnit);
  fprintf(stdout, "# King model:      rhoUnit is %le\n",    rhoUnit);
  fprintf(stdout, "# King model: velocityUnit is %le\n", lengthUnit * sqrt((double)newton * rhoUnit));
  fprintf(stdout, "#\n#\n");
  fflush(stdout);
#endif//KING_PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline int findIdx(const double rad, profile *prf)
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = 3 + NRADBIN;
  //-----------------------------------------------------------------------
  if( rad < prf[ll].rad + DBL_EPSILON ){    return (ll    );  }
  if( rad > prf[rr].rad - DBL_EPSILON ){    return (rr - 1);  }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)ll + (uint)rr) >> 1;
    //---------------------------------------------------------------------
    if( (prf[cc].rad - rad) * (prf[ll].rad - rad) <= 0.0 )      rr = (int)cc;
    else                                                        ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr )
      return (ll);
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline int bisection(const double val, const int num, double tab[], double *ratio)
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = num - 1;
  //-----------------------------------------------------------------------
  /* prohibit extraporation */
  if( val < tab[ll] + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( val > tab[rr] - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)(ll + rr)) >> 1;
    //---------------------------------------------------------------------
    if( (tab[cc] - val) * (tab[ll] - val) <= 0.0)      rr = (int)cc;
    else                                               ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr ){
      *ratio = (val - tab[ll]) / (tab[1 + ll] - tab[ll]);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void getDensityProfile(const int num, double *rad, double *rho, double *dr1, double *dr2, profile *prf)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate spherical averaged profile of density, enclosed mass and potential (internal part) */
  //-----------------------------------------------------------------------
  const int head = findIdx(rad[      0], prf);
  const int tail = findIdx(rad[num - 1], prf);
  //-----------------------------------------------------------------------
#if 0
  /* this part would have some bugs about least squared method */
  /* 2nd derivative is given by the least squared method */
  double SS, Sx, Sy, Sxx, Sxy;
  SS = Sx = Sy = Sxx = Sxy = 0.0;
  for(int ii = 0; ii < 8; ii++){
    //---------------------------------------------------------------------
    const double logx = log10(rho[ii]);
    const double logy = log10(dr2[ii]);
    SS  += 1.0;
    Sx  += logx;
    Sxx += logx * logx;
    Sy  +=        logy;
    Sxy += logx * logy;
    //---------------------------------------------------------------------
  }
  const double pp   = (SS * Sxy - Sx * Sy) / (SS * Sxx - Sx * Sx);
  const double bb   = pow(10.0, (Sy - pp * Sx) / SS);
  const double pinv = 1.0 / pp;
  const double binv = 1.0 / bb;
  if( (pp != -1.0) && (pp != -2.0) ){
    //---------------------------------------------------------------------
    const double p1inv = 1.0 / (1.0 + pp);
    const double p2inv = 1.0 / (2.0 + pp);
    //---------------------------------------------------------------------
    const double cc =               dr1[0] - bb      * pp      * p1inv         * pow(rad[0] * binv, 1.0 + pinv);
    const double dd = rho[0] - cc * rad[0] - bb * bb * pp * pp * p1inv * p2inv * pow(rad[0] * binv, 2.0 + pinv);
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < head; ii++){
      //-------------------------------------------------------------------
      const double xx = prf[ii].rad * binv;
      const double yy = pow(xx, pinv);
      //-------------------------------------------------------------------
      prf[ii].d2rho_dr2 =                                                                       yy;
      prf[ii]. drho_dr  =      cc               + bb      * pp      * p1inv         * xx      * yy;
      prf[ii].  rho     = dd + cc * prf[ii].rad + bb * bb * pp * pp * p1inv * p2inv * xx * xx * yy;
    }/* for(int ii = 0; ii < head; ii++){ */
    //---------------------------------------------------------------------
  }/* if( (pp != -1.0) && (pp != -2.0) ){ */
  else{
    if( pp == -2.0 ){
      //-------------------------------------------------------------------
      const double p1inv = -1.0;
      const double cc = dr1[0] - bb * pp * p1inv * pow(rad[0] * binv, 1.0 + pinv);
      const double dd = rho[0] - cc * rad[0] - bb * rad[0] * (log(rad[0]) - 1.0);
      //-------------------------------------------------------------------
#pragma omp parallel for
      for(int ii = 0; ii < head; ii++){
	//-----------------------------------------------------------------
	const double xx = prf[ii].rad * binv;
	const double yy = pow(xx, pinv);
	//-----------------------------------------------------------------
	prf[ii].d2rho_dr2 =                                                     yy;
	prf[ii]. drho_dr  =      cc               + bb      * pp * p1inv * xx * yy;
	prf[ii].  rho     = dd + cc * prf[ii].rad + bb * bb * pp * p1inv * log(prf[ii].rad);
      }/* for(int ii = 0; ii < head; ii++){ */
      //-----------------------------------------------------------------
    }/* if( pp == -2.0 ){ */
    else{
      //-------------------------------------------------------------------
      const double cc = dr1[0] - bb * log(rad[0]);
      const double dd = rho[0] - cc * rad[0] - bb * rad[0] * (log(rad[0]) - 1.0);
      //-------------------------------------------------------------------
#pragma omp parallel for
      for(int ii = 0; ii < head; ii++){
	//-----------------------------------------------------------------
	const double xx = prf[ii].rad * binv;
	const double yy = pow(xx, pinv);
	const double zz = log(prf[ii].rad);
	//-----------------------------------------------------------------
	prf[ii].d2rho_dr2 = yy;
	prf[ii]. drho_dr  = cc + bb * zz;
	prf[ii].  rho     = dd + prf[ii].rad * (cc + bb * (zz - 1.0));
      }/* for(int ii = 0; ii < head; ii++){ */
      //-------------------------------------------------------------------
    }/* else{ */
  }/* else{ */
#else
#pragma omp parallel for
  for(int ii = 0; ii < head; ii++){
    prf[ii].rho       = rho[0];
    prf[ii].drho_dr   = dr1[0];
    prf[ii].d2rho_dr2 = dr2[0];
  }
#endif
  //-----------------------------------------------------------------------
  /* interpolate: enc --> prf */
#pragma omp parallel for
  for(int ii = head; ii < tail; ii++){
    double alpha;
    const int idx = bisection(prf[ii].rad, num, rad, &alpha);
    prf[ii].rho       = (1.0 - alpha) * rho[idx] + alpha * rho[1 + idx];
#if 1
    prf[ii].drho_dr   = (1.0 - alpha) * dr1[idx] + alpha * dr1[1 + idx];
    prf[ii].d2rho_dr2 = (1.0 - alpha) * dr2[idx] + alpha * dr2[1 + idx];
#else
    prf[ii].drho_dr   = (rho[idx + 1] - rho[idx - 1]) / (rad[idx + 1] - rad[idx - 1]);
    prf[ii].d2rho_dr2 = (1.0 - alpha) * dr2[idx] + alpha * dr2[1 + idx];
#endif
  }/* for(int ii = head; ii < tail; ii++){ */
  //-----------------------------------------------------------------------
  /* fill the total enclosed mass in the outer region */
#pragma omp parallel for
  for(int ii = tail; ii < 4 + NRADBIN; ii++){
    prf[ii].rho       = 0.0;
    prf[ii].drho_dr   = 0.0;
    prf[ii].d2rho_dr2 = 0.0;
  }
  //-----------------------------------------------------------------------
#if 0
  const int ih = findIdx(1.0e-2, prf);
#pragma omp parallel for
  for(int ii = ih; ii < 3 + NRADBIN; ii++)
    prf[ii].drho_dr = (prf[ii + 1].rho - prf[ii - 1].rho) / (prf[ii + 1].rad - prf[ii - 1].rad);
#pragma omp parallel for
  for(int ii = 1 + ih; ii < 2 + NRADBIN; ii++)
    prf[ii].d2rho_dr2 = (prf[ii + 1].drho_dr - prf[ii - 1].drho_dr) / (prf[ii + 1].rad - prf[ii - 1].rad);
#endif
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < 4 + NRADBIN; ii++)
    fprintf(stderr, "%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].drho_dr, prf[ii].d2rho_dr2);
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileKing(profile *prf, profile_cfg *cfg)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const double W0   = cfg->king_W0;
  const double Mtot = cfg->Mtot;
  const double r0   = cfg->rs;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* memory allocation */
  //-----------------------------------------------------------------------
  int num = 0;
  int rem = NADD_KING;
  //-----------------------------------------------------------------------
  double *rad, *rho, *psi, *dr1, *dr2;
  allocateArray(rem, &rad, &rho, &psi, &dr1, &dr2);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* derive density profile of the King sphere */
  //-----------------------------------------------------------------------
  solvePoissonEqOfKingDF(W0, &rad, &psi, &rho, &dr1, &dr2, &num, &rem);
#if 0
  for(int ii = 1; ii < num - 1; ii += 100)
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", rad[ii], dr1[ii], (rho[ii + 1] - rho[ii - 1]) / (rad[ii + 1] - rad[ii - 1]), dr2[ii], (dr1[ii + 1] - dr1[ii - 1]) / (rad[ii + 1] - rad[ii - 1]));
  exit(0);
#endif
  double *enc;  enc = (double *)malloc(sizeof(double) * num);  if( enc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate enc");  }
  rescaleKingSphere(Mtot, r0, &cfg->king_rt, num, rad, rho, dr1, dr2, enc);
  //-----------------------------------------------------------------------
  cfg->king_c = log10(cfg->king_rt / cfg->rs);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* test cubic spline */
  //-----------------------------------------------------------------------
#ifdef  TEST_CUBIC_SPLINE
  const int ncs = num / 1024;
  double *xx;  xx = (double *)malloc(sizeof(double) * ncs);  if( xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xx");  }
  double *yy;  yy = (double *)malloc(sizeof(double) * ncs);  if( yy == NULL ){    __KILL__(stderr, "ERROR: failure to allocate yy");  }
  double *bp;  bp = (double *)malloc(sizeof(double) * ncs);  if( bp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate bp");  }
  double *y2;  y2 = (double *)malloc(sizeof(double) * ncs);  if( y2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate y2");  }
  for(int ii = 0; ii < ncs; ii++){
    xx[ii] = rad[ii * 1024];
    yy[ii] = rho[ii * 1024];
  }
  genCubicSpline1D(ncs, xx, yy, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, y2);
#if 1
  //-----------------------------------------------------------------------
  crs mat, ilu;
  double *mat_val, *ilu_val;
  int    *mat_col, *ilu_col, *mat_row, *ilu_row;
  mat_val = (double *)malloc(sizeof(double) * ncs * 3);  ilu_val = (double *)malloc(sizeof(double) * ncs * 3);  mat_row = (int *)malloc(sizeof(int) * (ncs + 1));
  mat_col = (   int *)malloc(sizeof(   int) * ncs * 3);  ilu_col = (   int *)malloc(sizeof(   int) * ncs * 3);  ilu_row = (int *)malloc(sizeof(int) * (ncs + 1));
  mat.val = mat_val;  mat.col = mat_col;  mat.row = mat_row;
  ilu.val = ilu_val;  ilu.col = ilu_col;  ilu.row = ilu_row;
  //-----------------------------------------------------------------------
  /* vector used in BiCGSTAB method */
  double *vec, *res, *sdw, *mid, *tmp, *sol;
  double *Api, *Ati, *Kri, *Kpi, *Kti;
  vec = (double *)malloc(sizeof(double) * ncs);  res = (double *)malloc(sizeof(double) * ncs);  sol = (double *)malloc(sizeof(double) * ncs);
  sdw = (double *)malloc(sizeof(double) * ncs);  mid = (double *)malloc(sizeof(double) * ncs);
  tmp = (double *)malloc(sizeof(double) * ncs);  Api = (double *)malloc(sizeof(double) * ncs);
  Ati = (double *)malloc(sizeof(double) * ncs);  Kri = (double *)malloc(sizeof(double) * ncs);
  Kpi = (double *)malloc(sizeof(double) * ncs);  Kti = (double *)malloc(sizeof(double) * ncs);
  //-----------------------------------------------------------------------
  /* assume natural cubic spline */
#pragma omp parallel for
  for(int ii = 0; ii < ncs; ii++)
    sol[ii] = 0.0;
  vec[0] = 0.0;
#pragma omp parallel for
  for(int ii = 1; ii < ncs - 1; ii++)
    vec[ii] = (yy[ii + 1] - yy[ii]) / (xx[ii + 1] - xx[ii]) - (yy[ii] - yy[ii - 1]) / (xx[ii] - xx[ii - 1]);
  vec[ncs - 1] = 0.0;
  int valIdx = 0;
  int rowIdx = 0;
  mat.row[rowIdx] = valIdx;
  mat.col[valIdx] = 0;  mat.val[valIdx] = 1.0;  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
  for(int ii = 1; ii < ncs - 1; ii++){
    mat.col[valIdx] = ii - 1;    mat.val[valIdx] = (xx[ii    ] - xx[ii - 1]) / 6.0;    valIdx++;
    mat.col[valIdx] = ii    ;    mat.val[valIdx] = (xx[ii + 1] - xx[ii - 1]) / 3.0;    valIdx++;
    mat.col[valIdx] = ii + 1;    mat.val[valIdx] = (xx[ii + 1] - xx[ii    ]) / 6.0;    valIdx++;
    rowIdx++;    mat.row[rowIdx] = valIdx;
  }
  mat.col[valIdx] = ncs - 1;  mat.val[valIdx] = 1.0;  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
  //-----------------------------------------------------------------------
  /* outer boundary condition of the potential field */
  getILU0(ncs, mat, ilu);
  pbicgstab(mat, ncs, vec, sol, res, sdw, mid, tmp, Api, Ati, ilu, Kri, Kpi, Kti, 1.0e-10);
  for(int ii = 0; ii < ncs; ii++)
    fprintf(stderr, "%e\t%e\t%e\t%e\n", xx[ii], yy[ii], y2[ii], sol[ii]);
  free(mat_val);  free(mat_col);  free(mat_row);
  free(ilu_val);  free(ilu_col);  free(ilu_row);
  free(vec);  free(res);  free(sdw);  free(mid);  free(tmp);  free(sol);
  free(Api);  free(Ati);  free(Kri);  free(Kpi);  free(Kti);
  exit(0);
#endif
  for(int ii = 1; ii < ncs - 1; ii++)
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", xx[ii],
  	    (yy[ii + 1] - yy[ii - 1]) / (xx[ii + 1] - xx[ii - 1]), getCubicSpline1stDifferential1D(xx[ii], ncs, xx, yy, y2),
  	    (yy[ii + 1] + yy[ii - 1] - 2.0 * yy[ii]) / ((xx[ii + 1] - xx[ii - 1]) * (xx[ii + 1] - xx[ii - 1])), getCubicSpline2ndDifferential1D(xx[ii], ncs, xx, y2));
  free(xx);
  free(yy);
  free(bp);
  free(y2);
  exit(0);
#endif//TEST_CUBIC_SPLINE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* return the derived density profile */
  //-----------------------------------------------------------------------
  getDensityProfile(num, rad, rho, dr1, dr2, prf);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  releaseArray(rad, rho, psi, dr1, dr2);
  free(enc);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
