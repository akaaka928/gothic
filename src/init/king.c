/**
 * @file king.c
 *
 * @brief Source code for King model (lowered isothermal model)
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/06/02 (Fri)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "macro.h"
#include "constants.h"

#include "magi.h"
#include "profile.h"
#include "king.h"


extern const real newton;

static const double inv3 = 1.0 / 3.0;
static const double inv6 = 1.0 / 6.0;

static const double convergence = 1.0e-4;
/* static const double convergence = 1.0e-5; */
static const double extreme     = 0.25;


/**
 * @fn allocateArray
 *
 * @brief Memory allocation for King model.
 *
 * @param (num) number of data points
 * @return (rad) array for position
 * @return (rho) array for density
 * @return (psi) array for potential
 * @return (dr1) array for 1st derivative of density w.r.t. radius
 * @return (dr2) array for 2nd derivative of density w.r.t. radius
 */
static inline void allocateArray(const int num, double **rad, double **rho, double **psi, double **dr1, double **dr2)
{
  __NOTE__("%s\n", "start");

  *rad = (double *)malloc(sizeof(double) * num);  if( *rad == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rad\n");  }
  *rho = (double *)malloc(sizeof(double) * num);  if( *rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rho\n");  }
  *psi = (double *)malloc(sizeof(double) * num);  if( *psi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate psi\n");  }
  *dr1 = (double *)malloc(sizeof(double) * num);  if( *dr1 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dr1\n");  }
  *dr2 = (double *)malloc(sizeof(double) * num);  if( *dr2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dr2\n");  }

  __NOTE__("%s\n", "end");
}

/**
 * @fn enlargeArray
 *
 * @brief Memory allocation for King model.
 *
 * @param (num) number of data points
 * @return (rad) array for position
 * @return (rho) array for density
 * @return (psi) array for potential
 * @return (dr1) array for 1st derivative of density w.r.t. radius
 * @return (dr2) array for 2nd derivative of density w.r.t. radius
 */
static inline void enlargeArray(const int num, double **rad, double **rho, double **psi, double **dr1, double **dr2)
{
  __NOTE__("%s\n", "start");

  *rad = realloc(*rad, sizeof(double) * num);  if( *rad == NULL ){    __KILL__(stderr, "ERROR: failure to rescale rad\n");  }
  *rho = realloc(*rho, sizeof(double) * num);  if( *rho == NULL ){    __KILL__(stderr, "ERROR: failure to rescale rho\n");  }
  *psi = realloc(*psi, sizeof(double) * num);  if( *psi == NULL ){    __KILL__(stderr, "ERROR: failure to rescale psi\n");  }
  *dr1 = realloc(*dr1, sizeof(double) * num);  if( *dr1 == NULL ){    __KILL__(stderr, "ERROR: failure to rescale dr1\n");  }
  *dr2 = realloc(*dr2, sizeof(double) * num);  if( *dr2 == NULL ){    __KILL__(stderr, "ERROR: failure to rescale dr2\n");  }

  __NOTE__("%s\n", "end");
}

/**
 * @fn releaseArray
 *
 * @brief Memory deallocation for King model.
 *
 * @param (rad) array for position
 * @param (rho) array for density
 * @param (psi) array for potential
 * @param (dr1) array for 1st derivative of density w.r.t. radius
 * @param (dr2) array for 2nd derivative of density w.r.t. radius
 */
static inline void  releaseArray(               double  *rad, double  *rho, double  *psi, double  *dr1, double  *dr2)
{
  __NOTE__("%s\n", "start");

  free(rad);
  free(rho);
  free(psi);
  free(dr1);
  free(dr2);

  __NOTE__("%s\n", "end");
}


/**
 * @fn kingErrFunc
 *
 * @brief Subfunction to solve Poisson's equation for King sphere.
 *
 * @param (psi) non-dimensional potential
 * @param (sqrt_psi) square root of psi
 * @return exp(psi) * erf(sqrt(psi))
 */
static inline double kingErrFunc(const double psi, const double sqrt_psi)
{
  return (exp(psi) * erf(sqrt_psi));
}

/**
 * @fn kingFunc0
 *
 * @brief Subfunction to solve Poisson's equation for King sphere.
 *
 * @param (sqrt_psi) square root of non-dimensional potential
 * @param (kingExpErrFunc) result of kingErrFunc
 * @return kingExpErrFunc - 2 * sqrt(psi / pi)
 */
static inline double kingFunc0(const double sqrt_psi, const double kingExpErrFunc)
{
  return (kingExpErrFunc - M_2_SQRTPI * sqrt_psi);
}

/**
 * @fn getKingDensity
 *
 * @brief Subfunction to solve Poisson's equation for King sphere.
 *
 * @param (psi) non-dimensional potential
 * @param (sqrt_psi) square root of psi
 * @param (kingExpErrFunc) result of kingErrFunc
 * @param (rho1) scale density
 * @return volume density at given non-dimensional potential
 */
static inline double getKingDensity(const double psi, const double sqrt_psi, const double kingExpErrFunc, const double rho1)
{
  return (rho1 * (kingExpErrFunc - M_2_SQRTPI * sqrt_psi * (1.0 + 2.0 * psi * inv3)));
}

/**
 * @fn calcKingDensity
 *
 * @brief Get density of King sphere.
 *
 * @param (psi) non-dimensional potential
 * @param (rho1) scale density
 * @return volume density at given non-dimensional potential
 */
static inline double calcKingDensity(double psi, double rho1)
{
  const double sqrt_psi = sqrt(psi);
  return (getKingDensity(psi, sqrt_psi, kingErrFunc(psi, sqrt_psi), rho1));
}


/**
 * @def kingFunc1(rad, uu, yy)
 *
 * @brief Wrapper for _kingFunc1(rad, yy)
 */
#define kingFunc1(rad, uu, yy) _kingFunc1(rad, yy)
/**
 * @fn _kingFunc1
 *
 * @brief Subfunction to solve Poisson's equation for King sphere.
 *
 * @param (rad) position of data point
 * @param (yy) r^2 * dW/dr
 * @return dW/dr
 */
static inline double _kingFunc1(const double rad, const double yy)
{
  return ((rad > DBL_EPSILON) ? (yy / (rad * rad)) : (0.0));
}


/**
 * @def kingFunc2(rad, psi, dWdr, rho0inv, rho1)
 *
 * @brief Wrapper for _kingFunc2(rad, psi, rho0inv, rho1)
 */
#define kingFunc2(rad, psi, dWdr, rho0inv, rho1) _kingFunc2(rad, psi, rho0inv, rho1)
/**
 * @fn _kingFunc2
 *
 * @brief Subfunction to solve Poisson's equation for King sphere.
 *
 * @param (rad) position of data point
 * @param (psi) non-dimensional potential at rad
 * @param (rho0inv) 1 / rho0
 * @param (rho1) rho1
 * @return right hand side of the Poisson's equation, see Equation (B3) in Miki & Umemura (in preparation)
 */
static inline double _kingFunc2(const double rad, const double psi, const double rho0inv, const double rho1)
{
  return  (-9.0 * rad * rad * rho0inv * calcKingDensity(psi, rho1));
}


/**
 * @fn rungeKutta4thForKing
 *
 * @brief Runge--Kutta integrator in 4th order.
 *
 * @param (rr) position of data point
 * @return (uu) W
 * @return (yy) r^2 * dW/dr
 * @param (hh) width of the grid
 * @param (rho0inv) 1 / rho0
 * @param (rho1) rho1
 *
 * @sa _kingFunc1
 * @sa _kingFunc2
 */
static inline void rungeKutta4thForKing(const double rr, double *uu, double *yy, const double hh, const double rho0inv, const double rho1)
{
  const double k1u = kingFunc1(rr,            *uu,                  *yy                                );
  const double k1y = kingFunc2(rr,            *uu,                  *yy,                  rho0inv, rho1);
  const double k2u = kingFunc1(rr + 0.5 * hh, *uu + 0.5 * hh * k1u, *yy + 0.5 * hh * k1y               );
  const double k2y = kingFunc2(rr + 0.5 * hh, *uu + 0.5 * hh * k1u, *yy + 0.5 * hh * k1y, rho0inv, rho1);
  const double k3u = kingFunc1(rr + 0.5 * hh, *uu + 0.5 * hh * k2u, *yy + 0.5 * hh * k2y               );
  const double k3y = kingFunc2(rr + 0.5 * hh, *uu + 0.5 * hh * k2u, *yy + 0.5 * hh * k2y, rho0inv, rho1);
  const double k4u = kingFunc1(rr +       hh, *uu +       hh * k3u, *yy +       hh * k3y               );
  const double k4y = kingFunc2(rr +       hh, *uu +       hh * k3u, *yy +       hh * k3y, rho0inv, rho1);

  *uu += (hh * inv6) * (k1u + 2.0 * k2u + 2.0 * k3u + k4u);
  *yy += (hh * inv6) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
}


/**
 * @fn solvePoissonEqOfKingDF
 *
 * @brief Solver of Poisson's equation for King sphere.
 *
 * @param (W0) dimensionless King parameter at the center
 * @return (Rad) position of data point
 * @return (Psi) non-dimensional potential
 * @return (Rho) density
 * @return (Dr1) drho/dr
 * @return (Dr2) d^2rho/dr^2
 * @return (num) number of data points
 * @return (rem) number of unused data points
 *
 * @sa calcKingDensity
 * @sa rungeKutta4thForKing
 * @sa kingErrFunc
 * @sa getKingDensity
 * @sa enlargeArray
 * @sa kingFunc0
 */
static inline void solvePoissonEqOfKingDF
(const double W0,
#ifdef  KING_CENTRAL_CUSP
 const double dWdx_0,
#endif//KING_CENTRAL_CUSP
 double **Rad, double **Psi, double **Rho, double **Dr1, double **Dr2, int *num, int *rem)
{
  __NOTE__("%s\n", "start");


  double *rad, *psi, *rho, *dr1, *dr2;

  /** set initial condition */
  double hh = MINRAD * 0.25;
  double hnew = hh;

  /** assumption is rcal := r / r0, r0 = 1.0, sigma = 1.0 */
  /** see Eq (4.106) in Galactic Dynamics, second edtion, page 305 */
  /** rho0 = (real)2.25 * sigma * sigma / (pi * newton * r0 * r0); */
  const double rho0 = 2.25 * M_1_PI;
  const double rho0inv = 1.0 / rho0;

  /** uu stores the dimensionless King parameter W, yy stores r^2 * dW/dr */
  /** W is the dimensionless King parameter Psi / sigma^2 */
  double uu = W0;
  /** boundary condition is dW/dr(r = 0) = 0: assume core profile at the center */
#ifdef  KING_CENTRAL_CUSP
  double yy = hh * hh * dWdx_0;
#else///KING_CENTRAL_CUSP
  double yy = 0.0;
#endif//KING_CENTRAL_CUSP

  const double rho1 = rho0 / calcKingDensity(W0, 1.0);

  int ii = 0;
  rad = *Rad;  rad[ii] = 0.0;
  psi = *Psi;  psi[ii] = uu;
  rho = *Rho;  rho[ii] = calcKingDensity(psi[ii], rho1);
#ifdef  KING_CENTRAL_CUSP
  dr1 = *Dr1;  dr1[ii] = rho1 * kingFunc0(sqrt(W0), kingErrFunc(W0, sqrt(W0))) * dWdx_0;
#else///KING_CENTRAL_CUSP
  dr1 = *Dr1;  dr1[ii] = 0.0;
#endif//KING_CENTRAL_CUSP
  dr2 = *Dr2;  dr2[ii] = 0.0;



  double rold = rho[ii];

  /** solve Poisson Equation using 4th-order Runge-Kutta method */
#ifdef  KING_PROGRESS_REPORT_ON
  fprintf(stdout, "# Poisson solver for King model start\n");
  fflush(stdout);
#endif//KING_PROGRESS_REPORT_ON
  while( true ){
    hh = hnew;
    const double uold = uu;
    const double yold = yy;
    const double uoldinv = 1.0 / uold;
    const double yoldinv = 1.0 / (1.0e-30 + fabs(yold));
    const double roldinv = 1.0 / rold;

    double dens, sqrtW, kingExpErrFunc;
    while( true ){
      uu = uold;
      yy = yold;
      rungeKutta4thForKing(rad[ii], &uu, &yy, hh, rho0inv, rho1);

      sqrtW = sqrt(uu);
      kingExpErrFunc = kingErrFunc(uu, sqrtW);
      dens = getKingDensity(uu, sqrtW, kingExpErrFunc, rho1);

      /** convergence tests */
      const double udiff = fabs(uu - uold) * uoldinv;
      const double ydiff = fabs(yy - yold) * yoldinv;
      double  diff = (udiff > ydiff) ? udiff : ydiff;
      const double rdiff = (convergence * dens > DBL_EPSILON * rho[0]) ? (fabs(dens - rold) * roldinv) : (0.0);
      if( rdiff > diff )
	diff = rdiff;
      if( diff < convergence ){
	if( diff < extreme * convergence )	  hnew = hh * 2.0;
	break;
      }/* if( diff < convergence ){ */
      else{
	hh *= 0.5;
	hnew = hh;
      }/* else{ */
    }/* while( true ){ */

    if( *rem == 0 ){
      enlargeArray((*num) + NADD_KING, Rad, Rho, Psi, Dr1, Dr2);
      rad = *Rad;
      psi = *Psi;
      rho = *Rho;
      dr1 = *Dr1;
      dr2 = *Dr2;
      *rem += NADD_KING;
    }/* if( *rem == 0 ){ */
    ii++;    *num += 1;    *rem -= 1;
    rad[ii] = rad[ii - 1] + hh;

    if( (uu > DBL_EPSILON) && (dens > DBL_EPSILON) ){
      psi[ii] = uu;
      rho[ii] = dens;
      rold = dens;

      const double tmp = kingFunc0(sqrtW, kingExpErrFunc);
      const double rinv = 1.0 / rad[ii];
      const double dW_dr = yy * rinv * rinv;
      dr1[ii] = rho1 *  tmp * dW_dr;
      const double d2W_dr2 = -9.0 * dens * rho0inv - 2.0 * dW_dr * rinv;
      dr2[ii] = rho1 * (tmp * d2W_dr2 + kingExpErrFunc * dW_dr * dW_dr);

#if 0
      if( ii < 32 )
	fprintf(stderr, "%e\t%e\t%e\t%e\n", rad[ii], uu, yy, dW_dr);
#endif
    }/* if( (uu > DBL_EPSILON) && (dens > DBL_EPSILON) ){ */
    else{
      psi[ii] = 0.0;
      rho[ii] = 0.0;
      dr1[ii] = 0.0;
      dr2[ii] = 0.0;
      break;
    }/* else{ */

#ifdef  KING_PROGRESS_REPORT_ON
    if( (ii % KING_PROGRESS_REPORT_ON) == 0 ){
      fprintf(stdout, "# Poisson solver for King model: %d steps finished: rad = %e, rho = %e, W = %e\n", ii, rad[ii], rho[ii], psi[ii]);
      fflush(stdout);
    }/* if( (ii % KING_PROGRESS_REPORT_ON) == 0 ){ */
#endif//KING_PROGRESS_REPORT_ON
  }/* while( true ){ */

#ifdef  KING_PROGRESS_REPORT_ON
  fprintf(stdout, "# Poisson solver for King model finish: %d elements used: rt / r0 = %e, c = %e\n#\n#\n", *num, rad[*num - 1], log10(rad[*num - 1]));
  fflush(stdout);
#endif//KING_PROGRESS_REPORT_ON

#ifdef  KING_PROGRESS_REPORT_ON
  fprintf(stdout, "# King model: rho1 is %le w/o scaling\n", rho1);
  fflush(stdout);
#endif//KING_PROGRESS_REPORT_ON


  __NOTE__("%s\n", "end");
}


/**
 * @fn rescaleKingSphere
 *
 * @brief Rescale the King profile.
 *
 * @param (Mtot) total mass of the King sphere
 * @param (r0) King radius
 * @return (rt) tidal radius
 * @param (num) number of data points
 * @return (rad) position of data point
 * @return (rho) density
 * @return (dr1) drho/dr
 * @return (dr2) d^2rho/dr^2
 * @return (enc) M(r)
 */
static void rescaleKingSphere(const double Mtot, const double r0, double *rt, const int num, double *rad, double *rho, double *dr1, double *dr2, double *enc)
{
  __NOTE__("%s\n", "start");


  double MencCalc = 2.0 * M_PI;/**< 4pi / 2 */
  enc[0] = MencCalc * rad[1] * rad[1] * rad[1] * (rho[0] + rho[1]) * inv3;
  for(int ii = 1; ii < num - 1; ii++)
    enc[ii] = enc[ii - 1] + MencCalc * rad[ii] * rad[ii + 1] * (rad[ii + 1] - rad[ii]) * (rho[ii] + rho[ii + 1]);
  MencCalc = enc[num - 1] = enc[num - 2];
  double r0Calc = 1.0;


  /** evaluate scaling factors */
  const double lengthUnit = r0   /   r0Calc;
  const double   massUnit = Mtot / MencCalc;

  const double    rhoUnit = massUnit / (lengthUnit * lengthUnit * lengthUnit);

  const double  drho_drUnit  = rhoUnit / lengthUnit;
  const double d2rho_dr2Unit = rhoUnit / (lengthUnit * lengthUnit);

  *rt = rad[num - 1] * lengthUnit;

#pragma omp parallel for
  for(int ii = 0; ii < num; ii++){
    rad[ii] *=    lengthUnit;
    rho[ii] *=       rhoUnit;
    enc[ii] *=      massUnit;
    dr1[ii] *=   drho_drUnit;
    dr2[ii] *= d2rho_dr2Unit;
  }/* for(int ii = 0; ii < num; ii++){ */

#ifdef  KING_PROGRESS_REPORT_ON
  fprintf(stdout, "# King model:   lengthUnit is %le\n", lengthUnit);
  fprintf(stdout, "# King model:     massUnit is %le\n",   massUnit);
  fprintf(stdout, "# King model:      rhoUnit is %le\n",    rhoUnit);
  fprintf(stdout, "# King model: velocityUnit is %le\n", lengthUnit * sqrt((double)newton * rhoUnit));
  fprintf(stdout, "#\n#\n");
  fflush(stdout);
#endif//KING_PROGRESS_REPORT_ON


  __NOTE__("%s\n", "end");
}


/**
 * @fn findIdx
 *
 * @brief Find a data element in the given array corresponding to the given value.
 *
 * @param (rad) radius
 * @param (prf) radial profile of the King sphere
 * @return the corresponding index to the given radius
 */
static inline int findIdx(const double rad, profile *prf)
{
  int ll =           0;
  int rr = NRADBIN - 1;

  if( rad < prf[ll].rad + DBL_EPSILON ){    return (ll    );  }
  if( rad > prf[rr].rad - DBL_EPSILON ){    return (rr - 1);  }

  while( true ){
    const uint cc = ((uint)ll + (uint)rr) >> 1;

    if( (prf[cc].rad - rad) * (prf[ll].rad - rad) <= 0.0 )      rr = (int)cc;
    else                                                        ll = (int)cc;

    if( (1 + ll) == rr )
      return (ll);
  }/* while( true ){ */
}

/**
 * @fn bisection
 *
 * @brief Execute bisection.
 *
 * @param (val) the target value
 * @param (num) number of data points
 * @param (tab) array contains data points
 * @return (ratio) parameter for linear interpolation
 * @return lower index of the corresponding data point
 */
static inline int bisection(const double val, const int num, double tab[], double *ratio)
{
  int ll = 0;
  int rr = num - 1;

  /** prohibit extraporation */
  if( val < tab[ll] + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( val > tab[rr] - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }

  while( true ){
    const uint cc = ((uint)(ll + rr)) >> 1;

    if( (tab[cc] - val) * (tab[ll] - val) <= 0.0)      rr = (int)cc;
    else                                               ll = (int)cc;

    if( (1 + ll) == rr ){
      *ratio = (val - tab[ll]) / (tab[1 + ll] - tab[ll]);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
  }/* while( true ){ */
}


/**
 * @fn getDensityProfile
 *
 * @brief Get density profile of the King model.
 *
 * @param (num) number of data points
 * @param (rad) position of data point
 * @param (rho) density
 * @param (dr1) drho/dr
 * @param (dr2) d^2rho/dr^2
 * @return (prf) radial profile of the component
 *
 * @sa findIdx
 * @sa bisection
 */
static inline void getDensityProfile(const int num, double *rad, double *rho, double *dr1, double *dr2, profile *prf)
{
  __NOTE__("%s\n", "start");


  /** calculate spherical averaged profile of density, enclosed mass and potential (internal part) */
  const int head = findIdx(rad[      0], prf);
  const int tail = findIdx(rad[num - 1], prf);

#pragma omp parallel for
  for(int ii = 0; ii < head; ii++){
    prf[ii].rho       = rho[0];
    prf[ii].drho_dr   = dr1[0];
    prf[ii].d2rho_dr2 = dr2[0];
  }/* for(int ii = 0; ii < head; ii++){ */

  /** interpolate: enc --> prf */
#pragma omp parallel for
  for(int ii = head; ii < tail; ii++){
    double alpha;
    const int idx = bisection(prf[ii].rad, num, rad, &alpha);
    prf[ii].rho       = (1.0 - alpha) * rho[idx] + alpha * rho[1 + idx];
    prf[ii].drho_dr   = (1.0 - alpha) * dr1[idx] + alpha * dr1[1 + idx];
    prf[ii].d2rho_dr2 = (1.0 - alpha) * dr2[idx] + alpha * dr2[1 + idx];
  }/* for(int ii = head; ii < tail; ii++){ */

  /** fill the total enclosed mass in the outer region */
#pragma omp parallel for
  for(int ii = tail; ii < NRADBIN; ii++){
    prf[ii].rho       = 0.0;
    prf[ii].drho_dr   = 0.0;
    prf[ii].d2rho_dr2 = 0.0;
  }/* for(int ii = tail; ii < NRADBIN; ii++){ */

#if 0
  /* for(int ii = 0; ii < num; ii += 1024) */
  for(int ii = 0; ii < num; ii++)
    fprintf(stderr, "%e\t%e\t%e\t%e\n", rad[ii], rho[ii], dr1[ii], dr2[ii]);
  fflush(NULL);
  exit(0);
#endif

  __NOTE__("%s\n", "end");
}

/**
 * @fn setDensityProfileKing
 *
 * @brief Get radial profile of the King model.
 *
 * @return (prf) radial profile of the component
 * @param (cfg) physical properties of the component
 *
 * @sa allocateArray
 * @sa solvePoissonEqOfKingDF
 * @sa rescaleKingSphere
 * @sa getDensityProfile
 * @sa releaseArray
 */
void setDensityProfileKing(profile *prf, profile_cfg *cfg)
{
  __NOTE__("%s\n", "start");


  const double W0   = cfg->king_W0;
  const double Mtot = cfg->Mtot;
  const double r0   = cfg->rs;

  /** memory allocation */
  int num = 0;
  int rem = NADD_KING;
  double *rad, *rho, *psi, *dr1, *dr2;
  allocateArray(rem, &rad, &rho, &psi, &dr1, &dr2);


  /** derive density profile of the King sphere */
  solvePoissonEqOfKingDF(W0,
#ifdef  KING_CENTRAL_CUSP
			 cfg->king_dWdx_0,
#endif//KING_CENTRAL_CUSP
			 &rad, &psi, &rho, &dr1, &dr2, &num, &rem);
  double *enc;  enc = (double *)malloc(sizeof(double) * num);  if( enc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate enc\n");  }
  rescaleKingSphere(Mtot, r0, &cfg->king_rt, num, rad, rho, dr1, dr2, enc);

  cfg->king_c = log10(cfg->king_rt / cfg->rs);

  /** return the derived density profile */
  getDensityProfile(num, rad, rho, dr1, dr2, prf);

  /** memory deallocation */
  releaseArray(rad, rho, psi, dr1, dr2);
  free(enc);


  __NOTE__("%s\n", "end");
}
