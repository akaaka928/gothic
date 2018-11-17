/**
 * @file eddington.c
 *
 * @brief Source code for Eddington's formula
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/11/13 (Tue)
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
#include <assert.h>

#include "macro.h"
#include "constants.h"

#include "magi.h"
#include "profile.h"
#include "eddington.h"

#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
#include "spline.h"
#else///ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
extern double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF

extern const real newton;


#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF

/**
 * @fn get_d2rho_dPsi2
 *
 * @brief Calculate second derivative of the mass density w.r.t. the relative potential.
 *
 * @param (skind) number of components
 * @param (prf) radial profile of the component
 * @return (ret) second derivative of the mass density w.r.t. the relative potential
 */
static inline void get_d2rho_dPsi2(const int skind, profile **prf, double *ret
#ifdef  USE_OSIPKOV_MERRITT_METHOD
				   , const double ra2inv
#endif//USE_OSIPKOV_MERRITT_METHOD
)
{
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const int jj = NRADBIN - 1 - ii;

    const double rad = prf[0][jj].rad;
    const double enc = prf[0][jj].enc_tot;
    const double rho = prf[0][jj].rho_tot;
    const double r2 = rad * rad;
    const double fac = 2.0 * (enc - 2.0 * M_PI * rho * r2 * rad) / (rad * enc);

    double common = r2 / enc;
    common *= common;

#ifdef  USE_OSIPKOV_MERRITT_METHOD
    const double fac1 = 2.0 * ra2inv;
    const double fac2 = 1.0 + r2 * ra2inv;
#endif//USE_OSIPKOV_MERRITT_METHOD

    for(int kk = 0; kk < skind; kk++){
      const double  drho_dr  = prf[kk][jj]. drho_dr;
      const double d2rho_dr2 = prf[kk][jj].d2rho_dr2;

#ifndef USE_OSIPKOV_MERRITT_METHOD
      ret[ii + kk * NRADBIN] = (d2rho_dr2 + fac * drho_dr) * common;
#else///USE_OSIPKOV_MERRITT_METHOD
      const double rhoi = fac1 * prf[kk][jj].rho;
      ret[ii + kk * NRADBIN] = (rhoi + 2.0 * rad * drho_dr + fac2 * d2rho_dr2 + fac * (rad * rhoi + fac2 * drho_dr)) * common;
#endif//USE_OSIPKOV_MERRITT_METHOD
    }/* for(int kk = 0; kk < skind; kk++){ */

  }/* for(int ii = 0; ii < NRADBIN; ii++){ */
}


/**
 * @fn get_integrand
 *
 * @brief Calculate integrand for Eddington formula.
 *
 * @param (ene) relative energy
 * @param (psi) relative potential
 * @param (skind) number of components
 * @param (xx) position of data points (psi)
 * @param (yy) value of data points (d2rho_dpsi2)
 * @param (y2) coefficients in cubic spline interpolation
 * @return (ret) integrand for Eddington formula
 */
static inline void get_integrand(const double ene, const double psi, const int skind, double * restrict xx, double * restrict yy, double * restrict y2, double * restrict ret)
{
  const double coe = 1.0 / sqrt(ene - psi);

  for(int ii = 0; ii < skind; ii++)
    ret[ii] = coe * getCubicSpline1D(psi, NRADBIN, xx, &yy[ii * NRADBIN], &y2[ii * NRADBIN]);
}


/**
 * @fn get_DEformula
 *
 * @brief Calculate integrand in double exponential formula.
 *
 * @param (tt) value of integration variable
 * @return (ret) integrand for Eddington formula
 * @param (max_pls_min) psi_max + psi_min
 * @param (max_mns_min) psi_max - psi_min
 * @param (xx) position of data points (psi)
 * @param (yy) value of data points (d2rho_dpsi2)
 * @param (y2) coefficients in cubic spline interpolation
 */
static inline void get_DEformula(const double tt, const int skind, double ret[restrict], const double max_pls_min, const double max_mns_min, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double sinh_t = M_PI_2 * sinh(tt);

  const double cosh_t = cosh(sinh_t);
  const double common = cosh(tt) * exp(0.5 * sinh_t) * sqrt(cosh_t) / (cosh_t * cosh_t);

  /* since psi_max is around unity and psi_min is around 0.01 to 0.1, round-off error would not be significant when calculating psi */
#if 0
  const double psi_max = 0.5 * (max_pls_min + max_mns_min);
  const double psi_min = 0.5 * (max_pls_min - max_mns_min);
  const double psi = 0.5 * (psi_min * exp(-sinh_t) + psi_max * exp(sinh_t)) / cosh_t;
#else
  const double psi = 0.5 * (max_pls_min + max_mns_min * tanh(sinh_t));
#endif

  for(int kk = 0; kk < skind; kk++)
    ret[kk] += common * getCubicSpline1D(psi, NRADBIN, xx, &yy[kk * NRADBIN], &y2[kk * NRADBIN]);
}


static inline void update_trapezoidal(const double hh, const double tmin, const double tmax, const int skind, double sum[restrict], double tmp[restrict], const double max_pls_min, const double max_mns_min, double * restrict xx, double * restrict yy, double * restrict y2)
{
  /** initialization */
  for(int kk = 0; kk < skind; kk++)
    tmp[kk] = 0.0;
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    get_DEformula(tt, skind, tmp, max_pls_min, max_mns_min, xx, yy, y2);
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  for(int kk = 0; kk < skind; kk++)
    sum[kk] = 0.5 * sum[kk] + hh * tmp[kk];
}


static inline void set_domain_boundary(const double hh, double * restrict tmin, double * restrict tmax, const int skind, double sum[restrict], double ffp[restrict], double ff0[restrict], double fft[restrict], const double max_pls_min, const double max_mns_min, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  for(int kk = 0; kk < skind; kk++)    ffp[kk] = 0.0;
  get_DEformula(tt, skind, ffp, max_pls_min, max_mns_min, xx, yy, y2);
  for(int kk = 0; kk < skind; kk++){
    const double tmp = ffp[kk];
    ff0[kk] = tmp;
    sum[kk] = tmp * hh;
  }/* for(int kk = 0; kk < skind; kk++){ */


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    for(int kk = 0; kk < skind; kk++){
      fft[kk] = ffp[kk];
      ffp[kk] = 0.0;
    }/* for(int kk = 0; kk < skind; kk++){ */

    tt += hh;    boundary = tt;
    get_DEformula(tt, skind, ffp, max_pls_min, max_mns_min, xx, yy, y2);

    damp = -1.0;
    for(int kk = 0; kk < skind; kk++){
      const double tmp = ffp[kk];
      sum[kk] += hh * tmp;

      damp = fmax(damp, fabs(fft[kk]) + fabs(tmp));
    }/* for(int kk = 0; kk < skind; kk++){ */
  }/* while( (damp > converge) && (boundary < maximum) ){ */
  *tmax = boundary;


  /** determine lower boundary */
  for(int kk = 0; kk < skind; kk++)    ffp[kk] = ff0[kk];
  tt = 0.0;
  boundary = 0.0;
  damp = 1.0;
  while( (damp > converge) && (boundary > -maximum) ){
    for(int kk = 0; kk < skind; kk++){
      fft[kk] = ffp[kk];
      ffp[kk] = 0.0;
    }/* for(int kk = 0; kk < skind; kk++){ */

    tt -= hh;    boundary = tt;
    get_DEformula(tt, skind, ffp, max_pls_min, max_mns_min, xx, yy, y2);

    damp = -1.0;
    for(int kk = 0; kk < skind; kk++){
      const double tmp = ffp[kk];
      sum[kk] += hh * tmp;

      damp = fmax(damp, fabs(fft[kk]) + fabs(tmp));
    }/* for(int kk = 0; kk < skind; kk++){ */
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;

  /* __NOTE__("*tmin = %e, *tmax = %e\n", *tmin, *tmax); */
}


static inline void integrate_DEformula(const int skind, double sum[restrict], double ffp[restrict], double ff0[restrict], double fft[restrict], const double psi_min, const double psi_max, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double criteria_abs = 1.0e-12;
  /* const double criteria_rel = 1.0e-10; */
  const double criteria_rel = 1.0e-8;
  /* const double criteria_rel = 1.0e-6; */
  /* const double criteria_rel = 1.0e-5; */
  /* const double criteria_rel = 1.0e-4; */

  const double max_pls_min = psi_max + psi_min;
  const double max_mns_min = psi_max - psi_min;

  double hh = 1.0;
  double tmin, tmax;
  set_domain_boundary(hh, &tmin, &tmax, skind, sum, ffp, ff0, fft, max_pls_min, max_mns_min, xx, yy, y2);


  while( true ){
    for(int kk = 0; kk < skind; kk++)
      ff0[kk] = sum[kk];

    hh *= 0.5;
    update_trapezoidal(hh, tmin, tmax, skind, sum, fft, max_pls_min, max_mns_min, xx, yy, y2);

    bool converge = true;
    for(int kk = 0; kk < skind; kk++)
      if( converge ){
	if( fabs(sum[kk]) > DBL_EPSILON ){
	  if( fabs(1.0 - ff0[kk] / sum[kk]) > criteria_rel )
	    converge = false;
	}
	else
	  if( fabs(sum[kk] - ff0[kk]) > criteria_abs )
	    converge = false;
      }/* if( converge ){ */

    if( converge )
      break;
  }/* while( true ){ */


  __NOTE__("tmin = %e, tmax = %e, hh = %e\n", tmin, tmax, hh);
}


/**
 * @fn integrateEddingtonFormula
 *
 * @brief Derive distribution function.
 *
 * @param (kind) number of spherical symmetric components
 * @param (prf) radial profile of the component
 * @return (fene) distribution function
 *
 * @sa getEddingtonFormula
 */
void integrateEddingtonFormula(const int skind, profile **prf, profile_cfg *cfg, dist_func **fene)
{
  __NOTE__("%s\n", "start");


  int iout = 0;
  for(int kk = 0; kk < skind; kk++)
    iout = (iout < cfg[kk].iout) ? cfg[kk].iout : iout;

  /** set integration range */
  const double Emax = prf[0][   0].psi_tot;
#ifndef USE_OSIPKOV_MERRITT_METHOD
  const double Emin = prf[0][iout].psi_tot;
#else///USE_OSIPKOV_MERRITT_METHOD
  /**# note: if the below implementation leads nan or inf, then set Emin as prf[0][iout].psi_tot; */
  const double Emin = 0.0;/**< Q = \mathcal{E} - L^2 / (2 r_a^2) */
#endif//USE_OSIPKOV_MERRITT_METHOD
  const double Ebin = (Emax - Emin) / (double)(NENEBIN - 1);

  /** memory allocation for cubic spline interpolation */
  double *xx;  xx = (double *)malloc(sizeof(double)         * NRADBIN);  if( xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xx\n");  }
  double *yy;  yy = (double *)malloc(sizeof(double) * skind * NRADBIN);  if( yy == NULL ){    __KILL__(stderr, "ERROR: failure to allocate yy\n");  }
  double *bp;  bp = (double *)malloc(sizeof(double) * skind * NRADBIN);  if( bp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate bp\n");  }
  double *y2;  y2 = (double *)malloc(sizeof(double) * skind * NRADBIN);  if( y2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate y2\n");  }


  /** preparation of data table */
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++)
    xx[ii] = prf[0][NRADBIN - 1 - ii].psi_tot;

  get_d2rho_dPsi2(skind, prf, yy
#ifdef  USE_OSIPKOV_MERRITT_METHOD
    , cfg[kk].ra2inv
#endif//USE_OSIPKOV_MERRITT_METHOD
    );

  /** execute cubic spline interpolation */
#pragma omp parallel for
  for(int kk = 0; kk < skind; kk++)
    genCubicSpline1D(NRADBIN, xx, &yy[kk * NRADBIN], &bp[kk * NRADBIN], NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, &y2[kk * NRADBIN]);


  /** calculate overall coefficient */
  double common = 1.0 / CAST_R2D(newton);
  common *= common * 0.125 * M_1_PI;


  const double Ecut = prf[0][iout].psi_tot;
  const int icut = (int)floor((Ecut - Emin) / Ebin);

  __NOTE__("iout = %d, Ecut = %e, icut = %d, Nbin = %d\n", iout, Ecut, icut, NENEBIN);

  for(int kk = 0; kk < skind; kk++)
#pragma omp parallel for
    for(int ii = 0; ii < icut + 1; ii++){
      fene[kk][ii].ene = CAST_D2R(Emin + Ebin * (double)ii);
      fene[kk][ii].val = ZERO;
    }/* for(int ii = 0; ii < icut + 1; ii++){ */

#pragma omp parallel
  {
    double sum[NKIND_MAX], ffp[NKIND_MAX], ff0[NKIND_MAX], fft[NKIND_MAX];

#pragma omp for schedule(dynamic, 4)
    for(int ii = icut + 1; ii < NENEBIN; ii++){
      const double ene = Emin + Ebin * (double)ii;
      integrate_DEformula(skind, sum, ffp, ff0, fft, Emin, ene, xx, yy, y2);

      const double factor = sqrt(ene - Emin) * common;
      assert( fpclassify(factor) == FP_NORMAL );
      for(int kk = 0; kk < skind; kk++){
	fene[kk][ii].ene = CAST_D2R(ene);
	fene[kk][ii].val = CAST_D2R(fmax(factor * sum[kk], 0.0));
      }/* for(int kk = 0; kk < skind; kk++){ */
    }/* for(int ii = icut + 1; ii < NENEBIN; ii++){ */
  }


  /** memory deallocation for cubic spline interpolation */
  free(xx);
  free(yy);
  free(bp);
  free(y2);


  __NOTE__("%s\n", "end");
}


#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
/**
 * @fn get_DEformula
 *
 * @brief Calculate integrand in double exponential formula.
 *
 * @param (tt) value of integration variable
 * @param (skind) number of components
 * @return (ret) integrand for Eddington formula
 * @param (max_pls_min) psi_max + psi_min
 * @param (max_mns_min) psi_max - psi_min
 * @param (xx) position of data points (psi)
 * @param (yy) value of data points (d2rho_dpsi2)
 * @param (y2) coefficients in cubic spline interpolation
 */
static inline void get_DEformula_vdisp(const double tt, const int skind, double v2f[restrict], double v4f[restrict], const double vesc2, const double psi, dist_func **df, const double Emin, const double invEbin)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double cosh_t = cosh(sinh_t);
  const double inv_cosh2_t = 1.0 / (cosh_t * cosh_t);

  const double v2 = 0.25 * vesc2 * exp(2.0 * sinh_t) * inv_cosh2_t;
  const double ene = psi - 0.5 * v2;

  const double common = cosh(tt) * inv_cosh2_t * v2;

  for(int kk = 0; kk < skind; kk++){
    const double val = common * getDF(ene, df[kk], Emin, invEbin);
    v2f[kk] += val;
    v4f[kk] += val * v2;
  }/* for(int kk = 0; kk < skind; kk++){ */
}


static inline void update_trapezoidal_vdisp(const double hh, const double tmin, const double tmax, const int skind, double v2f[restrict], double v4f[restrict], const double vesc2, const double psi, dist_func **df, const double Emin, const double invEbin, double ss2[restrict], double ss4[restrict])
{
  /** initialization */
  double tt = tmin + hh;
  for(int kk = 0; kk < skind; kk++){
    ss2[kk] = 0.0;
    ss4[kk] = 0.0;
  }/* for(int kk = 0; kk < skind; kk++){ */

  /** employ mid-point rule */
  while( tt < tmax ){
    get_DEformula_vdisp(tt, skind, ss2, ss4, vesc2, psi, df, Emin, invEbin);

    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  for(int kk = 0; kk < skind; kk++){
    v2f[kk] = 0.5 * v2f[kk] + hh * ss2[kk];
    v4f[kk] = 0.5 * v4f[kk] + hh * ss4[kk];
  }/* for(int kk = 0; kk < skind; kk++){ */
}


static inline void set_domain_boundary_vdisp(const double hh, double * restrict tmin, double * restrict tmax, const int skind, double v2f[restrict], double v4f[restrict], const double vesc2, const double psi, dist_func **df, const double Emin, const double invEbin, double fp2[restrict], double f02[restrict], double ft2[restrict], double fp4[restrict], double f04[restrict], double ft4[restrict])
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  for(int kk = 0; kk < skind; kk++)
    fp2[kk] = fp4[kk] = 0.0;
  get_DEformula_vdisp(tt, skind, fp2, fp4, vesc2, psi, df, Emin, invEbin);
  for(int kk = 0; kk < skind; kk++){
    const double pp2 = fp2[kk];
    const double pp4 = fp4[kk];
    f02[kk] = pp2;    v2f[kk] = hh * pp2;
    f04[kk] = pp4;    v4f[kk] = hh * pp4;
  }/* for(int kk = 0; kk < skind; kk++){ */


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    for(int kk = 0; kk < skind; kk++){
      ft2[kk] = fp2[kk];      fp2[kk] = 0.0;
      ft4[kk] = fp4[kk];      fp4[kk] = 0.0;
    }/* for(int kk = 0; kk < skind; kk++){ */

    tt += hh;    boundary = tt;
    get_DEformula_vdisp(tt, skind, fp2, fp4, vesc2, psi, df, Emin, invEbin);

    damp = -1.0;
    for(int kk = 0; kk < skind; kk++){
      const double tmp2 = fp2[kk];
      const double tmp4 = fp4[kk];

      v2f[kk] += hh * tmp2;      damp = fmax(damp, fabs(ft2[kk]) + fabs(tmp2));
      v4f[kk] += hh * tmp4;      damp = fmax(damp, fabs(ft4[kk]) + fabs(tmp4));
    }/* for(int kk = 0; kk < skind; kk++){ */
  }/* while( (damp > converge) && (boundary < maximum) ){ */
  *tmax = boundary;


  /** determine lower boundary */
  for(int kk = 0; kk < skind; kk++){
    fp2[kk] = f02[kk];
    fp4[kk] = f04[kk];
  }/* for(int kk = 0; kk < skind; kk++){ */
  tt = 0.0;
  boundary = 0.0;
  damp = 1.0;
  while( (damp > converge) && (boundary > -maximum) ){
    for(int kk = 0; kk < skind; kk++){
      ft2[kk] = fp2[kk];      fp2[kk] = 0.0;
      ft4[kk] = fp4[kk];      fp4[kk] = 0.0;
    }/* for(int kk = 0; kk < skind; kk++){ */

    tt -= hh;    boundary = tt;
    get_DEformula_vdisp(tt, skind, fp2, fp4, vesc2, psi, df, Emin, invEbin);

    damp = -1.0;
    for(int kk = 0; kk < skind; kk++){
      const double tmp2 = fp2[kk];
      const double tmp4 = fp4[kk];

      v2f[kk] += hh * tmp2;      damp = fmax(damp, fabs(ft2[kk]) + fabs(tmp2));
      v4f[kk] += hh * tmp4;      damp = fmax(damp, fabs(ft4[kk]) + fabs(tmp4));
    }/* for(int kk = 0; kk < skind; kk++){ */
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;
}


static inline void integrate_DEformula_vdisp(const int skind, double v2f[restrict], double v4f[restrict], const double vesc2, const double psi, dist_func **df, const double Emin, const double invEbin, double fp2[restrict], double f02[restrict], double ft2[restrict], double fp4[restrict], double f04[restrict], double ft4[restrict])
{
  const double criteria_abs = 1.0e-12;
  /* const double criteria_rel = 1.0e-8; */
  const double criteria_rel = 1.0e-6;
  /* const double criteria_rel = 1.0e-5; */
  /* const double criteria_rel = 1.0e-4; */

  double hh = 1.0;
  double tmin, tmax;
  set_domain_boundary_vdisp(hh, &tmin, &tmax, skind, v2f, v4f, vesc2, psi, df, Emin, invEbin, fp2, f02, ft2, fp4, f04, ft4);

  while( true ){
    for(int kk = 0; kk < skind; kk++){
      ft2[kk] = v2f[kk];
      ft4[kk] = v4f[kk];
    }/* for(int kk = 0; kk < skind; kk++){ */

    hh *= 0.5;
    update_trapezoidal_vdisp(hh, tmin, tmax, skind, v2f, v4f, vesc2, psi, df, Emin, invEbin, fp2, fp4);

    bool converge = true;
    for(int kk = 0; kk < skind; kk++)
      if( converge ){
	if( fabs(v2f[kk]) > DBL_EPSILON ){
	  if( fabs(1.0 - ft2[kk] / v2f[kk]) > criteria_rel )
	    converge = false;
	}
	else
	  if( fabs(v2f[kk] - ft2[kk]) > criteria_abs )
	    converge = false;

	if( converge ){
	  if( fabs(v4f[kk]) > DBL_EPSILON ){
	    if( fabs(1.0 - ft4[kk] / v4f[kk]) > criteria_rel )
	      converge = false;
	  }
	  else
	    if( fabs(v4f[kk] - ft4[kk]) > criteria_abs )
	      converge = false;
	}/* if( converge ){ */
      }/* if( converge ){ */

    if( converge )
      break;
  }/* while( true ){ */
}


/**
 * @fn calcVelocityDispersionProfile
 *
 * @brief Calculate velocity dispersion profile.
 *
 * @param (skind) number of spherical symmetric components
 * @return (prf) radial profile of the components
 * @param (fene) distribution function of the components
 *
 * @sa gaussQuadVelocity
 */
void calcVelocityDispersionProfile(const int skind, profile **prf, profile_cfg *cfg, dist_func **df)
{
  __NOTE__("%s\n", "start");


  const double    Emin = df[0][          0].ene;
  const double    Emax = df[0][NENEBIN - 1].ene;
  const double invEbin = (double)(NENEBIN - 1) / (Emax - Emin);


  int iout = 0;
  for(int kk = 0; kk < skind; kk++)
    iout = (iout < cfg[kk].iout) ? cfg[kk].iout : iout;


#pragma omp parallel
  {
    double v2f[NKIND_MAX], fp2[NKIND_MAX], f02[NKIND_MAX], ft2[NKIND_MAX];
    double v4f[NKIND_MAX], fp4[NKIND_MAX], f04[NKIND_MAX], ft4[NKIND_MAX];

#pragma omp for schedule(dynamic, 4) nowait
    for(int ii = 0; ii < iout + 1; ii += SKIP_INTERVAL_FOR_VELOCITY_DISPERSION){
      const double psi = prf[0][ii].psi_tot;
      const double vesc2 = 2.0 * (psi - Emin);

      /* call double exponential formula */
      integrate_DEformula_vdisp(skind, v2f, v4f, vesc2, psi, df, Emin, invEbin, fp2, f02, ft2, fp4, f04, ft4);

      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].sigr = sqrt(v4f[kk] / (DBL_MIN + 3.0 * v2f[kk]));
	prf[kk][ii].v2f  = v2f[kk];
	prf[kk][ii].v4f  = v4f[kk];
      }/* for(int kk = 0; kk < skind; kk++){ */
    }/* for(int ii = 0; ii < iout + 1; ii += SKIP_INTERVAL_FOR_VELOCITY_DISPERSION){ */

    for(int kk = 0; kk < skind; kk++)
#pragma omp for schedule(dynamic, 4) nowait
      for(int ii = iout + 1; ii < NRADBIN; ii++){
	prf[kk][ii].sigr = 0.0;
	prf[kk][ii].v2f  = 0.0;
	prf[kk][ii].v4f  = 0.0;
      }/* for(int ii = iout + 1; ii < NRADBIN; ii++){ */
  }


#   if  SKIP_INTERVAL_FOR_VELOCITY_DISPERSION != 1
#pragma omp parallel
  {

    for(int kk = 0; kk < skind; kk++){
      const int iout = cfg[kk].iout;

#pragma omp for schedule(dynamic, 4) nowait
      for(int ii = 0; ii < iout + 1; ii += SKIP_INTERVAL_FOR_VELOCITY_DISPERSION){
	const double sigr0 = prf[kk][ii                                        ].sigr;
	const double  v2f0 = prf[kk][ii                                        ]. v2f;
	const double  v4f0 = prf[kk][ii                                        ]. v4f;
	const double sigr1 = prf[kk][ii + SKIP_INTERVAL_FOR_VELOCITY_DISPERSION].sigr;
	const double  v2f1 = prf[kk][ii + SKIP_INTERVAL_FOR_VELOCITY_DISPERSION]. v2f;
	const double  v4f1 = prf[kk][ii + SKIP_INTERVAL_FOR_VELOCITY_DISPERSION]. v4f;

	double sigr_slope = (sigr1 - sigr0) / (double)SKIP_INTERVAL_FOR_VELOCITY_DISPERSION;
	double  v2f_slope = ( v2f1 -  v2f0) / (double)SKIP_INTERVAL_FOR_VELOCITY_DISPERSION;
	double  v4f_slope = ( v4f1 -  v4f0) / (double)SKIP_INTERVAL_FOR_VELOCITY_DISPERSION;

	for(int jj = 1; jj < SKIP_INTERVAL_FOR_VELOCITY_DISPERSION; jj++){
	  prf[kk][ii + jj].sigr = sigr0 + sigr_slope * (double)jj;
	  prf[kk][ii + jj]. v2f =  v2f0 +  v2f_slope * (double)jj;
	  prf[kk][ii + jj]. v4f =  v4f0 +  v4f_slope * (double)jj;
	}/* for(int jj = 1; jj < SKIP_INTERVAL_FOR_VELOCITY_DISPERSION; jj++){ */
      }/* for(int ii = 0; ii < iout + 1; ii += SKIP_INTERVAL_FOR_VELOCITY_DISPERSION){ */

    }/* for(int kk = 0; kk < skind; kk++){ */

  }
#endif//SKIP_INTERVAL_FOR_VELOCITY_DISPERSION != 1


  __NOTE__("%s\n", "end");
}
#endif//MAKE_VELOCITY_DISPERSION_PROFILE

#else///ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF

/**
 * @fn findIdx
 *
 * @brief Find a data element in the given array corresponding to the given value.
 *
 * @param (psi) relative potential
 * @param (prf) radial profile of the component
 * @return (ll) the corresponding lower index
 * @return (rr) the corresponding upper index
 */
static inline void findIdx(const double psi, profile *prf, int *ll, int *rr)
{
  bool bisection = true;
  *ll =           0;
  *rr = NRADBIN - 1;

  if( bisection == true )    if( fabs(prf[*ll].psi_tot - psi) / psi < DBL_EPSILON ){      bisection = false;      *rr = (*ll) + 1;    }
  if( bisection == true )    if( fabs(prf[*rr].psi_tot - psi) / psi < DBL_EPSILON ){      bisection = false;      *ll = (*rr) - 1;    }

  while( bisection ){
    const uint cc = ((uint)(*ll) + (uint)(*rr)) >> 1;

    if( (prf[cc].psi_tot - psi) * (prf[*ll].psi_tot - psi) <= 0.0 )      *rr = (int)cc;
    else                                                                 *ll = (int)cc;

    if( (1 + (*ll)) == (*rr) )      break;
  }/* while( bisection ){ */
}


/**
 * @fn getEddingtonFormula
 *
 * @brief Calculate distribution function using the Eddington's formula.
 *
 * @param (ene) relative energy
 * @param (psi) relative potential
 * @param (kind) number of components
 * @param (prf) radial profile of the component
 * @return (val) distribution function of each component
 *
 * @sa findIdx
 */
static inline void getEddingtonFormula(const double ene, const double psi, const int kind, profile **prf, double *val
#ifdef  USE_OSIPKOV_MERRITT_METHOD
				       , profile_cfg *cfg
#endif//USE_OSIPKOV_MERRITT_METHOD
				       )
{
  int ll, rr;
  findIdx(psi, prf[0], &ll, &rr);

  /** based on linear interpolation */
  const double r_dr = (psi - prf[0][ll].psi_tot) / (prf[0][rr].psi_tot - prf[0][ll].psi_tot);
  const double rad = (1.0 - r_dr) * prf[0][ll].rad     + r_dr * prf[0][rr].rad;
  const double enc = (1.0 - r_dr) * prf[0][ll].enc_tot + r_dr * prf[0][rr].enc_tot;
#if 1
  const double rho = (1.0 - r_dr) * prf[0][ll].rho_tot + r_dr * prf[0][rr].rho_tot;
  const double fac = 2.0 * (enc - 2.0 * M_PI * rho * rad * rad * rad) / (DBL_MIN + rad * enc);
#else
  const double fl = 2.0 * (prf[0][ll].enc_tot - 2.0 * M_PI * prf[0][ll].rho_tot * prf[0][ll].rad * prf[0][ll].rad * prf[0][ll].rad) / (DBL_MIN + prf[0][ll].rad * prf[0][ll].enc_tot);
  const double fr = 2.0 * (prf[0][rr].enc_tot - 2.0 * M_PI * prf[0][rr].rho_tot * prf[0][rr].rad * prf[0][rr].rad * prf[0][rr].rad) / (DBL_MIN + prf[0][rr].rad * prf[0][rr].enc_tot);
  const double fac = (1.0 - r_dr) * fl + r_dr * fr;
#endif

  double common = rad * rad / enc;
  common *= common / sqrt(ene - psi);

  for(int kk = 0; kk < kind; kk++){
    /** based on linear interpolation */
    const double  drho_dr  = (1.0 - r_dr) * prf[kk][ll]. drho_dr  + r_dr * prf[kk][rr]. drho_dr;
    const double d2rho_dr2 = (1.0 - r_dr) * prf[kk][ll].d2rho_dr2 + r_dr * prf[kk][rr].d2rho_dr2;

#ifndef USE_OSIPKOV_MERRITT_METHOD
    val[kk] = common * (d2rho_dr2 + drho_dr * fac);
#else///USE_OSIPKOV_MERRITT_METHOD
    const double ra2inv = cfg[kk].ra2inv;
    const double fac1 = 2.0 * ra2inv;
    const double fac2 = 1.0 + rad * rad * ra2inv;

    const double rhoi = fac1 * ((1.0 - r_dr) * prf[kk][ll].rho + r_dr * prf[kk][rr].rho);
    val[kk] = (rhoi + 2.0 * rad * drho_dr + fac2 * d2rho_dr2 + fac * (rad * rhoi + fac2 * drho_dr)) * common;
#endif//USE_OSIPKOV_MERRITT_METHOD
  }/* for(int kk = 0; kk < kind; kk++){ */
}


/**
 * @fn gaussQuad1dEddington
 *
 * @brief Execute Gaussian quadrature to get distribution function.
 *
 * @param (num) number of data points
 * @param (xmin) minimum of interval of integration
 * @param (xmax) maximum of interval of integration
 * @param (kind) number of components
 * @param (prf) radial profile of the component
 * @return (sum) distribution function of each component
 * @param (fm) temporary array
 * @param (fp) temporary array
 *
 * @sa getEddingtonFormula
 */
static inline void gaussQuad1dEddington
(const int num, const double xmin, const double xmax,
 const int kind, profile **prf, double * restrict sum, double * restrict fm, double * restrict fp
#ifdef  USE_OSIPKOV_MERRITT_METHOD
 , profile_cfg *cfg
#endif//USE_OSIPKOV_MERRITT_METHOD
 )
{
  const double mns = 0.5 * (xmax - xmin);
  const double pls = 0.5 * (xmax + xmin);

  for(int kk = 0; kk < kind; kk++)    sum[kk] = 0.0;
  if( num & 1 ){
    const double ww = gsl_gaussQD_weight[(num >> 1)];
    const double xx = pls + mns * gsl_gaussQD_pos[(num >> 1)];
    getEddingtonFormula(xmax, xx, kind, prf, fm
#ifdef  USE_OSIPKOV_MERRITT_METHOD
			, cfg
#endif//USE_OSIPKOV_MERRITT_METHOD
			);
    for(int kk = 0; kk < kind; kk++)
      sum[kk] = ww * fm[kk];
  }/* if( num & 1 ){ */

  for(int ii = (num >> 1) - 1; ii >= 0; ii--){
    const double ww = gsl_gaussQD_weight[ii];
    const double xp = pls + mns * gsl_gaussQD_pos[ii];
    const double xm = pls - mns * gsl_gaussQD_pos[ii];
    getEddingtonFormula(xmax, xp, kind, prf, fp
#ifdef  USE_OSIPKOV_MERRITT_METHOD
			, cfg
#endif//USE_OSIPKOV_MERRITT_METHOD
			);
    getEddingtonFormula(xmax, xm, kind, prf, fm
#ifdef  USE_OSIPKOV_MERRITT_METHOD
			, cfg
#endif//USE_OSIPKOV_MERRITT_METHOD
			);

    for(int kk = 0; kk < kind; kk++)
      sum[kk] += ww * (fm[kk] + fp[kk]);
  }/* for(int ii = (num >> 1) - 1; ii >= 0; ii--){ */

  for(int kk = 0; kk < kind; kk++)
    sum[kk] *= mns;
}


/**
 * @fn integrateEddingtonFormula
 *
 * @brief Derive distribution function.
 *
 * @param (kind) number of spherical symmetric components
 * @param (prf) radial profile of the component
 * @return (fene) distribution function
 *
 * @sa gaussQuad1dEddington
 */
void integrateEddingtonFormula(const int skind, profile **prf
#ifdef  USE_OSIPKOV_MERRITT_METHOD
			       , profile_cfg *cfg
#endif//USE_OSIPKOV_MERRITT_METHOD
			       , dist_func **fene)
{
  __NOTE__("%s\n", "start");


  /** set integration range */
  const double Emax = prf[0][          0].psi_tot;
#ifndef USE_OSIPKOV_MERRITT_METHOD
  const double Emin = prf[0][NRADBIN - 1].psi_tot;
#else///USE_OSIPKOV_MERRITT_METHOD
  /**# note: if the below implementation leads nan or inf, then set Emin as prf[0][NRADBIN - 1].psi_tot; */
  const double Emin = 0.0;/**< Q = \mathcal{E} - L^2 / (2 r_a^2) */
#endif//USE_OSIPKOV_MERRITT_METHOD
  const double Ebin = (Emax - Emin) / (double)(NENEBIN - 1);

  double common = M_1_PI / CAST_R2D(newton);
  common *= common * 0.5 * M_SQRT1_2;


/* #ifdef  NDIVIDE_GAUSSQD */
/*   const int nsub = NENEBIN / NDIVIDE_GAUSSQD; */
/*   int head = 0; */
/*   double Ein = 0.0; */
/*   static double sub[NKIND_MAX]; */
/*   for(int ii = 0; ii < skind; ii++)    sub[ii] = 0.0; */

/*   for(int iter = 0; iter < NDIVIDE_GAUSSQD; iter++){ */
/*     const int tail = head + nsub; */
/* #pragma omp parallel for schedule(dynamic, 16) */
/*     for(int ii = head; ii < tail; ii++){ */
/*       double sum[NKIND_MAX], fm[NKIND_MAX], fp[NKIND_MAX]; */
/*       const double ene = Emin + Ebin * (double)ii; */

/*       gaussQuad1dEddington(NINTBIN, Ein, ene, skind, prf, sum, fm, fp); */

/*       for(int kk = 0; kk < skind; kk++){ */
/* 	sum[kk] += sub[kk]; */
/* 	fene[kk][ii].ene = CAST_D2R(ene); */
/* 	fene[kk][ii].val = CAST_D2R(common * fmax(sum[kk], 0.0)); */
/*       }/\* for(int kk = 0; kk < skind; kk++){ *\/ */
/*     }/\* for(int ii = head; ii < tail; ii++){ *\/ */

/*     head = tail; */
/*     for(int ii = 0; ii < skind; ii++) */
/*       sub[ii] = fene[ii][head - 1].val; */
/*     Ein = fene[0][head - 1].ene; */
/*   }/\* for(int iter = 0; iter < NDIVIDE_GAUSSQD; iter++){ *\/ */
/* #else///NDIVIDE_GAUSSQD */
#pragma omp parallel
  {
    double sum[NKIND_MAX], fm[NKIND_MAX], fp[NKIND_MAX];

#pragma omp for schedule(dynamic, 4)
    for(int ii = 0; ii < NENEBIN; ii++){
      const double ene = Emin + Ebin * (double)ii;
      gaussQuad1dEddington(NINTBIN, 0.0, ene, skind, prf, sum, fm, fp
#ifdef  USE_OSIPKOV_MERRITT_METHOD
			   , cfg
#endif//USE_OSIPKOV_MERRITT_METHOD
			   );

      for(int kk = 0; kk < skind; kk++){
	fene[kk][ii].ene = CAST_D2R(ene);
	fene[kk][ii].val = CAST_D2R(fmax(common * sum[kk], 0.0));
      }/* for(int kk = 0; kk < skind; kk++){ */
    }/* for(int ii = 0; ii < NENEBIN; ii++){ */
  }
/* #endif//NDIVIDE_GAUSSQD */


  __NOTE__("%s\n", "end");
}


#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
/**
 * @fn gaussQuadVelocity
 *
 * @brief Gaussian quadrature for calculating velocity dispersion
 *
 * @param (num) number of data points
 * @param (psi) relative potential
 * @param (vmin) minimum of velocity
 * @param (vmax) maximum of velocity
 * @param (df) distribution function of the components
 * @param (Emin) minimum value of the relative energy per unit mass
 * @param (invEbin) inverse of energy bin in distribution function
 * @return (v2f) integrated value of v^2 * df(v)
 * @return (v4f) integrated value of v^4 * df(v)
 *
 * @sa getDF
 */
#ifndef USE_OSIPKOV_MERRITT_METHOD
void gaussQuadVelocity(const int num, const double psi, const double vmin, const double vmax, dist_func *df, const double Emin, const double invEbin, double * restrict v2f, double * restrict v4f);
void gaussQuadVelocity(const int num, const double psi, const double vmin, const double vmax, dist_func *df, const double Emin, const double invEbin, double * restrict v2f, double * restrict v4f)
{
  const double mns = 0.5 * (vmax - vmin);
  const double pls = 0.5 * (vmax + vmin);

  (*v2f) = (*v4f) = 0.0;

  if( num & 1 ){
    const double weight =             gsl_gaussQD_weight[(num >> 1)];
    const double  value = pls + mns * gsl_gaussQD_pos   [(num >> 1)];

    const double v2 = value * value;    const double v2df = v2 * getDF(psi - 0.5 * v2, df, Emin, invEbin);

    (*v2f) = weight * v2df;
    (*v4f) = weight * v2df * v2;
  }/* if( num & 1 ){ */

  for(int ii = (num >> 1) - 1; ii >= 0; ii--){
    const double weight = gsl_gaussQD_weight[ii];
    const double vp = pls + mns * gsl_gaussQD_pos[ii];
    const double vm = pls - mns * gsl_gaussQD_pos[ii];

    const double vp2 = vp * vp;    const double vp2df = vp2 * getDF(psi - 0.5 * vp2, df, Emin, invEbin);
    const double vm2 = vm * vm;    const double vm2df = vm2 * getDF(psi - 0.5 * vm2, df, Emin, invEbin);

    (*v2f) += weight * (      vp2df +       vm2df);
    (*v4f) += weight * (vp2 * vp2df + vm2 * vm2df);
  }/* for(int ii = (num >> 1) - 1; ii >= 0; ii--){ */

  (*v2f) *= mns;
  (*v4f) *= mns;
}
#else///USE_OSIPKOV_MERRITT_METHOD
void gaussQuadVelocity(const int num, const double psi, const double Emin, const double Emax, dist_func *df, const double Emin, const double invEbin, double * restrict v2f, double * restrict v4f)
{
  const double mns = 0.5 * (Emax - Emin);
  const double pls = 0.5 * (Emax + Emin);

  (*v2f) = (*v4f) = 0.0;

  if( num & 1 ){
    const double weight =             gsl_gaussQD_weight[(num >> 1)];
    const double  value = pls + mns * gsl_gaussQD_pos   [(num >> 1)];

    const double v2 = value * value;    const double v2df = v2 * getDF(psi - 0.5 * v2, df, Emin, invEbin);

    (*v2f) = weight * v2df;
    (*v4f) = weight * v2df * v2;
  }/* if( num & 1 ){ */

  for(int ii = (num >> 1) - 1; ii >= 0; ii--){
    const double weight = gsl_gaussQD_weight[ii];
    const double vp = pls + mns * gsl_gaussQD_pos[ii];
    const double vm = pls - mns * gsl_gaussQD_pos[ii];

    const double vp2 = vp * vp;    const double vp2df = vp2 * getDF(psi - 0.5 * vp2, df, Emin, invEbin);
    const double vm2 = vm * vm;    const double vm2df = vm2 * getDF(psi - 0.5 * vm2, df, Emin, invEbin);

    (*v2f) += weight * (      vp2df +       vm2df);
    (*v4f) += weight * (vp2 * vp2df + vm2 * vm2df);
  }/* for(int ii = (num >> 1) - 1; ii >= 0; ii--){ */

  (*v2f) *= mns;
  (*v4f) *= mns;
}
#endif//USE_OSIPKOV_MERRITT_METHOD


/**
 * @fn calcVelocityDispersionProfile
 *
 * @brief Calculate velocity dispersion profile.
 *
 * @param (skind) number of spherical symmetric components
 * @return (prf) radial profile of the components
 * @param (fene) distribution function of the components
 *
 * @sa gaussQuadVelocity
 */
void calcVelocityDispersionProfile(const int skind, profile **prf, dist_func **df)
{
  __NOTE__("%s\n", "start");


  static double Emin[NKIND_MAX], Emax[NKIND_MAX], invEbin[NKIND_MAX];
  for(int kk = 0; kk < skind; kk++){
    Emin[kk] = df[kk][          0].ene;
    Emax[kk] = df[kk][NENEBIN - 1].ene;
    invEbin[kk] = (double)(NENEBIN - 1) / (Emax[kk] - Emin[kk]);
  }/* for(int kk = 0; kk < skind; kk++){ */

#pragma omp parallel for schedule(dynamic, 4)
  for(int ii = 0; ii < NRADBIN; ii += SKIP_INTERVAL_FOR_VELOCITY_DISPERSION){
#ifndef USE_OSIPKOV_MERRITT_METHOD
    /** initialization */
    for(int kk = 0; kk < skind; kk++)
      prf[kk][ii].v2f = prf[kk][ii].v4f = 0.0;
#endif//USE_OSIPKOV_MERRITT_METHOD

    const double psi = prf[0][ii].psi_tot;
    for(int kk = 0; kk < skind; kk++){
#ifndef USE_OSIPKOV_MERRITT_METHOD
      /** numerical quadrature of v2f and v4f in vel [0, vesc] */
      double v2f, v4f;
      const double vesc = sqrt(2.0 * (psi - Emin[kk]));
      gaussQuadVelocity(NINTBIN, psi, 0.0, vesc, df[kk], Emin[kk], invEbin[kk], &v2f, &v4f);

      prf[kk][ii].sigr = sqrt(v4f / (3.0 * v2f));
      prf[kk][ii].v2f  = v2f;
      prf[kk][ii].v4f  = v4f;
#else///USE_OSIPKOV_MERRITT_METHOD
      /** numerical quadrature of sigr in Q [Emin, Psi] */


      prf[kk][ii].sigr = hoge;


#endif//USE_OSIPKOV_MERRITT_METHOD
    }/* for(int kk = 0; kk < skind; kk++){ */
  }/* for(int ii = 0; ii < NRADBIN; ii += SKIP_INTERVAL_FOR_VELOCITY_DISPERSION){ */


#   if  SKIP_INTERVAL_FOR_VELOCITY_DISPERSION != 1
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN - SKIP_INTERVAL_FOR_VELOCITY_DISPERSION; ii += SKIP_INTERVAL_FOR_VELOCITY_DISPERSION)
    for(int kk = 0; kk < skind; kk++){
      const double sigr0 = prf[kk][ii                                        ].sigr;
#ifndef USE_OSIPKOV_MERRITT_METHOD
      const double  v2f0 = prf[kk][ii                                        ]. v2f;
      const double  v4f0 = prf[kk][ii                                        ]. v4f;
#endif//USE_OSIPKOV_MERRITT_METHOD
      const double sigr1 = prf[kk][ii + SKIP_INTERVAL_FOR_VELOCITY_DISPERSION].sigr;
#ifndef USE_OSIPKOV_MERRITT_METHOD
      const double  v2f1 = prf[kk][ii + SKIP_INTERVAL_FOR_VELOCITY_DISPERSION]. v2f;
      const double  v4f1 = prf[kk][ii + SKIP_INTERVAL_FOR_VELOCITY_DISPERSION]. v4f;
#endif//USE_OSIPKOV_MERRITT_METHOD

      double sigr_slope = (sigr1 - sigr0) / (double)SKIP_INTERVAL_FOR_VELOCITY_DISPERSION;
#ifndef USE_OSIPKOV_MERRITT_METHOD
      double  v2f_slope = ( v2f1 -  v2f0) / (double)SKIP_INTERVAL_FOR_VELOCITY_DISPERSION;
      double  v4f_slope = ( v4f1 -  v4f0) / (double)SKIP_INTERVAL_FOR_VELOCITY_DISPERSION;
#endif//USE_OSIPKOV_MERRITT_METHOD

      for(int jj = 1; jj < SKIP_INTERVAL_FOR_VELOCITY_DISPERSION; jj++){
	prf[kk][ii + jj].sigr = sigr0 + sigr_slope * (double)jj;
#ifndef USE_OSIPKOV_MERRITT_METHOD
	prf[kk][ii + jj]. v2f =  v2f0 +  v2f_slope * (double)jj;
	prf[kk][ii + jj]. v4f =  v4f0 +  v4f_slope * (double)jj;
#endif//USE_OSIPKOV_MERRITT_METHOD
      }/* for(int jj = 1; jj < SKIP_INTERVAL_FOR_VELOCITY_DISPERSION; jj++){ */
    }/* for(int kk = 0; kk < skind; kk++){ */
#endif//SKIP_INTERVAL_FOR_VELOCITY_DISPERSION != 1
  __NOTE__("%s\n", "end");
}
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
