/**
 * @file profile.c
 *
 * @brief Source code for describing radial profile of spherical component(s)
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2019/09/19 (Thu)
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
#include "name.h"

#include "profile.h"

#ifdef  MAKE_COLUMN_DENSITY_PROFILE
#include "magi.h"
#endif//MAKE_COLUMN_DENSITY_PROFILE


#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
#include "spline.h"
#define NFIT_PROFILE (4)
#define NCAP_PROFILE (16)
#define NSPLINE_PROFILE (NRADBIN + NCAP_PROFILE)
#else///ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
extern double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE

#ifdef  ENABLE_GASEOUS_COMPONENT
#include "gas.h"
#endif//ENABLE_GASEOUS_COMPONENT

extern const real newton;
extern const double     mass_astro2com;
extern const double   length_astro2com;
extern const double velocity_astro2com;
#ifdef  ENABLE_GASEOUS_COMPONENT
extern const double temperature_astro2com;
#endif//ENABLE_GASEOUS_COMPONENT


/**
 * @fn setDensityProfilePlummer
 *
 * @brief Get radial profile of the Plummer model.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale radius of the component
 */
void setDensityProfilePlummer(profile *prf, const double rs)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad_rs = prf[ii].rad * rsinv;
    const double x2_pls_1 = 1.0 + rad_rs * rad_rs;
    const double     inv = 1.0 / x2_pls_1;
    const double sqrtinv = sqrt(inv);
    const double pow_m5_2 = inv * inv * sqrtinv;

    prf[ii].  rho     = pow_m5_2;
    prf[ii]. drho_dr  = - 5.0 * rad_rs * pow_m5_2 * inv * rsinv;
    prf[ii].d2rho_dr2 = 5.0 * pow_m5_2 * inv * (7.0 * rad_rs * rad_rs * inv - 1.0) * rsinv * rsinv;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setDensityProfileBurkert
 *
 * @brief Get radial profile of the Burkert model.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale radius of the component
 */
void setDensityProfileBurkert(profile *prf, const double rs)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad_rs      = prf[ii].rad * rsinv;
    const double rad_rs_p1   = rad_rs + 1.0;
    const double rad2_rs2    = rad_rs * rad_rs;
    const double rad2_rs2_p1 = 1.0 + rad2_rs2;
    const double inv = 1.0 / (rad_rs_p1 * rad2_rs2_p1);

    prf[ii].  rho     = inv;
    prf[ii]. drho_dr  = - (1.0 + rad_rs * (2.0 + 3.0 * rad_rs)) * inv * inv * rsinv;
    prf[ii].d2rho_dr2 = 4.0 * rad2_rs2 * (4.0 * rad_rs + 3.0 * rad2_rs2_p1) * inv * inv * inv * rsinv * rsinv;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setDensityProfileHernquist
 *
 * @brief Get radial profile of the Hernquist model.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale radius of the component
 */
void setDensityProfileHernquist(profile *prf, const double rs)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad_rs = prf[ii].rad * rsinv;
    const double rad_rs_p1 = rad_rs + 1.0;
    const double tmp = rad_rs * rad_rs_p1;
    const double inv = 1.0 / (tmp * rad_rs_p1);

    prf[ii].  rho     = 1.0 * inv * inv * tmp;
    prf[ii]. drho_dr  = - (1.0 + 4.0 * rad_rs) * inv * inv * rsinv;
    prf[ii].d2rho_dr2 = 2.0 * (1.0 + 5.0 * rad_rs * (1.0 + 2.0 * rad_rs)) * inv * inv * inv * rad_rs_p1 * rsinv * rsinv;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setDensityProfileNFW
 *
 * @brief Get radial profile of the NFW model.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale radius of the component
 */
void setDensityProfileNFW(profile *prf, const double rs)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad_rs = prf[ii].rad * rsinv;
    const double rad_rs_p1 = rad_rs + 1.0;
    const double tmp = rad_rs * rad_rs_p1;
    const double inv = 1.0 / (tmp * rad_rs_p1);

    prf[ii].  rho     = inv;
    prf[ii]. drho_dr  = - (1.0 + 3.0 * rad_rs) * inv * inv * rad_rs_p1 * rsinv;
    prf[ii].d2rho_dr2 = 2.0 * (2.0 * rad_rs * (3.0 * rad_rs + 2.0) + 1.0) * inv * inv * inv * rad_rs_p1 * rad_rs_p1 * rsinv * rsinv;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setDensityProfileMoore
 *
 * @brief Get radial profile of the Moore model.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale radius of the component
 */
void setDensityProfileMoore(profile *prf, const double rs)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad_rs   = prf[ii].rad * rsinv;
    const double     inv  = 1.0 / (rad_rs * (1.0 + rad_rs));
    const double sqrtinv  = sqrt(inv);
    const double pow_m5_2 = inv * inv * sqrtinv;

    prf[ii].  rho     =  1.0 * inv * sqrtinv;
    prf[ii]. drho_dr  = -1.5 * (1.0 + 2.0 * rad_rs) * pow_m5_2 * rsinv;
    prf[ii].d2rho_dr2 =  0.75 * (5.0 + 16.0 * rad_rs * (1.0 + rad_rs)) * pow_m5_2 * inv * rsinv * rsinv;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setDensityProfileEinasto
 *
 * @brief Get radial profile of the Einasto model.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale radius of the component
 * @param (alpha) parameter to determine the degree of transition steepness of the density slope
 */
void setDensityProfileEinasto(profile *prf, const double rs, const double alpha)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad_rs = prf[ii].rad * rsinv;
    const double xalpha = pow(rad_rs, alpha);
    const double xinv = 1.0 / rad_rs;

    const double rho = pow(M_E, (-2.0 / alpha) * (xalpha - 1.0));
    prf[ii].  rho     = rho;
    prf[ii]. drho_dr  = -2.0 * xalpha * xinv * rho * rsinv;
    prf[ii].d2rho_dr2 =  2.0 * xalpha * xinv * xinv * rho * rsinv * rsinv * (2.0 * xalpha + 1.0 - alpha);
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setDensityProfileAppKing
 *
 * @brief Get radial profile of the King model in empirical form.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale radius of the component
 * @param (rt) tidal radius of the component
 */
void setDensityProfileAppKing(profile *prf, const double rs, const double rt)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad_rs = prf[ii].rad * rsinv;
    const double CC = 1.0 / sqrt(1.0 + rt * rsinv * rt * rsinv);
    const double inv = 1.0 / (1.0 + rad_rs * rad_rs);
    const double sqrtinv = sqrt(inv);

    prf[ii].  rho     = (prf[ii].rad < rt) ? ((sqrtinv - CC) * (sqrtinv - CC)) : (0.0);
    prf[ii]. drho_dr  = (prf[ii].rad < rt) ? (-2.0 * (sqrtinv - CC) * rad_rs * inv * sqrtinv * rsinv) : (0.0);
    prf[ii].d2rho_dr2 = (prf[ii].rad < rt) ? ((2.0 * inv * inv * (4.0 * rad_rs * rad_rs * inv - 1.0) + 2.0 * CC * sqrtinv * inv * (1.0 - 3.0 * rad_rs * rad_rs * inv)) * rsinv * rsinv) : (0.0);
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setDensityProfileTwoPower
 *
 * @brief Get radial profile of the double power-law model.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale radius of the component
 * @param (alpha) inner power-law slope of the component
 * @param (beta) parameter to control the transition steepness from inner power-law slope to outer power-law slope
 * @param (gam) outer power-law slope of the component
 */
void setDensityProfileTwoPower(profile *prf, const double rs, const double alpha, const double beta, const double gam)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
  const double rsinv = 1.0 / rs;
  const double binv  = 1.0 / beta;
  const double agbi  = (alpha - gam) * binv;
  const double b1g2  = 1.0 + beta + 2.0 * gam;
  const double a1    = 1.0 + alpha;
  const double g1    = 1.0 + gam;
  const double b1    = 1.0 - beta;

#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double xx    = prf[ii].rad * rsinv;
    const double xinv  = 1.0 / xx;
    const double xa    = pow(xinv, alpha);
    const double xb    = pow(  xx, beta);
    const double xbp1  = 1.0 + xb;
    const double xbp1i = 1.0 / xbp1;
    const double base  = xa * pow(xbp1, agbi);

    prf[ii].  rho     =                  base;
    prf[ii]. drho_dr  = -rsinv         * base * xinv        * xbp1i         * (alpha                    + gam * xb);
    prf[ii].d2rho_dr2 =  rsinv * rsinv * base * xinv * xinv * xbp1i * xbp1i * (alpha * (a1 + b1g2 * xb) + gam * xb * (b1 + g1 * xb));
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setDensityProfileTriPower
 *
 * @brief Get radial profile of the triple power-law model.
 *
 * @return (prf) radial profile of the component
 * @param (rin) inner scale radius of the component
 * @param (rout) outer scale radius of the component
 * @param (alp) innermost power-law slope of the component
 * @param (bet) parameter to control the transition steepness from alp to gam
 * @param (gam) intermediate power-law slope of the component
 * @param (del) parameter to control the transition steepness from gam to eps
 * @param (eps) outermost power-law slope of the component
 */
void setDensityProfileTriPower(profile *prf, const double rin, const double rout, const double alp, const double bet, const double gam, const double del, const double eps)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
  /** At first, set alpha-beta-gamma profile */
  setDensityProfileTwoPower(prf, rin, alp, bet, gam);
  const double rsinv = 1.0 / rout;
  const double gmeod = (gam - eps) / del;

#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double yy   = prf[ii].rad * rsinv;
    const double yd   = pow(yy, del);
    const double ydp1 = 1.0 + yd;
    const double inv  = 1.0 / (yy * ydp1);

    const double base0 = pow(ydp1, gmeod);
    const double base1 = (gam - eps) * yd * rsinv * inv * base0;
    const double base2 =                    rsinv * inv * base1 * ((del - 1.0) + (gam - eps - 1.0) * yd);

    prf[ii].d2rho_dr2 = prf[ii].d2rho_dr2 * base0 + 2.0 * prf[ii].drho_dr * base1 + prf[ii].rho * base2;
    prf[ii]. drho_dr  = prf[ii]. drho_dr  * base0                                 + prf[ii].rho * base1;
    prf[ii].  rho     = prf[ii].  rho     * base0;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setDensityProfileAppLoweredEvans
 *
 * @brief Get radial profile of the approximated Lowered Evans model.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale radius of the component
 * @param (alpha) fitting parameter
 * @param (rc) fitting parameter
 * @param (beta) fitting parameter
 * @param (rt) fitting parameter
 * @param (invdelta) fitting parameter
 */
void setDensityProfileAppLoweredEvans(profile *prf, const double rs, const double alpha, const double rc, const double beta, const double rt, const double invdelta)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
  const double rsinv  = 1.0 / rs;
  const double rcinv  = 1.0 / rc;
  const double rc2 = rc * rc;
  const double  betax2 = 2.0 * beta;
  const double alphax2 = 2.0 * alpha;

#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad = prf[ii].rad;
    const double rad_rs = rad * rsinv;
    const double rad_rc = rad * rcinv;
    const double inv_1p_rad_rs  = 1.0 / (1.0 + rad_rs);
    const double inv_1p_rad_rc2 = 1.0 / (1.0 + rad_rc * rad_rc);
    const double inv_1p_rad_rs_alpha = pow(inv_1p_rad_rs , alpha);
    const double inv_1p_rad_rc2_beta = pow(inv_1p_rad_rc2, beta);
    const double rho = inv_1p_rad_rs_alpha * inv_1p_rad_rc2_beta;

    const double r2 = rad * rad;
    const double inv = 1.0 / ((rs + rad) * (rc2 + r2));
    const double val = alpha * (rc2 + r2) + betax2 * (rs + rad) * rad;
    const double drho = -inv * val * rho;

    const double sinv = 1.0 / (rs + rad);
    const double rinv = 1.0 / (rc2 + r2);
    const double d2rho = rho * ((1.0 + alpha) * alpha * sinv * sinv + betax2 * (rc2 + rad * ((1.0 + betax2) * rad + alphax2 * (rc2 + r2) * sinv)) * rinv * rinv);

    double smooth, dsmooth, d2smooth;
    if( rad < rt ){
      smooth   = 1.0;
      dsmooth  = 0.0;
      d2smooth = 0.0;
    }/* if( rad < rt ){ */
    else{
      smooth   = exp(-(rad - rt) * invdelta);
      dsmooth  = -invdelta            * smooth;
      d2smooth =  invdelta * invdelta * smooth;
    }/* else{ */

    prf[ii].  rho     = rho * smooth;
    prf[ii]. drho_dr  = rho * dsmooth + drho * smooth;
    prf[ii].d2rho_dr2 = rho * d2smooth + 2.0 * drho * dsmooth + d2rho * smooth;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn setContributionByCentralBH
 *
 * @brief Set central massive black hole.
 *
 * @return (prf) radial profile of the component
 * @param (cfg) physical properties of the component
 */
void setContributionByCentralBH(profile *prf, const profile_cfg cfg)
{
  __NOTE__("%s\n", "start");

  /** set density profile */
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    prf[ii].  rho     = 0.0;
    prf[ii]. drho_dr  = 0.0;
    prf[ii].d2rho_dr2 = 0.0;

    prf[ii].enc =                    cfg.Mtot;
    prf[ii].psi = CAST_R2D(newton) * cfg.Mtot / prf[ii].rad;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setIsothermalGas
 *
 * @brief Get radial profile of the isothermal gaseous component in hydrostatic.
 *
 * @return (prf) radial profile of the component
 * @param (Tgas) temperature of the gas
 * @param (mu) mean molecular weight of the gas
 * @param (rho0) scaling factor of the density
 */
void setIsothermalGas(profile *prf, const double Tgas, const double mu, const double rho0)
{
  __NOTE__("%s\n", "start");

  /** assumption: rho_tot, enc_tot, and psi_tot are already prepared */
  /* double rho_tot, enc_tot, psi_tot; */
  const double psi0 = prf[0].psi_tot;

  /** set density profile */
#if 1
  extern const double boltzmann, protonmass;
  const double coeff = (mu * protonmass) / (boltzmann * Tgas);
  /* __KILL__(stderr, "coeff = %e, psi0 = %e, exponent = %e, ret = %e\n", coeff, psi0, coeff * psi0, exp(coeff * psi0)); */
#else
  const double coeff = 1.0 / prf[0].psi_tot;
#endif

#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double psi = prf[ii].psi_tot;
    const double rho = rho0 * exp(-coeff * (psi0 - psi));/**< psi = -Phi(r) */

    /* const double rad = prf[ii].rad; */
    /* const double rinv = 1.0 / rad; */
    /* const double enc = prf[ii].enc_tot; */
    /* const double grav = newton * enc * rinv * rinv; */

    prf[ii].rho = rho;
    /* prf[ii]. drho_dr  = -coeff * grav * rho; */
    /* prf[ii].d2rho_dr2 =  coeff * grav * rho * (coeff * grav + 2.0 * (rinv - 2.0 * M_PI * rad * rad * prf[ii].rho_tot / enc)); */
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

#if 0
  for(int ii = 0; ii < NRADBIN; ii += (NRADBIN / 128))
    fprintf(stdout, "%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].psi_tot);
  __KILL__(stderr, "suspend for debugging\n");
#endif

  __NOTE__("%s\n", "end");
}


#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
/**
 * @fn leastSquaresMethod
 *
 * @brief Data fitting by the least squares method (assume power-law model).
 *
 * @param (num) number of data points
 * @param (xx) x
 * @param (yy) y = y(x)
 * @return (pp) the resultant power-law index
 * @return (bb) the resultant base
 */
static inline void leastSquaresMethod(const int num, double * restrict xx, double * restrict yy, double * restrict pp, double * restrict bb)
{
  double SS, Sx, Sy, Sxx, Sxy;
  SS = Sx = Sy = Sxx = Sxy = 0.0;

  for(int ii = 0; ii < num; ii++){
    const double logx = log10(xx[ii]);
    const double logy = log10(yy[ii]);
    SS  += 1.0;
    Sx  += logx;
    Sxx += logx * logx;
    Sy  +=        logy;
    Sxy += logx * logy;
  }/* for(int ii = 0; ii < num; ii++){ */

  *pp = (SS * Sxy - Sx * Sy) / (SS * Sxx - Sx * Sx);
  *bb = pow(10.0, (Sy - (*pp) * Sx) / SS);

#if 1
  /** tentative treatment for the case that y is almost zero and hence least squares method returns inf */
  if( isfinite(*bb) == 0 ){
    *pp = 1.0;
    *bb = 0.0;/* 0.5 * (yy[(num - 1) >> 1] + yy[num >> 1]);/\* this is the median of yy *\/ */
  }/* if( isfinite(*bb) == 0 ){ */
#endif
}


/**
 * @fn get_DEformula_internal
 *
 * @brief Calculate integrand in double exponential formula.
 *
 * @param (tt) value of integration variable
 * @param (rr) 
 * @param (xx) position of data points (psi)
 * @param (yy) value of data points (d2rho_dpsi2)
 * @param (y2) coefficients in cubic spline interpolation
 * @return integrand for 
 */
static inline double get_DEformula_internal(const double tt, const double rr, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double cosh_t = cosh(sinh_t);
  const double inv_cosh_t = 1.0 / cosh_t;

  const double one_pls_x = exp(sinh_t) * inv_cosh_t;

  return (cosh(tt) * inv_cosh_t * inv_cosh_t * one_pls_x * one_pls_x * getCubicSpline1D(0.5 * rr * one_pls_x, NSPLINE_PROFILE, xx, yy, y2));
}


/**
 * @fn get_DEformula_external
 *
 * @brief Calculate integrand in double exponential formula.
 *
 * @param (tt) value of integration variable
 * @param (rmin) 
 * @param (rmax) 
 * @param (xx) position of data points (psi)
 * @param (yy) value of data points (d2rho_dpsi2)
 * @param (y2) coefficients in cubic spline interpolation
 * @return integrand for 
 */
static inline double get_DEformula_external(const double tt, const double rmin, const double rmax, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double cosh_t = cosh(sinh_t);
  const double inv_cosh_t = 1.0 / cosh_t;

  const double one_pls_x = exp( sinh_t) * inv_cosh_t;
  const double one_mns_x = exp(-sinh_t) * inv_cosh_t;

  const double rr = 0.5 * (rmin * one_mns_x + rmax * one_pls_x);

  return (cosh(tt) * inv_cosh_t * inv_cosh_t * rr * getCubicSpline1D(rr, NSPLINE_PROFILE, xx, yy, y2));
}


static inline double update_trapezoidal_internal(const double hh, const double tmin, const double tmax, const double sum, const double rr, double * restrict xx, double * restrict yy, double * restrict y2)
{
  /** initialization */
  double sub = 0.0;
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    sub += get_DEformula_internal(tt, rr, xx, yy, y2);
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  return (0.5 * sum + hh * sub);
}


static inline double update_trapezoidal_external(const double hh, const double tmin, const double tmax, const double sum, const double rmin, const double rmax, double * restrict xx, double * restrict yy, double * restrict y2)
{
  /** initialization */
  double sub = 0.0;
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    sub += get_DEformula_external(tt, rmin, rmax, xx, yy, y2);
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  return (0.5 * sum + hh * sub);
}


static inline double set_domain_boundary_internal(const double rs, const double hh, double * restrict tmin, double * restrict tmax, const double rr, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  double fp = get_DEformula_internal(tt, rr, xx, yy, y2);
  const double f0 = fp;
  double sum = hh * fp;


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    const double ft = fp;

    tt += hh;    boundary = tt;
    fp = get_DEformula_internal(tt, rr, xx, yy, y2);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);
  }/* while( (damp > converge) && (boundary < maximum) ){ */
  *tmax = boundary;


  /** determine lower boundary */
  fp = f0;
  tt = 0.0;
  boundary = 0.0;
  damp = 1.0;
  while( (damp > converge) && (boundary > -maximum) ){
    const double ft = fp;

    tt -= hh;    boundary = tt;
    fp = get_DEformula_internal(tt, rr, xx, yy, y2);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);

    const double yy = M_PI_2 * sinh(tt);
    if( 0.5 * rr * exp(yy) / cosh(yy) > rs )
      damp = 1.0;
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;

  return (sum);
}


static inline double set_domain_boundary_external(const double hh, double * restrict tmin, double * restrict tmax, const double rmin, const double rmax, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  double fp = get_DEformula_external(tt, rmin, rmax, xx, yy, y2);
  const double f0 = fp;
  double sum = hh * fp;


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    const double ft = fp;

    tt += hh;    boundary = tt;
    fp = get_DEformula_external(tt, rmin, rmax, xx, yy, y2);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);

#if 0
    const double yy = M_PI_2 * sinh(tt);
    if( 0.5 * (rmin * exp(-yy) + rmax * exp(yy)) / cosh(yy) < 0.9 * rmax )
      damp = 1.0;
#endif
  }/* while( (damp > converge) && (boundary < maximum) ){ */
  *tmax = boundary;


  /** determine lower boundary */
  fp = f0;
  tt = 0.0;
  boundary = 0.0;
  damp = 1.0;
  while( (damp > converge) && (boundary > -maximum) ){
    const double ft = fp;

    tt -= hh;    boundary = tt;
    fp = get_DEformula_external(tt, rmin, rmax, xx, yy, y2);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);

    const double yy = M_PI_2 * sinh(tt);
    if( 0.5 * (rmin * exp(-yy) + rmax * exp(yy)) / cosh(yy) > 1.01 * rmin )
      damp = 1.0;
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;

  return (sum);
}


static inline double integrate_DEformula_internal(const double rr, const double rs, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double criteria_abs = 1.0e-12;
  /* const double criteria_rel = 1.0e-10; */
  const double criteria_rel = 1.0e-8;
  /* const double criteria_rel = 1.0e-7; */
  /* const double criteria_rel = 1.0e-6; */
  /* const double criteria_rel = 1.0e-5; */
  /* const double criteria_rel = 1.0e-4; */

  double hh = 1.0;
  double tmin, tmax;
  double sum = set_domain_boundary_internal(rs, hh, &tmin, &tmax, rr, xx, yy, y2);


  while( true ){
    const double f0 = sum;

    hh *= 0.5;
    sum = update_trapezoidal_internal(hh, tmin, tmax, sum, rr, xx, yy, y2);

    bool converge = true;
    if( fabs(sum) > DBL_EPSILON ){
      if( fabs(1.0 - f0 / sum) > criteria_rel )
	converge = false;
    }
    else
      if( fabs(sum - f0) > criteria_abs )
	converge = false;


    if( converge )
      break;
  }/* while( true ){ */

  return (rr * rr * rr * sum);
}


static inline double integrate_DEformula_external(const double rmin, const double rmax, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double criteria_abs = 1.0e-12;
  /* const double criteria_rel = 1.0e-10; */
  const double criteria_rel = 1.0e-8;
  /* const double criteria_rel = 1.0e-7; */
  /* const double criteria_rel = 1.0e-6; */
  /* const double criteria_rel = 1.0e-5; */
  /* const double criteria_rel = 1.0e-4; */

  double hh = 1.0;
  double tmin, tmax;
  double sum = set_domain_boundary_external(hh, &tmin, &tmax, rmin, rmax, xx, yy, y2);

  while( true ){
    const double f0 = sum;

    hh *= 0.5;
    sum = update_trapezoidal_external(hh, tmin, tmax, sum, rmin, rmax, xx, yy, y2);

    bool converge = true;
    if( fabs(sum) > DBL_EPSILON ){
      if( fabs(1.0 - f0 / sum) > criteria_rel )
	converge = false;
    }
    else
      if( fabs((sum) - f0) > criteria_abs )
	converge = false;


    if( converge )
      break;
  }/* while( true ){ */

  return ((rmax - rmin) * sum);
}
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE


/**
 * @fn integrateDensityProfile
 *
 * @brief Integrate the density profile.
 *
 * @return (prf) radial profile of the component
 * @param (logrbin) bin width in the logarithmic space
 */
void integrateDensityProfile(profile *prf, profile_cfg *cfg
#ifndef ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
			     , const double logrbin
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
			     )
{
  __NOTE__("%s\n", "start");


  const double Mtot = cfg->Mtot;
  const bool cutoff = cfg->cutoff;
  const double redge = cfg->rc;
  const double width = cfg->rc_width;

  /** set cutoff radius if necessary */
  if( cutoff == true ){
    const double invThick = 1.0 / width;
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
#ifdef  ERFC_SMOOTHING
      const double smooth_x = 0.5 * (prf[ii].rad - redge) * invThick;
      const double smooth = 0.5 * erfc(smooth_x);
      const double dhdx   = -0.25 * M_2_SQRTPI * invThick * exp(-smooth_x * smooth_x);
      const double d2hdx2 = -dhdx * invThick * smooth_x;
#else///ERFC_SMOOTHING
      const double smooth_tanh = tanh(0.5 * (prf[ii].rad - redge) * invThick);
      const double smooth = 0.5 * (1.0 - smooth_tanh);
      const double dhdx = -0.25 * invThick * (1.0 - smooth_tanh * smooth_tanh);
      const double d2hdx2 = -dhdx * invThick * smooth_tanh;
#endif//ERFC_SMOOTHING

      prf[ii].d2rho_dr2  = d2hdx2 * prf[ii].rho + 2.0 * dhdx * prf[ii].drho_dr + smooth * prf[ii].d2rho_dr2;
      prf[ii]. drho_dr   =  dhdx  * prf[ii].rho +     smooth * prf[ii].drho_dr;
      prf[ii].  rho     *= smooth;
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */
  }/* if( cutoff == true ){ */


#if 1
  int iout = NRADBIN - 1;
  for(int ii = NRADBIN - 3; ii >= 0; ii--)
    if( prf[ii].rho > DBL_MIN ){
      iout = ii + 2;
      break;
    }
  double rmax = prf[iout].rad;
#else
  int iout = NRADBIN - 1;
  double rmax = prf[NRADBIN - 1].rad;
#endif
  cfg->iout = iout;
  cfg->rmax = rmax;

/* #pragma omp parallel for */
/*   for(int ii = 0; ii < NRADBIN; ii++){ */
/*     prf[ii].enc = 0.0; */
/*     prf[ii].psi = 0.0; */
/*   } */

#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
  const double rs = cfg->rs;

/** memory allocation for cubic spline interpolation */
  double *xx;  xx = (double *)malloc(sizeof(double) * NSPLINE_PROFILE);  if( xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xx\n");  }
  double *bp;  bp = (double *)malloc(sizeof(double) * NSPLINE_PROFILE);  if( bp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate bp\n");  }
  double *yy;  yy = (double *)malloc(sizeof(double) * NSPLINE_PROFILE);  if( yy == NULL ){    __KILL__(stderr, "ERROR: failure to allocate yy\n");  }
  double *y2;  y2 = (double *)malloc(sizeof(double) * NSPLINE_PROFILE);  if( y2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate y2\n");  }

#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    xx[NCAP_PROFILE + ii] = prf[ii].rad;
    yy[NCAP_PROFILE + ii] = prf[ii].rho;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  double pp, bb;
  leastSquaresMethod(NFIT_PROFILE, &xx[NCAP_PROFILE], &yy[NCAP_PROFILE], &pp, &bb);
  double rbin = prf[0].rad / (double)NCAP_PROFILE;

  for(int ii = 1; ii < NCAP_PROFILE; ii++){
    xx[ii] = rbin * (double)ii;
    yy[ii] = bb * pow(xx[ii], pp);
  }/* for(int ii = 1; ii < NCAP_PROFILE; ii++){ */
  xx[0] = 0.0;
  yy[0] = yy[1];

#if 0
  for(int ii = 0; ii < NSPLINE_PROFILE; ii++)
    fprintf(stderr, "%e\t%e\n", xx[ii], yy[ii]);
  exit(0);
#endif


  /** execute cubic spline interpolation */
  genCubicSpline1D(NSPLINE_PROFILE, xx, yy, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, y2);

  const double M_PI_16 = 0.125 * M_PI_2;

#pragma omp parallel for schedule(dynamic, 8)
  for(int ii = 0; ii < iout + 1; ii++){
    const double rad = prf[ii].rad;
    const double enc = M_PI_16 * integrate_DEformula_internal(rad, rs, xx, yy, y2);
    const double ext = M_PI_4  * integrate_DEformula_external(rad, rmax, xx, yy, y2);

    prf[ii].enc = enc;
    prf[ii].psi = ext + enc / rad;
  }/* for(int ii = 0; ii < iout + 1; ii++){ */

#pragma omp parallel for
  for(int ii = iout + 1; ii < NRADBIN; ii++){
    prf[ii].enc = prf[iout].enc;
    prf[ii].psi = prf[ii].enc / prf[ii].rad;
  }/* for(int ii = iout + 1; ii < NRADBIN; ii++){ */

  free(xx);  free(bp);
  free(yy);  free(y2);

#else///ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE

  /** calculate enclosed mass and potential (internal part) */
  double Menc[2];
  Menc[0] = prf[0].rad * prf[0].rad * prf[0].rad * prf[0].rho;
  Menc[1] = prf[1].rad * prf[1].rad * prf[1].rad * prf[1].rho;
  prf[0].enc =              Menc[0] / 3.0;
  prf[1].enc = prf[0].enc + (prf[1].rad - prf[0].rad) * prf[1].rad * prf[0].rad * 0.5 * (prf[0].rho + prf[1].rho);
#if 1
  prf[0].psi = prf[0].enc / prf[0].rad;
  prf[1].psi = prf[1].enc / prf[1].rad;
#else
  prf[0].psi = prf[0].rad * prf[0].rad * prf[0].rho / 3.0;
  prf[1].psi = prf[1].rad * prf[1].rad * prf[1].rho / 3.0;
#endif
  Menc[0] += Menc[1] * 4.0;

  for(int ii = 2; ii < NRADBIN; ii++){
    const double mass = prf[ii].rad * prf[ii].rad * prf[ii].rad * prf[ii].rho;
    const int idx = ii & 1;

    prf[ii].enc = prf[idx].enc + (mass + Menc[idx]) * logrbin * M_LN10 / 3.0;
    prf[ii].psi = prf[ii].enc / prf[ii].rad;

    Menc[0] += mass * (double)(1 << (1 + (idx    )));
    Menc[1] += mass * (double)(1 << (1 + (idx ^ 1)));
  }/* for(int ii = 2; ii < NRADBIN; ii++){ */
#if 0
  if( cfg->kind == ISOTHERMAL_GAS ){
    fprintf(stdout, "internal: rho0 = %e, enc0 = %e, psi0 = %e, psi0_tot = %e\n", prf[0].rho, prf[0].enc, prf[0].psi, prf[0].psi_tot);
    fprintf(stdout, "internal: rho1 = %e, enc1 = %e, psi1 = %e, psi1_tot = %e\n", prf[1].rho, prf[1].enc, prf[1].psi, prf[1].psi_tot);
    fprintf(stdout, "internal: rho2 = %e, enc2 = %e, psi2 = %e, psi2_tot = %e\n", prf[2].rho, prf[2].enc, prf[2].psi, prf[2].psi_tot);
  }
#endif

  /** calculate potential (external part) */
  double Pext[2];
  Pext[0] = prf[NRADBIN - 1].rad * prf[NRADBIN - 1].rad * prf[NRADBIN - 1].rho;
  Pext[1] = prf[NRADBIN - 2].rad * prf[NRADBIN - 2].rad * prf[NRADBIN - 2].rho;
  double Pini[2];
  Pini[0] =           0.0;
  Pini[1] = Pini[0] + (prf[NRADBIN - 1].rad - prf[NRADBIN - 2].rad) * sqrt(prf[NRADBIN - 2].rad * prf[NRADBIN - 1].rad) * 0.5 * (prf[NRADBIN - 2].rho + prf[NRADBIN - 1].rho);
  Pext[0] += Pext[1] * 4.0;

  for(int ii = NRADBIN - 3; ii >= 0; ii--){
    const double psi = prf[ii].rad * prf[ii].rad * prf[ii].rho;
    const int idx = (int)((ii + 1) & 1);

    prf[ii].psi += Pini[idx] + (psi + Pext[idx]) * logrbin * M_LN10 / 3.0;

    Pext[0] += psi * (double)(1 << (1 + (idx    )));
    Pext[1] += psi * (double)(1 << (1 + (idx ^ 1)));
  }/* for(int ii = NRADBIN - 3; ii >= 0; ii--){ */
#if 0
  if( cfg->kind == ISOTHERMAL_GAS ){
    fprintf(stdout, "external: enc0 = %e, psi0 = %e\n", prf[0].enc, prf[0].psi);
    fprintf(stdout, "external: enc1 = %e, psi1 = %e\n", prf[1].enc, prf[1].psi);
    fprintf(stdout, "external: enc2 = %e, psi2 = %e\n", prf[2].enc, prf[2].psi);
  }
#endif

#if 1
  /** remove spike from the potential profile */
  double psi_in = prf[0].psi;
  for(int ii = 1; ii < NRADBIN - 1; ii++){
  if( prf[ii].psi > psi_in ){
  prf[ii].psi = fmin(psi_in, 0.5 * (psi_in + prf[ii + 1].psi));
}
  psi_in = prf[ii].psi;
}
#endif

#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE


  /** multiply overall factors */
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    prf[ii].enc *= 4.0 * M_PI;
    prf[ii].psi *= 4.0 * M_PI * CAST_R2D(newton);
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */


  /** set appropriate unit system */
  const double Mscale = Mtot / prf[NRADBIN - 1].enc;
#ifdef  ENABLE_GASEOUS_COMPONENT
  const bool initialize_rho0 = ((cfg->kind != KING) && (cfg->kind != ISOTHERMAL_GAS));
#else///ENABLE_GASEOUS_COMPONENT
  const bool initialize_rho0 =  (cfg->kind != KING);
#endif//ENABLE_GASEOUS_COMPONENT
  cfg->rho0 = initialize_rho0 ? (Mscale) : (cfg->rho0 * Mscale);
#if 0
  fprintf(stdout, "#Mtot = %e, Mcalc = %e, Mscale = %e, update = %d: rho0 = %e\n", Mtot, prf[NRADBIN - 1].enc, Mscale, !initialize_rho0, cfg->rho0);
  /* fprintf(stdout, "rho0 = %e, enc0 = %e, psi0 = %e\n", prf[0].rho, prf[0].enc, prf[0].psi); */
  /* fprintf(stdout, "rho1 = %e, enc1 = %e, psi1 = %e\n", prf[128].rho, prf[128].enc, prf[128].psi); */
#endif
#if 0
  if( isnan(Mscale) )
    for(int ii = 0; ii < NRADBIN; ii += 32)
      fprintf(stderr, "%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].enc, prf[ii].psi);
  fflush(NULL);
#endif
#if 0
  if( cfg->kind == ISOTHERMAL_GAS ){
    for(int ii = 0; ii < NRADBIN; ii++){
      fprintf(stdout, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].rho_tot, prf[ii].enc, prf[ii].enc_tot, prf[ii].psi, prf[ii].psi_tot);
      if( prf[ii].psi < 1.0e+4 )
	break;
    }
  }
#endif
#if 0
  if( (cfg->kind == ISOTHERMAL_GAS) && (Mscale < 10.0) ){
    for(int ii = 0; ii < NRADBIN; ii += (NRADBIN / 128))
      fprintf(stdout, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].rho_tot, prf[ii].enc, prf[ii].enc_tot, prf[ii].psi, prf[ii].psi_tot);
    __KILL__(stderr, "suspend for debugging\n");
  }
#endif
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    prf[ii].  rho     *= Mscale;
    prf[ii].  enc     *= Mscale;
    prf[ii].  psi     *= Mscale;
    prf[ii]. drho_dr  *= Mscale;
    prf[ii].d2rho_dr2 *= Mscale;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

#if 0
  if( cfg->kind == ISOTHERMAL_GAS ){
    fprintf(stdout, "#rescaled: rho0 = %e, enc0 = %e, psi0 = %e, psi0_tot = %e\n", prf[0].rho, prf[0].enc, prf[0].psi, prf[0].psi_tot);
    fprintf(stdout, "#rescaled: rho1 = %e, enc1 = %e, psi1 = %e, psi1_tot = %e\n", prf[1].rho, prf[1].enc, prf[1].psi, prf[1].psi_tot);
    fprintf(stdout, "#rescaled: rho2 = %e, enc2 = %e, psi2 = %e, psi2_tot = %e\n", prf[2].rho, prf[2].enc, prf[2].psi, prf[2].psi_tot);
  }
#endif

#if 0
  if( cfg->kind == ISOTHERMAL_GAS ){
    for(int ii = 0; ii < NRADBIN; ii += (NRADBIN / 128))
      fprintf(stdout, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].rho_tot, prf[ii].enc, prf[ii].enc_tot, prf[ii].psi, prf[ii].psi_tot);
    __KILL__(stderr, "suspend for debugging\n");
  }
#endif

#if 0
  for(int ii = iout - 3; ii >= 0; ii--)
    if( prf[ii].enc < Mtot ){
      iout = ii + 2;
      break;
    }
  rmax = prf[iout].rad;
  cfg->iout = iout;
  cfg->rmax = rmax;

  for(int ii = 0; ii < NRADBIN; ii++)
    fprintf(stderr, "%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].enc);
  exit(0);
#endif


#if 0
  for(int ii = 0; ii < NRADBIN; ii++)
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].enc, prf[ii].psi, prf[ii].drho_dr, prf[ii].d2rho_dr2);

  fprintf(stdout, "Mscale = %e\n", Mscale);
  fflush(NULL);

  exit(0);
#endif

  __NOTE__("%s\n", "end");
}


/**
 * @fn getAsymptoticSersicScale
 *
 * @brief Equation (25) in Ciotti & Bertin (1999), A&A, 352, 447--451
 *
 * @param (nn) Sersic index n
 * @return scale factor b_n
 */
/** based on : Eq.(25) */
static inline double getAsymptoticSersicScale(const double nn){  return (2.0 * nn + (-1.0 + 2.0 * (2.0 + 23.0 / (63.0 * nn)) / (135.0 * nn)) / 3.0);}


/**
 * @fn writeProfileCfgFormat
 *
 * @brief Print the expected format.
 *
 * @param (file) the specified file name
 * @param (cfg) physical properties of the component
 */
static inline void writeProfileCfgFormat(char *filename, const profile_cfg cfg)
{
  __NOTE__("%s\n", "start");


  fprintf(stderr, "ERROR: data written in \"%s\" does not match with format of specified model id \"%d\"\n", filename, cfg.kind);
  fprintf(stderr, "Expected format is below:\n");
  fprintf(stderr, "\tMtot<real>: total mass of the model in astrophysical units\n");
  if( cfg.kind != CENTRALBH ){
    fprintf(stderr, "\trs<real>: scale length of the model in astrophysical units\n");
#ifdef  USE_OSIPKOV_MERRITT_METHOD
    if( (cfg.kind != EXP_DISK) && (cfg.kind != SERSIC) && (cfg.kind != TBL_DISK) )
      fprintf(stderr, "\tra<real>: anisotropy radius of the model in astrophysical units (set 10000 rs as ra when a negative value is specified)\n");
#endif//USE_OSIPKOV_MERRITT_METHOD
  }/* if( cfg.kind != CENTRALBH ){ */

  /** some models require more information */
  if( cfg.kind == KING )
#ifdef  KING_CENTRAL_CUSP
    fprintf(stderr, "\tW0<real> dWdx_0<real>: dimensionless King parameter at the center\n");
#else///KING_CENTRAL_CUSP
    fprintf(stderr, "\tW0<real>: dimensionless King parameter at the center\n");
#endif//KING_CENTRAL_CUSP
  if( cfg.kind == APP_KING )
    fprintf(stderr, "\trt<real>: tidal radius of King model in empirical form in astrophysical units\n");
  if( cfg.kind == EINASTO )
    fprintf(stderr, "\talpha<real>: shape parameter of the Einasto model\n");
  if( cfg.kind == TWO_POWER ){
    fprintf(stderr, "\talpha<real>:    inner power-law slope of the two-power model\n");
    fprintf(stderr, "\t beta<real>: internal power-law slope of the two-power model\n");
    fprintf(stderr, "\tgamma<real>:    outer power-law slope of the two-power model\n");
  }/* if( cfg.kind == TWO_POWER ){ */
  if( cfg.kind == TRI_POWER ){
    fprintf(stderr, "\t   rout<real>: outer transition radius in astrophysical units\n");
    fprintf(stderr, "\t  alpha<real>:    innermost power-law slope of the three-power model\n");
    fprintf(stderr, "\t   beta<real>: transitional power-law slope of the three-power model\n");
    fprintf(stderr, "\t  gamma<real>: intermediate power-law slope of the three-power model\n");
    fprintf(stderr, "\t  delta<real>: transitional power-law slope of the three-power model\n");
    fprintf(stderr, "\tepsilon<real>:    outermost power-law slope of the three-power model\n");
  }/* if( cfg.kind == TRI_POWER ){ */
  if( cfg.kind == APP_EVANS ){
    fprintf(stderr, "\talpha<real>: inner power-law slope of the approximated Lowered Evans model\n");
    fprintf(stderr, "\t   rc<real>: second scale radius of the approximated Lowered Evans model\n");
    fprintf(stderr, "\t beta<real>: outer power-law slope of the approximated Lowered Evans model\n");
    fprintf(stderr, "\trt<real> wt<real>: exponential cutoff radius and width of the approximated Lowered Evans model, respectively\n");
  }/* if( cfg.kind == APP_EVANS ){ */
  if( cfg.kind == TABLE_RHO )
    fprintf(stderr, "\ttable<char *>: file name to be read density profile in table form\n");
  if( cfg.kind == TABLE_SIG )
    fprintf(stderr, "\ttable<char *>: file name to be read density profile in table form\n");
  if( cfg.kind == SPHSERSIC )
    fprintf(stderr, "\tn_sersic<real>: Sersic index\n");
  if( cfg.kind == SIGTWOPOW ){
    fprintf(stderr, "\talpha<real>:    inner power-law slope of the projected two-power model\n");
    fprintf(stderr, "\t beta<real>: internal power-law slope of the projected two-power model\n");
    fprintf(stderr, "\tgamma<real>:    outer power-law slope of the projected two-power model\n");
  }/* if( cfg.kind == SIGTWOPOW ){ */
  /* if( cfg.kind == ISOTHERMAL_GAS ){ */
  /*   fprintf(stderr, "\tTgas<real>: temperature of the gaseous component\n"); */
  /*   fprintf(stderr, "\t  mu<real>: mean molecular weight of the gaseous component\n"); */
  /* }/\* if( cfg.kind == ISOTHERMAL_GAS ){ *\/ */
  if( cfg.kind == TBL_DISK )
    fprintf(stderr, "\ttable<char *>: file name to be read column density profile in table form\n");
  if( (cfg.kind == EXP_DISK) || (cfg.kind == SERSIC) || (cfg.kind == TBL_DISK) ){
    if( cfg.kind == SERSIC )
      fprintf(stderr, "\tn_sersic<real>: Sersic index\n");
    fprintf(stderr, "\tRt<real> Rt_width<real>: cutoff radius and width of the disk mid-plane density in horizontal direction in astrophysical units, respectively\n");
    fprintf(stderr, "\tzd<real>: scale height of isothermal disk in the vertical direction in astrophysical units\n");
    fprintf(stderr, "\tsigmaR0<real> param<real>: velocity dispersion in radial direction at the center in astrophysical units (negative value indicates sigmaR0 = sigmaz0), perpendicular velocity dispersion over circular velocity\n");
#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
    fprintf(stderr, "\t\tparam is used as frac (scaling factor f in Miki & Umemura (2018), MNRAS, 475, 2269).\n");
    fprintf(stderr, "\t\tif the inputted sigmaR0 is negative, then the default value (= sigma_z(R = 0)) is substituted\n");
#else///ENFORCE_EPICYCLIC_APPROXIMATION
    fprintf(stderr, "\t\tparam is used as Toomre's Q-value at the scale length.\n");
    fprintf(stderr, "\t\tif the inputted sigmaR0 is negative, then sigmaR0 is automatically determined to obtain the inputted Toomre's Q-value.\n");
    fprintf(stderr, "\t\tif the inputted sigmaR0 and param is both negative, then the default value (= sigma_z(R = 0)) is substituted.\n");
#endif//ENFORCE_EPICYCLIC_APPROXIMATION
#ifndef USE_ZANG_HOHL_1978_EQ5
    fprintf(stderr, "\t\tif the inputted retrogradeFrac<real> is not zero ([0., 1.]), then rotation axis of particles with the given fraction is anti-parallel with normal component.\n");
#else///USE_ZANG_HOHL_1978_EQ5
    fprintf(stderr, "\t\tif the inputted retrogradeParam<real> is not zero ([0., 1.]), then particles satisfying equation (5) of Zang & Hohl (1978) are picked upped as retrograding component.\n");
#endif//USE_ZANG_HOHL_1978_EQ5
  }/* if( (cfg.kind == EXP_DISK) || (cfg.kind == SERSIC) || (cfg.kind == TBL_DISK) ){ */

  /** information about density cutoff */
  if( cfg.kind != CENTRALBH ){
    fprintf(stderr, "\tcutoff<bool>: set explicit density cutoff (1) or no (0)\n");
    if( cfg.cutoff )
      fprintf(stderr, "\trc<real> rc_width<real>: cutoff radius and width of the density cutoff in astrophysical units, respectively\n");
  }/* if( cfg.kind != CENTRALBH ) */

  __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);

  __NOTE__("%s\n", "end");
}

/**
 * @fn readProfileCfg
 *
 * @brief Read physical properties of the spherical component(s) with memory allocation.
 *
 * @param (fcfg) file name of the configuration
 * @return (unit) unit system
 * @return (kind) number of components
 * @return (cfg) physical properties of the component(s)
 *
 * @sa setPhysicalConstantsAndUnitSystem
 * @sa writeProfileCfgFormat
 */
void readProfileCfg(char *fcfg, int *unit, int *kind, profile_cfg **cfg)
{
  __NOTE__("%s\n", "start");


  /** read global settings */
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s", CFGFOLDER, fcfg);

  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  int checker = 1;

  /** read the specified unit system and set it */
  checker &= (1 == fscanf(fp, "%d", unit));
  setPhysicalConstantsAndUnitSystem(*unit, 1);
  /* setPhysicalConstantsAndUnitSystem(3, 1); */
  /* __KILL__(stderr, "suspend for confirmation\n"); */
  /** read # of components */
  checker &= (1 == fscanf(fp, "%d", kind));

  *cfg = (profile_cfg *)malloc(sizeof(profile_cfg) * (*kind));  if( *cfg == NULL ){    __KILL__(stderr, "ERROR: failure to allocate cfg\n");  }
  for(int ii = 0; ii < *kind; ii++)
    checker &= (4 == fscanf(fp, "%d\t%s\t%d\t%zu", &(*cfg)[ii].kind, (*cfg)[ii].file, &(*cfg)[ii].forceNum, &(*cfg)[ii].num));

  fclose(fp);
  if( !checker ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }


  /** read individual settings */
  for(int ii = 0; ii < *kind; ii++){
    sprintf(filename, "%s/%s", CFGFOLDER, (*cfg)[ii].file);
    fp = fopen(filename, "r");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    checker = 1;

    checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].Mtot));    (*cfg)[ii].Mtot *= mass_astro2com;
    if( (*cfg)[ii].kind != CENTRALBH ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].rs));      (*cfg)[ii].rs *= length_astro2com;
#ifdef  USE_OSIPKOV_MERRITT_METHOD
      if( ((*cfg)[ii].kind != EXP_DISK) && ((*cfg)[ii].kind != SERSIC) && ((*cfg)[ii].kind != TBL_DISK) ){
	(*cfg)[ii].ra = -1.0;
	checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].ra));
	(*cfg)[ii].ra *= length_astro2com;
	if( (*cfg)[ii].ra <= 0.0 )
	  (*cfg)[ii].ra  = 10000.0 * (*cfg)[ii].rs;
	(*cfg)[ii].ra2inv = 1.0 / (DBL_MIN + (*cfg)[ii].ra * (*cfg)[ii].ra);
      }/* if( ((*cfg)[ii].kind != EXP_DISK) && ((*cfg)[ii].kind != SERSIC) && ((*cfg)[ii].kind != TBL_DISK) ){ */
#endif//USE_OSIPKOV_MERRITT_METHOD
    }/* if( (*cfg)[ii].kind != CENTRALBH ){ */
    else
      if( ((*cfg)[ii].forceNum != 1) || ((*cfg)[ii].num != 1) ){
	__KILL__(stderr, "ERROR: number of BH particle must be specified to be unity\n");
      }/* if( ((*cfg)[ii].forceNum != 1) || ((*cfg)[ii].num != 1) ){ */

    /** parameter for King profile */
    if( (*cfg)[ii].kind ==     KING )
#ifdef  KING_CENTRAL_CUSP
      checker &= (2 == fscanf(fp, "%le %le", &(*cfg)[ii].king_W0, &(*cfg)[ii].king_dWdx_0));
#else///KING_CENTRAL_CUSP
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].king_W0));
#endif//KING_CENTRAL_CUSP
    if( (*cfg)[ii].kind == APP_KING ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].king_rt));      (*cfg)[ii].king_rt *= length_astro2com;
    }/* if( (*cfg)[ii].kind == APP_KING ){ */

    /** parameter for Einasto profile */
    if( (*cfg)[ii].kind == EINASTO )
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].einasto_alpha));

    /** parameters for two-power model */
    if( (*cfg)[ii].kind == TWO_POWER ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_alpha));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_beta));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_gamma));
    }/* if( (*cfg)[ii].kind == TWO_POWER ){ */

    /** parameters for three-power model */
    if( (*cfg)[ii].kind == TRI_POWER ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].tripower_rout));      (*cfg)[ii].tripower_rout *= length_astro2com;
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_alpha));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_beta));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_gamma));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].tripower_delta));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].tripower_epsilon));
    }/* if( (*cfg)[ii].kind == TRI_POWER ){ */

    /** parameters for approximated lowered Evans model */
    if( (*cfg)[ii].kind == APP_EVANS ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].alevans_alpha));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].alevans_rc));      (*cfg)[ii].alevans_rc *= length_astro2com;
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].alevans_beta));
      checker &= (2 == fscanf(fp, "%le %le", &(*cfg)[ii].alevans_rt, &(*cfg)[ii].alevans_wt));
      (*cfg)[ii].alevans_rt *= length_astro2com;
      (*cfg)[ii].alevans_wt *= length_astro2com;
    }/* if( (*cfg)[ii].kind == APP_EVANS ){ */

    /** parameter for density profile in table form */
    if( (*cfg)[ii].kind == TABLE_RHO )
      checker &= (1 == fscanf(fp, "%s", (*cfg)[ii].table));

    /** parameter for column density profile in table form */
    if( (*cfg)[ii].kind == TABLE_SIG )
      checker &= (1 == fscanf(fp, "%s", (*cfg)[ii].table));

    /** parameter for spherical Serisic profile */
    if( (*cfg)[ii].kind == SPHSERSIC ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].n_sersic));
      (*cfg)[ii].b_sersic = getAsymptoticSersicScale((*cfg)[ii].n_sersic);
    }/* if( (*cfg)[ii].kind == SPHSERSIC ){ */

    /** parameters for projected two-power model */
    if( (*cfg)[ii].kind == SIGTWOPOW ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_alpha));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_beta));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_gamma));
    }/* if( (*cfg)[ii].kind == SIGTWOPOW ){ */

    /** parameter for column density profile in table form */
    if( (*cfg)[ii].kind == TBL_DISK )
      checker &= (1 == fscanf(fp, "%s", (*cfg)[ii].table));

    /** parameters for disk component */
    if( ((*cfg)[ii].kind == EXP_DISK) || ((*cfg)[ii].kind == SERSIC) || ((*cfg)[ii].kind == TBL_DISK) ){
      if( (*cfg)[ii].kind == SERSIC ){
	checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].n_sersic));
	(*cfg)[ii].b_sersic = getAsymptoticSersicScale((*cfg)[ii].n_sersic);
      }/* if( (*cfg)[ii].kind == SERSIC ){ */
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].zd));      (*cfg)[ii].zd   *= length_astro2com;
      checker &= (2 == fscanf(fp, "%le %le", &(*cfg)[ii].vdispR0
#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
			      , &(*cfg)[ii].vdisp_frac)
#else///ENFORCE_EPICYCLIC_APPROXIMATION
			      , &(*cfg)[ii].toomre)
#endif//ENFORCE_EPICYCLIC_APPROXIMATION
		  );
      (*cfg)[ii].vdispR0 *= velocity_astro2com;
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].retrogradeFrac));
    }/* if( ((*cfg)[ii].kind == EXP_DISK) || ((*cfg)[ii].kind == SERSIC) || ((*cfg)[ii].kind == TBL_DISK) ){ */

  /** parameters for density cutoff */
  if( (*cfg)[ii].kind != CENTRALBH ){
      int input;
      checker &= (1 == fscanf(fp, "%d", &input));
      (*cfg)[ii].cutoff = (bool)input;
      if( (*cfg)[ii].cutoff ){
	checker &= (2 == fscanf(fp, "%le %le", &(*cfg)[ii].rc, &(*cfg)[ii].rc_width));
	(*cfg)[ii].rc       *= length_astro2com;
	(*cfg)[ii].rc_width *= length_astro2com;
      }/* if( (*cfg)[ii].cutoff ){ */
    }/* if( (*cfg)[ii].kind != CENTRALBH ){ */

    fclose(fp);
    if( !checker )
      writeProfileCfgFormat(filename, (*cfg)[ii]);
  }/* for(int ii = 0; ii < *kind; ii++){ */


  __NOTE__("%s\n", "end");
}

#ifdef  ENABLE_GASEOUS_COMPONENT

/**
 * @fn gas_writeProfileCfgFormat
 *
 * @brief Print the expected format.
 *
 * @param (file) the specified file name
 * @param (cfg) physical properties of the component
 */
static inline void gas_writeProfileCfgFormat(char *filename, const profile_cfg cfg)
{
  __NOTE__("%s\n", "start");


  fprintf(stderr, "ERROR: data written in \"%s\" does not match with format of specified model id \"%d\"\n", filename, cfg.kind);
  fprintf(stderr, "Expected format is below:\n");
  fprintf(stderr, "\tMtot<real>: total mass of the model in astrophysical units\n");

  /** some models require more information */
  if( cfg.kind == ISOTHERMAL_GAS ){
    fprintf(stderr, "\tTgas<real>: temperature of the gaseous component\n");
    fprintf(stderr, "\t  mu<real>: mean molecular weight of the gaseous component\n");

    fprintf(stderr, "\trc<real> rc_width<real>: cutoff radius and width of the density cutoff in astrophysical units, respectively\n");
  }/* if( cfg.kind == ISOTHERMAL_GAS ){ */

  __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);

  __NOTE__("%s\n", "end");
}

/**
 * @fn gas_readProfileCfg
 *
 * @brief Read physical properties of the spherical gaseous component(s) with memory allocation.
 *
 * @param (fcfg) file name of the configuration
 * @param (unit) unit system
 * @return (kind) number of components
 * @return (cfg) physical properties of the component(s)
 *
 * @sa setPhysicalConstantsAndUnitSystem
 * @sa writeProfileCfgFormat
 */
void gas_readProfileCfg(char *fcfg, const int unit, int *kind, profile_cfg **cfg)
{
  __NOTE__("%s\n", "start");


  /** read global settings */
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s", CFGFOLDER, fcfg);

  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  int checker = 1;

  /** read the specified unit system and set it */
  int unit_read;
  checker &= (1 == fscanf(fp, "%d", &unit_read));
  if( unit_read != unit ){
    __KILL__(stderr, "ERROR: unit system in \"%s\" is %d; however, it is inconsistent with %d for N-body particles\n", filename, unit_read, unit);
  }
  /** read # of components */
  checker &= (1 == fscanf(fp, "%d", kind));

  *cfg = (profile_cfg *)malloc(sizeof(profile_cfg) * (*kind));  if( *cfg == NULL ){    __KILL__(stderr, "ERROR: failure to allocate cfg\n");  }
  for(int ii = 0; ii < *kind; ii++)
    checker &= (4 == fscanf(fp, "%d\t%s\t%d\t%zu", &(*cfg)[ii].kind, (*cfg)[ii].file, &(*cfg)[ii].forceNum, &(*cfg)[ii].num));

  fclose(fp);
  if( !checker ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }


  /** read individual settings */
  for(int ii = 0; ii < *kind; ii++){
    sprintf(filename, "%s/%s", CFGFOLDER, (*cfg)[ii].file);
    fp = fopen(filename, "r");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    checker = 1;

    checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].Mtot));    (*cfg)[ii].Mtot *= mass_astro2com;

    /** parameters for isothermal gaseous component */
    if( (*cfg)[ii].kind == ISOTHERMAL_GAS ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].Tgas));      (*cfg)[ii].Tgas *= temperature_astro2com;
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].mu));

      (*cfg)[ii].cutoff = true;
      checker &= (2 == fscanf(fp, "%le %le", &(*cfg)[ii].rc, &(*cfg)[ii].rc_width));
      (*cfg)[ii].rc       *= length_astro2com;
      (*cfg)[ii].rc_width *= length_astro2com;
    }/* if( (*cfg)[ii].kind == ISOTHERMAL_GAS ){ */

    fclose(fp);
    if( !checker )
      gas_writeProfileCfgFormat(filename, (*cfg)[ii]);
  }/* for(int ii = 0; ii < *kind; ii++){ */


  __NOTE__("%s\n", "end");
}

#endif//ENABLE_GASEOUS_COMPONENT


#ifdef  MAKE_COLUMN_DENSITY_PROFILE
#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
/**
 * @fn get_DEformula
 *
 * @brief Calculate integrand in double exponential formula.
 *
 * @param (tt) value of integration variable
 * @param (rr) 
 * @return (ret) integrand for 
 * @param (xx) position of data points (psi)
 * @param (yy) value of data points (d2rho_dpsi2)
 * @param (y2) coefficients in cubic spline interpolation
 */
static inline void get_DEformula
(const double tt, const double R2, const double zmax, const int skind, double ret[restrict], double * restrict xx, double * restrict yy, double * restrict y2
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 , double v2f[restrict], double * restrict f2, double * restrict d2
 , double v4f[restrict], double * restrict f4, double * restrict d4
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double cosh_t = cosh(sinh_t);
  const double inv_cosh_t = 1.0 / cosh_t;

  const double common = cosh(tt) * inv_cosh_t * inv_cosh_t;

  const double zz = 0.5 * zmax * exp(sinh_t) * inv_cosh_t;
  const double rr = sqrt(R2 + zz * zz);

  for(int kk = 0; kk < skind; kk++){
    ret[kk] += common * getCubicSpline1D(rr, NRADBIN, xx, &yy[kk * NRADBIN], &y2[kk * NRADBIN]);
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    v2f[kk] += common * getCubicSpline1D(rr, NRADBIN, xx, &f2[kk * NRADBIN], &d2[kk * NRADBIN]);
    v4f[kk] += common * getCubicSpline1D(rr, NRADBIN, xx, &f4[kk * NRADBIN], &d4[kk * NRADBIN]);
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }/* for(int ii = 0; ii < skind; ii++){ */

}


static inline void update_trapezoidal
(const double hh, const double tmin, const double tmax, const int skind, double sum[restrict], const double R2, const double zmax, double * restrict xx, double * restrict yy, double * restrict y2, double sub[restrict]
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 , double v2f[restrict], double * restrict f2, double * restrict d2, double sub2[restrict]
 , double v4f[restrict], double * restrict f4, double * restrict d4, double sub4[restrict]
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
)
{
  /** initialization */
  for(int kk = 0; kk < skind; kk++){
    sub[kk] = 0.0;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    sub2[kk] = sub4[kk] = 0.0;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }/* for(int kk = 0; kk < skind; kk++){ */
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    get_DEformula(tt, R2, zmax, skind, sub, xx, yy, y2
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		  , sub2, f2, d2, sub4, f4, d4
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		  );
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  for(int kk = 0; kk < skind; kk++){
    sum[kk] = 0.5 * sum[kk] + hh * sub[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    v2f[kk] = 0.5 * v2f[kk] + hh * sub2[kk];
    v4f[kk] = 0.5 * v4f[kk] + hh * sub4[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }/* for(int kk = 0; kk < skind; kk++){ */
}


static inline void set_domain_boundary
(const double hh, double * restrict tmin, double * restrict tmax, const double R2, const double zmax, const int skind, double sum[restrict], double * restrict xx, double * restrict yy, double * restrict y2, double fp[restrict], double ft[restrict], double f0[restrict]
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 , double v2f[restrict], double * restrict f2, double * restrict d2, double v2fp[restrict], double v2ft[restrict], double v2f0[restrict]
 , double v4f[restrict], double * restrict f4, double * restrict d4, double v4fp[restrict], double v4ft[restrict], double v4f0[restrict]
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 )
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  for(int kk = 0; kk < skind; kk++){
    fp[kk] = 0.0;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    v2fp[kk] = v4fp[kk] = 0.0;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }/* for(int kk = 0; kk < skind; kk++){ */
  get_DEformula(tt, R2, zmax, skind, fp, xx, yy, y2
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		, v2fp, f2, d2, v4fp, f4, d4
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		);
  for(int kk = 0; kk < skind; kk++){
    const double tmp = fp[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    const double tmp2 = v2fp[kk];
    const double tmp4 = v4fp[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    f0[kk] = tmp;    sum[kk] = hh * tmp;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    v2f0[kk] = tmp2;    v2f[kk] = hh * tmp2;
    v4f0[kk] = tmp4;    v4f[kk] = hh * tmp4;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }/* for(int ii = 0; ii < skind; ii++){ */


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    for(int kk = 0; kk < skind; kk++){
      ft[kk] = fp[kk];      fp[kk] = 0.0;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      v2ft[kk] = v2fp[kk];      v2fp[kk] = 0.0;
      v4ft[kk] = v4fp[kk];      v4fp[kk] = 0.0;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    }/* for(int kk = 0; kk < skind; kk++){ */

    tt += hh;    boundary = tt;
    get_DEformula(tt, R2, zmax, skind, fp, xx, yy, y2
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		  , v2fp, f2, d2, v4fp, f4, d4
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		  );

    damp = -1.0;
    for(int kk = 0; kk < skind; kk++){
      sum[kk] += hh * fp[kk];      damp = fmax(damp, fabs(ft[kk]) + fabs(fp[kk]));
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      v2f[kk] += hh * v2fp[kk];      damp = fmax(damp, fabs(v2ft[kk]) + fabs(v2fp[kk]));
      v4f[kk] += hh * v4fp[kk];      damp = fmax(damp, fabs(v4ft[kk]) + fabs(v4fp[kk]));
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    }/* for(int ii = 0; ii < skind; ii++){ */

  }/* while( (damp > converge) && (boundary < maximum) ){ */
  *tmax = boundary;


  /** determine lower boundary */
  const double RR = sqrt(R2);
  for(int kk = 0; kk < skind; kk++){
    fp[kk] = f0[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    v2fp[kk] = v2f0[kk];
    v4fp[kk] = v4f0[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }/* for(int kk = 0; kk < skind; kk++){ */
  tt = 0.0;
  boundary = 0.0;
  damp = 1.0;
  while( (damp > converge) && (boundary > -maximum) ){
    for(int kk = 0; kk < skind; kk++){
      ft[kk] = fp[kk];      fp[kk] = 0.0;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      v2ft[kk] = v2fp[kk];      v2fp[kk] = 0.0;
      v4ft[kk] = v4fp[kk];      v4fp[kk] = 0.0;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    }/* for(int kk = 0; kk < skind; kk++){ */

    tt -= hh;    boundary = tt;
    get_DEformula(tt, R2, zmax, skind, fp, xx, yy, y2
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		  , v2fp, f2, d2, v4fp, f4, d4
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		  );

    damp = -1.0;
    for(int kk = 0; kk < skind; kk++){
      sum[kk] += hh * fp[kk];      damp = fmax(damp, fabs(ft[kk]) + fabs(fp[kk]));
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      v2f[kk] += hh * v2fp[kk];      damp = fmax(damp, fabs(v2ft[kk]) + fabs(v2fp[kk]));
      v4f[kk] += hh * v4fp[kk];      damp = fmax(damp, fabs(v4ft[kk]) + fabs(v4fp[kk]));
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    }/* for(int ii = 0; ii < skind; ii++){ */

    const double yy = M_PI_2 * sinh(tt);
    if( 0.5 * zmax * exp(yy) / cosh(yy) > 1.01 * RR )
      damp = 1.0;
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;
}


static inline void integrate_DEformula
(const double R2, const double zmax, const int skind, double sum[restrict], double * restrict xx, double * restrict yy, double * restrict y2, double sub[restrict], double ft[restrict], double f0[restrict]
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 , double v2f[restrict], double * restrict f2, double * restrict d2, double v2fs[restrict], double v2ft[restrict], double v2f0[restrict]
 , double v4f[restrict], double * restrict f4, double * restrict d4, double v4fs[restrict], double v4ft[restrict], double v4f0[restrict]
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
)
{
  const double criteria_abs = 1.0e-12;
  /* const double criteria_rel = 1.0e-8; */
  const double criteria_rel = 1.0e-5;
  /* const double criteria_rel = 1.0e-4; */

  double hh = 1.0;
  double tmin, tmax;
  set_domain_boundary(hh, &tmin, &tmax, R2, zmax, skind, sum, xx, yy, y2, sub, ft, f0
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		      , v2f, f2, d2, v2fs, v2ft, v2f0, v4f, f4, d4, v4fs, v4ft, v4f0
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		      );


  while( true ){
    for(int kk = 0; kk < skind; kk++){
      f0[kk] = sum[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      v2f0[kk] = v2f[kk];
      v4f0[kk] = v4f[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    }/* for(int kk = 0; kk < skind; kk++){ */

    hh *= 0.5;
    update_trapezoidal(hh, tmin, tmax, skind, sum, R2, zmax, xx, yy, y2, sub
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		       , v2f, f2, d2, v2fs, v4f, f4, d4, v4fs
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
		       );

    bool converge = true;
    for(int kk = 0; kk < skind; kk++){
      if( converge ){
	if( fabs(sum[kk]) > DBL_EPSILON ){
	  if( fabs(1.0 - f0[kk] / sum[kk]) > criteria_rel )
	    converge = false;
	}
	else
	  if( fabs(sum[kk] - f0[kk]) > criteria_abs )
	    converge = false;
      }

#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      if( converge ){
	if( fabs(v2f[kk]) > DBL_EPSILON ){
	  if( fabs(1.0 - v2f0[kk] / v2f[kk]) > criteria_rel )
	    converge = false;
	}
	else
	  if( fabs(v2f[kk] - v2f0[kk]) > criteria_abs )
	    converge = false;
      }

      if( converge ){
	if( fabs(v4f[kk]) > DBL_EPSILON ){
	  if( fabs(1.0 - v4f0[kk] / v4f[kk]) > criteria_rel )
	    converge = false;
	}
	else
	  if( fabs(v4f[kk] - v4f0[kk]) > criteria_abs )
	    converge = false;
      }
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    }/* for(int kk = 0; kk < skind; kk++){ */

    if( converge )
      break;
  }/* while( true ){ */

#if 0
  if( fpclassify(sum[0]) != FP_NORMAL ){
    fprintf(stderr, "tmin = %e, tmax = %e, hh = %e, R2 = %e\n", tmin, tmax, hh, R2);
    fflush(NULL);
    exit(0);
  }
#endif

}


/**
 * @fn calcColumnDensityProfile
 *
 * @brief Calculate column density profile.
 *
 * @param (skind) number of spherical symmetric components
 * @return (prf) radial profile of the components
 * @param (logrmax) rmax in logarithmic space
 * @param (cfg) physical properties of the component(s)
 */
void calcColumnDensityProfile(const int skind, profile **prf,
#ifndef ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
			      const double logrmax,
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
			      profile_cfg *cfg)
{
  __NOTE__("%s\n", "start");


  /** memory allocation for cubic spline interpolation */
  double *xx;  xx = (double *)malloc(sizeof(double)         * NRADBIN);  if( xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xx\n");  }
  double *bp;  bp = (double *)malloc(sizeof(double) * skind * NRADBIN);  if( bp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate bp\n");  }
  double *yy;  yy = (double *)malloc(sizeof(double) * skind * NRADBIN);  if( yy == NULL ){    __KILL__(stderr, "ERROR: failure to allocate yy\n");  }
  double *y2;  y2 = (double *)malloc(sizeof(double) * skind * NRADBIN);  if( y2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate y2\n");  }
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  double *f2;  f2 = (double *)malloc(sizeof(double) * skind * NRADBIN);  if( f2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate f2\n");  }
  double *d2;  d2 = (double *)malloc(sizeof(double) * skind * NRADBIN);  if( d2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate d2\n");  }
  double *f4;  f4 = (double *)malloc(sizeof(double) * skind * NRADBIN);  if( f4 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate f4\n");  }
  double *d4;  d4 = (double *)malloc(sizeof(double) * skind * NRADBIN);  if( d4 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate d4\n");  }
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)

#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++)
    xx[ii] = prf[0][ii].rad;

  for(int kk = 0; kk < skind; kk++)
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      yy[ii + kk * NRADBIN] = prf[kk][ii].rho;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      f2[ii + kk * NRADBIN] = prf[kk][ii].v2f;
      f4[ii + kk * NRADBIN] = prf[kk][ii].v4f;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  /** execute cubic spline interpolation */
#pragma omp parallel for
  for(int kk = 0; kk < skind; kk++){
    genCubicSpline1D(NRADBIN, xx, &yy[kk * NRADBIN], bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, &y2[kk * NRADBIN]);
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    genCubicSpline1D(NRADBIN, xx, &f2[kk * NRADBIN], bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, &d2[kk * NRADBIN]);
    genCubicSpline1D(NRADBIN, xx, &f4[kk * NRADBIN], bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, &d4[kk * NRADBIN]);
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }/* for(int kk = 0; kk < skind; kk++){ */


  int iout = 0;
  for(int ii = 0; ii < skind; ii++)
    if( cfg[ii].kind != CENTRALBH )
      iout = (iout < cfg[ii].iout) ? cfg[ii].iout : iout;
  const double rmax2 = prf[0][iout].rad * prf[0][iout].rad;

#pragma omp parallel
  {
    double Sig[NKIND_MAX], Sig0[NKIND_MAX], Sig1[NKIND_MAX], Sig2[NKIND_MAX];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    double v2f[NKIND_MAX], v2f0[NKIND_MAX], v2f1[NKIND_MAX], v2f2[NKIND_MAX];
    double v4f[NKIND_MAX], v4f0[NKIND_MAX], v4f1[NKIND_MAX], v4f2[NKIND_MAX];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)

#pragma omp for schedule(dynamic, 8) nowait
    for(int ii = 0; ii < iout + 1; ii += SKIP_INTERVAL_FOR_COLUMN_DENSITY){
      const double R2 = prf[0][ii].rad * prf[0][ii].rad;
      const double zmax = sqrt(rmax2 - R2);

      /* call DE formula */
      integrate_DEformula(R2, zmax, skind, Sig, xx, yy, y2, Sig0, Sig1, Sig2
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			  , v2f, f2, d2, v2f0, v2f1, v2f2
			  , v4f, f4, d4, v4f0, v4f1, v4f2
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			  );

      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].Sigma = M_PI_2 * zmax * Sig[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	prf[kk][ii].slos  = sqrt(v4f[kk] / (DBL_MIN + 3.0 * v2f[kk]));
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      }/* for(int kk = 0; kk < skind; kk++){ */

    }/* for(int ii = 0; ii < iout + 1; ii += SKIP_INTERVAL_FOR_COLUMN_DENSITY){ */

    for(int kk = 0; kk < skind; kk++)
#pragma omp for nowait
      for(int ii = iout + 1; ii < NRADBIN; ii++)
	prf[kk][ii].Sigma = 0.0;

  }



#if 0
  for(int ii = 0; ii < NRADBIN; ii++)
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", prf[0][ii].rad, prf[0][ii].rho, prf[0][ii].enc, prf[0][ii].psi, prf[0][ii].Sigma);

  fflush(NULL);

  exit(0);
#endif


#   if  SKIP_INTERVAL_FOR_COLUMN_DENSITY != 1
#pragma omp parallel
  {
    for(int kk = 0; kk < skind; kk++){
      const int iout = cfg[kk].iout;

#pragma omp for schedule(dynamic, 8) nowait
      for(int ii = 0; ii < iout + 1; ii += SKIP_INTERVAL_FOR_COLUMN_DENSITY){
	const double S0 = prf[kk][ii                                   ].Sigma;
	const double S1 = prf[kk][ii + SKIP_INTERVAL_FOR_COLUMN_DENSITY].Sigma;
	const double Slope = (S1 - S0) / (double)SKIP_INTERVAL_FOR_COLUMN_DENSITY;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	const double s0 = prf[kk][ii                                   ].slos;
	const double s1 = prf[kk][ii + SKIP_INTERVAL_FOR_COLUMN_DENSITY].slos;
	const double slope = (s1 - s0) / (double)SKIP_INTERVAL_FOR_COLUMN_DENSITY;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)

	for(int jj = 1; jj < SKIP_INTERVAL_FOR_COLUMN_DENSITY; jj++){
	  prf[kk][ii + jj].Sigma = S0 + Slope * (double)jj;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	  prf[kk][ii + jj].slos  = s0 + slope * (double)jj;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	}/* for(int jj = 1; jj < SKIP_INTERVAL_FOR_COLUMN_DENSITY; jj++){ */
      }/* for(int ii = 0; ii < iout + 1; ii += SKIP_INTERVAL_FOR_COLUMN_DENSITY){ */

    }/* for(int kk = 0; kk < skind; kk++){ */
  }
#endif//SKIP_INTERVAL_FOR_COLUMN_DENSITY != 1


  free(xx);  free(bp);
  free(yy);  free(y2);
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  free(f2);  free(d2);
  free(f4);  free(d4);
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)


  __NOTE__("%s\n", "end");
}







#else///ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE

/**
 * @fn findIdx
 *
 * @brief Find a data element in the given array corresponding to the given value.
 *
 * @param (rad) radius
 * @param (prf) radial profile of the component
 * @return (ll) the corresponding lower index
 * @return (rr) the corresponding upper index
 */
static inline void findIdx(const double rad, profile * restrict prf, int * restrict ll, int * restrict rr)
{
  bool bisection = true;
  *ll =           0;
  *rr = NRADBIN - 1;

  if( bisection == true )    if( fabs(prf[*ll].rad - rad) / rad < DBL_EPSILON ){      bisection = false;      *rr = (*ll) + 1;    }
  if( bisection == true )    if( fabs(prf[*rr].rad - rad) / rad < DBL_EPSILON ){      bisection = false;      *ll = (*rr) - 1;    }

  while( bisection ){
    const uint cc = ((uint)(*ll) + (uint)(*rr)) >> 1;

    if( (prf[cc].rad - rad) * (prf[*ll].rad - rad) <= 0.0 )      *rr = (int)cc;
    else                                                         *ll = (int)cc;

    if( (1 + (*ll)) == (*rr) )      break;

  }/* while( bisection ){ */
}


/**
 * @fn getInterpolatedDensity
 *
 * @brief Get linear interpolated density.
 *
 * @param (RR) position in R-direction
 * @param (zz) position in z-direction
 * @param (skind) number of spherical symmetric components
 * @param (prf) radial profile of the components
 * @return (Sig) the corresponding density of each component
 * @return (v2f) the corresponding integral of each component
 * @return (v4f) the corresponding integral of each component
 */
static inline void getInterpolatedDensity
(const double RR, const double zz, const int skind, profile * restrict * prf, double * restrict val
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 , double * restrict v2f, double * restrict v4f
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 )
{
  const double rad = sqrt(RR * RR + zz * zz);

  int ll, rr;
  findIdx(rad, prf[0], &ll, &rr);

  /** based on linear interpolation */
  const double alpha = (rad - prf[0][ll].rad) / (prf[0][rr].rad - prf[0][ll].rad);

  for(int kk = 0; kk < skind; kk++){
    val[kk] = (1.0 - alpha) * prf[kk][ll].rho + alpha * prf[kk][rr].rho;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    v2f[kk] = (1.0 - alpha) * prf[kk][ll].v2f + alpha * prf[kk][rr].v2f;
    v4f[kk] = (1.0 - alpha) * prf[kk][ll].v4f + alpha * prf[kk][rr].v4f;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }/* for(int kk = 0; kk < skind; kk++){ */
}

/**
 * @fn gaussQuadVertical
 *
 * @brief Gaussian quadrature for calculating column density.
 *
 * @param (RR) position in R-direction
 * @param (num) number of data points
 * @param (zmin) minimum of z
 * @param (zmax) maximum of z
 * @param (skind) number of spherical symmetric components
 * @param (prf) radial profile of the components
 * @return (sum) column density of each component
 * @param (fm) temporary array
 * @param (fp) temporary array
 * @return (v2f) integral of v2f(r) of each component
 * @param (v2fm) temporary array
 * @param (v2fp) temporary array
 * @return (v4f) integral of v4f(r) of each component
 * @param (v4fm) temporary array
 * @param (v4fp) temporary array
 *
 * @sa getInterpolatedDensity
 */
void gaussQuadVertical
(const double RR, const int num, const double zmin, const double zmax, const int skind, profile * restrict * prf, double * restrict sum, double * restrict fm, double * restrict fp
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 , double * restrict v2f, double * restrict v2fm, double * restrict v2fp
 , double * restrict v4f, double * restrict v4fm, double * restrict v4fp
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 );
void gaussQuadVertical
(const double RR, const int num, const double zmin, const double zmax, const int skind, profile * restrict * prf, double * restrict sum, double * restrict fm, double * restrict fp
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 , double * restrict v2f, double * restrict v2fm, double * restrict v2fp
 , double * restrict v4f, double * restrict v4fm, double * restrict v4fp
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
 )
{
  const double mns = 0.5 * (zmax - zmin);
  const double pls = 0.5 * (zmax + zmin);

  for(int kk = 0; kk < skind; kk++){
    sum[kk] = 0.0;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    v2f[kk] = v4f[kk] = 0.0;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }/* for(int kk = 0; kk < skind; kk++){ */

  if( num & 1 ){
    const double weight =             gsl_gaussQD_weight[(num >> 1)];
    const double  value = pls + mns * gsl_gaussQD_pos   [(num >> 1)];
    getInterpolatedDensity(RR, value, skind, prf, fm
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			   , v2fm, v4fm
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			   );
    for(int kk = 0; kk < skind; kk++){
      sum[kk] = weight * fm[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      v2f[kk] = weight * v2fm[kk];
      v4f[kk] = weight * v4fm[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    }/* for(int kk = 0; kk < skind; kk++){ */
  }/* if( num & 1 ){ */

  for(int ii = (num >> 1) - 1; ii >= 0; ii--){
    const double weight = gsl_gaussQD_weight[ii];
    const double zp = pls + mns * gsl_gaussQD_pos[ii];
    const double zm = pls - mns * gsl_gaussQD_pos[ii];
    getInterpolatedDensity(RR, zp, skind, prf, fp
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			   , v2fp, v4fp
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			   );
    getInterpolatedDensity(RR, zm, skind, prf, fm
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			   , v2fm, v4fm
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			   );

    for(int kk = 0; kk < skind; kk++){
      sum[kk] += weight * (fm[kk] + fp[kk]);
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      v2f[kk] += weight * (v2fm[kk] + v2fp[kk]);
      v4f[kk] += weight * (v4fm[kk] + v4fp[kk]);
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    }/* for(int kk = 0; kk < skind; kk++){ */
  }/* for(int ii = (num >> 1) - 1; ii >= 0; ii--){ */

  for(int kk = 0; kk < skind; kk++){
    sum[kk] *= mns;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    v2f[kk] *= mns;
    v4f[kk] *= mns;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }/* for(int kk = 0; kk < skind; kk++){ */
}

/**
 * @fn calcColumnDensityProfile
 *
 * @brief Calculate column density profile.
 *
 * @param (skind) number of spherical symmetric components
 * @return (prf) radial profile of the components
 * @param (logrmax) rmax in logarithmic space
 * @param (cfg) physical properties of the component(s)
 *
 * @sa gaussQuadVertical
 */
void calcColumnDensityProfile(const int skind, profile **prf, const double logrmax, profile_cfg *cfg)
{
  __NOTE__("%s\n", "start");


  double rs = DBL_MAX;
  for(int ii = 0; ii < skind; ii++)
    if( cfg[ii].kind != CENTRALBH )
      rs = fmin(rs, cfg[ii].rs);
  const double Rmax = pow(10.0, logrmax);

#pragma omp parallel
  {
    double *sum;    sum = (double *)malloc(skind * sizeof(double));    if( sum == NULL ){      __KILL__(stderr, "ERROR: failure to allocate sum\n");    }
    double *tfp;    tfp = (double *)malloc(skind * sizeof(double));    if( tfp == NULL ){      __KILL__(stderr, "ERROR: failure to allocate tfp\n");    }
    double *tfm;    tfm = (double *)malloc(skind * sizeof(double));    if( tfm == NULL ){      __KILL__(stderr, "ERROR: failure to allocate tfm\n");    }

#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    double *v2f;    v2f = (double *)malloc(skind * sizeof(double));    if( v2f == NULL ){      __KILL__(stderr, "ERROR: failure to allocate v2f\n");    }
    double *v2p;    v2p = (double *)malloc(skind * sizeof(double));    if( v2p == NULL ){      __KILL__(stderr, "ERROR: failure to allocate v2p\n");    }
    double *v2m;    v2m = (double *)malloc(skind * sizeof(double));    if( v2m == NULL ){      __KILL__(stderr, "ERROR: failure to allocate v2m\n");    }

    double *v4f;    v4f = (double *)malloc(skind * sizeof(double));    if( v4f == NULL ){      __KILL__(stderr, "ERROR: failure to allocate v4f\n");    }
    double *v4p;    v4p = (double *)malloc(skind * sizeof(double));    if( v4p == NULL ){      __KILL__(stderr, "ERROR: failure to allocate v4p\n");    }
    double *v4m;    v4m = (double *)malloc(skind * sizeof(double));    if( v4m == NULL ){      __KILL__(stderr, "ERROR: failure to allocate v4m\n");    }

    double *numer;    numer = (double *)malloc(skind * sizeof(double));    if( numer == NULL ){      __KILL__(stderr, "ERROR: failure to allocate numer\n");    }
    double *denom;    denom = (double *)malloc(skind * sizeof(double));    if( denom == NULL ){      __KILL__(stderr, "ERROR: failure to allocate denom\n");    }
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)

#pragma omp for schedule(dynamic, 8) nowait
    for(int ii = 0; ii < NRADBIN; ii += SKIP_INTERVAL_FOR_COLUMN_DENSITY){
      /** initialization */
      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].Sigma = 0.0;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	numer[kk] = denom[kk] = 0.0;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      }/* for(int kk = 0; kk < skind; kk++){ */

      const double RR = prf[0][ii].rad;
      const double R0 = fmin(RR, rs);
      const double R2 = fmax(RR, rs);
      const double R1 = 0.5 * (R0 + R2);

      gaussQuadVertical(RR, NINTBIN,  0.0     ,          R0, skind, prf, sum, tfm, tfp
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			, v2f, v2m, v2p, v4f, v4m, v4p
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			);
      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].Sigma += sum[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	numer[kk] += v4f[kk];
	denom[kk] += v2f[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      }/* for(int kk = 0; kk < skind; kk++){ */

      gaussQuadVertical(RR, NINTBIN,	    R0,		 R1, skind, prf, sum, tfm, tfp
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			, v2f, v2m, v2p, v4f, v4m, v4p
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			);
      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].Sigma += sum[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	numer[kk] += v4f[kk];
	denom[kk] += v2f[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      }/* for(int kk = 0; kk < skind; kk++){ */

      gaussQuadVertical(RR, NINTBIN,	    R1,		 R2, skind, prf, sum, tfm, tfp
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			, v2f, v2m, v2p, v4f, v4m, v4p
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			);
      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].Sigma += sum[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	numer[kk] += v4f[kk];
	denom[kk] += v2f[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      }/* for(int kk = 0; kk < skind; kk++){ */

      gaussQuadVertical(RR, NINTBIN,	    R2,	 2.0 *	 R2, skind, prf, sum, tfm, tfp
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			, v2f, v2m, v2p, v4f, v4m, v4p
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			);
      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].Sigma += sum[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	numer[kk] += v4f[kk];
	denom[kk] += v2f[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      }/* for(int kk = 0; kk < skind; kk++){ */

      gaussQuadVertical(RR, NINTBIN,  2.0 * R2, 10.0 *	 R2, skind, prf, sum, tfm, tfp
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			, v2f, v2m, v2p, v4f, v4m, v4p
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			);
      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].Sigma += sum[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	numer[kk] += v4f[kk];
	denom[kk] += v2f[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      }/* for(int kk = 0; kk < skind; kk++){ */

      gaussQuadVertical(RR, NINTBIN, 10.0 * R2,	 2.0 * Rmax, skind, prf, sum, tfm, tfp
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			, v2f, v2m, v2p, v4f, v4m, v4p
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
			);
      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].Sigma += sum[kk];
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	numer[kk] += v4f[kk];
	denom[kk] += v2f[kk];
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      }/* for(int kk = 0; kk < skind; kk++){ */


      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].Sigma *= 2.0;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	prf[kk][ii].slos = sqrt(numer[kk] / (3.0 * denom[kk]));
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      }/* for(int kk = 0; kk < skind; kk++){ */
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */

    free(sum);    free(tfp);    free(tfm);
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
    free(v2f);    free(v2p);    free(v2m);
    free(v4f);    free(v4p);    free(v4m);
    free(numer);    free(denom);
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
  }

#   if  SKIP_INTERVAL_FOR_COLUMN_DENSITY != 1
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN - SKIP_INTERVAL_FOR_COLUMN_DENSITY; ii += SKIP_INTERVAL_FOR_COLUMN_DENSITY)
    for(int kk = 0; kk < skind; kk++){
      const double S0 = prf[kk][ii                                   ].Sigma;
      const double S1 = prf[kk][ii + SKIP_INTERVAL_FOR_COLUMN_DENSITY].Sigma;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      const double s0 = prf[kk][ii                                   ].slos;
      const double s1 = prf[kk][ii + SKIP_INTERVAL_FOR_COLUMN_DENSITY].slos;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)

      double Slope = (S1 - S0) / (double)SKIP_INTERVAL_FOR_COLUMN_DENSITY;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      double slope = (s1 - s0) / (double)SKIP_INTERVAL_FOR_COLUMN_DENSITY;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      for(int jj = 1; jj < SKIP_INTERVAL_FOR_COLUMN_DENSITY; jj++){
	prf[kk][ii + jj].Sigma = S0 + Slope * (double)jj;
#   if  defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
	prf[kk][ii + jj].slos  = s0 + slope * (double)jj;
#endif//defined(MAKE_VELOCITY_DISPERSION_PROFILE) && !defined(USE_OSIPKOV_MERRITT_METHOD)
      }/* for(int jj = 1; jj < SKIP_INTERVAL_FOR_COLUMN_DENSITY; jj++){ */
    }/* for(int kk = 0; kk < skind; kk++){ */
#endif//SKIP_INTERVAL_FOR_COLUMN_DENSITY != 1


  __NOTE__("%s\n", "end");
}
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
#endif//MAKE_COLUMN_DENSITY_PROFILE
