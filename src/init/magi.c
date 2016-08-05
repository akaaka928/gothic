/*************************************************************************\
 *                                                                       *
                  last updated on 2016/08/02(Tue) 17:37:31
 *                                                                       *
 *    MAGI: "MAny-component Galactic Initial-conditions" generator       *
 *    Making Initial Condition Code of N-body Simulation                 *
 *       Based on distribution function given by Eddington's formula     *
 *       arbitrary spherical symmetric and isotropic DF                  *
 *       extension in multi-components is based on Kazantzidis+06        *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
/* conversion from physical unit to computational unit must be performed internally */
//-------------------------------------------------------------------------
#define USE_SZIP_COMPRESSION
//-------------------------------------------------------------------------
#define SHIFT_CENTER_PER_COMPONENT
//-------------------------------------------------------------------------
/* #define RESET_ROTATION_AXIS */
#   if  !defined(RESET_ROTATION_AXIS) && defined(SHIFT_CENTER_PER_COMPONENT)
#define          RESET_ROTATION_AXIS
#endif//!defined(RESET_ROTATION_AXIS) && defined(SHIFT_CENTER_PER_COMPONENT)
//-------------------------------------------------------------------------
#define ERFC_SMOOTHING
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
#       include <hdf5.h>
#       include <hdf5lib.h>
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <name.h>
#include <myutil.h>
#include <constants.h>
#include <timer.h>
#include <rotate.h>
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#include "../misc/allocate.h"
//-------------------------------------------------------------------------
#include "../file/io.h"
//-------------------------------------------------------------------------
#include "magi.h"
#include "king.h"
#include "table.h"
#include "abel.h"
#include "disk.h"
//-------------------------------------------------------------------------
extern const real newton;
extern const double     mass_astro2com,     mass2astro;extern const char        mass_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double   length_astro2com,   length2astro;extern const char      length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double     time_astro2com,     time2astro;extern const char        time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double velocity_astro2com, velocity2astro;extern const char    velocity_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                       energy2astro;extern const char      energy_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                  col_density2astro;extern const char col_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                      density2astro;extern const char     density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                      senergy2astro;extern const char     senergy_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
//-------------------------------------------------------------------------
gsl_rng *GSLRand;
#define UNIRAND_DBL ((double)gsl_rng_uniform(GSLRand))
#define UNIRAND     ((real)gsl_rng_uniform(GSLRand))
#define RANDVAL     (TWO * (UNIRAND) - UNITY)
//-------------------------------------------------------------------------
double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void setDensityProfilePlummer(profile *prf, const double rs)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double rad_rs = prf[ii].rad * rsinv;
    const double x2_pls_1 = 1.0 + rad_rs * rad_rs;
    const double     inv = 1.0 / x2_pls_1;
    const double sqrtinv = sqrt(inv);
    const double pow_m5_2 = inv * inv * sqrtinv;
    //---------------------------------------------------------------------
    prf[ii].  rho     = pow_m5_2;
    prf[ii]. drho_dr  = - 5.0 * rad_rs * pow_m5_2 * inv * rsinv;
    prf[ii].d2rho_dr2 = 5.0 * pow_m5_2 * inv * (7.0 * rad_rs * rad_rs * inv - 1.0) * rsinv * rsinv;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileBurkert(profile *prf, const double rs)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double rad_rs      = prf[ii].rad * rsinv;
    const double rad_rs_p1   = rad_rs + 1.0;
    const double rad2_rs2    = rad_rs * rad_rs;
    const double rad2_rs2_p1 = 1.0 + rad2_rs2;
    const double inv = 1.0 / (rad_rs_p1 * rad2_rs2_p1);
    //---------------------------------------------------------------------
    prf[ii].  rho     = inv;
    prf[ii]. drho_dr  = - (1.0 + rad_rs * (2.0 + 3.0 * rad_rs)) * inv * inv * rsinv;
    prf[ii].d2rho_dr2 = 4.0 * rad2_rs2 * (4.0 * rad_rs + 3.0 * rad2_rs2_p1) * inv * inv * inv * rsinv * rsinv;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileHernquist(profile *prf, const double rs)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double rad_rs = prf[ii].rad * rsinv;
    const double rad_rs_p1 = rad_rs + 1.0;
    const double tmp = rad_rs * rad_rs_p1;
    const double inv = 1.0 / (tmp * rad_rs_p1);
    //---------------------------------------------------------------------
    prf[ii].  rho     = 1.0 * inv * inv * tmp;
    prf[ii]. drho_dr  = - (1.0 + 4.0 * rad_rs) * inv * inv * rsinv;
    prf[ii].d2rho_dr2 = 2.0 * (1.0 + 5.0 * rad_rs * (1.0 + 2.0 * rad_rs)) * inv * inv * inv * rad_rs_p1 * rsinv * rsinv;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileNFW(profile *prf, const double rs)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double rad_rs = prf[ii].rad * rsinv;
    const double rad_rs_p1 = rad_rs + 1.0;
    const double tmp = rad_rs * rad_rs_p1;
    const double inv = 1.0 / (tmp * rad_rs_p1);
    //---------------------------------------------------------------------
    prf[ii].  rho     = inv;
    prf[ii]. drho_dr  = - (1.0 + 3.0 * rad_rs) * inv * inv * rad_rs_p1 * rsinv;
    prf[ii].d2rho_dr2 = 2.0 * (2.0 * rad_rs * (3.0 * rad_rs + 2.0) + 1.0) * inv * inv * inv * rad_rs_p1 * rad_rs_p1 * rsinv * rsinv;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileMoore(profile *prf, const double rs)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double rad_rs   = prf[ii].rad * rsinv;
    const double     inv  = 1.0 / (rad_rs * (1.0 + rad_rs));
    const double sqrtinv  = sqrt(inv);
    const double pow_m5_2 = inv * inv * sqrtinv;
    //---------------------------------------------------------------------
    prf[ii].  rho     =  1.0 * inv * sqrtinv;
    prf[ii]. drho_dr  = -1.5 * (1.0 + 2.0 * rad_rs) * pow_m5_2 * rsinv;
    prf[ii].d2rho_dr2 =  0.75 * (5.0 + 16.0 * rad_rs * (1.0 + rad_rs)) * pow_m5_2 * inv * rsinv * rsinv;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileEinasto(profile *prf, const double rs, const double alpha)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double rad_rs = prf[ii].rad * rsinv;
    const double xalpha = pow(rad_rs, alpha);
    const double xinv = 1.0 / rad_rs;
    //---------------------------------------------------------------------
    const double rho = pow(M_E, (-2.0 / alpha) * (xalpha - 1.0));
    prf[ii].  rho     = rho;
    prf[ii]. drho_dr  = -2.0 * xalpha * xinv * rho * rsinv;
    prf[ii].d2rho_dr2 =  2.0 * xalpha * xinv * xinv * rho * rsinv * rsinv * (2.0 * xalpha + 1.0 - alpha);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileAppKing(profile *prf, const double rs, const double rt)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
  const double rsinv = 1.0 / rs;
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double rad_rs = prf[ii].rad * rsinv;
    const double CC = 1.0 / sqrt(1.0 + rt * rsinv * rt * rsinv);
    const double inv = 1.0 / (1.0 + rad_rs * rad_rs);
    const double sqrtinv = sqrt(inv);
    //---------------------------------------------------------------------
    prf[ii].  rho     = (prf[ii].rad < rt) ? ((sqrtinv - CC) * (sqrtinv - CC)) : (0.0);
    prf[ii]. drho_dr  = (prf[ii].rad < rt) ? (-2.0 * (sqrtinv - CC) * rad_rs * inv * sqrtinv * rsinv) : (0.0);
    prf[ii].d2rho_dr2 = (prf[ii].rad < rt) ? ((2.0 * inv * inv * (4.0 * rad_rs * rad_rs * inv - 1.0) + 2.0 * CC * sqrtinv * inv * (1.0 - 3.0 * rad_rs * rad_rs * inv)) * rsinv * rsinv) : (0.0);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileTwoPower(profile *prf, const double rs, const double alpha, const double beta, const double gam)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
  const double rsinv = 1.0 / rs;
  const double binv  = 1.0 / beta;
  const double agbi  = (alpha - gam) * binv;
  const double b1g2  = 1.0 + beta + 2.0 * gam;
  const double a1    = 1.0 + alpha;
  const double g1    = 1.0 + gam;
  const double b1    = 1.0 - beta;
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double xx    = prf[ii].rad * rsinv;
    const double xinv  = 1.0 / xx;
    const double xa    = pow(xinv, alpha);
    const double xb    = pow(  xx, beta);
    const double xbp1  = 1.0 + xb;
    const double xbp1i = 1.0 / xbp1;
    const double base  = xa * pow(xbp1, agbi);
    //---------------------------------------------------------------------
    prf[ii].  rho     =                  base;
    prf[ii]. drho_dr  = -rsinv         * base * xinv        * xbp1i         * (alpha                    + gam * xb);
    prf[ii].d2rho_dr2 =  rsinv * rsinv * base * xinv * xinv * xbp1i * xbp1i * (alpha * (a1 + b1g2 * xb) + gam * xb * (b1 + g1 * xb));
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileTriPower(profile *prf, const double rin, const double rout, const double alp, const double bet, const double gam, const double del, const double eps)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
  /* At first, set alpha-beta-gamma profile */
  setDensityProfileTwoPower(prf, rin, alp, bet, gam);
  const double rsinv = 1.0 / rout;
  const double gmeod = (gam - eps) / del;
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double yy   = prf[ii].rad * rsinv;
    const double yd   = pow(yy, del);
    const double ydp1 = 1.0 + yd;
    const double inv  = 1.0 / (yy * ydp1);
    //---------------------------------------------------------------------
    const double base0 = pow(ydp1, gmeod);
    const double base1 = (gam - eps) * yd * rsinv * inv * base0;
    const double base2 =                    rsinv * inv * base1 * ((del - 1.0) + (gam - eps - 1.0) * yd);
    //---------------------------------------------------------------------
    prf[ii].d2rho_dr2 = prf[ii].d2rho_dr2 * base0 + 2.0 * prf[ii].drho_dr * base1 + prf[ii].rho * base2;
    prf[ii]. drho_dr  = prf[ii]. drho_dr  * base0                                 + prf[ii].rho * base1;
    prf[ii].  rho     = prf[ii].  rho     * base0;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileAppLoweredEvans(profile *prf, const double rs, const double alpha, const double rc, const double beta, const double rt, const double invdelta)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
  const double rsinv  = 1.0 / rs;
  const double rcinv  = 1.0 / rc;
  const double rc2 = rc * rc;
  const double  betax2 = 2.0 * beta;
  const double alphax2 = 2.0 * alpha;
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double rad = prf[ii].rad;
    const double rad_rs = rad * rsinv;
    const double rad_rc = rad * rcinv;
    const double inv_1p_rad_rs  = 1.0 / (1.0 + rad_rs);
    const double inv_1p_rad_rc2 = 1.0 / (1.0 + rad_rc * rad_rc);
    const double inv_1p_rad_rs_alpha = pow(inv_1p_rad_rs , alpha);
    const double inv_1p_rad_rc2_beta = pow(inv_1p_rad_rc2, beta);
    const double rho = inv_1p_rad_rs_alpha * inv_1p_rad_rc2_beta;
    //---------------------------------------------------------------------
    const double r2 = rad * rad;
    const double inv = 1.0 / ((rs + rad) * (rc2 + r2));
    const double val = alpha * (rc2 + r2) + betax2 * (rs + rad) * rad;
    const double drho = -inv * val * rho;
    //---------------------------------------------------------------------
    const double sinv = 1.0 / (rs + rad);
    const double rinv = 1.0 / (rc2 + r2);
    const double d2rho = rho * ((1.0 + alpha) * alpha * sinv * sinv + betax2 * (rc2 + rad * ((1.0 + betax2) * rad + alphax2 * (rc2 + r2) * sinv)) * rinv * rinv);
    //---------------------------------------------------------------------
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
    //---------------------------------------------------------------------
    prf[ii].  rho     = rho * smooth;
    prf[ii]. drho_dr  = rho * dsmooth + drho * smooth;
    prf[ii].d2rho_dr2 = rho * d2smooth + 2.0 * drho * dsmooth + d2rho * smooth;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < 4 + NRADBIN; ii += (NRADBIN / N_PRINT_LINES_ASCII))
    fprintf(stderr, "%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].drho_dr, prf[ii].d2rho_dr2);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setContributionByCentralBH(profile *prf, const profile_cfg cfg)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set density profile */
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    prf[ii].  rho     = 0.0;
    prf[ii]. drho_dr  = 0.0;
    prf[ii].d2rho_dr2 = 0.0;
    //---------------------------------------------------------------------
    prf[ii].enc =                  cfg.Mtot;
    prf[ii].psi = (double)newton * cfg.Mtot / prf[ii].rad;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void integrateDensityProfile(profile *prf, const double logrbin, const double Mtot, const bool cutoff, const double redge, const double width)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set cutoff radius if necessary */
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < 4 + NRADBIN; ii += (NRADBIN >> 6))
    fprintf(stderr, "%e\t%e\n", prf[ii].rad, prf[ii].rho);
#endif
  //-----------------------------------------------------------------------
  if( cutoff == true ){
    //---------------------------------------------------------------------
    const double invThick = 1.0 / width;
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < 4 + NRADBIN; ii++){
      //-------------------------------------------------------------------
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
      //-------------------------------------------------------------------
      prf[ii].d2rho_dr2  = d2hdx2 * prf[ii].rho + 2.0 * dhdx * prf[ii].drho_dr + smooth * prf[ii].d2rho_dr2;
      prf[ii]. drho_dr   =  dhdx  * prf[ii].rho +     smooth * prf[ii].drho_dr;
      prf[ii].  rho     *= smooth;
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
    //---------------------------------------------------------------------
  }/* if( cutoff == true ){ */
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < 4 + NRADBIN; ii += (NRADBIN >> 6))
    fprintf(stderr, "%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].drho_dr, prf[ii].d2rho_dr2);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate enclosed mass and potential (internal part) */
  //-----------------------------------------------------------------------
  double Menc[2];
  Menc[0] = prf[0].rad * prf[0].rad * prf[0].rad * prf[0].rho;
  Menc[1] = prf[1].rad * prf[1].rad * prf[1].rad * prf[1].rho;
  prf[0].enc =              Menc[0] / 3.0;
  prf[1].enc = prf[0].enc + (prf[1].rad - prf[0].rad) * prf[1].rad * prf[0].rad * 0.5 * (prf[0].rho + prf[1].rho);
  Menc[0] += Menc[1] * 4.0;
  //-----------------------------------------------------------------------
  for(int ii = 2; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    const double mass = prf[ii].rad * prf[ii].rad * prf[ii].rad * prf[ii].rho;
    const int idx = ii & 1;
    //---------------------------------------------------------------------
    prf[ii].enc = prf[idx].enc + (mass + Menc[idx]) * logrbin * M_LN10 / 3.0;
    prf[ii].psi = prf[ii].enc / prf[ii].rad;
    //---------------------------------------------------------------------
    Menc[0] += mass * (double)(1 << (1 + (idx    )));
    Menc[1] += mass * (double)(1 << (1 + (idx ^ 1)));
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate potential (external part) */
  //-----------------------------------------------------------------------
  double Pext[2];
  Pext[0] = prf[NRADBIN + 3].rad * prf[NRADBIN + 3].rad * prf[NRADBIN + 3].rho;
  Pext[1] = prf[NRADBIN + 2].rad * prf[NRADBIN + 2].rad * prf[NRADBIN + 2].rho;
  double Pini[2];
  Pini[0] =           0.0;
  Pini[1] = Pini[0] + (prf[NRADBIN + 3].rad - prf[NRADBIN + 2].rad) * sqrt(prf[NRADBIN + 2].rad * prf[NRADBIN + 3].rad) * 0.5 * (prf[NRADBIN + 2].rho + prf[NRADBIN + 3].rho);
  Pext[0] += Pext[1] * 4.0;
  //-----------------------------------------------------------------------
  for(int ii = NRADBIN + 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double psi = prf[ii].rad * prf[ii].rad * prf[ii].rho;
    const int idx = (int)((ii + 1) & 1);
    //---------------------------------------------------------------------
    prf[ii].psi += Pini[idx] + (psi + Pext[idx]) * logrbin * M_LN10 / 3.0;
    //---------------------------------------------------------------------
    Pext[0] += psi * (double)(1 << (1 + (idx    )));
    Pext[1] += psi * (double)(1 << (1 + (idx ^ 1)));
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* multiply overall factors */
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    prf[ii].enc *= 4.0 * M_PI;
    prf[ii].psi *= 4.0 * M_PI * (double)newton;
  }
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < 4 + NRADBIN; ii += (NRADBIN >> 6))
    fprintf(stderr, "%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].enc, prf[ii].psi);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set appropriate unit system */
  //-----------------------------------------------------------------------
  const double Mscale = Mtot / prf[NRADBIN + 3].enc;
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    prf[ii].  rho     *= Mscale;
    prf[ii].  enc     *= Mscale;
    prf[ii].  psi     *= Mscale;
    prf[ii]. drho_dr  *= Mscale;
    prf[ii].d2rho_dr2 *= Mscale;
  }
  //-----------------------------------------------------------------------
#if 0
  /* for(int ii = 0; ii < 4 + NRADBIN; ii += 1000) */
  for(int ii = 1; ii < 4 + NRADBIN - 1; ii++)
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].drho_dr, prf[ii].d2rho_dr2, (prf[ii + 1].rho - prf[ii - 1].rho) / (prf[ii + 1].rad - prf[ii - 1].rad), (prf[ii + 1].drho_dr - prf[ii - 1].drho_dr) / (prf[ii + 1].rad - prf[ii - 1].rad));
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void findIdx(const double psi, profile *prf, int *ll, int *rr)
{
  //-----------------------------------------------------------------------
  bool bisection = true;
  *ll = 2;
  *rr = 1 + NRADBIN;
  //-----------------------------------------------------------------------
  if( bisection == true )    if( fabs(prf[*ll].psi_tot - psi) / psi < DBL_EPSILON ){      bisection = false;      *rr = (*ll) + 1;    }
  if( bisection == true )    if( fabs(prf[*rr].psi_tot - psi) / psi < DBL_EPSILON ){      bisection = false;      *ll = (*rr) - 1;    }
  //-----------------------------------------------------------------------
  while( bisection ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)(*ll) + (uint)(*rr)) >> 1;
    //---------------------------------------------------------------------
    if( (prf[cc].psi_tot - psi) * (prf[*ll].psi_tot - psi) <= 0.0 )      *rr = (int)cc;
    else                                                                 *ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + (*ll)) == (*rr) )      break;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void getEddingtonFormula(const double ene, const double psi, const int kind, profile **prf, double *val)
{
  //-----------------------------------------------------------------------
  int ll, rr;
  findIdx(psi, prf[0], &ll, &rr);
  //-----------------------------------------------------------------------
  /* based on linear interpolation */
  const double r_dr = (psi - prf[0][ll].psi_tot) / (prf[0][rr].psi_tot - prf[0][ll].psi_tot);
  const double rad = (1.0 - r_dr) * prf[0][ll].rad     + r_dr * prf[0][rr].rad;
  const double rho = (1.0 - r_dr) * prf[0][ll].rho_tot + r_dr * prf[0][rr].rho_tot;
  const double enc = (1.0 - r_dr) * prf[0][ll].enc_tot + r_dr * prf[0][rr].enc_tot;
  const double fac = 2.0 / rad - 4.0 * M_PI * rho * rad * rad / enc;
  //-----------------------------------------------------------------------
  double common = M_1_PI * rad * rad / ((double)newton * enc);
  common *= common;
  common *= (0.5 * M_SQRT1_2 / sqrt(ene - psi));
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < kind; kk++){
    //---------------------------------------------------------------------
    /* based on linear interpolation */
    const double  drho_dr  = (1.0 - r_dr) * prf[kk][ll]. drho_dr  + r_dr * prf[kk][rr]. drho_dr;
    const double d2rho_dr2 = (1.0 - r_dr) * prf[kk][ll].d2rho_dr2 + r_dr * prf[kk][rr].d2rho_dr2;
    //---------------------------------------------------------------------
    val[kk] = common * (d2rho_dr2 + drho_dr * fac);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void gaussQuad1dEddington(const int num, const double xmin, const double xmax,
					const int kind, profile **prf, double *sum, double *fm, double *fp)
{
  //-----------------------------------------------------------------------
  const double mns = 0.5 * (xmax - xmin);
  const double pls = 0.5 * (xmax + xmin);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int kk = 0; kk < kind; kk++)    sum[kk] = 0.0;
  //-----------------------------------------------------------------------
  if( num & 1 ){
    const double weight =             gsl_gaussQD_weight[(num >> 1)];
    const double  value = pls + mns * gsl_gaussQD_pos   [(num >> 1)];
    getEddingtonFormula(xmax, value, kind, prf, fm);
    for(int kk = 0; kk < kind; kk++)
      sum[kk] = weight * fm[kk];
  }/* if( num & 1 ){ */
  //-----------------------------------------------------------------------
  for(int ii = (num >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double weight = gsl_gaussQD_weight[ii];
    const double xp = pls + mns * gsl_gaussQD_pos[ii];
    const double xm = pls - mns * gsl_gaussQD_pos[ii];
    getEddingtonFormula(xmax, xp, kind, prf, fp);
    getEddingtonFormula(xmax, xm, kind, prf, fm);
    //---------------------------------------------------------------------
    for(int kk = 0; kk < kind; kk++)
      sum[kk] += weight * (fm[kk] + fp[kk]);
    //---------------------------------------------------------------------
  }/* for(int ii = (num >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < kind; kk++)
    sum[kk] *= mns;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void integrateEddingtonFormula(const int skind, profile **prf, dist_func **fene)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set integration range */
  //-----------------------------------------------------------------------
  /* find maximum of r^2 rho and truncated radius for the particle distribution */
  double fmax = 0.0;
  int   iout = 1 + NRADBIN;
  for(int ii = 2; ii < 2 + NRADBIN; ii++){
    double floc = prf[0][ii].rad * prf[0][ii].rad * prf[0][ii].rho_tot;
    if(                                        floc > fmax ){      fmax = floc;    }
  }
  //-----------------------------------------------------------------------
  const double Emax = prf[0][2   ].psi_tot;
  const double Emin = prf[0][iout].psi_tot;
  const double Ebin = (Emax - Emin) / (double)(NENEBIN - 1);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  NDIVIDE_GAUSSQD
  //-----------------------------------------------------------------------
  const int nsub = NENEBIN / NDIVIDE_GAUSSQD;
  int head = 0;
  double Ein = 0.0;
  static double sub[NKIND_MAX];
  for(int ii = 0; ii < skind; ii++)    sub[ii] = 0.0;
  //-----------------------------------------------------------------------
  for(int iter = 0; iter < NDIVIDE_GAUSSQD; iter++){
    const int tail = head + nsub;
#pragma omp parallel for schedule(dynamic, 16)
    for(int ii = head; ii < tail; ii++){
      //-------------------------------------------------------------------
      double sum[NKIND_MAX], fm[NKIND_MAX], fp[NKIND_MAX];
      //-------------------------------------------------------------------
      const double ene = Emin + Ebin * (double)ii;
      gaussQuad1dEddington(NINTBIN, Ein, ene, skind, prf, sum, fm, fp);
      //-------------------------------------------------------------------
      for(int kk = 0; kk < skind; kk++){
	sum[kk] += sub[kk];
	fene[kk][ii].ene = (real)ene;
	fene[kk][ii].val = (sum[kk] >= 0.0) ? (real)sum[kk] : ZERO;
      }/* for(int kk = 0; kk < skind; kk++){ */
      //-------------------------------------------------------------------
    }/* for(int ii = head; ii < tail; ii++){ */
    //---------------------------------------------------------------------
    head = tail;
    for(int ii = 0; ii < skind; ii++)
      sub[ii] = fene[ii][head - 1].val;
    Ein = fene[0][head - 1].ene;
    //---------------------------------------------------------------------
  }/* for(int iter = 0; iter < NDIVIDE_GAUSSQD; iter++){ */
#pragma omp parallel for schedule(dynamic, 16)
  for(int ii = 0; ii < NENEBIN; ii++){
    //---------------------------------------------------------------------
    double sum[NKIND_MAX], fm[NKIND_MAX], fp[NKIND_MAX];
    //---------------------------------------------------------------------
    const double ene = Emin + Ebin * (double)ii;
    gaussQuad1dEddington(NINTBIN, 0.0, ene, skind, prf, sum, fm, fp);
    //---------------------------------------------------------------------
    for(int kk = 0; kk < skind; kk++){
      fene[kk][ii].ene = (real)ene;
      fene[kk][ii].val = (sum[kk] >= 0.0) ? (real)sum[kk] : ZERO;
    }/* for(int kk = 0; kk < skind; kk++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NENEBIN; ii++){ */
  //-----------------------------------------------------------------------
#else///NDIVIDE_GAUSSQD
  //-----------------------------------------------------------------------
#pragma omp parallel for schedule(dynamic, 16)
  for(int ii = 0; ii < NENEBIN; ii++){
    //---------------------------------------------------------------------
    double sum[NKIND_MAX], fm[NKIND_MAX], fp[NKIND_MAX];
    //---------------------------------------------------------------------
    const double ene = Emin + Ebin * (double)ii;
    gaussQuad1dEddington(NINTBIN, 0.0, ene, skind, prf, sum, fm, fp);
    //---------------------------------------------------------------------
    for(int kk = 0; kk < skind; kk++){
      fene[kk][ii].ene = (real)ene;
#if 1
      fene[kk][ii].val = (sum[kk] >= 0.0) ? (real)sum[kk] : ZERO;
#else
      fene[kk][ii].val = (real)sum[kk];
#endif
    }/* for(int kk = 0; kk < skind; kk++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NENEBIN; ii++){ */
  //-----------------------------------------------------------------------
#endif//NDIVIDE_GAUSSQD
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* based on Ciotti & Bertin (1999), A&A, 352, 447-451: Eq.(25) */
static inline double getAsymptoticSersicScale(const double nn){  return (2.0 * nn + (-1.0 + 2.0 * (2.0 + 23.0 / (63.0 * nn)) / (135.0 * nn)) / 3.0);}
//-------------------------------------------------------------------------
static inline void writeProfileCfgFormat(char *filename, const profile_cfg cfg)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  fprintf(stderr, "ERROR: data written in \"%s\" does not match with format of specified model id \"%d\"\n", filename, cfg.kind);
  fprintf(stderr, "Expected format is below:\n");
  fprintf(stderr, "\tMtot<real>: total mass of the model in astrophysical units\n");
  if( cfg.kind != CENTRALBH )
    fprintf(stderr, "\trs<real>: scale length of the model in astrophysical units\n");
  //-----------------------------------------------------------------------
  /* some models require more information */
  if( cfg.kind == KING )
    fprintf(stderr, "\tW0<real>: dimensionless King parameter at the center\n");
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
  if( (cfg.kind == EXP_DISK) || (cfg.kind == SERSIC) ){
    if( cfg.kind == SERSIC )
      fprintf(stderr, "\tn_sersic<real>: Sersic index\n");
    fprintf(stderr, "\tRt<real> Rt_width<real>: cutoff radius and width of the disk mid-plane density in horizontal direction in astrophysical units, respectively\n");
    fprintf(stderr, "\tzd<real>: scale height of isothermal disk in the vertical direction in astrophysical units\n");
    fprintf(stderr, "\tsigmaR0<real> frac<real>: velocity dispersion in radial direction at the center in the vertical direction in astrophysical units, perpendicular velocity dispersion over circular velocity\n");
    fprintf(stderr, "\t\tif the inputted sigmaR0 is negative, then the default value (= sigma_z(R = 0)) is substituted\n");
    fprintf(stderr, "\t\tif USE_ORIGINAL_VDISP_ESTIMATOR defined in src/init/disk_polar.h is ON, then frac is used; if that is OFF, then sigmaR0 is used.\n");
  }/* if( (cfg.kind == EXP_DISK) || (cfg.kind == SERSIC) ){ */
  //-----------------------------------------------------------------------
  /* information about density cutoff */
  if( cfg.kind != CENTRALBH ){
    fprintf(stderr, "\tcutoff<bool>: set explicit density cutoff (1) or no (0)\n");
    if( cfg.cutoff )
      fprintf(stderr, "\trc<real> rc_width<real>: cutoff radius and width of the density cutoff in astrophysical units, respectively\n");
  }/* if( cfg.kind != CENTRALBH ) */
  //-----------------------------------------------------------------------
  __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void readProfileCfg(char *fcfg, int *unit, int *kind, profile_cfg **cfg)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* read global settings */
  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s", CFGFOLDER, fcfg);
  //-----------------------------------------------------------------------
  fp = fopen(filename, "r");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  int checker = 1;
  //-----------------------------------------------------------------------
  /* read the specified unit system and set it */
  checker &= (1 == fscanf(fp, "%d", unit));
  setPhysicalConstantsAndUnitSystem(*unit, 1);
  /* read # of components */
  checker &= (1 == fscanf(fp, "%d", kind));
  //-----------------------------------------------------------------------
  *cfg = (profile_cfg *)malloc(sizeof(profile_cfg) * (*kind));  if( *cfg == NULL ){    __KILL__(stderr, "ERROR: failure to allocate cfg");  }
  for(int ii = 0; ii < *kind; ii++)
    checker &= (4 == fscanf(fp, "%d\t%s\t%d\t%zu", &(*cfg)[ii].kind, (*cfg)[ii].file, &(*cfg)[ii].forceNum, &(*cfg)[ii].num));
    /* checker &= (2 == fscanf(fp, "%d\t%s", &(*cfg)[ii].kind, (*cfg)[ii].file)); */
  //-----------------------------------------------------------------------
  fclose(fp);
  if( !checker ){
    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
  }/* if( !checker ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* read individual settings */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < *kind; ii++){
    //---------------------------------------------------------------------
    sprintf(filename, "%s/%s", CFGFOLDER, (*cfg)[ii].file);
    //---------------------------------------------------------------------
    fp = fopen(filename, "r");
    if( fp == NULL ){
      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
    }/* fp = fopen(filename, "r"); */
    checker = 1;
    //---------------------------------------------------------------------
    checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].Mtot));    (*cfg)[ii].Mtot *= mass_astro2com;
    if( (*cfg)[ii].kind != CENTRALBH ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].rs));      (*cfg)[ii].rs *= length_astro2com;
    }/* if( (*cfg)[ii].kind != CENTRALBH ){ */
    else
      if( ((*cfg)[ii].forceNum != 1) || ((*cfg)[ii].num != 1) ){
	__KILL__(stderr, "ERROR: number of BH particle must be specified to be unity\n");
      }/* if( ((*cfg)[ii].forceNum != 1) || ((*cfg)[ii].num != 1) ){ */
    //---------------------------------------------------------------------
    /* parameter for King profile */
    if( (*cfg)[ii].kind ==     KING )
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].king_W0));
    if( (*cfg)[ii].kind == APP_KING ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].king_rt));      (*cfg)[ii].king_rt *= length_astro2com;
    }/* if( (*cfg)[ii].kind == APP_KING ){ */
    //---------------------------------------------------------------------
    /* parameter for Einasto profile */
    if( (*cfg)[ii].kind == EINASTO )
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].einasto_alpha));
    //---------------------------------------------------------------------
    /* parameters for two-power model */
    if( (*cfg)[ii].kind == TWO_POWER ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_alpha));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_beta));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_gamma));
    }/* if( (*cfg)[ii].kind == TWO_POWER ){ */
    //---------------------------------------------------------------------
    /* parameters for three-power model */
    if( (*cfg)[ii].kind == TRI_POWER ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].tripower_rout));      (*cfg)[ii].tripower_rout *= length_astro2com;
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_alpha));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_beta));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_gamma));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].tripower_delta));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].tripower_epsilon));
    }/* if( (*cfg)[ii].kind == TRI_POWER ){ */
    //---------------------------------------------------------------------
    /* parameters for approximated lowered Evans model */
    if( (*cfg)[ii].kind == APP_EVANS ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].alevans_alpha));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].alevans_rc));      (*cfg)[ii].alevans_rc *= length_astro2com;
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].alevans_beta));
      checker &= (2 == fscanf(fp, "%le %le", &(*cfg)[ii].alevans_rt, &(*cfg)[ii].alevans_wt));
      (*cfg)[ii].alevans_rt *= length_astro2com;
      (*cfg)[ii].alevans_wt *= length_astro2com;
    }/* if( (*cfg)[ii].kind == APP_EVANS ){ */
    //---------------------------------------------------------------------
    /* parameter for density profile in table form */
    if( (*cfg)[ii].kind == TABLE_RHO )
      checker &= (1 == fscanf(fp, "%s", (*cfg)[ii].table));
    //---------------------------------------------------------------------
    /* parameter for column density profile in table form */
    if( (*cfg)[ii].kind == TABLE_SIG )
      checker &= (1 == fscanf(fp, "%s", (*cfg)[ii].table));
    //---------------------------------------------------------------------
    /* parameter for spherical Serisic profile */
    if( (*cfg)[ii].kind == SPHSERSIC ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].n_sersic));
      (*cfg)[ii].b_sersic = getAsymptoticSersicScale((*cfg)[ii].n_sersic);
    }/* if( (*cfg)[ii].kind == SPHSERSIC ){ */
    //---------------------------------------------------------------------
    /* parameters for projected two-power model */
    if( (*cfg)[ii].kind == SIGTWOPOW ){
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_alpha));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_beta));
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].twopower_gamma));
    }/* if( (*cfg)[ii].kind == SIGTWOPOW ){ */
    //---------------------------------------------------------------------
    /* parameters for disk component */
    if( ((*cfg)[ii].kind == EXP_DISK) || ((*cfg)[ii].kind == SERSIC) ){
      if( (*cfg)[ii].kind == SERSIC ){
	checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].n_sersic));
	(*cfg)[ii].b_sersic = getAsymptoticSersicScale((*cfg)[ii].n_sersic);
      }/* if( (*cfg)[ii].kind == SERSIC ){ */
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].zd));      (*cfg)[ii].zd   *= length_astro2com;
      checker &= (2 == fscanf(fp, "%le %le", &(*cfg)[ii].vdispR0, &(*cfg)[ii].vdisp_frac));      (*cfg)[ii].vdispR0 *= velocity_astro2com;
    }/* if( ((*cfg)[ii].kind == EXP_DISK) || ((*cfg)[ii].kind == SERSIC) ){ */
    //---------------------------------------------------------------------
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
    //---------------------------------------------------------------------
    fclose(fp);
    if( !checker )
      writeProfileCfgFormat(filename, (*cfg)[ii]);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void outputFundamentalInformation
(const int unit, const int kind, const int skind, profile_cfg *cfg, profile **prf, dist_func **df, const ulong Ntot, const real eps, const double snapshotInterval, const double ft, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  FILE *fp;
  char filename[256], date[64];
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* output useful information for multi-component analysis */
  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
  //-----------------------------------------------------------------------
  fprintf(fp, "%d\n", unit);
  fprintf(fp, "%d\t%d\n", kind, skind);
  for(int ii = 0; ii < kind; ii++)
    fprintf(fp, "%zu\n", cfg[ii].num);
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* output fundamental information of the particle distribution */
  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.info.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
  //-----------------------------------------------------------------------
  /* output global settings */
  getPresentDateInStrings(date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Fundamental information about the particle distribution\n");
  fprintf(fp, "\tgenerated on %s", date);
  fprintf(fp, "Physical quantities in Computational and Astrophysical units is listed\n");
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "#############################################################################\n");
  double Mtot = 0.0;
  for(int ii = 0; ii < kind; ii++)
    Mtot += cfg[ii].Mtot;
  fprintf(fp, "Total mass of the system  Mtot is %e (= %e %s)\n" , Mtot, Mtot * mass2astro, mass_astro_unit_name);
  fprintf(fp, "Total number of particles Ntot is %zu (= 2^%u)\n", Ntot, ilog2((uint)Ntot));
  fprintf(fp, "Number of components      kind is %d\n" , kind);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Length of Plummer softening  is %e (= %e %s)\n", eps, eps * length2astro, length_astro_unit_name);
  fprintf(fp, "Snapshot interval            is %e (= %e %s)\n", snapshotInterval, snapshotInterval * time2astro, time_astro_unit_name);
  fprintf(fp, "Final time of the simulation is %e (= %e %s)\n",               ft,               ft * time2astro, time_astro_unit_name);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "#############################################################################\n\n");
  //-----------------------------------------------------------------------
  ulong ihead = 0;
  for(int ii = 0; ii < kind; ii++){
    //---------------------------------------------------------------------
    /* output settings for individual component */
    fprintf(fp, "#############################################################################\n");
    fprintf(fp, "#############################################################################\n");
    fprintf(fp, "%d-th component: ", ii);
    switch( cfg[ii].kind ){
    case   PLUMMER:      fprintf(fp,                        "Plummer model\n");      break;
    case      KING:      fprintf(fp,                           "King model\n");      break;
    case   BURKERT:      fprintf(fp,                        "Burkert model\n");      break;
    case HERNQUIST:      fprintf(fp,                      "Hernquist model\n");      break;
    case       NFW:      fprintf(fp,                            "NFW model\n");      break;
    case     MOORE:      fprintf(fp,                          "Moore model\n");      break;
    case   EINASTO:      fprintf(fp,                        "Einasto model\n");      break;
    case  APP_KING:      fprintf(fp,         "King model in empirical form\n");      break;
    case TWO_POWER:      fprintf(fp,   "Two-power (alpha-beta-gamma) model\n");      break;
    case TRI_POWER:      fprintf(fp,                    "Three-power model\n");      break;
    case APP_EVANS:      fprintf(fp,     "Approximated lowered Evans model\n");      break;
    case TABLE_RHO:      fprintf(fp,        "Density profile in table form\n");      break;
    case TABLE_SIG:      fprintf(fp, "Column density profile in table form\n");      break;
    case SPHSERSIC:      fprintf(fp,           "Sersic profile (spherical)\n");      break;
    case SIGTWOPOW:      fprintf(fp,        "Two-power model in projection\n");      break;
    case  EXP_DISK:      fprintf(fp,                     "Exponential disk\n");      break;
    case    SERSIC:      fprintf(fp,                "Sersic profile (disk)\n");      break;
    case CENTRALBH:      fprintf(fp,           "Central massive black hole\n");      break;
    default:
      __KILL__(stderr, "ERROR: undefined profile ``%d'' specified\n", cfg[ii].kind);
      break;
    }/* switch( cfg[ii].kind ){ */
    fprintf(fp, "#############################################################################\n");
    fprintf(fp, "Total number of particles Ntot is %zu (about 2^%u)\n", cfg[ii].num, ilog2((int)cfg[ii].num));
    fprintf(fp, "Range of idx for the component is [%zu, %zu]\n", ihead, ihead + cfg[ii].num - 1);    ihead += cfg[ii].num;
    fprintf(fp, "Mass of each N-body particle m is %e (= %e %s)\n", cfg[ii].Mtot / (double)cfg[ii].num, cfg[ii].Mtot / (double)cfg[ii].num * mass2astro, mass_astro_unit_name);
    fprintf(fp, "#############################################################################\n");
    if( cfg[ii].kind == TABLE_RHO )
      fprintf(fp, "Given density profile is written in %s\n", cfg[ii].table);
    if( cfg[ii].kind == TABLE_SIG )
      fprintf(fp, "Given column density profile is written in %s\n", cfg[ii].table);
    fprintf(fp, "Total mass of the component Mtot is %e (= %e %s)\n", cfg[ii].Mtot, cfg[ii].Mtot * mass2astro, mass_astro_unit_name);
    if( cfg[ii].kind != CENTRALBH )
      fprintf(fp, "Scale radius of the component rs is %e (= %e %s)\n", cfg[ii].rs, cfg[ii].rs * length2astro, length_astro_unit_name);
    if( cfg[ii].kind == KING ){
      fprintf(fp, "Dimensionless King parameter  W0 is %e\n", cfg[ii].king_W0);
      fprintf(fp, "Tidal radius of the component rt is %e (= %e %s)\n", cfg[ii].king_rt, cfg[ii].king_rt * length2astro, length_astro_unit_name);
      fprintf(fp, "Concentration paramter         c is %e\n", cfg[ii].king_c);
    }/* if( cfg[ii].kind == KING ){ */
    if( cfg[ii].kind == APP_KING ){
      fprintf(fp, "Tidal radius of the component rt is %e (= %e %s)\n", cfg[ii].king_rt, cfg[ii].king_rt * length2astro, length_astro_unit_name);
      fprintf(fp, "Concentration paramter         c is %e\n", log10(cfg[ii].king_rt / cfg[ii].rs));
    }/* if( cfg[ii].kind == APP_KING ){ */
    if( cfg[ii].kind == EINASTO )
      fprintf(fp, "Shape parameter            alpha is %e\n", cfg[ii].einasto_alpha);
    if( cfg[ii].kind == TWO_POWER ){
      fprintf(fp, "   Inner power-law index   alpha is %e\n", cfg[ii].twopower_alpha);
      fprintf(fp, "Internal power-law index    beta is %e\n", cfg[ii].twopower_beta);
      fprintf(fp, "   Outer power-law index   gamma is %e\n", cfg[ii].twopower_gamma);
    }/* if( cfg[ii].kind == TWO_POWER ){ */
    if( cfg[ii].kind == TRI_POWER ){
      fprintf(fp, "Outer transition radius         rout is %e (= %e %s)\n", cfg[ii].tripower_rout, cfg[ii].tripower_rout * length2astro, length_astro_unit_name);
      fprintf(fp, "   Innermost power-law index   alpha is %e\n", cfg[ii].twopower_alpha);
      fprintf(fp, "Transitional power-law index    beta is %e\n", cfg[ii].twopower_beta);
      fprintf(fp, "Intermediate power-law index   gamma is %e\n", cfg[ii].twopower_gamma);
      fprintf(fp, "Transitional power-law index   delta is %e\n", cfg[ii].tripower_delta);
      fprintf(fp, "   Outermost power-law index epsilon is %e\n", cfg[ii].tripower_epsilon);
    }/* if( cfg[ii].kind == TRI_POWER ){ */
    if( cfg[ii].kind == APP_EVANS ){
      fprintf(fp, "Second scale radius           rc is %e (= %e %s)\n", cfg[ii].alevans_rc, cfg[ii].alevans_rc * length2astro, length_astro_unit_name);
      fprintf(fp, "Inner power-law index      alpha is %e\n", cfg[ii].alevans_alpha);
      fprintf(fp, "Outer power-law index       beta is %e\n", cfg[ii].alevans_beta);
      fprintf(fp, "Exponential cutoff radius     rt is %e (= %e %s)\n", cfg[ii].alevans_rt, cfg[ii].alevans_rt * length2astro, length_astro_unit_name);
      fprintf(fp, "Exponential cutoff width      wt is %e (= %e %s)\n", cfg[ii].alevans_wt, cfg[ii].alevans_wt * length2astro, length_astro_unit_name);
      cfg[ii].rs = cfg[ii].alevans_rc;
      fprintf(fp, "Re-defined scale radius       rs is %e (= %e %s)\n", cfg[ii].rs, cfg[ii].rs * length2astro, length_astro_unit_name);
    }/* if( cfg[ii].kind == APP_EVANS ){ */
    if( cfg[ii].kind == SPHSERSIC ){
      fprintf(fp, "Sersic index                   n is %e\n", cfg[ii].n_sersic);
      fprintf(fp, "Dimensionless scale factor     b is %e\n", cfg[ii].b_sersic);
    }/* if( cfg[ii].kind == SPHSERSIC ){ */
    if( cfg[ii].kind == SIGTWOPOW ){
      fprintf(fp, "   Inner power-law index   alpha is %e\n", cfg[ii].twopower_alpha);
      fprintf(fp, "Internal power-law index    beta is %e\n", cfg[ii].twopower_beta);
      fprintf(fp, "   Outer power-law index   gamma is %e\n", cfg[ii].twopower_gamma);
    }/* if( cfg[ii].kind == SIGTWOPOW ){ */
    if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) ){
      if( cfg[ii].kind == SERSIC ){
	fprintf(fp, "Sersic index                   n is %e\n", cfg[ii].n_sersic);
	fprintf(fp, "Dimensionless scale factor     b is %e\n", cfg[ii].b_sersic);
      }/* if( cfg[ii].kind == SERSIC ){ */
      fprintf(fp, "Scale height of the component zd is %e (= %e %s)\n", cfg[ii].zd, cfg[ii].zd * length2astro, length_astro_unit_name);
      fprintf(fp, "Central surface density   Sigma0 is %e (= %e %s)\n", cfg[ii].Sigma0, cfg[ii].Sigma0 * col_density2astro, col_density_astro_unit_name);
      fprintf(fp, "Circular speed at scale radius   is %e (= %e %s)\n", cfg[ii].vcirc_Rd , cfg[ii].vcirc_Rd  * velocity2astro, velocity_astro_unit_name);
      fprintf(fp, "Maximum circular speed           is %e (= %e %s)\n", cfg[ii].vcirc_max, cfg[ii].vcirc_max * velocity2astro, velocity_astro_unit_name);
      fprintf(fp, "Circular speed is maximized      at %e (= %e %s)\n", cfg[ii].vcirc_Rd , cfg[ii].vcirc_Rd  *   length2astro,   length_astro_unit_name);
#ifdef  USE_ORIGINAL_VDISP_ESTIMATOR
      fprintf(fp, "Horizontal velocity dispersion   is %e of circular velocity or vertical velocity dispersion (maximum is used)\n", cfg[ii].vdisp_frac);
#else///USE_ORIGINAL_VDISP_ESTIMATOR
      fprintf(fp, "Horizontal velocity dispersion   is %e (= %e %s)\n", cfg[ii].vdispR0  , cfg[ii].vdispR0   * velocity2astro, velocity_astro_unit_name);
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
      fprintf(fp, "Vertical   velocity dispersion   is %e (= %e %s)\n", cfg[ii].vdispz0  , cfg[ii].vdispz0   * velocity2astro, velocity_astro_unit_name);
      fprintf(fp, "Toomre's Q-value at scale radius is %e\n", cfg[ii].toomre);
      fprintf(fp, "Minimum of Toomre's Q-value      is %e\n", cfg[ii].Qmin);
    }/* if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) ){ */
    if( cfg[ii].kind != CENTRALBH ){
      if( cfg[ii].cutoff ){
	fprintf(fp, "Cutoff radius of the component   is %e (= %e %s)\n", cfg[ii].rc      , cfg[ii].rc       * length2astro, length_astro_unit_name);
	fprintf(fp, "Cutoff  width of the component   is %e (= %e %s)\n", cfg[ii].rc_width, cfg[ii].rc_width * length2astro, length_astro_unit_name);
	fprintf(fp, "Cutoff radius over scale radius  is %e\n", cfg[ii].rc / cfg[ii].rs);
	fprintf(fp, "Cutoff  width over scale radius  is %e\n", cfg[ii].rc_width / cfg[ii].rs);
      }/* if( cfg[ii].cutoff ){ */
      fprintf(fp, "Cutoff energy of the component   is %e (= %e %s) (= %e G Mtot / rs)\n", cfg[ii].Ecut, cfg[ii].Ecut * senergy2astro, senergy_astro_unit_name, cfg[ii].Ecut / (newton * cfg[ii].Mtot / cfg[ii].rs));
    }/* if( cfg[ii].kind != CENTRALBH ){ */
    fprintf(fp, "#############################################################################\n");
    //---------------------------------------------------------------------
    /* estimate typical timescale */
    if( cfg[ii].kind != CENTRALBH ){
      //-------------------------------------------------------------------
      /* estimate enclosed mass within rs */
      double Ms;
      {
	int jj = 2;
	while( true ){
	  if( prf[0][jj].rad > cfg[ii].rs ){	    jj--;	    break;	  }
	  jj++;
	  if( jj == (NRADBIN + 1) )	    break;
	}/* while( true ){ */
	Ms = prf[0][jj].enc_tot;
      }
      //-------------------------------------------------------------------
      /* estimate dynamical time at scale length for each component */
      const double Ns = (double)Ntot * (Ms / Mtot);
      const double tff = M_PI_2 * cfg[ii].rs * sqrt(cfg[ii].rs / (2.0 * (double)newton * Ms));
      const double t2r = tff * Ns / (32.0 * log(cfg[ii].rs / (double)eps));
      double trot;
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) )
	trot = 2.0 * M_PI * cfg[ii].rs / cfg[ii].vcirc_Rd;
      fprintf(fp, "Total number of particles within the scale length is       %e\n", Ns);
      fprintf(fp, "Enclosed mass of all components within the scale length is %e (= %e %s)\n",  Ms,  Ms * mass2astro, mass_astro_unit_name);
      fprintf(fp, "Free-fall time at the scale length                      is %e (= %e %s)\n", tff, tff * time2astro, time_astro_unit_name);
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) )
	fprintf(fp, "Rotation time scale at the scale length                 is %e (= %e x tff = %e %s)\n", trot, trot / tff, trot * time2astro, time_astro_unit_name);
      fprintf(fp, "Two-body relaxation time at the scale length            is %e (= %e %s)\n", t2r, t2r * time2astro, time_astro_unit_name);
      fprintf(fp, "#############################################################################\n");
      fprintf(fp, "Snapshot interval in the unit of free-fall time           is %e\n", (double)snapshotInterval / tff);
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) )
	fprintf(fp, "Snapshot interval in the unit of rotation time scale      is %e\n", (double)snapshotInterval / trot);
      fprintf(fp, "Snapshot interval in the unit of two-body relaxation time is %e\n", (double)snapshotInterval / t2r);
      fprintf(fp, "#############################################################################\n");
      fprintf(fp, "Final time of the simulation in the unit of free-fall time           is %e\n", (double)ft / tff);
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) )
	fprintf(fp, "Final time of the simulation in the unit of rotation time scale      is %e\n", (double)ft / trot);
      fprintf(fp, "Final time of the simulation in the unit of two-body relaxation time is %e\n", (double)ft / t2r);
      fprintf(fp, "#############################################################################\n");
      fprintf(fp, "#############################################################################\n");
      //-------------------------------------------------------------------
    }/* if( cfg[ii].kind != CENTRALBH ){ */
    //---------------------------------------------------------------------
    fprintf(fp, "\n");
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < kind; ii++){ */
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* output fundamental profile of the particle distribution */
  //-----------------------------------------------------------------------
  real *tmp_rad;  tmp_rad = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_rad == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_rad");  }
  real *tmp_rho;  tmp_rho = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_rho");  }
  real *tmp_enc;  tmp_enc = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_enc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_enc");  }
  real *tmp_psi;  tmp_psi = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_psi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_psi");  }
  real *tmp_tff;  tmp_tff = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_tff == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_tff");  }
  real *tmp_t2r;  tmp_t2r = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_t2r == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_t2r");  }
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "profile");
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#ifdef  DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  uint szip_pixels_per_block = 8;
  hsize_t cdims = 128 * szip_pixels_per_block;
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < kind; kk++){
    //---------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
    char grp[16];    sprintf(grp, "data%d", kk);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else///USE_HDF5_FORMAT
    sprintf(filename, "%s/%s.profile.%d.dat", DATAFOLDER, file, kk);
    fp = fopen(filename, "wb");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }
#endif//USE_HDF5_FORMAT
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      tmp_rad[ii] = (real)(prf[kk][ii].rad *  length2astro);
      tmp_rho[ii] = (real)(prf[kk][ii].rho * density2astro);
      tmp_enc[ii] = (real)(prf[kk][ii].enc *    mass2astro);
      tmp_psi[ii] = (real)(prf[kk][ii].psi * senergy2astro);
    }
    //---------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
    hsize_t dims = NRADBIN;
    dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
    /* write radius */
    dataset = H5Dcreate(group, "rad", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rad));
    chkHDF5err(H5Dclose(dataset));
    /* write density */
    dataset = H5Dcreate(group, "rho", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rho));
    chkHDF5err(H5Dclose(dataset));
    /* write enclosed mass */
    dataset = H5Dcreate(group, "enc", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_enc));
    chkHDF5err(H5Dclose(dataset));
    /* write potential */
    dataset = H5Dcreate(group, "Psi", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_psi));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------
    /* write attribute data */
    //---------------------------------------------------------------------
    /* create the data space for the attribute */
    hsize_t attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    //---------------------------------------------------------------------
    /* write # of arrays */
    int nradbin = NRADBIN;
    attribute = H5Acreate(group, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nradbin));
    chkHDF5err(H5Aclose(attribute));
    /* write scale radius */
    attribute = H5Acreate(group, "rs", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg[kk].rs));
    chkHDF5err(H5Aclose(attribute));
    /* write total mass */
    attribute = H5Acreate(group, "Mtot", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg[kk].Mtot));
    chkHDF5err(H5Aclose(attribute));
    /* profile ID */
    attribute = H5Acreate(group, "profile ID", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cfg[kk].kind));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
    //---------------------------------------------------------------------
#else///USE_HDF5_FORMAT
    int nradbin = NRADBIN;
    bool success = true;
    success &= (fwrite(&nradbin, sizeof( int),       1, fp) ==       1);
    success &= (fwrite( tmp_rad, sizeof(real), NRADBIN, fp) == NRADBIN);
    success &= (fwrite( tmp_rho, sizeof(real), NRADBIN, fp) == NRADBIN);
    success &= (fwrite( tmp_enc, sizeof(real), NRADBIN, fp) == NRADBIN);
    success &= (fwrite( tmp_psi, sizeof(real), NRADBIN, fp) == NRADBIN);
    if( !success ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
    fclose(fp);
#endif//USE_HDF5_FORMAT
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < kind; kk++){ */
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  /* evaluate typical timescale */
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad = prf[0][ii].rad;
    double enc = 0.0;
    for(int jj = 0; jj < kind; jj++)
      enc += prf[jj][ii].enc;
    const double tff = M_PI_2 * rad * sqrt(rad / (2.0 * (double)newton * enc));
    double t2r = tff * ((double)Ntot * (enc / Mtot)) / (32.0 * log(enc / (double)eps));
    if( t2r < 0.0 )
      t2r = 0.0;
    tmp_rad[ii] = (real)(rad * length2astro);
    tmp_tff[ii] = (real)(tff *   time2astro);
    tmp_t2r[ii] = (real)(t2r *   time2astro);
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */
  //-----------------------------------------------------------------------
  /* write typical timescale */
  hid_t timeScale = H5Gcreate(target, "time_scale", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t dims = NRADBIN;
  dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  chkHDF5err(H5Pset_chunk(property, 1, &cdims));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
  /* write radius */
  dataset = H5Dcreate(timeScale, "rad", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rad));
  chkHDF5err(H5Dclose(dataset));
  /* write free-fall time */
  dataset = H5Dcreate(timeScale, "free_fall_time", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_tff));
  chkHDF5err(H5Dclose(dataset));
  /* write relaxation time */
  dataset = H5Dcreate(timeScale, "relaxation_time", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_t2r));
  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  int nradbin = NRADBIN;
  hid_t attribute = H5Acreate(timeScale, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nradbin));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(timeScale));
  //-----------------------------------------------------------------------
  /* write attribute data */
  //-----------------------------------------------------------------------
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
  /* write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &kind));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
  attribute = H5Acreate(target,  "length_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,  length_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, density_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,    "mass_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,    mass_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "senergy_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, senergy_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,    "time_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,    time_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* close the file */
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
#endif//USE_HDF5_FORMAT
  free(tmp_rad);
  free(tmp_rho);
  free(tmp_enc);
  free(tmp_psi);
  free(tmp_tff);
  free(tmp_t2r);
  //-----------------------------------------------------------------------
#ifdef  OUTPUT_ASCII_PROFILE
  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.profile.txt", DATAFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
  //-----------------------------------------------------------------------
  fprintf(fp, "#r\trho(r)\tM(r)\tPhi(r)\n");
  fprintf(fp, "#\tgenerated on %s", date);
  fprintf(fp, "#number of conponents is %d\n", kind);
  fprintf(fp, "#format is total");
  for(int kk = 0; kk < kind; kk++)
    fprintf(fp, ", %d-th", kk);
  fprintf(fp, "\n");
  const int nskip = NRADBIN / N_PRINT_LINES_ASCII;
  for(int ii = 2; ii < 2 + NRADBIN; ii += nskip){
    fprintf(fp, "%e", prf[0][ii].rad);
    fprintf(fp, "\t%e",  prf[0][ii].rho_tot);    for(int kk = 0; kk < kind; kk++)      fprintf(fp, "\t%e",  prf[kk][ii].rho);
    fprintf(fp, "\t%e",  prf[0][ii].enc_tot);    for(int kk = 0; kk < kind; kk++)      fprintf(fp, "\t%e",  prf[kk][ii].enc);
    fprintf(fp, "\t%e", -prf[0][ii].psi_tot);    for(int kk = 0; kk < kind; kk++)      fprintf(fp, "\t%e", -prf[kk][ii].psi);
    fprintf(fp, "\n");
  }
  fclose(fp);
  //-----------------------------------------------------------------------
#endif//OUTPUT_ASCII_PROFILE
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* output distribution function of the particle distribution */
  //-----------------------------------------------------------------------
  real *tmp_ene;  tmp_ene = (real *)malloc(NENEBIN * sizeof(real));  if( tmp_ene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_ene");  }
  real *tmp_val;  tmp_val = (real *)malloc(NENEBIN * sizeof(real));  if( tmp_val == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_val");  }
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "df");
  target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < skind; kk++){
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NENEBIN; ii++){
      tmp_ene[ii] = (real)((double)df[kk][ii].ene * senergy2astro);
      tmp_val[ii] =                df[kk][ii].val;
    }
    //---------------------------------------------------------------------
    char grp[16];    sprintf(grp,  "series%d", kk);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //---------------------------------------------------------------------
    hsize_t dims = NENEBIN;
    dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
    /* write energy */
    dataset = H5Dcreate(group, "energy", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_ene));
    chkHDF5err(H5Dclose(dataset));
    /* write energy */
    dataset = H5Dcreate(group, "DF", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_val));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------
    /* write attribute data */
    //---------------------------------------------------------------------
    /* create the data space for the attribute */
    hsize_t attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    //---------------------------------------------------------------------
    /* write # of arrays */
    int nenebin = NENEBIN;
    attribute = H5Acreate(group, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nenebin));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < skind; kk++){ */
  //-----------------------------------------------------------------------
  /* write attribute data */
  //-----------------------------------------------------------------------
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
  /* write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skind));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  attribute = H5Acreate(target, "senergy_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, senergy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));
  //-----------------------------------------------------------------------
  /* write flag about DOUBLE_PRECISION */
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* close the file */
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
#else///USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < skind; kk++){
    //---------------------------------------------------------------------
    sprintf(filename, "%s/%s.df.%d.dat", DATAFOLDER, file, kk);
    fp = fopen(filename, "wb");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NENEBIN; ii++){
      tmp_ene[ii] = df[kk][ii].ene;
      tmp_val[ii] = df[kk][ii].val;
    }/* for(int ii = 0; ii < NENEBIN; ii++){ */
    //---------------------------------------------------------------------
    int nenebin = NENEBIN;
    bool success = true;
    success &= (fwrite(&nenebin, sizeof( int),       1, fp) ==       1);
    success &= (fwrite( tmp_ene, sizeof(real), NENEBIN, fp) == NENEBIN);
    success &= (fwrite( tmp_val, sizeof(real), NENEBIN, fp) == NENEBIN);
    if( !success ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
    fclose(fp);
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < skind; kk++){ */
  //-----------------------------------------------------------------------
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  free(tmp_ene);
  free(tmp_val);
  //-----------------------------------------------------------------------
#ifdef  OUTPUT_ASCII_PROFILE
  //-----------------------------------------------------------------------
  if( skind > 0 ){
    //---------------------------------------------------------------------
    sprintf(filename, "%s/%s.df.txt", DATAFOLDER, file);
    fp = fopen(filename, "w");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }
    //---------------------------------------------------------------------
    fprintf(fp, "#relE\tDF(relE)\n");
    fprintf(fp, "#\tgenerated on %s", date);
    fprintf(fp, "#number of spherical conponent(s) is %d\n", skind);
    fprintf(fp, "#format is energy\t%d-th", 0);
    for(int kk = 1; kk < skind; kk++)
      fprintf(fp, "\t%d-th", kk);
    fprintf(fp, "\n");
    const int nskip_ene = NENEBIN / N_PRINT_LINES_ASCII;
    for(int ii = 0; ii < NENEBIN; ii += nskip_ene){
      fprintf(fp, "%e", df[0][ii].ene);
      for(int kk = 0; kk < skind; kk++)
	fprintf(fp, "\t%e", df[kk][ii].val);
      fprintf(fp, "\n");
    }/* for(int ii = 0; ii < NENEBIN; ii += nskip_ene){ */
    fclose(fp);
    //---------------------------------------------------------------------
  }/* if( skind > 0 ){ */
  //-----------------------------------------------------------------------
#endif//OUTPUT_ASCII_PROFILE
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void isotropicDistribution(const real rad, real *vecx, real *vecy, real *vecz)
{
  //-----------------------------------------------------------------------
  const real proj = RANDVAL;
  *vecz = rad * proj;
  real Rproj = rad * SQRT(UNITY - proj * proj);
  //-----------------------------------------------------------------------
  real theta = TWO * (real)M_PI * UNIRAND;
  *vecx = Rproj * COS(theta);
  *vecy = Rproj * SIN(theta);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void shiftCenter(ulong num, nbody_particle *body)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate center-of-mass and bulk-motion */
  //-----------------------------------------------------------------------
  double com[3] = {0.0, 0.0, 0.0};
  double vel[3] = {0.0, 0.0, 0.0};
  double Mtot = 0.0;
  //-----------------------------------------------------------------------
  for(ulong ii = 0; ii < num; ii++){
    com[0] += (double)(body[ii].m * body[ii].x);    vel[0] += (double)(body[ii].m * body[ii].vx);
    com[1] += (double)(body[ii].m * body[ii].y);    vel[1] += (double)(body[ii].m * body[ii].vy);
    com[2] += (double)(body[ii].m * body[ii].z);    vel[2] += (double)(body[ii].m * body[ii].vz);
    Mtot   += (double)body[ii].m;
  }
  //-----------------------------------------------------------------------
  double Minv = 1.0 / Mtot;
  com[0] *= Minv;  vel[0] *= Minv;
  com[1] *= Minv;  vel[1] *= Minv;
  com[2] *= Minv;  vel[2] *= Minv;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* shift the coordinate system to the center-of-mass rest frame */
  //-----------------------------------------------------------------------
  const real rcom[3] = {(real)com[0], (real)com[1], (real)com[2]};
  const real rvel[3] = {(real)vel[0], (real)vel[1], (real)vel[2]};
#pragma omp parallel for
  for(ulong ii = 0; ii < num; ii++){
    body[ii].x -= rcom[0];    body[ii].vx -= rvel[0];
    body[ii].y -= rcom[1];    body[ii].vy -= rvel[1];
    body[ii].z -= rcom[2];    body[ii].vz -= rvel[2];
#ifdef  BLOCK_TIME_STEP
    body[ii].t0 = body[ii].t1 = 0.0;
    body[ii].dt = ZERO;
#endif//BLOCK_TIME_STEP
  }/* for(ulong ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
#if 0
  printf("#position shift = (%e, %e, %e)\n", com[0], com[1], com[2]);
  printf("#velocity shift = (%e, %e, %e)\n", vel[0], vel[1], vel[2]);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  RESET_ROTATION_AXIS
  //-----------------------------------------------------------------------
  /* calculate angular momentum vector */
  //-----------------------------------------------------------------------
  double amom[3] = {0.0, 0.0, 0.0};
  for(ulong ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    const double rx = (double)              body[ii]. x;
    const double ry = (double)              body[ii]. y;
    const double rz = (double)              body[ii]. z;
    const double px = (double)(body[ii].m * body[ii].vx);
    const double py = (double)(body[ii].m * body[ii].vy);
    const double pz = (double)(body[ii].m * body[ii].vz);
    //---------------------------------------------------------------------
    amom[0] += ry * pz - rz * py;
    amom[1] += rz * px - rx * pz;
    amom[2] += rx * py - ry * px;
    //---------------------------------------------------------------------
  }/* for(ulong ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
  /* rotate galaxy (if necessary) */
  //-----------------------------------------------------------------------
  const double L2 = amom[0] * amom[0] + amom[1] * amom[1] + amom[2] * amom[2];
  if( L2 > 1.0e-10 ){
    //---------------------------------------------------------------------
    real ini[3] = {(real)amom[0], (real)amom[1], (real)amom[2]};
    real fin[3] = {ZERO, ZERO, UNITY};
    //---------------------------------------------------------------------
    real rot[3][3], inv[3][3];
    initRotationMatrices(ini, fin, rot, inv);
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(ulong ii = 0; ii < num; ii++){
      //-------------------------------------------------------------------
      real bfr[3], aft[3];
      //-------------------------------------------------------------------
      /* rotate position */
      bfr[0] = body[ii].x;
      bfr[1] = body[ii].y;
      bfr[2] = body[ii].z;
      rotateVector(bfr, rot, aft);
      body[ii].x = aft[0];
      body[ii].y = aft[1];
      body[ii].z = aft[2];
      //-------------------------------------------------------------------
      /* rotate velocity */
      bfr[0] = body[ii].vx;
      bfr[1] = body[ii].vy;
      bfr[2] = body[ii].vz;
      rotateVector(bfr, rot, aft);
      body[ii].vx = aft[0];
      body[ii].vy = aft[1];
      body[ii].vz = aft[2];
      //-------------------------------------------------------------------
    }/* for(ulong ii = 0; ii < num; ii++){ */
    //---------------------------------------------------------------------
  }/* if( L2 > 1.0e-10 ){ */
  //-----------------------------------------------------------------------
#endif//RESET_ROTATION_AXIS
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline double getDF(const double ene, dist_func *df, const double Emin, const double invEbin)
{
  //-----------------------------------------------------------------------
#if 1
#if 0
  fprintf(stderr, "ene = %e, Emin = %e, invEbin = %e\n", ene, Emin, invEbin);
  fflush(stderr);
#endif
  const int ll = (int)((ene - Emin) * invEbin);
  return (df[ll].val + (df[ll + 1].val - df[ll].val) * (ene - df[ll].ene) / (df[ll + 1].ene - df[ll].ene));
#else
  int ll = (int)floor((ene - Emin) * invEbin);
  if( ll <           0 )    ll =           0;
  if( ll > NENEBIN - 2 )    ll = NENEBIN - 2;
  const int rr = ll + 1;
  const double alpha = (ene - df[ll].ene) * invEbin;
  return ((1.0 - alpha) * df[ll].val + alpha * df[rr].val);
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
double distributeSpheroidParticles(ulong *Nuse, nbody_particle *body, const real mass, profile_cfg cfg, profile *prf, dist_func *df)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const ulong num = cfg.num;
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# start distributing spherical component particles (%zu bodies: [%zu:%zu])\n", num, *Nuse, (*Nuse) + num - 1);
  fflush(stdout);
  const ulong nunit = (ulong)ceilf(0.1f * (float)num);
  ulong stage = 1;
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------
  const double    Emin = df[          0].ene;
  const double    Emax = df[NENEBIN - 1].ene;
  const double invEbin = (double)(NENEBIN - 1) / (Emax - Emin);
  //-----------------------------------------------------------------------
  double fmax = 0.0;
  int   iout = 0;
  for(int ii = 0; ii < NRADBIN; ii++){
    double floc = prf[ii].rad * prf[ii].rad * prf[ii].rho;
    if(                            floc > fmax ){      fmax = floc;    }
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */
  const double Ecut = (iout != 0) ? prf[iout].psi_tot : Emin;
#if 0
  fprintf(stderr, "Emin = %e, Emax = %e, fmax = %e, iout = %d, Ecut = %e\n", Emin, Emax, fmax, iout, Ecut);
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const double Mmax = cfg.Mtot;
  const double Mmin = prf[0].enc;
  for(ulong ii = *Nuse; ii < *Nuse + num; ii++){
    //---------------------------------------------------------------------
    /* set spatial distribution by table search */
    //---------------------------------------------------------------------
    const double tmp = Mmin + (Mmax - Mmin) * UNIRAND_DBL;
    int ll = 0;
    int rr = NRADBIN - 1;
    while( true ){
      //-------------------------------------------------------------------
      uint cc = (ll + rr) >> 1;
      //-------------------------------------------------------------------
      if( (prf[cc].enc - tmp) * (prf[ll].enc - tmp) <= 0.0 )	rr = cc;
      else	                                                ll = cc;
      //-------------------------------------------------------------------
      if( (ll + 1) == rr )	break;
      //-------------------------------------------------------------------
    }/* while( true ){ */
    //---------------------------------------------------------------------
    const double alpha = (tmp - prf[ll].enc) / (prf[rr].enc - prf[ll].enc);
    const double rad = (1.0 - alpha) * prf[ll].rad + alpha * prf[rr].rad;
    isotropicDistribution((real)rad, &(body[ii].x), &(body[ii].y), &(body[ii].z));
    __NOTE__("position of %zu-th particle determined: rad = %e\n", ii, rad);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* determine velocity distribution by rejection method */
    //---------------------------------------------------------------------
#if 1
    const double psi = (1.0 - alpha) * prf[ll].psi_tot + alpha * prf[rr].psi_tot;
#else
    const double psi = prf[ll].psi_tot + (prf[rr].psi_tot - prf[ll].psi_tot) * (tmp - prf[ll].enc) / (prf[rr].enc - prf[ll].enc);
#endif
    /* const double vesc = sqrt(2.0 * psi); */
    const double vesc = sqrt(2.0 * (psi - Ecut));
    const double v2Fmax = vesc * vesc * getDF(psi, df, Emin, invEbin);
    double vel;
    while( true ){
      //-------------------------------------------------------------------
      vel = vesc * UNIRAND_DBL;
      //-------------------------------------------------------------------
      const double ene = psi - 0.5 * vel * vel;
      const double val = vel * vel * getDF(ene, df, Emin, invEbin);
      const double try = v2Fmax * UNIRAND_DBL;
      if( val > try )	break;
      //-------------------------------------------------------------------
    }/* while( true ){ */
    //---------------------------------------------------------------------
    isotropicDistribution((real)vel, &(body[ii].vx), &(body[ii].vy), &(body[ii].vz));
    __NOTE__("velocity of %zu-th particle determined\n", ii);
    //---------------------------------------------------------------------
    body[ii].ax = body[ii].ay = body[ii].az = ZERO;
    body[ii].m   = mass;
    body[ii].pot = ZERO;
    //---------------------------------------------------------------------
    body[ii].idx = ii;
    //---------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
    if( (ii - (*Nuse)) == (stage * nunit) ){
      fprintf(stdout, "# ~%zu%% completed\n", stage * 10);
      fflush(stdout);
      stage++;
    }/* if( (ii - (*Nuse)) == (stage * nunit) ){ */
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------
  }/* for(ulong ii = *Nuse; ii < *Nuse + num; ii++){ */
  //-----------------------------------------------------------------------
  *Nuse += num;
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# finish distributing spherical component particles (%zu bodies)\n#\n#\n", num);
  fflush(stdout);
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (Ecut);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* real * restrict rhoFrac is an input array */
/* real * restrict _????? is output array */
//-------------------------------------------------------------------------
static void evaluateDiskProperties
(char *file, disk_data *disk_info, const int diskID, real * restrict rhoFrac,
 real * restrict _vcirc, real * restrict _sigmap, real * restrict _sigmaR, real * restrict _kappa, real * restrict _Omega, real * restrict _toomre, real * restrict _lambda)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifndef USE_ORIGINAL_VDISP_ESTIMATOR
  const double invRd = 1.0 / disk_info[diskID].cfg->rs;
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
  disk_info[diskID].cfg->vcirc_max   = -1.0;
  disk_info[diskID].cfg->vcirc_max_R = -1.0;
  disk_info[diskID].cfg->Qmin = DBL_MAX;
  bool passed = false;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NDISKBIN_HOR; ii++){
    //---------------------------------------------------------------------
    /* evaluate epicyclic quantities and circular speed */
#ifndef USE_POTENTIAL_SCALING_SCHEME
    const double Omega = sqrt(disk_info[diskID]. dPhidR [INDEX2D(NDISKBIN_RAD, NDISKBIN_VER, ii, 0)] / disk_info[diskID].hor[ii]);
    const double kappa = sqrt(disk_info[diskID].d2PhidR2[INDEX2D(NDISKBIN_RAD, NDISKBIN_VER, ii, 0)] + 3.0 * Omega * Omega);
#else///USE_POTENTIAL_SCALING_SCHEME
    const double Omega = sqrt(disk_info[diskID]. dPhidR [ii] / disk_info[diskID].hor[ii]);
    const double kappa = sqrt(disk_info[diskID].d2PhidR2[ii] + 3.0 * Omega * Omega);
#endif//USE_POTENTIAL_SCALING_SCHEME
    const double vcirc = disk_info[diskID].hor[ii] * Omega;
    //---------------------------------------------------------------------
    /* evaluate Toomre's Q-value */
#ifdef  USE_ORIGINAL_VDISP_ESTIMATOR
    const double sigmap = DISK_PERP_VDISP(disk_info[diskID].sigmaz[ii], vcirc, disk_info[diskID].cfg->vdisp_frac);
    /* const double gam2inv = 0.25 * (3.0 + disk_info[diskID].d2PhidR2[ii] / (1.0e-100 + Omega * Omega)); */
    /* const double sigmaR = sigmap / sqrt(gam2inv); */
    const double sigmaR = sigmap * 2.0 * Omega / (DBL_MIN + kappa);
#else///USE_ORIGINAL_VDISP_ESTIMATOR
    const double sigmaR = DISK_RADIAL_VDISP(disk_info[diskID].cfg->vdispR0, disk_info[diskID].hor[ii], invRd);
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
    const double toomre = sigmaR * kappa / (3.36 * newton * disk_info[diskID].Sigma[ii]);
    const double lambda = 4.0 * M_PI * M_PI * newton * disk_info[diskID].Sigma[ii] / (DBL_MIN + kappa * kappa);
    //---------------------------------------------------------------------
    /* find the maximum circular speed */
    if( vcirc > disk_info[diskID].cfg->vcirc_max ){
      disk_info[diskID].cfg->vcirc_max   = vcirc;
      disk_info[diskID].cfg->vcirc_max_R = disk_info[diskID].hor[ii];
    }/* if( vcirc > disk_info[diskID].cfg->vcirc_max ){ */
    //---------------------------------------------------------------------
    /* find the minimum Toomre's Q-value */
    if( (disk_info[diskID].Sigma[ii] > 1.0e-4 * disk_info[diskID].Sigma[0]) && (toomre > DBL_EPSILON) && (toomre < disk_info[diskID].cfg->Qmin) )
      disk_info[diskID].cfg->Qmin = toomre;
    //---------------------------------------------------------------------
    /* find the circular speed and Toomre's Q-value at the scale radius */
    if( !passed ){
      disk_info[diskID].cfg->vcirc_Rd = vcirc;
      disk_info[diskID].cfg->toomre   = toomre;
      if( disk_info[diskID].hor[ii] > disk_info[diskID].cfg->rs )
	passed = true;
    }/* if( !passed ){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* memorize calculated values */
    //---------------------------------------------------------------------
    _vcirc [ii] = (real)(vcirc  * velocity2astro);
    _sigmap[ii] = (real)(sigmap * velocity2astro);
    _sigmaR[ii] = (real)(sigmaR * velocity2astro);
    _kappa [ii] = (real)(kappa  /     time2astro);
    _Omega [ii] = (real)(Omega  /     time2astro);
    _lambda[ii] = (real)(lambda *   length2astro);
    _toomre[ii] = (real) toomre;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output profile data in ASCII */
  //-----------------------------------------------------------------------
  FILE *fp;
  char filename[256];
  sprintf(filename, "%s/%s.diskhor.%d.txt", DATAFOLDER, file, diskID);
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);
  }
  fprintf(fp, "#R (%s)\tSigma(R) (%s)\tfrac(R)\tv_c(R) (%s)\tsigma_p(R) (%s)\tsigma_R(R) (%s)\tsigma_z(R) (%s)\tkappa(R) (1 / %s)\tOmega(R) (1/%s)\tToomre-Q(R)\tlambda_crit(R) (%s)\n",
	  length_astro_unit_name, col_density_astro_unit_name, velocity_astro_unit_name, velocity_astro_unit_name, velocity_astro_unit_name, velocity_astro_unit_name, time_astro_unit_name, time_astro_unit_name, length_astro_unit_name);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NDISKBIN_HOR; ii++)
    fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", disk_info[diskID].hor[ii] * length2astro, disk_info[diskID].Sigma[ii] * col_density2astro, rhoFrac[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 0)], _vcirc[ii], _sigmap[ii], _sigmaR[ii], disk_info[diskID].sigmaz[ii] * velocity2astro, _kappa[ii], _Omega[ii], _toomre[ii], _lambda[ii]);
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static void writeDiskData(char *file, const int ndisk, disk_data *disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  real *tmp_hor;  tmp_hor = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
  real *tmp_ver;  tmp_ver = (real *)malloc(               NDISKBIN_VER * sizeof(real));
  real *tmp_sig;  tmp_sig = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
  real *tmp_Sig;  tmp_Sig = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
  real *_vcirc ;  _vcirc  = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
  real *_sigmap;  _sigmap = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
  real *_sigmaR;  _sigmaR = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
  real *_kappa ;  _kappa  = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
  real *_Omega ;  _Omega  = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
  real *_lambda;  _lambda = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
  real *_toomre;  _toomre = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  real *tmp__zd;  tmp__zd = (real *)malloc(NDISKBIN_HOR                * sizeof(real));
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  real *tmp_rho;  tmp_rho = (real *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(real));
  real *tmp_Phi;  tmp_Phi = (real *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(real));
  real *rhoFrac;  rhoFrac = (real *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(real));
  if( tmp_hor == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_hor");  }
  if( tmp_ver == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_ver");  }
  if( tmp_sig == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_sig");  }
  if( tmp_Sig == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_Sig");  }
  if( _vcirc  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _vcirc" );  }
  if( _sigmap == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _sigmap");  }
  if( _sigmaR == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _sigmaR");  }
  if( _kappa  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _kappa" );  }
  if( _Omega  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _Omega" );  }
  if( _lambda == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _lambda");  }
  if( _toomre == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _toomre");  }
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  if( tmp__zd == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp__zd");  }
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  if( tmp_rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_rho");  }
  if( tmp_Phi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_Phi");  }
  if( rhoFrac == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rhoFrac");  }
  //-----------------------------------------------------------------------
  if( ndisk > 1 )
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR * NDISKBIN_VER; ii++){
      disk[0].rhoTot[ii] = (*disk[0].rho)[ii];
      for(int jj = 1; jj < ndisk; jj++)
	disk[0].rhoTot[ii] += (*disk[jj].rho)[ii];
      disk[0].rhoTot[ii] = 1.0 / (DBL_MIN + disk[0].rhoTot[ii]);
    }/* for(int ii = 0; ii < NDISKBIN_HOR * NDISKBIN_VER; ii++){ */
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  char filename[128];
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "disk");
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#ifdef  DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  uint szip_pixels_per_block = 8;
  /* hsize_t cdims[2] = {128, 128 * szip_pixels_per_block}; */
  hsize_t cdims[2] = {128, 128};
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    /* output in HDF5 format */
    //---------------------------------------------------------------------
    char grp[16];    sprintf(grp, "data%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write two-dimensional data */
    //---------------------------------------------------------------------
    if( ndisk > 1 )
#pragma omp parallel for
      for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++)
	rhoFrac[jj] = (real)((*disk[ii].rho)[jj] * disk[ii].rhoTot[jj]);
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++){
      tmp_rho[jj] = (real)((*disk[ii].rho)[jj] * density2astro);
      tmp_Phi[jj] = (real)(  disk[ii].pot [jj] * senergy2astro);
    }/* for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++){ */
    //---------------------------------------------------------------------
    evaluateDiskProperties(file, disk, ii, rhoFrac, _vcirc, _sigmap, _sigmaR, _kappa, _Omega, _toomre, _lambda);
    //---------------------------------------------------------------------
    hsize_t dims[2] = {NDISKBIN_HOR, NDISKBIN_VER};
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 2, cdims));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
    /* write density */
    dataset = H5Dcreate(group, "rho", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rho));
    chkHDF5err(H5Dclose(dataset));
    /* write potential */
    dataset = H5Dcreate(group, "Phi", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_Phi));
    chkHDF5err(H5Dclose(dataset));
    /* write fraction */
    if( ndisk > 1 ){
      dataset = H5Dcreate(group, "fraction", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rhoFrac));
      chkHDF5err(H5Dclose(dataset));
    }/* if( ndisk > 1 ){ */
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write one-dimensional data */
    //---------------------------------------------------------------------
    for(int jj = 0; jj < NDISKBIN_HOR; jj++){
      tmp_hor[jj] = (real)(disk[ii].hor   [jj] *      length2astro);
      tmp_sig[jj] = (real)(disk[ii].sigmaz[jj] *    velocity2astro);
      tmp_Sig[jj] = (real)(disk[ii].Sigma [jj] * col_density2astro);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      tmp__zd[jj] = (real)(disk[ii].zd    [jj] *      length2astro);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    }/* for(int jj = 0; jj < NDISKBIN_HOR; jj++){ */
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
      tmp_ver[jj] = (real)(disk[ii].ver   [jj] *      length2astro);
    //---------------------------------------------------------------------
    /* horizontal variables */
    dataspace = H5Screate_simple(1, &dims[0], NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims[1]));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
    /* write horizontal position */
    dataset = H5Dcreate(group, "radius", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_hor));
    chkHDF5err(H5Dclose(dataset));
    /* write velocity dispersion in z-direction */
    dataset = H5Dcreate(group, "sigmaz", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_sig));
    chkHDF5err(H5Dclose(dataset));
    /* write column density */
    dataset = H5Dcreate(group, "Sigma", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_Sig));
    chkHDF5err(H5Dclose(dataset));
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
    /* write scale height */
    dataset = H5Dcreate(group, "zd", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp__zd));
    chkHDF5err(H5Dclose(dataset));
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    /* write velocity dispersion in R-direction */
    dataset = H5Dcreate(group, "sigmaR", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _sigmaR));
    chkHDF5err(H5Dclose(dataset));
    /* write velocity dispersion in tangential direction */
    dataset = H5Dcreate(group, "sigmap", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _sigmap));
    chkHDF5err(H5Dclose(dataset));
    /* write circular velocity */
    dataset = H5Dcreate(group, "vcirc", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _vcirc));
    chkHDF5err(H5Dclose(dataset));
    /* write kappa */
    dataset = H5Dcreate(group, "kappa", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _kappa));
    chkHDF5err(H5Dclose(dataset));
    /* write Omega */
    dataset = H5Dcreate(group, "Omega", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _Omega));
    chkHDF5err(H5Dclose(dataset));
    /* write lambda_critical */
    dataset = H5Dcreate(group, "lambda", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _lambda));
    chkHDF5err(H5Dclose(dataset));
    /* write Toomre's Q-value */
    dataset = H5Dcreate(group, "Toomre's Q", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _toomre));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------
    /* vertical variables */
    dataspace = H5Screate_simple(1, &dims[1], NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims[1]));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
    /* write vertical position */
    dataset = H5Dcreate(group, "height", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_ver));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write attribute data */
    //---------------------------------------------------------------------
    /* create the data space for the attribute */
    hsize_t attr_dims = 2;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    //---------------------------------------------------------------------
    /* write # of arrays */
    int narray[2] = {NDISKBIN_HOR, NDISKBIN_VER};
    attribute = H5Acreate(group, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, narray));
    chkHDF5err(H5Aclose(attribute));
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    narray[0] = 1;
    /* write scale radius */
    attribute = H5Acreate(group, "Rs", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &disk[ii].cfg->rs));
    chkHDF5err(H5Aclose(attribute));
    /* write scale height */
    attribute = H5Acreate(group, "zd", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &disk[ii].cfg->zd));
    chkHDF5err(H5Aclose(attribute));
    /* write total mass */
    attribute = H5Acreate(group, "Mtot", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &disk[ii].cfg->Mtot));
    chkHDF5err(H5Aclose(attribute));
    /* profile ID */
    attribute = H5Acreate(group, "profile ID", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &disk[ii].cfg->kind));
    chkHDF5err(H5Aclose(attribute));
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------
    chkHDF5err(H5Gclose(group));
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  //-----------------------------------------------------------------------
  /* write attribute data */
  //-----------------------------------------------------------------------
  /* create the data space for the attribute */
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  //-----------------------------------------------------------------------
  /* write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &ndisk));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
  attribute = H5Acreate(target,      "length_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,      length_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,     "density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,     density_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,     "senergy_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,     senergy_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,    "velocity_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,    velocity_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "col_density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, col_density_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,        "time_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,        time_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));
  //-----------------------------------------------------------------------
  /* write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* close the file */
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
#else///USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    FILE *fp;
    char filename[256];
    sprintf(filename, "%s/%s.diskdat.%d.dat", DATAFOLDER, file, ii);
    fp = fopen(filename, "wb");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }
    //---------------------------------------------------------------------
    bool success = true;
    const int nhorbin = NDISKBIN_HOR;
    const int nverbin = NDISKBIN_VER;
    /* const int nsphbin = NDISKBIN_RAD; */
    success &= (fwrite(&nhorbin, sizeof(int), 1, fp) == 1);
    success &= (fwrite(&nverbin, sizeof(int), 1, fp) == 1);
    /* success &= (fwrite(&nsphbin, sizeof(int), 1, fp) == 1); */
    //---------------------------------------------------------------------
    if( ndisk > 1 )
#pragma omp parallel for
      for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++)
	rhoFrac[jj] = (real)((*disk[ii].rho)[jj] * disk[ii].rhoTot[jj]);
    evaluateDiskProperties(file, disk, ii, rhoFrac, _vcirc, _sigmap, _sigmaR, _kappa, _Omega, _toomre, _lambda);
    //---------------------------------------------------------------------
    for(int jj = 0; jj < nhorbin; jj++){
      tmp_hor[jj] = (real)(disk[ii].hor   [jj] *      length2astro);
      tmp_sig[jj] = (real)(disk[ii].sigmaz[jj] *    velocity2astro);
      tmp_Sig[jj] = (real)(disk[ii].Sigma [jj] * col_density2astro);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      tmp__zd[jj] = (real)(disk[ii].zd    [jj] *      length2astro);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    }/* for(int jj = 0; jj < nhorbin; jj++){ */
    for(int jj = 0; jj < nverbin; jj++)
      tmp_ver[jj] = (real)(disk[ii].ver[jj] * length2astro);
    for(int jj = 0; jj < nhorbin * nverbin; jj++){
      tmp_rho[jj] = (real)((*disk[ii].rho)[jj] * density2astro);
      tmp_Phi[jj] = (real)(  disk[ii].pot [jj] * senergy2astro);
    }/* for(int jj = 0; jj < nhorbin * nverbin; jj++){ */
    //---------------------------------------------------------------------
    success &= (fwrite(&disk[ii].cfg->rs, sizeof(double),                 1, fp) ==                 1);
    success &= (fwrite(&disk[ii].cfg->zd, sizeof(double),                 1, fp) ==                 1);
    /* success &= (fwrite( disk[ii].hor    , sizeof(double), nhorbin          , fp) == nhorbin          ); */
    /* success &= (fwrite( disk[ii].ver    , sizeof(double),           nverbin, fp) ==           nverbin); */
    /* success &= (fwrite(*disk[ii].rho    , sizeof(double), nhorbin * nverbin, fp) == nhorbin * nverbin); */
    /* success &= (fwrite( disk[ii].phi    , sizeof(double), nhorbin * nverbin, fp) == nhorbin * nverbin); */
    /* success &= (fwrite( disk[ii].sigma  , sizeof(double), nhorbin          , fp) == nhorbin          ); */
    /* success &= (fwrite( disk[ii].Sigma  , sizeof(double), nhorbin          , fp) == nhorbin          ); */
    success &= (fwrite(tmp_hor, sizeof(real), nhorbin          , fp) == nhorbin          );
    success &= (fwrite(tmp_ver, sizeof(real),           nverbin, fp) ==           nverbin);
    success &= (fwrite(tmp_rho, sizeof(real), nhorbin * nverbin, fp) == nhorbin * nverbin);
    success &= (fwrite(tmp_Phi, sizeof(real), nhorbin * nverbin, fp) == nhorbin * nverbin);
    success &= (fwrite(tmp_sig, sizeof(real), nhorbin          , fp) == nhorbin          );
    success &= (fwrite(tmp_Sig, sizeof(real), nhorbin          , fp) == nhorbin          );
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
    success &= (fwrite(tmp__zd, sizeof(real), nhorbin          , fp) == nhorbin          );
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    /* success &= (fwrite(  sph_rad        , sizeof(double), nsphbin          , fp) == nsphbin          ); */
    /* success &= (fwrite(  sph_enc        , sizeof(double), nsphbin          , fp) == nsphbin          ); */
    /* success &= (fwrite(  sph_rho        , sizeof(double), nsphbin          , fp) == nsphbin          ); */
    //---------------------------------------------------------------------
    fclose(fp);
    //---------------------------------------------------------------------
    if( success != true ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  //-----------------------------------------------------------------------
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  free(tmp_hor);  free(tmp_ver);
  free(tmp_sig);  free(tmp_Sig);
  free(_vcirc);  free(_sigmap);  free(_sigmaR);
  free(_kappa);  free(_Omega);  free(_lambda);  free(_toomre);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  free(tmp__zd);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  free(tmp_rho);  free(tmp_Phi);  free(rhoFrac);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //-----------------------------------------------------------------------
  /* read input arguments */
  //-----------------------------------------------------------------------
  if( argc < 9 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 9);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -Ntot=<unsigned long int>\n");
    __FPRINTF__(stderr, "          -config=<char *>\n");
    __FPRINTF__(stderr, "          -eps=<real> -eta=<real>\n");
    __FPRINTF__(stderr, "          -ft=<real> -snapshotInterval=<real> -saveInterval=<real>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }
  //-----------------------------------------------------------------------
  /* read input arguments do not depend on the unit system adopted in the numerical simulation */
  char *file;  requiredCmdArg(getCmdArgStr( argc, (const char * const *)argv,   "file", &file));
  char *fcfg;  requiredCmdArg(getCmdArgStr( argc, (const char * const *)argv, "config", &fcfg));
  ulong Ntot;  requiredCmdArg(getCmdArgUlg( argc, (const char * const *)argv,   "Ntot", &Ntot));
  real   eta;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,    "eta", &eta));
  double saveInterval;  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "saveInterval", &saveInterval));
  //-----------------------------------------------------------------------
  /* set unit system by reading the configuration file about physical parameters of the initial distribution */
  int unit, kind;
  profile_cfg *cfg;
  readProfileCfg(fcfg, &unit, &kind, &cfg);
  if( kind > NKIND_MAX ){    __KILL__(stderr, "ERROR: kind(= %d) must be smaller than %d\n", kind, NKIND_MAX);  }
  //-----------------------------------------------------------------------
  /* read input arguments depend on the unit system adopted in the numerical simulation */
  double tmp;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "snapshotInterval", &tmp));
  double snapshotInterval = ldexp(1.0, (int)floor(log2(tmp * time_astro2com)));
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "eps", &tmp));  real   eps = (real)(tmp * length_astro2com);
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "ft",  &tmp));  double  ft =       (tmp *   time_astro2com);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set number of particles */
  //-----------------------------------------------------------------------
  ulong Nrem = 0;
  for(int ii = 0; ii < kind; ii++)
    Nrem += (cfg[ii].forceNum == 1) ? cfg[ii].num : 0;
  if( Nrem > Ntot ){
    __KILL__(stderr, "ERROR: the sum of number of particles for each component (%zu) exceeds the specified total number of particles (%zu).\n", Nrem, Ntot);
  }/* if( Nrem > Ntot ){ */
  Nrem = Ntot - Nrem;
  //-----------------------------------------------------------------------
  /* set number of particles to represent each profile */
  double Mtot = 0.0;
  for(int ii = 0; ii < kind; ii++)
    if( cfg[ii].forceNum != 1 )
      Mtot += cfg[ii].Mtot;
  const double Minv = 1.0 / Mtot;
  ulong Nuse = 0;
  int idx = kind;
  double Mmax = 0.0;
  for(int ii = 0; ii < kind; ii++){
    //---------------------------------------------------------------------
    if( cfg[ii].forceNum != 1 ){
      //-------------------------------------------------------------------
      /* number of particles is determined by mass fraction */
      cfg[ii].num = (ulong)(cfg[ii].Mtot * Minv * (double)Nrem);
      Nuse += cfg[ii].num;
      //-------------------------------------------------------------------
      if( cfg[ii].Mtot > Mmax ){
	//-----------------------------------------------------------------
	idx = ii;
	Mmax = cfg[ii].Mtot;
	//-----------------------------------------------------------------
      }/* if( cfg[ii].Mtot > Mmax ){ */
      //-------------------------------------------------------------------
    }/* if( cfg[ii].forceNum != 1 ){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < kind; ii++){ */
  if( (idx == kind) && (Nuse != Nrem) ){
    __KILL__(stderr, "ERROR: mismatch about number of particles detected (Nuse = %zu, Nrem = %zu) with %d components\n", Nuse, Nrem, kind);
  }/* if( (idx == kind) && (Nuse != Nrem) ){ */
  if( Nuse != Nrem )
    cfg[idx].num += (Nrem - Nuse);
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
  /* initialize the table for Gaussian Quadrature provided by GSL */
  for(int ii = 0; ii < NTBL_GAUSS_QD; ii++){
    gsl_gaussQD_pos   [ii] = 0.0;
    gsl_gaussQD_weight[ii] = 0.0;
  }/* for(int ii = 0; ii < NTBL_GAUSS_QD; ii++){ */
  gsl_integration_glfixed_table *tab;
  tab = gsl_integration_glfixed_table_alloc(NINTBIN);
  int max = (NINTBIN >> 1) + (NINTBIN & 1);
  for(int ii = 0; ii < max; ii++){
    gsl_gaussQD_pos   [ii] = (*tab).x[(max - 1) - ii];
    gsl_gaussQD_weight[ii] = (*tab).w[(max - 1) - ii];
  }/* for(int ii = 0; ii < max; ii++){ */
  gsl_integration_glfixed_table_free(tab);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set global settings */
  //-----------------------------------------------------------------------
  /* identifier for disk components is negative value (-1, -2) */
  /* identifier for spherical components is positive value (0, 1, 2, ...) */
  int skind = kind;
  for(int ii = 0; ii < kind; ii++)
    if( cfg[ii].kind < 0 )
      skind--;
  const int ndisk = kind - skind;
  const bool addDisk = (ndisk != 0) ? true : false;
  if( addDisk )
    for(int ii = skind; ii < kind; ii++)
      if( cfg[ii].kind >= 0 ){      	__KILL__(stderr, "ERROR: disk component must be last component(s).\n\tModify \"%s/%s\".\n", CFGFOLDER, fcfg);      }
  //-----------------------------------------------------------------------
  bool cutoff = false;
  double rmax = 0.0;
  for(int ii = 0; ii < kind; ii++){
    cutoff |= cfg[ii].cutoff;
    if( cfg[ii].cutoff ){      if( rmax < cfg[ii].rc )	rmax = cfg[ii].rc;    }
    else{                      if( rmax < cfg[ii].rs )	rmax = cfg[ii].rs;    }
  }/* for(int ii = 0; ii < kind; ii++){ */
  //-----------------------------------------------------------------------
  if( cutoff )    rmax *= (double)(NINTBIN * 10);/* <-- sufficiently greater than outermost radius times NINTBIN */
  else            rmax *= 1000.0;
  const double logrmin = log10(MINRAD);
  const double logrmax = log10(rmax);
  const double logrbin = (logrmax - logrmin) / (double)(4 + NRADBIN);
  const double invlogrbin = 1.0 / logrbin;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set distribution function */
  //-----------------------------------------------------------------------
  /* memory allocation for spherical component(s) */
  profile **prf, *_prf;
  /* 2 * 2 bins are added in the both edge */
  _prf = (profile  *)malloc(sizeof(profile  ) * kind * (4 + NRADBIN));  if( _prf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _prf");  }
  prf  = (profile **)malloc(sizeof(profile *) * kind                );  if(  prf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  prf");  }
  for(int ii = 0; ii < kind; ii++)
    prf[ii] = _prf + ii * (4 + NRADBIN);
#pragma omp parallel for
  for(int ii = 0; ii < kind; ii++)
    for(int jj = 0; jj < 4 + NRADBIN; jj++)
      prf[ii][jj].rad = pow(10.0, logrmin + logrbin * (double)jj);
  /* memory allocation for disk component */
  disk_data  *disk_info;
  double *disk_hor, *disk_ver, *disk_pot, *disk_dPhidR, *disk_d2PhidR2;
  double *sph_rad, *sph_enc, *sph_rho;
  double *disk_rho, *disk_rhoSum, *disk_rhoTot, *disk_Sigma, *disk_sigmaz, *disk_enc;
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  double *disk_zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  if( addDisk ){
    //---------------------------------------------------------------------
    /* allocate required arrays */
    //---------------------------------------------------------------------
    allocDiskProfile(ndisk, &disk_info, &disk_hor, &disk_ver, &disk_pot, &disk_dPhidR, &disk_d2PhidR2, &sph_rad, &sph_rho, &sph_enc,
		     &disk_rho, &disk_rhoSum, &disk_rhoTot, &disk_Sigma, &disk_sigmaz, &disk_enc
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
		     , &disk_zd
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
		     );
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* commit fundamental information */
    //---------------------------------------------------------------------
    for(int ii = 0; ii < ndisk; ii++){
      disk_info[ii].cfg        = &cfg[skind + ii];
      disk_info[ii].prf        =  prf[skind + ii];
      disk_info[ii].   logrbin =    logrbin;
      disk_info[ii].invlogrbin = invlogrbin;
    }/* for(int ii = 0; ii < ndisk; ii++){ */
    //---------------------------------------------------------------------
  }/* if( addDisk ){ */
  //-----------------------------------------------------------------------
  /* set density profile and mass profile for spherical component(s) */
  for(int ii = 0; ii < skind; ii++){
    //---------------------------------------------------------------------
    profile_abel_cfg dummy;    dummy.invRd = 1.0;
    //---------------------------------------------------------------------
#if 0
    fprintf(stdout, "%d-th component: kind = %d\n", ii, cfg[ii].kind);
    fflush(NULL);
#endif
    //---------------------------------------------------------------------
    switch( cfg[ii].kind ){
    case   PLUMMER:      setDensityProfilePlummer        (prf[ii],  cfg[ii].rs                                                                       );      break;
    case      KING:      setDensityProfileKing           (prf[ii], &cfg[ii]                                                                          );      break;
    case   BURKERT:      setDensityProfileBurkert        (prf[ii],  cfg[ii].rs                                                                       );      break;
    case HERNQUIST:      setDensityProfileHernquist      (prf[ii],  cfg[ii].rs                                                                       );      break;
    case       NFW:      setDensityProfileNFW            (prf[ii],  cfg[ii].rs                                                                       );      break;
    case     MOORE:      setDensityProfileMoore          (prf[ii],  cfg[ii].rs                                                                       );      break;
    case   EINASTO:      setDensityProfileEinasto        (prf[ii],  cfg[ii].rs, cfg[ii]. einasto_alpha                                               );      break;
    case  APP_KING:      setDensityProfileAppKing        (prf[ii],  cfg[ii].rs,                         cfg[ii].   king_rt                           );      break;
    case TWO_POWER:      setDensityProfileTwoPower       (prf[ii],  cfg[ii].rs, cfg[ii].twopower_alpha, cfg[ii].twopower_beta, cfg[ii].twopower_gamma);      break;
    case TRI_POWER:      setDensityProfileTriPower       (prf[ii],  cfg[ii].rs, cfg[ii].tripower_rout, cfg[ii].twopower_alpha, cfg[ii].twopower_beta, cfg[ii].twopower_gamma, cfg[ii].tripower_delta, cfg[ii].tripower_epsilon);      break;
    case APP_EVANS:      setDensityProfileAppLoweredEvans(prf[ii],  cfg[ii].rs, cfg[ii].alevans_alpha, cfg[ii].alevans_rc, cfg[ii].alevans_beta, cfg[ii].alevans_rt, 1.0 / cfg[ii].alevans_wt);      break;
    case TABLE_RHO:      setDensityProfileTable          (prf[ii],  cfg[ii].rs, cfg[ii].table                                                        );      break;
    case TABLE_SIG:      readColumnDensityProfileTable   (prf[ii],  cfg[ii].rs, cfg[ii].table, cfg[ii]                                               );      break;
    case SPHSERSIC:      execAbelTransform               (prf[ii],  cfg[ii]   , MINRAD, rmax, dummy                                                  );      break;
    case SIGTWOPOW:      execAbelTransform               (prf[ii],  cfg[ii]   , MINRAD, rmax, dummy                                                  );      break;
    case CENTRALBH:      setContributionByCentralBH      (prf[ii],  cfg[ii]                                                                          );      break;
    default:      __KILL__(stderr, "ERROR: undefined profile ``%d'' specified\n", cfg[ii].kind);      break;
    }/* switch( cfg[ii].kind ){ */
    //---------------------------------------------------------------------
    if( cfg[ii].kind != CENTRALBH ){
#if 1
      integrateDensityProfile(prf[ii], logrbin, cfg[ii].Mtot, cfg[ii].cutoff, cfg[ii].rc, cfg[ii].rc_width);
#else
      integrateDensityProfile(prf[ii], logrbin, cfg[ii].Mtot, cfg[ii].cutoff, cfg[ii].rc, (cfg[ii].rs < (cfg[ii].rc * 0.1)) ? (cfg[ii].rc * 0.1) : (cfg[ii].rs));
#endif
    }/* if( cfg[ii].kind != CENTRALBH ){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < skind; ii++){ */
  //-----------------------------------------------------------------------
  /* evaluate sum of density, enclosed mass and potential of all spherical component(s) */
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    double rho = 0.0;
    double enc = 0.0;
    double psi = 0.0;
    for(int kk = 0; kk < skind; kk++){
      rho += prf[kk][ii].rho;
      enc += prf[kk][ii].enc;
      psi += prf[kk][ii].psi;
    }/* for(int kk = 0; kk < skind; kk++){ */
    //---------------------------------------------------------------------
    for(int kk = 0; kk < kind; kk++){
      prf[kk][ii].rho_tot = rho;
      prf[kk][ii].enc_tot = enc;
      prf[kk][ii].psi_tot = psi;
    }/* for(int kk = 0; kk < kind; kk++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
#if 0
  for(int ii = 0; ii < 4 + NRADBIN; ii += (NRADBIN >> 6)){
    fprintf(stdout, "%e", prf[0][ii].rad);
    for(int jj = 0; jj < skind; jj++)
      fprintf(stdout, "\t%e\t%e\t%e", prf[jj][ii].rho, prf[jj][ii].enc, prf[jj][ii].psi);
    fprintf(stdout, "\n");
  }
  exit(0);
#endif
  //-----------------------------------------------------------------------
  /* set density profile for the disk component */
  if( addDisk ){
    //---------------------------------------------------------------------
    /* set disk_radius, disk_height, disk_pot */
    makeDiskPotentialTable(ndisk, disk_info);
#if 0
    printf("# potential-density pair of the disk component is calculated\n");
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	fprintf(stderr, "%e\t%e\t%e\t%e\n", disk_info[0].hor[ii], disk_info[0].ver[jj], (*disk_info[0].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)], disk_info[0].pot[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)]);
      fprintf(stderr, "\n");
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
    exit(0);
#endif
    //---------------------------------------------------------------------
    /* set profile of spherical averaged density, mass and potential */
    integrateSphericalDensityProfile(ndisk, disk_info);
#if 0
    printf("# spherical averaged quantities are calculated\n");
    for(int ii = 0; ii < 4 + NRADBIN; ii++)
      fprintf(stderr, "%e\t%e\t%e\t%e\n", prf[skind][ii].rad, prf[skind][ii].rho, prf[skind][ii].enc, prf[skind][ii].psi);
    exit(0);
#endif
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < 4 + NRADBIN; ii++){
      //-------------------------------------------------------------------
      double rho = 0.0;
      double enc = 0.0;
      double psi = 0.0;
      for(int kk = skind; kk < kind; kk++){
	rho += prf[kk][ii].rho;
	enc += prf[kk][ii].enc;
	psi += prf[kk][ii].psi;
      }
      //-------------------------------------------------------------------
      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].rho_tot += rho;
	prf[kk][ii].enc_tot += enc;
	prf[kk][ii].psi_tot += psi;
      }/* for(int kk = 0; kk < skind; kk++){ */
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
    //---------------------------------------------------------------------
  }/* if( addDisk ){ */
  //-----------------------------------------------------------------------
  /* integrate Eddington's formula numerically */
  dist_func **fene, *_fene;
  _fene = (dist_func  *)malloc(sizeof(dist_func  ) * skind * NENEBIN);  if( _fene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _fene");  }
  fene  = (dist_func **)malloc(sizeof(dist_func *) * skind          );  if(  fene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  fene");  }
  for(int ii = 0; ii < skind; ii++)
    fene[ii] = _fene + ii * NENEBIN;
  integrateEddingtonFormula(skind, prf, fene);
  //-----------------------------------------------------------------------
  if( addDisk ){
    /* differentiate potential along the radial direction on the equatorial plane */
    diffAxisymmetricPotential(disk_info[0]);
#if 0
    fprintf(stderr, "#R\tz\tPhi\tdPhidR\td2PhidR2\n");
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
#ifndef USE_POTENTIAL_SCALING_SCHEME
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", disk_info[0].hor[ii], disk_info[0].ver[jj], disk_info[0].pot[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)], disk_info[0].dPhidR[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)], disk_info[0].d2PhidR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)]);
      fprintf(stderr, "\n");
#else///USE_POTENTIAL_SCALING_SCHEME
	fprintf(stderr, "%e\t%e\t%e\t%e\n", disk_info[0].hor[ii], disk_info[0].pot[ii], disk_info[0].dPhidR[ii], disk_info[0].d2PhidR2[ii]);
#endif//USE_POTENTIAL_SCALING_SCHEME
    }/* for(int ii = 0; ii < NDBIN_RAD; ii++){ */
    exit(0);
#endif
#if 0
    printf("# differentiation along the radial direction on the equatorial plane are calculated\n");
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
#ifndef USE_POTENTIAL_SCALING_SCHEME
      for(int jj = 0; jj < NDISKBIN_VER; jj++){
	const double Omega2 = disk_info[0].dPhidR[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] / disk_info[0].hor[ii];
	const double kappa = sqrt(disk_info[0].d2PhidR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] + 3.0 * Omega2);
	const double gam   = 2.0 * sqrt(Omega2 / (disk_info[0].d2PhidR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] + 3.0 * Omega2));
	const double vcirc  = disk_info[0].hor[ii] * sqrt(Omega2);
	fprintf(stderr, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", disk_info[0].hor[ii], disk_info[0].ver[jj], Omega2, kappa, gam, 1.0 / (0.5 * sqrt(3.0 + disk_info[0].d2PhidR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] / (1.0e-100 + Omega2))), vcirc);
      }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
      fprintf(stderr, "\n");
#else///USE_POTENTIAL_SCALING_SCHEME
      const double Omega2 = disk_info[0].dPhidR[ii] / disk_info[0].hor[ii];
      const double kappa = sqrt(disk_info[0].d2PhidR2[ii] + 3.0 * Omega2);
      const double gam   = 2.0 * sqrt(Omega2 / (disk_info[0].d2PhidR2[ii] + 3.0 * Omega2));
      const double vcirc  = disk_info[0].hor[ii] * sqrt(Omega2);
      fprintf(stderr, "%e\t%e\t%e\t%e\t%e\t%e\n", disk_info[0].hor[ii], Omega2, kappa, gam, 1.0 / (0.5 * sqrt(3.0 + disk_info[0].d2PhidR2[ii] / (1.0e-100 + Omega2))), vcirc);
#endif//USE_POTENTIAL_SCALING_SCHEME
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
    exit(0);
#endif
    /* set velocity dispersion in vertical direction */
    calcVerticalVdisp(ndisk, disk_info);
    for(int ii = skind; ii < kind; ii++){
      /* cfg[ii].vdispz0 = disk_vdisp[0]; */
      cfg[ii].vdispz0 = disk_info[ii - skind].sigmaz[0];
      if( cfg[ii].vdispR0 < 0.0 )
	cfg[ii].vdispR0 = cfg[ii].vdispz0;
    }/* for(int ii = skind; ii < kind; ii++){ */
#if 0
    printf("# vertical velocity dispersion are calculated\n");
    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
	fprintf(stderr, "%e\t%e\n", disk_info[0].hor[ii], disk_info[0].sigmaz[ii]);
    exit(0);
#endif
    //---------------------------------------------------------------------
    /* output fundamental quantities of the disk component */
    writeDiskData(file, ndisk, disk_info);
    //---------------------------------------------------------------------

/*     //--------------------------------------------------------------------- */
/*     /\* evaluate and output kinematical properties of the disk component *\/ */
/*     //--------------------------------------------------------------------- */
/*     for(int hh = skind; hh < kind; hh++){ */
/*       //------------------------------------------------------------------- */
/*       const int diskID = hh - skind; */
/*       //------------------------------------------------------------------- */
/*       FILE *fp; */
/*       char filename[256]; */
/*       sprintf(filename, "%s/%s.diskhor.%d.txt", DATAFOLDER, file, diskID); */
/*       fp = fopen(filename, "w"); */
/*       if( fp == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);      } */
/*       fprintf(fp, "#R\tSigma(R)\tfrac(R)\tv_c(R)\tsigma_p(R)\tsigma_R(R)\tsigma_z(R)\tkappa(R)\tOmega(R)\tToomre-Q(R)\tlambda_crit(R)\n"); */
/*       //------------------------------------------------------------------- */
/* #ifndef USE_ORIGINAL_VDISP_ESTIMATOR */
/*       const double invRd = 1.0 / cfg[hh].rs; */
/* #endif//USE_ORIGINAL_VDISP_ESTIMATOR */
/*       cfg[hh].vcirc_max   = -1.0; */
/*       cfg[hh].vcirc_max_R = -1.0; */
/*       cfg[hh].Qmin = DBL_MAX; */
/*       bool passed = false; */
/*       //------------------------------------------------------------------- */
/* 	for(int ii = 0; ii < NDISKBIN_HOR; ii++) */
/* 	  { */
/* 	    //------------------------------------------------------------- */
/* 	    /\* evaluate epicyclic quantities and circular speed *\/ */
/* #ifndef USE_POTENTIAL_SCALING_SCHEME */
/* 	    const double Omega = sqrt(disk_info[diskID]. dPhidR [INDEX2D(NDISKBIN_RAD, NDISKBIN_VER, ii, 0)] / disk_info[diskID].hor[ii]); */
/* 	    const double kappa = sqrt(disk_info[diskID].d2PhidR2[INDEX2D(NDISKBIN_RAD, NDISKBIN_VER, ii, 0)] + 3.0 * Omega * Omega); */
/* #else///USE_POTENTIAL_SCALING_SCHEME */
/* 	    const double Omega = sqrt(disk_info[diskID]. dPhidR [ii] / disk_info[diskID].hor[ii]); */
/* 	    const double kappa = sqrt(disk_info[diskID].d2PhidR2[ii] + 3.0 * Omega * Omega); */
/* #endif//USE_POTENTIAL_SCALING_SCHEME */
/* 	    const double vcirc = disk_info[diskID].hor[ii] * Omega; */
/* 	    //------------------------------------------------------------- */
/* 	    /\* evaluate Toomre's Q-value *\/ */
/* #ifdef  USE_ORIGINAL_VDISP_ESTIMATOR */
/* 	    const double sigmap = DISK_PERP_VDISP(disk_info[diskID].sigmaz[ii], vcirc, cfg[hh].vdisp_frac); */
/* 	    const double gam2inv = 0.25 * (3.0 + disk_info[diskID].d2PhidR2[ii] / (1.0e-100 + Omega * Omega)); */
/* 	    const double sigmaR = sigmap / sqrt(gam2inv); */
/* #else///USE_ORIGINAL_VDISP_ESTIMATOR */
/* 	    const double sigmaR = DISK_RADIAL_VDISP(cfg[hh].vdispR0, disk_info[diskID].hor[ii], invRd); */
/* #endif//USE_ORIGINAL_VDISP_ESTIMATOR */
/* 	    const double toomre = sigmaR * kappa / (3.36 * newton * disk_info[diskID].Sigma[ii]); */
/* 	    const double lambda = 4.0 * M_PI * M_PI * newton * disk_info[diskID].Sigma[ii] / (kappa * kappa); */
/* 	    //------------------------------------------------------------- */
/* 	    /\* find the maximum circular speed *\/ */
/* 	    if( vcirc > cfg[hh].vcirc_max ){ */
/* 	      cfg[hh].vcirc_max   = vcirc; */
/* 	      cfg[hh].vcirc_max_R = disk_info[diskID].hor[ii]; */
/* 	    }/\* if( vcirc > cfg[hh].vcirc_max ){ *\/ */
/* 	    //------------------------------------------------------------- */
/* 	    /\* find the minimum Toomre's Q-value *\/ */
/* 	    if( (disk_info[diskID].Sigma[ii] > 1.0e-4 * disk_info[diskID].Sigma[0]) && (toomre > DBL_EPSILON) && (toomre < cfg[hh].Qmin) ) */
/* 	      cfg[hh].Qmin = toomre; */
/* 	    //------------------------------------------------------------- */
/* 	    /\* find the circular speed and Toomre's Q-value at the scale radius *\/ */
/* 	    if( !passed ){ */
/* 	      cfg[hh].vcirc_Rd = vcirc; */
/* 	      cfg[hh].toomre   = toomre; */
/* 	      if( disk_info[diskID].hor[ii] > cfg[hh].rs ) */
/* 		passed = true; */
/* 	    }/\* if( !passed ){ *\/ */
/* 	    //------------------------------------------------------------- */
/* 	    double rhoTot = DBL_MIN; */
/* 	    for(int jj = 0; jj < ndisk; jj++) */
/* 	      rhoTot += (*disk_info[jj].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 0)]; */
/* 	    //------------------------------------------------------------- */
/* 	    fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", disk_info[diskID].hor[ii], disk_info[diskID].Sigma[ii], (*disk_info[diskID].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 0)] / rhoTot, vcirc, sigmap, sigmaR, disk_info[diskID].sigmaz[ii], kappa, Omega, toomre, lambda); */
/* 	    //------------------------------------------------------------- */
/* 	  } */
/*       //------------------------------------------------------------------- */
/*       fclose(fp); */
/*       //------------------------------------------------------------------- */
/*     }/\* for(int hh = skind; hh < kind; hh++){ *\/ */
    //---------------------------------------------------------------------
  }/* if( addDisk ){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set particle distribution */
  //-----------------------------------------------------------------------
  nbody_particle *body;
  allocParticleDataAoS((int)Ntot, &body);
  /* const real mass = (real)(Mtot / (double)Ntot); */
  //-----------------------------------------------------------------------
  /* create spherical particle distribution */
  Nuse = 0;
  for(int ii = 0; ii < skind; ii++){
    //---------------------------------------------------------------------
    /* distribute spheroid particles */
    if( cfg[ii].kind != CENTRALBH )
      cfg[ii].Ecut = distributeSpheroidParticles(&Nuse, body, (real)(cfg[ii].Mtot / (double)cfg[ii].num), cfg[ii], &prf[ii][2], fene[ii]);
    else{
      body[Nuse]. x = body[Nuse]. y = body[Nuse]. z = ZERO;
      body[Nuse].vx = body[Nuse].vy = body[Nuse].vz = ZERO;
      body[Nuse].ax = body[Nuse].ay = body[Nuse].az = ZERO;
      body[Nuse].m   = (real)cfg[ii].Mtot;
      body[Nuse].pot = ZERO;
      body[Nuse].idx = Nuse;
      Nuse++;
    }/* else{ */
    //---------------------------------------------------------------------
    /* shift center-of-mass, remove bulk motion */
#ifdef  SHIFT_CENTER_PER_COMPONENT
    shiftCenter(cfg[ii].num, &body[Nuse - cfg[ii].num]);
#endif//SHIFT_CENTER_PER_COMPONENT
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < skind; ii++){ */
  //-----------------------------------------------------------------------
  /* add disk component if required */
  if( addDisk ){
    //---------------------------------------------------------------------
#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
    double Krand = 0.0;
    for(ulong ii = 0; ii < Nuse; ii++)
      Krand += body[ii].vx * body[ii].vx + body[ii].vy * body[ii].vy + body[ii].vz * body[ii].vz;
    for(int ii = 0; ii < ndisk; ii++)
      disk_info[ii].Krand_sph = Krand;
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION
    //---------------------------------------------------------------------
    for(int ii = 0; ii < ndisk; ii++){
      //-------------------------------------------------------------------
      /* distribute disk particles */
      distributeDiskParticles(&Nuse, body, (real)(disk_info[ii].cfg->Mtot / (double)disk_info[ii].cfg->num), disk_info[ii]);
      //-------------------------------------------------------------------
      /* shift center-of-mass, remove bulk motion */
#ifdef  SHIFT_CENTER_PER_COMPONENT
      shiftCenter(disk_info[ii].cfg->num, &body[Nuse - disk_info[ii].cfg->num]);
#endif//SHIFT_CENTER_PER_COMPONENT
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < ndisk; ii++){ */
    //---------------------------------------------------------------------
  }/* if( addDisk ){ */
  //-----------------------------------------------------------------------
#ifndef SHIFT_CENTER_PER_COMPONENT
  shiftCenter(Ntot, body);
#endif//SHIFT_CENTER_PER_COMPONENT
  //-----------------------------------------------------------------------
  if( skind == 0 )
    for(int ii = 0; ii < 4 + NRADBIN; ii++){
      prf[0][ii].enc_tot  = prf[ 0][ii].enc;
      for(int jj = 1; jj < ndisk; jj++)
      prf[0][ii].enc_tot += prf[jj][ii].enc;
    }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
  //-----------------------------------------------------------------------
  /* write fundamental information */
  outputFundamentalInformation(unit, kind, skind, cfg, prf, fene, Ntot, eps, snapshotInterval, ft, file);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  double time = 0.0;
  double  dt  = 0.0;
  int   last  = 1;
  ulong steps = 0;
  //-----------------------------------------------------------------------
  writeSettings(unit, Ntot, eps, eta, ft, snapshotInterval, saveInterval, file);
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
#ifdef  MONITOR_ENERGY_ERROR
  static energyError relEneErr = {1.0, DBL_MIN};
#endif//MONITOR_ENERGY_ERROR
  writeTentativeData(time, dt, steps, Ntot, body, file, &last, hdf5type
#ifdef  MONITOR_ENERGY_ERROR
		     , relEneErr
#endif//MONITOR_ENERGY_ERROR
		     );
  removeHDF5DataType(hdf5type);
#else///USE_HDF5_FORMAT
  writeTentativeData(time, dt, steps, Ntot, body, file, &last);
#endif//USE_HDF5_FORMAT
  updateConfigFile(last, file);
  //-----------------------------------------------------------------------
#ifdef  APPEND_ASCII_ICDATA
  FILE *fpascii;
  char asciifile[256];
  ulong Nhead = 0;
  for(int kk = 0; kk < kind; kk++){
    //---------------------------------------------------------------------
    sprintf(asciifile, "%s/%s%d_ascii.dat", DATAFOLDER, file, kk);
    fpascii = fopen(asciifile, "w");
    if( fpascii == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", asciifile);      }
    //---------------------------------------------------------------------
    for(ulong ii = 0; ii < cfg[kk].num; ii++)
      fprintf(fpascii, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", body[Nhead + ii].x, body[Nhead + ii].y, body[Nhead + ii].z, body[Nhead + ii].vx, body[Nhead + ii].vy, body[Nhead + ii].vz, body[Nhead + ii].m);
    Nhead += cfg[kk].num;
    //---------------------------------------------------------------------
    fclose(fpascii);
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < kind; kk++){ */
#endif//APPEND_ASCII_ICDATA
  //-----------------------------------------------------------------------
#ifdef  DUMPFILE_FOR_BONSAI
  FILE *fpbonsai;
  char bonsaifile[256];
  sprintf(bonsaifile, "%s/%s_bonsai.dat", DATAFOLDER, file);
  fpbonsai = fopen(bonsaifile, "wb");
  if( fpbonsai == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", bonsaifile);      }
  bool bonsaiSuccess = true;
  size_t bonsaiCount;
  /* header */
  struct dump {
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
  };
  typedef struct dump header;
  header bonsaiHeader;
  bonsaiHeader.time = time;
  bonsaiHeader.nbodies = (int)Ntot;
  bonsaiHeader.ndim = 3;
  bonsaiHeader.nsph = 0;
  bonsaiHeader.ndark = (int)Ntot;
  bonsaiHeader.nstar = 0;
  bonsaiCount = 1;  if( bonsaiCount != fwrite(&bonsaiHeader, sizeof(header), bonsaiCount, fpbonsai) )    bonsaiSuccess = false;
  /* main body */
  struct dark_particle {
    real mass;
    real pos[3];
    real vel[3];
    real eps;
    int idx;
  };
  for(ulong ii = 0; ii < Ntot; ii++){
    struct dark_particle tmp;
    tmp.mass   = body[ii]. m;
    tmp.pos[0] = body[ii]. x;
    tmp.pos[1] = body[ii]. y;
    tmp.pos[2] = body[ii]. z;
#if 1
    tmp.vel[0] = body[ii].vx;
    tmp.vel[1] = body[ii].vy;
    tmp.vel[2] = body[ii].vz;
#else
    tmp.vel[0] = 0.0f;
    tmp.vel[1] = 0.0f;
    tmp.vel[2] = 0.0f;
#endif
    tmp.eps    = eps;
    tmp.idx    = (int)body[ii].idx;
    bonsaiCount = 1;    if( bonsaiCount != fwrite(&tmp, sizeof(struct dark_particle), bonsaiCount, fpbonsai) )      bonsaiSuccess = false;
  }
  if( bonsaiSuccess != true ){
    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", bonsaifile);
  }
  fclose(fpbonsai);
#endif//DUMPFILE_FOR_BONSAI
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  freeParticleDataAoS(body);
  //-----------------------------------------------------------------------
  free(cfg);
  free(prf);
  free(_prf);
  free(fene);
  free(_fene);
  //-----------------------------------------------------------------------
  if( addDisk )
    freeDiskProfile(ndisk, disk_info, disk_hor, disk_ver, disk_pot, disk_dPhidR, disk_d2PhidR2, sph_rad, sph_rho, sph_enc,
		    disk_rho, disk_rhoSum, disk_rhoTot, disk_Sigma, disk_sigmaz, disk_enc
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
		    , disk_zd
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
		    );
  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
