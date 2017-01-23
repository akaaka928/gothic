/*************************************************************************\
 *                                                                       *
                  last updated on 2017/01/18(Wed) 18:00:56
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
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//-------------------------------------------------------------------------
#include "macro.h"
#include "constants.h"
//-------------------------------------------------------------------------
#include "profile.h"
//-------------------------------------------------------------------------
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
#include "magi.h"
#endif//MAKE_COLUMN_DENSITY_PROFILE
//-------------------------------------------------------------------------
extern const real newton;
extern const double     mass_astro2com;
extern const double   length_astro2com;
extern const double velocity_astro2com;
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 2; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = NRADBIN + 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* multiply overall factors */
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    prf[ii].enc *= 4.0 * M_PI;
    prf[ii].psi *= 4.0 * M_PI * (double)newton;
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
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
  if( cfg.kind == TBL_DISK )
    fprintf(stderr, "\ttable<char *>: file name to be read column density profile in table form\n");
  if( (cfg.kind == EXP_DISK) || (cfg.kind == SERSIC) || (cfg.kind == TBL_DISK) ){
    if( cfg.kind == SERSIC )
      fprintf(stderr, "\tn_sersic<real>: Sersic index\n");
    fprintf(stderr, "\tRt<real> Rt_width<real>: cutoff radius and width of the disk mid-plane density in horizontal direction in astrophysical units, respectively\n");
    fprintf(stderr, "\tzd<real>: scale height of isothermal disk in the vertical direction in astrophysical units\n");
    fprintf(stderr, "\tsigmaR0<real> frac<real>: velocity dispersion in radial direction at the center in the vertical direction in astrophysical units, perpendicular velocity dispersion over circular velocity\n");
    fprintf(stderr, "\t\tif the inputted sigmaR0 is negative, then the default value (= sigma_z(R = 0)) is substituted\n");
    fprintf(stderr, "\t\tif USE_ORIGINAL_VDISP_ESTIMATOR defined in src/init/disk_polar.h is ON, then frac is used; if that is OFF, then sigmaR0 is used.\n");
    fprintf(stderr, "\t\tif the inputted retrogradeFrac<real> is not zero ([0., 1.]), then rotation axis of particles with the given fraction is anti-parallel with normal component.\n");
  }/* if( (cfg.kind == EXP_DISK) || (cfg.kind == SERSIC) || (cfg.kind == TBL_DISK) ){ */
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
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  int checker = 1;
  //-----------------------------------------------------------------------
  /* read the specified unit system and set it */
  checker &= (1 == fscanf(fp, "%d", unit));
  setPhysicalConstantsAndUnitSystem(*unit, 1);
  /* read # of components */
  checker &= (1 == fscanf(fp, "%d", kind));
  //-----------------------------------------------------------------------
  *cfg = (profile_cfg *)malloc(sizeof(profile_cfg) * (*kind));  if( *cfg == NULL ){    __KILL__(stderr, "ERROR: failure to allocate cfg\n");  }
  for(int ii = 0; ii < *kind; ii++)
    checker &= (4 == fscanf(fp, "%d\t%s\t%d\t%zu", &(*cfg)[ii].kind, (*cfg)[ii].file, &(*cfg)[ii].forceNum, &(*cfg)[ii].num));
    /* checker &= (2 == fscanf(fp, "%d\t%s", &(*cfg)[ii].kind, (*cfg)[ii].file)); */
  //-----------------------------------------------------------------------
  fclose(fp);
  if( !checker ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }/* if( !checker ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* read individual settings */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < *kind; ii++){
    //---------------------------------------------------------------------
    sprintf(filename, "%s/%s", CFGFOLDER, (*cfg)[ii].file);
    //---------------------------------------------------------------------
    fp = fopen(filename, "r");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
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
    /* parameter for column density profile in table form */
    if( (*cfg)[ii].kind == TBL_DISK )
      checker &= (1 == fscanf(fp, "%s", (*cfg)[ii].table));
    //---------------------------------------------------------------------
    /* parameters for disk component */
    if( ((*cfg)[ii].kind == EXP_DISK) || ((*cfg)[ii].kind == SERSIC) || ((*cfg)[ii].kind == TBL_DISK) ){
      if( (*cfg)[ii].kind == SERSIC ){
	checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].n_sersic));
	(*cfg)[ii].b_sersic = getAsymptoticSersicScale((*cfg)[ii].n_sersic);
      }/* if( (*cfg)[ii].kind == SERSIC ){ */
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].zd));      (*cfg)[ii].zd   *= length_astro2com;
      checker &= (2 == fscanf(fp, "%le %le", &(*cfg)[ii].vdispR0, &(*cfg)[ii].vdisp_frac));      (*cfg)[ii].vdispR0 *= velocity_astro2com;
      checker &= (1 == fscanf(fp, "%le", &(*cfg)[ii].retrogradeFrac));
    }/* if( ((*cfg)[ii].kind == EXP_DISK) || ((*cfg)[ii].kind == SERSIC) || ((*cfg)[ii].kind == TBL_DISK) ){ */
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
  }/* for(int ii = 0; ii < *kind; ii++){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
//-------------------------------------------------------------------------
extern double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
//-------------------------------------------------------------------------
static inline void findIdx(const double rad, profile * restrict prf, int * restrict ll, int * restrict rr)
{
  //-----------------------------------------------------------------------
  bool bisection = true;
  *ll =               0;
  *rr = 4 + NRADBIN - 1;
  //-----------------------------------------------------------------------
  if( bisection == true )    if( fabs(prf[*ll].rad - rad) / rad < DBL_EPSILON ){      bisection = false;      *rr = (*ll) + 1;    }
  if( bisection == true )    if( fabs(prf[*rr].rad - rad) / rad < DBL_EPSILON ){      bisection = false;      *ll = (*rr) - 1;    }
  //-----------------------------------------------------------------------
  while( bisection ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)(*ll) + (uint)(*rr)) >> 1;
    //---------------------------------------------------------------------
    if( (prf[cc].rad - rad) * (prf[*ll].rad - rad) <= 0.0 )      *rr = (int)cc;
    else                                                         *ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + (*ll)) == (*rr) )      break;
    //---------------------------------------------------------------------
  }/* while( bisection ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void getInterpolatedDensity(const double RR, const double zz, const int skind, profile * restrict * prf, double *val)
{
  //-----------------------------------------------------------------------
  const double rad = SQRT(RR * RR + zz * zz);
  //-----------------------------------------------------------------------
  int ll, rr;
  findIdx(rad, prf[0], &ll, &rr);
  //-----------------------------------------------------------------------
  /* based on linear interpolation */
  const double alpha = (rad - prf[0][ll].rad) / (prf[0][rr].rad - prf[0][ll].rad);
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < skind; kk++)
    val[kk] = (UNITY - alpha) * prf[kk][ll].rho + alpha * prf[kk][rr].rho;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void gaussQuadVertical(const double RR, const int num, const double zmin, const double zmax, const int skind, profile * restrict * prf, double *sum, double *fm, double *fp);
void gaussQuadVertical(const double RR, const int num, const double zmin, const double zmax, const int skind, profile * restrict * prf, double *sum, double *fm, double *fp)
{
  //-----------------------------------------------------------------------
  const double mns = 0.5 * (zmax - zmin);
  const double pls = 0.5 * (zmax + zmin);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int kk = 0; kk < skind; kk++)    sum[kk] = 0.0;
  //-----------------------------------------------------------------------
  if( num & 1 ){
    const double weight =             gsl_gaussQD_weight[(num >> 1)];
    const double  value = pls + mns * gsl_gaussQD_pos   [(num >> 1)];
    getInterpolatedDensity(RR, value, skind, prf, fm);
    for(int kk = 0; kk < skind; kk++)
      sum[kk] = weight * fm[kk];
  }/* if( num & 1 ){ */
  //-----------------------------------------------------------------------
  for(int ii = (num >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double weight = gsl_gaussQD_weight[ii];
    const double zp = pls + mns * gsl_gaussQD_pos[ii];
    const double zm = pls - mns * gsl_gaussQD_pos[ii];
    getInterpolatedDensity(RR, zp, skind, prf, fp);
    getInterpolatedDensity(RR, zm, skind, prf, fm);
    //---------------------------------------------------------------------
    for(int kk = 0; kk < skind; kk++)
      sum[kk] += weight * (fm[kk] + fp[kk]);
    //---------------------------------------------------------------------
  }/* for(int ii = (num >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < skind; kk++)
    sum[kk] *= mns;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void calcColumnDensityProfile(const int skind, profile **prf, const double logrmax, profile_cfg *cfg)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  double rs = DBL_MAX;
  for(int ii = 0; ii < skind; ii++)
    if( cfg[ii].kind != CENTRALBH )
      rs = fmin(rs, cfg[ii].rs);
  //-----------------------------------------------------------------------
  const double Rmax = pow(10.0, logrmax);
  //-----------------------------------------------------------------------
#pragma omp parallel
  {
    //---------------------------------------------------------------------
    double *sum;    sum = (double *)malloc(skind * sizeof(double));    if( sum == NULL ){      __KILL__(stderr, "ERROR: failure to allocate sum\n");    }
    double *tfp;    tfp = (double *)malloc(skind * sizeof(double));    if( tfp == NULL ){      __KILL__(stderr, "ERROR: failure to allocate tfp\n");    }
    double *tfm;    tfm = (double *)malloc(skind * sizeof(double));    if( tfm == NULL ){      __KILL__(stderr, "ERROR: failure to allocate tfm\n");    }
    //---------------------------------------------------------------------
#pragma omp for schedule(auto) nowait
    for(int ii = 0; ii < 4 + NRADBIN; ii += SKIP_INTERVAL_FOR_COLUMN_DENSITY){
      //-------------------------------------------------------------------
      /* initialization */
      for(int kk = 0; kk < skind; kk++)
	prf[kk][ii].Sigma = 0.0;
      //-------------------------------------------------------------------
      const double RR = prf[0][ii].rad;
      const double R0 = fmin(RR, rs);
      const double R2 = fmax(RR, rs);
      const double R1 = 0.5 * (R0 + R2);
      //-------------------------------------------------------------------
      gaussQuadVertical(RR, NINTBIN,  0.0     ,          R0, skind, prf, sum, tfm, tfp);      for(int kk = 0; kk < skind; kk++)	prf[kk][ii].Sigma += sum[kk];
      gaussQuadVertical(RR, NINTBIN,	    R0,		 R1, skind, prf, sum, tfm, tfp);      for(int kk = 0; kk < skind; kk++)	prf[kk][ii].Sigma += sum[kk];
      gaussQuadVertical(RR, NINTBIN,	    R1,		 R2, skind, prf, sum, tfm, tfp);      for(int kk = 0; kk < skind; kk++)	prf[kk][ii].Sigma += sum[kk];
      gaussQuadVertical(RR, NINTBIN,	    R2,	 2.0 *	 R2, skind, prf, sum, tfm, tfp);      for(int kk = 0; kk < skind; kk++)	prf[kk][ii].Sigma += sum[kk];
      gaussQuadVertical(RR, NINTBIN,  2.0 * R2, 10.0 *	 R2, skind, prf, sum, tfm, tfp);      for(int kk = 0; kk < skind; kk++)	prf[kk][ii].Sigma += sum[kk];
      gaussQuadVertical(RR, NINTBIN, 10.0 * R2,	 2.0 * Rmax, skind, prf, sum, tfm, tfp);      for(int kk = 0; kk < skind; kk++)	prf[kk][ii].Sigma += sum[kk];
      //-------------------------------------------------------------------
      for(int kk = 0; kk < skind; kk++)
	prf[kk][ii].Sigma *= 2.0;
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
    //---------------------------------------------------------------------
    free(sum);    free(tfp);    free(tfm);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#   if  SKIP_INTERVAL_FOR_COLUMN_DENSITY != 4
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN - SKIP_INTERVAL_FOR_COLUMN_DENSITY; ii += SKIP_INTERVAL_FOR_COLUMN_DENSITY)
    for(int kk = 0; kk < skind; kk++){
      //-------------------------------------------------------------------
      const double S0 = prf[kk][ii                                   ].Sigma;
      const double S1 = prf[kk][ii + SKIP_INTERVAL_FOR_COLUMN_DENSITY].Sigma;
      //-------------------------------------------------------------------
      double slope = (S1 - S0) / (double)SKIP_INTERVAL_FOR_COLUMN_DENSITY;
      for(int jj = 1; jj < SKIP_INTERVAL_FOR_COLUMN_DENSITY; jj++)
	prf[kk][ii + jj].Sigma = S0 + slope * (double)jj;
      //-------------------------------------------------------------------
    }/* for(int kk = 0; kk < skind; kk++){ */
#else///SKIP_INTERVAL_FOR_COLUMN_DENSITY != 4
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii += 4)
    for(int kk = 0; kk < skind; kk++){
      //-------------------------------------------------------------------
      const double S0 = prf[kk][ii    ].Sigma;
      const double S1 = prf[kk][ii + 4].Sigma;
      //-------------------------------------------------------------------
      prf[kk][ii + 1].Sigma = 0.75 *  S0 + 0.25 * S1;
      prf[kk][ii + 2].Sigma = 0.5  * (S0 +        S1);
      prf[kk][ii + 3].Sigma = 0.25 *  S0 + 0.75 * S1;
      //-------------------------------------------------------------------
    }/* for(int kk = 0; kk < skind; kk++){ */
#endif//SKIP_INTERVAL_FOR_COLUMN_DENSITY != 4
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < NRADBIN; ii += 128)
    fprintf(stderr, "%e\t%e\t%e\n", prf[0][ii].rad, prf[0][ii].rho, prf[0][ii].Sigma);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//MAKE_COLUMN_DENSITY_PROFILE
//-------------------------------------------------------------------------
