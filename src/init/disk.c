/*************************************************************************\
 *                                                                       *
                  last updated on 2016/08/09(Tue) 18:01:27
 *                                                                       *
 *    Making Initial Condition Code of N-body Simulation                 *
 *       Assume balance of force in R and z direction                    *
 *          z-direction: gravity versus gradient of velocity dispersion  *
 *          R-direction: gravity versus centrifugal force                *
 *       Adopt epicyclic approximation                                   *
 *          connects velocity dispersion in R and phi direction          *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
/* #define CONFIRM_BUILDING_BLOCK */
//-------------------------------------------------------------------------
#define PROHIBIT_EXTRAPOLATION
//-------------------------------------------------------------------------
#define USE_ELLIPTIC_INTEGRAL
//-------------------------------------------------------------------------
#define SOPHISTICATED_REFINEMENT
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
//-------------------------------------------------------------------------
#include <gsl/gsl_rng.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#       ifdef  USE_ELLIPTIC_INTEGRAL
#include <gsl/gsl_sf_ellint.h>
#       endif//USE_ELLIPTIC_INTEGRAL
//-------------------------------------------------------------------------
#include <macro.h>
#include <constants.h>
//-------------------------------------------------------------------------
#include "blas.h"
#include "spline.h"
#include "magi.h"
#include "disk.h"
//-------------------------------------------------------------------------
extern const real newton;
//-------------------------------------------------------------------------
extern double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
//-------------------------------------------------------------------------
#define NMAX_GAUSS_QD_LOW (11)
#define NTBL_GAUSS_QD_LOW ((NMAX_GAUSS_QD_LOW >> 1) + (NMAX_GAUSS_QD_LOW & 1))
#define NINTBIN_LOW NMAX_GAUSS_QD_LOW
static double gsl_gaussQD_low_pos[NTBL_GAUSS_QD_LOW], gsl_gaussQD_low_weight[NTBL_GAUSS_QD_LOW];
//-------------------------------------------------------------------------
extern gsl_rng *GSLRand;
#define UNIRAND_DBL (      gsl_rng_uniform(GSLRand))
#define UNIRAND     ((real)gsl_rng_uniform(GSLRand))
#define RANDVAL_DBL (2.0 * (UNIRAND_DBL) - 1.0)
#define RANDVAL     (TWO * (UNIRAND    ) - UNITY)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef __ICC
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#     pragma warning (disable:869)
#endif//__ICC
//-------------------------------------------------------------------------
static inline double getSmoothCutoff (const double RR, const double Rt, const double invDelta){  return (0.5 * erfc(2.0 * (RR - Rt) * invDelta));}
#if 1
double getColumnDensityExp   (double RR, double invRd, const disk_util disk){
  return (exp(-                    RR * invRd                   ) * getSmoothCutoff(RR, disk.Rcutoff, disk.invRsmooth));}
double getColumnDensitySersic(double RR, double invRd, const disk_util disk){
  return (exp(-disk.sersic_b * pow(RR * invRd, disk.sersic_ninv)) * getSmoothCutoff(RR, disk.Rcutoff, disk.invRsmooth));}
#else
double getColumnDensityExp   (double RR, double invRd, double ninv, double bb, double Rt, double invDelta){
  return (exp(-         RR * invRd       ) * getSmoothCutoff(RR, Rt, invDelta));}
double getColumnDensitySersic(double RR, double invRd, double ninv, double bb, double Rt, double invDelta){
  return (exp(-bb * pow(RR * invRd, ninv)) * getSmoothCutoff(RR, Rt, invDelta));}
#endif
//-------------------------------------------------------------------------
static double invzd2invz0;
static inline void setz0inv4SechProfile(void){  invzd2invz0 = 2.0 * acosh(sqrt(M_E));}
#if 1
double getVerticalDensity(const double  zz, const double invzd, const disk_util disk){
  const double tmp = 1.0 / cosh(0.5 * zz * invzd * invzd2invz0);/* := sech() */
  return (tmp * tmp);
}
#else
double getVerticalDensity(const double  zz, const double invzd, const double ninv, const double bb, const double Rc, const double invRsmooth){
  const double tmp = 1.0 / cosh(0.5 * zz * invzd * invzd2invz0);/* := sech() */
  return (tmp * tmp);
}
#endif
//-------------------------------------------------------------------------
#ifdef __ICC
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#     pragma warning ( enable:869)
#endif//__ICC
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
static inline int bisection(const double val, const int num, double * restrict tab, const bool logtbl, const double invbin, double * restrict ratio)
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
      *ratio = (logtbl ? (log10(val / tab[ll])) : (val - tab[ll])) * invbin;
      return (ll);
    }/* if( (1 + ll) == rr ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline int findIdxSpherical(const double rad, profile *prf, const double invlogbin, double *ratio)
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = 3 + NRADBIN;
  //-----------------------------------------------------------------------
  if( rad < prf[ll].rad + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( rad > prf[rr].rad - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)ll + (uint)rr) >> 1;
    //---------------------------------------------------------------------
    if( (prf[cc].rad - rad) * (prf[ll].rad - rad) <= 0.0 )      rr = (int)cc;
    else                                                        ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr ){
      *ratio = log10(rad / prf[ll].rad) * invlogbin;
      return (ll);
    }/* if( (1 + ll) == rr ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double Phi_spherical(const double rad, profile *sph, const double invlogrbin_sph)
{
  //-----------------------------------------------------------------------
  double ratio;
  const int idx = findIdxSpherical(rad, sph, invlogrbin_sph, &ratio);
  //-----------------------------------------------------------------------
  return ((1.0 - ratio) * sph[idx].psi_tot + ratio * sph[1 + idx].psi_tot);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
//-------------------------------------------------------------------------
static inline int findIdxSphericalPsi(const double psi, profile *prf, double *ratio)
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = 3 + NRADBIN;
  //-----------------------------------------------------------------------
  if( psi > prf[ll].psi_tot - DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( psi < prf[rr].psi_tot + DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)ll + (uint)rr) >> 1;
    //---------------------------------------------------------------------
    if( (prf[cc].psi_tot - psi) * (prf[ll].psi_tot - psi) <= 0.0 )      rr = (int)cc;
    else                                                                ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr ){
      *ratio = (psi - prf[ll].psi_tot) / (prf[rr].psi_tot - prf[ll].psi_tot);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static const double diskDimmingHeight    = 16.0;
static const double diskDimmingHeightInv = 0.0625;
//-------------------------------------------------------------------------
static inline void getVariableDiskScaleHeight(const int ndisk, disk_data * restrict disk, const int NR, profile * restrict sph, const double invlogrbin_sph)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int hh = 0; hh < ndisk; hh++){
    //---------------------------------------------------------------------
    double *RR;    RR = disk[hh].hor;
    double *zd;    zd = disk[hh].zd;
    //---------------------------------------------------------------------
    const double zd0 = disk[hh].cfg->zd;
    const double zd2 = diskDimmingHeight * zd0 * diskDimmingHeight * zd0;
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NR; ii++){
      //-------------------------------------------------------------------
      const double R2 = RR[ii] * RR[ii];
      //-------------------------------------------------------------------
      /* get potential @ the reference points (mid plane and dimming scale) */
      const double PsiMid = Phi_spherical(RR[ii], sph, invlogrbin_sph);
      const double PsiDim = Phi_spherical(sqrt(R2 + zd2), sph, invlogrbin_sph);
      //-------------------------------------------------------------------
      const double Psi = PsiMid + diskDimmingHeightInv * (PsiDim - PsiMid);
      //-------------------------------------------------------------------
      double ratio;
      const int irr = findIdxSphericalPsi(Psi, sph, &ratio);
      const double rr = (1.0 - ratio) * sph[irr].rad + ratio * sph[1 + irr].rad;
      //-------------------------------------------------------------------
      /* fprintf(stderr, "%e\t%e\t%e\t%e\t%e\t%e\t%d\t%e\n", RR[ii], sqrt(rr * rr - R2), PsiMid, Psi, PsiDim, rr, irr, ratio); */
      const double zd1 = sqrt(rr * rr - R2);
      zd[ii] = (zd0 < zd1) ? zd0 : zd1;
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NR; ii++){ */
    //---------------------------------------------------------------------
  }/* for(int hh = 0; hh < ndisk; hh++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* static const double zd_scale_factor = 0.5; */
//-------------------------------------------------------------------------
/* static inline void initDiskThicknessControler(disk_data *disk) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   const double rs = disk->rs_spheroid; */
/*   //----------------------------------------------------------------------- */
/*   disk->zd_scale_min = 1.0; */
/*   disk->zd_scale_max = 1.0; */
/*   disk->zd_scale_R0  = 1.0; */
/*   disk->zd_scale_inv = 1.0; */
/*   //----------------------------------------------------------------------- */
/*   if( disk->cfg->zd > (zd_scale_factor * rs) ){ */
/*     //--------------------------------------------------------------------- */
/*     disk->zd_scale_min = zd_scale_factor * rs; */
/*     disk->zd_scale_max = disk->cfg->zd; */
/*     //--------------------------------------------------------------------- */
/*     disk->zd_scale_R0 = 0.25 * (rs + 2.0 * disk->cfg->zd); */
/*     disk->zd_scale_inv = 10.0 / (2.0 * disk->cfg->zd - rs); */
/*     //--------------------------------------------------------------------- */
/*   }/\* if( disk->cfg->zd > (zd_scale_factor * rs) ){ *\/ */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
/* static inline double modulateDiskThickness(const double RR, disk_data disk){ */
/*   //----------------------------------------------------------------------- */
/*   return (disk.zd_scale_min + 0.5 * (disk.zd_scale_max - disk.zd_scale_min) * (1.0 + tanh((RR - disk.zd_scale_R0) * disk.zd_scale_inv))); */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
/* void selectDominantComponent(const int ndisk, disk_data * restrict disk, const int skind, profile * restrict * prf, const double invlogbin, profile_cfg * restrict cfg) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "start"); */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   for(int ii = 0; ii < ndisk; ii++){ */
/*     //--------------------------------------------------------------------- */
/*     /\* pick up the dominant spherical component determines to determine gravitational force at r = scale height of the given disk component *\/ */
/*     //--------------------------------------------------------------------- */
/*     double ratio; */
/*     const int idx = findIdxSpherical(disk[ii].cfg->zd, prf[0], invlogbin, &ratio); */
/*     //--------------------------------------------------------------------- */
/*     int idxMax = 0; */
/*     double encMax = 0.0; */
/*     for(int jj = 0; jj < skind; jj++){ */
/*       //------------------------------------------------------------------- */
/*       const double enc = (1.0 - ratio) * prf[jj][idx].enc + ratio * prf[jj][1 + idx].enc; */
/*       if( enc > encMax ){	encMax = enc;	idxMax = jj;      } */
/*       //------------------------------------------------------------------- */
/*     }/\* for(int jj = 0; jj < skind; jj++){ *\/ */
/*     //--------------------------------------------------------------------- */
/*     // */
/*     //--------------------------------------------------------------------- */
/*     disk[ii].rs_spheroid = cfg[idxMax].rs; */
/*     //--------------------------------------------------------------------- */
/*   }/\* for(int ii = 0; ii < ndisk; ii++){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "end"); */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_GD_FORM_POTENTIAL
//-------------------------------------------------------------------------
static inline double func4GDpot(const double RR, const double a2, const disk_data disk)
{
  //-----------------------------------------------------------------------
  /* return (RR * disk.getColumnDensity(RR, disk.invRd, disk.sersic_ninv, disk.sersic_b, disk.Rcutoff, disk.invRsmooth) / sqrt(RR * RR - a2)); */
  return (RR * disk.getColumnDensity(RR, disk.invRd, disk.util) / sqrt(RR * RR - a2));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
double gaussQD4GDpot(const double min, const double max, const double a2, const disk_data disk)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  //-----------------------------------------------------------------------
  if( NINTBIN & 1 )
    sum = gsl_gaussQD_weight[(NINTBIN >> 1)] * func4GDpot(pls + mns * gsl_gaussQD_pos[(NINTBIN >> 1)], a2, disk);
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (func4GDpot(pls + mns * gsl_gaussQD_pos[ii], a2, disk) + func4GDpot(pls - mns * gsl_gaussQD_pos[ii], a2, disk));
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void calcGDpot
(const double zz, const int Nz, double * restrict ver, const double invdz, const double zp,
 const int ia, const double fa, const int iR, const double fR, const double aa, const double apR2, const double amR2,
 const int ndisk, disk_data * restrict disk, double * restrict invSigma, double * restrict ret)
{
  //-----------------------------------------------------------------------
  const double zpls2 = (zz - zp) * (zz - zp);  const double asinp = asin(2.0 * aa / (sqrt(zpls2 + apR2) + sqrt(zpls2 + amR2)));
  const double zmns2 = (zz + zp) * (zz + zp);  const double asinm = asin(2.0 * aa / (sqrt(zmns2 + apR2) + sqrt(zmns2 + amR2)));
  //-----------------------------------------------------------------------
#if 0
  int status = fpclassify(asinp);
  if( (status != FP_NORMAL) && (status != FP_ZERO) ){
    static char msg[64];
    switch( status ){
    case FP_NAN      :	sprintf(msg, "Not a Number"                                    );	break;
    case FP_INFINITE :	sprintf(msg, "either positive infinity or negative inifinity"  );	break;
    case FP_SUBNORMAL:	sprintf(msg, "too small to be represented in normalized format");	break;
    }/* switch( fpclassify(errMax) ){ */
    __KILL__(stderr, "ERROR: asinp is \"%s\".\n", msg);
  }/* if( fpclassify(errMax) != FP_NORMAL ){ */
  status = fpclassify(asinm);
  if( (status != FP_NORMAL) && (status != FP_ZERO) ){
    static char msg[64];
    switch( status ){
    case FP_NAN      :	sprintf(msg, "Not a Number"                                    );	break;
    case FP_INFINITE :	sprintf(msg, "either positive infinity or negative inifinity"  );	break;
    case FP_SUBNORMAL:	sprintf(msg, "too small to be represented in normalized format");	break;
    }/* switch( fpclassify(errMax) ){ */
    __KILL__(stderr, "ERROR: asinm is \"%s\".\n", msg);
  }/* if( fpclassify(errMax) != FP_NORMAL ){ */
#endif
  //-----------------------------------------------------------------------
  double fz;
  const int jz = bisection(zp, Nz, ver, false, invdz, &fz);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    const double zeta = invSigma[ii] *
      (((1.0 - fz) * (*disk[ii].rho)[INDEX2D(NR, Nz,     iR, jz)] + fz * (*disk[ii].rho)[INDEX2D(NR, Nz,     iR, 1 + jz)]) * (1.0 - fR) +
       ((1.0 - fz) * (*disk[ii].rho)[INDEX2D(NR, Nz, 1 + iR, jz)] + fz * (*disk[ii].rho)[INDEX2D(NR, Nz, 1 + iR, 1 + jz)]) *        fR   );
#if 0
    int status = fpclassify(zeta);
    if( (status != FP_NORMAL) && (status != FP_ZERO) ){
      static char msg[64];
      switch( status ){
      case FP_NAN      :	sprintf(msg, "Not a Number"                                    );	break;
      case FP_INFINITE :	sprintf(msg, "either positive infinity or negative inifinity"  );	break;
      case FP_SUBNORMAL:	sprintf(msg, "too small to be represented in normalized format");	break;
      }/* switch( fpclassify(errMax) ){ */
      __KILL__(stderr, "ERROR: zeta is \"%s\"; invSigma[%d] = %e, fR = %e, iR = %d, fz = %e, jz = %d; rho = %e, %e, %e, %e while DBL_MIN = %e.\n",
	       msg, ii, invSigma[ii], fR, iR, fz, jz,
	       (*disk[ii].rho)[INDEX2D(NR, Nz,     iR, jz)], (*disk[ii].rho)[INDEX2D(NR, Nz,     iR, 1 + jz)],
	       (*disk[ii].rho)[INDEX2D(NR, Nz, 1 + iR, jz)], (*disk[ii].rho)[INDEX2D(NR, Nz, 1 + iR, 1 + jz)], DBL_MIN);
    }/* if( fpclassify(errMax) != FP_NORMAL ){ */
#endif
    //---------------------------------------------------------------------
    const double diff = (1.0 - fa) * disk[ii].prf[ia].psi + fa * disk[ii].prf[1 + ia].psi;
    //---------------------------------------------------------------------
#if 0
    fprintf(stdout, "%d\t%e\t%e\t%e\t%e\t%e\n", iR, zp, zeta, diff, asinp, asinm);
#endif
    //---------------------------------------------------------------------
#if 0
    fprintf(stdout, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", sqrt(apR2) - aa, zz, aa, zp, asinp, asinm, 2.0 * aa / (sqrt(zpls2 + apR2) + sqrt(zpls2 + amR2)), 2.0 * aa / (sqrt(zmns2 + apR2) + sqrt(zmns2 + amR2)));
#endif
    //---------------------------------------------------------------------
    ret[ii] = zeta * diff * (asinp + asinm);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void _gaussQuad2d4GDpot
(const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 const int ia, const double fa, const int iR, const double fR, const double aa, const double apR2, const double amR2,
 const int ndisk, disk_data * restrict disk, double * restrict sum, double * invSigma, double * restrict tmp)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (zmax - zmin);
  const double pls = 0.5 * (zmax + zmin);
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] = 0.0;
  if( NINTBIN & 1 ){
    //---------------------------------------------------------------------
    const double ww =             gsl_gaussQD_weight[(NINTBIN >> 1)];
    const double zp = pls + mns * gsl_gaussQD_pos   [(NINTBIN >> 1)];
    calcGDpot(zz, Nz, ver, invdz, zp, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, invSigma, tmp);
    for(int ii = 0; ii < ndisk; ii++)      sum[ii] = ww * tmp[ii];
    //---------------------------------------------------------------------
  }/* if( NINTBIN & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double ww = gsl_gaussQD_weight[ii];
    //---------------------------------------------------------------------
    double zp = pls + mns * gsl_gaussQD_pos[ii];
    calcGDpot(zz, Nz, ver, invdz, zp, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, invSigma, tmp);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
    zp = pls - mns * gsl_gaussQD_pos[ii];
    calcGDpot(zz, Nz, ver, invdz, zp, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, invSigma, tmp);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
    //---------------------------------------------------------------------
  }/* for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  /* finalization */
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] *= mns;
  //-----------------------------------------------------------------------
#if 0
  fprintf(stdout, "%d\t%e\t%e\t%e\n", iR, zmin, zmax, sum[0]);
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline int findIdxSpherical4GDpot(const double rad, profile *prf, double *ratio)
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = NRADBIN - 1;
  //-----------------------------------------------------------------------
  if( rad < prf[ll].rho + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( rad > prf[rr].rho - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)ll + (uint)rr) >> 1;
    //---------------------------------------------------------------------
    if( (prf[cc].rho - rad) * (prf[ll].rho - rad) <= 0.0 )      rr = (int)cc;
    else                                                        ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr ){
      *ratio = (rad - prf[ll].rho) / (prf[rr].rho - prf[ll].rho);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void _setRdep4GDpot
(const double aa, int * restrict ia, double * restrict fa, const int ndisk, disk_data * restrict disk, double * restrict invSigma,
 const double RR, int * restrict iR, double * restrict fR, const int NR, double * restrict hor, const double invdR, double * restrict apR2, double * restrict amR2)
{
  //-----------------------------------------------------------------------
  *apR2 = (aa + RR) * (aa + RR);
  *amR2 = (aa - RR) * (aa - RR);
  //-----------------------------------------------------------------------
  *iR = bisection(aa, NR, hor, false, invdR, fR);
  *ia = findIdxSpherical4GDpot(aa, disk[0].prf, fa);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++)
    invSigma[ii] = 1.0 / (DBL_MIN + (1.0 - (*fR)) * disk[ii].Sigma[*iR] + (*fR) * disk[ii].Sigma[1 + (*iR)]);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void gaussQuad2d4GDpot
(const double RR, const int NR, double * restrict hor, const double invdR, const double Rmin, const double Rmax,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 const int ndisk, disk_data * restrict disk, double * restrict invSigma, double * restrict sub, double * restrict tmp, double * restrict sum)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (Rmax - Rmin);
  const double pls = 0.5 * (Rmax + Rmin);
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] = 0.0;
  if( NINTBIN & 1 ){
    //---------------------------------------------------------------------
    const double ww =             gsl_gaussQD_weight[(NINTBIN >> 1)];
    const double aa = pls + mns * gsl_gaussQD_pos   [(NINTBIN >> 1)];
    int ia, iR;
    double fa, fR, apR2, amR2;
    _setRdep4GDpot(aa, &ia, &fa, ndisk, disk, invSigma, RR, &iR, &fR, NR, hor, invdR, &apR2, &amR2);
    _gaussQuad2d4GDpot(zz, Nz, ver, invdz, zmin, zmax, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, tmp, invSigma, sub);
    for(int ii = 0; ii < ndisk; ii++)      sum[ii] = ww * tmp[ii];
    //---------------------------------------------------------------------
  }/* if( NINTBIN & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double ww = gsl_gaussQD_weight[ii];
    //---------------------------------------------------------------------
    double aa = pls + mns * gsl_gaussQD_pos[ii];
    int ia, iR;
    double fa, fR, apR2, amR2;
    _setRdep4GDpot(aa, &ia, &fa, ndisk, disk, invSigma, RR, &iR, &fR, NR, hor, invdR, &apR2, &amR2);
    _gaussQuad2d4GDpot(zz, Nz, ver, invdz, zmin, zmax, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, tmp, invSigma, sub);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
    //---------------------------------------------------------------------
    aa = pls - mns * gsl_gaussQD_pos[ii];
    _setRdep4GDpot(aa, &ia, &fa, ndisk, disk, invSigma, RR, &iR, &fR, NR, hor, invdR, &apR2, &amR2);
    _gaussQuad2d4GDpot(zz, Nz, ver, invdz, zmin, zmax, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, tmp, invSigma, sub);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
    //---------------------------------------------------------------------
  }/* for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  /* finalization */
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] *= mns;
  //-----------------------------------------------------------------------
#if 0
  fprintf(stdout, "%e\t%e\t%e\t%e\t%e\n", Rmin, Rmax, zmin, zmax, sum[0]);
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double integrateGDpot
(const double RR, const int NR, double * restrict hor, const double invdR, const double Rs, const double Rmax,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zs, const double zmax,
 const int ndisk, disk_data * restrict disk, double * restrict invSigma, double * restrict sub, double * restrict tmp, double * restrict sum)
{
  //-----------------------------------------------------------------------
  double Phi = 0.0;
  //-----------------------------------------------------------------------
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, 0.0     , 2.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, 2.0 * zs, 6.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, 6.0 * zs,       zmax, ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  //-----------------------------------------------------------------------
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 2.0 * Rs, 8.0 * Rs, zz, Nz, ver, invdz, 0.0     , 2.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 2.0 * Rs, 8.0 * Rs, zz, Nz, ver, invdz, 2.0 * zs, 6.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 2.0 * Rs, 8.0 * Rs, zz, Nz, ver, invdz, 6.0 * zs,       zmax, ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  //-----------------------------------------------------------------------
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 8.0 * Rs, Rmax, zz, Nz, ver, invdz, 0.0     , 2.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 8.0 * Rs, Rmax, zz, Nz, ver, invdz, 2.0 * zs, 6.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 8.0 * Rs, Rmax, zz, Nz, ver, invdz, 6.0 * zs,       zmax, ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  //-----------------------------------------------------------------------
  /* exit(0); */
  return (Phi);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_GD_FORM_POTENTIAL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef USE_ELLIPTIC_INTEGRAL
//-------------------------------------------------------------------------
static inline double rinv4pot(double d2, double dd, double phi){  return (1.0 / sqrt(DBL_MIN + d2 + dd * (1.0 - cos(phi))));}
/* static inline double rinv4pot(double d2, double dd, double phi){  return (1.0 / sqrt(d2 + dd * (1.0 - cos(phi))));} */
//-------------------------------------------------------------------------
double gaussQD4rinv4pot(const double min, const double max, const double r2, const double Rmean)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  //-----------------------------------------------------------------------
  if( NINTBIN & 1 )
    sum = gsl_gaussQD_weight[(NINTBIN >> 1)] * rinv4pot(r2, Rmean, pls + mns * gsl_gaussQD_pos[(NINTBIN >> 1)]);
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (rinv4pot(r2, Rmean, pls + mns * gsl_gaussQD_pos[ii]) + rinv4pot(r2, Rmean, pls - mns * gsl_gaussQD_pos[ii]));
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_ELLIPTIC_INTEGRAL
//-------------------------------------------------------------------------
static inline double calcPot
(const int iR, const double fR, const double Rp, const double R2, const double Rmean,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zp,
 double * restrict rho)
{
  //-----------------------------------------------------------------------
  const double zdiff = zz - zp;
#ifdef  USE_ELLIPTIC_INTEGRAL
  const double rinv = 1.0 / sqrt(R2 + zdiff * zdiff);
#endif//USE_ELLIPTIC_INTEGRAL
  //-----------------------------------------------------------------------
  double fz;
  const int jz = bisection(fabs(zp), Nz, ver, false, invdz, &fz);
  //-----------------------------------------------------------------------
  const double ret = Rp *
    (((1.0 - fz) * rho[INDEX2D(NR, Nz,     iR, jz)] + fz * rho[INDEX2D(NR, Nz,     iR, 1 + jz)]) * (1.0 - fR) +
     ((1.0 - fz) * rho[INDEX2D(NR, Nz, 1 + iR, jz)] + fz * rho[INDEX2D(NR, Nz, 1 + iR, 1 + jz)]) *        fR   )
#ifdef  USE_ELLIPTIC_INTEGRAL
    * 2.0 * rinv * gsl_sf_ellint_Kcomp(Rmean * rinv, GSL_PREC_DOUBLE)
#else///USE_ELLIPTIC_INTEGRAL
    * gaussQD4rinv4pot(0.0, M_PI, R2 + zdiff * zdiff, Rmean)
#endif//USE_ELLIPTIC_INTEGRAL
    ;
  //-----------------------------------------------------------------------
  return (ret);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
double _gaussQuad2d4calcPot
(const int iR, const double fR, const double Rp, const double R2, const double Rmean,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 double * restrict rho)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (zmax - zmin);
  const double pls = 0.5 * (zmax + zmin);
  double sum = 0.0;
  if( NINTBIN & 1 ){
    //---------------------------------------------------------------------
    const double weight =             gsl_gaussQD_weight[(NINTBIN >> 1)];
    const double     zp = pls + mns * gsl_gaussQD_pos   [(NINTBIN >> 1)];
    sum = weight * calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, zp, rho);
    //---------------------------------------------------------------------
  }/* if( NINTBIN & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double weight = gsl_gaussQD_weight[ii];
    //---------------------------------------------------------------------
    sum += weight * (calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, pls + mns * gsl_gaussQD_pos[ii], rho) +
		     calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, pls - mns * gsl_gaussQD_pos[ii], rho));
    //---------------------------------------------------------------------
  }/* for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void _setRdep4calcPot
(const double Rp, double * restrict R2, int * restrict iR, double * restrict fR, double * restrict Rmean,
 const double RR, const int NR, double * restrict hor, const double invdR)
{
  //-----------------------------------------------------------------------
#ifdef  USE_ELLIPTIC_INTEGRAL
  *Rmean = 2.0 * sqrt(RR * Rp);
  const double Rsum = RR + Rp;  *R2 = Rsum * Rsum;
#else///USE_ELLIPTIC_INTEGRAL
  *Rmean = 2.0 * RR * Rp;
  const double Rsub = RR - Rp;  *R2 = Rsub * Rsub;
#endif//USE_ELLIPTIC_INTEGRAL
  *iR = bisection(Rp, NR, hor, false, invdR, fR);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
double gaussQuad2d4calcPot
(const double RR, const int NR, double * restrict hor, const double invdR, const double Rmin, const double Rmax,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 double * restrict rho)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (Rmax - Rmin);
  const double pls = 0.5 * (Rmax + Rmin);
  double sum = 0.0;
  if( NINTBIN & 1 ){
    //---------------------------------------------------------------------
    const double weight =             gsl_gaussQD_weight[(NINTBIN >> 1)];
    const double     Rp = pls + mns * gsl_gaussQD_pos   [(NINTBIN >> 1)];
    double R2, fR, Rmean;
    int iR;
    _setRdep4calcPot(Rp, &R2, &iR, &fR, &Rmean, RR, NR, hor, invdR);
    sum = weight * _gaussQuad2d4calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, zmin, zmax, rho);
    //---------------------------------------------------------------------
  }/* if( NINTBIN & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double weight = gsl_gaussQD_weight[ii];
    //---------------------------------------------------------------------
    double Rp = pls + mns * gsl_gaussQD_pos[ii];
    double R2, fR, Rmean;
    int iR;
    _setRdep4calcPot(Rp, &R2, &iR, &fR, &Rmean, RR, NR, hor, invdR);
    sum += weight * _gaussQuad2d4calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, zmin, zmax, rho);
    //---------------------------------------------------------------------
    Rp = pls - mns * gsl_gaussQD_pos[ii];
    _setRdep4calcPot(Rp, &R2, &iR, &fR, &Rmean, RR, NR, hor, invdR);
    sum += weight * _gaussQuad2d4calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, zmin, zmax, rho);
    //---------------------------------------------------------------------
  }/* for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifndef USE_GD_FORM_POTENTIAL
static inline
#endif//USE_GD_FORM_POTENTIAL
double integratePotential
(const double RR, const int NR, double * restrict hor, const double invdR, const double Rs, const double Rmax,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zs, const double zmax, double * restrict rho)
{
  //-----------------------------------------------------------------------
  double Phi = 0.0;
  //-----------------------------------------------------------------------
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, -zmax, -zs  , rho);
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, -zs  ,  zs  , rho);
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz,  zs  ,  zmax, rho);
  //-----------------------------------------------------------------------
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 2.0 * Rs, Rmax, zz, Nz, ver, invdz, -zmax, -zs  , rho);
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 2.0 * Rs, Rmax, zz, Nz, ver, invdz, -zs  ,  zs  , rho);
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 2.0 * Rs, Rmax, zz, Nz, ver, invdz,  zs  ,  zmax, rho);
  //-----------------------------------------------------------------------
  return (Phi);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void calcOuterPotential
(const int NR, double * restrict hor, const double dR, const double invdR, const double Rmax,
 const int Nz, double * restrict ver, const double dz, const double invdz, const double zmax,
 const double Rs, const double zs, double * restrict Phi_NR, double * restrict Phi_Nz,
#ifdef  USE_GD_FORM_POTENTIAL
 const int ndisk, disk_data * restrict disk, double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum
#else///USE_GD_FORM_POTENTIAL
 double * restrict rho
#endif//USE_GD_FORM_POTENTIAL
#ifdef  CONFIRM_BUILDING_BLOCK
 , const double Mdisk
#endif//CONFIRM_BUILDING_BLOCK
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* \Phi_{  i, N_z} = potential @ R = (  i + 1/2) * R_max / N_R, z = (N_z + 1/2) * z_max / N_z; for i = 0, ..., N_R - 1 */
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < NR; ii++){
    for(int jj = 0; jj < Nz; jj++){
      fprintf(stderr, "%e\t%e", hor[ii], ver[jj]);
      for(int hh = 0; hh < ndisk; hh++)
	fprintf(stderr, "\t%e", (*disk[hh].rho)[INDEX2D(NR, Nz, ii, jj)]);
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }
  exit(0);
#endif
  //-----------------------------------------------------------------------
#ifdef  CONFIRM_BUILDING_BLOCK
  fprintf(stderr, "# boundary @ large z\n");
  fprintf(stderr, "# calc\tsphere\terror\n");
#else///CONFIRM_BUILDING_BLOCK
#pragma omp parallel for num_threads(CPUS_PER_PROCESS)
#endif//CONFIRM_BUILDING_BLOCK
  for(int jj = 0; jj < Nz; jj++){
    //---------------------------------------------------------------------
    const double RR = hor[NR - 1] + dR;
    const double zz = ver[jj];
    //---------------------------------------------------------------------
#ifdef  USE_GD_FORM_POTENTIAL
    const int tidx = omp_get_thread_num();
    double *inv = &stock_inv[ndisk * tidx];
    double *sub = &stock_sub[ndisk * tidx];
    double *tmp = &stock_tmp[ndisk * tidx];
    double *sum = &stock_sum[ndisk * tidx];
    Phi_NR[jj] = 4.0 * (double)newton * integrateGDpot(RR, NR, hor, invdR, Rs, Rmax, zz, Nz, ver, invdz, zs, zmax, ndisk, disk, inv, sub, tmp, sum);
#else///USE_GD_FORM_POTENTIAL
    Phi_NR[jj] = -2.0 * (double)newton * integratePotential(RR, NR, hor, invdR, Rs, Rmax, zz, Nz, ver, invdz, zs, zmax, rho);
#endif//USE_GD_FORM_POTENTIAL
    //---------------------------------------------------------------------
#ifdef  CONFIRM_BUILDING_BLOCK
    const double sphere = -(double)newton * Mdisk / sqrt(RR * RR + zz * zz);
    const double mndisk = -(double)newton * Mdisk / sqrt(RR * RR + (Rs + sqrt(zz * zz + zs * zs)) * (Rs + sqrt(zz * zz + zs * zs)));
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", Phi_NR[jj], sphere, 1.0 - Phi_NR[jj] / sphere, mndisk, 1.0 - Phi_NR[jj] / mndisk);
#endif//CONFIRM_BUILDING_BLOCK
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NR; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* \Phi_{N_R,   j} = potential @ R = (N_R + 1/2) * R_max / N_R, z = (  j + 1/2) * z_max / N_z; for j = 0, ..., N_z - 1 */
  //-----------------------------------------------------------------------
#ifdef  CONFIRM_BUILDING_BLOCK
  fprintf(stderr, "# boundary @ large R\n");
  fprintf(stderr, "# calc\tsphere\terror\n");
#else///CONFIRM_BUILDING_BLOCK
#pragma omp parallel for num_threads(CPUS_PER_PROCESS)
#endif//CONFIRM_BUILDING_BLOCK
  for(int ii = 0; ii < NR; ii++){
    //---------------------------------------------------------------------
    const double zz = ver[Nz - 1] + dz;
    const double RR = hor[ii];
    //---------------------------------------------------------------------
#ifdef  USE_GD_FORM_POTENTIAL
    const int tidx = omp_get_thread_num();
    double *inv = &stock_inv[ndisk * tidx];
    double *sub = &stock_sub[ndisk * tidx];
    double *tmp = &stock_tmp[ndisk * tidx];
    double *sum = &stock_sum[ndisk * tidx];
    Phi_Nz[ii] = 4.0 * (double)newton * integrateGDpot(RR, NR, hor, invdR, Rs, Rmax, zz, Nz, ver, invdz, zs, zmax, ndisk, disk, inv, sub, tmp, sum);
#else///USE_GD_FORM_POTENTIAL
    Phi_Nz[ii] = -2.0 * (double)newton * integratePotential(RR, NR, hor, invdR, Rs, Rmax, zz, Nz, ver, invdz, zs, zmax, rho);
#endif//USE_GD_FORM_POTENTIAL
    //---------------------------------------------------------------------
#ifdef  CONFIRM_BUILDING_BLOCK
    const double sphere = -(double)newton * Mdisk / sqrt(RR * RR + zz * zz);
    const double mndisk = -(double)newton * Mdisk / sqrt(RR * RR + (Rs + sqrt(zz * zz + zs * zs)) * (Rs + sqrt(zz * zz + zs * zs)));
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", Phi_Nz[ii], sphere, 1.0 - Phi_Nz[ii] / sphere, mndisk, 1.0 - Phi_Nz[ii] / mndisk);
#endif//CONFIRM_BUILDING_BLOCK
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#if 0
  bool failure = false;
  for(int ii = 0; ii < NR; ii++){
    int status = fpclassify(Phi_Nz[ii]);
    failure |= ((status != FP_NORMAL) && (status != FP_ZERO));
  }
  for(int jj = 0; jj < Nz; jj++){
    int status = fpclassify(Phi_NR[jj]);
    failure |= ((status != FP_NORMAL) && (status != FP_ZERO));
  }
  if( failure ){
    for(int ii = 0; ii < NR; ii++){
      for(int jj = 0; jj < Nz; jj++){
	fprintf(stderr, "%e\t%e", hor[ii], ver[jj]);
	for(int hh = 0; hh < ndisk; hh++)
	  fprintf(stderr, "\t%e", (*disk[hh].rho)[INDEX2D(NR, Nz, ii, jj)]);
	fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
    }

    /* for(int ii = 0; ii < NR; ii++){ */
    /*   fprintf(stderr, "%e\t%e", hor[ii], Phi_Nz[ii]); */
    /*   /\* for(int hh = 0; hh < ndisk; hh++) *\/ */
    /*   /\* 	fprintf(stderr, "\t%e", disk[hh].Sigma[ii]); *\/ */
    /*   fprintf(stderr, "\n"); */
    /* } */
    /* fprintf(stderr, "\n"); */
    /* for(int jj = 0; jj < Nz; jj++) */
    /*   fprintf(stderr, "%e\t%e\n", ver[jj], Phi_NR[jj]); */
    exit(0);
  }
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  CONFIRM_BUILDING_BLOCK
  fflush(NULL);
  exit(0);
#endif//CONFIRM_BUILDING_BLOCK
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void setSparseMatrix
(const int NR, const double dR, double * restrict RR, const int Nz, const double dz, double * restrict rho, double * restrict Phi_NR, double * restrict Phi_Nz,
 const crs mat, double * restrict vec)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* const int Nrow = NR * Nz; */
  const double dz2 = dz * dz;
  const double dR2 = dR * dR;
  int valIdx = 0;
  int rowIdx = 0;
  mat.row[rowIdx] = valIdx;
  //-----------------------------------------------------------------------
  const double Phi0 = 8.0 * M_PI * (double)newton * dR2 * dz2;
#pragma omp parallel for
  for(int ii = 0; ii < NR; ii++){
    //---------------------------------------------------------------------
    const double Ri = RR[ii];
    //---------------------------------------------------------------------
    for(int jj = 0; jj < Nz; jj++)
      vec[INDEX2D(NR, Nz, ii, jj)] = Phi0 * Ri * rho[INDEX2D(NR, Nz, ii, jj)];
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NR; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* inner boundary (ii = 0) */
  //-----------------------------------------------------------------------
  /* ii = 0; jj = 0; */
  mat.col[valIdx] = INDEX2D(NR, Nz, 0, 0);  mat.val[valIdx] = -dR * (dR2 + 2.0 * dz2);  valIdx++;
  mat.col[valIdx] = INDEX2D(NR, Nz, 0, 1);  mat.val[valIdx] =  dR *  dR2             ;  valIdx++;
  mat.col[valIdx] = INDEX2D(NR, Nz, 1, 0);  mat.val[valIdx] =  dR *        2.0 * dz2 ;  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
  //-----------------------------------------------------------------------
  /* ii = 0; 1 <= jj <= Nz - 2 */
  for(int jj = 1; jj < Nz - 1; jj++){
    //---------------------------------------------------------------------
    mat.col[valIdx] = INDEX2D(NR, Nz, 0, jj - 1);    mat.val[valIdx] =        dR *  dR2       ;    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, 0, jj    );    mat.val[valIdx] = -2.0 * dR * (dR2 + dz2);    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, 0, jj + 1);    mat.val[valIdx] =        dR *  dR2       ;    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, 1, jj    );    mat.val[valIdx] =  2.0 * dR *        dz2 ;    valIdx++;
    //---------------------------------------------------------------------
    rowIdx++;    mat.row[rowIdx] = valIdx;
    //---------------------------------------------------------------------
  }/* for(int jj = 1; jj < Nz - 1; jj++){ */
  //-----------------------------------------------------------------------
  /* ii = 0; jj = Nz - 1 */
  mat.col[valIdx] = INDEX2D(NR, Nz, 0, Nz - 2);  mat.val[valIdx] =        dR *  dR2       ;  valIdx++;
  mat.col[valIdx] = INDEX2D(NR, Nz, 0, Nz - 1);  mat.val[valIdx] = -2.0 * dR * (dR2 + dz2);  valIdx++;
  mat.col[valIdx] = INDEX2D(NR, Nz, 1, Nz - 1);  mat.val[valIdx] =  2.0 * dR *        dz2 ;  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
  vec[INDEX2D(NR, Nz, 0, Nz - 1)] -= dR * dR2 * Phi_Nz[0];
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* main domain (1 <= ii <= NR - 2) */
  //-----------------------------------------------------------------------
  for(int ii = 1; ii < NR - 1; ii++){
    //---------------------------------------------------------------------
    const double Ri = 2.0 * RR[ii];
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* inner boundary (jj = 0) */
    //---------------------------------------------------------------------
    mat.col[valIdx] = INDEX2D(NR, Nz, ii - 1, 0);    mat.val[valIdx] = (Ri - dR) *              dz2 ;    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, ii    , 0);    mat.val[valIdx] = -Ri       * (dR2 + 2.0 * dz2);    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, ii    , 1);    mat.val[valIdx] =  Ri       *  dR2             ;    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, ii + 1, 0);    mat.val[valIdx] = (Ri + dR) *              dz2 ;    valIdx++;
    //---------------------------------------------------------------------
    rowIdx++;    mat.row[rowIdx] = valIdx;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* main domain (1 <= jj <= Nz - 2) */
    //---------------------------------------------------------------------
    for(int jj = 1; jj < Nz - 1; jj++){
      //-------------------------------------------------------------------
      mat.col[valIdx] = INDEX2D(NR, Nz, ii - 1, jj    );      mat.val[valIdx] =       (Ri - dR) *        dz2 ;      valIdx++;
      mat.col[valIdx] = INDEX2D(NR, Nz, ii    , jj - 1);      mat.val[valIdx] =        Ri       *  dR2       ;      valIdx++;
      mat.col[valIdx] = INDEX2D(NR, Nz, ii    , jj    );      mat.val[valIdx] = -2.0 * Ri       * (dR2 + dz2);      valIdx++;
      mat.col[valIdx] = INDEX2D(NR, Nz, ii    , jj + 1);      mat.val[valIdx] =        Ri       *  dR2       ;      valIdx++;
      mat.col[valIdx] = INDEX2D(NR, Nz, ii + 1, jj    );      mat.val[valIdx] =       (Ri + dR) *        dz2 ;      valIdx++;
      //-------------------------------------------------------------------
      rowIdx++;      mat.row[rowIdx] = valIdx;
      //-------------------------------------------------------------------
    }/* for(int jj = 1; jj < Nz - 1; jj++){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* outer boundary (jj = Nz - 1) */
    //---------------------------------------------------------------------
    mat.col[valIdx] = INDEX2D(NR, Nz, ii - 1, Nz - 1);    mat.val[valIdx] =       (Ri - dR) *        dz2 ;    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, ii    , Nz - 2);    mat.val[valIdx] =        Ri       *  dR2       ;    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, ii    , Nz - 1);    mat.val[valIdx] = -2.0 * Ri       * (dR2 + dz2);    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, ii + 1, Nz - 1);    mat.val[valIdx] =       (Ri + dR) *        dz2 ;    valIdx++;
    //---------------------------------------------------------------------
    rowIdx++;    mat.row[rowIdx] = valIdx;
    //---------------------------------------------------------------------
    vec[INDEX2D(NR, Nz, ii, Nz - 1)] -= Ri * dR2 * Phi_Nz[ii];
    //---------------------------------------------------------------------
  }/* for(int ii = 1; ii < NR - 1; ii++){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* outer boundary (ii = NR - 1) */
  const double Rout = 2.0 * RR[NR - 1];
  //-----------------------------------------------------------------------
  /* ii = NR - 1; jj = 0; */
  mat.col[valIdx] = INDEX2D(NR, Nz, NR - 2, 0);  mat.val[valIdx] = (Rout - dR) *              dz2 ;  valIdx++;
  mat.col[valIdx] = INDEX2D(NR, Nz, NR - 1, 0);  mat.val[valIdx] = -Rout       * (dR2 + 2.0 * dz2);  valIdx++;
  mat.col[valIdx] = INDEX2D(NR, Nz, NR - 1, 1);  mat.val[valIdx] =  Rout       *  dR2             ;  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
  vec[INDEX2D(NR, Nz, NR - 1, 0)] -= (Rout + dR) * dz2 * Phi_NR[0];
  //-----------------------------------------------------------------------
  /* ii = NR - 1; 1 <= jj <= Nz - 2 */
  for(int jj = 1; jj < Nz - 1; jj++){
    //---------------------------------------------------------------------
    mat.col[valIdx] = INDEX2D(NR, Nz, NR - 2, jj    );    mat.val[valIdx] =       (Rout - dR) *        dz2 ;    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, NR - 1, jj - 1);    mat.val[valIdx] =        Rout       *  dR2       ;    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, NR - 1, jj    );    mat.val[valIdx] = -2.0 * Rout       * (dR2 + dz2);    valIdx++;
    mat.col[valIdx] = INDEX2D(NR, Nz, NR - 1, jj + 1);    mat.val[valIdx] =        Rout       *  dR2       ;    valIdx++;
    //---------------------------------------------------------------------
    rowIdx++;    mat.row[rowIdx] = valIdx;
    //---------------------------------------------------------------------
    vec[INDEX2D(NR, Nz, NR - 1, jj)] -= (Rout + dR) * dz2 * Phi_NR[jj];
    //---------------------------------------------------------------------
  }/* for(int jj = 1; jj < Nz - 1; jj++){ */
  //-----------------------------------------------------------------------
  /* ii = NR - 1; jj = Nz - 1 */
  mat.col[valIdx] = INDEX2D(NR, Nz, NR - 2, Nz - 1);  mat.val[valIdx] =       (Rout - dR) *        dz2 ;  valIdx++;
  mat.col[valIdx] = INDEX2D(NR, Nz, NR - 1, Nz - 2);  mat.val[valIdx] =        Rout       *  dR2       ;  valIdx++;
  mat.col[valIdx] = INDEX2D(NR, Nz, NR - 1, Nz - 1);  mat.val[valIdx] = -2.0 * Rout       * (dR2 + dz2);  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
  vec[INDEX2D(NR, Nz, NR - 1, Nz - 1)] -= (Rout + dR) * dz2 * Phi_NR[Nz - 1] + Rout * dR2 * Phi_Nz[NR - 1];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  if( mat.row[rowIdx] != (5 * NR * Nz - 2 * (NR + Nz)) ){
    //---------------------------------------------------------------------
    __KILL__(stderr, "ERROR: Number of non-zero elements is %d, while the expected is %d\n", mat.row[rowIdx], 5 * NR * Nz - 2 * (NR + Nz));
    //---------------------------------------------------------------------
  }/* if( mat.row[rowIdx] != (5 * NR * Nz - 2 * (NR + Nz)) ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* assume potential generated by a Miyamoto & Nagai disk for initial guess to the solution of the Poisson equation */
//-------------------------------------------------------------------------
void initPotentialField(const int NR, double * restrict RR, const int Nz, double * restrict zz, double * restrict Phi, const int ndisk, disk_data * restrict disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* Miyamoto & Nagai disk */
#pragma omp parallel for
  for(int ii = 0; ii < NR; ii++){
    //---------------------------------------------------------------------
    for(int jj = 0; jj < Nz; jj++)
      Phi[INDEX2D(NR, Nz, ii, jj)] = 0.0;
    //---------------------------------------------------------------------
    const double R2 = RR[ii] * RR[ii];
    //---------------------------------------------------------------------
    for(int jj = 0; jj < Nz; jj++){
      //-------------------------------------------------------------------
      const double z2 = zz[jj] * zz[jj];
      //-------------------------------------------------------------------
      for(int kk = 0; kk < ndisk; kk++){
	//-----------------------------------------------------------------
	const double Rs = disk[kk].cfg->rs + sqrt(z2 + disk[kk].cfg->zd * disk[kk].cfg->zd);
	Phi[INDEX2D(NR, Nz, ii, jj)] -= (double)newton * disk[kk].cfg->Mtot / sqrt(R2 + Rs * Rs);
	//-----------------------------------------------------------------
      }/* for(int kk = 0; kk < ndisk; kk++){ */
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < Nz; jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NR; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void getPotentialField
(const int ndisk, disk_data * restrict disk,
 const int NR, double * restrict RR, const double dR, const double invdR, const double Rmax,
 const int Nz, double * restrict zz, const double dz, const double invdz, const double zmax,
 double * restrict rho, double * restrict Phi, double * restrict Phi_NR, double * restrict Phi_Nz,
 const crs mat, const crs ilu, double * restrict vec, double * restrict res, double * restrict sdw, double * restrict mid, double * restrict tmp,
 double * restrict Api, double * restrict Ati, double * restrict Kri, double * restrict Kpi, double * restrict Kti
#ifdef  USE_GD_FORM_POTENTIAL
 , double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum
#endif//USE_GD_FORM_POTENTIAL
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* prepare Poisson equation in matrix form */
  //-----------------------------------------------------------------------
  double Mtot = 0.0;
  double Rs   = 0.0;
  double zs   = 0.0;
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    Mtot += disk[ii].cfg->Mtot;
    Rs   += disk[ii].cfg->Mtot * disk[ii].cfg->rs;
    zs   += disk[ii].cfg->Mtot * disk[ii].cfg->zd;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  const double Minv = 1.0 / Mtot;
  Rs *= Minv;
  zs *= Minv;
  //-----------------------------------------------------------------------
/* #ifdef  USE_GD_FORM_POTENTIAL */
/*   for(int ii = 0; ii < NR; ii++){ */
/*     for(int jj = 0; jj < Nz; jj++){ */
/*       fprintf(stderr, "%e\t%e", RR[ii], zz[jj]); */
/*       for(int hh = 0; hh < ndisk; hh++) */
/* 	fprintf(stderr, "\t%e", (*disk[hh].rho)[INDEX2D(NR, Nz, ii, jj)]); */
/*       fprintf(stderr, "\n"); */
/*     } */
/*     fprintf(stderr, "\n"); */
/*   } */
/*   exit(0); */
/* #endif//USE_GD_FORM_POTENTIAL */
  //-----------------------------------------------------------------------
#if 0
#ifdef  USE_GD_FORM_POTENTIAL
  static int printSteps = 0;
  if( printSteps == 1 ){
    for(int ii = 0; ii < NR; ii++){
      for(int jj = 0; jj < Nz; jj++){
	fprintf(stderr, "%e\t%e", RR[ii], zz[jj]);
	for(int hh = 0; hh < ndisk; hh++)
	  fprintf(stderr, "\t%e", (*disk[hh].rho)[INDEX2D(NR, Nz, ii, jj)]);
	fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
    }
    exit(0);
  }
  printSteps++;
#endif//USE_GD_FORM_POTENTIAL
#endif
  //-----------------------------------------------------------------------
  calcOuterPotential(NR, RR, dR, invdR, Rmax, Nz, zz, dz, invdz, zmax, Rs, zs, Phi_NR, Phi_Nz,
#ifdef  USE_GD_FORM_POTENTIAL
		     ndisk, disk, stock_inv, stock_sub, stock_tmp, stock_sum
#else///USE_GD_FORM_POTENTIAL
		     rho
#endif//USE_GD_FORM_POTENTIAL
#ifdef  CONFIRM_BUILDING_BLOCK
		     , Mtot
#endif//CONFIRM_BUILDING_BLOCK
		     );
  setSparseMatrix(NR, dR, RR, Nz, dz, rho, Phi_NR, Phi_Nz, mat, vec);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* solve Poisson equation using iterative method */
  //-----------------------------------------------------------------------
  const int Nrow = NR * Nz;
  getILU0(Nrow, mat, ilu);
  pbicgstab(mat, Nrow, vec, Phi, res, sdw, mid, tmp, Api, Ati, ilu, Kri, Kpi, Kti, CONVERGENCE_BICGSTAB);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* static inline double gaussQuad1d4Rho(double (*func)(double, double, double, double, double, double), const double min, const double max, const double xinv, const double ninv, const double bb, const double Rcutoff, const double invRsmooth) */
static inline double gaussQuad1d4Rho(double (*func)(double, double, disk_util), const double min, const double max, const double xinv, const disk_util disk)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  if( NINTBIN_LOW & 1 )
    sum = gsl_gaussQD_low_weight[(NINTBIN_LOW >> 1)] * func(pls + mns * gsl_gaussQD_low_pos[(NINTBIN_LOW >> 1)], xinv, disk);
    /* sum = gsl_gaussQD_low_weight[(NINTBIN_LOW >> 1)] * func(pls + mns * gsl_gaussQD_low_pos[(NINTBIN_LOW >> 1)], xinv, ninv, bb, Rcutoff, invRsmooth); */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_low_weight[ii] *
      (func(pls + mns * gsl_gaussQD_low_pos[ii], xinv, disk) +
       func(pls - mns * gsl_gaussQD_low_pos[ii], xinv, disk));
      /* (func(pls + mns * gsl_gaussQD_low_pos[ii], xinv, ninv, bb, Rcutoff, invRsmooth) + */
      /*  func(pls - mns * gsl_gaussQD_low_pos[ii], xinv, ninv, bb, Rcutoff, invRsmooth)); */
  //-----------------------------------------------------------------------
  /* finalization */
  return (sum * mns);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline double _setzdep4calcRho
(const double zz, const int Nz, double * restrict ver, const double invdz,
 double * restrict Phi, profile * restrict sph, const double invlogrbin_sph,
 const int iR, const double fR, const double R2, const double PsiR0, const double invPsi)
{
  //-----------------------------------------------------------------------
  double fz;
  const int jz = bisection(zz, Nz, ver, false, invdz, &fz);
  //-----------------------------------------------------------------------
  const double PsiRz =
    ((1.0 - fz) * (-Phi[INDEX2D(NR, Nz,     iR, jz)]) + fz * (-Phi[INDEX2D(NR, Nz,     iR, 1 + jz)])) * (1.0 - fR) +
    ((1.0 - fz) * (-Phi[INDEX2D(NR, Nz, 1 + iR, jz)]) + fz * (-Phi[INDEX2D(NR, Nz, 1 + iR, 1 + jz)])) *        fR  +
    Phi_spherical(sqrt(R2 + zz * zz), sph, invlogrbin_sph);
  //-----------------------------------------------------------------------
  return (exp(-(PsiRz - PsiR0) * invPsi));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void _setRdep4calcRho
(const double RR, double * restrict R2, double * restrict horDep, int * restrict iR, double * restrict fR, double * restrict PsiR0, double * restrict invPsi,
 const int NR, double * restrict horRad, const double invdR, const int Nz, double * restrict ver, const double invdz, double * restrict Phi,
 const disk_data disk, profile * restrict sph, const double invlogrbin_sph)
{
  //-----------------------------------------------------------------------
  *R2 = RR * RR;
  //-----------------------------------------------------------------------
  *iR = bisection(RR, NR, horRad, false, invdR, fR);
  *PsiR0 = (1.0 - (*fR)) * (-Phi[INDEX2D(NR, Nz, *iR, 0)]) + (*fR) * (-Phi[INDEX2D(NR, Nz, 1 + (*iR), 0)]) + Phi_spherical(RR, sph, invlogrbin_sph);
  //-----------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  /* const double value = modulateDiskThickness(RR, disk); */
  /* const double invzd = disk.invzd   / value; */
  /* const double    zd = disk.cfg->zd * value; */
  const double    zd = (1.0 - (*fR)) * disk.zd[*iR] + (*fR) * disk.zd[1 + (*iR)];
#else///ENABLE_VARIABLE_SCALE_HEIGHT
  const double    zd = disk.cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  const double invzd = 1.0 / zd;
  //-----------------------------------------------------------------------
  /* *horDep = invzd * disk.getColumnDensity(RR, disk.invRd, disk.sersic_ninv, disk.sersic_b, disk.Rcutoff, disk.invRsmooth); */
  *horDep = invzd * disk.getColumnDensity(RR, disk.invRd, disk.util);
  //-----------------------------------------------------------------------
  double fz;
  const int jz = bisection(zd, Nz, ver, false, invdz, &fz);
  const double Psi_Rzd =
    ((1.0 - fz) * (-Phi[INDEX2D(NR, Nz,      *iR , jz)]) + fz * (-Phi[INDEX2D(NR, Nz,      *iR , 1 + jz)])) * (1.0 - (*fR)) +
    ((1.0 - fz) * (-Phi[INDEX2D(NR, Nz, 1 + (*iR), jz)]) + fz * (-Phi[INDEX2D(NR, Nz, 1 + (*iR), 1 + jz)])) *        (*fR)  +
    Phi_spherical(sqrt((*R2) + zd * zd), sph, invlogrbin_sph);
  *invPsi = 1.0 / (Psi_Rzd - (*PsiR0));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double _gaussQuad2d4calcRho
(const double zmin, const double zmax, const int Nz, double * restrict ver, const double invdz, const int iR, const double fR, const double R2,
 const double PsiR0, const double invPsi, double * restrict Phi, profile * restrict sph, const double invlogrbin_sph)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (zmax - zmin);
  const double pls = 0.5 * (zmax + zmin);
  double sum = 0.0;
  if( NINTBIN_LOW & 1 ){
    //---------------------------------------------------------------------
    sum = gsl_gaussQD_low_weight[(NINTBIN_LOW >> 1)] *
      _setzdep4calcRho(pls + mns * gsl_gaussQD_low_pos[(NINTBIN_LOW >> 1)], Nz, ver, invdz, Phi, sph, invlogrbin_sph, iR, fR, R2, PsiR0, invPsi);
    //---------------------------------------------------------------------
  }/* if( NINTBIN_LOW & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_low_weight[ii] *
      (_setzdep4calcRho(pls + mns * gsl_gaussQD_low_pos[ii], Nz, ver, invdz, Phi, sph, invlogrbin_sph, iR, fR, R2, PsiR0, invPsi) +
       _setzdep4calcRho(pls - mns * gsl_gaussQD_low_pos[ii], Nz, ver, invdz, Phi, sph, invlogrbin_sph, iR, fR, R2, PsiR0, invPsi));
  //-----------------------------------------------------------------------
  /* finalization */
  return (sum * mns);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
double gaussQuad2d4calcRho
(const double Rmin, const double Rmax, const int NR, double * restrict hor, const double invdR,
 const double zmin, const double zmax, const int Nz, double * restrict ver, const double invdz,
 profile * restrict sph, const double invlogrbin_sph, double * restrict Phi, const disk_data disk)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (Rmax - Rmin);
  const double pls = 0.5 * (Rmax + Rmin);
  double sum = 0.0;
  if( NINTBIN_LOW & 1 ){
    //---------------------------------------------------------------------
    double R2, radVal, fR, PsiR0, invPsi;
    int iR;
    _setRdep4calcRho(pls + mns * gsl_gaussQD_low_pos[(NINTBIN_LOW >> 1)], &R2, &radVal, &iR, &fR, &PsiR0, &invPsi, NR, hor, invdR, Nz, ver, invdz, Phi, disk, sph, invlogrbin_sph);
    //---------------------------------------------------------------------
    sum = gsl_gaussQD_low_weight[(NINTBIN_LOW >> 1)] * radVal * _gaussQuad2d4calcRho(zmin, zmax, Nz, ver, invdz, iR, fR, R2, PsiR0, invPsi, Phi, sph, invlogrbin_sph);
    //---------------------------------------------------------------------
  }/* if( NINTBIN_LOW & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    double R2, radVal, fR, PsiR0, invPsi;
    int iR;
    //---------------------------------------------------------------------
    _setRdep4calcRho(pls + mns * gsl_gaussQD_low_pos[ii], &R2, &radVal, &iR, &fR, &PsiR0, &invPsi, NR, hor, invdR, Nz, ver, invdz, Phi, disk, sph, invlogrbin_sph);
    sum += gsl_gaussQD_low_weight[ii] * radVal * _gaussQuad2d4calcRho(zmin, zmax, Nz, ver, invdz, iR, fR, R2, PsiR0, invPsi, Phi, sph, invlogrbin_sph);
    //---------------------------------------------------------------------
    _setRdep4calcRho(pls - mns * gsl_gaussQD_low_pos[ii], &R2, &radVal, &iR, &fR, &PsiR0, &invPsi, NR, hor, invdR, Nz, ver, invdz, Phi, disk, sph, invlogrbin_sph);
    sum += gsl_gaussQD_low_weight[ii] * radVal * _gaussQuad2d4calcRho(zmin, zmax, Nz, ver, invdz, iR, fR, R2, PsiR0, invPsi, Phi, sph, invlogrbin_sph);
    //---------------------------------------------------------------------
  }/* for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  /* finalization */
  return (sum * mns);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* static inline double func4encSigma(const double RR, const disk_data disk){  return (RR * disk.getColumnDensity(RR, disk.invRd, disk.sersic_ninv, disk.sersic_b, disk.Rcutoff, disk.invRsmooth));} */
static inline double func4encSigma(const double RR, const disk_data disk){  return (RR * disk.getColumnDensity(RR, disk.invRd, disk.util));}
//-------------------------------------------------------------------------
static inline double gaussQD4encSigma(const double min, const double max, const disk_data disk)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  //-----------------------------------------------------------------------
  if( NINTBIN & 1 )
    sum = gsl_gaussQD_weight[(NINTBIN >> 1)] * func4encSigma(pls + mns * gsl_gaussQD_pos[(NINTBIN >> 1)], disk);
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (func4encSigma(pls + mns * gsl_gaussQD_pos[ii], disk) + func4encSigma(pls - mns * gsl_gaussQD_pos[ii], disk));
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void swapDblArrays(double **p0, double **p1)
{
  //-----------------------------------------------------------------------
  double *tmp;
  //-----------------------------------------------------------------------
  tmp = *p0;
  *p0 = *p1;
  *p1 = tmp;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void getPotDensPair
(const int ndisk, disk_data * restrict disk, const int NR, const double Rmax, const int Nz, const double zmax,
 double * restrict Phi_NR, double * restrict Phi_Nz,
 const crs mat, const crs ilu, double * restrict vec, double * restrict res, double * restrict sdw, double * restrict mid, double * restrict tmp,
 double * restrict Api, double * restrict Ati, double * restrict Kri, double * restrict Kpi, double * restrict Kti
#ifdef  USE_GD_FORM_POTENTIAL
 , double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum
#endif//USE_GD_FORM_POTENTIAL
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* load disk data (common settings for all components) */
  //-----------------------------------------------------------------------
  double * RR;  RR  = disk[0].hor;
  double * zz;  zz  = disk[0].ver;
  double *Phi;  Phi = disk[0].pot;
  //-----------------------------------------------------------------------
  profile *sph;  sph = disk[0].prf;
  const double invlogrbin_sph = disk[0].invlogrbin;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set grid points */
  //-----------------------------------------------------------------------
  /* assume regular grid in each direction */
  const double Rbin = Rmax / (double)NR;  const double invRbin = 1.0 / Rbin;
  const double zbin = zmax / (double)Nz;  const double invzbin = 1.0 / zbin;
  for(int ii = 0; ii < NR; ii++)    RR[ii] = (0.5 + (double)ii) * Rbin;
  for(int jj = 0; jj < Nz; jj++)    zz[jj] = (0.5 + (double)jj) * zbin;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    disk[ii].Rbin = Rbin;    disk[ii].invRbin = invRbin;
    disk[ii].zbin = zbin;    disk[ii].invzbin = invzbin;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  //-----------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  getVariableDiskScaleHeight(ndisk, disk, NR, sph, invlogrbin_sph);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set column density profile on the midplane */
  for(int hh = 0; hh < ndisk; hh++){
    //---------------------------------------------------------------------
    /* set column density profile */
#pragma omp parallel for
#ifndef PROHIBIT_EXTRAPOLATION
    for(int ii = 0; ii < NR; ii++)
#else///PROHIBIT_EXTRAPOLATION
    for(int ii = 0; ii < NR - 1; ii++)
#endif//PROHIBIT_EXTRAPOLATION
      disk[hh].Sigma[ii] = disk[hh].cfg->Sigma0 * gaussQuad1d4Rho(disk[hh].getColumnDensity, RR[ii] - 0.5 * Rbin, RR[ii] + 0.5 * Rbin, disk[hh].invRd, disk[hh].util) * invRbin;
      /* disk[hh].Sigma[ii] = disk[hh].cfg->Sigma0 * gaussQuad1d4Rho(disk[hh].getColumnDensity, RR[ii] - 0.5 * Rbin, RR[ii] + 0.5 * Rbin, disk[hh].invRd, disk[hh].sersic_ninv, disk[hh].sersic_b, disk[hh].Rcutoff, disk[hh].invRsmooth) * invRbin; */
#ifdef  PROHIBIT_EXTRAPOLATION
    disk[hh].Sigma[NR - 1] = 0.0;
#endif//PROHIBIT_EXTRAPOLATION
    //---------------------------------------------------------------------
  }/* for(int hh = 0; hh < ndisk; hh++){ */
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < NR; ii++){
    fprintf(stderr, "%e\t%e", RR[ii], disk[0].Sigma[ii]);
    for(int hh = 1; hh < ndisk; hh++)
      fprintf(stderr, "\t%e", disk[hh].Sigma[ii]);
    fprintf(stderr, "\n");
  }
  exit(0);
#endif
  //-----------------------------------------------------------------------
  /* set density field */
  /* initial vertical profile is assumed to be hyperbolic secant */
  static bool init = true;
  if( init ){
    //---------------------------------------------------------------------
    init = false;
    setz0inv4SechProfile();
    //---------------------------------------------------------------------
    for(int hh = 0; hh < ndisk; hh++){
#pragma omp parallel for
#ifndef PROHIBIT_EXTRAPOLATION
      for(int ii = 0; ii < NR; ii++)
#else///PROHIBIT_EXTRAPOLATION
      for(int ii = 0; ii < NR - 1; ii++)
#endif//PROHIBIT_EXTRAPOLATION
	{
	  //---------------------------------------------------------------
	  /* set scale height @ given R */
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
	  /* const double invzd = disk[hh].invzd / modulateDiskThickness(RR[ii], disk[hh]); */
	  const double invzd = 1.0 / disk[hh].zd[ii];
#else///ENABLE_VARIABLE_SCALE_HEIGHT
	  const double invzd = 1.0 / disk[hh].cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
	  const double rho0 = 0.5 * disk[hh].Sigma[ii] * invzd;
	  //---------------------------------------------------------------
#ifndef PROHIBIT_EXTRAPOLATION
	  for(int jj = 0; jj < Nz; jj++)
#else///PROHIBIT_EXTRAPOLATION
	  for(int jj = 0; jj < Nz - 1; jj++)
#endif//PROHIBIT_EXTRAPOLATION
	    (*disk[hh].rho)[INDEX2D(NR, Nz, ii, jj)] = rho0 * gaussQuad1d4Rho(getVerticalDensity, zz[jj] - 0.5 * zbin, zz[jj] + 0.5 * zbin, invzd, disk[hh].util);
	    /* (*disk[hh].rho)[INDEX2D(NR, Nz, ii, jj)] = rho0 * gaussQuad1d4Rho(getVerticalDensity, zz[jj] - 0.5 * zbin, zz[jj] + 0.5 * zbin, invzd, disk[hh].sersic_ninv, disk[hh].sersic_b, disk[hh].Rcutoff, disk[hh].invRsmooth); */
#ifdef  PROHIBIT_EXTRAPOLATION
	  (*disk[hh].rho)[INDEX2D(NR, Nz, ii, Nz - 1)] = 0.0;
#endif//PROHIBIT_EXTRAPOLATION
	  //---------------------------------------------------------------
	  /* calculate column density using Simpson's rule */
	  double Sigma = (*disk[hh].rho)[INDEX2D(NR, Nz, ii, 0)];
	  for(int jj = 1; jj < Nz - 1; jj++)
	    Sigma += (double)(1 << (1 + (jj & 1))) * (*disk[hh].rho)[INDEX2D(NR, Nz, ii, jj)];
	  Sigma += (*disk[hh].rho)[INDEX2D(NR, Nz, ii, Nz - 1)];
	  Sigma *= 2.0 * zbin / 3.0;/* 2.0 reflects plane symmetry about the equatorial plane (z = 0) */
	  /* calibrate column density */
	  const double Mscale = disk[hh].Sigma[ii] / (DBL_MIN + Sigma);
	  for(int jj = 0; jj < Nz; jj++)
	    (*disk[hh].rho)[INDEX2D(NR, Nz, ii, jj)] *= Mscale;
	  //---------------------------------------------------------------
	}
#ifdef  PROHIBIT_EXTRAPOLATION
    for(int jj = 0; jj < Nz; jj++)
      (*disk[hh].rho)[INDEX2D(NR, Nz, NR - 1, jj)] = 0.0;
#endif//PROHIBIT_EXTRAPOLATION
    }/* for(int hh = 0; hh < ndisk; hh++){ */
    //---------------------------------------------------------------------
    initPotentialField(NR, RR, Nz, zz, Phi, ndisk, disk);
    //---------------------------------------------------------------------
  }/* if( init ){ */
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < NR; ii++){
    for(int jj = 0; jj < Nz; jj++){
      fprintf(stderr, "%e\t%e", RR[ii], zz[jj]);
      for(int hh = 0; hh < ndisk; hh++)
	fprintf(stderr, "\t%e", (*disk[hh].rho)[INDEX2D(NR, Nz, ii, jj)]);
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }
  exit(0);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* iterative process to get the potential-density pair */
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# procedure start: NR = %d, Nz = %d\n", NR, Nz);
  fflush(stdout);
  int steps = 0;
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    /* calculate potential field from the given density field */
    //---------------------------------------------------------------------
    /* set density field */
    double *rhoTot;
    if( ndisk > 1 ){
      //-------------------------------------------------------------------
      rhoTot = disk[0].rhoTot;
      //-------------------------------------------------------------------
#pragma omp parallel for
      for(int ii = 0; ii < NR * Nz; ii++)
	rhoTot[ii] = 0.0;
      for(int hh = 0; hh < ndisk; hh++)
#pragma omp parallel for
	for(int ii = 0; ii < NR * Nz; ii++)
	  rhoTot[ii] += (*disk[hh].rho)[ii];
      //-------------------------------------------------------------------
    }/* if( ndisk > 1 ){ */
    else      rhoTot = *disk[0].rho;
    //---------------------------------------------------------------------
#if 0
    static int printSteps = 0;
    if( printSteps == 1 ){
      for(int ii = 0; ii < NR; ii++){
	for(int jj = 0; jj < Nz; jj++)
	  fprintf(stderr, "%e\t%e\t%e\n", RR[ii], zz[jj], rhoTot[INDEX2D(NR, Nz, ii, jj)]);
	fprintf(stderr, "\n");
      }
      exit(0);
    }
    printSteps++;
#endif
    //---------------------------------------------------------------------
    /* solve Poisson equation using ILU(0) preconditioned BiCGSTAB method */
    getPotentialField(ndisk, disk, NR, RR, Rbin, invRbin, Rmax, Nz, zz, zbin, invzbin, zmax,
		      rhoTot, Phi, Phi_NR, Phi_Nz, mat, ilu, vec, res, sdw, mid, tmp, Api, Ati, Kri, Kpi, Kti
#ifdef  USE_GD_FORM_POTENTIAL
		      , stock_inv, stock_sub, stock_tmp, stock_sum
#endif//USE_GD_FORM_POTENTIAL
		      );
    //---------------------------------------------------------------------
#if 0
    for(int ii = 0; ii < NR; ii++){
      for(int jj = 0; jj < Nz; jj++)
	fprintf(stderr, "%e\t%e\t%e\n", RR[ii], zz[jj], Phi[INDEX2D(NR, Nz, ii, jj)]);
      fprintf(stderr, "\n");
    }
    exit(0);
#endif
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* update density field from the derived potential field */
    //---------------------------------------------------------------------
    static double errMax;
    errMax = 0.0;
    //---------------------------------------------------------------------
    for(int hh = 0; hh < ndisk; hh++){
#pragma omp parallel for
#ifndef PROHIBIT_EXTRAPOLATION
      for(int ii = 0; ii < NR; ii++)
#else///PROHIBIT_EXTRAPOLATION
      for(int ii = 0; ii < NR - 1; ii++)
#endif//PROHIBIT_EXTRAPOLATION
	{
	  //---------------------------------------------------------------
	  const double Rmin = RR[ii] - 0.5 * Rbin;
	  const double Rmax = RR[ii] + 0.5 * Rbin;
	  /* const double rho0 = (*disk[hh].rho)[INDEX2D(NR, Nz, ii, 0)] / gaussQuad2d4calcRho(Rmin, Rmax, NR, RR, invRbin, 0.0, zbin, Nz, zz, invzbin, sph, invlogrbin_sph, Phi, disk[hh]); */
	  const double rho0 = (*disk[hh].rho)[INDEX2D(NR, Nz, ii, 0)] / (DBL_MIN + gaussQuad2d4calcRho(Rmin, Rmax, NR, RR, invRbin, 0.0, zbin, Nz, zz, invzbin, sph, invlogrbin_sph, Phi, disk[hh]));
#if 0
	  fprintf(stdout, "%d\t%e\t%e\n", hh, RR[ii], gaussQuad2d4calcRho(Rmin, Rmax, NR, RR, invRbin, 0.0, zbin, Nz, zz, invzbin, sph, invlogrbin_sph, Phi, disk[hh]));
#endif
	  //---------------------------------------------------------------
	  /* calculate vertical density profile */
	  (*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, 0)] = (*disk[hh].rho)[INDEX2D(NR, Nz, ii, 0)];
#ifndef PROHIBIT_EXTRAPOLATION
	  for(int jj = 1; jj < Nz; jj++)
#else///PROHIBIT_EXTRAPOLATION
	  for(int jj = 1; jj < Nz - 1; jj++)
#endif//PROHIBIT_EXTRAPOLATION
	    {
	      (*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, jj)] = rho0 * gaussQuad2d4calcRho
		(Rmin, Rmax, NR, RR, invRbin, zz[jj] - 0.5 * zbin, zz[jj] + 0.5 * zbin, Nz, zz, invzbin, sph, invlogrbin_sph, Phi, disk[hh]);
	    }
#ifdef  PROHIBIT_EXTRAPOLATION
	  (*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, Nz - 1)] = 0.0;
#endif//PROHIBIT_EXTRAPOLATION
	  //---------------------------------------------------------------
	  /* calculate column density using Simpson's rule */
	  double Sigma = (*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, 0)];
	  for(int jj = 1; jj < Nz - 1; jj++)
	    Sigma += (double)(1 << (1 + (jj & 1))) * (*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, jj)];
	  Sigma += (*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, Nz - 1)];
	  Sigma *= 2.0 * zbin / 3.0;/* 2.0 reflects plane symmetry about the equatorial plane (z = 0) */
	  //---------------------------------------------------------------
#if 0
	  if( ii == 55 )
	    for(int jj = 0; jj < Nz; jj++)
	      fprintf(stderr, "%e\t%e\n", disk[hh].ver[jj], (*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, jj)]);
#endif
	  //---------------------------------------------------------------
	  /* calibrate column density */
	  const double Mscale = disk[hh].Sigma[ii] / (DBL_MIN + Sigma);
	  for(int jj = 0; jj < Nz; jj++)
	    (*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, jj)] *= Mscale;
	  //---------------------------------------------------------------
#if 0
	  fprintf(stderr, "# ii = %d: Mscale = %e, Sigma[ii] = %e, Sigma = %e\n", ii, Mscale, disk[hh].Sigma[ii], Sigma);
#endif
	  //---------------------------------------------------------------
	  /* const double errVal = fabs(Mscale - 1.0); */
	  /* const double errVal = (disk[hh].Sigma[ii] != 0.0) ? fabs(Mscale - 1.0) : (0.0); */
	  const double errVal = (disk[hh].Sigma[ii] > NEGLECT_DENSITY_MINIMUM * disk[hh].Sigma[0]) ? fabs(Mscale - 1.0) : (0.0);
#pragma omp flush(errMax)
	  if( errVal > errMax )
#pragma omp critical
	    if( errVal > errMax )
	      errMax = errVal;
	  //---------------------------------------------------------------
	}
#ifdef  PROHIBIT_EXTRAPOLATION
      for(int jj = 0; jj < Nz; jj++)
	(*disk[hh].rhoSum)[INDEX2D(NR, Nz, NR - 1, jj)] = 0.0;
      disk[hh].Sigma[NR - 1] = 0.0;
#endif//PROHIBIT_EXTRAPOLATION
      //-------------------------------------------------------------------
    }/* for(int hh = 0; hh < ndisk; hh++){ */
    //---------------------------------------------------------------------
#if 0
    exit(0);
#endif
    //---------------------------------------------------------------------
#if 0
    const int status = fpclassify(errMax);
    if( (status != FP_NORMAL) && (status != FP_ZERO) ){
      static char msg[64];
      switch( status ){
      case FP_NAN      :	sprintf(msg, "Not a Number"                                    );	break;
      case FP_INFINITE :	sprintf(msg, "either positive infinity or negative inifinity"  );	break;
      case FP_SUBNORMAL:	sprintf(msg, "too small to be represented in normalized format");	break;
      }/* switch( fpclassify(errMax) ){ */
      __KILL__(stderr, "ERROR: errMax is \"%s\".\n", msg);
    }/* if( fpclassify(errMax) != FP_NORMAL ){ */
#endif
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* convergence tests */
    //---------------------------------------------------------------------
    /* convergence test for the column density */
    bool converge = true;
    if( errMax > CONVERGENCE_POTDENSPAIR )      converge = false;
#ifdef  PROGRESS_REPORT_ON
    fprintf(stdout, "# %d-th iteration: column density error is %e, procedure is %s\n", steps, errMax, (converge ? ("    converged") : ("not converged")));
    fflush(stdout);
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------
    /* convergence test for the density field */
    if( converge ){
      //-------------------------------------------------------------------
      static double worst;
      worst = DBL_EPSILON;
      //-------------------------------------------------------------------
      for(int hh = 0; hh < ndisk; hh++){
	//-----------------------------------------------------------------
#pragma omp parallel for
	for(int ii = 0; ii < NR; ii++){
	  //---------------------------------------------------------------
	  double errVal = 0.0;
	  //---------------------------------------------------------------
	  if( (*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, 0)] > NEGLECT_DENSITY_MINIMUM * (*disk[hh].rhoSum)[INDEX2D(NR, Nz, 0, 0)] )
	    for(int jj = 0; jj < Nz; jj++){
	      //-----------------------------------------------------------
	      const double err =
		((*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, jj)] > NEGLECT_DENSITY_MINIMUM * (*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, 0)]) ?
		(fabs(1.0 - (*disk[hh].rho)[INDEX2D(NR, Nz, ii, jj)] / ((*disk[hh].rhoSum)[INDEX2D(NR, Nz, ii, jj)]))) : (0.0);
	      errVal = (errVal > err) ? errVal : err;
	      //-----------------------------------------------------------
	    }/* for(int jj = 0; jj < Nz; jj++){ */
	  //---------------------------------------------------------------
#pragma omp flush(worst)
	  if( errVal > worst )
#pragma omp critical
	    if( errVal > worst )
	      worst = errVal;
	}/* for(int ii = 0; ii < NR; ii++){ */
	//-----------------------------------------------------------------
      }/* for(int hh = 0; hh < ndisk; hh++){ */
      //-------------------------------------------------------------------
      if( worst > CONVERGENCE_POTDENSPAIR )
	converge = false;
      //-------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
      fprintf(stdout, "# %d-th iteration:        density error is %e, procedure is %s\n", steps, worst, (converge ? ("    converged") : ("not converged")));
      fflush(stdout);
#endif//PROGRESS_REPORT_ON
      //-------------------------------------------------------------------
    }/* if( converge ){ */
    //---------------------------------------------------------------------
    for(int hh = 0; hh < ndisk; hh++)
      swapDblArrays(disk[hh].rho, disk[hh].rhoSum);
    //---------------------------------------------------------------------
    /* final confirmation */
    if( converge )
      break;
    //---------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
    steps++;
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# procedure finish after %d iteration(s): NR = %d, Nz = %d\n#\n#\n", 1 + steps, NR, Nz);
  fflush(stdout);
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void allocDiskProfile
(const int ndisk, disk_data **disk,
 double **hor, double **ver, double **pot, double **dPhidR, double **d2PhidR2, double **radSph, double **rhoSph, double **encSph,
 double **rho, double **rhoSum, double **rhoTot, double **Sigma, double **sigmaz, double **enc
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 , double **zd
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* allocate arrays to store physical quantities */
  //-----------------------------------------------------------------------
  *hor = (double *)malloc(NDISKBIN_HOR                * sizeof(double));
  *ver = (double *)malloc(               NDISKBIN_VER * sizeof(double));
  *pot = (double *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(double));
  if( *hor == NULL ){    __KILL__(stderr, "ERROR: failure to allocate hor");  }
  if( *ver == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ver");  }
  if( *pot == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pot");  }
  //-----------------------------------------------------------------------
#ifdef  USE_POTENTIAL_SCALING_SCHEME
  * dPhidR  = (double *)malloc(NDISKBIN_HOR                * sizeof(double));
  *d2PhidR2 = (double *)malloc(NDISKBIN_HOR                * sizeof(double));
#else///USE_POTENTIAL_SCALING_SCHEME
  * dPhidR  = (double *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(double));
  *d2PhidR2 = (double *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(double));
#endif//USE_POTENTIAL_SCALING_SCHEME
  if( * dPhidR  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dPhidR");  }
  if( *d2PhidR2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate d2PhidR2");  }
  //-----------------------------------------------------------------------
  *radSph = (double *)malloc(NDISKBIN_RAD * sizeof(double));
  *rhoSph = (double *)malloc(NDISKBIN_RAD * sizeof(double));
  *encSph = (double *)malloc(NDISKBIN_RAD * sizeof(double));
  if( *radSph == NULL ){    __KILL__(stderr, "ERROR: failure to allocate radSph");  }
  if( *rhoSph == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rhoSph");  }
  if( *encSph == NULL ){    __KILL__(stderr, "ERROR: failure to allocate encSph");  }
  //-----------------------------------------------------------------------
  *rho    = (double *)malloc(NDISKBIN_HOR * NDISKBIN_VER * ndisk * sizeof(double));
  *rhoSum = (double *)malloc(NDISKBIN_HOR * NDISKBIN_VER * ndisk * sizeof(double));
  if( ndisk > 1 )
    *rhoTot = (double *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(double));
  if( *rho    == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rho");  }
  if( *rhoSum == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rhoSum");  }
  if( ndisk > 1 )
    if( *rhoTot == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rhoTot");  }
  //-----------------------------------------------------------------------
  *Sigma  = (double *)malloc(NDISKBIN_HOR * ndisk * sizeof(double));
  *sigmaz = (double *)malloc(NDISKBIN_HOR * ndisk * sizeof(double));
  *enc    = (double *)malloc(NDISKBIN_HOR * ndisk * sizeof(double));
  if( *Sigma  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate Sigma");  }
  if( *sigmaz == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sigmaz");  }
  if( *enc    == NULL ){    __KILL__(stderr, "ERROR: failure to allocate enc");  }
  //-----------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  *zd = (double *)malloc(NDISKBIN_HOR * ndisk * sizeof(double));
  if( *zd == NULL ){    __KILL__(stderr, "ERROR: failure to allocate zd");  }
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* allocate utility structure and commit arrays */
  //-----------------------------------------------------------------------
  *disk = (disk_data *)malloc(ndisk * sizeof(disk_data));  if( *disk == NULL ){    __KILL__(stderr, "ERROR: failure to allocate disk" );  }
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    /* common arrays for all components */
    //---------------------------------------------------------------------
    (*disk)[ii].hor = *hor;
    (*disk)[ii].ver = *ver;
    (*disk)[ii].pot = *pot;
    //---------------------------------------------------------------------
    (*disk)[ii]. dPhidR  = * dPhidR ;
    (*disk)[ii].d2PhidR2 = *d2PhidR2;
    //---------------------------------------------------------------------
    (*disk)[ii].radSph = *radSph;
    (*disk)[ii].rhoSph = *rhoSph;
    (*disk)[ii].encSph = *encSph;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* individual arrays for each component */
    //---------------------------------------------------------------------
    (*disk)[ii].arr0 = &((*rho)   [ii * NDISKBIN_HOR * NDISKBIN_VER]);    (*disk)[ii].rho    = &(*disk)[ii].arr0;
    (*disk)[ii].arr1 = &((*rhoSum)[ii * NDISKBIN_HOR * NDISKBIN_VER]);    (*disk)[ii].rhoSum = &(*disk)[ii].arr1;
    //---------------------------------------------------------------------
    (*disk)[ii].Sigma  = &((*Sigma )[ii * NDISKBIN_HOR]);
    (*disk)[ii].sigmaz = &((*sigmaz)[ii * NDISKBIN_HOR]);
    (*disk)[ii].enc    = &((*enc   )[ii * NDISKBIN_HOR]);
    //---------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
    (*disk)[ii].zd = &((*zd)[ii * NDISKBIN_HOR]);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* additional setting in case of multiple component disk */
    (*disk)[ii].rhoTot = *rhoTot;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  freeDiskProfile
(const int ndisk, disk_data  *disk,
 double  *hor, double  *ver, double  *pot, double  *dPhidR, double  *d2PhidR2, double  *radSph, double  *rhoSph, double  *encSph,
 double  *rho, double  *rhoSum, double  *rhoTot, double  *Sigma, double  *sigmaz, double  *enc
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 , double  *zd
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  free(hor);  free(ver);  free(pot);
  //-----------------------------------------------------------------------
  free( dPhidR );
  free(d2PhidR2);
  //-----------------------------------------------------------------------
  free(radSph);  free(rhoSph);  free(encSph);
  //-----------------------------------------------------------------------
  free(rho);  free(rhoSum);
  if( ndisk > 1 )
    free(rhoTot);
  //-----------------------------------------------------------------------
  free(Sigma);  free(sigmaz);  free(enc);
  //-----------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  free(zd);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  //-----------------------------------------------------------------------
  free(disk);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void makeDiskPotentialTable(const int ndisk, disk_data * restrict disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  /* initialize the table for Gaussian Quadrature provided by GSL */
  for(int ii = 0; ii < NTBL_GAUSS_QD_LOW; ii++){
    gsl_gaussQD_low_pos   [ii] = 0.0;
    gsl_gaussQD_low_weight[ii] = 0.0;
  }
  gsl_integration_glfixed_table *tab;
  tab = gsl_integration_glfixed_table_alloc(NINTBIN_LOW);
  int max = (NINTBIN_LOW >> 1) + (NINTBIN_LOW & 1);
  for(int ii = 0; ii < max; ii++){
    gsl_gaussQD_low_pos   [ii] = (*tab).x[(max - 1) - ii];
    gsl_gaussQD_low_weight[ii] = (*tab).w[(max - 1) - ii];
  }
  gsl_integration_glfixed_table_free(tab);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* memory for the iterative solver of the Poisson equation */
  //-----------------------------------------------------------------------
  /* sparse matrix in CRS format */
  crs mat, ilu;
  static double mat_val[NNZ_CG], ilu_val[NNZ_CG];
  static    int mat_col[NNZ_CG], ilu_col[NNZ_CG], mat_row[NROW_CG + 1], ilu_row[NROW_CG + 1];
  mat.val = mat_val;  mat.col = mat_col;  mat.row = mat_row;
  ilu.val = ilu_val;  ilu.col = ilu_col;  ilu.row = ilu_row;
  //-----------------------------------------------------------------------
  /* vector used in BiCGSTAB method */
  static double vec[NCOL_CG], res[NCOL_CG], sdw[NCOL_CG], mid[NCOL_CG], tmp[NCOL_CG];
  static double Api[NCOL_CG], Ati[NCOL_CG], Kri[NCOL_CG], Kpi[NCOL_CG], Kti[NCOL_CG];
  //-----------------------------------------------------------------------
  /* outer boundary condition of the potential field */
  static double Phi_NR[NDISKBIN_VER], Phi_Nz[NDISKBIN_HOR];
  //-----------------------------------------------------------------------
  /* execute first touch for memories related to BiCGSTAB */
#pragma omp parallel
  {
#pragma omp for nowait
    for(int ii = 0; ii < NCOL_CG; ii++){
      vec[ii] = res[ii] = sdw[ii] = mid[ii] = tmp[ii] = 0.0;
      Api[ii] = Ati[ii] = Kri[ii] = Kpi[ii] = Kti[ii] = 0.0;
    }/* for(int ii = 0; ii < NCOL_CG; ii++){ */
#pragma omp for nowait
    for(int ii = 0; ii < NNZ_CG; ii++){
      mat_val[ii] = 0.0;
      mat_col[ii] = 0;
    }/* for(int ii = 0; ii < NNZ_CG; ii++){ */
#pragma omp for nowait
    for(int ii = 0; ii < NROW_CG; ii++)
      mat_row[ii] = 0;
  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* load disk data (individual settings for each component) */
  //-----------------------------------------------------------------------
  for(int hh = 0; hh < ndisk; hh++){
    //---------------------------------------------------------------------
    disk[hh].invRd = 1.0 / disk[hh].cfg->rs;
    /* disk[hh].invzd = 1.0 / disk[hh].cfg->zd; */
    //---------------------------------------------------------------------
    if( disk[hh].cfg->kind == SERSIC ){
      /* disk[hh].sersic_ninv = 1.0 / disk[hh].cfg->n_sersic; */
      /* disk[hh].sersic_b    =       disk[hh].cfg->b_sersic; */
      disk[hh].util.sersic_ninv = 1.0 / disk[hh].cfg->n_sersic;
      disk[hh].util.sersic_b    =       disk[hh].cfg->b_sersic;
    }/* if( disk[hh].cfg->kind == SERSIC ){ */
    //---------------------------------------------------------------------
    switch( disk[hh].cfg->kind ){
    case EXP_DISK:      disk[hh].getColumnDensity = getColumnDensityExp   ;      break;
    case   SERSIC:      disk[hh].getColumnDensity = getColumnDensitySersic;      break;
    default:      __KILL__(stderr, "ERROR: undefined model(%d) is specified as disk profile\n", disk[hh].cfg->kind);      break;
    }
    //---------------------------------------------------------------------
/* #ifdef  ENABLE_VARIABLE_SCALE_HEIGHT */
/*     initDiskThicknessControler(&disk[hh]); */
/* #endif//ENABLE_VARIABLE_SCALE_HEIGHT */
    //---------------------------------------------------------------------
  }/* for(int hh = 0; hh < ndisk; hh++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set box size */
  //-----------------------------------------------------------------------
  double zmax = disk[0].cfg->zd;
  double Rmax = disk[0].cfg->rs;
  //-----------------------------------------------------------------------
  for(int ii = 1; ii < ndisk; ii++){
    if( disk[ii].cfg->zd > zmax )      zmax = disk[ii].cfg->zd;
    if( disk[ii].cfg->rs > Rmax )      Rmax = disk[ii].cfg->rs;
  }/* for(int ii = 1; ii < ndisk; ii++){ */
  //-----------------------------------------------------------------------
  __NOTE__("zmax = %e\n", zmax);
  __NOTE__("Rmax = %e\n", Rmax);
  //-----------------------------------------------------------------------
  /* zmax = ldexp(1.0, (int)ceil(log2(zmax * diskMaxHeight))); */
  /* Rmax = ldexp(1.0, (int)ceil(log2(Rmax * diskMaxLength))); */
  /* zmax = 0.5 * Rmax; */
  //-----------------------------------------------------------------------
  zmax = ldexp(1.0, (int)ceil(log2((double)NHOR_OVER_NVER * zmax * diskMaxHeight)));
  Rmax = ldexp(1.0, (int)ceil(log2(                         Rmax * diskMaxLength)));
  if( zmax > Rmax )    Rmax = zmax;
  zmax = (diskMaxHeight / diskMaxLength) * Rmax / (double)NHOR_OVER_NVER;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set cutoff scale in the horizontal direction and set Sigma0 */
  //-----------------------------------------------------------------------
  for(int hh = 0; hh < ndisk; hh++){
    //---------------------------------------------------------------------
    /* additional setting about density cutoff in the horizontal direction */
    /* disk[hh].Rcutoff    = Rmax * 10.0;/\* i.e., the point at infinity *\/ */
    /* disk[hh].invRsmooth = disk[hh].invRd; */
    disk[hh].util.Rcutoff    = Rmax * 10.0;/* i.e., the point at infinity */
    disk[hh].util.invRsmooth = disk[hh].invRd;
    if( disk[hh].cfg->cutoff ){
      //-------------------------------------------------------------------
      /* disk[hh].Rcutoff    =       disk[hh].cfg->rc; */
      /* disk[hh].invRsmooth = 1.0 / disk[hh].cfg->rc_width; */
      disk[hh].util.Rcutoff    =       disk[hh].cfg->rc;
      disk[hh].util.invRsmooth = 1.0 / disk[hh].cfg->rc_width;
      //-------------------------------------------------------------------
    }/* if( disk[hh].cfg->cutoff ){ */
    //---------------------------------------------------------------------
#ifdef  NDIVIDE_GAUSSQD4DISK
    const double Rbin = Rmax / (double)NDIVIDE_GAUSSQD4DISK;
    double sum = 0.0;
    double Rin = 0.0;
    for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){
      sum += gaussQD4encSigma(Rin, Rin + Rbin, disk[hh]);
      Rin += Rbin;
    }/* for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){ */
    disk[hh].cfg->Sigma0 = disk[hh].cfg->Mtot / (2.0 * M_PI * sum);
#else///NDIVIDE_GAUSSQD4DISK
    disk[hh].cfg->Sigma0 = disk[hh].cfg->Mtot / (2.0 * M_PI * gaussQD4encSigma(0.0, Rmax, disk[hh]));
#endif//NDIVIDE_GAUSSQD4DISK
    //---------------------------------------------------------------------
  }/* for(int hh = 0; hh < ndisk; hh++){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate GD formalism */
  //-----------------------------------------------------------------------
#ifdef  USE_GD_FORM_POTENTIAL
  //-----------------------------------------------------------------------
  double *stock_inv;  stock_inv = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  double *stock_sub;  stock_sub = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  double *stock_tmp;  stock_tmp = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  double *stock_sum;  stock_sum = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  if( stock_inv == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_inv");  }
  if( stock_sub == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_sub");  }
  if( stock_tmp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_tmp");  }
  if( stock_sum == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_sum");  }
  //-----------------------------------------------------------------------
  const double    abin = Rmax / (double)NRADBIN;
  const double invabin = 1.0 / abin;
  for(int hh = 0; hh < ndisk; hh++){
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      //-------------------------------------------------------------------
      disk[hh].prf[ii].rho = (double)ii * abin;/* store the position ``a'' */
#ifdef  NDIVIDE_GAUSSQD4DISK
      double sum = 0.0;
      double Rin = disk[hh].prf[ii].rho;
      const double a2 = Rin * Rin;
      const double Rbin = (Rmax - Rin) / (double)NDIVIDE_GAUSSQD4DISK;
      for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){
	sum += gaussQD4GDpot(Rin, Rin + Rbin, a2, disk[hh]);
	Rin += Rbin;
      }/* for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){ */
      disk[hh].prf[ii].enc = disk[hh].cfg->Sigma0 * sum;/* store the integral */
#else///NDIVIDE_GAUSSQD4DISK
      disk[hh].prf[ii].enc = disk[hh].cfg->Sigma0 * gaussQD4GDpot(disk[hh].prf[ii].rho, Rmax, disk[hh].prf[ii].rho * disk[hh].prf[ii].rho, disk[hh]);/* store the integral */
#endif//NDIVIDE_GAUSSQD4DISK
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */
    //---------------------------------------------------------------------
    disk[hh].prf[0].psi = (disk[hh].prf[1].enc - disk[hh].prf[0].enc) * invabin;
#pragma omp parallel for
    for(int ii = 1; ii < NRADBIN - 1; ii++)
      disk[hh].prf[ii].psi = (disk[hh].prf[ii + 1].enc - disk[hh].prf[ii - 1].enc) * 0.5 * invabin;/* store the derivative */
    disk[hh].prf[NRADBIN - 1].psi = (disk[hh].prf[NRADBIN - 1].enc - disk[hh].prf[NRADBIN - 2].enc) * invabin;
    //---------------------------------------------------------------------
  }/* for(int hh = 0; hh < ndisk; hh++){ */
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < NRADBIN; ii++){
    fprintf(stderr, "%e\t%e", disk[0].prf[ii].rho, disk[0].prf[ii].psi);
    for(int hh = 1; hh < ndisk; hh++)
      fprintf(stderr, "\t%e", disk[hh].prf[ii].psi);
    fprintf(stderr, "\n");
  }
  exit(0);
#endif
  //-----------------------------------------------------------------------
#endif//USE_GD_FORM_POTENTIAL
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set reduced table size */
  int NR = NDISKBIN_HOR >> (INCREASE_TABLE_STEPS -1);
  int Nz = NDISKBIN_VER >> (INCREASE_TABLE_STEPS -1);
  for(int ii = 0; ii < INCREASE_TABLE_STEPS; ii++){
    //---------------------------------------------------------------------
    getPotDensPair(ndisk, disk, NR, Rmax, Nz, zmax, Phi_NR, Phi_Nz,
		   mat, ilu, vec, res, sdw, mid, tmp, Api, Ati, Kri, Kpi, Kti
#ifdef  USE_GD_FORM_POTENTIAL
		   , stock_inv, stock_sub, stock_tmp, stock_sum
#endif//USE_GD_FORM_POTENTIAL
		   );
    //---------------------------------------------------------------------
    /* increase table size for the next step */
    if( ii < (INCREASE_TABLE_STEPS - 1) ){
      //-------------------------------------------------------------------
      for(int hh = 0; hh < ndisk; hh++)
      	swapDblArrays(disk[hh].rho, disk[hh].rhoSum);
      //-------------------------------------------------------------------
      const int oldNR = NR;      NR <<= 1;
      const int oldNz = Nz;      Nz <<= 1;
      //-------------------------------------------------------------------
      for(int hh = 0; hh < ndisk; hh++){
	//-----------------------------------------------------------------
#ifdef  SOPHISTICATED_REFINEMENT
	const double rho_min = NEGLECT_DENSITY_MINIMUM * (*disk[hh].rhoSum)[INDEX2D(oldNR, oldNz, 0, 0)];
#endif//SOPHISTICATED_REFINEMENT
	//-----------------------------------------------------------------
	for(int jj = oldNR - 1; jj >= 0; jj--){
	  //---------------------------------------------------------------
#ifdef  SOPHISTICATED_REFINEMENT
	  const double R0 = disk[hh].hor[jj];
	  const double R0inv = 1.0 / R0;
	  const double Rpinv = 1.0 / (4.0 * R0 + disk[hh].Rbin);
	  const double Rminv = 1.0 / (4.0 * R0 - disk[hh].Rbin);
#endif//SOPHISTICATED_REFINEMENT
	  //---------------------------------------------------------------
	  for(int kk = oldNz - 1; kk >= 0; kk--){
	    //-------------------------------------------------------------
	    /* the array rhoSum[oldNR][oldNz] contains old density table */
	    const double rho_org = (*disk[hh].rhoSum)[INDEX2D(oldNR, oldNz, jj, kk)];
	    //-------------------------------------------------------------
#ifdef  SOPHISTICATED_REFINEMENT
	    //-------------------------------------------------------------
	    if( rho_org > rho_min ){
	      //-----------------------------------------------------------
	      const double rho_inv = 1.0 / (DBL_MIN + rho_org);
	      //-----------------------------------------------------------
	      double aa, bb;
	      if( (jj < oldNR - 1) && (jj > 0) )
		aa = ((*disk[hh].rhoSum)[INDEX2D(oldNR, oldNz, jj + 1, kk)] - (*disk[hh].rhoSum)[INDEX2D(oldNR, oldNz, jj - 1, kk)]) * rho_inv * 0.5 * disk[hh].invRbin;
	      else
		if( jj == 0 )		aa = ((*disk[hh].rhoSum)[INDEX2D(oldNR, oldNz, 1, kk)] * rho_inv - 1.0) * 0.5 * disk[hh].invRbin;
		else		aa = 0.0;
	      if( (kk < oldNz - 1) && (kk > 0) )
		bb = ((*disk[hh].rhoSum)[INDEX2D(oldNR, oldNz, jj, kk + 1)] - (*disk[hh].rhoSum)[INDEX2D(oldNR, oldNz, jj, kk - 1)]) * rho_inv * 0.5 * disk[hh].invzbin;
	      else
		if( kk == 0 )		bb = ((*disk[hh].rhoSum)[INDEX2D(oldNR, oldNz, jj, 1)] * rho_inv - 1.0) * 0.5 * disk[hh].invzbin;
		else		bb = 0.0;
	      //-----------------------------------------------------------
	      (*disk[hh].rho)[INDEX2D(NR, Nz, 1 + (jj << 1), 1 + (kk << 1))] = R0 * rho_org * Rpinv * (4.0 + (R0inv + aa) * disk[hh].Rbin + bb * disk[hh].zbin);
	      (*disk[hh].rho)[INDEX2D(NR, Nz, 1 + (jj << 1),     (kk << 1))] = R0 * rho_org * Rpinv * (4.0 + (R0inv + aa) * disk[hh].Rbin - bb * disk[hh].zbin);
	      (*disk[hh].rho)[INDEX2D(NR, Nz,     (jj << 1), 1 + (kk << 1))] = R0 * rho_org * Rminv * (4.0 - (R0inv + aa) * disk[hh].Rbin + bb * disk[hh].zbin);
	      (*disk[hh].rho)[INDEX2D(NR, Nz,     (jj << 1),     (kk << 1))] = R0 * rho_org * Rminv * (4.0 - (R0inv + aa) * disk[hh].Rbin - bb * disk[hh].zbin);
	      //-----------------------------------------------------------
	    }/* if( rho_org > rho_min ){ */
	    else{
	      //-----------------------------------------------------------
	      (*disk[hh].rho)[INDEX2D(NR, Nz, 1 + (jj << 1), 1 + (kk << 1))] = rho_org;
	      (*disk[hh].rho)[INDEX2D(NR, Nz, 1 + (jj << 1),     (kk << 1))] = rho_org;
	      (*disk[hh].rho)[INDEX2D(NR, Nz,     (jj << 1), 1 + (kk << 1))] = rho_org;
	      (*disk[hh].rho)[INDEX2D(NR, Nz,     (jj << 1),     (kk << 1))] = rho_org;
	      //-----------------------------------------------------------
	    }/* else{ */
	    //-------------------------------------------------------------
#else///SOPHISTICATED_REFINEMENT
	    //-------------------------------------------------------------
	    (*disk[hh].rho)[INDEX2D(NR, Nz, 1 + (jj << 1), 1 + (kk << 1))] = rho_org;
	    (*disk[hh].rho)[INDEX2D(NR, Nz, 1 + (jj << 1),     (kk << 1))] = rho_org;
	    (*disk[hh].rho)[INDEX2D(NR, Nz,     (jj << 1), 1 + (kk << 1))] = rho_org;
	    (*disk[hh].rho)[INDEX2D(NR, Nz,     (jj << 1),     (kk << 1))] = rho_org;
	    //-------------------------------------------------------------
#endif//SOPHISTICATED_REFINEMENT
	    //-------------------------------------------------------------
	  }/* for(int kk = oldNz - 1; kk >= 0; kk--){ */
	  //---------------------------------------------------------------
	}/* for(int jj = oldNR - 1; jj >= 0; jj--){ */
        //-----------------------------------------------------------------
#ifdef  PROHIBIT_EXTRAPOLATION
	for(int kk = oldNz - 1; kk >= 0; kk--)
	  (*disk[hh].rho)[INDEX2D(NR, Nz, NR - 2, kk)] = (*disk[hh].rho)[INDEX2D(NR, Nz, NR - 3, kk)];
#endif//PROHIBIT_EXTRAPOLATION
        //-----------------------------------------------------------------
      }/* for(int hh = 0; hh < ndisk; hh++){ */
      //-------------------------------------------------------------------
      for(int jj = oldNR - 1; jj >= 0; jj--)
	for(int kk = oldNz - 1; kk >= 0; kk--){
	  //---------------------------------------------------------------
	  const double pot_org = disk[0].pot[INDEX2D(oldNR, oldNz, jj, kk)];
	  //---------------------------------------------------------------
#ifdef  SOPHISTICATED_REFINEMENT
	  //---------------------------------------------------------------
	  double aa, bb;
	  if( (jj < oldNR - 1) && (jj > 0) )
	    aa = (disk[0].pot[INDEX2D(oldNR, oldNz, jj + 1, kk)] - disk[0].pot[INDEX2D(oldNR, oldNz, jj - 1, kk)]) * 0.5 * disk[0].invRbin;
	  else
	    if( jj == 0 )	      aa = (disk[0].pot[INDEX2D(oldNR, oldNz, 1, kk)] - pot_org) * 0.5 * disk[0].invRbin;
	    else		aa = 0.0;
	  if( (kk < oldNz - 1) && (kk > 0) )
	    bb = (disk[0].pot[INDEX2D(oldNR, oldNz, jj, kk + 1)] - disk[0].pot[INDEX2D(oldNR, oldNz, jj, kk - 1)]) * 0.5 * disk[0].invzbin;
	  else
	    if( kk == 0 )	      bb = (disk[0].pot[INDEX2D(oldNR, oldNz, jj, 1)] - pot_org) * 0.5 * disk[0].invzbin;
	    else		bb = 0.0;
	  //---------------------------------------------------------------
	  disk[0].pot[INDEX2D(NR, Nz, 1 + (jj << 1), 1 + (kk << 1))] = pot_org + 0.25 * disk[0].Rbin * aa + 0.25 * disk[0].zbin * bb;
	  disk[0].pot[INDEX2D(NR, Nz, 1 + (jj << 1),     (kk << 1))] = pot_org + 0.25 * disk[0].Rbin * aa - 0.25 * disk[0].zbin * bb;
	  disk[0].pot[INDEX2D(NR, Nz,     (jj << 1), 1 + (kk << 1))] = pot_org - 0.25 * disk[0].Rbin * aa + 0.25 * disk[0].zbin * bb;
	  disk[0].pot[INDEX2D(NR, Nz,     (jj << 1),     (kk << 1))] = pot_org - 0.25 * disk[0].Rbin * aa - 0.25 * disk[0].zbin * bb;
	  //---------------------------------------------------------------
#else///SOPHISTICATED_REFINEMENT
	  //---------------------------------------------------------------
	  disk[0].pot[INDEX2D(NR, Nz, 1 + (jj << 1), 1 + (kk << 1))] = pot_org;
	  disk[0].pot[INDEX2D(NR, Nz, 1 + (jj << 1),     (kk << 1))] = pot_org;
	  disk[0].pot[INDEX2D(NR, Nz,     (jj << 1), 1 + (kk << 1))] = pot_org;
	  disk[0].pot[INDEX2D(NR, Nz,     (jj << 1),     (kk << 1))] = pot_org;
	  //---------------------------------------------------------------
#endif//SOPHISTICATED_REFINEMENT
	  //---------------------------------------------------------------
	}/* for(int kk = oldNz - 1; kk >= 0; kk--){ */
      //-------------------------------------------------------------------
    }/* if( ii < (INCREASE_TABLE_STEPS - 1) ){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < INCREASE_TABLE_STEPS; ii++){ */
  //-----------------------------------------------------------------------
#ifdef  USE_GD_FORM_POTENTIAL
  free(stock_inv);
  free(stock_sub);
  free(stock_tmp);
  free(stock_sum);
#endif//USE_GD_FORM_POTENTIAL
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate 2 \pi \int_0^R dR' R' Sigma(R') */
  //-----------------------------------------------------------------------
#ifdef  NDIVIDE_GAUSSQD4DISK
  for(int hh = 0; hh < ndisk; hh++){
    const int nsub = NDISKBIN_HOR / NDIVIDE_GAUSSQD4DISK;
    int head = 0;
    double Min = 0.0;
    double Rin = 0.0;
    for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){
      const int tail = head + nsub;
#pragma omp parallel for
      for(int ii = head; ii < tail; ii++)
	disk[hh].enc[ii] = Min + gaussQD4encSigma(Rin, disk[hh].hor[ii], disk[hh]);
      head = tail;
      Rin = disk[hh].hor[head - 1];
      Min = disk[hh].enc[head - 1];
    }/* for(int iter = 0; iter < 4; iter++){ */
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
      disk[hh].enc[ii] *= 2.0 * M_PI * disk[hh].cfg->Sigma0;
  }/* for(int hh = 0; hh < ndisk; hh++){ */
#else///NDIVIDE_GAUSSQD4DISK
  for(int hh = 0; hh < ndisk; hh++)
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
      disk[hh].enc[ii] = 2.0 * M_PI * disk[hh].cfg->Sigma0 * gaussQD4encSigma(0.0, disk[hh].hor[ii], disk[hh]);
#endif//NDIVIDE_GAUSSQD4DISK
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < NDISKBIN_HOR; ii++){
    fprintf(stderr, "%e", disk[0].hor[ii]);
    for(int hh = 0; hh < ndisk; hh++)
      fprintf(stderr, "\t%.20e\t%.20e", disk[hh].Sigma[ii], disk[hh].enc[ii]);
    fprintf(stderr, "\n");
  }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate int_0^z dz' \rho(R, z') */
  //-----------------------------------------------------------------------
  static double xx[NDISKBIN_VER + 2], ff_ful[(NDISKBIN_VER + 2) * CPUS_PER_PROCESS], f2_ful[(NDISKBIN_VER + 2) * CPUS_PER_PROCESS], bp_ful[(NDISKBIN_VER + 2) * CPUS_PER_PROCESS];
  xx[0] = 0.0;
  for(int ii = 0; ii < NDISKBIN_VER; ii++)
    xx[1 + ii] = disk[0].ver[ii];
  xx[NDISKBIN_VER + 1] = xx[NDISKBIN_VER] + 0.5 * disk[0].zbin;
  const double z02 = disk[0].ver[0] * disk[0].ver[0];
  const double z12 = disk[0].ver[1] * disk[0].ver[1];
  const double zd2inv = 1.0 / (z12 - z02);

  for(int hh = 0; hh < ndisk; hh++)
#pragma omp parallel for num_threads(CPUS_PER_PROCESS)
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      //-------------------------------------------------------------------
      /* apply cubic spline interpolation */
      //-------------------------------------------------------------------
      const int tidx = omp_get_thread_num();
      double *ff = &ff_ful[(NDISKBIN_VER + 2) * tidx];
      double *f2 = &f2_ful[(NDISKBIN_VER + 2) * tidx];
      double *bp = &bp_ful[(NDISKBIN_VER + 2) * tidx];
      //-------------------------------------------------------------------
      ff[NDISKBIN_VER + 1] = (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)];
      for(int jj = NDISKBIN_VER - 1; jj >= 0; jj--)
	ff[1 + jj] = (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)];
      ff[0] = ((*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 0)] * z12 - (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 1)] * z02) * zd2inv;
      //-------------------------------------------------------------------
      genCubicSpline1D(NDISKBIN_VER + 2, xx, ff, bp, 0.0, NATURAL_CUBIC_SPLINE, f2);
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* numerical integral based on cubic spline interpolation */
      //-------------------------------------------------------------------
      double sum = getCubicSplineIntegral1D(0, xx, ff, f2);
      for(int jj = 0; jj < NDISKBIN_VER; jj++){
	//-----------------------------------------------------------------
	sum += getCubicSplineIntegral1D1stHalf(1 + jj, xx, ff, f2);
	(*disk[hh].rhoSum)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] = sum;
	//-----------------------------------------------------------------
	sum += getCubicSplineIntegral1D2ndHalf(1 + jj, xx, ff, f2);
	//-----------------------------------------------------------------
      }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
      //-------------------------------------------------------------------
      const double inv = 1.0 / (DBL_MIN + (*disk[hh].rhoSum)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)]);
      for(int jj = NDISKBIN_VER - 1; jj >= 0; jj--)
	(*disk[hh].rhoSum)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] *= inv;
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  //-----------------------------------------------------------------------
#if 0
  for(int jj = 0; jj < NDISKBIN_VER; jj++)
    fprintf(stderr, "%e\t%e\n", disk[0].ver[jj], (*disk[0].rhoSum)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj)]);
  exit(0);
#endif
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < NDISKBIN_HOR; ii++){
    for(int jj = 0; jj < NDISKBIN_VER; jj++){
      fprintf(stderr, "%e\t%e", disk[0].hor[ii], disk[0].ver[jj]);
      for(int hh = 0; hh < ndisk; hh++)
	fprintf(stderr, "\t%e", (*disk[hh].rhoSum)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)]);
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline double func4enc(const double rr, double * restrict rad, const double invdr, double * restrict rho)
{
  //-----------------------------------------------------------------------
  double fr;
  const int ir = bisection(rr, NDISKBIN_RAD, rad, false, invdr, &fr);
  //-----------------------------------------------------------------------
  return (rr * rr * ((1.0 - fr) * rho[ir] + fr * rho[1 + ir]));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
double gaussQD4enc(const double min, const double max, double * restrict rad, const double invdr, double * restrict rho)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  //-----------------------------------------------------------------------
  if( NINTBIN & 1 )
    sum = gsl_gaussQD_weight[(NINTBIN >> 1)] * func4enc(pls + mns * gsl_gaussQD_pos[(NINTBIN >> 1)], rad, invdr, rho);
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (func4enc(pls + mns * gsl_gaussQD_pos[ii], rad, invdr, rho) + func4enc(pls - mns * gsl_gaussQD_pos[ii], rad, invdr, rho));
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void integrateSphericalDensityProfile(const int ndisk, disk_data *disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#ifdef  BILINEAR_INTERPOLATION
  static double rho_grid[NDISKBIN_HOR + 1][NDISKBIN_VER + 1];
#endif//BILINEAR_INTERPOLATION
#ifdef  BICUBIC_INTERPOLATION
  static double rho_grid[(NDISKBIN_VER + 2) * (NDISKBIN_HOR + 2)], f2spline[(NDISKBIN_VER + 2) * (NDISKBIN_HOR + 2)];
  static double RR1d[NDISKBIN_HOR + 2], zz1d[NDISKBIN_VER + 2];
  static double bpspline[(NDISKBIN_HOR + 2) * CPUS_PER_PROCESS], ff_1dtot[(NDISKBIN_HOR + 2) * CPUS_PER_PROCESS], f2_1dtot[(NDISKBIN_HOR + 2) * CPUS_PER_PROCESS], bp_1dtot[(NDISKBIN_HOR + 2) * CPUS_PER_PROCESS];
  RR1d[0] = 0.0;  for(int ii = 0; ii < NDISKBIN_HOR; ii++)    RR1d[1 + ii] = disk[0].hor[ii];
  RR1d[NDISKBIN_HOR + 1] = RR1d[NDISKBIN_HOR] + disk[0].Rbin;
  zz1d[0] = 0.0;  for(int jj = 0; jj < NDISKBIN_VER; jj++)    zz1d[1 + jj] = disk[0].ver[jj];
  zz1d[NDISKBIN_VER + 1] = zz1d[NDISKBIN_VER] + disk[0].zbin;
#endif//BICUBIC_INTERPOLATION
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int hh = 0; hh < ndisk; hh++){
    //---------------------------------------------------------------------
    /* load disk data */
    //---------------------------------------------------------------------
    const double Mdisk = disk[hh].cfg->Mtot;
    //---------------------------------------------------------------------
    double *rad;    rad = disk[hh].radSph;
    double *rho;    rho = disk[hh].rhoSph;
    double *enc;    enc = disk[hh].encSph;
    //---------------------------------------------------------------------
    profile *prf;    prf = disk[hh].prf;
    const double logrbin_sph = disk[hh].logrbin;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* calculate spherical averaged profiles of (volume-)density and enclosed mass */
    //---------------------------------------------------------------------
    /* initialization */
    const double Rbin = disk[hh].Rbin;    const double Rmax = disk[hh].hor[NDISKBIN_HOR - 1] + 0.5 * Rbin;
    const double zbin = disk[hh].zbin;    const double zmax = disk[hh].ver[NDISKBIN_VER - 1] + 0.5 * zbin;
#ifndef BICUBIC_INTERPOLATION
    const double invRbin = disk[hh].invRbin;
    const double invzbin = disk[hh].invzbin;
#endif//BICUBIC_INTERPOLATION
    const double    rrbin = sqrt(Rmax * Rmax + zmax * zmax) / (double)NDISKBIN_RAD;
    const double invrrbin = 1.0 / rrbin;
    for(int ii = 0; ii < NDISKBIN_RAD; ii++){
      rad[ii] = (0.5 + (double)ii) * rrbin;
      rho[ii] = 0.0;
      enc[ii] = 0.0;
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
    //---------------------------------------------------------------------
    /* density assignment */
#ifdef  BILINEAR_INTERPOLATION
    rho_grid[0][0] = (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, 0)];
    for(int jj = 1; jj < NDISKBIN_VER; jj++)
      rho_grid[0][jj] = 0.5 * ((*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj - 1)] + (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj)]);
    rho_grid[0][NDISKBIN_VER] = (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, NDISKBIN_VER - 1)];

    for(int ii = 1; ii < NDISKBIN_HOR; ii++){
      rho_grid[ii][0] = 0.5 * ((*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, 0)] + (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 0)]);
      for(int jj = 1; jj < NDISKBIN_VER; jj++)
	rho_grid[ii][jj] = 0.25 * ((*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, jj - 1)] + (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, jj)] + (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj - 1)] + (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)]);
      rho_grid[ii][NDISKBIN_VER] = 0.5 * ((*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, NDISKBIN_VER - 1)] + (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)]);
    }

    rho_grid[NDISKBIN_HOR][0] = (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, 0)];
    for(int jj = 1; jj < NDISKBIN_VER; jj++)
      rho_grid[NDISKBIN_HOR][jj] = 0.5 * ((*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj - 1)] + (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)]);
    rho_grid[NDISKBIN_HOR][NDISKBIN_VER] = (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, NDISKBIN_VER - 1)];
#endif//BILINEAR_INTERPOLATION
    const double invRSub = 1.0 / (double)SUBDIVIDE_NUM_HOR;
    const double invzSub = 1.0 / (double)SUBDIVIDE_NUM_VER;
#ifdef  BICUBIC_INTERPOLATION
    /* copy main domain */
    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	rho_grid[INDEX2D(NDISKBIN_VER + 2, NDISKBIN_HOR + 2, jj + 1, ii + 1)] = (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)];
    /* set outer boundary */
    for(int jj = 0; jj < NDISKBIN_VER + 2; jj++)
      rho_grid[INDEX2D(NDISKBIN_VER + 2, NDISKBIN_HOR + 2, jj, NDISKBIN_HOR + 1)] = 0.0;
    for(int ii = 0; ii < NDISKBIN_HOR + 2; ii++)
      rho_grid[INDEX2D(NDISKBIN_VER + 2, NDISKBIN_HOR + 2, NDISKBIN_VER + 1, ii)] = 0.0;
    /* set inner boundary */
    const double z02 = disk[hh].ver[0] * disk[hh].ver[0];
    const double z12 = disk[hh].ver[1] * disk[hh].ver[1];
    const double zd2inv = 1.0 / (z12 - z02);
    for(int ii = 1; ii < NDISKBIN_HOR + 1; ii++)
      rho_grid[INDEX2D(NDISKBIN_VER + 2, NDISKBIN_HOR + 2, 0, ii)] = ((*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, 0)] * z12 - (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, 1)] * z02) * zd2inv;
    const double R02 = disk[hh].hor[0] * disk[hh].hor[0];
    const double R12 = disk[hh].hor[1] * disk[hh].hor[1];
    const double Rd2inv = 1.0 / (R12 - R02);
    for(int jj = 0; jj < NDISKBIN_VER + 1; jj++)
      rho_grid[INDEX2D(NDISKBIN_VER + 2, NDISKBIN_HOR + 2, jj, 0)] = (rho_grid[INDEX2D(NDISKBIN_VER + 2, NDISKBIN_HOR + 2, jj, 1)] * R12 - rho_grid[INDEX2D(NDISKBIN_VER + 2, NDISKBIN_HOR + 2, jj, 2)] * R02) * Rd2inv;
    genCubicSpline2D1st(NDISKBIN_VER + 2, NDISKBIN_HOR + 2, RR1d, rho_grid, f2spline, bpspline);
#pragma omp parallel for num_threads(CPUS_PER_PROCESS)
    for(int ii = 0; ii < NDISKBIN_HOR * SUBDIVIDE_NUM; ii++){
      //-------------------------------------------------------------------
      const int tidx = omp_get_thread_num();
      double *ff_1d = &ff_1dtot[(NDISKBIN_VER + 2) * tidx];
      double *f2_1d = &f2_1dtot[(NDISKBIN_VER + 2) * tidx];
      double *bp_1d = &bp_1dtot[(NDISKBIN_VER + 2) * tidx];
      //-------------------------------------------------------------------
      const double RR = (0.5 + (double)ii) * invRSub * Rbin;
      genCubicSpline2D2nd(RR, NDISKBIN_VER + 2, zz1d, NDISKBIN_HOR + 2, RR1d, rho_grid, f2spline, ff_1d, f2_1d, bp_1d);
      //-------------------------------------------------------------------
      for(int jj = 0; jj < NDISKBIN_VER * SUBDIVIDE_NUM; jj++){
	//-----------------------------------------------------------------
	const double zz = (0.5 + (double)jj) * invzSub * zbin;
	const int idx = (int)floor(sqrt(RR * RR + zz * zz) * invrrbin);
	//-----------------------------------------------------------------
#pragma omp atomic
	rho[idx] += getCubicSpline1D(zz, NDISKBIN_VER + 2, zz1d, ff_1d, f2_1d);
#pragma omp atomic
	enc[idx] += 1.0;
	//-----------------------------------------------------------------
      }/* for(int jj = 0; jj < NDISKBIN_VER * SUBDIVIDE_NUM; jj++){ */
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NDISKBIN_HOR * SUBDIVIDE_NUM; ii++){ */
#else///BICUBIC_INTERPOLATION
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      //-------------------------------------------------------------------
      const double Rmin = disk[hh].hor[ii] - 0.5 * Rbin;
#ifdef  BILINEAR_INTERPOLATION
      const int Ridx = ii;
#endif//BILINEAR_INTERPOLATION
      //-------------------------------------------------------------------
      for(int jj = 0; jj < NDISKBIN_VER; jj++){
	//-----------------------------------------------------------------
	const double zmin = disk[hh].ver[jj] - 0.5 * zbin;
#ifdef  BILINEAR_INTERPOLATION
	const int zidx = jj;
#endif//BILINEAR_INTERPOLATION
#ifndef BILINEAR_INTERPOLATION
	const double rhoSub = (*disk[hh].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)];
#endif//BILINEAR_INTERPOLATION
	//-----------------------------------------------------------------
	for(int kk = 0; kk < SUBDIVIDE_NUM_HOR; kk++){
	  //---------------------------------------------------------------
	  const double RR = Rmin + (0.5 + (double)kk) * invRSub * Rbin;
#ifdef  BILINEAR_INTERPOLATION
	  const double tt = (RR - ((double)Ridx * Rbin)) * invRbin;
#endif//BILINEAR_INTERPOLATION
	  //---------------------------------------------------------------
	  for(int ll = 0; ll < SUBDIVIDE_NUM_VER; ll++){
	    //-------------------------------------------------------------
	    const double zz = zmin + (0.5 + (double)ll) * invzSub * zbin;
	    const int idx = (int)floor(sqrt(RR * RR + zz * zz) * invrrbin);
#ifndef BILINEAR_INTERPOLATION
	    rho[idx] += rhoSub;
#else///BILINEAR_INTERPOLATION
	    const double uu = (zz - ((double)zidx * zbin)) * invzbin;
	    rho[idx] +=
	      ((1.0 - uu) * rho_grid[    Ridx][zidx] + uu * rho_grid[    Ridx][1 + zidx]) * (1.0 - tt) +
	      ((1.0 - uu) * rho_grid[1 + Ridx][zidx] + uu * rho_grid[1 + Ridx][1 + zidx]) *        tt;
#endif//BILINEAR_INTERPOLATION
	    enc[idx] += 1.0;
	    //-------------------------------------------------------------
	  }/* for(int ll = 0; ll < SUBDIVIDE_NUM; ll++){ */
	  //---------------------------------------------------------------
	}/* for(int kk = 0; kk < SUBDIVIDE_NUM; kk++){ */
	//-----------------------------------------------------------------
      }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
#endif//BICUBIC_INTERPOLATION
    //---------------------------------------------------------------------
    for(int ii = 0; ii < NDISKBIN_RAD; ii++){
      //-------------------------------------------------------------------
      rho[ii] /= enc[ii];
      enc[ii] = 0.0;
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NDISKBIN_RAD; ii++){ */
    //---------------------------------------------------------------------
#if 1
  static double rr[NDISKBIN_RAD + 2], ff[NDISKBIN_RAD + 2], f2[NDISKBIN_RAD + 2], bp[NDISKBIN_RAD + 2];
  rr[0] = 0.0;
  for(int ii = 0; ii < NDISKBIN_RAD; ii++)
    rr[1 + ii] = rad[ii];
  rr[NDISKBIN_RAD + 1] = rr[NDISKBIN_RAD] + 0.5 * rrbin;
  ff[NDISKBIN_RAD + 1] = rr[NDISKBIN_RAD + 1] * rr[NDISKBIN_RAD + 1] * rho[NDISKBIN_RAD - 1];
  for(int ii = NDISKBIN_RAD - 1; ii >= 0; ii--)
    ff[1 + ii] = rad[ii] * rad[ii] * rho[ii];
  ff[0] = 0.0;
  genCubicSpline1D(NDISKBIN_RAD + 2, rr, ff, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, f2);
  double sum = getCubicSplineIntegral1D(0, rr, ff, f2);
  for(int ii = 0; ii < NDISKBIN_RAD; ii++){
    //-----------------------------------------------------------------
    sum += getCubicSplineIntegral1D1stHalf(1 + ii, rr, ff, f2);
    enc[ii] = 4.0 * M_PI * sum;
    //-----------------------------------------------------------------
    sum += getCubicSplineIntegral1D2ndHalf(1 + ii, rr, ff, f2);
    //-----------------------------------------------------------------
  }/* for(int ii = 0; ii < NDISKBIN_RAD; jj++){ */
#if 0
  for(int ii = 0; ii < NDISKBIN_RAD; ii++)
    fprintf(stderr, "%e\t%e\n", rad[ii], enc[ii]);
  exit(0);
#endif
#else
    /* numerical integration to get enclosed mass profile using Simpson's rule */
    double Menc[2] = {rad[0] * rad[0] * rho[0], rad[1] * rad[1] * rho[1]};
    double Mini[2] = {0.0, gaussQD4enc(0.0, rad[0], rad, invrrbin, rho)};
    /* enc[0] = rad[0] * Menc[0] / 3.0; */
    /* enc[1] = enc[0] + 0.5 * rad[1] * rrbin * rrbin * 0.5 * (rho[0] + rho[1]); */
    enc[0] = Mini[1];
    enc[1] = gaussQD4enc(0.0, rad[1], rad, invrrbin, rho);
    Menc[0] += Menc[1] * 4.0;
    for(int ii = 2; ii < NDISKBIN_RAD; ii++){
      //-------------------------------------------------------------------
      const double mass = rad[ii] * rad[ii] * rho[ii];
      const int idx = ii & 1;
      //-------------------------------------------------------------------
      /* enc[ii] = enc[idx] + (Menc[idx] + mass) * rrbin / 3.0; */
      enc[ii] = Mini[idx] + (Menc[idx] + mass) * rrbin / 3.0;
      //-------------------------------------------------------------------
      Menc[0] += mass * (double)(1 << (1 + (idx    )));
      Menc[1] += mass * (double)(1 << (1 + (idx ^ 1)));
      //-------------------------------------------------------------------
    }/* for(int ii = 2; ii < NDISKBIN_RAD; ii++){ */
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_RAD; ii++)
      enc[ii] *= 4.0 * M_PI;
#endif
    //---------------------------------------------------------------------
    /* adjust numerical errors */
    const double Mscale = disk[hh].cfg->Mtot / enc[NDISKBIN_RAD - 1];
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_RAD; ii++){
      rho[ii] *= Mscale;
      enc[ii] *= Mscale;
    }/* for(int ii = 0; ii < NDISKBIN_RAD; ii++){ */
#ifdef  PROGRESS_REPORT_ON
    fprintf(stdout, "# spherical averaged density profile correction: multiplying %e for %d-th disk\n", Mscale, hh);
    fflush(stdout);
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* assign spherical averaged profile of density, enclosed mass and potential (internal part) */
    //---------------------------------------------------------------------
    double dummy;
    const int head = findIdxSpherical(rad[               0], prf, 1.0, &dummy);
    const int tail = findIdxSpherical(rad[NDISKBIN_RAD - 1], prf, 1.0, &dummy);
    //---------------------------------------------------------------------
    /* fit logr vs logM using 4 meshes */
    double SS, Sx, Sy[2], Sxx, Sxy[2];
    SS = Sx = Sxx = 0.0;
    for(int ii = 0; ii < 2; ii++)
      Sy[ii] = Sxy[ii] = 0.0;
    for(int ii = 0; ii < 4; ii++){
      const double lograd = log10(rad[ii]);
      const double logenc = log10(enc[ii]);
      const double logrho = log10(rho[ii]);
      SS     += 1.0;
      Sx     += lograd;
      Sxx    += lograd * lograd;
      Sy [0] +=          logenc;      Sy [1] +=          logrho;
      Sxy[0] += lograd * logenc;      Sxy[1] += lograd * logrho;
    }
    double power[2], basis[2];
    for(int ii = 0; ii < 2; ii++){
      power[ii] = (SS * Sxy[ii] - Sx * Sy[ii]) / (SS * Sxx - Sx * Sx);
      basis[ii] = pow(10.0, (Sy[ii] - power[ii] * Sx) / SS);
    }
    /* extrapolate: enc, rho --> prf */
#pragma omp parallel for
    for(int ii = 0; ii < head; ii++){
      prf[ii].enc = basis[0] * pow(prf[ii].rad, power[0]);
      prf[ii].rho = basis[1] * pow(prf[ii].rad, power[1]);
    }
    //---------------------------------------------------------------------
    /* interpolate: enc, rho --> prf */
#pragma omp parallel for
    for(int ii = head; ii < tail; ii++){
      double alpha;
      const int idx = bisection(prf[ii].rad, NDISKBIN_RAD, rad, false, invrrbin, &alpha);
      prf[ii].rho = (1.0 - alpha) * rho[idx] + alpha * rho[1 + idx];
      prf[ii].enc = (1.0 - alpha) * enc[idx] + alpha * enc[1 + idx];
    }/* for(int ii = head; ii < tail; ii++){ */
    //---------------------------------------------------------------------
    /* fill the total enclosed mass in the outer region */
#pragma omp parallel for
    for(int ii = tail; ii < 4 + NRADBIN; ii++){
      prf[ii].rho = 0.0;
      prf[ii].enc = Mdisk;
    }
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < 4 + NRADBIN; ii++)
      prf[ii].psi = prf[ii].enc / prf[ii].rad;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* calculate spherical averaged profile of potential (external part) */
    //---------------------------------------------------------------------
    double Pext[2];
    Pext[0] = prf[NRADBIN + 3].rad * prf[NRADBIN + 3].rad * prf[NRADBIN + 3].rho;
    Pext[1] = prf[NRADBIN + 2].rad * prf[NRADBIN + 2].rad * prf[NRADBIN + 2].rho;
    double Pini[2];
    Pini[0] =           0.0;
    Pini[1] = Pini[0] + (prf[NRADBIN + 3].rad - prf[NRADBIN + 2].rad) * sqrt(prf[NRADBIN + 2].rad * prf[NRADBIN + 3].rad) * 0.5 * (prf[NRADBIN + 2].rho + prf[NRADBIN + 3].rho);
    Pext[0] += Pext[1] * 4.0;
    //---------------------------------------------------------------------
    for(int ii = NRADBIN + 1; ii >= 0; ii--){
      //-------------------------------------------------------------------
      const double psi = prf[ii].rad * prf[ii].rad * prf[ii].rho;
      const int idx = (int)((ii + 1) & 1);
      //-------------------------------------------------------------------
      prf[ii].psi += 4.0 * M_PI * (Pini[idx] + (Pext[idx] + psi) * logrbin_sph * M_LN10 / 3.0);
      //-------------------------------------------------------------------
      Pext[0] += psi * (double)(1 << (1 + (idx    )));
      Pext[1] += psi * (double)(1 << (1 + (idx ^ 1)));
      //-------------------------------------------------------------------
    }
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* multiply overall factor */
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < 4 + NRADBIN; ii++)
      prf[ii].psi *= (double)newton;
    //---------------------------------------------------------------------
  }/* for(int hh = 0; hh < ndisk; hh++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define SMOOTH(a, b, c, d) ((((a) * (b) > 0.0) && ((b) * (c) > 0.0) && ((c) * (d) > 0.0)) ? (0.5 * ((b) + (c))) : (0.0)) */
//-------------------------------------------------------------------------
void diffAxisymmetricPotential(const disk_data disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  double *  Phi    ;    Phi     = disk.pot;
  double *       RR;        RR  = disk.hor;
  double *       zz;        zz  = disk.ver;
  double * dPhi_dR ;   dPhi_dR  = disk. dPhidR;
  double *d2Phi_dR2;  d2Phi_dR2 = disk.d2PhidR2;
  const double invRbin = disk.invRbin;
  profile *prf;  prf = disk.prf;
  const double invlogrbin = disk.invlogrbin;
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < NDISKBIN_HOR - 1; ii++){
    //---------------------------------------------------------------------
#ifdef  USE_POTENTIAL_SCALING_SCHEME
    const int jj = 0;
#else///USE_POTENTIAL_SCALING_SCHEME
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
#endif//USE_POTENTIAL_SCALING_SCHEME
      {
	//-----------------------------------------------------------------
	/* R-derivatives of potential given by the disk component */
	double _dPhidR__disk, d2PhidR2_disk;
	if( ii != 0 ){
	  _dPhidR__disk = 0.5     * invRbin * (Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, jj)] - Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, jj)]);
	  d2PhidR2_disk = invRbin * invRbin * (Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, jj)] + Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, jj)] - 2.0 * Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)]);
	}/* if( ii != 0 ){ */
	else{
	  _dPhidR__disk = 0.5     * invRbin * (Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, jj)] - Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj)]);
	  d2PhidR2_disk = invRbin * invRbin * (Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, jj)] - Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj)]);
	}/* else{ */
	//-----------------------------------------------------------------
	/* r-derivatives of potential given by the spherical component(s) */
	const double rad = sqrt(RR[ii] * RR[ii] + zz[jj] * zz[jj]);
	const double rinv = 1.0 / rad;
	const double drdR = RR[ii] * rinv;
	const double d2rdR2 = (1.0 - drdR * drdR) * rinv;
	/* find index corresponds to r = rad in prf */
	double alpha;
	const int idx = findIdxSpherical(rad, prf, invlogrbin, &alpha);
	/* Psi = -Phi + Phi0 = -Phi (i.e. Phi0 is assumed to be zero) */
	const double enc = (1.0 - alpha) * prf[idx].enc_tot + alpha * prf[1 + idx].enc_tot;
	const double rho = (1.0 - alpha) * prf[idx].rho_tot + alpha * prf[1 + idx].rho_tot;
	const double  dPhidr__sphe = (double)newton *                           enc * rinv * rinv;
	const double d2Phidr2_sphe = (double)newton * (4.0 * M_PI * rho - 2.0 * enc * rinv * rinv * rinv);
	//-----------------------------------------------------------------
	/* R-derivatives of potential given by the all components */
#ifndef USE_POTENTIAL_SCALING_SCHEME
	dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] = _dPhidR__disk + dPhidr__sphe * drdR;
	d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] = d2PhidR2_disk + dPhidr__sphe * drdR * d2rdR2 + d2Phidr2_sphe * drdR * drdR;
#else///USE_POTENTIAL_SCALING_SCHEME
	dPhi_dR  [ii] = _dPhidR__disk + dPhidr__sphe * drdR;
	d2Phi_dR2[ii] = d2PhidR2_disk + dPhidr__sphe * drdR * d2rdR2 + d2Phidr2_sphe * drdR * drdR;
#endif//USE_POTENTIAL_SCALING_SCHEME
	//-----------------------------------------------------------------
      }
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NDISKBIN_HOR - 1; ii++){ */
  //-----------------------------------------------------------------------
#ifdef  USE_POTENTIAL_SCALING_SCHEME
  dPhi_dR  [NDISKBIN_HOR - 1] =  dPhi_dR [NDISKBIN_HOR - 2];
  d2Phi_dR2[NDISKBIN_HOR - 1] = d2Phi_dR2[NDISKBIN_HOR - 2];
#else///USE_POTENTIAL_SCALING_SCHEME
  for(int jj = 0; jj < NDISKBIN_VER; jj++){
    //---------------------------------------------------------------------
    dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] =  dPhi_dR [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj)];
    d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] = d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj)];
    //---------------------------------------------------------------------
  }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
#endif//USE_POTENTIAL_SCALING_SCHEME
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void calcVerticalVdisp(const int ndisk, disk_data *disk_info)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int hh = 0; hh < ndisk; hh++){
    //---------------------------------------------------------------------
    /* load disk data */
    //---------------------------------------------------------------------
    const disk_data disk = disk_info[hh];
    //---------------------------------------------------------------------
    double *sig;  sig = disk.sigmaz;
    double *Phi;  Phi = disk.pot;
    double *hor;  hor = disk.hor;
    double *ver;  ver = disk.ver;
    //---------------------------------------------------------------------
    const double invzbin = disk.invzbin;
    //---------------------------------------------------------------------
    profile *sph;  sph = disk.prf;
    const double invlogrbin = disk.invlogrbin;
    //---------------------------------------------------------------------
#ifndef ENABLE_VARIABLE_SCALE_HEIGHT
    const double zd = disk.cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* calculate vertical velocity dispersion */
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      //-------------------------------------------------------------------
      int i0;
      double a0, Psi_s, Psi_d;
      //-------------------------------------------------------------------
      /* get potential @ (R, z = 0) */
      const double RR = hor[ii];
      i0 = findIdxSpherical(RR, sph, invlogrbin, &a0);
      Psi_s = (1.0 - a0) * sph[i0].psi_tot + a0 * sph[i0 + 1].psi_tot;
      Psi_d = -Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 0)];
      const double Psi_R_0 = Psi_d + Psi_s;
      //-------------------------------------------------------------------
      /* get potential @ (R, z = zd) */
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      /* double zd = disk.cfg->zd * modulateDiskThickness(RR, disk); */
      double zd = disk.zd[ii];
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      const double rad = sqrt(RR * RR + zd * zd);
      i0 = findIdxSpherical(rad, sph, invlogrbin, &a0);
      Psi_s = (1.0 - a0) * sph[i0].psi_tot + a0 * sph[i0 + 1].psi_tot;
      i0 = bisection(zd, NDISKBIN_VER, ver, false, invzbin, &a0);
      Psi_d = (1.0 - a0) * (-Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, i0)]) + a0 * (-Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 1 + i0)]);
      const double Psi_R_zd = Psi_d + Psi_s;
      //-------------------------------------------------------------------
      sig[ii] = sqrt(Psi_R_0 - Psi_R_zd);
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
    //---------------------------------------------------------------------
  }/* for(int hh = 0; hh < ndisk; hh++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void leastSquaredMethod(const int num, double * restrict xx, double * restrict yy, double * restrict pp, double * restrict bb)
{
  //-----------------------------------------------------------------------
  double SS, Sx, Sy, Sxx, Sxy;
  SS = Sx = Sy = Sxx = Sxy = 0.0;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    const double logx = log10(xx[ii]);
    const double logy = log10(yy[ii]);
    SS  += 1.0;
    Sx  += logx;
    Sxx += logx * logx;
    Sy  +=        logy;
    Sxy += logx * logy;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
  *pp = (SS * Sxy - Sx * Sy) / (SS * Sxx - Sx * Sx);
  *bb = pow(10.0, (Sy - (*pp) * Sx) / SS);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void distributeDiskParticles(ulong *Nuse, iparticle body, const real mass, const disk_data disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* load disk data */
  //-----------------------------------------------------------------------
  profile *prf;  prf = disk.prf;
  const double invlogrbin = disk.invlogrbin;
  //-----------------------------------------------------------------------
  double *pot;  pot =  disk.pot;
  double *hor;  hor =  disk.hor;
  double *ver;  ver =  disk.ver;
  double *enc;  enc =  disk.enc;
  double *rho;  rho = *disk.rhoSum;
  double *sig;  sig =  disk.sigmaz;
  //-----------------------------------------------------------------------
  double * dPhidR ;  dPhidR   =  disk. dPhidR;
  double *d2PhidR2;  d2PhidR2 =  disk.d2PhidR2;
  //-----------------------------------------------------------------------
  const ulong num = disk.cfg->num;
  const double frac = disk.cfg->vdisp_frac;
  const double Mmin = enc[0];
  const double Mmax = enc[NDISKBIN_HOR - 1];
  //-----------------------------------------------------------------------
  const ulong rand_half = ((gsl_rng_min(GSLRand) + gsl_rng_max(GSLRand)) >> 1);
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# start distributing disk particles (%zu bodies: [%zu:%zu])\n", num, *Nuse, (*Nuse) + num - 1);
  fflush(stdout);
  const ulong nunit = (ulong)ceilf(0.1f * (float)num);
  ulong stage = 1;
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------
#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
  double Krot = 0.0;
  double Krnd = 0.0;
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION
  //-----------------------------------------------------------------------
  for(ulong ii = *Nuse; ii < *Nuse + num; ii++){
    //---------------------------------------------------------------------
    /* determine Rg, phi, and z */
    //---------------------------------------------------------------------
    /* set Rg */
    const double Renc = Mmax * UNIRAND_DBL;
    double aRg, Rg;
    int iRg;
    if( Renc >= Mmin ){
      //-------------------------------------------------------------------
      /* use table search */
      iRg = bisection(Renc, NDISKBIN_HOR, enc, false, 1.0, &aRg);
      aRg /= (enc[1 + iRg] - enc[iRg]);
      Rg = (1.0 - aRg) * hor[iRg] + aRg * hor[1 + iRg];
      //-------------------------------------------------------------------
    }/* if( Renc >= Mmin ){ */
    else{
      //-------------------------------------------------------------------
      /* use least squared method */
      /* fit logR vs logM using 2 meshes */
      double pp, bb;
      leastSquaredMethod(2, hor, enc, &pp, &bb);
      Rg = pow(Renc / bb, 1.0 / pp);
      iRg = 0;
      aRg = 0.0;
      //-------------------------------------------------------------------
    }/* else{ */
    //---------------------------------------------------------------------
    /* set z using table search */
    const double zenc = UNIRAND_DBL;
    double azz, zz;
    int jzz;
    if( zenc <= rho[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, iRg, 0)] + DBL_EPSILON ){
      //-------------------------------------------------------------------
      /* use least squared method */
      /* fit logz vs logM using 2 meshes */
      double pp, bb;
      leastSquaredMethod(2, ver, &rho[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, iRg, 0)], &pp, &bb);
      zz = pow(zenc / bb, 1.0 / pp);
      jzz = 0;
      azz = 0.0;
      //-------------------------------------------------------------------
    }/* if( zenc < rho[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, iRg, 0)] ){ */
    else{
      //-------------------------------------------------------------------
      if( zenc <= (rho[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, iRg, NDISKBIN_VER - 1)] - DBL_EPSILON) ){
	/* use table search */
	jzz = bisection(zenc, NDISKBIN_VER, &rho[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, iRg, 0)], false, 1.0, &azz);
	azz /= (rho[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, iRg, 1 + jzz)] - rho[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, iRg, jzz)]);
	zz = (1.0 - azz) * ver[jzz] + azz * ver[1 + jzz];
      }/* if( zenc <= ([INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, iRg, NDISKBIN_VER)] - DBL_EPSILON) ){ */
      else{
	/* this means integrated mass along z-axis is zero */
	jzz = 0;
	azz = 0.0;
	zz = 0.0;
      }/* else{ */
      //-------------------------------------------------------------------
#if 0
      if( azz == 0.0 ){
	fprintf(stdout, "# ii = %zu, iRg = %d, ", ii - (*Nuse)
      }
#endif
      //-------------------------------------------------------------------
    }/* else{ */
    if( gsl_rng_get(GSLRand) < rand_half )
      zz *= -1.0;
    //---------------------------------------------------------------------
    __NOTE__("%zu-th particle: guiding center determined\n", ii - (*Nuse));
    //---------------------------------------------------------------------
    /* calculate physical quantities at R = Rg, z = zz */
    const double sigmaz = (1.0 - aRg) * sig[iRg] + aRg * sig[1 + iRg];
    const double diskpot =
      ((1.0 - azz) * pot[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER,     iRg, jzz)] + azz * pot[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER,     iRg, 1 + jzz)]) * (1.0 - aRg) +
      ((1.0 - azz) * pot[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1 + iRg, jzz)] + azz * pot[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1 + iRg, 1 + jzz)]) *        aRg;
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* calculate escape velocity at r^2 = Rg^2 + zz^2 */
    double delta;
    const int idx = findIdxSpherical(sqrt(Rg * Rg + zz * zz), prf, invlogrbin, &delta);
    const double sphepot = -((1.0 - delta) * prf[idx].psi_tot + delta * prf[1 + idx].psi_tot);
    const double PhiRz = diskpot + sphepot;
    const double vesc2 = -2.0 * PhiRz;
    //---------------------------------------------------------------------
    /* set a scale factor to include effects of disk thickness */
#ifdef  USE_POTENTIAL_SCALING_SCHEME
    double aR0;
    const int iR0 = findIdxSpherical(Rg, prf, invlogrbin, &aR0);
    const double PhiR0 =
      (1.0 - aRg) * pot[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, iRg, 0)] + aRg * pot[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1 + iRg, 0)] +
      (1.0 - aR0) * (-prf[iR0].psi_tot)                           + aR0 * (-prf[1 + iR0].psi_tot);
    double potScaling = PhiRz / PhiR0;
    if( potScaling > 1.0 )
      potScaling = 1.0;
#endif//USE_POTENTIAL_SCALING_SCHEME
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* calculate physical quantities at R = Rg, z = 0 */
    //---------------------------------------------------------------------
    /* determine sigmaR */
    /* calculate circular speed, epicyclic frequency, ... */
#ifndef USE_POTENTIAL_SCALING_SCHEME
    const double Omega2 =
      (((1.0 - azz) * dPhidR [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER,     iRg, jzz)] + azz *  dPhidR [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER,     iRg, 1 + jzz)]) * (1.0 - aRg) +
       ((1.0 - azz) * dPhidR [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1 + iRg, jzz)] + azz *  dPhidR [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1 + iRg, 1 + jzz)]) *        aRg) / Rg;
    const double d2Phi =
      ((1.0 - azz) * d2PhidR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER,     iRg, jzz)] + azz * d2PhidR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER,     iRg, 1 + jzz)]) * (1.0 - aRg) +
      ((1.0 - azz) * d2PhidR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1 + iRg, jzz)] + azz * d2PhidR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1 + iRg, 1 + jzz)]) *        aRg;
#else///USE_POTENTIAL_SCALING_SCHEME
    const double Omega2 = potScaling * ((1.0 - aRg) *  dPhidR [iRg] + aRg *  dPhidR [1 + iRg]) / Rg;
    const double d2Phi  = potScaling * ((1.0 - aRg) * d2PhidR2[iRg] + aRg * d2PhidR2[1 + iRg]);
#endif//USE_POTENTIAL_SCALING_SCHEME
    //---------------------------------------------------------------------
    double vcirc2 = Rg * Rg * Omega2;
    double vcirc  = sqrt(vcirc2);
    //---------------------------------------------------------------------
    const double gam2inv = 0.25 * (3.0 + d2Phi / (DBL_MIN + Omega2));
    const double gam2    = 1.0 / (DBL_MIN + gam2inv);
    const double sz2inv = 1.0 / (DBL_MIN + sigmaz * sigmaz);
#ifdef  USE_ORIGINAL_VDISP_ESTIMATOR
    const double sigmap = DISK_PERP_VDISP(sigmaz, vcirc, frac);
    const double sp2inv = 1.0 / (DBL_MIN + sigmap * sigmap);
    const double sR2inv = gam2inv * sp2inv;
    const double sigmaR = sqrt(gam2) * sigmap;
#else///USE_ORIGINAL_VDISP_ESTIMATOR
    const double sR2inv = 1.0 / DISK_RADIAL_VDISP2(sig[0] * sig[0], Rg, invRd);
    const double sp2inv = gam2 * sR2inv;
#ifdef  SPEEDUP_CONVERGENCE
    const double sigmaR = 1.0 / sqrt(DBL_MIN + sR2inv);
    const double sigmap = 1.0 / sqrt(DBL_MIN + sp2inv);
#endif//SPEEDUP_CONVERGENCE
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
    //---------------------------------------------------------------------
#if 0
    /* WARNING: THIS IS A TENTATIVE IMPLEMENTATION FOR DEBUGGING */
    if( vcirc2 > vesc2 ){
      //-------------------------------------------------------------------
      fprintf(stderr, "\tRg = %e\n", Rg);
      fprintf(stderr, "\tzz = %e\n", zz);
#ifdef  USE_POTENTIAL_SCALING_SCHEME
      fprintf(stderr, "\tpotScaling = %e\n", potScaling);
#endif//USE_POTENTIAL_SCALING_SCHEME
      fprintf(stderr, "\tOmega2 = %e\n", Omega2);
      fprintf(stderr, "\tkappa = %e\n", sqrt(d2Phi + 3.0 * Omega2));
      fprintf(stderr, "\tgam = %e\n", sqrt(gam2));
      fprintf(stderr, "\tPhiRz = %e\n", PhiRz);
      fprintf(stderr, "\tsphepot = %e\n", sphepot);
      fprintf(stderr, "\tvesc = %e\n", sqrt(vesc2));
      fprintf(stderr, "\tvcirc = %e\n", vcirc);
      fprintf(stderr, "\tsigmaR = %e\n", 1.0 / sqrt(sR2inv));
      fprintf(stderr, "\tsigmap = %e\n", 1.0 / sqrt(sp2inv));
      fprintf(stderr, "\tsigmaz = %e\n", sigmaz);
      //-------------------------------------------------------------------
      __KILL__(stderr, "ERROR: %zu-th particle: circular speed (%e) exceed escape velocity (%e)\n", ii - (*Nuse), vcirc, sqrt(vesc2));
      //-------------------------------------------------------------------
    }/* if( vcirc2 > vesc2 ){ */
#endif
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* determine particle position and velocity */
    double vx, vy, vz;
    double xx, yy;
    //---------------------------------------------------------------------
    const double vmax2 = vesc2 - vcirc2;
#ifndef  SPEEDUP_CONVERGENCE
    const double vmax  = sqrt(vmax2);
#endif//SPEEDUP_CONVERGENCE
    //---------------------------------------------------------------------
    double vR, vp;
    //---------------------------------------------------------------------
    while( true ){
      //-------------------------------------------------------------------
      while( true ){
	//-----------------------------------------------------------------
#ifdef  SPEEDUP_CONVERGENCE
	while( true ){
	  //---------------------------------------------------------------
	  vR = RANDVAL_DBL;
	  vp = RANDVAL_DBL;
	  vz = RANDVAL_DBL;
	  //---------------------------------------------------------------
	  if( vR * vR + vp * vp + vz * vz < 1.0 )
	    break;
	  //---------------------------------------------------------------
	}/* while( true ){ */
	/* 8x sigma --> DF = e-32 ~ 1.3e-14 is sufficiently small probability */
	vR *= 8.0 * sigmaR;
	vp *= 8.0 * sigmap;
	vz *= 8.0 * sigmaz;
#else///SPEEDUP_CONVERGENCE
      	vR = vmax * RANDVAL_DBL;
      	vp = vmax * RANDVAL_DBL;
      	vz = vmax * RANDVAL_DBL;
#endif//SPEEDUP_CONVERGENCE
	//-----------------------------------------------------------------
	if( vR * vR + vp * vp + vz * vz < vmax2 )
	  break;
	//-----------------------------------------------------------------
      }/* while( true ){ */
      //-------------------------------------------------------------------
      /* use Schwarzschild DF */
      const double val = exp(-0.5 * (vR * vR * sR2inv + vp * vp * sp2inv + vz * vz * sz2inv));
      const double try = UNIRAND_DBL;/* maximum of DF is unity */
      if( val > try )	break;
      //-------------------------------------------------------------------
    }/* while( true ){ */
    //---------------------------------------------------------------------
#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
    Krot += vcirc2;
    Krnd += vR * vR + vp * vp + vz * vz;
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION
    //---------------------------------------------------------------------
    /* uniform distribution in polar angle */
    vp += vcirc;
    const double phi = 2.0 * M_PI * UNIRAND_DBL;
    const double cosphi = cos(phi);
    const double sinphi = sin(phi);
    xx = Rg * cosphi;    vx = vR * cosphi - vp * sinphi;
    yy = Rg * sinphi;    vy = vR * sinphi + vp * cosphi;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    body.pos[ii].x = (real)xx;
    body.pos[ii].y = (real)yy;
    body.pos[ii].z = (real)zz;
    body.pos[ii].m = mass;
    //---------------------------------------------------------------------
    body.acc[ii].x = body.acc[ii].y = body.acc[ii].z = body.acc[ii].pot = ZERO;
    //---------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
    body.vel[ii].x = (real)vx;
    body.vel[ii].y = (real)vy;
    body.vel[ii].z = (real)vz;
#else///BLOCK_TIME_STEP
    body.vx[ii] = (real)vx;
    body.vy[ii] = (real)vy;
    body.vz[ii] = (real)vz;
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
    body.idx[ii] = ii;
    //---------------------------------------------------------------------
#if 0
    fprintf(stderr, "%e\t%e\t%e\t%e\t%e\t%e\n", xx, yy, zz, vx, vy, vz);
#endif
    //---------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
    if( (ii - (*Nuse)) == (stage * nunit) ){
      fprintf(stdout, "# ~%zu%% completed\n", stage * 10);
      fflush(stdout);
      stage++;
    }/* if( (ii - (*Nuse)) == (stage * nunit) ){ */
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  *Nuse += num;
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# finish distributing disk particles (%zu bodies)\n", num);
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------
#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
  fprintf(stdout, "# simple check based on Ostriker-Peebles criterion: Krot / (Krot + Krand) > ~0.28 is unstable to a bar-like mode\n");
  fprintf(stdout, "# \tw/o Ekin of spherical component: Krot / (Krot + Krand) = %e; i.e., %s to a bar-like mode\n", Krot / (Krot + Krnd), ((Krot / (Krot + Krnd)) > 0.28) ? "unstable" : "  stable");
  Krnd += disk.Krand_sph;
  fprintf(stdout, "# \tw/z Ekin of spherical component: Krot / (Krot + Krand) = %e; i.e., %s to a bar-like mode\n", Krot / (Krot + Krnd), ((Krot / (Krot + Krnd)) > 0.28) ? "unstable" : "  stable");
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "#\n#\n");
  fflush(stdout);
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
