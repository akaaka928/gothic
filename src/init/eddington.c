/*************************************************************************\
 *                                                                       *
                  last updated on 2016/08/12(Fri) 10:56:40
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
#include <macro.h>
#include <constants.h>
//-------------------------------------------------------------------------
#include "magi.h"
#include "profile.h"
#include "eddington.h"
//-------------------------------------------------------------------------
extern const real newton;
//-------------------------------------------------------------------------
extern double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
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
  }/* while( bisection ){ */
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
  }/* for(int kk = 0; kk < kind; kk++){ */
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
  }/* for(int ii = 2; ii < 2 + NRADBIN; ii++){ */
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
