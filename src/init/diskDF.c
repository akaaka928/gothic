/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/07(Wed) 16:13:52
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
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
//-------------------------------------------------------------------------
#include "macro.h"
#include "constants.h"
//-------------------------------------------------------------------------
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "profile.h"
#include "blas.h"
#include "spline.h"
#include "magi.h"
#include "potdens.h"
#include "diskDF.h"
//-------------------------------------------------------------------------
extern const real newton;
//-------------------------------------------------------------------------
extern double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
//-------------------------------------------------------------------------
#ifdef  USE_SFMT
#include "SFMT.h"
extern sfmt_t sfmt;
#ifdef  USE_SFMTJUMP
#pragma omp threadprivate(sfmt)
#endif//USE_SFMTJUMP
#define UNIRAND_DBL (sfmt_genrand_res53(&sfmt))
#else///USE_SFMT
#include <gsl/gsl_rng.h>
extern gsl_rng *GSLRand;
#define UNIRAND_DBL (gsl_rng_uniform(GSLRand))
#endif//USE_SFMT
#define UNIRAND     (CAST_D2R(UNIRAND_DBL))
#define RANDVAL_DBL (2.0 * (UNIRAND_DBL) - 1.0)
#define RANDVAL     (TWO * (UNIRAND    ) - UNITY)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef USE_SFMTJUMP
static inline
#else///USE_SFMTJUMP
int bisection(const double val, const int num, double * restrict tab, const bool logtbl, const double invbin, double * restrict ratio);
#endif//USE_SFMTJUMP
int bisection(const double val, const int num, double * restrict tab, const bool logtbl, const double invbin, double * restrict ratio)
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
#ifndef USE_SFMTJUMP
static inline
#else///USE_SFMTJUMP
int findIdxSpherical(const double rad, profile *prf, const double invlogbin, double *ratio);
#endif//USE_SFMTJUMP
int findIdxSpherical(const double rad, profile *prf, const double invlogbin, double *ratio)
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
double gaussQD4enc(const double min, const double max, double * restrict rad, const double invdr, double * restrict rho);
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
static inline void convertCylindrical2Spherical
(const int lev, const int ihead, const int itail, const int jhead, const int jtail,
 const double invrrbin, double * restrict rho, disk_data disk)
{
  //-----------------------------------------------------------------------
  const double fac  = 1.0 / (double)SUBDIVIDE_NUM;
  const double fac2 = fac * fac;
  //-----------------------------------------------------------------------
  const double dR = ldexp(disk.hh, -lev);
  const double dz = dR;
  //-----------------------------------------------------------------------
  for(int ii = ihead; ii < itail; ii++){
    //---------------------------------------------------------------------
    const double R0 = (0.5 + (double)ii) * dR;
    const double dV = R0 * dR * dz;
    const double Rmin = R0 - 0.5 * dR;
    //---------------------------------------------------------------------
    for(int jj = jhead; jj < jtail; jj++){
      //-------------------------------------------------------------------
      const double mass = (*disk.rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] * dV;
      const double zmin = (double)jj * dz;
      //-------------------------------------------------------------------
      for(int kk = 0; kk < SUBDIVIDE_NUM; kk++){
	//-----------------------------------------------------------------
	const double RR = Rmin + (0.5 + (double)kk) * dR * fac;
	const double R2 = RR * RR;
	//-----------------------------------------------------------------
	for(int ll = 0; ll < SUBDIVIDE_NUM; ll++){
	  //---------------------------------------------------------------
	  const double zz = zmin + (0.5 + (double)ll) * dz * fac;
	  const double r2 = R2 + zz * zz;
	  const int idx = (int)nearbyint(sqrt(r2) * invrrbin - 0.5);
	  //---------------------------------------------------------------
	  rho[idx] += fac2 * mass;
	  //---------------------------------------------------------------
	}/* for(int ll = 0; ll < SUBDIVIDE_NUM; ll++){ */
      }/* for(int kk = 0; kk < SUBDIVIDE_NUM; kk++){ */
      //-------------------------------------------------------------------
    }/* for(int jj = jhead; jj < jtail; jj++){ */
  }/* for(int ii = ihead; ii < itail; ii++){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void integrateSphericalDensityProfile(const int ndisk, const int maxLev, disk_data *disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int kk = 0; kk < ndisk; kk++){
    //---------------------------------------------------------------------
    /* load disk data */
    //---------------------------------------------------------------------
    const double Mdisk = disk[kk].cfg->Mtot;
    //---------------------------------------------------------------------
    double *rad;    rad = disk[kk].radSph;
    double *rho;    rho = disk[kk].rhoSph;
    double *enc;    enc = disk[kk].encSph;
    //---------------------------------------------------------------------
    profile *prf;    prf = disk[kk].prf;
    const double logrbin_sph = disk[kk].logrbin;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* calculate spherical averaged profiles of (volume-)density and enclosed mass */
    //---------------------------------------------------------------------
    /* initialization */
    const double    rrbin = sqrt(disk[kk].Rmax * disk[kk].Rmax + disk[kk].zmax * disk[kk].zmax) / (double)NDISKBIN_RAD;
    const double invrrbin = 1.0 / rrbin;
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_RAD; ii++){
      rad[ii] = (0.5 + (double)ii) * rrbin;
      rho[ii] = 0.0;
      enc[ii] = 0.0;
    }/* for(int ii = 0; ii < NDISKBIN_RAD; ii++){ */
    //---------------------------------------------------------------------
    /* density assignment */
    convertCylindrical2Spherical(maxLev - 1, 0, NDISKBIN_HOR >> 1, 0, NDISKBIN_VER >> 1, invrrbin, rho, disk[kk]);
    for(int ll = maxLev - 1; ll >= 0; ll--){
      convertCylindrical2Spherical(ll, 0, NDISKBIN_HOR >> 1, NDISKBIN_VER >> 1, NDISKBIN_VER, invrrbin, rho, disk[kk]);
      convertCylindrical2Spherical(ll, NDISKBIN_HOR >> 1, NDISKBIN_HOR, 0, NDISKBIN_VER, invrrbin, rho, disk[kk]);
    }/* for(int ll = maxLev - 1; ll >= 0; ll--){ */
    /* yield enclosed mass and volume density (by dividing by volume of pieces) */
    enc[0] = 4.0 * M_PI * rho[0];
    rho[0] *= 3.0 * invrrbin / (3.0 * rad[0] * rad[0] + 0.25 * rrbin * rrbin);
    for(int ii = 1; ii < NDISKBIN_RAD; ii++){
      enc[ii] = enc[ii - 1] + 4.0 * M_PI * rho[ii];
      rho[ii] *= 3.0 * invrrbin / (3.0 * rad[ii] * rad[ii] + 0.25 * rrbin * rrbin);
    }
    //---------------------------------------------------------------------
    /* adjust numerical errors */
    const double Mscale = disk[kk].cfg->Mtot / enc[NDISKBIN_RAD - 1];
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_RAD; ii++){
      rho[ii] *= Mscale;
      enc[ii] *= Mscale;
    }/* for(int ii = 0; ii < NDISKBIN_RAD; ii++){ */
#ifdef  PROGRESS_REPORT_ON
    fprintf(stdout, "# spherical averaged density profile correction: multiplying %e for %d-th disk\n", Mscale, kk);
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
    }/* for(int ii = 0; ii < 4; ii++){ */
    double power[2], basis[2];
    for(int ii = 0; ii < 2; ii++){
      power[ii] = (SS * Sxy[ii] - Sx * Sy[ii]) / (SS * Sxx - Sx * Sx);
      basis[ii] = pow(10.0, (Sy[ii] - power[ii] * Sx) / SS);
    }/* for(int ii = 0; ii < 2; ii++){ */
    /* extrapolate: enc, rho --> prf */
#pragma omp parallel for
    for(int ii = 0; ii < head; ii++){
      prf[ii].enc = basis[0] * pow(prf[ii].rad, power[0]);
      prf[ii].rho = basis[1] * pow(prf[ii].rad, power[1]);
    }/* for(int ii = 0; ii < head; ii++){ */
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
    }/* for(int ii = tail; ii < 4 + NRADBIN; ii++){ */
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
    }/* for(int ii = NRADBIN + 1; ii >= 0; ii--){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* multiply overall factor */
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < 4 + NRADBIN; ii++)
      prf[ii].psi *= CAST_R2D(newton);
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < ndisk; kk++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define SMOOTH(a, b, c, d) ((((a) * (b) > 0.0) && ((b) * (c) > 0.0) && ((c) * (d) > 0.0)) ? (0.5 * ((b) + (c))) : (0.0)) */
//-------------------------------------------------------------------------
/* #ifndef USE_POTENTIAL_SCALING_SCHEME */
/* //------------------------------------------------------------------------- */
/* static inline void diffResolutionBoundary(double * restrict ff, const int ii, const int jj, double * restrict fmm, double * restrict fmp, double * restrict fpm, double * restrict fpp) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   const int ic = ii >> 1;  const int im = ic - 1;  const int ip = ic + 1; */
/*   const int jc = jj >> 1;  const int jm = jc - 1;  const int jp = jc + 1; */
/*   //----------------------------------------------------------------------- */
/*   const double mm = ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, im, jc)] - 0.125 * (ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, im, jp)] - ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, im, jm)]); */
/*   const double cm = ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ic, jc)] - 0.125 * (ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ic, jp)] - ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ic, jm)]); */
/*   const double pm = ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ip, jc)] - 0.125 * (ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ip, jp)] - ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ip, jm)]); */
/*   //----------------------------------------------------------------------- */
/*   const double mp = ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, im, jc)] + 0.125 * (ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, im, jp)] - ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, im, jm)]); */
/*   const double cp = ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ic, jc)] + 0.125 * (ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ic, jp)] - ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ic, jm)]); */
/*   const double pp = ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ip, jc)] + 0.125 * (ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ip, jp)] - ff[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ip, jm)]); */
/*   //----------------------------------------------------------------------- */
/*   *fmm = cm - 0.125 * (mp - mm);/\* ii = i-, jj = j- *\/ */
/*   *fmp = cm + 0.125 * (mp - mm);/\* ii = i-, jj = j+ *\/ */
/*   *fpm = cp - 0.125 * (pp - pm);/\* ii = i+, jj = j- *\/ */
/*   *fpp = cp + 0.125 * (pp - pm);/\* ii = i+, jj = j+ *\/ */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
/* #endif//USE_POTENTIAL_SCALING_SCHEME */
//-------------------------------------------------------------------------
void diffAxisymmetricPotential(const int maxLev, const disk_data disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int lev = 0; lev < maxLev; lev++){
    //---------------------------------------------------------------------
    double *RR;    RR = &(disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)]);
    double *zz;    zz = &(disk.ver[INDEX2D(maxLev, NDISKBIN_VER, lev, 0)]);
    double *  Phi    ;      Phi     = &(disk.  pot   [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, 0)]);
#ifndef USE_POTENTIAL_SCALING_SCHEME
    double * dPhi_dR ;	   dPhi_dR  = &(disk. dPhidR [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, 0)]);
    double *d2Phi_dR2;	  d2Phi_dR2 = &(disk.d2PhidR2[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, 0)]);
    /* double *f1, *f2; */
    /* if( lev > 0 ){ */
    /*   f1 = &(disk. dPhidR [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev - 1, 0)]); */
    /*   f2 = &(disk.d2PhidR2[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev - 1, 0)]); */
    /* }/\* if( lev > 0 ){ *\/ */
#else///USE_POTENTIAL_SCALING_SCHEME
    double * dPhi_dR ;	   dPhi_dR  = &(disk. dPhidR [INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]);
    double *d2Phi_dR2;	  d2Phi_dR2 = &(disk.d2PhidR2[INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]);
#endif//USE_POTENTIAL_SCALING_SCHEME
    const double invRbin = 1.0 / ldexp(disk.hh, -lev);
    profile *prf;    prf = disk.prf;
    const double invlogrbin = disk.invlogrbin;
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR - 1; ii++){
      //-------------------------------------------------------------------
#ifdef  USE_POTENTIAL_SCALING_SCHEME
      const int jj = 0;
#else///USE_POTENTIAL_SCALING_SCHEME
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
#endif//USE_POTENTIAL_SCALING_SCHEME
	{
	  //---------------------------------------------------------------
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
	  //---------------------------------------------------------------
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
	  const double  dPhidr__sphe = CAST_R2D(newton) *                           enc * rinv * rinv;
	  const double d2Phidr2_sphe = CAST_R2D(newton) * (4.0 * M_PI * rho - 2.0 * enc * rinv * rinv * rinv);
	  //---------------------------------------------------------------
	  /* R-derivatives of potential given by the all components */
#ifndef USE_POTENTIAL_SCALING_SCHEME
	  dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] = _dPhidR__disk + dPhidr__sphe * drdR;
	  d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] = d2PhidR2_disk + dPhidr__sphe * drdR * d2rdR2 + d2Phidr2_sphe * drdR * drdR;
#else///USE_POTENTIAL_SCALING_SCHEME
	  dPhi_dR  [ii] = _dPhidR__disk + dPhidr__sphe * drdR;
	  d2Phi_dR2[ii] = d2PhidR2_disk + dPhidr__sphe * drdR * d2rdR2 + d2Phidr2_sphe * drdR * drdR;
#endif//USE_POTENTIAL_SCALING_SCHEME
	  //---------------------------------------------------------------
	}
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* modify values at the boundary of each resolution */
      //-------------------------------------------------------------------
/* #ifndef USE_POTENTIAL_SCALING_SCHEME */
/*       //------------------------------------------------------------------- */
/*       if( lev > 0 ){ */
/* 	//----------------------------------------------------------------- */
/* 	double mm1, mp1, pm1, pp1;	diffResolutionBoundary(f1, ii, NDISKBIN_VER - 1, &mm1, &mp1, &pm1, &pp1); */
/* 	double mm2, mp2, pm2, pp2;	diffResolutionBoundary(f2, ii, NDISKBIN_VER - 1, &mm2, &mp2, &pm2, &pp2); */
/* 	//----------------------------------------------------------------- */
/* 	if( ii & 1 ){ */
/* 	  //--------------------------------------------------------------- */
/* 	  /\* upper half *\/ */
/* 	  //--------------------------------------------------------------- */
/* 	  dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 2)] = pm1; */
/* 	  dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)] = pp1; */
/* 	  d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 2)] = pm2; */
/* 	  d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)] = pp2; */
/* 	  //--------------------------------------------------------------- */
/* 	}/\* if( ii & 1 ){ *\/ */
/* 	else{ */
/* 	  //--------------------------------------------------------------- */
/* 	  /\* lower half *\/ */
/* 	  //--------------------------------------------------------------- */
/* 	  dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 2)] = mm1; */
/* 	  dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)] = mp1; */
/* 	  d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 2)] = mm2; */
/* 	  d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)] = mp2; */
/* 	  //--------------------------------------------------------------- */
/* 	}/\* else{ *\/ */
/* 	//----------------------------------------------------------------- */
/*       }/\* if( lev > 0 ){ *\/ */
/*       //------------------------------------------------------------------- */
/* #endif//USE_POTENTIAL_SCALING_SCHEME */
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NDISKBIN_HOR - 1; ii++){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
#ifndef USE_POTENTIAL_SCALING_SCHEME
    //---------------------------------------------------------------------
    /* if( lev > 0 ){ */
    /*   //------------------------------------------------------------------- */
    /*   for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
    /* 	//----------------------------------------------------------------- */
    /* 	double mm1, mp1, pm1, pp1;	diffResolutionBoundary(f1, NDISKBIN_HOR - 1, jj, &mm1, &mp1, &pm1, &pp1); */
    /* 	double mm2, mp2, pm2, pp2;	diffResolutionBoundary(f2, NDISKBIN_HOR - 1, jj, &mm2, &mp2, &pm2, &pp2); */
    /* 	//----------------------------------------------------------------- */
    /* 	if( jj & 1 ){ */
    /* 	  //--------------------------------------------------------------- */
    /* 	  /\* upper half *\/ */
    /* 	  //--------------------------------------------------------------- */
    /* 	  dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj)] = mp1; */
    /* 	  dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] = pp1; */
    /* 	  d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj)] = mp2; */
    /* 	  d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] = pp2; */
    /* 	  //--------------------------------------------------------------- */
    /* 	}/\* if( jj & 1 ){ *\/ */
    /* 	//----------------------------------------------------------------- */
    /* 	else{ */
    /* 	  //--------------------------------------------------------------- */
    /* 	  /\* lower half *\/ */
    /* 	  //--------------------------------------------------------------- */
    /* 	  dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj)] = mm1; */
    /* 	  dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] = pm1; */
    /* 	  d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj)] = mm2; */
    /* 	  d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] = pm2; */
    /* 	  //--------------------------------------------------------------- */
    /* 	}/\* else{ *\/ */
    /* 	//----------------------------------------------------------------- */
    /*   }/\* for(int jj = 0; jj < NDISKBIN_VER; jj++){ *\/ */
    /*   //------------------------------------------------------------------- */
    /* }/\* if( lev > 0 ){ *\/ */
    /* //--------------------------------------------------------------------- */
    /* else */
      for(int jj = 0; jj < NDISKBIN_VER; jj++){
	//-----------------------------------------------------------------
	dPhi_dR  [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] =  dPhi_dR [INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj)];
	d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] = d2Phi_dR2[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj)];
	//-----------------------------------------------------------------
      }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
    //---------------------------------------------------------------------
#else///USE_POTENTIAL_SCALING_SCHEME
    //---------------------------------------------------------------------
    dPhi_dR  [NDISKBIN_HOR - 1] =  dPhi_dR [NDISKBIN_HOR - 2];
    d2Phi_dR2[NDISKBIN_HOR - 1] = d2Phi_dR2[NDISKBIN_HOR - 2];
    //---------------------------------------------------------------------
#endif//USE_POTENTIAL_SCALING_SCHEME
    //---------------------------------------------------------------------
  }/* for(int lev = 0; lev < maxLev; lev++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void calcVerticalVdisp(const int ndisk, const int maxLev, disk_data *disk_info)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int kk = 0; kk < ndisk; kk++){
    //---------------------------------------------------------------------
    /* load disk data */
    //---------------------------------------------------------------------
    const disk_data disk = disk_info[kk];
    const double invzbin0 = 1.0 / disk.hh;
    //---------------------------------------------------------------------
    profile *sph;    sph = disk.prf;
    const double invlogrbin = disk.invlogrbin;
    //---------------------------------------------------------------------
    for(int lev = 0; lev < maxLev; lev++){
      //---------------------------------------------------------------------
      double *sig;      sig = &(disk.sigmaz[INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]);
      double *Phi;      Phi = &(disk.pot   [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, 0)]);
      double *hor;      hor = &(disk.hor   [INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]);
      //-------------------------------------------------------------------
      const double invzbin = 1.0 / ldexp(disk.hh, -lev);
      //-------------------------------------------------------------------
#ifndef ENABLE_VARIABLE_SCALE_HEIGHT
      const double zd = disk.cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* calculate vertical velocity dispersion */
      //-------------------------------------------------------------------
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	//-----------------------------------------------------------------
	int i0;
	double a0, Psi_s, Psi_d;
	//-----------------------------------------------------------------
	/* get potential @ (R, z = 0) */
	const double RR = hor[ii];
	i0 = findIdxSpherical(RR, sph, invlogrbin, &a0);
	Psi_s = (1.0 - a0) * sph[i0].psi_tot + a0 * sph[i0 + 1].psi_tot;
	Psi_d = -Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 0)];
	const double Psi_R_0 = Psi_d + Psi_s;
	//-----------------------------------------------------------------
	/* get potential @ (R, z = zd) */
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
	/* double zd = disk.cfg->zd * modulateDiskThickness(RR, disk); */
	const double zd = disk.zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
	const double rad = sqrt(RR * RR + zd * zd);
	i0 = findIdxSpherical(rad, sph, invlogrbin, &a0);
	Psi_s = (1.0 - a0) * sph[i0].psi_tot + a0 * sph[i0 + 1].psi_tot;
	/* get i0 in the appropriate level of grid points (must be coarser than the current level) */
	i0 = bisection(zd, NDISKBIN_VER, disk.ver, false, invzbin0, &a0);
	int ll = (i0 != 0) ? ((int)ilog2((uint)(NDISKBIN_VER / i0))) : (lev);
	if( ll > lev )	  ll = lev;
	const int i1 = ii >> (lev - ll);
	i0 = bisection(zd, NDISKBIN_VER, &(disk.ver[INDEX2D(maxLev, NDISKBIN_VER, ll, 0)]), false, invzbin, &a0);
	Psi_d = (1.0 - a0) * (-Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, i1, i0)]) + a0 * (-Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, i1, 1 + i0)]);
	const double Psi_R_zd = Psi_d + Psi_s;
	//-----------------------------------------------------------------
	sig[ii] = sqrt(Psi_R_0 - Psi_R_zd);
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
      //-------------------------------------------------------------------
    }/* for(int lev = 0; lev < maxLev; lev++){ */
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < ndisk; kk++){ */
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
#ifndef USE_SFMTJUMP
static inline
#else///USE_SFMTJUMP
int bisection4nestedGrid
(const double val, const int num, double * restrict tab, const bool logtbl, const double invbin, double * restrict ratio, const int maxLev, int * restrict lev, double * restrict tab_lev);
#endif//USE_SFMTJUMP
int bisection4nestedGrid
(const double val, const int num, double * restrict tab, const bool logtbl, const double invbin, double * restrict ratio, const int maxLev, int * restrict lev, double * restrict tab_lev)
{
  //-----------------------------------------------------------------------
  /* 1. find the nested level */
  *lev = 0;
  for(int ii = maxLev - 1; ii >= 0; ii--)
    if( val <= tab_lev[ii] ){
      *lev = ii;
      break;
    }/* if( val <= tab_lev[ii] ){ */
  //-----------------------------------------------------------------------
  /* 2. find the index */
  int idx = bisection(val, num, &(tab[INDEX2D(maxLev, num, *lev, 0)]), logtbl, invbin, ratio);
  //-----------------------------------------------------------------------
  return (idx);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void distributeDiskParticles(ulong *Nuse, iparticle body, const real mass, const int maxLev, const disk_data disk)
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
  double *node_hor;  node_hor = disk.node_hor;
  double *node_ver;  node_ver = disk.node_ver;
  double *zone_hor;  zone_hor = disk.hor;
  double *zone_ver;  zone_ver = disk.ver;
  double *pot;  pot =  disk.pot;
  double *enc;  enc =  disk.enc;
  double *rho;  rho = *disk.rhoSum;
  double *sig;  sig =  disk.sigmaz;
  //-----------------------------------------------------------------------
  double * dPhidR ;  dPhidR   =  disk. dPhidR;
  double *d2PhidR2;  d2PhidR2 =  disk.d2PhidR2;
  //-----------------------------------------------------------------------
  const ulong num = disk.cfg->num;
#ifdef  USE_ORIGINAL_VDISP_ESTIMATOR
  const double frac = disk.cfg->vdisp_frac;
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
  const double Mmax = enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, 0, NDISKBIN_HOR)];
  //-----------------------------------------------------------------------
  /* const ulong rand_half = ((gsl_rng_min(GSLRand) + gsl_rng_max(GSLRand)) >> 1); */
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
#ifndef USE_SFMTJUMP
  const ulong nunit = (ulong)ceilf(0.1f * (float)num);
  ulong stage = 1;
#else///USE_SFMTJUMP
#pragma omp single nowait
#endif//USE_SFMTJUMP
  {
    fprintf(stdout, "#\n#\n# start distributing disk particles (%zu bodies: [%zu:%zu])\n", num, *Nuse, (*Nuse) + num - 1);
    fflush(stdout);
  }
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------
#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
  double Krot = 0.0;
  double Krnd = 0.0;
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION
  //-----------------------------------------------------------------------
  double *tab_lev;
  tab_lev = (double *)malloc(sizeof(double) * maxLev);
  if( tab_lev == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tab_lev\n");  }
  for(int lev = 0; lev < maxLev; lev++)
    tab_lev[lev] = enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, NDISKBIN_HOR)];
  //-----------------------------------------------------------------------
#ifdef  USE_SFMTJUMP
#pragma omp for
#endif//USE_SFMTJUMP
  for(ulong ii = *Nuse; ii < (*Nuse) + num; ii++){
    //---------------------------------------------------------------------
    /* determine Rg, phi, and z */
    //---------------------------------------------------------------------
    /* set Rg */
    const double Renc = Mmax * UNIRAND_DBL;
    double aRg, Rg;
    int iRg, lev;
    //---------------------------------------------------------------------
    /* use table search */
    iRg = bisection4nestedGrid(Renc, NDISKBIN_HOR + 1, enc, false, 1.0, &aRg, maxLev, &lev, tab_lev);
    aRg /= (enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, 1 + iRg)] - enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, iRg)]);
    Rg = (1.0 - aRg) * node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, iRg)] + aRg * node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, 1 + iRg)];
    /* correction about ``iRg'' and ``aRg'' for zone centeric coordinate from edge coordinate */
    if( Rg >= zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)] ){
      if( Rg <= zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, NDISKBIN_HOR - 1)] ){
	iRg = bisection(Rg, NDISKBIN_HOR, &(zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)]), false, 1.0, &aRg);
	aRg /= (zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + iRg)] - zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, iRg)]);
      }/* if( Rg <= zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, NDISKBIN_HOR - 1)] ){ */
      else{
	iRg = NDISKBIN_HOR - 2;
	aRg = 1.0;
      }/* else{ */
    }/* if( Rg >= zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)] ){ */
    else{
      iRg = 0;
      aRg = 0.0;
    }/* else{ */
    //---------------------------------------------------------------------
    /* calculate vertical velocity dispersion at R = Rg */
    const double sigmaz = (1.0 - aRg) * sig[INDEX2D(maxLev, NDISKBIN_HOR, lev, iRg)] + aRg * sig[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + iRg)];
    //---------------------------------------------------------------------
    /* set z using table search */
    const double zenc = UNIRAND_DBL;
    int jzz = 0;
    double azz = 0.0;
    double zz = 0.0;
    //---------------------------------------------------------------------
    while( true ){
      //-------------------------------------------------------------------
      if( zenc > (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev, iRg, NDISKBIN_VER)] + DBL_EPSILON) ){
	//-----------------------------------------------------------------
	lev--;
	if( lev < 0 ){	  lev = 0;	  break;	}
	/* reset iRg and aRg in coarser grid */
	iRg = bisection(Renc, NDISKBIN_HOR + 1, &(enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, 0)]), false, 1.0, &aRg);
	aRg /= (enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, 1 + iRg)] - enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, iRg)]);
	/* correction about ``iRg'' and ``aRg'' for zone centeric coordinate from edge coordinate */
	if( Rg >= zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)] ){
	  if( Rg <= zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, NDISKBIN_HOR - 1)] ){
	    iRg = bisection(Rg, NDISKBIN_HOR, &(zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)]), false, 1.0, &aRg);
	    aRg /= (zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + iRg)] - zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, iRg)]);
	  }/* if( Rg <= zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, NDISKBIN_HOR - 1)] ){ */
	  else{
	    iRg = NDISKBIN_HOR - 2;
	    aRg = 1.0;
	  }/* else{ */
	}/* if( Rg >= zone_hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)] ){ */
	else{
	  iRg = 0;
	  aRg = 0.0;
	}/* else{ */
	//-----------------------------------------------------------------
      }
      else
	break;
    }/* while( true ){ */
    //---------------------------------------------------------------------
    /* use table search */
    jzz = bisection(zenc, NDISKBIN_VER + 1, &rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev, iRg, 0)], false, 1.0, &azz);
    azz /= (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev, iRg, 1 + jzz)] - rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev, iRg, jzz)]);
    zz = (1.0 - azz) * node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev, jzz)] + azz * node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev, 1 + jzz)];
    /* correction about ``jzz'' and ``azz'' for zone centeric coordinate from edge coordinate */
    if( zz >= zone_ver[INDEX2D(maxLev, NDISKBIN_VER, lev, 0)] ){
      if( zz <= zone_ver[INDEX2D(maxLev, NDISKBIN_VER, lev, NDISKBIN_VER - 1)] ){
	jzz = bisection(zz, NDISKBIN_VER, &(zone_ver[INDEX2D(maxLev, NDISKBIN_VER, lev, 0)]), false, 1.0, &azz);
	azz /= (zone_ver[INDEX2D(maxLev, NDISKBIN_VER, lev, 1 + jzz)] - zone_ver[INDEX2D(maxLev, NDISKBIN_VER, lev, jzz)]);
      }/* if( zz <= zone_ver[INDEX2D(maxLev, NDISKBIN_VER, lev, NDISKBIN_VER - 1)] ){ */
      else{
	jzz = NDISKBIN_VER - 2;
	azz = 1.0;
      }/* else{ */
    }/* if( zz >= zone_ver[INDEX2D(maxLev, NDISKBIN_VER, lev, 0)] ){ */
    else{
      jzz = 0;
      azz = 0.0;
    }/* else{ */
    //---------------------------------------------------------------------
    if( UNIRAND_DBL < 0.5 )
      zz *= -1.0;
    //---------------------------------------------------------------------
    __NOTE__("%zu-th particle: guiding center determined\n", ii - (*Nuse));
    //---------------------------------------------------------------------
    /* calculate potential at R = Rg, z = zz */
    const double diskpot =
      ((1.0 - azz) * pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev,     iRg, jzz)] + azz * pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev,     iRg, 1 + jzz)]) * (1.0 - aRg) +
      ((1.0 - azz) * pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + iRg, jzz)] + azz * pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + iRg, 1 + jzz)]) *	       aRg;
    //---------------------------------------------------------------------
#if 0
    fprintf(stderr, "%d\t%e\t%e\n", lev, Rg, zz);
    fflush(stderr);
#endif
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
      (1.0 - aRg) * pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, iRg, 0)] + aRg * pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + iRg, 0)] +
      (1.0 - aR0) * (-prf[iR0].psi_tot)                                         + aR0 * (-prf[1 + iR0].psi_tot);
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
      (((1.0 - azz) * dPhidR [INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev,     iRg, jzz)] + azz *  dPhidR [INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev,     iRg, 1 + jzz)]) * (1.0 - aRg) +
       ((1.0 - azz) * dPhidR [INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + iRg, jzz)] + azz *  dPhidR [INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + iRg, 1 + jzz)]) *        aRg) / Rg;
    const double d2Phi =
      ((1.0 - azz) * d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev,     iRg, jzz)] + azz * d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev,     iRg, 1 + jzz)]) * (1.0 - aRg) +
      ((1.0 - azz) * d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + iRg, jzz)] + azz * d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + iRg, 1 + jzz)]) *        aRg;
#else///USE_POTENTIAL_SCALING_SCHEME
    const double Omega2 = potScaling * ((1.0 - aRg) *  dPhidR [INDEX2D(maxLev, NDISKBIN_HOR, lev, iRg)] + aRg *  dPhidR [INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + iRg)]) / Rg;
    const double d2Phi  = potScaling * ((1.0 - aRg) * d2PhidR2[INDEX2D(maxLev, NDISKBIN_HOR, lev, iRg)] + aRg * d2PhidR2[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + iRg)]);
#endif//USE_POTENTIAL_SCALING_SCHEME
    //---------------------------------------------------------------------
    double vcirc2 = Rg * Rg * Omega2;
    double vcirc  = sqrt(vcirc2);
    //---------------------------------------------------------------------
    const double gam2inv = 0.25 * (3.0 + d2Phi / (DBL_MIN + Omega2));
    const double gam2    = 1.0 / (DBL_MIN + gam2inv);
    const double sz2inv  = 1.0 / (DBL_MIN + sigmaz * sigmaz);
#ifdef  USE_ORIGINAL_VDISP_ESTIMATOR
    const double sigmap = DISK_PERP_VDISP(sigmaz, vcirc, frac);
    const double sp2inv = 1.0 / (DBL_MIN + sigmap * sigmap);
    const double sR2inv = gam2inv * sp2inv;
    const double sigmaR = sqrt(gam2) * sigmap;
#else///USE_ORIGINAL_VDISP_ESTIMATOR
    const double sR2inv = 1.0 / DISK_RADIAL_VDISP2(sig[INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, 0)] * sig[INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, 0)], Rg, disk.invRd);
    const double sp2inv = gam2 * sR2inv;
#ifdef  SPEEDUP_CONVERGENCE
    const double sigmaR = 1.0 / sqrt(DBL_MIN + sR2inv);
    const double sigmap = 1.0 / sqrt(DBL_MIN + sp2inv);
#endif//SPEEDUP_CONVERGENCE
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
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
    body.pos[ii].x = CAST_D2R(xx);
    body.pos[ii].y = CAST_D2R(yy);
    body.pos[ii].z = CAST_D2R(zz);
    body.pos[ii].m = mass;
    //---------------------------------------------------------------------
    body.acc[ii].x = body.acc[ii].y = body.acc[ii].z = body.acc[ii].pot = ZERO;
    //---------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
    body.vel[ii].x = CAST_D2R(vx);
    body.vel[ii].y = CAST_D2R(vy);
    body.vel[ii].z = CAST_D2R(vz);
#else///BLOCK_TIME_STEP
    body.vx[ii] = CAST_D2R(vx);
    body.vy[ii] = CAST_D2R(vy);
    body.vz[ii] = CAST_D2R(vz);
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
    body.idx[ii] = ii;
    //---------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
#ifndef USE_SFMTJUMP
    if( (ii - (*Nuse)) == (stage * nunit) ){
      fprintf(stdout, "# ~%zu%% completed\n", stage * 10);
      fflush(stdout);
      stage++;
    }/* if( (ii - (*Nuse)) == (stage * nunit) ){ */
#endif//USE_SFMTJUMP
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  *Nuse += num;
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
#ifdef  USE_SFMTJUMP
#pragma omp barrier
#pragma omp master
#endif//USE_SFMTJUMP
  fprintf(stdout, "# finish distributing disk particles (%zu bodies)\n", num);
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------
#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
#ifdef  USE_SFMTJUMP
#pragma omp single nowait
#endif//USE_SFMTJUMP
  {
    fprintf(stdout, "# simple check based on Ostriker-Peebles criterion: Krot / (Krot + Krand) > ~0.28 is unstable to a bar-like mode\n");
    fprintf(stdout, "# \tw/o Ekin of spherical component: Krot / (Krot + Krand) = %e; i.e., %s to a bar-like mode\n", Krot / (Krot + Krnd), ((Krot / (Krot + Krnd)) > 0.28) ? "unstable" : "  stable");
    Krnd += disk.Krand_sph;
    fprintf(stdout, "# \tw/z Ekin of spherical component: Krot / (Krot + Krand) = %e; i.e., %s to a bar-like mode\n", Krot / (Krot + Krnd), ((Krot / (Krot + Krnd)) > 0.28) ? "unstable" : "  stable");
  }
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION
  //-----------------------------------------------------------------------
  free(tab_lev);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
