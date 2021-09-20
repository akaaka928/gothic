/**
 * @file diskDF.c
 *
 * @brief Source code for generating initial condition of disk component(s)
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/04/23 (Mon)
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
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "macro.h"
#include "constants.h"
#include "rand.h"

#include "../misc/structure.h"

#include "profile.h"
#include "blas.h"
#include "spline.h"
#include "magi.h"
#include "potdens.h"
#include "diskDF.h"


extern const real newton;


/**
 * @fn bisection
 *
 * @brief Execute bisection.
 *
 * @param (val) the target value
 * @param (num) number of data points
 * @param (tab) array contains data points
 * @param (logrbl) true when data points are sampled in logarithmic space
 * @param (invbin) inverse of interval between grid points
 * @return (ratio) parameter for linear interpolation
 * @return lower index of the corresponding data point
 */
static inline int bisection(const double val, const int num, double * restrict tab, const bool logtbl, const double invbin, double * restrict ratio)
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
      *ratio = (logtbl ? (log10(val / tab[ll])) : (val - tab[ll])) * invbin;
      return (ll);
    }/* if( (1 + ll) == rr ){ */
  }/* while( true ){ */
}


/**
 * @fn findIdxSpherical
 *
 * @brief Find a data element in the given array corresponding to the given value.
 *
 * @param (rad) radius
 * @param (prf) radial profile of the component
 * @param (invlogbin) inverse of logarithmic interval between grid points
 * @return (ratio) parameter for linear interpolation
 * @return the corresponding index to the given radius
 */
static inline int findIdxSpherical(const double rad, profile *prf, const double invlogbin, double *ratio)
{
  int ll =           0;
  int rr = NRADBIN - 1;

  if( rad < prf[ll].rad + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( rad > prf[rr].rad - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }

  while( true ){
    const uint cc = ((uint)ll + (uint)rr) >> 1;

    if( (prf[cc].rad - rad) * (prf[ll].rad - rad) <= 0.0 )      rr = (int)cc;
    else                                                        ll = (int)cc;

    if( (1 + ll) == rr ){
      *ratio = log10(rad / prf[ll].rad) * invlogbin;
      return (ll);
    }/* if( (1 + ll) == rr ){ */
  }/* while( true ){ */
}


/**
 * @fn convertCylindrical2Spherical
 *
 * @brief Convert mass profile in cylindrical coordinates to that in polar coordinates.
 *
 * @param (lev) nested grid level
 * @param (ihead) head index for i-loop
 * @param (itail) tail index for i-loop
 * @param (jhead) head index for j-loop
 * @param (jtail) tail index for j-loop
 * @param (invrrbin) inverse of interval between grid points in polar coordinates
 * @return (rho) spherical averaged density
 * @param (disk) physical properties of the disk component
 */
static inline void convertCylindrical2Spherical
(const int lev, const int ihead, const int itail, const int jhead, const int jtail,
 const double invrrbin, double * restrict rho, disk_data disk)
{
  const double fac  = 1.0 / (double)SUBDIVIDE_NUM;
  const double fac2 = fac * fac;

  const double dR = ldexp(disk.hh, -lev);
  const double dz = dR;

  for(int ii = ihead; ii < itail; ii++){
    const double R0 = (0.5 + (double)ii) * dR;
    const double dV = R0 * dR * dz;
    const double Rmin = R0 - 0.5 * dR;

    for(int jj = jhead; jj < jtail; jj++){
      const double mass = (*disk.rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] * dV;
      const double zmin = (double)jj * dz;

      for(int kk = 0; kk < SUBDIVIDE_NUM; kk++){
	const double RR = Rmin + (0.5 + (double)kk) * dR * fac;
	const double R2 = RR * RR;

	for(int ll = 0; ll < SUBDIVIDE_NUM; ll++){
	  const double zz = zmin + (0.5 + (double)ll) * dz * fac;
	  const double r2 = R2 + zz * zz;
	  const int idx = (int)nearbyint(sqrt(r2) * invrrbin - 0.5);

	  rho[idx] += fac2 * mass;
	}/* for(int ll = 0; ll < SUBDIVIDE_NUM; ll++){ */
      }/* for(int kk = 0; kk < SUBDIVIDE_NUM; kk++){ */
    }/* for(int jj = jhead; jj < jtail; jj++){ */
  }/* for(int ii = ihead; ii < itail; ii++){ */
}


/**
 * @fn integrateSphericalDensityProfile
 *
 * @brief Integrate spherical averaged density profile.
 *
 * @param (ndisk) number of disk components
 * @param (maxLev) maximum level of nested grid
 * @return (disk) physical properties of the disk component
 *
 * @sa convertCylindrical2Spherical
 * @sa findIdxSpherical
 * @sa bisection
 */
void integrateSphericalDensityProfile(const int ndisk, const int maxLev, disk_data *disk)
{
  __NOTE__("%s\n", "start");


  for(int kk = 0; kk < ndisk; kk++){
    /** load disk data */
    const double Mdisk = disk[kk].cfg->Mtot;
    double *rad;    rad = disk[kk].radSph;
    double *rho;    rho = disk[kk].rhoSph;
    double *enc;    enc = disk[kk].encSph;
    profile *prf;    prf = disk[kk].prf;
    const double logrbin_sph = disk[kk].logrbin;

    /** calculate spherical averaged profiles of (volume-)density and enclosed mass */
    /** initialization */
    const double    rrbin = sqrt(disk[kk].Rmax * disk[kk].Rmax + disk[kk].zmax * disk[kk].zmax) / (double)NDISKBIN_RAD;
    const double invrrbin = 1.0 / rrbin;
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_RAD; ii++){
      rad[ii] = (0.5 + (double)ii) * rrbin;
      rho[ii] = 0.0;
      enc[ii] = 0.0;
    }/* for(int ii = 0; ii < NDISKBIN_RAD; ii++){ */

    /** density assignment */
    convertCylindrical2Spherical(maxLev - 1, 0, NDISKBIN_HOR >> 1, 0, NDISKBIN_VER >> 1, invrrbin, rho, disk[kk]);
    for(int ll = maxLev - 1; ll >= 0; ll--){
      convertCylindrical2Spherical(ll, 0, NDISKBIN_HOR >> 1, NDISKBIN_VER >> 1, NDISKBIN_VER, invrrbin, rho, disk[kk]);
      convertCylindrical2Spherical(ll, NDISKBIN_HOR >> 1, NDISKBIN_HOR, 0, NDISKBIN_VER, invrrbin, rho, disk[kk]);
    }/* for(int ll = maxLev - 1; ll >= 0; ll--){ */

    /** yield enclosed mass and volume density (by dividing by volume of pieces) */
    enc[0] = 4.0 * M_PI * rho[0];
    rho[0] *= 3.0 * invrrbin / (3.0 * rad[0] * rad[0] + 0.25 * rrbin * rrbin);
    for(int ii = 1; ii < NDISKBIN_RAD; ii++){
      enc[ii] = enc[ii - 1] + 4.0 * M_PI * rho[ii];
      rho[ii] *= 3.0 * invrrbin / (3.0 * rad[ii] * rad[ii] + 0.25 * rrbin * rrbin);
    }/* for(int ii = 1; ii < NDISKBIN_RAD; ii++){ */

    /** adjust numerical errors */
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


    /** assign spherical averaged profile of density, enclosed mass and potential (internal part) */
    double dummy;
    const int head = findIdxSpherical(rad[               0], prf, 1.0, &dummy);
    const int tail = findIdxSpherical(rad[NDISKBIN_RAD - 1], prf, 1.0, &dummy);

    /** fit logr vs logM using 4 meshes */
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

    /** extrapolate: enc, rho --> prf */
#pragma omp parallel for
    for(int ii = 0; ii < head; ii++){
      prf[ii].enc = basis[0] * pow(prf[ii].rad, power[0]);
      prf[ii].rho = basis[1] * pow(prf[ii].rad, power[1]);
    }/* for(int ii = 0; ii < head; ii++){ */

    /** interpolate: enc, rho --> prf */
#pragma omp parallel for
    for(int ii = head; ii < tail; ii++){
      double alpha;
      const int idx = bisection(prf[ii].rad, NDISKBIN_RAD, rad, false, invrrbin, &alpha);
      prf[ii].rho = (1.0 - alpha) * rho[idx] + alpha * rho[1 + idx];
      prf[ii].enc = (1.0 - alpha) * enc[idx] + alpha * enc[1 + idx];
    }/* for(int ii = head; ii < tail; ii++){ */

    /** fill the total enclosed mass in the outer region */
#pragma omp parallel for
    for(int ii = tail; ii < NRADBIN; ii++){
      prf[ii].rho = 0.0;
      prf[ii].enc = Mdisk;
    }/* for(int ii = tail; ii < NRADBIN; ii++){ */

#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++)
      prf[ii].psi = prf[ii].enc / prf[ii].rad;


    /** calculate spherical averaged profile of potential (external part) */
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

      prf[ii].psi += 4.0 * M_PI * (Pini[idx] + (Pext[idx] + psi) * logrbin_sph * M_LN10 / 3.0);

      Pext[0] += psi * (double)(1 << (1 + (idx    )));
      Pext[1] += psi * (double)(1 << (1 + (idx ^ 1)));
    }/* for(int ii = NRADBIN - 3; ii >= 0; ii--){ */


    /** multiply overall factor */
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++)
      prf[ii].psi *= CAST_R2D(newton);
  }/* for(int kk = 0; kk < ndisk; kk++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn diffAxisymmetricPotential
 *
 * @brief Calculate derivatives of potential field.
 *
 * @param (maxLev) maximum level of nested grid
 * @return (disk) physical properties of the disk component
 *
 * @sa findIdxSpherical
 */
void diffAxisymmetricPotential(const int maxLev, const disk_data disk)
{
  __NOTE__("%s\n", "start");


  profile *prf;
  prf = disk.prf;
  const double invlogrbin = disk.invlogrbin;

#ifndef USE_POTENTIAL_SCALING_SCHEME
  const int Nj = NDISKBIN_VER;
#else///USE_POTENTIAL_SCALING_SCHEME
  const int Nj = 1;
#endif//USE_POTENTIAL_SCALING_SCHEME


  for(int lev = 0; lev < maxLev; lev++){
    const double invRbin = 1.0 / ldexp(disk.hh, -lev);

    /** generate R-profile of \pdv{\Phi}{R} and \pdv[2]{\Phi}{R} */
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
#ifdef  USE_POTENTIAL_SCALING_SCHEME
      const int jj = 0;
#else///USE_POTENTIAL_SCALING_SCHEME
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
#endif//USE_POTENTIAL_SCALING_SCHEME
	{
	  const double Phi_m = disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, (ii > 0) ? (ii - 1) : (0), jj)];/**< symmetry @ R = 0 (ii == 0) */
	  const double Phi_0 = disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)];
	  double Phi_p = 0.0;

	  if( ii < (NDISKBIN_HOR - 1) )
	    Phi_p = disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii + 1, jj)];
	  else{
	    if( lev > 0 ){
	      const int j0 = jj >> 1;
	      const int im = ii >> 1;
	      const int ip = im + 1;

	      const int jm = (j0 > 0) ? (j0 - 1) : (0);/**< symmetry @ z = 0 */
	      const int jp = j0 + 1;

	      const double Phi_pm = 0.75 * disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, ip, jm)] + 0.25 * disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, im, jm)];
	      const double Phi_p0 = 0.75 * disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, ip, j0)] + 0.25 * disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, im, j0)];
	      const double Phi_pp = 0.75 * disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, ip, jp)] + 0.25 * disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, im, jp)];

	      Phi_p = 0.75 * Phi_p0 + 0.25 * ((jj & 1) ? (Phi_pp) : (Phi_pm));
	    }/* if( lev > 0 ){ */
	    else{
	      /* Phi_p = Phi_0;/\**< sufficiently far from the focusing domain *\/ */
	      Phi_p = 2.0 * Phi_0 - Phi_m;
	    }
	  }/* else{ */
	  disk. dPhidR [INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii, jj)] = 0.5     * invRbin * (Phi_p - Phi_m);
	  disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii, jj)] = invRbin * invRbin * (Phi_p + Phi_m - 2.0 * Phi_0);
	}
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */


    /** detect spikes of \pdv[2]{\Phi}{R} profile and smooth out by cubic spline interpolation */
    {
      const double hh = ldexp(disk.hh, -lev);
      static double spline_xx[NDISKBIN_HOR + 2], spline_ff[NDISKBIN_HOR + 2], spline_bp[NDISKBIN_HOR + 2], spline_f2[NDISKBIN_HOR + 2];
#ifdef  USE_POTENTIAL_SCALING_SCHEME
      const int jj = 0;
#else///USE_POTENTIAL_SCALING_SCHEME
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
#endif//USE_POTENTIAL_SCALING_SCHEME
	{
	  spline_xx[0] = 0.0;
	  spline_ff[0] = 1.25 * disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 0, jj)] - 0.25 * disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1, jj)];

	  spline_xx[1] = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)];
	  spline_ff[1] = disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 0, jj)];
	  int Nspline = 2;

	  double prev = disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, Nj, lev, 0, jj)];
	  double curr = disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, Nj, lev, 1, jj)];

	  for(int ii = 1; ii < NDISKBIN_HOR - 1; ii++){
	    double next = disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii + 1, jj)];

	    if( (prev * curr >= 0.0) || (curr * next >= 0.0) ){
	      spline_xx[Nspline] = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
	      spline_ff[Nspline] = disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)];
	      Nspline++;
	    }/* if( (prev * curr >= 0.0) || (curr * next >= 0.0) ){ */
	    else
	      disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] = DBL_MAX;

	    prev = curr;
	    curr = next;
	  }/* for(int ii = 1; ii < NDISKBIN_HOR - 1; ii++){ */

	  spline_xx[Nspline] = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, NDISKBIN_HOR - 1)];
	  spline_ff[Nspline] = disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, NDISKBIN_HOR - 1, jj)];
	  Nspline++;

	  spline_xx[Nspline] = hh + disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, NDISKBIN_HOR - 1)];

	  double S = 0.0, Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0, Syy = 0.0, Sxxx = 0.0, Sxxy = 0.0, Sxxxx = 0.0;
	  for(int ii = NDISKBIN_HOR - 3; ii < NDISKBIN_HOR; ii++){
	    const double xx = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
	    const double yy = disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)];
	    S += 1.0;
	    Sx += xx;
	    Sy += yy;
	    Sxx += xx * xx;
	    Sxy += xx * yy;
	    Syy += yy + yy;
	    Sxxx += xx * xx * xx;
	    Sxxy += xx * xx * yy;
	    Sxxxx += xx * xx * xx * xx;
	  }/* for(int ii = NDISKBIN_HOR - 3; ii < NDISKBIN_HOR; ii++){ */

	  const double aa = ((S * Sxxx - Sxx * Sx) * (S * Sxy - Sx * Sy) + (Sy * Sxx - S * Sxxy) * (S * Sxx - Sx * Sx)) / ((S * Sxxxx - Sxx * Sxx) * (S * Sxx - Sx * Sx) - (S * Sxxx - Sxx * Sx) * (S * Sxxx - Sxx * Sx));
	  const double bb = (S * Sxy - Sx * Sy - (S * Sxxx - Sxx * Sx) * aa) / (S * Sxx - Sx * Sx);
	  const double cc = (Sy - aa * Sxx - bb * Sx) / S;
	  spline_ff[Nspline] = cc + spline_xx[Nspline] * (bb + aa * spline_xx[Nspline]);
	  Nspline++;

	  if( Nspline < (NDISKBIN_HOR + 2) ){
	    for(int ii = Nspline; ii < NDISKBIN_HOR + 2; ii++){
	      spline_xx[Nspline] = spline_xx[Nspline - 1] + hh;
	      spline_ff[Nspline] = cc + spline_xx[Nspline] * (bb + aa * spline_xx[Nspline]);
	      Nspline++;
	    }/* for(int ii = Nspline; ii < NDISKBIN_HOR + 2; ii++){ */
	    double fpr = 2.0 * aa * (spline_ff[Nspline - 1] - bb);

	    genCubicSpline1D(Nspline, spline_xx, spline_ff, spline_bp, 0.0, fpr, spline_f2);

#pragma omp parallel for
	    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
	      if( disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii, jj)] > 0.5 * DBL_MAX ){
		const double RR = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
		disk. dPhidR [INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii, jj)] = getCubicSpline1stDifferential1D(RR, Nspline, spline_xx, spline_ff, spline_f2);
		disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii, jj)] = getCubicSpline2ndDifferential1D(RR, Nspline, spline_xx,            spline_f2);
	      }/* if( disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii, jj)] > 0.5 * DBL_MAX ){ */
	  }/* if( Nspline < (NDISKBIN_HOR + 2) ){ */
	}
    }

    /** superpose potential-field by spherical component(s) */
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      const double RR = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
      const double R2 = RR * RR;
#ifdef  USE_POTENTIAL_SCALING_SCHEME
      const int jj = 0;
#else///USE_POTENTIAL_SCALING_SCHEME
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
#endif//USE_POTENTIAL_SCALING_SCHEME
	{
	  /** r-derivatives of potential given by the spherical component(s) */
	  const double zz = disk.ver[INDEX2D(maxLev, NDISKBIN_VER, lev, jj)];
	  const double rad = sqrt(R2 + zz * zz);
	  const double rinv = 1.0 / rad;
	  const double drdR = RR * rinv;
	  const double d2rdR2 = (1.0 - drdR * drdR) * rinv;

	  /** find index corresponds to r = rad in prf */
	  double alpha;
	  const int idx = findIdxSpherical(rad, prf, invlogrbin, &alpha);
	  /** Psi = -Phi + Phi0 = -Phi (i.e. Phi0 is assumed to be zero) */
	  const double enc = (1.0 - alpha) * prf[idx].enc_tot + alpha * prf[1 + idx].enc_tot;
	  const double rho = (1.0 - alpha) * prf[idx].rho_tot + alpha * prf[1 + idx].rho_tot;
	  const double  dPhidr  = CAST_R2D(newton) *                           enc * rinv * rinv;
	  const double d2Phidr2 = CAST_R2D(newton) * (4.0 * M_PI * rho - 2.0 * enc * rinv * rinv * rinv);

	  /** R-derivatives of potential given by the all components */
	  disk. dPhidR [INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii, jj)] += dPhidr * drdR;
	  disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii, jj)] += dPhidr * drdR * d2rdR2 + d2Phidr2 * drdR * drdR;
	}
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  }/* for(int lev = 0; lev < maxLev; lev++){ */


#if 0
  for(int lev = 0; lev < maxLev; lev++){
    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
      fprintf(stderr, "%d\t%e\t%e\t%e\t%e\n", lev, disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], disk.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)], disk.dPhidR[INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii, 0)], disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, Nj, lev, ii, 0)]);
    fprintf(stderr, "\n");
  }
  exit(0);
#endif


  __NOTE__("%s\n", "end");
}


/**
 * @fn bisection4nestedGrid
 *
 * @brief Execute bisection for nested grid.
 *
 * @param (val) the target value
 * @param (num) number of data points
 * @param (tab) array contains data points
 * @param (logrbl) true when data points are sampled in logarithmic space
 * @param (invbin) inverse of interval between grid points
 * @return (ratio) parameter for linear interpolation
 * @param (maxLev) maximum level of nested grid
 * @return (lev) the corresponding level of nested grid
 * @param (tab_lev) array to determin lev
 * @return lower index of the corresponding data point
 *
 * @sa bisection
 */
static inline int bisection4nestedGrid
(const double val, const int num, double * restrict tab, const bool logtbl, const double invbin, double * restrict ratio, const int maxLev, int * restrict lev, double * restrict tab_lev)
{
  /** 1. find the nested level */
  *lev = 0;
  for(int ii = maxLev - 1; ii >= 0; ii--)
    if( val <= tab_lev[ii] ){
      *lev = ii;
      break;
    }/* if( val <= tab_lev[ii] ){ */

  /** 2. find the index */
  int idx = bisection(val, num, &(tab[INDEX2D(maxLev, num, *lev, 0)]), logtbl, invbin, ratio);

  return (idx);
}


/**
 * @fn Psi_spherical
 *
 * @brief Get relative potential of spherical symmetric components.
 *
 * @param (rad) radius
 * @param (sph) radial profile of the component
 * @param (invlogbin_sph) inverse of logarithmic interval between grid points
 * @return potential at r = rad
 *
 * @sa findIdxSpherical
 */
static inline double Psi_spherical(const double rad, profile *sph, const double invlogrbin_sph)
{
  double ratio;
  const int idx = findIdxSpherical(rad, sph, invlogrbin_sph, &ratio);

  return ((1.0 - ratio) * sph[idx].psi_tot + ratio * sph[1 + idx].psi_tot);
}


static inline double getPsi
(const double zz, const int maxLev, double * restrict node_ver, double * restrict tab_lev,
 const double RR, const double R2, const int lev, const int ii, double * restrict hor, double * restrict ver,
 double * restrict Phi, profile * restrict sph, const double invlogrbin_sph)
{
  int lev_z = 0;
  double azz;
  int jzz = bisection4nestedGrid(zz, NDISKBIN_VER + 1, node_ver, false, 1.0, &azz, maxLev, &lev_z, tab_lev);
  azz /= (node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, 1 + jzz)] - node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, jzz)]);

  int iiR = ii;
  double aaR = 0;
  if( lev > lev_z ){
    /** reset iiR and aaR in coarser grid */
    iiR = bisection(RR, NDISKBIN_HOR, &hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_z, 0)], false, 1.0, &aaR);
    aaR /= (hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_z, 1 + iiR)] - hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_z, iiR)]);
    /* aaR /= (enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev_z, 1 + iiR)] - enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev_z, iiR)]); */
  }/* if( lev > lev_z ){ */
  if( lev < lev_z ){
    /* reset jzz and azz in coarser grid */
    lev_z = lev;
    jzz = bisection(zz, NDISKBIN_VER + 1, &(node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, 0)]), false, 1.0, &azz);
    /* azz /= (node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, 1 + jzz)] - node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, jzz)]); */
  }/* if( lev < lev_z ){ */


  /** correction about ``jzz'' and ``azz'' for zone centeric coordinate from edge coordinate */
  if( zz >= ver[INDEX2D(maxLev, NDISKBIN_VER, lev_z, 0)] ){
    if( zz <= ver[INDEX2D(maxLev, NDISKBIN_VER, lev_z, NDISKBIN_VER - 1)] ){
      jzz = bisection(zz, NDISKBIN_VER, &ver[INDEX2D(maxLev, NDISKBIN_VER, lev_z, 0)], false, 1.0, &azz);
      azz /= (ver[INDEX2D(maxLev, NDISKBIN_VER, lev_z, 1 + jzz)] - ver[INDEX2D(maxLev, NDISKBIN_VER, lev_z, jzz)]);
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


  const double Phi_disk =
    ((1.0 - azz) * Phi[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_z,	    iiR, jzz)] + azz * Phi[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_z,     iiR, 1 + jzz)]) * (1.0 - aaR) +
    ((1.0 - azz) * Phi[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_z, 1 + iiR, jzz)] + azz * Phi[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_z, 1 + iiR, 1 + jzz)]) *        aaR;


  return (-Phi_disk + Psi_spherical(sqrt(R2 + zz * zz), sph, invlogrbin_sph));
}


/**
 * @fn calcVerticalVdisp
 *
 * @brief Calculate velocity dispersion in z-direction.
 *
 * @param (ndisk) number of disk components
 * @param (maxLev) maximum level of nested grid
 * @return (disk_info) physical properties of the disk component
 *
 * @sa findIdxSpherical
 * @sa bisection
 */
#if 1
void calcVerticalVdisp(const int ndisk, const int maxLev, disk_data *disk_info)
{
  __NOTE__("%s\n", "start");


  profile *sph;
  sph = disk_info[0].prf;
  const double invlogrbin = disk_info[0].invlogrbin;

  double *Phi;
  Phi = disk_info[0].pot;

  double *ver;
  ver = disk_info[0].ver;

  double *hor;
  hor = disk_info[0].hor;

  double *node_ver;
  node_ver = disk_info[0].node_ver;

  double *tab_lev;
  tab_lev = (double *)malloc(sizeof(double) * maxLev);
  if( tab_lev == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tab_lev\n");  }
  for(int ll = 0; ll < maxLev; ll++)
    tab_lev[ll] = node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, ll, NDISKBIN_VER)];


  for(int kk = 0; kk < ndisk; kk++){
    /** load disk data */
    const disk_data disk = disk_info[kk];


    for(int lev = 0; lev < maxLev; lev++){
      double *sig;      sig = &(disk.sigmaz[INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]);
      /* double *Phi;      Phi = &(disk.pot   [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, 0)]); */
      /* double *hor;      hor = &(disk.hor   [INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]); */

#ifndef ENABLE_VARIABLE_SCALE_HEIGHT
      const double zd = disk.cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT

      /** calculate vertical velocity dispersion */
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	/** get potential @ (R, z = 0) */
	const double RR = hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
	const double Psi_R_0 = -Phi[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] + Psi_spherical(RR, sph, invlogrbin);

	/** get potential @ (R, z = zd) */
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
	const double zd = disk.zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
	const double Psi_R_zd = getPsi(zd, maxLev, node_ver, tab_lev, RR, RR * RR, lev, ii, hor, ver, Phi, sph, invlogrbin);

	sig[ii] = sqrt(Psi_R_0 - Psi_R_zd);
      }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
    }/* for(int lev = 0; lev < maxLev; lev++){ */
  }/* for(int kk = 0; kk < ndisk; kk++){ */


  free(tab_lev);

  __NOTE__("%s\n", "end");
}
#else
void calcVerticalVdisp(const int ndisk, const int maxLev, disk_data *disk_info)
{
  __NOTE__("%s\n", "start");


  for(int kk = 0; kk < ndisk; kk++){
    /** load disk data */
    const disk_data disk = disk_info[kk];
    const double invzbin0 = 1.0 / disk.hh;
    profile *sph;    sph = disk.prf;
    const double invlogrbin = disk.invlogrbin;

    for(int lev = 0; lev < maxLev; lev++){
      double *sig;      sig = &(disk.sigmaz[INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]);
      double *Phi;      Phi = &(disk.pot   [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, 0)]);
      double *hor;      hor = &(disk.hor   [INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]);

      const double invzbin = 1.0 / ldexp(disk.hh, -lev);
#ifndef ENABLE_VARIABLE_SCALE_HEIGHT
      const double zd = disk.cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT

      /** calculate vertical velocity dispersion */
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	int i0;
	double a0, Psi_s, Psi_d;

	/** get potential @ (R, z = 0) */
	const double RR = hor[ii];
	i0 = findIdxSpherical(RR, sph, invlogrbin, &a0);
	Psi_s = (1.0 - a0) * sph[i0].psi_tot + a0 * sph[i0 + 1].psi_tot;
	Psi_d = -Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 0)];
	const double Psi_R_0 = Psi_d + Psi_s;

	/** get potential @ (R, z = zd) */
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
	const double zd = disk.zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
	const double rad = sqrt(RR * RR + zd * zd);
	i0 = findIdxSpherical(rad, sph, invlogrbin, &a0);
	Psi_s = (1.0 - a0) * sph[i0].psi_tot + a0 * sph[i0 + 1].psi_tot;
	/** get i0 in the appropriate level of grid points (must be coarser than the current level) */
	i0 = bisection(zd, NDISKBIN_VER, disk.ver, false, invzbin0, &a0);
	int ll = (i0 != 0) ? ((int)ilog2((uint)(NDISKBIN_VER / i0))) : (lev);
	if( ll > lev )	  ll = lev;
	const int i1 = ii >> (lev - ll);
	i0 = bisection(zd, NDISKBIN_VER, &(disk.ver[INDEX2D(maxLev, NDISKBIN_VER, ll, 0)]), false, invzbin, &a0);
	Psi_d = (1.0 - a0) * (-Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, i1, i0)]) + a0 * (-Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, i1, 1 + i0)]);
	const double Psi_R_zd = Psi_d + Psi_s;

	sig[ii] = sqrt(Psi_R_0 - Psi_R_zd);
      }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
    }/* for(int lev = 0; lev < maxLev; lev++){ */
  }/* for(int kk = 0; kk < ndisk; kk++){ */


  __NOTE__("%s\n", "end");
}
#endif


/**
 * @fn distributeDiskParticles
 *
 * @brief Distribute N-body particles in the disk component.
 *
 * @return (Nuse) total number of N-body particles distributed
 * @return (body) N-body particles
 * @param (mass) mass of N-body particles
 * @param (maxLev) maximum level of nested grid
 * @param (disk) physical properties of the disk component
 * @param (rand) state of pseudo random numbers
 *
 * @sa convertCylindrical2Spherical
 * @sa findIdxSpherical
 * @sa bisection
 * @sa bisection4nestedGrid
 */
void distributeDiskParticles(ulong *Nuse, iparticle body, const real mass, const int maxLev, const disk_data disk, rand_state *rand)
{
  __NOTE__("%s\n", "start");


  /** load disk data */
  profile *prf;  prf = disk.prf;
  const double invlogrbin = disk.invlogrbin;

  double *node_hor;  node_hor = disk.node_hor;
  double *node_ver;  node_ver = disk.node_ver;
  double *zone_hor;  zone_hor = disk.hor;
  double *zone_ver;  zone_ver = disk.ver;
  double *pot;  pot =  disk.pot;
  double *enc;  enc =  disk.enc;
  double *rho;  rho = *disk.rhoSum;
  double *sig;  sig =  disk.sigmaz;

  double * dPhidR ;  dPhidR   =  disk. dPhidR;
  double *d2PhidR2;  d2PhidR2 =  disk.d2PhidR2;

#ifdef  SWEEP_HIGH_ALTITUDE_COMPONENT
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  double *zd;  zd = disk.zd;
#else///ENABLE_VARIABLE_SCALE_HEIGHT
  const double zd = disk.cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
#endif//SWEEP_HIGH_ALTITUDE_COMPONENT

  const ulong num = disk.cfg->num;
  const double Mmax = enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, 0, NDISKBIN_HOR)];

#ifdef  PROGRESS_REPORT_ON
#ifdef  USE_SFMTJUMP
  const ulong nunit = (ulong)ceilf(0.1f * (float)num / (float)omp_get_num_threads());
#else///USE_SFMTJUMP
  const ulong nunit = (ulong)ceilf(0.1f * (float)num);
#endif//USE_SFMTJUMP
  ulong stage = 1;
  ulong Npart = 0;
#ifdef  USE_SFMTJUMP
#pragma omp single nowait
#endif//USE_SFMTJUMP
  {
    fprintf(stdout, "#\n#\n# start distributing disk particles (%zu bodies: [%zu:%zu])\n", num, *Nuse, (*Nuse) + num - 1);
    fflush(stdout);
  }
#endif//PROGRESS_REPORT_ON

#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
  const double mass_2 = 0.5 * CAST_R2D(mass);
  double Tdisk = 0.0;
  double Wdisk = 0.0;
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION

#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
  const double frac = disk.cfg->vdisp_frac;
  if( disk.cfg->vdispR0 < 0.0 )
    disk.cfg->vdispR0 = disk.cfg->vdispz0;
#else///ENFORCE_EPICYCLIC_APPROXIMATION
  if( disk.cfg->vdispR0 < 0.0 ){
    if( disk.cfg->toomre > 0.0 ){
      /* get physical quantities @ R = Rs to obtain Toomre's Q-value */
      const int lev = (int)floor(log2(disk.hh * ((double)NDISKBIN_HOR - 1.5) * disk.invRd));
      const int ii = (int)floor(ldexp(disk.cfg->rs / disk.hh, lev) - 0.5);
      const double aR = (disk.cfg->rs - disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]) / (disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + ii)] - disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);

      const double Sigma0 = (1.0 - aR) * disk.Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] + aR * disk.Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + ii)];

#ifndef USE_POTENTIAL_SCALING_SCHEME
      const double  dPhidR  = (1.0 - aR) * disk. dPhidR [INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] + aR * disk. dPhidR [INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + ii, 0)];
      const double d2PhidR2 = (1.0 - aR) * disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] + aR * disk.d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + ii, 0)];
#else///USE_POTENTIAL_SCALING_SCHEME
      const double  dPhidR  = (1.0 - aR) * disk. dPhidR [INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] + aR * disk. dPhidR [INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + ii)];
      const double d2PhidR2 = (1.0 - aR) * disk.d2PhidR2[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] + aR * disk.d2PhidR2[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + ii)];
#endif//USE_POTENTIAL_SCALING_SCHEME

      const double Omega0 = sqrt( dPhidR  * disk.invRd);
      const double kappa0 = sqrt(d2PhidR2 + 3.0 * Omega0 * Omega0);

#pragma omp single
      __NOTE__("Sigma0 = %e, kappa0 = %e, Q = %e, vdispR0 = %e\n", Sigma0, kappa0, disk.cfg->toomre, disk.cfg->vdispR0);
#pragma omp barrier
      disk.cfg->vdispR0 = disk.cfg->toomre * sqrt(M_E) * 3.36 * CAST_R2D(newton) * Sigma0 / (DBL_MIN + kappa0);
    }/* if( disk.cfg->toomre > 0.0 ){ */
    else
      disk.cfg->vdispR0 = disk.cfg->vdispz0;
  }/* if( disk.cfg->vdispR0 < 0.0 ){ */
  const double vdispR0_2 = disk.cfg->vdispR0 * disk.cfg->vdispR0;
  assert(vdispR0_2 > 0.0);
#endif//ENFORCE_EPICYCLIC_APPROXIMATION

#ifdef  USE_ZANG_HOHL_1978_EQ5
  const double Jstar_inv = disk.cfg->Jstar_inv;
#endif//USE_ZANG_HOHL_1978_EQ5


  double *tab_lev;
  tab_lev = (double *)malloc(sizeof(double) * maxLev);
  if( tab_lev == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tab_lev\n");  }
  for(int lev = 0; lev < maxLev; lev++)
    tab_lev[lev] = enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, NDISKBIN_HOR)];

#ifdef  SWEEP_HIGH_ALTITUDE_COMPONENT
  double *z_tab_lev;
  z_tab_lev = (double *)malloc(sizeof(double) * maxLev);
  if( z_tab_lev == NULL ){    __KILL__(stderr, "ERROR: failure to allocate z_tab_lev\n");  }
  for(int lev = 0; lev < maxLev; lev++)
    z_tab_lev[lev] = node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev, NDISKBIN_VER)];
#endif//SWEEP_HIGH_ALTITUDE_COMPONENT

  /** distribute N-body particles */
#ifdef  USE_SFMTJUMP
#pragma omp for
#endif//USE_SFMTJUMP
  for(ulong ii = *Nuse; ii < (*Nuse) + num; ii++){
    /** determine Rg, phi, and z */
    /** set Rg */
    const double Renc = Mmax * UNIRAND_DBL(rand);
    double aRg, Rg;
    int iRg, lev;

    /** use table search */
    iRg = bisection4nestedGrid(Renc, NDISKBIN_HOR + 1, enc, false, 1.0, &aRg, maxLev, &lev, tab_lev);
    aRg /= (enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, 1 + iRg)] - enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, iRg)]);
    Rg = (1.0 - aRg) * node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, iRg)] + aRg * node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, 1 + iRg)];
    /** correction about ``iRg'' and ``aRg'' for zone centeric coordinate from edge coordinate */
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

    /** calculate vertical velocity dispersion at R = Rg */
    const double sigmaz = (1.0 - aRg) * sig[INDEX2D(maxLev, NDISKBIN_HOR, lev, iRg)] + aRg * sig[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + iRg)];
    assert(sigmaz >= 0.0);


    double zenc_max = 1.0;
    int jzz = 0;
    double azz = 0.0;

#ifdef  SWEEP_HIGH_ALTITUDE_COMPONENT
    /** zz should be less than DISK_DIMMING_HEIGHT * zd */
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
    double zmax = DISK_DIMMING_HEIGHT * ((1.0 - aRg) * zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, iRg)] + aRg * zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + iRg)]);
#ifdef  ADDITIONAL_CONDITION_FOR_SCALE_HEIGHT
    zmax = fmin(zmax, DISK_DIMMING_SCALE * disk.cfg->rs);
#endif//ADDITIONAL_CONDITION_FOR_SCALE_HEIGHT
#else///ENABLE_VARIABLE_SCALE_HEIGHT
    double zmax = DISK_DIMMING_HEIGHT * zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    zmax = fmin(zmax, disk.cfg->rc);

    if( zmax < node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, 0, NDISKBIN_VER)] ){
      int lev_z = 0;
      int ii_R = iRg;
      double aa_R = aRg;
      jzz = bisection4nestedGrid(zmax, NDISKBIN_VER + 1, node_ver, false, 1.0, &azz, maxLev, &lev_z, z_tab_lev);
      azz /= (node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, 1 + jzz)] - node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, jzz)]);

      if( lev > lev_z ){
    	/** reset ii_R and aa_R in coarser grid */
    	ii_R = bisection(Renc, NDISKBIN_HOR + 1, &(enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev_z, 0)]), false, 1.0, &aa_R);
    	/* aa_R /= (enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev_z, 1 + ii_R)] - enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev_z, ii_R)]); */
      }/* if( lev > lev_z ){ */
      if( lev < lev_z ){
	/* reset jzz and azz in coarser grid */
    	lev_z = lev;
    	jzz = bisection(zmax, NDISKBIN_VER + 1, &(node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, 0)]), false, 1.0, &azz);
    	/* azz /= (node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, 1 + jzz)] - node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, jzz)]); */
      }/* if( lev < lev_z ){ */

      /* zenc_max = */
      /* 	(rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev_z,     ii_R,     jzz)] * (1.0 - azz) + */
      /* 	 rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev_z,     ii_R, 1 + jzz)] *        azz    ) * (1.0 - aa_R) + */
      /* 	(rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev_z, 1 + ii_R,     jzz)] * (1.0 - azz) + */
      /* 	 rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev_z, 1 + ii_R, 1 + jzz)] *        azz    ) *        aa_R; */
      /* zenc_max = */
      /* 	rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev_z, ii_R,     jzz)] * (1.0 - azz) + */
      /* 	rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev_z, ii_R, 1 + jzz)] *        azz; */
      /* zenc_max = rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev_z, ii_R, jzz)]; */
      zenc_max = fmin(rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev_z, ii_R, jzz)], rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev_z, 1 + ii_R, jzz)]);
    }/* if( zmax < node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, 0, NDISKBIN_VER)] ){ */

    zenc_max = fmin(zenc_max, 1.0);
#endif//SWEEP_HIGH_ALTITUDE_COMPONENT


    /** set z using table search */
    double zz = 0.0;
    jzz = 0;
    azz = 0.0;

#ifdef  SWEEP_HIGH_ALTITUDE_COMPONENT
    const int lev_stock = lev;
    const int iRg_stock = iRg;
    const double aRg_stock = aRg;
    zz = zmax * 2.0;
    while( zz > zmax ){
      /** recover the original version */
      lev = lev_stock;
      iRg = iRg_stock;
      aRg = aRg_stock;
#endif//SWEEP_HIGH_ALTITUDE_COMPONENT
      const double zenc = UNIRAND_DBL(rand) * zenc_max;

      while( true ){
	if( zenc > (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev, iRg, NDISKBIN_VER)] + DBL_EPSILON) ){
	  lev--;
	  if( lev < 0 ){	    lev = 0;	    break;	  }
	  /** reset iRg and aRg in coarser grid */
	  iRg = bisection(Renc, NDISKBIN_HOR + 1, &(enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, 0)]), false, 1.0, &aRg);
	  aRg /= (enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, 1 + iRg)] - enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, iRg)]);
	  /** correction about ``iRg'' and ``aRg'' for zone centeric coordinate from edge coordinate */
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
	}/* if( zenc > (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev, iRg, NDISKBIN_VER)] + DBL_EPSILON) ){ */
	else
	  break;
      }/* while( true ){ */

      /** use table search */
      jzz = bisection(zenc, NDISKBIN_VER + 1, &rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev, iRg, 0)], false, 1.0, &azz);
      azz /= (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev, iRg, 1 + jzz)] - rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, lev, iRg, jzz)]);
      zz = (1.0 - azz) * node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev, jzz)] + azz * node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev, 1 + jzz)];

      /** correction about ``jzz'' and ``azz'' for zone centeric coordinate from edge coordinate */
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

#ifdef  SWEEP_HIGH_ALTITUDE_COMPONENT
    }/* while( zz > zmax ){ */
#endif//SWEEP_HIGH_ALTITUDE_COMPONENT


    if( UNIRAND_DBL(rand) < 0.5 )
      zz *= -1.0;

    __NOTE__("%zu-th particle: guiding center determined\n", ii - (*Nuse));

    /** calculate potential at R = Rg, z = zz */
    const double diskpot =
      ((1.0 - azz) * pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev,     iRg, jzz)] + azz * pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev,     iRg, 1 + jzz)]) * (1.0 - aRg) +
      ((1.0 - azz) * pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + iRg, jzz)] + azz * pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 1 + iRg, 1 + jzz)]) *	       aRg;


    /** calculate escape velocity at r^2 = Rg^2 + zz^2 */
    double delta;
    const int idx = findIdxSpherical(sqrt(Rg * Rg + zz * zz), prf, invlogrbin, &delta);
    const double sphepot = -((1.0 - delta) * prf[idx].psi_tot + delta * prf[1 + idx].psi_tot);
    const double PhiRz = diskpot + sphepot;
    const double vesc2 = -2.0 * PhiRz;

#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
    Wdisk += mass_2 * (diskpot + sphepot);
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION

    /** set a scale factor to include effects of disk thickness */
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


    /** calculate physical quantities at R = Rg, z = 0 */
    /** determine sigmaR */
    /** calculate circular speed, epicyclic frequency, ... */
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

    double vcirc2 = Rg * Rg * Omega2;
    double vcirc  = sqrt(vcirc2);

    const double gam2inv = 0.25 * (3.0 + d2Phi / (DBL_MIN + Omega2));
    assert(gam2inv >= 0.0);
    const double gam2    = 1.0 / (DBL_MIN + gam2inv);
    assert(gam2 >= 0.0);
    const double sz2inv  = 1.0 / (DBL_MIN + sigmaz * sigmaz);
#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
    const double sigmap = DISK_PERP_VDISP(sigmaz, vcirc, frac);    assert(sigmap >= 0.0);
    const double sp2inv = 1.0 / (DBL_MIN + sigmap * sigmap);
    const double sR2inv = gam2inv * sp2inv;
    const double sigmaR = sqrt(gam2) * sigmap;    assert(sigmaR >= 0.0);
#else///ENFORCE_EPICYCLIC_APPROXIMATION
    const double sR2inv = 1.0 / DISK_RADIAL_VDISP2(vdispR0_2, Rg, disk.invRd);
    const double sp2inv = gam2 * sR2inv;
#ifdef  SPEEDUP_CONVERGENCE
    const double sigmaR = 1.0 / sqrt(DBL_MIN + sR2inv);    assert(sigmaR >= 0.0);
    const double sigmap = 1.0 / sqrt(DBL_MIN + sp2inv);    assert(sigmap >= 0.0);
    __NOTE__("i = %zu: sR2inv = %e, sp2inv = %e, vdispR0_2 = %e, Rg = %e, disk.invRd = %e, gam2 = %e\n", ii - (*Nuse), sR2inv, sp2inv, vdispR0_2, Rg, disk.invRd, gam2);
#endif//SPEEDUP_CONVERGENCE
#endif//ENFORCE_EPICYCLIC_APPROXIMATION


    /** determine particle position and velocity */
    double vx, vy, vz;
    double xx, yy;

    const double vmax2 = vesc2 - vcirc2;
    assert(vmax2 >= 0.0);
#ifndef  SPEEDUP_CONVERGENCE
    const double vmax  = sqrt(vmax2);
#endif//SPEEDUP_CONVERGENCE

    __NOTE__("%zu-th particle: examine velocity; %e, %e, %e, %e\n", ii - (*Nuse), sigmaR, sigmap, sigmaz, vmax2);
    double vR, vp;
    while( true ){
      while( true ){
#ifdef  SPEEDUP_CONVERGENCE
	while( true ){
	  vR = RANDVAL_DBL(rand);
	  vp = RANDVAL_DBL(rand);
	  vz = RANDVAL_DBL(rand);

	  if( vR * vR + vp * vp + vz * vz < 1.0 )
	    break;
	}/* while( true ){ */
	/** 8x sigma --> DF = e-32 ~ 1.3e-14 is sufficiently small probability */
	vR *= 8.0 * sigmaR;
	vp *= 8.0 * sigmap;
	vz *= 8.0 * sigmaz;
#else///SPEEDUP_CONVERGENCE
      	vR = vmax * RANDVAL_DBL(rand);
      	vp = vmax * RANDVAL_DBL(rand);
      	vz = vmax * RANDVAL_DBL(rand);
#endif//SPEEDUP_CONVERGENCE

	if( vR * vR + vp * vp + vz * vz < vmax2 )
	  break;
      }/* while( true ){ */

      /** use the Schwarzschild's DF */
      const double val = exp(-0.5 * (vR * vR * sR2inv + vp * vp * sp2inv + vz * vz * sz2inv));
      const double try = UNIRAND_DBL(rand);/**< maximum of DF is unity */
      if( val > try )	break;
    }/* while( true ){ */

#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
    Tdisk += mass_2 * vcirc2;
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION

#ifdef  USE_ZANG_HOHL_1978_EQ5
    /** flip circular velocity for retrograding particles */
    const double act = fmin(Rg * vcirc * Jstar_inv, 1.0);

    if( (disk.cfg->retrogradeFrac > 0.0) && (UNIRAND_DBL(rand) < (0.5 * (1.0 - 0.5 * act * (3.0 - act * act)))) )
      vcirc *= -1.0;
#endif//USE_ZANG_HOHL_1978_EQ5

    /** uniform distribution in polar angle */
    vp += vcirc;
    const double phi = 2.0 * M_PI * UNIRAND_DBL(rand);
    const double cosphi = cos(phi);
    const double sinphi = sin(phi);
    xx = Rg * cosphi;    vx = vR * cosphi - vp * sinphi;
    yy = Rg * sinphi;    vy = vR * sinphi + vp * cosphi;
    __NOTE__("%zu-th particle: velocity determined\n", ii - (*Nuse));

    body.pos[ii].x = CAST_D2R(xx);
    body.pos[ii].y = CAST_D2R(yy);
    body.pos[ii].z = CAST_D2R(zz);
    body.pos[ii].m = mass;

    body.acc[ii].x = body.acc[ii].y = body.acc[ii].z = body.acc[ii].pot = ZERO;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
    body.acc_ext[ii].x = body.acc_ext[ii].y = body.acc_ext[ii].z = body.acc_ext[ii].pot = ZERO;
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#ifdef  BLOCK_TIME_STEP
    body.vel[ii].x = CAST_D2R(vx);
    body.vel[ii].y = CAST_D2R(vy);
    body.vel[ii].z = CAST_D2R(vz);
#else///BLOCK_TIME_STEP
    body.vx[ii] = CAST_D2R(vx);
    body.vy[ii] = CAST_D2R(vy);
    body.vz[ii] = CAST_D2R(vz);
#endif//BLOCK_TIME_STEP

    body.idx[ii] = ii;

#ifdef  PROGRESS_REPORT_ON
    Npart++;
    if( Npart == (stage * nunit) ){
#ifdef  USE_SFMTJUMP
      fprintf(stdout, "# ~%zu%% on thread %d\n", stage * 10, omp_get_thread_num());
#else///USE_SFMTJUMP
      fprintf(stdout, "# ~%zu%%\n", stage * 10);
#endif//USE_SFMTJUMP
      fflush(stdout);
      stage++;
    }/* if( Npart == (stage * nunit) ){ */
#endif//PROGRESS_REPORT_ON
  }

  *Nuse += num;

#ifdef  PROGRESS_REPORT_ON
#ifdef  USE_SFMTJUMP
#pragma omp barrier
#pragma omp master
#endif//USE_SFMTJUMP
  fprintf(stdout, "# finish distributing disk particles (%zu bodies)\n", num);
#endif//PROGRESS_REPORT_ON

#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION

#ifdef  USE_SFMTJUMP
  static double Tdisk_tot, Wdisk_tot;
#pragma omp single
  Tdisk_tot = Wdisk_tot = 0.0;

#pragma omp atomic
  Tdisk_tot += Tdisk;
#pragma omp atomic
  Wdisk_tot += Wdisk;

#pragma omp barrier
#pragma omp single
  {
    disk.cfg->Tdisk = Tdisk_tot;
    disk.cfg->Wdisk = Wdisk_tot;
  }
#else///USE_SFMTJUMP
    disk.cfg->Tdisk = Tdisk;
    disk.cfg->Wdisk = Wdisk;
#endif//USE_SFMTJUMP
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION

  free(tab_lev);
#ifdef  SWEEP_HIGH_ALTITUDE_COMPONENT
  free(z_tab_lev);
#endif//SWEEP_HIGH_ALTITUDE_COMPONENT

  __NOTE__("%s\n", "end");
}


/**
 * @fn getEffectiveRadius
 *
 * @brief Evaluate effective radius of the disk component(s).
 *
 * @param (ndisk) number of disk components
 * @param (maxLev) maximum level of nested grid
 * @param (disk) physical properties of the disk component
 *
 * @sa bisection4nestedGrid
 */
void getEffectiveRadius(const int ndisk, const int maxLev, disk_data *disk)
{
  __NOTE__("%s\n", "start");


  double *tab_lev;
  tab_lev = (double *)malloc(sizeof(double) * maxLev);
  if( tab_lev == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tab_lev\n");  }

  for(int ii = 0; ii < ndisk; ii++){
    double *     enc;    enc      = disk[ii].enc;
    double *node_hor;    node_hor = disk[ii].node_hor;
    for(int ll = 0; ll < maxLev; ll++)
      tab_lev[ll] = enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, ll, NDISKBIN_HOR)];

    const double Mhalf = 0.5 * disk[ii].cfg->Mtot;

    double alp;
    int lev;
    int idx = bisection4nestedGrid(Mhalf, NDISKBIN_HOR + 1, enc, false, 1.0, &alp, maxLev, &lev, tab_lev);
    alp /= (enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, 1 + idx)] - enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, idx)]);

    disk[ii].cfg->Reff = (1.0 - alp) * node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, idx)] + alp * node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, 1 + idx)];
  }/* for(int ii = 0; ii < ndisk; ii++){ */

  free(tab_lev);


  __NOTE__("%s\n", "end");
}


/**
 * @fn findIdx4nestedGrid
 *
 * @brief Find a data element in the given nested array corresponding to the given radius.
 *
 * @param (RR) radius
 * @param (maxLev) maximum level of nested grid
 * @param (disk) physical properties of the disk component
 * @return (lev) the corresponding level of nested grid
 * @return (idx) the corresponding lower index
 * @return (alp) the corresponding coefficient for linear interpolation
 *
 * @sa bisection4nestedGrid
 */
void findIdx4nestedGrid(const double RR, const int maxLev, const disk_data disk, int * restrict lev, int * restrict idx, double * restrict alp)
{
  __NOTE__("%s\n", "start");


  double *tab_lev;
  tab_lev = (double *)malloc(sizeof(double) * maxLev);
  if( tab_lev == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tab_lev\n");  }

  double *node_hor;  node_hor = disk.node_hor;
  for(int ll = 0; ll < maxLev; ll++)
    tab_lev[ll] = node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, ll, NDISKBIN_HOR)];

  *idx = bisection4nestedGrid(RR, NDISKBIN_HOR + 1, node_hor, false, 1.0, alp, maxLev, lev, tab_lev);
  *alp /= (node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, *lev, 1 + (*idx))] - node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, *lev, *idx)]);

  free(tab_lev);


  __NOTE__("%s\n", "end");
}
