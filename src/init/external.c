/**
 * @file external.c
 *
 * @brief Source code for setting an external fixed potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/04 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "macro.h"

#include "../misc/structure.h"

#include "profile.h"
#include "spline.h"

#include "external.h"


#define NFIT_EXT (6)
const double tiny = 1.0e-20;

/**
 * @fn leastSquaresMethod
 *
 * @brief Data fitting by the least squares method (assume ax + b + c/x and sigma_i = 1).
 *
 * @param (num) number of data points
 * @param (xx) x
 * @param (yy) y = y(x)
 * @return (aa) the resultant a
 * @return (bb) the resultant b
 * @return (cc) the resultant c
 */
static inline void leastSquaresMethod(const int num, double * restrict xx, double * restrict yy, double * restrict aa, double * restrict bb, double * restrict cc)
{
  double Sxy = 0.0;
  double Sxx = 0.0;
  double Sx = 0.0;
  double S = 0.0;
  double Sy = 0.0;
  double Sp = 0.0;/**< S_{x'} */
  double Spy = 0.0;/**< S_{x'y} */
  double Spp = 0.0;/**< S_{x'x'} */

  for(int ii = 0; ii < num; ii++){
    const double x = xx[ii];
    const double y = yy[ii];
    const double xinv = 1.0 / (tiny + x);

    /** sigma_i^2 is assumed to be unity */
    S   += 1.0;
    Sx  += x;
    Sy  += y;
    Sxx += x * x;
    Sxy += x * y;
    Sp  += xinv;
    Spp += xinv * xinv;
    Spy += xinv * y;
  }/* for(int ii = 0; ii < num; ii++){ */

  const double a = ((Spp * Sx - S * Sp) * (Spp * Sy - Sp * Spy) + (Spy * S - Spp * Sxy) * (Spp * S - Sp * Sp)) / ((Spp * Sxx - S * S) * (Spp * S - Sp * Sp) - (Spp * Sx - S * Sp) * (Spp * Sx - S * Sp));
  const double b = (Spp * Sy - Sp * Spy - (Spp * Sx - S * Sp) * a) / (Spp * S - Sp * Sp);
  const double c = (Spy - a * S - b * Sp) / Spp;

  *aa = a;
  *bb = b;
  *cc = c;
}


/**
 * @fn genExtPotTbl1D
 *
 * @brief Generate table(s) for cubic spline interpolation to represent external fixed potential field(s).
 *
 * @param (kind) number of components
 * @param (prf) radial profile of component(s)
 * @return (pot) potential field for cubic spline interpolation
 */
void genExtPotTbl1D(const int kind, profile **prf, potential_field *pot)
{
  __NOTE__("%s\n", "start");


  static double xx[NFIT_EXT], ff[NFIT_EXT];
  static double rr[N_EXT_POT_SPHE], yy[N_EXT_POT_SPHE];
  static double bp[N_EXT_POT_SPHE], y2[N_EXT_POT_SPHE];


  potential table should be start from r = 0;


  for(int ii = 0; ii < kind; ii++){
    double aa, bb, cc;

    /** extrapolate for the innermost position by least squares method */
    for(int jj = 0; jj < NFIT_EXT; jj++){
      xx[jj] =  prf[ii][jj].rad;
      ff[jj] = -prf[ii][jj].psi;
    }/* for(int jj = 0; jj < NFIT_EXT; jj++){ */
    leastSquaresMethod(NFIT_EXT, xx, ff, &aa, &bb, &cc);
    rr[0] = 0.0;    pot[ii].rad[0] = ZERO;
    yy[0] = aa * rr[0] + bb + cc / (tiny + rr[0]);
    pot[ii].Phi[0].val = CAST_D2R(yy[0]);

    /** extrapolate for the outermost position by least squares method */
    for(int jj = 0; jj < NFIT_EXT; jj++){
      xx[jj] =  prf[ii][NRADBIN - NFIT_EXT + jj].rad;
      ff[jj] = -prf[ii][NRADBIN - NFIT_EXT + jj].psi;
    }/* for(int jj = 0; jj < NFIT_EXT; jj++){ */
    leastSquaresMethod(NFIT_EXT, xx, ff, &aa, &bb, &cc);
    rr[N_EXT_POT_SPHE - 1] = 10.0 * prf[ii][NRADBIN - 1].rad;    pot[ii].rad[N_EXT_POT_SPHE - 1] = CAST_D2R(rr[N_EXT_POT_SPHE - 1]);
    yy[N_EXT_POT_SPHE - 1] = aa * rr[N_EXT_POT_SPHE - 1] + bb + cc / (tiny + rr[N_EXT_POT_SPHE - 1]);
    pot[ii].Phi[N_EXT_POT_SPHE - 1].val = CAST_D2R(yy[N_EXT_POT_SPHE - 1]);

    /** apply cubic spline interpolation */
    const int skip = N_EXT_SPH / NRADBIN;
    for(int jj = 0; jj < N_EXT_SPH; jj++){
      const int kk = jj * skip;
      const double rad =  prf[ii][kk].rad;
      const double Phi = -prf[ii][kk].psi;

      const int ll = (N_EXT_CAP >> 1) + jj;
      rr[ll] = rad;    pot[ii].rad[ll]     = CAST_D2R(rad);
      yy[ll] = Phi;    pot[ii].Phi[ll].val = CAST_D2R(Phi);
    }/* for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){ */
    genCubicSpline1D(N_EXT_POT_SPHE, rr, yy, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, y2);

    for(int jj = 0; jj < N_EXT_POT_SPHE; jj++)
      pot[ii].Phi[jj].dr2 = CAST_D2R(y2[jj]);
  }/* for(int ii = 0; ii < kind; ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn superposePotFld1D
 *
 * @brief Superpose external fixed potential fields of multiple components.
 *
 * @param (kind) number of components
 * @param (skind) number of spherical symmetric components
 * @param (pot) potential field for cubic spline interpolation
 * @return (sphe) superposed potential field for cubic spline interpolation (only for spherical symmetric components)
 * @return (disk) superposed potential field for cubic spline interpolation (only for spherical averaged disk components)
 */
void superposePotFld1D(const int kind, const int skind, potential_field *pot, potential_field sphe, potential_field disk)
{
  __NOTE__("%s\n", "start");


  potential_field zero = {ZERO, ZERO};

  /** assume all potential field share the same radius */
#pragma omp parallel for
  for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){
    sphe.rad[ii] = pot[0].rad[ii];

    real val = ZERO;
    real dr2 = ZERO;
    for(int kk = 0; kk < skind; kk++){
      val += pot[kk].Phi[ii].val;
      dr2 += pot[kk].Phi[ii].dr2;
    }/* for(int kk = 0; kk < skind; kk++){ */
    sphe.Phi[ii].val = val;
    sphe.Phi[ii].dr2 = dr2;
  }/* for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){ */

  if( kind > skind ){
#pragma omp parallel for
    for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){
      disk.rad[ii] = pot[0].rad[ii];

      real val = ZERO;
      real dr2 = ZERO;
      for(int kk = skind; kk < kind; kk++){
	val += pot[kk].Phi[ii].val;
	dr2 += pot[kk].Phi[ii].dr2;
      }/* for(int kk = skind; kk < kind; kk++){ */
      disk.Phi[ii].val = val;
      disk.Phi[ii].dr2 = dr2;
    }/* for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){ */
  }/* if( kind > skind ){ */


  __NOTE__("%s\n", "end");
}
