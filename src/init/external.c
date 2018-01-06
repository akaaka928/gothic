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


#define NFIT_EXT (4)

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
    double pp, bb;

    /** extrapolate for the innermost position by least squares method */
    for(int jj = 0; jj < NFIT_EXT; jj++){
      xx[jj] =  prf[ii][jj].rad;
      ff[jj] = -prf[ii][jj].psi;
    }/* for(int jj = 0; jj < NFIT_EXT; jj++){ */
    leastSquaresMethod(NFIT_EXT, xx, ff, &pp, &bb);
    const double fpl = bb * pp * pow(xx[0], pp - 1.0);

    /** extrapolate for the outermost position by least squares method */
    for(int jj = 0; jj < NFIT_EXT; jj++){
      xx[jj] =  prf[ii][NRADBIN - NFIT_EXT + jj].rad;
      ff[jj] = -prf[ii][NRADBIN - NFIT_EXT + jj].psi;
    }/* for(int jj = 0; jj < NFIT_EXT; jj++){ */
    leastSquaresMethod(NFIT_EXT, xx, ff, &pp, &bb);
    const double fpr = bb * pp * pow(xx[NFIT_EXT - 1], pp - 1.0);

    /** apply cubic spline interpolation */
    const int skip = N_EXT_POT_SPHE / NRADBIN;
    for(int jj = 0; jj < N_EXT_POT_SPHE; jj++){
      const int kk = jj * skip;
      const double rad =  prf[ii][kk].rad;
      const double Phi = -prf[ii][kk].psi;

      rr[jj] = rad;    pot[ii].rad[jj]     = CAST_D2R(rad);
      yy[jj] = Phi;    pot[ii].Phi[jj].val = CAST_D2R(Phi);
    }/* for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){ */
    genCubicSpline1D(N_EXT_POT_SPHE, xx, yy, bp, ypl, ypr, y2);

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
