/**
 * @file external.c
 *
 * @brief Source code for setting an external fixed potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/18 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "macro.h"

#include "../misc/structure.h"

#include "profile.h"
#include "spline.h"

#include "external.h"


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


  static double rr[N_EXT_POT_SPHE], yy[N_EXT_POT_SPHE];
  static double bp[N_EXT_POT_SPHE], y2[N_EXT_POT_SPHE];


  for(int ii = 0; ii < kind; ii++){
    /** apply cubic spline interpolation */
    const int skip = NRADBIN / N_EXT_POT_SPHE;
#pragma omp parallel for
    for(int jj = 0; jj < N_EXT_POT_SPHE; jj++){
      const int kk = jj * skip;
      const double rad =  prf[ii][kk].rad;
      const double Phi = -prf[ii][kk].psi;

      rr[jj] = rad;    pot[ii].rad[jj]     = CAST_D2R(rad);
      yy[jj] = Phi;    pot[ii].Phi[jj].val = CAST_D2R(Phi);
    }/* for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){ */
    genCubicSpline1D(N_EXT_POT_SPHE, rr, yy, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, y2);

    for(int jj = 0; jj < N_EXT_POT_SPHE; jj++)
      pot[ii].Phi[jj].dr2 = CAST_D2R(y2[jj]);

    const real logrmin = LOG10(rr[                 0]);
    const real logrmax = LOG10(rr[N_EXT_POT_SPHE - 1]);
    pot[ii].logrmin = logrmin;
    pot[ii].logrbin = (logrmax - logrmin) / (real)(N_EXT_POT_SPHE - 1);
    pot[ii].num = N_EXT_POT_SPHE;
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
void superposePotFld1D(const int kind, const int skind, potential_field * restrict pot, potential_field * restrict sphe, potential_field * restrict disk)
{
  __NOTE__("%s\n", "start");

  /** assume all potential field share the same radius */
  sphe->logrmin = pot[0].logrmin;
  sphe->logrbin = pot[0].logrbin;
  sphe->num     = pot[0].num;
#pragma omp parallel for
  for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){
    sphe->rad[ii] = pot[0].rad[ii];

    real val = ZERO;
    real dr2 = ZERO;
    for(int kk = 0; kk < skind; kk++){
      val += pot[kk].Phi[ii].val;
      dr2 += pot[kk].Phi[ii].dr2;
    }/* for(int kk = 0; kk < skind; kk++){ */
    sphe->Phi[ii].val = val;
    sphe->Phi[ii].dr2 = dr2;
  }/* for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){ */

  if( kind > skind ){
    disk->logrmin = pot[0].logrmin;
    disk->logrbin = pot[0].logrbin;
    disk->num     = pot[0].num;
#pragma omp parallel for
    for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){
      disk->rad[ii] = pot[0].rad[ii];

      real val = ZERO;
      real dr2 = ZERO;
      for(int kk = skind; kk < kind; kk++){
	val += pot[kk].Phi[ii].val;
	dr2 += pot[kk].Phi[ii].dr2;
      }/* for(int kk = skind; kk < kind; kk++){ */
      disk->Phi[ii].val = val;
      disk->Phi[ii].dr2 = dr2;
    }/* for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){ */
  }/* if( kind > skind ){ */

  __NOTE__("%s\n", "end");
}
