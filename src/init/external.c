/**
 * @file external.c
 *
 * @brief Source code for setting an external fixed potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/24 (Wed)
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


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
/**
 * @fn genAdaptiveExtPotTbl1D
 *
 * @brief Generate table for cubic spline interpolation to represent external fixed potential field.
 *
 * @param (prf) radial profile of the component
 * @return (pot) potential field for cubic spline interpolation
 */
static inline void genAdaptiveExtPotTbl1D(profile *prf, potential_field *pot)
{
  __NOTE__("%s\n", "start");

  static double rr[N_EXT_POT_SPHE], yy[N_EXT_POT_SPHE];
  static double bp[N_EXT_POT_SPHE], y2[N_EXT_POT_SPHE];

  /** set sampling scale */
  const double dy = (-prf[NRADBIN - 1].psi + prf[0].psi) / (double)(N_EXT_POT_SPHE - 1);

  /** set innermost boundary */
  int num = 0;
  rr[num] =  prf[0].rad;
  double yold = -prf[0].psi;
  yy[num] = yold;
  num++;

  /** add sampling points */
  for(int kk = 1; kk < NRADBIN; kk++){
    double ynew = -prf[kk].psi;
    if( (ynew - yold) >= dy ){
      rr[num] = prf[kk].rad;
      yy[num] = ynew;
      yold = ynew;
      num++;
    }/* if( (ynew - yold) >= dy ){ */

    if( num >= N_EXT_POT_SPHE ){
      __KILL__(stderr, "ERROR: num (= %d) must be smaller than N_EXT_POT_SPHE (= %d)\n", num, N_EXT_POT_SPHE);
    }/* if( num >= N_EXT_POT_SPHE ){ */
  }/* for(int kk = 1; kk < NRADBIN; kk++){ */

  /** set outermost boundary */
  rr[num] =  prf[NRADBIN - 1].rad;
  yy[num] = -prf[NRADBIN - 1].psi;
  num++;

  /** apply cubic spline interpolation */
  pot->num = num;
  genCubicSpline1D(num, rr, yy, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, y2);

  for(int jj = 0; jj < num; jj++){
    pot->rad[jj]     = CAST_D2R(rr[jj]);
    pot->Phi[jj].val = CAST_D2R(yy[jj]);
    pot->Phi[jj].dr2 = CAST_D2R(y2[jj]);
  }/* for(int jj = 0; jj < num; jj++){ */

  __NOTE__("%s\n", "end");
}
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD


/**
 * @fn genExtPotTbl1D
 *
 * @brief Generate table(s) for cubic spline interpolation to represent external fixed potential field(s).
 *
 * @param (kind) number of components
 * @param (prf) radial profile of component(s)
 * @return (pot) potential field for cubic spline interpolation
 */
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
void genExtPotTbl1D(const int kind, profile **prf, potential_field *pot)
{
  __NOTE__("%s\n", "start");

  for(int ii = 0; ii < kind; ii++)
    genAdaptiveExtPotTbl1D(prf[ii], &pot[ii]);

  __NOTE__("%s\n", "end");
}
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
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

    /* const real logrmin = LOG10(rr[                 0]); */
    /* const real logrmax = LOG10(rr[N_EXT_POT_SPHE - 1]); */
    /* pot[ii].logrmin = logrmin; */
    /* pot[ii].logrbin = (logrmax - logrmin) / (real)(N_EXT_POT_SPHE - 1); */
    const double logrmin = log10(prf[ii][0].rad);
    const double logrmax = log10(prf[ii][(N_EXT_POT_SPHE - 1) * skip].rad);
    pot[ii].logrmin = CAST_D2R(logrmin);
    pot[ii].logrbin = CAST_D2R((logrmax - logrmin) / (double)(N_EXT_POT_SPHE - 1));
    pot[ii].num = N_EXT_POT_SPHE;
  }/* for(int ii = 0; ii < kind; ii++){ */


  __NOTE__("%s\n", "end");
}
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD


#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

/**
 * @fn genSuperposedPotFld1D
 *
 * @brief Generate superposed external, fixed potential-fields of multiple components.
 *
 * @param (kind) number of components
 * @param (skind) number of spherical symmetric components
 * @param (prf) radial profile of component(s)
 * @return (sphe) superposed potential field for cubic spline interpolation (only for spherical symmetric components)
 * @return (disk) superposed potential field for cubic spline interpolation (only for spherical averaged disk components)
 */
void genSuperposedPotFld1D(const int kind, const int skind, profile **prf, potential_field * restrict sphe, potential_field * restrict disk)
{
  __NOTE__("%s\n", "start");

  profile *tmp;
  tmp = (profile  *)malloc(sizeof(profile) * NRADBIN);
  if( tmp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp\n");  }

  /** generate superposed potential table for spherical component(s) */
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    tmp[ii].rad = prf[0][ii].rad;

    double psi = 0.0;
    for(int kk = 0; kk < skind; kk++)
      psi += prf[kk][ii].psi;

    tmp[ii].psi = psi;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */
  genAdaptiveExtPotTbl1D(tmp, sphe);

  /** generate superposed potential table for spherical averaged disk component(s) */
  if( kind > skind ){
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      double psi = 0.0;
      for(int kk = skind; kk < kind; kk++)
	psi += prf[kk][ii].psi;

      tmp[ii].psi = psi;
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */
    genAdaptiveExtPotTbl1D(tmp, disk);
  }/* if( kind > skind ){ */

  free(tmp);

  __NOTE__("%s\n", "end");
}

#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

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
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

#endif//SET_EXTERNAL_POTENTIAL_FIELD
