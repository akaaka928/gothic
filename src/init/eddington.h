/**
 * @file eddington.h
 *
 * @brief Header file for Eddington's formula
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/06/19 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef EDDINGTON_H
#define EDDINGTON_H


#include "macro.h"

#include "../init/profile.h"


#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA
#define ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA


#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
/* #define NENEBIN (262144) */
/* #define NENEBIN (16384) */
#define NENEBIN (8192)
#else///ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
#define NENEBIN (1048576)
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF


/**
 * @struct dist_func
 *
 * @brief structure for distribution function
 */
typedef struct
{
  real ene, val;
} dist_func;


/**
 * @fn getDF
 *
 * @brief Get distribution function of the component based on linear interpolation.
 *
 * @param (ene) a specified energy
 * @param (df) distribution function
 * @param (Emin) minimum of energy bin
 * @param (invEbin) inverse of energy bin width
 * @return DF at the specified energy
 */
static inline double getDF(const double ene, dist_func *df, const double Emin, const double invEbin)
{
  const int ll = (int)((ene - Emin) * invEbin);
  return (df[ll].val + (df[ll + 1].val - df[ll].val) * (ene - df[ll].ene) / (df[ll + 1].ene - df[ll].ene));
}


/* list of functions appeared in ``eddington.c'' */
void integrateEddingtonFormula(const int skind, profile **prf,
#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
			       profile_cfg *cfg,
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
			       dist_func **fene);

#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
void calcVelocityDispersionProfile(const int skind, profile **prf,
#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
				   profile_cfg *cfg,
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
				   dist_func **df);
#endif//MAKE_VELOCITY_DISPERSION_PROFILE


#endif//EDDINGTON_H
