/**
 * @file abel.h
 *
 * @brief Header file for Abel transformation to deproject density profile
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/06/22 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef ABEL_H
#define ABEL_H


#include "../init/profile.h"


#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA
#define ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_ABEL
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA


#define NABEL (4096)
#ifndef ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_ABEL
#define NDIVIDE_GAUSSQD4ABEL (32)
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_ABEL


/**
 * @struct profile_abel_cfg
 *
 * @brief structure for Abel transformation
 */
typedef struct
{
  double *xx, *yy, *y2;/**< pointers for column density profile in table form */
  double invRd;
  double ninv, bb;/**< parameters for Sersic profile */
  double alpha, beta, gam, del, eps;/**< parameters for two-power model (del := -1 - alpha; eps := (alpha - beta - gamma) / beta) */
  int num;/**< parameter for column density profile in table form */
} profile_abel_cfg;
/**
 * @struct abel_util
 *
 * @brief structure for Abel transformation
 */
typedef struct
{
  double (*getColumnDensityDerivative)(double, profile_abel_cfg);
  profile_abel_cfg cfg;
} abel_util;


/* list of functions appeared in ``able.c'' */
void execAbelTransform(profile *prf, const profile_cfg cfg, const double rmin, const double rmax, const profile_abel_cfg tmp);

void readColumnDensityProfileTable(profile *prf, const double rs, char *file, const profile_cfg cfg);
void readColumnDensityTable4Disk  (profile *prf, const double rs, char *file, int *num, double **xx, double **ff, double **f2, double **bp);


#endif//ABEL_H
