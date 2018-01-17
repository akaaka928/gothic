/**
 * @file table.c
 *
 * @brief Source code for generating table of f' and f'' from f using cubic spline interpolation
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/17 (Wed)
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
#include <math.h>

#include "macro.h"
#include "name.h"

#include "spline.h"
#include "profile.h"
#include "table.h"


/**
 * @fn writeDensityProfileTableFormat
 *
 * @brief Print the expected format.
 *
 * @param (file) the specified file name
 */
static inline void writeDensityProfileTableFormat(char *filename)
{
  fprintf(stderr, "ERROR: data written in \"%s\" does not match with the expected format\n", filename);
  fprintf(stderr, "Expected format is below:\n");
  fprintf(stderr, "\tnum<int>: number of data points\n");
  fprintf(stderr, "\txx<double>\tff<double>: radius normalized by the scale radius (must be sorted in the ascending order) and density in arbitrary unit, respectively; num lines\n");
  __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
}

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
 * @fn getInterpolatedDensityProfile
 *
 * @brief Get interpolated density profile from the machine-readable table format.
 *
 * @param (num) number of data points
 * @return (prf) radial profile of the component
 * @return (xx) position of data points
 * @return (ff) value of data points
 *
 * @sa leastSquaresMethod
 * @sa genCubicSpline1D
 * @sa getCubicSpline1D
 * @sa getCubicSpline1stDifferential1D
 * @sa getCubicSpline2ndDifferential1D
 */
void getInterpolatedDensityProfile(const int num, profile * restrict prf, double * restrict xx, double * restrict ff)
{
  __NOTE__("%s\n", "start");


  /** extrapolate for the innermost position by least squares method */
  double pp, bb;
  const double rmin = 0.5 * (prf[          0].rad < xx[          NPUT]) ? (prf[          0].rad) : (xx[          NPUT]);
  const double rmax = 2.0 * (prf[NRADBIN - 1].rad > xx[num - 1 - NPUT]) ? (prf[NRADBIN - 1].rad) : (xx[num - 1 - NPUT]);
  leastSquaresMethod(NFIT, &xx[                 NPUT], &ff[                 NPUT], &pp, &bb);
  const double logrmin = log10(rmin);
  double logrbin = (log10(xx[NPUT]) - logrmin) / (double)NPUT;
  for(int ii = 0; ii < NPUT; ii++){
    xx[ii] = pow(10.0, logrmin + logrbin * (double)ii);
    ff[ii] = bb * pow(xx[ii], pp);
  }/* for(int ii = 0; ii < NPUT; ii++){ */
  const double fpl = bb * pp * pow(rmin, pp - 1.0);

  /** extrapolate for the outermost position by least squares method */
  leastSquaresMethod(NFIT, &xx[num - 1 - NFIT - NPUT], &ff[num - 1 - NFIT - NPUT], &pp, &bb);
#if 0
  fprintf(stdout, "pp = %e, bb = %e\n", pp, bb);
#endif
  const double logrmax = log10(rmax);
  logrbin = (log10(xx[num - 1 - NPUT]) - logrmax) / (double)NPUT;
  for(int ii = 0; ii < NPUT; ii++){
    xx[num - 1 - ii] = pow(10.0, logrmax + logrbin * (double)ii);
    ff[num - 1 - ii] = bb * pow(xx[num - 1 - ii], pp);
#if 0
    fprintf(stderr, "xx[%d] = %e, ff[%d] = %e\n", num - 1 - ii, xx[num - 1 - ii], num - 1 - ii, ff[num - 1 - ii]);
    fflush(NULL);
#endif
  }/* for(int ii = 0; ii < NPUT; ii++){ */
  const double fpr = bb * pp * pow(rmax, pp - 1.0);


  /** generate table for cubic spline interpolation */
  double *bp;  bp = (double *)malloc(num * sizeof(double));  if( bp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate bp\n");  }
  double *f2;  f2 = (double *)malloc(num * sizeof(double));  if( f2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate f2\n");  }
  genCubicSpline1D(num, xx, ff, bp, fpl, fpr, f2);


  /** asign interpolated density profile */
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    prf[ii].  rho     = getCubicSpline1D               (prf[ii].rad, num, xx, ff, f2);
    prf[ii]. drho_dr  = getCubicSpline1stDifferential1D(prf[ii].rad, num, xx, ff, f2);
    prf[ii].d2rho_dr2 = getCubicSpline2ndDifferential1D(prf[ii].rad, num, xx,     f2);
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */


  /** memory deallocation */
  free(bp);  free(f2);


  __NOTE__("%s\n", "end");
}

/**
 * @fn setDensityProfileTable
 *
 * @brief Read data table for the spherical component with memory allocation.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale length of the component
 * @param (file) the specified file name
 *
 * @sa writeDensityProfileTableFormat
 * @sa getInterpolatedDensityProfile
 */
void setDensityProfileTable(profile *prf, const double rs, char *file)
{
  __NOTE__("%s\n", "start");


  /** read data table */
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s", CFGFOLDER, file);
  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  bool success = true;

  int num;
  if( 1 != fscanf(fp, "%d", &num) )    success = false;
  num += 2 * NPUT;

  double *xx;  xx = (double *)malloc(num * sizeof(double));  if( xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xx\n");  }
  double *ff;  ff = (double *)malloc(num * sizeof(double));  if( ff == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ff\n");  }

  for(int ii = NPUT; ii < num - NPUT; ii++)
    success &= (2 == fscanf(fp, "%le\t%le", &xx[ii], &ff[ii]));

  fclose(fp);
  if( success != true )
    writeDensityProfileTableFormat(filename);


  /** scale the length to the computational unit */
  /** assume unity corresponds to the scale radius in the computational unit */
  for(int ii = NPUT; ii < num - NPUT; ii++)
    xx[ii] *= rs;


  /** asign interpolated density profile */
  getInterpolatedDensityProfile(num, prf, xx, ff);


  /** memory deallocation */
  free(xx);  free(ff);


  __NOTE__("%s\n", "end");
}
