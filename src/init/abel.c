/**
 * @file abel.c
 *
 * @brief Source code for Abel transformation to deproject density profile
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/06/02 (Fri)
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

#include "magi.h"
#include "spline.h"
#include "profile.h"
#include "table.h"
#include "abel.h"


extern double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];


/**
 * @fn getColumnDensityDerivativeSersic
 *
 * @brief Calculate derivative of column density profile in Sersic law.
 *
 * @param (RR) position R
 * @param (cfg) physical properties of the component
 * @return 1st derivative of column density profile w.r.t. R
 */
double getColumnDensityDerivativeSersic(double RR, profile_abel_cfg cfg);
double getColumnDensityDerivativeSersic(double RR, profile_abel_cfg cfg)
{
  const double xx = cfg.bb * pow(RR * cfg.invRd, cfg.ninv);
  return (-cfg.ninv * xx * exp(-xx) / RR);
}

/**
 * @fn getColumnDensityDerivativeTwoPow
 *
 * @brief Calculate derivative of column density profile in double power-law.
 *
 * @param (RR) position R
 * @param (cfg) physical properties of the component
 * @return 1st derivative of column density profile w.r.t. R
 */
double getColumnDensityDerivativeTwoPow(double RR, profile_abel_cfg cfg);
double getColumnDensityDerivativeTwoPow(double RR, profile_abel_cfg cfg)
{
  const double xx = RR * cfg.invRd;
  const double xb = pow(xx, cfg.beta);
  return (-cfg.invRd * pow(xx, cfg.del) * pow(1.0 + xb, cfg.eps) * (cfg.alpha + cfg.gam * xb));
}

/**
 * @fn getColumnDensityDerivativeTable
 *
 * @brief Calculate derivative of column density profile in machine-readable table format.
 *
 * @param (RR) position R
 * @param (cfg) physical properties of the component
 * @return 1st derivative of column density profile w.r.t. R
 */
double getColumnDensityDerivativeTable(double RR, profile_abel_cfg cfg);
double getColumnDensityDerivativeTable(double RR, profile_abel_cfg cfg)
{
  return (getCubicSpline1stDifferential1D(RR, cfg.num, cfg.xx, cfg.yy, cfg.y2));
}


/**
 * @fn func4Abel
 *
 * @brief Function for Abel transformation.
 *
 * @param (RR) position R
 * @param (r2) r squared
 * @param (abel) physical properties of the component
 * @return 1st derivative of column density profile w.r.t. R over sqrt(R^2 - r^2)
 */
static inline double func4Abel(const double RR, const double r2, const abel_util abel)
{
  return (abel.getColumnDensityDerivative(RR, abel.cfg) / sqrt(RR * RR - r2));
}

/**
 * @fn gaussQD4Abel
 *
 * @brief Gaussian quadrature for Abel transformation.
 *
 * @param (min) minimum of R
 * @param (max) maximum of R
 * @param (r2) r squared
 * @param (abel) physical properties of the component
 * @return result of Gaussian quadrature
 *
 * @sa func4Abel
 */
double gaussQD4Abel(const double min, const double max, const double r2, const abel_util abel);
double gaussQD4Abel(const double min, const double max, const double r2, const abel_util abel)
{
  /* initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  if( NINTBIN & 1 )
    sum = gsl_gaussQD_weight[(NINTBIN >> 1)] * func4Abel(pls + mns * gsl_gaussQD_pos[(NINTBIN >> 1)], r2, abel);

  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (func4Abel(pls + mns * gsl_gaussQD_pos[ii], r2, abel) + func4Abel(pls - mns * gsl_gaussQD_pos[ii], r2, abel));

  /* finalization */
  return (-M_1_PI * mns * sum);
}


/**
 * @fn execAbelTransform
 *
 * @brief Execute Abel transformation.
 *
 * @return (prf) radial profile of the component
 * @param (cfg) physical properties of the component
 * @param (rmin) minimum of r
 * @param (rmax) maximum of r
 * @param (tmp) temporal data for machine-readable table format
 *
 * @sa gaussQD4Abel
 * @sa getInterpolatedDensityProfile
 */
void execAbelTransform(profile *prf, const profile_cfg cfg, const double rmin, const double rmax, const profile_abel_cfg tmp)
{
  __NOTE__("%s\n", "start");


  /* configure the specified profile */
  abel_util abel;
  abel.cfg.invRd = 1.0 / cfg.rs;
  switch( cfg.kind ){
  case SPHSERSIC:    /**< spherical symmetric profile which gives Sersic profile in the column density profile */
    abel.cfg.ninv = 1.0 / cfg.n_sersic;
    abel.cfg.bb   = cfg.b_sersic;
    abel.getColumnDensityDerivative = getColumnDensityDerivativeSersic;
    break;
  case SIGTWOPOW:    /**< spherical symmetric profile which gives Two-power profile in the column density profile */
    abel.cfg.alpha = cfg.twopower_alpha;
    abel.cfg.beta  = cfg.twopower_beta;
    abel.cfg.gam   = cfg.twopower_gamma;
    abel.cfg.del   = -1.0 - abel.cfg.alpha;
    abel.cfg.eps   = (abel.cfg.alpha - abel.cfg.beta - abel.cfg.gam) / abel.cfg.beta;
    abel.getColumnDensityDerivative = getColumnDensityDerivativeTwoPow;
    break;
  case TABLE_SIG:    /**< spherical symmetric profile which gives the specified column density profile in the table form */
    abel.cfg = tmp;
    abel.getColumnDensityDerivative = getColumnDensityDerivativeTable;
    break;
  default:
    __KILL__(stderr, "ERROR: inputted index %d is not implemented in this function.\n", cfg.kind);
    break;
  }/* switch( cfg.kind ){ */


  /* initialize temporary array and execute first touch */
  static double rad[NABEL], rho[NABEL];
  const double logrmin =  log10(rmin);
  const double logrbin = (log10(rmax) - logrmin) / (double)(NABEL - 1);
#pragma omp parallel for
  for(int ii = 0; ii < NABEL; ii++){
    rad[ii] = pow(10.0, logrmin + logrbin * (double)ii);
    rho[ii] = 0.0;
  }/* for(int ii = 0; ii < NABEL; ii++){ */


  /* execute Abel transform */
#pragma omp parallel for
  for(int ii = 0; ii < NABEL; ii++){
#ifdef  NDIVIDE_GAUSSQD4ABEL
    double sum = 0.0;
    double Rin = rad[ii];
    const double r2 = Rin * Rin;
    const double Rbin = (1.125 * rmax - Rin) / (double)NDIVIDE_GAUSSQD4ABEL;
    for(int iter = 0; iter < NDIVIDE_GAUSSQD4ABEL; iter++){
      sum += gaussQD4Abel(Rin, Rin + Rbin, r2, abel);
      Rin += Rbin;
    }/* for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){ */
    rho[ii] = sum;/* store the integral */
#else///NDIVIDE_GAUSSQD4ABEL
    rho[ii] = gaussQD4Abel(rad[ii], 1.125 * rmax, abel);
#endif//NDIVIDE_GAUSSQD4ABEL
  }/* for(int ii = 0; ii < NABEL; ii++){ */

#if 0
  for(int ii = 0; ii < NABEL; ii++)
    fprintf(stderr, "%e\t%e\n", rad[ii], rho[ii]);
  fflush(NULL);
#endif

  /* return deprojected density profile */
  getInterpolatedDensityProfile(NABEL, prf, rad, rho);

  __NOTE__("%s\n", "end");
}


/**
 * @fn writeColumnDensityProfileTableFormat
 *
 * @brief Print the expected format.
 *
 * @param (file) the specified file name
 */
static inline void writeColumnDensityProfileTableFormat(char *filename)
{
  fprintf(stderr, "ERROR: data written in \"%s\" does not match with the expected format\n", filename);
  fprintf(stderr, "Expected format is below:\n");
  fprintf(stderr, "\tnum<int>: number of data points\n");
  fprintf(stderr, "\txx<double>\tff<double>: radius normalized by the scale radius (must be sorted in the ascending order) and column density in arbitrary unit, respectively; num lines\n");
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
}


/**
 * @fn readColumnDensityProfileTable
 *
 * @brief Read data table for the spherical component with memory allocation.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale length of the component
 * @param (file) the specified file name
 * @param (cfg) physical properties of the component
 *
 * @sa writeColumnDensityProfileTableFormat
 * @sa leastSquaresMethod
 * @sa genCubicSpline1D
 * @sa execAbelTransform
 */
void readColumnDensityProfileTable(profile *prf, const double rs, char *file, const profile_cfg cfg)
{
  __NOTE__("%s\n", "start");


  /* read data table */
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s", CFGFOLDER, file);
  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  bool success = true;

  int num;
  if( 1 != fscanf(fp, "%d", &num) )    success = false;
  num += 2 * NPUT;

  double *xx;  xx = (double *)malloc(num * sizeof(double));  if( xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xx");  }
  double *ff;  ff = (double *)malloc(num * sizeof(double));  if( ff == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ff");  }

  for(int ii = NPUT; ii < num - NPUT; ii++)
    success &= (2 == fscanf(fp, "%le\t%le", &xx[ii], &ff[ii]));

  fclose(fp);
  if( success != true )
    writeColumnDensityProfileTableFormat(filename);


  /* scale the length to the computational unit */
  /* assume unity corresponds to the scale radius in the computational unit */
  for(int ii = NPUT; ii < num - NPUT; ii++)
    xx[ii] *= rs;


  /** extrapolate for the innermost position by least squares method */
  /** extrapolate for the outermost position by least squares method */
  double pp, bb;
  const double rmin = 0.5 * (prf[          0].rad < xx[          NPUT]) ? (prf[          0].rad) : (xx[          NPUT]);
  const double rmax = 2.0 * (prf[NRADBIN - 1].rad > xx[num - 1 - NPUT]) ? (prf[NRADBIN - 1].rad) : (xx[num - 1 - NPUT]);
  leastSquaresMethod(NFIT, &xx[NPUT], &ff[NPUT], &pp, &bb);
  const double logrmin = log10(rmin);
  double logrbin = (log10(xx[NPUT]) - logrmin) / (double)NPUT;
  for(int ii = 0; ii < NPUT; ii++){
    xx[ii] = pow(10.0, logrmin + logrbin * (double)ii);
    ff[ii] = bb * pow(xx[ii], pp);
  }/* for(int ii = 0; ii < NPUT; ii++){ */
  const double fpl = bb * pp * pow(rmin, pp - 1.0);
  leastSquaresMethod(NFIT, &xx[num - 1 - NFIT - NPUT], &ff[num - 1 - NFIT - NPUT], &pp, &bb);
  const double logrmax = log10(rmax);
  logrbin = (log10(xx[num - 1 - NPUT]) - logrmax) / (double)NPUT;
  for(int ii = 0; ii < NPUT; ii++){
    xx[num - 1 - ii] = pow(10.0, logrmax + logrbin * (double)ii);
    ff[num - 1 - ii] = bb * pow(xx[num - 1 - ii], pp);
  }/* for(int ii = 0; ii < NPUT; ii++){ */
  const double fpr = bb * pp * pow(rmax, pp - 1.0);


  /** execute cubic spline interpolation to derive 1st derivative of the column density profile */
  double *f2;  f2 = (double *)malloc(num * sizeof(double));  if( f2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate f2");  }
  double *bp;  bp = (double *)malloc(num * sizeof(double));  if( bp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate bp");  }
  genCubicSpline1D(num, xx, ff, bp, fpl, fpr, f2);


  /** execute Abel transformation */
  profile_abel_cfg abel;
  abel.num = num;
  abel.xx  = xx;
  abel.yy  = ff;
  abel.y2  = f2;
  execAbelTransform(prf, cfg, rmin, rmax, abel);


  /* memory deallocation */
  free(bp);  free(f2);
  free(xx);  free(ff);


  __NOTE__("%s\n", "end");
}


/**
 * @fn readColumnDensityTable4Disk
 *
 * @brief Read data table for the disk component with memory allocation.
 *
 * @return (prf) radial profile of the component
 * @param (rs) scale length of the component
 * @param (file) the specified file name
 * @return (num) number of data points
 * @return (xx) position of data points
 * @return (ff) column density at the position of data points
 * @return (f2) coefficients in cubic spline interpolation
 * @return (bp) coefficients in cubic spline interpolation
 *
 * @sa writeColumnDensityProfileTableFormat
 * @sa leastSquaresMethod
 * @sa genCubicSpline1D
 */
void readColumnDensityTable4Disk(profile *prf, const double rs, char *file, int *num, double **xx, double **ff, double **f2, double **bp)
{
  __NOTE__("%s\n", "start");


  /* read data table */
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s", CFGFOLDER, file);
  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  bool success = true;

  if( 1 != fscanf(fp, "%d", num) )    success = false;
  *num += 2 * NPUT;

  *xx = (double *)malloc((*num) * sizeof(double));  if( *xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate *xx");  }
  *ff = (double *)malloc((*num) * sizeof(double));  if( *ff == NULL ){    __KILL__(stderr, "ERROR: failure to allocate *ff");  }

  for(int ii = NPUT; ii < (*num) - NPUT; ii++)
    success &= (2 == fscanf(fp, "%le\t%le", &((*xx)[ii]), &((*ff)[ii])));

  fclose(fp);
  if( success != true )
    writeColumnDensityProfileTableFormat(filename);


  /* scale the length to the computational unit */
  /* assume unity corresponds to the scale radius in the computational unit */
  for(int ii = NPUT; ii < (*num) - NPUT; ii++)
    (*xx)[ii] *= rs;


  /** extrapolate for the innermost position by least squares method */
  /** extrapolate for the outermost position by least squares method */
  double pp, bb;
  const double rmin = 0.5 * (prf[          0].rad < (*xx)[          NPUT]) ? (prf[          0].rad) : ((*xx)[          NPUT]);
  const double rmax = 2.0 * (prf[NRADBIN - 1].rad > (*xx)[(*num) - 1 - NPUT]) ? (prf[NRADBIN - 1].rad) : ((*xx)[(*num) - 1 - NPUT]);
  leastSquaresMethod(NFIT, &((*xx)[NPUT]), &((*ff)[NPUT]), &pp, &bb);
  const double logrmin = log10(rmin);
  double logrbin = (log10((*xx)[NPUT]) - logrmin) / (double)NPUT;
  for(int ii = 0; ii < NPUT; ii++){
    (*xx)[ii] = pow(10.0, logrmin + logrbin * (double)ii);
    (*ff)[ii] = bb * pow((*xx)[ii], pp);
  }/* for(int ii = 0; ii < NPUT; ii++){ */
  const double fpl = bb * pp * pow(rmin, pp - 1.0);
  leastSquaresMethod(NFIT, &((*xx)[(*num) - 1 - NFIT - NPUT]), &((*ff)[(*num) - 1 - NFIT - NPUT]), &pp, &bb);
  const double logrmax = log10(rmax);
  logrbin = (log10((*xx)[(*num) - 1 - NPUT]) - logrmax) / (double)NPUT;
  for(int ii = 0; ii < NPUT; ii++){
    (*xx)[(*num) - 1 - ii] = pow(10.0, logrmax + logrbin * (double)ii);
    (*ff)[(*num) - 1 - ii] = bb * pow((*xx)[(*num) - 1 - ii], pp);
  }/* for(int ii = 0; ii < NPUT; ii++){ */
  const double fpr = bb * pp * pow(rmax, pp - 1.0);


  /** execute cubic spline interpolation to specify column density at arbitrary R */
  *f2 = (double *)malloc((*num) * sizeof(double));  if( *f2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate *f2");  }
  *bp = (double *)malloc((*num) * sizeof(double));  if( *bp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate *bp");  }
  genCubicSpline1D(*num, *xx, *ff, *bp, fpl, fpr, *f2);


  __NOTE__("%s\n", "end");
}
