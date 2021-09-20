/**
 * @file spline.c
 *
 * @brief Source code for cubic spline interpolation
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/10/26 (Thu)
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
#include <omp.h>

#include "macro.h"

#include "spline.h"


/**
 * @fn bisec
 *
 * @brief Execute bisection.
 *
 * @param (val) the target value
 * @param (num) number of data points
 * @param (tab) array contains data points
 * @return (ratio) parameter for linear interpolation
 * @return lower index of the corresponding data point
 */
static inline int bisec(const double val, const int num, double * restrict tab, double * restrict ratio)
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
      *ratio = (val - tab[ll]) / (tab[rr] - tab[ll]);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
  }/* while( true ){ */
}


/**
 * @fn genCubicSpline1D
 *
 * @brief Execute cubic spline interpolation in 1D.
 *
 * @param (num) number of data points
 * @param (xx) position of data points
 * @param (yy) value of data points
 * @return (bp) coefficients in cubic spline interpolation
 * @param (ypl) boundary condition at the lower limit
 * @param (ypr) boundary condition at the upper limit
 * @return (y2) coefficients in cubic spline interpolation
 */
void genCubicSpline1D(const int num, double * restrict xx, double * restrict yy, double * restrict bp, const double ypl, const double ypr, double * restrict y2)
{
  __NOTE__("%s\n", "start");

  /** set boundary condition at the left side */
  if( ypl >= NATURAL_CUBIC_SPLINE )
    y2[0] = bp[0] = 0.0;
  else{
    y2[0] = -0.5;
    bp[0] = (3.0 / (xx[1] - xx[0])) * ((yy[1] - yy[0]) / (xx[1] - xx[0]) - ypl);
  }/* else{ */


  /** forward sweep step in the tridiagonal matrix algorithm to get y2 */
  /** set -u' on y2[num], b' on bp[] */
  for(int ii = 1; ii < num - 1; ii++){
    const double xinv = 1.0 / (xx[ii + 1] - xx[ii - 1]);
    const double sig = (xx[ii] - xx[ii - 1]) * xinv;
    const double dinv = 1.0 / (2.0 + sig * y2[ii - 1]);
    y2[ii] = (sig - 1.0) * dinv;
    bp[ii] = (6.0 * ((yy[ii - 1] * (xx[ii + 1] - xx[ii]) + yy[ii + 1] * (xx[ii] - xx[ii - 1])) * xinv - yy[ii]) / ((xx[ii + 1] - xx[ii]) * (xx[ii] - xx[ii - 1])) - sig * bp[ii - 1]) * dinv;
  }/* for(int ii = 1; ii < num - 1; ii++){ */


  /** set boundary condition at the right side */
  /** default is natural cubic spline */
  double qr = 0.0;
  double ur = 0.0;
  if( ypr < NATURAL_CUBIC_SPLINE ){
    qr = 0.5;
    ur = (3.0 / (xx[num - 1] - xx[num - 2])) * (ypr - (yy[num - 1] - yy[num - 2]) / (xx[num - 1] - xx[num - 2]));
  }/* if( ypr < NATURAL_CUBIC_SPLINE ){ */
  y2[num - 1] = (ur - qr * bp[num - 2]) / (qr * y2[num - 2] + 1.0);


  /** backward substitution step in the tridiagonal matrix algorithm to get y2 */
  for(int ii = num - 2; ii >= 0; ii--)
    y2[ii] = bp[ii] + y2[ii] * y2[ii + 1];


  __NOTE__("%s\n", "end");
}

/**
 * @fn getCubicSpline1D
 *
 * @brief Get interpolated value based on cubic spline interpolation in 1D.
 *
 * @param (pp) position of the specified point
 * @param (num) number of data points
 * @param (xx) position of data points
 * @param (yy) value of data points
 * @param (y2) coefficients in cubic spline interpolation
 * @return interpolated value at pp
 *
 * @sa bisec
 */
double getCubicSpline1D(const double pp, const int num, double * restrict xx, double * restrict yy, double * restrict y2)
{
  double aa;
  const int ii = bisec(pp, num, xx, &aa);
  const double dx = xx[ii + 1] - xx[ii];

  return ((1.0 - aa) * yy[ii] + aa * (yy[ii + 1] + (aa - 1.0) * dx * dx * ((aa + 1.0) * y2[ii + 1] - (aa - 2.0) * y2[ii]) / 6.0));
}

/**
 * @fn getCubicSpline1stDifferential1D
 *
 * @brief Get 1st differential value at a specified point based on cubic spline interpolation in 1D.
 *
 * @param (pp) position of the specified point
 * @param (num) number of data points
 * @param (xx) position of data points
 * @param (yy) value of data points
 * @param (y2) coefficients in cubic spline interpolation
 * @return 1st differential value at pp
 *
 * @sa bisec
 */
double getCubicSpline1stDifferential1D(const double pp, const int num, double * restrict xx, double * restrict yy, double * restrict y2)
{
  double aa;
  const int ii = bisec(pp, num, xx, &aa);
  const double dx = xx[ii + 1] - xx[ii];

  return ((yy[ii + 1] - yy[ii]) / dx + dx * ((3.0 * aa * aa - 1.0) * y2[ii + 1] - (2.0 + 3.0 * aa * (aa - 2.0)) * y2[ii]) / 6.0);
}

/**
 * @fn getCubicSpline2ndDifferential1D
 *
 * @brief Get 2nd differential value at a specified point based on cubic spline interpolation in 1D.
 *
 * @param (pp) position of the specified point
 * @param (num) number of data points
 * @param (xx) position of data points
 * @param (y2) coefficients in cubic spline interpolation
 * @return 2nd differential value at pp
 *
 * @sa bisec
 */
double getCubicSpline2ndDifferential1D(const double pp, const int num, double * restrict xx, double * restrict y2)
{
  double aa;
  const int ii = bisec(pp, num, xx, &aa);

  return ((1.0 - aa) * y2[ii] + aa * y2[ii + 1]);
}

/**
 * @fn getCubicSplineIntegral1D
 *
 * @brief Get integrated value of a grid ii based on cubic spline interpolation in 1D.
 *
 * @param (ii) index of the specified grid
 * @param (xx) position of data points
 * @param (yy) value of data points
 * @param (y2) coefficients in cubic spline interpolation
 * @return integrated value of the grid ii
 */
double getCubicSplineIntegral1D(const int ii, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double dx = xx[ii + 1] - xx[ii];
  return (0.5 * dx * (yy[ii] + yy[ii + 1] - dx * dx * (y2[ii] + y2[ii + 1]) / 12.0));
}

/**
 * @fn getCubicSplineIntegral1D1stHalf
 *
 * @brief Get integrated value of a grid ii (first half only) based on cubic spline interpolation in 1D.
 *
 * @param (ii) index of the specified grid
 * @param (xx) position of data points
 * @param (yy) value of data points
 * @param (y2) coefficients in cubic spline interpolation
 * @return integrated value of the grid ii (first half only)
 */
double getCubicSplineIntegral1D1stHalf(const int ii, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double dx = xx[ii + 1] - xx[ii];
  return (0.0625 * dx * (7.0 * yy[ii] + yy[ii + 1] - dx * dx * (9.0 * y2[ii] + 7.0 * y2[ii + 1]) / 24.0));
}

/**
 * @fn getCubicSplineIntegral1D2ndHalf
 *
 * @brief Get integrated value of a grid ii (second half only) based on cubic spline interpolation in 1D.
 *
 * @param (ii) index of the specified grid
 * @param (xx) position of data points
 * @param (yy) value of data points
 * @param (y2) coefficients in cubic spline interpolation
 * @return integrated value of the grid ii (second half only)
 */
double getCubicSplineIntegral1D2ndHalf(const int ii, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double dx = xx[ii + 1] - xx[ii];
  return (0.0625 * dx * (yy[ii] + 7.0 * yy[ii + 1] - dx * dx * (7.0 * y2[ii] + 9.0 * y2[ii + 1]) / 24.0));
}


/**
 * @fn genCubicSpline2D1st
 *
 * @brief Execute cubic spline interpolation in 2D (first half).
 *
 * @param (nx) number of data points in x-direction
 * @param (ny) number of data points in y-direction
 * @param (yy) position of data points
 * @param (ff) value of data points
 * @return (f2) coefficients in cubic spline interpolation
 * @return (bp_ful) coefficients in cubic spline interpolation
 *
 * @sa genCubicSpline1D
 */
void genCubicSpline2D1st(const int nx, const int ny, double * restrict yy, double * restrict ff, double * restrict f2, double * restrict bp_ful)
{
  __NOTE__("%s\n", "start");

  /** assume natural cubic spline */
#pragma omp parallel for num_threads(CPUS_PER_PROCESS)
  for(int ii = 0; ii < nx; ii++){
    const int tidx = omp_get_thread_num();
    double *bp = &bp_ful[ny * tidx];
    genCubicSpline1D(ny, yy, &ff[INDEX2D(nx, ny, ii, 0)], bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, &f2[INDEX2D(nx, ny, ii, 0)]);
  }/* for(int ii = 0; ii < nx; ii++){ */

  __NOTE__("%s\n", "end");
}

/**
 * @fn genCubicSpline2D2nd
 *
 * @brief Execute cubic spline interpolation in 2D (second half).
 *
 * @param (py) y-position of the specified point
 * @param (nx) number of data points in x-direction
 * @param (xx) position of data points
 * @param (ny) number of data points in y-direction
 * @param (yy) position of data points
 * @param (ff) value of data points
 * @return (f2) coefficients in cubic spline interpolation
 * @return (ffx) coefficients in cubic spline interpolation
 * @return (f2x) coefficients in cubic spline interpolation
 * @return (bp) coefficients in cubic spline interpolation
 *
 * @sa bisec
 * @sa genCubicSpline1D
 */
void genCubicSpline2D2nd(const double py, const int nx, double * restrict xx, const int ny, double * restrict yy, double * restrict ff, double * restrict f2, double * restrict ffx, double * restrict f2x, double * restrict bp)
{
  __NOTE__("%s\n", "start");

  double aa;
  const int jj = bisec(py, ny, yy, &aa);
  const double dy = yy[jj + 1] - yy[jj];

#pragma omp parallel for
  for(int ii = 0; ii < nx; ii++)
    ffx[ii] = ((1.0 - aa) * ff[INDEX2D(nx, ny, ii, jj)] + aa * (ff[INDEX2D(nx, ny, ii, jj + 1)] + (aa - 1.0) * dy * dy * ((aa + 1.0) * f2[INDEX2D(nx, ny, ii, jj + 1)] - (aa - 2.0) * f2[INDEX2D(nx, ny, ii, jj)]) / 6.0));

  genCubicSpline1D(nx, xx, ffx, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, f2x);

  __NOTE__("%s\n", "end");
}


#if 0
#define NDATA (9)
#define NPLOT (256)
int main(void){
  static double xd[NDATA], yd[NDATA];
  xd[0] = 0.0;  yd[0] = 0.0;
  xd[1] = 0.1;  yd[1] = 0.0;
  xd[2] = 0.2;  yd[2] = 0.0;
  xd[3] = 0.4;  yd[3] = 0.2;
  xd[4] = 0.5;  yd[4] = 1.0;
  xd[5] = 0.6;  yd[5] = 0.2;
  xd[6] = 0.8;  yd[6] = 0.0;
  xd[7] = 0.9;  yd[7] = 0.0;
  xd[8] = 1.0;  yd[8] = 0.0;

  static double bp[NDATA], y2[NDATA];
#if 1
  genCubicSpline1D(NDATA, xd, yd, bp, 0.0, 0.0, y2);
#else
  genCubicSpline1D(NDATA, xd, yd, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, y2);
#endif

  const double xbin = 1.0 / (double)(NPLOT - 1);
  for(int ii = 0; ii < NPLOT; ii++){
    const double xp = xbin * (double)ii;
    printf("%e\t%e\t%e\t%e\n", xp,
	   getCubicSpline1D               (xp, NDATA, xd, yd, y2),
	   getCubicSpline1stDifferential1D(xp, NDATA, xd, yd, y2),
	   getCubicSpline2ndDifferential1D(xp, NDATA, xd,     y2));
  }

  return (0);
}
#endif
