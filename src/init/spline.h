/**
 * @file spline.h
 *
 * @brief Header file for cubic spline interpolation
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
#ifndef SPLINE_H
#define SPLINE_H


#include <float.h>


#define NATURAL_CUBIC_SPLINE DBL_MAX


/* list of functions appeared in ``spline.c'' */
void genCubicSpline1D(const int num, double * restrict xx, double * restrict yy, double * restrict bp, const double ypl, const double ypr, double * restrict y2);
double getCubicSpline1D(const double pp, const int num, double * restrict xx, double * restrict yy, double * restrict y2);

double getCubicSpline1stDifferential1D(const double pp, const int num, double * restrict xx, double * restrict yy, double * restrict y2);
double getCubicSpline2ndDifferential1D(const double pp, const int num, double * restrict xx, double * restrict y2);

double getCubicSplineIntegral1D       (const int ii, double * restrict xx, double * restrict yy, double * restrict y2);
double getCubicSplineIntegral1D1stHalf(const int ii, double * restrict xx, double * restrict yy, double * restrict y2);
double getCubicSplineIntegral1D2ndHalf(const int ii, double * restrict xx, double * restrict yy, double * restrict y2);

void genCubicSpline2D1st(const int nx, const int ny, double * restrict yy, double * restrict ff, double * restrict f2, double * restrict bp_ful);
void genCubicSpline2D2nd(const double py, const int nx, double * restrict xx, const int ny, double * restrict yy, double * restrict ff, double * restrict f2, double * restrict ffx, double * restrict f2x, double * restrict bp);


#endif//SPLINE_H
