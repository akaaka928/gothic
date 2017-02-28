/**
 * @file blas.h
 *
 * @brief Header file for BLAS (Basic Linear Algebra Subprograms)
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/02/21 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef BLAS_H
#define BLAS_H


#ifndef BLAS_PROGRESS_REPORT_ON
#define BLAS_PROGRESS_REPORT_ON (10)
#endif//BLAS_PROGRESS_REPORT_ON


/**
 * @struct crs
 *
 * @brief structure for sparse matrix in CRS format
 */
typedef struct
{
  double *val;
  int    *col;
  int    *row;
} crs;

/**
 * @struct soaBiCGSTAB
 *
 * @brief structure for BiCGSTAB method in CRS format
 */
typedef struct
{
  crs mat;
  double *vec, *res, *sdw, *mid, *tmp, *Api, *Ati;
} soaBiCGSTAB;

/**
 * @struct soaPreConditioning
 *
 * @brief structure for preconditioner with ILU(0) in CRS format
 */
typedef struct
{
  crs ilu;
  double *Kri, *Kpi, *Kti;
} soaPreConditioning;


/* list of functions appeared in ``blas.c'' */
void getILU0(const int num, crs mat, crs ilu);
void bicgstab
(const crs mat, const int num, double * restrict vec, double * restrict sol,
 double * restrict res, double * restrict sdw, double * restrict mid, double * restrict tmp,
 double * restrict Api, double * restrict Ati,
 const double tol);
void pbicgstab
(const crs mat, const int num, double * restrict vec, double * restrict sol,
 double * restrict res, double * restrict sdw, double * restrict mid, double * restrict tmp,
 double * restrict Api, double * restrict Ati,
 const crs ilu, double * restrict Kri, double * restrict Kpi, double * restrict Kti,
 const double tol);


#endif//BLAS_H
