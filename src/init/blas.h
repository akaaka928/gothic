/*************************************************************************\
 *                                                                       *
                  last updated on 2016/01/02(Sat) 16:25:54
 *                                                                       *
 *    Header File for BLAS (Basic Linear Algebra Subprograms)            *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef BLAS_H
#define BLAS_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef PROGRESS_REPORT_ON
#define PROGRESS_REPORT_ON
#endif//PROGRESS_REPORT_ON
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* structure for sparse matrix (in CRS format) */
//-------------------------------------------------------------------------
typedef struct
{
  double *val;
  int    *col;
  int    *row;
} crs;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "blas.c"
//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//BLAS_H
//-------------------------------------------------------------------------
