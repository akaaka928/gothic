/*************************************************************************\
 *                                                                       *
                  last updated on 2016/08/11(Thu) 17:12:31
 *                                                                       *
 *    BLAS (Basic Linear Algebra Subprograms)                            *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
/* #define DEBUG_MODE_FOR_BLAS_C */
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
//-------------------------------------------------------------------------
#include <macro.h>
//-------------------------------------------------------------------------
#include "blas.h"
//-------------------------------------------------------------------------
#ifdef  DEBUG_MODE_FOR_BLAS_C
#include <math.h>
#endif//DEBUG_MODE_FOR_BLAS_C
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline double dnrm2(const int head, const int tail, double * restrict vec)
{
  //-----------------------------------------------------------------------
  double ret = 0.0;
  for(int ii = head; ii < tail; ii++)
    ret += vec[ii] * vec[ii];
  //-----------------------------------------------------------------------
  return (ret);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double ddot(const int head, const int tail, double * restrict vec0, double * restrict vec1)
{
  //-----------------------------------------------------------------------
  double ret = 0.0;
  for(int ii = head; ii < tail; ii++)
    ret += vec0[ii] * vec1[ii];
  //-----------------------------------------------------------------------
  return (ret);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void spmv(crs mat, const int head, const int tail, double * restrict vec, double * restrict ret);
void spmv(crs mat, const int head, const int tail, double * restrict vec, double * restrict ret)
{
  //-----------------------------------------------------------------------
  for(int ii = head; ii < tail; ii++){
    //---------------------------------------------------------------------
    ret[ii] = 0.0;
    //---------------------------------------------------------------------
    const int jhead = mat.row[ii    ];
    const int jtail = mat.row[ii + 1];
    //---------------------------------------------------------------------
    for(int jj = jhead; jj < jtail; jj++)
      ret[ii] += mat.val[jj] * vec[mat.col[jj]];
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void getILU0(const int num, crs mat, crs ilu)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* initialization */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num + 1; ii++)
    ilu.row[ii] = mat.row[ii];
  for(int ii = 0; ii < ilu.row[num]; ii++){
    ilu.val[ii] = mat.val[ii];
    ilu.col[ii] = mat.col[ii];
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* execute Incomplete LU decomposition in the inner-product form */
  /* diagonal elements in L are fixed to unity */
  /* assume all diagonal elements are non-zero */
  //-----------------------------------------------------------------------
  for(int ii = 1; ii < num; ii++){
    //---------------------------------------------------------------------
    const int khead = ilu.row[ii    ];
    const int ktail = ilu.row[ii + 1];
    int kdiag = ktail;
    for(int kk = khead; kk < ktail; kk++){
      //-------------------------------------------------------------------
      const int idx = mat.col[kk];
      //-------------------------------------------------------------------
      if( idx == ii ){	kdiag = kk    ;	break;      }
      if( idx  > ii ){	__KILL__(stderr, "ERROR: at least one diagonal element is zero\n");      }
      //-------------------------------------------------------------------
    }/* for(int kk = khead; kk < ktail; kk++){ */
    //---------------------------------------------------------------------
    for(int mm = khead; mm < kdiag; mm++){
      //-------------------------------------------------------------------
      const int kk = ilu.col[mm];
      //-------------------------------------------------------------------
      const int lhead = ilu.row[kk    ];
      const int ltail = ilu.row[kk + 1];
      for(int ll = lhead; ll < ltail; ll++){
	//-----------------------------------------------------------------
	const int idx = ilu.col[ll];
	//-----------------------------------------------------------------
	if( idx == kk ){	  ilu.val[mm] /= ilu.val[ll];	  break;	}
	if( idx  > kk ){	  __KILL__(stderr, "ERROR: at least one diagonal element is zero\n");	}
	//-----------------------------------------------------------------
      }/* for(int ll = lhead; ll < ltail; ll++){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      for(int jj = mm + 1; jj < ktail; jj++){
	//-----------------------------------------------------------------
	const int jcol = ilu.col[jj];
	double Akj = 0.0;
	for(int ll = lhead; ll < ltail; ll++){
	  //---------------------------------------------------------------
	  const int idx = ilu.col[ll];
	  //---------------------------------------------------------------
	  if( idx == jcol ){	    Akj = ilu.val[ll];	    break;	  }
	  if( idx  > jcol )	    break;
	  //---------------------------------------------------------------
	}/* for(int ll = ldiag + 1; ll < ltail; ll++){ */
	//-----------------------------------------------------------------
	ilu.val[jj] -= ilu.val[mm] * Akj;
	//-----------------------------------------------------------------
      }/* for(int jj = mm + 1; jj < ktail; jj++){ */
      //-------------------------------------------------------------------
    }/* for(int mm = khead; mm < ktail; mm++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 1; ii < num; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void mulLUinv(const int num, crs ilu, double * restrict ini, double * restrict tmp, double * restrict ret);
void mulLUinv(const int num, crs ilu, double * restrict ini, double * restrict tmp, double * restrict ret)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* forward elimination */
  //-----------------------------------------------------------------------
  tmp[0] = ini[0];
  ret[0] = tmp[0];
  for(int ii = 1; ii < num; ii++){
    //---------------------------------------------------------------------
    tmp[ii] = ini[ii];
    //---------------------------------------------------------------------
    const int jhead = ilu.row[ii    ];
    const int jtail = ilu.row[ii + 1];
    for(int jj = jhead; jj < jtail; jj++){
      //-------------------------------------------------------------------
      const int kk = ilu.col[jj];
      if( kk == ii ){	ret[ii] = tmp[ii];	break;      }
      if( kk  > ii ){	__KILL__(stderr, "ERROR: the diagonal element would be skipped\n");      }
      //-------------------------------------------------------------------
      tmp[ii] -= ilu.val[jj] * tmp[kk];
      //-------------------------------------------------------------------
    }/* for(int jj = jhead; jj < jtail; jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 1; ii < num; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* backward substitution */
  //-----------------------------------------------------------------------
  for(int ii = num - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const int jhead = ilu.row[ii    ];
    const int jtail = ilu.row[ii + 1];
    for(int jj = jtail - 1; jj >= jhead; jj--){
      //-------------------------------------------------------------------
      const int kk = ilu.col[jj];
      if( kk == ii ){	ret[ii] /= ilu.val[jj];	break;      }
      if( kk  < ii ){	__KILL__(stderr, "ERROR: the diagonal element would be skipped\n");      }
      //-------------------------------------------------------------------
      ret[ii] -= ilu.val[jj] * ret[kk];
      //-------------------------------------------------------------------
    }/* for(int jj = jtail - 1; jj >= jhead; jj--){ */
    //---------------------------------------------------------------------
  }/* for(int ii = num - 2; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* mat :: input          :: matrix A in CRS format */
/* num :: input          :: # of elements in vector */
/* vec :: input          :: vector b */
/* sol :: input / output :: vector x, initialized before this function */
/* res :: temporary      :: residual vector r */
/* sdw :: temporary      :: shadow residual vector r0^* */
/* mid :: temporary      :: intermediate vector p */
/* tmp :: temporary      :: temporary vector t */
/* Api :: temporary      :: product of matrix A and vector pi */
/* Ati :: temporary      :: product of matrix A and vector ti */
//-------------------------------------------------------------------------
/* BiCGSTAB */
void bicgstab
(const crs mat, const int num, double * restrict vec, double * restrict sol,
 double * restrict res, double * restrict sdw, double * restrict mid, double * restrict tmp,
 double * restrict Api, double * restrict Ati,
 const double tol)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# BiCGSTAB start: # of columns is %d, tolerance value is %e\n", num, tol);
  fflush(stdout);
  int steps = 0;
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set convergence criterion */
  //-----------------------------------------------------------------------
  const double tol2 = tol * tol * dnrm2(0, num, vec);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* initialization: calculate residual vector r0, r0^*, and p0 */
  //-----------------------------------------------------------------------
  spmv(mat, 0, num, sol, sdw);
  for(int ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    res[ii] = vec[ii] - sdw[ii];
    sdw[ii] = res[ii];
    mid[ii] = res[ii];
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* iterative procedure */
  //-----------------------------------------------------------------------
  double sdw_res = ddot(0, num, sdw, res);
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    /* derive alpha */
    //---------------------------------------------------------------------
    spmv(mat, 0, num, mid, Api);
    const double alpha = sdw_res / ddot(0, num, sdw, Api);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* derive omega */
    //---------------------------------------------------------------------
    for(int ii = 0; ii < num; ii++)
      tmp[ii] = res[ii] - alpha * Api[ii];
    spmv(mat, 0, num, tmp, Ati);
    const double omega = ddot(0, num, tmp, Ati) / dnrm2(0, num, Ati);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* update vectors */
    //---------------------------------------------------------------------
    for(int ii = 0; ii < num; ii++){
      //-------------------------------------------------------------------
      sol[ii] += alpha * mid[ii] + omega * tmp[ii];
      res[ii]  =         tmp[ii] - omega * Ati[ii];
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < num; ii++){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* convergence test */
    //---------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
    const double error = dnrm2(0, num, res);
    fprintf(stdout, "#\t%d-th iteration: error is %e\n", steps, error);
    fflush(stdout);
    if( error < tol2 )      break;
    steps++;
#else///PROGRESS_REPORT_ON
    if( dnrm2(num, res) < tol2 )
      break;
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* update intermediate vector */
    //---------------------------------------------------------------------
    double beta = alpha / (omega * sdw_res);
    sdw_res = ddot(0, num, sdw, res);
    beta *= sdw_res;
    for(int ii = 0; ii < num; ii++)
      mid[ii] = res[ii] + beta * (mid[ii] - omega * Api[ii]);
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* modified preconditioned BiCGSTAB */
/* mat :: input          :: matrix A in CRS format */
/* num :: input          :: # of elements in vector */
/* vec :: input          :: vector b */
/* sol :: input / output :: vector x, initialized before this function */
/* res :: temporary      :: residual vector ri */
/* sdw :: temporary      :: shadow residual vector r0^* */
/* mid :: temporary      :: intermediate vector pi */
/* tmp :: temporary      :: temporary vector ti */
/* Api :: temporary      :: product of matrix A and vector      pi */
/* Ati :: temporary      :: product of matrix A and vector K^-1 ti */
/* ilu :: input          :: Incomplete LU factorized matrix K in CRS format */
/* Kri :: temporary      :: product of matrix K^-1 and vector   ri */
/* Kpi :: temporary      :: product of matrix K^-1 and vector A pi */
/* Kti :: temporary      :: product of matrix K^-1 and vector   ti */
//-------------------------------------------------------------------------
void pbicgstab
(const crs mat, const int num, double * restrict vec, double * restrict sol,
 double * restrict res, double * restrict sdw, double * restrict mid, double * restrict tmp,
 double * restrict Api, double * restrict Ati,
 const crs ilu, double * restrict Kri, double * restrict Kpi, double * restrict Kti,
 const double tol)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# preconditioned BiCGSTAB start: # of columns is %d, tolerance value is %e\n", num, tol);
  fflush(stdout);
  int steps = 0;
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* start OpenMP parallelized region */
  //-----------------------------------------------------------------------
  /* shared variables */
  double dnrm_ful, ddot_ful, sdw_res;
  //-----------------------------------------------------------------------
#pragma omp parallel
  {
    //---------------------------------------------------------------------
    /* configuration of OpenMP */
    //---------------------------------------------------------------------
    const int rank = omp_get_thread_num();
    const int size = omp_get_num_threads();
    //---------------------------------------------------------------------
    const int head = (num * (    rank)) / size;
    const int tail = (num * (1 + rank)) / size;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* set convergence criterion */
    //---------------------------------------------------------------------
    double dnrm_loc = dnrm2(head, tail, vec);
#pragma omp single
    dnrm_ful = 0.0;
#pragma omp atomic
    dnrm_ful += dnrm_loc;
#pragma omp barrier
    const double tol2 = tol * tol * dnrm_ful;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* initialization: calculate residual vector r0, r0^*, and p0 */
    //---------------------------------------------------------------------
    spmv(mat, head, tail, sol, sdw);
    for(int ii = head; ii < tail; ii++)
      res[ii] = vec[ii] - sdw[ii];
#pragma omp barrier
#pragma omp single
    mulLUinv(num, ilu, res, sdw, Kri);
    for(int ii = head; ii < tail; ii++){
      sdw[ii] = Kri[ii];
      mid[ii] = Kri[ii];
    }/* for(int ii = head; ii < tail; ii++){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* iterative procedure */
    //---------------------------------------------------------------------
    double ddot_loc = ddot(head, tail, sdw, Kri);
#pragma omp single
    sdw_res = 0.0;
#pragma omp atomic
    sdw_res += ddot_loc;
#pragma omp barrier
    //---------------------------------------------------------------------
    while( true ){
      //-------------------------------------------------------------------
      /* derive alpha */
      //-------------------------------------------------------------------
      spmv(mat, head, tail, mid, Api);
#pragma omp barrier
#pragma omp single
      mulLUinv(num, ilu, Api, Kti, Kpi);
      ddot_loc = ddot(head, tail, sdw, Kpi);
#pragma omp single
      ddot_ful = 0.0;
#pragma omp atomic
      ddot_ful += ddot_loc;
#pragma omp barrier
      const double alpha = sdw_res / ddot_ful;
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* derive omega */
      //-------------------------------------------------------------------
      for(int ii = head; ii < tail; ii++){
	tmp[ii] = res[ii] - alpha * Api[ii];
	Kti[ii] = Kri[ii] - alpha * Kpi[ii];
      }/* for(int ii = head; ii < tail; ii++){ */
#pragma omp barrier
      spmv(mat, head, tail, Kti, Ati);
      ddot_loc = ddot(head, tail, tmp, Ati);
      dnrm_loc = dnrm2(head, tail, Ati);
#pragma omp single nowait
      ddot_ful = 0.0;
#pragma omp single
      dnrm_ful = 0.0;
#pragma omp atomic
      ddot_ful += ddot_loc;
#pragma omp atomic
      dnrm_ful += dnrm_loc;
#pragma omp barrier
      const double omega = ddot_ful / dnrm_ful;
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* update vectors */
      //-------------------------------------------------------------------
      for(int ii = head; ii < tail; ii++){
	//-----------------------------------------------------------------
	sol[ii] += alpha * mid[ii] + omega * Kti[ii];
	res[ii]  =         tmp[ii] - omega * Ati[ii];
	//-----------------------------------------------------------------
      }/* for(int ii = head; ii < tail; ii++){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* convergence test */
      //-------------------------------------------------------------------
      dnrm_loc = dnrm2(head, tail, res);
#pragma omp barrier
#pragma omp single
      dnrm_ful = 0.0;
#pragma omp atomic
      dnrm_ful += dnrm_loc;
#pragma omp barrier
#ifdef  PROGRESS_REPORT_ON
#pragma omp single
      {
	fprintf(stdout, "#\t%d-th iteration: error is %e, goal is %e\n", steps, dnrm_ful, tol2);
	fflush(stdout);
      }
      if( dnrm_ful < tol2 )      break;
#pragma omp single nowait
      steps++;
#else///PROGRESS_REPORT_ON
      if( dnrm_ful < tol2 )      break;
#endif//PROGRESS_REPORT_ON
      //-------------------------------------------------------------------
#ifdef  DEBUG_MODE_FOR_BLAS_C
#pragma omp barrier
#pragma omp master
      {
	int status = fpclassify(dnrm_ful);
	if( (status != FP_NORMAL) && (status != FP_ZERO) ){
	  static char msg[64];
	  switch( status ){
	  case FP_NAN      :	sprintf(msg, "Not a Number"                                    );	break;
	  case FP_INFINITE :	sprintf(msg, "either positive infinity or negative inifinity"  );	break;
	  case FP_SUBNORMAL:	sprintf(msg, "too small to be represented in normalized format");	break;
	  }/* switch( fpclassify(errMax) ){ */
	  for(int ii = 0; ii < num; ii++)
	    fprintf(stdout, "%d\t%e\t%e\t%e\n", ii, vec[ii], sol[ii], res[ii]);
	  __KILL__(stderr, "ERROR: dnrm_ful is \"%s\".\n", msg);
	}/* if( fpclassify(errMax) != FP_NORMAL ){ */
      }
#endif//DEBUG_MODE_FOR_BLAS_C
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* update intermediate vector */
      //-------------------------------------------------------------------
      double beta = alpha / (omega * sdw_res);
#pragma omp single
      mulLUinv(num, ilu, res, Api, Kri);
      ddot_loc = ddot(head, tail, sdw, Kri);
#pragma omp single
      sdw_res = 0.0;
#pragma omp atomic
      sdw_res += ddot_loc;
#pragma omp barrier
      beta *= sdw_res;
      for(int ii = head; ii < tail; ii++)
	mid[ii] = Kri[ii] + beta * (mid[ii] - omega * Kpi[ii]);
#pragma omp barrier
      //-------------------------------------------------------------------
    }/* while( true ){ */
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  /* end OpenMP parallelized region */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#if 0
//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //-----------------------------------------------------------------------
  /* memory allocation */
  //-----------------------------------------------------------------------
  const int num = 12;
  //-----------------------------------------------------------------------
  /* allocate vectors */
  double *vec;  vec = (double *)malloc(num * sizeof(double));  if( vec == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vec");  }
  double *sol;  sol = (double *)malloc(num * sizeof(double));  if( sol == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sol");  }
  double *res;  res = (double *)malloc(num * sizeof(double));  if( res == NULL ){    __KILL__(stderr, "ERROR: failure to allocate res");  }
  double *sdw;  sdw = (double *)malloc(num * sizeof(double));  if( sdw == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sdw");  }
  double *mid;  mid = (double *)malloc(num * sizeof(double));  if( mid == NULL ){    __KILL__(stderr, "ERROR: failure to allocate mid");  }
  double *tmp;  tmp = (double *)malloc(num * sizeof(double));  if( tmp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp");  }
  double *Api;  Api = (double *)malloc(num * sizeof(double));  if( Api == NULL ){    __KILL__(stderr, "ERROR: failure to allocate Api");  }
  double *Ati;  Ati = (double *)malloc(num * sizeof(double));  if( Ati == NULL ){    __KILL__(stderr, "ERROR: failure to allocate Ati");  }
  //-----------------------------------------------------------------------
  const int nnz = 5 * num * num;/* estimation of possible maximum */
  double *val;  val = (double *)malloc( nnz      * sizeof(double));  if( val == NULL ){    __KILL__(stderr, "ERROR: failure to allocate val");  }
  int    *col;  col = (   int *)malloc( nnz      * sizeof(   int));  if( col == NULL ){    __KILL__(stderr, "ERROR: failure to allocate col");  }
  int    *row;  row = (   int *)malloc((num + 1) * sizeof(   int));  if( row == NULL ){    __KILL__(stderr, "ERROR: failure to allocate row");  }
  crs mat;
  mat.val = val;
  mat.col = col;
  mat.row = row;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* intialize matrix and vector */
  //-----------------------------------------------------------------------
  int idx = 0;
  row[0] = 0;
  for(int ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    if(  ii >  2                                         ){      val[idx] = 1.0;      col[idx] = ii - 3;      idx++;      }
    if( (ii >  0) && (ii != 3) && (ii != 6) && (ii != 9) ){      val[idx] = 1.0;      col[idx] = ii - 1;      idx++;      }
    {                                                            val[idx] = 4.0;      col[idx] = ii    ;      idx++;      }
    if( (ii < 11) && (ii != 2) && (ii != 5) && (ii != 8) ){      val[idx] = 1.0;      col[idx] = ii + 1;      idx++;      }
    if(  ii <  9                                         ){      val[idx] = 1.0;      col[idx] = ii + 3;      idx++;      }
    row[ii + 1] = idx;
  }/* for(int ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++)
    sol[ii] = (double)(ii + 1);
  //-----------------------------------------------------------------------
  spmv(mat, num, sol, vec);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++)
    sol[ii] = 0.0;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  bicgstab(mat, num, vec, sol, res, sdw, mid, tmp, Api, Ati, 1.0e-10);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++)
    fprintf(stdout, "x[%2d] = %e\n", ii, sol[ii]);
  /* solution is BiCGSTAB is confirmed */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  double *ival;  ival = (double *)malloc( nnz      * sizeof(double));  if( ival == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ival");  }
  int    *icol;  icol = (   int *)malloc( nnz      * sizeof(   int));  if( icol == NULL ){    __KILL__(stderr, "ERROR: failure to allocate icol");  }
  int    *irow;  irow = (   int *)malloc((num + 1) * sizeof(   int));  if( irow == NULL ){    __KILL__(stderr, "ERROR: failure to allocate irow");  }
  crs ilu;
  ilu.val = ival;
  ilu.col = icol;
  ilu.row = irow;
  //-----------------------------------------------------------------------
  getILU0(num, mat, ilu);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++)
    sol[ii] = 0.0;
  //-----------------------------------------------------------------------
  double * xi;  xi  = (double *)malloc(num * sizeof(double));  if(  xi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  xi");  }
  double *Kpi;  Kpi = (double *)malloc(num * sizeof(double));  if( Kpi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate Kpi");  }
  double *Kti;  Kti = (double *)malloc(num * sizeof(double));  if( Kti == NULL ){    __KILL__(stderr, "ERROR: failure to allocate Kti");  }
  //-----------------------------------------------------------------------
  pbicgstab(mat, num, vec, sol, res, sdw, mid, tmp, Api, Ati, ilu, xi, Kpi, Kti, 1.0e-10);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++)
    fprintf(stdout, "x[%2d] = %e\n", ii, sol[ii]);
  /* solution is preconditioned BiCGSTAB with updated is confirmed */
  //-----------------------------------------------------------------------
  free(xi);
  free(Kpi);
  free(Kti);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(vec);
  free(sol);
  free(res);
  free(sdw);
  free(mid);
  free(tmp);
  free(Api);
  free(Ati);
  //-----------------------------------------------------------------------
  free(val);
  free(col);
  free(row);
  //-----------------------------------------------------------------------
  free(ival);
  free(icol);
  free(irow);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif
//-------------------------------------------------------------------------
