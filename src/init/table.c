/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:26:52
 *                                                                       *
 *    Generate table of f'& f''from f using  Cubic Spline Interpolation  *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//-------------------------------------------------------------------------
#include "macro.h"
//-------------------------------------------------------------------------
#include "spline.h"
#include "profile.h"
#include "table.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* read data table with memory allocation */
//-------------------------------------------------------------------------
static inline void writeDensityProfileTableFormat(char *filename)
{
  //-----------------------------------------------------------------------
  fprintf(stderr, "ERROR: data written in \"%s\" does not match with the expected format\n", filename);
  fprintf(stderr, "Expected format is below:\n");
  fprintf(stderr, "\tnum<int>: number of data points\n");
  fprintf(stderr, "\txx<double>\tff<double>: radius normalized by the scale radius (must be sorted in the ascending order) and density in arbitrary unit, respectively; num lines\n");
  __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void leastSquaredMethod(const int num, double * restrict xx, double * restrict yy, double * restrict pp, double * restrict bb)
{
  //-----------------------------------------------------------------------
  double SS, Sx, Sy, Sxx, Sxy;
  SS = Sx = Sy = Sxx = Sxy = 0.0;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    const double logx = log10(xx[ii]);
#if 0
    const double logy = log10(yy[ii] + DBL_MIN);
#else
    const double logy = log10(yy[ii]);
#endif
    SS  += 1.0;
    Sx  += logx;
    Sxx += logx * logx;
    Sy  +=        logy;
    Sxy += logx * logy;
    //---------------------------------------------------------------------
#if 0
    fprintf(stderr, "xx[%d] = %e, yy[%d] = %e\n", ii, xx[ii], ii, yy[ii]);
#endif
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
  *pp = (SS * Sxy - Sx * Sy) / (SS * Sxx - Sx * Sx);
  *bb = pow(10.0, (Sy - (*pp) * Sx) / SS);
  //-----------------------------------------------------------------------
#if 0
  fprintf(stderr, "SS = %e, Sx = %e, Sxx = %e, Sy = %e, Sxy = %e\n", SS, Sx, Sxx, Sy, Sxy);
#endif
  //-----------------------------------------------------------------------
#if 1
  /* tentative treatment for the case that y is almost zero and hence least squared method returns inf */
  if( isfinite(*bb) == 0 )
    *bb = 0.0;
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void getInterpolatedDensityProfile(const int num, profile * restrict prf, double * restrict xx, double * restrict ff)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* extrapolate for the innermost position by least squared method */
  /* extrapolate for the outermost position by least squared method */
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  double pp, bb;
  const double rmin = 0.5 * (prf[          0].rad < xx[          NPUT]) ? (prf[          0].rad) : (xx[          NPUT]);
  const double rmax = 2.0 * (prf[NRADBIN + 3].rad > xx[num - 1 - NPUT]) ? (prf[NRADBIN + 3].rad) : (xx[num - 1 - NPUT]);
  leastSquaredMethod(NFIT, &xx[                 NPUT], &ff[                 NPUT], &pp, &bb);
  const double logrmin = log10(rmin);
  double logrbin = (log10(xx[NPUT]) - logrmin) / (double)NPUT;
  for(int ii = 0; ii < NPUT; ii++){
    xx[ii] = pow(10.0, logrmin + logrbin * (double)ii);
    ff[ii] = bb * pow(xx[ii], pp);
  }/* for(int ii = 0; ii < NPUT; ii++){ */
  const double fpl = bb * pp * pow(rmin, pp - 1.0);
  leastSquaredMethod(NFIT, &xx[num - 1 - NFIT - NPUT], &ff[num - 1 - NFIT - NPUT], &pp, &bb);
#if 0
  fprintf(stderr, "NFIT = %d, NPUT = %d, pp = %e, bb = %e\n", NFIT, NPUT, pp, bb);
  exit(0);
#endif
  const double logrmax = log10(rmax);
  logrbin = (log10(xx[num - 1 - NPUT]) - logrmax) / (double)NPUT;
  for(int ii = 0; ii < NPUT; ii++){
    xx[num - 1 - ii] = pow(10.0, logrmax + logrbin * (double)ii);
    ff[num - 1 - ii] = bb * pow(xx[num - 1 - ii], pp);
  }/* for(int ii = 0; ii < NPUT; ii++){ */
  const double fpr = bb * pp * pow(rmax, pp - 1.0);
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < num; ii++)
    fprintf(stderr, "%e\t%e\n", xx[ii], ff[ii]);
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* generate table for cubic spline interpolation */
  //-----------------------------------------------------------------------
  double *bp;  bp = (double *)malloc(num * sizeof(double));  if( bp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate bp\n");  }
  double *f2;  f2 = (double *)malloc(num * sizeof(double));  if( f2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate f2\n");  }
  genCubicSpline1D(num, xx, ff, bp, fpl, fpr, f2);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* asign interpolated density profile */
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    prf[ii].  rho     = getCubicSpline1D               (prf[ii].rad, num, xx, ff, f2);
    prf[ii]. drho_dr  = getCubicSpline1stDifferential1D(prf[ii].rad, num, xx, ff, f2);
    prf[ii].d2rho_dr2 = getCubicSpline2ndDifferential1D(prf[ii].rad, num, xx,     f2);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < 4 + NRADBIN; ii++)
    fprintf(stderr, "%e\t%e\t%e\t%e\n", prf[ii].rad, prf[ii].rho, prf[ii].drho_dr, prf[ii].d2rho_dr2);
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(bp);  free(f2);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setDensityProfileTable(profile *prf, const double rs, char *file)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* read data table */
  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s", CFGFOLDER, file);
  fp = fopen(filename, "r");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  bool success = true;
  //-----------------------------------------------------------------------
  int num;
  if( 1 != fscanf(fp, "%d", &num) )    success = false;
  num += 2 * NPUT;
  //-----------------------------------------------------------------------
  double *xx;  xx = (double *)malloc(num * sizeof(double));  if( xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xx\n");  }
  double *ff;  ff = (double *)malloc(num * sizeof(double));  if( ff == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ff\n");  }
  //-----------------------------------------------------------------------
  for(int ii = NPUT; ii < num - NPUT; ii++)
    success &= (2 == fscanf(fp, "%le\t%le", &xx[ii], &ff[ii]));
  //-----------------------------------------------------------------------
  fclose(fp);
  if( success != true )
    writeDensityProfileTableFormat(filename);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* scale length to the computational unit */
  /* assume unity corresponds to the scale radius in the computational unit */
  //-----------------------------------------------------------------------
  for(int ii = NPUT; ii < num - NPUT; ii++)
    xx[ii] *= rs;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* asign interpolated density profile */
  //-----------------------------------------------------------------------
  getInterpolatedDensityProfile(num, prf, xx, ff);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(xx);  free(ff);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
