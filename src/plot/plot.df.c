/*************************************************************************\
 *                                                                       *
                  last updated on 2017/01/20(Fri) 14:54:45
 *                                                                       *
 *    Plot Code of Elapsed time of Tree code                             *
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
#include <hdf5.h>
#include "hdf5lib.h"
//-------------------------------------------------------------------------
#include "macro.h"
#include "name.h"
#include "plplotlib.h"
#include "constants.h"
//-------------------------------------------------------------------------
#include "../init/magi.h"
#include "../init/eddington.h"
//-------------------------------------------------------------------------
/* extern const char senergy_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS]; */
/* extern double senergy2astro; */
/* extern double senergy2cgs; */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define LOGPLOT_HOR
#define LOGPLOT_VER
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define NMAX (4096)
/* #define NMAX (8192) */
/* #define NMAX (16384) */
/* #define NMAX (32768) */
/* #define NMAX (65536) */
//-------------------------------------------------------------------------
#define NDST (2)
static char datname[NDST][16] = {"hernquist", "plummer"};
static char  icname[NDST][16] = {"Hernquist", "Plummer"};
//-------------------------------------------------------------------------
/* #define NSER (2) */
#define NSER (3)
static char sername[NSER][2] = {"a", "t"
#   if  NSER > 2
				 , "x"
#endif//NSER > 2
};
static char synonym[NSER][16] = {"function", "table"
#   if  NSER > 2
				 , "analytic"
#endif//NSER > 2
};
//-------------------------------------------------------------------------
#define NPLT (2)
/* DF and relative error */
//-------------------------------------------------------------------------
#define NERR (1)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define NORMALIZE_VARIABLES
//-------------------------------------------------------------------------
#ifdef  NORMALIZE_VARIABLES
static double Enorm[NDST];
#endif//NORMALIZE_VARIABLES
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void readDF(char *file, PLFLT ene[restrict], PLFLT val[restrict], PLFLT rad[restrict]);
void compareDF(const bool energy, PLFLT ene[restrict][NSER][NMAX], PLFLT val[restrict][NSER][NMAX], PLFLT err[restrict][NERR][NMAX], PLplotPltRange box[restrict][NPLT], int argc, char **argv);
//-------------------------------------------------------------------------
void calcPlummerDF  (PLFLT ene[restrict], PLFLT val[restrict], PLFLT rad[restrict]);
void calcHernquistDF(PLFLT ene[restrict], PLFLT val[restrict], PLFLT rad[restrict]);
//-------------------------------------------------------------------------
#ifndef NORMALIZE_VARIABLES
void rescaleDF2CGS(PLFLT val[restrict][NSER][NMAX]);
#endif//NORMALIZE_VARIABLES
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //-----------------------------------------------------------------------
  /* initialization */
  //-----------------------------------------------------------------------
  if( argc < 1 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 1);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 1 ){ */
  //-----------------------------------------------------------------------
  modifyArgcArgv4PLplot(&argc, argv, 1);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* prepare dataset */
  static PLFLT ene[NDST][NSER][NMAX];
  static PLFLT val[NDST][NSER][NMAX];
  static PLFLT rad[NDST][NSER][NMAX];
  //-----------------------------------------------------------------------
  /* read DF generated by MAGI */
  int unit = -1;
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < ((NSER <= 2) ? NSER : 2); jj++){
      //-------------------------------------------------------------------
      char filename[128];
      sprintf(filename, "%s%s", sername[jj], datname[ii]);
      //-------------------------------------------------------------------
      /* verify unit system */
      char summary[128];
      sprintf(summary, "%s/%s.summary.txt", DOCUMENTFOLDER, filename);
      FILE *fp;
      fp = fopen(summary, "r");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", summary);      }
      int tmp;
      fscanf(fp, "%d\n", &tmp);
      fclose(fp);
      //-------------------------------------------------------------------
      if( (ii + jj) == 0 ){
	unit = tmp;
	setPhysicalConstantsAndUnitSystem(unit, 0);
      }/* if( (ii + jj) == 0 ){ */
      else
	if( tmp != unit ){
	  __KILL__(stderr, "ERROR: unit system does not match: tmp = %d while unit = %d\n", tmp, unit);
	}/* if( tmp != unit ){ */
      //-------------------------------------------------------------------
      readDF(filename, ene[ii][jj], val[ii][jj], rad[ii][jj]);
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < ((NSER <= 2) ? NSER : 2); jj++){ */
  //-----------------------------------------------------------------------
#   if  NSER > 2
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NMAX; jj++)
      ene[ii][NSER - 1][jj] = ene[ii][0][jj];
  calcHernquistDF(ene[0][NSER - 1], val[0][NSER - 1], rad[0][NSER - 1]);
  calcPlummerDF  (ene[1][NSER - 1], val[1][NSER - 1], rad[1][NSER - 1]);
#endif//NSER > 2
  //-----------------------------------------------------------------------
  static PLFLT err[NDST][NERR][NMAX];
  for(int ii = 0; ii < NDST; ii++)
    for(int kk = 0; kk < NMAX; kk++){
      //-------------------------------------------------------------------
      const PLFLT ans = val[ii][0][kk];
      const PLFLT inv = 1.0 / (DBL_MIN + ans);
      //-------------------------------------------------------------------
      for(int jj = 0; jj < NERR; jj++)
	err[ii][jj][kk] = fabs(val[ii][jj + 1][kk] - ans) * inv;
      //-------------------------------------------------------------------
    }/* for(int kk = 0; kk < NMAX; kk++){ */
  //-----------------------------------------------------------------------
#ifndef NORMALIZE_VARIABLES
  rescaleDF2CGS(val);
#endif//NORMALIZE_VARIABLES
  //-----------------------------------------------------------------------
#ifdef  NORMALIZE_VARIABLES
#if 0
  for(int ii = 0; ii < NDST; ii++)
    printf("Enorm[%d] = %e\n", ii, Enorm[ii]);
#endif
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NSER; jj++)
      for(int kk = 0; kk < NMAX; kk++)
	ene[ii][jj][kk] *= Enorm[ii];
#endif//NORMALIZE_VARIABLES
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* plot elapsed time as a function of total number of N-body particles */
  //-----------------------------------------------------------------------
  PLplotPltRange range[NDST][NPLT];
  for(int ii = 0; ii < NDST; ii++){
    //---------------------------------------------------------------------
    range[ii][0].xmin = ene[ii][0][       0];
    range[ii][0].xmax = ene[ii][0][NMAX - 1];
    //---------------------------------------------------------------------
    range[ii][0].ymin = DBL_MAX;
    range[ii][0].ymax = 0.0;
    for(int jj = 0; jj < NMAX; jj++)
      /* if( val[ii][0][jj] > DBL_EPSILON ){ */
      if( val[ii][0][jj] > 0.0 ){
	range[ii][0].ymax = fmax(val[ii][0][jj], range[ii][0].ymax);
	range[ii][0].ymin = fmin(val[ii][0][jj], range[ii][0].ymin);
      }/* if( val[ii][0][jj] > 0.0 ){ */
    range[ii][1].ymin = DBL_MAX;
    range[ii][1].ymax = 0.0;
    for(int jj = 0; jj < NERR; jj++)
      for(int kk = 0; kk < NMAX; kk++)
	if( err[ii][jj][kk] > DBL_EPSILON ){
	  range[ii][1].ymax = fmax(err[ii][jj][kk], range[ii][1].ymax);
	  range[ii][1].ymin = fmin(err[ii][jj][kk], range[ii][1].ymin);
	}/* if( err[ii][jj][kk] > DBL_EPSILON ){ */
    //---------------------------------------------------------------------
#ifdef  LOGPLOT_HOR
    range[ii][0].xmin = log10(range[ii][0].xmin);
    range[ii][0].xmax = log10(range[ii][0].xmax);
    range[ii][0].xlog = LOGARITHMIC_PLOT;
#else///LOGPLOT_HOR
    range[ii][0].xlog =      LINEAR_PLOT;
#endif//LOGPLOT_HOR
    //---------------------------------------------------------------------
#if 1
    range[ii][1].ymin = fmax(range[ii][1].ymin, FLT_EPSILON);
    range[ii][1].ymax = fmin(range[ii][1].ymax, 10.0);
#endif
    //---------------------------------------------------------------------
#ifdef  LOGPLOT_VER
    range[ii][0].ymin = log10(range[ii][0].ymin);    range[ii][1].ymin = log10(range[ii][1].ymin);
    range[ii][0].ymax = log10(range[ii][0].ymax);    range[ii][1].ymax = log10(range[ii][1].ymax);
#if 1
    range[ii][0].ymin = fmax(range[ii][0].ymin, range[ii][0].ymax - 12.0);
#endif
    range[ii][0].ylog = LOGARITHMIC_PLOT;
    range[ii][1].ylog = LOGARITHMIC_PLOT;
#else///LOGPLOT_VER
    range[ii][0].ylog =      LINEAR_PLOT;
    range[ii][1].ylog =      LINEAR_PLOT;
#endif//LOGPLOT_VER
    //---------------------------------------------------------------------
#if 1
    const PLFLT Lx = range[ii][0].xmax - range[ii][0].xmin;    range[ii][0].xmin -= 0.125 * Lx;    range[ii][0].xmax += 0.125 * Lx;
#ifdef  LOGPLOT_VER
    if( range[ii][0].ymax > log10(val[ii][0][NMAX - 1]) )
#else///LOGPLOT_VER
    if( range[ii][0].ymax >       val[ii][0][NMAX - 1]  )
#endif//LOGPLOT_VER
      {
	const PLFLT Ly = range[ii][0].ymax - range[ii][0].ymin;    range[ii][0].ymin -= 0.125 * Ly;    range[ii][0].ymax += 0.125 * Ly;
      }
    const PLFLT Ey = range[ii][1].ymax - range[ii][1].ymin;    range[ii][1].ymin -= 0.125 * Ey;    range[ii][1].ymax += 0.125 * Ey;
/* #ifndef LOGPLOT_HOR */
/*     range[ii][0].xmin = fmax(range[ii][0].xmin, 0.0); */
/* #endif//LOGPLOT_HOR */
#endif
    //---------------------------------------------------------------------
    range[ii][0].xgrd = range[ii][1].xgrd = range[ii][0].ygrd = range[ii][1].ygrd = true;
    //---------------------------------------------------------------------
    range[ii][1].xmin = range[ii][0].xmin;
    range[ii][1].xmax = range[ii][0].xmax;
    range[ii][1].xlog = range[ii][0].xlog;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NDST; ii++){ */
  //-----------------------------------------------------------------------
#if 1
  range[0][0].ymin -= 2.0;
  range[0][0].ymax -= 2.0;
  range[0][1].ymin = log10(FLT_EPSILON);
  range[1][1].ymin = log10(FLT_EPSILON);
  range[1][0].ymin = range[0][0].ymin;  range[1][0].ymax = range[0][0].ymax;
  range[0][1].ymin = range[1][1].ymin;  range[0][1].ymax = range[1][1].ymax;
#endif
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NPLT; jj++)
      fprintf(stdout, "range[%d][%d]: (%le, %le), (%le, %le)\n", ii, jj, range[ii][jj].xmin, range[ii][jj].xmax, range[ii][jj].ymin, range[ii][jj].ymax);
  fflush(stdout);
  exit(0);
#endif
  //-----------------------------------------------------------------------
#ifdef  LOGPLOT_HOR
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NSER; jj++)
      for(int kk = 0; kk < NMAX; kk++)
	ene[ii][jj][kk] = (ene[ii][jj][kk] > DBL_MIN) ? log10(ene[ii][jj][kk]) : DBL_MAX;
#endif//LOGPLOT_HOR
#ifdef  LOGPLOT_VER
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NSER; jj++)
      for(int kk = 0; kk < NMAX; kk++)
	val[ii][jj][kk] = (val[ii][jj][kk] > DBL_MIN) ? log10(val[ii][jj][kk]) : DBL_MAX;
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NERR; jj++)
      for(int kk = 0; kk < NMAX; kk++)
	err[ii][jj][kk] = (err[ii][jj][kk] > DBL_MIN) ? log10(err[ii][jj][kk]) : DBL_MAX;
#endif//LOGPLOT_VER
  //-----------------------------------------------------------------------
  compareDF(true, ene, val, err, range, argc, argv);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  LOGPLOT_HOR
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NPLT; jj++){
      //-------------------------------------------------------------------
      range[ii][jj].xmin = log10(1.0e-2);
      range[ii][jj].xmax = log10(9.9e+1);
      range[ii][jj].xlog = LOGARITHMIC_PLOT;
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < NPLT; jj++){ */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NSER; jj++)
      for(int kk = 0; kk < NMAX; kk++)
	rad[ii][jj][kk] = (rad[ii][jj][kk] > DBL_MIN) ? log10(rad[ii][jj][kk]) : DBL_MAX;
  //-----------------------------------------------------------------------
#else///LOGPLOT_HOR
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NPLT; jj++){
      //-------------------------------------------------------------------
      range[ii][jj].xmin = -5.0;
      range[ii][jj].xmax = 30.0;
      range[ii][jj].xlog = LINEAR_PLOT;
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < NPLT; jj++){ */
  //-----------------------------------------------------------------------
#endif//LOGPLOT_HOR
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NDST; ii++){
    range[ii][0].ymin = -4.0;
    range[ii][0].ymax =  4.0;
    range[ii][1].ymin = -6.0;
    range[ii][1].ymax = log10(9.9e-2);
  }/* for(int ii = 0; ii < NDST; ii++){ */
  //-----------------------------------------------------------------------
  compareDF(false, rad, val, err, range, argc, argv);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline int bisection(const real val, const int num, real * restrict tab, real * restrict ratio)
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = num - 1;
  //-----------------------------------------------------------------------
  /* prohibit extraporation */
  if( val > tab[ll] + EPSILON ){    *ratio =  ZERO;    return (ll    );  }
  if( val < tab[rr] - EPSILON ){    *ratio = UNITY;    return (rr - 1);  }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)(ll + rr)) >> 1;
    //---------------------------------------------------------------------
    if( (tab[cc] - val) * (tab[ll] - val) <= ZERO)      rr = (int)cc;
    else                                                ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr ){
      *ratio = (val - tab[ll]) / (tab[rr] - tab[ll]);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void readDF(char *file, PLFLT ene[restrict], PLFLT val[restrict], PLFLT rad[restrict])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  static real tmp_ene[NENEBIN], tmp_val[NENEBIN];
  static real tmp_rad[NRADBIN], tmp_psi[NRADBIN];
  //-----------------------------------------------------------------------
#ifdef  DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION
  //-----------------------------------------------------------------------
  /* open an existing file for DF */
  char filename[128];
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "df");
  hid_t dftarget = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  //-----------------------------------------------------------------------
  /* open an existing file for profile */
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "profile");
  hid_t rrtarget = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  //-----------------------------------------------------------------------
  /* read attributes */
  hid_t attribute = H5Aopen(dftarget, "kinds", H5P_DEFAULT);
  int kind;
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &kind));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* read data */
  for(int kk = 0; kk < kind; kk++){
    //---------------------------------------------------------------------
    char grp[16];    sprintf(grp, "series%d", kk);
    hid_t dfgroup = H5Gopen(dftarget, grp, H5P_DEFAULT);
    sprintf(grp, "data%d", kk);
    hid_t rrgroup = H5Gopen(rrtarget, grp, H5P_DEFAULT);
    //---------------------------------------------------------------------
    hid_t dataset;
    /* read energy */
    dataset = H5Dopen(dfgroup, "energy", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_ene));
    chkHDF5err(H5Dclose(dataset));
    /* read DF */
    dataset = H5Dopen(dfgroup, "DF", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_val));
    chkHDF5err(H5Dclose(dataset));
    /* read radius */
    dataset = H5Dopen(rrgroup, "rad", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rad));
    chkHDF5err(H5Dclose(dataset));
    /* read Psi */
    dataset = H5Dopen(rrgroup, "Psi", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_psi));
    chkHDF5err(H5Dclose(dataset));
    //---------------------------------------------------------------------
    /* read attributes */
    double rs;
    attribute = H5Aopen(rrgroup, "rs", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &rs));
    chkHDF5err(H5Aclose(attribute));
#ifdef  NORMALIZE_VARIABLES
    double Mtot;
    attribute = H5Aopen(rrgroup, "Mtot", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &Mtot));
    chkHDF5err(H5Aclose(attribute));
#endif//NORMALIZE_VARIABLES
    //---------------------------------------------------------------------
    /* close the group */
    chkHDF5err(H5Gclose(dfgroup));
    chkHDF5err(H5Gclose(rrgroup));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* save data */
    const int skip = NENEBIN / NMAX;
    //---------------------------------------------------------------------
    const PLFLT rsinv = 1.0 / rs;
    for(int ii = 0; ii < NMAX; ii++){
      real ratio;
      const int idx = bisection(tmp_ene[ii * skip], NRADBIN, tmp_psi, &ratio);
      rad[ii] = rsinv * (PLFLT)(((UNITY - ratio) * tmp_rad[idx] + ratio * tmp_rad[idx + 1]));
      /* fprintf(stdout, "%e\t%d\t%le\n", tmp_ene[ii], idx, rad[ii]); */
    }/* for(int ii = 0; ii < NMAX; ii++){ */
    /* exit(0); */
    //---------------------------------------------------------------------
#ifdef  NORMALIZE_VARIABLES
    extern real newton;
    const double G = (double)newton;
    const double scale = sqrt((4.0 * M_PI * G) * (4.0 * M_PI * G) * (4.0 * M_PI * G) * Mtot * rs * rs * rs);
#endif//NORMALIZE_VARIABLES
    extern double senergy2astro;
    extern double senergy2cgs;
    const double senergyAstro2CGS = senergy2cgs / senergy2astro;
#if 0
    fprintf(stdout, "senergyAstro2CGS = %e; senergy2cgs = %e, senergy2astro = %e\n", senergyAstro2CGS, senergy2cgs, senergy2astro);
#endif
    for(int ii = 0; ii < NMAX; ii++){
      //-------------------------------------------------------------------
      ene[ii] = (PLFLT)tmp_ene[ii * skip] * senergyAstro2CGS;
#ifdef  NORMALIZE_VARIABLES
      val[ii] = (PLFLT)tmp_val[ii * skip] * scale;
#else///NORMALIZE_VARIABLES
      val[ii] = (PLFLT)tmp_val[ii * skip];
#endif//NORMALIZE_VARIABLES
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NMAX; ii++){ */
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < kind; kk++){ */
  //-----------------------------------------------------------------------
  /* close the file */
  chkHDF5err(H5Fclose(dftarget));
  chkHDF5err(H5Fclose(rrtarget));
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void compareDF
(const bool energy, PLFLT ene[restrict][NSER][NMAX], PLFLT val[restrict][NSER][NMAX], PLFLT err[restrict][NERR][NMAX], PLplotPltRange box[restrict][NPLT], int argc, char **argv)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* global setting(s) */
  //-----------------------------------------------------------------------
  /* maximum number of data kinds */
  const PLINT pkind = NERR;
  const PLINT lkind = NSER;
  //-----------------------------------------------------------------------
  const PLBOOL puni = true;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* preparation for dataset */
  //-----------------------------------------------------------------------
  /* memory allocation for index */
  PLINT kind = (puni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  PLINT *num;  allocPLINT(&num, kind);
  //-----------------------------------------------------------------------
  /* set number of data points */
  for(PLINT ii = 0; ii < lkind; ii++)    num[        ii] = NMAX;
  for(PLINT ii = 0; ii < pkind; ii++)    num[lkind + ii] = NMAX;
  //-----------------------------------------------------------------------
  /* set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  ls[0].style = DOTTED_LINE;  ls[0].width = MIDDLE_LINE;  ls[0].color = RED;
  ls[1].style = SOLID_LINE;  ls[1].width = MIDDLE_LINE;  ls[1].color = BLACK;
#   if  NSER > 2
  ls[2].style = DASHED_LINE;  ls[2].width = MIDDLE_LINE;  ls[2].color = BLUE;
#endif//NSER > 2
  sprintf(pt[0].type, PLplotSymbolType[smallDot]);  pt[0].scale = PLplotSymbolSize[smallDot];  pt[0].color = BLACK;
#   if  NERR > 1
  sprintf(pt[1].type, PLplotSymbolType[smallDot]);  pt[1].scale = PLplotSymbolSize[smallDot];  pt[1].color = RED;
#endif//NERR > 1
  //-----------------------------------------------------------------------
  /* set labels */
  char hlab[PLplotCharWords], vlab[PLplotCharWords], elab[PLplotCharWords];
  if( energy ){
#ifdef  NORMALIZE_VARIABLES
    sprintf(hlab, "-#fiE#fr / (#fiG M#fr#dtot#u / #fir#fr#ds#u)");
#else///NORMALIZE_VARIABLES
    sprintf(hlab, "#fi#[0x2130]#fr (erg g#u-1#d)");
#endif//NORMALIZE_VARIABLES
  }/* if( energy ){ */
  else
    sprintf(hlab, "#fir#fr#dout#u#fi#fr(#fiE#fr) / #fir#fr#ds#u");
#ifdef  NORMALIZE_VARIABLES
  sprintf(vlab, "(#fi4#gpG r#fr#ds#u)#u3/2#d #fiM#fr#dtot#u#u1/2#d #fif#fr (#fiE#fr)");
#else///NORMALIZE_VARIABLES
  sprintf(vlab, "#fif#fr (g cm#u-5#d s#u2#d)");
#endif//NORMALIZE_VARIABLES
  sprintf(elab, "Relative difference");
  //-----------------------------------------------------------------------
  /* set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write = true;
  if( energy )    sprintf(basecap.side, "%s", "t");
  else            sprintf(basecap.side, "%s", "b");
  //-----------------------------------------------------------------------
  /* set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = true;
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
  for(int jj = 0; jj < NSER; jj++)
    sprintf(pleg.text[jj], "%s", synonym[jj]);
  if( energy )
    pleg.pos = PL_POSITION_RIGHT | PL_POSITION_BOTTOM | PL_POSITION_INSIDE;
  else
    pleg.pos = PL_POSITION_RIGHT | PL_POSITION_TOP    | PL_POSITION_INSIDE;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation to exploit multiple-plot function */
  //-----------------------------------------------------------------------
  /* maximum number of panels in each direction */
  const PLINT nxpanel = NDST;
  const PLINT nypanel = NPLT;
  //-----------------------------------------------------------------------
  /* memory allocation for data */
  PLFLT **hor;  allocPointer4PLFLT(&hor, nxpanel * nypanel * kind);
  PLFLT **ver;  allocPointer4PLFLT(&ver, nxpanel * nypanel * kind);
  //-----------------------------------------------------------------------
  /* specify plot or skip panel */
  PLBOOL *plot;  allocPLBOOL(&plot, nxpanel * nypanel);
  /* arrays related to draw line(s) and point(s) */
  PLINT *nlkind;  allocPLINT(&nlkind, nxpanel * nypanel);
  PLINT *npkind;  allocPLINT(&npkind, nxpanel * nypanel);
  PLINT **lnum;  allocPointer4PLINT(&lnum, nxpanel * nypanel);
  PLINT **pnum;  allocPointer4PLINT(&pnum, nxpanel * nypanel);
  PLFLT ***lx;  allocDoublePointer4PLFLT(&lx, nxpanel * nypanel);
  PLFLT ***ly;  allocDoublePointer4PLFLT(&ly, nxpanel * nypanel);
  PLFLT ***px;  allocDoublePointer4PLFLT(&px, nxpanel * nypanel);
  PLFLT ***py;  allocDoublePointer4PLFLT(&py, nxpanel * nypanel);
  PLplotLineStyle ** line;  allocPointer4PLplotLineStyle(& line, nxpanel * nypanel);
  PLplotPointType **point;  allocPointer4PLplotPointType(&point, nxpanel * nypanel);
  /* arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /* arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /* arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /* arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /* arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /* array to set figure name */
  char figfile[PLplotCharWords];
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel);
  assignChar4PLplot(nxpanel, figname, _figname);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* configure to plot DF */
  //-----------------------------------------------------------------------
  /* data preparation */
  for(int ii = 0; ii < nxpanel; ii++){
    for(int jj = 0; jj < lkind; jj++){
      hor[INDEX(nxpanel, nypanel, kind, ii, NPLT - 1, jj)] = ene[ii][jj];
      ver[INDEX(nxpanel, nypanel, kind, ii, NPLT - 1, jj)] = val[ii][jj];
    }/* for(int jj = 0; jj < NSER; jj++){ */
    for(int jj = 0; jj < pkind; jj++){
      hor[INDEX(nxpanel, nypanel, kind, ii,        0, jj)] = ene[ii][jj + 1];
      ver[INDEX(nxpanel, nypanel, kind, ii,        0, jj)] = err[ii][jj    ];
    }/* for(int jj = 0; jj < NERR; jj++){ */
  }/* for(int ii = 0; ii < nxpanel; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      //-------------------------------------------------------------------
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      //-------------------------------------------------------------------
      /* global setting(s) */
      plot[idx] = true;
      //-------------------------------------------------------------------
      /* errorbar setting(s) */
      errbar[idx] = false;
      //-------------------------------------------------------------------
      /* label setting(s) */
      xlabel[idx] = hlab;
      //-------------------------------------------------------------------
      /* misc(s) */
      leg[idx] = pleg;
      uni[idx] = false;
      cap[idx] = basecap;
      //-------------------------------------------------------------------
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */
  //-----------------------------------------------------------------------
  for(PLINT ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT jj = 0;
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* line setting(s) */
    nlkind[idx] = 0;
    line  [idx] = NULL;
    lnum  [idx] = NULL;
    lx    [idx] = NULL;
    ly    [idx] = NULL;
    //---------------------------------------------------------------------
    /* point setting(s) */
    npkind[idx] = pkind;
    point [idx] = pt;
    pnum  [idx] = &num[lkind];
    px    [idx] = &hor[INDEX2D(nxpanel * nypanel, kind, idx, 0)];
    py    [idx] = &ver[INDEX2D(nxpanel * nypanel, kind, idx, 0)];
    //---------------------------------------------------------------------
    /* plot area */
    range[idx] = box[ii][1];
    //---------------------------------------------------------------------
    /* label setting(s) */
    ylabel[idx] = elab;
    //---------------------------------------------------------------------
    /* misc(s) */
    leg[idx].write = false;
    //---------------------------------------------------------------------
  }/* for(PLINT ii = 0; ii < nxpanel; ii++){ */
  //-----------------------------------------------------------------------
  for(PLINT ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT jj = NPLT - 1;
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* line setting(s) */
    nlkind[idx] = lkind;
    line  [idx] = ls;
    lnum  [idx] = num;
    lx    [idx] = &hor[INDEX2D(nxpanel * nypanel, kind, idx, 0)];
    ly    [idx] = &ver[INDEX2D(nxpanel * nypanel, kind, idx, 0)];
    //---------------------------------------------------------------------
    /* point setting(s) */
    npkind[idx] = 0;
    point [idx] = NULL;
    pnum  [idx] = NULL;
    px    [idx] = NULL;
    py    [idx] = NULL;
    //---------------------------------------------------------------------
    /* plot area */
    range[idx] = box[ii][0];
    //---------------------------------------------------------------------
    /* label setting(s) */
    ylabel[idx] = vlab;
    //---------------------------------------------------------------------
  }/* for(PLINT ii = 0; ii < nxpanel; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* individual file names */
  if( energy ){
    sprintf(figfile, "%s", "df");
    for(int ii = 0; ii < nxpanel; ii++)
      sprintf(figname[ii], "%s_%s", datname[ii], "df");
  }/* if( energy ){ */
  else{
    sprintf(figfile, "%s", "df_rad");
    for(int ii = 0; ii < nxpanel; ii++)
      sprintf(figname[ii], "%s_%s", datname[ii], "df_rad");
  }/* else{ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create figure(s) */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 0);
    /* 97 is "a" in ASCII code */
    sprintf(cap[INDEX2D(nxpanel, nypanel, ii, 1)].text, "(%c) %s", 97, "DF");
    sprintf(cap[INDEX2D(nxpanel, nypanel, ii, 0)].text, "(%c) %s", 98, "Error");
    //---------------------------------------------------------------------
    plotData(1, nypanel, &plot[idx], true, false,
  	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
  	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
  	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
  	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
  	     &xlabel[idx], &ylabel[idx], "", figname[ii], argc, argv);
    //---------------------------------------------------------------------
  }/* for(int jj = 0; jj < nypanel; jj++){ */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      if( idx != INDEX2D(nxpanel, nypanel, 0, nypanel - 1) )
  	leg[idx].write = false;
      /* 97 is "a" in ASCII code */
      sprintf(cap[idx].text, "(%c) %s (%s)", 97 + INDEX2D(nypanel, nxpanel, nypanel - jj - 1, ii), (jj != 0) ? "DF" : "Error", icname[ii]);
    }/* for(int jj = 0; jj < nypanel; jj++){ */
  plotData(nxpanel, nypanel, plot, true, true,
  	   nlkind,  line, lnum, lx, ly,
  	   npkind, point, pnum, px, py,
  	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
  	   cap, leg, uni, range,
  	   xlabel, ylabel, "", figfile, argc, argv);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  free(figname);  free(_figname);
  //-----------------------------------------------------------------------
  free(hor);  free(ver);
  //-----------------------------------------------------------------------
  free(num);
  free(ls);  free(pt);
  free(plegTxt);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void calcPlummerDF(PLFLT ene[restrict], PLFLT val[restrict], PLFLT rad[restrict])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  extern real newton;
  extern double mass2com;
  extern double length2com;
  extern const double kpc;
  extern const double Msun;
  extern double senergy2com;
  //-----------------------------------------------------------------------
  const double Mtot = 1.0e+9 * Msun * mass2com;
  const double rs = 2.0 * kpc * length2com;
  const double G = (double)newton;
#ifdef  NORMALIZE_VARIABLES
  Enorm[1] = senergy2com * rs / (G * Mtot);
#endif//NORMALIZE_VARIABLES
  const double foo = G * G * Mtot * Mtot / (rs * rs);
  //-----------------------------------------------------------------------
#ifdef  NORMALIZE_VARIABLES
  const double factor = sqrt((4.0 * M_PI * G * rs) * (4.0 * M_PI * G * rs) * (4.0 * M_PI * G * rs) * Mtot) * 24.0 * M_SQRT2 * rs * rs / (7.0 * M_PI * M_PI * M_PI * G * G * G * G * G * Mtot * Mtot * Mtot * Mtot);
#else///NORMALIZE_VARIABLES
  const double factor = 24.0 * M_SQRT2 * rs * rs / (7.0 * M_PI * M_PI * M_PI * G * G * G * G * G * Mtot * Mtot * Mtot * Mtot);
#endif//NORMALIZE_VARIABLES
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NMAX; ii++){
    //---------------------------------------------------------------------
    const double tmp = ene[ii] * senergy2com;
    //---------------------------------------------------------------------
    val[ii] = factor * tmp * tmp * tmp * sqrt(tmp);
    rad[ii] = (foo > tmp * tmp) ? (sqrt(foo / (tmp * tmp) - 1.0)) : (0.0);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NMAX; ii++){ */
  //-----------------------------------------------------------------------
#if 0
  fprintf(stdout, "# Plummer profile\n");
  for(int ii = 0; ii < NMAX; ii++)
    fprintf(stdout, "%e\t%e\n", ene[ii], val[ii]);
  exit(0);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void calcHernquistDF(PLFLT ene[restrict], PLFLT val[restrict], PLFLT rad[restrict])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  extern real newton;
  extern double mass2com;
  extern double length2com;
  extern const double kpc;
  extern const double Msun;
  extern double senergy2com;
  //-----------------------------------------------------------------------
  const double Mtot = 1.0e+10 * Msun * mass2com;
  const double rs = 1.0 * kpc * length2com;
  const double G = (double)newton;
#ifdef  NORMALIZE_VARIABLES
  Enorm[0] = senergy2com * rs / (G * Mtot);
#endif//NORMALIZE_VARIABLES
  //-----------------------------------------------------------------------
  const double tmp = sqrt(rs / (G * Mtot));
#ifdef  NORMALIZE_VARIABLES
  const double factor = sqrt((4.0 * M_PI * G) * (4.0 * M_PI * G) * (4.0 * M_PI * G) * Mtot * rs * rs * rs) * Mtot * tmp * tmp * tmp / (8.0 * M_SQRT2 * M_PI * M_PI * M_PI * rs * rs * rs);
#else///NORMALIZE_VARIABLES
  const double factor = Mtot * tmp * tmp * tmp / (8.0 * M_SQRT2 * M_PI * M_PI * M_PI * rs * rs * rs);
#endif//NORMALIZE_VARIABLES
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NMAX; ii++){
    //---------------------------------------------------------------------
    const double qq = tmp * sqrt(senergy2com * ene[ii]);
    const double q2 = qq * qq;
    const double foo = sqrt(1.0 - q2);
    //---------------------------------------------------------------------
    val[ii] = (q2 < 1.0) ? (factor * (3.0 * asin(qq) + qq * foo * (1.0 - 2.0 * q2) * (8.0 * (q2 - 1.0) * q2 - 3.0)) / (foo * (1.0 - q2) * (1.0 - q2))) : (0.0);
    rad[ii] = (q2 < 1.0) ? (1.0 / q2 - 1.0) : 0.0;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NMAX; ii++){ */
  //-----------------------------------------------------------------------
#if 0
  fprintf(stdout, "# Hernquist profile\n");
  for(int ii = 0; ii < NMAX; ii++)
    fprintf(stdout, "%e\t%e\n", ene[ii], val[ii]);
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef NORMALIZE_VARIABLES
//-------------------------------------------------------------------------
void rescaleDF2CGS(PLFLT val[restrict][NSER][NMAX])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  extern double density2cgs;
  extern double senergy2cgs;
  const double df2cgs = density2cgs / senergy2cgs;
#if 0
  printf("density2cgs = %e, senergy2cgs = %e, df2cgs = %e\n", density2cgs, senergy2cgs, df2cgs);
  exit(0);
#endif
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NSER; jj++)
      for(int kk = 0; kk < NMAX; kk++)
	val[ii][jj][kk] *= df2cgs;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//NORMALIZE_VARIABLES
//-------------------------------------------------------------------------
