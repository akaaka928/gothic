/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:37:35
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
#include <string.h>
//-------------------------------------------------------------------------
#include "macro.h"
#include "myutil.h"
#include "plplotlib.h"
#include "name.h"
#include "constants.h"
//-------------------------------------------------------------------------
extern const double length2astro;
extern const char   length_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define MONOCHROME_PLOT */
//-------------------------------------------------------------------------
#define NORMALIZE_COL
//-------------------------------------------------------------------------
#define LOGPLOT_HOR
/* #define LOGPLOT_VER */
#define LOGPLOT_COL
//-------------------------------------------------------------------------
#define NKIND (5)
static char ballname[NKIND][16] = {"fischer03", "ritter90", "cartesian", "com", "smaller"};
/* static char synonym[NKIND - 1][32] = {"Ritter (1990)", "Rectangular cuboid", "Center of mass", "Cuboid with C.O.M."}; */
static char synonym[NKIND - 1][4] = {"EBS", "GEO", "COM", "CMP"};
//-------------------------------------------------------------------------
/* #define NX (32) */
/* #define NY (32) */
//-------------------------------------------------------------------------
#define NX (64)
#define NY (64)
//-------------------------------------------------------------------------
/* #define NX (128) */
/* #define NY (128) */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotElapsedTime
(const int Ntot, PLFLT * restrict rad, PLFLT * restrict rate, PLplotPltRange box, char file[], int argc, char **argv);
void plotDistribution
(PLFLT * restrict xmap, PLFLT * restrict ymap, PLFLT * restrict fmap, const PLFLT fmapMin, const PLFLT fmapMax, PLplotPltRange range, char file[], int argc, char **argv);
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //-----------------------------------------------------------------------
  /* initialization */
  //-----------------------------------------------------------------------
  if( argc < 2 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 2);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 2 ){ */
  //-----------------------------------------------------------------------
  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "file", &file));
  //-----------------------------------------------------------------------
  modifyArgcArgv4PLplot(&argc, argv, 2);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* prepare dataset */
  //-----------------------------------------------------------------------
  static int Nball[NKIND];
  {
    //---------------------------------------------------------------------
    /* read # of balls */
    //---------------------------------------------------------------------
    char filename[128];
    sprintf(filename, "%s/%s.ball.info.txt", LOGFOLDER, file);
    FILE *fp;
    fp = fopen(filename, "r");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    //---------------------------------------------------------------------
    int Ntot;
    char readname[128], testname[128];
    while( EOF != fscanf(fp, "%s\t%d", readname, &Ntot) ){
      //-------------------------------------------------------------------
      for(int ii = 0; ii < NKIND; ii++){
	sprintf(testname, "%s/%s.ball.%s.dat", DATAFOLDER, file, ballname[ii]);
	if( strcmp(testname, readname) == 0 )
	  Nball[ii] = Ntot;
      }/* for(int ii = 0; ii < NKIND; ii++){ */
      //-------------------------------------------------------------------
    }/* while( EOF != fscanf(fp, "%s\t%d", tmpname, &Nball) ){ */
    //---------------------------------------------------------------------
    fclose(fp);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* read and set unit system */
    //---------------------------------------------------------------------
    sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
    fp = fopen(filename, "r");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    int unit;
    fscanf(fp, "%d\n", &unit);
    fclose(fp);
    //---------------------------------------------------------------------
    setPhysicalConstantsAndUnitSystem(unit, 0);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#if 0
  fprintf(stdout, "length2astro = %e\n", length2astro);
  exit(0);
#endif
  //-----------------------------------------------------------------------
  const int Ntot = Nball[0];
  real  *ball;  ball = (real  *)malloc((size_t) NKIND      * (size_t)Ntot * (size_t)4 * sizeof( real));
  PLFLT * rad;  rad  = (PLFLT *)malloc((size_t) NKIND      * (size_t)Ntot             * sizeof(PLFLT));
  PLFLT *rate;  rate = (PLFLT *)malloc((size_t)(NKIND - 1) * (size_t)Ntot             * sizeof(PLFLT));
  if( ball == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ball\n");  }
  if(  rad == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  rad\n");  }
  if( rate == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rate\n");  }
  //-----------------------------------------------------------------------
  /* read the properties of enclosing balls */
  for(int ii = 0; ii < NKIND; ii++){
    //---------------------------------------------------------------------
    FILE *fp;
    char filename[128];
    sprintf(filename, "%s/%s.ball.%s.dat", DATAFOLDER, file, ballname[ii]);
    fp = fopen(filename, "rb");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    //---------------------------------------------------------------------
    fread(&ball[Ntot * 4 * ii], sizeof(real), Ntot * 4, fp);
    //---------------------------------------------------------------------
    fclose(fp);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NKIND; ii++){ */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NKIND * Ntot; ii++)
    rad[ii] = (PLFLT)SQRT(ball[INDEX2D(NKIND * Ntot, 4, ii, 3)]) * length2astro;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ntot; ii++){
    //---------------------------------------------------------------------
    if( rad[ii] > FLT_EPSILON ){
      //-------------------------------------------------------------------
      const PLFLT inv = 1.0 / rad[INDEX2D(NKIND, Ntot, 0, ii)];
      for(int jj = 1; jj < NKIND; jj++)
	rate[INDEX2D(NKIND - 1, Ntot, jj - 1, ii)] = rad[INDEX2D(NKIND, Ntot, jj, ii)] * inv;
      //-------------------------------------------------------------------
#if 0
      for(int jj = 1; jj < NKIND; jj++)
	if( rate[INDEX2D(NKIND - 1, Ntot, jj - 1, ii)] > 1.5 )
	  fprintf(stderr, "ALERT: %d-th group in \"%s\" has %le times bigger radius (%e) compared to \"%s\" (%e)\n",
		  jj, ballname[jj], rate[INDEX2D(NKIND - 1, Ntot, jj - 1, ii)], rad[INDEX2D(NKIND, Ntot, jj, ii)], ballname[0], rad[INDEX2D(NKIND, Ntot, 0, ii)]);
#endif
      //-------------------------------------------------------------------
    }/* if( rad[ii] > FLT_EPSILON ){ */
    //---------------------------------------------------------------------
    else{
      //-------------------------------------------------------------------
      rad[ii] = 0.0;
      for(int jj = 1; jj < NKIND; jj++)
	rate[INDEX2D(NKIND - 1, Ntot, jj - 1, ii)] = 1.0;
      //-------------------------------------------------------------------
    }/* else{ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ntot; ii++){ */
  //-----------------------------------------------------------------------
  /* confirmation */
  for(int ii = 0; ii < Ntot * (NKIND - 1); ii++)
    if( rate[ii] < (UNITY - EPSILON) )
      fprintf(stderr, "BUG: %d-th group in \"%s\" has %le times smaller radius (%e) compared to \"%s\" (%e)\n", ii / (NKIND - 1), ballname[(ii % (NKIND - 1)) + 1], rate[ii], rad[INDEX2D(NKIND, Ntot, 1 + (ii / Ntot), ii % Ntot)], ballname[0], rad[INDEX2D(NKIND, Ntot, 0, ii % Ntot)]);
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < Ntot * NKIND; ii++)
    printf("%e\n", rad[ii]);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* plot elapsed time as a function of accuracy */
  //-----------------------------------------------------------------------
  double xmax = 0.0;
#ifdef  LOGPLOT_HOR
  double xmin = DBL_MAX;
#endif//LOGPLOT_HOR
  for(int ii = 0; ii < Ntot; ii++){
    xmax = fmax(xmax, rad[ii]);
#ifdef  LOGPLOT_HOR
    if( rad[ii] > FLT_EPSILON )
      xmin = fmin(xmin, rad[ii]);
#endif//LOGPLOT_HOR
  }/* for(int ii = 0; ii < Ntot; ii++){ */
  double ymax = 1.0;
  for(int ii = 0; ii < Ntot * (NKIND - 1); ii++)
    ymax = fmax(ymax, rate[ii]);
#if 0
  printf("Ntot = %d\n", Ntot);
  printf("xmax = %e, ymax = %e\n", xmax, ymax);
#ifdef  LOGPLOT_HOR
  printf("xmin = %e\n", xmin);
#endif//LOGPLOT_HOR
  exit(0);
#endif
  //-----------------------------------------------------------------------
#ifdef  LOGPLOT_HOR
  for(int ii = 0; ii < Ntot; ii++)
    rad[ii] = log10(rad[ii]);
#endif//LOGPLOT_HOR
  //-----------------------------------------------------------------------
#ifdef  LOGPLOT_VER
  for(int ii = 0; ii < Ntot * (NKIND - 1); ii++)
    rate[ii] = log10(rate[ii]);
#endif//LOGPLOT_VER
  //-----------------------------------------------------------------------
  PLplotPltRange range;
#ifdef  LOGPLOT_HOR
  range.xmin = log10(xmin);
  range.xmax = log10(xmax);
  const double Lx = range.xmax - range.xmin;
  range.xmin -= 0.125 * Lx;
  range.xmax += 0.125 * Lx;
  range.xlog = LOGARITHMIC_PLOT;
#else//LOGPLOT_HOR
  range.xmin = 0.0;
  range.xmax = xmax;
  range.xmax *= 1.125;
  range.xlog = LINEAR_PLOT;
#endif//LOGPLOT_HOR
#ifdef  LOGPLOT_VER
  range.ymin = log10(1.0);
  range.ymax = log10(ymax);
  const double Ly = range.ymax - range.ymin;
  range.ymax += 0.125 * Ly;
  range.ylog = LOGARITHMIC_PLOT;
#else///LOGPLOT_VER
  range.ymin = 1.0;
  range.ymax = ymax;
  range.ymax *= 1.125;
  range.ylog = LINEAR_PLOT;
#endif//LOGPLOT_VER
  range.xgrd = range.ygrd = true;
#if 0
  printf("xmin = %e, xmax = %e, ymin = %e, ymax = %e\n", range.xmin, range.xmax, range.ymin, range.ymax);
  exit(0);
#endif
  //-----------------------------------------------------------------------
#if 0
  plotElapsedTime(Ntot, rad, rate, range, file, argc, argv);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* draw color maps */
  //-----------------------------------------------------------------------
  PLFLT *xmap;  allocPLFLT(&xmap,               NX * NY);
  PLFLT *ymap;  allocPLFLT(&ymap,               NX * NY);
  PLFLT *fmap;  allocPLFLT(&fmap, (NKIND - 1) * NX * NY);
  //-----------------------------------------------------------------------
  const PLFLT dx = (range.xmax - range.xmin) / (PLFLT)NX;
  const PLFLT dy = (range.ymax - range.ymin) / (PLFLT)NY;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NX; ii++){
    //---------------------------------------------------------------------
    const PLFLT xx = range.xmin + dx * (0.5 + (PLFLT)ii);
    //---------------------------------------------------------------------
    for(int jj = 0; jj < NY; jj++){
      //-------------------------------------------------------------------
      xmap[INDEX2D(NX, NY, ii, jj)] = xx;
      ymap[INDEX2D(NX, NY, ii, jj)] = range.ymin + dy * (0.5 + (PLFLT)jj);
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < NY; jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NX; ii++){ */
  //-----------------------------------------------------------------------
  const PLFLT frac = 1.0 / (PLFLT)Ntot;
  const PLFLT dxinv = 1.0 / dx;
  const PLFLT dyinv = 1.0 / dy;
  for(int ii = 0; ii < (NKIND - 1) * NX * NY; ii++)
    fmap[ii] = 0.0;
  for(int ii = 0; ii < Ntot; ii++){
    //---------------------------------------------------------------------
    const int ix = (int)floor((rad[ii] - range.xmin) * dxinv);
    //---------------------------------------------------------------------
    if( (ix >= 0) && (ix < NX) )
      for(int jj = 0; jj < NKIND - 1; jj++){
	//-----------------------------------------------------------------
	const int jy = (int)floor((rate[INDEX2D(NKIND - 1, Ntot, jj, ii)] - range.ymin) * dyinv);
	//-----------------------------------------------------------------
	if( (jy >= 0) && (jy < NY) )
	  fmap[INDEX(NKIND - 1, NX, NY, jj, ix, jy)] += frac;
	//-----------------------------------------------------------------
      }/* for(int jj = 0; jj < NKIND - 1; jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ntot; ii++){ */
  //-----------------------------------------------------------------------
  PLFLT fmapMin = 0.0;
  PLFLT fmapMax = -DBL_MAX;
#ifdef  NORMALIZE_COL
  for(int ii = 0; ii < NKIND - 1; ii++){
    //---------------------------------------------------------------------
    fmapMax = -DBL_MAX;
    //---------------------------------------------------------------------
    for(int jj = 0; jj < NX * NY; jj++)
      if( fmap[INDEX2D(NKIND - 1, NX * NY, ii, jj)] > fmapMax )
	fmapMax = fmap[INDEX2D(NKIND - 1, NX * NY, ii, jj)];
    //---------------------------------------------------------------------
    fmapMax = 1.0 / fmapMax;
    for(int jj = 0; jj < NX * NY; jj++)
      fmap[INDEX2D(NKIND - 1, NX * NY, ii, jj)] *= fmapMax;
#ifdef  LOGPLOT_COL
    for(int jj = 0; jj < NX * NY; jj++)
      if( fmap[INDEX2D(NKIND - 1, NX * NY, ii, jj)] > 0.0 ){
	fmap[INDEX2D(NKIND - 1, NX * NY, ii, jj)] = log10(fmap[INDEX2D(NKIND - 1, NX * NY, ii, jj)]);
	if( fmap[INDEX2D(NKIND - 1, NX * NY, ii, jj)] < fmapMin )
	  fmapMin = fmap[INDEX2D(NKIND - 1, NX * NY, ii, jj)];
      }/* if( fmap[INDEX2D(NKIND - 1, NX * NY, ii, jj)] > 0.0 ){ */
      else
	fmap[INDEX2D(NKIND - 1, NX * NY, ii, jj)] = -DBL_MAX;
#endif//LOGPLOT_COL
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NKIND - 1; ii++){ */
  fmapMax = 1.0;
#ifdef  LOGPLOT_COL
  fmapMax = 0.0;
#endif//LOGPLOT_COL
#else///NORMALIZE_COL
  for(int ii = 0; ii < (NKIND - 1) * NX * NY; ii++)
    if( fmap[ii] > 0.0 ){
#ifdef  LOGPLOT_COL
      fmap[ii] = log10(fmap[ii]);
#endif//LOGPLOT_COL
      if( fmap[ii] < fmapMin )	fmapMin = fmap[ii];
      if( fmap[ii] > fmapMax )	fmapMax = fmap[ii];
    }/* if( fmap[ii] > 0.0 ){ */
#endif//NORMALIZE_COL
  //-----------------------------------------------------------------------
#ifdef  LOGPLOT_HOR
  range.xmin += 0.125 * Lx;
  range.xmax -= 0.125 * Lx;
#else//LOGPLOT_HOR
  range.xmax /= 1.125;
#endif//LOGPLOT_HOR
#ifdef  LOGPLOT_VER
  range.ymax -= 0.125 * Ly;
#else///LOGPLOT_VER
  range.ymax /= 1.125;
#endif//LOGPLOT_VER
#if 0
  printf("xmin = %e, xmax = %e, ymin = %e, ymax = %e\n", range.xmin, range.xmax, range.ymin, range.ymax);
  exit(0);
#endif
  //-----------------------------------------------------------------------
#if 1
#ifndef LOGPLOT_VER
  range.ymax = fmin(range.ymax, 1.59);
#endif//LOGPLOT_VER
#endif
  //-----------------------------------------------------------------------
  plotDistribution(xmap, ymap, fmap, fmapMin, fmapMax, range, file, argc, argv);
  //-----------------------------------------------------------------------
  free(xmap);
  free(ymap);
  free(fmap);
  //-----------------------------------------------------------------------





  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(ball);
  free(rad);
  free(rate);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotElapsedTime
(const int Ntot, PLFLT * restrict rad, PLFLT * restrict rate, PLplotPltRange box, char file[], int argc, char **argv)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* global setting(s) */
  //-----------------------------------------------------------------------
  /* maximum number of data kinds */
  const PLINT pkind = NKIND - 2;
  const PLINT lkind = 0;
  //-----------------------------------------------------------------------
  const PLBOOL puni = false;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* preparation for dataset */
  //-----------------------------------------------------------------------
  /* memory allocation for index */
  PLINT kind = (puni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  PLINT *num;  allocPLINT(&num, kind);
  //-----------------------------------------------------------------------
  /* set number of data points */
  for(PLINT ii = 0; ii < kind; ii++)    num[ii] = Ntot;
  //-----------------------------------------------------------------------
  /* set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  for(PLINT ii = 0; ii < pkind; ii++){
    sprintf(pt[ii].type, PLplotSymbolType[smallDot]);
    pt[ii].scale = PLplotSymbolSize[smallDot];
  }/* for(PLINT ii = 0; ii < pkind; ii++){ */
  pt[0].color = RED;
  pt[1].color = BLUE;
  pt[2].color = BLACK;
  //-----------------------------------------------------------------------
  /* set labels */
  char hlab[PLplotCharWords];  sprintf(hlab, "#fir#fr#dSEB#u (%s)", length_astro_unit_name4plot);
  char vlab[PLplotCharWords];  sprintf(vlab, "#fir#fr#dball#u / #fir#fr#dSEB#u");
  //-----------------------------------------------------------------------
  /* set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write = false;
  sprintf(basecap.side, "%s", "b");
  //-----------------------------------------------------------------------
  /* set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = true;
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
  for(int ii = 0; ii < kind; ii++)
    sprintf(pleg.text[ii], synonym[ii]);
  pleg.pos = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;
  pleg.Tcolor = true;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation to exploit multiple-plot function */
  //-----------------------------------------------------------------------
  /* maximum number of panels in each direction */
  const PLINT nxpanel = 1;
  const PLINT nypanel = 1;
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
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* configure to plot elapsed time as a function of tree error by multiple processes */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < kind; ii++){
    hor[ii] = rad;
    ver[ii] = &rate[ii * Ntot];
  }/* for(int ii = 0; ii < kind; ii++){ */
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
      /* line setting(s) */
      nlkind[idx] = lkind;
      line  [idx] = ls;
      lnum  [idx] = num;
      lx    [idx] = &hor[INDEX2D(nxpanel * nypanel, lkind, idx, 0)];
      ly    [idx] = &ver[INDEX2D(nxpanel * nypanel, lkind, idx, 0)];
      //-------------------------------------------------------------------
      /* point setting(s) */
      npkind[idx] = pkind;
      point [idx] = pt;
      pnum  [idx] = num;
      px    [idx] = &hor[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];
      py    [idx] = &ver[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];
      //-------------------------------------------------------------------
      /* errorbar setting(s) */
      errbar[idx] = false;
      //-------------------------------------------------------------------
      /* plot area */
      range[idx] = box;
      //-------------------------------------------------------------------
      /* label setting(s) */
      xlabel[idx] = hlab;
      ylabel[idx] = vlab;
      //-------------------------------------------------------------------
      /* misc(s) */
      leg[idx] = pleg;
      uni[idx] = puni;
      cap[idx] = basecap;
      //-------------------------------------------------------------------
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* individual file names */
  sprintf(figfile, "%s_ball", file);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* create figure(s) */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  free(plegTxt);
  //-----------------------------------------------------------------------
  free(num);  free(hor);  free(ver);
  free(ls);  free(pt);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotDistribution
(PLFLT * restrict xmap, PLFLT * restrict ymap, PLFLT * restrict fmap, const PLFLT fmapMin, const PLFLT fmapMax, PLplotPltRange range, char file[], int argc, char **argv)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* preparation for dataset */
  //-----------------------------------------------------------------------
  /* set color scheme */
  PLplotColorScheme *color;
  setDefaultColorScheme(1, &color, false);
  for(PLINT ii = 0; ii < 1; ii++){
#ifdef  MONOCHROME_PLOT
    color[ii].type = 0;
    color[ii].reverse = true;
#else///MONOCHROME_PLOT
    /* color[ii].type = 11; */
    color[ii].type = 12;
    /* color[ii].type = 22; */
#endif//MONOCHROME_PLOT
  }/* for(PLINT ii = 0; ii < 1; ii++){ */
  //-----------------------------------------------------------------------
  /* set labels */
  char hlab[PLplotCharWords];  sprintf(hlab, "#fir#fr#dSEB#u (%s)", length_astro_unit_name4plot);
  char vlab[PLplotCharWords];  sprintf(vlab, "#fir#fr#dball#u / #fir#fr#dSEB#u");
#ifdef  NORMALIZE_COL
  char flab[PLplotCharWords];  sprintf(flab, "#fiN#fr / #fiN#fr#dmax#u");
#else///NORMALIZE_COL
  char flab[PLplotCharWords];  sprintf(flab, "Fraction");
#endif//NORMALIZE_COL
  //-----------------------------------------------------------------------
  /* set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write = true;
  sprintf(basecap.side, "%s", "t");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation to exploit multiple-plot function */
  //-----------------------------------------------------------------------
  /* maximum number of panels in each direction */
  const PLINT nxpanel = 2;
  const PLINT nypanel = 2;
  //-----------------------------------------------------------------------
  /* specify plot or skip panel */
  PLBOOL *draw;  allocPLBOOL(&draw, nxpanel * nypanel);
  /* arrays related to color map(s) */
  PLINT  *ndat;  allocPLINT        (&ndat, nxpanel * nypanel);
  PLINT  *nx  ;  allocPLINT        (&nx  , nxpanel * nypanel);
  PLINT  *ny  ;  allocPLINT        (&ny  , nxpanel * nypanel);
  PLFLT **xdat;  allocPointer4PLFLT(&xdat, nxpanel * nypanel);
  PLFLT **ydat;  allocPointer4PLFLT(&ydat, nxpanel * nypanel);
  PLFLT **fdat;  allocPointer4PLFLT(&fdat, nxpanel * nypanel);
  PLBOOL *flog;  allocPLBOOL       (&flog, nxpanel * nypanel);
  /* arrays related to colorbar(s) */
  PLplotColorScheme * col;  allocPLplotColorScheme(& col, nxpanel * nypanel);
  PLBOOL            * fix;  allocPLBOOL           (& fix, nxpanel * nypanel);
  PLFLT             *fmin;  allocPLFLT            (&fmin, nxpanel * nypanel);
  PLFLT             *fmax;  allocPLFLT            (&fmax, nxpanel * nypanel);
  /* specify overplot or skip vector map */
  PLBOOL *vmap;  allocPLBOOL(&vmap, nxpanel * nypanel);
  /* arrays related to vector map(s) */
  PLINT  *xskip;  allocPLINT        (&xskip, nxpanel * nypanel);
  PLINT  *yskip;  allocPLINT        (&yskip, nxpanel * nypanel);
  PLFLT **   vx;  allocPointer4PLFLT(&   vx, nxpanel * nypanel);
  PLFLT **   vy;  allocPointer4PLFLT(&   vy, nxpanel * nypanel);
  /* arrays related to vector style(s) */
  PLplotVectorStyle *vs_p;  allocPLplotVectorStyle(&vs_p, nxpanel * nypanel);
  /* arrays for caption(s) and plot range(s) */
  PLplotCaption  *cap;  allocPLplotCaption (&cap, nxpanel * nypanel);
  PLplotPltRange *box;  allocPLplotPltRange(&box, nxpanel * nypanel);
  /* arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  char **flabel;  allocPointer4Char4PLplot(&flabel, nxpanel * nypanel);
  /* array to set figure name */
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel * nypanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel * nypanel);
  assignChar4PLplot(nxpanel * nypanel, figname, _figname);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      //-------------------------------------------------------------------
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      //-------------------------------------------------------------------
      /* global setting(s) */
      draw[idx] = true;
      //-------------------------------------------------------------------
      /* number of arrays */
      ndat[idx] = NX * NY;
      nx  [idx] = NX;
      ny  [idx] = NY;
      //-------------------------------------------------------------------
      /* position of grid points */
      xdat[idx] = xmap;
      ydat[idx] = ymap;
      //-------------------------------------------------------------------
      /* plot settings for color map */
#ifdef  LOGPLOT_COL
      flog[idx] = LOGARITHMIC_PLOT;
#else///LOGPLOT_COL
      flog[idx] =      LINEAR_PLOT;
#endif//LOGPLOT_COL
      col [idx] = color[0];
      fix [idx] = true;
      fmin[idx] = fmapMin;
      fmax[idx] = fmapMax;
      //-------------------------------------------------------------------
      /* vector field */
      vmap[idx] = false;
      //-------------------------------------------------------------------
      /* label setting(s) */
      xlabel[idx] = hlab;
      ylabel[idx] = vlab;
      flabel[idx] = flab;
      //-------------------------------------------------------------------
      /* misc(s) */
      cap[idx] = basecap;
      box[idx] = range;
      //-------------------------------------------------------------------
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */
  //-----------------------------------------------------------------------
  /* individual setting(s) */
  fdat[INDEX2D(nxpanel, nypanel, 0, 1)] = &fmap[NX * NY * 0];  sprintf(cap[INDEX2D(nxpanel, nypanel, 0, 1)].text, "(%c) %s", 97 + INDEX2D(nypanel, nxpanel, 0, 0), synonym[0]);
  fdat[INDEX2D(nxpanel, nypanel, 1, 1)] = &fmap[NX * NY * 1];  sprintf(cap[INDEX2D(nxpanel, nypanel, 1, 1)].text, "(%c) %s", 97 + INDEX2D(nypanel, nxpanel, 0, 1), synonym[1]);
  fdat[INDEX2D(nxpanel, nypanel, 0, 0)] = &fmap[NX * NY * 2];  sprintf(cap[INDEX2D(nxpanel, nypanel, 0, 0)].text, "(%c) %s", 97 + INDEX2D(nypanel, nxpanel, 1, 0), synonym[2]);
  fdat[INDEX2D(nxpanel, nypanel, 1, 0)] = &fmap[NX * NY * 3];  sprintf(cap[INDEX2D(nxpanel, nypanel, 1, 0)].text, "(%c) %s", 97 + INDEX2D(nypanel, nxpanel, 1, 1), synonym[3]);
  //-----------------------------------------------------------------------
  /* individual file names */
  char figfile[PLplotCharWords];  sprintf(figfile, "%s_ball_rate", file);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 1)], "%s_%s", figfile, ballname[1 + 0]);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 1, 1)], "%s_%s", figfile, ballname[1 + 1]);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 0)], "%s_%s", figfile, ballname[1 + 2]);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 1, 0)], "%s_%s", figfile, ballname[1 + 3]);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* create figure(s) */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel - 1; jj++)
      col[INDEX2D(nxpanel, nypanel, ii, jj)].addColBar = false;
  drawMaps(nxpanel, nypanel, draw, true, true,
	   ndat, xdat, ydat, fdat, flog, nx, ny,
	   col, fix, fmin, fmax,
	   vmap, vs_p, xskip, vx, yskip, vy,
	   cap, box, xlabel, ylabel, flabel, "",
	   figfile, argc, argv);
  //-----------------------------------------------------------------------
#if 0
  for(PLINT idx = 0; idx < nxpanel * nypanel; idx++)
    drawMaps(1, 1, &draw[idx], false, false,
	     &ndat[idx], &xdat[idx], &ydat[idx], &fdat[idx], &flog[idx], &nx[idx], &ny[idx],
	     &col[idx], &fix[idx], &fmin[idx], &fmax[idx],
	     &vmap[idx], &vs_p[idx], &xskip[idx], &vx[idx], &yskip[idx], &vy[idx],
	     &cap[idx], &box[idx], &xlabel[idx], &ylabel[idx], &flabel[idx], "",
	     figname[idx], argc, argv);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(draw);
  free(ndat);  free(nx);  free(ny);
  free(xdat);  free(ydat);  free(fdat);
  free(flog);  free(col);  free(fix);  free(fmin);  free(fmax);
  free(vmap);  free(xskip);  free(yskip);  free(vx);  free(vy);  free(vs_p);
  free(cap);  free(box);
  free(xlabel);  free(ylabel);  free(flabel);
  free(figname);  free(_figname);
  //-----------------------------------------------------------------------
  free(color);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
