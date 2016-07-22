/*************************************************************************\
 *                                                                       *
                  last updated on 2016/06/11(Sat) 18:12:05
 *                                                                       *
 *    Plot Code of Breakdown of Tree code on K20X with acceleration MAC  *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#define USE_DISTINCT_PANELS
//-------------------------------------------------------------------------
#define HIGHLIGHTING_COLOR_MAP
#define HIGHLIGHTING_COLOR_MAP_VERTICAL
//-------------------------------------------------------------------------
#   if  !defined(HIGHLIGHTING_COLOR_MAP) && defined(HIGHLIGHTING_COLOR_MAP_VERTICAL)
#        define  HIGHLIGHTING_COLOR_MAP
#endif//!defined(HIGHLIGHTING_COLOR_MAP) && defined(HIGHLIGHTING_COLOR_MAP_VERTICAL)
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <plplot.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <plplotlib.h>
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define NMAX (1024) */
#define NMAX (45)
//-------------------------------------------------------------------------
#define NARC (3)
/* #define NARC (2) */
#define NDST (2)
//-------------------------------------------------------------------------
static char arcname[NARC][16] =
  {"comq", "tcag"
#   if  NARC > 2
   , "dm"
#endif//NARC > 2
};
static char icname[NDST][4] = {"nfw", "m31"};
//-------------------------------------------------------------------------
static char gpuname[NARC][16] =
  {"M2090", "K20X"
#   if  NARC > 2
   , "GTX TITAN X"
#endif//NARC > 2
};
static char dstname[NDST][4] = {"NFW", "M31"};
//-------------------------------------------------------------------------
#define NINTFOLDER "performance"
//-------------------------------------------------------------------------
#ifdef  HIGHLIGHTING_COLOR_MAP_VERTICAL
//-------------------------------------------------------------------------
#define NMAX_SLOW_SEQUENCE (6)
static PLINT slowSeqNum [NDST][NARC] = {{5, 4, 6}, {3, 4, 3}};
static PLFLT slowSeqStep[NDST][NARC][NMAX_SLOW_SEQUENCE] =
  {
    {{11.0, 19.0, 26.0, 34.0, 41.0, DBL_MAX}, {15.0, 23.0, 38.0, 42.0, DBL_MAX, DBL_MAX}, {0.0, 8.0, 15.0, 23.0, 30.0, 38.0}},
    {{0.0, 16.0, 32.0, DBL_MAX, DBL_MAX, DBL_MAX}, {0.0, 16.0, 32.0, 33.0, DBL_MAX, DBL_MAX}, {0.0, 16.0, 32.0, DBL_MAX, DBL_MAX, DBL_MAX}}
  };
//-------------------------------------------------------------------------
#endif//HIGHLIGHTING_COLOR_MAP_VERTICAL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* assumed Flops per interaction */
//-------------------------------------------------------------------------
static const PLFLT FlopsCount[2][NARC] =
  {
    /* based on flops of rsqrtf = (throughput of FP32 add) / (throughput of FP32 reciprocal square root) / (flops of FP32 add) = 8, 6, 4 Flops */
    {28.0, 26.0
#   if  NARC > 2
     , 24.0
#endif//NARC > 2
    },
    /* based on flops of rsqrtf = (throughput of FP32 FMA) / (throughput of FP32 reciprocal square root) / (flops of FP32 FMA) = 4, 3, 2 Flops */
    {24.0, 23.0
#   if  NARC > 2
     , 22.0
#endif//NARC > 2
    }
  };
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotPerformance
(PLFLT step[restrict][NARC][NMAX], PLFLT nips[restrict][NARC][NMAX], PLFLT nint[restrict][NARC][NMAX], PLplotPltRange box[restrict], int argc, char **argv);
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
  //-----------------------------------------------------------------------
  fprintf(stdout, "# IC on host\t# of interactions per sec\tGFlop/s(add based)\tGFlop/s(FMA based)\n");
  fprintf(stdout, "#\taveraged, maximum, minimum\n");
  //-----------------------------------------------------------------------
  /* read the measured data */
  static    int step[NDST][NARC][NMAX];
  static  ulong nint[NDST][NARC][NMAX];/* # of interactions */
  static double time[NDST][NARC][NMAX];/* elapsed time (second) */
  static PLFLT fstep[NDST][NARC][NMAX];
  static PLFLT fnint[NDST][NARC][NMAX];/* # of interactions */
  static PLFLT fnips[NDST][NARC][NMAX];/* # of interactions per second */
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NARC; jj++){
      //-------------------------------------------------------------------
      FILE *fp;
      char filename[128];
      sprintf(filename, "%s/%s/%s.log", NINTFOLDER, arcname[jj], icname[ii]);
      fp = fopen(filename, "r");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
      int ll = 0;
      /* #step	interactions	Ngroups	calcGrav_dev */
      while( (EOF != fscanf(fp, "%d\t%zu\t%*d\t%le", &step[ii][jj][ll], &nint[ii][jj][ll], &time[ii][jj][ll])) && (ll < NMAX) )
	ll++;
      fclose(fp);
      //-------------------------------------------------------------------
      double minnips = DBL_MAX;
      double maxnips = 0.0;
      double avgnips = 0.0;
      double    ttot = 0.0;
      for(int kk = 0; kk < ll; kk++){
	//-----------------------------------------------------------------
	fstep[ii][jj][kk] = (PLFLT)step[ii][jj][kk];
	double nips = (double)nint[ii][jj][kk];	ttot += time[ii][jj][kk];	avgnips += nips;
	fnint[ii][jj][kk] = log10(nips);	nips /= time[ii][jj][kk];
	fnips[ii][jj][kk] = log10(nips);
	//-----------------------------------------------------------------
	if( minnips > nips )	  minnips = nips;
	if( maxnips < nips )	  maxnips = nips;
	//-----------------------------------------------------------------
      }/* for(int kk = 0; kk < ll; kk++){ */
      //-------------------------------------------------------------------
      avgnips /= ttot;
      fprintf(stdout, "%s on %s:\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n", dstname[ii], arcname[jj],
	      avgnips, maxnips, minnips,
	      avgnips * FlopsCount[0][jj] * 1.0e-9, maxnips * FlopsCount[0][jj] * 1.0e-9, minnips * FlopsCount[0][jj] * 1.0e-9,
	      avgnips * FlopsCount[1][jj] * 1.0e-9, maxnips * FlopsCount[1][jj] * 1.0e-9, minnips * FlopsCount[1][jj] * 1.0e-9
	      );
      //-------------------------------------------------------------------
      for(int kk = ll; kk < NMAX; kk++){
	step[ii][jj][kk]  =  step[ii][jj][ll - 1];
	fstep[ii][jj][kk] = fstep[ii][jj][ll - 1];
	/* fnint[ii][jj][kk] = fnint[ii][jj][ll - 1]; */
	/* fnips[ii][jj][kk] = fnips[ii][jj][ll - 1]; */
	fnint[ii][jj][kk] = DBL_MAX;
	fnips[ii][jj][kk] = DBL_MAX;
      }/* for(int mm = ll; mm < NMAX; mm++){ */
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < NARC; jj++){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* plot # of interactions per second as a function of step */
  //-----------------------------------------------------------------------
  PLplotPltRange range[NDST];
  for(int ii = 0; ii < NDST; ii++){
    //---------------------------------------------------------------------
    range[ii].xmin =         0.0;
    range[ii].xlog = LINEAR_PLOT;
    range[ii].ymin = log10(5.0e+7);
    range[ii].ymax = log10(8.0e+11);
    range[ii].ylog = LOGARITHMIC_PLOT;
    range[ii].xgrd = range[ii].ygrd = true;
    //---------------------------------------------------------------------
    int smax = step[ii][0][NMAX - 1];
    for(int jj = 1; jj < NARC; jj++)
      if( smax < step[ii][jj][NMAX - 1] )
	smax = step[ii][jj][NMAX - 1];
    range[ii].xmax = (PLFLT)smax;
    //---------------------------------------------------------------------
#if 0
    range[ii].xmax = fmin(range[ii].xmax, 45.0);
#endif
    //---------------------------------------------------------------------
#if 1
    const PLFLT length = 2.0e-2 * (range[ii].xmax - range[ii].xmin);
    range[ii].xmin -= length;
    range[ii].xmax += length;
#endif
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NDST; ii++){ */
  //-----------------------------------------------------------------------
  plotPerformance(fstep, fnips, fnint, range, argc, argv);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotPerformance
(PLFLT step[restrict][NARC][NMAX], PLFLT nips[restrict][NARC][NMAX], PLFLT nint[restrict][NARC][NMAX], PLplotPltRange box[restrict], int argc, char **argv)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* global setting(s) */
  //-----------------------------------------------------------------------
  /* maximum number of data kinds */
  const PLINT pkind = 2 * NARC;
  const PLINT lkind = 0;
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
  for(PLINT ii = 0; ii < kind; ii++)    num[ii] = NMAX;
  //-----------------------------------------------------------------------
  /* set symbol and line style */
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
#ifdef  USE_DISTINCT_PANELS
#if 1
  pt[       0].color =   BLACK;  sprintf(pt[       0].type, PLplotSymbolType[openCircle ]);  pt[       0].scale = PLplotSymbolSize[openCircle ] * 0.6;
  pt[NARC    ].color =   BLACK;  sprintf(pt[NARC    ].type, PLplotSymbolType[openCircle ]);  pt[NARC    ].scale = PLplotSymbolSize[openCircle ] * 0.6;
  pt[       1].color =     RED;  sprintf(pt[       1].type, PLplotSymbolType[openDiamond]);  pt[       1].scale = PLplotSymbolSize[openDiamond] * 0.6;
  pt[NARC + 1].color =     RED;  sprintf(pt[NARC + 1].type, PLplotSymbolType[openDiamond]);  pt[NARC + 1].scale = PLplotSymbolSize[openDiamond] * 0.6;
#   if  NARC > 2
  pt[       2].color =    BLUE;  sprintf(pt[       2].type, PLplotSymbolType[openSquare ]);  pt[       2].scale = PLplotSymbolSize[openSquare ] * 0.6;
  pt[NARC + 2].color =    BLUE;  sprintf(pt[NARC + 2].type, PLplotSymbolType[openSquare ]);  pt[NARC + 2].scale = PLplotSymbolSize[openSquare ] * 0.6;
#endif//NARC > 2
#else
  pt[       0].color =   BLACK;  sprintf(pt[       0].type, PLplotSymbolType[fillCircle ]);  pt[       0].scale = PLplotSymbolSize[fillCircle ] * 0.6;
  pt[NARC    ].color =   BLACK;  sprintf(pt[NARC    ].type, PLplotSymbolType[openCircle ]);  pt[NARC    ].scale = PLplotSymbolSize[openCircle ] * 0.6;
  pt[       1].color =     RED;  sprintf(pt[       1].type, PLplotSymbolType[fillDiamond]);  pt[       1].scale = PLplotSymbolSize[fillDiamond] * 0.6;
  pt[NARC + 1].color =     RED;  sprintf(pt[NARC + 1].type, PLplotSymbolType[openDiamond]);  pt[NARC + 1].scale = PLplotSymbolSize[openDiamond] * 0.6;
#   if  NARC > 2
  pt[       2].color =    BLUE;  sprintf(pt[       2].type, PLplotSymbolType[fillSquare ]);  pt[       2].scale = PLplotSymbolSize[fillSquare ] * 0.6;
  pt[NARC + 2].color =    BLUE;  sprintf(pt[NARC + 2].type, PLplotSymbolType[openSquare ]);  pt[NARC + 2].scale = PLplotSymbolSize[openSquare ] * 0.6;
#endif//NARC > 2
#endif
#else///USE_DISTINCT_PANELS
  pt[0].color =   BLACK;  sprintf(pt[0].type, PLplotSymbolType[openCircle ]);  pt[0].scale = PLplotSymbolSize[openCircle ] * 0.6;
  pt[1].color =   BLACK;  sprintf(pt[1].type, PLplotSymbolType[fillCircle ]);  pt[1].scale = PLplotSymbolSize[fillCircle ] * 0.6;
  pt[2].color =     RED;  sprintf(pt[2].type, PLplotSymbolType[openDiamond]);  pt[2].scale = PLplotSymbolSize[openDiamond] * 0.6;
  pt[3].color =     RED;  sprintf(pt[3].type, PLplotSymbolType[fillDiamond]);  pt[3].scale = PLplotSymbolSize[fillDiamond] * 0.6;
#   if  NARC > 2
  pt[4].color =    BLUE;  sprintf(pt[4].type, PLplotSymbolType[openSquare ]);  pt[4].scale = PLplotSymbolSize[openSquare ] * 0.6;
  pt[5].color =    BLUE;  sprintf(pt[5].type, PLplotSymbolType[fillSquare ]);  pt[5].scale = PLplotSymbolSize[fillSquare ] * 0.6;
#endif//NARC > 2
#endif//USE_DISTINCT_PANELS
  //-----------------------------------------------------------------------
  /* set labels */
  char hlab[PLplotCharWords];  sprintf(hlab, "Step");
  char llab[PLplotCharWords];  sprintf(llab, "## of interactions per second");
  char rlab[PLplotCharWords];  sprintf(rlab, "## of interactions");
  //-----------------------------------------------------------------------
  /* set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write = true;
  sprintf(basecap.side, "%s", "t");
  static PLplotCaption dstcap[NDST];
  for(int ii = 0; ii < NDST; ii++){
    setDefaultCaption(&dstcap[ii]);
    dstcap[ii].write = true;
    sprintf(dstcap[ii].side, "%s", "t");
    sprintf(dstcap[ii].text, "%s", dstname[ii]);
  }/* for(int ii = 0; ii < NDST; ii++){ */
#ifdef  USE_DISTINCT_PANELS
  static PLplotCaption tgtcap[2];
  for(int ii = 0; ii < 2; ii++){
    setDefaultCaption(&tgtcap[ii]);
    tgtcap[ii].write = true;
    sprintf(tgtcap[ii].side, "%s", "t");
  }/* for(int ii = 0; ii < 2; ii++){ */
  sprintf(tgtcap[1].text, "%s", "Calculation amount");
  sprintf(tgtcap[0].text, "%s", "Calculation speed");
#endif//USE_DISTINCT_PANELS
  //-----------------------------------------------------------------------
  /* set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = false;
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
#ifdef  USE_DISTINCT_PANELS
  for(int ii = 0; ii < 2; ii++)
    for(int jj = 0; jj < NARC; jj++)
      sprintf(pleg.text[INDEX2D(2, NARC, ii, jj)], "%s", gpuname[jj]);
  pleg.pos = PL_POSITION_LEFT | PL_POSITION_BOTTOM;
#else///USE_DISTINCT_PANELS
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < 2; jj++)
      sprintf(pleg.text[INDEX2D(NARC, 2, ii, jj)], "%s on %s", (jj == 0) ? "nint" : "nips", gpuname[ii]);
  pleg.pos = PL_POSITION_RIGHT | PL_POSITION_BOTTOM | PL_POSITION_OUTSIDE;
#endif//USE_DISTINCT_PANELS
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation to exploit multiple-plot function */
  //-----------------------------------------------------------------------
  /* maximum number of panels in each direction */
  const PLINT nxpanel = NDST;
#ifdef  USE_DISTINCT_PANELS
  const PLINT nypanel = 2;
#else///USE_DISTINCT_PANELS
  const PLINT nypanel = 1;
#endif//USE_DISTINCT_PANELS
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
  /* configure to plot performance as a function of step */
  //-----------------------------------------------------------------------
  /* data preparation */
  //-----------------------------------------------------------------------
#ifdef  USE_DISTINCT_PANELS
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NARC; jj++){
      hor[INDEX(NDST, 2, NARC, ii, 0, jj)] = &step[ii][jj][0];
      hor[INDEX(NDST, 2, NARC, ii, 1, jj)] = &step[ii][jj][0];
      ver[INDEX(NDST, 2, NARC, ii, 0, jj)] = &nips[ii][jj][0];
      ver[INDEX(NDST, 2, NARC, ii, 1, jj)] = &nint[ii][jj][0];
    }/* for(int jj = 0; jj < NARC; jj++){ */
#else///USE_DISTINCT_PANELS
  for(int ii = 0; ii < NDST; ii++)
    for(int jj = 0; jj < NARC; jj++){
      hor[INDEX(NDST, NARC, 2, ii, jj, 0)] = &step[ii][jj][0];
      hor[INDEX(NDST, NARC, 2, ii, jj, 1)] = &step[ii][jj][0];
      ver[INDEX(NDST, NARC, 2, ii, jj, 0)] = &nint[ii][jj][0];
      ver[INDEX(NDST, NARC, 2, ii, jj, 1)] = &nips[ii][jj][0];
    }/* for(int jj = 0; jj < NARC; jj++){ */
#endif//USE_DISTINCT_PANELS
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
      line  [idx] = NULL;
      lnum  [idx] = NULL;
      lx    [idx] = NULL;
      ly    [idx] = NULL;
      //-------------------------------------------------------------------
      /* point setting(s) */
#ifdef  USE_DISTINCT_PANELS
      npkind[idx] = NARC;
      point [idx] = &pt[NARC * jj];
      pnum  [idx] = num;
      px    [idx] = &hor[INDEX2D(nxpanel * nypanel, NARC, idx, 0)];
      py    [idx] = &ver[INDEX2D(nxpanel * nypanel, NARC, idx, 0)];
#else///USE_DISTINCT_PANELS
      npkind[idx] = pkind;
      point [idx] = pt;
      pnum  [idx] = num;
      px    [idx] = &hor[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];
      py    [idx] = &ver[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];
#endif//USE_DISTINCT_PANELS
      //-------------------------------------------------------------------
      /* errorbar setting(s) */
      errbar[idx] = false;
      //-------------------------------------------------------------------
      /* plot area */
      range[idx] = box[ii];
#ifdef  USE_DISTINCT_PANELS
      if( jj == 0 ){	range[idx].ymax = log10(4.0e+11);	range[idx].ymin = log10(2.0e+9);      }
      if( jj == 1 ){	range[idx].ymax = log10(1.0e+11);	range[idx].ymin = log10(2.0e+7);      }
#endif//USE_DISTINCT_PANELS
      //-------------------------------------------------------------------
      /* label setting(s) */
      xlabel[idx] = hlab;
#ifdef  USE_DISTINCT_PANELS
      ylabel[idx] = (jj == 0) ? llab : rlab;
#else///USE_DISTINCT_PANELS
      ylabel[idx] = llab;
#endif//USE_DISTINCT_PANELS
      //-------------------------------------------------------------------
      /* misc(s) */
      leg[idx] = pleg;
      uni[idx] = puni;
      cap[idx] = basecap;
      /* 97 is "a" in ASCII code */
#ifdef  USE_DISTINCT_PANELS
      sprintf(cap[idx].text, "(%c) %s (%s)", 97 + INDEX2D(nypanel, nxpanel, nypanel - jj - 1, ii), tgtcap[jj].text, dstcap[ii].text);
#else///USE_DISTINCT_PANELS
      sprintf(cap[idx].text, "(%c) %s", 97 + INDEX2D(nypanel, nxpanel, nypanel - jj - 1, ii), dstcap[ii].text);
#endif//USE_DISTINCT_PANELS
      //-------------------------------------------------------------------
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */
  //-----------------------------------------------------------------------
#ifdef  USE_DISTINCT_PANELS
  leg[INDEX2D(nxpanel, nypanel, 0, 1)].write = true;
#endif//USE_DISTINCT_PANELS
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* individual file names */
  sprintf(figfile, "%s", "performance");
  for(int ii = 0; ii < NDST; ii++)
    sprintf(figname[ii], "%s_%s", icname[ii], "performance");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  HIGHLIGHTING_COLOR_MAP
  //-----------------------------------------------------------------------
  /* highlight fast, inetermediate, and slow sequences by background color */
  /* PLFLT colDat[2][2]; */
  PLFLT *_colDat = malloc(sizeof(PLFLT  ) * 2 * 2);
  PLFLT **colDat = malloc(sizeof(PLFLT *) * 2);
  for(int ii = 0; ii < 2; ii++)
    colDat[ii] = _colDat + ii * 2;
  for(int ii = 0; ii < 2; ii++)
    for(int jj = 0; jj < 2; jj++)
      colDat[ii][jj] = 0.5;
#ifndef HIGHLIGHTING_COLOR_MAP_VERTICAL
  PLFLT colMap[NDST][NARC][2][2];
  for(int ii = 0; ii < NDST; ii++){
    for(int jj = 0; jj < NARC; jj++){
      colMap[ii][jj][0][0] = box[ii].xmin;
      colMap[ii][jj][1][0] = box[ii].xmax;
    }/* for(int jj = 0; jj < NARC; jj++){ */
    colMap[ii][0][0][1] = log10(3.6e+9 );    colMap[ii][0][1][1] = log10(5.7e+9 );
    colMap[ii][1][0][1] = log10(3.6e+9 );    colMap[ii][1][1][1] = log10(6.2e+9 );
    colMap[ii][2][0][1] = log10(2.3e+10);    colMap[ii][2][1][1] = log10(3.4e+10);
  }/* for(int ii = 0; ii < NDST; ii++){ */
#endif//HIGHLIGHTING_COLOR_MAP_VERTICAL
  PLFLT pp[NARC][2], c0[NARC][2], c1[NARC][2], c2[NARC][2], aa[NARC][2];
  /* M2090: black (H = 0, L = 0, S = 1) */
  pp[0][0] = 0.0;  c0[0][0] =   0.0;  c1[0][0] = 0.0;  c2[0][0] = 1.0;  aa[0][0] = 0.1;
  pp[0][1] = 1.0;  c0[0][1] =   0.0;  c1[0][1] = 0.0;  c2[0][1] = 1.0;  aa[0][1] = 0.1;
  /* K20X: red (H = 0, L = 1/2, S = 1) */
  pp[1][0] = 0.0;  c0[1][0] =   0.0;  c1[1][0] = 0.5;  c2[1][0] = 1.0;  aa[1][0] = 0.1;
  pp[1][1] = 1.0;  c0[1][1] =   0.0;  c1[1][1] = 0.5;  c2[1][1] = 1.0;  aa[1][1] = 0.1;
  /* TITAN X: blue (H = 240, L = 1/2, S = 1) */
  pp[2][0] = 0.0;  c0[2][0] = 240.0;  c1[2][0] = 0.5;  c2[2][0] = 1.0;  aa[2][0] = 0.1;
  pp[2][1] = 1.0;  c0[2][1] = 240.0;  c1[2][1] = 0.5;  c2[2][1] = 1.0;  aa[2][1] = 0.1;
  PLBOOL rr[2] = {false, false};
  //-----------------------------------------------------------------------
#endif//HIGHLIGHTING_COLOR_MAP
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* create figure(s) */
  //-----------------------------------------------------------------------
  PLBOOL xjoin, yjoin;
  PLFLT scale;
  const PLFLT aspect = 1.0;
  //-----------------------------------------------------------------------
  xjoin = yjoin = true;
  scale = setScaleLength(nxpanel, nypanel);
  initPLplot(argc, argv, figfile);
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      //-------------------------------------------------------------------
      PLINT idx = nypanel * ii + jj;
      //-------------------------------------------------------------------
      /* unified control of lines and points */
      const PLBOOL useLS = false;
      const PLBOOL usePT =  true;
      //-------------------------------------------------------------------
      setViewPort((idx == 0) ? true : false, nxpanel, nypanel, idx, xjoin, yjoin, range[idx], aspect);
      //-------------------------------------------------------------------
#ifdef  HIGHLIGHTING_COLOR_MAP
#ifndef HIGHLIGHTING_COLOR_MAP_VERTICAL
      if( jj == 0 )
#endif//HIGHLIGHTING_COLOR_MAP_VERTICAL
	for(int ss = 0; ss < NARC; ss++){
	  plscmap1la(0, 2, pp[ss], c0[ss], c1[ss], c2[ss], aa[ss], rr);
#ifdef  HIGHLIGHTING_COLOR_MAP_VERTICAL
	  for(int kk = 0; kk < slowSeqNum[ii][ss]; kk++)
	    plimagefr((const PLFLT * const *)colDat, 2, 2, slowSeqStep[ii][ss][kk] - 0.5, slowSeqStep[ii][ss][kk] + 0.5, range[idx].ymin, range[idx].ymax, 0.0, 1.0, 0.0, 1.0, NULL, NULL);
#else///HIGHLIGHTING_COLOR_MAP_VERTICAL
	  plimagefr((const PLFLT * const *)colDat, 2, 2, colMap[ii][ss][0][0], colMap[ii][ss][1][0], colMap[ii][ss][0][1], colMap[ii][ss][1][1], 0.0, 1.0, 0.0, 1.0, NULL, NULL);
#endif//HIGHLIGHTING_COLOR_MAP_VERTICAL
	}/* for(int ss = 0; ss < NARC; ss++){ */
#endif//HIGHLIGHTING_COLOR_MAP
      //-------------------------------------------------------------------
      dotPoint(npkind[idx], pnum[idx], px[idx], py[idx], point[idx], scale);
      //-------------------------------------------------------------------
      /* make frame */
#ifdef  USE_DISTINCT_PANELS
      makeFrame(nxpanel, ii, range[idx].xlog, range[idx].xgrd, xjoin, xlabel[idx],
		nypanel, jj, range[idx].ylog, range[idx].ygrd, yjoin, ylabel[idx], "");
#else///USE_DISTINCT_PANELS
      applyDefaultSettings();
      plschr(0.0, scale);
      if( nxpanel == 1 )	yjoin = false;
      if( nypanel == 1 )	xjoin = false;
      //-------------------------------------------------------------------
      if( (xjoin != true) || (jj == 0) ){
	if( range[idx].xlog == LOGARITHMIC_PLOT )	  plbox("bclnst", 0.0, 0, "", 0.0, 0);
	else	                                          plbox("bcnst" , 0.0, 0, "", 0.0, 0);
	pllab(xlabel[idx], "", "");
      }/* if( (xjoin != true) || (jj == 0) ){ */
      else{
	if( range[idx].xlog == LOGARITHMIC_PLOT )	  plbox("bclst", 0.0, 0, "", 0.0, 0);
	else	                                          plbox("bcst" , 0.0, 0, "", 0.0, 0);
      }/* else{ */
      //-------------------------------------------------------------------
      if( (yjoin != true) || (ii == 0) ){
	/* left edge: write left axis with numbers and label */
	if( range[idx].ylog == LOGARITHMIC_PLOT )	  plbox("", 0.0, 0, "blnstv", 0.0, 0);
	else	                                          plbox("", 0.0, 0, "bnstv" , 0.0, 0);
	pllab("", ylabel[idx], "");
      }/* if( (yjoin != true) || (ii == 0) ){ */
      if( (yjoin != true) || (ii == (nxpanel - 1)) ){
	/* right edge: write right axis with numbers and label */
	if( range[idx].ylog == LOGARITHMIC_PLOT )	  plbox("", 0.0, 0, "clmstv", 0.0, 0);
	else	                                          plbox("", 0.0, 0, "cmstv" , 0.0, 0);
	plmtex("r", 4.0, 0.5, 0.5, rlab);
      }/* if( (yjoin != true) || (ii == (nxpanel - 1)) ){ */
      if( (yjoin == true) && (ii != 0) && (ii != (nxpanel - 1)) ){
	/* internal panels: write left and right axis */
	if( range[idx].ylog == LOGARITHMIC_PLOT ){	  plbox("", 0.0, 0, "blstv", 0.0, 0);	  plbox("", 0.0, 0, "clstv", 0.0, 0);	}
	else{	                                          plbox("", 0.0, 0, "bstv" , 0.0, 0);	  plbox("", 0.0, 0, "cstv" , 0.0, 0);	}
      }/* if( (yjoin == true) && (ii != 0) && (ii != (nxpanel - 1)) ){ */
      if( (yjoin == true) && (ii == 0) ){
	/* write right axis */
	if( range[idx].ylog == LOGARITHMIC_PLOT )	  plbox("", 0.0, 0, "clstv", 0.0, 0);
	else	                                          plbox("", 0.0, 0, "cstv" , 0.0, 0);
      }/* if( (yjoin == true) && (ii == 0) ){ */
      if( (yjoin == true) && (ii == (nxpanel - 1)) ){
	if( range[idx].ylog == LOGARITHMIC_PLOT )	  plbox("", 0.0, 0, "blstv", 0.0, 0);
	else	                                          plbox("", 0.0, 0, "bstv" , 0.0, 0);
      }/* if( (yjoin == true) && (ii == (nxpanel - 1)) ){ */
      //-------------------------------------------------------------------
      plwidth(THIN_LINE * scale * scale);
      setLineStyle(DOTTED_LINE, scale * scale);
      if( range[ii].xgrd == true )	plbox("g", 0.0, 0, " ", 0.0, 0);
      if( range[ii].ygrd == true )	plbox(" ", 0.0, 0, "g", 0.0, 0);
      plwidth(THIN_LINE);
      pllsty(1);
      //-------------------------------------------------------------------
      applyDefaultSettings();
#endif//USE_DISTINCT_PANELS
      //-------------------------------------------------------------------
      if( cap[idx].write == true )	makeCaption(cap[idx], scale);
      if( leg[idx].write == true )
	makeLegend (idx, aspect, scale, nxpanel, xjoin, nypanel, yjoin, uni[idx], leg[idx], useLS, nlkind[idx], line[idx], usePT, npkind[idx], point[idx]);
      //-------------------------------------------------------------------
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */
  exitPLplot();
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NDST; ii++){
    //---------------------------------------------------------------------
    cap[ii].write = false;
    //---------------------------------------------------------------------
    xjoin = false;
    yjoin = false;
    scale = setScaleLength(1, 1);
    //---------------------------------------------------------------------
    initPLplot(argc, argv, figname[ii]);
    //---------------------------------------------------------------------
    /* unified control of lines and points */
    const PLBOOL useLS = false;
    const PLBOOL usePT =  true;
    //---------------------------------------------------------------------
    setViewPort(true, 1, 1, 0, xjoin, yjoin, range[ii], aspect);
    //---------------------------------------------------------------------
    dotPoint(npkind[ii], pnum[ii], px[ii], py[ii], point[ii], scale);
    //---------------------------------------------------------------------
    /* make frame */
    applyDefaultSettings();
    plschr(0.0, scale);
    plbox("bcnst" , 0.0, 0, "", 0.0, 0);    pllab(xlabel[ii], "", "");
    plbox("", 0.0, 0, "blnstv", 0.0, 0);    pllab("", ylabel[ii], "");/*  left axis with numbers and label */
    plbox("", 0.0, 0, "clmstv", 0.0, 0);    plmtex("r", 4.0, 0.5, 0.5, rlab);/* right axis with numbers and label */
    plwidth(THIN_LINE * scale * scale);    setLineStyle(DOTTED_LINE, scale * scale);
    if( range[ii].xgrd == true )      plbox("g", 0.0, 0, " ", 0.0, 0);
    if( range[ii].xgrd == true )      plbox(" ", 0.0, 0, "g", 0.0, 0);
    plwidth(THIN_LINE);    pllsty(1);
    //---------------------------------------------------------------------
    if( cap[ii].write == true )      makeCaption(cap[ii], scale);
    if( leg[ii].write == true )
      makeLegend (ii, aspect, scale, 1, xjoin, 1, yjoin, uni[ii], leg[ii], useLS, nlkind[ii], line[ii], usePT, npkind[ii], point[ii]);
    //---------------------------------------------------------------------
    exitPLplot();
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NDST; ii++){ */
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
  free(plegTxt);
  //-----------------------------------------------------------------------
  free(num);  free(hor);  free(ver);
  free(pt);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
