/**
 * @file plot.ndep.c
 *
 * @brief Plot code of elapsed time of tree code
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2017/10/26 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "macro.h"
#include "plplotlib.h"


#define LINE_ONLY_PLOT
#define LOGPLOT_VER


#define NNUM (15)
#define NMAX (38)

#define NSUB (6)

#define NARC (3)
/* #define NARC (2) */
#define NDST (2)

static char arcname[NARC][16] =
  {"comq", "tcag"
#   if  NARC > 2
   , "dm"
#endif//NARC > 2
};
static char icname[NDST][4] = {"nfw", "m31"};

static char gpuname[NARC][16] =
  {"M2090", "K20X"
#   if  NARC > 2
   , "GTX TITAN X"
#endif//NARC > 2
};
static char dstname[NDST][4] = {"NFW", "M31"};

#define NTOTFOLDER "dependence/ntot"
#define   NIFOLDER "dependence/ni"
#define FRACFOLDER "dependence/breakdown"
/* #define FRACFOLDER "dependence/breakdown_2-06" */


void plotNtotDependence(PLINT ntot[][NDST][NNUM], PLFLT time[][NDST][NNUM], PLplotPltRange box, char file[], const bool averaged, int argc, char **argv);
void plotNtotBreakdown
(PLINT nbrk[restrict][NDST][NNUM], PLFLT tsum[restrict][NDST][NNUM], PLFLT walk[restrict][NDST][NNUM], PLFLT node[restrict][NDST][NNUM], PLFLT make[restrict][NDST][NNUM], PLFLT time[restrict][NDST][NNUM], PLplotPltRange box, int argc, char **argv);
void plotNiDependence(PLFLT inum[restrict][NDST][NSUB][NMAX], PLFLT time[restrict][NDST][NSUB][NMAX], PLplotPltRange box, int argc, char **argv);


int main(int argc, char **argv)
{
  /** initialization */
  if( argc < 1 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 1);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 1 ){ */
  modifyArgcArgv4PLplot(&argc, argv, 1);


  /** read the elapsed time (Ni dependence) */
  static PLINT NI[NARC][NDST][NSUB][NMAX];
  static PLFLT ni[NARC][NDST][NSUB][NMAX], ti[NARC][NDST][NSUB][NMAX];
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < NDST; jj++)
      for(int kk = 0; kk < NSUB; kk++)
	for(int ll = 0; ll < NMAX; ll++)
	  ti[ii][jj][kk][ll] = DBL_MAX;
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < NDST; jj++)
      for(int kk = 0; kk < NSUB; kk++){
	FILE *fp;
	char filename[128];
	sprintf(filename, "%s/%s/%s.%.2d.log", NIFOLDER, arcname[ii], icname[jj], 1 << kk);
	fp = fopen(filename, "r");
	if( fp != NULL ){
	  int ll = 0;
	  while( EOF != fscanf(fp, "%*d\t%d\t%le", &NI[ii][jj][kk][ll], &ti[ii][jj][kk][ll]) )
	    ll++;

	  fclose(fp);

	  for(int mm = 0; mm < ll; mm++){
	    ni[ii][jj][kk][mm] = log10(NI[ii][jj][kk][mm]);
#ifdef  LOGPLOT_VER
	    ti[ii][jj][kk][mm] = log10(ti[ii][jj][kk][mm]);
#endif//LOGPLOT_VER
	  }/* for(int mm = 0; mm < ll; mm++){ */

	  for(int mm = ll; mm < NMAX; mm++){
	    ni[ii][jj][kk][mm] = ni[ii][jj][kk][ll - 1];
/* #ifdef  LOGPLOT_VER */
/* 	    ti[ii][jj][kk][mm] = ti[ii][jj][kk][ll - 1]; */
/* #endif//LOGPLOT_VER */
	  }/* for(int mm = ll; mm < NMAX; mm++){ */
	}/* if( fp != NULL ){ */
      }/* for(int kk = 0; kk < NSUB; kk++){ */

  /** read the breakdown (Ntot dependence) */
  static PLINT nbrk[NARC][NDST][NNUM];
  static PLFLT walk[NARC][NDST][NNUM], make[NARC][NDST][NNUM], node[NARC][NDST][NNUM], time[NARC][NDST][NNUM], tsum[NARC][NDST][NNUM];
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < NDST; jj++){
      for(int kk = 0; kk < NNUM; kk++){
	nbrk[ii][jj][kk] = 0;
	walk[ii][jj][kk] = make[ii][jj][kk] = node[ii][jj][kk] = time[ii][jj][kk] = tsum[ii][jj][kk] = DBL_MAX;
      }/* for(int kk = 0; kk < NNUM; kk++){ */

      FILE *fp;
      char filename[128];
      sprintf(filename, "%s/%s.%s.log", FRACFOLDER, icname[jj], arcname[ii]);
      fp = fopen(filename, "r");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
      int ll = 0;
      PLFLT makeTree, sortPHcurve, pred, corr, exam, search;
      while( EOF != fscanf(fp, "%d\t%le\t%le\t%le\t%*e\t%le\t%*e\t%*e\t%*e\t%*e\t%le\t%le\t%*e\t%*e\t%le\t%le",
			   &nbrk[ii][jj][ll], &walk[ii][jj][ll], &node[ii][jj][ll], &makeTree, &sortPHcurve, &pred, &corr, &exam, &search) ){
	make[ii][jj][ll] = makeTree + sortPHcurve + exam + search;
	time[ii][jj][ll] = pred + corr;
	tsum[ii][jj][ll] = walk[ii][jj][ll] + node[ii][jj][ll] + make[ii][jj][ll] + time[ii][jj][ll];
	ll++;
      }
      fclose(fp);

      for(int kk = 0; kk < ll; kk++)
	if( walk[ii][jj][kk] < DBL_MAX ){
	  tsum[ii][jj][kk] = log10(tsum[ii][jj][kk]);
	  walk[ii][jj][kk] = log10(walk[ii][jj][kk]);
	  make[ii][jj][kk] = log10(make[ii][jj][kk]);
	  node[ii][jj][kk] = log10(node[ii][jj][kk]);
	  time[ii][jj][kk] = log10(time[ii][jj][kk]);
	}/* if( walk[ii][jj][kk] < DBL_MAX ){ */

      for(int kk = ll; kk < NNUM; kk++)
	nbrk[ii][jj][kk] = nbrk[ii][jj][ll - 1];
    }/* for(int jj = 0; jj < NDST; jj++){ */


  /** plot elapsed time as a function of total number of N-body particles */
  PLplotPltRange range;
  /* range.xmin = log10((double)(ntot[NARC - 1][NDST - 1][       0] >> 1)); */
  /* range.xmax = log10((double)(ntot[NARC - 1][NDST - 1][NNUM - 1] << 1)); */
  range.xmin = log10((double)(nbrk[NARC - 1][NDST - 1][       0] >> 1));
  range.xmax = log10((double)(nbrk[NARC - 1][NDST - 1][NNUM - 1] << 1));
  range.xlog = LOGARITHMIC_PLOT;
#ifdef  LOGPLOT_VER
  /* range.ymin = log10(1.0e-3); */
  /* range.ymax = log10(1.0e+0); */
  range.ylog = LOGARITHMIC_PLOT;
#else///LOGPLOT_VER
  /* range.ymin = 0.0; */
  /* range.ymax = 1.0; */
  range.ylog = LINEAR_PLOT;
#endif//LOGPLOT_VER
  range.xgrd = range.ygrd = true;

  /* plotNtotDependence(ntot, mean, range, "mean",  true, argc, argv); */

#ifdef  LOGPLOT_VER
  range.ymin = log10(1.0e-6);
  range.ymax = log10(9.9e-1);
#else///LOGPLOT_VER
  range.ymin = 0.0;
  range.ymax = 1.0;
#endif//LOGPLOT_VER
  plotNtotBreakdown(nbrk, tsum, walk, node, make, time, range, argc, argv);


  /** plot elapsed time as a function of number of i-particles under the fixed total number of N-body particles */
  range.xmin = log10(1.1e+1);
  range.xmax = log10(2.0e+7);
#ifdef  LOGPLOT_VER
  range.ymin = log10(1.0e-5);
  range.ymax = log10(1.0e+1);
  range.ylog = LOGARITHMIC_PLOT;
#else///LOGPLOT_VER
  range.ymin =  0.0;
  range.ymax = 10.0;
  range.ylog = LINEAR_PLOT;
#endif//LOGPLOT_VER

  plotNiDependence(ni, ti, range, argc, argv);


  return (0);
}


void plotNtotDependence(PLINT ntot[][NDST][NNUM], PLFLT time[][NDST][NNUM], PLplotPltRange box, char file[], const bool averaged, int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /** global setting(s) */
  /** maximum number of data kinds */
  const PLINT pkind = NARC;
  const PLINT lkind = NARC;

  const PLBOOL puni = true;


  /** preparation for dataset */
  /** memory allocation for index */
  PLINT kind = (puni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  if( (pkind != kind) || (lkind != kind) ){
    __KILL__(stderr, "ERROR: pkind(%d), lkind(%d), and kind(%d) must be the same\n", pkind, lkind, kind);
  }/* if( (pkind != kind) || (lkind != kind) ){ */
  PLINT *num;  allocPLINT(&num, kind);

  /** set number of data points */
  for(PLINT ii = 0; ii < kind; ii++)    num[ii] = NNUM;

  /** data preparation */
  static PLFLT nbody[NARC][NDST][NNUM];
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < NDST; jj++)
      for(int kk = 0; kk < NNUM; kk++)
	nbody[ii][jj][kk] = log10((double)ntot[ii][jj][kk]);

  /** set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  ls[0].color = BLACK;  ls[0].style = DOTTED_LINE;  ls[0].width = MIDDLE_LINE;
  ls[1].color = BLACK;  ls[1].style =  SOLID_LINE;  ls[1].width = MIDDLE_LINE;
  sprintf(pt[0].type, PLplotSymbolType[openCircle]);  pt[0].scale = PLplotSymbolSize[openCircle];
  sprintf(pt[1].type, PLplotSymbolType[fillSquare]);  pt[1].scale = PLplotSymbolSize[fillSquare];
#   if  NARC > 2
  ls[2].color = BLACK;  ls[2].style = DASHED_LINE;  ls[2].width = MIDDLE_LINE;
  sprintf(pt[2].type, PLplotSymbolType[openDiamond]);  pt[2].scale = PLplotSymbolSize[openDiamond];
#endif//NARC > 2
  for(int ii = 0; ii < kind; ii++)
    pt[ii].color = ls[ii].color;

  /** set labels */
  char hlab[PLplotCharWords];  sprintf(hlab, "#fiN#fr#dtot#u");
  char slab[PLplotCharWords];  sprintf(slab, "Elapsed time per step (s)");
  char tlab[PLplotCharWords];  sprintf(tlab, "Total elapsed time (s)");

  /** set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write = true;
  sprintf(basecap.side, "%s", "t");

  /** set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = true;
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
  for(int ii = 0; ii < NARC; ii++)
    sprintf(pleg.text[ii], gpuname[ii]);
  /* for(int ii = 0; ii < NARC; ii++) */
  /*   for(int jj = 0; jj < NDST; jj++) */
  /*     sprintf(pleg.text[INDEX2D(NARC, NDST, ii, jj)], "%s on %s", dstname[jj], gpuname[ii]); */
  pleg.pos = PL_POSITION_RIGHT | PL_POSITION_BOTTOM | PL_POSITION_INSIDE;


  /** memory allocation to exploit multiple-plot function */
  /** maximum number of panels in each direction */
  const PLINT nxpanel = NDST;
  const PLINT nypanel = 1;

  /** memory allocation for data */
  PLFLT **hor;  allocPointer4PLFLT(&hor, nxpanel * nypanel * kind);
  PLFLT **ver;  allocPointer4PLFLT(&ver, nxpanel * nypanel * kind);

  /** specify plot or skip panel */
  PLBOOL *plot;  allocPLBOOL(&plot, nxpanel * nypanel);
  /** arrays related to draw line(s) and point(s) */
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
  /** arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /** arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /** arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /** arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /** arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /** array to set figure name */
  char figfile[PLplotCharWords];
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel);
  assignChar4PLplot(nxpanel, figname, _figname);


  /** configure to plot elapsed time as a function of number of N-body particles */
  /** data preparation */
  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < kind; jj++){
      hor[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &nbody[jj][ii][0];
      ver[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = & time[jj][ii][0];
    }/* for(int jj = 0; jj < kind; jj++){ */


  /** common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

      /** global setting(s) */
      plot[idx] = true;

      /** line setting(s) */
      nlkind[idx] = lkind;
      line  [idx] = ls;
      lnum  [idx] = num;
      lx    [idx] = &hor[INDEX2D(nxpanel * nypanel, lkind, idx, 0)];
      ly    [idx] = &ver[INDEX2D(nxpanel * nypanel, lkind, idx, 0)];

      /** point setting(s) */
      npkind[idx] = pkind;
      point [idx] = pt;
      pnum  [idx] = num;
      px    [idx] = &hor[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];
      py    [idx] = &ver[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];

      /** errorbar setting(s) */
      errbar[idx] = false;

      /** plot area */
      range[idx] = box;

      /** label setting(s) */
      xlabel[idx] = hlab;
      ylabel[idx] = averaged ? slab : tlab;

      /** misc(s) */
      leg[idx] = pleg;
      uni[idx] = puni;
      cap[idx] = basecap;
      /** 97 is "a" in ASCII code */
      sprintf(cap[idx].text, "(%c) %s", 97 + INDEX2D(nypanel, nxpanel, nypanel - jj - 1, ii), dstname[ii]);
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */


  /** individual file names */
  sprintf(figfile, "%s_%s", file, "vs_ntot");
  for(int ii = 0; ii < nxpanel; ii++)
    sprintf(figname[ii], "%s_%s_%s", icname[ii], file, "vs_ntot");


  /** create figure(s) */
  for(int ii = 0; ii < nxpanel; ii++){
    const PLINT   idx = INDEX2D(nxpanel, nypanel, ii, 0);
    const PLBOOL ctmp = cap[idx].write;
    const PLINT  ltmp = leg[idx].pos;
    cap[idx].write = false;
    leg[idx].pos   = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;

    plotData(1, 1, &plot[idx], false, false,
	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
	     &xlabel[idx], &ylabel[idx], "", figname[ii], argc, argv);

    cap[idx].write = ctmp;
    leg[idx].pos   = ltmp;
  }/* for(int ii = 0; ii < nxpanel; ii++){ */

  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      if( idx != INDEX2D(nxpanel, nypanel, 0, nypanel - 1) )
	leg[idx].write = false;
    }/* for(int jj = 0; jj < nypanel; jj++){ */
  plotData(nxpanel, nypanel, plot, true, true,
	   nlkind,  line, lnum, lx, ly,
	   npkind, point, pnum, px, py,
	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
	   cap, leg, uni, range,
	   xlabel, ylabel, "", figfile, argc, argv);


  /** memory deallocation */
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  free(figname);  free(_figname);

  free(plegTxt);

  free(num);  free(hor);  free(ver);
  free(ls);  free(pt);


  __NOTE__("%s\n", "end");
}


void plotNtotBreakdown
(PLINT nbrk[restrict][NDST][NNUM], PLFLT tsum[restrict][NDST][NNUM], PLFLT walk[restrict][NDST][NNUM], PLFLT node[restrict][NDST][NNUM], PLFLT make[restrict][NDST][NNUM], PLFLT time[restrict][NDST][NNUM], PLplotPltRange box, int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /** global setting(s) */
  /** maximum number of data kinds */
  const PLINT pkind = 5;/**< sum, walk, node, make, and time */
  const PLINT lkind = 5;

  const PLBOOL puni = true;


  /** preparation for dataset */
  /** memory allocation for index */
  PLINT kind = (puni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  if( (pkind != kind) || (lkind != kind) ){
    __KILL__(stderr, "ERROR: pkind(%d), lkind(%d), and kind(%d) must be the same\n", pkind, lkind, kind);
  }/* if( (pkind != kind) || (lkind != kind) ){ */
  PLINT *num;  allocPLINT(&num, kind);

  /** set number of data points */
  for(PLINT ii = 0; ii < kind; ii++)    num[ii] = NNUM;

  /** data preparation */
  static PLFLT fbrk[NARC][NDST][NNUM];
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < NDST; jj++)
      for(int kk = 0; kk < NNUM; kk++)
	fbrk[ii][jj][kk] = log10((double)nbrk[ii][jj][kk]);

  /** set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  ls[0].color =   BLACK;  ls[0].style =             SOLID_LINE;  ls[0].width =   BOLD_LINE;
  ls[1].color =     RED;  ls[1].style =            DASHED_LINE;  ls[1].width = MIDDLE_LINE;
  ls[2].color =    BLUE;  ls[2].style =            DOTTED_LINE;  ls[2].width = MIDDLE_LINE;
  ls[3].color = MAGENTA;  ls[3].style =        DOT_DASHED_LINE;  ls[3].width = MIDDLE_LINE;
  ls[4].color =   GREEN;  ls[4].style = TRIPLE_DOT_DASHED_LINE;  ls[4].width = MIDDLE_LINE;
  sprintf(pt[0].type, PLplotSymbolType[smallDot    ]);  pt[0].scale = PLplotSymbolSize[smallDot    ];
  sprintf(pt[1].type, PLplotSymbolType[fillCircle  ]);  pt[1].scale = PLplotSymbolSize[fillSquare  ];
  sprintf(pt[2].type, PLplotSymbolType[fillSquare  ]);  pt[2].scale = PLplotSymbolSize[fillSquare  ];
  sprintf(pt[3].type, PLplotSymbolType[fillDiamond ]);  pt[3].scale = PLplotSymbolSize[fillDiamond ];
  sprintf(pt[4].type, PLplotSymbolType[fillTriangle]);  pt[4].scale = PLplotSymbolSize[fillTriangle];
  for(int ii = 0; ii < kind; ii++)
    pt[ii].color = ls[ii].color;

  /** set labels */
  char hlab[PLplotCharWords];  sprintf(hlab, "#fiN#fr#dtot#u");
  char slab[PLplotCharWords];  sprintf(slab, "Elapsed time per step (s)");

  /** set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write = true;
  sprintf(basecap.side, "%s", "t");

  /** set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = true;
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
  sprintf(pleg.text[0], "total");
  sprintf(pleg.text[1], "walk tree");
  sprintf(pleg.text[2], "calc node");
  sprintf(pleg.text[3], "make tree");
  sprintf(pleg.text[4], "pred/corr");
  pleg.pos = PL_POSITION_RIGHT | PL_POSITION_BOTTOM | PL_POSITION_INSIDE;


  /** memory allocation to exploit multiple-plot function */
  /** maximum number of panels in each direction */
  const PLINT nxpanel = NARC;
  const PLINT nypanel = NDST;

  /** memory allocation for data */
  PLFLT **hor;  allocPointer4PLFLT(&hor, nxpanel * nypanel * kind);
  PLFLT **ver;  allocPointer4PLFLT(&ver, nxpanel * nypanel * kind);

  /** specify plot or skip panel */
  PLBOOL *plot;  allocPLBOOL(&plot, nxpanel * nypanel);
  /** arrays related to draw line(s) and point(s) */
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
  /** arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /** arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /** arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /** arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /** arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /** array to set figure name */
  char figfile[PLplotCharWords];
  char *_figname;  allocChar4PLplot        (&_figname, NDST);
  char **figname;  allocPointer4Char4PLplot(& figname, NDST);
  assignChar4PLplot(NDST, figname, _figname);


  /** configure to plot elapsed tsum as a function of number of N-body particles */
  /** data preparation */
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < NDST; jj++){
      hor[INDEX(NARC, NDST, kind, ii, jj, 0)] = &fbrk[ii][NDST - jj - 1][0];      ver[INDEX(NARC, NDST, kind, ii, jj, 0)] = &tsum[ii][NDST - jj - 1][0];
      hor[INDEX(NARC, NDST, kind, ii, jj, 1)] = &fbrk[ii][NDST - jj - 1][0];      ver[INDEX(NARC, NDST, kind, ii, jj, 1)] = &walk[ii][NDST - jj - 1][0];
      hor[INDEX(NARC, NDST, kind, ii, jj, 2)] = &fbrk[ii][NDST - jj - 1][0];      ver[INDEX(NARC, NDST, kind, ii, jj, 2)] = &node[ii][NDST - jj - 1][0];
      hor[INDEX(NARC, NDST, kind, ii, jj, 3)] = &fbrk[ii][NDST - jj - 1][0];      ver[INDEX(NARC, NDST, kind, ii, jj, 3)] = &make[ii][NDST - jj - 1][0];
      hor[INDEX(NARC, NDST, kind, ii, jj, 4)] = &fbrk[ii][NDST - jj - 1][0];      ver[INDEX(NARC, NDST, kind, ii, jj, 4)] = &time[ii][NDST - jj - 1][0];
    }/* for(int jj = 0; jj < NDST; jj++){ */


  /** common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

      /** global setting(s) */
      plot[idx] = true;

      /** line setting(s) */
      nlkind[idx] = lkind;
      line  [idx] = ls;
      lnum  [idx] = num;
      lx    [idx] = &hor[INDEX2D(nxpanel * nypanel, lkind, idx, 0)];
      ly    [idx] = &ver[INDEX2D(nxpanel * nypanel, lkind, idx, 0)];

      /** point setting(s) */
      npkind[idx] = pkind;
      point [idx] = pt;
      pnum  [idx] = num;
      px    [idx] = &hor[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];
      py    [idx] = &ver[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];

      /** errorbar setting(s) */
      errbar[idx] = false;

      /** plot area */
      range[idx] = box;

      /** label setting(s) */
      xlabel[idx] = hlab;
      ylabel[idx] = slab;

      /** misc(s) */
      leg[idx] = pleg;
      uni[idx] = puni;
      cap[idx] = basecap;
      /** 97 is "a" in ASCII code */
      sprintf(cap[idx].text, "(%c) %s on %s", 97 + INDEX2D(nypanel, nxpanel, nypanel - jj - 1, ii), dstname[NDST - jj - 1], gpuname[ii]);
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */


  /** individual file names */
  sprintf(figfile, "%s", "breakdown_ndep");
  for(int ii = 0; ii < NDST; ii++)
    sprintf(figname[ii], "%s_%s", icname[ii], "breakdown_ndep");


  /** create figure(s) */
  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      if( idx != INDEX2D(nxpanel, nypanel, 0, nypanel - 1) )
	leg[idx].write = false;
    }/* for(int jj = 0; jj < nypanel; jj++){ */
  plotData(nxpanel, nypanel, plot, true, true,
	   nlkind,  line, lnum, lx, ly,
	   npkind, point, pnum, px, py,
	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
	   cap, leg, uni, range,
	   xlabel, ylabel, "", figfile, argc, argv);

  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      sprintf(cap[idx].text, "(%c) %s", 97 + ii, gpuname[ii]);
      if( ii == 0 )
	leg[idx].write = true;
    }/* for(int jj = 0; jj < nypanel; jj++){ */
  for(int ii = 0; ii < nxpanel; ii++)
    plot[INDEX2D(nxpanel, nypanel, ii, 1)] = false;
  plotData(nxpanel, nypanel, plot, true, true,
	   nlkind,  line, lnum, lx, ly,
	   npkind, point, pnum, px, py,
	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
	   cap, leg, uni, range,
	   xlabel, ylabel, "", figname[0], argc, argv);
  for(int ii = 0; ii < nxpanel; ii++){
    const PLINT src = INDEX2D(nxpanel, nypanel, ii, 1);
    const PLINT dst = INDEX2D(nxpanel, nypanel, ii, 0);

    nlkind[dst] = nlkind[src];
    line  [dst] =   line[src];
    lnum  [dst] =   lnum[src];
    lx    [dst] =     lx[src];
    ly    [dst] =     ly[src];

    npkind[dst] = npkind[src];
    point [dst] =  point[src];
    pnum  [dst] =   pnum[src];
    px    [dst] =     px[src];
    py    [dst] =     py[src];

    errbar[dst] = errbar[src];

    range[dst] = range[src];

    xlabel[dst] = xlabel[src];
    ylabel[dst] = ylabel[src];

    leg[dst] = leg[src];
    uni[dst] = uni[src];
    cap[dst] = cap[src];
  }/* for(int ii = 0; ii < nxpanel; ii++){ */
  plotData(nxpanel, nypanel, plot, true, true,
	   nlkind,  line, lnum, lx, ly,
	   npkind, point, pnum, px, py,
	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
	   cap, leg, uni, range,
	   xlabel, ylabel, "", figname[1], argc, argv);


  /** memory deallocation */
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  free(figname);  free(_figname);

  free(plegTxt);

  free(num);  free(hor);  free(ver);
  free(ls);  free(pt);


  __NOTE__("%s\n", "end");
}


void plotNiDependence(PLFLT inum[restrict][NDST][NSUB][NMAX], PLFLT time[restrict][NDST][NSUB][NMAX], PLplotPltRange box, int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /** global setting(s) */
  /** maximum number of data kinds */
#ifdef  LINE_ONLY_PLOT
  const PLINT pkind = 0;
#else///LINE_ONLY_PLOT
  const PLINT pkind = NSUB;
#endif//LINE_ONLY_PLOT
  const PLINT lkind = NSUB;

  const PLBOOL puni = true;


  /** preparation for dataset */
  /** memory allocation for index */
  PLINT kind = (puni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
#ifndef LINE_ONLY_PLOT
  if( (pkind != kind) || (lkind != kind) ){
    __KILL__(stderr, "ERROR: pkind(%d), lkind(%d), and kind(%d) must be the same\n", pkind, lkind, kind);
  }/* if( (pkind != kind) || (lkind != kind) ){ */
#endif//LINE_ONLY_PLOT
  PLINT *num;  allocPLINT(&num, kind);

  /** set number of data points */
  for(PLINT ii = 0; ii < kind; ii++)    num[ii] = NMAX;

  /** set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  ls[0].style =            DASHED_LINE;  ls[0].width = MIDDLE_LINE;
  ls[1].style =             SOLID_LINE;  ls[1].width = MIDDLE_LINE;
  ls[2].style =            DOTTED_LINE;  ls[2].width = MIDDLE_LINE;
  ls[3].style =        DOT_DASHED_LINE;  ls[3].width = MIDDLE_LINE;
  ls[4].style = TRIPLE_DOT_DASHED_LINE;  ls[4].width = MIDDLE_LINE;
  ls[5].style =             SOLID_LINE;  ls[5].width = MIDDLE_LINE;
#ifndef LINE_ONLY_PLOT
  for(int ii = 0; ii < kind; ii++)
    pt[ii].color = ls[ii].color;
#endif//LINE_ONLY_PLOT

  /** set labels */
  char hlab[PLplotCharWords];  sprintf(hlab, "#fiN#di#u#fr");
  char vlab[PLplotCharWords];  sprintf(vlab, "Elapsed time (s)");

  /** set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write = true;
  sprintf(basecap.side, "%s", "b");

  /** set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = true;
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
  for(int jj = 0; jj < NSUB; jj++)
    sprintf(pleg.text[jj], "#fiS#fr = %d", 1 << jj);
  pleg.pos = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;


  /** memory allocation to exploit multiple-plot function */
  /** maximum number of panels in each direction */
  const PLINT nxpanel = NARC;
  const PLINT nypanel = NDST;

  /** memory allocation for data */
  PLFLT **hor;  allocPointer4PLFLT(&hor, nxpanel * nypanel * kind);
  PLFLT **ver;  allocPointer4PLFLT(&ver, nxpanel * nypanel * kind);

  /** specify plot or skip panel */
  PLBOOL *plot;  allocPLBOOL(&plot, nxpanel * nypanel);
  /** arrays related to draw line(s) and point(s) */
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
  /** arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /** arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /** arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /** arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /** arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /** array to set figure name */
  char figfile[PLplotCharWords];
  char *_figname;  allocChar4PLplot        (&_figname, (nxpanel + 1) * nypanel);
  char **figname;  allocPointer4Char4PLplot(& figname, (nxpanel + 1) * nypanel);
  assignChar4PLplot((nxpanel + 1) * nypanel, figname, _figname);


  /** configure to plot elapsed time as a function of number of N-body particles */
  /** data preparation */
  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel; jj++)
      for(int kk = 0; kk < NSUB; kk++){
	hor[INDEX(nxpanel, nypanel, NSUB, ii, jj, kk)] = inum[ii][jj][kk];
	ver[INDEX(nxpanel, nypanel, NSUB, ii, jj, kk)] = time[ii][jj][kk];
      }/* for(int kk = 0; kk < NSUB; kk++){ */


  /** common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

      /** global setting(s) */
      plot[idx] = true;

      /** line setting(s) */
      nlkind[idx] = lkind;
      line  [idx] = ls;
      lnum  [idx] = num;
      lx    [idx] = &hor[INDEX2D(nxpanel * nypanel, lkind, idx, 0)];
      ly    [idx] = &ver[INDEX2D(nxpanel * nypanel, lkind, idx, 0)];

      /** point setting(s) */
      npkind[idx] = pkind;
#ifdef  LINE_ONLY_PLOT
      point [idx] = NULL;
      pnum  [idx] = NULL;
      px    [idx] = NULL;
      py    [idx] = NULL;
#else///LINE_ONLY_PLOT
      point [idx] = pt;
      pnum  [idx] = num;
      px    [idx] = &hor[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];
      py    [idx] = &ver[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];
#endif//LINE_ONLY_PLOT

      /** errorbar setting(s) */
      errbar[idx] = false;

      /** plot area */
      range[idx] = box;

      /** label setting(s) */
      xlabel[idx] = hlab;
      ylabel[idx] = vlab;

      /** misc(s) */
      leg[idx] = pleg;
      uni[idx] = puni;
      cap[idx] = basecap;
      /** 97 is "a" in ASCII code */
      sprintf(cap[idx].text, "(%c) %s on %s", 97 + INDEX2D(nypanel, nxpanel, nypanel - jj - 1, ii), dstname[jj], gpuname[ii]);
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */


  /** individual file names */
  sprintf(figfile, "%s", "nidep");
  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel; jj++)
      sprintf(figname[INDEX2D(nxpanel, nypanel, ii, jj)], "%s_%s_%s", icname[jj], "nidep", arcname[ii]);
  for(int jj = 0; jj < nypanel; jj++)
    sprintf(figname[nxpanel * nypanel + jj], "%s_%s", icname[jj], "nidep");


  /** create figure(s) */
  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      const PLBOOL ctmp = cap[idx].write;
      const PLINT  ltmp = leg[idx].pos;
      cap[idx].write = false;
      leg[idx].pos   = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;

      plotData(1, 1, &plot[idx], false, false,
	       &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
	       &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
	       &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
	       &cap[idx], &leg[idx], &uni[idx], &range[idx],
	       &xlabel[idx], &ylabel[idx], "", figname[idx], argc, argv);

      /* cap[idx].write = tmp; */
      cap[idx].write = ctmp;
      leg[idx].pos   = ltmp;
    }/* for(int jj = 0; jj < nypanel; jj++){ */

  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      if( idx != INDEX2D(nxpanel, nypanel, 0, nypanel - 1) )
	leg[idx].write = false;
    }/* for(int jj = 0; jj < nypanel; jj++){ */
  plotData(nxpanel, nypanel, plot, true, true,
	   nlkind,  line, lnum, lx, ly,
	   npkind, point, pnum, px, py,
	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
	   cap, leg, uni, range,
	   xlabel, ylabel, "", figfile, argc, argv);


  /** data preparation */
  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel; jj++)
      for(int kk = 0; kk < NSUB; kk++){
	hor[INDEX(nypanel, nxpanel, NSUB, jj, ii, kk)] = inum[ii][jj][kk];
	ver[INDEX(nypanel, nxpanel, NSUB, jj, ii, kk)] = time[ii][jj][kk];
      }/* for(int kk = 0; kk < NSUB; kk++){ */

  /** common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nypanel, nxpanel, jj, ii);

      /** misc(s) */
      leg[idx].write = (ii == 0) ? true : false;
      /** 97 is "a" in ASCII code */
      sprintf(cap[idx].text, "(%c) %s", 97 + ii, gpuname[ii]);
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */

  for(int jj = 0; jj < nypanel; jj++){
    const PLINT idx = INDEX2D(nypanel, nxpanel, jj, 0);
    plotData(nxpanel, 1, &plot[idx], true, true,
	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
	     &xlabel[idx], &ylabel[idx], "", figname[nxpanel * nypanel + jj], argc, argv);
  }/* for(int jj = 0; jj < nypanel; jj++){ */


  /** memory deallocation */
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  free(figname);  free(_figname);

  free(plegTxt);

  free(num);  free(hor);  free(ver);
  free(ls);  free(pt);


  __NOTE__("%s\n", "end");
}
