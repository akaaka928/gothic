/*************************************************************************\
 *                                                                       *
                  last updated on 2016/06/13(Mon) 12:52:37
 *                                                                       *
 *    Plot Code of Breakdown of Tree code on K20X with acceleration MAC  *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef EXEC_BENCHMARK
#define EXEC_BENCHMARK
#endif//EXEC_BENCHMARK
//-------------------------------------------------------------------------
#define SHOW_LEGEND_OF_SYMBOLS
//-------------------------------------------------------------------------
#define CONNECT_SYMBOLS_BY_LINE
#define HIGHLIGHTING_COLOR_MAP
//-------------------------------------------------------------------------
/* #define CUMULATIVE_PLOT */
//-------------------------------------------------------------------------
#ifdef  HIGHLIGHTING_COLOR_MAP
#define NSEQUENCES (3)
#endif//HIGHLIGHTING_COLOR_MAP
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <myutil.h>
#include <name.h>
#include <plplotlib.h>
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotBreakdown
(PLINT step[restrict],
 PLFLT walkTree[restrict], PLFLT makeTree[restrict], PLFLT sortPH[restrict], PLFLT calcMAC[restrict],
 PLFLT spltIgrp[restrict], PLFLT predJpar[restrict], PLFLT corrIpar[restrict],
#ifdef  CUMULATIVE_PLOT
 PLFLT walkTreeSum[restrict], PLFLT makeTreeSum[restrict], PLFLT sortPHSum[restrict], PLFLT calcMACSum[restrict],
 PLFLT spltIgrpSum[restrict], PLFLT predJparSum[restrict], PLFLT corrIparSum[restrict],
#endif//CUMULATIVE_PLOT
 PLplotPltRange box, char file[], int argc, char **argv);
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
  /* read the measured data */
  static PLINT step[BENCHMARK_STEPS];
  static PLFLT walkTree[BENCHMARK_STEPS], calcMAC[BENCHMARK_STEPS], makeTree[BENCHMARK_STEPS], sortPH[BENCHMARK_STEPS];
  static PLFLT getBody[BENCHMARK_STEPS], setBody[BENCHMARK_STEPS], copyNode[BENCHMARK_STEPS], copyCell[BENCHMARK_STEPS];
  static PLFLT predJpar[BENCHMARK_STEPS], corrIpar[BENCHMARK_STEPS], setGlbDt[BENCHMARK_STEPS], setLocDt[BENCHMARK_STEPS], adjustDt[BENCHMARK_STEPS];
  static PLFLT calcIsep[BENCHMARK_STEPS], spltIgrp[BENCHMARK_STEPS];
#ifdef  CUMULATIVE_PLOT
  static PLFLT walkTreeSum[BENCHMARK_STEPS], calcMACSum[BENCHMARK_STEPS], makeTreeSum[BENCHMARK_STEPS], sortPHSum[BENCHMARK_STEPS];
  static PLFLT getBodySum[BENCHMARK_STEPS], setBodySum[BENCHMARK_STEPS], copyNodeSum[BENCHMARK_STEPS], copyCellSum[BENCHMARK_STEPS];
  static PLFLT predJparSum[BENCHMARK_STEPS], corrIparSum[BENCHMARK_STEPS], setGlbDtSum[BENCHMARK_STEPS], setLocDtSum[BENCHMARK_STEPS], adjustDtSum[BENCHMARK_STEPS];
  static PLFLT calcIsepSum[BENCHMARK_STEPS], spltIgrpSum[BENCHMARK_STEPS];
#endif//CUMULATIVE_PLOT
  FILE *fp;
  char filename[128];
  sprintf(filename, "%s/%s.time.bare.dat", BENCH_LOG_FOLDER, file);
  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  //-----------------------------------------------------------------------
  int checker = 1;
  for(int ii = 0; ii < BENCHMARK_STEPS; ii++)
    checker &= (16 == fscanf(fp, "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le", &step[ii],
			     &walkTree[ii], &calcMAC[ii], &makeTree[ii], &setGlbDt[ii], &sortPH[ii],
			     &getBody[ii], &setBody[ii], &copyNode[ii], &copyCell[ii],
			     &predJpar[ii], &corrIpar[ii], &setLocDt[ii], &adjustDt[ii], &calcIsep[ii], &spltIgrp[ii]));
  //-----------------------------------------------------------------------
  fclose(fp);
  if( checker != 1 ){
    __KILL__(stderr, "ERROR: failure to read data from \"%s\"\n", filename);
  }/* if( checker != 1 ){ */
  //-----------------------------------------------------------------------
#ifdef  CUMULATIVE_PLOT
  walkTreeSum[0] = walkTree[0];
  calcMACSum [0] = calcMAC [0];
  makeTreeSum[0] = makeTree[0];
  setGlbDtSum[0] = setGlbDt[0];
  sortPHSum  [0] = sortPH  [0];
  getBodySum [0] =  getBody[0];
  setBodySum [0] =  setBody[0];
  copyNodeSum[0] = copyNode[0];
  copyCellSum[0] = copyCell[0];
  predJparSum[0] = predJpar[0];
  corrIparSum[0] = corrIpar[0];
  setLocDtSum[0] = setLocDt[0];
  adjustDtSum[0] = adjustDt[0];
  calcIsepSum[0] = calcIsep[0];
  spltIgrpSum[0] = spltIgrp[0];
  for(int ii = 1; ii < BENCHMARK_STEPS; ii++){
    walkTreeSum[ii] = walkTree[ii] + walkTreeSum[ii - 1];
    calcMACSum [ii] = calcMAC [ii] + calcMACSum [ii - 1];
    makeTreeSum[ii] = makeTree[ii] + makeTreeSum[ii - 1];
    setGlbDtSum[ii] = setGlbDt[ii] + setGlbDtSum[ii - 1];
    sortPHSum  [ii] = sortPH  [ii] + sortPHSum  [ii - 1];
    getBodySum [ii] =  getBody[ii] + getBodySum [ii - 1];
    setBodySum [ii] =  setBody[ii] + setBodySum [ii - 1];
    copyNodeSum[ii] = copyNode[ii] + copyNodeSum[ii - 1];
    copyCellSum[ii] = copyCell[ii] + copyCellSum[ii - 1];
    predJparSum[ii] = predJpar[ii] + predJparSum[ii - 1];
    corrIparSum[ii] = corrIpar[ii] + corrIparSum[ii - 1];
    setLocDtSum[ii] = setLocDt[ii] + setLocDtSum[ii - 1];
    adjustDtSum[ii] = adjustDt[ii] + adjustDtSum[ii - 1];
    calcIsepSum[ii] = calcIsep[ii] + calcIsepSum[ii - 1];
    spltIgrpSum[ii] = spltIgrp[ii] + spltIgrpSum[ii - 1];
  }/* for(int ii = 1; ii < BENCHMARK_STEPS; ii++){ */
#endif//CUMULATIVE_PLOT
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < BENCHMARK_STEPS; ii++){
    //---------------------------------------------------------------------
    walkTree[ii] =                         log10(walkTree[ii])          ;/* time to calculate gravity (calcGrav_dev) */
    makeTree[ii] = (makeTree[ii] != 0.0) ? log10(makeTree[ii]) : DBL_MAX;/* time to build tree structure (makeTree) */
    sortPH  [ii] = (sortPH  [ii] != 0.0) ? log10(sortPH  [ii]) : DBL_MAX;/* time to generate and sort PH-key and sort N-body particles (sortPHcurve) */
    calcMAC [ii] =                         log10(calcMAC [ii])          ;/* time to calculate mass, size and MAC for pseudo j-particles (calcMultipole) */
    spltIgrp[ii] = (spltIgrp[ii] != 0.0) ? log10(spltIgrp[ii]) : DBL_MAX;/* time to slice group of i-particles (searchNeighbor_dev) */
    predJpar[ii] =                         log10(predJpar[ii])          ;/* time to calculate position and velocity of j-particles (predoction_dev) */
    corrIpar[ii] = (corrIpar[ii] != 0.0) ? log10(corrIpar[ii]) : DBL_MAX;/* time to calculate position of i-particles (correction_dev) */
    calcIsep[ii] = (calcIsep[ii] != 0.0) ? log10(calcIsep[ii]) : DBL_MAX;/* time to calculate separation between i-particles within a group (examineNeighbor_dev) */
    setGlbDt[ii] = (setGlbDt[ii] != 0.0) ? log10(setGlbDt[ii]) : DBL_MAX;/* time to set global time step of the simulation (setTimeStep_dev) */
    setLocDt[ii] =                         log10(setLocDt[ii])          ;/* time to set  local time step of the simulation (setLaneTime_dev) */
    adjustDt[ii] =                         log10(adjustDt[ii])          ;/* time to adjust time step of N-body particles (adjustTime_dev) */
    setBody [ii] = (setBody [ii] != 0.0) ? log10( setBody[ii]) : DBL_MAX;/* time to copy N-body particles from host to device (cpBody_hst2dev) */
    getBody [ii] = (getBody [ii] != 0.0) ? log10( getBody[ii]) : DBL_MAX;/* time to copy N-body particles from device to host (cpBody_dev2hst) */
    copyNode[ii] = (copyNode[ii] != 0.0) ? log10(copyNode[ii]) : DBL_MAX;/* time to copy tree node from host to device (setTreeNode_dev) */
    copyCell[ii] = (copyCell[ii] != 0.0) ? log10(copyCell[ii]) : DBL_MAX;/* time to copy tree cell from host to device (setTreeNode_dev) */
    //---------------------------------------------------------------------
/* #ifdef  CUMULATIVE_PLOT */
/*     walkTreeSum[ii] = log10(walkTreeSum[ii]); */
/*     calcMACSum [ii] = log10(calcMACSum [ii]); */
/*     makeTreeSum[ii] = log10(makeTreeSum[ii]); */
/*     setGlbDtSum[ii] = log10(setGlbDtSum[ii]); */
/*     sortPHSum  [ii] = log10(sortPHSum  [ii]); */
/*     getBodySum [ii] = log10(getBodySum [ii]); */
/*     setBodySum [ii] = log10(setBodySum [ii]); */
/*     copyNodeSum[ii] = log10(copyNodeSum[ii]); */
/*     copyCellSum[ii] = log10(copyCellSum[ii]); */
/*     predJparSum[ii] = log10(predJparSum[ii]); */
/*     corrIparSum[ii] = log10(corrIparSum[ii]); */
/*     setLocDtSum[ii] = log10(setLocDtSum[ii]); */
/*     adjustDtSum[ii] = log10(adjustDtSum[ii]); */
/*     calcIsepSum[ii] = log10(calcIsepSum[ii]); */
/*     spltIgrpSum[ii] = log10(spltIgrpSum[ii]); */
/* #endif//CUMULATIVE_PLOT */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < BENCHMARK_STEPS; ii++){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* plot breakdown as a function of step */
  //-----------------------------------------------------------------------
  PLplotPltRange range;
  range.xmin =   0.0;
  /* range.xmax = 300.0; */
  range.xmax = 100.0;
  range.xlog = LINEAR_PLOT;
#if 1
  range.ymin = log10(1.0e-5);
#else
  range.ymin = log10(1.0e-3);
#endif
  range.ymax = log10(1.0e+1);
  range.ylog = LOGARITHMIC_PLOT;
  range.xgrd = range.ygrd = true;
#if 1
  /* const PLFLT length = 1.0; */
  const PLFLT length = 0.5;
  range.xmin -= length;
  range.xmax += length;
#endif
  //-----------------------------------------------------------------------
  plotBreakdown(step, walkTree, makeTree, sortPH, calcMAC, spltIgrp, predJpar, corrIpar,
#ifdef  CUMULATIVE_PLOT
		walkTreeSum, makeTreeSum, sortPHSum, calcMACSum, spltIgrpSum, predJparSum, corrIparSum,
#endif//CUMULATIVE_PLOT
		range, file, argc, argv);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotBreakdown
(PLINT step[restrict],
 PLFLT walkTree[restrict], PLFLT makeTree[restrict], PLFLT sortPH[restrict], PLFLT calcMAC[restrict],
 PLFLT spltIgrp[restrict], PLFLT predJpar[restrict], PLFLT corrIpar[restrict],
#ifdef  CUMULATIVE_PLOT
 PLFLT walkTreeSum[restrict], PLFLT makeTreeSum[restrict], PLFLT sortPHSum[restrict], PLFLT calcMACSum[restrict],
 PLFLT spltIgrpSum[restrict], PLFLT predJparSum[restrict], PLFLT corrIparSum[restrict],
#endif//CUMULATIVE_PLOT
 PLplotPltRange box, char file[], int argc, char **argv)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* global setting(s) */
  //-----------------------------------------------------------------------
  /* maximum number of data kinds */
  const PLINT pkind = 7;
#ifdef  CONNECT_SYMBOLS_BY_LINE
  /* const PLINT lkind = 1; */
  const PLINT lkind = pkind;
#else///CONNECT_SYMBOLS_BY_LINE
  const PLINT lkind = 0;
#endif//CONNECT_SYMBOLS_BY_LINE
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
  for(PLINT ii = 0; ii < kind; ii++)    num[ii] = BENCHMARK_STEPS;
  //-----------------------------------------------------------------------
  /* data preparation */
  static PLFLT xdat[BENCHMARK_STEPS];
  for(int ii = 0; ii < BENCHMARK_STEPS; ii++)
    xdat[ii] = (PLFLT)step[ii];
#ifdef  CONNECT_SYMBOLS_BY_LINE
  PLINT numTreeBuild = 0;
  for(int ii = 0; ii < BENCHMARK_STEPS; ii++)
    if( makeTree[ii] < DBL_MAX )
      numTreeBuild++;
  PLFLT *xtree;  allocPLFLT(&xtree, numTreeBuild);
  PLFLT *ymake;  allocPLFLT(&ymake, numTreeBuild);
  PLFLT *ysort;  allocPLFLT(&ysort, numTreeBuild);
  PLFLT *ysplt;  allocPLFLT(&ysplt, numTreeBuild);
#ifdef  CUMULATIVE_PLOT
  PLFLT *xtreeSum;  allocPLFLT(&xtreeSum, numTreeBuild);
  PLFLT *ymakeSum;  allocPLFLT(&ymakeSum, numTreeBuild);
  PLFLT *ysortSum;  allocPLFLT(&ysortSum, numTreeBuild);
  PLFLT *yspltSum;  allocPLFLT(&yspltSum, numTreeBuild);
#endif//CUMULATIVE_PLOT
  numTreeBuild = 0;
  for(int ii = 0; ii < BENCHMARK_STEPS; ii++)
    if( makeTree[ii] < DBL_MAX ){
      xtree[numTreeBuild] = xdat[ii];
      ymake[numTreeBuild] = makeTree[ii];
      ysort[numTreeBuild] = sortPH  [ii];
      ysplt[numTreeBuild] = spltIgrp[ii];
#ifdef  CUMULATIVE_PLOT
      xtreeSum[numTreeBuild] = xdat[ii];
      ymakeSum[numTreeBuild] = makeTreeSum[ii];
      ysortSum[numTreeBuild] =   sortPHSum[ii];
      yspltSum[numTreeBuild] = spltIgrpSum[ii];
#endif//CUMULATIVE_PLOT
      numTreeBuild++;
    }/* if( makeTree[ii] < DBL_MAX ){ */
#endif//CONNECT_SYMBOLS_BY_LINE
  //-----------------------------------------------------------------------
  /* set symbol and line style */
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  /* pt[0].color =     RED;  sprintf(pt[0].type, PLplotSymbolType[openCircle  ]);  pt[0].scale = PLplotSymbolSize[openCircle  ]; */
  pt[0].color =     RED;  sprintf(pt[0].type, PLplotSymbolType[fillCircle  ]);  pt[0].scale = PLplotSymbolSize[fillCircle  ];
  pt[1].color = MAGENTA;  sprintf(pt[1].type, PLplotSymbolType[fillCircle  ]);  pt[1].scale = PLplotSymbolSize[fillCircle  ];
  pt[2].color =    BLUE;  sprintf(pt[2].type, PLplotSymbolType[fillSquare  ]);  pt[2].scale = PLplotSymbolSize[fillSquare  ];
  pt[3].color =   BLACK;  sprintf(pt[3].type, PLplotSymbolType[    Plus    ]);  pt[3].scale = PLplotSymbolSize[    Plus    ];
  pt[4].color =    BLUE;  sprintf(pt[4].type, PLplotSymbolType[fillDiamond ]);  pt[4].scale = PLplotSymbolSize[fillDiamond ];
  pt[5].color =   GREEN;  sprintf(pt[5].type, PLplotSymbolType[    Plus    ]);  pt[5].scale = PLplotSymbolSize[    Plus    ];
  /* pt[6].color =   BLACK;  sprintf(pt[6].type, PLplotSymbolType[openTriangle]);  pt[6].scale = PLplotSymbolSize[openTriangle]; */
  pt[6].color =   BLACK;  sprintf(pt[6].type, PLplotSymbolType[fillTriangle]);  pt[6].scale = PLplotSymbolSize[fillTriangle];
  for(int ii = 0; ii < pkind; ii++)
    pt[ii].scale *= 0.2;
#ifdef  CONNECT_SYMBOLS_BY_LINE
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  for(int ii = 0; ii < lkind; ii++){
    ls[ii].color = pt[ii].color;
    ls[ii].width = THIN_LINE * 0.1;
    ls[ii].style = SOLID_LINE;
  }/* for(int ii = 0; ii < lkind; ii++){ */
#endif//CONNECT_SYMBOLS_BY_LINE
  //-----------------------------------------------------------------------
  /* set labels */
  char hlab[PLplotCharWords];  sprintf(hlab, "Step");
  char vlab[PLplotCharWords];  sprintf(vlab, "Execution time (s)");
#ifdef  CUMULATIVE_PLOT
  char clab[PLplotCharWords];  sprintf(clab, "Elapsed time (s)");
#endif//CUMULATIVE_PLOT
  //-----------------------------------------------------------------------
  /* set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write = false;
  sprintf(basecap.side, "%s", "t");
  //-----------------------------------------------------------------------
  /* set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = false;
#ifdef  SHOW_LEGEND_OF_SYMBOLS
  pleg.write = true;
#endif//SHOW_LEGEND_OF_SYMBOLS
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
  sprintf(pleg.text[0], "walk tree");
  sprintf(pleg.text[1], "make tree");
  sprintf(pleg.text[2], "PH-key");
  sprintf(pleg.text[3], "calc MAC");
  sprintf(pleg.text[4], "split #fii#fr-group");
  sprintf(pleg.text[5], "predictor");
  sprintf(pleg.text[6], "corrector");
  /* pleg.pos = PL_POSITION_RIGHT | PL_POSITION_TOP | PL_POSITION_OUTSIDE; */
  pleg.pos = PL_POSITION_RIGHT | PL_POSITION_TOP | PL_POSITION_BOTTOM | PL_POSITION_OUTSIDE;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation to exploit multiple-plot function */
  //-----------------------------------------------------------------------
  /* maximum number of panels in each direction */
  const PLINT nxpanel = 1;
#ifdef  CUMULATIVE_PLOT
  const PLINT nypanel = 2;
#else///CUMULATIVE_PLOT
  const PLINT nypanel = 1;
#endif//CUMULATIVE_PLOT
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
  /* configure to plot elapsed time as a function of number of N-body particles */
  //-----------------------------------------------------------------------
  /* data preparation */
  for(int ii = 0; ii < kind; ii++)
    hor[ii] = xdat;
  ver[0] = walkTree;
  ver[1] = makeTree;
  ver[2] = sortPH  ;
  ver[3] = calcMAC ;
  ver[4] = spltIgrp;
  ver[5] = predJpar;
  ver[6] = corrIpar;
#ifdef  CONNECT_SYMBOLS_BY_LINE
  num[1] = numTreeBuild;  hor[1] = xtree;  ver[1] = ymake;
  num[2] = numTreeBuild;  hor[2] = xtree;  ver[2] = ysort;
  num[4] = numTreeBuild;  hor[4] = xtree;  ver[4] = ysplt;
#endif//CONNECT_SYMBOLS_BY_LINE
#ifdef  CUMULATIVE_PLOT
  for(int ii = 0; ii < kind; ii++)
    hor[kind + ii] = xdat;
  ver[kind + 0] = walkTreeSum;
  ver[kind + 1] = makeTreeSum;
  ver[kind + 2] =   sortPHSum;
  ver[kind + 3] =  calcMACSum;
  ver[kind + 4] = spltIgrpSum;
  ver[kind + 5] = predJparSum;
  ver[kind + 6] = corrIparSum;
#ifdef  CONNECT_SYMBOLS_BY_LINE
  hor[kind + 1] = xtree;  ver[kind + 1] = ymakeSum;
  hor[kind + 2] = xtree;  ver[kind + 2] = ysortSum;
  hor[kind + 4] = xtree;  ver[kind + 4] = yspltSum;
#endif//CONNECT_SYMBOLS_BY_LINE
#endif//CUMULATIVE_PLOT
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
#ifdef  CONNECT_SYMBOLS_BY_LINE
      line  [idx] = ls;
      lnum  [idx] = num;
      lx    [idx] = &hor[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];
      ly    [idx] = &ver[INDEX2D(nxpanel * nypanel, pkind, idx, 0)];
#else///CONNECT_SYMBOLS_BY_LINE
      line  [idx] = NULL;
      lnum  [idx] = NULL;
      lx    [idx] = NULL;
      ly    [idx] = NULL;
#endif//CONNECT_SYMBOLS_BY_LINE
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
#ifdef  CUMULATIVE_PLOT
      if( idx != 0 ){
	range[idx].ylog = LINEAR_PLOT;
	range[idx].ymin = 0.0;
	range[idx].ymax = 40.0;
      }
#endif//CUMULATIVE_PLOT
      //-------------------------------------------------------------------
      /* label setting(s) */
      xlabel[idx] = hlab;
      ylabel[idx] = vlab;
#ifdef  CUMULATIVE_PLOT
      if( idx != 0 )
	ylabel[idx] = clab;
#endif//CUMULATIVE_PLOT
      //-------------------------------------------------------------------
      /* misc(s) */
      leg[idx] = pleg;
#ifdef  CUMULATIVE_PLOT
      if( idx != 0 )
	leg[idx].write = false;
#endif//CUMULATIVE_PLOT
      uni[idx] = puni;
      cap[idx] = basecap;
      //-------------------------------------------------------------------
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* individual file names */
#ifndef CUMULATIVE_PLOT
#ifdef  CONNECT_SYMBOLS_BY_LINE
  sprintf(figfile, "%s_%s", file, "breakdown_line");
#else///CONNECT_SYMBOLS_BY_LINE
  sprintf(figfile, "%s_%s", file, "breakdown");
#endif//CONNECT_SYMBOLS_BY_LINE
#else///CUMULATIVE_PLOT
#ifdef  CONNECT_SYMBOLS_BY_LINE
  sprintf(figfile, "%s_%s", file, "breakdown_cum_line");
#else///CONNECT_SYMBOLS_BY_LINE
  sprintf(figfile, "%s_%s", file, "breakdown_cum");
#endif//CONNECT_SYMBOLS_BY_LINE
#endif//CUMULATIVE_PLOT
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
  PLFLT colMap[NSEQUENCES][2][2];
  for(int ii = 0; ii < NSEQUENCES; ii++){
    colMap[ii][0][0] = box.xmin;
    colMap[ii][1][0] = box.xmax;
  }/* for(int ii = 0; ii < NSEQUENCES; ii++){ */
  colMap[0][0][1] = log10(1.5e+0);  colMap[0][1][1] = log10(3.2e+0);
  colMap[1][0][1] = log10(1.2e-1);  colMap[1][1][1] = log10(5.0e-1);
  colMap[2][0][1] = log10(3.0e-3);  colMap[2][1][1] = log10(1.6e-2);
  PLFLT pp[NSEQUENCES][2], c0[NSEQUENCES][2], c1[NSEQUENCES][2], c2[NSEQUENCES][2], aa[NSEQUENCES][2];
#if 1
  /* slow sequence: red (H = 0, L = 1/2, S = 1) */
  pp[0][0] = 0.0;  c0[0][0] = 0.0;  c1[0][0] = 0.5;  c2[0][0] = 1.0;  aa[0][0] = 0.4;
  pp[0][1] = 1.0;  c0[0][1] = 0.0;  c1[0][1] = 0.5;  c2[0][1] = 1.0;  aa[0][1] = 0.4;
  /* intermediate sequence: red (H = 0, L = 1/2, S = 1) */
  pp[1][0] = 0.0;  c0[1][0] = 0.0;  c1[1][0] = 0.5;  c2[1][0] = 1.0;  aa[1][0] = 0.2;
  pp[1][1] = 1.0;  c0[1][1] = 0.0;  c1[1][1] = 0.5;  c2[1][1] = 1.0;  aa[1][1] = 0.2;
  /* fast sequence: red (H = 0, L = 1/2, S = 1) */
  pp[2][0] = 0.0;  c0[2][0] = 0.0;  c1[2][0] = 0.5;  c2[2][0] = 1.0;  aa[2][0] = 0.1;
  pp[2][1] = 1.0;  c0[2][1] = 0.0;  c1[2][1] = 0.5;  c2[2][1] = 1.0;  aa[2][1] = 0.1;
#else
  /* slow sequence: cyan (H = 180, L = 1/2, S = 1) */
  pp[0][0] = 0.0;  c0[0][0] = 180.0;  c1[0][0] = 0.5;  c2[0][0] = 1.0;  aa[0][0] = 0.3;
  pp[0][1] = 1.0;  c0[0][1] = 180.0;  c1[0][1] = 0.5;  c2[0][1] = 1.0;  aa[0][1] = 0.3;
  /* intermediate sequence: green (H = 120, L = 1/2, S = 1) */
  pp[1][0] = 0.0;  c0[1][0] = 120.0;  c1[1][0] = 0.5;  c2[1][0] = 1.0;  aa[1][0] = 0.3;
  pp[1][1] = 1.0;  c0[1][1] = 120.0;  c1[1][1] = 0.5;  c2[1][1] = 1.0;  aa[1][1] = 0.3;
  /* fast sequence: yellow (H = 60, L = 1/2, S = 1) */
  pp[2][0] = 0.0;  c0[2][0] =  60.0;  c1[2][0] = 0.5;  c2[2][0] = 1.0;  aa[2][0] = 0.3;
  pp[2][1] = 1.0;  c0[2][1] =  60.0;  c1[2][1] = 0.5;  c2[2][1] = 1.0;  aa[2][1] = 0.3;
#endif
  PLBOOL rr[2] = {false, false};
  //-----------------------------------------------------------------------
#endif//HIGHLIGHTING_COLOR_MAP
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* create figure(s) */
  //-----------------------------------------------------------------------
  initPLplot(argc, argv, figfile);
  const PLBOOL xjoin = false;
  const PLBOOL yjoin = false;
  const PLFLT scale = 1.0;
  const PLFLT aspect = 0.25;
  //-----------------------------------------------------------------------
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      //-------------------------------------------------------------------
      PLINT idx = nypanel * ii + jj;
      if( plot[idx] ){
//-----------------------------------------------------------------
  	/* unified control of lines and points */
  	const PLBOOL useLS = false;
  	const PLBOOL usePT =  true;
  	//-----------------------------------------------------------------
  	setViewPort((idx == 0) ? true : false, nxpanel, nypanel, idx, xjoin, yjoin, range[idx], aspect);
  	//-----------------------------------------------------------------
#ifdef  HIGHLIGHTING_COLOR_MAP
#ifdef  CUMULATIVE_PLOT
	if( idx == 0 )
#endif//CUMULATIVE_PLOT
	for(int ss = 0; ss < NSEQUENCES; ss++){
	  plscmap1la(0, 2, pp[ss], c0[ss], c1[ss], c2[ss], aa[ss], rr);
	  plimagefr((const PLFLT * const *)colDat, 2, 2, colMap[ss][0][0], colMap[ss][1][0], colMap[ss][0][1], colMap[ss][1][1], 0.0, 1.0, 0.0, 1.0, NULL, NULL);
	}/* for(int ss = 0; ss < NSEQUENCES; ss++){ */
#endif//HIGHLIGHTING_COLOR_MAP
  	//-----------------------------------------------------------------
#if 1
	drawLine(1, lnum[idx], lx[idx], ly[idx],  line[idx], scale);
	drawLine(1, &lnum[idx][6], &lx[idx][6], &ly[idx][6],  &line[idx][6], scale);
#else
	drawLine(nlkind[idx], lnum[idx], lx[idx], ly[idx],  line[idx], scale);
#endif
  	dotPoint(npkind[idx], pnum[idx], px[idx], py[idx], point[idx], scale);
  	//-----------------------------------------------------------------
  	makeFrameAspect(nxpanel, ii, range[idx].xlog, range[idx].xgrd, xjoin, xlabel[idx], nypanel, jj, range[idx].ylog, range[idx].ygrd, yjoin, ylabel[idx], "", aspect);
#ifdef  HIGHLIGHTING_COLOR_MAP
#ifdef  CUMULATIVE_PLOT
	if( idx == 0 ){
#endif//CUMULATIVE_PLOT
	  plschr(0.0, cap[idx].scale * scale * sqrt(aspect) * 0.75);
	  plmtex("t", -1.8, 0.1, 0.0, "slow");
	  plmtex("t", -3.5, 0.1, 0.0, "intermediate");
	  plmtex("t", -8.0, 0.1, 0.0, "fast");
#ifdef  CUMULATIVE_PLOT
	}/* if( idx == 0 ){ */
#endif//CUMULATIVE_PLOT
#endif//HIGHLIGHTING_COLOR_MAP
  	//-----------------------------------------------------------------
  	if( cap[idx].write == true )	  makeCaption(cap[idx], scale);
  	//-----------------------------------------------------------------
  	if( leg[idx].write == true ){
	  PLFLT *bufPtSize;	  allocPLFLT(&bufPtSize, npkind[idx]);
	  for(int kk = 0; kk < npkind[idx]; kk++){
	    bufPtSize[kk] = point[idx][kk].scale;
	    point[idx][kk].scale = 0.8;
	  }/* for(int kk = 0; kk < npkind[idx]; kk++){ */
  	  makeLegend (idx, aspect, scale * 0.7, nxpanel, xjoin, nypanel, yjoin, uni[idx], leg[idx], useLS, nlkind[idx], line[idx], usePT, npkind[idx], point[idx]);
	  for(int kk = 0; kk < npkind[idx]; kk++)
	    point[idx][kk].scale = bufPtSize[kk];
	  free(bufPtSize);
	}/* if( leg[idx].write == true ){ */
  	//-----------------------------------------------------------------
      }/* if( plot[idx] ){ */
      //-------------------------------------------------------------------
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */
  //-----------------------------------------------------------------------
  exitPLplot();
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
#ifdef  HIGHLIGHTING_COLOR_MAP
  free(_colDat);
  free(colDat);
#endif//HIGHLIGHTING_COLOR_MAP
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
  free(pt);
#ifdef  CONNECT_SYMBOLS_BY_LINE
  free(ls);
#endif//CONNECT_SYMBOLS_BY_LINE
  //-----------------------------------------------------------------------
#ifdef  CONNECT_SYMBOLS_BY_LINE
  free(xtree);
  free(ymake);
  free(ysort);
  free(ysplt);
#ifdef  CUMULATIVE_PLOT
  free(ymakeSum);
  free(ysortSum);
  free(yspltSum);
#endif//CUMULATIVE_PLOT
#endif//CONNECT_SYMBOLS_BY_LINE
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
