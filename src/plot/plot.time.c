/**
 * @file plot.time.c
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

#define PLOT_SUMMARY_FIGURE_ONLY

#define COMPARE_FOCUS_VARIABLE

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "macro.h"
#include "myutil.h"
#include "mpilib.h"
#include "plplotlib.h"
#include "name.h"

#include "cdflib.h"


#define NMAX (18)

#define NARC (3)
/* #define NARC (2) */
#define NCHK (3)
#define NINT (2)
#define NMAC (4)
/* #define NMAC (3) */

static char arcname[NARC][ 8] =
  {"cc20", "cc35"
#   if  NARC > 2
   , "cc52"
#endif//NARC > 2
};
static char chkname[NCHK][ 8] = {"acc", "grv", "pot"};
static char macname[NMAC][16] =
  {"GADGET2_MAC", "Multipole_MAC", "Opening_angle"
#   if  NMAC > 3
   , "Bonsai2"
#endif//NMAC > 3
  };
static char synonym[NMAC][16] =
  {"acceleration", "multipole", "theta"
#   if  NMAC > 3
   , "bonsai2"
#endif//NMAC > 3
  };
static char intname[NINT][ 8] = {"block", "share"};


#define LOGPLOT_HOR
#define LOGPLOT_VER


void plotElapsedTime
(PLFLT err[restrict][NCHK][NARC][NINT][NMAC][NMAX], PLFLT time[restrict][NINT][NMAC][NMAX], PLplotPltRange box, char file[], MPIinfo mpi, int argc, char **argv);
void compareCheckTarget
(PLFLT err[restrict][NCHK][NARC][NINT][NMAC][NMAX], PLFLT time[restrict][NINT][NMAC][NMAX], PLplotPltRange box, char file[], MPIinfo mpi, int argc, char **argv);
void plotSpeedup
(PLFLT mac[restrict][NARC][NINT][NMAC][NMAX], PLFLT time[restrict][NINT][NMAC][NMAX], PLplotPltRange relativeMAC, PLplotPltRange openingAngle, char file[], int argc, char **argv);
void plotBenefit
(PLFLT mac[restrict][NARC][NINT][NMAC][NMAX], PLFLT time[restrict][NINT][NMAC][NMAX], PLplotPltRange relativeMAC, PLplotPltRange openingAngle, char file[], int argc, char **argv);


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);


  /** initialization */
  if( argc < 3 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 3);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -Ntot=<unsigned long int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 3 ){ */

  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "file", &file));
  ulong Ntot;  requiredCmdArg(getCmdArgUlg(argc, (const char * const *)argv, "Ntot", &Ntot));

  modifyArgcArgv4PLplot(&argc, argv, 3);


  /** prepare dataset */
  double *percentage;
  allocPercentile(&percentage);

  /** read the analyzed tree error for the initial snapshot */
  static double mac          [NCHK][NARC][NINT][NMAC][NMAX];
  static double err[NSUMMARY][NCHK][NARC][NINT][NMAC][NMAX];
  for(int jj = 0; jj < NCHK; jj++)
    for(int kk = 0; kk < NARC; kk++)
      for(int ll = 0; ll < NINT; ll++)
	for(int mm = 0; mm < NMAC; mm++)
	  for(int nn = 0; nn < NMAX; nn++)
	    mac[jj][kk][ll][mm][nn] = 0.0;
  for(int ii = 0; ii < NSUMMARY; ii++)
    for(int jj = 0; jj < NCHK; jj++)
      for(int kk = 0; kk < NARC; kk++)
	for(int ll = 0; ll < NINT; ll++)
	  for(int mm = 0; mm < NMAC; mm++)
	    for(int nn = 0; nn < NMAX; nn++)
	      err[ii][jj][kk][ll][mm][nn] = 0.0;
  for(int ii = 0; ii < NCHK; ii++)
    for(int jj = 0; jj < NARC; jj++)
      for(int kk = 0; kk < NMAC; kk++){
	if( ((((ii * NARC + jj) * NMAC) + kk) % mpi.size) == mpi.rank ){
	  FILE *fp;
	  char filename[128];
	  sprintf(filename, "%s/%s/%s.%s.%s.%.3d.txt", DATAFOLDER, arcname[jj], file, macname[kk], chkname[ii], 0);
	  fp = fopen(filename, "r");
	  /* if( fp == NULL ){	    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);	  } */
	  if( fp != NULL ){
	    int ll = 0;
	    while( EOF != fscanf(fp, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le",
				 &mac    [ii][jj][0][kk][ll],
				 &err[ 0][ii][jj][0][kk][ll], &err[ 1][ii][jj][0][kk][ll], &err[ 2][ii][jj][0][kk][ll], &err[ 3][ii][jj][0][kk][ll], &err[ 4][ii][jj][0][kk][ll],
				 &err[ 5][ii][jj][0][kk][ll], &err[ 6][ii][jj][0][kk][ll], &err[ 7][ii][jj][0][kk][ll], &err[ 8][ii][jj][0][kk][ll], &err[ 9][ii][jj][0][kk][ll],
				 &err[10][ii][jj][0][kk][ll], &err[11][ii][jj][0][kk][ll], &err[12][ii][jj][0][kk][ll], &err[13][ii][jj][0][kk][ll], &err[14][ii][jj][0][kk][ll],
				 &err[15][ii][jj][0][kk][ll], &err[16][ii][jj][0][kk][ll], &err[17][ii][jj][0][kk][ll], &err[18][ii][jj][0][kk][ll], &err[19][ii][jj][0][kk][ll],
				 &err[20][ii][jj][0][kk][ll], &err[21][ii][jj][0][kk][ll], &err[22][ii][jj][0][kk][ll], &err[23][ii][jj][0][kk][ll], &err[24][ii][jj][0][kk][ll],
				 &err[25][ii][jj][0][kk][ll], &err[26][ii][jj][0][kk][ll], &err[27][ii][jj][0][kk][ll])
		   )
	      ll++;
	    fclose(fp);

	    for(int mm = ll; mm < NMAX; mm++)
	      mac[ii][jj][0][kk][mm] = mac[ii][jj][0][kk][ll - 1];
	    for(int nn = 0; nn < NSUMMARY; nn++)
	      for(int mm = ll; mm < NMAX; mm++)
		err[nn][ii][jj][0][kk][mm] = err[nn][ii][jj][0][kk][ll - 1];

#ifdef  LOGPLOT_HOR
	    for(int nn = 0; nn < NSUMMARY; nn++)
	      for(int mm = 0; mm < NMAX; mm++)
		err[nn][ii][jj][0][kk][mm] = log10(err[nn][ii][jj][0][kk][mm]);
#endif//LOGPLOT_HOR

	    for(int oo = 0; oo < NINT; oo++)
	      for(int mm = 0; mm < NMAX; mm++)
		mac[ii][jj][oo][kk][mm] = mac[ii][jj][0][kk][mm];
	    for(int nn = 0; nn < NSUMMARY; nn++)
	      for(int oo = 0; oo < NINT; oo++)
		for(int mm = 0; mm < NMAX; mm++)
		  err[nn][ii][jj][oo][kk][mm] = err[nn][ii][jj][0][kk][mm];
	  }/* if( fp != NULL ){ */
	}/* if( ((((ii * NARC + jj) * NMAC) + kk) % mpi.size) == mpi.rank ){ */
      }/* for(int kk = 0; kk < NMAC; kk++){ */

  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE,  mac,            NCHK * NARC * NINT * NMAC * NMAX, MPI_DOUBLE, MPI_SUM, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE,  err, NSUMMARY * NCHK * NARC * NINT * NMAC * NMAX, MPI_DOUBLE, MPI_SUM, mpi.comm));

  /** read the elapsed time */
  static double time[NARC][NINT][NMAC][NMAX];
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < NINT; jj++)
      for(int kk = 0; kk < NMAC; kk++)
	for(int ll = 0; ll < NMAX; ll++)
	  time[ii][jj][kk][ll] = DBL_MAX;
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < NINT; jj++)
      for(int kk = 0; kk < NMAC; kk++){
	if( (((ii * NINT + jj) * NMAC + kk) % mpi.size) == mpi.rank ){
	  FILE *fp;
	  char filename[128];
	  sprintf(filename, "%s/%s.%s.%s.%s.%zu.time.log", LOGFOLDER, file, synonym[kk], intname[jj], arcname[ii], Ntot);
	  fp = fopen(filename, "r");
	  /* if( fp == NULL ){	    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);	  } */
	  if( fp != NULL ){
	    int ll = 0;
	    while( EOF != fscanf(fp, "%*le\t%*le\t%*le\t%le", &time[ii][jj][kk][ll] ) )
	      ll++;
	    fclose(fp);

#ifdef  LOGPLOT_VER
	    for(int mm = 0; mm < ll; mm++)
	      time[ii][jj][kk][mm] = log10(time[ii][jj][kk][mm]);
#endif//LOGPLOT_VER

	    for(int nn = 0; nn < NCHK; nn++)
	      for(int mm = ll; mm < NMAX; mm++)
		mac[nn][ii][jj][kk][mm] = mac[nn][ii][jj][kk][ll - 1];
	    for(int oo = 0; oo < NSUMMARY; oo++)
	      for(int nn = 0; nn < NCHK; nn++)
		for(int mm = ll; mm < NMAX; mm++)
		  err[oo][nn][ii][jj][kk][mm] = err[oo][nn][ii][jj][kk][ll - 1];
	  }/* if( fp != NULL ){ */
	}/* if( (((ii * NINT + jj) * NMAC + kk) % mpi.size) == mpi.rank ){ */
      }/* for(int kk = 0; kk < NMAC; kk++){ */

  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE,  mac,            NCHK * NARC * NINT * NMAC * NMAX, MPI_DOUBLE, MPI_MAX, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE,  err, NSUMMARY * NCHK * NARC * NINT * NMAC * NMAX, MPI_DOUBLE, MPI_MAX, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, time,                   NARC * NINT * NMAC * NMAX, MPI_DOUBLE, MPI_MIN, mpi.comm));


  /** plot elapsed time as a function of accuracy */
  PLplotPltRange range;
#ifdef  LOGPLOT_HOR
  range.xmin = log10(1.1e-6);
  range.xmax = log10(9.9e-1);
  range.xlog = LOGARITHMIC_PLOT;
#else//LOGPLOT_HOR
  range.xmin = 1.1e-6;
  range.xmax = 9.9e-1;
  range.xlog = LINEAR_PLOT;
#endif//LOGPLOT_HOR
#ifdef  LOGPLOT_VER
#   if  NARC > 2
  range.ymin = log10(2.0e-2);
  range.ymax = log10(2.0e+2);
#else///NARC > 2
  range.ymin = log10(1.0e-1);
  range.ymax = log10(9.9e+1);
#endif//NARC > 2
  range.ylog = LOGARITHMIC_PLOT;
#else///LOGPLOT_VER
  range.ymin =   0.0;
  range.ymax = 100.0;
  range.ylog = LINEAR_PLOT;
#endif//LOGPLOT_VER
  range.xgrd = range.ygrd = true;

  plotElapsedTime   (err, time, range, file, mpi, argc, argv);
#ifdef  COMPARE_FOCUS_VARIABLE
  compareCheckTarget(err, time, range, file, mpi, argc, argv);
#endif//COMPARE_FOCUS_VARIABLE


  /** plot speed-up as a function of accuracy controlling parameter */
  /** plot acceleration from Fermi architecture */
  PLplotPltRange relativeMAC, openingAngle;
  relativeMAC.xmin = log10(1.1e-6);
  relativeMAC.xmax = log10(9.9e-1);
  relativeMAC.xlog = LOGARITHMIC_PLOT;
#   if  NARC > 2
  relativeMAC.ymin = 0.4;
  relativeMAC.ymax = 4.9;
#else///NARC > 2
  relativeMAC.ymin = 0.9;
  relativeMAC.ymax = 1.7;
#endif//NARC > 2
  relativeMAC.ylog = LINEAR_PLOT;
  relativeMAC.xgrd = relativeMAC.ygrd = true;
  if( strcmp(file, "nfw") == 0 ){
#   if  NARC > 2
    relativeMAC.ymin = 0.0;
    relativeMAC.ymax = 4.0;
#else///NARC > 2
    relativeMAC.ymin = 0.8;
    relativeMAC.ymax = 2.1;
#endif//NARC > 2
  }/* if( strcmp(file, "nfw") == 0 ){ */
  openingAngle = relativeMAC;
  openingAngle.xmin = 0.0;
  openingAngle.xmax = 1.0;
  openingAngle.xlog = LINEAR_PLOT;
  /* for(int ii = 0; ii < NINT; ii++) */
  /*   for(int jj = 0; jj < NARC; jj++) */
  /*     for(int mm = 0; mm < NINT; mm++) */
  /* 	for(int kk = 0; kk < NMAC - 1; kk++) */
  /* 	  for(int ll = 0; ll < NMAX; ll++) */
  /* 	    mac[ii][jj][mm][kk][ll] = log10(mac[ii][jj][mm][kk][ll]); */
  for(int ii = 0; ii < NINT; ii++)
    for(int jj = 0; jj < NARC; jj++)
      for(int mm = 0; mm < NINT; mm++)
	for(int kk = 0; kk < ((NMAC <= 3) ? (NMAC) : (3)) - 1; kk++)
	  for(int ll = 0; ll < NMAX; ll++)
	    mac[ii][jj][mm][kk][ll] = log10(mac[ii][jj][mm][kk][ll]);

  if( mpi.rank == 0 )
    plotSpeedup(mac, time, relativeMAC, openingAngle, file, argc, argv);
  relativeMAC.ymin = openingAngle.ymin = 1.0;
  /* relativeMAC.ymax = openingAngle.ymax = 7.5; */
  relativeMAC.ymax = openingAngle.ymax = 6.5;
  if( strcmp(file, "nfw") == 0 ){
    relativeMAC.ymin = openingAngle.ymin = 1.0;
    relativeMAC.ymax = openingAngle.ymax = 4.5;
  }/* if( strcmp(file, "nfw") == 0 ){ */
  if( mpi.rank == mpi.size - 1 )
    plotBenefit(mac, time, relativeMAC, openingAngle, file, argc, argv);


  /** memory deallocation */
  freePercentile(percentage);


  exitMPI();
  return (0);
}


/** err [NSUMMARY][NCHK][NARC]      [NMAC][NMAX]; */
/** time                [NARC][NINT][NMAC][NMAX] */
void plotElapsedTime
(PLFLT err[restrict][NCHK][NARC][NINT][NMAC][NMAX], PLFLT time[restrict][NINT][NMAC][NMAX], PLplotPltRange box, char file[], MPIinfo mpi, int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /** global setting(s) */
  /** maximum number of data kinds */
#   if  NMAC <= 3
  const PLINT pkind = NINT * NMAC;
  const PLINT lkind = NINT * NMAC;
#else///NMAC <= 3
  const PLINT pkind = NINT * 3 + (NMAC - 3);
  const PLINT lkind = NINT * 3 + (NMAC - 3);
#endif//NMAC <= 3

  const PLBOOL puni = true;


  /** preparation for dataset */
  /** memory allocation for index */
  PLINT kind = (puni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  if( (pkind != kind) || (lkind != kind) ){
    __KILL__(stderr, "ERROR: pkind(%d), lkind(%d), and kind(%d) must be the same\n", pkind, lkind, kind);
  }/* if( (pkind != kind) || (lkind != kind) ){ */
  PLINT *num;  allocPLINT(&num, kind);

  /** set number of data points */
  for(PLINT ii = 0; ii < kind; ii++)    num[ii] = NMAX;

  /** set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  /** block time step */
  ls[0].color =   RED;  ls[0].style =  SOLID_LINE;  ls[0].width = MIDDLE_LINE;
  ls[1].color =  BLUE;  ls[1].style =  SOLID_LINE;  ls[1].width = MIDDLE_LINE;
  ls[2].color = BLACK;  ls[2].style =  SOLID_LINE;  ls[2].width = MIDDLE_LINE;
  sprintf(pt[0].type, PLplotSymbolType[fillCircle ]);  pt[0].scale = PLplotSymbolSize[fillCircle ];
  sprintf(pt[1].type, PLplotSymbolType[openSquare ]);  pt[1].scale = PLplotSymbolSize[openSquare ];
  sprintf(pt[2].type, PLplotSymbolType[fillDiamond]);  pt[2].scale = PLplotSymbolSize[fillDiamond];
  /** shared time step */
  ls[3].color =   RED;  ls[3].style = DOTTED_LINE;  ls[3].width = MIDDLE_LINE;
  ls[4].color =  BLUE;  ls[4].style = DOTTED_LINE;  ls[4].width = MIDDLE_LINE;
  ls[5].color = BLACK;  ls[5].style = DOTTED_LINE;  ls[5].width = MIDDLE_LINE;
  sprintf(pt[3].type, PLplotSymbolType[openCircle ]);  pt[3].scale = PLplotSymbolSize[openCircle ];
  sprintf(pt[4].type, PLplotSymbolType[fillSquare ]);  pt[4].scale = PLplotSymbolSize[fillSquare ];
  sprintf(pt[5].type, PLplotSymbolType[openDiamond]);  pt[5].scale = PLplotSymbolSize[openDiamond];
#   if  NMAC > 3
  /** Bonsai2 */
  ls[6].color = GREEN;  ls[6].style = DASHED_LINE;  ls[6].width = MIDDLE_LINE;
  sprintf(pt[6].type, PLplotSymbolType[fillTriangle]);  pt[6].scale = PLplotSymbolSize[fillTriangle];
#endif//NMAC > 3
  for(int ii = 0; ii < kind; ii++)
    pt[ii].color = ls[ii].color;

  /** set labels */
  char vlab[PLplotCharWords];  sprintf(vlab, "|#fi#<0x12>a#fr#di#u#utree#d - #fi#<0x12>a#fr#di#u#udirect#d| / #fia#fr#di#u#udirect#d");
  char slab[PLplotCharWords];  sprintf(slab, "|#fia#fr#di#u#utree#d / #fia#fr#di#u#udirect#d - 1|");
  char plab[PLplotCharWords];  sprintf(plab, "|#fi#gF#fr#di#u#utree#d / #fi#gF#fr#di#u#udirect#d -1|");
  char tlab[PLplotCharWords];  sprintf(tlab, "Elapsed time per step (s)");

  /** set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write =  true;
  sprintf(basecap.side, "%s", "b");
  PLplotCaption cc20cap;  setDefaultCaption(&cc20cap);  cc20cap.write = false;
  PLplotCaption cc35cap;  setDefaultCaption(&cc35cap);  cc35cap.write = false;
  PLplotCaption cc52cap;  setDefaultCaption(&cc52cap);  cc52cap.write = false;
  sprintf(cc20cap.side, "%s", "b");  sprintf(cc20cap.text, "%s", "M2090");
  sprintf(cc35cap.side, "%s", "b");  sprintf(cc35cap.text, "%s", "K20X");
  sprintf(cc52cap.side, "%s", "b");  sprintf(cc52cap.text, "%s", "GTX TITAN X");
  static PLplotCaption gpucap[NARC];
  gpucap[0] = cc20cap;
  gpucap[1] = cc35cap;
#   if  NARC > 2
  gpucap[2] = cc52cap;
#endif//NARC > 2
#   if  NSUMMARY == 28
  static PLplotCaption pcap[NSUMMARY];
  for(int ii = 0; ii < NSUMMARY; ii++){
    setDefaultCaption(&pcap[ii]);
    pcap[ii].write = false;
    sprintf(pcap[ii].side, "%s", "b");
  }/* for(int ii = 0; ii < NSUMMARY; ii++){ */
  sprintf(pcap[ 0].text, "%s", "100% Error");
  sprintf(pcap[ 1].text, "%s", "+5#fi#gs#fr Error");
  sprintf(pcap[ 2].text, "%s", "+4#fi#gs#fr Error");
  sprintf(pcap[ 3].text, "%s", "+3#fi#gs#fr Error");
  sprintf(pcap[ 4].text, "%s", "99% Error");
  sprintf(pcap[ 5].text, "%s", "+2#fi#gs#fr Error");
  sprintf(pcap[ 6].text, "%s", "95% Error");
  sprintf(pcap[ 7].text, "%s", "90% Error");
  sprintf(pcap[ 8].text, "%s", "80% Error");
  sprintf(pcap[ 9].text, "%s", "75% Error");
  sprintf(pcap[10].text, "%s", "70% Error");
  sprintf(pcap[11].text, "%s", "+1#fi#gs#fr Error");
  sprintf(pcap[12].text, "%s", "67% Error");
  sprintf(pcap[13].text, "%s", "60% Error");
  sprintf(pcap[14].text, "%s", "Median Error");
  sprintf(pcap[15].text, "%s", "40% Error");
  sprintf(pcap[16].text, "%s", "33% Error");
  sprintf(pcap[17].text, "%s", "-1#fi#gs#fr Error");
  sprintf(pcap[18].text, "%s", "30% Error");
  sprintf(pcap[19].text, "%s", "25% Error");
  sprintf(pcap[20].text, "%s", "20% Error");
  sprintf(pcap[21].text, "%s", "10% Error");
  sprintf(pcap[22].text, "%s", "5% Error");
  sprintf(pcap[23].text, "%s", "-2#fi#gs#fr Error");
  sprintf(pcap[24].text, "%s", "1%y Error");
  sprintf(pcap[25].text, "%s", "-3#fi#gs#fr Error");
  sprintf(pcap[26].text, "%s", "-4#fi#gs#fr Error");
  sprintf(pcap[27].text, "%s", "-5#fi#gs#fr Error");
#endif//NSUMMARY == 28

  /** set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = true;
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
  static char tmpname[3][16] = {"acc", "mul", "#gh"};
#   if  NMAC <= 3
  for(int ii = 0; ii < kind; ii++)
    sprintf(pleg.text[ii], "%s, %s", intname[ii / NMAC], tmpname[ii % NMAC]);
#else///NMAC <= 3
  for(int ii = 0; ii < NINT * 3; ii++)
    sprintf(pleg.text[ii], "%s, %s", intname[ii / 3], tmpname[ii % 3]);
  sprintf(pleg.text[6], "%s", "Bonsai");
#endif//NMAC <= 3
  pleg.pos = PL_POSITION_RIGHT | PL_POSITION_TOP | PL_POSITION_INSIDE;


  /** memory allocation to exploit multiple-plot function */
  /** maximum number of panels in each direction */
  const PLINT nxpanel = NARC;
  const PLINT nypanel = 2;

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


  /** configure to plot elapsed time as a function of tree error by multiple processes */
  for(int prec = mpi.rank; prec < NSUMMARY >> 1; prec += mpi.size){
    /** data preparation */
    for(int ii = 0; ii < nxpanel; ii++){
#   if  NMAC <= 3
      for(int jj = 0; jj < kind; jj++){
	hor[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &err [NSUMMARY >> 1][0][ii][jj / NMAC][jj % NMAC][0];/**< median error */
	hor[INDEX(nxpanel, nypanel, kind, ii, 1, jj)] = &err [         prec][0][ii][jj / NMAC][jj % NMAC][0];
	ver[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &time                  [ii][jj / NMAC][jj % NMAC][0];
	ver[INDEX(nxpanel, nypanel, kind, ii, 1, jj)] = &time                  [ii][jj / NMAC][jj % NMAC][0];
      }/* for(int jj = 0; jj < kind; jj++){ */
#else///NMAC <= 3
      /* GOTHIC */
      for(int jj = 0; jj < NINT * 3; jj++){
	hor[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &err [NSUMMARY >> 1][0][ii][jj / 3][jj % 3][0];/**< median error */
	hor[INDEX(nxpanel, nypanel, kind, ii, 1, jj)] = &err [         prec][0][ii][jj / 3][jj % 3][0];
	ver[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &time                  [ii][jj / 3][jj % 3][0];
	ver[INDEX(nxpanel, nypanel, kind, ii, 1, jj)] = &time                  [ii][jj / 3][jj % 3][0];
      }/* for(int jj = 0; jj < kind; jj++){ */
      /* Bonsai2 */
      hor[INDEX(nxpanel, nypanel, kind, ii, 0, 6)] = &err [NSUMMARY >> 1][0][ii][1][3][0];/**< median error */
      hor[INDEX(nxpanel, nypanel, kind, ii, 1, 6)] = &err [         prec][0][ii][1][3][0];
      ver[INDEX(nxpanel, nypanel, kind, ii, 0, 6)] = &time                  [ii][1][3][0];
      ver[INDEX(nxpanel, nypanel, kind, ii, 1, 6)] = &time                  [ii][1][3][0];
#endif//NMAC <= 3
    }/* for(int ii = 0; ii < nxpanel; ii++){ */

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
	xlabel[idx] = vlab;
	ylabel[idx] = tlab;

	/** misc(s) */
	leg[idx] = pleg;
	uni[idx] = puni;
	cap[idx] = basecap;
	/** 97 is "a" in ASCII code */
	sprintf(cap[idx].text, "(%c) %s on %s", 97 + INDEX2D(nypanel, nxpanel, nypanel - jj - 1, ii), pcap[(jj == 0) ? (NSUMMARY >> 1) : prec].text, gpucap[ii].text);
      }/* for(PLINT jj = 0; jj < nypanel; jj++){ */


    /** individual file names */
    sprintf(figfile, "%s_%s_%.2d", file, "acc", prec);
    for(int ii = 0; ii < nxpanel; ii++)
      sprintf(figname[ii], "%s_%s_%.2d_%s", file, "acc", prec, arcname[ii]);


    /** create figure(s) */
#ifndef PLOT_SUMMARY_FIGURE_ONLY
    for(int ii = 0; ii < nxpanel; ii++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 1);
      const PLBOOL tmp = cap[idx].write;
      cap[idx].write = false;

      plotData(1, 1, &plot[idx], false, false,
	       &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
	       &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
	       &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
	       &cap[idx], &leg[idx], &uni[idx], &range[idx],
	       &xlabel[idx], &ylabel[idx], "", figname[ii], argc, argv);

      cap[idx].write = tmp;
    }/* for(int ii = 0; ii < nxpanel; ii++){ */
#endif//PLOT_SUMMARY_FIGURE_ONLY

    for(int ii = 0; ii < nxpanel; ii++)
      for(int jj = 0; jj < nypanel; jj++){
	const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
	/* if( idx != INDEX2D(nxpanel, nypanel, 0, nypanel - 1) ) */
	if( idx != INDEX2D(nxpanel, nypanel, 0, 0) )
	  leg[idx].write = false;
      }/* for(int jj = 0; jj < nypanel; jj++){ */
    plotData(nxpanel, nypanel, plot, true, true,
	     nlkind,  line, lnum, lx, ly,
	     npkind, point, pnum, px, py,
	     errbar, NULL, NULL, NULL, NULL, NULL, NULL,
	     cap, leg, uni, range,
	     xlabel, ylabel, "", figfile, argc, argv);
  }/* for(int prec = mpi.rank; prec < NSUMMARY >> 1; prec += mpi.size){ */


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


/** err [NSUMMARY][NCHK][NARC]      [NMAC][NMAX]; */
/** time                [NARC][NINT][NMAC][NMAX] */
void compareCheckTarget
(PLFLT err[restrict][NCHK][NARC][NINT][NMAC][NMAX], PLFLT time[restrict][NINT][NMAC][NMAX], PLplotPltRange box, char file[], MPIinfo mpi, int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /** global setting(s) */
  /** maximum number of data kinds */
  const PLINT pkind = NINT * ((NMAC <= 3) ? NMAC : 3);
  const PLINT lkind = NINT * ((NMAC <= 3) ? NMAC : 3);

  const PLBOOL puni = true;


  /** preparation for dataset */
  /** memory allocation for index */
  PLINT kind = (puni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  if( (pkind != kind) || (lkind != kind) ){
    __KILL__(stderr, "ERROR: pkind(%d), lkind(%d), and kind(%d) must be the same\n", pkind, lkind, kind);
  }/* if( (pkind != kind) || (lkind != kind) ){ */
  PLINT *num;  allocPLINT(&num, kind);

  /** set number of data points */
  for(PLINT ii = 0; ii < kind; ii++)    num[ii] = NMAX;

  /** set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  /** block time step */
  ls[0].color =   RED;  ls[0].style =  SOLID_LINE;  ls[0].width =   BOLD_LINE;
  ls[1].color =  BLUE;  ls[1].style =  SOLID_LINE;  ls[1].width = MIDDLE_LINE;
  ls[2].color = BLACK;  ls[2].style =  SOLID_LINE;  ls[2].width = MIDDLE_LINE;
  sprintf(pt[0].type, PLplotSymbolType[fillCircle ]);  pt[0].scale = PLplotSymbolSize[fillCircle ] * 0.75;
  sprintf(pt[1].type, PLplotSymbolType[openSquare ]);  pt[1].scale = PLplotSymbolSize[openSquare ] * 0.75;
  sprintf(pt[2].type, PLplotSymbolType[fillDiamond]);  pt[2].scale = PLplotSymbolSize[fillDiamond] * 0.75;
  /** shared time step */
  ls[3].color =   RED;  ls[3].style = DOTTED_LINE;  ls[3].width =   BOLD_LINE;
  ls[4].color =  BLUE;  ls[4].style = DOTTED_LINE;  ls[4].width = MIDDLE_LINE;
  ls[5].color = BLACK;  ls[5].style = DOTTED_LINE;  ls[5].width = MIDDLE_LINE;
  sprintf(pt[3].type, PLplotSymbolType[openCircle ]);  pt[3].scale = PLplotSymbolSize[openCircle ] * 0.75;
  sprintf(pt[4].type, PLplotSymbolType[fillSquare ]);  pt[4].scale = PLplotSymbolSize[fillSquare ] * 0.75;
  sprintf(pt[5].type, PLplotSymbolType[openDiamond]);  pt[5].scale = PLplotSymbolSize[openDiamond] * 0.75;
  for(int ii = 0; ii < kind; ii++)
    pt[ii].color = ls[ii].color;

  /** set labels */
  char vlab[PLplotCharWords];  sprintf(vlab, "|#fi#<0x12>a#fr#di#u#utree#d - #fi#<0x12>a#fr#di#u#udirect#d| / #fia#fr#di#u#udirect#d");
  char slab[PLplotCharWords];  sprintf(slab, "|#fia#fr#di#u#utree#d / #fia#fr#di#u#udirect#d - 1|");
  char plab[PLplotCharWords];  sprintf(plab, "|#fi#gF#fr#di#u#utree#d / #fi#gF#fr#di#u#udirect#d -1|");
  char tlab[PLplotCharWords];  sprintf(tlab, "Elapsed time per step (s)");

  /** set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write =  true;
  sprintf(basecap.side, "%s", "b");
  PLplotCaption cc20cap;  setDefaultCaption(&cc20cap);  cc20cap.write = false;
  PLplotCaption cc35cap;  setDefaultCaption(&cc35cap);  cc35cap.write = false;
  PLplotCaption cc52cap;  setDefaultCaption(&cc52cap);  cc52cap.write = false;
  sprintf(cc20cap.side, "%s", "b");  sprintf(cc20cap.text, "%s", "M2090");
  sprintf(cc35cap.side, "%s", "b");  sprintf(cc35cap.text, "%s", "K20X");
  sprintf(cc52cap.side, "%s", "b");  sprintf(cc52cap.text, "%s", "GTX TITAN X");
#   if  NSUMMARY == 28
  static PLplotCaption pcap[NSUMMARY];
  for(int ii = 0; ii < NSUMMARY; ii++){
    setDefaultCaption(&pcap[ii]);
    pcap[ii].write = false;
    sprintf(pcap[ii].side, "%s", "b");
  }/* for(int ii = 0; ii < NSUMMARY; ii++){ */
  sprintf(pcap[ 0].text, "%s", "100% Error");
  sprintf(pcap[ 1].text, "%s", "+5#fi#gs#fr Error");
  sprintf(pcap[ 2].text, "%s", "+4#fi#gs#fr Error");
  sprintf(pcap[ 3].text, "%s", "+3#fi#gs#fr Error");
  sprintf(pcap[ 4].text, "%s", "99% Error");
  sprintf(pcap[ 5].text, "%s", "+2#fi#gs#fr Error");
  sprintf(pcap[ 6].text, "%s", "95% Error");
  sprintf(pcap[ 7].text, "%s", "90% Error");
  sprintf(pcap[ 8].text, "%s", "80% Error");
  sprintf(pcap[ 9].text, "%s", "75% Error");
  sprintf(pcap[10].text, "%s", "70% Error");
  sprintf(pcap[11].text, "%s", "+1#fi#gs#fr Error");
  sprintf(pcap[12].text, "%s", "67% Error");
  sprintf(pcap[13].text, "%s", "60% Error");
  sprintf(pcap[14].text, "%s", "Median Error");
  sprintf(pcap[15].text, "%s", "40% Error");
  sprintf(pcap[16].text, "%s", "33% Error");
  sprintf(pcap[17].text, "%s", "-1#fi#gs#fr Error");
  sprintf(pcap[18].text, "%s", "30% Error");
  sprintf(pcap[19].text, "%s", "25% Error");
  sprintf(pcap[20].text, "%s", "20% Error");
  sprintf(pcap[21].text, "%s", "10% Error");
  sprintf(pcap[22].text, "%s", "5% Error");
  sprintf(pcap[23].text, "%s", "-2#fi#gs#fr Error");
  sprintf(pcap[24].text, "%s", "1%y Error");
  sprintf(pcap[25].text, "%s", "-3#fi#gs#fr Error");
  sprintf(pcap[26].text, "%s", "-4#fi#gs#fr Error");
  sprintf(pcap[27].text, "%s", "-5#fi#gs#fr Error");
#endif//NSUMMARY == 28

  /** set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = true;
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
  static char tmpname[3][16] = {"acc", "mul", "#gh"};
  for(int ii = 0; ii < kind; ii++)
    sprintf(pleg.text[ii], "%s, %s", intname[ii / ((NMAC <= 3) ? NMAC : 3)], tmpname[ii % ((NMAC <= 3) ? NMAC : 3)]);
  pleg.pos = PL_POSITION_RIGHT | PL_POSITION_TOP | PL_POSITION_INSIDE;


  /** memory allocation to exploit multiple-plot function */
  /** maximum number of panels in each direction */
  const PLINT nxpanel = NCHK;
  const PLINT nypanel = 2;

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


  /** configure to plot elapsed time as a function of tree error by multiple processes */
  for(int prec = mpi.rank; prec < NSUMMARY >> 1; prec += mpi.size){
    /** data preparation */
    for(int ii = 0; ii < nxpanel; ii++)
      for(int jj = 0; jj < kind; jj++){
	hor[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &err [NSUMMARY >> 1][ii][0][jj / ((NMAC <= 3) ? NMAC : 3)][jj % ((NMAC <= 3) ? NMAC : 3)][0];/**< median error */
	hor[INDEX(nxpanel, nypanel, kind, ii, 1, jj)] = &err [         prec][ii][0][jj / ((NMAC <= 3) ? NMAC : 3)][jj % ((NMAC <= 3) ? NMAC : 3)][0];
	ver[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &time                   [0][jj / ((NMAC <= 3) ? NMAC : 3)][jj % ((NMAC <= 3) ? NMAC : 3)][0];
	ver[INDEX(nxpanel, nypanel, kind, ii, 1, jj)] = &time                   [0][jj / ((NMAC <= 3) ? NMAC : 3)][jj % ((NMAC <= 3) ? NMAC : 3)][0];
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
	switch( ii ){
	case 0:	  xlabel[idx] = vlab;	  break;
	case 1:	  xlabel[idx] = slab;	  break;
	case 2:	  xlabel[idx] = plab;	  break;
	default:	  __KILL__(stderr, "ERROR: ii(%d) must be 0, 1, or 2.\n", ii);	  break;
	}/* switch( ii ){ */
	ylabel[idx] = tlab;

	/** misc(s) */
	leg[idx] = pleg;
	uni[idx] = puni;
	cap[idx] = basecap;
	/** 97 is "a" in ASCII code */
	sprintf(cap[idx].text, "(%c) %s about %s", 97 + INDEX2D(nypanel, nxpanel, nypanel - jj - 1, ii), pcap[(jj == 0) ? (NSUMMARY >> 1) : prec].text, chkname[ii]);
      }/* for(PLINT jj = 0; jj < nypanel; jj++){ */


    /** individual file names */
    sprintf(figfile, "%s_%s_%.2d", file, arcname[0], prec);
    for(int ii = 0; ii < nxpanel; ii++)
      sprintf(figname[ii], "%s_%s_%.2d_%s", file, arcname[0], prec, chkname[ii]);


    /** create figure(s) */
#ifndef PLOT_SUMMARY_FIGURE_ONLY
    for(int ii = 0; ii < nxpanel; ii++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 1);
      const PLBOOL tmp = cap[idx].write;
      cap[idx].write = false;

      plotData(1, 1, &plot[idx], false, false,
	       &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
	       &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
	       &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
	       &cap[idx], &leg[idx], &uni[idx], &range[idx],
	       &xlabel[idx], &ylabel[idx], "", figname[ii], argc, argv);

      cap[idx].write = tmp;
    }/* for(int ii = 0; ii < nxpanel; ii++){ */
#endif//PLOT_SUMMARY_FIGURE_ONLY

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
  }/* for(int prec = mpi.rank; prec < NSUMMARY >> 1; prec += mpi.size){ */


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


/** mac [NCHK][NARC][NINT][NMAC][NMAX]; */
/** time      [NARC][NINT][NMAC][NMAX] */
void plotSpeedup
(PLFLT mac[restrict][NARC][NINT][NMAC][NMAX], PLFLT time[restrict][NINT][NMAC][NMAX], PLplotPltRange relativeMAC, PLplotPltRange openingAngle, char file[], int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /** global setting(s) */
  /** maximum number of data kinds */
  const PLINT pkind = NINT * (NARC - 1);
  const PLINT lkind = NINT * (NARC - 1);

  const PLBOOL puni = true;


  /** preparation for dataset */
  /** memory allocation for index */
  PLINT kind = (puni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  if( (pkind != kind) || (lkind != kind) ){
    __KILL__(stderr, "ERROR: pkind(%d), lkind(%d), and kind(%d) must be the same\n", pkind, lkind, kind);
  }/* if( (pkind != kind) || (lkind != kind) ){ */
  PLINT *num;  allocPLINT(&num, kind);

  /** set number of data points */
  for(PLINT ii = 0; ii < kind; ii++)    num[ii] = NMAX;

  /** data preparation and analysis */
  static PLFLT ss[((NMAC <= 3) ? NMAC : 3)][NINT][NARC - 1][NMAX];
#ifdef  LOGPLOT_VER
  for(int ii = 0; ii < NINT; ii++)
    for(int jj = 0; jj < ((NMAC <= 3) ? NMAC : 3); jj++)
      for(int kk = 0; kk < NMAX; kk++)
	for(int ll = 0; ll < NARC - 1; ll++)
	  ss[jj][ii][ll][kk] = pow(10.0, time[0][ii][jj][kk] - time[1 + ll][ii][jj][kk]);
#else///LOGPLOT_VER
  for(int ii = 0; ii < NINT; ii++)
    for(int jj = 0; jj < ((NMAC <= 3) ? NMAC : 3); jj++)
      for(int kk = 0; kk < NMAX; kk++){
	const PLFLT inv = 1.0 / time[0][ii][jj][kk];
	for(int ll = 0; ll < NARC - 1; ll++)
	  ss[jj][ii][ll][kk] = time[1 + ll][ii][jj][kk] * inv;
      }/* for(int kk = 0; kk < NMAX; kk++){ */
#endif//LOGPLOT_VER
  /** mac [NCHK][NARC][NINT][((NMAC <= 3) ? NMAC : 3)][NMAX]; */
  for(int ii = 0; ii < ((NMAC <= 3) ? NMAC : 3); ii++)
    for(int jj = 0; jj < NINT; jj++)
      for(int kk = 1; kk < NMAX; kk++)
	if( mac[0][0][jj][ii][kk] >= mac[0][0][jj][ii][kk - 1] )
	  for(int ll = 0; ll < NARC - 1; ll++)
	    ss[ii][jj][ll][kk] = DBL_MAX;

  /** set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  /** block time step on K20X */
  ls[0].color = BLACK;  ls[0].style =  SOLID_LINE;  ls[0].width = MIDDLE_LINE;
  sprintf(pt[0].type, PLplotSymbolType[fillCircle  ]);  pt[0].scale = PLplotSymbolSize[fillCircle  ] * 0.75;
  /** shared time step on K20X */
  ls[1].color = BLACK;  ls[1].style = DOTTED_LINE;  ls[1].width = MIDDLE_LINE;
  sprintf(pt[1].type, PLplotSymbolType[openCircle  ]);  pt[1].scale = PLplotSymbolSize[openCircle  ] * 0.75;
#   if  NARC > 2
  /** block time step on TITAN X */
  ls[2].color = BLACK;  ls[2].style =  SOLID_LINE;  ls[2].width = MIDDLE_LINE;
  sprintf(pt[2].type, PLplotSymbolType[fillSquare  ]);  pt[2].scale = PLplotSymbolSize[fillSquare  ] * 0.75;
  /** shared time step on TITAN X */
  ls[3].color = BLACK;  ls[3].style = DOTTED_LINE;  ls[3].width = MIDDLE_LINE;
  sprintf(pt[3].type, PLplotSymbolType[openSquare  ]);  pt[3].scale = PLplotSymbolSize[openSquare  ] * 0.75;
#endif//NARC > 2
  for(int ii = 0; ii < kind; ii++)
    pt[ii].color = ls[ii].color;

  /** set labels */
  /** char hlab[PLplotCharWords];  sprintf(hlab, "Accuracy controlling parameter"); */
  char alab[PLplotCharWords];  sprintf(alab, "#fi#gD#fr#dacc#u");
  char mlab[PLplotCharWords];  sprintf(mlab, "#fi#gD#fr#dmul#u");
  char tlab[PLplotCharWords];  sprintf(tlab, "#gh");
  char vlab[PLplotCharWords];  sprintf(vlab, "Speed up");

  /** set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write =  true;
  sprintf(basecap.side, "%s", "b");
  PLplotCaption   acccap;  setDefaultCaption(&  acccap);    acccap.write = false;
  PLplotCaption   mulcap;  setDefaultCaption(&  mulcap);    mulcap.write = false;
  PLplotCaption thetacap;  setDefaultCaption(&thetacap);  thetacap.write = false;
  sprintf(  acccap.side, "%s", "b");  sprintf(  acccap.text, "%s", "Acceleration MAC");
  sprintf(  mulcap.side, "%s", "b");  sprintf(  mulcap.text, "%s", "Multipole MAC");
  sprintf(thetacap.side, "%s", "b");  sprintf(thetacap.text, "%s", "Opening angle");
  static PLplotCaption maccap[((NMAC <= 3) ? NMAC : 3)];
  maccap[0] =   acccap;
  maccap[1] =   mulcap;
  maccap[2] = thetacap;

  /** set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = true;
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
  static char tmpname[NARC - 1][16] = {"K20X"
#   if  NARC > 2
				       , "GTX TITAN X"
#endif//NARC > 2
  };
  /* for(int ii = 0; ii < kind; ii++) */
  /*   sprintf(pleg.text[ii], "%s on %s", intname[ii / (NARC - 1)], tmpname[ii % (NARC - 1)]); */
  for(int ii = 0; ii < kind; ii++)
    sprintf(pleg.text[ii], "%s on %s", intname[ii % NINT], tmpname[ii / NINT]);
  pleg.pos = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;


  /** memory allocation to exploit multiple-plot function */
  /** maximum number of panels in each direction */
  const PLINT nxpanel = ((NMAC <= 3) ? NMAC : 3);
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


  /** configure to plot elapsed time as a function of tree error by multiple processes */
  /** data preparation */
  /** mac[NCHK][NARC][NINT][((NMAC <= 3) ? NMAC : 3)][NMAX]; */
  /** ss       [((NMAC <= 3) ? NMAC : 3)][NINT][NARC - 1][NMAX]; */
  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < kind; jj++){
      /* hor[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &mac[ 0][              0][jj / (NARC - 1)][ii             ][0]; */
      /* ver[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &ss [ii][jj / (NARC - 1)][jj % (NARC - 1)][0]; */
      hor[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &mac[ 0][0][jj % NINT][ii][0];
      ver[INDEX(nxpanel, nypanel, kind, ii, 0, jj)] = &ss [ii]   [jj % NINT][jj / NINT][0];
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
      range[idx] = (ii < 2) ? relativeMAC : openingAngle;

      /** label setting(s) */
      switch( ii ){
      case 0:	xlabel[idx] = alab;	break;
      case 1:	xlabel[idx] = mlab;	break;
      case 2:	xlabel[idx] = tlab;	break;
      default:
	__KILL__(stderr, "ERROR: ii(%d) must be smaller than 2\n", ii);
	break;
      }/* switch( ii ){ */
      ylabel[idx] = vlab;

      /** misc(s) */
      leg[idx] = pleg;
      uni[idx] = puni;
      cap[idx] = basecap;
      /* 97 is "a" in ASCII code */
      sprintf(cap[idx].text, "(%c) %s", 97 + INDEX2D(nypanel, nxpanel, nypanel - jj - 1, ii), maccap[ii].text);
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */


  /** individual file names */
  sprintf(figfile, "%s_%s", file, "speedup");
  for(int ii = 0; ii < nxpanel; ii++)
    sprintf(figname[ii], "%s_%s_%s", file, "speedup", synonym[ii]);


  /** create figure(s) */
/* #ifndef PLOT_SUMMARY_FIGURE_ONLY */
  for(int ii = 0; ii < nxpanel; ii++){
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 0);
    const PLBOOL tmp = cap[idx].write;
    cap[idx].write = false;

    plotData(1, 1, &plot[idx], false, false,
	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
	     &xlabel[idx], &ylabel[idx], "", figname[ii], argc, argv);
    cap[idx].write = tmp;
  }/* for(int ii = 0; ii < nxpanel; ii++){ */
/* #endif//PLOT_SUMMARY_FIGURE_ONLY */

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


/** mac [NCHK][NARC][NINT][NMAC][NMAX]; */
/** time      [NARC][NINT][NMAC][NMAX] */
void plotBenefit
(PLFLT mac[restrict][NARC][NINT][NMAC][NMAX], PLFLT time[restrict][NINT][NMAC][NMAX], PLplotPltRange relativeMAC, PLplotPltRange openingAngle, char file[], int argc, char **argv)
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
  for(PLINT ii = 0; ii < kind; ii++)    num[ii] = NMAX;

  /** data preparation and analysis */
  static PLFLT gain[NARC][((NMAC <= 3) ? NMAC : 3)][NMAX];
#ifdef  LOGPLOT_VER
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < ((NMAC <= 3) ? NMAC : 3); jj++)
      for(int kk = 0; kk < NMAX; kk++)
	gain[ii][jj][kk] = pow(10.0, time[ii][1][jj][kk] - time[ii][0][jj][kk]);
#else///LOGPLOT_VER
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < ((NMAC <= 3) ? NMAC : 3); jj++)
      for(int kk = 0; kk < NMAX; kk++)
	gain[ii][jj][kk] = time[ii][0][jj][kk] / time[ii][1][jj][kk];
#endif//LOGPLOT_VER
  for(int ii = 0; ii < NARC; ii++)
    for(int jj = 0; jj < ((NMAC <= 3) ? NMAC : 3); jj++)
      for(int kk = 1; kk < NMAX; kk++)
	if( mac[0][ii][1][jj][kk] >= mac[0][ii][1][jj][kk - 1] )
	  gain[ii][jj][kk] = DBL_MAX;

  /** set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  ls[0].color = BLACK;  ls[0].style = DASHED_LINE;  ls[0].width = MIDDLE_LINE;
  ls[1].color = BLACK;  ls[1].style =  SOLID_LINE;  ls[1].width = MIDDLE_LINE;
#   if  NARC > 2
  ls[2].color = BLACK;  ls[2].style = DOTTED_LINE;  ls[2].width = MIDDLE_LINE;
#endif//NARC > 2
  sprintf(pt[0].type, PLplotSymbolType[openCircle ]);  pt[0].scale = PLplotSymbolSize[openCircle ] * 0.75;
  sprintf(pt[1].type, PLplotSymbolType[openSquare ]);  pt[1].scale = PLplotSymbolSize[openSquare ] * 0.75;
#   if  NARC > 2
  sprintf(pt[2].type, PLplotSymbolType[fillDiamond]);  pt[2].scale = PLplotSymbolSize[fillDiamond] * 0.75;
#endif//NARC > 2
  for(int ii = 0; ii < kind; ii++)
    pt[ii].color = ls[ii].color;

  /** set labels */
  char alab[PLplotCharWords];  sprintf(alab, "#fi#gD#fr#dacc#u");
  char mlab[PLplotCharWords];  sprintf(mlab, "#fi#gD#fr#dmul#u");
  char tlab[PLplotCharWords];  sprintf(tlab, "#gh");
  char vlab[PLplotCharWords];  sprintf(vlab, "Speed up");

  /** set caption(s) */
  PLplotCaption basecap;  setDefaultCaption(&basecap);  basecap.write =  true;
  sprintf(basecap.side, "%s", "b");
  PLplotCaption   acccap;  setDefaultCaption(&  acccap);    acccap.write = false;
  PLplotCaption   mulcap;  setDefaultCaption(&  mulcap);    mulcap.write = false;
  PLplotCaption thetacap;  setDefaultCaption(&thetacap);  thetacap.write = false;
  sprintf(  acccap.side, "%s", "b");  sprintf(  acccap.text, "%s", "Acceleration MAC");
  sprintf(  mulcap.side, "%s", "b");  sprintf(  mulcap.text, "%s", "Multipole MAC");
  sprintf(thetacap.side, "%s", "b");  sprintf(thetacap.text, "%s", "Opening angle");
  static PLplotCaption maccap[((NMAC <= 3) ? NMAC : 3)];
  maccap[0] =   acccap;
  maccap[1] =   mulcap;
  maccap[2] = thetacap;

  /** set legends */
  PLplotLegend pleg;  setDefaultLegend(&pleg, false);  pleg.write = true;
  char *plegTxt;
  allocChar4PLplot(&plegTxt, kind);
  allocPointer4Char4PLplot(&(pleg.text), kind);
  assignChar4PLplot(kind, pleg.text, plegTxt);
  sprintf(pleg.text[0], "%s", "M2090");
  sprintf(pleg.text[1], "%s", "K20X");
#   if  NARC > 2
  sprintf(pleg.text[2], "%s", "GTX TITAN X");
#endif//NARC > 2
  pleg.pos = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;


  /** memory allocation to exploit multiple-plot function */
  /** maximum number of panels in each direction */
  const PLINT nxpanel = ((NMAC <= 3) ? NMAC : 3);
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


  /** configure to plot elapsed time as a function of tree error by multiple processes */
  /** data preparation */
  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < kind; jj++){
      hor[INDEX2D(nxpanel, kind, ii, jj)] = &mac [0][jj][0][ii][0];
      ver[INDEX2D(nxpanel, kind, ii, jj)] = &gain   [jj]   [ii][0];
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
      range[idx] = (ii < 2) ? relativeMAC : openingAngle;

      /** label setting(s) */
      switch( ii ){
      case 0:	xlabel[idx] = alab;	break;
      case 1:	xlabel[idx] = mlab;	break;
      case 2:	xlabel[idx] = tlab;	break;
      default:
	__KILL__(stderr, "ERROR: ii(%d) must be smaller than 2\n", ii);
	break;
      }/* switch( ii ){ */
      ylabel[idx] = vlab;

      /** misc(s) */
      leg[idx] = pleg;
      uni[idx] = puni;
      cap[idx] = basecap;
      /** 97 is "a" in ASCII code */
      sprintf(cap[idx].text, "(%c) %s", 97 + INDEX2D(nypanel, nxpanel, nypanel - jj - 1, ii), maccap[ii].text);
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */


  /** individual file names */
  sprintf(figfile, "%s_%s", file, "gain");
  for(int ii = 0; ii < nxpanel; ii++)
    sprintf(figname[ii], "%s_%s_%s", file, "gain", synonym[ii]);


  /** create figure(s) */
/* #ifndef PLOT_SUMMARY_FIGURE_ONLY */
  for(int ii = 0; ii < nxpanel; ii++){
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 0);
    const PLBOOL tmp = cap[idx].write;
    cap[idx].write = false;

    plotData(1, 1, &plot[idx], false, false,
	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
	     &xlabel[idx], &ylabel[idx], "", figname[ii], argc, argv);

    cap[idx].write = tmp;
  }/* for(int ii = 0; ii < nxpanel; ii++){ */
/* #endif//PLOT_SUMMARY_FIGURE_ONLY */

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
