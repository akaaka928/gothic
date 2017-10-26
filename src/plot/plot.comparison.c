/**
 * @file plot.comparison.c
 *
 * @brief Plot code of cumulative distribution function for tree code
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
#include <mpi.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "plplotlib.h"
#include "name.h"


#define TYPICAL_ERROR_FOR_GALAXY_SCALE

#define INVERSE_VER
#define LOGPLOT_HOR
/* #define LOGPLOT_VER */

#   if  defined(INVERSE_VER) && !defined(LOGPLOT_VER)
#undef  INVERSE_VER
#endif//defined(INVERSE_VER) && !defined(LOGPLOT_VER)


#define NTOT (8388608)
#define PROJECTNAME "m31"


void compareCDF
(const int num, PLFLT bhcdf[restrict], PLFLT bhacc[restrict], PLFLT wscdf[restrict], PLFLT wsacc[restrict], PLFLT g2cdf[restrict], PLFLT g2acc[restrict],
 PLplotPltRange box, char file[], int argc, char **argv);


int main(int argc, char **argv)
{
  if( argc < 1 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 1);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 1 ){ */

  modifyArgcArgv4PLplot(&argc, argv, 1);


  /** read analyzed dataset */
  static PLFLT bhcdf[NTOT], bhacc[NTOT];
  static PLFLT wscdf[NTOT], wsacc[NTOT];
  static PLFLT g2cdf[NTOT], g2acc[NTOT];

#ifdef  USE_HDF5_FORMAT
  /** open HDF5 files */
#ifdef  TYPICAL_ERROR_FOR_GALAXY_SCALE
  /** typical accuracy for galactic scale N-body simulations, initial condition is M31 */
  char bhfile[128];  sprintf(bhfile, "%s/%s/%s.%s.%s.%.2d.%s.%.3d.h5", DATAFOLDER, "cc35/compare", PROJECTNAME, "tree", "bh", 6, "err.cdf", 0);
  char wsfile[128];  sprintf(wsfile, "%s/%s/%s.%s.%s.%.2d.%s.%.3d.h5", DATAFOLDER, "cc35/compare", PROJECTNAME, "tree", "ws", 2, "err.cdf", 0);
  char g2file[128];  sprintf(g2file, "%s/%s/%s.%s.%s.%.2d.%s.%.3d.h5", DATAFOLDER, "cc35/compare", PROJECTNAME, "tree", "g2", 7, "err.cdf", 0);
#else///TYPICAL_ERROR_FOR_GALAXY_SCALE
  char bhfile[128];  sprintf(bhfile, "%s/%s/%s.%s.%s.%.2d.%s.%.3d.h5", DATAFOLDER, "cc35/compare", PROJECTNAME, "tree", "bh",  1, "err.cdf", 0);
  char wsfile[128];  sprintf(wsfile, "%s/%s/%s.%s.%s.%.2d.%s.%.3d.h5", DATAFOLDER, "cc35/compare", PROJECTNAME, "tree", "ws", 14, "err.cdf", 0);
  char g2file[128];  sprintf(g2file, "%s/%s/%s.%s.%s.%.2d.%s.%.3d.h5", DATAFOLDER, "cc35/compare", PROJECTNAME, "tree", "g2", 17, "err.cdf", 0);
#endif//TYPICAL_ERROR_FOR_GALAXY_SCALE
  /** typical accuracy for cosmological N-body simulations, initial condition is M31 */
  hid_t bhtarget = H5Fopen(bhfile, H5F_ACC_RDONLY, H5P_DEFAULT);  hid_t bhgroup = H5Gopen(bhtarget, "relErr", H5P_DEFAULT);
  hid_t wstarget = H5Fopen(wsfile, H5F_ACC_RDONLY, H5P_DEFAULT);  hid_t wsgroup = H5Gopen(wstarget, "relErr", H5P_DEFAULT);
  hid_t g2target = H5Fopen(g2file, H5F_ACC_RDONLY, H5P_DEFAULT);  hid_t g2group = H5Gopen(g2target, "relErr", H5P_DEFAULT);
  /** read CDF */
  hid_t dataset;
  dataset = H5Dopen(bhgroup, "cdf", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, bhcdf));  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(bhgroup, "acc", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, bhacc));  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(wsgroup, "cdf", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, wscdf));  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(wsgroup, "acc", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, wsacc));  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(g2group, "cdf", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, g2cdf));  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(g2group, "acc", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, g2acc));  chkHDF5err(H5Dclose(dataset));
  /** close HDF5 files */
  chkHDF5err(H5Gclose(bhgroup));  chkHDF5err(H5Fclose(bhtarget));
  chkHDF5err(H5Gclose(wsgroup));  chkHDF5err(H5Fclose(wstarget));
  chkHDF5err(H5Gclose(g2group));  chkHDF5err(H5Fclose(g2target));
#endif//USE_HDF5_FORMAT

#   if  defined(LOGPLOT_HOR) || defined(LOGPLOT_VER)
  for(int ii = 0; ii < NTOT; ii++){
#ifdef  LOGPLOT_HOR
    bhacc[ii] = log10(bhacc[ii]);
    wsacc[ii] = log10(wsacc[ii]);
    g2acc[ii] = log10(g2acc[ii]);
#endif//LOGPLOT_HOR
#ifdef  LOGPLOT_VER
#ifdef  INVERSE_VER
    bhcdf[ii] = log10(1.0 - bhcdf[ii]);
    wscdf[ii] = log10(1.0 - wscdf[ii]);
    g2cdf[ii] = log10(1.0 - g2cdf[ii]);
#else///INVERSE_VER
    bhcdf[ii] = log10(bhcdf[ii]);
    wscdf[ii] = log10(wscdf[ii]);
    g2cdf[ii] = log10(g2cdf[ii]);
#endif//INVERSE_VER
#endif//LOGPLOT_VER
  }/* for(int ii = 0; ii < NTOT; ii++){ */
#endif//defined(LOGPLOT_HOR) || defined(LOGPLOT_VER)


  /** plot the cumulative distribution function */
  PLplotPltRange range;
#ifdef  LOGPLOT_HOR
#ifdef  TYPICAL_ERROR_FOR_GALAXY_SCALE
  range.xmin = log10(1.0e-4);
  range.xmax = log10(1.0e-2);
#else///TYPICAL_ERROR_FOR_GALAXY_SCALE
  range.xmin = log10(7.0e-7);
  range.xmax = log10(3.0e-4);
#endif//TYPICAL_ERROR_FOR_GALAXY_SCALE
  range.xlog = LOGARITHMIC_PLOT;
#else//LOGPLOT_HOR
  range.xmin = 0.0;
  range.xmax = 1.0;
  range.xlog = LINEAR_PLOT;
#endif//LOGPLOT_HOR
#ifdef  LOGPLOT_VER
  range.ymin = log10(1.0e-6);
  range.ymax = log10(1.0e-0);
  range.ylog = LOGARITHMIC_PLOT;
#else///LOGPLOT_VER
  range.ymin = 0.0;
  range.ymax = 1.0;
  range.ylog = LINEAR_PLOT;
#endif//LOGPLOT_VER
  range.xgrd = range.ygrd = true;


  compareCDF(NTOT, bhcdf, bhacc, wscdf, wsacc, g2cdf, g2acc, range, PROJECTNAME, argc, argv);


  return (0);
}


void compareCDF
(const int num, PLFLT bhcdf[restrict], PLFLT bhacc[restrict], PLFLT wscdf[restrict], PLFLT wsacc[restrict], PLFLT g2cdf[restrict], PLFLT g2acc[restrict],
 PLplotPltRange box, char file[], int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /** global setting(s) */
  /** maximum number of data kinds */
  const PLINT pkind = 0;
  const PLINT lkind = 3;

  const PLBOOL Auni = true;/**< A stands for acceleration */


  /** preparation for dataset */
  /** memory allocation for index */
  PLINT Akind = (Auni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  PLINT *Anum;  allocPLINT(&Anum, Akind);

  /** set number of data points */
  for(PLINT ii = 0; ii < Akind; ii++)    Anum[ii] = num;

  /** memory allocation for data */
  PLFLT **hor;  allocPointer4PLFLT(&hor, Akind);
  PLFLT **ver;  allocPointer4PLFLT(&ver, Akind);
  /** data preparation */
  if( Auni ){
    hor[0] = g2acc;    ver[0] = g2cdf;
    hor[1] = wsacc;    ver[1] = wscdf;
    hor[2] = bhacc;    ver[2] = bhcdf;
  }/* if( Auni ){ */
  else{
    __KILL__(stderr, "ERROR: not supported value for Auni(%d)\n", Auni);
  }/* else{ */

  /** set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  ls[0].color =   RED;  ls[0].style =  SOLID_LINE;  ls[0].width = MIDDLE_LINE;
  ls[1].color =  BLUE;  ls[1].style = DASHED_LINE;  ls[1].width = MIDDLE_LINE;
  ls[2].color = BLACK;  ls[2].style = DOTTED_LINE;  ls[2].width = MIDDLE_LINE;

  /** set labels */
  char Alab[PLplotCharWords];  sprintf(Alab, "|#fi#<0x12>a#fr#di#u#utree#d - #fi#<0x12>a#fr#di#u#udirect#d| / #fia#fr#di#u#udirect#d");
#ifdef  INVERSE_VER
  char ylab[PLplotCharWords];  sprintf(ylab, "1 - #fiF#fr");
#else///INVERSE_VER
  char ylab[PLplotCharWords];  sprintf(ylab, "#fiF#fr");
#endif//INVERSE_VER

  /** set caption(s) */
  PLplotCaption Acap;  setDefaultCaption(&Acap);  Acap.write = false;
  sprintf(Acap.side, "%s", "t");
  sprintf(Acap.text, "%s", "CDF for acceleration error");

  /** set legends */
  PLplotLegend Aleg;  setDefaultLegend(&Aleg, false);  Aleg.write = true;
  char *AlegTxt;
  {
    allocChar4PLplot(&AlegTxt, Akind);
    allocPointer4Char4PLplot(&(Aleg.text), Akind);
    assignChar4PLplot(Akind, Aleg.text, AlegTxt);
  }
  sprintf(Aleg.text[0], "%s", "acc");
  sprintf(Aleg.text[1], "%s", "mul");
  sprintf(Aleg.text[2], "%s", "#gh");
#ifdef  INVERSE_VER
  Aleg.pos = PL_POSITION_LEFT | PL_POSITION_BOTTOM | PL_POSITION_INSIDE;
#else///INVERSE_VER
  Aleg.pos = PL_POSITION_LEFT | PL_POSITION_TOP    | PL_POSITION_INSIDE;
#endif//INVERSE_VER


  /** memory allocation to exploit multiple-plot function */
  /** maximum number of panels in each direction */
  const PLINT nxpanel = 1;
  const PLINT nypanel = 1;

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


  /** configure to cumulative distribution function of relative error */
  /** common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

      /** global setting(s) */
      plot[idx] = true;

      /** line setting(s) */
      nlkind[idx] = IMIN(Akind, lkind);
      line  [idx] = ls;
      lnum  [idx] = Anum;
      lx    [idx] = hor;
      ly    [idx] = ver;

      /** point setting(s) */
      npkind[idx] = IMIN(Akind, pkind);
      point [idx] = pt;
      pnum  [idx] = Anum;
      px    [idx] = hor;
      py    [idx] = ver;

      /** errorbar setting(s) */
      errbar[idx] = false;

      /** plot area */
      range[idx] = box;

      /** label setting(s) */
      xlabel[idx] = Alab;
      ylabel[idx] = ylab;

      /** miscs */
      cap[idx] = Acap;
      leg[idx] = Aleg;
      uni[idx] = Auni;
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */

  /** individual file names */
  sprintf(figfile, "%s_%s_%s", file, "cdf", "comparison");


  /** create figure(s) */
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

  free(Anum);  free(hor);  free(ver);
  free(ls);  free(pt);

  free(AlegTxt);


  __NOTE__("%s\n", "end");
}
