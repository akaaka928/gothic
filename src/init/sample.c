/**
 * @file sample.c
 *
 * @brief Source code for generating sample of machine-readable table for MAGI
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/02/24 (Fri)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "macro.h"


/** CFGFOLDER must same with CFGFOLDER defined in magi.h */
#define CFGFOLDER "cfg"
#define DSTDIR "table"

/* /\** ~10% error *\/ */
/* #define NRAD (64) */
/* #define NHOR (64) */

/** ~1% error */
#define NRAD (128)
#define NHOR (128)

/* /\** ~0.1% error *\/ */
/* #define NRAD (256) */
/* #define NHOR (256) */


/** based on Ciotti & Bertin (1999), A&A, 352, 447-451: Eq.(25) */
static inline double getAsymptoticSersicScale(const double nn){  return (2.0 * nn + (-1.0 + 2.0 * (2.0 + 23.0 / (63.0 * nn)) / (135.0 * nn)) / 3.0);}


int main(void)
{
  /** set table range */
  /* const real rmin = 1.0e-6; */
  const real rmin = 1.0e-5;
  /* const real rmin = 1.0e-4; */
  /* const real rmin = 1.0e-3; */
  /* const real rmax = 1.0e+3; */
  /* const real rmin = 1.0e-2; */
  /* const real rmax = 1.0e+2; */
  /* const real rmin = 1.0e-3; */
  const real rmax = 5.0e+1;

  const real logrmin = LOG10(rmin);
  const real logrmax = LOG10(rmax);
  const real logrbin = (logrmax - logrmin) / (real)(NRAD - 1);
  const real logRbin = (logrmax - logrmin) / (real)(NHOR - 1);


  /** set density profile in table form */
  static real xx[NRAD], core[NRAD], cusp[NRAD], twop[NRAD];
  static real RR[NHOR], proj[NHOR], disk[NHOR];

  for(int ii = 0; ii < NRAD; ii++)
    xx[ii] = POW10(logrmin + logrbin * (real)ii);

  /** Plummer model */
  for(int ii = 0; ii < NRAD; ii++)
    core[ii] = POW(UNITY + xx[ii] * xx[ii], -HALF * FIVE);

  /** Hernquist model */
  for(int ii = 0; ii < NRAD; ii++)
    cusp[ii] = UNITY / (xx[ii] * (UNITY + xx[ii]) * (UNITY + xx[ii]) * (UNITY + xx[ii]));

  /** Two-power model */
  const real alpha = TWO;
  const real  beta = TWO;  const real binv = UNITY / TWO;
  const real   gam = FOUR;
  for(int ii = 0; ii < NRAD; ii++)
    twop[ii] = POW(xx[ii], -alpha) * POW(UNITY + POW(xx[ii], beta), (alpha - gam) * binv);


  /** set column density profile in table form */
  for(int ii = 0; ii < NHOR; ii++)
    RR[ii] = POW10(logrmin + logRbin * (real)ii);

  /** de Vaucouleurs profile */
  const real bb = (real)getAsymptoticSersicScale(4.0);
  for(int ii = 0; ii < NHOR; ii++)
    proj[ii] = EXP(-bb * POW(RR[ii], QUARTER));

  /** exponential disk */
  for(int ii = 0; ii < NHOR; ii++)
    disk[ii] = EXP(-RR[ii]);


  /** output ascii file */
  FILE *fp;
  char filename[128];

  sprintf(filename, "%s/%s/%s.dat", CFGFOLDER, DSTDIR, "tplummer");
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fprintf(fp, "%d\n", NRAD);
  for(int ii = 0; ii < NRAD; ii++)
    fprintf(fp, "%e\t%e\n", xx[ii], core[ii]);
  fclose(fp);

  sprintf(filename, "%s/%s/%s.dat", CFGFOLDER, DSTDIR, "thernquist");
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fprintf(fp, "%d\n", NRAD);
  for(int ii = 0; ii < NRAD; ii++)
    fprintf(fp, "%e\t%e\n", xx[ii], cusp[ii]);
  fclose(fp);

  sprintf(filename, "%s/%s/%s.dat", CFGFOLDER, DSTDIR, "dblpower");
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fprintf(fp, "%d\n", NRAD);
  for(int ii = 0; ii < NRAD; ii++)
    fprintf(fp, "%e\t%e\n", xx[ii], twop[ii]);
  fclose(fp);

  sprintf(filename, "%s/%s/%s.dat", CFGFOLDER, DSTDIR, "bulgetbl");
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fprintf(fp, "%d\n", NHOR);
  for(int ii = 0; ii < NHOR; ii++)
    fprintf(fp, "%e\t%e\n", RR[ii], proj[ii]);
  fclose(fp);

  sprintf(filename, "%s/%s/%s.dat", CFGFOLDER, DSTDIR, "disktbl");
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fprintf(fp, "%d\n", NHOR);
  for(int ii = 0; ii < NHOR; ii++)
    fprintf(fp, "%e\t%e\n", RR[ii], disk[ii]);
  fclose(fp);


  return (0);
}
