/**
 * @file merge.c
 *
 * @brief Source code for MAGI (MAny-component Galaxy Initializer)
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/05/31 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


int main(int argc, char **argv)
{
  /** read input arguments */
  if( argc < 9 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 9);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -list=<char *>\n");
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -eps=<real> -eta=<real>\n");
    __FPRINTF__(stderr, "          -ft=<real> -snapshotInterval=<real> -saveInterval=<real>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 9 ){ */

  /** read input arguments do not depend on the unit system adopted in the numerical simulation */
  char *file;  requiredCmdArg(getCmdArgStr( argc, (const char * const *)argv, "file", &file));
  char *fcfg;  requiredCmdArg(getCmdArgStr( argc, (const char * const *)argv, "list", &fcfg));
  ulong Ntot;  requiredCmdArg(getCmdArgUlg( argc, (const char * const *)argv, "Ntot", &Ntot));
  real   eta;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,  "eta", &eta));
  double saveInterval;  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "saveInterval", &saveInterval));


  /* open the list, get # of files, make list of individual files, get total # of particles */
  /* allocate memory */
  /* read individual files, check unit system */
  /* how about shift vector */




  /** below are taken from gothic.c */

  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv,  "file", &file));
  double tmp;



  /** read settings about the simulation */
  static ulong Ntot;
  static real eps, eta;
  static double ft, snapshotInterval, saveInterval;
  static int unit;
#ifdef  SERIALIZED_EXECUTION
  readSettings        (&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file);
#else///SERIALIZED_EXECUTION
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file, mpi);
#endif//SERIALIZED_EXECUTION

  /** read setting to dump tentative results of the simulation */
  static int last;
#ifdef  SERIALIZED_EXECUTION
  readConfigFile        (&last, file);
#else///SERIALIZED_EXECUTION
  readConfigFileParallel(&last, file, mpi);
#endif//SERIALIZED_EXECUTION

}
