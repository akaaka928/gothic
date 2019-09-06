/**
 * @file extend.c
 *
 * @brief Source code for extending final time of N-body simulations
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2019/09/03 (Tue)
 *
 * Copyright (C) 2018 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "myutil.h"
#include "constants.h"

#include "../file/io.h"


int main(int argc, char **argv)
{
  /** configure the details of the numerical simulation */
  /** read command line arguments */
  if( argc < 3 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 3);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -ft=<double>\n");
    __FPRINTF__(stderr, "          -saveInterval=<double> (optional)\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }
  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv,  "file", &file));


  /** read settings about the simulation */
  static ulong Ntot;
  static real eps, eta;
  static double ft, snapshotInterval, saveInterval;
  static int unit;
  readSettings(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file);


  /** read input arguments depend on the unit system adopted in the numerical simulation */
  extern const double time_astro2com;
  double tmp;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "ft",  &tmp));
  ft = tmp * time_astro2com;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "saveInterval", &tmp)) == myUtilAvail )
    saveInterval = tmp;


  writeSettings(unit, Ntot, eps, eta, ft, snapshotInterval, saveInterval, file);


  return (0);
}
