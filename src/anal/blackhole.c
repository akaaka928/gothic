/**
 * @file blackhole.c
 *
 * @brief Estimator for black hole mass using black hole mass--bulge mass relation by Kormendy & Ho (2013, ARA&A, 51, 511)
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2018/12/13 (Thu)
 *
 * Copyright (C) 2018 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "myutil.h"


static inline double get_blackhole_mass(const double Mbulge)
{
  /** obtain black hole mass (5.29) */
  const double Mbh = 0.49e+9 * pow(Mbulge * 1.0e-11, 1.17);

  return (Mbh);
}


int main(int argc, char **argv)
{
  /** initialization */
  if( argc < 2 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 2);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -Mbulge=<double> (in units of solar mass)\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 2 ){ */


  /** read input argument */
  static double Mbulge;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "Mbulge", &Mbulge));

  /** calculate black hole mass */
  const double Mbh = get_blackhole_mass(Mbulge);


  /** output calculated results */
  fprintf(stdout, "Bulge mass Mbulge = %e Msun\n", Mbulge);
  fprintf(stdout, "Black hole mass Mbh = %e Msun\n", Mbh);


  return (0);
}
