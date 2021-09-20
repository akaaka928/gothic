/**
 * @file spheroid.c
 *
 * @brief Estimator for mass and effective radius within a dark matter halo stellar mass--halo mass relation by Behroozi et al. (2013, ApJ, 770, 57) or Shankar et al. (2017, ApJ, 840, 34) and size--stellar mass relation by Bernardi et al. (2014, MNRAS, 443, 874)
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2018/12/14 (Fri)
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


static inline double get_stellar_mass_Shankar2017(const double M200)
{
  const double alp = 2.15;
  const double bet = 1.77;
  const double Ms0 = 8.511380382e+10;/**< log(Ms0) = 10.93 */
  const double Mh0 = 6.918309709e+11;/**< log(Mh0) = 11.84 */

  const double xx = M200 / Mh0;

  /** obtain stellar mass (5.22) */
  const double Ms = Ms0 * pow(xx, alp) / (1.0 + pow(xx, bet));

  return (Ms);
}


static inline double support_function_Behroozi2013(const double xx, const double redshift, const double am1, const double nu)
{
  const double alp = -1.412 + 0.731 * am1 * nu;
  const double del = 3.508 + nu * (2.608 * am1 - 0.043 * redshift);
  const double gam = 0.316 + nu * (1.319 * am1 + 0.279 * redshift);

  const double val = -log10(1.0 + pow(10.0, alp * xx)) + del * pow(log10(1.0 + exp(xx)), gam) / (1.0 + exp(pow(10.0, -xx)));

  return (val);
}
static inline double get_stellar_mass_Behroozi2013(const double M200, const double redshift)
{
  /** set scale factor */
  const double scale_factor = 1.0 / (1.0 + redshift);
  const double am1 = scale_factor - 1.0;

  const double nu = exp(-4.0 * scale_factor * scale_factor);

  const double logeps = -1.777 - am1 * (0.119 + 0.006 * nu);
  const double logM1 = 11.514 - nu * (1.793 * am1 + 0.251 * redshift);

  const double logMs = logeps + logM1 + support_function_Behroozi2013(log10(M200) - logM1, redshift, am1, nu) - support_function_Behroozi2013(0.0, redshift, am1, nu);
  const double Mstar = pow(10.0, logMs);

  return (Mstar);
}


static inline double get_effective_radius(const double Mstar)
{
  const double p0 = 13.4131;/**< parameter for SerExp (early types) */
  const double p1 = -2.9324;/**< parameter for SerExp (early types) */
  const double p2 = 0.1607;/**< parameter for SerExp (early types) */

  /* const double p0 = 11.2699;/\**< parameter for SerExp (early-type bulges) *\/ */
  /* const double p1 = -2.3026;/\**< parameter for SerExp (early-type bulges) *\/ */
  /* const double p2 = 0.1227;/\**< parameter for SerExp (early-type bulges) *\/ */

  /** obtain effective radius (5.24) */
  const double logM = log10(Mstar);
  const double logR = p0 + logM * (p1 + p2 * logM);

  const double Reff = pow(10.0, logR);

  return (Reff);
}


int main(int argc, char **argv)
{
  /** initialization */
  if( argc < 2 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 2);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -M200=<double> (in units of solar mass)\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 2 ){ */


  /** read input argument */
  static double M200;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "M200", &M200));


  /** calculate stellar mass */
  const double Mstar_S17 = get_stellar_mass_Shankar2017(M200);
  const double Mstar_B13 = get_stellar_mass_Behroozi2013(M200, 0.0);

  /** calculate effective radius */
  const double Reff_S17 = get_effective_radius(Mstar_S17);
  const double Reff_B13 = get_effective_radius(Mstar_B13);


  /** output calculated results */
  fprintf(stdout, "Virial mass M200 = %e Msun\n", M200);
  fprintf(stdout, "Shankar+17: Stellar mass Mstar = %e Msun, effective radius Reff = %e kpc\n", Mstar_S17, Reff_S17);
  fprintf(stdout, "Behroozi+13: Stellar mass Mstar = %e Msun, effective radius Reff = %e kpc\n", Mstar_B13, Reff_B13);


  return (0);
}
