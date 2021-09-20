/**
 * @file dmhalo.c
 *
 * @brief Estimator for concentration parameter using c-M relation by Prada et al. (2012, MNRAS, 423, 3018)
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


double newton_cgs;
double Msun;
double kpc;

double hubble_nondim;/**< H_0 / 100 */
double Omega_Lambda;
double Omega_matter;

double hubble;
double rho_crit;


static inline double cmin(const double xx)
{
  const double c0 = 3.681;
  const double c1 = 5.033;
  const double AA = 6.948;
  const double x0 = 0.424;

  return (c0 + (c1 - c0) * (0.5 + M_1_PI * atan(AA * (xx - x0))));
}
static inline double smin(const double xx)
{
  const double s0 = 1.047;
  const double s1 = 1.646;
  const double BB = 7.386;
  const double x1 = 0.526;

  return (s0 + (s1 - s0) * (0.5 + M_1_PI * atan(BB * (xx - x1))));
}


static inline double get_DEformula(const double tt, const double xx)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double cosh_t = cosh(sinh_t);
  const double  exp_t =  exp(sinh_t) * xx;

  const double tmp = sqrt(exp_t / (8.0 * cosh_t * cosh_t * cosh_t + exp_t * exp_t * exp_t));

  return (cosh(tt) * cosh_t * tmp * tmp * tmp);
}
static inline double update_trapezoidal(const double hh, const double tmin, const double tmax, const double sum, const double xx)
{
  /** initialization */
  double sub = 0.0;
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    sub += get_DEformula(tt, xx);
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  return (0.5 * sum + hh * sub);
}
static inline double set_domain_boundary(const double hh, double * restrict tmin, double * restrict tmax, const double xx)
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  double fp = get_DEformula(tt, xx);
  const double f0 = fp;
  double sum = hh * fp;


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    const double ft = fp;

    tt += hh;    boundary = tt;
    fp = get_DEformula(tt, xx);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);
  }/* while( (damp > converge) && (boundary < maximum) ){ */
  *tmax = boundary;


  /** determine lower boundary */
  fp = f0;
  tt = 0.0;
  boundary = 0.0;
  damp = 1.0;
  while( (damp > converge) && (boundary > -maximum) ){
    const double ft = fp;

    tt -= hh;    boundary = tt;
    fp = get_DEformula(tt, xx);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;

  return (sum);
}
static inline double integrate_DEformula(const double xx)
{
  const double criteria_abs = 1.0e-16;
  /* const double criteria_abs = 1.0e-14; */
  /* const double criteria_abs = 1.0e-12; */
  const double criteria_rel = 1.0e-12;
  /* const double criteria_rel = 1.0e-10; */
  /* const double criteria_rel = 1.0e-8; */
  /* const double criteria_rel = 1.0e-6; */

  double hh = 1.0;
  double tmin, tmax;
  double sum = set_domain_boundary(hh, &tmin, &tmax, xx);

  while( true ){
    const double f0 = sum;

    hh *= 0.5;
    sum = update_trapezoidal(hh, tmin, tmax, sum, xx);

    bool converge = true;
    if( fabs(sum) > DBL_EPSILON ){
      if( fabs(1.0 - f0 / sum) > criteria_rel )
	converge = false;
    }
    else
      if( fabs(sum - f0) > criteria_abs )
	converge = false;

    if( converge )
      break;
  }/* while( true ){ */

  return (2.0 * M_PI * xx * sum);
}


static inline double get_concentration(const double M200, const double redshift)
{
  /** set scale factor */
  const double scale_factor = 1.0 / (1.0 + redshift);

  /** calculate x (5.1) */
  const double x0 = cbrt(Omega_Lambda / Omega_matter);
  const double xx = scale_factor * x0;

  /** calculate D(a) (5.2) */
  const double Da = 2.5 * cbrt(Omega_matter / Omega_Lambda) * sqrt(1.0 + 1.0 / (xx * xx * xx)) * integrate_DEformula(xx);
  const double D1 = 2.5 * cbrt(Omega_matter / Omega_Lambda) * sqrt(1.0 + 1.0 / (x0 * x0 * x0)) * integrate_DEformula(x0);

  /** set y (5.3) */
  const double yy = 1.0e+12 / (hubble_nondim * M200);

  /** calculate root-mean-square density fluctuation sigma(M, a) (5.3) */
  const double sigma = (Da / D1) * 16.9 * pow(yy, 0.41) / (1.0 + 1.102 * pow(yy, 0.20) + 6.22 * pow(yy, 0.333));
  /* const double sigma = Da * 16.9 * pow(yy, 0.41) / (1.0 + 1.102 * pow(yy, 0.20) + 6.22 * pow(yy, 0.333)); */
  /* const double sigma = 16.9 * pow(yy, 0.41) / (1.0 + 1.102 * pow(yy, 0.20) + 6.22 * pow(yy, 0.333)); */

  /** determine parameters B0(x) and B1(x) (5.4) and (5.5) */
  const double B0 = cmin(xx) / cmin(x0);
  const double B1 = smin(xx) / smin(x0);

  /** obtain sigma' (5.8) */
  const double sp = B1 * sigma;

  /** obtain \mathcal{C} (5.8) */
  const double alp = 2.881;
  const double bet = 1.257;
  const double gam = 1.022;
  const double del = 0.060;
  const double c0 = alp * ((1.0 + pow(sp / bet, gam)) * exp(del / (sp * sp)));/**< concentration at z = 0 */

  /** derive concentration parameter (5.10) */
  const double concentration = B0 * c0;

  return (concentration);
}


static inline double get_Virial_radius(const double M200)
{
  const double r200 = cbrt(3.0 * M200 * Msun / (800.0 * M_PI * rho_crit)) / kpc;/**< in units of kpc */

  return (r200);
}


static inline double get_scale_radius(const double r200, const double concentration)
{
  const double rs = r200 / concentration;

  return (rs);
}


int main(int argc, char **argv)
{
  /** initialization */
  if( argc < 3 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 3);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -M200=<double> (in units of solar mass)\n");
    __FPRINTF__(stderr, "          -redshift=<double>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 3 ){ */


  /** read input arguments */
  static double M200;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "M200", &M200));
  static double redshift;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "redshift", &redshift));


  /** set physical and astrophysical constants */
  newton_cgs = 6.67408e-8;/**< 2014 CODATA */
  Msun = 1.98847e+33;/**< gram */
  kpc = 3.08567758149137e+21;/**< cm */

  /** set cosmological parameters (adopt Planck Collaboration 2018: TT,TE,EE+lowE+lensing+BAO 68% limits) */
  hubble_nondim = 0.6766;/**< H_0 = (67.66 \pm 0.42) km s^-1 Mpc^-1 */
  Omega_Lambda  = 0.6889;
  Omega_matter  = 0.3111;
  /* hubble_nondim = 0.70; */
  /* Omega_Lambda = 0.73; */
  /* Omega_matter = 0.27; */

  /** set critical density */
  hubble = 100.0 * hubble_nondim * (1.0e+5 / (1.0e+3 * kpc));/**< Hubble constant in units of cgs */
  hubble *= sqrt(Omega_Lambda + Omega_matter * (1.0 + redshift) * (1.0 + redshift) * (1.0 + redshift));/**< Hubble parameter at redshift z */
  rho_crit = 3.0 * hubble * hubble / (8.0 * M_PI * newton_cgs);/**< in units of cgs */


  /** calculate concentration parameter */
  const double concentration = get_concentration(M200, redshift);

  /** calculate r200 (5.11) */
  const double r200 = get_Virial_radius(M200);

  /** obtain scale length (5.12) */
  const double rs = get_scale_radius(r200, concentration);


  /** output calculated results */
  fprintf(stdout, "Redshift z = %e\n", redshift);
  fprintf(stdout, "Virial mass M200 = %e Msun\n", M200);
  fprintf(stdout, "Concentration c = %e\n", concentration);
  fprintf(stdout, "Virial radius r200 = %e kpc\n", r200);
  fprintf(stdout, "Scale radius rs = %e kpc\n", rs);

  fprintf(stdout, "log M200 (h^-1 Msun) = %e\n", log10(M200 * hubble_nondim));
  fprintf(stdout, "log c = %e\n", log10(concentration));


  return (0);
}
