/**
 * @file plot.distribution.c
 *
 * @brief Plot code for radial profiles of multiple components
 *
 * @author Yohei Miki (University of Tsukuba)
 *
 * @date 2017/08/14 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#define HEAD_GROUP_IDX (1)

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "myutil.h"
#include "name.h"
#include "constants.h"
#include "mpilib.h"
#include "plplotlib.h"

#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../file/io.h"

#include "../init/spline.h"


extern const double      length2astro;extern const char      length_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double        time2astro;extern const char        time_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double        mass2astro;extern const char        mass_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double    velocity2astro;extern const char    velocity_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];

static const int Nminimum = 16;



typedef struct
{
  ulong idx;
  real x, y, z;
  real vx, vy, vz;
  real m, pot;
  real rad, hor;
  int gid;
} nbody_particle;





void getCenter(const int num, nbody_particle *body, double com[restrict], double vel[restrict]);


#ifdef __ICC
/* Disable ICC's remark #161: unrecognized #pragma */
#     pragma warning (disable:161)
#endif//__ICC

int radAscendingOrder(const void *a, const void *b);
int radAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((nbody_particle *)a)->rad > ((nbody_particle *)b)->rad ){    return ( 1);  }
  else{    if( ((nbody_particle *)a)->rad < ((nbody_particle *)b)->rad ){    return (-1);  }
    else{                                                                    return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

int gidAscendingOrder(const void *a, const void *b);
int gidAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((nbody_particle *)a)->gid > ((nbody_particle *)b)->gid ){    return ( 1);  }
  else{    if( ((nbody_particle *)a)->gid < ((nbody_particle *)b)->gid ){    return (-1);  }
    else{                                                                    return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

#ifdef __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC


static inline double getEffectivePotential(const double rr, const double L2, const int Nspline, double * restrict xx, double * restrict yy, double * restrict y2)
{
  /* const double rinv = 1.0 / (DBL_MIN + rr); */
  const double rinv = 1.0 / (1.0e-100 + rr);
  return (getCubicSpline1D(rr, Nspline, xx, yy, y2) + 0.5 * L2 * rinv * rinv);
}


static inline double getRadialKineticEnergy(const double rr, const double Etot, const double L2, const int Nspline, double * restrict xx, double * restrict yy, double * restrict y2)
{
  return (2.0 * (Etot - getEffectivePotential(rr, L2, Nspline, xx, yy, y2)));
}


/**
 * @fn sign
 *
 * @brief Returns sign(b) * abs(a).
 *
 * @param (aa) variable a
 * @param (bb) variable b
 * @return sign(bb) * fabs(aa)
 */
static inline double sign(const double aa, const double bb){  return ((bb >= 0.0) ? (fabs(aa)) : (-fabs(aa)));}


/**
 * @fn brent_root
 *
 * @brief Find the root of the target function.
 *
 */
static inline double brent_root
(const double rmin, const double rmax, const double Etot, const double L2, const int Nspline, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const int Niter_max = 16384;
  const double tol = 1.0e-10;

  /** find root using Brent's method */
  double dd =  0.0;/**< displacement in the previous step */
  double ee =  0.0;/**< displacement in the step prior to the previous one*/
  double aa = rmin;  double fa = getRadialKineticEnergy(aa, Etot, L2, Nspline, xx, yy, y2);
  double bb = rmax;  double fb = getRadialKineticEnergy(bb, Etot, L2, Nspline, xx, yy, y2);
  double cc =   bb;  double fc = fb;

  if( fa * fb > 0.0 ){
    /* for(int ii = 0; ii < Nspline; ii++) */
    /*   fprintf(stdout, "%e\t%e\n", xx[ii], getEffectivePotential(xx[ii], L2, Nspline, xx, yy, y2)); */
    __KILL__(stderr, "ERROR: root may not be bracketed in [%e, %e]; fa = %e, fb = %e; Etot = %e, L2 = %e\n", aa, bb, fa, fb, Etot, L2);
  }/* if( fa * fb > 0.0 ){ */


  int iter = 0;
  while( true ){

    if( (fb * fc) > 0.0 ){
      cc = aa;      fc = fa;
      ee = dd = bb - aa;
    }/* if( (fb * fc) > 0.0 ){ */

    if( fabs(fc) < fabs(fb) ){
      aa = bb;      bb = cc;      cc = aa;
      fa = fb;      fb = fc;      fc = fa;
    }/* if( fabs(fc) < fabs(fb) ){ */

    const double tol1 = 2.0 * DBL_EPSILON * fabs(bb) + 0.5 * tol;
    const double mm = 0.5 * (cc - bb);

    if( (fabs(mm) <= tol1) || (fb == 0.0) )
      return (bb);

    if( (fabs(ee) >= tol1) && (fabs(fa) > fabs(fb)) ){
      /** try inverse quadratic interpolation */
      const double ss = fb / fa;
      double pp, qq;
      if( aa == cc ){
	pp = 2.0 * mm * ss;
	qq = 1.0 - ss;
      }/* if( za == zc ){ */
      else{
	const double inv = 1.0 / fc;
	const double rr = fb * inv;
	qq = fa * inv;
	pp = ss * (2.0 * mm * qq * (qq - rr) - (bb - aa) * (rr - 1.0));
	qq = (qq - 1.0) * (rr - 1.0) * (ss - 1.0);
      }/* else{ */

      if( pp > 0.0 )
	qq = -qq;
      pp = fabs(pp);

      /** validate the result of the inverse quadratic interpolation */
      if( (2.0 * pp) < fmin(3.0 * mm * qq - fabs(tol1 * qq), fabs(ee * qq)) ){
	/** accept the inverse quadratic interpolation */
	ee = dd;
	dd = pp / qq;
      }
      else{
	/** reject the inverse quadratic interpolation and adopt the bisection method */
	dd = mm;
	ee = dd;
      }/* else{ */
    }/* if( (fabs(ee) >= tol1) && (fabs(fa) > fabs(fb)) ){ */
    else{
      /** adopt the bisection method */
      dd = mm;
      ee = dd;
    }/* else{ */

    aa = bb;
    fa = fb;

    bb += ((fabs(dd) > tol1) ? (dd) : (sign(tol1, mm)));
    fb = getRadialKineticEnergy(bb, Etot, L2, Nspline, xx, yy, y2);

    iter++;
    if( iter > Niter_max ){
      __KILL__(stderr, "ERROR: Brent's method was not converged in %d steps.\n", iter);
    }/* if( iter > Niter_max ){ */
  }/* while( true ){ */
}


/**
 * @fn brent_minimum
 *
 * @brief Find the minimum of the target function.
 *
 */
static inline double brent_minimum
(double * restrict f_min, const double rmin, const double rmax, const double L2, const int Nspline, double * restrict rr, double * restrict yy, double * restrict y2)
{
  const int Niter_max = 128;
  const double tol = 1.0e-10;

  const double gold = 0.5 * (3.0 - sqrt(5.0));/**< the golden ratio */

  /** find minimum using Brent's method */
  double dd =  0.0;/**< displacement in the previous step */
  double ee =  0.0;/**< displacement in the step prior to the previous one*/

  double aa = fmin(rmin, rmax);
  double bb = fmax(rmin, rmax);

  double xx, ww, vv;
  double fx, fw, fv;
  xx = ww = vv = bb;
  fx = fw = fv = getEffectivePotential(xx, L2, Nspline, rr, yy, y2);


  int iter = 0;
  while( true ){
    const double xm = 0.5 * (aa + bb);
    const double tol0 = tol * fabs(xx) + DBL_EPSILON;
    const double tol1 = 2.0 * tol0;


    /* convergence check */
    if( fabs(xx - xm) <= (tol1 - 0.5 * (bb - aa)) ){
      /* fprintf(stdout, "# iter = %d, xx = %e, xm = %e, aa = %e, bb = %e, left = %e, right = %e, tol1 = %e, sub = %e\n", iter, xx, xm, aa, bb, fabs(xx - xm), tol1 - 0.5 * (bb - aa), tol1, 0.5 * (bb - aa)); */
      *f_min = fx;
      return xx;
    }/* if( fabs(xx - xm) <= (tol1 - 0.5 * (bb - aa)) ){ */
    if( iter > Niter_max ){
      __KILL__(stderr, "ERROR: Brent's method was not converged in %d steps (current error is %e; rmin = %e, rmax = %e, a = %e, b = %e).\n", iter, fabs(xx - xm), rmin, rmax, aa, bb);
    }/* if( iter > Niter_max ){ */


    /* execute interpolation */
    if( fabs(ee) > tol0 ){
      /* try parabolic interpolation */
      double rr = (xx - ww) * (fx - fv);
      double qq = (xx - vv) * (fx - fw);
      double pp = (xx - vv) * qq - (xx - ww) * rr;
      qq = 2.0 * (qq - rr);
      if( qq > 0.0 )	pp = -pp;
      else	        qq = -qq;
      const double tmp = ee;
      ee = dd;

      /* validate the result of the parabolic interpolation */
      if( (fabs(pp) >= fabs(0.5 * qq * tmp)) || (pp <= qq * (aa - xx)) || (pp >= qq * (bb - xx)) ){
	/* reject the parabolic interpolation and adopt the golden section search */
	ee = ((xx >= xm) ? (aa - xx) : (bb - xx));
	dd = gold * ee;
      }/* if( (fabs(pp) >= fabs(0.5 * qq * tmp)) || (pp <= qq * (aa - xx)) || (pp >= qq * (bb - xx)) ){ */
      else{
	/* accept the parabolic interpolation */
	dd = pp / qq;
	double tt = xx + dd;
	if( ((tt - aa) < tol1) || ((bb - tt) < tol1) )
	  dd = sign(tol0, xm - xx);
      }/* else{ */
    }/* if( fabs(ee) > tol0 ){ */
    else{
      /* adopt the golden section search */
      ee = ((xx >= xm) ? (aa - xx) : (bb - xx));
      dd = gold * ee;
    }/* else{ */


    const double uu = ((fabs(dd) >= tol0) ? (xx + dd) : (xx + sign(tol0, dd)));
    const double fu = getEffectivePotential(uu, L2, Nspline, rr, yy, y2);


    if( fu <= fx ){
      if( uu >= xx )	aa = xx;
      else	        bb = xx;

      vv = ww;      fv = fw;
      ww = xx;      fw = fx;
      xx = uu;      fx = fu;
    }/* if( fu <= fx ){ */
    else{
      if( uu < xx )	aa = uu;
      else	        bb = uu;

      if( (fu <= fw) || (ww == xx) ){
	vv = ww;	fv = fw;
	ww = uu;	fw = fu;
      }/* if( (fu <= fw) || (ww == xx) ){ */
      else
	if( (fu <= fv) || (vv == xx) || (vv == ww) ){
	  vv = uu;
	  fv = fu;
	}/* if( (fu <= fv) || (vv == xx) || (vv == ww) ){ */
    }/* else{ */

    iter++;
  }/* while( true ){ */
}






/**
 * @fn get_DEformula
 *
 * @brief Calculate integrand in double exponential formula.
 *
 * @param (tt) value of integration variable
 * @return (ret) integrand for Eddington formula
 * @param (max_pls_min) psi_max + psi_min
 * @param (max_mns_min) psi_max - psi_min
 * @param (xx) position of data points (psi)
 * @param (yy) value of data points (d2rho_dpsi2)
 * @param (y2) coefficients in cubic spline interpolation
 */
static inline double get_DEformula(const double tt, const double rmin, const double rmax, const double Etot, const double L2, const int Nspline, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double cosh_t = cosh(sinh_t);
  const double inv_cosh_t = 1.0 / cosh_t;

  const double one_pls_x = exp( sinh_t) * inv_cosh_t;
  const double one_mns_x = exp(-sinh_t) * inv_cosh_t;

  const double rr = 0.5 * (rmin * one_mns_x + rmax * one_pls_x);

#if 1
  /* workaround */
  const double ene = fmax(0.0, getRadialKineticEnergy(rr, Etot, L2, Nspline, xx, yy, y2));
#else
  const double ene = getRadialKineticEnergy(rr, Etot, L2, Nspline, xx, yy, y2);
#endif
#ifndef NDEBUG
  if( ene < 0.0 ){
    __KILL__(stderr, "val = %e, rmin = %e, rmax = %e, rr = %e, Etot = %e, L2 = %e, Phi = %e\n", ene, rmin, rmax, rr, Etot, L2, getCubicSpline1D(rr, Nspline, xx, yy, y2));
  }
#endif//NDEBUG

  return (cosh(tt) * inv_cosh_t * inv_cosh_t * sqrt(ene));
}


static inline double update_trapezoidal(const double hh, const double tmin, const double tmax, const double sum, const double rmin, const double rmax, const double Etot, const double L2, const int Nspline, double * restrict xx, double * restrict yy, double * restrict y2)
{
  /** initialization */
  double sub = 0.0;
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    sub += get_DEformula(tt, rmin, rmax, Etot, L2, Nspline, xx, yy, y2);
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

#ifndef NDEBUG
  if( fpclassify(0.5 * sum + hh * sub) != FP_NORMAL ){
    __KILL__(stderr, "sum = %e, hh = %e, sub = %e\n", sum, hh, sub);
  }
#endif//NDEBUG

  return (0.5 * sum + hh * sub);
}


static inline double set_domain_boundary(const double hh, double * restrict tmin, double * restrict tmax, const double rmin, const double rmax, const double Etot, const double L2, const int Nspline, double * restrict xx, double * restrict yy, double * restrict y2)
{
  /* const double converge = 1.0e-16; */
  const double converge = 1.0e-12;
  const double maximum = 128.0;

  double tt = 0.0;
  double fp = get_DEformula(tt, rmin, rmax, Etot, L2, Nspline, xx, yy, y2);
  const double f0 = fp;
  double sum = hh * fp;


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    const double ft = fp;

    tt += hh;    boundary = tt;
    fp = get_DEformula(tt, rmin, rmax, Etot, L2, Nspline, xx, yy, y2);

/* #ifndef NDEBUG */
/*     if( fpclassify(fp) != FP_NORMAL ){ */
/*       __KILL__(stderr, "fp = %e, tt = %e\n", fp, tt); */
/*     } */
/* #endif//NDEBUG */

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
    fp = get_DEformula(tt, rmin, rmax, Etot, L2, Nspline, xx, yy, y2);

/* #ifndef NDEBUG */
/*     if( fpclassify(fp) != FP_NORMAL ){ */
/*       __KILL__(stderr, "fp = %e, tt = %e\n", fp, tt); */
/*     } */
/* #endif//NDEBUG */

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;

#ifndef NDEBUG
  if( fpclassify(sum) != FP_NORMAL ){
    __KILL__(stderr, "sum = %e, tmin = %e, tmax = %e, hh = %e, f0 = %e\n", sum, *tmin, *tmax, hh, f0);
  }
#endif//NDEBUG

  return (sum);
}


static inline double integrate_DEformula(const double rmin, const double rmax, const double Etot, const double L2, const int Nspline, double * restrict xx, double * restrict yy, double * restrict y2)
{
  const double criteria_abs = 1.0e-12;
  /* const double criteria_rel = 1.0e-10; */
  /* const double criteria_rel = 1.0e-8; */
  const double criteria_rel = 1.0e-6;
  /* const double criteria_rel = 1.0e-5; */
  /* const double criteria_rel = 1.0e-4; */
  /* const double criteria_rel = 1.0e-3; */

  double hh = 1.0;
  double tmin, tmax;
  double sum = set_domain_boundary(hh, &tmin, &tmax, rmin, rmax, Etot, L2, Nspline, xx, yy, y2);


  while( true ){
    const double f0 = sum;

    hh *= 0.5;
    sum = update_trapezoidal(hh, tmin, tmax, sum, rmin, rmax, Etot, L2, Nspline, xx, yy, y2);

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

  sum *= 0.25 * (rmax - rmin);

#ifndef NDEBUG
  if( fpclassify(sum) != FP_NORMAL ){
    __KILL__(stderr, "sum = %e, rmax = %e, rmin = %e\n", sum, rmax, rmin);
  }
#endif//NDEBUG

  return (sum);
}


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);

  /** initialization */
  if( argc < 6 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 6);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -ncrit=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 6 ){ */

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));
  int    ncrit;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "ncrit", &ncrit));

  modifyArgcArgv4PLplot(&argc, argv, 6);


  /** load global settings of particle distribution */
  int last;
  readConfigFile(&last, file);
  int unit;
  ulong Ntot;
  real eta, eps;
  double ft, snapshotInterval, saveInterval;
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file, mpi);
  nbody_particle *body;
  body = (nbody_particle *)malloc(sizeof(nbody_particle) * Ntot);
  if( body == NULL ){    __KILL__(stderr, "ERROR: failure to allocate body\n");  }
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
  nbody_hdf5 hdf5;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
  ulong *hdf5_idx;
  allocSnapshotArray(&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx, (int)Ntot, &hdf5);
#else///USE_HDF5_FORMAT
  iparticle ibody;
  ulong *idx;
  position *pos;
  acceleration *acc;
#ifdef  BLOCK_TIME_STEP
  velocity *vel;
  ibody_time *ti;
#else///BLOCK_TIME_STEP
  real *vx, *vy, *vz;
#endif//BLOCK_TIME_STEP
  allocParticleData((int)Ntot, &ibody, &idx, &pos, &acc,
#ifdef  BLOCK_TIME_STEP
		    &vel, &ti
#else///BLOCK_TIME_STEP
		    &vx, &vy, &vz
#endif//BLOCK_TIME_STEP
		    );
#endif//USE_HDF5_FORMAT


  int num = (int)ceil((double)Ntot / (double)ncrit);
  double *rad;  rad = (double *)malloc(sizeof(double) * (2 + num));  if( rad == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }
  double *pot;  pot = (double *)malloc(sizeof(double) * (2 + num));  if( pot == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }
  double *bp;  bp = (double *)malloc(sizeof(double) * (2 + num));  if( bp == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }
  double *y2;  y2 = (double *)malloc(sizeof(double) * (2 + num));  if( y2 == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }


  int kind = 1;
  int skind = 1;
  FILE *fp;
  char filename[256];
  sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

  fscanf(fp, "%*d");/**< skip reading unit */
  fscanf(fp, "%d\t%d", &kind, &skind);

  ulong *gnum;  gnum = (ulong *)malloc(sizeof(ulong) * kind);  if( gnum == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }
  for(int ii = 0; ii < kind; ii++)
    fscanf(fp, "%zu", &gnum[ii]);

  fclose(fp);

  ulong *ghead;  ghead = (ulong *)malloc(sizeof(ulong) * kind);  if( ghead == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }
  ulong *gtail;  gtail = (ulong *)malloc(sizeof(ulong) * kind);  if( gtail == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }

  ghead[0] = 0;
  for(int ii = 1; ii < kind; ii++){
    const ulong idx = ghead[ii - 1] + gnum[ii - 1];
    ghead[ii    ] = idx;
    gtail[ii - 1] = idx;
  }/* for(int ii = 1; ii < kind; ii++){ */
  gtail[kind - 1] = Ntot;


  const int Nmax = num + kind * Nminimum;
  double *rr;  rr = (double *)malloc(sizeof(double) * Nmax);  if( rr == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }
  double *Jr;  Jr = (double *)malloc(sizeof(double) * Nmax);  if( Jr == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }
  double *Jt;  Jt = (double *)malloc(sizeof(double) * Nmax);  if( Jt == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }
  double *Jp;  Jp = (double *)malloc(sizeof(double) * Nmax);  if( Jp == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }

  int *ahead;  ahead = (int *)malloc(sizeof(int) * kind);  if( ahead == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }
  int *anum ;  anum  = (int *)malloc(sizeof(int) * kind);  if( anum  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");  }
  for(int ii = 0; ii < kind; ii++){
    ahead[ii] = 0;
    anum [ii] = 0;
  }/* for(int ii = 0; ii < kind; ii++){ */


  static double action2astro;
  action2astro = length2astro * velocity2astro;


  /** read particle distribution and analyze */
  int ifile = (int)(start + mpi.rank);
  for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){
    double time;
    ulong steps;
    int unit_read;
#ifdef  USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file, (uint)filenum, &hdf5, hdf5type);
    for(int ii = 0; ii < (int)Ntot; ii++){
      body[ii].x   = hdf5.pos[ii * 3    ];
      body[ii].y   = hdf5.pos[ii * 3 + 1];
      body[ii].z   = hdf5.pos[ii * 3 + 2];
      body[ii].vx  = hdf5.vel[ii * 3    ];
      body[ii].vy  = hdf5.vel[ii * 3 + 1];
      body[ii].vz  = hdf5.vel[ii * 3 + 2];
      body[ii].m   = hdf5.m[ii];
      body[ii].pot = hdf5.pot[ii];
      body[ii].idx = hdf5.idx[ii];
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#else///USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file, ibody, (uint)filenum);
    for(int ii = 0; ii < (int)Ntot; ii++){
      body[ii].x   = ibody.pos[ii].x;
      body[ii].y   = ibody.pos[ii].y;
      body[ii].z   = ibody.pos[ii].z;
      body[ii].m   = ibody.pos[ii].m;
      body[ii].vx  = ibody.vel[ii].x;
      body[ii].vy  = ibody.vel[ii].y;
      body[ii].vz  = ibody.vel[ii].z;
      body[ii].pot = ibody.acc[ii].pot;
      body[ii].idx = ibody.idx[ii];
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#endif//USE_HDF5_FORMAT
    if( unit_read != unit ){
      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
    }/* if( unit_read != unit ){ */

    for(int ii = 0; ii < (int)Ntot; ii++){
      int gid = 0;
      for(int jj = 0; jj < kind; jj++){
	if( body[ii].idx < gtail[jj] ){
	  gid = jj;
	  break;
	}/* if( body[ii].idx < gtail ){ */
      }/* for(int jj = 0; jj < kind; jj++){ */
      body[ii].gid = gid;
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */

    /** find center-of-mass and bulk velocity */
    double com[3], vel[3];
    getCenter((int)Ntot, body, com, vel);


    /** smoothing of gravitational potential profile */
    const double inv_ncrit = 1.0 / (double)ncrit;
    for(int ii = 0; ii < num; ii++){
      const int head =  ii      * ncrit;
      const int tail = (ii + 1) * ncrit;

      const int midIdx = head + (ncrit >> 1);
      rad[1 + ii] = CAST_R2D((ncrit % 1) ? body[midIdx].rad : (HALF * (body[midIdx - 1].rad + body[midIdx].rad)));/**< median of radius */

      double sum = 0.0;
      for(int jj = head; jj < tail; jj++)
	sum += CAST_R2D(body[jj].pot);
      pot[1 + ii] = sum * inv_ncrit;
    }/* for(int ii = 0; ii < num; ii++){ */
    rad[0] = 0.0;
    pot[0] = pot[1];
    rad[1 + num] = 1.0e+3 * rad[1 + num - 1];
    pot[1 + num] = 0.0;


    /** cubic spline interpolation */
    genCubicSpline1D(num + 2, rad, pot, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, y2);



    /** evaluate angle-action variables for each component in spherical gravitational potential */
    qsort(body, Ntot, sizeof(nbody_particle), gidAscendingOrder);
#if 0
    for(ulong ii = 0; ii < Ntot; ii++)
      printf("%d: %zu\n", body[ii].gid, body[ii].idx);
#endif

    for(int kk = HEAD_GROUP_IDX; kk < kind; kk++){

      if( kk >= skind )
	for(ulong ii = ghead[kk]; ii < gtail[kk]; ii++)
	  body[ii].rad = body[ii].hor;

      qsort(&body[ghead[kk]], gnum[kk], sizeof(nbody_particle), radAscendingOrder);

      const ulong ncrit_analysis = ((gnum[kk] / (ulong)ncrit) > (ulong)Nminimum) ? ((ulong)ncrit) : (gnum[kk] / (ulong)Nminimum);
      const double inv_ncrit_analysis = 1.0 / (double)ncrit_analysis;

      int idx = ahead[kk];
      for(ulong head = ghead[kk]; head < gtail[kk]; head += ncrit_analysis){
	const ulong tail = head + ncrit_analysis;
	if( tail >= gtail[kk] )
	  break;

	double Jr_sum = 0.0;
	double Jt_sum = 0.0;
	double Jp_sum = 0.0;

	for(ulong ii = head; ii < tail; ii++){
	  const double xx = CAST_R2D(body[ii].x) - com[0];
	  const double yy = CAST_R2D(body[ii].y) - com[1];
	  const double zz = CAST_R2D(body[ii].z) - com[2];
	  const double vx = CAST_R2D(body[ii].vx) - vel[0];
	  const double vy = CAST_R2D(body[ii].vy) - vel[1];
	  const double vz = CAST_R2D(body[ii].vz) - vel[2];

	  const double Lx = yy * vz - zz * vy;
	  const double Ly = zz * vx - xx * vz;
	  const double Lz = xx * vy - yy * vx;
	  const double L2 = Lx * Lx + Ly * Ly + Lz * Lz;

	  const double rr = sqrt(xx * xx + yy * yy + zz * zz);
#if 0
	  const double EE = fmax(CAST_R2D(body[ii].pot), getCubicSpline1D(rr, num + 2, rad, pot, y2)) + 0.5 * (vx * vx + vy * vy + vz * vz);
#else
	  const double EE = getCubicSpline1D(rr, num + 2, rad, pot, y2) + 0.5 * (vx * vx + vy * vy + vz * vz);
#endif


#if 0
	  for(int ll = 0; ll < num + 2; ll++)
	    fprintf(stdout, "%e\t%e\n", rad[ll], getEffectivePotential(rad[ll], L2, num + 2, rad, pot, y2));
	  exit(0);
#endif


	  /* find pericenter and apocenter */
#if 1
	  double f_min;
	  /* const double req = brent_minimum(&f_min, rad[0], rad[num + 1], L2, num + 2, rad, pot, y2); */
	  double req = brent_minimum(&f_min, rad[0], rad[num], L2, num + 2, rad, pot, y2);

	  /* tentative workaround for the case of f_min > EE (smearing potential table is the appropriate remedy) */
	  double rm = 1.0;
	  double rp = 1.0;
	  int try = 0;
	  while( f_min > EE ){
	    rm *= 0.99;
	    rp *= 1.01;
	    req = brent_minimum(&f_min, rm * rr, rp * rr, L2, num + 2, rad, pot, y2);
	    try++;
	    if( try > 300 ){
	      for(int ll = 0; ll < num + 2; ll++)
		fprintf(stdout, "%e\t%e\n", rad[ll], getEffectivePotential(rad[ll], L2, num + 2, rad, pot, y2));
	      __FPRINTF__(stderr, "# E_eff(r = %e) = %e, E_tot(r = %e) = %e\n", req, f_min, rr, EE);
	      __KILL__(stderr, "ERROR: iteration failed even after %d trials\n", try);
	    }/* if( try > 300 ){ */
	  }
	  /* if( EE < f_min ){ */
	  /*   for(int ll = 0; ll < num + 2; ll++) */
	  /*     fprintf(stdout, "%e\t%e\n", rad[ll], getEffectivePotential(rad[ll], L2, num + 2, rad, pot, y2)); */
	  /*   __KILL__(stderr, "# E_eff(r = %e) = %e, E_tot(r = %e) = %e\n", req, f_min, rr, EE); */
	  /* } */
	  assert(EE >= f_min);
	  /* fprintf(stdout, "# req = %e, fmin = %e\n", req, f_min); */
	  /* fflush(NULL); */
	  const double rperi = brent_root(rad[0], req, EE, L2, num + 2, rad, pot, y2);
	  const double rapo  = brent_root(req, rad[num + 1], EE, L2, num + 2, rad, pot, y2);

#else
#if 0
	  const double rperi = brent_root(1.0e-30, 1.125 * rr, EE, L2, num + 2, rad, pot, y2);
	  const double rapo  = brent_root(0.875 * rr, rad[num + 1], EE, L2, num + 2, rad, pot, y2);
#else
	  /* const double rperi = brent_root(1.0e-100, rr, EE, L2, num + 2, rad, pot, y2); */
	  /* const double rapo  = brent_root(rr, rad[num + 1], EE, L2, num + 2, rad, pot, y2); */
	  const double rperi = brent_root(1.0e-100, rr * (1.0 + DBL_EPSILON), EE, L2, num + 2, rad, pot, y2);
	  const double rapo  = brent_root(rr * (1.0 - DBL_EPSILON), rad[num + 1], EE, L2, num + 2, rad, pot, y2);
#endif
#endif
	  assert(rapo >= rperi);

	  /* calculate Jr */
	  Jr_sum += integrate_DEformula(rperi, rapo, EE, L2, num + 2, rad, pot, y2);

	  Jt_sum += sqrt(L2) - fabs(Lz);
	  Jp_sum += Lz;
	}/* for(ulong ii = head; ii < tail; ii++){ */

	const ulong midIdx = head + (ncrit_analysis >> 1);
	rr[idx] = length2astro * CAST_R2D((ncrit_analysis % 1) ? body[midIdx].rad : (HALF * (body[midIdx - 1].rad + body[midIdx].rad)));/**< median of radius */;

	Jr[idx] = action2astro * Jr_sum * inv_ncrit_analysis;
	Jt[idx] = action2astro * Jt_sum * inv_ncrit_analysis;
	Jp[idx] = action2astro * Jp_sum * inv_ncrit_analysis;

	idx++;
      }/* for(ulong head = ghead[kk]; head < gtail[kk]; head += ncrit_analysis){ */

      anum[kk] = idx - ahead[kk];
      if( kk != (kind - 1) )
	ahead[kk + 1] = idx;
    }/* for(int kk = HEAD_GROUP_IDX; kk < kind; kk++){ */


    /** print out the result */
#ifdef  USE_HDF5_FORMAT
    sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "action", filenum);
    hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* write attribute data */
    /* create the data space for the attribute */
    hsize_t attr_dims = 1;
    hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    /* write current time */
    double wtime = time * time2astro;
    attribute = H5Acreate(target, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &wtime));
    chkHDF5err(H5Aclose(attribute));
    /* write # of components */
    attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &kind));
    chkHDF5err(H5Aclose(attribute));
    /* write # of spherical components */
    attribute = H5Acreate(target, "skinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skind));
    chkHDF5err(H5Aclose(attribute));
    /* write # of skipped components */
    int skipped = HEAD_GROUP_IDX;
    attribute = H5Acreate(target, "skipped kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skipped));
    chkHDF5err(H5Aclose(attribute));

    extern const char length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
    extern const char   time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
    hid_t str4format = H5Tcopy(H5T_C_S1);
    chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
    attribute = H5Acreate(target, "length_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, str4format, length_astro_unit_name));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(target, "time_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, str4format, time_astro_unit_name));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Tclose(str4format));

    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    for(int kk = HEAD_GROUP_IDX; kk < kind; kk++){
      char grp[16];      sprintf(grp, "data%d", kk);
      hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      hsize_t dims = anum[kk];

      dataspace = H5Screate_simple(1, &dims, NULL);
      hid_t dataset;

      /* write r or R */
      if( kk < skind )
	dataset = H5Dcreate(group, "r", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      else
	dataset = H5Dcreate(group, "R", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rr[ahead[kk]]));
      chkHDF5err(H5Dclose(dataset));

      dataset = H5Dcreate(group, "Jr", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Jr[ahead[kk]]));
      chkHDF5err(H5Dclose(dataset));

      dataset = H5Dcreate(group, "Jtheta", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Jt[ahead[kk]]));
      chkHDF5err(H5Dclose(dataset));

      dataset = H5Dcreate(group, "Jphi", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Jp[ahead[kk]]));
      chkHDF5err(H5Dclose(dataset));

      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));


      /* write attribute data */
      /* create the data space for the attribute */
      hsize_t attr_dims = 1;
      dataspace = H5Screate_simple(1, &attr_dims, NULL);
      hid_t attribute;
      /* write # of data elements */
      attribute = H5Acreate(group, "number", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &anum[kk]));
      chkHDF5err(H5Aclose(attribute));
      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));

      chkHDF5err(H5Gclose(group));
    }/* for(int kk = HEAD_GROUP_IDX; kk < kind; kk++){ */

    /* close the file */
    chkHDF5err(H5Fclose(target));
#else///USE_HDF5_FORMAT
    sprintf(filename, "%s/%s.%s%.3u.txt", DATAFOLDER, file, "action", filenum);
    fp = fopen(filename, "w");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }

    fprintf(fp, "# r or R\tJ_r\tJ_theta\tJ_phi\n");
    for(int kk = HEAD_GROUP_IDX; kk < kind; kk++){
      for(ulong ii = ahead[kk]; ii < ahead[kk] + anum[kk]; ii++)
	fprintf(fp, "%e\t%e\t%e\t%e\n", rr[ii], Jr[ii], Jt[ii], Jp[ii]);

      fprintf(fp, "\n");
    }/* for(int kk = HEAD_GROUP_IDX; kk < kind; kk++){ */

    fclose(fp);
#endif//USE_HDF5_FORMAT

    ifile += (int)(interval * mpi.size);
  }/* for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){ */


#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
  freeSnapshotArray(hdf5_pos, hdf5_vel, hdf5_acc, hdf5_m, hdf5_pot, hdf5_idx);
#else///USE_HDF5_FORMAT
  freeParticleData(idx, pos, acc,
#ifdef  BLOCK_TIME_STEP
		    vel, ti
#else///BLOCK_TIME_STEP
		    vx, vy, vz
#endif//BLOCK_TIME_STEP
		    );
#endif//USE_HDF5_FORMAT
  free(body);

  free(rad);  free(pot);  free(bp);  free(y2);
  free(gnum);  free(ghead);  free(gtail);

  free(rr);  free(Jr);  free(Jt);  free(Jp);
  free(ahead);  free(anum);


  exitMPI();

  return (0);
}


void getCenter(const int num, nbody_particle *body, double com[restrict], double vel[restrict])
{
  __NOTE__("%s\n", "start");

    /** sort by particle position */

  double mtot = 0.0;
  double comx = 0.0;
  double comy = 0.0;
  double comz = 0.0;
  double velx = 0.0;
  double vely = 0.0;
  double velz = 0.0;

  bool converge = false;
  while( true ){

    for(int ii = 0; ii < num; ii++){
      const real xx = body[ii].x - CAST_D2R(comx);
      const real yy = body[ii].y - CAST_D2R(comy);
      const real zz = body[ii].z - CAST_D2R(comz);
      const real R2 = 1.0e-30f + xx * xx + yy * yy;
      const real r2 = R2 + zz * zz;
      body[ii].hor = R2 * RSQRT(R2);
      body[ii].rad = r2 * RSQRT(r2);
    }/* for(int ii = 0; ii < num; ii++){ */

    qsort(body, num, sizeof(nbody_particle), radAscendingOrder);

    if( converge )
      break;


    mtot = 0.0;
    velx = 0.0;
    vely = 0.0;
    velz = 0.0;

    double newx = 0.0;
    double newy = 0.0;
    double newz = 0.0;

    for(int ii = 0; ii < (num >> 1); ii++){
      const double xx = CAST_R2D(body[ii].x) - comx;
      const double yy = CAST_R2D(body[ii].y) - comy;
      const double zz = CAST_R2D(body[ii].z) - comz;

      const double mass = CAST_R2D(body[ii].m);
      mtot += mass;

      const double vx = CAST_R2D(body[ii].vx);
      const double vy = CAST_R2D(body[ii].vy);
      const double vz = CAST_R2D(body[ii].vz);

      newx += mass * xx;
      newy += mass * yy;
      newz += mass * zz;

      velx += mass * vx;
      vely += mass * vy;
      velz += mass * vz;
    }/* for(int ii = 0; ii < (num >> 1); ii++){ */

    mtot = 1.0 / mtot;
    newx *= mtot;
    newy *= mtot;
    newz *= mtot;
    velx *= mtot;
    vely *= mtot;
    velz *= mtot;

    const double dx = newx - comx;
    const double dy = newy - comy;
    const double dz = newz - comz;

    converge = ( dx * dx + dy * dy + dz * dz < 1.0e-6 );

    comx = newx;
    comy = newy;
    comz = newz;
  }/* while( true ){ */

  com[0] = comx;
  com[1] = comy;
  com[2] = comz;

  vel[0] = velx;
  vel[1] = vely;
  vel[2] = velz;


  __NOTE__("%s\n", "end");
}
