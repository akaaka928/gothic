/**
 * @file potdens.c
 *
 * @brief Source code for calculating potential-density pair of disk component(s)
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/04/24 (Tue)
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
#include <omp.h>

#include <gsl/gsl_integration.h>

#ifdef  USE_LIS
#include "lis.h"
#endif//USE_LIS

#include "macro.h"
#include "constants.h"

#include "profile.h"
#include "abel.h"
#include "blas.h"
#include "spline.h"
#include "magi.h"
#include "potdens.h"


extern const real newton;

extern double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
#define NMAX_GAUSS_QD_LOW (11)
#define NTBL_GAUSS_QD_LOW ((NMAX_GAUSS_QD_LOW >> 1) + (NMAX_GAUSS_QD_LOW & 1))
#define NINTBIN_LOW NMAX_GAUSS_QD_LOW
static double gsl_gaussQD_low_pos[NTBL_GAUSS_QD_LOW], gsl_gaussQD_low_weight[NTBL_GAUSS_QD_LOW];


#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#ifdef __ICC
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#     pragma warning (disable:869)
#endif//__ICC

/**
 * @fn getSmoothCutoff
 *
 * @brief Calculate complementary error function based smoother.
 *
 * @param (RR) position R
 * @param (Rc) cutoff radius
 * @param (invDelta) inverse of cutoff width
 * @return smoothing function
 */
static inline double getSmoothCutoff (const double RR, const double Rc, const double invDelta){  return (0.5 * erfc(2.0 * (RR - Rc) * invDelta));}

/**
 * @fn ColumnDensityExp
 *
 * @brief Calculate column density profile of exponential disk.
 *
 * @param (RR) position R
 * @param (invRd) inverse of scale length
 * @param (disk) physical properties of the disk component
 * @return column density at the specified radius
 */
double getColumnDensityExp   (double RR, double invRd, const disk_util disk);
double getColumnDensityExp   (double RR, double invRd, const disk_util disk){
  return (exp(-                    RR * invRd                   )   * getSmoothCutoff(RR, disk.Rcutoff, disk.invRsmooth));}

/**
 * @fn ColumnDensitySersic
 *
 * @brief Calculate column density profile of Sersic disk.
 *
 * @param (RR) position R
 * @param (invRd) inverse of scale length
 * @param (disk) physical properties of the disk component
 * @return column density at the specified radius
 */
double getColumnDensitySersic(double RR, double invRd, const disk_util disk);
double getColumnDensitySersic(double RR, double invRd, const disk_util disk){
  return (exp(-disk.sersic_b * pow(RR * invRd, disk.sersic_ninv))   * getSmoothCutoff(RR, disk.Rcutoff, disk.invRsmooth));}

/**
 * @fn ColumnDensitySpline
 *
 * @brief Calculate column density profile of machine-readable table format (using cubic spline interpolation).
 *
 * @param (RR) position R
 * @param (invRd) inverse of scale length
 * @param (disk) physical properties of the disk component
 * @return column density at the specified radius
 *
 * @sa getCubicSpline1D
 */
double getColumnDensitySpline(double RR, double invRd, const disk_util disk);
double getColumnDensitySpline(double RR, double invRd, const disk_util disk){
  return (getCubicSpline1D(RR, disk.num, disk.xx, disk.ff, disk.f2) * getSmoothCutoff(RR, disk.Rcutoff, disk.invRsmooth));}

/**
 * @var invzd2invz0
 * global constant to set hyperbolic secant profile in the vertical direction
 */
static double invzd2invz0;
/**
 * @fn setz0inv4SechProfile
 *
 * @brief Set invzd2invz0.
 */
static inline void setz0inv4SechProfile(void){  invzd2invz0 = 2.0 * acosh(sqrt(M_E));}

/**
 * @fn getVerticalDensity
 *
 * @brief Get vertical density of hyperbolic secant profile.
 *
 * @param (zz) height z
 * @param (invzd) inverse of scale height
 * @param (disk) physical properties of the disk component
 * @return volume density at the specified height
 */
double getVerticalDensity(const double  zz, const double invzd, const disk_util disk);
double getVerticalDensity(const double  zz, const double invzd, const disk_util disk){
  const double tmp = 1.0 / cosh(0.5 * zz * invzd * invzd2invz0);/**< := sech() */
  return (tmp * tmp);
}

#ifdef  __ICC
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#     pragma warning ( enable:869)
#endif//__ICC
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)


/**
 * @fn func4encSigma
 *
 * @brief Subfunction to execute Gaussian quadrature to get cumulative column density profile
 *
 * @param (RR) position R
 * @param (disk) physical properties of the disk component
 * @return R * Sigma(R)
 */
static inline double func4encSigma(const double RR, const disk_data disk){  return (RR * disk.getColumnDensity(RR, disk.invRd, disk.util));}

/**
 * @fn gaussQD4encSigma
 *
 * @brief Execute Gaussian quadrature to get cumulative column density profile
 *
 * @param (min) minimum of R
 * @param (max) maximum of R
 * @param (disk) physical properties of the disk component
 * @return integrated value
 */
static inline double gaussQD4encSigma(const double min, const double max, const disk_data disk)
{
  /** initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;

  if( NINTBIN & 1 )
    sum = gsl_gaussQD_weight[(NINTBIN >> 1)] * func4encSigma(pls + mns * gsl_gaussQD_pos[(NINTBIN >> 1)], disk);

  /** numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (func4encSigma(pls + mns * gsl_gaussQD_pos[ii], disk) + func4encSigma(pls - mns * gsl_gaussQD_pos[ii], disk));

  /** finalization */
  return (mns * sum);
}


/**
 * @fn freeDiskProfile
 *
 * @brief Deallocate memory for disk component(s)
 *
 * @param (ndisk) number of disk components
 * @param (disk) physical properties of the disk component
 * @param (hor) RR, zone center value
 * @param (ver) zz, zone center value
 * @param (node_hor) RR, node center value
 * @param (node_ver) zz, node center value
 * @param (pot) potential field
 * @param (rho0) density field
 * @param (rho1) density field
 * @param (rhoTot) superposed density field of multiple disk components
 * @param (dPhidR) dPhi/dR
 * @param (d2PhidR2) d^2Phi/dR^2
 * @param (Sigma) column density profile
 * @param (vsigz) vertical velocity dispersion profile
 * @param (enc) enclosed mass profile
 * @param (zd) scale height of the disk component(s)
 * @param (radSph) radius of spherical averaged profile
 * @param (rhoSph) density of spherical averaged profile
 * @param (encSph) enclosed mass of spherical averaged profile
 * @param (spline_xx) data points for cubic spline interpolation
 * @param (spline_ff) values of data points for cubic spline interpolation
 * @param (spline_f2) coefficients for cubic spline interpolation
 * @param (spline_bp) coefficients for cubic spline interpolation
 */
void   freeDiskProfile
(const int ndisk, disk_data  *disk,
 double  *hor, double  *ver, double  *node_hor, double  *node_ver,
 double  *pot, double  *rho0, double  *rho1, double  *rhoTot, double  *dPhidR, double  *d2PhidR2,
 double  *Sigma, double  *vsigz, double  *enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 double  *zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 double  *radSph, double  *rhoSph, double  *encSph,
 double  *spline_xx, double  *spline_ff, double  *spline_f2, double  *spline_bp)
{
  __NOTE__("%s\n", "start");

  free(disk);
  free(hor);  free(ver);  free(node_hor);  free(node_ver);
  free(pot);  free(rho0);  free(rho1);  free(dPhidR);  free(d2PhidR2);
  if( ndisk > 1 )
    free(rhoTot);
  free(Sigma);  free(vsigz);  free(enc);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  free(zd);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  free(radSph);  free(rhoSph);  free(encSph);
  free(spline_xx);  free(spline_ff);  free(spline_f2);  free(spline_bp);

  __NOTE__("%s\n", "end");
}


/**
 * @fn allocDiskProfile
 *
 * @brief Allocate memory for disk component(s)
 *
 * @param (ndisk) number of disk components
 * @return (disk) physical quantities of the disk component
 * @param (disk_cfg) physical properties of the disk component
 * @return (maxLev) maximum level of nested grid
 * @return (disk_prf) profile of the disk component
 * @param (skind) number of spherical symmetric components
 * @param (logrbin) rbin in logarithmic space
 * @param (invlogrbin) inverse of logrbin
 * @return (hor) RR, zone center value
 * @return (ver) zz, zone center value
 * @return (node_hor) RR, node center value
 * @return (node_ver) zz, node center value
 * @return (pot) potential field
 * @return (rho0) density field
 * @return (rho1) density field
 * @return (rhoTot) superposed density field of multiple disk components
 * @return (dPhidR) dPhi/dR
 * @return (d2PhidR2) d^2Phi/dR^2
 * @return (Sigma) column density profile
 * @return (vsigz) vertical velocity dispersion profile
 * @return (enc) enclosed mass profile
 * @return (zd) scale height of the disk component(s)
 * @return (radSph) radius of spherical averaged profile
 * @return (rhoSph) density of spherical averaged profile
 * @return (encSph) enclosed mass of spherical averaged profile
 * @return (spline_xx) data points for cubic spline interpolation
 * @return (spline_ff) values of data points for cubic spline interpolation
 * @return (spline_f2) coefficients for cubic spline interpolation
 * @return (spline_bp) coefficients for cubic spline interpolation
 */
void allocDiskProfile
(const int ndisk, disk_data **disk, profile_cfg *disk_cfg, int *maxLev, profile **disk_prf, const int skind, const double logrbin, const double invlogrbin,
 double **hor, double **ver, double **node_hor, double **node_ver,
 double **pot, double **rho0, double **rho1, double **rhoTot, double **dPhidR, double **d2PhidR2,
 double **Sigma, double **vsigz, double **enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 double **zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 double **radSph, double **rhoSph, double **encSph,
 double **spline_xx, double **spline_ff, double **spline_f2, double **spline_bp)
{
  __NOTE__("%s\n", "start");


  /** determine maximum level of the nested grid */
  /** configuration of coarsest grid */
  double Rmax = 0.0;
  double zmax = 0.0;
  for(int ii = 0; ii < ndisk; ii++){
    if( disk_cfg[ii].cutoff ){      if( Rmax <  disk_cfg[ii].rc                    )	Rmax = disk_cfg[ii].rc                  ;    }
    else{                           if( Rmax < (disk_cfg[ii].rs * DISK_MAX_LENGTH) )	Rmax = disk_cfg[ii].rs * DISK_MAX_LENGTH;    }
    if(                                 zmax < (disk_cfg[ii].zd * DISK_MAX_LENGTH) )    zmax = disk_cfg[ii].zd * DISK_MAX_LENGTH;
  }/* for(int ii = 0; ii < ndisk; ii++){ */

  int log2hmax;
  double maxLz, maxLR;
  if( Rmax > ldexp(zmax, NHOR_OVER_NVER) ){    log2hmax = (int)ceil(log2(DISK_MAX_SAFETY * Rmax));    maxLR = ldexp(1.0, log2hmax);    maxLz = ldexp(maxLR, -NHOR_OVER_NVER);  }
  else{                                        log2hmax = (int)ceil(log2(DISK_MAX_SAFETY * zmax));    maxLz = ldexp(1.0, log2hmax);    maxLR = ldexp(maxLz,  NHOR_OVER_NVER);  }
  const double hh = maxLz / (double)NDISKBIN_VER;
  log2hmax = (int)nearbyint(log2(hh));

  /** configuration of finest grid */
  double Rmin = DBL_MAX;
  double zmin = DBL_MAX;
  for(int ii = 0; ii < ndisk; ii++){
    if( Rmin > disk_cfg[ii].rs )      Rmin = disk_cfg[ii].rs;
    if( zmin > disk_cfg[ii].zd )      zmin = disk_cfg[ii].zd;
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  const int log2hmin = (int)floor(log2(DISK_MIN_LENGTH * fmin(Rmin, zmin)));

  /** hmin corresponds to 2^(1-Lmax) h */
  *maxLev = 1 - log2hmin + log2hmax;


  /** allocate arrays to store physical quantities */
  /** horizontal and vertical axes */
  *hor      = (double *)malloc((*maxLev) *  NDISKBIN_HOR      * sizeof(double));  if( *     hor == NULL ){    __KILL__(stderr, "ERROR: failure to allocate hor\n");  }
  *ver      = (double *)malloc((*maxLev) *  NDISKBIN_VER      * sizeof(double));  if( *     ver == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ver\n");  }
  *node_hor = (double *)malloc((*maxLev) * (NDISKBIN_HOR + 1) * sizeof(double));  if( *node_hor == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node_hor\n");  }
  *node_ver = (double *)malloc((*maxLev) * (NDISKBIN_VER + 1) * sizeof(double));  if( *node_ver == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node_ver\n");  }

  /** potential-density pair */
  *pot  = (double *)malloc(        (*maxLev) * NDISKBIN_HOR *  NDISKBIN_VER      * sizeof(double));  if( *pot  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pot\n" );  }
  *rho0 = (double *)malloc(ndisk * (*maxLev) * NDISKBIN_HOR * (NDISKBIN_VER + 1) * sizeof(double));  if( *rho0 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rho0\n");  }
  *rho1 = (double *)malloc(ndisk * (*maxLev) * NDISKBIN_HOR * (NDISKBIN_VER + 1) * sizeof(double));  if( *rho1 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rho1\n");  }
  /** multiple component disk */
  if( ndisk > 1 ){
    *rhoTot = (double *)malloc((*maxLev) * NDISKBIN_HOR * NDISKBIN_VER * sizeof(double));    if( *rhoTot == NULL ){      __KILL__(stderr, "ERROR: failure to allocate rhoTot\n");    }
  }/* if( ndisk > 1 ){ */

#ifndef USE_POTENTIAL_SCALING_SCHEME
  * dPhidR  = (double *)malloc((*maxLev) * NDISKBIN_HOR * NDISKBIN_VER * sizeof(double));
  *d2PhidR2 = (double *)malloc((*maxLev) * NDISKBIN_HOR * NDISKBIN_VER * sizeof(double));
#else///USE_POTENTIAL_SCALING_SCHEME
  * dPhidR  = (double *)malloc((*maxLev) * NDISKBIN_HOR                * sizeof(double));
  *d2PhidR2 = (double *)malloc((*maxLev) * NDISKBIN_HOR                * sizeof(double));
#endif//USE_POTENTIAL_SCALING_SCHEME
  if( * dPhidR  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dPhidR\n");  }
  if( *d2PhidR2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate d2PhidR2\n");  }

  /** horizontal profile */
  *Sigma  = (double *)malloc(ndisk * (*maxLev) *  NDISKBIN_HOR      * sizeof(double));  if( *Sigma == NULL ){    __KILL__(stderr, "ERROR: failure to allocate Sigma\n");  }
  *vsigz  = (double *)malloc(ndisk * (*maxLev) *  NDISKBIN_HOR      * sizeof(double));  if( *vsigz == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vsigz\n");  }
  *enc    = (double *)malloc(ndisk * (*maxLev) * (NDISKBIN_HOR + 1) * sizeof(double));  if( * enc  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate enc\n");  }

#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  *zd = (double *)malloc(ndisk * (*maxLev) * NDISKBIN_HOR * sizeof(double));
  if( *zd == NULL ){    __KILL__(stderr, "ERROR: failure to allocate zd\n");  }
#endif//ENABLE_VARIABLE_SCALE_HEIGHT

  /** for spherical averaged profile */
  *radSph = (double *)malloc(NDISKBIN_RAD * sizeof(double));  if( *radSph == NULL ){    __KILL__(stderr, "ERROR: failure to allocate radSph\n");  }
  *rhoSph = (double *)malloc(NDISKBIN_RAD * sizeof(double));  if( *rhoSph == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rhoSph\n");  }
  *encSph = (double *)malloc(NDISKBIN_RAD * sizeof(double));  if( *encSph == NULL ){    __KILL__(stderr, "ERROR: failure to allocate encSph\n");  }

  /** for spline fitting */
  /** note: NDISKBIN_HOR > NDISKBIN_VER */
  int Nspline = 2 + ((*maxLev) + 1) * (NDISKBIN_HOR >> 1);
  *spline_xx  = (double *)malloc(Nspline * CPUS_PER_PROCESS * sizeof(double));
  *spline_ff  = (double *)malloc(Nspline * CPUS_PER_PROCESS * sizeof(double));
  *spline_f2  = (double *)malloc(Nspline * CPUS_PER_PROCESS * sizeof(double));
  *spline_bp  = (double *)malloc(Nspline * CPUS_PER_PROCESS * sizeof(double));
  if( *spline_xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate spline_xx\n");  }
  if( *spline_ff == NULL ){    __KILL__(stderr, "ERROR: failure to allocate spline_ff\n");  }
  if( *spline_f2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate spline_f2\n");  }
  if( *spline_bp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate spline_bp\n");  }


  /** set horizontal and vertical axes */
  for(int ii = 0; ii < *maxLev; ii++){
    const double bin = ldexp(hh, -ii);

    (*node_hor)[INDEX2D(*maxLev, NDISKBIN_HOR + 1, ii, 0)] = 0.0;
#pragma omp parallel for
    for(int jj = 1; jj < 1 + NDISKBIN_HOR; jj++)      (*node_hor)[INDEX2D(*maxLev, 1 + NDISKBIN_HOR, ii, jj)] = bin *        (double)jj;
#pragma omp parallel for
    for(int jj = 0; jj <     NDISKBIN_HOR; jj++)      (*     hor)[INDEX2D(*maxLev,     NDISKBIN_HOR, ii, jj)] = bin * (0.5 + (double)jj);
    (*node_ver)[INDEX2D(*maxLev, NDISKBIN_VER + 1, ii, 0)] = 0.0;
#pragma omp parallel for
    for(int jj = 1; jj < 1 + NDISKBIN_VER; jj++)      (*node_ver)[INDEX2D(*maxLev, 1 + NDISKBIN_VER, ii, jj)] = bin *        (double)jj;
#pragma omp parallel for
    for(int jj = 0; jj <     NDISKBIN_VER; jj++)      (*     ver)[INDEX2D(*maxLev,     NDISKBIN_VER, ii, jj)] = bin * (0.5 + (double)jj);
  }/* for(int ii = 0; ii < *maxLev; ii++){ */


  /** allocate utility structure and commit arrays */
  *disk = (disk_data *)malloc(ndisk * sizeof(disk_data));  if( *disk == NULL ){    __KILL__(stderr, "ERROR: failure to allocate disk\n");  }

  for(int ii = 0; ii < ndisk; ii++){
    /** commit disk properties */
    (*disk)[ii].cfg   = &disk_cfg[ii];
    (*disk)[ii].Rmax  = maxLR;
    (*disk)[ii].zmax  = maxLz;
    (*disk)[ii].hh    = hh;
    (*disk)[ii].invRd = 1.0 / disk_cfg[ii].rs;
    (*disk)[ii].prf   = disk_prf[skind + ii];
    (*disk)[ii].logrbin = logrbin;
    (*disk)[ii].invlogrbin = invlogrbin;

    /** settings for Sersic profile */
    if( disk_cfg[ii].kind == SERSIC ){
      (*disk)[ii].util.sersic_ninv = 1.0 / disk_cfg[ii].n_sersic;
      (*disk)[ii].util.sersic_b    =       disk_cfg[ii].b_sersic;
    }/* if( disk_cfg[ii].kind == SERSIC ){ */

    /** additional setting about density cutoff in the horizontal direction */
    (*disk)[ii].util.Rcutoff    = Rmax * 10.0;/**< i.e., the point at infinity */
    (*disk)[ii].util.invRsmooth = (*disk)[ii].invRd;
    if( disk_cfg[ii].cutoff ){
      (*disk)[ii].util.Rcutoff    =       disk_cfg[ii].rc;
      (*disk)[ii].util.invRsmooth = 1.0 / disk_cfg[ii].rc_width;
    }/* if( disk_cfg[ii].cutoff ){ */

    /** initialize Toomre's Q-value analysis */
    (*disk)[ii].cfg->vcirc_max   = -1.0;
    (*disk)[ii].cfg->vcirc_max_R = -1.0;
    (*disk)[ii].cfg->Qmin0 = DBL_MAX;
    (*disk)[ii].cfg->Qmin1 = DBL_MAX;
    (*disk)[ii].cfg->Qmin2 = DBL_MAX;
    (*disk)[ii].cfg->qminR0 = -1.0;
    (*disk)[ii].cfg->qminR1 = -1.0;
    (*disk)[ii].cfg->qminR2 = -1.0;
    (*disk)[ii].cfg->passed = false;

    /** set function pointer */
    switch( disk_cfg[ii].kind ){
    case EXP_DISK:      (*disk)[ii].getColumnDensity = getColumnDensityExp   ;      break;
    case   SERSIC:      (*disk)[ii].getColumnDensity = getColumnDensitySersic;      break;
    case TBL_DISK:      (*disk)[ii].getColumnDensity = getColumnDensitySpline;
      readColumnDensityTable4Disk((*disk)[ii].prf, disk_cfg[ii].rs, disk_cfg[ii].file, &((*disk)[ii].util.num), &((*disk)[ii].util.xx), &((*disk)[ii].util.ff), &((*disk)[ii].util.f2), &((*disk)[ii].util.bp));
      break;
    default:      __KILL__(stderr, "ERROR: undefined model(%d) is specified as disk profile\n", disk_cfg[ii].kind);      break;
    }/* switch( disk_cfg[ii].kind ){ */

    /** set Sigma0 */
#ifdef  NDIVIDE_GAUSSQD4DISK
    const double Rbin = Rmax / (double)NDIVIDE_GAUSSQD4DISK;
    double sum = 0.0;
    double Rin = 0.0;
    for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){
      sum += gaussQD4encSigma(Rin, Rin + Rbin, (*disk)[ii]);
      Rin += Rbin;
    }/* for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){ */
    disk_cfg[ii].Sigma0 = disk_cfg[ii].Mtot / (2.0 * M_PI * sum);
#else///NDIVIDE_GAUSSQD4DISK
    disk_cfg[ii].Sigma0 = disk_cfg[ii].Mtot / (2.0 * M_PI * gaussQD4encSigma(0.0, Rmax, (*disk)[ii]));
#endif//NDIVIDE_GAUSSQD4DISK


    /** common arrays for all components */
    (*disk)[ii].hor = *hor;
    (*disk)[ii].ver = *ver;
    (*disk)[ii].node_hor = *node_hor;
    (*disk)[ii].node_ver = *node_ver;
    (*disk)[ii].pot = *pot;

    (*disk)[ii]. dPhidR  = * dPhidR ;
    (*disk)[ii].d2PhidR2 = *d2PhidR2;

    /** spherical average profile */
    (*disk)[ii].radSph = *radSph;
    (*disk)[ii].rhoSph = *rhoSph;
    (*disk)[ii].encSph = *encSph;

    /** spline fitting */
    (*disk)[ii].spline_xx = *spline_xx;
    (*disk)[ii].spline_ff = *spline_ff;
    (*disk)[ii].spline_f2 = *spline_f2;
    (*disk)[ii].spline_bp = *spline_bp;


    /** individual arrays for each component */
    (*disk)[ii].rho0 = &((*rho0)[INDEX(ndisk, *maxLev, NDISKBIN_HOR * (NDISKBIN_VER + 1), ii, 0, 0)]);    (*disk)[ii].rho    = &(*disk)[ii].rho0;
    (*disk)[ii].rho1 = &((*rho1)[INDEX(ndisk, *maxLev, NDISKBIN_HOR * (NDISKBIN_VER + 1), ii, 0, 0)]);    (*disk)[ii].rhoSum = &(*disk)[ii].rho1;

    (*disk)[ii].Sigma  = &((*Sigma)[INDEX(ndisk, *maxLev, NDISKBIN_HOR    , ii, 0, 0)]);
    (*disk)[ii].sigmaz = &((*vsigz)[INDEX(ndisk, *maxLev, NDISKBIN_HOR    , ii, 0, 0)]);
    (*disk)[ii].enc    = &((* enc )[INDEX(ndisk, *maxLev, NDISKBIN_HOR + 1, ii, 0, 0)]);

#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
    (*disk)[ii].zd = &((*zd)[INDEX(ndisk, *maxLev, NDISKBIN_HOR, ii, 0, 0)]);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT


    /** additional setting in case of multiple component disk */
    if( ndisk > 1 )
      (*disk)[ii].rhoTot = *rhoTot;


    /** initialization as a first touch */
    for(int ll = 0; ll < *maxLev; ll++)
#pragma omp parallel for
      for(int jj = 0; jj < NDISKBIN_HOR; jj++){
	(*disk)[ii].Sigma [INDEX2D(maxLev, NDISKBIN_HOR    , ll, jj)] = 0.0;
	(*disk)[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR    , ll, jj)] = 0.0;
	(*disk)[ii].enc   [INDEX2D(maxLev, NDISKBIN_HOR + 1, ll, jj)] = 0.0;
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
	(*disk)[ii].zd    [INDEX2D(maxLev, NDISKBIN_HOR, ll, jj)] = 0.0;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
	for(int kk = 0; kk < NDISKBIN_VER; kk++){
	  (*disk)[ii].rho0[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
	  (*disk)[ii].rho1[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
	}/* for(int kk = 0; kk < NDISKBIN_VER; kk++){ */
      }/* for(int jj = 0; jj < NDISKBIN_HOR; jj++){ */
  }/* for(int ii = 0; ii < ndisk; ii++){ */


  /** initialization as a first touch */
  for(int ll = 0; ll < *maxLev; ll++)
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_HOR; jj++)
      for(int kk = 0; kk < NDISKBIN_VER; kk++){
	(*pot)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
#ifndef USE_POTENTIAL_SCALING_SCHEME
	(* dPhidR )[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
	(*d2PhidR2)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
#endif//USE_POTENTIAL_SCALING_SCHEME
      }/* for(int kk = 0; kk < NDISKBIN_VER; kk++){ */

#pragma omp parallel num_threads(CPUS_PER_PROCESS)
  {
    const int tidx = omp_get_thread_num();
    const int num = 2 + ((*maxLev) + 1) * (NDISKBIN_HOR >> 1);
    for(int ii = 0; ii < num; ii++){
      (*spline_ff)[ii + tidx * num] = 0.0;
      (*spline_f2)[ii + tidx * num] = 0.0;
      (*spline_bp)[ii + tidx * num] = 0.0;
    }/* for(int ii = 0; ii < num; ii++){ */
  }

  if( ndisk > 1 )
    for(int ll = 0; ll < *maxLev; ll++)
#pragma omp parallel for
      for(int jj = 0; jj < NDISKBIN_HOR; jj++)
	for(int kk = 0; kk < NDISKBIN_VER; kk++)
	  (*rhoTot)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;


  __NOTE__("%s\n", "end");
}


/**
 * @fn bisection
 *
 * @brief Execute bisection.
 *
 * @param (val) the target value
 * @param (num) number of data points
 * @param (tab) array contains data points
 * @param (logrbl) true when data points are sampled in logarithmic space
 * @param (invbin) inverse of interval between grid points
 * @return (ratio) parameter for linear interpolation
 * @return lower index of the corresponding data point
 */
static inline int bisection(const double val, const int num, double * restrict tab, const bool logtbl, const double invbin, double * restrict ratio)
{
  int ll =       0;
  int rr = num - 1;

  /** prohibit extraporation */
  if( val < tab[ll] + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( val > tab[rr] - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }

  while( true ){
    const uint cc = ((uint)(ll + rr)) >> 1;

    if( (tab[cc] - val) * (tab[ll] - val) <= 0.0)      rr = (int)cc;
    else                                               ll = (int)cc;

    if( (1 + ll) == rr ){
      *ratio = (logtbl ? (log10(val / tab[ll])) : (val - tab[ll])) * invbin;
      return (ll);
    }/* if( (1 + ll) == rr ){ */
  }/* while( true ){ */
}


/**
 * @fn findIdxSpherical
 *
 * @brief Find a data element in the given array corresponding to the given value.
 *
 * @param (rad) radius
 * @param (prf) radial profile of the component
 * @param (invlogbin) inverse of logarithmic interval between grid points
 * @return (ratio) parameter for linear interpolation
 * @return the corresponding index to the given radius
 */
static inline int findIdxSpherical(const double rad, profile *prf, const double invlogbin, double *ratio)
{
  int ll =           0;
  int rr = NRADBIN - 1;

  if( rad < prf[ll].rad + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( rad > prf[rr].rad - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }

  while( true ){
    const uint cc = ((uint)ll + (uint)rr) >> 1;

    if( (prf[cc].rad - rad) * (prf[ll].rad - rad) <= 0.0 )      rr = (int)cc;
    else                                                        ll = (int)cc;

    if( (1 + ll) == rr ){
      *ratio = log10(rad / prf[ll].rad) * invlogbin;
      return (ll);
    }/* if( (1 + ll) == rr ){ */
  }/* while( true ){ */
}


/**
 * @fn Psi_spherical
 *
 * @brief Get relative potential of spherical symmetric components.
 *
 * @param (rad) radius
 * @param (sph) radial profile of the component
 * @param (invlogbin_sph) inverse of logarithmic interval between grid points
 * @return potential at r = rad
 *
 * @sa findIdxSpherical
 */
static inline double Psi_spherical(const double rad, profile *sph, const double invlogrbin_sph)
{
  double ratio;
  const int idx = findIdxSpherical(rad, sph, invlogrbin_sph, &ratio);

  return ((1.0 - ratio) * sph[idx].psi_tot + ratio * sph[1 + idx].psi_tot);
}


#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT

/**
 * @fn findIdxSphericalPsi
 *
 * @brief Find a data element in the given array corresponding to the given value.
 *
 * @param (psi) relative potential
 * @param (prf) radial profile of the component
 * @return (ratio) parameter for linear interpolation
 * @return the corresponding index to the given potential
 */
static inline int findIdxSphericalPsi(const double psi, profile *prf, double *ratio)
{
  int ll =           0;
  int rr = NRADBIN - 1;

  if( psi > prf[ll].psi_tot - DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( psi < prf[rr].psi_tot + DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }

  while( true ){
    const uint cc = ((uint)ll + (uint)rr) >> 1;

    if( (prf[cc].psi_tot - psi) * (prf[ll].psi_tot - psi) <= 0.0 )      rr = (int)cc;
    else                                                                ll = (int)cc;

    if( (1 + ll) == rr ){
      *ratio = (psi - prf[ll].psi_tot) / (prf[rr].psi_tot - prf[ll].psi_tot);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
  }/* while( true ){ */
}


/**
 * @fn bisection4nestedGrid
 *
 * @brief Execute bisection for nested grid.
 *
 * @param (val) the target value
 * @param (num) number of data points
 * @param (tab) array contains data points
 * @param (logrbl) true when data points are sampled in logarithmic space
 * @param (invbin) inverse of interval between grid points
 * @return (ratio) parameter for linear interpolation
 * @param (maxLev) maximum level of nested grid
 * @return (lev) the corresponding level of nested grid
 * @param (tab_lev) array to determin lev
 * @return lower index of the corresponding data point
 *
 * @sa bisection
 */
static inline int bisection4nestedGrid
(const double val, const int num, double * restrict tab, const bool logtbl, const double invbin, double * restrict ratio, const int maxLev, int * restrict lev, double * restrict tab_lev)
{
  /** 1. find the nested level */
  *lev = 0;
  for(int ii = maxLev - 1; ii >= 0; ii--)
    if( val <= tab_lev[ii] ){
      *lev = ii;
      break;
    }/* if( val <= tab_lev[ii] ){ */

  /** 2. find the index */
  int idx = bisection(val, num, &(tab[INDEX2D(maxLev, num, *lev, 0)]), logtbl, invbin, ratio);

  return (idx);
}


static inline double getPsi
(const double zz, const int maxLev, double * restrict node_ver, double * restrict tab_lev,
 const double RR, const double R2, const int lev, const int ii, double * restrict hor, double * restrict ver,
 double * restrict Phi, profile * restrict sph, const double invlogrbin_sph)
{
  int lev_z = 0;
  double azz;
  int jzz = bisection4nestedGrid(zz, NDISKBIN_VER + 1, node_ver, false, 1.0, &azz, maxLev, &lev_z, tab_lev);
  azz /= (node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, 1 + jzz)] - node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, jzz)]);

  int iiR = ii;
  double aaR = 0;
  if( lev > lev_z ){
    /** reset iiR and aaR in coarser grid */
    iiR = bisection(RR, NDISKBIN_HOR, &hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_z, 0)], false, 1.0, &aaR);
    aaR /= (hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_z, 1 + iiR)] - hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_z, iiR)]);
    /* aaR /= (enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev_z, 1 + iiR)] - enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev_z, iiR)]); */
  }/* if( lev > lev_z ){ */
  if( lev < lev_z ){
    /* reset jzz and azz in coarser grid */
    lev_z = lev;
    jzz = bisection(zz, NDISKBIN_VER + 1, &(node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, 0)]), false, 1.0, &azz);
    /* azz /= (node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, 1 + jzz)] - node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev_z, jzz)]); */
  }/* if( lev < lev_z ){ */


  /** correction about ``jzz'' and ``azz'' for zone centeric coordinate from edge coordinate */
  if( zz >= ver[INDEX2D(maxLev, NDISKBIN_VER, lev_z, 0)] ){
    if( zz <= ver[INDEX2D(maxLev, NDISKBIN_VER, lev_z, NDISKBIN_VER - 1)] ){
      jzz = bisection(zz, NDISKBIN_VER, &ver[INDEX2D(maxLev, NDISKBIN_VER, lev_z, 0)], false, 1.0, &azz);
      azz /= (ver[INDEX2D(maxLev, NDISKBIN_VER, lev_z, 1 + jzz)] - ver[INDEX2D(maxLev, NDISKBIN_VER, lev_z, jzz)]);
    }/* if( zz <= zone_ver[INDEX2D(maxLev, NDISKBIN_VER, lev, NDISKBIN_VER - 1)] ){ */
    else{
      jzz = NDISKBIN_VER - 2;
      azz = 1.0;
    }/* else{ */
  }/* if( zz >= zone_ver[INDEX2D(maxLev, NDISKBIN_VER, lev, 0)] ){ */
  else{
    jzz = 0;
    azz = 0.0;
  }/* else{ */


  const double Phi_disk =
    ((1.0 - azz) * Phi[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_z,	    iiR, jzz)] + azz * Phi[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_z,     iiR, 1 + jzz)]) * (1.0 - aaR) +
    ((1.0 - azz) * Phi[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_z, 1 + iiR, jzz)] + azz * Phi[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_z, 1 + iiR, 1 + jzz)]) *        aaR;


  return (-Phi_disk + Psi_spherical(sqrt(R2 + zz * zz), sph, invlogrbin_sph));
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
 * @fn brent4potential
 *
 * @brief Find the scale height of the disk component.
 *
 * @param (Psi) relative potential at the target z
 * @param (zd0) scale height as an input value
 * @param (PsiMid) relative potential at the mid-plane
 * @param (ver) z-axis of the disk table
 * @param (Phi) potential table of the disk component(s)
 * @param (sph) profile of the spherical symmetric component(s)
 * @param (invlogrbin_sph) inverse of interval of data points in logarithmic space for spherical symmetric component(s)
 * @param (maxLev) maximum level of nested grid
 * @param (RR) R-position of the specified location
 * @param (R2) R squared
 * @param (lev) the corresponding level of nested grid
 * @param (ii) index corresponding RR in the current level of nested grid
 * @param (hor) z-axis of the disk table
 * @param (node_ver) z-axis of the disk table, as position of boundaries of cells
 * @param (tab_lev) array to determin lev
 * @return corresponding z
 */
static inline double brent4potential
(const double Psi, const double zd0, const double PsiMid, double * restrict ver, double * restrict Phi, profile * restrict sph, const double invlogrbin_sph,
 const int maxLev, const double RR, const double R2, const int lev, const int ii, double * restrict hor, double * restrict node_ver, double * restrict tab_lev)
{
  /** evaluate potential @ z = zd0 */
  const double Psi_zd = getPsi(zd0, maxLev, node_ver, tab_lev, RR, R2, lev, ii, hor, ver, Phi, sph, invlogrbin_sph);

  /** if zd0 <= zdim, then return zd0 */
  if( Psi_zd >= Psi )
    return zd0;

  const int Niter_max = 100;
  const double tol = 1.0e-3;

  /** find zdim in (0, zd0) using Brent's method */
  double dd = 0.0;/**< displacement in the previous step */
  double ee = 0.0;/**< displacement in the step prior to the previous one*/
  double aa = 0.0;  double fa = PsiMid - Psi;
  double bb = zd0;  double fb = Psi_zd - Psi;
  double cc =  bb;  double fc = fb;

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

    if( fabs(dd) > tol1 )
      bb += dd;
    else
      bb += sign(tol1, mm);

    fb = getPsi(bb, maxLev, node_ver, tab_lev, RR, R2, lev, ii, hor, ver, Phi, sph, invlogrbin_sph) - Psi;

    iter++;
    if( iter > Niter_max ){
      __KILL__(stderr, "ERROR: Brent's method was not converged in %d steps.\n", iter);
    }/* if( iter > Niter_max ){ */
  }/* while( true ){ */
}


/**
 * @fn getVariableDiskScaleHeight
 *
 * @brief Reset scale height of disk component(s) to remove the needle-like structure.
 *
 * @param (ndisk) number of disk components
 * @param (maxLev) maximum level of nested grid
 * @return (disk) physical quantities of the disk component(s)
 * @param (prf) profile of the spherical symmetric component(s)
 * @param (invlogrbin_sph) inverse of interval of data points in logarithmic space for spherical symmetric component(s)
 *
 * @sa Phi_spherical
 * @sa findIdxSphericalPsi
 */
static inline void getVariableDiskScaleHeight(const int ndisk, const int maxLev, const int lev, disk_data * restrict disk, profile * restrict sph, const double invlogrbin_sph)
{
  __NOTE__("%s\n", "start");


  double *Phi;
  Phi = disk[0].pot;

  double *ver;
  ver = disk[0].ver;

  double *hor;
  hor = disk[0].hor;

  double *node_ver;
  node_ver = disk[0].node_ver;

  double *tab_lev;
  tab_lev = (double *)malloc(sizeof(double) * maxLev);
  if( tab_lev == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tab_lev\n");  }
  for(int ll = 0; ll < maxLev; ll++)
    tab_lev[ll] = node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, ll, NDISKBIN_VER)];


  for(int kk = 0; kk < ndisk; kk++){
    const double zd0 = disk[kk].cfg->zd;
#ifdef  ADDITIONAL_CONDITION_FOR_SCALE_HEIGHT
    double zdim = fmin(DISK_DIMMING_HEIGHT * zd0, DISK_DIMMING_SCALE * disk[kk].cfg->rs);
#else///ADDITIONAL_CONDITION_FOR_SCALE_HEIGHT
    double zdim = DISK_DIMMING_HEIGHT * zd0;
#endif//ADDITIONAL_CONDITION_FOR_SCALE_HEIGHT
    zdim = fmin(zdim, disk[kk].cfg->rc);


    double * zd;      zd  = &(disk[kk]. zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)]);

#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      const double RR = disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
      const double R2 = RR * RR;


      /** get potential of disk component(s) @ the reference points (mid plane and dimming scale) */



      /** get potential @ the reference points (mid plane and dimming scale) */
      const double PsiMid = Psi_spherical(RR, sph, invlogrbin_sph) - Phi[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)];
      const double PsiDim = getPsi(zdim, maxLev, node_ver, tab_lev, RR, R2, lev, ii, hor, ver, Phi, sph, invlogrbin_sph);

      const double Psi = PsiMid + DISK_DIMMING_HEIGHT_INV * (PsiDim - PsiMid);

#if 1
      zd[ii] = brent4potential(Psi, zd0, PsiMid, ver, Phi, sph, invlogrbin_sph, maxLev, RR, R2, lev, ii, hor, node_ver, tab_lev);
#else
      double ratio;
      const int irr = findIdxSphericalPsi(Psi, sph, &ratio);
      const double rr = (1.0 - ratio) * sph[irr].rad + ratio * sph[1 + irr].rad;
      const double zd1 = sqrt(rr * rr - R2);
      zd[ii] = fmin(zd0, zd1);
#endif
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  }/* for(int kk = 0; kk < ndisk; kk++){ */

  free(tab_lev);


  __NOTE__("%s\n", "end");
}

#endif//ENABLE_VARIABLE_SCALE_HEIGHT


/**
 * @fn gaussQuad1d4Rho
 *
 * @brief Execute Gaussian quadrature to obtain density.
 *
 * @param (func) function pointer
 * @param (min) minimum of interval of integration
 * @param (max) maximum of interval of integration
 * @param (xinv) inverse of scale length
 * @param (disk) physical properties of the disk component
 * @return integrated value
 */
static inline double gaussQuad1d4Rho(double (*func)(double, double, disk_util), const double min, const double max, const double xinv, const disk_util disk)
{
  /** initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  if( NINTBIN_LOW & 1 )
    sum = gsl_gaussQD_low_weight[(NINTBIN_LOW >> 1)] * func(pls + mns * gsl_gaussQD_low_pos[(NINTBIN_LOW >> 1)], xinv, disk);

  /** numerical integration */
  for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_low_weight[ii] *
      (func(pls + mns * gsl_gaussQD_low_pos[ii], xinv, disk) +
       func(pls - mns * gsl_gaussQD_low_pos[ii], xinv, disk));

  /** finalization */
  return (sum * mns);
}


/**
 * @fn initPotentialField
 *
 * @brief Initialize the potential field of the disk component(s).
 *
 * @param (RR) radius R
 * @param (zz) height z
 * @return (Phi) the potential field
 * @param (ndisk) number of disk components
 * @param (maxLev) maximum level of nested grid
 * @param (disk) physical properties of the disk component
 */
static inline void initPotentialField(const int ndisk, const int maxLev, disk_data *disk)
{
  __NOTE__("%s\n", "start");


  /** assume potential generated by a Miyamoto & Nagai disk for initial guess to the solution of the Poisson's equation */
  for(int lev = 0; lev < maxLev; lev++){
    double * RR;    RR  = &(disk[0].hor[INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]);
    double * zz;    zz  = &(disk[0].ver[INDEX2D(maxLev,                NDISKBIN_VER, lev, 0)]);
    double *Phi;    Phi = &(disk[0].pot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, 0)]);

#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] = 0.0;

      const double R2 = RR[ii] * RR[ii];

      for(int jj = 0; jj < NDISKBIN_VER; jj++){
	const double z2 = zz[jj] * zz[jj];

	double sum = 0.0;
	for(int kk = 0; kk < ndisk; kk++){
	  const double Rs = disk[kk].cfg->rs + sqrt(z2 + disk[kk].cfg->zd * disk[kk].cfg->zd);
	  sum += disk[kk].cfg->Mtot / sqrt(R2 + Rs * Rs);
	}/* for(int kk = 0; kk < ndisk; kk++){ */
	Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] -= CAST_R2D(newton) * sum;
      }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  }/* for(int lev = 0; lev < maxLev; lev++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn setColumnDensityProfile
 *
 * @brief Initialize the density field of the disk component(s).
 *
 * @param (ndisk) number of disk components
 * @param (maxLev) maximum level of nested grid
 * @return (disk) physical properties of the disk component(s)
 * @param (sph) radial profile of the spherical symmetric component(s)
 * @param (invlogrbin_sph) inverse of interval of data points in logarithmic space for spherical symmetric component(s)
 *
 * @sa gaussQuad1d4Rho
 * @sa setz0inv4SechProfile
 * @sa getVerticalDensity
 * @sa initPotentialField
 */
static inline void setColumnDensityProfile(const int ndisk, const int maxLev, disk_data * restrict disk, profile * restrict sph, const double invlogrbin_sph)
{
  __NOTE__("%s\n", "start");



  initPotentialField(ndisk, maxLev, disk);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  for(int ll = 0; ll < maxLev; ll++)
    getVariableDiskScaleHeight(ndisk, maxLev, ll, disk, sph, invlogrbin_sph);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT

  /** set column density profile on the midplane */
  for(int kk = 0; kk < ndisk; kk++)
    for(int ll = 0; ll < maxLev; ll++){
      const double hh = ldexp(disk[kk].hh, -ll);

      const double hinv = 1.0 / hh;
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++)
	disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)] = disk[kk].cfg->Sigma0 * gaussQuad1d4Rho(disk[kk].getColumnDensity, disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)] - 0.5 * hh, disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)] + 0.5 * hh, disk[kk].invRd, disk[kk].util) * hinv;
    }/* for(int ll = 0; ll < maxLev; ll++){ */


  /** set a tentative density distribution as an initial value in the coarsest grid */
  /** initial vertical profile is assumed to be hyperbolic secant */
  setz0inv4SechProfile();
  for(int kk = 0; kk < ndisk; kk++){
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      /** set scale height @ given R */
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      const double invzd = 1.0 / disk[kk].zd[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)];
#else///ENABLE_VARIABLE_SCALE_HEIGHT
      const double invzd = 1.0 / disk[kk].cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      const double rho0 = 0.5 * disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)] * invzd;

      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	(*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, jj)] = rho0 * gaussQuad1d4Rho(getVerticalDensity, disk[kk].ver[INDEX2D(maxLev, NDISKBIN_VER, 0, jj)] - 0.5 * disk[kk].hh, disk[kk].ver[INDEX2D(maxLev, NDISKBIN_VER, 0, jj)] + 0.5 * disk[kk].hh, invzd, disk[kk].util);


      /** calculate surface density */
#if 1
      double Sigma = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, 0)];
      for(int jj = 1; jj < NDISKBIN_VER - 1; jj++)
	Sigma += (double)(2 << (jj & 1)) * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, jj)];
      Sigma += 2.0 * disk[kk].hh / 3.0;/**< 2.0 reflects plane symmetry about the equatorial plane (z = 0) */
#else
      double Sigma = 0.0;
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	Sigma += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, jj)];
      Sigma *= 2.0 * disk[kk].hh;/**< 2.0 reflects plane symmetry about the equatorial plane (z = 0) */
#endif

      /** calibrate surface density */
      const double Mscale = disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)] / (DBL_MIN + Sigma);
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	(*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, jj)] *= Mscale;
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  }/* for(int kk = 0; kk < ndisk; kk++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn coarsePotential4boundary
 *
 * @brief Set outer boundary condition of the potential field.
 *
 * @param (lev_lrs) nested level of lower resolution grids
 * @return (pot) outer boundary condition
 */
static inline void coarsePotential4boundary(const int lev_lrs, double * restrict pot)
{
  __NOTE__("%s\n", "start");


  const int lev_hrs = lev_lrs + 1;

#pragma omp parallel for
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
    const int il = ii << 1;
    const int ih = il + 1;

    const int jj = (NDISKBIN_VER >> 1) - 1;

    const int jl = jj << 1;
    const int jh = jl + 1;

    pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] =
      (pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jh)] +
       pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jh)]) * 0.25;
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */

  {
    const int ii = (NDISKBIN_HOR >> 1) - 1;

    const int il = ii << 1;
    const int ih = il + 1;

#pragma omp parallel for
    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){
      const int jl = jj << 1;
      const int jh = jl + 1;

      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] =
	(pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jh)] +
	 pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jh)]) * 0.25;
    }/* for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){ */
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn coarsePotential
 *
 * @brief Update potential field in the coarser level.
 *
 * @param (lev_lrs) nested level of lower resolution grids
 * @return (pot) potential field
 */
static inline void coarsePotential(const int lev_lrs, double * restrict pot)
{
  __NOTE__("%s\n", "start");


  const int lev_hrs = lev_lrs + 1;

#pragma omp parallel for
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
    const int il = ii << 1;
    const int ih = il + 1;

    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){
      const int jl = jj << 1;
      const int jh = jl + 1;

      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] =
	(pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jh)] +
	 pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jh)]) * 0.25;
    }/* for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){ */
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn finePotential
 *
 * @brief Update potential field in the finer level.
 *
 * @param (lev_hrs) nested level of higher resolution grids
 * @return (pot) potential field
 */
static inline void   finePotential(const int lev_hrs, double *pot)
{
  __NOTE__("%s\n", "start");


  const int lev_lrs = lev_hrs - 1;

#pragma omp parallel for
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){
      double dfdR, dfdz;
      dfdR = 0.125 * (pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii + 1, jj)] - pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, (ii > 0) ? (ii - 1) : (0), jj)]);
      dfdz = 0.125 * (pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj + 1)] - pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, (jj > 0) ? (jj - 1) : (0))]);

      const double pot00 = pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)];

      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs,      ii << 1 ,      jj << 1 )] = pot00 - dfdR - dfdz;
      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs,      ii << 1 , 1 + (jj << 1))] = pot00 - dfdR + dfdz;
      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, 1 + (ii << 1),      jj << 1 )] = pot00 + dfdR - dfdz;
      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, 1 + (ii << 1), 1 + (jj << 1))] = pot00 + dfdR + dfdz;
    }/* for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){ */
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn coarseDensity
 *
 * @brief Update density field in the coarser level.
 *
 * @param (lev_lrs) nested level of lower resolution grids
 * @param (disk) physical quantities of the disk component
 * @return (rho) density field
 */
static inline void coarseDensity(const int lev_lrs, const disk_data disk, double * restrict rho)
{
  __NOTE__("%s\n", "start");


  const int lev_hrs = lev_lrs + 1;

#pragma omp parallel for
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
    const int il = ii << 1;
    const int ih = il + 1;

    const double R0 = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_lrs, ii)];
    const double Rl = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_hrs, il)];
    const double Rh = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_hrs, ih)];
    const double Rinv = 0.25 / R0;

    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){
      const int jl = jj << 1;
      const int jh = jl + 1;

      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] =
	((rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jl)] + rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jh)]) * Rl +
	 (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jl)] + rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jh)]) * Rh) * Rinv;
    }/* for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){ */
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn fineDensity
 *
 * @brief Update density field in the finer level.
 *
 * @param (lev_hrs) nested level of higher resolution grids
 * @return (rho) density field
 */
static inline void   fineDensity(const int lev_hrs, const disk_data disk, double *rho)
{
  __NOTE__("%s\n", "start");


  const int lev_lrs = lev_hrs - 1;
  const double hh_lrs = ldexp(disk.hh, -lev_lrs);

#pragma omp parallel for
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
    const double tmp = -hh_lrs / (3.0 * disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_lrs, ii)]);
    const double inner = (double)( 2 + 12 * ii) / (double)(3 + 12 * ii);
    const double outer = (double)(10 + 12 * ii) / (double)(9 + 12 * ii);

    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){
      double dfdR, dfdz;
      dfdR = 0.125 * (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii + 1, jj)] - rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, (ii > 0) ? (ii - 1) : (0), jj)]);
      dfdz = 0.125 * (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj + 1)] - rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, (jj > 0) ? (jj - 1) : (0))]);
      const double rhoAvg = rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] + tmp * dfdR;

      double rho00 = rhoAvg - inner * dfdR - dfdz;
      double rho01 = rhoAvg - inner * dfdR + dfdz;
      double rho10 = rhoAvg + outer * dfdR - dfdz;
      double rho11 = rhoAvg + outer * dfdR + dfdz;
      if( (rho00 < 0.0) || (rho01 < 0.0) || (rho10 < 0.0) || (rho11 < 0.0) )
      	rho00 = rho01 = rho10 = rho11 = rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)];
      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs,      ii << 1 ,      jj << 1 )] = rho00;
      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs,      ii << 1 , 1 + (jj << 1))] = rho01;
      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, 1 + (ii << 1),      jj << 1 )] = rho10;
      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, 1 + (ii << 1), 1 + (jj << 1))] = rho11;
    }/* for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){ */
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn findIdxSpherical4GDpot
 *
 * @brief Find a data element in the given array corresponding to the given value.
 *
 * @param (rad) radius
 * @param (prf) radial profile of the component
 * @return (ratio) parameter for linear interpolation
 * @return the corresponding index to the given radius
 */
static inline int findIdxSpherical4GDpot(const double rad, profile *prf, double *ratio)
{
  int ll =           0;
  int rr = NRADBIN - 1;

  if( rad < prf[ll].rho + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( rad > prf[rr].rho - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }

  while( true ){
    const uint cc = ((uint)ll + (uint)rr) >> 1;

    if( (prf[cc].rho - rad) * (prf[ll].rho - rad) <= 0.0 )      rr = (int)cc;
    else                                                        ll = (int)cc;

    if( (1 + ll) == rr ){
      *ratio = (rad - prf[ll].rho) / (prf[rr].rho - prf[ll].rho);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
  }/* while( true ){ */
}


/**
 * @fn func4GDpot
 *
 * @brief Subfunction to execute Gaussian quadrature to get potential field.
 *
 * @param (RR) position R
 * @param (a2) a squared
 * @param (disk) physical quantities of the disk component
 * @return R * Sigma(R) / sqrt(R^2 - a^2)
 */
static inline double func4GDpot(const double RR, const double a2, const disk_data disk)
{
  return (RR * disk.getColumnDensity(RR, disk.invRd, disk.util) / sqrt(RR * RR - a2));
}

/**
 * @fn gaussQD4GDpot
 *
 * @brief Execute Gaussian quadrature to get potential field.
 *
 * @param (min) minimum of R
 * @param (max) maximum of R
 * @param (a2) a squared
 * @param (disk) physical properties of the disk component
 * @return integrated value
 *
 * @sa func4GDpot
 */
double gaussQD4GDpot(const double min, const double max, const double a2, const disk_data disk);
double gaussQD4GDpot(const double min, const double max, const double a2, const disk_data disk)
{
  /** initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;

  if( NINTBIN & 1 )
    sum = gsl_gaussQD_weight[(NINTBIN >> 1)] * func4GDpot(pls + mns * gsl_gaussQD_pos[(NINTBIN >> 1)], a2, disk);

  /** numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (func4GDpot(pls + mns * gsl_gaussQD_pos[ii], a2, disk) + func4GDpot(pls - mns * gsl_gaussQD_pos[ii], a2, disk));

  /** finalization */
  return (mns * sum);
}

/**
 * @fn calcGDpot
 *
 * @brief Calculate potential field.
 *
 * @param (zz) a specified height z
 * @param (Nz) number of grid points in z-direction
 * @param (ver) height z
 * @param (invdz) inverse of interval of grid points in z-direction
 * @param (zp) a parameter
 * @param (ia) a parameter
 * @param (fa) a parameter
 * @param (iR) a parameter
 * @param (fR) a parameter
 * @param (aa) a parameter
 * @param (apR2) a parameter
 * @param (amR2) a parameter
 * @param (ndisk) number of disk components
 * @param (disk) physical properties of the disk component
 * @param (invSigma) inverse of the column density
 * @return (ret) calculated value
 *
 * @sa bisection
 */
static inline void calcGDpot
(const double zz, const int Nz, double * restrict ver, const double invdz, const double zp,
 const int ia, const double fa, const int iR, const double fR, const double aa, const double apR2, const double amR2,
 const int ndisk, disk_data * restrict disk, double * restrict invSigma, double * restrict ret)
{
  const double zpls2 = (zz - zp) * (zz - zp);  const double asinp = asin(2.0 * aa / (sqrt(zpls2 + apR2) + sqrt(zpls2 + amR2)));
  const double zmns2 = (zz + zp) * (zz + zp);  const double asinm = asin(2.0 * aa / (sqrt(zmns2 + apR2) + sqrt(zmns2 + amR2)));

  double fz;
  const int jz = bisection(zp, Nz, ver, false, invdz, &fz);

  for(int ii = 0; ii < ndisk; ii++){
    const double zeta = invSigma[ii] *
      (((1.0 - fz) * (*disk[ii].rho)[INDEX2D(NR, Nz,     iR, jz)] + fz * (*disk[ii].rho)[INDEX2D(NR, Nz,     iR, 1 + jz)]) * (1.0 - fR) +
       ((1.0 - fz) * (*disk[ii].rho)[INDEX2D(NR, Nz, 1 + iR, jz)] + fz * (*disk[ii].rho)[INDEX2D(NR, Nz, 1 + iR, 1 + jz)]) *        fR   );

    const double diff = (1.0 - fa) * disk[ii].prf[ia].psi + fa * disk[ii].prf[1 + ia].psi;

    ret[ii] = zeta * diff * (asinp + asinm);
  }/* for(int ii = 0; ii < ndisk; ii++){ */
}

/**
 * @fn _gaussQuad2d4GDpot
 *
 * @brief Execute Gaussian quadrature for calculating potential field.
 *
 * @param (zz) a specified height z
 * @param (Nz) number of grid points in z-direction
 * @param (ver) height z
 * @param (invdz) inverse of interval of grid points in z-direction
 * @param (zmin) minimum of z-range
 * @param (zmax) maximum of z-range
 * @param (ia) a parameter
 * @param (fa) a parameter
 * @param (iR) a parameter
 * @param (fR) a parameter
 * @param (aa) a parameter
 * @param (apR2) a parameter
 * @param (amR2) a parameter
 * @param (ndisk) number of disk components
 * @param (disk) physical properties of the disk component
 * @return (sum) calculated value
 * @param (invSigma) inverse of the column density
 * @param (tmp) temporary array
 *
 * @sa bisection
 */
static inline void _gaussQuad2d4GDpot
(const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 const int ia, const double fa, const int iR, const double fR, const double aa, const double apR2, const double amR2,
 const int ndisk, disk_data * restrict disk, double * restrict sum, double * invSigma, double * restrict tmp)
{
  /** initialization */
  const double mns = 0.5 * (zmax - zmin);
  const double pls = 0.5 * (zmax + zmin);
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] = 0.0;

  if( NINTBIN & 1 ){
    const double ww =             gsl_gaussQD_weight[(NINTBIN >> 1)];
    const double zp = pls + mns * gsl_gaussQD_pos   [(NINTBIN >> 1)];
    calcGDpot(zz, Nz, ver, invdz, zp, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, invSigma, tmp);
    for(int ii = 0; ii < ndisk; ii++)      sum[ii] = ww * tmp[ii];
  }/* if( NINTBIN & 1 ){ */

  /** numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){
    const double ww = gsl_gaussQD_weight[ii];

    double zp = pls + mns * gsl_gaussQD_pos[ii];
    calcGDpot(zz, Nz, ver, invdz, zp, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, invSigma, tmp);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
    zp = pls - mns * gsl_gaussQD_pos[ii];
    calcGDpot(zz, Nz, ver, invdz, zp, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, invSigma, tmp);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
  }/* for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){ */

  /** finalization */
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] *= mns;
}

/**
 * @fn _setRdep4GDpot
 *
 * @brief Set R-dependence of potential field for Gaussian quadrature.
 *
 * @param (aa) a parameter
 * @return (ia) a parameter
 * @return (fa) a parameter
 * @param (ndisk) number of disk components
 * @param (disk) physical properties of the disk component
 * @return (invSigma) inverse of the column density
 * @param (RR) a specified radius R
 * @return (iR) a parameter
 * @return (fR) a parameter
 * @param (NR) number of grid points in R-direction
 * @param (hor) radius R
 * @param (invdR) inverse of interval of grid points in R-direction
 * @return (apR2) a parameter
 * @return (amR2) a parameter
 *
 * @sa bisection
 * @sa findIdxSpherical4GDpot
 */
static inline void _setRdep4GDpot
(const double aa, int * restrict ia, double * restrict fa, const int ndisk, disk_data * restrict disk, double * restrict invSigma,
 const double RR, int * restrict iR, double * restrict fR, const int NR, double * restrict hor, const double invdR, double * restrict apR2, double * restrict amR2)
{
  *apR2 = (aa + RR) * (aa + RR);
  *amR2 = (aa - RR) * (aa - RR);

  *iR = bisection(aa, NR, hor, false, invdR, fR);
  *ia = findIdxSpherical4GDpot(aa, disk[0].prf, fa);

  for(int ii = 0; ii < ndisk; ii++)
    invSigma[ii] = 1.0 / (DBL_MIN + (1.0 - (*fR)) * disk[ii].Sigma[*iR] + (*fR) * disk[ii].Sigma[1 + (*iR)]);
}

/**
 * @fn gaussQuad2d4GDpot
 *
 * @brief Execute Gaussian quadrature to calculate potential field.
 *
 * @param (RR) a specified radius R
 * @param (NR) number of grid points in R-direction
 * @param (hor) radius R
 * @param (invdR) inverse of interval of grid points in R-direction
 * @param (Rmin) minimum of R-range
 * @param (Rmax) maximum of R-range
 * @param (zz) a specified height z
 * @param (Nz) number of grid points in z-direction
 * @param (ver) height z
 * @param (invdz) inverse of interval of grid points in z-direction
 * @param (zmin) minimum of z-range
 * @param (zmax) maximum of z-range
 * @param (ndisk) number of disk components
 * @param (disk) physical properties of the disk component
 * @param (invSigma) inverse of the column density
 * @param (sub) temporary array
 * @param (tmp) temporary array
 * @return (sum) calculated results
 *
 * @sa _setRdep4GDpot
 * @sa _gaussQuad2d4GDpot
 */
static inline void gaussQuad2d4GDpot
(const double RR, const int NR, double * restrict hor, const double invdR, const double Rmin, const double Rmax,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 const int ndisk, disk_data * restrict disk, double * restrict invSigma, double * restrict sub, double * restrict tmp, double * restrict sum)
{
  /** initialization */
  const double mns = 0.5 * (Rmax - Rmin);
  const double pls = 0.5 * (Rmax + Rmin);
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] = 0.0;
  if( NINTBIN & 1 ){
    const double ww =             gsl_gaussQD_weight[(NINTBIN >> 1)];
    const double aa = pls + mns * gsl_gaussQD_pos   [(NINTBIN >> 1)];
    int ia, iR;
    double fa, fR, apR2, amR2;
    _setRdep4GDpot(aa, &ia, &fa, ndisk, disk, invSigma, RR, &iR, &fR, NR, hor, invdR, &apR2, &amR2);
    _gaussQuad2d4GDpot(zz, Nz, ver, invdz, zmin, zmax, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, tmp, invSigma, sub);
    for(int ii = 0; ii < ndisk; ii++)      sum[ii] = ww * tmp[ii];
  }/* if( NINTBIN & 1 ){ */

  /** numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){
    const double ww = gsl_gaussQD_weight[ii];

    double aa = pls + mns * gsl_gaussQD_pos[ii];
    int ia, iR;
    double fa, fR, apR2, amR2;
    _setRdep4GDpot(aa, &ia, &fa, ndisk, disk, invSigma, RR, &iR, &fR, NR, hor, invdR, &apR2, &amR2);
    _gaussQuad2d4GDpot(zz, Nz, ver, invdz, zmin, zmax, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, tmp, invSigma, sub);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];

    aa = pls - mns * gsl_gaussQD_pos[ii];
    _setRdep4GDpot(aa, &ia, &fa, ndisk, disk, invSigma, RR, &iR, &fR, NR, hor, invdR, &apR2, &amR2);
    _gaussQuad2d4GDpot(zz, Nz, ver, invdz, zmin, zmax, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, tmp, invSigma, sub);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
  }/* for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){ */

  /** finalization */
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] *= mns;
}

/**
 * @fn integrateGDpot
 *
 * @brief Calculate potential field based on formalism in Binney & Tremaine (2008).
 *
 * @param (RR) a specified radius R
 * @param (NR) number of grid points in R-direction
 * @param (hor) radius R
 * @param (invdR) inverse of interval of grid points in R-direction
 * @param (Rs) scale radius of the disk component
 * @param (Rmax) maximum of R-range
 * @param (zz) a specified height z
 * @param (Nz) number of grid points in z-direction
 * @param (ver) height z
 * @param (invdz) inverse of interval of grid points in z-direction
 * @param (zs) scale height of the disk component
 * @param (zmax) maximum of z-range
 * @param (ndisk) number of disk components
 * @param (disk) physical properties of the disk component
 * @param (invSigma) inverse of the column density
 * @param (sub) temporary array
 * @param (tmp) temporary array
 * @param (sum) temporary array
 * @return value of potential
 *
 * @sa gaussQuad2d4GDpot
 */
static inline double integrateGDpot
(const double RR, const int NR, double * restrict hor, const double invdR, const double Rs, const double Rmax,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zs, const double zmax,
 const int ndisk, disk_data * restrict disk, double * restrict invSigma, double * restrict sub, double * restrict tmp, double * restrict sum)
{
  double Phi = 0.0;

  gaussQuad2d4GDpot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, 0.0     , 2.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, 2.0 * zs, 6.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, 6.0 * zs,       zmax, ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];

  gaussQuad2d4GDpot(RR, NR, hor, invdR, 2.0 * Rs, 8.0 * Rs, zz, Nz, ver, invdz, 0.0     , 2.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 2.0 * Rs, 8.0 * Rs, zz, Nz, ver, invdz, 2.0 * zs, 6.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 2.0 * Rs, 8.0 * Rs, zz, Nz, ver, invdz, 6.0 * zs,       zmax, ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];

  gaussQuad2d4GDpot(RR, NR, hor, invdR, 8.0 * Rs, Rmax, zz, Nz, ver, invdz, 0.0     , 2.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 8.0 * Rs, Rmax, zz, Nz, ver, invdz, 2.0 * zs, 6.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 8.0 * Rs, Rmax, zz, Nz, ver, invdz, 6.0 * zs,       zmax, ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];

  return (Phi);
}


/**
 * @fn calcOuterPotential
 *
 * @brief Calculate potential outside the boundary based on formalism in Binney & Tremaine (2008).
 *
 * @param (hor) radius R
 * @param (Rmax) maximum of R-range
 * @param (ver) height z
 * @param (zmax) maximum of z-range
 * @param (hh) interval of grid points (dR = dz)
 * @param (invhh) inverse of hh
 * @param (Rs) scale radius of the disk component
 * @param (zs) scale height of the disk component
 * @return (Phi_NR) potential outside the boundary
 * @return (Phi_Nz) potential outside the boundary
 * @param (ndisk) number of disk components
 * @param (disk) physical properties of the disk component
 * @param (stock_inv) tenporary array
 * @param (stock_sub) temporary array
 * @param (stock_tmp) temporary array
 * @param (stock_sum) temporary array
 *
 * @sa integrateGDpot
 */
static inline void calcOuterPotential
(double * restrict hor, const double Rmax, double * restrict ver, const double zmax,
 const double hh, const double invhh,
 const double Rs, const double zs, double * restrict Phi_NR, double * restrict Phi_Nz,
 const int ndisk, disk_data * restrict disk, double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum
 )
{
  __NOTE__("%s\n", "start");


  /** \Phi_{  i, N_z} = potential @ R = (  i + 1/2) * R_max / N_R, z = (N_z + 1/2) * z_max / N_z; for i = 0, ..., N_R - 1 */
#pragma omp parallel for num_threads(CPUS_PER_PROCESS)
  for(int jj = 0; jj < NDISKBIN_VER; jj++){
    const double RR = hor[NDISKBIN_HOR - 1] + hh;
    const double zz = ver[jj];

    const int tidx = omp_get_thread_num();
    double *inv = &stock_inv[ndisk * tidx];
    double *sub = &stock_sub[ndisk * tidx];
    double *tmp = &stock_tmp[ndisk * tidx];
    double *sum = &stock_sum[ndisk * tidx];
    Phi_NR[jj] =  4.0 * CAST_R2D(newton) * integrateGDpot    (RR, NDISKBIN_HOR, hor, invhh, Rs, Rmax, zz, NDISKBIN_VER, ver, invhh, zs, zmax, ndisk, disk, inv, sub, tmp, sum);
  }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */

  /** \Phi_{N_R,   j} = potential @ R = (N_R + 1/2) * R_max / N_R, z = (  j + 1/2) * z_max / N_z; for j = 0, ..., N_z - 1 */
#pragma omp parallel for num_threads(CPUS_PER_PROCESS)
  for(int ii = 0; ii < NDISKBIN_HOR; ii++){
    const double zz = ver[NDISKBIN_VER - 1] + hh;
    const double RR = hor[ii];

    const int tidx = omp_get_thread_num();
    double *inv = &stock_inv[ndisk * tidx];
    double *sub = &stock_sub[ndisk * tidx];
    double *tmp = &stock_tmp[ndisk * tidx];
    double *sum = &stock_sum[ndisk * tidx];
    Phi_Nz[ii] =  4.0 * CAST_R2D(newton) * integrateGDpot    (RR, NDISKBIN_HOR, hor, invhh, Rs, Rmax, zz, NDISKBIN_VER, ver, invhh, zs, zmax, ndisk, disk, inv, sub, tmp, sum);
  }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn setResultantVector
 *
 * @brief Set external potential field as outer boundary condition.
 *
 * @param (lev) level of nested grid
 * @param (hh) interval of grid points (dR = dz)
 * @param (RR) radius R
 * @param (rho) density field
 * @param (Phi_NR) potential outside the boundary
 * @param (Phi_Nz) potential outside the boundary
 * @return (vec) resultant vector
 */
static inline void setResultantVector
(const int lev, const double hh, double * restrict RR, double * restrict rho, double * restrict Phi_NR, double * restrict Phi_Nz,
#ifdef  USE_LIS
 LIS_VECTOR vec
#else///USE_LIS
 double * restrict vec
#endif//USE_LIS
)
{
  __NOTE__("%s\n", "start");


  const double Phi0 = 8.0 * M_PI * CAST_R2D(newton) * hh;
#pragma omp parallel for
  for(int ii = 0; ii < NDISKBIN_HOR; ii++){
    const double Ri = RR[ii];

    for(int jj = 0; jj < NDISKBIN_VER; jj++){
#ifdef  USE_LIS
      lis_vector_set_value(LIS_INS_VALUE, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj), Phi0 * Ri * rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)], vec);
#else///USE_LIS
      vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] = Phi0 * Ri * rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)];
#endif//USE_LIS
    }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */

    /** boundary condition @ z-edge */
#ifdef  USE_LIS
    lis_vector_set_value(LIS_ADD_VALUE, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1), -(double)(1 + (ii << 1)) * Phi_Nz[ii], vec);
#else///USE_LIS
    vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)] -= (double)(1 + (ii << 1)) * Phi_Nz[ii];
#endif//USE_LIS
  }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */

  /** boundary condition @ R-edge */
  for(int jj = 0; jj < NDISKBIN_VER; jj++){
#ifdef  USE_LIS
    lis_vector_set_value(LIS_ADD_VALUE, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj), -(double)(NDISKBIN_HOR << 1) * Phi_NR[jj], vec);
#else///USE_LIS
    vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] -= (double)(NDISKBIN_HOR << 1) * Phi_NR[jj];
#endif//USE_LIS
  }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn setSparseMatrix
 *
 * @brief Set sparse matrix to solve the Poisson's equation.
 */
#ifdef  USE_LIS
static inline void setSparseMatrix(LIS_MATRIX mat)
#else///USE_LIS
static inline void setSparseMatrix(const crs mat)
#endif//USE_LIS
{
  __NOTE__("%s\n", "start");


  int rowIdx = 0;
#ifndef USE_LIS
  int valIdx = 0;
  mat.row[rowIdx] = valIdx;
#endif//USE_LIS


  /** inner boundary (ii = 0) */
  /** ii = 0; jj = 0; */
#ifdef  USE_LIS
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, 0), -3.0, mat);
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, 1),  1.0, mat);
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, 0),  2.0, mat);
  rowIdx++;
#else///USE_LIS
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, 0);  mat.val[valIdx] = -3.0;  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, 1);  mat.val[valIdx] =  1.0;  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, 0);  mat.val[valIdx] =  2.0;  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
#endif//USE_LIS

  /** ii = 0; 1 <= jj <= NDISKBIN_VER - 2 */
  for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){
#ifdef  USE_LIS
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj - 1),  1.0, mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj    ), -4.0, mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj + 1),  1.0, mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, jj    ),  2.0, mat);
    rowIdx++;
#else///USE_LIS
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj - 1);    mat.val[valIdx] =  1.0;    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj    );    mat.val[valIdx] = -4.0;    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj + 1);    mat.val[valIdx] =  1.0;    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, jj    );    mat.val[valIdx] =  2.0;    valIdx++;
    rowIdx++;    mat.row[rowIdx] = valIdx;
#endif//USE_LIS
  }/* for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){ */

  /** ii = 0; jj = NDISKBIN_VER - 1 */
#ifdef  USE_LIS
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, NDISKBIN_VER - 2),  1.0, mat);
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, NDISKBIN_VER - 1), -4.0, mat);
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, NDISKBIN_VER - 1),  2.0, mat);
  rowIdx++;
#else///USE_LIS
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, NDISKBIN_VER - 2);  mat.val[valIdx] =  1.0;  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, NDISKBIN_VER - 1);  mat.val[valIdx] = -4.0;  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, NDISKBIN_VER - 1);  mat.val[valIdx] =  2.0;  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
#endif//USE_LIS


  /** main domain (1 <= ii <= NDISKBIN_HOR - 2) */
  for(int ii = 1; ii < NDISKBIN_HOR - 1; ii++){
    /** inner boundary (jj = 0) */
#ifdef  USE_LIS
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, 0),  (double)(      ii << 1      ), mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , 0), -(double)((1 + (ii << 1)) * 3), mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , 1),  (double)( 1 + (ii << 1)     ), mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, 0),  (double)((1 +  ii     ) << 1), mat);
    rowIdx++;
#else///USE_LIS
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, 0);    mat.val[valIdx] =  (double)(      ii << 1      );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , 0);    mat.val[valIdx] = -(double)((1 + (ii << 1)) * 3);    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , 1);    mat.val[valIdx] =  (double)( 1 + (ii << 1)     );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, 0);    mat.val[valIdx] =  (double)((1 +  ii     ) << 1);    valIdx++;
    rowIdx++;    mat.row[rowIdx] = valIdx;
#endif//USE_LIS

    /** main domain (1 <= jj <= Nz - 2) */
    for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){
#ifdef  USE_LIS
      lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, jj    ),  (double)(      ii << 1       ), mat);
      lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , jj - 1),  (double)( 1 + (ii << 1)      ), mat);
      lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , jj    ), -(double)((1 + (ii << 1)) << 2), mat);
      lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , jj + 1),  (double)( 1 + (ii << 1)      ), mat);
      lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, jj    ),  (double)((1 +  ii      ) << 1), mat);
      rowIdx++;
#else///USE_LIS
      mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, jj    );      mat.val[valIdx] =  (double)(      ii << 1       );      valIdx++;
      mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , jj - 1);      mat.val[valIdx] =  (double)( 1 + (ii << 1)      );      valIdx++;
      mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , jj    );      mat.val[valIdx] = -(double)((1 + (ii << 1)) << 2);      valIdx++;
      mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , jj + 1);      mat.val[valIdx] =  (double)( 1 + (ii << 1)      );      valIdx++;
      mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, jj    );      mat.val[valIdx] =  (double)((1 +  ii      ) << 1);      valIdx++;
      rowIdx++;      mat.row[rowIdx] = valIdx;
#endif//USE_LIS
    }/* for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){ */

    /** outer boundary (jj = Nz - 1) */
#ifdef  USE_LIS
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, NDISKBIN_VER - 1),  (double)(      ii << 1       ), mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , NDISKBIN_VER - 2),  (double)( 1 + (ii << 1)      ), mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , NDISKBIN_VER - 1), -(double)((1 + (ii << 1)) << 2), mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, NDISKBIN_VER - 1),  (double)((1 +  ii      ) << 1), mat);
    rowIdx++;
#else///USE_LIS
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, NDISKBIN_VER - 1);    mat.val[valIdx] =  (double)(      ii << 1       );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , NDISKBIN_VER - 2);    mat.val[valIdx] =  (double)( 1 + (ii << 1)      );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , NDISKBIN_VER - 1);    mat.val[valIdx] = -(double)((1 + (ii << 1)) << 2);    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, NDISKBIN_VER - 1);    mat.val[valIdx] =  (double)((1 +  ii      ) << 1);    valIdx++;
    rowIdx++;    mat.row[rowIdx] = valIdx;
#endif//USE_LIS
  }/* for(int ii = 1; ii < NDISKBIN_HOR - 1; ii++){ */


  /** outer boundary (ii = NDISKBIN_HOR - 1) */
  /* ii = NDISKBIN_HOR - 1; jj = 0; */
#ifdef  USE_LIS
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, 0),  (double)(      (NDISKBIN_HOR - 1) << 1      ), mat);
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, 0), -(double)((1 + ((NDISKBIN_HOR - 1) << 1)) * 3), mat);
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, 1),  (double)( 1 + ((NDISKBIN_HOR - 1) << 1)     ), mat);
  rowIdx++;
#else///USE_LIS
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, 0);  mat.val[valIdx] =  (double)(      (NDISKBIN_HOR - 1) << 1      );  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, 0);  mat.val[valIdx] = -(double)((1 + ((NDISKBIN_HOR - 1) << 1)) * 3);  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, 1);  mat.val[valIdx] =  (double)( 1 + ((NDISKBIN_HOR - 1) << 1)     );  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
#endif//USE_LIS

  /** ii = NDISKBIN_HOR - 1; 1 <= jj <= NDISKBIN_VER - 2 */
  for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){
#ifdef  USE_LIS
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj    ),  (double)(      (NDISKBIN_HOR - 1) << 1       ), mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj - 1),  (double)( 1 + ((NDISKBIN_HOR - 1) << 1)      ), mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj    ), -(double)((1 + ((NDISKBIN_HOR - 1) << 1)) << 2), mat);
    lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj + 1),  (double)( 1 + ((NDISKBIN_HOR - 1) << 1)      ), mat);
    rowIdx++;
#else///USE_LIS
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj    );    mat.val[valIdx] =  (double)(      (NDISKBIN_HOR - 1) << 1       );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj - 1);    mat.val[valIdx] =  (double)( 1 + ((NDISKBIN_HOR - 1) << 1)      );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj    );    mat.val[valIdx] = -(double)((1 + ((NDISKBIN_HOR - 1) << 1)) << 2);    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj + 1);    mat.val[valIdx] =  (double)( 1 + ((NDISKBIN_HOR - 1) << 1)      );    valIdx++;
    rowIdx++;    mat.row[rowIdx] = valIdx;
#endif//USE_LIS
  }/* for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){ */

  /** ii = NR - 1; jj = NDISKBIN_VER - 1 */
#ifdef  USE_LIS
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, NDISKBIN_VER - 1),  (double)(      (NDISKBIN_HOR - 1 ) << 1       ), mat);
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, NDISKBIN_VER - 2),  (double)( 1 + ((NDISKBIN_HOR - 1 ) << 1)      ), mat);
  lis_matrix_set_value(LIS_INS_VALUE, rowIdx, INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, NDISKBIN_VER - 1), -(double)((1 + ((NDISKBIN_HOR - 1 ) << 1)) << 2), mat);
  rowIdx++;
#else///USE_LIS
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, NDISKBIN_VER - 1);  mat.val[valIdx] =  (double)(      (NDISKBIN_HOR - 1 ) << 1       );  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, NDISKBIN_VER - 2);  mat.val[valIdx] =  (double)( 1 + ((NDISKBIN_HOR - 1 ) << 1)      );  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, NDISKBIN_VER - 1);  mat.val[valIdx] = -(double)((1 + ((NDISKBIN_HOR - 1 ) << 1)) << 2);  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
#endif//USE_LIS


#ifdef  USE_LIS
  /* lis_matrix_set_type(mat, LIS_MATRIX_CSR); */
  lis_matrix_set_type(mat, LIS_MATRIX_MSR);
  lis_matrix_assemble(mat);
#else///USE_LIS
  if( mat.row[rowIdx] != (5 * NDISKBIN_HOR * NDISKBIN_VER - 2 * (NDISKBIN_HOR + NDISKBIN_VER)) ){
    __KILL__(stderr, "ERROR: Number of non-zero elements is %d, while the expected is %d\n", mat.row[rowIdx], 5 * NDISKBIN_HOR * NDISKBIN_VER - 2 * (NDISKBIN_HOR + NDISKBIN_VER));
  }/* if( mat.row[rowIdx] != (5 * NDISKBIN_HOR * NDISKBIN_VER - 2 * (NDISKBIN_HOR + NDISKBIN_VER)) ){ */
#endif//USE_LIS


  __NOTE__("%s\n", "end");
}


/**
 * @fn getPotentialField
 *
 * @brief Solve the Poisson's equation.
 *
 * @param (ndisk) number of disk components
 * @param (disk) physical properties of the disk component
 * @param (lev) level of nested grid
 * @param (hh) interval of grid points (dR = dz)
 * @param (invhh) inverse of hh
 * @param (RR) radius R
 * @param (zz) height z
 * @param (rho) density field
 * @return (Phi) potential field
 * @param (Phi_NR) potential outside the boundary
 * @param (Phi_Nz) potential outside the boundary
 * @param (mat) sparse matrix
 * @param (pre) preconditioner
 * @param (stock_inv) tenporary array
 * @param (stock_sub) temporary array
 * @param (stock_tmp) temporary array
 * @param (stock_sum) temporary array
 *
 * @sa coarsePotential4boundary
 * @sa calcOuterPotential
 * @sa setResultantVector
 * @sa pbicgstab
 */
static inline void getPotentialField
(const int ndisk, disk_data * restrict disk, const int lev,
 const double hh, const double invhh, double * restrict RR, double * restrict zz,
 double * restrict rho, double * restrict Phi, double * restrict Phi_NR, double * restrict Phi_Nz,
#ifdef  USE_LIS
 LIS_MATRIX lis_mat, LIS_VECTOR lis_b, LIS_VECTOR lis_x, LIS_SOLVER lis_solver,
#else///USE_LIS
 soaBiCGSTAB mat, soaPreConditioning pre,
#endif//USE_LIS
 double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum
 )
{
  __NOTE__("%s\n", "start");


  /** prepare Poisson equation in matrix form (vector part) */
  if( lev > 0 ){
    const int lev_lrs = lev - 1;
    coarsePotential4boundary(lev_lrs, disk[0].pot);

#if 1
    /** set z-edge[NR] */
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii += 2){
      const int il = ii >> 1;
      const int jl = (NDISKBIN_VER - 1) >> 1;

      const double mm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, (il > 0) ? (il - 1) : (0), jl)];
      const double mp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, (il > 0) ? (il - 1) : (0), jl + 1)];

      const double cm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, il, jl)];
      const double cp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, il, jl + 1)];

      const double pm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, il + 1, jl)];
      const double pp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, il + 1, jl + 1)];

      const double mz = 0.75 * mp + 0.25 * mm;
      const double cz = 0.75 * cp + 0.25 * cm;
      const double pz = 0.75 * pp + 0.25 * pm;

      Phi_Nz[ii    ] = 0.75 * cz + 0.25 * mz;
      Phi_Nz[ii + 1] = 0.75 * cz + 0.25 * pz;
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii += 2){ */


    /** set R-edge[Nz] */
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_VER; jj += 2){
      const int il = (NDISKBIN_HOR - 1) >> 1;
      const int jl = jj >> 1;

      const double mm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, il, (jl > 0) ? (jl - 1) : (0))];
      const double mc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, il, jl)];
      const double mp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, il, jl + 1)];

      const double pm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, il + 1, (jl > 0) ? (jl - 1) : (0))];
      const double pc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, il + 1, jl)];
      const double pp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, il + 1, jl + 1)];

      const double Rm = 0.75 * pm + 0.25 * mm;
      const double Rc = 0.75 * pc + 0.25 * mc;
      const double Rp = 0.75 * pp + 0.25 * mp;

      Phi_NR[jj    ] = 0.75 * Rc + 0.25 * Rm;
      Phi_NR[jj + 1] = 0.75 * Rc + 0.25 * Rp;
    }/* for(int jj = 0; jj < NDISKBIN_VER; jj += 2){ */
#else
    /** set z-edge[NR] */
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      /* const double mm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii > 1) ? ((ii >> 1) - 1) : (0), (NDISKBIN_VER >> 1) - 1)]; */
      /* const double mc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii > 1) ? ((ii >> 1) - 1) : (0),  NDISKBIN_VER >> 1)	 ]; */
      /* const double mp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii > 1) ? ((ii >> 1) - 1) : (0), (NDISKBIN_VER >> 1) + 1)]; */
      const double mm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii > 2) ? ((ii >> 1) - 1) : (0), (NDISKBIN_VER >> 1) - 1)];
      const double mc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii > 2) ? ((ii >> 1) - 1) : (0),  NDISKBIN_VER >> 1)	 ];
      const double mp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii > 2) ? ((ii >> 1) - 1) : (0), (NDISKBIN_VER >> 1) + 1)];

      const double cm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		    ii >> 1	       , (NDISKBIN_VER >> 1) - 1)];
      const double cc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		    ii >> 1	       ,  NDISKBIN_VER >> 1)	 ];
      const double cp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		    ii >> 1	       , (NDISKBIN_VER >> 1) + 1)];

      const double pm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		   (ii >> 1) + 1       , (NDISKBIN_VER >> 1) - 1)];
      const double pc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		   (ii >> 1) + 1       ,  NDISKBIN_VER >> 1)	 ];
      const double pp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		   (ii >> 1) + 1       , (NDISKBIN_VER >> 1) + 1)];

      const double zmm = mc - 0.125 * (mp - mm);
      const double zcm = cc - 0.125 * (cp - cm);
      const double zpm = pc - 0.125 * (pp - pm);

      Phi_Nz[ii] = zcm + (double)(((ii & 1) << 1) - 1) * 0.125 * (zpm - zmm);
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */

    /** set R-edge[Nz] */
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_VER; jj++){
      /* const double mm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) - 1, (jj > 1 ) ? ((jj >> 1) - 1) : (0))]; */
      /* const double cm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,	NDISKBIN_HOR >> 1     , (jj > 1 ) ? ((jj >> 1) - 1) : (0))]; */
      /* const double pm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) + 1, (jj > 1 ) ? ((jj >> 1) - 1) : (0))]; */
      const double mm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) - 1, (jj > 2 ) ? ((jj >> 1) - 1) : (0))];
      const double cm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,	NDISKBIN_HOR >> 1     , (jj > 2 ) ? ((jj >> 1) - 1) : (0))];
      const double pm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) + 1, (jj > 2 ) ? ((jj >> 1) - 1) : (0))];

      const double mc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) - 1,		      jj >> 1		 )];
      const double cc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,	NDISKBIN_HOR >> 1     ,		      jj >> 1		 )];
      const double pc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) + 1,		      jj >> 1		 )];

      const double mp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) - 1,		     (jj >> 1) + 1	 )];
      const double cp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,	NDISKBIN_HOR >> 1     ,		     (jj >> 1) + 1	 )];
      const double pp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) + 1,		     (jj >> 1) + 1	 )];

      const double Rmm = cm - 0.125 * (pm - mm);
      const double Rmc = cc - 0.125 * (pc - mc);
      const double Rmp = cp - 0.125 * (pp - mp);

      Phi_NR[jj] = Rmc + (double)(((jj & 1) << 1) - 1) * 0.125 * (Rmp - Rmm);
    }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
#endif
  }/* if( lev > 0 ){ */
  else{
    double Mtot = 0.0;
    double Rs   = 0.0;
    double zs   = 0.0;
    for(int ii = 0; ii < ndisk; ii++){
      Mtot += disk[ii].cfg->Mtot;
      Rs   += disk[ii].cfg->Mtot * disk[ii].cfg->rs;
      zs   += disk[ii].cfg->Mtot * disk[ii].cfg->zd;
    }/* for(int ii = 0; ii < ndisk; ii++){ */
    const double Minv = 1.0 / Mtot;
    Rs *= Minv;
    zs *= Minv;

    calcOuterPotential(RR, disk[0].Rmax, zz, disk[0].zmax, hh, invhh, Rs, zs, Phi_NR, Phi_Nz,
		       ndisk, disk, stock_inv, stock_sub, stock_tmp, stock_sum
#ifdef  CONFIRM_BUILDING_BLOCK
		       , Mtot
#endif//CONFIRM_BUILDING_BLOCK
		       );
  }/* else{ */

  setResultantVector(lev, hh, RR, rho, Phi_NR, Phi_Nz,
#ifdef  USE_LIS
		     lis_b
#else///USE_LIS
		     mat.vec
#endif//USE_LIS

		     );


  /** solve Poisson equation using iterative method */
#ifdef  USE_LIS
  lis_solve(lis_mat, lis_b, lis_x, lis_solver);
#pragma omp parallel for
  for(int ii = 0; ii < NCOL_CG; ii++)
    lis_vector_get_value(lis_x, ii, &Phi[ii]);
#else///USE_LIS
  pbicgstab(mat.mat, NDISKBIN_HOR * NDISKBIN_VER, mat.vec, Phi, mat.res, mat.sdw, mat.mid, mat.tmp, mat.Api, mat.Ati, pre.ilu, pre.Kri, pre.Kpi, pre.Kti, CONVERGENCE_BICGSTAB);
#endif//USE_LIS


  __NOTE__("%s\n", "end");
}


/**
 * @fn _setzdep4calcRho
 *
 * @brief Set z-dependence of density field.
 *
 * @param (zz) a specified height z
 * @param (Nz) number of grid points in z-direction
 * @param (ver) height z
 * @param (invdz) inverse of interval of grid points in z-direction
 * @param (Phi) potential field
 * @param (sph) radial profile of the spherical symmetric component
 * @param (invlogbin_sph) inverse of logarithmic interval between grid points
 * @param (iR) a parameter
 * @param (fR) a parameter
 * @param (R2) R squared
 * @param (PsiR0) Psi @ R = R, z = 0
 * @param (invPsi) 1 / (Psi(R, zd) - Psi(R, 0))
 * @return z-dependence
 *
 * @sa bisection
 * @sa Phi_spherical
 */
static inline double _setzdep4calcRho
(const double zz, const int Nz, double * restrict ver, const double invdz,
 double * restrict Phi, profile * restrict sph, const double invlogrbin_sph,
 const int iR, const double fR, const double R2, const double PsiR0, const double invPsi)
{
  double fz;
  const int jz = bisection(zz, Nz, ver, false, invdz, &fz);

  const double PsiRz =
    ((1.0 - fz) * (-Phi[INDEX2D(NR, Nz,     iR, jz)]) + fz * (-Phi[INDEX2D(NR, Nz,     iR, 1 + jz)])) * (1.0 - fR) +
    ((1.0 - fz) * (-Phi[INDEX2D(NR, Nz, 1 + iR, jz)]) + fz * (-Phi[INDEX2D(NR, Nz, 1 + iR, 1 + jz)])) *        fR  +
    Psi_spherical(sqrt(R2 + zz * zz), sph, invlogrbin_sph);

#if 0
  if( (fpclassify(exp(-(PsiRz - PsiR0) * invPsi)) != FP_NORMAL) && (fpclassify(exp(-(PsiRz - PsiR0) * invPsi)) != FP_ZERO) && (fpclassify(exp(-(PsiRz - PsiR0) * invPsi)) != FP_SUBNORMAL)){
    __FPRINTF__(stderr, "iR = %d, zz = %e, PsiRz = %e, PsiR0 = %e, invPsi = %e, exponent = %e\n", iR, zz, PsiRz, PsiR0, invPsi, -(PsiRz - PsiR0) * invPsi);
  }
#endif

  return (exp(-(PsiRz - PsiR0) * invPsi));
}

/**
 * @fn _setRdep4calcRho
 *
 * @brief Set R-dependence of density field.
 *
 * @param (RR) a specified radius R
 * @return (R2) R squared
 * @return (horDep) R-dependence
 * @return (iR) a parameter
 * @return (fR) a parameter
 * @return (PsiR0) Psi @ R = R, z = 0
 * @return (invPsi) 1 / (Psi(R, zd) - Psi(R, 0))
 * @param (NR) number of grid points in R-direction
 * @param (horRad) radius R
 * @param (invdR) inverse of interval of grid points in R-direction
 * @param (Nz) number of grid points in z-direction
 * @param (ver) height z
 * @param (invdz) inverse of interval of grid points in z-direction
 * @param (Phi) potential field
 * @param (lev) level of nested grid
 * @param (disk) physical properties of the disk component
 * @param (sph) radial profile of the spherical symmetric component
 * @param (invlogbin_sph) inverse of logarithmic interval between grid points
 *
 * @sa bisection
 * @sa Phi_spherical
 */
static inline void _setRdep4calcRho
(const double RR, double * restrict R2, double * restrict horDep, int * restrict iR, double * restrict fR, double * restrict PsiR0, double * restrict invPsi,
 const int NR, double * restrict horRad, const double invdR, const int Nz, double * restrict Phi,
 const int lev, const disk_data disk, profile * restrict sph, const double invlogrbin_sph)
{
  *R2 = RR * RR;

  *iR = bisection(RR, NR, horRad, false, invdR, fR);
#if 0
  *PsiR0 = Psi_spherical(RR, sph, invlogrbin_sph) - 0.25 * ((1.0 - (*fR)) * (5.0 * Phi[INDEX2D(NR, Nz, *iR, 0)] - Phi[INDEX2D(NR, Nz, *iR, 1)]) + (*fR) * (5.0 * Phi[INDEX2D(NR, Nz, 1 + (*iR), 0)] - Phi[INDEX2D(NR, Nz, 1 + (*iR), 1)]));
#else
  *PsiR0 = (1.0 - (*fR)) * (-Phi[INDEX2D(NR, Nz, *iR, 0)]) + (*fR) * (-Phi[INDEX2D(NR, Nz, 1 + (*iR), 0)]) + Psi_spherical(RR, sph, invlogrbin_sph);
#endif

#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  const double    zd = (1.0 - (*fR)) * disk.zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, *iR)] + (*fR) * disk.zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + (*iR))];
#else///ENABLE_VARIABLE_SCALE_HEIGHT
  const double    zd = disk.cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  const double invzd = 1.0 / zd;

  *horDep = invzd * disk.getColumnDensity(RR, disk.invRd, disk.util);

  /** find appropriate grid level */
  const int levz = (int)floor(log2(disk.hh * ((double)Nz - 0.5) * invzd));
  const int lev_zd = (levz > lev) ? lev : levz;
  /** find the corresponding grid location */
  const double hinv_zd = 1.0 / ldexp(disk.hh, -lev_zd);
  double fz;
  int jz = bisection(zd, Nz, &disk.ver[INDEX2D(maxLev, NDISKBIN_VER, lev_zd, 0)], false, hinv_zd, &fz);
  double aaR = *fR;
  int iiR = *iR;
  if( lev_zd != lev )
    iiR = bisection(RR, NR, &disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_zd, 0)], false, hinv_zd, &aaR);
  const double Psi_Rzd =
    ((1.0 - fz) * (-disk.pot[INDEX(maxLev, NR, Nz, lev_zd,     iiR, jz)]) + fz * (-disk.pot[INDEX(maxLev, NR, Nz, lev_zd,     iiR, 1 + jz)])) * (1.0 - aaR) +
    ((1.0 - fz) * (-disk.pot[INDEX(maxLev, NR, Nz, lev_zd, 1 + iiR, jz)]) + fz * (-disk.pot[INDEX(maxLev, NR, Nz, lev_zd, 1 + iiR, 1 + jz)])) *        aaR  +
    Psi_spherical(sqrt((*R2) + zd * zd), sph, invlogrbin_sph);

  *invPsi = 1.0 / (Psi_Rzd - (*PsiR0));
}


/**
 * @fn _gaussQuad2d4calcRho
 *
 * @brief Execute Gaussian quadrature to obtain density field.
 *
 * @param (zmin) minimum of z-range
 * @param (zmax) maximum of z-range
 * @param (Nz) number of grid points in z-direction
 * @param (ver) height z
 * @param (invdz) inverse of interval of grid points in z-direction
 * @param (iR) a parameter
 * @param (fR) a parameter
 * @param (R2) R squared
 * @param (PsiR0) Psi @ R = R, z = 0
 * @param (invPsi) 1 / (Psi(R, zd) - Psi(R, 0))
 * @param (Phi) potential field
 * @param (sph) radial profile of the spherical symmetric component
 * @param (invlogbin_sph) inverse of logarithmic interval between grid points
 *
 * @sa _setzdep4calcRho
 */
static inline double _gaussQuad2d4calcRho
(const double zmin, const double zmax, const int Nz, double * restrict ver, const double invdz, const int iR, const double fR, const double R2,
 const double PsiR0, const double invPsi, double * restrict Phi, profile * restrict sph, const double invlogrbin_sph)
{
  /** initialization */
  const double mns = 0.5 * (zmax - zmin);
  const double pls = 0.5 * (zmax + zmin);
  double sum = 0.0;
  if( NINTBIN_LOW & 1 ){
    sum = gsl_gaussQD_low_weight[(NINTBIN_LOW >> 1)] *
      _setzdep4calcRho(pls + mns * gsl_gaussQD_low_pos[(NINTBIN_LOW >> 1)], Nz, ver, invdz, Phi, sph, invlogrbin_sph, iR, fR, R2, PsiR0, invPsi);
  }/* if( NINTBIN_LOW & 1 ){ */

  /** numerical integration */
  for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_low_weight[ii] *
      (_setzdep4calcRho(pls + mns * gsl_gaussQD_low_pos[ii], Nz, ver, invdz, Phi, sph, invlogrbin_sph, iR, fR, R2, PsiR0, invPsi) +
       _setzdep4calcRho(pls - mns * gsl_gaussQD_low_pos[ii], Nz, ver, invdz, Phi, sph, invlogrbin_sph, iR, fR, R2, PsiR0, invPsi));

  /** finalization */
  return (sum * mns);
}

/**
 * @fn gaussQuad2d4calcRho
 *
 * @brief Execute Gaussian quadrature to obtain density field.
 *
 * @param (Rmin) minimum of R-range
 * @param (Rmax) maximum of R-range
 * @param (NR) number of grid points in R-direction
 * @param (hor) radius R
 * @param (invdR) inverse of interval of grid points in R-direction
 * @param (zmin) minimum of z-range
 * @param (zmax) maximum of z-range
 * @param (Nz) number of grid points in z-direction
 * @param (ver) height z
 * @param (invdz) inverse of interval of grid points in z-direction
 * @param (sph) radial profile of the spherical symmetric component
 * @param (invlogbin_sph) inverse of logarithmic interval between grid points
 * @param (lev) level of nested grid
 * @param (disk) physical properties of the disk component
 *
 * @sa _setRdep4calcRho
 * @sa _gaussQuad2d4calcRho
 */
static inline double gaussQuad2d4calcRho
(const double Rmin, const double Rmax, const int NR, double * restrict hor, const double invdR,
 const double zmin, const double zmax, const int Nz, double * restrict ver, const double invdz,
 profile * restrict sph, const double invlogrbin_sph, double * restrict Phi, const int lev, const disk_data disk)
{
  /** initialization */
  const double mns = 0.5 * (Rmax - Rmin);
  const double pls = 0.5 * (Rmax + Rmin);
  double sum = 0.0;
  if( NINTBIN_LOW & 1 ){
    double R2, radVal, fR, PsiR0, invPsi;
    int iR;
    _setRdep4calcRho(pls + mns * gsl_gaussQD_low_pos[(NINTBIN_LOW >> 1)],
		     &R2, &radVal, &iR, &fR, &PsiR0, &invPsi, NR, hor, invdR, Nz, Phi, lev, disk, sph, invlogrbin_sph);

    sum = gsl_gaussQD_low_weight[(NINTBIN_LOW >> 1)] * radVal * _gaussQuad2d4calcRho(zmin, zmax, Nz, ver, invdz, iR, fR, R2, PsiR0, invPsi, Phi, sph, invlogrbin_sph);
  }/* if( NINTBIN_LOW & 1 ){ */

  /** numerical integration */
  for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--){
    double R2, radVal, fR, PsiR0, invPsi;
    int iR;
    _setRdep4calcRho(pls + mns * gsl_gaussQD_low_pos[ii], &R2, &radVal, &iR, &fR, &PsiR0, &invPsi, NR, hor, invdR, Nz, Phi, lev, disk, sph, invlogrbin_sph);
    sum += gsl_gaussQD_low_weight[ii] * radVal * _gaussQuad2d4calcRho(zmin, zmax, Nz, ver, invdz, iR, fR, R2, PsiR0, invPsi, Phi, sph, invlogrbin_sph);

    _setRdep4calcRho(pls - mns * gsl_gaussQD_low_pos[ii], &R2, &radVal, &iR, &fR, &PsiR0, &invPsi, NR, hor, invdR, Nz, Phi, lev, disk, sph, invlogrbin_sph);
    sum += gsl_gaussQD_low_weight[ii] * radVal * _gaussQuad2d4calcRho(zmin, zmax, Nz, ver, invdz, iR, fR, R2, PsiR0, invPsi, Phi, sph, invlogrbin_sph);
  }/* for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--){ */

  /** finalization */
  return (sum * mns);
}


/**
 * @fn swapDblArrays
 *
 * @brief Swap two arrays.
 *
 * @return (p0) data array
 * @return (p1) data array
 */
static inline void swapDblArrays(double **p0, double **p1)
{
  double *tmp;
  tmp = *p0;
  *p0 = *p1;
  *p1 = tmp;
}


/* #define USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE */
#ifdef  USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
static inline double get_disk_rho(double * restrict rho, const double RR, const double Rinv, const double zz, const int maxLev, const double h0, const double h0inv)
{

  const int levR = (int)floor(log2(h0 * ((double)NDISKBIN_HOR - 0.5) * Rinv));
  const int levz = (int)floor(log2(h0 * ((double)NDISKBIN_VER - 0.5) / zz));
  int lev = (levz > levR) ? levR : levz;
  lev = (lev < (maxLev - 1)) ? lev : (maxLev - 1);

  const int ii = (int)fmax(floor(ldexp(RR * h0inv, lev) - 0.5), 0.0);
  const int jj = (int)fmax(floor(ldexp(zz * h0inv, lev) - 0.5), 0.0);

  return (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)]);
}


/* DE formula for zmin <= z <= zmax */
static inline double get_DEformula_part(const double tt, const double zmin, const double zmax, double * restrict rho, const double RR, const double Rinv, const int maxLev, const double h0, const double h0inv)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double inv_cosh_t = 1.0 / cosh(sinh_t);

  const double zz = 0.5 * (zmin * exp(-sinh_t) + zmax * exp(sinh_t)) * inv_cosh_t;

  return (cosh(tt) * inv_cosh_t * inv_cosh_t * get_disk_rho(rho, RR, Rinv, zz, maxLev, h0, h0inv));
}
/* DE formula for 0 <= z <= zmax */
static inline double get_DEformula_full(const double tt, const double zmax, double * restrict rho, const double RR, const double Rinv, const int maxLev, const double h0, const double h0inv)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double inv_cosh_t = 1.0 / cosh(sinh_t);

  const double zz = 0.5 * zmax * exp(sinh_t) * inv_cosh_t;

  return (cosh(tt) * inv_cosh_t * inv_cosh_t * get_disk_rho(rho, RR, Rinv, zz, maxLev, h0, h0inv));
}


static inline double update_trapezoidal_part(const double hh, const double tmin, const double tmax, const double sum, const double zmin, const double zmax, double * restrict rho, const double RR, const double Rinv, const int maxLev, const double h0, const double h0inv)
{
  /** initialization */
  double tmp = 0.0;
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    tmp += get_DEformula_part(tt, zmin, zmax, rho, RR, Rinv, maxLev, h0, h0inv);
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  return (0.5 * sum + hh * tmp);
}
static inline double update_trapezoidal_full(const double hh, const double tmin, const double tmax, const double sum, const double zmax, double * restrict rho, const double RR, const double Rinv, const int maxLev, const double h0, const double h0inv)
{
  /** initialization */
  double tmp = 0.0;
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    tmp += get_DEformula_full(tt, zmax, rho, RR, Rinv, maxLev, h0, h0inv);
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  return (0.5 * sum + hh * tmp);
}


static inline double set_domain_boundary_part(const double hh, double * restrict tmin, double * restrict tmax, const double zmin, const double zmax, double * restrict rho, const double RR, const double Rinv, const int maxLev, const double h0, const double h0inv)
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  double fp = get_DEformula_part(tt, zmin, zmax, rho, RR, Rinv, maxLev, h0, h0inv);
  const double f0 = fp;
  double sum = hh * fp;


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    const double ft = fp;

    tt += hh;    boundary = tt;
    fp = get_DEformula_part(tt, zmin, zmax, rho, RR, Rinv, maxLev, h0, h0inv);

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
    fp = get_DEformula_part(tt, zmin, zmax, rho, RR, Rinv, maxLev, h0, h0inv);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;

  return (sum);
}
static inline double set_domain_boundary_full(const double hh, double * restrict tmin, double * restrict tmax, const double zmax, double * restrict rho, const double RR, const double Rinv, const int maxLev, const double h0, const double h0inv)
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  double fp = get_DEformula_full(tt, zmax, rho, RR, Rinv, maxLev, h0, h0inv);
  const double f0 = fp;
  double sum = hh * fp;


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    const double ft = fp;

    tt += hh;    boundary = tt;
    fp = get_DEformula_full(tt, zmax, rho, RR, Rinv, maxLev, h0, h0inv);

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
    fp = get_DEformula_full(tt, zmax, rho, RR, Rinv, maxLev, h0, h0inv);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;

  return (sum);
}


static inline double integrate_DEformula_part(const double zmin, const double zmax, double * restrict rho, const double RR, const double Rinv, const int maxLev, const double h0, const double h0inv)
{
  const double criteria_abs = 1.0e-12;
  /* const double criteria_rel = 1.0e-10; */
  /* const double criteria_rel = 1.0e-8; */
  /* const double criteria_rel = 1.0e-7; */
  /* const double criteria_rel = 1.0e-6; */
  /* const double criteria_rel = 1.0e-5; */
  const double criteria_rel = 1.0e-4;

  double hh = 1.0;
  double tmin, tmax;
  double sum = set_domain_boundary_part(hh, &tmin, &tmax, zmin, zmax, rho, RR, Rinv, maxLev, h0, h0inv);

  while( true ){
    const double f0 = sum;

    hh *= 0.5;
    sum = update_trapezoidal_part(hh, tmin, tmax, sum, zmin, zmax, rho, RR, Rinv, maxLev, h0, h0inv);

    if( (fabs(sum) > DBL_EPSILON) ? (fabs(1.0 - f0 / sum) <= criteria_rel) : (fabs(sum - f0) <= criteria_abs) )
      break;
  }/* while( true ){ */

  return (sum);
}
static inline double integrate_DEformula_full(const double zmax, double * restrict rho, const double RR, const double Rinv, const int maxLev, const double h0, const double h0inv)
{
  const double criteria_abs = 1.0e-12;
  /* const double criteria_rel = 1.0e-10; */
  /* const double criteria_rel = 1.0e-8; */
  /* const double criteria_rel = 1.0e-7; */
  /* const double criteria_rel = 1.0e-6; */
  /* const double criteria_rel = 1.0e-5; */
  const double criteria_rel = 1.0e-4;

  double hh = 1.0;
  double tmin, tmax;
  double sum = set_domain_boundary_full(hh, &tmin, &tmax, zmax, rho, RR, Rinv, maxLev, h0, h0inv);

  while( true ){
    const double f0 = sum;

    hh *= 0.5;
    sum = update_trapezoidal_full(hh, tmin, tmax, sum, zmax, rho, RR, Rinv, maxLev, h0, h0inv);

    if( (fabs(sum) > DBL_EPSILON) ? (fabs(1.0 - f0 / sum) <= criteria_rel) : (fabs(sum - f0) <= criteria_abs) )
      break;
  }/* while( true ){ */

  return (sum);
}
#endif//USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE


#define SCALE_BY_SURFACE_DENSITY
#define ASSIGN_COARSER_PATCH_FOR_SURFACE_DENSITY

/**
 * @fn getPotDensPair
 *
 * @brief Obtain potential-density pair of the disk component(s).
 *
 * @param (ndisk) number of disk components
 * @param (maxLev) maximum level of nested grid
 * @param (lev) level of nested grid
 * @param (levOld) level of nested grid in the previous iteration
 * @return (disk) physical properties of the disk component
 * @param (Phi_NR) potential outside the boundary
 * @param (Phi_Nz) potential outside the boundary
 * @param (mat) sparse matrix
 * @param (pre) preconditioner
 * @param (stock_inv) tenporary array
 * @param (stock_sub) temporary array
 * @param (stock_tmp) temporary array
 * @param (stock_sum) temporary array
 *
 * @sa fineDensity
 * @sa finePotential
 * @sa coarseDensity
 * @sa coarsePotential
 * @sa getPotentialField
 * @sa gaussQuad2d4calcRho
 */
void getPotDensPair
(const int ndisk,
#   if  defined(ITERATE_VARIABLE_SCALE_HEIGHT) && defined(ENABLE_VARIABLE_SCALE_HEIGHT)
 const int maxLev,
#endif//defined(ITERATE_VARIABLE_SCALE_HEIGHT) && defined(ENABLE_VARIABLE_SCALE_HEIGHT)
#ifdef  USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
 const int highestLev,
#endif//USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
 const int lev, const int levOld, disk_data * restrict disk,
 double * restrict Phi_NR, double * restrict Phi_Nz,
#ifdef  USE_LIS
 LIS_MATRIX lis_mat, LIS_VECTOR lis_b, LIS_VECTOR lis_x, LIS_SOLVER lis_solver,
#else///USE_LIS
 soaBiCGSTAB mat, soaPreConditioning pre,
#endif//USE_LIS
 double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum
 );
void getPotDensPair
(const int ndisk,
#   if  defined(ITERATE_VARIABLE_SCALE_HEIGHT) && defined(ENABLE_VARIABLE_SCALE_HEIGHT)
 const int maxLev,
#endif//defined(ITERATE_VARIABLE_SCALE_HEIGHT) && defined(ENABLE_VARIABLE_SCALE_HEIGHT)
#ifdef  USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
 const int highestLev,
#endif//USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
 const int lev, const int levOld, disk_data * restrict disk, double * restrict Phi_NR, double * restrict Phi_Nz,
#ifdef  USE_LIS
 LIS_MATRIX lis_mat, LIS_VECTOR lis_b, LIS_VECTOR lis_x, LIS_SOLVER lis_solver,
#else///USE_LIS
 soaBiCGSTAB mat, soaPreConditioning pre,
#endif//USE_LIS
 double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum)
{
  __NOTE__("%s\n", "start");


  /** load disk data (common settings for all components) */
  double * RR;  RR  = &(disk[0].hor[INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]);
  double * zz;  zz  = &(disk[0].ver[INDEX2D(maxLev,                NDISKBIN_VER, lev, 0)]);
  double *Phi;  Phi = &(disk[0].pot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, 0)]);

  /** load spherical components(s) data (common settings for all components) */
  profile *sph;  sph = disk[0].prf;
  const double invlogrbin_sph = disk[0].invlogrbin;

  /** set grid points */
  const double    hh = ldexp(disk[0].hh, -lev);
  const double invhh = 1.0 / hh;


  /** set density field and guess potential field */
  if( lev != levOld ){
    if( lev > levOld ){        finePotential(lev, disk[0].pot);      for(int ii = 0; ii < ndisk; ii++)	  fineDensity(lev, disk[ii], *(disk[ii].rho));    }
    if( lev < levOld ){      coarsePotential(lev, disk[0].pot);      for(int ii = 0; ii < ndisk; ii++)	coarseDensity(lev, disk[ii], *(disk[ii].rho));    }
  }/* if( levOld != KICKOFF_POISSON_SOLVER ){ */


  /** iterative process to get the potential-density pair */
  /** modify surface density profile if necessary */
#ifdef  USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
  if( lev > 0 ){
#if 1
    /** adopt Simpson rule for numerical integral */

    if( lev > levOld ){
      /** previous step is coarser grid ==>> evaluate mass above the region */
      for(int kk = 0; kk < ndisk; kk++){
#pragma omp parallel for
	for(int ii = 0; ii < NDISKBIN_HOR; ii++)
	  disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] = 0.5 * disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];

	for(int ll = lev - 1; ll >= 0; ll--){
	  const int ldiff = lev - ll;
	  const double dh = ldexp(disk[0].hh, -ll) / 3.0;

#pragma omp parallel for
	  for(int ii = 0; ii < (NDISKBIN_HOR >> ldiff); ii++){
	    double mass = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, NDISKBIN_VER >> 1)];
	    for(int jj = (NDISKBIN_VER >> 1) + 1; jj < NDISKBIN_VER - 1; jj++)
	      mass += (double)(2 << (jj & 1)) * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)];
	    mass += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, NDISKBIN_VER - 1)];
	    mass *= dh;

	    for(int mm = (ii << ldiff); mm < (ii << ldiff) + (1 << ldiff); mm++)
	      disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, mm)] -= mass;
	  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> ldiff); ii++){ */
	}/* for(int ll = lev - 1; ll >= 0; ll--){ */
      }/* for(int kk = 0; kk < ndisk; kk++){ */
    }/* if( lev > levOld ){ */
/*     else{ */
/*       /\** previous step is finer grid ==>> evaluate mass below the region *\/ */
/*       for(int kk = 0; kk < ndisk; kk++){ */
/* #pragma omp parallel for */
/* 	for(int ii = 0; ii < NDISKBIN_HOR; ii++) */
/* 	  disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] = 0.0; */

/* 	for(int ll = lev; ll < highestLev + 1; ll++){ */
/* 	  const int ldiff = ll - lev; */
/* 	  const double dh = ldexp(disk[0].hh, -ll) / 3.0; */

/* 	  for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
/* 	    int j0 = (ll != highestLev) ? ((ii < (NDISKBIN_HOR >> 1)) ? (NDISKBIN_VER >> 1) : 0) : 0; */
/* 	    double mass = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, j0)]; */
/* 	    for(int jj = j0 + 1; jj < NDISKBIN_VER - 1; jj++) */
/* 	      mass += (double)(2 << (jj & 1)) * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)]; */
/* 	    mass += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, NDISKBIN_VER - 1)]; */
/* 	    mass *= dh; */

/* 	    for(int mm = (ii >> ldiff); mm < (ii >> ldiff) + (1 << ldiff); mm++) */
/* 	      disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, mm)] += mass; */
/* 	  }/\* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ *\/ */
/* 	}/\* for(int ll = lev; ll < highestLev + 1; ll++){ *\/ */
/*       }/\* for(int kk = 0; kk < ndisk; kk++){ *\/ */
/*     }/\* else{ *\/ */






#else
    const double h0 = disk[0].hh;
    const double h0inv = 1.0 / h0;
    const double zmax = h0 * (double)NDISKBIN_VER;
    const double zmid = hh * (double)NDISKBIN_VER;
    const double zmid_inv = 1.0 / zmid_inv;

    for(int kk = 0; kk < ndisk; kk++)
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	const double Ri = RR[ii];
	const double Rinv = 1.0 / Ri;
	const double full = M_PI_4 *  zmid         * integrate_DEformula_full(      zmid, *(disk[kk].rho), Ri, Rinv, highestLev + 1, h0, h0inv);
#if 1
	const double part = 0.5 * disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] - M_PI_4 * (zmax - zmid) * integrate_DEformula_part(zmid, zmax, *(disk[kk].rho), Ri, Rinv, highestLev + 1, h0, h0inv);

	disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] = (full >= part) ? (full) : (part);
	/* if( disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] < 0.0 ){ */
	/*   __KILL__(stderr, "lev = %d, maxLev = %d, R = %e, ans = %e, tot = %e, full = %e, part = %e, zmax = %e, zmid = %e\n", lev, highestLev + 1, Ri, disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], 0.5 * disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], full, part, zmax, zmid); */
	/* } */
#else
	const double part = M_PI_4 * (zmax - zmid) * integrate_DEformula_part(zmid, zmax, *(disk[kk].rho), Ri, Rinv, highestLev + 1, h0, h0inv);

	disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] = (full >= part) ? (full) : (0.5 * disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] - part);
	/* if( disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] < 0.0 ){ */
	/*   __KILL__(stderr, "lev = %d, maxLev = %d, R = %e, ans = %e, tot = %e, full = %e, part = %e, zmax = %e, zmid = %e\n", lev, highestLev + 1, Ri, disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], 0.5 * disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], full, part, zmax, zmid); */
	/* } */
#endif
      }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
#endif
  }/* if( lev > 0 ){ */

#else///USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
#if 1
  /** adopt Simpson rule for numerical integral */
  if( lev > levOld ){
    /** previous step is coarser grid ==>> evaluate mass above the region */
    for(int kk = 0; kk < ndisk; kk++){
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++)
	disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] = 0.5 * disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];

      for(int ll = lev - 1; ll >= 0; ll--){
	const int ldiff = lev - ll;
	const double dh = ldexp(disk[0].hh, -ll) / 3.0;

#pragma omp parallel for
	for(int ii = 0; ii < (NDISKBIN_HOR >> ldiff); ii++){
	  double mass = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, NDISKBIN_VER >> 1)];
	  for(int jj = (NDISKBIN_VER >> 1) + 1; jj < NDISKBIN_VER - 1; jj++)
	    mass += (double)(2 << (jj & 1)) * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)];
	  mass += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, NDISKBIN_VER - 1)];
	  mass *= dh;

	  for(int mm = (ii << ldiff); mm < (ii << ldiff) + (1 << ldiff); mm++)
	    disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, mm)] -= mass;
	}/* for(int ii = 0; ii < (NDISKBIN_HOR >> ldiff); ii++){ */
      }/* for(int ll = lev - 1; ll >= 0; ll--){ */
    }/* for(int kk = 0; kk < ndisk; kk++){ */
  }/* if( lev > levOld ){ */
#else
  if( lev > 0 ){
#ifdef  ASSIGN_COARSER_PATCH_FOR_SURFACE_DENSITY
    const int levCoarse = lev - 1;
    const double hh_tmp = ldexp(disk[0].hh, -levCoarse);
    for(int kk = 0; kk < ndisk; kk++){
#pragma omp parallel for
      for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
	/** evaluate mass to be assigned */
	double mass = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levCoarse, ii, 0)];
	for(int jj = 1; jj < ((NDISKBIN_VER >> 1) - 1); jj++)
	  mass += (double)(2 << (jj & 1)) * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levCoarse, ii, jj)];
	mass += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levCoarse, ii, (NDISKBIN_VER >> 1) - 1)];
	mass *= hh_tmp / 3.0;
	mass *= disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, levCoarse, ii)] * hh_tmp;/**< consider mass within the corresponding shell */

	const int il = ii << 1;
	const int ir = il + 1;
	double ml = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, il, 0)];
	double mr = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ir, 0)];
	for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){
	  const double coeff = (double)(2 << (jj & 1));
	  ml += coeff * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, il, jj)];
	  mr += coeff * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ir, jj)];
	}/* for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){ */
	ml += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, il, NDISKBIN_VER - 1)];
	mr += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ir, NDISKBIN_VER - 1)];
	ml *= disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, il)];/**< consider mass within the corresponding shell */
	mr *= disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ir)];/**< consider mass within the corresponding shell */
	const double inv = 1.0 / (ml + mr);
	/** assign mass to fine grids */
	disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, il)] = ml * inv * mass / (disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, il)] * hh);
	disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ir)] = mr * inv * mass / (disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ir)] * hh);
      }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */
    }/* for(int ii = 0; ii < ndisk; ii++){ */
#else///ASSIGN_COARSER_PATCH_FOR_SURFACE_DENSITY

    /** subtract mass above the domain in this level */
    for(int ll = lev - 1; ll >= 0; ll--){
      const int ldiff = lev - ll;
      const double hh_tmp = ldexp(disk[0].hh, -ll);
#pragma omp parallel for
      for(int ii = 0; ii < (NDISKBIN_HOR >> ldiff); ii++)
	for(int kk = 0; kk < ndisk; kk++){
	  /** evaluate mass to be subtracted */
#if 1
	  double mass = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, NDISKBIN_VER >> 1)];
	  for(int jj = (NDISKBIN_VER >> 1) + 1; jj < NDISKBIN_VER - 1; jj++)
	    mass += (double)(2 << (jj & 1)) * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)];
	  mass += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, NDISKBIN_VER - 1)];
	  mass *= hh_tmp / 3.0;
#else
	  double mass = 0.0;
	  for(int jj = (NDISKBIN_VER >> 1); jj < NDISKBIN_VER; jj++)
	    mass += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)];
	  mass *= hh_tmp;
#endif
	  for(int mm = (ii << ldiff); mm < (ii << ldiff) + (1 << ldiff); mm++){
	    disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, mm)] = 0.5 * disk[ii].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] - mass;
	  }/* for(int mm = 0; mm < (1 << ldiff); mm++){ */
	}/* for(int kk = 0; kk < ndisk; kk++){ */
    }/* for(int ll = lev - 1; ll >= 0; ll--){ */
#endif//ASSIGN_COARSER_PATCH_FOR_SURFACE_DENSITY
  }/* if( lev > 0 ){ */
#endif
#endif//USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE


#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# procedure start: NR = %d, Nz = %d\n", NDISKBIN_HOR, NDISKBIN_VER);
  fflush(stdout);
  int steps = 0;
#endif//PROGRESS_REPORT_ON
  while( true ){
    /** calculate potential field from the given density field */
    /** set density field */
    double *rhoTot;
    if( ndisk > 1 ){
      rhoTot = disk[0].rhoTot;

#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR * NDISKBIN_VER; ii++)
	rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)] = 0.0;
      for(int kk = 0; kk < ndisk; kk++)
#pragma omp parallel for
	for(int ii = 0; ii < NDISKBIN_HOR * NDISKBIN_VER; ii++)
	  rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)] += (*disk[kk].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)];
    }/* if( ndisk > 1 ){ */
    else      rhoTot = *disk[0].rho;

    /** solve Poisson equation using ILU(0) preconditioned BiCGSTAB method */
    getPotentialField(ndisk, disk, lev, hh, invhh, RR, zz, rhoTot, Phi, Phi_NR, Phi_Nz,
#ifdef  USE_LIS
		      lis_mat, lis_b, lis_x, lis_solver,
#else///USE_LIS
		      mat, pre,
#endif//USE_LIS
		      stock_inv, stock_sub, stock_tmp, stock_sum);


#   if  defined(ITERATE_VARIABLE_SCALE_HEIGHT) && defined(ENABLE_VARIABLE_SCALE_HEIGHT)
    getVariableDiskScaleHeight(ndisk, maxLev, lev, disk, sph, invlogrbin_sph);
#endif//defined(ITERATE_VARIABLE_SCALE_HEIGHT) && defined(ENABLE_VARIABLE_SCALE_HEIGHT)


    /** update density field from the derived potential field */
#ifndef SCALE_BY_SURFACE_DENSITY
    static double errMax;
    errMax = 0.0;
#endif//SCALE_BY_SURFACE_DENSITY
    for(int kk = 0; kk < ndisk; kk++){
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	const double Rmin = RR[ii] - 0.5 * hh;
	const double Rmax = RR[ii] + 0.5 * hh;

#ifdef  SCALE_BY_SURFACE_DENSITY
	/** calculate vertical density profile */
	for(int jj = 0; jj < NDISKBIN_VER; jj++)
	  (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] =
	    gaussQuad2d4calcRho(Rmin, Rmax, NDISKBIN_HOR, RR, invhh, zz[jj] - 0.5 * hh, zz[jj] + 0.5 * hh, NDISKBIN_VER, zz, invhh, sph, invlogrbin_sph, Phi, lev, disk[kk]);
#else///SCALE_BY_SURFACE_DENSITY
	const double rho0 = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] /
	  (DBL_MIN + gaussQuad2d4calcRho(Rmin, Rmax, NDISKBIN_HOR, RR, invhh, 0.0, hh, NDISKBIN_VER, zz, invhh, sph, invlogrbin_sph, Phi, lev, disk[kk]));

	/** calculate vertical density profile */
	(*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)];
	for(int jj = 1; jj < NDISKBIN_VER; jj++)
	  (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] =
	    rho0 * gaussQuad2d4calcRho(Rmin, Rmax, NDISKBIN_HOR, RR, invhh, zz[jj] - 0.5 * hh, zz[jj] + 0.5 * hh, NDISKBIN_VER, zz, invhh, sph, invlogrbin_sph, Phi, lev, disk[kk]);
#endif//SCALE_BY_SURFACE_DENSITY

	/** calculate surface density */
#if 1
	double Sigma = (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)];
	for(int jj = 1; jj < NDISKBIN_VER - 1; jj++)
	  Sigma += (double)(2 << (jj & 1)) * (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)];
	Sigma += (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, NDISKBIN_VER - 1)];
	Sigma *= hh / 3.0;
#else
	double Sigma = 0.0;
	for(int jj = 0; jj < NDISKBIN_VER; jj++)
	  Sigma += (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)];
	Sigma *= hh;
#endif

	/** calibrate surface density */
	const double Mscale = disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] / (DBL_MIN + Sigma);
	for(int jj = 0; jj < NDISKBIN_VER; jj++)
	  (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] *= Mscale;

#if 0
	if( fpclassify(Mscale) != FP_NORMAL ){
	  __KILL__(stderr, "ERROR: lev = %d, ii = %d, Mscale = %e, Sigma = %e, assigned = %e, hh = %e, zd = %e\n", lev, ii, Mscale, Sigma, disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], hh, disk[kk].zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);
	}
#endif

#ifndef SCALE_BY_SURFACE_DENSITY
	const double errVal = (disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] > NEGLECT_DENSITY_MINIMUM * disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)]) ? fabs(Mscale - 1.0) : (0.0);
#pragma omp flush(errMax)
	if( errVal > errMax )
#pragma omp critical
	  if( errVal > errMax )
	    errMax = errVal;
#endif//SCALE_BY_SURFACE_DENSITY
      }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
    }/* for(int kk = 0; kk < ndisk; kk++){ */


    /** convergence tests */
    /** convergence test for the surface density */
    bool converge = true;
#ifndef SCALE_BY_SURFACE_DENSITY
    if( errMax > CONVERGENCE_POTDENSPAIR )      converge = false;
#ifdef  PROGRESS_REPORT_ON
    fprintf(stdout, "# %d-th iteration: surface density error is %e, procedure is %s\n", steps, errMax, (converge ? ("    converged") : ("not converged")));
    fflush(stdout);
#endif//PROGRESS_REPORT_ON
#endif//SCALE_BY_SURFACE_DENSITY

    /** convergence test for the volume-density field */
    if( converge ){
      static double worst;
      worst = DBL_EPSILON;
      for(int kk = 0; kk < ndisk; kk++){
#pragma omp parallel for
	for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	  double errVal = 0.0;
	  if( (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] > NEGLECT_DENSITY_MINIMUM * (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 0, 0)] )
	    for(int jj = 0; jj < NDISKBIN_VER; jj++){
	      const double err =
		((*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] > NEGLECT_DENSITY_MINIMUM * (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)]) ?
		(fabs(1.0 - (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] / ((*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)]))) : (0.0);
	      errVal = (errVal > err) ? errVal : err;
	    }/* for(int jj = 0; jj < Nz; jj++){ */
#pragma omp flush(worst)
	  if( errVal > worst )
#pragma omp critical
	    if( errVal > worst )
	      worst = errVal;
	}/* for(int ii = 0; ii < NR; ii++){ */
      }/* for(int kk = 0; kk < ndisk; kk++){ */
      if( worst > CONVERGENCE_POTDENSPAIR )
	converge = false;

#ifdef  PROGRESS_REPORT_ON
      fprintf(stdout, "# %d-th iteration:        density error is %e, procedure is %s\n", steps, worst, (converge ? ("    converged") : ("not converged")));
      fflush(stdout);
#endif//PROGRESS_REPORT_ON
    }/* if( converge ){ */

    for(int kk = 0; kk < ndisk; kk++)
      swapDblArrays(disk[kk].rho, disk[kk].rhoSum);

    /** final confirmation */
    if( converge )
      break;

#ifdef  PROGRESS_REPORT_ON
    steps++;
#endif//PROGRESS_REPORT_ON
  }/* while( true ){ */


#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# procedure finish after %d iteration(s): NR = %d, Nz = %d\n#\n#\n", 1 + steps, NDISKBIN_HOR, NDISKBIN_VER);
  fflush(stdout);
#endif//PROGRESS_REPORT_ON


  __NOTE__("%s\n", "end");
}


/**
 * @fn makeDiskPotentialTable
 *
 * @brief Obtain the density and potential field of the disk component(s).
 *
 * @param (ndisk) number of disk components
 * @param (maxLev) maximum level of nested grid
 * @return (disk) physical properties of the disk component
 *
 * @sa setSparseMatrix
 * @sa getILU0
 * @sa setColumnDensityProfile
 * @sa gaussQD4GDpot
 * @sa getPotDensPair
 */
void makeDiskPotentialTable(const int ndisk, const int maxLev, disk_data * restrict disk)
{
  __NOTE__("%s\n", "start");


  /** initialize the table for Gaussian Quadrature provided by GSL */
  for(int ii = 0; ii < NTBL_GAUSS_QD_LOW; ii++){
    gsl_gaussQD_low_pos   [ii] = 0.0;
    gsl_gaussQD_low_weight[ii] = 0.0;
  }/* for(int ii = 0; ii < NTBL_GAUSS_QD_LOW; ii++){ */
  gsl_integration_glfixed_table *tab;
  tab = gsl_integration_glfixed_table_alloc(NINTBIN_LOW);
  int max = (NINTBIN_LOW >> 1) + (NINTBIN_LOW & 1);
  for(int ii = 0; ii < max; ii++){
    gsl_gaussQD_low_pos   [ii] = (*tab).x[(max - 1) - ii];
    gsl_gaussQD_low_weight[ii] = (*tab).w[(max - 1) - ii];
  }/* for(int ii = 0; ii < max; ii++){ */
  gsl_integration_glfixed_table_free(tab);

  /** outer boundary condition of the potential field */
  static double Phi_NR[NDISKBIN_VER], Phi_Nz[NDISKBIN_HOR];

  /** allocate memory for the iterative solver of the Poisson equation */
  /** sparse matrix in CRS format */
#ifdef  USE_LIS
  LIS_MATRIX lis_mat;
  lis_matrix_create(0, &lis_mat);
  lis_matrix_set_size(lis_mat, 0, NROW_CG);
  setSparseMatrix(lis_mat);

  LIS_VECTOR lis_b, lis_x;
  lis_vector_create(0, &lis_b);
  lis_vector_set_size(lis_b, 0, NCOL_CG);
  lis_vector_duplicate(lis_b, &lis_x);

  LIS_SOLVER lis_solver;
  lis_solver_create(&lis_solver);
  /* lis_solver_set_option("-i bicgstab", lis_solver);/\**< 2.564329e+02 sec for diskTbl *\/ */
  /* lis_solver_set_option("-i bicgstab -p ilu", lis_solver);/\**< 1.304378e+02 sec for diskTbl *\/ */
  lis_solver_set_option("-i bicgstab -p ilut", lis_solver);/**< 1.136713e+02 sec for diskTbl */
  /* lis_solver_set_option("-i minres", lis_solver);/\**< 2.602234e+02 sec for diskTbl *\/ */
  /* lis_solver_set_option("-i minres -p ilu", lis_solver); /\**< too long elapsed time for diskTbl *\/ */
  /* lis_solver_set_option("-i minres -p ilut", lis_solver); /\**< 2.602234e+02 sec for diskTbl *\/ */
  lis_solver_set_option("-tol 1.0e-10 -maxiter 16384", lis_solver);
#else///USE_LIS
  crs mat, ilu;
  static double mat_val[NNZ_CG], ilu_val[NNZ_CG];
  static    int mat_col[NNZ_CG], ilu_col[NNZ_CG], mat_row[NROW_CG + 1], ilu_row[NROW_CG + 1];
  mat.val = mat_val;  mat.col = mat_col;  mat.row = mat_row;
  ilu.val = ilu_val;  ilu.col = ilu_col;  ilu.row = ilu_row;
  /** vector used in BiCGSTAB method */
  static double vec[NCOL_CG], res[NCOL_CG], sdw[NCOL_CG], mid[NCOL_CG], tmp[NCOL_CG];
  static double Api[NCOL_CG], Ati[NCOL_CG], Kri[NCOL_CG], Kpi[NCOL_CG], Kti[NCOL_CG];

  /** execute first touch for memories related to BiCGSTAB */
#pragma omp parallel
  {
#pragma omp for nowait
    for(int ii = 0; ii < NCOL_CG; ii++){
      vec[ii] = res[ii] = sdw[ii] = mid[ii] = tmp[ii] = 0.0;
      Api[ii] = Ati[ii] = Kri[ii] = Kpi[ii] = Kti[ii] = 0.0;
    }/* for(int ii = 0; ii < NCOL_CG; ii++){ */
#pragma omp for nowait
    for(int ii = 0; ii < NNZ_CG; ii++){
      mat_val[ii] = 0.0;
      mat_col[ii] = 0;
    }/* for(int ii = 0; ii < NNZ_CG; ii++){ */
#pragma omp for nowait
    for(int ii = 0; ii < NROW_CG; ii++)
      mat_row[ii] = 0;
  }

  soaBiCGSTAB smat;
  smat.mat = mat;  smat.vec = vec;  smat.res = res;  smat.sdw = sdw;
  smat.mid = mid;  smat.tmp = tmp;  smat.Api = Api;  smat.Ati = Ati;
  soaPreConditioning pre;
  pre.ilu = ilu;  pre.Kri = Kri;  pre.Kpi = Kpi;  pre.Kti = Kti;

  /** prepare Poisson equation in matrix form (matrix part) */
  setSparseMatrix(smat.mat);
  getILU0(NDISKBIN_HOR * NDISKBIN_VER, smat.mat, pre.ilu);
#endif//USE_LIS


  /** set column density profile */
  setColumnDensityProfile(ndisk, maxLev, disk, disk[0].prf, disk[0].invlogrbin);

  /** temporal value for generating potential-density pair of the disk component(s) */
  for(int kk = 0; kk < ndisk; kk++)
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
      disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)] = 0.5 * disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)];


  /** calculate GD formalism */
  double *stock_inv;  stock_inv = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  double *stock_sub;  stock_sub = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  double *stock_tmp;  stock_tmp = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  double *stock_sum;  stock_sum = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  if( stock_inv == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_inv\n");  }
  if( stock_sub == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_sub\n");  }
  if( stock_tmp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_tmp\n");  }
  if( stock_sum == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_sum\n");  }

  const double Rmax = disk[0].Rmax;
  const double    abin = Rmax / (double)NRADBIN;
  const double invabin = 1.0 / abin;
  for(int kk = 0; kk < ndisk; kk++){
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      disk[kk].prf[ii].rho = (double)ii * abin;/**< store the position ``a'' */
#ifdef  NDIVIDE_GAUSSQD4DISK
      double sum = 0.0;
      double Rin = disk[kk].prf[ii].rho;
      const double a2 = Rin * Rin;
      const double Rbin = (Rmax - Rin) / (double)NDIVIDE_GAUSSQD4DISK;
      for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){
	sum += gaussQD4GDpot(Rin, Rin + Rbin, a2, disk[kk]);
	Rin += Rbin;
      }/* for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){ */
      disk[kk].prf[ii].enc = disk[kk].cfg->Sigma0 * sum;/**< store the integral */
#else///NDIVIDE_GAUSSQD4DISK
      disk[kk].prf[ii].enc = disk[kk].cfg->Sigma0 * gaussQD4GDpot(disk[kk].prf[ii].rho, Rmax, disk[kk].prf[ii].rho * disk[kk].prf[ii].rho, disk[kk]);/**< store the integral */
#endif//NDIVIDE_GAUSSQD4DISK
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */

    disk[kk].prf[0].psi = (disk[kk].prf[1].enc - disk[kk].prf[0].enc) * invabin;
#pragma omp parallel for
    for(int ii = 1; ii < NRADBIN - 1; ii++)
      disk[kk].prf[ii].psi = (disk[kk].prf[ii + 1].enc - disk[kk].prf[ii - 1].enc) * 0.5 * invabin;/* store the derivative */
    disk[kk].prf[NRADBIN - 1].psi = (disk[kk].prf[NRADBIN - 1].enc - disk[kk].prf[NRADBIN - 2].enc) * invabin;
  }/* for(int kk = 0; kk < ndisk; kk++){ */


  /** iterative procedure adopting Full Multigrid Algorithm (gamma = 2) */
  int lev = 0;
  int old = lev;
  int inc = 1;
  int top = 1;
  int num = 0;
#ifdef  USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
  int highestLev = lev;
#endif//USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
  while( true ){
#ifdef  PROGRESS_REPORT_ON
    fprintf(stdout, "# grid level: lev = %d, oldLev = %d\n", lev, old);
    fflush(stdout);
#endif//PROGRESS_REPORT_ON

    getPotDensPair(ndisk,
#   if  defined(ITERATE_VARIABLE_SCALE_HEIGHT) && defined(ENABLE_VARIABLE_SCALE_HEIGHT)
		   maxLev,
#endif//defined(ITERATE_VARIABLE_SCALE_HEIGHT) && defined(ENABLE_VARIABLE_SCALE_HEIGHT)
#ifdef  USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
		   highestLev,
#endif//USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
		   lev, old, disk, Phi_NR, Phi_Nz,
#ifdef  USE_LIS
		   lis_mat, lis_b, lis_x, lis_solver,
#else///USE_LIS
		   smat, pre,
#endif//USE_LIS
		   stock_inv, stock_sub, stock_tmp, stock_sum);

#ifdef  USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE
    highestLev = (lev > highestLev) ? lev : highestLev;
#endif//USE_DEFORMULA_FOR_SURFACE_DENSITY_ESTIMATE

    /** procedure proposed by Press & Teukolsky (1991), Computers in Physics 5, 514, Multigrid Methods for Boundary Value Problems. I */
    if( lev == top ){   num++;  inc = -1;
      if( num == 2 ){	top++;	num =  0;      }
    }/* if( lev == top ){ */
    /** in case of the coarsest grid */
    if( (lev == 0) && (inc == -1) )      inc = 1;
    /** exit the loop if we found the solution in the finest grid */
    if( lev == (maxLev - 1) )      break;

    /** go to the next level */
    old  = lev;
    lev += inc;
  }/* while( true ){ */

#if 0
  for(int ll = maxLev - 2; ll >= 0; ll--)
    coarsePotential(ll, disk[0].pot );
#endif

  for(int kk = 0; kk < ndisk; kk++)
    for(int ll = maxLev - 2; ll >= 0; ll--)
      coarseDensity  (ll, disk[kk], *(disk[kk].rho));

  free(stock_inv);
  free(stock_sub);
  free(stock_tmp);
  free(stock_sum);


  /** calculate 2 \pi \int_0^R dR' R' Sigma(R') */
#ifdef  NDIVIDE_GAUSSQD4DISK
  for(int kk = 0; kk < ndisk; kk++){
    double Min = 0.0;
    double Rin = 0.0;

    for(int lev = maxLev - 1; lev >= 0; lev--){
      const int nsub = (lev != (maxLev - 1)) ? ((NDISKBIN_HOR >> 1) / NDIVIDE_GAUSSQD4DISK) : (NDISKBIN_HOR / NDIVIDE_GAUSSQD4DISK);
      int head = (lev != (maxLev - 1)) ?  (NDISKBIN_HOR >> 1)                         : (0);
      disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head)] = Min;

      for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){
	const int tail = head + nsub;
#pragma omp parallel for
	for(int ii = head; ii < tail; ii++)
	  disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii + 1)] = Min + gaussQD4encSigma(Rin, disk[kk].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii + 1)], disk[kk]);

	head = tail;
	Rin = disk[kk].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head)];
	Min = disk[kk].	    enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head)];
      }/* for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){ */

#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR + 1; ii++)
	disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii)] *= 2.0 * M_PI * disk[kk].cfg->Sigma0;
    }/* for(int lev = maxLev - 1; lev >= 0; lev--){ */
  }/* for(int kk = 0; kk < ndisk; kk++){ */
#else///NDIVIDE_GAUSSQD4DISK
  for(int kk = 0; kk < ndisk; kk++){
    double Min = 0.0;
    double Rin = 0.0;

    for(int lev = maxLev - 1; lev >= 0; lev--){
      const int nsub = (lev != (maxLev - 1)) ? (NDISKBIN_HOR >> 1) : (NDISKBIN_HOR);
      const int head = (lev != (maxLev - 1)) ? (NDISKBIN_HOR >> 1) : (0);
      disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head)] = Min;

#pragma omp parallel for
      for(int ii = head; ii < head + nsub; ii++)
	disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii + 1)] = Min + gaussQD4encSigma(Rin, disk[kk].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii + 1)], disk[kk]);

      Rin = disk[kk].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head + nsub)];
      Min = disk[kk].	  enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head + nsub)];

#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR + 1; ii++)
	disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii)] *= 2.0 * M_PI * disk[kk].cfg->Sigma0;
    }/* for(int lev = maxLev - 1; lev >= 0; lev--){ */
  }/* for(int kk = 0; kk < ndisk; kk++){ */
#endif//NDIVIDE_GAUSSQD4DISK
  for(int kk = 0; kk < ndisk; kk++){
    for(int lev = maxLev - 2; lev >= 0; lev--){
      for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++)
	disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii)] = disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev + 1, ii << 1)];
    }/* for(int lev = maxLev - 2; lev >= 0; lev--){ */
  }/* for(int kk = 0; kk < ndisk; kk++){ */


  /** calculate int_0^z dz' \rho(R, z') */
  for(int kk = 0; kk < ndisk; kk++){
    for(int ll = 0; ll < maxLev; ll++){
      const double hh = ldexp(disk[kk].hh, -ll);

#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	const double dV = hh * hh * disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)];
	double sum = 0.0;
	(*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii, 0)] = sum;

	for(int jj = 0; jj < NDISKBIN_VER; jj++){
	  sum += dV * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)];
	  (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii, jj + 1)] = sum;
	}/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
      }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
    }/* for(int ll = 0; ll < maxLev; ll++){ */

    /** normalization */
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)] = 1.0 / (DBL_MIN + (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, 0, ii, NDISKBIN_VER)]);
      for(int ll = 1; ll < maxLev; ll++){
	disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)] = disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, ll - 1, ii >> 1)] * ((*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii, NDISKBIN_VER)] + (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii ^ 1, NDISKBIN_VER)]) / (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii, NDISKBIN_VER)];
      }/* for(int ll = 1; ll < maxLev; ll++){ */
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */

    for(int ll = 0; ll < maxLev; ll++)
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++)
	for(int jj = 0; jj < NDISKBIN_VER + 1; jj++)
	  (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii, jj)] *= disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)];
  }/* for(int kk = 0; kk < ndisk; kk++){ */

#ifdef  USE_LIS
  lis_matrix_destroy(lis_mat);
  lis_vector_destroy(lis_b);
  lis_vector_destroy(lis_x);
  lis_solver_destroy(lis_solver);
#endif//USE_LIS


  __NOTE__("%s\n", "end");
}
