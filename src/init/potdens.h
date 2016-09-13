/*************************************************************************\
 *                                                                       *
                  last updated on 2016/09/13(Tue) 11:47:50
 *                                                                       *
 *    Header File for Definition to generate initial condition of disk   *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef POTDENS_H
#define POTDENS_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MACRO_H
#include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef PROFILE_H
#include "../init/profile.h"
#endif//PROFILE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define CHECK_OSTRIKER_PEEBLES_CRITERION
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define NDISKBIN_VER (32) */
#define NDISKBIN_VER (64)
/* NHOR_OVER_NVER is defined as an exponent of two */
#define NHOR_OVER_NVER (2)
#define NDISKBIN_HOR (NDISKBIN_VER << NHOR_OVER_NVER)
//-------------------------------------------------------------------------
#define NDISKBIN_RAD (16384)
//-------------------------------------------------------------------------
#define  NNZ_CG (5 * NDISKBIN_HOR * NDISKBIN_VER - 2 * (NDISKBIN_HOR + NDISKBIN_VER))
#define NROW_CG (NDISKBIN_HOR * NDISKBIN_VER)
#define NCOL_CG (NDISKBIN_HOR * NDISKBIN_VER)
//-------------------------------------------------------------------------
#define DISK_MAX_SAFETY (2.0)
#define DISK_MAX_LENGTH (10.0)
#define DISK_MIN_LENGTH (3.90625e-3)
//-------------------------------------------------------------------------
#define NDIVIDE_GAUSSQD4DISK (4)
//-------------------------------------------------------------------------
/* #define PROHIBIT_EXTRAPOLATION *//* <-- currently, not implemented */
#define ENABLE_VARIABLE_SCALE_HEIGHT
//-------------------------------------------------------------------------
#define USE_GD_FORM_POTENTIAL
#ifdef  USE_GD_FORM_POTENTIAL
#define NGDPOT (1048576)
#endif//USE_GD_FORM_POTENTIAL
//-------------------------------------------------------------------------
#define USE_ELLIPTIC_INTEGRAL
//-------------------------------------------------------------------------
#define CONVERGENCE_BICGSTAB    (1.0e-10)
#define CONVERGENCE_POTDENSPAIR (1.0e-4)
#define NEGLECT_DENSITY_MINIMUM (1.0e-10)
//-------------------------------------------------------------------------
#define INDEX4D(n0, n1, n2, n3, ii, jj, kk, ll) ((ll) + (n3) * ((kk) + (n2) * ((jj) + (n1) * (ii))))
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
  double sersic_ninv, sersic_b;
  double Rcutoff, invRsmooth;
} disk_util;
//-------------------------------------------------------------------------
typedef struct
{
  /* physical properties of disk component */
  profile_cfg *cfg;
  /* physical properties of spherical component(s) */
  profile *prf;
  /* arrays */
  double *hor;/* [nest level][NR] array */
  double *ver;/* [nest level][Nz] array */
  double *pot, *rhoTot;/* [nest level][NR][Nz] arrays */
  double **rho, **rhoSum, *rho0, *rho1;/* [Ndisk][nest level][NR][Nz] arrays */
  double *dPhidR, *d2PhidR2;/* [nest level][NR][Nz] arrays */
  double *Sigma, *sigmaz, *enc;/* [Ndisk][nest level][NR] array */
  disk_util util;
  double (*getColumnDensity)(double, double, disk_util);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  double *zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  /* arrays for spherical averaged profile */
  double *radSph, *rhoSph, *encSph;/* [NDISKBIN_RAD] array */
  /* arrays for spline fit */
  double *spline_xx, *spline_ff, *spline_f2, *spline_bp;
  double Rmax, zmax, hh;/* configuration of domain */
  double invRd;/* configuration of disk component */
  double logrbin, invlogrbin;/* configuration of spherical component(s) */
#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
  double Krand_sph;
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION
} disk_data;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "potdens.c"
//-------------------------------------------------------------------------
void   freeDiskProfile
(const int ndisk, disk_data  *disk,
 double  *hor, double  *ver,
 double  *pot, double  *rho0, double  *rho1, double  *rhoTot, double  *dPhidR, double  *d2PhidR2,
 double  *Sigma, double  *vsigz, double  *enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 double  *zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 double  *radSph, double  *rhoSph, double  *encSph,
 double  *spline_xx, double  *spline_ff, double  *spline_f2, double  *spline_bp);
//-------------------------------------------------------------------------
void allocDiskProfile
(const int ndisk, disk_data **disk, profile_cfg *disk_cfg, int *maxLev,
 double **hor, double **ver,
 double **pot, double **rho0, double **rho1, double **rhoTot,
 double **dPhidR, double **d2PhidR2,
 double **Sigma, double **vsigz, double **enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 double **zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 double **radSph, double **rhoSph, double **encSph,
 double **spline_xx, double **spline_ff, double **spline_f2, double **spline_bp);
//-------------------------------------------------------------------------
void makeDiskPotentialTable(const int ndisk, const int maxLev, disk_data * restrict disk);
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//POTDENS_H
//-------------------------------------------------------------------------
