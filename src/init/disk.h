/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:24:08
 *                                                                       *
 *    Header File for Definition to generate initial condition of disk   *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef DISK_H
#define DISK_H
//-------------------------------------------------------------------------
#include "macro.h"
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#include "../init/profile.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define NDIVIDE_GAUSSQD4DISK (4)
//-------------------------------------------------------------------------
/* this option is introduced to ensure z_d < r_s in the innermost region; this is problem on thick disk models */
#define ENABLE_VARIABLE_SCALE_HEIGHT
#define CHECK_OSTRIKER_PEEBLES_CRITERION
#define USE_ORIGINAL_VDISP_ESTIMATOR
#define SPEEDUP_CONVERGENCE
/* #define USE_POTENTIAL_SCALING_SCHEME */
//-------------------------------------------------------------------------
#define USE_GD_FORM_POTENTIAL
#ifdef  USE_GD_FORM_POTENTIAL
#define NGDPOT (1048576)
#endif//USE_GD_FORM_POTENTIAL
//-------------------------------------------------------------------------
#define CONVERGENCE_BICGSTAB    (1.0e-10)
#define CONVERGENCE_POTDENSPAIR (1.0e-4)
#define NEGLECT_DENSITY_MINIMUM (1.0e-10)
//-------------------------------------------------------------------------
#ifndef RUN_ON_PC
#define INCREASE_TABLE_STEPS (4)
#define NDISKBIN_VER (512)
#else///RUN_ON_PC
#define INCREASE_TABLE_STEPS (3)
#define NDISKBIN_VER (256)
#endif//RUN_ON_PC
//-------------------------------------------------------------------------
#define NHOR_OVER_NVER (8)
#define NDISKBIN_HOR (NHOR_OVER_NVER * NDISKBIN_VER)
//-------------------------------------------------------------------------
#define NCOL_CG (NDISKBIN_HOR * NDISKBIN_VER)
#define NROW_CG (NDISKBIN_HOR * NDISKBIN_VER)
#define  NNZ_CG (5 * NDISKBIN_HOR * NDISKBIN_VER - 2 * (NDISKBIN_HOR + NDISKBIN_VER))
//-------------------------------------------------------------------------
#define NDISKBIN_RAD (16384)
#define SUBDIVIDE_NUM_HOR (NDISKBIN_RAD / NDISKBIN_HOR)
#define SUBDIVIDE_NUM_VER (NDISKBIN_RAD / NDISKBIN_VER)
//-------------------------------------------------------------------------
#define BILINEAR_INTERPOLATION
/* #define BICUBIC_INTERPOLATION */
//-------------------------------------------------------------------------
#   if  defined(BILINEAR_INTERPOLATION) && defined(BICUBIC_INTERPOLATION)
#undef          BILINEAR_INTERPOLATION
#endif//defined(BILINEAR_INTERPOLATION) && defined(BICUBIC_INTERPOLATION)
//-------------------------------------------------------------------------
static const double diskMaxLength = 10.0;
static const double diskMaxHeight = 20.0;
/* static const double diskMaxLength =  8.0; */
/* static const double diskMaxHeight = 16.0; */
//-------------------------------------------------------------------------
typedef struct
{
  double sersic_ninv, sersic_b;
  double Rcutoff, invRsmooth;
} disk_util;
//-------------------------------------------------------------------------
typedef struct
{
  /* common arrays */
  double *pot;
  double *hor, Rbin, invRbin;
  double *ver, zbin, invzbin;
  double *dPhidR, *d2PhidR2;
  double *radSph, *rhoSph, *encSph;/* NDISKBIN_RAD elements */
  profile *prf;
  double logrbin, invlogrbin;
  /* individual arrays */
  double **rho, **rhoSum, *arr0, *arr1, *rhoTot;
  double *Sigma, *sigmaz, *enc;/* enc は，2\pi \int_0^R dR' R' Sigma(R') の結果を格納 */
  /* physical properties */
  profile_cfg *cfg;
  double invRd;/* set @ makeDiskPotentialTable */
  /* double sersic_ninv, sersic_b; */
  /* double Rcutoff, invRsmooth; */
  /* double (*getColumnDensity)(double, double, double, double, double, double); */
  disk_util util;
  double (*getColumnDensity)(double, double, disk_util);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  double *zd;
  /* double zd_scale_min, zd_scale_max, zd_scale_R0, zd_scale_inv; */
  /* double rs_spheroid; */
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
  double Krand_sph;
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION
} disk_data;
//-------------------------------------------------------------------------
#ifdef  USE_ORIGINAL_VDISP_ESTIMATOR
//-------------------------------------------------------------------------
/* static double VEL0, VEL1; */
/* #define DISK_PERP_VDISP(sigmaz, vcirc, frac) (VEL0=(sigmaz), VEL1=((frac) * (vcirc)), (((VEL0) < (VEL1)) ? (VEL0) : (VEL1))) */
#define DISK_PERP_VDISP(sigmaz, vcirc, frac) (fmin(sigmaz, (frac) * (vcirc)))
//-------------------------------------------------------------------------
#else///USE_ORIGINAL_VDISP_ESTIMATOR
//-------------------------------------------------------------------------
/* same method with GalactICS */
#define DISK_RADIAL_VDISP2(sz0_2, RR, invRd) ((sz0_2) *      exp(-(RR) * (invRd)))
#define DISK_RADIAL_VDISP( sz0  , RR, invRd) ((sz0  ) * sqrt(exp(-(RR) * (invRd))))
//-------------------------------------------------------------------------
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "disk.c"
//-------------------------------------------------------------------------
/* #ifdef  ENABLE_VARIABLE_SCALE_HEIGHT */
/* void selectDominantComponent(const int ndisk, disk_data * restrict disk, const int skind, profile * restrict * prf, const double invlogbin, profile_cfg * restrict cfg); */
/* #endif//ENABLE_VARIABLE_SCALE_HEIGHT */
//-------------------------------------------------------------------------
void allocDiskProfile
(const int ndisk, disk_data **disk,
 double **hor, double **ver, double **pot, double **dPhidR, double **d2PhidR2, double **radSph, double **rhoSph, double **encSph,
 double **rho, double **rhoSum, double **rhoTot, double **Sigma, double **sigmaz, double **enc
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 , double **zd
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 );
void  freeDiskProfile
(const int ndisk, disk_data  *disk,
 double  *hor, double  *ver, double  *pot, double  *dPhidR, double  *d2PhidR2, double  *radSph, double  *rhoSph, double  *encSph,
 double  *rho, double  *rhoSum, double  *rhoTot, double  *Sigma, double  *sigmaz, double  *enc
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 , double  *zd
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 );
//-------------------------------------------------------------------------
void makeDiskPotentialTable(const int ndisk, disk_data * restrict disk);
//-------------------------------------------------------------------------
void integrateSphericalDensityProfile(const int ndisk, disk_data *disk);
//-------------------------------------------------------------------------
void diffAxisymmetricPotential(const disk_data disk);
void calcVerticalVdisp(const int ndisk, disk_data *disk_info);
//-------------------------------------------------------------------------
void distributeDiskParticles(ulong *Nuse, iparticle body, const real mass, const disk_data disk);
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//DISK_H
//-------------------------------------------------------------------------
