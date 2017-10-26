/**
 * @file potdens.h
 *
 * @brief Header file for calculating potential-density pair of disk component(s)
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/10/26 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef POTDENS_H
#define POTDENS_H


#include "macro.h"

#include "../init/profile.h"


/**
 * @def NDISKBIN_VER
 *
 * @brief number of grid points of the density and potential fields in z-direction
 */
#define NDISKBIN_VER (64)

/**
 * @def NHOR_OVER_NVER
 *
 * @brief NR / Nz, defined as an exponent of two
 */
#define NHOR_OVER_NVER (2)

/**
 * @def NDISKBIN_HOR
 *
 * @brief number of grid points of the density and potential fields in R-direction
 */
#define NDISKBIN_HOR (NDISKBIN_VER << NHOR_OVER_NVER)


/** number of arrays for BiCGSTAB method */
#define  NNZ_CG (5 * NDISKBIN_HOR * NDISKBIN_VER - 2 * (NDISKBIN_HOR + NDISKBIN_VER))
#define NROW_CG (NDISKBIN_HOR * NDISKBIN_VER)
#define NCOL_CG (NDISKBIN_HOR * NDISKBIN_VER)


/** macros for setting physical scale of the density and potential fields of the disk component(s) */
#define DISK_MAX_SAFETY (2.0)
#define DISK_MAX_LENGTH (10.0)
#define DISK_MIN_LENGTH (3.90625e-3)


/* #define NDIVIDE_GAUSSQD4DISK (4) */
#define NDIVIDE_GAUSSQD4DISK (8)
/* #define NDIVIDE_GAUSSQD4DISK (16) */


/**
 * @def ENABLE_VARIABLE_SCALE_HEIGHT
 *
 * @brief enable to remove the needle-like structure
 * @detail Equation (8) in Miki & Umemura (in preparation)
 */
#define ENABLE_VARIABLE_SCALE_HEIGHT

/**
 * @def ITERATE_VARIABLE_SCALE_HEIGHT
 *
 * @brief determine the scale height of disk component(s) by iterations
 */
#define ITERATE_VARIABLE_SCALE_HEIGHT

/**
 * @def DISK_DIMMING_HEIGHT
 *
 * @brief parameter to remove the needle-like structure
 * @detail see Equation (8) in Miki & Umemura (in preparation)
 */
#define DISK_DIMMING_HEIGHT (16.0)

/**
 * @def DISK_DIMMING_HEIGHT_INV
 *
 * @brief inverse of diskDimmingHeight
 */
#define DISK_DIMMING_HEIGHT_INV (6.25e-2)

/**
 * @def ADDITIONAL_CONDITION_FOR_SCALE_HEIGHT
 *
 * @brief enable to remove the needle-like structure
 */
#define ADDITIONAL_CONDITION_FOR_SCALE_HEIGHT

/**
 * @def DISK_DIMMING_SCALE
 *
 * @brief parameter to remove the needle-like structure for very thick disk
 */
#define DISK_DIMMING_SCALE (5.0)


/**
 * @def CONVERGENCE_BICGSTAB
 *
 * @brief tolerance value for BiCGSTAB method
 */
#define CONVERGENCE_BICGSTAB    (1.0e-10)

/**
 * @def CONVERGENCE_POTDENSPAIR
 *
 * @brief tolerance value for iterating potential-density pair of the disk component(s)
 */
#define CONVERGENCE_POTDENSPAIR (1.0e-4)
/* #define CONVERGENCE_POTDENSPAIR (3.1e-4) */

/**
 * @def NEGLECT_DENSITY_MINIMUM
 *
 * @brief tolerance value for low density regions
 */
#define NEGLECT_DENSITY_MINIMUM (1.0e-10)


#define INDEX4D(n0, n1, n2, n3, ii, jj, kk, ll) ((ll) + (n3) * ((kk) + (n2) * ((jj) + (n1) * (ii))))

/* #define NDISKBIN_RAD (16384) */
#define NDISKBIN_RAD (131072)
/* #define NDISKBIN_RAD (1048576) */


/**
 * @struct disk_util
 *
 * @brief structure for disk component(s)
 */
typedef struct
{
  double *xx, *ff, *f2, *bp;/**< arrays for spline fit (column density profile in table form) */
  double sersic_ninv, sersic_b;
  double Rcutoff, invRsmooth;
  int num;/**< # of arrays for spline fit */
} disk_util;

/**
 * @struct disk_data
 *
 * @brief structure for disk component(s)
 */
typedef struct
{
  profile_cfg *cfg;  /**< physical properties of disk component */
  profile *prf;  /**< physical properties of spherical component(s) */
  double *hor;/**< [nest level][NR] array */
  double *ver;/**< [nest level][Nz] array */
  double *node_hor;/**< [nest level][NR + 1] array */
  double *node_ver;/**< [nest level][Nz + 1] array */
  double *pot, *rhoTot;/**< [nest level][NR][Nz] arrays */
  double **rho, **rhoSum, *rho0, *rho1;/**< [Ndisk][nest level][NR][Nz + 1] arrays: Nz + 1 is for rhoSum, to exploit in bisection to determine vertical position */
  double *dPhidR, *d2PhidR2;/**< [nest level][NR][Nz] arrays */
  double *Sigma, *sigmaz, *enc;/**< [Ndisk][nest level][NR] array */
  disk_util util;
  double (*getColumnDensity)(double, double, disk_util);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  double *zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  double *radSph, *rhoSph, *encSph;/**< arrays for spherical averaged profile, [NDISKBIN_RAD] array */
  double *spline_xx, *spline_ff, *spline_f2, *spline_bp;  /**< arrays for spline fit */
  double Rmax, zmax, hh;/**< configuration of domain */
  double invRd;/**< configuration of disk component */
  double logrbin, invlogrbin;/**< configuration of spherical component(s) */
} disk_data;


/* list of functions appeared in ``potdens.c'' */
void   freeDiskProfile
(const int ndisk, disk_data  *disk,
 double  *hor, double  *ver, double  *node_hor, double  *node_ver,
 double  *pot, double  *rho0, double  *rho1, double  *rhoTot, double  *dPhidR, double  *d2PhidR2,
 double  *Sigma, double  *vsigz, double  *enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 double  *zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 double  *radSph, double  *rhoSph, double  *encSph,
 double  *spline_xx, double  *spline_ff, double  *spline_f2, double  *spline_bp);

void allocDiskProfile
(const int ndisk, disk_data **disk, profile_cfg *disk_cfg, int *maxLev, profile **disk_prf, const int skind, const double logrbin, const double invlogrbin,
 double **hor, double **ver, double **node_hor, double **node_ver,
 double **pot, double **rho0, double **rho1, double **rhoTot,
 double **dPhidR, double **d2PhidR2,
 double **Sigma, double **vsigz, double **enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 double **zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 double **radSph, double **rhoSph, double **encSph,
 double **spline_xx, double **spline_ff, double **spline_f2, double **spline_bp);

void makeDiskPotentialTable(const int ndisk, const int maxLev, disk_data * restrict disk);


#endif//POTDENS_H
