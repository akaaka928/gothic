/**
 * @file plot.distribution.c
 *
 * @brief Plot code for radial profiles of multiple components
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2018/02/12 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#define DUMP_DISKHEIGHT_PROFILE

/* #define PLOT_VELOCITY_DISPERSION */
/* /\* #define PLOT_TOOMRE_Q *\/ */

#define DUMP_SPLITTED_SNAPSHOT
/* #define PLOT_VERTICAL_PROFILE */
#define USE_RMS_THICKNESS
#define USE_MEDIAN_FOR_POSITION
#define USE_ITERATED_THICKNESS

#define OUTPUT_BULK_MOTION

#define OVERPLOT_INITIAL_DISKHEIGHT


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

#ifdef  DUMP_SPLITTED_SNAPSHOT
#include <unistd.h>
#endif//DUMP_SPLITTED_SNAPSHOT

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include <gsl/gsl_integration.h>

#include "macro.h"
#include "myutil.h"
#include "name.h"
#include "constants.h"
#include "mpilib.h"
#include "plplotlib.h"

#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../file/io.h"


#define NMAX_GAUSS_QD (51)
#define NTBL_GAUSS_QD ((NMAX_GAUSS_QD >> 1) + (NMAX_GAUSS_QD & 1))
#define NINTBIN NMAX_GAUSS_QD
real gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];

/* #define PLOT_INDIVIDUAL_PROFILES */
/* #define PLOT_SPLITTED_DISTRIBUTION_MAPS */

#define REDUCE_PARTICLE_DISTRIBUTION_MAP (8192)

extern const double      length2astro;extern const char      length_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double        time2astro;extern const char        time_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double        mass2astro;extern const char        mass_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double     density2astro;extern const char     density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double col_density2astro;extern const char col_density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
#ifdef  PLOT_VELOCITY_DISPERSION
extern const double    velocity2astro;extern const char    velocity_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
#endif//PLOT_VELOCITY_DISPERSION

static int ncrit_base;
static const int Nminimum = 16;

static const int NAllocUnit = 32;



/** this value must be the same with NDISKBIN in ../init/disk.h */
#define NUM_ANALYTIC (256)
#define NUM_HORIZONTAL_BIN (10)

typedef struct
{
  ulong head, num;/**< data for N-body particles */
  real *rad, *rho, *enc, *Sig;/**< tables represent the initial profile */
#ifdef  PLOT_VELOCITY_DISPERSION
  real *sigr, *siglos;
#endif//PLOT_VELOCITY_DISPERSION
#ifndef USE_HDF5_FORMAT
  real *psi;/**< tables represent the initial profile */
#endif//USE_HDF5_FORMAT
  double disk_Rd, disk_zd;
  real *disk_radius, *disk_height, *disk_azimuth, *disk_rho, *disk_vol, *disk_rad1d, *disk_Sigma;/**< tables represent the initial profile of the disk component */
#ifdef  PLOT_VELOCITY_DISPERSION
  real *disk_sigR, *disk_sigp, *disk_sigz, *disk_vcirc;
#ifdef  PLOT_TOOMRE_Q
  real *disk_toomre;
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
#ifndef USE_HDF5_FORMAT
  real *disk_pot, *disk_sig;/**< tables represent the initial profile of the disk component */
#endif//USE_HDF5_FORMAT
  real *hor, *ver, *Sigma, *Menc;/**< tables represent the decomposed profile */
  real *ver_pos[NUM_HORIZONTAL_BIN], *ver_rho[NUM_HORIZONTAL_BIN];
  int nrad;/**< number of elements used in tables */
  int disk_nrad, disk_nazi, disk_nsph;/**< number of elements used in tables for the disk component */
  int prf_head, prf_num;/**< data for radial profile of N-body particles */
  int prf_hor_head, prf_hor_num;/**< data for horizontal profile of N-body particles */
  int prf_ver_head[NUM_HORIZONTAL_BIN], prf_ver_num[NUM_HORIZONTAL_BIN];/**< data for vertical profile of N-body particles */
} model;

typedef struct
{
  ulong idx;
  real  x,  y,  z;
  real vx, vy, vz;
  real ax, ay, az;
  real m, pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  real ax_ext, ay_ext, az_ext, pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
} nbody_particle;


static inline void allocProfileArray(int num, real **rad, real **rho, real **enc
#ifdef  PLOT_VELOCITY_DISPERSION
				     , real **sig
#endif//PLOT_VELOCITY_DISPERSION
				     )
{
  __NOTE__("%s\n", "start");

  *rad = (real *)malloc(sizeof(real) * num);
  *rho = (real *)malloc(sizeof(real) * num);
  *enc = (real *)malloc(sizeof(real) * num);
  if( (*rad == NULL) || (*rho == NULL) || (*enc == NULL) ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }
#ifdef  PLOT_VELOCITY_DISPERSION
  *sig = (real *)malloc(sizeof(real) * num);
  if( *sig == NULL ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }
#endif//PLOT_VELOCITY_DISPERSION

  __NOTE__("%s\n", "end");
}

static inline void enlargeProfileArray(int num, real **rad, real **rho, real **enc
#ifdef  PLOT_VELOCITY_DISPERSION
				     , real **sig
#endif//PLOT_VELOCITY_DISPERSION
				       )
{
  __NOTE__("%s\n", "start");

  *rad = realloc(*rad, sizeof(real) * num);
  *rho = realloc(*rho, sizeof(real) * num);
  *enc = realloc(*enc, sizeof(real) * num);
  if( (*rad == NULL) || (*rho == NULL) || (*enc == NULL) ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }
#ifdef  PLOT_VELOCITY_DISPERSION
  *sig = realloc(*sig, sizeof(real) * num);
  if( *sig == NULL ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }
#endif//PLOT_VELOCITY_DISPERSION

  __NOTE__("%s\n", "end");
}


static inline void allocDecomposedProfileArray
(int num, real **pos, real **rho, const bool height, real **dsp
#ifdef  PLOT_VELOCITY_DISPERSION
    , const bool vdisp, real **sigR, real **sigp, real **sigz
#ifdef  PLOT_TOOMRE_Q
    , real **toomre
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
)
{
  __NOTE__("%s\n", "start");

  *pos = (real *)malloc(sizeof(real) * num);
  *rho = (real *)malloc(sizeof(real) * num);
  if( (*pos == NULL) || (*rho == NULL) ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }

  if( height ){
    *dsp = (real *)malloc(sizeof(real) * num);
    if( *dsp == NULL ){
      __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
    }
  }/* if( height ){ */

#ifdef  PLOT_VELOCITY_DISPERSION
  if( vdisp ){
    *sigR   = (real *)malloc(sizeof(real) * num);
    *sigp   = (real *)malloc(sizeof(real) * num);
    *sigz   = (real *)malloc(sizeof(real) * num);
    if( (*sigR == NULL) || (*sigp == NULL) || (*sigz == NULL) ){
      __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
    }
#ifdef  PLOT_TOOMRE_Q
    *toomre = (real *)malloc(sizeof(real) * num);
    if( *toomre == NULL ){
      __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
    }
#endif//PLOT_TOOMRE_Q
  }
#endif//PLOT_VELOCITY_DISPERSION

  __NOTE__("%s\n", "end");
}

static inline void enlargeDecomposedProfileArray
  (int num, real **pos, real **rho, const bool height, real **dsp
#ifdef  PLOT_VELOCITY_DISPERSION
    , const bool vdisp, real **sigR, real **sigp, real **sigz
#ifdef  PLOT_TOOMRE_Q
    , real **toomre
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
    )
{
  __NOTE__("%s\n", "start");

  *pos = realloc(*pos, sizeof(real) * num);
  *rho = realloc(*rho, sizeof(real) * num);
  if( (*pos == NULL) || (*rho == NULL) ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }

  if( height ){
    *dsp = realloc(*dsp, sizeof(real) * num);
    if( *dsp == NULL ){
      __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
    }
  }/* if( height ){ */

#ifdef  PLOT_VELOCITY_DISPERSION
  if( vdisp ){
    *sigR   = realloc(*  sigR, sizeof(real) * num);
    *sigp   = realloc(*  sigp, sizeof(real) * num);
    *sigz   = realloc(*  sigz, sizeof(real) * num);
    if( (*sigR == NULL) || (*sigp == NULL) || (*sigz == NULL) ){
      __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
    }
#ifdef  PLOT_TOOMRE_Q
    *toomre = realloc(*toomre, sizeof(real) * num);
    if( *toomre == NULL ){
      __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
    }
#endif//PLOT_TOOMRE_Q
  }
#endif//PLOT_VELOCITY_DISPERSION

  __NOTE__("%s\n", "end");
}


#ifdef __ICC
/* Disable ICC's remark #161: unrecognized #pragma */
#     pragma warning (disable:161)
#endif//__ICC

int idxAscendingOrder(const void *a, const void *b);
int idxAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((const nbody_particle *)a)->idx > ((const nbody_particle *)b)->idx ){    return ( 1);  }
  else{    if( ((const nbody_particle *)a)->idx < ((const nbody_particle *)b)->idx ){    return (-1);  }
    else{                                                                    return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

int radAscendingOrder(const void *a, const void *b);
int radAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((const nbody_particle *)a)->ax > ((const nbody_particle *)b)->ax ){    return ( 1);  }
  else{    if( ((const nbody_particle *)a)->ax < ((const nbody_particle *)b)->ax ){    return (-1);  }
    else{                                                                  return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

int horAscendingOrder(const void *a, const void *b);
int horAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((const nbody_particle *)a)->ay > ((const nbody_particle *)b)->ay ){    return ( 1);  }
  else{    if( ((const nbody_particle *)a)->ay < ((const nbody_particle *)b)->ay ){    return (-1);  }
    else{                                                                  return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

int verAscendingOrder(const void *a, const void *b);
int verAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((const nbody_particle *)a)->az > ((const nbody_particle *)b)->az ){    return ( 1);  }
  else{    if( ((const nbody_particle *)a)->az < ((const nbody_particle *)b)->az ){    return (-1);  }
    else{                                                                  return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

#ifdef __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC


void plotDistributionMaps
(ulong Np, const int ngroup, model *group, nbody_particle *body, PLFLT time,
 PLplotPltRange xybox, PLplotPltRange xzbox, PLplotPltRange zybox,
 char file[], const int filenum, int argc, char **argv);

void plotRadialProfile
(const int ntot, const int ngroup, const bool overlay_initial, model *group, PLFLT time,
 real *rad, real *rho, real *enc,
#ifdef  PLOT_VELOCITY_DISPERSION
 real *sig, PLplotPltRange velbox,
#endif//PLOT_VELOCITY_DISPERSION
 PLplotPltRange rhobox, PLplotPltRange encbox,
 char *file, const int filenum, int argc, char **argv);

void plotHorizontalProfile
(const int ntot, const int ngroup, const bool overlay_initial, model *group, PLFLT time,
 real *rad, real *rho, PLplotPltRange rhobox,
 const int nspheroids, real *zdisp, PLplotPltRange zbox,
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
 const int ntot_t0, int * restrict hor_head_t0, int * restrict hor_num_t0, real * restrict rad_t0, real * restrict zdisp_t0,
#endif//OVERPLOT_INITIAL_DISKHEIGHT
#ifdef  PLOT_VELOCITY_DISPERSION
 real * restrict sigR, real * restrict sigp, real * restrict sigz, PLplotPltRange Rvelrange,
#ifdef  PLOT_TOOMRE_Q
 real * restrict toomre, PLplotPltRange Qrange,
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
 char *file, const int filenum, int argc, char **argv);
void plotVerticalProfile
(const int ntot, const int ngroup, const bool overlay_initial, model *group, PLFLT time,
 real *rad, real *rho, PLplotPltRange rhobox,
 char *file, const int filenum, int argc, char **argv);

void analyzeRadialProfile
(int Np, nbody_particle *body_tot, const int kind, model *group, int *num, int *rem, real **rad, real **rho, real **enc,
#ifdef  PLOT_VELOCITY_DISPERSION
 real **sigr,
#endif//PLOT_VELOCITY_DISPERSION
 const real r2max);
void analyzeDecomposedProfile
(nbody_particle *body_tot, const int kind, model *group,
 int *numHor, int *remHor, real **hor_pos, real **hor_rho, real **hor_zdisp,
 int *numVer, int *remVer, real **ver_pos, real **ver_rho
#ifdef  PLOT_VELOCITY_DISPERSION
 , real **sigR, real **sigp, real **sigz
#ifdef  PLOT_TOOMRE_Q
 , real **toomre
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
);

void integrateColumnDensity(const int num, const real logRbin, real *rad, real *Sigma, real *mass);


static inline int bisec4nested(const real rad, const int num, real *val, const int levMin, int *lev, real *ratio)
{
  int ll = 0;
  int rr = num - 1;

  *lev = 0;
  for(int ii = levMin; ii >= 0; ii--)
    if( rad <= (val[INDEX2D(maxLev, num, ii, 0)] + val[INDEX2D(maxLev, num, ii, num - 1)]) ){
      *lev = ii;
      break;
    }

  if( rad < val[INDEX2D(maxLev, num, *lev, ll)] ){    *ratio =  ZERO;    return (ll    );  }
  if( rad > val[INDEX2D(maxLev, num, *lev, rr)] ){    *ratio = UNITY;    return (rr - 1);  }

  while( true ){
    const uint cc = ((uint)ll + (uint)rr) >> 1;

    if( (val[INDEX2D(maxLev, num, *lev, cc)] - rad) * (val[INDEX2D(maxLev, num, *lev, ll)] - rad) <= ZERO )      rr = (int)cc;
    else                                                 ll = (int)cc;

    if( (1 + ll) == rr ){
      *ratio = (rad - val[INDEX2D(maxLev, num, *lev, ll)]) / (val[INDEX2D(maxLev, num, *lev, ll + 1)] - val[INDEX2D(maxLev, num, *lev, ll)]);
      return (ll);
    }
  }
}

static inline int bisec(const real rad, const int nmax, real *val)
{
  int ll = 0;
  int rr = nmax - 1;

  if( FABS(val[ll] - rad) / rad < EPSILON )    return (ll    );
  if( FABS(val[rr] - rad) / rad < EPSILON )    return (rr - 1);

  while( true ){
    const uint cc = ((uint)ll + (uint)rr) >> 1;

    if( (val[cc] - rad) * (val[ll] - rad) <= ZERO )      rr = (int)cc;
    else                                                 ll = (int)cc;
    if( (1 + ll) == rr )
      return (ll);
  }
}

static inline void findIdx(const real rad, model prf, int *ll, int *rr)
{
  bool bisection = true;
  *ll = 0;
  *rr = prf.nrad - 1;

  if( bisection == true )    if( fabs(prf.rad[*ll] - rad) / rad < EPSILON ){      bisection = false;      *rr = (*ll) + 1;    }
  if( bisection == true )    if( fabs(prf.rad[*rr] - rad) / rad < EPSILON ){      bisection = false;      *ll = (*rr) - 1;    }

  while( bisection ){
    const uint cc = ((uint)(*ll) + (uint)(*rr)) >> 1;

    if( (prf.rad[cc] - rad) * (prf.rad[*ll] - rad) <= ZERO )      *rr = (int)cc;
    else                                                          *ll = (int)cc;

    if( (1 + (*ll)) == (*rr) )      break;
  }
}

static inline void getInterpolatedDensity(const real RR, const real zz, const int skind, model *group, real *val)
{
  const real rad = SQRT(RR * RR + zz * zz);

  int ll, rr;
  findIdx(rad, group[0], &ll, &rr);

  /* based on linear interpolation */
  const real alpha = (rad - group[0].rad[ll]) / (group[0].rad[rr] - group[0].rad[ll]);
  for(int kk = 0; kk < skind; kk++)
    val[kk] = (UNITY - alpha) * group[kk].rho[ll] + alpha * group[kk].rho[rr];
}

void gaussQuadVertical(const real RR, const int num, const real zmin, const real zmax, const int skind, model *group, real *sum, real *fm, real *fp);
void gaussQuadVertical(const real RR, const int num, const real zmin, const real zmax, const int skind, model *group, real *sum, real *fm, real *fp)
{
  const real mns = HALF * (zmax - zmin);
  const real pls = HALF * (zmax + zmin);

  for(int kk = 0; kk < skind; kk++)    sum[kk] = ZERO;

  if( num & 1 ){
    const real weight =             gsl_gaussQD_weight[(num >> 1)];
    const real  value = pls + mns * gsl_gaussQD_pos   [(num >> 1)];
    getInterpolatedDensity(RR, value, skind, group, fm);
    for(int kk = 0; kk < skind; kk++)
      sum[kk] = weight * fm[kk];
  }

  for(int ii = (num >> 1) - 1; ii >= 0; ii--){
    const real weight = gsl_gaussQD_weight[ii];
    const real zp = pls + mns * gsl_gaussQD_pos[ii];
    const real zm = pls - mns * gsl_gaussQD_pos[ii];
    getInterpolatedDensity(RR, zp, skind, group, fp);
    getInterpolatedDensity(RR, zm, skind, group, fm);

    for(int kk = 0; kk < skind; kk++)
      sum[kk] += weight * (fm[kk] + fp[kk]);
  }

  for(int kk = 0; kk < skind; kk++)
    sum[kk] *= mns;
}


real R_rhor(real RR, real zz, model prf, real invlogrbin);
real R_rhor(real RR, real zz, model prf, real invlogrbin)
{
  const real rad = SQRT(RR * RR + zz * zz);
  int ll, rr;
  findIdx(rad, prf, &ll, &rr);

  const real alpha = (LOG10(rad / prf.rad[ll])) * invlogrbin;
  return (RR * ((UNITY - alpha) * prf.rho[ll] + alpha * prf.rho[rr]));
}

static inline real _gaussQuad2d(real (*func)(real, real, model, real), const real invlogrbin, const model prf, const real xval, const int ny, const real ymin, const real ymax)
{
  const real mns = HALF * (ymax - ymin);
  const real pls = HALF * (ymax + ymin);

  real sum = ZERO;
  if( ny & 1 )    sum = gsl_gaussQD_weight[(ny >> 1)] * func(xval, pls + mns * gsl_gaussQD_pos[(ny >> 1)], prf, invlogrbin);
  for(int ii = (ny >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (func(xval, pls - mns * gsl_gaussQD_pos[ii], prf, invlogrbin) +
				     func(xval, pls + mns * gsl_gaussQD_pos[ii], prf, invlogrbin));
  sum *= mns;

  return (sum);
}

real gaussQuad2d(real (*func)(real, real, model, real), const real invlogrbin, const model prf,
		 const int nx, const real xmin, const real xmax,
		 const int ny, const real ymin, const real ymax);
real gaussQuad2d(real (*func)(real, real, model, real), const real invlogrbin, const model prf,
		 const int nx, const real xmin, const real xmax,
		 const int ny, const real ymin, const real ymax)
{
  const real mns = HALF * (xmax - xmin);
  const real pls = HALF * (xmax + xmin);

  real sum = ZERO;
  if( nx & 1 )    sum = gsl_gaussQD_weight[(nx >> 1)] * _gaussQuad2d(func, invlogrbin, prf, pls + mns * gsl_gaussQD_pos[(nx >> 1)], ny, ymin, ymax);
  for(int ii = (nx >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (_gaussQuad2d(func, invlogrbin, prf, pls - mns * gsl_gaussQD_pos[ii], ny, ymin, ymax) +
				     _gaussQuad2d(func, invlogrbin, prf, pls + mns * gsl_gaussQD_pos[ii], ny, ymin, ymax));
  sum *= mns;

  return (sum);
}


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);

  /** initialization */
  if( argc < 7 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 7);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -problem=<int>\n");
    __FPRINTF__(stderr, "          -ncrit=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 7 ){ */

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));
  int  problem;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,  "problem", &problem));
  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "ncrit", &ncrit_base));

  modifyArgcArgv4PLplot(&argc, argv, 7);

  /* initialize the table for Gaussian Quadrature provided by GSL */
  for(int ii = 0; ii < NTBL_GAUSS_QD; ii++){
    gsl_gaussQD_pos   [ii] = 0.0;
    gsl_gaussQD_weight[ii] = 0.0;
  }/* for(int ii = 0; ii < NTBL_GAUSS_QD; ii++){ */
  gsl_integration_glfixed_table *tab;
  tab = gsl_integration_glfixed_table_alloc(NINTBIN);
  int max = (NINTBIN >> 1) + (NINTBIN & 1);
  for(int ii = 0; ii < max; ii++){
    gsl_gaussQD_pos   [ii] = (real)(*tab).x[(max - 1) - ii];
    gsl_gaussQD_weight[ii] = (real)(*tab).w[(max - 1) - ii];
  }/* for(int ii = 0; ii < max; ii++){ */
  gsl_integration_glfixed_table_free(tab);


  /** load global settings of particle distribution */
  int last;
  readConfigFile(&last, file);
  int unit;
  ulong Ntot;
  real eta, eps;
  double ft, snapshotInterval, saveInterval;
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file, mpi);
  nbody_particle *body;
  /* allocParticleDataAoS((int)Ntot, &body); */
  body = (nbody_particle *)malloc(sizeof(nbody_particle) * Ntot);
  if( body == NULL ){    __KILL__(stderr, "ERROR: failure to allocate body\n");  }
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
  nbody_hdf5 hdf5;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  real *hdf5_acc_ext, *hdf5_pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  ulong *hdf5_idx;
  allocSnapshotArray
    (&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     &hdf5_acc_ext, &hdf5_pot_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     (int)Ntot, &hdf5);
#else///USE_HDF5_FORMAT
  iparticle ibody;
  ulong *idx;
  position *pos;
  acceleration *acc;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  acceleration *ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  velocity *vel;
  ibody_time *ti;
#else///BLOCK_TIME_STEP
  real *vx, *vy, *vz;
#endif//BLOCK_TIME_STEP
  allocParticleData
    ((int)Ntot, &ibody, &idx, &pos, &acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     &ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
     &vel, &ti
#else///BLOCK_TIME_STEP
     &vx, &vy, &vz
#endif//BLOCK_TIME_STEP
     );
#endif//USE_HDF5_FORMAT


  /** set plot range */
  PLFLT radius;
  double radmin, rhomin, encmin, Sigmamin;  double Rmin =  0.0, zmin = 0.0;
  double radmax, rhomax, encmax, Sigmamax;  double Rmax = 10.0, zmax = 2.0;
#ifdef  PLOT_VELOCITY_DISPERSION
  double velmin = 0.0;
  double velmax = 3.0e+2;
#ifdef  PLOT_TOOMRE_Q
  double Qmin = 0.0;
  double Qmax = 5.0;
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
  switch( problem ){
  case  0:    /* cold collapse of a uniform sphere */
    radius = 3.0;
    radmin = 1.0e-2;    rhomin = 1.0e-6;    encmin = 3.2e-6;    Sigmamin = 1.0e-6;
    radmax = 2.0e+1;    rhomax = 2.0e+2;    encmax = 3.2e+0;    Sigmamax = 2.0e+2;
    break;
  case  1:    /* king sphere with W0 = 3 */
    radius = 5.0;
    radmin = 1.0e-2;    rhomin = 1.0e-6;    encmin = 3.2e-6;    Sigmamin = 1.0e-6;
    radmax = 1.0e+1;    rhomax = 1.0e+0;    encmax = 3.2e+0;    Sigmamax = 1.0e+0;
    break;
  case  2:    /* Hernquist sphere with C = 10 */
    radius = 10.0;
    radmin = 1.0e-2;    rhomin = 1.0e-9;    encmin = 3.2e-5;    Sigmamin = 1.0e-8;
    radmax = 5.0e+1;    rhomax = 1.0e+2;    encmax = 3.2e+0;    Sigmamax = 1.0e+2;
    break;
  case  3:    /* NFW sphere with C = 5 */
    radius = 5.0;
    radmin = 1.0e-2;    rhomin = 1.0e-8;    encmin = 3.2e-5;    Sigmamin = 1.0e-6;
    radmax = 5.0e+1;    rhomax = 1.0e+2;    encmax = 3.2e+0;    Sigmamax = 1.0e+1;
    break;
  case  4:    /* Einasto sphere with C = 10 */
    radius = 10.0;
    radmin = 1.0e-2;    rhomin = 1.0e-9;    encmin = 3.2e-5;    Sigmamin = 1.0e-7;
    radmax = 5.0e+1;    rhomax = 1.0e+2;    encmax = 3.2e+0;    Sigmamax = 1.0e+1;
    break;
  case  5:    /* Plummer sphere with C = 10 */
    radius = 10.0;
    radmin = 1.0e-2;    rhomin = 1.0e-9;    encmin = 3.2e-5;    Sigmamin = 1.0e-7;
    radmax = 5.0e+1;    rhomax = 1.0e+0;    encmax = 3.2e+0;    Sigmamax = 1.0e+1;
    break;
  case  6:    /* Burkert sphere with C = 10 */
    radius = 10.0;
    radmin = 1.0e-2;    rhomin = 1.0e-9;    encmin = 3.2e-5;    Sigmamin = 1.0e-7;
    radmax = 5.0e+1;    rhomax = 1.0e+0;    encmax = 3.2e+0;    Sigmamax = 1.0e+0;
    break;
  case  7:    /* Moore sphere with C = 10 */
    radius = 10.0;
    radmin = 1.0e-2;    rhomin = 1.0e-9;    encmin = 3.2e-5;    Sigmamin = 1.0e-7;
    radmax = 5.0e+1;    rhomax = 1.0e+3;    encmax = 3.2e+0;    Sigmamax = 1.0e+1;
    break;
  case  8:    /* Two-power sphere with C = 10 */
    radius = 10.0;
    radmin = 1.0e-2;    rhomin = 1.0e-7;    encmin = 1.0e-3;    Sigmamin = 1.0e-6;
    radmax = 5.0e+1;    rhomax = 1.0e+4;    encmax = 3.2e+0;    Sigmamax = 1.0e+2;
    break;
  case 10:    /* king sphere with W0 = 3 within an Einasto sphere with C = 10 */
    radius = 10.0;
    radmin = 1.0e-2;    rhomin = 1.0e-9;    encmin = 3.2e-5;    Sigmamin = 1.0e-7;
    radmax = 5.0e+1;    rhomax = 1.0e+2;    encmax = 3.2e+0;    Sigmamax = 1.0e+1;
    break;
  case 11:    /* A galaxy with multiple components */
    radius = 5.0e+1;
    radmin = 1.0e-2;    rhomin = 2.0e-9;    encmin = 3.2e-6;    Sigmamin = 1.0e-7;    Rmin =  0.0;    zmin = 0.0;
    radmax = 4.0e+2;    rhomax = 1.0e+0;    encmax = 3.2e+0;    Sigmamax = 1.0e-1;    Rmax = 30.0;    zmax = 1.0;
    break;
  case 12:    /* A galaxy with multiple components */
    radius = 5.0e+1;
    radmin = 1.0e-2;    rhomin = 2.0e-9;    encmin = 3.2e-6;    Sigmamin = 1.0e-7;    Rmin =  0.0;    zmin = 0.0;
    radmax = 4.0e+2;    rhomax = 1.0e+0;    encmax = 3.2e+0;    Sigmamax = 1.0e-1;    Rmax = 30.0;    zmax = 1.6;
    break;
  case 13:    /* A galaxy with multiple components */
    radius = 5.0e+1;
    radmin = 1.0e-2;    rhomin = 2.0e-9;    encmin = 3.2e-5;    Sigmamin = 1.0e-6;    Rmin =  0.0;    zmin = 0.0;
    radmax = 4.0e+2;    rhomax = 1.0e+0;    encmax = 3.2e+1;    Sigmamax = 1.0e-0;    Rmax = 30.0;    zmax = 2.0;
    break;
  case 20:    /* M31 model determined by Fardal et al. (2007) */
    radius = 1.0e+1;
    radmin = 1.0e-2;    rhomin = 1.0e-7;    encmin = 2.0e-1;    Sigmamin = 4.0e-6;    Rmin =  0.0;    zmin = 0.0;
    radmax = 1.0e+3;    rhomax = 2.0e+4;    encmax = 2.0e+4;    Sigmamax = 5.0e+3;    Rmax = 50.0;    zmax = 1.0;
    break;
  case 22:    /* A trial multi components galaxy model */
    radius = 1.0e+1;
    radmin = 1.0e-2;    rhomin = 1.0e-6;    encmin = 2.0e-1;    Sigmamin = 5.0e-5;    Rmin =  0.0;    zmin = 0.0;
    radmax = 1.0e+3;    rhomax = 2.0e+3;    encmax = 8.0e+4;    Sigmamax = 5.0e+2;    Rmax = 50.0;    zmax = 2.0;
    break;
  case 23:    /* MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008) */
    radius = 1.0e+1;
    radmin = 1.0e-2;    rhomin = 1.0e-6;    encmin = 2.0e-1;    Sigmamin = 5.0e-5;    Rmin =  0.0;    zmin = 0.0;
    radmax = 1.0e+3;    rhomax = 2.0e+3;    encmax = 8.0e+4;    Sigmamax = 5.0e+2;    Rmax = 50.0;    zmax = 2.0;
    break;
  case 24:    /* M31 model (Sofue 2015; Gilbert et al. 2012) */
    radius = 1.0e+1;
    radmin = 1.0e-3;    rhomin = 1.0e-7;    encmin = 2.0e-1;    Sigmamin = 4.0e-6;    Rmin =  0.0;    zmin = 0.0;
    radmax = 1.0e+3;    rhomax = 2.0e+3;    encmax = 2.0e+5;    Sigmamax = 5.0e+1;    Rmax = 50.0;    zmax = 1.0;
    break;
  case 25:    /* A trial multi components galaxy model */
    radius = 1.0e+1;
    radmin = 1.0e-2;    rhomin = 1.0e-6;    encmin = 2.0e-1;    Sigmamin = 5.0e-5;    Rmin =  0.0;    zmin = 0.0;
    radmax = 1.0e+3;    rhomax = 2.0e+3;    encmax = 8.0e+4;    Sigmamax = 5.0e+2;    Rmax = 50.0;    zmax = 2.0;
    break;
  case 26:    /* A trial multi components galaxy model (spherical model) */
    radius = 1.0e+1;
    radmin = 1.0e-1;    rhomin = 1.0e-8;    encmin = 2.0e-1;    Sigmamin = 5.0e-6;    Rmin =  0.0;    zmin = 0.0;
    radmax = 1.0e+3;    rhomax = 2.0e+2;    encmax = 4.0e+4;    Sigmamax = 1.0e+2;    Rmax = 50.0;    zmax = 2.0;
    break;
  case 27:    /* M31 model determined by Fardal et al. (2007) with stellar halo */
    radius = 1.0e+1;
    radmin = 1.0e-2;    rhomin = 1.0e-7;    encmin = 2.0e-1;    Sigmamin = 4.0e-6;    Rmin =  0.0;    zmin = 0.0;
    radmax = 1.0e+3;    rhomax = 2.0e+4;    encmax = 2.0e+4;    Sigmamax = 5.0e+3;    Rmax = 50.0;    zmax = 1.0;
    break;
  case 28:    /* A trial multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk) */
    radius = 1.0e+1;
    radmin = 1.0e-2;    rhomin = 1.0e-6;    encmin = 2.0e-1;    Sigmamin = 5.0e-5;    Rmin =  0.0;    zmin = 0.0;
    radmax = 1.0e+3;    rhomax = 2.0e+3;    encmax = 8.0e+4;    Sigmamax = 5.0e+2;    Rmax = 40.0;    zmax = 1.5;
    break;
  case 29:    /* Multi components galaxy model by Vasiliev & Athanassoula (2015) */
    radius = 1.0e+1;
    radmin = 1.0e-2;    rhomin = 1.0e-7;    encmin = 1.0e-3;    Sigmamin = 1.0e-5;    Rmin = 0.0;    zmin = 0.0;
    radmax = 1.0e+2;    rhomax = 1.0e+2;    encmax = 9.9e+1;    Sigmamax = 4.0e-0;    Rmax = 8.0;    zmax = 0.2;
    break;
  case 30:    /* time evolution of MW/A defined in Kuijken & Dubinski (1995) */
    radius = 1.0e+1;
    radmin = 1.0e-3;    rhomin = 1.0e-5;    encmin = 1.0e-3;    Sigmamin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    radmax = 5.0e+1;    rhomax = 1.0e+4;    encmax = 1.0e+2;    Sigmamax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 31:    /* time evolution of MW/B defined in Kuijken & Dubinski (1995) */
    radius = 1.0e+1;
    radmin = 1.0e-3;    rhomin = 1.0e-5;    encmin = 1.0e-3;    Sigmamin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    radmax = 5.0e+1;    rhomax = 1.0e+4;    encmax = 1.0e+2;    Sigmamax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 32:    /* time evolution of MW/C defined in Kuijken & Dubinski (1995) */
    radius = 1.0e+1;
    radmin = 1.0e-3;    rhomin = 1.0e-5;    encmin = 1.0e-3;    Sigmamin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    radmax = 5.0e+1;    rhomax = 1.0e+3;    encmax = 1.0e+3;    Sigmamax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 33:    /* time evolution of MW/D defined in Kuijken & Dubinski (1995) */
    radius = 1.0e+1;
    radmin = 1.0e-3;    rhomin = 1.0e-6;    encmin = 1.0e-3;    Sigmamin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    radmax = 1.0e+2;    rhomax = 1.0e+3;    encmax = 1.0e+3;    Sigmamax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 34:    /* time evolution of M31/A defined in Widrow et al. (2003) */
    radius = 1.0e+1;
    radmin = 1.0e-3;    rhomin = 1.0e-5;    encmin = 1.0e-2;    Sigmamin = 1.0e-4;    Rmin = 0.0;    zmin = 0.00;
    radmax = 1.0e+2;    rhomax = 1.0e+4;    encmax = 1.0e+2;    Sigmamax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 35:    /* time evolution of M31/D defined in Widrow et al. (2003) */
    radius = 1.0e+1;
    radmin = 1.0e-3;    rhomin = 1.0e-5;    encmin = 1.0e-2;    Sigmamin = 1.0e-4;    Rmin = 0.0;    zmin = 0.0;
    radmax = 1.0e+2;    rhomax = 1.0e+4;    encmax = 1.0e+2;    Sigmamax = 1.0e+3;    Rmax = 2.0;    zmax = 0.2;
    break;
  case 36:    /* time evolution of MWa defined in Widrow & Dubinski (2005) */
    radius = 5.0e+1;
    radmin = 1.0e-3;    rhomin = 1.0e-6;    encmin = 1.0e-3;    Sigmamin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    radmax = 1.0e+2;    rhomax = 1.0e+5;    encmax = 1.0e+3;    Sigmamax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 37:    /* time evolution of MWb defined in Widrow & Dubinski (2005) */
    radius = 5.0e+1;
    radmin = 1.0e-3;    rhomin = 1.0e-6;    encmin = 1.0e-3;    Sigmamin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    radmax = 1.0e+2;    rhomax = 1.0e+5;    encmax = 1.0e+3;    Sigmamax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 38:    /* MW like galaxy (based on Widrow et al. 2008, Bedorf et al. 2014) */
    radius = 5.0e+1;
    radmin = 1.0e-2;    rhomin = 1.0e-6;    encmin = 1.0e-2;    Sigmamin = 1.0e-5;    Rmin =  0.0;    zmin = 0.0;
    radmax = 1.0e+3;    rhomax = 1.0e+4;    encmax = 1.0e+5;    Sigmamax = 1.0e+3;    Rmax = 28.0;    zmax = 0.8;
    break;
  case 40:    /* Plummer sphere in table form with C = 20 */
    radius = 20.0;
    radmin = 1.0e-2;    rhomin = 1.0e-7;    encmin = 3.2e-3;    Sigmamin = 1.0e-7;
    radmax = 5.0e+1;    rhomax = 1.0e+1;    encmax = 3.2e+1;    Sigmamax = 1.0e+1;
    break;
  case 41:    /* Two-power sphere in table form with C = 20 */
    radius = 20.0;
    radmin = 1.0e-2;    rhomin = 1.0e-7;    encmin = 1.0e-2;    Sigmamin = 1.0e-6;
    radmax = 5.0e+1;    rhomax = 1.0e+4;    encmax = 3.2e+1;    Sigmamax = 1.0e+3;
    break;
  case 42:    /* de Vaucouleurs sphere in table form */
    radius = 20.0;
    radmin = 1.0e-2;    rhomin = 1.0e-5;    encmin = 1.0e-2;    Sigmamin = 1.0e-4;
    radmax = 5.0e+1;    rhomax = 1.0e+1;    encmax = 3.2e+3;    Sigmamax = 1.0e+1;
    break;
  case 50:    /* de Vaucouleurs sphere */
    radius = 20.0;
    radmin = 1.0e-2;    rhomin = 1.0e-5;    encmin = 1.0e-2;    Sigmamin = 1.0e-4;
    radmax = 5.0e+1;    rhomax = 1.0e+1;    encmax = 3.2e+3;    Sigmamax = 1.0e+1;
    break;
  case 51:    /* projected Two-power model */
    radius = 20.0;
    radmin = 1.0e-2;    rhomin = 1.0e-5;    encmin = 3.2e-3;    Sigmamin = 1.0e-4;
    radmax = 5.0e+1;    rhomax = 1.0e+3;    encmax = 3.2e+1;    Sigmamax = 1.0e+3;
    break;
  case 80:    /* A trial multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk) */
    radius = 1.0e+1;
    radmin = 1.0e-2;    rhomin = 1.0e-6;    encmin = 2.0e-1;    Sigmamin = 5.0e-5;    Rmin =  0.0;    zmin = 0.0;
    radmax = 1.0e+3;    rhomax = 2.0e+3;    encmax = 8.0e+4;    Sigmamax = 1.0e+2;    Rmax = 25.0;    zmax = 1.5;
    break;
  default:
    radius = 3.0;
    radmin = 1.0e-2;    rhomin = 1.0e-6;    encmin = 3.2e-6;    Sigmamin = 1.0e-6;
    radmax = 2.0e+1;    rhomax = 2.0e+2;    encmax = 3.2e+0;    Sigmamax = 2.0e+2;
    break;
  }
#   if  !defined(DUMP_SPLITTED_SNAPSHOT) && !defined(USE_HDF5_FORMAT)
  PLplotPltRange xybox, xzbox, zybox;
  xybox.xmin = xzbox.xmin = -radius * (PLFLT)length2astro;  xybox.xmax = xzbox.xmax =  radius * (PLFLT)length2astro;
  xybox.ymin = zybox.ymin = -radius * (PLFLT)length2astro;  xybox.ymax = zybox.ymax =  radius * (PLFLT)length2astro;
  xzbox.ymin = zybox.xmin = -radius * (PLFLT)length2astro;  xzbox.ymax = zybox.xmax =  radius * (PLFLT)length2astro;
  xybox.xlog = xzbox.xlog = zybox.xlog = LINEAR_PLOT;  xybox.xgrd = xzbox.xgrd = zybox.xgrd = false;
  xybox.ylog = xzbox.ylog = zybox.ylog = LINEAR_PLOT;  xybox.ygrd = xzbox.ygrd = zybox.ygrd = false;
#endif//!defined(DUMP_SPLITTED_SNAPSHOT) && !defined(USE_HDF5_FORMAT)

  int num = 0;
  int rem = NAllocUnit;
  real *rad, *rho, *enc;
#ifdef  PLOT_VELOCITY_DISPERSION
  real *sigr;
#endif//PLOT_VELOCITY_DISPERSION
  allocProfileArray(rem, &rad, &rho, &enc
#ifdef  PLOT_VELOCITY_DISPERSION
		    , &sigr
#endif//PLOT_VELOCITY_DISPERSION
		    );
  int num_hor = 0;
  int rem_hor = NAllocUnit;
  real *hor_pos, *hor_rho, *hor_zdisp;
#ifdef  PLOT_VELOCITY_DISPERSION
  real *sigR, *sigp, *sigz;
#ifdef  PLOT_TOOMRE_Q
  real *toomre;
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
  allocDecomposedProfileArray(rem_hor, &hor_pos, &hor_rho,  true, &hor_zdisp
#ifdef  PLOT_VELOCITY_DISPERSION
			      , true, &sigR, &sigp, &sigz
#ifdef  PLOT_TOOMRE_Q
			      , &toomre
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
			      );
  int num_ver = 0;
  int rem_ver = NAllocUnit;
  real *ver_pos, *ver_rho;
  allocDecomposedProfileArray(rem_ver, &ver_pos, &ver_rho, false, NULL
#ifdef  PLOT_VELOCITY_DISPERSION
			      , false, NULL, NULL, NULL
#ifdef  PLOT_TOOMRE_Q
			      , NULL
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
			      );
  PLplotPltRange rhorange, encrange;
  rhorange.xmin = (PLFLT)log10(radmin *  length2astro);  encrange.xmin = (PLFLT)log10(radmin * length2astro);
  rhorange.xmax = (PLFLT)log10(radmax *  length2astro);  encrange.xmax = (PLFLT)log10(radmax * length2astro);
  rhorange.ymin = (PLFLT)log10(rhomin * density2astro);  encrange.ymin = (PLFLT)log10(encmin *   mass2astro);
  rhorange.ymax = (PLFLT)log10(rhomax * density2astro);  encrange.ymax = (PLFLT)log10(encmax *   mass2astro);
  rhorange.xlog = encrange.xlog = LOGARITHMIC_PLOT;  rhorange.xgrd = encrange.xgrd = true;
  rhorange.ylog = encrange.ylog = LOGARITHMIC_PLOT;  rhorange.ygrd = encrange.ygrd = true;
  PLplotPltRange Sigmarange;
  Sigmarange.xmin = (PLFLT)log10(  radmin *      length2astro);
  Sigmarange.xmax = (PLFLT)log10(  radmax *      length2astro);
  Sigmarange.ymin = (PLFLT)log10(Sigmamin * col_density2astro);
  Sigmarange.ymax = (PLFLT)log10(Sigmamax * col_density2astro);
  Sigmarange.xlog = LOGARITHMIC_PLOT;  Sigmarange.xgrd = true;
  Sigmarange.ylog = LOGARITHMIC_PLOT;  Sigmarange.ygrd = true;
  PLplotPltRange heightrange;
  heightrange.xmin = (PLFLT)Rmin * (PLFLT)length2astro;  heightrange.ymin = (PLFLT)zmin * (PLFLT)length2astro;
  heightrange.xmax = (PLFLT)Rmax * (PLFLT)length2astro;  heightrange.ymax = (PLFLT)zmax * (PLFLT)length2astro;
  heightrange.xlog = heightrange.ylog = LINEAR_PLOT;  heightrange.xgrd = heightrange.ygrd = true;
#ifdef  PLOT_VELOCITY_DISPERSION
  PLplotPltRange rvelrange;
  rvelrange.xmin = (PLFLT)log10(radmin * length2astro);
  rvelrange.xmax = (PLFLT)log10(radmax * length2astro);
  rvelrange.ymin = (PLFLT)velmin * (PLFLT)velocity2astro;
  rvelrange.ymax = (PLFLT)velmax * (PLFLT)velocity2astro;
  rvelrange.xlog = LOGARITHMIC_PLOT;  rvelrange.xgrd = true;
  rvelrange.ylog =      LINEAR_PLOT;  rvelrange.ygrd = true;
  PLplotPltRange Rvelrange;
  Rvelrange.xmin = (PLFLT)Rmin * (PLFLT)length2astro;
  Rvelrange.xmax = (PLFLT)Rmax * (PLFLT)length2astro;
  Rvelrange.ymin = (PLFLT)velmin * (PLFLT)velocity2astro;
  Rvelrange.ymax = (PLFLT)velmax * (PLFLT)velocity2astro;
  Rvelrange.xlog = LINEAR_PLOT;  Rvelrange.xgrd = true;
  Rvelrange.ylog = LINEAR_PLOT;  Rvelrange.ygrd = true;
#ifdef  PLOT_TOOMRE_Q
  PLplotPltRange Qrange;
  Qrange.xmin = (PLFLT)Rmin * (PLFLT)length2astro;
  Qrange.xmax = (PLFLT)Rmax * (PLFLT)length2astro;
  Qrange.ymin = (PLFLT)Qmin;
  Qrange.ymax = (PLFLT)Qmax;
  Qrange.xlog = LINEAR_PLOT;  Qrange.xgrd = true;
  Qrange.ylog = LINEAR_PLOT;  Qrange.ygrd = true;
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION


  /** read analytic profile and analyze */
#ifndef SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  USE_HDF5_FORMAT
#ifdef  DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION
#endif//USE_HDF5_FORMAT
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  int kind = 1;
  int skind = 1;
  model *group;
  bool multi_group = false;
  if( problem >= 1 ){
    multi_group = true;

    FILE *fp;
    char filename[256];
    sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
    fp = fopen(filename, "r");
    if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

    fscanf(fp, "%*d");/* skip reading unit */
    fscanf(fp, "%d\t%d", &kind, &skind);
    group = (model *)malloc(kind * sizeof(model));
    if( group == NULL ){      __KILL__(stderr, "ERROR: failure to allocate group\n");    }
    for(int ii = 0; ii < kind; ii++)
      fscanf(fp, "%zu", &group[ii].num);

    fclose(fp);
  }/* if( problem >= 1 ){ */

  if( !multi_group ){
    group = (model *)malloc(kind * sizeof(model));
    if( group == NULL ){      __KILL__(stderr, "ERROR: failure to allocate group\n");    }

    group[0].num = Ntot;
  }/* if( !multi_group ){ */

  group[0].head = 0;
  for(int ii = 1; ii < kind; ii++)
    group[ii].head = group[ii - 1].head + group[ii - 1].num;

  /** load initial profile */
  bool overlay_initial = false;
#ifndef SET_EXTERNAL_POTENTIAL_FIELD
  if( problem >= 1 ){
    overlay_initial = true;

    for(int ii = 0; ii < kind; ii++){
      /* read radial profile of volume density, enclosed mass and potential */
      char filename[256];

#ifdef  USE_HDF5_FORMAT
      /* open the target file and group */
      sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "profile");
      hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      /* read flag about DOUBLE_PRECISION */
      int useDP;
      hid_t attribute = H5Aopen(target, "useDP", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
      chkHDF5err(H5Aclose(attribute));

      char grp[16];      sprintf(grp, "data%d", ii);
      hid_t h5group = H5Gopen(target, grp, H5P_DEFAULT);

      /* read attribute data */
      attribute = H5Aopen(h5group, "num", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &group[ii].nrad));
      chkHDF5err(H5Aclose(attribute));
      /* memory allocation */
      group[ii].rad = (real *)malloc(group[ii].nrad * sizeof(real));
      group[ii].rho = (real *)malloc(group[ii].nrad * sizeof(real));
      group[ii].enc = (real *)malloc(group[ii].nrad * sizeof(real));
      group[ii].Sig = (real *)malloc(group[ii].nrad * sizeof(real));
      if( group[ii].rad == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].rad\n", ii);      }
      if( group[ii].rho == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].rho\n", ii);      }
      if( group[ii].enc == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].enc\n", ii);      }
      if( group[ii].Sig == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].Sig\n", ii);      }
#ifdef  PLOT_VELOCITY_DISPERSION
      group[ii].sigr   = (real *)malloc(group[ii].nrad * sizeof(real));
      group[ii].siglos = (real *)malloc(group[ii].nrad * sizeof(real));
      if( group[ii].sigr   == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].sigr\n", ii);      }
      if( group[ii].siglos == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].siglos\n", ii);      }
#endif//PLOT_VELOCITY_DISPERSION

      /* read the profile data */
      hid_t dataset;
      /* read radius */
      dataset = H5Dopen(h5group, "rad", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].rad));
      chkHDF5err(H5Dclose(dataset));
      /* read density */
      dataset = H5Dopen(h5group, "rho", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].rho));
      chkHDF5err(H5Dclose(dataset));
      /* read enclosed mass */
      dataset = H5Dopen(h5group, "enc", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].enc));
      chkHDF5err(H5Dclose(dataset));
      /* read column density */
      dataset = H5Dopen(h5group, "Sigma", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].Sig));
      chkHDF5err(H5Dclose(dataset));
#ifdef  PLOT_VELOCITY_DISPERSION
      /* read 1-dimensional velocity dispersion */
      dataset = H5Dopen(h5group, "sigma_r", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].sigr));
      chkHDF5err(H5Dclose(dataset));
      /* read velocity dispersion along the line-of-sight */
      dataset = H5Dopen(h5group, "sigma_los", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].siglos));
      chkHDF5err(H5Dclose(dataset));
#endif//PLOT_VELOCITY_DISPERSION

      /* close the file */
      chkHDF5err(H5Gclose(h5group));
      chkHDF5err(H5Fclose(target));

#ifdef  DOUBLE_PRECISION
      if( useDP != 1 ){
	__KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);
      }/* if( useDP != 1 ){ */
#else///DOUBLE_PRECISION
      if( useDP != 0 ){
	__KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);
      }/* if( useDP != 0 ){ */
#endif//DOUBLE_PRECISION

#else///USE_HDF5_FORMAT

      FILE *fp;
      sprintf(filename, "%s/%s.profile.%d.dat", DATAFOLDER, file, ii);
      fp = fopen(filename, "rb");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);      }

      bool success = true;
      success &= (fread(&group[ii].nrad, sizeof(int), 1, fp) == 1);
      group[ii].rad = (real *)malloc(group[ii].nrad * sizeof(real));
      group[ii].rho = (real *)malloc(group[ii].nrad * sizeof(real));
      group[ii].enc = (real *)malloc(group[ii].nrad * sizeof(real));
      group[ii].psi = (real *)malloc(group[ii].nrad * sizeof(real));
      if( group[ii].rad   == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].rad\n", ii);      }
      if( group[ii].rho   == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].rho\n", ii);      }
      if( group[ii].enc   == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].enc\n", ii);      }
      if( group[ii].psi   == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].psi\n", ii);      }

      success &= (fread(group[ii].rad, sizeof(real), group[ii].nrad, fp) == group[ii].nrad);
      success &= (fread(group[ii].rho, sizeof(real), group[ii].nrad, fp) == group[ii].nrad);
      success &= (fread(group[ii].enc, sizeof(real), group[ii].nrad, fp) == group[ii].nrad);
      success &= (fread(group[ii].psi, sizeof(real), group[ii].nrad, fp) == group[ii].nrad);

      if( !success ){	__KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);      }
      fclose(fp);

#endif//USE_HDF5_FORMAT

      group[ii].hor   = (real *)malloc(NUM_ANALYTIC * sizeof(real));
      group[ii].ver   = (real *)malloc(NUM_ANALYTIC * sizeof(real));
      group[ii].Sigma = (real *)malloc(NUM_ANALYTIC * sizeof(real));
      group[ii].Menc  = (real *)malloc(NUM_ANALYTIC * sizeof(real));
      if( group[ii].hor   == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].hor\n"  , ii);      }
      if( group[ii].ver   == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].ver\n"  , ii);      }
      if( group[ii].Sigma == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].Sigma\n", ii);      }
      if( group[ii].Menc  == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group[%d].Menc\n" , ii);      }

      for(int jj = 0; jj < NUM_HORIZONTAL_BIN; jj++){
	group[ii].ver_pos[jj] = (real *)malloc(NUM_ANALYTIC * sizeof(real));
	group[ii].ver_rho[jj] = (real *)malloc(NUM_ANALYTIC * sizeof(real));
	if( group[ii].ver_pos[jj] == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[%d].ver_pos[%d]\n", ii, jj);	}
	if( group[ii].ver_rho[jj] == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[%d].ver_rho[%d]\n", ii, jj);	}
      }/* for(int jj = 0; jj < NUM_HORIZONTAL_BIN; jj++){ */


      /** read 2D-distribution of volume density and potential for the disk component */
      if( ii >= skind ){
#ifdef  USE_HDF5_FORMAT

	/* read HDF5 file */
	char filename[128];
	sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "disk");
	hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	/* read flag about DOUBLE_PRECISION */
	int _useDP;
	attribute = H5Aopen(target, "useDP", H5P_DEFAULT);
	chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &_useDP));
	chkHDF5err(H5Aclose(attribute));
	int maxLev;
	attribute = H5Aopen(target, "maxLev", H5P_DEFAULT);
	chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &maxLev));
	chkHDF5err(H5Aclose(attribute));

	char grp[16];	sprintf(grp, "data%d", ii - skind);
	hid_t h5group = H5Gopen(target, grp, H5P_DEFAULT);

	/* read attributes */
	hid_t attribute;
	/* read scale radius and scale height */
	attribute = H5Aopen(h5group, "Rs", H5P_DEFAULT);
	chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &group[ii].disk_Rd));
	chkHDF5err(H5Aclose(attribute));
	attribute = H5Aopen(h5group, "zd", H5P_DEFAULT);
	chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &group[ii].disk_zd));
	chkHDF5err(H5Aclose(attribute));

	sprintf(grp, "patch%d", 0);
	hid_t patch = H5Gopen(h5group, grp, H5P_DEFAULT);
	int dims[2];
	/* read # of arrays */
	attribute = H5Aopen(patch, "num", H5P_DEFAULT);
	chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, dims));
	chkHDF5err(H5Aclose(attribute));
	chkHDF5err(H5Gclose(patch));
	group[ii].disk_nrad = dims[0];
	group[ii].disk_nazi = dims[1];

	/* memory allocation */
	group[ii].disk_radius = (real *)malloc( maxLev      *  group[ii].disk_nrad                             * sizeof(real));
	group[ii].disk_height = (real *)malloc( maxLev      *                              group[ii].disk_nazi * sizeof(real));
	group[ii].disk_rho    = (real *)malloc( maxLev      *  group[ii].disk_nrad       * group[ii].disk_nazi * sizeof(real));
	group[ii].disk_rad1d  = (real *)malloc((maxLev + 1) * (group[ii].disk_nrad >> 1)                       * sizeof(real));
	group[ii].disk_Sigma  = (real *)malloc((maxLev + 1) * (group[ii].disk_nrad >> 1)                       * sizeof(real));
	if( group[ii].disk_radius == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_radius\n");	}
	if( group[ii].disk_height == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_height\n");	}
	if( group[ii].disk_rho	  == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_rho\n"   );	}
	if( group[ii].disk_rad1d  == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_rad1d\n" );	}
	if( group[ii].disk_Sigma  == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_Sigma\n" );	}
#ifdef  PLOT_VELOCITY_DISPERSION
	group[ii].disk_sigR   = (real *)malloc((maxLev + 1) * (group[ii].disk_nrad >> 1) * sizeof(real));
	group[ii].disk_sigp   = (real *)malloc((maxLev + 1) * (group[ii].disk_nrad >> 1) * sizeof(real));
	group[ii].disk_sigz   = (real *)malloc((maxLev + 1) * (group[ii].disk_nrad >> 1) * sizeof(real));
	group[ii].disk_vcirc  = (real *)malloc((maxLev + 1) * (group[ii].disk_nrad >> 1) * sizeof(real));
	if( group[ii].disk_sigR   == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_sigR\n");	}
	if( group[ii].disk_sigp   == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_sigp\n");	}
	if( group[ii].disk_sigz   == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_sigz\n");	}
	if( group[ii].disk_vcirc  == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_vcirc\n");	}
#ifdef  PLOT_TOOMRE_Q
	group[ii].disk_toomre = (real *)malloc((maxLev + 1) * (group[ii].disk_nrad >> 1) * sizeof(real));
	if( group[ii].disk_toomre == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_toomre\n");	}
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION

	/* read disk data */
	for(int lev = 0; lev < maxLev; lev++){
	  sprintf(grp, "patch%d", lev);
	  patch = H5Gopen(h5group, grp, H5P_DEFAULT);

	  /* read horizontal position */
	  dataset = H5Dopen(patch, "radius", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(group[ii].disk_radius[INDEX2D(maxLev, group[ii].disk_nrad, lev, 0)])));
	  chkHDF5err(H5Dclose(dataset));
	  /* read vertical position */
	  dataset = H5Dopen(patch, "height", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(group[ii].disk_height[INDEX2D(maxLev, group[ii].disk_nazi, lev, 0)])));
	  chkHDF5err(H5Dclose(dataset));
	  /* read density distribution */
	  dataset = H5Dopen(patch, "rho", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(group[ii].disk_rho[INDEX2D(maxLev, group[ii].disk_nrad * group[ii].disk_nazi, lev, 0)])));
	  chkHDF5err(H5Dclose(dataset));
	  chkHDF5err(H5Gclose(patch));
	}/* for(int lev = 0; lev < maxLev; lev++){ */

	sprintf(grp, "1D_data");
	patch = H5Gopen(h5group, grp, H5P_DEFAULT);
	/* read radius */
	dataset = H5Dopen(patch, "radius", H5P_DEFAULT);
	chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].disk_rad1d));
	chkHDF5err(H5Dclose(dataset));
	/* read column density profile */
	dataset = H5Dopen(patch, "Sigma", H5P_DEFAULT);
	chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].disk_Sigma));
	chkHDF5err(H5Dclose(dataset));
#ifdef  PLOT_VELOCITY_DISPERSION
	/* read velocity dispersion in R-direction */
	dataset = H5Dopen(patch, "sigmaR", H5P_DEFAULT);
	chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].disk_sigR));
	chkHDF5err(H5Dclose(dataset));
	/* read velocity dispersion in tangential direction */
	dataset = H5Dopen(patch, "sigmap", H5P_DEFAULT);
	chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].disk_sigp));
	chkHDF5err(H5Dclose(dataset));
	/* read velocity dispersion in z-direction */
	dataset = H5Dopen(patch, "sigmaz", H5P_DEFAULT);
	chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].disk_sigz));
	chkHDF5err(H5Dclose(dataset));
	/* read circular velocity */
	dataset = H5Dopen(patch, "vcirc", H5P_DEFAULT);
	chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].disk_vcirc));
	chkHDF5err(H5Dclose(dataset));
#ifdef  PLOT_TOOMRE_Q
	/* read Toomre's Q-value */
	dataset = H5Dopen(patch, "Toomre's Q", H5P_DEFAULT);
	chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, group[ii].disk_toomre));
	chkHDF5err(H5Dclose(dataset));
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
	chkHDF5err(H5Gclose(patch));

	/* close the file */
	chkHDF5err(H5Gclose(h5group));
	chkHDF5err(H5Fclose(target));
#ifdef  DOUBLE_PRECISION
	if( _useDP != 1 ){
	  __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", _useDP, true);
	}/* if( useDP != 1 ){ */
#else///DOUBLE_PRECISION
	if( _useDP != 0 ){
	  __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", _useDP, false);
	}/* if( useDP != 0 ){ */
#endif//DOUBLE_PRECISION

#else///USE_HDF5_FORMAT

	sprintf(filename, "%s/%s.diskdat.%d.dat", DATAFOLDER, file, ii - skind);
	fp = fopen(filename, "rb");
	if( fp == NULL ){	  __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);	}

	success = true;
	success &= (fread(&group[ii].disk_nrad, sizeof(   int), 1, fp) == 1);
	success &= (fread(&group[ii].disk_nazi, sizeof(   int), 1, fp) == 1);
	success &= (fread(&group[ii].disk_nsph, sizeof(   int), 1, fp) == 1);
	success &= (fread(&group[ii].disk_Rd  , sizeof(double), 1, fp) == 1);
	success &= (fread(&group[ii].disk_zd  , sizeof(double), 1, fp) == 1);
	group[ii].disk_radius  = (real *)malloc(group[ii].disk_nrad                       * sizeof(real));
	group[ii].disk_height  = (real *)malloc(                      group[ii].disk_nazi * sizeof(real));
	group[ii].disk_rho     = (real *)malloc(group[ii].disk_nrad * group[ii].disk_nazi * sizeof(real));
	group[ii].disk_pot     = (real *)malloc(group[ii].disk_nrad * group[ii].disk_nazi * sizeof(real));
	group[ii].disk_sig     = (real *)malloc(group[ii].disk_nrad                       * sizeof(real));
	group[ii].disk_Sigma   = (real *)malloc(group[ii].disk_nrad                       * sizeof(real));
	if( group[ii].disk_radius  == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_radius\n");	}
	if( group[ii].disk_height  == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_height\n");	}
	if( group[ii].disk_rho	   == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_rho\n"   );	}
	if( group[ii].disk_pot	   == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_pot\n"   );	}
	if( group[ii].disk_sig	   == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_sig\n"   );	}
	if( group[ii].disk_Sigma   == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_Sigma\n" );	}

	success &= (fread(group[ii].disk_radius , sizeof(real), group[ii].disk_nrad                      , fp) == group[ii].disk_nrad                      );
	success &= (fread(group[ii].disk_height , sizeof(real),                       group[ii].disk_nazi, fp) ==                       group[ii].disk_nazi);
	success &= (fread(group[ii].disk_rho    , sizeof(real), group[ii].disk_nrad * group[ii].disk_nazi, fp) == group[ii].disk_nrad * group[ii].disk_nazi);
	success &= (fread(group[ii].disk_pot    , sizeof(real), group[ii].disk_nrad * group[ii].disk_nazi, fp) == group[ii].disk_nrad * group[ii].disk_nazi);
	success &= (fread(group[ii].disk_sig    , sizeof(real), group[ii].disk_nrad                      , fp) == group[ii].disk_nrad                      );
	success &= (fread(group[ii].disk_Sigma  , sizeof(real), group[ii].disk_nrad                      , fp) == group[ii].disk_nrad                      );

	if( !success ){	  __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);	}
	fclose(fp);
#endif//USE_HDF5_FORMAT


	/** create volume density table and column density table of the disk component */
	const real    Rmax = group[ii].disk_radius[group[ii].disk_nrad - 1];
	const real    zmax = group[ii].disk_height[group[ii].disk_nazi - 1];
	const real logRmax = LOG10(Rmax);	const real logRmin = LOG10(group[ii].disk_radius[INDEX2D(maxLev, group[ii].disk_nrad, maxLev - 1, 0)]);
	const real logzmax = LOG10(zmax);	const real logzmin = LOG10(group[ii].disk_height[INDEX2D(maxLev, group[ii].disk_nazi, maxLev - 1, 0)]);
	const real logRbin = (logRmax - logRmin) / (real)(NUM_ANALYTIC - 1);
	const real logzbin = (logzmax - logzmin) / (real)(NUM_ANALYTIC - 1);
	for(int jj = 0; jj < NUM_ANALYTIC; jj++){
	  group[ii].hor[jj] = POW(TEN, logRmin + logRbin * (real)jj);
	  group[ii].ver[jj] = POW(TEN, logzmin + logzbin * (real)jj);
	}/* for(int jj = 0; jj < NUM_ANALYTIC; jj++){ */
	group[ii].disk_vol = (real *)malloc(NUM_ANALYTIC * NUM_ANALYTIC * sizeof(real));
	if( group[ii].disk_vol == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_vol\n");	}

	for(int jj = 0; jj < NUM_ANALYTIC; jj++){
	  const real RR = group[ii].hor[jj];

	  int idx = bisec(RR, (maxLev + 1) * (group[ii].disk_nrad >> 1), group[ii].disk_rad1d);
	  real alpha = ((idx != 0) && (idx != ((maxLev + 1) * (group[ii].disk_nrad >> 1) - 1))) ? ((RR - group[ii].disk_rad1d[idx]) / (group[ii].disk_rad1d[1 + idx] - group[ii].disk_rad1d[idx])) : ((idx == 0) ? (ZERO) : (UNITY));
	  group[ii].Sigma[jj] = (UNITY - alpha) * group[ii].disk_Sigma[idx] + alpha * group[ii].disk_Sigma[idx + 1];

	  int rlev;
	  real arad;
	  int irad = bisec4nested(RR, group[ii].disk_nrad, group[ii].disk_radius, maxLev - 1, &rlev, &arad);

	  for(int kk = 0; kk < NUM_ANALYTIC; kk++){
	    const real zz = group[ii].ver[kk];

	    /* 2D-interpolation */
	    int zlev;
	    real aazi;
	    const int iazi = bisec4nested(zz, group[ii].disk_nazi, group[ii].disk_height, rlev, &zlev, &aazi);
	    if( rlev != zlev ){
	      irad >>= (rlev - zlev);
	      rlev = zlev;
	      arad = (RR - group[ii].disk_radius[INDEX2D(maxLev, group[ii].disk_nrad, rlev, irad)]) / (group[ii].disk_radius[INDEX2D(maxLev, group[ii].disk_nrad, rlev, 1 + irad)] - group[ii].disk_radius[INDEX2D(maxLev, group[ii].disk_nrad, rlev, irad)]);
	    }/* if( rlev != zlev ){ */
	    group[ii].disk_vol[INDEX2D(NUM_ANALYTIC, NUM_ANALYTIC, jj, kk)] =
	      (group[ii].disk_rho[INDEX(maxLev, group[ii].disk_nrad, group[ii].disk_nazi, zlev,     irad,     iazi)] * (UNITY - aazi) +
	       group[ii].disk_rho[INDEX(maxLev, group[ii].disk_nrad, group[ii].disk_nazi, zlev,     irad, 1 + iazi)] *	        aazi) * (UNITY - arad) +
	      (group[ii].disk_rho[INDEX(maxLev, group[ii].disk_nrad, group[ii].disk_nazi, zlev, 1 + irad,     iazi)] * (UNITY - aazi) +
	       group[ii].disk_rho[INDEX(maxLev, group[ii].disk_nrad, group[ii].disk_nazi, zlev, 1 + irad, 1 + iazi)] *	        aazi) *	         arad;
	  }/* for(int kk = 0; kk < NUM_ANALYTIC; kk++){ */
	}/* for(int jj = 0; jj < NUM_ANALYTIC; jj++){ */


	/** create vertical density table of the disk component */
	/** divide into NUM_HORIZONTAL_BIN bins */
	const real alpha = POW(TEN, logRbin);
	const real factor = (UNITY + alpha) * (UNITY + alpha) * (UNITY + alpha) * (alpha - UNITY) / (FOUR * alpha * alpha);
	integrateColumnDensity(NUM_ANALYTIC, logRbin, group[ii].hor, group[ii].Sigma, group[ii].Menc);
	const real Menc_bin = group[ii].Menc[NUM_ANALYTIC - 1] / (real)NUM_HORIZONTAL_BIN;
	int head = 0;
	for(int jj = 0; jj < NUM_HORIZONTAL_BIN; jj++){
	  /* set tail */
	  const real Mmin = group[ii].Menc[head];
	  int tail;
	  for(tail = head + 1; tail < NUM_ANALYTIC; tail++)
	    if( group[ii].Menc[tail] > Mmin + Menc_bin )
	      break;

	  const real Rinner = group[ii].hor[head    ];
	  const real Router = group[ii].hor[tail - 1];
	  const real R2inv  = UNITY / (Router * Router - Rinner * Rinner);

	  for(int kk = 0; kk < NUM_ANALYTIC; kk++){
	    group[ii].ver_pos[jj][kk] = group[ii].ver[kk];
	    group[ii].ver_rho[jj][kk] = 0.0;
	  }
	  for(int kk = head; kk < tail; kk++){
	    const real R2 = factor * R2inv * group[ii].hor[kk] * group[ii].hor[kk];
	    for(int mm = 0; mm < NUM_ANALYTIC; mm++)
	      group[ii].ver_rho[jj][mm] += R2 * group[ii].disk_vol[INDEX2D(NUM_ANALYTIC, NUM_ANALYTIC, kk, mm)];
	  }

	  head = tail;
	}/* for(int jj = 0; jj < NUM_HORIZONTAL_BIN; jj++){ */
      }/* if( ii >= skind ){ */
    }/* for(int ii = 0; ii < kind; ii++){ */


    /** analyze column density profile of spherical components */
    const real logRmin = (real)log10(radmin);
    const real logRmax = (real)log10(radmax);
    const real logRbin = (logRmax - logRmin) / (real)(NUM_ANALYTIC - 1);
    for(int ii = 0; ii < NUM_ANALYTIC; ii++){
      const real RR = POW(TEN, logRmin + logRbin * (real)ii);
      int ll, rr;
      findIdx(RR, group[0], &ll, &rr);
      const real ratio = (RR - group[0].rad[ll]) / (group[0].rad[rr] - group[0].rad[ll]);
      for(int kk = 0; kk < skind; kk++){
    	group[kk].hor  [ii] = RR;
    	group[kk].Sigma[ii] = (UNITY - ratio) * group[kk].Sig[ll] + ratio * group[kk].Sig[rr];
      }/* for(int kk = 0; kk < skind; kk++){ */
    }/* for(int ii = 0; ii < NUM_ANALYTIC; ii++){ */

    for(int kk = 0; kk < skind; kk++)
      integrateColumnDensity(NUM_ANALYTIC, logRbin, group[kk].hor, group[kk].Sigma, group[kk].Menc);


    /** analyze vertical density profile of spherical components */
    for(int kk = 0; kk < skind; kk++){
      const real    logrmin = LOG10(group[kk].rad[0]);
      const real    logrmax = LOG10(group[kk].rad[group[kk].nrad - 1]);
      const real invlogrbin = (real)(group[kk].nrad - 1) / (logrmax - logrmin);

      const real logzmin = (real)log10(radmin);
      const real logzmax = (real)log10(radmax);
      const real logzbin = (logzmax - logzmin) / (real)(NUM_ANALYTIC - 1);
      const real alpha = POW(TEN, logzbin);
      const real factor = (alpha * alpha - UNITY) / (TWO * alpha);

      /** divide into NUM_HORIZONTAL_BIN bins */
      const real Menc_bin = group[kk].Menc[NUM_ANALYTIC - 1] / (real)NUM_HORIZONTAL_BIN;
      int head = 0;
      for(int ii = 0; ii < NUM_HORIZONTAL_BIN; ii++){
	/* set tail */
	const real Mmin = group[kk].Menc[head];
	int tail;
	for(tail = head + 1; tail < NUM_ANALYTIC; tail++)
	  if( group[kk].Menc[tail] > Mmin + Menc_bin )
	    break;

	const real Rinner = group[kk].hor[head    ];
	const real Router = group[kk].hor[tail - 1];
	const real R2inv  = UNITY / (Router * Router - Rinner * Rinner);

	for(int jj = 0; jj < NUM_ANALYTIC; jj++){
	  const real zz = POW(TEN, logzmin + logzbin * (real)jj);
	  const real dz = factor * zz;

	  group[kk].ver_pos[ii][jj] = zz;
	  group[kk].ver_rho[ii][jj] = (TWO * R2inv / dz) * gaussQuad2d(R_rhor, invlogrbin, group[kk], NINTBIN, Rinner, Router, NINTBIN, zz - HALF * dz, zz + HALF * dz);
	}/* for(int jj = 0; jj < NUM_ANALYTIC; jj++){ */
	head = tail;

      }/* for(int ii = 0; ii < NUM_HORIZONTAL_BIN; ii++){ */
    }/* for(int kk = 0; kk < skind; kk++){ */
  }/* if( problem >= 1 ){ */
#endif//SET_EXTERNAL_POTENTIAL_FIELD


  /** read particle distribution and analyze */
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
  int num_hor_t0;
  real *hor_pos_t0, *hor_zdisp_t0;
  int *prf_hor_head_t0, *prf_hor_num_t0;
  prf_hor_head_t0 = (int *)malloc(sizeof(int) * kind);
  prf_hor_num_t0  = (int *)malloc(sizeof(int) * kind);
  if( prf_hor_head_t0 == NULL ){	__KILL__(stderr, "%s\n", "ERROR: failure to allocate prf_hor_head_t0.");      }
  if( prf_hor_num_t0  == NULL ){	__KILL__(stderr, "%s\n", "ERROR: failure to allocate prf_hor_num_t0.");      }
#endif//OVERPLOT_INITIAL_DISKHEIGHT

#ifdef  OUTPUT_BULK_MOTION
  const int nfile = (end - start + 1) / interval;
  real *comPos;  comPos = (real *)malloc(sizeof(real) * 3 * nfile * kind);  if( comPos  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate comPos");  }
  real *comVel;  comVel = (real *)malloc(sizeof(real) * 3 * nfile * kind);  if( comVel  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate comVel");  }
  for(int ii = 0; ii < nfile * 3 * kind; ii++)
    comPos[ii] = comVel[ii] = ZERO;
#endif//OUTPUT_BULK_MOTION

  int ifile = (int)(start + mpi.rank);
  for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){
    double time;
    ulong steps;
    int unit_read;
#ifdef  USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file, (uint)filenum, &hdf5, hdf5type);
    for(int ii = 0; ii < (int)Ntot; ii++){
      body[ii]. x  = hdf5.pos[ii * 3];      body[ii]. y = hdf5.pos[ii * 3 + 1];      body[ii].z   = hdf5.pos[ii * 3 + 2];
      body[ii].vx  = hdf5.vel[ii * 3];      body[ii].vy = hdf5.vel[ii * 3 + 1];      body[ii].vz  = hdf5.vel[ii * 3 + 2];
      body[ii].ax  = hdf5.acc[ii * 3];      body[ii].ay = hdf5.acc[ii * 3 + 1];      body[ii].az  = hdf5.acc[ii * 3 + 2];
      body[ii].idx = hdf5.idx[ii    ];      body[ii]. m = hdf5.  m[ii        ];      body[ii].pot = hdf5.pot[ii        ];
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      body[ii].ax_ext = hdf5.acc_ext[ii * 3];      body[ii].ay_ext = hdf5.acc_ext[ii * 3 + 1];      body[ii].az_ext = hdf5.acc_ext[ii * 3 + 2];
      body[ii].pot_ext = hdf5.pot_ext[ii];
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#else///USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file, ibody, (uint)filenum);
    for(int ii = 0; ii < (int)Ntot; ii++){
      body[ii]. x  = ibody.pos[ii].x;      body[ii]. y = ibody.pos[ii].y;      body[ii]. z  = ibody.pos[ii].z;
      body[ii].vx  = ibody.vel[ii].x;      body[ii].vy = ibody.vel[ii].y;      body[ii].vz  = ibody.vel[ii].z;
      body[ii].ax  = ibody.acc[ii].x;      body[ii].ay = ibody.acc[ii].y;      body[ii].az  = ibody.acc[ii].z;
      body[ii].idx = ibody.idx[ii]  ;      body[ii]. m = ibody.pos[ii].m;      body[ii].pot = ibody.acc[ii].pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      body[ii].ax_ext = ibody.acc_ext[ii].x;      body[ii].ay_ext = ibody.acc_ext[ii].y;      body[ii].az_ext = ibody.acc_ext[ii].z;
      body[ii].pot_ext = ibody.acc_ext[ii].pot;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#endif//USE_HDF5_FORMAT
    if( unit_read != unit ){
      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
    }/* if( unit_read != unit ){ */
    /* sort by particle index */
    qsort(body, Ntot, sizeof(nbody_particle), idxAscendingOrder);

#ifdef  OUTPUT_BULK_MOTION
    for(int kk = 0; kk < kind; kk++){
      double Mtot = 0.0;
      double com[3] = {0.0, 0.0, 0.0};
      double vel[3] = {0.0, 0.0, 0.0};

      for(ulong ii = group[kk].head; ii < group[kk].head + group[kk].num; ii++){
	const double mass = CAST_R2D(body[ii].m);

	Mtot += mass;
	com[0] += mass * CAST_R2D(body[ii].x);	vel[0] += mass * CAST_R2D(body[ii].vx);
	com[1] += mass * CAST_R2D(body[ii].y);	vel[1] += mass * CAST_R2D(body[ii].vy);
	com[2] += mass * CAST_R2D(body[ii].z);	vel[2] += mass * CAST_R2D(body[ii].vz);
      }/* for(ulong ii = group[kk].head; ii < group[kk].head + group[kk].num; ii++){ */

      Mtot = 1.0 / Mtot;
#pragma unroll
      for(int ii = 0; ii < 3; ii++){
	com[ii] *= Mtot;
	vel[ii] *= Mtot;

	comPos[INDEX(nfile, kind, 3, (filenum - start) / interval, kk, ii)] = CAST_D2R(com[ii]);
	comVel[INDEX(nfile, kind, 3, (filenum - start) / interval, kk, ii)] = CAST_D2R(vel[ii]);
      }/* for(int ii = 0; ii < 3; ii++){ */
    }/* for(int kk = 0; kk < kind; kk++){ */
#endif//OUTPUT_BULK_MOTION

#ifdef  DUMP_SPLITTED_SNAPSHOT
#ifdef  USE_HDF5_FORMAT
    if( multi_group ){
      char filename[128];
      sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "split", filenum);
      if( 0 != access(filename, F_OK) ){
	for(int ii = 0; ii < (int)Ntot; ii++){
	  hdf5.idx[ii    ] = body[ii].idx;	hdf5.  m[ii        ] = body[ii]. m;	hdf5.pot[ii        ] = body[ii].pot;
	  hdf5.pos[ii * 3] = body[ii]. x ;	hdf5.pos[ii * 3 + 1] = body[ii]. y;	hdf5.pos[ii * 3 + 2] = body[ii]. z ;
	  hdf5.vel[ii * 3] = body[ii].vx ;	hdf5.vel[ii * 3 + 1] = body[ii].vy;	hdf5.vel[ii * 3 + 2] = body[ii].vz ;
	  hdf5.acc[ii * 3] = body[ii].ax ;	hdf5.acc[ii * 3 + 1] = body[ii].ay;	hdf5.acc[ii * 3 + 2] = body[ii].az ;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
	  hdf5.acc_ext[ii * 3] = body[ii].ax_ext;	  hdf5.acc_ext[ii * 3 + 1] = body[ii].ay_ext;	  hdf5.acc_ext[ii * 3 + 2] = body[ii].az_ext;
	  hdf5.pot_ext[ii] = body[ii].pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
	}/* for(int ii = 0; ii < (int)Ntot; ii++){ */

	int *group_head;      group_head = (int *)malloc(sizeof(int) * kind);
	int *group_num ;      group_num  = (int *)malloc(sizeof(int) * kind);
	if( group_head == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group_head\n");      }
	if( group_num  == NULL ){	__KILL__(stderr, "ERROR: failure to allocate group_num\n" );      }

	for(int ii = 0; ii < kind; ii++){
	  group_head[ii] = (int)group[ii].head;
	  group_num [ii] = (int)group[ii].num;
	}/* for(int ii = 0; ii < kind; ii++){ */

	writeSnapshotMultiGroups(time, steps, &hdf5, file, filenum, hdf5type, kind, group_head, group_num);
	free(group_head);
	free(group_num);
      }/* if( 0 != access(filename, F_OK) ){ */
      else{
	fprintf(stderr, "# \"%s\" was not updated for reducing the elapsed time.\n", filename);
	fflush(stderr);
      }/* else{ */
    }/* if( multi_group ){ */
#endif//USE_HDF5_FORMAT
#endif//DUMP_SPLITTED_SNAPSHOT
#   if  !defined(DUMP_SPLITTED_SNAPSHOT) && !defined(USE_HDF5_FORMAT)
    /* plot spatial distribution */
    plotDistributionMaps(Ntot, kind, group, body, time, xybox, xzbox, zybox, file, ifile, argc, argv);
#endif//!defined(DUMP_SPLITTED_SNAPSHOT) && !defined(USE_HDF5_FORMAT)

    /* plot spherical averaged profile */
    analyzeRadialProfile(Ntot, body, kind, group, &num, &rem, &rad, &rho, &enc,
#ifdef  PLOT_VELOCITY_DISPERSION
			 &sigr,
#endif//PLOT_VELOCITY_DISPERSION
			 (const real)(radius * radius));
    plotRadialProfile(num,  kind, overlay_initial, group, time, rad, rho, enc,
#ifdef  PLOT_VELOCITY_DISPERSION
		      sigr, rvelrange,
#endif//PLOT_VELOCITY_DISPERSION
		      rhorange, encrange, file, ifile, argc, argv);
    rem += num;
    num = 0;

    /* plot horizontal (R := sqrt(x^2 + y^2)) and vertical (|z|) profile */
    analyzeDecomposedProfile(body, kind, group,
			     &num_hor, &rem_hor, &hor_pos, &hor_rho, &hor_zdisp,
			     &num_ver, &rem_ver, &ver_pos, &ver_rho
#ifdef  PLOT_VELOCITY_DISPERSION
			     , &sigR, &sigp, &sigz
#ifdef  PLOT_TOOMRE_Q
			     , &toomre
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
			     );

#ifdef  DUMP_DISKHEIGHT_PROFILE
    for(int ii = skind; ii < kind; ii++){
      char datfile[256];
      sprintf(datfile, "%s/%s.disk%d.snp%.3d.dat", DATAFOLDER, file, ii - skind, filenum);
      FILE *fp;
      fp = fopen(datfile, "wb");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", datfile);      }

      bool success = true;
      size_t tmp;

      tmp = group[ii].prf_hor_num;      if( tmp != fwrite(&(hor_pos  [group[ii].prf_hor_head]), sizeof(real), tmp, fp) )	success = false;
      tmp = group[ii].prf_hor_num;      if( tmp != fwrite(&(hor_zdisp[group[ii].prf_hor_head]), sizeof(real), tmp, fp) )	success = false;

      if( success != true ){	__KILL__(stderr, "ERROR: failure to write \"%s\"\n", datfile);      }
      fclose(fp);
    }/* for(int ii = skind; ii < kind; ii++){ */
#endif//DUMP_DISKHEIGHT_PROFILE

#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
    if( filenum == (start + mpi.rank * interval) ){
      num_hor_t0 = num_hor;
      chkMPIerr(MPI_Bcast(&num_hor_t0, 1, MPI_INT, 0, mpi.comm));
      hor_pos_t0   = (real *)malloc(sizeof(real) * num_hor_t0);
      hor_zdisp_t0 = (real *)malloc(sizeof(real) * num_hor_t0);
      if( (hor_pos_t0 == NULL) || (hor_zdisp_t0 == NULL) ){	__KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");      }
      if( mpi.rank == 0 ){
	for(int ii = 0; ii < num_hor_t0; ii++){
	  hor_pos_t0  [ii] = hor_pos  [ii];
	  hor_zdisp_t0[ii] = hor_zdisp[ii];
	}/* for(int ii = 0; ii < num_hor_t0; ii++){ */
	for(int ii = 0; ii < kind; ii++){
	  prf_hor_head_t0[ii] = group[ii].prf_hor_head;
	  prf_hor_num_t0 [ii] = group[ii].prf_hor_num;
	}/* for(int ii = 0; ii < kind; ii++){ */
      }/* if( mpi.rank == 0 ){ */
      chkMPIerr(MPI_Bcast(  hor_pos_t0, num_hor_t0, MPI_REALDAT, 0, mpi.comm));
      chkMPIerr(MPI_Bcast(hor_zdisp_t0, num_hor_t0, MPI_REALDAT, 0, mpi.comm));
      chkMPIerr(MPI_Bcast(prf_hor_head_t0, kind, MPI_INT, 0, mpi.comm));
      chkMPIerr(MPI_Bcast( prf_hor_num_t0, kind, MPI_INT, 0, mpi.comm));
    }/* if( filenum == (start + mpi.rank * interval) ){ */
#endif//OVERPLOT_INITIAL_DISKHEIGHT
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
    plotHorizontalProfile(num_hor, kind, overlay_initial, group, time, hor_pos, hor_rho, Sigmarange, skind, hor_zdisp, heightrange,
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
			  num_hor_t0, prf_hor_head_t0, prf_hor_num_t0, hor_pos_t0, hor_zdisp_t0,
#endif//OVERPLOT_INITIAL_DISKHEIGHT
#ifdef  PLOT_VELOCITY_DISPERSION
			  sigR, sigp, sigz, Rvelrange,
#ifdef  PLOT_TOOMRE_Q
			  toomre, Qrange,
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
			  file, ifile, argc, argv);
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#ifdef  PLOT_VERTICAL_PROFILE
    plotVerticalProfile(num_ver, kind, overlay_initial, group, time, ver_pos, ver_rho, zrhorange, file, ifile, argc, argv);
#endif//PLOT_VERTICAL_PROFILE
    rem_hor += num_hor;    num_hor = 0;
    rem_ver += num_ver;    num_ver = 0;

    ifile += (int)(interval * mpi.size);
  }/* for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){ */

#ifdef  OUTPUT_BULK_MOTION
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, comPos, nfile * kind * 3, MPI_REALDAT, MPI_SUM, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, comVel, nfile * kind * 3, MPI_REALDAT, MPI_SUM, mpi.comm));
  if( mpi.rank == 0 ){
    fprintf(stdout, "# position of center-of-mass\n");
    fprintf(stdout, "# file");
    for(int ii = 0; ii < kind; ii++)      fprintf(stdout, "\tcom(%d)_x\tcom(%d)_y\tcom(%d)_z", ii, ii, ii);
    fprintf(stdout, "\n");
    for(int ii = 0; ii < nfile; ii++){
      fprintf(stdout, "%d", start + ii * interval);
      for(int kk = 0; kk < kind; kk++)
	fprintf(stdout, "\t%e\t%e\t%e", comPos[INDEX(nfile, kind, 3, ii, kk, 0)], comPos[INDEX(nfile, kind, 3, ii, kk, 1)], comPos[INDEX(nfile, kind, 3, ii, kk, 2)]);
      fprintf(stdout, "\n");
    }/* for(int ii = 0; ii < nfile; ii++){ */

    fprintf(stdout, "# velocity of bulk motion\n");
    fprintf(stdout, "# file");
    for(int ii = 0; ii < kind; ii++)      fprintf(stdout, "\tvel(%d)_x\tvel(%d)_y\tvel(%d)_z", ii, ii, ii);
    fprintf(stdout, "\n");
    for(int ii = 0; ii < nfile; ii++){
      fprintf(stdout, "%d", start + ii * interval);
      for(int kk = 0; kk < kind; kk++)
	fprintf(stdout, "\t%e\t%e\t%e", comVel[INDEX(nfile, kind, 3, ii, kk, 0)], comVel[INDEX(nfile, kind, 3, ii, kk, 1)], comVel[INDEX(nfile, kind, 3, ii, kk, 2)]);
      fprintf(stdout, "\n");
    }/* for(int ii = 0; ii < nfile; ii++){ */
  }/* if( mpi.rank == 0 ){ */
  free(comPos);
  free(comVel);
#endif//OUTPUT_BULK_MOTION

#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
  free(prf_hor_head_t0);
  free(prf_hor_num_t0);
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if( hor_pos_t0 != NULL )    free(hor_pos_t0);
  if( hor_zdisp_t0 != NULL )  free(hor_zdisp_t0);
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#endif//OVERPLOT_INITIAL_DISKHEIGHT


#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
  freeSnapshotArray
    (hdf5_pos, hdf5_vel, hdf5_acc, hdf5_m, hdf5_pot, hdf5_idx
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , hdf5_acc_ext, hdf5_pot_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
#else///USE_HDF5_FORMAT
  freeParticleData
    (idx, pos, acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
     vel, ti
#else///BLOCK_TIME_STEP
     vx, vy, vz
#endif//BLOCK_TIME_STEP
     );
#endif//USE_HDF5_FORMAT
  free(body);
  free(rad);  free(rho);  free(enc);
  free(hor_pos);  free(hor_rho);  free(hor_zdisp);
  free(ver_pos);  free(ver_rho);

  if( overlay_initial )
    for(int ii = 0; ii < kind; ii++){
      free(group[ii].rad);
      free(group[ii].rho);
      free(group[ii].enc);
      free(group[ii].Sig);
#ifndef USE_HDF5_FORMAT
      free(group[ii].psi);
#endif//USE_HDF5_FORMAT

      free(group[ii].hor);
      free(group[ii].ver);
      free(group[ii].Sigma);
      free(group[ii].Menc);

      for(int jj = 0; jj < NUM_HORIZONTAL_BIN; jj++){
	free(group[ii].ver_pos[jj]);
	free(group[ii].ver_rho[jj]);
      }/* for(int jj = 0; jj < NUM_HORIZONTAL_BIN; jj++){ */

      if( ii >= skind ){
	free(group[ii].disk_radius );
	free(group[ii].disk_height );
	free(group[ii].disk_rho    );
#ifndef USE_HDF5_FORMAT
	free(group[ii].disk_pot    );
	free(group[ii].disk_sig    );
#endif//USE_HDF5_FORMAT
	free(group[ii].disk_Sigma  );
	free(group[ii].disk_rad1d  );
	free(group[ii].disk_vol    );
      }/* if( ii >= skind ){ */
    }/* for(int ii = 0; ii < kind; ii++){ */
  free(group);


  exitMPI();

  return (0);
}


void plotDistributionMaps
(ulong Np, const int ngroup, model *group, nbody_particle *body, PLFLT time,
 PLplotPltRange xybox, PLplotPltRange xzbox, PLplotPltRange zybox,
 char *file, const int filenum, int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /* global setting(s) */
  /* maximum number of data kinds */
  const PLINT pkind = (PLINT)ngroup;
  const PLINT lkind = 0;

  /* preparation for dataset */
  /* memory allocation for index */
  PLINT * num;  allocPLINT(& num, pkind);
  /* set number of data points */
#ifdef  REDUCE_PARTICLE_DISTRIBUTION_MAP
  const PLINT Ntot = (Np > REDUCE_PARTICLE_DISTRIBUTION_MAP) ? REDUCE_PARTICLE_DISTRIBUTION_MAP : Np;
  const PLINT Ninc = (PLINT)Np / Ntot;
#else///REDUCE_PARTICLE_DISTRIBUTION_MAP
  const PLINT Ntot = (PLINT)Np;
  const PLINT Ninc = 1;
#endif//REDUCE_PARTICLE_DISTRIBUTION_MAP
  for(PLINT ii = 0; ii < pkind; ii++)    num[ii] = (PLINT)group[ii].num / Ninc;
  PLINT Ntmp = num[0];
  for(PLINT ii = 1; ii < pkind; ii++)    Ntmp += num[ii];
  if( Ntmp != Ntot )
    num[0] += Ntot - Ntmp;

  /* memory allocation for data */
  PLFLT *_xpt, *_ypt, *_zpt;
  {
    PLINT tot = 0;
    for(PLINT ii = 0; ii < pkind; ii++)      tot += num[ii];
    allocPLFLT(&_xpt, tot);
    allocPLFLT(&_ypt, tot);
    allocPLFLT(&_zpt, tot);
  }
  PLFLT **xpt;  allocPointer4PLFLT(&xpt, pkind);
  PLFLT **ypt;  allocPointer4PLFLT(&ypt, pkind);
  PLFLT **zpt;  allocPointer4PLFLT(&zpt, pkind);
  xpt[0] = _xpt;
  ypt[0] = _ypt;
  zpt[0] = _zpt;
  for(PLINT ii = 1; ii < pkind; ii++){
    xpt[ii] = xpt[ii - 1] + num[ii - 1];
    ypt[ii] = ypt[ii - 1] + num[ii - 1];
    zpt[ii] = zpt[ii - 1] + num[ii - 1];
  }

  /* data preparation */
  for(PLINT jj = 0; jj < num[0]; jj++){
    /* const PLINT kk = jj * Ninc; */
    const PLFLT xpos = (PLFLT)body[jj].x * (PLFLT)length2astro;
    const PLFLT ypos = (PLFLT)body[jj].y * (PLFLT)length2astro;
    const PLFLT zpos = (PLFLT)body[jj].z * (PLFLT)length2astro;
    xpt[0][jj] = xpos;
    ypt[0][jj] = ypos;
    zpt[0][jj] = zpos;
  }
  for(PLINT kk = pkind - 1; kk >= 1; kk--){
    for(PLINT jj = 0; jj < num[kk]; jj++){
      const PLINT ll = (PLINT)group[kk].head + jj;

      const PLFLT xpos = (PLFLT)body[ll].x * (PLFLT)length2astro;
      const PLFLT ypos = (PLFLT)body[ll].y * (PLFLT)length2astro;
      const PLFLT zpos = (PLFLT)body[ll].z * (PLFLT)length2astro;
      xpt[kk][jj] = xpos;
      ypt[kk][jj] = ypos;
      zpt[kk][jj] = zpos;
    }
  }

  /* set symbol style */
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  for(PLINT ii = 0; ii < pkind; ii++){
    /* pt[ii].type  = (char *)PLplotSymbolType[smallDot]; */
    sprintf(pt[ii].type, "%s", PLplotSymbolType[smallDot]);
    pt[ii].scale =       PLplotSymbolSize[smallDot];
  }
  if( pkind - 1 >= 1 )    pt[pkind - 1].color = RED;
  if( pkind - 2 >= 1 )    pt[pkind - 2].color = BLUE;
  if( pkind - 3 >= 1 )    pt[pkind - 3].color = MAGENTA;

  /* set labels */
  char xlab[PLplotCharWords];  sprintf(xlab, "#fix #fr(%s)", length_astro_unit_name4plot);
  char ylab[PLplotCharWords];  sprintf(ylab, "#fiy #fr(%s)", length_astro_unit_name4plot);
  char zlab[PLplotCharWords];  sprintf(zlab, "#fiz #fr(%s)", length_astro_unit_name4plot);

  /* set caption(s) */
  PLplotCaption xycap;  setDefaultCaption(&xycap);  xycap.write = false;
  PLplotCaption xzcap;  setDefaultCaption(&xzcap);  xzcap.write = false;
  PLplotCaption zycap;  setDefaultCaption(&zycap);  zycap.write = false;
  sprintf(xycap.side, "t");
  sprintf(xzcap.side, "t");
  sprintf(zycap.side, "t");
  sprintf(xycap.text, "#fit#fr = %e (%s)", (double)time * time2astro, time_astro_unit_name4plot);
  sprintf(xzcap.text, "#fit#fr = %e (%s)", (double)time * time2astro, time_astro_unit_name4plot);
  sprintf(zycap.text, "#fit#fr = %e (%s)", (double)time * time2astro, time_astro_unit_name4plot);

  /* set legends */
  PLplotLegend xyleg;  setDefaultLegend(&xyleg, false);  xyleg.write = false;
  PLplotLegend xzleg;  setDefaultLegend(&xzleg, false);  xzleg.write = false;
  PLplotLegend zyleg;  setDefaultLegend(&zyleg, false);  zyleg.write = false;
  char *xylegTxt;
  {
    allocChar4PLplot(&xylegTxt, pkind);
    allocPointer4Char4PLplot(&(xyleg.text), pkind);
    assignChar4PLplot(pkind, xyleg.text, xylegTxt);
  }
  sprintf(xyleg.text[0], "position");
  char *xzlegTxt;
  {
    allocChar4PLplot(&xzlegTxt, pkind);
    allocPointer4Char4PLplot(&(xzleg.text), pkind);
    assignChar4PLplot(pkind, xzleg.text, xzlegTxt);
  }
  sprintf(xzleg.text[0], "position");
  char *zylegTxt;
  {
    allocChar4PLplot(&zylegTxt, pkind);
    allocPointer4Char4PLplot(&(zyleg.text), pkind);
    assignChar4PLplot(pkind, zyleg.text, zylegTxt);
  }
  sprintf(zyleg.text[0], "position");
  PLBOOL xyuni = false;
  PLBOOL xzuni = false;
  PLBOOL zyuni = false;

  char title[PLplotCharWords];
  sprintf(title, "#fit#fr = %6.2f (%s)", (PLFLT)((double)time * time2astro), time_astro_unit_name4plot);



  /* memory allocation to exploit multiple-plot function */
  /* maximum number of panels in each direction */
  const PLINT nxpanel = 2;
  const PLINT nypanel = 2;

  /* specify plot or skip panel */
  PLBOOL *plot;
  allocPLBOOL(&plot, nxpanel * nypanel);
  /* arrays related to draw line(s) and point(s) */
  PLINT *nlkind;  allocPLINT(&nlkind, nxpanel * nypanel);
  PLINT *npkind;  allocPLINT(&npkind, nxpanel * nypanel);
  PLINT **lnum;  allocPointer4PLINT(&lnum, nxpanel * nypanel);
  PLINT **pnum;  allocPointer4PLINT(&pnum, nxpanel * nypanel);
  PLFLT ***lx;  allocDoublePointer4PLFLT(&lx, nxpanel * nypanel);
  PLFLT ***ly;  allocDoublePointer4PLFLT(&ly, nxpanel * nypanel);
  PLFLT ***px;  allocDoublePointer4PLFLT(&px, nxpanel * nypanel);
  PLFLT ***py;  allocDoublePointer4PLFLT(&py, nxpanel * nypanel);
  PLplotLineStyle ** line;  allocPointer4PLplotLineStyle(& line, nxpanel * nypanel);
  PLplotPointType **point;  allocPointer4PLplotPointType(&point, nxpanel * nypanel);
  /* arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /* arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /* arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /* arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /* arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /* array to set figure name */
  char figfile[PLplotCharWords];
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel * nypanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel * nypanel);
  assignChar4PLplot(nxpanel * nypanel, figname, _figname);


  /* configure to time evolution of energy */
  /* common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

      /* global setting(s) */
      plot[idx] = true;

      /* line setting(s) */
      nlkind[idx] = lkind;
      lnum  [idx] = NULL;
      line  [idx] = NULL;
      lx    [idx] = NULL;
      ly    [idx] = NULL;

      /* point setting(s) */
      npkind[idx] = pkind;
      pnum  [idx] = num;
      point [idx] = pt;

      /* errorbar setting(s) */
      errbar[idx] = 0;
    }
  }

  /* a NULL panel */
  {
    const PLINT idx = INDEX2D(nxpanel, nypanel, 1, 1);
    plot[idx] = false;
  }

  /* setting(s) for position in x-direction shown along horizontal axis */
  for(PLINT jj = 0; jj < nypanel; jj++){
    const PLINT idx = INDEX2D(nxpanel, nypanel, 0, jj);
    /* point setting(s) */
    px[idx] = xpt;
    /* label setting(s) */
    xlabel[idx] = xlab;
  }

  /* setting(s) for position in y-direction shown along vertical axis */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 0);

    /* point setting(s) */
    py[idx] = ypt;
    /* label setting(s) */
    ylabel[idx] = ylab;
  }

  /* individual setting(s) */
  /* data arrays and labes */
  py[INDEX2D(nxpanel, nypanel, 0, 1)] = zpt;  ylabel[INDEX2D(nxpanel, nypanel, 0, 1)] = zlab;
  px[INDEX2D(nxpanel, nypanel, 1, 0)] = zpt;  xlabel[INDEX2D(nxpanel, nypanel, 1, 0)] = zlab;
  /* plot area */
  range[INDEX2D(nxpanel, nypanel, 0, 0)] = xybox;
  range[INDEX2D(nxpanel, nypanel, 0, 1)] = xzbox;
  range[INDEX2D(nxpanel, nypanel, 1, 0)] = zybox;
  /* captions */
  cap[INDEX2D(nxpanel, nypanel, 0, 0)] = xycap;
  cap[INDEX2D(nxpanel, nypanel, 0, 1)] = xzcap;
  cap[INDEX2D(nxpanel, nypanel, 1, 0)] = zycap;
  /* legends */
  leg[INDEX2D(nxpanel, nypanel, 0, 0)] = xyleg;  uni[INDEX2D(nxpanel, nypanel, 0, 0)] = xyuni;
  leg[INDEX2D(nxpanel, nypanel, 0, 1)] = xzleg;  uni[INDEX2D(nxpanel, nypanel, 0, 1)] = xzuni;
  leg[INDEX2D(nxpanel, nypanel, 1, 0)] = zyleg;  uni[INDEX2D(nxpanel, nypanel, 1, 0)] = zyuni;
  /* individual file names */
  sprintf(figfile,                                  "%s_%s_%.3d", file, "posmap", filenum);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 0)], "%s_%s_%.3d", file, "xyproj", filenum);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 1)], "%s_%s_%.3d", file, "xzproj", filenum);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 1, 0)], "%s_%s_%.3d", file, "zyproj", filenum);


  /* create figure(s) */
  plotData(nxpanel, nypanel, plot, true, true,
  	   nlkind,  line, lnum, lx, ly,
  	   npkind, point, pnum, px, py,
  	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
  	   cap, leg, uni, range,
  	   xlabel, ylabel, "", figfile, argc, argv);
#ifdef  PLOT_SPLITTED_DISTRIBUTION_MAPS
  for(PLINT idx = 0; idx < nxpanel * nypanel; idx++)
    if( plot[idx] == true )
      plotData(1, 1, &plot[idx], false, false,
  	       &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
  	       &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
  	       &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
  	       &cap[idx], &leg[idx], &uni[idx], &range[idx],
  	       &xlabel[idx], &ylabel[idx], "", figname[idx], argc, argv);
#endif//PLOT_SPLITTED_DISTRIBUTION_MAPS


  /* memory deallocation */
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  free(figname);  free(_figname);
  free(_xpt);  free(xpt);
  free(_ypt);  free(ypt);
  free(_zpt);  free(zpt);
  free(pt);
  free(xylegTxt);
  free(xzlegTxt);
  free(zylegTxt);
  free(num);


  __NOTE__("%s\n", "end");
}


void plotRadialProfile
(const int ntot, const int ngroup, const bool overlay_initial, model *group, PLFLT time,
 real *rad, real *rho, real *enc,
#ifdef  PLOT_VELOCITY_DISPERSION
 real *sig, PLplotPltRange velbox,
#endif//PLOT_VELOCITY_DISPERSION
 PLplotPltRange rhobox, PLplotPltRange encbox,
 char *file, const int filenum, int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /* global setting(s) */
  int skip = 0;
  for(int ii = 0; ii < ngroup; ii++)
    if( group[ii].num == 1 )
      skip++;

  /* maximum number of data kinds */
  const PLINT  pkind = ngroup - skip;
  const PLINT  lkind = overlay_initial ? (ngroup - skip) : 0;
  const PLBOOL unify = false;

  PLINT *grpIdx;  allocPLINT(&grpIdx, pkind);
  {
    int jj = 0;
    for(int ii = 0; ii < ngroup; ii++)
      if( group[ii].num > 1 ){
	grpIdx[jj] = ii;
	jj++;
      }/* if( group[ii].num > 1 ){ */
  }


  /* preparation for dataset */
  /* memory allocation for index */
  PLINT *npt;  allocPLINT(&npt, pkind);
  PLINT *nln;  allocPLINT(&nln, lkind);

  /* set number of data points */
  for(PLINT ii = 0; ii < pkind; ii++)    npt[ii] = (PLINT)group[grpIdx[ii]].prf_num;
  for(PLINT ii = 0; ii < lkind; ii++)    nln[ii] = (PLINT)NUM_ANALYTIC;

  /* memory allocation for data */
  PLFLT *_rad_p;  allocPLFLT(&_rad_p,                 ntot);  PLFLT **rad_p;  allocPointer4PLFLT(&rad_p, pkind);
  PLFLT *_rho_p;  allocPLFLT(&_rho_p,                 ntot);  PLFLT **rho_p;  allocPointer4PLFLT(&rho_p, pkind);
  PLFLT *_enc_p;  allocPLFLT(&_enc_p,                 ntot);  PLFLT **enc_p;  allocPointer4PLFLT(&enc_p, pkind);
  PLFLT *_rad_l;  allocPLFLT(&_rad_l, NUM_ANALYTIC * lkind);  PLFLT **rad_l;  allocPointer4PLFLT(&rad_l, lkind);
  PLFLT *_rho_l;  allocPLFLT(&_rho_l, NUM_ANALYTIC * lkind);  PLFLT **rho_l;  allocPointer4PLFLT(&rho_l, lkind);
  PLFLT *_enc_l;  allocPLFLT(&_enc_l, NUM_ANALYTIC * lkind);  PLFLT **enc_l;  allocPointer4PLFLT(&enc_l, lkind);
#ifdef  PLOT_VELOCITY_DISPERSION
  PLFLT *_sig_p;  allocPLFLT(&_sig_p,                 ntot);  PLFLT **sig_p;  allocPointer4PLFLT(&sig_p, pkind);
  PLFLT *_sig_l;  allocPLFLT(&_sig_l, NUM_ANALYTIC * lkind);  PLFLT **sig_l;  allocPointer4PLFLT(&sig_l, lkind);
#endif//PLOT_VELOCITY_DISPERSION

  /* set pointers */
  rad_p[0] = _rad_p;
  rho_p[0] = _rho_p;
  enc_p[0] = _enc_p;
#ifdef  PLOT_VELOCITY_DISPERSION
  sig_p[0] = _sig_p;
#endif//PLOT_VELOCITY_DISPERSION
  for(PLINT ii = 1; ii < pkind; ii++){
    rad_p[ii] = rad_p[ii - 1] + npt[ii - 1];
    rho_p[ii] = rho_p[ii - 1] + npt[ii - 1];
    enc_p[ii] = enc_p[ii - 1] + npt[ii - 1];
#ifdef  PLOT_VELOCITY_DISPERSION
    sig_p[ii] = sig_p[ii - 1] + npt[ii - 1];
#endif//PLOT_VELOCITY_DISPERSION
  }/* for(PLINT ii = 1; ii < pkind; ii++){ */
  for(PLINT ii = 0; ii < lkind; ii++){
    rad_l[ii] = _rad_l + ii * NUM_ANALYTIC;
    rho_l[ii] = _rho_l + ii * NUM_ANALYTIC;
    enc_l[ii] = _enc_l + ii * NUM_ANALYTIC;
#ifdef  PLOT_VELOCITY_DISPERSION
    sig_l[ii] = _sig_l + ii * NUM_ANALYTIC;
#endif//PLOT_VELOCITY_DISPERSION
  }/* for(PLINT ii = 0; ii < lkind; ii++){ */

  /* data preparation for particle distribution */
  for(PLINT ii = 0; ii < pkind; ii++)
    for(PLINT jj = 0; jj < npt[ii]; jj++){
      const int kk = group[grpIdx[ii]].prf_head + jj;

      rad_p[ii][jj] = (PLFLT)log10((double)rad[kk] *  length2astro);
      rho_p[ii][jj] = (PLFLT)log10((double)rho[kk] * density2astro);
      enc_p[ii][jj] = (PLFLT)log10((double)enc[kk] *    mass2astro);
#ifdef  PLOT_VELOCITY_DISPERSION
      sig_p[ii][jj] = (PLFLT)sig[kk] * (PLFLT)velocity2astro;
#endif//PLOT_VELOCITY_DISPERSION
    }/* for(PLINT jj = 0; jj < npt[ii]; jj++){ */

  /* data preparation for analytic distribution */
  for(PLINT ii = 0; ii < lkind; ii++){
    const PLINT jskip = group[grpIdx[ii]].nrad / nln[ii];
    for(PLINT jj = 0; jj < nln[ii]; jj++){
      const int kk = jj * jskip;

      rad_l[ii][jj] = (PLFLT)log10((double)group[grpIdx[ii]].rad[kk]);
      rho_l[ii][jj] = (PLFLT)log10((double)group[grpIdx[ii]].rho[kk]);
      enc_l[ii][jj] = (PLFLT)log10((double)group[grpIdx[ii]].enc[kk]);
#ifdef  PLOT_VELOCITY_DISPERSION
      sig_l[ii][jj] = (PLFLT)group[grpIdx[ii]].sigr[kk];
#endif//PLOT_VELOCITY_DISPERSION
    }/* for(PLINT jj = 0; jj < nln[ii]; jj++){ */
  }/* for(PLINT ii = 0; ii < lkind; ii++){ */

  /* set symbol style */
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);

  /* set labels */
  char radlab[PLplotCharWords];  sprintf(radlab, "#fir #fr(%s)", length_astro_unit_name4plot);
  char rholab[PLplotCharWords];  sprintf(rholab, "#fi#gr #fr(%s)", density_astro_unit_name4plot);
  char enclab[PLplotCharWords];  sprintf(enclab, "#fiM#fr#denc#u (%s)", mass_astro_unit_name4plot);
#ifdef  PLOT_VELOCITY_DISPERSION
  char siglab[PLplotCharWords];  sprintf(enclab, "#fi#gs#dr#u #fr(%s)", velocity_astro_unit_name4plot);
#endif//PLOT_VELOCITY_DISPERSION

  /* set caption(s) */
  PLplotCaption rhocap;  setDefaultCaption(&rhocap);  rhocap.write = false;
  PLplotCaption enccap;  setDefaultCaption(&enccap);  enccap.write = false;
  sprintf(rhocap.side, "%s", "t");
  sprintf(enccap.side, "%s", "t");
  sprintf(rhocap.text, "%s = %e (%s)", "#fit#fr", (double)time * time2astro, time_astro_unit_name4plot);
  sprintf(enccap.text, "%s = %e (%s)", "#fit#fr", (double)time * time2astro, time_astro_unit_name4plot);
#ifdef  PLOT_VELOCITY_DISPERSION
  PLplotCaption sigcap;  setDefaultCaption(&sigcap);  sigcap.write = false;
  sprintf(sigcap.side, "%s", "t");
  sprintf(sigcap.text, "%s = %e (%s)", "#fit#fr", (double)time * time2astro, time_astro_unit_name4plot);
#endif//PLOT_VELOCITY_DISPERSION

  /* set legends */
  PLplotLegend rholeg;  setDefaultLegend(&rholeg, false);  rholeg.write = false;
  PLplotLegend encleg;  setDefaultLegend(&encleg, false);  encleg.write = false;
  char *rholegTxt;
  char *enclegTxt;
#ifdef  PLOT_VELOCITY_DISPERSION
  PLplotLegend sigleg;  setDefaultLegend(&sigleg, false);  sigleg.write = false;
  char *siglegTxt;
#endif//PLOT_VELOCITY_DISPERSION
  {
    PLINT dkind = (unify) ? (IMAX(pkind, lkind)) : (pkind + lkind);
    allocChar4PLplot(&rholegTxt, dkind);    allocPointer4Char4PLplot(&(rholeg.text), dkind);    assignChar4PLplot(dkind, rholeg.text, rholegTxt);
    allocChar4PLplot(&enclegTxt, dkind);    allocPointer4Char4PLplot(&(encleg.text), dkind);    assignChar4PLplot(dkind, encleg.text, enclegTxt);
#ifdef  PLOT_VELOCITY_DISPERSION
    allocChar4PLplot(&siglegTxt, dkind);    allocPointer4Char4PLplot(&(sigleg.text), dkind);    assignChar4PLplot(dkind, sigleg.text, siglegTxt);
#endif//PLOT_VELOCITY_DISPERSION
  }

  char title[PLplotCharWords];
  sprintf(title, "#fit#fr = %6.2f (%s)", (PLFLT)((double)time * time2astro), time_astro_unit_name4plot);


  /* memory allocation to exploit multiple-plot function */
  /* maximum number of panels in each direction */
  const PLINT nxpanel = 1;
#ifdef  PLOT_VELOCITY_DISPERSION
  const PLINT nypanel = 3;
#else///PLOT_VELOCITY_DISPERSION
  const PLINT nypanel = 2;
#endif//PLOT_VELOCITY_DISPERSION

  /* specify plot or skip panel */
  PLBOOL *plot;
  allocPLBOOL(&plot, nxpanel * nypanel);
  /* arrays related to draw line(s) and point(s) */
  PLINT *nlkind;  allocPLINT(&nlkind, nxpanel * nypanel);
  PLINT *npkind;  allocPLINT(&npkind, nxpanel * nypanel);
  PLINT **lnum;  allocPointer4PLINT(&lnum, nxpanel * nypanel);
  PLINT **pnum;  allocPointer4PLINT(&pnum, nxpanel * nypanel);
  PLFLT ***lx;  allocDoublePointer4PLFLT(&lx, nxpanel * nypanel);
  PLFLT ***ly;  allocDoublePointer4PLFLT(&ly, nxpanel * nypanel);
  PLFLT ***px;  allocDoublePointer4PLFLT(&px, nxpanel * nypanel);
  PLFLT ***py;  allocDoublePointer4PLFLT(&py, nxpanel * nypanel);
  PLplotLineStyle ** line;  allocPointer4PLplotLineStyle(& line, nxpanel * nypanel);
  PLplotPointType **point;  allocPointer4PLplotPointType(&point, nxpanel * nypanel);
  /* arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /* arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /* arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /* arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /* arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /* array to set figure name */
  char figfile[PLplotCharWords];
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel * nypanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel * nypanel);
  assignChar4PLplot(nxpanel * nypanel, figname, _figname);


  /* configuration about radial profile of particle distribution  */
  /* common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

      /* global setting(s) */
      plot[idx] = true;
      /* line setting(s) */
      nlkind[idx] = lkind;
      lnum  [idx] = nln;
      line  [idx] = ls;
      lx    [idx] = rad_l;
      /* point setting(s) */
      npkind[idx] = pkind;
      pnum  [idx] = npt;
      point [idx] = pt;
      px    [idx] = rad_p;
      /* errorbar setting(s) */
      errbar[idx] = 0;
      /* legend(s) */
      uni[idx] = unify;
      /* label setting(s) */
      xlabel[idx] = radlab;
    }
  }

  /* setting(s) for enclosed mass profile */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 0);

    /* line setting(s) */
    ly    [idx] = enc_l;
    /* point setting(s) */
    py    [idx] = enc_p;
    /* plot area */
    range[idx] = encbox;
    /* label setting(s) */
    ylabel[idx] = enclab;
    /* miscs */
    cap[idx] = enccap;
    leg[idx] = encleg;
  }

  /* setting(s) for volume density profile */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 1);

    /* line setting(s) */
    ly    [idx] = rho_l;
    /* point setting(s) */
    py    [idx] = rho_p;
    /* plot area */
    range[idx] = rhobox;
    /* label setting(s) */
    ylabel[idx] = rholab;
    /* miscs */
    cap[idx] = rhocap;
    leg[idx] = rholeg;
  }

#ifdef  PLOT_VELOCITY_DISPERSION
  /* setting(s) for radial velocity dispersion profile */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 2);

    /* line setting(s) */
    ly    [idx] = sig_l;
    /* point setting(s) */
    py    [idx] = sig_p;
    /* plot area */
    range[idx] = velbox;
    /* label setting(s) */
    ylabel[idx] = siglab;
    /* miscs */
    cap[idx] = sigcap;
    leg[idx] = sigleg;
  }
#endif//PLOT_VELOCITY_DISPERSION

  /* individual file names */
  sprintf(figfile                                 , "%s_%s_%.3d", file, "profile", filenum);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 0)], "%s_%s_%.3d", file, "encmass", filenum);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 1)], "%s_%s_%.3d", file, "density", filenum);
#ifdef  PLOT_VELOCITY_DISPERSION
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 2)], "%s_%s_%.3d", file, "sigma_r", filenum);
#endif//PLOT_VELOCITY_DISPERSION


  /* create figure(s) */
#ifdef  PLOT_INDIVIDUAL_PROFILES
  for(PLINT idx = 0; idx < nxpanel * nypanel; idx++)
    plotData(1, 1, &plot[idx], false, false,
  	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
  	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
  	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
  	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
  	     &xlabel[idx], &ylabel[idx], "", figname[idx], argc, argv);
#endif//PLOT_INDIVIDUAL_PROFILES
  plotData(nxpanel, nypanel, plot, true, true,
	   nlkind,  line, lnum, lx, ly,
	   npkind, point, pnum, px, py,
	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
	   cap, leg, uni, range,
	   xlabel, ylabel, "", figfile, argc, argv);


  /* memory deallocation */
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  free(_rad_p);  free(rad_p);  free(_rad_l);  free(rad_l);
  free(_rho_p);  free(rho_p);  free(_rho_l);  free(rho_l);
  free(_enc_p);  free(enc_p);  free(_enc_l);  free(enc_l);
#ifdef  PLOT_VELOCITY_DISPERSION
  free(_sig_p);  free(sig_p);  free(_sig_l);  free(sig_l);
#endif//PLOT_VELOCITY_DISPERSION
  free(pt);  free(ls);
  free(rholegTxt);
  free(enclegTxt);
#ifdef  PLOT_VELOCITY_DISPERSION
  free(siglegTxt);
#endif//PLOT_VELOCITY_DISPERSION
  free(grpIdx);


  __NOTE__("%s\n", "end");
}


void plotHorizontalProfile
(const int ntot, const int ngroup, const bool overlay_initial, model *group, PLFLT time,
 real *rad, real *rho, PLplotPltRange rhobox,
 const int nspheroids, real *zdisp, PLplotPltRange zbox,
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
 const int ntot_t0, int * restrict hor_head_t0, int * restrict hor_num_t0, real * restrict rad_t0, real * restrict zdisp_t0,
#endif//OVERPLOT_INITIAL_DISKHEIGHT
#ifdef  PLOT_VELOCITY_DISPERSION
 real * restrict sigR, real * restrict sigp, real * restrict sigz, PLplotPltRange Rvelrange,
#ifdef  PLOT_TOOMRE_Q
 real * restrict toomre, PLplotPltRange Qrange,
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
 char *file, const int filenum, int argc, char **argv)
{
  __NOTE__("%s\n", "start");


#ifdef  PLOT_VELOCITY_DISPERSION
#ifdef  PLOT_TOOMRE_Q
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION


  /* global setting(s) */
  int skip = 0;
  for(int ii = 0; ii < ngroup; ii++)
    if( group[ii].num == 1 )
      skip++;

  /* maximum number of data kinds */
  const PLINT  pkind = ngroup - skip;
  const PLINT  lkind = overlay_initial ? (ngroup - skip) : 0;
  const PLBOOL unify = false;

  const PLINT  ndisk = ngroup - nspheroids;
  const bool    disk = (ndisk != 0) ? true : false;

  PLINT *grpIdx;  allocPLINT(&grpIdx, pkind);
  {
    int jj = 0;
    for(int ii = 0; ii < ngroup; ii++)
      if( group[ii].num > 1 ){
	grpIdx[jj] = ii;
	jj++;
      }/* if( group[ii].num > 1 ){ */
  }


  /* preparation for dataset */
  /* memory allocation for index */
  PLINT *npt;  allocPLINT(&npt, pkind);
  PLINT *nln;  allocPLINT(&nln, lkind);

  /* set number of data points */
  for(PLINT ii = 0; ii < pkind; ii++)    npt[ii] = (PLINT)group[grpIdx[ii]].prf_hor_num;
  for(PLINT ii = 0; ii < lkind; ii++)    nln[ii] = (PLINT)NUM_ANALYTIC;

  /* memory allocation for data */
  PLFLT *_rad_p;  allocPLFLT(&_rad_p,                 ntot);  PLFLT **rad_p;  allocPointer4PLFLT(&rad_p, pkind);
  PLFLT *_rho_p;  allocPLFLT(&_rho_p,                 ntot);  PLFLT **rho_p;  allocPointer4PLFLT(&rho_p, pkind);
  PLFLT *_rad_l;  allocPLFLT(&_rad_l, NUM_ANALYTIC * lkind);  PLFLT **rad_l;  allocPointer4PLFLT(&rad_l, lkind);
  PLFLT *_rho_l;  allocPLFLT(&_rho_l, NUM_ANALYTIC * lkind);  PLFLT **rho_l;  allocPointer4PLFLT(&rho_l, lkind);
  PLFLT *_hor_p;  allocPLFLT(&_hor_p,                 ntot);  PLFLT **hor_p;  allocPointer4PLFLT(&hor_p, ndisk);
  PLFLT *_ver_p;  allocPLFLT(&_ver_p,                 ntot);  PLFLT **ver_p;  allocPointer4PLFLT(&ver_p, ndisk);
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
  PLINT *nln_t0;  allocPLINT(&nln_t0, ndisk * 2);
  PLFLT *_hor_l;  allocPLFLT(&_hor_l,       ntot_t0 + ntot);  PLFLT **hor_l;  allocPointer4PLFLT(&hor_l, ndisk * 2);
  PLFLT *_ver_l;  allocPLFLT(&_ver_l,       ntot_t0 + ntot);  PLFLT **ver_l;  allocPointer4PLFLT(&ver_l, ndisk * 2);
  for(int ii = 0; ii < ndisk; ii++)    nln_t0[        ii] = group[grpIdx[nspheroids - skip + ii]].prf_hor_num;
  for(int ii = 0; ii < ndisk; ii++)    nln_t0[ndisk + ii] = hor_num_t0[nspheroids - skip + ii];
#endif//OVERPLOT_INITIAL_DISKHEIGHT
#ifdef  PLOT_VELOCITY_DISPERSION
  PLFLT *_sigR_p;  allocPLFLT(&_sigR_p,                 ntot);  PLFLT **sigR_p;  allocPointer4PLFLT(&sigR_p, ndisk);
  PLFLT *_sigR_l;  allocPLFLT(&_sigR_l, NUM_ANALYTIC * ndisk);  PLFLT **sigR_l;  allocPointer4PLFLT(&sigR_l, ndisk);
  PLFLT *_sigp_p;  allocPLFLT(&_sigp_p,                 ntot);  PLFLT **sigp_p;  allocPointer4PLFLT(&sigp_p, ndisk);
  PLFLT *_sigp_l;  allocPLFLT(&_sigp_l, NUM_ANALYTIC * ndisk);  PLFLT **sigp_l;  allocPointer4PLFLT(&sigp_l, ndisk);
  PLFLT *_sigz_p;  allocPLFLT(&_sigz_p,                 ntot);  PLFLT **sigz_p;  allocPointer4PLFLT(&sigz_p, ndisk);
  PLFLT *_sigz_l;  allocPLFLT(&_sigz_l, NUM_ANALYTIC * ndisk);  PLFLT **sigz_l;  allocPointer4PLFLT(&sigz_l, ndisk);
#ifdef  PLOT_TOOMRE_Q
  PLFLT *_toomre_p;  allocPLFLT(&_toomre_p,                 ntot);  PLFLT **toomre_p;  allocPointer4PLFLT(&toomre_p, ndisk);
  PLFLT *_toomre_l;  allocPLFLT(&_toomre_l, NUM_ANALYTIC * ndisk);  PLFLT **toomre_l;  allocPointer4PLFLT(&toomre_l, ndisk);
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION

  /* set pointers */
  rad_p[0] = _rad_p;  rho_p[0] = _rho_p;
  hor_p[0] = _hor_p;  ver_p[0] = _ver_p;
  for(PLINT ii = 1; ii < pkind; ii++){
    rad_p[ii] = rad_p[ii - 1] + npt[ii - 1];
    rho_p[ii] = rho_p[ii - 1] + npt[ii - 1];
  }/* for(PLINT ii = 1; ii < pkind; ii++){ */
#ifdef  PLOT_VELOCITY_DISPERSION
  sigR_p[0] = _sigR_p;
  sigp_p[0] = _sigp_p;
  sigz_p[0] = _sigz_p;
#ifdef  PLOT_TOOMRE_Q
  toomre_p[0] = _toomre_p;
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
  for(PLINT ii = 1; ii < ndisk; ii++){
    /* hor_p[ii] = hor_p[ii - 1] + npt[ii - 1]; */
    /* ver_p[ii] = ver_p[ii - 1] + npt[ii - 1]; */
    hor_p[ii] = hor_p[ii - 1] + npt[nspheroids - skip + ii - 1];
    ver_p[ii] = ver_p[ii - 1] + npt[nspheroids - skip + ii - 1];
#ifdef  PLOT_VELOCITY_DISPERSION
    sigR_p[ii] = sigR_p[ii - 1] + npt[nspheroids - skip + ii - 1];
    sigp_p[ii] = sigp_p[ii - 1] + npt[nspheroids - skip + ii - 1];
    sigz_p[ii] = sigz_p[ii - 1] + npt[nspheroids - skip + ii - 1];
#ifdef  PLOT_TOOMRE_Q
    toomre_p[ii] = toomre_p[ii - 1] + npt[nspheroids - skip + ii - 1];
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
  }/* for(PLINT ii = 1; ii < ndisk; ii++){ */

#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
  hor_l[0] = _hor_l;
  ver_l[0] = _ver_l;
  for(PLINT ii = 1; ii < ndisk * 2; ii++){
    hor_l[ii] = hor_l[ii - 1] + nln_t0[ii - 1];
    ver_l[ii] = ver_l[ii - 1] + nln_t0[ii - 1];
  }/* for(PLINT ii = 1; ii < ndisk * 2; ii++){ */
#endif//OVERPLOT_INITIAL_DISKHEIGHT
#ifdef  PLOT_VELOCITY_DISPERSION
  for(PLINT ii = 0; ii < ndisk; ii++){
    sigR_l[ii] = _sigR_l + ii * NUM_ANALYTIC;
    sigp_l[ii] = _sigp_l + ii * NUM_ANALYTIC;
    sigz_l[ii] = _sigz_l + ii * NUM_ANALYTIC;
#ifdef  PLOT_TOOMRE_Q
    toomre_l[ii] = _toomre_l + ii * NUM_ANALYTIC;
#endif//PLOT_TOOMRE_Q
  }/* for(PLINT ii = 1; ii < ndisk; ii++){ */
#endif//PLOT_VELOCITY_DISPERSION
  for(PLINT ii = 0; ii < lkind; ii++){
    rad_l[ii] = _rad_l + ii * NUM_ANALYTIC;
    rho_l[ii] = _rho_l + ii * NUM_ANALYTIC;
  }/* for(PLINT ii = 0; ii < lkind; ii++){ */

  /* data preparation for particle distribution */
  for(PLINT ii = 0; ii < pkind; ii++)
    for(PLINT jj = 0; jj < npt[ii]; jj++){
      /* const int kk = group[grpIdx[ii]].prf_head + jj; */
      const int kk = group[grpIdx[ii]].prf_hor_head + jj;

      rad_p[ii][jj] = (PLFLT)log10((double)rad[kk] *      length2astro);
      rho_p[ii][jj] = (PLFLT)log10((double)rho[kk] * col_density2astro);
    }/* for(PLINT jj = 0; jj < npt[ii]; jj++){ */
  for(PLINT ii = 0; ii < ndisk; ii++){
    const int ll = nspheroids + ii;
    for(PLINT jj = 0; jj < npt[ll - skip]; jj++){
      /* const int kk = group[ll].prf_head + jj; */
      const int kk = group[ll].prf_hor_head + jj;

      hor_p[ii][jj] = (PLFLT)rad  [kk] * (PLFLT)length2astro;
      ver_p[ii][jj] = (PLFLT)zdisp[kk] * (PLFLT)length2astro;
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
      hor_l[ii][jj] = (PLFLT)rad  [kk] * (PLFLT)length2astro;
      ver_l[ii][jj] = (PLFLT)zdisp[kk] * (PLFLT)length2astro;
#endif//OVERPLOT_INITIAL_DISKHEIGHT

#ifdef  PLOT_VELOCITY_DISPERSION
      sigR_p[ii][jj] = (PLFLT)sigR[kk] * (PLFLT)velocity2astro;
      sigp_p[ii][jj] = (PLFLT)sigp[kk] * (PLFLT)velocity2astro;
      sigz_p[ii][jj] = (PLFLT)sigz[kk] * (PLFLT)velocity2astro;
#ifdef  PLOT_TOOMRE_Q
      toomre_p[ii][jj] = (PLFLT)toomre[kk];
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
    }/* for(PLINT jj = 0; jj < npt[ii]; jj++){ */
  }/* for(PLINT ii = 0; ii < ndisk; ii++){ */
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
  for(PLINT ii = 0; ii < ndisk; ii++)
    for(PLINT jj = 0; jj < hor_num_t0[nspheroids + ii]; jj++){
      const int kk = hor_head_t0[nspheroids + ii] + jj;

      hor_l[ndisk + ii][jj] =   rad_t0[kk] * (PLFLT)length2astro;
      ver_l[ndisk + ii][jj] = zdisp_t0[kk] * (PLFLT)length2astro;
    }/* for(PLINT jj = 0; jj < hor_num_t0[ll - skip]; jj++){ */
#endif//OVERPLOT_INITIAL_DISKHEIGHT
  /* data preparation for analytic distribution */
  for(PLINT ii = 0; ii < lkind; ii++)
    for(PLINT jj = 0; jj < nln[ii]; jj++){
      rad_l[ii][jj] = (PLFLT)log10((double)group[grpIdx[ii]].hor  [jj]);
      rho_l[ii][jj] = (PLFLT)log10((double)group[grpIdx[ii]].Sigma[jj]);
    }/* for(PLINT jj = 0; jj < nln[ii]; jj++){ */
#ifdef  PLOT_VELOCITY_DISPERSION
  for(PLINT ii = 0; ii < ndisk; ii++){
      const PLINT jj = ii + nspheroids;
      for(PLINT kk = 0; kk < NUM_ANALYTIC; kk++){
      sigR_l[ii][kk] = (PLFLT)group[grpIdx[jj]].disk_sigR[kk];
      sigp_l[ii][kk] = (PLFLT)group[grpIdx[jj]].disk_sigp[kk];
      sigz_l[ii][kk] = (PLFLT)group[grpIdx[jj]].disk_sigz[kk];
#ifdef  PLOT_TOOMRE_Q
      toomre_l[ii][kk] = (PLFLT)group[grpIdx[jj]].disk_toomre[kk];
#endif//PLOT_TOOMRE_Q
    }/* for(PLINT kk = 0; kk < NUM_ANALYTIC; kk++){ */
    }
#endif//PLOT_VELOCITY_DISPERSION

  /* set symbol style */
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
  PLplotLineStyle *ls_t0;  setDefaultLineStyle(ndisk * 2, &ls_t0);
  if( ndisk == 2 ){
    ls_t0[1].style = DASHED_LINE;
    ls_t0[2].style = DOT_DASHED_LINE;
    ls_t0[3].style = DOTTED_LINE;
  }/* if( ndisk == 2 ){ */
  for(int ii = 0; ii < ndisk; ii++){
    ls_t0[        ii].width =   BOLD_LINE;
    ls_t0[ndisk + ii].width = MIDDLE_LINE;
    ls_t0[ndisk + ii].color = ls_t0[ii].color;
  }/* for(int ii = 0; ii < ndisk; ii++){ */
#endif//OVERPLOT_INITIAL_DISKHEIGHT

  /* set labels */
  char radlab[PLplotCharWords];  sprintf(radlab, "#fiR #fr(%s)", length_astro_unit_name4plot);
  char rholab[PLplotCharWords];  sprintf(rholab, "#gS #fr(%s)", col_density_astro_unit_name4plot);
#ifdef  USE_RMS_THICKNESS
  char verlab[PLplotCharWords];  sprintf(verlab, "RMS(#fiz#fr) (%s)", length_astro_unit_name4plot);
#else///USE_RMS_THICKNESS
  char verlab[PLplotCharWords];  sprintf(verlab, "(<#fiz#fr#u2#d>-<#fiz#fr>#u2#d)#u1/2#d (%s)", length_astro_unit_name4plot);
#endif//USE_RMS_THICKNESS
#ifdef  PLOT_VELOCITY_DISPERSION
  char sigRlab[PLplotCharWords];  sprintf(sigRlab, "#gs#dR#u #fr(%s)", velocity_astro_unit_name4plot);
  char sigplab[PLplotCharWords];  sprintf(sigplab, "#gs#dp#u #fr(%s)", velocity_astro_unit_name4plot);
  char sigzlab[PLplotCharWords];  sprintf(sigzlab, "#gs#dz#u #fr(%s)", velocity_astro_unit_name4plot);
#ifdef  PLOT_TOOMRE_Q
  char toomreab[PLplotCharWords];  sprintf(toomrelab, "#fiQ");
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION

  /* set caption(s) */
  PLplotCaption rhocap;  setDefaultCaption(&rhocap);  rhocap.write = false;
  sprintf(rhocap.side, "%s", "t");
  sprintf(rhocap.text, "%s = %e (%s)", "#fit#fr", (double)time * time2astro, time_astro_unit_name4plot);
  PLplotCaption vercap;  setDefaultCaption(&vercap);  vercap.write = false;
  sprintf(vercap.side, "%s", "t");
  sprintf(vercap.text, "%s = %e (%s)", "#fit#fr", (double)time * time2astro, time_astro_unit_name4plot);
#ifdef  PLOT_VELOCITY_DISPERSION
  PLplotCaption sigcap;  setDefaultCaption(&sigcap);  sigcap.write = false;
  sprintf(sigcap.side, "%s", "t");
  sprintf(sigcap.text, "%s = %e (%s)", "#fit#fr", (double)time * time2astro, time_astro_unit_name4plot);
#ifdef  PLOT_TOOMRE_Q
  PLplotCaption Qcap;  setDefaultCaption(&Qcap);  Qcap.write = false;
  sprintf(Qcap.side, "%s", "t");
  sprintf(Qcap.text, "%s = %e (%s)", "#fit#fr", (double)time * time2astro, time_astro_unit_name4plot);
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION

  /* set legends */
  PLplotLegend rholeg;  setDefaultLegend(&rholeg, false);  rholeg.write = false;
  char *rholegTxt;
  {
    PLINT dkind = (unify) ? (IMAX(pkind, lkind)) : (pkind + lkind);
    allocChar4PLplot(&rholegTxt, dkind);    allocPointer4Char4PLplot(&(rholeg.text), dkind);    assignChar4PLplot(dkind, rholeg.text, rholegTxt);
  }
  PLplotLegend verleg;  setDefaultLegend(&verleg, false);  verleg.write = false;
  char *verlegTxt;
  {
    PLINT dkind = ndisk;
    allocChar4PLplot(&verlegTxt, dkind);    allocPointer4Char4PLplot(&(verleg.text), dkind);    assignChar4PLplot(dkind, verleg.text, verlegTxt);
  }
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
  PLplotLegend zt0leg;  setDefaultLegend(&zt0leg, false);  zt0leg.write = true;
  char *zt0legTxt;
  {
    PLINT dkind = ndisk * 2;
    allocChar4PLplot(&zt0legTxt, dkind);    allocPointer4Char4PLplot(&(zt0leg.text), dkind);    assignChar4PLplot(dkind, zt0leg.text, zt0legTxt);
    if( ndisk == 1 )
      for(int ii = 0; ii < dkind; ii++)
	sprintf(zt0leg.text[ii], "#fit#fr=%.1f (%s)", (((ii / ndisk) == 0) ? ((double)time) : 0.0) * time2astro, time_astro_unit_name4plot);
    else
      for(int ii = 0; ii < dkind; ii++)
	sprintf(zt0leg.text[ii], "disk%d (#fit#fr=%.1f %s)", ii % ndisk, (((ii / ndisk) == 0) ? ((double)time) : 0.0) * time2astro, time_astro_unit_name4plot);
  }
  zt0leg.pos = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;
#endif//OVERPLOT_INITIAL_DISKHEIGHT
#ifdef  PLOT_VELOCITY_DISPERSION
  PLplotLegend sigRleg;  setDefaultLegend(&sigRleg, false);  sigRleg.write = false;
  PLplotLegend sigpleg;  setDefaultLegend(&sigpleg, false);  sigpleg.write = false;
  PLplotLegend sigzleg;  setDefaultLegend(&sigzleg, false);  sigzleg.write = false;
  char *sigRlegTxt;
  char *sigplegTxt;
  char *sigzlegTxt;
  {
    PLINT dkind = (unify) ? (IMAX(pkind, lkind)) : (pkind + lkind);
    allocChar4PLplot(&sigRlegTxt, dkind);    allocPointer4Char4PLplot(&(sigRleg.text), dkind);    assignChar4PLplot(dkind, sigRleg.text, sigRlegTxt);
    allocChar4PLplot(&sigplegTxt, dkind);    allocPointer4Char4PLplot(&(sigpleg.text), dkind);    assignChar4PLplot(dkind, sigpleg.text, sigplegTxt);
    allocChar4PLplot(&sigzlegTxt, dkind);    allocPointer4Char4PLplot(&(sigzleg.text), dkind);    assignChar4PLplot(dkind, sigzleg.text, sigzlegTxt);
  }
#ifdef  PLOT_TOOMRE_Q
  PLplotLegend toomreleg;  setDefaultLegend(&toomreleg, false);  toomreleg.write = false;
  char *toomrelegTxt;
  {
    PLINT dkind = (unify) ? (IMAX(pkind, lkind)) : (pkind + lkind);
    allocChar4PLplot(&toomrelegTxt, dkind);    allocPointer4Char4PLplot(&(toomreleg.text), dkind);    assignChar4PLplot(dkind, toomreleg.text, toomrelegTxt);
  }
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION

  char title[PLplotCharWords];
  sprintf(title, "#fit#fr = %6.2f (%s)", (PLFLT)((double)time * time2astro), time_astro_unit_name4plot);


  /* memory allocation to exploit multiple-plot function */
  /* maximum number of panels in each direction */
#ifdef  PLOT_VELOCITY_DISPERSION
#ifndef PLOT_TOOMRE_Q
  const PLINT nxpanel = 2;
#else///PLOT_TOOMRE_Q
  const PLINT nxpanel = 3;
#endif//PLOT_TOOMRE_Q
#else///PLOT_VELOCITY_DISPERSION
  const PLINT nxpanel = 1;
#endif//PLOT_VELOCITY_DISPERSION
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
  const PLINT nypanel = 1 + (disk ? 2 : 0);
#else///OVERPLOT_INITIAL_DISKHEIGHT
  const PLINT nypanel = 1 + (disk ? 1 : 0);
#endif//OVERPLOT_INITIAL_DISKHEIGHT

  /* specify plot or skip panel */
  PLBOOL *plot;
  allocPLBOOL(&plot, nxpanel * nypanel);
  /* arrays related to draw line(s) and point(s) */
  PLINT *nlkind;  allocPLINT(&nlkind, nxpanel * nypanel);
  PLINT *npkind;  allocPLINT(&npkind, nxpanel * nypanel);
  PLINT **lnum;  allocPointer4PLINT(&lnum, nxpanel * nypanel);
  PLINT **pnum;  allocPointer4PLINT(&pnum, nxpanel * nypanel);
  PLFLT ***lx;  allocDoublePointer4PLFLT(&lx, nxpanel * nypanel);
  PLFLT ***ly;  allocDoublePointer4PLFLT(&ly, nxpanel * nypanel);
  PLFLT ***px;  allocDoublePointer4PLFLT(&px, nxpanel * nypanel);
  PLFLT ***py;  allocDoublePointer4PLFLT(&py, nxpanel * nypanel);
  PLplotLineStyle ** line;  allocPointer4PLplotLineStyle(& line, nxpanel * nypanel);
  PLplotPointType **point;  allocPointer4PLplotPointType(&point, nxpanel * nypanel);
  /* arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /* arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /* arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /* arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /* arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /* array to set figure name */
  /* char figfile[PLplotCharWords]; */
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel * nypanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel * nypanel);
  assignChar4PLplot(nxpanel * nypanel, figname, _figname);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* configuration about radial profile of particle distribution  */
  //-----------------------------------------------------------------------
  /* common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    for(PLINT jj = 0; jj < nypanel; jj++){
      //-------------------------------------------------------------------
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      //-------------------------------------------------------------------
      /* global setting(s) */
      plot[idx] = true;
      //-------------------------------------------------------------------
      /* point setting(s) */
      point [idx] = pt;
      //-------------------------------------------------------------------
      /* errorbar setting(s) */
      errbar[idx] = 0;
      //-------------------------------------------------------------------
      /* legend(s) */
      uni[idx] = unify;
      //-------------------------------------------------------------------
      /* label setting(s) */
      xlabel[idx] = radlab;
      //-------------------------------------------------------------------
    }
  }
  //-----------------------------------------------------------------------
  /* setting(s) for column density profile */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 0);
    //---------------------------------------------------------------------
    /* line setting(s) */
    nlkind[idx] = lkind;
    lnum  [idx] = nln;
    line  [idx] = ls;
    lx    [idx] = rad_l;
    ly    [idx] = rho_l;
    //---------------------------------------------------------------------
    /* point setting(s) */
    npkind[idx] = pkind;
    pnum  [idx] = npt;
    px    [idx] = rad_p;
    py    [idx] = rho_p;
    //---------------------------------------------------------------------
    /* plot area */
    range[idx] = rhobox;
    //---------------------------------------------------------------------
    /* label setting(s) */
    ylabel[idx] = rholab;
    //---------------------------------------------------------------------
    /* miscs */
    cap[idx] = rhocap;
    leg[idx] = rholeg;
    //---------------------------------------------------------------------
    /* file name */
    sprintf(figname[idx], "%s_%s_%.3d", file, "horizontal", filenum);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  /* setting(s) for scale height profile */
  if( disk ){
    for(PLINT ii = 0; ii < nxpanel; ii++){
      //-------------------------------------------------------------------
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 1);
      //-------------------------------------------------------------------
      /* line setting(s) */
      nlkind[idx] = 0;
      //-------------------------------------------------------------------
      /* point setting(s) */
      npkind[idx] = ndisk;
      pnum  [idx] = &npt[nspheroids - skip];
      px    [idx] = hor_p;
      py    [idx] = ver_p;
      //-------------------------------------------------------------------
      /* plot area */
      range[idx] = zbox;
      //-------------------------------------------------------------------
      /* label setting(s) */
      ylabel[idx] = verlab;
      //-------------------------------------------------------------------
      /* miscs */
      cap[idx] = vercap;
      leg[idx] = verleg;
      //-------------------------------------------------------------------
      /* file name */
      sprintf(figname[idx], "%s_%s_%.3d", file, "diskheight", filenum);
      //-------------------------------------------------------------------
    }
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
    for(PLINT ii = 0; ii < nxpanel; ii++){
      //-------------------------------------------------------------------
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 2);
      //-------------------------------------------------------------------
      /* point setting(s) */
      npkind[idx] = 0;
      //-------------------------------------------------------------------
      /* line setting(s) */
      nlkind[idx] = ndisk * 2;
      line  [idx] = ls_t0;
      lnum  [idx] = nln_t0;
      lx    [idx] = hor_l;
      ly    [idx] = ver_l;
      //-------------------------------------------------------------------
      /* plot area */
      range[idx] = zbox;
      //-------------------------------------------------------------------
      /* label setting(s) */
      ylabel[idx] = verlab;
      //-------------------------------------------------------------------
      /* miscs */
      cap[idx] = vercap;
      leg[idx] = zt0leg;
      //-------------------------------------------------------------------
      /* file name */
      sprintf(figname[idx], "%s_%s_%.3d", file, "diskheight_l", filenum);
      //-------------------------------------------------------------------
    }
#endif//OVERPLOT_INITIAL_DISKHEIGHT
  }/* if( disk ){ */
  //-----------------------------------------------------------------------
  /* individual file names */
  /* sprintf(figfile                                 , "%s.%s.%.3d", file, "profile", filenum); */
  /* sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 0)], "%s.%s.%.3d", file, "encmass", filenum); */
  /* sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 1)], "%s.%s.%.3d", file, "density", filenum); */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create figure(s) */
  //-----------------------------------------------------------------------
  for(PLINT idx = 0; idx < nxpanel * nypanel; idx++)
    plotData(1, 1, &plot[idx], false, false,
  	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
  	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
  	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
  	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
  	     &xlabel[idx], &ylabel[idx], "", figname[idx], argc, argv);
  /* plotData(nxpanel, nypanel, plot, true, true, */
  /* 	   nlkind,  line, lnum, lx, ly, */
  /* 	   npkind, point, pnum, px, py, */
  /* 	   errbar, NULL, NULL, NULL, NULL, NULL, NULL, */
  /* 	   cap, leg, uni, range, */
  /* 	   xlabel, ylabel, "", figfile, argc, argv); */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  //-----------------------------------------------------------------------
  free(_rad_p);  free(rad_p);  free(_rad_l);  free(rad_l);
  free(_rho_p);  free(rho_p);  free(_rho_l);  free(rho_l);
  free(pt);  free(ls);
  //-----------------------------------------------------------------------
  free(_hor_p);  free(hor_p);
  free(_ver_p);  free(ver_p);
  //-----------------------------------------------------------------------
  free(rholegTxt);
  free(verlegTxt);
  //-----------------------------------------------------------------------
  free(grpIdx);
  //-----------------------------------------------------------------------
#ifdef  OVERPLOT_INITIAL_DISKHEIGHT
  free(nln_t0);  free(_hor_l);  free(hor_l);  free(_ver_l);  free(ver_l);
  free(ls_t0);
  free(zt0legTxt);
#endif//OVERPLOT_INITIAL_DISKHEIGHT
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotVerticalProfile
(const int ntot, const int ngroup, const bool overlay_initial, model *group, PLFLT time,
 real *rad, real *rho, PLplotPltRange rhobox,
 char *file, const int filenum, int argc, char **argv)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* global setting(s) */
  //-----------------------------------------------------------------------
  const PLINT num_horizontal_bin = NUM_HORIZONTAL_BIN - 3;
  //-----------------------------------------------------------------------
  int skip = 0;
  for(int ii = 0; ii < ngroup; ii++)
    if( group[ii].num == 1 )
      skip++;
  //-----------------------------------------------------------------------
  /* maximum number of data kinds */
  const PLINT  pkind = (ngroup - skip) * num_horizontal_bin;
  const PLINT  lkind = overlay_initial ? ((ngroup - skip) * num_horizontal_bin) : 0;
  const PLBOOL unify = false;
  //-----------------------------------------------------------------------
  PLINT *grpIdx;  allocPLINT(&grpIdx, pkind);
  {
    int jj = 0;
    for(int ii = 0; ii < ngroup; ii++)
      if( group[ii].num > 1 ){
	grpIdx[jj] = ii;
	jj++;
      }/* if( group[ii].num > 1 ){ */
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* preparation for dataset */
  //-----------------------------------------------------------------------
  /* memory allocation for index */
  PLINT *npt;  allocPLINT(&npt, pkind);
  PLINT *nln;  allocPLINT(&nln, lkind);
  //-----------------------------------------------------------------------
  /* set number of data points */
  for(int ii = 0; ii < ngroup - skip; ii++)
    for(int jj = 0; jj < num_horizontal_bin; jj++)
      npt[ii * num_horizontal_bin + jj] = (PLINT)group[grpIdx[ii]].prf_ver_num[num_horizontal_bin - 1 - jj];
      /* npt[ii * num_horizontal_bin + jj] = (PLINT)group[ii].prf_ver_num[jj]; */
  for(PLINT ii = 0; ii < lkind; ii++)    nln[ii] = (PLINT)NUM_ANALYTIC;
  //-----------------------------------------------------------------------
  /* memory allocation for data */
  PLFLT *_rad_p;  allocPLFLT(&_rad_p,                 ntot);  PLFLT **rad_p;  allocPointer4PLFLT(&rad_p, pkind);
  PLFLT *_rho_p;  allocPLFLT(&_rho_p,                 ntot);  PLFLT **rho_p;  allocPointer4PLFLT(&rho_p, pkind);
  PLFLT *_rad_l;  allocPLFLT(&_rad_l, NUM_ANALYTIC * lkind);  PLFLT **rad_l;  allocPointer4PLFLT(&rad_l, lkind);
  PLFLT *_rho_l;  allocPLFLT(&_rho_l, NUM_ANALYTIC * lkind);  PLFLT **rho_l;  allocPointer4PLFLT(&rho_l, lkind);
  /* set pointers */
  rad_p[0] = _rad_p;
  rho_p[0] = _rho_p;
  for(PLINT ii = 1; ii < pkind; ii++){
    rad_p[ii] = rad_p[ii - 1] + npt[ii - 1];
    rho_p[ii] = rho_p[ii - 1] + npt[ii - 1];
  }/* for(PLINT ii = 1; ii < pkind; ii++){ */
  for(PLINT ii = 0; ii < lkind; ii++){
    rad_l[ii] = _rad_l + ii * NUM_ANALYTIC;
    rho_l[ii] = _rho_l + ii * NUM_ANALYTIC;
  }/* for(PLINT ii = 0; ii < lkind; ii++){ */
  /* data preparation */
  for(int ii = 0; ii < ngroup - skip; ii++)
    for(int jj = 0; jj < num_horizontal_bin; jj++){
      const int kk = ii * num_horizontal_bin + jj;
      const int nn = num_horizontal_bin - 1 - jj;
      /* data preparation for particle distribution */
      for(PLINT ll = 0; ll < npt[kk]; ll++){
	const int mm = group[grpIdx[ii]].prf_ver_head[nn] + ll;
	rad_p[kk][ll] = (PLFLT)log10((double)rad[mm] *  length2astro);
	rho_p[kk][ll] = (PLFLT)log10((double)rho[mm] * density2astro);
      }/* for(PLINT ll = 0; ll < npt[kk]; ll++){ */
      /* data preparation for analytic distribution */
      if( overlay_initial )
	for(PLINT ll = 0; ll < nln[kk]; ll++){
	  rad_l[kk][ll] = (PLFLT)log10((double)group[grpIdx[ii]].ver_pos[nn][ll]);
	  rho_l[kk][ll] = (PLFLT)log10((double)group[grpIdx[ii]].ver_rho[nn][ll]);
	}/* for(PLINT ll = 0; ll < nln[kk]; ll++){ */
    }/* for(int jj = 0; jj < num_horizontal_bin; jj++){ */
  //-----------------------------------------------------------------------
  /* set symbol style */
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  for(int ii = 0; ii < ngroup - skip; ii++)
    for(int jj = 0; jj < num_horizontal_bin; jj++){
      const int kk = ii * num_horizontal_bin + jj;
      pt[kk].color =       PLplotColor     [(num_horizontal_bin - 1 - jj) % 15];
      pt[kk].scale =       PLplotSymbolSize[(num_horizontal_bin - 1 - jj) % 14] * 0.5;
      sprintf(pt[kk].type, "%s", PLplotSymbolType[(num_horizontal_bin - 1 - jj) % 14]);
      /* pt[kk].type  = (char *)PLplotSymbolType[(num_horizontal_bin - 1 - jj) % 14]; */
      if( overlay_initial ){
	ls[kk].color = PLplotColor[(num_horizontal_bin - 1 - jj) % 15];
	ls[kk].style =             (num_horizontal_bin - 1 - jj) %  5;
	ls[kk].width = THIN_LINE;
      }/* if( overlay_initial ){ */
    }/* for(int jj = 0; jj < num_horizontal_bin; jj++){ */
  //-----------------------------------------------------------------------
  /* set labels */
  char radlab[PLplotCharWords];  sprintf(radlab, "#fiz #fr(%s)", length_astro_unit_name4plot);
  char rholab[PLplotCharWords];  sprintf(rholab, "#fi#gr #fr(%s)", density_astro_unit_name4plot);
  //-----------------------------------------------------------------------
  /* set caption(s) */
  PLplotCaption rhocap;  setDefaultCaption(&rhocap);  rhocap.write = false;
  sprintf(rhocap.side, "%s", "t");
  sprintf(rhocap.text, "%s = %e (%s)", "#fit#fr", (double)time * time2astro, time_astro_unit_name4plot);
  //-----------------------------------------------------------------------
  /* set legends */
  PLplotLegend rholeg;  setDefaultLegend(&rholeg, false);  rholeg.write = false;
  char *rholegTxt;
  {
    PLINT dkind = (unify) ? (IMAX(pkind, lkind)) : (pkind + lkind);
    allocChar4PLplot(&rholegTxt, dkind);    allocPointer4Char4PLplot(&(rholeg.text), dkind);    assignChar4PLplot(dkind, rholeg.text, rholegTxt);
  }
  /* sprintf(rholeg.text[0], "%s", "halo"); */
  /* sprintf(rholeg.text[1], "%s", "bulge"); */
  /* sprintf(rholeg.text[2], "%s", "disk"); */
  //-----------------------------------------------------------------------
  char title[PLplotCharWords];
  sprintf(title, "#fit#fr = %6.2f (%s)", (PLFLT)((double)time * time2astro), time_astro_unit_name4plot);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation to exploit multiple-plot function */
  //-----------------------------------------------------------------------
  /* maximum number of panels in each direction */
  const PLINT nxpanel = 1;
  const PLINT nypanel = ngroup - skip;
  //-----------------------------------------------------------------------
  /* specify plot or skip panel */
  PLBOOL *plot;
  allocPLBOOL(&plot, nxpanel * nypanel);
  /* arrays related to draw line(s) and point(s) */
  PLINT *nlkind;  allocPLINT(&nlkind, nxpanel * nypanel);
  PLINT *npkind;  allocPLINT(&npkind, nxpanel * nypanel);
  PLINT **lnum;  allocPointer4PLINT(&lnum, nxpanel * nypanel);
  PLINT **pnum;  allocPointer4PLINT(&pnum, nxpanel * nypanel);
  PLFLT ***lx;  allocDoublePointer4PLFLT(&lx, nxpanel * nypanel);
  PLFLT ***ly;  allocDoublePointer4PLFLT(&ly, nxpanel * nypanel);
  PLFLT ***px;  allocDoublePointer4PLFLT(&px, nxpanel * nypanel);
  PLFLT ***py;  allocDoublePointer4PLFLT(&py, nxpanel * nypanel);
  PLplotLineStyle ** line;  allocPointer4PLplotLineStyle(& line, nxpanel * nypanel);
  PLplotPointType **point;  allocPointer4PLplotPointType(&point, nxpanel * nypanel);
  /* arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /* arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /* arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /* arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /* arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /* array to set figure name */
  /* char figfile[PLplotCharWords]; */
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel * nypanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel * nypanel);
  assignChar4PLplot(nxpanel * nypanel, figname, _figname);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* configuration about radial profile of particle distribution  */
  //-----------------------------------------------------------------------
  /* common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    for(PLINT jj = 0; jj < nypanel; jj++){
      //-------------------------------------------------------------------
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      //-------------------------------------------------------------------
      /* global setting(s) */
      plot[idx] = true;
      //-------------------------------------------------------------------
      /* line setting(s) */
      nlkind[idx] = overlay_initial ? num_horizontal_bin : 0;
      lnum  [idx] = &nln  [overlay_initial ? (idx * num_horizontal_bin) : 0];
      line  [idx] = &ls   [overlay_initial ? (idx * num_horizontal_bin) : 0];
      lx    [idx] = &rad_l[overlay_initial ? (idx * num_horizontal_bin) : 0];
      ly    [idx] = &rho_l[overlay_initial ? (idx * num_horizontal_bin) : 0];
      //-------------------------------------------------------------------
      /* point setting(s) */
      npkind[idx] = num_horizontal_bin;
      pnum  [idx] = &npt  [idx * num_horizontal_bin];
      point [idx] = &pt   [idx * num_horizontal_bin];
      px    [idx] = &rad_p[idx * num_horizontal_bin];
      py    [idx] = &rho_p[idx * num_horizontal_bin];
      //-------------------------------------------------------------------
      /* errorbar setting(s) */
      errbar[idx] = 0;
      //-------------------------------------------------------------------
      /* plot area */
      range[idx] = rhobox;
      //-------------------------------------------------------------------
      /* legend(s) */
      uni[idx] = unify;
      //-------------------------------------------------------------------
      /* label setting(s) */
      xlabel[idx] = radlab;
      ylabel[idx] = rholab;
      //-------------------------------------------------------------------
      /* miscs */
      cap[idx] = rhocap;
      leg[idx] = rholeg;
      //-------------------------------------------------------------------
      /* file name */
      sprintf(figname[idx], "%s_%s_%s%d_%.3d", file, "vertical", "series", idx, filenum);
      //-------------------------------------------------------------------
    }
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create figure(s) */
  //-----------------------------------------------------------------------
  for(PLINT idx = 0; idx < nxpanel * nypanel; idx++)
    plotData(1, 1, &plot[idx], false, false,
  	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
  	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
  	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
  	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
  	     &xlabel[idx], &ylabel[idx], "", figname[idx], argc, argv);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  //-----------------------------------------------------------------------
  free(_rad_p);  free(rad_p);  free(_rad_l);  free(rad_l);
  free(_rho_p);  free(rho_p);  free(_rho_l);  free(rho_l);
  free(pt);  free(ls);
  //-----------------------------------------------------------------------
  free(rholegTxt);
  //-----------------------------------------------------------------------
  free(grpIdx);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


void analyzeRadialProfile
(int Np, nbody_particle *body_tot, const int kind, model *group, int *num, int *rem, real **rad, real **rho, real **enc,
#ifdef  PLOT_VELOCITY_DISPERSION
 real **sigr,
#endif//PLOT_VELOCITY_DISPERSION
 const real r2max)
{
  __NOTE__("%s\n", "start");


  /** shift the center to the center-of-mass of the system */
  /* calculate the center-of-mass of the system */
  double Mtot_dbl = 0.0;
  double comx_dbl = 0.0, comy_dbl = 0.0, comz_dbl = 0.0;
#ifdef  PLOT_VELOCITY_DISPERSION
  double velx_dbl = 0.0, vely_dbl = 0.0, velz_dbl = 0.0;
#endif//PLOT_VELOCITY_DISPERSION

  for(int ii = 0; ii < Np; ii++){
    const double xx = CAST_R2D(body_tot[ii].x);
    const double yy = CAST_R2D(body_tot[ii].y);
    const double zz = CAST_R2D(body_tot[ii].z);
    const double r2 = xx * xx + yy * yy + zz * zz;
    if( r2 < 4.0 * r2max ){
      const double mass = CAST_R2D(body_tot[ii].m);
      Mtot_dbl += mass;

      comx_dbl += mass * xx;
      comy_dbl += mass * yy;
      comz_dbl += mass * zz;

#ifdef  PLOT_VELOCITY_DISPERSION
      velx_dbl += mass * CAST_R2D(body_tot[ii].vx);
      vely_dbl += mass * CAST_R2D(body_tot[ii].vy);
      velz_dbl += mass * CAST_R2D(body_tot[ii].vz);
#endif//PLOT_VELOCITY_DISPERSION
    }
  }
  Mtot_dbl = 1.0 / Mtot_dbl;
  real com[3] = {CAST_D2R(comx_dbl * Mtot_dbl), CAST_D2R(comy_dbl * Mtot_dbl), CAST_D2R(comz_dbl * Mtot_dbl)};
#ifdef  PLOT_VELOCITY_DISPERSION
  real vel[3] = {CAST_D2R(velx_dbl * Mtot_dbl), CAST_D2R(vely_dbl * Mtot_dbl), CAST_D2R(velz_dbl * Mtot_dbl)};
#endif//PLOT_VELOCITY_DISPERSION
  /* calculate the distance from the center-of-mass of the system */
  for(int ii = 0; ii < Np; ii++){
    body_tot[ii].x -= com[0];
    body_tot[ii].y -= com[1];
    body_tot[ii].z -= com[2];

    const real R2 = body_tot[ii].x * body_tot[ii].x + body_tot[ii].y * body_tot[ii].y;
    const real z2 = body_tot[ii].z * body_tot[ii].z;
    body_tot[ii].ax = R2 + z2;
    body_tot[ii].ay = R2;
    body_tot[ii].az = z2;

#ifdef  PLOT_VELOCITY_DISPERSION
    body_tot[ii].vx -= vel[0];
    body_tot[ii].vy -= vel[1];
    body_tot[ii].vz -= vel[2];
#endif//PLOT_VELOCITY_DISPERSION
  }


  /** make radial profile */
  *num = 0;
  for(int kk = 0; kk < kind; kk++){
    group[kk].prf_head = *num;

    /* sort the array in ascending distance from the center */
    nbody_particle *body;
    body = &body_tot[group[kk].head];
    qsort(body, group[kk].num, sizeof(nbody_particle), radAscendingOrder);

    real *prad, *prho, *penc;
#ifdef  PLOT_VELOCITY_DISPERSION
    real *psig;
#endif//PLOT_VELOCITY_DISPERSION
    real inner = ZERO;
    real Menc  = ZERO;
    prad = *rad;
    prho = *rho;
    penc = *enc;
#ifdef  PLOT_VELOCITY_DISPERSION
    psig = *sigr;
#endif//PLOT_VELOCITY_DISPERSION
    const int ncrit = ((int)(group[kk].num / ncrit_base) > Nminimum) ? ncrit_base : (int)(group[kk].num / Nminimum);
#ifdef  PLOT_VELOCITY_DISPERSION
    const real inv_ncrit = UNITY / (real)ncrit;
#endif//PLOT_VELOCITY_DISPERSION
    for(int head = 0; head < (int)group[kk].num; head += ncrit){
      /** check # of unused elements */
      if( *rem == 0 ){
	enlargeProfileArray(*num + NAllocUnit, rad, rho, enc
#ifdef  PLOT_VELOCITY_DISPERSION
			    , sigr
#endif//PLOT_VELOCITY_DISPERSION
			    );
	*rem += NAllocUnit;
	prad = *rad;
	prho = *rho;
	penc = *enc;
#ifdef  PLOT_VELOCITY_DISPERSION
	psig = *sigr;
#endif//PLOT_VELOCITY_DISPERSION
      }

      int tail = head + ncrit;
      if( (ulong)tail - 1 >= group[kk].num )	break;

      real outer = SQRT(body[tail - 1].ax);
#ifdef  USE_MEDIAN_FOR_POSITION
      if( ncrit % 1 ){
	/* if ncrit is an odd number */
	prad[*num] = SQRT(body[head + (ncrit >> 1)].ax);
      }
      else{
	/* if ncrit is an even number */
	prad[*num] = HALF * (SQRT(body[head + (ncrit >> 1) - 1].ax) + SQRT(body[head + (ncrit >> 1)].ax));
      }
#else///USE_MEDIAN_FOR_POSITION
      prad[*num] = HALF * (inner + outer);
#endif//USE_MEDIAN_FOR_POSITION

#ifdef  PLOT_VELOCITY_DISPERSION
      real vr2 = ZERO;
      real vrm = ZERO;
      for(int ii = head; ii < tail; ii++){
	const real vr = (body_tot[ii].vx * body_tot[ii].x + body_tot[ii].vy * body_tot[ii].y + body_tot[ii].vz * body_tot[ii].z) / sqrt(body_tot[ii].ax);	/**< ax = r2; */
	vrm += vr;
	vr2 += vr * vr;
      }/* for(int ii = head; ii < tail; ii++){ */
      vrm *= inv_ncrit;
      vr2 *= inv_ncrit;
      psig[*num] = SQRT(vr2 - vrm * vrm);
#endif//PLOT_VELOCITY_DISPERSION

      real mass = ZERO;
      for(int ii = head; ii < tail; ii++)
	mass += body[ii].m;
      prho[*num] = mass / ((real)(4.0 * M_PI / 3.0) * (outer * outer * outer - inner * inner * inner));
      penc[*num] = Menc + HALF * mass;

      Menc += mass;
      inner = outer;

      *num += 1;
      *rem -= 1;
    }/* for(int head = 0; head < (int)group[kk].num; head += ncrit){ */

    group[kk].prf_num = *num - (group[kk].prf_head);
  }


  __NOTE__("%s\n", "end");
}


void analyzeDecomposedProfile(nbody_particle *body_tot, const int kind, model *group,
			      int *numHor, int *remHor, real **hor_pos, real **hor_rho, real **hor_zdisp,
			      int *numVer, int *remVer, real **ver_pos, real **ver_rho
#ifdef  PLOT_VELOCITY_DISPERSION
			      , real **sigR, real **sigp, real **sigz
#ifdef  PLOT_TOOMRE_Q
			      , real **toomre
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
			      )
{
  __NOTE__("%s\n", "start");


  *numHor = 0;
  *numVer = 0;

  for(int kk = 0; kk < kind; kk++){
    real *ppos, *prho, *pzdisp;
#ifdef  PLOT_VELOCITY_DISPERSION
    real *psigR, *psigp, *psigz;
#ifdef  PLOT_TOOMRE_Q
    real *ptoomre;
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION


    /** make horizontal profile on the midplane */
    /* sort the array in ascending distance from the center on the midplane */
    nbody_particle *body;
    body = &body_tot[group[kk].head];
    qsort(body, group[kk].num, sizeof(nbody_particle), horAscendingOrder);

    group[kk].prf_hor_head = *numHor;
    real inner = ZERO;
    ppos   = *hor_pos;
    prho   = *hor_rho;
    pzdisp = *hor_zdisp;
#ifdef  PLOT_VELOCITY_DISPERSION
    psigR   = *sigR;
    psigp   = *sigp;
    psigz   = *sigz;
#ifdef  PLOT_TOOMRE_Q
    ptoomre = *toomre;
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION

    const int ncrit = ((int)(group[kk].num / ncrit_base) > Nminimum) ? ncrit_base : (int)(group[kk].num / Nminimum);
#ifdef  PLOT_VELOCITY_DISPERSION
    const real inv_ncrit = UNITY / (real)ncrit;
#endif//PLOT_VELOCITY_DISPERSION
    for(int head = 0; head < (int)group[kk].num; head += ncrit){
      /* check # of unused elements */
      if( *remHor == 0 ){
      	enlargeDecomposedProfileArray(*numHor + NAllocUnit, hor_pos, hor_rho, true, hor_zdisp
#ifdef  PLOT_VELOCITY_DISPERSION
				      , true, sigR, sigp, sigz
#ifdef  PLOT_TOOMRE_Q
				      , toomre
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
				      );
      	*remHor += NAllocUnit;
      	ppos   = *hor_pos;
      	prho   = *hor_rho;
	pzdisp = *hor_zdisp;
#ifdef  PLOT_VELOCITY_DISPERSION
	psigR   = *sigR;
	psigp   = *sigp;
	psigz   = *sigz;
#ifdef  PLOT_TOOMRE_Q
	ptoomre = *toomre;
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
      }

      int tail = head + ncrit;
      if( (ulong)tail - 1 >= group[kk].num )	break;

      real outer = SQRT(body[tail - 1].ay);
#ifdef  USE_MEDIAN_FOR_POSITION
      if( ncrit % 1 ){
	/* if ncrit is an odd number */
	ppos[*numHor] = SQRT(body[head + (ncrit >> 1)].ay);
      }
      else{
	/* if ncrit is an even number */
	ppos[*numHor] = HALF * (SQRT(body[head + (ncrit >> 1) - 1].ay) + SQRT(body[head + (ncrit >> 1)].ay));
      }
#else///USE_MEDIAN_FOR_POSITION
      ppos[*numHor] = HALF * (inner + outer);
#endif//USE_MEDIAN_FOR_POSITION

      real mass = ZERO;
      real zvar = ZERO;
#ifndef USE_RMS_THICKNESS
      real mean = ZERO;
#endif//USE_RMS_THICKNESS
      for(int ii = head; ii < tail; ii++){
      	mass += body[ii].m;
	zvar += body[ii].z * body[ii].z;
#ifndef USE_RMS_THICKNESS
	mean += body[ii].z;
#endif//USE_RMS_THICKNESS
      }/* for(int ii = head; ii < tail; ii++){ */
      prho  [*numHor] = mass / ((real)(M_PI) * (outer * outer - inner * inner));

#ifndef USE_RMS_THICKNESS
      mass = UNITY / (real)(tail - head);
      mean *= mass;
      pzdisp[*numHor] = SQRT(zvar * mass - mean * mean);
#else///USE_RMS_THICKNESS
      pzdisp[*numHor] = SQRT(zvar / (real)(tail - head));
#endif//USE_RMS_THICKNESS


#ifdef  USE_ITERATED_THICKNESS
      int numOld = tail - head;
      while( true ){
	/* set upper limit as 5 sigma */
	const real trial  = FIVE * pzdisp[*numHor];
	const real trial2 = trial * trial;

	zvar = ZERO;
#ifndef USE_RMS_THICKNESS
	const real mold = mean;
	mean = ZERO;
#endif//USE_RMS_THICKNESS
	int numValid = 0;
	for(int ii = head; ii < tail; ii++){
#ifdef  USE_RMS_THICKNESS
	  if(   ( body[ii].z         *  body[ii].z        ) < trial2 )
#else///USE_RMS_THICKNESS
	    if( ((body[ii].z - mold) * (body[ii].z - mold)) < trial2 )
#endif//USE_RMS_THICKNESS
	      {
		numValid++;
		zvar += body[ii].z * body[ii].z;
#ifndef USE_RMS_THICKNESS
		mean += body[ii].z;
#endif//USE_RMS_THICKNESS
	      }
	}/* for(int ii = head; ii < tail; ii++){ */

#ifndef USE_RMS_THICKNESS
	mass = UNITY / (real)numValid;
	mean *= mass;
	pzdisp[*numHor] = SQRT(zvar * mass - mean * mean);
#else///USE_RMS_THICKNESS
	pzdisp[*numHor] = SQRT(zvar / (real)numValid);
#endif//USE_RMS_THICKNESS

	if( numValid == numOld )
	  break;
	else
	  numOld = numValid;
      }
#endif//USE_ITERATED_THICKNESS

#ifdef  PLOT_VELOCITY_DISPERSION
      real vR2 = ZERO;      real vRm = ZERO;
      real vp2 = ZERO;      real vpm = ZERO;
      real vz2 = ZERO;      real vzm = ZERO;
      for(int ii = head; ii < tail; ii++){
	const real vz = body_tot[ii].vz;
	const real invR = SQRT(body_tot[ii].ay);	/**< ay = R2; */
	const real xx = body_tot[ii].x;
	const real yy = body_tot[ii].y;
	const real vx = body_tot[ii].vx;
	const real vy = body_tot[ii].vy;
	vzm += vz;
	vz2 += vz * vz;
	const real vR = ( xx * vx + yy * vy) * invR;
	const real vp = (-yy * vx + xx * vy) * invR;
	vRm += vR;
	vR2 += vR * vR;
	vpm += vp;
	vp2 += vp * vp;
      }/* for(int ii = hea; ii < tail; ii++){ */
      vzm *= inv_ncrit;      vz2 *= inv_ncrit;
      vRm *= inv_ncrit;      vR2 *= inv_ncrit;
      vpm *= inv_ncrit;      vp2 *= inv_ncrit;
      psigz[*numHor] = SQRT(vz2 - vzm * vzm);
      psigR[*numHor] = SQRT(vR2 - vRm * vRm);
      psigp[*numHor] = SQRT(vp2 - vpm * vpm);

#ifdef  PLOT_TOOMRE_Q
      /* const double toomre = sigmaR * kappa / (DBL_MIN + 3.36 * newton * Sigma); */
      ptoomre[*numHor] = ;
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION

      inner = outer;
      *numHor += 1;
      *remHor -= 1;
    }/* for(int head = 0; head < (int)group[kk].num; head += ncrit){ */
    group[kk].prf_hor_num = *numHor - (group[kk].prf_hor_head);


    /** make vertical profile in each R bin */
    const int numSub = (int)group[kk].num / NUM_HORIZONTAL_BIN;
    for(int ii = 0; ii < NUM_HORIZONTAL_BIN; ii++){
      /** sort the array in ascending distance from the the midplane */
      nbody_particle *data;
      data = &body_tot[numSub * ii + (int)group[kk].head];
      const real R2inner = data[         0].x * data[         0].x + data[         0].y * data[         0].y;
      const real R2outer = data[numSub - 1].x * data[numSub - 1].x + data[numSub - 1].y * data[numSub - 1].y;
      const real piR2 = (real)(M_PI * (double)(R2outer - R2inner));
      qsort(data, numSub, sizeof(nbody_particle), verAscendingOrder);

      group[kk].prf_ver_head[ii] = *numVer;
      real lower = ZERO;
      ppos = *ver_pos;
      prho = *ver_rho;

      for(int head = 0; head < numSub; head += ncrit){
	/** check # of unused elements */
	if( *remVer == 0 ){
	  enlargeDecomposedProfileArray(*numVer + NAllocUnit, ver_pos, ver_rho, false, NULL
#ifdef  PLOT_VELOCITY_DISPERSION
					, false, NULL, NULL, NULL
#ifdef  PLOT_TOOMRE_Q
					, NULL
#endif//PLOT_TOOMRE_Q
#endif//PLOT_VELOCITY_DISPERSION
					);
	  *remVer += NAllocUnit;
	  ppos = *ver_pos;
	  prho = *ver_rho;
	}

	int tail = head + ncrit;
	if( tail - 1 >= numSub )	  break;

	real upper = FABS(data[tail - 1].z);
#ifdef  USE_MEDIAN_FOR_POSITION
	if( ncrit % 1 ){
	  /* if ncrit is an odd number */
	  ppos[*numVer] = SQRT(body[head + (ncrit >> 1)].z);
	}
	else{
	  /* if ncrit is an even number */
	  ppos[*numVer] = HALF * (SQRT(body[head + (ncrit >> 1) - 1].z) + SQRT(body[head + (ncrit >> 1)].z));
	}
#else///USE_MEDIAN_FOR_POSITION
	ppos[*numVer] = HALF * (lower + upper);
#endif//USE_MEDIAN_FOR_POSITION
	real mass = ZERO;
	for(int jj = head; jj < tail; jj++)
	  mass += data[jj].m;
	prho[*numVer] = HALF * mass / (piR2 * (upper - lower));
	/* both z > 0 and z < 0 components are summed, symmetry about |z| is assumed */
	lower = upper;
	*numVer += 1;
	*remVer -= 1;
      }/* for(int head = 0; head < numSub; head += ncrit){ */
      group[kk].prf_ver_num[ii] = *numVer - (group[kk].prf_ver_head[ii]);
    }/* for(int ii = 0; ii < NUM_HORIZONTAL_BIN; ii++){ */
  }/* for(int kk = 0; kk < kind; kk++){ */


  __NOTE__("%s\n", "end");
}


void integrateColumnDensity(const int num, const real logRbin, real *rad, real *Sigma, real *mass)
{
  __NOTE__("%s\n", "start");


  /* calculate enclosed mass within a cylinder */
  real Menc[2];
  Menc[0] = rad[0] * rad[0] * Sigma[0];
  Menc[1] = rad[1] * rad[1] * Sigma[1];
  mass[0] = Menc[0];
  mass[1] = mass[0] + (rad[1] - rad[0]) * (rad[1] - rad[0]) * TWO * (Sigma[0] + Sigma[1]);
  Menc[0] += Menc[1] * FOUR;

  for(int ii = 2; ii < num; ii++){
    const real mtmp = rad[ii] * rad[ii] * Sigma[ii];
    const int idx = ii & 1;

    mass[ii] = mass[idx] + (mtmp + Menc[idx]) * logRbin * (real)M_LN10 * ONE_THIRD;

    Menc[0] += mtmp * (real)(1 << (1 + (idx    )));
    Menc[1] += mtmp * (real)(1 << (1 + (idx ^ 1)));
  }

  /* multiply overall factor */
  for(int ii = 0; ii < num; ii++)
    mass[ii] *= TWO * (real)M_PI;


  __NOTE__("%s\n", "end");
}
