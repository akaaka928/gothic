/**
 * @file accretion.c
 *
 * @brief Analyzer for mass accretion related quantities (mass, timescale, and accretion rate)
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2019/07/31 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

/**
 * @def USE_SZIP_COMPRESSION
 *
 * @brief On to enable Szip compression for HDF5 files (default is ON).
 */
/* #define USE_SZIP_COMPRESSION */

/**
 * @def USE_GZIP_COMPRESSION
 *
 * @brief On to enable gzip compression for HDF5 files (default is ON).
 *
 * @detail currently, h5py does not accept Szip compression in default.
 */
#define USE_GZIP_COMPRESSION

#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#define USE_FILE_COMPRESSION
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
/* The maximum number of elements in a chunk is 2^32-1 which is equal to 4,294,967,295 */
/* The maximum size for any chunk is 4GB */
#define MAXIMUM_CHUNK_SIZE      ((hsize_t)1 << 31)
#define MAXIMUM_CHUNK_SIZE_4BIT ((hsize_t)1 << 30)
#define MAXIMUM_CHUNK_SIZE_8BIT ((hsize_t)1 << 29)
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "myutil.h"
#include "name.h"
#include "constants.h"

#include "../init/profile.h"


/* const double criteria_abs = 1.0e-30; */
/* const double criteria_abs = 1.0e-20; */
/* const double criteria_abs = 1.0e-16; */
/* const double criteria_abs = 1.0e-14; */
/* const double criteria_abs = 1.0e-12; */
const double criteria_abs = 1.0e-10;
/* const double criteria_rel = 1.0e-12; */
/* const double criteria_rel = 1.0e-10; */
/* const double criteria_rel = 1.0e-8; */
/* const double criteria_rel = 1.0e-6; */
const double criteria_rel = 1.0e-5;
/* const double criteria_rel = 1.0e-4; */


static inline double getDF(const double QQ, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  if( QQ < Qmin )
    return (0.0);
  const int ii = (int)((QQ - Qmin) * invQbin);
  return (CAST_R2D(df[ii] + (df[ii + 1] - df[ii]) * (QQ - ene[ii]) / (ene[ii + 1] - ene[ii])));
}
static inline double getPsi(const double rad, const double logrmin, const double invlogrbin, real * restrict Psi)
{
  const double rr = (log10(rad) - logrmin) * invlogrbin;
  const int ii = (int)rr;
  if( ii < 0 )
    return (CAST_R2D(Psi[0]));
  return (CAST_R2D(Psi[ii]) + (CAST_R2D(Psi[ii + 1]) - CAST_R2D(Psi[ii])) * (rr - (double)ii));
}


static inline double get_DEformula_azimuthal(const double tt, const double Psi, const double r2_ra2, const double v2_2, const double alp_2, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double cosh_t = cosh(sinh_t);
  const double inv_cosh_t = 1.0 / cosh_t;

  const double sineta = sin(alp_2 * exp(sinh_t) * inv_cosh_t);
  const double QQ = Psi - v2_2 * (1.0 + r2_ra2 * sineta * sineta);

  return (cosh(tt) * inv_cosh_t * inv_cosh_t * getDF(QQ, ene, df, Qmin, invQbin));
}
static inline double update_trapezoidal_azimuthal(const double hh, const double tmin, const double tmax, const double sum, const double Psi, const double r2_ra2, const double v2_2, const double alp_2, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  /** initialization */
  double sub = 0.0;
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    sub += get_DEformula_azimuthal(tt, Psi, r2_ra2, v2_2, alp_2, ene, df, Qmin, invQbin);
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  return (0.5 * sum + hh * sub);
}
static inline double set_domain_boundary_azimuthal(const double hh, double * restrict tmin, double * restrict tmax, const double Psi, const double r2_ra2, const double v2_2, const double alp_2, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  double fp = get_DEformula_azimuthal(tt, Psi, r2_ra2, v2_2, alp_2, ene, df, Qmin, invQbin);
  const double f0 = fp;
  double sum = hh * fp;


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    const double ft = fp;

    tt += hh;    boundary = tt;
    fp = get_DEformula_azimuthal(tt, Psi, r2_ra2, v2_2, alp_2, ene, df, Qmin, invQbin);

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
    fp = get_DEformula_azimuthal(tt, Psi, r2_ra2, v2_2, alp_2, ene, df, Qmin, invQbin);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;

  return (sum);
}
static inline double integrate_DEformula_azimuthal(const double lmax, const double rad, const double Psi, const double r2_ra2, const double vel, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  const double alp_2 = 0.5 * asin(lmax / fmax(rad * vel, lmax));
#ifndef NDEBUG
  if( fpclassify(alp_2) != FP_NORMAL ){
    __KILL__(stderr, "ERROR: alp_2 = %e: lmax = %e, rad = %e, vel = %e; x = %e\n", alp_2, lmax, rad, vel, lmax / (rad * vel));
  }/* if( fpclassify(alp_2) != FP_NORMAL ){ */
#endif//NDEBUG
  const double v2_2 = 0.5 * vel * vel;

  double hh = 1.0;
  double tmin, tmax;
  double sum = set_domain_boundary_azimuthal(hh, &tmin, &tmax, Psi, r2_ra2, v2_2, alp_2, ene, df, Qmin, invQbin);

  while( true ){
    const double f0 = sum;

    hh *= 0.5;
    sum = update_trapezoidal_azimuthal(hh, tmin, tmax, sum, Psi, r2_ra2, v2_2, alp_2, ene, df, Qmin, invQbin);

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

  return (alp_2 * sum);
}


static inline double get_DEformula_velocity(const double tt, const double vesc_2, const double lmax, const double rad, const double Psi, const double r2_ra2, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double cosh_t = cosh(sinh_t);
  const double inv_cosh_t = 1.0 / cosh_t;

  const double vel = vesc_2 * exp(sinh_t) * inv_cosh_t;

  return (cosh(tt) * inv_cosh_t * inv_cosh_t * vel * vel * integrate_DEformula_azimuthal(lmax, rad, Psi, r2_ra2, vel, ene, df, Qmin, invQbin));
}
static inline double update_trapezoidal_velocity(const double hh, const double tmin, const double tmax, const double sum, const double vesc_2, const double lmax, const double rad, const double Psi, const double r2_ra2, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  /** initialization */
  double sub = 0.0;
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    sub += get_DEformula_velocity(tt, vesc_2, lmax, rad, Psi, r2_ra2, ene, df, Qmin, invQbin);
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  return (0.5 * sum + hh * sub);
}
static inline double set_domain_boundary_velocity(const double hh, double * restrict tmin, double * restrict tmax, const double vesc_2, const double lmax, const double rad, const double Psi, const double r2_ra2, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  double fp = get_DEformula_velocity(tt, vesc_2, lmax, rad, Psi, r2_ra2, ene, df, Qmin, invQbin);
  const double f0 = fp;
  double sum = hh * fp;


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    const double ft = fp;

    tt += hh;    boundary = tt;
    fp = get_DEformula_velocity(tt, vesc_2, lmax, rad, Psi, r2_ra2, ene, df, Qmin, invQbin);

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
    fp = get_DEformula_velocity(tt, vesc_2, lmax, rad, Psi, r2_ra2, ene, df, Qmin, invQbin);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;

  return (sum);
}
static inline double integrate_DEformula_velocity(const double rad, const double logrmin, const double invlogrbin, real * restrict Psi_tot, const double lmax, const double ra2inv, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  const double r2_ra2 = rad * rad * ra2inv;
  const double Psi = getPsi(rad, logrmin, invlogrbin, Psi_tot);
  const double vesc_2 = 0.5 * sqrt(2.0 * Psi);

  double hh = 1.0;
  double tmin, tmax;
  double sum = set_domain_boundary_velocity(hh, &tmin, &tmax, vesc_2, lmax, rad, Psi, r2_ra2, ene, df, Qmin, invQbin);

  while( true ){
    const double f0 = sum;

    hh *= 0.5;
    sum = update_trapezoidal_velocity(hh, tmin, tmax, sum, vesc_2, lmax, rad, Psi, r2_ra2, ene, df, Qmin, invQbin);

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

  return (vesc_2 * sum);
}


static inline double get_DEformula_radial(const double tt, const double rr, const double drp_4, const double drm_4, const double logrmin, const double invlogrbin, real * restrict Psi_tot, const double lmax, const double ra2inv, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  const double sinh_t = M_PI_2 * sinh(tt);
  const double cosh_t = cosh(sinh_t);
  const double inv_cosh_t = 1.0 / cosh_t;

  const double rad = rr + inv_cosh_t * (drp_4 * exp(sinh_t) - drm_4 * exp(-sinh_t));

  return (cosh(tt) * inv_cosh_t * inv_cosh_t * rad * rad * integrate_DEformula_velocity(rad, logrmin, invlogrbin, Psi_tot, lmax, ra2inv, ene, df, Qmin, invQbin));
}
static inline double update_trapezoidal_radial(const double hh, const double tmin, const double tmax, const double sum, const double rr, const double drp_4, const double drm_4, const double logrmin, const double invlogrbin, real * restrict Psi_tot, const double lmax, const double ra2inv, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  /** initialization */
  double sub = 0.0;
  double tt = tmin + hh;

  /** employ mid-point rule */
  while( tt < tmax ){
    sub += get_DEformula_radial(tt, rr, drp_4, drm_4, logrmin, invlogrbin, Psi_tot, lmax, ra2inv, ene, df, Qmin, invQbin);
    tt += 2.0 * hh;
  }/* while( tt < tmax ){ */

  return (0.5 * sum + hh * sub);
}
static inline double set_domain_boundary_radial(const double hh, double * restrict tmin, double * restrict tmax, const double rr, const double drp_4, const double drm_4, const double logrmin, const double invlogrbin, real * restrict Psi_tot, const double lmax, const double ra2inv, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  const double converge = 1.0e-16;
  const double maximum = 128.0;

  double tt = 0.0;
  double fp = get_DEformula_radial(tt, rr, drp_4, drm_4, logrmin, invlogrbin, Psi_tot, lmax, ra2inv, ene, df, Qmin, invQbin);
  const double f0 = fp;
  double sum = hh * fp;


  /** determine upper boundary */
  double boundary = 0.0;
  double damp = 1.0;
  while( (damp > converge) && (boundary < maximum) ){
    const double ft = fp;

    tt += hh;    boundary = tt;
    fp = get_DEformula_radial(tt, rr, drp_4, drm_4, logrmin, invlogrbin, Psi_tot, lmax, ra2inv, ene, df, Qmin, invQbin);

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
    fp = get_DEformula_radial(tt, rr, drp_4, drm_4, logrmin, invlogrbin, Psi_tot, lmax, ra2inv, ene, df, Qmin, invQbin);

    sum += hh * fp;
    damp = fabs(ft) + fabs(fp);
  }/* while( (damp > converge) && (boundary > -maximum) ){ */
  *tmin = boundary;

  return (sum);
}
static inline double integrate_DEformula_radial(const double rr, const double drp_4, const double drm_4, const double logrmin, const double invlogrbin, real * restrict Psi_tot, const double lmax, const double ra2inv, real * restrict ene, real * restrict df, const double Qmin, const double invQbin)
{
  double hh = 1.0;
  double tmin, tmax;
  double sum = set_domain_boundary_radial(hh, &tmin, &tmax, rr, drp_4, drm_4, logrmin, invlogrbin, Psi_tot, lmax, ra2inv, ene, df, Qmin, invQbin);

  while( true ){
    const double f0 = sum;

    hh *= 0.5;
    sum = update_trapezoidal_radial(hh, tmin, tmax, sum, rr, drp_4, drm_4, logrmin, invlogrbin, Psi_tot, lmax, ra2inv, ene, df, Qmin, invQbin);

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

  return ((drp_4 + drm_4) * sum);
}


int main(int argc, char **argv)
{
  /** initialization */
  if( argc < 4 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 4);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -config=<char *>\n");
    __FPRINTF__(stderr, "          -BH_mass=<double>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 4 ){ */

  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,   "file", &file));
  char *fcfg;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "config", &fcfg));

  /** set unit system by reading the configuration file about physical parameters of the initial distribution */
  int unit, kind;
  profile_cfg *cfg;
  readProfileCfg(fcfg, &unit, &kind, &cfg);
  if( kind > NKIND_MAX ){    __KILL__(stderr, "ERROR: kind(= %d) must be smaller than %d\n", kind, NKIND_MAX);  }


#ifdef  USE_HDF5_FORMAT
  /** read DF */
  char filename[256];
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "df");
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  /** read attributes */
  hid_t attribute;
  int skind;
  attribute = H5Aopen(target, "kinds", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &skind));
  chkHDF5err(H5Aclose(attribute));
  int anisotropy;
  attribute = H5Aopen(target, "anisotropy", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &anisotropy));
  chkHDF5err(H5Aclose(attribute));
  int useDP;
  attribute = H5Aopen(target, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
#ifdef  USE_DOUBLE_PRECISION
  if( useDP != 1 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);
  }/* if( useDP != 1 ){ */
  const hid_t hdf5_real = H5T_NATIVE_DOUBLE;
#else///USE_DOUBLE_PRECISION
  if( useDP != 0 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);
  }/* if( useDP != 0 ){ */
  const hid_t hdf5_real = H5T_NATIVE_FLOAT;
#endif//USE_DOUBLE_PRECISION

  hid_t tmpgrp = H5Gopen(target, "series0", H5P_DEFAULT);
  int num;
  attribute = H5Aopen(tmpgrp, "num", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &num));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(tmpgrp));

  /** memory allocation */
  real **ene, *_ene;
  _ene = (real  *)malloc(sizeof(real) * skind * num);  if( _ene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _ene\n");  }
  ene  = (real **)malloc(sizeof(real) * skind      );  if(  ene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  ene\n");  }
  for(int ii = 0; ii < skind; ii++)
    ene[ii] = _ene + ii * num;
  real **df, *_df;
  _df = (real  *)malloc(sizeof(real) * skind * num);  if( _df == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _df\n");  }
  df  = (real **)malloc(sizeof(real) * skind      );  if(  df == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  df\n");  }
  for(int ii = 0; ii < skind; ii++)
    df[ii] = _df + ii * num;

  extern const double senergy_astro2com;
  for(int ii = 0; ii < skind; ii++){
    __NOTE__("read DF for %d-th component\n", ii);
    /** open an existing group */
    char grp[16];
    sprintf(grp, "series%d", ii);
    hid_t group = H5Gopen(target, grp, H5P_DEFAULT);

    /** read data */
    hid_t dataset;
    dataset = H5Dopen(group, "energy", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, ene[ii]));
    chkHDF5err(H5Dclose(dataset));
#pragma omp parallel for
    for(int jj = 0; jj < num; jj++)
      ene[ii][jj] = CAST_D2R(CAST_R2D(ene[ii][jj]) * senergy_astro2com);
    dataset = H5Dopen(group, "DF", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT,  df[ii]));
    chkHDF5err(H5Dclose(dataset));

    chkHDF5err(H5Gclose(group));
    __NOTE__("set DF for %d-th component\n", ii);
  }/* for(int ii = 0; ii < skind; ii++){ */

  /** close the file */
  chkHDF5err(H5Fclose(target));


  /** obtain radial profiles of Psi(r) and M(r) */
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "profile");
  target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  attribute = H5Aopen(target, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
#ifdef  USE_DOUBLE_PRECISION
  if( useDP != 1 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);
  }/* if( useDP != 1 ){ */
#else///USE_DOUBLE_PRECISION
  if( useDP != 0 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);
  }/* if( useDP != 0 ){ */
#endif//USE_DOUBLE_PRECISION

  tmpgrp = H5Gopen(target, "data0", H5P_DEFAULT);
  int nrad;
  attribute = H5Aopen(tmpgrp, "num", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &nrad));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(tmpgrp));

  /** memory allocation */
  real *rad;
  rad = (real *)malloc(sizeof(real) * nrad);  if( rad == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rad\n");  }
  real **enc, *_enc;
  _enc = (real  *)malloc(sizeof(real) * kind * nrad);  if( _enc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _enc\n");  }
  enc  = (real **)malloc(sizeof(real) * kind       );  if(  enc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  enc\n");  }
  for(int ii = 0; ii < kind; ii++)
    enc[ii] = _enc + ii * nrad;
  real **Psi, *_Psi;
  _Psi = (real  *)malloc(sizeof(real) * kind * nrad);  if( _Psi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _Psi\n");  }
  Psi  = (real **)malloc(sizeof(real) * kind       );  if(  Psi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  Psi\n");  }
  for(int ii = 0; ii < kind; ii++)
    Psi[ii] = _Psi + ii * nrad;
  real *enc_tot;
  enc_tot = (real *)malloc(sizeof(real) * nrad);  if( enc_tot == NULL ){    __KILL__(stderr, "ERROR: failure to allocate enc_tot\n");  }
  real *Psi_tot;
  Psi_tot = (real *)malloc(sizeof(real) * nrad);  if( Psi_tot == NULL ){    __KILL__(stderr, "ERROR: failure to allocate Psi_tot\n");  }
  real **rho, *_rho;
  _rho = (real  *)malloc(sizeof(real) * kind * nrad);  if( _rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _rho\n");  }
  rho  = (real **)malloc(sizeof(real) * kind       );  if(  rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  rho\n");  }
  for(int ii = 0; ii < kind; ii++)
    rho[ii] = _rho + ii * nrad;
  double *ra;
  ra = (double *)malloc(sizeof(double) * kind);  if( ra == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ra\n");  }
  double *ra2inv;
  ra2inv = (double *)malloc(sizeof(double) * kind);  if( ra2inv == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ra2inv\n");  }
  double *rs;  rs = (double *)malloc(sizeof(double) * kind);  if( rs == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rs\n");  }
  double *rc;  rc = (double *)malloc(sizeof(double) * kind);  if( rc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rc\n");  }


  for(int ii = 0; ii < kind; ii++){
    __NOTE__("read profile for %d-th component\n", ii);
    char grp[16];
    sprintf(grp, "data%d", ii);
    hid_t group = H5Gopen(target, grp, H5P_DEFAULT);

    /** read data */
    hid_t dataset;
    if( ii == 0 ){
      dataset = H5Dopen(group, "rad", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rad));
      chkHDF5err(H5Dclose(dataset));
    }/* if( ii == 0 ){ */
    dataset = H5Dopen(group, "enc", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, enc[ii]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dopen(group, "Psi", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, Psi[ii]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dopen(group, "rho", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho[ii]));
    chkHDF5err(H5Dclose(dataset));

    attribute = H5Aopen(group, "rs", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &rs[ii]));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Aopen(group, "rc", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &rc[ii]));
    chkHDF5err(H5Aclose(attribute));
    if( anisotropy ){
      attribute = H5Aopen(group, "ra", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &ra[ii]));
      chkHDF5err(H5Aclose(attribute));
      ra2inv[ii] = 1.0 / (DBL_MIN + ra[ii] * ra[ii]);
    }/* if( anisotropy ){ */
    else
      ra[ii] = ra2inv[ii] = 0.0;

    chkHDF5err(H5Gclose(group));
    __NOTE__("set profile for %d-th component\n", ii);
  }/* for(int ii = 0; ii < kind; ii++){ */


  /** close the file */
  chkHDF5err(H5Fclose(target));


  /** care unit system of rad, enc and Psi in HDF5 file (perhaps, it is in astrophysical unit system) */
  extern const double length_astro2com;
  extern const double mass_astro2com;
  extern const double density_astro2com;
#pragma omp parallel for
  for(int ii = 0; ii < nrad; ii++){
    rad[ii] = CAST_D2R(CAST_R2D(rad[ii]) * length_astro2com);

    double enc_tmp = 0.0;
    double Psi_tmp = 0.0;
    for(int jj = 0; jj < kind; jj++){
      enc[jj][ii] = CAST_D2R(CAST_R2D(enc[jj][ii]) *    mass_astro2com);      enc_tmp += CAST_R2D(enc[jj][ii]);
      Psi[jj][ii] = CAST_D2R(CAST_R2D(Psi[jj][ii]) * senergy_astro2com);      Psi_tmp += CAST_R2D(Psi[jj][ii]);
      rho[jj][ii] = CAST_D2R(CAST_R2D(rho[jj][ii]) * density_astro2com);
    }/* for(int jj = 0; jj < kind; jj++){ */

    enc_tot[ii] = CAST_D2R(enc_tmp);
    Psi_tot[ii] = CAST_D2R(Psi_tmp);
  }/* for(int ii = 0; ii < nrad; ii++){ */
  const double logrmin = log10(CAST_R2D(rad[0]));
  const double logrbin = (log10(CAST_R2D(rad[nrad - 1])) - logrmin) / (double)(nrad - 1);
  const double invlogrbin = 1.0 / logrbin;
  const double fac_drp = pow(10.0, logrbin) - 1.0;
  const double fac_drm = 1.0 - 1.0 / pow(10.0, logrbin);
  __NOTE__("nrad = %d, logrmin = %e, logrbin = %e, invlogrbin = %e\n", nrad, logrmin, logrbin, invlogrbin);
#if 0
  __FPRINTF__(stdout, "rad[1] - rad[0] = %e; %e, %e\n", CAST_R2D(rad[1]) - CAST_R2D(rad[0]), CAST_R2D(rad[0]) * (pow(10.0, logrbin) - 1.0), CAST_R2D(rad[1]) * (1.0 - 1.0 / pow(10.0, logrbin)));
  __FPRINTF__(stdout, "rad[1001] - rad[1000] = %e; %e, %e\n", CAST_R2D(rad[1001]) - CAST_R2D(rad[1000]), CAST_R2D(rad[1000]) * (pow(10.0, logrbin) - 1.0), CAST_R2D(rad[1001]) * (1.0 - 1.0 / pow(10.0, logrbin)));
  exit(0);
#endif


#else///USE_HDF5_FORMAT
  __KILL__(stderr, "ERROR: this program supports HDF5 format only.\n");
#endif//USE_HDF5_FORMAT


  /** set central BH and calculate maximum value of specific angular momentum */
  static double M_bh = 0.0;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "BH_mass", &M_bh));
  M_bh *= mass_astro2com;
  extern const real newton;
  extern const double lightspeed;
  const double lmax = 2.0 * sqrt(3.0) * CAST_R2D(newton) * M_bh / lightspeed;
  extern const double mass2astro;
  extern const double length2astro;
  extern const double velocity2astro;
  fprintf(stdout, "M_bh = %e Msun, lmax = %e kpc km/s\n", M_bh * mass2astro, lmax * length2astro * velocity2astro);


  /** integrate distribution function and evaluate accreted mass onto central MBH as a function of radius */
  real **loss, *_loss;
  _loss = (real  *)malloc(sizeof(real) * kind * nrad);  if( _loss == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _loss\n");  }
  loss  = (real **)malloc(sizeof(real) * kind       );  if(  loss == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  loss\n");  }
  for(int ii = 0; ii < skind; ii++)
    loss[ii] = _loss + ii * nrad;
  int *head;  head = (int *)malloc(sizeof(int) * kind);  if( head == NULL ){    __KILL__(stderr, "ERROR: failure to allocate head\n");  }
  int *tail;  tail = (int *)malloc(sizeof(int) * kind);  if( tail == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tail\n");  }

  const double factor = M_PI * M_PI * 2.0 * M_PI * M_PI * M_PI;
  for(int ii = 0; ii < skind; ii++){
    const double Qmin = CAST_R2D(ene[ii][0]);
    const double invQbin = (double)(num - 1) / (ene[ii][num - 1] - Qmin);

    const double rmin = rs[ii] * 9.765625e-4 * 0.5;
    const double rmax = rc[ii] * 1.5;
    int jhead = 0;
    for(int jj = 0; jj < nrad; jj++){
      const double rr = CAST_R2D(rad[jj]);
      if( (rr >= rmin) && (rr <= rmax) && (rho[ii][jj] > ZERO) ){
	jhead = jj;
	break;
      }
      loss[ii][jj] = ZERO;
    }/* for(int jj = 0; jj < nrad; jj++){ */
    int jtail = nrad - 1;
    for(int jj = jhead; jj < nrad; jj++){
      const double rr = CAST_R2D(rad[jj]);
      if( (rr < rmin) || (rr > rmax) || (rho[ii][jj] <= ZERO) ){
	jtail = jj - 1;
	break;
      }/* if( (rr < rmin) || (rr > rmax) || (rho[ii][jj] <= ZERO) ){ */
    }/* for(int jj = jhead; jj < nrad; jj++){ */
    for(int jj = jtail + 1; jj < nrad; jj++)
      loss[ii][jj] = ZERO;
    head[ii] = jhead;
    tail[ii] = jtail;

#pragma omp parallel for schedule(dynamic,4)
    for(int jj = jhead; jj < jtail + 1; jj++){
      const double rr = CAST_R2D(rad[jj]);
      /* calculate mass in loss cone */
      const double drp_4 = rr * fac_drp * 0.125;
      const double drm_4 = (jj > 0) ? (rr * fac_drm * 0.125) : (0.0);
      loss[ii][jj] = CAST_D2R(factor * integrate_DEformula_radial(rr, drp_4, drm_4, logrmin, invlogrbin, Psi_tot, lmax, ra2inv[ii], ene[ii], df[ii], Qmin, invQbin) * mass2astro);
    }/* for(int jj = jhead; jj < jtail + 1; jj++){ */
  }/* for(int ii = 0; ii < skind; ii++){ */


  /** evaluate timescale to generate loss cone (is free-fall time sufficient? because, angular momentum is small) */
  real *tff;
  tff = (real *)malloc(sizeof(real) * nrad);  if( tff == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tff\n");  }
  extern const double time2astro;
#pragma omp parallel for
  for(int ii = 0; ii < nrad; ii++)
    tff[ii] = CAST_D2R(M_PI_2 * rad[ii] * sqrt(rad[ii] / (2.0 * newton * enc_tot[ii])) * time2astro);


  /** evaluate total accreted mass and output loss[skind][nrad], Msum[skind][nrad], frac[skind][nrad], and tff[nrad] */
  real **Msum, *_Msum;
  _Msum = (real  *)malloc(sizeof(real) * kind * nrad);  if( _Msum == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _Msum\n");  }
  Msum  = (real **)malloc(sizeof(real) * kind       );  if(  Msum == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  Msum\n");  }
  for(int ii = 0; ii < kind; ii++)
    Msum[ii] = _Msum + ii * nrad;
  real **frac, *_frac;
  _frac = (real  *)malloc(sizeof(real) * kind * nrad);  if( _frac == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _frac\n");  }
  frac  = (real **)malloc(sizeof(real) * kind       );  if(  frac == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  frac\n");  }
  for(int ii = 0; ii < kind; ii++)
    frac[ii] = _frac + ii * nrad;
  extern const char   mass_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  extern const char length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  extern const char   time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];

  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "losscone");
  target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t attr_dims = 1;
  hid_t attrspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(target, "kind", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skind));
  chkHDF5err(H5Aclose(attribute));
  int minhead = head[0];
  int maxtail = tail[0];
  for(int ii = 1; ii < skind; ii++){
    if( minhead > head[ii] )      minhead = head[ii];
    if( maxtail < tail[ii] )      maxtail = tail[ii];
  }/* for(int ii = 1; ii < skind; ii++){ */
  const int ntot = maxtail + 1 - minhead;
  attribute = H5Acreate(target, "nrad", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &ntot));
  chkHDF5err(H5Aclose(attribute));

  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
  attribute = H5Acreate(target, "length_astro_unit_name", str4format, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, length_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "mass_astro_unit_name", str4format, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, mass_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "time_astro_unit_name", str4format, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, time_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));

  /* create the data space for the dataset */
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_FILE_COMPRESSION
  const hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  hsize_t cdims_loc;
#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  const uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  const uint szip_pixels_per_block = 8;
  hsize_t szip_cdims = 128 * szip_pixels_per_block;
#else///USE_SZIP_COMPRESSION
  /* compression using gzip */
  const uint gzip_compress_level = 9;/**< 9 is the maximum compression ratio */
  hsize_t gzip_cdims = 1024;
#endif//USE_SZIP_COMPRESSION
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#else///USE_FILE_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_FILE_COMPRESSION

  hsize_t dims = ntot;
  dataspace = H5Screate_simple(1, &dims, NULL);
#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
  property = H5Pcreate(H5P_DATASET_CREATE);
#ifdef  USE_SZIP_COMPRESSION
  szip_cdims = 128 * szip_pixels_per_block;
  cdims_loc = szip_cdims;
  if( dims < cdims_loc )
    cdims_loc = dims;
  if( dims > (hsize_t)szip_pixels_per_block ){
    if( cdims_loc > cdims_max )
      cdims_loc = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
  }/* if( dims > (hsize_t)szip_pixels_per_block ){ */
  else
    property = H5P_DEFAULT;
#else///USE_SZIP_COMPRESSION
  gzip_cdims = 1024;
  cdims_loc = gzip_cdims;
  if( dims < cdims_loc )
    cdims_loc = dims;
  if( cdims_loc > cdims_max )
    cdims_loc = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
  chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_SZIP_COMPRESSION
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)

  dataset = H5Dcreate(target, "rad", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rad[minhead]));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(target, "tff", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tff[minhead]));
  chkHDF5err(H5Dclose(dataset));
  for(int ii = 0; ii < skind; ii++){
    for(int jj = 0; jj < nrad; jj++)
      enc[ii][jj] = CAST_D2R(CAST_R2D(enc[ii][jj]) * mass2astro);

    char tag[16];
    sprintf(tag, "Menc_%d", ii);
    dataset = H5Dcreate(target, tag, hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &enc[ii][minhead]));
    chkHDF5err(H5Dclose(dataset));
  }/* for(int ii = 0; ii < skind; ii++){ */

#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#ifdef  USE_SZIP_COMPRESSION
  if( dims > (hsize_t)szip_pixels_per_block )
#endif//USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
  chkHDF5err(H5Sclose(dataspace));


  for(int ii = 0; ii < skind; ii++){
    const int jhead = head[ii];
    const int jtail = tail[ii];
    double Mloss = 0.0;
    double Min = CAST_R2D(enc[ii][jhead - 1]);
    for(int jj = jhead; jj < jtail + 1; jj++){
      const double Mtmp = CAST_R2D(loss[ii][jj]);
      Mloss += Mtmp;
      Msum[ii][jj] = CAST_D2R(Mloss);
      const double Mhere = CAST_R2D(enc[ii][jj]);
      frac[ii][jj] = Mtmp / (DBL_MIN + Mhere - Min);
      Min = Mhere;
    }/* for(int jj = jhead; jj < jtail + 1; jj++){ */

    const real Mhalf = CAST_D2R(0.5 * Mloss);
    real rhalf = ZERO;
    real thalf = ZERO;
    real Mdot = ZERO;
    for(int jj = jhead; jj < jtail + 1; jj++)
      if( Msum[ii][jj] > Mhalf ){
	rhalf = CAST_D2R(CAST_R2D(rad[jj]) * length2astro);
	thalf = tff[jj];
	Mdot = Mhalf / thalf;
	fprintf(stdout, "%d-th component: Mloss = %e %s, Mhalf = %e %s, rhalf = %e %s, thalf = %e %s; mean accretion rate = %e %s/%s\n", ii, Mloss, mass_astro_unit_name, Mhalf, mass_astro_unit_name, rhalf, length_astro_unit_name, thalf, time_astro_unit_name, Mdot, mass_astro_unit_name, time_astro_unit_name);
	break;
      }/* if( Msum[ii][jj] > Mhalf ){ */


    /** output */
    char grp[16];
    sprintf(grp, "data%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    const int jnum = jtail + 1 - jhead;
    dims = jnum;
    dataspace = H5Screate_simple(1, &dims, NULL);
#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
    property = H5Pcreate(H5P_DATASET_CREATE);
#ifdef  USE_SZIP_COMPRESSION
    szip_cdims = 128 * szip_pixels_per_block;
    cdims_loc = szip_cdims;
    if( dims < cdims_loc )
      cdims_loc = dims;
    if( dims > (hsize_t)szip_pixels_per_block ){
      if( cdims_loc > cdims_max )
	cdims_loc = cdims_max;
      chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( dims > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#else///USE_SZIP_COMPRESSION
    gzip_cdims = 1024;
    cdims_loc = gzip_cdims;
    if( dims < cdims_loc )
      cdims_loc = dims;
    if( cdims_loc > cdims_max )
      cdims_loc = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_SZIP_COMPRESSION
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)


    /** 1D (nrad) arrays */
    dataset = H5Dcreate(group, "rad", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rad[jhead]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "Mloc", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &loss[ii][jhead]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "Msum", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Msum[ii][jhead]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "frac", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &frac[ii][jhead]));
    chkHDF5err(H5Dclose(dataset));

    /** attributes */
    attribute = H5Acreate(group, "nrad", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &jnum));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "ra", H5T_NATIVE_DOUBLE, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &ra[ii]));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "Mloss", H5T_NATIVE_DOUBLE, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Mloss));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "Mhalf", hdf5_real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, hdf5_real, &Mhalf));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "rhalf", hdf5_real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, hdf5_real, &rhalf));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "thalf", hdf5_real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, hdf5_real, &thalf));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "Mdot", hdf5_real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, hdf5_real, &Mdot));
    chkHDF5err(H5Aclose(attribute));

#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#ifdef  USE_SZIP_COMPRESSION
    if( dims > (hsize_t)szip_pixels_per_block )
#endif//USE_SZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)

    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
  }/* for(int ii = 0; ii < skind; ii++){ */

  /** close the dataspace */
  chkHDF5err(H5Sclose(attrspace));

  /** close the file */
  chkHDF5err(H5Fclose(target));


  free(ene);  free(_ene);
  free(df);  free(_df);

  free(rad);
  free(enc);  free(_enc);
  free(Psi);  free(_Psi);
  free(rho);  free(_rho);
  free(enc_tot);
  free(Psi_tot);
  free(ra);  free(ra2inv);
  free(rs);  free(rc);

  free(loss);  free(_loss);
  free(head);  free(tail);

  free(tff);
  free(Msum);  free(_Msum);
  free(frac);  free(_frac);


  return (0);
}
