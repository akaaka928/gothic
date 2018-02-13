/**
 * @file external.c
 *
 * @brief Source code for setting an external fixed potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/02/13 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "macro.h"

#include "../misc/structure.h"

#include "profile.h"
#include "spline.h"

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
#include "potdens.h"
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK

#include "external.h"


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
/**
 * @fn genAdaptiveExtPotTbl1D
 *
 * @brief Generate table for cubic spline interpolation to represent external fixed potential field.
 *
 * @param (prf) radial profile of the component
 * @return (pot) potential field for cubic spline interpolation
 */
static inline void genAdaptiveExtPotTbl1D(profile *prf, potential_field *pot)
{
  __NOTE__("%s\n", "start");

  static double rr[N_EXT_POT_SPHE], yy[N_EXT_POT_SPHE];
  static double bp[N_EXT_POT_SPHE], y2[N_EXT_POT_SPHE];

  /** set sampling scale */
  const double dy = (-prf[NRADBIN - 1].psi + prf[0].psi) / (double)(N_EXT_POT_SPHE - 1);

  /** set innermost boundary */
  int num = 0;
  rr[num] =  prf[0].rad;
  double yold = -prf[0].psi;
  yy[num] = yold;
  num++;

  /** add sampling points */
  for(int kk = 1; kk < NRADBIN; kk++){
    double ynew = -prf[kk].psi;
    if( (ynew - yold) >= dy ){
      rr[num] = prf[kk].rad;
      yy[num] = ynew;
      yold = ynew;
      num++;
    }/* if( (ynew - yold) >= dy ){ */

    if( num >= N_EXT_POT_SPHE ){
      __KILL__(stderr, "ERROR: num (= %d) must be smaller than N_EXT_POT_SPHE (= %d)\n", num, N_EXT_POT_SPHE);
    }/* if( num >= N_EXT_POT_SPHE ){ */
  }/* for(int kk = 1; kk < NRADBIN; kk++){ */

  /** set outermost boundary */
  rr[num] =  prf[NRADBIN - 1].rad;
  yy[num] = -prf[NRADBIN - 1].psi;
  num++;

  /** apply cubic spline interpolation */
  pot->num = num;
  genCubicSpline1D(num, rr, yy, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, y2);

  for(int jj = 0; jj < num; jj++){
    pot->rad[jj]     = CAST_D2R(rr[jj]);
    pot->Phi[jj].val = CAST_D2R(yy[jj]);
    pot->Phi[jj].dr2 = CAST_D2R(y2[jj]);
  }/* for(int jj = 0; jj < num; jj++){ */

  __NOTE__("%s\n", "end");
}
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD


/**
 * @fn genExtPotTbl1D
 *
 * @brief Generate table(s) for cubic spline interpolation to represent external fixed potential field(s).
 *
 * @param (kind) number of components
 * @param (prf) radial profile of component(s)
 * @return (pot) potential field for cubic spline interpolation
 */
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
void genExtPotTbl1D(const int kind, profile **prf, potential_field *pot)
{
  __NOTE__("%s\n", "start");

  for(int ii = 0; ii < kind; ii++)
    genAdaptiveExtPotTbl1D(prf[ii], &pot[ii]);

  __NOTE__("%s\n", "end");
}
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
void genExtPotTbl1D(const int kind, profile **prf, potential_field *pot, const double invlogrbin)
{
  __NOTE__("%s\n", "start");

  const double cvt = M_LOG10E * 0.5 * invlogrbin;/**< log(e) / (2 log(rbin)) */

  for(int ii = 0; ii < kind; ii++){
    const int skip = NRADBIN / N_EXT_POT_SPHE;

    {
      /* const int jj = 0; */
      const double rad =  prf[ii][0].rad;
      const double Phi = -prf[ii][0].psi;
      const double Fr  = (prf[ii][1].psi - prf[ii][0].psi) * 2.0 * cvt / rad;

      pot[ii].Phi[0].Phi = CAST_D2R(Phi);
      pot[ii].Phi[0].Fr  = CAST_D2R(Fr);
    }

#pragma omp parallel for
    for(int jj = 1; jj < N_EXT_POT_SPHE - 1; jj++){
      const int kk = jj * skip;
      const double rad =  prf[ii][kk].rad;
      const double Phi = -prf[ii][kk].psi;
      const double Fr  = (prf[ii][kk + 1].psi - prf[ii][kk - 1].psi) * cvt / rad;

      pot[ii].Phi[jj].Phi = CAST_D2R(Phi);
      pot[ii].Phi[jj].Fr  = CAST_D2R(Fr);
    }/* for(int jj = 0; jj < N_EXT_POT_SPHE; jj++){ */


    {
      const int jj = N_EXT_POT_SPHE - 1;
      const int kk = jj * skip;
      const double rad =  prf[ii][kk].rad;
      const double Phi = -prf[ii][kk].psi;
      const double Fr  = ((kk + 1) <= (NRADBIN - 1)) ? ((prf[ii][kk + 1].psi - prf[ii][kk - 1].psi) * cvt / rad) : ((prf[ii][NRADBIN - 1].psi - prf[ii][NRADBIN - 2].psi) * 2.0 * cvt / rad);

      pot[ii].Phi[jj].Phi = CAST_D2R(Phi);
      pot[ii].Phi[jj].Fr  = CAST_D2R(Fr);
    }

    const double logrmin = log10(prf[ii][0].rad);
    const double logrmax = log10(prf[ii][(N_EXT_POT_SPHE - 1) * skip].rad);
    pot[ii].logrmin = CAST_D2R(logrmin);
    pot[ii].logrbin = CAST_D2R((logrmax - logrmin) / (double)(N_EXT_POT_SPHE - 1));
    pot[ii].num = N_EXT_POT_SPHE;
  }/* for(int ii = 0; ii < kind; ii++){ */


  __NOTE__("%s\n", "end");
}
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD


#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

/**
 * @fn genSuperposedPotFld1D
 *
 * @brief Generate superposed external, fixed potential-fields of multiple components.
 *
 * @param (kind) number of components
 * @param (skind) number of spherical symmetric components
 * @param (prf) radial profile of component(s)
 * @return (sphe) superposed potential field for cubic spline interpolation (only for spherical symmetric components)
 * @return (disk) superposed potential field for cubic spline interpolation (only for spherical averaged disk components)
 */
void genSuperposedPotFld1D(const int kind, const int skind, profile **prf, potential_field * restrict sphe, potential_field * restrict disk)
{
  __NOTE__("%s\n", "start");

  profile *tmp;
  tmp = (profile  *)malloc(sizeof(profile) * NRADBIN);
  if( tmp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp\n");  }

  /** generate superposed potential table for spherical component(s) */
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    tmp[ii].rad = prf[0][ii].rad;

    double psi = 0.0;
    for(int kk = 0; kk < skind; kk++)
      psi += prf[kk][ii].psi;

    tmp[ii].psi = psi;
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */
  genAdaptiveExtPotTbl1D(tmp, sphe);

  /** generate superposed potential table for spherical averaged disk component(s) */
  if( kind > skind ){
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      double psi = 0.0;
      for(int kk = skind; kk < kind; kk++)
	psi += prf[kk][ii].psi;

      tmp[ii].psi = psi;
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */
    genAdaptiveExtPotTbl1D(tmp, disk);
  }/* if( kind > skind ){ */

  free(tmp);

  __NOTE__("%s\n", "end");
}

#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

/**
 * @fn superposePotFld1D
 *
 * @brief Superpose external fixed potential fields of multiple components.
 *
 * @param (kind) number of components
 * @param (skind) number of spherical symmetric components
 * @param (pot) potential field for cubic spline interpolation
 * @return (sphe) superposed potential field for cubic spline interpolation (only for spherical symmetric components)
 * @return (disk) superposed potential field for cubic spline interpolation (only for spherical averaged disk components)
 */
void superposePotFld1D(const int kind, const int skind, potential_field * restrict pot, potential_field * restrict sphe, potential_field * restrict disk)
{
  __NOTE__("%s\n", "start");

  /** all potential field share the same radius */
  sphe->logrmin = pot[0].logrmin;
  sphe->logrbin = pot[0].logrbin;
  sphe->num     = pot[0].num;
#pragma omp parallel for
  for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){
    real Phi = ZERO;
    real Fr  = ZERO;
    for(int kk = 0; kk < skind; kk++){
      Phi += pot[kk].Phi[ii].Phi;
      Fr  += pot[kk].Phi[ii].Fr;
    }/* for(int kk = 0; kk < skind; kk++){ */
    sphe->Phi[ii].Phi = Phi;
    sphe->Phi[ii].Fr  = Fr;
  }/* for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){ */

  if( kind > skind ){
    disk->logrmin = pot[0].logrmin;
    disk->logrbin = pot[0].logrbin;
    disk->num     = pot[0].num;
#pragma omp parallel for
    for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){
      real Phi = ZERO;
      real Fr  = ZERO;
      for(int kk = skind; kk < kind; kk++){
	Phi += pot[kk].Phi[ii].Phi;
	Fr  += pot[kk].Phi[ii].Fr;
      }/* for(int kk = skind; kk < kind; kk++){ */
      disk->Phi[ii].Phi = Phi;
      disk->Phi[ii].Fr  = Fr;
    }/* for(int ii = 0; ii < N_EXT_POT_SPHE; ii++){ */
  }/* if( kind > skind ){ */

  __NOTE__("%s\n", "end");
}
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
/**
 * @fn extractDiskPotential
 *
 * @brief Extract potential field by disk components.
 *
 * @param (maxLev) maximum level of nested grid
 * @param (data) physical quantities of the disk component
 * @param (sphe) superposed potential field for cubic spline interpolation (only for spherical averaged disk components)
 * @return (disk) required set by GOTHIC
 */
void extractDiskPotential(const int maxLev, const disk_data data, const potential_field sphe, disk_potential *disk)
{
  __NOTE__("%s\n", "start");

  /** copy spherically averaged potential profile for long-range force */
  disk->sphe.rad = sphe.rad;
  disk->sphe.Phi = sphe.Phi;
  disk->sphe.num = sphe.num;

  disk->maxLev = maxLev;
  disk->NR = NDISKBIN_HOR;
  disk->Nz = NDISKBIN_VER;
  disk->hh = data.hh;

#pragma omp parallel for
  for(int ii = 0; ii < maxLev * NDISKBIN_HOR; ii++)
    disk->RR[ii] = CAST_D2R(data.hor[ii]);

#pragma omp parallel for
  for(int ii = 0; ii < maxLev * NDISKBIN_VER; ii++)
    disk->zz[ii] = CAST_D2R(data.ver[ii]);

  /* ii = 0 or jj = 0 are reserved for expressing symmetry over R = 0 or z = 0, respectively */
  for(int lev = 0; lev < maxLev; lev++){
    /* ii = 0; jj = 0 */
    disk->Phi[INDEX(maxLev, NDISKBIN_HOR + 1, NDISKBIN_VER + 1, lev, 0, 0)] = CAST_D2R(data.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 0, 0)]);

    /* ii = 0; jj >= 1 */
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
      disk->Phi[INDEX(maxLev, NDISKBIN_HOR + 1, NDISKBIN_VER + 1, lev, 0, 1 + jj)] = CAST_D2R(data.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 0, jj)]);

    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      disk->Phi[INDEX(maxLev, NDISKBIN_HOR + 1, NDISKBIN_VER + 1, lev, 1 + ii, 0)] = CAST_D2R(data.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)]);

      for(int jj = 0; jj < NDISKBIN_VER; jj++)
      disk->Phi[INDEX(maxLev, NDISKBIN_HOR + 1, NDISKBIN_VER + 1, lev, 1 + ii, 1 + jj)] = CAST_D2R(data.pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)]);
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  }/* for(int lev = 0; lev < maxLev; lev++){ */

  __NOTE__("%s\n", "end");
}

#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

/**
 * @fn extractDiskPotential
 *
 * @brief Extract potential field by disk components.
 *
 * @param (maxLev) maximum level of nested grid
 * @param (data) physical quantities of the disk component
 * @param (sphe) superposed potential field for cubic spline interpolation (only for spherical averaged disk components)
 * @return (disk) required set by GOTHIC
 */
void extractDiskPotential(const int maxLev, const int ndisk, disk_data *data, const potential_field sphe, disk_potential *disk)
{
  __NOTE__("%s\n", "start");


  /** copy spherically averaged potential profile for long-range force */
  disk->sphe.Phi = sphe.Phi;
  disk->sphe.num = sphe.num;
  disk->sphe.logrmin = sphe.logrmin;
  disk->sphe.logrbin = sphe.logrbin;

  disk->NR = NR_EXT_POT_DISK;
  disk->Nz = NZ_EXT_POT_DISK;

  double Rmin = DBL_MAX;
  double zmin = DBL_MAX;
  for(int ii = 0; ii < ndisk; ii++){
    Rmin = fmin(Rmin, data[ii].cfg->rs);
    zmin = fmin(zmin, data[ii].cfg->zd);
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  const int log2hmin = (int)floor(log2(EXT_POT_DISK_MIN_LENGTH * fmin(Rmin, zmin)));
  const real hh   = disk->hh   = LDEXP(UNITY,  log2hmin);
  const real hinv = disk->hinv = LDEXP(UNITY, -log2hmin);
  const int levTarget = (int)nearbyint(log2(data[0].hh * disk->hinv));

  if( ((hh * (double)NR_EXT_POT_DISK) >= data[0].Rmax) || ((hh * (double)NZ_EXT_POT_DISK) >= data[0].zmax) ){
    __KILL__(stderr, "ERROR: shrink domain size (Rmax = %e, zmax = %e) by reducing NZ_EXT_POT_DISK(%d) or EXT_POT_DISK_MIN_LENGTH(%e) to fit the size of numerical potential field (Rmax = %e, zmax = %e)\n", hh * (double)NR_EXT_POT_DISK, hh * (double)NZ_EXT_POT_DISK, NZ_EXT_POT_DISK, EXT_POT_DISK_MIN_LENGTH, data[0].Rmax, data[0].zmax);
  }/* if( ((hh * (double)NR_EXT_POT_DISK) >= data[0].Rmax) || ((hh * (double)NZ_EXT_POT_DISK) >= data[0].zmax) ){ */
  const real h0 = data[0].hh;
  const real h0inv = UNITY / h0;


  /** set gravitational potential field */
  /* ii = 0 (R-symmetry); jj = 0 (z-symmetry) */
  disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, 0, 0)] = CAST_D2R(data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, maxLev - 1, 0, 0)]);/**< (Phi[-0, -0] + Phi[-0, 0] + Phi[0, -0] + Phi[0, 0]) / 4 */

  /* ii = 0 (R-symmetry); jj >= 1 */
  for(int jj = 1; jj < NZ_EXT_POT_DISK + 1; jj++){
    const real zz = hh * (real)jj;
    const int levz = (int)FLOOR(LOG2(h0 * ((real)NDISKBIN_VER - HALF) / zz));
    const int jm = (int)FMAX(FLOOR(LDEXP(zz * h0inv, levz) - HALF), ZERO);
    if( levz >= levTarget ){
      disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, 0, jj)] =
	CAST_D2R(0.5 * (data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levz, 0, jm)] + data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levz, 0, jm + 1)]));/**< (Phi[-0, jm] + Phi[-0, jp] + Phi[0, jm] + Phi[0, jp]) / 4 */
    }/* if( levz >= levTarget ){ */
    else{
      const real dh = LDEXP(h0, -levz);
      const real dhinv = UNITY / dh;
      const real dPhidR = CAST_D2R(data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levz, 1, jm)] - data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levz, 0, jm)]) * HALF * dhinv;/**< Phi[-0, jm] = Phi[0, jm] */
      const real dPhidz = CAST_D2R(data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levz, 0, jm + 1)] - data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levz, 0, (jm > 0) ? (jm - 1) : 0)]) * HALF * dhinv;

      const real R0 = dh * HALF;
      const real z0 = dh * (HALF + (real)jm);

      disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, 0, jj)] = CAST_D2R(data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levz, 0, jm)]) - dPhidR * R0 + dPhidz * (zz - z0);/**< RR = 0 */
    }/* else{ */
  }/* for(int jj = 1; jj < NZ_EXT_POT_DISK + 1; jj++){ */

  /* ii >= 1 */
  for(int ii = 1; ii < NR_EXT_POT_DISK + 1; ii++){
    const real RR = hh * (real)ii;
    const int levR = (int)FLOOR(LOG2(h0 * ((real)NDISKBIN_HOR - HALF) / RR));

    /* jj = 0 (z-symmetry) */
    {
      const int im = (int)FMAX(FLOOR(LDEXP(RR * h0inv, levR) - HALF), ZERO);
      if( levR >= levTarget ){
	disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, 0)] =
	  CAST_D2R(0.5 * (data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, im, 0)] + data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, im + 1, 0)]));/**< (Phi[im, -0] + Phi[im, 0] + Phi[ip, -0] + Phi[ip, 0]) / 4 */
      }/* if( levR >= levTarget ){ */
      else{
	const real dh = LDEXP(h0, -levR);
	const real dhinv = UNITY / dh;
	const real dPhidR = CAST_D2R(data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, im + 1, 0)] - data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, (im > 0) ? (im - 1) : 0, 0)]) * HALF * dhinv;
	const real dPhidz = CAST_D2R(data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, im, 1)] - data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, im, 0)]) * HALF * dhinv;/**< Phi[im, -0] = Phi[im, 0] */

	const real R0 = dh * (HALF + (real)im);
	const real z0 = dh * HALF;

	disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, 0)] = CAST_D2R(data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, im, 0)]) + dPhidR * (RR - R0) - dPhidz * z0;/**< zz = 0 */
      }/* else{ */
    }

    /* jj >= 1 */
    for(int jj = 1; jj < NZ_EXT_POT_DISK + 1; jj++){
      const real zz = hh * (real)jj;
      const int levz = (int)FLOOR(LOG2(h0 * ((real)NDISKBIN_VER - HALF) / zz));

      int lev = (levz < levR) ? levz : levR;
      lev = (lev < (maxLev - 1)) ? lev : (maxLev - 1);

      const int im = (int)FMAX(FLOOR(LDEXP(RR * h0inv, lev) - HALF), ZERO);
      const int jm = (int)FMAX(FLOOR(LDEXP(zz * h0inv, lev) - HALF), ZERO);

      if( lev >= levTarget ){
	disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, jj)] =
	  CAST_D2R(0.25 * (data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, im, jm)] + data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, im, jm + 1)] + data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, im + 1, jm)] + data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, levR, im + 1, jm + 1)]));
      }/* if( lev >= levTarget ){ */
      else{
	const real dh = LDEXP(h0, -lev);
	const real dhinv = UNITY / dh;
	const real dPhidR = CAST_D2R(data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, im + 1, jm)] - data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, (im > 0) ? (im - 1) : 0, jm)]) * HALF * dhinv;
	const real dPhidz = CAST_D2R(data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, im, jm + 1)] - data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, im, (jm > 0) ? (jm - 1) : 0)]) * HALF * dhinv;

	const real R0 = dh * (HALF + (real)im);
	const real z0 = dh * (HALF + (real)jm);

	disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, jj)] = CAST_D2R(data[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, im, jm)]) + dPhidR * (RR - R0) + dPhidz * (zz - z0);
      }/* else{ */
    }/* for(int jj = 1; jj < NZ_EXT_POT_DISK + 1; jj++){ */
  }/* for(int ii = 1; ii < NR_EXT_POT_DISK + 1; ii++){ */


  /** set gravitational force field */
  /* ii = 0 (R-symmetry); jj = 0 (z-symmetry) */
  disk_grav FRz = {ZERO, ZERO};
  disk->FRz[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, 0, 0)] = FRz;

  /* ii = 0 (R-symmetry); jj >= 1 */
  for(int jj = 1; jj < NZ_EXT_POT_DISK; jj++){
    FRz.z = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, 0, jj + 1)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, 0, jj - 1)]) * HALF * hinv;
    disk->FRz[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, 0, jj)] = FRz;
  }/* for(int jj = 1; jj < NZ_EXT_POT_DISK; jj++){ */

  /* ii = 0 (R-symmetry); jj = Nz (z-edge) */
  FRz.z = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, 0, NZ_EXT_POT_DISK)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, 0, NZ_EXT_POT_DISK - 1)]) * hinv;
  disk->FRz[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, 0, NZ_EXT_POT_DISK)] = FRz;

  /* ii >= 1 */
  for(int ii = 1; ii < NR_EXT_POT_DISK; ii++){
    /* jj = 0 (z-symmetry) */
    FRz.R = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii + 1, 0)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii - 1, 0)]) * HALF * hinv;
    FRz.z = ZERO;
    disk->FRz[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, 0)] = FRz;

    /* jj >= 1 */
    for(int jj = 1; jj < NZ_EXT_POT_DISK; jj++){
      FRz.R = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii + 1, jj)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii - 1, jj)]) * HALF * hinv;
      FRz.z = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, jj + 1)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, jj - 1)]) * HALF * hinv;
      disk->FRz[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, jj)] = FRz;
    }/* for(int jj = 1; jj < NZ_EXT_POT_DISK; jj++){ */

    /* jj = Nz (z-edge) */
    FRz.R = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii + 1, NZ_EXT_POT_DISK)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii - 1, NZ_EXT_POT_DISK)]) * HALF * hinv;
    FRz.z = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, NZ_EXT_POT_DISK)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, NZ_EXT_POT_DISK - 1)]) * hinv;
    disk->FRz[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, ii, NZ_EXT_POT_DISK)] = FRz;
  }/* for(int ii = 1; ii < NR_EXT_POT_DISK; ii++){ */

  /* ii = NR (R-edge); jj = 0 (z-symmetry) */
  FRz.R = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK, 0)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK - 1, 0)]) * hinv;
  FRz.z = ZERO;
  disk->FRz[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK, 0)] = FRz;

  /* ii = NR (R-edge); jj >= 1 */
  for(int jj = 1; jj < NZ_EXT_POT_DISK; jj++){
    FRz.R = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK, jj)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK - 1, jj)]) * hinv;
    FRz.z = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK, jj + 1)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK, jj - 1)]) * HALF * hinv;
    disk->FRz[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK, jj)] = FRz;
  }/* for(int jj = 1; jj < NZ_EXT_POT_DISK; jj++){ */

  /* ii = NR (R-edge); jj = Nz (z-edge) */
  FRz.R = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK, NZ_EXT_POT_DISK)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK - 1, NZ_EXT_POT_DISK)]) * hinv;
  FRz.z = -(disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK, NZ_EXT_POT_DISK)] - disk->Phi[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK, NZ_EXT_POT_DISK - 1)]) * hinv;
  disk->FRz[INDEX2D(NR_EXT_POT_DISK + 1, NZ_EXT_POT_DISK + 1, NR_EXT_POT_DISK, NZ_EXT_POT_DISK)] = FRz;


  __NOTE__("%s\n", "end");
}
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK

#endif//SET_EXTERNAL_POTENTIAL_FIELD
