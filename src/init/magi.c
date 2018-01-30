/**
 * @file magi.c
 *
 * @brief Source code for MAGI (MAny-component Galaxy Initializer)
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/30 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

/** conversion from physical unit to computational unit must be performed internally */

/**
 * @def USE_SZIP_COMPRESSION
 *
 * @brief On to enable Szip compression for HDF5 files (default is OFF).
 */
/* #define USE_SZIP_COMPRESSION */

/**
 * @def RESET_ROTATION_AXIS
 *
 * @brief On to reset rotation axis to the z-axis (default is OFF).
 */
/* #define RESET_ROTATION_AXIS */


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "name.h"
#include "myutil.h"
#include "constants.h"
#include "timer.h"
#include "rotate.h"
#include "rand.h"

#ifdef  USE_SFMTJUMP
#include "SFMT-jump.h"
#include "sfmtjump_polynomial.h"
#include <omp.h>
#endif//USE_SFMTJUMP

#include "../misc/structure.h"
#include "../misc/allocate.h"

#ifndef RUN_WITHOUT_GOTHIC
#include "../misc/tune.h"
#include "../misc/brent.h"
#endif//RUN_WITHOUT_GOTHIC

#include "../file/io.h"

#include "magi.h"
#include "king.h"
#include "profile.h"
#include "eddington.h"
#include "table.h"
#include "abel.h"
#include "potdens.h"
#include "diskDF.h"

#ifdef  MAKE_COLUMN_DENSITY_PROFILE
#include "spline.h"
#endif//MAKE_COLUMN_DENSITY_PROFILE

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
#include "external.h"
#endif//SET_EXTERNAL_POTENTIAL_FIELD


/* global constants to set unit system, defined in constants.c */
extern const real newton;
extern const double     mass_astro2com,     mass2astro;extern const char        mass_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double   length_astro2com,   length2astro;extern const char      length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double     time_astro2com,     time2astro;extern const char        time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                     velocity2astro;extern const char    velocity_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                  col_density2astro;extern const char col_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                      density2astro;extern const char     density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                      senergy2astro;extern const char     senergy_astro_unit_name[CONSTANTS_H_CHAR_WORDS];


#include <gsl/gsl_integration.h>
double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];


/** measured by void initBenchmark_cpu(void) and void stopBenchmark_cpu(double *result) */
typedef struct
{
  double bodyAlloc, file;
  double spheAlloc, spheProf, spheDist, spheInfo, eddington;
  double diskAlloc, diskProf, diskDist, diskInfo, diskTbl, diskVel;
  double observe;
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
  double vdisp;
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  double column;
#endif//MAKE_COLUMN_DENSITY_PROFILE
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  double external;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
} breakdown;


/**
 * @fn isotropicDistribution
 *
 * @brief Generate isotropic distribution.
 *
 * @param (rad) norm of the vector
 * @return (vecx) x-component of the resultant vector
 * @return (vecy) y-component of the resultant vector
 * @return (vecz) z-component of the resultant vector
 * @param (rand) state of the pseudo random number generator
 *
 * @sa RANDVAL
 */
static inline void isotropicDistribution(const real rad, real *vecx, real *vecy, real *vecz, rand_state *rand)
{
  const real proj = RANDVAL(rand);
  *vecz = rad * proj;
  real Rproj = rad * SQRT(UNITY - proj * proj);

  real theta = TWO * CAST_D2R(M_PI) * UNIRAND(rand);
  *vecx = Rproj * COS(theta);
  *vecy = Rproj * SIN(theta);
}


#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#ifdef  __ICC
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#     pragma warning (disable:869)
#endif//__ICC

/**
 * @fn isInnerParticle4spherical
 *
 * @brief Detect particles locate in the inner regions.
 *
 * @param (x2) x squared
 * @param (y2) y squared
 * @param (z2) z squared
 * @param (rmax2) rmax squared
 * @param (zmax2) zmax squared
 * @return true when the particle locates in the inner region
 */
bool isInnerParticle4spherical(const double x2, const double y2, const double z2, const double rmax2, const double zmax2);
bool isInnerParticle4spherical(const double x2, const double y2, const double z2, const double rmax2, const double zmax2){  return ((x2 + y2 + z2) < rmax2);}

/**
 * @fn isInnerParticle4disk
 *
 * @brief Detect particles locate in the inner regions.
 *
 * @param (x2) x squared
 * @param (y2) y squared
 * @param (z2) z squared
 * @param (Rmax2) Rmax squared
 * @param (zmax2) zmax squared
 * @return true when the particle locates in the inner region
 */
bool isInnerParticle4disk     (const double x2, const double y2, const double z2, const double Rmax2, const double zmax2);
bool isInnerParticle4disk     (const double x2, const double y2, const double z2, const double Rmax2, const double zmax2){  return (((x2 + y2) < Rmax2) && (z2 < zmax2));}

#ifdef  __ICC
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#     pragma warning ( enable:869)
#endif//__ICC
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)


/**
 * @fn pickRetrogradingParticles
 *
 * @brief Select retrograding particles.
 *
 * @param (num) total number of particles in the component
 * @param (head) index of the head particle in the component
 * @param (frac) retrograding fraction
 * @return (retroNum) number of retrograding particles
 * @return (retroHead) head index of retrograding particles
 */
static inline void pickRetrogradingParticles(const ulong num, const ulong head, const double frac, ulong *retroNum, ulong *retroHead)
{
  *retroNum  = (ulong)floor((double)num * frac);
  *retroHead = head + (num - (*retroNum));
}
/**
 * @fn retrograder
 *
 * @brief Flip rotating velocity of particles.
 *
 * @param (num) total number of particles in the component
 * @param (head) index of the head particle in the component
 * @return (body) physical quantities of N-body particles
 * @param (frac) retrograding fraction
 *
 * @sa pickRetrogradingParticles
 */
void retrograder(const ulong num, const ulong head, iparticle body, const double frac);
void retrograder(const ulong num, const ulong head, iparticle body, const double frac)
{
  __NOTE__("%s\n", "start");

  /** pick up retrograding particles */
  ulong retroNum, retroHead;
  pickRetrogradingParticles(num, head, frac, &retroNum, &retroHead);

  /** flip rotating velocity */
#ifndef USE_SFMTJUMP
#pragma omp parallel
#endif//USE_SFMTJUMP
#pragma omp for
  for(ulong ii = retroHead; ii < retroHead + retroNum; ii++){
#ifdef  BLOCK_TIME_STEP
    body.vel[ii].x *= -UNITY;
    body.vel[ii].y *= -UNITY;
#else///BLOCK_TIME_STEP
    body.vx[ii]    *= -UNITY;
    body.vy[ii]    *= -UNITY;
#endif//BLOCK_TIME_STEP
  }/* for(ulong ii = retroHead; ii < retroHead + retroNum; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn shiftCenter
 *
 * @brief Shift center to set the center-of-mass at the coordinate origin with null bulk velocity.
 *
 * @param (num) total number of particles in the component
 * @param (head) index of the head particle in the component
 * @return (body) physical quantities of N-body particles
 * @param (cfg) physical properties of the component
 *
 * @sa isInnerParticle4spherical
 * @sa isInnerParticle4disk
 * @sa initRotationMatrices
 * @sa rotateVector
 */
void shiftCenter(const ulong num, const ulong head, iparticle body, const profile_cfg cfg);
void shiftCenter(const ulong num, const ulong head, iparticle body, const profile_cfg cfg)
{
  __NOTE__("%s\n", "start");


#ifndef USE_SFMTJUMP
#pragma omp parallel
#endif//USE_SFMTJUMP
  {
    /** calculate center-of-mass and bulk-motion */
    double comx_loc = 0.0, comy_loc = 0.0, comz_loc = 0.0;
    double velx_loc = 0.0, vely_loc = 0.0, velz_loc = 0.0;
    double Mtot_loc = 0.0;
    static double comx, comy, comz, velx, vely, velz, Mtot;
#pragma omp single
    comx = comy = comz = velx = vely = velz = Mtot = 0.0;

    /** use particles (r < 3 rs            ) for spherical components */
    /** use particles (R < 3 Rd, |z| < 3 zd) for      disk components */
    const double rmax2 = 9.0 * cfg.rs * cfg.rs;
    const double zmax2 = 9.0 * ((cfg.kind < 0) ? (cfg.zd * cfg.zd) : (cfg.rs * cfg.rs));
    bool (*isInnerParticle)(double, double, double, double, double) = (cfg.kind < 0) ? isInnerParticle4disk : isInnerParticle4spherical;

#pragma omp for
    for(ulong ii = head; ii < head + num; ii++){
      const double xx = CAST_R2D(body.pos[ii].x);
      const double yy = CAST_R2D(body.pos[ii].y);
      const double zz = CAST_R2D(body.pos[ii].z);
      if( isInnerParticle(xx * xx, yy * yy, zz * zz, rmax2, zmax2) ){
	const double mass = CAST_R2D(body.pos[ii].m);
	Mtot_loc += mass;

	comx_loc += mass * xx;
	comy_loc += mass * yy;
	comz_loc += mass * zz;
#ifdef  BLOCK_TIME_STEP
	velx_loc += mass * CAST_R2D(body.vel[ii].x);
	vely_loc += mass * CAST_R2D(body.vel[ii].y);
	velz_loc += mass * CAST_R2D(body.vel[ii].z);
#else///BLOCK_TIME_STEP
	velx_loc += mass * CAST_R2D(body.vx[ii]);
	vely_loc += mass * CAST_R2D(body.vy[ii]);
	velz_loc += mass * CAST_R2D(body.vz[ii]);
#endif//BLOCK_TIME_STEP
      }/* if( isInnerParticle(xx * xx, yy * yy, zz * zz, rmax2, zmax2) ){ */
    }/* for(ulong ii = head; ii < head + num; ii++){ */

#pragma omp atomic
    Mtot += Mtot_loc;
#pragma omp atomic
    comx += comx_loc;
#pragma omp atomic
    comy += comy_loc;
#pragma omp atomic
    comz += comz_loc;
#pragma omp atomic
    velx += velx_loc;
#pragma omp atomic
    vely += vely_loc;
#pragma omp atomic
    velz += velz_loc;
#pragma omp barrier
#pragma omp single
    {
      double Minv = 1.0 / (DBL_MIN + Mtot);
      comx *= Minv;      comy *= Minv;      comz *= Minv;
      velx *= Minv;      vely *= Minv;      velz *= Minv;
    }

#ifdef  PROGRESS_REPORT_ON
#pragma omp single nowait
    {
      fprintf(stdout, "# center-of-mass shift: %e, %e, %e\n", comx, comy, comz);
      fprintf(stdout, "#    bulk motion shift: %e, %e, %e\n", velx, vely, velz);
      fflush(stdout);
    }
#endif//PROGRESS_REPORT_ON


    /** shift the coordinate system to the center-of-mass rest frame */
#pragma omp for
    for(ulong ii = head; ii < head + num; ii++){
      body.pos[ii].x = CAST_D2R(CAST_R2D(body.pos[ii].x) - comx);
      body.pos[ii].y = CAST_D2R(CAST_R2D(body.pos[ii].y) - comy);
      body.pos[ii].z = CAST_D2R(CAST_R2D(body.pos[ii].z) - comz);
#ifdef  BLOCK_TIME_STEP
      body.time[ii].t0 = body.time[ii].t1 = 0.0;
      body.vel[ii].x = CAST_D2R(CAST_R2D(body.vel[ii].x) - velx);
      body.vel[ii].y = CAST_D2R(CAST_R2D(body.vel[ii].y) - vely);
      body.vel[ii].z = CAST_D2R(CAST_R2D(body.vel[ii].z) - velz);
      body.vel[ii].dt = ZERO;
#else///BLOCK_TIME_STEP
      body.vx[ii] = CAST_D2R(CAST_R2D(body.vx[ii]) - velx);
      body.vy[ii] = CAST_D2R(CAST_R2D(body.vy[ii]) - vely);
      body.vz[ii] = CAST_D2R(CAST_R2D(body.vz[ii]) - velz);
#endif//BLOCK_TIME_STEP
    }/* for(ulong ii = head; ii < head + num; ii++){ */


#ifdef  RESET_ROTATION_AXIS
    /** calculate angular momentum vector */
    double Lx_loc = 0.0, Ly_loc = 0.0, Lz_loc = 0.0;
    static double Lx, Ly, Lz;
#pragma omp single
    Lx = Ly = Lz = 0.0;
#pragma omp for
    for(ulong ii = head; ii < head + num; ii++){
      const double rx = CAST_R2D(body.pos[ii].x);
      const double ry = CAST_R2D(body.pos[ii].y);
      const double rz = CAST_R2D(body.pos[ii].z);
      if( isInnerParticle(rx * rx, ry * ry, rz * rz, rmax2, zmax2) ){
	const double mass = CAST_R2D(body.pos[ii].m);
#ifdef  BLOCK_TIME_STEP
	const double px = CAST_R2D(body.vel[ii].x) * mass;
	const double py = CAST_R2D(body.vel[ii].y) * mass;
	const double pz = CAST_R2D(body.vel[ii].z) * mass;
#else///BLOCK_TIME_STEP
	const double px = CAST_R2D(body.vx[ii]) * mass;
	const double py = CAST_R2D(body.vy[ii]) * mass;
	const double pz = CAST_R2D(body.vz[ii]) * mass;
#endif//BLOCK_TIME_STEP

	Lx_loc += ry * pz - rz * py;
	Ly_loc += rz * px - rx * pz;
	Lz_loc += rx * py - ry * px;
      }/* if( isInnerParticle(rx * rx, ry * ry, rz * rz, rmax2, zmax2) ){ */
    }/* for(ulong ii = head; ii < head + num; ii++){ */

    /** rotate galaxy (if necessary) */
#pragma omp atomic
    Lx += Lx_loc;
#pragma omp atomic
    Ly += Ly_loc;
#pragma omp atomic
    Lz += Lz_loc;
#pragma omp barrier
    const double L2 = Lx * Lx + Ly * Ly + Lz * Lz;
    if( L2 > 1.0e-6 ){
      real ini[3] = {CAST_D2R(Lx), CAST_D2R(Ly), CAST_D2R(Lz)};
      real fin[3] = {ZERO, ZERO, UNITY};

      real rot[3][3], inv[3][3];
      initRotationMatrices(ini, fin, rot, inv);

#ifndef USE_SFMTJUMP
#pragma omp parallel
#endif//USE_SFMTJUMP
#pragma omp for
      for(ulong ii = head; ii < head + num; ii++){
	real bfr[3], aft[3];
	/** rotate position */
	bfr[0] = body.pos[ii].x;
	bfr[1] = body.pos[ii].y;
	bfr[2] = body.pos[ii].z;
	rotateVector(bfr, rot, aft);
	body.pos[ii].x = aft[0];
	body.pos[ii].y = aft[1];
	body.pos[ii].z = aft[2];

	/** rotate velocity */
#ifdef  BLOCK_TIME_STEP
	bfr[0] = body.vel[ii].x;
	bfr[1] = body.vel[ii].y;
	bfr[2] = body.vel[ii].z;
#else///BLOCK_TIME_STEP
	bfr[0] = body.vx[ii];
	bfr[1] = body.vy[ii];
	bfr[2] = body.vz[ii];
#endif//BLOCK_TIME_STEP
	rotateVector(bfr, rot, aft);
#ifdef  BLOCK_TIME_STEP
	body.vel[ii].x = aft[0];
	body.vel[ii].y = aft[1];
	body.vel[ii].z = aft[2];
#else///BLOCK_TIME_STEP
	body.vx[ii] = aft[0];
	body.vy[ii] = aft[1];
	body.vz[ii] = aft[2];
#endif//BLOCK_TIME_STEP
      }/* for(ulong ii = head; ii < head + num; ii++){ */
    }/* if( L2 > 1.0e-6 ){ */
#endif//RESET_ROTATION_AXIS
  }


  __NOTE__("%s\n", "end");
}


/**
 * @fn distributeSpheroidParticles
 *
 * @brief Distribute N-body particles in the spherical symmetric component.
 *
 * @return (Nuse) total number of N-body particles distributed
 * @return (body) N-body particles
 * @param (mass) mass of N-body particles
 * @param (cfg) physical properties of the component
 * @param (prf) radial profile of the component
 * @param (df) distribution function of the component
 * @param (rand) state of pseudo random numbers
 * @return cutoff energy corresponding to the escape velocity
 *
 * @sa UNIRAND_DBL
 * @sa isotropicDistribution
 * @sa getDF
 */
double distributeSpheroidParticles(ulong *Nuse, iparticle body, const real mass, profile_cfg cfg, profile *prf, dist_func *df, rand_state *rand);
double distributeSpheroidParticles(ulong *Nuse, iparticle body, const real mass, profile_cfg cfg, profile *prf, dist_func *df, rand_state *rand)
{
  __NOTE__("%s\n", "start");


  const ulong num = cfg.num;

#ifdef  PROGRESS_REPORT_ON
#ifdef  USE_SFMTJUMP
  const ulong nunit = (ulong)floorf(0.1f * (float)num / (float)omp_get_num_threads());
#else///USE_SFMTJUMP
  const ulong nunit = (ulong)floorf(0.1f * (float)num);
#endif//USE_SFMTJUMP
  ulong stage = 1;
  ulong Npart = 0;
#ifdef  USE_SFMTJUMP
#pragma omp single nowait
#endif//USE_SFMTJUMP
  {
    fprintf(stdout, "#\n#\n# start distributing spherical component particles (%zu bodies: [%zu:%zu])\n", num, *Nuse, (*Nuse) + num - 1);
    fflush(stdout);
  }
#endif//PROGRESS_REPORT_ON

  const double    Emin = df[          0].ene;
  const double    Emax = df[NENEBIN - 1].ene;
  const double invEbin = (double)(NENEBIN - 1) / (Emax - Emin);

#if 1
  const int iout = cfg.iout;
  const double Ecut = prf[iout].psi_tot;
#else
  /* no more required or something wrong... */
  double fmax = 0.0;
  int   iout = 0;
  for(int ii = 0; ii < NRADBIN; ii++){
    double floc = prf[ii].rad * prf[ii].rad * prf[ii].rho;
    if(                            floc > fmax ){      fmax = floc;    }
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */
  const double Ecut = (iout != 0) ? prf[iout].psi_tot : Emin;
#endif

  const double Mmax = cfg.Mtot;
  const double Mmin = prf[0].enc;

#if 0
#pragma omp single
  {
    fprintf(stdout, "Mmin = %e, Mmax = %e; Emin = %e, Emax = %e; Ecut = %e, iout = %d\n", Mmin, Mmax, Emin, Emax, Ecut, iout);
    fflush(NULL);
    exit(0);
  }
#endif

#ifdef  USE_SFMTJUMP
#pragma omp for
#endif//USE_SFMTJUMP
  for(ulong ii = *Nuse; ii < *Nuse + num; ii++){
    /** set spatial distribution by table search */
    const double tmp = Mmin + (Mmax - Mmin) * UNIRAND_DBL(rand);
    int ll = 0;
    /* int rr = NRADBIN - 1; */
    int rr = iout;
    while( true ){
      uint cc = (ll + rr) >> 1;

      if( (prf[cc].enc - tmp) * (prf[ll].enc - tmp) <= 0.0 )	rr = cc;
      else	                                                ll = cc;

      if( (ll + 1) == rr )	break;
    }/* while( true ){ */

    const double alpha = (tmp - prf[ll].enc) / (prf[rr].enc - prf[ll].enc);
    const double rad = (1.0 - alpha) * prf[ll].rad + alpha * prf[rr].rad;
    isotropicDistribution(CAST_D2R(rad), &(body.pos[ii].x), &(body.pos[ii].y), &(body.pos[ii].z), rand);

#if 0
    fprintf(stderr, "ii = %zu, rad = %e\n", ii, rad);
    fflush(NULL);
#endif

    /** determine velocity distribution by rejection method */
    const double psi = (1.0 - alpha) * prf[ll].psi_tot + alpha * prf[rr].psi_tot;
    const double vesc = sqrt(2.0 * (psi - Ecut));

#if 0
    if( fpclassify(vesc) != FP_NORMAL ){
      fprintf(stderr, "rad = %e, psi = %e, Ecut = %e, Emin = %e, rout = %e, vesc = %e\n", rad, psi, Ecut, Emin, cfg.rmax, vesc);
      fflush(NULL);
      exit(0);
    }
#endif

    const double v2Fmax = vesc * vesc * getDF(psi, df, Emin, invEbin);
    double vel;
    while( true ){
      vel = vesc * UNIRAND_DBL(rand);

      const double ene = psi - 0.5 * vel * vel;
      const double val = vel * vel * getDF(ene, df, Emin, invEbin);
      const double try = v2Fmax * UNIRAND_DBL(rand);

#if 0
      if( ii == 655364 ){
	fprintf(stderr, "vesc = %e, vel = %e, rad = %e, psi = %e, ene = %e, v2Fmax = %e, val = %e, try = %e\n", vesc, vel, rad, psi, ene, v2Fmax, val, try);
	fflush(stderr);
      }
#endif

      if( val > try )	break;
    }/* while( true ){ */

#if 0
    fprintf(stderr, "ii = %zu, vel = %e\n", ii, vel);
    fflush(NULL);
#endif


#ifdef  BLOCK_TIME_STEP
    isotropicDistribution(CAST_D2R(vel), &(body.vel[ii].x), &(body.vel[ii].y), &(body.vel[ii].z), rand);
#else///BLOCK_TIME_STEP
    isotropicDistribution(CAST_D2R(vel), &(body.vx[ii]), &(body.vy[ii]), &(body.vz[ii]), rand);
#endif//BLOCK_TIME_STEP

    body.acc[ii].x   = body.acc[ii].y = body.acc[ii].z = ZERO;
    body.pos[ii].m   = mass;
    body.acc[ii].pot = -CAST_D2R(psi);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
    body.acc_ext[ii].x = body.acc_ext[ii].y = body.acc_ext[ii].z = body.acc_ext[ii].pot = ZERO;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    body.idx[ii] = ii;

#ifdef  PROGRESS_REPORT_ON
    Npart++;
    if( Npart == (stage * nunit) ){
#ifdef  USE_SFMTJUMP
      fprintf(stdout, "# ~%zu%% on thread %d\n", stage * 10, omp_get_thread_num());
#else///USE_SFMTJUMP
      fprintf(stdout, "# ~%zu%%\n", stage * 10);
#endif//USE_SFMTJUMP
      fflush(stdout);
      stage++;
    }/* if( Npart == (stage * nunit) ){ */
#endif//PROGRESS_REPORT_ON
  }/* for(ulong ii = *Nuse; ii < *Nuse + num; ii++){ */

  *Nuse += num;

#ifdef  PROGRESS_REPORT_ON
#ifdef  USE_SFMTJUMP
#pragma omp barrier
#pragma omp master
#endif//USE_SFMTJUMP
  {
    fprintf(stdout, "# finish distributing spherical component particles (%zu bodies)\n", num);
    fflush(stdout);
  }
#endif//PROGRESS_REPORT_ON


  __NOTE__("%s\n", "end");
  return (Ecut);
}


void outputFundamentalInformation
(const int unit, const int kind, const int skind, profile_cfg *cfg, profile **prf, const ulong Ntot, const real eps, const double snapshotInterval, const double ft, char file[]);
void outputRadialProfiles
(const int kind, profile **prf,
#ifdef  USE_HDF5_FORMAT
 profile_cfg *cfg, const real eps,
#endif//USE_HDF5_FORMAT
 char file[]);

void outputDistributionFunction(const int skind, dist_func **df, char file[]);

#ifdef  USE_HDF5_FORMAT
void outputRepresentativeQuantities
(const int kind, profile_cfg *cfg, const int skind, profile **prf, const int ndisk, const int maxLev, disk_data * disk, char file[]);
#endif//USE_HDF5_FORMAT

void evaluateObservables(const int kind, const int skind, profile_cfg *cfg, profile **prf);
void writeDiskData(char *file, const int ndisk, const int maxLev, disk_data *disk);


int main(int argc, char **argv)
{
  /** read input arguments */
  if( argc < 9 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 9);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -Ntot=<unsigned long int>\n");
    __FPRINTF__(stderr, "          -config=<char *>\n");
    __FPRINTF__(stderr, "          -eps=<real> -eta=<real>\n");
    __FPRINTF__(stderr, "          -ft=<real> -snapshotInterval=<real> -saveInterval=<real>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 9 ){ */

  /** read input arguments do not depend on the unit system adopted in the numerical simulation */
  char *file;  requiredCmdArg(getCmdArgStr( argc, (const char * const *)argv,   "file", &file));
  char *fcfg;  requiredCmdArg(getCmdArgStr( argc, (const char * const *)argv, "config", &fcfg));
  ulong Ntot;  requiredCmdArg(getCmdArgUlg( argc, (const char * const *)argv,   "Ntot", &Ntot));
  real   eta;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,    "eta", &eta));
  double saveInterval;  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "saveInterval", &saveInterval));

  /** set unit system by reading the configuration file about physical parameters of the initial distribution */
  int unit, kind;
  profile_cfg *cfg;
  readProfileCfg(fcfg, &unit, &kind, &cfg);
  if( kind > NKIND_MAX ){    __KILL__(stderr, "ERROR: kind(= %d) must be smaller than %d\n", kind, NKIND_MAX);  }

  /** read input arguments depend on the unit system adopted in the numerical simulation */
  double tmp;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "snapshotInterval", &tmp));
  double snapshotInterval = ldexp(1.0, (int)floor(log2(tmp * time_astro2com)));
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "eps", &tmp));  real   eps = CAST_D2R(tmp * length_astro2com);
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "ft",  &tmp));  double  ft =         (tmp *   time_astro2com);

  static breakdown execTime;


  /** set number of particles */
  ulong Nrem = 0;
  for(int ii = 0; ii < kind; ii++)
    Nrem += (cfg[ii].forceNum == 1) ? cfg[ii].num : 0;
  if( Nrem > Ntot ){
    __KILL__(stderr, "ERROR: the sum of number of particles for each component (%zu) exceeds the specified total number of particles (%zu).\n", Nrem, Ntot);
  }/* if( Nrem > Ntot ){ */
  Nrem = Ntot - Nrem;

  /** set number of particles to represent each profile */
  double Mtot = 0.0;
  for(int ii = 0; ii < kind; ii++)
    if( cfg[ii].forceNum != 1 )
      Mtot += cfg[ii].Mtot;
  const double Minv = 1.0 / Mtot;
  ulong Nuse = 0;
  int kidx = kind;
  double Mmax = 0.0;
  for(int ii = 0; ii < kind; ii++){
    if( cfg[ii].forceNum != 1 ){
      /** number of particles is determined by mass fraction */
      cfg[ii].num = (ulong)(cfg[ii].Mtot * Minv * (double)Nrem);
      Nuse += cfg[ii].num;

      if( cfg[ii].Mtot > Mmax ){
	kidx = ii;
	Mmax = cfg[ii].Mtot;
      }/* if( cfg[ii].Mtot > Mmax ){ */
    }/* if( cfg[ii].forceNum != 1 ){ */
  }/* for(int ii = 0; ii < kind; ii++){ */

  if( (kidx == kind) && (Nuse != Nrem) ){
    __KILL__(stderr, "ERROR: mismatch about number of particles detected (Nuse = %zu, Nrem = %zu) with %d components\n", Nuse, Nrem, kind);
  }/* if( (kidx == kind) && (Nuse != Nrem) ){ */
  if( Nuse != Nrem )
    cfg[kidx].num += (Nrem - Nuse);


  /** initialize the table for Gaussian Quadrature provided by GSL */
  for(int ii = 0; ii < NTBL_GAUSS_QD; ii++){
    gsl_gaussQD_pos   [ii] = 0.0;
    gsl_gaussQD_weight[ii] = 0.0;
  }/* for(int ii = 0; ii < NTBL_GAUSS_QD; ii++){ */
  gsl_integration_glfixed_table *tab;
  tab = gsl_integration_glfixed_table_alloc(NINTBIN);
  int max = (NINTBIN >> 1) + (NINTBIN & 1);
  for(int ii = 0; ii < max; ii++){
    gsl_gaussQD_pos   [ii] = (*tab).x[(max - 1) - ii];
    gsl_gaussQD_weight[ii] = (*tab).w[(max - 1) - ii];
  }/* for(int ii = 0; ii < max; ii++){ */
  gsl_integration_glfixed_table_free(tab);


  /** set global settings */
  /** identifier for disk components is negative value (-1, -2) */
  /** identifier for spherical components is positive value (0, 1, 2, ...) */
  int skind = kind;
  for(int ii = 0; ii < kind; ii++)
    if( cfg[ii].kind < 0 )
      skind--;
  const int ndisk = kind - skind;
  const bool addDisk = (ndisk != 0) ? true : false;
  if( addDisk )
    for(int ii = skind; ii < kind; ii++)
      if( cfg[ii].kind >= 0 ){      	__KILL__(stderr, "ERROR: disk component must be last component(s).\n\tModify \"%s/%s\".\n", CFGFOLDER, fcfg);      }

  double rmax = 0.0;
#if 1
  bool cutoff = false;
  for(int ii = 0; ii < kind; ii++){
    cutoff |= cfg[ii].cutoff;
    if( cfg[ii].cutoff ){      if( rmax < cfg[ii].rc )	rmax = cfg[ii].rc;    }
    else{                      if( rmax < cfg[ii].rs )	rmax = cfg[ii].rs;    }
  }/* for(int ii = 0; ii < kind; ii++){ */

  if( cutoff )    rmax *= (double)(NINTBIN * 10);/**< sufficiently greater than outermost radius times NINTBIN */
  else            rmax *= 1000.0;
#else
  for(int ii = 0; ii < kind; ii++){
    if( cfg[ii].cutoff )
      rmax = fmax(rmax,  10.0 * cfg[ii].rc);
    else
      rmax = fmax(rmax, 100.0 * cfg[ii].rs);
  }/* for(int ii = 0; ii < kind; ii++){ */
#endif

  const double logrmin = log10(MINRAD);
  const double logrmax = log10(rmax);
  const double logrbin = (logrmax - logrmin) / (double)NRADBIN;
  const double invlogrbin = 1.0 / logrbin;


  /** set distribution function */
  /** memory allocation for spherical component(s) */
  initBenchmark_cpu();
  profile **prf, *_prf;
  /** 2 * 2 bins are added in the both edge */
  _prf = (profile  *)malloc(sizeof(profile  ) * kind * NRADBIN);  if( _prf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _prf\n");  }
  prf  = (profile **)malloc(sizeof(profile *) * kind          );  if(  prf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  prf\n");  }
  for(int ii = 0; ii < kind; ii++)
    prf[ii] = _prf + ii * NRADBIN;
#pragma omp parallel for
  for(int ii = 0; ii < kind; ii++)
    for(int jj = 0; jj < NRADBIN; jj++)
      prf[ii][jj].rad = pow(10.0, logrmin + logrbin * (double)jj);
  stopBenchmark_cpu(&execTime.spheAlloc);
  /** memory allocation for disk component(s) */
  initBenchmark_cpu();
  int maxLev;
  disk_data  *disk_info;
  double *disk_hor, *disk_ver, *disk_node_hor, *disk_node_ver, *disk_pot, *disk_dPhidR, *disk_d2PhidR2;
  double *sph_rad, *sph_enc, *sph_rho;
  double *disk_rho, *disk_rhoSum, *disk_rhoTot, *disk_Sigma, *disk_sigmaz, *disk_enc;
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  double *disk_zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  double *spline_xx, *spline_ff, *spline_f2, *spline_bp;
  if( addDisk )
    allocDiskProfile(ndisk, &disk_info, &cfg[skind], &maxLev, prf, skind, logrbin, invlogrbin,
		     &disk_hor, &disk_ver, &disk_node_hor, &disk_node_ver,
		     &disk_pot, &disk_rho, &disk_rhoSum, &disk_rhoTot,
		     &disk_dPhidR, &disk_d2PhidR2,
		     &disk_Sigma, &disk_sigmaz, &disk_enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
		     &disk_zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
		     &sph_rad, &sph_rho, &sph_enc,
		     &spline_xx, &spline_ff, &spline_f2, &spline_bp);
  stopBenchmark_cpu(&execTime.diskAlloc);

  /** set density profile and mass profile for spherical component(s) */
  int nsphere = skind;
  initBenchmark_cpu();
  for(int ii = 0; ii < skind; ii++){
    profile_abel_cfg dummy;    dummy.invRd = 1.0;
    switch( cfg[ii].kind ){
    case   PLUMMER:      setDensityProfilePlummer        (prf[ii],  cfg[ii].rs                                                                       );      break;
    case      KING:      setDensityProfileKing           (prf[ii], &cfg[ii]                                                                          );      break;
    case   BURKERT:      setDensityProfileBurkert        (prf[ii],  cfg[ii].rs                                                                       );      break;
    case HERNQUIST:      setDensityProfileHernquist      (prf[ii],  cfg[ii].rs                                                                       );      break;
    case       NFW:      setDensityProfileNFW            (prf[ii],  cfg[ii].rs                                                                       );      break;
    case     MOORE:      setDensityProfileMoore          (prf[ii],  cfg[ii].rs                                                                       );      break;
    case   EINASTO:      setDensityProfileEinasto        (prf[ii],  cfg[ii].rs, cfg[ii]. einasto_alpha                                               );      break;
    case  APP_KING:      setDensityProfileAppKing        (prf[ii],  cfg[ii].rs,                         cfg[ii].   king_rt                           );      break;
    case TWO_POWER:      setDensityProfileTwoPower       (prf[ii],  cfg[ii].rs, cfg[ii].twopower_alpha, cfg[ii].twopower_beta, cfg[ii].twopower_gamma);      break;
    case TRI_POWER:      setDensityProfileTriPower       (prf[ii],  cfg[ii].rs, cfg[ii].tripower_rout, cfg[ii].twopower_alpha, cfg[ii].twopower_beta, cfg[ii].twopower_gamma, cfg[ii].tripower_delta, cfg[ii].tripower_epsilon);      break;
    case APP_EVANS:      setDensityProfileAppLoweredEvans(prf[ii],  cfg[ii].rs, cfg[ii].alevans_alpha, cfg[ii].alevans_rc, cfg[ii].alevans_beta, cfg[ii].alevans_rt, 1.0 / cfg[ii].alevans_wt);      break;
    case TABLE_RHO:      setDensityProfileTable          (prf[ii],  cfg[ii].rs, cfg[ii].table                                                        );      break;
    case TABLE_SIG:      readColumnDensityProfileTable   (prf[ii],  cfg[ii].rs, cfg[ii].table, cfg[ii]                                               );      break;
    case SPHSERSIC:      execAbelTransform               (prf[ii],  cfg[ii]   , MINRAD, rmax, dummy                                                  );      break;
    case SIGTWOPOW:      execAbelTransform               (prf[ii],  cfg[ii]   , MINRAD, rmax, dummy                                                  );      break;
    case CENTRALBH:      setContributionByCentralBH      (prf[ii],  cfg[ii]                                                                          );      nsphere--;      break;
    default:      __KILL__(stderr, "ERROR: undefined profile ``%d'' specified\n", cfg[ii].kind);      break;
    }/* switch( cfg[ii].kind ){ */

    if( cfg[ii].kind != CENTRALBH )
      integrateDensityProfile(prf[ii], &cfg[ii]
#ifndef ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
			      , logrbin
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
			      );
  }/* for(int ii = 0; ii < skind; ii++){ */

  if( nsphere != skind ){
    for(int kk = nsphere; kk < skind; kk++)
      if( cfg[kk].kind != CENTRALBH ){	__KILL__(stderr, "ERROR: central BH must be last component(s) in spherical components.\n\tModify \"%s/%s\".\n", CFGFOLDER, fcfg);      }
  }/* if( nsphere != skind ){ */

  /** evaluate sum of density, enclosed mass and potential of all spherical component(s) */
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    double rho = 0.0;
    double enc = 0.0;
    double psi = 0.0;
    for(int kk = 0; kk < skind; kk++){
      rho += prf[kk][ii].rho;
      enc += prf[kk][ii].enc;
      psi += prf[kk][ii].psi;
    }/* for(int kk = 0; kk < skind; kk++){ */

    for(int kk = 0; kk < kind; kk++){
      prf[kk][ii].rho_tot = rho;
      prf[kk][ii].enc_tot = enc;
      prf[kk][ii].psi_tot = psi;
    }/* for(int kk = 0; kk < kind; kk++){ */
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */
  stopBenchmark_cpu(&execTime.spheProf);

  /* set density profile for the disk component(s) */
  if( addDisk ){
    /** set disk_radius, disk_height, disk_pot */
    initBenchmark_cpu();
    makeDiskPotentialTable(ndisk, maxLev, disk_info);
    stopBenchmark_cpu(&execTime.diskTbl);

    /** set profile of spherical averaged density, mass and potential */
    initBenchmark_cpu();
    integrateSphericalDensityProfile(ndisk, maxLev, disk_info);

#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      double rho = 0.0;
      double enc = 0.0;
      double psi = 0.0;
      for(int kk = skind; kk < kind; kk++){
	rho += prf[kk][ii].rho;
	enc += prf[kk][ii].enc;
	psi += prf[kk][ii].psi;
      }/* for(int kk = skind; kk < kind; kk++){ */

      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].rho_tot += rho;
	prf[kk][ii].enc_tot += enc;
	prf[kk][ii].psi_tot += psi;
      }/* for(int kk = 0; kk < skind; kk++){ */
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */
    stopBenchmark_cpu(&execTime.diskProf);
  }/* if( addDisk ){ */

  /** integrate Eddington's formula numerically */
  initBenchmark_cpu();
  dist_func **fene, *_fene;
  _fene = (dist_func  *)malloc(sizeof(dist_func  ) * nsphere * NENEBIN);  if( _fene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _fene\n");  }
  fene  = (dist_func **)malloc(sizeof(dist_func *) * nsphere          );  if(  fene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  fene\n");  }
  for(int ii = 0; ii < nsphere; ii++)
    fene[ii] = _fene + ii * NENEBIN;
  integrateEddingtonFormula(nsphere, prf,
#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
			    cfg,
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
			    fene);
  stopBenchmark_cpu(&execTime.eddington);

#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
  initBenchmark_cpu();
  calcVelocityDispersionProfile(nsphere, prf,
#ifdef  ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
				cfg,
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_DF
				fene);
  stopBenchmark_cpu(&execTime.vdisp);
#endif//MAKE_VELOCITY_DISPERSION_PROFILE

#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  initBenchmark_cpu();
  calcColumnDensityProfile(nsphere, prf,
#ifndef ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
			   logrmax,
#endif//ADOPT_DOUBLE_EXPONENTIAL_FORMULA_FOR_PROFILE
			   cfg);
  stopBenchmark_cpu(&execTime.column);
#endif//MAKE_COLUMN_DENSITY_PROFILE

  if( addDisk ){
    initBenchmark_cpu();
    /** differentiate potential along the radial direction on the equatorial plane */
    diffAxisymmetricPotential(maxLev, disk_info[0]);
    /** set velocity dispersion in vertical direction */
    calcVerticalVdisp(ndisk, maxLev, disk_info);

    for(int ii = skind; ii < kind; ii++)
      cfg[ii].vdispz0 = disk_info[ii - skind].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, 0)];

    stopBenchmark_cpu(&execTime.diskVel);
  }/* if( addDisk ){ */


  /** set particle distribution */
  initBenchmark_cpu();
  iparticle body;
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
    ((int)Ntot, &body, &idx, &pos, &acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     &ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
     &vel, &ti
#else///BLOCK_TIME_STEP
     &vx, &vy, &vz
#endif//BLOCK_TIME_STEP
     );
  stopBenchmark_cpu(&execTime.bodyAlloc);

  /** parallel region for OpenMP */
#ifdef  USE_SFMTJUMP
#pragma omp parallel private(Nuse)
#endif//USE_SFMTJUMP
  {
    /** initialize pseudo random number generator */
    rand_state *rand;
    initRandNum(&rand);
#ifdef  USE_SFMTJUMP
    for(int ii = 0; ii < omp_get_thread_num(); ii++)
      SFMT_jump(rand, SFMTJUMP_10_100);
#endif//USE_SFMTJUMP

    /** create spherical particle distribution */
#ifdef  USE_SFMTJUMP
#pragma omp barrier
#pragma omp master
#endif//USE_SFMTJUMP
    initBenchmark_cpu();
    Nuse = 0;
    for(int ii = 0; ii < skind; ii++){
      /** distribute spheroid particles */
      if( cfg[ii].kind != CENTRALBH )
	cfg[ii].Ecut = distributeSpheroidParticles(&Nuse, body, CAST_D2R(cfg[ii].Mtot / (double)cfg[ii].num), cfg[ii], &prf[ii][2], fene[ii], rand);
      else{
#ifdef  USE_SFMTJUMP
#pragma omp single nowait
#endif//USE_SFMTJUMP
	{
	  fprintf(stdout, "#\n#\n# start distributing spherical component particles (%zu bodies: [%zu:%zu])\n", cfg[ii].num, Nuse, Nuse + cfg[ii].num - 1);
	  fflush(stdout);
	  body.pos[Nuse].x = body.pos[Nuse].y = body.pos[Nuse].z = ZERO;	  body.pos[Nuse].m   = CAST_D2R(cfg[ii].Mtot);
	  body.acc[Nuse].x = body.acc[Nuse].y = body.acc[Nuse].z = ZERO;	  body.acc[Nuse].pot = ZERO;
#ifdef  BLOCK_TIME_STEP
	  body.vel[Nuse].x = body.vel[Nuse].y = body.vel[Nuse].z = ZERO;
#else///BLOCK_TIME_STEP
	  body.vx[Nuse] = body.vy[Nuse] = body.vz[Nuse] = ZERO;
#endif//BLOCK_TIME_STEP
	  body.idx[Nuse] = Nuse;
	  fprintf(stdout, "# finish distributing spherical component particles (%zu bodies)\n", cfg[ii].num);
	  fflush(stdout);
	}
	Nuse++;
      }/* else{ */

      /** shift center-of-mass, remove bulk motion */
      shiftCenter(cfg[ii].num, Nuse - cfg[ii].num, body, cfg[ii]);
    }/* for(int ii = 0; ii < skind; ii++){ */
#ifdef  USE_SFMTJUMP
#pragma omp barrier
#pragma omp master
#endif//USE_SFMTJUMP
    stopBenchmark_cpu(&execTime.spheDist);

    /** add disk component if required */
    if( addDisk ){
#ifdef  USE_SFMTJUMP
#pragma omp barrier
#pragma omp master
#endif//USE_SFMTJUMP
      initBenchmark_cpu();

#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
      static double Wsphe;
#pragma omp single
      Wsphe = 0.0;
      double Wsphe_loc = 0.0;
#pragma omp for nowait
      for(ulong ii = 0; ii < Nuse; ii++)
	Wsphe_loc += HALF * body.pos[ii].m * body.acc[ii].pot;
#pragma omp atomic
      Wsphe += Wsphe_loc;
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION

      for(int ii = 0; ii < ndisk; ii++){
	/** distribute disk particles */
	distributeDiskParticles(&Nuse, body, CAST_D2R(disk_info[ii].cfg->Mtot / (double)disk_info[ii].cfg->num), maxLev, disk_info[ii], rand);

	/** flip velocity vector to generate retrograding particles */
	if( disk_info[ii].cfg->retrogradeFrac > 0.0 )
	  retrograder(disk_info[ii].cfg->num, Nuse - disk_info[ii].cfg->num, body, disk_info[ii].cfg->retrogradeFrac);

	/** shift center-of-mass, remove bulk motion */
	shiftCenter(disk_info[ii].cfg->num, Nuse - disk_info[ii].cfg->num, body, *disk_info[ii].cfg);
      }/* for(int ii = 0; ii < ndisk; ii++){ */

#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
#pragma omp single
      {
	double Tdisk = 0.0;

	fprintf(stdout, "# simple check based on Ostriker--Peebles criterion: T / |W| > ~0.14 is unstable to a bar-like mode\n");

	for(int ii = 0; ii < ndisk; ii++){
	  const double tt = disk_info[ii].cfg->Tdisk / (-disk_info[ii].cfg->Wdisk);
	  fprintf(stdout, "#\tT_disk%d / |W| = %e; i.e., %s to a bar-like mode (W = \\int dV \\rho_disk%d \\Phi / 2)\n", ii, tt, (tt > 0.14) ? "unstable" : "  stable", ii);

	  Tdisk += disk_info[ii].cfg->Tdisk;
	  Wsphe += disk_info[ii].cfg->Wdisk;
	}/* for(int ii = 0; ii < ndisk; ii++){ */

#if 0
	extern const double energy2astro;
	fprintf(stdout, "# T = %e, W = %e\n", Tdisk * energy2astro, -Wsphe * energy2astro);
#endif

	Tdisk /= (-Wsphe);
	fprintf(stdout, "#\tT / |W| = %e; i.e., %s to a bar-like mode (W = \\int dV \\rho \\Phi / 2)\n", Tdisk, (Tdisk > 0.14) ? "unstable" : "  stable");
      }
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION

#ifdef  USE_SFMTJUMP
#pragma omp barrier
#pragma omp master
#endif//USE_SFMTJUMP
      stopBenchmark_cpu(&execTime.diskDist);
    }/* if( addDisk ){ */

    /** release pseudo random number generator */
    freeRandNum(rand);
  }

  /** generate tables for fixed potential field */
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  initBenchmark_cpu();
  real *rad_pot;
  pot2 *Phi_pot;
  potential_field *pot_tbl, pot_tbl_sphe, pot_tbl_disk;
  allocPotentialField(&rad_pot, &Phi_pot, &pot_tbl, N_EXT_POT_SPHE, kind, &pot_tbl_sphe, skind, &pot_tbl_disk);

  real *RR_diskpot, *zz_diskpot, *Phi_diskpot;
  disk_potential diskpot;
  if( addDisk )
    allocDiskPotential(&RR_diskpot, &zz_diskpot, &Phi_diskpot, maxLev, NDISKBIN_HOR, NDISKBIN_VER, &diskpot);

  genExtPotTbl1D(kind, prf, pot_tbl);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  genSuperposedPotFld1D(kind, skind, prf, &pot_tbl_sphe, &pot_tbl_disk);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  superposePotFld1D(kind, skind, pot_tbl, &pot_tbl_sphe, &pot_tbl_disk);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  if( addDisk )
    extractDiskPotential(maxLev, disk_info[0], pot_tbl_disk, &diskpot);

  stopBenchmark_cpu(&execTime.external);
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  if( skind == 0 )
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      prf[0][ii].enc_tot  = prf[ 0][ii].enc;
      for(int jj = 1; jj < ndisk; jj++)
      prf[0][ii].enc_tot += prf[jj][ii].enc;
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */


  /** evaluate observable quantities */
  initBenchmark_cpu();
  evaluateObservables(kind, skind, cfg, prf);
  if( addDisk )
    getEffectiveRadius(ndisk, maxLev, disk_info);
  stopBenchmark_cpu(&execTime.observe);


  /** output fundamental quantities of the disk component */
  if( addDisk ){
    initBenchmark_cpu();
    writeDiskData(file, ndisk, maxLev, disk_info);
    stopBenchmark_cpu(&execTime.diskInfo);
  }/* if( addDisk ){ */

  /** write fundamental information */
  initBenchmark_cpu();
  outputFundamentalInformation(unit, kind, skind, cfg, prf, Ntot, eps, snapshotInterval, ft, file);
  outputRadialProfiles(kind, prf,
#ifdef  USE_HDF5_FORMAT
		       cfg, eps,
#endif//USE_HDF5_FORMAT
		       file);
  outputDistributionFunction(nsphere, fene, file);
#ifdef  USE_HDF5_FORMAT
  outputRepresentativeQuantities(kind, cfg, skind, prf, ndisk, maxLev, disk_info, file);
#endif//USE_HDF5_FORMAT
  stopBenchmark_cpu(&execTime.spheInfo);

  /** write particle data */
  double time = 0.0;
  initBenchmark_cpu();

#   if  !defined(WRITE_IN_TIPSY_FORMAT) && !defined(WRITE_IN_GALACTICS_FORMAT)
  double  dt  = 0.0;
  int   last  = 1;
  ulong steps = 0;

  writeSettings(unit, Ntot, eps, eta, ft, snapshotInterval, saveInterval, file);
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
#ifndef RUN_WITHOUT_GOTHIC
#ifdef  MONITOR_ENERGY_ERROR
  static energyError relEneErr = {1.0, DBL_MIN};
#endif//MONITOR_ENERGY_ERROR
  static rebuildTree rebuild;
  static measuredTime measured;
  static autoTuningParam rebuildParam;
  static brentStatus status;
  static brentMemory memory;
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT

  writeTentativeData(time, dt, steps, Ntot, body, file, &last
#ifdef  USE_HDF5_FORMAT
		     , hdf5type
#ifndef RUN_WITHOUT_GOTHIC
		     , rebuild, measured, rebuildParam, status, memory
#ifdef  MONITOR_ENERGY_ERROR
		     , relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
		     );
  updateConfigFile(last, file);

#else///!defined(WRITE_IN_TIPSY_FORMAT) && !defined(WRITE_IN_GALACTICS_FORMAT)

#ifdef  WRITE_IN_TIPSY_FORMAT
  writeTipsyFile(time, CAST_R2F(eps), (int)Ntot, body, file);
#endif//WRITE_IN_TIPSY_FORMAT

#ifdef  WRITE_IN_GALACTICS_FORMAT
  int hidx = 0;
  for(int kk = 0; kk < kind; kk++){
    writeGalactICSFile(time, hidx, cfg[kk].num, body, file, kk);
    hidx += cfg[kk].num;
  }/* for(int kk = 0; kk < kind; kk++){ */
#endif//WRITE_IN_GALACTICS_FORMAT

#endif//!defined(WRITE_IN_TIPSY_FORMAT) && !defined(WRITE_IN_GALACTICS_FORMAT)

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  writeFixedPotentialTable
    (unit, pot_tbl_sphe, skind, pot_tbl
#ifdef  USE_HDF5_FORMAT
     , hdf5type
#else///USE_HDF5_FORMAT
#ifdef  WRITE_IN_GALACTICS_FORMAT
     , false
#else///WRITE_IN_GALACTICS_FORMAT
     , true
#endif//WRITE_IN_GALACTICS_FORMAT
#endif//USE_HDF5_FORMAT
     , file);

  if( addDisk )
    writeFixedDiskPotential
      (unit, diskpot
#ifdef  USE_HDF5_FORMAT
       , hdf5type
#else///USE_HDF5_FORMAT
#ifdef  WRITE_IN_GALACTICS_FORMAT
       , false
#else///WRITE_IN_GALACTICS_FORMAT
       , true
#endif//WRITE_IN_GALACTICS_FORMAT
#endif//USE_HDF5_FORMAT
       , file);

#endif//SET_EXTERNAL_POTENTIAL_FIELD

#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
#endif//USE_HDF5_FORMAT

  stopBenchmark_cpu(&execTime.file);


  /** write result of measured breakdown */
  double ttot = 0.0;
  ttot += execTime.bodyAlloc + execTime.file;
  ttot += execTime.spheAlloc + execTime.spheProf + execTime.spheDist + execTime.spheInfo + execTime.eddington;
  ttot += execTime.diskAlloc + execTime.diskProf + execTime.diskDist + execTime.diskInfo + execTime.diskTbl + execTime.diskVel;
  ttot += execTime.observe;
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
  ttot += execTime.vdisp;
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  ttot += execTime.column;
#endif//MAKE_COLUMN_DENSITY_PROFILE
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  ttot += execTime.external;
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  FILE *fp;
  char filename[128];
  sprintf(filename, "%s/%s.N%zu.%s.%s.txt", LOGFOLDER, file, Ntot, "magi", "breakdown");
  if( 0 != access(filename, F_OK) ){
    fp = fopen(filename, "w");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	    "ttot", "bodyAlloc", "file",
	    "spheAlloc", "spheProf", "spheDist", "spheInfo", "eddington",
	    "diskAlloc", "diskProf", "diskDist", "diskInfo", "diskTbl", "diskVel",
	    "observe",
    "vdisp", "column",
    "external");
    fclose(fp);
  }/* if( 0 != access(filename, F_OK) ){ */
  fp = fopen(filename, "a");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e",
	  ttot, execTime.bodyAlloc, execTime.file,
	  execTime.spheAlloc, execTime.spheProf, execTime.spheDist, execTime.spheInfo, execTime.eddington,
	  execTime.diskAlloc, execTime.diskProf, execTime.diskDist, execTime.diskInfo, execTime.diskTbl, execTime.diskVel,
	  execTime.observe);
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
  fprintf(fp, "\t%e", execTime.vdisp);
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  fprintf(fp, "\t%e", execTime.column);
#endif//MAKE_COLUMN_DENSITY_PROFILE
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  fprintf(fp, "\t%e", execTime.external);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  fprintf(fp, "\n");
  fclose(fp);

  fprintf(stdout, "#\n#\n# benchmark result:\n");
  fprintf(stdout, "# total elapsed time is %e s\n", ttot);
  ttot = 100.0 / ttot;
  if( addDisk )    fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.diskTbl  , execTime.diskTbl   * ttot, "generating potential--density pair of disk component(s)");
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.spheDist , execTime.spheDist  * ttot, "distributing N-body particles of spherical component(s)");
  if( addDisk )    fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.diskDist , execTime.diskDist  * ttot, "distributing N-body particles of disk component(s)");
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.eddington, execTime.eddington * ttot, "calculating DF(s) of spherical component(s) using Eddington formula");
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.spheProf , execTime.spheProf  * ttot, "generating radial profile of spherical component(s)");
  if( addDisk )    fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.diskProf , execTime.diskProf  * ttot, "calculating spherical averaged profile of disk component(s)");
  if( addDisk )    fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.diskVel  , execTime.diskVel   * ttot, "calculating vertical velocity dispersion of disk component(s)");
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.file    , execTime.file     * ttot, "writing particle data");
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.spheInfo, execTime.spheInfo * ttot, "writing fundamental data of spherical component(s)");
  if( addDisk )    fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.diskInfo, execTime.diskInfo * ttot, "writing fundamental data of disk component(s)");
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.vdisp , execTime.vdisp  * ttot, "calculating velocity dispersion profile of spherical component(s)");
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.column, execTime.column * ttot, "calculating column density profile of spherical component(s)");
#endif//MAKE_COLUMN_DENSITY_PROFILE
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.observe, execTime.observe * ttot, "evaluating observables of spherical component(s)");
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.spheAlloc, execTime.spheAlloc * ttot, "memory allocation for spherical component(s)");
  if( addDisk )    fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.diskAlloc, execTime.diskAlloc * ttot, "memory allocation for disk component(s)");
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.bodyAlloc, execTime.bodyAlloc * ttot, "memory allocation for N-body particles");
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  fprintf(stdout, "# %e s (%5.2f%%) for %s\n", execTime.external, execTime.external * ttot, "execute cubic spline interpolation to generate external potential field");
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#ifdef  PRINT_OUT_ASCII_DATA_FOR_QUICK_CHECK
  /** print out particle distribution in ASCII format for quick check */
  sprintf(filename, "%s/%s.%s.txt", DATAFOLDER, file, "ascii");
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fprintf(fp, "#x\ty\tz\tvx\tvy\tvz\n");
  ulong head = 0;
  for(int kk = 0; kk < kind; kk++){
    for(ulong ii = head; ii < head + IMIN(cfg[kk].num, 1024); ii++)
      fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\n",
	      body.pos[ii].x, body.pos[ii].y, body.pos[ii].z,
#ifdef  BLOCK_TIME_STEP
	      body.vel[ii].x, body.vel[ii].y, body.vel[ii].z
#else///BLOCK_TIME_STEP
	      body.vx[ii], body.vy[ii], body.vz[ii]
#endif//BLOCK_TIME_STEP
	      );
    fprintf(fp, "\n");
    head += cfg[kk].num;
  }/* for(int kk = 0; kk < kind; kk++){ */
  fclose(fp);

  sprintf(filename, "%s/ascii.gp", DATAFOLDER);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fprintf(fp, "se si sq\n");
  fprintf(fp, "se st da d\n");
  fprintf(fp, "se xl 'x'\n");
  fprintf(fp, "se yl 'z'\n");
  fprintf(fp, "p ");
  for(int ii = 0; ii < kind; ii++){
    fprintf(fp, "'%s.ascii.txt' u 1:3 ev :::%d::%d", file, ii, ii);
    if( ii != (kind - 1) )
      fprintf(fp, ",\\\n  ");
  }/* for(int ii = 0; ii < kind; ii++){ */
  fprintf(fp, "\n");
  fclose(fp);
#endif//PRINT_OUT_ASCII_DATA_FOR_QUICK_CHECK


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
  free(cfg);
  free(prf);
  free(_prf);
  free(fene);
  free(_fene);
  if( addDisk )
    freeDiskProfile
      (ndisk, disk_info,
       disk_hor, disk_ver, disk_node_hor, disk_node_ver,
       disk_pot, disk_rho, disk_rhoSum, disk_rhoTot, disk_dPhidR, disk_d2PhidR2,
       disk_Sigma, disk_sigmaz, disk_enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
       disk_zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
       sph_rad, sph_rho, sph_enc,
       spline_xx, spline_ff, spline_f2, spline_bp);

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  freePotentialField(rad_pot, Phi_pot, pot_tbl);
  if( addDisk )
    freeDiskPotential(RR_diskpot, zz_diskpot, Phi_diskpot);
#endif//SET_EXTERNAL_POTENTIAL_FIELDp

  return (0);
}


/* for evaluating relaxation time of 2-body relaxation */
/**
 * @fn findIdx_rad
 *
 * @brief Find a data element in the given array corresponding to the given radius.
 *
 * @param (rad) rad
 * @param (prf) radial profile of the component
 * @return (ll) the corresponding lower index
 * @return (rr) the corresponding upper index
 */
static inline void findIdx_rad(const double rad, profile * restrict prf, int * restrict ll, int * restrict rr)
{
  bool bisection = true;
  *ll =           0;
  *rr = NRADBIN - 1;

  if( bisection == true )    if( fabs(prf[*ll].rad - rad) / rad < DBL_EPSILON ){      bisection = false;      *rr = (*ll) + 1;    }
  if( bisection == true )    if( fabs(prf[*rr].rad - rad) / rad < DBL_EPSILON ){      bisection = false;      *ll = (*rr) - 1;    }

  while( bisection ){
    const uint cc = ((uint)(*ll) + (uint)(*rr)) >> 1;

    if( (prf[cc].rad - rad) * (prf[*ll].rad - rad) <= 0.0 )      *rr = (int)cc;
    else                                                         *ll = (int)cc;

    if( (1 + (*ll)) == (*rr) )
      break;
  }/* while( bisection ){ */
}


/**
 * @fn outputFundamentalInformation
 *
 * @brief Print out fundamental information on the generated system.
 *
 * @param (unit) unit system
 * @param (kind) number of components
 * @param (skind) number of spherical symmetric components
 * @param (cfg) physical properties of the components
 * @param (prf) radial profile of the components
 * @param (Ntot) total number of N-body particles
 * @return (eps) value of Plummer softening length
 * @return (SnapshotInterval) time interval to write snapshot files
 * @return (ft) finish time of the simulation
 * @param (file) name of the simulation
 */
void outputFundamentalInformation
(const int unit, const int kind, const int skind, profile_cfg *cfg, profile **prf, const ulong Ntot, const real eps, const double snapshotInterval, const double ft, char file[])
{
  __NOTE__("%s\n", "start");

  FILE *fp;
  char filename[256], date[64];

  /** output useful information for multi-component analysis */
  sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

  fprintf(fp, "%d\n", unit);
  fprintf(fp, "%d\t%d\n", kind, skind);
  for(int ii = 0; ii < kind; ii++)
    fprintf(fp, "%zu\n", cfg[ii].num);

  fclose(fp);


  /** output fundamental information of the particle distribution */
  sprintf(filename, "%s/%s.info.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

  /** output global settings */
  getPresentDateInStrings(date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Fundamental information about the particle distribution\n");
  fprintf(fp, "\tgenerated on %s", date);
  fprintf(fp, "Physical quantities in Computational and Astrophysical units is listed\n");
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "#############################################################################\n");
  double Mtot = 0.0;
  for(int ii = 0; ii < kind; ii++)
    Mtot += cfg[ii].Mtot;
  fprintf(fp, "Total mass of the system  Mtot is %e (= %e %s)\n" , Mtot, Mtot * mass2astro, mass_astro_unit_name);
  fprintf(fp, "Total number of particles Ntot is %zu (= 2^%u)\n", Ntot, ilog2((uint)Ntot));
  fprintf(fp, "Number of components      kind is %d\n" , kind);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Length of Plummer softening  is %e (= %e %s)\n", eps, eps * length2astro, length_astro_unit_name);
  fprintf(fp, "Snapshot interval            is %e (= %e %s)\n", snapshotInterval, snapshotInterval * time2astro, time_astro_unit_name);
  fprintf(fp, "Final time of the simulation is %e (= %e %s)\n",               ft,               ft * time2astro, time_astro_unit_name);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "#############################################################################\n\n");

  ulong ihead = 0;
  for(int ii = 0; ii < kind; ii++){
    /** output settings for individual component */
    fprintf(fp, "#############################################################################\n");
    fprintf(fp, "#############################################################################\n");
    fprintf(fp, "%d-th component: ", ii);
    switch( cfg[ii].kind ){
    case   PLUMMER:      fprintf(fp,                        "Plummer model\n");      break;
    case      KING:      fprintf(fp,                           "King model\n");      break;
    case   BURKERT:      fprintf(fp,                        "Burkert model\n");      break;
    case HERNQUIST:      fprintf(fp,                      "Hernquist model\n");      break;
    case       NFW:      fprintf(fp,                            "NFW model\n");      break;
    case     MOORE:      fprintf(fp,                          "Moore model\n");      break;
    case   EINASTO:      fprintf(fp,                        "Einasto model\n");      break;
    case  APP_KING:      fprintf(fp,         "King model in empirical form\n");      break;
    case TWO_POWER:      fprintf(fp,   "Two-power (alpha-beta-gamma) model\n");      break;
    case TRI_POWER:      fprintf(fp,                    "Three-power model\n");      break;
    case APP_EVANS:      fprintf(fp,     "Approximated lowered Evans model\n");      break;
    case TABLE_RHO:      fprintf(fp,        "Density profile in table form\n");      break;
    case TABLE_SIG:      fprintf(fp, "Column density profile in table form\n");      break;
    case SPHSERSIC:      fprintf(fp,           "Sersic profile (spherical)\n");      break;
    case SIGTWOPOW:      fprintf(fp,        "Two-power model in projection\n");      break;
    case  EXP_DISK:      fprintf(fp,                     "Exponential disk\n");      break;
    case    SERSIC:      fprintf(fp,                "Sersic profile (disk)\n");      break;
    case  TBL_DISK:      fprintf(fp,         "Disk component in table form\n");      break;
    case CENTRALBH:      fprintf(fp,           "Central massive black hole\n");      break;
    default:
      __KILL__(stderr, "ERROR: undefined profile ``%d'' specified\n", cfg[ii].kind);
      break;
    }/* switch( cfg[ii].kind ){ */
    fprintf(fp, "#############################################################################\n");
    fprintf(fp, "Total number of particles Ntot is %zu (about 2^%u)\n", cfg[ii].num, ilog2((int)cfg[ii].num));
    fprintf(fp, "Range of idx for the component is [%zu, %zu]\n", ihead, ihead + cfg[ii].num - 1);
    if( cfg[ii].kind < 0 )
      if( cfg[ii].retrogradeFrac > 0.0 ){
	ulong retroNum, retroHead;
	pickRetrogradingParticles(cfg[ii].num, ihead, cfg[ii].retrogradeFrac, &retroNum, &retroHead);
	fprintf(fp, "\tRange of idx for   prograding particles is [%zu, %zu]\n", ihead, ihead + cfg[ii].num - retroNum - 1);
	fprintf(fp, "\tRange of idx for retrograding particles is [%zu, %zu]\n", retroHead, retroHead + retroNum - 1);
      }/* if( cfg[ii].retrogradeFrac > 0.0 ){ */
    ihead += cfg[ii].num;
    fprintf(fp, "Mass of each N-body particle m is %e (= %e %s)\n", cfg[ii].Mtot / (double)cfg[ii].num, cfg[ii].Mtot / (double)cfg[ii].num * mass2astro, mass_astro_unit_name);
    fprintf(fp, "#############################################################################\n");
    if( cfg[ii].kind == TABLE_RHO )
      fprintf(fp, "Given density profile is written in %s\n", cfg[ii].table);
    if( cfg[ii].kind == TABLE_SIG )
      fprintf(fp, "Given column density profile is written in %s\n", cfg[ii].table);
    fprintf(fp, "Total mass of the component Mtot is %e (= %e %s)\n", cfg[ii].Mtot, cfg[ii].Mtot * mass2astro, mass_astro_unit_name);
    if( cfg[ii].kind != CENTRALBH )
      fprintf(fp, "Scale radius of the component rs is %e (= %e %s)\n", cfg[ii].rs, cfg[ii].rs * length2astro, length_astro_unit_name);
    if( cfg[ii].kind == KING ){
      fprintf(fp, "Dimensionless King parameter  W0 is %e\n", cfg[ii].king_W0);
#ifdef  KING_CENTRAL_CUSP
      fprintf(fp, "Dimensionless slope      dW/dx_0 is %e\n", cfg[ii].king_dWdx_0);
#endif//KING_CENTRAL_CUSP
      fprintf(fp, "Tidal radius of the component rt is %e (= %e %s)\n", cfg[ii].king_rt, cfg[ii].king_rt * length2astro, length_astro_unit_name);
      fprintf(fp, "Concentration paramter         c is %e\n", cfg[ii].king_c);
    }/* if( cfg[ii].kind == KING ){ */
    if( cfg[ii].kind == APP_KING ){
      fprintf(fp, "Tidal radius of the component rt is %e (= %e %s)\n", cfg[ii].king_rt, cfg[ii].king_rt * length2astro, length_astro_unit_name);
      fprintf(fp, "Concentration paramter         c is %e\n", log10(cfg[ii].king_rt / cfg[ii].rs));
    }/* if( cfg[ii].kind == APP_KING ){ */
    if( cfg[ii].kind == EINASTO )
      fprintf(fp, "Shape parameter            alpha is %e\n", cfg[ii].einasto_alpha);
    if( cfg[ii].kind == TWO_POWER ){
      fprintf(fp, "   Inner power-law index   alpha is %e\n", cfg[ii].twopower_alpha);
      fprintf(fp, "Internal power-law index    beta is %e\n", cfg[ii].twopower_beta);
      fprintf(fp, "   Outer power-law index   gamma is %e\n", cfg[ii].twopower_gamma);
    }/* if( cfg[ii].kind == TWO_POWER ){ */
    if( cfg[ii].kind == TRI_POWER ){
      fprintf(fp, "Outer transition radius         rout is %e (= %e %s)\n", cfg[ii].tripower_rout, cfg[ii].tripower_rout * length2astro, length_astro_unit_name);
      fprintf(fp, "   Innermost power-law index   alpha is %e\n", cfg[ii].twopower_alpha);
      fprintf(fp, "Transitional power-law index    beta is %e\n", cfg[ii].twopower_beta);
      fprintf(fp, "Intermediate power-law index   gamma is %e\n", cfg[ii].twopower_gamma);
      fprintf(fp, "Transitional power-law index   delta is %e\n", cfg[ii].tripower_delta);
      fprintf(fp, "   Outermost power-law index epsilon is %e\n", cfg[ii].tripower_epsilon);
    }/* if( cfg[ii].kind == TRI_POWER ){ */
    if( cfg[ii].kind == APP_EVANS ){
      fprintf(fp, "Second scale radius           rc is %e (= %e %s)\n", cfg[ii].alevans_rc, cfg[ii].alevans_rc * length2astro, length_astro_unit_name);
      fprintf(fp, "Inner power-law index      alpha is %e\n", cfg[ii].alevans_alpha);
      fprintf(fp, "Outer power-law index       beta is %e\n", cfg[ii].alevans_beta);
      fprintf(fp, "Exponential cutoff radius     rt is %e (= %e %s)\n", cfg[ii].alevans_rt, cfg[ii].alevans_rt * length2astro, length_astro_unit_name);
      fprintf(fp, "Exponential cutoff width      wt is %e (= %e %s)\n", cfg[ii].alevans_wt, cfg[ii].alevans_wt * length2astro, length_astro_unit_name);
      cfg[ii].rs = cfg[ii].alevans_rc;
      fprintf(fp, "Re-defined scale radius       rs is %e (= %e %s)\n", cfg[ii].rs, cfg[ii].rs * length2astro, length_astro_unit_name);
    }/* if( cfg[ii].kind == APP_EVANS ){ */
    if( cfg[ii].kind == SPHSERSIC ){
      fprintf(fp, "Sersic index                   n is %e\n", cfg[ii].n_sersic);
      fprintf(fp, "Dimensionless scale factor     b is %e\n", cfg[ii].b_sersic);
    }/* if( cfg[ii].kind == SPHSERSIC ){ */
    if( cfg[ii].kind == SIGTWOPOW ){
      fprintf(fp, "   Inner power-law index   alpha is %e\n", cfg[ii].twopower_alpha);
      fprintf(fp, "Internal power-law index    beta is %e\n", cfg[ii].twopower_beta);
      fprintf(fp, "   Outer power-law index   gamma is %e\n", cfg[ii].twopower_gamma);
    }/* if( cfg[ii].kind == SIGTWOPOW ){ */
    if( (ii < skind) && (cfg[ii].kind != CENTRALBH) ){
      fprintf(fp, "Scale density of the component rho0 is %e (= %e %s)\n", cfg[ii].rho0, cfg[ii].rho0 * density2astro, density_astro_unit_name);
    }/* if( (ii < skind) && (cfg[ii].kind != CENTRALBH) ){ */
    fprintf(fp, "#############################################################################\n");

    if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) || (cfg[ii].kind == TBL_DISK) ){
      if( cfg[ii].kind == SERSIC ){
	fprintf(fp, "Sersic index                   n is %e\n", cfg[ii].n_sersic);
	fprintf(fp, "Dimensionless scale factor     b is %e\n", cfg[ii].b_sersic);
      }/* if( cfg[ii].kind == SERSIC ){ */
      fprintf(fp, "Scale height of the component zd is %e (= %e %s)\n", cfg[ii].zd, cfg[ii].zd * length2astro, length_astro_unit_name);
      fprintf(fp, "``ENABLE_VARIABLE_SCALE_HEIGHT'' is %s.\n",
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
	      "on"
#else///ENABLE_VARIABLE_SCALE_HEIGHT
	      "off"
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
	      );
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      fprintf(fp, "Dimming height of the component is %f times zd\n", DISK_DIMMING_HEIGHT);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      fprintf(fp, "Retrograding fraction            is %e\n", cfg[ii].retrogradeFrac);
      fprintf(fp, "Central surface density   Sigma0 is %e (= %e %s)\n", cfg[ii].Sigma0, cfg[ii].Sigma0 * col_density2astro, col_density_astro_unit_name);
      fprintf(fp, "Circular speed at scale radius   is %e (= %e %s)\n", cfg[ii].vcirc_Rd   , cfg[ii].vcirc_Rd    * velocity2astro, velocity_astro_unit_name);
      fprintf(fp, "Maximum circular speed           is %e (= %e %s)\n", cfg[ii].vcirc_max  , cfg[ii].vcirc_max   * velocity2astro, velocity_astro_unit_name);
      fprintf(fp, "Circular speed is maximized      at %e (= %e %s)\n", cfg[ii].vcirc_max_R, cfg[ii].vcirc_max_R *   length2astro,   length_astro_unit_name);
      fprintf(fp, "``ENFORCE_EPICYCLIC_APPROXIMATION'' is %s.\n",
#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
	      "on"
#else///ENFORCE_EPICYCLIC_APPROXIMATION
	      "off"
#endif//ENFORCE_EPICYCLIC_APPROXIMATION
	      );
#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
      fprintf(fp, "Horizontal velocity dispersion   is %e of circular velocity or vertical velocity dispersion (maximum is used)\n", cfg[ii].vdisp_frac);
#else///ENFORCE_EPICYCLIC_APPROXIMATION
      fprintf(fp, "Horizontal velocity dispersion   is %e (= %e %s)\n", cfg[ii].vdispR0  , cfg[ii].vdispR0   * velocity2astro, velocity_astro_unit_name);
#endif//ENFORCE_EPICYCLIC_APPROXIMATION
      fprintf(fp, "Vertical   velocity dispersion   is %e (= %e %s)\n", cfg[ii].vdispz0  , cfg[ii].vdispz0   * velocity2astro, velocity_astro_unit_name);
      fprintf(fp, "Toomre's Q-value at scale radius is %e\n", cfg[ii].toomre);
#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
      fprintf(fp, "Minimum of Toomre's Q-value      is %e at R = %e (= %e %s), when excluding central region (estimation from velocity dispersion profile)\n", cfg[ii].Qmin2, cfg[ii].qminR2, cfg[ii].qminR2 * length2astro, length_astro_unit_name);
#endif//ENFORCE_EPICYCLIC_APPROXIMATION
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      fprintf(fp, "Minimum of Toomre's Q-value      is %e at R = %e (= %e %s), when excluding central region (estimation from scale height profile)\n", cfg[ii].Qmin1, cfg[ii].qminR1, cfg[ii].qminR1 * length2astro, length_astro_unit_name);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      fprintf(fp, "Minimum of Toomre's Q-value      is %e at R = %e (= %e %s)\n", cfg[ii].Qmin0, cfg[ii].qminR0, cfg[ii].qminR0 * length2astro, length_astro_unit_name);
    }/* if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) || (cfg[ii].kind == TBL_DISK) ){ */

    if( cfg[ii].kind != CENTRALBH ){
      fprintf(fp, "#############################################################################\n");
      fprintf(fp, "Representative scales:\n");
      fprintf(fp, "Half-mass radius         r_{1/2} is %e (= %e %s)\n", cfg[ii].rhalf, cfg[ii].rhalf * length2astro, length_astro_unit_name);
      fprintf(fp, "Effective radius         R_{eff} is %e (= %e %s)\n", cfg[ii].Reff, cfg[ii].Reff * length2astro, length_astro_unit_name);

      if( ii < skind ){/* for spherical component(s) */
	fprintf(fp, "#########\n");
	fprintf(fp, "Representative quantities at the center:\n");
	fprintf(fp, "Volume density               rho is %e (= %e %s)\n", prf[ii][0].rho, prf[ii][0].rho * density2astro, density_astro_unit_name);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
	fprintf(fp, "Column density             Sigma is %e (= %e %s)\n", prf[ii][0].Sigma, prf[ii][0].Sigma * col_density2astro, col_density_astro_unit_name);
#endif//MAKE_COLUMN_DENSITY_PROFILE
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
	fprintf(fp, "Velocity dispersion      sigma_r is %e (= %e %s)\n", prf[ii][0].sigr, prf[ii][0].sigr * velocity2astro, velocity_astro_unit_name);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
	fprintf(fp, "LoS velocity dispersion  sig_los is %e (= %e %s)\n", prf[ii][0].slos, prf[ii][0].slos * velocity2astro, velocity_astro_unit_name);
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE

	int ll, rr;
	double alpha, val;

	/* physical quantities @ r = rs */
	findIdx_rad(cfg[ii].rs, prf[ii], &ll, &rr);
	alpha = (cfg[ii].rs - prf[ii][ll].rad) / (prf[ii][rr].rad - prf[ii][ll].rad);

	fprintf(fp, "#########\n");
	fprintf(fp, "Representative quantities at the scale radius:\n");

	val = (1.0 - alpha) * prf[ii][ll].rho + alpha * prf[ii][rr].rho;
	fprintf(fp, "Volume density               rho is %e (= %e %s)\n", val, val * density2astro, density_astro_unit_name);

#ifdef  MAKE_COLUMN_DENSITY_PROFILE
	val = (1.0 - alpha) * prf[ii][ll].Sigma + alpha * prf[ii][rr].Sigma;
	fprintf(fp, "Column density             Sigma is %e (= %e %s)\n", val, val * col_density2astro, col_density_astro_unit_name);
#endif//MAKE_COLUMN_DENSITY_PROFILE

	val = (1.0 - alpha) * prf[ii][ll].enc + alpha * prf[ii][rr].enc;
	fprintf(fp, "Enclosed mass of this component  is %e (= %e %s)\n", val, val * mass2astro, mass_astro_unit_name);

	val = (1.0 - alpha) * prf[ii][ll].enc_tot + alpha * prf[ii][rr].enc_tot;
	fprintf(fp, "Enclosed mass of  all components is %e (= %e %s)\n", val, val * mass2astro, mass_astro_unit_name);

#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
	val = (1.0 - alpha) * prf[ii][ll].sigr + alpha * prf[ii][rr].sigr;
	fprintf(fp, "Velocity dispersion      sigma_r is %e (= %e %s)\n", val, val * velocity2astro, velocity_astro_unit_name);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
	val = (1.0 - alpha) * prf[ii][ll].slos + alpha * prf[ii][rr].slos;
	fprintf(fp, "LoS velocity dispersion  sig_los is %e (= %e %s)\n", val, val * velocity2astro, velocity_astro_unit_name);
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE

      /* physical quantities @ r = r_{1/2} */
	findIdx_rad(cfg[ii].rhalf, prf[ii], &ll, &rr);
	alpha = (cfg[ii].rhalf - prf[ii][ll].rad) / (prf[ii][rr].rad - prf[ii][ll].rad);

	fprintf(fp, "#########\n");
	fprintf(fp, "Representative quantities at the half-mass radius:\n");

	val = (1.0 - alpha) * prf[ii][ll].rho + alpha * prf[ii][rr].rho;
	fprintf(fp, "Volume density               rho is %e (= %e %s)\n", val, val * density2astro, density_astro_unit_name);

#ifdef  MAKE_COLUMN_DENSITY_PROFILE
	val = (1.0 - alpha) * prf[ii][ll].Sigma + alpha * prf[ii][rr].Sigma;
	fprintf(fp, "Column density             Sigma is %e (= %e %s)\n", val, val * col_density2astro, col_density_astro_unit_name);
#endif//MAKE_COLUMN_DENSITY_PROFILE

	val = (1.0 - alpha) * prf[ii][ll].enc + alpha * prf[ii][rr].enc;
	fprintf(fp, "Enclosed mass of this component  is %e (= %e %s)\n", val, val * mass2astro, mass_astro_unit_name);

	val = (1.0 - alpha) * prf[ii][ll].enc_tot + alpha * prf[ii][rr].enc_tot;
	fprintf(fp, "Enclosed mass of  all components is %e (= %e %s)\n", val, val * mass2astro, mass_astro_unit_name);

#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
	val = (1.0 - alpha) * prf[ii][ll].sigr + alpha * prf[ii][rr].sigr;
	fprintf(fp, "Velocity dispersion      sigma_r is %e (= %e %s)\n", val, val * velocity2astro, velocity_astro_unit_name);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
	val = (1.0 - alpha) * prf[ii][ll].slos + alpha * prf[ii][rr].slos;
	fprintf(fp, "LoS velocity dispersion  sig_los is %e (= %e %s)\n", val, val * velocity2astro, velocity_astro_unit_name);
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE

      /* physical quantities @ r = R_{eff} */
	findIdx_rad(cfg[ii].Reff, prf[ii], &ll, &rr);
	alpha = (cfg[ii].Reff - prf[ii][ll].rad) / (prf[ii][rr].rad - prf[ii][ll].rad);

	fprintf(fp, "#########\n");
	fprintf(fp, "Representative quantities at the effective radius:\n");

	val = (1.0 - alpha) * prf[ii][ll].rho + alpha * prf[ii][rr].rho;
	fprintf(fp, "Volume density               rho is %e (= %e %s)\n", val, val * density2astro, density_astro_unit_name);

#ifdef  MAKE_COLUMN_DENSITY_PROFILE
	val = (1.0 - alpha) * prf[ii][ll].Sigma + alpha * prf[ii][rr].Sigma;
	fprintf(fp, "Column density             Sigma is %e (= %e %s)\n", val, val * col_density2astro, col_density_astro_unit_name);
#endif//MAKE_COLUMN_DENSITY_PROFILE

	val = (1.0 - alpha) * prf[ii][ll].enc + alpha * prf[ii][rr].enc;
	fprintf(fp, "Enclosed mass of this component  is %e (= %e %s)\n", val, val * mass2astro, mass_astro_unit_name);

	val = (1.0 - alpha) * prf[ii][ll].enc_tot + alpha * prf[ii][rr].enc_tot;
	fprintf(fp, "Enclosed mass of  all components is %e (= %e %s)\n", val, val * mass2astro, mass_astro_unit_name);

#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
	val = (1.0 - alpha) * prf[ii][ll].sigr + alpha * prf[ii][rr].sigr;
	fprintf(fp, "Velocity dispersion      sigma_r is %e (= %e %s)\n", val, val * velocity2astro, velocity_astro_unit_name);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
	val = (1.0 - alpha) * prf[ii][ll].slos + alpha * prf[ii][rr].slos;
	fprintf(fp, "LoS velocity dispersion  sig_los is %e (= %e %s)\n", val, val * velocity2astro, velocity_astro_unit_name);
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
      }/* if( ii < skind ){ */

      fprintf(fp, "#############################################################################\n");
      if( cfg[ii].cutoff ){
	fprintf(fp, "Cutoff radius of the component   is %e (= %e %s)\n", cfg[ii].rc      , cfg[ii].rc       * length2astro, length_astro_unit_name);
	fprintf(fp, "Cutoff  width of the component   is %e (= %e %s)\n", cfg[ii].rc_width, cfg[ii].rc_width * length2astro, length_astro_unit_name);
	fprintf(fp, "Cutoff radius over scale radius  is %e\n", cfg[ii].rc / cfg[ii].rs);
	fprintf(fp, "Cutoff  width over scale radius  is %e\n", cfg[ii].rc_width / cfg[ii].rs);
      }/* if( cfg[ii].cutoff ){ */
      fprintf(fp, "Cutoff energy of the component   is %e (= %e %s) (= %e G Mtot / rs)\n", cfg[ii].Ecut, cfg[ii].Ecut * senergy2astro, senergy_astro_unit_name, cfg[ii].Ecut / (CAST_R2D(newton) * cfg[ii].Mtot / cfg[ii].rs));
    }/* if( cfg[ii].kind != CENTRALBH ){ */
    fprintf(fp, "#############################################################################\n");

    /** estimate typical timescale */
    if( cfg[ii].kind != CENTRALBH ){
      /** estimate enclosed mass within rs */
      double Ms = 0.0;
      double Ns = 0.0;
      {
	int ll, rr;
	findIdx_rad(cfg[ii].rs, prf[ii], &ll, &rr);
	const double alpha = (cfg[ii].rs - prf[ii][ll].rad) / (prf[ii][rr].rad - prf[ii][ll].rad);

	for(int kk = 0; kk < kind; kk++){
	  Ms = (1.0 - alpha) * prf[kk][ll].enc + alpha * prf[kk][rr].enc;
	  Ns += nearbyint((double)cfg[kk].num * (Ms / cfg[kk].Mtot));
	}/* for(int kk = 0; kk < kind; kk++){ */

	Ms = (1.0 - alpha) * prf[0][ll].enc_tot + alpha * prf[0][rr].enc_tot;
      }

      /** estimate dynamical time at scale length for each component */
      const double tff = M_PI_2 * cfg[ii].rs * sqrt(cfg[ii].rs / (2.0 * CAST_R2D(newton) * Ms));
      const double t2r = tff * Ns / (32.0 * log(cfg[ii].rs / CAST_R2D(eps)));
      double trot = 0.0;
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) || (cfg[ii].kind == TBL_DISK) )
	trot = 2.0 * M_PI * cfg[ii].rs / cfg[ii].vcirc_Rd;
      fprintf(fp, "Total number of particles within the scale length is       %e\n", Ns);
      fprintf(fp, "Enclosed mass of all components within the scale length is %e (= %e %s)\n",  Ms,  Ms * mass2astro, mass_astro_unit_name);
      fprintf(fp, "Free-fall time at the scale length                      is %e (= %e %s)\n", tff, tff * time2astro, time_astro_unit_name);
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) || (cfg[ii].kind == TBL_DISK) )
	fprintf(fp, "Rotation time scale at the scale length                 is %e (= %e x tff = %e %s)\n", trot, trot / tff, trot * time2astro, time_astro_unit_name);
      fprintf(fp, "Two-body relaxation time at the scale length            is %e (= %e %s)\n", t2r, t2r * time2astro, time_astro_unit_name);
      fprintf(fp, "#############################################################################\n");
      fprintf(fp, "Snapshot interval in the unit of free-fall time           is %e\n", CAST_R2D(snapshotInterval) / tff);
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) || (cfg[ii].kind == TBL_DISK) )
	fprintf(fp, "Snapshot interval in the unit of rotation time scale      is %e\n", CAST_R2D(snapshotInterval) / trot);
      fprintf(fp, "Snapshot interval in the unit of two-body relaxation time is %e\n", CAST_R2D(snapshotInterval) / t2r);
      fprintf(fp, "#############################################################################\n");
      fprintf(fp, "Final time of the simulation in the unit of free-fall time           is %e\n", CAST_R2D(ft) / tff);
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) || (cfg[ii].kind == TBL_DISK) )
	fprintf(fp, "Final time of the simulation in the unit of rotation time scale      is %e\n", CAST_R2D(ft) / trot);
      fprintf(fp, "Final time of the simulation in the unit of two-body relaxation time is %e\n", CAST_R2D(ft) / t2r);
      fprintf(fp, "#############################################################################\n");
      fprintf(fp, "#############################################################################\n");
    }/* if( cfg[ii].kind != CENTRALBH ){ */
    fprintf(fp, "\n");
  }/* for(int ii = 0; ii < kind; ii++){ */
  fclose(fp);


  __NOTE__("%s\n", "end");
}


/**
 * @fn outputRadialProfiles
 *
 * @brief Print out radial profiles of the generated system.
 *
 * @param (kind) number of components
 * @param (prf) radial profile of the components
 * @param (cfg) physical properties of the components
 * @return (eps) value of Plummer softening length
 * @param (file) name of the simulation
 */
void outputRadialProfiles
(const int kind, profile **prf,
#ifdef  USE_HDF5_FORMAT
 profile_cfg *cfg, const real eps,
#endif//USE_HDF5_FORMAT
 char file[])
{
  __NOTE__("%s\n", "start");


  /** output fundamental profile of the particle distribution */
  real *tmp_rad;  tmp_rad = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_rad == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_rad\n");  }
  real *tmp_rho;  tmp_rho = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_rho\n");  }
  real *tmp_enc;  tmp_enc = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_enc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_enc\n");  }
  real *tmp_psi;  tmp_psi = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_psi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_psi\n");  }
  real *tmp_tff;  tmp_tff = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_tff == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_tff\n");  }
  real *tmp_t2r;  tmp_t2r = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_t2r == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_t2r\n");  }
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
  real *tmp_sig;  tmp_sig = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_sig == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_sig\n");  }
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  real *tmp_los;  tmp_los = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_los == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_los\n");  }
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  real *tmp_Sig;  tmp_Sig = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_Sig == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_Sig\n");  }
#endif//MAKE_COLUMN_DENSITY_PROFILE


  char filename[128];
#ifdef  USE_HDF5_FORMAT
  /** create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "profile");
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#ifdef  DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  uint szip_pixels_per_block = 8;
  hsize_t cdims = 128 * szip_pixels_per_block;
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
#endif//USE_HDF5_FORMAT

  for(int kk = 0; kk < kind; kk++){
#ifdef  USE_HDF5_FORMAT
    char grp[16];    sprintf(grp, "data%d", kk);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else///USE_HDF5_FORMAT
    FILE *fp;
#ifndef WRITE_IN_GALACTICS_FORMAT
    sprintf(filename, "%s/%s.profile.%d.dat", DATAFOLDER, file, kk);
    fp = fopen(filename, "wb");
#else///WRITE_IN_GALACTICS_FORMAT
    sprintf(filename, "%s/%s.profile.%d.txt", DATAFOLDER, file, kk);
    fp = fopen(filename, "w");
#endif//WRITE_IN_GALACTICS_FORMAT
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }
#endif//USE_HDF5_FORMAT

#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      tmp_rad[ii] = CAST_D2R(prf[kk][ii].rad *  length2astro);
      tmp_rho[ii] = CAST_D2R(prf[kk][ii].rho * density2astro);
      tmp_enc[ii] = CAST_D2R(prf[kk][ii].enc *    mass2astro);
      tmp_psi[ii] = CAST_D2R(prf[kk][ii].psi * senergy2astro);
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
      tmp_sig[ii] = CAST_D2R(prf[kk][ii].sigr * velocity2astro);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
      tmp_los[ii] = CAST_D2R(prf[kk][ii].slos * velocity2astro);
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
      tmp_Sig[ii] = CAST_D2R(prf[kk][ii].Sigma * col_density2astro);
#endif//MAKE_COLUMN_DENSITY_PROFILE
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */

#ifdef  USE_HDF5_FORMAT
    hsize_t dims = NRADBIN;
    dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION

    /** write radius */
    dataset = H5Dcreate(group, "rad", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rad));
    chkHDF5err(H5Dclose(dataset));
    /** write density */
    dataset = H5Dcreate(group, "rho", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rho));
    chkHDF5err(H5Dclose(dataset));
    /** write enclosed mass */
    dataset = H5Dcreate(group, "enc", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_enc));
    chkHDF5err(H5Dclose(dataset));
    /** write potential */
    dataset = H5Dcreate(group, "Psi", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_psi));
    chkHDF5err(H5Dclose(dataset));
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
    /** write 1-dimensional velocity dispersion */
    dataset = H5Dcreate(group, "sigma_r", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_sig));
    chkHDF5err(H5Dclose(dataset));
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
    /** write velocity dispersion along the line-of-sight */
    dataset = H5Dcreate(group, "sigma_los", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_los));
    chkHDF5err(H5Dclose(dataset));
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
    /** write column density */
    dataset = H5Dcreate(group, "Sigma", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_Sig));
    chkHDF5err(H5Dclose(dataset));
#endif//MAKE_COLUMN_DENSITY_PROFILE
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** write attribute data */
    /** create the data space for the attribute */
    hsize_t attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;

    /** write # of arrays */
    int nradbin = NRADBIN;
    attribute = H5Acreate(group, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nradbin));
    chkHDF5err(H5Aclose(attribute));
    /** write scale radius */
    attribute = H5Acreate(group, "rs", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg[kk].rs));
    chkHDF5err(H5Aclose(attribute));
    /** write total mass */
    attribute = H5Acreate(group, "Mtot", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg[kk].Mtot));
    chkHDF5err(H5Aclose(attribute));
    /** profile ID */
    attribute = H5Acreate(group, "profile ID", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cfg[kk].kind));
    chkHDF5err(H5Aclose(attribute));

    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
#else///USE_HDF5_FORMAT
#ifndef WRITE_IN_GALACTICS_FORMAT
    int nradbin = NRADBIN;
    bool success = true;
    success &= (fwrite(&nradbin, sizeof( int),       1, fp) ==       1);
    success &= (fwrite( tmp_rad, sizeof(real), NRADBIN, fp) == NRADBIN);
    success &= (fwrite( tmp_rho, sizeof(real), NRADBIN, fp) == NRADBIN);
    success &= (fwrite( tmp_enc, sizeof(real), NRADBIN, fp) == NRADBIN);
    success &= (fwrite( tmp_psi, sizeof(real), NRADBIN, fp) == NRADBIN);
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
    success &= (fwrite( tmp_sig, sizeof(real), NRADBIN, fp) == NRADBIN);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
    success &= (fwrite( tmp_los, sizeof(real), NRADBIN, fp) == NRADBIN);
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
    success &= (fwrite( tmp_Sig, sizeof(real), NRADBIN, fp) == NRADBIN);
#endif//MAKE_COLUMN_DENSITY_PROFILE
    if( !success ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
#else///WRITE_IN_GALACTICS_FORMAT
    fprintf(fp, "#r\trho(r)\tM(r)\tPsi(r)");
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
    fprintf(fp, "\tsigma_r");
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
    fprintf(fp, "\tsigma_los");
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
    fprintf(fp, "\tSigma");
#endif//MAKE_COLUMN_DENSITY_PROFILE
    fprintf(fp, "\n");
    for(int ii = 0; ii < NRADBIN; ii += (NRADBIN / 1024)){
      fprintf(fp, "%e\t%e\t%e\t%e", tmp_rad[ii], tmp_rho[ii], tmp_enc[ii], tmp_psi[ii]);
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
      fprintf(fp, "\t%e", tmp_sig[ii]);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
      fprintf(fp, "\t%e", tmp_los[ii]);
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
      fprintf(fp, "\t%e", tmp_Sig[ii]);
#endif//MAKE_COLUMN_DENSITY_PROFILE
      fprintf(fp, "\n");
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */
#endif//WRITE_IN_GALACTICS_FORMAT
    fclose(fp);
#endif//USE_HDF5_FORMAT
  }/* for(int kk = 0; kk < kind; kk++){ */


  /** evaluate typical timescale */
#ifdef  USE_HDF5_FORMAT
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad = prf[0][ii].rad;
    const double enc = prf[0][ii].enc_tot;

    double Ns = 0.0;
    for(int kk = 0; kk < kind; kk++)
      Ns += nearbyint((double)cfg[kk].num * (prf[kk][ii].enc / cfg[kk].Mtot));

    const double tff = M_PI_2 * rad * sqrt(rad / (2.0 * CAST_R2D(newton) * enc));
    const double t2r = fmax(tff * (double)Ns / (32.0 * log(enc / CAST_R2D(eps))), 0.0);

    tmp_rad[ii] = CAST_D2R(rad * length2astro);
    tmp_tff[ii] = CAST_D2R(tff *   time2astro);
    tmp_t2r[ii] = CAST_D2R(t2r *   time2astro);
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */

  /** write typical timescale */
  hid_t timeScale = H5Gcreate(target, "time_scale", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t dims = NRADBIN;
  dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  chkHDF5err(H5Pset_chunk(property, 1, &cdims));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION

  /** write radius */
  dataset = H5Dcreate(timeScale, "rad", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rad));
  chkHDF5err(H5Dclose(dataset));
  /** write free-fall time */
  dataset = H5Dcreate(timeScale, "free_fall_time", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_tff));
  chkHDF5err(H5Dclose(dataset));
  /** write relaxation time */
  dataset = H5Dcreate(timeScale, "relaxation_time", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_t2r));
  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));

  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  int nradbin = NRADBIN;
  hid_t attribute = H5Acreate(timeScale, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nradbin));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(timeScale));

  /** write attribute data */
  /** create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);

  /** write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &kind));
  chkHDF5err(H5Aclose(attribute));

  /** write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));

  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
  attribute = H5Acreate(target,      "length_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,     length_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,     "density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,     density_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,        "mass_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,        mass_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,     "senergy_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,     senergy_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,        "time_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,        time_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "col_density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, col_density_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,    "velocity_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,    velocity_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));

  /** close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  /** close the file */
  chkHDF5err(H5Fclose(target));
#endif//USE_HDF5_FORMAT

  free(tmp_rad);
  free(tmp_rho);
  free(tmp_enc);
  free(tmp_psi);
  free(tmp_tff);
  free(tmp_t2r);
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
  free(tmp_sig);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  free(tmp_los);
#endif//MAKE_COLUMN_DENSITY_PROFILE
#endif//MAKE_VELOCITY_DISPERSION_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  free(tmp_Sig);
#endif//MAKE_COLUMN_DENSITY_PROFILE


  __NOTE__("%s\n", "end");
}


/**
 * @fn outputDistributionFunction
 *
 * @brief Print out distribution function of the generated system.
 *
 * @param (skind) number of spherical symmetric components
 * @param (df) distribution function of the component
 * @param (file) name of the simulation
 */
void outputDistributionFunction(const int skind, dist_func **df, char file[])
{
  __NOTE__("%s\n", "start");

  char filename[256];


  /** output distribution function of the particle distribution */
  real *tmp_ene;  tmp_ene = (real *)malloc(NENEBIN * sizeof(real));  if( tmp_ene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_ene\n");  }
  real *tmp_val;  tmp_val = (real *)malloc(NENEBIN * sizeof(real));  if( tmp_val == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_val\n");  }

#ifdef  USE_HDF5_FORMAT
  /** create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "df");
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hid_t dataspace;
  hsize_t attr_dims;
  hid_t attribute;

  for(int kk = 0; kk < skind; kk++){
#pragma omp parallel for
    for(int ii = 0; ii < NENEBIN; ii++){
      tmp_ene[ii] = CAST_D2R(CAST_R2D(df[kk][ii].ene) * senergy2astro);
      tmp_val[ii] =                   df[kk][ii].val;
    }/* for(int ii = 0; ii < NENEBIN; ii++){ */

    char grp[16];    sprintf(grp,  "series%d", kk);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

#ifdef  USE_SZIP_COMPRESSION
    hid_t property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#else///USE_SZIP_COMPRESSION
    hid_t property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION

    hsize_t dims = NENEBIN;
    dataspace = H5Screate_simple(1, &dims, NULL);

#ifdef  DOUBLE_PRECISION
    hid_t hdf5_real = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
    hid_t hdf5_real = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION

    /** write energy */
    hid_t dataset;
    dataset = H5Dcreate(group, "energy", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_ene));
    chkHDF5err(H5Dclose(dataset));
    /** write energy */
    dataset = H5Dcreate(group, "DF", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_val));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION

    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** write attribute data */
    /** create the data space for the attribute */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);

    /** write # of arrays */
    int nenebin = NENEBIN;
    attribute = H5Acreate(group, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nenebin));
    chkHDF5err(H5Aclose(attribute));

    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
  }/* for(int kk = 0; kk < skind; kk++){ */

  /** write attribute data */
  /** create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);

  /** write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skind));
  chkHDF5err(H5Aclose(attribute));

  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
  attribute = H5Acreate(target, "senergy_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, senergy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));

  /** write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));

  /** close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  /** close the file */
  chkHDF5err(H5Fclose(target));
#else///USE_HDF5_FORMAT
  FILE *fp;

  for(int kk = 0; kk < skind; kk++){
    sprintf(filename, "%s/%s.df.%d.dat", DATAFOLDER, file, kk);
    fp = fopen(filename, "wb");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }

#pragma omp parallel for
    for(int ii = 0; ii < NENEBIN; ii++){
      tmp_ene[ii] = CAST_D2R(CAST_R2D(df[kk][ii].ene) * senergy2astro);
      tmp_val[ii] =                   df[kk][ii].val;
    }/* for(int ii = 0; ii < NENEBIN; ii++){ */

    int nenebin = NENEBIN;
    bool success = true;
    success &= (fwrite(&nenebin, sizeof( int),       1, fp) ==       1);
    success &= (fwrite( tmp_ene, sizeof(real), NENEBIN, fp) == NENEBIN);
    success &= (fwrite( tmp_val, sizeof(real), NENEBIN, fp) == NENEBIN);
    if( !success ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
    fclose(fp);
  }/* for(int kk = 0; kk < skind; kk++){ */
#endif//USE_HDF5_FORMAT

  free(tmp_ene);
  free(tmp_val);


  __NOTE__("%s\n", "end");
}



#ifdef  USE_HDF5_FORMAT
/**
 * @fn writeVals_at_givenRad
 *
 * @brief Write physical properties at the given radius.
 */
static inline void writeVals_at_givenRad
(const hid_t group, const double rad, char subgrp[], double * restrict attr,
 const int kind, const int skind, profile **prf,
 const int maxLev, disk_data * disk)
{
  __NOTE__("%s\n", "start");

  /** extract physical quantities */
  double alp;
  int ll, rr;
  findIdx_rad(rad, prf[0], &ll, &rr);
  alp = (rad - prf[0][ll].rad) / (prf[0][rr].rad - prf[0][ll].rad);

  double bet = 0.0;
  int lev = 0;
  int idx = 0;
  if( kind != skind )
    findIdx4nestedGrid(rad, maxLev, disk[0], &lev, &idx, &bet);


  /** open a group */
  hid_t sub = H5Gcreate(group, subgrp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


  /** write physical quantities of the representative component */
  hsize_t attr_dims;
  hid_t dataspace, attribute;
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  double astro;

  /** write the radius */
  astro = rad * length2astro;
  attribute = H5Acreate(sub, "r", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &astro));
  chkHDF5err(H5Aclose(attribute));

  /** write total mass */
  astro = ((1.0 - alp) * prf[0][ll].enc_tot + alp * prf[0][rr].enc_tot) * mass2astro;
  attribute = H5Acreate(sub, "Mtot", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &astro));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));


  /** write physical quantities of all components */
  attr_dims = kind;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);

  /** write enclosed mass */
  for(int jj = 0; jj < kind; jj++)
    attr[jj] = ((1.0 - alp) * prf[jj][ll].enc + alp * prf[jj][rr].enc) * mass2astro;
  attribute = H5Acreate(sub, "M_enc", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, attr));
  chkHDF5err(H5Aclose(attribute));

  /** write column density */
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  for(int jj = 0; jj < skind; jj++)
    attr[jj] = ((1.0 - alp) * prf[jj][ll].Sigma + alp * prf[jj][rr].Sigma) * col_density2astro;
  for(int jj = skind; jj < kind; jj++)
    attr[jj] = ((1.0 - bet) * disk[jj - skind].Sigma[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, idx)] + bet * disk[jj - skind].Sigma[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, idx + 1)]) * col_density2astro;

  attribute = H5Acreate(sub, "Sigma", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, attr));
  chkHDF5err(H5Aclose(attribute));
#endif//MAKE_COLUMN_DENSITY_PROFILE

  chkHDF5err(H5Sclose(dataspace));


  /** write one-dimensional velocity dispersion */
#ifdef  MAKE_VELOCITY_DISPERSION_PROFILE
  attr_dims = skind;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);

  /** write velocity dispersion of the r-component */
  for(int jj = 0; jj < skind; jj++)
    attr[jj] = ((1.0 - alp) * prf[jj][ll].sigr + alp * prf[jj][rr].sigr) * velocity2astro;
  attribute = H5Acreate(sub, "sigma_r", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, attr));
  chkHDF5err(H5Aclose(attribute));

  /** write velocity dispersion along the line-of-sight */
  for(int jj = 0; jj < skind; jj++)
    attr[jj] = ((1.0 - alp) * prf[jj][ll].slos + alp * prf[jj][rr].slos) * velocity2astro;
  attribute = H5Acreate(sub, "sigma_los", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, attr));
  chkHDF5err(H5Aclose(attribute));

  chkHDF5err(H5Sclose(dataspace));
#endif//MAKE_VELOCITY_DISPERSION_PROFILE


  /** close the group */
  chkHDF5err(H5Gclose(sub));


  __NOTE__("%s\n", "end");
}
#endif//USE_HDF5_FORMAT


#ifdef  USE_HDF5_FORMAT
/**
 * @fn outputRepresentativeQuantities
 *
 * @brief Print out fundamental information on the generated system.
 *
 * @param (kind) number of components
 * @param (cfg) physical properties of the components
 * @param (skind) number of spherical symmetric components
 * @param (prf) radial profile of the components
 * @param (ndisk) number of disk components
 * @param (maxLev) maximum level of nested grid
 * @param (disk) physical properties of the disk component
 * @param (file) name of the simulation
 */
void outputRepresentativeQuantities
(const int kind, profile_cfg *cfg, const int skind, profile **prf, const int ndisk, const int maxLev, disk_data * disk, char file[])
{
  __NOTE__("%s\n", "start");


  /** output fundamental physical properties of each component */
  double *tmp_attr;  tmp_attr = (double *)malloc(kind * sizeof(double));
  if( tmp_attr == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_attr\n");  }

  /** create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "property");
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t dataspace, attribute;
  hsize_t attr_dims;

  /** write physical properties as attribute */
  for(int kk = 0; kk < kind; kk++){
    char grp[16];    sprintf(grp, "grp%d", kk);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /** write physical quantities at characteristic scales */
    if( cfg[kk].kind != CENTRALBH ){
      writeVals_at_givenRad(group, cfg[kk].rs   , "rs"   , tmp_attr, kind, skind, prf, maxLev, disk);
      writeVals_at_givenRad(group, cfg[kk].rhalf, "rhalf", tmp_attr, kind, skind, prf, maxLev, disk);
      writeVals_at_givenRad(group, cfg[kk].Reff , "Reff" , tmp_attr, kind, skind, prf, maxLev, disk);
    }/* if( cfg[kk].kind != CENTRALBH ){ */

    if( kk >= skind ){
      const int diskID = kk - skind;

      attr_dims = 1;
      dataspace = H5Screate_simple(1, &attr_dims, NULL);

      double astro;
      /** write scale height */
      astro = cfg[kk].zd * length2astro;
      attribute = H5Acreate(group, "zd", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &astro));
      chkHDF5err(H5Aclose(attribute));

      /** write column density @ R = 0 */
      astro = cfg[kk].Sigma0 * col_density2astro;
      attribute = H5Acreate(group, "Sigma0", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &astro));
      chkHDF5err(H5Aclose(attribute));

      /** write circular velocity @ R = Rd */
      astro = cfg[kk].vcirc_Rd * velocity2astro;
      attribute = H5Acreate(group, "vcirc(Rd)", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &astro));
      chkHDF5err(H5Aclose(attribute));

      /** write maximum circular velocity */
      astro = cfg[kk].vcirc_max * velocity2astro;
      attribute = H5Acreate(group, "vcirc_max", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &astro));
      chkHDF5err(H5Aclose(attribute));

      /** write radius where the circular velocity is maximized */
      astro = cfg[kk].vcirc_max_R * length2astro;
      attribute = H5Acreate(group, "vcirc_max_R", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &astro));
      chkHDF5err(H5Aclose(attribute));

      /** write Toomre's Q-value @ R = Rd */
      attribute = H5Acreate(group, "Q(R = Rd)", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg[kk].toomre));
      chkHDF5err(H5Aclose(attribute));

      /** write velocity dispersion in z-direction  @ R = 0 */
      astro = disk[diskID].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, 0)] * velocity2astro;
      attribute = H5Acreate(group, "sigma_z(R = 0)", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &astro));
      chkHDF5err(H5Aclose(attribute));


      /** calculate velocity dispersion in R-direction @ R = 0 */
      const int lev = maxLev - 1;
      const int ii = 0;
#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
#ifndef USE_POTENTIAL_SCALING_SCHEME
      const double Omega = sqrt(disk[diskID]. dPhidR [INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] / disk[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);
      const double kappa = sqrt(disk[diskID].d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] + 3.0 * Omega * Omega);
#else///USE_POTENTIAL_SCALING_SCHEME
      const double Omega = sqrt(disk[diskID]. dPhidR [INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] / disk[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);
      const double kappa = sqrt(disk[diskID].d2PhidR2[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] + 3.0 * Omega * Omega);
#endif//USE_POTENTIAL_SCALING_SCHEME
      const double vcirc = disk[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] * Omega;
      const double sigmap = DISK_PERP_VDISP(disk[diskID].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], vcirc, disk[diskID].cfg->vdisp_frac);
      const double sigmaR = sigmap * 2.0 * Omega / (DBL_MIN + kappa);
#else///ENFORCE_EPICYCLIC_APPROXIMATION
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      const double sigmaR = DISK_RADIAL_VDISP(disk[diskID].cfg->vdispR0, disk[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], 1.0 / disk[diskID].zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);
#else///ENABLE_VARIABLE_SCALE_HEIGHT
      const double sigmaR = DISK_RADIAL_VDISP(disk[diskID].cfg->vdispR0, disk[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], disk[diskID].invRd);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
#endif//ENFORCE_EPICYCLIC_APPROXIMATION

      /** write velocity dispersion in R-direction  @ R = 0 */
      astro = sigmaR * velocity2astro;
      attribute = H5Acreate(group, "sigma_R(R = 0)", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &astro));
      chkHDF5err(H5Aclose(attribute));

      /** write velocity dispersion in p-direction  @ R = 0 */
#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
      astro = sigmap * velocity2astro;
      attribute = H5Acreate(group, "sigma_p(R = 0)", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &astro));
      chkHDF5err(H5Aclose(attribute));
#endif//ENFORCE_EPICYCLIC_APPROXIMATION

      chkHDF5err(H5Sclose(dataspace));
    }/* if( kk >= skind ){ */

    /** close the group */
    chkHDF5err(H5Gclose(group));
  }/* for(int kk = 0; kk < kind; kk++){ */

  /** write attribute data */
  /** create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);

  /** write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &kind));
  chkHDF5err(H5Aclose(attribute));

  /** write # of spherical components */
  attribute = H5Acreate(target, "skinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skind));
  chkHDF5err(H5Aclose(attribute));

  /** write # of disk components */
  attribute = H5Acreate(target, "ndisk", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &ndisk));
  chkHDF5err(H5Aclose(attribute));

  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
  attribute = H5Acreate(target,      "length_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,      length_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,        "mass_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,        mass_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,    "velocity_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,    velocity_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "col_density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, col_density_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));

  /** close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  /** close the file */
  chkHDF5err(H5Fclose(target));

  free(tmp_attr);


  __NOTE__("%s\n", "end");
}
#endif//USE_HDF5_FORMAT


/**
 * @fn findIdx_enc
 *
 * @brief Find a data element in the given array corresponding to the given enclosed mass.
 *
 * @param (enc) enclosed mass
 * @param (prf) radial profile of the component
 * @return (ll) the corresponding lower index
 * @return (rr) the corresponding upper index
 */
static inline void findIdx_enc(const double enc, profile * restrict prf, int * restrict ll, int * restrict rr)
{
  bool bisection = true;
  *ll =           0;
  *rr = NRADBIN - 1;

  if( bisection == true )    if( fabs(prf[*ll].enc - enc) / enc < DBL_EPSILON ){      bisection = false;      *rr = (*ll) + 1;    }
  if( bisection == true )    if( fabs(prf[*rr].enc - enc) / enc < DBL_EPSILON ){      bisection = false;      *ll = (*rr) - 1;    }

  while( bisection ){
    const uint cc = ((uint)(*ll) + (uint)(*rr)) >> 1;

    if( (prf[cc].enc - enc) * (prf[*ll].enc - enc) <= 0.0 )      *rr = (int)cc;
    else                                                         *ll = (int)cc;

    if( (1 + (*ll)) == (*rr) )
      break;
  }/* while( bisection ){ */
}


/**
 * @fn evaluateObservables
 *
 * @brief Evaluate observable quantities of spherical component(s)
 *
 * @param (kind) number of components
 * @param (skind) number of spherical symmetric components
 * @param (cfg) physical properties of the components
 * @param (prf) radial profile of the components
 *
 * @sa findIdx_enc
 * @sa genCubicSpline1D
 * @sa getCubicSplineIntegral1D
 */
void evaluateObservables(const int kind, const int skind, profile_cfg *cfg, profile **prf)
{
  __NOTE__("%s\n", "start");


  /** evaluate half-mass radius r_{1/2}: the radius containing half the total mass */
  for(int kk = 0; kk < kind; kk++)
    if( cfg[kk].kind != CENTRALBH ){
      const double Mhalf = 0.5 * cfg[kk].Mtot;
      int ll, rr;
      findIdx_enc(Mhalf, prf[kk], &ll, &rr);

      const double alpha = (Mhalf - prf[kk][ll].enc) / (prf[kk][rr].enc - prf[kk][ll].enc);

      cfg[kk].rhalf = (1.0 - alpha) * prf[kk][ll].rad + alpha * prf[kk][rr].rad;
    }/* if( cfg[ii].kind != CENTRALBH ){ */


#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  /** evaluate effective radius R_{eff}: the radius containing half the total mass in the projected plane */
#pragma omp parallel for
  for(int kk = 0; kk < skind; kk++)
    if( cfg[kk].kind != CENTRALBH ){
      const double Mhalf = 0.5 * 0.5 * cfg[kk].Mtot;

      const int num = NRADBIN / SKIP_INTERVAL_FOR_EFFECTIVE_RADIUS;

      double *xx;      xx = (double *)malloc(num * sizeof(double));
      double *ff;      ff = (double *)malloc(num * sizeof(double));
      double *bp;      bp = (double *)malloc(num * sizeof(double));
      double *f2;      f2 = (double *)malloc(num * sizeof(double));

      if( xx == NULL ){	__KILL__(stderr, "ERROR: failure to allocate xx\n");      }
      if( ff == NULL ){	__KILL__(stderr, "ERROR: failure to allocate ff\n");      }
      if( bp == NULL ){	__KILL__(stderr, "ERROR: failure to allocate bp\n");      }
      if( f2 == NULL ){	__KILL__(stderr, "ERROR: failure to allocate f2\n");      }

      for(int ii = 0; ii < num; ii++){
	const int jj = ii * SKIP_INTERVAL_FOR_EFFECTIVE_RADIUS;
	xx[ii] = prf[kk][jj].rad;/**< R */
	ff[ii] = 2.0 * M_PI * prf[kk][jj].rad * prf[kk][jj].Sigma;/**< 2 \pi R \Sigma(R) */
      }/* for(int ii = 0; ii < num; ii++){ */

      genCubicSpline1D(num, xx, ff, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, f2);

      double R0 = 0.0;
      double M0 = 0.0;
      for(int ii = 0; ii < num; ii++){
	const double R1 = xx[ii];
	const double M1 = M0 + getCubicSplineIntegral1D(ii, xx, ff, f2);

	if( M1 >= Mhalf ){
	  const double alpha = (Mhalf - M0) / (M1 - M0);
	  cfg[kk].Reff = (1.0 - alpha) * R0 + alpha * R1;
	  break;
	}/* if( M1 >= Mhalf ){ */

	R0 = R1;
	M0 = M1;
      }/* for(int ii = 0; ii < num; ii++){ */

      free(xx);      free(ff);      free(bp);      free(f2);
    }/* if( cfg[kk].kind != CENTRALBH ){ */
#endif//MAKE_COLUMN_DENSITY_PROFILE


  __NOTE__("%s\n", "end");
}


/**
 * @fn evaluateDiskProperties
 *
 * @brief Evaluate physical properties of the disk component.
 *
 * @param (disk_info) physical quantities of the disk component
 * @param (diskID) index of the disk component
 * @param (lev) nested grid level
 * @param (ihead) head index
 * @param (itail) tail index
 * @return (_vcirc) profile of circular velocity
 * @return (_sigmap) profile of velocity dispersion in the azimuthal direction
 * @return (_sigmaR) profile of velocity dispersion in the horizontal direction
 * @return (_kappa) profile of epicycle frequency
 * @return (_Omega) profile of circular frequency
 * @return (_toomre) profile of Toomre's Q-value
 * @return (_lambda) profile of critical wavelength in the Toomre's analysis (wavelength greater than this is unstable)
 */
static void evaluateDiskProperties
(disk_data *disk_info, const int diskID, const int lev, const int ihead, const int itail,
 real * restrict _vcirc, real * restrict _sigmap, real * restrict _sigmaR, real * restrict _kappa, real * restrict _Omega, real * restrict _toomre, real * restrict _lambda)
{
  __NOTE__("%s\n", "start");


#ifndef ENFORCE_EPICYCLIC_APPROXIMATION
  const double invRd = 1.0 / disk_info[diskID].cfg->rs;
#endif//ENFORCE_EPICYCLIC_APPROXIMATION
  bool passed = disk_info[diskID].cfg->passed;

  for(int ii = ihead; ii < itail + 1; ii++){
    /** evaluate epicyclic quantities and circular speed */
#ifndef USE_POTENTIAL_SCALING_SCHEME
    const double  dPhidR  = disk_info[diskID]. dPhidR [INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)];
    const double d2PhidR2 = disk_info[diskID].d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)];
#else///USE_POTENTIAL_SCALING_SCHEME
    const double  dPhidR  = disk_info[diskID]. dPhidR [INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
    const double d2PhidR2 = disk_info[diskID].d2PhidR2[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
#endif//USE_POTENTIAL_SCALING_SCHEME
    const double Omega = sqrt( dPhidR  / disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);
    const double kappa = sqrt(d2PhidR2 + 3.0 * Omega * Omega);
    const double vcirc = disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] * Omega;

    /** evaluate Toomre's Q-value */
#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
    const double sigmap = DISK_PERP_VDISP(disk_info[diskID].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], vcirc, disk_info[diskID].cfg->vdisp_frac);
    const double sigmaR = sigmap * 2.0 * Omega / (DBL_MIN + kappa);
#else///ENFORCE_EPICYCLIC_APPROXIMATION
    const double sigmaR = DISK_RADIAL_VDISP(disk_info[diskID].cfg->vdispR0, disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], invRd);
    const double sigmap = sigmaR * kappa * 0.5 / (DBL_MIN + Omega);
#endif//ENFORCE_EPICYCLIC_APPROXIMATION
    const double Sigma = disk_info[diskID].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
    const double toomre = sigmaR * kappa / (DBL_MIN + 3.36 * CAST_R2D(newton) * Sigma);
    const double lambda = 4.0 * M_PI * M_PI * CAST_R2D(newton) * Sigma / (DBL_MIN + kappa * kappa);

    /** find the maximum circular speed */
    if( vcirc > disk_info[diskID].cfg->vcirc_max ){
      disk_info[diskID].cfg->vcirc_max   = vcirc;
      disk_info[diskID].cfg->vcirc_max_R = disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
    }/* if( vcirc > disk_info[diskID].cfg->vcirc_max ){ */

    /** find the minimum Toomre's Q-value */
    if( toomre < disk_info[diskID].cfg->Qmin0 ){
      disk_info[diskID].cfg->Qmin0  = toomre;
      disk_info[diskID].cfg->qminR0 = disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
    }/* if( toomre < disk_info[diskID].cfg->Qmin0 ){ */

#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
    if( (toomre < disk_info[diskID].cfg->Qmin1) && ((1.0 + DBL_EPSILON) * disk_info[diskID].zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] >= disk_info[diskID].cfg->zd ) ){
      disk_info[diskID].cfg->Qmin1  = toomre;
      disk_info[diskID].cfg->qminR1 = disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
    }
#endif//ENABLE_VARIABLE_SCALE_HEIGHT

#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
    if( (toomre < disk_info[diskID].cfg->Qmin2) && ((1.0 + DBL_EPSILON) * sigmap >= disk_info[diskID].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]) ){
      disk_info[diskID].cfg->Qmin2  = toomre;
      disk_info[diskID].cfg->qminR2 = disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
    }
#endif//ENFORCE_EPICYCLIC_APPROXIMATION

    /** find the circular speed and Toomre's Q-value at the scale radius */
    if( !passed ){
      disk_info[diskID].cfg->vcirc_Rd = vcirc;
      disk_info[diskID].cfg->toomre   = toomre;
      if( disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] > disk_info[diskID].cfg->rs )
	passed = true;
    }/* if( !passed ){ */


    /** memorize calculated values */
    const int jj = ii - ihead;
    _vcirc [jj] = CAST_D2R(vcirc  * velocity2astro);    /* if( isinf(_vcirc [jj]) == 1 )	    _vcirc [jj] = REAL_MAX; */
    _sigmap[jj] = CAST_D2R(sigmap * velocity2astro);	/* if( isinf(_sigmap[jj]) == 1 )      _sigmap[jj] = REAL_MAX; */
    _sigmaR[jj] = CAST_D2R(sigmaR * velocity2astro);	if( isinf(_sigmaR[jj]) == 1 )	   _sigmaR[jj] = REAL_MAX;
    _kappa [jj] = CAST_D2R(kappa  /     time2astro);    /* if( isinf(_kappa [jj]) == 1 )	    _kappa [jj] = REAL_MAX; */
    _Omega [jj] = CAST_D2R(Omega  /     time2astro);    /* if( isinf(_Omega [jj]) == 1 )	    _Omega [jj] = REAL_MAX; */
    _lambda[jj] = CAST_D2R(lambda *   length2astro);	if( isinf(_lambda[jj]) == 1 )	   _lambda[jj] = REAL_MAX;
    _toomre[jj] = CAST_D2R(toomre                 );	if( isinf(_toomre[jj]) == 1 )	   _toomre[jj] = REAL_MAX;
  }/* for(int ii = ihead; ii < itail + 1; ii++){ */

  disk_info[diskID].cfg->passed = passed;


  __NOTE__("%s\n", "end");
}


/**
 * @fn writeDiskData
 *
 * @brief Print out fundamental information on the disk component(s).
 *
 * @param (file) name of the simulation
 * @param (ndisk) number of disk components
 * @param (maxLev) maximum level of nested grid
 * @return (disk) physical properties of the disk component
 *
 * @sa evaluateDiskProperties
 */
void writeDiskData(char *file, const int ndisk, const int maxLev, disk_data *disk)
{
  __NOTE__("%s\n", "start");


  /** arrays for 2D-plot */
#ifdef  USE_HDF5_FORMAT
  real *node_RR;  node_RR = (real *)malloc((NDISKBIN_HOR + 1) * sizeof(real));  if( node_RR == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node_RR\n");  }
  real *node_zz;  node_zz = (real *)malloc((NDISKBIN_VER + 1) * sizeof(real));  if( node_zz == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node_zz\n");  }
#endif//USE_HDF5_FORMAT
  real *tmp_rho;  tmp_rho = (real *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(real));  if( tmp_rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_rho\n");  }
  real *rhoFrac;  rhoFrac = (real *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(real));  if( rhoFrac == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rhoFrac\n");  }
  real *tmp_Phi;  tmp_Phi = (real *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(real));  if( tmp_Phi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_Phi\n");  }

  /** arrays for 1D-plot */
  real *tmp_ver;  tmp_ver = (real *)malloc((NDISKBIN_VER >> 1) * (maxLev + 1) * sizeof(real));  if( tmp_ver == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_ver\n");  }
  real *tmp_hor;  tmp_hor = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( tmp_hor == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_hor\n");  }
  real *tmp_sig;  tmp_sig = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( tmp_sig == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_sig\n");  }
  real *tmp_Sig;  tmp_Sig = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( tmp_Sig == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_Sig\n");  }
  real *_vcirc ;  _vcirc  = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( _vcirc  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _vcirc\n" );  }
  real *_sigmap;  _sigmap = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( _sigmap == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _sigmap\n");  }
  real *_sigmaR;  _sigmaR = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( _sigmaR == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _sigmaR\n");  }
  real *_kappa ;  _kappa  = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( _kappa  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _kappa\n" );  }
  real *_Omega ;  _Omega  = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( _Omega  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _Omega\n" );  }
  real *_lambda;  _lambda = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( _lambda == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _lambda\n");  }
  real *_toomre;  _toomre = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( _toomre == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _toomre\n");  }
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  real *tmp__zd;  tmp__zd = (real *)malloc((NDISKBIN_HOR >> 1) * (maxLev + 1) * sizeof(real));  if( tmp__zd == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp__zd\n");  }
#endif//ENABLE_VARIABLE_SCALE_HEIGHT

  if( ndisk > 1 )
    for(int lev = 0; lev < maxLev; lev++)
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR * NDISKBIN_VER; ii++){
	disk[0].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)] = (*disk[0].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)];
	for(int jj = 1; jj < ndisk; jj++)
	  disk[0].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)] += (*disk[jj].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)];
	disk[0].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)] = 1.0 / (DBL_MIN + disk[0].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)]);
      }/* for(int ii = 0; ii < NDISKBIN_HOR * NDISKBIN_VER; ii++){ */


#ifdef  USE_HDF5_FORMAT
  char filename[128];
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "disk");
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#ifdef  DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION

  /** preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_SZIP_COMPRESSION
  /** compression using szip */
  uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  uint szip_pixels_per_block = 8;
  hsize_t cdims[2];
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION

  for(int ii = 0; ii < ndisk; ii++){
    /** output in HDF5 format */
    char grp[16];    sprintf(grp, "data%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /** write two-dimensional data */
    for(int lev = 0; lev < maxLev; lev++){
      char subgrp[16];      sprintf(subgrp, "patch%d", lev);
      hid_t patch = H5Gcreate(group, subgrp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      /** write data */
      if( ndisk > 1 )
#pragma omp parallel for
	for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++)
	  rhoFrac[jj] = CAST_D2R((*disk[ii].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * disk[ii].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)]);
#pragma omp parallel for
      for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++){
	tmp_rho[jj] = CAST_D2R((*disk[ii].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * density2astro);
	tmp_Phi[jj] = CAST_D2R(	 disk[ii].pot [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * senergy2astro);
      }/* for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++){ */

      hsize_t dims[2] = {NDISKBIN_HOR, NDISKBIN_VER};
      dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
      property = H5Pcreate(H5P_DATASET_CREATE);
      cdims[0] = (dims[0] <  128                          ) ? (dims[0]) : (128);
      cdims[1] = (dims[1] < (128 * szip_pixels_per_block) ) ? (dims[1]) : (128 * szip_pixels_per_block);
      chkHDF5err(H5Pset_chunk(property, 2, cdims));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION

      /** write density */
      dataset = H5Dcreate(patch, "rho", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rho));
      chkHDF5err(H5Dclose(dataset));
      /** write potential */
      dataset = H5Dcreate(patch, "Phi", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_Phi));
      chkHDF5err(H5Dclose(dataset));
      /** write fraction */
      if( ndisk > 1 ){
	dataset = H5Dcreate(patch, "fraction", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
	chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rhoFrac));
	chkHDF5err(H5Dclose(dataset));
      }/* if( ndisk > 1 ){ */
#ifdef  USE_SZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
      /** close the dataspace */
      chkHDF5err(H5Sclose(dataspace));

      /** prepare horizontal and vertical axes */
      double width[2];
      width[0] = ldexp(disk[ii].hh, -lev) * length2astro;
      width[1] = ldexp(disk[ii].hh, -lev) * length2astro;

      /** note: VisIt requires position of edges of the grid */
#pragma omp parallel
      {
#pragma omp for nowait
	for(int jj = 0; jj < NDISKBIN_HOR; jj++)
	  tmp_hor[jj] = CAST_D2R(disk[ii].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] * length2astro);
#pragma omp for nowait
	for(int jj = 0; jj < NDISKBIN_VER; jj++)
	  tmp_ver[jj] = CAST_D2R(disk[ii].ver[INDEX2D(maxLev, NDISKBIN_VER, lev, jj)] * length2astro);
#pragma omp for nowait
	for(int jj = 0; jj < NDISKBIN_HOR + 1; jj++)
	  node_RR[jj] = CAST_D2R(disk[ii].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, jj)] * length2astro);
#pragma omp for nowait
	for(int jj = 0; jj < NDISKBIN_VER + 1; jj++)
	  node_zz[jj] = CAST_D2R(disk[ii].node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev, jj)] * length2astro);
      }

      /** write horizontal axis */
      dataspace = H5Screate_simple(1, &dims[0], NULL);
#ifdef  USE_SZIP_COMPRESSION
      property = H5Pcreate(H5P_DATASET_CREATE);
      cdims[0] = (dims[0] < (128 * szip_pixels_per_block) ) ? (dims[0]) : (128 * szip_pixels_per_block);
      chkHDF5err(H5Pset_chunk(property, 1, &cdims[0]));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
      /** write horizontal position */
      dataset = H5Dcreate(patch, "radius", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_hor));
      chkHDF5err(H5Dclose(dataset));
      /** close the dataspace */
      chkHDF5err(H5Sclose(dataspace));

      dims[0]++;
      dataspace = H5Screate_simple(1, &dims[0], NULL);
#ifdef  USE_SZIP_COMPRESSION
      property = H5Pcreate(H5P_DATASET_CREATE);
      cdims[0] = (dims[0] < (128 * szip_pixels_per_block) ) ? (dims[0]) : (128 * szip_pixels_per_block);
      chkHDF5err(H5Pset_chunk(property, 1, &cdims[0]));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION

      /** write horizontal position */
      dataset = H5Dcreate(patch, "node_RR", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, node_RR));
      chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
      /** close the dataspace */
      chkHDF5err(H5Sclose(dataspace));

      /** write vertical axis */
      dataspace = H5Screate_simple(1, &dims[1], NULL);
#ifdef  USE_SZIP_COMPRESSION
      property = H5Pcreate(H5P_DATASET_CREATE);
      chkHDF5err(H5Pset_chunk(property, 1, &cdims[1]));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION

      /** write vertical position */
      dataset = H5Dcreate(patch, "height", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_ver));
      chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
      /** close the dataspace */
      chkHDF5err(H5Sclose(dataspace));

      dims[1]++;
      dataspace = H5Screate_simple(1, &dims[1], NULL);
#ifdef  USE_SZIP_COMPRESSION
      property = H5Pcreate(H5P_DATASET_CREATE);
      cdims[1] = (dims[1] < (128 * szip_pixels_per_block) ) ? (dims[1]) : (128 * szip_pixels_per_block);
      chkHDF5err(H5Pset_chunk(property, 1, &cdims[1]));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION

      /** write vertical position */
      dataset = H5Dcreate(patch, "node_zz", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, node_zz));
      chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
      /** close the dataspace */
      chkHDF5err(H5Sclose(dataspace));

      /** write attribute data */
      /** create the data space for the attribute */
      hsize_t attr_dims = 2;
      dataspace = H5Screate_simple(1, &attr_dims, NULL);
      hid_t attribute;

      /** write # of arrays */
      int narray[2] = {NDISKBIN_HOR, NDISKBIN_VER};
      attribute = H5Acreate(patch, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, narray));
      chkHDF5err(H5Aclose(attribute));
      /** write the head index which corresponding to the regular grid having the same resolution with this level (i.e., (NDISKBIN_HOR << lev) * (NDISKBIN_VER << lev) grid points) */
      narray[0] = 0;
      narray[1] = 0;
      attribute = H5Acreate(patch, "head", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, narray));
      chkHDF5err(H5Aclose(attribute));
      /** write mesh sizes */
      attribute = H5Acreate(patch, "width", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, width));
      chkHDF5err(H5Aclose(attribute));

      chkHDF5err(H5Sclose(dataspace));
      attr_dims = 1;
      dataspace = H5Screate_simple(1, &attr_dims, NULL);

      /** Nested level */
      attribute = H5Acreate(patch, "level", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &lev));
      chkHDF5err(H5Aclose(attribute));
      /** Parent patch */
      int parent = (lev > 0) ? (lev - 1) : 0;
      attribute = H5Acreate(patch, "parent", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &parent));
      chkHDF5err(H5Aclose(attribute));
      /** close the dataspace */
      chkHDF5err(H5Sclose(dataspace));

      chkHDF5err(H5Gclose(patch));
    }/* for(int lev = 0; lev < maxLev; lev++){ */


    /** prepare one-dimensional data */
    /** data preparation in the finest grid */
    evaluateDiskProperties(disk, ii, maxLev - 1, 0, NDISKBIN_HOR - 1, _vcirc, _sigmap, _sigmaR, _kappa, _Omega, _toomre, _lambda);
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_HOR; jj++){
      tmp_hor[jj] = CAST_D2R(disk[ii].hor   [INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] *      length2astro);
      tmp_sig[jj] = CAST_D2R(disk[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] *    velocity2astro);
      tmp_Sig[jj] = CAST_D2R(disk[ii].Sigma [INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] * col_density2astro);
#ifdef	ENABLE_VARIABLE_SCALE_HEIGHT
      tmp__zd[jj] = CAST_D2R(disk[ii].zd    [INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] *      length2astro);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    }/* for(int jj = 0; jj < NDISKBIN_HOR; jj++){ */
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
      tmp_ver[jj] = CAST_D2R(disk[ii].ver   [INDEX2D(maxLev, NDISKBIN_VER, maxLev - 1, jj)] *      length2astro);

    /** data preparation in coarser grids */
    for(int lev = maxLev - 2; lev >= 0; lev--)
      evaluateDiskProperties(disk, ii, lev, NDISKBIN_HOR >> 1, NDISKBIN_HOR - 1,
			     &(_vcirc [(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_sigmap[(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_sigmaR[(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_kappa [(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_Omega [(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_toomre[(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_lambda[(NDISKBIN_HOR >> 1) * (maxLev - lev)]));
#pragma omp parallel
    for(int lev = maxLev - 2; lev >= 0; lev--){
#pragma omp for nowait
      for(int jj = 0; jj < (NDISKBIN_HOR >> 1); jj++){
	tmp_hor[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = CAST_D2R(disk[ii].hor   [INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] *      length2astro);
	tmp_sig[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = CAST_D2R(disk[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] *    velocity2astro);
	tmp_Sig[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = CAST_D2R(disk[ii].Sigma [INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] * col_density2astro);
#ifdef	ENABLE_VARIABLE_SCALE_HEIGHT
	tmp__zd[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = CAST_D2R(disk[ii].zd    [INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] *      length2astro);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      }/* for(int jj = 0; jj < NDISKBIN_HOR; jj++){ */
#pragma omp for nowait
      for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++)
	tmp_ver[(NDISKBIN_VER >> 1) * (maxLev - lev) + jj] = CAST_D2R(disk[ii].ver   [INDEX2D(maxLev, NDISKBIN_VER, lev, (NDISKBIN_VER >> 1) + jj)] *      length2astro);
    }/* for(int lev = maxLev - 2; lev >= 0; lev--){ */

    /** write one-dimensional data */
    char subgrp[16];
    sprintf(subgrp, "1D_data");
    hid_t patch = H5Gcreate(group, subgrp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims[2] = {(NDISKBIN_HOR >> 1) * (maxLev + 1), (NDISKBIN_VER >> 1) * (maxLev + 1)};

    /** write horizontal variables */
    dataspace = H5Screate_simple(1, &dims[0], NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    cdims[0] = (dims[0] < (128 * szip_pixels_per_block) ) ? (dims[0]) : (128 * szip_pixels_per_block);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims[0]));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION

    /** write horizontal position */
    dataset = H5Dcreate(patch, "radius", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_hor));
    chkHDF5err(H5Dclose(dataset));
    /** write velocity dispersion in z-direction */
    dataset = H5Dcreate(patch, "sigmaz", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_sig));
    chkHDF5err(H5Dclose(dataset));
    /** write column density */
    dataset = H5Dcreate(patch, "Sigma", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_Sig));
    chkHDF5err(H5Dclose(dataset));
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
    /** write scale height */
    dataset = H5Dcreate(patch, "zd", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp__zd));
    chkHDF5err(H5Dclose(dataset));
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    /** write velocity dispersion in R-direction */
    dataset = H5Dcreate(patch, "sigmaR", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _sigmaR));
    chkHDF5err(H5Dclose(dataset));
    /** write velocity dispersion in tangential direction */
    dataset = H5Dcreate(patch, "sigmap", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _sigmap));
    chkHDF5err(H5Dclose(dataset));
    /** write circular velocity */
    dataset = H5Dcreate(patch, "vcirc", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _vcirc));
    chkHDF5err(H5Dclose(dataset));
    /** write kappa */
    dataset = H5Dcreate(patch, "kappa", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _kappa));
    chkHDF5err(H5Dclose(dataset));
    /** write Omega */
    dataset = H5Dcreate(patch, "Omega", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _Omega));
    chkHDF5err(H5Dclose(dataset));
    /** write lambda_critical */
    dataset = H5Dcreate(patch, "lambda", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _lambda));
    chkHDF5err(H5Dclose(dataset));
    /** write Toomre's Q-value */
    dataset = H5Dcreate(patch, "Toomre's Q", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _toomre));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** write vertical variables */
    dataspace = H5Screate_simple(1, &dims[1], NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims[1]));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION

    /** write vertical position */
    dataset = H5Dcreate(patch, "height", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_ver));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** write attribute data */
    /** create the data space for the attribute */
    hsize_t attr_dims = 2;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    /** write # of arrays */
    attribute = H5Acreate(patch, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, dims));
    chkHDF5err(H5Aclose(attribute));
    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(patch));


    /** write attribute data */
    /** create the data space for the attribute */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);

    /** write scale radius */
    attribute = H5Acreate(group, "Rs", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &disk[ii].cfg->rs));
    chkHDF5err(H5Aclose(attribute));
    /** write scale height */
    attribute = H5Acreate(group, "zd", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &disk[ii].cfg->zd));
    chkHDF5err(H5Aclose(attribute));
    /** write total mass */
    attribute = H5Acreate(group, "Mtot", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &disk[ii].cfg->Mtot));
    chkHDF5err(H5Aclose(attribute));
    /** profile ID */
    attribute = H5Acreate(group, "profile ID", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &disk[ii].cfg->kind));
    chkHDF5err(H5Aclose(attribute));
    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    chkHDF5err(H5Gclose(group));
  }/* for(int ii = 0; ii < ndisk; ii++){ */

  /** write attribute data */
  /** create the data space for the attribute */
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  /** write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &ndisk));
  chkHDF5err(H5Aclose(attribute));
  /** write maximum # of nested level */
  attribute = H5Acreate(target, "maxLev", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &maxLev));
  chkHDF5err(H5Aclose(attribute));

  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
  attribute = H5Acreate(target,      "length_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,      length_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,     "density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,     density_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,     "senergy_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,     senergy_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,    "velocity_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,    velocity_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "col_density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, col_density_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,        "time_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,        time_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));

  /** write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));

  /** close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  /** close the file */
  chkHDF5err(H5Fclose(target));
#else///USE_HDF5_FORMAT
  for(int ii = 0; ii < ndisk; ii++){
    FILE *fp;
    char filename[256];
    sprintf(filename, "%s/%s.diskdat.%d.dat", DATAFOLDER, file, ii);
    fp = fopen(filename, "wb");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }

    bool success = true;
    const int nhorbin = NDISKBIN_HOR;
    const int nverbin = NDISKBIN_VER;
    success &= (fwrite(&nhorbin, sizeof(int), 1, fp) == 1);
    success &= (fwrite(&nverbin, sizeof(int), 1, fp) == 1);
    success &= (fwrite(& maxLev, sizeof(int), 1, fp) == 1);
    success &= (fwrite(&disk[ii].cfg->rs, sizeof(double), 1, fp) == 1);
    success &= (fwrite(&disk[ii].cfg->zd, sizeof(double), 1, fp) == 1);

    for(int lev = 0; lev < maxLev; lev++){
      if( ndisk > 1 )
#pragma omp parallel for
	for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++)
	  rhoFrac[jj] = CAST_D2R((*disk[ii].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * disk[ii].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)]);
      evaluateDiskProperties(disk, ii, lev, 0, NDISKBIN_HOR - 1, _vcirc, _sigmap, _sigmaR, _kappa, _Omega, _toomre, _lambda);

      for(int jj = 0; jj < nhorbin; jj++){
	tmp_hor[jj] = CAST_D2R(disk[ii].hor   [INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] *      length2astro);
	tmp_sig[jj] = CAST_D2R(disk[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] *    velocity2astro);
	tmp_Sig[jj] = CAST_D2R(disk[ii].Sigma [INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] * col_density2astro);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
	tmp__zd[jj] = CAST_D2R(disk[ii].zd    [INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] *      length2astro);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      }/* for(int jj = 0; jj < nhorbin; jj++){ */
      for(int jj = 0; jj < nverbin; jj++)
	tmp_ver[jj] = CAST_D2R(disk[ii].ver[INDEX2D(maxLev, NDISKBIN_VER, lev, jj)] * length2astro);
      for(int jj = 0; jj < nhorbin * nverbin; jj++){
	tmp_rho[jj] = CAST_D2R((*disk[ii].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * density2astro);
	tmp_Phi[jj] = CAST_D2R(  disk[ii].pot [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * senergy2astro);
      }/* for(int jj = 0; jj < nhorbin * nverbin; jj++){ */

      success &= (fwrite(tmp_hor, sizeof(real), nhorbin          , fp) == (size_t)nhorbin                  );
      success &= (fwrite(tmp_ver, sizeof(real),           nverbin, fp) ==                   (size_t)nverbin);
      success &= (fwrite(tmp_rho, sizeof(real), nhorbin * nverbin, fp) == (size_t)nhorbin * (size_t)nverbin);
      success &= (fwrite(tmp_Phi, sizeof(real), nhorbin * nverbin, fp) == (size_t)nhorbin * (size_t)nverbin);
      success &= (fwrite(tmp_sig, sizeof(real), nhorbin          , fp) == (size_t)nhorbin                  );
      success &= (fwrite(tmp_Sig, sizeof(real), nhorbin          , fp) == (size_t)nhorbin                  );
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      success &= (fwrite(tmp__zd, sizeof(real), nhorbin          , fp) == (size_t)nhorbin                  );
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    }/* for(int lev = 0; lev < maxLev; lev++){ */

    fclose(fp);
    if( success != true ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }


    sprintf(filename, "%s/%s.disk_hor.%d.txt", DATAFOLDER, file, ii);
    fp = fopen(filename, "w");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }

    /** prepare one-dimensional data */
    /** data preparation in the finest grid */
    evaluateDiskProperties(disk, ii, maxLev - 1, 0, NDISKBIN_HOR - 1, _vcirc, _sigmap, _sigmaR, _kappa, _Omega, _toomre, _lambda);
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_HOR; jj++){
      tmp_hor[jj] = CAST_D2R(disk[ii].hor   [INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] *      length2astro);
      tmp_sig[jj] = CAST_D2R(disk[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] *    velocity2astro);
      tmp_Sig[jj] = CAST_D2R(disk[ii].Sigma [INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] * col_density2astro);
#ifdef	ENABLE_VARIABLE_SCALE_HEIGHT
      tmp__zd[jj] = CAST_D2R(disk[ii].zd    [INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] *      length2astro);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    }/* for(int jj = 0; jj < NDISKBIN_HOR; jj++){ */

    /** data preparation in coarser grids */
    for(int lev = maxLev - 2; lev >= 0; lev--)
      evaluateDiskProperties(disk, ii, lev, NDISKBIN_HOR >> 1, NDISKBIN_HOR - 1,
			     &(_vcirc [(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_sigmap[(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_sigmaR[(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_kappa [(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_Omega [(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_toomre[(NDISKBIN_HOR >> 1) * (maxLev - lev)]),
			     &(_lambda[(NDISKBIN_HOR >> 1) * (maxLev - lev)]));
#pragma omp parallel
    for(int lev = maxLev - 2; lev >= 0; lev--){
#pragma omp for nowait
      for(int jj = 0; jj < (NDISKBIN_HOR >> 1); jj++){
	tmp_hor[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = CAST_D2R(disk[ii].hor   [INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] *      length2astro);
	tmp_sig[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = CAST_D2R(disk[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] *    velocity2astro);
	tmp_Sig[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = CAST_D2R(disk[ii].Sigma [INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] * col_density2astro);
#ifdef	ENABLE_VARIABLE_SCALE_HEIGHT
	tmp__zd[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = CAST_D2R(disk[ii].zd    [INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] *      length2astro);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      }/* for(int jj = 0; jj < NDISKBIN_HOR; jj++){ */
    }/* for(int lev = maxLev - 2; lev >= 0; lev--){ */


    fprintf(fp, "#R\tSigma\tsigmaR\tsigmap\tsigmaz\tvcirc\tkappa\tOmega\tlambda\tQ-val");
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
    fprintf(fp, "\tzd");
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    fprintf(fp, "\n");

    for(int jj = 0; jj < (NDISKBIN_HOR >> 1) * (maxLev + 1); jj++){
      fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e", tmp_hor[jj], tmp_Sig[jj], _sigmaR[jj], _sigmap[jj], tmp_sig[jj], _vcirc[jj], _kappa[jj], _Omega[jj], _lambda[jj], _toomre[jj]);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      fprintf(fp, "\t%e", tmp__zd[jj]);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    fprintf(fp, "\n");
    }/* for(int jj = 0; jj < (NDISKBIN_HOR >> 1) * (maxLev + 1); jj++){ */

    fclose(fp);
  }/* for(int ii = 0; ii < ndisk; ii++){ */
#endif//USE_HDF5_FORMAT


  free(tmp_hor);  free(tmp_ver);
  free(tmp_sig);  free(tmp_Sig);
  free(_vcirc);  free(_sigmap);  free(_sigmaR);
  free(_kappa);  free(_Omega);  free(_lambda);  free(_toomre);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  free(tmp__zd);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  free(tmp_rho);  free(tmp_Phi);  free(rhoFrac);

  __NOTE__("%s\n", "end");
}
