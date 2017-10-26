/**
 * @file uniformsphere.c
 *
 * @brief Source code for generating uniform sphere with Gaussian velocity dispersion
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2017/10/26 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

/**< conversion from physical unit to computational unit must be performed internally */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "myutil.h"
#include "constants.h"
#include "timer.h"
#include "name.h"
#include "rand.h"

#include "../misc/structure.h"
#include "../misc/allocate.h"

#ifndef RUN_WITHOUT_GOTHIC
#include "../misc/tune.h"
#include "../misc/brent.h"
#endif//RUN_WITHOUT_GOTHIC

#include "../file/io.h"


extern const double mass_astro2com, velocity_astro2com, length_astro2com, time_astro2com;
extern const real newton;


/**
 * @fn isotropicDistribution
 *
 * @brief Generate isotropic distribution.
 *
 * @param (rad) norm of the vector
 * @return (body) array for N-body particles
 * @param (idx) index of the N-body particle
 * @param (rand) state of the pseudo random number generator
 *
 * @sa RANDVAL
 */
static inline void isotropicDistribution(const real rad, iparticle body, const int idx, rand_state *rand)
{
  const real proj = RANDVAL(rand);
  body.pos[idx].z = rad * proj;
  real Rproj  = rad * SQRT(UNITY - proj * proj);

  real theta = TWO * CAST_D2R(M_PI) * UNIRAND(rand);
  body.pos[idx].x = Rproj * COS(theta);
  body.pos[idx].y = Rproj * SIN(theta);
}


/**
 * @fn gaussian
 *
 * @brief Generate Guassian with mean = 0.0, dispersion = 1.0 by Box--Muller method.
 *
 * @param (rand) state of the pseudo random number generator
 * @return resultant value
 *
 * @sa RANDVAL
 */
static inline real gaussian(rand_state *rand)
{
  real xx;
  real r2 = ZERO;
  while( (r2 < EPSILON) || (r2 + EPSILON > UNITY) ){
    xx = RANDVAL(rand);
    real yy = RANDVAL(rand);
    r2 = xx * xx + yy * yy;
  }

  return (SQRT(-TWO * LOG(r2) / r2) * xx);
}


/**
 * @fn makeUniformSphere
 *
 * @brief Generate uniform sphere.
 *
 * @param (num) total number of N-body particles
 * @return (body) array for N-body particles
 * @param (mtot) total mass of the sphere
 * @param (length) radius of the sphere
 * @param (sigma) velocity dispersion of the particles
 * @param (rand) state of the pseudo random number generator
 *
 * @sa isotropicDistribution
 * @sa gaussian
 */
void makeUniformSphere(ulong num, iparticle body, real mtot, real length, real sigma, rand_state *rand);
void makeUniformSphere(ulong num, iparticle body, real mtot, real length, real sigma, rand_state *rand)
{
  __NOTE__("%s\n", "start");

  real mass = mtot / (real)num;
  const real sigma1D = sigma / (real)M_SQRT3;
  for(ulong ii = 0; ii < num; ii++){
    isotropicDistribution(length * POW(UNIRAND(rand), ONE_THIRD), body, ii, rand);
    body.pos[ii].m = mass;

    body.acc[ii].x   = ZERO;
    body.acc[ii].y   = ZERO;
    body.acc[ii].z   = ZERO;
    body.acc[ii].pot = ZERO;

#ifdef  BLOCK_TIME_STEP
    body.vel[ii].x = sigma1D * gaussian(rand);
    body.vel[ii].y = sigma1D * gaussian(rand);
    body.vel[ii].z = sigma1D * gaussian(rand);
#else///BLOCK_TIME_STEP
    body.vx[ii] = sigma1D * gaussian();
    body.vy[ii] = sigma1D * gaussian();
    body.vz[ii] = sigma1D * gaussian();
#endif//BLOCK_TIME_STEP

    body.idx[ii] = ii;
  }/* for(ulong ii = 0; ii < num; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn shiftCenter
 *
 * @brief Shift center to set the center-of-mass at the coordinate origin with null bulk velocity.
 *
 * @param (num) total number of particles in the component
 * @return (body) physical quantities of N-body particles
 */
void shiftCenter(ulong num, iparticle body);
void shiftCenter(ulong num, iparticle body)
{
  __NOTE__("%s\n", "start");


  double com[3] = {0.0, 0.0, 0.0};
  double vel[3] = {0.0, 0.0, 0.0};
  double Mtot = 0.0;

  for(ulong ii = 0; ii < num; ii++){
    const double mass = CAST_R2D(body.pos[ii].m);
    Mtot   += mass;

    com[0] += mass * CAST_R2D(body.pos[ii].x);
    com[1] += mass * CAST_R2D(body.pos[ii].y);
    com[2] += mass * CAST_R2D(body.pos[ii].z);
#ifdef  BLOCK_TIME_STEP
    vel[0] += mass * CAST_R2D(body.vel[ii].x);
    vel[1] += mass * CAST_R2D(body.vel[ii].y);
    vel[2] += mass * CAST_R2D(body.vel[ii].z);
#else///BLOCK_TIME_STEP
    vel[0] += mass * CAST_R2D(body.vx[ii]);
    vel[1] += mass * CAST_R2D(body.vy[ii]);
    vel[2] += mass * CAST_R2D(body.vz[ii]);
#endif//BLOCK_TIME_STEP
  }/* for(ulong ii = 0; ii < num; ii++){ */

  double Minv = 1.0 / Mtot;
  com[0] *= Minv;  vel[0] *= Minv;
  com[1] *= Minv;  vel[1] *= Minv;
  com[2] *= Minv;  vel[2] *= Minv;
  const real rcom[3] = {CAST_D2R(com[0]), CAST_D2R(com[1]), CAST_D2R(com[2])};
  const real rvel[3] = {CAST_D2R(vel[0]), CAST_D2R(vel[1]), CAST_D2R(vel[2])};

  for(ulong ii = 0; ii < num; ii++){
    body.pos[ii].x -= rcom[0];
    body.pos[ii].y -= rcom[1];
    body.pos[ii].z -= rcom[2];

#ifdef  BLOCK_TIME_STEP
    body.time[ii].t0 = body.time[ii].t1 = 0.0;
    body.vel[ii].x -= rvel[0];
    body.vel[ii].y -= rvel[1];
    body.vel[ii].z -= rvel[2];
    body.vel[ii].dt = ZERO;
#else///BLOCK_TIME_STEP
    body.vx[ii] -= rvel[0];
    body.vy[ii] -= rvel[1];
    body.vz[ii] -= rvel[2];
#endif//BLOCK_TIME_STEP
  }/* for(ulong ii = 0; ii < num; ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn outputFundamentalInformationOfColdSphere
 *
 * @brief Print out fundamental information on the generated system.
 *
 * @param (Mtot) total mass of the sphere
 * @param (rad) radius of the sphere
 * @param (virial) Virial ratio of the sphere
 * @param (sigma) velocity dispersion of the particles
 * @param (Ntot) total number of N-body particles
 * @return (eps) value of Plummer softening length
 * @return (SnapshotInterval) time interval to write snapshot files
 * @return (ft) finish time of the simulation
 * @param (file) name of the simulation
 *
 * @sa getPresentDateInStrings
 */
void outputFundamentalInformationOfColdSphere
(const real Mtot, const real rad, const real virial, const real sigma, const ulong Ntot, const real eps, const double snapshotInterval, const double ft, char file[]);
void outputFundamentalInformationOfColdSphere
(const real Mtot, const real rad, const real virial, const real sigma, const ulong Ntot, const real eps, const double snapshotInterval, const double ft, char file[])
{
  __NOTE__("%s\n", "start");

  FILE *fp;
  char filename[256], date[64];


  /**< output fundamental information of the Cold sphere */
  sprintf(filename, "%s/%s.info.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "%s\n", "ERROR: fundamental information file of the Cold sphere couldn't open.");  }

  getPresentDateInStrings(date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Fundamental information about the Cold sphere\n");
  fprintf(fp, "\tgenerated on %s", date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Total mass of the Cold sphere Mtot is                    %e\n", Mtot);
  fprintf(fp, "Radius of the Cold sphere rad is                         %e\n", rad);
  fprintf(fp, "Virial ratio of the Cold sphere is                       %e\n", virial);
  fprintf(fp, "Gaussian velocity dispersion of the Cold sphere sigma is %e\n", sigma);
  fprintf(fp, "#############################################################################\n");

  const real Ms = Mtot;
  const real Ns = (real)Ntot;
  const real tff = (real)M_PI_2 * rad * SQRTRATIO(rad, TWO * newton * Ms);
  const real t2r = tff * Ns / ((real)32.0 * LOG(rad / eps));
  fprintf(fp, "Number of N-body particles to represent the Cold sphere is  %zu (= 2^%u)\n", Ntot, ilog2((uint)Ntot));
  fprintf(fp, "Length of Plummer softening is                              %e\n", eps);
  fprintf(fp, "Number of particles within the scale length is              %e\n", Ns);
  fprintf(fp, "Enclosed mass within the scale length is                    %e\n", Ms);
  fprintf(fp, "Free-fall time at the scale length in computational unit is %e\n", tff);
  fprintf(fp, "Two-body relaxation time at the scale length is             %e\n", t2r);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Snapshot interval in the computational unit               is %e\n", snapshotInterval);
  fprintf(fp, "Snapshot interval in the unit of free-fall time           is %e\n", snapshotInterval / tff);
  fprintf(fp, "Snapshot interval in the unit of two-body relaxation time is %e\n", snapshotInterval / t2r);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Final time of the simulation in the computational unit               is %e\n", ft);
  fprintf(fp, "Final time of the simulation in the unit of free-fall time           is %e\n", ft / tff);
  fprintf(fp, "Final time of the simulation in the unit of two-body relaxation time is %e\n", ft / t2r);
  fprintf(fp, "#############################################################################\n");
  fclose(fp);


  __NOTE__("%s\n", "end");
}


int main(int argc, char **argv)
{
  if( argc < 12 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 12);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -unit=<int>\n");
    __FPRINTF__(stderr, "          -Ntot=<unsigned long int>\n");
    __FPRINTF__(stderr, "          -Mtot=<real> -virial=<real> -rad=<real>\n");
    __FPRINTF__(stderr, "          -eps=<real> -ft=<real> -eta=<real>\n");
    __FPRINTF__(stderr, "          -ft=<real> -snapshotInterval=<real> -saveInterval=<real>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }

  char *file;  requiredCmdArg(getCmdArgStr( argc, (const char * const *)argv, "file", &file));
  ulong Ntot;  requiredCmdArg(getCmdArgUlg( argc, (const char * const *)argv, "Ntot", &Ntot));
  int   unit;  requiredCmdArg(getCmdArgInt( argc, (const char * const *)argv, "unit", &unit));
  setPhysicalConstantsAndUnitSystem(unit, 1);
  double tmp;
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "Mtot",   &tmp));  real Mtot   = CAST_D2R(tmp *   mass_astro2com);
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "virial", &tmp));  real virial = CAST_D2R(tmp);
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "rad",    &tmp));  real rad    = CAST_D2R(tmp * length_astro2com);
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "eps",    &tmp));  real eps    = CAST_D2R(tmp * length_astro2com);
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "ft",     &tmp));  double ft   =         (tmp *   time_astro2com);
  real eta;
  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "eta", &eta));
  double saveInterval;  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "saveInterval", &saveInterval));
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "snapshotInterval", &tmp));
  double snapshotInterval = ldexp(1.0, (int)floor(log2(tmp * time_astro2com)));


  /**< calculate velocity dispersion to realize the given Virial ratio */
  real sigma = SQRT((real)1.2 * newton * Mtot * virial / rad);

  /**< initialize pseudo random number generator */
  rand_state *rand;
  initRandNum(&rand);


  writeSettings(unit, Ntot, eps, eta, ft, snapshotInterval, saveInterval, file);
  iparticle body;
  ulong *idx;
  position *pos;
  acceleration *acc;
#ifdef  BLOCK_TIME_STEP
  velocity *vel;
  ibody_time *ti;
#else///BLOCK_TIME_STEP
  real *vx, *vy, *vz;
#endif//BLOCK_TIME_STEP
  allocParticleData((int)Ntot, &body, &idx, &pos, &acc,
#ifdef  BLOCK_TIME_STEP
		    &vel, &ti
#else///BLOCK_TIME_STEP
		    &vx, &vy, &vz
#endif//BLOCK_TIME_STEP
		    );

  makeUniformSphere(Ntot, body, Mtot, rad, sigma, rand);
  shiftCenter(Ntot, body);

  /**< release pseudo random number generator */
  freeRandNum(rand);

  outputFundamentalInformationOfColdSphere(Mtot, rad, virial, sigma, Ntot, eps, snapshotInterval, ft, file);


  double  time  = 0.0;
  double   dt   = 0.0;
  int   last  = 1;
  ulong steps = 0;

#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
#ifdef  MONITOR_ENERGY_ERROR
  static energyError relEneErr = {1.0, DBL_MIN};
#endif//MONITOR_ENERGY_ERROR
  static rebuildTree rebuild;
  static measuredTime measured;
  static autoTuningParam rebuildParam;
  static brentStatus status;
  static brentMemory memory;
  writeTentativeData(time, dt, steps, Ntot, body, file, &last, hdf5type
		     , rebuild, measured, rebuildParam, status, memory
#ifdef  MONITOR_ENERGY_ERROR
		     , relEneErr
#endif//MONITOR_ENERGY_ERROR
		     );
  removeHDF5DataType(hdf5type);
#else///USE_HDF5_FORMAT
  writeTentativeData(time, dt, steps, Ntot, body, file, &last);
#endif//USE_HDF5_FORMAT
  updateConfigFile(last, file);


  freeParticleData(idx, pos, acc,
#ifdef  BLOCK_TIME_STEP
		    vel, ti
#else///BLOCK_TIME_STEP
		    vx, vy, vz
#endif//BLOCK_TIME_STEP
		    );
  return (0);
}
