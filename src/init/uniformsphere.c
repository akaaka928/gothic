/*************************************************************************\
 *                                                                       *
                  last updated on 2017/01/25(Wed) 16:38:35
 *                                                                       *
 *    Making Initial Condition Code of N-body Simulation                 *
 *       Uniform sphere w/ Gaussian velocity dispersion                  *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
/* conversion from physical unit to computational unit must be performed internally */
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
#include "macro.h"
#include "myutil.h"
#include "constants.h"
#include "timer.h"
#include "name.h"
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#include "../misc/allocate.h"
//-------------------------------------------------------------------------
#include "../misc/tune.h"
#           if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#include "../misc/brent.h"
#        endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
//-------------------------------------------------------------------------
#include "../file/io.h"
//-------------------------------------------------------------------------
extern const double mass_astro2com, velocity_astro2com, length_astro2com, time_astro2com;
extern const real newton;
//-------------------------------------------------------------------------
#ifdef  USE_SFMT
#include "SFMT.h"
sfmt_t sfmt;
#define RANDPOS (CAST_D2R(sfmt_genrand_res53(&sfmt)))
#else///USE_SFMT
#include <gsl/gsl_rng.h>
gsl_rng *GSLRand;
#define RANDPOS ((real)gsl_rng_uniform(GSLRand))
#endif//USE_SFMT
#define RANDVAL (TWO * (RANDPOS) - UNITY)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void isotropicDistribution(const real rad, iparticle body, const int idx)
{
  //-----------------------------------------------------------------------
  const real proj = RANDVAL;
  body.pos[idx].z = rad * proj;
  real Rproj  = rad * SQRT(UNITY - proj * proj);
  //-----------------------------------------------------------------------
  real theta = TWO * (real)M_PI * RANDPOS;
  body.pos[idx].x = Rproj * COS(theta);
  body.pos[idx].y = Rproj * SIN(theta);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* Gaussian with mean = 0.0, dispersion = 1.0 by Box-Muller Method */
static inline real gaussian(void)
{
  //-----------------------------------------------------------------------
  real x, y, z, r2;
  //-----------------------------------------------------------------------
  r2 = ZERO;
  while( (r2 < EPSILON) || (r2 + EPSILON > (real)1.0) ){
    x = RANDVAL;
    y = RANDVAL;
    r2 = x * x + y * y;
  }
  //-----------------------------------------------------------------------
  z = SQRT(-(real)2.0 * LOG(r2) / r2) * x;
  //-----------------------------------------------------------------------
  return (z);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void makeUniformSphere(ulong num, iparticle body, real mtot, real length, real sigma);
void makeUniformSphere(ulong num, iparticle body, real mtot, real length, real sigma)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  real mass = mtot / (real)num;
  //-----------------------------------------------------------------------
  for(ulong ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    isotropicDistribution(length * POW(RANDPOS, ONE_THIRD), body, ii);
    body.pos[ii].m = mass;
    //---------------------------------------------------------------------
    body.acc[ii].x   = ZERO;
    body.acc[ii].y   = ZERO;
    body.acc[ii].z   = ZERO;
    body.acc[ii].pot = ZERO;
    //---------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
    body.vel[ii].x = sigma * gaussian();
    body.vel[ii].y = sigma * gaussian();
    body.vel[ii].z = sigma * gaussian();
#else///BLOCK_TIME_STEP
    body.vx[ii] = sigma * gaussian();
    body.vy[ii] = sigma * gaussian();
    body.vz[ii] = sigma * gaussian();
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
    body.idx[ii] = ii;
    //---------------------------------------------------------------------
  }/* for(ulong ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void shiftCenter(ulong num, iparticle body);
void shiftCenter(ulong num, iparticle body)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  double com[3] = {0.0, 0.0, 0.0};
  double vel[3] = {0.0, 0.0, 0.0};
  double Mtot = 0.0;
  //-----------------------------------------------------------------------
  for(ulong ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    const double mass = (double)body.pos[ii].m;
    Mtot   += mass;
    //---------------------------------------------------------------------
    com[0] += mass * (double)body.pos[ii].x;
    com[1] += mass * (double)body.pos[ii].y;
    com[2] += mass * (double)body.pos[ii].z;
#ifdef  BLOCK_TIME_STEP
    vel[0] += mass * (double)body.vel[ii].x;
    vel[1] += mass * (double)body.vel[ii].y;
    vel[2] += mass * (double)body.vel[ii].z;
#else///BLOCK_TIME_STEP
    vel[0] += mass * (double)body.vx[ii];
    vel[1] += mass * (double)body.vy[ii];
    vel[2] += mass * (double)body.vz[ii];
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
  }/* for(ulong ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  double Minv = 1.0 / Mtot;
  com[0] *= Minv;  vel[0] *= Minv;
  com[1] *= Minv;  vel[1] *= Minv;
  com[2] *= Minv;  vel[2] *= Minv;
  const real rcom[3] = {(real)com[0], (real)com[1], (real)com[2]};
  const real rvel[3] = {(real)vel[0], (real)vel[1], (real)vel[2]};
  //-----------------------------------------------------------------------
  for(ulong ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    body.pos[ii].x -= rcom[0];
    body.pos[ii].y -= rcom[1];
    body.pos[ii].z -= rcom[2];
    //---------------------------------------------------------------------
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
    //---------------------------------------------------------------------
  }/* for(ulong ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
#if 0
  printf("position shift = (%e, %e, %e)\n", com[0], com[1], com[2]);
  printf("velocity shift = (%e, %e, %e)\n", vel[0], vel[1], vel[2]);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void outputFundamentalInformationOfColdSphere
(const real Mtot, const real rad, const real sigma, const ulong Ntot, const real eps, const double snapshotInterval, const double ft, char file[]);
void outputFundamentalInformationOfColdSphere
(const real Mtot, const real rad, const real sigma, const ulong Ntot, const real eps, const double snapshotInterval, const double ft, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  FILE *fp;
  char filename[256], date[64];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output fundamental information of the Cold sphere */
  sprintf(filename, "%s/%s.info.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "%s\n", "ERROR: fundamental information file of the Cold sphere couldn't open.");
  }
  //-----------------------------------------------------------------------
  getPresentDateInStrings(date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Fundamental information about the Cold sphere\n");
  fprintf(fp, "\tgenerated on %s", date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Total mass of the Cold sphere Mtot is                    %e\n", Mtot);
  fprintf(fp, "Radius of the Cold sphere rad is                         %e\n", rad);
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
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //-----------------------------------------------------------------------
  if( argc < 12 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 12);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -unit=<int>\n");
    __FPRINTF__(stderr, "          -Ntot=<unsigned long int>\n");
    __FPRINTF__(stderr, "          -Mtot=<real> -sigma=<real> -rad=<real>\n");
    __FPRINTF__(stderr, "          -eps=<real> -ft=<real> -eta=<real>\n");
    __FPRINTF__(stderr, "          -ft=<real> -snapshotInterval=<real> -saveInterval=<real>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }
  //-----------------------------------------------------------------------
  /* setPhysicalConstantsAndUnitSystem(UNITSYSTEM, 1); */
  //-----------------------------------------------------------------------
  char *file;  requiredCmdArg(getCmdArgStr( argc, (const char * const *)argv, "file", &file));
  ulong Ntot;  requiredCmdArg(getCmdArgUlg( argc, (const char * const *)argv, "Ntot", &Ntot));
  int   unit;  requiredCmdArg(getCmdArgInt( argc, (const char * const *)argv, "unit", &unit));
  setPhysicalConstantsAndUnitSystem(unit, 1);
  double tmp;
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "Mtot",  &tmp));  real Mtot  = (real)(tmp *     mass_astro2com);
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "sigma", &tmp));  real sigma = (real)(tmp * velocity_astro2com);
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "rad",   &tmp));  real rad   = (real)(tmp *   length_astro2com);
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "eps",   &tmp));  real eps   = (real)(tmp *   length_astro2com);
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "ft",    &tmp));  double ft  =       (tmp *     time_astro2com);
  real eta;
  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "eta", &eta));
  double saveInterval;  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "saveInterval", &saveInterval));
  requiredCmdArg(getCmdArgDbl( argc, (const char * const *)argv, "snapshotInterval", &tmp));
  double snapshotInterval = ldexp(1.0, (int)floor(log2(tmp * time_astro2com)));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* settings about random numbers */
#ifdef  USE_SFMT
  sfmt_init_gen_rand(&sfmt, 5489);
  /* sfmt_init_gen_rand(&sfmt, 19650218); */
  /* 5489 for 32 bit Mersenne Twister */
  /* 19650218 for 64 bit Mersenne Twister */
#else///USE_SFMT
  const gsl_rng_type *RandType;
  gsl_rng_env_setup();
  RandType = gsl_rng_mt19937;
  GSLRand = gsl_rng_alloc(RandType);
  gsl_rng_set(GSLRand, 5489);
#endif//USE_SFMT
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  writeSettings(unit, Ntot, eps, eta, ft, snapshotInterval, saveInterval, file);
  //-----------------------------------------------------------------------
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
  makeUniformSphere(Ntot, body, Mtot, rad, sigma);
  shiftCenter(Ntot, body);
  //-----------------------------------------------------------------------
  outputFundamentalInformationOfColdSphere(Mtot, rad, sigma, Ntot, eps, snapshotInterval, ft, file);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  double  time  = 0.0;
  double   dt   = 0.0;
  int   last  = 1;
  ulong steps = 0;
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
#ifdef  MONITOR_ENERGY_ERROR
  static energyError relEneErr = {1.0, DBL_MIN};
#endif//MONITOR_ENERGY_ERROR
  static rebuildTree rebuild;
  static measuredTime measured;
#ifdef  WALK_TREE_COMBINED_MODEL
  static autoTuningParam rebuildParam;
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
  static brentStatus status;
  static brentMemory memory;
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
  writeTentativeData(time, dt, steps, Ntot, body, file, &last, hdf5type
		     , rebuild, measured
#ifdef  WALK_TREE_COMBINED_MODEL
		     , rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
		     , status, memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
		     , relEneErr
#endif//MONITOR_ENERGY_ERROR
		     );
  removeHDF5DataType(hdf5type);
#else///USE_HDF5_FORMAT
  writeTentativeData(time, dt, steps, Ntot, body, file, &last);
#endif//USE_HDF5_FORMAT
  updateConfigFile(last, file);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  freeParticleData(idx, pos, acc,
#ifdef  BLOCK_TIME_STEP
		    vel, ti
#else///BLOCK_TIME_STEP
		    vx, vy, vz
#endif//BLOCK_TIME_STEP
		    );
  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
