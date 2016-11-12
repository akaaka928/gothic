/*************************************************************************\
 *                                                                       *
                  last updated on 2016/11/11(Fri) 10:56:26
 *                                                                       *
 *    MAGI: "MAny-component Galactic Initial-conditions" generator       *
 *    Making Initial Condition Code of N-body Simulation                 *
 *       Based on distribution function given by Eddington's formula     *
 *       arbitrary spherical symmetric and isotropic DF                  *
 *       extension in multi-components is based on Kazantzidis+06        *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
/* conversion from physical unit to computational unit must be performed internally */
//-------------------------------------------------------------------------
#define USE_SZIP_COMPRESSION
//-------------------------------------------------------------------------
/* #define RESET_ROTATION_AXIS */
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
#       include <hdf5.h>
#       include <hdf5lib.h>
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <name.h>
#include <myutil.h>
#include <constants.h>
#include <timer.h>
#include <rotate.h>
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
#include "magi.h"
#include "king.h"
#include "profile.h"
#include "eddington.h"
#include "table.h"
#include "abel.h"
/* #include "disk.h" */
#include "potdens.h"
#include "diskDF.h"
//-------------------------------------------------------------------------
extern const real newton;
extern const double     mass_astro2com,     mass2astro;extern const char        mass_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double   length_astro2com,   length2astro;extern const char      length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double     time_astro2com,     time2astro;extern const char        time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                     velocity2astro;extern const char    velocity_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                  col_density2astro;extern const char col_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                      density2astro;extern const char     density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                      senergy2astro;extern const char     senergy_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
//-------------------------------------------------------------------------
gsl_rng *GSLRand;
#define UNIRAND_DBL ((double)gsl_rng_uniform(GSLRand))
#define UNIRAND     (  (real)gsl_rng_uniform(GSLRand))
#define RANDVAL     (TWO * (UNIRAND) - UNITY)
//-------------------------------------------------------------------------
double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void isotropicDistribution(const real rad, real *vecx, real *vecy, real *vecz)
{
  //-----------------------------------------------------------------------
  const real proj = RANDVAL;
  *vecz = rad * proj;
  real Rproj = rad * SQRT(UNITY - proj * proj);
  //-----------------------------------------------------------------------
  real theta = TWO * (real)M_PI * UNIRAND;
  *vecx = Rproj * COS(theta);
  *vecy = Rproj * SIN(theta);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifdef __ICC
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#     pragma warning (disable:869)
#endif//__ICC
bool isInnerParticle4spherical(const double x2, const double y2, const double z2, const double rmax2, const double zmax2);
bool isInnerParticle4disk     (const double x2, const double y2, const double z2, const double Rmax2, const double zmax2);
bool isInnerParticle4spherical(const double x2, const double y2, const double z2, const double rmax2, const double zmax2){  return ((x2 + y2 + z2) < rmax2);}
bool isInnerParticle4disk     (const double x2, const double y2, const double z2, const double Rmax2, const double zmax2){  return (((x2 + y2) < Rmax2) && (z2 < zmax2));}
#ifdef __ICC
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#     pragma warning ( enable:869)
#endif//__ICC
//-------------------------------------------------------------------------
void shiftCenter(const ulong num, const ulong head, iparticle body, const profile_cfg cfg);
void shiftCenter(const ulong num, const ulong head, iparticle body, const profile_cfg cfg)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate center-of-mass and bulk-motion */
  //-----------------------------------------------------------------------
  double com[3] = {0.0, 0.0, 0.0};
  double vel[3] = {0.0, 0.0, 0.0};
  double Mtot = 0.0;
  //-----------------------------------------------------------------------
  /* use particles (r < 3 rs            ) for spherical components */
  /* use particles (R < 3 Rd, |z| < 3 zd) for      disk components */
  const double rmax2 = 9.0 * cfg.rs * cfg.rs;
  const double zmax2 = 9.0 * ((cfg.kind < 0) ? (cfg.zd * cfg.zd) : (cfg.rs * cfg.rs));
  bool (*isInnerParticle)(double, double, double, double, double) = (cfg.kind < 0) ? isInnerParticle4disk : isInnerParticle4spherical;
  //-----------------------------------------------------------------------
  for(ulong ii = head; ii < head + num; ii++){
    //---------------------------------------------------------------------
    const double xx = (double)body.pos[ii].x;
    const double yy = (double)body.pos[ii].y;
    const double zz = (double)body.pos[ii].z;
    //---------------------------------------------------------------------
    if( isInnerParticle(xx * xx, yy * yy, zz * zz, rmax2, zmax2) ){
      //-------------------------------------------------------------------
      const double mass = (double)body.pos[ii].m;
      Mtot   += mass;
      //-------------------------------------------------------------------
      com[0] += mass * xx;
      com[1] += mass * yy;
      com[2] += mass * zz;
#ifdef  BLOCK_TIME_STEP
      vel[0] += mass * (double)body.vel[ii].x;
      vel[1] += mass * (double)body.vel[ii].y;
      vel[2] += mass * (double)body.vel[ii].z;
#else///BLOCK_TIME_STEP
      vel[0] += mass * (double)body.vx[ii];
      vel[1] += mass * (double)body.vy[ii];
      vel[2] += mass * (double)body.vz[ii];
#endif//BLOCK_TIME_STEP
      //-------------------------------------------------------------------
    }/* if( isInnerParticle(xx * xx, yy * yy, zz * zz, rmax2, zmax2) ){ */
    //---------------------------------------------------------------------
  }/* for(ulong ii = head; ii < head + num; ii++){ */
  //-----------------------------------------------------------------------
  double Minv = 1.0 / (DBL_MIN + Mtot);
  com[0] *= Minv;  vel[0] *= Minv;
  com[1] *= Minv;  vel[1] *= Minv;
  com[2] *= Minv;  vel[2] *= Minv;
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# center-of-mass shift: %e, %e, %e\n", com[0], com[1], com[2]);
  fprintf(stdout, "#    bulk motion shift: %e, %e, %e\n", vel[0], vel[1], vel[2]);
  fflush(stdout);
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* shift the coordinate system to the center-of-mass rest frame */
  //-----------------------------------------------------------------------
#if 1
  const real rcom[3] = {(real)com[0], (real)com[1], (real)com[2]};
  const real rvel[3] = {(real)vel[0], (real)vel[1], (real)vel[2]};
#else
  const real rcom[3] = {ZERO, ZERO, ZERO};
  const real rvel[3] = {ZERO, ZERO, ZERO};
#endif
#pragma omp parallel for
  for(ulong ii = head; ii < head + num; ii++){
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
  }/* for(ulong ii = head; ii < head + num; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  RESET_ROTATION_AXIS
  //-----------------------------------------------------------------------
  /* calculate angular momentum vector */
  //-----------------------------------------------------------------------
  double amom[3] = {0.0, 0.0, 0.0};
  for(ulong ii = head; ii < head + num; ii++){
    //---------------------------------------------------------------------
    const double rx = (double)body.pos[ii].x;
    const double ry = (double)body.pos[ii].y;
    const double rz = (double)body.pos[ii].z;
    //---------------------------------------------------------------------
    if( isInnerParticle(rx * rx, ry * ry, rz * rz, rmax2, zmax2) ){
      //-------------------------------------------------------------------
      const double mass = (double)body.pos[ii].m;
      //-------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
      const double px = (double)body.vel[ii].x * mass;
      const double py = (double)body.vel[ii].y * mass;
      const double pz = (double)body.vel[ii].z * mass;
#else///BLOCK_TIME_STEP
      const double px = (double)body.vx[ii] * mass;
      const double py = (double)body.vy[ii] * mass;
      const double pz = (double)body.vz[ii] * mass;
#endif//BLOCK_TIME_STEP
      //-------------------------------------------------------------------
      amom[0] += ry * pz - rz * py;
      amom[1] += rz * px - rx * pz;
      amom[2] += rx * py - ry * px;
      //-------------------------------------------------------------------
    }/* if( isInnerParticle(rx * rx, ry * ry, rz * rz, rmax2, zmax2) ){ */
    //---------------------------------------------------------------------
  }/* for(ulong ii = head; ii < head + num; ii++){ */
  //-----------------------------------------------------------------------
  /* rotate galaxy (if necessary) */
  //-----------------------------------------------------------------------
#if 1
  const double L2 = amom[0] * amom[0] + amom[1] * amom[1] + amom[2] * amom[2];
#else
  const double L2 = 1.0;
#endif
  if( L2 > 1.0e-6 ){
    //---------------------------------------------------------------------
#if 1
    real ini[3] = {(real)amom[0], (real)amom[1], (real)amom[2]};
    real fin[3] = {ZERO, ZERO, UNITY};
#else
    real ini[3] = {ZERO, ZERO, UNITY};
    real fin[3] = {ZERO, UNITY, ZERO};
#endif
    //---------------------------------------------------------------------
    real rot[3][3], inv[3][3];
    initRotationMatrices(ini, fin, rot, inv);
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(ulong ii = head; ii < head + num; ii++){
      //-------------------------------------------------------------------
      real bfr[3], aft[3];
      //-------------------------------------------------------------------
      /* rotate position */
      bfr[0] = body.pos[ii].x;
      bfr[1] = body.pos[ii].y;
      bfr[2] = body.pos[ii].z;
      rotateVector(bfr, rot, aft);
      body.pos[ii].x = aft[0];
      body.pos[ii].y = aft[1];
      body.pos[ii].z = aft[2];
      //-------------------------------------------------------------------
      /* rotate velocity */
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
      //-------------------------------------------------------------------
    }/* for(ulong ii = head; ii < head + num; ii++){ */
    //---------------------------------------------------------------------
  }/* if( L2 > 1.0e-6 ){ */
  //-----------------------------------------------------------------------
#endif//RESET_ROTATION_AXIS
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline double getDF(const double ene, dist_func *df, const double Emin, const double invEbin)
{
  //-----------------------------------------------------------------------
#if 1
  const int ll = (int)((ene - Emin) * invEbin);
  return (df[ll].val + (df[ll + 1].val - df[ll].val) * (ene - df[ll].ene) / (df[ll + 1].ene - df[ll].ene));
#else
  int ll = (int)floor((ene - Emin) * invEbin);
  if( ll <           0 )    ll =           0;
  if( ll > NENEBIN - 2 )    ll = NENEBIN - 2;
  const int rr = ll + 1;
  const double alpha = (ene - df[ll].ene) * invEbin;
  return ((1.0 - alpha) * df[ll].val + alpha * df[rr].val);
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
double distributeSpheroidParticles(ulong *Nuse, iparticle body, const real mass, profile_cfg cfg, profile *prf, dist_func *df);
double distributeSpheroidParticles(ulong *Nuse, iparticle body, const real mass, profile_cfg cfg, profile *prf, dist_func *df)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const ulong num = cfg.num;
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "#\n#\n# start distributing spherical component particles (%zu bodies: [%zu:%zu])\n", num, *Nuse, (*Nuse) + num - 1);
  fflush(stdout);
  const ulong nunit = (ulong)ceilf(0.1f * (float)num);
  ulong stage = 1;
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------
  const double    Emin = df[          0].ene;
  const double    Emax = df[NENEBIN - 1].ene;
  const double invEbin = (double)(NENEBIN - 1) / (Emax - Emin);
  //-----------------------------------------------------------------------
  double fmax = 0.0;
  int   iout = 0;
  for(int ii = 0; ii < NRADBIN; ii++){
    double floc = prf[ii].rad * prf[ii].rad * prf[ii].rho;
    if(                            floc > fmax ){      fmax = floc;    }
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */
  const double Ecut = (iout != 0) ? prf[iout].psi_tot : Emin;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const double Mmax = cfg.Mtot;
  const double Mmin = prf[0].enc;
  //-----------------------------------------------------------------------
  for(ulong ii = *Nuse; ii < *Nuse + num; ii++){
    //---------------------------------------------------------------------
    /* set spatial distribution by table search */
    //---------------------------------------------------------------------
    const double tmp = Mmin + (Mmax - Mmin) * UNIRAND_DBL;
    int ll = 0;
    int rr = NRADBIN - 1;
    while( true ){
      //-------------------------------------------------------------------
      uint cc = (ll + rr) >> 1;
      //-------------------------------------------------------------------
      if( (prf[cc].enc - tmp) * (prf[ll].enc - tmp) <= 0.0 )	rr = cc;
      else	                                                ll = cc;
      //-------------------------------------------------------------------
      if( (ll + 1) == rr )	break;
      //-------------------------------------------------------------------
    }/* while( true ){ */
    //---------------------------------------------------------------------
    const double alpha = (tmp - prf[ll].enc) / (prf[rr].enc - prf[ll].enc);
    const double rad = (1.0 - alpha) * prf[ll].rad + alpha * prf[rr].rad;
    isotropicDistribution((real)rad, &(body.pos[ii].x), &(body.pos[ii].y), &(body.pos[ii].z));
    __NOTE__("position of %zu-th particle determined: rad = %e\n", ii, rad);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* determine velocity distribution by rejection method */
    //---------------------------------------------------------------------
#if 1
    const double psi = (1.0 - alpha) * prf[ll].psi_tot + alpha * prf[rr].psi_tot;
#else
    const double psi = prf[ll].psi_tot + (prf[rr].psi_tot - prf[ll].psi_tot) * (tmp - prf[ll].enc) / (prf[rr].enc - prf[ll].enc);
#endif
    /* const double vesc = sqrt(2.0 * psi); */
    const double vesc = sqrt(2.0 * (psi - Ecut));
    const double v2Fmax = vesc * vesc * getDF(psi, df, Emin, invEbin);
    double vel;
    while( true ){
      //-------------------------------------------------------------------
      vel = vesc * UNIRAND_DBL;
      //-------------------------------------------------------------------
      const double ene = psi - 0.5 * vel * vel;
      const double val = vel * vel * getDF(ene, df, Emin, invEbin);
      const double try = v2Fmax * UNIRAND_DBL;
      if( val > try )	break;
      //-------------------------------------------------------------------
    }/* while( true ){ */
    //---------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
    isotropicDistribution((real)vel, &(body.vel[ii].x), &(body.vel[ii].y), &(body.vel[ii].z));
#else///BLOCK_TIME_STEP
    isotropicDistribution((real)vel, &(body.vx[ii]), &(body.vy[ii]), &(body.vz[ii]));
#endif//BLOCK_TIME_STEP
    __NOTE__("velocity of %zu-th particle determined\n", ii);
    //---------------------------------------------------------------------
    body.acc[ii].x   = body.acc[ii].y = body.acc[ii].z = ZERO;
    body.pos[ii].m   = mass;
    body.acc[ii].pot = ZERO;
    //---------------------------------------------------------------------
    body.idx[ii] = ii;
    //---------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
    if( (ii - (*Nuse)) == (stage * nunit) ){
      fprintf(stdout, "# ~%zu%% completed\n", stage * 10);
      fflush(stdout);
      stage++;
    }/* if( (ii - (*Nuse)) == (stage * nunit) ){ */
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------
  }/* for(ulong ii = *Nuse; ii < *Nuse + num; ii++){ */
  //-----------------------------------------------------------------------
  *Nuse += num;
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# finish distributing spherical component particles (%zu bodies)\n", num);
  fflush(stdout);
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (Ecut);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void outputFundamentalInformation
(const int unit, const int kind, const int skind, profile_cfg *cfg, profile **prf, dist_func **df, const ulong Ntot, const real eps, const double snapshotInterval, const double ft, char file[]);
void writeDiskData(char *file, const int ndisk, const int maxLev, disk_data *disk);
//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //-----------------------------------------------------------------------
  /* read input arguments */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* read input arguments do not depend on the unit system adopted in the numerical simulation */
  char *file;  requiredCmdArg(getCmdArgStr( argc, (const char * const *)argv,   "file", &file));
  char *fcfg;  requiredCmdArg(getCmdArgStr( argc, (const char * const *)argv, "config", &fcfg));
  ulong Ntot;  requiredCmdArg(getCmdArgUlg( argc, (const char * const *)argv,   "Ntot", &Ntot));
  real   eta;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,    "eta", &eta));
  double saveInterval;  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "saveInterval", &saveInterval));
  //-----------------------------------------------------------------------
  /* set unit system by reading the configuration file about physical parameters of the initial distribution */
  int unit, kind;
  profile_cfg *cfg;
  readProfileCfg(fcfg, &unit, &kind, &cfg);
  if( kind > NKIND_MAX ){    __KILL__(stderr, "ERROR: kind(= %d) must be smaller than %d\n", kind, NKIND_MAX);  }
  //-----------------------------------------------------------------------
  /* read input arguments depend on the unit system adopted in the numerical simulation */
  double tmp;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "snapshotInterval", &tmp));
  double snapshotInterval = ldexp(1.0, (int)floor(log2(tmp * time_astro2com)));
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "eps", &tmp));  real   eps = (real)(tmp * length_astro2com);
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "ft",  &tmp));  double  ft =       (tmp *   time_astro2com);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set number of particles */
  //-----------------------------------------------------------------------
  ulong Nrem = 0;
  for(int ii = 0; ii < kind; ii++)
    Nrem += (cfg[ii].forceNum == 1) ? cfg[ii].num : 0;
  if( Nrem > Ntot ){
    __KILL__(stderr, "ERROR: the sum of number of particles for each component (%zu) exceeds the specified total number of particles (%zu).\n", Nrem, Ntot);
  }/* if( Nrem > Ntot ){ */
  Nrem = Ntot - Nrem;
  //-----------------------------------------------------------------------
  /* set number of particles to represent each profile */
  double Mtot = 0.0;
  for(int ii = 0; ii < kind; ii++)
    if( cfg[ii].forceNum != 1 )
      Mtot += cfg[ii].Mtot;
  const double Minv = 1.0 / Mtot;
  ulong Nuse = 0;
  int kidx = kind;
  double Mmax = 0.0;
  for(int ii = 0; ii < kind; ii++){
    //---------------------------------------------------------------------
    if( cfg[ii].forceNum != 1 ){
      //-------------------------------------------------------------------
      /* number of particles is determined by mass fraction */
      cfg[ii].num = (ulong)(cfg[ii].Mtot * Minv * (double)Nrem);
      Nuse += cfg[ii].num;
      //-------------------------------------------------------------------
      if( cfg[ii].Mtot > Mmax ){
	//-----------------------------------------------------------------
	kidx = ii;
	Mmax = cfg[ii].Mtot;
	//-----------------------------------------------------------------
      }/* if( cfg[ii].Mtot > Mmax ){ */
      //-------------------------------------------------------------------
    }/* if( cfg[ii].forceNum != 1 ){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < kind; ii++){ */
  if( (kidx == kind) && (Nuse != Nrem) ){
    __KILL__(stderr, "ERROR: mismatch about number of particles detected (Nuse = %zu, Nrem = %zu) with %d components\n", Nuse, Nrem, kind);
  }/* if( (kidx == kind) && (Nuse != Nrem) ){ */
  if( Nuse != Nrem )
    cfg[kidx].num += (Nrem - Nuse);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* GSL initialization */
  //-----------------------------------------------------------------------
  /* initialize random number provided by GSL */
  const gsl_rng_type *RandType;
  gsl_rng_env_setup();
  RandType = gsl_rng_mt19937;
  GSLRand  = gsl_rng_alloc(RandType);
  gsl_rng_set(GSLRand, 5489);
  //-----------------------------------------------------------------------
  /* initialize the table for Gaussian Quadrature provided by GSL */
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
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set global settings */
  //-----------------------------------------------------------------------
  /* identifier for disk components is negative value (-1, -2) */
  /* identifier for spherical components is positive value (0, 1, 2, ...) */
  int skind = kind;
  for(int ii = 0; ii < kind; ii++)
    if( cfg[ii].kind < 0 )
      skind--;
  const int ndisk = kind - skind;
  const bool addDisk = (ndisk != 0) ? true : false;
  if( addDisk )
    for(int ii = skind; ii < kind; ii++)
      if( cfg[ii].kind >= 0 ){      	__KILL__(stderr, "ERROR: disk component must be last component(s).\n\tModify \"%s/%s\".\n", CFGFOLDER, fcfg);      }
  //-----------------------------------------------------------------------
  bool cutoff = false;
  double rmax = 0.0;
  for(int ii = 0; ii < kind; ii++){
    cutoff |= cfg[ii].cutoff;
    if( cfg[ii].cutoff ){      if( rmax < cfg[ii].rc )	rmax = cfg[ii].rc;    }
    else{                      if( rmax < cfg[ii].rs )	rmax = cfg[ii].rs;    }
  }/* for(int ii = 0; ii < kind; ii++){ */
  //-----------------------------------------------------------------------
  if( cutoff )    rmax *= (double)(NINTBIN * 10);/* <-- sufficiently greater than outermost radius times NINTBIN */
  else            rmax *= 1000.0;
  const double logrmin = log10(MINRAD);
  const double logrmax = log10(rmax);
  const double logrbin = (logrmax - logrmin) / (double)(4 + NRADBIN);
  const double invlogrbin = 1.0 / logrbin;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set distribution function */
  //-----------------------------------------------------------------------
  /* memory allocation for spherical component(s) */
  profile **prf, *_prf;
  /* 2 * 2 bins are added in the both edge */
  _prf = (profile  *)malloc(sizeof(profile  ) * kind * (4 + NRADBIN));  if( _prf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _prf\n");  }
  prf  = (profile **)malloc(sizeof(profile *) * kind                );  if(  prf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  prf\n");  }
  for(int ii = 0; ii < kind; ii++)
    prf[ii] = _prf + ii * (4 + NRADBIN);
#pragma omp parallel for
  for(int ii = 0; ii < kind; ii++)
    for(int jj = 0; jj < 4 + NRADBIN; jj++)
      prf[ii][jj].rad = pow(10.0, logrmin + logrbin * (double)jj);
  /* memory allocation for disk component */
  int maxLev;
  disk_data  *disk_info;
  double *disk_hor, *disk_ver, *disk_node_hor, *disk_node_ver, *disk_pot, *disk_dPhidR, *disk_d2PhidR2;
  double *sph_rad, *sph_enc, *sph_rho;
  double *disk_rho, *disk_rhoSum, *disk_rhoTot, *disk_Sigma, *disk_sigmaz, *disk_enc;
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  double *disk_zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  double *spline_xx, *spline_ff, *spline_f2, *spline_bp;
  if( addDisk ){
    //---------------------------------------------------------------------
    /* allocate required arrays */
    //---------------------------------------------------------------------
    allocDiskProfile(ndisk, &disk_info, &cfg[skind], &maxLev,
		     &disk_hor, &disk_ver, &disk_node_hor, &disk_node_ver,
		     &disk_pot, &disk_rho, &disk_rhoSum, &disk_rhoTot,
		     &disk_dPhidR, &disk_d2PhidR2,
		     &disk_Sigma, &disk_sigmaz, &disk_enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
		     &disk_zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
		     &sph_rad, &sph_rho, &sph_enc,
		     &spline_xx, &spline_ff, &spline_f2, &spline_bp);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* commit fundamental information */
    //---------------------------------------------------------------------
    for(int ii = 0; ii < ndisk; ii++){
      //-------------------------------------------------------------------
      /* disk_info[ii].cfg        = &cfg[skind + ii]; */
      disk_info[ii].prf        =  prf[skind + ii];
      disk_info[ii].   logrbin =    logrbin;
      disk_info[ii].invlogrbin = invlogrbin;
      //-------------------------------------------------------------------
      /* initialize for Toomre's Q-value analysis */
      disk_info[ii].cfg->vcirc_max   = -1.0;
      disk_info[ii].cfg->vcirc_max_R = -1.0;
      disk_info[ii].cfg->Qmin = DBL_MAX;
      disk_info[ii].cfg->passed = false;
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < ndisk; ii++){ */
    //---------------------------------------------------------------------
  }/* if( addDisk ){ */
  //-----------------------------------------------------------------------
  /* set density profile and mass profile for spherical component(s) */
  for(int ii = 0; ii < skind; ii++){
    //---------------------------------------------------------------------
    profile_abel_cfg dummy;    dummy.invRd = 1.0;
    //---------------------------------------------------------------------
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
    case CENTRALBH:      setContributionByCentralBH      (prf[ii],  cfg[ii]                                                                          );      break;
    default:      __KILL__(stderr, "ERROR: undefined profile ``%d'' specified\n", cfg[ii].kind);      break;
    }/* switch( cfg[ii].kind ){ */
    //---------------------------------------------------------------------
    if( cfg[ii].kind != CENTRALBH ){
#if 1
      integrateDensityProfile(prf[ii], logrbin, cfg[ii].Mtot, cfg[ii].cutoff, cfg[ii].rc, cfg[ii].rc_width);
#else
      integrateDensityProfile(prf[ii], logrbin, cfg[ii].Mtot, cfg[ii].cutoff, cfg[ii].rc, (cfg[ii].rs < (cfg[ii].rc * 0.1)) ? (cfg[ii].rc * 0.1) : (cfg[ii].rs));
#endif
    }/* if( cfg[ii].kind != CENTRALBH ){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < skind; ii++){ */
  //-----------------------------------------------------------------------
  /* evaluate sum of density, enclosed mass and potential of all spherical component(s) */
#pragma omp parallel for
  for(int ii = 0; ii < 4 + NRADBIN; ii++){
    //---------------------------------------------------------------------
    double rho = 0.0;
    double enc = 0.0;
    double psi = 0.0;
    for(int kk = 0; kk < skind; kk++){
      rho += prf[kk][ii].rho;
      enc += prf[kk][ii].enc;
      psi += prf[kk][ii].psi;
    }/* for(int kk = 0; kk < skind; kk++){ */
    //---------------------------------------------------------------------
    for(int kk = 0; kk < kind; kk++){
      prf[kk][ii].rho_tot = rho;
      prf[kk][ii].enc_tot = enc;
      prf[kk][ii].psi_tot = psi;
    }/* for(int kk = 0; kk < kind; kk++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
  //-----------------------------------------------------------------------
  /* set density profile for the disk component */
  if( addDisk ){
    //---------------------------------------------------------------------
    /* set disk_radius, disk_height, disk_pot */
    makeDiskPotentialTable(ndisk, maxLev, disk_info);
#if 0
    writeDiskData(file, ndisk, maxLev, disk_info);
    exit(0);
#endif
    //---------------------------------------------------------------------
    /* set profile of spherical averaged density, mass and potential */
    integrateSphericalDensityProfile(ndisk, maxLev, disk_info);
#if 0
    writeDiskData(file, ndisk, maxLev, disk_info);
    exit(0);
#endif
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < 4 + NRADBIN; ii++){
      //-------------------------------------------------------------------
      double rho = 0.0;
      double enc = 0.0;
      double psi = 0.0;
      for(int kk = skind; kk < kind; kk++){
	rho += prf[kk][ii].rho;
	enc += prf[kk][ii].enc;
	psi += prf[kk][ii].psi;
      }/* for(int kk = skind; kk < kind; kk++){ */
      //-------------------------------------------------------------------
      for(int kk = 0; kk < skind; kk++){
	prf[kk][ii].rho_tot += rho;
	prf[kk][ii].enc_tot += enc;
	prf[kk][ii].psi_tot += psi;
      }/* for(int kk = 0; kk < skind; kk++){ */
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
    //---------------------------------------------------------------------
  }/* if( addDisk ){ */
  //-----------------------------------------------------------------------
  /* integrate Eddington's formula numerically */
  dist_func **fene, *_fene;
  _fene = (dist_func  *)malloc(sizeof(dist_func  ) * skind * NENEBIN);  if( _fene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _fene\n");  }
  fene  = (dist_func **)malloc(sizeof(dist_func *) * skind          );  if(  fene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  fene\n");  }
  for(int ii = 0; ii < skind; ii++)
    fene[ii] = _fene + ii * NENEBIN;
  integrateEddingtonFormula(skind, prf, fene);
  //-----------------------------------------------------------------------
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  calcColumnDensityProfile(skind, prf, logrmax, cfg);
#endif//MAKE_COLUMN_DENSITY_PROFILE
  //-----------------------------------------------------------------------
  if( addDisk ){
    /* differentiate potential along the radial direction on the equatorial plane */
    diffAxisymmetricPotential(maxLev, disk_info[0]);
#if 0
    writeDiskData(file, ndisk, maxLev, disk_info);
    exit(0);
#endif
    /* set velocity dispersion in vertical direction */
    calcVerticalVdisp(ndisk, maxLev, disk_info);
    for(int ii = skind; ii < kind; ii++){
      /* cfg[ii].vdispz0 = disk_vdisp[0]; */
      cfg[ii].vdispz0 = disk_info[ii - skind].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, 0)];
      if( cfg[ii].vdispR0 < 0.0 )
	cfg[ii].vdispR0 = cfg[ii].vdispz0;
    }/* for(int ii = skind; ii < kind; ii++){ */
    //---------------------------------------------------------------------
    /* output fundamental quantities of the disk component */
    writeDiskData(file, ndisk, maxLev, disk_info);
    //---------------------------------------------------------------------
  }/* if( addDisk ){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set particle distribution */
  //-----------------------------------------------------------------------
  /* nbody_particle *body; */
  /* allocParticleDataAoS((int)Ntot, &body); */
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
  //-----------------------------------------------------------------------
  /* create spherical particle distribution */
  Nuse = 0;
  for(int ii = 0; ii < skind; ii++){
    //---------------------------------------------------------------------
    /* distribute spheroid particles */
    if( cfg[ii].kind != CENTRALBH )
      cfg[ii].Ecut = distributeSpheroidParticles(&Nuse, body, (real)(cfg[ii].Mtot / (double)cfg[ii].num), cfg[ii], &prf[ii][2], fene[ii]);
    else{
      body.pos[Nuse].x = body.pos[Nuse].y = body.pos[Nuse].z = ZERO;      body.pos[Nuse].m   = (real)cfg[ii].Mtot;
      body.acc[Nuse].x = body.acc[Nuse].y = body.acc[Nuse].z = ZERO;      body.acc[Nuse].pot = ZERO;
#ifdef  BLOCK_TIME_STEP
      body.vel[Nuse].x = body.vel[Nuse].y = body.vel[Nuse].z = ZERO;
#else///BLOCK_TIME_STEP
      body.vx[Nuse] = body.vy[Nuse] = body.vz[Nuse] = ZERO;
#endif//BLOCK_TIME_STEP
      body.idx[Nuse] = Nuse;
      Nuse++;
    }/* else{ */
    //---------------------------------------------------------------------
    /* shift center-of-mass, remove bulk motion */
    shiftCenter(cfg[ii].num, Nuse - cfg[ii].num, body, cfg[ii]);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < skind; ii++){ */
  //-----------------------------------------------------------------------
  /* add disk component if required */
  if( addDisk ){
    //---------------------------------------------------------------------
#ifdef  CHECK_OSTRIKER_PEEBLES_CRITERION
    double Krand = 0.0;
    for(ulong ii = 0; ii < Nuse; ii++){
#ifdef  BLOCK_TIME_STEP
      Krand += body.vel[ii].x * body.vel[ii].x + body.vel[ii].y * body.vel[ii].y + body.vel[ii].z * body.vel[ii].z;
#else///BLOCK_TIME_STEP
      Krand += body.vx[ii] * body.vx[ii] + body.vy[ii] * body.vy[ii] + body.vz[ii] * body.vz[ii];
#endif//BLOCK_TIME_STEP
    }/* for(ulong ii = 0; ii < Nuse; ii++){ */
    for(int ii = 0; ii < ndisk; ii++)
      disk_info[ii].Krand_sph = Krand;
#endif//CHECK_OSTRIKER_PEEBLES_CRITERION
    //---------------------------------------------------------------------
    for(int ii = 0; ii < ndisk; ii++){
      //-------------------------------------------------------------------
      /* distribute disk particles */
      distributeDiskParticles(&Nuse, body, (real)(disk_info[ii].cfg->Mtot / (double)disk_info[ii].cfg->num), maxLev, disk_info[ii]);
      //-------------------------------------------------------------------
      /* shift center-of-mass, remove bulk motion */
      shiftCenter(disk_info[ii].cfg->num, Nuse - disk_info[ii].cfg->num, body, *disk_info[ii].cfg);
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < ndisk; ii++){ */
    //---------------------------------------------------------------------
  }/* if( addDisk ){ */
  //-----------------------------------------------------------------------
  if( skind == 0 )
    for(int ii = 0; ii < 4 + NRADBIN; ii++){
      prf[0][ii].enc_tot  = prf[ 0][ii].enc;
      for(int jj = 1; jj < ndisk; jj++)
      prf[0][ii].enc_tot += prf[jj][ii].enc;
    }/* for(int ii = 0; ii < 4 + NRADBIN; ii++){ */
  //-----------------------------------------------------------------------
  /* write fundamental information */
  outputFundamentalInformation(unit, kind, skind, cfg, prf, fene, Ntot, eps, snapshotInterval, ft, file);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  double time = 0.0;
  double  dt  = 0.0;
  int   last  = 1;
  ulong steps = 0;
  //-----------------------------------------------------------------------
  writeSettings(unit, Ntot, eps, eta, ft, snapshotInterval, saveInterval, file);
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
#ifdef  APPEND_ASCII_ICDATA
  FILE *fpascii;
  char asciifile[256];
  ulong Nhead = 0;
  for(int kk = 0; kk < kind; kk++){
    //---------------------------------------------------------------------
    sprintf(asciifile, "%s/%s%d_ascii.dat", DATAFOLDER, file, kk);
    fpascii = fopen(asciifile, "w");
    if( fpascii == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", asciifile);      }
    //---------------------------------------------------------------------
    for(ulong ii = 0; ii < cfg[kk].num; ii++)
      fprintf(fpascii, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", body[Nhead + ii].x, body[Nhead + ii].y, body[Nhead + ii].z, body[Nhead + ii].vx, body[Nhead + ii].vy, body[Nhead + ii].vz, body[Nhead + ii].m);
    Nhead += cfg[kk].num;
    //---------------------------------------------------------------------
    fclose(fpascii);
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < kind; kk++){ */
#endif//APPEND_ASCII_ICDATA
  //-----------------------------------------------------------------------
#ifdef  DUMPFILE_FOR_BONSAI
  FILE *fpbonsai;
  char bonsaifile[256];
  sprintf(bonsaifile, "%s/%s_bonsai.dat", DATAFOLDER, file);
  fpbonsai = fopen(bonsaifile, "wb");
  if( fpbonsai == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", bonsaifile);      }
  bool bonsaiSuccess = true;
  size_t bonsaiCount;
  /* header */
  struct dump {
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
  };
  typedef struct dump header;
  header bonsaiHeader;
  bonsaiHeader.time = time;
  bonsaiHeader.nbodies = (int)Ntot;
  bonsaiHeader.ndim = 3;
  bonsaiHeader.nsph = 0;
  bonsaiHeader.ndark = (int)Ntot;
  bonsaiHeader.nstar = 0;
  bonsaiCount = 1;  if( bonsaiCount != fwrite(&bonsaiHeader, sizeof(header), bonsaiCount, fpbonsai) )    bonsaiSuccess = false;
  /* main body */
  struct dark_particle {
    real mass;
    real pos[3];
    real vel[3];
    real eps;
    int idx;
  };
  for(ulong ii = 0; ii < Ntot; ii++){
    struct dark_particle tmp;
    tmp.mass   = body.pos[ii].m;
    tmp.pos[0] = body.pos[ii].x;
    tmp.pos[1] = body.pos[ii].y;
    tmp.pos[2] = body.pos[ii].z;
    tmp.vel[0] = body.vel[ii].x;
    tmp.vel[1] = body.vel[ii].y;
    tmp.vel[2] = body.vel[ii].z;
    tmp.eps    = eps;
    tmp.idx    = (int)body.idx[ii];
    bonsaiCount = 1;    if( bonsaiCount != fwrite(&tmp, sizeof(struct dark_particle), bonsaiCount, fpbonsai) )      bonsaiSuccess = false;
  }
  if( bonsaiSuccess != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", bonsaifile);  }
  fclose(fpbonsai);
#endif//DUMPFILE_FOR_BONSAI
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
  free(cfg);
  free(prf);
  free(_prf);
  free(fene);
  free(_fene);
  //-----------------------------------------------------------------------
  if( addDisk )
    freeDiskProfile(ndisk, disk_info,
		    disk_hor, disk_ver, disk_node_hor, disk_node_ver,
		    disk_pot, disk_rho, disk_rhoSum, disk_rhoTot, disk_dPhidR, disk_d2PhidR2,
		     disk_Sigma, disk_sigmaz, disk_enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
		     disk_zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
		     sph_rad, sph_rho, sph_enc,
		     spline_xx, spline_ff, spline_f2, spline_bp);
  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void outputFundamentalInformation
(const int unit, const int kind, const int skind, profile_cfg *cfg, profile **prf, dist_func **df, const ulong Ntot, const real eps, const double snapshotInterval, const double ft, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  FILE *fp;
  char filename[256], date[64];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output useful information for multi-component analysis */
  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
  //-----------------------------------------------------------------------
  fprintf(fp, "%d\n", unit);
  fprintf(fp, "%d\t%d\n", kind, skind);
  for(int ii = 0; ii < kind; ii++)
    fprintf(fp, "%zu\n", cfg[ii].num);
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output fundamental information of the particle distribution */
  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.info.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
  //-----------------------------------------------------------------------
  /* output global settings */
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
  //-----------------------------------------------------------------------
  ulong ihead = 0;
  for(int ii = 0; ii < kind; ii++){
    //---------------------------------------------------------------------
    /* output settings for individual component */
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
    case CENTRALBH:      fprintf(fp,           "Central massive black hole\n");      break;
    default:
      __KILL__(stderr, "ERROR: undefined profile ``%d'' specified\n", cfg[ii].kind);
      break;
    }/* switch( cfg[ii].kind ){ */
    fprintf(fp, "#############################################################################\n");
    fprintf(fp, "Total number of particles Ntot is %zu (about 2^%u)\n", cfg[ii].num, ilog2((int)cfg[ii].num));
    fprintf(fp, "Range of idx for the component is [%zu, %zu]\n", ihead, ihead + cfg[ii].num - 1);    ihead += cfg[ii].num;
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
    if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) ){
      if( cfg[ii].kind == SERSIC ){
	fprintf(fp, "Sersic index                   n is %e\n", cfg[ii].n_sersic);
	fprintf(fp, "Dimensionless scale factor     b is %e\n", cfg[ii].b_sersic);
      }/* if( cfg[ii].kind == SERSIC ){ */
      fprintf(fp, "Scale height of the component zd is %e (= %e %s)\n", cfg[ii].zd, cfg[ii].zd * length2astro, length_astro_unit_name);
      fprintf(fp, "Central surface density   Sigma0 is %e (= %e %s)\n", cfg[ii].Sigma0, cfg[ii].Sigma0 * col_density2astro, col_density_astro_unit_name);
      fprintf(fp, "Circular speed at scale radius   is %e (= %e %s)\n", cfg[ii].vcirc_Rd , cfg[ii].vcirc_Rd  * velocity2astro, velocity_astro_unit_name);
      fprintf(fp, "Maximum circular speed           is %e (= %e %s)\n", cfg[ii].vcirc_max, cfg[ii].vcirc_max * velocity2astro, velocity_astro_unit_name);
      fprintf(fp, "Circular speed is maximized      at %e (= %e %s)\n", cfg[ii].vcirc_Rd , cfg[ii].vcirc_Rd  *   length2astro,   length_astro_unit_name);
#ifdef  USE_ORIGINAL_VDISP_ESTIMATOR
      fprintf(fp, "Horizontal velocity dispersion   is %e of circular velocity or vertical velocity dispersion (maximum is used)\n", cfg[ii].vdisp_frac);
#else///USE_ORIGINAL_VDISP_ESTIMATOR
      fprintf(fp, "Horizontal velocity dispersion   is %e (= %e %s)\n", cfg[ii].vdispR0  , cfg[ii].vdispR0   * velocity2astro, velocity_astro_unit_name);
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
      fprintf(fp, "Vertical   velocity dispersion   is %e (= %e %s)\n", cfg[ii].vdispz0  , cfg[ii].vdispz0   * velocity2astro, velocity_astro_unit_name);
      fprintf(fp, "Toomre's Q-value at scale radius is %e\n", cfg[ii].toomre);
      fprintf(fp, "Minimum of Toomre's Q-value      is %e\n", cfg[ii].Qmin);
    }/* if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) ){ */
    if( cfg[ii].kind != CENTRALBH ){
      if( cfg[ii].cutoff ){
	fprintf(fp, "Cutoff radius of the component   is %e (= %e %s)\n", cfg[ii].rc      , cfg[ii].rc       * length2astro, length_astro_unit_name);
	fprintf(fp, "Cutoff  width of the component   is %e (= %e %s)\n", cfg[ii].rc_width, cfg[ii].rc_width * length2astro, length_astro_unit_name);
	fprintf(fp, "Cutoff radius over scale radius  is %e\n", cfg[ii].rc / cfg[ii].rs);
	fprintf(fp, "Cutoff  width over scale radius  is %e\n", cfg[ii].rc_width / cfg[ii].rs);
      }/* if( cfg[ii].cutoff ){ */
      fprintf(fp, "Cutoff energy of the component   is %e (= %e %s) (= %e G Mtot / rs)\n", cfg[ii].Ecut, cfg[ii].Ecut * senergy2astro, senergy_astro_unit_name, cfg[ii].Ecut / (newton * cfg[ii].Mtot / cfg[ii].rs));
    }/* if( cfg[ii].kind != CENTRALBH ){ */
    fprintf(fp, "#############################################################################\n");
    //---------------------------------------------------------------------
    /* estimate typical timescale */
    if( cfg[ii].kind != CENTRALBH ){
      //-------------------------------------------------------------------
      /* estimate enclosed mass within rs */
      double Ms;
      {
	int jj = 2;
	while( true ){
	  if( prf[0][jj].rad > cfg[ii].rs ){	    jj--;	    break;	  }
	  jj++;
	  if( jj == (NRADBIN + 1) )	    break;
	}/* while( true ){ */
	Ms = prf[0][jj].enc_tot;
      }
      //-------------------------------------------------------------------
      /* estimate dynamical time at scale length for each component */
      const double Ns = (double)Ntot * (Ms / Mtot);
      const double tff = M_PI_2 * cfg[ii].rs * sqrt(cfg[ii].rs / (2.0 * (double)newton * Ms));
      const double t2r = tff * Ns / (32.0 * log(cfg[ii].rs / (double)eps));
      double trot;
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) )
	trot = 2.0 * M_PI * cfg[ii].rs / cfg[ii].vcirc_Rd;
      fprintf(fp, "Total number of particles within the scale length is       %e\n", Ns);
      fprintf(fp, "Enclosed mass of all components within the scale length is %e (= %e %s)\n",  Ms,  Ms * mass2astro, mass_astro_unit_name);
      fprintf(fp, "Free-fall time at the scale length                      is %e (= %e %s)\n", tff, tff * time2astro, time_astro_unit_name);
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) )
	fprintf(fp, "Rotation time scale at the scale length                 is %e (= %e x tff = %e %s)\n", trot, trot / tff, trot * time2astro, time_astro_unit_name);
      fprintf(fp, "Two-body relaxation time at the scale length            is %e (= %e %s)\n", t2r, t2r * time2astro, time_astro_unit_name);
      fprintf(fp, "#############################################################################\n");
      fprintf(fp, "Snapshot interval in the unit of free-fall time           is %e\n", (double)snapshotInterval / tff);
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) )
	fprintf(fp, "Snapshot interval in the unit of rotation time scale      is %e\n", (double)snapshotInterval / trot);
      fprintf(fp, "Snapshot interval in the unit of two-body relaxation time is %e\n", (double)snapshotInterval / t2r);
      fprintf(fp, "#############################################################################\n");
      fprintf(fp, "Final time of the simulation in the unit of free-fall time           is %e\n", (double)ft / tff);
      if( (cfg[ii].kind == EXP_DISK) || (cfg[ii].kind == SERSIC) )
	fprintf(fp, "Final time of the simulation in the unit of rotation time scale      is %e\n", (double)ft / trot);
      fprintf(fp, "Final time of the simulation in the unit of two-body relaxation time is %e\n", (double)ft / t2r);
      fprintf(fp, "#############################################################################\n");
      fprintf(fp, "#############################################################################\n");
      //-------------------------------------------------------------------
    }/* if( cfg[ii].kind != CENTRALBH ){ */
    //---------------------------------------------------------------------
    fprintf(fp, "\n");
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < kind; ii++){ */
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output fundamental profile of the particle distribution */
  //-----------------------------------------------------------------------
  real *tmp_rad;  tmp_rad = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_rad == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_rad\n");  }
  real *tmp_rho;  tmp_rho = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_rho\n");  }
  real *tmp_enc;  tmp_enc = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_enc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_enc\n");  }
  real *tmp_psi;  tmp_psi = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_psi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_psi\n");  }
  real *tmp_tff;  tmp_tff = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_tff == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_tff\n");  }
  real *tmp_t2r;  tmp_t2r = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_t2r == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_t2r\n");  }
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  real *tmp_Sig;  tmp_Sig = (real *)malloc(NRADBIN * sizeof(real));  if( tmp_Sig == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_Sig\n");  }
#endif//MAKE_COLUMN_DENSITY_PROFILE
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
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
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < kind; kk++){
    //---------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
    char grp[16];    sprintf(grp, "data%d", kk);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else///USE_HDF5_FORMAT
    sprintf(filename, "%s/%s.profile.%d.dat", DATAFOLDER, file, kk);
    fp = fopen(filename, "wb");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }
#endif//USE_HDF5_FORMAT
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      tmp_rad[ii] = (real)(prf[kk][ii].rad *  length2astro);
      tmp_rho[ii] = (real)(prf[kk][ii].rho * density2astro);
      tmp_enc[ii] = (real)(prf[kk][ii].enc *    mass2astro);
      tmp_psi[ii] = (real)(prf[kk][ii].psi * senergy2astro);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
      tmp_Sig[ii] = (real)(prf[kk][ii].Sigma * col_density2astro);
#endif//MAKE_COLUMN_DENSITY_PROFILE
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */
    //---------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
    hsize_t dims = NRADBIN;
    dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
    /* write radius */
    dataset = H5Dcreate(group, "rad", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rad));
    chkHDF5err(H5Dclose(dataset));
    /* write density */
    dataset = H5Dcreate(group, "rho", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rho));
    chkHDF5err(H5Dclose(dataset));
    /* write enclosed mass */
    dataset = H5Dcreate(group, "enc", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_enc));
    chkHDF5err(H5Dclose(dataset));
    /* write potential */
    dataset = H5Dcreate(group, "Psi", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_psi));
    chkHDF5err(H5Dclose(dataset));
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
    /* write column density */
    dataset = H5Dcreate(group, "Sigma", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_Sig));
    chkHDF5err(H5Dclose(dataset));
#endif//MAKE_COLUMN_DENSITY_PROFILE
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------
    /* write attribute data */
    //---------------------------------------------------------------------
    /* create the data space for the attribute */
    hsize_t attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    //---------------------------------------------------------------------
    /* write # of arrays */
    int nradbin = NRADBIN;
    attribute = H5Acreate(group, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nradbin));
    chkHDF5err(H5Aclose(attribute));
    /* write scale radius */
    attribute = H5Acreate(group, "rs", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg[kk].rs));
    chkHDF5err(H5Aclose(attribute));
    /* write total mass */
    attribute = H5Acreate(group, "Mtot", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg[kk].Mtot));
    chkHDF5err(H5Aclose(attribute));
    /* profile ID */
    attribute = H5Acreate(group, "profile ID", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cfg[kk].kind));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
    //---------------------------------------------------------------------
#else///USE_HDF5_FORMAT
    int nradbin = NRADBIN;
    bool success = true;
    success &= (fwrite(&nradbin, sizeof( int),       1, fp) ==       1);
    success &= (fwrite( tmp_rad, sizeof(real), NRADBIN, fp) == NRADBIN);
    success &= (fwrite( tmp_rho, sizeof(real), NRADBIN, fp) == NRADBIN);
    success &= (fwrite( tmp_enc, sizeof(real), NRADBIN, fp) == NRADBIN);
    success &= (fwrite( tmp_psi, sizeof(real), NRADBIN, fp) == NRADBIN);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
    success &= (fwrite( tmp_Sig, sizeof(real), NRADBIN, fp) == NRADBIN);
#endif//MAKE_COLUMN_DENSITY_PROFILE
    if( !success ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
    fclose(fp);
#endif//USE_HDF5_FORMAT
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < kind; kk++){ */
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  /* evaluate typical timescale */
#pragma omp parallel for
  for(int ii = 0; ii < NRADBIN; ii++){
    const double rad = prf[0][ii].rad;
    double enc = 0.0;
    for(int jj = 0; jj < kind; jj++)
      enc += prf[jj][ii].enc;
    const double tff = M_PI_2 * rad * sqrt(rad / (2.0 * (double)newton * enc));
    double t2r = tff * ((double)Ntot * (enc / Mtot)) / (32.0 * log(enc / (double)eps));
    if( t2r < 0.0 )
      t2r = 0.0;
    tmp_rad[ii] = (real)(rad * length2astro);
    tmp_tff[ii] = (real)(tff *   time2astro);
    tmp_t2r[ii] = (real)(t2r *   time2astro);
  }/* for(int ii = 0; ii < NRADBIN; ii++){ */
  //-----------------------------------------------------------------------
  /* write typical timescale */
  hid_t timeScale = H5Gcreate(target, "time_scale", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t dims = NRADBIN;
  dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  chkHDF5err(H5Pset_chunk(property, 1, &cdims));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
  /* write radius */
  dataset = H5Dcreate(timeScale, "rad", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rad));
  chkHDF5err(H5Dclose(dataset));
  /* write free-fall time */
  dataset = H5Dcreate(timeScale, "free_fall_time", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_tff));
  chkHDF5err(H5Dclose(dataset));
  /* write relaxation time */
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
  //-----------------------------------------------------------------------
  /* write attribute data */
  //-----------------------------------------------------------------------
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
  /* write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &kind));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
  attribute = H5Acreate(target,  "length_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,  length_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, density_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,    "mass_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,    mass_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "senergy_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, senergy_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target,    "time_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format,    time_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "col_density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, col_density_astro_unit_name));  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* close the file */
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
#endif//USE_HDF5_FORMAT
  free(tmp_rad);
  free(tmp_rho);
  free(tmp_enc);
  free(tmp_psi);
  free(tmp_tff);
  free(tmp_t2r);
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  free(tmp_Sig);
#endif//MAKE_COLUMN_DENSITY_PROFILE
  //-----------------------------------------------------------------------
#ifdef  OUTPUT_ASCII_PROFILE
  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.profile.txt", DATAFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
  //-----------------------------------------------------------------------
  fprintf(fp, "#r\trho(r)\tM(r)\tPhi(r)\n");
  fprintf(fp, "#\tgenerated on %s", date);
  fprintf(fp, "#number of conponents is %d\n", kind);
  fprintf(fp, "#format is total");
  for(int kk = 0; kk < kind; kk++)
    fprintf(fp, ", %d-th", kk);
  fprintf(fp, "\n");
  const int nskip = NRADBIN / N_PRINT_LINES_ASCII;
  for(int ii = 2; ii < 2 + NRADBIN; ii += nskip){
    fprintf(fp, "%e", prf[0][ii].rad);
    fprintf(fp, "\t%e",  prf[0][ii].rho_tot);    for(int kk = 0; kk < kind; kk++)      fprintf(fp, "\t%e",  prf[kk][ii].rho);
    fprintf(fp, "\t%e",  prf[0][ii].enc_tot);    for(int kk = 0; kk < kind; kk++)      fprintf(fp, "\t%e",  prf[kk][ii].enc);
    fprintf(fp, "\t%e", -prf[0][ii].psi_tot);    for(int kk = 0; kk < kind; kk++)      fprintf(fp, "\t%e", -prf[kk][ii].psi);
    fprintf(fp, "\n");
  }/* for(int ii = 2; ii < 2 + NRADBIN; ii += nskip){ */
  fclose(fp);
  //-----------------------------------------------------------------------
#endif//OUTPUT_ASCII_PROFILE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output distribution function of the particle distribution */
  //-----------------------------------------------------------------------
  real *tmp_ene;  tmp_ene = (real *)malloc(NENEBIN * sizeof(real));  if( tmp_ene == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_ene\n");  }
  real *tmp_val;  tmp_val = (real *)malloc(NENEBIN * sizeof(real));  if( tmp_val == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_val\n");  }
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "df");
  target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < skind; kk++){
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NENEBIN; ii++){
      tmp_ene[ii] = (real)((double)df[kk][ii].ene * senergy2astro);
      tmp_val[ii] =                df[kk][ii].val;
    }/* for(int ii = 0; ii < NENEBIN; ii++){ */
    //---------------------------------------------------------------------
    char grp[16];    sprintf(grp,  "series%d", kk);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //---------------------------------------------------------------------
    hsize_t dims = NENEBIN;
    dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
    /* write energy */
    dataset = H5Dcreate(group, "energy", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_ene));
    chkHDF5err(H5Dclose(dataset));
    /* write energy */
    dataset = H5Dcreate(group, "DF", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_val));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------
    /* write attribute data */
    //---------------------------------------------------------------------
    /* create the data space for the attribute */
    hsize_t attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    //---------------------------------------------------------------------
    /* write # of arrays */
    int nenebin = NENEBIN;
    attribute = H5Acreate(group, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nenebin));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < skind; kk++){ */
  //-----------------------------------------------------------------------
  /* write attribute data */
  //-----------------------------------------------------------------------
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
  /* write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skind));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  attribute = H5Acreate(target, "senergy_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, senergy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));
  //-----------------------------------------------------------------------
  /* write flag about DOUBLE_PRECISION */
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* close the file */
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
#else///USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < skind; kk++){
    //---------------------------------------------------------------------
    sprintf(filename, "%s/%s.df.%d.dat", DATAFOLDER, file, kk);
    fp = fopen(filename, "wb");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NENEBIN; ii++){
      tmp_ene[ii] = df[kk][ii].ene;
      tmp_val[ii] = df[kk][ii].val;
    }/* for(int ii = 0; ii < NENEBIN; ii++){ */
    //---------------------------------------------------------------------
    int nenebin = NENEBIN;
    bool success = true;
    success &= (fwrite(&nenebin, sizeof( int),       1, fp) ==       1);
    success &= (fwrite( tmp_ene, sizeof(real), NENEBIN, fp) == NENEBIN);
    success &= (fwrite( tmp_val, sizeof(real), NENEBIN, fp) == NENEBIN);
    if( !success ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
    fclose(fp);
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < skind; kk++){ */
  //-----------------------------------------------------------------------
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  free(tmp_ene);
  free(tmp_val);
  //-----------------------------------------------------------------------
#ifdef  OUTPUT_ASCII_PROFILE
  //-----------------------------------------------------------------------
  if( skind > 0 ){
    //---------------------------------------------------------------------
    sprintf(filename, "%s/%s.df.txt", DATAFOLDER, file);
    fp = fopen(filename, "w");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }
    //---------------------------------------------------------------------
    fprintf(fp, "#relE\tDF(relE)\n");
    fprintf(fp, "#\tgenerated on %s", date);
    fprintf(fp, "#number of spherical conponent(s) is %d\n", skind);
    fprintf(fp, "#format is energy\t%d-th", 0);
    for(int kk = 1; kk < skind; kk++)
      fprintf(fp, "\t%d-th", kk);
    fprintf(fp, "\n");
    const int nskip_ene = NENEBIN / N_PRINT_LINES_ASCII;
    for(int ii = 0; ii < NENEBIN; ii += nskip_ene){
      fprintf(fp, "%e", df[0][ii].ene);
      for(int kk = 0; kk < skind; kk++)
	fprintf(fp, "\t%e", df[kk][ii].val);
      fprintf(fp, "\n");
    }/* for(int ii = 0; ii < NENEBIN; ii += nskip_ene){ */
    fclose(fp);
    //---------------------------------------------------------------------
  }/* if( skind > 0 ){ */
  //-----------------------------------------------------------------------
#endif//OUTPUT_ASCII_PROFILE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* real * restrict _????? are output arrays */
//-------------------------------------------------------------------------
static void evaluateDiskProperties
(disk_data *disk_info, const int diskID, const int lev, const int ihead, const int itail,
 real * restrict _vcirc, real * restrict _sigmap, real * restrict _sigmaR, real * restrict _kappa, real * restrict _Omega, real * restrict _toomre, real * restrict _lambda)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifndef USE_ORIGINAL_VDISP_ESTIMATOR
  const double invRd = 1.0 / disk_info[diskID].cfg->rs;
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
  /* /\* disk_info[diskID].cfg->vcirc_max   = -1.0; *\/ */
  /* /\* disk_info[diskID].cfg->vcirc_max_R = -1.0; *\/ */
  /* /\* disk_info[diskID].cfg->Qmin = DBL_MAX; *\/ */
  /* bool passed = false; */
  bool passed = disk_info[diskID].cfg->passed;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int ii = ihead; ii < itail + 1; ii++){
    //---------------------------------------------------------------------
    /* evaluate epicyclic quantities and circular speed */
#ifndef USE_POTENTIAL_SCALING_SCHEME
    const double Omega = sqrt(disk_info[diskID]. dPhidR [INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] / disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);
    const double kappa = sqrt(disk_info[diskID].d2PhidR2[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] + 3.0 * Omega * Omega);
#else///USE_POTENTIAL_SCALING_SCHEME
    const double Omega = sqrt(disk_info[diskID]. dPhidR [INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] / disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);
    const double kappa = sqrt(disk_info[diskID].d2PhidR2[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] + 3.0 * Omega * Omega);
#endif//USE_POTENTIAL_SCALING_SCHEME
    const double vcirc = disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] * Omega;
    //---------------------------------------------------------------------
    /* evaluate Toomre's Q-value */
#ifdef  USE_ORIGINAL_VDISP_ESTIMATOR
    const double sigmap = DISK_PERP_VDISP(disk_info[diskID].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], vcirc, disk_info[diskID].cfg->vdisp_frac);
    const double sigmaR = sigmap * 2.0 * Omega / (DBL_MIN + kappa);
#else///USE_ORIGINAL_VDISP_ESTIMATOR
    const double sigmaR = DISK_RADIAL_VDISP(disk_info[diskID].cfg->vdispR0, disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], invRd);
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
    const double Sigma = disk_info[diskID].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
    const double toomre = sigmaR * kappa / (DBL_MIN + 3.36 * newton * Sigma);
    const double lambda = 4.0 * M_PI * M_PI * newton * Sigma / (DBL_MIN + kappa * kappa);
    //---------------------------------------------------------------------
    /* find the maximum circular speed */
    if( vcirc > disk_info[diskID].cfg->vcirc_max ){
      disk_info[diskID].cfg->vcirc_max   = vcirc;
      disk_info[diskID].cfg->vcirc_max_R = disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)];
    }/* if( vcirc > disk_info[diskID].cfg->vcirc_max ){ */
    //---------------------------------------------------------------------
    /* find the minimum Toomre's Q-value */
    /* if( (disk_info[diskID].Sigma[ii] > 1.0e-4 * disk_info[diskID].Sigma[0]) && (toomre > DBL_EPSILON) && (toomre < disk_info[diskID].cfg->Qmin) ) */
    if( toomre < disk_info[diskID].cfg->Qmin )
      disk_info[diskID].cfg->Qmin = toomre;
    //---------------------------------------------------------------------
    /* find the circular speed and Toomre's Q-value at the scale radius */
    if( !passed ){
      disk_info[diskID].cfg->vcirc_Rd = vcirc;
      disk_info[diskID].cfg->toomre   = toomre;
      if( disk_info[diskID].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] > disk_info[diskID].cfg->rs )
	passed = true;
    }/* if( !passed ){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* memorize calculated values */
    //---------------------------------------------------------------------
    const int jj = ii - ihead;
    _vcirc [jj] = (real)(vcirc	* velocity2astro);    /* if( isinf(_vcirc [jj]) == 1 )      _vcirc [jj] = REAL_MAX; */
    _sigmap[jj] = (real)(sigmap * velocity2astro);    /* if( isinf(_sigmap[jj]) == 1 )      _sigmap[jj] = REAL_MAX; */
    _sigmaR[jj] = (real)(sigmaR * velocity2astro);    if( isinf(_sigmaR[jj]) == 1 )      _sigmaR[jj] = REAL_MAX;
    _kappa [jj] = (real)(kappa	/     time2astro);    /* if( isinf(_kappa [jj]) == 1 )      _kappa [jj] = REAL_MAX; */
    _Omega [jj] = (real)(Omega	/     time2astro);    /* if( isinf(_Omega [jj]) == 1 )      _Omega [jj] = REAL_MAX; */
    _lambda[jj] = (real)(lambda *   length2astro);    if( isinf(_lambda[jj]) == 1 )      _lambda[jj] = REAL_MAX;
    _toomre[jj] = (real) toomre                  ;    if( isinf(_toomre[jj]) == 1 )      _toomre[jj] = REAL_MAX;
    //---------------------------------------------------------------------
  }/* for(int ii = ihead; ii < itail + 1; ii++){ */
  //-----------------------------------------------------------------------
  disk_info[diskID].cfg->passed = passed;
  //-----------------------------------------------------------------------

  /* //----------------------------------------------------------------------- */
  /* /\* output profile data in ASCII *\/ */
  /* //----------------------------------------------------------------------- */
  /* FILE *fp; */
  /* char filename[256]; */
  /* sprintf(filename, "%s/%s.diskhor.%d.txt", DATAFOLDER, file, diskID); */
  /* fp = fopen(filename, "w"); */
  /* if( fp == NULL ){ */
  /*   __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename); */
  /* } */
  /* fprintf(fp, "#R (%s)\tSigma(R) (%s)\tfrac(R)\tv_c(R) (%s)\tsigma_p(R) (%s)\tsigma_R(R) (%s)\tsigma_z(R) (%s)\tkappa(R) (1 / %s)\tOmega(R) (1/%s)\tToomre-Q(R)\tlambda_crit(R) (%s)\n", */
  /* 	  length_astro_unit_name, col_density_astro_unit_name, velocity_astro_unit_name, velocity_astro_unit_name, velocity_astro_unit_name, velocity_astro_unit_name, time_astro_unit_name, time_astro_unit_name, length_astro_unit_name); */
  /* //----------------------------------------------------------------------- */
  /* for(int ii = 0; ii < NDISKBIN_HOR; ii++) */
  /*   fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", disk_info[diskID].hor[ii] * length2astro, disk_info[diskID].Sigma[ii] * col_density2astro, rhoFrac[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, 0)], _vcirc[ii], _sigmap[ii], _sigmaR[ii], disk_info[diskID].sigmaz[ii] * velocity2astro, _kappa[ii], _Omega[ii], _toomre[ii], _lambda[ii]); */
  /* //----------------------------------------------------------------------- */
  /* fclose(fp); */
  /* //----------------------------------------------------------------------- */


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeDiskData(char *file, const int ndisk, const int maxLev, disk_data *disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* arrays for 2D-plot */
#ifdef  USE_HDF5_FORMAT
  real *node_RR;  node_RR = (real *)malloc((NDISKBIN_HOR + 1) * sizeof(real));  if( node_RR == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node_RR\n");  }
  real *node_zz;  node_zz = (real *)malloc((NDISKBIN_VER + 1) * sizeof(real));  if( node_zz == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node_zz\n");  }
#endif//USE_HDF5_FORMAT
  real *tmp_rho;  tmp_rho = (real *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(real));  if( tmp_rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_rho\n");  }
  real *rhoFrac;  rhoFrac = (real *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(real));  if( rhoFrac == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rhoFrac\n");  }
  real *tmp_Phi;  tmp_Phi = (real *)malloc(NDISKBIN_HOR * NDISKBIN_VER * sizeof(real));  if( tmp_Phi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_Phi\n");  }
  //-----------------------------------------------------------------------
  /* arrays for 1D-plot */
#ifdef  USE_HDF5_FORMAT
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
#else///USE_HDF5_FORMAT
  real *tmp_ver;  tmp_ver = (real *)malloc(NDISKBIN_VER * sizeof(real));  if( tmp_ver == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_ver\n");  }
  real *tmp_hor;  tmp_hor = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( tmp_hor == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_hor\n");  }
  real *tmp_sig;  tmp_sig = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( tmp_sig == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_sig\n");  }
  real *tmp_Sig;  tmp_Sig = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( tmp_Sig == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_Sig\n");  }
  real *_vcirc ;  _vcirc  = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( _vcirc  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _vcirc\n" );  }
  real *_sigmap;  _sigmap = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( _sigmap == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _sigmap\n");  }
  real *_sigmaR;  _sigmaR = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( _sigmaR == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _sigmaR\n");  }
  real *_kappa ;  _kappa  = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( _kappa  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _kappa\n" );  }
  real *_Omega ;  _Omega  = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( _Omega  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _Omega\n" );  }
  real *_lambda;  _lambda = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( _lambda == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _lambda\n");  }
  real *_toomre;  _toomre = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( _toomre == NULL ){    __KILL__(stderr, "ERROR: failure to allocate _toomre\n");  }
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  real *tmp__zd;  tmp__zd = (real *)malloc(NDISKBIN_HOR * sizeof(real));  if( tmp__zd == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp__zd\n");  }
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  if( ndisk > 1 )
    for(int lev = 0; lev < maxLev; lev++)
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR * NDISKBIN_VER; ii++){
	disk[0].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)] = (*disk[0].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)];
	for(int jj = 1; jj < ndisk; jj++)
	  disk[0].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)] += (*disk[jj].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)];
	disk[0].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)] = 1.0 / (DBL_MIN + disk[0].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)]);
      }/* for(int ii = 0; ii < NDISKBIN_HOR * NDISKBIN_VER; ii++){ */
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  char filename[128];
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "disk");
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
  hsize_t cdims[2];
  /* hsize_t cdims[2] = {128, 128 * szip_pixels_per_block}; */
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    /* output in HDF5 format */
    //---------------------------------------------------------------------
    char grp[16];    sprintf(grp, "data%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write two-dimensional data */
    //---------------------------------------------------------------------
    for(int lev = 0; lev < maxLev; lev++){
      //-------------------------------------------------------------------
      char subgrp[16];      sprintf(subgrp, "patch%d", lev);
      hid_t patch = H5Gcreate(group, subgrp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* write data */
      //-------------------------------------------------------------------
      if( ndisk > 1 )
#pragma omp parallel for
	for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++)
	  rhoFrac[jj] = (real)((*disk[ii].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * disk[ii].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)]);
#pragma omp parallel for
      for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++){
	tmp_rho[jj] = (real)((*disk[ii].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * density2astro);
	tmp_Phi[jj] = (real)(  disk[ii].pot [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * senergy2astro);
      }/* for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++){ */
      //-------------------------------------------------------------------
      hsize_t dims[2] = {NDISKBIN_HOR, NDISKBIN_VER};
      dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
      property = H5Pcreate(H5P_DATASET_CREATE);
      cdims[0] = (dims[0] <  128                          ) ? (dims[0]) : (128);
      cdims[1] = (dims[1] < (128 * szip_pixels_per_block) ) ? (dims[1]) : (128 * szip_pixels_per_block);
      chkHDF5err(H5Pset_chunk(property, 2, cdims));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
      /* write density */
      dataset = H5Dcreate(patch, "rho", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_rho));
      chkHDF5err(H5Dclose(dataset));
      /* write potential */
      dataset = H5Dcreate(patch, "Phi", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_Phi));
      chkHDF5err(H5Dclose(dataset));
      /* write fraction */
      if( ndisk > 1 ){
	dataset = H5Dcreate(patch, "fraction", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
	chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rhoFrac));
	chkHDF5err(H5Dclose(dataset));
      }/* if( ndisk > 1 ){ */
#ifdef  USE_SZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* prepare horizontal and vertical axes */
      //-------------------------------------------------------------------
      double width[2];
      width[0] = ldexp(disk[ii].hh, -lev) * length2astro;
      width[1] = ldexp(disk[ii].hh, -lev) * length2astro;
      //-------------------------------------------------------------------
      /* note: VisIt requires position of edges of the grid */
#pragma omp parallel
      {
#pragma omp for nowait
	for(int jj = 0; jj < NDISKBIN_HOR; jj++)
	  tmp_hor[jj] = (real)(disk[ii].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] * length2astro);
#pragma omp for nowait
	for(int jj = 0; jj < NDISKBIN_VER; jj++)
	  tmp_ver[jj] = (real)(disk[ii].ver[INDEX2D(maxLev, NDISKBIN_VER, lev, jj)] * length2astro);
#pragma omp for nowait
	for(int jj = 0; jj < NDISKBIN_HOR + 1; jj++)
	  node_RR[jj] = (real)(disk[ii].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, jj)] * length2astro);
#pragma omp for nowait
	for(int jj = 0; jj < NDISKBIN_VER + 1; jj++)
	  node_zz[jj] = (real)(disk[ii].node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev, jj)] * length2astro);
      }
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* write horizontal axis */
      //-------------------------------------------------------------------
      dataspace = H5Screate_simple(1, &dims[0], NULL);
#ifdef  USE_SZIP_COMPRESSION
      property = H5Pcreate(H5P_DATASET_CREATE);
      cdims[0] = (dims[0] < (128 * szip_pixels_per_block) ) ? (dims[0]) : (128 * szip_pixels_per_block);
      chkHDF5err(H5Pset_chunk(property, 1, &cdims[0]));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
      /* write horizontal position */
      dataset = H5Dcreate(patch, "radius", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_hor));
      chkHDF5err(H5Dclose(dataset));
      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));
      //-------------------------------------------------------------------
      dims[0]++;
      dataspace = H5Screate_simple(1, &dims[0], NULL);
#ifdef  USE_SZIP_COMPRESSION
      property = H5Pcreate(H5P_DATASET_CREATE);
      cdims[0] = (dims[0] < (128 * szip_pixels_per_block) ) ? (dims[0]) : (128 * szip_pixels_per_block);
      chkHDF5err(H5Pset_chunk(property, 1, &cdims[0]));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
      /* write horizontal position */
      dataset = H5Dcreate(patch, "node_RR", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, node_RR));
      chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* write vertical axis */
      //-------------------------------------------------------------------
      dataspace = H5Screate_simple(1, &dims[1], NULL);
#ifdef  USE_SZIP_COMPRESSION
      property = H5Pcreate(H5P_DATASET_CREATE);
      chkHDF5err(H5Pset_chunk(property, 1, &cdims[1]));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
      /* write vertical position */
      dataset = H5Dcreate(patch, "height", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_ver));
      chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));
      //-------------------------------------------------------------------
      dims[1]++;
      dataspace = H5Screate_simple(1, &dims[1], NULL);
#ifdef  USE_SZIP_COMPRESSION
      property = H5Pcreate(H5P_DATASET_CREATE);
      cdims[1] = (dims[1] < (128 * szip_pixels_per_block) ) ? (dims[1]) : (128 * szip_pixels_per_block);
      chkHDF5err(H5Pset_chunk(property, 1, &cdims[1]));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
      /* write vertical position */
      dataset = H5Dcreate(patch, "node_zz", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, node_zz));
      chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* write attribute data */
      //-------------------------------------------------------------------
      /* create the data space for the attribute */
      hsize_t attr_dims = 2;
      dataspace = H5Screate_simple(1, &attr_dims, NULL);
      hid_t attribute;
      //-------------------------------------------------------------------
      /* write # of arrays */
      int narray[2] = {NDISKBIN_HOR, NDISKBIN_VER};
      attribute = H5Acreate(patch, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, narray));
      chkHDF5err(H5Aclose(attribute));
      /* write the head index which corresponding to the regular grid having the same resolution with this level (i.e., (NDISKBIN_HOR << lev) * (NDISKBIN_VER << lev) grid points) */
      narray[0] = 0;
      narray[1] = 0;
      attribute = H5Acreate(patch, "head", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, narray));
      chkHDF5err(H5Aclose(attribute));
      /* write mesh sizes */
      attribute = H5Acreate(patch, "width", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, width));
      chkHDF5err(H5Aclose(attribute));
      //-------------------------------------------------------------------
      chkHDF5err(H5Sclose(dataspace));
      attr_dims = 1;
      dataspace = H5Screate_simple(1, &attr_dims, NULL);
      //-------------------------------------------------------------------
      /* Nested level */
      attribute = H5Acreate(patch, "level", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &lev));
      chkHDF5err(H5Aclose(attribute));
      /* Parent patch */
      int parent = (lev > 0) ? (lev - 1) : 0;
      attribute = H5Acreate(patch, "parent", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &parent));
      chkHDF5err(H5Aclose(attribute));
      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      chkHDF5err(H5Gclose(patch));
      //-------------------------------------------------------------------
    }/* for(int lev = 0; lev < maxLev; lev++){ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* prepare one-dimensional data */
    //---------------------------------------------------------------------
    /* data preparation in the finest grid */
    evaluateDiskProperties(disk, ii, maxLev - 1, 0, NDISKBIN_HOR - 1, _vcirc, _sigmap, _sigmaR, _kappa, _Omega, _toomre, _lambda);
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_HOR; jj++){
      //-------------------------------------------------------------------
      tmp_hor[jj] = (real)(disk[ii].hor   [INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] *      length2astro);
      tmp_sig[jj] = (real)(disk[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] *    velocity2astro);
      tmp_Sig[jj] = (real)(disk[ii].Sigma [INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] * col_density2astro);
#ifdef	ENABLE_VARIABLE_SCALE_HEIGHT
      tmp__zd[jj] = (real)(disk[ii].zd    [INDEX2D(maxLev, NDISKBIN_HOR, maxLev - 1, jj)] *      length2astro);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < NDISKBIN_HOR; jj++){ */
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
      tmp_ver[jj] = (real)(disk[ii].ver   [INDEX2D(maxLev, NDISKBIN_VER, maxLev - 1, jj)] *      length2astro);
    //---------------------------------------------------------------------
    /* data preparation in coarser grids */
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
      //-------------------------------------------------------------------
#pragma omp for nowait
      for(int jj = 0; jj < (NDISKBIN_HOR >> 1); jj++){
	tmp_hor[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = (real)(disk[ii].hor   [INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] *      length2astro);
	tmp_sig[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = (real)(disk[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] *    velocity2astro);
	tmp_Sig[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = (real)(disk[ii].Sigma [INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] * col_density2astro);
#ifdef	ENABLE_VARIABLE_SCALE_HEIGHT
	tmp__zd[(NDISKBIN_HOR >> 1) * (maxLev - lev) + jj] = (real)(disk[ii].zd    [INDEX2D(maxLev, NDISKBIN_HOR, lev, (NDISKBIN_HOR >> 1) + jj)] *      length2astro);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      }/* for(int jj = 0; jj < NDISKBIN_HOR; jj++){ */
#pragma omp for nowait
      for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++)
	tmp_ver[(NDISKBIN_VER >> 1) * (maxLev - lev) + jj] = (real)(disk[ii].ver   [INDEX2D(maxLev, NDISKBIN_VER, lev, (NDISKBIN_VER >> 1) + jj)] *      length2astro);
      //-------------------------------------------------------------------
    }/* for(int lev = maxLev - 2; lev >= 0; lev--){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write one-dimensional data */
    //---------------------------------------------------------------------
    char subgrp[16];
    sprintf(subgrp, "1D_data");
    hid_t patch = H5Gcreate(group, subgrp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //---------------------------------------------------------------------
    hsize_t dims[2] = {(NDISKBIN_HOR >> 1) * (maxLev + 1), (NDISKBIN_VER >> 1) * (maxLev + 1)};
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write horizontal variables */
    //---------------------------------------------------------------------
    dataspace = H5Screate_simple(1, &dims[0], NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    cdims[0] = (dims[0] < (128 * szip_pixels_per_block) ) ? (dims[0]) : (128 * szip_pixels_per_block);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims[0]));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
      /* write horizontal position */
    dataset = H5Dcreate(patch, "radius", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_hor));
    chkHDF5err(H5Dclose(dataset));
    /* write velocity dispersion in z-direction */
    dataset = H5Dcreate(patch, "sigmaz", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_sig));
    chkHDF5err(H5Dclose(dataset));
    /* write column density */
    dataset = H5Dcreate(patch, "Sigma", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_Sig));
    chkHDF5err(H5Dclose(dataset));
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
    /* write scale height */
    dataset = H5Dcreate(patch, "zd", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp__zd));
    chkHDF5err(H5Dclose(dataset));
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      /* write velocity dispersion in R-direction */
    dataset = H5Dcreate(patch, "sigmaR", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _sigmaR));
    chkHDF5err(H5Dclose(dataset));
    /* write velocity dispersion in tangential direction */
    dataset = H5Dcreate(patch, "sigmap", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _sigmap));
    chkHDF5err(H5Dclose(dataset));
    /* write circular velocity */
    dataset = H5Dcreate(patch, "vcirc", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _vcirc));
    chkHDF5err(H5Dclose(dataset));
    /* write kappa */
    dataset = H5Dcreate(patch, "kappa", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _kappa));
    chkHDF5err(H5Dclose(dataset));
    /* write Omega */
    dataset = H5Dcreate(patch, "Omega", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _Omega));
    chkHDF5err(H5Dclose(dataset));
    /* write lambda_critical */
    dataset = H5Dcreate(patch, "lambda", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _lambda));
    chkHDF5err(H5Dclose(dataset));
    /* write Toomre's Q-value */
    dataset = H5Dcreate(patch, "Toomre's Q", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, _toomre));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
      /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write vertical variables */
    //---------------------------------------------------------------------
    dataspace = H5Screate_simple(1, &dims[1], NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    chkHDF5err(H5Pset_chunk(property, 1, &cdims[1]));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
      /* write vertical position */
    dataset = H5Dcreate(patch, "height", hdf5_real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_ver));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
      /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write attribute data */
    //---------------------------------------------------------------------
    /* create the data space for the attribute */
    hsize_t attr_dims = 2;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    //---------------------------------------------------------------------
    /* write # of arrays */
    attribute = H5Acreate(patch, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, dims));
    chkHDF5err(H5Aclose(attribute));
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------
    chkHDF5err(H5Gclose(patch));
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* write attribute data */
    //---------------------------------------------------------------------
    /* create the data space for the attribute */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    //---------------------------------------------------------------------
    /* write scale radius */
    attribute = H5Acreate(group, "Rs", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &disk[ii].cfg->rs));
    chkHDF5err(H5Aclose(attribute));
    /* write scale height */
    attribute = H5Acreate(group, "zd", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &disk[ii].cfg->zd));
    chkHDF5err(H5Aclose(attribute));
    /* write total mass */
    attribute = H5Acreate(group, "Mtot", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &disk[ii].cfg->Mtot));
    chkHDF5err(H5Aclose(attribute));
    /* profile ID */
    attribute = H5Acreate(group, "profile ID", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &disk[ii].cfg->kind));
    chkHDF5err(H5Aclose(attribute));
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------
    chkHDF5err(H5Gclose(group));
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  //-----------------------------------------------------------------------
  /* write attribute data */
  //-----------------------------------------------------------------------
  /* create the data space for the attribute */
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  /* write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &ndisk));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write maximum # of nested level */
  attribute = H5Acreate(target, "maxLev", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &maxLev));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* close the file */
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
#else///USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    FILE *fp;
    char filename[256];
    sprintf(filename, "%s/%s.diskdat.%d.dat", DATAFOLDER, file, ii);
    fp = fopen(filename, "wb");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);    }
    //---------------------------------------------------------------------
    bool success = true;
    //---------------------------------------------------------------------
    const int nhorbin = NDISKBIN_HOR;
    const int nverbin = NDISKBIN_VER;
    success &= (fwrite(&nhorbin, sizeof(int), 1, fp) == 1);
    success &= (fwrite(&nverbin, sizeof(int), 1, fp) == 1);
    success &= (fwrite(& maxLev, sizeof(int), 1, fp) == 1);
    //---------------------------------------------------------------------
    success &= (fwrite(&disk[ii].cfg->rs, sizeof(double), 1, fp) == 1);
    success &= (fwrite(&disk[ii].cfg->zd, sizeof(double), 1, fp) == 1);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    for(int lev = 0; lev < maxLev; lev++){
      //-------------------------------------------------------------------
      if( ndisk > 1 )
#pragma omp parallel for
	for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++)
	  rhoFrac[jj] = (real)((*disk[ii].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * disk[ii].rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)]);
      evaluateDiskProperties(disk, ii, lev, 0, NDISKBIN_HOR - 1, _vcirc, _sigmap, _sigmaR, _kappa, _Omega, _toomre, _lambda);
      //-------------------------------------------------------------------
      for(int jj = 0; jj < nhorbin; jj++){
	tmp_hor[jj] = (real)(disk[ii].hor   [INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] *     length2astro);
	tmp_sig[jj] = (real)(disk[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] *    velocity2astro);
	tmp_Sig[jj] = (real)(disk[ii].Sigma [INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] * col_density2astro);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
	tmp__zd[jj] = (real)(disk[ii].zd    [INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] *      length2astro);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      }/* for(int jj = 0; jj < nhorbin; jj++){ */
      for(int jj = 0; jj < nverbin; jj++)
	tmp_ver[jj] = (real)(disk[ii].ver[INDEX2D(maxLev, NDISKBIN_VER, lev, jj)] * length2astro);
      for(int jj = 0; jj < nhorbin * nverbin; jj++){
	tmp_rho[jj] = (real)((*disk[ii].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * density2astro);
	tmp_Phi[jj] = (real)(  disk[ii].pot [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * senergy2astro);
      }/* for(int jj = 0; jj < nhorbin * nverbin; jj++){ */
      //---------------------------------------------------------------------
      success &= (fwrite(tmp_hor, sizeof(real), nhorbin          , fp) == nhorbin          );
      success &= (fwrite(tmp_ver, sizeof(real),           nverbin, fp) ==           nverbin);
      success &= (fwrite(tmp_rho, sizeof(real), nhorbin * nverbin, fp) == nhorbin * nverbin);
      success &= (fwrite(tmp_Phi, sizeof(real), nhorbin * nverbin, fp) == nhorbin * nverbin);
      success &= (fwrite(tmp_sig, sizeof(real), nhorbin          , fp) == nhorbin          );
      success &= (fwrite(tmp_Sig, sizeof(real), nhorbin          , fp) == nhorbin          );
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      success &= (fwrite(tmp__zd, sizeof(real), nhorbin          , fp) == nhorbin          );
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      //-------------------------------------------------------------------
    }/* for(int lev = 0; lev < maxLev; lev++){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    fclose(fp);
    //---------------------------------------------------------------------
    if( success != true ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  //-----------------------------------------------------------------------
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  free(tmp_hor);  free(tmp_ver);
  free(tmp_sig);  free(tmp_Sig);
  free(_vcirc);  free(_sigmap);  free(_sigmaR);
  free(_kappa);  free(_Omega);  free(_lambda);  free(_toomre);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  free(tmp__zd);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  free(tmp_rho);  free(tmp_Phi);  free(rhoFrac);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
