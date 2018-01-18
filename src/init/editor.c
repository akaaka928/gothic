/**
 * @file editor.c
 *
 * @brief Source code for manipulating output of MAGI (MAny-component Galaxy Initializer)
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2018/01/18 (Thu)
 *
 * Copyright (C) 2018 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

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

#include "../misc/structure.h"
#include "../misc/allocate.h"

#ifndef RUN_WITHOUT_GOTHIC
#include "../misc/tune.h"
#include "../misc/brent.h"
#endif//RUN_WITHOUT_GOTHIC

#include "../file/io.h"


/* global constants to set unit system, defined in constants.c */
extern const double length_astro2com,   length2astro;extern const char   length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double   time_astro2com,     time2astro;extern const char     time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double                   velocity2astro;extern const char velocity_astro_unit_name[CONSTANTS_H_CHAR_WORDS];


/**
 * @struct component
 *
 * @brief structure for each component
 */
typedef struct
{
  ulong  num;/**< number of particles for the specified component */
  ulong head;/**< head index of particles for the specified component */
  int skip;/**< whether to skip the component (1) or not (0) */
  int disk;/**< whether the component is disk (1) or not (0) */
} component;


/**
 * @struct object
 *
 * @brief structure for each object
 */
typedef struct
{
  char file[128];/**< file name of particle distribution */
  char cfg[128];/**< file name of configuration */
  double angle;/**< rotation angle of the system (in radian) */
  double ax, ay, az;/**< rotation axis of the system */
  double xx, yy, zz;/**< position of center-of-mass */
  double vx, vy, vz;/**< bulk velocity of the system */
  int kind;/**< number of components */
  int head;/**< head index of components */
} object;


/**
 * @fn writeSystemCfgFormat
 *
 * @brief Print the expected format.
 *
 * @param (filename) the specified file name
 */
static inline void writeSystemCfgFormat(char *filename)
{
  __NOTE__("%s\n", "start");

  fprintf(stderr, "ERROR: information written in \"%s\" does not match with the expected format\n", filename);
  fprintf(stderr, "Expected format is below:\n");
  fprintf(stderr, "\tangle<real>: rotation angle in units of radian\n");
  fprintf(stderr, "\tax<real> ay<real> az<real>: rotation axis of the system\n");
  fprintf(stderr, "\tx<real> y<real> z<real>: initial position of the system\n");
  fprintf(stderr, "\tvx<real> vy<real> vz<real>: initial velocity of the system\n");
  fprintf(stderr, "\tskip<int>: number of components to be removed\n");
  fprintf(stderr, "\tidx0<int>: index of the 1st component to be removed\n");
  fprintf(stderr, "\tidx1<int>: index of the 2nd component to be removed\n");
  fprintf(stderr, "\tidx2<int>: index of the 3rd component to be removed\n");

  __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);

  __NOTE__("%s\n", "end");
}

/**
 * @fn writeEditorCfgFormat
 *
 * @brief Print the expected format.
 *
 * @param (filename) the specified file name
 */
static inline void writeEditorCfgFormat(char *filename)
{
  __NOTE__("%s\n", "start");

  fprintf(stderr, "ERROR: information written in \"%s\" does not match with the expected format\n", filename);
  fprintf(stderr, "Expected format is below:\n");
  fprintf(stderr, "\tnum<int>: number of input files\n");
  fprintf(stderr, "\tfile0<char *> cfg0<char *>: filename of particle data and configuration, respectively, for the 1st object\n");
  fprintf(stderr, "\tfile1<char *> cfg1<char *>: filename of particle data and configuration, respectively, for the 2nd object\n");
  fprintf(stderr, "\tfile2<char *> cfg2<char *>: filename of particle data and configuration, respectively, for the 3rd object\n");

  __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);

  __NOTE__("%s\n", "end");
}

/**
 * @fn readEditorCfg
 *
 * @brief Read summary of object(s) to be edited.
 *
 * @param (cfg) file name of the configuration
 * @return (unit) unit system of the simulation
 * @return (Nobj) number of objects
 * @return (obj) summary of the object(s)
 * @return (Ncmp) number of components
 * @return (cmp) summary of component(s) in the object(s)
 *
 * @sa writeEditorCfgFormat
 * @sa writeSystemCfgFormat
 */
static inline void readEditorCfg(char *cfg, int *unit, int *Nobj, object **obj, int *Ncmp, component **cmp)
{
  __NOTE__("%s\n", "start");


  /** read global settings */
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s", CFGFOLDER, cfg);

  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  int checker = 1;

  /** read the number of input files */
  checker &= (1 == fscanf(fp, "%d", Nobj));
  *obj = (object *)malloc(sizeof(object) * (*Nobj));  if( *obj == NULL ){    __KILL__(stderr, "ERROR: failure to allocate obj\n");  }

  for(int ii = 0; ii < *Nobj; ii++)
    checker &= (2 == fscanf(fp, "%s %s", (*obj)[ii].file, (*obj)[ii].cfg));

  fclose(fp);
  if( !checker )
    writeEditorCfgFormat(filename);


  int kind = 0;
  for(int ii = 0; ii < *Nobj; ii++){
    sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, (*obj)[ii].file);
    fp = fopen(filename, "r");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }

    int unit_tmp;
    checker = 1;
    checker &= (1 == fscanf(fp, "%d", &unit_tmp));
    checker &= (1 == fscanf(fp, "%d\t%*d", &(*obj)[ii].kind));

    fclose(fp);
    if( !checker ){      __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);    }

    (*obj)[ii].head = kind;
    kind += (*obj)[ii].kind;
    if( ii > 0 ){
      if( unit_tmp != *unit ){
	__KILL__(stderr, "ERROR: unit system of all objects listed in \"%s/%s\" must be identical.\n", CFGFOLDER, cfg);
      }/* if( unit_tmp != unit ){ */
    }/* if( ii > 0 ){ */
    else
      *unit = unit_tmp;
  }/* for(int ii = 0; ii < *Nobj; ii++){ */

  *Ncmp = kind;
  *cmp = (component *)malloc(sizeof(component) * kind);  if( *cmp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate cmp\n");  }


  /** read individual objects */
  for(int ii = 0; ii < *Nobj; ii++){
    sprintf(filename, "%s/%s", CFGFOLDER, (*obj)[ii].cfg);
    fp = fopen(filename, "r");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    checker = 1;

    checker &= (1 == fscanf(fp, "%le", &(*obj)[ii].angle));
    checker &= (3 == fscanf(fp, "%le %le %le", &(*obj)[ii].ax, &(*obj)[ii].ay, &(*obj)[ii].az));
    checker &= (3 == fscanf(fp, "%le %le %le", &(*obj)[ii].xx, &(*obj)[ii].yy, &(*obj)[ii].zz));
    checker &= (3 == fscanf(fp, "%le %le %le", &(*obj)[ii].vx, &(*obj)[ii].vy, &(*obj)[ii].vz));

    int skip = 0;
    checker &= (1 == fscanf(fp, "%d", &skip));

    FILE *sfp;
    char sfile[128];
    sprintf(sfile, "%s/%s.summary.txt", DOCUMENTFOLDER, (*obj)[ii].file);
    sfp = fopen(sfile, "r");
    if( sfp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", sfile);    }
    int success = 1;
    int ndisk, skind;
    success &= (0 == fscanf(sfp, "%*d"));
    success &= (1 == fscanf(sfp, "%*d\t%d", &ndisk));
    skind = (*obj)[ii].kind - ndisk;

    ulong head = 0;
    for(int jj = 0; jj < (*obj)[ii].kind; jj++){
      success &= (1 == fscanf(sfp, "%zu", &(*cmp)[(*obj)[ii].head + jj].num));

      (*cmp)[(*obj)[ii].head + jj].head = head;
      head += (*cmp)[(*obj)[ii].head + jj].num;

      (*cmp)[(*obj)[ii].head + jj].skip = 0;

      if( jj >= skind )
	(*cmp)[(*obj)[ii].head + jj].disk = 1;
    }/* for(int jj = 0; jj < (*obj)[ii].kind; jj++){ */

    fclose(sfp);
    if( !success ){      __KILL__(stderr, "ERROR: failure to read \"%s\"\n", sfile);    }

    for(int jj = 0; jj < skip; jj++){
      int idx;
      checker &= (1 == fscanf(fp, "%d", &idx));

      (*cmp)[(*obj)[ii].head + idx].skip = 1;
    }/* for(int jj = 0; jj < skip; jj++){ */

    fclose(fp);
    if( !checker )
      writeSystemCfgFormat(filename);
  }

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
  if(          (ulong *)a > (ulong *)b ){    return ( 1);  }
  else{    if( (ulong *)a < (ulong *)b ){    return (-1);  }
    else{                                                                    return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

#ifdef __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC


static inline void addSystem(object obj, component *cmp, ulong *head, const iparticle body
#ifdef  USE_HDF5_FORMAT
			     , hdf5struct hdf5type
#endif//USE_HDF5_FORMAT
			     )
{
  __NOTE__("%s\n", "start");


  /**< count up number of N-body particles in the object */
  ulong num = 0;
  for(int ii = obj.head; ii < obj.head + obj.kind; ii++)
    num += cmp[ii].num;

  /** memory allocation for tentative space */
  iparticle tmp;
  ulong *idx;
  position *pos;
  acceleration *acc;
#ifdef  BLOCK_TIME_STEP
  velocity *vel;
  ibody_time *ti;
#else///BLOCK_TIME_STEP
  real *vx, *vy, *vz;
#endif//BLOCK_TIME_STEP
  allocParticleData((int)num, &tmp, &idx, &pos, &acc,
#ifdef  BLOCK_TIME_STEP
		    &vel, &ti
#else///BLOCK_TIME_STEP
		    &vx, &vy, &vz
#endif//BLOCK_TIME_STEP
		    );

  /**< read distribution of N-body particles */
  double time, dt;
  ulong steps;
#ifndef RUN_WITHOUT_GOTHIC
#ifdef  MONITOR_ENERGY_ERROR
  static energyError relEneErr = {1.0, DBL_MIN};
#endif//MONITOR_ENERGY_ERROR
  static int dropPrevTune;
  static rebuildTree rebuild;
  static measuredTime measured;
  static autoTuningParam rebuildParam;
  static brentStatus status;
  static brentMemory memory;
#endif//RUN_WITHOUT_GOTHIC
  readTentativeData(&time, &dt, &steps, (int)num, tmp, obj.file, 0
#ifdef  USE_HDF5_FORMAT
		    , hdf5type
#ifndef RUN_WITHOUT_GOTHIC
		    , &dropPrevTune, &rebuild, &measured, &rebuildParam, &status, &memory
#ifdef  MONITOR_ENERGY_ERROR
		    , &relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
		    );

  /**< set rotation matrix */
  static real axis[3];
  axis[0] = CAST_D2R(obj.ax);
  axis[1] = CAST_D2R(obj.ay);
  axis[2] = CAST_D2R(obj.az);
  real invnorm = RSQRT(CAST_D2R(1.0e-20) + axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  for(int ii = 0; ii < 3; ii++)
    axis[ii] *= invnorm;
  static real rot[3][3], inv[3][3];
  setRodriguesRotationMatrix(axis, CAST_D2R(sin(obj.angle)), CAST_D2R(cos(obj.angle)), rot, inv);

  /**< set shift vector */
  static real shift[3];
  shift[0] = CAST_D2R(obj.xx);
  shift[1] = CAST_D2R(obj.yy);
  shift[2] = CAST_D2R(obj.zz);

  /**< set drift vector */
  static real drift[3];
  drift[0] = CAST_D2R(obj.vx);
  drift[1] = CAST_D2R(obj.vy);
  drift[2] = CAST_D2R(obj.vz);

#pragma omp parallel for
  for(ulong ii = 0; ii < num; ii++){
    /** rotate and shift particle position */
    const real xx = tmp.pos[ii].x;
    const real yy = tmp.pos[ii].y;
    const real zz = tmp.pos[ii].z;
    tmp.pos[ii].x = shift[0] + rot[0][0] * xx + rot[0][1] * yy + rot[0][2] * zz;
    tmp.pos[ii].y = shift[1] + rot[1][0] * xx + rot[1][1] * yy + rot[1][2] * zz;
    tmp.pos[ii].z = shift[2] + rot[2][0] * xx + rot[2][1] * yy + rot[2][2] * zz;

    /** rotate and drift particle velocity */
#ifdef  BLOCK_TIME_STEP
    const real vx = tmp.vel[ii].x;
    const real vy = tmp.vel[ii].y;
    const real vz = tmp.vel[ii].z;
    tmp.vel[ii].x = drift[0] + rot[0][0] * vx + rot[0][1] * vy + rot[0][2] * vz;
    tmp.vel[ii].y = drift[1] + rot[1][0] * vx + rot[1][1] * vy + rot[1][2] * vz;
    tmp.vel[ii].z = drift[2] + rot[2][0] * vx + rot[2][1] * vy + rot[2][2] * vz;
#else///BLOCK_TIME_STEP
    const real vx = tmp.vx[ii];
    const real vy = tmp.vy[ii];
    const real vz = tmp.vz[ii];
    tmp.vx[ii] = drift[0] + rot[0][0] * vx + rot[0][1] * vy + rot[0][2] * vz;
    tmp.vy[ii] = drift[1] + rot[1][0] * vx + rot[1][1] * vy + rot[1][2] * vz;
    tmp.vz[ii] = drift[2] + rot[2][0] * vx + rot[2][1] * vy + rot[2][2] * vz;
#endif//BLOCK_TIME_STEP
  }/* for(ulong ii = 0; ii < num; ii++){ */

  /**< sort N-body particles */
  qsort(tmp.idx, num, sizeof(ulong), idxAscendingOrder);

  /**< copy particle position and velocity */
  for(int ii = obj.head; ii < obj.head + obj.kind; ii++)
    if( cmp[ii].skip != 1 ){
      for(ulong jj = 0; jj < cmp[ii].num; jj++){
	const ulong src = tmp.idx[cmp[ii].head + jj];
	const ulong dst = (*head) + jj;

	body.pos[dst] = tmp.pos[src];
#ifdef  BLOCK_TIME_STEP
	body.vel[dst] = tmp.vel[src];
#else///BLOCK_TIME_STEP
	body.vx[dst] = tmp.vx[src];
	body.vy[dst] = tmp.vy[src];
	body.vz[dst] = tmp.vz[src];
#endif//BLOCK_TIME_STEP
      }/* for(ulong jj = 0; jj < cmp[ii].num; jj++){ */
      *head += cmp[ii].num;
    }/* if( cmp[ii].skip != 1 ){ */

  /** memory deallocation for tentative space */
  freeParticleData(idx, pos, acc,
#ifdef  BLOCK_TIME_STEP
		    vel, ti
#else///BLOCK_TIME_STEP
		    vx, vy, vz
#endif//BLOCK_TIME_STEP
		    );


  __NOTE__("%s\n", "end");
}


int main(int argc, char **argv)
{
  /** read input arguments */
  if( argc < 9 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 9);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -list=<char *>\n");
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -eps=<real> -eta=<real>\n");
    __FPRINTF__(stderr, "          -ft=<real> -snapshotInterval=<real> -saveInterval=<real>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 9 ){ */

  /** read input arguments do not depend on the unit system adopted in the numerical simulation */
  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "file", &file));
  char *fcfg;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "list", &fcfg));
  double tmp;
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv,  "eta", &tmp));
  const real eta = CAST_D2R(tmp);
  double saveInterval;  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "saveInterval", &saveInterval));

  /** read system configuration and set unit system */
  int unit;
  int Nobj, Ncmp;
  object *obj;
  component *cmp;
  readEditorCfg(fcfg, &unit, &Nobj, &obj, &Ncmp, &cmp);
  setPhysicalConstantsAndUnitSystem(unit, 1);

  /** read input arguments depend on the unit system adopted in the numerical simulation */
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "eps", &tmp));
  const real eps = CAST_D2R(tmp * length_astro2com);
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv,  "ft", &tmp));
  const double ft = (tmp * time_astro2com);
  requiredCmdArg(getCmdArgDbl(argc, (const char * const *)argv, "snapshotInterval", &tmp));
  const double snapshotInterval = ldexp(1.0, (int)floor(log2(tmp * time_astro2com)));


  /** count up total number of N-body particles in the resultant system */
  static ulong Ntot = 0;
  for(int ii = 0; ii < Nobj; ii++)
    for(int jj = obj[ii].head; jj < obj[ii].head + obj[ii].kind; jj++)
      if( cmp[jj].skip == 0 )
	Ntot += cmp[jj].num;


  /** memory allocation */
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

  /** read particle data and modify */
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
#endif//USE_HDF5_FORMAT
  ulong head = 0;
  for(int ii = 0; ii < Nobj; ii++)
    addSystem(obj[ii], cmp, &head, body
#ifdef  USE_HDF5_FORMAT
	      , hdf5type
#endif//USE_HDF5_FORMAT
	      );

  /** initialize particle acceleration, time, and index */
  static acceleration acc_zero = {ZERO, ZERO, ZERO, ZERO};
#ifdef  BLOCK_TIME_STEP
  static ibody_time time_zero = {0.0, 0.0};
#endif//BLOCK_TIME_STEP
  for(ulong ii = 0; ii < Ntot; ii++){
    body.acc[ii] = acc_zero;
#ifdef  BLOCK_TIME_STEP
    body.time[ii] = time_zero;
#endif//BLOCK_TIME_STEP
    body.idx[ii] = ii;
  }/* for(ulong ii = 0; ii < Ntot; ii++){ */

  /** write particle data */
  double time = 0.0;
#   if  !defined(WRITE_IN_TIPSY_FORMAT) && !defined(WRITE_IN_GALACTICS_FORMAT)
  double  dt  = 0.0;
  int   last  = 1;
  ulong steps = 0;
  writeSettings(unit, Ntot, eps, eta, ft, snapshotInterval, saveInterval, file);

#ifdef  USE_HDF5_FORMAT
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

#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
#endif//USE_HDF5_FORMAT


  /** output useful information for multi-component analysis */
  FILE *fp;
  char filename[256];
  sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

  fprintf(fp, "%d\n", unit);

  int kind = 0;
  int skind = 0;
  for(int ii = 0; ii < Nobj; ii++)
    for(int jj = obj[ii].head; jj < obj[ii].head + obj[ii].kind; jj++)
      if( cmp[jj].skip != 1 ){
	kind++;
	skind += (cmp[jj].disk != 1);
      }/* if( cmp[jj].skip != 1 ){ */
  fprintf(fp, "%d\t%d\n", kind, skind);

  for(int ii = 0; ii < Nobj; ii++)
    for(int jj = obj[ii].head; jj < obj[ii].head + obj[ii].kind; jj++)
      if( cmp[jj].skip != 1 )
	fprintf(fp, "%zu\n", cmp[jj].num);

  fclose(fp);


  /** output fundamental information of the particle distribution */
  sprintf(filename, "%s/%s.info.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

  /** output global settings */
  char date[64];
  getPresentDateInStrings(date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Fundamental information about the particle distribution\n");
  fprintf(fp, "\tgenerated on %s", date);
  fprintf(fp, "Physical quantities in Computational and Astrophysical units is listed\n");
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Total number of particles Ntot is %zu (= 2^%u)\n", Ntot, ilog2((uint)Ntot));
  fprintf(fp, "Number of components      kind is %d\n", kind);
  fprintf(fp, "Number of original objects     is %d\n", Nobj);
  fprintf(fp, "Number of original components  is %d\n", Ncmp);
  for(int ii = 0; ii < Nobj; ii++)
    fprintf(fp, "\t%d-th object: contains %d components\n", ii, obj[ii].kind);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Length of Plummer softening  is %e (= %e %s)\n", eps, eps * length2astro, length_astro_unit_name);
  fprintf(fp, "Snapshot interval            is %e (= %e %s)\n", snapshotInterval, snapshotInterval * time2astro, time_astro_unit_name);
  fprintf(fp, "Final time of the simulation is %e (= %e %s)\n",               ft,               ft * time2astro, time_astro_unit_name);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "#############################################################################\n\n");

  head = 0;
  for(int ii = 0; ii < Nobj; ii++){
    fprintf(fp, "#############################################################################\n");
    fprintf(fp, "#############################################################################\n");
    fprintf(fp, "%d-th object: %s (%d components)\n", ii, obj[ii].file, obj[ii].kind);
    fprintf(fp, "#############################################################################\n");
    fprintf(fp, "rotation angle is %e radian\n", obj[ii].angle);
    fprintf(fp, "rotation axis is (%e, %e, %e)\n", obj[ii].ax, obj[ii].ay, obj[ii].az);
    fprintf(fp, "initial position is (%e, %e, %e) = (%e %s, %e %s, %e %s)\n", obj[ii].xx, obj[ii].yy, obj[ii].zz, obj[ii].xx * length2astro, length_astro_unit_name, obj[ii].yy * length2astro, length_astro_unit_name, obj[ii].zz * length2astro, length_astro_unit_name);
    fprintf(fp, "initial velocity is (%e, %e, %e) = (%e %s, %e %s, %e %s)\n", obj[ii].vx, obj[ii].vy, obj[ii].vz, obj[ii].vx * velocity2astro, velocity_astro_unit_name, obj[ii].vy * velocity2astro, velocity_astro_unit_name, obj[ii].vz * velocity2astro, velocity_astro_unit_name);
    fprintf(fp, "#############################################################################\n");

    for(int jj = obj[ii].head; jj < obj[ii].head + obj[ii].kind; jj++)
      if( cmp[jj].skip != 1 ){
	fprintf(fp, "%d-th component in original file:\n", jj - obj[ii].head);
	fprintf(fp, "\tnumber of particles is %zu (about 2^%u)\n", cmp[jj].num, ilog2((int)cmp[jj].num));
	fprintf(fp, "\trange of idx for the component is [%zu, %zu]\n", head, head + cmp[jj].num - 1);
	head += cmp[jj].num;
	fprintf(fp, "\tthis is a %s component\n", cmp[jj].disk ? "disk" : "spherical");
      }/* if( cmp[jj].skip != 1 ){ */

    fprintf(fp, "#############################################################################\n\n");
  }/* for(int ii = 0; ii < Nobj; ii++){ */


  /** memory deallocation */
  free(obj);
  free(cmp);
  freeParticleData(idx, pos, acc,
#ifdef  BLOCK_TIME_STEP
		    vel, ti
#else///BLOCK_TIME_STEP
		    vx, vy, vz
#endif//BLOCK_TIME_STEP
		    );


  return (0);
}
