/**
 * @file extract.c
 *
 * @brief Analyzer for radial profiles of multiple components
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2018/07/10 (Tue)
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
#define USE_SZIP_COMPRESSION

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
#include <mpi.h>

#include <unistd.h>
#include <sys/stat.h>

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
#include "mpilib.h"

#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../file/io.h"


extern const double      length2astro;extern const char      length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double        time2astro;extern const char        time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double        mass2astro;extern const char        mass_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double    velocity2astro;extern const char    velocity_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double     density2astro;extern const char     density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double col_density2astro;extern const char col_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];

static int ncrit_base;
static const int Nminimum = 16;

static const int NAllocUnit = 32;



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
  real rad, hor;
} nbody_particle;


static inline void allocProfileArray(int num, real **rad, real **rho, real **enc, real **sig)
{
  __NOTE__("%s\n", "start");

  *rad = (real *)malloc(sizeof(real) * num);
  *rho = (real *)malloc(sizeof(real) * num);
  *enc = (real *)malloc(sizeof(real) * num);
  *sig = (real *)malloc(sizeof(real) * num);
  if( (*rad == NULL) || (*rho == NULL) || (*enc == NULL) || (*sig == NULL) ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }

  __NOTE__("%s\n", "end");
}

static inline void enlargeProfileArray(int num, real **rad, real **rho, real **enc, real **sig)
{
  __NOTE__("%s\n", "start");

  *rad = realloc(*rad, sizeof(real) * num);
  *rho = realloc(*rho, sizeof(real) * num);
  *enc = realloc(*enc, sizeof(real) * num);
  *sig = realloc(*sig, sizeof(real) * num);
  if( (*rad == NULL) || (*rho == NULL) || (*enc == NULL) || (*sig == NULL) ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }

  __NOTE__("%s\n", "end");
}


static inline void allocHorizontalProfileArray(int num, real **pos, real **Sigma, real **height, real **sigR, real **sigp, real **sigz)
{
  __NOTE__("%s\n", "start");

  *pos    = (real *)malloc(sizeof(real) * num);
  *Sigma  = (real *)malloc(sizeof(real) * num);
  *height = (real *)malloc(sizeof(real) * num);
  *sigR   = (real *)malloc(sizeof(real) * num);
  *sigp   = (real *)malloc(sizeof(real) * num);
  *sigz   = (real *)malloc(sizeof(real) * num);
  if( (*pos == NULL) || (*Sigma == NULL) || (*height == NULL) || (*sigR == NULL) || (*sigp == NULL) || (*sigz == NULL) ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }

  __NOTE__("%s\n", "end");
}

static inline void enlargeHorizontalProfileArray(int num, real **pos, real **Sigma, real **height, real **sigR, real **sigp, real **sigz)
{
  __NOTE__("%s\n", "start");

  *pos    = realloc(*   pos, sizeof(real) * num);
  *Sigma  = realloc(* Sigma, sizeof(real) * num);
  *height = realloc(*height, sizeof(real) * num);
  *sigR   = realloc(*  sigR, sizeof(real) * num);
  *sigp   = realloc(*  sigp, sizeof(real) * num);
  *sigz   = realloc(*  sigz, sizeof(real) * num);
  if( (*pos == NULL) || (*Sigma == NULL) || (*height == NULL) || (*sigR == NULL) || (*sigp == NULL) || (*sigz == NULL) ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
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
  if(          ((const nbody_particle *)a)->rad > ((const nbody_particle *)b)->rad ){    return ( 1);  }
  else{    if( ((const nbody_particle *)a)->rad < ((const nbody_particle *)b)->rad ){    return (-1);  }
    else{                                                                                return ( 0);  }  }
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
  if(          ((const nbody_particle *)a)->hor > ((const nbody_particle *)b)->hor ){    return ( 1);  }
  else{    if( ((const nbody_particle *)a)->hor < ((const nbody_particle *)b)->hor ){    return (-1);  }
    else{                                                                                return ( 0);  }  }
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
  if(          (const real *)a > (const real *)b ){    return ( 1);  }
  else{    if( (const real *)a < (const real *)b ){    return (-1);  }
    else{                                              return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

#ifdef __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC


real getCenter(const int num, nbody_particle *body, real com[restrict], real vel[restrict]);
void analyzeRadialProfile(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body_tot, int *num, int *rem, int * restrict prfHead, int * restrict prfNum, real **rad, real **rho, real **enc, real **sig, real com_tot[restrict], real vel_tot[restrict]);
void analyzeHorizontalProfile(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body_tot, int *num, int *rem, int * restrict prfHead, int * restrict prfNum, real **pos, real **Sigma, real **height, real **sigR, real **sigp, real **sigz, real com_tot[restrict], real vel_tot[restrict]);
void generateMassDistributionMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body_tot, const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz, real * restrict rho_map);
void generateSurfaceDensityMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body, const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx, const int nv, const real vmin, const real dv, real * restrict f_xv, real * restrict f_yv, real * restrict f_zv);

#ifdef  USE_HDF5_FORMAT
void writeAnalyzedProfiles
(const double time, const ulong steps, char file[], const uint id, hdf5struct type, const int kind, int * restrict bodyNum,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nz, real * restrict zz, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx,
 const int nv, real * restrict vv, real * restrict f_xv, real * restrict f_yv, real * restrict f_zv,
 const int nx3D, real * restrict rho_xx, const int ny3D, real * restrict rho_yy, const int nz3D, real * restrict rho_zz, real * restrict rho_map,
 int * restrict prfHead, int * restrict prfNum, real * restrict rad, real * restrict rho, real * restrict enc, real * restrict sig,
   int * restrict prfHead_hor, int * restrict prfNum_hor, real * restrict hor, real * restrict Sigma, real * restrict height, real * restrict sigR, real * restrict sigp, real * restrict sigz,
   real * restrict rhalf, real * restrict com, real * restrict vel);

#ifdef  HDF5_FOR_ZINDAIJI
void writeZindaijiFile(const int Ntot, nbody_hdf5 hdf5, const real eps, const int kind, int *bodyHead, int *bodyNum, int *type, const double time, char file[], const uint id);
#endif//HDF5_FOR_ZINDAIJI
#endif//USE_HDF5_FORMAT


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);

  /** initialization */
  if( argc < 21 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 21);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -ncrit=<int>\n");
    __FPRINTF__(stderr, "          -nx=<int> -xmin=<real> -xmax=<real>\n");
    __FPRINTF__(stderr, "          -ny=<int> -ymin=<real> -ymax=<real>\n");
    __FPRINTF__(stderr, "          -nz=<int> -zmin=<real> -zmax=<real>\n");
    __FPRINTF__(stderr, "          -nv=<int> -vmin=<real> -vmax=<real>\n");
    __FPRINTF__(stderr, "          -nx3D=<int> -ny3D=<int> -nz3D=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 21 ){ */

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));
  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "ncrit", &ncrit_base));

  int nx;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nx", &nx));
  int ny;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "ny", &ny));
  int nz;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nz", &nz));
  int nx3D;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nx3D", &nx3D));
  int ny3D;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "ny3D", &ny3D));
  int nz3D;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nz3D", &nz3D));
  real xmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "xmin", &xmin));
  real ymin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "ymin", &ymin));
  real zmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "zmin", &zmin));
  real xmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "xmax", &xmax));
  real ymax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "ymax", &ymax));
  real zmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "zmax", &zmax));

  int nv;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nv", &nv));
  real vmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vmin", &vmin));
  real vmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vmax", &vmax));

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

  setPhysicalConstantsAndUnitSystem(unit, 1);

  /** read number of components */
  int kind = 0;
  int skind = 0;
  int *bodyHead, *bodyNum;
  FILE *fp;
  char filename[128];
  sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "r");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  int unit_tmp;
  bool checker = true;
  checker &= (1 == fscanf(fp, "%d", &unit_tmp));
  checker &= (2 == fscanf(fp, "%d\t%d", &kind, &skind));
  bodyHead = (int *)malloc(sizeof(int) * kind);  if( bodyHead == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate bodyHead");  }
  bodyNum  = (int *)malloc(sizeof(int) * kind);  if( bodyNum  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate bodyNum");  }
#ifdef  HDF5_FOR_ZINDAIJI
  int *bodyType;
  bodyType = (int *)malloc(sizeof(int) * kind);  if( bodyType == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate bodyType");  }
#endif//HDF5_FOR_ZINDAIJI
  for(int ii = 0; ii < kind; ii++)
    checker &= (1 == fscanf(fp, "%d", &bodyNum[ii]));
#ifdef  HDF5_FOR_ZINDAIJI
  for(int ii = 0; ii < kind; ii++)
    checker &= (1 == fscanf(fp, "%d", &bodyType[ii]));
#endif//HDF5_FOR_ZINDAIJI
  fclose(fp);
  if( !checker ){
    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
  }/* if( !checker ){ */
  bodyHead[0] = 0;
  for(int ii = 1; ii < kind; ii++)
    bodyHead[ii] = bodyHead[ii - 1] + bodyNum[ii - 1];
#ifdef  HDF5_FOR_ZINDAIJI
  for(int ii = 0; ii < kind; ii++)
    bodyType[ii] &= 3;
#endif//HDF5_FOR_ZINDAIJI


  const int nfile = (end - start + 1) / interval;
  real *comPos;  comPos = (real *)malloc(sizeof(real) * 3 * nfile * kind);  if( comPos == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate comPos");  }
  real *comVel;  comVel = (real *)malloc(sizeof(real) * 3 * nfile * kind);  if( comVel == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate comVel");  }
  real * rhalf;  rhalf  = (real *)malloc(sizeof(real)     * nfile * kind);  if(  rhalf == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rhalf");  }
  for(int ii = 0; ii < nfile * 3 * kind; ii++)
    comPos[ii] = comVel[ii] = ZERO;
  for(int ii = 0; ii < nfile * kind; ii++)
    rhalf[ii] = ZERO;

  int num = 0;
  int rem = NAllocUnit;
  real *rad, *rho, *enc, *sig;
  allocProfileArray(rem, &rad, &rho, &enc, &sig);
  int *prfHead, *prfNum;
  prfHead = (int *)malloc(sizeof(int) * kind);  if( prfHead == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate prfHead");  }
  prfNum  = (int *)malloc(sizeof(int) * kind);  if( prfNum  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate prfNum");  }

  int num_hor = 0;
  int rem_hor = NAllocUnit;
  real *hor, *Sigma, *height;
  real *sigR, *sigp, *sigz;
  allocHorizontalProfileArray(rem_hor, &hor, &Sigma, &height, &sigR, &sigp, &sigz);
  int *prfHead_hor, *prfNum_hor;
  prfHead_hor = (int *)malloc(sizeof(int) * kind);  if( prfHead_hor == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate prfHead_hor");  }
  prfNum_hor  = (int *)malloc(sizeof(int) * kind);  if( prfNum_hor  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate prfNum_hor");  }


  /** assume zone-centered mapping */
  xmin = CAST_D2R(CAST_R2D(xmin) / length2astro);
  ymin = CAST_D2R(CAST_R2D(ymin) / length2astro);
  zmin = CAST_D2R(CAST_R2D(zmin) / length2astro);
  xmax = CAST_D2R(CAST_R2D(xmax) / length2astro);
  ymax = CAST_D2R(CAST_R2D(ymax) / length2astro);
  zmax = CAST_D2R(CAST_R2D(zmax) / length2astro);
  real *xx;  xx = (real *)malloc(sizeof(real) * (nx + 1));  if( xx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xx");  }
  real *yy;  yy = (real *)malloc(sizeof(real) * (ny + 1));  if( yy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate yy");  }
  real *zz;  zz = (real *)malloc(sizeof(real) * (nz + 1));  if( zz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate zz");  }
  const real dx = (xmax - xmin) / (real)nx;  for(int ii = 0; ii < nx + 1; ii++)    xx[ii] = xmin + dx * (real)ii;
  const real dy = (ymax - ymin) / (real)ny;  for(int jj = 0; jj < ny + 1; jj++)    yy[jj] = ymin + dy * (real)jj;
  const real dz = (zmax - zmin) / (real)nz;  for(int kk = 0; kk < nz + 1; kk++)    zz[kk] = zmin + dz * (real)kk;
  real *Sigma_xy;  Sigma_xy = (real *)malloc(sizeof(real) * kind * nx * ny);  if( Sigma_xy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Sigma_xy");  }
  real *Sigma_yz;  Sigma_yz = (real *)malloc(sizeof(real) * kind * ny * nz);  if( Sigma_yz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Sigma_yz");  }
  real *Sigma_zx;  Sigma_zx = (real *)malloc(sizeof(real) * kind * nz * nx);  if( Sigma_zx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Sigma_zx");  }

  real *rho_xx ;  rho_xx  = (real *)malloc(sizeof(real) * (nx3D + 1));  if( rho_xx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rho_xx");  }
  real *rho_yy ;  rho_yy  = (real *)malloc(sizeof(real) * (ny3D + 1));  if( rho_yy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rho_yy");  }
  real *rho_zz ;  rho_zz  = (real *)malloc(sizeof(real) * (nz3D + 1));  if( rho_zz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rho_zz");  }
  const real rho_dx = (xmax - xmin) / (real)nx3D;  for(int ii = 0; ii < nx3D + 1; ii++)    rho_xx[ii] = xmin + rho_dx * (real)ii;
  const real rho_dy = (ymax - ymin) / (real)ny3D;  for(int jj = 0; jj < ny3D + 1; jj++)    rho_yy[jj] = ymin + rho_dy * (real)jj;
  const real rho_dz = (zmax - zmin) / (real)nz3D;  for(int kk = 0; kk < nz3D + 1; kk++)    rho_zz[kk] = zmin + rho_dz * (real)kk;
  real *rho_map;
  rho_map = (real *)malloc(sizeof(real) * kind * nx3D * ny3D * nz3D);  if( rho_map == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rho_map");  }

  vmin = CAST_D2R(CAST_R2D(vmin) / velocity2astro);
  vmax = CAST_D2R(CAST_R2D(vmax) / velocity2astro);
  real *vv;  vv = (real *)malloc(sizeof(real) * (nv + 1));  if( vv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vv");  }
  const real dv = (vmax - vmin) / (real)nv;  for(int ii = 0; ii < nv + 1; ii++)    vv[ii] = vmin + dv * (real)ii;
  real *f_xv;  f_xv = (real *)malloc(sizeof(real) * kind * nx * nv);  if( f_xv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_xv");  }
  real *f_yv;  f_yv = (real *)malloc(sizeof(real) * kind * ny * nv);  if( f_yv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_yv");  }
  real *f_zv;  f_zv = (real *)malloc(sizeof(real) * kind * nz * nv);  if( f_zv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_zv");  }



#ifdef  USE_HDF5_FORMAT
  /** obtain minimum and maximum values for matplotlib */
  double radmin = DBL_MAX, radmax = 0.0;
  double rhomin = DBL_MAX, rhomax = 0.0;
  double encmin = DBL_MAX, encmax = 0.0;
  double                   sigmax = 0.0;
  double hormin = DBL_MAX, hormax = 0.0;
  double Sigmin = DBL_MAX, Sigmax = 0.0;
  double sigRmax = 0.0, sigpmax = 0.0, sigzmax = 0.0;
  double diskmax = 0.0;
#endif//USE_HDF5_FORMAT


  for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){
    /** read snapshot */
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


    /** sort particle data by index */
    qsort(body, Ntot, sizeof(nbody_particle), idxAscendingOrder);


#ifdef  USE_HDF5_FORMAT
#ifdef  HDF5_FOR_ZINDAIJI
    sprintf(filename, "%s/%s_%s%.3u.h5", DATAFOLDER, file, "zindaiji", filenum);
#else///HDF5_FOR_ZINDAIJI
    sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "split", filenum);
#endif//HDF5_FOR_ZINDAIJI
    bool dump_file = (0 != access(filename, F_OK));
    if( !dump_file ){
      struct stat stat_file;
      stat(filename, &stat_file);
      char tmpname[128];
      sprintf(tmpname, "%s/%s.%s%.3u.h5", DATAFOLDER, file, SNAPSHOT, filenum);
      struct stat stat_snap;
      stat(tmpname, &stat_snap);
      if( stat_snap.st_ctime > stat_file.st_ctime )
	dump_file = true;
    }/* if( !dump_file ){ */
    if( dump_file ){
      for(int ii = 0; ii < (int)Ntot; ii++){
	hdf5.idx[ii    ] = body[ii].idx;	hdf5.  m[ii        ] = body[ii]. m;	hdf5.pot[ii        ] = body[ii].pot;
	hdf5.pos[ii * 3] = body[ii]. x ;	hdf5.pos[ii * 3 + 1] = body[ii]. y;	hdf5.pos[ii * 3 + 2] = body[ii]. z ;
	hdf5.vel[ii * 3] = body[ii].vx ;	hdf5.vel[ii * 3 + 1] = body[ii].vy;	hdf5.vel[ii * 3 + 2] = body[ii].vz ;
	hdf5.acc[ii * 3] = body[ii].ax ;	hdf5.acc[ii * 3 + 1] = body[ii].ay;	hdf5.acc[ii * 3 + 2] = body[ii].az ;
#   if  defined(SET_EXTERNAL_POTENTIAL_FIELD) && !defined(HDF5_FOR_ZINDAIJI)
	hdf5.acc_ext[ii * 3] = body[ii].ax_ext;	  hdf5.acc_ext[ii * 3 + 1] = body[ii].ay_ext;	  hdf5.acc_ext[ii * 3 + 2] = body[ii].az_ext;
	hdf5.pot_ext[ii] = body[ii].pot_ext;
#endif//defined(SET_EXTERNAL_POTENTIAL_FIELD) && !defined(HDF5_FOR_ZINDAIJI)
      }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#ifdef  HDF5_FOR_ZINDAIJI
      writeZindaijiFile(Ntot, hdf5, eps, kind, bodyHead, bodyNum, bodyType, time, file, (uint)filenum);
#else///HDF5_FOR_ZINDAIJI
      writeSnapshotMultiGroups(time, steps, &hdf5, file, filenum, hdf5type, kind, bodyHead, bodyNum);
#endif//HDF5_FOR_ZINDAIJI
    }/* if( dump_file ){ */
    else{
      fprintf(stdout, "# \"%s\" was not updated for reducing the elapsed time.\n", filename);
      fflush(stdout);
    }/* else{ */
#endif//USE_HDF5_FORMAT


    /** obtain the location of the center-of-mass and its bulk velocity for particles within half-mass radius */
    for(int kk = 0; kk < kind; kk++)
      rhalf[INDEX2D(nfile, kind, (filenum - start) / interval, kk)] = getCenter(bodyNum[kk], &body[bodyHead[kk]], &comPos[INDEX(nfile, kind, 3, (filenum - start) / interval, kk, 0)], &comVel[INDEX(nfile, kind, 3, (filenum - start) / interval, kk, 0)]);


    /** obtain radial profile */
    analyzeRadialProfile(kind, bodyHead, bodyNum, body, &num, &rem, prfHead, prfNum, &rad, &rho, &enc, &sig, &comPos[INDEX(nfile, kind, 3, (filenum - start) / interval, 0, 0)], &comVel[INDEX(nfile, kind, 3, (filenum - start) / interval, 0, 0)]);
    rem += num;
    num = 0;

    /** obtain horizontal profile */
    analyzeHorizontalProfile(kind, bodyHead, bodyNum, body, &num_hor, &rem_hor, prfHead_hor, prfNum_hor, &hor, &Sigma, &height, &sigR, &sigp, &sigz, &comPos[INDEX(nfile, kind, 3, (filenum - start) / interval, 0, 0)], &comVel[INDEX(nfile, kind, 3, (filenum - start) / interval, 0, 0)]);
    rem_hor += num_hor;
    num_hor = 0;


    /** obtain mass distribution map for volume rendering */
    generateMassDistributionMaps(kind, bodyHead, bodyNum, body, eps, nx3D, xmin, rho_dx, ny3D, ymin, rho_dy, nz3D, zmin, rho_dz, rho_map);
    generateSurfaceDensityMaps(kind, bodyHead, bodyNum, body, eps, nx, xmin, dx, ny, ymin, dy, nz, zmin, dz, Sigma_xy, Sigma_yz, Sigma_zx, nv, vmin, dv, f_xv, f_yv, f_zv);


#ifdef  USE_HDF5_FORMAT
    /** dump analyzed results for matplotlib and/or VisIt */
    writeAnalyzedProfiles
      (time, steps, file, filenum, hdf5type, kind, bodyNum,
       nx, xx, ny, yy, nz, zz, Sigma_xy, Sigma_yz, Sigma_zx,
       nv, vv, f_xv, f_yv, f_zv,
       nx3D, rho_xx, ny3D, rho_yy, nz3D, rho_zz, rho_map,
       prfHead, prfNum, rad, rho, enc, sig,
       prfHead_hor, prfNum_hor, hor, Sigma, height, sigR, sigp, sigz,
       &rhalf[INDEX2D(nfile, kind, (filenum - start) / interval, 0)], &comPos[INDEX(nfile, kind, 3, (filenum - start) / interval, 0, 0)], &comVel[INDEX(nfile, kind, 3, (filenum - start) / interval, 0, 0)]);


    /** obtain minimum and maximum values for matplotlib */
    for(int kk = 0; kk < kind; kk++)
      if( bodyNum[kk] > Nminimum ){
	radmin = fmin(radmin, rad[prfHead[kk]]);	radmax = fmax(radmax, rad[prfHead[kk] + prfNum[kk] - 1]);
	encmin = fmin(encmin, enc[prfHead[kk]]);	encmax = fmax(encmax, enc[prfHead[kk] + prfNum[kk] - 1]);

	hormin = fmin(hormin, hor[prfHead_hor[kk]                     ]);
	hormax = fmax(hormax, hor[prfHead_hor[kk] + prfNum_hor[kk] - 1]);
      }/* if( bodyNum[kk] > Nminimum ){ */

    for(int ii = prfHead[0]; ii < prfHead[kind - 1] + prfNum[kind - 1]; ii++){
      rhomin = fmin(rhomin, rho[ii]);
      rhomax = fmax(rhomax, rho[ii]);
      sigmax = fmax(sigmax, sig[ii]);
    }/* for(int ii = prfHead[0]; ii < prfHead[kind - 1] + prfNum[kind - 1]; ii++){ */

    for(int ii = prfHead_hor[0]; ii < prfHead_hor[kind - 1] + prfNum_hor[kind - 1]; ii++){
      Sigmin = fmin(Sigmin, Sigma[ii]);
      Sigmax = fmax(Sigmax, Sigma[ii]);
    }/* for(int ii = prfHead_hor[0]; ii < prfHead_hor[kind - 1] + prfNum_hor[kind - 1]; ii++){ */

    for(int ii = prfHead_hor[skind]; ii < prfHead_hor[kind - 1] + prfNum_hor[kind - 1]; ii++){
      sigRmax = fmax(sigRmax, sigR[ii]);
      sigpmax = fmax(sigpmax, sigp[ii]);
      sigzmax = fmax(sigzmax, sigz[ii]);

      diskmax = fmax(diskmax, height[ii]);
    }/* for(int ii = prfHead_hor[0]; ii < prfHead_hor[kind - 1] + prfNum_hor[kind - 1]; ii++){ */
#endif//USE_HDF5_FORMAT
  }/* for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){ */


  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? comPos : MPI_IN_PLACE, comPos, nfile * kind * 3, MPI_REALDAT, MPI_SUM, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? comVel : MPI_IN_PLACE, comVel, nfile * kind * 3, MPI_REALDAT, MPI_SUM, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ?  rhalf : MPI_IN_PLACE,  rhalf, nfile * kind    , MPI_REALDAT, MPI_SUM, 0, mpi.comm));
#ifdef  USE_HDF5_FORMAT
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & radmin: MPI_IN_PLACE, & radmin, 1, MPI_DOUBLE, MPI_MIN, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & radmax: MPI_IN_PLACE, & radmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & hormin: MPI_IN_PLACE, & hormin, 1, MPI_DOUBLE, MPI_MIN, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & hormax: MPI_IN_PLACE, & hormax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & encmin: MPI_IN_PLACE, & encmin, 1, MPI_DOUBLE, MPI_MIN, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & encmax: MPI_IN_PLACE, & encmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & rhomin: MPI_IN_PLACE, & rhomin, 1, MPI_DOUBLE, MPI_MIN, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & rhomax: MPI_IN_PLACE, & rhomax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & Sigmin: MPI_IN_PLACE, & Sigmin, 1, MPI_DOUBLE, MPI_MIN, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & Sigmax: MPI_IN_PLACE, & Sigmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & sigmax: MPI_IN_PLACE, & sigmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? &sigRmax: MPI_IN_PLACE, &sigRmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? &sigpmax: MPI_IN_PLACE, &sigpmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? &sigzmax: MPI_IN_PLACE, &sigzmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? &diskmax: MPI_IN_PLACE, &diskmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
#endif//USE_HDF5_FORMAT

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

    fprintf(stdout, "# half-mass radius\n");
    fprintf(stdout, "# file");
    for(int ii = 0; ii < kind; ii++)      fprintf(stdout, "\trhalf(%d)", ii);
    fprintf(stdout, "\n");
    for(int ii = 0; ii < nfile; ii++){
      fprintf(stdout, "%d", start + ii * interval);
      for(int kk = 0; kk < kind; kk++)
	fprintf(stdout, "\t%e", rhalf[INDEX2D(nfile, kind, ii, kk)]);
      fprintf(stdout, "\n");
    }/* for(int ii = 0; ii < nfile; ii++){ */

    sprintf(filename, "%s/%s.minmax.txt", DATAFOLDER, file);
    fp = fopen(filename, "w");
    if( fp == NULL ){
      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
    }
    fprintf(fp, "%e\t%e\n", radmin, radmax);
    fprintf(fp, "%e\t%e\n", rhomin, rhomax);
    fprintf(fp, "%e\t%e\n", encmin, encmax);
    fprintf(fp, "%e\t%e\n",    0.0, sigmax);
    fprintf(fp, "%e\t%e\n", hormin, hormax);
    fprintf(fp, "%e\t%e\n", Sigmin, Sigmax);
    fprintf(fp, "%e\t%e\n",    0.0, sigRmax);
    fprintf(fp, "%e\t%e\n",    0.0, sigpmax);
    fprintf(fp, "%e\t%e\n",    0.0, sigzmax);
    fprintf(fp, "%e\t%e\n",    0.0, diskmax);
    fclose(fp);
  }/* if( mpi.rank == 0 ){ */


  free(comPos);
  free(comVel);
  free(rhalf);


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
  free(bodyHead);
  free(bodyNum);
#ifdef  HDF5_FOR_ZINDAIJI
  free(bodyType);
#endif//HDF5_FOR_ZINDAIJI

  free(prfHead);  free(prfNum);
  free(rad);  free(rho);  free(enc);  free(sig);

  free(prfHead_hor);  free(prfNum_hor);
  free(hor);  free(Sigma);  free(height);  free(sigR);  free(sigp);  free(sigz);

  free(xx);  free(yy);  free(zz);
  free(Sigma_xy);  free(Sigma_yz);  free(Sigma_zx);

  free(rho_xx);  free(rho_yy);  free(rho_zz);
  free(rho_map);

  free(vv);  free(f_xv);  free(f_yv);  free(f_zv);

  exitMPI();

  return (0);
}


real getCenter(const int num, nbody_particle *body, real com[restrict], real vel[restrict])
{
  __NOTE__("%s for %d particles\n", "start", num);


  real rhalf = ZERO;

  if( num > 1 ){
    double comx = 0.0;
    double comy = 0.0;
    double comz = 0.0;
    double velx = 0.0;
    double vely = 0.0;
    double velz = 0.0;

    int Npart = num >> 1;

    bool converge = false;
    int steps = 0;
    while( true ){
      for(int ii = 0; ii < num; ii++){
	const real xx = CAST_D2R(CAST_R2D(body[ii].x) - comx);
	const real yy = CAST_D2R(CAST_R2D(body[ii].y) - comy);
	const real zz = CAST_D2R(CAST_R2D(body[ii].z) - comz);
	/* is recipe for correcting precession required?? */
	const real R2 = 1.0e-30f + xx * xx + yy * yy;
	const real r2 = R2 + zz * zz;
	body[ii].hor = R2 * RSQRT(R2);
	body[ii].rad = r2 * RSQRT(r2);
      }/* for(int ii = 0; ii < num; ii++){ */

      /** sort by particle position */
      qsort(body, num, sizeof(nbody_particle), radAscendingOrder);
      rhalf = (num & 1) ? (body[num >> 1].rad) : (HALF * (body[num >> 1].rad + body[(num >> 1) + 1].rad));

      if( converge || (Npart < 1) )
	break;

      velx = 0.0;
      vely = 0.0;
      velz = 0.0;

      double mtot = 0.0;
      double newx = 0.0;
      double newy = 0.0;
      double newz = 0.0;

      /** find center-of-mass and bulk velocity of particles within the half-mass radius (assumption: equal-mass particles) */
      for(int ii = 0; ii < Npart; ii++){
	const double mass = CAST_R2D(body[ii].m);
	mtot += mass;

	const double xx = CAST_R2D(body[ii].x);
	const double yy = CAST_R2D(body[ii].y);
	const double zz = CAST_R2D(body[ii].z);

	const double vx = CAST_R2D(body[ii].vx);
	const double vy = CAST_R2D(body[ii].vy);
	const double vz = CAST_R2D(body[ii].vz);

	newx += mass * xx;
	newy += mass * yy;
	newz += mass * zz;

	velx += mass * vx;
	vely += mass * vy;
	velz += mass * vz;
      }/* for(int ii = 0; ii < Npart; ii++){ */

      const double minv = 1.0 / mtot;
      newx *= minv;
      newy *= minv;
      newz *= minv;
      velx *= minv;
      vely *= minv;
      velz *= minv;

      const double dx = newx - comx;
      const double dy = newy - comy;
      const double dz = newz - comz;

      converge = (dx * dx + dy * dy + dz * dz) < 1.0e-6;

      if( (steps > 1024) && !converge ){
	__FPRINTF__(stderr, "Warning: does not converged after %d iterations (%d particles; error are %e): shrink Npart(%d) to %d.\n", steps, num, dx * dx + dy * dy + dz * dz, Npart, Npart >> 1);
	Npart >>= 1;
	steps = -1;
      }/* if( (steps > 1024) && !converge ){ */

      comx = newx;
      comy = newy;
      comz = newz;

      steps++;
    }/* while( true ){ */
    com[0] = CAST_D2R(comx);
    com[1] = CAST_D2R(comy);
    com[2] = CAST_D2R(comz);
    vel[0] = CAST_D2R(velx);
    vel[1] = CAST_D2R(vely);
    vel[2] = CAST_D2R(velz);
  }/* if( num > 1 ){ */
  else{
    com[0] = body[0].x;
    com[1] = body[0].y;
    com[2] = body[0].z;
    vel[0] = body[0].vx;
    vel[1] = body[0].vy;
    vel[2] = body[0].vz;
  }/* else{ */

  __NOTE__("%s\n", "end");
  return (rhalf);
}


void analyzeRadialProfile(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body_tot, int *num, int *rem, int * restrict prfHead, int * restrict prfNum, real **rad, real **rho, real **enc, real **sig, real com_tot[restrict], real vel_tot[restrict])
{
  __NOTE__("%s\n", "start");


  /** make radial profile */
  *num = 0;
  for(int kk = 0; kk < kind; kk++){
    nbody_particle *body;
    body = &body_tot[bodyHead[kk]];

    const real com[3] = {com_tot[INDEX2D(kind, 3, kk, 0)], com_tot[INDEX2D(kind, 3, kk, 1)], com_tot[INDEX2D(kind, 3, kk, 2)]};
    const real vel[3] = {vel_tot[INDEX2D(kind, 3, kk, 0)], vel_tot[INDEX2D(kind, 3, kk, 1)], vel_tot[INDEX2D(kind, 3, kk, 2)]};

    prfHead[kk] = *num;
    if( bodyNum[kk] > Nminimum ){
      /* /\* sort the array in ascending distance from the center *\/ */
      /* qsort(body, bodyNum[kk], sizeof(nbody_particle), radAscendingOrder); */

      real *prad, *prho, *penc, *psig;
      real inner = ZERO;
      real Menc  = ZERO;
      prad = *rad;
      prho = *rho;
      penc = *enc;
      psig = *sig;
      const int ncrit = ((int)(bodyNum[kk] / ncrit_base) > Nminimum) ? ncrit_base : (int)(bodyNum[kk] / Nminimum);
      const real inv_ncrit = UNITY / (real)ncrit;

      for(int head = 0; head < bodyNum[kk]; head += ncrit){
	/** check # of unused elements */
	if( *rem == 0 ){
	  enlargeProfileArray(*num + NAllocUnit, rad, rho, enc, sig);
	  *rem += NAllocUnit;
	  prad = *rad;
	  prho = *rho;
	  penc = *enc;
	  psig = *sig;
	}/* if( *rem == 0 ){ */

	int tail = head + ncrit;
	if( tail - 1 >= bodyNum[kk] )
	  break;

	real outer = body[tail - 1].rad;
	prad[*num] = (ncrit & 1) ? (body[head + (ncrit >> 1)].rad) : (HALF * (body[head + (ncrit >> 1) - 1].rad + body[head + (ncrit >> 1)].rad));/**< use median of particle location */

	real vr2 = ZERO;
	real vrm = ZERO;
	for(int ii = head; ii < tail; ii++){
	  const real vr = ((body[ii].vx - vel[0]) * (body[ii].x - com[0]) + (body[ii].vy - vel[1]) * (body[ii].y - com[1]) + (body[ii].vz - vel[2]) * (body[ii].z - com[2])) / body[ii].rad;
	  vrm += vr;
	  vr2 += vr * vr;
	}/* for(int ii = head; ii < tail; ii++){ */
	vrm *= inv_ncrit;
	vr2 *= inv_ncrit;
	psig[*num] = SQRT(vr2 - vrm * vrm);

	real mass = ZERO;
	for(int ii = head; ii < tail; ii++)
	  mass += body[ii].m;
	prho[*num] = mass / (CAST_D2R(4.0 * M_PI / 3.0) * (outer * outer * outer - inner * inner * inner));
	penc[*num] = Menc + HALF * mass;

	Menc += mass;
	inner = outer;

	*num += 1;
	*rem -= 1;
      }/* for(int head = 0; head < bodyNum[kk]; head += ncrit){ */
    }/* if( bodyNum[kk] > Nminimum ){ */

    prfNum[kk] = *num - prfHead[kk];
  }/* for(int kk = 0; kk < kind; kk++){ */


  __NOTE__("%s\n", "end");
}


void analyzeHorizontalProfile(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body_tot, int *num, int *rem, int * restrict prfHead, int * restrict prfNum, real **pos, real **Sigma, real **height, real **sigR, real **sigp, real **sigz, real com_tot[restrict], real vel_tot[restrict])
{
  __NOTE__("%s\n", "start");


  /** make horizontal profile on the midplane */
  *num = 0;
  for(int kk = 0; kk < kind; kk++){
    nbody_particle *body;
    body = &body_tot[bodyHead[kk]];

    const real com[3] = {com_tot[INDEX2D(kind, 3, kk, 0)], com_tot[INDEX2D(kind, 3, kk, 1)], com_tot[INDEX2D(kind, 3, kk, 2)]};
    const real vel[3] = {vel_tot[INDEX2D(kind, 3, kk, 0)], vel_tot[INDEX2D(kind, 3, kk, 1)], vel_tot[INDEX2D(kind, 3, kk, 2)]};

    prfHead[kk] = *num;
    if( bodyNum[kk] > Nminimum ){
      /* sort the array in ascending distance from the center */
      qsort(body, bodyNum[kk], sizeof(nbody_particle), horAscendingOrder);

      real *ppos, *pSigma, *pheight, *psigR, *psigp, *psigz;
      real inner = ZERO;
      ppos    = *pos;
      pSigma  = *Sigma;
      pheight = *height;
      psigR   = *sigR;
      psigp   = *sigp;
      psigz   = *sigz;

      const int ncrit = ((int)(bodyNum[kk] / ncrit_base) > Nminimum) ? ncrit_base : (int)(bodyNum[kk] / Nminimum);
      const real inv_ncrit = UNITY / (real)ncrit;
      real *ver;
      ver = (real *)malloc(sizeof(real) * ncrit);
      if( ver == NULL ){
	__KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
      }/* if( ver == NULL ){ */

      for(int head = 0; head < bodyNum[kk]; head += ncrit){
	/* check # of unused elements */
	if( *rem == 0 ){
	  enlargeHorizontalProfileArray(*num + NAllocUnit, pos, Sigma, height, sigR, sigp, sigz);
	  *rem += NAllocUnit;
	  ppos    = *pos;
	  pSigma  = *Sigma;
	  pheight = *height;
	  psigR   = *sigR;
	  psigp   = *sigp;
	  psigz   = *sigz;
	}/* if( *rem == 0 ){ */

	int tail = head + ncrit;
	if( tail - 1 >= bodyNum[kk] )
	  break;

	real outer = body[tail - 1].hor;
	ppos[*num] = (ncrit & 1) ? (body[head + (ncrit >> 1)].hor) : (HALF * (body[head + (ncrit >> 1) - 1].hor + body[head + (ncrit >> 1)].hor));/**< use median of particle location */

	real mass = ZERO;
	real vR2 = ZERO;      real vRm = ZERO;
	real vp2 = ZERO;      real vpm = ZERO;
	real vz2 = ZERO;      real vzm = ZERO;
	for(int ii = head; ii < tail; ii++){
	  mass += body[ii].m;

	  ver[ii - head] = FABS(body[ii].z - com[2]);

	  const real invR = UNITY / body[ii].hor;
	  const real xx = body[ii].x - com[0];
	  const real yy = body[ii].y - com[1];
	  const real vx = body[ii].vx - vel[0];
	  const real vy = body[ii].vy - vel[1];
	  const real vz = body[ii].vz - vel[2];
	  vzm += vz;
	  vz2 += vz * vz;
	  const real vR = ( xx * vx + yy * vy) * invR;
	  const real vp = (-yy * vx + xx * vy) * invR;
	  vRm += vR;
	  vR2 += vR * vR;
	  vpm += vp;
	  vp2 += vp * vp;
	}/* for(int ii = hea; ii < tail; ii++){ */

	pSigma[*num] = mass / (CAST_D2R(M_PI) * (outer * outer - inner * inner));

	qsort(ver, ncrit, sizeof(real), verAscendingOrder);
	pheight[*num] = (ncrit & 1) ? (ver[(ncrit >> 1)]) : (HALF * (ver[(ncrit >> 1) - 1] + ver[(ncrit >> 1)]));/**< use median of height */

	vzm *= inv_ncrit;      vz2 *= inv_ncrit;
	vRm *= inv_ncrit;      vR2 *= inv_ncrit;
	vpm *= inv_ncrit;      vp2 *= inv_ncrit;
	psigz[*num] = SQRT(vz2 - vzm * vzm);
	psigR[*num] = SQRT(vR2 - vRm * vRm);
	psigp[*num] = SQRT(vp2 - vpm * vpm);

	inner = outer;
	*num += 1;
	*rem -= 1;
      }/* for(int head = 0; head < (int)group[kk].num; head += ncrit){ */
      free(ver);
    }/* if( bodyNum[kk] > Nminimum ){ */

    prfNum[kk] = *num - prfHead[kk];
  }/* for(int kk = 0; kk < kind; kk++){ */


  __NOTE__("%s\n", "end");
}


/* 3 sigma means that neglecting component of 0.26% */
/* 5 sigma means that neglecting component of 6e-7 */
/* #define SPREAD (3) */
#define SPREAD (5)
void generateMassDistributionMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body_tot, const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz, real * restrict rho_map)
{
  __NOTE__("%s\n", "start");

  const real dxinv = UNITY / dx;
  const real dyinv = UNITY / dy;
  const real dzinv = UNITY / dz;

  const real sig = HALF * eps;/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real invsig = UNITY / sig;

  const int nx_smooth = (int)CEIL(SPREAD * dxinv * sig);
  const int ny_smooth = (int)CEIL(SPREAD * dyinv * sig);
  const int nz_smooth = (int)CEIL(SPREAD * dzinv * sig);
  real *erfx;  erfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 2));  if( erfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfx");  }
  real *erfy;  erfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 2));  if( erfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfy");  }
  real *erfz;  erfz = (real *)malloc(sizeof(real) * (2 * nz_smooth + 2));  if( erfz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfz");  }
  for(int ii = 0; ii < 2 * nx_smooth + 2; ii++)    erfx[ii] = ERF((dx * ((real)(ii - nx_smooth) - HALF)) * invsig);
  for(int ii = 0; ii < 2 * ny_smooth + 2; ii++)    erfy[ii] = ERF((dy * ((real)(ii - ny_smooth) - HALF)) * invsig);
  for(int ii = 0; ii < 2 * nz_smooth + 2; ii++)    erfz[ii] = ERF((dz * ((real)(ii - nz_smooth) - HALF)) * invsig);

  real *psfx;  psfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 1));  if( psfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfx");  }
  real *psfy;  psfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 1));  if( psfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfy");  }
  real *psfz;  psfz = (real *)malloc(sizeof(real) * (2 * nz_smooth + 1));  if( psfz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfz");  }
  for(int ii = 0; ii < 2 * nx_smooth + 1; ii++)    psfx[ii] = HALF * (erfx[ii + 1] - erfx[ii]);
  for(int ii = 0; ii < 2 * ny_smooth + 1; ii++)    psfy[ii] = HALF * (erfy[ii + 1] - erfy[ii]);
  for(int ii = 0; ii < 2 * nz_smooth + 1; ii++)    psfz[ii] = HALF * (erfz[ii + 1] - erfz[ii]);


  const real dVinv = dxinv * dyinv * dzinv;

  for(int kk = 0; kk < kind; kk++){
    for(int ii = 0; ii < nx * ny * nz; ii++)
      rho_map[INDEX2D(kind, nx * ny * nz, kk, ii)] = ZERO;

    nbody_particle *body;
    body = &body_tot[bodyHead[kk]];

    for(int ii = 0; ii < bodyNum[kk]; ii++){
      const real mi = body[ii].m;
      const real xi = body[ii].x;
      const real yi = body[ii].y;
      const real zi = body[ii].z;

      const int l0 = (int)FLOOR((xi - xmin) * dxinv);
      const int m0 = (int)FLOOR((yi - ymin) * dyinv);
      const int n0 = (int)FLOOR((zi - zmin) * dzinv);

      for(int sx = 0; sx < 2 * nx_smooth + 1; sx++){
	const int ll = l0 + sx - nx_smooth;
	if( (ll >= 0) && (ll < nx) ){
	  const real mx = mi * psfx[sx];

	  for(int sy = 0; sy < 2 * ny_smooth + 1; sy++){
	    const int mm = m0 + sy - ny_smooth;
	    if( (mm >= 0) && (mm < ny) ){
	      const real my = mx * psfy[sy];

	      for(int sz = 0; sz < 2 * nz_smooth + 1; sz++){
		const int nn = n0 + sz - nz_smooth;
		if( (nn >= 0) && (nn < nz) )
		  rho_map[INDEX4D(kind, nx, ny, nz, kk, ll, mm, nn)] += my * psfz[sz];
	      }/* for(int sz = 0; sz < 2 * nz_smooth + 1; sz++){ */

	    }/* if( (mm >= 0) && (mm < ny) ){ */
	  }/* for(int sy = 0; sy < 2 * ny_smooth + 1; sy++){ */

	}/* if( (ll >= 0) && (ll < nx) ){ */
      }/* for(int sx = 0; sx < 2 * nx_smooth + 1; sx++){ */

    }/* for(int ii = 0; ii < bodyNum[kk]; ii++){ */

    for(int ii = 0; ii < nx * ny * nz; ii++)
      rho_map[INDEX2D(kind, nx * ny * nz, kk, ii)] *= dVinv;
  }/* for(int cmp = 0; cmp < kind; cmp++){ */

  free(erfx);  free(erfy);  free(erfz);
  free(psfx);  free(psfy);  free(psfz);


  __NOTE__("%s\n", "end");
}


#define SMOOTHING_FOR_VISUALIZATION TWO
/* #define SMOOTHING_FOR_VISUALIZATION THREE */
void generateSurfaceDensityMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body, const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx, const int nv, const real vmin, const real dv, real * restrict f_xv, real * restrict f_yv, real * restrict f_zv)
{
  __NOTE__("%s\n", "start");

  const real dxinv = UNITY / dx;
  const real dyinv = UNITY / dy;
  const real dzinv = UNITY / dz;
  const real dvinv = UNITY / dv;

  const real xsig = FMAX(HALF * CAST_D2R(CAST_R2D(eps) / length2astro), SMOOTHING_FOR_VISUALIZATION * FABS(dx));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real ysig = FMAX(HALF * CAST_D2R(CAST_R2D(eps) / length2astro), SMOOTHING_FOR_VISUALIZATION * FABS(dy));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real zsig = FMAX(HALF * CAST_D2R(CAST_R2D(eps) / length2astro), SMOOTHING_FOR_VISUALIZATION * FABS(dz));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real vsig =                                                     SMOOTHING_FOR_VISUALIZATION * FABS(dv) ;

  const real invxsig = UNITY / xsig;
  const real invysig = UNITY / ysig;
  const real invzsig = UNITY / zsig;
  const real invvsig = UNITY / vsig;

  const int nx_smooth = (int)CEIL(SPREAD * FABS(dxinv) * xsig);
  const int ny_smooth = (int)CEIL(SPREAD * FABS(dyinv) * ysig);
  const int nz_smooth = (int)CEIL(SPREAD * FABS(dzinv) * zsig);
  const int nv_smooth = (int)CEIL(SPREAD * FABS(dvinv) * vsig);
  real *erfx;  erfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 2));  if( erfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfx");  }
  real *erfy;  erfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 2));  if( erfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfy");  }
  real *erfz;  erfz = (real *)malloc(sizeof(real) * (2 * nz_smooth + 2));  if( erfz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfz");  }
  real *erfv;  erfv = (real *)malloc(sizeof(real) * (2 * nv_smooth + 2));  if( erfv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfv");  }
  for(int ii = 0; ii < 2 * nx_smooth + 2; ii++)    erfx[ii] = ERF((dx * ((real)(ii - nx_smooth) - HALF)) * invxsig);
  for(int ii = 0; ii < 2 * ny_smooth + 2; ii++)    erfy[ii] = ERF((dy * ((real)(ii - ny_smooth) - HALF)) * invysig);
  for(int ii = 0; ii < 2 * nz_smooth + 2; ii++)    erfz[ii] = ERF((dz * ((real)(ii - nz_smooth) - HALF)) * invzsig);
  for(int ii = 0; ii < 2 * nv_smooth + 2; ii++)    erfv[ii] = ERF((dv * ((real)(ii - nv_smooth) - HALF)) * invvsig);

  real *psfx;  psfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 1));  if( psfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfx");  }
  real *psfy;  psfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 1));  if( psfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfy");  }
  real *psfz;  psfz = (real *)malloc(sizeof(real) * (2 * nz_smooth + 1));  if( psfz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfz");  }
  real *psfv;  psfv = (real *)malloc(sizeof(real) * (2 * nv_smooth + 1));  if( psfv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfv");  }
  for(int ii = 0; ii < 2 * nx_smooth + 1; ii++)    psfx[ii] = HALF * (erfx[ii + 1] - erfx[ii]);
  for(int ii = 0; ii < 2 * ny_smooth + 1; ii++)    psfy[ii] = HALF * (erfy[ii + 1] - erfy[ii]);
  for(int ii = 0; ii < 2 * nz_smooth + 1; ii++)    psfz[ii] = HALF * (erfz[ii + 1] - erfz[ii]);
  for(int ii = 0; ii < 2 * nv_smooth + 1; ii++)    psfv[ii] = HALF * (erfv[ii + 1] - erfv[ii]);


  const real dSxyinv = dxinv * dyinv;
  const real dSyzinv =         dyinv * dzinv;
  const real dSzxinv = dxinv	     * dzinv;

  const real df_xvinv = dxinv * dvinv;
  const real df_yvinv = dyinv * dvinv;
  const real df_zvinv = dzinv * dvinv;
#if 0
  for(int ii = 0; ii < 2 * nx_smooth + 2; ii++)    fprintf(stdout, "erfx[%d] = %e\n", ii, erfx[ii]);
  for(int ii = 0; ii < 2 * ny_smooth + 2; ii++)    fprintf(stdout, "erfy[%d] = %e\n", ii, erfy[ii]);
  for(int ii = 0; ii < 2 * nz_smooth + 2; ii++)    fprintf(stdout, "erfz[%d] = %e\n", ii, erfz[ii]);
  for(int ii = 0; ii < 2 * nx_smooth + 1; ii++)	   fprintf(stdout, "psfx[%d] = %e\n", ii, psfx[ii]);
  for(int ii = 0; ii < 2 * ny_smooth + 1; ii++)	   fprintf(stdout, "psfy[%d] = %e\n", ii, psfy[ii]);
  for(int ii = 0; ii < 2 * nz_smooth + 1; ii++)	   fprintf(stdout, "psfz[%d] = %e\n", ii, psfz[ii]);
  __FPRINTF__(stdout, "dSxyinv = %e, dSyzinv = %e, dSzxinv = %e\n", dSxyinv, dSyzinv, dSzxinv);
#endif

  for(int kk = 0; kk < kind; kk++){
    for(int ii = 0; ii < nx * ny; ii++)
      Sigma_xy[INDEX2D(kind, nx * ny, kk, ii)] = ZERO;
    for(int ii = 0; ii < ny * nz; ii++)
      Sigma_yz[INDEX2D(kind, ny * nz, kk, ii)] = ZERO;
    for(int ii = 0; ii < nz * nx; ii++)
      Sigma_zx[INDEX2D(kind, nz * nx, kk, ii)] = ZERO;
    for(int ii = 0; ii < nx * nv; ii++)
      f_xv[INDEX2D(kind, nx * nv, kk, ii)] = ZERO;
    for(int ii = 0; ii < ny * nv; ii++)
      f_yv[INDEX2D(kind, ny * nv, kk, ii)] = ZERO;
    for(int ii = 0; ii < nz * nv; ii++)
      f_zv[INDEX2D(kind, nz * nv, kk, ii)] = ZERO;

    for(int ii = bodyHead[kk]; ii < bodyHead[kk] + bodyNum[kk]; ii++){
      const real mi = CAST_D2R(CAST_R2D(body[ii].m) * mass2astro);
      const real xi = body[ii].x;
      const real yi = body[ii].y;
      const real zi = body[ii].z;
      const real vx = body[ii].vx;
      const real vy = body[ii].vy;
      const real vz = body[ii].vz;

      const int l0 = (int)FLOOR((xi - xmin) * dxinv);
      const int m0 = (int)FLOOR((yi - ymin) * dyinv);
      const int n0 = (int)FLOOR((zi - zmin) * dzinv);
      const int o0 = (int)FLOOR((vx - vmin) * dvinv);
      const int p0 = (int)FLOOR((vy - vmin) * dvinv);
      const int q0 = (int)FLOOR((vz - vmin) * dvinv);


      for(int sx = 0; sx < 2 * nx_smooth + 1; sx++){
	const int ll = l0 + sx - nx_smooth;
	if( (ll >= 0) && (ll < nx) ){
	  const real mx = mi * psfx[sx];

	  for(int sv = 0; sv < 2 * nv_smooth + 1; sv++){
	    const int oo = o0 + sv - nv_smooth;
	    if( (oo >= 0) && (oo < nv) )
	      f_xv[INDEX3D(kind, nx, nv, kk, ll, oo)] += mx * psfv[sv];
	  }/* for(int sv = 0; sv < 2 * nv_smooth + 1; sv++){ */

	  for(int sy = 0; sy < 2 * ny_smooth + 1; sy++){
	    const int mm = m0 + sy - ny_smooth;
	    if( (mm >= 0) && (mm < ny) )
	      Sigma_xy[INDEX(kind, nx, ny, kk, ll, mm)] += mx * psfy[sy];
	  }/* for(int sy = 0; sy < 2 * ny_smooth + 1; sy++){ */

	  for(int sz = 0; sz < 2 * nz_smooth + 1; sz++){
	    const int nn = n0 + sz - nz_smooth;
	    if( (nn >= 0) && (nn < nz) )
	      Sigma_zx[INDEX3D(kind, nz, nx, kk, nn, ll)] += mx * psfz[sz];
	  }/* for(int sz = 0; sz < 2 * nz_smooth + 1; sz++){ */
	}/* if( (ll >= 0) && (ll < nx) ){ */
      }/* for(int sx = 0; sx < 2 * nx_smooth + 1; sx++){ */

      for(int sy = 0; sy < 2 * ny_smooth + 1; sy++){
	const int mm = m0 + sy - ny_smooth;
	if( (mm >= 0) && (mm < ny) ){
	  const real my = mi * psfy[sy];

	  for(int sv = 0; sv < 2 * nv_smooth + 1; sv++){
	    const int pp = p0 + sv - nv_smooth;
	    if( (pp >= 0) && (pp < nv) )
	      f_yv[INDEX3D(kind, ny, nv, kk, mm, pp)] += my * psfv[sv];
	  }/* for(int sv = 0; sv < 2 * nv_smooth + 1; sv++){ */

	  for(int sz = 0; sz < 2 * nz_smooth + 1; sz++){
	    const int nn = n0 + sz - nz_smooth;
	    if( (nn >= 0) && (nn < nz) )
	      Sigma_yz[INDEX3D(kind, ny, nz, kk, mm, nn)] += my * psfz[sz];
	  }/* for(int sz = 0; sz < 2 * nz_smooth + 1; sz++){ */
	}/* if( (mm >= 0) && (mm < ny) ){ */
      }/* for(int sy = 0; sy < 2 * ny_smooth + 1; sy++){ */

      for(int sz = 0; sz < 2 * nz_smooth + 1; sz++){
	const int nn = n0 + sz - nz_smooth;
	if( (nn >= 0) && (nn < nz) ){
	  const real mz = mi * psfz[sz];
	  for(int sv = 0; sv < 2 * nv_smooth + 1; sv++){
	    const int qq = q0 + sv - nv_smooth;
	    if( (qq >= 0) && (qq < nv) )
	      f_zv[INDEX3D(kind, nz, nv, kk, nn, qq)] += mz * psfv[sv];
	  }/* for(int sv = 0; sv < 2 * nv_smooth + 1; sv++){ */
	}/* if( (nn >= 0) && (nn < nz) ){ */
      }/* for(int sz = 0; sz < 2 * nz_smooth + 1; sz++){ */

    }/* for(int ii = 0; ii < bodyNum[kk]; ii++){ */

    for(int ii = 0; ii < nx * ny; ii++)
      Sigma_xy[INDEX2D(kind, nx * ny, kk, ii)] *= dSxyinv;
    for(int ii = 0; ii < ny * nz; ii++)
      Sigma_yz[INDEX2D(kind, ny * nz, kk, ii)] *= dSyzinv;
    for(int ii = 0; ii < nz * nx; ii++)
      Sigma_zx[INDEX2D(kind, nz * nx, kk, ii)] *= dSzxinv;

    for(int ii = 0; ii < nx * nv; ii++)
      f_xv[INDEX2D(kind, nx * nv, kk, ii)] *= df_xvinv;
    for(int ii = 0; ii < ny * nv; ii++)
      f_yv[INDEX2D(kind, ny * nv, kk, ii)] *= df_yvinv;
    for(int ii = 0; ii < nz * nv; ii++)
      f_zv[INDEX2D(kind, nz * nv, kk, ii)] *= df_zvinv;
  }/* for(int kk = 0; kk < kind; kk++){ */

  free(erfx);  free(erfy);  free(erfz);  free(erfv);
  free(psfx);  free(psfy);  free(psfz);  free(psfv);


  __NOTE__("%s\n", "end");
}


#ifdef  USE_HDF5_FORMAT
#ifdef  PREPARE_XDMF_FILES
static inline void writeXdmf(char file[], const uint id, const int nx, const int ny, const int nz, const int kind)
{
  __NOTE__("%s\n", "start");

  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.xmf", DATAFOLDER, file, "plt", id);
  FILE *fp;
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  sprintf(filename, "%s.%s%.3u.h5", file, "plt", id);

  fprintf(fp, "<?xml version=\"1.0\" ?>\n");
  fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(fp, "<Xdmf Version=\"3.0\">\n");
  fprintf(fp, "<Domain>\n");
  for(int kk = 0; kk < kind; kk++){
    fprintf(fp, "<Grid Name=\"field%d\" GridType=\"Uniform\">\n", kk);
    fprintf(fp, "<Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n", nx + 1, ny + 1, nz + 1);
    fprintf(fp, "<Geometry GeometryType=\"VXVYVZ\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nx + 1, sizeof(real));
    fprintf(fp, "%s:/field%d/x\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", ny + 1, sizeof(real));
    fprintf(fp, "%s:/field%d/y\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nz + 1, sizeof(real));
    fprintf(fp, "%s:/field%d/z\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "</Geometry>\n");
    fprintf(fp, "<Attribute Name=\"rho\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nx, ny, nz, sizeof(real));
    fprintf(fp, "%s:/field%d/rho\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "</Attribute>\n");
    fprintf(fp, "</Grid>\n");

    fprintf(fp, "<Grid Name=\"field%d_xy\" GridType=\"Uniform\">\n", kk);
    fprintf(fp, "<Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", nx + 1, ny + 1);
    /**# Is VXVY available in XDMF? */
    fprintf(fp, "<Geometry GeometryType=\"VXVY\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nx + 1, sizeof(real));
    fprintf(fp, "%s:/field%d/x\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", ny + 1, sizeof(real));
    fprintf(fp, "%s:/field%d/y\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    /* fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nz + 1, sizeof(real)); */
    /* fprintf(fp, "%s:/field%d/z\n", filename, kk); */
    /* fprintf(fp, "</DataItem>\n"); */
    fprintf(fp, "</Geometry>\n");
    fprintf(fp, "<Attribute Name=\"Sigma\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nx, ny, sizeof(real));
    fprintf(fp, "%s:/field%d/Sigma_xy\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "</Attribute>\n");
    fprintf(fp, "</Grid>\n");

    fprintf(fp, "<Grid Name=\"field%d_yz\" GridType=\"Uniform\">\n", kk);
    fprintf(fp, "<Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny + 1, nz + 1);
    /**# Is VXVY available in XDMF? */
    fprintf(fp, "<Geometry GeometryType=\"VXVY\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", ny + 1, sizeof(real));
    fprintf(fp, "%s:/field%d/y\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nz + 1, sizeof(real));
    fprintf(fp, "%s:/field%d/z\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    /* fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nx + 1, sizeof(real)); */
    /* fprintf(fp, "%s:/field%d/x\n", filename, kk); */
    /* fprintf(fp, "</DataItem>\n"); */
    fprintf(fp, "</Geometry>\n");
    fprintf(fp, "<Attribute Name=\"Sigma\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", ny, nz, sizeof(real));
    fprintf(fp, "%s:/field%d/Sigma_yz\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "</Attribute>\n");
    fprintf(fp, "</Grid>\n");

    fprintf(fp, "<Grid Name=\"field%d_zx\" GridType=\"Uniform\">\n", kk);
    fprintf(fp, "<Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", nz + 1, nx + 1);
    /**# Is VXVY available in XDMF? */
    fprintf(fp, "<Geometry GeometryType=\"VXVY\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nz + 1, sizeof(real));
    fprintf(fp, "%s:/field%d/z\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nx + 1, sizeof(real));
    fprintf(fp, "%s:/field%d/x\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    /* fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", ny + 1, sizeof(real)); */
    /* fprintf(fp, "%s:/field%d/y\n", filename, kk); */
    /* fprintf(fp, "</DataItem>\n"); */
    fprintf(fp, "</Geometry>\n");
    fprintf(fp, "<Attribute Name=\"Sigma\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", nz, nx, sizeof(real));
    fprintf(fp, "%s:/field%d/Sigma_zx\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "</Attribute>\n");
    fprintf(fp, "</Grid>\n");
  }/* for(int kk = 0; kk < kind; kk++){ */
  fprintf(fp, "</Domain>\n");
  fprintf(fp, "</Xdmf>\n");

  fclose(fp);

  __NOTE__("%s\n", "end");
}
#endif//PREPARE_XDMF_FILES


/**
 * @fn writeAnalyzedProfiles
 *
 * @brief Write analyzed profiles of the N-body simulation.
 */
void writeAnalyzedProfiles
(const double time, const ulong steps, char file[], const uint id, hdf5struct type, const int kind, int * restrict bodyNum,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nz, real * restrict zz, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx,
 const int nv, real * restrict vv, real * restrict f_xv, real * restrict f_yv, real * restrict f_zv,
 const int nx3D, real * restrict rho_xx, const int ny3D, real * restrict rho_yy, const int nz3D, real * restrict rho_zz, real * restrict rho_map,
 int * restrict prfHead, int * restrict prfNum, real * restrict rad, real * restrict rho, real * restrict enc, real * restrict sig,
 int * restrict prfHead_hor, int * restrict prfNum_hor, real * restrict hor, real * restrict Sigma, real * restrict height, real * restrict sigR, real * restrict sigp, real * restrict sigz,
 real * restrict rhalf, real * restrict com, real * restrict vel)
{
  __NOTE__("%s\n", "start");

  static bool firstCall = true;
  const double phase_space_density2astro = mass2astro / (length2astro * velocity2astro);

  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "plt", id);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


  /* create the data space for the dataset */
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_FILE_COMPRESSION
  const hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  hsize_t cdims_loc[3];
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  const uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  const uint szip_pixels_per_block = 8;
  const hsize_t szip_cdims[3] = {1, 1, 128 * szip_pixels_per_block};
#endif//USE_SZIP_COMPRESSION
#ifdef  USE_GZIP_COMPRESSION
  /* compression using gzip */
  const uint gzip_compress_level = 9;/**< 9 is the maximum compression ratio */
  const hsize_t gzip_cdims[3] = {1, 1, 1024};
#endif//USE_GZIP_COMPRESSION
#else///USE_FILE_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_FILE_COMPRESSION



  /* write attribute data */
  /* create the data space for the attribute */
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  /* write current time */
  double wtime = time * time2astro;
  attribute = H5Acreate(target, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &wtime));
  chkHDF5err(H5Aclose(attribute));
  /* write # of steps */
  attribute = H5Acreate(target, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  /* write # of N-body particles */
  ulong num_ulong = 0;
  for(int ii = 0; ii < kind; ii++)
    num_ulong += (ulong)bodyNum[ii];
  attribute = H5Acreate(target, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &num_ulong));
  chkHDF5err(H5Aclose(attribute));
  /* write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &kind));
  chkHDF5err(H5Aclose(attribute));
  /* write # of grid points */
  attribute = H5Acreate(target, "nx", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nx));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "ny", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &ny));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nz", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nz));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nv", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nv));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nx3D", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nx3D));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "ny3D", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &ny3D));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nz3D", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nz3D));
  chkHDF5err(H5Aclose(attribute));
  /* write flag about DOUBLE_PRECISION */
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
  attribute = H5Acreate(target, "length_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, length_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "mass_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, mass_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "time_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, time_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "velocity_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, velocity_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, density_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "col_density_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, col_density_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));


  for(int ii = 0; ii < kind; ii++){
    char grp[16];    sprintf(grp, "field%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t attr_dims = 1;
    hid_t attribute;

    /** 3D (nx * ny * nz) array */
    hsize_t dims[3] = {nx3D, ny3D, nz3D};
    dataspace = H5Screate_simple(3, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    cdims_loc[0] = szip_cdims[0];
    cdims_loc[1] = szip_cdims[1];
    cdims_loc[2] = szip_cdims[2];
    if( (dims[0] * dims[1] * dims[2]) > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[2] < cdims_loc[2] ){
	cdims_loc[2] = dims[2];
	cdims_loc[1] = dims[2] / cdims_loc[2];
      }/* if( dims[2] < cdims_loc[2] ){ */
      if( dims[1] < cdims_loc[1] ){
	cdims_loc[1] = dims[1];
	cdims_loc[0] = dims[1] / cdims_loc[1];
      }/* if( dims[1] < cdims_loc[1] ){ */
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];
      if( cdims_loc[0] * cdims_loc[1] * cdims_loc[2] > cdims_max )
	cdims_loc[0] = cdims_max / (cdims_loc[1] * cdims_loc[2]);
      chkHDF5err(H5Pset_chunk(property, 3, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( (dims[0] * dims[1] * dims[2]) > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#else///USE_SZIP_COMPRESSION
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[0];
    cdims_loc[1] = gzip_cdims[1];
    cdims_loc[2] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[2] < cdims_loc[2] ){
      cdims_loc[2] = dims[2];
      cdims_loc[1] = dims[2] / cdims_loc[2];
    }/* if( dims[1] < cdims_loc[1] ){ */
    if( dims[1] < cdims_loc[1] ){
      cdims_loc[1] = dims[1];
      cdims_loc[0] = dims[1] / cdims_loc[1];
    }/* if( dims[1] < cdims_loc[1] ){ */
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] * cdims_loc[2] > cdims_max )
      cdims_loc[0] = cdims_max / (cdims_loc[1] * cdims_loc[2]);
    chkHDF5err(H5Pset_chunk(property, 3, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
#endif//USE_SZIP_COMPRESSION
    for(int jj = INDEX2D(kind, nx3D * ny3D * nz3D, ii, 0); jj < INDEX2D(kind, nx3D * ny3D * nz3D, ii + 1, 0); jj++)
      rho_map[jj] = CAST_D2R(CAST_R2D(rho_map[jj]) * density2astro);
    dataset = H5Dcreate(group, "rho", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rho_map[INDEX2D(kind, nx3D * ny3D * nz3D, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    if( (dims[0] * dims[1] * dims[2]) > (hsize_t)szip_pixels_per_block )
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (nx * ny) array */
    dims[0] = nx;
    dims[1] = ny;
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    cdims_loc[1] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[1] < cdims_loc[1] ){
      cdims_loc[1] = dims[1];
      cdims_loc[0] = dims[1] / cdims_loc[1];
    }/* if( dims[1] < cdims_loc[1] ){ */
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] > cdims_max )
      cdims_loc[0] = cdims_max / cdims_loc[1];
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    for(int jj = INDEX2D(kind, nx * ny, ii, 0); jj < INDEX2D(kind, nx * ny, ii + 1, 0); jj++)
      Sigma_xy[jj] = CAST_D2R(CAST_R2D(Sigma_xy[jj]) * col_density2astro);
    dataset = H5Dcreate(group, "Sigma_xy", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Sigma_xy[INDEX2D(kind, nx * ny, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (ny * nz) array */
    dims[0] = ny;
    dims[1] = nz;
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    cdims_loc[1] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[1] < cdims_loc[1] ){
      cdims_loc[1] = dims[1];
      cdims_loc[0] = dims[1] / cdims_loc[1];
    }/* if( dims[1] < cdims_loc[1] ){ */
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] > cdims_max )
      cdims_loc[0] = cdims_max / cdims_loc[1];
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    for(int jj = INDEX2D(kind, ny * nz, ii, 0); jj < INDEX2D(kind, ny * nz, ii + 1, 0); jj++)
      Sigma_yz[jj] = CAST_D2R(CAST_R2D(Sigma_yz[jj]) * col_density2astro);
    dataset = H5Dcreate(group, "Sigma_yz", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Sigma_yz[INDEX2D(kind, ny * nz, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (nz * nx) array */
    dims[0] = nz;
    dims[1] = nx;
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    cdims_loc[1] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[1] < cdims_loc[1] ){
      cdims_loc[1] = dims[1];
      cdims_loc[0] = dims[1] / cdims_loc[1];
    }/* if( dims[1] < cdims_loc[1] ){ */
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] > cdims_max )
      cdims_loc[0] = cdims_max / cdims_loc[1];
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    for(int jj = INDEX2D(kind, nz * nx, ii, 0); jj < INDEX2D(kind, nz * nx, ii + 1, 0); jj++)
      Sigma_zx[jj] = CAST_D2R(CAST_R2D(Sigma_zx[jj]) * col_density2astro);
    dataset = H5Dcreate(group, "Sigma_zx", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Sigma_zx[INDEX2D(kind, nz * nx, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (nx * nv) array */
    dims[0] = nx;
    dims[1] = nv;
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    cdims_loc[1] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[1] < cdims_loc[1] ){
      cdims_loc[1] = dims[1];
      cdims_loc[0] = dims[1] / cdims_loc[1];
    }/* if( dims[1] < cdims_loc[1] ){ */
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] > cdims_max )
      cdims_loc[0] = cdims_max / cdims_loc[1];
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    for(int jj = INDEX2D(kind, nx * nv, ii, 0); jj < INDEX2D(kind, nx * nv, ii + 1, 0); jj++)
      f_xv[jj] = CAST_D2R(CAST_R2D(f_xv[jj]) * phase_space_density2astro);
    dataset = H5Dcreate(group, "f_xv", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &f_xv[INDEX2D(kind, nx * nv, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (ny * nv) array */
    dims[0] = ny;
    dims[1] = nv;
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    cdims_loc[1] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[1] < cdims_loc[1] ){
      cdims_loc[1] = dims[1];
      cdims_loc[0] = dims[1] / cdims_loc[1];
    }/* if( dims[1] < cdims_loc[1] ){ */
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] > cdims_max )
      cdims_loc[0] = cdims_max / cdims_loc[1];
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    for(int jj = INDEX2D(kind, ny * nv, ii, 0); jj < INDEX2D(kind, ny * nv, ii + 1, 0); jj++)
      f_yv[jj] = CAST_D2R(CAST_R2D(f_yv[jj]) * phase_space_density2astro);
    dataset = H5Dcreate(group, "f_yv", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &f_yv[INDEX2D(kind, ny * nv, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (nz * nv) array */
    dims[0] = nz;
    dims[1] = nv;
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    cdims_loc[1] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[1] < cdims_loc[1] ){
      cdims_loc[1] = dims[1];
      cdims_loc[0] = dims[1] / cdims_loc[1];
    }/* if( dims[1] < cdims_loc[1] ){ */
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] > cdims_max )
      cdims_loc[0] = cdims_max / cdims_loc[1];
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    for(int jj = INDEX2D(kind, nz * nv, ii, 0); jj < INDEX2D(kind, nz * nv, ii + 1, 0); jj++)
      f_zv[jj] = CAST_D2R(CAST_R2D(f_zv[jj]) * phase_space_density2astro);
    dataset = H5Dcreate(group, "f_zv", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &f_zv[INDEX2D(kind, nz * nv, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nx) array */
    dims[0] = nx + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    if( firstCall )
      for(int jj = 0; jj < nx + 1; jj++)
	xx[jj] = CAST_D2R(CAST_R2D(xx[jj]) * length2astro);
    dataset = H5Dcreate(group, "x", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, xx));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (ny) array */
    dims[0] = ny + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    if( firstCall )
      for(int jj = 0; jj < ny + 1; jj++)
	yy[jj] = CAST_D2R(CAST_R2D(yy[jj]) * length2astro);
    dataset = H5Dcreate(group, "y", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, yy));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nz) array */
    dims[0] = nz + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    if( firstCall )
      for(int jj = 0; jj < nz + 1; jj++)
	zz[jj] = CAST_D2R(CAST_R2D(zz[jj]) * length2astro);
    dataset = H5Dcreate(group, "z", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, zz));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nv) array */
    dims[0] = nv + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    if( firstCall )
      for(int jj = 0; jj < nv + 1; jj++)
	vv[jj] = CAST_D2R(CAST_R2D(vv[jj]) * velocity2astro);
    dataset = H5Dcreate(group, "v", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, vv));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nx3D) array */
    dims[0] = nx3D + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    if( firstCall )
      for(int jj = 0; jj < nx3D + 1; jj++)
	rho_xx[jj] = CAST_D2R(CAST_R2D(rho_xx[jj]) * length2astro);
    dataset = H5Dcreate(group, "rho_x", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho_xx));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (ny3D) array */
    dims[0] = ny3D + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    if( firstCall )
      for(int jj = 0; jj < ny3D + 1; jj++)
	rho_yy[jj] = CAST_D2R(CAST_R2D(rho_yy[jj]) * length2astro);
    dataset = H5Dcreate(group, "rho_y", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho_yy));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nz3D) array */
    dims[0] = nz3D + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[2];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    if( firstCall )
      for(int jj = 0; jj < nz3D + 1; jj++)
	rho_zz[jj] = CAST_D2R(CAST_R2D(rho_zz[jj]) * length2astro);
    dataset = H5Dcreate(group, "rho_z", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho_zz));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));


    if( bodyNum[ii] > Nminimum ){
      chkHDF5err(H5Gclose(group));
      sprintf(grp, "rad%d", ii);
      group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      /** 1D (num) arrays */
      dims[0] = prfNum[ii];
      dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
      cdims_loc[0] = gzip_cdims[2];
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];
      if( cdims_loc[0] > cdims_max )
	cdims_loc[0] = cdims_max;
      chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
      chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
      for(int jj = prfHead[ii]; jj < prfHead[ii] + prfNum[ii]; jj++){
	rad[jj] = CAST_D2R(CAST_R2D(rad[jj]) *   length2astro);
	rho[jj] = CAST_D2R(CAST_R2D(rho[jj]) *  density2astro);
	enc[jj] = CAST_D2R(CAST_R2D(enc[jj]) *     mass2astro);
	sig[jj] = CAST_D2R(CAST_R2D(sig[jj]) * velocity2astro);
      }/* for(int jj = prfHead[ii]; jj < prfHead[ii] + prfNum[ii]; jj++){ */
      dataset = H5Dcreate(group, "rad", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rad[prfHead[ii]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dcreate(group, "rho", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rho[prfHead[ii]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dcreate(group, "enc", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &enc[prfHead[ii]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dcreate(group, "sig", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sig[prfHead[ii]]));
      chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));
      /* write attribute data */
      attr_dims = 1;
      dataspace = H5Screate_simple(1, &attr_dims, NULL);
      /* write # of grid points */
      attribute = H5Acreate(group, "Nrad", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &prfNum[ii]));
      chkHDF5err(H5Aclose(attribute));
      chkHDF5err(H5Sclose(dataspace));


      chkHDF5err(H5Gclose(group));
      sprintf(grp, "hor%d", ii);
      group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      /** 1D (num_hor) arrays */
      dims[0] = prfNum_hor[ii];
      dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
      cdims_loc[0] = gzip_cdims[2];
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];
      if( cdims_loc[0] > cdims_max )
	cdims_loc[0] = cdims_max;
      chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
      chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
      for(int jj = prfHead_hor[ii]; jj < prfHead_hor[ii] + prfNum_hor[ii]; jj++){
	hor   [jj] = CAST_D2R(CAST_R2D(hor   [jj]) *      length2astro);
	Sigma [jj] = CAST_D2R(CAST_R2D(Sigma [jj]) * col_density2astro);
	height[jj] = CAST_D2R(CAST_R2D(height[jj]) *      length2astro);
	sigR  [jj] = CAST_D2R(CAST_R2D(sigR  [jj]) *    velocity2astro);
	sigp  [jj] = CAST_D2R(CAST_R2D(sigp  [jj]) *    velocity2astro);
	sigz  [jj] = CAST_D2R(CAST_R2D(sigz  [jj]) *    velocity2astro);
      }/* for(int jj = prfHead[ii]; jj < prfHead[ii] + prfNum[ii]; jj++){ */
      dataset = H5Dcreate(group, "hor", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &hor[prfHead_hor[ii]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dcreate(group, "Sigma", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Sigma[prfHead_hor[ii]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dcreate(group, "height", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &height[prfHead_hor[ii]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dcreate(group, "sigR", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sigR[prfHead_hor[ii]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dcreate(group, "sigp", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sigp[prfHead_hor[ii]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dcreate(group, "sigz", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sigz[prfHead_hor[ii]]));
      chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));
      /* write attribute data */
      attr_dims = 1;
      dataspace = H5Screate_simple(1, &attr_dims, NULL);
      /* write # of grid points */
      attribute = H5Acreate(group, "Nhor", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &prfNum_hor[ii]));
      chkHDF5err(H5Aclose(attribute));
      chkHDF5err(H5Sclose(dataspace));
    }/* if( bodyNum[ii] > Nminimum ){ */

    /* write attribute data */
    chkHDF5err(H5Gclose(group));
    sprintf(grp, "attr%d", ii);
    group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    /* write # of N-body particles */
    int num = bodyNum[ii];
    attribute = H5Acreate(group, "number", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &num));
    chkHDF5err(H5Aclose(attribute));
    /* write half-mass radius */
    rhalf[ii] = CAST_D2R(CAST_R2D(rhalf[ii]) * length2astro);
    attribute = H5Acreate(group, "rhalf", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, type.real, &rhalf[ii]));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Sclose(dataspace));
    attr_dims = 3;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    for(int jj = INDEX2D(kind, 3, ii, 0); jj < INDEX2D(kind, 3, ii + 1, 0); jj++){
      com[jj] = CAST_D2R(CAST_R2D(com[jj]) *   length2astro);
      vel[jj] = CAST_D2R(CAST_R2D(vel[jj]) * velocity2astro);
    }/* for(int jj = INDEX2D(kind, 3, ii, 0); jj < INDEX2D(kind, 3, ii + 1, 0); jj++){ */
    attribute = H5Acreate(group, "com", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, type.real, &com[INDEX2D(kind, 3, ii, 0)]));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "vel", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, type.real, &vel[INDEX2D(kind, 3, ii, 0)]));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Sclose(dataspace));


    chkHDF5err(H5Gclose(group));
  }/* for(int ii = 0; ii < kind; ii++){ */

  /* close the file */
  chkHDF5err(H5Fclose(target));
#ifdef  PREPARE_XDMF_FILES
  writeXdmf(file, id, nx, ny, nz, kind);
#endif//PREPARE_XDMF_FILES

  firstCall = false;

  __NOTE__("%s\n", "end");
}


#ifdef  HDF5_FOR_ZINDAIJI
void writeZindaijiFile(const int Ntot, nbody_hdf5 hdf5, const real eps, const int kind, int *bodyHead, int *bodyNum, int *type, const double time, char file[], const uint id)
{
  __NOTE__("%s\n", "start");


  float *tmp_flt;  tmp_flt = (float *)malloc(Ntot * 3 * sizeof(float));  if( tmp_flt == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_flt");  }
  int   *tmp_int;  tmp_int = (int   *)malloc(Ntot     * sizeof(int)  );  if( tmp_int == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp_int");  }
  const double vel2zin = length2astro / time2astro;
  const double acc2zin = 0.5 * vel2zin / time2astro;



  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
  sprintf(filename, "%s/%s_%s%.3u.h5", DATAFOLDER, file, "zindaiji", id);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


  /* create the data space for the dataset */
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_FILE_COMPRESSION
  const hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  hsize_t cdims_loc[2];
#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  const uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  const uint szip_pixels_per_block = 8;
  hsize_t szip_cdims[2] = {32 * szip_pixels_per_block, 3};
#else///USE_SZIP_COMPRESSION
  /* compression using gzip */
  const uint gzip_compress_level = 9;/**< 9 is the maximum compression ratio */
  hsize_t gzip_cdims[2] = {512, 3};
#endif//USE_SZIP_COMPRESSION
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#else///USE_FILE_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_FILE_COMPRESSION


  /* write attribute data */
  /* create the data space for the attribute */
  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  static const char format_ver[CONSTANTS_H_CHAR_WORDS] = "0.0";
  attribute = H5Acreate(target, "Format Version", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, format_ver));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));


  char grp[16];  sprintf(grp, "snap%d", 0);
  hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* write attribute data */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  /* write # of N-body particles */
  attribute = H5Acreate(group, "number", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &Ntot));
  chkHDF5err(H5Aclose(attribute));
  /* write current time */
  float time_float = (float)(time * time2astro);
  attribute = H5Acreate(group, "time", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_FLOAT, &time_float));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));

  /* write particle data */
  /* 2D (num * 3) array */
  hsize_t dims[2] = {Ntot, 3};
  dataspace = H5Screate_simple(2, dims, NULL);
#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#ifdef  USE_SZIP_COMPRESSION
  cdims_loc[0] = szip_cdims[0];
  cdims_loc[1] = szip_cdims[1];
  if( dims[0] * dims[1] > (hsize_t)szip_pixels_per_block ){
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] > cdims_max )
      cdims_loc[0] = cdims_max / cdims_loc[1];
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
  }/* if( dims[0] * dims[1] > (hsize_t)szip_pixels_per_block ){ */
  else
    property = H5P_DEFAULT;
#else///USE_SZIP_COMPRESSION
  cdims_loc[0] = gzip_cdims[0];
  cdims_loc[1] = gzip_cdims[1];
  property = H5Pcreate(H5P_DATASET_CREATE);
  if( dims[0] < cdims_loc[0] )
    cdims_loc[0] = dims[0];
  if( cdims_loc[0] * cdims_loc[1] > cdims_max )
    cdims_loc[0] = cdims_max / cdims_loc[1];
  chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
  chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_SZIP_COMPRESSION
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)

  /* write particle position */
  for(int ii = 0; ii < Ntot * 3; ii++)
    tmp_flt[ii] = (float)(CAST_R2D(hdf5.pos[ii]) * length2astro);
  dataset = H5Dcreate(group, "xyz", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_flt));
  chkHDF5err(H5Dclose(dataset));
  /* write particle velocity */
  for(int ii = 0; ii < Ntot * 3; ii++)
    tmp_flt[ii] = (float)(CAST_R2D(hdf5.vel[ii]) * vel2zin);
  dataset = H5Dcreate(group, "vxvyvz", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_flt));
  chkHDF5err(H5Dclose(dataset));
  /* write particle acceleration */
  for(int ii = 0; ii < Ntot * 3; ii++)
    tmp_flt[ii] = (float)(CAST_R2D(hdf5.acc[ii]) * acc2zin);
  dataset = H5Dcreate(group, "axayaz", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_flt));
  chkHDF5err(H5Dclose(dataset));

#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#ifdef  USE_SZIP_COMPRESSION
  if( dims[0] * dims[1] > (hsize_t)szip_pixels_per_block )
#endif//USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));


  /* 1D (num) arrays */
  dataspace = H5Screate_simple(1, dims, NULL);
#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
  property = H5Pcreate(H5P_DATASET_CREATE);
#ifdef  USE_SZIP_COMPRESSION
  szip_cdims[0] = 128 * szip_pixels_per_block;
  cdims_loc[0] = szip_cdims[0];
  if( dims[0] < cdims_loc[0] )
    cdims_loc[0] = dims[0];
  if( dims[0] > (hsize_t)szip_pixels_per_block ){
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
  }/* if( (hsize_t)num[ii] > (hsize_t)szip_pixels_per_block ){ */
  else
    property = H5P_DEFAULT;
#else///USE_SZIP_COMPRESSION
  gzip_cdims[0] = 1024;
  cdims_loc[0] = gzip_cdims[0];
  if( dims[0] < cdims_loc[0] )
    cdims_loc[0] = dims[0];
  if( cdims_loc[0] > cdims_max )
    cdims_loc[0] = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
  chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_SZIP_COMPRESSION
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)

  /* write particle index */
  for(int ii = 0; ii < Ntot; ii++)
    tmp_int[ii] = (int)hdf5.idx[ii];
  dataset = H5Dcreate(group, "ID", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_int));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "index", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_int));
  chkHDF5err(H5Dclose(dataset));
  /* write particle radius */
  for(int ii = 0; ii < Ntot; ii++)
    tmp_flt[ii] = CAST_R2F(eps);
  dataset = H5Dcreate(group, "radius", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_int));
  chkHDF5err(H5Dclose(dataset));
  /* write particle type */
  for(int kk = 0; kk < kind; kk++)
    for(int ii = bodyHead[kk]; ii < bodyHead[kk] + bodyNum[kk]; ii++)
      tmp_int[ii] = type[kk];
  dataset = H5Dcreate(group, "type", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_int));
  chkHDF5err(H5Dclose(dataset));

#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#ifdef  USE_SZIP_COMPRESSION
  if( dims[0] * dims[1] > (hsize_t)szip_pixels_per_block )
#endif//USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));


  chkHDF5err(H5Gclose(group));
  chkHDF5err(H5Fclose(target));

  free(tmp_flt);
  free(tmp_int);

  __NOTE__("%s\n", "end");
}
#endif//HDF5_FOR_ZINDAIJI


#endif//USE_HDF5_FORMAT
