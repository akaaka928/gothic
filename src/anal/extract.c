/**
 * @file extract.c
 *
 * @brief Analyzer for radial profiles of multiple components
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2018/03/08 (Thu)
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
void generateMassDistributionMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body_tot, const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz, real * restrict rho_map, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx);

#ifdef  USE_HDF5_FORMAT
void writeAnalyzedProfiles
(const double time, const ulong steps, char file[], const uint id, hdf5struct type, const int kind, int * restrict bodyNum,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nz, real * restrict zz, real * restrict rho_map, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx,
 int * restrict prfHead, int * restrict prfNum, real * restrict rad, real * restrict rho, real * restrict enc, real * restrict sig,
   int * restrict prfHead_hor, int * restrict prfNum_hor, real * restrict hor, real * restrict Sigma, real * restrict height, real * restrict sigR, real * restrict sigp, real * restrict sigz,
   real * restrict rhalf, real * restrict com, real * restrict vel);
#endif//USE_HDF5_FORMAT


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);

  /** initialization */
  if( argc < 15 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 15);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -ncrit=<int>\n");
    __FPRINTF__(stderr, "          -nx=<int> -xmin=<real> -xmax=<real>\n");
    __FPRINTF__(stderr, "          -ny=<int> -ymin=<real> -ymax=<real>\n");
    __FPRINTF__(stderr, "          -nz=<int> -zmin=<real> -zmax=<real>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 15 ){ */

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));
  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "ncrit", &ncrit_base));

  int nx;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nx", &nx));
  int ny;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "ny", &ny));
  int nz;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nz", &nz));
  real xmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "xmin", &xmin));
  real ymin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "ymin", &ymin));
  real zmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "zmin", &zmin));
  real xmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "xmax", &xmax));
  real ymax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "ymax", &ymax));
  real zmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "zmax", &zmax));


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


  /** read number of components */
  int kind = 0;
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
  checker &= (1 == fscanf(fp, "%d\t%*d", &kind));
  bodyHead = (int *)malloc(sizeof(int) * kind);  if( bodyHead == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate bodyHead");  }
  bodyNum  = (int *)malloc(sizeof(int) * kind);  if( bodyNum  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate bodyNum");  }
  for(int ii = 0; ii < kind; ii++)
    checker &= (1 == fscanf(fp, "%d", &bodyNum[ii]));
  fclose(fp);
  if( !checker ){
    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
  }
  bodyHead[0] = 0;
  for(int ii = 1; ii < kind; ii++)
    bodyHead[ii] = bodyHead[ii - 1] + bodyNum[ii - 1];


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
  real *xx;  xx = (real *)malloc(sizeof(real) * nx);  if( xx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xx");  }
  real *yy;  yy = (real *)malloc(sizeof(real) * ny);  if( yy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate yy");  }
  real *zz;  zz = (real *)malloc(sizeof(real) * nz);  if( zz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate zz");  }
  const real dx = (xmax - xmin) / (real)nx;  for(int ii = 0; ii < nx; ii++)    xx[ii] = xmin + dx * (HALF + (real)ii);
  const real dy = (ymax - ymin) / (real)ny;  for(int jj = 0; jj < ny; jj++)    yy[jj] = ymin + dy * (HALF + (real)jj);
  const real dz = (zmax - zmin) / (real)nz;  for(int kk = 0; kk < nz; kk++)    zz[kk] = zmin + dz * (HALF + (real)kk);

  real *rho_map;
  rho_map = (real *)malloc(sizeof(real) * kind * nx * ny * nz);  if( rho_map == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rho_map");  }
  real *Sigma_xy;  Sigma_xy = (real *)malloc(sizeof(real) * kind * nx * ny);  if( Sigma_xy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Sigma_xy");  }
  real *Sigma_yz;  Sigma_yz = (real *)malloc(sizeof(real) * kind * ny * nz);  if( Sigma_yz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Sigma_yz");  }
  real *Sigma_zx;  Sigma_zx = (real *)malloc(sizeof(real) * kind * nz * nx);  if( Sigma_zx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Sigma_zx");  }



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
    sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "split", filenum);
    bool dump_splitted_snapshot = (0 != access(filename, F_OK));
    if( !dump_splitted_snapshot ){
      struct stat stat_split;
      stat(filename, &stat_split);
      char tmpname[128];
      sprintf(tmpname, "%s/%s.%s%.3u.h5", DATAFOLDER, file, SNAPSHOT, filenum);
      struct stat stat_snap;
      stat(tmpname, &stat_snap);
      if( stat_snap.st_ctime > stat_split.st_ctime )
	dump_splitted_snapshot = true;
    }/* if( !dump_splitted_snapshot ){ */
    if( dump_splitted_snapshot ){
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

      writeSnapshotMultiGroups(time, steps, &hdf5, file, filenum, hdf5type, kind, bodyHead, bodyNum);
    }/* if( dump_splitted_snapshot ){ */
    else{
      fprintf(stderr, "# \"%s\" was not updated for reducing the elapsed time.\n", filename);
      fflush(stderr);
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
    generateMassDistributionMaps(kind, bodyHead, bodyNum, body, eps, nx, xmin, dx, ny, ymin, dy, nz, zmin, dz, rho_map, Sigma_xy, Sigma_yz, Sigma_zx);


    /** dump analyzed results for matplotlib and/or VisIt */
#ifdef  USE_HDF5_FORMAT
    writeAnalyzedProfiles
      (time, steps, file, filenum, hdf5type, kind, bodyNum,
       nx, xx, ny, yy, nz, zz, rho_map, Sigma_xy, Sigma_yz, Sigma_zx,
       prfHead, prfNum, rad, rho, enc, sig,
       prfHead_hor, prfNum_hor, hor, Sigma, height, sigR, sigp, sigz,
       &rhalf[INDEX2D(nfile, kind, (filenum - start) / interval, 0)], &comPos[INDEX(nfile, kind, 3, (filenum - start) / interval, 0, 0)], &comVel[INDEX(nfile, kind, 3, (filenum - start) / interval, 0, 0)]);
#endif//USE_HDF5_FORMAT
  }/* for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){ */


  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, comPos, nfile * kind * 3, MPI_REALDAT, MPI_SUM, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, comVel, nfile * kind * 3, MPI_REALDAT, MPI_SUM, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE,  rhalf, nfile * kind    , MPI_REALDAT, MPI_SUM, mpi.comm));
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

  free(prfHead);  free(prfNum);
  free(rad);  free(rho);  free(enc);  free(sig);

  free(prfHead_hor);  free(prfNum_hor);
  free(hor);  free(Sigma);  free(height);  free(sigR);  free(sigp);  free(sigz);

  free(xx);  free(yy);  free(zz);
  free(rho_map);
  free(Sigma_xy);  free(Sigma_yz);  free(Sigma_zx);

  exitMPI();

  return (0);
}


real getCenter(const int num, nbody_particle *body, real com[restrict], real vel[restrict])
{
  __NOTE__("%s\n", "start");


  double mtot = 0.0;
  double comx = 0.0;
  double comy = 0.0;
  double comz = 0.0;
  double velx = 0.0;
  double vely = 0.0;
  double velz = 0.0;

  bool converge = false;
  while( true ){

    for(int ii = 0; ii < num; ii++){
      const real xx = body[ii].x - CAST_D2R(comx);
      const real yy = body[ii].y - CAST_D2R(comy);
      const real zz = body[ii].z - CAST_D2R(comz);
      /* is recipe for correcting precession required?? */
      const real R2 = 1.0e-30f + xx * xx + yy * yy;
      const real r2 = R2 + zz * zz;
      body[ii].hor = R2 * RSQRT(R2);
      body[ii].rad = r2 * RSQRT(r2);
    }/* for(int ii = 0; ii < num; ii++){ */

    /** sort by particle position */
    qsort(body, num, sizeof(nbody_particle), radAscendingOrder);

    if( converge )
      break;


    mtot = 0.0;
    velx = 0.0;
    vely = 0.0;
    velz = 0.0;

    double newx = 0.0;
    double newy = 0.0;
    double newz = 0.0;

    /** find center-of-mass and bulk velocity of particles within the half-mass radius (assumption: equal-mass particles) */
    for(int ii = 0; ii < (num >> 1); ii++){
      const double xx = CAST_R2D(body[ii].x) - comx;
      const double yy = CAST_R2D(body[ii].y) - comy;
      const double zz = CAST_R2D(body[ii].z) - comz;

      const double mass = CAST_R2D(body[ii].m);
      mtot += mass;

      const double vx = CAST_R2D(body[ii].vx);
      const double vy = CAST_R2D(body[ii].vy);
      const double vz = CAST_R2D(body[ii].vz);

      newx += mass * xx;
      newy += mass * yy;
      newz += mass * zz;

      velx += mass * vx;
      vely += mass * vy;
      velz += mass * vz;
    }/* for(int ii = 0; ii < (num >> 1); ii++){ */

    mtot = 1.0 / mtot;
    newx *= mtot;
    newy *= mtot;
    newz *= mtot;
    velx *= mtot;
    vely *= mtot;
    velz *= mtot;

    const double dx = newx - comx;
    const double dy = newy - comy;
    const double dz = newz - comz;

    converge = ( (dx * dx + dy * dy + dz * dz) < 1.0e-6 );

    comx = newx;
    comy = newy;
    comz = newz;
  }/* while( true ){ */

  com[0] = CAST_D2R(comx);
  com[1] = CAST_D2R(comy);
  com[2] = CAST_D2R(comz);
  vel[0] = CAST_D2R(velx);
  vel[1] = CAST_D2R(vely);
  vel[2] = CAST_D2R(velz);


  __NOTE__("%s\n", "end");
  return (HALF * (body[num >> 1].rad + body[(num >> 1) + 1].rad));
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

    /* /\* sort the array in ascending distance from the center *\/ */
    /* qsort(body, bodyNum[kk], sizeof(nbody_particle), radAscendingOrder); */

    prfHead[kk] = *num;
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
      if( tail - 1 >= bodyNum[kk] )	break;

      real outer = body[tail - 1].rad;
      prad[*num] = (ncrit & 1) ? (body[head + (ncrit >> 1)].rad) : (HALF * (body[head + (ncrit >> 1) - 1].rad + body[head + (ncrit >> 1)].rad));/**< use median of particle location */

      real vr2 = ZERO;
      real vrm = ZERO;
      for(int ii = head; ii < tail; ii++){
	const real vr = ((body[ii].vx - vel[0]) * (body[ii].x - com[0]) + (body[ii].vy - vel[1]) * (body[ii].y - com[1]) + (body[ii].vz - vel[2]) * (body[ii].z - com[2])) / body_tot[ii].rad;
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

    /* sort the array in ascending distance from the center */
    qsort(body, bodyNum[kk], sizeof(nbody_particle), horAscendingOrder);

    prfHead[kk] = *num;
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
    }


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
      if( tail - 1 >= bodyNum[kk] )	break;

      real outer = body[tail - 1].hor;
      ppos[*num] = (ncrit & 1) ? (body[head + (ncrit >> 1)].hor) : (HALF * (body[head + (ncrit >> 1) - 1].hor + body[head + (ncrit >> 1)].hor));/**< use median of particle location */

      real mass = ZERO;
      real vR2 = ZERO;      real vRm = ZERO;
      real vp2 = ZERO;      real vpm = ZERO;
      real vz2 = ZERO;      real vzm = ZERO;
      for(int ii = head; ii < tail; ii++){
      	mass += body[ii].m;

	ver[ii - head] = body[ii].z - com[2];

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

    prfNum[kk] = *num - prfHead[kk];
  }/* for(int kk = 0; kk < kind; kk++){ */


  __NOTE__("%s\n", "end");
}


/* 3 sigma means that neglecting component of 0.26% */
/* 5 sigma means that neglecting component of 6e-7 */
/* #define SPREAD (3) */
#define SPREAD (5)
void generateMassDistributionMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body_tot, const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz, real * restrict rho_map, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx)
{
  __NOTE__("%s\n", "start");

  const real sig = HALF * eps;/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real invsig = UNITY / sig;

  const int nx_smooth = (int)CEIL(SPREAD * dx * invsig);
  const int ny_smooth = (int)CEIL(SPREAD * dy * invsig);
  const int nz_smooth = (int)CEIL(SPREAD * dz * invsig);
  real *erfx;  erfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 2));  if( erfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfx");  }
  real *erfy;  erfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 2));  if( erfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfy");  }
  real *erfz;  erfz = (real *)malloc(sizeof(real) * (2 * nz_smooth + 2));  if( erfz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfz");  }
  for(int ii = 0; ii < 2 * nx_smooth + 2; ii++)    erfx[ii] = ERF((dx * (real)(ii - nx_smooth)) * invsig);
  for(int ii = 0; ii < 2 * ny_smooth + 2; ii++)    erfy[ii] = ERF((dy * (real)(ii - ny_smooth)) * invsig);
  for(int ii = 0; ii < 2 * nz_smooth + 2; ii++)    erfz[ii] = ERF((dz * (real)(ii - nz_smooth)) * invsig);

  real *psfx;  psfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 1));  if( psfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfx");  }
  real *psfy;  psfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 1));  if( psfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfy");  }
  real *psfz;  psfz = (real *)malloc(sizeof(real) * (2 * nz_smooth + 1));  if( psfz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfz");  }
  for(int ii = 0; ii < 2 * nx_smooth + 1; ii++)    psfx[ii] = HALF * (erfx[ii + 1] - erfx[ii]);
  for(int ii = 0; ii < 2 * ny_smooth + 1; ii++)    psfy[ii] = HALF * (erfy[ii + 1] - erfy[ii]);
  for(int ii = 0; ii < 2 * nz_smooth + 1; ii++)    psfz[ii] = HALF * (erfz[ii + 1] - erfz[ii]);

  const real dxinv = UNITY / dx;
  const real dyinv = UNITY / dy;
  const real dzinv = UNITY / dz;

  const real dVinv = dxinv * dyinv * dzinv;
  const real dSxyinv = dxinv * dyinv;
  const real dSyzinv = dyinv * dzinv;
  const real dSzxinv = dzinv * dxinv;

  for(int kk = 0; kk < kind; kk++){
    for(int ii = 0; ii < nx * ny * nz; ii++)
      rho_map[INDEX2D(kind, nx * ny * nz, kk, ii)] = ZERO;
    for(int ii = 0; ii < nx * ny; ii++)
      Sigma_xy[INDEX2D(kind, nx * ny, kk, ii)] = ZERO;
    for(int ii = 0; ii < ny * nz; ii++)
      Sigma_yz[INDEX2D(kind, ny * nz, kk, ii)] = ZERO;
    for(int ii = 0; ii < nz * nx; ii++)
      Sigma_zx[INDEX2D(kind, nz * nx, kk, ii)] = ZERO;

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

      for(int ll = l0 - nx_smooth; ll < l0 + nx_smooth + 1; ll++)
	if( (ll >= 0) && (ll < nx) ){
	  const real mx = mi * psfx[ll - l0 - nx_smooth];
	  for(int mm = m0 - ny_smooth; mm < m0 + ny_smooth + 1; mm++)
	    if( (mm >= 0) && (mm < ny) ){
	      const real my = mx * psfy[mm - m0 - ny_smooth];
	      Sigma_xy[INDEX(kind, nx, ny, kk, ll, mm)] += my;
	      for(int nn = n0 - nz_smooth; nn < n0 + nz_smooth + 1; nn++)
		if( (nn >= 0) && (nn < nz) )
		  rho_map[INDEX4D(kind, nx, ny, nz, kk, ll, mm, nn)] += my * psfz[nn - n0 - nz_smooth];
	    }/* if( (mm >= 0) && (mm < ny) ){ */
	  for(int nn = n0 - nz_smooth; nn < n0 + nz_smooth + 1; nn++)
	    if( (nn >= 0) && (nn < nz) )
	      Sigma_zx[INDEX(kind, nz, nx, kk, nn, ll)] += mx * psfz[nn - n0 - nz_smooth];
	}/* if( (ll >= 0) && (ll < nx) ){ */

      for(int mm = m0 - ny_smooth; mm < m0 + ny_smooth + 1; mm++)
	if( (mm >= 0) && (mm < ny) ){
	  const real my = mi * psfy[mm - m0 - ny_smooth];
	  for(int nn = n0 - nz_smooth; nn < n0 + nz_smooth + 1; nn++)
	    if( (nn >= 0) && (nn < nz) )
	      Sigma_yz[INDEX(kind, ny, nz, kk, mm, nn)] += my * psfz[nn - n0 - nz_smooth];
	}/* if( (mm >= 0) && (mm < ny) ){ */
    }/* for(int ii = 0; ii < bodyNum[kk]; ii++){ */

    for(int ii = 0; ii < nx * ny * nz; ii++)
      rho_map[INDEX2D(kind, nx * ny * nz, kk, ii)] *= dVinv;
    for(int ii = 0; ii < nx * ny; ii++)
      Sigma_xy[INDEX2D(kind, nx * ny, kk, ii)] *= dSxyinv;
    for(int ii = 0; ii < ny * nz; ii++)
      Sigma_yz[INDEX2D(kind, ny * nz, kk, ii)] *= dSyzinv;
    for(int ii = 0; ii < nz * nx; ii++)
      Sigma_zx[INDEX2D(kind, nz * nx, kk, ii)] *= dSzxinv;
  }/* for(int cmp = 0; cmp < kind; cmp++){ */

  free(erfx);  free(erfy);  free(erfz);
  free(psfx);  free(psfy);  free(psfz);


  __NOTE__("%s\n", "end");
}


#ifdef  USE_HDF5_FORMAT
/**
 * @fn writeAnalyzedProfiles
 *
 * @brief Write analyzed profiles of the N-body simulation.
 */
void writeAnalyzedProfiles
(const double time, const ulong steps, char file[], const uint id, hdf5struct type, const int kind, int * restrict bodyNum,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nz, real * restrict zz, real * restrict rho_map, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx,
 int * restrict prfHead, int * restrict prfNum, real * restrict rad, real * restrict rho, real * restrict enc, real * restrict sig,
 int * restrict prfHead_hor, int * restrict prfNum_hor, real * restrict hor, real * restrict Sigma, real * restrict height, real * restrict sigR, real * restrict sigp, real * restrict sigz,
 real * restrict rhalf, real * restrict com, real * restrict vel)
{
  __NOTE__("%s\n", "start");

  static bool firstCall = true;

  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "plot", id);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


  /* create the data space for the dataset */
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  const uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  const uint szip_pixels_per_block = 8;
  const hsize_t cdims[3] = {1, 1, 128 * szip_pixels_per_block};
  const hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION


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
    char grp[16];    sprintf(grp, "data%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /** 3D (nx * ny * nz) array */
    hsize_t dims[3] = {nx, ny, nz};
    dataspace = H5Screate_simple(3, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    hsize_t cdims_loc[3] = {cdims[0], cdims[1], cdims[2]};
    if( (dims[0] * dims[1] * dims[2]) > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[2] < cdims_loc[2] ){
	cdims_loc[2] = dims[2];
	cdims_loc[1] = cdims[2] / cdims_loc[2];
      }/* if( dims[2] < cdims_loc[2] ){ */
      if( dims[1] < cdims_loc[1] ){
	cdims_loc[1] = dims[1];
	cdims_loc[0] = cdims[1] / cdims_loc[1];
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
#endif//USE_SZIP_COMPRESSION
    for(int jj = INDEX2D(kind, nx * ny * nz, ii, 0); jj < INDEX2D(kind, nx * ny * nz, ii + 1, 0); jj++)
      rho_map[jj] = CAST_D2R(CAST_R2D(rho_map[jj]) * density2astro);
    dataset = H5Dcreate(group, "rho_map", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rho_map[INDEX2D(kind, nx * ny * nz, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (nx * ny) array */
    dims[0] = nx;
    dims[1] = ny;
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    cdims_loc[0] = cdims[1];
    cdims_loc[1] = cdims[2];
    if( (dims[0] * dims[1]) > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[1] < cdims_loc[1] ){
	cdims_loc[1] = dims[1];
	cdims_loc[0] = cdims[1] / cdims_loc[1];
      }/* if( dims[1] < cdims_loc[1] ){ */
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];

      if( cdims_loc[0] * cdims_loc[1] > cdims_max )
	cdims_loc[0] = cdims_max / cdims_loc[1];
      chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( (dims[0] * dims[1]) > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
    for(int jj = INDEX2D(kind, nx * ny, ii, 0); jj < INDEX2D(kind, nx * ny, ii + 1, 0); jj++)
      Sigma_xy[jj] = CAST_D2R(CAST_R2D(Sigma_xy[jj]) * col_density2astro);
    dataset = H5Dcreate(group, "Sigma_xy", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Sigma_xy[INDEX2D(kind, nx * ny, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (ny * nz) array */
    dims[0] = ny;
    dims[1] = nz;
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    cdims_loc[0] = cdims[1];
    cdims_loc[1] = cdims[2];
    if( (dims[0] * dims[1]) > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[1] < cdims_loc[1] ){
	cdims_loc[1] = dims[1];
	cdims_loc[0] = cdims[1] / cdims_loc[1];
      }/* if( dims[1] < cdims_loc[1] ){ */
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];

      if( cdims_loc[0] * cdims_loc[1] > cdims_max )
	cdims_loc[0] = cdims_max / cdims_loc[1];
      chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( (dims[0] * dims[1]) > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
    for(int jj = INDEX2D(kind, ny * nz, ii, 0); jj < INDEX2D(kind, ny * nz, ii + 1, 0); jj++)
      Sigma_yz[jj] = CAST_D2R(CAST_R2D(Sigma_yz[jj]) * col_density2astro);
    dataset = H5Dcreate(group, "Sigma_yz", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Sigma_yz[INDEX2D(kind, ny * nz, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (nz * nx) array */
    dims[0] = nz;
    dims[1] = nx;
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    cdims_loc[0] = cdims[1];
    cdims_loc[1] = cdims[2];
    if( (dims[0] * dims[1]) > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[1] < cdims_loc[1] ){
	cdims_loc[1] = dims[1];
	cdims_loc[0] = cdims[1] / cdims_loc[1];
      }/* if( dims[1] < cdims_loc[1] ){ */
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];

      if( cdims_loc[0] * cdims_loc[1] > cdims_max )
	cdims_loc[0] = cdims_max / cdims_loc[1];
      chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( (dims[0] * dims[1]) > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
    for(int jj = INDEX2D(kind, nz * nx, ii, 0); jj < INDEX2D(kind, nz * nx, ii + 1, 0); jj++)
      Sigma_zx[jj] = CAST_D2R(CAST_R2D(Sigma_zx[jj]) * col_density2astro);
    dataset = H5Dcreate(group, "Sigma_zx", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Sigma_zx[INDEX2D(kind, nz * nx, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nx) array */
    dims[0] = nx;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    cdims_loc[0] = cdims[2];
    if( dims[0] > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];

      if( cdims_loc[0] > cdims_max )
	cdims_loc[0] = cdims_max;
      chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( dims[0] > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
    if( firstCall )
      for(int jj = 0; jj < nx; jj++)
	xx[jj] = CAST_D2R(CAST_R2D(xx[jj]) * length2astro);
    dataset = H5Dcreate(group, "x", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, xx));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (ny) array */
    dims[0] = ny;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    cdims_loc[0] = cdims[2];
    if( dims[0] > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];

      if( cdims_loc[0] > cdims_max )
	cdims_loc[0] = cdims_max;
      chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( dims[0] > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
    if( firstCall )
      for(int jj = 0; jj < ny; jj++)
	yy[jj] = CAST_D2R(CAST_R2D(yy[jj]) * length2astro);
    dataset = H5Dcreate(group, "y", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, yy));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nz) array */
    dims[0] = nz;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    cdims_loc[0] = cdims[2];
    if( dims[0] > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];

      if( cdims_loc[0] > cdims_max )
	cdims_loc[0] = cdims_max;
      chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( dims[0] > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
    if( firstCall )
      for(int jj = 0; jj < nz; jj++)
	zz[jj] = CAST_D2R(CAST_R2D(zz[jj]) * length2astro);
    dataset = H5Dcreate(group, "z", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, zz));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (num) arrays */
    dims[0] = prfNum[ii];
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    cdims_loc[0] = cdims[2];
    if( dims[0] > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];

      if( cdims_loc[0] > cdims_max )
	cdims_loc[0] = cdims_max;
      chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( dims[0] > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
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
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (num_hor) arrays */
    dims[0] = prfNum_hor[ii];
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    cdims_loc[0] = cdims[2];
    if( dims[0] > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( dims[0] < cdims_loc[0] )
	cdims_loc[0] = dims[0];

      if( cdims_loc[0] > cdims_max )
	cdims_loc[0] = cdims_max;
      chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( dims[0] > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
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
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));


    /* write attribute data */
    /* create the data space for the attribute */
    hsize_t attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    /* write # of N-body particles */
    int num = bodyNum[ii];
    attribute = H5Acreate(group, "number", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &num));
    chkHDF5err(H5Aclose(attribute));
    /* write # of grid points */
    attribute = H5Acreate(group, "Nrad", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &prfNum[ii]));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "Nhor", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &prfNum_hor[ii]));
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

  firstCall = false;

  __NOTE__("%s\n", "end");
}
#endif//USE_HDF5_FORMAT
