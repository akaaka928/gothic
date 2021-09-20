/**
 * @file subaru.c
 *
 * @brief Generate surface density maps in the observation fields of Subaru/HSC
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2019/12/05 (Thu)
 *
 * Copyright (C) 2019 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

/**
 * @def USE_DEGREE_FOR_SURFACE_DENSITY_MAP
 *
 * @brief On to use degree to set degree as a length scale in surface density maps instead of kpc (default is ON).
 */
#define USE_DEGREE_FOR_SURFACE_DENSITY_MAP

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
#include <mpi.h>

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
#include "rotate.h"

#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../file/io.h"
#include "../anal/m31coord.h"

#include "rand.h"
#ifdef  USE_SFMTJUMP
#include "SFMT-jump.h"
#include "sfmtjump_polynomial.h"
#endif//USE_SFMTJUMP


/** 5 fields in Subaru/HSC observations by Komiyama et al. (2018): copied from py/m31.py, 4 additional fields proposed for HSC or PFS observations: copied from py/radec.py */
const real hsc_fov = 0.75;/**< FoV is 0.75 degree in radius (HSC spec) */
const real hsc_fov2 = 0.75 * 0.75;
/* #define NFIELD (9) */
/* const real  xi0[NFIELD] = {-4.653744156858361, -5.7227710294896745, -5.530826768687439, -4.842739589150247, -5.912437203950868, -4.05, -3.2, -6.1, -6.3}; */
/* const real eta0[NFIELD] = { 3.635683362752261,  4.47091006650441  ,  5.816784796173433,  2.290423176281251,  3.125050968157858,  1.25,  0.2,  6.9,  8.2}; */
/* char *field_name[NFIELD] = {"f003", "f004", "f009", "f022", "f023", "S0", "S1", "N0", "N1"}; */
#define NFIELD (5)
const real  xi0[NFIELD] = {-4.653744156858361, -5.7227710294896745, -5.530826768687439, -4.842739589150247, -5.912437203950868/* , -4.05, -3.2, -6.1, -6.3 */};
const real eta0[NFIELD] = { 3.635683362752261,  4.47091006650441  ,  5.816784796173433,  2.290423176281251,  3.125050968157858/* ,  1.25,  0.2,  6.9,  8.2 */};
char *field_name[NFIELD] = {"f003", "f004", "f009", "f022", "f023"/* , "S0", "S1", "N0", "N1" */};


/** snapshots without DM sub-halo interaction, M = 10^7, 10^7.5, 10^8, 10^8.5, 10^9, and 10^9.5 Msun */
#define NMODEL (7)
/** all snapshot must be located in dat/ -->> creating symbolic link is required */
/* char *modelTag[NMODEL] = {"nws-continue", "nws-test-m7_0-orbit0", "nws-test-m7_5-orbit0", "nws-test-m8_0-orbit0", "nws-test-m8_5-orbit0", "nws-test-m9_0-orbit0", "nws-test-m9_5-orbit0"}; */
char *modelTag[NMODEL] = {"nws-continue", "nws-test-m7_0-orbit1", "nws-test-m7_5-orbit1", "nws-test-m8_0-orbit1", "nws-test-m8_5-orbit1", "nws-test-m9_0-orbit1", "nws-test-m9_5-orbit1"};
/* char *modelTag[NMODEL] = {"nws-continue", "nws-test-m7_0-orbit6", "nws-test-m7_5-orbit6", "nws-test-m8_0-orbit6", "nws-test-m8_5-orbit6", "nws-test-m9_0-orbit6", "nws-test-m9_5-orbit6"}; */
/* char *modelTag[NMODEL] = {"nws-continue", "nws-test-m7_0-orbit7", "nws-test-m7_5-orbit7", "nws-test-m8_0-orbit7", "nws-test-m8_5-orbit7", "nws-test-m9_0-orbit7", "nws-test-m9_5-orbit7"}; */


/** noise patterns: white noise (S/N = 1, 2, 3, 4, 5, 10), stellar halo in M31, and MW foreground contamination (stellar halo and MW are in preparation) */
#define NNOISE (6)
const int sn_val[NNOISE] = {1, 2, 3, 4, 5, 10};
const real sn_inv[NNOISE] = {UNITY, HALF, ONE_THIRD, QUARTER, ONE_FIFTH, HALF * ONE_FIFTH};

/** total number of noise particles */
#define NDIST (16384)



extern const double   length2astro;
extern const double     time2astro;
extern const double     mass2astro;
extern const double velocity2astro;


#ifdef  __ICC
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
  if(          ((const nbody_aos *)a)->idx > ((const nbody_aos *)b)->idx ){    return ( 1);  }
  else{    if( ((const nbody_aos *)a)->idx < ((const nbody_aos *)b)->idx ){    return (-1);  }
    else{                                                                      return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}
#ifdef  __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC


void read_model
(const bool without_subhalo, const uint filenum,
 const int modelID, nbody_aos *body, const int unit, const real eps,
 const ulong Ntot, const int kind, int * restrict bodyHead, int * restrict bodyNum,
 real rot[][3], real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nv, real * restrict vv,
 real * restrict map, real * restrict vmap,
 rand_state *rand
#ifdef  USE_HDF5_FORMAT
 , nbody_hdf5 hdf5, const hdf5struct hdf5type
#endif//USE_HDF5_FORMAT
 );

#ifdef  USE_HDF5_FORMAT
void writeSubaruMaps
(char file[], const uint id, hdf5struct type, const int kind,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nv, real * restrict vv,
 real * restrict map, real * restrict vmap);
#endif//USE_HDF5_FORMAT


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);

  /** initialization */
  if( argc < 8 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 8);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -nx=<int> -ny=<int> -nv=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));/**< this would be "nws-continue" */
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));

  int nx;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nx", &nx));
  int ny;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "ny", &ny));
  int nv;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nv", &nv));

  /** set coordinate transformation from M31's disk coordinate to observerd frame */
  static real rot[3][3], inv[3][3];/**< rot: M31 disk frame (X, Y, Z) ==>> observed frame (x, y, z) */
  setRotationMatrix(inv, rot);

  /** load global settings of particle distribution */
  int last;
  readConfigFile(&last, file);
  int unit;
  ulong Ntot;
  real eta_tmp, eps;
  double ft, snapshotInterval, saveInterval;
  readSettingsParallel(&unit, &Ntot, &eps, &eta_tmp, &ft, &snapshotInterval, &saveInterval, file, mpi);
  setPhysicalConstantsAndUnitSystem(unit, 1);
  eps = CAST_D2R(CAST_R2D(eps) * length2astro);

  nbody_aos *body;
  /* allocParticleDataAoS((int)Ntot, &body); */
  body = (nbody_aos *)malloc(sizeof(nbody_aos) * (Ntot + 1));/**< 1 is test-particle corresponding to DM sub-halo */
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
     (int)(Ntot + 1), &hdf5);/**< 1 is test-particle corresponding to DM sub-halo */
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
  real *  xi;  xi   = (real *)malloc(sizeof(real) * Ntot);  if(   xi == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xi.");  }
  real * eta;  eta  = (real *)malloc(sizeof(real) * Ntot);  if(  eta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate eta.");  }
  real *dist;  dist = (real *)malloc(sizeof(real) * Ntot);  if( dist == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate dist.");  }
  real * vxi;  vxi  = (real *)malloc(sizeof(real) * Ntot);  if( vxi  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vxi.");  }
  real *veta;  veta = (real *)malloc(sizeof(real) * Ntot);  if( veta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate veta.");  }
  real *vlos;  vlos = (real *)malloc(sizeof(real) * Ntot);  if( vlos == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vlos.");  }



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



  real *map;  map = (real *)malloc(sizeof(real) * (kind * (NMODEL + NNOISE)) * NFIELD * nx * ny);  if( map == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate map");  }
  real *vmap;  vmap = (real *)malloc(sizeof(real) * (kind * (NMODEL + NNOISE)) * NFIELD * nv * ny);  if( vmap == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vmap");  }

  real *xx;  xx = (real *)malloc(sizeof(real) * NFIELD * (nx + 1));  if( xx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xx");  }
  real *yy;  yy = (real *)malloc(sizeof(real) * NFIELD * (ny + 1));  if( yy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate yy");  }
  real *vv;  vv = (real *)malloc(sizeof(real) * NFIELD * (nv + 1));  if( vv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vv");  }



  /** initialize pseudo random number generator */
  rand_state *rand;
  initRandNum(&rand);
#ifdef  USE_SFMTJUMP
  for(int ii = 0; ii < mpi.rank; ii++)
    SFMT_jump(rand, SFMTJUMP_10_100);
#endif//USE_SFMTJUMP


  for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){
    /** analyze snapshot without DM sub-halo and create noise maps*/
    /** analyze snapshot with various mass of DM sub-halo */
    for(int ii = 0; ii < NMODEL; ii++)
      read_model(ii == 0, filenum, ii, body, unit, eps,
		 Ntot, kind, bodyHead, bodyNum,
		 rot, xi, eta, dist, vxi, veta, vlos,
		 nx, xx, ny, yy, nv, vv, map, vmap,
		 rand
#ifdef  USE_HDF5_FORMAT
		 , hdf5, hdf5type
#endif//USE_HDF5_FORMAT
		 );

    writeSubaruMaps(file, filenum, hdf5type, kind, nx, xx, ny, yy, nv, vv, map, vmap);
  }


  free(body);

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

  free(xi);  free(eta);  free(dist);
  free(vxi);  free(veta);  free(vlos);

  free(bodyHead);  free(bodyNum);

  free(map);  free(vmap);
  free(xx);  free(yy);  free(vv);


  exitMPI();

  return (0);
}


/* 3 sigma means that neglecting component of 0.26% */
/* 5 sigma means that neglecting component of 6e-7 */
/* #define SPREAD (3) */
#define SPREAD (5)

#define SMOOTHING_FOR_VISUALIZATION TWO
/* #define SMOOTHING_FOR_VISUALIZATION THREE */

void read_model
(const bool without_subhalo, const uint filenum,
 const int modelID, nbody_aos *body, const int unit, const real eps,
 const ulong Ntot, const int kind, int * restrict bodyHead, int * restrict bodyNum,
 real rot[][3], real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nv, real * restrict vv,
 real * restrict map, real * restrict vmap,
 rand_state *rand
#ifdef  USE_HDF5_FORMAT
 , nbody_hdf5 hdf5, const hdf5struct hdf5type
#endif//USE_HDF5_FORMAT
 )
{
  __NOTE__("%s\n", "start");

  /** read snapshot */
  double time;
  ulong steps;
  int unit_read;
  const ulong Ntot_sim = without_subhalo ? Ntot : (Ntot + 1);
#ifdef  USE_HDF5_FORMAT
  readSnapshot(&unit_read, &time, &steps, Ntot_sim, modelTag[modelID], filenum, &hdf5, hdf5type);
  for(int ii = 0; ii < (int)Ntot_sim; ii++){
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
  readSnapshot(&unit_read, &time, &steps, Ntot_sim, modelTag[modelID], ibody, filenum);
  for(int ii = 0; ii < (int)Ntot_sim; ii++){
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
  __NOTE__("%s is read\n", modelTag[modelID]);


  /** sort particle data by index */
  qsort(body, Ntot_sim, sizeof(nbody_aos), idxAscendingOrder);

  /** obtain surface density maps for visualization */
  standard_coordinate(Ntot, body, rot, xi, eta, dist, vxi, veta, vlos);


  /** generate list of particle index in each field */
  int *list;  list = (int *)malloc(sizeof(int) * Ntot         );  if( list == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate list");  }
  int * num;  num  = (int *)malloc(sizeof(int) * NFIELD * kind);  if(  num == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate num");  }
  int *head;  head = (int *)malloc(sizeof(int) * NFIELD * kind);  if( head == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate head");  }
  for(int ff = 0; ff < NFIELD * kind; ff++){
    num [ff] = 0;
    head[ff] = 0;
  }
  for(int ii = 0; ii < (int)Ntot; ii++)
    list[ii] = (int)Ntot;

  int mem = 0;
  for(int ff = 0; ff < NFIELD; ff++)
    for(int kk = 0; kk < kind; kk++){
      for(int ii = bodyHead[kk]; ii < bodyHead[kk] + bodyNum[kk]; ii++){
	const real dx =  xi[ii] -  xi0[ff];
	const real dy = eta[ii] - eta0[ff];

	if( dx * dx + dy * dy <= hsc_fov2 ){
	  list[mem] = ii;
	  mem++;
	}
      }
      num[INDEX2D(NFIELD, kind, ff, kk)] = mem - head[INDEX2D(NFIELD, kind, ff, kk)];
      if( INDEX2D(NFIELD, kind, ff, kk) + 1 < NFIELD * kind )
	head[INDEX2D(NFIELD, kind, ff, kk) + 1] = mem;
    }


  /** initialize focusing regions */
  if( without_subhalo ){
    for(int ff = 0; ff < NFIELD; ff++){
      /** vmin and vmax will be determined by snapshot without DM sub-halo */
      real vmin =  REAL_MAX;
      real vmax = -REAL_MAX;

      for(int ii = head[INDEX2D(NFIELD, kind, ff, 0)]; ii < head[INDEX2D(NFIELD, kind, ff, kind - 1)] + num[INDEX2D(NFIELD, kind, ff, kind - 1)]; ii++){
	const real vel = vlos[list[ii]];

	vmin = FMIN(vmin, vel);
	vmax = FMAX(vmax, vel);
      }

      /** set grids for the line-of-sight velocity in each fields */
      /** set xx, yy, and vv */
      const real xmin =  xi0[ff] - hsc_fov;      const real xmax =  xi0[ff] + hsc_fov;      const real dx = (xmax - xmin) / (real)nx;
      const real ymin = eta0[ff] - hsc_fov;      const real ymax = eta0[ff] + hsc_fov;      const real dy = (ymax - ymin) / (real)ny;

      for(int ii = 0; ii < nx + 1; ii++)	xx[INDEX2D(NFIELD, nx + 1, ff, ii)] = xmin + dx * (real)ii;
      for(int ii = 0; ii < ny + 1; ii++)	yy[INDEX2D(NFIELD, ny + 1, ff, ii)] = ymin + dy * (real)ii;

      const real dv = (vmax - vmin) / (real)nv;
      for(int ii = 0; ii < nv + 1; ii++)
	vv[INDEX2D(NFIELD, nv + 1, ff, ii)] = vmin + dv * (real)ii;
    }
  }
  __NOTE__("%s is initialized\n", modelTag[modelID]);


  /** generate maps in each field */
  const real deg2kpc = zm31 * CAST_D2R(tan(1.0 * M_PI / 180.0));
  const real kpc2deg = UNITY / deg2kpc;
  for(int ff = 0; ff < NFIELD; ff++){
    const real xmin = xx[INDEX2D(NFIELD, nx + 1, ff, 0)];
    const real ymin = yy[INDEX2D(NFIELD, ny + 1, ff, 0)];
    const real vmin = vv[INDEX2D(NFIELD, nv + 1, ff, 0)];

    const real xmax = xx[INDEX2D(NFIELD, nx + 1, ff, nx)];
    const real ymax = yy[INDEX2D(NFIELD, ny + 1, ff, ny)];
    const real vmax = vv[INDEX2D(NFIELD, nv + 1, ff, nv)];

    const real dx = (xmax - xmin) / (real)nx;
    const real dy = (ymax - ymin) / (real)ny;
    const real dv = (vmax - vmin) / (real)nv;

    const real dxinv = UNITY / dx;
    const real dyinv = UNITY / dy;
    const real dvinv = UNITY / dv;

    const real xsig = FMAX(HALF * eps * kpc2deg, SMOOTHING_FOR_VISUALIZATION * FABS(dx));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
    const real ysig = FMAX(HALF * eps * kpc2deg, SMOOTHING_FOR_VISUALIZATION * FABS(dy));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
    const real vsig = FMAX(HALF                , SMOOTHING_FOR_VISUALIZATION * FABS(dv));/**< sigma in velocity space is set to be 1 km/s (velocity precision is 5-10 km/s in Chiba et al. 2016: http://dx.doi.org/10.1017/S1743921315006912) */

    const real invxsig = UNITY / xsig;
    const real invysig = UNITY / ysig;
    const real invvsig = UNITY / vsig;

    const int nx_smooth = (int)CEIL(SPREAD * FABS(dxinv) * xsig);
    const int ny_smooth = (int)CEIL(SPREAD * FABS(dyinv) * ysig);
    const int nv_smooth = (int)CEIL(SPREAD * FABS(dvinv) * vsig);

    real *erfx;  erfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 2));  if( erfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfx");  }
    real *erfy;  erfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 2));  if( erfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfy");  }
    real *erfv;  erfv = (real *)malloc(sizeof(real) * (2 * nv_smooth + 2));  if( erfv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfv");  }
    for(int ii = 0; ii < 2 * nx_smooth + 2; ii++)    erfx[ii] = ERF((dx * ((real)(ii - nx_smooth) - HALF)) * invxsig);
    for(int ii = 0; ii < 2 * ny_smooth + 2; ii++)    erfy[ii] = ERF((dy * ((real)(ii - ny_smooth) - HALF)) * invysig);
    for(int ii = 0; ii < 2 * nv_smooth + 2; ii++)    erfv[ii] = ERF((dv * ((real)(ii - nv_smooth) - HALF)) * invvsig);

    real *psfx;  psfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 1));  if( psfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfx");  }
    real *psfy;  psfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 1));  if( psfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfy");  }
    real *psfv;  psfv = (real *)malloc(sizeof(real) * (2 * nv_smooth + 1));  if( psfv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfv");  }
    for(int ii = 0; ii < 2 * nx_smooth + 1; ii++)    psfx[ii] = HALF * (erfx[ii + 1] - erfx[ii]);
    for(int ii = 0; ii < 2 * ny_smooth + 1; ii++)    psfy[ii] = HALF * (erfy[ii + 1] - erfy[ii]);
    for(int ii = 0; ii < 2 * nv_smooth + 1; ii++)    psfv[ii] = HALF * (erfv[ii + 1] - erfv[ii]);

    const real dSxyinv = dxinv * dyinv;
    const real dfyvinv = dyinv * dvinv;

    for(int kk = 0; kk < kind; kk++){
      __NOTE__("%s: map initialization: %d/%d for %s\n", field_name[ff], kk, kind, modelTag[modelID]);
      for(int ii = 0; ii < nx * ny; ii++)
	map [INDEX(NFIELD, kind * (NMODEL + NNOISE), nx * ny, ff, INDEX2D(kind, NMODEL + NNOISE, kk, modelID), ii)] = ZERO;
      for(int ii = 0; ii < ny * nv; ii++)
	vmap[INDEX(NFIELD, kind * (NMODEL + NNOISE), ny * nv, ff, INDEX2D(kind, NMODEL + NNOISE, kk, modelID), ii)] = ZERO;
      real mass = ZERO;

      for(int ii = head[INDEX2D(NFIELD, kind, ff, kk)]; ii < head[INDEX2D(NFIELD, kind, ff, kk)] + num[INDEX2D(NFIELD, kind, ff, kk)]; ii++){
	const int idx = list[ii];
	const real mi = CAST_D2R(CAST_R2D(body[idx].m) * mass2astro);
	const real xp =   xi[idx];
	const real yp =  eta[idx];
	const real vp = vlos[idx];
	mass += mi;

	const int l0 = (int)FLOOR((xp - xmin) * dxinv);
	const int m0 = (int)FLOOR((yp - ymin) * dyinv);
	const int n0 = (int)FLOOR((vp - vmin) * dvinv);

	for(int sy = 0; sy < 2 * ny_smooth + 1; sy++){
	  const int mm = m0 + sy - ny_smooth;
	  if( (mm >= 0) && (mm < ny) ){
	    const real my = mi * psfy[sy];

	    for(int sx = 0; sx < 2 * nx_smooth + 1; sx++){
	      const int ll = l0 + sx - nx_smooth;
	      if( (ll >= 0) && (ll < nx) )
		map[INDEX(NFIELD * kind * (NMODEL + NNOISE), nx, ny, INDEX(NFIELD, kind, NMODEL + NNOISE, ff, kk, modelID), ll, mm)] += my * psfx[sx];

	    }

	    for(int sv = 0; sv < 2 * nv_smooth + 1; sv++){
	      const int nn = n0 + sv - nv_smooth;
	      if( (nn >= 0) && (nn < nv) )
		vmap[INDEX(NFIELD * kind * (NMODEL + NNOISE), ny, nv, INDEX(NFIELD, kind, NMODEL + NNOISE, ff, kk, modelID), mm, nn)] += my * psfv[sv];
	    }
	  }
	}

      }


      for(int ii = 0; ii < nx * ny; ii++)
	map [INDEX(NFIELD, kind * (NMODEL + NNOISE), nx * ny, ff, INDEX2D(kind, NMODEL + NNOISE, kk, modelID), ii)] *= dSxyinv;
      for(int ii = 0; ii < ny * nv; ii++)
	vmap[INDEX(NFIELD, kind * (NMODEL + NNOISE), ny * nv, ff, INDEX2D(kind, NMODEL + NNOISE, kk, modelID), ii)] *= dfyvinv;
      __NOTE__("%s: map generated: %d/%d for %s\n", field_name[ff], kk, kind, modelTag[modelID]);


      /** generate noise map */
      if( without_subhalo ){
	/* field: ff, kind: kk */
	for(int sn = 0; sn < NNOISE; sn++){
	  for(int ii = 0; ii < nx * ny; ii++)
	    map [INDEX(NFIELD, kind * (NMODEL + NNOISE), nx * ny, ff, INDEX2D(kind, NMODEL + NNOISE, kk, NMODEL + sn), ii)] = ZERO;
	  for(int ii = 0; ii < ny * nv; ii++)
	    vmap[INDEX(NFIELD, kind * (NMODEL + NNOISE), ny * nv, ff, INDEX2D(kind, NMODEL + NNOISE, kk, NMODEL + sn), ii)] = ZERO;
	}
	__NOTE__("%s: noise map initialized: %d/%d for %s\n", field_name[ff], kk, kind, modelTag[modelID]);

	const real m_noise = mass / (real)NDIST;


	for(int ii = 0; ii < NDIST; ii++){
	  real xn, yn;
	  while( true ){
	    xn = xmin + (xmax - xmin) * UNIRAND(rand);
	    yn = ymin + (ymax - ymin) * UNIRAND(rand);

	    const real dx = xn -  xi0[ff];
	    const real dy = yn - eta0[ff];
	    if( dx * dx + dy * dy <= hsc_fov2 )
	      break;
	  }
	  const real vn = vmin + (vmax - vmin) * UNIRAND(rand);
	  __NOTE__("%d/%d noise particle\n", ii, NDIST);

	  const int l0 = (int)FLOOR((xn - xmin) * dxinv);
	  const int m0 = (int)FLOOR((yn - ymin) * dyinv);
	  const int n0 = (int)FLOOR((vn - vmin) * dvinv);

	  for(int sy = 0; sy < 2 * ny_smooth + 1; sy++){
	    const int mm = m0 + sy - ny_smooth;
	    if( (mm >= 0) && (mm < ny) ){
	      const real my = m_noise * psfy[sy];

	      for(int sx = 0; sx < 2 * nx_smooth + 1; sx++){
		const int ll = l0 + sx - nx_smooth;
		if( (ll >= 0) && (ll < nx) ){
		  const real mset = my * psfx[sx];

		  for(int sn = 0; sn < NNOISE; sn++)
		    map[INDEX(NFIELD * kind * (NMODEL + NNOISE), nx, ny, INDEX(NFIELD, kind, NMODEL + NNOISE, ff, kk, NMODEL + sn), ll, mm)] += mset * sn_inv[sn];
		}

	      }

	      for(int sv = 0; sv < 2 * nv_smooth + 1; sv++){
		const int nn = n0 + sv - nv_smooth;
		if( (nn >= 0) && (nn < nv) ){
		  const real mset = my * psfv[sv];

		  for(int sn = 0; sn < NNOISE; sn++)
		    vmap[INDEX(NFIELD * kind * (NMODEL + NNOISE), ny, nv, INDEX(NFIELD, kind, NMODEL + NNOISE, ff, kk, NMODEL + sn), mm, nn)] += mset * sn_inv[sn];
		}
	      }
	    }
	  }
	}

	for(int sn = 0; sn < NNOISE; sn++){
	  for(int ii = 0; ii < nx * ny; ii++)
	    map [INDEX(NFIELD, kind * (NMODEL + NNOISE), nx * ny, ff, INDEX2D(kind, NMODEL + NNOISE, kk, NMODEL + sn), ii)] *= dSxyinv;
	  for(int ii = 0; ii < ny * nv; ii++)
	    vmap[INDEX(NFIELD, kind * (NMODEL + NNOISE), ny * nv, ff, INDEX2D(kind, NMODEL + NNOISE, kk, NMODEL + sn), ii)] *= dfyvinv;
	}
	__NOTE__("%s: noise map generated: %d/%d for %s\n", field_name[ff], kk, kind, modelTag[modelID]);
      }
    }

    free(erfx);    free(erfy);    free(erfv);
    free(psfx);    free(psfy);    free(psfv);
  }

  free(list);
  free(num);
  free(head);


  __NOTE__("%s\n", "end");
}


#ifdef  USE_HDF5_FORMAT
/**
 * @fn writeSubaruMaps
 *
 * @brief Write analyzed profiles of the N-body simulation.
 */
void writeSubaruMaps
(char file[], const uint id, hdf5struct type, const int kind,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nv, real * restrict vv,
 real * restrict map, real * restrict vmap)
{
  __NOTE__("%s\n", "start");

  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "subaru", id);
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
  attribute = H5Acreate(target, "nv", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nv));
  chkHDF5err(H5Aclose(attribute));
  /* write flag about USE_DOUBLE_PRECISION */
#ifdef  USE_DOUBLE_PRECISION
  const int useDP = 1;
#else///USE_DOUBLE_PRECISION
  const int useDP = 0;
#endif//USE_DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));


  hid_t group = H5Gcreate(target, "subaru", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(group, "FoV", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &hsc_fov));
  chkHDF5err(H5Aclose(attribute));
  int tmp = NFIELD;
  attribute = H5Acreate(group, "Nfield", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &tmp));
  chkHDF5err(H5Aclose(attribute));
  tmp = NMODEL;
  attribute = H5Acreate(group, "Nmodel", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &tmp));
  chkHDF5err(H5Aclose(attribute));
  tmp = NNOISE;
  attribute = H5Acreate(group, "Nnoise", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &tmp));
  chkHDF5err(H5Aclose(attribute));
  tmp = NDIST;
  attribute = H5Acreate(group, "Ndist", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &tmp));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));

  hid_t H5T_VARIABLE_STRING = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(H5T_VARIABLE_STRING, H5T_VARIABLE));

  attr_dims = NFIELD;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(group, "xi", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, xi0));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "eta", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, eta0));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "field", H5T_VARIABLE_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_VARIABLE_STRING, field_name));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));

  attr_dims = NMODEL;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(group, "model", H5T_VARIABLE_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_VARIABLE_STRING, modelTag));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));

  chkHDF5err(H5Tclose(H5T_VARIABLE_STRING));

  attr_dims = NNOISE;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(group, "SN", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, sn_val));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));

  chkHDF5err(H5Gclose(group));


  for(int ii = 0; ii < kind; ii++){
    char grp[16];
    sprintf(grp, "component%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /** 2D (nx * ny) array */
    hsize_t dims[2] = {nx, ny};
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
    for(int ff = 0; ff < NFIELD; ff++){
      char dat[64];
      for(int model = 0; model < NMODEL; model++){
	sprintf(dat, "%s-%s-map", field_name[ff], modelTag[model]);
	dataset = H5Dcreate(group, dat, type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
	chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &map[INDEX(NFIELD, kind * (NMODEL + NNOISE), nx * ny, ff, INDEX2D(kind, NMODEL + NNOISE, ii, model), 0)]));
	chkHDF5err(H5Dclose(dataset));
      }
      for(int noise = 0; noise < NNOISE; noise++){
	sprintf(dat, "%s-SN%02d-map", field_name[ff], sn_val[noise]);
	dataset = H5Dcreate(group, dat, type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
	chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &map[INDEX(NFIELD, kind * (NMODEL + NNOISE), nx * ny, ff, INDEX2D(kind, NMODEL + NNOISE, ii, NMODEL + noise), 0)]));
	chkHDF5err(H5Dclose(dataset));
      }
    }
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
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
    for(int ff = 0; ff < NFIELD; ff++){
      char dat[64];
      for(int model = 0; model < NMODEL; model++){
	sprintf(dat, "%s-%s-vel", field_name[ff], modelTag[model]);
	dataset = H5Dcreate(group, dat, type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
	chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vmap[INDEX(NFIELD, kind * (NMODEL + NNOISE), ny * nv, ff, INDEX2D(kind, NMODEL + NNOISE, ii, model), 0)]));
	chkHDF5err(H5Dclose(dataset));
      }
      for(int noise = 0; noise < NNOISE; noise++){
	sprintf(dat, "%s-SN%02d-vel", field_name[ff], sn_val[noise]);
	dataset = H5Dcreate(group, dat, type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
	chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vmap[INDEX(NFIELD, kind * (NMODEL + NNOISE), ny * nv, ff, INDEX2D(kind, NMODEL + NNOISE, ii, NMODEL + noise), 0)]));
	chkHDF5err(H5Dclose(dataset));
      }
    }
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));


    /** 1D (nx + 1) array */
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
    for(int ff = 0; ff < NFIELD; ff++){
      char dat[32];
      sprintf(dat, "%s-xi", field_name[ff]);
      dataset = H5Dcreate(group, dat, type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &xx[INDEX2D(NFIELD, nx + 1, ff, 0)]));
      chkHDF5err(H5Dclose(dataset));
    }
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (ny + 1) array */
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
    for(int ff = 0; ff < NFIELD; ff++){
      char dat[32];
      sprintf(dat, "%s-eta", field_name[ff]);
      dataset = H5Dcreate(group, dat, type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yy[INDEX2D(NFIELD, nx + 1, ff, 0)]));
      chkHDF5err(H5Dclose(dataset));
    }
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nv + 1) array */
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
    for(int ff = 0; ff < NFIELD; ff++){
      char dat[32];
      sprintf(dat, "%s-vlos", field_name[ff]);
      dataset = H5Dcreate(group, dat, type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vv[INDEX2D(NFIELD, nv + 1, ff, 0)]));
      chkHDF5err(H5Dclose(dataset));
    }
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));


    chkHDF5err(H5Gclose(group));
  }





  /* close the file */
  chkHDF5err(H5Fclose(target));




  __NOTE__("%s\n", "end");
}
#endif//USE_HDF5_FORMAT
