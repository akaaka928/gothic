/**
 * @file m31obs.c
 *
 * @brief Generate surface density maps in the observed coordinate (for N-body runs in M31's disk coordinate)
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

#ifdef  USE_HDF5_FORMAT
#ifdef  HDF5_FOR_ZINDAIJI
#include <unistd.h>
#include <sys/stat.h>
#endif//HDF5_FOR_ZINDAIJI
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



const real zm31 = CAST_D2R(776.0);/**< same with Komiyama et al. (2018); median of NED */
/* const real zm31 = CAST_D2R(773.0);/\**< M31 distance in Conn et al. (2016) *\/ */

const real vm31x = CAST_D2R(0.0);
const real vm31y = CAST_D2R(0.0);
const real vm31z = CAST_D2R(-300.0);


extern const double      length2astro;
extern const double        time2astro;
extern const double        mass2astro;
extern const double    velocity2astro;


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
} nbody_particle;


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
#ifdef __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC



void setRotationMatrix(real rot[restrict][3], real inv[restrict][3]);
void standard_coordinate(const int num, nbody_particle *body, real rot[restrict][3], real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos);

void generateMassDistributionMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body, real rot[restrict][3], const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz, real * restrict rho_map);
void generateSurfaceDensityMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body, real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vlos, const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx, const int nv, const real vmin, const real dv, real * restrict f_xv, real * restrict f_yv, real * restrict f_zv);

#ifdef  USE_HDF5_FORMAT
void writeM31coordinateData
(const double time, const ulong steps, char file[], const uint id, hdf5struct type, const int kind, int * restrict bodyHead, int * restrict bodyNum,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nz, real * restrict zz, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx,
 const int nv, real * restrict vl, real * restrict f_xv, real * restrict f_yv, real * restrict f_zv,
 const int nx3D, real * restrict rho_xx, const int ny3D, real * restrict rho_yy, const int nz3D, real * restrict rho_zz, real * restrict rho_map,
 real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos
#ifdef  HDF5_FOR_ZINDAIJI
 , real rot[restrict][3], nbody_hdf5 hdf5, nbody_particle *body, int * restrict bodyType, const real eps
#endif//HDF5_FOR_ZINDAIJI
);

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
  if( argc < 20 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 20);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -nxi=<int> -ximin=<real> -ximax=<real>\n");
    __FPRINTF__(stderr, "          -neta=<int> -etamin=<real> -etamax=<real>\n");
    __FPRINTF__(stderr, "          -nD=<int> -Dmin=<real> -Dmax=<real>\n");
    __FPRINTF__(stderr, "          -nvlos=<int> -vlosmin=<real> -vlosmax=<real>\n");
    __FPRINTF__(stderr, "          -nx3D=<int> -ny3D=<int> -nz3D=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 20 ){ */

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));

  int nx;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nxi" , &nx));
  int ny;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "neta", &ny));
  int nz;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nD"  , &nz));
  int nx3D;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nx3D", &nx3D));
  int ny3D;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "ny3D", &ny3D));
  int nz3D;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nz3D", &nz3D));
  real xmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,  "ximax", &xmin));/**< x = -xi */
  real ymin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "etamin", &ymin));/**< y = eta */
  real zmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,   "Dmin", &zmin));/**< z = D */
  real xmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,  "ximin", &xmax));/**< x = -xi */
  real ymax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "etamax", &ymax));/**< y = eta */
  real zmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,   "Dmax", &zmax));/**< z = D */

  int nvlos;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nvlos", &nvlos));
  real vlosmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vlosmin", &vlosmin));
  real vlosmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vlosmax", &vlosmax));

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
  }
  bodyHead[0] = 0;
  for(int ii = 1; ii < kind; ii++)
    bodyHead[ii] = bodyHead[ii - 1] + bodyNum[ii - 1];
#ifdef  HDF5_FOR_ZINDAIJI
  for(int ii = 0; ii < kind; ii++)
    bodyType[ii] &= 3;
#endif//HDF5_FOR_ZINDAIJI

  /** assume zone-centered mapping */
  real *rho_xx ;  rho_xx  = (real *)malloc(sizeof(real) * (nx3D + 1));  if( rho_xx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rho_xx");  }
  real *rho_yy ;  rho_yy  = (real *)malloc(sizeof(real) * (ny3D + 1));  if( rho_yy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rho_yy");  }
  real *rho_zz ;  rho_zz  = (real *)malloc(sizeof(real) * (nz3D + 1));  if( rho_zz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rho_zz");  }
  real *rho_map;  rho_map = (real *)malloc(sizeof(real) * kind * nx3D * ny3D * nz3D);  if( rho_map == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rho_map");  }
  const real deg2kpc = zm31 * CAST_D2R(tan(1.0 * M_PI / 180.0));
  const real rho_xmin = deg2kpc * xmax;  const real rho_xmax = deg2kpc * xmin;
  const real rho_ymin = deg2kpc * ymin;  const real rho_ymax = deg2kpc * ymax;
  const real rho_zmin =           zmin;  const real rho_zmax =           zmax;
  const real rho_dx = (rho_xmax - rho_xmin) / (real)nx3D;  for(int ii = 0; ii < nx3D + 1; ii++)    rho_xx[ii] = rho_xmin + rho_dx * (real)ii;
  const real rho_dy = (rho_ymax - rho_ymin) / (real)ny3D;  for(int jj = 0; jj < ny3D + 1; jj++)    rho_yy[jj] = rho_ymin + rho_dy * (real)jj;
  const real rho_dz = (rho_zmax - rho_zmin) / (real)nz3D;  for(int kk = 0; kk < nz3D + 1; kk++)    rho_zz[kk] = rho_zmin + rho_dz * (real)kk;

  real *xx;  xx = (real *)malloc(sizeof(real) * (nx + 1));  if( xx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xx");  }
  real *yy;  yy = (real *)malloc(sizeof(real) * (ny + 1));  if( yy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate yy");  }
  real *zz;  zz = (real *)malloc(sizeof(real) * (nz + 1));  if( zz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate zz");  }
  zmin += zm31;
  zmax += zm31;
  const real dx = (xmax - xmin) / (real)nx;  for(int ii = 0; ii < nx + 1; ii++)    xx[ii] = xmin + dx * (real)ii;
  const real dy = (ymax - ymin) / (real)ny;  for(int jj = 0; jj < ny + 1; jj++)    yy[jj] = ymin + dy * (real)jj;
  const real dz = (zmax - zmin) / (real)nz;  for(int kk = 0; kk < nz + 1; kk++)    zz[kk] = zmin + dz * (real)kk;
  real *Sigma_xy;  Sigma_xy = (real *)malloc(sizeof(real) * kind * nx * ny);  if( Sigma_xy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Sigma_xy");  }
  real *Sigma_yz;  Sigma_yz = (real *)malloc(sizeof(real) * kind * ny * nz);  if( Sigma_yz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Sigma_yz");  }
  real *Sigma_zx;  Sigma_zx = (real *)malloc(sizeof(real) * kind * nz * nx);  if( Sigma_zx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Sigma_zx");  }

  real *vl;  vl = (real *)malloc(sizeof(real) * (nvlos + 1));  if( vl == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vl");  }
  const real dv = (vlosmax - vlosmin) / (real)nvlos;  for(int ii = 0; ii < nvlos + 1; ii++)    vl[ii] = vlosmin + dv * (real)ii;
  real *f_xv;  f_xv = (real *)malloc(sizeof(real) * kind * nx * nvlos);  if( f_xv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_xv");  }
  real *f_yv;  f_yv = (real *)malloc(sizeof(real) * kind * ny * nvlos);  if( f_yv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_yv");  }
  real *f_zv;  f_zv = (real *)malloc(sizeof(real) * kind * nz * nvlos);  if( f_zv == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_zv");  }


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

    /** obtain mass distribution map for volume rendering */
    generateMassDistributionMaps(kind, bodyHead, bodyNum, body, rot, eps, nx3D, rho_xmin, rho_dx, ny3D, rho_ymin, rho_dy, nz3D, rho_zmin, rho_dz, rho_map);

    /** obtain surface density maps for visualization */
    standard_coordinate(Ntot, body, rot, xi, eta, dist, vxi, veta, vlos);
    generateSurfaceDensityMaps(kind, bodyHead, bodyNum, body, xi, eta, dist, vlos, eps, nx, xmin, dx, ny, ymin, dy, nz, zmin, dz, Sigma_xy, Sigma_yz, Sigma_zx, nvlos, vlosmin, dv, f_xv, f_yv, f_zv);


    /** dump analyzed results for matplotlib and/or VisIt */
#ifdef  USE_HDF5_FORMAT
    writeM31coordinateData
      (time, steps, file, filenum, hdf5type, kind, bodyHead, bodyNum,
       nx, xx, ny, yy, nz, zz, Sigma_xy, Sigma_yz, Sigma_zx,
       nvlos, vl, f_xv, f_yv, f_zv,
       nx3D, rho_xx, ny3D, rho_yy, nz3D, rho_zz, rho_map,
       xi, eta, dist, vxi, veta, vlos
#ifdef  HDF5_FOR_ZINDAIJI
       , rot, hdf5, body, bodyType, eps
#endif//HDF5_FOR_ZINDAIJI
       );
#endif//USE_HDF5_FORMAT
  }/* for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){ */


  if( mpi.rank == 0 ){
    fprintf(stdout, "assumed value in the analysis:\n");
    fprintf(stdout, "distance to M31 is %e kpc\n", zm31);
    fprintf(stdout, "vM31_x = %e km/s\n", vm31x);
    fprintf(stdout, "vM31_y = %e km/s\n", vm31y);
    fprintf(stdout, "vM31_z = %e km/s\n", vm31z);
  }/* if( mpi.rank == 0 ){ */


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

  free(xi);  free(eta);  free(dist);
  free(vxi);  free(veta);  free(vlos);

  free(xx);  free(yy);  free(zz);
  free(Sigma_xy);  free(Sigma_yz);  free(Sigma_zx);

  free(vl);
  free(f_xv);  free(f_yv);  free(f_zv);

  free(rho_xx);  free(rho_yy);  free(rho_zz);
  free(rho_map);

  exitMPI();

  return (0);
}


/**
 * @fn scalarProd
 *
 * @brief Calculate scalar product.
 *
 * @param (aa) vector A
 * @param (bb) vector B
 * @return scalar product of A and B
 */
static inline real scalarProd(real aa[], real bb[])
{
  return (aa[0] * bb[0] + aa[1] * bb[1] + aa[2] * bb[2]);
}
/**
 * @fn vectorProd
 *
 * @brief Calculate vector product.
 *
 * @param (aa) vector A
 * @param (bb) vector B
 * @return (ans) vector product of A and B
 */
static inline void vectorProd(real aa[], real bb[], real ans[])
{
  ans[0] = aa[1] * bb[2] - aa[2] * bb[1];
  ans[1] = aa[2] * bb[0] - aa[0] * bb[2];
  ans[2] = aa[0] * bb[1] - aa[1] * bb[0];
}

static inline void rgemm(real aa[restrict][3], real bb[restrict][3], real ans[restrict][3])
{
  for(int ii = 0; ii < 3; ii++)
    for(int jj = 0; jj < 3; jj++)
      ans[ii][jj] = ZERO;

  for(int ii = 0; ii < 3; ii++)
    for(int kk = 0; kk < 3; kk++)
      for(int jj = 0; jj < 3; jj++)
	ans[ii][jj] += aa[ii][kk] * bb[kk][jj];
}

static inline void transpose(real ini[restrict][3], real fin[restrict][3])
{
  for(int ii = 0; ii < 3; ii++)
    for(int jj = 0; jj < 3; jj++)
      fin[jj][ii] = ini[ii][jj];
}


void setRotationMatrix(real rot[restrict][3], real inv[restrict][3])
{
  __NOTE__("%s\n", "start");

  const real theta = -37.0 * CAST_D2R(M_PI / 180.0);/**< position angle */
  const real   phi = -77.0 * CAST_D2R(M_PI / 180.0);/**< inclination */

  const real cost = COS(theta);  const real cosp = COS(phi);
  const real sint = SIN(theta);  const real sinp = SIN(phi);

  /* 1st rotation: rotation axis is z-axis, rotation angle is theta */
  static real rot1[3][3], inv1[3][3];
  static real axis1[3] = {ZERO, ZERO, UNITY};
  setRodriguesRotationMatrix(axis1, sint, cost, rot1, inv1);

  /* 2nd rotation: rotation axis is rotated y-axis (rot1 * (0, 1, 0)), rotation angle is phi */
  static real rot2[3][3], inv2[3][3];
  static real axis2[3];
  axis1[0] = ZERO;
  axis1[1] = UNITY;
  axis1[2] = ZERO;
  rotateVector(axis1, rot1, axis2);
  setRodriguesRotationMatrix(axis2, sinp, cosp, rot2, inv2);

  /* get rotation axis of M31's disk in observed frame */
  static real axis3[3], spin[3];
  spin[0] = ZERO;
  spin[1] = ZERO;
  spin[2] = -UNITY;
  rotateVector(spin, rot1, axis3);
  rotateVector(axis3, rot2, spin);

  /* rotate spin axis of M31's disk */
  axis3[0] = ZERO;
  axis3[1] = ZERO;
  axis3[2] = UNITY;
  static real rot3[3][3], inv3[3][3];
  initRotationMatrices(spin, axis3, rot3, inv3);

  /* major axis in M31's disk frame */
  static real major_disk[3] = {ZERO, UNITY, ZERO};
  static real major_tmp[3];
  rotateVector(major_disk, rot1, major_tmp);
  static real major_obs[3];
  rotateVector(major_tmp, rot2, major_obs);
  rotateVector(major_disk, inv3, major_tmp);
  static real rot4[3][3], inv4[3][3];
  static real axis[3];
  vectorProd(major_tmp, major_obs, axis);
  const real sintheta = -SQRT(scalarProd(axis, axis));
  const real costheta = scalarProd(major_tmp, major_obs);
  setRodriguesRotationMatrix(spin, sintheta, costheta, rot4, inv4);

  rgemm(rot4, inv3, inv);
  transpose(inv, rot);

  __NOTE__("%s\n", "end");
}


void standard_coordinate(const int num, nbody_particle *body, real rot[restrict][3], real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos)
{
  __NOTE__("%s\n", "start");

  const real rad2deg = CAST_D2R(180.0 * M_1_PI);
  static real ini[3], fin[3];

  for(int ii = 0; ii < num; ii++){
    /* coordinate rotation of particle position */
    ini[0] = CAST_D2R(CAST_R2D(body[ii].x) * length2astro);
    ini[1] = CAST_D2R(CAST_R2D(body[ii].y) * length2astro);
    ini[2] = CAST_D2R(CAST_R2D(body[ii].z) * length2astro);
    rotateVector(ini, rot, fin);
    const real xx = fin[0];
    const real yy = fin[1];
    const real zz = fin[2] + zm31;

    /* coordinate rotation of particle velocity */
    ini[0] = CAST_D2R(CAST_R2D(body[ii].vx) * velocity2astro);
    ini[1] = CAST_D2R(CAST_R2D(body[ii].vy) * velocity2astro);
    ini[2] = CAST_D2R(CAST_R2D(body[ii].vz) * velocity2astro);
    rotateVector(ini, rot, fin);
    const real vx = fin[0] + vm31x;
    const real vy = fin[1] + vm31y;
    const real vz = fin[2] + vm31z;

    const real tmp = UNITY / zz;
    xi [ii] = rad2deg * ATAN(xx * tmp);
    eta[ii] = rad2deg * ATAN(yy * tmp);
    const real d2 = xx * xx + yy * yy + zz * zz;
    const real dinv = RSQRT(d2);
    dist[ii] = d2 * dinv;

    vxi [ii] = (zz * vx - xx * vz) * dinv;
    veta[ii] = (zz * vy - yy * vz) * dinv;
    vlos[ii] = (xx * vx + yy * vy + zz * vz) * dinv;
  }/* for(int ii = 0; ii < num; ii++){ */

  __NOTE__("%s\n", "end");
}


/* 3 sigma means that neglecting component of 0.26% */
/* 5 sigma means that neglecting component of 6e-7 */
/* #define SPREAD (3) */
#define SPREAD (5)
void generateMassDistributionMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body, real rot[restrict][3], const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz, real * restrict rho_map)
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
  for(int ii = 0; ii < 2 * ny_smooth + 2; ii++)	   erfy[ii] = ERF((dy * ((real)(ii - ny_smooth) - HALF)) * invsig);
  for(int ii = 0; ii < 2 * nz_smooth + 2; ii++)	   erfz[ii] = ERF((dz * ((real)(ii - nz_smooth) - HALF)) * invsig);

  real *psfx;  psfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 1));  if( psfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfx");  }
  real *psfy;  psfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 1));  if( psfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfy");  }
  real *psfz;  psfz = (real *)malloc(sizeof(real) * (2 * nz_smooth + 1));  if( psfz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfz");  }
  for(int ii = 0; ii < 2 * nx_smooth + 1; ii++)    psfx[ii] = HALF * (erfx[ii + 1] - erfx[ii]);
  for(int ii = 0; ii < 2 * ny_smooth + 1; ii++)    psfy[ii] = HALF * (erfy[ii + 1] - erfy[ii]);
  for(int ii = 0; ii < 2 * nz_smooth + 1; ii++)    psfz[ii] = HALF * (erfz[ii + 1] - erfz[ii]);

  const real dVinv = dxinv * dyinv * dzinv;

  static real ini[3], fin[3];
  for(int kk = 0; kk < kind; kk++){
    for(int ii = 0; ii < nx * ny * nz; ii++)
      rho_map[INDEX2D(kind, nx * ny * nz, kk, ii)] = ZERO;

    for(int ii = bodyHead[kk]; ii < bodyHead[kk] + bodyNum[kk]; ii++){
      const real mi = CAST_D2R(CAST_R2D(body[ii].m) *   mass2astro);

      /* coordinate rotation of particle position */
      ini[0] = CAST_D2R(CAST_R2D(body[ii].x) * length2astro);
      ini[1] = CAST_D2R(CAST_R2D(body[ii].y) * length2astro);
      ini[2] = CAST_D2R(CAST_R2D(body[ii].z) * length2astro);
      rotateVector(ini, rot, fin);
      const real xx = fin[0];
      const real yy = fin[1];
      const real zz = fin[2];

      const int l0 = (int)FLOOR((xx - xmin) * dxinv);
      const int m0 = (int)FLOOR((yy - ymin) * dyinv);
      const int n0 = (int)FLOOR((zz - zmin) * dzinv);

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
  }/* for(int kk = 0; kk < kind; kk++){ */

  free(erfx);  free(erfy);  free(erfz);
  free(psfx);  free(psfy);  free(psfz);


  __NOTE__("%s\n", "end");
}


#define SMOOTHING_FOR_VISUALIZATION TWO
/* #define SMOOTHING_FOR_VISUALIZATION THREE */
void generateSurfaceDensityMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body, real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vlos, const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx, const int nv, const real vmin, const real dv, real * restrict f_xv, real * restrict f_yv, real * restrict f_zv)
{
  __NOTE__("%s\n", "start");

  const real deg2kpc = zm31 * CAST_D2R(tan(1.0 * M_PI / 180.0));
  const real kpc2deg = UNITY / deg2kpc;

  const real dxinv = UNITY / dx;
  const real dyinv = UNITY / dy;
  const real dzinv = UNITY / dz;
  const real dvinv = UNITY / dv;

  const real xsig = FMAX(HALF * eps * kpc2deg, SMOOTHING_FOR_VISUALIZATION * FABS(dx));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real ysig = FMAX(HALF * eps * kpc2deg, SMOOTHING_FOR_VISUALIZATION * FABS(dy));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real zsig = FMAX(HALF * eps          , SMOOTHING_FOR_VISUALIZATION * FABS(dz));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real vsig =                            SMOOTHING_FOR_VISUALIZATION * FABS(dv) ;

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



#ifdef  USE_DEGREE_FOR_SURFACE_DENSITY_MAP
  /** show Msun deg^-2 */
  /* const real inv_kpc2deg = deg2kpc; */
  const real dSxyinv = dxinv * dyinv;
  /* const real dSyzinv =         dyinv * dzinv * inv_kpc2deg; */
  /* const real dSzxinv = dxinv	     * dzinv * inv_kpc2deg; */
  const real dSyzinv =         dyinv * dzinv;
  const real dSzxinv = dxinv	     * dzinv;

  const real df_xvinv = dxinv               * dvinv;
  const real df_yvinv = dyinv               * dvinv;
  /* const real df_zvinv = dzinv * inv_kpc2deg * dvinv; */
  const real df_zvinv = dzinv * dvinv;
#else///USE_DEGREE_FOR_SURFACE_DENSITY_MAP
  /** show Msun kpc^-2 */
  const real inv_deg2kpc = kpc2deg;
  const real dSxyinv = dxinv * inv_deg2kpc * dyinv * inv_deg2kpc;
  const real dSyzinv =	                     dyinv * inv_deg2kpc * dzinv;
  const real dSzxinv = dxinv * inv_deg2kpc			 * dzinv;

  const real df_xvinv = dxinv * inv_deg2kpc * dvinv;
  const real df_yvinv = dyinv * inv_deg2kpc * dvinv;
  const real df_zvinv = dzinv               * dvinv;
#endif//USE_DEGREE_FOR_SURFACE_DENSITY_MAP
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
      const real xx =   xi[ii];
      const real yy =  eta[ii];
      const real zz = dist[ii];
      const real vv = vlos[ii];

      const int l0 = (int)FLOOR((xx - xmin) * dxinv);
      const int m0 = (int)FLOOR((yy - ymin) * dyinv);
      const int n0 = (int)FLOOR((zz - zmin) * dzinv);
      const int o0 = (int)FLOOR((vv - vmin) * dvinv);


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
	    const int oo = o0 + sv - nv_smooth;
	    if( (oo >= 0) && (oo < nv) )
	      f_yv[INDEX3D(kind, ny, nv, kk, mm, oo)] += my * psfv[sv];
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
	    const int oo = o0 + sv - nv_smooth;
	    if( (oo >= 0) && (oo < nv) )
	      f_zv[INDEX3D(kind, nz, nv, kk, nn, oo)] += mz * psfv[sv];
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
/**
 * @fn writeM31coordinateData
 *
 * @brief Write analyzed profiles of the N-body simulation.
 */
void writeM31coordinateData
(const double time, const ulong steps, char file[], const uint id, hdf5struct type, const int kind, int * restrict bodyHead, int * restrict bodyNum,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nz, real * restrict zz, real * restrict Sigma_xy, real * restrict Sigma_yz, real * restrict Sigma_zx,
 const int nv, real * restrict vl, real * restrict f_xv, real * restrict f_yv, real * restrict f_zv,
 const int nx3D, real * restrict rho_xx, const int ny3D, real * restrict rho_yy, const int nz3D, real * restrict rho_zz, real * restrict rho_map,
 real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos
#ifdef  HDF5_FOR_ZINDAIJI
 , real rot[restrict][3], nbody_hdf5 hdf5, nbody_particle *body, int * restrict bodyType, const real eps
#endif//HDF5_FOR_ZINDAIJI
)
{
  __NOTE__("%s\n", "start");

  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "m31obs", id);
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
  attribute = H5Acreate(target, "nx3D", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nx3D));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "ny3D", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &ny3D));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nz3D", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nz3D));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nv", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nv));
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
  /* write flag about USE_DEGREE_FOR_SURFACE_DENSITY_MAP */
#ifdef  USE_DEGREE_FOR_SURFACE_DENSITY_MAP
  const int useDegree = 1;
#else///USE_DEGREE_FOR_SURFACE_DENSITY_MAP
  const int useDegree = 0;
#endif//USE_DEGREE_FOR_SURFACE_DENSITY_MAP
  attribute = H5Acreate(target, "useDegree", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDegree));
  chkHDF5err(H5Aclose(attribute));
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));


  for(int ii = 0; ii < kind; ii++){
    char grp[16];    sprintf(grp, "field%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t attr_dims = 1;
    hid_t attribute;

    /** 3D (nx3D * ny3D * nz3D) array */
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
#endif//USE_SZIP_COMPRESSION
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
    dataset = H5Dcreate(group, "xi", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, xx));
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
    dataset = H5Dcreate(group, "x", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho_xx));
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
    dataset = H5Dcreate(group, "eta", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, yy));
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
    dataset = H5Dcreate(group, "y", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho_yy));
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
    dataset = H5Dcreate(group, "D", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, zz));
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
    dataset = H5Dcreate(group, "z", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho_zz));
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
    dataset = H5Dcreate(group, "vlos", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));


    sprintf(grp, "obs%d", ii);
    group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    /** 1D (num) arrays */
    dims[0] = bodyNum[ii];
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
    dataset = H5Dcreate(group, "xi", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &xi[bodyHead[ii]]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "eta", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &eta[bodyHead[ii]]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "dist", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dist[bodyHead[ii]]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "vxi", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vxi[bodyHead[ii]]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "veta", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &veta[bodyHead[ii]]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "vlos", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vlos[bodyHead[ii]]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
      /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    /* write attribute data */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    /* write # of particles */
    attribute = H5Acreate(group, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &bodyNum[ii]));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
  }/* for(int ii = 0; ii < kind; ii++){ */

  /* close the file */
  chkHDF5err(H5Fclose(target));


#ifdef  HDF5_FOR_ZINDAIJI
  sprintf(filename, "%s/%s_%s%.3u.h5", DATAFOLDER, file, "zindaiji", id);
  bool dump_file = (0 != access(filename, F_OK));
  if( !dump_file ){
    struct stat stat_file;
    stat(filename, &stat_file);
    char tmpname[128];
    sprintf(tmpname, "%s/%s.%s%.3u.h5", DATAFOLDER, file, SNAPSHOT, id);
    struct stat stat_snap;
    stat(tmpname, &stat_snap);
    if( stat_snap.st_ctime > stat_file.st_ctime )
      dump_file = true;
  }/* if( !dump_file ){ */
  if( dump_file ){
    int Ntot = 0;
    for(int kk = 0; kk < kind; kk++)
      Ntot += bodyNum[kk];
    static real ini[3], fin[3];
    /* execute coordinate rotation */
    const double velocity2com = 1.0 / velocity2astro;
    const double vm31x_com = vm31x * velocity2com;
    const double vm31y_com = vm31y * velocity2com;
    const double vm31z_com = vm31z * velocity2com;
    for(int ii = 0; ii < Ntot; ii++){
      hdf5.idx[ii    ] = body[ii].idx;	hdf5.  m[ii        ] = body[ii]. m;	hdf5.pot[ii        ] = body[ii].pot;
      ini[0] = CAST_D2R(CAST_R2D(body[ii].x) * length2astro);      ini[1] = CAST_D2R(CAST_R2D(body[ii].y) * length2astro);      ini[2] = CAST_D2R(CAST_R2D(body[ii].z) * length2astro);
      rotateVector(ini, rot, fin);
      fin[2] += zm31;
      hdf5.pos[ii * 3] = fin[0];	hdf5.pos[ii * 3 + 1] = fin[1];	hdf5.pos[ii * 3 + 2] = fin[2];
      ini[0] = body[ii].vx;      ini[1] = body[ii].vy;      ini[2] = body[ii].vz;
      fin[0] = CAST_D2R(CAST_R2D(fin[0]) + vm31x_com);
      fin[1] = CAST_D2R(CAST_R2D(fin[1]) + vm31y_com);
      fin[2] = CAST_D2R(CAST_R2D(fin[2]) + vm31z_com);
      hdf5.vel[ii * 3] = fin[0];	hdf5.vel[ii * 3 + 1] = fin[1];	hdf5.vel[ii * 3 + 2] = fin[2];
      ini[0] = body[ii].ax;      ini[1] = body[ii].ay;      ini[2] = body[ii].az;
      hdf5.acc[ii * 3] = fin[0];	hdf5.acc[ii * 3 + 1] = fin[1];	hdf5.acc[ii * 3 + 2] = fin[2];
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
    writeZindaijiFile(Ntot, hdf5, eps, kind, bodyHead, bodyNum, bodyType, time, file, id);
  }/* if( dump_file ){ */
  else{
    fprintf(stdout, "# \"%s\" was not updated for reducing the elapsed time.\n", filename);
    fflush(stdout);
  }/* else{ */
#endif//HDF5_FOR_ZINDAIJI


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
