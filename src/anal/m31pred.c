/**
 * @file m31pred.c
 *
 * @brief Generate surface density maps for future observations in the observed coordinate (for N-body runs in M31's disk coordinate)
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2020/11/26 (Thu)
 *
 * Copyright (C) 2020 Yohei Miki
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
#include "../anal/m31coord.h"


extern const double      length2astro;
extern const double        time2astro;
extern const double        mass2astro;
extern const double    velocity2astro;


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



void generateSurfaceDensityMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_aos *body,
  real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta,
  const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz,
  const int nvxi, const real vximin, const real dvxi, real * restrict f_xi_vxi, real * restrict f_eta_vxi, real * restrict f_D_vxi,
  const int nveta, const real vetamin, const real dveta, real * restrict f_xi_veta, real * restrict f_eta_veta, real * restrict f_D_veta);


#ifdef  USE_HDF5_FORMAT
void readM31coordinateData
(char file[], const uint id, double * restrict time, hdf5struct type,
 int * restrict kind, int * restrict bodyHead, int * restrict bodyNum,
 real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos);

void writeM31predictionData
(const double time, char file[], const uint id, hdf5struct type, const int kind,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nz, real * restrict zz,
 const int nvxi, real * restrict vxi, real * restrict f_xi_vxi, real * restrict f_eta_vxi, real * restrict f_D_vxi,
 const int nveta, real * restrict veta, real * restrict f_xi_veta, real * restrict f_eta_veta, real * restrict f_D_veta);
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
    __FPRINTF__(stderr, "          -nvxi=<int> -vximin=<real> -vximax=<real>\n");
    __FPRINTF__(stderr, "          -nveta=<int> -vetamin=<real> -vetamax=<real>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 20 ){ */

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));

  int nx;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nxi" , &nx));
  int ny;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "neta", &ny));
  int nz;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nD"  , &nz));
  real xmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,  "ximax", &xmin));/**< x = -xi */
  real ymin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "etamin", &ymin));/**< y = eta */
  real zmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,   "Dmin", &zmin));/**< z = D */
  real xmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,  "ximin", &xmax));/**< x = -xi */
  real ymax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "etamax", &ymax));/**< y = eta */
  real zmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv,   "Dmax", &zmax));/**< z = D */

  int nvxi;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nvxi", &nvxi));
  real vximin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vximin", &vximin));
  real vximax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vximax", &vximax));
  int nveta;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nveta", &nveta));
  real vetamin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vetamin", &vetamin));
  real vetamax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vetamax", &vetamax));

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
  body = (nbody_aos *)malloc(sizeof(nbody_aos) * Ntot);
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
  for(int ii = 0; ii < kind; ii++)
    checker &= (1 == fscanf(fp, "%d", &bodyNum[ii]));
  fclose(fp);
  if( !checker ){
    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
  }
  bodyHead[0] = 0;
  for(int ii = 1; ii < kind; ii++)
    bodyHead[ii] = bodyHead[ii - 1] + bodyNum[ii - 1];


  /** assume zone-centered mapping */
  /* const real deg2kpc = zm31 * CAST_D2R(tan(1.0 * M_PI / 180.0)); */
  real *xx;  xx = (real *)malloc(sizeof(real) * (nx + 1));  if( xx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xx");  }
  real *yy;  yy = (real *)malloc(sizeof(real) * (ny + 1));  if( yy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate yy");  }
  real *zz;  zz = (real *)malloc(sizeof(real) * (nz + 1));  if( zz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate zz");  }
  zmin += zm31;
  zmax += zm31;
  const real dx = (xmax - xmin) / (real)nx;  for(int ii = 0; ii < nx + 1; ii++)    xx[ii] = xmin + dx * (real)ii;
  const real dy = (ymax - ymin) / (real)ny;  for(int jj = 0; jj < ny + 1; jj++)    yy[jj] = ymin + dy * (real)jj;
  const real dz = (zmax - zmin) / (real)nz;  for(int kk = 0; kk < nz + 1; kk++)    zz[kk] = zmin + dz * (real)kk;

  real *vx;  vx = (real *)malloc(sizeof(real) * (nvxi  + 1));  if( vx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vx");  }
  real *vy;  vy = (real *)malloc(sizeof(real) * (nveta + 1));  if( vy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vy");  }
  const real dvx = ( vximax -  vximin) / (real)nvxi ;  for(int ii = 0; ii < nvxi  + 1; ii++)    vx[ii] = vximin  + dvx * (real)ii;
  const real dvy = (vetamax - vetamin) / (real)nveta;  for(int ii = 0; ii < nveta + 1; ii++)    vy[ii] = vetamin + dvy * (real)ii;

  real *f_xi_vxi ;  f_xi_vxi  = (real *)malloc(sizeof(real) * kind * nx * nvxi );  if( f_xi_vxi  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_xi_vxi");  }
  real *f_xi_veta;  f_xi_veta = (real *)malloc(sizeof(real) * kind * nx * nveta);  if( f_xi_veta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_xi_veta");  }
  real *f_eta_vxi ;  f_eta_vxi  = (real *)malloc(sizeof(real) * kind * ny * nvxi );  if( f_eta_vxi  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_eta_vxi");  }
  real *f_eta_veta;  f_eta_veta = (real *)malloc(sizeof(real) * kind * ny * nveta);  if( f_eta_veta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_eta_veta");  }
  real *f_D_vxi ;  f_D_vxi  = (real *)malloc(sizeof(real) * kind * nz * nvxi );  if( f_D_vxi  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_D_vxi");  }
  real *f_D_veta;  f_D_veta = (real *)malloc(sizeof(real) * kind * nz * nveta);  if( f_D_veta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_D_veta");  }


  for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){
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
    qsort(body, Ntot, sizeof(nbody_aos), idxAscendingOrder);

#ifdef  USE_HDF5_FORMAT
    readM31coordinateData(file, filenum, &time, hdf5type, &kind, bodyHead, bodyNum, xi, eta, dist, vxi, veta, vlos);
#endif//USE_HDF5_FORMAT

    /** obtain surface density maps for visualization */
    generateSurfaceDensityMaps(kind, bodyHead, bodyNum, body, xi, eta, dist, vxi, veta, eps, nx, xmin, dx, ny, ymin, dy, nz, zmin, dz, nvxi, vximin, dvx, f_xi_vxi, f_eta_vxi, f_D_vxi, nveta, vetamin, dvy, f_xi_veta, f_eta_veta, f_D_veta);

    /** dump analyzed results for matplotlib and/or VisIt */
#ifdef  USE_HDF5_FORMAT
    writeM31predictionData
      (time, file, filenum, hdf5type, kind,
       nx, xx, ny, yy, nz, zz,
       nvxi, vx, f_xi_vxi, f_eta_vxi, f_D_vxi,
       nveta, vy, f_xi_veta, f_eta_veta, f_D_veta);
#endif//USE_HDF5_FORMAT
  }/* for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){ */


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
  free(vx);  free(vy);

  free(f_xi_vxi );  free(f_eta_vxi );  free(f_D_vxi );
  free(f_xi_veta);  free(f_eta_veta);  free(f_D_veta);

  exitMPI();

  return (0);
}


/* 3 sigma means that neglecting component of 0.26% */
/* 5 sigma means that neglecting component of 6e-7 */
/* #define SPREAD (3) */
#define SPREAD (5)
#define SMOOTHING_FOR_VISUALIZATION TWO
/* #define SMOOTHING_FOR_VISUALIZATION THREE */
void generateSurfaceDensityMaps(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_aos *body,
  real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta,
  const real eps, const int nx, const real xmin, const real dx, const int ny, const real ymin, const real dy, const int nz, const real zmin, const real dz,
  const int nvxi, const real vximin, const real dvxi, real * restrict f_xi_vxi, real * restrict f_eta_vxi, real * restrict f_D_vxi,
  const int nveta, const real vetamin, const real dveta, real * restrict f_xi_veta, real * restrict f_eta_veta, real * restrict f_D_veta)
{
  __NOTE__("%s\n", "start");

  const real deg2kpc = zm31 * CAST_D2R(tan(1.0 * M_PI / 180.0));
  const real kpc2deg = UNITY / deg2kpc;

  const real dxinv = UNITY / dx;
  const real dyinv = UNITY / dy;
  const real dzinv = UNITY / dz;
  const real dvxiinv = UNITY / dvxi;
  const real dvetainv = UNITY / dveta;

  const real xsig = FMAX(HALF * eps * kpc2deg, SMOOTHING_FOR_VISUALIZATION * FABS(dx));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real ysig = FMAX(HALF * eps * kpc2deg, SMOOTHING_FOR_VISUALIZATION * FABS(dy));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real zsig = FMAX(HALF * eps          , SMOOTHING_FOR_VISUALIZATION * FABS(dz));/**< sigma in Gaussian is roughly corresponding to eps of Plummer softening */
  const real vxisig = SMOOTHING_FOR_VISUALIZATION * FABS(dvxi);
  const real vetasig = SMOOTHING_FOR_VISUALIZATION * FABS(dveta);

  const real invxsig = UNITY / xsig;
  const real invysig = UNITY / ysig;
  const real invzsig = UNITY / zsig;
  const real invvxisig = UNITY / vxisig;
  const real invvetasig = UNITY / vetasig;

  const int nx_smooth = (int)CEIL(SPREAD * FABS(dxinv) * xsig);
  const int ny_smooth = (int)CEIL(SPREAD * FABS(dyinv) * ysig);
  const int nz_smooth = (int)CEIL(SPREAD * FABS(dzinv) * zsig);
  const int nvxi_smooth = (int)CEIL(SPREAD * FABS(dvxiinv) * vxisig);
  const int nveta_smooth = (int)CEIL(SPREAD * FABS(dvetainv) * vetasig);

  real *erfx;  erfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 2));  if( erfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfx");  }
  real *erfy;  erfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 2));  if( erfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfy");  }
  real *erfz;  erfz = (real *)malloc(sizeof(real) * (2 * nz_smooth + 2));  if( erfz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfz");  }
  real *erfvxi;  erfvxi = (real *)malloc(sizeof(real) * (2 * nvxi_smooth + 2));  if( erfvxi == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfvxi");  }
  real *erfveta;  erfveta = (real *)malloc(sizeof(real) * (2 * nveta_smooth + 2));  if( erfveta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfveta");  }
  for(int ii = 0; ii < 2 * nx_smooth + 2; ii++)    erfx[ii] = ERF((dx * ((real)(ii - nx_smooth) - HALF)) * invxsig);
  for(int ii = 0; ii < 2 * ny_smooth + 2; ii++)    erfy[ii] = ERF((dy * ((real)(ii - ny_smooth) - HALF)) * invysig);
  for(int ii = 0; ii < 2 * nz_smooth + 2; ii++)    erfz[ii] = ERF((dz * ((real)(ii - nz_smooth) - HALF)) * invzsig);
  for(int ii = 0; ii < 2 * nvxi_smooth + 2; ii++)    erfvxi[ii] = ERF((dvxi * ((real)(ii - nvxi_smooth) - HALF)) * invvxisig);
  for(int ii = 0; ii < 2 * nveta_smooth + 2; ii++)    erfveta[ii] = ERF((dveta * ((real)(ii - nveta_smooth) - HALF)) * invvetasig);

  real *psfx;  psfx = (real *)malloc(sizeof(real) * (2 * nx_smooth + 1));  if( psfx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfx");  }
  real *psfy;  psfy = (real *)malloc(sizeof(real) * (2 * ny_smooth + 1));  if( psfy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfy");  }
  real *psfz;  psfz = (real *)malloc(sizeof(real) * (2 * nz_smooth + 1));  if( psfz == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfz");  }
  real *psfvxi;  psfvxi = (real *)malloc(sizeof(real) * (2 * nvxi_smooth + 1));  if( psfvxi == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfvxi");  }
  real *psfveta;  psfveta = (real *)malloc(sizeof(real) * (2 * nveta_smooth + 1));  if( psfveta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfveta");  }
  for(int ii = 0; ii < 2 * nx_smooth + 1; ii++)    psfx[ii] = HALF * (erfx[ii + 1] - erfx[ii]);
  for(int ii = 0; ii < 2 * ny_smooth + 1; ii++)    psfy[ii] = HALF * (erfy[ii + 1] - erfy[ii]);
  for(int ii = 0; ii < 2 * nz_smooth + 1; ii++)    psfz[ii] = HALF * (erfz[ii + 1] - erfz[ii]);
  for(int ii = 0; ii < 2 * nvxi_smooth + 1; ii++)    psfvxi[ii] = HALF * (erfvxi[ii + 1] - erfvxi[ii]);
  for(int ii = 0; ii < 2 * nveta_smooth + 1; ii++)    psfveta[ii] = HALF * (erfveta[ii + 1] - erfveta[ii]);



#ifdef  USE_DEGREE_FOR_SURFACE_DENSITY_MAP
  /** show Msun deg^-2 */
  const real dfinv_xi_vxi   = dxinv * dvxiinv;
  const real dfinv_xi_veta  = dxinv * dvetainv;
  const real dfinv_eta_vxi  = dyinv * dvxiinv;
  const real dfinv_eta_veta = dyinv * dvetainv;
  const real dfinv_D_vxi    = dzinv * dvxiinv;
  const real dfinv_D_veta   = dzinv * dvetainv;
#else///USE_DEGREE_FOR_SURFACE_DENSITY_MAP
  /** show Msun kpc^-2 */
  const real inv_deg2kpc = kpc2deg;
  const real dfinv_xi_vxi   = dxinv * inv_deg2kpc * dvxiinv;
  const real dfinv_xi_veta  = dxinv * inv_deg2kpc * dvetainv;
  const real dfinv_eta_vxi  = dyinv * inv_deg2kpc * dvxiinv;
  const real dfinv_eta_veta = dyinv * inv_deg2kpc * dvetainv;
  const real dfinv_D_vxi    = dzinv * dvxiinv;
  const real dfinv_D_veta   = dzinv * dvetainv;
#endif//USE_DEGREE_FOR_SURFACE_DENSITY_MAP

  for(int kk = 0; kk < kind; kk++){
    for(int ii = 0; ii < nx * nvxi ; ii++)      f_xi_vxi [INDEX2D(kind, nx * nvxi , kk, ii)] = ZERO;
    for(int ii = 0; ii < nx * nveta; ii++)      f_xi_veta[INDEX2D(kind, nx * nveta, kk, ii)] = ZERO;
    for(int ii = 0; ii < ny * nvxi ; ii++)      f_eta_vxi [INDEX2D(kind, ny * nvxi , kk, ii)] = ZERO;
    for(int ii = 0; ii < ny * nveta; ii++)      f_eta_veta[INDEX2D(kind, ny * nveta, kk, ii)] = ZERO;
    for(int ii = 0; ii < nz * nvxi ; ii++)      f_D_vxi [INDEX2D(kind, nz * nvxi , kk, ii)] = ZERO;
    for(int ii = 0; ii < nz * nveta; ii++)      f_D_veta[INDEX2D(kind, nz * nveta, kk, ii)] = ZERO;

    for(int ii = bodyHead[kk]; ii < bodyHead[kk] + bodyNum[kk]; ii++){
      const real mi = CAST_D2R(CAST_R2D(body[ii].m) * mass2astro);
      const real xx =   xi[ii];
      const real yy =  eta[ii];
      const real zz = dist[ii];
      const real vx = vxi[ii];
      const real vy = veta[ii];

      const int l0 = (int)FLOOR((xx - xmin) * dxinv);
      const int m0 = (int)FLOOR((yy - ymin) * dyinv);
      const int n0 = (int)FLOOR((zz - zmin) * dzinv);
      const int o0 = (int)FLOOR((vx - vximin) * dvxiinv);
      const int p0 = (int)FLOOR((vy - vetamin) * dvetainv);



      for(int svx = 0; svx < 2 * nvxi_smooth + 1; svx++){
        const int oo = o0 + svx - nvxi_smooth;
        if((oo >= 0) && (oo < nvxi)){
          const real mvxi = mi * psfvxi[svx];

          for(int sx = 0; sx < 2 * nx_smooth + 1; sx++){
      	    const int ll = l0 + sx - nx_smooth;
	          if((ll >= 0) && (ll < nx))
	            f_xi_vxi[INDEX3D(kind, nx, nvxi, kk, ll, oo)] += mvxi * psfx[sx];
          }

          for(int sy = 0; sy < 2 * ny_smooth + 1; sy++){
      	    const int mm = m0 + sy - ny_smooth;
	          if((mm >= 0) && (mm < ny))
	            f_eta_vxi[INDEX3D(kind, ny, nvxi, kk, mm, oo)] += mvxi * psfy[sy];
          }

          for(int sz = 0; sz < 2 * nz_smooth + 1; sz++){
      	    const int nn = n0 + sz - nz_smooth;
	          if((nn >= 0) && (nn < nz))
	            f_D_vxi[INDEX3D(kind, nz, nvxi, kk, nn, oo)] += mvxi * psfz[sz];
          }

        }
      }


      for(int svy = 0; svy < 2 * nveta_smooth + 1; svy++){
        const int pp = p0 + svy - nveta_smooth;
        if((pp >= 0) && (pp < nveta)){
          const real mveta = mi * psfveta[svy];

          for(int sx = 0; sx < 2 * nx_smooth + 1; sx++){
      	    const int ll = l0 + sx - nx_smooth;
	          if((ll >= 0) && (ll < nx))
	            f_xi_veta[INDEX3D(kind, nx, nveta, kk, ll, pp)] += mveta * psfx[sx];
          }

          for(int sy = 0; sy < 2 * ny_smooth + 1; sy++){
      	    const int mm = m0 + sy - ny_smooth;
	          if((mm >= 0) && (mm < ny))
	            f_eta_veta[INDEX3D(kind, ny, nveta, kk, mm, pp)] += mveta * psfy[sy];
          }

          for(int sz = 0; sz < 2 * nz_smooth + 1; sz++){
      	    const int nn = n0 + sz - nz_smooth;
	          if((nn >= 0) && (nn < nz))
	            f_D_veta[INDEX3D(kind, nz, nveta, kk, nn, pp)] += mveta * psfz[sz];
          }

        }
      }

    }/* for(int ii = 0; ii < bodyNum[kk]; ii++){ */

    for(int ii = 0; ii < nx * nvxi ; ii++)      f_xi_vxi [INDEX2D(kind, nx * nvxi , kk, ii)] *= dfinv_xi_vxi;
    for(int ii = 0; ii < nx * nveta; ii++)      f_xi_veta[INDEX2D(kind, nx * nveta, kk, ii)] *= dfinv_xi_veta;
    for(int ii = 0; ii < ny * nvxi ; ii++)      f_eta_vxi [INDEX2D(kind, ny * nvxi , kk, ii)] *= dfinv_eta_vxi;
    for(int ii = 0; ii < ny * nveta; ii++)      f_eta_veta[INDEX2D(kind, ny * nveta, kk, ii)] *= dfinv_eta_veta;
    for(int ii = 0; ii < nz * nvxi ; ii++)      f_D_vxi [INDEX2D(kind, nz * nvxi , kk, ii)] *= dfinv_D_vxi;
    for(int ii = 0; ii < nz * nveta; ii++)      f_D_veta[INDEX2D(kind, nz * nveta, kk, ii)] *= dfinv_D_veta;
  }/* for(int kk = 0; kk < kind; kk++){ */

  free(erfx);  free(erfy);  free(erfz);  free(erfvxi);  free(erfveta);
  free(psfx);  free(psfy);  free(psfz);  free(psfvxi);  free(psfveta);


  __NOTE__("%s\n", "end");
}


#ifdef  USE_HDF5_FORMAT
/**
 * @fn readM31coordinateData
 *
 * @brief Read analyzed profiles of the N-body simulation.
 */
void readM31coordinateData
(char file[], const uint id, double * restrict time, hdf5struct type,
 int * restrict kind, int * restrict bodyHead, int * restrict bodyNum,
 real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos)
{
  __NOTE__("%s\n", "start");

  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "m31obs", id);
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);


  /** read attribute data */
  hid_t attribute;
  /* read current time in units of astrophysical unit */
  attribute = H5Aopen(target, "time", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, time));
  chkHDF5err(H5Aclose(attribute));
  /* /\* read # of steps *\/ */
  /* attribute = H5Aopen(target, "steps", H5P_DEFAULT); */
  /* chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, steps)); */
  /* chkHDF5err(H5Aclose(attribute)); */
  /* /\* read # of N-body particles *\/ */
  /* ulong num_ulong = 0; */
  /* attribute = H5Aopen(target, "number", H5P_DEFAULT); */
  /* chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, &num_ulong)); */
  /* chkHDF5err(H5Aclose(attribute)); */
  /* read flag about USE_DOUBLE_PRECISION */
  attribute = H5Aopen(target, "kinds", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, kind));
  chkHDF5err(H5Aclose(attribute));
  int useDP;
  attribute = H5Aopen(target, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
#ifdef  USE_DOUBLE_PRECISION
  if( useDP != 1 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);
  }/* if( useDP != 1 ){ */
#else///USE_DOUBLE_PRECISION
  if( useDP != 0 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);
  }/* if( useDP != 0 ){ */
#endif//USE_DOUBLE_PRECISION


  bodyHead[0] = 0;
  for(int kk = 0; kk < *kind; kk++){
    /* open an existing group */
    char groupname[16];
    sprintf(groupname, "obs%d", kk);

    bodyNum[kk] = 0;
    if(H5Lexists(target, groupname, H5P_DEFAULT)){
      hid_t group = H5Gopen(target, groupname, H5P_DEFAULT);

      attribute = H5Aopen(group, "num", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &bodyNum[kk]));
      chkHDF5err(H5Aclose(attribute));

      hid_t dataset;
      dataset = H5Dopen(group, "xi", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &xi[bodyHead[kk]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dopen(group, "eta", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &eta[bodyHead[kk]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dopen(group, "dist", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dist[bodyHead[kk]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dopen(group, "vxi", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vxi[bodyHead[kk]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dopen(group, "veta", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &veta[bodyHead[kk]]));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dopen(group, "vlos", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vlos[bodyHead[kk]]));
      chkHDF5err(H5Dclose(dataset));

      chkHDF5err(H5Gclose(group));
    }
    if((kk + 1) < (*kind))
      bodyHead[kk + 1] = bodyHead[kk] + bodyNum[kk];
  }

  /** close the file */
  chkHDF5err(H5Fclose(target));


  __NOTE__("%s\n", "end");
}


/**
 * @fn writeM31predictionData
 *
 * @brief Write analyzed profiles of the N-body simulation.
 */
void writeM31predictionData
(const double time, char file[], const uint id, hdf5struct type, const int kind,
 const int nx, real * restrict xx, const int ny, real * restrict yy, const int nz, real * restrict zz,
 const int nvxi, real * restrict vxi, real * restrict f_xi_vxi, real * restrict f_eta_vxi, real * restrict f_D_vxi,
 const int nveta, real * restrict veta, real * restrict f_xi_veta, real * restrict f_eta_veta, real * restrict f_D_veta)
{
  __NOTE__("%s\n", "start");

  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "m31prd", id);
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
  attribute = H5Acreate(target, "nvxi", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nvxi));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nveta", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nveta));
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

    /** 2D (nx * nvxi) array */
    hsize_t dims[2] = {nx, nvxi};
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
    dataset = H5Dcreate(group, "xi-vxi", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &f_xi_vxi[INDEX2D(kind, nx * nvxi, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (nx * nveta) array */
    dims[0] = nx;
    dims[1] = nveta;
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
    dataset = H5Dcreate(group, "xi-veta", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &f_xi_veta[INDEX2D(kind, nx * nveta, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (ny * nvxi) array */
    dims[0] = ny;
    dims[1] = nvxi;
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
    dataset = H5Dcreate(group, "eta-vxi", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &f_eta_vxi[INDEX2D(kind, ny * nvxi, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (ny * nveta) array */
    dims[0] = ny;
    dims[1] = nveta;
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
    dataset = H5Dcreate(group, "eta-veta", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &f_eta_veta[INDEX2D(kind, ny * nveta, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (nz * nvxi) array */
    dims[0] = nz;
    dims[1] = nvxi;
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
    dataset = H5Dcreate(group, "D-vxi", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &f_D_vxi[INDEX2D(kind, nz * nvxi, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));

    /** 2D (nz * nveta) array */
    dims[0] = nz;
    dims[1] = nveta;
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
    dataset = H5Dcreate(group, "D-veta", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &f_D_veta[INDEX2D(kind, nz * nveta, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
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
    chkHDF5err(H5Sclose(dataspace));


    /** 1D (nvxi) array */
    dims[0] = nvxi + 1;
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
    dataset = H5Dcreate(group, "vxi", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, vxi));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nveta) array */
    dims[0] = nveta + 1;
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
    dataset = H5Dcreate(group, "veta", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, veta));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    chkHDF5err(H5Sclose(dataspace));

    chkHDF5err(H5Gclose(group));
  }/* for(int ii = 0; ii < kind; ii++){ */

  /* close the file */
  chkHDF5err(H5Fclose(target));

  __NOTE__("%s\n", "end");
}


#endif//USE_HDF5_FORMAT
