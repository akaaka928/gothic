/**
 * @file hsc.c
 *
 * @brief Generate surface density maps in the observation fields of Subaru/HSC
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2019/11/09 (Sat)
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


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);

  /** initialization */
  if( argc < 5 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 5);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 5 ){ */

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));

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
