/**
 * @file io.c
 *
 * @brief Source code for Input/Output functions in GOTHIC and MAGI
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/12/17 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

/**
 * @def USE_SZIP_COMPRESSION
 *
 * @brief On to enable Szip compression for HDF5 files (default is OFF).
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
#include <math.h>

#ifndef DISABLE_MPI
#include <mpi.h>
#include "mpilib.h"
#endif//DISABLE_MPI

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
#include "name.h"
#include "constants.h"

#include "../misc/structure.h"
#ifndef RUN_WITHOUT_GOTHIC
#include "../misc/benchmark.h"
#include "../misc/tune.h"
#include "../misc/brent.h"
#endif//RUN_WITHOUT_GOTHIC

#include "io.h"

#ifndef RUN_WITHOUT_GOTHIC
#ifdef  EXEC_BENCHMARK
#       include <unistd.h>
#endif//EXEC_BENCHMARK
#ifdef  HUNT_WALK_PARAMETER
#       include "../tree/walk_dev.h"
#       include "../tree/neighbor_dev.h"
#       include "../tree/shrink_dev.h"
#endif//HUNT_WALK_PARAMETER
#ifdef  HUNT_MAKE_PARAMETER
#       include "../sort/peano_dev.h"
#       include "../tree/make_dev.h"
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_NODE_PARAMETER
#       include "../tree/make_dev.h"
#endif//HUNT_NODE_PARAMETER
#ifdef  HUNT_TIME_PARAMETER
#       include "../time/adv_dev.h"
#endif//HUNT_TIME_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
#       include "../tree/neighbor_dev.h"
#endif//HUNT_FIND_PARAMETER
#endif//RUN_WITHOUT_GOTHIC


/* global constants to set unit system, defined in constants.c */
extern const double        time2astro;extern const char        time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
#ifdef  USE_HDF5_FORMAT
extern const real newton;
extern const double      length2astro;extern const char      length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double        mass2astro;extern const char        mass_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double     density2astro;extern const char     density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double col_density2astro;extern const char col_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double      energy2astro;extern const char      energy_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double     senergy2astro;extern const char     senergy_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double    velocity2astro;extern const char    velocity_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double       accel2astro;extern const char       accel_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
#endif//USE_HDF5_FORMAT


#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
/**
 * @fn createMPIcfg_dataio
 *
 * @brief Configure MPI communicator for data I/O.
 *
 * @return (cfg) MPI communicator for data I/O
 * @param (mpi) current MPI communicator
 */
void createMPIcfg_dataio(MPIcfg_dataio *cfg, const MPIinfo mpi)
{
  __NOTE__("%s\n", "start");

  cfg->comm = mpi.comm;
  cfg->info = mpi.info;
  cfg->size = mpi.size;
  cfg->rank = mpi.rank;

  __NOTE__("%s\n", "end");
}
/**
 * @fn updateMPIcfg_dataio
 *
 * @brief Update MPI communicator for data I/O.
 *
 * @return (cfg) MPI communicator for data I/O
 * @param (num) number of N-body particles contained in this MPI process
 *
 * @detail The partition of particles' array is updated using MPI_Scan.
 */
static inline void updateMPIcfg_dataio(MPIcfg_dataio *cfg, const int num)
{
  __NOTE__("%s\n", "start");

  static ulong unit, psum;
  unit = (ulong)num;
  chkMPIerr(MPI_Scan(&unit, &psum, 1, MPI_UNSIGNED_LONG, MPI_SUM, cfg->comm));
  cfg->head = (MPI_Offset)(psum - unit);

  __NOTE__("%s\n", "end");
}
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)


static inline void writeConfigFileFormat(char file[])
{
  __NOTE__("%s\n", "start");
  fprintf(stderr, "Name of the configuration file should be \"%s\" (argv[1])\n", file);
  fprintf(stderr, "\tline 0: index of last updated tentative file series\n");
  __KILL__(stderr, "%s\n", "ERROR: format of the configuration file wrong.");
  __NOTE__("%s\n", "end");
}

/**
 * @fn updateConfigFile
 *
 * @brief Update the file index of the restarter file.
 *
 * @param (last) index of the latest restarter file
 * @param (file) name of the restarter file
 */
void updateConfigFile(int last, char file[])
{
  __NOTE__("%s\n", "start");

  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s", DATAFOLDER, file, CONFIG);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fprintf(fp, "%d\n", last);
  fclose(fp);

  __NOTE__("%s\n", "end");
}
/**
 * @fn readConfigFile
 *
 * @brief Read the file index of the latest restarter file.
 *
 * @return (last) index of the latest restarter file
 * @param (file) name of the restarter file
 */
void readConfigFile(int *last, char file[])
{
  __NOTE__("%s\n", "start");

  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s", DATAFOLDER, file, CONFIG);
  fp = fopen(filename, "r");
  if( fp == NULL ){
    writeConfigFileFormat(filename);
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }/* if( fp == NULL ){ */
  int checker = 1;
  checker &= (1 == fscanf(fp, "%d", last));
  fclose(fp);
  if( !checker )    writeConfigFileFormat(filename);

  __NOTE__("%s\n", "end");
}


#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
/**
 * @fn readConfigFileParallel
 *
 * @brief Read the file index of the latest restarter file with MPI-IO.
 *
 * @return (last) index of the latest restarter file
 * @param (file) name of the restarter file
 * @param (mpi) information on MPI process
 *
 * @sa readConfigFile
 */
void readConfigFileParallel(int *last, char file[], MPIinfo mpi)
{
  __NOTE__("%s\n", "start");

  /* only the root process read the specified file */
  if( mpi.rank == 0 )
    readConfigFile(last, file);

  /* broadcast the read parameters from the root process */
  chkMPIerr(MPI_Bcast(last, 1, MPI_INT, 0, mpi.comm));

  __NOTE__("%s\n", "end");
}
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)


/**
 * @fn readSettings
 *
 * @brief Read the settings about the N-body simulation.
 *
 * @return (unit) unit system of the simulation
 * @return (Ntot) total number of N-body particles
 * @return (eps) value of Plummer softening length
 * @return (eta) value of safety parameter to determine the time step
 * @return (ft) finish time of the simulation
 * @return (SnapshotInterval) time interval to write snapshot files
 * @return (SaveInterval) time interval to dump the tentative results of the simulation
 * @param (file) name of the simulation
 */
void readSettings(int *unit, ulong *Ntot, real *eps, real *eta, double *ft, double *SnapshotInterval, double *SaveInterval, char file[])
{
  __NOTE__("%s\n", "start");

  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s", DATAFOLDER, file, SETTINGS);
  fp = fopen(filename, "rb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  /* read values from the file */
  bool success = true;
  size_t tmp;
  tmp = 1;  if( tmp != fread(            Ntot, sizeof( ulong), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(             eps, sizeof(  real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(             eta, sizeof(  real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(              ft, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(SnapshotInterval, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(    SaveInterval, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(            unit, sizeof(   int), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }
  fclose(fp);

  /* set the unit system of the simulation */
  setPhysicalConstantsAndUnitSystem(*unit, 0);

  __NOTE__("%s\n", "end");
}


#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
/**
 * @fn readSettingsParallel
 *
 * @brief Read the settings about the N-body simulation with MPI-IO.
 *
 * @return (unit) unit system of the simulation
 * @return (Ntot) total number of N-body particles
 * @return (eps) value of Plummer softening length
 * @return (eta) value of safety parameter to determine the time step
 * @return (ft) finish time of the simulation
 * @return (SnapshotInterval) time interval to write snapshot files
 * @return (SaveInterval) time interval to dump the tentative results of the simulation
 * @param (file) name of the simulation
 * @param (mpi) information on MPI process
 *
 * @sa readSettings
 */
void readSettingsParallel(int *unit, ulong *Ntot, real *eps, real *eta, double *ft, double *SnapshotInterval, double *SaveInterval, char file[], MPIinfo mpi)
{
  __NOTE__("%s\n", "start");

  /* only the root process read the specified file */
  if( mpi.rank == 0 )
    readSettings(unit, Ntot, eps, eta, ft, SnapshotInterval, SaveInterval, file);

  /* broadcast the read parameters from the root process */
  chkMPIerr(MPI_Bcast(            Ntot, 1, MPI_UNSIGNED_LONG, 0, mpi.comm));
  chkMPIerr(MPI_Bcast(             eps, 1, MPI_REALDAT,       0, mpi.comm));
  chkMPIerr(MPI_Bcast(             eta, 1, MPI_REALDAT,       0, mpi.comm));
  chkMPIerr(MPI_Bcast(              ft, 1, MPI_DOUBLE,        0, mpi.comm));
  chkMPIerr(MPI_Bcast(SnapshotInterval, 1, MPI_DOUBLE,        0, mpi.comm));
  chkMPIerr(MPI_Bcast(    SaveInterval, 1, MPI_DOUBLE,        0, mpi.comm));
  chkMPIerr(MPI_Bcast(            unit, 1, MPI_INT,           0, mpi.comm));

  /* set the unit system of the simulation */
  setPhysicalConstantsAndUnitSystem(*unit, 0);

  __NOTE__("%s\n", "end");
}
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)


/**
 * @fn writeSettings
 *
 * @brief Write the settings about the N-body simulation.
 *
 * @param (unit) unit system of the simulation
 * @param (Ntot) total number of N-body particles
 * @param (eps) value of Plummer softening length
 * @param (eta) value of safety parameter to determine the time step
 * @param (ft) finish time of the simulation
 * @param (SnapshotInterval) time interval to write snapshot files
 * @param (SaveInterval) time interval to dump the tentative results of the simulation
 * @param (file) name of the simulation
 */
void writeSettings(int unit, ulong Ntot, real eps, real eta, double ft, double SnapshotInterval, double SaveInterval, char file[])
{
  __NOTE__("%s\n", "start");

  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s", DATAFOLDER, file, SETTINGS);
  fp = fopen(filename, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  bool success = true;
  size_t tmp;
  tmp = 1;  if( tmp != fwrite(&            Ntot, sizeof( ulong), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&             eps, sizeof(  real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&             eta, sizeof(  real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&              ft, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&SnapshotInterval, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&    SaveInterval, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&            unit, sizeof(   int), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }
  fclose(fp);

  __NOTE__("%s\n", "end");
}


#ifdef  USE_HDF5_FORMAT
/**
 * @fn createHDF5DataType
 *
 * @brief Create data type for HDF5.
 *
 * @return (type) data type for HDF5
 */
void createHDF5DataType(hdf5struct *type)
{
  __NOTE__("%s\n", "start");

  /* commit a data type to write unit name as strings */
  type->str4unit = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(type->str4unit, CONSTANTS_H_CHAR_WORDS));/* memo: sizeof(char) is unity */

  /* commit a data type to switch float and double */
#ifdef  USE_DOUBLE_PRECISION
  type->real = H5T_NATIVE_DOUBLE;
#else///USE_DOUBLE_PRECISION
  type->real = H5T_NATIVE_FLOAT;
#endif//USE_DOUBLE_PRECISION

  /* commit a data type of rebuildTree for GOTHIC */
#ifndef RUN_WITHOUT_GOTHIC
  type->rebuildTree = H5Tcreate(H5T_COMPOUND, sizeof(rebuildTree));
  chkHDF5err(H5Tinsert(type->rebuildTree, "interval", HOFFSET(rebuildTree, interval), H5T_NATIVE_DOUBLE));
#   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
  chkHDF5err(H5Tinsert(type->rebuildTree, "avg", HOFFSET(rebuildTree, avg), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->rebuildTree, "var", HOFFSET(rebuildTree, var), H5T_NATIVE_DOUBLE));
#endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
  chkHDF5err(H5Tinsert(type->rebuildTree, "reuse", HOFFSET(rebuildTree, reuse), H5T_NATIVE_INT));
#ifdef  BLOCK_TIME_STEP
  chkHDF5err(H5Tinsert(type->rebuildTree, "adjust", HOFFSET(rebuildTree, adjust), H5T_NATIVE_BOOL));
#endif//BLOCK_TIME_STEP
#endif//RUN_WITHOUT_GOTHIC

  /* commit a data type of measuredTime for GOTHIC */
#ifndef RUN_WITHOUT_GOTHIC
  type->measuredTime = H5Tcreate(H5T_COMPOUND, sizeof(measuredTime));
  chkHDF5err(H5Tinsert(type->measuredTime, "walkTree[0]", HOFFSET(measuredTime, walkTree[0]), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->measuredTime, "walkTree[1]", HOFFSET(measuredTime, walkTree[1]), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->measuredTime, "makeTree"   , HOFFSET(measuredTime, makeTree   ), H5T_NATIVE_DOUBLE));
#ifndef SERIALIZED_EXECUTION
  chkHDF5err(H5Tinsert(type->measuredTime,  "sum_excg"   , HOFFSET(measuredTime, sum_excg   ), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->measuredTime,  "sum_rebuild", HOFFSET(measuredTime, sum_rebuild), H5T_NATIVE_DOUBLE));
#ifdef  MONITOR_LETGEN_TIME
  chkHDF5err(H5Tinsert(type->measuredTime,  "excg"       , HOFFSET(measuredTime, excg       ), H5T_NATIVE_DOUBLE));
#endif//MONITOR_LETGEN_TIME
#endif//SERIALIZED_EXECUTION
#endif//RUN_WITHOUT_GOTHIC


  /* commit data types for auto-tuning in GOTHIC */
#ifndef RUN_WITHOUT_GOTHIC
  /* commit a data type of statVal */
  type->statVal = H5Tcreate(H5T_COMPOUND, sizeof(statVal));
  chkHDF5err(H5Tinsert(type->statVal, "S"  , HOFFSET(statVal, S  ), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->statVal, "Sx" , HOFFSET(statVal, Sx ), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->statVal, "Sy" , HOFFSET(statVal, Sy ), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->statVal, "Sxx", HOFFSET(statVal, Sxx), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->statVal, "Sxy", HOFFSET(statVal, Sxy), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->statVal, "Syy", HOFFSET(statVal, Syy), H5T_NATIVE_DOUBLE));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  chkHDF5err(H5Tinsert(type->statVal, "Sxxxx", HOFFSET(statVal, Sxxxx), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->statVal, "Sxxx" , HOFFSET(statVal, Sxxx ), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->statVal, "Sxxy" , HOFFSET(statVal, Sxxy ), H5T_NATIVE_DOUBLE));
#endif//USE_PARABOLIC_GROWTH_MODEL

  /* commit a data type of guessTime */
  type->guessTime = H5Tcreate(H5T_COMPOUND, sizeof(guessTime));
  chkHDF5err(H5Tinsert(type->guessTime,  "slope", HOFFSET(guessTime,  slope), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->guessTime,  "icept", HOFFSET(guessTime,  icept), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->guessTime, "rchisq", HOFFSET(guessTime, rchisq), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->guessTime,   "time", HOFFSET(guessTime,   time), H5T_NATIVE_DOUBLE));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  chkHDF5err(H5Tinsert(type->guessTime, "second", HOFFSET(guessTime, second), H5T_NATIVE_DOUBLE));
#endif//USE_PARABOLIC_GROWTH_MODEL

  /* commit a data type of brentFunc */
  type->brentFunc = H5Tcreate(H5T_COMPOUND, sizeof(brentFunc));
  chkHDF5err(H5Tinsert(type->brentFunc, "pos", HOFFSET(brentFunc, pos), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->brentFunc, "val", HOFFSET(brentFunc, val), H5T_NATIVE_DOUBLE));

  /* commit a data type of brentStatus */
  type->brentStatus = H5Tcreate(H5T_COMPOUND, sizeof(brentStatus));
  chkHDF5err(H5Tinsert(type->brentStatus, "x", HOFFSET(brentStatus, x), type->brentFunc));
  chkHDF5err(H5Tinsert(type->brentStatus, "w", HOFFSET(brentStatus, w), type->brentFunc));
  chkHDF5err(H5Tinsert(type->brentStatus, "v", HOFFSET(brentStatus, v), type->brentFunc));
  chkHDF5err(H5Tinsert(type->brentStatus, "u", HOFFSET(brentStatus, u), type->brentFunc));
  chkHDF5err(H5Tinsert(type->brentStatus, "a", HOFFSET(brentStatus, a), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->brentStatus, "b", HOFFSET(brentStatus, b), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->brentStatus, "d", HOFFSET(brentStatus, d), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->brentStatus, "e", HOFFSET(brentStatus, e), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->brentStatus, "gold", HOFFSET(brentStatus, gold), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->brentStatus, "initialized", HOFFSET(brentStatus, initialized), H5T_NATIVE_BOOL));

  /* commit a data type of brentMemory */
  type->brentMemory = H5Tcreate(H5T_COMPOUND, sizeof(brentMemory));
  chkHDF5err(H5Tinsert(type->brentMemory, "previous", HOFFSET(brentMemory, previous), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->brentMemory,   "totNum", HOFFSET(brentMemory,   totNum), H5T_NATIVE_INT));
  chkHDF5err(H5Tinsert(type->brentMemory, "degraded", HOFFSET(brentMemory, degraded), H5T_NATIVE_INT));
  chkHDF5err(H5Tinsert(type->brentMemory, "interval", HOFFSET(brentMemory, interval), H5T_NATIVE_INT));
#endif//RUN_WITHOUT_GOTHIC

  /* commit a data type of pot2 */
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  type->pot2 = H5Tcreate(H5T_COMPOUND, sizeof(pot2));
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  chkHDF5err(H5Tinsert(type->pot2, "val", HOFFSET(pot2, val), type->real));
  chkHDF5err(H5Tinsert(type->pot2, "dr2", HOFFSET(pot2, dr2), type->real));
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  chkHDF5err(H5Tinsert(type->pot2, "Phi", HOFFSET(pot2, Phi), type->real));
  chkHDF5err(H5Tinsert(type->pot2,  "Fr", HOFFSET(pot2,  Fr), type->real));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
  type->disk_grav = H5Tcreate(H5T_COMPOUND, sizeof(disk_grav));
  chkHDF5err(H5Tinsert(type->disk_grav, "R", HOFFSET(disk_grav, R), type->real));
  chkHDF5err(H5Tinsert(type->disk_grav, "z", HOFFSET(disk_grav, z), type->real));
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  /* commit a data type of deviceMonitors */
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
  type->gpu_clock = H5Tcreate(H5T_COMPOUND, sizeof(gpu_clock));
  chkHDF5err(H5Tinsert(type->gpu_clock, "elapsed", HOFFSET(gpu_clock, elapsed), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->gpu_clock, "devClock", HOFFSET(gpu_clock, devClock), H5T_NATIVE_UINT));
  chkHDF5err(H5Tinsert(type->gpu_clock, "temperature", HOFFSET(gpu_clock, temperature), H5T_NATIVE_UINT));
  chkHDF5err(H5Tinsert(type->gpu_clock, "power", HOFFSET(gpu_clock, power), H5T_NATIVE_UINT));
  chkHDF5err(H5Tinsert(type->gpu_clock, "grpNum", HOFFSET(gpu_clock, grpNum), H5T_NATIVE_INT));
#endif//REPORT_GPU_CLOCK_FREQUENCY

  __NOTE__("%s\n", "end");
}
/**
 * @fn removeHDF5DataType
 *
 * @brief Remove data type for HDF5.
 *
 * @param (type) data type for HDF5
 */
void removeHDF5DataType(hdf5struct  type)
{
  __NOTE__("%s\n", "start");

#ifndef RUN_WITHOUT_GOTHIC
  chkHDF5err(H5Tclose(type.brentMemory));
  chkHDF5err(H5Tclose(type.brentStatus));
  chkHDF5err(H5Tclose(type.brentFunc));
  chkHDF5err(H5Tclose(type.guessTime));
  chkHDF5err(H5Tclose(type.statVal));
  chkHDF5err(H5Tclose(type.measuredTime));
  chkHDF5err(H5Tclose(type.rebuildTree));
#endif//RUN_WITHOUT_GOTHIC

  chkHDF5err(H5Tclose(type.str4unit));

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  chkHDF5err(H5Tclose(type.pot2));
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  chkHDF5err(H5Tclose(type.disk_grav));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#ifdef  REPORT_GPU_CLOCK_FREQUENCY
  chkHDF5err(H5Tclose(type.gpu_clock));
#endif//REPORT_GPU_CLOCK_FREQUENCY

  __NOTE__("%s\n", "end");
}
#endif//USE_HDF5_FORMAT


#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
static inline void appendGPUclockInfo
(const int monitor_step, gpu_clock *deviceMonitors,
#ifdef  USE_HDF5_FORMAT
 const hid_t target, const hdf5struct type
#else///USE_HDF5_FORMAT
 char *filename, FILE *fp
#endif//USE_HDF5_FORMAT
 )
{
  __NOTE__("%s\n", "start");


  const int monitor = (CLOCK_RECORD_STEPS < monitor_step) ? CLOCK_RECORD_STEPS : monitor_step;

  /* data reordering */
  static gpu_clock record[CLOCK_RECORD_STEPS];
  const int origin = monitor_step & (CLOCK_RECORD_STEPS - 1);
  const int turnaround = monitor - origin;
  for(int ii = 0; ii < turnaround; ii++)
    record[ii] = deviceMonitors[origin + ii];
  for(int ii = turnaround; ii < monitor; ii++)
    record[ii] = deviceMonitors[ii - turnaround];

#ifdef  USE_HDF5_FORMAT
  hid_t group = H5Gcreate(target, "GPUinfo", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* write attribute */
  hsize_t attr_dims = 1;
  hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute = H5Acreate(group, "steps", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &monitor));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));

  /* write measured data */
  hsize_t clock_dims = (hsize_t)monitor;
  dataspace = H5Screate_simple(1, &clock_dims, NULL);
  hid_t dataset = H5Dcreate(group, "record", type.gpu_clock, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.gpu_clock, H5S_ALL, H5S_ALL, H5P_DEFAULT, record));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(dataspace));

  chkHDF5err(H5Gclose(group));
#else///USE_HDF5_FORMAT
  size_t tmp;
  bool success = true;
  tmp = 1;  if( tmp != fwrite(&monitor, sizeof(int), tmp, fp) )    success = false;
  tmp = monitor;  if( tmp != fwrite(record, sizeof(gpu_clock), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }
#endif//USE_HDF5_FORMAT


  __NOTE__("%s\n", "end");
}

#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
static inline void appendGPUclockInfoParallel
(const int monitor_step, gpu_clock *deviceMonitors, MPIcfg_dataio mpi,
#ifdef  USE_HDF5_FORMAT
 const hid_t target, const hdf5struct type
#else///USE_HDF5_FORMAT
 MPI_File fh, MPI_Offset *disp
#endif//USE_HDF5_FORMAT
 )
{
  __NOTE__("%s\n", "start");


  const int monitor = (CLOCK_RECORD_STEPS < monitor_step) ? CLOCK_RECORD_STEPS : monitor_step;

  if( monitor > 0 ){

    /* data reordering */
    static gpu_clock record[CLOCK_RECORD_STEPS];
    const int origin = monitor_step & (CLOCK_RECORD_STEPS - 1);
    const int turnaround = monitor - origin;
    for(int ii = 0; ii < turnaround; ii++)
      record[ii] = deviceMonitors[origin + ii];
    for(int ii = turnaround; ii < monitor; ii++)
      record[ii] = deviceMonitors[ii - turnaround];

#ifdef  USE_HDF5_FORMAT
    hid_t group = H5Gcreate(target, "GPUinfo", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* write attribute */
    hsize_t attr_dims = 1;
    hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute = H5Acreate(group, "steps", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &monitor));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Sclose(dataspace));

    /* create (distributed) dataset */
    /* create dataspace */
    hsize_t dims_ful[2] = {mpi.size, monitor};  hid_t fulSpace = H5Screate_simple(2, dims_ful, NULL);
    hsize_t dims_mem[2] = {       1, monitor};
    hsize_t dims_loc[2] = {       1, monitor};  hid_t locSpace = H5Screate_simple(2, dims_loc, NULL);
    /* create chunked dataset */
    hid_t data_create = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t dims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
    if( dims_mem[0] * dims_mem[1] > dims_max )
      dims_mem[0] = dims_max / dims_mem[1];
    chkHDF5err(H5Pset_chunk(data_create, 2, dims_mem));

    /* configuration about distributed dataset */
    hsize_t  count[2] = {1, 1};
    hsize_t stride[2] = {1, 1};
    hsize_t  block[2] = {dims_loc[0], dims_loc[1]};
    hsize_t offset[2] = {mpi.rank, 0};

    /* set up the collective transfer properties function */
    hid_t w_property = H5Pcreate(H5P_DATASET_XFER);
    chkHDF5err(H5Pset_dxpl_mpio(w_property, H5FD_MPIO_COLLECTIVE));

    /* write measured data */
    hid_t dataset = H5Dcreate(group, "record", type.gpu_clock, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hid_t hyperslab = H5Dget_space(dataset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
    chkHDF5err(H5Dwrite(dataset, type.gpu_clock, locSpace, hyperslab, w_property, record));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dataset));

    chkHDF5err(H5Pclose(w_property));
    chkHDF5err(H5Pclose(data_create));
    chkHDF5err(H5Sclose(locSpace));
    chkHDF5err(H5Sclose(fulSpace));

    chkHDF5err(H5Gclose(group));
#else///USE_HDF5_FORMAT
    MPI_Status status;

    /* the root process writes steps, an unsigned int value */
    chkMPIerr(MPI_File_set_view(fh, *disp, MPI_INT, MPI_INT, "native", MPI_INFO_NULL));
    if( mpi.rank == 0 )
      chkMPIerr(MPI_File_write(fh, &monitor, 1, MPI_INT, &status));
    chkMPIerr(MPI_File_sync(fh));
    (*disp) += 1 * (MPI_Offset)sizeof(int);

    /* the whole processes write measured clock frequency of GPUs */
    chkMPIerr(MPI_File_set_view(fh, (*disp) + mpi.rank * (MPI_Offset)monitor * (MPI_Offset)sizeof(gpu_clock), MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL));
    chkMPIerr(MPI_File_write(fh, record, monitor * sizeof(gpu_clock), MPI_BYTE, &status));
    chkMPIerr(MPI_File_sync(fh));
    (*disp) += mpi.size * (MPI_Offset)monitor * (MPI_Offset)sizeof(gpu_clock);
#endif//USE_HDF5_FORMAT

  }/* if( monitor > 0 ){ */

  __NOTE__("%s\n", "end");
}
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)


/**
 * @fn readTentativeData
 *
 * @brief Read data file of N-body particles.
 *
 * @return (time) current time of the simulation
 * @return (dt) current time step of the simulation
 * @return (steps) number of time steps calculated in the simulation
 * @param (num) number of N-body particles
 * @return (body) the particle data
 * @param (file) name of the simulation
 * @param (last) index of the latest restarter file
 * @param (type) data type for HDF5 (only for HDF5 enabled runs)
 * @return (dropPrevTune) a flag to ignore the settings of auto-tuning in the previous run (only for GOTHIC with HDF5)
 * @return (rebuild) tree rebuild information in the previous run (only for GOTHIC with HDF5)
 * @return (measured) measured execution time in the previous run (only for GOTHIC with HDF5)
 * @return (rebuildParam) parameters for auto-tuning of tree rebuild interval in the previous run (only for GOTHIC with HDF5)
 * @return (status) parameters for auto-tuning based on Brent method in the previous run (only for GOTHIC with HDF5)
 * @return (memory) parameters for auto-tuning based on Brent method in the previous run (only for GOTHIC with HDF5)
 * @return (relEneErr) energy error in the previous run (only for GOTHIC with HDF5)
 */
void  readTentativeData(double *time, double *dt, ulong *steps, double *elapsed, int num, iparticle body, char file[], int  last
#ifdef  USE_HDF5_FORMAT
			, hdf5struct type
#ifndef RUN_WITHOUT_GOTHIC
			, int *dropPrevTune, rebuildTree *rebuild, measuredTime *measured, autoTuningParam *rebuildParam, brentStatus *status, brentMemory *memory
#ifdef  MONITOR_ENERGY_ERROR
			, energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
#ifdef  ONLINE_ANALYSIS
			, real * restrict score_all
#endif//ONLINE_ANALYSIS
			)
{
  __NOTE__("%s\n", "start");


  /* open an existing file with read only option */
  char filename[128];
#ifndef USE_HDF5_FORMAT
  FILE *fp;
  sprintf(filename, "%s/%s.%s%d.dat", DATAFOLDER, file, TENTATIVE, last);
  fp = fopen(filename, "rb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fread( time, sizeof(double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(   dt, sizeof(double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(steps, sizeof( ulong), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(elapsed, sizeof(double), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.pos, sizeof(    position), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.acc, sizeof(acceleration), tmp, fp) )    success = false;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  tmp = num;  if( tmp != fread(body.acc_ext, sizeof(acceleration), tmp, fp) )    success = false;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  tmp = num;  if( tmp != fread(body. vel, sizeof(  velocity), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.time, sizeof(ibody_time), tmp, fp) )    success = false;
#else///BLOCK_TIME_STEP
  tmp = num;  if( tmp != fread(body.vx, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.vy, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.vz, sizeof(real), tmp, fp) )    success = false;
#endif//BLOCK_TIME_STEP
  tmp = num;  if( tmp != fread(body.idx, sizeof(ulong), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }

  fclose(fp);
#else///USE_HDF5_FORMAT
  sprintf(filename, "%s/%s.%s%d.h5", DATAFOLDER, file, TENTATIVE, last);
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  /* open an existing group */
  hid_t group = H5Gopen(target, "nbody", H5P_DEFAULT);


  /* read attribute data */
  hid_t attribute;
  /* read current time */
  attribute = H5Aopen(group, "time", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, time));
  chkHDF5err(H5Aclose(attribute));
  /* read time step */
  attribute = H5Aopen(group, "dt", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, dt));
  chkHDF5err(H5Aclose(attribute));
  /* read # of steps */
  attribute = H5Aopen(group, "steps", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, steps));
  chkHDF5err(H5Aclose(attribute));
  /* read elapsed time */
  attribute = H5Aopen(group, "elapsed", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, elapsed));
  chkHDF5err(H5Aclose(attribute));
  /* read # of N-body particles */
  ulong num_ulong = 0;
  attribute = H5Aopen(group, "number", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, &num_ulong));
  chkHDF5err(H5Aclose(attribute));
  /* read inverse of total energy at the initial condition */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Aopen(group, "E0inv", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &relEneErr->E0inv));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  /* read maximum value of the relative error of the total energy */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Aopen(group, "errMax", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &relEneErr->errMax));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  /* read flag about USE_DOUBLE_PRECISION */
  int useDP;
  attribute = H5Aopen(group, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* read flag about BLOCK_TIME_STEP */
  int blockTimeStep;
  attribute = H5Aopen(group, "blockTimeStep", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &blockTimeStep));
  chkHDF5err(H5Aclose(attribute));


  /* simple error checks */
  if( num_ulong != (ulong)num ){    __KILL__(stderr, "ERROR: number of N-body particles in the file (%zu) differs with that in the code (%d)\n", num_ulong, num);  }
#ifdef  USE_DOUBLE_PRECISION
  if( useDP != 1 ){    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);  }
#else///USE_DOUBLE_PRECISION
  if( useDP != 0 ){    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);  }
#endif//USE_DOUBLE_PRECISION
#ifdef  BLOCK_TIME_STEP
  if( blockTimeStep != 1 ){    __KILL__(stderr, "ERROR: blockTimeStep (%d) differs with that in the code (%d)\n", blockTimeStep, true);  }
#else///BLOCK_TIME_STEP
  if( blockTimeStep != 0 ){    __KILL__(stderr, "ERROR: blockTimeStep (%d) differs with that in the code (%d)\n", blockTimeStep, false);  }
#endif//BLOCK_TIME_STEP


  /* read particle data */
  hid_t dataset;
  /* read particle position */
  dataset = H5Dopen(group, "position", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.pos));
  chkHDF5err(H5Dclose(dataset));
  /* read particle acceleration */
  dataset = H5Dopen(group, "acceleration", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.acc));
  chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* read particle acceleration by external potential field */
  dataset = H5Dopen(group, "acceleration_external", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.acc_ext));
  chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  /* read particle velocity */
  dataset = H5Dopen(group, "velocity", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.vel));
  chkHDF5err(H5Dclose(dataset));
  /* read particle time */
  dataset = H5Dopen(group, "time", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.time));
  chkHDF5err(H5Dclose(dataset));
#else///BLOCK_TIME_STEP
  /* read particle velocity */
  dataset = H5Dopen(group, "vx", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.vx));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "vy", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.vy));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "vz", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.vz));
  chkHDF5err(H5Dclose(dataset));
#endif//BLOCK_TIME_STEP
  /* read particle index */
  dataset = H5Dopen(group, "index", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.idx));
  chkHDF5err(H5Dclose(dataset));
  /* end access to the dataset and release the corresponding resource */
  chkHDF5err(H5Gclose(group));


#ifndef  RUN_WITHOUT_GOTHIC
  /* read parameters for auto-tuning in GOTHIC */
  if( (*steps != 0) && (*dropPrevTune == 0) ){
    group = H5Gopen(target, "parameters in auto-tuning", H5P_DEFAULT);

    /* read number of MPI processes in the previous run */
    int procs;
    attribute = H5Aopen(group, "MPI_procs", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &procs));
    chkHDF5err(H5Aclose(attribute));

    if( procs != 1 ){
      *dropPrevTune = 1;
      chkHDF5err(H5Gclose(group));
    }/* if( procs != 1 ){ */
  }/* if( (*steps != 0) && (*dropPrevTune == 0) ){ */

  if( (*steps != 0) && (*dropPrevTune == 0) ){
    /* read attributes */
    /* read flag about FORCE_ADJUSTING_PARTICLE_TIME_STEPS */
    int forceAdjust;
    attribute = H5Aopen(group, "FORCE_ADJUSTING_PARTICLE_TIME_STEPS", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &forceAdjust));
    chkHDF5err(H5Aclose(attribute));
    /* read flag about USE_PARABOLIC_GROWTH_MODEL */
    int parabolic;
    attribute = H5Aopen(group, "USE_PARABOLIC_GROWTH_MODEL", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &parabolic));
    chkHDF5err(H5Aclose(attribute));

    /* read parameters */
    /* read rebuildTree */
#ifdef  FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    if( forceAdjust == 1 )
#else///FORCE_ADJUSTING_PARTICLE_TIME_STEPS
      if( forceAdjust == 0 )
#endif//FORCE_ADJUSTING_PARTICLE_TIME_STEPS
	{
	  dataset = H5Dopen(group, "rebuild tree", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, type.rebuildTree, H5S_ALL, H5S_ALL, H5P_DEFAULT, rebuild));
	  chkHDF5err(H5Dclose(dataset));
	}
    /* read measuredTime */
    dataset = H5Dopen(group, "measured time", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, type.measuredTime, H5S_ALL, H5S_ALL, H5P_DEFAULT, measured));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
    if( parabolic == 1 )
#else///USE_PARABOLIC_GROWTH_MODEL
      if( parabolic == 0 )
#endif//USE_PARABOLIC_GROWTH_MODEL
	{
	  /* read statVal for linear growth model */
	  dataset = H5Dopen(group, "stats (linear)", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, type.statVal, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam->linearStats)));
	  chkHDF5err(H5Dclose(dataset));
	  /* read guessTime for linear growth model */
	  dataset = H5Dopen(group, "guess (linear)", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, type.guessTime, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam->linearGuess)));
	  chkHDF5err(H5Dclose(dataset));
	  /* read statVal for power-law growth model */
	  dataset = H5Dopen(group, "stats (power)", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, type.statVal, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam->powerStats)));
	  chkHDF5err(H5Dclose(dataset));
	  /* read guessTime for power-law growth model */
	  dataset = H5Dopen(group, "guess (power)", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, type.guessTime, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam->powerGuess)));
	  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
	  /* read statVal for parabolic growth model */
	  dataset = H5Dopen(group, "stats (parabolic)", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, type.statVal, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam->parabolicStats)));
	  chkHDF5err(H5Dclose(dataset));
	  /* read guessTime for parabolic growth model */
	  dataset = H5Dopen(group, "guess (parabolic)", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, type.guessTime, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam->parabolicGuess)));
	  chkHDF5err(H5Dclose(dataset));
#endif//USE_PARABOLIC_GROWTH_MODEL
	}
    /* read brentStatus */
    dataset = H5Dopen(group, "Brent status", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, type.brentStatus, H5S_ALL, H5S_ALL, H5P_DEFAULT, status));
    chkHDF5err(H5Dclose(dataset));
    /* read brentMemory */
    dataset = H5Dopen(group, "Brent memory", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, type.brentMemory, H5S_ALL, H5S_ALL, H5P_DEFAULT, memory));
    chkHDF5err(H5Dclose(dataset));

    chkHDF5err(H5Gclose(group));
  }/* if( (*steps != 0) && (*dropPrevTune == 0) ){ */
  else
    *dropPrevTune = 1;
#endif//RUN_WITHOUT_GOTHIC


#ifdef  ONLINE_ANALYSIS
  if( H5Lexists(target, "NWstream", H5P_DEFAULT) ){
    group = H5Gopen(target, "NWstream", H5P_DEFAULT);

    dataset = H5Dopen(group, "score", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, score_all));
    chkHDF5err(H5Dclose(dataset));

    chkHDF5err(H5Gclose(group));
  }
#endif//ONLINE_ANALYSIS


  /* close the file */
  chkHDF5err(H5Fclose(target));
#endif//USE_HDF5_FORMAT

  __NOTE__("%s\n", "end");
}
/**
 * @fn writeTentativeData
 *
 * @brief Write data file of N-body particles.
 *
 * @param (time) current time of the simulation
 * @param (dt) current time step of the simulation
 * @param (steps) number of time steps calculated in the simulation
 * @param (num) number of N-body particles
 * @param (body) the particle data
 * @param (file) name of the simulation
 * @return (last) index of the latest restarter file
 * @param (type) data type for HDF5 (only for HDF5 enabled runs)
 * @param (rebuild) tree rebuild information in the previous run (only for GOTHIC with HDF5)
 * @param (measured) measured execution time in the previous run (only for GOTHIC with HDF5)
 * @param (rebuildParam) parameters for auto-tuning of tree rebuild interval in the previous run (only for GOTHIC with HDF5)
 * @param (status) parameters for auto-tuning based on Brent method in the previous run (only for GOTHIC with HDF5)
 * @param (memory) parameters for auto-tuning based on Brent method in the previous run (only for GOTHIC with HDF5)
 * @param (relEneErr) energy error in the previous run (only for GOTHIC with HDF5)
 */
void writeTentativeData
(double  time, double  dt, ulong  steps, double elapsed, ulong num, iparticle body, char file[], int *last
#ifdef  USE_HDF5_FORMAT
 , hdf5struct type
#ifndef RUN_WITHOUT_GOTHIC
 , rebuildTree rebuild, measuredTime measured, autoTuningParam rebuildParam, brentStatus status, brentMemory memory
#ifdef  MONITOR_ENERGY_ERROR
 , energyError relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
 , const bool dumpGPUclock, gpu_clock *deviceMonitors, const int monitor_step
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
#ifdef  ONLINE_ANALYSIS
 , const bool dumpNWstream, const int Nscore, real * restrict score_all
#endif//ONLINE_ANALYSIS
 )
{
  __NOTE__("%s\n", "start");


  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
#ifndef USE_HDF5_FORMAT
  FILE *fp;
  sprintf(filename, "%s/%s.%s%d.dat", DATAFOLDER, file, TENTATIVE, (*last) ^ 1);
  fp = fopen(filename, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fwrite(& time, sizeof(        double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(&   dt, sizeof(        double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(&steps, sizeof(         ulong), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(&elapsed, sizeof(      double), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.pos, sizeof(    position), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.acc, sizeof(acceleration), tmp, fp) )    success = false;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  tmp = num;  if( tmp != fwrite(body.acc_ext, sizeof(acceleration), tmp, fp) )    success = false;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body. vel, sizeof(  velocity), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.time, sizeof(ibody_time), tmp, fp) )    success = false;
#else///BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body.vx, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.vy, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.vz, sizeof(real), tmp, fp) )    success = false;
#endif//BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body.idx, sizeof(ulong), tmp, fp) )    success = false;

#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
  if( dumpGPUclock )
    appendGPUclockInfo(monitor_step, deviceMonitors, filename, fp);
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)

  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }
  fclose(fp);
  *last ^= 1;
#else///USE_HDF5_FORMAT
  sprintf(filename, "%s/%s.%s%d.h5", DATAFOLDER, file, TENTATIVE, (*last) ^ 1);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t group = H5Gcreate(target, "nbody", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


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

  /* write num * real4 arrays */
  hsize_t dims = num * 4;
  dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  hsize_t cdims_loc = cdims;
  if( dims < cdims_loc )
    cdims_loc = (hsize_t)1 << llog2((ulong)dims);
  hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  if( cdims_loc > cdims_max )
    cdims_loc = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
  /* write particle position */
  dataset = H5Dcreate(group, "position", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.pos));
  chkHDF5err(H5Dclose(dataset));
  /* write particle acceleration */
  dataset = H5Dcreate(group, "acceleration", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.acc));
  chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* write particle acceleration by external potential field */
  dataset = H5Dcreate(group, "acceleration_external", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.acc_ext));
  chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  /* write particle velocity */
  dataset = H5Dcreate(group, "velocity", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.vel));
  chkHDF5err(H5Dclose(dataset));
#endif//BLOCK_TIME_STEP
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));

#ifdef  BLOCK_TIME_STEP
  /* write num * double2 array */
  dims = num * 2;
  dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  cdims_loc = cdims;
  if( dims < cdims_loc )
    cdims_loc = (hsize_t)1 << llog2((ulong)dims);
  cdims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
  if( cdims_loc > cdims_max )
    cdims_loc = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
  /* write particle time */
  dataset = H5Dcreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.time));
  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));
#endif//BLOCK_TIME_STEP


  /* write num array */
  dims = num;
  dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  cdims_loc = cdims;
  if( dims < cdims_loc )
    cdims_loc = (hsize_t)1 << llog2((ulong)dims);
  cdims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
  if( cdims_loc > cdims_max )
    cdims_loc = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
#ifndef BLOCK_TIME_STEP
  /* write particle velocity */
  dataset = H5Dcreate(group, "vx", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.vx));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "vy", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.vy));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "vz", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.vz));
  chkHDF5err(H5Dclose(dataset));
#endif//BLOCK_TIME_STEP
  /* write particle index */
  dataset = H5Dcreate(group, "index", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.idx));
  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));


  /* write attribute data */
  /* create the data space for the attribute */
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  /* write current time */
  attribute = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time));
  chkHDF5err(H5Aclose(attribute));
  /* write time step */
  attribute = H5Acreate(group, "dt", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &dt));
  chkHDF5err(H5Aclose(attribute));
  /* write # of steps */
  attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  /* write elapsed time */
  attribute = H5Acreate(group, "elapsed", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &elapsed));
  chkHDF5err(H5Aclose(attribute));
  /* write # of N-body particles */
  ulong num_ulong = (ulong)num;
  attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &num_ulong));
  chkHDF5err(H5Aclose(attribute));
  /* write inverse of total energy at the initial condition */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Acreate(group, "E0inv", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &relEneErr.E0inv));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  /* write maximum value of the relative error of the total energy */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Acreate(group, "errMax", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &relEneErr.errMax));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  /* write flag about USE_DOUBLE_PRECISION */
#ifdef  USE_DOUBLE_PRECISION
  const int useDP = 1;
#else///USE_DOUBLE_PRECISION
  const int useDP = 0;
#endif//USE_DOUBLE_PRECISION
  attribute = H5Acreate(group, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* write flag about BLOCK_TIME_STEP */
#ifdef  BLOCK_TIME_STEP
  const int blockTimeStep = 1;
#else///BLOCK_TIME_STEP
  const int blockTimeStep = 0;
#endif//BLOCK_TIME_STEP
  attribute = H5Acreate(group, "blockTimeStep", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &blockTimeStep));
  chkHDF5err(H5Aclose(attribute));
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(group));


  /* output parameters for auto-tuning (when steps > 0) */
#ifndef  RUN_WITHOUT_GOTHIC
  if( steps != 0 ){
    group = H5Gcreate(target, "parameters in auto-tuning", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* write parameters */
    dims = 1;
    dataspace = H5Screate_simple(1, &dims, NULL);
    /* output rebuildTree */
    dataset = H5Dcreate(group, "rebuild tree", type.rebuildTree, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.rebuildTree, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rebuild));
    chkHDF5err(H5Dclose(dataset));
    /* output measuredTime */
    dataset = H5Dcreate(group, "measured time", type.measuredTime, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.measuredTime, H5S_ALL, H5S_ALL, H5P_DEFAULT, &measured));
    chkHDF5err(H5Dclose(dataset));
    /* output statVal for linear growth model */
    dataset = H5Dcreate(group, "stats (linear)", type.statVal, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.statVal, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam.linearStats)));
    chkHDF5err(H5Dclose(dataset));
    /* output guessTime for linear growth model */
    dataset = H5Dcreate(group, "guess (linear)", type.guessTime, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.guessTime, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam.linearGuess)));
    chkHDF5err(H5Dclose(dataset));
    /* output statVal for power-law growth model */
    dataset = H5Dcreate(group, "stats (power)", type.statVal, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.statVal, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam.powerStats)));
    chkHDF5err(H5Dclose(dataset));
    /* output guessTime for power-law growth model */
    dataset = H5Dcreate(group, "guess (power)", type.guessTime, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.guessTime, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam.powerGuess)));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
    /* output statVal for parabolic growth model */
    dataset = H5Dcreate(group, "stats (parabolic)", type.statVal, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.statVal, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam.parabolicStats)));
    chkHDF5err(H5Dclose(dataset));
    /* output guessTime for parabolic growth model */
    dataset = H5Dcreate(group, "guess (parabolic)", type.guessTime, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.guessTime, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(rebuildParam.parabolicGuess)));
    chkHDF5err(H5Dclose(dataset));
#endif//USE_PARABOLIC_GROWTH_MODEL
    /* output brentStatus */
    dataset = H5Dcreate(group, "Brent status", type.brentStatus, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.brentStatus, H5S_ALL, H5S_ALL, H5P_DEFAULT, &status));
    chkHDF5err(H5Dclose(dataset));
    /* output brentMemory */
    dataset = H5Dcreate(group, "Brent memory", type.brentMemory, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.brentMemory, H5S_ALL, H5S_ALL, H5P_DEFAULT, &memory));
    chkHDF5err(H5Dclose(dataset));
    chkHDF5err(H5Sclose(dataspace));

    /* write attributes */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    int flag;
    /* write number of MPI processes in the current run */
    flag = 1;
    attribute = H5Acreate(group, "MPI_procs", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about FORCE_ADJUSTING_PARTICLE_TIME_STEPS */
#ifdef  FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    flag = 1;
#else///FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    flag = 0;
#endif//FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    attribute = H5Acreate(group, "FORCE_ADJUSTING_PARTICLE_TIME_STEPS", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about USE_PARABOLIC_GROWTH_MODEL */
#ifdef  USE_PARABOLIC_GROWTH_MODEL
    flag = 1;
#else///USE_PARABOLIC_GROWTH_MODEL
    flag = 0;
#endif//USE_PARABOLIC_GROWTH_MODEL
    attribute = H5Acreate(group, "USE_PARABOLIC_GROWTH_MODEL", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
  }/* if( steps != 0 ){ */
#endif//RUN_WITHOUT_GOTHIC

#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
  if( dumpGPUclock )
    appendGPUclockInfo(monitor_step, deviceMonitors, target, type);
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)


#ifdef  ONLINE_ANALYSIS
  if( dumpNWstream ){
    group = H5Gcreate(target, "NWstream", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dims = Nscore;
    dataspace = H5Screate_simple(1, &dims, NULL);
    dataset = H5Dcreate(group, "score", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, score_all));
    chkHDF5err(H5Dclose(dataset));
    chkHDF5err(H5Sclose(dataspace));

    chkHDF5err(H5Gclose(group));
  }
#endif//ONLINE_ANALYSIS


  /* close the file */
  chkHDF5err(H5Fclose(target));
  *last ^= 1;
#endif//USE_HDF5_FORMAT


  __NOTE__("%s\n", "end");
}


#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
/**
 * @fn readTentativeDataParallel
 *
 * @brief Read data file of N-body particles.
 *
 * @return (time) current time of the simulation
 * @return (dt) current time step of the simulation
 * @return (steps) number of time steps calculated in the simulation
 * @return (num) number of N-body particles contained in this MPI process
 * @return (body) the particle data
 * @param (file) name of the simulation
 * @param (last) index of the latest restarter file
 * @return (mpi) MPI communicator for data I/O
 * @param (Ntot) total number of N-body particles
 * @param (type) data type for HDF5 (only for HDF5 enabled runs)
 * @return (dropPrevTune) a flag to ignore the settings of auto-tuning in the previous run (only for GOTHIC with HDF5)
 * @return (rebuild) tree rebuild information in the previous run (only for GOTHIC with HDF5)
 * @return (measured) measured execution time in the previous run (only for GOTHIC with HDF5)
 * @return (rebuildParam) parameters for auto-tuning of tree rebuild interval in the previous run (only for GOTHIC with HDF5)
 * @return (status) parameters for auto-tuning based on Brent method in the previous run (only for GOTHIC with HDF5)
 * @return (memory) parameters for auto-tuning based on Brent method in the previous run (only for GOTHIC with HDF5)
 * @return (relEneErr) energy error in the previous run (only for GOTHIC with HDF5)
 */
void  readTentativeDataParallel(double *time, double *dt, ulong *steps, double *elapsed, int *num, iparticle body, char file[], int  last, MPIcfg_dataio *mpi, ulong Ntot
#ifdef  USE_HDF5_FORMAT
				, hdf5struct type
#ifndef RUN_WITHOUT_GOTHIC
				, int *dropPrevTune, rebuildTree *rebuild, measuredTime *measured, autoTuningParam *rebuildParam, brentStatus *status, brentMemory *memory
#ifdef  MONITOR_ENERGY_ERROR
				, energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
				)
{
  __NOTE__("%s\n", "start");


  /* open an existing file with read only option */
  char filename[128];
#ifndef USE_HDF5_FORMAT
  /* reset MPI_Offset */
  updateMPIcfg_dataio(mpi, *num);

  sprintf(filename, "%s/%s.%s%d.dat", DATAFOLDER, file, TENTATIVE, last);
  MPI_File fh;
  chkMPIerr(MPI_File_open(mpi->comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh));

  MPI_Status status;
  MPI_Offset disp = 0;
  /* the root process reads and broadcasts time, a double value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_read(fh, time, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_Bcast(time, 1, MPI_DOUBLE, 0, mpi->comm));
  disp += 1 * (MPI_Offset)sizeof(double);
  /* the root process reads and broadcasts dt, a double value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_read(fh, dt, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_Bcast(dt, 1, MPI_DOUBLE, 0, mpi->comm));
  disp += 1 * (MPI_Offset)sizeof(double);
  /* the root process reads and broadcasts # of steps, an unsigned long value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_read(fh, steps, 1, MPI_UNSIGNED_LONG, &status));
  chkMPIerr(MPI_Bcast(steps, 1, MPI_UNSIGNED_LONG, 0, mpi->comm));
  disp += 1 * (MPI_Offset)sizeof(ulong);
  /* the root process reads and broadcasts elapsed time, a double value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_read(fh, elapsed, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_Bcast(elapsed, 1, MPI_DOUBLE, 0, mpi->comm));
  disp += 1 * (MPI_Offset)sizeof(double);
  /* the whole processes read position */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(position), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.pos, (*num) * 4, MPI_REALDAT, &status));
  disp += Ntot * (MPI_Offset)sizeof(position);
  /* the whole processes read acceleration */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(acceleration), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.acc, (*num) * 4, MPI_REALDAT, &status));
  disp += Ntot * (MPI_Offset)sizeof(acceleration);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* the whole processes read acceleration by external potential field */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(acceleration), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.acc_ext, (*num) * 4, MPI_REALDAT, &status));
  disp += Ntot * (MPI_Offset)sizeof(acceleration);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  /* the whole processes read velocity and time */
#ifdef  BLOCK_TIME_STEP
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(velocity), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.vel, (*num) * 4, MPI_REALDAT, &status));
  disp += Ntot * (MPI_Offset)sizeof(velocity);
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(ibody_time), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.time, (*num) * 2, MPI_DOUBLE, &status));
  disp += Ntot * (MPI_Offset)sizeof(ibody_time);
#else///BLOCK_TIME_STEP
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(real), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.vx, *num, MPI_REALDAT, &status));
  disp += Ntot * (MPI_Offset)sizeof(real);
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(real), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.vy, *num, MPI_REALDAT, &status));
  disp += Ntot * (MPI_Offset)sizeof(real);
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(real), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.vz, *num, MPI_REALDAT, &status));
  disp += Ntot * (MPI_Offset)sizeof(real);
#endif//BLOCK_TIME_STEP
  /* the whole processes read index */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(ulong), MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.idx, *num, MPI_UNSIGNED_LONG, &status));
  disp += Ntot * (MPI_Offset)sizeof(ulong);

  /* close the target file */
  chkMPIerr(MPI_File_close(&fh));
#else///USE_HDF5_FORMAT
  sprintf(filename, "%s/%s.%s%d.h5", DATAFOLDER, file, TENTATIVE, last);
  hid_t f_property = H5Pcreate(H5P_FILE_ACCESS);
  chkHDF5err(H5Pset_fapl_mpio(f_property, mpi->comm, mpi->info));
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, f_property);
  chkHDF5err(H5Pclose(f_property));
  /* open an existing group */
  hid_t group = H5Gopen(target, "nbody", H5P_DEFAULT);
  /* create property list for collective dataset read */
  hid_t r_property = H5Pcreate(H5P_DATASET_XFER);
  chkHDF5err(H5Pset_dxpl_mpio(r_property, H5FD_MPIO_COLLECTIVE));


  /* read attribute data */
  ulong Nread = 0;
  int useDP;
  int blockTimeStep;
  if( mpi->rank == 0 ){
    hid_t attribute;
    /* read current time */
    attribute = H5Aopen(group, "time", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, time));
    chkHDF5err(H5Aclose(attribute));
    /* read time step */
    attribute = H5Aopen(group, "dt", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, dt));
    chkHDF5err(H5Aclose(attribute));
    /* read # of steps */
    attribute = H5Aopen(group, "steps", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, steps));
    chkHDF5err(H5Aclose(attribute));
    /* read elapsed time */
    attribute = H5Aopen(group, "elapsed", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, elapsed));
    chkHDF5err(H5Aclose(attribute));
    /* read # of N-body particles */
    attribute = H5Aopen(group, "number", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, &Nread));
    chkHDF5err(H5Aclose(attribute));
    /* read inverse of total energy at the initial condition */
#ifdef  MONITOR_ENERGY_ERROR
    attribute = H5Aopen(group, "E0inv", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &relEneErr->E0inv));
    chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
    /* read maximum value of the relative error of the total energy */
#ifdef  MONITOR_ENERGY_ERROR
    attribute = H5Aopen(group, "errMax", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &relEneErr->errMax));
    chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
    /* read flag about USE_DOUBLE_PRECISION */
    attribute = H5Aopen(group, "useDP", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
    chkHDF5err(H5Aclose(attribute));
    /* read flag about BLOCK_TIME_STEP */
    attribute = H5Aopen(group, "blockTimeStep", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &blockTimeStep));
    chkHDF5err(H5Aclose(attribute));
  }/* if( mpi->rank == 0 ){ */

  /* broadcast the read attributes from the root process */
  chkMPIerr(MPI_Bcast(  time, 1, MPI_DOUBLE       , 0, mpi->comm));
  chkMPIerr(MPI_Bcast(    dt, 1, MPI_DOUBLE       , 0, mpi->comm));
  chkMPIerr(MPI_Bcast( steps, 1, MPI_UNSIGNED_LONG, 0, mpi->comm));
  chkMPIerr(MPI_Bcast(elapsed, 1, MPI_DOUBLE, 0, mpi->comm));
  chkMPIerr(MPI_Bcast(&Nread, 1, MPI_UNSIGNED_LONG, 0, mpi->comm));
#ifdef  MONITOR_ENERGY_ERROR
  chkMPIerr(MPI_Bcast(&relEneErr->E0inv , 1, MPI_DOUBLE, 0, mpi->comm));
  chkMPIerr(MPI_Bcast(&relEneErr->errMax, 1, MPI_DOUBLE, 0, mpi->comm));
#endif//MONITOR_ENERGY_ERROR
  chkMPIerr(MPI_Bcast(&useDP, 1, MPI_INT, 0, mpi->comm));
  chkMPIerr(MPI_Bcast(&blockTimeStep, 1, MPI_INT, 0, mpi->comm));


  /* simple error checks */
  if( Nread != Ntot ){
    __KILL__(stderr, "ERROR: number of N-body particles in the file (%zu) differs with that in the code (%zu)\n", Nread, Ntot);
  }/* if( Nread != Ntot ){ */
#ifdef  USE_DOUBLE_PRECISION
  if( useDP != 1 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);
  }/* if( useDP != 1 ){ */
#else///USE_DOUBLE_PRECISION
  if( useDP != 0 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);
  }/* if( useDP != 0 ){ */
#endif//USE_DOUBLE_PRECISION
#ifdef  BLOCK_TIME_STEP
  if( blockTimeStep != 1 ){
    __KILL__(stderr, "ERROR: blockTimeStep (%d) differs with that in the code (%d)\n", blockTimeStep, true);
  }/* if( blockTimeStep != 1 ){ */
#else///BLOCK_TIME_STEP
  if( blockTimeStep != 0 ){
    __KILL__(stderr, "ERROR: blockTimeStep (%d) differs with that in the code (%d)\n", blockTimeStep, false);
  }/* if( blockTimeStep != 0 ){ */
#endif//BLOCK_TIME_STEP


#ifndef  RUN_WITHOUT_GOTHIC
  /* read parameters for auto-tuning in GOTHIC */
  if( (*steps != 0) && (*dropPrevTune == 0) ){
    chkHDF5err(H5Gclose(group));
    group = H5Gopen(target, "parameters in auto-tuning", H5P_DEFAULT);

    /* read number of MPI processes in the previous run */
    int procs;
    if( mpi->rank == 0 ){
      hid_t attribute = H5Aopen(group, "MPI_procs", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &procs));
      chkHDF5err(H5Aclose(attribute));
    }/* if( mpi->rank == 0 ){ */
    chkMPIerr(MPI_Bcast(&procs, 1, MPI_INT, 0, mpi->comm));

    if( procs != mpi->size ){
      *dropPrevTune = 1;
      chkHDF5err(H5Gclose(group));
      group = H5Gopen(target, "nbody", H5P_DEFAULT);
    }/* if( procs != mpi->size ){ */
  }/* if( (*steps != 0) && (*dropPrevTune == 0) ){ */

  if( (*steps != 0) && (*dropPrevTune == 0) ){
    /* read attributes */
    int forceAdjust;
    int monitor, parabolic;
    if( mpi->rank == 0 ){
      hid_t attribute;
      /* read flag about FORCE_ADJUSTING_PARTICLE_TIME_STEPS */
      attribute = H5Aopen(group, "FORCE_ADJUSTING_PARTICLE_TIME_STEPS", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &forceAdjust));
      chkHDF5err(H5Aclose(attribute));
      /* read flag about MONITOR_LETGEN_TIME */
      attribute = H5Aopen(group, "MONITOR_LETGEN_TIME", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &monitor));
      chkHDF5err(H5Aclose(attribute));
      /* read flag about USE_PARABOLIC_GROWTH_MODEL */
      attribute = H5Aopen(group, "USE_PARABOLIC_GROWTH_MODEL", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &parabolic));
      chkHDF5err(H5Aclose(attribute));
    }/* if( mpi->rank == 0 ){ */

    /* broadcast the read attributes from the root process */
    chkMPIerr(MPI_Bcast(&forceAdjust, 1, MPI_INT, 0, mpi->comm));
    chkMPIerr(MPI_Bcast(&    monitor, 1, MPI_INT, 0, mpi->comm));
    chkMPIerr(MPI_Bcast(&  parabolic, 1, MPI_INT, 0, mpi->comm));


    /* read parameters */
    hsize_t dims_loc = 1;
    hid_t locSpace = H5Screate_simple(1, &dims_loc, NULL);
    /* configuration about domain decomposition */
    hsize_t  count = 1;
    hsize_t stride = 1;
    hsize_t  block = dims_loc;
    hsize_t offset = mpi->rank;
    hid_t dataset, hyperslab;

    /* read rebuildTree */
#ifdef  FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    if( forceAdjust == 1 )
#else///FORCE_ADJUSTING_PARTICLE_TIME_STEPS
      if( forceAdjust == 0 )
#endif//FORCE_ADJUSTING_PARTICLE_TIME_STEPS
	{
	  dataset = H5Dopen(group, "rebuild tree", H5P_DEFAULT);
	  hyperslab = H5Dget_space(dataset);
	  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
	  chkHDF5err(H5Dread(dataset, type.rebuildTree, locSpace, hyperslab, r_property, rebuild));
	  chkHDF5err(H5Sclose(hyperslab));
	  chkHDF5err(H5Dclose(dataset));
	}
    /* read measuredTime */
#ifdef  MONITOR_LETGEN_TIME
    if( monitor == 1 )
#else///MONITOR_LETGEN_TIME
      if( monitor == 0 )
#endif//MONITOR_LETGEN_TIME
	{
	  dataset = H5Dopen(group, "measured time", H5P_DEFAULT);
	  hyperslab = H5Dget_space(dataset);
	  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
	  chkHDF5err(H5Dread(dataset, type.measuredTime, locSpace, hyperslab, r_property, measured));
	  chkHDF5err(H5Sclose(hyperslab));
	  chkHDF5err(H5Dclose(dataset));
	}
#ifdef  USE_PARABOLIC_GROWTH_MODEL
    if( parabolic == 1 )
#else///USE_PARABOLIC_GROWTH_MODEL
      if( parabolic == 0 )
#endif//USE_PARABOLIC_GROWTH_MODEL
	{
	  /* read statVal for linear growth model */
	  dataset = H5Dopen(group, "stats (linear)", H5P_DEFAULT);
	  hyperslab = H5Dget_space(dataset);
	  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
	  chkHDF5err(H5Dread(dataset, type.statVal, locSpace, hyperslab, r_property, &(rebuildParam->linearStats)));
	  chkHDF5err(H5Sclose(hyperslab));
	  chkHDF5err(H5Dclose(dataset));
	  /* read guessTime for linear growth model */
	  dataset = H5Dopen(group, "guess (linear)", H5P_DEFAULT);
	  hyperslab = H5Dget_space(dataset);
	  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
	  chkHDF5err(H5Dread(dataset, type.guessTime, locSpace, hyperslab, r_property, &(rebuildParam->linearGuess)));
	  chkHDF5err(H5Sclose(hyperslab));
	  chkHDF5err(H5Dclose(dataset));
	  /* read statVal for power-law growth model */
	  dataset = H5Dopen(group, "stats (power)", H5P_DEFAULT);
	  hyperslab = H5Dget_space(dataset);
	  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
	  chkHDF5err(H5Dread(dataset, type.statVal, locSpace, hyperslab, r_property, &(rebuildParam->powerStats)));
	  chkHDF5err(H5Sclose(hyperslab));
	  chkHDF5err(H5Dclose(dataset));
	  /* read guessTime for power-law growth model */
	  dataset = H5Dopen(group, "guess (power)", H5P_DEFAULT);
	  hyperslab = H5Dget_space(dataset);
	  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
	  chkHDF5err(H5Dread(dataset, type.guessTime, locSpace, hyperslab, r_property, &(rebuildParam->powerGuess)));
	  chkHDF5err(H5Sclose(hyperslab));
	  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
	  /* read statVal for parabolic growth model */
	  dataset = H5Dopen(group, "stats (parabolic)", H5P_DEFAULT);
	  hyperslab = H5Dget_space(dataset);
	  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
	  chkHDF5err(H5Dread(dataset, type.statVal, locSpace, hyperslab, r_property, &(rebuildParam->parabolicStats)));
	  chkHDF5err(H5Sclose(hyperslab));
	  chkHDF5err(H5Dclose(dataset));
	  /* read guessTime for parabolic growth model */
	  dataset = H5Dopen(group, "guess (parabolic)", H5P_DEFAULT);
	  hyperslab = H5Dget_space(dataset);
	  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
	  chkHDF5err(H5Dread(dataset, type.guessTime, locSpace, hyperslab, r_property, &(rebuildParam->parabolicGuess)));
	  chkHDF5err(H5Sclose(hyperslab));
	  chkHDF5err(H5Dclose(dataset));
#endif//USE_PARABOLIC_GROWTH_MODEL
	}

    /* read brentStatus */
    dataset = H5Dopen(group, "Brent status", H5P_DEFAULT);
    hyperslab = H5Dget_space(dataset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dread(dataset, type.brentStatus, locSpace, hyperslab, r_property, status));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dataset));
    /* read brentMemory */
    dataset = H5Dopen(group, "Brent memory", H5P_DEFAULT);
    hyperslab = H5Dget_space(dataset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dread(dataset, type.brentMemory, locSpace, hyperslab, r_property, memory));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dataset));

    /* read # of particles in each process */
    dataset = H5Dopen(group, "num", H5P_DEFAULT);
    hyperslab = H5Dget_space(dataset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dread(dataset, H5T_NATIVE_INT, locSpace, hyperslab, r_property, num));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dataset));
    chkHDF5err(H5Sclose(locSpace));

    chkHDF5err(H5Gclose(group));
    group = H5Gopen(target, "nbody", H5P_DEFAULT);
  }/* if( (*steps != 0) && (*dropPrevTune == 0) ){ */
  else
    *dropPrevTune = 1;
#endif//RUN_WITHOUT_GOTHIC


  /* reset MPI_Offset */
  updateMPIcfg_dataio(mpi, *num);


  /* read particle data */
  /* dataset for real4 arrays */
  /* dataspace */
  hsize_t dims_loc = 4 * (*num);
  hid_t locSpace = H5Screate_simple(1, &dims_loc, NULL);
  /* configuration about domain decomposition */
  hsize_t  count = 1;
  hsize_t stride = 1;
  hsize_t  block = dims_loc;
  hsize_t offset = mpi->head * 4;

  hid_t dataset, hyperslab;
  /* read particle position */
  dataset = H5Dopen(group, "position", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, type.real, locSpace, hyperslab, r_property, body.pos));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* read particle acceleration */
  dataset = H5Dopen(group, "acceleration", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, type.real, locSpace, hyperslab, r_property, body.acc));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* read particle acceleration by external potential field */
  dataset = H5Dopen(group, "acceleration_external", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, type.real, locSpace, hyperslab, r_property, body.acc_ext));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  dataset = H5Dopen(group, "velocity", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, type.real, locSpace, hyperslab, r_property, body.vel));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
#endif//BLOCK_TIME_STEP
  chkHDF5err(H5Sclose(locSpace));

#ifdef  BLOCK_TIME_STEP
  /* dataset for double2 array */
  /* dataspace */
  dims_loc = 2 * (*num);
  locSpace = H5Screate_simple(1, &dims_loc, NULL);
  /* configuration about domain decomposition */
  count  = 1;
  stride = 1;
  block  = dims_loc;
  offset = mpi->head * 2;

  /* read particle time */
  dataset = H5Dopen(group, "time", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, locSpace, hyperslab, r_property, body.time));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(locSpace));
#endif//BLOCK_TIME_STEP


  /* dataset for num array */
  /* dataspace */
  dims_loc = *num;
  locSpace = H5Screate_simple(1, &dims_loc, NULL);
  /* configuration about domain decomposition */
  count  = 1;
  stride = 1;
  block  = dims_loc;
  offset = mpi->head;

#ifndef BLOCK_TIME_STEP
  /* read particle velocity */
  dataset = H5Dopen(group, "vx", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, type.real, locSpace, hyperslab, r_property, body.vx));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "vy", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, type.real, locSpace, hyperslab, r_property, body.vy));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "vz", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, type.real, locSpace, hyperslab, r_property, body.vz));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
#endif//BLOCK_TIME_STEP
  dataset = H5Dopen(group, "index", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_ULONG, locSpace, hyperslab, r_property, body.idx));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(locSpace));
  /* finish collective dataset read */
  chkHDF5err(H5Pclose(r_property));


  /* close the file */
  /* end access to the dataset and release the corresponding resource */
  chkHDF5err(H5Gclose(group));
  chkHDF5err(H5Fclose(target));
#endif//USE_HDF5_FORMAT


  __NOTE__("%s\n", "end");
}
/**
 * @fn writeTentativeDataParallel
 *
 * @brief Write data file of N-body particles.
 *
 * @param (time) current time of the simulation
 * @param (dt) current time step of the simulation
 * @param (steps) number of time steps calculated in the simulation
 * @param (num) number of N-body particles contained in this MPI process
 * @param (body) the particle data
 * @param (file) name of the simulation
 * @return (last) index of the latest restarter file
 * @return (mpi) MPI communicator for data I/O
 * @param (Ntot) total number of N-body particles
 * @param (type) data type for HDF5 (only for HDF5 enabled runs)
 * @param (rebuild) tree rebuild information in the previous run (only for GOTHIC with HDF5)
 * @param (measured) measured execution time in the previous run (only for GOTHIC with HDF5)
 * @param (rebuildParam) parameters for auto-tuning of tree rebuild interval in the previous run (only for GOTHIC with HDF5)
 * @param (status) parameters for auto-tuning based on Brent method in the previous run (only for GOTHIC with HDF5)
 * @param (memory) parameters for auto-tuning based on Brent method in the previous run (only for GOTHIC with HDF5)
 * @param (relEneErr) energy error in the previous run (only for GOTHIC with HDF5)
 */
void writeTentativeDataParallel
(double  time, double  dt, ulong  steps, double elapsed, int num, iparticle body, char file[], int *last, MPIcfg_dataio *mpi, ulong Ntot
#ifdef  USE_HDF5_FORMAT
 , hdf5struct type
#ifndef RUN_WITHOUT_GOTHIC
 , rebuildTree rebuild, measuredTime measured, autoTuningParam rebuildParam, brentStatus status, brentMemory memory
#ifdef  MONITOR_ENERGY_ERROR
 , energyError relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
 , const bool dumpGPUclock, gpu_clock *deviceMonitors, const int monitor_step
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
 )
{
  __NOTE__("%s\n", "start");


  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
#ifndef USE_HDF5_FORMAT
  /* reset MPI_Offset */
  updateMPIcfg_dataio(mpi, num);

  /* open the target file */
  sprintf(filename, "%s/%s.%s%d.dat", DATAFOLDER, file, TENTATIVE, (*last) ^ 1);
  MPI_File fh;
  chkMPIerr(MPI_File_open(mpi->comm, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh));
  chkMPIerr(MPI_File_sync(fh));
  chkMPIerr(MPI_File_set_size(fh, 0));
  chkMPIerr(MPI_File_sync(fh));

  MPI_Status status;
  MPI_Offset disp = 0;
  /* the root process writes time, a double value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &time, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(double);
  /* the root process writes dt, a double value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &dt, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(double);
  /* the root process writes steps, an unsigned long value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &steps, 1, MPI_UNSIGNED_LONG, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(ulong);
  /* the root process writes dt, a double value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &elapsed, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(double);
  /* the whole processes write position */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(position), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.pos, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(position);
  /* the whole processes write acceleration */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(acceleration), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.acc, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(acceleration);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* the whole processes write acceleration by external potential field */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(acceleration), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.acc_ext, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(acceleration);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  /* the whole processes write velocity and time */
#ifdef  BLOCK_TIME_STEP
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(velocity), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.vel, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(velocity);
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(ibody_time), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.time, num * 2, MPI_DOUBLE, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(ibody_time);
#else///BLOCK_TIME_STEP
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(real), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.vx, num, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(real);
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(real), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.vy, num, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(real);
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(real), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.vz, num, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(real);
#endif//BLOCK_TIME_STEP
  /* the whole processes write index */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(ulong), MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.idx, num, MPI_UNSIGNED_LONG, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(ulong);

#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
  if( dumpGPUclock )
    appendGPUclockInfoParallel(monitor_step, deviceMonitors, *mpi, fh, &disp);
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)

  /* close the output file */
  chkMPIerr(MPI_File_close(&fh));
  *last ^= 1;
#else///USE_HDF5_FORMAT
  hid_t f_property = H5Pcreate(H5P_FILE_ACCESS);
  chkHDF5err(H5Pset_fapl_mpio(f_property, mpi->comm, mpi->info));
  sprintf(filename, "%s/%s.%s%d.h5", DATAFOLDER, file, TENTATIVE, (*last) ^ 1);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, f_property);
  chkHDF5err(H5Pclose(f_property));
  hid_t group = H5Gcreate(target, "nbody", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create property list for collective dataset write */
  hid_t w_property = H5Pcreate(H5P_DATASET_XFER);
  chkHDF5err(H5Pset_dxpl_mpio(w_property, H5FD_MPIO_COLLECTIVE));
  /* configuration about domain decomposition */
  updateMPIcfg_dataio(mpi, num);


  /* create (distributed) dataset for real4 arrays */
  /* create dataspace */
  hsize_t dims_ful = 4 * Ntot            ;  hid_t fulSpace = H5Screate_simple(1, &dims_ful, NULL);
  hsize_t dims_loc = 4 *  num            ;  hid_t locSpace = H5Screate_simple(1, &dims_loc, NULL);
  hsize_t dims_mem = 4 * Ntot / mpi->size;
  /* configuration about domain decomposition */
  hsize_t  count = 1;
  hsize_t stride = 1;
  hsize_t  block = dims_loc;
  hsize_t offset = mpi->head * 4;
  /* create chunked dataset */
  hid_t data_create = H5Pcreate(H5P_DATASET_CREATE);
  hsize_t dims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  if( dims_mem > dims_max )
    dims_mem = dims_max;
  chkHDF5err(H5Pset_chunk(data_create, 1, &dims_mem));
  hid_t dset, hyperslab;

  /* write particle position */
  dset = H5Dcreate(group, "position", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, type.real, locSpace, hyperslab, w_property, body.pos));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));
  /* write particle acceleration */
  dset = H5Dcreate(group, "acceleration", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, type.real, locSpace, hyperslab, w_property, body.acc));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* write particle acceleration by external potential field */
  dset = H5Dcreate(group, "acceleration_external", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, type.real, locSpace, hyperslab, w_property, body.acc_ext));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  /* write particle velocity */
  dset = H5Dcreate(group, "velocity", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, type.real, locSpace, hyperslab, w_property, body.vel));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));
#endif//BLOCK_TIME_STEP

  /* close/release resources */
  chkHDF5err(H5Pclose(data_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));


#ifdef  BLOCK_TIME_STEP
  /* create (distributed) dataset for double2 array */
  /* create dataspace */
  dims_ful = 2 * Ntot            ;  fulSpace = H5Screate_simple(1, &dims_ful, NULL);
  dims_mem = 2 * Ntot / mpi->size;
  dims_loc = 2 *  num            ;  locSpace = H5Screate_simple(1, &dims_loc, NULL);
  /* configuration about domain decomposition */
  count  = 1;
  stride = 1;
  block  = dims_loc;
  offset = mpi->head * 2;
  /* create chunked dataset */
  data_create = H5Pcreate(H5P_DATASET_CREATE);
  dims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
  if( dims_mem > dims_max )
    dims_mem = dims_max;
  chkHDF5err(H5Pset_chunk(data_create, 1, &dims_mem));

  /* write particle time */
  dset = H5Dcreate(group, "time", H5T_NATIVE_DOUBLE, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, H5T_NATIVE_DOUBLE, locSpace, hyperslab, w_property, body.time));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));

  /* close/release resources */
  chkHDF5err(H5Pclose(data_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));
#endif//BLOCK_TIME_STEP


  /* create (distributed) dataset for num arrays */
  /* create dataspace */
  dims_ful = Ntot            ;  fulSpace = H5Screate_simple(1, &dims_ful, NULL);
  dims_mem = Ntot / mpi->size;
  dims_loc =  num            ;  locSpace = H5Screate_simple(1, &dims_loc, NULL);
  /* configuration about domain decomposition */
  count  = 1;
  stride = 1;
  block  = dims_loc;
  offset = mpi->head;
  /* create chunked dataset */
  data_create = H5Pcreate(H5P_DATASET_CREATE);
  dims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
  if( dims_mem > dims_max )
    dims_mem = dims_max;
  chkHDF5err(H5Pset_chunk(data_create, 1, &dims_mem));

#ifndef BLOCK_TIME_STEP
  /* write particle velocity */
  dset = H5Dcreate(group, "vx", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, type.real, locSpace, hyperslab, w_property, body.vx));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));
  dset = H5Dcreate(group, "vy", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, type.real, locSpace, hyperslab, w_property, body.vy));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));
  dset = H5Dcreate(group, "vz", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, type.real, locSpace, hyperslab, w_property, body.vz));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));
#endif//BLOCK_TIME_STEP
  dset = H5Dcreate(group, "index", H5T_NATIVE_ULONG, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, H5T_NATIVE_ULONG, locSpace, hyperslab, w_property, body.idx));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));

  /* close/release resources */
  chkHDF5err(H5Pclose(data_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));


  /* write attribute data */
  /* create the data space for the attribute */
  hsize_t attr_dims = 1;
  hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  /* write current time */
  attribute = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time));
  chkHDF5err(H5Aclose(attribute));
  /* write time step */
  attribute = H5Acreate(group, "dt", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &dt));
  chkHDF5err(H5Aclose(attribute));
  /* write # of steps */
  attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  /* write # of steps */
  attribute = H5Acreate(group, "elapsed", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &elapsed));
  chkHDF5err(H5Aclose(attribute));
  /* write # of N-body particles */
  attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &Ntot));
  chkHDF5err(H5Aclose(attribute));
  /* write inverse of total energy at the initial condition */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Acreate(group, "E0inv", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &relEneErr.E0inv));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  /* write maximum value of the relative error of the total energy */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Acreate(group, "errMax", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &relEneErr.errMax));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  /* write flag about USE_DOUBLE_PRECISION */
#ifdef  USE_DOUBLE_PRECISION
  const int useDP = 1;
#else///USE_DOUBLE_PRECISION
  const int useDP = 0;
#endif//USE_DOUBLE_PRECISION
  attribute = H5Acreate(group, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* write flag about BLOCK_TIME_STEP */
#ifdef  BLOCK_TIME_STEP
  const int blockTimeStep = 1;
#else///BLOCK_TIME_STEP
  const int blockTimeStep = 0;
#endif//BLOCK_TIME_STEP
  attribute = H5Acreate(group, "blockTimeStep", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &blockTimeStep));
  chkHDF5err(H5Aclose(attribute));

  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(group));


#ifndef  RUN_WITHOUT_GOTHIC
  /* write parameters for auto-tuning (when steps > 0) */
  if( steps != 0 ){
    group = H5Gcreate(target, "parameters in auto-tuning", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* create dataspace */
    dims_ful = mpi->size;    fulSpace = H5Screate_simple(1, &dims_ful, NULL);
    dims_loc =         1;    locSpace = H5Screate_simple(1, &dims_loc, NULL);
    dims_mem =         1;
    /* configuration about domain decomposition */
    count  = 1;
    stride = 1;
    block  = dims_loc;
    offset = mpi->rank;
    /* create chunked dataset */
    data_create = H5Pcreate(H5P_DATASET_CREATE);
    dims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
    if( dims_mem > dims_max )
      dims_mem = dims_max;
    chkHDF5err(H5Pset_chunk(data_create, 1, &dims_mem));

    /* output rebuildTree */
    dset = H5Dcreate(group, "rebuild tree", type.rebuildTree, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, type.rebuildTree, locSpace, hyperslab, w_property, &rebuild));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));
    /* output measuredTime */
    dset = H5Dcreate(group, "measured time", type.measuredTime, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, type.measuredTime, locSpace, hyperslab, w_property, &measured));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));

    /* output statVal for linear growth model */
    dset = H5Dcreate(group, "stats (linear)", type.statVal, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, type.statVal, locSpace, hyperslab, w_property, &(rebuildParam.linearStats)));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));
    /* output guessTime for linear growth model */
    dset = H5Dcreate(group, "guess (linear)", type.guessTime, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, type.guessTime, locSpace, hyperslab, w_property, &(rebuildParam.linearGuess)));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));
    /* output statVal for power-law growth model */
    dset = H5Dcreate(group, "stats (power)", type.statVal, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, type.statVal, locSpace, hyperslab, w_property, &(rebuildParam.powerStats)));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));
    /* output guessTime for power-law growth model */
    dset = H5Dcreate(group, "guess (power)", type.guessTime, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, type.guessTime, locSpace, hyperslab, w_property, &(rebuildParam.powerGuess)));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
    /* output statVal for parabolic growth model */
    dset = H5Dcreate(group, "stats (parabolic)", type.statVal, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, type.statVal, locSpace, hyperslab, w_property, &(rebuildParam.parabolicStats)));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));
    /* output guessTime for parabolic growth model */
    dset = H5Dcreate(group, "guess (parabolic)", type.guessTime, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, type.guessTime, locSpace, hyperslab, w_property, &(rebuildParam.parabolicGuess)));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));
#endif//USE_PARABOLIC_GROWTH_MODEL

    /* output brentStatus */
    dset = H5Dcreate(group, "Brent status", type.brentStatus, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, type.brentStatus, locSpace, hyperslab, w_property, &status));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));
    /* output brentMemory */
    dset = H5Dcreate(group, "Brent memory", type.brentMemory, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, type.brentMemory, locSpace, hyperslab, w_property, &memory));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));

    /* output # of particles in each process */
    dset = H5Dcreate(group, "num", H5T_NATIVE_INT, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, H5T_NATIVE_INT, locSpace, hyperslab, w_property, &num));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));

    /* close/release resources */
    chkHDF5err(H5Pclose(data_create));
    chkHDF5err(H5Sclose(locSpace));
    chkHDF5err(H5Sclose(fulSpace));


    /* write attributes */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    int flag;
    /* write number of MPI processes in the current run */
    attribute = H5Acreate(group, "MPI_procs", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &mpi->size));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about FORCE_ADJUSTING_PARTICLE_TIME_STEPS */
#ifdef  FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    flag = 1;
#else///FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    flag = 0;
#endif//FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    attribute = H5Acreate(group, "FORCE_ADJUSTING_PARTICLE_TIME_STEPS", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about MONITOR_LETGEN_TIME */
#ifdef  MONITOR_LETGEN_TIME
    flag = 1;
#else///MONITOR_LETGEN_TIME
    flag = 0;
#endif//MONITOR_LETGEN_TIME
    attribute = H5Acreate(group, "MONITOR_LETGEN_TIME", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about USE_PARABOLIC_GROWTH_MODEL */
#ifdef  USE_PARABOLIC_GROWTH_MODEL
    flag = 1;
#else///USE_PARABOLIC_GROWTH_MODEL
    flag = 0;
#endif//USE_PARABOLIC_GROWTH_MODEL
    attribute = H5Acreate(group, "USE_PARABOLIC_GROWTH_MODEL", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));

    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
  }/* if( steps != 0 ){ */
#endif//RUN_WITHOUT_GOTHIC
  chkHDF5err(H5Pclose(w_property));

  /* write GPU information */
#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
  if( dumpGPUclock )
    appendGPUclockInfoParallel(monitor_step, deviceMonitors, *mpi, target, type);
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)

  /* close the file */
  /* finish collective dataset write */
  chkHDF5err(H5Fclose(target));
  *last ^= 1;
#endif//USE_HDF5_FORMAT


  __NOTE__("%s\n", "end");
}
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)


/**
 * @fn writeTipsyFile
 *
 * @brief Write particle data in TIPSY format.
 *
 * @param (time) current time of the simulation
 * @param (eps) softening length
 * @param (num) number of N-body particles
 * @param (body) the particle data
 * @param (file) name of the simulation
 */
void writeTipsyFile
(double time, float eps, int num, iparticle body, char file[]
#ifdef  ENABLE_GASEOUS_COMPONENT
 , const int Nsph, sph_particle sph
#endif//ENABLE_GASEOUS_COMPONENT
)
{
  __NOTE__("%s\n", "start");

#ifndef ENABLE_GASEOUS_COMPONENT
const int Nsph = 0;
#endif//ENABLE_GASEOUS_COMPONENT

  char filename[256];
  sprintf(filename, "%s/%s.tipsy", DATAFOLDER, file);
  FILE *fp;
  fp = fopen(filename, "wb");
  if( fp == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);      }

  bool success = true;
  size_t count;

  /** header */
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
  bonsaiHeader.nbodies = num + Nsph;
  bonsaiHeader.ndim = 3;
  bonsaiHeader.nsph = Nsph;
  bonsaiHeader.ndark = num;
  bonsaiHeader.nstar = 0;
  count = 1;  if( count != fwrite(&bonsaiHeader, sizeof(header), count, fp) )    success = false;

  /** main body */
#ifdef  ENABLE_GASEOUS_COMPONENT
  struct gas_particle {
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float hsmooth;
    float metals;
    /* float phi; */
    int idx;
  };
  struct gas_particle *gas;
  gas = (struct gas_particle *)malloc(sizeof(struct gas_particle) * Nsph);
  if( gas == NULL ){    __KILL__(stderr, "ERROR: failure to allocate gas\n");  }
  for(int ii = 0; ii < Nsph; ii++){
    gas[ii].mass = CAST_R2F(sph.pos[ii].m);
    gas[ii].pos[0] = CAST_R2F(sph.pos[ii].x);
    gas[ii].pos[1] = CAST_R2F(sph.pos[ii].y);
    gas[ii].pos[2] = CAST_R2F(sph.pos[ii].z);
    gas[ii].vel[0] = CAST_R2F(sph.vx[ii]);
    gas[ii].vel[1] = CAST_R2F(sph.vy[ii]);
    gas[ii].vel[2] = CAST_R2F(sph.vz[ii]);
    gas[ii].rho = CAST_R2F(sph.rho[ii]);
    gas[ii].temp = CAST_R2F(sph.T[ii]);
    gas[ii].hsmooth = 0.0f;
    gas[ii].metals = 0.0f;
    gas[ii].idx = (int)sph.idx[ii];
#if 0
    if( (ii & 1023) == 0 )
      fprintf(stdout, "m = %e, x = %e, y = %e, z = %e, vx = %e, vy = %e, vz = %e, rho = %e, T = %e, idx = %d\n", gas[ii].mass, gas[ii].pos[0], gas[ii].pos[1], gas[ii].pos[2], gas[ii].vel[0], gas[ii].vel[1], gas[ii].vel[2], gas[ii].rho, gas[ii].temp, gas[ii].idx);
#endif
  }
  count = Nsph;
  if( count != fwrite(gas, sizeof(struct gas_particle), count, fp) )
    success = false;
#endif//ENABLE_GASEOUS_COMPONENT

  struct dark_particle {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    int idx;
  };
  struct dark_particle *dark;
  dark = (struct dark_particle *)malloc(sizeof(struct dark_particle) * num);
  if( dark == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dark\n");  }
  for(int ii = 0; ii < num; ii++){
    dark[ii].mass   = CAST_R2F(body.pos[ii].m);
    dark[ii].pos[0] = CAST_R2F(body.pos[ii].x);
    dark[ii].pos[1] = CAST_R2F(body.pos[ii].y);
    dark[ii].pos[2] = CAST_R2F(body.pos[ii].z);
#ifdef  BLOCK_TIME_STEP
    dark[ii].vel[0] = CAST_R2F(body.vel[ii].x);
    dark[ii].vel[1] = CAST_R2F(body.vel[ii].y);
    dark[ii].vel[2] = CAST_R2F(body.vel[ii].z);
#else///BLOCK_TIME_STEP
    dark[ii].vel[0] = CAST_R2F(body.vx[ii]);
    dark[ii].vel[1] = CAST_R2F(body.vy[ii]);
    dark[ii].vel[2] = CAST_R2F(body.vz[ii]);
#endif//BLOCK_TIME_STEP
    dark[ii].eps    = eps;
    dark[ii].idx    = (int)body.idx[ii];
  }
  count = num;
  if( count != fwrite(dark, sizeof(struct dark_particle), count, fp) )
    success = false;

  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }

  fclose(fp);

  __NOTE__("%s\n", "end");
}


/**
 * @fn writeGalactICSFile
 *
 * @brief Write particle data in GalactICS format.
 *
 * @param (time) current time of the simulation
 * @param (head) head index of N-body particles
 * @param (num) number of N-body particles
 * @param (body) the particle data
 * @param (file) name of the file
 * @param (series) index of the component
 */
void writeGalactICSFile(double time, int head, int num, iparticle body, char file[], int series)
{
  __NOTE__("%s\n", "start");

  char filename[256];
  sprintf(filename, "%s/%s.%d", DATAFOLDER, file, series);
  FILE *fp;
  fp = fopen(filename, "w");
  if( fp == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);      }

  fprintf(fp, "%d  % .7e\n", num, time);

  for(int ii = head; ii < head + num; ii++)
    fprintf(fp, "% .7e % .7e % .7e % .7e % .7e % .7e % .7e\n",
	    body.pos[ii].m, body.pos[ii].x, body.pos[ii].y, body.pos[ii].z,
#ifdef  BLOCK_TIME_STEP
	    body.vel[ii].x, body.vel[ii].y, body.vel[ii].z
#else///BLOCK_TIME_STEP
	    body.vx[ii], body.vy[ii], body.vz[ii]
#endif//BLOCK_TIME_STEP
	    );

  fclose(fp);

  __NOTE__("%s\n", "end");
}


/**
 * @fn writeGADGETFile
 *
 * @brief Write particle data in GADGET format.
 */
void writeGADGETFile
(const int Ntot, double time, int kind, int skind, int * head, int * num, iparticle body, char file[]
#ifdef  USE_HDF5_FORMAT
 , hdf5struct type
#endif//USE_HDF5_FORMAT
)
{
  __NOTE__("%s\n", "start");

  /** header for GADGET */
  struct io_header{
    int npart[6];
    double mass[6];
    double time;
    double redshift;
    int flag_sfr;
    int flag_feedback;
    uint npartTotal[6];
    int flag_cooling;
    int num_files;
    double BoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    int flag_stellarage;
    int flag_metals;
    uint npartTotalHighWord[6];
    int flag_entropy_instead_u;
    char fill[60];
  };
  static struct io_header header;
  header.npart[0] = 0;/**< Gas */
  header.npart[1] = 0;/**< Halo */
  header.npart[2] = 0;/**< Disk */
  header.npart[3] = 0;/**< Bulge */
  header.npart[4] = 0;/**< Stars */
  header.npart[5] = 0;/**< Bndry */
  header.time = time;
  header.num_files = 1;
  header.flag_entropy_instead_u = 0;

  real *pos;  pos = (real *)malloc(Ntot * 3 * sizeof(real));  if( pos == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pos");  }
  real *vel;  vel = (real *)malloc(Ntot * 3 * sizeof(real));  if( vel == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vel");  }
  uint *idx;  idx = (uint *)malloc(Ntot     * sizeof(uint));  if( idx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate idx");  }
  real *mss;  mss = (real *)malloc(Ntot     * sizeof(real));  if( mss == NULL ){    __KILL__(stderr, "ERROR: failure to allocate mss");  }

  /* set 0-th component as halo particles */
  for(int ii = 0; ii < num[0]; ii++){
    pos[INDEX2D(Ntot, 3, ii, 0)] = body.pos[ii].x;
    pos[INDEX2D(Ntot, 3, ii, 1)] = body.pos[ii].y;
    pos[INDEX2D(Ntot, 3, ii, 2)] = body.pos[ii].z;
    mss[                 ii    ] = body.pos[ii].m;
#ifdef  BLOCK_TIME_STEP
    vel[INDEX2D(Ntot, 3, ii, 0)] = body.vel[ii].x;
    vel[INDEX2D(Ntot, 3, ii, 1)] = body.vel[ii].y;
    vel[INDEX2D(Ntot, 3, ii, 2)] = body.vel[ii].z;
#else///BLOCK_TIME_STEP
    vel[INDEX2D(Ntot, 3, ii, 0)] = body.vx[ii];
    vel[INDEX2D(Ntot, 3, ii, 1)] = body.vy[ii];
    vel[INDEX2D(Ntot, 3, ii, 2)] = body.vz[ii];
#endif//BLOCK_TIME_STEP
  }/* for(int ii = 0; ii < num[0]; ii++){ */
  header.npart[1] += num[0];
  int tail = num[0];

  /* set disk particles */
  for(int kk = skind; kk < kind; kk++){
    for(int ii = 0; ii < num[kk]; ii++){
      const int src = head[kk] + ii;
      const int dst = tail + ii;

      pos[INDEX2D(Ntot, 3, dst, 0)] = body.pos[src].x;
      pos[INDEX2D(Ntot, 3, dst, 1)] = body.pos[src].y;
      pos[INDEX2D(Ntot, 3, dst, 2)] = body.pos[src].z;
      mss[                 dst    ] = body.pos[src].m;
#ifdef  BLOCK_TIME_STEP
      vel[INDEX2D(Ntot, 3, dst, 0)] = body.vel[src].x;
      vel[INDEX2D(Ntot, 3, dst, 1)] = body.vel[src].y;
      vel[INDEX2D(Ntot, 3, dst, 2)] = body.vel[src].z;
#else///BLOCK_TIME_STEP
      vel[INDEX2D(Ntot, 3, dst, 0)] = body.vx[src];
      vel[INDEX2D(Ntot, 3, dst, 1)] = body.vy[src];
      vel[INDEX2D(Ntot, 3, dst, 2)] = body.vz[src];
#endif//BLOCK_TIME_STEP
    }/* for(int ii = 0; ii < num[kk]; ii++){ */
    header.npart[2] += num[kk];
    tail += num[kk];
  }/* for(int kk = skind; kk < kind; kk++){ */

  /* set remained spherical components as bulge particles */
  for(int kk = 1; kk < skind; kk++){
    for(int ii = 0; ii < num[kk]; ii++){
      const int src = head[kk] + ii;
      const int dst = tail + ii;

      pos[INDEX2D(Ntot, 3, dst, 0)] = body.pos[src].x;
      pos[INDEX2D(Ntot, 3, dst, 1)] = body.pos[src].y;
      pos[INDEX2D(Ntot, 3, dst, 2)] = body.pos[src].z;
      mss[                 dst    ] = body.pos[src].m;
#ifdef  BLOCK_TIME_STEP
      vel[INDEX2D(Ntot, 3, dst, 0)] = body.vel[src].x;
      vel[INDEX2D(Ntot, 3, dst, 1)] = body.vel[src].y;
      vel[INDEX2D(Ntot, 3, dst, 2)] = body.vel[src].z;
#else///BLOCK_TIME_STEP
      vel[INDEX2D(Ntot, 3, dst, 0)] = body.vx[src];
      vel[INDEX2D(Ntot, 3, dst, 1)] = body.vy[src];
      vel[INDEX2D(Ntot, 3, dst, 2)] = body.vz[src];
#endif//BLOCK_TIME_STEP
    }/* for(int ii = 0; ii < num[kk]; ii++){ */
    header.npart[3] += num[kk];
    tail += num[kk];
  }/* for(int kk = 1; kk < skind; kk++) */

  for(uint ii = 0; ii < (uint)Ntot; ii++)
    idx[ii] = ii;

  for(int kk = 0; kk < 6; kk++){
    header.mass[kk] = 0.0;
    header.npartTotal[kk] = (uint)header.npart[kk];
    header.npartTotalHighWord[kk] = 0;/**< most significant word of 64-bit total particle numbers (for N >= 2^23) */
  }/* for(int kk = 0; kk < 6; kk++){ */


  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
#ifndef USE_HDF5_FORMAT
  FILE *fp;
  sprintf(filename, "%s/%s.gadget", DATAFOLDER, file);
  fp = fopen(filename, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  int dummy = 256;
  bool success = true;
  size_t tmp;
  tmp =        1;  if( tmp != fwrite(& dummy, sizeof( dummy), tmp, fp) )    success = false;
  tmp =        1;  if( tmp != fwrite(&header, sizeof(header), tmp, fp) )    success = false;
  tmp =        1;  if( tmp != fwrite(& dummy, sizeof( dummy), tmp, fp) )    success = false;
  tmp = Ntot * 3;  if( tmp != fwrite(pos, sizeof(real), tmp, fp) )    success = false;
  tmp = Ntot * 3;  if( tmp != fwrite(vel, sizeof(real), tmp, fp) )    success = false;
  tmp = Ntot    ;  if( tmp != fwrite(idx, sizeof(uint), tmp, fp) )    success = false;
  tmp = Ntot    ;  if( tmp != fwrite(mss, sizeof(real), tmp, fp) )    success = false;

  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }
  fclose(fp);

#else///USE_HDF5_FORMAT

  sprintf(filename, "%s/%s.hdf5", DATAFOLDER, file);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t group = H5Gcreate(target, "Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* write attributes */
  hsize_t attr_dims = 6;
  hid_t dataspace;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  /* write number of particles */
  attribute = H5Acreate(group, "NumPart_ThisFile", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &header.npart));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "NumPart_Total", H5T_NATIVE_UINT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_UINT, &header.npartTotal));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "NumPart_Total_HighWord", H5T_NATIVE_UINT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_UINT, &header.npartTotalHighWord));
  chkHDF5err(H5Aclose(attribute));
  /* write particle mass */
  attribute = H5Acreate(group, "MassTable", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &header.mass));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));

  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  /* write current time */
  attribute = H5Acreate(group, "Time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &header.time));
  chkHDF5err(H5Aclose(attribute));
  /* write configuration */
  attribute = H5Acreate(group, "NumFilesPerSnapshot", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &header.num_files));
  chkHDF5err(H5Aclose(attribute));
  /* write misc */
  attribute = H5Acreate(group, "Flag_Entropy_ICs", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &header.flag_entropy_instead_u));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(group));

  tail = 0;
  for(int kk = 0; kk < 6; kk++)
    if( header.npart[kk] > 0 ){
      char grp[16];      sprintf(grp, "PartType%d", kk);
      group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      hid_t dataset;
      hsize_t dims = header.npart[kk] * 3;
      dataspace = H5Screate_simple(1, &dims, NULL);
      /* write particle position */
      dataset = H5Dcreate(group, "Coordinates", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pos[tail * 3]));
      chkHDF5err(H5Dclose(dataset));
      /* write particle velocity */
      dataset = H5Dcreate(group, "Velocities", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vel[tail * 3]));
      chkHDF5err(H5Dclose(dataset));
      chkHDF5err(H5Sclose(dataspace));

      dims = header.npart[kk];
      dataspace = H5Screate_simple(1, &dims, NULL);
      /* write particle index */
      dataset = H5Dcreate(group, "ParticleIDs", H5T_NATIVE_UINT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &idx[tail]));
      chkHDF5err(H5Dclose(dataset));
      /* write particle mass */
      dataset = H5Dcreate(group, "Masses", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mss[tail]));
      chkHDF5err(H5Dclose(dataset));
      chkHDF5err(H5Sclose(dataspace));

      chkHDF5err(H5Gclose(group));
      tail += header.npart[kk];
    }/* if( header.npart[kk] > 0 ){ */

  chkHDF5err(H5Fclose(target));
#endif//USE_HDF5_FORMAT

  free(pos);
  free(vel);
  free(idx);
  free(mss);

  __NOTE__("%s\n", "end");
}


/**
 * @fn writeAnthemFile
 *
 * @brief Write particle data in Anthem format.
 */
void writeAnthemFile(const ulong Ntot, iparticle body, const int unit, char file[])
{
  __NOTE__("%s\n", "start");


  /** open the target file */
  char filename[128];
  sprintf(filename, "%s/%s_dmp0.anthem", DATAFOLDER, file);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t group = H5Gcreate(target, "info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /** write attributes */
  hsize_t attr_dims = 1;
  hid_t dataspace;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  /** write attributes about configuration */
#ifdef  USE_DOUBLE_PRECISION
  const hid_t HDF5_FPTYPE_REAL = H5T_NATIVE_DOUBLE;
  int FP_PREQ = 64;
  double *write_buf;
  write_buf = (double *)malloc(sizeof(double) * Ntot);
  if( write_buf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate write_buf\n");  }
#else///USE_DOUBLE_PRECISION
  const hid_t HDF5_FPTYPE_REAL = H5T_NATIVE_FLOAT;
  int FP_PREQ = 32;
  float *write_buf;
  write_buf = (float *)malloc(sizeof(float) * Ntot);
  if( write_buf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate write_buf\n");  }
#endif//USE_DOUBLE_PRECISION
  attribute = H5Acreate(group, "FP_M", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &FP_PREQ));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "FP_L", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  FP_PREQ = 32;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &FP_PREQ));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "FP_H", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  FP_PREQ = 64;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &FP_PREQ));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));


  /** write attributes about physical properties of the system */
  group = H5Gcreate(target, "system", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  attribute = H5Acreate(group, "num", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &Ntot));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));


  /** write unit system */
  group = H5Gcreate(target, "unit", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  attribute = H5Acreate(group, "ID", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  int anthem_unit;
  switch( unit ){
  case ANTHEM_UNIT_UNITY:
    anthem_unit = -1;
    break;
  case ANTHEM_UNIT_GALAXY:
    anthem_unit = 0;
    break;
  case ANTHEM_UNIT_DWARF:
    anthem_unit = 1;
    break;
  case ANTHEM_UNIT_GLOBULAR_CLUSTER:
    anthem_unit = 2;
    break;
  case ANTHEM_UNIT_GALACTICS:
    anthem_unit = 10;
    break;
  default:
    __KILL__(stderr, "ERROR: undefined unit (%d) is passed\n", unit);
    break;
  }
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &anthem_unit));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "newton", HDF5_FPTYPE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, HDF5_FPTYPE_REAL, &newton));
  chkHDF5err(H5Aclose(attribute));
  extern const double lightspeed;
  const real lightspeed_real = CAST_D2R(lightspeed);
  attribute = H5Acreate(group, "speed_of_light", HDF5_FPTYPE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, HDF5_FPTYPE_REAL, &lightspeed_real));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_mass", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &mass2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_length", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &length2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_velocity", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &velocity2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &energy2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_specific_energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &senergy2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_momentum", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const double momentum2astro;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &momentum2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_angular_momentum", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const double angular_momentum2astro;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &angular_momentum2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_specific_angular_momentum", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const double specific_angular_momentum2astro;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &specific_angular_momentum2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_density", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &density2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_surface_density", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &col_density2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_phase_space_density_6D", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const double phase_space_density_6D2astro;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &phase_space_density_6D2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_phase_space_density", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const double phase_space_density2astro;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &phase_space_density2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_angle_phase_space_density", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const double angle_phase_space_density2astro;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &angle_phase_space_density2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_EJ_phase_space_density", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const double EJ_phase_space_density2astro;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &EJ_phase_space_density2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_number_density", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const double ndensity2astro;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &ndensity2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "com2astr_column_density", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const double colndensity2astro;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &colndensity2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "max_length", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  const int max_length = CONSTANTS_H_PLOT_WORDS;
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &max_length));
  chkHDF5err(H5Aclose(attribute));



  /* attr_dims = 1; */
  /* dataspace = H5Screate_simple(1, &attr_dims, NULL); */
  /* attribute = H5Acreate(subgroup, "length_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT); */
  /* chkHDF5err(H5Awrite(attribute, type.str4unit, length_astro_unit_name)); */
  /* chkHDF5err(H5Aclose(attribute)); */
  /* type->str4unit = H5Tcopy(H5T_C_S1); */
  /* chkHDF5err(H5Tset_size(type->str4unit, CONSTANTS_H_CHAR_WORDS));/\* memo: sizeof(char) is unity *\/ */



  // write texts
  hid_t H5T_ANTHEM_STRING = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(H5T_ANTHEM_STRING, H5T_VARIABLE));/* memo: sizeof(char) is unity */
  attribute = H5Acreate(group, "tex_time", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char time_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  char *tmp_str = (char *)malloc(sizeof(char) * CONSTANTS_H_PLOT_WORDS);
  sprintf(tmp_str, time_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_mass", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char mass_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, mass_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_length", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char length_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, length_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_velocity", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char velocity_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, velocity_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_energy", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char energy_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, energy_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_specific_energy", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char senergy_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, senergy_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_momentum", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char momentum_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, momentum_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_angular_momentum", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char angular_momentum_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, angular_momentum_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_specific_angular_momentum", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char specific_angular_momentum_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, specific_angular_momentum_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, density_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_surface_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char col_density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, col_density_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_phase_space_density_6D", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char phase_space_density_6D_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, phase_space_density_6D_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_phase_space_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char phase_space_density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, phase_space_density_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_angle_phase_space_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char angle_phase_space_density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, angle_phase_space_density_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_EJ_phase_space_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char EJ_phase_space_density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, EJ_phase_space_density_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_number_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char ndensity_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, ndensity_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "tex_column_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char colndensity_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
  sprintf(tmp_str, colndensity_astro_unit_name4plot);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_time", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  sprintf(tmp_str, time_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_mass", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  sprintf(tmp_str, mass_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_length", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  sprintf(tmp_str, length_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_velocity", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  sprintf(tmp_str, velocity_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_energy", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  sprintf(tmp_str, energy_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_specific_energy", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  sprintf(tmp_str, senergy_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_momentum", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char momentum_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, momentum_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_angular_momentum", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char angular_momentum_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, angular_momentum_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_specific_angular_momentum", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char specific_angular_momentum_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, specific_angular_momentum_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, density_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_surface_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char col_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, col_density_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_phase_space_density_6D", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char phase_space_density_6D_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, phase_space_density_6D_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_phase_space_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char phase_space_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, phase_space_density_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_angle_phase_space_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char angle_phase_space_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, angle_phase_space_density_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_EJ_phase_space_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char EJ_phase_space_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, EJ_phase_space_density_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_number_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char ndensity_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, ndensity_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ascii_column_density", H5T_ANTHEM_STRING, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  extern const char colndensity_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
  sprintf(tmp_str, colndensity_astro_unit_name);
  chkHDF5err(H5Awrite(attribute, H5T_ANTHEM_STRING, &tmp_str));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(H5T_ANTHEM_STRING));
  free(tmp_str);

  chkHDF5err(H5Gclose(group));
  chkHDF5err(H5Sclose(dataspace));


  /** write particle data */
  group = H5Gcreate(target, "nbody", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t dims = Ntot;
  dataspace = H5Screate_simple(1, &dims, NULL);
  hid_t dataset;
  dataset = H5Dcreate(group, "id", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.idx));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "mass", HDF5_FPTYPE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for(ulong ii = 0; ii < Ntot; ii++)
    write_buf[ii] = body.pos[ii].m;
  chkHDF5err(H5Dwrite(dataset, HDF5_FPTYPE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, write_buf));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "posx", HDF5_FPTYPE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for(ulong ii = 0; ii < Ntot; ii++)
    write_buf[ii] = body.pos[ii].x;
  chkHDF5err(H5Dwrite(dataset, HDF5_FPTYPE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, write_buf));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "posy", HDF5_FPTYPE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for(ulong ii = 0; ii < Ntot; ii++)
    write_buf[ii] = body.pos[ii].y;
  chkHDF5err(H5Dwrite(dataset, HDF5_FPTYPE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, write_buf));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "posz", HDF5_FPTYPE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for(ulong ii = 0; ii < Ntot; ii++)
    write_buf[ii] = body.pos[ii].z;
  chkHDF5err(H5Dwrite(dataset, HDF5_FPTYPE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, write_buf));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "velx", HDF5_FPTYPE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for(ulong ii = 0; ii < Ntot; ii++){
#ifdef  BLOCK_TIME_STEP
    write_buf[ii] = body.vel[ii].x;
#else///BLOCK_TIME_STEP
    write_buf[ii] = body.vx[ii];
#endif//BLOCK_TIME_STEP
  }
  chkHDF5err(H5Dwrite(dataset, HDF5_FPTYPE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, write_buf));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "vely", HDF5_FPTYPE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for(ulong ii = 0; ii < Ntot; ii++){
#ifdef  BLOCK_TIME_STEP
    write_buf[ii] = body.vel[ii].y;
#else///BLOCK_TIME_STEP
    write_buf[ii] = body.vy[ii];
#endif//BLOCK_TIME_STEP
  }
  chkHDF5err(H5Dwrite(dataset, HDF5_FPTYPE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, write_buf));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "velz", HDF5_FPTYPE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for(ulong ii = 0; ii < Ntot; ii++){
#ifdef  BLOCK_TIME_STEP
    write_buf[ii] = body.vel[ii].z;
#else///BLOCK_TIME_STEP
    write_buf[ii] = body.vz[ii];
#endif//BLOCK_TIME_STEP
  }
  chkHDF5err(H5Dwrite(dataset, HDF5_FPTYPE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, write_buf));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(group));

  /** close the file */
  chkHDF5err(H5Fclose(target));


  FILE *fp;
  sprintf(filename, "%s/%s.txt", DATAFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fprintf(fp, "_dmp0\n");
  fprintf(fp, "%d\n", 0);
  fclose(fp);


  __NOTE__("%s\n", "end");
}


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
/**
 * @fn writeFixedPotentialTable
 *
 * @brief Write fixed potential field represented in cubic spline interpolation.
 *
 * @param (unit) unit system of the potential field
 * @param (pot_tbl_sphe) superposed potential field for cubic spline interpolation (only for spherical symmetric components)
 * @param (skind) number of spherical symmetric components in MAGI
 * @param (pot_tbl) potential field for cubic spline interpolation
 * @param (binary) write in binary format when binary is true; otherwise, write in ASCII format
 */
void writeFixedPotentialTable
(const int unit, potential_field pot_tbl_sphe, const int skind, potential_field *pot_tbl,
#ifdef  USE_HDF5_FORMAT
 hdf5struct type,
#else///USE_HDF5_FORMAT
 const bool binary,
#endif//USE_HDF5_FORMAT
 char file[])
{
  __NOTE__("%s\n", "start");


  char cfgfile[128];
  sprintf(cfgfile, "%s/%s.%s.cfg", DATAFOLDER, file, "ext_pot_sphe");
  FILE *fp_cfg;
  fp_cfg = fopen(cfgfile, "w");
  if( fp_cfg == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", cfgfile);  }
  fprintf(fp_cfg, "%d\n", 1);

  char filename[256];

#ifdef  USE_HDF5_FORMAT

  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "ext_pot");
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t dims = pot_tbl_sphe.num;
  hid_t dataspace = H5Screate_simple(1, &dims, NULL);

  /* preparation for data compression */
  hid_t dataset, property;
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  uint szip_pixels_per_block = 8;
  hsize_t cdims = 128 * szip_pixels_per_block;
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION

#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  hsize_t cdims_loc = cdims;
  if( dims < cdims_loc )
    cdims_loc = (hsize_t)1 << llog2((ulong)dims);
  hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  if( cdims_loc > cdims_max )
    cdims_loc = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION


  /* write potential table for superposed spherical components */
  fprintf(fp_cfg, "%s\n", "spherical");
  hid_t group = H5Gcreate(target, "spherical", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* write radius */
  dataset = H5Dcreate(group, "r", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, pot_tbl_sphe.rad));
  chkHDF5err(H5Dclose(dataset));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* write potential */
  dataset = H5Dcreate(group, "Phi(r)", type.pot2, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.pot2, H5S_ALL, H5S_ALL, H5P_DEFAULT, pot_tbl_sphe.Phi));
  chkHDF5err(H5Dclose(dataset));

  /* write attribute data */
  hsize_t attr_dims = 1;
  hid_t attrspace = H5Screate_simple(1, &attr_dims, NULL);
  /* write # of data points */
  hid_t attribute = H5Acreate(group, "num", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &(pot_tbl_sphe.num)));
  chkHDF5err(H5Aclose(attribute));
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* write log_10(r_min) */
  attribute = H5Acreate(group, "logrmin", type.real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &(pot_tbl_sphe.logrmin)));
  chkHDF5err(H5Aclose(attribute));
  /* write logrbin */
  attribute = H5Acreate(group, "logrbin", type.real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &(pot_tbl_sphe.logrbin)));
  chkHDF5err(H5Aclose(attribute));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  chkHDF5err(H5Gclose(group));


  /* write potential table for each component */
  if( skind > 1 ){
    group = H5Gcreate(target, "individual", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int ii = 0; ii < skind; ii++){
      if( (hsize_t)pot_tbl[ii].num != dims ){
#ifdef  USE_SZIP_COMPRESSION
	chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
	chkHDF5err(H5Sclose(dataspace));
	dims = pot_tbl[ii].num;
	dataspace = H5Screate_simple(1, &dims, NULL);

#ifdef  USE_SZIP_COMPRESSION
	property = H5Pcreate(H5P_DATASET_CREATE);
	cdims_loc = cdims;
	if( dims < cdims_loc )
	  cdims_loc = (hsize_t)1 << llog2((ulong)dims);
	if( cdims_loc > cdims_max )
	  cdims_loc = cdims_max;
	chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
	chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
      }/* if( (hsize_t)pot_tbl[ii].num != dims ){ */

      char grp[16];
      sprintf(grp, "data%d", ii);
      fprintf(fp_cfg, "%s/%s\n", "individual", grp);
      hid_t sub = H5Gcreate(group, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
      /* write radius */
      dataset = H5Dcreate(sub, "r", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, pot_tbl[ii].rad));
      chkHDF5err(H5Dclose(dataset));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
      /* write potential */
      dataset = H5Dcreate(sub, "Phi(r)", type.pot2, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, type.pot2, H5S_ALL, H5S_ALL, H5P_DEFAULT, pot_tbl[ii].Phi));
      chkHDF5err(H5Dclose(dataset));

      /* write attribute data */
      /* write # of data points */
      attribute = H5Acreate(sub, "num", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &(pot_tbl[ii].num)));
      chkHDF5err(H5Aclose(attribute));
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
      /* write log_10(r_min) */
      attribute = H5Acreate(sub, "logrmin", type.real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, type.real, &(pot_tbl[ii].logrmin)));
      chkHDF5err(H5Aclose(attribute));
      /* write logrbin */
      attribute = H5Acreate(sub, "logrbin", type.real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, type.real, &(pot_tbl[ii].logrbin)));
      chkHDF5err(H5Aclose(attribute));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

      chkHDF5err(H5Gclose(sub));
    }/* for(int ii = 0; ii < skind; ii++){ */
    chkHDF5err(H5Gclose(group));
  }/* if( skind > 1 ){ */

#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION

  /* write attribute data */
  /* write # of spherical components */
  attribute = H5Acreate(target, "skind", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skind));
  chkHDF5err(H5Aclose(attribute));
  /* write flag about USE_DOUBLE_PRECISION */
#ifdef  USE_DOUBLE_PRECISION
  const int useDP = 1;
#else///USE_DOUBLE_PRECISION
  const int useDP = 0;
#endif//USE_DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));

  /* write unit system */
  group = H5Gcreate(target, "unit_system", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* write index of unit system */
  attribute = H5Acreate(group, "unit", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &unit));
  chkHDF5err(H5Aclose(attribute));
  /* write gravitational constant G */
  attribute = H5Acreate(group, "newton", type.real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &newton));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));

  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Sclose(attrspace));

  /* close the file */
  chkHDF5err(H5Fclose(target));

#else///USE_HDF5_FORMAT

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  int Ndat = pot_tbl_sphe.num;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  const int Ndat = pot_tbl_sphe.num;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  /* write numeric table for superposed spherical components */
  sprintf(filename, "%s/%s.%s.dat", DATAFOLDER, file, "ext_pot");
  FILE *fp;
  fp = fopen(filename, binary ? "wb" : "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
  if( binary ){
    bool success = true;
    size_t tmp;
    tmp = 1;    if( tmp != fwrite(&unit, sizeof(int), tmp, fp) )      success = false;
    tmp = 1;    if( tmp != fwrite(&Ndat, sizeof(int), tmp, fp) )      success = false;
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = 1;    if( tmp != fwrite(&(pot_tbl_sphe.logrmin), sizeof(real), tmp, fp) )      success = false;
    tmp = 1;    if( tmp != fwrite(&(pot_tbl_sphe.logrbin), sizeof(real), tmp, fp) )      success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

    fprintf(fp_cfg, "%d\n", READ_SUPERPOSED_TABLE_SPHE);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = Ndat;    if( tmp != fwrite(pot_tbl_sphe.rad, sizeof(real), tmp, fp) )      success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = Ndat;    if( tmp != fwrite(pot_tbl_sphe.Phi, sizeof(pot2), tmp, fp) )      success = false;

    if( success != true ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
  }/* if( binary ){ */
  else{
    fprintf(fp, "%d\n", unit);
    fprintf(fp, "%d\n", Ndat);
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    fprintf(fp, "%e\t%e\n", pot_tbl_sphe.logrmin, pot_tbl_sphe.logrbin);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

    fprintf(fp_cfg, "%d\n", READ_SUPERPOSED_TABLE_SPHE);
    for(int ii = 0; ii < Ndat; ii++){
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
      fprintf(fp, "%e\t%e\t%e\n", pot_tbl_sphe.rad[ii], pot_tbl_sphe.Phi[ii].val, pot_tbl_sphe.Phi[ii].dr2);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
      fprintf(fp, "%e\t%e\n", pot_tbl_sphe.Phi[ii].Phi, pot_tbl_sphe.Phi[ii].Fr);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    }/* for(int ii = 0; ii < Ndat; ii++){ */
  }/* else{ */
  fclose(fp);


  /* write numeric table for each component */
  if( skind >= 1 )
    for(int kk = 0; kk < skind; kk++){
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
      Ndat = pot_tbl[kk].num;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

      fprintf(fp_cfg, "%d\n", kk);
      sprintf(filename, "%s/%s.%s.%d", DATAFOLDER, file, "pot", kk);
      fp = fopen(filename, binary ? "wb" : "w");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);      }
      if( binary ){
	bool success = true;
	size_t tmp;
	tmp =    1;	if( tmp != fwrite(&unit, sizeof(int), tmp, fp) )	  success = false;
	tmp =    1;	if( tmp != fwrite(&Ndat, sizeof(int), tmp, fp) )	  success = false;
	tmp = Ndat;	if( tmp != fwrite(pot_tbl[kk].Phi, sizeof(pot2), tmp, fp) )	  success = false;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	tmp = Ndat;	if( tmp != fwrite(pot_tbl[kk].rad, sizeof(real), tmp, fp) )	  success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

	if( success != true ){	  __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);	}
      }/* if( binary ){ */
      else{
	fprintf(fp, "%d\n", unit);
	fprintf(fp, "%d\n", Ndat);
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
        fprintf(fp, "%e\t%e\n", pot_tbl[kk].logrmin, pot_tbl[kk].logrbin);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

	for(int ii = 0; ii < Ndat; ii++)
	  fprintf(fp, "%e\t%e\t%e\n",
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
		  pot_tbl[kk].rad[ii],
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	          POW10(pot_tbl[kk].logrmin + pot_tbl[kk].logrbin * (real)ii),
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
		  pot_tbl[kk].Phi[ii].Phi, pot_tbl[kk].Phi[ii].Fr);
      }/* else{ */
      fclose(fp);
    }/* for(int kk = 0; kk < skind; kk++){ */

#endif//USE_HDF5_FORMAT

  fclose(fp_cfg);


  __NOTE__("%s\n", "end");
}


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
/**
 * @fn writeFixedDiskPotential
 *
 * @brief Write fixed potential field by disk components
 *
 * @param (unit) unit system of the potential field
 * @param (disk) potential field by the disk components
 * @param (binary) write in binary format when binary is true; otherwise, write in ASCII format
 */
void writeFixedDiskPotential
(const int unit, const disk_potential disk
#ifdef  USE_HDF5_FORMAT
 , hdf5struct type
#else///USE_HDF5_FORMAT
 , const bool binary
#endif//USE_HDF5_FORMAT
 , char *file)
{
  __NOTE__("%s\n", "start");
#if 0
  for(int ii = 0; ii < disk.NR + 1; ii++){
    for(int jj = 0; jj < disk.Nz + 1; jj++)
      fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", disk.hh * (real)ii, disk.hh * (real)jj, disk.Phi[INDEX2D(disk.NR + 1, disk.Nz + 1, ii, jj)], disk.FRz[INDEX2D(disk.NR + 1, disk.Nz + 1, ii, jj)].R, disk.FRz[INDEX2D(disk.NR + 1, disk.Nz + 1, ii, jj)].z);
    fprintf(stderr, "\n");
  }/* for(int ii = 0; ii < disk.NR + 1; ii++){ */
#endif

  char filename[128];
#ifdef  USE_HDF5_FORMAT
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "ext_disk");
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  hsize_t dims = (disk.NR + 1) * (disk.Nz + 1) * disk.maxLev;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  hsize_t dims = (disk.NR + 1) * (disk.Nz + 1);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  hid_t dataspace = H5Screate_simple(1, &dims, NULL);

  /* preparation for data compression */
  hid_t dataset, property;
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  uint szip_pixels_per_block = 8;
  hsize_t cdims = 128 * szip_pixels_per_block;
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION

#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  hsize_t cdims_loc = cdims;
  if( dims < cdims_loc )
    cdims_loc = (hsize_t)1 << llog2((ulong)dims);
  hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  if( cdims_loc > cdims_max )
    cdims_loc = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION


  /* write potential table for superposed spherical components */
  hid_t group = H5Gcreate(target, "2D", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* write \Phi(R, z) */
  dataset = H5Dcreate(group, "Phi", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk.Phi));
  chkHDF5err(H5Dclose(dataset));
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* write F(R, z) */
  dataset = H5Dcreate(group, "FRz", type.disk_grav, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.disk_grav, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk.FRz));
  chkHDF5err(H5Dclose(dataset));
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* write R */
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));
  dims = disk.maxLev * disk.NR;
  dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  cdims_loc = cdims;
  if( dims < cdims_loc )
    cdims_loc = (hsize_t)1 << llog2((ulong)dims);
  if( cdims_loc > cdims_max )
    cdims_loc = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
  dataset = H5Dcreate(group, "R", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk.RR));
  chkHDF5err(H5Dclose(dataset));
  /* write z */
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));
  dims = disk.maxLev * disk.Nz;
  dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  cdims_loc = cdims;
  if( dims < cdims_loc )
    cdims_loc = (hsize_t)1 << llog2((ulong)dims);
  if( cdims_loc > cdims_max )
    cdims_loc = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
  dataset = H5Dcreate(group, "z", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk.zz));
  chkHDF5err(H5Dclose(dataset));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

  /* write attribute data */
  hsize_t attr_dims = 1;
  hid_t attrspace = H5Screate_simple(1, &attr_dims, NULL);
  /* write hh */
  hid_t attribute = H5Acreate(group, "hh", type.real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &(disk.hh)));
  chkHDF5err(H5Aclose(attribute));
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* write # of nested levels */
  attribute = H5Acreate(group, "maxLev", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &(disk.maxLev)));
  chkHDF5err(H5Aclose(attribute));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* write # of R grids */
  attribute = H5Acreate(group, "NR", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &(disk.NR)));
  chkHDF5err(H5Aclose(attribute));
  /* write # of z grids */
  attribute = H5Acreate(group, "Nz", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &(disk.Nz)));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));

  /* write spherical averaged potential profile */
  group = H5Gcreate(target, "spherical", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));
  dims = disk.sphe.num;
  dataspace = H5Screate_simple(1, &dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  cdims_loc = cdims;
  if( dims < cdims_loc )
    cdims_loc = (hsize_t)1 << llog2((ulong)dims);
  if( cdims_loc > cdims_max )
    cdims_loc = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* write radius */
  dataset = H5Dcreate(group, "r", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk.sphe.rad));
  chkHDF5err(H5Dclose(dataset));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* write potential */
  dataset = H5Dcreate(group, "Phi(r)", type.pot2, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.pot2, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk.sphe.Phi));
  chkHDF5err(H5Dclose(dataset));

  /* write attribute data */
  attr_dims = 1;
  attrspace = H5Screate_simple(1, &attr_dims, NULL);
  /* write # of data points */
  attribute = H5Acreate(group, "num", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &(disk.sphe.num)));
  chkHDF5err(H5Aclose(attribute));
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  /* write log_10(r_min) */
  attribute = H5Acreate(group, "logrmin", type.real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &(disk.sphe.logrmin)));
  chkHDF5err(H5Aclose(attribute));
  /* write logrbin */
  attribute = H5Acreate(group, "logrbin", type.real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &(disk.sphe.logrbin)));
  chkHDF5err(H5Aclose(attribute));
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  chkHDF5err(H5Gclose(group));

#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION

  /* write attribute data */
  /* write flag about USE_DOUBLE_PRECISION */
#ifdef  USE_DOUBLE_PRECISION
  const int useDP = 1;
#else///USE_DOUBLE_PRECISION
  const int useDP = 0;
#endif//USE_DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));

  /* write unit system */
  group = H5Gcreate(target, "unit_system", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* write index of unit system */
  attribute = H5Acreate(group, "unit", H5T_NATIVE_INT, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &unit));
  chkHDF5err(H5Aclose(attribute));
  /* write gravitational constant G */
  attribute = H5Acreate(group, "newton", type.real, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &newton));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));

  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Sclose(attrspace));

  /* close the file */
  chkHDF5err(H5Fclose(target));

#else///USE_HDF5_FORMAT

  sprintf(filename, "%s/%s.%s.dat", DATAFOLDER, file, "ext_disk");
  FILE *fp;
  fp = fopen(filename, binary ? "wb" : "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
  if( binary ){
    bool success = true;
    size_t tmp;
    tmp = 1;    if( tmp != fwrite(&unit, sizeof(int), tmp, fp) )      success = false;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = 1;    if( tmp != fwrite(&disk.maxLev, sizeof(int), tmp, fp) )      success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = 1;    if( tmp != fwrite(&disk.NR, sizeof(int), tmp, fp) )      success = false;
    tmp = 1;    if( tmp != fwrite(&disk.Nz, sizeof(int), tmp, fp) )      success = false;
    tmp = 1;    if( tmp != fwrite(&disk.hh, sizeof(real), tmp, fp) )      success = false;

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = disk.maxLev *  disk.NR                     ;    if( tmp != fwrite(disk.RR , sizeof(real), tmp, fp) )      success = false;
    tmp = disk.maxLev *                  disk.Nz     ;    if( tmp != fwrite(disk.zz , sizeof(real), tmp, fp) )      success = false;
    tmp = disk.maxLev * (disk.NR + 1) * (disk.Nz + 1);    if( tmp != fwrite(disk.Phi, sizeof(real), tmp, fp) )      success = false;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = (disk.NR + 1) * (disk.Nz + 1);    if( tmp != fwrite(disk.Phi, sizeof(real), tmp, fp) )      success = false;
    tmp = (disk.NR + 1) * (disk.Nz + 1);    if( tmp != fwrite(disk.FRz, sizeof(disk_grav), tmp, fp) )      success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

    tmp = 1;    if( tmp != fwrite(&(disk.sphe.num), sizeof(int), tmp, fp) )      success = false;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = disk.sphe.num;    if( tmp != fwrite(disk.sphe.rad, sizeof(real), tmp, fp) )      success = false;
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = 1;    if( tmp != fwrite(&(disk.sphe.logrmin), sizeof(real), tmp, fp) )      success = false;
    tmp = 1;    if( tmp != fwrite(&(disk.sphe.logrbin), sizeof(real), tmp, fp) )      success = false;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    tmp = disk.sphe.num;    if( tmp != fwrite(disk.sphe.Phi, sizeof(pot2), tmp, fp) )      success = false;

    if( success != true ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
  }/* if( binary ){ */
  else{
    fprintf(fp, "%d\n", unit);
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    fprintf(fp, "%d\t%d\t%d\n", disk.maxLev, disk.NR, disk.Nz);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    fprintf(fp, "%d\t%d\n", disk.NR, disk.Nz);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    fprintf(fp, "%e\n", disk.hh);

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    for(int ii = 0; ii < disk.maxLev * disk.NR; ii++)
      fprintf(fp, "%e\n", disk.RR[ii]);
    for(int ii = 0; ii < disk.maxLev * disk.Nz; ii++)
      fprintf(fp, "%e\n", disk.zz[ii]);
    for(int ii = 0; ii < disk.maxLev * (disk.NR + 1) * (disk.Nz + 1); ii++)
      fprintf(fp, "%e\n", disk.Phi[ii]);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    for(int ii = 0; ii < disk.NR + 1; ii++){
      const real RR = disk.hh * (real)ii;
      for(int jj = 0; jj < disk.Nz + 1; jj++){
	const int idx = INDEX2D(disk.NR + 1, disk.Nz + 1, ii, jj);
	fprintf(fp, "%e\t%e\t%e\t%e\t%e\n", RR, disk.hh * (real)jj, disk.Phi[idx], disk.FRz[idx].R, disk.FRz[idx].z);
      }/* for(int jj = 0; jj < disk.Nz + 1; jj++){ */
    }/* for(int ii = 0; ii < disk.NR + 1; ii++){ */
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

    fprintf(fp, "%d\n", disk.sphe.num);
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
    fprintf(fp, "%e\t%e\n", disk.sphe.logrmin, disk.sphe.logrbin);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

    for(int ii = 0; ii < disk.sphe.num; ii++)
      fprintf(fp, "%e\t%e\t%e\n",
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	      disk.sphe.rad[ii],
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	      LOG10(disk.sphe.logrmin + disk.sphe.logrbin * (real)ii),
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
	      disk.sphe.Phi[ii].Phi, disk.sphe.Phi[ii].Fr);
  }/* if( binary ){ */
  fclose(fp);

#endif//USE_HDF5_FORMAT

  __NOTE__("%s\n", "end");
}
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK

#endif//SET_EXTERNAL_POTENTIAL_FIELD


/**
 * @fn readSnapshot
 *
 * @brief Read snapshot file of the N-body simulation.
 *
 * @return (unit) unit system of the simulation
 * @return (time) current time of the simulation
 * @return (steps) number of time steps calculated in the simulation
 * @param (num) number of N-body particles
 * @param (file) name of the simulation
 * @param (id) ID of the snapshot
 * @return (body) the particle data
 * @param (type) data type for HDF5 (only for HDF5 enabled runs)
 */
void  readSnapshot(int *unit, double *time, ulong *steps, int num, char file[], uint id
#ifdef  USE_HDF5_FORMAT
		   , nbody_hdf5 *body, hdf5struct type
#else///USE_HDF5_FORMAT
		   , iparticle body
#endif//USE_HDF5_FORMAT
		   )
{
  __NOTE__("%s\n", "start");


  /* open an existing file with read only option */
  char filename[128];
#ifndef USE_HDF5_FORMAT
  FILE *fp;
  sprintf(filename, "%s/%s.%s%.3u.dat", DATAFOLDER, file, SNAPSHOT, id);
  fp = fopen(filename, "rb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fread( unit, sizeof(           int), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread( time, sizeof(        double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(steps, sizeof(         ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.pos, sizeof(    position), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.acc, sizeof(acceleration), tmp, fp) )    success = false;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  tmp = num;  if( tmp != fread(body.acc_ext, sizeof(acceleration), tmp, fp) )    success = false;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  tmp = num;  if( tmp != fread(body.vel, sizeof(velocity), tmp, fp) )    success = false;
#else///BLOCK_TIME_STEP
  tmp = num;  if( tmp != fread(body.vx, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.vy, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.vz, sizeof(real), tmp, fp) )    success = false;
#endif//BLOCK_TIME_STEP
  tmp = num;  if( tmp != fread(body.idx, sizeof(ulong), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }

  fclose(fp);
#else///USE_HDF5_FORMAT
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, SNAPSHOT, id);
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  /* open an existing group */
  char groupname[16];
  sprintf(groupname, SNAPSHOT);
  hid_t group = H5Gopen(target, groupname, H5P_DEFAULT);


  /* read particle data */
  hid_t dataset;
  /* read particle position */
  dataset = H5Dopen(group, "position", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->pos));
  chkHDF5err(H5Dclose(dataset));
  /* read particle velocity */
  dataset = H5Dopen(group, "velocity", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->vel));
  chkHDF5err(H5Dclose(dataset));
  /* read particle acceleration */
  dataset = H5Dopen(group, "acceleration", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->acc));
  chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* read particle acceleration by external potential field */
  dataset = H5Dopen(group, "acceleration_external", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->acc_ext));
  chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  /* read particle mass */
  dataset = H5Dopen(group, "mass", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->m));
  chkHDF5err(H5Dclose(dataset));
  /* read particle potential */
  dataset = H5Dopen(group, "potential", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->pot));
  chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* read particle potential by external potential field */
  dataset = H5Dopen(group, "potential_external", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->pot_ext));
  chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  /* read particle index */
  dataset = H5Dopen(group, "index", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->idx));
  chkHDF5err(H5Dclose(dataset));


  /* read attribute data */
  hid_t attribute;
  /* read current time */
  attribute = H5Aopen(group, "time", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, time));
  chkHDF5err(H5Aclose(attribute));
  /* read # of steps */
  attribute = H5Aopen(group, "steps", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, steps));
  chkHDF5err(H5Aclose(attribute));
  /* read # of N-body particles */
  ulong num_ulong = 0;
  attribute = H5Aopen(group, "number", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, &num_ulong));
  chkHDF5err(H5Aclose(attribute));
  /* read flag about USE_DOUBLE_PRECISION */
  int useDP;
  attribute = H5Aopen(group, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));

  chkHDF5err(H5Gclose(group));



  /* read unit system */
  group = H5Gopen(target, "unit_system", H5P_DEFAULT);
  attribute = H5Aopen(group, "unit", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, unit));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));



  /* close the file */
  chkHDF5err(H5Fclose(target));


  /* simple error check */
  if( num_ulong != (ulong)num ){
    __KILL__(stderr, "ERROR: number of N-body particles in the file (%zu) differs with that in the code (%d)\n", num_ulong, num);
  }
#ifdef  USE_DOUBLE_PRECISION
  if( useDP != 1 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);
  }/* if( useDP != 1 ){ */
#else///USE_DOUBLE_PRECISION
  if( useDP != 0 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);
  }/* if( useDP != 0 ){ */
#endif//USE_DOUBLE_PRECISION
#endif//USE_HDF5_FORMAT


  __NOTE__("%s\n", "end");
}


#ifdef  REPORT_COMPUTE_RATE
static inline void writeComputeRate(char file[], const uint id, const double time, const ulong steps, const double speed, const double speed_run, const double complete, double guess, const double brent_avg, const double rebuild_interval)
{
  char filename[128];
  sprintf(filename, "%s/%s.%s.log", LOGFOLDER, file, "speed");
  FILE *fp;
  if( id == 0 ){
    fp = fopen(filename, "w");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    fprintf(fp, "#time(%s)\tsteps\tfile\trate(s/step)\tspeed(%s/s)\tbrent_rate\tmake_interval\tcomplete\tremain\n", time_astro_unit_name, time_astro_unit_name);
  }/* if( id == 0 ){ */
  else{
    fp = fopen(filename, "a");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
  }/* else{ */

  if( steps == 0 )
    guess = 0.0;
  const int day = (int)floor(guess * 1.15740740740e-5);
  guess -= (double)day * 86400.0;/**< 1 day = 24h * 3600s = 84600s */
  const int hour = (int)floor(guess * 2.77777777778e-4);
  guess -= (double)hour * 3600.0;
  const int minute = (int)floor(guess * 1.66666666667e-2);
  guess -= (double)minute * 60.0;

  fprintf(fp, "%e\t%zu\t%u\t%e\t%e\t%e\t%e\t%7.3lf%%\t%dd %2dh %2dm %6.3lfs\n", time * time2astro, steps, id, speed, speed_run, brent_avg, rebuild_interval, complete, day, hour, minute, guess);
  fclose(fp);
}
#endif//REPORT_COMPUTE_RATE


#ifdef  PREPARE_XDMF_FILES
static inline void writeXdmf(const int num, char file[], const uint id)
{
  __NOTE__("%s\n", "start");

  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.xmf", DATAFOLDER, file, SNAPSHOT, id);
  FILE *fp;
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  sprintf(filename, "%s.%s%.3u.h5", file, SNAPSHOT, id);

  fprintf(fp, "<?xml version=\"1.0\" ?>\n");
  fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(fp, "<Xdmf Version=\"3.0\">\n");
  fprintf(fp, "<Domain>\n");
  fprintf(fp, "<Grid Name=\"N-body\" GridType=\"Uniform\">\n");
  fprintf(fp, "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n", num);
  fprintf(fp, "<Geometry GeometryType=\"XYZ\">\n");
  fprintf(fp, "<DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", num, sizeof(real));
  fprintf(fp, "%s:/%s/position\n", filename, SNAPSHOT);
  fprintf(fp, "</DataItem>\n");
  fprintf(fp, "</Geometry>\n");
  fprintf(fp, "<Attribute Name=\"velocity\" AttributeType=\"Vector\" Center=\"Node\">\n");
  fprintf(fp, "<DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", num, sizeof(real));
  fprintf(fp, "%s:/%s/velocity\n", filename, SNAPSHOT);
  fprintf(fp, "</DataItem>\n");
  fprintf(fp, "</Attribute>\n");
  fprintf(fp, "<Attribute Name=\"acceleration\" AttributeType=\"Vector\" Center=\"Node\">\n");
  fprintf(fp, "<DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", num, sizeof(real));
  fprintf(fp, "%s:/%s/acceleration\n", filename, SNAPSHOT);
  fprintf(fp, "</DataItem>\n");
  fprintf(fp, "</Attribute>\n");
  fprintf(fp, "<Attribute Name=\"mass\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", num, sizeof(real));
  fprintf(fp, "%s:/%s/mass\n", filename, SNAPSHOT);
  fprintf(fp, "</DataItem>\n");
  fprintf(fp, "</Attribute>\n");
  fprintf(fp, "<Attribute Name=\"potential\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", num, sizeof(real));
  fprintf(fp, "%s:/%s/potential\n", filename, SNAPSHOT);
  fprintf(fp, "</DataItem>\n");
  fprintf(fp, "</Attribute>\n");
  fprintf(fp, "<Attribute Name=\"index\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"UInt\" Precision=\"%zu\" Format=\"HDF\">\n", num, sizeof(ulong));
  fprintf(fp, "%s:/%s/index\n", filename, SNAPSHOT);
  fprintf(fp, "</DataItem>\n");
  fprintf(fp, "</Attribute>\n");
  fprintf(fp, "</Grid>\n");
  fprintf(fp, "</Domain>\n");
  fprintf(fp, "</Xdmf>\n");

  fclose(fp);

  __NOTE__("%s\n", "end");
}
#endif//PREPARE_XDMF_FILES


/**
 * @fn writeSnapshot
 *
 * @brief Write snapshot file of the N-body simulation.
 *
 * @param (unit) unit system of the simulation
 * @param (time) current time of the simulation
 * @param (steps) number of time steps calculated in the simulation
 * @param (num) number of N-body particles
 * @param (file) name of the simulation
 * @param (id) ID of the snapshot
 * @param (body) the particle data
 * @param (type) data type for HDF5 (only for HDF5 enabled runs)
 */
void writeSnapshot
(int  unit, double  time, ulong  steps, int num, char file[], uint id
#ifdef  USE_HDF5_FORMAT
 , nbody_hdf5 *body, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
 , energyError *relEneErr, const bool update
#endif//MONITOR_ENERGY_ERROR
#else///USE_HDF5_FORMAT
 , iparticle body
#endif//USE_HDF5_FORMAT
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
 , gpu_clock *deviceMonitors, const int monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  REPORT_COMPUTE_RATE
 , const double speed, const double speed_run, const double complete, const double guess, const double brent_avg, const double rebuild_interval
#endif//REPORT_COMPUTE_RATE
 )
{
  __NOTE__("%s\n", "start");


  char filename[128];

  /* calculate total energy */
#   if  defined(USE_HDF5_FORMAT) && defined(MONITOR_ENERGY_ERROR)
  double Ekin = 0.0;
  double Epot = 0.0;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  double Eext = 0.0;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  if( update ){
    for(int ii = 0; ii < num; ii++){
      const double mass = CAST_R2D(body->m  [ii	     ]);
      const double velx = CAST_R2D(body->vel[ii * 3    ]);
      const double vely = CAST_R2D(body->vel[ii * 3 + 1]);
      const double velz = CAST_R2D(body->vel[ii * 3 + 2]);
      Ekin += mass * (velx * velx + vely * vely + velz * velz);
      Epot += mass * CAST_R2D(body->pot[ii]);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      Eext += mass * CAST_R2D(body->pot_ext[ii]);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* for(int ii = 0; ii < num; ii++){ */

    Ekin *= 0.5 * energy2astro;
#ifndef SET_EXTERNAL_POTENTIAL_FIELD
    Epot *= 0.5 * energy2astro;
#else///SET_EXTERNAL_POTENTIAL_FIELD
    Epot = Eext + 0.5 * Epot;
    Epot *= energy2astro;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  }
  const double Etot = Ekin + Epot;

  FILE *fp;
  if( update ){
    sprintf(filename, "%s/%s.%s.log", LOGFOLDER, file, "energy");
    if( id == 0 ){
      relEneErr->E0inv  = 1.0 / Etot;
      relEneErr->errMax = DBL_MIN;
      fp = fopen(filename, "w");
      if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
      fprintf(fp, "#time(%s)\tsteps\tfile\trelative_error\tEtot(%s)\tEkin(%s)\tEpot(%s)\n", time_astro_unit_name, energy_astro_unit_name, energy_astro_unit_name, energy_astro_unit_name);
    }/* if( id == 0 ){ */
    else{
      fp = fopen(filename, "a");
      if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    }/* else{ */
  }

  const double Eerr = Etot * relEneErr->E0inv - 1.0;
  if( fabs(Eerr) > fabs(relEneErr->errMax) )
    relEneErr->errMax = Eerr;

  if( update ){
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
    fprintf(fp, "%e\t%zu\t%u\t% e\t%e\t%e\t%e\n", time * time2astro, steps, id, Eerr, Etot, Ekin, Epot);
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
    fclose(fp);
  }
#endif//defined(USE_HDF5_FORMAT) && defined(MONITOR_ENERGY_ERROR)

#ifdef  REPORT_COMPUTE_RATE
  writeComputeRate(file, id, time, steps, speed, speed_run, complete, guess, brent_avg, rebuild_interval);
#endif//REPORT_COMPUTE_RATE

  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
#ifndef USE_HDF5_FORMAT
  FILE *fp;
  sprintf(filename, "%s/%s.%s%.3u.dat", DATAFOLDER, file, SNAPSHOT, id);
  fp = fopen(filename, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fwrite(& unit, sizeof(           int), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(& time, sizeof(        double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(&steps, sizeof(         ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.pos, sizeof(    position), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.acc, sizeof(acceleration), tmp, fp) )    success = false;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  tmp = num;  if( tmp != fwrite(body.acc_ext, sizeof(acceleration), tmp, fp) )    success = false;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body.vel, sizeof(velocity), tmp, fp) )    success = false;
#else///BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body.vx, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.vy, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.vz, sizeof(real), tmp, fp) )    success = false;
#endif//BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body.idx, sizeof(ulong), tmp, fp) )    success = false;

#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
  appendGPUclockInfo(monitor_step, deviceMonitors, filename, fp);
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)

  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }

  fclose(fp);
#else///USE_HDF5_FORMAT
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, SNAPSHOT, id);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* create a group and data space */
  char groupname[16];
  sprintf(groupname, SNAPSHOT);
  hid_t group = H5Gcreate(target, groupname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_FILE_COMPRESSION
  hsize_t cdims_max;
  hsize_t cdims_loc[2];
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  const uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  const uint szip_pixels_per_block = 8;
  /* cdims[2] = {128 * szip_pixels_per_block, 3}; */
  hsize_t szip_cdims[2] = {32 * szip_pixels_per_block, 3};
#endif//USE_SZIP_COMPRESSION
#ifdef  USE_GZIP_COMPRESSION
  /* compression using gzip */
  const uint gzip_compress_level = 9;/**< 9 is the maximum compression ratio */
  /* const hsize_t gzip_cdims[2] = {1024, 1}; */
  hsize_t gzip_cdims[2] = {256, 3};
#endif//USE_GZIP_COMPRESSION
#else///USE_FILE_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_FILE_COMPRESSION


  /* 2D (num * 3) array */
  hsize_t dims[2] = {num, 3};
  dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_FILE_COMPRESSION
  cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
#ifdef  USE_SZIP_COMPRESSION
  cdims_loc[0] = szip_cdims[0];
  cdims_loc[1] = szip_cdims[1];
  if( (dims[0] * dims[1]) > (hsize_t)szip_pixels_per_block ){
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] > cdims_max )
      cdims_loc[0] = cdims_max / cdims_loc[1];
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
  }
  else
    property = H5P_DEFAULT;
#else///USE_SZIP_COMPRESSION
#ifdef  USE_GZIP_COMPRESSION
  cdims_loc[0] = gzip_cdims[0];
  cdims_loc[1] = gzip_cdims[1];
  property = H5Pcreate(H5P_DATASET_CREATE);
  if( dims[0] < cdims_loc[0] )
    cdims_loc[0] = dims[0];
  if( cdims_loc[0] * cdims_loc[1] > cdims_max )
    cdims_loc[0] = cdims_max / cdims_loc[1];
  chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
  chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
#endif//USE_SZIP_COMPRESSION
#endif//USE_FILE_COMPRESSION


  /* write particle position */
  dataset = H5Dcreate(group, "position", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->pos));
  chkHDF5err(H5Dclose(dataset));
  /* write particle velocity */
  dataset = H5Dcreate(group, "velocity", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->vel));
  chkHDF5err(H5Dclose(dataset));
  /* write particle acceleration */
  dataset = H5Dcreate(group, "acceleration", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->acc));
  chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* write particle acceleration by external potential field */
  dataset = H5Dcreate(group, "acceleration_external", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->acc_ext));
  chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  USE_FILE_COMPRESSION
#ifdef  USE_SZIP_COMPRESSION
  if( (dims[0] * dims[1]) > (hsize_t)szip_pixels_per_block )
#endif//USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_FILE_COMPRESSION
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));

  /* 1D (num) arrays */
  dataspace = H5Screate_simple(1, dims, NULL);

#ifdef  USE_FILE_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  cdims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
#ifdef  USE_SZIP_COMPRESSION
  szip_cdims[0] = 128 * szip_pixels_per_block;
  cdims_loc[0] = szip_cdims[0];
  if( dims[0] < cdims_loc[0] )
    cdims_loc[0] = dims[0];
  if( cdims_loc[0] > cdims_max )
    cdims_loc[0] = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#else///USE_SZIP_COMPRESSION
#ifdef  USE_GZIP_COMPRESSION
  gzip_cdims[0] = 1024;
  if( dims[0] < cdims_loc[0] )
    cdims_loc[0] = dims[0];
  if( cdims_loc[0] > cdims_max )
    cdims_loc[0] = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
  chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
#endif//USE_SZIP_COMPRESSION
#endif//USE_FILE_COMPRESSION
  /* write particle mass */
  dataset = H5Dcreate(group, "mass", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->m));
  chkHDF5err(H5Dclose(dataset));
  /* write particle potential */
  dataset = H5Dcreate(group, "potential", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->pot));
  chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* write particle potential by external potential field */
  dataset = H5Dcreate(group, "potential_external", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->pot_ext));
  chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  /* write particle index */
  dataset = H5Dcreate(group, "index", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->idx));
  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_FILE_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_FILE_COMPRESSION
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));

  /* write attribute data */
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  /* write current time */
  attribute = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time));
  chkHDF5err(H5Aclose(attribute));
  /* write # of steps */
  attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  /* write # of N-body particles */
  ulong num_ulong = (ulong)num;
  attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &num_ulong));
  chkHDF5err(H5Aclose(attribute));
  /* write flag about USE_DOUBLE_PRECISION */
#ifdef  USE_DOUBLE_PRECISION
  const int useDP = 1;
#else///USE_DOUBLE_PRECISION
  const int useDP = 0;
#endif//USE_DOUBLE_PRECISION
  attribute = H5Acreate(group, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));

  /* close the group */
  chkHDF5err(H5Gclose(group));


  /* write unit system */
  group = H5Gcreate(target, "unit_system", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(group, "unit", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &unit));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "newton", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &newton));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));

  /* write conversion factors */
  hid_t subgroup = H5Gcreate(group, "conversion factors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(subgroup, "length2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &length2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "time2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "mass2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &mass2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "density2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &density2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "col_density2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &col_density2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "energy2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &energy2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "senergy2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &senergy2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "velocity2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &velocity2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "accel2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &accel2astro));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(subgroup));

  /* write axis labels */
  subgroup = H5Gcreate(group, "axis labels", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(subgroup, "length_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, length_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "time_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, time_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "mass_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, mass_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "density_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, density_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "col_density_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, col_density_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "energy_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, energy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "senergy_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, senergy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "velocity_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, velocity_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "accel_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, accel_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(subgroup));

  /* close the group for unit system */
  chkHDF5err(H5Gclose(group));


#ifdef  MONITOR_ENERGY_ERROR
  if( update ){
    /* write energy conservation */
    group = H5Gcreate(target, "energy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    /* create the data space for the attribute */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);

    attribute = H5Acreate(group, "total energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Etot));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "kinetic energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Ekin));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "potential energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Epot));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(group, "relative error", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Eerr));
    chkHDF5err(H5Aclose(attribute));

    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
  }
#endif//MONITOR_ENERGY_ERROR

  /* write GPU information */
#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
  appendGPUclockInfo(monitor_step, deviceMonitors, target, type);
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)


  /* close the file */
  chkHDF5err(H5Fclose(target));

#ifdef  PREPARE_XDMF_FILES
  writeXdmf(num, file, id);
#endif//PREPARE_XDMF_FILES
#endif//USE_HDF5_FORMAT

#   if  defined(USE_HDF5_FORMAT) && defined(MONITOR_ENERGY_ERROR) && !defined(SET_EXTERNAL_POTENTIAL_FIELD)
  if( update )
    if( Eerr > 0.1 ){
      __KILL__(stderr, "detect too large energy error: Eerr = %e\n", Eerr);
    }/* if( Eerr > 0.1 ){ */
#endif//defined(USE_HDF5_FORMAT) && defined(MONITOR_ENERGY_ERROR) && !defined(SET_EXTERNAL_POTENTIAL_FIELD)


  __NOTE__("%s\n", "end");
}


#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
/**
 * @fn writeSnapshotParallel
 *
 * @brief Write snapshot file of the N-body simulation.
 *
 * @param (unit) unit system of the simulation
 * @param (time) current time of the simulation
 * @param (steps) number of time steps calculated in the simulation
 * @param (num) number of N-body particles contained in this MPI process
 * @param (file) name of the simulation
 * @param (id) ID of the snapshot
 * @return (mpi) MPI communicator for data I/O
 * @param (Ntot) total number of N-body particles
 * @param (body) the particle data
 * @param (type) data type for HDF5 (only for HDF5 enabled runs)
 */
void writeSnapshotParallel
(int  unit, double  time, ulong  steps, int num, char file[], uint id, MPIcfg_dataio *mpi, ulong Ntot
#ifdef  USE_HDF5_FORMAT
 , nbody_hdf5 *body, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
 , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
#else///USE_HDF5_FORMAT
 , iparticle body
#endif//USE_HDF5_FORMAT
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
 , gpu_clock *deviceMonitors, const int monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  REPORT_COMPUTE_RATE
 , const double speed, const double speed_run, const double complete, const double guess, const double brent_avg, const double rebuild_interval
#endif//REPORT_COMPUTE_RATE
 )
{
  __NOTE__("%s\n", "start");


  char filename[128];

  /* calculate total energy */
#   if  defined(USE_HDF5_FORMAT) && defined(MONITOR_ENERGY_ERROR)
  double Ekin = 0.0;
  double Epot = 0.0;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  double Eext = 0.0;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  for(int ii = 0; ii < num; ii++){
    const double mass = CAST_R2D(body->m  [ii	     ]);
    const double velx = CAST_R2D(body->vel[ii * 3    ]);
    const double vely = CAST_R2D(body->vel[ii * 3 + 1]);
    const double velz = CAST_R2D(body->vel[ii * 3 + 2]);
    Ekin += mass * (velx * velx + vely * vely + velz * velz);
    Epot += mass * CAST_R2D(body->pot[ii]);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
    Eext += mass * CAST_R2D(body->pot_ext[ii]);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  }/* for(int ii = 0; ii < num; ii++){ */

  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &Ekin, 1, MPI_DOUBLE, MPI_SUM, mpi->comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &Epot, 1, MPI_DOUBLE, MPI_SUM, mpi->comm));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &Eext, 1, MPI_DOUBLE, MPI_SUM, mpi->comm));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  Ekin *= 0.5 * energy2astro;
#ifndef SET_EXTERNAL_POTENTIAL_FIELD
  Epot *= 0.5 * energy2astro;
#else///SET_EXTERNAL_POTENTIAL_FIELD
  Epot = Eext + 0.5 * Epot;
  Epot *= energy2astro;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  const double Etot = Ekin + Epot;

  if( id == 0 ){
    relEneErr->E0inv  = 1.0 / Etot;
    relEneErr->errMax = DBL_MIN;
  }/* if( id == 0 ){ */

  const double Eerr = Etot * relEneErr->E0inv - 1.0;
  if( fabs(Eerr) > fabs(relEneErr->errMax) )
    relEneErr->errMax = Eerr;

  if( mpi->rank == 0 ){
    FILE *fp;
    sprintf(filename, "%s/%s.%s.log", LOGFOLDER, file, "energy");

    if( id == 0 ){
      fp = fopen(filename, "w");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
      fprintf(fp, "#time(%s)\tsteps\tfile\trelative_error\tEtot(%s)\tEkin(%s)\tEpot(%s)\n",
	      time_astro_unit_name, energy_astro_unit_name, energy_astro_unit_name, energy_astro_unit_name);
    }/* if( id == 0 ){ */
    else{
      fp = fopen(filename, "a");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
    }/* else{ */

    fprintf(fp, "%e\t%zu\t%u\t% e\t%e\t%e\t%e\n", time * time2astro, steps, id, Eerr, Etot, Ekin, Epot);
    fclose(fp);
  }/* if( mpi->rank == 0 ){ */
#endif//defined(USE_HDF5_FORMAT) && defined(MONITOR_ENERGY_ERROR)

#ifdef  REPORT_COMPUTE_RATE
  if( mpi->rank == 0 )
    writeComputeRate(file, id, time, steps, speed, speed_run, complete, guess, brent_avg, rebuild_interval);
#endif//REPORT_COMPUTE_RATE


  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
#ifndef USE_HDF5_FORMAT
  /* reset MPI_Offset */
  updateMPIcfg_dataio(mpi, num);

  /* open the target file */
  sprintf(filename, "%s/%s.%s%.3u.dat", DATAFOLDER, file, SNAPSHOT, id);
  MPI_File fh;
  chkMPIerr(MPI_File_open(mpi->comm, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh));
  chkMPIerr(MPI_File_sync(fh));
  chkMPIerr(MPI_File_set_size(fh, 0));
  chkMPIerr(MPI_File_sync(fh));

  MPI_Status status;
  MPI_Offset disp = 0;

  /* the root process writes time, a real value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_INT, MPI_INT, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &unit, 1, MPI_INT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(int);
  /* the root process writes time, a real value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &time, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(double);
  /* the root process writes steps, an unsigned long value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &steps, 1, MPI_UNSIGNED_LONG, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(ulong);
  /* the whole processes write position */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(position), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.pos, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(position);
  /* the whole processes write acceleration */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(acceleration), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.acc, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(acceleration);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* the whole processes write acceleration by external potential field */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(acceleration), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.acc_ext, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(acceleration);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  /* the whole processes write velocity */
#ifdef  BLOCK_TIME_STEP
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(velocity), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.vel, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(velocity);
#else///BLOCK_TIME_STEP
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(real), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.vx, num, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(real);
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(real), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.vy, num, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(real);
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(real), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.vz, num, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(real);
#endif//BLOCK_TIME_STEP
  /* the whole processes write index */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(ulong), MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.idx, num, MPI_UNSIGNED_LONG, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(ulong);

#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
  appendGPUclockInfoParallel(monitor_step, deviceMonitors, *mpi, fh, &disp);
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)

  /* close the target file */
  chkMPIerr(MPI_File_close(&fh));
#else///USE_HDF5_FORMAT
  hid_t f_property = H5Pcreate(H5P_FILE_ACCESS);
  chkHDF5err(H5Pset_fapl_mpio(f_property, mpi->comm, mpi->info));
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, SNAPSHOT, id);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, f_property);
  chkHDF5err(H5Pclose(f_property));

  /* create a group and data space */
  hid_t group = H5Gcreate(target, SNAPSHOT, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


  /* create (distributed) dataset */
  /* create dataspace */
  hsize_t dims_ful[2] = {Ntot            , 3};  hid_t fulSpace = H5Screate_simple(2, dims_ful, NULL);
  hsize_t dims_mem[2] = {Ntot / mpi->size, 3};
  hsize_t dims_loc[2] = { num            , 3};  hid_t locSpace = H5Screate_simple(2, dims_loc, NULL);
  /* create chunked dataset */
  hid_t data_create = H5Pcreate(H5P_DATASET_CREATE);
  hsize_t dims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  if( dims_mem[0] * dims_mem[1] > dims_max )
    dims_mem[0] = dims_max / dims_mem[1];
  chkHDF5err(H5Pset_chunk(data_create, 2, dims_mem));
  hid_t dataset, hyperslab;
  hid_t data_access = H5P_DEFAULT;

  /* configuration about domain decomposition */
  updateMPIcfg_dataio(mpi, num);
  hsize_t  count[2] = {1, 1};
  hsize_t stride[2] = {1, 1};
  hsize_t  block[2] = {dims_loc[0], dims_loc[1]};
  hsize_t offset[2] = {mpi->head, 0};

  /* set up the collective transfer properties function */
  hid_t w_property = H5Pcreate(H5P_DATASET_XFER);
  chkHDF5err(H5Pset_dxpl_mpio(w_property, H5FD_MPIO_COLLECTIVE));

  /* 2D (num, 3) arrays */
  /* write particle position */
  dataset = H5Dcreate(group, "position", type.real, fulSpace, H5P_DEFAULT, data_create, data_access);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->pos));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* write particle velocity */
  dataset = H5Dcreate(group, "velocity", type.real, fulSpace, H5P_DEFAULT, data_create, data_access);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->vel));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* write particle acceleration */
  dataset = H5Dcreate(group, "acceleration", type.real, fulSpace, H5P_DEFAULT, data_create, data_access);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->acc));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* write particle acceleration by external potential field */
  dataset = H5Dcreate(group, "acceleration_external", type.real, fulSpace, H5P_DEFAULT, data_create, data_access);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->acc_ext));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  /* close the dataspace */
  chkHDF5err(H5Pclose(data_access));
  chkHDF5err(H5Pclose(data_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));

  /* 1D (num) arrays */
  fulSpace = H5Screate_simple(1, dims_ful, NULL);
  locSpace = H5Screate_simple(1, dims_loc, NULL);
  data_create = H5Pcreate(H5P_DATASET_CREATE);
  dims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
  if( dims_mem[0] > dims_max )
    dims_mem[0] = dims_max;
  chkHDF5err(H5Pset_chunk(data_create, 1, dims_mem));
  /* write particle mass */
  dataset = H5Dcreate(group, "mass", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->m));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* write particle potential */
  dataset = H5Dcreate(group, "potential", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->pot));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* write particle potential by external potential field */
  dataset = H5Dcreate(group, "potential_external", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->pot_ext));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  /* write particle index */
  dataset = H5Dcreate(group, "index", H5T_NATIVE_ULONG, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_ULONG, locSpace, hyperslab, w_property, body->idx));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* close the dataspace */
  chkHDF5err(H5Pclose(data_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));
  chkHDF5err(H5Pclose(w_property));


  /* write attribute data */
  hsize_t attr_dims = 1;
  hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
  /* write current time */
  hid_t attribute = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time));
  chkHDF5err(H5Aclose(attribute));
  /* write # of steps */
  attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  /* write # of N-body particles */
  attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &Ntot));
  chkHDF5err(H5Aclose(attribute));
  /* write flag about USE_DOUBLE_PRECISION */
#ifdef  USE_DOUBLE_PRECISION
  const int useDP = 1;
#else///USE_DOUBLE_PRECISION
  const int useDP = 0;
#endif//USE_DOUBLE_PRECISION
  attribute = H5Acreate(group, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));

  /* close the group */
  chkHDF5err(H5Gclose(group));


  /* write unit system */
  group = H5Gcreate(target, "unit_system", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(group, "unit", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &unit));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "newton", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &newton));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Sclose(dataspace));

  /* write conversion factors */
  hid_t subgroup = H5Gcreate(group, "conversion factors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(subgroup, "length2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &length2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "time2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "mass2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &mass2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "density2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &density2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "col_density2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &col_density2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "energy2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &energy2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "senergy2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &senergy2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "velocity2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &velocity2astro));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "accel2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &accel2astro));
  chkHDF5err(H5Aclose(attribute));

  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(subgroup));


  /* write axis labels */
  subgroup = H5Gcreate(group, "axis labels", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(subgroup, "length_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, length_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "time_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, time_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "mass_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, mass_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "density_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, density_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "col_density_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, col_density_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "energy_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, energy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "senergy_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, senergy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "velocity_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, velocity_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(subgroup, "accel_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, accel_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));

  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(subgroup));
  /* close the group for unit system */
  chkHDF5err(H5Gclose(group));


  /* write energy conservation */
#ifdef  MONITOR_ENERGY_ERROR
  group = H5Gcreate(target, "energy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);

  attribute = H5Acreate(group, "total energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Etot));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "kinetic energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Ekin));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "potential energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Epot));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "relative error", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Eerr));
  chkHDF5err(H5Aclose(attribute));

  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(group));
#endif//MONITOR_ENERGY_ERROR

  /* write GPU information */
#   if  defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)
  appendGPUclockInfoParallel(monitor_step, deviceMonitors, *mpi, target, type);
#endif//defined(REPORT_GPU_CLOCK_FREQUENCY) && !defined(RUN_WITHOUT_GOTHIC)


  /* close the file */
  chkHDF5err(H5Fclose(target));

#ifdef  PREPARE_XDMF_FILES
  if( mpi->rank == 0 )
    writeXdmf(Ntot, file, id);
#endif//PREPARE_XDMF_FILES
#endif//USE_HDF5_FORMAT

#   if  defined(USE_HDF5_FORMAT) && defined(MONITOR_ENERGY_ERROR) && !defined(SET_EXTERNAL_POTENTIAL_FIELD)
  if( Eerr > 0.1 ){
    __KILL__(stderr, "detect too large energy error: Eerr = %e\n", Eerr);
  }/* if( Eerr > 0.1 ){ */
#endif//defined(USE_HDF5_FORMAT) && defined(MONITOR_ENERGY_ERROR) && !defined(SET_EXTERNAL_POTENTIAL_FIELD)


  __NOTE__("%s\n", "end");
}
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)


#ifdef  USE_HDF5_FORMAT
#ifdef  PREPARE_XDMF_FILES
static inline void writeXdmf_split(char file[], const uint id, const int kind, int *num)
{
  __NOTE__("%s\n", "start");

  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.xmf", DATAFOLDER, file, "split", id);
  FILE *fp;
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  sprintf(filename, "%s.%s%.3u.h5", file, "split", id);

  fprintf(fp, "<?xml version=\"1.0\" ?>\n");
  fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(fp, "<Xdmf Version=\"3.0\">\n");
  fprintf(fp, "<Domain>\n");
  for(int kk = 0; kk < kind; kk++){
    fprintf(fp, "<Grid Name=\"data%d\" GridType=\"Uniform\">\n", kk);
    fprintf(fp, "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n", num[kk]);
    fprintf(fp, "<Geometry GeometryType=\"XYZ\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", num[kk], sizeof(real));
    fprintf(fp, "%s:/data%d/position\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "</Geometry>\n");
    fprintf(fp, "<Attribute Name=\"velocity\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", num[kk], sizeof(real));
    fprintf(fp, "%s:/data%d/velocity\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "</Attribute>\n");
    fprintf(fp, "<Attribute Name=\"acceleration\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", num[kk], sizeof(real));
    fprintf(fp, "%s:/data%d/acceleration\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "</Attribute>\n");
    fprintf(fp, "<Attribute Name=\"mass\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", num[kk], sizeof(real));
    fprintf(fp, "%s:/data%d/mass\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "</Attribute>\n");
    fprintf(fp, "<Attribute Name=\"potential\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", num[kk], sizeof(real));
    fprintf(fp, "%s:/data%d/potential\n", filename, kk);
    fprintf(fp, "</DataItem>\n");
    fprintf(fp, "</Attribute>\n");
    fprintf(fp, "<Attribute Name=\"index\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(fp, "<DataItem Dimensions=\"%d\" NumberType=\"UInt\" Precision=\"%zu\" Format=\"HDF\">\n", num[kk], sizeof(ulong));
    fprintf(fp, "%s:/data%d/index\n", filename, kk);
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
 * @fn writeSnapshotMultiGroups
 *
 * @brief Write snapshot file of the N-body simulation, which is consisting of multiple components.
 *
 * @param (time) current time of the simulation
 * @param (steps) number of time steps calculated in the simulation
 * @param (*body) the particle data
 * @param (file) name of the simulation
 * @param (id) ID of the snapshot
 * @param (type) data type for HDF5
 * @param (kind) number of particle groups
 * @return (head) index of the head particle in each group
 * @return (num) number of N-body particles in each group
 */
#   if  !defined(USE_SZIP_COMPRESSION) && !defined(USE_GZIP_COMPRESSION)
#define USE_GZIP_COMPRESSION
#define TENTATIVE_USE_GZIP_COMPRESSION
#endif//!defined(USE_SZIP_COMPRESSION) && !defined(USE_GZIP_COMPRESSION)
void writeSnapshotMultiGroups(double  time, ulong  steps, nbody_hdf5 *body, char file[], uint id, hdf5struct type, int kind, int *head, int *num)
{
  __NOTE__("%s\n", "start");


  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "split", id);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


  /* create the data space for the dataset */
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
  const hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  hsize_t cdims_loc[2];
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  const uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  const uint szip_pixels_per_block = 8;
  hsize_t szip_cdims[2] = {32 * szip_pixels_per_block, 3};
#endif//USE_SZIP_COMPRESSION
#ifdef  USE_GZIP_COMPRESSION
  /* compression using gzip */
  const uint gzip_compress_level = 9;/**< 9 is the maximum compression ratio */
  /* const hsize_t gzip_cdims[2] = {1, 1024}; */
  hsize_t gzip_cdims[2] = {512, 3};
#endif//USE_GZIP_COMPRESSION
#else///defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
  property = H5P_DEFAULT;
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)


  /* write attribute data */
  /* create the data space for the attribute */
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
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
    num_ulong += (ulong)num[ii];
  attribute = H5Acreate(target, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &num_ulong));
  chkHDF5err(H5Aclose(attribute));
  /* write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &kind));
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

  attribute = H5Acreate(target, "length_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, length_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "time_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, time_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "mass_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, mass_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "senergy_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, senergy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "velocity_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, velocity_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "accel_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, accel_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));


  chkHDF5err(H5Tclose(str4format));
  chkHDF5err(H5Sclose(dataspace));


  /* write particle data */
  for(int ii = 0; ii < kind; ii++){
    char grp[16];    sprintf(grp, "data%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* 2D (num * 3) array */
    hsize_t dims[2] = {num[ii], 3};
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
    for(int jj = 3 * head[ii]; jj < 3 * (head[ii] + num[ii]); jj++)
      body->pos[jj] = CAST_D2R(CAST_R2D(body->pos[jj]) * length2astro);
    /* coordinate transformation */
    dataset = H5Dcreate(group, "position", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->pos[head[ii] * 3]));
    chkHDF5err(H5Dclose(dataset));
    /* write particle velocity */
    /* coordinate transformation */
    for(int jj = 3 * head[ii]; jj < 3 * (head[ii] + num[ii]); jj++)
      body->vel[jj] = CAST_D2R(CAST_R2D(body->vel[jj]) * velocity2astro);
    dataset = H5Dcreate(group, "velocity", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->vel[head[ii] * 3]));
    chkHDF5err(H5Dclose(dataset));
    /* write particle acceleration */
    /* coordinate transformation */
    for(int jj = 3 * head[ii]; jj < 3 * (head[ii] + num[ii]); jj++)
      body->acc[jj] = CAST_D2R(CAST_R2D(body->acc[jj]) * accel2astro);
    dataset = H5Dcreate(group, "acceleration", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->acc[head[ii] * 3]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
    /* coordinate transformation */
    for(int jj = 3 * head[ii]; jj < 3 * (head[ii] + num[ii]); jj++)
      body->acc_ext[jj] = CAST_D2R(CAST_R2D(body->acc_ext[jj]) * accel2astro);
    dataset = H5Dcreate(group, "acceleration_external", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->acc_ext[head[ii] * 3]));
    chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
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
    /* write particle mass */
    for(int jj = head[ii]; jj < head[ii] + num[ii]; jj++)
      body->m[jj] = CAST_D2R(CAST_R2D(body->m[jj]) * mass2astro);
    dataset = H5Dcreate(group, "mass", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->m[head[ii]]));
    chkHDF5err(H5Dclose(dataset));
    /* write particle potential */
    for(int jj = head[ii]; jj < head[ii] + num[ii]; jj++)
      body->pot[jj] = CAST_D2R(CAST_R2D(body->pot[jj]) * senergy2astro);
    dataset = H5Dcreate(group, "potential", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->pot[head[ii]]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
    /* write particle potential by external potential field */
    for(int jj = head[ii]; jj < head[ii] + num[ii]; jj++)
      body->pot_ext[jj] = CAST_D2R(CAST_R2D(body->pot_ext[jj]) * senergy2astro);
    dataset = H5Dcreate(group, "potential_external", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->pot_ext[head[ii]]));
    chkHDF5err(H5Dclose(dataset));
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    /* write particle index */
    dataset = H5Dcreate(group, "index", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->idx[head[ii]]));
    chkHDF5err(H5Dclose(dataset));
#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#ifdef  USE_SZIP_COMPRESSION
    if( dims[0] > (hsize_t)szip_pixels_per_block )
#endif//USE_SZIP_COMPRESSION
      chkHDF5err(H5Pclose(property));
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));


    /* write attribute data */
    /* create the data space for the attribute */
    hsize_t attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    /* write # of N-body particles */
    ulong num_ulong = (ulong)num[ii];
    attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &num_ulong));
    chkHDF5err(H5Aclose(attribute));
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
  }/* for(int ii = 0; ii < kind; ii++){ */


  /* close the file */
  chkHDF5err(H5Fclose(target));
#ifdef  PREPARE_XDMF_FILES
  writeXdmf_split(file, id, kind, num);
#endif//PREPARE_XDMF_FILES

  __NOTE__("%s\n", "end");
}
#   if  defined(TENTATIVE_USE_GZIP_COMPRESSION) && defined(USE_GZIP_COMPRESSION)
#undef TENTATIVE_USE_GZIP_COMPRESSION
#undef USE_GZIP_COMPRESSION
#endif//defined(TENTATIVE_USE_GZIP_COMPRESSION) && defined(USE_GZIP_COMPRESSION)
#endif//USE_HDF5_FORMAT


/**
 * @fn readApproxAccel
 *
 * @brief Read particles' acceleration calculated by GOTHIC
 *
 * @return (time) current time of the simulation
 * @return (step) number of time steps calculated in the simulation
 * @param (num) number of N-body particles
 * @return (idx) index of the particles
 * @return (acc) acceleration of the particles
 * @param (file) name of the simulation
 */
void  readApproxAccel(double *time, ulong *step, int num, ulong *idx, acceleration *acc, char file[])
{
  __NOTE__("%s\n", "start");


  FILE *fp;
  fp = fopen(file, "rb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", file);  }

  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fread(time, sizeof(      double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(step, sizeof(       ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread( idx, sizeof(       ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread( acc, sizeof(acceleration), tmp, fp) )    success = false;

  if( success != true ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", file);  }
  fclose(fp);


  __NOTE__("%s\n", "end");
}
/**
 * @fn writeApproxAccel
 *
 * @brief Write particles' acceleration calculated by GOTHIC
 *
 * @param (time) current time of the simulation
 * @param (step) number of time steps calculated in the simulation
 * @param (num) number of N-body particles
 * @param (idx) index of the particles
 * @param (acc) acceleration of the particles
 * @param (file) name of the simulation
 */
void writeApproxAccel(double  time, ulong  step, int num, ulong *idx, acceleration *acc, char file[])
{
  __NOTE__("%s\n", "start");


  FILE *fp;
  fp = fopen(file, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", file);  }

  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fwrite(&time, sizeof(      double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(&step, sizeof(       ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(  idx, sizeof(       ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(  acc, sizeof(acceleration), tmp, fp) )    success = false;

  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", file);  }
  fclose(fp);


  __NOTE__("%s\n", "end");
}


#ifdef  EXEC_BENCHMARK
/**
 * @fn dumpBenchmark
 *
 * @brief Write results of measured execution time of GOTHIC.
 *
 * @param (jobID) index of the simulation run
 * @param (file) name of the simulation
 * @param (steps) number of time steps calculated in the simulation
 * @param (dat) measured execution time
 */
void dumpBenchmark(int jobID, char file[], int steps, wall_clock_time *dat)
{
  __NOTE__("%s\n", "start");


  char filename[128];
  FILE *fp;

#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  wall_clock_time mean = {0};
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)


  sprintf(filename, "%s/%s.%s%.8d.%s.dat", BENCH_LOG_FOLDER, file, WALLCLOCK, jobID, "bare");
  int newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  if( newFile ){
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", "step",
	    "calcGrav_dev", "calcMom_dev", "makeTree", "setTimeStep_dev", "sortPHcurve",
	    "cpBody_dev2hst", "cpBody_hst2dev");
#ifdef  BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s\t%s\t%s", "prediction_dev", "correction_dev", "setLaneTime_dev", "adjustTime_dev");
#else///BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s", "advPos_dev", "advVel_dev");
#endif//BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s", "examineNeighbor_dev", "searchNeighbor_dev");
#ifdef  HUNT_MAKE_PARAMETER
    fprintf(fp, "\t%s\t%s\t%s\t%s\t%s\t%s", "genPHkey_kernel", "rsortKey_library", "sortBody_kernel", "makeTree_kernel", "linkTree_kernel", "trimTree_kernel");
    fprintf(fp, "\t%s\t%s\t%s\t%s\t%s", "initTreeLink_kernel", "initTreeCell_kernel", "initTreeNode_kernel", "initTreeBody_kernel", "copyRealBody_kernel");
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
    fprintf(fp, "\t%s\t%s\t%s\t%s", "searchNeighbor_kernel", "sortNeighbors", "countNeighbors_kernel", "commitNeighbors");
#endif//HUNT_FIND_PARAMETER
    fprintf(fp, "\n");
  }/* if( newFile ){ */
  for(int ii = 0; ii < steps; ii++){
    fprintf(fp, "%4d\t%e\t%e\t%e\t%e\t%e\t%e\t%e", ii,
	    dat[ii].calcGravity_dev, dat[ii].calcMultipole_dev, dat[ii].makeTree, dat[ii].setTimeStep_dev, dat[ii].sortParticlesPHcurve,
	    dat[ii].copyParticle_dev2hst, dat[ii].copyParticle_hst2dev);
#ifdef  BLOCK_TIME_STEP
    fprintf(fp, "\t%e\t%e\t%e\t%e", dat[ii].prediction_dev, dat[ii].correction_dev, dat[ii].setLaneTime_dev, dat[ii].adjustParticleTime_dev);
#else///BLOCK_TIME_STEP
    fprintf(fp, "\t%e\t%e", dat[ii].advPos_dev, dat[ii].advVel_dev);
#endif//BLOCK_TIME_STEP
    fprintf(fp, "\t%e\t%e", dat[ii].examineNeighbor_dev, dat[ii].searchNeighbor_dev);
#ifdef  HUNT_MAKE_PARAMETER
    fprintf(fp, "\t%e\t%e\t%e\t%e\t%e\t%e", dat[ii].genPHkey_kernel, dat[ii].rsortKey_library, dat[ii].sortBody_kernel, dat[ii].makeTree_kernel, dat[ii].linkTree_kernel, dat[ii].trimTree_kernel);
    fprintf(fp, "\t%e\t%e\t%e\t%e\t%e", dat[ii].initTreeLink_kernel, dat[ii].initTreeCell_kernel, dat[ii].initTreeNode_kernel, dat[ii].initTreeBody_kernel, dat[ii].copyRealBody_kernel);
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
    fprintf(fp, "\t%e\t%e\t%e\t%e", dat[ii].searchNeighbor_kernel, dat[ii].sortNeighbors, dat[ii].countNeighbors_kernel, dat[ii].commitNeighbors);
#endif//HUNT_FIND_PARAMETER
    fprintf(fp, "\n");
    mean.calcGravity_dev      += dat[ii].calcGravity_dev;
    mean.calcMultipole_dev    += dat[ii].calcMultipole_dev;
    mean.makeTree             += dat[ii].makeTree;
#ifdef  BLOCK_TIME_STEP
    mean.prediction_dev         += dat[ii].prediction_dev;
    mean.correction_dev         += dat[ii].correction_dev;
    mean.setLaneTime_dev        += dat[ii].setLaneTime_dev;
    mean.adjustParticleTime_dev += dat[ii].adjustParticleTime_dev;
#else///BLOCK_TIME_STEP
    mean.advPos_dev           += dat[ii].advPos_dev;
    mean.advVel_dev           += dat[ii].advVel_dev;
#endif//BLOCK_TIME_STEP
    mean.setTimeStep_dev      += dat[ii].setTimeStep_dev;
    mean.sortParticlesPHcurve += dat[ii].sortParticlesPHcurve;
    mean.copyParticle_dev2hst += dat[ii].copyParticle_dev2hst;
    mean.copyParticle_hst2dev += dat[ii].copyParticle_hst2dev;
    mean.examineNeighbor_dev += dat[ii].examineNeighbor_dev;
    mean. searchNeighbor_dev += dat[ii]. searchNeighbor_dev;
#ifdef  HUNT_MAKE_PARAMETER
    mean.genPHkey_kernel     += dat[ii].genPHkey_kernel;
    mean.rsortKey_library    += dat[ii].rsortKey_library;
    mean.sortBody_kernel     += dat[ii].sortBody_kernel;
    mean.makeTree_kernel     += dat[ii].makeTree_kernel;
    mean.linkTree_kernel     += dat[ii].linkTree_kernel;
    mean.trimTree_kernel     += dat[ii].trimTree_kernel;
    mean.initTreeLink_kernel += dat[ii].initTreeLink_kernel;
    mean.initTreeCell_kernel += dat[ii].initTreeCell_kernel;
    mean.initTreeNode_kernel += dat[ii].initTreeNode_kernel;
    mean.initTreeBody_kernel += dat[ii].initTreeBody_kernel;
    mean.copyRealBody_kernel += dat[ii].copyRealBody_kernel;
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
    mean.searchNeighbor_kernel += dat[ii].searchNeighbor_kernel;
    mean.sortNeighbors         += dat[ii].sortNeighbors;
    mean.countNeighbors_kernel += dat[ii].countNeighbors_kernel;
    mean.commitNeighbors       += dat[ii].commitNeighbors;
#endif//HUNT_FIND_PARAMETER
  }
  fclose(fp);


  sprintf(filename, "%s/%s.%s%.8d.%s.dat", BENCH_LOG_FOLDER, file, WALLCLOCK, jobID, "mean");
  newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  if( newFile ){
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s",
	    "calcGrav_dev", "calcMom_dev", "makeTree", "setTimeStep_dev", "sortPHcurve",
	    "cpBody_dev2hst", "cpBody_hst2dev");
#ifdef  BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s\t%s\t%s", "prediction_dev", "correction_dev", "setLaneTime_dev", "adjustTime_dev");
#else///BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s", "advPos_dev", "advVel_dev");
#endif//BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s", "examineNeighbor_dev", "searchNeighbor_dev");
#ifdef  HUNT_MAKE_PARAMETER
    fprintf(fp, "\t%s\t%s\t%s\t%s\t%s\t%s", "genPHkey_kernel", "rsortKey_library", "sortBody_kernel", "makeTree_kernel", "linkTree_kernel", "trimTree_kernel");
    fprintf(fp, "\t%s\t%s\t%s\t%s\t%s", "initTreeLink_kernel", "initTreeCell_kernel", "initTreeNode_kernel", "initTreeBody_kernel", "copyRealBody_kernel");
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
    fprintf(fp, "\t%s\t%s\t%s\t%s", "searchNeighbor_kernel", "sortNeighbors", "countNeighbors_kernel", "commitNeighbors");
#endif//HUNT_FIND_PARAMETER
    fprintf(fp, "\n");
  }
  const double inv = 1.0 / (double)steps;
  fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e",
	  inv * mean.calcGravity_dev, inv * mean.calcMultipole_dev, inv * mean.makeTree, inv * mean.setTimeStep_dev, inv * mean.sortParticlesPHcurve,
	  inv * mean.copyParticle_dev2hst, inv * mean.copyParticle_hst2dev);
#ifdef  BLOCK_TIME_STEP
  fprintf(fp, "\t%e\t%e\t%e\t%e", inv * mean.prediction_dev, inv * mean.correction_dev, inv * mean.setLaneTime_dev, inv * mean.adjustParticleTime_dev);
#else///BLOCK_TIME_STEP
  fprintf(fp, "\t%e\t%e", inv * mean.advPos_dev, inv * mean.advVel_dev);
#endif//BLOCK_TIME_STEP
  fprintf(fp, "\t%e\t%e", inv * mean.examineNeighbor_dev, inv * mean.searchNeighbor_dev);
#ifdef  HUNT_MAKE_PARAMETER
  fprintf(fp, "\t%e\t%e\t%e\t%e\t%e\t%e", inv * mean.genPHkey_kernel, inv * mean.rsortKey_library, inv * mean.sortBody_kernel, inv * mean.makeTree_kernel, inv * mean.linkTree_kernel, inv * mean.trimTree_kernel);
  fprintf(fp, "\t%e\t%e\t%e\t%e\t%e", inv * mean.initTreeLink_kernel, inv * mean.initTreeCell_kernel, inv * mean.initTreeNode_kernel, inv * mean.initTreeBody_kernel, inv * mean.copyRealBody_kernel);
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
  fprintf(fp, "\t%e\t%e\t%e\t%e", inv * mean.searchNeighbor_kernel, inv * mean.sortNeighbors, inv * mean.countNeighbors_kernel, inv * mean.commitNeighbors);
#endif//HUNT_FIND_PARAMETER
  fprintf(fp, "\n");
  fclose(fp);

#ifdef  HUNT_WALK_PARAMETER
  sprintf(filename, "%s/%s.%s.dat", BENCH_LOG_FOLDER, file, "walk");
  newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  if( newFile ){
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", "calcGravity_dev", "NTHREADS", "TSUB", "NWARP", "NLOOP", "NUNROLL",
	    "USE_WARP_SHUFFLE_FUNC", "GET_BUFID",
	    "NBLOCKS_PER_SM", "USE_WARP_REDUCE_FUNC", "USE_L2_SET_ASIDE_POLICY", "NLEVEL_TREE_NODE_L2_PERSISTING");
#ifdef  NEIGHBOR_PHKEY_LEVEL
    fprintf(fp, "\t%s", "NEIGHBOR_PHKEY_LEVEL");
#endif//NEIGHBOR_PHKEY_LEVEL
    fprintf(fp, "\t%s", "NEIGHBOR_LENGTH_SHRINK_FACTOR");
    fprintf(fp, "\n");
  }/* if( newFile ){ */
  fprintf(fp, "%e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
	  inv * mean.calcGravity_dev,
	  NTHREADS, TSUB, NWARP, NLOOP, NUNROLL
	  /* use warp shuffle instruction or not */
#ifdef  USE_WARP_SHUFFLE_FUNC
	  , 1
#else///USE_WARP_SHUFFLE_FUNC
	  , 0
#endif//USE_WARP_SHUFFLE_FUNC
	  /* how to occupy buffer on global memory */
#ifdef  USE_SMID_TO_GET_BUFID
	  , 1
#else///USE_SMID_TO_GET_BUFID
#ifdef  TRY_MODE_ABOUT_BUFFER
	  , 2
#else///TRY_MODE_ABOUT_BUFFER
	  , 0
#endif//TRY_MODE_ABOUT_BUFFER
#endif//USE_SMID_TO_GET_BUFID
	  , NBLOCKS_PER_SM
#ifdef  USE_WARP_REDUCE_FUNC
	  , 1
#else///USE_WARP_REDUCE_FUNC
	  , 0
#endif//USE_WARP_REDUCE_FUNC
#ifdef  USE_L2_SET_ASIDE_POLICY
	  , 1, NLEVEL_TREE_NODE_L2_PERSISTING
#else///USE_L2_SET_ASIDE_POLICY
	  , 0, 0
#endif//USE_L2_SET_ASIDE_POLICY
	  );
  /* criterion about PH-key level to split i-particles */
#ifdef  NEIGHBOR_PHKEY_LEVEL
  fprintf(fp, "\t%d", NEIGHBOR_PHKEY_LEVEL);
#endif//NEIGHBOR_PHKEY_LEVEL
  /* criterion about neighbor length to split i-particles */
  fprintf(fp, "\t%e", NEIGHBOR_LENGTH_SHRINK_FACTOR);
  fprintf(fp, "\n");
  fclose(fp);
#endif//HUNT_WALK_PARAMETER

#   if  defined(HUNT_MAKE_PARAMETER) || defined(HUNT_FIND_PARAMETER)
  extern int treeBuildCalls;
  const double invSteps = 1.0 / (double)treeBuildCalls;
#endif//defined(HUNT_MAKE_PARAMETER) || defined(HUNT_FIND_PARAMETER)

#ifdef  HUNT_MAKE_PARAMETER
  sprintf(filename, "%s/%s.%s.dat", BENCH_LOG_FOLDER, file, "make");
  newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  if( newFile ){
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s", "makeTree", "sortParticlesPHcurve", "genPHkey_kernel", "rsortKey_library", "NTHREADS_PH", "sortBody_kernel", "NTHREADS_PHSORT");
  fprintf(fp, "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
	  "makeTree_kernel", "linkTree_kernel", "trimTree_kernel",
	  "NTHREADS_MAKE_TREE", "TSUB_MAKE_TREE", "NTHREADS_LINK_TREE", "NTHREADS_TRIM_TREE",
	  "USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE", "USE_WARP_REDUCE_FUNC_MAKE_TREE_STRUCTURE");
  fprintf(fp, "\t%s\t%s\t%s\t%s\t%s\t%s", "initTreeLink_kernel", "initTreeCell_kernel", "initTreeNode_kernel", "initTreeBody_kernel", "copyRealBody_kernel", "Ttot");
  fprintf(fp, "\n");
  }/* if( newFile ){ */
  fprintf(fp, "%e\t%e\t%e\t%e\t%d\t%e\t%d", invSteps * mean.makeTree, invSteps * mean.sortParticlesPHcurve,
	  invSteps * mean.genPHkey_kernel, invSteps * mean.rsortKey_library, NTHREADS_PH, invSteps * mean.sortBody_kernel, NTHREADS_PHSORT);
  fprintf(fp, "\t%e\t%e\t%e\t%d\t%d\t%d\t%d\t%d\t%d",
	  invSteps * mean.makeTree_kernel, invSteps * mean.linkTree_kernel, invSteps * mean.trimTree_kernel,
	  NTHREADS_MAKE_TREE, TSUB_MAKE_TREE, NTHREADS_LINK_TREE, NTHREADS_TRIM_TREE
	  /* use warp shuffle for making tree structure instruction or not */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
	  , 1
#else///USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
	  , 0
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
#ifdef  USE_WARP_REDUCE_FUNC_MAKE_TREE_STRUCTURE
	  , 1
#else///USE_WARP_REDUCE_FUNC_MAKE_TREE_STRUCTURE
	  , 0
#endif//USE_WARP_REDUCE_FUNC_MAKE_TREE_STRUCTURE
	  );
  fprintf(fp, "\t%e\t%e\t%e\t%e\t%e\t%d", invSteps * mean.initTreeLink_kernel, invSteps * mean.initTreeCell_kernel, invSteps * mean.initTreeNode_kernel, inv * mean.initTreeBody_kernel, inv * mean.copyRealBody_kernel, NTHREADS_INIT_BODY);
  fprintf(fp, "\n");
  fclose(fp);
#endif//HUNT_MAKE_PARAMETER

#ifdef  HUNT_NODE_PARAMETER
  sprintf(filename, "%s/%s.%s.dat", BENCH_LOG_FOLDER, file, "node");
  newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  if( newFile )
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\n",
	    "calcMultipole", "NTHREADS_MAC", "TSUB_MAC", "USE_WARP_SHUFFLE_FUNC_MAC", "USE_WARP_REDUCE_FUNC_MAC");
  fprintf(fp, "%e\t%d\t%d\t%d\t%d\n",
	  inv * mean.calcMultipole_dev,
	  NTHREADS_MAC, TSUB_MAC
	  /* use warp shuffle instruction or not */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
	  , 1
#else///USE_WARP_SHUFFLE_FUNC_MAC
	  , 0
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#ifdef  USE_WARP_REDUCE_FUNC_MAC
	  , 1
#else///USE_WARP_REDUCE_FUNC_MAC
	  , 0
#endif//USE_WARP_REDUCE_FUNC_MAC
	  );
  fclose(fp);
#endif//HUNT_NODE_PARAMETER

#ifdef  HUNT_TIME_PARAMETER
  sprintf(filename, "%s/%s.%s.dat", BENCH_LOG_FOLDER, file, "time");
  newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  if( newFile ){
#ifdef  BLOCK_TIME_STEP
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s", "prediction_dev", "correction_dev", "setLaneTime_dev", "adjustParticleTime_dev", "setTimeStep_dev");
#else///BLOCK_TIME_STEP
    fprintf(fp, "#%s\t%s\t%s", "advPos_dev", "advVel_dev", "setTimeStep_dev");
#endif//BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s\t%s\n", "NTHREADS_TIME", "USE_WARP_SHUFFLE_FUNC_TIME", "USE_WARP_REDUCE_FUNC_TIME");
  }/* if( newFile ){ */
#ifdef  BLOCK_TIME_STEP
  fprintf(fp, "%e\t%e\t%e\t%e\t%e",
	  inv * mean.prediction_dev, inv * mean.correction_dev, inv * mean.setLaneTime_dev,
	  inv * mean.adjustParticleTime_dev, inv * mean.setTimeStep_dev);
#else///BLOCK_TIME_STEP
  fprintf(fp, "%e\t%e\t%e", inv * mean.advPos_dev, inv * mean.advVel_dev, inv * mean.setTimeStep_dev);
#endif//BLOCK_TIME_STEP
  fprintf(fp, "\t%d\t%d\t%d\n", NTHREADS_TIME
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
	  , 1
#else///USE_WARP_SHUFFLE_FUNC_TIME
	  , 0
#endif//USE_WARP_SHUFFLE_FUNC_TIME
#ifdef  USE_WARP_REDUCE_FUNC_TIME
	  , 1
#else///USE_WARP_REDUCE_FUNC_TIME
	  , 0
#endif//USE_WARP_REDUCE_FUNC_TIME
	  );
  fclose(fp);
#endif//HUNT_TIME_PARAMETER

#ifdef  HUNT_FIND_PARAMETER
  sprintf(filename, "%s/%s.%s.dat", BENCH_LOG_FOLDER, file, "find");
  newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  if( newFile )
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
	    , "examineNeighbor_dev", "searchNeighbor_dev"
	    , "facileNeighborSearching_kernel", "NTHREADS_FACILE_NS", "NEIGHBOR_NUM", "SMEM_PREF_FOR_NEIGHBOR_SEARCH"
	    , "sortNeighbors"
	    , "countContinuousNeighbors_kernel", "NTHREADS_SHRINK"
	    , "commitNeighbors");
  fprintf(fp, "%e\t%e\t%e\t%d\t%d\t%d\t%e\t%e\t%d\t%e\n"
	  , invSteps * mean.examineNeighbor_dev, invSteps * mean.searchNeighbor_dev
	  , invSteps * mean.searchNeighbor_kernel, NTHREADS_FACILE_NS, NEIGHBOR_NUM
#ifdef  SMEM_PREF_FOR_NEIGHBOR_SEARCH
	  , 1
#else///SMEM_PREF_FOR_NEIGHBOR_SEARCH
	  , 0
#endif//SMEM_PREF_FOR_NEIGHBOR_SEARCH
	  , invSteps * mean.sortNeighbors
	  , invSteps * mean.countNeighbors_kernel, NTHREADS_SHRINK
	  , invSteps * mean.commitNeighbors);
  fclose(fp);
#endif//HUNT_FIND_PARAMETER


  __NOTE__("%s\n", "end");
}
#endif//EXEC_BENCHMARK


#ifdef  COUNT_INTERACTIONS
/**
 * @fn dumpStatistics
 *
 * @brief Write # of interactions in GOTHIC.
 *
 * @param (jobID) index of the simulation run
 * @param (file) name of the simulation
 * @param (steps) number of time steps calculated in the simulation
 * @param (PHlevel) level of the PH-key hierarchy in the simulation
 * @param (dat) measured # of interactions
 */
void dumpStatistics(int jobID, char file[], int steps, int PHlevel, tree_metrics *dat)
{
  __NOTE__("%s\n", "start");


  char filename[128];
  FILE *fp;

  sprintf(filename, "%s/%s.%s%.8d.%s.dat", BENCH_LOG_FOLDER, file, WALK_STAT, jobID, "Nj");
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	  "step",
	  "interactions",
	  "Ngroups",
	  "minimum",
	  "maximum",
	  "mean",
	  "dispersion");
  for(int ii = 0; ii < steps; ii++)
    fprintf(fp, "%d\t%zu\t%d\t%d\t%d\t%e\t%e\n",
	    ii,
	    dat[ii].Ninteractions,
	    dat[ii].Nj.Ngroup,
	    dat[ii].Nj.min,
	    dat[ii].Nj.max,
	    dat[ii].Nj.mean,
	    dat[ii].Nj.sdev);

  fclose(fp);


  sprintf(filename, "%s/%s.%s%.8d.%s.dat", BENCH_LOG_FOLDER, file, WALK_STAT, jobID, "Nbuf");
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	  "step",
	  "bufSize",
	  "Ngroups",
	  "minimum",
	  "maximum",
	  "mean",
	  "dispersion");
  for(int ii = 0; ii < steps; ii++)
    fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%e\t%e\n",
	    ii,
	    dat[ii].bufSize,
	    dat[ii].Nbuf.Ngroup,
	    dat[ii].Nbuf.min,
	    dat[ii].Nbuf.max,
	    dat[ii].Nbuf.mean,
	    dat[ii].Nbuf.sdev);
  fclose(fp);


  sprintf(filename, "%s/%s.%s%.8d.%s.dat", BENCH_LOG_FOLDER, file, TREE_STAT, jobID, "total");
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	  "step",
	  "cellNum",
	  "nodeNum",
	  "mjMean",
	  "mjSdev",
	  "r2Mean",
	  "r2Sdev");
  for(int ii = 0; ii < steps; ii++)
    fprintf(fp, "%d\t%d\t%d\t%e\t%e\t%e\t%e\n",
	    ii,
	    dat[ii].total.cellNum,
	    dat[ii].total.nodeNum,
	    dat[ii].total.mjMean,
	    dat[ii].total.mjSdev,
	    dat[ii].total.r2Mean,
	    dat[ii].total.r2Sdev);
  fclose(fp);


  for(int ll = 0; ll < PHlevel; ll++){
    sprintf(filename, "%s/%s.%s%.8d.%s%.2d.dat", BENCH_LOG_FOLDER, file, TREE_STAT, jobID, "level", ll);
    fp = fopen(filename, "w");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }

    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	    "step",
	    "cellNum",
	    "nodeNum",
	    "mjMean",
	    "mjSdev",
	    "r2Mean",
	    "r2Sdev");
    for(int ii = 0; ii < steps; ii++)
      fprintf(fp, "%d\t%d\t%d\t%e\t%e\t%e\t%e\n",
	      ii,
	      dat[ii].level[ll].cellNum,
	      dat[ii].level[ll].nodeNum,
	      dat[ii].level[ll].mjMean,
	      dat[ii].level[ll].mjSdev,
	      dat[ii].level[ll].r2Mean,
	      dat[ii].level[ll].r2Sdev);

    fclose(fp);
  }/* for(int ll = 0; ll < PHlevel; ll++){ */


  __NOTE__("%s\n", "end");
}
#endif//COUNT_INTERACTIONS
