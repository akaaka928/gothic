/*************************************************************************\
 *                                                                       *
                  last updated on 2017/01/17(Tue) 20:10:24
 *                                                                       *
 *    Input/Output Code of N-body calculation                            *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
/* #define USE_SZIP_COMPRESSION */
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
//-------------------------------------------------------------------------
#ifdef  MONITOR_ENERGY_ERROR
#include <math.h>
#endif//MONITOR_ENERGY_ERROR
//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
/* The maximum number of elements in a chunk is 2^32-1 which is equal to 4,294,967,295 */
/* The maximum size for any chunk is 4GB */
#define MAXIMUM_CHUNK_SIZE      ((hsize_t)1 << 31)
#define MAXIMUM_CHUNK_SIZE_4BIT ((hsize_t)1 << 30)
#define MAXIMUM_CHUNK_SIZE_8BIT ((hsize_t)1 << 29)
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
#include "macro.h"
#include "name.h"
#include "mpilib.h"
#include "constants.h"
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/tune.h"
#           if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#include "../misc/brent.h"
#        endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
//-------------------------------------------------------------------------
#include "io.h"
//-------------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
#       include <unistd.h>
#endif//EXEC_BENCHMARK
//-------------------------------------------------------------------------
#ifdef  HUNT_WALK_PARAMETER
#       include "../tree/walk_dev.h"
#ifdef  BRUTE_FORCE_LOCALIZATION
#ifdef  LOCALIZE_I_PARTICLES
#       include "../tree/neighbor_dev.h"
#endif//LOCALIZE_I_PARTICLES
#       include "../tree/shrink_dev.h"
#else///BRUTE_FORCE_LOCALIZATION
#       include "../tree/stat_dev.h"
#endif//BRUTE_FORCE_LOCALIZATION
#endif//HUNT_WALK_PARAMETER
//-------------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
#       include "../sort/peano_dev.h"
#       include "../tree/make_dev.h"
#endif//HUNT_MAKE_PARAMETER
//-------------------------------------------------------------------------
#ifdef  HUNT_NODE_PARAMETER
#       include "../tree/make_dev.h"
#endif//HUNT_NODE_PARAMETER
//-------------------------------------------------------------------------
#ifdef  HUNT_TIME_PARAMETER
#       include "../time/adv_dev.h"
#endif//HUNT_TIME_PARAMETER
//-------------------------------------------------------------------------
#   if  defined(HUNT_FIND_PARAMETER) && defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#       include "../tree/neighbor_dev.h"
#endif//defined(HUNT_FIND_PARAMETER) && defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
extern const real newton;
extern const double      length2astro;extern const char      length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double        time2astro;extern const char        time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double        mass2astro;extern const char        mass_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double     density2astro;extern const char     density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double col_density2astro;extern const char col_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double      energy2astro;extern const char      energy_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double     senergy2astro;extern const char     senergy_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double    velocity2astro;extern const char    velocity_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double       accel2astro;extern const char       accel_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void createMPIcfg_dataio(MPIcfg_dataio *cfg, MPIinfo mpi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  cfg->comm = mpi.comm;
  cfg->info = mpi.info;
  cfg->size = mpi.size;
  cfg->rank = mpi.rank;
  /* commitMPIbyte(&(cfg->body), (int)sizeof(nbody_particle)); */
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void updateMPIcfg_dataio(MPIcfg_dataio *cfg, int num)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  static ulong unit, psum;
  //-----------------------------------------------------------------------
  unit = (ulong)num;
  chkMPIerr(MPI_Scan(&unit, &psum, 1, MPI_UNSIGNED_LONG, MPI_SUM, cfg->comm));
  cfg->head = (MPI_Offset)(psum - unit);
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void writeConfigFileFormat(char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  fprintf(stderr, "Name of the configuration file should be \"%s\" (argv[1])\n", file);
  fprintf(stderr, "\tline 0: index of last updated tentative file series\n");
  __KILL__(stderr, "%s\n", "ERROR: format of the configuration file wrong.");
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void updateConfigFile(int last, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.%s", DOCUMENTFOLDER, file, CONFIG);
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  //-----------------------------------------------
  fprintf(fp, "%d\n", last);
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void readConfigFile(int *last, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.%s", DOCUMENTFOLDER, file, CONFIG);
  fp = fopen(filename, "r");
  if( fp == NULL ){
    writeConfigFileFormat(filename);
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  int checker = 1;
  //-----------------------------------------------------------------------
  checker &= (1 == fscanf(fp, "%d", last));
  //-----------------------------------------------------------------------
  fclose(fp);
  if( !checker )    writeConfigFileFormat(filename);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void readConfigFileParallel(int *last, char file[], MPIinfo mpi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* the root process read the specified file */
  if( mpi.rank == 0 ){
    //---------------------------------------------------------------------
    char filename[128];
    FILE *fp;
    //---------------------------------------------------------------------
    sprintf(filename, "%s/%s.%s", DOCUMENTFOLDER, file, CONFIG);
    fp = fopen(filename, "r");
    if( fp == NULL ){
      writeConfigFileFormat(filename);
      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
    }
    int checker = 1;
    //---------------------------------------------------------------------
    checker &= (1 == fscanf(fp, "%d", last));
    //---------------------------------------------------------------------
    fclose(fp);
    if( !checker )
      writeConfigFileFormat(filename);
    //---------------------------------------------------------------------
  }/* if( mpi.rank == 0 ){ */
  //-----------------------------------------------------------------------
  /* broadcast the read parameters from the root process */
  //-----------------------------------------------------------------------
  chkMPIerr(MPI_Bcast(last, 1, MPI_INT, 0, mpi.comm));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void readSettings(int *unit, ulong *Ntot, real *eps, real *eta, double *ft, double *SnapshotInterval, double *SaveInterval, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s", DATAFOLDER, file, CONFIG);
  fp = fopen(filename, "rb");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  bool success = true;
  size_t tmp;
  tmp = 1;  if( tmp != fread(            Ntot, sizeof( ulong), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(             eps, sizeof(  real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(             eta, sizeof(  real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(              ft, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(SnapshotInterval, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(    SaveInterval, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(            unit, sizeof(   int), tmp, fp) )    success = false;
  if( success != true ){
    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------
  setPhysicalConstantsAndUnitSystem(*unit, 0);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void readSettingsParallel(int *unit, ulong *Ntot, real *eps, real *eta, double *ft, double *SnapshotInterval, double *SaveInterval, char file[], MPIinfo mpi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* the root process read the specified file */
  if( mpi.rank == 0 ){
    //---------------------------------------------------------------------
    char filename[128];
    FILE *fp;
    sprintf(filename, "%s/%s.%s", DATAFOLDER, file, CONFIG);
    fp = fopen(filename, "rb");
    if( fp == NULL ){
      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
    }
    //---------------------------------------------------------------------
    bool success = true;
    size_t tmp;
    tmp = 1;  if( tmp != fread(            Ntot, sizeof( ulong), tmp, fp) )    success = false;
    tmp = 1;  if( tmp != fread(             eps, sizeof(  real), tmp, fp) )    success = false;
    tmp = 1;  if( tmp != fread(             eta, sizeof(  real), tmp, fp) )    success = false;
    tmp = 1;  if( tmp != fread(              ft, sizeof(double), tmp, fp) )    success = false;
    tmp = 1;  if( tmp != fread(SnapshotInterval, sizeof(double), tmp, fp) )    success = false;
    tmp = 1;  if( tmp != fread(    SaveInterval, sizeof(double), tmp, fp) )    success = false;
    tmp = 1;  if( tmp != fread(            unit, sizeof(   int), tmp, fp) )    success = false;
    if( success != true ){
      __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
    }
    //---------------------------------------------------------------------
    fclose(fp);
    //---------------------------------------------------------------------
  }/* if( mpi.rank == 0 ){ */
  //-----------------------------------------------------------------------
  /* broadcast the read parameters from the root process */
  chkMPIerr(MPI_Bcast(            Ntot, 1, MPI_UNSIGNED_LONG, 0, mpi.comm));
  chkMPIerr(MPI_Bcast(             eps, 1, MPI_REALDAT,       0, mpi.comm));
  chkMPIerr(MPI_Bcast(             eta, 1, MPI_REALDAT,       0, mpi.comm));
  chkMPIerr(MPI_Bcast(              ft, 1, MPI_DOUBLE,        0, mpi.comm));
  chkMPIerr(MPI_Bcast(SnapshotInterval, 1, MPI_DOUBLE,        0, mpi.comm));
  chkMPIerr(MPI_Bcast(    SaveInterval, 1, MPI_DOUBLE,        0, mpi.comm));
  chkMPIerr(MPI_Bcast(            unit, 1, MPI_INT,           0, mpi.comm));
  //-----------------------------------------------------------------------
  setPhysicalConstantsAndUnitSystem(*unit, 0);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void writeSettings(int unit, ulong Ntot, real eps, real eta, double ft, double SnapshotInterval, double SaveInterval, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s", DATAFOLDER, file, CONFIG);
  fp = fopen(filename, "wb");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  bool success = true;
  size_t tmp;
  tmp = 1;  if( tmp != fwrite(&            Ntot, sizeof( ulong), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&             eps, sizeof(  real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&             eta, sizeof(  real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&              ft, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&SnapshotInterval, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&    SaveInterval, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&            unit, sizeof(   int), tmp, fp) )    success = false;
  if( success != true ){
    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
//-------------------------------------------------------------------------
void createHDF5DataType(hdf5struct *type)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit a data type to write unit name as strings */
  //-----------------------------------------------------------------------
  type->str4unit = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(type->str4unit, CONSTANTS_H_CHAR_WORDS));/* memo: sizeof(char) is unity */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit a data type to switch float and double */
  //-----------------------------------------------------------------------
#ifdef  DOUBLE_PRECISION
  type->real = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
  type->real = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit a data type of rebuildTree */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit a data type of measuredTime */
  //-----------------------------------------------------------------------
  type->measuredTime = H5Tcreate(H5T_COMPOUND, sizeof(measuredTime));
  chkHDF5err(H5Tinsert(type->measuredTime, "walkTree[0]", HOFFSET(measuredTime, walkTree[0]), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->measuredTime, "walkTree[1]", HOFFSET(measuredTime, walkTree[1]), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->measuredTime, "makeTree"   , HOFFSET(measuredTime, makeTree   ), H5T_NATIVE_DOUBLE));
#ifdef  WALK_TREE_TOTAL_SUM_MODEL
  chkHDF5err(H5Tinsert(type->measuredTime,   "incSum"   , HOFFSET(measuredTime,  incSum    ), H5T_NATIVE_DOUBLE));
#endif//WALK_TREE_TOTAL_SUM_MODEL
#ifndef SERIALIZED_EXECUTION
  chkHDF5err(H5Tinsert(type->measuredTime,  "sum_excg"   , HOFFSET(measuredTime, sum_excg   ), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->measuredTime,  "sum_rebuild", HOFFSET(measuredTime, sum_rebuild), H5T_NATIVE_DOUBLE));
#ifdef  MONITOR_LETGEN_TIME
  chkHDF5err(H5Tinsert(type->measuredTime,  "excg"       , HOFFSET(measuredTime, excg       ), H5T_NATIVE_DOUBLE));
#endif//MONITOR_LETGEN_TIME
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
  //-----------------------------------------------------------------------
  /* commit a data type of statVal */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit a data type of guessTime */
  //-----------------------------------------------------------------------
  type->guessTime = H5Tcreate(H5T_COMPOUND, sizeof(guessTime));
  chkHDF5err(H5Tinsert(type->guessTime,  "slope", HOFFSET(guessTime,  slope), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->guessTime,  "icept", HOFFSET(guessTime,  icept), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->guessTime, "rchisq", HOFFSET(guessTime, rchisq), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->guessTime,   "time", HOFFSET(guessTime,   time), H5T_NATIVE_DOUBLE));
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  chkHDF5err(H5Tinsert(type->guessTime, "second", HOFFSET(guessTime, second), H5T_NATIVE_DOUBLE));
#endif//USE_PARABOLIC_GROWTH_MODEL
  //-----------------------------------------------------------------------
#endif//WALK_TREE_COMBINED_MODEL
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  USE_BRENT_METHOD
  //-----------------------------------------------------------------------
  /* commit a data type of brentFunc */
  //-----------------------------------------------------------------------
  type->brentFunc = H5Tcreate(H5T_COMPOUND, sizeof(brentFunc));
  chkHDF5err(H5Tinsert(type->brentFunc, "pos", HOFFSET(brentFunc, pos), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->brentFunc, "val", HOFFSET(brentFunc, val), H5T_NATIVE_DOUBLE));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit a data type of brentStatus */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit a data type of brentMemory */
  //-----------------------------------------------------------------------
  type->brentMemory = H5Tcreate(H5T_COMPOUND, sizeof(brentMemory));
  chkHDF5err(H5Tinsert(type->brentMemory, "previous", HOFFSET(brentMemory, previous), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->brentMemory,   "totNum", HOFFSET(brentMemory,   totNum), H5T_NATIVE_INT));
  chkHDF5err(H5Tinsert(type->brentMemory, "degraded", HOFFSET(brentMemory, degraded), H5T_NATIVE_INT));
  chkHDF5err(H5Tinsert(type->brentMemory, "interval", HOFFSET(brentMemory, interval), H5T_NATIVE_INT));
  //-----------------------------------------------------------------------
#endif//USE_BRENT_METHOD
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void removeHDF5DataType(hdf5struct  type)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  USE_BRENT_METHOD
  chkHDF5err(H5Tclose(type.brentMemory));
  chkHDF5err(H5Tclose(type.brentStatus));
  chkHDF5err(H5Tclose(type.brentFunc));
#endif//USE_BRENT_METHOD
  //-----------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
  chkHDF5err(H5Tclose(type.guessTime));
  chkHDF5err(H5Tclose(type.statVal));
#endif//WALK_TREE_COMBINED_MODEL
  chkHDF5err(H5Tclose(type.measuredTime));
  chkHDF5err(H5Tclose(type.rebuildTree));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Tclose(type.str4unit));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  readTentativeData(double *time, double *dt, ulong *steps, int num, iparticle body, char file[], int  last, hdf5struct type
			, int *dropPrevTune
			, rebuildTree *rebuild, measuredTime *measured
#ifdef  WALK_TREE_COMBINED_MODEL
			, autoTuningParam *rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
			, brentStatus *status, brentMemory *memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
			, energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
			)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* open an existing file with read only option */
  //-----------------------------------------------------------------------
  char filename[128];
  sprintf(filename, "%s/%s.%s%d.h5", DATAFOLDER, file, TENTATIVE, last);
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  //-----------------------------------------------------------------------
  /* open an existing group */
  hid_t group = H5Gopen(target, "nbody", H5P_DEFAULT);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* read attribute data */
  //-----------------------------------------------------------------------
  hid_t attribute;
  //-----------------------------------------------------------------------
  /* read current time */
  attribute = H5Aopen(group, "time", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, time));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* read time step */
  attribute = H5Aopen(group, "dt", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, dt));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* read # of steps */
  attribute = H5Aopen(group, "steps", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, steps));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* read # of N-body particles */
  ulong num_ulong = 0;
  attribute = H5Aopen(group, "number", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, &num_ulong));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* read inverse of total energy at the initial condition */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Aopen(group, "E0inv", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &relEneErr->E0inv));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  /* read maximum value of the relative error of the total energy */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Aopen(group, "errMax", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &relEneErr->errMax));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  /* read flag about DOUBLE_PRECISION */
  int useDP;
  attribute = H5Aopen(group, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* read flag about BLOCK_TIME_STEP */
  int blockTimeStep;
  attribute = H5Aopen(group, "blockTimeStep", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &blockTimeStep));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* simple error check */
  //-----------------------------------------------------------------------
  if( num_ulong != (ulong)num ){    __KILL__(stderr, "ERROR: number of N-body particles in the file (%zu) differs with that in the code (%d)\n", num_ulong, num);  }
  //-----------------------------------------------------------------------
#ifdef  DOUBLE_PRECISION
  if( useDP != 1 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);
  }/* if( useDP != 1 ){ */
#else///DOUBLE_PRECISION
  if( useDP != 0 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);
  }/* if( useDP != 0 ){ */
#endif//DOUBLE_PRECISION
  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
  if( blockTimeStep != 1 ){
    __KILL__(stderr, "ERROR: blockTimeStep (%d) differs with that in the code (%d)\n", blockTimeStep, true);
  }/* if( blockTimeStep != 1 ){ */
#else///BLOCK_TIME_STEP
  if( blockTimeStep != 0 ){
    __KILL__(stderr, "ERROR: blockTimeStep (%d) differs with that in the code (%d)\n", blockTimeStep, false);
  }/* if( blockTimeStep != 0 ){ */
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* read particle data */
  //-----------------------------------------------------------------------
  hid_t dataset;
  /* read particle position */
  dataset = H5Dopen(group, "position", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.pos));
  chkHDF5err(H5Dclose(dataset));
  /* read particle acceleration */
  dataset = H5Dopen(group, "acceleration", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.acc));
  chkHDF5err(H5Dclose(dataset));
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
  //-----------------------------------------------------------------------
  /* end access to the dataset and release the corresponding resource */
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  if( (*steps != 0) && (*dropPrevTune == 0) ){
    //---------------------------------------------------------------------
    group = H5Gopen(target, "parameters in auto-tuning", H5P_DEFAULT);
    //---------------------------------------------------------------------
    /* read attributes */
    //---------------------------------------------------------------------
    /* read flag about FORCE_ADJUSTING_PARTICLE_TIME_STEPS */
    int forceAdjust;
    attribute = H5Aopen(group, "FORCE_ADJUSTING_PARTICLE_TIME_STEPS", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &forceAdjust));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* read flag about WALK_TREE_COMBINED_MODEL */
    int combined;
    attribute = H5Aopen(group, "WALK_TREE_COMBINED_MODEL", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &combined));
    chkHDF5err(H5Aclose(attribute));
    /* read flag about WALK_TREE_TOTAL_SUM_MODEL */
    int totSum;
    attribute = H5Aopen(group, "WALK_TREE_TOTAL_SUM_MODEL", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &totSum));
    chkHDF5err(H5Aclose(attribute));
    /* read flag about USE_PARABOLIC_GROWTH_MODEL */
    int parabolic;
    attribute = H5Aopen(group, "USE_PARABOLIC_GROWTH_MODEL", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &parabolic));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* read flag about LOCALIZE_I_PARTICLES */
    int localize;
    attribute = H5Aopen(group, "LOCALIZE_I_PARTICLES", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &localize));
    chkHDF5err(H5Aclose(attribute));
    /* read flag about USE_BRENT_METHOD */
    int useBrent;
    attribute = H5Aopen(group, "USE_BRENT_METHOD", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useBrent));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* read parameters */
    //---------------------------------------------------------------------
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
#ifdef  WALK_TREE_TOTAL_SUM_MODEL
    if( totSum == 1 )
#else///WALK_TREE_TOTAL_SUM_MODEL
      if( totSum == 0 )
#endif//WALK_TREE_TOTAL_SUM_MODEL
	{
	  dataset = H5Dopen(group, "measured time", H5P_DEFAULT);
	  chkHDF5err(H5Dread(dataset, type.measuredTime, H5S_ALL, H5S_ALL, H5P_DEFAULT, measured));
	  chkHDF5err(H5Dclose(dataset));
	}
    //---------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
    if( combined == 1 ){
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
    }/* if( combined == 1 ){ */
#endif//WALK_TREE_COMBINED_MODEL
    //---------------------------------------------------------------------
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
    if( (localize == 1) && (useBrent == 1) ){
      /* read brentStatus */
      dataset = H5Dopen(group, "Brent status", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, type.brentStatus, H5S_ALL, H5S_ALL, H5P_DEFAULT, status));
      chkHDF5err(H5Dclose(dataset));
      /* read brentMemory */
      dataset = H5Dopen(group, "Brent memory", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, type.brentMemory, H5S_ALL, H5S_ALL, H5P_DEFAULT, memory));
      chkHDF5err(H5Dclose(dataset));
    }/* if( (localize == 1) && (useBrent == 1) ){ */
    else
      *dropPrevTune = 1;
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    chkHDF5err(H5Gclose(group));
    //---------------------------------------------------------------------
  }/* if( (*steps != 0) && (*dropPrevTune == 0) ){ */
  else
    *dropPrevTune = 1;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* close the file */
  //-----------------------------------------------------------------------
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  readTentativeDataParallel(double *time, double *dt, ulong *steps, int *num, iparticle body, char file[], int  last, MPIcfg_dataio *mpi, ulong Ntot, hdf5struct type
				, int *dropPrevTune, rebuildTree *rebuild, measuredTime *measured
#ifdef  WALK_TREE_COMBINED_MODEL
				, autoTuningParam *rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
				, brentStatus *status, brentMemory *memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
				, energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
				)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* open an existing file with read only option */
  //-----------------------------------------------------------------------
  char filename[128];
  sprintf(filename, "%s/%s.%s%d.h5", DATAFOLDER, file, TENTATIVE, last);
  hid_t f_property = H5Pcreate(H5P_FILE_ACCESS);
  chkHDF5err(H5Pset_fapl_mpio(f_property, mpi->comm, mpi->info));
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, f_property);
  chkHDF5err(H5Pclose(f_property));
  //-----------------------------------------------------------------------
  /* open an existing group */
  hid_t group = H5Gopen(target, "nbody", H5P_DEFAULT);
  //-----------------------------------------------------------------------
  /* create property list for collective dataset read */
  hid_t r_property = H5Pcreate(H5P_DATASET_XFER);
  chkHDF5err(H5Pset_dxpl_mpio(r_property, H5FD_MPIO_COLLECTIVE));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* read attribute data */
  //-----------------------------------------------------------------------
  ulong Nread = 0;
  int useDP;
  int blockTimeStep;
  //-----------------------------------------------------------------------
  if( mpi->rank == 0 ){
    //---------------------------------------------------------------------
    hid_t attribute;
    //---------------------------------------------------------------------
    /* read current time */
    attribute = H5Aopen(group, "time", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, time));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* read time step */
    attribute = H5Aopen(group, "dt", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, dt));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* read # of steps */
    attribute = H5Aopen(group, "steps", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, steps));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* read # of N-body particles */
    attribute = H5Aopen(group, "number", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, &Nread));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* read inverse of total energy at the initial condition */
#ifdef  MONITOR_ENERGY_ERROR
    attribute = H5Aopen(group, "E0inv", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &relEneErr->E0inv));
    chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
    //---------------------------------------------------------------------
    /* read maximum value of the relative error of the total energy */
#ifdef  MONITOR_ENERGY_ERROR
    attribute = H5Aopen(group, "errMax", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &relEneErr->errMax));
    chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
    //---------------------------------------------------------------------
    /* read flag about DOUBLE_PRECISION */
    attribute = H5Aopen(group, "useDP", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* read flag about BLOCK_TIME_STEP */
    attribute = H5Aopen(group, "blockTimeStep", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &blockTimeStep));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
  }/* if( mpi->rank == 0 ){ */
  //-----------------------------------------------------------------------
  chkMPIerr(MPI_Bcast(  time, 1, MPI_DOUBLE       , 0, mpi->comm));
  chkMPIerr(MPI_Bcast(    dt, 1, MPI_DOUBLE       , 0, mpi->comm));
  chkMPIerr(MPI_Bcast( steps, 1, MPI_UNSIGNED_LONG, 0, mpi->comm));
  chkMPIerr(MPI_Bcast(&Nread, 1, MPI_UNSIGNED_LONG, 0, mpi->comm));
#ifdef  MONITOR_ENERGY_ERROR
  chkMPIerr(MPI_Bcast(&relEneErr->E0inv , 1, MPI_DOUBLE, 0, mpi->comm));
  chkMPIerr(MPI_Bcast(&relEneErr->errMax, 1, MPI_DOUBLE, 0, mpi->comm));
#endif//MONITOR_ENERGY_ERROR
  chkMPIerr(MPI_Bcast(&useDP, 1, MPI_INT, 0, mpi->comm));
  chkMPIerr(MPI_Bcast(&blockTimeStep, 1, MPI_INT, 0, mpi->comm));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* simple error check */
  //-----------------------------------------------------------------------
  if( Nread != Ntot ){
    __KILL__(stderr, "ERROR: number of N-body particles in the file (%zu) differs with that in the code (%zu)\n", Nread, Ntot);
  }/* if( Nread != Ntot ){ */
  //-----------------------------------------------------------------------
#ifdef  DOUBLE_PRECISION
  if( useDP != 1 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);
  }/* if( useDP != 1 ){ */
#else///DOUBLE_PRECISION
  if( useDP != 0 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);
  }/* if( useDP != 0 ){ */
#endif//DOUBLE_PRECISION
  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
  if( blockTimeStep != 1 ){
    __KILL__(stderr, "ERROR: blockTimeStep (%d) differs with that in the code (%d)\n", blockTimeStep, true);
  }/* if( blockTimeStep != 1 ){ */
#else///BLOCK_TIME_STEP
  if( blockTimeStep != 0 ){
    __KILL__(stderr, "ERROR: blockTimeStep (%d) differs with that in the code (%d)\n", blockTimeStep, false);
  }/* if( blockTimeStep != 0 ){ */
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  if( (*steps != 0) && (*dropPrevTune == 0) ){
    //---------------------------------------------------------------------
    chkHDF5err(H5Gclose(group));
    group = H5Gopen(target, "parameters in auto-tuning", H5P_DEFAULT);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* read attributes */
    //---------------------------------------------------------------------
    int forceAdjust;
    int monitor, combined, totSum, parabolic;
    int localize, useBrent;
    //---------------------------------------------------------------------
    if( mpi->rank == 0 ){
      //-------------------------------------------------------------------
      hid_t attribute;
      //-------------------------------------------------------------------
      /* read flag about FORCE_ADJUSTING_PARTICLE_TIME_STEPS */
      attribute = H5Aopen(group, "FORCE_ADJUSTING_PARTICLE_TIME_STEPS", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &forceAdjust));
      chkHDF5err(H5Aclose(attribute));
      //-------------------------------------------------------------------
      /* read flag about MONITOR_LETGEN_TIME */
      attribute = H5Aopen(group, "MONITOR_LETGEN_TIME", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &monitor));
      chkHDF5err(H5Aclose(attribute));
      /* read flag about WALK_TREE_COMBINED_MODEL */
      attribute = H5Aopen(group, "WALK_TREE_COMBINED_MODEL", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &combined));
      chkHDF5err(H5Aclose(attribute));
      /* read flag about WALK_TREE_TOTAL_SUM_MODEL */
      attribute = H5Aopen(group, "WALK_TREE_TOTAL_SUM_MODEL", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &totSum));
      chkHDF5err(H5Aclose(attribute));
      /* read flag about USE_PARABOLIC_GROWTH_MODEL */
      attribute = H5Aopen(group, "USE_PARABOLIC_GROWTH_MODEL", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &parabolic));
      chkHDF5err(H5Aclose(attribute));
      //-------------------------------------------------------------------
      /* read flag about LOCALIZE_I_PARTICLES */
      attribute = H5Aopen(group, "LOCALIZE_I_PARTICLES", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &localize));
      chkHDF5err(H5Aclose(attribute));
      /* read flag about USE_BRENT_METHOD */
      attribute = H5Aopen(group, "USE_BRENT_METHOD", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useBrent));
      chkHDF5err(H5Aclose(attribute));
      //-------------------------------------------------------------------
    }/* if( mpi->rank == 0 ){ */
    //---------------------------------------------------------------------
    chkMPIerr(MPI_Bcast(&forceAdjust, 1, MPI_INT, 0, mpi->comm));
    chkMPIerr(MPI_Bcast(&    monitor, 1, MPI_INT, 0, mpi->comm));
    chkMPIerr(MPI_Bcast(&   combined, 1, MPI_INT, 0, mpi->comm));
    chkMPIerr(MPI_Bcast(&     totSum, 1, MPI_INT, 0, mpi->comm));
    chkMPIerr(MPI_Bcast(&  parabolic, 1, MPI_INT, 0, mpi->comm));
    chkMPIerr(MPI_Bcast(&   localize, 1, MPI_INT, 0, mpi->comm));
    chkMPIerr(MPI_Bcast(&   useBrent, 1, MPI_INT, 0, mpi->comm));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* read parameters */
    //---------------------------------------------------------------------
    hsize_t dims_loc = 1;
    hid_t locSpace = H5Screate_simple(1, &dims_loc, NULL);
    /* configuration about domain decomposition */
    hsize_t  count = 1;
    hsize_t stride = 1;
    hsize_t  block = dims_loc;
    hsize_t offset = mpi->rank;
    hid_t dataset, hyperslab;
    //---------------------------------------------------------------------
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
    /* WALK_TREE_TOTAL_SUM_MODEL, MONITOR_LETGEN_TIME */
#ifdef  WALK_TREE_TOTAL_SUM_MODEL
    if( totSum == 1 )
#else///WALK_TREE_TOTAL_SUM_MODEL
      if( totSum == 0 )
#endif//WALK_TREE_TOTAL_SUM_MODEL
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
    //---------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
    if( combined == 1 ){
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
    }/* if( combined == 1 ){ */
#endif//WALK_TREE_COMBINED_MODEL
    //---------------------------------------------------------------------
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
    if( (localize == 1) && (useBrent == 1) ){
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
    }/* if( (localize == 1) && (useBrent == 1) ){ */
    else
      *dropPrevTune = 1;
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
    //---------------------------------------------------------------------
    /* read # of particles in each process */
    dataset = H5Dopen(group, "num", H5P_DEFAULT);
    hyperslab = H5Dget_space(dataset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dread(dataset, H5T_NATIVE_INT, locSpace, hyperslab, r_property, num));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dataset));
    //---------------------------------------------------------------------
    chkHDF5err(H5Sclose(locSpace));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    chkHDF5err(H5Gclose(group));
    group = H5Gopen(target, "nbody", H5P_DEFAULT);
    //---------------------------------------------------------------------
  }/* if( (*steps != 0) && (*dropPrevTune == 0) ){ */
  else
    *dropPrevTune = 1;
  //-----------------------------------------------------------------------
  /* reset MPI_Offset */
  updateMPIcfg_dataio(mpi, *num);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* read particle data */
  //-----------------------------------------------------------------------
  /* dataset for real4 arrays */
  //-----------------------------------------------------------------------
  /* dataspace */
  hsize_t dims_loc = 4 * (*num);
  hid_t locSpace = H5Screate_simple(1, &dims_loc, NULL);
  //-----------------------------------------------------------------------
  /* configuration about domain decomposition */
  hsize_t  count = 1;
  hsize_t stride = 1;
  hsize_t  block = dims_loc;
  hsize_t offset = mpi->head * 4;
  //-----------------------------------------------------------------------
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
#ifdef  BLOCK_TIME_STEP
  dataset = H5Dopen(group, "velocity", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, type.real, locSpace, hyperslab, r_property, body.vel));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(locSpace));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
  /* dataset for double2 array */
  //-----------------------------------------------------------------------
  /* dataspace */
  dims_loc = 2 * (*num);
  locSpace = H5Screate_simple(1, &dims_loc, NULL);
  //-----------------------------------------------------------------------
  /* configuration about domain decomposition */
  count  = 1;
  stride = 1;
  block  = dims_loc;
  offset = mpi->head * 2;
  //-----------------------------------------------------------------------
  /* read particle time */
  dataset = H5Dopen(group, "time", H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, locSpace, hyperslab, r_property, body.time));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(locSpace));
  //-----------------------------------------------------------------------
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* dataset for num array */
  //-----------------------------------------------------------------------
  /* dataspace */
  dims_loc = *num;
  locSpace = H5Screate_simple(1, &dims_loc, NULL);
  //-----------------------------------------------------------------------
  /* configuration about domain decomposition */
  count  = 1;
  stride = 1;
  block  = dims_loc;
  offset = mpi->head;
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(locSpace));
  //-----------------------------------------------------------------------
  /* finish collective dataset read */
  chkHDF5err(H5Pclose(r_property));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* close the file */
  //-----------------------------------------------------------------------
  /* end access to the dataset and release the corresponding resource */
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
#if 0
  sprintf(filename, "%s/%s.%s%d_%d.txt", DATAFOLDER, file, "rank", mpi->rank, mpi->size);
  FILE *fp;
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  for(int ii = 0; ii < *num; ii++)
    fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%zu\n",
	    body.pos[ii].m,
	    body.pos[ii].x, body.pos[ii].y, body.pos[ii].z,
	    body.vel[ii].x, body.vel[ii].y, body.vel[ii].z,
	    body.idx[ii]);
  fclose(fp);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeTentativeData(double  time, double  dt, ulong  steps, ulong num, iparticle body, char file[], int *last, hdf5struct type
			, rebuildTree rebuild, measuredTime measured
#ifdef  WALK_TREE_COMBINED_MODEL
			, autoTuningParam rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
			, brentStatus status, brentMemory memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
			, energyError relEneErr
#endif//MONITOR_ENERGY_ERROR
			)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  //-----------------------------------------------------------------------
  char filename[128];
  sprintf(filename, "%s/%s.%s%d.h5", DATAFOLDER, file, TENTATIVE, (*last) ^ 1);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //-----------------------------------------------------------------------
  hid_t group = H5Gcreate(target, "nbody", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* preparation for data compression */
  //-----------------------------------------------------------------------
  hid_t dataset, dataspace, property;
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  uint szip_pixels_per_block = 8;
  hsize_t cdims = 128 * szip_pixels_per_block;
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* num * real4 arrays */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* write particle position */
  dataset = H5Dcreate(group, "position", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.pos));
  chkHDF5err(H5Dclose(dataset));
  /* write particle acceleration */
  dataset = H5Dcreate(group, "acceleration", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.acc));
  chkHDF5err(H5Dclose(dataset));
#ifdef  BLOCK_TIME_STEP
  /* write particle velocity */
  dataset = H5Dcreate(group, "velocity", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.vel));
  chkHDF5err(H5Dclose(dataset));
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
  /* num * double2 array */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* write particle time */
  dataset = H5Dcreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.time));
  chkHDF5err(H5Dclose(dataset));
  //-----------------------------------------------------------------------
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* num array */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* write particle index */
  dataset = H5Dcreate(group, "index", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, body.idx));
  chkHDF5err(H5Dclose(dataset));
  //-----------------------------------------------------------------------
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* write attribute data */
  //-----------------------------------------------------------------------
  /* create the data space for the attribute */
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  //-----------------------------------------------------------------------
  /* write current time */
  attribute = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write time step */
  attribute = H5Acreate(group, "dt", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &dt));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write # of steps */
  attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write # of N-body particles */
  ulong num_ulong = (ulong)num;
  attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &num_ulong));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write inverse of total energy at the initial condition */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Acreate(group, "E0inv", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &relEneErr.E0inv));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  /* write maximum value of the relative error of the total energy */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Acreate(group, "errMax", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &relEneErr.errMax));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  /* write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
  attribute = H5Acreate(group, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write flag about BLOCK_TIME_STEP */
#ifdef  BLOCK_TIME_STEP
  const int blockTimeStep = 1;
#else///BLOCK_TIME_STEP
  const int blockTimeStep = 0;
#endif//BLOCK_TIME_STEP
  attribute = H5Acreate(group, "blockTimeStep", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &blockTimeStep));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output parameters for auto-tuning (when steps > 0) */
  //-----------------------------------------------------------------------
  if( steps != 0 ){
    //---------------------------------------------------------------------
    group = H5Gcreate(target, "parameters in auto-tuning", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write parameters */
    //---------------------------------------------------------------------
    dims = 1;
    dataspace = H5Screate_simple(1, &dims, NULL);
    //---------------------------------------------------------------------
    /* output rebuildTree */
    dataset = H5Dcreate(group, "rebuild tree", type.rebuildTree, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.rebuildTree, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rebuild));
    chkHDF5err(H5Dclose(dataset));
    /* output measuredTime */
    dataset = H5Dcreate(group, "measured time", type.measuredTime, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.measuredTime, H5S_ALL, H5S_ALL, H5P_DEFAULT, &measured));
    chkHDF5err(H5Dclose(dataset));
    //---------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
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
#endif//WALK_TREE_COMBINED_MODEL
    //---------------------------------------------------------------------
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
    /* output brentStatus */
    dataset = H5Dcreate(group, "Brent status", type.brentStatus, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.brentStatus, H5S_ALL, H5S_ALL, H5P_DEFAULT, &status));
    chkHDF5err(H5Dclose(dataset));
    /* output brentMemory */
    dataset = H5Dcreate(group, "Brent memory", type.brentMemory, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.brentMemory, H5S_ALL, H5S_ALL, H5P_DEFAULT, &memory));
    chkHDF5err(H5Dclose(dataset));
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
    //---------------------------------------------------------------------
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write attributes */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    int flag;
    //---------------------------------------------------------------------
    /* write flag about FORCE_ADJUSTING_PARTICLE_TIME_STEPS */
#ifdef  FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    flag = 1;
#else///FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    flag = 0;
#endif//FORCE_ADJUSTING_PARTICLE_TIME_STEPS
    attribute = H5Acreate(group, "FORCE_ADJUSTING_PARTICLE_TIME_STEPS", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about WALK_TREE_COMBINED_MODEL */
#ifdef  WALK_TREE_COMBINED_MODEL
    flag = 1;
#else///WALK_TREE_COMBINED_MODEL
    flag = 0;
#endif//WALK_TREE_COMBINED_MODEL
    attribute = H5Acreate(group, "WALK_TREE_COMBINED_MODEL", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about WALK_TREE_TOTAL_SUM_MODEL */
#ifdef  WALK_TREE_TOTAL_SUM_MODEL
    flag = 1;
#else///WALK_TREE_TOTAL_SUM_MODEL
    flag = 0;
#endif//WALK_TREE_TOTAL_SUM_MODEL
    attribute = H5Acreate(group, "WALK_TREE_TOTAL_SUM_MODEL", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
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
    /* write flag about LOCALIZE_I_PARTICLES */
#ifdef  LOCALIZE_I_PARTICLES
    flag = 1;
#else///LOCALIZE_I_PARTICLES
    flag = 0;
#endif//LOCALIZE_I_PARTICLES
    attribute = H5Acreate(group, "LOCALIZE_I_PARTICLES", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about USE_BRENT_METHOD */
#ifdef  USE_BRENT_METHOD
    flag = 1;
#else///USE_BRENT_METHOD
    flag = 0;
#endif//USE_BRENT_METHOD
    attribute = H5Acreate(group, "USE_BRENT_METHOD", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about TEST_BRENT_METHOD */
#ifdef  TEST_BRENT_METHOD
    flag = 1;
#else///TEST_BRENT_METHOD
    flag = 0;
#endif//TEST_BRENT_METHOD
    attribute = H5Acreate(group, "TEST_BRENT_METHOD", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
    //---------------------------------------------------------------------
  }/* if( steps != 0 ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* close the file */
  //-----------------------------------------------------------------------
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
  *last ^= 1;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeTentativeDataParallel(double  time, double  dt, ulong  steps, int num, iparticle body, char file[], int *last, MPIcfg_dataio *mpi, ulong Ntot, hdf5struct type
				, rebuildTree rebuild, measuredTime measured
#ifdef  WALK_TREE_COMBINED_MODEL
				, autoTuningParam rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
				, brentStatus status, brentMemory memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
				, energyError relEneErr
#endif//MONITOR_ENERGY_ERROR
				)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  //-----------------------------------------------------------------------
  char filename[128];
  hid_t f_property = H5Pcreate(H5P_FILE_ACCESS);
  chkHDF5err(H5Pset_fapl_mpio(f_property, mpi->comm, mpi->info));
  sprintf(filename, "%s/%s.%s%d.h5", DATAFOLDER, file, TENTATIVE, (*last) ^ 1);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, f_property);
  chkHDF5err(H5Pclose(f_property));
  //-----------------------------------------------------------------------
  hid_t group = H5Gcreate(target, "nbody", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //-----------------------------------------------------------------------
  /* create property list for collective dataset write */
  hid_t w_property = H5Pcreate(H5P_DATASET_XFER);
  chkHDF5err(H5Pset_dxpl_mpio(w_property, H5FD_MPIO_COLLECTIVE));
  //-----------------------------------------------------------------------
  /* configuration about domain decomposition */
  updateMPIcfg_dataio(mpi, num);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create (distributed) dataset for real4 arrays */
  //-----------------------------------------------------------------------
  /* create dataspace */
  hsize_t dims_ful = 4 * Ntot            ;  hid_t fulSpace = H5Screate_simple(1, &dims_ful, NULL);
  hsize_t dims_loc = 4 *  num            ;  hid_t locSpace = H5Screate_simple(1, &dims_loc, NULL);
  hsize_t dims_mem = 4 * Ntot / mpi->size;
  //-----------------------------------------------------------------------
  /* configuration about domain decomposition */
  hsize_t  count = 1;
  hsize_t stride = 1;
  hsize_t  block = dims_loc;
  hsize_t offset = mpi->head * 4;
  //-----------------------------------------------------------------------
  /* create chunked dataset */
  hid_t data_create = H5Pcreate(H5P_DATASET_CREATE);
  hsize_t dims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  if( dims_mem > dims_max )
    dims_mem = dims_max;
  chkHDF5err(H5Pset_chunk(data_create, 1, &dims_mem));
  hid_t dset, hyperslab;
  //-----------------------------------------------------------------------
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
#ifdef  BLOCK_TIME_STEP
  /* write particle velocity */
  dset = H5Dcreate(group, "velocity", type.real, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, type.real, locSpace, hyperslab, w_property, body.vel));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
  /* close/release resources */
  chkHDF5err(H5Pclose(data_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
  /* create (distributed) dataset for double2 array */
  //-----------------------------------------------------------------------
  /* create dataspace */
  dims_ful = 2 * Ntot            ;  fulSpace = H5Screate_simple(1, &dims_ful, NULL);
  dims_mem = 2 * Ntot / mpi->size;
  dims_loc = 2 *  num            ;  locSpace = H5Screate_simple(1, &dims_loc, NULL);
  //-----------------------------------------------------------------------
  /* configuration about domain decomposition */
  count  = 1;
  stride = 1;
  block  = dims_loc;
  offset = mpi->head * 2;
  //-----------------------------------------------------------------------
  /* create chunked dataset */
  data_create = H5Pcreate(H5P_DATASET_CREATE);
  dims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
  if( dims_mem > dims_max )
    dims_mem = dims_max;
  chkHDF5err(H5Pset_chunk(data_create, 1, &dims_mem));
  //-----------------------------------------------------------------------
  /* write particle time */
  dset = H5Dcreate(group, "time", H5T_NATIVE_DOUBLE, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  chkHDF5err(H5Dwrite(dset, H5T_NATIVE_DOUBLE, locSpace, hyperslab, w_property, body.time));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));
  //-----------------------------------------------------------------------
  /* close/release resources */
  chkHDF5err(H5Pclose(data_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));
  //-----------------------------------------------------------------------
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create (distributed) dataset for num arrays */
  //-----------------------------------------------------------------------
  /* create dataspace */
  dims_ful = Ntot            ;  fulSpace = H5Screate_simple(1, &dims_ful, NULL);
  dims_mem = Ntot / mpi->size;
  dims_loc =  num            ;  locSpace = H5Screate_simple(1, &dims_loc, NULL);
  //-----------------------------------------------------------------------
  /* configuration about domain decomposition */
  count  = 1;
  stride = 1;
  block  = dims_loc;
  offset = mpi->head;
  //-----------------------------------------------------------------------
  /* create chunked dataset */
  data_create = H5Pcreate(H5P_DATASET_CREATE);
  dims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
  if( dims_mem > dims_max )
    dims_mem = dims_max;
  chkHDF5err(H5Pset_chunk(data_create, 1, &dims_mem));
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* close/release resources */
  chkHDF5err(H5Pclose(data_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));
  //-----------------------------------------------------------------------
  /* /\* finish collective dataset write *\/ */
  /* chkHDF5err(H5Pclose(w_property)); */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* write attribute data */
  //-----------------------------------------------------------------------
  /* create the data space for the attribute */
  hsize_t attr_dims = 1;
  hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  //-----------------------------------------------------------------------
  /* write current time */
  attribute = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write time step */
  attribute = H5Acreate(group, "dt", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &dt));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write # of steps */
  attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write # of N-body particles */
  attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &Ntot));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write inverse of total energy at the initial condition */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Acreate(group, "E0inv", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &relEneErr.E0inv));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  /* write maximum value of the relative error of the total energy */
#ifdef  MONITOR_ENERGY_ERROR
  attribute = H5Acreate(group, "errMax", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &relEneErr.errMax));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  /* write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
  attribute = H5Acreate(group, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write flag about BLOCK_TIME_STEP */
#ifdef  BLOCK_TIME_STEP
  const int blockTimeStep = 1;
#else///BLOCK_TIME_STEP
  const int blockTimeStep = 0;
#endif//BLOCK_TIME_STEP
  attribute = H5Acreate(group, "blockTimeStep", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &blockTimeStep));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output parameters for auto-tuning (when steps > 0) */
  //-----------------------------------------------------------------------
  if( steps != 0 ){
    //---------------------------------------------------------------------
    group = H5Gcreate(target, "parameters in auto-tuning", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //---------------------------------------------------------------------
    /* /\* create property list for collective dataset write *\/ */
    /* w_property = H5Pcreate(H5P_DATASET_XFER); */
    /* chkHDF5err(H5Pset_dxpl_mpio(w_property, H5FD_MPIO_COLLECTIVE)); */
    //---------------------------------------------------------------------
    /* create dataspace */
    dims_ful = mpi->size;    fulSpace = H5Screate_simple(1, &dims_ful, NULL);
    dims_loc =         1;    locSpace = H5Screate_simple(1, &dims_loc, NULL);
    dims_mem =         1;
    //---------------------------------------------------------------------
    /* configuration about domain decomposition */
    count  = 1;
    stride = 1;
    block  = dims_loc;
    offset = mpi->rank;
    //---------------------------------------------------------------------
    /* create chunked dataset */
    data_create = H5Pcreate(H5P_DATASET_CREATE);
    dims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
    if( dims_mem > dims_max )
      dims_mem = dims_max;
    chkHDF5err(H5Pset_chunk(data_create, 1, &dims_mem));
    //---------------------------------------------------------------------
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
    //---------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
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
#endif//WALK_TREE_COMBINED_MODEL
    //---------------------------------------------------------------------
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
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
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
    //---------------------------------------------------------------------
    /* output # of particles in each process */
    dset = H5Dcreate(group, "num", H5T_NATIVE_INT, fulSpace, H5P_DEFAULT, data_create, H5P_DEFAULT);
    hyperslab = H5Dget_space(dset);
    chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
    chkHDF5err(H5Dwrite(dset, H5T_NATIVE_INT, locSpace, hyperslab, w_property, &num));
    chkHDF5err(H5Sclose(hyperslab));
    chkHDF5err(H5Dclose(dset));
    //---------------------------------------------------------------------
    /* close/release resources */
    chkHDF5err(H5Pclose(data_create));
    chkHDF5err(H5Sclose(locSpace));
    chkHDF5err(H5Sclose(fulSpace));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write attributes */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    int flag;
    //---------------------------------------------------------------------
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
    /* write flag about WALK_TREE_COMBINED_MODEL */
#ifdef  WALK_TREE_COMBINED_MODEL
    flag = 1;
#else///WALK_TREE_COMBINED_MODEL
    flag = 0;
#endif//WALK_TREE_COMBINED_MODEL
    attribute = H5Acreate(group, "WALK_TREE_COMBINED_MODEL", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about WALK_TREE_TOTAL_SUM_MODEL */
#ifdef  WALK_TREE_TOTAL_SUM_MODEL
    flag = 1;
#else///WALK_TREE_TOTAL_SUM_MODEL
    flag = 0;
#endif//WALK_TREE_TOTAL_SUM_MODEL
    attribute = H5Acreate(group, "WALK_TREE_TOTAL_SUM_MODEL", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
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
    /* write flag about LOCALIZE_I_PARTICLES */
#ifdef  LOCALIZE_I_PARTICLES
    flag = 1;
#else///LOCALIZE_I_PARTICLES
    flag = 0;
#endif//LOCALIZE_I_PARTICLES
    attribute = H5Acreate(group, "LOCALIZE_I_PARTICLES", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about USE_BRENT_METHOD */
#ifdef  USE_BRENT_METHOD
    flag = 1;
#else///USE_BRENT_METHOD
    flag = 0;
#endif//USE_BRENT_METHOD
    attribute = H5Acreate(group, "USE_BRENT_METHOD", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    /* write flag about TEST_BRENT_METHOD */
#ifdef  TEST_BRENT_METHOD
    flag = 1;
#else///TEST_BRENT_METHOD
    flag = 0;
#endif//TEST_BRENT_METHOD
    attribute = H5Acreate(group, "TEST_BRENT_METHOD", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &flag));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));
    //---------------------------------------------------------------------
  }/* if( steps != 0 ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* close the file */
  //-----------------------------------------------------------------------
  /* finish collective dataset write */
  chkHDF5err(H5Pclose(w_property));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
  *last ^= 1;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///USE_HDF5_FORMAT
//-------------------------------------------------------------------------
void readTentativeData(double *time, double *dt, ulong *steps, int num, iparticle body, char file[], int last)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s%d.dat", DATAFOLDER, file, TENTATIVE, last);
  fp = fopen(filename, "rb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  //-----------------------------------------------------------------------
  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fread( time, sizeof(double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(   dt, sizeof(double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(steps, sizeof( ulong), tmp, fp) )    success = false;
  /* tmp = num;  if( tmp != fread( body, sizeof(nbody_particle), tmp, fp) )    success = false; */
  tmp = num;  if( tmp != fread(body.pos, sizeof(    position), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.acc, sizeof(acceleration), tmp, fp) )    success = false;
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
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void readTentativeDataParallel(double *time, double *dt, ulong *steps, int *num, iparticle body, char file[], int last, MPIcfg_dataio *mpi, ulong Ntot)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* reset MPI_Offset */
  updateMPIcfg_dataio(mpi, *num);
  //-----------------------------------------------------------------------
  /* open the target file */
  char filename[128];
  sprintf(filename, "%s/%s.%s%d.dat", DATAFOLDER, file, TENTATIVE, last);
  MPI_File fh;
  chkMPIerr(MPI_File_open(mpi->comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  MPI_Status status;
  MPI_Offset disp = 0;
  //-----------------------------------------------------------------------
  /* the root process reads and broadcasts time, a real value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_read(fh, time, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_Bcast(time, 1, MPI_DOUBLE, 0, mpi->comm));
  disp += 1 * (MPI_Offset)sizeof(double);
  //-----------------------------------------------------------------------
  /* the root process reads and broadcasts dt, a real value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_read(fh, dt, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_Bcast(dt, 1, MPI_DOUBLE, 0, mpi->comm));
  disp += 1 * (MPI_Offset)sizeof(double);
  //-----------------------------------------------------------------------
  /* the root process reads and broadcasts write steps, an unsigned long value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_read(fh, steps, 1, MPI_UNSIGNED_LONG, &status));
  chkMPIerr(MPI_Bcast(steps, 1, MPI_UNSIGNED_LONG, 0, mpi->comm));
  disp += 1 * (MPI_Offset)sizeof(ulong);
  //-----------------------------------------------------------------------
  /* the whole processes read body, an iparticle array */
  /* chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(nbody_particle), mpi->body, mpi->body, "native", MPI_INFO_NULL)); */
  /* chkMPIerr(MPI_File_read(fh, body, num, mpi->body, &status)); */
  //-----------------------------------------------------------------------
  /* the whole processes read position */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(position), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.pos, (*num) * 4, MPI_REALDAT, &status));
  disp += Ntot * (MPI_Offset)sizeof(position);
  //-----------------------------------------------------------------------
  /* the whole processes read acceleration */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(acceleration), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.acc, (*num) * 4, MPI_REALDAT, &status));
  disp += Ntot * (MPI_Offset)sizeof(acceleration);
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* the whole processes read index */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(ulong), MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body.idx, *num, MPI_UNSIGNED_LONG, &status));
  disp += Ntot * (MPI_Offset)sizeof(ulong);
  //-----------------------------------------------------------------------
  /* close the target file */
  chkMPIerr(MPI_File_close(&fh));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void writeTentativeData(double time, double dt, ulong steps, ulong num, iparticle body, char file[], int *last)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s%d.dat", DATAFOLDER, file, TENTATIVE, (*last) ^ 1);
  fp = fopen(filename, "wb");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fwrite(& time, sizeof(        double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(&   dt, sizeof(        double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(&steps, sizeof(         ulong), tmp, fp) )    success = false;
  /* tmp = num;  if( tmp != fwrite(  body, sizeof(nbody_particle), tmp, fp) )    success = false; */
  tmp = num;  if( tmp != fwrite(body.pos, sizeof(    position), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.acc, sizeof(acceleration), tmp, fp) )    success = false;
#ifdef  BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body. vel, sizeof(  velocity), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.time, sizeof(ibody_time), tmp, fp) )    success = false;
#else///BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body.vx, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.vy, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.vz, sizeof(real), tmp, fp) )    success = false;
#endif//BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body.idx, sizeof(ulong), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }
  //-----------------------------------------------------------------------
  fclose(fp);
  *last ^= 1;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeTentativeDataParallel(double time, double dt, ulong steps, int num, iparticle body, char file[], int *last, MPIcfg_dataio *mpi, ulong Ntot)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* reset MPI_Offset */
  updateMPIcfg_dataio(mpi, num);
  //-----------------------------------------------------------------------
  /* open the target file */
  char filename[128];
  sprintf(filename, "%s/%s.%s%d.dat", DATAFOLDER, file, TENTATIVE, (*last) ^ 1);
  MPI_File fh;
  chkMPIerr(MPI_File_open(mpi->comm, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh));
  chkMPIerr(MPI_File_sync(fh));
  chkMPIerr(MPI_File_set_size(fh, 0));
  chkMPIerr(MPI_File_sync(fh));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  MPI_Status status;
  MPI_Offset disp = 0;
  //-----------------------------------------------------------------------
  /* the root process writes time, a real value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &time, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(double);
  //-----------------------------------------------------------------------
  /* the root process writes dt, a real value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &dt, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(double);
  //-----------------------------------------------------------------------
  /* the root process writes steps, an unsigned long value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &steps, 1, MPI_UNSIGNED_LONG, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(ulong);
  //-----------------------------------------------------------------------
  /* /\* the whole processes write body, an iparticle array *\/ */
  /* chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(nbody_particle), mpi->body, mpi->body, "native", MPI_INFO_NULL)); */
  /* chkMPIerr(MPI_File_write(fh, body, num, mpi->body, &status)); */
  /* chkMPIerr(MPI_File_sync(fh)); */
  //-----------------------------------------------------------------------
  /* the whole processes write position */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(position), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.pos, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(position);
  //-----------------------------------------------------------------------
  /* the whole processes write acceleration */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(acceleration), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.acc, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(acceleration);
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* the whole processes write index */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(ulong), MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.idx, num, MPI_UNSIGNED_LONG, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(ulong);
  //-----------------------------------------------------------------------
  /* close the output file */
  chkMPIerr(MPI_File_close(&fh));
  *last ^= 1;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
//-------------------------------------------------------------------------
void  readSnapshot(int *unit, double *time, ulong *steps, int num, nbody_hdf5 *body, char file[], uint id, hdf5struct type)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* open an existing file with read only option */
  //-----------------------------------------------------------------------
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, SNAPSHOT, id);
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  //-----------------------------------------------------------------------
  /* open an existing group */
  char groupname[16];
  sprintf(groupname, SNAPSHOT);
  hid_t group = H5Gopen(target, groupname, H5P_DEFAULT);
  //-----------------------------------------------------------------------
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
  /* read particle mass */
  dataset = H5Dopen(group, "mass", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->m));
  chkHDF5err(H5Dclose(dataset));
  /* read particle potential */
  dataset = H5Dopen(group, "potential", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->pot));
  chkHDF5err(H5Dclose(dataset));
  /* read particle index */
  dataset = H5Dopen(group, "index", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->idx));
  chkHDF5err(H5Dclose(dataset));
  //-----------------------------------------------------------------------
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
  /* read flag about DOUBLE_PRECISION */
  int useDP;
  attribute = H5Aopen(group, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* read unit system */
  //-----------------------------------------------------------------------
  group = H5Gopen(target, "unit system", H5P_DEFAULT);
  //-----------------------------------------------------------------------
  attribute = H5Aopen(group, "UNITSYSTEM", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, unit));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* close the file */
  //-----------------------------------------------------------------------
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* simple error check */
  //-----------------------------------------------------------------------
  if( num_ulong != (ulong)num ){
    __KILL__(stderr, "ERROR: number of N-body particles in the file (%zu) differs with that in the code (%d)\n", num_ulong, num);
  }
  //-----------------------------------------------------------------------
#ifdef  DOUBLE_PRECISION
  if( useDP != 1 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);
  }/* if( useDP != 1 ){ */
#else///DOUBLE_PRECISION
  if( useDP != 0 ){
    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);
  }/* if( useDP != 0 ){ */
#endif//DOUBLE_PRECISION
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeSnapshot(int  unit, double  time, ulong  steps, int num, nbody_hdf5 *body, char file[], uint id, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
		   , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
		   )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  char filename[128];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate total energy */
  //-----------------------------------------------------------------------
#ifdef  MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  double Ekin = 0.0;
  double Epot = 0.0;
  for(int ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    const double mass = (double)body->m  [ii        ];
    const double velx = (double)body->vel[ii * 3    ];
    const double vely = (double)body->vel[ii * 3 + 1];
    const double velz = (double)body->vel[ii * 3 + 2];
    Ekin += mass * (velx * velx + vely * vely + velz * velz);
    Epot += mass * (double)body->pot[ii];
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
  Ekin *= 0.5 * energy2astro;
  Epot *= 0.5 * energy2astro;
  const double Etot = Ekin + Epot;
  //-----------------------------------------------------------------------
  FILE *fp;
  sprintf(filename, "%s/%s.%s.log", LOGFOLDER, file, "energy");
  //-----------------------------------------------------------------------
  if( id == 0 ){
    //---------------------------------------------------------------------
    relEneErr->E0inv  = 1.0 / Etot;
    relEneErr->errMax = DBL_MIN;
    //---------------------------------------------------------------------
    fp = fopen(filename, "w");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    fprintf(fp, "#time(%s)\tsteps\tfile\trelative_error\tEtot(%s)\tEkin(%s)\tEpot(%s)\n", time_astro_unit_name, energy_astro_unit_name, energy_astro_unit_name, energy_astro_unit_name);
    //---------------------------------------------------------------------
  }/* if( id == 0 ){ */
  else{
    fp = fopen(filename, "a");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
  }/* else{ */
  //-----------------------------------------------------------------------
  const double Eerr = Etot * relEneErr->E0inv - 1.0;
  if( fabs(Eerr) > fabs(relEneErr->errMax) )
    relEneErr->errMax = Eerr;
  //-----------------------------------------------------------------------
  fprintf(fp, "%e\t%zu\t%u\t% e\t%e\t%e\t%e\n", time * time2astro, steps, id, Eerr, Etot, Ekin, Epot);
  fclose(fp);
  //-----------------------------------------------------------------------
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, SNAPSHOT, id);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create a group and data space */
  char groupname[16];
  sprintf(groupname, SNAPSHOT);
  hid_t group = H5Gcreate(target, groupname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  uint szip_pixels_per_block = 8;
  /* hsize_t cdims[2] = {128 * szip_pixels_per_block, 3}; */
  hsize_t cdims[2] = {32 * szip_pixels_per_block, 3};
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
  //-----------------------------------------------------------------------
  /* 2D (num * 3) array */
  hsize_t dims[2] = {num, 3};
  dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  hsize_t cdims_loc[2] = {cdims[0], cdims[1]};
  if( (hsize_t)num < cdims_loc[0] )
    cdims_loc[0] = (hsize_t)num;
  hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  if( cdims_loc[0] * cdims_loc[1] > cdims_max )
    cdims_loc[0] = cdims_max / cdims_loc[1];
  chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
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
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* 1D (num) arrays */
  dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
  property = H5Pcreate(H5P_DATASET_CREATE);
  cdims[0] = 128 * szip_pixels_per_block;
  cdims_loc[0] = cdims[0];
  if( (hsize_t)num < cdims_loc[0] )
    cdims_loc[0] = (hsize_t)num;
  cdims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
  if( cdims_loc[0] > cdims_max )
    cdims_loc[0] = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
  chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
  /* write particle mass */
  dataset = H5Dcreate(group, "mass", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->m));
  chkHDF5err(H5Dclose(dataset));
  /* write particle potential */
  dataset = H5Dcreate(group, "potential", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->pot));
  chkHDF5err(H5Dclose(dataset));
  /* write particle index */
  dataset = H5Dcreate(group, "index", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, body->idx));
  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
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
  /* write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
  attribute = H5Acreate(group, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* close the group */
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* write unit system */
  //-----------------------------------------------------------------------
  group = H5Gcreate(target, "unit system", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
  /* int unit = UNITSYSTEM; */
  attribute = H5Acreate(group, "UNITSYSTEM", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &unit));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "newton", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &newton));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* write conversion factors */
  hid_t subgroup = H5Gcreate(group, "conversion factors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(subgroup));
  //-----------------------------------------------------------------------
  /* write axis labels */
  subgroup = H5Gcreate(group, "axis labels", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(subgroup));
  //-----------------------------------------------------------------------
  /* close the group for unit system */
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  /* write energy conservation */
  //-----------------------------------------------------------------------
  group = H5Gcreate(target, "energy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* close the file */
  //-----------------------------------------------------------------------
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeSnapshotParallel(int  unit, double  time, ulong  steps, int num, nbody_hdf5 *body, char file[], uint id, MPIcfg_dataio *mpi, ulong Ntot, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
		   , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
			   )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  char filename[128];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate total energy */
  //-----------------------------------------------------------------------
#ifdef  MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  double Ekin = 0.0;
  double Epot = 0.0;
  for(int ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    const double mass = (double)body->m  [ii        ];
    const double velx = (double)body->vel[ii * 3    ];
    const double vely = (double)body->vel[ii * 3 + 1];
    const double velz = (double)body->vel[ii * 3 + 2];
    Ekin += mass * (velx * velx + vely * vely + velz * velz);
    Epot += mass * (double)body->pot[ii];
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &Ekin, 1, MPI_DOUBLE, MPI_SUM, mpi->comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &Epot, 1, MPI_DOUBLE, MPI_SUM, mpi->comm));
  //-----------------------------------------------------------------------
  Ekin *= 0.5 * energy2astro;
  Epot *= 0.5 * energy2astro;
  const double Etot = Ekin + Epot;
  //-----------------------------------------------------------------------
  if( id == 0 ){
    //---------------------------------------------------------------------
    relEneErr->E0inv  = 1.0 / Etot;
    relEneErr->errMax = DBL_MIN;
    //---------------------------------------------------------------------
  }/* if( id == 0 ){ */
  //-----------------------------------------------------------------------
  const double Eerr = Etot * relEneErr->E0inv - 1.0;
  if( fabs(Eerr) > fabs(relEneErr->errMax) )
    relEneErr->errMax = Eerr;
  //-----------------------------------------------------------------------
  if( mpi->rank == 0 ){
    //---------------------------------------------------------------------
    FILE *fp;
    sprintf(filename, "%s/%s.%s.log", LOGFOLDER, file, "energy");
    //---------------------------------------------------------------------
    if( id == 0 ){
      //-------------------------------------------------------------------
      fp = fopen(filename, "w");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
      fprintf(fp, "#time(%s)\tsteps\tfile\trelative_error\tEtot(%s)\tEkin(%s)\tEpot(%s)\n",
	      time_astro_unit_name, energy_astro_unit_name, energy_astro_unit_name, energy_astro_unit_name);
      //-------------------------------------------------------------------
    }/* if( id == 0 ){ */
    else{
      fp = fopen(filename, "a");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
    }/* else{ */
    //-----------------------------------------------------------------------
    fprintf(fp, "%e\t%zu\t%u\t% e\t%e\t%e\t%e\n", time * time2astro, steps, id, Eerr, Etot, Ekin, Epot);
    fclose(fp);
    //-----------------------------------------------------------------------
  }/* if( mpi->rank == 0 ){ */
  //-----------------------------------------------------------------------
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  //-----------------------------------------------------------------------
  hid_t f_property = H5Pcreate(H5P_FILE_ACCESS);
  chkHDF5err(H5Pset_fapl_mpio(f_property, mpi->comm, mpi->info));
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, SNAPSHOT, id);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, f_property);
  chkHDF5err(H5Pclose(f_property));
  //-----------------------------------------------------------------------
  /* create a group and data space */
  hid_t group = H5Gcreate(target, SNAPSHOT, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create (distributed) dataset */
  //-----------------------------------------------------------------------
  /* create dataspace */
  hsize_t dims_ful[2] = {Ntot            , 3};  hid_t fulSpace = H5Screate_simple(2, dims_ful, NULL);
  hsize_t dims_mem[2] = {Ntot / mpi->size, 3};
  hsize_t dims_loc[2] = { num            , 3};  hid_t locSpace = H5Screate_simple(2, dims_loc, NULL);
  //-----------------------------------------------------------------------
  /* create chunked dataset */
  hid_t data_create = H5Pcreate(H5P_DATASET_CREATE);
  hsize_t dims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  if( dims_mem[0] * dims_mem[1] > dims_max )
    dims_mem[0] = dims_max / dims_mem[1];
  chkHDF5err(H5Pset_chunk(data_create, 2, dims_mem));
  hid_t dataset, hyperslab;
  //-----------------------------------------------------------------------
#if 0
  /* slower than original implementation on dm, 8M bodies by 2 GPUs */
  /* reset size of chunk cache (the default size is 1 MiB) */
  hid_t data_access = H5Pcreate(H5P_DATASET_ACCESS);
  size_t chunk_size = dims_mem[0] * dims_mem[1] * sizeof(real);
  if( chunk_size > 1048576 )
    chkHDF5err(H5Pset_chunk_cache(data_access, H5D_CHUNK_CACHE_NSLOTS_DEFAULT, chunk_size, 1));
#else
  hid_t data_access = H5P_DEFAULT;
#endif
  //-----------------------------------------------------------------------
  /* configuration about domain decomposition */
  updateMPIcfg_dataio(mpi, num);
  hsize_t  count[2] = {1, 1};
  hsize_t stride[2] = {1, 1};
  hsize_t  block[2] = {dims_loc[0], dims_loc[1]};
  hsize_t offset[2] = {mpi->head, 0};
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: block = (%Lu, %Lu), offset = (%Lu, %Lu)\n", mpi->rank, block[0], block[1], offset[0], offset[1]);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  //-----------------------------------------------------------------------
  /* set up the collective transfer properties function */
  hid_t w_property = H5Pcreate(H5P_DATASET_XFER);
  chkHDF5err(H5Pset_dxpl_mpio(w_property, H5FD_MPIO_COLLECTIVE));
  //-----------------------------------------------------------------------
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
  /* close the dataspace */
  chkHDF5err(H5Pclose(data_access));
  chkHDF5err(H5Pclose(data_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));
  //-----------------------------------------------------------------------
  /* 1D (num) arrays */
  fulSpace = H5Screate_simple(1, dims_ful, NULL);
  locSpace = H5Screate_simple(1, dims_loc, NULL);
  data_create = H5Pcreate(H5P_DATASET_CREATE);
  /* chkHDF5err(H5Pset_chunk(data_create, 1, dims_loc)); */
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
  //-----------------------------------------------------------------------
  /* write attribute data */
  hsize_t attr_dims = 1;
  hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
  /* write current time */
#ifdef  DBG_PARALLEL_HDF5
  union {double d; ulong i;} buf_d;
  buf_d.d = time;
  fprintf(stdout, "rank %d: time = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  hid_t attribute = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time));
  chkHDF5err(H5Aclose(attribute));
  /* write # of steps */
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: steps = %lu\n", mpi->rank, steps);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  /* write # of N-body particles */
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: Ntot = %lu\n", mpi->rank, Ntot);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &Ntot));
  chkHDF5err(H5Aclose(attribute));
  /* write flag about DOUBLE_PRECISION */
#ifdef  DOUBLE_PRECISION
  const int useDP = 1;
#else///DOUBLE_PRECISION
  const int useDP = 0;
#endif//DOUBLE_PRECISION
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: useDP = %d\n", mpi->rank, useDP);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* close the group */
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* write unit system */
  //-----------------------------------------------------------------------
  group = H5Gcreate(target, "unit system", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
  attribute = H5Acreate(group, "UNITSYSTEM", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &unit));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
#ifdef  DOUBLE_PRECISION
  union {double r; ulong i;} buf_r;  buf_r.r = newton;
  fprintf(stdout, "rank %d: newton = %lx\n", mpi->rank, buf_r.i);
#else///DOUBLE_PRECISION
  union {float  r; uint  i;} buf_r;  buf_r.r = newton;
  fprintf(stdout, "rank %d: newton = %x\n", mpi->rank, buf_r.i);
#endif//DOUBLE_PRECISION
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "newton", type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.real, &newton));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* write conversion factors */
  hid_t subgroup = H5Gcreate(group, "conversion factors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = length2astro;
  fprintf(stdout, "rank %d: length2astro = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "length2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &length2astro));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = time2astro;
  fprintf(stdout, "rank %d: time2astro = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "time2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time2astro));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = mass2astro;
  fprintf(stdout, "rank %d: mass2astro = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "mass2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &mass2astro));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = density2astro;
  fprintf(stdout, "rank %d: density2astro = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "density2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &density2astro));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = col_density2astro;
  fprintf(stdout, "rank %d: col_density2astro = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "col_density2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &col_density2astro));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = energy2astro;
  fprintf(stdout, "rank %d: energy2astro = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "energy2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &energy2astro));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = senergy2astro;
  fprintf(stdout, "rank %d: senergy2astro = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "senergy2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &senergy2astro));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = velocity2astro;
  fprintf(stdout, "rank %d: velocity2astro = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "velocity2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &velocity2astro));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = accel2astro;
  fprintf(stdout, "rank %d: accel2astro = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "accel2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &accel2astro));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(subgroup));
  //-----------------------------------------------------------------------
  /* write axis labels */
  subgroup = H5Gcreate(group, "axis labels", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: length_astro_unit_name = \"%s\"\n", mpi->rank, length_astro_unit_name);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "length_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, length_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: time_astro_unit_name = \"%s\"\n", mpi->rank, time_astro_unit_name);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "time_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, time_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: mass_astro_unit_name = \"%s\"\n", mpi->rank, mass_astro_unit_name);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "mass_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, mass_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: density_astro_unit_name = \"%s\"\n", mpi->rank, density_astro_unit_name);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "density_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, density_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: col_density_astro_unit_name = \"%s\"\n", mpi->rank, col_density_astro_unit_name);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "col_density_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, col_density_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: energy_astro_unit_name = \"%s\"\n", mpi->rank, energy_astro_unit_name);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "energy_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, energy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: senergy_astro_unit_name = \"%s\"\n", mpi->rank, senergy_astro_unit_name);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "senergy_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, senergy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: velocity_astro_unit_name = \"%s\"\n", mpi->rank, velocity_astro_unit_name);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "velocity_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, velocity_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: accel_astro_unit_name = \"%s\"\n", mpi->rank, accel_astro_unit_name);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(subgroup, "accel_astro_unit_name", type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, type.str4unit, accel_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(subgroup));
  //-----------------------------------------------------------------------
  /* close the group for unit system */
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  /* write energy conservation */
  //-----------------------------------------------------------------------
  group = H5Gcreate(target, "energy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  //-----------------------------------------------------------------------
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = Etot;
  fprintf(stdout, "rank %d: Etot = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "total energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Etot));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = Ekin;
  fprintf(stdout, "rank %d: Ekin = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "kinetic energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Ekin));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = Epot;
  fprintf(stdout, "rank %d: Epot = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "potential energy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Epot));
  chkHDF5err(H5Aclose(attribute));
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = Eerr;
  fprintf(stdout, "rank %d: Eerr = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "relative error", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &Eerr));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* close the file */
  //-----------------------------------------------------------------------
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeSnapshotMultiGroups(double  time, ulong  steps, nbody_hdf5 *body, char file[], uint id, hdf5struct type, int kind, int *head, int *num)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  //-----------------------------------------------------------------------
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "split", id);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create the data space for the dataset */
  //-----------------------------------------------------------------------
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  uint szip_pixels_per_block = 8;
  /* hsize_t cdims[2] = {128 * szip_pixels_per_block, 3}; */
  hsize_t cdims[2] = {32 * szip_pixels_per_block, 3};
#else///USE_SZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* write attribute data */
  //-----------------------------------------------------------------------
  /* create the data space for the attribute */
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  //-----------------------------------------------------------------------
  /* write current time */
  double wtime = time * time2astro;
  attribute = H5Acreate(target, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &wtime));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write # of steps */
  attribute = H5Acreate(target, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write # of N-body particles */
  ulong num_ulong = 0;
  for(int ii = 0; ii < kind; ii++)
    num_ulong += (ulong)num[ii];
  attribute = H5Acreate(target, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &num_ulong));
  chkHDF5err(H5Aclose(attribute));
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
  //-----------------------------------------------------------------------
#ifdef  HDF5_FOR_ZINDAIJI
  //-----------------------------------------------------------------------
  /* write additional attribute for Zindaiji */
  static const char format_ver[CONSTANTS_H_CHAR_WORDS] = "0.0";
  attribute = H5Acreate(target, "Format Version", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, format_ver));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  real *tmp;  tmp = (real *)malloc(num_ulong * 3 * sizeof(real));  if( tmp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tmp");  }
  const double vel2zin = length2astro / time2astro;
  const double acc2zin = 0.5 * vel2zin / time2astro;
  //-----------------------------------------------------------------------
#endif//HDF5_FOR_ZINDAIJI
  //-----------------------------------------------------------------------
  chkHDF5err(H5Tclose(str4format));
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* write particle data */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < kind; ii++){
    //---------------------------------------------------------------------
    char grp[16];    sprintf(grp, "data%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //---------------------------------------------------------------------
    /* 2D (num * 3) array */
    hsize_t dims[2] = {num[ii], 3};
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    hsize_t cdims_loc[2] = {cdims[0], cdims[1]};
    if( (hsize_t)(num[ii] * 3) > (hsize_t)szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( (hsize_t)num[ii] < cdims_loc[0] )
	cdims_loc[0] = (hsize_t)num[ii];
      hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
      if( cdims_loc[0] * cdims_loc[1] > cdims_max )
	cdims_loc[0] = cdims_max / cdims_loc[1];
      chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( (hsize_t)(num[ii] * 3) > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
    /* write particle position */
    for(int jj = 3 * head[ii]; jj < 3 * (head[ii] + num[ii]); jj++)
      body->pos[jj] = (real)((double)body->pos[jj] * length2astro);
#ifdef  HDF5_FOR_ZINDAIJI
    dataset = H5Dcreate(group,      "xyz", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->pos[head[ii] * 3]));
    chkHDF5err(H5Dclose(dataset));
#endif//HDF5_FOR_ZINDAIJI
    /* coordinate transformation */
    dataset = H5Dcreate(group, "position", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->pos[head[ii] * 3]));
    chkHDF5err(H5Dclose(dataset));
    /* write particle velocity */
#ifdef  HDF5_FOR_ZINDAIJI
    for(int jj = 0; jj < 3 * num[ii]; jj++)
      tmp[jj] = (real)((double)body->vel[jj + 3 * head[ii]] * vel2zin);
    dataset = H5Dcreate(group,   "vxvyvz", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp));
    chkHDF5err(H5Dclose(dataset));
#endif//HDF5_FOR_ZINDAIJI
    /* coordinate transformation */
    for(int jj = 3 * head[ii]; jj < 3 * (head[ii] + num[ii]); jj++)
      body->vel[jj] = (real)((double)body->vel[jj] * velocity2astro);
    dataset = H5Dcreate(group, "velocity", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->vel[head[ii] * 3]));
    chkHDF5err(H5Dclose(dataset));
    /* write particle acceleration */
#ifdef  HDF5_FOR_ZINDAIJI
    /* write particle acceleration for Zindaiji */
    for(int jj = 0; jj < 3 * num[ii]; jj++)
      tmp[jj] = (real)((double)body->acc[jj + 3 * head[ii]] * acc2zin);
    dataset = H5Dcreate(group, "axayaz", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp));
    chkHDF5err(H5Dclose(dataset));
#endif//HDF5_FOR_ZINDAIJI
    /* coordinate transformation */
    for(int jj = 3 * head[ii]; jj < 3 * (head[ii] + num[ii]); jj++)
      body->acc[jj] = (real)((double)body->acc[jj] * accel2astro);
    dataset = H5Dcreate(group, "acceleration", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->acc[head[ii] * 3]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    if( (hsize_t)(num[ii] * 3) > (hsize_t)szip_pixels_per_block )
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------
    /* 1D (num) arrays */
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    cdims[0] = 128 * szip_pixels_per_block;
    cdims_loc[0] = cdims[0];
    if( (hsize_t)num[ii] < cdims_loc[0] )
      cdims_loc[0] = (hsize_t)num[ii];
    if( (hsize_t)num[ii] > (hsize_t)szip_pixels_per_block ){
      hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_8BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_8BIT : MAXIMUM_CHUNK_SIZE;
      if( cdims_loc[0] > cdims_max )
	cdims_loc[0] = cdims_max;
      chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( (hsize_t)num[ii] > (hsize_t)szip_pixels_per_block ){ */
    else
      property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
    /* write particle mass */
    for(int jj = head[ii]; jj < head[ii] + num[ii]; jj++)
      body->m[jj] = (real)((double)body->m[jj] * mass2astro);
    dataset = H5Dcreate(group, "mass", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->m[head[ii]]));
    chkHDF5err(H5Dclose(dataset));
    /* write particle potential */
    for(int jj = head[ii]; jj < head[ii] + num[ii]; jj++)
      body->pot[jj] = (real)((double)body->pot[jj] * senergy2astro);
    dataset = H5Dcreate(group, "potential", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->pot[head[ii]]));
    chkHDF5err(H5Dclose(dataset));
    /* write particle index */
    dataset = H5Dcreate(group, "index", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &body->idx[head[ii]]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    if( (hsize_t)num[ii] > (hsize_t)szip_pixels_per_block )
      chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* write attribute data */
    //---------------------------------------------------------------------
    /* create the data space for the attribute */
    hsize_t attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    //---------------------------------------------------------------------
    /* write # of N-body particles */
    ulong num_ulong = (ulong)num[ii];
    attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &num_ulong));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
#ifdef  HDF5_FOR_ZINDAIJI
    //---------------------------------------------------------------------
    /* write current time */
    float time_float = (float)(time * time2astro);
    attribute = H5Acreate(group, "time", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_FLOAT, &time_float));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
    /* write # of steps */
    attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
    chkHDF5err(H5Aclose(attribute));
    //---------------------------------------------------------------------
#endif//HDF5_FOR_ZINDAIJI
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    chkHDF5err(H5Gclose(group));
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < kind; ii++){ */
  //-----------------------------------------------------------------------
#ifdef  HDF5_FOR_ZINDAIJI
  free(tmp);
#endif//HDF5_FOR_ZINDAIJI
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  /* close the file */
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///USE_HDF5_FORMAT
//-------------------------------------------------------------------------
void readSnapshot(int *unit, double *time, ulong *steps, int num, iparticle body, char file[], uint id)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* read binary file */
  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s%.3u.dat", DATAFOLDER, file, SNAPSHOT, id);
  fp = fopen(filename, "rb");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fread( unit, sizeof(           int), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread( time, sizeof(        double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(steps, sizeof(         ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.pos, sizeof(    position), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.acc, sizeof(acceleration), tmp, fp) )    success = false;
#ifdef  BLOCK_TIME_STEP
  tmp = num;  if( tmp != fread(body.vel, sizeof(velocity), tmp, fp) )    success = false;
#else///BLOCK_TIME_STEP
  tmp = num;  if( tmp != fread(body.vx, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.vy, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread(body.vz, sizeof(real), tmp, fp) )    success = false;
#endif//BLOCK_TIME_STEP
  tmp = num;  if( tmp != fread(body.idx, sizeof(ulong), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeSnapshot(int  unit, double  time, ulong  steps, int  num, iparticle body, char file[], uint id)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s%.3u.dat", DATAFOLDER, file, SNAPSHOT, id);
  fp = fopen(filename, "wb");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fwrite(& unit, sizeof(           int), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(& time, sizeof(        double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(&steps, sizeof(         ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.pos, sizeof(    position), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.acc, sizeof(acceleration), tmp, fp) )    success = false;
#ifdef  BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body.vel, sizeof(velocity), tmp, fp) )    success = false;
#else///BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body.vx, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.vy, sizeof(real), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(body.vz, sizeof(real), tmp, fp) )    success = false;
#endif//BLOCK_TIME_STEP
  tmp = num;  if( tmp != fwrite(body.idx, sizeof(ulong), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
void writeSnapshotParallel(int  unit, double  time, ulong  steps, int  num, iparticle body, char file[], uint id, MPIcfg_dataio *mpi, ulong Ntot)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* reset MPI_Offset */
  updateMPIcfg_dataio(mpi, num);
  //-----------------------------------------------------------------------
  /* open the target file */
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.dat", DATAFOLDER, file, SNAPSHOT, id);
  MPI_File fh;
  chkMPIerr(MPI_File_open(mpi->comm, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh));
  chkMPIerr(MPI_File_sync(fh));
  chkMPIerr(MPI_File_set_size(fh, 0));
  chkMPIerr(MPI_File_sync(fh));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  MPI_Status status;
  MPI_Offset disp = 0;
  //-----------------------------------------------------------------------
  /* the root process writes time, a real value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_INT, MPI_INT, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &unit, 1, MPI_INT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(int);
  //-----------------------------------------------------------------------
  /* the root process writes time, a real value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &time, 1, MPI_DOUBLE, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(double);
  //-----------------------------------------------------------------------
  /* the root process writes steps, an unsigned long value */
  chkMPIerr(MPI_File_set_view(fh, disp, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  if( mpi->rank == 0 )
    chkMPIerr(MPI_File_write(fh, &steps, 1, MPI_UNSIGNED_LONG, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += 1 * (MPI_Offset)sizeof(ulong);
  //-----------------------------------------------------------------------
  /* the whole processes write position */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(position), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.pos, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(position);
  //-----------------------------------------------------------------------
  /* the whole processes write acceleration */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(acceleration), MPI_REALDAT, MPI_REALDAT, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.acc, num * 4, MPI_REALDAT, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(acceleration);
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* the whole processes write index */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(ulong), MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body.idx, num, MPI_UNSIGNED_LONG, &status));
  chkMPIerr(MPI_File_sync(fh));
  disp += Ntot * (MPI_Offset)sizeof(ulong);
  //-----------------------------------------------------------------------
  /* close the target file */
  chkMPIerr(MPI_File_close(&fh));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
//-------------------------------------------------------------------------
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void  readApproxAccel(double *time, ulong *step, int num, ulong *idx, acceleration *acc, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  FILE *fp;
  fp = fopen(file, "rb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", file);  }
  //-----------------------------------------------------------------------
  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fread(time, sizeof(      double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(step, sizeof(       ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread( idx, sizeof(       ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread( acc, sizeof(acceleration), tmp, fp) )    success = false;
  if( success != true ){
    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", file);
  }
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeApproxAccel(double  time, ulong  step, int num, ulong *idx, acceleration *acc, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  FILE *fp;
  fp = fopen(file, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", file);  }
  //-----------------------------------------------------------------------
  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fwrite(&time, sizeof(      double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fwrite(&step, sizeof(       ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(  idx, sizeof(       ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fwrite(  acc, sizeof(acceleration), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", file);  }
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
//-------------------------------------------------------------------------
void dumpBenchmark(int jobID, char file[], int steps, wall_clock_time *dat)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  //-----------------------------------------------------------------------
  wall_clock_time mean =
    {0.0, 0.0, 0.0,/* calcGravity_dev, calcMultipole, makeTree */
#ifdef  BLOCK_TIME_STEP
     0.0, 0.0, 0.0, 0.0,/* prediction_dev, correction_dev, setLaneTime_dev, adjustParticleTime_dev */
#else///BLOCK_TIME_STEP
     0.0, 0.0,/* advPos_dev, advVel_dev */
#endif//BLOCK_TIME_STEP
     0.0, 0.0,/* setTimeStep_dev, sortParticlesPHcurve */
     0.0, 0.0,/* copyParticle_dev2hst, copyParticle_hst2dev */
     0.0, 0.0/* setTreeNode_dev, setTreeCell_dev */
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
     , 0.0, 0.0/* examineNeighbor_dev, searchNeighbor_dev */
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  HUNT_MAKE_PARAMETER
     , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/* genPHkey_kernel, rsortKey_library, sortBody_kernel, makeTree_kernel, linkTree_kernel, trimTree_kernel */
     , 0.0, 0.0, 0.0, 0.0, 0.0/* initTreeLink_kernel, initTreeCell_kernel, initTreeNode_kernel, initTreeBody_kernel, copyRealBody_kernel */
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
     , 0.0, 0.0, 0.0, 0.0/* searchNeighbor_kernel, sortNeighbors, countNeighbors_kernel, commitNeighbors */
#endif//HUNT_FIND_PARAMETER
    };
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.%s%.8d.%s.dat", BENCH_LOG_FOLDER, file, WALLCLOCK, jobID, "bare");
  //-----------------------------------------------------------------------
  int newFile = (0 != access(filename, F_OK));
  //-----------------------------------------------------------------------
  fp = fopen(filename, "a");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  if( newFile ){
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", "step",
	    "calcGrav_dev", "calcMultipole", "makeTree", "setTimeStep_dev", "sortPHcurve",
	    "cpBody_dev2hst", "cpBody_hst2dev", "setTreeNode_dev", "setTreeCell_dev");
#ifdef  BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s\t%s\t%s", "prediction_dev", "correction_dev", "setLaneTime_dev", "adjustTime_dev");
#else///BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s", "advPos_dev", "advVel_dev");
#endif//BLOCK_TIME_STEP
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
    fprintf(fp, "\t%s\t%s", "examineNeighbor_dev", "searchNeighbor_dev");
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
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
    fprintf(fp, "%4d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e", ii,
	    dat[ii].calcGravity_dev, dat[ii].calcMultipole, dat[ii].makeTree, dat[ii].setTimeStep_dev, dat[ii].sortParticlesPHcurve,
	    dat[ii].copyParticle_dev2hst, dat[ii].copyParticle_hst2dev, dat[ii].setTreeNode_dev, dat[ii].setTreeCell_dev);
#ifdef  BLOCK_TIME_STEP
    fprintf(fp, "\t%e\t%e\t%e\t%e", dat[ii].prediction_dev, dat[ii].correction_dev, dat[ii].setLaneTime_dev, dat[ii].adjustParticleTime_dev);
#else///BLOCK_TIME_STEP
    fprintf(fp, "\t%e\t%e", dat[ii].advPos_dev, dat[ii].advVel_dev);
#endif//BLOCK_TIME_STEP
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
    fprintf(fp, "\t%e\t%e", dat[ii].examineNeighbor_dev, dat[ii].searchNeighbor_dev);
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  HUNT_MAKE_PARAMETER
    fprintf(fp, "\t%e\t%e\t%e\t%e\t%e\t%e", dat[ii].genPHkey_kernel, dat[ii].rsortKey_library, dat[ii].sortBody_kernel, dat[ii].makeTree_kernel, dat[ii].linkTree_kernel, dat[ii].trimTree_kernel);
    fprintf(fp, "\t%e\t%e\t%e\t%e\t%e", dat[ii].initTreeLink_kernel, dat[ii].initTreeCell_kernel, dat[ii].initTreeNode_kernel, dat[ii].initTreeBody_kernel, dat[ii].copyRealBody_kernel);
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
    fprintf(fp, "\t%e\t%e\t%e\t%e", dat[ii].searchNeighbor_kernel, dat[ii].sortNeighbors, dat[ii].countNeighbors_kernel, dat[ii].commitNeighbors);
#endif//HUNT_FIND_PARAMETER
    fprintf(fp, "\n");
    mean.calcGravity_dev      += dat[ii].calcGravity_dev;
    mean.calcMultipole        += dat[ii].calcMultipole;
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
    mean.setTreeNode_dev      += dat[ii].setTreeNode_dev;
    mean.setTreeCell_dev      += dat[ii].setTreeCell_dev;
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
    mean.examineNeighbor_dev += dat[ii].examineNeighbor_dev;
    mean. searchNeighbor_dev += dat[ii]. searchNeighbor_dev;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
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
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.%s%.8d.%s.dat", BENCH_LOG_FOLDER, file, WALLCLOCK, jobID, "mean");
  newFile = (0 != access(filename, F_OK));
  //-----------------------------------------------------------------------
  fp = fopen(filename, "a");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  //-----------------------------------------------
  if( newFile ){
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
	    "calcGrav_dev", "calcMultipole", "makeTree", "setTimeStep_dev", "sortPHcurve",
	    "cpBody_dev2hst", "cpBody_hst2dev", "setTreeNode_dev", "setTreeCell_dev");
#ifdef  BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s\t%s\t%s", "prediction_dev", "correction_dev", "setLaneTime_dev", "adjustTime_dev");
#else///BLOCK_TIME_STEP
    fprintf(fp, "\t%s\t%s", "advPos_dev", "advVel_dev");
#endif//BLOCK_TIME_STEP
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
    fprintf(fp, "\t%s\t%s", "examineNeighbor_dev", "searchNeighbor_dev");
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
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
  fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e",
	  inv * mean.calcGravity_dev, inv * mean.calcMultipole, inv * mean.makeTree, inv * mean.setTimeStep_dev, inv * mean.sortParticlesPHcurve,
	  inv * mean.copyParticle_dev2hst, inv * mean.copyParticle_hst2dev, inv * mean.setTreeNode_dev, inv * mean.setTreeCell_dev);
#ifdef  BLOCK_TIME_STEP
  fprintf(fp, "\t%e\t%e\t%e\t%e", inv * mean.prediction_dev, inv * mean.correction_dev, inv * mean.setLaneTime_dev, inv * mean.adjustParticleTime_dev);
#else///BLOCK_TIME_STEP
  fprintf(fp, "\t%e\t%e", inv * mean.advPos_dev, inv * mean.advVel_dev);
#endif//BLOCK_TIME_STEP
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
  fprintf(fp, "\t%e\t%e", inv * mean.examineNeighbor_dev, inv * mean.searchNeighbor_dev);
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  HUNT_MAKE_PARAMETER
  fprintf(fp, "\t%e\t%e\t%e\t%e\t%e\t%e", inv * mean.genPHkey_kernel, inv * mean.rsortKey_library, inv * mean.sortBody_kernel, inv * mean.makeTree_kernel, inv * mean.linkTree_kernel, inv * mean.trimTree_kernel);
  fprintf(fp, "\t%e\t%e\t%e\t%e\t%e", inv * mean.initTreeLink_kernel, inv * mean.initTreeCell_kernel, inv * mean.initTreeNode_kernel, inv * mean.initTreeBody_kernel, inv * mean.copyRealBody_kernel);
#endif//HUNT_MAKE_PARAMETER
#ifdef  HUNT_FIND_PARAMETER
  fprintf(fp, "\t%e\t%e\t%e\t%e", inv * mean.searchNeighbor_kernel, inv * mean.sortNeighbors, inv * mean.countNeighbors_kernel, inv * mean.commitNeighbors);
#endif//HUNT_FIND_PARAMETER
  fprintf(fp, "\n");
  fclose(fp);
  //-----------------------------------------------------------------------
#ifdef  HUNT_WALK_PARAMETER
  sprintf(filename, "%s/%s.%s.dat", BENCH_LOG_FOLDER, file, "walk");
  newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  if( newFile ){
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", "calcGravity_dev", "NTHREADS", "TSUB", "NWARP", "NLOOP",
	    "USE_WARP_SHUFFLE_FUNC", "GET_BUFID", "LOCALIZE_I_PARTICLES");
#ifdef  NEIGHBOR_PHKEY_LEVEL
    fprintf(fp, "\t%s", "NEIGHBOR_PHKEY_LEVEL");
#endif//NEIGHBOR_PHKEY_LEVEL
#ifdef  BRUTE_FORCE_LOCALIZATION
#       ifdef  USE_BRENT_METHOD
    fprintf(fp, "\t%s", "NEIGHBOR_LENGTH_SHRINK_FACTOR");
#       else///USE_BRENT_METHOD
    fprintf(fp, "\t%s\t%s", "NEIGHBOR_LENGTH_UPPER_FRACTION", "NEIGHBOR_LENGTH_EXTEND_FRACTION");
#       endif//USE_BRENT_METHOD
#endif//BRUTE_FORCE_LOCALIZATION
    fprintf(fp, "\n");
  }/* if( newFile ){ */
  fprintf(fp, "%e\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
	  inv * mean.calcGravity_dev,
	  NTHREADS, TSUB, NWARP, NLOOP
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
	  /* localize i-particles or not */
#ifdef  LOCALIZE_I_PARTICLES
	  , 1
#else///LOCALIZE_I_PARTICLES
	  , 0
#endif//LOCALIZE_I_PARTICLES
	  );
  /* criterion about PH-key level to split i-particles */
#ifdef  NEIGHBOR_PHKEY_LEVEL
  fprintf(fp, "\t%d", NEIGHBOR_PHKEY_LEVEL);
#endif//NEIGHBOR_PHKEY_LEVEL
  /* criterion about neighbor length to split i-particles */
#ifdef  BRUTE_FORCE_LOCALIZATION
#       ifdef  USE_BRENT_METHOD
  fprintf(fp, "\t%e", NEIGHBOR_LENGTH_SHRINK_FACTOR);
#       else///USE_BRENT_METHOD
  fprintf(fp, "\t%e\t%e", NEIGHBOR_LENGTH_UPPER_FRACTION, NEIGHBOR_LENGTH_EXTEND_FRACTION);
#       endif//USE_BRENT_METHOD
#endif//BRUTE_FORCE_LOCALIZATION
  fprintf(fp, "\n");
  fclose(fp);
#endif//HUNT_WALK_PARAMETER
  //-----------------------------------------------------------------------
#   if  defined(HUNT_MAKE_PARAMETER) || (defined(HUNT_FIND_PARAMETER) && defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES))
  extern int treeBuildCalls;
  const double invSteps = 1.0 / (double)treeBuildCalls;
#endif//defined(HUNT_MAKE_PARAMETER) || (defined(HUNT_FIND_PARAMETER) && defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES))
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  sprintf(filename, "%s/%s.%s.dat", BENCH_LOG_FOLDER, file, "make");
  newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  if( newFile ){
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s", "makeTree", "sortParticlesPHcurve", "genPHkey_kernel", "rsortKey_library", "NTHREADS_PH", "sortBody_kernel", "NTHREADS_PHSORT");
#ifdef  MAKE_TREE_ON_DEVICE
  fprintf(fp, "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
	  "makeTree_kernel", "linkTree_kernel", "trimTree_kernel",
	  "NTHREADS_MAKE_TREE", "TSUB_MAKE_TREE", "NTHREADS_LINK_TREE", "NTHREADS_TRIM_TREE",
	  "USE_WARP_SHUFFLE_FUNC_MAKE_TREE", "USE_WARP_SHUFFLE_FUNC_LINK_TREE");
#endif//MAKE_TREE_ON_DEVICE
  fprintf(fp, "\t%s\t%s\t%s\t%s\t%s\t%s", "initTreeLink_kernel", "initTreeCell_kernel", "initTreeNode_kernel", "initTreeBody_kernel", "copyRealBody_kernel", "Ttot");
  fprintf(fp, "\n");
  }/* if( newFile ){ */
  fprintf(fp, "%e\t%e\t%e\t%e\t%d\t%e\t%d", invSteps * mean.makeTree, invSteps * mean.sortParticlesPHcurve,
	  invSteps * mean.genPHkey_kernel, invSteps * mean.rsortKey_library, NTHREADS_PH, invSteps * mean.sortBody_kernel, NTHREADS_PHSORT);
#ifdef  MAKE_TREE_ON_DEVICE
  fprintf(fp, "\t%e\t%e\t%e\t%d\t%d\t%d\t%d\t%d\t%d",
	  invSteps * mean.makeTree_kernel, invSteps * mean.linkTree_kernel, invSteps * mean.trimTree_kernel,
	  NTHREADS_MAKE_TREE, TSUB_MAKE_TREE, NTHREADS_LINK_TREE, NTHREADS_TRIM_TREE
	  /* use warp shuffle for making tree structure instruction or not */
#       ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_TREE
	  , 1
#       else///USE_WARP_SHUFFLE_FUNC_MAKE_TREE
	  , 0
#       endif//USE_WARP_SHUFFLE_FUNC_MAKE_TREE
	  /* use warp shuffle for linking tree structure instruction or not */
#       ifdef  USE_WARP_SHUFFLE_FUNC_LINK_TREE
	  , 1
#       else///USE_WARP_SHUFFLE_FUNC_LINK_TREE
	  , 0
#       endif//USE_WARP_SHUFFLE_FUNC_LINK_TREE
	  );
#endif//MAKE_TREE_ON_DEVICE
  fprintf(fp, "\t%e\t%e\t%e\t%e\t%e\t%d", invSteps * mean.initTreeLink_kernel, invSteps * mean.initTreeCell_kernel, invSteps * mean.initTreeNode_kernel, inv * mean.initTreeBody_kernel, inv * mean.copyRealBody_kernel, NTHREADS_INIT_BODY);
  fprintf(fp, "\n");
  fclose(fp);
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------
#ifdef  HUNT_NODE_PARAMETER
  sprintf(filename, "%s/%s.%s.dat", BENCH_LOG_FOLDER, file, "node");
  newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  if( newFile )
    fprintf(fp, "#%s\t%s\t%s\t%s\n",
	    "calcMultipole", "NTHREADS_MAC", "TSUB_MAC", "USE_WARP_SHUFFLE_FUNC_MAC");
  fprintf(fp, "%e\t%d\t%d\t%d\n",
	  inv * mean.calcMultipole,
	  NTHREADS_MAC, TSUB_MAC
	  /* use warp shuffle instruction or not */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
	  , 1
#else///USE_WARP_SHUFFLE_FUNC_MAC
	  , 0
#endif//USE_WARP_SHUFFLE_FUNC_MAC
	  );
  fclose(fp);
#endif//HUNT_NODE_PARAMETER
  //-----------------------------------------------------------------------
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
    fprintf(fp, "\t%s\t%s\n", "NTHREADS_TIME", "USE_WARP_SHUFFLE_FUNC_TIME");
  }/* if( newFile ){ */
#ifdef  BLOCK_TIME_STEP
  fprintf(fp, "%e\t%e\t%e\t%e\t%e",
	  inv * mean.prediction_dev, inv * mean.correction_dev, inv * mean.setLaneTime_dev,
	  inv * mean.adjustParticleTime_dev, inv * mean.setTimeStep_dev);
#else///BLOCK_TIME_STEP
  fprintf(fp, "%e\t%e\t%e", inv * mean.advPos_dev, inv * mean.advVel_dev, inv * mean.setTimeStep_dev);
#endif//BLOCK_TIME_STEP
  fprintf(fp, "\t%d\t%d\n", NTHREADS_TIME
#ifdef  USE_WARP_SHUFFLE_FUNC_TIME
	  , 1
#else///USE_WARP_SHUFFLE_FUNC_TIME
	  , 0
#endif//USE_WARP_SHUFFLE_FUNC_TIME
	  );
  fclose(fp);
#endif//HUNT_TIME_PARAMETER
  //-----------------------------------------------------------------------
#   if  defined(HUNT_FIND_PARAMETER) && defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
  sprintf(filename, "%s/%s.%s.dat", BENCH_LOG_FOLDER, file, "find");
  newFile = (0 != access(filename, F_OK));
  fp = fopen(filename, "a");  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  if( newFile )
    fprintf(fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	    "examineNeighbor_dev", "searchNeighbor_dev", "searchNeighbor_kernel", "sortNeighbors", "countNeighbors_kernel", "commitNeighbors",
	    "NTHREADS_SHRINK", "NTHREADS_FACILE_NS", "NTHREADS_NEIGHBOR", "TSUB_NEIGHBOR", "NEIGHBOR_NUM",
	    "SMEM_PREF_FOR_NEIGHBOR_SEARCH", "USE_WARP_SHUFFLE_FUNC_NEIGHBOR");
  fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
	  invSteps * mean.examineNeighbor_dev, invSteps * mean.searchNeighbor_dev, invSteps * mean.searchNeighbor_kernel, invSteps * mean.sortNeighbors, invSteps * mean.countNeighbors_kernel, invSteps * mean.commitNeighbors,
	  NTHREADS_SHRINK, NTHREADS_FACILE_NS, NTHREADS_NEIGHBOR, TSUB_NEIGHBOR, NEIGHBOR_NUM
#ifdef  SMEM_PREF_FOR_NEIGHBOR_SEARCH
	  , 1
#else///SMEM_PREF_FOR_NEIGHBOR_SEARCH
	  , 0
#endif//SMEM_PREF_FOR_NEIGHBOR_SEARCH
#ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  , 1
#else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  , 0
#endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
	  );
  fclose(fp);
#endif//defined(HUNT_FIND_PARAMETER) && defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//EXEC_BENCHMARK
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
//-------------------------------------------------------------------------
void dumpStatistics(int jobID, char file[], int steps, int PHlevel, tree_metrics *dat)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.%s%.8d.%s.dat", BENCH_LOG_FOLDER, file, WALK_STAT, jobID, "Nj");
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }/* if( fp == NULL ){ */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.%s%.8d.%s.dat", BENCH_LOG_FOLDER, file, WALK_STAT, jobID, "Nbuf");
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }/* if( fp == NULL ){ */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  sprintf(filename, "%s/%s.%s%.8d.%s.dat", BENCH_LOG_FOLDER, file, TREE_STAT, jobID, "total");
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }/* if( fp == NULL ){ */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int ll = 0; ll < PHlevel; ll++){
    //---------------------------------------------------------------------
    sprintf(filename, "%s/%s.%s%.8d.%s%.2d.dat", BENCH_LOG_FOLDER, file, TREE_STAT, jobID, "level", ll);
    fp = fopen(filename, "w");
    if( fp == NULL ){
      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
    }/* if( fp == NULL ){ */
    //---------------------------------------------------------------------
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
    //---------------------------------------------------------------------
    fclose(fp);
    //---------------------------------------------------------------------
  }/* for(int ll = 0; ll < PHlevel; ll++){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//COUNT_INTERACTIONS
//-------------------------------------------------------------------------
