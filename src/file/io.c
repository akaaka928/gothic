/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/22(Fri) 17:21:34
 *                                                                       *
 *    Input/Output Code of N-body calculation                            *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#define USE_SZIP_COMPRESSION
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
#       include <hdf5.h>
#       include <hdf5lib.h>
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
#include <macro.h>
#include <name.h>
#include <mpilib.h>
#include <constants.h>
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
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
  commitMPIbyte(&(cfg->body), (int)sizeof(nbody_particle));
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
  chkHDF5err(H5Tset_size(type->str4unit, CONSTANTS_H_CHAR_WORDS));
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
  /* commit data type of nbody_particle */
  //-----------------------------------------------------------------------
  type->nbody = H5Tcreate(H5T_COMPOUND, sizeof(nbody_particle));
  chkHDF5err(H5Tinsert(type->nbody, "index", HOFFSET(nbody_particle, idx), H5T_NATIVE_ULONG));
#ifdef  BLOCK_TIME_STEP
  chkHDF5err(H5Tinsert(type->nbody, "time", HOFFSET(nbody_particle, t0), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->nbody, "t + dt", HOFFSET(nbody_particle, t1), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->nbody, "dt", HOFFSET(nbody_particle, dt), type->real));
#endif//BLOCK_TIME_STEP
  chkHDF5err(H5Tinsert(type->nbody,   "x", HOFFSET(nbody_particle,   x), type->real));
  chkHDF5err(H5Tinsert(type->nbody,   "y", HOFFSET(nbody_particle,   y), type->real));
  chkHDF5err(H5Tinsert(type->nbody,   "z", HOFFSET(nbody_particle,   z), type->real));
  chkHDF5err(H5Tinsert(type->nbody,  "vx", HOFFSET(nbody_particle,  vx), type->real));
  chkHDF5err(H5Tinsert(type->nbody,  "vy", HOFFSET(nbody_particle,  vy), type->real));
  chkHDF5err(H5Tinsert(type->nbody,  "vz", HOFFSET(nbody_particle,  vz), type->real));
  chkHDF5err(H5Tinsert(type->nbody,  "ax", HOFFSET(nbody_particle,  ax), type->real));
  chkHDF5err(H5Tinsert(type->nbody,  "ay", HOFFSET(nbody_particle,  ay), type->real));
  chkHDF5err(H5Tinsert(type->nbody,  "az", HOFFSET(nbody_particle,  az), type->real));
  chkHDF5err(H5Tinsert(type->nbody,   "m", HOFFSET(nbody_particle,   m), type->real));
  chkHDF5err(H5Tinsert(type->nbody, "pot", HOFFSET(nbody_particle, pot), type->real));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* commit data type of position, acceleration, velocity, and ibody_time */
  //-----------------------------------------------------------------------
  type->position     = H5Tcreate(H5T_COMPOUND, sizeof(position));
  type->acceleration = H5Tcreate(H5T_COMPOUND, sizeof(acceleration));
#ifdef  BLOCK_TIME_STEP
  type->velocity     = H5Tcreate(H5T_COMPOUND, sizeof(velocity));
  type->ibody_time   = H5Tcreate(H5T_COMPOUND, sizeof(ibody_time));
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
  chkHDF5err(H5Tinsert(type->    position,   "x", HOFFSET(    position,   x), type->real));
  chkHDF5err(H5Tinsert(type->    position,   "y", HOFFSET(    position,   y), type->real));
  chkHDF5err(H5Tinsert(type->    position,   "z", HOFFSET(    position,   z), type->real));
  chkHDF5err(H5Tinsert(type->    position,   "m", HOFFSET(    position,   m), type->real));
  chkHDF5err(H5Tinsert(type->acceleration,   "x", HOFFSET(acceleration,   x), type->real));
  chkHDF5err(H5Tinsert(type->acceleration,   "y", HOFFSET(acceleration,   y), type->real));
  chkHDF5err(H5Tinsert(type->acceleration,   "z", HOFFSET(acceleration,   z), type->real));
  chkHDF5err(H5Tinsert(type->acceleration, "pot", HOFFSET(acceleration, pot), type->real));
#ifdef  BLOCK_TIME_STEP
  chkHDF5err(H5Tinsert(type->    velocity,   "x", HOFFSET(    velocity,   x), type->real));
  chkHDF5err(H5Tinsert(type->    velocity,   "y", HOFFSET(    velocity,   y), type->real));
  chkHDF5err(H5Tinsert(type->    velocity,   "z", HOFFSET(    velocity,   z), type->real));
  chkHDF5err(H5Tinsert(type->    velocity,  "dt", HOFFSET(    velocity,  dt), type->real));
  chkHDF5err(H5Tinsert(type->ibody_time,   "time", HOFFSET(ibody_time, t0), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(type->ibody_time, "t + dt", HOFFSET(ibody_time, t1), H5T_NATIVE_DOUBLE));
#endif//BLOCK_TIME_STEP
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
  chkHDF5err(H5Tclose(type.str4unit));
  chkHDF5err(H5Tclose(type.nbody));
  chkHDF5err(H5Tclose(type.position));
  chkHDF5err(H5Tclose(type.acceleration));
#ifdef  BLOCK_TIME_STEP
  chkHDF5err(H5Tclose(type.velocity));
  chkHDF5err(H5Tclose(type.ibody_time));
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  readTentativeData(double *time, double *dt, ulong *steps, int num, nbody_particle *body, char file[], int  last, hdf5struct type
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
  /* read particle data */
  hid_t dataset = H5Dopen(group, "data", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, type.nbody, H5S_ALL, H5S_ALL, H5P_DEFAULT, body));
  chkHDF5err(H5Dclose(dataset));
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
  /* close the file */
  //-----------------------------------------------------------------------
  /* end access to the dataset and release the corresponding resource */
  chkHDF5err(H5Gclose(group));
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
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  readTentativeDataParallel(double *time, double *dt, ulong *steps, int num, nbody_particle *body, char file[], int  last, MPIcfg_dataio *mpi, ulong Ntot, hdf5struct type
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
  /* reset MPI_Offset */
  updateMPIcfg_dataio(mpi, num);
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
  hid_t dataset = H5Dopen(group, "data", H5P_DEFAULT);
  //-----------------------------------------------------------------------
  /* configuration about domain decomposition */
  hsize_t dims_loc = num;
  hid_t locSpace = H5Screate_simple(1, &dims_loc, NULL);
  hsize_t  count = 1;
  hsize_t stride = 1;
  hsize_t  block = dims_loc;
  hsize_t offset = mpi->head;
  //-----------------------------------------------------------------------
  /* select hyperslab in the file */
  hid_t hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  //-----------------------------------------------------------------------
  /* create property list for collective dataset read */
  hid_t r_property = H5Pcreate(H5P_DATASET_XFER);
  chkHDF5err(H5Pset_dxpl_mpio(r_property, H5FD_MPIO_COLLECTIVE));
  //-----------------------------------------------------------------------
  /* read particle data */
  chkHDF5err(H5Dread(dataset, type.nbody, locSpace, hyperslab, r_property, body));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Pclose(r_property));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(locSpace));
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
  /* close the file */
  //-----------------------------------------------------------------------
  /* end access to the dataset and release the corresponding resource */
  chkHDF5err(H5Gclose(group));
  //-----------------------------------------------------------------------
  chkHDF5err(H5Fclose(target));
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
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeTentativeData(double  time, double  dt, ulong  steps, int num, nbody_particle *body, char file[], int *last, hdf5struct type
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
  /* create the data space for the dataset */
  //-----------------------------------------------------------------------
  hsize_t dims = num;
  hid_t dataspace = H5Screate_simple(1, &dims, NULL);
  //-----------------------------------------------------------------------
  /* SZIP compression can only be used with atomic datatypes that are integer, float, or char. */
  /* It cannot be applied to compound, array, variable-length, enumerations, or other user-defined datatypes. */
  /* The call to H5Dcreate will fail if attempting to create an SZIP compressed dataset with a non-allowed datatype. */
  /* The conflict can only be detected when the property list is used. */
  //-----------------------------------------------------------------------
  /* create the dataset */
  hid_t dataset = H5Dcreate(group, "data", type.nbody, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------
  /* write particle data */
  chkHDF5err(H5Dwrite(dataset, type.nbody, H5S_ALL, H5S_ALL, H5P_DEFAULT, body));
  chkHDF5err(H5Dclose(dataset));
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
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* close the file */
  //-----------------------------------------------------------------------
  chkHDF5err(H5Gclose(group));
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
  *last ^= 1;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeTentativeDataParallel(double  time, double  dt, ulong  steps, int num, nbody_particle *body, char file[], int *last, MPIcfg_dataio *mpi, ulong Ntot, hdf5struct type
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

  //-----------------------------------------------------------------------
  /* create (distributed) dataset */
  //-----------------------------------------------------------------------
  /* create dataspace */
  hsize_t dims_ful = Ntot            ;  hid_t fulSpace = H5Screate_simple(1, &dims_ful, NULL);
  hsize_t dims_mem = Ntot / mpi->size;
  hsize_t dims_loc =  num            ;  hid_t locSpace = H5Screate_simple(1, &dims_loc, NULL);
  //-----------------------------------------------------------------------
  /* create chunked dataset */
  hid_t date_create = H5Pcreate(H5P_DATASET_CREATE);
  /* chkHDF5err(H5Pset_chunk(date_create, 1, &dims_loc)); */
  chkHDF5err(H5Pset_chunk(date_create, 1, &dims_mem));
  hid_t dset = H5Dcreate(group, "data", type.nbody, fulSpace, H5P_DEFAULT, date_create, H5P_DEFAULT);
  //-----------------------------------------------------------------------
  /* configuration about domain decomposition */
  updateMPIcfg_dataio(mpi, num);
  hsize_t  count = 1;
  hsize_t stride = 1;
  hsize_t  block = dims_loc;
  hsize_t offset = mpi->head;
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: block = %Lu, offset = %Lu\n", mpi->rank, block, offset);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  //-----------------------------------------------------------------------
  /* select hyperslab in the file */
  hid_t hyperslab = H5Dget_space(dset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &offset, &stride, &count, &block));
  //-----------------------------------------------------------------------
  /* create property list for collective dataset write */
  hid_t w_property = H5Pcreate(H5P_DATASET_XFER);
  chkHDF5err(H5Pset_dxpl_mpio(w_property, H5FD_MPIO_COLLECTIVE));
  chkHDF5err(H5Dwrite(dset, type.nbody, locSpace, hyperslab, w_property, body));
  //-----------------------------------------------------------------------
  /* close/release resources */
  chkHDF5err(H5Pclose(w_property));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dset));
  chkHDF5err(H5Pclose(date_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));
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
#ifdef  DBG_PARALLEL_HDF5
  union {double d; ulong i;} buf_d;
  buf_d.d = time;
  fprintf(stdout, "rank %d: time = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write time step */
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = dt;
  fprintf(stdout, "rank %d: dt = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "dt", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &dt));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write # of steps */
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: steps = %lu\n", mpi->rank, steps);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write # of N-body particles */
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: Ntot = %lu\n", mpi->rank, Ntot);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &Ntot));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* write inverse of total energy at the initial condition */
#ifdef  MONITOR_ENERGY_ERROR
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = relEneErr.E0inv;
  fprintf(stdout, "rank %d: E0inv = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "E0inv", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &relEneErr.E0inv));
  chkHDF5err(H5Aclose(attribute));
#endif//MONITOR_ENERGY_ERROR
  //-----------------------------------------------------------------------
  /* write maximum value of the relative error of the total energy */
#ifdef  MONITOR_ENERGY_ERROR
#ifdef  DBG_PARALLEL_HDF5
  buf_d.d = relEneErr.errMax;
  fprintf(stdout, "rank %d: errMax = %lx\n", mpi->rank, buf_d.i);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
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
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: useDP = %d\n", mpi->rank, useDP);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
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
#ifdef  DBG_PARALLEL_HDF5
  fprintf(stdout, "rank %d: blockTimeStep = %d\n", mpi->rank, blockTimeStep);
  fflush(stdout);
#endif//DBG_PARALLEL_HDF5
  attribute = H5Acreate(group, "blockTimeStep", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &blockTimeStep));
  chkHDF5err(H5Aclose(attribute));
  //-----------------------------------------------------------------------
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* close the file */
  //-----------------------------------------------------------------------
  chkHDF5err(H5Gclose(group));
  chkHDF5err(H5Fclose(target));
  //-----------------------------------------------------------------------
  *last ^= 1;
  //-----------------------------------------------------------------------
#ifdef  DBG_PARALLEL_HDF5
#if 1
  MPI_Finalize();
  exit(0);
#endif
#endif//DBG_PARALLEL_HDF5
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///USE_HDF5_FORMAT
//-------------------------------------------------------------------------
void readTentativeData(double *time, double *dt, ulong *steps, int num, nbody_particle *body, char file[], int last)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  char filename[128];
  FILE *fp;
  sprintf(filename, "%s/%s.%s%d.dat", DATAFOLDER, file, TENTATIVE, last);
  fp = fopen(filename, "rb");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  bool success = true;
  size_t tmp;
  tmp =   1;  if( tmp != fread( time, sizeof(        double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(   dt, sizeof(        double), tmp, fp) )    success = false;
  tmp =   1;  if( tmp != fread(steps, sizeof(         ulong), tmp, fp) )    success = false;
  tmp = num;  if( tmp != fread( body, sizeof(nbody_particle), tmp, fp) )    success = false;
  if( success != true ){
    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void readTentativeDataParallel(double *time, double *dt, ulong *steps, int num, nbody_particle *body, char file[], int last, MPIcfg_dataio *mpi)
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
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(nbody_particle), mpi->body, mpi->body, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_read(fh, body, num, mpi->body, &status));
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
void writeTentativeData(double time, double dt, ulong steps, int num, nbody_particle *body, char file[], int *last)
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
  tmp = num;  if( tmp != fwrite(  body, sizeof(nbody_particle), tmp, fp) )    success = false;
  if( success != true ){
    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  fclose(fp);
  *last ^= 1;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeTentativeDataParallel(double time, double dt, ulong steps, int num, nbody_particle *body, char file[], int *last, MPIcfg_dataio *mpi)
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
  /* the whole processes write body, an iparticle array */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(nbody_particle), mpi->body, mpi->body, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body, num, mpi->body, &status));
  chkMPIerr(MPI_File_sync(fh));
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
  hid_t date_create = H5Pcreate(H5P_DATASET_CREATE);
  chkHDF5err(H5Pset_chunk(date_create, 2, dims_mem));
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
  dataset = H5Dcreate(group, "position", type.real, fulSpace, H5P_DEFAULT, date_create, data_access);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->pos));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* write particle velocity */
  dataset = H5Dcreate(group, "velocity", type.real, fulSpace, H5P_DEFAULT, date_create, data_access);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->vel));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* write particle acceleration */
  dataset = H5Dcreate(group, "acceleration", type.real, fulSpace, H5P_DEFAULT, date_create, data_access);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->acc));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* close the dataspace */
  chkHDF5err(H5Pclose(data_access));
  chkHDF5err(H5Pclose(date_create));
  chkHDF5err(H5Sclose(locSpace));
  chkHDF5err(H5Sclose(fulSpace));
  //-----------------------------------------------------------------------
  /* 1D (num) arrays */
  fulSpace = H5Screate_simple(1, dims_ful, NULL);
  locSpace = H5Screate_simple(1, dims_loc, NULL);
  date_create = H5Pcreate(H5P_DATASET_CREATE);
  /* chkHDF5err(H5Pset_chunk(date_create, 1, dims_loc)); */
  chkHDF5err(H5Pset_chunk(date_create, 1, dims_mem));
  /* write particle mass */
  dataset = H5Dcreate(group, "mass", type.real, fulSpace, H5P_DEFAULT, date_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->m));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* write particle potential */
  dataset = H5Dcreate(group, "potential", type.real, fulSpace, H5P_DEFAULT, date_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, type.real, locSpace, hyperslab, w_property, body->pot));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* write particle index */
  dataset = H5Dcreate(group, "index", H5T_NATIVE_ULONG, fulSpace, H5P_DEFAULT, date_create, H5P_DEFAULT);
  hyperslab = H5Dget_space(dataset);
  chkHDF5err(H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset, stride, count, block));
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_ULONG, locSpace, hyperslab, w_property, body->idx));
  chkHDF5err(H5Sclose(hyperslab));
  chkHDF5err(H5Dclose(dataset));
  /* close the dataspace */
  chkHDF5err(H5Pclose(date_create));
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
#ifdef  DBG_PARALLEL_HDF5
#if 1
  MPI_Finalize();
  exit(0);
#endif
#endif//DBG_PARALLEL_HDF5
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
    if( (uint)(num[ii] * 3) > szip_pixels_per_block ){
      property = H5Pcreate(H5P_DATASET_CREATE);
      if( (hsize_t)num[ii] < cdims_loc[0] )
	cdims_loc[0] = (hsize_t)num[ii];
      chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( (uint)(num[ii] * 3) > szip_pixels_per_block ){ */
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
    if( (uint)(num[ii] * 3) > szip_pixels_per_block )
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
    if( (uint)num[ii] > szip_pixels_per_block ){
      chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
    }/* if( (uint)num[ii] > szip_pixels_per_block ){ */
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
    if( (uint)num[ii] > szip_pixels_per_block )
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
void readSnapshot(int *unit, double *time, ulong *steps, int num, nbody_particle *body, char file[], uint id)
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
  tmp = num;  if( tmp != fread( body, sizeof(nbody_particle), tmp, fp) )    success = false;
  if( success != true ){
    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------
  /* setPhysicalConstantsAndUnitSystem(*unit, 0); */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void writeSnapshot(int  unit, double  time, ulong  steps, int  num, nbody_particle *body, char file[], uint id)
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
  tmp = num;  if( tmp != fwrite(  body, sizeof(nbody_particle), tmp, fp) )    success = false;
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
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
void writeSnapshotParallel(int  unit, double  time, ulong  steps, int  num, nbody_particle *body, char file[], uint id, MPIcfg_dataio *mpi)
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
  /* the whole processes write body, an iparticle array */
  chkMPIerr(MPI_File_set_view(fh, disp + mpi->head * (MPI_Offset)sizeof(nbody_particle), mpi->body, mpi->body, "native", MPI_INFO_NULL));
  chkMPIerr(MPI_File_write(fh, body, num, mpi->body, &status));
  chkMPIerr(MPI_File_sync(fh));
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
