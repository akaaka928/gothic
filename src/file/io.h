/**
 * @file io.h
 *
 * @brief Header file for Input/Output functions in GOTHIC and MAGI
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/05/09 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef IO_H
#define IO_H


#include <stdbool.h>
#include "macro.h"
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#include <mpilib.h>
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)

#include "../misc/structure.h"
#ifndef RUN_WITHOUT_GOTHIC
#include "../misc/benchmark.h"
#include "../misc/tune.h"
#include "../misc/brent.h"
#endif//RUN_WITHOUT_GOTHIC

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#endif//USE_HDF5_FORMAT


/** macros to specify the file name */
#define CONFIG    "cfg.dat"
#define SETTINGS  "run.dat"
#define TENTATIVE "tmp"
#define SNAPSHOT  "snp"
#define BRUTE "direct"
#define APPBH "tree.bh"
#define APPWS "tree.ws"
#define APPG2 "tree.g2"
#define ERRFILE_NUM "err.num"
#define ERRFILE_LST "err.lst"
#define ERRFILE_CDF "err.cdf"
#define ACCERR "accerr"
#define TREE_STAT "stat.tree"
#define WALK_STAT "stat.walk"

#ifndef USE_HDF5_FORMAT
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
#define READ_SUPERPOSED_TABLE_SPHE (-1)
#define READ_SUPERPOSED_TABLE_DISK (-2)
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#endif//USE_HDF5_FORMAT


/**
 * @struct MPIcfg_dataio
 *
 * @brief structure for MPI-IO
 */
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
typedef struct
{
  MPI_Offset head;
  MPI_Comm comm;
  MPI_Info info;
  int rank, size;
} MPIcfg_dataio;
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)

#ifdef  USE_HDF5_FORMAT
/**
 * @struct hdf5struct
 *
 * @brief structure for HDF5
 */
typedef struct
{
  hid_t real;
  hid_t str4unit;
#ifndef RUN_WITHOUT_GOTHIC
  hid_t rebuildTree, measuredTime, statVal, guessTime;/**< to take over the current status of auto-tuning about tree rebuild interval */
  hid_t brentFunc, brentStatus, brentMemory;/**< to take over the current status of auto-tuning to remove low-dense regions */
#endif//RUN_WITHOUT_GOTHIC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  hid_t pot2;
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  hid_t disk_grav;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
  hid_t gpu_clock;
#endif//REPORT_GPU_CLOCK_FREQUENCY
} hdf5struct;
#endif//USE_HDF5_FORMAT

/* list of functions appeared in ``io.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void createMPIcfg_dataio(MPIcfg_dataio *cfg, const MPIinfo mpi);
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)

  void updateConfigFile        (int  last, char file[]);
  void   readConfigFile        (int *last, char file[]);
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void   readConfigFileParallel(int *last, char file[], MPIinfo mpi);
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)

  void readSettings        (int *unit, ulong *Ntot, real *eps, real *eta, double *ft, double *SnapshotInterval, double *SaveInterval, char file[]);
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void readSettingsParallel(int *unit, ulong *Ntot, real *eps, real *eta, double *ft, double *SnapshotInterval, double *SaveInterval, char file[], MPIinfo mpi);
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void writeSettings       (int  unit, ulong  Ntot, real  eps, real  eta, double  ft, double  SnapshotInterval, double  SaveInterval, char file[]);

#ifdef  USE_HDF5_FORMAT
  void createHDF5DataType(hdf5struct *type);
  void removeHDF5DataType(hdf5struct  type);
#endif//USE_HDF5_FORMAT

  /* read/write particle data */
  void  readTentativeData(double *time, double *dt, ulong *steps, int num, iparticle body, char file[], int  last
#ifdef  USE_HDF5_FORMAT
			  , hdf5struct type
#ifndef RUN_WITHOUT_GOTHIC
			  , int *dropPrevTune, rebuildTree *rebuild, measuredTime *measured, autoTuningParam *rebuildParam, brentStatus *status, brentMemory *memory
#ifdef  MONITOR_ENERGY_ERROR
			  , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
			  );
  void writeTentativeData
  (double  time, double  dt, ulong  steps, ulong num, iparticle body, char file[], int *last
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
   );

  /* read/write particle data with MPI-IO */
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void  readTentativeDataParallel(double *time, double *dt, ulong *steps, int *num, iparticle body, char file[], int  last, MPIcfg_dataio *mpi, ulong Ntot
#ifdef  USE_HDF5_FORMAT
				  , hdf5struct type
#ifndef RUN_WITHOUT_GOTHIC
				  , int *dropPrevTune, rebuildTree *rebuild, measuredTime *measured, autoTuningParam *rebuildParam, brentStatus *status, brentMemory *memory
#ifdef  MONITOR_ENERGY_ERROR
				  , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
				  );
  void writeTentativeDataParallel
  (double  time, double  dt, ulong  steps, int num, iparticle body, char file[], int *last, MPIcfg_dataio *mpi, ulong Ntot
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
   );
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)

  /* write particle data in TIPSY format */
  void writeTipsyFile(double time, float eps, int num, iparticle body, char file[]);

  /* write particle data in GalactICS format */
  void writeGalactICSFile(double time, int head, int num, iparticle body, char file[], int series);

  /* write particle data in GADGET format */
  void writeGADGETFile
  (const int Ntot, double time, int kind, int skind, int * head, int * num, iparticle body, char file[]
#ifdef  USE_HDF5_FORMAT
   , hdf5struct type
#endif//USE_HDF5_FORMAT
   );


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  /* write fixed potential field */
  void writeFixedPotentialTable
  (const int unit, potential_field pot_tbl_sphe, const int skind, potential_field *pot_tbl,
#ifdef  USE_HDF5_FORMAT
   hdf5struct type,
#else///USE_HDF5_FORMAT
   const bool binary,
#endif//USE_HDF5_FORMAT
   char file[]);

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
  void writeFixedDiskPotential
  (const int unit, const disk_potential disk
#ifdef  USE_HDF5_FORMAT
   , hdf5struct type
#else///USE_HDF5_FORMAT
   , const bool binary
#endif//USE_HDF5_FORMAT
   , char *file);
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  /* read/write snapshot generated by GOTHIC */
  void  readSnapshot
  (int *unit, double *time, ulong *steps, int num, char file[], uint id
#ifdef  USE_HDF5_FORMAT
   , nbody_hdf5 *body, hdf5struct type
#else///USE_HDF5_FORMAT
   , iparticle body
#endif//USE_HDF5_FORMAT
   );
  void writeSnapshot
  (int  unit, double  time, ulong  steps, int num, char file[], uint id
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
   );

#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
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
   );

#ifdef  USE_HDF5_FORMAT
  void writeSnapshotMultiGroups(double  time, ulong  steps, nbody_hdf5 *body, char file[], uint id, hdf5struct type, int kind, int *head, int *num);
#endif//USE_HDF5_FORMAT
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)

  /* read/write particles' acceleration calculated by GOTHIC */
  void  readApproxAccel(double *time, ulong *step, int num, ulong *idx, acceleration *acc, char file[]);
  void writeApproxAccel(double  time, ulong  step, int num, ulong *idx, acceleration *acc, char file[]);

#ifdef  EXEC_BENCHMARK
  /* write results of measured execution time of GOTHIC */
  void dumpBenchmark(int jobID, char file[], int steps, wall_clock_time *dat);
#endif//EXEC_BENCHMARK

#ifdef  COUNT_INTERACTIONS
  /* write # of interactions in GOTHIC */
  void dumpStatistics(int jobID, char file[], int steps, int PHlevel, tree_metrics *dat);
#endif//COUNT_INTERACTIONS
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//IO_H
