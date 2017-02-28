/**
 * @file io.h
 *
 * @brief Header file for Input/Output functions in GOTHIC and MAGI
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/02/27 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef IO_H
#define IO_H


#include "macro.h"
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#include <mpilib.h>
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)

#include "../misc/structure.h"
#ifndef RUN_WITHOUT_GOTHIC
#include "../misc/benchmark.h"
#include "../misc/tune.h"
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD) && !defined(BRENT_H)
#include "../misc/brent.h"
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD) && !defined(BRENT_H)
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
} hdf5struct;
#endif//USE_HDF5_FORMAT


/* list of functions appeared in ``io.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void createMPIcfg_dataio(MPIcfg_dataio *cfg, MPIinfo mpi);
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
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
			  );
  void writeTentativeData(double  time, double  dt, ulong  steps, ulong num, iparticle body, char file[], int *last
#ifdef  USE_HDF5_FORMAT
			  , hdf5struct type
#ifndef RUN_WITHOUT_GOTHIC
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
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
			  );

  /* read/write particle data with MPI-IO */
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void  readTentativeDataParallel(double *time, double *dt, ulong *steps, int *num, iparticle body, char file[], int  last, MPIcfg_dataio *mpi, ulong Ntot
#ifdef  USE_HDF5_FORMAT
				  , hdf5struct type
#ifndef RUN_WITHOUT_GOTHIC
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
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
				  );
  void writeTentativeDataParallel(double  time, double  dt, ulong  steps, int num, iparticle body, char file[], int *last, MPIcfg_dataio *mpi, ulong Ntot
#ifdef  USE_HDF5_FORMAT
				  , hdf5struct type
#ifndef RUN_WITHOUT_GOTHIC
				  , rebuildTree rebuild, measuredTime measured
#ifdef  WALK_TREE_COMBINED_MODEL
				  , autoTuningParam rebuildParam
#endif//WALK_TREE_COMBINED_MODEL
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
				  , brentStatus status, brentMemory memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  MONITOR_ENERGY_ERROR
				  , energyError  relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//RUN_WITHOUT_GOTHIC
#endif//USE_HDF5_FORMAT
				  );
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)

  /* write particle data in TIPSY format */
  void writeTipsyFile(double time, float eps, int num, iparticle body, char file[]);

  /* read/write snapshot generated by GOTHIC */
  void  readSnapshot(int *unit, double *time, ulong *steps, int num, char file[], uint id
#ifdef  USE_HDF5_FORMAT
		     , nbody_hdf5 *body, hdf5struct type
#else///USE_HDF5_FORMAT
		     , iparticle body
#endif//USE_HDF5_FORMAT
		     );
  void writeSnapshot(int  unit, double  time, ulong  steps, int num, char file[], uint id
#ifdef  USE_HDF5_FORMAT
		     , nbody_hdf5 *body, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
		     , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
#else///USE_HDF5_FORMAT
		     , iparticle body
#endif//USE_HDF5_FORMAT
		     );

#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void writeSnapshotParallel(int  unit, double  time, ulong  steps, int num, char file[], uint id, MPIcfg_dataio *mpi, ulong Ntot
#ifdef  USE_HDF5_FORMAT
			     , nbody_hdf5 *body, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
			     , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
#else///USE_HDF5_FORMAT
			     , iparticle body
#endif//USE_HDF5_FORMAT
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
