/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/05(Tue) 17:03:26
 *                                                                       *
 *    Header File for Input/Output Code of N-body simulation             *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef IO_H
#define IO_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#       ifndef MPILIB_H
#include <mpilib.h>
#       endif//MPILIB_H
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
//-------------------------------------------------------------------------
#ifndef BENCHMARK_H
#       include "../misc/benchmark.h"
#endif//BENCHMARK_H
//-------------------------------------------------------------------------
#ifndef STRUCTURE_H
#       include "../misc/structure.h"
#endif//STRUCTURE_H
//-------------------------------------------------------------------------
#   if  defined(USE_HDF5_FORMAT) && !defined(_HDF5_H)
#       include <hdf5.h>
#endif//defined(USE_HDF5_FORMAT) && !defined(_HDF5_H)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define IORETRY (16)
//-------------------------------------------------------------------------
#define CONFIG    "cfg.dat"
#define TENTATIVE "tmp"
#define SNAPSHOT  "snp"
/* #define APPROX    "app" */
#define BRUTE "direct"
#define APPBH "tree.bh"
#define APPWS "tree.ws"
#define APPG2 "tree.g2"
#define ERRFILE_NUM "err.num"
#define ERRFILE_LST "err.lst"
#define ERRFILE_CDF "err.cdf"
//-------------------------------------------------------------------------
#define TREE_STAT "stat.tree"
#define WALK_STAT "stat.walk"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
typedef struct
{
  MPI_Datatype body;
  MPI_Offset   head;
  MPI_Comm comm;
  MPI_Info info;
  int rank, size;
} MPIcfg_dataio;
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
typedef struct
{
  hid_t nbody;
  hid_t position;
  hid_t acceleration;
#ifdef  BLOCK_TIME_STEP
  hid_t velocity;
  hid_t ibody_time;
#endif//BLOCK_TIME_STEP
  hid_t real;
  hid_t str4unit;
} hdf5struct;
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "io.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void createMPIcfg_dataio(MPIcfg_dataio *cfg, MPIinfo mpi);
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  //-----------------------------------------------------------------------
  void updateConfigFile        (int  last, char file[]);
  void   readConfigFile        (int *last, char file[]);
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void   readConfigFileParallel(int *last, char file[], MPIinfo mpi);
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  //-----------------------------------------------------------------------
  void readSettings(int *unit, ulong *Ntot, real *eps, real *eta, double *ft, double *SnapshotInterval, double *SaveInterval, char file[]);
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void readSettingsParallel(int *unit, ulong *Ntot, real *eps, real *eta, double *ft, double *SnapshotInterval, double *SaveInterval, char file[], MPIinfo mpi);
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void writeSettings(int unit, ulong Ntot, real eps, real eta, double ft, double SnapshotInterval, double SaveInterval, char file[]);
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  void createHDF5DataType(hdf5struct *type);
  void removeHDF5DataType(hdf5struct  type);
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  void  readTentativeData        (double *time, double *dt, ulong *steps, int num, nbody_particle *body, char file[], int  last
#ifdef  USE_HDF5_FORMAT
				  , hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
				  , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
				  );
  void writeTentativeData        (double  time, double  dt, ulong  steps, int num, nbody_particle *body, char file[], int *last
#ifdef  USE_HDF5_FORMAT
				  , hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
				  , energyError relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
				  );
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void  readTentativeDataParallel(double *time, double *dt, ulong *steps, int num, nbody_particle *body, char file[], int  last, MPIcfg_dataio *mpi
#ifdef  USE_HDF5_FORMAT
				  , ulong Ntot, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
				  , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
				  );
  void writeTentativeDataParallel(double  time, double  dt, ulong  steps, int num, nbody_particle *body, char file[], int *last, MPIcfg_dataio *mpi
#ifdef  USE_HDF5_FORMAT
				  , ulong Ntot, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
				  , energyError  relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
				  );
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  void  readSnapshot(int *unit, double *time, ulong *steps, int num, nbody_hdf5 *body, char file[], uint id, hdf5struct type);
  void writeSnapshot(int  unit, double  time, ulong  steps, int num, nbody_hdf5 *body, char file[], uint id, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
		     , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
		     );
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void writeSnapshotParallel(int  unit, double  time, ulong  steps, int num, nbody_hdf5 *body, char file[], uint id, MPIcfg_dataio *mpi, ulong Ntot, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
			     , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
			     );
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void writeSnapshotMultiGroups(double  time, ulong  steps, nbody_hdf5 *body, char file[], uint id, hdf5struct type, int kind, int *head, int *num);
#else///USE_HDF5_FORMAT
  void  readSnapshot(int *unit, double *time, ulong *steps, int num, nbody_particle *body, char file[], uint id);
  void writeSnapshot(int  unit, double  time, ulong  steps, int num, nbody_particle *body, char file[], uint id);
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
  void writeSnapshotParallel(int  unit, double  time, ulong  steps, int num, nbody_particle *body, char file[], uint id, MPIcfg_dataio *mpi
#ifdef  USE_HDF5_FORMAT
			     , ulong Ntot, hdf5struct type
#ifdef  MONITOR_ENERGY_ERROR
			     , energyError *relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
			     );
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  /* void  readApproxAccel(double *time, ulong *steps, int num, ulong *idx, acceleration *direct, acceleration *octree, char file[], uint id); */
  /* void writeApproxAccel(double  time, ulong  steps, int num, ulong *idx, acceleration *direct, acceleration *octree, char file[], uint id); */
  void  readApproxAccel(double *time, ulong *step, int num, ulong *idx, acceleration *acc, char file[]);
  void writeApproxAccel(double  time, ulong  step, int num, ulong *idx, acceleration *acc, char file[]);
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  void dumpBenchmark(int jobID, char file[], int steps, wall_clock_time *dat);
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
  void dumpStatistics(int jobID, char file[], int steps, int PHlevel, tree_metrics *dat);
#endif//COUNT_INTERACTIONS
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//IO_H
//-------------------------------------------------------------------------
