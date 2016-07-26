/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/26(Tue) 16:39:45
 *                                                                       *
 *    Header File for N-body calculation with MPI parallelization        *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef EXCHANGE_H
#define EXCHANGE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  !defined(MPI_INCLUDED) && !defined(OMPI_MPI_H)
#       include <mpi.h>
#endif//!defined(MPI_INCLUDED) && !defined(OMPI_MPI_H)
//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef MPILIB_H
#       include <mpilib.h>
#endif//MPILIB_H
//-------------------------------------------------------------------------
#ifndef STRUCTURE_H
#       include "../misc/structure.h"
#endif//STRUCTURE_H
//-------------------------------------------------------------------------
#ifndef MPICFG_H
#       include "../para/mpicfg.h"
#endif//MPICFG_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* constants to set particle sampling rate */
//-------------------------------------------------------------------------
#define DEFAULT_SAMPLING_RATE   (1.0e-4f)
/* #define DEFAULT_SAMPLING_NUMBER (512.0f) */
#define DEFAULT_SAMPLING_NUMBER (64.0f)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
  real min[3], max[3];
  MPI_Request req;
} domainCfg;
//-------------------------------------------------------------------------
typedef struct
{
  real xmin, xmax, ymin, ymax, zmin, zmax;
  int rank, num, head;
  MPI_Request pos, idx;
#ifdef  GADGET_MAC
  MPI_Request acc;
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
  MPI_Request vel, time;
#else///BLOCK_TIME_STEP
  MPI_Request vx, vy, vz;
#endif//BLOCK_TIME_STEP
} sendCfg;
//-------------------------------------------------------------------------
typedef struct
{
  int num, head;
  MPI_Request req, pos, idx;
#ifdef  GADGET_MAC
  MPI_Request acc;
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
  MPI_Request vel, time;
#else///BLOCK_TIME_STEP
  MPI_Request vx, vy, vz;
#endif//BLOCK_TIME_STEP
} recvCfg;
//-------------------------------------------------------------------------
typedef struct
{
  int dstRank, bodyIdx;
} domainDecomposeKey;
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
//-- List of functions appeared in "exchange.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  void configORBtopology
  (int *ndim, MPIinfo **orb, domainCfg **box, sendCfg **sendBuf, recvCfg **recvBuf, MPIcfg_tree *mpi,
   const ulong Ntot, real *samplingRate, int *sampleNumMax, real **sampleLoc, real **sampleFul,
   int **recvNum, int **recvDsp, real **boxMin, real **boxMax);
  void removeORBtopology
  (const int ndim, MPIinfo  *orb, domainCfg  *box, sendCfg  *sendBuf, recvCfg  *recvBuf,
   real  *sampleLoc, real  *sampleFul, int  *recvNum, int  *recvDsp, real  *boxMin, real  *boxMax);
  //-----------------------------------------------------------------------
  void exchangeParticles
  (const int numOld,              iparticle src,
#   if  defined(GENERATE_PHKEY_ON_DEVICE) || !defined(CALC_MULTIPOLE_ON_DEVICE)
   iparticle ibody_dev,
#endif//defined(GENERATE_PHKEY_ON_DEVICE) || !defined(CALC_MULTIPOLE_ON_DEVICE)
   const int numMax, int *numNew, iparticle dst,
   domainDecomposeKey * restrict key,
   sendCfg *sendBuf, recvCfg *recvBuf,
   const float samplingRate, const int sampleNumMax, const double tloc, real *sampleLoc, real *sampleFul, int *recvNum, int *recvDsp, real *boxMin, real *boxMax,
   const int ndim, MPIinfo *orb, domainCfg *domain, MPIcfg_tree mpi
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
   );
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//EXCHANGE_H
//-------------------------------------------------------------------------
