/*************************************************************************\
 *                                                                       *
                  last updated on 2016/11/05(Sat) 15:38:37
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
#include <mpi.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <mpilib.h>
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#include "../misc/tune.h"
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#include "../misc/brent.h"
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
//-------------------------------------------------------------------------
#include "../para/mpicfg.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* constants to set particle sampling rate */
#define DEFAULT_SAMPLING_RATE   (1.0e-4f)
#define DEFAULT_SAMPLING_NUMBER (64.0f)
/* #define DEFAULT_SAMPLING_NUMBER (512.0f) */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define NDIM (3)
//-------------------------------------------------------------------------
typedef struct
{
  real min[NDIM], max[NDIM];
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
#ifdef  REDUCE_EXCHANGING
//-------------------------------------------------------------------------
/* constant to detect load imbalance */
static const double loadImbalanceCrit = 0.9;
#define SLOW_DOWN_PROCS_CRIT(tot) ((tot) >> 4)
//-------------------------------------------------------------------------
typedef struct
{
  MPI_Comm comm;
  double invsize;
  double tmax, tmin;
  bool enable, execute;
} loadImbalanceDetector;
//-------------------------------------------------------------------------
#endif//REDUCE_EXCHANGING
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
  (domainCfg **box, sendCfg **sendBuf, recvCfg **recvBuf, MPIcfg_tree *mpi,
#ifdef  REDUCE_EXCHANGING
   MPIinfo orb[],
#else///REDUCE_EXCHANGING
   int *ndim, MPIinfo **orb,
#endif//REDUCE_EXCHANGING
   const ulong Ntot, real *samplingRate, int *sampleNumMax, real **sampleLoc, real **sampleFul,
   int **recvNum, int **recvDsp, real **boxMin, real **boxMax);
  void removeORBtopology
  (domainCfg  *box, sendCfg  *sendBuf, recvCfg  *recvBuf,
#ifndef REDUCE_EXCHANGING
   const int ndim, MPIinfo  *orb,
#endif//REDUCE_EXCHANGING
   real  *sampleLoc, real  *sampleFul, int  *recvNum, int  *recvDsp, real  *boxMin, real  *boxMax);
  //-----------------------------------------------------------------------
  void exchangeParticles
  (const int numOld, const int numMax, int *numNew, iparticle ibody_dev, iparticle src, iparticle dst,
   domainDecomposeKey * restrict key,
   sendCfg *sendBuf, recvCfg *recvBuf,
   const float samplingRate, const int sampleNumMax, const double tloc, real *sampleLoc, real *sampleFul, int *recvNum, int *recvDsp, real *boxMin, real *boxMax,
#ifdef  REDUCE_EXCHANGING
   double *exchangeInterval, loadImbalanceDetector balancer[], MPIinfo orb[],
#else///REDUCE_EXCHANGING
   const int ndim, MPIinfo *orb,
#endif//REDUCE_EXCHANGING
   domainCfg *domain, MPIcfg_tree mpi, measuredTime *measured
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
   , brentStatus *status, brentMemory *memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
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
