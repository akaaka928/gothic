/*************************************************************************\
 *                                                                       *
                  last updated on 2017/01/17(Tue) 12:08:00
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
#include "macro.h"
#include "mpilib.h"
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#ifdef  EXCHANGE_USING_GPUS
#include "../misc/tune.h"
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#include "../misc/brent.h"
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#endif//EXCHANGE_USING_GPUS
//-------------------------------------------------------------------------
#include "../para/mpicfg.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* constants to set particle sampling rate */
#define DEFAULT_SAMPLING_RATE   (1.0e-4f)
/* #define DEFAULT_SAMPLING_NUMBER (64.0f) */
#define DEFAULT_SAMPLING_NUMBER (128.0f)
/* #define DEFAULT_SAMPLING_NUMBER (256.0f) */
/* #define DEFAULT_SAMPLING_NUMBER (512.0f) */
//-------------------------------------------------------------------------
typedef struct
{
  int *rnum, *disp;/* former recvNum and recvDsp */
  float *xmin, *xmax, *ymin, *ymax, *zmin, *zmax;
  float rate;/* former samplingRate */
  int Nmax;/* former sampleNumMax */
} sampling;
//-------------------------------------------------------------------------
typedef struct
{
  float *x_hst, *y_hst, *z_hst;
  float *x_dev, *y_dev, *z_dev;
  int *i_hst, *i_dev;
} samplePos;
//-------------------------------------------------------------------------
typedef struct
{
  float *x, *y, *z;
} particlePos;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
#ifdef  EXCHANGE_USING_GPUS
  float *xmin, *xmax, *ymin, *ymax, *zmin, *zmax;
  MPI_Request *req;
#else///EXCHANGE_USING_GPUS
  real min[3], max[3];
  MPI_Request req;
#endif//EXCHANGE_USING_GPUS
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
#ifdef  EXCHANGE_USING_GPUS
  int *dstRank_hst, *dstRank_dev, *bodyIdx_dev;
#else///EXCHANGE_USING_GPUS
  int dstRank, bodyIdx;
#endif//EXCHANGE_USING_GPUS
} domainDecomposeKey;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* constant to detect load imbalance */
/* static const double loadImbalanceCrit = 0.9; */
static const double loadImbalanceCrit = 0.95;
#define SLOW_DOWN_PROCS_CRIT(tot) ((tot) >> 4)
//-------------------------------------------------------------------------
typedef struct
{
  double tmin, tmax;
  bool enable, execute;
} loadImbalanceDetector;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "exchange.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
#ifdef  EXCHANGE_USING_GPUS
  //-----------------------------------------------------------------------
  muse allocateORMtopology(float **dxmin, float **dxmax, float **dymin, float **dymax, float **dzmin, float **dzmax, MPI_Request **dmreq,
			   float **sxmin, float **sxmax, float **symin, float **symax, float **szmin, float **szmax,
			   sendCfg **sendBuf, recvCfg **recvBuf, int **rnum, int **disp,
			   MPIinfo orm[], MPIinfo rep[],
			   const int Nx, const int Ny, const int Nz, MPIcfg_tree *mpi,
			   domainCfg *dom, sampling *sample, const ulong Ntot);
  void  releaseORMtopology(float  *dxmin, float  *dxmax, float  *dymin, float  *dymax, float  *dzmin, float  *dzmax, MPI_Request  *dmreq,
			   float  *sxmin, float  *sxmax, float  *symin, float  *symax, float  *szmin, float  *szmax,
			   sendCfg  *sendBuf, recvCfg  *recvBuf, int  *rnum, int  *disp,
			   MPIinfo orm[], MPIinfo rep[], const int rank);
  //-----------------------------------------------------------------------
#else///EXCHANGE_USING_GPUS
  //-----------------------------------------------------------------------
  void configORBtopology
  (domainCfg **box, sendCfg **sendBuf, recvCfg **recvBuf, MPIcfg_tree *mpi, int *ndim, MPIinfo **orb,
   const ulong Ntot, real *samplingRate, int *sampleNumMax, real **sampleLoc, real **sampleFul,
   int **recvNum, int **recvDsp, real **boxMin, real **boxMax);
  void removeORBtopology
  (domainCfg  *box, sendCfg  *sendBuf, recvCfg  *recvBuf, const int ndim, MPIinfo  *orb,
   real  *sampleLoc, real  *sampleFul, int  *recvNum, int  *recvDsp, real  *boxMin, real  *boxMax);
  //-----------------------------------------------------------------------
  void exchangeParticles
  (const int numOld, const int numMax, int *numNew, iparticle ibody_dev, iparticle src, iparticle dst,
   domainDecomposeKey *key,
   sendCfg *sendBuf, recvCfg *recvBuf,
   const float samplingRate, const int sampleNumMax, const double tloc, real *sampleLoc, real *sampleFul, int *recvNum, int *recvDsp, real *boxMin, real *boxMax,
   const int ndim, MPIinfo *orb,
   domainCfg *domain, MPIcfg_tree mpi, measuredTime *measured
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
   , brentStatus *status, brentMemory *memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
   );
  //-----------------------------------------------------------------------
#endif//EXCHANGE_USING_GPUS
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//EXCHANGE_H
//-------------------------------------------------------------------------
