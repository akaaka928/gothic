/*************************************************************************\
 *                                                                       *
                  last updated on 2016/11/11(Fri) 11:47:18
 *                                                                       *
 *    Header File for domain decomposition using GPUs                    *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef EXCHANGE_DEV_H
#define EXCHANGE_DEV_H
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "../para/exchange.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef NTHREADS_BOX
#          if  (GPUGEN >= 30)
#define NTHREADS_BOX (1024)
#       else///(GPUGEN >= 30)
#define NTHREADS_BOX (256)
#       endif//(GPUGEN >= 30)
#endif//NTHREADS_BOX
//-------------------------------------------------------------------------
/* tentative treatment for GTX TITAN X: (# of SM = 24) * (# of blocks per SM = 8) = 192 exceeds NTHREADS_BOX */
#   if  GPUGEN == 52
#   if  NTHREADS_BOX < 256
#undef  NTHREADS_BOX
#define NTHREADS_BOX  (256)
#endif//NTHREADS_BOX < 256
#endif//GPUGEN == 52
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef NTHREADS_SETIDX
#define NTHREADS_SETIDX (1024)
#endif//NTHREADS_SETIDX
//-------------------------------------------------------------------------
#ifndef NTHREADS_SORTYZ
#define NTHREADS_SORTYZ (1024)
#endif//NTHREADS_SORTYZ
//-------------------------------------------------------------------------
#ifndef NTHREADS_SETPOS
#define NTHREADS_SETPOS (1024)
#endif//NTHREADS_SETPOS
//-------------------------------------------------------------------------
#ifndef NTHREADS_DDSORT
#define NTHREADS_DDSORT (1024)
#endif//NTHREADS_DDSORT
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* must be determined by benchmarks */
//-------------------------------------------------------------------------
#ifndef NCRIT_XPOS_SORT
#define NCRIT_XPOS_SORT (65536)
#endif//NCRIT_XPOS_SORT
//-------------------------------------------------------------------------
#ifndef NCRIT_YPOS_SORT
#define NCRIT_YPOS_SORT (65536)
#endif//NCRIT_YPOS_SORT
//-------------------------------------------------------------------------
#ifndef NCRIT_ZPOS_SORT
#define NCRIT_ZPOS_SORT (65536)
#endif//NCRIT_ZPOS_SORT
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
  float4 *min_dev, *max_dev, *min_hst, *max_hst;
  int *gsync0, *gsync1;
} soaBoxSize;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "exchange_dev.cu"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  muse allocateBoxSize_dev(float4 **min_hst, float4 **max_hst, float4 **min_dev, float4 **max_dev, int **gsync0, int **gsync1, soaBoxSize *soa, const deviceProp devProp);
  void  releaseBoxSize_dev(float4  *min_hst, float4  *max_hst, float4  *min_dev, float4  *max_dev, int  *gsync0, int  *gsync1);
  //-----------------------------------------------------------------------
  void getBoxSize_dev(const int num, position * RESTRICT ipos, soaBoxSize soa, const deviceProp devProp, cudaStream_t stream);
  //-----------------------------------------------------------------------
  muse allocateSamplePos
  (float **x0hst, float **x1hst, float **y0hst, float **y1hst, float **z0hst, float **z1hst, int **idhst,
   float **x0dev, float **x1dev, float **y0dev, float **y1dev, float **z0dev, float **z1dev, int **iddev,
   samplePos * RESTRICT pos0, samplePos * RESTRICT pos1, const sampling sample);
  void  releaseSamplePos
  (float  *x0hst, float  *x1hst, float  *y0hst, float  *y1hst, float  *z0hst, float  *z1hst, int  *idhst,
   float  *x0dev, float  *x1dev, float  *y0dev, float  *y1dev, float  *z0dev, float  *z1dev, int  *iddev);
  //-----------------------------------------------------------------------
  muse allocateParticlePosition(float **xhst, float **yhst, float **zhst, particlePos *hst,
				float **xdev, float **ydev, float **zdev, particlePos *dev, int **rank_hst, int **rank_dev, int **idx_dev, domainDecomposeKey *key, const ulong Ntot);
  void  releaseParticlePosition(float  *xhst, float  *yhst, float  *zhst,
				float  *xdev, float  *ydev, float  *zdev,	            int  *rank_hst, int  *rank_dev, int  *idx_dev);
  //-----------------------------------------------------------------------
  void exchangeParticles_dev
  (const int numOld, const ulong Ntot, const int numMax, int *numNew,
   iparticle * RESTRICT src_dev, iparticle * RESTRICT dst_dev, iparticle * RESTRICT src_hst, iparticle * RESTRICT dst_hst,
   particlePos pos_hst, particlePos pos_dev, domainDecomposeKey key, sendCfg *sendBuf, recvCfg *recvBuf,
   MPIinfo orm[], MPIinfo rep[], domainCfg domain, MPIcfg_tree mpi,
   const double tloc, sampling sample, samplePos loc, samplePos ful,
   soaBoxSize soa, const deviceProp devProp, const deviceInfo devInfo,
   double *exchangeInterval, measuredTime *measured
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
#endif//EXCHANGE_DEV_H
//-------------------------------------------------------------------------
