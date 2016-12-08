/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:36:43
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
#include "cudalib.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX */
//-------------------------------------------------------------------------
/* #ifdef  ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX */
/* #define  */
/* #endif//ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX */
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
/* in L1 cache preferred configuration, capacity of shared memory is 16KiB per SM */
/* real4 smem[NTHREADS_BOX] corresponds 16 * NTHREADS_BOX bytes */
#define NBLOCKS_PER_SM_BOX (1024 / NTHREADS_BOX)
//-------------------------------------------------------------------------
#define REGISTERS_PER_THREAD_BOX (30)
//-------------------------------------------------------------------------
/* limitation from number of registers */
#   if  NBLOCKS_PER_SM_BOX > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_BOX * NTHREADS_BOX))
#undef  NBLOCKS_PER_SM_BOX
#define NBLOCKS_PER_SM_BOX   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_BOX * NTHREADS_BOX))
#endif//NBLOCKS_PER_SM_BOX > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_BOX * NTHREADS_BOX))
//-------------------------------------------------------------------------
/* maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_BOX > (MAX_THREADS_PER_SM / NTHREADS_BOX)
#undef  NBLOCKS_PER_SM_BOX
#define NBLOCKS_PER_SM_BOX   (MAX_THREADS_PER_SM / NTHREADS_BOX)
#endif//NBLOCKS_PER_SM_BOX > (MAX_THREADS_PER_SM / NTHREADS_BOX)
//-------------------------------------------------------------------------
/* maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_BOX > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_BOX
#define NBLOCKS_PER_SM_BOX   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_BOX > MAX_BLOCKS_PER_SM
//-------------------------------------------------------------------------
/* maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_BOX > ((MAX_WARPS_PER_SM * 32) / NTHREADS_BOX)
#undef  NBLOCKS_PER_SM_BOX
#define NBLOCKS_PER_SM_BOX   ((MAX_WARPS_PER_SM * 32) / NTHREADS_BOX)
#endif//NBLOCKS_PER_SM_BOX > ((MAX_WARPS_PER_SM * 32) / NTHREADS_BOX)
//-------------------------------------------------------------------------
/* # of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_BOX < 1
#undef  NBLOCKS_PER_SM_BOX
#define NBLOCKS_PER_SM_BOX  (1)
#endif//NBLOCKS_PER_SM_BOX < 1
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
#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
//-------------------------------------------------------------------------
//-- List of functions appeared in "exchange_dev.cu"
//-------------------------------------------------------------------------
#include "macro.h"
#include "../misc/structure.h"
#include "../sort/peano_dev.h"
#include "../para/exchange.h"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  void checkBoxSize_dev(const deviceProp devProp);
  /* muse allocateBoxSize_dev(float4 **min_hst, float4 **max_hst, float4 **min_dev, float4 **max_dev, int **gsync0, int **gsync1, soaBoxSize *soa, const deviceProp devProp); */
  /* void  releaseBoxSize_dev(float4  *min_hst, float4  *max_hst, float4  *min_dev, float4  *max_dev, int  *gsync0, int  *gsync1); */
  //-----------------------------------------------------------------------
  void getBoxSize_dev(const int num, position * RESTRICT ipos, soaPHsort soa, const deviceProp devProp, cudaStream_t stream);
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
   soaPHsort soa, const deviceProp devProp, const deviceInfo devInfo,
   double *exchangeInterval, measuredTime *measured
#ifdef  ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
   , const float epsinv
#endif//ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
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
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//EXCHANGE_DEV_H
//-------------------------------------------------------------------------
