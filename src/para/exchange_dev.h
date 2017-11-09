/**
 * @file exchange_dev.h
 *
 * @brief Header file for domain decomposition using GPUs with MPI
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/11/09 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef EXCHANGE_DEV_H
#define EXCHANGE_DEV_H


#include "cudalib.h"


#define SHARE_PH_BOX_BOUNDARY

/* #define ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX */

#   if  !defined(SHARE_PH_BOX_BOUNDARY) && defined(ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX)
#define SHARE_PH_BOX_BOUNDARY
#endif//!defined(SHARE_PH_BOX_BOUNDARY) && defined(ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX)


/**
 * @def NTHREADS_BOX
 *
 * @brief number of threads per block for getBoxSize_kernel
 */
#ifndef NTHREADS_BOX
#   if  (GPUGEN >= 30)
#define NTHREADS_BOX (1024)
#else///(GPUGEN >= 30)
#define NTHREADS_BOX (256)
#endif//(GPUGEN >= 30)
#endif//NTHREADS_BOX


/** tentative treatment for GTX TITAN X: (# of SM = 24) * (# of blocks per SM = 8) = 192 exceeds NTHREADS_BOX */
#   if  GPUGEN == 52
#   if  NTHREADS_BOX < 256
#undef  NTHREADS_BOX
#define NTHREADS_BOX  (256)
#endif//NTHREADS_BOX < 256
#endif//GPUGEN == 52

#   if  GPUGEN >= 60
/** capacity of shared memory is 64KiB per SM on newer GPUs */
/** real4 smem[NTHREADS_BOX] corresponds 16 * NTHREADS_BOX bytes */
#define NBLOCKS_PER_SM_BOX (4096 / NTHREADS_BOX)
#else///GPUGEN >= 60
/** in L1 cache preferred configuration, capacity of shared memory is 16KiB per SM on older GPUs */
/** real4 smem[NTHREADS_BOX] corresponds 16 * NTHREADS_BOX bytes */
#define NBLOCKS_PER_SM_BOX (1024 / NTHREADS_BOX)
#endif//GPUGEN >= 60

#define REGISTERS_PER_THREAD_BOX (32)

/** limitation from number of registers */
#   if  NBLOCKS_PER_SM_BOX > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_BOX * NTHREADS_BOX))
#undef  NBLOCKS_PER_SM_BOX
#define NBLOCKS_PER_SM_BOX   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_BOX * NTHREADS_BOX))
#endif//NBLOCKS_PER_SM_BOX > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_BOX * NTHREADS_BOX))

/** maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_BOX > (MAX_THREADS_PER_SM / NTHREADS_BOX)
#undef  NBLOCKS_PER_SM_BOX
#define NBLOCKS_PER_SM_BOX   (MAX_THREADS_PER_SM / NTHREADS_BOX)
#endif//NBLOCKS_PER_SM_BOX > (MAX_THREADS_PER_SM / NTHREADS_BOX)

/** maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_BOX > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_BOX
#define NBLOCKS_PER_SM_BOX   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_BOX > MAX_BLOCKS_PER_SM

/** maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_BOX > ((MAX_WARPS_PER_SM * 32) / NTHREADS_BOX)
#undef  NBLOCKS_PER_SM_BOX
#define NBLOCKS_PER_SM_BOX   ((MAX_WARPS_PER_SM * 32) / NTHREADS_BOX)
#endif//NBLOCKS_PER_SM_BOX > ((MAX_WARPS_PER_SM * 32) / NTHREADS_BOX)

/** # of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_BOX < 1
#undef  NBLOCKS_PER_SM_BOX
#define NBLOCKS_PER_SM_BOX  (1)
#endif//NBLOCKS_PER_SM_BOX < 1


/**
 * @def NTHREADS_SETIDX
 *
 * @brief number of threads per block for setIndex_kernel
 */
#ifndef NTHREADS_SETIDX
#define NTHREADS_SETIDX (1024)
#endif//NTHREADS_SETIDX

/**
 * @def NTHREADS_PICKUP
 *
 * @brief number of threads per block for pickupSamples_kernel
 */
#ifndef NTHREADS_PICKUP
#define NTHREADS_PICKUP (256)
#endif//NTHREADS_PICKUP

/**
 * @def NTHREADS_SORTYZ
 *
 * @brief number of threads per block for sortSamplePos_yz_kernel
 */
#ifndef NTHREADS_SORTYZ
#define NTHREADS_SORTYZ (1024)
#endif//NTHREADS_SORTYZ

/**
 * @def NTHREADS_SETPOS
 *
 * @brief number of threads per block for setParticlePosition_kernel
 */
#ifndef NTHREADS_SETPOS
#define NTHREADS_SETPOS (1024)
#endif//NTHREADS_SETPOS

/**
 * @def NTHREADS_DDSORT
 *
 * @brief number of threads per block for sortParticlesDDkey_kernel
 */
#ifndef NTHREADS_DDSORT
#define NTHREADS_DDSORT (1024)
#endif//NTHREADS_DDSORT


/**
 * @def NTHREADS_ASSIGN
 *
 * @brief number of threads per block for assignNewDomain_kernel
 */
#ifndef NTHREADS_ASSIGN
#define NTHREADS_ASSIGN (256)
#endif//NTHREADS_ASSIGN


#   if  defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
/* list of functions appeared in ``exchange_dev.cu'' */

#include "macro.h"
#include "../misc/structure.h"
#include "../sort/peano_dev.h"
#include "../para/exchange.h"


/**
 * @struct sendDom
 *
 * @brief structure for sending domain
 */
typedef struct
{
  float *xmin_dev, *xmax_dev, *ymin_dev, *ymax_dev, *zmin_dev, *zmax_dev;
  float *xmin_hst, *xmax_hst, *ymin_hst, *ymax_hst, *zmin_hst, *zmax_hst;
  int *rank, *rank_hst;
  int *numNew, *numNew_hst;
  int4 *gmem;
  int *gsync0, *gsync1;
} sendDom;


#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  void checkBoxSize_dev(const deviceProp devProp);

  void getBoxSize_dev(const int num, position * RESTRICT ipos, soaPHsort soa, const deviceProp devProp, cudaStream_t stream);

  muse allocateSamplePos
  (float **x0hst, float **x1hst, float **y0hst, float **y1hst, float **z0hst, float **z1hst, int **idhst,
   float **x0dev, float **x1dev, float **y0dev, float **y1dev, float **z0dev, float **z1dev, int **iddev,
   samplePos * RESTRICT pos0, samplePos * RESTRICT pos1, const sampling sample);
  void  releaseSamplePos
  (float  *x0hst, float  *x1hst, float  *y0hst, float  *y1hst, float  *z0hst, float  *z1hst, int  *idhst,
   float  *x0dev, float  *x1dev, float  *y0dev, float  *y1dev, float  *z0dev, float  *z1dev, int  *iddev);

  muse allocateParticlePosition(float **xhst, float **yhst, float **zhst, particlePos *hst,
				float **xdev, float **ydev, float **zdev, particlePos *dev, int **rank_hst, int **rank_dev, int **idx_dev, domainDecomposeKey *key, const ulong Ntot);
  void  releaseParticlePosition(float  *xhst, float  *yhst, float  *zhst,
				float  *xdev, float  *ydev, float  *zdev,	            int  *rank_hst, int  *rank_dev, int  *idx_dev);

  muse allocateDomainPos(float **xmin_dev, float **xmax_dev, float **ymin_dev, float **ymax_dev, float **zmin_dev, float **zmax_dev,
			 float **xmin_hst, float **xmax_hst, float **ymin_hst, float **ymax_hst, float **zmin_hst, float **zmax_hst,
			 int **rank, int **rank_hst, int **numNew, int **numNew_hst, int4 **gmem, int **gsync0, int **gsync1,
			 sendDom *dom, const int Ngpu, const deviceProp devProp);
  void  releaseDomainPos(float  *xmin_dev, float  *xmax_dev, float  *ymin_dev, float  *ymax_dev, float  *zmin_dev, float  *zmax_dev,
			 float  *xmin_hst, float  *xmax_hst, float  *ymin_hst, float  *ymax_hst, float  *zmin_hst, float  *zmax_hst,
			 int  *rank, int  *rank_hst, int  *numNew, int  *numNew_hst, int4  *gmem, int  *gsync0, int  *gsync1);

  void exchangeParticles_dev
  (const int numOld, const ulong Ntot, const int numMax, int *numNew,
   iparticle * RESTRICT src_dev, iparticle * RESTRICT dst_dev, iparticle * RESTRICT src_hst, iparticle * RESTRICT dst_hst,
   sendDom domBoundary, particlePos pos_hst, particlePos pos_dev, domainDecomposeKey key, sendCfg *sendBuf, recvCfg *recvBuf,
   MPIinfo orm[], MPIinfo rep[], domainCfg domain, MPIcfg_tree mpi,
   const double tloc, sampling sample, samplePos loc, samplePos ful,
   soaPHsort soa, const deviceProp devProp, const deviceInfo devInfo,
   measuredTime *measured
#ifdef  CARE_EXTERNAL_PARTICLES
   , domainLocation *location
#endif//CARE_EXTERNAL_PARTICLES
   , brentStatus *status, brentMemory *memory
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
   );
#ifdef  __CUDACC__
}
#endif//__CUDACC__
#endif//defined(MPI_INCLUDED) || defined(OMPI_MPI_H)


#endif//EXCHANGE_DEV_H
