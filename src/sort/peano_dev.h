/**
 * @file peano_dev.h
 *
 * @brief Header file for calculating Peano--Hilbert space-filling curve
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/01 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef PEANO_DEV_H
#define PEANO_DEV_H


#include <sys/time.h>

#include "macro.h"
#include "cudalib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "../sort/peano.h"


/**
 * @def NTHREADS_PH
 *
 * @brief number of threads per block for calcPHkey_kernel
 */
#ifndef NTHREADS_PH
#   if  (GPUGEN >= 30)
#define NTHREADS_PH (1024)
#else///(GPUGEN >= 30)
#define NTHREADS_PH (256)
#endif//(GPUGEN >= 30)
#endif//NTHREADS_PH


/** tentative treatment for GTX TITAN X: (# of SM = 24) * (# of blocks per SM = 8) = 192 exceeds NTHREADS_PH */
#   if  GPUGEN == 52
#   if  NTHREADS_PH < 256
#undef  NTHREADS_PH
#define NTHREADS_PH  (256)
#endif//NTHREADS_PH < 256
#endif//GPUGEN == 52


/**
 * @def NTHREADS_PHSORT
 *
 * @brief number of threads per block for sortParticlesPHcurve_kernel
 */
#ifndef NTHREADS_PHSORT
#   if  (GPUGEN >= 52)
#define NTHREADS_PHSORT (1024)
#else///(GPUGEN >= 52)
#   if  (GPUGEN >= 30)
#define NTHREADS_PHSORT (256)
#else///(GPUGEN >= 30)
#define NTHREADS_PHSORT (1024)
#endif//(GPUGEN >= 30)
#endif//(GPUGEN >= 52)
#endif//NTHREADS_PHSORT


/**
 * @struct soaPHsort
 *
 * @brief structure for PH-key
 */
typedef struct
{
  PHint *key;
  int *idx;
  float4 *min, *max;
#ifndef SERIALIZED_EXECUTION
  float4 *box_min, *box_max;
  float4 *min_hst, *max_hst;
#endif//SERIALIZED_EXECUTION
  int *gsync0, *gsync1;
#ifdef  CUB_AVAILABLE
  void *temp_storage;
  size_t temp_storage_size;
#endif//CUB_AVAILABLE
} soaPHsort;


#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
/* #define TIME_BASED_MODIFICATION */
/**
 * @struct domainLocation
 *
 * @brief structure for PH-key with domain decomposition
 */
typedef struct
{
  float3 boxmin, boxmax;
  float *diameter_dev, *diameter_hst;
#ifdef  TIME_BASED_MODIFICATION
  float elapsed, dtmin, dtinv;
  float linv;
#else///TIME_BASED_MODIFICATION
  float step;
  float dL_L;
#endif//TIME_BASED_MODIFICATION
  float eps, eta;
} domainLocation;
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)


/* list of functions appeared in ``peano_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  muse allocPeanoHilbertKey_dev
  (const int num, int **idx_dev, PHint **key_dev, PHint **key_hst, float4 **minall, float4 **maxall, int **gsync0, int **gsync1,
#ifndef SERIALIZED_EXECUTION
   float4 **box_min, float4 **box_max, float4 **min_hst, float4 **max_hst,
#endif//SERIALIZED_EXECUTION
   soaPHsort *dev, soaPHsort *hst,
#ifdef  CUB_AVAILABLE
   soaPHsort *pre, void **temp_storage, int **idx_pre, PHint **key_pre,
#endif//CUB_AVAILABLE
   const deviceProp devProp);
  void  freePeanoHilbertKey_dev
  (int  *idx_dev, PHint  *key_dev, PHint  *key_hst, float4  *minall, float4  *maxall, int  *gsync0, int  *gsync1
#ifndef SERIALIZED_EXECUTION
   , float4  *box_min, float4  *box_max, float4  *min_hst, float4  *max_hst
#endif//SERIALIZED_EXECUTION
#ifdef  CUB_AVAILABLE
   , void  *temp_storage, int  *idx_pre, PHint  *key_pre
#endif//CUB_AVAILABLE
   );

  void sortParticlesPHcurve_dev(const int num, iparticle * RESTRICT src, iparticle * RESTRICT dst, soaPHsort dev, const deviceProp devProp
#ifdef  CUB_AVAILABLE
				, soaPHsort pre
#endif//CUB_AVAILABLE
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
				, domainLocation *location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
				, struct timeval *start
#ifdef  EXEC_BENCHMARK
				, wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
				);

#ifdef  BLOCK_TIME_STEP
  void copySortedParticles_dev(const int Ni, const iparticle pi);
#endif//BLOCK_TIME_STEP

  void initPHinfo_dev(PHinfo *info);
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//PEANO_DEV_H
