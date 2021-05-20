/**
 * @file peano_dev.h
 *
 * @brief Header file for calculating Peano--Hilbert space-filling curve
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/12/08 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef PEANO_DEV_H
#define PEANO_DEV_H


#include <time.h>

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
#   if  (GPUVER >= 80)
#define NTHREADS_PH (1024)
#else///(GPUVER >= 80)
#   if  (GPUVER >= 70)
#define NTHREADS_PH (512)
#else///(GPUVER >= 70)
#   if  (GPUVER >= 30)
#define NTHREADS_PH (1024)
#else///(GPUVER >= 30)
#define NTHREADS_PH (256)
#endif//(GPUVER >= 30)
#endif//(GPUVER >= 70)
#endif//(GPUVER >= 80)
#endif//NTHREADS_PH


/** tentative treatment for GTX TITAN X: (# of SM = 24) * (# of blocks per SM = 8) = 192 exceeds NTHREADS_PH */
#   if  GPUVER >= 52
#   if  NTHREADS_PH < 256
#undef  NTHREADS_PH
#define NTHREADS_PH  (256)
#endif//NTHREADS_PH < 256
#endif//GPUVER >= 52


/** tentative treatment for Tesla P100: (# of SM = 56) * (# of blocks per SM = 6) = 336 exceeds NTHREADS_PH */
#   if  GPUVER >= 60
#   if  NTHREADS_PH < 512
#undef  NTHREADS_PH
#define NTHREADS_PH  (512)
#endif//NTHREADS_PH < 512
#endif//GPUVER >= 60


/**
 * @def NTHREADS_PHSORT
 *
 * @brief number of threads per block for sortParticlesPHcurve_kernel
 */
#ifndef NTHREADS_PHSORT
#   if  GPUVER >= 80
#define NTHREADS_PHSORT (256)
#else///GPUVER >= 80
#   if  GPUVER == 52
#define NTHREADS_PHSORT (1024)
#else///GPUVER == 52
#define NTHREADS_PHSORT (256)
#endif//GPUVER == 52
#endif//GPUVER >= 80
#endif//NTHREADS_PHSORT


#ifdef  RESET_CENTER_OF_MASS
#define NTHREADS_RESET_COM (1024)
#endif//RESET_CENTER_OF_MASS


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
#ifdef  RESET_CENTER_OF_MASS
  float *r2;
  double4 *com_all;
  double3 *vel_all;
  bool *converge, *converge_hst;
  int *gsync0_com, *gsync1_com;
#endif//RESET_CENTER_OF_MASS
#ifdef  CUB_AVAILABLE
  void *temp_storage;
  size_t temp_storage_size;
#endif//CUB_AVAILABLE
#ifdef  USE_OCCUPANCY_CALCULATOR
  int numBlocksPerSM_peano;
#ifdef  RESET_CENTER_OF_MASS
  int numBlocksPerSM_reset;
#endif//RESET_CENTER_OF_MASS
#ifndef SERIALIZED_EXECUTION
  int numBlocksPerSM_box;
#endif//SERIALIZED_EXECUTION
#endif//USE_OCCUPANCY_CALCULATOR
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
  void sortParticlesPHcurve_dev(const int num, iparticle * RESTRICT src, iparticle * RESTRICT dst, soaPHsort dev, const deviceProp devProp
#ifdef  CUB_AVAILABLE
				, soaPHsort pre
#endif//CUB_AVAILABLE
#ifdef  SET_SINK_PARTICLES
				, const int Nbh, const sinkparticle bh
#endif//SET_SINK_PARTICLES
#ifdef  RESET_CENTER_OF_MASS
				, const int com_group_head, const int com_group_num
#endif//RESET_CENTER_OF_MASS
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
				, domainLocation *location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
				, struct timespec *start
#ifdef  EXEC_BENCHMARK
				, wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
				);

#ifdef  BLOCK_TIME_STEP
  void copySortedParticles_dev(const int Ni, const iparticle pi);
#endif//BLOCK_TIME_STEP

  void initPHinfo_dev(PHinfo *info);

  muse allocPeanoHilbertKey_dev
  (const int num, int **idx_dev, PHint **key_dev, PHint **key_hst, float4 **minall, float4 **maxall, int **gsync0, int **gsync1,
#ifdef  RESET_CENTER_OF_MASS
   float **r2_dev, double4 **com_all, double3 **vel_all, int **gsync0_com, int **gsync1_com, bool **converge, bool **converge_hst,
#endif//RESET_CENTER_OF_MASS
#ifndef SERIALIZED_EXECUTION
   float4 **box_min, float4 **box_max, float4 **min_hst, float4 **max_hst,
#endif//SERIALIZED_EXECUTION
   soaPHsort *dev, soaPHsort *hst,
#ifdef  CUB_AVAILABLE
   soaPHsort *pre, void **temp_storage, int **idx_pre, PHint **key_pre,
#ifdef  RESET_CENTER_OF_MASS
    float **r2_pre,
#endif//RESET_CENTER_OF_MASS
#endif//CUB_AVAILABLE
   const deviceProp devProp);
  void  freePeanoHilbertKey_dev
  (int  *idx_dev, PHint  *key_dev, PHint  *key_hst, float4  *minall, float4  *maxall, int  *gsync0, int  *gsync1
#ifdef  RESET_CENTER_OF_MASS
   , float  *r2_dev, double4  *com_all, double3  *vel_all, int  *gsync0_com, int  *gsync1_com, bool  *converge, bool  *converge_hst
#endif//RESET_CENTER_OF_MASS
#ifndef SERIALIZED_EXECUTION
   , float4  *box_min, float4  *box_max, float4  *min_hst, float4  *max_hst
#endif//SERIALIZED_EXECUTION
#ifdef  CUB_AVAILABLE
   , void  *temp_storage, int  *idx_pre, PHint  *key_pre
#ifdef  RESET_CENTER_OF_MASS
   , float  *r2_pre
#endif//RESET_CENTER_OF_MASS
#endif//CUB_AVAILABLE
   );
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//PEANO_DEV_H
