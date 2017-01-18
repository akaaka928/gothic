/*************************************************************************\
 *                                                                       *
                  last updated on 2017/01/18(Wed) 11:12:04
 *                                                                       *
 *    Header File for constructing octree structure                      *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef PEANO_DEV_H
#define PEANO_DEV_H
//-------------------------------------------------------------------------
#include <sys/time.h>
//-------------------------------------------------------------------------
#include "macro.h"
#include "cudalib.h"
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "../sort/peano.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef NTHREADS_PH
#          if  (GPUGEN >= 30)
#define NTHREADS_PH (1024)
#       else///(GPUGEN >= 30)
#define NTHREADS_PH (256)
#       endif//(GPUGEN >= 30)
#endif//NTHREADS_PH
//-------------------------------------------------------------------------
/* tentative treatment for GTX TITAN X: (# of SM = 24) * (# of blocks per SM = 8) = 192 exceeds NTHREADS_PH */
#   if  GPUGEN == 52
#   if  NTHREADS_PH < 256
#undef  NTHREADS_PH
#define NTHREADS_PH  (256)
#endif//NTHREADS_PH < 256
#endif//GPUGEN == 52
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef NTHREADS_PHSORT
#          if  (GPUGEN >= 52)
#define NTHREADS_PHSORT (1024)
#       else///(GPUGEN >= 52)
#          if  (GPUGEN >= 30)
#define NTHREADS_PHSORT (256)
#       else///(GPUGEN >= 30)
#define NTHREADS_PHSORT (1024)
#       endif//(GPUGEN >= 30)
#       endif//(GPUGEN >= 52)
#endif//NTHREADS_PHSORT
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
//-------------------------------------------------------------------------
/* #define TIME_BASED_MODIFICATION */
//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  GENERATE_PHKEY_ON_DEVICE
//-------------------------------------------------------------------------
//-- List of functions appeared in "peano_dev.cu"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  muse allocPeanoHilbertKey_dev
  (const int num, int **idx_dev, PHint **key_dev, PHint **key_hst, float4 **minall, float4 **maxall, int **gsync0, int **gsync1,
#ifndef SERIALIZED_EXECUTION
   float4 **box_min, float4 **box_max, float4 **min_hst, float4 **max_hst,
#endif//SERIALIZED_EXECUTION
#ifndef CALC_MULTIPOLE_ON_DEVICE
   PHinfo **info_hst,
#endif//CALC_MULTIPOLE_ON_DEVICE
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
#ifndef CALC_MULTIPOLE_ON_DEVICE
   , PHinfo  *info_hst
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef  CUB_AVAILABLE
   , void  *temp_storage, int  *idx_pre, PHint  *key_pre
#endif//CUB_AVAILABLE
   );
  //-----------------------------------------------------------------------
  void sortParticlesPHcurve_dev(const int num, iparticle * RESTRICT src, iparticle * RESTRICT dst
				, soaPHsort dev, const deviceProp devProp
#ifdef  CUB_AVAILABLE
				, soaPHsort pre
#endif//CUB_AVAILABLE
#ifndef MAKE_TREE_ON_DEVICE
				, const soaPHsort hst
#endif//MAKE_TREE_ON_DEVICE
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
				, domainLocation *location
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
				, struct timeval *start
#ifdef  EXEC_BENCHMARK
				, wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
				);
  //-----------------------------------------------------------------------
#   if  defined(BLOCK_TIME_STEP) && defined(LOCALIZE_I_PARTICLES) && defined(BRUTE_FORCE_LOCALIZATION)
  void copySortedParticles_dev(const int Ni, const iparticle pi);
#endif//defined(BLOCK_TIME_STEP) && defined(LOCALIZE_I_PARTICLES) && defined(BRUTE_FORCE_LOCALIZATION)
  //-----------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
  void initPHinfo_dev(PHinfo *info);
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------
#endif//GENERATE_PHKEY_ON_DEVICE
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//PEANO_DEV_H
//-------------------------------------------------------------------------
