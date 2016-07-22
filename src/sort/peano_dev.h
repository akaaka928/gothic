/*************************************************************************\
 *                                                                       *
                  last updated on 2016/03/08(Tue) 09:45:41
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


//-------------------------------------------------------------------------
#ifndef _SYS_TIME_H
#      include <sys/time.h>
#endif//_SYS_TIME_H
//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef CUDALIB_H
#       include <cudalib.h>
#endif//CUDALIB_H
//-------------------------------------------------------------------------
#ifndef BENCHMARK_H
#       include "../misc/benchmark.h"
#endif//BENCHMARK_H
//-------------------------------------------------------------------------
#ifndef STRUCTURE_H
#       include "../misc/structure.h"
#endif//STRUCTURE_H
//-------------------------------------------------------------------------
#ifndef PEANO_H
#       include "../sort/peano.h"
#endif//PEANO_H
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
  real4 *min, *max;
  int *gsync0, *gsync1;
#ifdef  CUB_AVAILABLE
  void *temp_storage;
  size_t temp_storage_size;
#endif//CUB_AVAILABLE
} soaPHsort;
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
/*   muse allocPeanoHilbertKey_dev */
/*   (const int num, int **idx_dev, PHint **key_dev, PHint **key_hst, real4 **minall, real4 **maxall, int **gsync0, int **gsync1, const deviceProp devProp */
/* #ifndef CALC_MULTIPOLE_ON_DEVICE */
/*    , PHinfo **info_hst */
/* #endif//CALC_MULTIPOLE_ON_DEVICE */
/*    ); */
  muse allocPeanoHilbertKey_dev
  (const int num, int **idx_dev, PHint **key_dev, PHint **key_hst, real4 **minall, real4 **maxall, int **gsync0, int **gsync1,
#ifndef CALC_MULTIPOLE_ON_DEVICE
   PHinfo **info_hst,
#endif//CALC_MULTIPOLE_ON_DEVICE
   soaPHsort *dev, soaPHsort *hst,
#ifdef  CUB_AVAILABLE
   soaPHsort *pre, void **temp_storage, int **idx_pre, PHint **key_pre,
#endif//CUB_AVAILABLE
   const deviceProp devProp);
  void  freePeanoHilbertKey_dev
  (int  *idx_dev, PHint  *key_dev, PHint  *key_hst, real4  *minall, real4  *maxall, int  *gsync0, int  *gsync1
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
