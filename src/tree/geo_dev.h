/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/15(Fri) 13:48:03
 *                                                                       *
 *    Header File for generating enclosing ball as GEO                   *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef GEO_DEV_H
#define GEO_DEV_H
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
#ifndef MAKE_H
#       include "../tree/make.h"
#endif//MAKE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
  real *r2;
  int Nblock;
/* #ifdef  CUB_AVAILABLE */
/*   void *temp_storage; */
/*   size_t temp_storage_size; */
/* #endif//CUB_AVAILABLE */
} soaGEO;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "geo_dev.cu"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  void calc_r2max_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, iparticle *pi, soaGEO dev
#ifdef  EXEC_BENCHMARK
		      , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		      );
  //-----------------------------------------------------------------------
  muse allocGeometricEnclosingBall_dev
  (real **r2_dev
/* #ifdef  CUB_AVAILABLE */
/*    , void **temp_storage */
/* #endif//CUB_AVAILABLE */
   , soaGEO *dev, const int num_max);
  void  freeGeometricEnclosingBall_dev
  (real  *r2_dev
/* #ifdef  CUB_AVAILABLE */
/*    , void  *temp_storage */
/* #endif//CUB_AVAILABLE */
   );
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//GEO_DEV_H
//-------------------------------------------------------------------------
