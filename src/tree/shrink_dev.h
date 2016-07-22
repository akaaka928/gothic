/*************************************************************************\
 *                                                                       *
                  last updated on 2016/03/08(Tue) 09:43:05
 *                                                                       *
 *    Header File to examine tree statistics                             *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef SHRINK_DEV_H
#define SHRINK_DEV_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef CUDALIB_H
#       include <cudalib.h>
#endif//CUDALIB_H
//-------------------------------------------------------------------------
#ifndef MAKE_H
#       include "../tree/make.h"
#endif//MAKE_H
//-------------------------------------------------------------------------
#ifndef MAKE_DEV_H
#       include "../tree/make_dev.h"
#endif//MAKE_DEV_H
//-------------------------------------------------------------------------
#   if  !defined(NEIGHBOR_DEV_H) && defined(LOCALIZE_I_PARTICLES)
#       include "../tree/neighbor_dev.h"
#endif//!defined(NEIGHBOR_DEV_H) && defined(LOCALIZE_I_PARTICLES)
//-------------------------------------------------------------------------
#   if  defined(USE_BRENT_METHOD) && !defined(BRENT_H)
#       include "../misc/brent.h"
#endif//defined(USE_BRENT_METHOD) && !defined(BRENT_H)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define NUM_IGROUP_SAFETY_FACTOR (8) */
#define NUM_IGROUP_SAFETY_FACTOR (16)
//-------------------------------------------------------------------------
#ifndef NTHREADS_SHRINK
#          if  (GPUGEN >= 52)
#define NTHREADS_SHRINK (128)
#       else///(GPUGEN >= 52)
#          if  (GPUGEN >= 30)
#define NTHREADS_SHRINK (1024)
#       else///(GPUGEN >= 30)
#define NTHREADS_SHRINK (128)
#       endif//(GPUGEN >= 30)
#       endif//(GPUGEN >= 52)
#endif//NTHREADS_SHRINK
//-------------------------------------------------------------------------
/* #ifndef SHRINK_FRACTION */
/* /\* 1/2^2, 1/4^2, 1/8^2, 1/16^2, 1/32^2 *\/ */
/* /\* #define SHRINK_FRACTION (2.5e-1f) *\/ */
/* /\* #define SHRINK_FRACTION (6.25e-2f) *\/ */
/* #define SHRINK_FRACTION (1.5625e-2f) */
/* /\* #define SHRINK_FRACTION (3.90625e-3f) *\/ */
/* /\* #define SHRINK_FRACTION (9.765625e-4f) *\/ */
/* #endif//SHRINK_FRACTION */
//-------------------------------------------------------------------------
#ifdef  USE_BRENT_METHOD
//-------------------------------------------------------------------------
#ifndef NEIGHBOR_LENGTH_SHRINK_FACTOR
#define NEIGHBOR_LENGTH_SHRINK_FACTOR (0.8f)
#endif//NEIGHBOR_LENGTH_SHRINK_FACTOR
//-------------------------------------------------------------------------
#else///USE_BRENT_METHOD
//-------------------------------------------------------------------------
#ifndef NEIGHBOR_LENGTH_UPPER_FRACTION
/* #define NEIGHBOR_LENGTH_UPPER_FRACTION (1.0e-4f) */
/* #define NEIGHBOR_LENGTH_UPPER_FRACTION (5.0e-4f) */
#define NEIGHBOR_LENGTH_UPPER_FRACTION (1.0e-3f)
/* #define NEIGHBOR_LENGTH_UPPER_FRACTION (5.0e-3f) */
/* #define NEIGHBOR_LENGTH_UPPER_FRACTION (1.0e-2f) */
#endif//NEIGHBOR_LENGTH_UPPER_FRACTION
//-------------------------------------------------------------------------
#ifndef NEIGHBOR_LENGTH_EXTEND_FRACTION
#define NEIGHBOR_LENGTH_EXTEND_FRACTION (0.1f)
/* #define NEIGHBOR_LENGTH_EXTEND_FRACTION (0.2f) */
/* #define NEIGHBOR_LENGTH_EXTEND_FRACTION (0.25f) */
/* #define NEIGHBOR_LENGTH_EXTEND_FRACTION (0.0625f) */
/* #define NEIGHBOR_LENGTH_EXTEND_FRACTION (0.0f) */
#endif//NEIGHBOR_LENGTH_EXTEND_FRACTION
//-------------------------------------------------------------------------
#endif//USE_BRENT_METHOD
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  CUB_AVAILABLE
typedef struct
{
  void *temp_storage;
  real *out;
  size_t temp_storage_size;
} soaCUBreal;
#endif//CUB_AVAILABLE
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "shrink_dev.cu"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  void freeParticleGroups
  (laneinfo  *laneInfo_hst, laneinfo  *laneInfo_dev, double  *laneTime_dev
#ifdef  LOCALIZE_I_PARTICLES
   , int  *inum_hst, int  *inum_dev
#   if  defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
   , void  *temp_storage, real  *outCub
#endif//defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
#endif//LOCALIZE_I_PARTICLES
   );
  muse allocParticleGroups
  (laneinfo **laneInfo_hst, laneinfo **laneInfo_dev, double **laneTime_dev
#ifdef  LOCALIZE_I_PARTICLES
   , int **inum_hst, int **inum_dev
#   if  defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
   , soaCUBreal *util, void **temp_storage, real **outCub, iparticle body_dev
#endif//defined(CUB_AVAILABLE) && !defined(USE_BRENT_METHOD)
#endif//LOCALIZE_I_PARTICLES
   , int *inumPerLane, int *maxNgrp, const int num_max, deviceProp devProp);
  //-----------------------------------------------------------------------
#ifdef  LOCALIZE_I_PARTICLES
  void examineParticleSeparation
  (const int Ni, iparticle body_dev
#ifdef  USE_BRENT_METHOD
   , brentStatus *brent
#else///USE_BRENT_METHOD
#ifdef  CUB_AVAILABLE
   , soaCUBreal util
#endif//CUB_AVAILABLE
   , real *rmax
#endif//USE_BRENT_METHOD
#ifndef FACILE_NEIGHBOR_SEARCH
   , const soaTreeCell cell, const soaTreeNode node, const soaMakeTreeBuf makeBuf, const soaNeighborSearchBuf searchBuf, deviceProp devProp
#endif//FACILE_NEIGHBOR_SEARCH
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
   );
#endif///LOCALIZE_I_PARTICLES
  //-----------------------------------------------------------------------
  void updateParticleGroups
  (const int Ni, laneinfo *laneInfo, const int inumPerLane, const int maxNgrp, int *Ngrp
#ifdef  LOCALIZE_I_PARTICLES
   , const iparticle body_dev, int *inum_dev, int *inum_hst
#ifdef  USE_BRENT_METHOD
   , const real rmax
#endif//USE_BRENT_METHOD
#   if  !defined(FACILE_NEIGHBOR_SEARCH) && !defined(USE_BRENT_METHOD)
   , const soaTreeCell cell, const soaTreeNode node, const soaMakeTreeBuf makeBuf, const soaNeighborSearchBuf searchBuf, deviceProp devProp
#endif//!defined(FACILE_NEIGHBOR_SEARCH) && !defined(USE_BRENT_METHOD)
#endif//LOCALIZE_I_PARTICLES
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
   );
  //-----------------------------------------------------------------------
  void commitParticleGroups(const int Ngrp, laneinfo *laneInfo_hst, laneinfo *laneInfo_dev);
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//SHRINK_DEV_H
//-------------------------------------------------------------------------
