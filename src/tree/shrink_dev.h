/**
 * @file shrink_dev.h
 *
 * @brief Header file to split i-particle groups on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/08/29 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef SHRINK_DEV_H
#define SHRINK_DEV_H


#include "macro.h"
#include "cudalib.h"

#include "../tree/make.h"
#include "../tree/make_dev.h"
#include "../tree/neighbor_dev.h"
#include "../misc/brent.h"


/**
 * @def NUM_IGROUP_SAFETY_FACTOR
 *
 * @brief a parameter to guess maximum number of particle groups
 */
/* #define NUM_IGROUP_SAFETY_FACTOR (8) */
#define NUM_IGROUP_SAFETY_FACTOR (16)


/**
 * @def NTHREADS_SHRINK
 *
 * @brief number of threads per block for countContinuousNeighbor_kernel
 */
#ifndef NTHREADS_SHRINK
#   if  (GPUGEN >= 52)
#define NTHREADS_SHRINK (128)
#else///(GPUGEN >= 52)
#   if  (GPUGEN >= 30)
#define NTHREADS_SHRINK (1024)
#else///(GPUGEN >= 30)
#define NTHREADS_SHRINK (128)
#endif//(GPUGEN >= 30)
#endif//(GPUGEN >= 52)
#endif//NTHREADS_SHRINK


/**
 * @def NEIGHBOR_LENGTH_SHRINK_FACTOR
 *
 * @brief a key parameter to set minimum length in Brent's method
 */
#ifndef NEIGHBOR_LENGTH_SHRINK_FACTOR
#define NEIGHBOR_LENGTH_SHRINK_FACTOR (0.8f)
#endif//NEIGHBOR_LENGTH_SHRINK_FACTOR


/* list of functions appeared in ``shrink_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  void freeParticleGroups
  (laneinfo  *laneInfo_hst, laneinfo  *laneInfo_dev, double  *laneTime_dev, int  *inum_hst, int  *inum_dev);
  muse allocParticleGroups
  (laneinfo **laneInfo_hst, laneinfo **laneInfo_dev, double **laneTime_dev, int **inum_hst, int **inum_dev,
   int *inumPerLane, int *maxNgrp, const int num_max);

  void examineParticleSeparation(const int Ni, iparticle body_dev, brentStatus *brent
#ifdef  EXEC_BENCHMARK
				 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
				 );

  void updateParticleGroups
  (const int Ni, laneinfo *laneInfo, const int inumPerLane, const int maxNgrp, int *Ngrp,
   const iparticle body_dev, int *inum_dev, int *inum_hst, const real rmax
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
   );

  void commitParticleGroups(const int Ngrp, laneinfo *laneInfo_hst, laneinfo *laneInfo_dev);
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//SHRINK_DEV_H
