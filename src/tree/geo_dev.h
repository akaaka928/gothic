/**
 * @file geo_dev.h
 *
 * @brief Header file for generating enclosing ball
 * whose center is the geometric on of the enclosing rectangular cuboid
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/10/26 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef GEO_DEV_H
#define GEO_DEV_H


#include <sys/time.h>

#include "macro.h"
#include "cudalib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "../tree/make.h"


/**
 * @struct soaGEO
 *
 * @brief structure for geometric enclosing ball
 */
typedef struct
{
  real *r2;
  int Nblock;
} soaGEO;


/* list of functions appeared in ``geo_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  void calc_r2max_dev(const int Ngrp, laneinfo * RESTRICT laneInfo, iparticle *pi, soaGEO dev
#ifdef  EXEC_BENCHMARK
		      , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		      );

  muse allocGeometricEnclosingBall_dev(real **r2_dev, soaGEO *dev, const int num_max);
  void  freeGeometricEnclosingBall_dev(real  *r2_dev);
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//GEO_DEV_H
