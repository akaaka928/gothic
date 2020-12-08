/**
 * @file neighbor_dev.h
 *
 * @brief Header file for neighbor search using breadth-first octree on GPU
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
#ifndef NEIGHBOR_DEV_H
#define NEIGHBOR_DEV_H


#include <stdbool.h>

#include "macro.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "../tree/make.h"
#include "../tree/walk_dev.h"


#ifndef HUNT_FIND_PARAMETER
/* the below macro is enabled in the default option; switched off in the parameter survey mode to use -D and -U from Makefile */
#define SMEM_PREF_FOR_NEIGHBOR_SEARCH
#endif//HUNT_FIND_PARAMETER


#define NEIGHBOR_NUM (TSUB / NWARP)
#define NEIGHBOR_NUM_INC (NEIGHBOR_NUM + 1)


/**
 * @def NTHREADS_FACILE_NS
 *
 * @brief number of threads per block for facileNeighborSearching_kernel
 */
#ifndef NTHREADS_FACILE_NS
#   if  GPUVER >= 80
#define NTHREADS_FACILE_NS (256)
#else///GPUVER >= 80
#   if  GPUVER >= 60
#define NTHREADS_FACILE_NS (512)
#else///GPUVER >= 60
#define NTHREADS_FACILE_NS (256)
#endif//GPUVER >= 60
#endif//GPUVER >= 80
#endif//NTHREADS_FACILE_NS


/* list of functions appeared in ``neighbor_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  void facileNeighborSearching_dev(const int Ni, const iparticle pi);

#ifdef  SMEM_PREF_FOR_NEIGHBOR_SEARCH
  void setGlobalConstants_neighbor_dev_cu(void);
#endif//SMEM_PREF_FOR_NEIGHBOR_SEARCH
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//NEIGHBOR_DEV_H
