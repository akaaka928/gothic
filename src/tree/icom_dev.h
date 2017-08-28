/**
 * @file icom_dev.h
 *
 * @brief Header file for generating Enclosing Ball containing all N-body particles on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/08/28 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef ICOM_DEV_H
#define ICOM_DEV_H


/**
 * @def OCTREE_BASED_SEARCH
 *
 * @brief activates octree based traverse while finding the farthest particle from the specified point
 */
#define OCTREE_BASED_SEARCH


/**
 * @def ADOPT_EBS_FOR_LET
 *
 * @brief adopt efficient bounding sphere instead of geometric approximation for pseudo i-particle in LET generation
 */
#define ADOPT_EBS_FOR_LET


/**
 * @def NTHREADS_EB
 *
 * @brief number of threads per block for enclosing ball generator
 */
/**# ここの default 値は性能測定の結果決定する */
/**# また，オプションごとに最適値が変わるはずなのでそこも考慮すること */
#ifndef NTHREADS_EB
#define NTHREADS_EB (256)
#endif//NTHREADS_EB


#define USE_WARP_SHUFFLE_FUNC_EB
#   if  defined(USE_WARP_SHUFFLE_FUNC_EB) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_EB
#endif//defined(USE_WARP_SHUFFLE_FUNC_EB) && (GPUGEN < 30)


#ifdef  OCTREE_BASED_SEARCH
#define NI_R2MAX_ESTIMATE (1024)
#define NBUF_EB (((NI_R2MAX_ESTIMATE) + (NTHREADS_EB) - 1) / (NTHREADS_EB))
/** maximum number of NBUF_EB is 4 to use float4 in union */
#   if  NBUF_EB > 4
#undef  NBUF_EB
#define NBUF_EB (4)
#undef  NI_R2MAX_ESTIMATE
#define NI_R2MAX_ESTIMATE ((NTHREADS_EB) * (NBUF_EB))
#endif//NBUF_EB > 4
#endif//OCTREE_BASED_SEARCH


#include "../misc/benchmark.h"
#include "../misc/structure.h"

#ifdef  OCTREE_BASED_SEARCH
#include "../tree/make.h"
#include "../tree/make_dev.h"
#else///OCTREE_BASED_SEARCH
#include "cudalib.h"
#include "../sort/peano_dev.h"
#endif//OCTREE_BASED_SEARCH

/* list of functions appeared in ``icom_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
#ifndef OCTREE_BASED_SEARCH
  muse allocApproxEnclosingBall_dev(void **dev, const deviceProp devProp);
  void  freeApproxEnclosingBall_dev(void  *dev);
#endif//OCTREE_BASED_SEARCH

  void getApproxEnclosingBall_dev
  (const int num, const iparticle body
#ifdef  OCTREE_BASED_SEARCH
   , const soaTreeCell cell, const soaTreeNode node, const soaMakeTreeBuf buf
#else///OCTREE_BASED_SEARCH
   , void *gmem, const soaPHsort soa, const deviceProp devProp
#endif//OCTREE_BASED_SEARCH
   , const cudaStream_t stream
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
   );
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//ICOM_DEV_H
