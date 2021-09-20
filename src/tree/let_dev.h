/**
 * @file let_dev.h
 *
 * @brief Header file for building locally essential tree (LET)
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/11/30 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef LET_DEV_H
#define LET_DEV_H


#include <mpi.h>

#include "macro.h"
#include "cudalib.h"

#include "../sort/peano.h"
#include "../tree/make.h"
#include "../tree/buf_inc.h"
#include "../tree/let.h"
#include "../tree/walk_dev.h"


#ifdef  USE_ENCLOSING_BALL_FOR_LET
/**
 * @def SKIP_LET_GENERATOR_FOR_NEARBY_NODE
 *
 * @brief disable LET generator (only) for nearby node
 */
#define SKIP_LET_GENERATOR_FOR_NEARBY_NODE
#endif//USE_ENCLOSING_BALL_FOR_LET


#ifdef  SKIP_LET_GENERATOR_FOR_NEARBY_NODE
/**
 * @def THRESHOLD_TO_SKIP_LET_GENERATOR
 *
 * @brief a parameter to discriminate between nearby and distant nodes
 */
#define THRESHOLD_TO_SKIP_LET_GENERATOR (1.25e-1f)
/* #define THRESHOLD_TO_SKIP_LET_GENERATOR (6.25e-2f) */
#endif//SKIP_LET_GENERATOR_FOR_NEARBY_NODE


/**
 * @def NTHREADS_MAKE_LET
 *
 * @brief number of threads per block for makeLET_kernel
 */
#ifndef NTHREADS_MAKE_LET
/* #define NTHREADS_MAKE_LET (128) */
/* #define NTHREADS_MAKE_LET (256) */
/* #define NTHREADS_MAKE_LET (512) */
#define NTHREADS_MAKE_LET (1024)
#endif//NTHREADS_MAKE_LET

/** maximum value of NTHREADS_MAKE_LET is 1024 to calculate prefix sum in 2 stages (32^2) */
#   if  NTHREADS_MAKE_LET > 1024
#undef  NTHREADS_MAKE_LET
#define NTHREADS_MAKE_LET  (1024)
#endif//NTHREADS_MAKE_LET > 1024

/** maximum value of NTHREADS_MAKE_LET on Fermi generation GPUs is 512 */
#   if  (NTHREADS_MAKE_LET > 512) && (GPUVER < 30)
#undef   NTHREADS_MAKE_LET
#define  NTHREADS_MAKE_LET  (512)
#endif//(NTHREADS_MAKE_LET > 512) && (GPUVER < 30)


#define USE_WARP_SHUFFLE_FUNC_MAKE_LET
#   if  defined(USE_WARP_SHUFFLE_FUNC_MAKE_LET) && (GPUVER < 30)
#undef          USE_WARP_SHUFFLE_FUNC_MAKE_LET
#endif//defined(USE_WARP_SHUFFLE_FUNC_MAKE_LET) && (GPUVER < 30)

#define USE_WARP_REDUCE_FUNC_MAKE_LET
#   if  defined(USE_WARP_REDUCE_FUNC_MAKE_LET) && !defined(ENABLE_WARP_REDUCE_FUNC)
#undef          USE_WARP_REDUCE_FUNC_MAKE_LET
#endif//defined(USE_WARP_REDUCE_FUNC_MAKE_LET) && !defined(ENABLE_WARP_REDUCE_FUNC)


#ifdef  DIV_NTHREADS_MAKE_LET
#undef  DIV_NTHREADS_MAKE_LET
#endif//DIV_NTHREADS_MAKE_LET
#   if  NTHREADS_MAKE_LET == 1024
#define DIV_NTHREADS_MAKE_LET(a) ((a) >> 10)
#endif//NTHREADS_MAKE_LET == 1024
#   if  NTHREADS_MAKE_LET ==  512
#define DIV_NTHREADS_MAKE_LET(a) ((a) >>  9)
#endif//NTHREADS_MAKE_LET ==  512
#   if  NTHREADS_MAKE_LET ==  256
#define DIV_NTHREADS_MAKE_LET(a) ((a) >>  8)
#endif//NTHREADS_MAKE_LET ==  256
#   if  NTHREADS_MAKE_LET ==  128
#define DIV_NTHREADS_MAKE_LET(a) ((a) >>  7)
#endif//NTHREADS_MAKE_LET ==  128
#   if  NTHREADS_MAKE_LET ==   64
#define DIV_NTHREADS_MAKE_LET(a) ((a) >>  6)
#endif//NTHREADS_MAKE_LET ==   64
#   if  NTHREADS_MAKE_LET ==   32
#define DIV_NTHREADS_MAKE_LET(a) ((a) >>  5)
#endif//NTHREADS_MAKE_LET ==   32


/** # of blocks per SM is same with tree walk */
/** uint is 4B */
/** SM size is 64KiB on Pascal GPUs */
/** SM size is 48KiB or 16KiB on Fermi or Kepler GPUs */
/** # of elements in SM is 16K on newer GPUs */
/** # of elements in SM is 12K or 4K on older GPUs */
/** SM usage is NTHREADS_MAKE_LET * (NQUEUE_LET + 1) */
#ifndef NQUEUE_LET
#   if  SMPREF_LET == 1
#   if  NBLOCKS_PER_SM == 2
#define NQUEUE_LET (DIV_NTHREADS_MAKE_LET(SMEM_SIZE_SM_PREF >> 3) - 1)
#else///NBLOCKS_PER_SM == 2
#define NQUEUE_LET (DIV_NTHREADS_MAKE_LET((SMEM_SIZE_SM_PREF >> 2) / NBLOCKS_PER_SM) - 1)
#endif//NBLOCKS_PER_SM == 2
#else///SMPREF_LET == 1
#   if  NBLOCKS_PER_SM == 2
#define NQUEUE_LET (DIV_NTHREADS_MAKE_LET(SMEM_SIZE_L1_PREF >> 3) - 1)
#else///NBLOCKS_PER_SM == 2
#define NQUEUE_LET (DIV_NTHREADS_MAKE_LET((SMEM_SIZE_L1_PREF >> 2) / NBLOCKS_PER_SM) - 1)
#endif//NBLOCKS_PER_SM == 2
#endif//SMPREF_LET == 1
/* #   if  GPUVER >= 60 */
/* #define NQUEUE_LET (DIV_NTHREADS_MAKE_LET(16384 / NBLOCKS_PER_SM) - 1) */
/* #else///GPUVER >= 60 */
/* #   if  SMPREF_LET == 1 */
/* #define NQUEUE_LET (DIV_NTHREADS_MAKE_LET(12288 / NBLOCKS_PER_SM) - 1) */
/* #else///SMPREF_LET == 1 */
/* #define NQUEUE_LET (DIV_NTHREADS_MAKE_LET( 4096 / NBLOCKS_PER_SM) - 1) */
/* #endif//SMPREF_LET == 1 */
/* #endif//GPUVER >= 60 */
#endif//NQUEUE_LET


/* list of functions appeared in ``let_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  muse configLETtopology
  (domainInfo **info,
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
   position **min, position **max,
#endif//USE_RECTANGULAR_BOX_FOR_LET
   position **ipos,
#ifdef  GADGET_MAC
   real **amin,
#endif//GADGET_MAC
   int **numSend_hst, int **numSend_dev, cudaStream_t **stream, int *Nstream, const deviceProp gpu, MPIcfg_tree mpi);
  void releaseLETtopology
  (domainInfo  *info,
#ifdef  USE_RECTANGULAR_BOX_FOR_LET
   position  *min, position  *max,
#endif//USE_RECTANGULAR_BOX_FOR_LET
   position  *ipos,
#ifdef  GADGET_MAC
   real  *amin,
#endif//GADGET_MAC
   int  *numSend_hst, int  *numSend_dev, cudaStream_t  *stream, int  Nstream);

  void callGenLET
  (const cudaStream_t stream, domainInfo *let, const soaTreeNode tree, const soaTreeWalkBuf buf
#ifdef  SKIP_LET_GENERATOR_FOR_NEARBY_NODE
   , const position src
#endif//SKIP_LET_GENERATOR_FOR_NEARBY_NODE
#ifdef  MONITOR_LETGEN_TIME
   , unsigned long long int * RESTRICT cycles
#endif//MONITOR_LETGEN_TIME
   );

  void setGlobalConstants_let_dev_cu
  (
#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
   const real theta2_hst
#else///!defined(GADGET_MAC) && !defined(WS93_MAC)
   void
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
   );
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//LET_DEV_H
