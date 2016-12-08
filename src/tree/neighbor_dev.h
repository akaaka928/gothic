/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:56:03
 *                                                                       *
 *    Header File for neighbor searching using breadth-first tree        *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef NEIGHBOR_DEV_H
#define NEIGHBOR_DEV_H
//-------------------------------------------------------------------------
#include <stdbool.h>
//-------------------------------------------------------------------------
#include "macro.h"
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "../tree/make.h"
#include "../tree/walk_dev.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef HUNT_FIND_PARAMETER
/* the below macro is enabled in the default option; switched off in the parameter survey mode to use -D and -U from Makefile */
#define SMEM_PREF_FOR_NEIGHBOR_SEARCH
#endif//HUNT_FIND_PARAMETER
//-------------------------------------------------------------------------
#define NEIGHBOR_NUM (TSUB / NWARP)
#define NEIGHBOR_NUM_INC (NEIGHBOR_NUM + 1)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  FACILE_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------
#ifndef NTHREADS_FACILE_NS
#define NTHREADS_FACILE_NS (256)
#endif//NTHREADS_FACILE_NS
//-------------------------------------------------------------------------
#endif//FACILE_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef FACILE_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------
#define USE_SMID_TO_GET_BUFID_NEIGHBOR
#define TRY_MODE_ABOUT_BUFFER_NEIGHBOR
//-------------------------------------------------------------------------
#ifndef HUNT_FIND_PARAMETER
/* the below macro is enabled in the default option; switched off in the parameter survey mode to use -D and -U from Makefile */
#define USE_WARP_SHUFFLE_FUNC_NEIGHBOR
#endif//HUNT_FIND_PARAMETER
//-------------------------------------------------------------------------
#   if  defined(USE_WARP_SHUFFLE_FUNC_NEIGHBOR) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_NEIGHBOR
#endif//defined(USE_WARP_SHUFFLE_FUNC_NEIGHBOR) && (GPUGEN < 30)
//-------------------------------------------------------------------------
#ifndef NTHREADS_NEIGHBOR
/* #define NTHREADS_NEIGHBOR (512) */
#define NTHREADS_NEIGHBOR (256)
/* #define NTHREADS_NEIGHBOR (128) */
#endif//NTHREADS_NEIGHBOR
/* NTHREADS_NEIGHBOR must be equal or smaller than 512 due to the capacity of shared memory */
#   if  NTHREADS_NEIGHBOR > 512
#undef  NTHREADS_NEIGHBOR
#define NTHREADS_NEIGHBOR  (512)
#endif//NTHREADS_NEIGHBOR > 512
//-------------------------------------------------------------------------
#ifndef TSUB_NEIGHBOR
#define TSUB_NEIGHBOR (32)
#endif//TSUB_NEIGHBOR
/* TSUB_NEIGHBOR must be equal or smaller than NTHREADS_NEIGHBOR */
#   if  TSUB_NEIGHBOR > NTHREADS_NEIGHBOR
#undef  TSUB_NEIGHBOR
#define TSUB_NEIGHBOR   NTHREADS_NEIGHBOR
#endif//TSUB_NEIGHBOR > NTHREADS_NEIGHBOR
//-------------------------------------------------------------------------
#   if  NTHREADS_NEIGHBOR >= 128
/* TSUB_NEIGHBOR must be equal or greater than 4 */
#          if  TSUB_NEIGHBOR < 4
#       undef  TSUB_NEIGHBOR
#       define TSUB_NEIGHBOR  (4)
#       endif//TSUB_NEIGHBOR < 4
#   if  NTHREADS_NEIGHBOR >= 256
/* TSUB_NEIGHBOR must be equal or greater than 8 */
#          if  TSUB_NEIGHBOR < 8
#       undef  TSUB_NEIGHBOR
#       define TSUB_NEIGHBOR  (8)
#       endif//TSUB_NEIGHBOR < 8
#   if  NTHREADS_NEIGHBOR == 512
/* TSUB_NEIGHBOR must be equal to 32 */
#          if  TSUB_NEIGHBOR != 32
#       undef  TSUB_NEIGHBOR
#       define TSUB_NEIGHBOR   (32)
#       endif//TSUB_NEIGHBOR != 32
#endif//NTHREADS_NEIGHBOR == 512
#endif//NTHREADS_NEIGHBOR >= 256
#endif//NTHREADS_NEIGHBOR >= 128
//-------------------------------------------------------------------------
#define NGROUPS_NEIGHBOR (NTHREADS_NEIGHBOR / TSUB_NEIGHBOR)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define NI_NEIGHBOR_ESTIMATE (NEIGHBOR_NUM << 3) */
#define NI_NEIGHBOR_ESTIMATE (NEIGHBOR_NUM << 2)
/* #define NI_NEIGHBOR_ESTIMATE (NEIGHBOR_NUM << 1) */
//-------------------------------------------------------------------------
/* minimum number of NI_NEIGHBOR_ESTIMATE is TSUB_NEIGHBOR to guarantee NBUF_NEIGHBOR is greater than unity */
#   if  NI_NEIGHBOR_ESTIMATE < TSUB_NEIGHBOR
#undef  NI_NEIGHBOR_ESTIMATE
#define NI_NEIGHBOR_ESTIMATE  (TSUB_NEIGHBOR)
#endif//NI_NEIGHBOR_ESTIMATE < TSUB_NEIGHBOR
//-------------------------------------------------------------------------
#define NBUF_NEIGHBOR (NI_NEIGHBOR_ESTIMATE / TSUB_NEIGHBOR)
/* #define NBUF_NEIGHBOR (4) */
//-------------------------------------------------------------------------
/* minimum number of NBUF_NEIGHBOR is 2 to satisfy NBUF_NEIGHBOR > 1 */
#   if  NBUF_NEIGHBOR < 2
#undef  NBUF_NEIGHBOR
#define NBUF_NEIGHBOR  (2)
#undef  NI_NEIGHBOR_ESTIMATE
#define NI_NEIGHBOR_ESTIMATE (TSUB_NEIGHBOR * NBUF_NEIGHBOR)
#endif//NBUF_NEIGHBOR < 2
//-------------------------------------------------------------------------
/* maximum number of NBUF_NEIGHBOR is 4 to use float4 in union */
#   if  NBUF_NEIGHBOR > 4
#undef  NBUF_NEIGHBOR
#define NBUF_NEIGHBOR  (4)
#undef  NI_NEIGHBOR_ESTIMATE
#define NI_NEIGHBOR_ESTIMATE (TSUB_NEIGHBOR * NBUF_NEIGHBOR)
#endif//NBUF_NEIGHBOR > 4
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MAKE_DEV_H
#       include "../tree/make_dev.h"
#endif//MAKE_DEV_H
//-------------------------------------------------------------------------
/* mycudaMalloc((void **)more0Buf, devProp.numSM * NBLOCKS_PER_SM_MAC      * (NTHREADS_MAC      / TSUB_MAC     ) * NUM_ALLOC_MACBUF      * sizeof(int)); */
/* mycudaMalloc((void **)more0Buf, devProp.numSM * NBLOCKS_PER_SM_NEIGHBOR * (NTHREADS_NEIGHBOR / TSUB_NEIGHBOR) * NUM_ALLOC_NEIGHBORBUF * sizeof(int)); */
//-------------------------------------------------------------------------
#define NUM_ALLOC_NEIGHBORBUF (NUM_ALLOC_MACBUF * (NBLOCKS_PER_SM_MAC * NGROUPS_MAC) / (NBLOCKS_PER_SM_NEIGHBOR * NGROUPS_NEIGHBOR))
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* options related to radix sorting */
//-------------------------------------------------------------------------
#define CHECK_2BITS
#define NTHREADS_SORT NTHREADS_NEIGHBOR
#define TSUB_SORT TSUB_NEIGHBOR
//-------------------------------------------------------------------------
#define SORT_ELEMENTS_PER_THREAD NBUF_NEIGHBOR
//-------------------------------------------------------------------------
#define KEY_BITS (32)
#include "../sort/radix_dev.h"
/* #define SORT_ONLY */
/* #include "../sort/radix_dev.h" */
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
/* limitation from capacity of shared memory */
//-------------------------------------------------------------------------
/* capacity of shared memory */
/* 48 KiB = 48 * 1024 B = 49152 B */
/* 16 KiB = 16 * 1024 B = 16384 B */
#ifdef  SMEM_PREF_FOR_NEIGHBOR_SEARCH
#       define NS_SHARED_MEM_CAPACITY (49152)
#else///SMEM_PREF_FOR_NEIGHBOR_SEARCH
#       define NS_SHARED_MEM_CAPACITY (16384)
#endif//SMEM_PREF_FOR_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------
#define REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
//-------------------------------------------------------------------------
/* usage of shared memory */
#       ifdef  REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
#define USAGE_NS_COMMON (12 * NTHREADS_NEIGHBOR * NBUF_NEIGHBOR)
#       else///REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
#define USAGE_NS_COMMON (16 * NTHREADS_NEIGHBOR * NBUF_NEIGHBOR)
#       endif//REDUCE_SM_USAGE_IN_NEIGHBOR_DEV_CU
/* #define USAGE_NS_SORT (8 * (NTHREADS_NEIGHBOR * SORT_ELEMENTS_PER_THREAD + NEIGHBOR_NUM * NGROUPS_NEIGHBOR)) */
#define USAGE_NS_SORT (8 * (NTHREADS_NEIGHBOR * SORT_ELEMENTS_PER_THREAD + NEIGHBOR_NUM_INC * NGROUPS_NEIGHBOR))
#       ifdef  USE_WARP_SHUFFLE_FUNC_NEIGHBOR
#define USAGE_NS_SHARE (0)
#       else///USE_WARP_SHUFFLE_FUNC_NEIGHBOR
#define USAGE_NS_SHARE (4 * NTHREADS_NEIGHBOR)
#       endif//USE_WARP_SHUFFLE_FUNC_NEIGHBOR
#       ifdef  USE_WARP_SHUFFLE_FUNC_SORT
#define USAGE_NS_TEMP (0)
#       else///USE_WARP_SHUFFLE_FUNC_SORT
#define USAGE_NS_TEMP (4 * NTHREADS_NEIGHBOR * SORT_ELEMENTS_PER_THREAD)
#       endif//USE_WARP_SHUFFLE_FUNC_SORT
#          if  RADIX_SORT_CHECK_BITS > 2
#define USAGE_NS_BUFFER (16 * NTHREADS_SORT)
#       else///RADIX_SORT_CHECK_BITS > 2
#define USAGE_NS_BUFFER (0)
#       endif//RADIX_SORT_CHECK_BITS > 2
//-------------------------------------------------------------------------
#define NBLOCKS_PER_SM_NEIGHBOR (NS_SHARED_MEM_CAPACITY / (USAGE_NS_COMMON + USAGE_NS_SORT + USAGE_NS_SHARE + USAGE_NS_TEMP + USAGE_NS_BUFFER))
//-------------------------------------------------------------------------
/* #   if  NBLOCKS_PER_SM_NEIGHBOR > 2 */
/* #undef  NBLOCKS_PER_SM_NEIGHBOR */
/* #define NBLOCKS_PER_SM_NEIGHBOR  (2) */
/* #endif//NBLOCKS_PER_SM_NEIGHBOR > 2 */
//-------------------------------------------------------------------------
/* #   if  NBLOCKS_PER_SM_NEIGHBOR > 1 */
/* #undef  NBLOCKS_PER_SM_NEIGHBOR */
/* #define NBLOCKS_PER_SM_NEIGHBOR  (1) */
/* #endif//NBLOCKS_PER_SM_NEIGHBOR > 1 */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
  int *gsync0, *gsync1;
  uint *freeLst;
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
  uint *freeNum;
  int *active;
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
} soaNeighborSearchBuf;
//-------------------------------------------------------------------------
#endif//FACILE_NEIGHBOR_SEARCH
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "neighbor_dev.cu"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
#ifdef  FACILE_NEIGHBOR_SEARCH
  //-----------------------------------------------------------------------
  void facileNeighborSearching_dev(const int Ni, const iparticle pi);
  //-----------------------------------------------------------------------
#else///FACILE_NEIGHBOR_SEARCH
  //-----------------------------------------------------------------------
  /* muse allocNeighborSearch_dev(int **gsync0, int **gsync1, deviceProp devProp); */
  muse allocNeighborSearch_dev(int **gsync0, int **gsync1, uint **freeLst,
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
			       uint **freeNum, int **active,
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
			       soaNeighborSearchBuf *buf, deviceProp devProp);
  void  freeNeighborSearch_dev(int  *gsync0, int  *gsync1, uint  *freeLst
#   if  !defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
			       , uint  *freeNum, int  *active
#endif//!defined(USE_SMID_TO_GET_BUFID_NEIGHBOR) && !defined(TRY_MODE_ABOUT_BUFFER_NEIGHBOR)
			       );
  //-----------------------------------------------------------------------
  void searchNeighbors_dev
  (const int Ni, const iparticle pi, const soaTreeCell cell, const soaTreeNode node,
   const soaMakeTreeBuf makeBuf, const soaNeighborSearchBuf searchBuf, deviceProp devProp);
  //-----------------------------------------------------------------------
#endif//FACILE_NEIGHBOR_SEARCH
  //-----------------------------------------------------------------------
#ifdef  SMEM_PREF_FOR_NEIGHBOR_SEARCH
  void setGlobalConstants_neighbor_dev_cu(void);
#endif//SMEM_PREF_FOR_NEIGHBOR_SEARCH
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//NEIGHBOR_DEV_H
//-------------------------------------------------------------------------
