/*************************************************************************\
 *                                                                       *
                  last updated on 2016/11/01(Tue) 10:23:50
 *                                                                       *
 *    Header File for constructing octree structure                      *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef LET_DEV_H
#define LET_DEV_H
//-------------------------------------------------------------------------
#include <mpi.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../sort/peano.h"
//-------------------------------------------------------------------------
#include "../para/mpicfg.h"
//-------------------------------------------------------------------------
#include "../tree/macutil.h"
#include "../tree/make.h"
#include "../tree/buf_inc.h"
#include "../tree/let.h"
#include "../tree/walk_dev.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef NTHREADS_MAKE_LET
/* #define NTHREADS_MAKE_LET (128) */
/* #define NTHREADS_MAKE_LET (256) */
/* #define NTHREADS_MAKE_LET (512) */
#define NTHREADS_MAKE_LET (1024)
#endif//NTHREADS_MAKE_LET
//-------------------------------------------------------------------------
/* maximum value of NTHREADS_MAKE_LET is 1024 to calculate prefix sum in 2 stages (32^2) */
#   if  NTHREADS_MAKE_LET > 1024
#undef  NTHREADS_MAKE_LET
#define NTHREADS_MAKE_LET  (1024)
#endif//NTHREADS_MAKE_LET > 1024
//-------------------------------------------------------------------------
#define USE_WARP_SHUFFLE_FUNC_MAKE_LET
//-------------------------------------------------------------------------
#   if  defined(USE_WARP_SHUFFLE_FUNC_MAKE_LET) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_MAKE_LET
#endif//defined(USE_WARP_SHUFFLE_FUNC_MAKE_LET) && (GPUGEN < 30)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  DIV_NTHREADS_MAKE_LET
#undef  DIV_NTHREADS_MAKE_LET
#endif//DIV_NTHREADS_MAKE_LET
//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* # of blocks per SM is unity */
/* SM size is 48KiB or 16KiB */
/* uint is 4B */
/* # of elements in SM is 12K or 4K */
/* SM usage is NTHREADS_MAKE_LET * (NQUEUE_LET + 1) */
//-------------------------------------------------------------------------
#ifndef NQUEUE_LET
#          if  SMPREF_LET == 1
#define NQUEUE_LET (DIV_NTHREADS_MAKE_LET(12288 / NBLOCKS_PER_SM) - 1)
#       else///SMPREF_LET == 1
#define NQUEUE_LET (DIV_NTHREADS_MAKE_LET( 4096 / NBLOCKS_PER_SM) - 1)
#       endif//SMPREF_LET == 1
#endif//NQUEUE_LET
//-------------------------------------------------------------------------
/* #ifndef NQUEUE_LET */
/* #          if  SMPREF_LET == 1 */
/* #define NQUEUE_LET (DIV_NTHREADS_MAKE_LET(12288) - 1) */
/* #       else///SMPREF_LET == 1 */
/* #define NQUEUE_LET (DIV_NTHREADS_MAKE_LET( 4096) - 1) */
/* #       endif//SMPREF_LET == 1 */
/* #endif//NQUEUE_LET */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "let_dev.cu"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  muse  configLETtopology(domainInfo **info, position **ipos,
#ifdef  GADGET_MAC
			  real **amin,
#endif//GADGET_MAC
#ifdef  BUILD_LET_ON_DEVICE
			  int **numSend_hst, int **numSend_dev,
#endif//BUILD_LET_ON_DEVICE
#ifdef  ALLOCATE_LETBUFFER
			  uint **buf, soaTreeWalkBuf *treebuf,
#endif//ALLOCATE_LETBUFFER
/* #ifdef  USE_CUDA_EVENT */
/* 			  cudaEvent_t **iniEvent, cudaEvent_t **finEvent, */
/* #endif//USE_CUDA_EVENT */
			  cudaStream_t **stream, int *Nstream, const deviceProp gpu, MPIcfg_tree mpi);
  void releaseLETtopology(domainInfo  *info, position  *ipos,
#ifdef  GADGET_MAC
			  real  *amin,
#endif//GADGET_MAC
#ifdef  BUILD_LET_ON_DEVICE
			  int  *numSend_hst, int  *numSend_dev,
#endif//BUILD_LET_ON_DEVICE
#ifdef  ALLOCATE_LETBUFFER
			  uint  *buf,
#endif//ALLOCATE_LETBUFFER
/* #ifdef  USE_CUDA_EVENT */
/* 			  cudaEvent_t  *iniEvent, cudaEvent_t  *finEvent, */
/* #endif//USE_CUDA_EVENT */
			  cudaStream_t  *stream, int  Nstream);
  //-----------------------------------------------------------------------
#ifdef  DBG_LETGEN_ON_GPU
  void printTreeNode(const int Nj, uint * RESTRICT more, jparticle * RESTRICT jpos, real * RESTRICT mj);
#endif//DBG_LETGEN_ON_GPU
  void callGenLET
  (const cudaStream_t stream, domainInfo *let, MPIcfg_tree mpi, const soaTreeNode tree, const soaTreeWalkBuf buf
#ifdef  MONITOR_LETGEN_TIME
#ifdef  USE_CUDA_EVENT
   , const cudaEvent_t iniEvent, const cudaEvent_t finEvent
#else///USE_CUDA_EVENT
   , unsigned long long int * RESTRICT cycles
#endif//USE_CUDA_EVENT
#endif//MONITOR_LETGEN_TIME
   );
  //-----------------------------------------------------------------------
  void setGlobalConstants_let_dev_cu
  (
#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
   const real theta2_hst
#else///endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
   void
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
   );
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//LET_DEV_H
//-------------------------------------------------------------------------
