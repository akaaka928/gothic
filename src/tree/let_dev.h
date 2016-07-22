/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/21(Thu) 11:54:57
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


//-------------------------------------------------------------------------
#   if  !defined(MPI_INCLUDED) && !defined(OMPI_MPI_H)
#       include <mpi.h>
#endif//!defined(MPI_INCLUDED) && !defined(OMPI_MPI_H)
//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef CUDALIB_H
#       include <cudalib.h>
#endif//CUDALIB_H
//-------------------------------------------------------------------------
#ifndef PEANO_H
#       include "../sort/peano.h"
#endif//PEANO_H
//-------------------------------------------------------------------------
#ifndef MACUTIL_H
#       include "../tree/macutil.h"
#endif//MACUTIL_H
//-------------------------------------------------------------------------
#ifndef MAKE_H
#       include "../tree/make.h"
#endif//MAKE_H
//-------------------------------------------------------------------------
#ifndef BUF_INC_H
#       include "../tree/buf_inc.h"
#endif//BUF_INC_H
//-------------------------------------------------------------------------
#ifndef MPICFG_H
#       include "../para/mpicfg.h"
#endif//MPICFG_H
//-------------------------------------------------------------------------
#ifndef LET_H
#       include "../tree/let.h"
#endif//LET_H
//-------------------------------------------------------------------------
#ifndef WALK_DEV_H
#       include "../tree/walk_dev.h"
#endif//WALK_DEV_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef NTHREADS_MAKE_LET
#define NTHREADS_MAKE_LET (128)
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
			  cudaStream_t **stream, int *Nstream, const deviceProp gpu, MPIcfg_tree mpi);
  void releaseLETtopology(domainInfo  *info, position  *ipos,
#ifdef  GADGET_MAC
			  real  *amin,
#endif//GADGET_MAC
#ifdef  BUILD_LET_ON_DEVICE
			  int  *numSend_hst, int  *numSend_dev,
#endif//BUILD_LET_ON_DEVICE
			  cudaStream_t  *stream, int  Nstream);
  //-----------------------------------------------------------------------
  void callGenLET
  (const cudaStream_t stream, domainInfo *let, MPIcfg_tree mpi,
   const soaTreeNode tree, const int numSendGuess, const soaTreeWalkBuf buf
#ifdef  MONITOR_LETGEN_TIME
   , unsigned long long int * RESTRICT cycles
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
