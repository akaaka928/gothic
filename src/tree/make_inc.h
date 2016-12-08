/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:55:39
 *                                                                       *
 *    Header File for duplicating parallel prefix sum library on GPU     *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef MAKE_DEL_H
#include "../tree/make_del.h"
#endif//MAKE_DEL_H
//-------------------------------------------------------------------------
#ifndef MAKE_INC_H
#define MAKE_INC_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define USE_WARP_SHUFFLE_FUNC_MAKE_INC
//-------------------------------------------------------------------------
#   if  defined(USE_WARP_SHUFFLE_FUNC_MAKE_INC) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_MAKE_INC
#endif//defined(USE_WARP_SHUFFLE_FUNC_MAKE_INC) && (GPUGEN < 30)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  NTHREADS_MAKE_INC ==  128
#define PREFIX_SUM_MAKE_INC_BLCK prefixSumMakeBlck0128
#define PREFIX_SUM_MAKE_INC_GRID prefixSumMakeGrid0128
#endif//NTHREADS_MAKE_INC ==  128
//-------------------------------------------------------------------------
#   if  NTHREADS_MAKE_INC ==  256
#define PREFIX_SUM_MAKE_INC_BLCK prefixSumMakeBlck0256
#define PREFIX_SUM_MAKE_INC_GRID prefixSumMakeGrid0256
#endif//NTHREADS_MAKE_INC ==  256
//-------------------------------------------------------------------------
#   if  NTHREADS_MAKE_INC ==  512
#define PREFIX_SUM_MAKE_INC_BLCK prefixSumMakeBlck0512
#define PREFIX_SUM_MAKE_INC_GRID prefixSumMakeGrid0512
#endif//NTHREADS_MAKE_INC ==  512
//-------------------------------------------------------------------------
#   if  NTHREADS_MAKE_INC == 1024
#define PREFIX_SUM_MAKE_INC_BLCK prefixSumMakeBlck1024
#define PREFIX_SUM_MAKE_INC_GRID prefixSumMakeGrid1024
#endif//NTHREADS_MAKE_INC == 1024
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  MAKE_DEL_H
#undef  MAKE_DEL_H
#endif//MAKE_DEL_H
//-------------------------------------------------------------------------
#endif//MAKE_INC_H
//-------------------------------------------------------------------------
