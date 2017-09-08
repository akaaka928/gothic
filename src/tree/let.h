/**
 * @file let.h
 *
 * @brief Header file for building locally essential tree (LET)
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/09/07 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef LET_H
#define LET_H


#ifndef SERIALIZED_EXECUTION
#include <mpi.h>
#endif//SERIALIZED_EXECUTION

#include "macro.h"
#include "cudalib.h"

#include "../sort/peano.h"
#include "../tree/macutil.h"
#include "../tree/make.h"

#ifndef SERIALIZED_EXECUTION
#include "../para/mpicfg.h"
#endif//SERIALIZED_EXECUTION


/**
 * @def EXTEND_NUM_TREE_NODE
 *
 * @brief extend number of arrays for tree node
 */
#define EXTEND_NUM_TREE_NODE (2.0f)


#ifdef  SERIALIZED_EXECUTION
#undef  EXTEND_NUM_TREE_NODE
#define EXTEND_NUM_TREE_NODE (1.0f)
#else///SERIALIZED_EXECUTION


/**
 * @def LETSIZE_OVERESTIMATION_STEPS
 *
 * @brief a parameter to shrink size of LET buffer
 */
#define LETSIZE_OVERESTIMATION_STEPS (4)

/**
 * @def LETSIZE_REDUCE_CRITERION
 *
 * @brief a parameter to shrink size of LET buffer
 */
#define LETSIZE_REDUCE_CRITERION (0.25f)

/**
 * @def LETSIZE_REDUCE_FACTOR
 *
 * @brief a parameter to shrink size of LET buffer
 */
#define LETSIZE_REDUCE_FACTOR (0.75f)


/**
 * @def LETSIZE_OVERESTIMATION_FACTOR
 *
 * @brief a parameter to predict size of LET buffer
 */
/* #define LETSIZE_OVERESTIMATION_FACTOR (2.0f) */
/* #define LETSIZE_OVERESTIMATION_FACTOR (1.5f) */
#define LETSIZE_OVERESTIMATION_FACTOR (1.25f)
/* #define LETSIZE_OVERESTIMATION_FACTOR (1.125f) */


#define LET_INCNUM_UNIT (32)
#define ALIGN_BUF_FOR_LET(ini) ((ini) + (((ini) & (LET_INCNUM_UNIT - 1)) != 0) * (LET_INCNUM_UNIT - ((ini) & (LET_INCNUM_UNIT - 1))))


/**
 * @struct domainInfo
 *
 * @brief structure for LET
 */
typedef struct
{
  position icom;
  int *numSend_hst, *numSend_dev;
  MPI_Request reqSendInfo, reqRecvInfo;
#ifdef  MPI_ONE_SIDED_FOR_LET
  MPI_Request reqSendHead, reqRecvHead;
  int headDisp;
#else///MPI_ONE_SIDED_FOR_LET
  MPI_Request reqSendMore, reqRecvMore;
  MPI_Request reqSendJpos, reqRecvJpos;
  MPI_Request reqSendMass, reqRecvMass;
#endif//MPI_ONE_SIDED_FOR_LET
  int headSend, headRecv, numSend, numRecv, numRecvGuess;
  int maxSend, maxRecv, numFull, numNode;/**< numFull is # of nodes in local full tree, numNode is # of nodes in distant full tree */
  int overEstimateSend, overEstimateRecv;
  int rank;
#ifdef  GADGET_MAC
  real amin;
#endif//GADGET_MAC
} domainInfo;


/* list of functions appeared in ``let.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  void initMask(const int numNode, int2 *list, uint *mask);

  void shareNodePosition(const int Ndomain, domainInfo *info, position *pos_ful, const position pos_loc,
#ifdef  GADGET_MAC
			 real *acc_ful, const real acc_loc,
#endif//GADGET_MAC
			 MPIcfg_tree mpi);

  void   setLETpartition(const int Ndomain, domainInfo *info);
  void guessLETpartition(const int Ndomain, domainInfo *info, const int numNode, const position icom, MPIcfg_tree mpi);
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//SERIALIZED_EXECUTION
#endif//LET_H
