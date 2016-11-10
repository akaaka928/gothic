/*************************************************************************\
 *                                                                       *
                  last updated on 2016/10/28(Fri) 16:31:42
 *                                                                       *
 *    Header File for constructing octree structure                      *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef LET_H
#define LET_H
//-------------------------------------------------------------------------
#include <mpi.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../sort/peano.h"
//-------------------------------------------------------------------------
#include "../tree/macutil.h"
#include "../tree/make.h"
//-------------------------------------------------------------------------
#include "../para/mpicfg.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define EXTEND_NUM_TREE_NODE (2.0f)
//-------------------------------------------------------------------------
#ifdef  SERIALIZED_EXECUTION
#       undef  EXTEND_NUM_TREE_NODE
#       define EXTEND_NUM_TREE_NODE (1.0f)
#else///SERIALIZED_EXECUTION
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define LETSIZE_OVERESTIMATION_STEPS (4)
//-------------------------------------------------------------------------
/* #define LETSIZE_OVERESTIMATION_FACTOR (2.0f) */
/* #define LETSIZE_OVERESTIMATION_FACTOR (1.5f) */
#define LETSIZE_OVERESTIMATION_FACTOR (1.25f)
/* #define LETSIZE_OVERESTIMATION_FACTOR (1.125f) */
//-------------------------------------------------------------------------
#define LETSIZE_REDUCE_CRITERION (0.25f)
#define LETSIZE_REDUCE_FACTOR (0.75f)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define TGRP (1)
//-------------------------------------------------------------------------
/* /\* used capacity of shared memory *\/ */
/* TGRP * (1 + NBUF_LET) * (1 + 4 + 1) * 4 bytes */
/* TGRP * (1 + NBUF_MSK)               * 4 bytes */
/* TGRP                                * 4 bytes */
/* sum: TGRP * 4 * (1 + (1 + NBUF_MSK) + 6 * (1 + NBUF_LET)) */
/*    = TGRP * 4 * (8 + NBUF_MSK + 6 * NBUF_LET) bytes */
/* NBUF_MSK = 2 * NBUS_MSK ==> sum = TGRP * 4 * 8 * (1 + NBUF_LET) */
/* therefore, NBUF_LET = 2^n - 1 is the condition */
//-------------------------------------------------------------------------
#define NBUF_LET (7)
#define NBUF_MSK (2 * NBUF_LET)
//-------------------------------------------------------------------------
#define LET_INCNUM_UNIT (32)
/* #define ALIGN_BUF_FOR_LET(ini) ((ini) + (((ini) %  LET_INCNUM_UNIT     ) != 0) * (LET_INCNUM_UNIT - ((ini) & (LET_INCNUM_UNIT - 1)))) */
#define ALIGN_BUF_FOR_LET(ini) ((ini) + (((ini) & (LET_INCNUM_UNIT - 1)) != 0) * (LET_INCNUM_UNIT - ((ini) & (LET_INCNUM_UNIT - 1))))
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
  position icom;
#ifdef  BUILD_LET_ON_DEVICE
  int *numSend_hst, *numSend_dev;
#endif//BUILD_LET_ON_DEVICE
  MPI_Request reqSendInfo, reqSendMore, reqSendJpos, reqSendMass;
  MPI_Request reqRecvInfo, reqRecvMore, reqRecvJpos, reqRecvMass;
  int headSend, headRecv, numSend, numRecv, numRecvGuess;
  int maxSend, maxRecv, numFull, numNode;/* numFull is # of nodes in local full tree, numNode is # of nodes in distant full tree */
  int overEstimateSend, overEstimateRecv;
  int rank;
#ifdef  GADGET_MAC
  real amin;
#endif//GADGET_MAC
} domainInfo;
//-------------------------------------------------------------------------
/* #ifndef SERIALIZED_EXECUTION */
/* typedef struct */
/* { */
/*   domainInfo *let; */
/*   int2 *list; */
/*   uint *mask; */
/* } soaLET; */
/* #endif//SERIALIZED_EXECUTION */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "let.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  /* muse allocArrays2BuildLET(uint **mask, int2 **list); */
  /* void  freeArrays2BuildLET(uint  *mask, int2  *list); */
  //-----------------------------------------------------------------------
  void initMask(const int numNode, int2 *list, uint *mask);
  //-----------------------------------------------------------------------
/* #ifndef BUILD_LET_ON_DEVICE */
/*   muse  configLETtopology(domainInfo **info, position **ipos, MPIcfg_tree mpi); */
/*   void releaseLETtopology(domainInfo  *info, position  *ipos); */
/* #endif//BUILD_LET_ON_DEVICE */
/*   /\* void  commitLETnodes   (const int numNode, PHinfo *level, uint *ptag, int2 *list, int *leafLevel); *\/ */
  //-----------------------------------------------------------------------
  void shareNodePosition(const int Ndomain, domainInfo *info, position *pos_ful, const position pos_loc,
#ifdef  GADGET_MAC
			 real *acc_ful, const real acc_loc,
#endif//GADGET_MAC
			 MPIcfg_tree mpi);
  //-----------------------------------------------------------------------
  void   setLETpartition(const int Ndomain, domainInfo *info);
  void guessLETpartition(const int Ndomain, domainInfo *info, const int numNode, const position icom, MPIcfg_tree mpi);
  //-----------------------------------------------------------------------
  /* void guessLETpartition(const int Ndomain, domainInfo *info, const int numNode, const jparticle pj_root_hst, const real bmax2_root_hst, MPIcfg_tree mpi); */
  //-----------------------------------------------------------------------
  /* void makeLET */
  /* (const position icom, */
  /*  const int leafLevel, int2 * RESTRICT list, uint * RESTRICT mask, */
  /*  const uint * RESTRICT more_org, const jparticle * RESTRICT jpos_org, const real * RESTRICT mj_org, */
  /*  uint       * RESTRICT more_let,       jparticle * RESTRICT jpos_let,       real * RESTRICT mj_let, */
  /*  int *numLETnode); */
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//SERIALIZED_EXECUTION
//-------------------------------------------------------------------------
#endif//LET_H
//-------------------------------------------------------------------------
