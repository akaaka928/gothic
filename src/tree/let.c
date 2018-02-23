/**
 * @file let.c
 *
 * @brief Source code for building locally essential tree (LET)
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/02/22 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

#include "macro.h"
#include "mpilib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../para/mpicfg.h"

#include "make.h"
#include "let.h"


#ifndef SERIALIZED_EXECUTION


#define MPI_TAG_LET(rank, size) ((rank) + ((size) * 20))


/**
 * @fn shareNodePosition
 *
 * @brief Exchange position and size of the root cell.
 * @detail position contains (x, y, z, and bmax2)
 */
void shareNodePosition(const int Ndomain, domainInfo *info, position *pos_ful, const position pos_loc,
#ifdef  GADGET_MAC
		       real *acc_ful, const real acc_loc,
#endif//GADGET_MAC
		       MPIcfg_tree mpi)
{
  __NOTE__("%s\n", "start");


  chkMPIerr(MPI_Allgather(&pos_loc, 1,    mpi.ipos, pos_ful, 1,    mpi.ipos, mpi.comm));
#ifdef  GADGET_MAC
  chkMPIerr(MPI_Allgather(&acc_loc, 1, MPI_REALDAT, acc_ful, 1, MPI_REALDAT, mpi.comm));
#endif//GADGET_MAC

  for(int ii = 0; ii < Ndomain - 1; ii++){
    info[ii].icom = pos_ful[info[ii].rank];
#ifdef  GADGET_MAC
    info[ii].amin = acc_ful[info[ii].rank];
#endif//GADGET_MAC
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn setLETpartition
 *
 * @brief Set memory alignment for send and receive buffers.
 */
void setLETpartition(const int Ndomain, domainInfo *info)
{
  __NOTE__("%s\n", "start");

  for(int ii = 0; ii < Ndomain - 1; ii++){
    info[ii].maxSend = ALIGN_BUF_FOR_LET(info[ii].maxSend);
    info[ii].maxRecv = ALIGN_BUF_FOR_LET(info[ii].maxRecv);
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn guessLETpartition
 *
 * @brief Set partition for LET buffer based on guess of the required size of send and receive buffers.
 *
 * @param (icom) physical properties of enclosing ball (x, y, z, and r2)
 */
void guessLETpartition(const int Ndomain, domainInfo *info, const int numNode, const position icom, MPIcfg_tree mpi)
{
  __NOTE__("%s (icom = (%e, %e, %e, %e) @ rank %d)\n", "start", icom.x, icom.y, icom.z, icom.m, mpi.rank);


  /** exchange # of nodes in full tree */
  const int numMax = ALIGN_BUF_FOR_LET(numNode);

  for(int ii = 0; ii < Ndomain - 1; ii++){
    /** initialization */
    info[ii].overEstimateSend = 0;
    info[ii].overEstimateRecv = 0;
    info[ii].numFull = numNode;

    chkMPIerr(MPI_Isend(&          numMax  , 1, MPI_INT, info[ii].rank, MPI_TAG_LET(     mpi.rank, mpi.size), mpi.comm, &(info[ii].reqSendInfo)));
    chkMPIerr(MPI_Irecv(&(info[ii].numNode), 1, MPI_INT, info[ii].rank, MPI_TAG_LET(info[ii].rank, mpi.size), mpi.comm, &(info[ii].reqRecvInfo)));
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */

  for(int ii = 0; ii < Ndomain - 1; ii++){
    MPI_Status statusRecv;    chkMPIerr(MPI_Wait(&(info[ii].reqRecvInfo), &statusRecv));
    MPI_Status statusSend;    chkMPIerr(MPI_Wait(&(info[ii].reqSendInfo), &statusSend));
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */


  /** guess required depth to get enough accuracy */
  for(int ii = 0; ii < Ndomain - 1; ii++){
    const real rx = icom.x - info[ii].icom.x;
    const real ry = icom.y - info[ii].icom.y;
    const real rz = icom.z - info[ii].icom.z;
    const real r2 = FLT_MIN + rx * rx + ry * ry + rz * rz;

    const real lambda = FMAX(UNITY - SQRTRATIO(info[ii].icom.m, r2), ZERO);
    /* const real theta = SQRTRATIO(icom.m, EPSILON + (lambda * lambda) * r2);/\**< rough estimated value of theta *\/ */

    /** assumption: LET corresponds to full tree in theta <= thetamin */
    /** rough assumption: minimum value of theta is 0.0625 = 2^-4 -->> 1/thetamin^2 = 2^8 = 256*/
    /** rough assumption: computaional amount scales as thetamin^-2 */

    /** numLET = numNode * (theta^2 / thetamin^2); */
    /** int numLET = (int)CEIL((real)numNode * FMIN((real)256.0 * icom.m / (EPSILON + (lambda * lambda) * r2), UNITY)); */
    real factor = LDEXP(UNITY, FMIN(CEIL(LOG2((real)256.0 * icom.m / (EPSILON + (lambda * lambda) * r2))), ZERO));
    factor *= LETSIZE_OVERESTIMATION_FACTOR;
    const int numLET = ALIGN_BUF_FOR_LET((int)CEIL((real)numNode * factor));
    info[ii].maxSend = (numLET < numMax) ? numLET : numMax;

    chkMPIerr(MPI_Isend(&(info[ii].maxSend), 1, MPI_INT, info[ii].rank,      mpi.rank, mpi.comm, &(info[ii].reqSendInfo)));
    chkMPIerr(MPI_Irecv(&(info[ii].maxRecv), 1, MPI_INT, info[ii].rank, info[ii].rank, mpi.comm, &(info[ii].reqRecvInfo)));
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */

  for(int ii = 0; ii < Ndomain - 1; ii++){
    MPI_Status statusSend;    chkMPIerr(MPI_Wait(&(info[ii].reqSendInfo), &statusSend));
    MPI_Status statusRecv;    chkMPIerr(MPI_Wait(&(info[ii].reqRecvInfo), &statusRecv));
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */

  setLETpartition(Ndomain, info);


  /** count total number of the roughly estimated LET nodes */
  /** if the estimated values are small, then enlarge values for safety */
  int numSend = 0;
  int numRecv = 0;
  for(int ii = 0; ii < Ndomain - 1; ii++){
    numSend += info[ii].maxSend;
    numRecv += info[ii].maxRecv;
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */

#if 1
  const int upper = (int)floorf(0.9f * EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - ALIGN_BUF_FOR_LET(numNode);
  if( 2 * (numSend + numRecv) < upper ){
    const float factor = floorf((float)upper / (float)(numSend + numRecv));

    for(int ii = 0; ii < Ndomain - 1; ii++){
      if( info[ii].maxSend !=          numMax  )	info[ii].maxSend = ALIGN_BUF_FOR_LET((int)floorf(factor * (float)info[ii].maxSend));
      if( info[ii].maxRecv != info[ii].numNode )	info[ii].maxRecv = ALIGN_BUF_FOR_LET((int)floorf(factor * (float)info[ii].maxRecv));
    }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */
  }/* if( 2 * (numSend + numRecv) < upper ){ */
#endif

  if( info[0].maxSend + info[0].maxRecv > ((int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - ALIGN_BUF_FOR_LET(numNode)) ){
    __KILL__(stderr, "ERROR: rough estimation routine predicts # of LET nodes(%d) succeed the capacity of the array(%d - %d); communication with first partner would fail due to lack of MPI buffer.\n", info[0].maxSend + info[0].maxRecv, (int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE), ALIGN_BUF_FOR_LET(numNode));
  }/* if( numSend + numRecv > (int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - ALIGN_BUF_FOR_LET(numNode) ){ */


  __NOTE__("%s\n", "end");
}


#endif//SERIALIZED_EXECUTION
