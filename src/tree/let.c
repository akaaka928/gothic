/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:48:46
 *                                                                       *
 *    Octree N-body calculation for collisionless systems on NVIDIA GPUs *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
//-------------------------------------------------------------------------
#include "macro.h"
#include "mpilib.h"
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../para/mpicfg.h"
//-------------------------------------------------------------------------
#include "make.h"
#include "let.h"
//-------------------------------------------------------------------------


/* //------------------------------------------------------------------------- */
/* /\* arrays to construct LET nodes *\/ */
/* //------------------------------------------------------------------------- */
/* muse allocArrays2BuildLET(uint **mask, int2 **list) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "start"); */
/*   //----------------------------------------------------------------------- */
/*   muse alloc = {0, 0}; */
/*   *mask = (uint *)malloc(NUM_ALLOC_TREE_NODE * sizeof(uint));  if( *mask == NULL ){    __KILL__(stderr, "ERROR: failure to allocate mask\n");  } */
/*   alloc.host +=          NUM_ALLOC_TREE_NODE * sizeof(uint) ; */
/*   *list = (int2 *)malloc(NUM_PHKEY_LEVEL     * sizeof(int2));  if( *list == NULL ){    __KILL__(stderr, "ERROR: failure to allocate list\n");  } */
/*   alloc.host +=          NUM_ALLOC_TREE_NODE * sizeof(int2) ; */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "end"); */
/*   //----------------------------------------------------------------------- */
/*   return (alloc); */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
/* void  freeArrays2BuildLET(uint  *mask, int2  *list) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "start"); */
/*   //----------------------------------------------------------------------- */
/*   free(mask); */
/*   free(list); */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "end"); */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */


/* //------------------------------------------------------------------------- */
/* void initMask(const int numNode, int2 *list, uint *mask) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "start"); */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* initialize ``mask'' array *\/ */
/*   //----------------------------------------------------------------------- */
/*   const int idxHead = list[0].x; */
/*   const int numHead = list[0].y; */
/*   //----------------------------------------------------------------------- */
/*   /\* the values of ``mask'' belongs in the top level must be unity *\/ */
/*   for(int ii = idxHead; ii < idxHead + numHead; ii++) */
/*     mask[ii] = 1; */
/*   //----------------------------------------------------------------------- */
/*   /\* all the values of ``mask'' must be zero *\/ */
/*   for(int ii = idxHead + numHead; ii < idxHead + numNode; ii++) */
/*     mask[ii] = 0; */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "end"); */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */


//-------------------------------------------------------------------------
/* #ifndef BUILD_LET_ON_DEVICE */
/* //------------------------------------------------------------------------- */
/* muse configLETtopology(domainInfo **info, position **ipos, MPIcfg_tree mpi) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "start"); */
/*   //----------------------------------------------------------------------- */
/*   muse alloc = {0, 0}; */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   *info = (domainInfo *)malloc(mpi.size * sizeof(domainInfo));  alloc.host +=                mpi.size * sizeof(domainInfo) ; */
/*   *ipos = (position   *)malloc(mpi.size * sizeof(position  ));  alloc.host +=                mpi.size * sizeof(position  ) ; */
/*   if( *info == NULL ){    __KILL__(stderr, "ERROR: failure to allocate info\n");  } */
/*   if( *ipos == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ipos\n");  } */
/*   //----------------------------------------------------------------------- */
/*   for(int ii = 0; ii < mpi.size - 1; ii++) */
/*     (*info)[ii].rank = mpi.rank ^ (1 + ii); */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "end"); */
/*   //----------------------------------------------------------------------- */
/*   return (alloc); */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
/* void releaseLETtopology(domainInfo  *info, position  *ipos) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "start"); */
/*   //----------------------------------------------------------------------- */
/*   free(info); */
/*   free(ipos); */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "end"); */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
/* #endif//BUILD_LET_ON_DEVICE */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* void commitLETnodes(const int numNode, PHinfo *level, uint *ptag, int2 *list, int *leafLevel) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "start"); */
/*   //----------------------------------------------------------------------- */
/*   /\* commit # of local tree nodes in each tree level *\/ */
/*   *leafLevel = NUM_PHKEY_LEVEL; */
/*   for(int ii = 0; ii < NUM_PHKEY_LEVEL; ii++){ */
/*     list[ii].x = ptag[level[ii].head] & IDXMASK;/\* head index: *\/ */
/*     if( list[ii].x == NULL_NODE ){ */
/*       *leafLevel = ii; */
/*       break; */
/*     } */
/*   } */
/*   for(int ii = 0; ii < (*leafLevel) - 1; ii++) */
/*     list[ii].y = list[ii + 1].x - list[ii].x;/\* number *\/ */
/*   list[(*leafLevel) - 1].y = numNode - list[(*leafLevel) - 1].x; */
/*   //----------------------------------------------------------------------- */
/* #if 0 */
/*   fprintf(stderr, "leafLevel = %d, numNode = %d\n", *leafLevel, numNode); */
/*   for(int ii = 0; ii < *leafLevel; ii++) */
/*     fprintf(stderr, "ii = %d: x = %d, y = %d\n", ii, list[ii].x, list[ii].y); */
/*   fflush(stderr); */
/*   MPI_Barrier(MPI_COMM_WORLD); */
/*   MPI_Finalize(); */
/*   exit(0); */
/* #endif */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "end"); */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* /\* send position, r2 of the root cell *\/ */
/* position mypos = {pj_root_hst.x, pj_root_hst.y, pj_root_hst.z, bmax2_root_hst}; */
//-------------------------------------------------------------------------
void shareNodePosition(const int Ndomain, domainInfo *info, position *pos_ful, const position pos_loc,
#ifdef  GADGET_MAC
		       real *acc_ful, const real acc_loc,
#endif//GADGET_MAC
		       MPIcfg_tree mpi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  chkMPIerr(MPI_Allgather(&pos_loc, 1,    mpi.ipos, pos_ful, 1,    mpi.ipos, mpi.comm));
#ifdef  GADGET_MAC
  chkMPIerr(MPI_Allgather(&acc_loc, 1, MPI_REALDAT, acc_ful, 1, MPI_REALDAT, mpi.comm));
#endif//GADGET_MAC
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ndomain - 1; ii++){
    info[ii].icom = pos_ful[info[ii].rank];
#ifdef  GADGET_MAC
    info[ii].amin = acc_ful[info[ii].rank];
#endif//GADGET_MAC
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* static inline int increaseNumber(const int ini) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   /\* int fin = ini + (ini / 4); *\/ */
/*   /\* fin += ((fin % LET_INCNUM_UNIT) != 0) * (LET_INCNUM_UNIT - (fin & (LET_INCNUM_UNIT - 1))); *\/ */
/*   /\* //----------------------------------------------------------------------- *\/ */
/*   /\* return (fin); *\/ */
/*   //----------------------------------------------------------------------- */
/*   return (ini + ((ini % LET_INCNUM_UNIT) != 0) * (LET_INCNUM_UNIT - (ini & (LET_INCNUM_UNIT - 1)))); */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------
void setLETpartition(const int Ndomain, domainInfo *info)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ndomain - 1; ii++){
    //---------------------------------------------------------------------
    info[ii].maxSend = ALIGN_BUF_FOR_LET(info[ii].maxSend);
    info[ii].maxRecv = ALIGN_BUF_FOR_LET(info[ii].maxRecv);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* position {x, y, z, m} を enclosing ball 格納用（変数としては，icom に置いてあげている）に使う */
/* x, y, z --> position */
/* m --> radius squared */
void guessLETpartition(const int Ndomain, domainInfo *info, const int numNode, const position icom, MPIcfg_tree mpi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s (icom = (%e, %e, %e, %e) @ rank %d)\n", "start", icom.x, icom.y, icom.z, icom.m, mpi.rank);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* exchange # of nodes in full tree */
  //-----------------------------------------------------------------------
  const int numMax = ALIGN_BUF_FOR_LET(numNode);
  __NOTE__("numMax = %d, numNode = %d @ rank %d\n", numMax, numNode, mpi.rank);
  for(int ii = 0; ii < Ndomain - 1; ii++){
    //---------------------------------------------------------------------
    /* initialization */
    info[ii].overEstimateSend = 0;
    info[ii].overEstimateRecv = 0;
    info[ii].numFull = numNode;
    //---------------------------------------------------------------------
    chkMPIerr(MPI_Isend(&          numMax  , 1, MPI_INT, info[ii].rank,      mpi.rank, mpi.comm, &(info[ii].reqSendInfo)));
    chkMPIerr(MPI_Irecv(&(info[ii].numNode), 1, MPI_INT, info[ii].rank, info[ii].rank, mpi.comm, &(info[ii].reqRecvInfo)));
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ndomain - 1; ii++){
    //---------------------------------------------------------------------
    MPI_Status statusRecv;    chkMPIerr(MPI_Wait(&(info[ii].reqRecvInfo), &statusRecv));
    MPI_Status statusSend;    chkMPIerr(MPI_Wait(&(info[ii].reqSendInfo), &statusSend));
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* guess required depth to get enough accuracy */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ndomain - 1; ii++){
    //---------------------------------------------------------------------
    const real rx = icom.x - info[ii].icom.x;
    const real ry = icom.y - info[ii].icom.y;
    const real rz = icom.z - info[ii].icom.z;
    const real r2 = FLT_MIN + rx * rx + ry * ry + rz * rz;
    //---------------------------------------------------------------------
    const real lambda = FMAX(UNITY - SQRTRATIO(info[ii].icom.m, r2), ZERO);
    /* const real theta = SQRTRATIO(icom.m, EPSILON + (lambda * lambda) * r2);/\* rough estimated value of theta *\/ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* assumption: LET corresponds to full tree in theta <= thetamin */
    /* rough assumption: minimum value of theta is 0.0625 = 2^-4 -->> 1/thetamin^2 = 2^8 = 256*/
    /* rough assumption: computaional amount scales as thetamin^-2 */
    //---------------------------------------------------------------------
    /* numLET = numNode * (theta^2 / thetamin^2); */
    /* int numLET = (int)CEIL((real)numNode * FMIN((real)256.0 * icom.m / (EPSILON + (lambda * lambda) * r2), UNITY)); */
    real factor = LDEXP(UNITY, FMIN(CEIL(LOG2((real)256.0 * icom.m / (EPSILON + (lambda * lambda) * r2))), ZERO));
    factor *= LETSIZE_OVERESTIMATION_FACTOR;
    const int numLET = ALIGN_BUF_FOR_LET((int)CEIL((real)numNode * factor));
    __NOTE__("factor = %e, numLET = %d, numMax = %d @ rank %d\n", factor, numLET, numMax, mpi.rank);
    info[ii].maxSend = (numLET < numMax) ? numLET : numMax;
    //---------------------------------------------------------------------
    chkMPIerr(MPI_Isend(&(info[ii].maxSend), 1, MPI_INT, info[ii].rank,      mpi.rank, mpi.comm, &(info[ii].reqSendInfo)));
    chkMPIerr(MPI_Irecv(&(info[ii].maxRecv), 1, MPI_INT, info[ii].rank, info[ii].rank, mpi.comm, &(info[ii].reqRecvInfo)));
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ndomain - 1; ii++){
    //---------------------------------------------------------------------
    MPI_Status statusSend;    chkMPIerr(MPI_Wait(&(info[ii].reqSendInfo), &statusSend));
    MPI_Status statusRecv;    chkMPIerr(MPI_Wait(&(info[ii].reqRecvInfo), &statusRecv));
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */
  //-----------------------------------------------------------------------
  setLETpartition(Ndomain, info);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* count total number of the roughly estimated LET nodes */
  /* if the estimated values are small, then enlarge values for safety */
  //-----------------------------------------------------------------------
  int numSend = 0;
  int numRecv = 0;
  for(int ii = 0; ii < Ndomain - 1; ii++){
    numSend += info[ii].maxSend;
    numRecv += info[ii].maxRecv;
  }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */
  //-----------------------------------------------------------------------
  __NOTE__("numSend = %d, numRecv = %d, numNode = %d @ rank %d\n", numSend, numRecv, numNode, mpi.rank);
  //-----------------------------------------------------------------------
#if 1
  const int upper = (int)floorf(0.9f * EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - ALIGN_BUF_FOR_LET(numNode);
  if( 2 * (numSend + numRecv) < upper ){
    //---------------------------------------------------------------------
    const float factor = floorf((float)upper / (float)(numSend + numRecv));
    //---------------------------------------------------------------------
    for(int ii = 0; ii < Ndomain - 1; ii++){
      if( info[ii].maxSend !=          numMax  )	info[ii].maxSend = ALIGN_BUF_FOR_LET((int)floorf(factor * (float)info[ii].maxSend));
      if( info[ii].maxRecv != info[ii].numNode )	info[ii].maxRecv = ALIGN_BUF_FOR_LET((int)floorf(factor * (float)info[ii].maxRecv));
    }/* for(int ii = 0; ii < Ndomain - 1; ii++){ */
    //---------------------------------------------------------------------
  }/* if( 2 * (numSend + numRecv) < upper ){ */
#endif
  //-----------------------------------------------------------------------
  if( info[0].maxSend + info[0].maxRecv > ((int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - ALIGN_BUF_FOR_LET(numNode)) ){
    __KILL__(stderr, "ERROR: rough estimation routine predicts # of LET nodes(%d) succeed the capacity of the array(%d - %d); communication with first partner would fail due to lack of MPI buffer.\n", info[0].maxSend + info[0].maxRecv, (int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE), ALIGN_BUF_FOR_LET(numNode));
  }/* if( numSend + numRecv > (int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - ALIGN_BUF_FOR_LET(numNode) ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
/* void guessLETpartition(const int Ndomain, domainInfo *info, const int numNode, const jparticle pj_root_hst, const real bmax2_root_hst, int2 * list, MPIcfg_tree mpi) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "start"); */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* guess required depth to get enough accuracy *\/ */
/*   //----------------------------------------------------------------------- */
/*   /\* send position, r2 of the root cell *\/ */
/*   const position icom = {pj_root_hst.x, pj_root_hst.y, pj_root_hst.z, bmax2_root_hst}; */
/*   //----------------------------------------------------------------------- */
/*   for(int ii = 0; ii < Ndomain - 1; ii++){ */
/*     //--------------------------------------------------------------------- */
/*     const real rx = info[ii].icom.x - icom.x; */
/*     const real ry = info[ii].icom.y - icom.y; */
/*     const real rz = info[ii].icom.z - icom.z; */
/*     const real r2 = rx * rx + ry * ry + rz * rz; */
/*     //--------------------------------------------------------------------- */
/*     const real lambda = UNITY - SQRTRATIO(icom.m, r2); */
/*     const real theta = SQRTRATIO(icom.m, lambda * lambda * r2); */
/*     //--------------------------------------------------------------------- */
/*     // */
/*     //--------------------------------------------------------------------- */
/*     /\* rough assumption: minimum value of theta is 0.03125 = 2^{-5} *\/ */
/*     /\* typical value of theta in cosmological simulation is ~0.1 *\/ */
/*     /\* here, we set # of steps required to theta becomes lower than 2^{-5} *\/ */
/*     //--------------------------------------------------------------------- */
/*     info[ii].numRecv = 5 + (int)CEIL(LOG2(theta)); */
/*     if( info[ii].numRecv <               0 )      info[ii].numRecv =                   2;/\* exchange tree data at least level 2 *\/ */
/*     if( info[ii].numRecv > NUM_PHKEY_LEVEL )      info[ii].numRecv = NUM_PHKEY_LEVEL - 1;/\* full tree data is required *\/ */
/*     //--------------------------------------------------------------------- */
/*   }/\* for(int ii = 0; ii < Ndomain - 1; ii++){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* exchange information on the estimated depth of LETs *\/ */
/*   //----------------------------------------------------------------------- */
/*   for(int ii = 0; ii < Ndomain - 1; ii++){ */
/*     //--------------------------------------------------------------------- */
/*     chkMPIerr(MPI_Isend(&(info[ii].numRecv), 1, MPI_INT, info[ii].rank,      mpi.rank, mpi.comm, &(info[ii].reqSendInfo))); */
/*     chkMPIerr(MPI_Irecv(&(info[ii].numSend), 1, MPI_INT, info[ii].rank, info[ii].rank, mpi.comm, &(info[ii].reqRecvInfo))); */
/*     //--------------------------------------------------------------------- */
/*   }/\* for(int ii = 0; ii < Ndomain - 1; ii++){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   for(int ii = 0; ii < Ndomain - 1; ii++){ */
/*     //--------------------------------------------------------------------- */
/*     MPI_Status statusSend;    chkMPIerr(MPI_Wait(&(info[ii].reqSendInfo), &statusSend)); */
/*     MPI_Status statusRecv;    chkMPIerr(MPI_Wait(&(info[ii].reqRecvInfo), &statusRecv)); */
/*     //--------------------------------------------------------------------- */
/*   }/\* for(int ii = 0; ii < Ndomain - 1; ii++){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   for(int ii = 0; ii < Ndomain - 1; ii++){ */
/*     //--------------------------------------------------------------------- */
/*     /\* estimate # of LET nodes *\/ */
/*     //--------------------------------------------------------------------- */
/*     int num = 0; */
/*     for(int jj = 0; jj < info[ii].numSend; jj++) */
/*       num += list[jj].y; */
/*     //--------------------------------------------------------------------- */
/* #if 1 */
/*     info[ii].numSend = increaseNumber(num); */
/* #else */
/*     info[ii].numSend = num; */
/* #endif */
/*     //--------------------------------------------------------------------- */
/*     chkMPIerr(MPI_Isend(&(info[ii].numSend), 1, MPI_INT, info[ii].rank,      mpi.rank, mpi.comm, &(info[ii].reqSendInfo))); */
/*     chkMPIerr(MPI_Irecv(&(info[ii].numRecv), 1, MPI_INT, info[ii].rank, info[ii].rank, mpi.comm, &(info[ii].reqRecvInfo))); */
/*     //--------------------------------------------------------------------- */
/*   }/\* for(int ii = 0; ii < Ndomain - 1; ii++){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   for(int ii = 0; ii < Ndomain - 1; ii++){ */
/*     //--------------------------------------------------------------------- */
/*     MPI_Status statusSend;    chkMPIerr(MPI_Wait(&(info[ii].reqSendInfo), &statusSend)); */
/*     MPI_Status statusRecv;    chkMPIerr(MPI_Wait(&(info[ii].reqRecvInfo), &statusRecv)); */
/*     //--------------------------------------------------------------------- */
/*   }/\* for(int ii = 0; ii < Ndomain - 1; ii++){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   setLETpartition(Ndomain, info); */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* count total number of the roughly estimated LET nodes *\/ */
/*   /\* if the estimated values are small, then enlarge values for safety *\/ */
/*   //----------------------------------------------------------------------- */
/*   int numSend = 0; */
/*   int numRecv = 0; */
/*   for(int ii = 0; ii < Ndomain - 1; ii++){ */
/*     numSend += info[ii].numSend; */
/*     numRecv += info[ii].numRecv; */
/*   } */
/*   //----------------------------------------------------------------------- */
/*   const int upper = (int)floorf(0.9f * EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - increaseNumber(numNode); */
/*   if( 2 * (numSend + numRecv) < upper ){ */
/*     //--------------------------------------------------------------------- */
/*     const float factor = floorf((float)upper / (float)(numSend + numRecv)); */
/*     //--------------------------------------------------------------------- */
/*     for(int ii = 0; ii < Ndomain - 1; ii++){ */
/* #if 0 */
/*       info[ii].numSend = (int)floorf(factor * (float)info[ii].numSend); */
/*       info[ii].numRecv = (int)floorf(factor * (float)info[ii].numRecv); */
/* #else */
/*       info[ii].numSend = increaseNumber((int)floorf(factor * (float)info[ii].numSend)); */
/*       info[ii].numRecv = increaseNumber((int)floorf(factor * (float)info[ii].numRecv)); */
/* #endif */
/*     }/\* for(int ii = 0; ii < Ndomain - 1; ii++){ *\/ */
/*     //--------------------------------------------------------------------- */
/*   }/\* if( 2 * (numSend + numRecv) < upper ){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   if( info[0].numSend + info[0].numRecv > ((int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - increaseNumber(numNode)) ){ */
/*     __KILL__(stderr, "ERROR: rough estimation routine predicts # of LET nodes(%d) succeed the capacity of the array(%d - %d); communication with first partner would fail due to lack of MPI buffer.\n", info[0].numSend + info[0].numRecv, (int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE), increaseNumber(numNode)); */
/*   }/\* if( numSend + numRecv > (int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - increaseNumber(numNode) ){ *\/ */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "end"); */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
/* /\* parallel prefix sum within a group of TGRP threads (TGRP <= 32 to use implicit synchronization) *\/ */
/* /\* type of prefix sum is inclusive *\/ */
/* //------------------------------------------------------------------------- */
/* __device__ __forceinline__ void prefixSumTsub(volatile uint_real * smem, const int tidx, const int lane) */
/* { */
/*   //----------------------------------------------------------------------- */
/* #ifndef BUILD_LET_ON_DEVICE */
/*   //----------------------------------------------------------------------- */
/* #pragma omp barrier */
/*   if( lane == 0 ){ */
/*     //--------------------------------------------------------------------- */
/*     int sum = 0; */
/*     for(int ii = 0; ii < TGRP; ii++){ */
/*       sum += smem[tidx + ii].i; */
/*       smem[tidx + ii].i = sum; */
/*     } */
/*     //--------------------------------------------------------------------- */
/*   } */
/* #pragma omp barrier */
/*   //----------------------------------------------------------------------- */
/* #else///BUILD_LET_ON_DEVICE */
/*   //----------------------------------------------------------------------- */
/* #if TGRP >=  2 */
/*   if( lane >=  1 )    smem[tidx].i += smem[tidx -  1].i; */
/* #if TGRP >=  4 */
/*   if( lane >=  2 )    smem[tidx].i += smem[tidx -  2].i; */
/* #if TGRP >=  8 */
/*   if( lane >=  4 )    smem[tidx].i += smem[tidx -  4].i; */
/* #if TGRP >= 16 */
/*   if( lane >=  8 )    smem[tidx].i += smem[tidx -  8].i; */
/* #if TGRP == 32 */
/*   if( lane >= 16 )    smem[tidx].i += smem[tidx - 16].i; */
/* #endif//TGRP == 32 */
/* #endif//TGRP >= 16 */
/* #endif//TGRP >=  8 */
/* #endif//TGRP >=  4 */
/* #endif//TGRP >=  2 */
/*   //----------------------------------------------------------------------- */
/* #endif//BUILD_LET_ON_DEVICE */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */


/* //------------------------------------------------------------------------- */
/* /\* make width-first LET (Locally Essential Tree) *\/ */
/* //------------------------------------------------------------------------- */
/* /\* icom     :: input          :: position and squared radius of a pseudo i-particle corresponding to N-body particles in a different domain *\/ */
/* /\* list     :: input          :: head index and number of elements for ``more'' belongs to the corresponding PH-key level *\/ */
/* /\* mask     :: input / output :: a flag indicating skip or read the corresponding tree node (0 or 1) *\/ */
/* /\* more_org :: input          :: head index and number of child particles of the corresponding j-particle (full tree data; i.e., local data) *\/ */
/* /\* jpos_org :: input          :: position and squared radius of pseudo N-body particle as j-particles (full tree data; i.e., local data) *\/ */
/* /\*   mj_org :: input          :: mass of pseudo N-body particle as j-particles (full tree data; i.e., local data) *\/ */
/* /\* more_let ::         output :: head index and number of child particles of the corresponding j-particle (subtracted tree data; i.e., LET) *\/ */
/* /\* jpos_let ::         output :: position and squared radius of pseudo N-body particle as j-particles (subtracted tree data; i.e., LET) *\/ */
/* /\*   mj_let ::         output :: mass of pseudo N-body particle as j-particles (subtracted tree data; i.e., LET) *\/ */
/* //------------------------------------------------------------------------- */
/* /\* assumption: maximum number of threads is TGRP <= 32 --> variables such as ``lane'' are not required *\/ */
/* /\* upgrade to multi thread version is possible if utilize prefix sum *\/ */
/* /\* LET generation for multiple node is possible if icom is array and multiple mask, more_let are set *\/ */
/* void makeLET */
/* (const position icom, */
/*  const int leafLevel, int2 * RESTRICT list, uint * RESTRICT mask, */
/*  const uint * RESTRICT more_org, const jparticle * RESTRICT jpos_org, const real * RESTRICT mj_org, */
/*  uint       * RESTRICT more_let,       jparticle * RESTRICT jpos_let,       real * RESTRICT mj_let, */
/*  int *numLETnode) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   /\* identify thread properties *\/ */
/*   //----------------------------------------------------------------------- */
/* #ifndef BUILD_LET_ON_DEVICE */
/*   const int tidx = 0; */
/* #else///BUILD_LET_ON_DEVICE */
/*   const int tidx = THREADIDX_X1D; */
/* #endif//BUILD_LET_ON_DEVICE */
/*   //----------------------------------------------------------------------- */
/*   /\* const uint head = tidx - lane; <-- always 0 *\/ */
/*   /\* const uint tail = head + (TGRP - 1); <-- always TGRP - 1 *\/ */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* shared quantities in the thread parallelized version *\/ */
/*   //----------------------------------------------------------------------- */
/*   /\* buffers to reduce the number of small packet data transfers *\/ */
/*   __shared__      uint more_buf[TGRP * (1 + NBUF_LET)]; */
/*   __shared__ jparticle jpos_buf[TGRP * (1 + NBUF_LET)]; */
/*   __shared__      real   mj_buf[TGRP * (1 + NBUF_LET)]; */
/*   __shared__      uint mask_buf[TGRP * (1 + NBUF_MSK)];/\* contain index of the mask array set to be unity *\/ */
/*   /\* //----------------------------------------------------------------------- *\/ */
/*   /\* /\\* head index of the shared buffer for more, jpos, and mj within a thread group *\\/ *\/ */
/*   /\* const uint hbuf = (head / TGRP) * TGRP * (1 + NBUF); *\/ */
/*   //----------------------------------------------------------------------- */
/*   /\* initialize buffers *\/ */
/* #pragma unroll */
/*   for(int ii = 0; ii < 1 + NBUF_LET; ii++) */
/*     more_buf[tidx + TGRP * ii] = NULL_NODE; */
/* #pragma unroll */
/*   for(int ii = 0; ii < 1 + NBUF_MSK; ii++) */
/*     mask_buf[tidx + TGRP * ii] = 0; */
/*   //----------------------------------------------------------------------- */
/*   int Nbuf_let = 0; */
/*   int Nbuf_msk = 0; */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* shared values within the threads *\/ */
/*   //----------------------------------------------------------------------- */
/*   /\* to store prefix sum *\/ */
/*   __shared__ uint_real smem[TGRP]; */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* total number of skipped tree nodes in upper and the same tree levels *\/ */
/*   int Nskip_upper = 0; */
/*   //----------------------------------------------------------------------- */
/*   int dstLetHead = 0; */
/*   int    LetHead = 0; */
/*   //----------------------------------------------------------------------- */
/*   for(int ii = 0; ii < leafLevel; ii++){ */
/*     //--------------------------------------------------------------------- */
/*     const int srcLetHead = list[ii].x; */
/*     const int srcLetNum  = list[ii].y; */
/*     //--------------------------------------------------------------------- */
/*     const int Niter = (srcLetNum + (TGRP - 1)) / TGRP; */
/*     //--------------------------------------------------------------------- */
/*     /\* number of skipped tree nodes in the lower tree level *\/ */
/*     int Nskip_lower = 0; */
/*     //--------------------------------------------------------------------- */
/*     // */
/*     //--------------------------------------------------------------------- */
/*     for(int jj = 0; jj < Niter; jj++){ */
/*       //------------------------------------------------------------------- */
/*       const int arrIdx = tidx + jj * TGRP; */
/*       int  read = 0; */
/*       uint more = NULL_NODE; */
/*       int  cnum = 0; */
/*       jparticle jpos; */
/*       real mj; */
/*       //------------------------------------------------------------------- */
/*       bool use = true; */
/*       if( arrIdx >= srcLetNum )	use = false; */
/*       //------------------------------------------------------------------- */
/*       // */
/*       //------------------------------------------------------------------- */
/*       /\* read mask value corresponding the original tree node *\/ */
/*       if( use ){ */
/* 	//----------------------------------------------------------------- */
/* 	read = mask    [srcLetHead + arrIdx]; */
/* 	more = more_org[srcLetHead + arrIdx]; */
/* 	cnum = 1 + (more >> IDXBITS); */
/* 	//----------------------------------------------------------------- */
/* 	/\* if the cell is a leaf cell in the full tree, then there is no child cell *\/ */
/* 	if( (more & IDXMASK) == (srcLetHead + arrIdx) ) */
/* 	  cnum = 0; */
/* 	//----------------------------------------------------------------- */
/*       } */
/*       //------------------------------------------------------------------- */
/*       int node = read; */
/*       smem[tidx].i = read; */
/*       prefixSumTsub(smem, tidx, tidx); */
/*       const int LetIdx = smem[tidx    ].i - read; */
/*       const int LetNum = smem[TGRP - 1].i; */
/*       //------------------------------------------------------------------- */
/*       // */
/*       //------------------------------------------------------------------- */
/*       /\* if the original tree node is copied to LET node *\/ */
/*       //------------------------------------------------------------------- */
/*       if( read ){ */
/* 	jpos = jpos_org[srcLetHead + arrIdx]; */
/* 	mj   =   mj_org[srcLetHead + arrIdx]; */
/* 	//----------------------------------------------------------------- */
/* 	// */
/* 	//----------------------------------------------------------------- */
/* 	/\* calculate distance between the pseudo i-particle (a representative particle for all particles in another another) and the candidate tree node *\/ */
/* 	//----------------------------------------------------------------- */
/* 	/\* set a pseudo i-particle *\/ */
/* 	const real rx = jpos.x - icom.x; */
/* 	const real ry = jpos.y - icom.y; */
/* 	const real rz = jpos.z - icom.z; */
/* 	const real r2 = rx * rx + ry * ry + rz * rz; */
/* 	real lambda = UNITY - SQRTRATIO(icom.m, r2);/\* icom.m is the square of the radius *\/ */
/* 	if( lambda < EPSILON )	  lambda = ZERO; */
/* 	/\* calculate distance between the pseudo i-particle and the candidate j-particle *\/ */
/* 	//----------------------------------------------------------------- */
/* #ifdef  WS93_MAC */
/* 	if(   jpos.w < lambda * lambda * r2 ) */
/* #else///WS93_MAC */
/* 	  if( jpos.w < lambda * lambda * r2 * theta2 ) */
/* #endif//WS93_MAC */
/* 	    { */
/* 	      //----------------------------------------------------------- */
/* 	      /\* distant node ==>> child cells are not included in the LET *\/ */
/* 	      //----------------------------------------------------------- */
/* 	      node = 0; */
/* 	      more = LetHead + LetIdx; */
/* 	      jpos.w = -UNITY;/\* squared size for the distant node is set to be negative *\/ */
/* 	      //----------------------------------------------------------- */
/* 	    } */
/* 	//----------------------------------------------------------------- */
/*       } */
/*       //------------------------------------------------------------------- */
/*       int Nskip = (!node) * cnum; */
/*       smem[tidx].i = Nskip; */
/*       prefixSumTsub(smem, tidx, tidx); */
/*       //------------------------------------------------------------------- */
/*       // */
/*       //------------------------------------------------------------------- */
/*       /\* commit LET node to the (small) shared memory *\/ */
/*       //------------------------------------------------------------------- */
/*       more_buf[Nbuf_let + LetIdx] = more - node * (Nskip_upper + Nskip_lower + smem[tidx].i); */
/*       jpos_buf[Nbuf_let + LetIdx] = jpos; */
/*       mj_buf  [Nbuf_let + LetIdx] =   mj; */
/*       Nbuf_let += LetNum; */
/*       //------------------------------------------------------------------- */
/*       Nskip_lower += smem[TGRP - 1].i; */
/*       LetHead += LetNum; */
/*       //------------------------------------------------------------------- */
/*       /\* write back to the main memory (CPU implementation), the global memory (GPU implementation), or communication with another process (FPGA implementation) *\/ */
/*       if( Nbuf_let >= TGRP * NBUF_LET ){ */
/* 	//----------------------------------------------------------------- */
/* #pragma unroll */
/* 	for(int ll = tidx; ll < TGRP * NBUF_LET; ll += TGRP){ */
/* 	  //--------------------------------------------------------------- */
/* 	  more_let[dstLetHead + ll] = more_buf[ll]; */
/* 	  jpos_let[dstLetHead + ll] = jpos_buf[ll]; */
/* 	  mj_let  [dstLetHead + ll] =   mj_buf[ll]; */
/* 	  //--------------------------------------------------------------- */
/* 	} */
/* 	//----------------------------------------------------------------- */
/* 	dstLetHead += TGRP * NBUF_LET; */
/* 	Nbuf_let   -= TGRP * NBUF_LET; */
/* 	//----------------------------------------------------------------- */
/*       } */
/*       //------------------------------------------------------------------- */
/*       // */
/*       // */
/*       //------------------------------------------------------------------- */
/*       /\* save mask values *\/ */
/*       //------------------------------------------------------------------- */
/*       int Nmask = node * cnum; */
/*       smem[tidx].i = Nmask; */
/*       prefixSumTsub(smem, tidx, tidx); */
/*       int mbuf_head = smem[tidx].i - Nmask; */
/*       const int mbuf_nadd = smem[TGRP - 1].i; */
/*       //------------------------------------------------------------------- */
/*       if( (Nbuf_msk + mbuf_nadd) < (TGRP * (1 + NBUF_MSK)) ){ */
/* 	//----------------------------------------------------------------- */
/* 	/\* if shared buffer for ``mask'' has enough space *\/ */
/* 	//----------------------------------------------------------------- */
/* 	for(int cc = 0; cc < Nmask; cc++) */
/* 	  mask_buf[Nbuf_msk + mbuf_head + cc] = (more & IDXMASK) + cc; */
/* 	Nbuf_msk += mbuf_nadd; */
/* 	//----------------------------------------------------------------- */
/*       } */
/*       else{ */
/* 	//----------------------------------------------------------------- */
/* 	/\* if shared buffer for ``mask'' does not have enough space *\/ */
/* 	//----------------------------------------------------------------- */
/* 	const int Nsave = (Nbuf_msk + mbuf_nadd) / (TGRP * (1 + NBUF_MSK)); */
/* 	const int Nrem  = (Nbuf_msk + mbuf_nadd) % (TGRP * (1 + NBUF_MSK)); */
/* 	for(int kk = 0; kk < Nsave; kk++){ */
/* 	  //--------------------------------------------------------------- */
/* 	  int Nopen = (TGRP * (1 + NBUF_MSK)) - Nbuf_msk; */
/* 	  //--------------------------------------------------------------- */
/* 	  /\* fill the shared buffer *\/ */
/* 	  if( mbuf_head < Nopen ){ */
/* 	    //------------------------------------------------------------- */
/* 	    Nopen -= mbuf_head; */
/* 	    const int Ncopy = (Nmask < Nopen) ? (Nmask) : (Nopen); */
/* 	    //------------------------------------------------------------- */
/* 	    /\* implicit condition: below statement is executed if Ncopy > 0 *\/ */
/* 	    for(int cc = 0; cc < Ncopy; cc++) */
/* 	      mask_buf[Nbuf_msk + mbuf_head + cc] = (more & IDXMASK) + cc; */
/* 	    //------------------------------------------------------------- */
/* 	    Nmask -= Ncopy; */
/* 	    more  += Ncopy; */
/* 	    //------------------------------------------------------------- */
/* 	  } */
/* 	  //--------------------------------------------------------------- */
/* 	  // */
/* 	  //--------------------------------------------------------------- */
/* 	  /\* send mask index from the shared buffer to the global memory *\/ */
/* 	  //--------------------------------------------------------------- */
/* #pragma unroll */
/* 	  for(int ll = tidx; ll < (TGRP * (1 + NBUF_MSK)); ll += TGRP) */
/* 	    mask[mask_buf[ll]] = 1; */
/* 	  //--------------------------------------------------------------- */
/* 	  Nbuf_msk   = 0; */
/* 	  mbuf_head -= IMIN(mbuf_head, Nopen); */
/* 	  //--------------------------------------------------------------- */
/* 	} */
/* 	//----------------------------------------------------------------- */
/* 	// */
/* 	//----------------------------------------------------------------- */
/* 	/\* store unsaved mask values to the shared buffer *\/ */
/* 	//----------------------------------------------------------------- */
/* 	/\* implicit condition: below statement is executed if Nmask > 0 *\/ */
/* 	for(int cc = 0; cc < Nmask; cc++) */
/* 	  mask_buf[mbuf_head + cc] = (more & IDXMASK) + cc; */
/* 	Nbuf_msk = Nrem; */
/* 	//----------------------------------------------------------------- */
/*       } */
/*       //------------------------------------------------------------------- */
/*     } */
/*     //--------------------------------------------------------------------- */
/*     // */
/*     //--------------------------------------------------------------------- */
/*     /\* move mask values of the lower level cells on the shared buffer to the global memory *\/ */
/*     //--------------------------------------------------------------------- */
/*     /\* implicit condition: below statement is executed if Nbuf_msk > 0 *\/ */
/*     for(int ll = tidx; ll < Nbuf_msk; ll += TGRP) */
/*       mask[mask_buf[ll]] = 1; */
/*     //--------------------------------------------------------------------- */
/*     Nbuf_msk = 0; */
/*     //--------------------------------------------------------------------- */

/*     //--------------------------------------------------------------------- */
/*     Nskip_upper += Nskip_lower; */
/*     //--------------------------------------------------------------------- */
/* #if 1 */
/*     if( ii < (leafLevel - 1) ) */
/*       if( Nskip_lower == list[ii + 1].y ) */
/* 	break; */
/* #endif */
/*     //--------------------------------------------------------------------- */
/*   } */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* write back to the main memory (CPU implementation), the global memory (GPU implementation), or communication with another process (FPGA implementation) *\/ */
/*   //----------------------------------------------------------------------- */
/*   /\* implicit condition: below statement is executed if Nbuf_let > 0 *\/ */
/*   for(int ll = tidx; ll < Nbuf_let; ll += TGRP){ */
/*     //--------------------------------------------------------------------- */
/*     more_let[dstLetHead + ll] = more_buf[ll]; */
/*     jpos_let[dstLetHead + ll] = jpos_buf[ll]; */
/*     mj_let  [dstLetHead + ll] =   mj_buf[ll]; */
/*     //--------------------------------------------------------------------- */
/*   } */
/*   //----------------------------------------------------------------------- */
/*   if( tidx == 0 ) */
/*     *numLETnode = dstLetHead + Nbuf_let; */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
