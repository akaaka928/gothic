/**
 * @file let_dev.cu
 *
 * @brief Source code for building locally essential tree (LET)
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/27 (Mon)
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
#include <helper_cuda.h>
#include <sys/time.h>

#include "macro.h"
#include "cudalib.h"
#include "timer.h"
#include "mpilib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/device.h"

#include "macutil.h"
#include "make.h"
#include "buf_inc.h"
#include "../para/mpicfg.h"
#include "let.h"
#include "walk_dev.h"
#include "let_dev.h"
#include "buf_inc.cu"


#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
__constant__  real theta2;
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)


/**
 * @fn configLETtopology
 *
 * @brief Memory allocation to configure MPI topology for exchanging LETs.
 */
muse configLETtopology(domainInfo **info, position **ipos,
#ifdef  GADGET_MAC
		       real **amin,
#endif//GADGET_MAC
		       int **numSend_hst, int **numSend_dev, cudaStream_t **stream, int *Nstream, const deviceProp gpu, MPIcfg_tree mpi)
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  *info = (domainInfo *)malloc(mpi.size * sizeof(domainInfo));  if( *info == NULL ){    __KILL__(stderr, "ERROR: failure to allocate info\n");  }
  alloc.host                += mpi.size * sizeof(domainInfo);
  *ipos = (position   *)malloc(mpi.size * sizeof(position  ));  if( *ipos == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ipos\n");  }
  alloc.host                += mpi.size * sizeof(position  );
#ifdef  GADGET_MAC
  *amin = (real       *)malloc(mpi.size * sizeof(real      ));  if( *amin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate amin\n");  }
  alloc.host                += mpi.size * sizeof(real      );
#endif//GADGET_MAC

  for(int ii = 0; ii < mpi.size - 1; ii++)
    (*info)[ii].rank = mpi.rank ^ (1 + ii);

  mycudaMalloc    ((void **)numSend_dev, mpi.size * sizeof(int));  alloc.device += mpi.size * sizeof(int);
  mycudaMallocHost((void **)numSend_hst, mpi.size * sizeof(int));  alloc.host   += mpi.size * sizeof(int);
  for(int ii = 0; ii < mpi.size; ii++){
    (*info)[ii].numSend_hst = &((*numSend_hst)[ii]);
    (*info)[ii].numSend_dev = &((*numSend_dev)[ii]);
  }/* for(int ii = 0; ii < mpi.size; ii++){ */

  *Nstream = NBLOCKS_PER_SM * gpu.numSM;
  *stream = (cudaStream_t *)malloc((*Nstream) * sizeof(cudaStream_t));
  alloc.host +=                    (*Nstream) * sizeof(cudaStream_t);
  if( *stream == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stream\n");  }
#   if  GPUVER >= 35
  int highest, lowest;
  checkCudaErrors(cudaDeviceGetStreamPriorityRange(&lowest, &highest));
  for(int ii = 0; ii < *Nstream; ii++)
    checkCudaErrors(cudaStreamCreateWithPriority(&((*stream)[ii]), cudaStreamDefault, highest));
    /* checkCudaErrors(cudaStreamCreateWithPriority(&((*stream)[ii]), cudaStreamNonBlocking, highest)); */
#else///GPUVER >= 35
#pragma unroll
  for(int ii = 0; ii < *Nstream; ii++)
    checkCudaErrors(cudaStreamCreate(&((*stream)[ii])));
#endif//GPUVER >= 35


  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn releaseLETtopology
 *
 * @brief Memory deallocation to configure MPI topology for exchanging LETs.
 */
void releaseLETtopology(domainInfo  *info, position  *ipos,
#ifdef  GADGET_MAC
			real  *amin,
#endif//GADGET_MAC
			int  *numSend_hst, int  *numSend_dev, cudaStream_t  *stream, int  Nstream)
{
  __NOTE__("%s\n", "start");

  free(info);
  free(ipos);
#ifdef  GADGET_MAC
  free(amin);
#endif//GADGET_MAC
  mycudaFree    (numSend_dev);
  mycudaFreeHost(numSend_hst);

  for(int ii = 0; ii < Nstream; ii++)
    mycudaStreamDestroy(stream[ii]);
  free(stream);

  __NOTE__("%s\n", "end");
}


#define NTHREADS_SCAN_INC NTHREADS_MAKE_LET
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_LET
#define USE_WARP_SHUFFLE_FUNC_SCAN_INC
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_LET
#include "../util/scan_inc.cu"


/**
 * @fn copyData_s2s
 *
 * @brief Move data from shared memory to shared memory.
 * @detail implicit synchronization within 32 threads in a warp is assumed
 */
__device__ __forceinline__ void copyData_s2s
(volatile uint * src, int sidx,
 volatile uint * dst, int didx, const int num, const int tidx)
{
  const int iter = DIV_NTHREADS_MAKE_LET(num);
  const int frac = num & (NTHREADS_MAKE_LET - 1);/**< Nload % NTHREADS_MAKE_LET */

  union {uint4 i; uint a[4];} tmp;
  for(int kk = 0; kk < (iter >> 2); kk++){
    /** load */
    tmp.i.x = src[sidx                        ];
    tmp.i.y = src[sidx +     NTHREADS_MAKE_LET];
    tmp.i.z = src[sidx + 2 * NTHREADS_MAKE_LET];
    tmp.i.w = src[sidx + 3 * NTHREADS_MAKE_LET];
    sidx += 4 * NTHREADS_MAKE_LET;
    __syncthreads();
    /** store; */
    dst[didx                        ] = tmp.i.x;
    dst[didx +     NTHREADS_MAKE_LET] = tmp.i.y;
    dst[didx + 2 * NTHREADS_MAKE_LET] = tmp.i.z;
    dst[didx + 3 * NTHREADS_MAKE_LET] = tmp.i.w;
    didx += 4 * NTHREADS_MAKE_LET;
  }/* for(int kk = 0; kk < (iter >> 2); kk++){ */

  const int loop = iter & 3;/**< 0, 1, 2, 3 */
  /** improved implementation, without Local memory */
  switch( loop ){
  case 0:    break;
  case 1:
    tmp.i.x = src[sidx];    sidx += NTHREADS_MAKE_LET;    __syncthreads();
    dst[didx] = tmp.i.x;    didx += NTHREADS_MAKE_LET;
    break;
  case 2:
    tmp.i.x = src[sidx];    tmp.i.y = src[sidx + NTHREADS_MAKE_LET];    sidx += 2 * NTHREADS_MAKE_LET;    __syncthreads();
    dst[didx] = tmp.i.x;    dst[didx + NTHREADS_MAKE_LET] = tmp.i.y;    didx += 2 * NTHREADS_MAKE_LET;
    break;
  case 3:
    tmp.i.x = src[sidx];    tmp.i.y = src[sidx + NTHREADS_MAKE_LET];    tmp.i.z = src[sidx + 2 * NTHREADS_MAKE_LET];    sidx += 3 * NTHREADS_MAKE_LET;    __syncthreads();
    dst[didx] = tmp.i.x;    dst[didx + NTHREADS_MAKE_LET] = tmp.i.y;    dst[didx + 2 * NTHREADS_MAKE_LET] = tmp.i.z;    didx += 3 * NTHREADS_MAKE_LET;
    break;
  }/* switch( loop ){ */

  if( frac > 0 ){
    if( tidx < frac )      tmp.i.x = src[sidx];
    __syncthreads();
    if( tidx < frac )      dst[didx] = tmp.i.x;
  }/* if( frac > 0 ){ */

  __syncthreads();
}

/**
 * @fn copyData_g2s
 *
 * @brief Move data from global memory to shared memory.
 */
__device__ __forceinline__ void copyData_g2s
(uint * RESTRICT gbuf, size_t srcHead, uint * RESTRICT sbuf, int dstHead, int numCopy, const int tidx)
{
  /** fraction processing at loading from the head of destination array */
  const int numTemp = NTHREADS_MAKE_LET - (int)(srcHead & (NTHREADS_MAKE_LET - 1));/**< NTHREADS_MAKE_LET - (srcHead % NTHREADS_MAKE_LET) */
  const int numHead = (numTemp < numCopy) ? numTemp : numCopy;
  if( tidx < numHead )
    sbuf[dstHead + tidx] = gbuf[srcHead + tidx];
  dstHead += numHead;
  srcHead += numHead;
  numCopy -= numHead;

  /** sequential load from source on the global memory and store to destination on the shared memory */
  for(int ii = tidx; ii < numCopy; ii += NTHREADS_MAKE_LET)
    sbuf[dstHead + ii] = gbuf[srcHead + ii];

  __syncthreads();
}

/**
 * @fn copyData_s2g
 *
 * @brief Move data from shared memory to global memory.
 */
__device__ __forceinline__ void copyData_s2g
(uint * RESTRICT sbuf, int srcHead, uint * RESTRICT gbuf, size_t dstHead, int numCopy, const int tidx)
{
  /** fraction processing at storing to the head of destination array */
  const int numTemp = NTHREADS_MAKE_LET - (int)(dstHead & (NTHREADS_MAKE_LET - 1));/**< NTHREADS_MAKE_LET - (dstHead % NTHREADS_MAKE_LET) */
  const int numHead = (numTemp < numCopy) ? numTemp : numCopy;
  if( tidx < numHead )
    gbuf[dstHead + tidx] = sbuf[srcHead + tidx];
  dstHead += numHead;
  srcHead += numHead;
  numCopy -= numHead;

  /** sequential load from source on the shared memory and store to destination on the global memory */
  for(int ii = tidx; ii < numCopy; ii += NTHREADS_MAKE_LET)
    gbuf[dstHead + ii] = sbuf[srcHead + ii];

  __syncthreads();
}

/**
 * @fn copyData_g2g
 *
 * @brief Move data from global memory to global memory.
 */
__device__ __forceinline__ void copyData_g2g
(uint * gbuf, size_t srcHead, size_t dstHead, int Ncopy, const int Ndisp, const int tidx)
{
  /** configure the settings */
  const int Nfirst = Ndisp & (NTHREADS_MAKE_LET - 1);/**< Ndisp % NTHREADS_MAKE_LET */
  const int  ldIdx = (tidx + Nfirst) & (NTHREADS_MAKE_LET - 1);  /**< ldIdx is Nfirst, Nfirst + 1, ..., NTHREADS_MAKE_LET - 1, 0, 1, ..., Nfirst - 1 for tidx of 0, 1, 2, ..., NTHREADS_MAKE_LET - 1 */
  const int grpIdx = (ldIdx < Nfirst) ? 0 : 1;

  srcHead += Ndisp - Nfirst;/**< hereafter, srcHead is NTHREADS_MAKE_LET elements aligned */

  /** fraction processing at loading from the head of source array */
  uint temp = gbuf[srcHead + ldIdx];
  srcHead += NTHREADS_MAKE_LET;

  /** sequential load and store from source to destination on the global memory */
  const int Niter = BLOCKSIZE(Ncopy, NTHREADS_MAKE_LET);
  for(int iter = 0; iter < Niter; iter += 4){
    const int Nmove = (Ncopy > (4 * NTHREADS_MAKE_LET)) ? (4 * NTHREADS_MAKE_LET) : (Ncopy);

    /** load from the source array on the global memory */
    /** load from temp (fraction processing) as initialization */
    union {uint4 i; uint a[4];} local;
    const int Nloop = BLOCKSIZE(Nmove, NTHREADS_MAKE_LET);/**< 1, 2, 3, 4 */

    /** improved implementation, without Local memory */
    switch( Nloop ){
    case 4:
      if(  grpIdx )	        local.i.x = temp;
      temp = gbuf[srcHead +  ldIdx                         ];
      if( !grpIdx )      	local.i.x = temp;
      else	                local.i.y = temp;
      temp = gbuf[srcHead + (ldIdx +     NTHREADS_MAKE_LET)];
      if( !grpIdx )      	local.i.y = temp;
      else	                local.i.z = temp;
      temp = gbuf[srcHead + (ldIdx + 2 * NTHREADS_MAKE_LET)];
      if( !grpIdx )      	local.i.z = temp;
      else	                local.i.w = temp;
      temp = gbuf[srcHead + (ldIdx + 3 * NTHREADS_MAKE_LET)];
      if( !grpIdx )      	local.i.w = temp;
      __syncthreads();
      gbuf[dstHead + (tidx			  )] = local.i.x;
      gbuf[dstHead + (tidx +	 NTHREADS_MAKE_LET)] = local.i.y;
      gbuf[dstHead + (tidx + 2 * NTHREADS_MAKE_LET)] = local.i.z;
      gbuf[dstHead + (tidx + 3 * NTHREADS_MAKE_LET)] = local.i.w;
      break;
    case 3:
      if(  grpIdx )      	local.i.x = temp;
      temp = gbuf[srcHead +  ldIdx                         ];
      if( !grpIdx )      	local.i.x = temp;
      else              	local.i.y = temp;
      temp = gbuf[srcHead + (ldIdx +     NTHREADS_MAKE_LET)];
      if( !grpIdx )      	local.i.y = temp;
      else              	local.i.z = temp;
      temp = gbuf[srcHead + (ldIdx + 2 * NTHREADS_MAKE_LET)];
      if( !grpIdx )      	local.i.z = temp;
      __syncthreads();
      gbuf[dstHead + (tidx			  )] = local.i.x;
      gbuf[dstHead + (tidx +	 NTHREADS_MAKE_LET)] = local.i.y;
      gbuf[dstHead + (tidx + 2 * NTHREADS_MAKE_LET)] = local.i.z;
      break;
    case 2:
      if(  grpIdx )      	local.i.x = temp;
      temp = gbuf[srcHead +  ldIdx                     ];
      if( !grpIdx )      	local.i.x = temp;
      else               	local.i.y = temp;
      temp = gbuf[srcHead + (ldIdx + NTHREADS_MAKE_LET)];
      if( !grpIdx )      	local.i.y = temp;
      __syncthreads();
      gbuf[dstHead + (tidx		      )] = local.i.x;
      gbuf[dstHead + (tidx + NTHREADS_MAKE_LET)] = local.i.y;
      break;
    case 1:
      if(  grpIdx )      	local.i.x = temp;
      temp = gbuf[srcHead + ldIdx];
      if( !grpIdx )      	local.i.x = temp;
      __syncthreads();
      gbuf[dstHead + tidx] = local.i.x;
      break;
    }/* switch( Nloop ){ */

    Ncopy   -= Nmove;
    srcHead += Nmove;
    dstHead += Nmove;
  }/* for(int iter = 0; iter < Niter; iter += 4){ */

  __syncthreads();
}


/**
 * @fn enqueueChildNodes
 *
 * @brief Enqueue child nodes.
 */
__device__ __forceinline__ void enqueueChildNodes
(const int tidx, const int lane, int * RESTRICT smem, const int leaf, const uint subNode,
 uint * RESTRICT smbuf,                  int *rem_sm, int *num_sm,
 uint * RESTRICT gmbuf, const size_t hb, int *rem_gm, int *num_gm, int *head_gm, int *tail_gm)
{
  /** 1. compact the given sparse tree nodes */
  int add = PREFIX_SUM_BLCK(leaf, smem, lane, tidx) - leaf;/**< exclusive prefix sum of leaf */
  int Ntot = smem[NTHREADS_MAKE_LET - 1];
  __syncthreads();
  smem[(leaf) ? (add) : (Ntot + tidx - add)] = (leaf) ? subNode : NULL_NODE;
  __syncthreads();

  /** 2. copy tree nodes to the shared memory */
  const int Nsm = (Ntot < *rem_sm) ? (Ntot) : (*rem_sm);
  copyData_s2s((uint *)smem, tidx, smbuf, tidx + (*num_sm), Nsm, tidx);
  *num_sm += Nsm;
  *rem_sm -= Nsm;
  Ntot    -= Nsm;

  /** 3. move tree nodes on the global memory, if necessary */
  if( Ntot > *rem_gm ){
    copyData_g2g(gmbuf, hb, hb, *num_gm, *head_gm, tidx);
    * rem_gm += *head_gm;
    *tail_gm -= *head_gm;
    *head_gm  = 0;
  }/* if( Ntot > *rem_gm ){ */

  /** 4. copy tree nodes to the global memory */
  copyData_s2g((uint *)smem, Nsm, gmbuf, hb + (*tail_gm), Ntot, tidx);
  * rem_gm -= Ntot;
  * num_gm += Ntot;
  *tail_gm += Ntot;
}


/**
 * @fn makeLET_kernel
 *
 * @brief Generate locally essential tree (LET) based on the width-first tree traversal.
 *
 * @param (icom) position and squared radius of a pseudo i-particle corresponding to N-body particles in a different domain
 * @return (numLETnode) the total number of LET nodes
 * @param (more_org) head index and number of child particles of the corresponding j-particle (full tree data; i.e., local data)
 * @param (jpos_org) position and squared radius of pseudo N-body particle as j-particles (full tree data; i.e., local data)
 * @param (  mj_org) mass of pseudo N-body particle as j-particles (full tree data; i.e., local data)
 * @return (more_let) head index and number of child particles of the corresponding j-particle (subtracted tree data; i.e., LET)
 * @return (jpos_let) position and squared radius of pseudo N-body particle as j-particles (subtracted tree data; i.e., LET)
 * @return (  mj_let) mass of pseudo N-body particle as j-particles (subtracted tree data; i.e., LET)
 * @param (active) a shared value to lock the shared quantities (freeNum, freeLst) to control usage of buffer
 * @param (freeNum) an unsigned integer represents # of unused bufferes
 * @param (freeLst) a list of unused bufferes
 * @param (buffer) tentative memory space to store tree cells which does not fit within the limited space of the shared memory
 * @param (bufSize) size of the buffer
 * @return (overflow) a variable to detect buffer overflow
 */
__global__ void __launch_bounds__(NTHREADS_MAKE_LET, NBLOCKS_PER_SM) makeLET_kernel
(READ_ONLY position icom,
#ifdef  GADGET_MAC
 READ_ONLY real amin,
#endif//GADGET_MAC
 int * RESTRICT numLETnode,
 READ_ONLY uint * RESTRICT more_org, READ_ONLY jparticle * RESTRICT jpos_org, READ_ONLY real * RESTRICT mj_org,
           uint * RESTRICT more_let,           jparticle * RESTRICT jpos_let,           real * RESTRICT mj_let,
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
 int * RESTRICT active, uint * RESTRICT freeNum,
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  !defined(USE_SMID_TO_GET_BUFID) &&  defined(TRY_MODE_ABOUT_BUFFER)
 const int freeNum,
#endif//!defined(USE_SMID_TO_GET_BUFID) &&  defined(TRY_MODE_ABOUT_BUFFER)
 uint * RESTRICT freeLst, uint * RESTRICT buffer, const int bufSize, int * RESTRICT overflow
#   if  !defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
 , unsigned long long int * RESTRICT cycles
#endif//!defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
)
{
  /** start stop watch */
#   if  !defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
  const long long int initCycle = clock64();
#endif//!defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)

  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  const int lane = tidx & (warpSize - 1);/**< index of the thread within a thread group */

  /** shared values within the threads */
  __shared__ uint queue[NTHREADS_MAKE_LET * NQUEUE_LET];
  __shared__  int  smem[NTHREADS_MAKE_LET];

#ifdef  USE_SMID_TO_GET_BUFID
  const int bufTarget = occupyBuffer(tidx, freeLst, queue);
#else///USE_SMID_TO_GET_BUFID
#ifdef  TRY_MODE_ABOUT_BUFFER
  const int bufTarget = occupyBuffer(tidx, BLOCKIDX_X1D, freeNum, freeLst, queue);
#else///TRY_MODE_ABOUT_BUFFER
  occupyBuffer(tidx, freeNum, freeLst, queue, active);
#endif//TRY_MODE_ABOUT_BUFFER
#endif//USE_SMID_TO_GET_BUFID
  const int bufIdx = (int)queue[0];
  __syncthreads();
  size_t buf0Head = (size_t)bufIdx * (size_t)bufSize;


  /** sweep all tree nodes by executing tree-traversal */
  /** initialize queue for tree nodes */
#pragma unroll
  for(int jj = 0; jj < NQUEUE_LET; jj++)
    queue[tidx + NTHREADS_MAKE_LET * jj] = NULL_NODE;/**< size >= NTHREADS_MAKE_LET * NQUEUE_LET */

  /** set child j-cells in queue on the shared memory */
  int rem = 1;
  uint jcell = 0;/**< 0 means that the head index is 0 and the number of tree nodes is 1; i.e., it is the root node */
  if( tidx == 0 )
    queue[0] = jcell;

  /** initialize queue for j-cells and interaction list by a representative thread */
  int sendNum = 0;/**< number of LET nodes already stored in the global memory */
  int  totNum = 1;/**< total number of nodes stored to the queue */
  int bufHead = 0;
  int bufTail = 0;
  int bufOpen = bufSize;
  int bufUsed = 0;


  /** tree traversal in a width-first manner */
  int fail = 0;
  while( true ){
    /** if the queue becomes empty, then exit the while loop */
    __syncthreads();
    if( rem == 0 )
      break;

    /** pick up a queue from stack */
    /** tentative load from the stack */
    int cnum = 0;
    jcell = NULL_NODE;
    if( tidx < rem ){
      jcell = queue[tidx];
      cnum = 1 + (int)(jcell >> IDXBITS);
    }/* if( lane < rem ){ */
    jcell &= IDXMASK;

    /** predict the head index on the shared memory by parallel prefix sum */
    int hidx = PREFIX_SUM_BLCK(cnum, smem, lane, tidx) - cnum;/**< exclusive prefix sum of cnum */

    smem[tidx] = NULL_NODE;
    __syncthreads();

    int remove = 0;
    if( (cnum != 0) && (hidx < NTHREADS_MAKE_LET) ){
      /** local data can be uploaded to the shared memory */
      int unum = NTHREADS_MAKE_LET - hidx;
      if( cnum < unum )	  unum = cnum;

      /** upload local data */
      for(int jj = 0; jj < unum; jj++){
	smem[hidx] = (int)jcell;/**< because hidx < NTHREADS_MAKE_LET */
	hidx++;
	jcell++;
      }/* for(int jj = 0; jj < unum; jj++){ */

      /** eliminate stocked j-cells from the queue */
      if( unum == cnum )
	remove = 1;
      else{
	jcell += ((uint)(cnum - unum - 1) << IDXBITS);
	queue[tidx] = jcell;
      }/* else{ */
    }/* if( (cnum != 0) && (hidx < NTHREADS_MAKE_LET) ){ */

    /** set an index of j-cell */
    __syncthreads();
    const int target = smem[tidx];

    /** remove scanned j-cells if possible */
    PREFIX_SUM_BLCK(remove, smem, lane, tidx);
    remove = smem[NTHREADS_MAKE_LET - 1];

    if( remove != 0 ){
      rem -= remove;
      copyData_s2s(queue, tidx + remove, queue, tidx, rem, tidx);
    }/* if( remove != 0 ){ */
    else
      __syncthreads();

    /** pick up pseudo particles */
    /** prefixSum to submit an LET node */
    int returnLET = (target != NULL_NODE) ? 1 : 0;
    hidx = sendNum + PREFIX_SUM_BLCK(returnLET, smem, lane, tidx) - returnLET;/** index of the corresponding LET node, which is based on exclusive prefix sum of calc */
    sendNum += smem[NTHREADS_MAKE_LET - 1];

    /** only the active threads pick up a j-cell from the global memory */
    jparticle jpos_tmp;
    uint      more_tmp;
    int childNum = 0;
    int hasChild = 0;
    if( returnLET ){
      jpos_tmp     = jpos_org[target];      /**< get position of pseudo j-particle */
      mj_let[hidx] =   mj_org[target];      /**< send mj of an LET node */

      /** set a pseudo i-particle */
      const real rx = jpos_tmp.x - icom.x;
      const real ry = jpos_tmp.y - icom.y;
      const real rz = jpos_tmp.z - icom.z;
      const real r2 = FLT_MIN + rx * rx + ry * ry + rz * rz;
      real lambda = FMAX(UNITY - SQRTRATIO(icom.m, r2), ZERO);

      /** calculate distance between the pseudo i-particle and the candidate j-particle */
      lambda *= lambda * r2;
#ifdef  GADGET_MAC
      /** alpha * |a| * r^4 > G * M * l^2 */
      if( jpos_tmp.w < lambda * lambda * amin )
#else///GADGET_MAC
#ifdef  WS93_MAC
	  if( jpos_tmp.w < lambda )
#else///WS93_MAC
	    /** (l / r) < theta */
	    if( jpos_tmp.w < lambda * theta2 )
#endif//WS93_MAC
#endif//GADGET_MAC
	      {
		/** distant node ==>> child cells are not included in the LET */
		more_tmp = hidx;
		jpos_tmp.w = -UNITY;/**< squared size for the distant node is set to be negative */
	      }
	    else{
	      /** near node ==> child cells are included in the LET */
	      /** add child-cells of near tree-cells to the tentative stack */
	      more_tmp = more_org[target];
	      childNum = 1 + (int)(more_tmp >> IDXBITS);
	      hasChild = 1;
	    }/* else{ */
    }/* if( returnLET ){ */


    /** if the shared memory has open space and some tree cells are stored on the global memory, then load tree-cells from the global memory to the shared memory */
    /** evaluate available size of the queue on the shared memory */
    int Nsm_rem = NQUEUE_LET * NTHREADS_MAKE_LET - rem;

    if( (bufUsed != 0) && (Nsm_rem > 0) ){
      const int Nload = (Nsm_rem < bufUsed) ? (Nsm_rem) : (bufUsed);
      copyData_g2s(buffer, buf0Head + bufHead, queue, rem, Nload, tidx);      /**< hq is tidx */
      rem     += Nload;
      Nsm_rem -= Nload;
      bufUsed -= Nload;
      bufHead += Nload;

      if( bufUsed == 0 ){
	bufHead = 0;
	bufTail = 0;
	bufOpen = bufSize;
      }/* if( bufUsed == 0 ){ */
    }/* if( (bufUsed != 0) && (Nsm_rem > 0) ){ */
    else
      __syncthreads();

    /** copy child-cells of near tree-cells stored in the tentative stack to the stack on the shared memory and/or the global memory */
    enqueueChildNodes(tidx, lane, smem, hasChild, more_tmp, queue, &Nsm_rem, &rem, buffer, buf0Head, &bufOpen, &bufUsed, &bufHead, &bufTail);
    fail += (bufTail > bufSize);

    /** if current node has child nodes in LET, then head index of more_tmp must be rewritten */
    /** prefixSum to extend LET */
    int leafHead = PREFIX_SUM_BLCK(childNum, smem, lane, tidx) - childNum;/**< exclusive prefix sum of nchild */
    /** modify more pointer using leafHead */
    if( childNum > 0 )
      more_tmp = ((uint)(childNum - 1) << IDXBITS) + (uint)(totNum + leafHead);
    totNum += smem[NTHREADS_MAKE_LET - 1];

    /** add tree nodes to LET (mj_tmp is already stored) */
    if( returnLET ){
      jpos_let[hidx] = jpos_tmp;
      more_let[hidx] = more_tmp;
    }/* if( returnLET ){ */
  }/* while( true ){ */


  /* finalizing LET generator */
  if( tidx == 0 ){
    *numLETnode = sendNum;
    atomicAdd(overflow, fail);
  }/* if( tidx == 0 ){ */

#ifdef  USE_SMID_TO_GET_BUFID
  releaseBuffer(tidx, freeLst, (uint)bufIdx, bufTarget);
#else///USE_SMID_TO_GET_BUFID
#ifdef  TRY_MODE_ABOUT_BUFFER
  releaseBuffer(tidx, freeLst, (uint)bufIdx, bufTarget);
#else///TRY_MODE_ABOUT_BUFFER
  releaseBuffer(tidx, freeNum, freeLst, bufIdx, active);
#endif//TRY_MODE_ABOUT_BUFFER
#endif//USE_SMID_TO_GET_BUFID

#   if  !defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
  long long int exitCycle = clock64();
  if( tidx == 0 ){
    unsigned long long int elapsed = (unsigned long long int)(exitCycle - initCycle);
    atomicAdd(cycles, elapsed);
  }/* if( tidx == 0 ){ */
#endif//!defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
}


/**
 * @fn callGenLET
 *
 * @brief Generate locally essential tree (LET) based on the width-first tree traversal.
 */
extern "C"
void callGenLET
  (const cudaStream_t stream, domainInfo *let, MPIcfg_tree mpi, const soaTreeNode tree, const soaTreeWalkBuf buf
#ifdef  MONITOR_LETGEN_TIME
#ifdef  USE_CUDA_EVENT
   , const cudaEvent_t iniEvent, const cudaEvent_t finEvent
#else///USE_CUDA_EVENT
   , unsigned long long int * RESTRICT cycles
#endif//USE_CUDA_EVENT
#endif//MONITOR_LETGEN_TIME
   )
{
  __NOTE__("%s\n", "start");


#   if  defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
  checkCudaErrors(cudaEventRecord(iniEvent, 0));
#endif//defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)

  makeLET_kernel<<<1, NTHREADS_MAKE_LET, SMEM_SIZE, stream>>>
    ((*let).icom,
#ifdef  GADGET_MAC
     (*let).amin,
#endif//GADGET_MAC
     (*let).numSend_dev,
     tree.more, tree.jpos, tree.mj,
     &(tree.more[(*let).headSend]), &(tree.jpos[(*let).headSend]), &(tree.mj[(*let).headSend]),
#ifndef USE_SMID_TO_GET_BUFID
#ifndef TRY_MODE_ABOUT_BUFFER
     buf.active,
#endif//TRY_MODE_ABOUT_BUFFER
     buf.freeNum,
#endif//USE_SMID_TO_GET_BUFID
     buf.freeLst, buf.buffer, NGROUPS * buf.bufSize, buf.fail
#   if  !defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
     , cycles
#endif//!defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
     );

#   if  defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
  checkCudaErrors(cudaEventRecord(finEvent, 0));
#endif//defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)


  __NOTE__("%s\n", "end");
}


/**
 * @fn setGlobalConstants_let_dev_cu
 *
 * @brief Set global constants for let_dev.cu and initialize kernel functions.
 */
extern "C"
void setGlobalConstants_let_dev_cu
(
#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
 const real theta2_hst
#else///endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
 void
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
 )
{
  __NOTE__("%s\n", "start");


#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
#   if  CUDART_VERSION >= 5000
  cudaMemcpyToSymbol( theta2 , &theta2_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#else//CUDART_VERSION >= 5000
  cudaMemcpyToSymbol("theta2", &theta2_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#endif//CUDART_VERSION >= 5000
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)

#   if  SMPREF_LET == 1
  checkCudaErrors(cudaFuncSetCacheConfig(makeLET_kernel, cudaFuncCachePreferShared));
#endif//SMPREF_LET == 1


  __NOTE__("%s\n", "end");
}
