/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/22(Fri) 10:48:19
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
#include <helper_cuda.h>
#include <sys/time.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
#include <timer.h>
#include <mpilib.h>
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/device.h"
//-------------------------------------------------------------------------
#include "macutil.h"
#include "make.h"
#include "buf_inc.h"
//-------------------------------------------------------------------------
#include "../para/mpicfg.h"
#include "let.h"
#include "walk_dev.h"
#include "let_dev.h"
//-------------------------------------------------------------------------
#include "buf_inc.cu"
//-------------------------------------------------------------------------
#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
__constant__  real theta2;
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
muse configLETtopology(domainInfo **info, position **ipos,
#ifdef  GADGET_MAC
		       real **amin,
#endif//GADGET_MAC
#ifdef  BUILD_LET_ON_DEVICE
		       int **numSend_hst, int **numSend_dev,
#endif//BUILD_LET_ON_DEVICE
		       cudaStream_t **stream, int *Nstream, const deviceProp gpu, MPIcfg_tree mpi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  *info = (domainInfo *)malloc(mpi.size * sizeof(domainInfo));  if( *info == NULL ){    __KILL__(stderr, "ERROR: failure to allocate info\n");  }
  alloc.host                += mpi.size * sizeof(domainInfo);
  *ipos = (position   *)malloc(mpi.size * sizeof(position  ));  if( *ipos == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ipos\n");  }
  alloc.host                += mpi.size * sizeof(position  );
#ifdef  GADGET_MAC
  *amin = (real       *)malloc(mpi.size * sizeof(real      ));  if( *amin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate amin\n");  }
  alloc.host                += mpi.size * sizeof(real      );
#endif//GADGET_MAC
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < mpi.size - 1; ii++)
    (*info)[ii].rank = mpi.rank ^ (1 + ii);
  //-----------------------------------------------------------------------
#ifdef  BUILD_LET_ON_DEVICE
  mycudaMalloc    ((void **)numSend_dev, mpi.size * sizeof(int));  alloc.device += mpi.size * sizeof(int);
  mycudaMallocHost((void **)numSend_hst, mpi.size * sizeof(int));  alloc.host   += mpi.size * sizeof(int);
  for(int ii = 0; ii < mpi.size; ii++){
    (*info)[ii].numSend_hst = &((*numSend_hst)[ii]);
    (*info)[ii].numSend_dev = &((*numSend_dev)[ii]);
  }/* for(int ii = 0; ii < mpi.size; ii++){ */
#endif//BUILD_LET_ON_DEVICE
  //-----------------------------------------------------------------------
  *Nstream = NBLOCKS_PER_SM * gpu.numSM;
  *stream = (cudaStream_t *)malloc((*Nstream) * sizeof(cudaStream_t));
  alloc.host +=                    (*Nstream) * sizeof(cudaStream_t);
  if( *stream == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stream");  }
#pragma unroll
  for(int ii = 0; ii < *Nstream; ii++)
    checkCudaErrors(cudaStreamCreate(&((*stream)[ii])));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void releaseLETtopology(domainInfo  *info, position  *ipos,
#ifdef  GADGET_MAC
			real  *amin,
#endif//GADGET_MAC
#ifdef  BUILD_LET_ON_DEVICE
			int  *numSend_hst, int  *numSend_dev,
#endif//BUILD_LET_ON_DEVICE
			cudaStream_t  *stream, int  Nstream
			)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  free(info);
  free(ipos);
#ifdef  GADGET_MAC
  free(amin);
#endif//GADGET_MAC
#ifdef  BUILD_LET_ON_DEVICE
  mycudaFree    (numSend_dev);
  mycudaFreeHost(numSend_hst);
#endif//BUILD_LET_ON_DEVICE
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Nstream; ii++)
    mycudaStreamDestroy(stream[ii]);
  free(stream);
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* parallel prefix sum within a block */
/* type of prefix sum is inclusive */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
__device__ __forceinline__ int prefixSum(int val, const int tidx, const int lane, volatile int * smem)
{
  //-----------------------------------------------------------------------
  /* 1. prefix sum within a warp */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_LET
  //-----------------------------------------------------------------------
  /* load index */
  int tmp;
  /* calculate inclusive prefix sum */
  tmp = __shfl_up(val,  1, warpSize);  if( lane >=  1 )    val += tmp;
  tmp = __shfl_up(val,  2, warpSize);  if( lane >=  2 )    val += tmp;
  tmp = __shfl_up(val,  4, warpSize);  if( lane >=  4 )    val += tmp;
  tmp = __shfl_up(val,  8, warpSize);  if( lane >=  8 )    val += tmp;
  tmp = __shfl_up(val, 16, warpSize);  if( lane >= 16 )    val += tmp;
  /* return calculated inclusive prefix sum */
  smem[tidx] = val;
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_MAKE_LET
  //-----------------------------------------------------------------------
  smem[tidx] = val;
  if( lane >=  1 ){    val += smem[tidx -  1];    smem[tidx] = val;  }
  if( lane >=  2 ){    val += smem[tidx -  2];    smem[tidx] = val;  }
  if( lane >=  4 ){    val += smem[tidx -  4];    smem[tidx] = val;  }
  if( lane >=  8 ){    val += smem[tidx -  8];    smem[tidx] = val;  }
  if( lane >= 16 ){    val += smem[tidx - 16];    smem[tidx] = val;  }
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_LET
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#   if  NTHREADS_MAKE_LET >= 64
  //-----------------------------------------------------------------------
  /* 2. prefix sum about the tail of each warp */
  //-----------------------------------------------------------------------
  int scan = val;
  __syncthreads();
  /* warpSize = 32 = 2^5 */
  if( tidx < (NTHREADS_MAKE_LET >> 5) ){
    //---------------------------------------------------------------------
    val = smem[tidx * warpSize + warpSize - 1];
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_LET
#   if  NTHREADS_MAKE_LET >=   64
    const int groupSize = NTHREADS_MAKE_LET >> 5;
    tmp = __shfl_up(val,  1, groupSize);    if( lane >=  1 )      val += tmp;
#   if  NTHREADS_MAKE_LET >=  128
    tmp = __shfl_up(val,  2, groupSize);    if( lane >=  2 )      val += tmp;
#   if  NTHREADS_MAKE_LET >=  256
    tmp = __shfl_up(val,  4, groupSize);    if( lane >=  4 )      val += tmp;
#   if  NTHREADS_MAKE_LET >=  512
    tmp = __shfl_up(val,  8, groupSize);    if( lane >=  8 )      val += tmp;
#   if  NTHREADS_MAKE_LET == 1024
    tmp = __shfl_up(val, 16, groupSize);    if( lane >= 16 )      val += tmp;
#endif//NTHREADS_MAKE_LET == 1024
#endif//NTHREADS_MAKE_LET >=  512
#endif//NTHREADS_MAKE_LET >=  256
#endif//NTHREADS_MAKE_LET >=  128
#endif//NTHREADS_MAKE_LET >=   64
    smem[tidx] = val;
#else///USE_WARP_SHUFFLE_FUNC_MAKE_LET
    smem[tidx] = val;
#   if  NTHREADS_MAKE_LET >=   64
    if( lane >=  1 )      smem[tidx] += smem[tidx -  1];
#   if  NTHREADS_MAKE_LET >=  128
    if( lane >=  2 )      smem[tidx] += smem[tidx -  2];
#   if  NTHREADS_MAKE_LET >=  256
    if( lane >=  4 )      smem[tidx] += smem[tidx -  4];
#   if  NTHREADS_MAKE_LET >=  512
    if( lane >=  8 )      smem[tidx] += smem[tidx -  8];
#   if  NTHREADS_MAKE_LET == 1024
    if( lane >= 16 )      smem[tidx] += smem[tidx - 16];
#endif//NTHREADS_MAKE_LET == 1024
#endif//NTHREADS_MAKE_LET >=  512
#endif//NTHREADS_MAKE_LET >=  256
#endif//NTHREADS_MAKE_LET >=  128
#endif//NTHREADS_MAKE_LET >=   64
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_LET
    //---------------------------------------------------------------------
  }/* if( tidx < (NTHREADS_MAKE_LET >> 5) ){ */
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 3. prefix sum within a block */
  //-----------------------------------------------------------------------
  /* warpSize = 32 = 2^5 */
  if( tidx >= warpSize )
    scan += smem[(tidx >> 5) - 1];
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 4. upload calculate prefix sum */
  //-----------------------------------------------------------------------
  smem[tidx] = scan;
  val = scan;
  __syncthreads();
  //-----------------------------------------------------------------------
#endif//NTHREADS_MAKE_LET >= 64
  //-----------------------------------------------------------------------
  return (val);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_s2s
(volatile uint * src, int sidx,
 volatile uint * dst, int didx, const int num, const int tidx)
{
  //-----------------------------------------------------------------------
  const int iter = DIV_NTHREADS_MAKE_LET(num);
  const int frac = num & (NTHREADS_MAKE_LET - 1);/* := Nload % NTHREADS_MAKE_LET */
  //-----------------------------------------------------------------------
  union {uint4 i; uint a[4];} tmp;
  for(int kk = 0; kk < (iter >> 2); kk++){
    /* load */
    tmp.i.x = src[sidx                        ];
    tmp.i.y = src[sidx +     NTHREADS_MAKE_LET];
    tmp.i.z = src[sidx + 2 * NTHREADS_MAKE_LET];
    tmp.i.w = src[sidx + 3 * NTHREADS_MAKE_LET];
    sidx += 4 * NTHREADS_MAKE_LET;
    __syncthreads();
    /* store; */
    dst[didx                        ] = tmp.i.x;
    dst[didx +     NTHREADS_MAKE_LET] = tmp.i.y;
    dst[didx + 2 * NTHREADS_MAKE_LET] = tmp.i.z;
    dst[didx + 3 * NTHREADS_MAKE_LET] = tmp.i.w;
    didx += 4 * NTHREADS_MAKE_LET;
  }/* for(int kk = 0; kk < (iter >> 2); kk++){ */
  //-----------------------------------------------------------------------
  const int loop = iter & 3;
  /* load */
#pragma unroll
  for(int ii = 0; ii < loop; ii++){
    tmp.a[ii] = src[sidx];
    sidx += NTHREADS_MAKE_LET;
  }/* for(int ii = 0; ii < loop; ii++){ */
  if( loop != 0 )
    __syncthreads();
  /* store; */
#pragma unroll
  for(int ii = 0; ii < loop; ii++){
    dst[didx] = tmp.a[ii];
    didx += NTHREADS_MAKE_LET;
  }/* for(int ii = 0; ii < loop; ii++){ */
  //-----------------------------------------------------------------------
  if( frac > 0 ){
    if( tidx < frac )      tmp.i.x = src[sidx];
    __syncthreads();
    if( tidx < frac )      dst[didx] = tmp.i.x;
  }/* if( frac > 0 ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_g2s
(uint * RESTRICT gbuf, int srcHead, uint * RESTRICT sbuf, int dstHead, int numCopy, const int tidx)
{
  //-----------------------------------------------------------------------
  /* fraction processing at loading from the head of destination array */
  //-----------------------------------------------------------------------
  const int numTemp = NTHREADS_MAKE_LET - (srcHead & (NTHREADS_MAKE_LET - 1));/* := NTHREADS_MAKE_LET - (srcHead % NTHREADS_MAKE_LET) */
  const int numHead = (numTemp < numCopy) ? numTemp : numCopy;
  if( tidx < numHead )
    sbuf[dstHead + tidx] = gbuf[srcHead + tidx];
  dstHead += numHead;
  srcHead += numHead;
  numCopy -= numHead;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* sequential load from source on the global memory and store to destination on the shared memory */
  //-----------------------------------------------------------------------
  for(int ii = tidx; ii < numCopy; ii += NTHREADS_MAKE_LET)
    sbuf[dstHead + ii] = gbuf[srcHead + ii];
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_s2g
(uint * RESTRICT sbuf, int srcHead, uint * RESTRICT gbuf, int dstHead, int numCopy, const int tidx)
{
  //-----------------------------------------------------------------------
  /* fraction processing at storing to the head of destination array */
  //-----------------------------------------------------------------------
  const int numTemp = NTHREADS_MAKE_LET - (dstHead & (NTHREADS_MAKE_LET - 1));/* := NTHREADS_MAKE_LET - (dstHead % NTHREADS_MAKE_LET) */
  const int numHead = (numTemp < numCopy) ? numTemp : numCopy;
  if( tidx < numHead )
    gbuf[dstHead + tidx] = sbuf[srcHead + tidx];
  dstHead += numHead;
  srcHead += numHead;
  numCopy -= numHead;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* sequential load from source on the shared memory and store to destination on the global memory */
  //-----------------------------------------------------------------------
  for(int ii = tidx; ii < numCopy; ii += NTHREADS_MAKE_LET)
    gbuf[dstHead + ii] = sbuf[srcHead + ii];
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_g2g
(uint * gbuf, int srcHead, int dstHead, int Ncopy, const int Ndisp, const int tidx)
{
  //-----------------------------------------------------------------------
  /* configure the settings */
  //-----------------------------------------------------------------------
  const int Nfirst = Ndisp & (NTHREADS_MAKE_LET - 1);/* := Ndisp % NTHREADS_MAKE_LET */
  /* ldIdx is Nfirst, Nfirst + 1, ..., NTHREADS_MAKE_LET - 1, 0, 1, ..., Nfirst - 1 for tidx of 0, 1, 2, ..., NTHREADS_MAKE_LET - 1 */
  const int  ldIdx = (tidx + Nfirst) & (NTHREADS_MAKE_LET - 1);/* := (tidx + Nfirst) % NTHREADS_MAKE_LET */
  const int grpIdx = (ldIdx < Nfirst) ? 0 : 1;
  //-----------------------------------------------------------------------
  srcHead += Ndisp - Nfirst;/* hereafter, srcHead is NTHREADS_MAKE_LET elements aligned */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* fraction processing at loading from the head of source array */
  //-----------------------------------------------------------------------
  uint temp = gbuf[srcHead + ldIdx];
  srcHead += NTHREADS_MAKE_LET;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* sequential load and store from source to destination on the global memory */
  //-----------------------------------------------------------------------
  const int Niter = BLOCKSIZE(Ncopy, NTHREADS_MAKE_LET);
  /* for(int iter = 0; iter < Niter; iter++){ */
  for(int iter = 0; iter < Niter; iter += 4){
    //---------------------------------------------------------------------
    /* const int Nmove = (Ncopy > NTHREADS_MAKE_LET) ? (NTHREADS_MAKE_LET) : (Ncopy); */
    const int Nmove = (Ncopy > (4 * NTHREADS_MAKE_LET)) ? (4 * NTHREADS_MAKE_LET) : (Ncopy);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* load from the source array on the global memory */
    //---------------------------------------------------------------------
    /* load from temp (fraction processing) as initialization */
    /* uint local = temp; */
    union {uint4 i; uint a[4];} local;
    //---------------------------------------------------------------------
    /* load from global memory, store to shared memory or temp (fraction processing) */
    /* temp = gbuf[srcHead + ldIdx]; */
    /* if( !grpIdx ) */
    /*   local = temp; */
    const int Nloop = BLOCKSIZE(Nmove, NTHREADS_MAKE_LET);
#pragma unroll
    for(int ii = 0; ii < Nloop; ii++){
      if(  grpIdx )      	local.a[ii] = temp;
      temp = gbuf[srcHead + ldIdx + ii * NTHREADS_MAKE_LET];
      if( !grpIdx )      	local.a[ii] = temp;
    }/* for(int ii = 0; ii < Nloop; ii++){ */
    //---------------------------------------------------------------------
    __syncthreads();
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* store to the destination array on the global memory */
    //---------------------------------------------------------------------
    /* gbuf[dstHead + tidx] = local; */
#pragma unroll
    for(int ii = 0; ii < Nloop; ii++)
      gbuf[dstHead + tidx + ii * NTHREADS_MAKE_LET] = local.a[ii];
    //---------------------------------------------------------------------
    Ncopy   -= Nmove;
    srcHead += Nmove;
    dstHead += Nmove;
    //---------------------------------------------------------------------
  }/* for(int iter = 0; iter < Niter; iter += 4){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* uint smem[NTHREADS_MAKE_LET]; */
/* uint node[NTHREADS_MAKE_LET]; */
/* sizes of smem and node are the same, do not used at the same time ==>> use smem as node */
//-------------------------------------------------------------------------
__device__ __forceinline__ void enqueueChildNodes
(const int tidx, const int lane, int * RESTRICT smem, const int leaf, const uint subNode,
 uint * RESTRICT smbuf,                  int *rem_sm, int *num_sm,
 uint * RESTRICT gmbuf, const size_t hb, int *rem_gm, int *num_gm, int *head_gm, int *tail_gm
)
{
  //-----------------------------------------------------------------------
  /* 1. compact the given sparse tree nodes */
  //-----------------------------------------------------------------------
  int add = prefixSum(leaf, tidx, lane, smem) - leaf;/* exclusive prefix sum of leaf */
  //-----------------------------------------------------------------------
  __syncthreads();
  int Ntot = smem[NTHREADS_MAKE_LET - 1];
  /* node[(leaf) ? (add) : (Ntot + tidx - add)] = (leaf) ? subNode : NULL_NODE; */
  smem[(leaf) ? (add) : (Ntot + tidx - add)] = (leaf) ? subNode : NULL_NODE;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 2. copy tree nodes to the shared memory */
  //-----------------------------------------------------------------------
  const int Nsm = (Ntot < *rem_sm) ? (Ntot) : (*rem_sm);
  /* copyData_s2s(node, tidx, smbuf, tidx + (*num_sm), Nsm, tidx); */
  copyData_s2s((uint *)smem, tidx, smbuf, tidx + (*num_sm), Nsm, tidx);
  //-----------------------------------------------------------------------
  *num_sm += Nsm;
  *rem_sm -= Nsm;
  Ntot    -= Nsm;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 3. move tree nodes on the global memory, if necessary */
  //-----------------------------------------------------------------------
  if( Ntot > *rem_gm ){
    //---------------------------------------------------------------------
    copyData_g2g(gmbuf, hb, hb, *num_gm, *head_gm, tidx);
/* #pragma unroll */
/*     for(int ii = tidx; ii < *num_gm; ii += NTHREADS_MAKE_LET) */
/*       gmbuf[hb + (size_t)ii] = gmbuf[hb + (size_t)((*head_gm) + ii)]; */
/*       /\* gmbuf[hb + (size_t)((*head_gm) + ii)] = gmbuf[hb + (size_t)ii]; *\/ */
    //---------------------------------------------------------------------
    * rem_gm += *head_gm;
    *tail_gm -= *head_gm;
    *head_gm  = 0;
    //---------------------------------------------------------------------
  }/* if( Ntot > *rem_gm ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 4. copy tree nodes to the global memory */
  //-----------------------------------------------------------------------
#if 0
  memcpy((void *)&gmbuf[hb + (size_t)(*tail_gm)], (const void *)&smem[Nsm], sizeof(uint) * Ntot);
#else
  /* copyData_s2g(node, Nsm, gmbuf, hb + (size_t)(*tail_gm), Ntot, tidx); */
  copyData_s2g((uint *)smem, Nsm, gmbuf, hb + (size_t)(*tail_gm), Ntot, tidx);
#endif
  //-----------------------------------------------------------------------
  * rem_gm -= Ntot;
  * num_gm += Ntot;
  *tail_gm += Ntot;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* make width-first LET (Locally Essential Tree) */
//-------------------------------------------------------------------------
/* icom       :: input          :: position and squared radius of a pseudo i-particle corresponding to N-body particles in a different domain */
/* numLETnode ::         output :: the total number of LET nodes */
/* more_org   :: input          :: head index and number of child particles of the corresponding j-particle (full tree data; i.e., local data) */
/* jpos_org   :: input          :: position and squared radius of pseudo N-body particle as j-particles (full tree data; i.e., local data) */
/*   mj_org   :: input          :: mass of pseudo N-body particle as j-particles (full tree data; i.e., local data) */
/* more_let   ::         output :: head index and number of child particles of the corresponding j-particle (subtracted tree data; i.e., LET) */
/* jpos_let   ::         output :: position and squared radius of pseudo N-body particle as j-particles (subtracted tree data; i.e., LET) */
/*   mj_let   ::         output :: mass of pseudo N-body particle as j-particles (subtracted tree data; i.e., LET) */
/* active     ::                :: a shared value to lock the shared quantities (freeNum, freeLst) to control usage of buffer */
/* freeNum    ::                :: an unsigned integer represents # of unused bufferes */
/* freeLst    ::                :: a list of unused bufferes */
/* buffer     ::                :: tentative memory space to store tree cells which does not fit within the limited space of the shared memory */
/* bufSize    :: input          :: size of the buffer */
/* overflow   ::         output :: a variable to detect buffer overflow */
//-------------------------------------------------------------------------
__global__ void makeLET_kernel
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
#ifdef  MONITOR_LETGEN_TIME
 , unsigned long long int * RESTRICT cycles
#endif//MONITOR_LETGEN_TIME
)
{
  //-----------------------------------------------------------------------
  /* start stop watch */
  //-----------------------------------------------------------------------
#ifdef  MONITOR_LETGEN_TIME
  const long long int initCycle = clock64();
#endif//MONITOR_LETGEN_TIME
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int lane = tidx & (warpSize - 1);/* index of the thread within a thread group */
  //-----------------------------------------------------------------------
  /* const int head = tidx - lane; */
  /* const int tail = head + (warpSize - 1); */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* shared values within the threads */
  //-----------------------------------------------------------------------
  __shared__ uint queue[NTHREADS_MAKE_LET * NQUEUE_LET];
  __shared__  int  smem[NTHREADS_MAKE_LET];
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID
  const int target = occupyBuffer(tidx, freeLst, queue);
#else///USE_SMID_TO_GET_BUFID
#ifdef  TRY_MODE_ABOUT_BUFFER
  const int target = occupyBuffer(tidx, BLOCKIDX_X1D, freeNum, freeLst, queue);
#else///TRY_MODE_ABOUT_BUFFER
  occupyBuffer(tidx, freeNum, freeLst, queue, active);
#endif//TRY_MODE_ABOUT_BUFFER
#endif//USE_SMID_TO_GET_BUFID
  const int bufIdx = (int)queue[0];
  __syncthreads();
  size_t buf0Head = (size_t)bufIdx * (size_t)bufSize;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* sweep all tree nodes by executing tree-traversal */
  //-----------------------------------------------------------------------
  /* initialize queue for tree nodes */
#pragma unroll
  for(int jj = 0; jj < NQUEUE_LET; jj++)
    queue[tidx + NTHREADS_MAKE_LET * jj] = NULL_NODE;/* size >= NTHREADS_MAKE_LET * NQUEUE_LET */
  //-----------------------------------------------------------------------
  /* initialize queue for j-cells and interaction list by a representative thread */
  int sendNum = 0;/* # of LET nodes already stored in the global memory */
  int letTail = 0;/* the tail index of child cells for LET nodes already stored in the global memory */
  int bufHead = 0;
  int bufTail = 0;
  int bufOpen = bufSize;
  int bufUsed = 0;
  /* set child j-cells in queue on the shared memory */
  const int root = 0;
  uint jcell = more_org[root];
  int rem = 1 + (int)(jcell >> IDXBITS);
  jcell &= IDXMASK;
  //-----------------------------------------------------------------------
  if( rem > NTHREADS_MAKE_LET ){
    //---------------------------------------------------------------------
    /* if rem exceeds NTHREADS_MAKE_LET, then number of child j-cells must be shrunk */
    //---------------------------------------------------------------------
    queue[tidx] = jcell + (uint)tidx;/* size >= NTHREADS_MAKE_LET */
    //---------------------------------------------------------------------
    if( tidx == (NTHREADS_MAKE_LET - 1) )
      queue[tidx] += (uint)((rem - NTHREADS_MAKE_LET) << IDXBITS);/* size >= NTHREADS_MAKE_LET */
    //---------------------------------------------------------------------
    rem = NTHREADS_MAKE_LET;
    //---------------------------------------------------------------------
  }/* if( rem > NTHREADS_MAKE_LET ){ */
  else{
    //---------------------------------------------------------------------
    /* upload rem (<= NTHREADS_MAKE_LET) child j-cells to the shared memory */
    //---------------------------------------------------------------------
    if( tidx < rem )
      queue[tidx] = more_org[jcell + tidx];/* size >= NTHREADS_MAKE_LET */
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* tree traversal in a width-first manner */
  //-----------------------------------------------------------------------
  int fail = 0;
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    /* if the queue becomes empty, then exit the while loop */
    //---------------------------------------------------------------------
    __syncthreads();
    if( rem == 0 )
      break;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* pick up a queue from stack */
    //---------------------------------------------------------------------
    /* tentative load from the stack */
    int cnum = 0;
    jcell = NULL_NODE;
    if( tidx < rem ){
      jcell = queue[tidx];
      cnum = 1 + (int)(jcell >> IDXBITS);
    }/* if( lane < rem ){ */
    jcell &= IDXMASK;
    //---------------------------------------------------------------------
    /* predict the head index on the shared memory by parallel prefix sum */
    int hidx = prefixSum(cnum, tidx, lane, smem) - cnum;/* exclusive prefix sum of cnum */
    //---------------------------------------------------------------------
    smem[tidx] = NULL_NODE;
    __syncthreads();
    //---------------------------------------------------------------------
    int remove = 0;
    if( (cnum != 0) && (hidx < NTHREADS_MAKE_LET) ){
      //-------------------------------------------------------------------
      /* local data can be uploaded to the shared memory */
      int unum = NTHREADS_MAKE_LET - hidx;
      if( cnum < unum )	  unum = cnum;
      //-------------------------------------------------------------------
      /* upload local data */
      for(int jj = 0; jj < unum; jj++){
	/* list[hidx & (NTHREADS_MAKE_LET - 1)] = jcell; */
	smem[hidx] = (int)jcell;/* because hidx < NTHREADS_MAKE_LET */
	hidx++;
	jcell++;
      }/* for(int jj = 0; jj < unum; jj++){ */
      //-------------------------------------------------------------------
      /* eliminate stocked j-cells from the queue */
      if( unum == cnum )
	remove = 1;
      else{
	jcell += ((uint)(cnum - unum - 1) << IDXBITS);
	queue[tidx] = jcell;
      }/* else{ */
      //-------------------------------------------------------------------
    }/* if( (cnum != 0) && (hidx < NTHREADS_MAKE_LET) ){ */
    //---------------------------------------------------------------------
    /* set an index of j-cell */
    __syncthreads();
    const int target = smem[tidx];
    //---------------------------------------------------------------------
    /* remove scanned j-cells if possible */
    prefixSum(remove, tidx, lane, smem);
    remove = smem[NTHREADS_MAKE_LET - 1];
    //---------------------------------------------------------------------
    if( remove != 0 ){
      rem -= remove;
      copyData_s2s(queue, tidx + remove, queue, tidx, rem, tidx);
    }/* if( remove != 0 ){ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* pick up pseudo particles */
    //---------------------------------------------------------------------
    /* prefixSum to submit an LET node */
    int returnLET = (target != NULL_NODE) ? 1 : 0;
    hidx = sendNum + prefixSum(returnLET, tidx, lane, smem) - returnLET;/* index of the corresponding LET node, which is based on exclusive prefix sum of calc */
    sendNum += smem[NTHREADS_MAKE_LET - 1];
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* only the active threads pick up a j-cell from the global memory */
    //---------------------------------------------------------------------
    jparticle jpos_tmp;
    uint      more_tmp;
    int childNum = 0;
    int hasChild = 0;
    //---------------------------------------------------------------------
    if( returnLET ){
      //-------------------------------------------------------------------
      jpos_tmp     = jpos_org[target];      /* get position of pseudo j-particle */
      mj_let[hidx] =   mj_org[target];      /* send mj of an LET node */
      //-------------------------------------------------------------------
      /* set a pseudo i-particle */
      const real rx = jpos_tmp.x - icom.x;
      const real ry = jpos_tmp.y - icom.y;
      const real rz = jpos_tmp.z - icom.z;
      const real r2 = 1.0e-30f + rx * rx + ry * ry + rz * rz;
#if 1
      real lambda = FMAX(UNITY - SQRTRATIO(icom.m, r2), ZERO);
#else
      real lambda = UNITY - SQRTRATIO(icom.m, r2);
      if( lambda < EPSILON )	lambda = ZERO;
#endif
      /* calculate distance between the pseudo i-particle and the candidate j-particle */
      //-------------------------------------------------------------------
#ifdef  GADGET_MAC
      /* alpha * |a| * r^4 > G * M * l^2 */
#if 1
      lambda *= lambda * r2;
      if(   jpos_tmp.w < lambda               * lambda               * amin )
#else
	if( jpos_tmp.w < lambda * lambda * r2 * lambda * lambda * r2 * amin )
#endif
#else///GADGET_MAC
#ifdef  WS93_MAC
	  if(   jpos_tmp.w < lambda * lambda * r2 )
#else///WS93_MAC
	    /* (l / r) < theta */
	    if( jpos_tmp.w < lambda * lambda * r2 * theta2 )
#endif//WS93_MAC
#endif//GADGET_MAC
	      {
		//---------------------------------------------------------
		/* distant node ==>> child cells are not included in the LET */
		//---------------------------------------------------------
		more_tmp = hidx;
		jpos_tmp.w = -UNITY;/* squared size for the distant node is set to be negative */
		//---------------------------------------------------------
	      }
	    else
	      {
		//---------------------------------------------------------
		/* near node ==> child cells are included in the LET */
		//---------------------------------------------------------
		/* add child-cells of near tree-cells to the tentative stack */
		more_tmp = more_org[target];
		childNum = 1 + (int)(more_tmp >> IDXBITS);
		hasChild = 1;
		//---------------------------------------------------------
	      }
      //-------------------------------------------------------------------
    }/* if( returnLET ){ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* if the shared memory has open space and some tree cells are stored on the global memory, then load tree-cells from the global memory to the shared memory */
    //---------------------------------------------------------------------
    /* evaluate available size of the queue on the shared memory */
    int Nsm_rem = NQUEUE_LET * NTHREADS_MAKE_LET - rem;
    //---------------------------------------------------------------------
    if( (bufUsed != 0) && (Nsm_rem > 0) ){
      //-------------------------------------------------------------------
      /* hq is tidx */
      const int Nload = (Nsm_rem < bufUsed) ? (Nsm_rem) : (bufUsed);
#if 0
      memcpy((void *)&queue[rem], (const void *)&buffer[buf0Head + bufHead], sizeof(uint) * Nload);
#else
      copyData_g2s(buffer, buf0Head + bufHead, queue, rem, Nload, tidx);
#endif
      //-------------------------------------------------------------------
      rem     += Nload;
      Nsm_rem -= Nload;
      bufUsed -= Nload;
      bufHead += Nload;
      //-------------------------------------------------------------------
      if( bufUsed == 0 ){
	bufHead = 0;
	bufTail = 0;
	bufOpen = bufSize;
      }/* if( bufUsed == 0 ){ */
      //-------------------------------------------------------------------
    }/* if( (bufUsed != 0) && (Nsm_rem > 0) ){ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* copy child-cells of near tree-cells stored in the tentative stack to the stack on the shared memory and/or the global memory */
    //---------------------------------------------------------------------
    enqueueChildNodes(tidx, lane, smem, hasChild, more_tmp, queue, &Nsm_rem, &rem, buffer, buf0Head, &bufOpen, &bufUsed, &bufHead, &bufTail);
    fail += (bufOpen < 0);
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* if current node has child nodes in LET, then head index of more_tmp must be rewritten */
    //---------------------------------------------------------------------
    /* prefixSum to extend LET */
    int leafHead = prefixSum(childNum, tidx, lane, smem) - childNum;/* exclusive prefix sum of nchild */
    //---------------------------------------------------------------------
    /* modify more pointer using leafHead */
    if( childNum > 0 )
      more_tmp = ((uint)(childNum - 1) << IDXBITS) + (uint)(letTail + leafHead);
    letTail += smem[NTHREADS_MAKE_LET - 1];
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* add tree nodes to LET (mj_tmp is already stored) */
    //---------------------------------------------------------------------
    if( returnLET ){
      jpos_let[hidx] = jpos_tmp;
      more_let[hidx] = more_tmp;
    }/* if( returnLET ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* finalizing LET generator */
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    *numLETnode = sendNum;
    atomicAdd(overflow, fail);
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID
  releaseBuffer(tidx, freeLst, (uint)bufIdx, target);
#else///USE_SMID_TO_GET_BUFID
#ifdef  TRY_MODE_ABOUT_BUFFER
  releaseBuffer(tidx, freeLst, (uint)bufIdx, target);
#else///TRY_MODE_ABOUT_BUFFER
  releaseBuffer(tidx, freeNum, freeLst, bufIdx, active);
#endif//TRY_MODE_ABOUT_BUFFER
#endif//USE_SMID_TO_GET_BUFID
  //-----------------------------------------------------------------------
#ifdef  MONITOR_LETGEN_TIME
  long long int exitCycle = clock64();
  if( tidx == 0 ){
    unsigned long long int elapsed = (unsigned long long int)(exitCycle - initCycle);
    atomicAdd(cycles, elapsed);
  }/* if( tidx == 0 ){ */
#endif//MONITOR_LETGEN_TIME
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
 /* generate LET (Locally Essential Tree) */
//-------------------------------------------------------------------------
extern "C"
void callGenLET
  (const cudaStream_t stream, domainInfo *let, MPIcfg_tree mpi,
   const soaTreeNode tree, const int numSendGuess, const soaTreeWalkBuf buf
#ifdef  MONITOR_LETGEN_TIME
   , unsigned long long int * RESTRICT cycles
#endif//MONITOR_LETGEN_TIME
   )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  /* checkCudaErrors(cudaStreamSynchronize(stream)); */
  /* makeLET_kernel<<<1, NTHREADS_MAKE_LET>>> */
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
#ifdef  MONITOR_LETGEN_TIME
     , cycles
#endif//MONITOR_LETGEN_TIME
     );
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaMemcpy((*let).numSend_hst, (*let).numSend_dev, sizeof(int), cudaMemcpyDeviceToHost));
  /* checkCudaErrors(cudaMemcpyAsync((*let).numSend_hst, (*let).numSend_dev, sizeof(int), cudaMemcpyDeviceToHost, stream)); */
  /* checkCudaErrors(cudaStreamSynchronize(stream)); */
  let->numSend = *((*let).numSend_hst);
  if( let->numSend > numSendGuess ){
    __KILL__(stderr, "ERROR: predicted size of send buffer(%d) is not sufficient for true size of that(%d) @ rank %d for rand %d.\n\tsuggestion: consider increasing \"LETSIZE_OVERESTIMATION_FACTOR\" defined in src/tree/let.h (current value is %f).\n", numSendGuess, let->numSend, mpi.rank, let->rank, LETSIZE_OVERESTIMATION_FACTOR);
  }/* if( *numSend_hst > numSendGuess ){ */
  //-----------------------------------------------------------------------
#ifdef  DBG_LETGEN_ON_GPU
  fprintf(stdout, "numSend = %d @ rank %d\n", (*let).numSend, mpi.rank);
  /* checkCudaErrors(cudaDeviceSynchronize()); */
  /* MPI_Finalize(); */
  /* exit(0); */
#endif//DBG_LETGEN_ON_GPU
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#   if  !defined(GADGET_MAC) && !defined(WS93_MAC)
#   if  CUDART_VERSION >= 5000
  cudaMemcpyToSymbol( theta2 , &theta2_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#else//CUDART_VERSION >= 5000
  cudaMemcpyToSymbol("theta2", &theta2_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#endif//CUDART_VERSION >= 5000
#endif//!defined(GADGET_MAC) && !defined(WS93_MAC)
  //-----------------------------------------------------------------------
#   if  SMPREF_LET == 1
  checkCudaErrors(cudaFuncSetCacheConfig(makeLET_kernel, cudaFuncCachePreferShared));
#endif//SMPREF_LET == 1
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
