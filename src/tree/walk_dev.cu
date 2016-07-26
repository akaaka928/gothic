/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/26(Tue) 16:34:33
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
#include <helper_cuda.h>
#include <sys/time.h>
#ifndef SERIALIZED_EXECUTION
#       include <mpi.h>
#endif//SERIALIZED_EXECUTION
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
#       include <unistd.h>
#endif//PRINT_PSEUDO_PARTICLE_INFO
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
#include <timer.h>
#ifndef SERIALIZED_EXECUTION
#       include <mpilib.h>
#endif//SERIALIZED_EXECUTION
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
#       include <name.h>
#endif//PRINT_PSEUDO_PARTICLE_INFO
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/device.h"
//-------------------------------------------------------------------------
#include "macutil.h"
#include "make.h"
#include "buf_inc.h"
//-------------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
#       include "../para/mpicfg.h"
#       include "let.h"
#       include "let_dev.h"
#endif//SERIALIZED_EXECUTION
//-------------------------------------------------------------------------
#include "walk_dev.h"
//-------------------------------------------------------------------------
#define PARTIAL_SUM_ACCELERATION
//-------------------------------------------------------------------------
#ifdef  PARTIAL_SUM_ACCELERATION
#define ACCURATE_ACCUMULATION
/* #define ACCURATE_PARTIAL_SUM */
#endif//PARTIAL_SUM_ACCELERATION
//-------------------------------------------------------------------------
#   if  defined(ACCURATE_PARTIAL_SUM) && defined(ACCURATE_ACCUMULATION)
#undef          ACCURATE_PARTIAL_SUM
#endif//defined(ACCURATE_PARTIAL_SUM) && defined(ACCURATE_ACCUMULATION)
//-------------------------------------------------------------------------
/* #define MERGE_QUEUED_TREE_NODES */
//-------------------------------------------------------------------------
__constant__  real newton;
__constant__  real epsinv;
#ifndef INDIVIDUAL_GRAVITATIONAL_SOFTENING
__constant__  real eps2;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifndef WS93_MAC
__constant__  real theta2;
#endif//WS93_MAC
__constant__ jnode jnode0;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* set CUDA streams */
//-------------------------------------------------------------------------
extern "C"
muse setCUDAstreams_dev(cudaStream_t **stream, kernelStream *sinfo, deviceInfo *info, deviceProp *prop)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* determine # of CUDA streams */
  sinfo->idx = 0;
  sinfo->num = 2;
  //-----------------------------------------------------------------------
  /* allocate array for CUDA streams */
  *stream = (cudaStream_t *)malloc((size_t)(sinfo->num) * sizeof(cudaStream_t));  if( *stream == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stream");  }
  alloc.host +=                    (size_t)(sinfo->num) * sizeof(cudaStream_t) ;
  sinfo->stream = *stream;
  //-----------------------------------------------------------------------
  /* set CUDA streams */
  for(int ii = 0; ii < 2; ii++)
    sinfo->stream[ii] = info->stream[ii];
  for(int ii = 2; ii < sinfo->num; ii++)
    checkCudaErrors(cudaStreamCreate(&(sinfo->stream[ii])));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* arrays to store properties of tree cells (allocated on the global memory) */
//-------------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID
//-------------------------------------------------------------------------
/* complicated treatments is a remedy for ``not contiguous'' case of smid */
//-------------------------------------------------------------------------
__global__ void initFreeLst(const int numLanes, uint * RESTRICT freeLst, const int numFul, READ_ONLY int * RESTRICT smid)
{
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------
  if( tidx < numFul )
    freeLst[tidx] = INT_MAX;
  //-----------------------------------------------------------------------
  if( tidx < numLanes ){
    //---------------------------------------------------------------------
    const int target = (tidx % NBLOCKS_PER_SM) + smid[tidx / NBLOCKS_PER_SM] * NBLOCKS_PER_SM;
    //---------------------------------------------------------------------
    freeLst[target] = (uint)tidx;
    //---------------------------------------------------------------------
  }/* if( tidx < numLanes ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///USE_SMID_TO_GET_BUFID
//-------------------------------------------------------------------------
__global__ void initFreeLst
(const int numLanes, uint * RESTRICT freeLst
#ifndef TRY_MODE_ABOUT_BUFFER
 , uint * RESTRICT freeNum, int * RESTRICT active
#endif//TRY_MODE_ABOUT_BUFFER
 )
{
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------
  if( tidx < numLanes ){
    //---------------------------------------------------------------------
#ifdef  TRY_MODE_ABOUT_BUFFER
    freeLst[tidx] = (uint)tidx;
#else///TRY_MODE_ABOUT_BUFFER
    freeLst[tidx] = (uint)(numLanes - (tidx + 1));
#endif//TRY_MODE_ABOUT_BUFFER
    //---------------------------------------------------------------------
#ifndef TRY_MODE_ABOUT_BUFFER
    if( tidx == 0 ){
      *freeNum = (uint)numLanes;
      *active  = 1;
    }
#endif//TRY_MODE_ABOUT_BUFFER
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_SMID_TO_GET_BUFID
//-------------------------------------------------------------------------
extern "C"
void  freeTreeBuffer_dev
(int  *failure, uint  *buffer, uint  *freeLst
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
 , uint  *freeNum, int  *active
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
 , unsigned long long int  *cycles_hst, unsigned long long int  *cycles_dev
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
 , unsigned long long int  *cycles_let_hst, unsigned long long int  *cycles_let_dev
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  mycudaFree(failure);
  mycudaFree(buffer);
  mycudaFree(freeLst);
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
  mycudaFree(freeNum);
  mycudaFree(active);
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  mycudaFree    (cycles_dev);
  mycudaFreeHost(cycles_hst);
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  mycudaFree    (cycles_let_dev);
  mycudaFreeHost(cycles_let_hst);
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
muse allocTreeBuffer_dev
(int **failure, uint **buffer, uint **freeLst,
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
 uint **freeNum, int **active,
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
 unsigned long long int **cycles_hst, unsigned long long int **cycles_dev,
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
 unsigned long long int **cycles_let_hst, unsigned long long int **cycles_let_dev,
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
 soaTreeWalkBuf *buf, const int num_max, const muse used, const deviceProp gpu)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  mycudaMalloc((void **)failure, 1 * sizeof(int));
  alloc.device +=                1 * sizeof(int);
  const int fail_hst = 0;
  checkCudaErrors(cudaMemcpy(*failure, &fail_hst, sizeof(int), cudaMemcpyHostToDevice));
  //-----------------------------------------------------------------------
  const int nblocks = NBLOCKS_PER_SM * gpu.numSM;
  //-----------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID
  int last = 0;
  int num = 0;
  int *smid_dev;  mycudaMalloc    ((void **)&smid_dev, sizeof(int) * gpu.numSM);
  int *smid_hst;  mycudaMallocHost((void **)&smid_hst, sizeof(int) * gpu.numSM);
  for(int ii = 0; ii < 64; ii++)
    if( gpu.smid[ii] != -1 ){
      smid_hst[num] = gpu.smid[ii];      num++;
      last = ii;
    }
  last++;
  mycudaMalloc((void **)freeLst, (NBLOCKS_PER_SM * last) * sizeof(uint));  alloc.device += (NBLOCKS_PER_SM * last) * sizeof(uint);
#else///USE_SMID_TO_GET_BUFID
  mycudaMalloc((void **)freeLst, nblocks * sizeof(uint));  alloc.device += nblocks * sizeof(uint);
#endif//USE_SMID_TO_GET_BUFID
  //-----------------------------------------------------------------------
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
  mycudaMalloc((void **)freeNum,           sizeof(uint));  alloc.device +=           sizeof(uint);
  mycudaMalloc((void **) active,           sizeof( int));  alloc.device +=           sizeof( int);
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
  //-----------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID
  checkCudaErrors(cudaMemcpy(smid_dev, smid_hst, sizeof(int) * gpu.numSM, cudaMemcpyHostToDevice));
  initFreeLst<<<1, NBLOCKS_PER_SM * last>>>(nblocks, *freeLst, NBLOCKS_PER_SM * last, smid_dev);
  mycudaFree    (smid_dev);
  mycudaFreeHost(smid_hst);
#else///USE_SMID_TO_GET_BUFID
  initFreeLst<<<1, nblocks>>>(nblocks, *freeLst
#ifndef TRY_MODE_ABOUT_BUFFER
			      , *freeNum, *active
#endif//TRY_MODE_ABOUT_BUFFER
			      );
#endif//USE_SMID_TO_GET_BUFID
  //-----------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  mycudaMalloc    ((void **)cycles_dev, sizeof(unsigned long long int));  alloc.device += sizeof(unsigned long long int);
  mycudaMallocHost((void **)cycles_hst, sizeof(unsigned long long int));  alloc.host   += sizeof(unsigned long long int);
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  //-----------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  mycudaMalloc    ((void **)cycles_let_dev, sizeof(unsigned long long int));  alloc.device += sizeof(unsigned long long int);
  mycudaMallocHost((void **)cycles_let_hst, sizeof(unsigned long long int));  alloc.host   += sizeof(unsigned long long int);
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#if 1
  size_t unused, total;
  queryFreeDeviceMemory(&unused, &total);
#else
  const size_t unused = gpu.gmemSize - used.device - alloc.device;
#endif
#ifdef  CUB_AVAILABLE
  const size_t safety = GLOBAL_MEMORY_SYSBUF;
#else///CUB_AVAILABLE
  /* latters are pessimistic guess about device memory for CUDA thrust (PH-key sort, time step sort) */
  const size_t safety = GLOBAL_MEMORY_SYSBUF + (size_t)num_max * (sizeof(PHint) + sizeof(real));
#endif//CUB_AVAILABLE
  const size_t booked = (unused > safety) ? (unused - safety) : (unused >> 1);
  int bufUnit = (int)(booked / ((size_t)(NGROUPS * nblocks) * (sizeof(uint))));
  /* *bufUnit should be aligned in 32 bytes order (= 128 bits) --> 8 or 4 elements for single or double precision, respectively */
  bufUnit -= (bufUnit & 7);
  //-----------------------------------------------------------------------
  const size_t walkBufSize = (size_t)(NGROUPS * nblocks) * (size_t)bufUnit * sizeof(uint);
  mycudaMalloc((void **)buffer, walkBufSize);
  alloc.device +=               walkBufSize ;
  //-----------------------------------------------------------------------
  /* alert if the size for the walk buffer is smaller than 64 Ni B (512MiB @ N = 8M) */
  if( walkBufSize < ((size_t)num_max << 6) ){
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning:\tthe size for the walk buffer is %zu B (= %zu KiB = %zu MiB = %zu GiB), might be too small\n", walkBufSize, walkBufSize >> 10, walkBufSize >> 20, walkBufSize >> 30);
    fprintf(stderr, "suggestion:\tconsider decreasing \"TREE_SAFETY_VAL\" defined in src/tree/make.h (current value is %f)\n", TREE_SAFETY_VAL);
    fflush(stderr);
  }/* if( walkBufSize < ((size_t)num_max << 6) ){ */
  //-----------------------------------------------------------------------
  buf->fail    = *failure;
  buf->freeLst = *freeLst;
  buf->buffer  = *buffer;
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
  buf->freeNum = *freeNum;
  buf->active  = *active;
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  !defined(USE_SMID_TO_GET_BUFID) &&  defined(TRY_MODE_ABOUT_BUFFER)
  buf->freeNum = NBLOCKS_PER_SM * devProp.numSM;
#endif//!defined(USE_SMID_TO_GET_BUFID) &&  defined(TRY_MODE_ABOUT_BUFFER)
  buf->bufSize = bufUnit;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
//-------------------------------------------------------------------------
/* initialize count of Nj and Nbuf */
//-------------------------------------------------------------------------
__global__ void initCounter_kernel(int * RESTRICT Nj, int * RESTRICT Nb)
{
  //-----------------------------------------------------------------------
  Nj[GLOBALIDX_X1D] = 0;
  Nb[GLOBALIDX_X1D] = 0;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//COUNT_INTERACTIONS
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* initialize acceleration and potential */
//-------------------------------------------------------------------------
/* acc ::         output :: acceleration and potential of N-body particles */
//-------------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
//-------------------------------------------------------------------------
__global__ void initAcc_kernel
(acceleration * RESTRICT acc, const int laneNum, const laneinfo * RESTRICT laneInfo
#ifdef  GADGET_MAC
 , acceleration * RESTRICT old
#endif//GADGET_MAC
 )
{
  //-----------------------------------------------------------------------
#if 0
  const int tidx = THREADIDX_X1D;
  const int lane = tidx & (TSUB - 1);/* index of the thread within a thread group */
  /* const int laneIdx = GLOBALIDX_X1D / TSUB; */
  const int laneIdx = DIV_TSUB(GLOBALIDX_X1D);
#else
  /* const int lane    = THREADIDX_X1D & ((TSUB / NWARP) - 1); */
  /* const int laneIdx = GLOBALIDX_X1D /  (TSUB / NWARP); */
  const int lane    = THREADIDX_X1D & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);
#endif
#if 0
  const laneinfo info = laneInfo[laneIdx];
#else
  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];
#endif
  //-----------------------------------------------------------------------
  if( lane < info.num ){
    //---------------------------------------------------------------------
#ifdef  GADGET_MAC
    old[info.head + lane] = acc[info.head + lane];
#endif//GADGET_MAC
    //---------------------------------------------------------------------
    const acceleration ai = {ZERO, ZERO, ZERO, ZERO};
    acc[info.head + lane] = ai;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///BLOCK_TIME_STEP
//-------------------------------------------------------------------------
__global__ void initAcc_kernel
(acceleration *acc
#ifdef  GADGET_MAC
 , acceleration * RESTRICT old
#endif//GADGET_MAC
 )
{
  //-----------------------------------------------------------------------
  const acceleration ai = {ZERO, ZERO, ZERO, ZERO};
  //-----------------------------------------------------------------------
#ifdef  GADGET_MAC
  old[GLOBALIDX_X1D] = acc[GLOBALIDX_X1D];
#endif//GADGET_MAC
  //-----------------------------------------------------------------------
  acc[GLOBALIDX_X1D] = ai;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//BLOCK_TIME_STEP
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* multiply Gravitational constant and subtract self-interaction */
//-------------------------------------------------------------------------
/* acc :: input / output :: acceleration and potential of N-body particles */
/* pos :: input          :: position and mass of N-body particles */
//-------------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
//-------------------------------------------------------------------------
__global__ void trimAcc_kernel(acceleration * RESTRICT acc, READ_ONLY position * RESTRICT pos, const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo)
{
  //-----------------------------------------------------------------------
#if 0
  const int tidx = THREADIDX_X1D;
  const int lane = tidx & (TSUB - 1);/* index of the thread within a thread group */
  /* const int laneIdx = GLOBALIDX_X1D / TSUB; */
  const int laneIdx = DIV_TSUB(GLOBALIDX_X1D);
#else
  /* const int lane    = THREADIDX_X1D & ((TSUB / NWARP) - 1); */
  /* const int laneIdx = GLOBALIDX_X1D /  (TSUB / NWARP); */
  const int lane    = THREADIDX_X1D & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);
#endif
#if 0
  const laneinfo info = laneInfo[laneIdx];
#else
  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];
#endif
  //-----------------------------------------------------------------------
  if( lane < info.num ){
    //---------------------------------------------------------------------
    const int ii = info.head + lane;
    //---------------------------------------------------------------------
    /* load acceleration and mass */
    acceleration ai = acc[ii];
    const real   mi = pos[ii].m;
    //---------------------------------------------------------------------
    /* eliminate self-interaction */
    ai.pot -= epsinv * mi;
    //---------------------------------------------------------------------
    /* multiply Gravitational constant */
    ai.x   *=  newton;
    ai.y   *=  newton;
    ai.z   *=  newton;
    ai.pot *= -newton;
    //---------------------------------------------------------------------
    /* store acceleration */
    acc[ii] = ai;
    //---------------------------------------------------------------------
  }/* if( lane < info.num ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///BLOCK_TIME_STEP
//-------------------------------------------------------------------------
__global__ void trimAcc_kernel(acceleration * RESTRICT acc, READ_ONLY position * RESTRICT pos)
{
  //-----------------------------------------------------------------------
  const int ii = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  /* load acceleration and mass */
  acceleration ai = acc[ii];
  const real   mi = pos[ii].m;
  //-----------------------------------------------------------------------
  /* eliminate self-interaction */
  ai.pot -= epsinv * mi;
  //-----------------------------------------------------------------------
  /* multiply Gravitational constant */
  ai.x   *=  newton;
  ai.y   *=  newton;
  ai.z   *=  newton;
  ai.pot *= -newton;
  //-----------------------------------------------------------------------
  /* store acceleration */
  acc[ii] = ai;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//BLOCK_TIME_STEP
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* parallel prefix sum within a group of TSUB threads (TSUB <= 32 to use implicit synchronization) */
/* type of prefix sum is inclusive */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__  int prefixSumTsub(const int psum,                                            const int lane)
#else///USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ void prefixSumTsub(const int psum, volatile uint_real * smem, const int tidx, const int lane)
#endif//USE_WARP_SHUFFLE_FUNC
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  int val = psum;
  int tmp;
#   if  TSUB >=  2
  tmp = __shfl_up(val,  1, TSUB);  if( lane >=  1 )    val += tmp;
#   if  TSUB >=  4
  tmp = __shfl_up(val,  2, TSUB);  if( lane >=  2 )    val += tmp;
#   if  TSUB >=  8
  tmp = __shfl_up(val,  4, TSUB);  if( lane >=  4 )    val += tmp;
#   if  TSUB >= 16
  tmp = __shfl_up(val,  8, TSUB);  if( lane >=  8 )    val += tmp;
#   if  TSUB == 32
  tmp = __shfl_up(val, 16, TSUB);  if( lane >= 16 )    val += tmp;
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
  return (val);
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx].i = psum;
#   if  TSUB >=  2
  if( lane >=  1 )    smem[tidx].i += smem[tidx -  1].i;
#   if  TSUB >=  4
  if( lane >=  2 )    smem[tidx].i += smem[tidx -  2].i;
#   if  TSUB >=  8
  if( lane >=  4 )    smem[tidx].i += smem[tidx -  4].i;
#   if  TSUB >= 16
  if( lane >=  8 )    smem[tidx].i += smem[tidx -  8].i;
#   if  TSUB == 32
  if( lane >= 16 )    smem[tidx].i += smem[tidx - 16].i;
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
#ifdef  USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__  int prefixSumTsubMultiple(int psum,                                            const int lane, const int Niter)
#else///USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ void prefixSumTsubMultiple(int psum, volatile uint_real * smem, const int tidx, const int lane, const int Niter, const int tail)
#endif//USE_WARP_SHUFFLE_FUNC
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  int smem = prefixSumTsub(psum, lane);
#else///USE_WARP_SHUFFLE_FUNC
  prefixSumTsub(psum, smem, tidx, lane);
#endif//USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  for(int iter = 1; iter < Niter; iter++){
#ifdef  USE_WARP_SHUFFLE_FUNC
    const uint inc = (__shfl(smem, TSUB - 1, TSUB) >> (IDX_SHIFT_BITS * (iter - 1))) & IDX_SHIFT_MASK;
    smem         += (inc << (IDX_SHIFT_BITS * iter));
#else///USE_WARP_SHUFFLE_FUNC
    const uint inc = (smem[tail].i                 >> (IDX_SHIFT_BITS * (iter - 1))) & IDX_SHIFT_MASK;
    smem[tidx].i += (inc << (IDX_SHIFT_BITS * iter));
#endif//USE_WARP_SHUFFLE_FUNC
  }
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  return (smem);
#endif//USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* maximum value within a group of TSUB threads (TSUB <= 32 to use implicit synchronization) */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
/* continuous NWARP threads have the same value as input */
//-------------------------------------------------------------------------
__device__ __forceinline__ real getMaximumRealTsub
(
#ifdef  USE_WARP_SHUFFLE_FUNC
 const real max
#else///USE_WARP_SHUFFLE_FUNC
 real max, volatile uint_real * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
 )
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  real val = max;
#if 1
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = __shfl_xor(val,      NWARP, TSUB);  if( tmp > val )    val = tmp;
#   if  TSUB >= ( 4 * NWARP)
  tmp = __shfl_xor(val,  2 * NWARP, TSUB);  if( tmp > val )    val = tmp;
#   if  TSUB >= ( 8 * NWARP)
  tmp = __shfl_xor(val,  4 * NWARP, TSUB);  if( tmp > val )    val = tmp;
#   if  TSUB >= (16 * NWARP)
  tmp = __shfl_xor(val,  8 * NWARP, TSUB);  if( tmp > val )    val = tmp;
#   if  TSUB == (32 * NWARP)
  tmp = __shfl_xor(val, 16 * NWARP, TSUB);  if( tmp > val )    val = tmp;
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
#else
#   if  TSUB >=  2
  real tmp;
  tmp = __shfl_xor(val,  1, TSUB);  if( tmp > val )    val = tmp;
#   if  TSUB >=  4
  tmp = __shfl_xor(val,  2, TSUB);  if( tmp > val )    val = tmp;
#   if  TSUB >=  8
  tmp = __shfl_xor(val,  4, TSUB);  if( tmp > val )    val = tmp;
#   if  TSUB >= 16
  tmp = __shfl_xor(val,  8, TSUB);  if( tmp > val )    val = tmp;
#   if  TSUB == 32
  tmp = __shfl_xor(val, 16, TSUB);  if( tmp > val )    val = tmp;
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
#endif
  return (__shfl(val, 0, TSUB));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx].r = max;
  //-----------------------------------------------------------------------
#if 1
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = smem[tidx ^ (     NWARP)].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max;
#   if  TSUB >= ( 4 * NWARP)
  tmp = smem[tidx ^ ( 2 * NWARP)].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max;
#   if  TSUB >= ( 8 * NWARP)
  tmp = smem[tidx ^ ( 4 * NWARP)].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max;
#   if  TSUB >= (16 * NWARP)
  tmp = smem[tidx ^ ( 8 * NWARP)].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max;
#   if  TSUB == (32 * NWARP)
  tmp = smem[tidx ^ (16 * NWARP)].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max;
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
#else
#   if  TSUB >=  2
  real tmp;
  tmp = smem[tidx ^  1].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max;
#   if  TSUB >=  4
  tmp = smem[tidx ^  2].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max;
#   if  TSUB >=  8
  tmp = smem[tidx ^  4].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max;
#   if  TSUB >= 16
  tmp = smem[tidx ^  8].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max;
#   if  TSUB == 32
  tmp = smem[tidx ^ 16].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max;
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
#endif
  //-----------------------------------------------------------------------
  return (smem[head].r);
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* minimum value within a group of TSUB threads (TSUB <= 32 to use implicit synchronization) */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
/* continuous NWARP threads have the same value as input */
//-------------------------------------------------------------------------
__device__ __forceinline__ real getMinimumRealTsub
(
#ifdef  USE_WARP_SHUFFLE_FUNC
 const real min
#else///USE_WARP_SHUFFLE_FUNC
 real min, volatile uint_real * smem, const int tidx, const int head
#endif//USE_WARP_SHUFFLE_FUNC
 )
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  real val = min;
#if 1
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = __shfl_xor(val,      NWARP, TSUB);  if( tmp < val )    val = tmp;
#   if  TSUB >= ( 4 * NWARP)
  tmp = __shfl_xor(val,  2 * NWARP, TSUB);  if( tmp < val )    val = tmp;
#   if  TSUB >= ( 8 * NWARP)
  tmp = __shfl_xor(val,  4 * NWARP, TSUB);  if( tmp < val )    val = tmp;
#   if  TSUB >= (16 * NWARP)
  tmp = __shfl_xor(val,  8 * NWARP, TSUB);  if( tmp < val )    val = tmp;
#   if  TSUB == (32 * NWARP)
  tmp = __shfl_xor(val, 16 * NWARP, TSUB);  if( tmp < val )    val = tmp;
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
#else
#   if  TSUB >=  2
  real tmp;
  tmp = __shfl_xor(val,  1, TSUB);  if( tmp < val )    val = tmp;
#   if  TSUB >=  4
  tmp = __shfl_xor(val,  2, TSUB);  if( tmp < val )    val = tmp;
#   if  TSUB >=  8
  tmp = __shfl_xor(val,  4, TSUB);  if( tmp < val )    val = tmp;
#   if  TSUB >= 16
  tmp = __shfl_xor(val,  8, TSUB);  if( tmp < val )    val = tmp;
#   if  TSUB == 32
  tmp = __shfl_xor(val, 16, TSUB);  if( tmp < val )    val = tmp;
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
#endif
  return (__shfl(val, 0, TSUB));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx].r = min;
  //-----------------------------------------------------------------------
#if 1
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = smem[tidx ^ (     NWARP)].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min;
#   if  TSUB >= ( 4 * NWARP)
  tmp = smem[tidx ^ ( 2 * NWARP)].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min;
#   if  TSUB >= ( 8 * NWARP)
  tmp = smem[tidx ^ ( 4 * NWARP)].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min;
#   if  TSUB >= (16 * NWARP)
  tmp = smem[tidx ^ ( 8 * NWARP)].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min;
#   if  TSUB == (32 * NWARP)
  tmp = smem[tidx ^ (16 * NWARP)].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min;
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
#else
#   if  TSUB >=  2
  real tmp;
  tmp = smem[tidx ^  1].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min;
#   if  TSUB >=  4
  tmp = smem[tidx ^  2].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min;
#   if  TSUB >=  8
  tmp = smem[tidx ^  4].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min;
#   if  TSUB >= 16
  tmp = smem[tidx ^  8].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min;
#   if  TSUB == 32
  tmp = smem[tidx ^ 16].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min;
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
#endif
  //-----------------------------------------------------------------------
  return (smem[head].r);
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_s2s
(uint *src, int sidx,
 uint *dst, int didx, const int num, const int lane)
{
  //-----------------------------------------------------------------------
  const int iter = DIV_TSUB(num);
  const int frac = num & (TSUB - 1);/* := Nload % TSUB */
  //-----------------------------------------------------------------------
  /* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
  for(int kk = 0; kk < iter; kk++){
    dst[didx] = src[sidx];
    sidx += TSUB;
    didx += TSUB;
  }
  //-----------------------------------------------------------------------
  if( lane < frac )
    dst[didx] = src[sidx];
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_g2s
(uint * RESTRICT gbuf, size_t srcHead, uint * RESTRICT sbuf, int dstHead, int numCopy, const int lane)
{
  //-----------------------------------------------------------------------
  /* fraction processing at loading from the head of destination array */
  //-----------------------------------------------------------------------
  const int numTemp = TSUB - (int)(srcHead & (TSUB - 1));/* := TSUB - (srcHead % TSUB) */
  const int numHead = (numTemp < numCopy) ? numTemp : numCopy;
  if( lane < numHead )
    sbuf[dstHead + lane] = gbuf[srcHead + (size_t)lane];
  dstHead +=         numHead;
  srcHead += (size_t)numHead;
  numCopy -=         numHead;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* sequential load from source on the global memory and store to destination on the shared memory */
  //-----------------------------------------------------------------------
  for(int ii = lane; ii < numCopy; ii += TSUB)
    sbuf[dstHead + ii] = gbuf[srcHead + (size_t)ii];
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_s2g
(uint * RESTRICT sbuf, int srcHead, uint * RESTRICT gbuf, size_t dstHead, int numCopy, const int lane)
{
  //-----------------------------------------------------------------------
  /* fraction processing at storing to the head of destination array */
  //-----------------------------------------------------------------------
  const int numTemp = TSUB - (int)(dstHead & (TSUB - 1));/* := TSUB - (dstHead % TSUB) */
  const int numHead = (numTemp < numCopy) ? numTemp : numCopy;
  if( lane < numHead )
    gbuf[dstHead + (size_t)lane] = sbuf[srcHead + lane];
  dstHead += (size_t)numHead;
  srcHead +=         numHead;
  numCopy -=         numHead;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* sequential load from source on the shared memory and store to destination on the global memory */
  //-----------------------------------------------------------------------
  for(int ii = lane; ii < numCopy; ii += TSUB)
    gbuf[dstHead + (size_t)ii] = sbuf[srcHead + ii];
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_g2g
(uint * RESTRICT gbuf, int srcHead, int dstHead, int Ncopy, const int Ndisp, const int lane)
{
  //-----------------------------------------------------------------------
  /* configure the settings */
  //-----------------------------------------------------------------------
  const int Nfirst = Ndisp & (TSUB - 1);/* := Ndisp % TSUB */
  /* ldIdx is Nfirst, Nfirst + 1, ..., TSUB - 1, 0, 1, ..., Nfirst - 1 for lane of 0, 1, 2, ..., TSUB - 1 */
  const int  ldIdx = (lane + Nfirst) & (TSUB - 1);/* := (lane + Nfirst) % TSUB */
  const int grpIdx = (ldIdx < Nfirst) ? 0 : 1;
  //-----------------------------------------------------------------------
  srcHead += Ndisp - Nfirst;/* hereafter, srcHead is TSUB elements aligned */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* fraction processing at loading from the head of source array */
  //-----------------------------------------------------------------------
  uint temp = gbuf[srcHead + ldIdx];
  srcHead += TSUB;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#if 0
  //-----------------------------------------------------------------------
  /* sequential load and store from source to destination on the global memory */
  /* this implementation employes the local memory */
  //-----------------------------------------------------------------------
  const int Niter = BLOCKSIZE(Ncopy, TSUB * NSTOCK);
  for(int iter = 0; iter < Niter; iter++){
    //---------------------------------------------------------------------
    const int Nmove = (Ncopy > (TSUB * NSTOCK)) ? (TSUB * NSTOCK) : (Ncopy);
    const int Nloop = BLOCKSIZE(Nmove, TSUB);
    //---------------------------------------------------------------------
    //
    //---------------------------------------------------------------------
    /* load from the source array on the global memory */
    //---------------------------------------------------------------------
    /* load from temp (fraction processing) as initialization */
    union {uint4 i; uint a[4];} local = {temp, NULL_NODE, NULL_NODE, NULL_NODE};
    //---------------------------------------------------------------------
    /* load from global memory */
#pragma unroll
    for(int ii = 0; ii < Nloop - 1; ii++)
      local.a[grpIdx + ii] = gbuf[srcHead + TSUB * ii + ldIdx];
    //---------------------------------------------------------------------
    /* load from global memory, store to shared memory or temp (fraction processing) */
    temp = gbuf[srcHead + TSUB * (Nloop - 1) + ldIdx];
    if( !grpIdx ){
      local.a[Nloop - 1] = temp;
    }/* if( !grpIdx ){ */
    //---------------------------------------------------------------------
    //
    //---------------------------------------------------------------------
    /* store to the destination array on the global memory */
    //---------------------------------------------------------------------
#pragma unroll
    for(int ii = 0; ii < Nloop; ii++)
      gbuf[dstHead + TSUB * ii + lane] = local.a[ii];
    //---------------------------------------------------------------------
    Ncopy   -= Nmove;
    srcHead += Nmove;
    dstHead += Nmove;
    //---------------------------------------------------------------------
  }/* for(int iter = 0; iter < Niter; iter++){ */
  //-----------------------------------------------------------------------
#else
  //-----------------------------------------------------------------------
  /* sequential load and store from source to destination on the global memory */
  //-----------------------------------------------------------------------
  const int Niter = BLOCKSIZE(Ncopy, TSUB);
  for(int iter = 0; iter < Niter; iter++){
    //---------------------------------------------------------------------
    const int Nmove = (Ncopy > TSUB) ? (TSUB) : (Ncopy);
    //---------------------------------------------------------------------
    //
    //---------------------------------------------------------------------
    /* load from the source array on the global memory */
    //---------------------------------------------------------------------
    /* load from temp (fraction processing) as initialization */
    uint local = temp;
    //---------------------------------------------------------------------
    /* load from global memory, store to shared memory or temp (fraction processing) */
    temp = gbuf[srcHead + ldIdx];
    if( !grpIdx )
      local = temp;
    //---------------------------------------------------------------------
    //
    //---------------------------------------------------------------------
    /* store to the destination array on the global memory */
    //---------------------------------------------------------------------
    gbuf[dstHead + lane] = local;
    //---------------------------------------------------------------------
    Ncopy   -= Nmove;
    srcHead += Nmove;
    dstHead += Nmove;
    //---------------------------------------------------------------------
  }/* for(int iter = 0; iter < Niter; iter++){ */
  //-----------------------------------------------------------------------
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* merge continuous tree nodes */
//-------------------------------------------------------------------------
/* uint smem[TSUB]; */
/* uint node[TSUB * NSTOCK]; */
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ void cpChildNodes
(                           uint * RESTRICT node, jnode jidx,
 uint leaf,                 const int lane,
 uint * RESTRICT smbuf, const    int hq, int *rem_sm, int *num_sm,
 uint * RESTRICT gmbuf, const size_t hb, int *rem_gm, int *num_gm, int *head_gm, int *tail_gm
			    )
#else///USE_WARP_SHUFFLE_FUNC
__device__ __forceinline__ void cpChildNodes
(uint_real * RESTRICT smem, uint * RESTRICT node, jnode jidx,
 uint leaf, const int tidx, const int lane, const int tail,
 uint * RESTRICT smbuf, const    int hq, int *rem_sm, int *num_sm,
 uint * RESTRICT gmbuf, const size_t hb, int *rem_gm, int *num_gm, int *head_gm, int *tail_gm
 )
#endif//USE_WARP_SHUFFLE_FUNC
{
  //-----------------------------------------------------------------------
  int iter;
  //-----------------------------------------------------------------------
  /* 1. compact the given sparse tree nodes */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  int smem = prefixSumTsubMultiple(leaf, lane, NSTOCK);
  uint nadd = smem         - leaf;/* exclusive prefix sum of leaf */
#else///USE_WARP_SHUFFLE_FUNC
  prefixSumTsubMultiple(leaf, smem, tidx, lane, NSTOCK, tail);
  uint nadd = smem[tidx].i - leaf;/* exclusive prefix sum of leaf */
#endif//USE_WARP_SHUFFLE_FUNC
#pragma unroll
  for(iter = 0; iter < NSTOCK; iter++){
    if( (leaf >> (IDX_SHIFT_BITS * iter)) & 1 ){
      const uint hidx = (nadd >> (IDX_SHIFT_BITS * iter)) & IDX_SHIFT_MASK;
      node[hidx] = jidx.idx[iter];
    }/* if( (leaf >> (IDX_SHIFT_BITS * iter)) & 1 ){ */
    jidx.idx[iter] = NULL_NODE;
  }/* for(iter = 0; iter < NSTOCK; iter++){ */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
  const int nold = (int)((__shfl(smem, TSUB - 1, TSUB) >> (IDX_SHIFT_BITS * (NSTOCK - 1))) & IDX_SHIFT_MASK);
#else///USE_WARP_SHUFFLE_FUNC
  const int nold = (int)((       smem[tail].i          >> (IDX_SHIFT_BITS * (NSTOCK - 1))) & IDX_SHIFT_MASK);
#endif//USE_WARP_SHUFFLE_FUNC
  for(int ii = nold + lane; ii < NSTOCK * TSUB; ii += TSUB)
    node[ii] = NULL_NODE;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  int Ntot;
  //-----------------------------------------------------------------------
#ifdef  MERGE_QUEUED_TREE_NODES
  const int Niter = BLOCKSIZE(nold, TSUB);
  if( Niter != NSTOCK )
#endif//MERGE_QUEUED_TREE_NODES
#ifdef  USE_WARP_SHUFFLE_FUNC
    Ntot = (int)((__shfl(smem, TSUB - 1, TSUB) >> (IDX_SHIFT_BITS * (NSTOCK - 1))) & IDX_SHIFT_MASK);
#else///USE_WARP_SHUFFLE_FUNC
    Ntot = (int)((       smem[tail].i          >> (IDX_SHIFT_BITS * (NSTOCK - 1))) & IDX_SHIFT_MASK);
#endif//USE_WARP_SHUFFLE_FUNC
#ifdef  MERGE_QUEUED_TREE_NODES
  else{
    //---------------------------------------------------------------------
#if 1
    /* partial, faster version */
    //---------------------------------------------------------------------
    /* 2. examine continuity of the given tree nodes */
    /* 3. construct merged tree nodes */
    //---------------------------------------------------------------------
    iter = 0;
    leaf = 0;
#pragma unroll
    for(int ii = 2 * lane; ii < nold; ii += 2 * TSUB){
      //-------------------------------------------------------------------
      jidx.idx[2 * iter    ] = node[ii    ];
      jidx.idx[2 * iter + 1] = node[ii + 1];
      //-------------------------------------------------------------------
      const uint  numFormer = (jidx.idx[2 * iter    ] >> IDXBITS) + 1;
      const uint  numLatter = (jidx.idx[2 * iter + 1] >> IDXBITS) + 1;
      const uint tailFormer = (jidx.idx[2 * iter    ] &  IDXMASK) + numFormer;/* := tail index + 1 */
      const uint headLatter =  jidx.idx[2 * iter + 1] &  IDXMASK;
      //-------------------------------------------------------------------
      if( (tailFormer == headLatter) && ((numFormer + numLatter) <= NLEAF) ){
	jidx.idx[2 * iter    ] += (numLatter << IDXBITS);
	jidx.idx[2 * iter + 1]  = NULL_NODE;
      }/* if( (tailFormer == headLatter) && ((numFormer + numLatter) <= NLEAF) ){ */
      //-------------------------------------------------------------------
      uint numNodes = 0;
#pragma unroll
      for(int jj = 0; jj < 2; jj++)
	numNodes += (jidx.idx[2 * iter + jj] != NULL_NODE);
      leaf += (numNodes << (IDX_SHIFT_BITS * iter));
      //-------------------------------------------------------------------
      iter++;
      //-------------------------------------------------------------------
    }/* for(int ii = 2 * lane; ii < nold; ii += 2 * TSUB){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* 4. count up number of reconstructed tree nodes */
    //---------------------------------------------------------------------
    const int Nloop = BLOCKSIZE(nold, 2 * TSUB);
#ifdef  USE_WARP_SHUFFLE_FUNC
    smem = prefixSumTsubMultiple(leaf, lane, Nloop);
#else///USE_WARP_SHUFFLE_FUNC
    prefixSumTsubMultiple(leaf, smem, tidx, lane, Nloop, tail);
#endif//USE_WARP_SHUFFLE_FUNC
    //---------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
    Ntot = (int)((__shfl(smem, TSUB - 1, TSUB) >> (IDX_SHIFT_BITS * (Nloop - 1))) & IDX_SHIFT_MASK);
#else///USE_WARP_SHUFFLE_FUNC
    Ntot = (int)((       smem[tail].i          >> (IDX_SHIFT_BITS * (Nloop - 1))) & IDX_SHIFT_MASK);
#endif//USE_WARP_SHUFFLE_FUNC
    //---------------------------------------------------------------------
    for(int ii = Ntot + lane; ii < nold; ii += TSUB)
      node[ii] = NULL_NODE;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* 5. set the reconstructed tree nodes on the shared memory */
    //---------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
    smem         -= leaf;/* exclusive prefix sum */
#else///USE_WARP_SHUFFLE_FUNC
    smem[tidx].i -= leaf;/* exclusive prefix sum */
#endif//USE_WARP_SHUFFLE_FUNC
#pragma unroll
    for(int ii = 0; ii < Nloop; ii++){
      //-------------------------------------------------------------------
      const int  numNodes = (int)((leaf         >> (IDX_SHIFT_BITS * ii)) & IDX_SHIFT_MASK);
#ifdef  USE_WARP_SHUFFLE_FUNC
      const int headNodes = (int)((smem         >> (IDX_SHIFT_BITS * ii)) & IDX_SHIFT_MASK);
#else///USE_WARP_SHUFFLE_FUNC
      const int headNodes = (int)((smem[tidx].i >> (IDX_SHIFT_BITS * ii)) & IDX_SHIFT_MASK);
#endif//USE_WARP_SHUFFLE_FUNC
      //-------------------------------------------------------------------
#pragma unroll
      for(int jj = 0; jj < numNodes; jj++)
	node[headNodes + jj] = jidx.idx[2 * ii + jj];
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < Nloop; ii++){ */
    //---------------------------------------------------------------------
#else
    /* complete, slower version */
    //---------------------------------------------------------------------
    /* 2. examine continuity of the given tree nodes */
    //---------------------------------------------------------------------
    iter = 0;
    leaf = 0;
#pragma unroll
    for(int ii = lane; ii < nold; ii += TSUB){
      //-------------------------------------------------------------------
      jidx.idx[iter] = node[ii];
      //-------------------------------------------------------------------
      if( ii != 0 ){
	//-----------------------------------------------------------------
	uint tail_id = node[ii - 1];
	const uint num  = 1 + (tail_id >> IDXBITS);
	tail_id = (tail_id & IDXMASK) + num;/* := tail index + 1 */
	//-----------------------------------------------------------------
	leaf += ((tail_id == (jidx.idx[iter] & IDXMASK)) << (IDX_SHIFT_BITS * iter));
	//-----------------------------------------------------------------
      }/* if( ii != 0 ){ */
      //-------------------------------------------------------------------
      iter++;
      //-------------------------------------------------------------------
    }/* for(int ii = lane; ii < nold; ii += TSUB){ */
    //---------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
    smem = prefixSumTsubMultiple(leaf, lane, Niter);
#else///USE_WARP_SHUFFLE_FUNC
    prefixSumTsubMultiple(leaf, smem, tidx, lane, Niter, tail);
#endif//USE_WARP_SHUFFLE_FUNC
    //---------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
    const int nnew = nold - (int)((__shfl(smem, TSUB - 1, TSUB) >> (IDX_SHIFT_BITS * (Niter - 1))) & IDX_SHIFT_MASK);
#else///USE_WARP_SHUFFLE_FUNC
    const int nnew = nold - (int)((       smem[tail].i          >> (IDX_SHIFT_BITS * (Niter - 1))) & IDX_SHIFT_MASK);
#endif//USE_WARP_SHUFFLE_FUNC
    const int Nloop = BLOCKSIZE(nnew, TSUB);
    //---------------------------------------------------------------------
    for(int ii = nnew + lane; ii < nold; ii += TSUB)
      node[ii] = NULL_NODE;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* 3. construct merged tree nodes */
    //---------------------------------------------------------------------
    /* store head index to the shared memory */
    iter = 0;
#pragma unroll
    for(int ii = lane; ii < nold; ii += TSUB){
      if( !((leaf >> (IDX_SHIFT_BITS * iter)) & 1) )
#ifdef  USE_WARP_SHUFFLE_FUNC
	node[ii - (int)((smem         >> (IDX_SHIFT_BITS * iter)) & IDX_SHIFT_MASK)] = jidx.idx[iter] & IDXMASK;
#else///USE_WARP_SHUFFLE_FUNC
	node[ii - (int)((smem[tidx].i >> (IDX_SHIFT_BITS * iter)) & IDX_SHIFT_MASK)] = jidx.idx[iter] & IDXMASK;
#endif//USE_WARP_SHUFFLE_FUNC
      iter++;
    }/* for(int ii = lane; ii < nold; ii += TSUB){ */
    //---------------------------------------------------------------------
    /* load head index from the shared memory */
    jnode head;
    iter = 0;
#pragma unroll
    for(int ii = lane; ii < Nloop * TSUB; ii += TSUB){
      head.idx[iter] = node[ii];
      node[ii] = 0;
      iter++;
    }/* for(int ii = lane; ii < Nloop * TSUB; ii += TSUB){ */
    //---------------------------------------------------------------------
    /* sum up number of child nodes */
    /* TENTATIVE IMPLEMENTATION: in future update, atomic operation should be removed */
    iter = 0;
#pragma unroll
    for(int ii = lane; ii < nold; ii += TSUB){
      if( jidx.idx[iter] != NULL_NODE )
#ifdef  USE_WARP_SHUFFLE_FUNC
	atomicAdd(&node[ii - (int)((smem         >> (IDX_SHIFT_BITS * iter)) & IDX_SHIFT_MASK)], 1 + (jidx.idx[iter] >> IDXBITS));
#else///USE_WARP_SHUFFLE_FUNC
	atomicAdd(&node[ii - (int)((smem[tidx].i >> (IDX_SHIFT_BITS * iter)) & IDX_SHIFT_MASK)], 1 + (jidx.idx[iter] >> IDXBITS));
#endif//USE_WARP_SHUFFLE_FUNC
      iter++;
    }/* for(int ii = lane; ii < nold; ii += TSUB){ */
    //---------------------------------------------------------------------
    /* load number of child nodes from the shared memory */
    jnode num;
    iter = 0;
#pragma unroll
    for(int ii = lane; ii < Nloop * TSUB; ii += TSUB){
      num.idx[iter] = node[ii];
      iter++;
    }/* for(int ii = lane; ii < Nloop * TSUB; ii += TSUB){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* 4. count up number of reconstructed tree nodes */
    //---------------------------------------------------------------------
    nadd = 0;
#pragma unroll
    for(iter = 0; iter < Nloop; iter++)
      nadd += ((uint)(BLOCKSIZE(num.idx[iter], NLEAF)) << (IDX_SHIFT_BITS * iter));
#ifdef  USE_WARP_SHUFFLE_FUNC
    smem = prefixSumTsubMultiple(nadd, lane, Nloop);
#else///USE_WARP_SHUFFLE_FUNC
    prefixSumTsubMultiple(nadd, smem, tidx, lane, Nloop, tail);
#endif//USE_WARP_SHUFFLE_FUNC
    //---------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
    Ntot = (int)((__shfl(smem, TSUB - 1, TSUB) >> (IDX_SHIFT_BITS * (Nloop - 1))) & IDX_SHIFT_MASK);
#else///USE_WARP_SHUFFLE_FUNC
    Ntot = (int)((       smem[tail].i          >> (IDX_SHIFT_BITS * (Nloop - 1))) & IDX_SHIFT_MASK);
#endif//USE_WARP_SHUFFLE_FUNC
    //---------------------------------------------------------------------
#pragma unroll
    for(int ii = Ntot + lane; ii < Nloop * TSUB; ii += TSUB)
      node[ii] = NULL_NODE;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* 5. set the reconstructed tree nodes on the shared memory */
    //---------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
    smem         -= nadd;/* exclusive prefix sum */
#else///USE_WARP_SHUFFLE_FUNC
    smem[tidx].i -= nadd;/* exclusive prefix sum */
#endif//USE_WARP_SHUFFLE_FUNC
    iter = 0;
    for(int ii = lane; ii < nnew; ii += TSUB){
      //-------------------------------------------------------------------
      const int nsplit = (int)((nadd         >> (IDX_SHIFT_BITS * iter)) & IDX_SHIFT_MASK);
#ifdef  USE_WARP_SHUFFLE_FUNC
      const int hidx   = (int)((smem         >> (IDX_SHIFT_BITS * iter)) & IDX_SHIFT_MASK);
#else///USE_WARP_SHUFFLE_FUNC
      const int hidx   = (int)((smem[tidx].i >> (IDX_SHIFT_BITS * iter)) & IDX_SHIFT_MASK);
#endif//USE_WARP_SHUFFLE_FUNC
      //-------------------------------------------------------------------
      uint idx_node = head.idx[iter];
      uint rem_node =  num.idx[iter];
      //-------------------------------------------------------------------
      for(int jj = 0; jj < nsplit; jj++){
	//-----------------------------------------------------------------
	const uint num_node = (rem_node < NLEAF) ? rem_node : NLEAF;
	node[hidx + jj] = ((num_node - 1) << IDXBITS) + idx_node;
	//-----------------------------------------------------------------
	idx_node += num_node;
	rem_node -= num_node;
	//-----------------------------------------------------------------
      }/* for(int jj = 0; jj < nsplit; jj++){ */
      //-------------------------------------------------------------------
      iter++;
      //-------------------------------------------------------------------
    }/* for(int ii = lane; ii < nnew; ii += TSUB){ */
    //---------------------------------------------------------------------
#endif
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
#endif//MERGE_QUEUED_TREE_NODES
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* 6. copy merged tree nodes to the shared memory */
  //-----------------------------------------------------------------------
  const int Nsm = (Ntot < *rem_sm) ? (Ntot) : (*rem_sm);
  copyData_s2s(node, lane, smbuf, hq + (*num_sm), Nsm, lane);
  //-----------------------------------------------------------------------
  *num_sm += Nsm;
  *rem_sm -= Nsm;
  Ntot    -= Nsm;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 7. move tree nodes on the global memory, if necessary */
  //-----------------------------------------------------------------------
  if( Ntot > *rem_gm ){
    //---------------------------------------------------------------------
#if 1
    copyData_g2g(gmbuf, hb, hb, *num_gm, *head_gm, lane);
#else
    /* probably, below implementation is wrong */
#pragma unroll
    for(int ii = lane; ii < *num_gm; ii += TSUB)
      gmbuf[hb + (size_t)((*head_gm) + ii)] = gmbuf[hb + (size_t)ii];
#endif
    //---------------------------------------------------------------------
    * rem_gm += *head_gm;
    *tail_gm -= *head_gm;
    *head_gm  = 0;
    //---------------------------------------------------------------------
  }/* if( Ntot > *rem_gm ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 8. copy merged tree nodes to the global memory */
  //-----------------------------------------------------------------------
#if 0
  memcpy((void *)&gmbuf[hb + (size_t)(*tail_gm)], (const void *)&node[Nsm], sizeof(uint) * Ntot);
#else
  copyData_s2g(node, Nsm, gmbuf, hb + (size_t)(*tail_gm), Ntot, lane);
#endif
  //-----------------------------------------------------------------------
  * rem_gm -= Ntot;
  * num_gm += Ntot;
  *tail_gm += Ntot;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* calculate body--body interaction based on direct summation */
//-------------------------------------------------------------------------
/* pi   :: input          :: position and eps2 of N-body particles */
/* ai   :: input / output :: acceleration and potential of N-body particles */
/* jpos :: input          :: position and mass of N-body particles */
//-------------------------------------------------------------------------
__device__ __forceinline__ void calc_interaction
(const position pi, acceleration * RESTRICT ai, jparticle * RESTRICT jpos
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
 , real * RESTRICT eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  ACCURATE_ACCUMULATION
 , acceleration * RESTRICT res
#endif//ACCURATE_ACCUMULATION
#ifdef  IJ_PARALLELIZATION
 , const int lane
#endif//IJ_PARALLELIZATION
)
{
  //-----------------------------------------------------------------------
#ifdef  PARTIAL_SUM_ACCELERATION
#       ifdef  ACCURATE_PARTIAL_SUM
  acceleration res_loc = {ZERO, ZERO, ZERO, ZERO};
#       endif//ACCURATE_PARTIAL_SUM
  acceleration acc     = {ZERO, ZERO, ZERO, ZERO};
#endif//PARTIAL_SUM_ACCELERATION
  //-----------------------------------------------------------------------
#if 0
#       if DIV_NWARP(NLOOP * TSUB) == 128
#pragma unroll 64
#       endif
#       if DIV_NWARP(NLOOP * TSUB) == 96
#pragma unroll 48
#       endif
#       if DIV_NWARP(NLOOP * TSUB) == 64
#pragma unroll 32
#       endif
#       if DIV_NWARP(NLOOP * TSUB) == 48
#pragma unroll 24
#       endif
#       if DIV_NWARP(NLOOP * TSUB) == 32
#pragma unroll 16
#       endif
#       if DIV_NWARP(NLOOP * TSUB) == 24
#pragma unroll 12
#       endif
#       if DIV_NWARP(NLOOP * TSUB) == 16
#pragma unroll  8
#       endif
#       if DIV_NWARP(NLOOP * TSUB) == 12
#pragma unroll  6
#       endif
#       if DIV_NWARP(NLOOP * TSUB) ==  8
#pragma unroll  4
#       endif
#       if DIV_NWARP(NLOOP * TSUB) ==  6
#pragma unroll  3
#       endif
#       if DIV_NWARP(NLOOP * TSUB) ==  4
#pragma unroll  2
#       endif
#else
#pragma unroll
#endif
#ifdef  IJ_PARALLELIZATION
  for(int jj = lane; jj < NLOOP * TSUB; jj += NWARP)
#else///IJ_PARALLELIZATION
  for(int jj = 0; jj < NLOOP * TSUB; jj++)
#endif//IJ_PARALLELIZATION
    {
      //-------------------------------------------------------------------
      /* load j-particle from shared memory */
      jparticle pj = jpos[jj];
      //-------------------------------------------------------------------
      /* calculate distance between j-particel and i-particle */
      const real rx = pj.x - pi.x;
      const real ry = pj.y - pi.y;
      const real rz = pj.z - pi.z;
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
      const real r2 = eps2[jj] + rx * rx + ry * ry + rz * rz;
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
      const real r2 = pi.m     + rx * rx + ry * ry + rz * rz;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
      real rinv = RSQRT(r2);
      //-------------------------------------------------------------------
      /* calculate common factor for all direction */
      pj.w *= rinv;/* mj / r */
      rinv *= rinv;/* 1  / r^2 */
      rinv *= pj.w;/* mj / r^3 */
      //-------------------------------------------------------------------
      /* calculate gravitational acceleration of i-particle */
#ifdef  PARTIAL_SUM_ACCELERATION
#       ifdef  ACCURATE_PARTIAL_SUM
      /* R := R + x_i */
      res_loc.x   += rx * rinv;
      res_loc.y   += ry * rinv;
      res_loc.z   += rz * rinv;
      res_loc.pot += r2 * rinv;
      /* T := S */
      acceleration tmp_loc = acc;
      /* S := S + R */
      acc.x   += res_loc.x;
      acc.y   += res_loc.y;
      acc.z   += res_loc.z;
      acc.pot += res_loc.pot;
      /* T := S - T */
      tmp_loc.x   = acc.x   - tmp_loc.x;
      tmp_loc.y   = acc.y   - tmp_loc.y;
      tmp_loc.z   = acc.z   - tmp_loc.z;
      tmp_loc.pot = acc.pot - tmp_loc.pot;
      /* R := R - T */
      res_loc.x   -= tmp_loc.x;
      res_loc.y   -= tmp_loc.y;
      res_loc.z   -= tmp_loc.z;
      res_loc.pot -= tmp_loc.pot;
#       else///ACCURATE_PARTIAL_SUM
      acc.x   += rx * rinv;
      acc.y   += ry * rinv;
      acc.z   += rz * rinv;
      acc.pot += r2 * rinv;/* if necessary */
#       endif//ACCURATE_PARTIAL_SUM
#else///PARTIAL_SUM_ACCELERATION
      ai->x   += rx * rinv;
      ai->y   += ry * rinv;
      ai->z   += rz * rinv;
      ai->pot += r2 * rinv;/* if necessary */
#endif//PARTIAL_SUM_ACCELERATION
      //-------------------------------------------------------------------
    }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  PARTIAL_SUM_ACCELERATION
  //-----------------------------------------------------------------------
#ifdef  ACCURATE_ACCUMULATION
  //-----------------------------------------------------------------------
  /* R := R + x_i */
  res->x   += acc.x;
  res->y   += acc.y;
  res->z   += acc.z;
  res->pot += acc.pot;
  /* T := S */
  acceleration tmp = *ai;
  /* S := S + R */
  ai->x   += res->x;
  ai->y   += res->y;
  ai->z   += res->z;
  ai->pot += res->pot;
  /* T := S - T */
  tmp.x   = ai->x   - tmp.x;
  tmp.y   = ai->y   - tmp.y;
  tmp.z   = ai->z   - tmp.z;
  tmp.pot = ai->pot - tmp.pot;
  /* R := R - T */
  res->x   -= tmp.x;
  res->y   -= tmp.y;
  res->z   -= tmp.z;
  res->pot -= tmp.pot;
  //-----------------------------------------------------------------------
#else///ACCURATE_ACCUMULATION
  //-----------------------------------------------------------------------
#ifdef  ACCURATE_PARTIAL_SUM
  acc.x   += res_loc.x;
  acc.y   += res_loc.y;
  acc.z   += res_loc.z;
  acc.pot += res_loc.pot;
#endif//ACCURATE_PARTIAL_SUM
  ai->x   += acc.x;
  ai->y   += acc.y;
  ai->z   += acc.z;
  ai->pot += acc.pot;
  //-----------------------------------------------------------------------
#endif//ACCURATE_ACCUMULATION
  //-----------------------------------------------------------------------
#endif//PARTIAL_SUM_ACCELERATION
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  COMPARE_WITH_DIRECT_SOLVER
__global__ void calcAccDirect_kernel
(position *ipos, acceleration * RESTRICT iacc, position *jpos, const int Nj
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
 , const real eps2_val
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
)
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* load poisition of an i-particle */
  //-----------------------------------------------------------------------
  const int idx = GLOBALIDX_X1D;
  position     pi = ipos[idx];
#ifndef INDIVIDUAL_GRAVITATIONAL_SOFTENING
  pi.m = eps2;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
  acceleration ai = {ZERO, ZERO, ZERO, ZERO};
  //-----------------------------------------------------------------------
#ifdef  ACCURATE_ACCUMULATION
  acceleration res = {ZERO, ZERO, ZERO, ZERO};
#endif//ACCURATE_ACCUMULATION
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const position massless = {ZERO, ZERO, ZERO, ZERO};
  //-----------------------------------------------------------------------
  __shared__ jparticle pj[NTHREADS * NLOOP];
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
  __shared__      real eps2[NTHREADS * NLOOP];
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
  //-----------------------------------------------------------------------
  for(int jj = 0; jj < Nj; jj += NTHREADS * NLOOP){
    //---------------------------------------------------------------------
    __syncthreads();
    for(int ll = 0; ll < NLOOP; ll++){
      //-------------------------------------------------------------------
      position pj_loc = (jj + NTHREADS * ll + tidx < Nj) ? jpos[jj + NTHREADS * ll + tidx] : massless;
      //-------------------------------------------------------------------
      jparticle pj_tmp;
      pj_tmp.x = pj_loc.x;
      pj_tmp.y = pj_loc.y;
      pj_tmp.z = pj_loc.z;
      pj_tmp.w = pj_loc.m;
      //-------------------------------------------------------------------
      pj  [NTHREADS * ll + tidx] = pj_tmp;
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
      eps2[NTHREADS * ll + tidx] = eps2_val;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
      //-------------------------------------------------------------------
    }
    __syncthreads();
    //---------------------------------------------------------------------
#pragma unroll
    for(int kk = 0; kk < NTHREADS * NLOOP; kk += TSUB * NLOOP)
#ifdef  IJ_PARALLELIZATION
#pragma unroll
      for(int ll = 0; ll < NWARP; ll++)
#endif//IJ_PARALLELIZATION
      calc_interaction
	(pi, &ai, &pj[kk]
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	 , &eps2[kk]
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  ACCURATE_ACCUMULATION
	 , &res
#endif//ACCURATE_ACCUMULATION
#ifdef  IJ_PARALLELIZATION
	 , ll
#endif//IJ_PARALLELIZATION
	 );
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#ifdef  ACCURATE_ACCUMULATION
  ai.x   += res.x;
  ai.y   += res.y;
  ai.z   += res.z;
  ai.pot += res.pot;
#endif//ACCURATE_ACCUMULATION
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* store acceleration of an i-particle from each thread */
  //-----------------------------------------------------------------------
#if 1
  iacc[idx] = ai;
#else
  atomicAdd(&(iacc[idx].x  ), ai.x  );
  atomicAdd(&(iacc[idx].y  ), ai.y  );
  atomicAdd(&(iacc[idx].z  ), ai.z  );
  atomicAdd(&(iacc[idx].pot), ai.pot);
#endif
  //-----------------------------------------------------------------------
}
#endif//COMPARE_WITH_DIRECT_SOLVER
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#if 1
//-------------------------------------------------------------------------
#include "buf_inc.cu"
//-------------------------------------------------------------------------
#else
//-------------------------------------------------------------------------
/*    tidx: thread index within a block */
/* freeNum: number of free buffers on global memory */
/* freeLst: list of index for the buffer on global memory */
/*  bufIdx: head index of buffer allocated on shared memory */
//-------------------------------------------------------------------------
#ifdef  USE_SMID_TO_GET_BUFID
//-------------------------------------------------------------------------
__device__ __forceinline__ uint getSMidx()
{
  //-----------------------------------------------------------------------
  uint rr;
  asm("mov.u32 %0, %%smid;" : "=r"(rr));
  return (rr);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ int occupyBuffer(const int tidx, uint *freeLst, uint *bufIdx)
{
  //-----------------------------------------------------------------------
  int target = 0;
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
#   if  NBLOCKS_PER_SM == 2
    target         = getSMidx() * NBLOCKS_PER_SM;
#else///NBLOCKS_PER_SM == 2
    const int head = getSMidx() * NBLOCKS_PER_SM;
#endif//NBLOCKS_PER_SM == 2
    //---------------------------------------------------------------------
    uint tmp = UINT_MAX;
    while( tmp == UINT_MAX ){
      //-------------------------------------------------------------------
#   if  NBLOCKS_PER_SM == 2
      //-------------------------------------------------------------------
      target ^= 1;
      tmp = atomicExch(&freeLst[target], UINT_MAX);
      //-------------------------------------------------------------------
#else///NBLOCKS_PER_SM == 2
      //-------------------------------------------------------------------
      for(int ii = 0; ii < NBLOCKS_PER_SM; ii++){
	//-----------------------------------------------------------------
	tmp = atomicExch(&freeLst[head + ii], UINT_MAX);
	//-----------------------------------------------------------------
	if( tmp != UINT_MAX ){
	  target = head + ii;
	  break;
	}/* if( tmp != UINT_MAX ){ */
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < NBLOCKS_PER_SM; ii++){ */
      //-------------------------------------------------------------------
#endif//NBLOCKS_PER_SM == 2
      //-------------------------------------------------------------------
    }/* while( tmp == UINT_MAX ){ */
    //---------------------------------------------------------------------
    bufIdx[0] = tmp;
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  return (target);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void releaseBuffer(const int tidx, uint *freeLst, uint bufIdx, const int target)
{
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    atomicExch(&freeLst[target], bufIdx);
    //---------------------------------------------------------------------
/* #   if  NBLOCKS_PER_SM == 2 */
/*     int target = getSMidx() * NBLOCKS_PER_SM; */
/* #else///NBLOCKS_PER_SM == 2 */
/*     const int head = getSMidx() * NBLOCKS_PER_SM; */
/* #endif//NBLOCKS_PER_SM == 2 */
/*     //--------------------------------------------------------------------- */
/*     while( bufIdx != UINT_MAX ){ */
/*       //------------------------------------------------------------------- */
/* #   if  NBLOCKS_PER_SM == 2 */
/*       //------------------------------------------------------------------- */
/*       bufIdx = atomicExch(&freeLst[target], bufIdx); */
/*       target ^= 1; */
/*       //------------------------------------------------------------------- */
/* #else///NBLOCKS_PER_SM == 2 */
/*       //------------------------------------------------------------------- */
/*       for(int ii = 0; ii < NBLOCKS_PER_SM; ii++){ */
/* 	//----------------------------------------------------------------- */
/* 	bufIdx = atomicExch(&freeLst[head + ii], bufIdx); */
/* 	//----------------------------------------------------------------- */
/* 	if( bufIdx == UINT_MAX ) */
/* 	  break; */
/* 	//----------------------------------------------------------------- */
/*       }/\* for(int ii = 0; ii < NBLOCKS_PER_SM; ii++){ *\/ */
/*       //------------------------------------------------------------------- */
/* #endif//NBLOCKS_PER_SM == 2 */
/*       //------------------------------------------------------------------- */
/*     }/\* while( bufIdx != UINT_MAX ){ *\/ */
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///USE_SMID_TO_GET_BUFID
//-------------------------------------------------------------------------
#ifdef  TRY_MODE_ABOUT_BUFFER
//-------------------------------------------------------------------------
__device__ __forceinline__ int occupyBuffer(const int tidx, const int bidx, const int bufNum, uint *freeLst, uint *bufIdx)
{
  //-----------------------------------------------------------------------
  int target = 0;
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    target = bidx % bufNum;
    //---------------------------------------------------------------------
    while( true ){
      //-------------------------------------------------------------------
      const uint tmp = atomicExch(&freeLst[target], UINT_MAX);
      //-------------------------------------------------------------------
      if( tmp != UINT_MAX ){
	bufIdx[0] = tmp;
	break;
      }/* if( tmp != UINT_MAX ){ */
      //-------------------------------------------------------------------
      target++;
      target %= bufNum;
      //-------------------------------------------------------------------
    }/* while( true ){ */
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  return (target);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void releaseBuffer(const int tidx, uint *freeLst, uint bufIdx, const int target)
{
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    atomicExch(&freeLst[target], bufIdx);
    //---------------------------------------------------------------------
    /* int target = bidx % bufNum; */
    /* //--------------------------------------------------------------------- */
    /* while( true ){ */
    /*   //------------------------------------------------------------------- */
    /*   const uint tmp = atomicExch(&freeLst[target], bufIdx); */
    /*   //------------------------------------------------------------------- */
    /*   if( tmp == UINT_MAX )	break; */
    /*   //------------------------------------------------------------------- */
    /*   bufIdx = tmp; */
    /*   target++; */
    /*   target %= bufNum; */
    /*   //------------------------------------------------------------------- */
    /* }/\* while( true ){ *\/ */
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///TRY_MODE_ABOUT_BUFFER
//-------------------------------------------------------------------------
__device__ __forceinline__ void  occupyBuffer(const int tidx, uint *freeNum, uint *freeLst, uint *bufIdx, int *active)
{
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    /* lock the shared array */
    while( true )
      if( atomicAnd(active, 0) )
	break;
    //---------------------------------------------------------------------
    /* pick up a free buffer */
    const uint target = atomicDec(freeNum, UINT_MAX) - 1;
    /* bufIdx [0]      = freeLst[target]; */
    /* freeLst[target] = UINT_MAX; */
    bufIdx[0] = atomicExch(&freeLst[target], UINT_MAX);
#ifdef  DBG_TREE_WALK
    printf("%u buffers are free, index of %u is in use\n", target + 1, bufIdx[0]);
#endif//DBG_TREE_WALK
    //---------------------------------------------------------------------
    /* release the shared array */
    atomicOr(active, 1);
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void releaseBuffer(const int tidx, uint *freeNum, uint *freeLst, const int bufIdx, int *active)
{
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  if( tidx == 0 ){
    //---------------------------------------------------------------------
    /* lock the shared array */
    while( true )
      if( atomicAnd(active, 0) )
	break;
    //---------------------------------------------------------------------
    /* release the buffer */
    const uint target = atomicInc(freeNum, UINT_MAX);
    freeLst[target] = (uint)bufIdx;
#ifdef  DBG_TREE_WALK
    printf("%u buffers are free, index of %u is released\n", target + 1, bufIdx);
#endif//DBG_TREE_WALK
    //---------------------------------------------------------------------
    /* release the shared array */
    atomicOr(active, 1);
    //---------------------------------------------------------------------
  }/* if( tidx == 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//TRY_MODE_ABOUT_BUFFER
//-------------------------------------------------------------------------
#endif//USE_SMID_TO_GET_BUFID
//-------------------------------------------------------------------------
#endif
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* calculate gravitational acceleration based on the width-first tree traversal */
//-------------------------------------------------------------------------
/* laneInfo :: input          :: head index and number of ``active'' i-particles */
/* ipos     :: input          :: position and mass of N-body particles */
/* iacc     :: input / output :: acceleration and potential of N-body particles */
/* more     :: input          :: head index and number of child particles of the corresponding j-particle */
/* root     :: input          :: index of the root tree node */
/* jpos     :: input          :: position and squared radius of pseudo N-body particle as j-particles */
/* mj       :: input          :: mass of pseudo N-body particle as j-particles */
/* active   ::                :: a shared value to lock the shared quantities (freeNum, freeLst) to control usage of buffer */
/* freeNum  ::                :: an unsigned integer represents # of unused bufferes */
/* freeLst  ::                :: a list of unused bufferes */
/* buffer   ::                :: tentative memory space to store tree cells which does not fit within the limited space of the shared memory */
/* bufSize  :: input          :: size of the buffer */
/* overflow ::         output :: a variable to detect buffer overflow */
//-------------------------------------------------------------------------
#   if  defined(ADOPT_SMALLEST_ENCLOSING_BALL) || defined(ADOPT_APPROXIMATED_ENCLOSING_BALL)
#include "../tree/seb_dev.cu"
#endif//defined(ADOPT_SMALLEST_ENCLOSING_BALL) || defined(ADOPT_APPROXIMATED_ENCLOSING_BALL)
//-------------------------------------------------------------------------
__global__ void __launch_bounds__(NTHREADS, NBLOCKS_PER_SM) calcAcc_kernel
     (READ_ONLY laneinfo * RESTRICT laneInfo, READ_ONLY position * RESTRICT ipos, jnode * RESTRICT iacc,
#ifdef  GADGET_MAC
      READ_ONLY acceleration * RESTRICT iacc_old,
#endif//GADGET_MAC
      const int root, READ_ONLY uint * RESTRICT more, READ_ONLY jparticle * RESTRICT jpos, READ_ONLY jmass * RESTRICT mj,
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
      int * RESTRICT active, uint * RESTRICT freeNum,
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  !defined(USE_SMID_TO_GET_BUFID) &&  defined(TRY_MODE_ABOUT_BUFFER)
      const int freeNum,
#endif//!defined(USE_SMID_TO_GET_BUFID) &&  defined(TRY_MODE_ABOUT_BUFFER)
      uint * RESTRICT freeLst, uint * RESTRICT buffer, const int bufSize, int * RESTRICT overflow
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
      , unsigned long long int * RESTRICT cycles
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  COUNT_INTERACTIONS
      , int * RESTRICT stockNj, int * RESTRICT stockNbuf
#endif//COUNT_INTERACTIONS
      )
{
  //-----------------------------------------------------------------------
  /* start stop watch */
  //-----------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) && !defined(PRINT_PSEUDO_PARTICLE_INFO)
  const long long int initCycle = clock64();
#endif//!defined(SERIALIZED_EXECUTION) && !defined(PRINT_PSEUDO_PARTICLE_INFO)
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int lane = tidx & (TSUB - 1);/* index of the thread within a thread group */
  //-----------------------------------------------------------------------
  const int head = tidx - lane;
  const int tail = head + (TSUB - 1);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* shared quantities in the thread parallelized version */
  //-----------------------------------------------------------------------
  __shared__ jnode   pj[NTHREADS * (NLOOP + 1)];
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
  __shared__ real  eps2[NTHREADS * (NLOOP + 1)];
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
  __shared__ uint queue[NTHREADS * NQUEUE];
  //-----------------------------------------------------------------------
  /* const int hq = lane + (head / TSUB) * TSUB * NQUEUE;/\* head index of the shared array close and queue within a thread group *\/ */
  /* const int hp =        (head / TSUB) * TSUB * (NLOOP + 1);/\* head index of the shared array pj within a thread group *\/ */
  const int hq = lane + DIV_TSUB(head) * TSUB * NQUEUE;/* head index of the shared array close and queue within a thread group */
  const int hp =        DIV_TSUB(head) * TSUB * (NLOOP + 1);/* head index of the shared array pj within a thread group */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* shared values within the threads */
  //-----------------------------------------------------------------------
  /* to store prefix sum */
#ifdef  USE_WARP_SHUFFLE_FUNC
  int smem;
#else///USE_WARP_SHUFFLE_FUNC
  __shared__ uint_real smem[NTHREADS];
#endif//USE_WARP_SHUFFLE_FUNC
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
  size_t buf0Head = (size_t)(bufIdx * NGROUPS + DIV_TSUB(head)) * (size_t)bufSize;
#ifdef  DBG_TREE_WALK
  if( lane == 0 )
    printf("buffer: %ld--%ld\n", buf0Head, buf0Head + bufSize - 1);
#endif//DBG_TREE_WALK
  //-----------------------------------------------------------------------
#ifdef  DBG_TREE_WALK
  real mjtot = ZERO;
#endif//DBG_TREE_WALK
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate gravitational force using hierarchical tree-structure */
  //-----------------------------------------------------------------------
  /* const int laneIdx = GLOBALIDX_X1D / TSUB; */
  const int laneIdx = DIV_TSUB(GLOBALIDX_X1D);
  const laneinfo info = laneInfo[laneIdx];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* load poisition of an i-particle */
  //-----------------------------------------------------------------------
#ifdef  IJ_PARALLELIZATION
  /* const bool skip = ((lane / NWARP) < info.num) ? (false) : (true); */
  const bool skip = (DIV_NWARP(lane) < info.num) ? (false) : (true);
  int      jtag =              lane & (NWARP - 1);
  /* const int idx = info.head + (lane /  NWARP); */
  const int idx = info.head + DIV_NWARP(lane);
#ifdef  DBG_TREE_WALK
  if( lane == 0 )
    printf("%d\t%d\t%d\t%d\t%d\n", GLOBALIDX_X1D, tidx, laneIdx, info.head, info.num);
#endif//DBG_TREE_WALK
#else///IJ_PARALLELIZATION
  const bool skip = (lane < info.num) ? (false) : (true);
  const int idx = info.head + lane;
#endif//IJ_PARALLELIZATION
  //-----------------------------------------------------------------------
  position   pi = {ZERO, ZERO, ZERO, UNITY};/* x, y, z, m */
  position icom = {ZERO, ZERO, ZERO, UNITY};/* x, y, z, m; m contains r2max */
  if( !skip ){
    /* load position and mass of i-particle from global memory */
    pi = ipos[idx];
#ifndef ADOPT_ENCLOSING_BALL
    icom = pi;
#endif//ADOPT_ENCLOSING_BALL
  }/* if( !skip ){ */
  //-----------------------------------------------------------------------
  int fail = hp + lane;
  jnode jidx;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
  const long long int initCycle = clock64();
#endif//PRINT_PSEUDO_PARTICLE_INFO
  //-----------------------------------------------------------------------
  /* set an enclosing sphere contains whole i-particles within TSUB threads */
  //-----------------------------------------------------------------------
#ifdef  ADOPT_ENCLOSING_BALL
  //-----------------------------------------------------------------------
  pj[hp + lane].pi = pi;
  icom = pj[hp].pi;
  if( skip )    pi = icom;
  //-----------------------------------------------------------------------
#ifdef  ADOPT_SMALLEST_ENCLOSING_BALL
  //-----------------------------------------------------------------------
  /* adopt the smallest enclosing ball */
  {
    pos4seb sebPos = {pi.x, pi.y, pi.z, false};
    real4 sebCen;
    findSEB(lane, &pj[hp], &sebPos, &sebCen, (real *)&pj[hp + TSUB], (real *)&pj[hp + TSUB + NDIM_SEB], (int *)&pj[hp + TSUB + 2 * NDIM_SEB], (real *)&pj[hp + TSUB + 2 * NDIM_SEB + 1]
#ifndef USE_WARP_SHUFFLE_FUNC
	    , smem, (int *)&queue[hq - tidx], tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
	    );
    icom.x = sebCen.x;
    icom.y = sebCen.y;
    icom.z = sebCen.z;
  }
  //-----------------------------------------------------------------------
#endif//ADOPT_SMALLEST_ENCLOSING_BALL
#ifdef  ADOPT_APPROXIMATED_ENCLOSING_BALL
  //-----------------------------------------------------------------------
  /* adopt the approximated enclosing ball proposed by Ritter (1990) */
  approxSEB(lane, &pj[hp], pi, &icom
#ifndef USE_WARP_SHUFFLE_FUNC
	    , smem, (int *)&queue[hq - tidx], tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
	    );
  //-----------------------------------------------------------------------
#endif//ADOPT_APPROXIMATED_ENCLOSING_BALL
#   if  !defined(ADOPT_SMALLEST_ENCLOSING_BALL) && !defined(ADOPT_APPROXIMATED_ENCLOSING_BALL)
  //-----------------------------------------------------------------------
  /* adopt a simple estimation of enclosing ball using the minimum bounding box in Cartesian coordinates */
  {
#ifdef  USE_WARP_SHUFFLE_FUNC
    const real xmin = getMinimumRealTsub(pi.x                  );    const real xmax = getMaximumRealTsub(pi.x                  );
    const real ymin = getMinimumRealTsub(pi.y                  );    const real ymax = getMaximumRealTsub(pi.y                  );
    const real zmin = getMinimumRealTsub(pi.z                  );    const real zmax = getMaximumRealTsub(pi.z                  );
#else///USE_WARP_SHUFFLE_FUNC
    const real xmin = getMinimumRealTsub(pi.x, smem, tidx, head);    const real xmax = getMaximumRealTsub(pi.x, smem, tidx, head);
    const real ymin = getMinimumRealTsub(pi.y, smem, tidx, head);    const real ymax = getMaximumRealTsub(pi.y, smem, tidx, head);
    const real zmin = getMinimumRealTsub(pi.z, smem, tidx, head);    const real zmax = getMaximumRealTsub(pi.z, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC
    icom.x = HALF * (xmin + xmax);
    icom.y = HALF * (ymin + ymax);
    icom.z = HALF * (zmin + zmax);
  }
#ifdef  COMPARE_ENCLOSING_BALLS
  position ball = icom;
  icom = pi;
#endif//COMPARE_ENCLOSING_BALLS
  //-----------------------------------------------------------------------
#endif//!defined(ADOPT_SMALLEST_ENCLOSING_BALL) && !defined(ADOPT_APPROXIMATED_ENCLOSING_BALL)
  //-----------------------------------------------------------------------
#endif//ADOPT_ENCLOSING_BALL
#   if  !defined(ADOPT_ENCLOSING_BALL) || defined(COMPARE_ENCLOSING_BALLS)
  //-----------------------------------------------------------------------
  /* calculate center-of-mass of a group of i-particles as an enclosing sphere */
  icom.x *= icom.m;
  icom.y *= icom.m;
  icom.z *= icom.m;
  pj[fail].pi = icom;
  /* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
#   if  TSUB >=  2
#                    if  NWARP <  2
  jidx = pj[fail ^  1];  icom.x += jidx.pi.x;  icom.y += jidx.pi.y;  icom.z += jidx.pi.z;  icom.m += jidx.pi.m;  pj[fail].pi = icom;
#                 endif//NWARP <  2
#   if  TSUB >=  4
#                    if  NWARP <  4
  jidx = pj[fail ^  2];  icom.x += jidx.pi.x;  icom.y += jidx.pi.y;  icom.z += jidx.pi.z;  icom.m += jidx.pi.m;  pj[fail].pi = icom;
#                 endif//NWARP <  4
#   if  TSUB >=  8
#                    if  NWARP <  8
  jidx = pj[fail ^  4];  icom.x += jidx.pi.x;  icom.y += jidx.pi.y;  icom.z += jidx.pi.z;  icom.m += jidx.pi.m;  pj[fail].pi = icom;
#                 endif//NWARP <  8
#   if  TSUB >= 16
#                    if  NWARP < 16
  jidx = pj[fail ^  8];  icom.x += jidx.pi.x;  icom.y += jidx.pi.y;  icom.z += jidx.pi.z;  icom.m += jidx.pi.m;  pj[fail].pi = icom;
#                 endif//NWARP < 16
#   if  TSUB == 32
#                    if  NWARP < 32
  jidx = pj[fail ^ 16];  icom.x += jidx.pi.x;  icom.y += jidx.pi.y;  icom.z += jidx.pi.z;  icom.m += jidx.pi.m;  pj[fail].pi = icom;
#                 endif//NWARP < 32
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
  icom = pj[hp].pi;/* icom.m = Mtot */
  icom.m = UNITY / icom.m;/* tentative use as minv */
  icom.x *= icom.m;
  icom.y *= icom.m;
  icom.z *= icom.m;
  //-----------------------------------------------------------------------
#ifdef  COMPARE_ENCLOSING_BALLS
  {
    real dx, dy, dz;
    /* calculate size of geometrical estimated enclosing ball */
    dx = pi.x - ball.x;
    dy = pi.y - ball.y;
    dz = pi.z - ball.z;
    ball.m = getMaximumRealTsub(dx * dx + dy * dy + dz * dz
#ifndef USE_WARP_SHUFFLE_FUNC
				, smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
				);
    /* calculate size of mass-weighted estimated enclosing ball */
    dx = pi.x - icom.x;
    dy = pi.y - icom.y;
    dz = pi.z - icom.z;
    icom.m = getMaximumRealTsub(dx * dx + dy * dy + dz * dz
#ifndef USE_WARP_SHUFFLE_FUNC
				, smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
				);
    /* adopt smaller one */
    if( ball.m < icom.m )
      icom = ball;
  }
#endif//COMPARE_ENCLOSING_BALLS
  //-----------------------------------------------------------------------
#endif//!defined(ADOPT_ENCLOSING_BALL) || defined(COMPARE_ENCLOSING_BALLS)
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate radius of the sphere which include whole i-particles within the group centered on (icom.x, icom.y, icom.z) */
  //-----------------------------------------------------------------------
#ifdef  GADGET_MAC
  real amin;
#endif//GADGET_MAC
  {
    //---------------------------------------------------------------------
    /* calculate displacement of i-particle and center-of-mass */
    const real rx = pi.x - icom.x;
    const real ry = pi.y - icom.y;
    const real rz = pi.z - icom.z;
    //---------------------------------------------------------------------
    /* calculate maximum of r squared */
    icom.m = getMaximumRealTsub(1.0e-30f + rx * rx + ry * ry + rz * rz
#ifndef USE_WARP_SHUFFLE_FUNC
				, smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
				);
    //---------------------------------------------------------------------
#ifdef  GADGET_MAC
    //---------------------------------------------------------------------
    acceleration ai_old = {ZERO, ZERO, ZERO, ZERO};
    if( !skip )
      ai_old = iacc_old[idx];
    //---------------------------------------------------------------------
    /* calculate minimum of a squared */
    amin = getMinimumRealTsub(1.0e-30f + ai_old.x * ai_old.x + ai_old.y * ai_old.y + ai_old.z * ai_old.z
#ifndef USE_WARP_SHUFFLE_FUNC
			      , smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
			      );
    amin *= RSQRT(amin);
    //---------------------------------------------------------------------
#endif//GADGET_MAC
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#ifndef INDIVIDUAL_GRAVITATIONAL_SOFTENING
  /* set square of softening length at pi.m */
  pi.m = eps2;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
  //-----------------------------------------------------------------------
#ifdef  DBG_TREE_WALK
  if( lane == 0 )
    printf("%d: %f, %f, %f w/ r2max = %e\n", GLOBALIDX_X1D, icom.x, icom.y, icom.z, icom.m);
#endif//DBG_TREE_WALK
  //-----------------------------------------------------------------------
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
  long long int exitCycle = clock64();
  if( lane == 0 ){
    unsigned long long int elapsed = (unsigned long long int)(exitCycle - initCycle);
    atomicAdd(cycles, elapsed);
    iacc[DIV_TSUB(GLOBALIDX_X1D)].pi = icom;
  }/* if( tidx == 0 ){ */
#endif//PRINT_PSEUDO_PARTICLE_INFO
  //-----------------------------------------------------------------------



  //-----------------------------------------------------------------------
#ifndef PRINT_PSEUDO_PARTICLE_INFO
  //-----------------------------------------------------------------------
  /* sweep all j-cells by executing tree-traversal */
  //-----------------------------------------------------------------------
  /* initialize queue for j-cells */
#pragma unroll
  for(int jj = 0; jj < NQUEUE; jj++)
    queue[hq + TSUB * jj] = NULL_NODE;
  //-----------------------------------------------------------------------
  /* initialize queue for j-cells and interaction list by a representative thread */
  int Nj = 0;
  int bufHead = 0;
  int bufTail = 0;
  int bufOpen = bufSize;
  int bufUsed = 0;
  /* set child j-cells in queue on the shared memory */
  uint jcell = more[root];
  int rem = 1 + (jcell >> IDXBITS);
  jcell &= IDXMASK;
  //-----------------------------------------------------------------------
  if( rem > TSUB ){
    //---------------------------------------------------------------------
    /* if rem exceeds TSUB, then number of child j-cells must be shrunk */
    //---------------------------------------------------------------------
    queue[hq] = jcell + lane;
    //---------------------------------------------------------------------
    if( tidx == tail )
      queue[hq] += ((rem - TSUB) << IDXBITS);
    //---------------------------------------------------------------------
    rem = TSUB;
    //---------------------------------------------------------------------
  }/* if( rem > TSUB ){ */
  else{
    //---------------------------------------------------------------------
    /* upload rem (<= TSUB) child j-cells to the shared memory */
    //---------------------------------------------------------------------
    if( lane < rem )
      queue[hq] = more[jcell + lane];
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
  if( info.num == 0 )
    rem = 0;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* tree traversal in a width-first manner */
  //-----------------------------------------------------------------------
  acceleration  ai = {ZERO, ZERO, ZERO, ZERO};/* ax, ay, az, pot */
#ifdef  ACCURATE_ACCUMULATION
  acceleration res = {ZERO, ZERO, ZERO, ZERO};
#endif//ACCURATE_ACCUMULATION
  //-----------------------------------------------------------------------
  fail = 0;
  //-----------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
  int Nj_tot = 0;
  int Nb_max = 0;
#endif//COUNT_INTERACTIONS
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    /* if the queue becomes empty, then exit the while loop */
    //---------------------------------------------------------------------
    if( rem == 0 )
      break;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* pick up a queue from stack */
    //---------------------------------------------------------------------
    /* initialize the shared memory */
    uint leaf = 0;
#pragma unroll
    for(int iter = 0; iter < NSTOCK; iter++)
      jidx.idx[iter] = NULL_NODE;
    pj[hp + lane + NLOOP * TSUB] = jidx;
    //---------------------------------------------------------------------
    /* tentative load from the stack */
    int cnum = 0;
    jcell = NULL_NODE;
    if( lane < rem ){
      jcell = queue[hq];
      cnum = 1 + (jcell >> IDXBITS);
    }/* if( lane < rem ){ */
    //---------------------------------------------------------------------
    /* predict the head index on the shared memory by parallel prefix sum */
#ifdef  USE_WARP_SHUFFLE_FUNC
    smem = prefixSumTsub(cnum, lane);
    int hidx = smem         - cnum;/* exclusive prefix sum of cnum */
#else///USE_WARP_SHUFFLE_FUNC
    prefixSumTsub(cnum, smem, tidx, lane);
    int hidx = smem[tidx].i - cnum;/* exclusive prefix sum of cnum */
#endif//USE_WARP_SHUFFLE_FUNC
    //---------------------------------------------------------------------
    int remove = 0;
    if( (cnum != 0) && (hidx < TSUB * NSTOCK) ){
      //-------------------------------------------------------------------
      /* local data can be uploaded to the shared memory */
      int unum = TSUB * NSTOCK - hidx;
      if( cnum < unum )	  unum = cnum;
      //-------------------------------------------------------------------
      /* upload local data */
      jcell &= IDXMASK;
      for(int jj = 0; jj < unum; jj++){
	pj[hp + NLOOP * TSUB + (hidx & (TSUB - 1))].idx[DIV_TSUB(hidx)] = jcell;/* assumes TSUB is a power of 2 */
	hidx++;
	jcell++;
      }/* for(int jj = 0; jj < unum; jj++){ */
      //-------------------------------------------------------------------
      /* eliminate stocked j-cells from the queue */
      if( unum == cnum )
	remove = 1;
      else{
	jcell += ((cnum - unum - 1) << IDXBITS);
	queue[hq] = jcell;
      }/* else{ */
      //-------------------------------------------------------------------
    }/* if( (cnum != 0) && (hidx < TSUB * NSTOCK) ){ */
    //---------------------------------------------------------------------
    /* remove scanned j-cells if possible */
#ifdef  USE_WARP_SHUFFLE_FUNC
    smem = prefixSumTsub(remove, lane);
    remove = __shfl(smem, TSUB - 1, TSUB);
#else///USE_WARP_SHUFFLE_FUNC
    prefixSumTsub(remove, smem, tidx, lane);
    remove = smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC
    //---------------------------------------------------------------------
    if( remove != 0 ){
      rem -= remove;
      copyData_s2s(queue, hq + remove, queue, hq, rem, lane);
    }/* if( remove != 0 ){ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* pick up pseudo particles from NSTOCK buffers */
    //---------------------------------------------------------------------
    jidx = pj[hp + lane + NLOOP * TSUB];
    //---------------------------------------------------------------------
#pragma unroll
    for(int iter = 0; iter < NSTOCK; iter++){
      //-------------------------------------------------------------------
      /* set an index of j-cell */
      const uint target = jidx.idx[iter];
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* only the active threads pick up a j-cell from the global memory */
      //-------------------------------------------------------------------
      jparticle jcnd;
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
      real      jeps2;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
      int calc = 0;
      //-------------------------------------------------------------------
      if( target != NULL_NODE ){
	//-----------------------------------------------------------------
	jcnd = jpos[target];
	//-----------------------------------------------------------------
	/* set a pseudo i-particle */
	const real rx = jcnd.x - icom.x;
	const real ry = jcnd.y - icom.y;
	const real rz = jcnd.z - icom.z;
	const real r2 = rx * rx + ry * ry + rz * rz;
#if 1
	real lambda = FMAX(UNITY - SQRTRATIO(icom.m, r2), ZERO);
#else
	real lambda = UNITY - SQRTRATIO(icom.m, r2);
	if( lambda < EPSILON )	    lambda = ZERO;
#endif
	/* calculate distance between the pseudo i-particle and the candidate j-particle */
	//-----------------------------------------------------------------
#ifdef  GADGET_MAC
	/* alpha * |a| * r^4 > G * M * l^2 */
#if 1
	lambda *= lambda * r2;
	if( jcnd.w < lambda               * lambda               * amin )
#else
	if( jcnd.w < lambda * lambda * r2 * lambda * lambda * r2 * amin )
#endif
#else///GADGET_MAC
#ifdef  WS93_MAC
	  if(   jcnd.w < lambda * lambda * r2 )
#else///WS93_MAC
	    /* (l / r) < theta */
	    if( jcnd.w < lambda * lambda * r2 * theta2 )
#endif//WS93_MAC
#endif//GADGET_MAC
	      {
		/* add the candidate j-particle to the interaction list */
		const jmass mj_tmp = mj[target];
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
		jcnd.w = mj_tmp.mass;
		jeps2  = mj_tmp.eps2;
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
		jcnd.w = mj_tmp;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
		calc = 1;
	      }
	    else
	      {
		/* add child-cells of near tree-cells to the tentative stack */
		leaf += (1 << (IDX_SHIFT_BITS * iter));
		jidx.idx[iter] = more[target];
	      }
	//-----------------------------------------------------------------
      }/* if( target != NULL_NODE ){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* prefixSum to build a local interaction list */
      //-------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC
      smem = prefixSumTsub(calc, lane);
      hidx = smem         - calc;/* exclusive prefix sum of calc */
#else///USE_WARP_SHUFFLE_FUNC
      prefixSumTsub(calc, smem, tidx, lane);
      hidx = smem[tidx].i - calc;/* exclusive prefix sum of calc */
#endif//USE_WARP_SHUFFLE_FUNC
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* add distant tree-cells to the interaction list */
      if( calc ){
	pj  [hp + Nj + hidx].pos = jcnd;
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	eps2[hp + Nj + hidx] = jeps2;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
      }
#ifdef  USE_WARP_SHUFFLE_FUNC
      Nj += __shfl(smem, TSUB - 1, TSUB);/* inclusive prefix sum of calc */
#else///USE_WARP_SHUFFLE_FUNC
      Nj += smem[tail].i;/* inclusive prefix sum of calc */
#endif//USE_WARP_SHUFFLE_FUNC
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* calculate body--body interaction if sufficient size of interaction list is available */
      //-------------------------------------------------------------------
      if( Nj >= NLOOP * TSUB ){
	//-----------------------------------------------------------------
	calc_interaction
	  (pi, &ai, (jparticle *)&pj[hp]
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	   , &eps2[hp]
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  ACCURATE_ACCUMULATION
	   , &res
#endif//ACCURATE_ACCUMULATION
#ifdef  IJ_PARALLELIZATION
	   , jtag
#endif//IJ_PARALLELIZATION
	   );
	//-----------------------------------------------------------------
#ifdef  DBG_TREE_WALK
	for(int ll = 0; ll < NLOOP * TSUB; ll++)
	  mjtot += pj[hp + ll].pos.w;
#endif//DBG_TREE_WALK
	//-----------------------------------------------------------------
	pj  [hp + lane] = pj  [hp + lane + NLOOP * TSUB];
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	eps2[hp + lane] = eps2[hp + lane + NLOOP * TSUB];
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
	Nj -= NLOOP * TSUB;
	//-----------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
	Nj_tot += NLOOP * TSUB;
#endif//COUNT_INTERACTIONS
	//-----------------------------------------------------------------
      }/* if( Nj >= NLOOP * TSUB ){ */
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < NSTOCK; iter++){ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* if the shared memory has open space and some tree cells are stored on the global memory, then load tree-cells from the global memory to the shared memory */
    //---------------------------------------------------------------------
    /* evaluate available size of the queue on the shared memory */
    int Nsm_rem = NQUEUE * TSUB - rem;
    //---------------------------------------------------------------------
    if(  (bufUsed != 0) && (Nsm_rem > 0) ){
      //-------------------------------------------------------------------
      const int Nload = (Nsm_rem < bufUsed) ? (Nsm_rem) : (bufUsed);
#if 0
      memcpy((void *)&queue[hq - lane + rem], (const void *)&buffer[buf0Head + (size_t)bufHead], sizeof(uint) * Nload);
#else
      copyData_g2s(buffer, buf0Head + (size_t)bufHead, queue, hq - lane + rem, Nload, lane);
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
#ifdef  USE_WARP_SHUFFLE_FUNC
    cpChildNodes(      (uint *)(&pj[hp + NLOOP * TSUB]), jidx, leaf,       lane,       queue, hq, &Nsm_rem, &rem, buffer, buf0Head, &bufOpen, &bufUsed, &bufHead, &bufTail);
#else///USE_WARP_SHUFFLE_FUNC
    cpChildNodes(smem, (uint *)(&pj[hp + NLOOP * TSUB]), jidx, leaf, tidx, lane, tail, queue, hq, &Nsm_rem, &rem, buffer, buf0Head, &bufOpen, &bufUsed, &bufHead, &bufTail);
#endif//USE_WARP_SHUFFLE_FUNC
    //---------------------------------------------------------------------
    fail += (bufOpen < 0);
#ifdef  COUNT_INTERACTIONS
    if( bufUsed > Nb_max )
      Nb_max = bufUsed;
#endif//COUNT_INTERACTIONS
    //---------------------------------------------------------------------
#if 0
    if( fail > 0 )
      break;
#endif
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate body--body interaction for remained j-particles */
  //-----------------------------------------------------------------------
  if( Nj != 0 ){
    //---------------------------------------------------------------------
    /* add massless particles at the tail of the interaction list */
    const int Ndummy = NLOOP * TSUB - Nj;
    const jparticle massless = {ZERO, ZERO, ZERO, ZERO};
    for(int jj = 0; jj < NLOOP; jj++){
      const int addr = lane + jj * TSUB;
      if( addr < Ndummy ){
	pj  [hp + Nj + addr].pos = massless;
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	eps2[hp + Nj + addr]     = UNITY;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
      }/* if( addr < Ndummy ){ */
    }/* for(int jj = 0; jj < NLOOP; jj++){ */
    //---------------------------------------------------------------------
    calc_interaction
      (pi, &ai, (jparticle *)&pj[hp]
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
       , &eps2[hp]
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  ACCURATE_ACCUMULATION
       , &res
#endif//ACCURATE_ACCUMULATION
#ifdef  IJ_PARALLELIZATION
       , jtag
#endif//IJ_PARALLELIZATION
       );
#ifdef  DBG_TREE_WALK
    for(int ll = 0; ll < NLOOP * TSUB; ll++)
      mjtot += pj[hp + ll].pos.w;
#endif//DBG_TREE_WALK
#ifdef  COUNT_INTERACTIONS
    Nj_tot += NLOOP * TSUB;
#endif//COUNT_INTERACTIONS
    //---------------------------------------------------------------------
  }/* if( Nj != 0 ){ */
  //-----------------------------------------------------------------------
  /* force accumulation */
#ifdef  ACCURATE_ACCUMULATION
  ai.x   += res.x;
  ai.y   += res.y;
  ai.z   += res.z;
  ai.pot += res.pot;
#endif//ACCURATE_ACCUMULATION
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* store acceleration of an i-particle from each thread */
  //-----------------------------------------------------------------------
  if( !skip ){
    //---------------------------------------------------------------------
    /* NOTE: implicit synchronization within 32 threads (a warp) is assumed for NWARP = 8, 16, 32 */
#   if  defined(SERIALIZED_EXECUTION) && (NWARP == 1)
    iacc[idx].ai = ai;
#else///defined(SERIALIZED_EXECUTION) && (NWARP == 1)
#          if  NWARP > 1
    pj[hp + lane].ai = ai;
    /* set index to accumulate acceleration without atomic operations */
#   if  NWARP > 4
    const int gtag = jtag >> 2;
    jtag &= 3;
#endif//NWARP > 4
    const int itag = hp + lane - jtag;
#   if  NWARP == 2
    atomicAdd(&(iacc[idx].val[jtag]), pj[itag].val[jtag] + pj[itag + 1].val[jtag]);    jtag ^= 2;
    atomicAdd(&(iacc[idx].val[jtag]), pj[itag].val[jtag] + pj[itag + 1].val[jtag]);
#endif//NWARP == 2
#   if  NWARP == 4
    atomicAdd(&(iacc[idx].val[jtag]), pj[itag].val[jtag] + pj[itag + 1].val[jtag] + pj[itag + 2].val[jtag] + pj[itag + 3].val[jtag]);
#endif//NWARP == 4
#   if  NWARP >= 8
    pj[itag].val[jtag] += pj[itag + 1].val[jtag] + pj[itag + 2].val[jtag] + pj[itag + 3].val[jtag];
#   if  NWARP == 32
    if( gtag < 4 )
      pj[itag].val[jtag] += pj[itag + 16].val[jtag];
#endif//NWARP == 32
#   if  NWARP >= 16
    if( gtag < 2 )
      pj[itag].val[jtag] += pj[itag + 8].val[jtag];
#endif//NWARP >= 16
    if( gtag == 0 )
      atomicAdd(&(iacc[idx].val[jtag]), pj[itag].val[jtag] + pj[itag + 4].val[jtag]);
#endif//NWARP >= 8
#       else///NWARP > 1
    atomicAdd(&(iacc[idx].ai.x  ), ai.x  );
    atomicAdd(&(iacc[idx].ai.y  ), ai.y  );
    atomicAdd(&(iacc[idx].ai.z  ), ai.z  );
    atomicAdd(&(iacc[idx].ai.pot), ai.pot);
#       endif//NWARP > 1
#endif//defined(SERIALIZED_EXECUTION) && (NWARP == 1)
    //---------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
    atomicAdd(&(stockNj  [idx]), Nj_tot);
    atomicAdd(&(stockNbuf[idx]), Nb_max);
#endif//COUNT_INTERACTIONS
    //---------------------------------------------------------------------
    if( tidx == head )
      atomicAdd(overflow, fail);
    //---------------------------------------------------------------------
#ifdef  DBG_TREE_WALK
    if( tidx == head )
      printf("mjtot is %e for %d-th thread\n", mjtot, tidx);
#endif//DBG_TREE_WALK
    //---------------------------------------------------------------------
  }/* if( !skip ){ */
  //-----------------------------------------------------------------------
#endif//PRINT_PSEUDO_PARTICLE_INFO
  //-----------------------------------------------------------------------


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
#   if  !defined(SERIALIZED_EXECUTION) && !defined(PRINT_PSEUDO_PARTICLE_INFO)
  long long int exitCycle = clock64();
  if( tidx == 0 ){
    unsigned long long int elapsed = (unsigned long long int)(exitCycle - initCycle);
    atomicAdd(cycles, elapsed);
  }/* if( tidx == 0 ){ */
#endif//!defined(SERIALIZED_EXECUTION) && !defined(PRINT_PSEUDO_PARTICLE_INFO)
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* calculate gravitational acceleration and potential */
//-------------------------------------------------------------------------
/* Ni     :: input          :: total number of N-body particles stored in this process */
/* ipos   :: input          :: position and mass of N-body particles */
/* iacc   ::         output :: acceleration and potential of N-body particles */
/* more   :: input          :: head index and number of child particles of the corresponding j-particle */
/* jpos   :: input          :: position and squared radius of pseudo N-body particle as j-particles */
/* mj     :: input          :: mass of pseudo N-body particle as j-particles */
/* buffer ::                :: tentative memory space to store tree cells which does not fit within the limited space of the shared memory */
//-------------------------------------------------------------------------
static inline void callCalcGravityFunc
(const dim3 blck, const dim3 thrd, kernelStream *sinfo, int *sidx,
 laneinfo * RESTRICT laneInfo, const iparticle pi, const int rootIdx, const soaTreeNode tree
#ifndef SERIALIZED_EXECUTION
 , const int jhead
#endif//SERIALIZED_EXECUTION
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
 , unsigned long long int * RESTRICT cycles
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
 , const soaTreeWalkBuf buf
#ifdef  COUNT_INTERACTIONS
 , iparticle_treeinfo treeInfo
#endif//COUNT_INTERACTIONS
 )
{
  //-----------------------------------------------------------------------
#   if  defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
  if( blck.x != 0 ){
#endif//defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
    //---------------------------------------------------------------------
    if( blck.x <= MAX_BLOCKS_PER_GRID )
      calcAcc_kernel<<<blck, thrd, SMEM_SIZE, sinfo->stream[*sidx]>>>
	(laneInfo,
#ifdef  BLOCK_TIME_STEP
	 pi.jpos,
#else///BLOCK_TIME_STEP
	 pi.pos,
#endif//BLOCK_TIME_STEP
	 (jnode *)pi.acc,
#ifdef  GADGET_MAC
	 pi.acc_old,
#endif//GADGET_MAC
	 rootIdx,
#ifdef  SERIALIZED_EXECUTION
	 tree.more, tree.jpos, tree.mj,
#else///SERIALIZED_EXECUTION
	 &tree.more[jhead], &tree.jpos[jhead], &tree.mj[jhead],
#endif//SERIALIZED_EXECUTION
#ifndef USE_SMID_TO_GET_BUFID
#ifndef TRY_MODE_ABOUT_BUFFER
	 buf.active,
#endif//TRY_MODE_ABOUT_BUFFER
	 buf.freeNum,
#endif//USE_SMID_TO_GET_BUFID
	 buf.freeLst, buf.buffer, buf.bufSize, buf.fail
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
	 , cycles
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  COUNT_INTERACTIONS
	 , treeInfo.Nj, treeInfo.Nbuf
#endif//COUNT_INTERACTIONS
	 );
    //---------------------------------------------------------------------
    else{
      //-------------------------------------------------------------------
      int Nrem = blck.x;
      const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
      int hidx = 0;
      //-------------------------------------------------------------------
      for(int iter = 0; iter < Niter; iter++){
	//-----------------------------------------------------------------
	int Nblck = MAX_BLOCKS_PER_GRID;
	if( Nblck > Nrem )	  Nblck = Nrem;
	//-----------------------------------------------------------------
	int Nsub = Nblck * NGROUPS;
	calcAcc_kernel<<<Nblck, thrd.x, SMEM_SIZE, sinfo->stream[*sidx]>>>
	  (&laneInfo[hidx],
#ifdef  BLOCK_TIME_STEP
	   pi.jpos,
#else///BLOCK_TIME_STEP
	   pi.pos,
#endif//BLOCK_TIME_STEP
	   (jnode *)pi.acc,
#ifdef  GADGET_MAC
	   pi.acc_old,
#endif//GADGET_MAC
	   rootIdx,
#ifdef  SERIALIZED_EXECUTION
	   tree.more, tree.jpos, tree.mj,
#else///SERIALIZED_EXECUTION
	   &tree.more[jhead], &tree.jpos[jhead], &tree.mj[jhead],
#endif//SERIALIZED_EXECUTION
#ifndef USE_SMID_TO_GET_BUFID
#ifndef TRY_MODE_ABOUT_BUFFER
	   buf.active,
#endif//TRY_MODE_ABOUT_BUFFER
	   buf.freeNum,
#endif//USE_SMID_TO_GET_BUFID
	   buf.freeLst, buf.buffer, buf.bufSize, buf.fail
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
	   , cycles
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  COUNT_INTERACTIONS
	   , treeInfo.Nj, treeInfo.Nbuf
#endif//COUNT_INTERACTIONS
	   );
	//-----------------------------------------------------------------
	hidx += Nsub;
	Nrem -= Nblck;
	//-----------------------------------------------------------------
      }/* for(int iter = 0; iter < Niter; iter++){ */
      //-------------------------------------------------------------------
    }/* else{ */
    //---------------------------------------------------------------------
/* #ifndef SERIALIZED_EXECUTION */
/*     /\* evaluate GPU time based on clock cycle counter *\/ */
/*     checkCudaErrors(cudaMemcpyAsync(&clockCycles, *cycles, sizeof(unsigned long long int), cudaMemcpyDeviceToHost, sinfo->stream[*sidx])); */
/* #endif//SERIALIZED_EXECUTION */
    //---------------------------------------------------------------------
    *sidx ^= 1;
    //---------------------------------------------------------------------
#   if  defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
  }/* if( blck.x != 0 ){ */
#endif//defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
void calcGravity_dev
(const int grpNum
#ifdef  BLOCK_TIME_STEP
 , double *reduce, const int totNum
#endif//BLOCK_TIME_STEP
 , laneinfo * RESTRICT laneInfo, const int Ni, const iparticle pi, const soaTreeNode tree, const soaTreeWalkBuf buf
 , kernelStream *sinfo, deviceProp devProp, double *time
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
 , char *file
#endif//PRINT_PSEUDO_PARTICLE_INFO
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
 , unsigned long long int *cycles_hst, unsigned long long int *cycles_dev
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifndef SERIALIZED_EXECUTION
 , double *twalk, const int pjNum
#ifdef  LET_COMMUNICATION_VIA_HOST
 , const soaTreeNode tree_hst
#endif//LET_COMMUNICATION_VIA_HOST
 , const int Nlet, domainInfo *let, const int Nstream_let, cudaStream_t stream_let[], MPIcfg_tree mpi
#ifdef  MONITOR_LETGEN_TIME
 , double *tlet, unsigned long long int *cycles_let_hst, unsigned long long int *cycles_let_dev
#endif//MONITOR_LETGEN_TIME
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
 , iparticle_treeinfo treeInfo
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
#ifdef  COMPARE_WITH_DIRECT_SOLVER
 , const bool approxGravity
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
 , const real eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#endif//COMPARE_WITH_DIRECT_SOLVER
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  int Nrem;
  //-----------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
  /* initialize count of Nj and Nbuf */
  initCounter_kernel<<<BLOCKSIZE(Ni, NTHREADS), NTHREADS>>>(treeInfo.Nj, treeInfo.Nbuf);
  getLastCudaError("initCounter_kernel");
#endif//COUNT_INTERACTIONS
  //-----------------------------------------------------------------------
#if 0
  fprintf(stderr, "grpNum = %d\n", grpNum);
  exit(0);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set thread-block configuration */
  static dim3 thrd, blck;
  thrd.x =                           NTHREADS;  thrd.y = 1;  thrd.z = 1;
  /* blck.x = BLOCKSIZE(grpNum * NWARP, NGROUPS);  blck.y = 1;  blck.z = 1; */
  blck.x = BLOCKSIZE(grpNum, NGROUPS);  blck.y = 1;  blck.z = 1;
  //-----------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  /* initialize clock cycle counter */
  *cycles_hst = 0;
  checkCudaErrors(cudaMemcpy(cycles_dev, cycles_hst, sizeof(unsigned long long int), cudaMemcpyHostToDevice));
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  *cycles_let_hst = 0;
  checkCudaErrors(cudaMemcpy(cycles_let_dev, cycles_let_hst, sizeof(unsigned long long int), cudaMemcpyHostToDevice));
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  //-----------------------------------------------------------------------
  static struct timeval start;
  checkCudaErrors(cudaDeviceSynchronize());
  gettimeofday(&start, NULL);
  /* initStopwatch_dev(devInfo); */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* initialize acceleration and potential */
  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
  Nrem = BLOCKSIZE(grpNum, NWARP * NGROUPS);
#else///BLOCK_TIME_STEP
  Nrem = BLOCKSIZE(Ni, NTHREADS);
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
  /* when grid splitting is not required... */
  if( Nrem <= MAX_BLOCKS_PER_GRID ){
#ifdef  BLOCK_TIME_STEP
#ifndef SERIALIZED_EXECUTION
    if( grpNum != 0 )
#endif//SERIALIZED_EXECUTION
      initAcc_kernel<<<Nrem, thrd>>>
	(pi.acc, BLOCKSIZE(grpNum, NGROUPS) * NGROUPS, laneInfo
#ifdef  GADGET_MAC
	 , pi.acc_old
#endif//GADGET_MAC
	 );
#else///BLOCK_TIME_STEP
    initAcc_kernel<<<Nrem, NTHREADS>>>
      (pi.acc
#ifdef  GADGET_MAC
       , pi.acc_old
#endif//GADGET_MAC
       );
#endif//BLOCK_TIME_STEP
  }/* if( Nrem <= MAX_BLOCKS_PER_GRID ){ */
  //-----------------------------------------------------------------------
  /* when grid splitting is required... */
  else{
    //---------------------------------------------------------------------
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;
    //---------------------------------------------------------------------
    for(int iter = 0; iter < Niter; iter++){
      //-------------------------------------------------------------------
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;
      //-------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
      int Nsub = Nblck * NWARP * NGROUPS;
      initAcc_kernel<<<Nblck, thrd.x>>>
	(pi.acc, BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx]
#ifdef  GADGET_MAC
	 , pi.acc_old
#endif//GADGET_MAC
	 );
#else///BLOCK_TIME_STEP
      int Nsub = Nblck * NTHREADS;
      initAcc_kernel<<<Nblck, NTHREADS>>>
	(&pi.acc[hidx]
#ifdef  GADGET_MAC
	 , &pi.acc_old[hidx]
#endif//GADGET_MAC
	 );
#endif//BLOCK_TIME_STEP
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
  getLastCudaError("initAcc_kernel");
  //-----------------------------------------------------------------------
  /* calculate gravitational acceleration based on the width-first tree traversal */
  //-----------------------------------------------------------------------
#ifdef  COMPARE_WITH_DIRECT_SOLVER
  if( approxGravity )
#endif//COMPARE_WITH_DIRECT_SOLVER
    {
      //-------------------------------------------------------------------
      /* gravity from j-particles within local process */
      //-------------------------------------------------------------------
      /* set CUDA streams */
      int sidx = sinfo->idx;
      //-------------------------------------------------------------------
      callCalcGravityFunc(blck, thrd, sinfo, &sidx, laneInfo, pi, 0, tree
#ifndef SERIALIZED_EXECUTION
			  , 0
#endif//SERIALIZED_EXECUTION
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
			  , cycles_dev
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
			  , buf
#ifdef  COUNT_INTERACTIONS
			  , treeInfo
#endif//COUNT_INTERACTIONS
			  );
      //-------------------------------------------------------------------
      /* estimate performance indicator of block time step */
#ifdef  BLOCK_TIME_STEP
      const double block = (double)BLOCKSIZE(BLOCKSIZE(grpNum, NGROUPS), NBLOCKS_PER_SM * devProp.numSM);
      const double share = (double)BLOCKSIZE(BLOCKSIZE(totNum, NGROUPS), NBLOCKS_PER_SM * devProp.numSM);
#ifdef  WALK_TREE_COMBINED_MODEL
      *reduce = share / block;
#else///WALK_TREE_COMBINED_MODEL
      *reduce = block / share;
#endif//WALK_TREE_COMBINED_MODEL
#endif//BLOCK_TIME_STEP
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
      //-------------------------------------------------------------------
      /* gravity from j-particles within other process(es) */
      //-------------------------------------------------------------------
      const int headLETsend = pjNum + ((pjNum % LET_INCNUM_UNIT) != 0) * (LET_INCNUM_UNIT - (pjNum & (LET_INCNUM_UNIT - 1)));
      const int remLETbuf = (int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - headLETsend;
      int idxProcs = 0;
      int remProcs = Nlet - 1;
      while( true ){
	//-----------------------------------------------------------------
	/* get maximum number of processes which possible to communicate by limitation of memory capacity */
	//-----------------------------------------------------------------
	int remBuf = remLETbuf;
	int numProcs = remProcs;
	for(int ii = 0; ii < remProcs; ii++){
	  //---------------------------------------------------------------
	  remBuf -= let[idxProcs + ii].maxSend + let[idxProcs + ii].maxRecv;
	  //---------------------------------------------------------------
	  if( remBuf < 0 ){
	    numProcs = ii - 1;
	    break;
	  }/* if( remBuf < 0 ){ */
	  //---------------------------------------------------------------
	}/* for(int ii = 0; ii < remProcs; ii++){ */
	//-----------------------------------------------------------------
	chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &numProcs, 1, MPI_INT, MPI_MIN, mpi.comm));
	if( numProcs < 1 ){
	  __KILL__(stderr, "ERROR: numProcs is %d, due to lack of remLETbuf(%d) while 0-th target requires numSend(%d) and numRecv(%d).\n\tIncrease EXTEND_NUM_TREE_NODE(%f) defined in src/tree/let.h and/or TREE_SAFETY_VAL(%f) defined in src/tree/make.h.\n", numProcs, remLETbuf, let[idxProcs].maxSend, let[idxProcs].maxRecv, EXTEND_NUM_TREE_NODE, TREE_SAFETY_VAL);
	}/* if( numProcs < 1 ){ */
	//-----------------------------------------------------------------
	let[idxProcs].headSend = headLETsend;
	int numSend = 0;
	int numRecv = 0;
	for(int ii = 0; ii < numProcs - 1; ii++){
	  let[idxProcs + ii    ].headRecv = let[idxProcs + ii].headSend + let[idxProcs + ii].maxSend;	  numSend += let[idxProcs + ii].maxSend;
	  let[idxProcs + ii + 1].headSend = let[idxProcs + ii].headRecv + let[idxProcs + ii].maxRecv;	  numRecv += let[idxProcs + ii].maxRecv;
	}/* for(int ii = 0; ii < numProcs - 1; ii++){ */
	let[idxProcs + numProcs - 1].headRecv = let[idxProcs + numProcs - 1].headSend + let[idxProcs + numProcs - 1].maxSend;
	numSend += let[idxProcs + numProcs - 1].maxSend;
	numRecv += let[idxProcs + numProcs - 1].maxRecv;
	//-----------------------------------------------------------------
#if 1
	if( numSend + numRecv > remLETbuf ){
	  __KILL__(stderr, "ERROR: there is no sufficient size of buffers to complete communications.\n\tnumSend = %d, numRecv = %d, remLETbuf = %d @ rank %d\n", numSend, numRecv, remLETbuf, mpi.rank);
	}
#endif
	//-----------------------------------------------------------------


	//-----------------------------------------------------------------
	/* FUTURE UPDATE: divide below procedure to overlap communication and calculation */
	//-----------------------------------------------------------------
	/* generate numProcs LET(s) and send them to other process(es) */
	for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){
	  //---------------------------------------------------------------
	  const int streamIdxLET = ii % Nstream_let;
	  //---------------------------------------------------------------
#ifndef DBG_PARALLEL_HDF5
	  callGenLET(stream_let[streamIdxLET], let, mpi,
		     tree, let[ii].maxSend, buf
#ifdef  MONITOR_LETGEN_TIME
		     , cycles_let_dev
#endif//MONITOR_LETGEN_TIME
		     );
#ifndef NDEBUG
	  printf("numSend[%d -->> %d] is %d\n", mpi.rank, let[ii].rank, let[ii].numSend);
	  fflush(NULL);
#endif//NDEBUG
#else///DBG_PARALLEL_HDF5
	  let[ii].headSend = 0;
	  let[ii]. numSend = pjNum;
#endif//DBG_PARALLEL_HDF5

#ifdef  LET_COMMUNICATION_VIA_HOST
	  checkCudaErrors(cudaMemcpyAsync(&(tree_hst.more[let[ii].headSend]), &(tree.more[let[ii].headSend]), sizeof(     uint) * let[ii].numSend, cudaMemcpyDeviceToHost, stream_let[streamIdxLET]));
	  checkCudaErrors(cudaMemcpyAsync(&(tree_hst.jpos[let[ii].headSend]), &(tree.jpos[let[ii].headSend]), sizeof(jparticle) * let[ii].numSend, cudaMemcpyDeviceToHost, stream_let[streamIdxLET]));
	  checkCudaErrors(cudaMemcpyAsync(&(tree_hst.mj  [let[ii].headSend]), &(tree.mj  [let[ii].headSend]), sizeof(     real) * let[ii].numSend, cudaMemcpyDeviceToHost, stream_let[streamIdxLET]));
#endif//LET_COMMUNICATION_VIA_HOST
	  /* const int numSendGuess = let[ii].numSend; */
	  /* makeLET(let.let[ii].icom, leafLevel, let.list, let.mask, tree_hst.more, tree_hst.jpos, tree_hst.mj, */
	  /* 	  &(tree_hst.more[let.let[ii].headSend]), &(tree_hst.jpos[let.let[ii].headSend]), &(tree_hst.mj[let.let[ii].headSend]), */
	  /* 	  &(let.let[ii].numSend)); */
	  /* if( let[ii].numSend > numSendGuess ){ */
	  /*   __KILL__(stderr, "ERROR: predicted size of send buffer(%d) is not sufficient for true size of that(%d) @ rank %d with %d-th partner.\n", numSendGuess, let[ii].numSend, mpi.rank, ii); */
	  /* } */
	  //---------------------------------------------------------------
	  /* send # of LET nodes */
	  chkMPIerr(MPI_Isend(&(let[ii].numSend), 1, MPI_INT, let[ii].rank,     mpi.rank, mpi.comm, &(let[ii].reqSendInfo)));
	  let[ii].numRecvGuess = let[ii].maxRecv;
	  chkMPIerr(MPI_Irecv(&(let[ii].numRecv), 1, MPI_INT, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvInfo)));
#ifdef  DBG_LETGEN_ON_GPU
	  fprintf(stdout, "rank = %d: ii = %d: LET(target is %d): numSend = %d out of %d nodes\n", mpi.rank, ii, let[ii].rank, let[ii].numSend, pjNum);
	  fflush(stdout);
#endif//DBG_LETGEN_ON_GPU
	  //---------------------------------------------------------------
	  /* MPI_Isend */
#ifdef  LET_COMMUNICATION_VIA_HOST
	  checkCudaErrors(cudaStreamSynchronize(stream_let[streamIdxLET]));
	  chkMPIerr(MPI_Isend(&(tree_hst.more[let[ii].headSend]), let[ii].numSend, mpi.more, let[ii].rank, mpi.rank, mpi.comm, &(let[ii].reqSendMore)));
	  chkMPIerr(MPI_Isend(&(tree_hst.jpos[let[ii].headSend]), let[ii].numSend, mpi.jpos, let[ii].rank, mpi.rank, mpi.comm, &(let[ii].reqSendJpos)));
	  chkMPIerr(MPI_Isend(&(tree_hst.mj  [let[ii].headSend]), let[ii].numSend, mpi.mass, let[ii].rank, mpi.rank, mpi.comm, &(let[ii].reqSendMass)));
#else///LET_COMMUNICATION_VIA_HOST
	  chkMPIerr(MPI_Isend(&(tree.more[let[ii].headSend]), let[ii].numSend, mpi.more, let[ii].rank, mpi.rank, mpi.comm, &(let[ii].reqSendMore)));
	  chkMPIerr(MPI_Isend(&(tree.jpos[let[ii].headSend]), let[ii].numSend, mpi.jpos, let[ii].rank, mpi.rank, mpi.comm, &(let[ii].reqSendJpos)));
	  chkMPIerr(MPI_Isend(&(tree.mj  [let[ii].headSend]), let[ii].numSend, mpi.mass, let[ii].rank, mpi.rank, mpi.comm, &(let[ii].reqSendMass)));
#endif//LET_COMMUNICATION_VIA_HOST
	  //---------------------------------------------------------------
	}
	//-----------------------------------------------------------------
	for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){
	  //---------------------------------------------------------------
	  MPI_Status status;
	  chkMPIerr(MPI_Wait(&(let[ii].reqRecvInfo), &status));
#ifdef  DBG_LETGEN_ON_GPU
	  fprintf(stdout, "rank = %d: ii = %d: LET(origin is %d): numRecv = %d\n",
		  mpi.rank, ii, let[ii].rank, let[ii].numRecv);
	  fflush(stdout);
#endif//DBG_LETGEN_ON_GPU
	  if( let[ii].numRecv > let[ii].numRecvGuess ){
	    __KILL__(stderr, "ERROR: predicted size of receive buffer(%d) is not sufficient for true size of that(%d) @ rank %d with %d-th partner.\n\tsuggestion: consider increasing \"LETSIZE_REDUCE_FACTOR\" defined in src/tree/let.h (current value is %f) to at least %f.\n", let[ii].numRecvGuess, let[ii].numRecv, mpi.rank, ii, LETSIZE_REDUCE_FACTOR, LETSIZE_REDUCE_FACTOR * (float)let[ii].numRecv / (float)let[ii].numRecvGuess);
	  }/* if( let[ii].numRecv > let[ii].numRecvGuess ){ */
	  //---------------------------------------------------------------
	  /* receive # of LET nodes */
#ifdef  LET_COMMUNICATION_VIA_HOST
	  chkMPIerr(MPI_Irecv(&(tree_hst.more[let[ii].headRecv]), let[ii].numRecv, mpi.more, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvMore)));
	  chkMPIerr(MPI_Irecv(&(tree_hst.jpos[let[ii].headRecv]), let[ii].numRecv, mpi.jpos, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvJpos)));
	  chkMPIerr(MPI_Irecv(&(tree_hst.mj  [let[ii].headRecv]), let[ii].numRecv, mpi.mass, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvMass)));
#else///LET_COMMUNICATION_VIA_HOST
	  chkMPIerr(MPI_Irecv(&(tree.more[let[ii].headRecv]), let[ii].numRecv, mpi.more, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvMore)));
	  chkMPIerr(MPI_Irecv(&(tree.jpos[let[ii].headRecv]), let[ii].numRecv, mpi.jpos, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvJpos)));
	  chkMPIerr(MPI_Irecv(&(tree.mj  [let[ii].headRecv]), let[ii].numRecv, mpi.mass, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvMass)));
#endif//LET_COMMUNICATION_VIA_HOST
	  //---------------------------------------------------------------
	}
	//-----------------------------------------------------------------
	/* receive numProcs LET(s) and calculate gravity from them */
	for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){
	  //---------------------------------------------------------------
	  /* copy LET nodes from host to device */
	  //---------------------------------------------------------------
	  const int streamIdxLET = ii % Nstream_let;
	  MPI_Status statusMore;	  chkMPIerr(MPI_Wait(&(let[ii].reqRecvJpos), &statusMore));
#ifdef  LET_COMMUNICATION_VIA_HOST
	  checkCudaErrors(cudaMemcpyAsync(&(tree.jpos[let[ii].headRecv]), &(tree_hst.jpos[let[ii].headRecv]), sizeof(jparticle) * let[ii].numRecv, cudaMemcpyHostToDevice, stream_let[streamIdxLET]));
#endif//LET_COMMUNICATION_VIA_HOST
	  MPI_Status statusMass;	  chkMPIerr(MPI_Wait(&(let[ii].reqRecvMass), &statusMass));
#ifdef  LET_COMMUNICATION_VIA_HOST
	  checkCudaErrors(cudaMemcpyAsync(&(tree.mj  [let[ii].headRecv]), &(tree_hst.mj  [let[ii].headRecv]), sizeof(     real) * let[ii].numRecv, cudaMemcpyHostToDevice, stream_let[streamIdxLET]));
#endif//LET_COMMUNICATION_VIA_HOST
	  MPI_Status statusJpos;	  chkMPIerr(MPI_Wait(&(let[ii].reqRecvMore), &statusJpos));
#ifdef  LET_COMMUNICATION_VIA_HOST
	  checkCudaErrors(cudaMemcpyAsync(&(tree.more[let[ii].headRecv]), &(tree_hst.more[let[ii].headRecv]), sizeof(     uint) * let[ii].numRecv, cudaMemcpyHostToDevice, stream_let[streamIdxLET]));
#endif//LET_COMMUNICATION_VIA_HOST
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaStreamSynchronize(stream_let[streamIdxLET]));
#ifdef  DBG_LETGEN_ON_GPU
	  fprintf(stdout, "received LET from rank %d (%d-th partner)\n", let[ii].rank, ii);
	  fflush(stdout);
	  /* MPI_Finalize(); */
	  /* exit(0); */
#endif//DBG_LETGEN_ON_GPU
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* calculate gravity from LET */
	  //---------------------------------------------------------------
	  callCalcGravityFunc(blck, thrd, sinfo, &sidx, laneInfo, pi, 0, tree, let[ii].headRecv, cycles_dev, buf
#ifdef  COUNT_INTERACTIONS
			      , treeInfo
#endif//COUNT_INTERACTIONS
			      );
	  //---------------------------------------------------------------
	}/* for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){ */
	//-----------------------------------------------------------------
	for(int ii = 0; ii < numProcs; ii++){
	  //---------------------------------------------------------------
	  MPI_Status statusInfo;	  chkMPIerr(MPI_Wait(&(let[ii].reqSendInfo), &statusInfo));
	  MPI_Status statusMore;	  chkMPIerr(MPI_Wait(&(let[ii].reqSendMore), &statusMore));
	  MPI_Status statusJpos;	  chkMPIerr(MPI_Wait(&(let[ii].reqSendJpos), &statusJpos));
	  MPI_Status statusMass;	  chkMPIerr(MPI_Wait(&(let[ii].reqSendMass), &statusMass));
	  //---------------------------------------------------------------
	}/* for(int ii = 0; ii < numProcs; ii++){ */
	//-----------------------------------------------------------------
	idxProcs += numProcs;
	remProcs -= numProcs;
	//-----------------------------------------------------------------


	//-----------------------------------------------------------------
	if( remProcs <= 0 )	  break;
	//-----------------------------------------------------------------
      }/* while( true ){ */
      //-------------------------------------------------------------------
      /* preparation for communication in the next step */
#ifdef  BLOCK_TIME_STEP
      const float letsize_scaler = (float)(share / block);
#else///BLOCK_TIME_STEP
      const float letsize_scaler = UNITY;
#endif//BLOCK_TIME_STEP
      for(int ii = 0; ii < Nlet - 1; ii++){
	//-----------------------------------------------------------------
	/* /\* guess the minimum size of the buffer *\/ */
	/* int minSend = (int)ceilf(letsize_scaler * (float)let[ii].numSend);	minSend += 32 - (minSend & 31); */
	/* int minRecv = (int)ceilf(letsize_scaler * (float)let[ii].numRecv);	minRecv += 32 - (minRecv & 31); */
	//-----------------------------------------------------------------
	if( ceilf(letsize_scaler * (float)let[ii].numSend) < (LETSIZE_REDUCE_CRITERION * (float)let[ii].maxSend) )	  let[ii].overEstimateSend++;
	if( ceilf(letsize_scaler * (float)let[ii].numRecv) < (LETSIZE_REDUCE_CRITERION * (float)let[ii].maxRecv) )	  let[ii].overEstimateRecv++;
	//-----------------------------------------------------------------
	if( let[ii].overEstimateSend >= LETSIZE_OVERESTIMATION_STEPS ){
	  let[ii].maxSend = (int)ceilf(LETSIZE_REDUCE_FACTOR * (float)let[ii].maxSend);	  let[ii].maxSend += 32 - (let[ii].maxSend & 31);
	  let[ii].overEstimateSend = 0;
	}/* if( let[ii].overEstimateSend >= LETSIZE_OVERESTIMATION_STEPS ){ */
	//-----------------------------------------------------------------
	if( let[ii].overEstimateRecv >= LETSIZE_OVERESTIMATION_STEPS ){
	  let[ii].maxRecv = (int)ceilf(LETSIZE_REDUCE_FACTOR * (float)let[ii].maxRecv);	  let[ii].maxRecv += 32 - (let[ii].maxRecv & 31);
	  let[ii].overEstimateRecv = 0;
	}/* if( let[ii].overEstimateRecv >= LETSIZE_OVERESTIMATION_STEPS ){ */
	//-----------------------------------------------------------------
	/* let[ii].maxSend = (int)ceilf(fminf((float)let[ii].maxSend, letsize_scaler * (float)let[ii].numSend));	let[ii].maxSend += 32 - (let[ii].maxSend & 31); */
	/* let[ii].maxRecv = (int)ceilf(fminf((float)let[ii].maxRecv, letsize_scaler * (float)let[ii].numRecv));	let[ii].maxRecv += 32 - (let[ii].maxRecv & 31); */
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < Nlet - 1; ii++){ */
      setLETpartition(Nlet, let);
#if 0
      fprintf(stderr, "maxSend = %d while pjNum = %d @ rank %d\n", let[0].maxSend, pjNum, mpi.rank);
      fflush(stderr);
#endif
      //-------------------------------------------------------------------
      /* /\* check the size does not exceed the LET buffer *\/ */
      /* if( let[0].maxSend + let[0].maxRecv > ((int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - headLETsend) ){ */
      /* 	__KILL__(stderr, "ERROR: rough estimation routine predicts # of LET nodes(%d) succeed the capacity of the array(%d - %d); communication with first partner would fail due to lack of MPI buffer.\n", let[0].maxSend + let[0].maxRecv, (int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE), headLETsend); */
      /* }/\* if( let[0].maxSend + let[0].maxRecv > ((int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - headLETsend) ){ *\/ */
      //-------------------------------------------------------------------
#endif//SERIALIZED_EXECUTION
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
      sinfo->idx = sidx;
      checkCudaErrors(cudaDeviceSynchronize());
      /* checkCudaErrors(cudaStreamSynchronize(sinfo->stream[sidx    ])); */
      /* checkCudaErrors(cudaStreamSynchronize(sinfo->stream[sidx ^ 1])); */
      //-------------------------------------------------------------------
#ifdef  DBG_LETGEN_ON_GPU
      fprintf(stdout, "force calculation finished on rank %d\n", mpi.rank);
      fflush(stdout);
#endif//DBG_LETGEN_ON_GPU
      //-------------------------------------------------------------------
      int fail_hst;
      checkCudaErrors(cudaMemcpy(&fail_hst, buf.fail, sizeof(int), cudaMemcpyDeviceToHost));
      if( fail_hst != 0 ){
#ifdef  SERIALIZED_EXECUTION
	__KILL__(stderr, "ERROR: bufUsed exceeds bufSize of %d at least %u times.\nPLEASE re-simulate after decreasing NUM_BODY_MAX(%d) or GLOBAL_MEMORY_SYSBUF(%zu) defined in src/misc/structure.h or TREE_SAFETY_VAL(%f) defined in src/tree/make.h.\n", buf.bufSize, fail_hst, NUM_BODY_MAX, (size_t)GLOBAL_MEMORY_SYSBUF, TREE_SAFETY_VAL);
#else///SERIALIZED_EXECUTION
	__KILL__(stderr, "ERROR: bufUsed exceeds bufSize of %d at least %u times.\nPLEASE re-simulate after decreasing NUM_BODY_MAX(%d) or GLOBAL_MEMORY_SYSBUF(%zu) defined in src/misc/structure.h or TREE_SAFETY_VAL(%f) defined in src/tree/make.h, or EXTEND_NUM_TREE_NODE(%f) defined in src/tree/let.h.\n", buf.bufSize, fail_hst, NUM_BODY_MAX, (size_t)GLOBAL_MEMORY_SYSBUF, TREE_SAFETY_VAL, EXTEND_NUM_TREE_NODE);
#endif//SERIALIZED_EXECUTION
      }/* if( fail_hst != 0 ){ */
      //-------------------------------------------------------------------
    }
  //-----------------------------------------------------------------------
#ifdef  COMPARE_WITH_DIRECT_SOLVER
  else{
    //---------------------------------------------------------------------
    Nrem = BLOCKSIZE(Ni, NTHREADS);
    if( Nrem <= MAX_BLOCKS_PER_GRID )
      calcAccDirect_kernel<<<Nrem, NTHREADS>>>
	(pi.pos, pi.acc, pi.pos, Ni
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	 , eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
	 );
    //---------------------------------------------------------------------
    else{
      //-------------------------------------------------------------------
      const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
      int hidx = 0;
      //-------------------------------------------------------------------
      for(int iter = 0; iter < Niter; iter++){
	//-----------------------------------------------------------------
	int Nblck = MAX_BLOCKS_PER_GRID;
	if( Nblck > Nrem )	  Nblck = Nrem;
	//-----------------------------------------------------------------
	int Nsub = Nblck * NTHREADS;
	calcAccDirect_kernel<<<Nblck, NTHREADS>>>
	  (&pi.pos[hidx], &pi.acc[hidx], &pi.pos[hidx], Nsub
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	   , &eps2[hidx]
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
	   );
	//-----------------------------------------------------------------
	hidx += Nsub;
	Nrem -= Nblck;
	//-----------------------------------------------------------------
      }/* for(int iter = 0; iter < Niter; iter++){ */
      //-------------------------------------------------------------------
    }/* else{ */
    //---------------------------------------------------------------------
/*     calcAccDirect_kernel<<<BLOCKSIZE(Ni, NTHREADS), NTHREADS>>> */
/*       (pi.pos, pi.acc, pi.pos, Ni */
/* #ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING */
/*        , eps2 */
/* #endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING */
/*        ); */
    //---------------------------------------------------------------------
    getLastCudaError("calcAccDirect");
    //---------------------------------------------------------------------
  }
#endif//COMPARE_WITH_DIRECT_SOLVER
  //-----------------------------------------------------------------------
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
  checkCudaErrors(cudaDeviceSynchronize());
  /* get total clock cycles to compute enclosing ball */
  checkCudaErrors(cudaMemcpy(cycles_hst, cycles_dev, sizeof(unsigned long long int), cudaMemcpyDeviceToHost));
  /* get information on enclosing ball */
  acceleration *seb;
  const int Nseb = blck.x * NGROUPS;
  mycudaMallocHost((void **)&seb, (size_t)Nseb * sizeof(acceleration));
  checkCudaErrors(cudaMemcpy(seb, pi.acc, (size_t)Nseb * sizeof(acceleration), cudaMemcpyDeviceToHost));
  /* set file tag */
  FILE *fp;
#ifndef ADOPT_ENCLOSING_BALL
  char sebfile[] = "com";
#else///ADOPT_ENCLOSING_BALL
#ifdef  ADOPT_SMALLEST_ENCLOSING_BALL
  char sebfile[] = "fischer03";
#endif//ADOPT_SMALLEST_ENCLOSING_BALL
#ifdef  ADOPT_APPROXIMATED_ENCLOSING_BALL
  char sebfile[] = "ritter90";
#endif//ADOPT_APPROXIMATED_ENCLOSING_BALL
#ifdef  COMPARE_ENCLOSING_BALLS
  char sebfile[] = "smaller";
#endif//COMPARE_ENCLOSING_BALLS
#   if  !defined(ADOPT_SMALLEST_ENCLOSING_BALL) && !defined(ADOPT_APPROXIMATED_ENCLOSING_BALL) && !defined(COMPARE_ENCLOSING_BALLS)
  char sebfile[] = "cartesian";
#endif//!defined(ADOPT_SMALLEST_ENCLOSING_BALL) && !defined(ADOPT_APPROXIMATED_ENCLOSING_BALL) && !defined(COMPARE_ENCLOSING_BALLS)
#endif//ADOPT_ENCLOSING_BALL
  char filename[128];
  /* output computing cost of enclosing ball */
  sprintf(filename, "%s/%s.ball.clock.%s.txt", LOGFOLDER, file, sebfile);
  fp = fopen(filename, "a");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  char date[64], hostname[64];
  getPresentDateInStrings(date);
  gethostname(hostname, sizeof(hostname));
  fprintf(fp, "\twith %s equipped on %s\n", devProp.name, devProp.host);
  fprintf(fp, "\tmeasured on %s", date);
  fprintf(fp, "Nseb = %d, Ntot = %d, enclosing ball of continuous %d particles\n", Nseb, Ni, DIV_NWARP(TSUB));
  fprintf(fp, "%Lu cycles @ Ttot = %d, Tsub = %d, Nb_sm = %d, Nsm = %d\n", *cycles_hst, NTHREADS, TSUB, NBLOCKS_PER_SM, devProp.numSM);
  *cycles_hst /= (unsigned long long int)((NTHREADS >> 5) * (devProp.numSM * NBLOCKS_PER_SM));/* divide by product of (# of warps within a thread) and (# of concurrent blocks) */
  fprintf(fp, "%Lu cycles after divided by the product of (# of warps within a thread) and (# of concurrent blocks)\n", *cycles_hst);
  fprintf(fp, "%le seconds @ %lf GHz\n", (double)(*cycles_hst) / (devProp.coreClk * 1.0e+9), devProp.coreClk);
  fclose(fp);
  /* output properties of enclosing ball */
  sprintf(filename, "%s/%s.ball.%s.dat", DATAFOLDER, file, sebfile);
  fp = fopen(filename, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  fwrite(seb, sizeof(acceleration), Nseb, fp);
  fclose(fp);
  /* output summary of enclosing ball for future analysis */
  sprintf(date, "%s/%s.ball.info.txt", LOGFOLDER, file);
  fp = fopen(date, "a");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", date);  }
  fprintf(fp, "%s\t%d\n", filename, Nseb);
  fclose(fp);
  /* finalize the computation */
  mycudaFreeHost(seb);
  exit(0);
#endif//PRINT_PSEUDO_PARTICLE_INFO
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* mutiply Gravitational constant and subtract self-interaction for potential */
  //-----------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
  Nrem = BLOCKSIZE(grpNum, NWARP * NGROUPS);
#else///BLOCK_TIME_STEP
  Nrem = BLOCKSIZE(Ni, NTHREADS);
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
  /* when grid splitting is not required... */
  if( Nrem <= MAX_BLOCKS_PER_GRID ){
#ifdef  BLOCK_TIME_STEP
#ifndef SERIALIZED_EXECUTION
    if( grpNum != 0 )
#endif//SERIALIZED_EXECUTION
      trimAcc_kernel<<<Nrem, thrd>>>(pi.acc, pi.pos, BLOCKSIZE(grpNum, NGROUPS) * NGROUPS, laneInfo);
#else///BLOCK_TIME_STEP
    trimAcc_kernel<<<Nrem, NTHREADS>>>(pi.acc, pi.pos);
#endif//BLOCK_TIME_STEP
  }/* if( Nrem <= MAX_BLOCKS_PER_GRID ){ */
  //-----------------------------------------------------------------------
  /* when grid splitting is required... */
  else{
    //---------------------------------------------------------------------
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;
    //---------------------------------------------------------------------
    for(int iter = 0; iter < Niter; iter++){
      //-------------------------------------------------------------------
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;
      //-------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
      int Nsub = Nblck * NWARP * NGROUPS;
      trimAcc_kernel<<<Nblck, thrd.x>>>(pi.acc, pi.pos, BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx]);
#else///BLOCK_TIME_STEP
      int Nsub = Nblck * NTHREADS;
      trimAcc_kernel<<<Nblck, NTHREADS>>>(&pi.acc[hidx], &pi.pos[hidx]);
#endif//BLOCK_TIME_STEP
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
/* #ifdef  BLOCK_TIME_STEP */
/* #ifndef SERIALIZED_EXECUTION */
/*   if( grpNum != 0 ) */
/* #endif//SERIALIZED_EXECUTION */
/*     { */
/*       /\* trimAcc_kernel<<<                              blck, thrd>>>(pi.acc, pi.pos, laneInfo); *\/ */
/*       trimAcc_kernel<<<BLOCKSIZE(grpNum, NWARP * NGROUPS), thrd>>>(pi.acc, pi.pos, BLOCKSIZE(grpNum, NGROUPS) * NGROUPS, laneInfo); */
/* #endif */
/*     } */
/* #else///BLOCK_TIME_STEP */
/*   trimAcc_kernel<<<BLOCKSIZE(Ni, NTHREADS), NTHREADS>>>(pi.acc, pi.pos); */
/* #endif//BLOCK_TIME_STEP */
  //-----------------------------------------------------------------------
  getLastCudaError("trimAcc_kernel");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  static struct timeval finish;
  checkCudaErrors(cudaDeviceSynchronize());
  gettimeofday(&finish, NULL);
  *time = calcElapsedTimeInSec(start, finish);
  //-----------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
  /* evaluate GPU time based on clock cycle counter */
  checkCudaErrors(cudaMemcpy(cycles_hst, cycles_dev, sizeof(unsigned long long int), cudaMemcpyDeviceToHost));
  *twalk = (double)(*cycles_hst) * devProp.coreClk / (double)(devProp.numSM * NBLOCKS_PER_SM);
#ifdef  MONITOR_LETGEN_TIME
  checkCudaErrors(cudaMemcpy(cycles_let_hst, cycles_let_dev, sizeof(unsigned long long int), cudaMemcpyDeviceToHost));
  *tlet = (double)(*cycles_let_hst) * devProp.coreClk / (double)(devProp.numSM * NBLOCKS_PER_SM);
#endif//MONITOR_LETGEN_TIME
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  elapsed->calcGravity_dev += *time;
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
extern "C"
void setGlobalConstants_walk_dev_cu
(const real newton_hst, const real eps2_hst
#ifndef WS93_MAC
 , const real theta2_hst
#endif//WS93_MAC
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const real epsinv_hst = RSQRT(eps2_hst);
  //-----------------------------------------------------------------------
  jnode jnode0_hst;
#pragma unroll
  for(int ii = 0; ii < NSTOCK; ii++)
    jnode0_hst.idx[ii] = 0;
  //-----------------------------------------------------------------------
#   if  CUDART_VERSION >= 5000
  cudaMemcpyToSymbol( newton , &newton_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#ifndef INDIVIDUAL_GRAVITATIONAL_SOFTENING
  cudaMemcpyToSymbol( eps2   , &  eps2_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
  cudaMemcpyToSymbol( epsinv , &epsinv_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#ifndef WS93_MAC
  cudaMemcpyToSymbol( theta2 , &theta2_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#endif//WS93_MAC
  cudaMemcpyToSymbol( jnode0 , &jnode0_hst, sizeof(jnode), 0, cudaMemcpyHostToDevice);
#else//CUDART_VERSION >= 5000
  cudaMemcpyToSymbol("newton", &newton_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#ifndef INDIVIDUAL_GRAVITATIONAL_SOFTENING
  cudaMemcpyToSymbol("eps2"  , &  eps2_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
  cudaMemcpyToSymbol("epsinv", &epsinv_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#ifndef WS93_MAC
  cudaMemcpyToSymbol("theta2", &theta2_hst, sizeof( real), 0, cudaMemcpyHostToDevice);
#endif//WS93_MAC
  cudaMemcpyToSymbol("jnode0", &jnode0_hst, sizeof(jnode), 0, cudaMemcpyHostToDevice);
#endif//CUDART_VERSION >= 5000
  //-----------------------------------------------------------------------
#   if  SMPREF == 1
  checkCudaErrors(cudaFuncSetCacheConfig(calcAcc_kernel, cudaFuncCachePreferShared));
#endif//SMPREF == 1
#   if  WIDEBANK == 1
  checkCudaErrors(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte));
#endif//WIDEBANK == 1
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* error checking before running the kernel */
  //-----------------------------------------------------------------------
  struct cudaFuncAttributes funcAttr;
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, calcAcc_kernel));
  int regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS);
  int memLimit = ((SMPREF == 1) ? (48 * 1024) : (16 * 1024)) / funcAttr.sharedSizeBytes;
  int Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck != NBLOCKS_PER_SM ){
    //---------------------------------------------------------------------
    __KILL__(stderr, "ERROR: # of blocks per SM for calcAcc_kernel is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM defined in src/tree/walk_dev.cu is %d\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM);
    //---------------------------------------------------------------------
  }/* if( Nblck != NBLOCKS_PER_SM ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
