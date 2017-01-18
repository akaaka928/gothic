/*************************************************************************\
 *                                                                       *
                  last updated on 2017/01/17(Tue) 20:11:35
 *                                                                       *
 *    Octree N-body calculation for collisionless systems on NVIDIA GPUs *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
/* #define DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
#define DOUBLE_BUFFER_FOR_LET
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <helper_cuda.h>
#include <sys/time.h>
#ifndef SERIALIZED_EXECUTION
#include <mpi.h>
#endif//SERIALIZED_EXECUTION
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
#include <unistd.h>
#endif//PRINT_PSEUDO_PARTICLE_INFO
//-------------------------------------------------------------------------
#include "macro.h"
#include "cudalib.h"
#include "timer.h"
#ifndef SERIALIZED_EXECUTION
#include "mpilib.h"
#endif//SERIALIZED_EXECUTION
#ifdef  PRINT_PSEUDO_PARTICLE_INFO
#include "name.h"
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
#include "../misc/tune.h"
#include "../para/mpicfg.h"
#include "let.h"
#include "let_dev.h"
#endif//SERIALIZED_EXECUTION
//-------------------------------------------------------------------------
#include "walk_dev.h"
//-------------------------------------------------------------------------
#   if  !defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
#define USE_GPU_BASE_CLOCK_FREQ
#if 1
#   if  (__CUDACC_VER_MINOR__ + 10 * __CUDACC_VER_MAJOR__) >= 80
#include <nvml.h>
#undef  USE_GPU_BASE_CLOCK_FREQ
#define USE_MEASURED_CLOCK_FREQ
nvmlDevice_t deviceHandler;
#endif//(__CUDACC_VER_MINOR__ + 10 * __CUDACC_VER_MAJOR__) >= 80
#endif
#endif//!defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
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
muse setCUDAstreams_dev(cudaStream_t **stream, kernelStream *sinfo, deviceInfo *info, deviceProp *prop
/* #   if  defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)) */
/* 			, cudaEvent_t **iniEvent, cudaEvent_t **finEvent */
/* #endif//defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)) */
			)
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
  *stream = (cudaStream_t *)malloc((size_t)(sinfo->num) * sizeof(cudaStream_t));  if( *stream == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stream\n");  }
  alloc.host +=                    (size_t)(sinfo->num) * sizeof(cudaStream_t) ;
  sinfo->stream = *stream;
  //-----------------------------------------------------------------------
  /* set CUDA streams */
  for(int ii = 0; ii < 2; ii++)
    sinfo->stream[ii] = info->stream[ii];
  for(int ii = 2; ii < sinfo->num; ii++)
    checkCudaErrors(cudaStreamCreate(&(sinfo->stream[ii])));
  //-----------------------------------------------------------------------
#if 0
  int priority;
  for(int ii = 0; ii < sinfo->num; ii++){
    checkCudaErrors(cudaStreamGetPriority(sinfo->stream[ii], &priority));
    fprintf(stdout, "priority of stream[%d] is %d\n", ii, priority);
  }/* for(int ii = 0; ii < sinfo->num; ii++){ */
  fflush(stdout);
#endif
  //-----------------------------------------------------------------------
/* #if 1 */
/*   int highest, lowest; */
/*   checkCudaErrors(cudaDeviceGetStreamPriorityRange(&lowest, &highest)); */
/*   for(int ii = 0; ii < *Nstream; ii++) */
/*     checkCudaErrors(cudaStreamCreateWithPriority(&((*stream)[ii]), cudaStreamDefault, highest)); */
/*     /\* checkCudaErrors(cudaStreamCreateWithPriority(&((*stream)[ii]), cudaStreamNonBlocking, highest)); *\/ */
/* #else */
/* #pragma unroll */
/*   for(int ii = 0; ii < *Nstream; ii++) */
/*     checkCudaErrors(cudaStreamCreate(&((*stream)[ii]))); */
/* #endif */
  //-----------------------------------------------------------------------
/*   /\* allocate and set CUDA events *\/ */
/* #   if  defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)) */
/*   *iniEvent = (cudaEvent_t *)malloc((size_t)(sinfo->num) * sizeof(cudaEvent_t));  if( *iniEvent == NULL ){    __KILL__(stderr, "ERROR: failure to allocate iniEvent\n");  } */
/*   *finEvent = (cudaEvent_t *)malloc((size_t)(sinfo->num) * sizeof(cudaEvent_t));  if( *finEvent == NULL ){    __KILL__(stderr, "ERROR: failure to allocate finEvent\n");  } */
/*   alloc.host +=                     (size_t)(sinfo->num) * sizeof(cudaEvent_t); */
/*   alloc.host +=                     (size_t)(sinfo->num) * sizeof(cudaEvent_t); */
/*   for(int ii = 0; ii < sinfo->num; ii++){ */
/*     checkCudaErrors(cudaEventCreate(&((*iniEvent)[ii]))); */
/*     checkCudaErrors(cudaEventCreate(&((*finEvent)[ii]))); */
/*   }/\* for(int ii = 0; ii < sinfo->num; ii++){ *\/ */
/* #endif//defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)) */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* set CUDA streams */
//-------------------------------------------------------------------------
#   if  defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
//-------------------------------------------------------------------------
extern "C"
muse allocateCUDAevents_dev
(cudaEvent_t **iniWalk, cudaEvent_t **finWalk
#ifdef  MONITOR_LETGEN_TIME
 , cudaEvent_t **iniMake, cudaEvent_t **finMake
#endif//MONITOR_LETGEN_TIME
 , const int Ngpu)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* allocate array for CUDA events */
  *iniWalk = (cudaEvent_t *)malloc((size_t)Ngpu * sizeof(cudaEvent_t));  if( *iniWalk == NULL ){    __KILL__(stderr, "ERROR: failure to allocate iniWalk\n");  }
  *finWalk = (cudaEvent_t *)malloc((size_t)Ngpu * sizeof(cudaEvent_t));  if( *finWalk == NULL ){    __KILL__(stderr, "ERROR: failure to allocate finWalk\n");  }
  alloc.host +=                    (size_t)Ngpu * sizeof(cudaEvent_t) ;
  alloc.host +=                    (size_t)Ngpu * sizeof(cudaEvent_t) ;
#ifdef  MONITOR_LETGEN_TIME
  *iniMake = (cudaEvent_t *)malloc((size_t)(Ngpu - 1) * sizeof(cudaEvent_t));	 if( *iniMake == NULL ){    __KILL__(stderr, "ERROR: failure to allocate iniMake\n");  }
  *finMake = (cudaEvent_t *)malloc((size_t)(Ngpu - 1) * sizeof(cudaEvent_t));	 if( *finMake == NULL ){    __KILL__(stderr, "ERROR: failure to allocate finMake\n");  }
  alloc.host +=                    (size_t)(Ngpu - 1) * sizeof(cudaEvent_t) ;
  alloc.host +=                    (size_t)(Ngpu - 1) * sizeof(cudaEvent_t) ;
#endif//MONITOR_LETGEN_TIME
  //-----------------------------------------------------------------------
  /* set CUDA events */
  for(int ii = 0; ii < Ngpu; ii++){
    checkCudaErrors(cudaEventCreate(&((*iniWalk)[ii])));
    checkCudaErrors(cudaEventCreate(&((*finWalk)[ii])));
  }/* for(int ii = 0; ii < Ngpu; ii++){ */
#ifdef  MONITOR_LETGEN_TIME
  for(int ii = 0; ii < Ngpu - 1; ii++){
    checkCudaErrors(cudaEventCreate(&((*iniMake)[ii])));
    checkCudaErrors(cudaEventCreate(&((*finMake)[ii])));
  }/* for(int ii = 0; ii < Ngpu; ii++){ */
#endif//MONITOR_LETGEN_TIME
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
void  releaseCUDAevents_dev
(cudaEvent_t  *iniWalk, cudaEvent_t  *finWalk
#ifdef  MONITOR_LETGEN_TIME
 , cudaEvent_t  *iniMake, cudaEvent_t  *finMake
#endif//MONITOR_LETGEN_TIME
 , const int Ngpu)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* destroy CUDA events */
  for(int ii = 0; ii < Ngpu; ii++){
    mycudaEventDestroy(iniWalk[ii]);
    mycudaEventDestroy(finWalk[ii]);
  }/* for(int ii = 0; ii < Ngpu; ii++){ */
#ifdef  MONITOR_LETGEN_TIME
  for(int ii = 0; ii < Ngpu - 1; ii++){
    mycudaEventDestroy(iniMake[ii]);
    mycudaEventDestroy(finMake[ii]);
  }/* for(int ii = 0; ii < Ngpu - 1; ii++){ */
#endif//MONITOR_LETGEN_TIME
  //-----------------------------------------------------------------------
  /* deallocate CUDA events */
  free(iniWalk);
  free(finWalk);
#ifdef  MONITOR_LETGEN_TIME
  free(iniMake);
  free(finMake);
#endif//MONITOR_LETGEN_TIME
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
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
#ifndef USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
 , unsigned long long int  *cycles_hst, unsigned long long int  *cycles_dev
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
 , unsigned long long int  *cycles_let_hst, unsigned long long int  *cycles_let_dev
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#endif//USE_CUDA_EVENT
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
#ifndef USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  mycudaFree    (cycles_dev);
  mycudaFreeHost(cycles_hst);
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  mycudaFree    (cycles_let_dev);
  mycudaFreeHost(cycles_let_hst);
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#endif//USE_CUDA_EVENT
  //-----------------------------------------------------------------------
#ifdef  USE_MEASURED_CLOCK_FREQ
  nvmlShutdown();
#endif//USE_MEASURED_CLOCK_FREQ
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
#ifndef USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
 unsigned long long int **cycles_hst, unsigned long long int **cycles_dev,
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
 unsigned long long int **cycles_let_hst, unsigned long long int **cycles_let_dev,
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#endif//USE_CUDA_EVENT
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
#ifndef USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  mycudaMalloc    ((void **)cycles_dev, sizeof(unsigned long long int));  alloc.device += sizeof(unsigned long long int);
  mycudaMallocHost((void **)cycles_hst, sizeof(unsigned long long int));  alloc.host   += sizeof(unsigned long long int);
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  //-----------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  mycudaMalloc    ((void **)cycles_let_dev, sizeof(unsigned long long int));  alloc.device += sizeof(unsigned long long int);
  mycudaMallocHost((void **)cycles_let_hst, sizeof(unsigned long long int));  alloc.host   += sizeof(unsigned long long int);
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#endif//USE_CUDA_EVENT
  //-----------------------------------------------------------------------
#ifdef  USE_MEASURED_CLOCK_FREQ
  nvmlInit();
#if 1
  nvmlDeviceGetHandleByIndex(gpu.idx, &deviceHandler);
#else
  nvmlReturn_t nvmlMsg = nvmlDeviceGetHandleByIndex(gpu.idx, &deviceHandler);
  printf("nvmlMsg = %d\n", nvmlMsg);
  MPI_Finalize();
  exit(0);
#endif
#endif//USE_MEASURED_CLOCK_FREQ
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
  if( (booked / ((size_t)(NGROUPS * nblocks) * (sizeof(uint)))) > INT_MAX ){
    __KILL__(stderr, "ERROR: expected size for bufUnit (%zu) exceeds INT_MAX\n\trewrite \"calcAcc_kernel()\" in \"src/tree/walk_dev.cu\"\n", (booked / ((size_t)(NGROUPS * nblocks) * (sizeof(uint)))));
  }/* if( (booked / ((size_t)(NGROUPS * nblocks) * (sizeof(uint)))) > INT_MAX ){ */
  int bufUnit = (int)(booked / ((size_t)(NGROUPS * nblocks) * (sizeof(uint))));
  /* *bufUnit should be aligned in 32 bytes order (= 128 bits) --> 8 or 4 elements for single or double precision, respectively */
  bufUnit -= (bufUnit & 7);
#ifndef SERIALIZED_EXECUTION
  if( ((size_t)bufUnit * (size_t)NGROUPS) > INT_MAX ){
    __KILL__(stderr, "ERROR: expected size for bufUnit for LET (%zu) exceeds INT_MAX\n\trewrite \"makeLET_kernel()\" in \"src/tree/let_dev.cu\"\n", ((size_t)bufUnit * (size_t)NGROUPS));
  }/* if( ((size_t)bufUnit * (size_t)NGROUPS) > INT_MAX ){ */
#endif//SERIALIZED_EXECUTION
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
  buf->freeNum = NBLOCKS_PER_SM * gpu.numSM;
#endif//!defined(USE_SMID_TO_GET_BUFID) &&  defined(TRY_MODE_ABOUT_BUFFER)
  buf->bufSize = bufUnit;
  //-----------------------------------------------------------------------
#if 0
  fprintf(stdout, "bufUnit = %d, bufTot = %zu, bufSize = %zu\n", bufUnit, walkBufSize / sizeof(uint), walkBufSize);
  fflush(stdout);
#endif
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
#ifdef  DPADD_FOR_ACC
 , DPacc * RESTRICT tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
 , acceleration * RESTRICT res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
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
#ifdef  DPADD_FOR_ACC
    const DPacc dac = {0.0, 0.0, 0.0, 0.0};
    tmp[info.head + lane] = dac;
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
    res[info.head + lane] = ai;
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
    //---------------------------------------------------------------------
  }/* if( lane < info.num ){ */
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
#ifdef  DPADD_FOR_ACC
 , DPacc * RESTRICT tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
 , acceleration * RESTRICT res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
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
#ifdef  DPADD_FOR_ACC
    const DPacc dac = {0.0, 0.0, 0.0, 0.0};
    tmp[GLOBALIDX_X1D] = dac;
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
  res[GLOBALIDX_X1D] = ai;
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
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
__global__ void trimAcc_kernel(acceleration * RESTRICT acc, READ_ONLY position * RESTRICT pos, const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo
#ifdef  DPADD_FOR_ACC
 , READ_ONLY DPacc * RESTRICT tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
			       , READ_ONLY acceleration * RESTRICT res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
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
    const int ii = info.head + lane;
    //---------------------------------------------------------------------
#ifndef DPADD_FOR_ACC
    /* load acceleration */
    acceleration ai = acc[ii];
    /* eliminate self-interaction */
    ai.pot -= epsinv * pos[ii].m;
#endif//DPADD_FOR_ACC
    //---------------------------------------------------------------------
#ifdef  DPADD_FOR_ACC
    DPacc dacc = tmp[ii];
    acceleration ai;
    ai.x   = CAST_D2R(dacc.x);
    ai.y   = CAST_D2R(dacc.y);
    ai.z   = CAST_D2R(dacc.z);
    ai.pot = CAST_D2R(dacc.pot - CAST_R2D(epsinv * pos[ii].m));
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
    acceleration corr = res[ii];
    ai.x   = CAST_D2R(CAST_R2D(ai.x  ) + CAST_R2D(corr.x  ));
    ai.y   = CAST_D2R(CAST_R2D(ai.y  ) + CAST_R2D(corr.y  ));
    ai.z   = CAST_D2R(CAST_R2D(ai.z  ) + CAST_R2D(corr.z  ));
    ai.pot = CAST_D2R(CAST_R2D(ai.pot) + CAST_R2D(corr.pot));
#if 0
    printf("res(%e) = %e, %e, %e, %e\n", pos[ii].x, corr.x, corr.y, corr.z, corr.pot);
#endif
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
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
__global__ void trimAcc_kernel(acceleration * RESTRICT acc, READ_ONLY position * RESTRICT pos
#ifdef  DPADD_FOR_ACC
			       , READ_ONLY DPacc * RESTRICT tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
			       , READ_ONLY acceleration * RESTRICT res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
			       )
{
  //-----------------------------------------------------------------------
  const int ii = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
#ifndef DPADD_FOR_ACC
  /* load acceleration and mass */
  acceleration ai = acc[ii];
  /* eliminate self-interaction */
  ai.pot -= epsinv * pos[ii].m;
#endif//DPADD_FOR_ACC
  //-----------------------------------------------------------------------
#ifdef  DPADD_FOR_ACC
  DPacc dacc = tmp[ii];
  acceleration ai;
  ai.x   = CAST_D2R(dacc.x);
  ai.y   = CAST_D2R(dacc.y);
  ai.z   = CAST_D2R(dacc.z);
  ai.pot = CAST_D2R(dacc.pot - CAST_R2D(epsinv * pos[ii].m));
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
  acceleration corr = res[ii];
  ai.x   = CAST_D2R(CAST_R2D(ai.x  ) + CAST_R2D(corr.x  ));
  ai.y   = CAST_D2R(CAST_R2D(ai.y  ) + CAST_R2D(corr.y  ));
  ai.z   = CAST_D2R(CAST_R2D(ai.z  ) + CAST_R2D(corr.z  ));
  ai.pot = CAST_D2R(CAST_R2D(ai.pot) + CAST_R2D(corr.pot));
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
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
/* #   if  TSUB >= ( 2 * NWARP) */
/*   real tmp; */
/*   tmp = __shfl_xor(val,      NWARP, TSUB);  if( tmp > val )    val = tmp; */
/* #   if  TSUB >= ( 4 * NWARP) */
/*   tmp = __shfl_xor(val,  2 * NWARP, TSUB);  if( tmp > val )    val = tmp; */
/* #   if  TSUB >= ( 8 * NWARP) */
/*   tmp = __shfl_xor(val,  4 * NWARP, TSUB);  if( tmp > val )    val = tmp; */
/* #   if  TSUB >= (16 * NWARP) */
/*   tmp = __shfl_xor(val,  8 * NWARP, TSUB);  if( tmp > val )    val = tmp; */
/* #   if  TSUB == (32 * NWARP) */
/*   tmp = __shfl_xor(val, 16 * NWARP, TSUB);  if( tmp > val )    val = tmp; */
/* #endif//TSUB == (32 * NWARP) */
/* #endif//TSUB >= (16 * NWARP) */
/* #endif//TSUB >= ( 8 * NWARP) */
/* #endif//TSUB >= ( 4 * NWARP) */
/* #endif//TSUB >= ( 2 * NWARP) */
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = __shfl_xor(val,      NWARP, TSUB);  val = FMAX(val, tmp);
#   if  TSUB >= ( 4 * NWARP)
  tmp = __shfl_xor(val,  2 * NWARP, TSUB);  val = FMAX(val, tmp);
#   if  TSUB >= ( 8 * NWARP)
  tmp = __shfl_xor(val,  4 * NWARP, TSUB);  val = FMAX(val, tmp);
#   if  TSUB >= (16 * NWARP)
  tmp = __shfl_xor(val,  8 * NWARP, TSUB);  val = FMAX(val, tmp);
#   if  TSUB == (32 * NWARP)
  tmp = __shfl_xor(val, 16 * NWARP, TSUB);  val = FMAX(val, tmp);
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
  return (__shfl(val, 0, TSUB));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx].r = max;
  //-----------------------------------------------------------------------
/* #   if  TSUB >= ( 2 * NWARP) */
/*   real tmp; */
/*   tmp = smem[tidx ^ (     NWARP)].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max; */
/* #   if  TSUB >= ( 4 * NWARP) */
/*   tmp = smem[tidx ^ ( 2 * NWARP)].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max; */
/* #   if  TSUB >= ( 8 * NWARP) */
/*   tmp = smem[tidx ^ ( 4 * NWARP)].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max; */
/* #   if  TSUB >= (16 * NWARP) */
/*   tmp = smem[tidx ^ ( 8 * NWARP)].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max; */
/* #   if  TSUB == (32 * NWARP) */
/*   tmp = smem[tidx ^ (16 * NWARP)].r;  if( tmp > max ){    max = tmp;  }  smem[tidx].r = max; */
/* #endif//TSUB == (32 * NWARP) */
/* #endif//TSUB >= (16 * NWARP) */
/* #endif//TSUB >= ( 8 * NWARP) */
/* #endif//TSUB >= ( 4 * NWARP) */
/* #endif//TSUB >= ( 2 * NWARP) */
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = smem[tidx ^ (     NWARP)].r;  max = FMAX(max, tmp);  smem[tidx].r = max;
#   if  TSUB >= ( 4 * NWARP)
  tmp = smem[tidx ^ ( 2 * NWARP)].r;  max = FMAX(max, tmp);  smem[tidx].r = max;
#   if  TSUB >= ( 8 * NWARP)
  tmp = smem[tidx ^ ( 4 * NWARP)].r;  max = FMAX(max, tmp);  smem[tidx].r = max;
#   if  TSUB >= (16 * NWARP)
  tmp = smem[tidx ^ ( 8 * NWARP)].r;  max = FMAX(max, tmp);  smem[tidx].r = max;
#   if  TSUB == (32 * NWARP)
  tmp = smem[tidx ^ (16 * NWARP)].r;  max = FMAX(max, tmp);  smem[tidx].r = max;
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
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
/* #   if  TSUB >= ( 2 * NWARP) */
/*   real tmp; */
/*   tmp = __shfl_xor(val,      NWARP, TSUB);  if( tmp < val )    val = tmp; */
/* #   if  TSUB >= ( 4 * NWARP) */
/*   tmp = __shfl_xor(val,  2 * NWARP, TSUB);  if( tmp < val )    val = tmp; */
/* #   if  TSUB >= ( 8 * NWARP) */
/*   tmp = __shfl_xor(val,  4 * NWARP, TSUB);  if( tmp < val )    val = tmp; */
/* #   if  TSUB >= (16 * NWARP) */
/*   tmp = __shfl_xor(val,  8 * NWARP, TSUB);  if( tmp < val )    val = tmp; */
/* #   if  TSUB == (32 * NWARP) */
/*   tmp = __shfl_xor(val, 16 * NWARP, TSUB);  if( tmp < val )    val = tmp; */
/* #endif//TSUB == (32 * NWARP) */
/* #endif//TSUB >= (16 * NWARP) */
/* #endif//TSUB >= ( 8 * NWARP) */
/* #endif//TSUB >= ( 4 * NWARP) */
/* #endif//TSUB >= ( 2 * NWARP) */
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = __shfl_xor(val,      NWARP, TSUB);  val = FMIN(val, tmp);
#   if  TSUB >= ( 4 * NWARP)
  tmp = __shfl_xor(val,  2 * NWARP, TSUB);  val = FMIN(val, tmp);
#   if  TSUB >= ( 8 * NWARP)
  tmp = __shfl_xor(val,  4 * NWARP, TSUB);  val = FMIN(val, tmp);
#   if  TSUB >= (16 * NWARP)
  tmp = __shfl_xor(val,  8 * NWARP, TSUB);  val = FMIN(val, tmp);
#   if  TSUB == (32 * NWARP)
  tmp = __shfl_xor(val, 16 * NWARP, TSUB);  val = FMIN(val, tmp);
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
  return (__shfl(val, 0, TSUB));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
  smem[tidx].r = min;
  //-----------------------------------------------------------------------
/* #   if  TSUB >= ( 2 * NWARP) */
/*   real tmp; */
/*   tmp = smem[tidx ^ (     NWARP)].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min; */
/* #   if  TSUB >= ( 4 * NWARP) */
/*   tmp = smem[tidx ^ ( 2 * NWARP)].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min; */
/* #   if  TSUB >= ( 8 * NWARP) */
/*   tmp = smem[tidx ^ ( 4 * NWARP)].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min; */
/* #   if  TSUB >= (16 * NWARP) */
/*   tmp = smem[tidx ^ ( 8 * NWARP)].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min; */
/* #   if  TSUB == (32 * NWARP) */
/*   tmp = smem[tidx ^ (16 * NWARP)].r;  if( tmp < min ){    min = tmp;  }  smem[tidx].r = min; */
/* #endif//TSUB == (32 * NWARP) */
/* #endif//TSUB >= (16 * NWARP) */
/* #endif//TSUB >= ( 8 * NWARP) */
/* #endif//TSUB >= ( 4 * NWARP) */
/* #endif//TSUB >= ( 2 * NWARP) */
#   if  TSUB >= ( 2 * NWARP)
  real tmp;
  tmp = smem[tidx ^ (     NWARP)].r;  min = FMIN(min, tmp);  smem[tidx].r = min;
#   if  TSUB >= ( 4 * NWARP)
  tmp = smem[tidx ^ ( 2 * NWARP)].r;  min = FMIN(min, tmp);  smem[tidx].r = min;
#   if  TSUB >= ( 8 * NWARP)
  tmp = smem[tidx ^ ( 4 * NWARP)].r;  min = FMIN(min, tmp);  smem[tidx].r = min;
#   if  TSUB >= (16 * NWARP)
  tmp = smem[tidx ^ ( 8 * NWARP)].r;  min = FMIN(min, tmp);  smem[tidx].r = min;
#   if  TSUB == (32 * NWARP)
  tmp = smem[tidx ^ (16 * NWARP)].r;  min = FMIN(min, tmp);  smem[tidx].r = min;
#endif//TSUB == (32 * NWARP)
#endif//TSUB >= (16 * NWARP)
#endif//TSUB >= ( 8 * NWARP)
#endif//TSUB >= ( 4 * NWARP)
#endif//TSUB >= ( 2 * NWARP)
  //-----------------------------------------------------------------------
  return (smem[head].r);
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_s2s(uint *src, int sidx, uint *dst, int didx, const int num, const int lane)
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
__device__ __forceinline__ void copyData_g2s(uint * RESTRICT gbuf, size_t srcHead, uint * RESTRICT sbuf, int dstHead, int numCopy, const int lane)
{
  //-----------------------------------------------------------------------
  /* fraction processing at loading from the head of destination array */
  //-----------------------------------------------------------------------
  const int numTemp = TSUB - (int)(srcHead & (TSUB - 1));/* := TSUB - (srcHead % TSUB) */
  const int numHead = (numTemp < numCopy) ? numTemp : numCopy;
  if( lane < numHead )
    sbuf[dstHead + lane] = gbuf[srcHead + lane];
  dstHead += numHead;
  srcHead += numHead;
  numCopy -= numHead;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* sequential load from source on the global memory and store to destination on the shared memory */
  //-----------------------------------------------------------------------
  for(int ii = lane; ii < numCopy; ii += TSUB)
    sbuf[dstHead + ii] = gbuf[srcHead + ii];
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_s2g(uint * RESTRICT sbuf, int srcHead, uint * RESTRICT gbuf, size_t dstHead, int numCopy, const int lane)
{
  //-----------------------------------------------------------------------
  /* fraction processing at storing to the head of destination array */
  //-----------------------------------------------------------------------
  const int numTemp = TSUB - (int)(dstHead & (TSUB - 1));/* := TSUB - (dstHead % TSUB) */
  const int numHead = (numTemp < numCopy) ? numTemp : numCopy;
  if( lane < numHead )
    gbuf[dstHead + lane] = sbuf[srcHead + lane];
  dstHead += numHead;
  srcHead += numHead;
  numCopy -= numHead;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* sequential load from source on the shared memory and store to destination on the global memory */
  //-----------------------------------------------------------------------
  for(int ii = lane; ii < numCopy; ii += TSUB)
    gbuf[dstHead + ii] = sbuf[srcHead + ii];
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void copyData_g2g(uint * RESTRICT gbuf, size_t srcHead, size_t dstHead, int Ncopy, const int Ndisp, const int lane)
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
    copyData_g2g(gmbuf, hb, hb, *num_gm, *head_gm, lane);
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
  copyData_s2g(node, Nsm, gmbuf, hb + (*tail_gm), Ntot, lane);
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
#   if  defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 600)
/* after Pascal generation, native atomicAdd for FP64 is provided */
__device__ __forceinline__ double atomicAdd(double* address, double val)
{
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif//defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 600)
/* __device__ __forceinline__ float atomicAdd(float* addr, double val) */
/* { */
/*   unsigned int* addr_as_ui = (unsigned int*)addr; */
/*   unsigned int old = *addr_as_ui, assumed; */
/*   do { */
/*     assumed = old; */
/*     old = atomicCAS(addr_as_ui, assumed, __float_as_int((float)(val + (double)__int_as_float(assumed)))); */
/*   } while (assumed != old); */
/*   return __int_as_float(old); */
/* } */
/* __device__ __forceinline__ float atomicPrecAdd(float* addr, double val) */
/* { */
/*   int* addr_as_i = (int*)addr; */
/*   int old = *addr_as_i, assumed; */
/*   do { */
/*     assumed = old; */
/*     old = atomicCAS(addr_as_i, assumed, __float_as_int((float)(val + (double)__int_as_float(assumed)))); */
/*   } while (assumed != old); */
/*   return __int_as_float(old); */
/* } */
//-------------------------------------------------------------------------
__global__ void __launch_bounds__(NTHREADS, NBLOCKS_PER_SM) calcAcc_kernel
     (READ_ONLY laneinfo * RESTRICT laneInfo, READ_ONLY position * RESTRICT ipos, jnode * RESTRICT iacc,
#ifdef  GADGET_MAC
      READ_ONLY acceleration * RESTRICT iacc_old,
#endif//GADGET_MAC
      const int root, READ_ONLY uint * RESTRICT more, READ_ONLY jparticle * RESTRICT jpos, READ_ONLY jmass * RESTRICT mj,
#ifdef  DPADD_FOR_ACC
      DPacc * RESTRICT dacc,
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
      jnode * RESTRICT ires,
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
#   if  !defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
      int * RESTRICT active, uint * RESTRICT freeNum,
#endif//!defined(USE_SMID_TO_GET_BUFID) && !defined(TRY_MODE_ABOUT_BUFFER)
#   if  !defined(USE_SMID_TO_GET_BUFID) &&  defined(TRY_MODE_ABOUT_BUFFER)
      const int freeNum,
#endif//!defined(USE_SMID_TO_GET_BUFID) &&  defined(TRY_MODE_ABOUT_BUFFER)
      uint * RESTRICT freeLst, uint * RESTRICT buffer, const int bufSize, int * RESTRICT overflow
#   if  !defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
      , unsigned long long int * RESTRICT cycles
#endif//!defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
#ifdef  COUNT_INTERACTIONS
      , int * RESTRICT stockNj, int * RESTRICT stockNbuf
#endif//COUNT_INTERACTIONS
      )
{
  //-----------------------------------------------------------------------
  /* start stop watch */
  //-----------------------------------------------------------------------
#   if  !defined(USE_CUDA_EVENT) && !defined(SERIALIZED_EXECUTION) && !defined(PRINT_PSEUDO_PARTICLE_INFO)
  const long long int initCycle = clock64();
#endif//!defined(USE_CUDA_EVENT) && !defined(SERIALIZED_EXECUTION) && !defined(PRINT_PSEUDO_PARTICLE_INFO)
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
#if 0
  for(int ii = tidx; ii < NTHREADS * (NLOOP + 1); ii += NTHREADS){
    const position massless = {ZERO, ZERO, ZERO, ZERO};
    pj[ii].pi = massless;
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
    eps2[ii] = ZERO;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
  }
  for(int ii = tidx; ii < NTHREADS * NQUEUE; ii += NTHREADS)
    queue[ii] = NULL_NODE;
  __syncthreads();
#endif
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
  size_t buf0Head = (size_t)(bufIdx * NGROUPS + DIV_TSUB(head)) * (size_t)bufSize;
#ifdef  DBG_TREE_WALK
  if( lane == 0 )
    printf("buffer: %ld--%ld\n", buf0Head, buf0Head + bufSize - 1);
#endif//DBG_TREE_WALK
  //-----------------------------------------------------------------------
#if 0
  if( tidx == 0 )
    printf("walk: SM %d, tag %d\n", bufTarget / NBLOCKS_PER_SM, bufTarget);
#endif
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
#ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
  bool writeWholeData = false;
  if( !skip )
    if(     (-10.270200f < pi.x) && (pi.x < -10.270000f) )
      if(   (  0.646161f < pi.y) && (pi.y <   0.646163f) )
	if( (- 0.884657f < pi.z) && (pi.z < - 0.884655f) )
	  writeWholeData = true;
#if 1
#   if  TSUB >=  2
  writeWholeData |= __shfl_xor(writeWholeData,  1, TSUB);
#   if  TSUB >=  4
  writeWholeData |= __shfl_xor(writeWholeData,  2, TSUB);
#   if  TSUB >=  8
  writeWholeData |= __shfl_xor(writeWholeData,  4, TSUB);
#   if  TSUB >= 16
  writeWholeData |= __shfl_xor(writeWholeData,  8, TSUB);
#   if  TSUB == 32
  writeWholeData |= __shfl_xor(writeWholeData, 16, TSUB);
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
#else
  if( jtag != 0 )
    writeWholeData = false;
#endif
#if 0
  if( lane != 0 )
    writeWholeData = false;
#endif
#endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
  //-----------------------------------------------------------------------
#if 0
  if( (idx == 0) && (tidx == 0) )
    printf("pi_x = %e; pj_x = %e\n", pi.x, jpos[0].x);
  /* if( (idx == 0) && (tidx == 0) ) */
  /*   printf("G = %e, eps2 = %e\n", newton, eps2); */
#endif
  //-----------------------------------------------------------------------
#if 0
  if( (idx == 0) && (tidx == 0) )
    printf("(%e @ %e, %e, %e)\n", mj[0], jpos[0].x, jpos[0].y, jpos[0].z);
#endif
  //-----------------------------------------------------------------------
  int fail = hp + lane;
  jnode jidx;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#   if  !defined(USE_CUDA_EVENT) && defined(PRINT_PSEUDO_PARTICLE_INFO)
  const long long int initCycle = clock64();
#endif//!defined(USE_CUDA_EVENT) && defined(PRINT_PSEUDO_PARTICLE_INFO)
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
	    /* , smem, (int *)&queue[hq - tidx], tidx, head */
	    , smem, (int *)&queue[hq - lane], tidx, head
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
	    /* , smem, (int *)&queue[hq - tidx], tidx, head */
	    , smem, (int *)&queue[hq - lane], tidx, head
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
    ball.m = getMaximumRealTsub(FLT_MIN + dx * dx + dy * dy + dz * dz
#ifndef USE_WARP_SHUFFLE_FUNC
				, smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
				);
    /* calculate size of mass-weighted estimated enclosing ball */
    dx = pi.x - icom.x;
    dy = pi.y - icom.y;
    dz = pi.z - icom.z;
    icom.m = getMaximumRealTsub(FLT_MIN + dx * dx + dy * dy + dz * dz
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
#ifdef  YMIKI_MAC
  acceleration amin;
#else///YMIKI_MAC
  real amin;
#endif//YMIKI_MAC
#endif//GADGET_MAC
  {
    //---------------------------------------------------------------------
    /* calculate displacement of i-particle and center-of-mass */
    const real rx = pi.x - icom.x;
    const real ry = pi.y - icom.y;
    const real rz = pi.z - icom.z;
    //---------------------------------------------------------------------
    /* calculate maximum of r squared */
    icom.m = getMaximumRealTsub(FLT_MIN + rx * rx + ry * ry + rz * rz
#ifndef USE_WARP_SHUFFLE_FUNC
				, smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
				);
    //---------------------------------------------------------------------
#ifdef  GADGET_MAC
    //---------------------------------------------------------------------
    acceleration ai_old = {ZERO, ZERO, ZERO, ZERO};
    if( !skip ){
      ai_old = iacc_old[idx];
/* #ifdef  YMIKI_MAC */
/*       ai_old.pot = UNITY; */
/* #endif//YMIKI_MAC */
    }/* if( !skip ){ */
    //---------------------------------------------------------------------
    /* calculate minimum of a squared */
    const real tmp = getMinimumRealTsub(FLT_MIN + ai_old.x * ai_old.x + ai_old.y * ai_old.y + ai_old.z * ai_old.z
#ifndef USE_WARP_SHUFFLE_FUNC
					, smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC
					);
#ifndef YMIKI_MAC
    amin = tmp * RSQRT(tmp);
#endif//YMIKI_MAC
    //---------------------------------------------------------------------
#ifdef  YMIKI_MAC
    /* calculate bulk acceleration of a group of i-particles */
    pj[fail].ai = ai_old;
    /* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
#   if  TSUB >=  2
#                    if  NWARP <  2
    jidx = pj[fail ^  1];    ai_old.x += jidx.ai.x;    ai_old.y += jidx.ai.y;    ai_old.z += jidx.ai.z;    /* ai_old.pot += jidx.ai.pot; */    pj[fail].ai = ai_old;
#                 endif//NWARP <  2
#   if  TSUB >=  4
#                    if  NWARP <  4
    jidx = pj[fail ^  2];    ai_old.x += jidx.ai.x;    ai_old.y += jidx.ai.y;    ai_old.z += jidx.ai.z;    /* ai_old.pot += jidx.ai.pot; */    pj[fail].ai = ai_old;
#                 endif//NWARP <  4
#   if  TSUB >=  8
#                    if  NWARP <  8
    jidx = pj[fail ^  4];    ai_old.x += jidx.ai.x;    ai_old.y += jidx.ai.y;    ai_old.z += jidx.ai.z;    /* ai_old.pot += jidx.ai.pot; */    pj[fail].ai = ai_old;
#                 endif//NWARP <  8
#   if  TSUB >= 16
#                    if  NWARP < 16
    jidx = pj[fail ^  8];    ai_old.x += jidx.ai.x;    ai_old.y += jidx.ai.y;    ai_old.z += jidx.ai.z;    /* ai_old.pot += jidx.ai.pot; */    pj[fail].ai = ai_old;
#                 endif//NWARP < 16
#   if  TSUB == 32
#                    if  NWARP < 32
    jidx = pj[fail ^ 16];    ai_old.x += jidx.ai.x;    ai_old.y += jidx.ai.y;    ai_old.z += jidx.ai.z;    /* ai_old.pot += jidx.ai.pot; */    pj[fail].ai = ai_old;
#                 endif//NWARP < 32
#endif//TSUB == 32
#endif//TSUB >= 16
#endif//TSUB >=  8
#endif//TSUB >=  4
#endif//TSUB >=  2
    amin = pj[hp].ai;
#if 1
    /* amin.pot = UNITY / amin.pot; */
    /* amin.x *= amin.pot; */
    /* amin.y *= amin.pot; */
    /* amin.z *= amin.pot; */
    amin.pot = SQRTRATIO(tmp, FLT_MIN + amin.x * amin.x + amin.y * amin.y + amin.z * amin.z);
    amin.x *= amin.pot;
    amin.y *= amin.pot;
    amin.z *= amin.pot;
#else
    amin.pot = SQRTRATIO(tmp * amin.pot * amin.pot, FLT_MIN + amin.x * amin.x + amin.y * amin.y + amin.z * amin.z);
    amin.x *= amin.pot;
    amin.y *= amin.pot;
    amin.z *= amin.pot;
#endif
#endif//YMIKI_MAC
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
#   if  !defined(USE_CUDA_EVENT) && defined(PRINT_PSEUDO_PARTICLE_INFO)
  long long int exitCycle = clock64();
  if( lane == 0 ){
    unsigned long long int elapsed = (unsigned long long int)(exitCycle - initCycle);
    atomicAdd(cycles, elapsed);
    iacc[DIV_TSUB(GLOBALIDX_X1D)].pi = icom;
  }/* if( tidx == 0 ){ */
#endif//!defined(USE_CUDA_EVENT) && defined(PRINT_PSEUDO_PARTICLE_INFO)
  //-----------------------------------------------------------------------
/* #if 0 */
/*   amin = ZERO; */
/* #endif */
/* #ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
/*   if( writeWholeData && (jtag == 0) ) */
/*     printf("amin = %e, root = %d, Mj = %e\n", amin, root, mj[root]); */
/* #endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
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
#if 0
  int bufTailMax = bufTail;
#endif
  /* set child j-cells in queue on the shared memory */
#if 1
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
#else
  uint jcell;
  int rem = 1;
  if( lane == 0 )
    queue[hq] = 0;
#endif
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
/* #ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
/*     if( writeWholeData ){ */
/*       printf("%u+%u", queue[hq] & IDXMASK, 1 + (queue[hq] >> IDXBITS)); */
/*       for(int ii = 1; ii < rem; ii++) */
/* 	printf("\t%u+%u", queue[hq + ii] & IDXMASK, 1 + (queue[hq + ii] >> IDXBITS)); */
/*       printf("\n"); */
/*     } */
/* #endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
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
	const real r2 = FLT_MIN + rx * rx + ry * ry + rz * rz;
	real lambda = FMAX(UNITY - SQRTRATIO(icom.m, r2), ZERO);
	/* calculate distance between the pseudo i-particle and the candidate j-particle */
	//-----------------------------------------------------------------
#ifndef YMIKI_MAC
	lambda *= lambda * r2;
#endif//YMIKI_MAC
#ifdef  GADGET_MAC
#ifndef YMIKI_MAC
	/* alpha * |a| * r^4 > G * M * l^2 */
	if( jcnd.w < lambda * lambda * amin )
#else///YMIKI_MAC
	/* alpha * |(a, r)| * r^3 > G * M * l^2 */
	lambda *= lambda;	lambda *= lambda;	          /* lambda := lambda^4 */
	lambda *= (amin.x * rx + amin.y * ry + amin.z * rz);      /* lambda := lambda^4 * (ai, rij) */
	lambda *= r2;                                             /* lambda := lambda^4 * (ai, rij) * d^2 */
	if( jcnd.w < lambda * lambda * r2 )
#endif//YMIKI_MAC
#else///GADGET_MAC
#ifdef  WS93_MAC
	  if(   jcnd.w < lambda )
#else///WS93_MAC
	    /* (l / r) < theta */
	    if( jcnd.w < lambda * theta2 )
#endif//WS93_MAC
#endif//GADGET_MAC
	      {
		/* add the candidate j-particle to the interaction list */
#ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
		if( writeWholeData )
		  printf("%e\t%e\t%e\t%e\t%e\n", mj[target], rx, ry, rz, r2);
		  /* printf("%e\t%e\t%e\t%e\t%e\t%d\n", mj[target], jcnd.x, jcnd.y, jcnd.z, jcnd.w, target); */
#endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
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
/* #ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
/* 		if( writeWholeData ) */
/* 		  printf("%e\t%e\t%e\t%e\t%e\n", mj[target], jcnd.x, jcnd.y, jcnd.z, jcnd.w); */
/* #endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
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
/* #ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
/* 	if( writeWholeData ) */
/* 	  for(int ll = 0; ll < NLOOP * TSUB; ll++) */
/* 	    printf("%e\t%e\t%e\t%e\n", pj[hp + ll].pos.w, pj[hp + ll].pos.x, pj[hp + ll].pos.y, pj[hp + ll].pos.z); */
/* #endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
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
      copyData_g2s(buffer, buf0Head + bufHead, queue, hq - lane + rem, Nload, lane);
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
    /* fail += (bufOpen < 0); */
    fail += (bufTail > bufSize);
#if 0
    if( fail != 0 )
      buffer[ULONG_MAX] = NULL_NODE;
#endif
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
#if 0
    if( bufTail > bufTailMax )
      bufTailMax = bufTail;
#endif
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
#if 0
  if( lane == 0 )
    printf("bufTailMax = %d\n", bufTailMax);
#endif
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
    //---------------------------------------------------------------------
/* #ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
/*     if( writeWholeData ) */
/*       for(int ll = 0; ll < NLOOP * TSUB; ll++) */
/* 	printf("%e\t%e\t%e\t%e\n", pj[hp + ll].pos.w, pj[hp + ll].pos.x, pj[hp + ll].pos.y, pj[hp + ll].pos.z); */
/* #endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION */
    //---------------------------------------------------------------------
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
  /* accumulation of residuals in Kahan summation */
#   if  defined(ACCURATE_ACCUMULATION) && (NWARP > 1)
  /* set index to accumulate acceleration without atomic operations */
#   if  NWARP > 4
  const int gtag = jtag >> 2;
  jtag &= 3;
#endif//NWARP > 4
  const int itag = hp + lane - jtag;
  pj[hp + lane].ai = res;
#   if  NWARP == 2
  pj[itag].val[jtag] += pj[itag + 1].val[jtag];  jtag ^= 2;
  pj[itag].val[jtag] += pj[itag + 1].val[jtag];
#endif//NWARP == 2
#   if  NWARP == 4
  pj[itag].val[jtag] += pj[itag + 1].val[jtag] + pj[itag + 2].val[jtag] + pj[itag + 3].val[jtag];
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
    pj[itag].val[jtag] += pj[itag + 4].val[jtag];
#endif//NWARP >= 8
  res = pj[itag].ai;
#endif//defined(ACCURATE_ACCUMULATION) && (NWARP > 1)
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
#   if  NWARP > 1
#ifndef ACCURATE_ACCUMULATION
    /* set index to accumulate acceleration without atomic operations */
#   if  NWARP > 4
    const int gtag = jtag >> 2;
    jtag &= 3;
#endif//NWARP > 4
    const int itag = hp + lane - jtag;
#endif//ACCURATE_ACCUMULATION
    pj[hp + lane].ai = ai;
#   if  NWARP == 2
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
    /* T := S, S := S + R; T := S - T; R := R - T */
    real sum0 = pj[itag].val[jtag] + pj[itag + 1].val[jtag];    real tmp0 = atomicAdd(&(iacc[idx].val[jtag]), sum0);    tmp0 = (tmp0 + sum0) - tmp0;    sum0 -= tmp0;    jtag ^= 2;
    real sum1 = pj[itag].val[jtag] + pj[itag + 1].val[jtag];    real tmp1 = atomicAdd(&(iacc[idx].val[jtag]), sum1);    tmp1 = (tmp1 + sum1) - tmp1;    sum1 -= tmp1;
    pj[hp + lane].ai = res;
    atomicAdd(&(ires[idx].val[jtag]), sum1 + pj[itag].val[jtag]);    jtag ^= 2;
    atomicAdd(&(ires[idx].val[jtag]), sum0 + pj[itag].val[jtag]);
#else///defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
    atomicAdd(&(iacc[idx].val[jtag]), pj[itag].val[jtag] + pj[itag + 1].val[jtag]);    jtag ^= 2;
    atomicAdd(&(iacc[idx].val[jtag]), pj[itag].val[jtag] + pj[itag + 1].val[jtag]);
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
#endif//NWARP == 2
#   if  NWARP == 4
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
    /* T := S, S := S + R; T := S - T; R := R - T */
    real sum = pj[itag].val[jtag] + pj[itag + 1].val[jtag] + pj[itag + 2].val[jtag] + pj[itag + 3].val[jtag];
#if 1
    real tmp = atomicAdd(&(iacc[idx].val[jtag]), sum);    tmp = (tmp + sum) - tmp;    sum -= tmp;
#else
    real tmp = atomicAdd(&(iacc[idx].val[jtag]), sum);    tmp = CAST_D2R(CAST_R2D(tmp + sum) - CAST_R2D(tmp));    sum = CAST_D2R(CAST_R2D(sum) - CAST_R2D(tmp));
#endif
    pj[hp + lane].ai = res;
    atomicAdd(&(ires[idx].val[jtag]), sum + pj[itag].val[jtag]);
#else///defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
#if 0
    if( jtag == 0 )
      printf("%e\t%e\t%e\t%e\t%e\n", pi.x,
	     pj[itag].ai.x   + pj[itag + 1].ai.x   + pj[itag + 2].ai.x   + pj[itag + 3].ai.x,
	     pj[itag].ai.y   + pj[itag + 1].ai.y   + pj[itag + 2].ai.y   + pj[itag + 3].ai.y,
	     pj[itag].ai.z   + pj[itag + 1].ai.z   + pj[itag + 2].ai.z   + pj[itag + 3].ai.z,
	     pj[itag].ai.pot + pj[itag + 1].ai.pot + pj[itag + 2].ai.pot + pj[itag + 3].ai.pot);
#endif
    atomicAdd(&(iacc[idx].val[jtag]), pj[itag].val[jtag] + pj[itag + 1].val[jtag] + pj[itag + 2].val[jtag] + pj[itag + 3].val[jtag]);
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
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
    if( gtag == 0 ){
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
      /* T := S, S := S + R; T := S - T; R := R - T */
      real sum = pj[itag].val[jtag] + pj[itag + 4].val[jtag];
      real tmp = atomicAdd(&(iacc[idx].val[jtag]), sum);      tmp = (tmp + sum) - tmp;      sum -= tmp;
      pj[hp + lane].ai = res;
      atomicAdd(&(ires[idx].val[jtag]), sum + pj[itag].val[jtag]);
#else///defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
      atomicAdd(&(iacc[idx].val[jtag]), pj[itag].val[jtag] + pj[itag + 4].val[jtag]);
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
    }/* if( gtag == 0 ){ */
#endif//NWARP >= 8
#else///NWARP > 1
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
    /* /\* R := R + x_i *\/ */
    /* res.x += ai.x;    res.y += ai.y;    res.z += ai.z;    res.pot += ai.pot; */
    /* T := S, S := S + R; T := S - T; R := R - T */
    acceleration tmp;
    tmp.x   = atomicAdd(&(iacc[idx].ai.x  ), ai.x  );    tmp.x   = (tmp.x   + ai.x  ) - tmp.x  ;    ai.x   -= tmp.x  ;    atomicAdd(&(ires[idx].ai.x  ), ai.x   + res.x  );
    tmp.y   = atomicAdd(&(iacc[idx].ai.y  ), ai.y  );    tmp.y   = (tmp.y   + ai.y  ) - tmp.y  ;    ai.y   -= tmp.y  ;    atomicAdd(&(ires[idx].ai.y  ), ai.y   + res.y  );
    tmp.z   = atomicAdd(&(iacc[idx].ai.z  ), ai.z  );    tmp.z   = (tmp.z   + ai.z  ) - tmp.z  ;    ai.z   -= tmp.z  ;    atomicAdd(&(ires[idx].ai.z  ), ai.z   + res.z  );
    tmp.pot = atomicAdd(&(iacc[idx].ai.pot), ai.pot);    tmp.pot = (tmp.pot + ai.pot) - tmp.pot;    ai.pot -= tmp.pot;    atomicAdd(&(ires[idx].ai.pot), ai.pot + res.pot);
    /* acceleration tmp; */
    /* tmp.x   = atomicAdd(&(iacc[idx].ai.x  ), ai.x  );    tmp.x   = iacc[idx].ai.x   - tmp.x  ;    ai.x   -= tmp.x  ;    atomicAdd(&(ires[idx].ai.x  ), ai.x   + res.x  ); */
    /* tmp.y   = atomicAdd(&(iacc[idx].ai.y  ), ai.y  );    tmp.y   = iacc[idx].ai.y   - tmp.y  ;    ai.y   -= tmp.y  ;    atomicAdd(&(ires[idx].ai.y  ), ai.y   + res.y  ); */
    /* tmp.z   = atomicAdd(&(iacc[idx].ai.z  ), ai.z  );    tmp.z   = iacc[idx].ai.z   - tmp.z  ;    ai.z   -= tmp.z  ;    atomicAdd(&(ires[idx].ai.z  ), ai.z   + res.z  ); */
    /* tmp.pot = atomicAdd(&(iacc[idx].ai.pot), ai.pot);    tmp.pot = iacc[idx].ai.pot - tmp.pot;    ai.pot -= tmp.pot;    atomicAdd(&(ires[idx].ai.pot), ai.pot + res.pot); */
#else///defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
#ifndef DPADD_FOR_ACC
    atomicAdd(&(iacc[idx].ai.x  ), ai.x  );
    atomicAdd(&(iacc[idx].ai.y  ), ai.y  );
    atomicAdd(&(iacc[idx].ai.z  ), ai.z  );
    atomicAdd(&(iacc[idx].ai.pot), ai.pot);
#else///DPADD_FOR_ACC
    atomicAdd(&(dacc[idx].x  ), CAST_R2D(ai.x  ) + CAST_R2D(res.x  ));
    atomicAdd(&(dacc[idx].y  ), CAST_R2D(ai.y  ) + CAST_R2D(res.y  ));
    atomicAdd(&(dacc[idx].z  ), CAST_R2D(ai.z  ) + CAST_R2D(res.z  ));
    atomicAdd(&(dacc[idx].pot), CAST_R2D(ai.pot) + CAST_R2D(res.pot));
#endif//DPADD_FOR_ACC
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION)
#endif//NWARP > 1
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
  releaseBuffer(tidx, freeLst, (uint)bufIdx, bufTarget);
#else///USE_SMID_TO_GET_BUFID
#ifdef  TRY_MODE_ABOUT_BUFFER
  releaseBuffer(tidx, freeLst, (uint)bufIdx, bufTarget);
#else///TRY_MODE_ABOUT_BUFFER
  releaseBuffer(tidx, freeNum, freeLst, bufIdx, active);
#endif//TRY_MODE_ABOUT_BUFFER
#endif//USE_SMID_TO_GET_BUFID
  //-----------------------------------------------------------------------
#   if  !defined(USE_CUDA_EVENT) && !defined(SERIALIZED_EXECUTION) && !defined(PRINT_PSEUDO_PARTICLE_INFO)
  long long int exitCycle = clock64();
  if( tidx == 0 ){
    unsigned long long int elapsed = (unsigned long long int)(exitCycle - initCycle);
    atomicAdd(cycles, elapsed);
  }/* if( tidx == 0 ){ */
#endif//!defined(USE_CUDA_EVENT) && !defined(SERIALIZED_EXECUTION) && !defined(PRINT_PSEUDO_PARTICLE_INFO)
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#if 0
//-------------------------------------------------------------------------
__global__ void printTreeNode_kernel(const int num, READ_ONLY jparticle * RESTRICT pj, READ_ONLY jmass * RESTRICT mj)
{
  for(int ii = 0; ii < num; ii++){
    const jparticle jpos = pj[ii];
    const jmass     mass = mj[ii];
    printf("%e\t%e\t%e\t%e\n", mass, jpos.x, jpos.y, jpos.z);
  }/* for(int ii = 0; ii < num; ii++){ */
}
//-------------------------------------------------------------------------
void printTreeNode_dev(const int num, const soaTreeNode node, MPIinfo mpi)
{
  for(int ii = 0; ii < mpi.size; ii++){
    if( ii == mpi.rank )
      printTreeNode_kernel<<<1, 1>>>(num, node.jpos, node.mj);
    checkCudaErrors(cudaDeviceSynchronize());
    fflush(stdout);
    MPI_Barrier(mpi.comm);
    if( ii == mpi.rank )
      printf("# rank %d / %d\n", mpi.rank, mpi.size);
    fflush(stdout);
    MPI_Barrier(mpi.comm);
  }
  MPI_Finalize();
  exit(0);
}
//-------------------------------------------------------------------------
void printFullTree_dev(int n0, int n1, int h1, const soaTreeNode node, MPIinfo mpi)
{
  for(int ii = 0; ii < mpi.size; ii++){
    if( ii == mpi.rank ){
      printTreeNode_kernel<<<1, 1>>>(n0,   node.jpos     ,   node.mj     );
      printTreeNode_kernel<<<1, 1>>>(n1, &(node.jpos[h1]), &(node.mj[h1]));
    }
    checkCudaErrors(cudaDeviceSynchronize());
    fflush(stdout);
    MPI_Barrier(mpi.comm);
    if( ii == mpi.rank )
      printf("# rank %d / %d\n", mpi.rank, mpi.size);
    fflush(stdout);
    MPI_Barrier(mpi.comm);
  }
  MPI_Finalize();
  exit(0);
}
//-------------------------------------------------------------------------
__global__ void printTreeLink_kernel(READ_ONLY uint * RESTRICT more, READ_ONLY jparticle * RESTRICT pj, READ_ONLY jmass * RESTRICT mj, uint * RESTRICT buffer)
{
  buffer[0] = 0;
  int rem = 1;
  int head = rem;
  int ii = 0;

  while( true ){
    if( rem == 0 )
      break;

    uint jj = buffer[ii];
    /* printf("%e\t%e\t%e\t%e\t%e\n", mj[jj], pj[jj].x, pj[jj].y, pj[jj].z, pj[jj].w); */
    printf("%e\t%e\t%e\t%e\n", mj[jj], pj[jj].x, pj[jj].y, pj[jj].z);
    rem--;

    uint jcell = more[jj];
    int num = 1 + (jcell >> IDXBITS);
    jcell &= IDXMASK;
    if( jcell == jj )
      num = 0;
    for(int kk = 0; kk < num; kk++)
      buffer[head + kk] = jcell + kk;

    head += num;
    rem += num;
    ii++;
  }
}
//-------------------------------------------------------------------------
 void printTreeLink_dev(int head, const soaTreeNode node, uint * buf, MPIinfo mpi)
{
  checkCudaErrors(cudaDeviceSynchronize());
  MPI_Barrier(mpi.comm);
  for(int ii = 0; ii < mpi.size; ii++){
    if( ii == mpi.rank ){
      printTreeLink_kernel<<<1, 1>>>(  node.more       ,   node.jpos     ,     node.mj       , buf);
      printTreeLink_kernel<<<1, 1>>>(&(node.more[head]), &(node.jpos[head]), &(node.mj[head]), buf);
    }
    checkCudaErrors(cudaDeviceSynchronize());
    fflush(stdout);
    MPI_Barrier(mpi.comm);
    if( ii == mpi.rank )
      printf("# rank %d / %d\n", mpi.rank, mpi.size);
    fflush(stdout);
    MPI_Barrier(mpi.comm);
  }
  MPI_Finalize();
  exit(0);
}
//-------------------------------------------------------------------------
#endif
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
 , const int grpNum, const int jhead
#endif//SERIALIZED_EXECUTION
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  USE_CUDA_EVENT
 , int *Nwalk, cudaEvent_t *iniEvent, cudaEvent_t *finEvent
#else///USE_CUDA_EVENT
 , unsigned long long int * RESTRICT cycles
#endif//USE_CUDA_EVENT
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
 , const soaTreeWalkBuf buf
#ifdef  COUNT_INTERACTIONS
 , iparticle_treeinfo treeInfo
#endif//COUNT_INTERACTIONS
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s (grpNum = %d, jhead = %d)\n", "start", grpNum, jhead);
  //-----------------------------------------------------------------------
#ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
  checkCudaErrors(cudaDeviceSynchronize());
  chkMPIerr(MPI_Barrier(MPI_COMM_WORLD));
#endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
  //-----------------------------------------------------------------------
#if 0
  int deviceID;
  checkCudaErrors(cudaGetDevice(&deviceID));
  fprintf(stdout, "jhead = %d on device %d\n", jhead, deviceID);
  fflush(stdout);
#endif
  //-----------------------------------------------------------------------
#   if  defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
  checkCudaErrors(cudaEventRecord(iniEvent[*Nwalk], 0));
#endif//defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
  //-----------------------------------------------------------------------
#   if  defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
  if( grpNum != 0 ){
#endif//defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
    //---------------------------------------------------------------------
    if( blck.x <= MAX_BLOCKS_PER_GRID ){
      //-------------------------------------------------------------------
#if 1
      calcAcc_kernel<<<blck, thrd, SMEM_SIZE, sinfo->stream[*sidx]>>>
#else
      calcAcc_kernel<<<blck, thrd>>>
#endif
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
	 &(tree.more[jhead]), &(tree.jpos[jhead]), &(tree.mj[jhead]),
#endif//SERIALIZED_EXECUTION
#ifdef  DPADD_FOR_ACC
	 pi.tmp,
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 (jnode *)pi.res,
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
#ifndef USE_SMID_TO_GET_BUFID
#ifndef TRY_MODE_ABOUT_BUFFER
	 buf.active,
#endif//TRY_MODE_ABOUT_BUFFER
	 buf.freeNum,
#endif//USE_SMID_TO_GET_BUFID
	 buf.freeLst, buf.buffer, buf.bufSize, buf.fail
#   if  !defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
	 , cycles
#endif//!defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
#ifdef  COUNT_INTERACTIONS
	 , treeInfo.Nj, treeInfo.Nbuf
#endif//COUNT_INTERACTIONS
	 );
      //-------------------------------------------------------------------
      *sidx ^= 1;
      //-------------------------------------------------------------------
    }/* if( blck.x <= MAX_BLOCKS_PER_GRID ){ */
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
	   &(tree.more[jhead]), &(tree.jpos[jhead]), &(tree.mj[jhead]),
#endif//SERIALIZED_EXECUTION
#ifdef  DPADD_FOR_ACC
	   pi.tmp,
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	   (jnode *)pi.res,
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
#ifndef USE_SMID_TO_GET_BUFID
#ifndef TRY_MODE_ABOUT_BUFFER
	   buf.active,
#endif//TRY_MODE_ABOUT_BUFFER
	   buf.freeNum,
#endif//USE_SMID_TO_GET_BUFID
	   buf.freeLst, buf.buffer, buf.bufSize, buf.fail
#   if  !defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
	   , cycles
#endif//!defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
#ifdef  COUNT_INTERACTIONS
	   , treeInfo.Nj, treeInfo.Nbuf
#endif//COUNT_INTERACTIONS
	   );
	//-----------------------------------------------------------------
	hidx += Nsub;
	Nrem -= Nblck;
	//-----------------------------------------------------------------
	*sidx ^= 1;
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
#   if  defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
  }/* if( grpNum != 0 ){ */
#endif//defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
  //-----------------------------------------------------------------------
#   if  defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
  checkCudaErrors(cudaEventRecord(finEvent[*Nwalk], 0));
  *Nwalk += 1;
#endif//defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO))
  //-----------------------------------------------------------------------
#ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
  checkCudaErrors(cudaDeviceSynchronize());
  chkMPIerr(MPI_Barrier(MPI_COMM_WORLD));
  fflush(NULL);
  printf("calcAcc finish\n");
  fflush(NULL);
  chkMPIerr(MPI_Barrier(MPI_COMM_WORLD));
#endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
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
#ifdef  USE_CUDA_EVENT
 , cudaEvent_t *iniCalcAcc, cudaEvent_t *finCalcAcc
#else///USE_CUDA_EVENT
 , unsigned long long int *cycles_hst, unsigned long long int *cycles_dev
#endif//USE_CUDA_EVENT
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifndef SERIALIZED_EXECUTION
 , measuredTime *measured, const int pjNum
#ifdef  LET_COMMUNICATION_VIA_HOST
 , const soaTreeNode tree_hst
#endif//LET_COMMUNICATION_VIA_HOST
 , const int Nlet, domainInfo *let, const int Nstream_let, cudaStream_t stream_let[], MPIcfg_tree mpi
#ifdef  MONITOR_LETGEN_TIME
#ifdef  USE_CUDA_EVENT
 , cudaEvent_t *iniMakeLET, cudaEvent_t *finMakeLET
#else///USE_CUDA_EVENT
 , unsigned long long int *cycles_let_hst, unsigned long long int *cycles_let_dev
#endif//USE_CUDA_EVENT
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
  fprintf(stdout, "rank %d: grpNum = %d\n", mpi.rank, grpNum);
  fflush(stdout);
#endif
  //-----------------------------------------------------------------------
#if 0
  jparticle *pj_hst;
  mycudaMallocHost((void **)&pj_hst, pjNum * sizeof(jparticle));
  checkCudaErrors(cudaMemcpy(pj_hst, tree.jpos, pjNum * sizeof(jparticle), cudaMemcpyDeviceToHost));
  jmass *mj_hst;
  mycudaMallocHost((void **)&mj_hst, pjNum * sizeof(jmass));
  checkCudaErrors(cudaMemcpy(mj_hst, tree.mj, pjNum * sizeof(jmass), cudaMemcpyDeviceToHost));

  FILE *fp;
  char filename[128];
  sprintf(filename, "%s/%s.%d_%d.txt", "dat", "jpos", mpi.rank, mpi.size);
  fp = fopen(filename, "w");
  for(int ii = 0; ii < pjNum; ii++)
    fprintf(fp, "%e\t%e\t%e\t%e\n", pj_hst[ii].x, pj_hst[ii].y, pj_hst[ii].z, mj_hst[ii]);

  fclose(fp);
  mycudaFreeHost(pj_hst);
  mycudaFreeHost(mj_hst);
  MPI_Finalize();
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
#if 0
  fprintf(stdout, "rank %d: grpNum = %d with Nsm is %d (%d loops); amin = %e, pos = (%e, %e, %e), r = %e; bufSize is %d\n",
	  mpi.rank, grpNum, devProp.numSM, BLOCKSIZE(blck.x, NBLOCKS_PER_SM * devProp.numSM),
	  let[0].amin, let[0].icom.x, let[0].icom.y, let[0].icom.z, SQRT(let[0].icom.m),
	  buf.bufSize);
  fflush(stdout);
#endif
  //-----------------------------------------------------------------------
  /* initialize measurement counters */
#ifdef  USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  int Nwalk = 0;
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  int Nmake = 0;
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#else///USE_CUDA_EVENT
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
  *cycles_hst = 0;
  checkCudaErrors(cudaMemcpy(cycles_dev, cycles_hst, sizeof(unsigned long long int), cudaMemcpyHostToDevice));
#endif//!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#   if  !defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
  *cycles_let_hst = 0;
  checkCudaErrors(cudaMemcpy(cycles_let_dev, cycles_let_hst, sizeof(unsigned long long int), cudaMemcpyHostToDevice));
#endif//!defined(SERIALIZED_EXECUTION) && defined(MONITOR_LETGEN_TIME)
#endif//USE_CUDA_EVENT
  //-----------------------------------------------------------------------
#ifdef  USE_MEASURED_CLOCK_FREQ
  uint clockWalk;/* in MHz */
#endif//USE_MEASURED_CLOCK_FREQ
  //-----------------------------------------------------------------------
#   if  defined(SERIALIZED_EXECUTION) || defined(EXEC_BENCHMARK)
  static struct timeval start;
  checkCudaErrors(cudaDeviceSynchronize());
  gettimeofday(&start, NULL);
#endif//defined(SERIALIZED_EXECUTION) || defined(EXEC_BENCHMARK)
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
#ifdef  DPADD_FOR_ACC
	 , pi.tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 , pi.res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 );
#else///BLOCK_TIME_STEP
    initAcc_kernel<<<Nrem, NTHREADS>>>
      (pi.acc
#ifdef  GADGET_MAC
       , pi.acc_old
#endif//GADGET_MAC
#ifdef  DPADD_FOR_ACC
       , pi.tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
       , pi.res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
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
#ifdef  DPADD_FOR_ACC
	 , pi.tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 , pi.res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 );
#else///BLOCK_TIME_STEP
      int Nsub = Nblck * NTHREADS;
      initAcc_kernel<<<Nblck, NTHREADS>>>
	(&pi.acc[hidx]
#ifdef  GADGET_MAC
	 , &pi.acc_old[hidx]
#endif//GADGET_MAC
#ifdef  DPADD_FOR_ACC
	 , &pi.tmp[hidx]
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 , &pi.res[hidx]
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
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
#if 0
  checkCudaErrors(cudaDeviceSynchronize());
#endif
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
			  , grpNum, 0
#endif//SERIALIZED_EXECUTION
#   if  !defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)
#ifdef  USE_CUDA_EVENT
			  , &Nwalk, iniCalcAcc, finCalcAcc
#else///USE_CUDA_EVENT
			  , cycles_dev
#endif//USE_CUDA_EVENT
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
#ifdef  USE_MEASURED_CLOCK_FREQ
      /* measure clock frequency as a reference value */
      nvmlDeviceGetClock(deviceHandler, NVML_CLOCK_SM, NVML_CLOCK_ID_CURRENT, &clockWalk);
#endif//USE_MEASURED_CLOCK_FREQ
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
#ifndef SERIALIZED_EXECUTION
      //-------------------------------------------------------------------
      /* gravity from j-particles within other process(es) */
      //-------------------------------------------------------------------
      /* rewrite from MPI_Isend/MPI_Irecv to MPI_Put may accelerate the simulation */
      //-------------------------------------------------------------------
      int idxProcs = 0;
      int remProcs = Nlet - 1;
#ifdef  DOUBLE_BUFFER_FOR_LET
      static int headLETsend[2], headLETrecv[2], sizeLETbuf[2], sizeLETsend[2], sizeLETrecv[2];
      /* 1st half */
      headLETsend[0] = ALIGN_BUF_FOR_LET(pjNum);
      sizeLETbuf [0] = ((int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - headLETsend[0]) >> 1;
      headLETrecv[0] = ALIGN_BUF_FOR_LET(headLETsend[0] + (sizeLETbuf[0] >> 1));
      sizeLETsend[0] = headLETrecv[0] - headLETsend[0];
      sizeLETrecv[0] = sizeLETbuf [0] - sizeLETsend[0];
      /* 2nd half */
      headLETsend[1] = ALIGN_BUF_FOR_LET(headLETrecv[0] + sizeLETrecv[0]);
      sizeLETbuf [1] = (int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - headLETsend[1];
      headLETrecv[1] = ALIGN_BUF_FOR_LET(headLETsend[1] + (sizeLETbuf[1] >> 1));
      sizeLETsend[1] = headLETrecv[1] - headLETsend[1];
      sizeLETrecv[1] = sizeLETbuf [1] - sizeLETsend[1];
#else///DOUBLE_BUFFER_FOR_LET
      int LETsteps = 0;
      const int headLETsend = ALIGN_BUF_FOR_LET(pjNum);
      const int  remLETbuf  = (int)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE) - headLETsend;
      const int  remLETsend = ALIGN_BUF_FOR_LET(remLETbuf >> 1);
      const int  remLETrecv = remLETbuf - remLETsend;
      const int headLETrecv = headLETsend + remLETsend;
#endif//DOUBLE_BUFFER_FOR_LET
#   if  defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
      int prevLETstreams = 0;
#endif//defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
      while( true ){
	//-----------------------------------------------------------------
	/* get maximum number of processes which possible to communicate by limitation of memory capacity */
	//-----------------------------------------------------------------
#ifdef  DOUBLE_BUFFER_FOR_LET
	int remSend = sizeLETsend[sidx];
	int remRecv = sizeLETrecv[sidx];
#else///DOUBLE_BUFFER_FOR_LET
	int remSend = remLETsend;
	int remRecv = remLETrecv;
#endif//DOUBLE_BUFFER_FOR_LET
	int numProcs = remProcs;
	for(int ii = 0; ii < remProcs; ii++){
	  //---------------------------------------------------------------
	  remSend -= let[idxProcs + ii].maxSend;
	  remRecv -= let[idxProcs + ii].maxRecv;
	  //---------------------------------------------------------------
	  if( (remSend < 0) || (remRecv < 0) ){
	    numProcs = ii - 1;
	    break;
	  }/* if( (remSend < 0) || (remRecv < 0) ){ */
	  //---------------------------------------------------------------
	}/* for(int ii = 0; ii < remProcs; ii++){ */
	//-----------------------------------------------------------------
	if( (numProcs < 1) && (mpi.size > 1) ){
	  __KILL__(stderr, "ERROR: numProcs is %d, due to lack of sizeLETsend(%d) or sizeLETrecv(%d) while 0-th target requires numSend(%d) and numRecv(%d) @ rank %d.\n\tIncrease EXTEND_NUM_TREE_NODE(%f) defined in src/tree/let.h and/or TREE_SAFETY_VAL(%f) defined in src/tree/make.h.\n", numProcs,
#ifdef  DOUBLE_BUFFER_FOR_LET
		   sizeLETsend[sidx], sizeLETrecv[sidx],
#else///DOUBLE_BUFFER_FOR_LET
		   remLETsend, remLETrecv,
#endif//DOUBLE_BUFFER_FOR_LET
		   let[idxProcs].maxSend, let[idxProcs].maxRecv, mpi.rank, EXTEND_NUM_TREE_NODE, TREE_SAFETY_VAL);
	}/* if( (numProcs < 1) && (mpi.size > 1) ){ */
	chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, &numProcs, 1, MPI_INT, MPI_MIN, mpi.comm));
	//-----------------------------------------------------------------
	/* set send buffer for LET on device */
#ifdef  DOUBLE_BUFFER_FOR_LET
	let[idxProcs].headSend = headLETsend[sidx];
#else///DOUBLE_BUFFER_FOR_LET
	let[idxProcs].headSend = headLETsend;
#endif//DOUBLE_BUFFER_FOR_LET
	for(int ii = 0; ii < numProcs - 1; ii++)
	  let[idxProcs + ii + 1].headSend = let[idxProcs + ii].headSend + ALIGN_BUF_FOR_LET(let[idxProcs + ii].maxSend);
	//-----------------------------------------------------------------


	//-----------------------------------------------------------------
	/* FUTURE UPDATE: divide below procedure to overlap communication and calculation */
	//-----------------------------------------------------------------
	/* generate numProcs LET(s) */
	for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){
	  //---------------------------------------------------------------
	  const int streamIdxLET = ii % Nstream_let;
	  //---------------------------------------------------------------
	  callGenLET(stream_let[streamIdxLET], &let[ii], mpi, tree, buf
#ifdef  MONITOR_LETGEN_TIME
#ifdef  USE_CUDA_EVENT
		     , iniMakeLET[Nmake], finMakeLET[Nmake]
#else///USE_CUDA_EVENT
		     , cycles_let_dev
#endif//USE_CUDA_EVENT
#endif//MONITOR_LETGEN_TIME
		     );
#   if  defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
	  Nmake++;
#endif//defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME)
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaMemcpyAsync(let[ii].numSend_hst, let[ii].numSend_dev, sizeof(int), cudaMemcpyDeviceToHost, stream_let[streamIdxLET]));
	  //---------------------------------------------------------------
#if 0
	  printf("rank %d: grpNum = %d, Nlet = %d\n", mpi.rank, grpNum, *(let[ii].numSend_hst));
	  fflush(stdout);
#endif
	  //---------------------------------------------------------------
	}/* for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){ */
	//-----------------------------------------------------------------
	/* share # of LET nodes */
	for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){
	  //---------------------------------------------------------------
	  const int streamIdxLET = ii % Nstream_let;
	  //---------------------------------------------------------------
	  checkCudaErrors(cudaStreamSynchronize(stream_let[streamIdxLET]));
	  let[ii].numSend = *(let[ii].numSend_hst);
	  if( let[ii].numSend > let[ii].maxSend ){
	    __KILL__(stderr, "ERROR: predicted size of send buffer (%d) is not sufficient for true size of that (%d) @ rank %d for rand %d.\n\tsuggestion: consider increasing \"LETSIZE_REDUCE_FACTOR\" defined in src/tree/let.h (current value is %f) to at least %f.\n", let[ii].maxSend, let[ii].numSend, mpi.rank, let[ii].rank, LETSIZE_REDUCE_FACTOR, LETSIZE_REDUCE_FACTOR * (float)let[ii].numSend / (float)let[ii].maxSend);
	  }/* if( let[ii].numSend > let[ii].maxSend ){ */
	  __NOTE__("numSend = %d, numFull = %d @ rank %d\n", let[ii].numSend, let[ii].numFull, mpi.rank);
	  //---------------------------------------------------------------
	  /* send # of LET nodes */
	  chkMPIerr(MPI_Isend(&(let[ii].numSend), 1, MPI_INT, let[ii].rank,  mpi.rank, mpi.comm, &(let[ii].reqSendInfo)));
	  let[ii].numRecvGuess = let[ii].maxRecv;
	  chkMPIerr(MPI_Irecv(&(let[ii].numRecv), 1, MPI_INT, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvInfo)));
#ifdef  DBG_LETGEN_ON_GPU
	  fprintf(stdout, "rank = %d: ii = %d: LET(target is %d): numSend = %d out of %d nodes\n", mpi.rank, ii, let[ii].rank, let[ii].numSend, pjNum);
	  fflush(stdout);
#endif//DBG_LETGEN_ON_GPU
	  //---------------------------------------------------------------
	  /* copy LET nodes from device to host */
#ifdef  LET_COMMUNICATION_VIA_HOST
	  checkCudaErrors(cudaMemcpyAsync(&(tree_hst.more[let[ii].headSend]), &(tree.more[let[ii].headSend]), sizeof(     uint) * let[ii].numSend, cudaMemcpyDeviceToHost, stream_let[streamIdxLET]));
	  checkCudaErrors(cudaMemcpyAsync(&(tree_hst.jpos[let[ii].headSend]), &(tree.jpos[let[ii].headSend]), sizeof(jparticle) * let[ii].numSend, cudaMemcpyDeviceToHost, stream_let[streamIdxLET]));
	  checkCudaErrors(cudaMemcpyAsync(&(tree_hst.mj  [let[ii].headSend]), &(tree.mj  [let[ii].headSend]), sizeof(    jmass) * let[ii].numSend, cudaMemcpyDeviceToHost, stream_let[streamIdxLET]));
#endif//LET_COMMUNICATION_VIA_HOST
	  //---------------------------------------------------------------
	}/* for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){ */
	//-----------------------------------------------------------------
	/* send numProcs LET(s) to other process(es) */
	for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){
	  //---------------------------------------------------------------
	  const int streamIdxLET = ii % Nstream_let;
	  //---------------------------------------------------------------
	  /* send LET nodes using MPI_Isend */
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
	}/* for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){ */
	//-----------------------------------------------------------------
	/* chkMPIerr(MPI_Barrier(mpi.comm)); */
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* receive # of LET nodes */
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
	  //---------------------------------------------------------------
	}/* for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){ */
	//-----------------------------------------------------------------
	/* chkMPIerr(MPI_Barrier(mpi.comm)); */
	//-----------------------------------------------------------------
	/* set receive buffer for LET on device */
#ifdef  DOUBLE_BUFFER_FOR_LET
	let[idxProcs].headRecv = headLETrecv[sidx];
#else///DOUBLE_BUFFER_FOR_LET
	let[idxProcs].headRecv = headLETrecv;
#endif//DOUBLE_BUFFER_FOR_LET
	int numRecv = 0;
	for(int ii = 0; ii < numProcs - 1; ii++){
	  const int numRecvBuf = ALIGN_BUF_FOR_LET(let[idxProcs + ii].numRecv);
	  let[idxProcs + ii + 1].headRecv = numRecvBuf + let[idxProcs + ii].headRecv;
	  numRecv                        += numRecvBuf;
	}/* for(int ii = 0; ii < numProcs - 1; ii++){ */
	numRecv += ALIGN_BUF_FOR_LET(let[idxProcs + numProcs - 1].numRecv);
	//-----------------------------------------------------------------
#ifdef  DOUBLE_BUFFER_FOR_LET
	if( numRecv > sizeLETrecv[sidx] )
#else///DOUBLE_BUFFER_FOR_LET
	if( numRecv > remLETrecv )
#endif//DOUBLE_BUFFER_FOR_LET
	  {
	    __KILL__(stderr, "ERROR: lack of remLETrecv(%d) to store numRecv(%d) LET nodes.\n\tIncrease EXTEND_NUM_TREE_NODE(%f) defined in src/tree/let.h and/or TREE_SAFETY_VAL(%f) defined in src/tree/make.h.\n",
#ifdef  DOUBLE_BUFFER_FOR_LET
		     sizeLETrecv[sidx],
#else///DOUBLE_BUFFER_FOR_LET
		     remLETrecv,
#endif//DOUBLE_BUFFER_FOR_LET
		     numRecv, EXTEND_NUM_TREE_NODE, TREE_SAFETY_VAL);
	  }
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* receive LET nodes */
	//-----------------------------------------------------------------
	/* before receiving LET nodes, gravity calculation using LET nodes stored in the receive buffer in the previous loop must be finished */
#ifndef DOUBLE_BUFFER_FOR_LET
	if( LETsteps > 0 ){
#   if  defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
	  if( grpNum != 0 ){
#endif//defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
	    checkCudaErrors(cudaStreamSynchronize(sinfo->stream[sidx ^ 1]));
	    if( blck.x > MAX_BLOCKS_PER_GRID )
	      checkCudaErrors(cudaStreamSynchronize(sinfo->stream[sidx ]));
#   if  defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
	  }/* if( grpNum != 0 ){ */
#endif//defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
	}/* if( LETsteps > 0 ){ */
#endif//DOUBLE_BUFFER_FOR_LET
	//-----------------------------------------------------------------
	for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){
	  //---------------------------------------------------------------
#ifdef  LET_COMMUNICATION_VIA_HOST
	  chkMPIerr(MPI_Irecv(&(tree_hst.more[let[ii].headRecv]), let[ii].numRecv, mpi.more, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvMore)));
	  chkMPIerr(MPI_Irecv(&(tree_hst.jpos[let[ii].headRecv]), let[ii].numRecv, mpi.jpos, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvJpos)));
	  chkMPIerr(MPI_Irecv(&(tree_hst.mj  [let[ii].headRecv]), let[ii].numRecv, mpi.mass, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvMass)));
#else///LET_COMMUNICATION_VIA_HOST
	  chkMPIerr(MPI_Irecv(&(tree    .more[let[ii].headRecv]), let[ii].numRecv, mpi.more, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvMore)));
	  chkMPIerr(MPI_Irecv(&(tree    .jpos[let[ii].headRecv]), let[ii].numRecv, mpi.jpos, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvJpos)));
	  chkMPIerr(MPI_Irecv(&(tree    .mj  [let[ii].headRecv]), let[ii].numRecv, mpi.mass, let[ii].rank, let[ii].rank, mpi.comm, &(let[ii].reqRecvMass)));
#endif//LET_COMMUNICATION_VIA_HOST
	  //---------------------------------------------------------------
	}/* for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){ */
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* receive numProcs LET(s) and calculate gravity from them */
	//-----------------------------------------------------------------
	for(int ii = idxProcs; ii < idxProcs + numProcs; ii++){
	  //---------------------------------------------------------------
	  /* copy LET nodes from host to device */
	  //---------------------------------------------------------------
	  MPI_Status statusMore;	  chkMPIerr(MPI_Wait(&(let[ii].reqRecvJpos), &statusMore));
#ifdef  LET_COMMUNICATION_VIA_HOST
	  checkCudaErrors(cudaMemcpyAsync(&(tree.jpos[let[ii].headRecv]), &(tree_hst.jpos[let[ii].headRecv]), sizeof(jparticle) * let[ii].numRecv, cudaMemcpyHostToDevice, sinfo->stream[sidx]));
#endif//LET_COMMUNICATION_VIA_HOST
	  MPI_Status statusMass;	  chkMPIerr(MPI_Wait(&(let[ii].reqRecvMass), &statusMass));
#ifdef  LET_COMMUNICATION_VIA_HOST
	  checkCudaErrors(cudaMemcpyAsync(&(tree.mj  [let[ii].headRecv]), &(tree_hst.mj  [let[ii].headRecv]), sizeof(    jmass) * let[ii].numRecv, cudaMemcpyHostToDevice, sinfo->stream[sidx]));
#endif//LET_COMMUNICATION_VIA_HOST
	  MPI_Status statusJpos;	  chkMPIerr(MPI_Wait(&(let[ii].reqRecvMore), &statusJpos));
#ifdef  LET_COMMUNICATION_VIA_HOST
	  checkCudaErrors(cudaMemcpyAsync(&(tree.more[let[ii].headRecv]), &(tree_hst.more[let[ii].headRecv]), sizeof(     uint) * let[ii].numRecv, cudaMemcpyHostToDevice, sinfo->stream[sidx]));
#endif//LET_COMMUNICATION_VIA_HOST
	  //---------------------------------------------------------------
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
#ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
	  printf("# %d-th LET\n", ii);
	  checkCudaErrors(cudaDeviceSynchronize());
	  chkMPIerr(MPI_Barrier(MPI_COMM_WORLD));
	  fflush(NULL);
#endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
	  callCalcGravityFunc(blck, thrd, sinfo, &sidx, laneInfo, pi, 0, tree, grpNum, let[ii].headRecv
#ifdef  USE_CUDA_EVENT
			      , &Nwalk, iniCalcAcc, finCalcAcc
#else///USE_CUDA_EVENT
			      , cycles_dev
#endif//USE_CUDA_EVENT
			      , buf
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
#ifndef DOUBLE_BUFFER_FOR_LET
	LETsteps++;
#endif//DOUBLE_BUFFER_FOR_LET
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
#endif//SERIALIZED_EXECUTION
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
      sinfo->idx = sidx;
#if 0
      checkCudaErrors(cudaDeviceSynchronize());
      /* checkCudaErrors(cudaStreamSynchronize(sinfo->stream[sidx    ])); */
      /* checkCudaErrors(cudaStreamSynchronize(sinfo->stream[sidx ^ 1])); */
#endif
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
	__KILL__(stderr, "ERROR: bufUsed exceeds bufSize of %d at least %d times.\nPLEASE re-simulate after decreasing NUM_BODY_MAX(%d) or GLOBAL_MEMORY_SYSBUF(%zu) defined in src/misc/structure.h or TREE_SAFETY_VAL(%f) defined in src/tree/make.h.\n", buf.bufSize, fail_hst, NUM_BODY_MAX, (size_t)GLOBAL_MEMORY_SYSBUF, TREE_SAFETY_VAL);
#else///SERIALIZED_EXECUTION
	__KILL__(stderr, "ERROR: bufUsed exceeds bufSize of %d at least %d times.\nPLEASE re-simulate after decreasing NUM_BODY_MAX(%d) or GLOBAL_MEMORY_SYSBUF(%zu) defined in src/misc/structure.h or TREE_SAFETY_VAL(%f) defined in src/tree/make.h, or EXTEND_NUM_TREE_NODE(%f) defined in src/tree/let.h.\n", buf.bufSize, fail_hst, NUM_BODY_MAX, (size_t)GLOBAL_MEMORY_SYSBUF, TREE_SAFETY_VAL, EXTEND_NUM_TREE_NODE);
#endif//SERIALIZED_EXECUTION
      }/* if( fail_hst != 0 ){ */
      //-------------------------------------------------------------------
/* #   if  !defined(SERIALIZED_EXECUTION) && defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME) */
/*       for(int ii = 0; ii < prevLETstreams; ii++){ */
/* 	checkCudaErrors(cudaEventSynchronize(finMakeLET[ii])); */
/* 	checkCudaErrors(cudaEventElapsedTime(&makeLET_ms, iniCalcAcc[ii], finCalcAcc[ii])); */
/* 	makeLET += (double)makeLET_ms * 1.0e-3; */
/*       }/\* for(int jj = 0; jj < prevLETstreams; jj++){ *\/ */
/* #endif//!defined(SERIALIZED_EXECUTION) && defined(USE_CUDA_EVENT) && defined(MONITOR_LETGEN_TIME) */
      //-------------------------------------------------------------------
/* #   if  defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)) */
/*       checkCudaErrors(cudaEventSynchronize(finCalcAcc[sidx ^ 1])); */
/*       checkCudaErrors(cudaEventElapsedTime(&calcAcc_ms, iniCalcAcc[sidx ^ 1], finCalcAcc[sidx ^ 1])); */
/*       calcAcc += (double)calcAcc_ms * 1.0e-3; */
/* #endif//defined(USE_CUDA_EVENT) && (!defined(SERIALIZED_EXECUTION) || defined(PRINT_PSEUDO_PARTICLE_INFO)) */
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
      trimAcc_kernel<<<Nrem, thrd>>>
	(pi.acc, pi.pos, BLOCKSIZE(grpNum, NGROUPS) * NGROUPS, laneInfo
#ifdef  DPADD_FOR_ACC
	 , pi.tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 , pi.res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 );
#else///BLOCK_TIME_STEP
    trimAcc_kernel<<<Nrem, NTHREADS>>>
      (pi.acc, pi.pos
#ifdef  DPADD_FOR_ACC
       , pi.tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
       , pi.res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
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
      trimAcc_kernel<<<Nblck, thrd.x>>>
	(pi.acc, pi.pos, BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx]
#ifdef  DPADD_FOR_ACC
	 , pi.tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 , pi.res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 );
#else///BLOCK_TIME_STEP
      int Nsub = Nblck * NTHREADS;
      trimAcc_kernel<<<Nblck, NTHREADS>>>
	(&pi.acc[hidx], &pi.pos[hidx]
#ifdef  DPADD_FOR_ACC
	 , &pi.tmp[hidx]
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
	 , &pi.res[hidx]
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
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
  getLastCudaError("trimAcc_kernel");
  //-----------------------------------------------------------------------
#if 0
  checkCudaErrors(cudaDeviceSynchronize());
  MPI_Finalize();
  exit(0);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#   if  defined(SERIALIZED_EXECUTION) || defined(EXEC_BENCHMARK)
  static struct timeval finish;
  checkCudaErrors(cudaDeviceSynchronize());
  gettimeofday(&finish, NULL);
  *time = calcElapsedTimeInSec(start, finish);
#ifdef  EXEC_BENCHMARK
  elapsed->calcGravity_dev += *time;
#endif//EXEC_BENCHMARK
#endif//defined(SERIALIZED_EXECUTION) || defined(EXEC_BENCHMARK)
  //-----------------------------------------------------------------------
  /* evaluate GPU time */
#ifndef SERIALIZED_EXECUTION
  checkCudaErrors(cudaDeviceSynchronize());
#ifdef  USE_CUDA_EVENT
  double calcAcc = 0.0;
  for(int ii = 0; ii < Nwalk; ii++){
    float tmp_ms;
    checkCudaErrors(cudaEventElapsedTime(&tmp_ms, iniCalcAcc[ii], finCalcAcc[ii]));
    calcAcc += (double)tmp_ms;
  }/* for(int ii = 0; ii < Nwalk; ii++){ */
  calcAcc *= 1.0e-3;
#else///USE_CUDA_EVENT
  checkCudaErrors(cudaMemcpy(cycles_hst, cycles_dev, sizeof(unsigned long long int), cudaMemcpyDeviceToHost));
  /* # of launched blocks for tree traversal = # of blocks per kernel function * (# of local tree + # of LETs) = # of blocks per kernel function * # of GPUs */
#ifdef  USE_MEASURED_CLOCK_FREQ
  double devClock = (double)clockWalk * 1.0e+6;
#if 0
  printf("%e Hz on rank %d\n", devClock, mpi.rank);
  MPI_Finalize();
  exit(0);
#endif
#endif//USE_MEASURED_CLOCK_FREQ
#ifdef  USE_GPU_BASE_CLOCK_FREQ
  const double devClock = devProp.coreClk * 1.0e+9;
#endif//USE_GPU_BASE_CLOCK_FREQ
  /* const double calcAcc = ((double)(*cycles_hst) / (devClock * (double)(blck.x * mpi.size))) * (double)BLOCKSIZE(blck.x * mpi.size, devProp.numSM * NBLOCKS_PER_SM); */
  const double calcAcc = ((double)(*cycles_hst) / (devClock * (double)(blck.x * mpi.size))) * (double)BLOCKSIZE(blck.x * mpi.size, devProp.numSM);
#endif//USE_CUDA_EVENT
  measured->sum_excg    += calcAcc;
  measured->sum_rebuild += calcAcc;
  *time   = calcAcc;
#ifdef  MONITOR_LETGEN_TIME
#ifdef  USE_CUDA_EVENT
  double makeLET = 0.0;
  for(int ii = 0; ii < Nmake; ii++){
    float tmp_ms;
    checkCudaErrors(cudaEventElapsedTime(&tmp_ms, iniMakeLET[ii], finMakeLET[ii]));
    makeLET += (double)tmp_ms;
  }/* for(int ii = 0; ii < Nwalk; ii++){ */
  makeLET *= 1.0e-3;
#else///USE_CUDA_EVENT
  checkCudaErrors(cudaMemcpy(cycles_let_hst, cycles_let_dev, sizeof(unsigned long long int), cudaMemcpyDeviceToHost));
  /* # of launched blocks for LET generator = # of LETs = # of GPUs - 1 = # of MPI processes - 1 */
  /* const double makeLET = ((double)(*cycles_let_hst) / (devClock * (double)(mpi.size - 1))) * BLOCKSIZE(mpi.size - 1, devProp.numSM * NBLOCKS_PER_SM); */
  const double makeLET = ((double)(*cycles_let_hst) / (devClock * (double)(mpi.size - 1))) * BLOCKSIZE(mpi.size - 1, devProp.numSM);
#endif//USE_CUDA_EVENT
  measured->sum_excg    += makeLET;
  measured->sum_rebuild += makeLET;
#if 0
  static struct timeval finish;
  checkCudaErrors(cudaDeviceSynchronize());
  gettimeofday(&finish, NULL);
  fprintf(stdout, "rank %d: %e + %e | %e\n", mpi.rank, calcAcc, makeLET, calcElapsedTimeInSec(start, finish));
  fflush(stdout);
#endif
#endif//MONITOR_LETGEN_TIME
#endif//SERIALIZED_EXECUTION
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
#ifdef  DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
  MPI_Finalize();
  exit(0);
#endif//DEBUG_PRINT_FOR_PARTICLE_ACCELERATION
  //-----------------------------------------------------------------------
#if 0
  MPIinfo mpi_tmp;
  mpi_tmp.rank = mpi.rank;
  mpi_tmp.size = mpi.size;
  mpi_tmp.comm = mpi.comm;
  printTreeNode_dev(pjNum, tree, mpi_tmp);
#endif
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
