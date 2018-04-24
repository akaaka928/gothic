/**
 * @file exchange_dev.cu
 *
 * @brief Source code for domain decomposition using GPUs with MPI
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/04/13 (Fri)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#define OMIT_INPLACE_IN_SCATTERV

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <helper_cuda.h>

#include <thrust/sort.h>
#include <thrust/device_ptr.h>

#include "macro.h"
#include "mpilib.h"
#include "cudalib.h"
#include "timer.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/tune.h"
#include "../misc/brent.h"

#include "../util/gsync_dev.cu"

#include "../time/adv_dev.h"
#include "../sort/peano_dev.h"

#include "mpicfg.h"
#include "exchange.h"
#include "exchange_dev.h"

#ifdef  USE_RECTANGULAR_BOX_FOR_LET
#include "../tree/make.h"/**< to read laneinfo */
#endif//USE_RECTANGULAR_BOX_FOR_LET

#ifndef SERIALIZED_EXECUTION
#define MPI_TAG_BODY(rank, size) ((rank) + ((size) * 10))
#define MPI_TAG_IPOS(rank, size) ((rank) + ((size) * 11))
#define MPI_TAG_IIDX(rank, size) ((rank) + ((size) * 12))
#define MPI_TAG_IVEL(rank, size) ((rank) + ((size) * 13))
#define MPI_TAG_TIME(rank, size) ((rank) + ((size) * 14))
#define MPI_TAG_VELX(rank, size) ((rank) + ((size) * 15))
#define MPI_TAG_VELY(rank, size) ((rank) + ((size) * 16))
#define MPI_TAG_VELZ(rank, size) ((rank) + ((size) * 17))
#define MPI_TAG_IACC(rank, size) ((rank) + ((size) * 18))
#define MPI_TAG_IEXT(rank, size) ((rank) + ((size) * 19))
#endif//SERIALIZED_EXECUTION


/** tentative treatment for GTX TITAN X: (# of SM = 24) * (# of blocks per SM = 8) = 192 exceeds NTHREADS_ASSIGN */
#ifndef RUN_ON_PC
/* #   if  GPUGEN == 52 */
#   if  NTHREADS_ASSIGN < 256
#undef  NTHREADS_ASSIGN
#define NTHREADS_ASSIGN  (256)
#endif//NTHREADS_ASSIGN < 256
/* #endif//GPUGEN == 52 */
#endif//RUN_ON_PC


#ifndef USE_OCCUPANCY_CALCULATOR
#   if  GPUGEN >= 60
/** capacity of shared memory is 64KiB per SM on newer GPUs */
/** int4 num_sm[NTHREADS_ASSIGN] corresponds 16 * NTHREADS_ASSIGN bytes */
#define NBLOCKS_PER_SM_ASSIGN (4096 / NTHREADS_ASSIGN)
#else///GPUGEN >= 60
/** in L1 cache preferred configuration, capacity of shared memory is 16KiB per SM on older GPUs */
/** int4 num_sm[NTHREADS_ASSIGN] corresponds 16 * NTHREADS_ASSIGN bytes */
#define NBLOCKS_PER_SM_ASSIGN (1024 / NTHREADS_ASSIGN)
#endif//GPUGEN >= 60

#define REGISTERS_PER_THREAD_ASSIGN (32)

/** limitation from number of registers */
#   if  NBLOCKS_PER_SM_ASSIGN > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_ASSIGN * NTHREADS_ASSIGN))
#undef  NBLOCKS_PER_SM_ASSIGN
#define NBLOCKS_PER_SM_ASSIGN   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_ASSIGN * NTHREADS_ASSIGN))
#endif//NBLOCKS_PER_SM_ASSIGN > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_ASSIGN * NTHREADS_ASSIGN))

/** maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_ASSIGN > (MAX_THREADS_PER_SM / NTHREADS_ASSIGN)
#undef  NBLOCKS_PER_SM_ASSIGN
#define NBLOCKS_PER_SM_ASSIGN   (MAX_THREADS_PER_SM / NTHREADS_ASSIGN)
#endif//NBLOCKS_PER_SM_ASSIGN > (MAX_THREADS_PER_SM / NTHREADS_ASSIGN)

/** maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_ASSIGN > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_ASSIGN
#define NBLOCKS_PER_SM_ASSIGN   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_ASSIGN > MAX_BLOCKS_PER_SM

/** maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_ASSIGN > ((MAX_WARPS_PER_SM * 32) / NTHREADS_ASSIGN)
#undef  NBLOCKS_PER_SM_ASSIGN
#define NBLOCKS_PER_SM_ASSIGN   ((MAX_WARPS_PER_SM * 32) / NTHREADS_ASSIGN)
#endif//NBLOCKS_PER_SM_ASSIGN > ((MAX_WARPS_PER_SM * 32) / NTHREADS_ASSIGN)

/** # of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_ASSIGN < 1
#undef  NBLOCKS_PER_SM_ASSIGN
#define NBLOCKS_PER_SM_ASSIGN  (1)
#endif//NBLOCKS_PER_SM_ASSIGN < 1
#endif//USE_OCCUPANCY_CALCULATOR


/**
 * @fn pickupSamples_kernel
 *
 * @brief Pick up sample particles on GPU.
 */
__global__ void pickupSamples_kernel(const int num, const int nskip, READ_ONLY position * ipos, float * RESTRICT xi, float * RESTRICT yi, float * RESTRICT zi, const int offset)
{
  const int ii = offset + GLOBALIDX_X1D;
  if( ii < num ){
    const position pi = ipos[ii * nskip];
    xi[ii] = CAST_R2F(pi.x);
    yi[ii] = CAST_R2F(pi.y);
    zi[ii] = CAST_R2F(pi.z);
  }/* if( ii < num ){ */
}


#define NTHREADS_COMPARE_VEC3_INC NTHREADS_BOX
#include "../util/compare_vec3_inc.cu"
/**
 * @fn getBoxSize_kernel
 *
 * @brief Calculate box size to include all N-body particles in the local domain.
 */
__global__ void
#ifndef USE_OCCUPANCY_CALCULATOR
__launch_bounds__(NTHREADS_BOX, NBLOCKS_PER_SM_BOX)
#endif//USE_OCCUPANCY_CALCULATOR
getBoxSize_kernel
     (const int num, READ_ONLY position * RESTRICT ipos, float4 * RESTRICT min_all, float4 * RESTRICT max_all,
      int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int bnum =   GRIDDIM_X1D;

  const int hidx = tidx - (tidx & (warpSize - 1));

  const int ihead = (num *      bidx ) / bnum;
  const int itail = (num * (1 + bidx)) / bnum;


  /** calculate required box size to contain all N-body particles in the local domain */
  float4 min = { FLT_MAX,  FLT_MAX,  FLT_MAX,  FLT_MAX};
  float4 max = {-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX};
  for(int ih = ihead; ih < itail; ih += NTHREADS_BOX){
    const int ii = ih + tidx;
    if( ii < num ){
      const position pi = ipos[ii];
      const float xi = CAST_R2F(pi.x);
      const float yi = CAST_R2F(pi.y);
      const float zi = CAST_R2F(pi.z);

      min.x = fminf(min.x, xi);      max.x = fmaxf(max.x, xi);
      min.y = fminf(min.y, yi);      max.y = fmaxf(max.y, yi);
      min.z = fminf(min.z, zi);      max.z = fmaxf(max.z, zi);
    }/* if( ii < num ){ */
  }/* for(int ih = ihead; ih < itail; ih += NTHREADS_BOX){ */


  /** get the box size within a warp */
  __shared__ float4 smem[NTHREADS_BOX];

  GET_VEC3_MIN_MAX_GRID(&min, &max, smem, tidx, hidx, min_all, max_all, bidx, bnum, gsync0, gsync1);
  if( (bnum == 1) && ((bidx + tidx) == 0) ){
    min_all[0] = min;
    max_all[0] = max;
  }/* if( (bnum == 1) && ((bidx + tidx) == 0) ){ */
}


/**
 * @fn checkBoxSize_dev
 *
 * @brief Error detection before running the kernel
 *
 * @sa getBoxSize_kernel
 */
extern "C"
void checkBoxSize_dev
(
#ifdef  USE_OCCUPANCY_CALCULATOR
 soaPHsort *soa
#else///USE_OCCUPANCY_CALCULATOR
 const deviceProp devProp
#endif//USE_OCCUPANCY_CALCULATOR
)
{
  __NOTE__("%s\n", "start");

#ifdef  USE_OCCUPANCY_CALCULATOR
  checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&soa->numBlocksPerSM_box, getBoxSize_kernel, NTHREADS_BOX, 0));
#else///USE_OCCUPANCY_CALCULATOR
  struct cudaFuncAttributes funcAttr;
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, getBoxSize_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_BOX ){
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in getBoxSize_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_BOX);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_BOX = %d.\n", GPUGEN, GPUVER, NTHREADS_BOX);
    fflush (stderr);
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAC ){ */

  int regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS_BOX);
  if( regLimit > (MAX_REGISTERS_PER_SM / NTHREADS_BOX) )
    regLimit = (MAX_REGISTERS_PER_SM / NTHREADS_BOX);
  int memLimit = SMEM_SIZE_L1_PREF / funcAttr.sharedSizeBytes;
  int Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / NTHREADS_BOX) )    Nblck = MAX_THREADS_PER_SM / NTHREADS_BOX;
  if( Nblck >   MAX_BLOCKS_PER_SM                      )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / NTHREADS_BOX) )    Nblck = ((MAX_WARPS_PER_SM * 32) / NTHREADS_BOX);

  if( Nblck != NBLOCKS_PER_SM_BOX ){
    __KILL__(stderr, "ERROR: # of blocks per SM for getBoxSize_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_BOX defined in src/para/box_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_BOX, funcAttr.numRegs, REGISTERS_PER_THREAD_BOX, GPUGEN, GPUVER);
  }/* if( Nblck != NBLOCKS_PER_SM ){ */

  if( (devProp.numSM * NBLOCKS_PER_SM_BOX) > NTHREADS_BOX ){
    __KILL__(stderr, "ERROR: product (%d) of devProp.numSM(%d) * NBLOCKS_PER_SM_BOX(%d) must be smaller than NTHREADS_BOX(%d) to use shared memory.\n", devProp.numSM * NBLOCKS_PER_SM_BOX, devProp.numSM, NBLOCKS_PER_SM_BOX, NTHREADS_BOX);
  }/* if( (devProp.numSM * NBLOCKS_PER_SM_BOX) > NTHREADS_BOX ){ */
#endif//USE_OCCUPANCY_CALCULATOR

  __NOTE__("%s\n", "end");
}


/**
 * @fn getBoxSize_dev
 *
 * @brief Calculate box size to include all N-body particles in the local domain.
 */
extern "C"
void getBoxSize_dev(const int num, position * RESTRICT ipos, soaPHsort soa, const deviceProp devProp, cudaStream_t stream)
{
  __NOTE__("%s\n", "start");


#ifdef  USE_OCCUPANCY_CALCULATOR
  const int NBLOCKS_PER_SM_BOX = soa.numBlocksPerSM_box;
#endif//USE_OCCUPANCY_CALCULATOR

  getBoxSize_kernel<<<devProp.numSM * NBLOCKS_PER_SM_BOX, NTHREADS_BOX, SMEM_SIZE, stream>>>(num, ipos, soa.min, soa.max, soa.gsync0, soa.gsync1);
  getLastCudaError("getBoxSize_kernel");

  checkCudaErrors(cudaMemcpyAsync(soa.min_hst, soa.min, sizeof(float4), cudaMemcpyDeviceToHost, stream));
  checkCudaErrors(cudaMemcpyAsync(soa.max_hst, soa.max, sizeof(float4), cudaMemcpyDeviceToHost, stream));

  checkCudaErrors(cudaDeviceSynchronize());


  __NOTE__("%s\n", "end");
}


#ifdef  USE_RECTANGULAR_BOX_FOR_LET
extern "C"
void getEnclosingBox_dev(const int grpNum, const int Ngrp, laneinfo *laneInfo_hst, iparticle pi, soaPHsort soa, const deviceProp devProp)
{
  __NOTE__("%s\n", "start");
  __NOTE__("grpNum = %d, Ngrp = %d\n", grpNum, Ngrp);

#ifdef  USE_OCCUPANCY_CALCULATOR
  const int NBLOCKS_PER_SM_BOX = soa.numBlocksPerSM_box;
#endif//USE_OCCUPANCY_CALCULATOR

  if( grpNum > 0 ){
    /* getBoxSize_kernel<<<devProp.numSM * NBLOCKS_PER_SM_BOX, NTHREADS_BOX>>>(laneInfo_hst[(grpNum < Ngrp) ? (grpNum) : (Ngrp - 1)].head, pi.pos, soa.min, soa.max, soa.gsync0, soa.gsync1); */
    getBoxSize_kernel<<<devProp.numSM * NBLOCKS_PER_SM_BOX, NTHREADS_BOX>>>(laneInfo_hst[grpNum - 1].head + laneInfo_hst[grpNum - 1].num, pi.pos, soa.min, soa.max, soa.gsync0, soa.gsync1);
    getLastCudaError("getBoxSize_kernel");

    checkCudaErrors(cudaMemcpy(soa.min_hst, soa.min, sizeof(float4), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(soa.max_hst, soa.max, sizeof(float4), cudaMemcpyDeviceToHost));

    pi.min_hst->x = CAST_F2R(soa.min_hst->x);    pi.min_hst->y = CAST_F2R(soa.min_hst->y);    pi.min_hst->z = CAST_F2R(soa.min_hst->z);
    pi.max_hst->x = CAST_F2R(soa.max_hst->x);    pi.max_hst->y = CAST_F2R(soa.max_hst->y);    pi.max_hst->z = CAST_F2R(soa.max_hst->z);

    pi.icom_hst->x = HALF * (pi.min_hst->x + pi.max_hst->x);
    pi.icom_hst->y = HALF * (pi.min_hst->y + pi.max_hst->y);
    pi.icom_hst->z = HALF * (pi.min_hst->z + pi.max_hst->z);

    const real dx = HALF * (pi.max_hst->x - pi.min_hst->x);
    const real dy = HALF * (pi.max_hst->y - pi.min_hst->y);
    const real dz = HALF * (pi.max_hst->z - pi.min_hst->z);
    const real r2 = 1.0e-30f + dx * dx + dy * dy + dz * dz;
    pi.icom_hst->m = pi.min_hst->m = pi.max_hst->m = r2;
  }/* if( grpNum > 0 ){ */
  else{
#ifdef  SKIP_UNUSED_LET_GENERATION
    pi.icom_hst->m = pi.min_hst->m = pi.max_hst->m = -UNITY;
#else///SKIP_UNUSED_LET_GENERATION
    pi.icom_hst->m = pi.min_hst->m = pi.max_hst->m = ZERO;
#endif//SKIP_UNUSED_LET_GENERATION
  }/* else{ */

  __NOTE__("xl = %e, yl = %e, zl = %e, r2 = %e\n", pi. min_hst->x, pi. min_hst->y, pi. min_hst->z, pi. min_hst->m);
  __NOTE__("x0 = %e, y0 = %e, z0 = %e, r2 = %e\n", pi.icom_hst->x, pi.icom_hst->y, pi.icom_hst->z, pi.icom_hst->m);
  __NOTE__("xu = %e, yu = %e, zu = %e, r2 = %e\n", pi. max_hst->x, pi. max_hst->y, pi. max_hst->z, pi. max_hst->m);

  __NOTE__("%s\n", "end");
}
#endif//USE_RECTANGULAR_BOX_FOR_LET


/**
 * @fn setIndex_kernel
 *
 * @brief Set index for sorting on GPU.
 */
__global__ void setIndex_kernel(const int num, int *idx, const int offset)
{
  const int ii = offset + GLOBALIDX_X1D;
  if( ii < num )
    idx[ii] = ii;
}

/**
 * @fn sortSamplePos_yz_kernel
 *
 * @brief Sort position arrays (y and z) on GPU.
 */
__global__ void sortSamplePos_yz_kernel
(const int num, READ_ONLY int * RESTRICT old, float * RESTRICT dy, READ_ONLY float * RESTRICT sy, float * RESTRICT dz, READ_ONLY float * RESTRICT sz, const int offset)
{
  const int ii = offset + GLOBALIDX_X1D;

  if( ii < num ){
    const int jj = old[ii];
    dy[ii] = sy[jj];
    dz[ii] = sz[jj];
  }/* if( ii < num ){ */
}

/**
 * @fn sort_xpos_dev
 *
 * @brief Sort position array by x on GPU.
 */
static inline void sort_xpos_dev(const int num, samplePos src, samplePos dst)
{
  __NOTE__("%s\n", "start");

  int Nrem = BLOCKSIZE(num, NTHREADS_SETIDX);
  int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
  int hidx = 0;
  for(int iter = 0; iter < Niter; iter++){
    int Nblck = (Nrem < MAX_BLOCKS_PER_GRID) ? Nrem : MAX_BLOCKS_PER_GRID;
    int Nsub = Nblck * NTHREADS_SETIDX;
    setIndex_kernel<<<Nblck, NTHREADS_SETIDX>>>(num, src.i_dev, hidx);

    hidx += Nsub;
    Nrem -= Nblck;
  }/* for(int iter = 0; iter < Niter; iter++){ */
  getLastCudaError("setIndex_kernel");

  /** sort using thrust */
  thrust::stable_sort_by_key((thrust::device_ptr<float>)(src.x_dev), (thrust::device_ptr<float>)((src.x_dev) + num), (thrust::device_ptr<int>)(src.i_dev));
  checkCudaErrors(cudaMemcpy(dst.x_dev, src.x_dev, num * sizeof(float), cudaMemcpyDeviceToDevice));

  Nrem = BLOCKSIZE(num, NTHREADS_SORTYZ);
  Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
  hidx = 0;
  for(int iter = 0; iter < Niter; iter++){
    int Nblck = (Nrem < MAX_BLOCKS_PER_GRID) ? Nrem : MAX_BLOCKS_PER_GRID;
    int Nsub = Nblck * NTHREADS_SORTYZ;
    sortSamplePos_yz_kernel<<<Nblck, NTHREADS_SORTYZ>>>(num, src.i_dev, dst.y_dev, src.y_dev, dst.z_dev, src.z_dev, hidx);

    hidx += Nsub;
    Nrem -= Nblck;
  }/* for(int iter = 0; iter < Niter; iter++){ */
  getLastCudaError("sortSamplePos_yz_kernel");

  __NOTE__("%s\n", "end");
}

/**
 * @fn sort_xpos_hst
 *
 * @brief Sort position array (x) on CPU.
 */
static inline void sort_xpos_hst(const int num, samplePos RESTRICT src, samplePos RESTRICT dst)
{
  __NOTE__("%s\n", "start");

  for(int ii = 0; ii < num; ii++)
    src.i_hst[ii] = ii;

  /** sort using thrust */
  thrust::stable_sort_by_key(src.x_hst, (src.x_hst) + num, src.i_hst);

  for(int ii = 0; ii < num; ii++){
    const int jj  = src.i_hst[ii];
    dst.x_hst[ii] = src.x_hst[ii];
    dst.y_hst[ii] = src.y_hst[jj];
    dst.z_hst[ii] = src.z_hst[jj];
  }/* for(int ii = 0; ii < num; ii++){ */

#ifndef NDEBUG
#if 0
  __PRINTF__("src:\n");
  for(int ii = 0; ii < num; ii++)
    __PRINTF__("%d: x = %e, y = %e, z = %e\n", ii, src.x_hst[ii], src.y_hst[ii], src.z_hst[ii]);

  __PRINTF__("\ndst:\n");
  for(int ii = 0; ii < num; ii++)
    __PRINTF__("%d: x = %e, y = %e, z = %e\n", ii, dst.x_hst[ii], dst.y_hst[ii], dst.z_hst[ii]);
  MPI_Abort(MPI_COMM_WORLD, 0);
  exit(0);
#endif
#endif//NDEBUG

  __NOTE__("%s\n", "end");
}

/**
 * @fn sort_xpos
 *
 * @brief Sort position array by x.
 */
static inline void sort_xpos(const int num, samplePos * RESTRICT src, samplePos * RESTRICT dst, const bool device)
{
  __NOTE__("%s\n", "start");

  if( device )    sort_xpos_dev(num, *src, *dst);
  else            sort_xpos_hst(num, *src, *dst);

  /** swap the list structure */
  samplePos _tmp;
  _tmp = *src;
  *src = *dst;
  *dst = _tmp;

  __NOTE__("%s\n", "end");
}


/**
 * @fn sort_ypos_dev
 *
 * @brief Sort position array by y on GPU.
 */
static inline void sort_ypos_dev(const int num, samplePos src)
{
  __NOTE__("%s\n", "start");

  /** sort using thrust */
  thrust::stable_sort_by_key((thrust::device_ptr<float>)(src.y_dev), (thrust::device_ptr<float>)((src.y_dev) + num), (thrust::device_ptr<float>)(src.z_dev));

#ifndef NDEBUG
#if 0
  checkCudaErrors(cudaMemcpy(src.y_hst, src.y_dev, sizeof(float) * num, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(src.z_hst, src.z_dev, sizeof(float) * num, cudaMemcpyDeviceToHost));
  for(int ii = 0; ii < num; ii++)
    __PRINTF__("%d: y = %e, z = %e\n", ii, src.y_hst[ii], src.z_hst[ii]);
  MPI_Abort(MPI_COMM_WORLD, 0);
  exit(0);
#endif
#endif//NDEBUG

  __NOTE__("%s\n", "end");
}

/**
 * @fn sort_ypos_hst
 *
 * @brief Sort position array (y) on CPU.
 */
static inline void sort_ypos_hst(const int num, samplePos src)
{
  __NOTE__("%s\n", "start");

  /** sort using thrust */
  thrust::stable_sort_by_key(src.y_hst, (src.y_hst) + num, src.z_hst);

#ifndef NDEBUG
#if 0
  for(int ii = 0; ii < num; ii++)
    __PRINTF__("%d: y = %e, z = %e\n", ii, src.y_hst[ii], src.z_hst[ii]);
  MPI_Abort(MPI_COMM_WORLD, 0);
  exit(0);
#endif
#endif//NDEBUG

  __NOTE__("%s\n", "end");
}

/**
 * @fn sort_ypos
 *
 * @brief Sort position array by y.
 */
static inline void sort_ypos(const int num, samplePos * RESTRICT src, const bool device)
{
  __NOTE__("%s\n", "start");

  if( device )    sort_ypos_dev(num, *src);
  else            sort_ypos_hst(num, *src);

  __NOTE__("%s\n", "end");
}


/**
 * @fn sort_zpos_dev
 *
 * @brief Sort position array by z on GPU.
 */
static inline void sort_zpos_dev(const int num, samplePos src)
{
  __NOTE__("%s\n", "start");

  /** sort using thrust */
  thrust::stable_sort((thrust::device_ptr<float>)(src.z_dev), (thrust::device_ptr<float>)((src.z_dev) + num));

  __NOTE__("%s\n", "end");
}

/**
 * @fn sort_zpos_hst
 *
 * @brief Sort position array (z) on CPU.
 */
static inline void sort_zpos_hst(const int num, samplePos src)
{
  __NOTE__("%s\n", "start");

  /** sort using thrust */
  thrust::stable_sort(src.z_hst, (src.z_hst) + num);

#ifndef NDEBUG
#if 0
  for(int ii = 0; ii < num; ii++)
    __PRINTF__("%d: z = %e\n", ii, src.z_hst[ii]);
  MPI_Abort(MPI_COMM_WORLD, 0);
  exit(0);
#endif
#endif//NDEBUG

  __NOTE__("%s\n", "end");
}

/**
 * @fn sort_zpos
 *
 * @brief Sort position array by z.
 */
static inline void sort_zpos(const int num, samplePos * RESTRICT src, const bool device)
{
  __NOTE__("%s\n", "start");

  if( device )    sort_zpos_dev(num, *src);
  else            sort_zpos_hst(num, *src);

  __NOTE__("%s\n", "end");
}


/**
 * @fn allocateSamplePos
 *
 * @brief Memory allocation for sampling particles.
 */
extern "C"
muse allocateSamplePos
(float **x0hst, float **x1hst, float **y0hst, float **y1hst, float **z0hst, float **z1hst, int **idhst,
 float **x0dev, float **x1dev, float **y0dev, float **y1dev, float **z0dev, float **z1dev, int **iddev,
 samplePos * RESTRICT pos0, samplePos * RESTRICT pos1, const sampling sample)
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  mycudaMallocHost((void **)x0hst, sample.Nmax * sizeof(float));  alloc.host   += sample.Nmax * sizeof(float);
  mycudaMallocHost((void **)x1hst, sample.Nmax * sizeof(float));  alloc.host   += sample.Nmax * sizeof(float);
  mycudaMallocHost((void **)y0hst, sample.Nmax * sizeof(float));  alloc.host   += sample.Nmax * sizeof(float);
  mycudaMallocHost((void **)y1hst, sample.Nmax * sizeof(float));  alloc.host   += sample.Nmax * sizeof(float);
  mycudaMallocHost((void **)z0hst, sample.Nmax * sizeof(float));  alloc.host   += sample.Nmax * sizeof(float);
  mycudaMallocHost((void **)z1hst, sample.Nmax * sizeof(float));  alloc.host   += sample.Nmax * sizeof(float);
  mycudaMallocHost((void **)idhst, sample.Nmax * sizeof(  int));  alloc.host   += sample.Nmax * sizeof(  int);
  mycudaMalloc    ((void **)x0dev, sample.Nmax * sizeof(float));  alloc.device += sample.Nmax * sizeof(float);
  mycudaMalloc    ((void **)x1dev, sample.Nmax * sizeof(float));  alloc.device += sample.Nmax * sizeof(float);
  mycudaMalloc    ((void **)y0dev, sample.Nmax * sizeof(float));  alloc.device += sample.Nmax * sizeof(float);
  mycudaMalloc    ((void **)y1dev, sample.Nmax * sizeof(float));  alloc.device += sample.Nmax * sizeof(float);
  mycudaMalloc    ((void **)z0dev, sample.Nmax * sizeof(float));  alloc.device += sample.Nmax * sizeof(float);
  mycudaMalloc    ((void **)z1dev, sample.Nmax * sizeof(float));  alloc.device += sample.Nmax * sizeof(float);
  mycudaMalloc    ((void **)iddev, sample.Nmax * sizeof(  int));  alloc.device += sample.Nmax * sizeof(  int);

  pos0->x_hst = *x0hst;  pos0->y_hst = *y0hst;  pos0->z_hst = *z0hst;  pos0->i_hst = *idhst;
  pos0->x_dev = *x0dev;  pos0->y_dev = *y0dev;  pos0->z_dev = *z0dev;  pos0->i_dev = *iddev;
  pos1->x_hst = *x1hst;  pos1->y_hst = *y1hst;  pos1->z_hst = *z1hst;  pos1->i_hst = *idhst;
  pos1->x_dev = *x1dev;  pos1->y_dev = *y1dev;  pos1->z_dev = *z1dev;  pos1->i_dev = *iddev;


  __NOTE__("%s\n", "end");
  return (alloc);
}

/**
 * @fn releaseSamplePos
 *
 * @brief Memory deallocation for sampling particles.
 */
extern "C"
void  releaseSamplePos
(float  *x0hst, float  *x1hst, float  *y0hst, float  *y1hst, float  *z0hst, float  *z1hst, int  *idhst,
 float  *x0dev, float  *x1dev, float  *y0dev, float  *y1dev, float  *z0dev, float  *z1dev, int  *iddev)
{
  __NOTE__("%s\n", "start");

  mycudaFreeHost(x0hst);  mycudaFreeHost(y0hst);  mycudaFreeHost(z0hst);
  mycudaFreeHost(x1hst);  mycudaFreeHost(y1hst);  mycudaFreeHost(z1hst);  mycudaFreeHost(idhst);
  mycudaFree    (x0dev);  mycudaFree    (y0dev);  mycudaFree    (z0dev);
  mycudaFree    (x1dev);  mycudaFree    (y1dev);  mycudaFree    (z1dev);  mycudaFree    (iddev);

  __NOTE__("%s\n", "end");
}


#if 0
/**
 * @fn setParticlePosition_kernel
 *
 * @brief Set particle position.
 */
__global__ void setParticlePosition_kernel(const int num, READ_ONLY position * ipos, float * RESTRICT xi, float * RESTRICT yi, float * RESTRICT zi, const int offset)
{
  const int ii = offset + GLOBALIDX_X1D;

  if( ii < num ){
    const position pi = ipos[ii];
    xi[ii] = CAST_R2F(pi.x);
    yi[ii] = CAST_R2F(pi.y);
    zi[ii] = CAST_R2F(pi.z);
  }/* if( ii < num ){ */
}

/**
 * @fn copyParticlePositionAsync_dev2hst
 *
 * @brief Copy particle position.
 */
static inline void copyParticlePositionAsync_dev2hst(const int Ni, position * RESTRICT ipos, particlePos dev, particlePos hst, cudaStream_t stream)
{
  __NOTE__("%s\n", "start");

  int Nrem = BLOCKSIZE(Ni, NTHREADS_SETPOS);
  const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
  int hidx = 0;
  for(int iter = 0; iter < Niter; iter++){
    int Nblck = (Nrem < MAX_BLOCKS_PER_GRID) ? Nrem : MAX_BLOCKS_PER_GRID;
    int Nsub = Nblck * NTHREADS_SETPOS;
    setParticlePosition_kernel<<<Nblck, NTHREADS_SETPOS, SMEM_SIZE, stream>>>(Ni, ipos, dev.x, dev.y, dev.z, hidx);

    hidx += Nsub;
    Nrem -= Nblck;
  }/* for(int iter = 0; iter < Niter; iter++){ */
  getLastCudaError("setParticlePosition_kernel");

  checkCudaErrors(cudaMemcpyAsync(hst.x, dev.x, sizeof(float) * Ni, cudaMemcpyDeviceToHost, stream));
  checkCudaErrors(cudaMemcpyAsync(hst.y, dev.y, sizeof(float) * Ni, cudaMemcpyDeviceToHost, stream));
  checkCudaErrors(cudaMemcpyAsync(hst.z, dev.z, sizeof(float) * Ni, cudaMemcpyDeviceToHost, stream));

  __NOTE__("%s\n", "end");
}
#endif


#define UNKNOWN_DST_RANK (-1)
__global__ void initDstRank(const int num, int * dstRank)
{
  const int ii = GLOBALIDX_X1D;
  if( ii < num )
    dstRank[ii] = UNKNOWN_DST_RANK;
}

/**
 * @struct union_i4_f4
 *
 * @brief union for reducing usage of shared memory
 */
typedef union __align__(16)
{
  int4 i4;
  float4 f4;
} union_i4_f4;
#define NTHREADS_SCAN_VEC4_INC NTHREADS_ASSIGN
#include "../util/scan_vec4_inc.cu"
/**
 * @fn assignNewDomain_kernel
 *
 * @brief Assign N-body particles to new computational domain on GPU.
 */
__global__ void
#ifndef USE_OCCUPANCY_CALCULATOR
__launch_bounds__(NTHREADS_ASSIGN, NBLOCKS_PER_SM_ASSIGN)
#endif//USE_OCCUPANCY_CALCULATOR
assignNewDomain_kernel
     (const int num, int * RESTRICT numNew_gm, READ_ONLY position * RESTRICT ipos,
      const int domHead, const int domRem, READ_ONLY int * RESTRICT domRank, int * RESTRICT dstRank, READ_ONLY float * RESTRICT xmin, READ_ONLY float * RESTRICT xmax, READ_ONLY float * RESTRICT ymin, READ_ONLY float * RESTRICT ymax, READ_ONLY float * RESTRICT zmin, READ_ONLY float * RESTRICT zmax,
      int4 * RESTRICT gmem, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int bnum =   GRIDDIM_X1D;

  const int hidx = tidx - (tidx & (warpSize - 1));

  const int ihead = (num *      bidx ) / bnum;
  const int itail = (num * (1 + bidx)) / bnum;

  /* __shared__ int4 num_sm[NTHREADS_ASSIGN]; */
  __shared__ union_i4_f4 num_sm[NTHREADS_ASSIGN];

  /** load boundaries of domains */
  const int domNum = (domRem >= 4) ? 4 : domRem;
  /* __shared__ float3 boxmin_sm[4], boxmax_sm[4]; */
  /* __shared__ int4 rank_sm;/\* corresponds to sendBuf[jj].rank *\/ */
  if( tidx < domNum ){
    /* boxmin_sm[tidx].x = xmin[domHead + tidx]; */
    /* boxmax_sm[tidx].x = xmax[domHead + tidx]; */
    /* boxmin_sm[tidx].y = ymin[domHead + tidx]; */
    /* boxmax_sm[tidx].y = ymax[domHead + tidx]; */
    /* boxmin_sm[tidx].z = zmin[domHead + tidx]; */
    /* boxmax_sm[tidx].z = zmax[domHead + tidx]; */
    num_sm[    tidx].f4.x = xmin[domHead + tidx];
    num_sm[4 + tidx].f4.x = xmax[domHead + tidx];
    num_sm[    tidx].f4.y = ymin[domHead + tidx];
    num_sm[4 + tidx].f4.y = ymax[domHead + tidx];
    num_sm[    tidx].f4.z = zmin[domHead + tidx];
    num_sm[4 + tidx].f4.z = zmax[domHead + tidx];

    switch( tidx ){
    /* case 0:      rank_sm.x = domRank[domHead    ];      break; */
    /* case 1:      rank_sm.y = domRank[domHead + 1];      break; */
    /* case 2:      rank_sm.z = domRank[domHead + 2];      break; */
    /* case 3:      rank_sm.w = domRank[domHead + 3];      break; */
    case 0:      num_sm[8].i4.x = domRank[domHead    ];      break;
    case 1:      num_sm[8].i4.y = domRank[domHead + 1];      break;
    case 2:      num_sm[8].i4.z = domRank[domHead + 2];      break;
    case 3:      num_sm[8].i4.w = domRank[domHead + 3];      break;
    }/* switch( tidx ){ */
  }/* if( tidx < domNum ){ */
  else{
    if( tidx < 4 ){
      /* boxmin_sm[tidx].x =  0.25f * FLT_MAX; */
      /* boxmax_sm[tidx].x = -0.25f * FLT_MAX; */
      /* boxmin_sm[tidx].y =  0.25f * FLT_MAX; */
      /* boxmax_sm[tidx].y = -0.25f * FLT_MAX; */
      /* boxmin_sm[tidx].z =  0.25f * FLT_MAX; */
      /* boxmax_sm[tidx].z = -0.25f * FLT_MAX; */
      num_sm[    tidx].f4.x =  0.25f * FLT_MAX;
      num_sm[4 + tidx].f4.x = -0.25f * FLT_MAX;
      num_sm[    tidx].f4.y =  0.25f * FLT_MAX;
      num_sm[4 + tidx].f4.y = -0.25f * FLT_MAX;
      num_sm[    tidx].f4.z =  0.25f * FLT_MAX;
      num_sm[4 + tidx].f4.z = -0.25f * FLT_MAX;

      switch( tidx ){
      /* case 3:	rank_sm.x = -1;	break; */
      /* case 2:	rank_sm.x = -1;	break; */
      /* case 1:	rank_sm.x = -1;	break; */
      /* case 0:	rank_sm.x = -1;	break; */
      case 3:	num_sm[8].i4.x = -1;	break;
      case 2:	num_sm[8].i4.x = -1;	break;
      case 1:	num_sm[8].i4.x = -1;	break;
      case 0:	num_sm[8].i4.x = -1;	break;
      }/* switch( tidx ){ */
    }/* if( tidx < 4 ){ */
  }/* else{ */
  __syncthreads();


  /** determine process rank for each particle to belong */
  /* const int4 rank = rank_sm; */
  const int4 rank = num_sm[8].i4;
  int4 numNew = {0, 0, 0, 0};
  for(int ih = ihead; ih < itail; ih += NTHREADS_ASSIGN){
    const int ii = ih + tidx;
    if( (ii < itail) && (dstRank[ii] == UNKNOWN_DST_RANK) ){
      const position pi = ipos[ii];
      const float xi = CAST_R2F(pi.x);
      const float yi = CAST_R2F(pi.y);
      const float zi = CAST_R2F(pi.z);

      for(int jj = 0; jj < domNum; jj++){
	/* if( (xi >= boxmin_sm[jj].x) && (xi <= boxmax_sm[jj].x) && */
	/*     (yi >= boxmin_sm[jj].y) && (yi <= boxmax_sm[jj].y) && */
	/*     (zi >= boxmin_sm[jj].z) && (zi <= boxmax_sm[jj].z) ){ */
	if( (xi >= num_sm[jj].f4.x) && (xi <= num_sm[4 + jj].f4.x) &&
	    (yi >= num_sm[jj].f4.y) && (yi <= num_sm[4 + jj].f4.y) &&
	    (zi >= num_sm[jj].f4.z) && (zi <= num_sm[4 + jj].f4.z) ){

	  switch( jj ){
	  case 0:	    dstRank[ii] = rank.x;	    numNew.x++;	    break;
	  case 1:	    dstRank[ii] = rank.y;	    numNew.y++;	    break;
	  case 2:	    dstRank[ii] = rank.z;	    numNew.z++;	    break;
	  case 3:	    dstRank[ii] = rank.w;	    numNew.w++;	    break;
	  }/* switch(jj){ */

	  break;
	}
      }/* for(int jj = 0; jj < domNum; jj++){ */
    }/* if( (ii < itail) && (dstRank[ii] == UNKNOWN_DST_RANK) ){ */
  }/* for(int ih = ihead; ih < itail; ih += NTHREADS_ASSIGN){ */
  __syncthreads();/* additional synchronization to use union */


  /** reduction about numNew */
  numNew = TOTAL_SUM_VEC4_GRID(numNew, (int4 *)num_sm, tidx, hidx, gmem, bidx, bnum, gsync0, gsync1);


  /** upload number of N-body particles to be assigned to new computational domain */
  if( (bidx + tidx) == 0 ){
    numNew_gm[domHead    ] = numNew.x;
    numNew_gm[domHead + 1] = numNew.y;
    numNew_gm[domHead + 2] = numNew.z;
    numNew_gm[domHead + 3] = numNew.w;
  }/* if( (bidx + tidx) == 0 ){ */
}


/**
 * @fn sortParticlesDDkey_kernel
 *
 * @brief Sort particle data by domain index on GPU.
 */
__global__ void sortParticlesDDkey_kernel
(const int num, READ_ONLY int * RESTRICT old,
 ulong        * RESTRICT didx , READ_ONLY ulong        * RESTRICT sidx ,
 position     * RESTRICT dpos , READ_ONLY position     * RESTRICT spos ,
#ifdef  GADGET_MAC
 acceleration * RESTRICT dacc , READ_ONLY acceleration * RESTRICT sacc ,
#endif//GADGET_MAC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 acceleration * RESTRICT dext , READ_ONLY acceleration * RESTRICT sext ,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
 velocity     * RESTRICT dvel , READ_ONLY velocity     * RESTRICT svel ,
 ibody_time   * RESTRICT dtime, READ_ONLY ibody_time   * RESTRICT stime
#else///BLOCK_TIME_STEP
 real         * RESTRICT dvx  , READ_ONLY real         * RESTRICT svx  ,
 real         * RESTRICT dvy  , READ_ONLY real         * RESTRICT svy  ,
 real         * RESTRICT dvz  , READ_ONLY real         * RESTRICT svz
#endif//BLOCK_TIME_STEP
 , const int offset)
{
  const int ii = offset + GLOBALIDX_X1D;

  if( ii < num ){
    /** load old tag */
    const int jj = old[ii];

    didx [ii] = sidx [jj];
    dpos [ii] = spos [jj];
#ifdef  GADGET_MAC
    dacc [ii] = sacc [jj];
#endif//GADGET_MAC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
    dext [ii] = sext [jj];
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
    dvel [ii] = svel [jj];
    dtime[ii] = stime[jj];
#else///BLOCK_TIME_STEP
    dvx  [ii] = svx  [jj];
    dvy  [ii] = svy  [jj];
    dvz  [ii] = svz  [jj];
#endif//BLOCK_TIME_STEP

#ifndef NDEBUG
    if( spos[jj].m < EPSILON )
      printf("src = %d, dst = %d: massless (%e)\n", jj, ii, spos[jj].m);
#endif//NDEBUG
  }/* if( ii < num ){ */
}

/**
 * @fn sortDomainDecomposeKey
 *
 * @brief Sort particle data by domain index on GPU.
 *
 * @sa setIndex_kernel
 * @sa sortParticlesDDkey_kernel
 */
static inline void sortDomainDecomposeKey(const int num, domainDecomposeKey key, iparticle * RESTRICT src, iparticle * RESTRICT dst)
{
  __NOTE__("%s\n", "start");


  checkCudaErrors(cudaMemcpy(key.dstRank_dev, key.dstRank_hst, num * sizeof(int), cudaMemcpyHostToDevice));
  int Nrem = BLOCKSIZE(num, NTHREADS_SETIDX);
  int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
  int hidx = 0;
  for(int iter = 0; iter < Niter; iter++){
    int Nblck = (Nrem < MAX_BLOCKS_PER_GRID) ? Nrem : MAX_BLOCKS_PER_GRID;
    int Nsub = Nblck * NTHREADS_SETIDX;
    setIndex_kernel<<<Nblck, NTHREADS_SETIDX>>>(num, key.bodyIdx_dev, hidx);

    hidx += Nsub;
    Nrem -= Nblck;
  }/* for(int iter = 0; iter < Niter; iter++){ */
  getLastCudaError("setIndex_kernel");

  thrust::stable_sort_by_key((thrust::device_ptr<int>)(key.dstRank_dev), (thrust::device_ptr<int>)(key.dstRank_dev + num), (thrust::device_ptr<int>)(key.bodyIdx_dev));

  Nrem = BLOCKSIZE(num, NTHREADS_DDSORT);
  Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
  hidx = 0;
  for(int iter = 0; iter < Niter; iter++){
    int Nblck = (Nrem < MAX_BLOCKS_PER_GRID) ? Nrem : MAX_BLOCKS_PER_GRID;
    int Nsub = Nblck * NTHREADS_DDSORT;
    sortParticlesDDkey_kernel<<<Nblck, NTHREADS_DDSORT>>>
      (num, key.bodyIdx_dev,
       (*dst).idx , (*src).idx,
       (*dst).pos , (*src).pos,
#ifdef  GADGET_MAC
       (*dst).acc , (*src).acc,
#endif//GADGET_MAC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
       (*dst).acc_ext, (*src).acc_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
       (*dst).vel , (*src).vel,
       (*dst).time, (*src).time
#else///BLOCK_TIME_STEP
       (*dst).vx  , (*src).vx,
       (*dst).vy  , (*src).vy,
       (*dst).vz  , (*src).vz
#endif//BLOCK_TIME_STEP
       , hidx);

    hidx += Nsub;
    Nrem -= Nblck;
  }/* for(int iter = 0; iter < Niter; iter++){ */
  getLastCudaError("sortParticlesDDkey_kernel");


  /** swap the list structure */
  iparticle _tmp;
  _tmp = *src;
  *src = *dst;
  *dst = _tmp;


  __NOTE__("%s\n", "end");
}


/**
 * @fn allocateParticlePosition
 *
 * @brief Memory allocation particle position.
 */
extern "C"
muse allocateParticlePosition(float **xhst, float **yhst, float **zhst, particlePos *hst,
			      float **xdev, float **ydev, float **zdev, particlePos *dev,
			      int **rank_hst, int **rank_dev, int **idx_dev, domainDecomposeKey *key, const size_t num)
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  /* size_t num = (size_t)((float)Ntot * MAX_FACTOR_FROM_EQUIPARTITION); */
  size_t size = num;
  if( (num % NTHREADS) != 0 )
    size += NTHREADS - (num % NTHREADS);

  mycudaMallocHost((void **)xhst, size * sizeof(float));  alloc.host   += size * sizeof(float);
  mycudaMallocHost((void **)yhst, size * sizeof(float));  alloc.host   += size * sizeof(float);
  mycudaMallocHost((void **)zhst, size * sizeof(float));  alloc.host   += size * sizeof(float);
  mycudaMalloc    ((void **)xdev, size * sizeof(float));  alloc.device += size * sizeof(float);
  mycudaMalloc    ((void **)ydev, size * sizeof(float));  alloc.device += size * sizeof(float);
  mycudaMalloc    ((void **)zdev, size * sizeof(float));  alloc.device += size * sizeof(float);

  hst->x = *xhst;  hst->y = *yhst;  hst->z = *zhst;
  dev->x = *xdev;  dev->y = *ydev;  dev->z = *zdev;

  mycudaMallocHost((void **)rank_hst, size * sizeof(int));  alloc.host   += size * sizeof(int);
  mycudaMalloc    ((void **)rank_dev, size * sizeof(int));  alloc.device += size * sizeof(int);
  mycudaMalloc    ((void **) idx_dev, size * sizeof(int));  alloc.device += size * sizeof(int);

  key->dstRank_hst = *rank_hst;
  key->dstRank_dev = *rank_dev;
  key->bodyIdx_dev = * idx_dev;


  __NOTE__("%s\n", "end");
  return (alloc);
}

/**
 * @fn releaseParticlePosition
 *
 * @brief Memory deallocation particle position.
 */
extern "C"
void  releaseParticlePosition(float  *xhst, float  *yhst, float  *zhst,
			      float  *xdev, float  *ydev, float	 *zdev,
			      int  *rank_hst, int  *rank_dev, int  *idx_dev)
{
  __NOTE__("%s\n", "start");

  mycudaFreeHost(xhst);  mycudaFreeHost(yhst);  mycudaFreeHost(zhst);
  mycudaFree    (xdev);  mycudaFree    (ydev);  mycudaFree    (zdev);

  mycudaFreeHost(rank_hst);
  mycudaFree    (rank_dev);  mycudaFree(idx_dev);

  __NOTE__("%s\n", "end");
}


#ifndef USE_OCCUPANCY_CALCULATOR
/**
 * @fn checkDomainPos_dev
 *
 * @brief Error detection before running the kernel
 *
 * @sa assignNewDomain_kernel
 */
void checkDomainPos_dev(const deviceProp devProp)
{
  __NOTE__("%s\n", "start");

  struct cudaFuncAttributes funcAttr;
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, assignNewDomain_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_ASSIGN ){
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in assignNewDomain_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_ASSIGN);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_ASSIGN = %d.\n", GPUGEN, GPUVER, NTHREADS_ASSIGN);
    fflush (stderr);
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_ASSIGN ){ */

  int regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS_ASSIGN);
  if( regLimit > (MAX_REGISTERS_PER_SM / NTHREADS_ASSIGN) )
    regLimit = (MAX_REGISTERS_PER_SM / NTHREADS_ASSIGN);
  int memLimit = SMEM_SIZE_L1_PREF / funcAttr.sharedSizeBytes;
  int Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / NTHREADS_ASSIGN) )    Nblck = MAX_THREADS_PER_SM       / NTHREADS_ASSIGN;
  if( Nblck >   MAX_BLOCKS_PER_SM                      )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / NTHREADS_ASSIGN) )    Nblck = ((MAX_WARPS_PER_SM * 32) / NTHREADS_ASSIGN);

  if( Nblck != NBLOCKS_PER_SM_ASSIGN ){
    __KILL__(stderr, "ERROR: # of blocks per SM for assignNewDomain_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_ASSIGN defined in src/para/exchange_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_ASSIGN, funcAttr.numRegs, REGISTERS_PER_THREAD_ASSIGN, GPUGEN, GPUVER);
  }/* if( Nblck != NBLOCKS_PER_SM ){ */

  if( (devProp.numSM * NBLOCKS_PER_SM_ASSIGN) > NTHREADS_ASSIGN ){
    __KILL__(stderr, "ERROR: product (%d) of devProp.numSM(%d) * NBLOCKS_PER_SM_ASSIGN(%d) must be smaller than NTHREADS_ASSIGN(%d) to use shared memory.\n", devProp.numSM * NBLOCKS_PER_SM_ASSIGN, devProp.numSM, NBLOCKS_PER_SM_ASSIGN, NTHREADS_ASSIGN);
  }/* if( (devProp.numSM * NBLOCKS_PER_SM_ASSIGN) > NTHREADS_ASSIGN ){ */


  __NOTE__("%s\n", "end");
}
#endif//USE_OCCUPANCY_CALCULATOR


/**
 * @fn allocateDomainPos
 *
 * @brief Memory allocation to send domain boundary toward device.
 */
extern "C"
muse allocateDomainPos(float **xmin_dev, float **xmax_dev, float **ymin_dev, float **ymax_dev, float **zmin_dev, float **zmax_dev,
		       float **xmin_hst, float **xmax_hst, float **ymin_hst, float **ymax_hst, float **zmin_hst, float **zmax_hst,
		       int **rank, int **rank_hst, int **numNew, int **numNew_hst, int4 **gmem, int **gsync0, int **gsync1,
		       sendDom *dom, const int Ngpu, const deviceProp devProp)
{
  __NOTE__("%s\n", "start");


#   if  GPUGEN >= 70
  /* remove shared memory if __global__ function does not use */
  checkCudaErrors(cudaFuncSetCacheConfig(     pickupSamples_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, 0));
  checkCudaErrors(cudaFuncSetCacheConfig(          setIndex_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, 0));
  checkCudaErrors(cudaFuncSetCacheConfig(  sortSamplePos_yz_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, 0));
  checkCudaErrors(cudaFuncSetCacheConfig(              initDstRank, cudaFuncAttributeMaxDynamicSharedMemorySize, 0));
  checkCudaErrors(cudaFuncSetCacheConfig(sortParticlesDDkey_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, 0));
#endif//GPUGEN >= 70


  muse alloc = {0, 0};

  mycudaMalloc    ((void **)xmin_dev, Ngpu * sizeof(float));  alloc.device += Ngpu * sizeof(float);  mycudaMalloc    ((void **)xmax_dev, Ngpu * sizeof(float));  alloc.device += Ngpu * sizeof(float);
  mycudaMalloc	  ((void **)ymin_dev, Ngpu * sizeof(float));  alloc.device += Ngpu * sizeof(float);  mycudaMalloc    ((void **)ymax_dev, Ngpu * sizeof(float));  alloc.device += Ngpu * sizeof(float);
  mycudaMalloc	  ((void **)zmin_dev, Ngpu * sizeof(float));  alloc.device += Ngpu * sizeof(float);  mycudaMalloc    ((void **)zmax_dev, Ngpu * sizeof(float));  alloc.device += Ngpu * sizeof(float);
  mycudaMallocHost((void **)xmin_hst, Ngpu * sizeof(float));  alloc.host   += Ngpu * sizeof(float);  mycudaMallocHost((void **)xmax_hst, Ngpu * sizeof(float));  alloc.host   += Ngpu * sizeof(float);
  mycudaMallocHost((void **)ymin_hst, Ngpu * sizeof(float));  alloc.host   += Ngpu * sizeof(float);  mycudaMallocHost((void **)ymax_hst, Ngpu * sizeof(float));  alloc.host   += Ngpu * sizeof(float);
  mycudaMallocHost((void **)zmin_hst, Ngpu * sizeof(float));  alloc.host   += Ngpu * sizeof(float);  mycudaMallocHost((void **)zmax_hst, Ngpu * sizeof(float));  alloc.host   += Ngpu * sizeof(float);

  dom->xmin_dev = *xmin_dev;  dom->xmax_dev = *xmax_dev;  dom->xmin_hst = *xmin_hst;  dom->xmax_hst = *xmax_hst;
  dom->ymin_dev = *ymin_dev;  dom->ymax_dev = *ymax_dev;  dom->ymin_hst = *ymin_hst;  dom->ymax_hst = *ymax_hst;
  dom->zmin_dev = *zmin_dev;  dom->zmax_dev = *zmax_dev;  dom->zmin_hst = *zmin_hst;  dom->zmax_hst = *zmax_hst;


  mycudaMalloc    ((void **)rank    , Ngpu * sizeof(int));  alloc.device += Ngpu * sizeof(int);
  mycudaMallocHost((void **)rank_hst, Ngpu * sizeof(int));  alloc.host   += Ngpu * sizeof(int);

  /* number of elements for *numNew must be multiple of 4 */
  const size_t num = ((Ngpu % 4) == 0) ? Ngpu : (Ngpu + (4 - (Ngpu % 4)));
  mycudaMalloc    ((void **)numNew    , num * sizeof(int));  alloc.device += num * sizeof(int);
  mycudaMallocHost((void **)numNew_hst, num * sizeof(int));  alloc.host   += num * sizeof(int);

#ifdef  USE_OCCUPANCY_CALCULATOR
  checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&dom->numBlocksPerSM_assign, assignNewDomain_kernel, NTHREADS_ASSIGN, 0));
  const int NBLOCKS_PER_SM_ASSIGN = dom->numBlocksPerSM_assign;
#endif//USE_OCCUPANCY_CALCULATOR
  const size_t size = devProp.numSM * NBLOCKS_PER_SM_ASSIGN;
  mycudaMalloc((void **)gmem  , size * sizeof(int4));  alloc.device += size * sizeof(int4);
  mycudaMalloc((void **)gsync0, size * sizeof(int ));  alloc.device += size * sizeof(int);
  mycudaMalloc((void **)gsync1, size * sizeof(int ));  alloc.device += size * sizeof(int);

  dom->rank   = *rank  ;  dom->  rank_hst = *  rank_hst;
  dom->numNew = *numNew;  dom->numNew_hst = *numNew_hst;
  dom->gmem   = *gmem;
  dom->gsync0 = *gsync0;
  dom->gsync1 = *gsync1;

  /* check # of blocks launched in assignNewDomain_kernel(); */
#ifndef USE_OCCUPANCY_CALCULATOR
  checkDomainPos_dev(devProp);
#endif//USE_OCCUPANCY_CALCULATOR


  __NOTE__("%s\n", "end");
  return (alloc);
}

/**
 * @fn releaseDomainPos
 *
 * @brief Memory deallocation.
 */
extern "C"
void  releaseDomainPos(float  *xmin_dev, float  *xmax_dev, float  *ymin_dev, float  *ymax_dev, float  *zmin_dev, float  *zmax_dev,
		       float  *xmin_hst, float  *xmax_hst, float  *ymin_hst, float  *ymax_hst, float  *zmin_hst, float  *zmax_hst,
		       int  *rank, int  *rank_hst, int  *numNew, int  *numNew_hst, int4  *gmem, int  *gsync0, int  *gsync1)
{
  __NOTE__("%s\n", "start");

  mycudaFree    (xmin_dev);  mycudaFree    (xmax_dev);  mycudaFree    (ymin_dev);  mycudaFree    (ymax_dev);  mycudaFree    (zmin_dev);  mycudaFree    (zmax_dev);
  mycudaFreeHost(xmin_hst);  mycudaFreeHost(xmax_hst);  mycudaFreeHost(ymin_hst);  mycudaFreeHost(ymax_hst);  mycudaFreeHost(zmin_hst);  mycudaFreeHost(zmax_hst);

  mycudaFree(rank);  mycudaFreeHost(rank_hst);
  mycudaFree(numNew);  mycudaFreeHost(numNew_hst);
  mycudaFree(gmem);
  mycudaFree(gsync0);  mycudaFree(gsync1);

  __NOTE__("%s\n", "end");
}


#ifdef  ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
static inline float cutting(float ll, float rr, float min, float dL, float dLinv){  return (min + dL * nearbyintf(((0.5f * (ll + rr)) - min) * dLinv));}
#else///ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
static inline float cutting(float ll, float rr                                  ){  return (                        0.5f * (ll + rr)                 );}
#endif//ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX


/**
 * @fn exchangeParticles_dev
 *
 * @brief Exchange N-body particles among multiple GPUs with MPI.
 *
 * @sa copyParticlePositionAsync_dev2hst
 * @sa getBoxSize_dev
 * @sa sort_xpos
 * @sa sort_ypos
 * @sa sort_zpos
 * @sa cutting
 * @sa sortDomainDecomposeKey
 */
extern "C"
void exchangeParticles_dev
(const int numOld, const ulong Ntot, const int numMax, int *numNew,
 iparticle * RESTRICT src_dev, iparticle * RESTRICT dst_dev, iparticle * RESTRICT src_hst, iparticle * RESTRICT dst_hst,
 sendDom domBoundary, particlePos pos_hst, particlePos pos_dev, domainDecomposeKey key, sendCfg *sendBuf, recvCfg *recvBuf,
 MPIinfo orm[], MPIinfo rep[], domainCfg domain, MPIcfg_tree mpi,
 const double tloc, sampling sample, samplePos loc, samplePos ful,
 soaPHsort soa, const deviceProp devProp, const deviceInfo devInfo,
 measuredTime *measured
#ifdef  CARE_EXTERNAL_PARTICLES
 , domainLocation *location
#endif//CARE_EXTERNAL_PARTICLES
 , brentStatus *status, brentMemory *memory
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
 )
{
  __NOTE__("%s\n", "start");


  static struct timespec start;
  checkCudaErrors(cudaDeviceSynchronize());
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);

  const bool device = true;
  const bool   host = false;

  /** pick up sample particles */
  /** weight is determined using elapsed time by each process */
  __NOTE__("rank %d: tloc = %e, numOld = %d\n", mpi.rank, tloc, numOld);
  double ttot = tloc;
  chkMPIerr(MPI_Allreduce(&tloc, &ttot, 1, MPI_DOUBLE, MPI_SUM, mpi.comm));
  const float frac = fminf((float)(tloc / ttot), MAX_FACTOR_INCREASE / (float)mpi.size);
  const int Nsub = (int)ceilf((float)Ntot * sample.rate * frac);
  __NOTE__("rank %d: tloc = %e, ttot = %e, frac = %e, numOld = %d, Nsub = %d\n", mpi.rank, tloc, ttot, tloc / ttot, numOld, Nsub);

  const int iskip = (Nsub < numOld) ? (numOld / Nsub) : (1);
  int sendNum = numOld / iskip;
  int Nrem = BLOCKSIZE(sendNum, NTHREADS_PICKUP);
  int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
  int hidx = 0;
  for(int iter = 0; iter < Niter; iter++){
    int Nblck = (Nrem < MAX_BLOCKS_PER_GRID) ? Nrem : MAX_BLOCKS_PER_GRID;
    int Nsub = Nblck * NTHREADS_PICKUP;
    pickupSamples_kernel<<<Nblck, NTHREADS_PICKUP>>>(sendNum, iskip, (*src_dev).pos, loc.x_dev, loc.y_dev, loc.z_dev, hidx);

    hidx += Nsub;
    Nrem -= Nblck;
  }/* for(int iter = 0; iter < Niter; iter++){ */
  getLastCudaError("pickupSamples_kernel");
  __NOTE__("rank %d: iskip = %d, sendNum = %d\n", mpi.rank, iskip, sendNum);

  /** sort sample particles in each direction */
  if(             mpi.dim[0] != 1 )     sort_xpos(sendNum, &loc, &ful, device);
  else{	   if(	  mpi.dim[1] != 1 )	sort_ypos(sendNum, &loc      , device);
    else      if( mpi.dim[2] != 1 )	sort_zpos(sendNum, &loc      , device);
  }
#ifdef  MPI_VIA_HOST
  checkCudaErrors(cudaMemcpy(loc.x_hst, loc.x_dev, sendNum * sizeof(float), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(loc.y_hst, loc.y_dev, sendNum * sizeof(float), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(loc.z_hst, loc.z_dev, sendNum * sizeof(float), cudaMemcpyDeviceToHost));
  __NOTE__("rank %d: cudaMemcpy from device to host\n", mpi.rank);
#endif//MPI_VIA_HOST


  /** gather sampling points to the root process */
  chkMPIerr(MPI_Gather(&sendNum, 1, MPI_INT, sample.rnum, 1, MPI_INT, 0, mpi.comm));
  __NOTE__("rank %d: MPI_Gather\n", mpi.rank);

  /** set receive displacements (root process only) */
  int recvNum = 0;
  if( mpi.rank == 0 ){
    sample.disp[0] = 0;
    for(int jj = 1; jj < mpi.size; jj++)
      sample.disp[jj] = sample.disp[jj - 1] + sample.rnum[jj - 1];
    recvNum = sample.disp[mpi.size - 1] + sample.rnum[mpi.size - 1];
  }/* if( mpi.rank == 0 ){ */
  __NOTE__("rank %d: recvNum = %d\n", mpi.rank, recvNum);


  /** gather particle data to the root process */
#ifdef  MPI_ONE_SIDED_FOR_EXCG

  __NOTE__("MPI_Get start\n");
  if( mpi.rank == 0 ){
    chkMPIerr(MPI_Win_lock_all(0, loc.win_x));
    chkMPIerr(MPI_Win_lock_all(0, loc.win_y));
    chkMPIerr(MPI_Win_lock_all(0, loc.win_z));

    for(int ii = 0; ii < mpi.size; ii++){
      chkMPIerr(MPI_Get(ful.x_hst, sample.rnum[ii], MPI_FLOAT, ii, sample.disp[ii], sample.rnum[ii], MPI_FLOAT, loc.win_x));
      chkMPIerr(MPI_Get(ful.y_hst, sample.rnum[ii], MPI_FLOAT, ii, sample.disp[ii], sample.rnum[ii], MPI_FLOAT, loc.win_y));
      chkMPIerr(MPI_Get(ful.z_hst, sample.rnum[ii], MPI_FLOAT, ii, sample.disp[ii], sample.rnum[ii], MPI_FLOAT, loc.win_z));
    }/* for(int ii = 0; ii < mpi.size; ii++){ */

    chkMPIerr(MPI_Win_unlock_all(loc.win_x));
    chkMPIerr(MPI_Win_unlock_all(loc.win_y));
    chkMPIerr(MPI_Win_unlock_all(loc.win_z));
  }/* if( mpi.rank == 0 ){ */
  __NOTE__("MPI_Get finish\n");

#else///MPI_ONE_SIDED_FOR_EXCG

#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Gatherv(loc.x_dev, sendNum, MPI_FLOAT, ful.x_hst, sample.rnum, sample.disp, MPI_REALDAT, 0, mpi.comm));
  chkMPIerr(MPI_Gatherv(loc.y_dev, sendNum, MPI_FLOAT, ful.y_hst, sample.rnum, sample.disp, MPI_REALDAT, 0, mpi.comm));
  chkMPIerr(MPI_Gatherv(loc.z_dev, sendNum, MPI_FLOAT, ful.z_hst, sample.rnum, sample.disp, MPI_REALDAT, 0, mpi.comm));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Gatherv(loc.x_hst, sendNum, MPI_FLOAT, ful.x_hst, sample.rnum, sample.disp, MPI_REALDAT, 0, mpi.comm));
  chkMPIerr(MPI_Gatherv(loc.y_hst, sendNum, MPI_FLOAT, ful.y_hst, sample.rnum, sample.disp, MPI_REALDAT, 0, mpi.comm));
  chkMPIerr(MPI_Gatherv(loc.z_hst, sendNum, MPI_FLOAT, ful.z_hst, sample.rnum, sample.disp, MPI_REALDAT, 0, mpi.comm));
#endif//MPI_VIA_HOST

#endif//MPI_ONE_SIDED_FOR_EXCG


  /** set current (local) distribution */
  /** get (current) box size for the local distribution of N-body particles */
  getBoxSize_dev(numOld, (*src_dev).pos, soa, devProp, devInfo.stream[1]);
  float3 min, max;
  min.x = soa.min_hst->x;
  min.y = soa.min_hst->y;
  min.z = soa.min_hst->z;
  max.x = soa.max_hst->x;
  max.y = soa.max_hst->y;
  max.z = soa.max_hst->z;
  __NOTE__("xmin = %e, xmax = %e, ymin = %e, ymax = %e, zmin = %e, zmax = %e\n", min.x, max.x, min.y, max.y, min.z, max.z);

#ifdef  SHARE_PH_BOX_BOUNDARY
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, soa.min_hst, 3, MPI_FLOAT, MPI_MIN, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, soa.max_hst, 3, MPI_FLOAT, MPI_MAX, mpi.comm));
  checkCudaErrors(cudaMemcpyAsync(soa.box_min, soa.min_hst, sizeof(float4), cudaMemcpyHostToDevice, devInfo.stream[1]));
  checkCudaErrors(cudaMemcpyAsync(soa.box_max, soa.max_hst, sizeof(float4), cudaMemcpyHostToDevice, devInfo.stream[1]));
#endif//SHARE_PH_BOX_BOUNDARY

#ifdef  ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
  float diameter = 0.0f;
  diameter = fmaxf(diameter, max.x - min.x);
  diameter = fmaxf(diameter, max.y - min.y);
  diameter = fmaxf(diameter, max.z - min.z);
  diameter = ldexpf(1.0f, (int)ceilf(log2f(diameter)));
  soa.min_hst->x = 0.5f * (soa.min_hst->x + soa.max_hst->x - diameter);
  soa.min_hst->y = 0.5f * (soa.min_hst->y + soa.max_hst->y - diameter);
  soa.min_hst->z = 0.5f * (soa.min_hst->z + soa.max_hst->z - diameter);
  /* Nlevel ~ log_8(Ntot) = log_2(Ntot) / 3 */
  const float dL    = diameter * ldexpf(1.0f, -(((int)ceilf(log2f((float)Ntot) / 3.0f)) >> 1));
  const float dLinv = 1.0f / dL;
#endif//ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX


  /** determine local domain */
  float local_xmin = -0.5f * FLT_MAX;  float local_ymin = -0.5f * FLT_MAX;  float local_zmin = -0.5f * FLT_MAX;
  float local_xmax =  0.5f * FLT_MAX;  float local_ymax =  0.5f * FLT_MAX;  float local_zmax =  0.5f * FLT_MAX;

  /** domain decomposition in x-direction */
  if( mpi.dim[0] != 1 ){
    if( orm[0].rank == 0 ){
      __NOTE__("x-decomposition by mpi.rank = %d\n", mpi.rank);
      sendNum = recvNum;
      /** the root process determine the partition */
      if( rep[0].rank == 0 ){
	sort_xpos(recvNum, &ful, &loc, host);
	sample.xmin[0] = -0.5f * FLT_MAX;
	for(int ii = 0; ii < rep[0].size; ii++){
	  int Nini = (sendNum * (    ii)) / rep[0].size;
	  int Nfin = (sendNum * (1 + ii)) / rep[0].size;
	  sample.rnum[ii] = Nfin - Nini;

	  if( ii != (rep[0].size - 1) ){
	    const float middle = cutting(ful.x_hst[Nfin], ful.x_hst[Nfin + 1]
#ifdef  ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
					 , soa.min_hst->x, dL, dLinv
#endif//ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
					 );
	    sample.xmax[ii    ] = middle;
	    sample.xmin[ii + 1] = middle;
	    __NOTE__("sendNum = %d, size = %d, ii = %d, Nini = %d, Nfin = %d, middle = %e\n", sendNum, rep[0].size, ii, Nini, Nfin, middle);
	  }/* if( ii != (rep[0].size - 1) ){ */
	}/* for(int ii = 0; ii < rep[0].size; ii++){ */
	sample.xmax[rep[0].size - 1] = 0.5f * FLT_MAX;
      }/* if( rep[0].rank == 0 ){ */
      chkMPIerr(MPI_Scatter(sample.xmin, 1, MPI_FLOAT, &local_xmin, 1, MPI_FLOAT, 0, rep[0].comm));
      chkMPIerr(MPI_Scatter(sample.xmax, 1, MPI_FLOAT, &local_xmax, 1, MPI_FLOAT, 0, rep[0].comm));

      /** scatter sample particles if necessary */
      if( (mpi.dim[1] != 1) || (mpi.dim[2] != 1) ){
	/** set send displacements (root process only) */
	if( rep[0].rank == 0 ){
	  sample.disp[0] = 0;
	  for(int ii = 1; ii < rep[0].size; ii++)
	    sample.disp[ii] = sample.disp[ii - 1] + sample.rnum[ii - 1];
	}/* if( rep[0].rank == 0 ){ */

	/** scatter particle data from the root process */
	chkMPIerr(MPI_Scatter(sample.rnum, 1, MPI_INT, &recvNum, 1, MPI_INT, 0, rep[0].comm));
#ifndef OMIT_INPLACE_IN_SCATTERV
	chkMPIerr(MPI_Scatterv(ful.y_hst, sample.rnum, sample.disp, MPI_FLOAT, MPI_IN_PLACE, recvNum, MPI_FLOAT, 0, rep[0].comm));
	chkMPIerr(MPI_Scatterv(ful.z_hst, sample.rnum, sample.disp, MPI_FLOAT, MPI_IN_PLACE, recvNum, MPI_FLOAT, 0, rep[0].comm));
#else///OMIT_INPLACE_IN_SCATTERV
	chkMPIerr(MPI_Scatterv(ful.y_hst, sample.rnum, sample.disp, MPI_FLOAT, (rep[0].rank != 0) ? ful.y_hst : loc.y_hst, recvNum, MPI_FLOAT, 0, rep[0].comm));
	chkMPIerr(MPI_Scatterv(ful.z_hst, sample.rnum, sample.disp, MPI_FLOAT, (rep[0].rank != 0) ? ful.z_hst : loc.z_hst, recvNum, MPI_FLOAT, 0, rep[0].comm));
#endif//OMIT_INPLACE_IN_SCATTERV
      }/* if( (mpi.dim[1] != 1) || (mpi.dim[2] != 1) ){ */
    }/* if( orm[0].rank == 0 ){ */

    /** MPI_Bcast in orm[0].comm */
    if( (mpi.dim[1] != 1) || (mpi.dim[2] != 1) ){
      chkMPIerr(MPI_Bcast(&local_xmin, 1, MPI_FLOAT, 0, orm[0].comm));
      chkMPIerr(MPI_Bcast(&local_xmax, 1, MPI_FLOAT, 0, orm[0].comm));
    }/* if( (mpi.dim[1] != 1) || (mpi.dim[2] != 1) ){ */
  }/* if( mpi.dim[0] != 1 ){ */

  /** domain decomposition in y-direction */
  if( mpi.dim[1] != 1 ){
    if( orm[1].rank == 0 ){
      __NOTE__("y-decomposition by mpi.rank = %d\n", mpi.rank);
      sendNum = recvNum;
      /** the root process determine the partition */
      if( rep[1].rank == 0 ){
	sort_ypos(recvNum, &ful, host);
	sample.ymin[0] = -0.5f * FLT_MAX;
	for(int ii = 0; ii < rep[1].size; ii++){
	  int Nini = (sendNum * (    ii)) / rep[1].size;
	  int Nfin = (sendNum * (1 + ii)) / rep[1].size;
	  sample.rnum[ii] = Nfin - Nini;
	  if( ii != (rep[1].size - 1) ){
	    const float middle = cutting(ful.y_hst[Nfin], ful.y_hst[Nfin + 1]
#ifdef  ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
					 , soa.min_hst->y, dL, dLinv
#endif//ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
					 );
	    sample.ymax[ii    ] = middle;
	    sample.ymin[ii + 1] = middle;
	  }/* if( ii != (rep[1].size - 1) ){ */
	}/* for(int ii = 0; ii < rep[1].size; ii++){ */
	sample.ymax[rep[1].size - 1] = 0.5f * FLT_MAX;
      }/* if( rep[1].rank == 0 ){ */
      chkMPIerr(MPI_Scatter(sample.ymin, 1, MPI_FLOAT, &local_ymin, 1, MPI_FLOAT, 0, rep[1].comm));
      chkMPIerr(MPI_Scatter(sample.ymax, 1, MPI_FLOAT, &local_ymax, 1, MPI_FLOAT, 0, rep[1].comm));

      /** scatter sample particles if necessary */
      if( mpi.dim[2] != 1 ){
	/** set send displacements (root process only) */
	if( rep[1].rank == 0 ){
	  sample.disp[0] = 0;
	  for(int ii = 1; ii < rep[1].size; ii++)
	    sample.disp[ii] = sample.disp[ii - 1] + sample.rnum[ii - 1];
	}/* if( rep[1].rank == 0 ){ */

	/** scatter particle data from the root process */
	chkMPIerr(MPI_Scatter(sample.rnum, 1, MPI_INT, &recvNum, 1, MPI_INT, 0, rep[1].comm));
#ifndef OMIT_INPLACE_IN_SCATTERV
	chkMPIerr(MPI_Scatterv(ful.z_hst, sample.rnum, sample.disp, MPI_FLOAT, MPI_IN_PLACE, recvNum, MPI_FLOAT, 0, rep[1].comm));
#else///OMIT_INPLACE_IN_SCATTERV
	chkMPIerr(MPI_Scatterv(ful.z_hst, sample.rnum, sample.disp, MPI_FLOAT, (rep[1].rank != 0) ? ful.z_hst : loc.z_hst, recvNum, MPI_FLOAT, 0, rep[1].comm));
#endif//OMIT_INPLACE_IN_SCATTERV
      }/* if( mpi.dim[2] != 1 ){ */
    }/* if( orm[1].rank == 0 ){ */

    /** MPI_Bcast in orm[1].comm */
    if( mpi.dim[2] != 1 ){
      chkMPIerr(MPI_Bcast(&local_ymin, 1, MPI_FLOAT, 0, orm[1].comm));
      chkMPIerr(MPI_Bcast(&local_ymax, 1, MPI_FLOAT, 0, orm[1].comm));
    }/* if( mpi.dim[2] != 1 ){ */
  }/* if( mpi.dim[1] != 1 ){ */

  /** domain decomposition in z-direction */
  if( mpi.dim[2] != 1 ){
    if( orm[2].rank == 0 ){
      __NOTE__("z-decomposition by mpi.rank = %d\n", mpi.rank);
      sendNum = recvNum;
      /** the root process determine the partition */
      if( rep[2].rank == 0 ){
	sort_zpos(recvNum, &ful, host);
	sample.zmin[0] = -0.5f * FLT_MAX;
	for(int ii = 0; ii < rep[2].size; ii++){
	  int Nini = (sendNum * (    ii)) / rep[2].size;
	  int Nfin = (sendNum * (1 + ii)) / rep[2].size;
	  sample.rnum[ii] = Nfin - Nini;

	  if( ii != (rep[2].size - 1) ){
	    const float middle = cutting(ful.z_hst[Nfin], ful.z_hst[Nfin + 1]
#ifdef  ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
					 , soa.min_hst->z, dL, dLinv
#endif//ALIGN_DOMAIN_BOUNDARY_TO_PH_BOX
					 );
	    sample.zmax[ii    ] = middle;
	    sample.zmin[ii + 1] = middle;
	  }/* if( ii != (rep[2].size - 1) ){ */
	}/* for(int ii = 0; ii < rep[2].size; ii++){ */
	sample.zmax[rep[2].size - 1] = 0.5f * FLT_MAX;
      }/* if( rep[2].rank == 0 ){ */
      chkMPIerr(MPI_Scatter(sample.zmin, 1, MPI_FLOAT, &local_zmin, 1, MPI_FLOAT, 0, rep[2].comm));
      chkMPIerr(MPI_Scatter(sample.zmax, 1, MPI_FLOAT, &local_zmax, 1, MPI_FLOAT, 0, rep[2].comm));
    }/* if( orm[2].rank == 0 ){ */
  }/* if( mpi.dim[2] != 1 ){ */


  /** share the decomposed domain */
  chkMPIerr(MPI_Allgather(&local_xmin, 1, MPI_FLOAT, domain.xmin, 1, MPI_FLOAT, mpi.comm));
  chkMPIerr(MPI_Allgather(&local_xmax, 1, MPI_FLOAT, domain.xmax, 1, MPI_FLOAT, mpi.comm));
  chkMPIerr(MPI_Allgather(&local_ymin, 1, MPI_FLOAT, domain.ymin, 1, MPI_FLOAT, mpi.comm));
  chkMPIerr(MPI_Allgather(&local_ymax, 1, MPI_FLOAT, domain.ymax, 1, MPI_FLOAT, mpi.comm));
  chkMPIerr(MPI_Allgather(&local_zmin, 1, MPI_FLOAT, domain.zmin, 1, MPI_FLOAT, mpi.comm));
  chkMPIerr(MPI_Allgather(&local_zmax, 1, MPI_FLOAT, domain.zmax, 1, MPI_FLOAT, mpi.comm));
  __NOTE__("[%e, %e]x[%e, %e]x[%e, %e]\n",
	   domain.xmin[mpi.rank], domain.xmax[mpi.rank],
	   domain.ymin[mpi.rank], domain.ymax[mpi.rank],
	   domain.zmin[mpi.rank], domain.zmax[mpi.rank]);

#ifdef  CARE_EXTERNAL_PARTICLES
  location->boxmin.x = local_xmin;
  location->boxmin.y = local_ymin;
  location->boxmin.z = local_zmin;
  location->boxmax.x = local_xmax;
  location->boxmax.y = local_ymax;
  location->boxmax.z = local_zmax;
#ifdef  TIME_BASED_MODIFICATION
  location->elapsed = 0.0f;
  location->dtmin   = 0.5f * FLT_MAX;
  location->dtinv   = 1.0f / location->dtmin;
#else///TIME_BASED_MODIFICATION
  location->step = 0;
#endif//TIME_BASED_MODIFICATION
#endif//CARE_EXTERNAL_PARTICLES


  /** exchange N-body particles */
  const int numProcs = mpi.size;
  for(int ii = 0; ii < numProcs; ii++){
#ifdef  MPI_ONE_SIDED_FOR_EXCG
    chkMPIerr(MPI_Irecv(&(recvBuf[ii].body), 1, mpi.send, ii, MPI_TAG_BODY(ii, mpi.size), mpi.comm, &(recvBuf[ii].req)));
#else///MPI_ONE_SIDED_FOR_EXCG
    chkMPIerr(MPI_Irecv(&(recvBuf[ii].num ), 1, MPI_INT , ii, MPI_TAG_BODY(ii, mpi.size), mpi.comm, &(recvBuf[ii].req)));
#endif//MPI_ONE_SIDED_FOR_EXCG
  }/* for(int ii = 0; ii < numProcs; ii++){ */

  int overlapNum = 0;
#ifdef  MPI_ONE_SIDED_FOR_EXCG
  const sendBody nullSend = {0};
#else///MPI_ONE_SIDED_FOR_EXCG
  const int nullSend = 0;
#endif//MPI_ONE_SIDED_FOR_EXCG
  /* this function should be performed on host, use pinned memory */
  for(int ii = 0; ii < numProcs; ii++){
    if( (min.x <= domain.xmax[ii]) && (max.x >= domain.xmin[ii]) &&
	(min.y <= domain.ymax[ii]) && (max.y >= domain.ymin[ii]) &&
	(min.z <= domain.zmax[ii]) && (max.z >= domain.zmin[ii]) ){
      /** if spatial overlap is detected, ... */
      sendBuf[overlapNum].rank = ii;
      domBoundary.rank_hst[overlapNum] = ii;

      domBoundary.xmin_hst[overlapNum] = (domain.xmin[ii] < -0.25f * FLT_MAX) ? (domain.xmin[ii]) : (fmaxf(min.x, domain.xmin[ii]));
      domBoundary.ymin_hst[overlapNum] = (domain.ymin[ii] < -0.25f * FLT_MAX) ? (domain.ymin[ii]) : (fmaxf(min.y, domain.ymin[ii]));
      domBoundary.zmin_hst[overlapNum] = (domain.zmin[ii] < -0.25f * FLT_MAX) ? (domain.zmin[ii]) : (fmaxf(min.z, domain.zmin[ii]));

      domBoundary.xmax_hst[overlapNum] = (domain.xmax[ii] >  0.25f * FLT_MAX) ? (domain.xmax[ii]) : (fminf(max.x, domain.xmax[ii]));
      domBoundary.ymax_hst[overlapNum] = (domain.ymax[ii] >  0.25f * FLT_MAX) ? (domain.ymax[ii]) : (fminf(max.y, domain.ymax[ii]));
      domBoundary.zmax_hst[overlapNum] = (domain.zmax[ii] >  0.25f * FLT_MAX) ? (domain.zmax[ii]) : (fminf(max.z, domain.zmax[ii]));

      overlapNum++;
    }
    else{      /* if covered areas do not overlap, ... */
      __NOTE__("nullSend for rank %d from rank %d\n", ii, mpi.rank);
#ifdef  MPI_ONE_SIDED_FOR_EXCG
      chkMPIerr(MPI_Isend(&nullSend, 1, mpi.send, ii, MPI_TAG_BODY(mpi.rank, mpi.size), mpi.comm, &(domain.req[ii])));
#else///MPI_ONE_SIDED_FOR_EXCG
      chkMPIerr(MPI_Isend(&nullSend, 1, MPI_INT , ii, MPI_TAG_BODY(mpi.rank, mpi.size), mpi.comm, &(domain.req[ii])));
#endif//MPI_ONE_SIDED_FOR_EXCG
    }/* else{ */
  }/* for(int ii = 0; ii < numProcs; ii++){ */


  /** send domBoundary from host to device */
  checkCudaErrors(cudaMemcpy(domBoundary.rank    , domBoundary.rank_hst, sizeof(int)   * overlapNum, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(domBoundary.xmin_dev, domBoundary.xmin_hst, sizeof(float) * overlapNum, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(domBoundary.xmax_dev, domBoundary.xmax_hst, sizeof(float) * overlapNum, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(domBoundary.ymin_dev, domBoundary.ymin_hst, sizeof(float) * overlapNum, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(domBoundary.ymax_dev, domBoundary.ymax_hst, sizeof(float) * overlapNum, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(domBoundary.zmin_dev, domBoundary.zmin_hst, sizeof(float) * overlapNum, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(domBoundary.zmax_dev, domBoundary.zmax_hst, sizeof(float) * overlapNum, cudaMemcpyHostToDevice));


  /** determine process rank for each particle to belong */
  initDstRank<<<BLOCKSIZE(numOld, 1024), 1024>>>(numOld, key.dstRank_dev);

#ifdef  USE_OCCUPANCY_CALCULATOR
  const int NBLOCKS_PER_SM_ASSIGN = domBoundary.numBlocksPerSM_assign;
#endif//USE_OCCUPANCY_CALCULATOR
  for(int ii = 0; ii < overlapNum; ii += 4)
    assignNewDomain_kernel<<<devProp.numSM * NBLOCKS_PER_SM_ASSIGN, NTHREADS_ASSIGN>>>(numOld, domBoundary.numNew, (*src_dev).pos, ii, overlapNum - ii, domBoundary.rank, key.dstRank_dev, domBoundary.xmin_dev, domBoundary.xmax_dev, domBoundary.ymin_dev, domBoundary.ymax_dev, domBoundary.zmin_dev, domBoundary.zmax_dev, domBoundary.gmem, domBoundary.gsync0, domBoundary.gsync1);
  getLastCudaError("assignNewDomain_kernel");

  checkCudaErrors(cudaMemcpy(domBoundary.numNew_hst, domBoundary.numNew, overlapNum * sizeof(int), cudaMemcpyDeviceToHost));


#ifdef  MPI_ONE_SIDED_FOR_EXCG
  sendBuf[0].body.head = 0;
#else///MPI_ONE_SIDED_FOR_EXCG
  sendBuf[0].     head = 0;
#endif//MPI_ONE_SIDED_FOR_EXCG
  for(int ii = 0; ii < overlapNum; ii++){
#ifdef  MPI_ONE_SIDED_FOR_EXCG
    sendBuf[ii].body.num = domBoundary.numNew_hst[ii];
#else///MPI_ONE_SIDED_FOR_EXCG
    sendBuf[ii].     num = domBoundary.numNew_hst[ii];
#endif//MPI_ONE_SIDED_FOR_EXCG

    if( ii > 0 ){
#ifdef  MPI_ONE_SIDED_FOR_EXCG
      sendBuf[ii].body.head = sendBuf[ii - 1].body.head + sendBuf[ii - 1].body.num;
#else///MPI_ONE_SIDED_FOR_EXCG
      sendBuf[ii].     head = sendBuf[ii - 1].     head + sendBuf[ii - 1].     num;
#endif//MPI_ONE_SIDED_FOR_EXCG
    }/* if( ii > 0 ){ */

#ifdef  MPI_ONE_SIDED_FOR_EXCG
    chkMPIerr(MPI_Isend(&(sendBuf[ii].body), 1, mpi.send, sendBuf[ii].rank, MPI_TAG_BODY(mpi.rank, mpi.size), mpi.comm, &(domain.req[sendBuf[ii].rank])));
#else///MPI_ONE_SIDED_FOR_EXCG
    chkMPIerr(MPI_Isend(&(sendBuf[ii].num ), 1, MPI_INT , sendBuf[ii].rank, MPI_TAG_BODY(mpi.rank, mpi.size), mpi.comm, &(domain.req[sendBuf[ii].rank])));
#endif//MPI_ONE_SIDED_FOR_EXCG
  }/* for(int ii = 0; ii < overlapNum; ii++){ */

#ifdef  MPI_ONE_SIDED_FOR_EXCG
  if( (sendBuf[overlapNum - 1].body.head + sendBuf[overlapNum - 1].body.num) != numOld ){
    __KILL__(stderr, "ERROR: total number of scattered particles (%d) is differ from that of local particles (%d)\n", sendBuf[overlapNum - 1].body.head + sendBuf[overlapNum - 1].body.num, numOld);
  }/* if( (sendBuf[overlapNum - 1].body.head + sendBuf[overlapNum - 1].body.num) != numOld ){ */
#else///MPI_ONE_SIDED_FOR_EXCG
  if( (sendBuf[overlapNum - 1].head + sendBuf[overlapNum - 1].num) != numOld ){
    __KILL__(stderr, "ERROR: total number of scattered particles (%d) is differ from that of local particles (%d)\n", sendBuf[overlapNum - 1].head + sendBuf[overlapNum - 1].num, numOld);
  }/* if( (sendBuf[overlapNum - 1].head + sendBuf[overlapNum - 1].num) != numOld ){ */
#endif//MPI_ONE_SIDED_FOR_EXCG

  sortDomainDecomposeKey(numOld, key, src_dev, dst_dev);


  /** send particle data */
#ifdef  MPI_VIA_HOST
  copyParticle_dev2hst(numOld, *src_dev, *src_hst
#ifdef  EXEC_BENCHMARK
		       , elapsed
#endif//EXEC_BENCHMARK
		       );
#endif//MPI_VIA_HOST


#ifdef  MPI_ONE_SIDED_FOR_EXCG
  /** receive particle data */
#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_lock_all(0, (*src_dev).win_ipos));
  chkMPIerr(MPI_Win_lock_all(0, (*src_dev).win_idx));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_lock_all(0, (*src_hst).win_ipos));
  chkMPIerr(MPI_Win_lock_all(0, (*src_hst).win_idx));
#endif//MPI_VIA_HOST

#ifdef  BLOCK_TIME_STEP
#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_lock_all(0, (*src_dev).win_ivel));
  chkMPIerr(MPI_Win_lock_all(0, (*src_dev).win_time));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_lock_all(0, (*src_hst).win_ivel));
  chkMPIerr(MPI_Win_lock_all(0, (*src_hst).win_time));
#endif//MPI_VIA_HOST
#else///BLOCK_TIME_STEP
#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_lock_all(0, (*src_dev).win_vx));
  chkMPIerr(MPI_Win_lock_all(0, (*src_dev).win_vy));
  chkMPIerr(MPI_Win_lock_all(0, (*src_dev).win_vz));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_lock_all(0, (*src_hst).win_vx));
  chkMPIerr(MPI_Win_lock_all(0, (*src_hst).win_vy));
  chkMPIerr(MPI_Win_lock_all(0, (*src_hst).win_vz));
#endif//MPI_VIA_HOST
#endif//BLOCK_TIME_STEP

#ifdef  GADGET_MAC
#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_lock_all(0, (*src_dev).win_iacc));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_lock_all(0, (*src_hst).win_iacc));
#endif//MPI_VIA_HOST
#endif//GADGET_MAC

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_lock_all(0, (*src_dev).win_iext));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_lock_all(0, (*src_hst).win_iext));
#endif//MPI_VIA_HOST
#endif//SET_EXTERNAL_POTENTIAL_FIELD


  *numNew = 0;
  int recvHead = 0;
  for(int ii = 0; ii < numProcs; ii++){
    /** receive recvNum */
    MPI_Status status;
    chkMPIerr(MPI_Wait(&(recvBuf[ii].req), &status));

    /** if recvNum > 0, then set receive buffer */
    if( recvBuf[ii].body.num > 0 ){
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Get(&((*dst_dev).pos[recvHead]), recvBuf[ii].body.num, mpi.ipos         , ii, recvBuf[ii].body.head, recvBuf[ii].body.num, mpi.ipos         , (*src_dev).win_ipos));
      chkMPIerr(MPI_Get(&((*dst_dev).idx[recvHead]), recvBuf[ii].body.num, MPI_UNSIGNED_LONG, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, MPI_UNSIGNED_LONG, (*src_dev).win_idx));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Get(&((*dst_hst).pos[recvHead]), recvBuf[ii].body.num, mpi.ipos         , ii, recvBuf[ii].body.head, recvBuf[ii].body.num, mpi.ipos         , (*src_hst).win_ipos));
      chkMPIerr(MPI_Get(&((*dst_hst).idx[recvHead]), recvBuf[ii].body.num, MPI_UNSIGNED_LONG, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, MPI_UNSIGNED_LONG, (*src_hst).win_idx));
#endif//MPI_VIA_HOST

#ifdef  BLOCK_TIME_STEP
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Get(&((*dst_dev).vel [recvHead]), recvBuf[ii].body.num, mpi.ivel, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, mpi.ivel, (*src_dev).win_ivel));
      chkMPIerr(MPI_Get(&((*dst_dev).time[recvHead]), recvBuf[ii].body.num, mpi.time, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, mpi.time, (*src_dev).win_time));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Get(&((*dst_hst).vel [recvHead]), recvBuf[ii].body.num, mpi.ivel, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, mpi.ivel, (*src_hst).win_ivel));
      chkMPIerr(MPI_Get(&((*dst_hst).time[recvHead]), recvBuf[ii].body.num, mpi.time, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, mpi.time, (*src_hst).win_time));
#endif//MPI_VIA_HOST
#else///BLOCK_TIME_STEP
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Get(&((*dst_dev).vx[recvHead]), recvBuf[ii].body.num, MPI_REALDAT, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, MPI_REALDAT, (*src_dev).win_vx));
      chkMPIerr(MPI_Get(&((*dst_dev).vy[recvHead]), recvBuf[ii].body.num, MPI_REALDAT, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, MPI_REALDAT, (*src_dev).win_vy));
      chkMPIerr(MPI_Get(&((*dst_dev).vz[recvHead]), recvBuf[ii].body.num, MPI_REALDAT, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, MPI_REALDAT, (*src_dev).win_vz));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Get(&((*dst_hst).vx[recvHead]), recvBuf[ii].body.num, MPI_REALDAT, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, MPI_REALDAT, (*src_hst).win_vx));
      chkMPIerr(MPI_Get(&((*dst_hst).vy[recvHead]), recvBuf[ii].body.num, MPI_REALDAT, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, MPI_REALDAT, (*src_hst).win_vy));
      chkMPIerr(MPI_Get(&((*dst_hst).vz[recvHead]), recvBuf[ii].body.num, MPI_REALDAT, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, MPI_REALDAT, (*src_hst).win_vz));
#endif//MPI_VIA_HOST
#endif//BLOCK_TIME_STEP

#ifdef  GADGET_MAC
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Get(&((*dst_dev).acc[recvHead]), recvBuf[ii].body.num, mpi.iacc, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, mpi.iacc, (*src_dev).win_iacc));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Get(&((*dst_hst).acc[recvHead]), recvBuf[ii].body.num, mpi.iacc, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, mpi.iacc, (*src_hst).win_iacc));
#endif//MPI_VIA_HOST
#endif//GADGET_MAC

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Get(&((*dst_dev).acc_ext[recvHead]), recvBuf[ii].body.num, mpi.iacc, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, mpi.iacc, (*src_dev).win_iext));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Get(&((*dst_hst).acc_ext[recvHead]), recvBuf[ii].body.num, mpi.iacc, ii, recvBuf[ii].body.head, recvBuf[ii].body.num, mpi.iacc, (*src_hst).win_iext));
#endif//MPI_VIA_HOST
#endif//SET_EXTERNAL_POTENTIAL_FIELD


      *numNew  += recvBuf[ii].body.num;
      recvHead += recvBuf[ii].body.num;
    }/* if( recvBuf[ii].body.num > 0 ){ */
  }/* for(int ii = 0; ii < numProcs; ii++){ */


  /* complete MPI communications */
  for(int ii = 0; ii < numProcs; ii++)
    if( recvBuf[ii].body.num > 0 ){
      __NOTE__("MPI_Win_flush_local from rank %d (%d particles)\n", ii, recvBuf[ii].body.num);
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Win_flush_local(ii, (*src_dev).win_ipos));
      chkMPIerr(MPI_Win_flush_local(ii, (*src_dev).win_idx));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Win_flush_local(ii, (*src_hst).win_ipos));
      chkMPIerr(MPI_Win_flush_local(ii, (*src_hst).win_idx));
#endif//MPI_VIA_HOST

#ifdef  BLOCK_TIME_STEP
#ifdef  MPI_VIA_HOST
      chkMPIerr(MPI_Win_flush_local(ii, (*src_dev).win_ivel));
      chkMPIerr(MPI_Win_flush_local(ii, (*src_dev).win_time));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Win_flush_local(ii, (*src_hst).win_ivel));
      chkMPIerr(MPI_Win_flush_local(ii, (*src_hst).win_time));
#endif//MPI_VIA_HOST
#else///BLOCK_TIME_STEP
#ifdef  MPI_VIA_HOST
      chkMPIerr(MPI_Win_flush_local(ii, (*src_dev).win_vx));
      chkMPIerr(MPI_Win_flush_local(ii, (*src_dev).win_vy));
      chkMPIerr(MPI_Win_flush_local(ii, (*src_dev).win_vz));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Win_flush_local(ii, (*src_hst).win_vx));
      chkMPIerr(MPI_Win_flush_local(ii, (*src_hst).win_vy));
      chkMPIerr(MPI_Win_flush_local(ii, (*src_hst).win_vz));
#endif//MPI_VIA_HOST
#endif//BLOCK_TIME_STEP

#ifdef  GADGET_MAC
#ifdef  MPI_VIA_HOST
      chkMPIerr(MPI_Win_flush_local(ii, (*src_dev).win_iacc));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Win_flush_local(ii, (*src_hst).win_iacc));
#endif//MPI_VIA_HOST
#endif//GADGET_MAC

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  MPI_VIA_HOST
      chkMPIerr(MPI_Win_flush_local(ii, (*src_dev).win_iext));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Win_flush_local(ii, (*src_hst).win_iext));
#endif//MPI_VIA_HOST
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* if( recvBuf[ii].body.num > 0 ){ */


#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_unlock_all((*src_dev).win_ipos));
  chkMPIerr(MPI_Win_unlock_all((*src_dev).win_idx));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_unlock_all((*src_hst).win_ipos));
  chkMPIerr(MPI_Win_unlock_all((*src_hst).win_idx));
#endif//MPI_VIA_HOST

#ifdef  BLOCK_TIME_STEP
#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_unlock_all((*src_dev).win_ivel));
  chkMPIerr(MPI_Win_unlock_all((*src_dev).win_time));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_unlock_all((*src_hst).win_ivel));
  chkMPIerr(MPI_Win_unlock_all((*src_hst).win_time));
#endif//MPI_VIA_HOST
#else///BLOCK_TIME_STEP
#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_unlock_all((*src_dev).win_vx));
  chkMPIerr(MPI_Win_unlock_all((*src_dev).win_vy));
  chkMPIerr(MPI_Win_unlock_all((*src_dev).win_vz));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_unlock_all((*src_hst).win_vx));
  chkMPIerr(MPI_Win_unlock_all((*src_hst).win_vy));
  chkMPIerr(MPI_Win_unlock_all((*src_hst).win_vz));
#endif//MPI_VIA_HOST
#endif//BLOCK_TIME_STEP

#ifdef  GADGET_MAC
#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_unlock_all((*src_dev).win_iacc));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_unlock_all((*src_hst).win_iacc));
#endif//MPI_VIA_HOST
#endif//GADGET_MAC

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
#ifndef MPI_VIA_HOST
  chkMPIerr(MPI_Win_unlock_all((*src_dev).win_iext));
#else///MPI_VIA_HOST
  chkMPIerr(MPI_Win_unlock_all((*src_hst).win_iext));
#endif//MPI_VIA_HOST
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#else///MPI_ONE_SIDED_FOR_EXCG

  /** receive particle data */
  *numNew = 0;
  int recvHead = 0;
  for(int ii = 0; ii < numProcs; ii++){
    /** receive recvNum */
    MPI_Status status;
    chkMPIerr(MPI_Wait(&(recvBuf[ii].req), &status));
    __NOTE__("recvBuf[%d].req: head = %d, num = %d, src = rank %d\n", ii, recvHead, recvBuf[ii].num, ii);

    /** if recvNum > 0, then set receive buffer */
    if( recvBuf[ii].num > 0 ){
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Irecv(&((*dst_dev).pos[recvHead]), recvBuf[ii].num, mpi.ipos         , ii, MPI_TAG_IPOS(ii, mpi.size), mpi.comm, &(recvBuf[ii].pos)));
      chkMPIerr(MPI_Irecv(&((*dst_dev).idx[recvHead]), recvBuf[ii].num, MPI_UNSIGNED_LONG, ii, MPI_TAG_IIDX(ii, mpi.size), mpi.comm, &(recvBuf[ii].idx)));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Irecv(&((*dst_hst).pos[recvHead]), recvBuf[ii].num, mpi.ipos         , ii, MPI_TAG_IPOS(ii, mpi.size), mpi.comm, &(recvBuf[ii].pos)));
      chkMPIerr(MPI_Irecv(&((*dst_hst).idx[recvHead]), recvBuf[ii].num, MPI_UNSIGNED_LONG, ii, MPI_TAG_IIDX(ii, mpi.size), mpi.comm, &(recvBuf[ii].idx)));
#endif//MPI_VIA_HOST

#ifdef  BLOCK_TIME_STEP
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Irecv(&((*dst_dev).vel [recvHead]), recvBuf[ii].num, mpi.ivel, ii, MPI_TAG_IVEL(ii, mpi.size), mpi.comm, &(recvBuf[ii].vel)));
      chkMPIerr(MPI_Irecv(&((*dst_dev).time[recvHead]), recvBuf[ii].num, mpi.time, ii, MPI_TAG_TIME(ii, mpi.size), mpi.comm, &(recvBuf[ii].time)));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Irecv(&((*dst_hst).vel [recvHead]), recvBuf[ii].num, mpi.ivel, ii, MPI_TAG_IVEL(ii, mpi.size), mpi.comm, &(recvBuf[ii].vel)));
      chkMPIerr(MPI_Irecv(&((*dst_hst).time[recvHead]), recvBuf[ii].num, mpi.time, ii, MPI_TAG_TIME(ii, mpi.size), mpi.comm, &(recvBuf[ii].time)));
#endif//MPI_VIA_HOST
#else///BLOCK_TIME_STEP
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Irecv(&((*dst_dev).vx[recvHead]), recvBuf[ii].num, MPI_REALDAT, ii, MPI_TAG_VELX(ii, mpi.size), mpi.comm, &(recvBuf[ii].vx )));
      chkMPIerr(MPI_Irecv(&((*dst_dev).vy[recvHead]), recvBuf[ii].num, MPI_REALDAT, ii, MPI_TAG_VELY(ii, mpi.size), mpi.comm, &(recvBuf[ii].vy )));
      chkMPIerr(MPI_Irecv(&((*dst_dev).vz[recvHead]), recvBuf[ii].num, MPI_REALDAT, ii, MPI_TAG_VELZ(ii, mpi.size), mpi.comm, &(recvBuf[ii].vz )));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Irecv(&((*dst_hst).vx[recvHead]), recvBuf[ii].num, MPI_REALDAT, ii, MPI_TAG_VELX(ii, mpi.size), mpi.comm, &(recvBuf[ii].vx)));
      chkMPIerr(MPI_Irecv(&((*dst_hst).vy[recvHead]), recvBuf[ii].num, MPI_REALDAT, ii, MPI_TAG_VELY(ii, mpi.size), mpi.comm, &(recvBuf[ii].vy)));
      chkMPIerr(MPI_Irecv(&((*dst_hst).vz[recvHead]), recvBuf[ii].num, MPI_REALDAT, ii, MPI_TAG_VELZ(ii, mpi.size), mpi.comm, &(recvBuf[ii].vz)));
#endif//MPI_VIA_HOST
#endif//BLOCK_TIME_STEP

#ifdef  GADGET_MAC
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Irecv(&((*dst_dev).acc[recvHead]), recvBuf[ii].num, mpi.iacc, ii, MPI_TAG_IACC(ii, mpi.size), mpi.comm, &(recvBuf[ii].acc)));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Irecv(&((*dst_hst).acc[recvHead]), recvBuf[ii].num, mpi.iacc, ii, MPI_TAG_IACC(ii, mpi.size), mpi.comm, &(recvBuf[ii].acc)));
#endif//MPI_VIA_HOST
#endif//GADGET_MAC

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Irecv(&((*dst_dev).acc_ext[recvHead]), recvBuf[ii].num, mpi.iacc, ii, MPI_TAG_IEXT(ii, mpi.size), mpi.comm, &(recvBuf[ii].ext)));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Irecv(&((*dst_hst).acc_ext[recvHead]), recvBuf[ii].num, mpi.iacc, ii, MPI_TAG_IEXT(ii, mpi.size), mpi.comm, &(recvBuf[ii].ext)));
#endif//MPI_VIA_HOST
#endif//SET_EXTERNAL_POTENTIAL_FIELD

      *numNew  += recvBuf[ii].num;
      recvHead += recvBuf[ii].num;
    }/* if( recvBuf[ii].num > 0 ){ */
  }/* for(int ii = 0; ii < numProcs; ii++){ */

  /** send particle data */
  for(int ii = 0; ii < overlapNum; ii++){
    /** if sendBuf[ii].num > 0, then send N-body particles */
    if( sendBuf[ii].num > 0 ){
      __NOTE__("sendBuf[%d]: head = %d, num = %d, dst = rank %d\n", ii, sendBuf[ii].head, sendBuf[ii].num, sendBuf[ii].rank);

#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Isend(&((*src_dev).pos[sendBuf[ii].head]), sendBuf[ii].num, mpi.ipos         , sendBuf[ii].rank, MPI_TAG_IPOS(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].pos)));
      chkMPIerr(MPI_Isend(&((*src_dev).idx[sendBuf[ii].head]), sendBuf[ii].num, MPI_UNSIGNED_LONG, sendBuf[ii].rank, MPI_TAG_IIDX(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].idx)));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Isend(&((*src_hst).pos[sendBuf[ii].head]), sendBuf[ii].num, mpi.ipos         , sendBuf[ii].rank, MPI_TAG_IPOS(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].pos)));
      chkMPIerr(MPI_Isend(&((*src_hst).idx[sendBuf[ii].head]), sendBuf[ii].num, MPI_UNSIGNED_LONG, sendBuf[ii].rank, MPI_TAG_IIDX(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].idx)));
#endif//MPI_VIA_HOST

#ifdef  BLOCK_TIME_STEP
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Isend(&((*src_dev).vel [sendBuf[ii].head]), sendBuf[ii].num, mpi.ivel, sendBuf[ii].rank, MPI_TAG_IVEL(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].vel)));
      chkMPIerr(MPI_Isend(&((*src_dev).time[sendBuf[ii].head]), sendBuf[ii].num, mpi.time, sendBuf[ii].rank, MPI_TAG_TIME(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].time)));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Isend(&((*src_hst).vel [sendBuf[ii].head]), sendBuf[ii].num, mpi.ivel, sendBuf[ii].rank, MPI_TAG_IVEL(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].vel)));
      chkMPIerr(MPI_Isend(&((*src_hst).time[sendBuf[ii].head]), sendBuf[ii].num, mpi.time, sendBuf[ii].rank, MPI_TAG_TIME(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].time)));
#endif//MPI_VIA_HOST
#else///BLOCK_TIME_STEP
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Isend(&((*src_dev).vx[sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT, sendBuf[ii].rank, MPI_TAG_VELX(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].vx)));
      chkMPIerr(MPI_Isend(&((*src_dev).vy[sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT, sendBuf[ii].rank, MPI_TAG_VELY(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].vy)));
      chkMPIerr(MPI_Isend(&((*src_dev).vz[sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT, sendBuf[ii].rank, MPI_TAG_VELZ(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].vz)));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Isend(&((*src_hst).vx[sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT, sendBuf[ii].rank, MPI_TAG_VELX(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].vx)));
      chkMPIerr(MPI_Isend(&((*src_hst).vy[sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT, sendBuf[ii].rank, MPI_TAG_VELY(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].vy)));
      chkMPIerr(MPI_Isend(&((*src_hst).vz[sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT, sendBuf[ii].rank, MPI_TAG_VELZ(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].vz)));
#endif//MPI_VIA_HOST
#endif//BLOCK_TIME_STEP

#ifdef  GADGET_MAC
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Isend(&((*src_dev).acc[sendBuf[ii].head]), sendBuf[ii].num, mpi.iacc, sendBuf[ii].rank, MPI_TAG_IACC(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].acc)));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Isend(&((*src_hst).acc[sendBuf[ii].head]), sendBuf[ii].num, mpi.iacc, sendBuf[ii].rank, MPI_TAG_IACC(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].acc)));
#endif//MPI_VIA_HOST
#endif//GADGET_MAC

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
#ifndef MPI_VIA_HOST
      chkMPIerr(MPI_Isend(&((*src_dev).acc_ext[sendBuf[ii].head]), sendBuf[ii].num, mpi.iacc, sendBuf[ii].rank, MPI_TAG_IEXT(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].ext)));
#else///MPI_VIA_HOST
      chkMPIerr(MPI_Isend(&((*src_hst).acc_ext[sendBuf[ii].head]), sendBuf[ii].num, mpi.iacc, sendBuf[ii].rank, MPI_TAG_IEXT(mpi.rank, mpi.size), mpi.comm, &(sendBuf[ii].ext)));
#endif//MPI_VIA_HOST
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#ifndef NDEBUG
#ifdef  MPI_VIA_HOST
      for(int jj = 0; jj < sendBuf[ii].num; jj++){
	const int kk = sendBuf[ii].head + jj;
	if( (*src_hst).pos[kk].m < EPSILON )
	  __NOTE__("pi[%d]: sendmass = %e, x = %e, y = %e, z = %e\n", kk, (*src_hst).pos[kk].m, (*src_hst).pos[kk].x, (*src_hst).pos[kk].y, (*src_hst).pos[kk].z);
      }
#endif//MPI_VIA_HOST
#endif//NDEBUG
    }/* if( sendBuf[ii].num > 0 ){ */
  }/* for(int ii = 0; ii < overlapNum; ii++){ */


  /** complete MPI communications */
  for(int ii = 0; ii < overlapNum; ii++){
    if( sendBuf[ii].num > 0 ){
      __NOTE__("send %d bodies to rank %d (%d/%d)\n", sendBuf[ii].num, sendBuf[ii].rank, ii, overlapNum);
      MPI_Status  pos;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. pos), &pos));    __NOTE__("sendBuf[%d].pos\n", ii);
      MPI_Status  idx;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. idx), &idx));    __NOTE__("sendBuf[%d].idx\n", ii);
#ifdef  BLOCK_TIME_STEP
      MPI_Status  vel;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. vel), &vel));    __NOTE__("sendBuf[%d].vel\n", ii);
      MPI_Status time;    chkMPIerr(MPI_Wait(&(sendBuf[ii].time), &time));    __NOTE__("sendBuf[%d].time\n", ii);
#else///BLOCK_TIME_STEP
      MPI_Status   vx;    chkMPIerr(MPI_Wait(&(sendBuf[ii].  vx), & vx));    __NOTE__("sendBuf[%d].vx\n", ii);
      MPI_Status   vy;    chkMPIerr(MPI_Wait(&(sendBuf[ii].  vy), & vy));    __NOTE__("sendBuf[%d].vy\n", ii);
      MPI_Status   vz;    chkMPIerr(MPI_Wait(&(sendBuf[ii].  vz), & vz));    __NOTE__("sendBuf[%d].vz\n", ii);
#endif//BLOCK_TIME_STEP
#ifdef  GADGET_MAC
      MPI_Status  acc;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. acc), &acc));    __NOTE__("sendBuf[%d].acc\n", ii);
#endif//GADGET_MAC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      MPI_Status  ext;    chkMPIerr(MPI_Wait(&(sendBuf[ii].acc_ext), &ext));    __NOTE__("sendBuf[%d].ext\n", ii);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* if( sendBuf[ii].num > 0 ){ */
  }/* for(int ii = 0; ii < overlapNum; ii++){ */

  for(int ii = 0; ii < numProcs; ii++)
    if( recvBuf[ii].num > 0 ){
      __NOTE__("recv %d bodies from rank %d\n", recvBuf[ii].num, ii);
      MPI_Status  pos;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. pos), &pos));      __NOTE__("recvBuf[%d].pos\n", ii);
      MPI_Status  idx;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. idx), &idx));      __NOTE__("recvBuf[%d].idx\n", ii);
#ifdef  BLOCK_TIME_STEP
      MPI_Status  vel;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. vel), &vel));      __NOTE__("recvBuf[%d].vel\n", ii);
      MPI_Status time;      chkMPIerr(MPI_Wait(&(recvBuf[ii].time), &time));      __NOTE__("recvBuf[%d].time\n", ii);
#else///BLOCK_TIME_STEP
      MPI_Status   vx;      chkMPIerr(MPI_Wait(&(recvBuf[ii].  vx), & vx));      __NOTE__("recvBuf[%d].vx\n", ii);
      MPI_Status   vy;      chkMPIerr(MPI_Wait(&(recvBuf[ii].  vy), & vy));      __NOTE__("recvBuf[%d].vy\n", ii);
      MPI_Status   vz;      chkMPIerr(MPI_Wait(&(recvBuf[ii].  vz), & vz));      __NOTE__("recvBuf[%d].vz\n", ii);
#endif//BLOCK_TIME_STEP
#ifdef  GADGET_MAC
      MPI_Status  acc;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. acc), &acc));      __NOTE__("recvBuf[%d].acc\n", ii);
#endif//GADGET_MAC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      MPI_Status  ext;      chkMPIerr(MPI_Wait(&(recvBuf[ii].acc_ext), &ext));      __NOTE__("recvBuf[%d].ext\n", ii);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* if( recvBuf[ii].num > 0 ){ */

#endif//MPI_ONE_SIDED_FOR_EXCG


  for(int ii = 0; ii < numProcs; ii++){
    MPI_Status dom;
    chkMPIerr(MPI_Wait(&(domain.req[ii]), &dom));
  }/* for(int ii = 0; ii < numProcs; ii++){ */


  /** confirmation */
  if( *numNew > numMax ){
    __KILL__(stderr, "ERROR: # of required receive buffer (%d) exceeds the maximum number of particles per process (%d).\n\tsuggestion: consider increasing \"MAX_FACTOR_INCREASE\" or \"MAX_FACTOR_SAFETY\" defined in src/para/mpicfg.h (current values are %f and %f, respectively) at least %e times.\n", *numNew, numMax, MAX_FACTOR_INCREASE, MAX_FACTOR_SAFETY, (float)(*numNew) / (float)numMax);
  }/* if( *numNew > numMax ){ */

  const int diff = (*numNew) - numOld;
  int diff_sum;
  chkMPIerr(MPI_Reduce(&diff, &diff_sum, 1, MPI_INT, MPI_SUM, 0, mpi.comm));
  if( mpi.rank == 0 )
    if( diff_sum != 0 ){
      __KILL__(stderr, "ERROR: domain decomposition cause some error (duplication of %d particles)\n", diff_sum);
    }/* if( diff_sum != 0 ){ */


#ifdef  MPI_VIA_HOST
  /** copy N-body particles from host to device */
  copyParticle_hst2dev(*numNew, *dst_hst, *src_dev
#ifdef  EXEC_BENCHMARK
		       , elapsed
#endif//EXEC_BENCHMARK
		       );

#ifndef NDEBUG
  for(int ii = 0; ii < *numNew; ii++)
    if( (*dst_hst).pos[ii].m < EPSILON )
      __NOTE__("pi[%d]: recvmass = %e, x = %e, y = %e, z = %e\n", ii, (*dst_hst).pos[ii].m, (*dst_hst).pos[ii].x, (*dst_hst).pos[ii].y, (*dst_hst).pos[ii].z);
#endif//NDEBUG

#endif//MPI_VIA_HOST
  __NOTE__("numOld = %d, numNew = %d\n", numOld, *numNew);

  iparticle _tmp_hst;
  _tmp_hst = *src_hst;
  *src_hst = *dst_hst;
  *dst_hst = _tmp_hst;
#ifndef MPI_VIA_HOST
  _tmp_hst = *src_dev;
  *src_dev = *dst_dev;
  *dst_dev = _tmp_hst;
#endif//MPI_VIA_HOST


  /** modify parameters related to auto-tuning */
  const double scale = (double)(*numNew) / (double)numOld;
  measured->walkTree[0] *= scale;
  measured->walkTree[1] *= scale;
  measured->makeTree    *= scale;
  status->x.val *= scale;
  status->w.val *= scale;
  status->v.val *= scale;
  status->u.val *= scale;
  memory->previous *= scale;
#if 0
  if( memory->interval >= BRENT_METHOD_LAUNCH )
    memory->interval = BRENT_METHOD_LAUNCH - 1;
#endif

  /** reset counter */
  measured->sum_excg = 0.0;


  static struct timespec finish;
  checkCudaErrors(cudaDeviceSynchronize());
  clock_gettime(CLOCK_MONOTONIC_RAW, &finish);
  measured->excg = calcElapsedTimeInSec(start, finish);
#ifdef  EXEC_BENCHMARK
  elapsed->excgBody_dev = measured->excg;
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}
