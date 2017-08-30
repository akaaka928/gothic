/**
 * @file exchange_dev.cu
 *
 * @brief Source code for domain decomposition using GPUs with MPI
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/08/30 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
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



#define NTHREADS_COMPARE_VEC3_INC NTHREADS_BOX
#include "../util/compare_vec3_inc.cu"
/**
 * @fn getBoxSize_kernel
 *
 * @brief Calculate box size to include all N-body particles in the local domain.
 */
__global__ void __launch_bounds__(NTHREADS_BOX, NBLOCKS_PER_SM_BOX) getBoxSize_kernel
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
    if( ii < itail ){
      const position pi = ipos[ii];
      const float xi = CAST_R2F(pi.x);
      const float yi = CAST_R2F(pi.y);
      const float zi = CAST_R2F(pi.z);

      min.x = fminf(min.x, xi);      max.x = fmaxf(max.x, xi);
      min.y = fminf(min.y, yi);      max.y = fmaxf(max.y, yi);
      min.z = fminf(min.z, zi);      max.z = fmaxf(max.z, zi);
    }/* if( ii < itail ){ */
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
void checkBoxSize_dev(const deviceProp devProp)
{
  __NOTE__("%s\n", "start");


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


  getBoxSize_kernel<<<devProp.numSM * NBLOCKS_PER_SM_BOX, NTHREADS_BOX, SMEM_SIZE, stream>>>(num, ipos, soa.min, soa.max, soa.gsync0, soa.gsync1);
  getLastCudaError("getBoxSize_kernel");

  checkCudaErrors(cudaMemcpyAsync(soa.min_hst, soa.min, sizeof(float4), cudaMemcpyDeviceToHost, stream));
  checkCudaErrors(cudaMemcpyAsync(soa.max_hst, soa.max, sizeof(float4), cudaMemcpyDeviceToHost, stream));

  checkCudaErrors(cudaDeviceSynchronize());


  __NOTE__("%s\n", "end");
}


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


  checkCudaErrors(cudaMemcpy(src.x_dev, src.x_hst, num * sizeof(float), cudaMemcpyHostToDevice));
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
  checkCudaErrors(cudaMemcpy(dst.x_hst, src.x_dev, num * sizeof(float), cudaMemcpyDeviceToHost));

  checkCudaErrors(cudaMemcpy(src.y_dev, src.y_hst, num * sizeof(float), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(src.z_dev, src.z_hst, num * sizeof(float), cudaMemcpyHostToDevice));
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
  checkCudaErrors(cudaMemcpy(dst.y_hst, dst.y_dev, num * sizeof(float), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(dst.z_hst, dst.z_dev, num * sizeof(float), cudaMemcpyDeviceToHost));


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

  __NOTE__("%s\n", "end");
}

/**
 * @fn sort_xpos
 *
 * @brief Sort position array by x.
 */
static inline void sort_xpos(const int num, samplePos * RESTRICT src, samplePos * RESTRICT dst)
{
  __NOTE__("%s\n", "start");

  if( num >= NCRIT_XPOS_SORT )    sort_xpos_dev(num, *src, *dst);
  else                            sort_xpos_hst(num, *src, *dst);

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


  checkCudaErrors(cudaMemcpy(src.y_dev, src.y_hst, num * sizeof(float), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(src.z_dev, src.z_hst, num * sizeof(float), cudaMemcpyHostToDevice));

  /** sort using thrust */
  thrust::stable_sort_by_key((thrust::device_ptr<float>)(src.y_dev), (thrust::device_ptr<float>)((src.y_dev) + num), (thrust::device_ptr<float>)(src.z_dev));

  checkCudaErrors(cudaMemcpy(src.y_hst, src.y_dev, num * sizeof(float), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(src.z_hst, src.z_dev, num * sizeof(float), cudaMemcpyDeviceToHost));


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

  __NOTE__("%s\n", "end");
}

/**
 * @fn sort_ypos
 *
 * @brief Sort position array by y.
 */
static inline void sort_ypos(const int num, samplePos * RESTRICT src)
{
  __NOTE__("%s\n", "start");

  if( num >= NCRIT_YPOS_SORT )    sort_ypos_dev(num, *src);
  else                            sort_ypos_hst(num, *src);

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

  checkCudaErrors(cudaMemcpy(src.z_dev, src.z_hst, num * sizeof(float), cudaMemcpyHostToDevice));

  /** sort using thrust */
  thrust::stable_sort((thrust::device_ptr<float>)(src.z_dev), (thrust::device_ptr<float>)((src.z_dev) + num));

  checkCudaErrors(cudaMemcpy(src.z_hst, src.z_dev, num * sizeof(float), cudaMemcpyDeviceToHost));

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

  __NOTE__("%s\n", "end");
}

/**
 * @fn sort_zpos
 *
 * @brief Sort position array by z.
 */
static inline void sort_zpos(const int num, samplePos * RESTRICT src)
{
  __NOTE__("%s\n", "start");

  if( num >= NCRIT_ZPOS_SORT )    sort_zpos_dev(num, *src);
  else                            sort_zpos_hst(num, *src);

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
#ifdef  BLOCK_TIME_STEP
 velocity     * RESTRICT dvel , READ_ONLY velocity     * RESTRICT svel ,
 ibody_time   * RESTRICT dtime, READ_ONLY ibody_time   * RESTRICT stime
#else///BLOCK_TIME_STEP
 real         * RESTRICT dvx  , READ_ONLY real         * RESTRICT svx  ,
 real         * RESTRICT dvy  , READ_ONLY real         * RESTRICT svy  ,
 real         * RESTRICT dvz  , READ_ONLY real         * RESTRICT svz
#endif//BLOCK_TIME_STEP
 , const int offset
)
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
#ifdef  BLOCK_TIME_STEP
    dvel [ii] = svel [jj];
    dtime[ii] = stime[jj];
#else///BLOCK_TIME_STEP
    dvx  [ii] = svx  [jj];
    dvy  [ii] = svy  [jj];
    dvz  [ii] = svz  [jj];
#endif//BLOCK_TIME_STEP
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
			      int **rank_hst, int **rank_dev, int **idx_dev, domainDecomposeKey *key, const ulong Ntot)
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  size_t num = (size_t)((float)Ntot * MAX_FACTOR_FROM_EQUIPARTITION);
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
 particlePos pos_hst, particlePos pos_dev, domainDecomposeKey key, sendCfg *sendBuf, recvCfg *recvBuf,
 MPIinfo orm[], MPIinfo rep[], domainCfg domain, MPIcfg_tree mpi,
 const double tloc, sampling sample, samplePos loc, samplePos ful,
 soaPHsort soa, const deviceProp devProp, const deviceInfo devInfo,
 measuredTime *measured
#ifdef  CARE_EXTERNAL_PARTICLES
 , domainLocation *location
#endif//CARE_EXTERNAL_PARTICLES
 , brentStatus *status, brentMemory *memory
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
 )
{
  __NOTE__("%s\n", "start");


  static struct timespec start;
  checkCudaErrors(cudaDeviceSynchronize());
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);

  /** set current (local) distribution */
  /** copy particle position from device to host */
  copyParticlePositionAsync_dev2hst(numOld, (*src_dev).pos, pos_dev, pos_hst, devInfo.stream[0]);

  /** get (current) box size for the local distribution of N-body particles */
  getBoxSize_dev(numOld, (*src_dev).pos, soa, devProp, devInfo.stream[1]);
  float3 min, max;
  min.x = soa.min_hst->x;
  min.y = soa.min_hst->y;
  min.z = soa.min_hst->z;
  max.x = soa.max_hst->x;
  max.y = soa.max_hst->y;
  max.z = soa.max_hst->z;

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


  /** pick up sample particles */
  /** weight is determined using elapsed time by each process */
  __NOTE__("rank %d: tloc = %e, numOld = %d, xmin = %e, xmax = %e, ymin = %e, ymax = %e, zmin = %e, zmax = %e\n", mpi.rank, tloc, numOld, min.x, max.x, min.y, max.y, min.z, max.z);
  double ttot = tloc;
  chkMPIerr(MPI_Allreduce(&tloc, &ttot, 1, MPI_DOUBLE, MPI_SUM, mpi.comm));
  const float frac = fminf((float)(tloc / ttot), MAX_FACTOR_INCREASE / (float)mpi.size);
  const int Nsub = (int)ceilf((float)Ntot * sample.rate * frac);
  __NOTE__("rank %d: tloc = %e, ttot = %e, frac = %e, numOld = %d, Nsub = %d\n", mpi.rank, tloc, ttot, tloc / ttot, numOld, Nsub);

  const int iskip = (Nsub < numOld) ? (numOld / Nsub) : (1);
  int sendNum = 0;
  for(int ii = 0; ii < numOld; ii += iskip){
    loc.x_hst[sendNum] = pos_hst.x[ii];
    loc.y_hst[sendNum] = pos_hst.y[ii];
    loc.z_hst[sendNum] = pos_hst.z[ii];
    sendNum++;
  }/* for(int ii = 0; ii < numOld; ii += iskip){ */
  __NOTE__("rank %d: iskip = %d, sendNum = %d\n", mpi.rank, iskip, sendNum);

  /** sort sample particles in each direction */
  if(             mpi.dim[0] != 1 )     sort_xpos(sendNum, &loc, &ful);
  else{	   if(	  mpi.dim[1] != 1 )	sort_ypos(sendNum, &loc);
    else      if( mpi.dim[2] != 1 )	sort_zpos(sendNum, &loc);
  }
  __NOTE__("rank %d\n", mpi.rank);


  /** gather sampling points to the root process */
  chkMPIerr(MPI_Gather(&sendNum, 1, MPI_INT, sample.rnum, 1, MPI_INT, 0, mpi.comm));

  /** set receive displacements (root process only) */
  int recvNum = 0;
  if( mpi.rank == 0 ){
    sample.disp[0] = 0;
    for(int jj = 1; jj < mpi.size; jj++)
      sample.disp[jj] = sample.disp[jj - 1] + sample.rnum[jj - 1];
    recvNum = sample.disp[mpi.size - 1] + sample.rnum[mpi.size - 1];
  }/* if( orm[ii].rank == 0 ){ */

  /** gather particle data to the root process */
  chkMPIerr(MPI_Gatherv(loc.x_hst, sendNum, MPI_FLOAT, ful.x_hst, sample.rnum, sample.disp, MPI_REALDAT, 0, mpi.comm));
  chkMPIerr(MPI_Gatherv(loc.y_hst, sendNum, MPI_FLOAT, ful.y_hst, sample.rnum, sample.disp, MPI_REALDAT, 0, mpi.comm));
  chkMPIerr(MPI_Gatherv(loc.z_hst, sendNum, MPI_FLOAT, ful.z_hst, sample.rnum, sample.disp, MPI_REALDAT, 0, mpi.comm));


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
	sort_xpos(recvNum, &ful, &loc);
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
	chkMPIerr(MPI_Scatterv(ful.y_hst, sample.rnum, sample.disp, MPI_FLOAT, MPI_IN_PLACE, recvNum, MPI_FLOAT, 0, rep[0].comm));
	chkMPIerr(MPI_Scatterv(ful.z_hst, sample.rnum, sample.disp, MPI_FLOAT, MPI_IN_PLACE, recvNum, MPI_FLOAT, 0, rep[0].comm));
      }/* if( (mpi.dim[1] != 1) || (mpi.dim[2] != 1) ){ */
    }/* if( orm[0].rank == 0 ){ */

    /** MPI_Bcast in orm[0].comm */
    chkMPIerr(MPI_Bcast(&local_xmin, 1, MPI_FLOAT, 0, orm[0].comm));
    chkMPIerr(MPI_Bcast(&local_xmax, 1, MPI_FLOAT, 0, orm[0].comm));
  }/* if( mpi.dim[0] != 1 ){ */


  /** domain decomposition in y-direction */
  if( mpi.dim[1] != 1 ){
    if( orm[1].rank == 0 ){
      __NOTE__("y-decomposition by mpi.rank = %d\n", mpi.rank);
      sendNum = recvNum;
      /** the root process determine the partition */
      if( rep[1].rank == 0 ){
	sort_ypos(recvNum, &ful);
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
	chkMPIerr(MPI_Scatterv(ful.z_hst, sample.rnum, sample.disp, MPI_FLOAT, MPI_IN_PLACE, recvNum, MPI_FLOAT, 0, rep[1].comm));
      }/* if( mpi.dim[2] != 1 ){ */
    }/* if( orm[1].rank == 0 ){ */

    /** MPI_Bcast in orm[1].comm */
    chkMPIerr(MPI_Bcast(&local_ymin, 1, MPI_FLOAT, 0, orm[1].comm));
    chkMPIerr(MPI_Bcast(&local_ymax, 1, MPI_FLOAT, 0, orm[1].comm));
  }/* if( mpi.dim[1] != 1 ){ */


  /** domain decomposition in z-direction */
  if( mpi.dim[2] != 1 ){
    if( orm[2].rank == 0 ){
      __NOTE__("z-decomposition by mpi.rank = %d\n", mpi.rank);
      sendNum = recvNum;
      /** the root process determine the partition */
      if( rep[2].rank == 0 ){
	sort_zpos(recvNum, &ful);
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
  __NOTE__("rank %d: [%e, %e]x[%e, %e]x[%e, %e]\n", mpi.rank,
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
  for(int ii = 0; ii < numProcs; ii++)
    chkMPIerr(MPI_Irecv(&(recvBuf[ii].num), 1, MPI_INT, ii, ii, mpi.comm, &(recvBuf[ii].req)));

  const int nullSend = 0;
  int overlapNum = 0;
  for(int ii = 0; ii < numProcs; ii++){
    if( (min.x <= domain.xmax[ii]) && (max.x >= domain.xmin[ii]) &&
	(min.y <= domain.ymax[ii]) && (max.y >= domain.ymin[ii]) &&
	(min.z <= domain.zmax[ii]) && (max.z >= domain.zmin[ii]) ){
      /** if spatial overlap is detected, ... */
      sendBuf[overlapNum].rank = ii;
      sendBuf[overlapNum].num  =  0;

      sendBuf[overlapNum].xmin = (domain.xmin[ii] < -0.25f * FLT_MAX) ? (domain.xmin[ii]) : ((min.x > domain.xmin[ii]) ? (min.x) : (domain.xmin[ii]));
      sendBuf[overlapNum].ymin = (domain.ymin[ii] < -0.25f * FLT_MAX) ? (domain.ymin[ii]) : ((min.y > domain.ymin[ii]) ? (min.y) : (domain.ymin[ii]));
      sendBuf[overlapNum].zmin = (domain.zmin[ii] < -0.25f * FLT_MAX) ? (domain.zmin[ii]) : ((min.z > domain.zmin[ii]) ? (min.z) : (domain.zmin[ii]));

      sendBuf[overlapNum].xmax = (domain.xmax[ii] >  0.25f * FLT_MAX) ? (domain.xmax[ii]) : ((max.x < domain.xmax[ii]) ? (max.x) : (domain.xmax[ii]));
      sendBuf[overlapNum].ymax = (domain.ymax[ii] >  0.25f * FLT_MAX) ? (domain.ymax[ii]) : ((max.y < domain.ymax[ii]) ? (max.y) : (domain.ymax[ii]));
      sendBuf[overlapNum].zmax = (domain.zmax[ii] >  0.25f * FLT_MAX) ? (domain.zmax[ii]) : ((max.z < domain.zmax[ii]) ? (max.z) : (domain.zmax[ii]));

      overlapNum++;
    }
    else{
      /* if covered areas do not overlap, ... */
      chkMPIerr(MPI_Isend(&nullSend, 1, MPI_INT, ii, mpi.rank, mpi.comm, &(domain.req[ii])));
    }/* else{ */

    __NOTE__("rank %d, dst = %d\n", mpi.rank, ii);
  }/* for(int ii = 0; ii < numProcs; ii++){ */


  /** determine process rank for each particle to belong */
  for(int ii = 0; ii < numOld; ii++){
#ifndef NDEBUG
    bool find = false;
#endif//NDEBUG

    for(int jj = 0; jj < overlapNum; jj++){
      if( (pos_hst.x[ii] >= sendBuf[jj].xmin) && (pos_hst.x[ii] <= sendBuf[jj].xmax) &&
	  (pos_hst.y[ii] >= sendBuf[jj].ymin) && (pos_hst.y[ii] <= sendBuf[jj].ymax) &&
	  (pos_hst.z[ii] >= sendBuf[jj].zmin) && (pos_hst.z[ii] <= sendBuf[jj].zmax) ){
	key.dstRank_hst[ii] = sendBuf[jj].rank;
	sendBuf[jj].num++;
#ifndef NDEBUG
	find = true;
#endif//NDEBUG
	break;
      }
    }/* for(int jj = 0; jj < overlapNum; jj++){ */

#ifndef NDEBUG
    if( !find ){
#if 0
      fprintf(stderr, "numOld = %d, numProcs = %d, overlapNum = %d @ rank %d\n", numOld, numProcs, overlapNum, mpi.rank);
      fprintf(stderr, "local: xmin = %f, xmax = %f, ymin = %f, ymax = %f, zmin = %f, zmax = %f @ rank %d\n",
	      min.x, max.x, min.y, max.y, min.z, max.z, mpi.rank);
      for(int jj = 0; jj < numProcs; jj++){
	fprintf(stderr, "domain[%d]: xmin = %f, xmax = %f, ymin = %f, ymax = %f, zmin = %f, zmax = %f @ rank %d\n",
		jj, domain.xmin[jj], domain.xmax[jj], domain.ymin[jj], domain.ymax[jj], domain.zmin[jj], domain.zmax[jj], mpi.rank);
      }
      for(int jj = 0; jj < overlapNum; jj++){
	fprintf(stderr, "sendBuf[%d]: xmin = %f, xmax = %f, ymin = %f, ymax = %f, zmin = %f, zmax = %f @ rank %d\n",
		jj, sendBuf[jj].xmin, sendBuf[jj].xmax, sendBuf[jj].ymin, sendBuf[jj].ymax, sendBuf[jj].zmin, sendBuf[jj].zmax, mpi.rank);
      }
      fprintf(stderr, "ii = %d: x = %f, y = %f, z = %f @ rank %d\n", ii, pos_hst.x[ii], pos_hst.y[ii], pos_hst.z[ii], mpi.rank);
#endif
      __KILL__(stderr, "ERROR: target MPI rank is missing\n");
    }
#endif//NDEBUG
  }/* for(int ii = 0; ii < numOld; ii++){ */

  sendBuf[0].head = 0;
  for(int ii = 0; ii < overlapNum; ii++){
    chkMPIerr(MPI_Isend(&(sendBuf[ii].num), 1, MPI_INT, sendBuf[ii].rank, mpi.rank, mpi.comm, &(domain.req[sendBuf[ii].rank])));
    if( ii > 0 )
      sendBuf[ii].head = sendBuf[ii - 1].head + sendBuf[ii - 1].num;
  }/* for(int ii = 0; ii < overlapNum; ii++){ */

  if( (sendBuf[overlapNum - 1].head + sendBuf[overlapNum - 1].num) != numOld ){
    __KILL__(stderr, "ERROR: total number of scattered particles (%d) is differ from that of local particles (%d)\n",
	     sendBuf[overlapNum - 1].head + sendBuf[overlapNum - 1].num, numOld);
  }

  sortDomainDecomposeKey(numOld, key, src_dev, dst_dev);

  for(int ii = 0; ii < numProcs; ii++){
    MPI_Status status;
    chkMPIerr(MPI_Wait(&(domain.req[ii]), &status));
  }/* for(int ii = 0; ii < numProcs; ii++){ */


  /** if GPU-GPU direct communication is not available */
  copyParticle_dev2hst(numOld, *src_dev, *src_hst
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );
  for(int ii = 0; ii < overlapNum; ii++){
    chkMPIerr(MPI_Isend(&((*src_hst).pos [sendBuf[ii].head]), sendBuf[ii].num, mpi.ipos         , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].pos)));
#ifdef  GADGET_MAC
    chkMPIerr(MPI_Isend(&((*src_hst).acc [sendBuf[ii].head]), sendBuf[ii].num, mpi.iacc         , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].acc)));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
    chkMPIerr(MPI_Isend(&((*src_hst).vel [sendBuf[ii].head]), sendBuf[ii].num, mpi.ivel         , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].vel)));
    chkMPIerr(MPI_Isend(&((*src_hst).time[sendBuf[ii].head]), sendBuf[ii].num, mpi.time         , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].time)));
#else///BLOCK_TIME_STEP
    chkMPIerr(MPI_Isend(&((*src_hst).vx  [sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT      , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].vx )));
    chkMPIerr(MPI_Isend(&((*src_hst).vy  [sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT      , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].vy )));
    chkMPIerr(MPI_Isend(&((*src_hst).vz  [sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT      , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].vz )));
#endif//BLOCK_TIME_STEP
    chkMPIerr(MPI_Isend(&((*src_hst).idx [sendBuf[ii].head]), sendBuf[ii].num, MPI_UNSIGNED_LONG, sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].idx)));
  }/* for(int ii = 0; ii < overlapNum; ii++){ */


  /** src is receive buffer */
  *numNew = 0;
  for(int ii = 0; ii < numProcs; ii++)
    recvBuf[ii].head = 0;
  for(int ii = 0; ii < numProcs; ii++){
    /** receive recvNum */
    MPI_Status status;
    chkMPIerr(MPI_Wait(&(recvBuf[ii].req), &status));

    /** if recvNum != 0, then set receive buffer */
    if( recvBuf[ii].num != 0 ){
      chkMPIerr(MPI_Irecv(&((*dst_hst).pos [recvBuf[ii].head]), recvBuf[ii].num, mpi.ipos         , ii, ii, mpi.comm, &(recvBuf[ii].pos)));
#ifdef  GADGET_MAC
      chkMPIerr(MPI_Irecv(&((*dst_hst).acc [recvBuf[ii].head]), recvBuf[ii].num, mpi.iacc         , ii, ii, mpi.comm, &(recvBuf[ii].acc)));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
      chkMPIerr(MPI_Irecv(&((*dst_hst).vel [recvBuf[ii].head]), recvBuf[ii].num, mpi.ipos         , ii, ii, mpi.comm, &(recvBuf[ii].vel)));
      chkMPIerr(MPI_Irecv(&((*dst_hst).time[recvBuf[ii].head]), recvBuf[ii].num, mpi.ipos         , ii, ii, mpi.comm, &(recvBuf[ii].time)));
#else///BLOCK_TIME_STEP
      chkMPIerr(MPI_Irecv(&((*dst_hst).vx  [recvBuf[ii].head]), recvBuf[ii].num, MPI_REALDAT      , ii, ii, mpi.comm, &(recvBuf[ii].vx )));
      chkMPIerr(MPI_Irecv(&((*dst_hst).vy  [recvBuf[ii].head]), recvBuf[ii].num, MPI_REALDAT      , ii, ii, mpi.comm, &(recvBuf[ii].vy )));
      chkMPIerr(MPI_Irecv(&((*dst_hst).vz  [recvBuf[ii].head]), recvBuf[ii].num, MPI_REALDAT      , ii, ii, mpi.comm, &(recvBuf[ii].vz )));
#endif//BLOCK_TIME_STEP
      chkMPIerr(MPI_Irecv(&((*dst_hst).idx [recvBuf[ii].head]), recvBuf[ii].num, MPI_UNSIGNED_LONG, ii, ii, mpi.comm, &(recvBuf[ii].idx)));

      *numNew += recvBuf[ii].num;

      if( ii + 1 < numProcs )
	recvBuf[ii + 1].head = recvBuf[ii].head + recvBuf[ii].num;
    }/* if( recvBuf[ii].num != 0 ){ */
  }/* for(int ii = 0; ii < numProcs; ii++){ */

  if( *numNew > numMax ){
    __KILL__(stderr, "ERROR: # of required receive buffer (%d) exceeds the maximum number of particles per process (%d).\n\tsuggestion: consider increasing \"MAX_FACTOR_INCREASE\" or \"MAX_FACTOR_SAFETY\" defined in src/para/mpicfg.h (current values are %f and %f, respectively) at least %e times.\n", *numNew, numMax, MAX_FACTOR_INCREASE, MAX_FACTOR_SAFETY, (float)(*numNew) / (float)numMax);
  }/* if( *numNew > numMax ){ */


  /** complete MPI communications */
  for(int ii = 0; ii < overlapNum; ii++){
    MPI_Status  pos;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. pos), &pos));
#ifdef  GADGET_MAC
    MPI_Status  acc;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. acc), &acc));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
    MPI_Status  vel;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. vel), &vel));
    MPI_Status time;    chkMPIerr(MPI_Wait(&(sendBuf[ii].time), &time));
#else///BLOCK_TIME_STEP
    MPI_Status   vx;    chkMPIerr(MPI_Wait(&(sendBuf[ii].  vx), & vx));
    MPI_Status   vy;    chkMPIerr(MPI_Wait(&(sendBuf[ii].  vy), & vy));
    MPI_Status   vz;    chkMPIerr(MPI_Wait(&(sendBuf[ii].  vz), & vz));
#endif//BLOCK_TIME_STEP
    MPI_Status  idx;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. idx), &idx));
  }/* for(int ii = 0; ii < overlapNum; ii++){ */

  for(int ii = 0; ii < numProcs; ii++)
    if( recvBuf[ii].num != 0 ){
      MPI_Status  pos;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. pos), &pos));
#ifdef  GADGET_MAC
      MPI_Status  acc;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. acc), &acc));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
      MPI_Status  vel;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. vel), &vel));
      MPI_Status time;      chkMPIerr(MPI_Wait(&(recvBuf[ii].time), &time));
#else///BLOCK_TIME_STEP
      MPI_Status   vx;      chkMPIerr(MPI_Wait(&(recvBuf[ii].  vx), & vx));
      MPI_Status   vy;      chkMPIerr(MPI_Wait(&(recvBuf[ii].  vy), & vy));
      MPI_Status   vz;      chkMPIerr(MPI_Wait(&(recvBuf[ii].  vz), & vz));
#endif//BLOCK_TIME_STEP
      MPI_Status  idx;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. idx), &idx));
    }/* if( recvBuf[ii].num != 0 ){ */

  /** confirmation */
  const int diff = (*numNew) - numOld;
  int diff_sum;
  chkMPIerr(MPI_Reduce(&diff, &diff_sum, 1, MPI_INT, MPI_SUM, 0, mpi.comm));
  if( mpi.rank == 0 )
    if( diff_sum != 0 ){
      __KILL__(stderr, "ERROR: domain decomposition cause some error (duplication of %d particles)\n", diff_sum);
    }/* if( diff_sum != 0 ){ */


  /** copy N-body particles from host to device */
  copyParticle_hst2dev(*numNew, *dst_hst, *src_dev
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );
  __NOTE__("numOld = %d, numNew = %d @ rank %d\n", numOld, *numNew, mpi.rank);

  iparticle _tmp_hst;
  _tmp_hst = *src_hst;
  *src_hst = *dst_hst;
  *dst_hst = _tmp_hst;


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

  /** reset counter */
  measured->sum_excg = 0.0;


  static struct timespec finish;
  checkCudaErrors(cudaDeviceSynchronize());
  clock_gettime(CLOCK_MONOTONIC_RAW, &finish);
  measured->excg = calcElapsedTimeInSec(start, finish);
#ifdef  EXEC_BENCHMARK
  elapsed->excgBody = measured->excg;
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}
