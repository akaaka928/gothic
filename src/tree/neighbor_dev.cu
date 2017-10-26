/**
 * @file neighbor_dev.cu
 *
 * @brief Source code for neighbor search using breadth-first octree on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/10/26 (Thu)
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
#include <helper_cuda.h>

#include "macro.h"
#include "cudalib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "../util/gsync_dev.cu"

#include "make.h"
#include "make_dev.h"/**< to read NBLOCKS_PER_SM_MAC */
#include "neighbor_dev.h"


/**
 * @fn facileNeighborSearching_kernel
 *
 * @brief Execute simplified neighbor search.
 */
__global__ void facileNeighborSearching_kernel(const int Ni, READ_ONLY position * RESTRICT ibody, real * RESTRICT neighbor_length)
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;

  __shared__ position ipos[NTHREADS_FACILE_NS + NEIGHBOR_NUM];


  if( gidx < Ni ){
    /** load position of i-particles and neighbor candidates */
    ipos[tidx] = ibody[gidx];
    if( tidx < NEIGHBOR_NUM ){
      int idx = gidx + NTHREADS_FACILE_NS;      if( idx > (Ni - 1) )      	idx = Ni - 1;
      ipos[NTHREADS_FACILE_NS + tidx] = ibody[idx];
    }/* if( tidx < NEIGHBOR_NUM ){ */
    __syncthreads();

    /** calculate distance with NEIGHBOR_NUM particles and remember the maximum */
    const position pi = ipos[tidx];
    real r2max = ZERO;
#pragma unroll
    for(int ii = 0; ii < NEIGHBOR_NUM; ii++){
      const position pj = ipos[tidx + 1 + ii];

      const real dx = pj.x - pi.x;
      const real dy = pj.y - pi.y;
      const real dz = pj.z - pi.z;

      r2max = FMAX(r2max, FLT_MIN + dx * dx + dy * dy + dz * dz);
    }/* for(int ii = 0; ii < NEIGHBOR_NUM; ii++){ */

    /** store the derived guess about the length of neighbor arm */
    neighbor_length[gidx] = r2max * RSQRT(r2max);
  }/* if( gidx < Ni ){ */
}


/**
 * @fn facileNeighborSearching_dev
 *
 * @brief Execute simplified neighbor search.
 *
 * @sa facileNeighborSearching_kernel
 */
extern "C"
void facileNeighborSearching_dev(const int Ni, const iparticle pi)
{
  __NOTE__("%s\n", "start");


  int Nrem = BLOCKSIZE(Ni, NTHREADS_FACILE_NS);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    facileNeighborSearching_kernel<<<Nrem, NTHREADS_FACILE_NS>>>(Ni, pi.pos, pi.neighbor);
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_FACILE_NS;
      facileNeighborSearching_kernel<<<Nblck, NTHREADS_FACILE_NS>>>(Nsub, &pi.pos[hidx], &pi.neighbor[hidx]);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("facileNeighborSearching_kernel");


  __NOTE__("%s\n", "end");
}


#ifdef  SMEM_PREF_FOR_NEIGHBOR_SEARCH
/**
 * @fn setGlobalConstants_neighbor_dev_cu
 *
 * @brief Initialize kernel function in neighbor_dev.cu.
 */
extern "C"
void setGlobalConstants_neighbor_dev_cu(void)
{
  __NOTE__("%s\n", "start");

  checkCudaErrors(cudaFuncSetCacheConfig(facileNeighborSearching_kernel, cudaFuncCachePreferShared));

  __NOTE__("%s\n", "end");
}
#endif//SMEM_PREF_FOR_NEIGHBOR_SEARCH
