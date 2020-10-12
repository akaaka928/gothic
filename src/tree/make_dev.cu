/**
 * @file make_dev.cu
 *
 * @brief Source code for constructing octree structure on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/09/14 (Mon)
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

#   if  defined(USE_COOPERATIVE_GROUPS) || (GPUGEN >= 70)
#include <cooperative_groups.h>
using namespace cooperative_groups;
#endif//defined(USE_COOPERATIVE_GROUPS) || (GPUGEN >= 70)

#include "macro.h"
#include "cudalib.h"
#ifndef SERIALIZED_EXECUTION
#include <time.h>
#include "timer.h"
#endif//SERIALIZED_EXECUTION

#include "../misc/benchmark.h"
#include "../misc/structure.h"
#include "../misc/device.h"

#include "../util/gsync_dev.cu"

#include "../sort/peano.h"

#include "make.h"
#include "make_dev.h"
#include "let.h"/**< necessary to read EXTEND_NUM_TREE_NODE */

#include "../sort/peano_dev.h"


#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
#define MODIFY_INTERNAL_DOMAIN
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)


#ifndef USE_OCCUPANCY_CALCULATOR
/** limitation from capacity of shared memory */
/** capacity of shared memory is 64KiB on newer GPUs */
/** in shared memory preferred configuration, capacity of shared memory is 48KiB per SM on older GPUs */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#define NBLOCKS_PER_SM_MAC ( (SMEM_SIZE_SM_PREF >> 4) / (NTHREADS_MAC * (1 +     NBUF_MAC)))
/* #   if  GPUVER >= 60 */
/* /\* 4096 = 64 * 1024 / 16 *\/ */
/* #define NBLOCKS_PER_SM_MAC ( 4096 / (NTHREADS_MAC * (1 +     NBUF_MAC))) */
/* #else///GPUVER >= 60 */
/* /\* 3072 = 48 * 1024 / 16 *\/ */
/* #define NBLOCKS_PER_SM_MAC ( 3072 / (NTHREADS_MAC * (1 +     NBUF_MAC))) */
/* #endif//GPUVER >= 60 */
#else///USE_WARP_SHUFFLE_FUNC_MAC
#define NBLOCKS_PER_SM_MAC ((SMEM_SIZE_SM_PREF >> 2) / (NTHREADS_MAC * (5 + 4 * NBUF_MAC)))
/* #   if  GPUVER >= 60 */
/* /\* 16384 = 64 * 1024 / 4 *\/ */
/* #define NBLOCKS_PER_SM_MAC (16384 / (NTHREADS_MAC * (5 + 4 * NBUF_MAC))) */
/* #else///GPUVER >= 60 */
/* /\* 12288 = 48 * 1024 / 4 *\/ */
/* #define NBLOCKS_PER_SM_MAC (12288 / (NTHREADS_MAC * (5 + 4 * NBUF_MAC))) */
/* #endif//GPUVER >= 60 */
#endif//USE_WARP_SHUFFLE_FUNC_MAC

#define REGISTERS_PER_THREAD_MAC (80)
/* calcMultipole_kernel uses 60 registers @ Tesla M2090 */
#   if  GPUVER == 20
#undef  REGISTERS_PER_THREAD_MAC
#ifdef  GADGET_MAC
#define REGISTERS_PER_THREAD_MAC (59)
#else///GADGET_MAC
#define REGISTERS_PER_THREAD_MAC (60)
#endif//GADGET_MAC
#endif//GPUVER == 20
/* calcMultipole_kernel uses 57 registers @ Tesla K20X w/z warp shuffle */
/* calcMultipole_kernel uses 62 registers @ Tesla K20X w/o warp shuffle */
#   if  GPUVER == 35
#undef  REGISTERS_PER_THREAD_MAC
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#ifdef  GADGET_MAC
#define REGISTERS_PER_THREAD_MAC (55)
#else///GADGET_MAC
#ifdef  WS93_MAC
#define REGISTERS_PER_THREAD_MAC (57)
#else///WS93_MAC
#define REGISTERS_PER_THREAD_MAC (49)
#endif//WS93_MAC
#endif//GADGET_MAC
#else///USE_WARP_SHUFFLE_FUNC_MAC
#define REGISTERS_PER_THREAD_MAC (62)
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#endif//GPUVER == 35
/* calcMultipole_kernel uses 64 registers @ GTX 750 Ti w/z warp shuffle */
/* calcMultipole_kernel uses 72 registers @ GTX 750 Ti w/o warp shuffle */
#   if  GPUVER == 50
#undef  REGISTERS_PER_THREAD_MAC
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#define REGISTERS_PER_THREAD_MAC (64)
#else///USE_WARP_SHUFFLE_FUNC_MAC
#define REGISTERS_PER_THREAD_MAC (72)
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#endif//GPUVER == 50
/* calcMultipole_kernel uses 64 registers @ GTX 970 w/z warp shuffle */
/* calcMultipole_kernel uses 72 registers @ GTX 970 w/o warp shuffle */
#   if  GPUVER == 52
#undef  REGISTERS_PER_THREAD_MAC
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#define REGISTERS_PER_THREAD_MAC (64)
#else///USE_WARP_SHUFFLE_FUNC_MAC
#define REGISTERS_PER_THREAD_MAC (72)
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#endif//GPUVER == 52
/* calcMultipole_kernel uses 80 registers @ Tesla P100 */
#   if  GPUVER == 60
#undef  REGISTERS_PER_THREAD_MAC
#define REGISTERS_PER_THREAD_MAC (80)
#endif//GPUVER == 60

/** limitation from number of registers */
#   if  NBLOCKS_PER_SM_MAC > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAC * NTHREADS_MAC))
#undef  NBLOCKS_PER_SM_MAC
#define NBLOCKS_PER_SM_MAC   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAC * NTHREADS_MAC))
#endif//NBLOCKS_PER_SM_MAC > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAC * NTHREADS_MAC))

/** maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_MAC > (DIV_NTHREADS_MAC(MAX_THREADS_PER_SM))
#undef  NBLOCKS_PER_SM_MAC
#define NBLOCKS_PER_SM_MAC   (DIV_NTHREADS_MAC(MAX_THREADS_PER_SM))
#endif//NBLOCKS_PER_SM_MAC > (DIV_NTHREADS_MAC(MAX_THREADS_PER_SM))

/** maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_MAC > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_MAC
#define NBLOCKS_PER_SM_MAC   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_MAC > MAX_BLOCKS_PER_SM

/** maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_MAC > (DIV_NTHREADS_MAC(MAX_WARPS_PER_SM * 32))
#undef  NBLOCKS_PER_SM_MAC
#define NBLOCKS_PER_SM_MAC   (DIV_NTHREADS_MAC(MAX_WARPS_PER_SM * 32))
#endif//NBLOCKS_PER_SM_MAC > (DIV_NTHREADS_MAC(MAX_WARPS_PER_SM * 32))

/** # of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_MAC < 1
#undef  NBLOCKS_PER_SM_MAC
#define NBLOCKS_PER_SM_MAC  (1)
#endif//NBLOCKS_PER_SM_MAC < 1


/** limitations from capacity of shared memory (assuming L1 cache preferred) */
/** SM usage is 20 * NTHREADS_MAKE_TREE + 16 * NUM_PHKEY_LEVEL bytes */
#define NBLOCKS_PER_SM_MAKE_TREE   (SMEM_SIZE_L1_PREF / (20 * NTHREADS_MAKE_TREE + 16 * NUM_PHKEY_LEVEL))
/** SM usage is  4 * NTHREADS_LINK_TREE bytes */
#define NBLOCKS_PER_SM_LINK_TREE   (SMEM_SIZE_L1_PREF / (4 * NTHREADS_LINK_TREE))

#define REGISTERS_PER_THREAD_MAKE_TREE (47)
#define REGISTERS_PER_THREAD_LINK_TREE (25)
/* makeTree_kernel and linkTree_kernel use 53 and 27 registers, respectively @ Tesla M2090 */
#   if  GPUVER == 20
#undef  REGISTERS_PER_THREAD_MAKE_TREE
#define REGISTERS_PER_THREAD_MAKE_TREE (53)
#undef  REGISTERS_PER_THREAD_LINK_TREE
#define REGISTERS_PER_THREAD_LINK_TREE (27)
#endif//GPUVER == 20
/* makeTree_kernel and linkTree_kernel use 49 and 27 registers, respectively @ Tesla K20X */
#   if  GPUVER == 35
#undef  REGISTERS_PER_THREAD_MAKE_TREE
#define REGISTERS_PER_THREAD_MAKE_TREE (49)
#undef  REGISTERS_PER_THREAD_LINK_TREE
#define REGISTERS_PER_THREAD_LINK_TREE (27)
#endif//GPUVER == 35
/* makeTree_kernel and linkTree_kernel use 64 and 23 registers, respectively @ GTX 750 Ti w/z Ttot = 128, 512, 1024 */
/* makeTree_kernel and linkTree_kernel use 66 and 23 registers, respectively @ GTX 750 Ti w/z Ttot = 256 */
#   if  GPUVER == 50
#undef  REGISTERS_PER_THREAD_MAKE_TREE
#   if  NTHREADS_MAKE_TREE == 256
#define REGISTERS_PER_THREAD_MAKE_TREE (66)
#else///NTHREADS_MAKE_TREE == 256
#define REGISTERS_PER_THREAD_MAKE_TREE (64)
#endif//NTHREADS_MAKE_TREE == 256
#undef  REGISTERS_PER_THREAD_LINK_TREE
#define REGISTERS_PER_THREAD_LINK_TREE (23)
#endif//GPUVER == 50
/* makeTree_kernel and linkTree_kernel use 64 and 23 registers, respectively @ GTX 970 w/z Ttot = 128, 512, 1024 @ CUDA 7.5 */
/* makeTree_kernel and linkTree_kernel use 66 and 23 registers, respectively @ GTX 970 w/z Ttot = 256            @ CUDA 7.5 */
/* makeTree_kernel use 70 registers @ GTX TITAN X w/z Ttot = 128 @ CUDA 8.0 */
/* linkTree_kernel use 24 registers @ GTX TITAN X w/z Ttot = 256 @ CUDA 8.0 */
#   if  GPUVER == 52
#undef  REGISTERS_PER_THREAD_MAKE_TREE
#   if  NTHREADS_MAKE_TREE == 256
#define REGISTERS_PER_THREAD_MAKE_TREE (66)
#else///NTHREADS_MAKE_TREE == 256
#define REGISTERS_PER_THREAD_MAKE_TREE (62)
#endif//NTHREADS_MAKE_TREE == 256
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
#undef  REGISTERS_PER_THREAD_LINK_TREE
#define REGISTERS_PER_THREAD_LINK_TREE (23)
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
#endif//GPUVER == 52

/* makeTree_kernel use 47 registers @ Tesla P100 w/z Ttot = 256 @ CUDA 8.0 */
/* linkTree_kernel use 25 registers @ Tesla P100 w/z Ttot = 256 @ CUDA 8.0 */
#   if  GPUVER == 60
#undef  REGISTERS_PER_THREAD_MAKE_TREE
#define REGISTERS_PER_THREAD_MAKE_TREE (47)
#undef  REGISTERS_PER_THREAD_LINK_TREE
#define REGISTERS_PER_THREAD_LINK_TREE (25)
#endif//GPUVER == 60

/** limitations from number of registers */
#   if  NBLOCKS_PER_SM_MAKE_TREE > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAKE_TREE * NTHREADS_MAKE_TREE))
#undef  NBLOCKS_PER_SM_MAKE_TREE
#define NBLOCKS_PER_SM_MAKE_TREE   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAKE_TREE * NTHREADS_MAKE_TREE))
#endif//NBLOCKS_PER_SM_MAKE_TREE > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAKE_TREE * NTHREADS_MAKE_TREE))
#   if  NBLOCKS_PER_SM_LINK_TREE > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_LINK_TREE * NTHREADS_LINK_TREE))
#undef  NBLOCKS_PER_SM_LINK_TREE
#define NBLOCKS_PER_SM_LINK_TREE   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_LINK_TREE * NTHREADS_LINK_TREE))
#endif//NBLOCKS_PER_SM_LINK_TREE > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_LINK_TREE * NTHREADS_LINK_TREE))

/** maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_MAKE_TREE > (DIV_NTHREADS_MAKE_TREE(MAX_THREADS_PER_SM))
#undef  NBLOCKS_PER_SM_MAKE_TREE
#define NBLOCKS_PER_SM_MAKE_TREE   (DIV_NTHREADS_MAKE_TREE(MAX_THREADS_PER_SM))
#endif//NBLOCKS_PER_SM_MAKE_TREE > (DIV_NTHREADS_MAKE_TREE(MAX_THREADS_PER_SM))
#   if  NBLOCKS_PER_SM_LINK_TREE > (DIV_NTHREADS_LINK_TREE(MAX_THREADS_PER_SM))
#undef  NBLOCKS_PER_SM_LINK_TREE
#define NBLOCKS_PER_SM_LINK_TREE   (DIV_NTHREADS_LINK_TREE(MAX_THREADS_PER_SM))
#endif//NBLOCKS_PER_SM_LINK_TREE > (DIV_NTHREADS_LINK_TREE(MAX_THREADS_PER_SM))

/** maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_MAKE_TREE > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_MAKE_TREE
#define NBLOCKS_PER_SM_MAKE_TREE   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_MAKE_TREE > MAX_BLOCKS_PER_SM
#   if  NBLOCKS_PER_SM_LINK_TREE > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_LINK_TREE
#define NBLOCKS_PER_SM_LINK_TREE   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_LINK_TREE > MAX_BLOCKS_PER_SM

/** maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_MAKE_TREE > (DIV_NTHREADS_MAKE_TREE(MAX_WARPS_PER_SM * 32))
#undef  NBLOCKS_PER_SM_MAKE_TREE
#define NBLOCKS_PER_SM_MAKE_TREE   (DIV_NTHREADS_MAKE_TREE(MAX_WARPS_PER_SM * 32))
#endif//NBLOCKS_PER_SM_MAKE_TREE > (DIV_NTHREADS_MAKE_TREE(MAX_WARPS_PER_SM * 32))
#   if  NBLOCKS_PER_SM_LINK_TREE > (DIV_NTHREADS_LINK_TREE(MAX_WARPS_PER_SM * 32))
#undef  NBLOCKS_PER_SM_LINK_TREE
#define NBLOCKS_PER_SM_LINK_TREE   (DIV_NTHREADS_LINK_TREE(MAX_WARPS_PER_SM * 32))
#endif//NBLOCKS_PER_SM_LINK_TREE > (DIV_NTHREADS_LINK_TREE(MAX_WARPS_PER_SM * 32))

/** # of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_MAKE_TREE < 1
#undef  NBLOCKS_PER_SM_MAKE_TREE
#define NBLOCKS_PER_SM_MAKE_TREE  (1)
#endif//NBLOCKS_PER_SM_MAKE_TREE < 1
#   if  NBLOCKS_PER_SM_LINK_TREE < 1
#undef  NBLOCKS_PER_SM_LINK_TREE
#define NBLOCKS_PER_SM_LINK_TREE  (1)
#endif//NBLOCKS_PER_SM_LINK_TREE < 1


#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)

#define NBLOCKS_PER_SM_OUTFLOW (2)

#ifdef  TIME_BASED_MODIFICATION
#define REGISTERS_PER_THREAD_OUTFLOW (32)
#else///TIME_BASED_MODIFICATION
#ifdef  USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES
#define REGISTERS_PER_THREAD_OUTFLOW (35)
#else///USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES
#define REGISTERS_PER_THREAD_OUTFLOW (31)
#endif//USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES
#endif//TIME_BASED_MODIFICATION

/** limitations from number of registers */
#   if  NBLOCKS_PER_SM_OUTFLOW > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_OUTFLOW * NTHREADS_OUTFLOW))
#undef  NBLOCKS_PER_SM_OUTFLOW
#define NBLOCKS_PER_SM_OUTFLOW   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_OUTFLOW * NTHREADS_OUTFLOW))
#endif//NBLOCKS_PER_SM_OUTFLOW > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_OUTFLOW * NTHREADS_OUTFLOW))

/** maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_OUTFLOW > (DIV_NTHREADS_OUTFLOW(MAX_THREADS_PER_SM))
#undef  NBLOCKS_PER_SM_OUTFLOW
#define NBLOCKS_PER_SM_OUTFLOW   (DIV_NTHREADS_OUTFLOW(MAX_THREADS_PER_SM))
#endif//NBLOCKS_PER_SM_OUTFLOW > (DIV_NTHREADS_OUTFLOW(MAX_THREADS_PER_SM))

/** maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_OUTFLOW > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_OUTFLOW
#define NBLOCKS_PER_SM_OUTFLOW   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_OUTFLOW > MAX_BLOCKS_PER_SM

/** maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_OUTFLOW > (DIV_NTHREADS_OUTFLOW(MAX_WARPS_PER_SM * 32))
#undef  NBLOCKS_PER_SM_OUTFLOW
#define NBLOCKS_PER_SM_OUTFLOW   (DIV_NTHREADS_OUTFLOW(MAX_WARPS_PER_SM * 32))
#endif//NBLOCKS_PER_SM_OUTFLOW > (DIV_NTHREADS_OUTFLOW(MAX_WARPS_PER_SM * 32))

/** number of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_OUTFLOW < 1
#undef  NBLOCKS_PER_SM_OUTFLOW
#define NBLOCKS_PER_SM_OUTFLOW  (1)
#endif//NBLOCKS_PER_SM_OUTFLOW < 1
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)

#endif//USE_OCCUPANCY_CALCULATOR


#   if  defined(GADGET_MAC) || defined(WS93_MAC)
__constant__ real mac_delta;/**< G / alpha or sqrt(3 / accErr) */
#endif//defined(GADGET_MAC) || defined(WS93_MAC)
__constant__  treecell null_cell_device;
__constant__ jparticle zero_pj;
__constant__ jmass     zero_mj;


/**
 * @union alignedInt
 *
 * @brief union for store multiple integers in registers
 */
/**
 * @union alignedFlt
 *
 * @brief union for store multiple floats in registers
 */
#   if  NBUF_MAC == 4
typedef union __align__(16)
{
  int4 i4;
  int  ia[4];
} alignedInt;
typedef union __align__(16)
{
  real4 r4;
  real  ra[4];
} alignedFlt;
#endif//NBUF_MAC == 4
#   if  NBUF_MAC == 2
typedef union __align__(8)
{
  int2 i4;
  int  ia[2];
} alignedInt;
typedef union __align__(8)
{
  real2 r4;
  real  ra[2];
} alignedFlt;
#endif//NBUF_MAC == 2


/**
 * @union int_real
 *
 * @brief union for switching integer and float
 */
typedef union
{
  int  i;
  real r;
} int_real;


/**
 * @fn initPHhierarchy_kernel
 *
 * @brief Initialize information on hierarchy of Peano--Hilbert key.
 */
__global__ void initPHhierarchy_kernel(PHinfo *level)
{
  const int ii = GLOBALIDX_X1D;

  if( ii < NUM_PHKEY_LEVEL ){
    /** load current information on PH-key hierarchy */
    PHinfo phLev = level[ii];

    /** initialize PH-key information */
    phLev.head = NULL_CELL;
    phLev.num  = 0;

    /** set a root cell */
    if( ii == 0 ){
      phLev.head = 0;
      phLev.num  = 1;
    }/* if( ii == 0 ){ */

    /** store initialized information on PH-key hierarchy */
    level[ii] = phLev;
  }/* if( ii < NUM_PHKEY_LEVEL ){ */
}


/**
 * @fn initTreeCellOffset_kernel
 *
 * @brief Initialize information on tree cells.
 */
__global__ void initTreeCellOffset_kernel
(const int cellHead, const int cellNum, treecell * RESTRICT cell, PHint * RESTRICT hkey, uint * RESTRICT parent, uint * RESTRICT children, bool * RESTRICT leaf, uint * RESTRICT ptag, const int piNum)
{
  const int ii = cellHead + GLOBALIDX_X1D;

  /** initialize tree-cell information */
  treecell tmp_null_cell = {0, NULL_CELL};
  int tmp_hkey = -1;

  /** set a root cell */
  if( ii == 0 ){
    tmp_null_cell.head = 0;
    tmp_null_cell.num  = piNum;
    tmp_hkey           = 0;
  }/* if( ii == 0 ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

  if( ii < cellNum ){
    /** store initialized information on tree-cell */
    cell    [ii] = tmp_null_cell;
    hkey    [ii] = tmp_hkey;
    parent  [ii] = NULL_CELL;
    children[ii] = NULL_CELL;
    leaf    [ii] = true;
    ptag    [ii] = NULL_NODE;
  }/* if( ii < cellNum ){ */
}


/**
 * @fn initTreeNode_kernel
 *
 * @brief Initialize information on pseudo particles.
 */
/* __global__ void initTreeNode_kernel(const int pjNum, uint * RESTRICT more, int * RESTRICT node2cell) */
__global__ void initTreeNode_kernel(uint * RESTRICT more, int * RESTRICT node2cell)
{
  const int jj = GLOBALIDX_X1D;
  /* if( jj < pjNum ){ */
  /*   more     [jj] = NULL_NODE; */
  /*   node2cell[jj] = NULL_CELL; */
  /* }/\* if( jj < pjNum ){ *\/ */
  more     [jj] = NULL_NODE;
  node2cell[jj] = NULL_CELL;
}


/**
 * @fn initTreeLink_kernel
 *
 * @brief Initialize information on relation between i-particle and real j-particle.
 */
/* __global__ void initTreeLink_kernel(const int piNum, int *jtag) */
__global__ void initTreeLink_kernel(int *jtag)
{
  /* const int ii = GLOBALIDX_X1D; */
  /* if( ii < piNum ) */
  /*   jtag[ii] = NULL_NODE; */
  jtag[GLOBALIDX_X1D] = NULL_NODE;
}


/**
 * @fn makeTree_kernel
 *
 * @brief Make octree structure in a width-first manner.
 */
#define NTHREADS_SCAN_INC NTHREADS_MAKE_TREE
#ifdef  USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
#define USE_WARP_SHUFFLE_FUNC_SCAN_INC
#endif//USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
#include "../util/scan_inc.cu"
__global__ void
#ifndef USE_OCCUPANCY_CALCULATOR
__launch_bounds__(NTHREADS_MAKE_TREE, NBLOCKS_PER_SM_MAKE_TREE)
#endif//USE_OCCUPANCY_CALCULATOR
makeTree_kernel
     (int * RESTRICT leafLev_gm, PHinfo * RESTRICT level,
      int * RESTRICT numCell_gm, treecell * RESTRICT cell, PHint * RESTRICT hkey, bool * RESTRICT leaf, uint * RESTRICT children, uint * RESTRICT parent,
      READ_ONLY PHint * RESTRICT peano,
      int * RESTRICT gmem, volatile int * RESTRICT scanNum_gm
      , int * RESTRICT gsync0Ful, int * RESTRICT gsync1Ful, int * RESTRICT gsync0Loc, int * RESTRICT gsync1Loc
      )
{
  /** identify thread properties */
  const int bnum =   GRIDDIM_X1D;
  const int tidx = THREADIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int gnum = NGROUPS_MAKE_TREE * bnum;
  const int lane = tidx & (TSUB_MAKE_TREE - 1);/**< index of the thread within a thread group */
  const int scanLane = tidx & (warpSize - 1);

  const int ghead = tidx - lane;
  const int gtail = ghead + (TSUB_MAKE_TREE - 1);

  const int  idx = DIV_TSUB_MAKE_TREE(tidx);

#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  thread_block_tile<TSUB_MAKE_TREE> tile = tiled_partition<TSUB_MAKE_TREE>(this_thread_block());
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

  /** initialize number of tree cells */
  int numCell = 1;

  /** shared quantities in the thread parallelized version */
  __shared__ int smem[NTHREADS_MAKE_TREE];
  __shared__ treecell cell_sm[NTHREADS_MAKE_TREE];
  __shared__ PHint    hkey_sm[NTHREADS_MAKE_TREE];

  __shared__ PHinfo lev_sm[NUM_PHKEY_LEVEL];
  if( tidx < NUM_PHKEY_LEVEL )
    lev_sm[tidx] = level[tidx];
  __syncthreads();


  /** make tree structure in a width-first manner */
  for(int levelIdx = 0; levelIdx < MAXIMUM_PHKEY_LEVEL; levelIdx++){
    /** set level of tree cells to examine in this procedure */
    const PHinfo lev = lev_sm[levelIdx];
    const int cellLev  = lev.level;
    const int cellHead = lev.head;
    int cellNum  = lev.num;

    /** set keylevel of Peano-Hilbert key hierarchy about child cell(s) */
    /** PH key level is common within the threads */
    const   int leafLevel = cellLev - 1;
    const PHint leafScale = (PHint)1 << (leafLevel * 3);

    PHinfo daughter = lev_sm[levelIdx + 1];
#ifdef  SPLIT_CHILD_CELL
    const   int leafNmax = daughter.nmax;
#endif//SPLIT_CHILD_CELL

    /** global properties of tree cells belong to the same level hierarchy */
    int cellTail;  /**< equal to cellHead of tree cells belong to the lower level */
    int     tail;  /**< equal to the head index of the write buffer for child cells */
    tail = cellTail = cellHead + cellNum;


    /** examine tree cells */
    const int Niter = BLOCKSIZE(cellNum, gnum);
    for(int iter = 0; iter < Niter; iter++){
      const int cnumSub = (gnum < cellNum) ? gnum : cellNum;
      const int bnumSub = BLOCKSIZE(cnumSub, NGROUPS_MAKE_TREE);/**< number of active blocks (required to global synchronization) */
      if( bidx < bnumSub ){
	/** load a tree cell and evaluate either the cell has child cell(s) or no */
	/** coalesced load from global memory to shared memory */
	int cidx = cellHead + gnum * iter + bidx * NGROUPS_MAKE_TREE + tidx;
	if( (tidx < NGROUPS_MAKE_TREE) && (cidx < cellTail) ){
	  cell_sm[tidx] = cell[cidx];
	  hkey_sm[tidx] = hkey[cidx];
	}/* if( (tidx < NGROUPS_MAKE_TREE) && (cidx < cellTail) ){ */
	__syncthreads();

	cidx = cellHead + gnum * iter + bidx * NGROUPS_MAKE_TREE + idx;/**< common in TSUB_MAKE_TREE (= CELL_UNIT = 8) threads */
	const treecell root = (cidx < cellTail) ? (cell_sm[idx]) : (null_cell_device);/**< common in TSUB_MAKE_TREE (= CELL_UNIT = 8) threads */
	const PHint keyHead = (cidx < cellTail) ? (hkey_sm[idx]) : ((PHint)(-1));/**< common in TSUB_MAKE_TREE (= CELL_UNIT = 8) threads */

	const bool node = (root.num > NCRIT) ? (true) : (false);/**< common in TSUB_MAKE_TREE (= CELL_UNIT = 8) threads */

	/** divide the responsible tree cell if the cell is not a leaf cell */
	__syncthreads();
	cell_sm[tidx].head = 0;
	cell_sm[tidx].num = root.num;
	int addChild = 0;
	if( node ){
	  int ihead = 0;
	  int itail = root.num - 1;

	  /** the responsible tree cell is a node cell */
	  if( lane == 0 ){
	    leaf[cidx] = false;
	    hkey_sm[tidx] = peano[root.head + itail];
	  }/* if( lane == 0 ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  tile.sync();/**< tile.sync() to reduce warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

	  /** find CELL_UNIT - 1 (== TSUB_MAKE_TREE - 1 == 7) boundaries (in maximum) */
	  if( lane < (CELL_UNIT - 1) ){
	    /** properties of the PH-key of the focused child cell */
	    /** NOTE: the child cell may not exist */
	    const PHint ktail = keyHead + (1 + lane) * leafScale;

	    /** check whether the specified region contains keys between khead and ktail */
	    if( ktail <= hkey_sm[ghead] ){
	      /** find the tail of the PH-key for the child cell using bisection method */
	      while( true ){
		/** when the procedure finds the tail of the PH-key... */
		if( itail <= (1 + ihead) )		  break;

		uint ihalf = (uint)(ihead + itail) >> 1;
		if( peano[root.head + ihalf] <= ktail )		  ihead = (int)ihalf;
		else		                                  itail = (int)ihalf;
	      }/* while( true ){ */

	      itail = ihead;
	    }/* if( ktail <= hkey_sm[ghead] ){ */

	    if( ((peano[root.head] - keyHead) / leafScale) <= lane )
	      itail++;
	    cell_sm[tidx].num = itail;
	  }/* if( lane < (CELL_UNIT - 1) ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  tile.sync();/**< tile.sync() for consistency */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

	  /** assume implicit synchronization within a warp */
	  if( lane > 0 )
	    cell_sm[tidx].head = cell_sm[tidx - 1].num;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  tile.sync();/**< tile.sync() to reduce warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  cell_sm[tidx].num -= cell_sm[tidx].head;
	  cell_sm[tidx].head += root.head;

#ifdef  SPLIT_CHILD_CELL
	  addChild = (cell_sm[tidx].num > 0) ? BLOCKSIZE(cell_sm[tidx].num, leafNmax) : 0;
#else///SPLIT_CHILD_CELL
	  if( cell_sm[tidx].num > 0 )
	    addChild = 1;
#endif//SPLIT_CHILD_CELL
	}/* if( node ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	__syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

	/** calculate prefix sum about addChild */
	int headIdx, scanNum;
	int targetIdx = PREFIX_SUM_GRID_WITH_PARTITION(addChild, smem, scanLane, tidx, &headIdx, &scanNum, gmem, bidx, bnumSub, gsync0Loc, gsync1Loc);
	if( bidx + tidx == 0 )
	  *scanNum_gm = scanNum;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	__syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

	targetIdx += tail - addChild;/**< this must be an exclusive scan */

	const int tag = smem[gtail] - ((ghead > 0) ? (smem[ghead - 1]) : (headIdx));


	/** store child cells to the global memory */
	/** NOTE: tentative implementation (uncoalesced store to global memory) */
	/** if performance is too low due to uncoalesced store, then try to implement coalesced version (probably, use additional shared memory to stock splitted child cells) */
	if( (node) && (lane == 0) )
	  children[cidx] = ((tag - 1) << IDXBITS) + targetIdx;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	__syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

#ifdef  SPLIT_CHILD_CELL
	treecell childCell = cell_sm[tidx];
	int nrem = childCell.num;
	const int nchild = BLOCKSIZE(nrem, addChild);
	for(int childIdx = 0; childIdx < addChild; childIdx++){
	  /** estimate the number of particles per child cell */
	  childCell.num = (nchild < nrem) ? nchild : nrem;

	  /** commit the new child cell */
	  cell  [targetIdx + childIdx] = childCell;
	  hkey  [targetIdx + childIdx] = keyHead + lane * leafScale;
	  parent[targetIdx + childIdx] = cidx;

	  /** preparation to the next child cell */
	  childCell.head += childCell.num;
	}/* for(int childIdx = 0; childIdx < addChild; childIdx++){ */
#else///SPLIT_CHILD_CELL
	/** commit the new child cell */
	if( addChild > 0 ){
	  cell  [targetIdx] = cell_sm[tidx];
	  hkey  [targetIdx] = keyHead + lane * leafScale;
	  parent[targetIdx] = cidx;
	}/* if( addChild > 0 ){ */
#endif//SPLIT_CHILD_CELL
      }/* if( bidx < bnumSub ){ */

      /** global synchronization within bnum blocks */
      globalSync(tidx, bidx, bnum, gsync0Ful, gsync1Ful);

      if( tidx == 0 )
	smem[0] = *scanNum_gm;
      __syncthreads();
      const  int addCellNum = smem[0];
      tail    += addCellNum;
      numCell += addCellNum;
    }/* for(int iter = 0; iter < Niter; iter++){ */

    daughter.head =        cellTail;
    daughter.num  = tail - cellTail;
    if( tidx == 0 )
      lev_sm[levelIdx + 1] = daughter;

    if( daughter.num == 0 ){
      __syncthreads();
      break;
    }/* if( daughter.num == 0 ){ */

    globalSync(tidx, bidx, bnum, gsync0Ful, gsync1Ful);
  }/* for(int levelIdx = 0; levelIdx < MAXIMUM_PHKEY_LEVEL; levelIdx++){ */


  if( bidx == 0 ){
    if( tidx < NUM_PHKEY_LEVEL )
      level[tidx] = lev_sm[tidx];

    if( tidx == NUM_PHKEY_LEVEL )
      *numCell_gm = numCell;

    if( tidx == 0 ){
      int leafLev = NUM_PHKEY_LEVEL - 1;
      for(int levelIdx = (NUM_PHKEY_LEVEL - 1); levelIdx >= 0; levelIdx--)
	if( lev_sm[levelIdx].num != 0 ){
	  leafLev = levelIdx;
	  break;
	}/* if( lev_sm[levelIdx].num != 0 ){ */
      *leafLev_gm = leafLev;
    }/* if( tidx == 0 ){ */
  }/* if( bidx == 0 ){ */
}


/**
 * @fn linkNode
 *
 * @brief Extend the pseudo-particle chain.
 */
#   if  NTHREADS_SCAN_INC != NTHREADS_LINK_TREE
#undef	NTHREADS_SCAN_INC
#define NTHREADS_SCAN_INC    NTHREADS_LINK_TREE
#include "../util/scan_inc.cu"
#endif//NTHREADS_SCAN_INC != NTHREADS_LINK_TREE
__device__ __forceinline__ void linkNode
(const treecell root, const bool leaf_cell, uint * RESTRICT ptag,
 uint * RESTRICT more, int * RESTRICT jtag, int * RESTRICT phead,
 const int bnum, const int bidx, const int tidx, const int lane,
 volatile int * RESTRICT smem, int * RESTRICT gmem, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  /** load fundamental information on the focusing tree cell */
  const int head = root.head;
  /** below procedure is switched on if head != NULL_CELL */
  /** if the cell is a leaf cell, then root.num pseudo particles are set. otherwise, a pseudo particle is set.  */
  const int nadd = (head != NULL_CELL) ? ((leaf_cell) ? (root.num) : (1)) : (0);

  /** calculate prefix sum about addChild */
  int scanNum;
  int headIdx = PREFIX_SUM_GRID(nadd, smem, lane, tidx, &scanNum, gmem, bidx, bnum, gsync0, gsync1);
  headIdx += (*phead) - nadd;/**< this must be an exclusive scan */

  if( nadd > 0 )
    *ptag = ((nadd - 1) << IDXBITS) + headIdx;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

  /** when the tree cell is a leaf cell, then more index specifies the cell itself */
  if( (head != NULL_CELL) && leaf_cell ){
    for(int jj = 0; jj < nadd; jj++){
      /** commit the new particle to the j-particle array */
      const int jidx = headIdx + jj;
      more[jidx] = jidx;

      /** connect an i-particle with the corresponding real j-particle */
      jtag[head + jj] = jidx;
    }/* for(int jj = 0; jj < nadd; jj++){ */
  }/* if( (head != NULL_CELL) && lear_cell ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

  *phead += scanNum;
}


/**
 * @fn linkTree_kernel
 *
 * @brief Link tree nodes.
 */
__global__ void
#ifndef USE_OCCUPANCY_CALCULATOR
__launch_bounds__(NTHREADS_LINK_TREE, NBLOCKS_PER_SM_LINK_TREE)
#endif//USE_OCCUPANCY_CALCULATOR
linkTree_kernel
     (const int cellTail, READ_ONLY treecell * RESTRICT cell, READ_ONLY bool * RESTRICT leaf, uint * RESTRICT ptag,
      int * RESTRICT pjNum, uint * RESTRICT more, int * RESTRICT jtag,
      int * RESTRICT gmem, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  /** identify thread properties */
  const int bnum =   GRIDDIM_X1D;
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int lane = tidx & (warpSize - 1);

  /** shared quantities in the thread parallelized version */
  __shared__ int smem[NTHREADS_LINK_TREE];

  /** extend the pseudo particle chain */
  /** examine tree cells */
  int cellNum = cellTail;
  int numNode = 0;

  const int chkNumStep = NTHREADS_LINK_TREE * bnum;
  const int Nloop = BLOCKSIZE(cellNum, chkNumStep);
  for(int loop = 0; loop < Nloop; loop++){
    const int chkNum = (chkNumStep < cellNum) ? chkNumStep : cellNum;

    /** load a tree cell and evaluate either the cell has child cell(s) or no */
    const int cidx = gidx + chkNumStep * loop;
    treecell root = (cidx < cellTail) ? (cell[cidx]) : (null_cell_device);

    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    /** extend the pseudo particle chain */
    linkNode(root, leaf[cidx], &ptag[cidx], more, jtag, &numNode,
	     bnum, bidx, tidx, lane, smem, gmem, gsync0, gsync1);

    cellNum -= chkNum;
  }/* for(int loop = 0; loop < Nloop; loop++){ */

  /** return # of pseudo j-particles */
  if( gidx == 0 )
    *pjNum = numNode;
}


/**
 * @fn trimTree_kernel
 *
 * @brief Extend the pseudo-particle chain.
 */
__global__ void trimTree_kernel
(const int cellHead, const int cellTail, READ_ONLY treecell * RESTRICT cell, READ_ONLY bool * RESTRICT leaf, READ_ONLY uint * RESTRICT children, READ_ONLY uint * RESTRICT ptag,
 int * RESTRICT node2cell, uint * RESTRICT more)
{
  const int cidx = cellHead + GLOBALIDX_X1D;

  /** connect tree nodes in a different PH-level */
  if( (cidx < cellTail) && (cell[cidx].head != NULL_CELL) ){
    /** relate the tree cell and corresponding tree nodes */
    const uint node = ptag[cidx];
    const int pidx = node & IDXMASK;
#ifdef  SPLIT_CHILD_CELL
    const int pnum = 1 + (int)(node >> IDXBITS);
    for(int jj = 0; jj < pnum; jj++)
      node2cell[pidx + jj] = cidx;
#else///SPLIT_CHILD_CELL
    node2cell[pidx] = cidx;
#endif//SPLIT_CHILD_CELL

    if( leaf[cidx] == false ){
      /** count up number of child nodes */
      const uint child = children[cidx];
      const  int chead =               child  & IDXMASK;
      const  int ctail = chead + (1 + (child >> IDXBITS));
      const uint phead = ptag[chead];
      int childNum = 1 + (phead >> IDXBITS);
      for(int jj = chead + 1; jj < ctail; jj++)
	childNum += (1 + (ptag[jj] >> IDXBITS));

      /** commit child particles as pseudo particle */
      more[pidx] = ((childNum - 1) << IDXBITS) + (phead & IDXMASK);
    }/* if( leaf[cidx] == false ){ */
  }/* if( (cidx < cellTail) && (cell[cidx].head != NULL_CELL) ){ */
}


/**
 * @fn makeTreeStructure_dev
 *
 * @brief Make octree structure in a width-first manner.
 */
extern "C"
void makeTreeStructure_dev
(const int piNum, PHint * RESTRICT peano,
 int * RESTRICT leafLev, int * RESTRICT leafLev_dev, int * RESTRICT scanNum_dev,
 int * RESTRICT cellNum, int * RESTRICT cellNum_dev, const soaTreeCell cell,
 int * RESTRICT nodeNum, int * RESTRICT nodeNum_dev, const soaTreeNode node,
 const soaMakeTreeBuf buf, deviceProp devProp
#ifndef SERIALIZED_EXECUTION
 , const int piNum_prev
#endif//SERIALIZED_EXECUTION
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
)
{
  __NOTE__("%s\n", "start");


  /** start stop watch */
#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK

  /** initialize tree structure */
#ifdef  SERIALIZED_EXECUTION
  const int piNum_prev = piNum;
#endif//SERIALIZED_EXECUTION
  int Nrem = BLOCKSIZE(piNum_prev, NTHREADS_INIT_LINK);
  if( Nrem <= MAX_BLOCKS_PER_GRID ){
    /* initTreeLink_kernel<<<Nrem, NTHREADS_INIT_LINK>>>(piNum_prev, node.jtag); */
    initTreeLink_kernel<<<Nrem, NTHREADS_INIT_LINK>>>(node.jtag);
  }/* if( Nrem <= MAX_BLOCKS_PER_GRID ){ */
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_INIT_LINK;
      /* initTreeLink_kernel<<<Nblck, NTHREADS_INIT_LINK>>>(Nsub, &node.jtag[hidx]); */
      initTreeLink_kernel<<<Nblck, NTHREADS_INIT_LINK>>>(&node.jtag[hidx]);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->initTreeLink_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER

  Nrem = BLOCKSIZE(*cellNum, NTHREADS_INIT_CELL);
  __NOTE__("cellNum in the previous step: %d\n", *cellNum);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    initTreeCellOffset_kernel<<<Nrem, NTHREADS_INIT_CELL>>>(0, *cellNum, cell.cell, cell.hkey, cell.parent, cell.children, cell.leaf, cell.ptag, piNum);
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_INIT_CELL;
      initTreeCellOffset_kernel<<<Nblck, NTHREADS_INIT_CELL>>>(hidx, *cellNum, cell.cell, cell.hkey, cell.parent, cell.children, cell.leaf, cell.ptag, piNum);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->initTreeCell_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER

  Nrem = BLOCKSIZE(*nodeNum, NTHREADS_INIT_NODE);
  __NOTE__("launch initTreeNode_kernel\n");
  if( Nrem <= MAX_BLOCKS_PER_GRID ){
    /* initTreeNode_kernel<<<Nrem, NTHREADS_INIT_NODE>>>(*nodeNum, node.more, node.node2cell); */
    initTreeNode_kernel<<<Nrem, NTHREADS_INIT_NODE>>>(node.more, node.node2cell);
  }/* if( Nrem <= MAX_BLOCKS_PER_GRID ){ */
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_INIT_NODE;
      /* initTreeNode_kernel<<<Nblck, NTHREADS_INIT_NODE>>>(Nsub, &node.more[hidx], &node.node2cell[hidx]); */
      initTreeNode_kernel<<<Nblck, NTHREADS_INIT_NODE>>>(&node.more[hidx], &node.node2cell[hidx]);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->initTreeNode_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER

  __NOTE__("launch initPHhierarchy_kernel\n");
  initPHhierarchy_kernel<<<1, NUM_PHKEY_LEVEL>>>(cell.level);

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->makeTree));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER


  /** make tree structure */
#ifdef  USE_OCCUPANCY_CALCULATOR
  const int NBLOCKS_PER_SM_MAKE_TREE = buf.numBlocksPerSM_make_tree;
#endif//USE_OCCUPANCY_CALCULATOR
  __NOTE__("launch makeTree_kernel: %d blocks per SM, %d SMs\n", NBLOCKS_PER_SM_MAKE_TREE, devProp.numSM);
  makeTree_kernel<<<devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE, NTHREADS_MAKE_TREE>>>
    (leafLev_dev, cell.level, cellNum_dev, cell.cell, cell.hkey, cell.leaf, cell.children, cell.parent, peano,
     buf.gmem_make_tree, scanNum_dev, buf.gsync0_make_tree, buf.gsync1_make_tree, buf.gsync2_make_tree, buf.gsync3_make_tree);
  getLastCudaError("makeTree_kernel");

  checkCudaErrors(cudaMemcpy(leafLev, leafLev_dev, sizeof(int), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cellNum, cellNum_dev, sizeof(int), cudaMemcpyDeviceToHost));
  if( *cellNum > NUM_ALLOC_TREE_CELL ){
    __KILL__(stderr, "ERROR: allocated vector length for tree-cell arrays are too short, increase TREE_SAFETY_VAL(%f) defined in src/tree/make.h or use more MPI processes!\n", TREE_SAFETY_VAL);
  }/* if( *cellNum > NUM_ALLOC_TREE_CELL ){ */
#if 0
  fprintf(stdout, "*cellNum = %d\n", *cellNum);
  fflush(NULL);
  exit(0);
#endif
  __NOTE__("leafLev = %d\n", *leafLev);
  __NOTE__("cellNum in the current step: %d\n", *cellNum);

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->makeTree_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER


  /** link tree structure */
#ifdef  USE_OCCUPANCY_CALCULATOR
  const int NBLOCKS_PER_SM_LINK_TREE = buf.numBlocksPerSM_link_tree;
#endif//USE_OCCUPANCY_CALCULATOR
  __NOTE__("launch linkTree_kernel: %d blocks per SM, %d SMs\n", NBLOCKS_PER_SM_LINK_TREE, devProp.numSM);
  linkTree_kernel<<<devProp.numSM * NBLOCKS_PER_SM_LINK_TREE, NTHREADS_LINK_TREE>>>
    (*cellNum, cell.cell, cell.leaf, cell.ptag,
     nodeNum_dev, node.more, node.jtag,
     buf.gmem_link_tree, buf.gsync0_link_tree, buf.gsync1_link_tree);
  getLastCudaError("linkTree_kernel");
  checkCudaErrors(cudaMemcpy(nodeNum, nodeNum_dev, sizeof(int), cudaMemcpyDeviceToHost));

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->linkTree_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER


  /** trim tree structure */
  Nrem = BLOCKSIZE(*cellNum, NTHREADS_TRIM_TREE);
#if 0
  fprintf(stdout, "Nrem = %d, *cellNum = %d\n", Nrem, *cellNum);
  fflush(NULL);
  exit(0);
#endif
  __NOTE__("launch trimTree_kernel\n");
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    trimTree_kernel<<<Nrem, NTHREADS_TRIM_TREE>>>(0, *cellNum, cell.cell, cell.leaf, cell.children, cell.ptag, node.node2cell, node.more);
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_TRIM_TREE;
      trimTree_kernel<<<Nblck, NTHREADS_TRIM_TREE>>>(hidx, *cellNum, cell.cell, cell.leaf, cell.children, cell.ptag, node.node2cell, node.more);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */
  getLastCudaError("trimTree_kernel");

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->trimTree_kernel));
#endif//HUNT_MAKE_PARAMETER


  /** stop stop watch */
#ifdef  EXEC_BENCHMARK
#ifndef HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->makeTree));
#else///HUNT_MAKE_PARAMETER
  elapsed->makeTree += elapsed->initTreeLink_kernel + elapsed->initTreeCell_kernel + elapsed->initTreeNode_kernel + elapsed->makeTree_kernel + elapsed->linkTree_kernel + elapsed->trimTree_kernel;
#endif//HUNT_MAKE_PARAMETER
#endif//EXEC_BENCHMARK


  __NOTE__("%s\n", "end");
}


/**
 * @fn allocTreeCell_dev
 *
 * @brief Allocate arrays to store properties of tree cells.
 */
extern "C"
muse allocTreeCell_dev
(soaTreeCell *dev, treecell **cell_dev, bool **leaf_dev, uint **node_dev, PHinfo **info_dev,
 PHint **hkey_dev, uint **parent_dev, uint **children_dev, int **leafLev_dev, int **numCell_dev, int **numNode_dev, int **scanNum_dev
#ifdef  COUNT_INTERACTIONS
 , soaTreeCell *hst, treecell **cell_hst, bool **leaf_hst, uint **node_hst, PHinfo **info_hst
#endif//COUNT_INTERACTIONS
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
 , deviceProp devProp
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
 )
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  const size_t num = NUM_ALLOC_TREE_CELL;

  mycudaMalloc((void **)    cell_dev, num * sizeof(treecell));  alloc.device += num * sizeof(treecell);
  mycudaMalloc((void **)    leaf_dev, num * sizeof(    bool));  alloc.device += num * sizeof(    bool);
  mycudaMalloc((void **)    node_dev, num * sizeof(    uint));  alloc.device += num * sizeof(    uint);
  mycudaMalloc((void **)    hkey_dev, num * sizeof(   PHint));  alloc.device += num * sizeof(   PHint);
  mycudaMalloc((void **)  parent_dev, num * sizeof(    uint));  alloc.device += num * sizeof(    uint);
  mycudaMalloc((void **)children_dev, num * sizeof(    uint));  alloc.device += num * sizeof(    uint);
  mycudaMalloc((void **) leafLev_dev, num * sizeof(     int));  alloc.device +=       sizeof(     int);
  mycudaMalloc((void **) numCell_dev, num * sizeof(     int));  alloc.device +=       sizeof(     int);
  mycudaMalloc((void **) numNode_dev, num * sizeof(     int));  alloc.device +=       sizeof(     int);
  mycudaMalloc((void **) scanNum_dev, num * sizeof(     int));  alloc.device +=       sizeof(     int);

#ifdef  COUNT_INTERACTIONS
  mycudaMallocHost((void **)cell_hst, num * sizeof(treecell));  alloc.host += num * sizeof(treecell);
  mycudaMallocHost((void **)leaf_hst, num * sizeof(    bool));  alloc.host += num * sizeof(    bool);
  mycudaMallocHost((void **)node_hst, num * sizeof(    uint));  alloc.host += num * sizeof(    uint);
#endif//COUNT_INTERACTIONS

#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  mycudaMalloc    ((void **)info_dev, (NUM_PHKEY_LEVEL + 1) * sizeof(PHinfo));  alloc.device += (NUM_PHKEY_LEVEL + 1) * sizeof(PHinfo);
  /** enforce BLOCKSIZE(clev.num, bnum * NGROUPS_OUTFLOW) becomes zero */
  PHinfo info_tmp = {0, 1 - devProp.numSM * NBLOCKS_PER_SM_OUTFLOW * NGROUPS_OUTFLOW, 0, 0};
  checkCudaErrors(cudaMemcpy(&((*info_dev)[NUM_PHKEY_LEVEL]), &info_tmp, sizeof(PHinfo), cudaMemcpyHostToDevice));
#else///!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  mycudaMalloc    ((void **)info_dev,  NUM_PHKEY_LEVEL      * sizeof(PHinfo));  alloc.device +=  NUM_PHKEY_LEVEL      * sizeof(PHinfo);
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)

  initPHinfo_dev(*info_dev);

#ifdef  COUNT_INTERACTIONS
  mycudaMallocHost((void **)info_hst, NUM_PHKEY_LEVEL * sizeof(PHinfo));  alloc.host   += NUM_PHKEY_LEVEL * sizeof(PHinfo);
#endif//COUNT_INTERACTIONS


  dev->cell  = *cell_dev;
  dev->leaf  = *leaf_dev;
  dev->ptag  = *node_dev;
  dev->level = *info_dev;
  dev->hkey     = *    hkey_dev;
  dev->parent   = *  parent_dev;
  dev->children = *children_dev;
#ifdef  COUNT_INTERACTIONS
  hst->cell  = *cell_hst;
  hst->leaf  = *leaf_hst;
  hst->ptag  = *node_hst;
  hst->level = *info_hst;
#endif//COUNT_INTERACTIONS


  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freeTreeCell_dev
 *
 * @brief Deallocate arrays to store properties of tree cells.
 */
extern "C"
void  freeTreeCell_dev
(treecell  *cell_dev, bool  *leaf_dev, uint  *node_dev, PHinfo  *info_dev,
 PHint  *hkey_dev, uint  *parent_dev, uint  *children_dev, int  *leafLev_dev, int  *numCell_dev, int  *numNode_dev, int  *scanNum_dev
#ifdef  COUNT_INTERACTIONS
 , treecell  *cell_hst, bool  *leaf_hst, uint  *node_hst, PHinfo  *info_hst
#endif//COUNT_INTERACTIONS
)
{
  __NOTE__("%s\n", "start");


  mycudaFree(cell_dev);
  mycudaFree(leaf_dev);
  mycudaFree(node_dev);
  mycudaFree(info_dev);
  mycudaFree(    hkey_dev);
  mycudaFree(  parent_dev);
  mycudaFree(children_dev);
  mycudaFree( leafLev_dev);
  mycudaFree( numCell_dev);
  mycudaFree( numNode_dev);
  mycudaFree( scanNum_dev);

#ifdef  COUNT_INTERACTIONS
  mycudaFreeHost(cell_hst);
  mycudaFreeHost(leaf_hst);
  mycudaFreeHost(node_hst);
  mycudaFreeHost(info_hst);
#endif//COUNT_INTERACTIONS


  __NOTE__("%s\n", "end");
}


#ifdef  GADGET_MAC
/**
 * @fn enforceBarnesHutMAC
 *
 * @brief Rewrite MAC for initial step (GADGET-MAC requires pre-estimated acceleration).
 *
 * @param (Ni) number of i-particles
 * @return (ai) acceleration and potential of i-particles
 * @param (Nj) number of j-particles
 * @return (pj) position and MAC of j-particles
 * @return (mac_bak) tentative stock of GADGET-MAC of all j-particles
 * @param (bmax) size of distribution of j-particles
 */
__global__ void enforceBarnesHutMAC
(const int Ni, acceleration * RESTRICT ai,
 const int Nj, jparticle * RESTRICT pj, real * RESTRICT mac_bak, READ_ONLY real * RESTRICT bmax)
{
  const int gidx = GLOBALIDX_X1D;

  if( gidx < Ni ){
    /** set old acceleration as |a| is unity */
    const acceleration iacc = {UNITY, ZERO, ZERO, ZERO};
    ai[gidx] = iacc;
  }/* if( gidx < Ni ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

  if( gidx < Nj ){
    /** load position and MAC of j-particle */
    jparticle jpos = pj[gidx];

    /** save the GADGET-MAC */
    mac_bak[gidx] = jpos.w;

    /** modify MAC for pseudo particles */
    /** squared size for the real particle is set to be negative (-UNITY) */
    if( jpos.w > -HALF ){
      /** new MAC for initial step is Barnes-Hut criterion */
      /** jcnd.w = bmax^4 / theta^4 in this case */
      const real bjmax = bmax[gidx];
      /** theta is assumed to be 0.5 for simplicity */
      jpos.w = TWO * bjmax;
      jpos.w *= jpos.w;
      jpos.w *= jpos.w;
      pj[gidx] = jpos;
    }/* if( jpos.w > -HALF ){ */
  }/* if( gidx < Nj ){ */
}


/**
 * @fn recoverGADGET_MAC
 *
 * @brief Recover MAC for further time steps.
 *
 * @param (Nj) number of j-particles
 * @return (pj) position and MAC of j-particles
 * @param (mac_bak) tentative stock of GADGET-MAC of all j-particles
 */
__global__ void recoverGADGET_MAC
(const int Nj, jparticle * RESTRICT pj, READ_ONLY real * RESTRICT mac_bak)
{
  const int jj = GLOBALIDX_X1D;

  if( jj < Nj ){
    /** load position and MAC of j-particle */
    jparticle jpos = pj[jj];
    const real mac = mac_bak[jj];

    /** recover the original MAC */
    jpos.w = mac;

    /** return the original position and MAC of j-particle */
    pj[jj] = jpos;
  }/* if( jj < Nj ){ */
}


/**
 * @fn enforceBarnesHutMAC_dev
 *
 * @brief Rewrite MAC for initial step (GADGET-MAC requires pre-estimated acceleration).
 *
 * @param (Ni) number of i-particles
 * @return (pi) i-particles
 * @param (Nj) number of j-particles
 * @return (pj) j-particles
 *
 * @sa enforceBarnesHutMAC
 */
extern "C"
void enforceBarnesHutMAC_dev(const int Ni, const iparticle pi, const int Nj, const soaTreeNode pj)
{
  __NOTE__("%s\n", "start");

  /** NOTE: Nj is always greater than Ni */
  int Njrem = BLOCKSIZE(Nj, 1024);
  int Nirem = BLOCKSIZE(Ni, 1024);
  if( Njrem <= MAX_BLOCKS_PER_GRID )
    enforceBarnesHutMAC<<<Njrem, 1024>>>(Ni, pi.acc, Nj, pj.jpos, pj.mac, pj.bmax);
  else{
    const int Niter = BLOCKSIZE(Njrem, MAX_BLOCKS_PER_GRID);
    int hjidx = 0;
    int hiidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Njblck = MAX_BLOCKS_PER_GRID;      if( Njblck > Njrem )	Njblck = Njrem;
      int Niblck = MAX_BLOCKS_PER_GRID;      if( Niblck > Nirem )	Niblck = Nirem;

      int Nisub = Niblck * 1024;
      int Njsub = Njblck * 1024;
      enforceBarnesHutMAC<<<Njblck, 1024>>>(Nisub, &pi.acc[hiidx], Njsub, &pj.jpos[hjidx], &pj.mac[hjidx], &pj.bmax[hjidx]);

      hjidx += Njsub;      Njrem -= Njblck;
      hiidx += Nisub;      Nirem -= Niblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("enforceBarnesHutMAC");

  __NOTE__("%s\n", "end");
}


/**
 * @fn recoverGADGET_MAC_dev
 *
 * @brief Recover MAC for further time steps.
 *
 * @param (Nj) number of j-particles
 * @return (pj) j-particles
 *
 * @sa recoverGADGET_MAC
 */
extern "C"
void recoverGADGET_MAC_dev(const int Nj, const soaTreeNode pj)
{
  __NOTE__("%s\n", "start");

  int Nrem = BLOCKSIZE(Nj, 1024);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    recoverGADGET_MAC<<<Nrem, 1024>>>(Nj, pj.jpos, pj.mac);
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * 1024;
      recoverGADGET_MAC<<<Nblck, 1024>>>(Nsub, &pj.jpos[hidx], &pj.mac[hidx]);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("recoverGADGET_MAC");

  __NOTE__("%s\n", "end");
}
#endif//GADGET_MAC



/**
 * @fn initTreeBody_kernel
 *
 * @brief Initialize mass and position of pseudo j-particles.
 *
 * @return (pj) position and MAC of j-particles
 * @return (mj) mass of j-particles
 * @return (bmax) size of j-particles
 * @return (mr2) quadrupole moment of j-particles
 */
__global__ void initTreeBody_kernel
(jparticle * RESTRICT pj, jmass * RESTRICT mj, real * RESTRICT bmax
#ifdef  WS93_MAC
 , real * RESTRICT mr2
#endif//WS93_MAC
)
{
  const int gidx = GLOBALIDX_X1D;

  pj  [gidx] = zero_pj;
  mj  [gidx] = zero_mj;
  bmax[gidx] = ZERO;
#ifdef  WS93_MAC
  mr2 [gidx] = ZERO;
#endif//WS93_MAC
}


/**
 * @fn copyRealBody_kernel
 *
 * @brief Set N-body particles as real j-particles.
 */
__global__ void copyRealBody_kernel
(const int Ni, READ_ONLY int * RESTRICT jtag, READ_ONLY position * RESTRICT pi,
 jparticle * RESTRICT pj, jmass * RESTRICT mj
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
 , const real eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  WS93_MAC
 , real * RESTRICT mr2
#endif//WS93_MAC
 )
{
  const int gidx = GLOBALIDX_X1D;

  if( gidx < Ni ){
    /** load an N-body particle */
    const position ipos = pi  [gidx];
    const int      jidx = jtag[gidx];

    /** set the N-body particle as a real particle */
    jparticle jpos;
    jpos.x = ipos.x;
    jpos.y = ipos.y;
    jpos.z = ipos.z;
    jpos.w = -UNITY;/**< squared size for the real particle is set to be negative */

#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
    const jmass mj_tmp = {ipos.m, eps2};
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
    const jmass mj_tmp =  ipos.m;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING

    /** commit the new particle to the j-particle array */
    pj[jidx] = jpos;
    mj[jidx] = mj_tmp;
#ifdef  WS93_MAC
    mr2[jidx] = ipos.m * (jpos.x * jpos.x + jpos.y * jpos.y + jpos.z * jpos.z);
#endif//WS93_MAC
  }/* if( gidx < Ni ){ */
}


#define TSUB_SCAN_INC TSUB_MAC
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#define USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#include "../util/scan_tsub_inc.cu"


#define TSUB_COMPARE_INC TSUB_MAC
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#define USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#include "../util/compare_tsub_inc.cu"


/**
 * @fn calcMultipole_kernel
 *
 * @brief Calculate multipole moment(s) of tree cells based on the width-first search.
 *
 * @param (bottomLev) bottom level of the PH-hierarchy
 * @param (level) head index and number of tree cells contained in the corresponding hierarchy
 * @param (cell) head index and number of N-body particles contained in the corresponding tree cell
 * @param (leaf) a flag to remember the corresponding tree cell is either leaf(true) of node(false)
 * @param (pi) position and mass of N-body particles
 * @param (node) head index and number of pseudo particles of the corresponding tree cell (a leaf cell contains up to NCRIT particles)
 * @param (more) head index and number of child particles of the corresponding j-particle (a leaf cell contains up to NCRIT particles)
 * @param (node2cell) index of the tree cell corresponding a pseudo particle
 * @return (pj) position and squared radius of pseudo N-body particle as j-particles
 * @return (mj) mass of pseudo N-body particle as j-particles
 * @return (mr2) mass times r squared of pseudo N-body particle as j-particles
 * @return (bmax) size of pseudo N-body particle as j-particles
 */
__global__ void
#ifndef USE_OCCUPANCY_CALCULATOR
__launch_bounds__(NTHREADS_MAC, NBLOCKS_PER_SM_MAC)
#endif//USE_OCCUPANCY_CALCULATOR
calcMultipole_kernel
     (const int bottomLev, READ_ONLY PHinfo * RESTRICT level,
      READ_ONLY treecell * RESTRICT cell, READ_ONLY bool * RESTRICT leaf, READ_ONLY position * RESTRICT pi,
      READ_ONLY uint * RESTRICT node, READ_ONLY uint * RESTRICT more, READ_ONLY int * RESTRICT node2cell,
      jparticle * RESTRICT pj, jmass * RESTRICT mj, real * RESTRICT bmax,
      int * RESTRICT more0Buf, int * RESTRICT more1Buf, real * RESTRICT rjmaxBuf, int * RESTRICT overflow
#ifndef USE_COOPERATIVE_GROUPS
      , int * RESTRICT gsync0, int * RESTRICT gsync1
#endif//USE_COOPERATIVE_GROUPS
#ifdef  WS93_MAC
      , real * RESTRICT mr2
#endif//WS93_MAC
#ifdef  COUNT_INTERACTIONS
      , tree_stats * RESTRICT stats
#endif//COUNT_INTERACTIONS
      )
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int bnum =   GRIDDIM_X1D;

  const int lane = tidx & (TSUB_MAC - 1);
  const int head = tidx - lane;
#ifndef USE_WARP_SHUFFLE_FUNC_MAC
  const int tail = head + TSUB_MAC - 1;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
  const int hbuf = DIV_TSUB_MAC(head) * TSUB_MAC * NBUF_MAC;/**< head index of the shared array close and queue within a thread group */

#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_MAC < 32)
  thread_block_tile<TSUB_MAC> tile = tiled_partition<TSUB_MAC>(this_thread_block());
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_MAC < 32)

#ifdef  USE_COOPERATIVE_GROUPS
  grid_group grid = this_grid();
#endif//USE_COOPERATIVE_GROUPS

#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
  int smem;
#else///USE_WARP_SHUFFLE_FUNC_MAC
  __shared__  int_real  smem[NTHREADS_MAC];
#endif//USE_WARP_SHUFFLE_FUNC_MAC
  __shared__ jparticle pj_sm[NTHREADS_MAC];
  __shared__      real rjbuf[NTHREADS_MAC * NBUF_MAC];
  __shared__  int      list0[NTHREADS_MAC * NBUF_MAC];
  __shared__  int      list1[NTHREADS_MAC * NBUF_MAC];
  __shared__  int      pjidx[NTHREADS_MAC * NBUF_MAC];  /**< pjidx would be able to be removed by integrating with list; however, reducing SM usage from 24KB to 20KB has no effect @ Ttot = 512 */

  /** head index of remote buffers */
  const int bufHead = (DIV_TSUB_MAC(head) + bidx * NGROUPS_MAC) * NUM_ALLOC_MACBUF;


  /** calculate multipole moment(s) of pseudo j-particles */
  for(int levelIdx = bottomLev; levelIdx >= 0; levelIdx--){
#ifdef  COUNT_INTERACTIONS
    tree_stats local = {ZERO, ZERO, ZERO, ZERO, 0, 0};
#endif//COUNT_INTERACTIONS

    /** load a tree cell and evaluate either the cell is a node or a leaf */
    /** if cidx is not less than tail, then the thread has a null leaf-cell */
    const PHinfo clev = level[levelIdx];

/* #ifndef NDEBUG */
/*     if( gidx == 0 ) */
/*       printf("l%d: levelIdx = %d, bsize = %d, bnum = %d\n", __LINE__, levelIdx, BLOCKSIZE(clev.num, bnum * NGROUPS_MAC), bnum); */
/* #endif//NDEBUG */

#ifdef  USE_COOPERATIVE_GROUPS
    grid.sync();
#else///USE_COOPERATIVE_GROUPS
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
#endif//USE_COOPERATIVE_GROUPS

    /** loop to set a maximum number for # of blocks */
    for(int ii = 0; ii < BLOCKSIZE(clev.num, bnum * NGROUPS_MAC); ii++){
      const      int cidx =  clev.head + DIV_TSUB_MAC(gidx + ii * bnum * NTHREADS_MAC);/**< common in TSUB_MAC (<= 32) threads */
      const treecell root = (cidx < clev.head + clev.num) ? (cell[cidx]) : (null_cell_device);/**< common in TSUB_MAC (<= 32) threads */
/* #ifndef NDEBUG */
/*       if( gidx == 0 ) */
/* 	printf("l%d: cidx = %d\n", __LINE__, cidx); */
/* #endif//NDEBUG */

      /** extend the pseudo particle chain */
      if( root.head != NULL_CELL ){/**< common in TSUB_MAC (<= 32) threads */
#ifdef  COUNT_INTERACTIONS
	local.cellNum++;
#endif//COUNT_INTERACTIONS

	if( !leaf[cidx] ){/**< common in TSUB_MAC (<= 32) threads */
/* #ifndef NDEBUG */
/* 	  if( (lane == 0) && (levelIdx == 15)  && (ii > 3) && (bidx == 305) ) */
/* 	    printf("l%d: %d; %d, %d: cidx = %d\n", __LINE__, ii, bidx, tidx, cidx); */
/* #ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP */
/* 	  /\** synchronization to reduce warp divergence *\/ */
/* #   if  TSUB_MAC == 32 */
/* 	  __syncwarp(); */
/* #else///TSUB_MAC == 32 */
/* 	  tile.sync(); */
/* #endif//TSUB_MAC == 32 */
/* #endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP */
/* #endif//NDEBUG */
	  /** when the tree cell is a node cell, then calculate multipole moment(s) of the cell */
	  jparticle jcom = {ZERO, ZERO, ZERO, ZERO};
#ifdef  WS93_MAC
	  real mjrj2 = FLT_MIN;
#endif//WS93_MAC

	  /** sum up multipole moment(s) of child cells */
	  /** only leaf cells can point multiple tree nodes */
	  const  int jidx  = node[cidx] & IDXMASK;
	  uint more_tmp = more[jidx];
	  int cnum  = (more_tmp >> IDXBITS) + 1;
	  int chead = (more_tmp  & IDXMASK);

	  /** calculate multipole moments of the pseudo particle group */
	  for(int jj = lane; jj < cnum; jj += TSUB_MAC){
	    const int       sidx = chead + jj;
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const real      mass = mj[sidx].mass;
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const real      mass = mj[sidx];
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const jparticle jpos = pj[sidx];

	    /** calculate total mass */
	    jcom.w += mass;
	    /** calculate center-of-mass */
	    jcom.x += mass * jpos.x;
	    jcom.y += mass * jpos.y;
	    jcom.z += mass * jpos.z;
	    /** calculate trace of quadrupole moment */
#ifdef  WS93_MAC
	    mjrj2 += mr2[sidx];
#endif//WS93_MAC
	  }/* for(int jj = lane; jj < cnum; jj += TSUB_MAC){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
	  __syncwarp();
#else///TSUB_MAC == 32
	  tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

	  /** sum up partial sums within TSUB_MAC threads */
	  pj_sm[tidx] = jcom;
	  jparticle jtmp;
#   if  TSUB_MAC >=  2
	  jtmp =  pj_sm[tidx ^  1];	  jcom.x += jtmp.x;	  jcom.y += jtmp.y;	  jcom.z += jtmp.z;	  jcom.w += jtmp.w;	  pj_sm[tidx] = jcom;
#   if  TSUB_MAC >=  4
	  jtmp =  pj_sm[tidx ^  2];	  jcom.x += jtmp.x;	  jcom.y += jtmp.y;	  jcom.z += jtmp.z;	  jcom.w += jtmp.w;	  pj_sm[tidx] = jcom;
#   if  TSUB_MAC >=  8
	  jtmp =  pj_sm[tidx ^  4];	  jcom.x += jtmp.x;	  jcom.y += jtmp.y;	  jcom.z += jtmp.z;	  jcom.w += jtmp.w;	  pj_sm[tidx] = jcom;
#   if  TSUB_MAC >= 16
	  jtmp =  pj_sm[tidx ^  8];	  jcom.x += jtmp.x;	  jcom.y += jtmp.y;	  jcom.z += jtmp.z;	  jcom.w += jtmp.w;	  pj_sm[tidx] = jcom;
#   if  TSUB_MAC == 32
	  jtmp =  pj_sm[tidx ^ 16];	  jcom.x += jtmp.x;	  jcom.y += jtmp.y;	  jcom.z += jtmp.z;	  jcom.w += jtmp.w;	  pj_sm[tidx] = jcom;
#endif//TSUB_MAC == 32
#endif//TSUB_MAC >= 16
#endif//TSUB_MAC >=  8
#endif//TSUB_MAC >=  4
#endif//TSUB_MAC >=  2
	  jcom = pj_sm[head];

#ifdef  WS93_MAC
#   if  defined(USE_WARP_SHUFFLE_FUNC_MAC) && (TSUB_MAC < 32)
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  const uint SHFL_MASK_TSUB_MAC = __activemask();/**< multiple groups of lanes may call the below warp shuffle instructions */
#else///ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  const uint SHFL_MASK_TSUB_MAC = SHFL_MASK_32;
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#endif//defined(USE_WARP_SHUFFLE_FUNC_MAC) && (TSUB_MAC < 32)
	  mjrj2 = TOTAL_SUM_TSUB(mjrj2
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
				 , SHFL_MASK_TSUB_MAC
#else///USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
				 , smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
				 );
	  if( lane == 0 )
	    mr2[jidx] = mjrj2;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
	  __syncwarp();
#else///TSUB_MAC == 32
	  tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#endif//WS93_MAC


	  /** calculate multipole moments of the pseudo particle group */
	  /** calculate center-of-mass */
	  real mtot = jcom.w;
	  real minv = UNITY / (FLT_MIN + mtot);
	  jcom.x *= minv;
	  jcom.y *= minv;
	  jcom.z *= minv;
#ifdef  WS93_MAC
	  /** calculate trace of quadrupole moment */
	  real B2 = mjrj2 - mtot * (jcom.x * jcom.x + jcom.y * jcom.y + jcom.z * jcom.z);
#endif//WS93_MAC

	  /** commit a pseudo particle */
	  if( lane == 0 ){
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const jmass mj_loc = {mtot, ZERO};
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const jmass mj_loc =  mtot;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    mj[jidx] =  mj_loc;
	  }/* if( lane == 0 ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
	  __syncwarp();
#else///TSUB_MAC == 32
	  tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

#ifdef  COUNT_INTERACTIONS
	  local.nodeNum++;
	  local.mjMean +=        mtot;
	  local.mjSdev += mtot * mtot;
#endif//COUNT_INTERACTIONS


	  /** estimate size of particle distribution */
	  /** initialize list of examined tree nodes and related variables */
	  int inum = root.num;
	  int Ntry = 1;
	  if( lane == 0 )	    list0[hbuf] = cidx;
/* #ifndef NDEBUG */
/* 	  /\* if( inum > NI_BMAX_ESTIMATE ) *\/ */
/* 	  /\*   printf("l%d: inum = %d @ (%d, %d)\n", __LINE__, inum, bidx, tidx); *\/ */
/* 	  if( (cidx == 3504557) || (cidx == 3504558) ) */
/* 	    printf("l%d: inum = %d(%d) @ (%d, %d)\n", __LINE__, inum, NI_BMAX_ESTIMATE, bidx, tidx); */
/* #endif//NDEBUG */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
	  __syncwarp();
#else///TSUB_MAC == 32
	  tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

	  /** pick up NI_BMAX_ESTIMATE i-particles in maximum to estimate bmax */
	  while( inum > NI_BMAX_ESTIMATE ){
/* #ifndef NDEBUG */
/* 	    /\* if( lane == 0 ) *\/ */
/* 	    if( (lane == 0) && (levelIdx == 15) && (bidx == 305) && (ii > 3) && ((tidx == 80) || (tidx == 96)) ) */
/* 	      printf("l%d: inum = %d @ (%d, %d)\n", __LINE__, inum, bidx, tidx); */
/* #ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP */
/* 	    /\** synchronization to reduce warp divergence *\/ */
/* #   if  TSUB_MAC == 32 */
/* 	    __syncwarp(); */
/* #else///TSUB_MAC == 32 */
/* 	    tile.sync(); */
/* #endif//TSUB_MAC == 32 */
/* #endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP */
/* #endif//NDEBUG */
	    real rmin = ZERO;
	    int Nloc = 0;
	    int Nbuf = 0;

	    int Niter = BLOCKSIZE(Ntry, TSUB_MAC * NBUF_MAC);
	    for(int iter = 0; iter < Niter; iter++){
	      const int Nsweep = (Ntry > (TSUB_MAC * NBUF_MAC)) ? (TSUB_MAC * NBUF_MAC) : Ntry;
	      const int ibuf_loop = BLOCKSIZE(Nsweep, TSUB_MAC);
	      for(int ibuf = 0; ibuf < ibuf_loop; ibuf++){
		cnum = 0;
		if( (lane + ibuf * TSUB_MAC) < Nsweep ){
		  /** load a tree node corresponding the tree cell */
		  more_tmp = node[list0[hbuf + lane + ibuf * TSUB_MAC]];
		  const int nodenum  = 1 + (more_tmp >> IDXBITS);
		  const int nodehead =      more_tmp  & IDXMASK;

		  /** load all child nodes of the tree cell */
		  more_tmp = more[nodehead];
		  cnum  = 1 + (more_tmp >> IDXBITS);
		  chead =      more_tmp  & IDXMASK;
		  for(int jj = 1; jj < nodenum; jj++)
		    cnum += (1 + (more[nodehead + jj] >> IDXBITS));
		}/* if( (lane + ibuf * TSUB_MAC) < Nsweep ){ */
/* #ifndef NDEBUG */
/* 		/\* if( lane == 0 ) *\/ */
/* 		if( (cidx == 3504557) || (cidx == 3504558) ) */
/* 		  printf("l%d: cnum = %d (%d/%d) @ (%d, %d)\n", __LINE__, cnum, iter, Niter, bidx, tidx); */
/* #endif//NDEBUG */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		/** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
		__syncwarp();
#else///TSUB_MAC == 32
		tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP


#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#   if  TSUB_MAC < 32
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		uint SHFL_MASK_TSUB_MAC = __activemask();/**< multiple groups of lanes may call the below warp shuffle instructions */
#else///ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		const uint SHFL_MASK_TSUB_MAC = SHFL_MASK_32;
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#endif//TSUB_MAC < 32
		smem = PREFIX_SUM_TSUB(cnum, lane, SHFL_MASK_TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		PREFIX_SUM_TSUB(cnum, lane, (int *)smem, tidx);
#endif//USE_WARP_SHUFFLE_FUNC_MAC

/* #ifndef NDEBUG */
/* 		if( (cidx == 3504557) || (cidx == 3504558) ) */
/* 		  printf("l%d: smem = %d @ (%d, %d)\n", __LINE__, smem, bidx, tidx); */
/* 		/\* if( lane == 0 ) *\/ */
/* 		/\*   printf("l%d: cnum = %d @ (%d, %d)\n", __LINE__, cnum, bidx, tidx); *\/ */
/* #endif//NDEBUG */

#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		const int lend = BLOCKSIZE(__SHFL(SHFL_MASK_TSUB_MAC, smem, TSUB_MAC - 1, TSUB_MAC), NBUF_MAC * TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		const int lend = BLOCKSIZE(                           smem[tail].i,                  NBUF_MAC * TSUB_MAC);
#endif//USE_WARP_SHUFFLE_FUNC_MAC

/* #ifndef NDEBUG */
/* 		if( (cidx == 3504557) || (cidx == 3504558) ) */
/* 		  printf("l%d: lend = %d @ (%d, %d)\n", __LINE__, lend, bidx, tidx); */
/* #endif//NDEBUG */

		for(int ll = 0; ll < lend; ll++){
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  const int unum =
		    (smem         <= (TSUB_MAC * NBUF_MAC)) ? cnum :
		    ((smem         >= (cnum + TSUB_MAC * NBUF_MAC)) ? (0) : (cnum + TSUB_MAC * NBUF_MAC - smem));
		  const int shead = hbuf + smem         - cnum;
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  const int unum =
		    (smem[tidx].i <= (TSUB_MAC * NBUF_MAC)) ? cnum :
		    ((smem[tidx].i >= (cnum + TSUB_MAC * NBUF_MAC)) ? (0) : (cnum + TSUB_MAC * NBUF_MAC - smem[tidx].i));
		  const int shead = hbuf + smem[tidx].i - cnum;
#endif//USE_WARP_SHUFFLE_FUNC_MAC

		  for(int jj = 0; jj < unum; jj++){
		    pjidx[shead + jj] = chead;
		    chead++;
		  }/* for(int jj = 0; jj < unum; jj++){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		  /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
		  __syncwarp();
#else///TSUB_MAC == 32
		  tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		  cnum -= unum;
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  const int Ntmp = smem         - (NBUF_MAC * TSUB_MAC);/**< Ntmp is a temporal buffer */
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  const int Ntmp = smem[tidx].i - (NBUF_MAC * TSUB_MAC);/**< Ntmp is a temporal buffer */
#endif//USE_WARP_SHUFFLE_FUNC_MAC

/* #ifndef NDEBUG */
/* 		  if( (cidx == 3504557) || (cidx == 3504558) ) */
/* 		    printf("l%d: cnum = %d (%d/%d) @ (%d, %d)\n", __LINE__, cnum, ll, lend, bidx, tidx); */
/* #endif//NDEBUG */

		  /** pick up candidate tree nodes */
#   if  NBUF_MAC == 4
		  alignedFlt rjmax_loc = { REAL_MIN,  REAL_MIN,  REAL_MIN,  REAL_MIN};
		  alignedInt pjidx_loc = {NULL_NODE, NULL_NODE, NULL_NODE, NULL_NODE};
#endif//NBUF_MAC == 4
#   if  NBUF_MAC == 2
		  alignedFlt rjmax_loc = { REAL_MIN,  REAL_MIN};
		  alignedInt pjidx_loc = {NULL_NODE, NULL_NODE};
#endif//NBUF_MAC == 2

#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_MAC < 32)
		  SHFL_MASK_TSUB_MAC = __activemask();/**< mask may be change from the previous one */
#endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_MAC < 32)
		  int stail = __SHFL(SHFL_MASK_TSUB_MAC, smem, TSUB_MAC - 1, TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  int stail =                            smem[tail].i                 ;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		  stail = (stail < (TSUB_MAC * NBUF_MAC)) ? stail : (TSUB_MAC * NBUF_MAC);
/* #ifndef NDEBUG */
/* 		  if( (cidx == 3504557) || (cidx == 3504558) ) */
/* 		    printf("l%d: stail = %d (%d/%d) @ (%d, %d)\n", __LINE__, stail, ll, lend, bidx, tidx); */
/* #endif//NDEBUG */

#pragma unroll
		  for(int kk = 0; kk < NBUF_MAC; kk++){
		    const int jj = lane + kk * TSUB_MAC;
		    if( jj >= stail )
		      break;

		    const int kidx = pjidx[hbuf + jj];
		    const jparticle jpos = pj[kidx];

		    const real dx = jpos.x - jcom.x;
		    const real dy = jpos.y - jcom.y;
		    const real dz = jpos.z - jcom.z;
		    const real d2 = FLT_MIN + dx * dx + dy * dy + dz * dz;
		    const real dr = d2 * RSQRT(d2);

		    const real rjmax = bmax[kidx] + dr;
		    const real rjmin = -rjmax + (TWO * (UNITY - EPSILON)) * dr;
		    rmin = FMAX(rjmin, rmin);

		    if( rjmax > rmin ){
		      pjidx_loc.ia[kk] = kidx;
		      rjmax_loc.ra[kk] = rjmax;
		    }/* if( rjmax > rmin ){ */
		  }/* for(int kk = 0; kk < NBUF_MAC; kk++){ */
/* #ifndef NDEBUG */
/* 		  if( (cidx == 3504557) || (cidx == 3504558) ) */
/* 		    printf("l%d: rmin = %e (%d/%d) @ (%d, %d)\n", __LINE__, rmin, ll, lend, bidx, tidx); */
/* #endif//NDEBUG */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		  /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
		  __syncwarp();
#else///TSUB_MAC == 32
		  tile.sync();
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
		  SHFL_MASK_TSUB_MAC = __activemask();/**< mask may be change from the previous one */
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

		  /** share rmin within TSUB_MAC threads */
		  rmin = GET_MAX_TSUB(rmin
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
				      , SHFL_MASK_TSUB_MAC
#else///USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
				      , (real *)smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
				      );
/* #ifndef NDEBUG */
/* 		  if( (cidx == 3504557) || (cidx == 3504558) ) */
/* 		    printf("l%d: rmin = %e (%d/%d) @ (%d, %d)\n", __LINE__, rmin, ll, lend, bidx, tidx); */
/* #endif//NDEBUG */


		  /** recheck local buffer (is really rjmax greater than rmin ?) */
#pragma unroll
		  for(int jj = 0; jj < NBUF_MAC; jj++){
		    const int share = ( rjmax_loc.ra[jj] > rmin ) ? 1 : 0;
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		    smem = PREFIX_SUM_TSUB(share, lane, SHFL_MASK_TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		    PREFIX_SUM_TSUB(share, lane, (int *)smem, tidx);
#endif//USE_WARP_SHUFFLE_FUNC_MAC

		    if( share ){
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		      const int dst = hbuf + Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_MAC
		      const int dst = hbuf + Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		      list1[dst] = pjidx_loc.ia[jj];
		      rjbuf[dst] = rjmax_loc.ra[jj];
		    }/* if( share ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		    /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
		    __syncwarp();
#else///TSUB_MAC == 32
		    tile.sync();
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
		    SHFL_MASK_TSUB_MAC = __activemask();/**< mask may be change from the previous one */
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		    Nloc += __SHFL(SHFL_MASK_TSUB_MAC, smem, TSUB_MAC - 1, TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		    Nloc +=                            smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_MAC

		    if( Nloc > ((NBUF_MAC - 1) * TSUB_MAC) ){
		      for(int kk = lane; kk < Nloc; kk += TSUB_MAC){
			more1Buf[bufHead + Nbuf + kk] = list1[hbuf + kk];
			rjmaxBuf[bufHead + Nbuf + kk] = rjbuf[hbuf + kk];
		      }/* for(int kk = lane; kk < Nloc; kk += TSUB_MAC){ */

		      Nbuf += Nloc;
		      Nloc = 0;
		    }/* if( Nloc > ((NBUF_MAC - 1) * TSUB_MAC) ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		    /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
		    __syncwarp();
#else///TSUB_MAC == 32
		    tile.sync();
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
		    SHFL_MASK_TSUB_MAC = __activemask();/**< mask may be change from the previous one */
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		  }/* for(int jj = 0; jj < NBUF_MAC; jj++){ */

#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  smem         = Ntmp;/**< Ntmp is a temporal buffer */
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  smem[tidx].i = Ntmp;/**< Ntmp is a temporal buffer */
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		}/* for(int ll = 0; ll < lend; ll++){ */
	      }/* for(int ibuf = 0; ibuf < NBUF_MAC; ibuf++){ */

	      Ntry -= Nsweep;

	      /** copy data from global memory to shared memory */
	      const int Ncopy = (Ntry < (TSUB_MAC * NBUF_MAC)) ? (Ntry) : (TSUB_MAC * NBUF_MAC);
	      for(int jj = lane; jj < Ncopy; jj += TSUB_MAC)
		list0[hbuf + jj] = more0Buf[(bufHead + (TSUB_MAC * NBUF_MAC) * (iter + 1)) + jj];
	    }/* for(int iter = 0; iter < Niter; iter++){ */

	    if( Nbuf != 0 ){
	      for(int ll = lane; ll < Nloc; ll += TSUB_MAC){
		more1Buf[bufHead + Nbuf + ll] = list1[hbuf + ll];
		rjmaxBuf[bufHead + Nbuf + ll] = rjbuf[hbuf + ll];
	      }/* for(int ll = lane; ll < Nloc; ll += TSUB_MAC){ */

	      for(int ll = lane; ll < TSUB_MAC * NBUF_MAC; ll += TSUB_MAC){
		list1[hbuf + ll] = more1Buf[bufHead + ll];
		rjbuf[hbuf + ll] = rjmaxBuf[bufHead + ll];
	      }/* for(int ll = lane; ll < TSUB_MAC * NBUF_MAC; ll += TSUB_MAC){ */
	    }/* if( Nbuf != 0 ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	    /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
	    __syncwarp();
#else///TSUB_MAC == 32
	    tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

	    Ntry = Nbuf + Nloc;
	    if( (lane == 0) && (Ntry > NUM_ALLOC_MACBUF) )
	      atomicAdd(overflow, 1);
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	    /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
	    __syncwarp();
#else///TSUB_MAC == 32
	    tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP


	    /** list up all child nodes that satisfy rjmax > rmin */
	    inum = 0;
	    Nloc = 0;
	    Nbuf = 0;

	    Niter = BLOCKSIZE(Ntry, NBUF_MAC * TSUB_MAC);
	    for(int iter = 0; iter < Niter; iter++){
	      const int krem = (Ntry < (NBUF_MAC * TSUB_MAC)) ? Ntry : (NBUF_MAC * TSUB_MAC);
	      const int knum = BLOCKSIZE(krem, TSUB_MAC);

	      for(int ki = 0; ki < knum; ki++){
		int cellIdx = lane + ki * TSUB_MAC;
		int  add = 0;
		int iadd = 0;

		/** select distant tree cells */
		if( cellIdx < krem ){
		  /** when the current node must be taken into account */
		  cellIdx += hbuf;
		  if( rjbuf[cellIdx] > rmin ){
		    /** count up total number of contained i-particles */
		    cellIdx = node2cell[list1[cellIdx]];
		    iadd = cell[cellIdx].num;

		    add = 1;
		  }/* if( rjbuf[cellIdx] > rmin ){ */
		}/* if( cellIdx < krem ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		/** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
		__syncwarp();
#else///TSUB_MAC == 32
		tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

#   if  defined(USE_WARP_SHUFFLE_FUNC_MAC) && (TSUB_MAC < 32)
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		uint SHFL_MASK_TSUB_MAC = __activemask();/**< multiple groups of lanes may call the below warp shuffle instructions */
#else///ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		const uint SHFL_MASK_TSUB_MAC = SHFL_MASK_32;
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#endif//defined(USE_WARP_SHUFFLE_FUNC_MAC) && (TSUB_MAC < 32)

		/** remove duplicated tree cells */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		smem = PREFIX_SUM_TSUB(add, lane, SHFL_MASK_TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		PREFIX_SUM_TSUB(add, lane, (int *)smem, tidx);
#endif//USE_WARP_SHUFFLE_FUNC_MAC

		if( add ){
		  /** test uploading... */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  const int smidx = Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  const int smidx = Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		  list0[hbuf + smidx] = cellIdx;

		  /** if detect duplication, upload flag is turned off */
		  if( ((smidx > 0) && (list0[hbuf + smidx - 1] == cellIdx)) || ((Nbuf > 0) && (smidx == 0) && (more0Buf[bufHead + Nbuf - 1] == cellIdx)) ){
		    add  = 0;
		    iadd = 0;
		  }/* if( ((smidx > 0) && (list0[hbuf + smidx - 1] == cellIdx)) || ((Nbuf > 0) && (smidx == 0) && (more0Buf[bufHead + Nbuf - 1] == cellIdx)) ){ */
		}/* if( add ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		/** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
		__syncwarp();
#else///TSUB_MAC == 32
		tile.sync();
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		SHFL_MASK_TSUB_MAC = __activemask();/**< multiple groups of lanes may call the below warp shuffle instructions */
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

		/** save tree cells on the local buffer */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		smem = PREFIX_SUM_TSUB(add, lane, SHFL_MASK_TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		PREFIX_SUM_TSUB(add, lane, (int *)smem, tidx);
#endif//USE_WARP_SHUFFLE_FUNC_MAC

		if( add ){
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  const int smidx = Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  const int smidx = Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		  list0[hbuf + smidx] = cellIdx;
		}/* if( add ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		/** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
		__syncwarp();
#else///TSUB_MAC == 32
		tile.sync();
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		SHFL_MASK_TSUB_MAC = __activemask();/**< multiple groups of lanes may call the below warp shuffle instructions */
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		Nloc += __SHFL(SHFL_MASK_TSUB_MAC, smem, TSUB_MAC - 1, TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		Nloc +=                            smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_MAC

		/** move data to the remote buffer if necessary */
		if( Nloc > ((NBUF_MAC - 1) * TSUB_MAC) ){
		  for(int ll = lane; ll < Nloc; ll += TSUB_MAC)
		    more0Buf[bufHead + Nbuf + ll] = list0[hbuf + ll];
		  Nbuf += Nloc;
		  Nloc  = 0;
		}/* if( Nloc > ((NBUF_MAC - 1) * TSUB_MAC) ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		/** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
		__syncwarp();
#else///TSUB_MAC == 32
		tile.sync();
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		SHFL_MASK_TSUB_MAC = __activemask();/**< multiple groups of lanes may call the below warp shuffle instructions */
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

		/** sum up iadd within TSUB_MAC threads */
		iadd = TOTAL_SUM_TSUB(iadd
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
				      , SHFL_MASK_TSUB_MAC
#else///USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
				      , (int *)smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
				      );

		inum += iadd;
	      }/* for(int ki = 0; ki < knum; ki++){ */

	      Ntry -= krem;

	      /** copy data from remote buffer to local buffer */
	      const int Ncopy = (Ntry < (TSUB_MAC * NBUF_MAC)) ? (Ntry) : (TSUB_MAC * NBUF_MAC);
	      for(int jj = lane; jj < Ncopy; jj += TSUB_MAC){
		rjbuf[hbuf + jj] = rjmaxBuf[bufHead + NBUF_MAC * TSUB_MAC * (iter + 1) + jj];
		list1[hbuf + jj] = more1Buf[bufHead + NBUF_MAC * TSUB_MAC * (iter + 1) + jj];
	      }/* for(int jj = lane; jj < Ncopy; jj += TSUB_MAC){ */
	    }/* for(int iter = 0; iter < Niter1; iter++){ */

	    if( Nbuf != 0 ){
	      for(int ll = lane; ll < Nloc; ll += TSUB_MAC)
		more0Buf[bufHead + Nbuf + ll] = list0[hbuf + ll];

	      for(int ll = lane; ll < NBUF_MAC * TSUB_MAC; ll += TSUB_MAC)
		list0[hbuf + ll] = more0Buf[bufHead + ll];
	    }/* if( Nbuf != 0 ){ */

	    Ntry = Nloc + Nbuf;
	    if( (lane == 0) && (Ntry > NUM_ALLOC_MACBUF) )
	      atomicAdd(overflow, 1);

/* #ifndef NDEBUG */
/* 	    /\* if( lane == 0 ) *\/ */
/* 	    if( (levelIdx == 15) && (ii > 3) && (bidx == 305) && ((head == 80) || (head == 96)) ) */
/* 	      printf("l%d: inum = %d(%d) @ (%d, %d)\n", __LINE__, inum, NI_BMAX_ESTIMATE, bidx, tidx); */
/* #endif//NDEBUG */

#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	    /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
	    __syncwarp();
#else///TSUB_MAC == 32
	    tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  }/* while( inum > NI_BMAX_ESTIMATE ){ */


	  /** check positions of all the pick upped i-particles */
	  real jbmax = ZERO;
	  /** Since NI_BMAX_ESTIMATE <= TSUB_MAC * NBUF_MAC, Ntry1 is less than TSUB_MAC * NBUF_MAC */
	  /** load index of the pick upped i-particles to list1 */
	  const int Niter = BLOCKSIZE(Ntry, TSUB_MAC);
	  int Ncand = 0;
	  for(int iter = 0; iter < Niter; iter++){
	    treecell cand;
	    int pnum = 0;
	    if( lane < Ntry ){
	      cand = cell[list0[hbuf + iter * TSUB_MAC + lane]];
	      pnum = cand.num;
	    }/* if( lane < Ntry ){ */
/* #ifndef NDEBUG */
/* 	  /\* if( lane == 0 ) *\/ */
/* 	  if( (cidx == 3504557) || (cidx == 3504558) ) */
/* 	    printf("l%d: pnum = %d @ (%d, %d)\n", __LINE__, pnum, bidx, tidx); */
/* #endif//NDEBUG */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	    /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
	    __syncwarp();
#else///TSUB_MAC == 32
	    tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#   if  TSUB_MAC < 32
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		uint SHFL_MASK_TSUB_MAC = __activemask();/**< multiple groups of lanes may call the below warp shuffle instructions */
#else///ENABLE_IMPLICIT_SYNC_WITHIN_WARP
		const uint SHFL_MASK_TSUB_MAC = SHFL_MASK_32;
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#endif//TSUB_MAC < 32
	    smem = PREFIX_SUM_TSUB(pnum, lane, SHFL_MASK_TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
	    PREFIX_SUM_TSUB(pnum, lane, (int *)smem, tidx);
#endif//USE_WARP_SHUFFLE_FUNC_MAC

#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
	    for(int jj = 0; jj < pnum; jj++)
	      list1[hbuf + Ncand + smem         - pnum + jj] = cand.head + jj;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	    /** synchronization for consistency */
#   if  TSUB_MAC == 32
	    __syncwarp();
#else///TSUB_MAC == 32
	    tile.sync();
	    SHFL_MASK_TSUB_MAC = __activemask();/**< multiple groups of lanes may call the below warp shuffle instructions */
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	    Ncand += __SHFL(SHFL_MASK_TSUB_MAC, smem, TSUB_MAC - 1, TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
	    for(int jj = 0; jj < pnum; jj++)
	      list1[hbuf + Ncand + smem[tidx].i - pnum + jj] = cand.head + jj;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	    /** synchronization for consistency */
#   if  TSUB_MAC == 32
	    __syncwarp();
#else///TSUB_MAC == 32
	    tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	    Ncand += smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
/* #ifndef NDEBUG */
/* 	    if( lane == 0 ) */
/* 	      printf("l%d: Ncand = %d @ (%d, %d)\n", __LINE__, Ncand, bidx, tidx); */
/* #endif//NDEBUG */

	    Ntry -= TSUB_MAC;
	  }/* for(int iter = 0; iter < Niter; iter++){ */

	  for(int jj = lane; jj < Ncand; jj += TSUB_MAC){
	    const position ipos = pi[list1[hbuf + jj]];

	    const real dx = ipos.x - jcom.x;
	    const real dy = ipos.y - jcom.y;
	    const real dz = ipos.z - jcom.z;
	    const real r2 = FLT_MIN + dx * dx + dy * dy + dz * dz;

	    jbmax = FMAX(jbmax, r2);
	  }/* for(int jj = lane; jj < Ncand; jj += TSUB_MAC){ */
/* #ifndef NDEBUG */
/* 	  /\* if( lane == 0 ) *\/ */
/* 	  if( (lane == 0) && (levelIdx == 15) && (ii > 3) && (bidx == 305) && ((tidx == 80) || (tidx == 96)) ) */
/* 	    printf("l%d: jbmax = %e @ (%d, %d)\n", __LINE__, jbmax, bidx, tidx); */
/* #endif//NDEBUG */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
	  __syncwarp();
#else///TSUB_MAC == 32
	  tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

#   if  defined(USE_WARP_SHUFFLE_FUNC_MAC) && (TSUB_MAC < 32)
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  uint SHFL_MASK_TSUB_MAC = __activemask();/**< multiple groups of lanes may call the below warp shuffle instructions */
#else///ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  const uint SHFL_MASK_TSUB_MAC = SHFL_MASK_32;
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
#endif//defined(USE_WARP_SHUFFLE_FUNC_MAC) && (TSUB_MAC < 32)
	  jbmax = GET_MAX_TSUB(jbmax
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
			       , SHFL_MASK_TSUB_MAC
#else///USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
			       , (real *)smem, tidx, head
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
			       );

	  jbmax *= RSQRT(jbmax);
	  if( lane == 0 )
	    bmax[jidx] = jbmax;
/* #ifndef NDEBUG */
/* 	  /\* if( lane == 0 ) *\/ */
/* 	  if( (cidx == 3504557) || (cidx == 3504558) ) */
/* 	    printf("l%d: jbmax = %e @ (%d, %d)\n", __LINE__, jbmax, bidx, tidx); */
/* #endif//NDEBUG */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  /** synchronization to reduce warp divergence */
#   if  TSUB_MAC == 32
	  __syncwarp();
#else///TSUB_MAC == 32
	  tile.sync();
#endif//TSUB_MAC == 32
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

#ifdef  GADGET_MAC
	  jcom.w  = mac_delta * mtot * jbmax * jbmax;
#ifdef  YMIKI_MAC
	  jcom.w *= jcom.w;
#endif//YMIKI_MAC
#else///GADGET_MAC
#ifdef  WS93_MAC
	  const real bmax_2 = HALF * jbmax;
	  B2 *= RSQRT(B2);
	  B2  = (bmax_2 * bmax_2) + mac_delta * B2;
	  jcom.w  = bmax_2 + B2 * RSQRT(B2);
	  jcom.w *= jcom.w;
#else///WS93_MAC
	  jcom.w  = jbmax * jbmax;
#endif//WS93_MAC
#endif//GADGET_MAC

	  if( lane == 0 )
	    pj[jidx] = jcom;
/* #ifndef NDEBUG */
/* 	  if( (lane == 0) && (levelIdx == 15) && (ii > 3) && (bidx == 305) ) */
/* 	    printf("l%d: %d; %d, %d\n", __LINE__, ii, bidx, tidx); */
/* #endif//NDEBUG */

#ifdef  COUNT_INTERACTIONS
	  local.nodeNum++;
	  local.r2Mean +=          jcom.w;
	  local.r2Sdev += jcom.w * jcom.w;
#endif//COUNT_INTERACTIONS
	}/* if( !leaf[cidx] ){ */
      }/* if( root.head != NULL_CELL ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
      __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    }/* for(int ii = 0; ii < BLOCKSIZE(clev.num, bnum * NGROUPS_MAC); ii++){ */


#ifdef  COUNT_INTERACTIONS
    if( lane == 0 ){
      atomicAdd(&(stats[levelIdx].nodeNum), local.nodeNum);
      atomicAdd(&(stats[levelIdx].cellNum), local.cellNum);
      atomicAdd(&(stats[levelIdx].mjMean), local. mjMean);
      atomicAdd(&(stats[levelIdx].mjSdev), local. mjSdev);
      atomicAdd(&(stats[levelIdx].r2Mean), local. r2Mean);
      atomicAdd(&(stats[levelIdx].r2Sdev), local. r2Sdev);
    }/* if( lane == 0 ){ */
#endif//COUNT_INTERACTIONS
  }/* for(int levelIdx = bottomLev; levelIdx >= 0; levelIdx--){ */
}


#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
/**
 * @fn detectOuterParticle
 *
 * @brief Detect particles locate exterior of the local domain.
 */
#   if  NTHREADS_SCAN_INC != NTHREADS_OUTFLOW
#undef	NTHREADS_SCAN_INC
#define NTHREADS_SCAN_INC    NTHREADS_OUTFLOW
#include "../util/scan_inc.cu"
#endif//NTHREADS_SCAN_INC != NTHREADS_OUTFLOW
__device__ __forceinline__ bool detectOuterParticle(const jparticle jpos, const float rad, const float3 bmin, const float3 bmax)
{
  float sep = 0.5f * FLT_MAX;
  float pjx = CAST_R2F(jpos.x);  pjx = fminf((pjx - rad) - bmin.x, bmax.x - (pjx + rad));
  float pjy = CAST_R2F(jpos.y);  pjy = fminf((pjy - rad) - bmin.y, bmax.y - (pjy + rad));
  float pjz = CAST_R2F(jpos.z);  pjz = fminf((pjz - rad) - bmin.z, bmax.z - (pjz + rad));
  sep = fminf(sep, pjx);
  pjy = fminf(pjy, pjz);
  sep = fminf(pjy, sep);
  return ((sep > 0.0f) ? false : true);
}


/**
 * @fn copyData_g2g
 *
 * @brief Move data from global memory to global memory
 */
__device__ __forceinline__ void copyData_g2g
  (uint * RESTRICT ubuf,
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
   float * RESTRICT fbuf,
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
   size_t srcHead, size_t dstHead, int Ncopy, const int Ndisp, const int gidx, const int Nthread, const int tidx, const int bidx, const int bnum, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  /** configure the settings */
  const int Nfirst = Ndisp % Nthread;
  const int  ldIdx = (gidx + Nfirst) % Nthread;  /**< ldIdx is Nfirst, Nfirst + 1, ..., Nthread - 1, 0, 1, ..., Nfirst - 1 for gidx of 0, 1, 2, ..., Nthread - 1 */
  const int grpIdx = (ldIdx < Nfirst) ? 0 : 1;

  srcHead += Ndisp - Nfirst;/**< hereafter, srcHead is warpSize elements aligned */

  /** fraction processing at loading from the head of source array */
  uint  utmp = ubuf[srcHead + ldIdx];
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
  float ftmp = fbuf[srcHead + ldIdx];
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
  srcHead += Nthread;

  /** sequential load and store from source to destination on the global memory */
  const int Niter = BLOCKSIZE(Ncopy, Nthread);
  for(int iter = 0; iter < Niter; iter++){
    const int Nmove = (Ncopy > Nthread) ? (Nthread) : (Ncopy);

    /** load from the source array on the global memory */
    /** load from temp (fraction processing) as initialization */
    uint  uloc = utmp;
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
    float floc = ftmp;
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)

    /** load from global memory, store to shared memory or temp (fraction processing) */
    utmp = ubuf[srcHead + ldIdx];
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
    ftmp = fbuf[srcHead + ldIdx];
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
    if( !grpIdx ){
      uloc = utmp;
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
      floc = ftmp;
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
    }/* if( !grpIdx ){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

    /** store to the destination array on the global memory */
    ubuf[dstHead + gidx] = uloc;
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
    fbuf[dstHead + gidx] = floc;
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)

    Ncopy   -= Nmove;
    srcHead += Nmove;
    dstHead += Nmove;

    globalSync(tidx, bidx, bnum, gsync0, gsync1);
  }/* for(int iter = 0; iter < Niter; iter++){ */
}


/**
 * @fn checkOutflow_kernel
 *
 * @brief Detect particles locate exterior of the local domain.
 */
__global__ void
#ifndef USE_OCCUPANCY_CALCULATOR
__launch_bounds__(NTHREADS_OUTFLOW, NBLOCKS_PER_SM_OUTFLOW)
#endif//USE_OCCUPANCY_CALCULATOR
checkOutflow_kernel
  (const int topLev, READ_ONLY PHinfo * level, READ_ONLY uint * RESTRICT ptag,
   READ_ONLY uint * RESTRICT more, jparticle * RESTRICT pj, real * RESTRICT bmax,
   const float3 boxmin, const float3 boxmax,
   uint * RESTRICT ubuf,
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
#ifdef  TIME_BASED_MODIFICATION
   const float linv,
#endif//TIME_BASED_MODIFICATION
   float * RESTRICT fbuf,
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
   int * RESTRICT overflow, const size_t bufSize,
   int * RESTRICT gmem, int * RESTRICT gsync0, int * RESTRICT gsync1)
{
  /** identify thread properties */
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int bnum =   GRIDDIM_X1D;

  const int lane = tidx & (TSUB_OUTFLOW - 1);  /**< for local summation */
  const int lane32 = tidx & (warpSize - 1);  /**< for globalPrefixSum */

/* #   if  !defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_OUTFLOW < 32) */
/*   thread_block_tile<TSUB_OUTFLOW> tile = tiled_partition<TSUB_OUTFLOW>(this_thread_block()); */
/* #endif//!defined(ENABLE_IMPLICIT_SYNC_WITHIN_WARP) && (TSUB_OUTFLOW < 32) */

  __shared__ int smem[NTHREADS_OUTFLOW];


  /** initialization */
  int rem = 0;
  {
    const PHinfo clev = level[topLev];

    for(int ii = 0; ii < BLOCKSIZE(clev.num, bnum * NGROUPS_OUTFLOW); ii++){
      /** pick up a tree cell */
      const  int cidx = DIV_TSUB_OUTFLOW(gidx + ii * bnum * NTHREADS_OUTFLOW);/**< common in TSUB_OUTFLOW (<= 32) threads */
      const uint node = (cidx < clev.num) ? (ptag[clev.head + cidx]) : (NULL_NODE);

      /** pick up tree nodes and check the position */
      const uint jnum  = (node >> IDXBITS) + 1;
      const uint jhead = (node	& IDXMASK);
      int num = 0;
      int list = 0;

      if( node != NULL_NODE )
	for(uint jj = lane; jj < jnum; jj += TSUB_OUTFLOW){
	  if( more[jhead + jj] != (jhead + jj) )
	    if( detectOuterParticle(pj[jhead + jj], CAST_R2F(bmax[jhead + jj]), boxmin, boxmax) ){
	      num++;
	      list |= 1 << (DIV_TSUB_OUTFLOW(jj));
	    }/* if( detectOuterParticle(pj[jhead + jj], CAST_R2F(bmax[jhead + jj])) ){ */
	}/* for(uint jj = lane; jj < jnum; jj += TSUB_OUTFLOW){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
      __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

      /** exclusive scan within a grid */
      int scanNum;
      int headIdx = PREFIX_SUM_GRID(num, smem, lane32, tidx, &scanNum, gmem, bidx, bnum, gsync0, gsync1);
      headIdx -= num;/**< this must be an exclusive scan */

      /** store the picked tree nodes on the global memory */
#if 1
      the below if block contains a bug: every threads have individual num; therefore, lane should not be used in the for block;
      careful revision is required;
#endif
      if( num != 0 ){
	int kk = 0;
	for(uint jj = lane; jj < jnum; jj += TSUB_OUTFLOW){
	  if( (list >> DIV_TSUB_OUTFLOW(jj)) & 1 ){
	    ubuf[rem + headIdx + kk] =          more[jhead + jj];
#ifdef  USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES
	    fbuf[rem + headIdx + kk] = CAST_R2F(  pj[jhead + jj].w);
#endif//USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES
#ifdef  TIME_BASED_MODIFICATION
	    fbuf[rem + headIdx + kk] = ldexpf(linv, topLev);
#endif//TIME_BASED_MODIFICATION
	    kk++;
	  }/* if( (list >> DIV_TSUB_OUTFLOW(jj)) & 1 ){ */
	  if( kk == num )
	    break;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	  __syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	}/* for(int jj = 0; jj < DIV_TSUB_OUTFLOW(NLEAF); jj++){ */
      }/* if( num != 0 ){ */






      rem += scanNum;
    }/* for(int ii = 0; ii < BLOCKSIZE(clev.num, bnum * NTHREADS_OUTFLOW); ii++){ */
  }


  /** detect particles locate outside the preset domain */
  size_t hbuf = 0;
  size_t tbuf = rem;

  while( true ){
    /** if the queue becomes empty, then exit the while loop */
    if( rem == 0 )
      break;

    /** pick up a tree node per TSUB_OUTFLOW threads */
    const size_t nidx = hbuf + DIV_TSUB_OUTFLOW(gidx);
    const uint   root = (nidx < tbuf) ? (ubuf[nidx]) : (NULL_NODE);
    const int checked = (rem < (bnum * NGROUPS_OUTFLOW)) ? (rem) : (bnum * NGROUPS_OUTFLOW);
    rem  -= checked;
    hbuf += checked;
    if( (rem == 0) && (hbuf > (checked + bnum * NGROUPS_OUTFLOW)) )
      hbuf = tbuf = 0;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

    /** pick up tree nodes and check the position */
    const uint jnum  = (root >> IDXBITS) + 1;
    const uint jhead = (root  & IDXMASK);
    int num = 0;
    int list = 0;
    if( root != NULL_NODE )
      for(uint jj = lane; jj < jnum; jj += TSUB_OUTFLOW){
	if( more[jhead + jj] != (jhead + jj) )
	  if( detectOuterParticle(pj[jhead + jj], CAST_R2F(bmax[jhead + jj]), boxmin, boxmax) ){
	    num++;
	    list |= 1 << (DIV_TSUB_OUTFLOW(jj));
	  }/* if( detectOuterParticle(pj[jhead + jj], CAST_R2F(bmax[jhead + jj])) ){ */
      }/* for(uint jj = lane; jj < jnum; jj += TSUB_OUTFLOW){ */
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

    /** exclusive scan within a grid */
    int scanNum;
    int headIdx = PREFIX_SUM_GRID(num, smem, lane32, tidx, &scanNum, gmem, bidx, bnum, gsync0, gsync1);
    headIdx -= num;/**< this must be an exclusive scan */

    /** edit MAC of the external tree nodes and commit their child nodes as next candidates */
#if 1
    the below if block contains a bug: every threads have individual num; therefore, lane should not be used in the for block;
    careful revision is required;
#endif
    if( num != 0 ){
#ifdef  TIME_BASED_MODIFICATION
      const float base_inv = fbuf[nidx];
#endif//TIME_BASED_MODIFICATION
#ifdef  USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES
      const float fmac = fbuf[nidx];
#endif//USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES
      int kk = 0;
      for(uint jj = lane; jj < jnum; jj += TSUB_OUTFLOW){
	if( (list >> DIV_TSUB_OUTFLOW(jj)) & 1 ){
	  ubuf[tbuf + (size_t)(headIdx + kk)] =          more[jhead + jj];
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
#ifdef  TIME_BASED_MODIFICATION
	  fbuf[tbuf + (size_t)(headIdx + kk)] = 2.0f * base_inv;
	  /** we want calculate 2^(5/4)^n */
#ifdef  GADGET_MAC
	  pj[jhead + jj].w *= exp2f(5.0f * fmaxf(ceilf(log2(CAST_R2F(bmax[jhead + jj]) * base_inv)), 1.0f));
#else///GADGET_MAC
	  pj[jhead + jj].w *= exp2f(2.5f * fmaxf(ceilf(log2(CAST_R2F(bmax[jhead + jj]) * base_inv)), 1.0f));
#endif//GADGET_MAC
#endif//TIME_BASED_MODIFICATION
#ifdef  USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES
	  fbuf[tbuf + (size_t)(headIdx + kk)] = CAST_R2F(  pj[jhead + jj].w);
	  pj[jhead + jj].w = CAST_F2R(fmac);
#endif//USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES
#else///defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
	  /** in the worst case, 8 processes can share the same location */
	  /** then, m_J approaches to m_J / 8 */
#ifdef  GADGET_MAC
	  /** diJ^4 is \propto \propto m_J */
	  pj[jhead + jj].w *= 8.0f;
	  /* pj[jhead + jj].w *= 32.0f;/\* 32 = 2^5 *\/ */
#else///GADGET_MAC
	  /** in limit of b_J is 0, diJ^2 is \propto m_J^1/2 */
	  /** pj[jhead + jj].w *= 2.0f; */
	  pj[jhead + jj].w *= 5.656854249f;/* 5.656854249 = 2^2.5 */
#endif//GADGET_MAC
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
	  kk++;
	}/* if( (list >> DIV_TSUB_OUTFLOW(jj)) & 1 ){ */
	if( kk == num )
	  break;
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
	__syncwarp();
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
      }/* for(uint jj = lane; jj < jnum; jj += TSUB_OUTFLOW){ */
    }/* if( num != 0 ){ */

    tbuf += scanNum;
    rem  += scanNum;
    if( (gidx == 0) && (tbuf > bufSize) )
      atomicAdd(overflow, 1);
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP

    globalSync(tidx, bidx, bnum, gsync0, gsync1);

    if( bufSize - tbuf < DIV_TSUB_OUTFLOW(bnum * NTHREADS_OUTFLOW) )
      copyData_g2g(ubuf,
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
		   fbuf,
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
		   0, 0, tbuf - hbuf, hbuf, gidx, bnum * NTHREADS_OUTFLOW, tidx, bidx, bnum, gsync0, gsync1);
#ifndef ENABLE_IMPLICIT_SYNC_WITHIN_WARP
    __syncwarp();/**< __syncwarp() to remove warp divergence */
#endif//ENABLE_IMPLICIT_SYNC_WITHIN_WARP
  }/* while( true ){ */
}
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)


#ifdef  COUNT_INTERACTIONS
/**
 * @fn initStatistics_kernel
 *
 * @brief Initialize counter for Nj and Nbuf.
 */
__global__ void initStatistics_kernel(tree_stats *stats)
{
  const tree_stats zero_stats = {ZERO, ZERO, ZERO, ZERO, 0, 0};

  for(int ll = 0; ll < MAXIMUM_PHKEY_LEVEL; ll++)
    stats[ll] = zero_stats;
}
#endif//COUNT_INTERACTIONS


/**
 * @fn calcMultipole_dev
 *
 * @brief Calculate multipole moment(s) of tree cells based on the width-first search.
 */
extern "C"
void calcMultipole_dev
(const int bottomLev, const soaTreeCell cell,
 const int piNum, const iparticle pi, const int pjNum, const soaTreeNode node,
 const soaMakeTreeBuf buf, deviceProp devProp
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
 , const real eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifndef SERIALIZED_EXECUTION
#ifdef  CARE_EXTERNAL_PARTICLES
 , domainLocation *location
#endif//CARE_EXTERNAL_PARTICLES
 , double *tmac
#endif//SERIALIZED_EXECUTION
#ifdef  COUNT_INTERACTIONS
 , const soaTreeCell cell_hst, tree_stats * RESTRICT stats_hst
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
 )
{
  __NOTE__("%s\n", "start");


#ifdef  COUNT_INTERACTIONS
  tree_stats *stats_dev;
  mycudaMalloc((void **)&stats_dev, MAXIMUM_PHKEY_LEVEL * sizeof(tree_stats));
  initStatistics_kernel<<<1, 1>>>(stats_dev);
#endif//COUNT_INTERACTIONS

#ifdef  EXEC_BENCHMARK
  initStopwatch();
#else///EXEC_BENCHMARK
#ifndef SERIALIZED_EXECUTION
  static struct timespec start;
  checkCudaErrors(cudaDeviceSynchronize());
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
#endif//SERIALIZED_EXECUTION
#endif//EXEC_BENCHMARK


  /** initialize position and mass of pseudo j-particles */
  int Nrem = BLOCKSIZE(pjNum, NTHREADS_INIT_BODY);
  __NOTE__("initTreeBody_kernel for Nrem = %d (pjNum = %d)\n", Nrem, pjNum);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    initTreeBody_kernel<<<Nrem, NTHREADS_INIT_BODY>>>(node.jpos, node.mj, node.bmax
#ifdef  WS93_MAC
						      , node.mr2
#endif//WS93_MAC
						      );
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_INIT_BODY;
      initTreeBody_kernel<<<Nblck, NTHREADS_INIT_BODY>>>(&node.jpos[hidx], &node.mj[hidx], &node.bmax[hidx]
#ifdef  WS93_MAC
							 , &node.mr2[hidx]
#endif//WS93_MAC
							 );

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("initTreeBody_kernel");

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->initTreeBody_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER


  /** set position and mass of j-particles */
  Nrem = BLOCKSIZE(piNum, NTHREADS_COPY_BODY);
  __NOTE__("copyRealBody_kernel for Nrem = %d (piNum = %d)\n", Nrem, piNum);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    copyRealBody_kernel<<<Nrem, NTHREADS_COPY_BODY>>>(piNum, node.jtag,
#ifdef  BLOCK_TIME_STEP
						      pi.jpos,
#else///BLOCK_TIME_STEP
						      pi. pos,
#endif//BLOCK_TIME_STEP
						      node.jpos, node.mj
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
						      , eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  WS93_MAC
						      , node.mr2
#endif//WS93_MAC
						      );
  else{
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

      int Nsub = Nblck * NTHREADS_COPY_BODY;
      copyRealBody_kernel<<<Nblck, NTHREADS_COPY_BODY>>>(piNum, &node.jtag[hidx],
#ifdef  BLOCK_TIME_STEP
							 &pi.jpos[hidx],
#else///BLOCK_TIME_STEP
							 &pi. pos[hidx],
#endif//BLOCK_TIME_STEP
							 node.jpos, node.mj
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
							 , eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  WS93_MAC
							 , node.mr2
#endif//WS93_MAC
							 );

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("copyRealBody_kernel");

#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->copyRealBody_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER


  /** set pseudo j-particles */
#ifdef  USE_OCCUPANCY_CALCULATOR
  const int NBLOCKS_PER_SM_MAC = buf.numBlocksPerSM_mac;
#endif//USE_OCCUPANCY_CALCULATOR
  __NOTE__("launch calcMultipole_kernel: %d blocks per SM, %d SMs; bottomLev = %d\n", NBLOCKS_PER_SM_MAC, devProp.numSM, bottomLev);
#ifdef  USE_COOPERATIVE_GROUPS
  const int bottomLev_tmp = bottomLev - 1;
#ifdef  COUNT_INTERACTIONS
  tree_stats *stats_tmp = &stats_dev[0];
#endif//COUNT_INTERACTIONS
  void *kernelArgs[] = {
    (void *)&bottomLev_tmp,
    (void *)&cell.level,
    (void *)&cell.cell,
    (void *)&cell.leaf,
#ifdef  BLOCK_TIME_STEP
    (void *)&pi.jpos,
#else///BLOCK_TIME_STEP
    (void *)&pi. pos,
#endif//BLOCK_TIME_STEP
    (void *)&cell.ptag,
    (void *)&node.more,
    (void *)&node.node2cell,
    (void *)&node.jpos,
    (void *)&node.mj,
    (void *)&node.bmax,
    (void *)&buf.more0,
    (void *)&buf.more1,
    (void *)&buf.rjmax,
    (void *)&buf.fail
#ifdef  WS93_MAC
    , (void *)&node.mr2
#endif//WS93_MAC
#ifdef  COUNT_INTERACTIONS
    , (void *)&stats_tmp
#endif//COUNT_INTERACTIONS
  };
  checkCudaErrors(cudaLaunchCooperativeKernel((void *)calcMultipole_kernel, devProp.numSM * NBLOCKS_PER_SM_MAC, NTHREADS_MAC, kernelArgs));
#else///USE_COOPERATIVE_GROUPS
  calcMultipole_kernel<<<devProp.numSM * NBLOCKS_PER_SM_MAC, NTHREADS_MAC>>>
    (bottomLev - 1, cell.level, cell.cell, cell.leaf,
#ifdef  BLOCK_TIME_STEP
     pi.jpos,
#else///BLOCK_TIME_STEP
     pi. pos,
#endif//BLOCK_TIME_STEP
     cell.ptag, node.more, node.node2cell, node.jpos, node.mj, node.bmax,
     buf.more0, buf.more1, buf.rjmax, buf.fail, buf.gsync0, buf.gsync1
#ifdef  WS93_MAC
     , node.mr2
#endif//WS93_MAC
#ifdef  COUNT_INTERACTIONS
     , &stats_dev[0]
#endif//COUNT_INTERACTIONS
     );
#endif//USE_COOPERATIVE_GROUPS

  getLastCudaError("calcMultipole_kernel");

  int fail_hst;
  checkCudaErrors(cudaMemcpy(&fail_hst, buf.fail, sizeof(int), cudaMemcpyDeviceToHost));
  if( fail_hst != 0 ){
    __KILL__(stderr, "ERROR: buffer (%d elements per %d threads group) overflow at least %d times.\nPLEASE re-simulate after increasing NUM_ALLOC_MACBUF defined in src/tree/make.h.\n",
	     NUM_ALLOC_MACBUF, TSUB_MAC, fail_hst);
  }/* if( fail_hst != 0 ){ */

#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  float3 boxmin = location->boxmin;
  float3 boxmax = location->boxmax;

#ifdef  MODIFY_INTERNAL_DOMAIN
#ifdef  TIME_BASED_MODIFICATION
  const float skin = location->elapsed * location->dtinv * location->eps * location->eta;
#else///TIME_BASED_MODIFICATION
  const float skin =                     location->step  * location->eps * location->eta;
#endif//TIME_BASED_MODIFICATION
  boxmin.x += skin;  boxmin.y += skin;  boxmin.z += skin;
  boxmax.x -= skin;  boxmax.y -= skin;  boxmax.z -= skin;
#endif//MODIFY_INTERNAL_DOMAIN

  int topLev = NUM_PHKEY_LEVEL;
#ifdef  TIME_BASED_MODIFICATION
  if( location->elapsed > 0.5f * location->dtmin ){
    topLev = -(int)ceilf(log2f(location->elapsed * location->dtinv * location->eta * location->eps * location->linv)) - 1;/**< commit parent tree cells */
    if( topLev < 0 )
      topLev = 0;
    if( topLev >= bottomLev )
      topLev = NUM_PHKEY_LEVEL;
  }/* if( location->elapsed > 0.5f * location->dtmin ){ */
#else///TIME_BASED_MODIFICATION
  if( location->step > 0.5f ){
    topLev = -(int)ceilf(log2f(location->step * location->dL_L)) - 1;/**< commit parent tree cells */
    if( topLev < 0 )
      topLev = 0;
    if( topLev >= bottomLev )
      topLev = NUM_PHKEY_LEVEL;
  }/* if( location->step > 0.5f ){ */
#endif//TIME_BASED_MODIFICATION

  __NOTE__("checkOutflow_kernel\n");
#ifdef  USE_OCCUPANCY_CALCULATOR
  const int NBLOCKS_PER_SM_OUTFLOW = buf.numBlocksPerSM_outflow;
#endif//USE_OCCUPANCY_CALCULATOR
  checkOutflow_kernel<<<devProp.numSM * NBLOCKS_PER_SM_OUTFLOW, NTHREADS_OUTFLOW>>>
    (topLev, cell.level, cell.ptag, node.more, node.jpos, node.bmax, boxmin, boxmax,
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
     (uint *)buf.ibuf_external,
#ifdef  TIME_BASED_MODIFICATION
     location->linv,
#endif//TIME_BASED_MODIFICATION
     (float *)buf.rbuf_external,
#else///defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
     buf.ubuf_external,
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
     buf.fail, buf.Nbuf_external, buf.gmem_external, buf.gsync0_external, buf.gsync1_external);

  getLastCudaError("checkOutflow_kernel");
  checkCudaErrors(cudaMemcpy(&fail_hst, buf.fail, sizeof(int), cudaMemcpyDeviceToHost));
  if( fail_hst != 0 ){
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
    __KILL__(stderr, "ERROR: buffer (%zu elements in total) overflow at least %d times.\nPLEASE re-simulate after increasing NUM_ALLOC_MACBUF(%d) defined in src/tree/make.h.\n",
	     buf.Nbuf_external, fail_hst, NUM_ALLOC_MACBUF);
#else///defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
    __KILL__(stderr, "ERROR: buffer (%zu elements in total) overflow at least %d times.\nPLEASE re-simulate after decreasing NUM_BODY_MAX(%d) or GLOBAL_MEMORY_SYSBUF(%zu) defined in src/misc/structure.h or TREE_SAFETY_VAL(%f) defined in src/tree/make.h, or EXTEND_NUM_TREE_NODE(%f) defined in src/tree/let.h.\n",
      buf.Nbuf_external, fail_hst, NUM_BODY_MAX, (size_t)GLOBAL_MEMORY_SYSBUF, TREE_SAFETY_VAL, EXTEND_NUM_TREE_NODE);
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
  }/* if( fail_hst != 0 ){ */

#ifndef TIME_BASED_MODIFICATION
  location->step += 1.0f;
#endif//TIME_BASED_MODIFICATION
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)


#ifdef  EXEC_BENCHMARK
  double time = 0.0;
  stopStopwatch(&time);
  elapsed->calcMultipole_dev = time;
  /* *tmac += time; */
#else///EXEC_BENCHMARK
#ifndef SERIALIZED_EXECUTION
  static struct timespec finish;
  checkCudaErrors(cudaDeviceSynchronize());
  clock_gettime(CLOCK_MONOTONIC_RAW, &finish);
  *tmac += calcElapsedTimeInSec(start, finish);
#endif//SERIALIZED_EXECUTION
#endif//EXEC_BENCHMARK


#ifdef  COUNT_INTERACTIONS
  checkCudaErrors(cudaMemcpy(cell_hst.level, cell.level, sizeof(PHinfo) * MAXIMUM_PHKEY_LEVEL, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(stats_hst, stats_dev, MAXIMUM_PHKEY_LEVEL * sizeof(tree_stats), cudaMemcpyDeviceToHost));
  mycudaFree(stats_dev);
#endif//COUNT_INTERACTIONS


  __NOTE__("%s\n", "end");
}


/**
 * @fn allocTreeNode_dev
 *
 * @brief Allocate arrays to store properties of tree nodes (tree cells and N-body particles).
 */
extern "C"
muse allocTreeNode_dev
(soaTreeNode *dev, uint **more_dev, jparticle **pj_dev, jmass **mj_dev, real **bmax_dev, int **n2c_dev, int **gsync0, int **gsync1,
#ifdef  WS93_MAC
 real **mr2_dev,
#endif//WS93_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
 soaTreeNode *hst, uint **more_hst, jparticle **pj_hst, jmass **mj_hst,
#endif//!defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
 int **gmem_make_tree, int **gsync0_make_tree, int **gsync1_make_tree, int **gsync2_make_tree, int **gsync3_make_tree,
 int **gmem_link_tree, int **gsync0_link_tree, int **gsync1_link_tree,
#ifdef  GADGET_MAC
 real **mac_dev,
#endif//GADGET_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
 int **gmem_external, int **gsync0_external, int **gsync1_external, float **diameter_dev, float **diameter_hst, domainLocation *location, const float eps, const float eta,
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
 int **more0Buf, int **more1Buf, real **rjmaxBuf, int **fail_dev, soaMakeTreeBuf *buf, const deviceProp devProp)
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

#ifdef  USE_OCCUPANCY_CALCULATOR
/* #ifdef  WS93_MAC */
/*   __NOTE__("WS93_MAC is defined\n"); */
/* #endif//WS93_MAC */
/* #ifdef  GADGET_MAC */
/*   __NOTE__("GADGET_MAC is defined\n"); */
/* #endif//GADGET_MAC */
/* #ifdef  SERIALIZED_EXECUTION */
/*   __NOTE__("SERIALIZED_EXECUTION is defined\n"); */
/* #endif//SERIALIZED_EXECUTION */
/* #ifdef  MPI_VIA_HOST */
/*   __NOTE__("MPI_VIA_HOST is defined\n"); */
/* #endif//MPI_VIA_HOST */
/* #ifdef  CARE_EXTERNAL_PARTICLES */
/*   __NOTE__("CARE_EXTERNAL_PARTICLES is defined\n"); */
/* #endif//CARE_EXTERNAL_PARTICLES */
  /* __NOTE__("address of dev is %p\n", dev); */
  /* __NOTE__("address of more_dev is %p\n", more_dev); */
  /* __NOTE__("address of pj_dev is %p\n", pj_dev); */
  /* __NOTE__("address of mj_dev is %p\n", mj_dev); */
  /* __NOTE__("address of bmax_dev is %p\n", bmax_dev); */
  /* __NOTE__("address of n2c_dev is %p\n", n2c_dev); */
  /* __NOTE__("address of gsync0 is %p\n", gsync0); */
  /* __NOTE__("address of gsync1 is %p\n", gsync1); */
  /* __NOTE__("address of gmem_make_tree is %p\n", gmem_make_tree); */
  /* __NOTE__("address of gsync0_make_tree is %p\n", gsync0_make_tree); */
  /* __NOTE__("address of gsync1_make_tree is %p\n", gsync1_make_tree); */
  /* __NOTE__("address of gsync2_make_tree is %p\n", gsync2_make_tree); */
  /* __NOTE__("address of gsync3_make_tree is %p\n", gsync3_make_tree); */
  /* __NOTE__("address of gmem_link_tree is %p\n", gmem_link_tree); */
  /* __NOTE__("address of gsync0_link_tree is %p\n", gsync0_link_tree); */
  /* __NOTE__("address of gsync1_link_tree is %p\n", gsync1_link_tree); */
  /* __NOTE__("address of mac_dev is %p\n", mac_dev); */
  /* __NOTE__("address of more0Buf is %p\n", more0Buf); */
  /* __NOTE__("address of more1Buf is %p\n", more1Buf); */
  /* __NOTE__("address of rjmaxBuf is %p\n", rjmaxBuf); */
  /* __NOTE__("address of fail_dev is %p\n", fail_dev); */
  /* __NOTE__("address of buf is %p\n", buf); */
  checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&buf->numBlocksPerSM_mac, calcMultipole_kernel, NTHREADS_MAC, 0));
  const int NBLOCKS_PER_SM_MAC = buf->numBlocksPerSM_mac;
  __NOTE__("NBLOCKS_PER_SM_MAC = %d\n", NBLOCKS_PER_SM_MAC);
  checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&buf->numBlocksPerSM_make_tree, makeTree_kernel, NTHREADS_MAKE_TREE, 0));
  const int NBLOCKS_PER_SM_MAKE_TREE = buf->numBlocksPerSM_make_tree;
  __NOTE__("NBLOCKS_PER_SM_MAKE_TREE = %d\n", NBLOCKS_PER_SM_MAKE_TREE);
  checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&buf->numBlocksPerSM_link_tree, linkTree_kernel, NTHREADS_LINK_TREE, 0));
  const int NBLOCKS_PER_SM_LINK_TREE = buf->numBlocksPerSM_link_tree;
  __NOTE__("NBLOCKS_PER_SM_LINK_TREE = %d\n", NBLOCKS_PER_SM_LINK_TREE);
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&buf->numBlocksPerSM_outflow, checkOutflow_kernel, NTHREADS_OUTFLOW, 0));
  const int NBLOCKS_PER_SM_OUTFLOW = buf->numBlocksPerSM_outflow;
  __NOTE__("NBLOCKS_PER_SM_OUTFLOW = %d\n", NBLOCKS_PER_SM_OUTFLOW);
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
#endif//USE_OCCUPANCY_CALCULATOR


  const size_t num = (size_t)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE);
  __NOTE__("num = %zu\n", num);
  mycudaMalloc    ((void **) more_dev, num * sizeof(uint     ));  alloc.device += num * sizeof(uint     );  dev->more      = * more_dev;
  mycudaMalloc    ((void **)   pj_dev, num * sizeof(jparticle));  alloc.device += num * sizeof(jparticle);  dev->jpos      = *   pj_dev;
  mycudaMalloc    ((void **)   mj_dev, num * sizeof(jmass    ));  alloc.device += num * sizeof(jmass    );  dev->mj        = *   mj_dev;
  mycudaMalloc    ((void **)  n2c_dev, num * sizeof(int      ));  alloc.device += num * sizeof(int      );  dev->node2cell = *  n2c_dev;

#   if  !defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
  __NOTE__("%zu MiB for more_hst\n", num * sizeof(uint) >> 20);
  mycudaMallocHost((void **) more_hst, num * sizeof(uint     ));  alloc.host += num * sizeof(uint     );  hst->more      = * more_hst;
  __NOTE__("%zu MiB for pj_hst\n", num * sizeof(jparticle) >> 20);
  mycudaMallocHost((void **)   pj_hst, num * sizeof(jparticle));  alloc.host += num * sizeof(jparticle);  hst->jpos      = *   pj_hst;
  __NOTE__("%zu MiB for mj_hst\n", num * sizeof(jmass) >> 20);
  mycudaMallocHost((void **)   mj_hst, num * sizeof(jmass    ));  alloc.host += num * sizeof(jmass    );  hst->mj        = *   mj_hst;
#endif//!defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)

#ifdef  GADGET_MAC
  mycudaMalloc((void **)mac_dev, num * sizeof(real));  alloc.device += num * sizeof(real);  dev->mac = *mac_dev;
#endif//GADGET_MAC

#ifdef  WS93_MAC
  mycudaMalloc    ((void **) mr2_dev, num * sizeof(real));  alloc.device += num * sizeof(real);  dev->mr2  = * mr2_dev;
#endif//WS93_MAC
  mycudaMalloc    ((void **)bmax_dev, num * sizeof(real));  alloc.device += num * sizeof(real);  dev->bmax = *bmax_dev;

  mycudaMalloc    ((void **)more0Buf, devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(int ));
  alloc.device +=                     devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(int ) ;
  mycudaMalloc    ((void **)more1Buf, devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(int ));
  alloc.device +=                     devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(int ) ;
  mycudaMalloc    ((void **)rjmaxBuf, devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(real));
  alloc.device +=                     devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(real) ;

#ifndef USE_COOPERATIVE_GROUPS
  mycudaMalloc((void **)gsync0, devProp.numSM * NBLOCKS_PER_SM_MAC * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_MAC * sizeof(int);
  mycudaMalloc((void **)gsync1, devProp.numSM * NBLOCKS_PER_SM_MAC * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_MAC * sizeof(int);
  initGsync_kernel<<<1, devProp.numSM * NBLOCKS_PER_SM_MAC>>>(devProp.numSM * NBLOCKS_PER_SM_MAC, *gsync0, *gsync1);
  getLastCudaError("initGsync_kernel");
#endif//USE_COOPERATIVE_GROUPS

  mycudaMalloc((void **)fail_dev, 1 * sizeof(int));
  alloc.device +=                 1 * sizeof(int);
  const int fail_hst = 0;
  checkCudaErrors(cudaMemcpy(*fail_dev, &fail_hst, sizeof(int), cudaMemcpyHostToDevice));

  mycudaMalloc((void **)  gmem_make_tree, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int);
  mycudaMalloc((void **)gsync0_make_tree, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int);
  mycudaMalloc((void **)gsync1_make_tree, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int);
  mycudaMalloc((void **)gsync2_make_tree, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int);
  mycudaMalloc((void **)gsync3_make_tree, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int);
  mycudaMalloc((void **)  gmem_link_tree, devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int);
  mycudaMalloc((void **)gsync0_link_tree, devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int);
  mycudaMalloc((void **)gsync1_link_tree, devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int);
  initGsync_kernel<<<1, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE>>>(devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE, *gsync0_make_tree, *gsync1_make_tree);
  getLastCudaError("initGsync_kernel");
  initGsync_kernel<<<1, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE>>>(devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE, *gsync2_make_tree, *gsync3_make_tree);
  getLastCudaError("initGsync_kernel");
  initGsync_kernel<<<1, devProp.numSM * NBLOCKS_PER_SM_LINK_TREE>>>(devProp.numSM * NBLOCKS_PER_SM_LINK_TREE, *gsync0_link_tree, *gsync1_link_tree);
  getLastCudaError("initGsync_kernel");

#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  mycudaMalloc((void **)  gmem_external, devProp.numSM * NBLOCKS_PER_SM_OUTFLOW * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_OUTFLOW * sizeof(int);
  mycudaMalloc((void **)gsync0_external, devProp.numSM * NBLOCKS_PER_SM_OUTFLOW * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_OUTFLOW * sizeof(int);
  mycudaMalloc((void **)gsync1_external, devProp.numSM * NBLOCKS_PER_SM_OUTFLOW * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_OUTFLOW * sizeof(int);
  initGsync_kernel<<<1, devProp.numSM * NBLOCKS_PER_SM_OUTFLOW>>>(devProp.numSM * NBLOCKS_PER_SM_OUTFLOW, *gsync0_external, *gsync1_external);
  getLastCudaError("initGsync_kernel");
  mycudaMalloc    ((void **)diameter_dev, sizeof(float));  alloc.device += sizeof(float);
  mycudaMallocHost((void **)diameter_hst, sizeof(float));  alloc.host   += sizeof(float);
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)

  buf->Nbuf   = devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF;

  buf->more0  = *more0Buf;
  buf->more1  = *more1Buf;
  buf->rjmax  = *rjmaxBuf;
  buf->fail   = *fail_dev;
#ifndef USE_COOPERATIVE_GROUPS
  buf->gsync0 = *gsync0;
  buf->gsync1 = *gsync1;
#endif//USE_COOPERATIVE_GROUPS

  buf->  gmem_make_tree = *  gmem_make_tree;
  buf->  gmem_link_tree = *  gmem_link_tree;
  buf->gsync0_make_tree = *gsync0_make_tree;
  buf->gsync1_make_tree = *gsync1_make_tree;
  buf->gsync2_make_tree = *gsync2_make_tree;
  buf->gsync3_make_tree = *gsync3_make_tree;
  buf->gsync0_link_tree = *gsync0_link_tree;
  buf->gsync1_link_tree = *gsync1_link_tree;

#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
  buf->  Nbuf_external = devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF;
  buf->  ibuf_external = *more0Buf;/**< cast required */
  buf->  rbuf_external = *rjmaxBuf;/**< cast required */
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
  buf->  gmem_external = *  gmem_external;
  buf->gsync0_external = *gsync0_external;
  buf->gsync1_external = *gsync1_external;
  location->eps = eps;
  location->eta = eta;
#ifdef  TIME_BASED_MODIFICATION
  location->elapsed = 0.0f;
  location->dtmin   = 0.5f * FLT_MAX;
  location->dtinv   = 1.0f / location->dtmin;
#else///TIME_BASED_MODIFICATION
  location->step = 0.0f;
#endif//TIME_BASED_MODIFICATION
  location->diameter_dev = *diameter_dev;
  location->diameter_hst = *diameter_hst;
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)


  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freeTreeNode_dev
 *
 * @brief Deallocate arrays to store properties of tree nodes (tree cells and N-body particles).
 */
extern "C"
void  freeTreeNode_dev
(uint  *more_dev, jparticle  *pj_dev, jmass  *mj_dev, real  *bmax_dev, int  *n2c_dev, int  *gsync0, int  *gsync1,
#ifdef  WS93_MAC
 real  *mr2_dev,
#endif//WS93_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
 uint  *more_hst, jparticle  *pj_hst, jmass  *mj_hst,
#endif//!defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
 int  *gmem_make_tree, int  *gsync0_make_tree, int  *gsync1_make_tree, int  *gsync2_make_tree, int  *gsync3_make_tree,
 int  *gmem_link_tree, int  *gsync0_link_tree, int  *gsync1_link_tree,
#ifdef  GADGET_MAC
 real  *mac_dev,
#endif//GADGET_MAC
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
 int  *gmem_external, int  *gsync0_external, int  *gsync1_external, float  *diameter_dev, float  *diameter_hst,
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
 int  *more0Buf, int  *more1Buf, real  *rjmaxBuf, int  *fail_dev)
{
  __NOTE__("%s\n", "start");


  mycudaFree    ( more_dev);
  mycudaFree    (   pj_dev);
  mycudaFree    (   mj_dev);
  mycudaFree    (  n2c_dev);

#   if  !defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)
  mycudaFreeHost( more_hst);
  mycudaFreeHost(   pj_hst);
  mycudaFreeHost(   mj_hst);
#endif//!defined(SERIALIZED_EXECUTION) && defined(MPI_VIA_HOST)

#ifdef  GADGET_MAC
  mycudaFree(mac_dev);
#endif//GADGET_MAC

#ifdef  WS93_MAC
  mycudaFree( mr2_dev);
#endif//WS93_MAC
  mycudaFree(bmax_dev);

  mycudaFree(more0Buf);
  mycudaFree(more1Buf);
  mycudaFree(rjmaxBuf);
#ifndef USE_COOPERATIVE_GROUPS
  mycudaFree(gsync0);
  mycudaFree(gsync1);
#endif//USE_COOPERATIVE_GROUPS

  mycudaFree(fail_dev);

  mycudaFree(gmem_make_tree);  mycudaFree(gsync0_make_tree);  mycudaFree(gsync1_make_tree);  mycudaFree(gsync2_make_tree);  mycudaFree(gsync3_make_tree);
  mycudaFree(gmem_link_tree);
  mycudaFree(gsync0_link_tree);  mycudaFree(gsync1_link_tree);

#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  mycudaFree(gmem_external);
  mycudaFree(gsync0_external);
  mycudaFree(gsync1_external);
  mycudaFree(diameter_dev);  mycudaFreeHost(diameter_hst);
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)


  __NOTE__("%s\n", "end");
}


/**
 * @fn setGlobalConstants_make_dev_cu
 *
 * @brief Set global constants for make_dev.cu and initialize kernel functions.
 */
extern "C"
void setGlobalConstants_make_dev_cu
(
#ifdef  GADGET_MAC
 const real accErr, const real newton
#else///GADGET_MAC
#ifdef  WS93_MAC
 const real accErr
#else///WS93_MAC
 void
#endif//WS93_MAC
#endif//GADGET_MAC
 )
{
  __NOTE__("%s\n", "start");


#ifdef  GADGET_MAC
  const real      host_val = newton / accErr;
#else///GADGET_MAC
#ifdef  WS93_MAC
  const real      host_val = SQRTRATIO(THREE, accErr);
#endif//WS93_MAC
#endif//GADGET_MAC
  const jparticle host_pj  = {ZERO, ZERO, -UNITY, ZERO};
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
  const jmass     host_mj  = {ZERO, ZERO};
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
  const jmass     host_mj  =  ZERO;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING

#   if  CUDART_VERSION >= 5000
#   if  defined(GADGET_MAC) || defined(WS93_MAC)
  cudaMemcpyToSymbol( mac_delta        , &host_val , sizeof(real)     , 0, cudaMemcpyHostToDevice);
#endif//defined(GADGET_MAC) || defined(WS93_MAC)
  cudaMemcpyToSymbol( null_cell_device , &null_cell, sizeof(treecell) , 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol( zero_pj          , &host_pj  , sizeof(jparticle), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol( zero_mj          , &host_mj  , sizeof(jmass)    , 0, cudaMemcpyHostToDevice);
#else///CUDART_VERSION >= 5000
#   if  defined(GADGET_MAC) || defined(WS93_MAC)
  cudaMemcpyToSymbol("mac_delta"       , &host_val , sizeof(real)     , 0, cudaMemcpyHostToDevice);
#endif//defined(GADGET_MAC) || defined(WS93_MAC)
  cudaMemcpyToSymbol("null_cell_device", &null_cell, sizeof(treecell) , 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("zero_pj"         , &host_pj  , sizeof(jparticle), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("zero_mj"         , &host_mj  , sizeof(jmass)    , 0, cudaMemcpyHostToDevice);
#endif//CUDART_VERSION >= 5000


#   if  GPUVER < 70
  checkCudaErrors(cudaFuncSetCacheConfig(calcMultipole_kernel, cudaFuncCachePreferShared));
#else///GPUVER < 70
  /* checkCudaErrors(cudaFuncSetAttribute(calcMultipole_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_SM)); */
  checkCudaErrors(cudaFuncSetAttribute(calcMultipole_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_SM_64K));

  /* remove shared memory if __global__ function does not use */
  checkCudaErrors(cudaFuncSetAttribute(initPHhierarchy_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(initTreeCellOffset_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(initTreeNode_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(initTreeLink_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(makeTree_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_SM_32K));
  checkCudaErrors(cudaFuncSetAttribute(linkTree_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_SM_08K));
  checkCudaErrors(cudaFuncSetAttribute(trimTree_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(initTreeBody_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(copyRealBody_kernel, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
#ifdef  GADGET_MAC
  checkCudaErrors(cudaFuncSetAttribute(enforceBarnesHutMAC, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
  checkCudaErrors(cudaFuncSetAttribute(recoverGADGET_MAC, cudaFuncAttributePreferredSharedMemoryCarveout, CARVEOUT_MAX_L1));
#endif//GADGET_MAC

#endif//GPUVER < 70


#ifndef USE_OCCUPANCY_CALCULATOR
  /** error checking before running the kernel */
  struct cudaFuncAttributes funcAttr;

  /** check calcMultipole_kernel() */
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, calcMultipole_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAC ){
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in calcMultipole_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_MAC);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_MAC = %d.\n", GPUGEN, GPUVER, NTHREADS_MAC);
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
    fprintf(stderr, "note: warp shuffle instruction is  enabled\n");
#else///USE_WARP_SHUFFLE_FUNC_MAC
    fprintf(stderr, "note: warp shuffle instruction is disabled\n");
#endif//USE_WARP_SHUFFLE_FUNC_MAC
    fflush (stderr);
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAC ){ */
  int regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS_MAC);
  if( regLimit > (MAX_REGISTERS_PER_SM / NTHREADS_MAC) )
    regLimit = (MAX_REGISTERS_PER_SM / NTHREADS_MAC);
  int memLimit = SMEM_SIZE_SM_PREF / funcAttr.sharedSizeBytes;
  int Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / NTHREADS_MAC) )    Nblck = MAX_THREADS_PER_SM / NTHREADS_MAC;
  if( Nblck >   MAX_BLOCKS_PER_SM                       )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / NTHREADS_MAC) )    Nblck = ((MAX_WARPS_PER_SM * 32) / NTHREADS_MAC);
  if( Nblck != NBLOCKS_PER_SM_MAC ){
    __KILL__(stderr, "ERROR: # of blocks per SM for calcMultipole_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_MAC defined in src/tree/make_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_MAC, funcAttr.numRegs, REGISTERS_PER_THREAD_MAC, GPUGEN, GPUVER);
  }/* if( Nblck != NBLOCKS_PER_SM ){ */

  /** check makeTree_kernel() */
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, makeTree_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAKE_TREE ){
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in makeTree_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_MAKE_TREE);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_MAKE_TREE = %d.\n", GPUGEN, GPUVER, NTHREADS_MAKE_TREE);
    fflush (stderr);
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAKE_TREE ){ */
  regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS_MAKE_TREE);
  if( regLimit > (MAX_REGISTERS_PER_SM / NTHREADS_MAKE_TREE) )
    regLimit = (MAX_REGISTERS_PER_SM / NTHREADS_MAKE_TREE);
  memLimit = SMEM_SIZE_L1_PREF / funcAttr.sharedSizeBytes;
  Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / NTHREADS_MAKE_TREE) )    Nblck = MAX_THREADS_PER_SM / NTHREADS_MAKE_TREE;
  if( Nblck >   MAX_BLOCKS_PER_SM                             )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / NTHREADS_MAKE_TREE) )    Nblck = ((MAX_WARPS_PER_SM * 32) / NTHREADS_MAKE_TREE);
  if( Nblck != NBLOCKS_PER_SM_MAKE_TREE ){
    __KILL__(stderr, "ERROR: # of blocks per SM for makeTree_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_MAKE_TREE defined in src/tree/make_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_MAKE_TREE, funcAttr.numRegs, REGISTERS_PER_THREAD_MAKE_TREE, GPUGEN, GPUVER);
  }/* if( Nblck != NBLOCKS_PER_SM ){ */

  /** check linkTree_kernel() */
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, linkTree_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_LINK_TREE ){
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in linkTree_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_LINK_TREE);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_LINK_TREE = %d.\n", GPUGEN, GPUVER, NTHREADS_LINK_TREE);
    fflush (stderr);
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_LINK_TREE ){ */
  regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS_LINK_TREE);
  if( regLimit > (MAX_REGISTERS_PER_SM / NTHREADS_LINK_TREE) )
    regLimit = (MAX_REGISTERS_PER_SM / NTHREADS_LINK_TREE);
  memLimit = SMEM_SIZE_L1_PREF / funcAttr.sharedSizeBytes;
  Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / NTHREADS_LINK_TREE) )    Nblck = MAX_THREADS_PER_SM / NTHREADS_LINK_TREE;
  if( Nblck >   MAX_BLOCKS_PER_SM                             )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / NTHREADS_LINK_TREE) )    Nblck = ((MAX_WARPS_PER_SM * 32) / NTHREADS_LINK_TREE);
  if( Nblck != NBLOCKS_PER_SM_LINK_TREE ){
    __KILL__(stderr, "ERROR: # of blocks per SM for linkTree_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_LINK_TREE defined in src/tree/make_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_LINK_TREE, funcAttr.numRegs, REGISTERS_PER_THREAD_LINK_TREE, GPUGEN, GPUVER);
  }/* if( Nblck != NBLOCKS_PER_SM ){ */

#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  /** check checkOutflow_kernel() */
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, checkOutflow_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_OUTFLOW ){
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in checkOutflow_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_OUTFLOW);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_OUTFLOW = %d.\n", GPUGEN, GPUVER, NTHREADS_OUTFLOW);
    fflush (stderr);
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_OUTFLOW ){ */
  regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS_OUTFLOW);
  if( regLimit > (MAX_REGISTERS_PER_SM / NTHREADS_OUTFLOW) )
    regLimit = (MAX_REGISTERS_PER_SM / NTHREADS_OUTFLOW);
  memLimit = SMEM_SIZE_L1_PREF / funcAttr.sharedSizeBytes;
  Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / NTHREADS_OUTFLOW) )    Nblck = MAX_THREADS_PER_SM / NTHREADS_OUTFLOW;
  if( Nblck >   MAX_BLOCKS_PER_SM                           )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / NTHREADS_OUTFLOW) )    Nblck = ((MAX_WARPS_PER_SM * 32) / NTHREADS_OUTFLOW);
  if( Nblck != NBLOCKS_PER_SM_OUTFLOW ){
    __KILL__(stderr, "ERROR: # of blocks per SM for checkOutflow_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_OUTFLOW defined in src/tree/make_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_OUTFLOW, funcAttr.numRegs, REGISTERS_PER_THREAD_OUTFLOW, GPUGEN, GPUVER);
  }/* if( Nblck != NBLOCKS_PER_SM ){ */
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
#endif//USE_OCCUPANCY_CALCULATOR


  __NOTE__("%s\n", "end");
}


#ifndef NDEBUG
/**
 * @fn setGlobalConstants_make_dev_cu
 *
 * @brief Set global constants for make_dev.cu and initialize kernel functions.
 */
extern "C"
void printPHkey_location(const int piNum, PHint * RESTRICT peano, const iparticle pi)
{
  __NOTE__("%s\n", "start");

  PHint *key_hst;
  mycudaMallocHost((void **)&key_hst, sizeof(PHint) * piNum);
  checkCudaErrors(cudaMemcpy(key_hst, peano, sizeof(PHint) * piNum, cudaMemcpyDeviceToHost));

  position *pos_hst;
  mycudaMallocHost((void **)&pos_hst, sizeof(position) * piNum);
  checkCudaErrors(cudaMemcpy(pos_hst, pi.pos, sizeof(position) * piNum, cudaMemcpyDeviceToHost));

  ulong *idx_hst;
  mycudaMallocHost((void **)&idx_hst, sizeof(ulong) * piNum);
  checkCudaErrors(cudaMemcpy(idx_hst, pi.idx, sizeof(ulong) * piNum, cudaMemcpyDeviceToHost));

  PHint old = key_hst[0];
  for(int ii = 1; ii < piNum; ii++){
    PHint key = key_hst[ii];
    if( key == old )
      fprintf(stderr, "%d\t%zu\t%e\t%e\t%e\t%e\t%zu\n", ii, key, pos_hst[ii].x, pos_hst[ii].y, pos_hst[ii].z, pos_hst[ii].m, idx_hst[ii]);
    old = key;
  }/* for(int ii = 1; ii < piNum; ii++){ */

  mycudaFreeHost(key_hst);
  mycudaFreeHost(pos_hst);
  mycudaFreeHost(idx_hst);

  __NOTE__("%s\n", "end");
}
#endif//NDEBUG
