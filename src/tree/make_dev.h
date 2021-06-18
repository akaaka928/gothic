/**
 * @file make_dev.h
 *
 * @brief Header file for constructing octree structure on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/12/08 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef MAKE_DEV_H
#define MAKE_DEV_H


#include <stdbool.h>

#include "macro.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "../sort/peano.h"
#include "../tree/make.h"

#ifndef SERIALIZED_EXECUTION
#include "../sort/peano_dev.h"
#endif//SERIALIZED_EXECUTION


#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && !defined(TIME_BASED_MODIFICATION)
/* #define USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES */
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES) && !defined(TIME_BASED_MODIFICATION)


#ifndef HUNT_NODE_PARAMETER
/* the below macro is enabled in the default option; switched off in the parameter survey mode to use -D and -U from Makefile */
#define USE_WARP_SHUFFLE_FUNC_MAC
#define USE_WARP_REDUCE_FUNC_MAC

/** for better performance on V100 GPUs */
#   if  (GPUVER == 70) || (GPUVER == 75)
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#undef  USE_WARP_SHUFFLE_FUNC_MAC
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#endif//(GPUVER == 70) || (GPUVER == 75)

/** for better performance on A100 GPUs */
#   if  (GPUVER >= 80)
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#undef  USE_WARP_SHUFFLE_FUNC_MAC
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#ifdef  USE_WARP_REDUCE_FUNC_MAC
#undef  USE_WARP_REDUCE_FUNC_MAC
#endif//USE_WARP_REDUCE_FUNC_MAC
#endif//(GPUVER >= 80)

#endif//HUNT_NODE_PARAMETER


#   if  defined(USE_WARP_SHUFFLE_FUNC_MAC) && (GPUVER < 30)
#undef          USE_WARP_SHUFFLE_FUNC_MAC
#endif//defined(USE_WARP_SHUFFLE_FUNC_MAC) && (GPUVER < 30)

#   if  defined(USE_WARP_REDUCE_FUNC_MAC) && !defined(ENABLE_WARP_REDUCE_FUNC)
#undef          USE_WARP_REDUCE_FUNC_MAC
#endif//defined(USE_WARP_REDUCE_FUNC_MAC) && !defined(ENABLE_WARP_REDUCE_FUNC)


/**
 * @def NTHREADS_MAC
 *
 * @brief number of threads per block for calcMultipole_kernel
 */
#ifndef NTHREADS_MAC
#   if  GPUVER >= 80
#   if  GPUGEN <= 60
/** Pascal mode */
#define NTHREADS_MAC (512)
#else///GPUGEN <= 60
/** Ampere mode */
#define NTHREADS_MAC (128)
#endif//GPUGEN <= 60
#else///GPUVER >= 80
#   if  (GPUVER == 52) || (GPUVER == 60) || (GPUVER == 61)
#define NTHREADS_MAC (256)
#else///(GPUVER == 52) || (GPUVER == 60) || (GPUVER == 61)
#define NTHREADS_MAC (128)
#endif//(GPUVER == 52) || (GPUVER == 60) || (GPUVER == 61)
#endif//GPUVER >= 80
#endif//NTHREADS_MAC

/** NTHREADS_MAC must be equal or smaller than 512 due to the capacity of shared memory */
#   if  NTHREADS_MAC > 512
#undef  NTHREADS_MAC
#define NTHREADS_MAC  (512)
#endif//NTHREADS_MAC > 512

#ifdef  DIV_NTHREADS_MAC
#undef  DIV_NTHREADS_MAC
#endif//DIV_NTHREADS_MAC
#   if  NTHREADS_MAC == 512
#define DIV_NTHREADS_MAC(a) ((a) >> 9)
#endif//NTHREADS_MAC == 512
#   if  NTHREADS_MAC == 256
#define DIV_NTHREADS_MAC(a) ((a) >> 8)
#endif//NTHREADS_MAC == 256
#   if  NTHREADS_MAC == 128
#define DIV_NTHREADS_MAC(a) ((a) >> 7)
#endif//NTHREADS_MAC == 128


/**
 * @def TSUB_MAC
 *
 * @brief number of threads that share a common tree node for calcMultipole_kernel
 */
#ifndef TSUB_MAC
#   if  (GPUVER == 60) || (GPUVER == 61)
#define TSUB_MAC (16)
#else///(GPUVER == 60) || (GPUVER == 61)
#define TSUB_MAC (32)
#endif//(GPUVER == 60) || (GPUVER == 61)
#endif//TSUB_MAC

/** TSUB_MAC must be equal or smaller than NTHREADS_MAC */
#   if  TSUB_MAC > NTHREADS_MAC
#undef  TSUB_MAC
#define TSUB_MAC   NTHREADS_MAC
#endif//TSUB_MAC > NTHREADS_MAC

#   if  NTHREADS_MAC >= 128
/** TSUB_MAC must be equal or greater than 4 */
#   if  TSUB_MAC < 4
#undef  TSUB_MAC
#define TSUB_MAC  (4)
#endif//TSUB_MAC < 4
#   if  NTHREADS_MAC >= 256
/** TSUB_MAC must be equal or greater than 8 */
#   if  TSUB_MAC < 8
#undef  TSUB_MAC
#define TSUB_MAC  (8)
#endif//TSUB_MAC < 8
#   if  NTHREADS_MAC == 512
/** TSUB_MAC must be equal to 32 */
#   if  TSUB_MAC != 32
#undef  TSUB_MAC
#define TSUB_MAC   (32)
#endif//TSUB_MAC != 32
#endif//NTHREADS_MAC == 512
#endif//NTHREADS_MAC >= 256
#endif//NTHREADS_MAC >= 128


#define NBUF_MAC (((NJ_BMAX_ESTIMATE) + (TSUB_MAC) - 1) / TSUB_MAC)

/** maximum number of NBUF_MAC is 4 to use float4 in union */
#   if  NBUF_MAC > 4
#undef  NBUF_MAC
#define NBUF_MAC (4)
#undef  TSUB_MAC
#define TSUB_MAC ((NJ_BMAX_ESTIMATE + (NBUF_MAC) - 1) >> 2)
#endif//NBUF_MAC > 4

#ifdef  DIV_TSUB_MAC
#undef  DIV_TSUB_MAC
#endif//DIV_TSUB_MAC
#   if  TSUB_MAC == 32
#define DIV_TSUB_MAC(a) ((a) >> 5)
#define SHFL_MASK_TSUB_MAC SHFL_MASK_32
#endif//TSUB_MAC == 32
#   if  TSUB_MAC == 16
#define DIV_TSUB_MAC(a) ((a) >> 4)
#endif//TSUB_MAC == 16
#   if  TSUB_MAC ==  8
#define DIV_TSUB_MAC(a) ((a) >> 3)
#endif//TSUB_MAC ==  8
#   if  TSUB_MAC ==  4
#define DIV_TSUB_MAC(a) ((a) >> 2)
#endif//TSUB_MAC ==  4
#   if  TSUB_MAC ==  2
#define DIV_TSUB_MAC(a) ((a) >> 1)
#endif//TSUB_MAC ==  2
#   if  TSUB_MAC ==  1
#define DIV_TSUB_MAC(a) (a)
#endif//TSUB_MAC ==  1


#define NGROUPS_MAC (DIV_TSUB_MAC(NTHREADS_MAC))


#   if  NI_BMAX_ESTIMATE > (TSUB_MAC * NBUF_MAC)
#undef  NI_BMAX_ESTIMATE
#define NI_BMAX_ESTIMATE   (TSUB_MAC * NBUF_MAC)
#endif//NI_BMAX_ESTIMATE > (TSUB_MAC * NBUF_MAC)


#ifndef HUNT_MAKE_PARAMETER
/* the below macro is disabled in the default option for better performance; switched off in the parameter survey mode to use -D from Makefile */
/* #define USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE */
// #define USE_WARP_REDUCE_FUNC_MAKE_TREE_STRUCTURE

/** for better performance on P100, V100, and A100 GPUs */
#   if  !defined(USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE) && (GPUVER >= 60)
#define          USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
#endif//!defined(USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE) && (GPUVER >= 60)
#endif//HUNT_MAKE_PARAMETER


#   if  defined(USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE) && (GPUVER < 30)
#undef          USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE
#endif//defined(USE_WARP_SHUFFLE_FUNC_MAKE_TREE_STRUCTURE) && (GPUVER < 30)

#   if  defined(USE_WARP_REDUCE_FUNC_MAKE_TREE_STRUCTURE) && !defined(ENABLE_WARP_REDUCE_FUNC)
#undef          USE_WARP_REDUCE_FUNC_MAKE_TREE_STRUCTURE
#endif//defined(USE_WARP_REDUCE_FUNC_MAKE_TREE_STRUCTURE) && !defined(ENABLE_WARP_REDUCE_FUNC)


/**
 * @def NTHREADS_MAKE_TREE
 *
 * @brief number of threads per block for makeTree_kernel
 */
#ifndef NTHREADS_MAKE_TREE
#   if  (GPUVER >= 60)
#define NTHREADS_MAKE_TREE (512)
#else///(GPUVER >= 60)
#define NTHREADS_MAKE_TREE (128)
#endif//(GPUVER >= 60)
#endif//NTHREADS_MAKE_TREE


/**
 * @def NTHREADS_LINK_TREE
 *
 * @brief number of threads per block for linkTree_kernel
 */
#ifndef NTHREADS_LINK_TREE
#   if  (GPUVER >= 80)
#   if  GPUGEN <= 60
// Pascal mode
#define NTHREADS_LINK_TREE (512)
#else///GPUGEN <= 60
// Ampere mode
#define NTHREADS_LINK_TREE (256)
#endif//GPUGEN <= 60
#else///(GPUVER >= 80)
#   if  (GPUVER >= 30)
#define NTHREADS_LINK_TREE (256)
#else///(GPUVER >= 30)
#define NTHREADS_LINK_TREE (128)
#endif//(GPUVER >= 30)
#endif//(GPUVER >= 80)
#endif//NTHREADS_LINK_TREE


/**
 * @def NTHREADS_TRIM_TREE
 *
 * @brief number of threads per block for trimTree_kernel
 */
#ifndef NTHREADS_TRIM_TREE
#   if  (GPUVER >= 80)
#define NTHREADS_TRIM_TREE (128)
#else///(GPUVER >= 80)
#   if  (GPUVER >= 70)
#define NTHREADS_TRIM_TREE (256)
#else///(GPUVER >= 70)
#define NTHREADS_TRIM_TREE (128)
#endif//(GPUVER >= 70)
#endif//(GPUVER >= 80)
#endif//NTHREADS_TRIM_TREE


/** NTHREADS_MAKE_TREE must be equal or smaller than 512 due to the capacity of shared memory */
#   if  NTHREADS_MAKE_TREE > 512
#undef  NTHREADS_MAKE_TREE
#define NTHREADS_MAKE_TREE  (512)
#endif//NTHREADS_MAKE_TREE > 512

#ifdef  DIV_NTHREADS_MAKE_TREE
#undef  DIV_NTHREADS_MAKE_TREE
#endif//DIV_NTHREADS_MAKE_TREE
#   if  NTHREADS_MAKE_TREE == 512
#define DIV_NTHREADS_MAKE_TREE(a) ((a) >> 9)
#endif//NTHREADS_MAKE_TREE == 512
#   if  NTHREADS_MAKE_TREE == 256
#define DIV_NTHREADS_MAKE_TREE(a) ((a) >> 8)
#endif//NTHREADS_MAKE_TREE == 256
#   if  NTHREADS_MAKE_TREE == 128
#define DIV_NTHREADS_MAKE_TREE(a) ((a) >> 7)
#endif//NTHREADS_MAKE_TREE == 128

#ifdef  DIV_NTHREADS_LINK_TREE
#undef  DIV_NTHREADS_LINK_TREE
#endif//DIV_NTHREADS_LINK_TREE
#   if  NTHREADS_LINK_TREE == 1024
#define DIV_NTHREADS_LINK_TREE(a) ((a) >> 10)
#endif//NTHREADS_LINK_TREE == 1024
#   if  NTHREADS_LINK_TREE ==  512
#define DIV_NTHREADS_LINK_TREE(a) ((a) >>  9)
#endif//NTHREADS_LINK_TREE ==  512
#   if  NTHREADS_LINK_TREE ==  256
#define DIV_NTHREADS_LINK_TREE(a) ((a) >>  8)
#endif//NTHREADS_LINK_TREE ==  256
#   if  NTHREADS_LINK_TREE ==  128
#define DIV_NTHREADS_LINK_TREE(a) ((a) >>  7)
#endif//NTHREADS_LINK_TREE ==  128


/**
 * @def TSUB_MAKE_TREE
 *
 * @brief number of threads that share a common tree cell for makeTree_kernel
 */
#ifndef TSUB_MAKE_TREE
#define TSUB_MAKE_TREE   CELL_UNIT
#endif//TSUB_MAKE_TREE
#   if  TSUB_MAKE_TREE > CELL_UNIT
#undef  TSUB_MAKE_TREE
#define TSUB_MAKE_TREE   CELL_UNIT
#endif//TSUB_MAKE_TREE > CELL_UNIT

#ifdef  DIV_TSUB_MAKE_TREE
#undef  DIV_TSUB_MAKE_TREE
#endif//DIV_TSUB_MAKE_TREE
#   if  TSUB_MAKE_TREE == 32
#define DIV_TSUB_MAKE_TREE(a) ((a) >> 5)
#endif//TSUB_MAKE_TREE == 32
#   if  TSUB_MAKE_TREE == 16
#define DIV_TSUB_MAKE_TREE(a) ((a) >> 4)
#endif//TSUB_MAKE_TREE == 16
#   if  TSUB_MAKE_TREE ==  8
#define DIV_TSUB_MAKE_TREE(a) ((a) >> 3)
#endif//TSUB_MAKE_TREE ==  8
#   if  TSUB_MAKE_TREE ==  4
#define DIV_TSUB_MAKE_TREE(a) ((a) >> 2)
#endif//TSUB_MAKE_TREE ==  4
#   if  TSUB_MAKE_TREE ==  2
#define DIV_TSUB_MAKE_TREE(a) ((a) >> 1)
#endif//TSUB_MAKE_TREE ==  2
#   if  TSUB_MAKE_TREE ==  1
#define DIV_TSUB_MAKE_TREE(a) (a)
#endif//TSUB_MAKE_TREE ==  1


#define NGROUPS_MAKE_TREE (DIV_TSUB_MAKE_TREE(NTHREADS_MAKE_TREE))


/**
 * @def NTHREADS_INIT_LINK
 *
 * @brief number of threads per block for initTreeLink_kernel
 */
#ifndef NTHREADS_INIT_LINK
#   if  (GPUVER >= 80)
#define NTHREADS_INIT_LINK (512)
#else///(GPUVER >= 80)
#   if  (GPUVER >= 30)
#   if  (GPUVER == 60) || (GPUVER == 61)
#define NTHREADS_INIT_LINK (256)
#else///(GPUVER == 60) || (GPUVER == 61)
#define NTHREADS_INIT_LINK (512)
#endif//(GPUVER == 60) || (GPUVER == 61)
#else///(GPUVER >= 30)
#define NTHREADS_INIT_LINK (256)
#endif//(GPUVER >= 30)
#endif//(GPUVER >= 80)
#endif//NTHREADS_INIT_LINK


/**
 * @def NTHREADS_INIT_CELL
 *
 * @brief number of threads per block for initTreeCell_kernel
 */
#ifndef NTHREADS_INIT_CELL
#   if  (GPUVER >= 80)
#define NTHREADS_INIT_CELL (512)
#else///(GPUVER >= 80)
#   if  (GPUVER >= 70)
#define NTHREADS_INIT_CELL (512)
#else///(GPUVER >= 70)
#   if  (GPUVER == 60) || (GPUVER == 61)
#define NTHREADS_INIT_CELL (128)
#else///(GPUVER == 60) || (GPUVER == 61)
#   if  (GPUVER >= 52)
#define NTHREADS_INIT_CELL (256)
#else///(GPUVER >= 52)
#   if  (GPUVER >= 30)
#define NTHREADS_INIT_CELL (128)
#else///(GPUVER >= 30)
#define NTHREADS_INIT_CELL (512)
#endif//(GPUVER >= 30)
#endif//(GPUVER >= 52)
#endif//(GPUVER == 60) || (GPUVER == 61)
#endif//(GPUVER >= 70)
#endif//(GPUVER >= 80)
#endif//NTHREADS_INIT_CELL


/**
 * @def NTHREADS_INIT_NODE
 *
 * @brief number of threads per block for initTreeNode_kernel
 */
#ifndef NTHREADS_INIT_NODE
#   if  GPUVER >= 70
#define NTHREADS_INIT_NODE (512)
#else///GPUVER >= 70
#   if  GPUVER >= 30
#define NTHREADS_INIT_NODE (256)
#else///GPUVER >= 30
#define NTHREADS_INIT_NODE (128)
#endif//GPUVER >= 30
#endif//GPUVER >= 70
#endif//NTHREADS_INIT_NODE


/**
 * @def NTHREADS_INIT_BODY
 *
 * @brief number of threads per block for initTreeBody_kernel
 */
#ifndef NTHREADS_INIT_BODY
#   if  GPUVER >= 80
#define NTHREADS_INIT_BODY (256)
#else///GPUVER >= 80
#define NTHREADS_INIT_BODY (128)
#endif//GPUVER >= 80
#endif//NTHREADS_INIT_BODY


/**
 * @def NTHREADS_COPY_BODY
 *
 * @brief number of threads per block for copyRealBody_kernel
 */
#ifndef NTHREADS_COPY_BODY
#   if  (GPUVER >= 80)
#define NTHREADS_COPY_BODY (1024)
#else///(GPUVER >= 80)
#   if  (GPUVER >= 70)
#define NTHREADS_COPY_BODY (512)
#else///(GPUVER >= 70)
#   if  (GPUVER >= 60)
#define NTHREADS_COPY_BODY (256)
#else///(GPUVER >= 60)
#   if  (GPUVER >= 52)
#define NTHREADS_COPY_BODY (1024)
#else///(GPUVER >= 52)
#   if  (GPUVER >= 30)
#define NTHREADS_COPY_BODY (512)
#else///(GPUVER >= 30)
#define NTHREADS_COPY_BODY (128)
#endif//(GPUVER >= 30)
#endif//(GPUVER >= 52)
#endif//(GPUVER >= 60)
#endif//(GPUVER >= 70)
#endif//(GPUVER >= 80)
#endif//NTHREADS_COPY_BODY


#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
/**
 * @def NTHREADS_OUTFLOW
 *
 * @brief number of threads per block for checkOutflow_kernel
 */
#define NTHREADS_OUTFLOW (1024)

#ifdef  DIV_NTHREADS_OUTFLOW
#undef  DIV_NTHREADS_OUTFLOW
#endif//DIV_NTHREADS_OUTFLOW
#   if  NTHREADS_OUTFLOW == 1024
#define DIV_NTHREADS_OUTFLOW(a) ((a) >> 10)
#endif//NTHREADS_OUTFLOW == 1024
#   if  NTHREADS_OUTFLOW ==  512
#define DIV_NTHREADS_OUTFLOW(a) ((a) >>  9)
#endif//NTHREADS_OUTFLOW ==  512
#   if  NTHREADS_OUTFLOW ==  256
#define DIV_NTHREADS_OUTFLOW(a) ((a) >>  8)
#endif//NTHREADS_OUTFLOW ==  256
#   if  NTHREADS_OUTFLOW ==  128
#define DIV_NTHREADS_OUTFLOW(a) ((a) >>  7)
#endif//NTHREADS_OUTFLOW ==  128

/**
 * @def TSUB_OUTFLOW
 *
 * @brief number of threads that share a common tree cell for checkOutflow_kernel
 */
#define TSUB_OUTFLOW (8)

/** TSUB_OUTFLOW must be equal or greater than NLEAF / 32 to store NLEAF / TSUB_OUTFLOW flags in an integer (32bit) */
#   if  TSUB_OUTFLOW < (NLEAF >> 5)
#undef  TSUB_OUTFLOW
#define TSUB_OUTFLOW   (NLEAF >> 5)
#endif//TSUB_OUTFLOW < (NLEAF >> 5)

#ifdef  DIV_TSUB_OUTFLOW
#undef  DIV_TSUB_OUTFLOW
#endif//DIV_TSUB_OUTFLOW
#   if  TSUB_OUTFLOW == 32
#define DIV_TSUB_OUTFLOW(a) ((a) >> 5)
#endif//TSUB_OUTFLOW == 32
#   if  TSUB_OUTFLOW == 16
#define DIV_TSUB_OUTFLOW(a) ((a) >> 4)
#endif//TSUB_OUTFLOW == 16
#   if  TSUB_OUTFLOW ==  8
#define DIV_TSUB_OUTFLOW(a) ((a) >> 3)
#endif//TSUB_OUTFLOW ==  8
#   if  TSUB_OUTFLOW ==  4
#define DIV_TSUB_OUTFLOW(a) ((a) >> 2)
#endif//TSUB_OUTFLOW ==  4
#   if  TSUB_OUTFLOW ==  2
#define DIV_TSUB_OUTFLOW(a) ((a) >> 1)
#endif//TSUB_OUTFLOW ==  2
#   if  TSUB_OUTFLOW ==  1
#define DIV_TSUB_OUTFLOW(a) (a)
#endif//TSUB_OUTFLOW ==  1

#define NGROUPS_OUTFLOW (DIV_TSUB_OUTFLOW(NTHREADS_OUTFLOW))
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)


/**
 * @struct soaMakeTreeBuf
 *
 * @brief structure for building octree structure (SoA)
 */
typedef struct
{
  size_t Nbuf;
  real *rjmax;
  int *more0, *more1;
  int *fail;
  int *gsync0, *gsync1;
  int *gmem_make_tree, *gmem_link_tree;
  int *gsync0_make_tree, *gsync1_make_tree, *gsync2_make_tree, *gsync3_make_tree;
  int *gsync0_link_tree, *gsync1_link_tree;
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
#   if  defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
  real *rbuf_external;
  int  *ibuf_external;
#else///defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
  uint *ubuf_external;
#endif//defined(USE_PARENT_MAC_FOR_EXTERNAL_PARTICLES) || defined(TIME_BASED_MODIFICATION)
  size_t Nbuf_external;
  int *gmem_external, *gsync0_external, *gsync1_external;
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
#ifdef  USE_OCCUPANCY_CALCULATOR
  int numBlocksPerSM_mac;
  int numBlocksPerSM_make_tree;
  int numBlocksPerSM_link_tree;
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
  int numBlocksPerSM_outflow;
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
#endif//USE_OCCUPANCY_CALCULATOR
} soaMakeTreeBuf;


/* list of functions appeared in ``make_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
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
   );

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
   int **more0Buf, int **more1Buf, real **rjmaxBuf, int **fail_dev, soaMakeTreeBuf *buf, const deviceProp devProp);

  void  freeTreeNode_dev
  (uint  *more_dev, jparticle  *pj_dev, jmass  *mj_dev, real  *bmax_dev, int  *n2c_dev, int  *gsync0, int  *gsync1,
#       ifdef  WS93_MAC
   real  *mr2_dev,
#       endif//WS93_MAC
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
   int  *more0Buf, int  *more1Buf, real  *rjmaxBuf, int  *fail_dev);

  muse allocTreeCell_dev
  (soaTreeCell *dev, treecell **cell_dev, bool **leaf_dev, uint **node_dev, PHinfo **info_dev,
   PHint **hkey_dev, uint **parent_dev, uint **children_dev, int **leafLev_dev, int **numCell_dev, int **numNode_dev, int **scanNum_dev
#ifdef  COUNT_INTERACTIONS
   , soaTreeCell *hst, treecell **cell_hst, bool **leaf_hst, uint **node_hst, PHinfo **info_hst
#endif//COUNT_INTERACTIONS
#   if  !defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
   , deviceProp devProp
#endif//!defined(SERIALIZED_EXECUTION) && defined(CARE_EXTERNAL_PARTICLES)
   );

  void  freeTreeCell_dev
  (treecell  *cell_dev, bool  *leaf_dev, uint  *node_dev, PHinfo  *info_dev,
   PHint  *hkey_dev, uint  *parent_dev, uint  *children_dev, int  *leafLev_dev, int  *numCell_dev, int  *numNode_dev, int  *scanNum_dev
#ifdef  COUNT_INTERACTIONS
   , treecell  *cell_hst, bool  *leaf_hst, uint  *node_hst, PHinfo  *info_hst
#endif//COUNT_INTERACTIONS
   );

#ifdef  GADGET_MAC
  void enforceBarnesHutMAC_dev(const int Ni, const iparticle pi, const int Nj, const soaTreeNode pj);
  void recoverGADGET_MAC_dev(const int Nj, const soaTreeNode pj);
#endif//GADGET_MAC

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
   );

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
   );

#ifndef NDEBUG
  void printPHkey_location(const int piNum, PHint * RESTRICT peano, const iparticle pi);
#endif//NDEBUG

#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//MAKE_DEV_H
