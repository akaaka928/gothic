/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/21(Thu) 16:21:55
 *                                                                       *
 *    Header File for constructing octree structure                      *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef MAKE_DEV_H
#define MAKE_DEV_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  !defined(_STDBOOL_H) && !defined(_STDBOOL)
#       include <stdbool.h>
#endif//!defined(_STDBOOL_H) && !defined(_STDBOOL)
//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef BENCHMARK_H
#       include "../misc/benchmark.h"
#endif//BENCHMARK_H
//-------------------------------------------------------------------------
#ifndef STRUCTURE_H
#       include "../misc/structure.h"
#endif//STRUCTURE_H
//-------------------------------------------------------------------------
#ifndef PEANO_H
#       include "../sort/peano.h"
#endif//PEANO_H
//-------------------------------------------------------------------------
#ifndef MACUTIL_H
#       include "../tree/macutil.h"
#endif//MACUTIL_H
//-------------------------------------------------------------------------
#ifndef MAKE_H
#       include "../tree/make.h"
#endif//MAKE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef HUNT_NODE_PARAMETER
/* the below macro is enabled in the default option; switched off in the parameter survey mode to use -D and -U from Makefile */
#define USE_WARP_SHUFFLE_FUNC_MAC
#endif//HUNT_NODE_PARAMETER
//-------------------------------------------------------------------------
#   if  defined(USE_WARP_SHUFFLE_FUNC_MAC) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_MAC
#endif//defined(USE_WARP_SHUFFLE_FUNC_MAC) && (GPUGEN < 30)
//-------------------------------------------------------------------------
#ifndef NTHREADS_MAC
#          if  (GPUGEN >= 52)
#define NTHREADS_MAC (256)
#       else///(GPUGEN >= 52)
#define NTHREADS_MAC (128)
#       endif//(GPUGEN >= 52)
#endif//NTHREADS_MAC
/* NTHREADS_MAC must be equal or smaller than 512 due to the capacity of shared memory */
#   if  NTHREADS_MAC > 512
#undef  NTHREADS_MAC
#define NTHREADS_MAC  (512)
#endif//NTHREADS_MAC > 512
//-------------------------------------------------------------------------
#ifndef TSUB_MAC
#define TSUB_MAC (32)
#endif//TSUB_MAC
/* TSUB_MAC must be equal or smaller than NTHREADS_MAC */
#   if  TSUB_MAC > NTHREADS_MAC
#undef  TSUB_MAC
#define TSUB_MAC   NTHREADS_MAC
#endif//TSUB_MAC > NTHREADS_MAC
//-------------------------------------------------------------------------
#   if  NTHREADS_MAC >= 128
/* TSUB_MAC must be equal or greater than 4 */
#          if  TSUB_MAC < 4
#       undef  TSUB_MAC
#       define TSUB_MAC  (4)
#       endif//TSUB_MAC < 4
#   if  NTHREADS_MAC >= 256
/* TSUB_MAC must be equal or greater than 8 */
#          if  TSUB_MAC < 8
#       undef  TSUB_MAC
#       define TSUB_MAC  (8)
#       endif//TSUB_MAC < 8
#   if  NTHREADS_MAC == 512
/* TSUB_MAC must be equal to 32 */
#          if  TSUB_MAC != 32
#       undef  TSUB_MAC
#       define TSUB_MAC   (32)
#       endif//TSUB_MAC != 32
#endif//NTHREADS_MAC == 512
#endif//NTHREADS_MAC >= 256
#endif//NTHREADS_MAC >= 128
//-------------------------------------------------------------------------
#define NGROUPS_MAC (NTHREADS_MAC / TSUB_MAC)
//-------------------------------------------------------------------------
#define NBUF_MAC (((NJ_BMAX_ESTIMATE) + (TSUB_MAC) - 1) / TSUB_MAC)
//-------------------------------------------------------------------------
/* maximum number of NBUF_MAC is 4 to use float4 in union */
#   if  NBUF_MAC > 4
#undef  NBUF_MAC
#define NBUF_MAC (4)
#undef  TSUB_MAC
#define TSUB_MAC ((NJ_BMAX_ESTIMATE + (NBUF_MAC) - 1) / NBUF_MAC)
#endif//NBUF_MAC > 4
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  NI_BMAX_ESTIMATE > (TSUB_MAC * NBUF_MAC)
#undef  NI_BMAX_ESTIMATE
#define NI_BMAX_ESTIMATE   (TSUB_MAC * NBUF_MAC)
#endif//NI_BMAX_ESTIMATE > (TSUB_MAC * NBUF_MAC)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------
#ifndef HUNT_MAKE_PARAMETER
/* the below macro is disabled in the default option for better performance; switched off in the parameter survey mode to use -D from Makefile */
/* #define USE_WARP_SHUFFLE_FUNC_MAKE_TREE */
/* the below macro is disabled in the default option for better performance; switched off in the parameter survey mode to use -D from Makefile */
/* #define USE_WARP_SHUFFLE_FUNC_LINK_TREE */
#endif//HUNT_MAKE_PARAMETER
//-------------------------------------------------------------------------
#   if  defined(USE_WARP_SHUFFLE_FUNC_MAKE_TREE) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_MAKE_TREE
#endif//defined(USE_WARP_SHUFFLE_FUNC_MAKE_TREE) && (GPUGEN < 30)
//-------------------------------------------------------------------------
#   if  defined(USE_WARP_SHUFFLE_FUNC_LINK_TREE) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_LINK_TREE
#endif//defined(USE_WARP_SHUFFLE_FUNC_LINK_TREE) && (GPUGEN < 30)
//-------------------------------------------------------------------------
#ifndef NTHREADS_MAKE_TREE
#define NTHREADS_MAKE_TREE (128)
#endif//NTHREADS_MAKE_TREE
#ifndef NTHREADS_LINK_TREE
#          if  (GPUGEN >= 30)
#define NTHREADS_LINK_TREE (256)
#       else///(GPUGEN >= 30)
#define NTHREADS_LINK_TREE (128)
#       endif//(GPUGEN >= 30)
#endif//NTHREADS_LINK_TREE
#ifndef NTHREADS_TRIM_TREE
#define NTHREADS_TRIM_TREE (128)
#endif//NTHREADS_TRIM_TREE
//-------------------------------------------------------------------------
/* NTHREADS_MAKE_TREE must be equal or smaller than 512 due to the capacity of shared memory */
#   if  NTHREADS_MAKE_TREE > 512
#undef  NTHREADS_MAKE_TREE
#define NTHREADS_MAKE_TREE  (512)
#endif//NTHREADS_MAKE_TREE > 512
//-------------------------------------------------------------------------
#ifndef TSUB_MAKE_TREE
#define TSUB_MAKE_TREE   CELL_UNIT
#endif//TSUB_MAKE_TREE
#   if  TSUB_MAKE_TREE > CELL_UNIT
#undef  TSUB_MAKE_TREE
#define TSUB_MAKE_TREE   CELL_UNIT
#endif//TSUB_MAKE_TREE > CELL_UNIT
#define NGROUPS_MAKE_TREE (NTHREADS_MAKE_TREE / TSUB_MAKE_TREE)
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
#ifndef NTHREADS_INIT_LINK
#          if  (GPUGEN >= 30)
#define NTHREADS_INIT_LINK (512)
#       else///(GPUGEN >= 30)
#define NTHREADS_INIT_LINK (256)
#       endif//(GPUGEN >= 30)
#endif//NTHREADS_INIT_LINK
//-------------------------------------------------------------------------
#ifndef NTHREADS_INIT_CELL
#          if  (GPUGEN >= 52)
#define NTHREADS_INIT_CELL (256)
#       else///(GPUGEN >= 52)
#          if  (GPUGEN >= 30)
#define NTHREADS_INIT_CELL (128)
#       else///(GPUGEN >= 30)
#define NTHREADS_INIT_CELL (512)
#       endif//(GPUGEN >= 30)
#       endif//(GPUGEN >= 52)
#endif//NTHREADS_INIT_CELL
//-------------------------------------------------------------------------
#ifndef NTHREADS_INIT_NODE
#          if  (GPUGEN >= 30)
#define NTHREADS_INIT_NODE (256)
#       else///(GPUGEN >= 30)
#define NTHREADS_INIT_NODE (128)
#       endif//(GPUGEN >= 30)
#endif//NTHREADS_INIT_NODE
//-------------------------------------------------------------------------
#endif//MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef NTHREADS_INIT_BODY
#define NTHREADS_INIT_BODY (128)
#endif//NTHREADS_INIT_BODY
//-------------------------------------------------------------------------
#ifndef NTHREADS_COPY_BODY
#          if  (GPUGEN >= 52)
#define NTHREADS_COPY_BODY (1024)
#       else///(GPUGEN >= 52)
#          if  (GPUGEN >= 30)
#define NTHREADS_COPY_BODY (512)
#       else///(GPUGEN >= 30)
#define NTHREADS_COPY_BODY (128)
#       endif//(GPUGEN >= 30)
#       endif//(GPUGEN >= 52)
#endif//NTHREADS_COPY_BODY
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
  int *more0, *more1;
  real *rjmax;
  int *fail;
  int *gsync0, *gsync1;
#ifdef  MAKE_TREE_ON_DEVICE
  int *gmem_make_tree, *gmem_link_tree;
  int *gsync0_make_tree, *gsync1_make_tree, *gsync2_make_tree, *gsync3_make_tree;
  int *gsync0_link_tree, *gsync1_link_tree;
#endif//MAKE_TREE_ON_DEVICE
} soaMakeTreeBuf;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "make_dev.cu"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
  void makeTreeStructure_dev
  (const int piNum, PHint * RESTRICT peano,
   int * RESTRICT leafLev, int * RESTRICT leafLev_dev, int * RESTRICT scanNum_dev,
   int * RESTRICT cellNum, int * RESTRICT cellNum_dev, const soaTreeCell cell,
   int * RESTRICT nodeNum, int * RESTRICT nodeNum_dev, const soaTreeNode node,
   const soaMakeTreeBuf buf, deviceProp devProp
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
   );
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
/*   muse allocTreeNode_dev */
/*   (uint **more_dev, jparticle **pj_dev, jmass **mj_dev, */
/* #ifndef MAKE_TREE_ON_DEVICE */
/*    uint **more_hst, int **n2c_hst, */
/* #endif//MAKE_TREE_ON_DEVICE */
/* #ifdef  GADGET_MAC */
/*    real **mac_dev, */
/* #endif//GADGET_MAC */
/* #ifdef  CALC_MULTIPOLE_ON_DEVICE */
/* #       ifdef  WS93_MAC */
/*    real **mr2_dev, */
/* #       endif//WS93_MAC */
/*    real **bmax_dev, int **n2c_dev, int **gsync0, int **gsync1, deviceProp devProp, */
/* #else///CALC_MULTIPOLE_ON_DEVICE */
/* #       ifdef  WS93_MAC */
/*    real **mr2_hst, */
/* #       endif//WS93_MAC */
/*    real **bmax_hst, */
/* #endif//CALC_MULTIPOLE_ON_DEVICE */
/* #   if  !defined(CALC_MULTIPOLE_ON_DEVICE) || (!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)) */
/*    jparticle **pj_hst, jmass **mj_hst, */
/* #endif//!defined(CALC_MULTIPOLE_ON_DEVICE) || (!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)) */
/* #   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH) */
/*    int **niSub_dev, */
/* #ifndef MAKE_TREE_ON_DEVICE */
/*    int **niSub_hst, */
/* #endif//MAKE_TREE_ON_DEVICE */
/* #endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH) */
/* #ifdef  MAKE_TREE_ON_DEVICE */
/*    int **gmem_make_tree, int **gsync0_make_tree, int **gsync1_make_tree, int **gsync2_make_tree, int **gsync3_make_tree, */
/*    int **gmem_link_tree, int **gsync0_link_tree, int **gsync1_link_tree, */
/* #endif//MAKE_TREE_ON_DEVICE */
/*    int **more0Buf, int **more1Buf, real **rjmaxBuf, int **fail_dev, soaTreeNode *dev, soaTreeNode *hst, soaMakeTreeBuf *buf); */
/*   void  freeTreeNode_dev */
/*   (uint  *more_dev, jparticle  *pj_dev, jmass  *mj_dev, */
/* #ifndef MAKE_TREE_ON_DEVICE */
/*    uint  *more_hst, int  *n2c_hst, */
/* #endif//MAKE_TREE_ON_DEVICE */
/* #ifdef  GADGET_MAC */
/*    real  *mac_dev, */
/* #endif//GADGET_MAC */
/* #ifdef  CALC_MULTIPOLE_ON_DEVICE */
/* #       ifdef  WS93_MAC */
/*    real  *mr2_dev, */
/* #       endif//WS93_MAC */
/*    real  *bmax_dev, int  *n2c_dev, int  *gsync0, int  *gsync1, */
/* #else///CALC_MULTIPOLE_ON_DEVICE */
/* #       ifdef  WS93_MAC */
/*    real  *mr2_hst, */
/* #       endif//WS93_MAC */
/*    real  *bmax_hst, */
/* #endif//CALC_MULTIPOLE_ON_DEVICE */
/* #   if  !defined(CALC_MULTIPOLE_ON_DEVICE) || (!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)) */
/*    jparticle  *pj_hst, jmass  *mj_hst, */
/* #endif//!defined(CALC_MULTIPOLE_ON_DEVICE) || (!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)) */
/* #   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH) */
/*    int  *niSub_dev, */
/* #ifndef MAKE_TREE_ON_DEVICE */
/*    int  *niSub_hst, */
/* #endif//MAKE_TREE_ON_DEVICE */
/* #endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH) */
/* #ifdef  MAKE_TREE_ON_DEVICE */
/*    int  *gmem_make_tree, int  *gsync0_make_tree, int  *gsync1_make_tree, int  *gsync2_make_tree, int  *gsync3_make_tree, */
/*    int  *gmem_link_tree, int  *gsync0_link_tree, int  *gsync1_link_tree, */
/* #endif//MAKE_TREE_ON_DEVICE */
/*    int  *more0Buf, int  *more1Buf, real  *rjmaxBuf, int  *fail_dev); */
  muse allocTreeNode_dev
  (uint **more_dev, jparticle **pj_dev, jmass **mj_dev,
#ifdef  CALC_MULTIPOLE_ON_DEVICE
   real **bmax_dev, int **n2c_dev, int **gsync0, int **gsync1, deviceProp devProp,
#       ifdef  WS93_MAC
   real **mr2_dev,
#       endif//WS93_MAC
#else///CALC_MULTIPOLE_ON_DEVICE
   real **bmax_hst,
#       ifdef  WS93_MAC
   real **mr2_hst,
#       endif//WS93_MAC
#endif//CALC_MULTIPOLE_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
   int **niSub_dev,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
   uint **more_hst,
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
   jparticle **pj_hst, jmass **mj_hst,
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
#ifdef  MAKE_TREE_ON_DEVICE
   int **gmem_make_tree, int **gsync0_make_tree, int **gsync1_make_tree, int **gsync2_make_tree, int **gsync3_make_tree,
   int **gmem_link_tree, int **gsync0_link_tree, int **gsync1_link_tree,
#else///MAKE_TREE_ON_DEVICE
   int **n2c_hst,
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
   int **niSub_hst,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#endif//MAKE_TREE_ON_DEVICE
#ifdef  GADGET_MAC
   real **mac_dev,
#endif//GADGET_MAC
   int **more0Buf, int **more1Buf, real **rjmaxBuf, int **fail_dev, soaTreeNode *dev, soaTreeNode *hst, soaMakeTreeBuf *buf);
  void  freeTreeNode_dev
  (uint  *more_dev, jparticle  *pj_dev, jmass  *mj_dev,
#ifdef  CALC_MULTIPOLE_ON_DEVICE
   real  *bmax_dev, int  *n2c_dev, int  *gsync0, int  *gsync1,
#       ifdef  WS93_MAC
   real  *mr2_dev,
#       endif//WS93_MAC
#else///CALC_MULTIPOLE_ON_DEVICE
   real  *bmax_hst,
#       ifdef  WS93_MAC
   real  *mr2_hst,
#       endif//WS93_MAC
#endif//CALC_MULTIPOLE_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
   int  *niSub_dev,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
   uint  *more_hst,
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
   jparticle  *pj_hst, jmass  *mj_hst,
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
#ifdef  MAKE_TREE_ON_DEVICE
   int  *gmem_make_tree, int  *gsync0_make_tree, int  *gsync1_make_tree, int  *gsync2_make_tree, int  *gsync3_make_tree,
   int  *gmem_link_tree, int  *gsync0_link_tree, int  *gsync1_link_tree,
#else///MAKE_TREE_ON_DEVICE
   int  *n2c_hst,
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
   int  *niSub_hst,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#endif//MAKE_TREE_ON_DEVICE
#ifdef  GADGET_MAC
   real  *mac_dev,
#endif//GADGET_MAC
   int  *more0Buf, int  *more1Buf, real  *rjmaxBuf, int  *fail_dev);
  //-----------------------------------------------------------------------
  void setTreeNode_dev
  (const size_t Nj, const soaTreeNode dev, const soaTreeNode hst
#ifdef  CALC_MULTIPOLE_ON_DEVICE
   , const size_t Ni
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
   );
  //-----------------------------------------------------------------------
#   if  defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
  muse allocTreeCell_dev
  (treecell **cell_dev, bool **leaf_dev, uint **node_dev, PHinfo **info_dev,
#ifdef  MAKE_TREE_ON_DEVICE
   PHint **hkey_dev, uint **parent_dev, uint **children_dev, int **leafLev_dev, int **numCell_dev, int **numNode_dev, int **scanNum_dev,
#else///MAKE_TREE_ON_DEVICE
   treecell **cell_hst, bool **leaf_hst, uint **node_hst, PHinfo **info_hst,
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
   treecell **cell_hst, bool **leaf_hst, uint **node_hst, PHinfo **info_hst,
#endif//defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
   soaTreeCell *dev, soaTreeCell *hst);
  void  freeTreeCell_dev
  (treecell  *cell_dev, bool  *leaf_dev, uint  *node_dev, PHinfo  *info_dev,
#ifdef  MAKE_TREE_ON_DEVICE
   PHint  *hkey_dev, uint  *parent_dev, uint  *children_dev, int  *leafLev_dev, int  *numCell_dev, int  *numNode_dev, int  *scanNum_dev
#else///MAKE_TREE_ON_DEVICE
   treecell  *cell_hst, bool  *leaf_hst, uint  *node_hst, PHinfo  *info_hst
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
   , treecell  *cell_hst, bool  *leaf_hst, uint  *node_hst, PHinfo  *info_hst
#endif//defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
   );
#endif//defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
  //-----------------------------------------------------------------------
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifndef MAKE_TREE_ON_DEVICE
/*   void setTreeCell_dev */
/*   (const size_t num, */
/*    treecell *cell_dev, bool *leaf_dev, uint *node_dev, PHinfo *info_dev, */
/*    treecell *cell_hst, bool *leaf_hst, uint *node_hst, PHinfo *info_hst */
/* #ifdef  EXEC_BENCHMARK */
/*    , wall_clock_time *elapsed */
/* #endif//EXEC_BENCHMARK */
/*    ); */
  void setTreeCell_dev(const size_t num, const soaTreeCell dev, const soaTreeCell hst
#ifdef  EXEC_BENCHMARK
		       , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		       );
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifdef  GADGET_MAC
  /* void enforceBarnesHutMAC_dev */
  /* (const int Ni, acceleration * RESTRICT ai_dev, */
  /*  const int Nj, jparticle * RESTRICT pj_dev, real * RESTRICT mac_dev, real * RESTRICT bmax_dev); */
  /* void recoverGADGET_MAC_dev */
  /* (const int Nj, jparticle * RESTRICT pj_dev, real * RESTRICT mac_dev); */
  void enforceBarnesHutMAC_dev(const int Ni, const iparticle pi, const int Nj, const soaTreeNode pj);
  void recoverGADGET_MAC_dev(const int Nj, const soaTreeNode pj);
#endif//GADGET_MAC
  //-----------------------------------------------------------------------
/*   void calcMultipole_dev */
/*   (const int bottomLev, PHinfo * RESTRICT level, PHinfo * RESTRICT level_dev, treecell * RESTRICT cell, bool * RESTRICT leaf, */
/*    const int piNum, position * RESTRICT pi, int * RESTRICT jtag, */
/*    const int pjNum, uint * RESTRICT node, uint * RESTRICT more, int * RESTRICT node2cell, jparticle * RESTRICT pj, jmass * RESTRICT mj, real * RESTRICT bmax, */
/*    int * RESTRICT more0Buf, int * RESTRICT more1Buf, real * RESTRICT rjmaxBuf, int * RESTRICT fail_dev, */
/*    int * RESTRICT gsync0, int * RESTRICT gsync1, deviceProp devProp */
/* #ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING */
/*    , const real eps2 */
/* #endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING */
/* #ifdef  WS93_MAC */
/*    , real * RESTRICT mr2 */
/* #endif//WS93_MAC */
/* #   if  !defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE) */
/*    , jparticle * RESTRICT pj_hst, jmass * RESTRICT mj_hst, real * RESTRICT bmax_root_hst */
/* #endif//!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE) */
/* #ifdef  COUNT_INTERACTIONS */
/*    , tree_stats * RESTRICT stats_hst */
/* #endif//COUNT_INTERACTIONS */
/* #ifdef  EXEC_BENCHMARK */
/*    , wall_clock_time *elapsed */
/* #endif//EXEC_BENCHMARK */
/*    ); */
  void calcMultipole_dev
  (const int bottomLev, const soaTreeCell cell,
   const int piNum, const iparticle pi, const int pjNum, const soaTreeNode node,
   const soaMakeTreeBuf buf, deviceProp devProp
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
   , const real eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#   if  !defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
   , const soaTreeNode node_hst, real * RESTRICT bmax_root_hst
#endif//!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
#ifdef  COUNT_INTERACTIONS
   , const soaTreeCell cell_hst, tree_stats * RESTRICT stats_hst
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
   );
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//MAKE_DEV_H
//-------------------------------------------------------------------------
