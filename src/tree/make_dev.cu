/*************************************************************************\
 *                                                                       *
                  last updated on 2016/10/16(Sun) 16:40:38
 *                                                                       *
 *    Constructing octree structure for collisionless systems            *
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
#include <helper_cuda.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "../sort/peano.h"
//-------------------------------------------------------------------------
#include "macutil.h"
#include "make.h"
#include "make_dev.h"
#include "let.h"/* <-- necessary to read EXTEND_NUM_TREE_NODE */
//-------------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
#include "../sort/peano_dev.h"
#endif//MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------
#include "../misc/gsync_dev.cu"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* limitation from capacity of shared memory */
/* in shared memory preferred configuration, capacity of shared memory is 48KiB per SM */
#       ifdef  USE_WARP_SHUFFLE_FUNC_MAC
/* 3072 = 48 * 1024 / 16 */
#define NBLOCKS_PER_SM_MAC ( 3072 / (NTHREADS_MAC * (1 +     NBUF_MAC)))
#       else///USE_WARP_SHUFFLE_FUNC_MAC
/* 12288 = 48 * 1024 / 4 */
#define NBLOCKS_PER_SM_MAC (12288 / (NTHREADS_MAC * (5 + 4 * NBUF_MAC)))
#       endif//USE_WARP_SHUFFLE_FUNC_MAC
//-------------------------------------------------------------------------
#define REGISTERS_PER_THREAD_MAC (64)
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
#       ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#ifdef  GADGET_MAC
#define REGISTERS_PER_THREAD_MAC (55)
#else///GADGET_MAC
#ifdef  WS93_MAC
#define REGISTERS_PER_THREAD_MAC (57)
#else///WS93_MAC
#define REGISTERS_PER_THREAD_MAC (49)
#endif//WS93_MAC
#endif//GADGET_MAC
#       else///USE_WARP_SHUFFLE_FUNC_MAC
#define REGISTERS_PER_THREAD_MAC (62)
#       endif//USE_WARP_SHUFFLE_FUNC_MAC
#endif//GPUVER == 35
/* calcMultipole_kernel uses 64 registers @ GTX 750 Ti w/z warp shuffle */
/* calcMultipole_kernel uses 72 registers @ GTX 750 Ti w/o warp shuffle */
#   if  GPUVER == 50
#undef  REGISTERS_PER_THREAD_MAC
#       ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#define REGISTERS_PER_THREAD_MAC (64)
#       else///USE_WARP_SHUFFLE_FUNC_MAC
#define REGISTERS_PER_THREAD_MAC (72)
#       endif//USE_WARP_SHUFFLE_FUNC_MAC
#endif//GPUVER == 50
/* calcMultipole_kernel uses 64 registers @ GTX 970 w/z warp shuffle */
/* calcMultipole_kernel uses 72 registers @ GTX 970 w/o warp shuffle */
#   if  GPUVER == 52
#undef  REGISTERS_PER_THREAD_MAC
#       ifdef  USE_WARP_SHUFFLE_FUNC_MAC
#define REGISTERS_PER_THREAD_MAC (64)
#       else///USE_WARP_SHUFFLE_FUNC_MAC
#define REGISTERS_PER_THREAD_MAC (72)
#       endif//USE_WARP_SHUFFLE_FUNC_MAC
#endif//GPUVER == 52
//-------------------------------------------------------------------------
/* limitation from number of registers */
#   if  NBLOCKS_PER_SM_MAC > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAC * NTHREADS_MAC))
#undef  NBLOCKS_PER_SM_MAC
#define NBLOCKS_PER_SM_MAC   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAC * NTHREADS_MAC))
#endif//NBLOCKS_PER_SM_MAC > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAC * NTHREADS_MAC))
//-------------------------------------------------------------------------
/* maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_MAC > (MAX_THREADS_PER_SM / NTHREADS_MAC)
#undef  NBLOCKS_PER_SM_MAC
#define NBLOCKS_PER_SM_MAC   (MAX_THREADS_PER_SM / NTHREADS_MAC)
#endif//NBLOCKS_PER_SM_MAC > (MAX_THREADS_PER_SM / NTHREADS_MAC)
//-------------------------------------------------------------------------
/* maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_MAC > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_MAC
#define NBLOCKS_PER_SM_MAC   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_MAC > MAX_BLOCKS_PER_SM
//-------------------------------------------------------------------------
/* maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_MAC > ((MAX_WARPS_PER_SM * 32) / NTHREADS_MAC)
#undef  NBLOCKS_PER_SM_MAC
#define NBLOCKS_PER_SM_MAC   ((MAX_WARPS_PER_SM * 32) / NTHREADS_MAC)
#endif//NBLOCKS_PER_SM_MAC > ((MAX_WARPS_PER_SM * 32) / NTHREADS_MAC)
//-------------------------------------------------------------------------
/* # of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_MAC < 1
#undef  NBLOCKS_PER_SM_MAC
#define NBLOCKS_PER_SM_MAC  (1)
#endif//NBLOCKS_PER_SM_MAC < 1
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------
/* limitations from capacity of shared memory (assuming L1 cache preferred) */
/* SM usage is 20 * NTHREADS_MAKE_TREE + 16 * NUM_PHKEY_LEVEL bytes */
/* #   if  NBLOCKS_PER_SM_MAKE_TREE > (16384 / (20 * NTHREADS_MAKE_TREE + 16 * NUM_PHKEY_LEVEL)) */
/* #undef  NBLOCKS_PER_SM_MAKE_TREE */
#define NBLOCKS_PER_SM_MAKE_TREE   (16384 / (20 * NTHREADS_MAKE_TREE + 16 * NUM_PHKEY_LEVEL))
/* #endif//NBLOCKS_PER_SM_MAKE_TREE > (16384 / (20 * NTHREADS_MAKE_TREE + 16 * NUM_PHKEY_LEVEL)) */
/* SM usage is  4 * NTHREADS_LINK_TREE bytes */
/* #   if  NBLOCKS_PER_SM_LINK_TREE > (16384 / (4 * NTHREADS_LINK_TREE)) */
/* #undef  NBLOCKS_PER_SM_LINK_TREE */
#define NBLOCKS_PER_SM_LINK_TREE   (16384 / (4 * NTHREADS_LINK_TREE))
/* #endif//NBLOCKS_PER_SM_LINK_TREE > (16384 / (4 * NTHREADS_LINK_TREE)) */
//-------------------------------------------------------------------------
#define REGISTERS_PER_THREAD_MAKE_TREE (64)
#define REGISTERS_PER_THREAD_LINK_TREE (32)
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
#          if  NTHREADS_MAKE_TREE == 256
#define REGISTERS_PER_THREAD_MAKE_TREE (66)
#       else///NTHREADS_MAKE_TREE == 256
#define REGISTERS_PER_THREAD_MAKE_TREE (64)
#       endif//NTHREADS_MAKE_TREE == 256
#undef  REGISTERS_PER_THREAD_LINK_TREE
#define REGISTERS_PER_THREAD_LINK_TREE (23)
#endif//GPUVER == 50
/* makeTree_kernel and linkTree_kernel use 64 and 23 registers, respectively @ GTX 970 w/z Ttot = 128, 512, 1024 @ CUDA 7.5 */
/* makeTree_kernel and linkTree_kernel use 66 and 23 registers, respectively @ GTX 970 w/z Ttot = 256            @ CUDA 7.5 */
/* makeTree_kernel use 70 registers @ GTX TITAN X w/z Ttot = 128 @ CUDA 8.0 */
/* linkTree_kernel use 24 registers @ GTX TITAN X w/z Ttot = 256 @ CUDA 8.0 */
#   if  GPUVER == 52
#undef  REGISTERS_PER_THREAD_MAKE_TREE
#          if  NTHREADS_MAKE_TREE == 256
#define REGISTERS_PER_THREAD_MAKE_TREE (66)
#       else///NTHREADS_MAKE_TREE == 256
/* #define REGISTERS_PER_THREAD_MAKE_TREE (64) */
#define REGISTERS_PER_THREAD_MAKE_TREE (70)
#       endif//NTHREADS_MAKE_TREE == 256
#undef  REGISTERS_PER_THREAD_LINK_TREE
/* #define REGISTERS_PER_THREAD_LINK_TREE (23) */
#define REGISTERS_PER_THREAD_LINK_TREE (24)
#endif//GPUVER == 52
//-------------------------------------------------------------------------
/* limitations from number of registers */
#   if  NBLOCKS_PER_SM_MAKE_TREE > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAKE_TREE * NTHREADS_MAKE_TREE))
#undef  NBLOCKS_PER_SM_MAKE_TREE
#define NBLOCKS_PER_SM_MAKE_TREE   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAKE_TREE * NTHREADS_MAKE_TREE))
#endif//NBLOCKS_PER_SM_MAKE_TREE > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_MAKE_TREE * NTHREADS_MAKE_TREE))
#   if  NBLOCKS_PER_SM_LINK_TREE > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_LINK_TREE * NTHREADS_LINK_TREE))
#undef  NBLOCKS_PER_SM_LINK_TREE
#define NBLOCKS_PER_SM_LINK_TREE   (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_LINK_TREE * NTHREADS_LINK_TREE))
#endif//NBLOCKS_PER_SM_LINK_TREE > (MAX_REGISTERS_PER_SM / (REGISTERS_PER_THREAD_LINK_TREE * NTHREADS_LINK_TREE))
//-------------------------------------------------------------------------
/* maximum # of registers per SM */
#   if  NBLOCKS_PER_SM_MAKE_TREE > (MAX_THREADS_PER_SM / NTHREADS_MAKE_TREE)
#undef  NBLOCKS_PER_SM_MAKE_TREE
#define NBLOCKS_PER_SM_MAKE_TREE   (MAX_THREADS_PER_SM / NTHREADS_MAKE_TREE)
#endif//NBLOCKS_PER_SM_MAKE_TREE > (MAX_THREADS_PER_SM / NTHREADS_MAKE_TREE)
#   if  NBLOCKS_PER_SM_LINK_TREE > (MAX_THREADS_PER_SM / NTHREADS_LINK_TREE)
#undef  NBLOCKS_PER_SM_LINK_TREE
#define NBLOCKS_PER_SM_LINK_TREE   (MAX_THREADS_PER_SM / NTHREADS_LINK_TREE)
#endif//NBLOCKS_PER_SM_LINK_TREE > (MAX_THREADS_PER_SM / NTHREADS_LINK_TREE)
//-------------------------------------------------------------------------
/* maximum # of resident blocks per SM */
#   if  NBLOCKS_PER_SM_MAKE_TREE > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_MAKE_TREE
#define NBLOCKS_PER_SM_MAKE_TREE   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_MAKE_TREE > MAX_BLOCKS_PER_SM
#   if  NBLOCKS_PER_SM_LINK_TREE > MAX_BLOCKS_PER_SM
#undef  NBLOCKS_PER_SM_LINK_TREE
#define NBLOCKS_PER_SM_LINK_TREE   MAX_BLOCKS_PER_SM
#endif//NBLOCKS_PER_SM_LINK_TREE > MAX_BLOCKS_PER_SM
//-------------------------------------------------------------------------
/* maximum # of resident warps per SM */
#   if  NBLOCKS_PER_SM_MAKE_TREE > ((MAX_WARPS_PER_SM * 32) / NTHREADS_MAKE_TREE)
#undef  NBLOCKS_PER_SM_MAKE_TREE
#define NBLOCKS_PER_SM_MAKE_TREE   ((MAX_WARPS_PER_SM * 32) / NTHREADS_MAKE_TREE)
#endif//NBLOCKS_PER_SM_MAKE_TREE > ((MAX_WARPS_PER_SM * 32) / NTHREADS_MAKE_TREE)
#   if  NBLOCKS_PER_SM_LINK_TREE > ((MAX_WARPS_PER_SM * 32) / NTHREADS_LINK_TREE)
#undef  NBLOCKS_PER_SM_LINK_TREE
#define NBLOCKS_PER_SM_LINK_TREE   ((MAX_WARPS_PER_SM * 32) / NTHREADS_LINK_TREE)
#endif//NBLOCKS_PER_SM_LINK_TREE > ((MAX_WARPS_PER_SM * 32) / NTHREADS_LINK_TREE)
//-------------------------------------------------------------------------
/* # of blocks per SM must not be zero */
#   if  NBLOCKS_PER_SM_MAKE_TREE < 1
#undef  NBLOCKS_PER_SM_MAKE_TREE
#define NBLOCKS_PER_SM_MAKE_TREE  (1)
#endif//NBLOCKS_PER_SM_MAKE_TREE < 1
#   if  NBLOCKS_PER_SM_LINK_TREE < 1
#undef  NBLOCKS_PER_SM_LINK_TREE
#define NBLOCKS_PER_SM_LINK_TREE  (1)
#endif//NBLOCKS_PER_SM_LINK_TREE < 1
//-------------------------------------------------------------------------
#endif//MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  GADGET_MAC
__constant__      real mac_delta;/* := G / alpha */
#endif//GADGET_MAC
#ifdef  WS93_MAC
__constant__      real mac_delta;/* := sqrt(3 / accErr) */
#endif//WS93_MAC
__constant__  treecell null_cell_device;
__constant__ jparticle zero_pj;
__constant__ jmass     zero_mj;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------
typedef union
{
  int  i;
  real r;
} int_real;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------
/* initialize information on hierarchy of Peano--Hilbert key */
__global__ void initPHhierarchy_kernel(PHinfo *level)
{
  //-----------------------------------------------------------------------
  const int ii = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( ii < NUM_PHKEY_LEVEL ){
    //---------------------------------------------------------------------
    /* load current information on PH-key hierarchy */
    PHinfo phLev = level[ii];
    //---------------------------------------------------------------------
    /* initialize PH-key information */
    phLev.head = NULL_CELL;
    phLev.num  = 0;
    //---------------------------------------------------------------------
    /* set a root cell */
    if( ii == 0 ){
      //-------------------------------------------------------------------
      phLev.head = 0;
      phLev.num  = 1;
      //-------------------------------------------------------------------
    }/* if( ii == 0 ){ */
    //---------------------------------------------------------------------
    /* store initialized information on PH-key hierarchy */
    level[ii] = phLev;
    //---------------------------------------------------------------------
  }/* if( ii < NUM_PHKEY_LEVEL ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* initialize information on tree cells */
__global__ void initTreeCell_kernel
(const int cellNum, treecell * RESTRICT cell, PHint * RESTRICT hkey, uint * RESTRICT parent, uint * RESTRICT children, bool * RESTRICT leaf, uint * RESTRICT ptag, const int piNum)
{
  //-----------------------------------------------------------------------
  const int ii = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( ii < cellNum ){
    //---------------------------------------------------------------------
    /* initialize tree-cell information */
    treecell tmp_null_cell = {0, NULL_CELL};
    int tmp_hkey = -1;
    //---------------------------------------------------------------------
    /* set a root cell */
    if( ii == 0 ){
      //-------------------------------------------------------------------
      tmp_null_cell.head = 0;
      tmp_null_cell.num  = piNum;
      tmp_hkey           = 0;
      //-------------------------------------------------------------------
    }/* if( ii == 0 ){ */
    //---------------------------------------------------------------------
    /* store initialized information on tree-cell */
    cell    [ii] = tmp_null_cell;
    hkey    [ii] = tmp_hkey;
    parent  [ii] = NULL_CELL;
    children[ii] = NULL_CELL;
    leaf    [ii] = true;
    ptag    [ii] = NULL_NODE;
    //---------------------------------------------------------------------
  }/* if( ii < cellNum ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* initialize information on tree cells */
__global__ void initTreeCellOffset_kernel
(const int cellHead, const int cellNum, treecell * RESTRICT cell, PHint * RESTRICT hkey, uint * RESTRICT parent, uint * RESTRICT children, bool * RESTRICT leaf, uint * RESTRICT ptag, const int piNum)
{
  //-----------------------------------------------------------------------
  const int ii = cellHead + GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( ii < cellNum ){
    //---------------------------------------------------------------------
    /* initialize tree-cell information */
    treecell tmp_null_cell = {0, NULL_CELL};
    int tmp_hkey = -1;
    //---------------------------------------------------------------------
    /* set a root cell */
    if( ii == 0 ){
      //-------------------------------------------------------------------
      tmp_null_cell.head = 0;
      tmp_null_cell.num  = piNum;
      tmp_hkey           = 0;
      //-------------------------------------------------------------------
    }/* if( ii == 0 ){ */
    //---------------------------------------------------------------------
    /* store initialized information on tree-cell */
    cell    [ii] = tmp_null_cell;
    hkey    [ii] = tmp_hkey;
    parent  [ii] = NULL_CELL;
    children[ii] = NULL_CELL;
    leaf    [ii] = true;
    ptag    [ii] = NULL_NODE;
    //---------------------------------------------------------------------
  }/* if( ii < cellNum ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* initialize information on pseudo particles */
__global__ void initTreeNode_kernel
(const int pjNum, uint * RESTRICT more, int * RESTRICT node2cell
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
 , int * RESTRICT niSub
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
)
{
  //-----------------------------------------------------------------------
  const int jj = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( jj < pjNum ){
    //---------------------------------------------------------------------
    more     [jj] = NULL_NODE;
    node2cell[jj] = NULL_CELL;
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
    niSub    [jj] = 0;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
    //---------------------------------------------------------------------
  }/* if( jj < pjNum ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* initialize information on relation between i-particle and real j-particle */
__global__ void initTreeLink_kernel(const int piNum, int *jtag)
{
  //-----------------------------------------------------------------------
  const int ii = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( ii < piNum )
    jtag[ii] = NULL_NODE;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define NTHREADS_MAKE_INC NTHREADS_MAKE_TREE
#include "../tree/make_inc.cu"
//-------------------------------------------------------------------------
__global__ void __launch_bounds__(NTHREADS_MAKE_TREE, NBLOCKS_PER_SM_MAKE_TREE) makeTree_kernel
     (int * RESTRICT leafLev_gm, PHinfo * RESTRICT level,
      int * RESTRICT numCell_gm, treecell * RESTRICT cell, PHint * RESTRICT hkey, bool * RESTRICT leaf, uint * RESTRICT children, uint * RESTRICT parent,
      READ_ONLY PHint * RESTRICT peano,
      int * RESTRICT gmem, volatile int * RESTRICT scanNum_gm,
      int * RESTRICT gsync0Ful, int * RESTRICT gsync1Ful, int * RESTRICT gsync0Loc, int * RESTRICT gsync1Loc)
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int bnum =   GRIDDIM_X1D;
  const int tidx = THREADIDX_X1D;
  /* const int gidx = GLOBALIDX_X1D; */
  const int bidx =  BLOCKIDX_X1D;
  /* const int gnum =  BLOCKDIM_X1D * bnum; */
  const int gnum = NGROUPS_MAKE_TREE * bnum;
  const int lane = tidx & (TSUB_MAKE_TREE - 1);/* index of the thread within a thread group */
  const int scanLane = tidx & (warpSize - 1);
  //-----------------------------------------------------------------------
  const int ghead = tidx - lane;
  const int gtail = ghead + (TSUB_MAKE_TREE - 1);
  //-----------------------------------------------------------------------
  const int  idx = tidx / TSUB_MAKE_TREE;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* shared quantities in the thread parallelized version */
  //-----------------------------------------------------------------------
  __shared__ PHinfo lev_sm[NUM_PHKEY_LEVEL];
  if( tidx < NUM_PHKEY_LEVEL )
    lev_sm[tidx] = level[tidx];
  __shared__ treecell cell_sm[NTHREADS_MAKE_TREE];
  __shared__ PHint    hkey_sm[NTHREADS_MAKE_TREE];
  //-----------------------------------------------------------------------
  __shared__ int smem[NTHREADS_MAKE_TREE];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* initialize # of cells and nodes */
  //-----------------------------------------------------------------------
  /* int cellNum  = 1; */
  /* int cellRem  = NUM_ALLOC_TREE_CELL - 1; */
  /* int nodeHead = 0; */
  /* int nodeNum  = 0; */
  int numCell = 1;
  /* int numNode = 0; */
  /* int hidNode = 0;/\* head index *\/ */
  /* /\* static int phead;/\\* the head index of arrays to store pseudo j-particles (pj, mj) *\\/ *\/ */
  /* /\* phead = 0; *\/ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* make tree structure in a width-first manner */
  //-----------------------------------------------------------------------
  __syncthreads();
  for(int levelIdx = 0; levelIdx < MAXIMUM_PHKEY_LEVEL; levelIdx++){
    //---------------------------------------------------------------------
    /* set level of tree cells to examine in this procedure */
    //---------------------------------------------------------------------
    /* const PHinfo lev = level[levelIdx]; */
    const PHinfo lev = lev_sm[levelIdx];
    const int cellLev  = lev.level;
    const int cellHead = lev.head;
    int cellNum  = lev.num;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* set keylevel of Peano-Hilbert key hierarchy about child cell(s) */
    //---------------------------------------------------------------------
    /* PH key level is common within the threads */
    const   int leafLevel = cellLev - 1;
    const PHint leafScale = (PHint)1 << (leafLevel * 3);
    //---------------------------------------------------------------------
    /* PHinfo daughter = level[levelIdx + 1]; */
    PHinfo daughter = lev_sm[levelIdx + 1];
#ifdef  SPLIT_CHILD_CELL
    const   int leafNmax = daughter.nmax;
#endif//SPLIT_CHILD_CELL
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* global properties of tree cells belong to the same level hierarchy */
    //---------------------------------------------------------------------
    int cellTail;  /* equal to cellHead of tree cells belong to the lower level */
    int     tail;  /* equal to the head index of the write buffer for child cells */
    //---------------------------------------------------------------------
    tail = cellTail = cellHead + cellNum;
    //---------------------------------------------------------------------
    /* if( (tail % CELL_UNIT) != 0 ){ */
    /*   tail += (CELL_UNIT - (cellTail % CELL_UNIT)); */
    /*   cellTail = tail; */
    /*   numCell = tail; */
    /* }/\* if( (tail % CELL_UNIT) != 0 ){ *\/ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* examine tree cells */
    //---------------------------------------------------------------------
    const int Niter = BLOCKSIZE(cellNum, gnum);
    for(int iter = 0; iter < Niter; iter++){
      //-------------------------------------------------------------------
      const int cnumSub = (gnum < cellNum) ? gnum : cellNum;
      /* const int bnumSub = BLOCKSIZE(cnumSub, NTHREADS_MAKE_TREE);/\* # of active blocks (required to global synchronization) *\/ */
      const int bnumSub = BLOCKSIZE(cnumSub, NGROUPS_MAKE_TREE);/* # of active blocks (required to global synchronization) */
      //-------------------------------------------------------------------
      if( bidx < bnumSub ){
	//-----------------------------------------------------------------
	/* load a tree cell and evaluate either the cell has child cell(s) or no */
	//-----------------------------------------------------------------
	/* coalesced load from global memory to shared memory */
	int cidx = cellHead + gnum * iter + bidx * NGROUPS_MAKE_TREE + tidx;
	if( (tidx < NGROUPS_MAKE_TREE) && (cidx < cellTail) ){
	  //---------------------------------------------------------------
	  cell_sm[tidx] = cell[cidx];
	  hkey_sm[tidx] = hkey[cidx];
	  //---------------------------------------------------------------
	}/* if( (tidx < NGROUPS_MAKE_TREE) && (cidx < cellTail) ){ */
	__syncthreads();
	cidx = cellHead + gnum * iter + bidx * NGROUPS_MAKE_TREE + idx;
	const treecell root = (cidx < cellTail) ? (cell_sm[idx]) : (null_cell_device);
	const PHint keyHead = (cidx < cellTail) ? (hkey_sm[idx]) : ((PHint)(-1));
	//-----------------------------------------------------------------
	const bool node = (root.num > NCRIT) ? (true) : (false);
	//-----------------------------------------------------------------


	//-----------------------------------------------------------------
	/* divide the responsible tree cell if the cell is not a leaf cell */
	//-----------------------------------------------------------------
	__syncthreads();
	/* cell_sm[tidx] = root; */
	cell_sm[tidx].head = 0;
	/* cell_sm[tidx].num  = 0; */
	/* if( lane == TSUB_MAKE_TREE - 1 ) */
	  cell_sm[tidx].num = root.num;
	int addChild = 0;
	//-----------------------------------------------------------------
	if( node ){
	  //---------------------------------------------------------------
	  /* the responsible tree cell is a node cell */
	  if( lane == 0 )
	    leaf[cidx] = false;
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* find CELL_UNIT - 1 (== TSUB_MAKE_TREE - 1 == 7) boundaries (in maximum) */
	  //---------------------------------------------------------------
	  if( lane < (CELL_UNIT - 1) ){
	    //-------------------------------------------------------------
	    /* properties of the PH-key of the focused child cell */
	    /* NOTE: the child cell may not exist */
	    //-------------------------------------------------------------
	    /* const PHint khead = keyHead +      lane  * leafScale; */
	    const PHint ktail = keyHead + (1 + lane) * leafScale;
	    //-------------------------------------------------------------
	    int ihead = 0;
	    int itail = root.num - 1;
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    /* check whether the specified region contains keys between khead and ktail */
	    //-------------------------------------------------------------
	    /* if( lane < 2 ) */
	    /*   hkey_sm[tidx] = peano[root.head + lane * itail];/\* lane = 0 or 1 *\/ */
	    /* if( (hkey_sm[ghead] < ktail) && (ktail <= hkey_sm[ghead + 1]) ){ */
	    if( lane == 0 )
	      hkey_sm[tidx] = peano[root.head + itail];
	    if( ktail <= hkey_sm[ghead] ){
	      //-----------------------------------------------------------
#ifdef  USE_BRENT_SCHEME_FOR_TREE_BUILD
	      /* //----------------------------------------------------------- */
	      /* ihead = rootFindBrentIdx(ihead, itail, &peano[root.head], ktail, UNITY); */
	      /* //----------------------------------------------------------- */
#else///USE_BRENT_SCHEME_FOR_TREE_BUILD
	      //-----------------------------------------------------------
	      /* find the tail of the PH-key for the child cell using bisection method */
	      //-----------------------------------------------------------
	      while( true ){
		//---------------------------------------------------------
		/* when the procedure finds the tail of the PH-key... */
		if( itail <= (1 + ihead) )		  break;
		//---------------------------------------------------------
		uint ihalf = (uint)(ihead + itail) >> 1;
		if( peano[root.head + ihalf] <= ktail )		  ihead = (int)ihalf;
		else		                                  itail = (int)ihalf;
		//---------------------------------------------------------
	      }/* while( true ){ */
	      //-----------------------------------------------------------
#endif//USE_BRENT_SCHEME_FOR_TREE_BUILD
	      //-----------------------------------------------------------
	      itail = ihead;
	      //-----------------------------------------------------------
	    }/* if( ktail <= hkey_sm[ghead] ){ */
	    //-------------------------------------------------------------
	    if( ((peano[root.head] - keyHead) / leafScale) <= lane )
	      itail++;
	    cell_sm[tidx].num = itail;
	    //-------------------------------------------------------------
	  }/* if( lane < (CELL_UNIT - 1) ){ */
	  //---------------------------------------------------------------
#if 0
	  __syncthreads();
	  if( tidx == 0 ){
	    for(int ll = 0; ll < cnumSub * TSUB_MAKE_TREE; ll++)
	      printf("cellLev = %d, cidx = %d, head = %d, itail = %d\n", cellLev, cellHead + gnum * iter + (ll / TSUB_MAKE_TREE), cell_sm[ll / TSUB_MAKE_TREE].head, cell_sm[ll].num);
	  }
	  __syncthreads();
#endif
	  //---------------------------------------------------------------
	  /* assume implicit synchronization within a warp */
	  if( lane > 0 )
	    cell_sm[tidx].head = cell_sm[tidx - 1].num;
	  cell_sm[tidx].num -= cell_sm[tidx].head;

	  /* if( lane > 0 ) */
	  /*   cell_sm[tidx].head = 1 + cell_sm[tidx - 1].num; */
	  /* cell_sm[tidx].num -= cell_sm[tidx].head; */

	  /* __syncthreads(); */
	  /* cell_sm[tidx].num -= cell_sm[tidx].head; */

	  /* if( lane == (CELL_UNIT - 1) ) */
	  /*   cell_sm[tidx].num = root.num - cell_sm[tidx].head; */

	  /* if( lane > 0 ){ */
	  /* /\* if( lane > (int)((hkey_sm[ghead] - keyHead) / leafScale) ){ *\/ */
	  /*   cell_sm[tidx].head = cell_sm[tidx - 1].num; */
	  /*   cell_sm[tidx].num -= cell_sm[tidx    ].head; */
	  /* }/\* if( lane > 0 ){ *\/ */
	  /* /\* else *\/ */
	  /* /\*   cell_sm[tidx].head = 0; *\/ */
	  cell_sm[tidx].head += root.head;
#if 0
	  __syncthreads();
	  if( tidx == 0 ){
	    for(int cc = 0; cc < cnumSub; cc++){
	      printf("cidx = %d: (head, num) = (%d, %d)", cellHead + cc, cell_sm[tidx + cc * CELL_UNIT].head, cell_sm[tidx + cc * CELL_UNIT].num);
	      for(int ll = 1; ll < CELL_UNIT; ll++)
		printf(", (%d, %d)", cell_sm[tidx + ll + cc * CELL_UNIT].head, cell_sm[tidx + ll + cc * CELL_UNIT].num);
	      printf("\n");
	    }
	  }
	  __syncthreads();
#endif
	  //---------------------------------------------------------------
#ifdef  SPLIT_CHILD_CELL
	  /* /\* split child cell if the number of N-body particles exceeds the upper limit *\/ */
	  /* const int nsplit = (itail < leafNmax) ? (1) : (((itail) + (leafNmax - 1)) / (leafNmax)); */
	  /* /\* __NOTE__("jj = %d, lone = %d, cellLev = %d, leafNmax = %d, cidx = %d, itail = %d, nsplit = %d\n", jj, lone, cellLev, leafNmax, cidx, itail, nsplit); *\/ */
	  /* /\* const uint nchild = itail / nsplit; *\/ */
	  /* const int nchild = (itail + (nsplit - 1)) / nsplit; */

	  /* split child cell if the number of N-body particles exceeds the upper limit */
	  /* /\* addChild = BLOCKSIZE(cell_sm[tidx].num, leafNmax); *\/ */
	  /* addChild = (cell_sm[tidx].num < leafNmax) ? (1) : BLOCKSIZE(cell_sm[tidx].num, leafNmax); */
	  /* addChild = (cell_sm[tidx].num < leafNmax) ? ((cell_sm[tidx].num > 0) ? 1 : 0) : BLOCKSIZE(cell_sm[tidx].num, leafNmax); */
	  addChild = (cell_sm[tidx].num > 0) ? BLOCKSIZE(cell_sm[tidx].num, leafNmax) : 0;
#else///SPLIT_CHILD_CELL
	  if( cell_sm[tidx].num > 0 )
	    addChild = 1;
#endif//SPLIT_CHILD_CELL
	  //---------------------------------------------------------------
	}/* if( node ){ */
	//-----------------------------------------------------------------
#if 0
	if( addChild > 0 )
	  printf("levelIdx = %d: cidx = %d: bidx = %d: tidx = %d: addChild = %d, cell_sm.num = %d, leafNmax = %d\n", levelIdx, cidx, bidx, tidx, addChild, cell_sm[tidx].num, leafNmax);
#endif
	//-----------------------------------------------------------------
	/* cellLev = 21, cidx = 0, tag = 8 */
	/* cellLev = 20, cidx = 1, tag = 8 */
	/* cellLev = 20, cidx = 2, tag = 6 */
	/* cellLev = 20, cidx = 3, tag = 3 */
	/* cellLev = 20, cidx = 4, tag = 5 */
	/* cellLev = 20, cidx = 5, tag = 4 */
	/* cellLev = 20, cidx = 6, tag = 3 */
	/* cellLev = 20, cidx = 7, tag = 4 */
	/* cellLev = 20, cidx = 8, tag = 6 */
	//-----------------------------------------------------------------


	//-----------------------------------------------------------------
	/* calculate prefix sum about addChild */
	//-----------------------------------------------------------------
	/* 1. scan within a block */
	PREFIX_SUM_MAKE_INC_BLCK(addChild, tidx, scanLane, smem);
	/* 2. save result of local prefix sum */
	const int scanIdx = smem[tidx] - addChild;/* this must be an exclusive scan */
	int scanNum = smem[NTHREADS_MAKE_TREE - 1];
#if 0
	if( tidx == 0 ){
	  printf("levelIdx = %d: bidx = %d: smem = %d", levelIdx, bidx, smem[0]);
	  for(int ii = 1; ii < NTHREADS_MAKE_TREE; ii++)
	    printf(", %d", smem[ii]);
	  printf("\n");
	}
#endif
	int headIdx = 0;
	const int tag = smem[gtail] - ((gtail >= TSUB_MAKE_TREE) ? smem[gtail - TSUB_MAKE_TREE] : 0);
	/* 3. scan within a grid */
	PREFIX_SUM_MAKE_INC_GRID(&headIdx, &scanNum, bnumSub, bidx, tidx, scanLane, smem, gmem, gsync0Loc, gsync1Loc);
	/* 4. the thread (gidx = 0) store the total number of child cells added in this step */
	if( bidx + tidx == 0 )
	  *scanNum_gm = scanNum;
#if 0
	if( tidx == 0 )
	  printf("levelIdx = %d: cellNum = %d, bnum = %d, gnum = %d: iter = %d/%d: bidx = %d/%d: scanNum(ful) = %d, scanNum(loc) = %d\n", levelIdx, cellNum, bnum, gnum, iter, Niter, bidx, bnumSub, scanNum, smem[NTHREADS_MAKE_TREE - 1]);
#endif
	//-----------------------------------------------------------------
/* #if 1 */
/* 	if( tag > NLEAF ){ */
/* 	  __KILL__(stderr, "ERROR: # of child cells (%d) exceeds NLEAF(%d). Enlarge NLEAF and/or MAXIMUM_PHKEY_LEVEL(%d)\n", tag, NLEAF, MAXIMUM_PHKEY_LEVEL); */
/* 	} */
/* #endif */
	//-----------------------------------------------------------------


	//-----------------------------------------------------------------
	/* store child cells to the global memory */
	/* NOTE: tentative implementation (uncoalesced store to global memory) */
	/* if performance is too low due to uncoalesced store, then try to implement coalesced version (probably, use additional shared memory to stock splitted child cells) */
	//-----------------------------------------------------------------
	int targetIdx = tail + headIdx + scanIdx;
	/* if( lane == 0 ) */
	if( (node) && (lane == 0) )
	  children[cidx] = ((tag - 1) << IDXBITS) + targetIdx;
	//-----------------------------------------------------------------
#ifdef  SPLIT_CHILD_CELL
	//-----------------------------------------------------------------
	treecell childCell = cell_sm[tidx];
	int nrem = childCell.num;
	const int nchild = BLOCKSIZE(nrem, addChild);
	/* const int nchild = (itail + (addChild - 1)) / addChild; */
	/* int nrem = itail; */
	for(int childIdx = 0; childIdx < addChild; childIdx++){
	  //---------------------------------------------------------------
	  /* estimate the number of particles per child cell */
	  //---------------------------------------------------------------
	  childCell.num = (nchild < nrem) ? nchild : nrem;
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* commit the new child cell */
	  //---------------------------------------------------------------
	  cell  [targetIdx + childIdx] = childCell;
	  hkey  [targetIdx + childIdx] = keyHead + lane * leafScale;
	  parent[targetIdx + childIdx] = cidx;
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* preparation to the next child cell */
	  //---------------------------------------------------------------
	  childCell.head += childCell.num;
	  //---------------------------------------------------------------
	}/* for(int childIdx = 0; childIdx < addChild; childIdx++){ */
	//-----------------------------------------------------------------
#else///SPLIT_CHILD_CELL
	//-----------------------------------------------------------------
	/* commit the new child cell */
	//-----------------------------------------------------------------
	if( addChild > 0 ){
	  //---------------------------------------------------------------
	  cell  [targetIdx] = cell_sm[tidx];
	  hkey  [targetIdx] = keyHead + lane * leafScale;
	  parent[targetIdx] = cidx;
	  //---------------------------------------------------------------
	}/* if( addChild > 0 ){ */
	//-----------------------------------------------------------------
#endif//SPLIT_CHILD_CELL
	//-----------------------------------------------------------------
/* 	/\* error check *\/ */
/* 	if( lone != 0 ){ */
/* 	  __NOTE__("the cell has %d child-cells\n", tag); */
/* #   if  MAXIMUM_PHKEY_LEVEL <= 10 */
/* 	  __NOTE__("head of the key of the first child cell is %u\n", keyHead); */
/* 	  __NOTE__("tail of the key of the  last child cell is %u\n", ktail); */
/* 	  __NOTE__("key scale of child cells is %u\n", leafScale); */
/* #else///MAXIMUM_PHKEY_LEVEL <= 10 */
/* 	  __NOTE__("head of the key of the first child cell is %zu\n", keyHead); */
/* 	  __NOTE__("tail of the key of the  last child cell is %zu\n", ktail); */
/* 	  __NOTE__("key scale of child cells is %zu\n", leafScale); */
/* #endif//MAXIMUM_PHKEY_LEVEL <= 10 */
/* 	  __KILL__(stderr, "ERROR: tree construction failed. # of left particle(s) is %d(/%d).\n", lone, root.num); */
/* 	}/\* if( lone != 0 ){ *\/ */
	//-----------------------------------------------------------------
      }/* if( bidx < bnumSub ){ */
      //-------------------------------------------------------------------
      /* global synchronization within bnum blocks */
      globalSync(tidx, bidx, bnum, gsync0Ful, gsync1Ful);
      //-------------------------------------------------------------------
      if( tidx == 0 )
	smem[0] = *scanNum_gm;
      __syncthreads();
      /* const int increment = BLOCKSIZE(scanNum, CELL_UNIT) * CELL_UNIT; */
      /* tail += increment; */
      /* numCell += increment; */
      const  int addCellNum = smem[0];
      tail    += addCellNum;
      numCell += addCellNum;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
    daughter.head =        cellTail;
    daughter.num  = tail - cellTail;
    if( tidx == 0 )
      lev_sm[levelIdx + 1] = daughter;
    //---------------------------------------------------------------------
#if 0
    if( (bidx + tidx) == 0 )
      for(int ii = 0; ii < MAXIMUM_PHKEY_LEVEL; ii++)
	printf("levelIdx = %d: ii = %d: head = %d, num = %d\n", levelIdx, ii, lev_sm[ii].head, lev_sm[ii].num);
#endif
    //---------------------------------------------------------------------
    if( daughter.num == 0 )
      break;
    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0Ful, gsync1Ful);
    //---------------------------------------------------------------------
  }/* for(int levelIdx = 0; levelIdx < MAXIMUM_PHKEY_LEVEL; levelIdx++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  if( bidx == 0 ){
    //---------------------------------------------------------------------
    if( tidx < NUM_PHKEY_LEVEL )
      level[tidx] = lev_sm[tidx];
    //---------------------------------------------------------------------
    if( tidx == NUM_PHKEY_LEVEL )
      *numCell_gm = numCell;
    //---------------------------------------------------------------------
    if( tidx == 0 ){
      //-------------------------------------------------------------------
      int leafLev = NUM_PHKEY_LEVEL - 1;
      for(int levelIdx = (NUM_PHKEY_LEVEL - 1); levelIdx >= 0; levelIdx--)
	if( lev_sm[levelIdx].num != 0 ){
	  leafLev = levelIdx;
	  break;
	}/* if( lev_sm[levelIdx].num != 0 ){ */
      *leafLev_gm = leafLev;
      //-------------------------------------------------------------------
    }/* if( tidx == 0 ){ */
    //---------------------------------------------------------------------
  }/* if( bidx == 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* extend the pseudo particle chain */
//-------------------------------------------------------------------------
#   if  NTHREADS_MAKE_INC != NTHREADS_LINK_TREE
#undef  NTHREADS_MAKE_INC
#define NTHREADS_MAKE_INC    NTHREADS_LINK_TREE
#include "../tree/make_inc.cu"
#endif//NTHREADS_MAKE_INC != NTHREADS_LINK_TREE
//-------------------------------------------------------------------------
__device__ __forceinline__ void linkNode
(const treecell root, const bool leaf_cell, uint * RESTRICT ptag,
 uint * RESTRICT more, int * RESTRICT jtag, int * RESTRICT phead,
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
 int * RESTRICT niSub,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
 const int bnum, const int bidx, const int tidx, const int lane,
 volatile int * RESTRICT smem, int * RESTRICT gmem, int * RESTRICT gsync0, int * RESTRICT gsync1
)
{
  //-----------------------------------------------------------------------
  /* load fundamental information on the focusing tree cell */
  //-----------------------------------------------------------------------
  const int head = root.head;
  /* below procedure is switched on if head != NULL_CELL */
  /* if the cell is a leaf cell, then root.num pseudo particles are set. otherwise, a pseudo particle is set.  */
  const int nadd = (head != NULL_CELL) ? ((leaf_cell) ? (root.num) : (1)) : (0);
  /* if( nadd > NLEAF ){ */
  /*   __KILL__(stderr, "ERROR: nadd (%d) exceeds NLEAF(%d). Enlarge NLEAF and/or MAXIMUM_PHKEY_LEVEL(%d)\n", nadd, NLEAF, MAXIMUM_PHKEY_LEVEL); */
  /* }/\* if( nadd > NLEAF ){ *\/ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate prefix sum about addChild */
  //-----------------------------------------------------------------------
  /* 1. scan within a block */
  PREFIX_SUM_MAKE_INC_BLCK(nadd, tidx, lane, smem);
  /* 2. save result of local prefix sum */
  const int scanIdx = smem[tidx] - nadd;/* this must be an exclusive scan */
  int scanNum = smem[NTHREADS_LINK_TREE - 1];
  int headIdx = 0;
  /* 3. scan within a grid */
  PREFIX_SUM_MAKE_INC_GRID(&headIdx, &scanNum, bnum, bidx, tidx, lane, smem, gmem, gsync0, gsync1);
  /* /\* 4. the thread (gidx = 0) store the total number of child cells added in this step *\/ */
  /* if( bidx + tidx == 0 ) */
  /*   *scanNum_gm = scanNum; */
  headIdx += scanIdx + (*phead);
  if( nadd > 0 )    *ptag = ((nadd - 1) << IDXBITS) + headIdx;
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  else    niSub[*phead] = root.num;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* when the tree cell is a leaf cell, then more index specifies the cell itself */
  //-----------------------------------------------------------------------
  if( (head != NULL_CELL) && leaf_cell ){
    //---------------------------------------------------------------------
    for(int jj = 0; jj < nadd; jj++){
      //-------------------------------------------------------------------
      /* commit the new particle to the j-particle array */
      const int jidx = headIdx + jj;
      more[jidx] = jidx;
      //-------------------------------------------------------------------
      /* connect an i-particle with the corresponding real j-particle */
      jtag[head + jj] = jidx;
      //-------------------------------------------------------------------
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
      niSub[jidx] = 1;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < nadd; jj++){ */
    //---------------------------------------------------------------------
  }/* if( (head != NULL_CELL) && lear_cell ){ */
  //-----------------------------------------------------------------------
  *phead += scanNum;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__global__ void __launch_bounds__(NTHREADS_LINK_TREE, NBLOCKS_PER_SM_LINK_TREE) linkTree_kernel
     (const int cellTail, READ_ONLY treecell * RESTRICT cell, READ_ONLY bool * RESTRICT leaf, uint * RESTRICT ptag,
      int * RESTRICT pjNum, uint * RESTRICT more, int * RESTRICT jtag,
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
      int * RESTRICT niSub,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
      int * RESTRICT gmem, int * RESTRICT gsync0, int * RESTRICT gsync1
      )
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int bnum =   GRIDDIM_X1D;
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int lane = tidx & (warpSize - 1);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* shared quantities in the thread parallelized version */
  //-----------------------------------------------------------------------
  __shared__ int smem[NTHREADS_LINK_TREE];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* extend the pseudo particle chain */
  /* examine tree cells */
  //-----------------------------------------------------------------------
  int cellNum = cellTail;
  int numNode = 0;
  //-----------------------------------------------------------------------
  const int chkNumStep = NTHREADS_LINK_TREE * bnum;
  const int Nloop = BLOCKSIZE(cellNum, chkNumStep);
  for(int loop = 0; loop < Nloop; loop++){
    //---------------------------------------------------------------------
    const int chkNum = (chkNumStep < cellNum) ? chkNumStep : cellNum;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* load a tree cell and evaluate either the cell has child cell(s) or no */
    //---------------------------------------------------------------------
    const int cidx = gidx + chkNumStep * loop;
    treecell root = (cidx < cellTail) ? (cell[cidx]) : (null_cell_device);
    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* extend the pseudo particle chain */
    //---------------------------------------------------------------------
    linkNode(root, leaf[cidx], &ptag[cidx], more, jtag, &numNode,
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	     niSub,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	     bnum, bidx, tidx, lane, smem, gmem, gsync0, gsync1);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    cellNum -= chkNum;
    //---------------------------------------------------------------------
  }/* for(int loop = 0; loop < Nloop; loop++){ */
  //-----------------------------------------------------------------------
  /* return # of pseudo j-particles */
  if( gidx == 0 )
    *pjNum = numNode;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* extend the pseudo particle chain */
//-------------------------------------------------------------------------
__global__ void trimTree_kernel
(const int cellHead, const int cellTail, READ_ONLY treecell * RESTRICT cell, READ_ONLY bool * RESTRICT leaf, READ_ONLY uint * RESTRICT children, READ_ONLY uint * RESTRICT ptag,
 int * RESTRICT node2cell, uint * RESTRICT more)
{
  //-----------------------------------------------------------------------
  const int cidx = cellHead + GLOBALIDX_X1D;
  /* int cellNum = cellTail; */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* connect tree nodes in a different PH-level */
  //-----------------------------------------------------------------------
  if( (cidx < cellTail) && (cell[cidx].head != NULL_CELL) ){
    //---------------------------------------------------------------------
    /* relate the tree cell and corresponding tree nodes */
    const uint node = ptag[cidx];
    const int pidx = node & IDXMASK;
#ifdef  SPLIT_CHILD_CELL
    const int pnum = 1 + (int)(node >> IDXBITS);
    for(int jj = 0; jj < pnum; jj++)
      node2cell[pidx + jj] = cidx;
#else///SPLIT_CHILD_CELL
    node2cell[pidx] = cidx;
#endif//SPLIT_CHILD_CELL
    //---------------------------------------------------------------------
    if( leaf[cidx] == false ){
      //-------------------------------------------------------------------
      /* count up number of child nodes */
      //-------------------------------------------------------------------
      const uint child = children[cidx];
      const  int chead =               child  & IDXMASK;
      const  int ctail = chead + (1 + (child >> IDXBITS));
      const uint phead = ptag[chead];
      int childNum = 1 + (phead >> IDXBITS);
      for(int jj = chead + 1; jj < ctail; jj++)
	childNum += (1 + (ptag[jj] >> IDXBITS));
      //-------------------------------------------------------------------
/* #if 1 */
/*       if( childNum > NLEAF ){ */
/* 	__KILL__(stderr, "ERROR: childNum(%d) exceeds NLEAF(%d) for cidx = %d @ levelIdx = %d. Enlarge NLEAF and/or MAXIMUM_PHKEY_LEVEL(%d)\n", childNum, NLEAF, cidx, levelIdx, MAXIMUM_PHKEY_LEVEL); */
/*       } */
/* #endif */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* commit child particles as pseudo particle */
      //-------------------------------------------------------------------
      more[pidx] = ((childNum - 1) << IDXBITS) + (phead & IDXMASK);
      //-------------------------------------------------------------------
    }/* if( leaf[cidx] == false ){ */
    //---------------------------------------------------------------------
  }/* if( (cidx < cellTail) && (cell[cidx].head != NULL_CELL) ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
extern "C"
void makeTreeStructure_dev
(const int piNum, PHint * RESTRICT peano,
 int * RESTRICT leafLev, int * RESTRICT leafLev_dev, int * RESTRICT scanNum_dev,
 int * RESTRICT cellNum, int * RESTRICT cellNum_dev, const soaTreeCell cell,
 int * RESTRICT nodeNum, int * RESTRICT nodeNum_dev, const soaTreeNode node,
 const soaMakeTreeBuf buf, deviceProp devProp
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* start stop watch */
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* initialize tree structure */
  //-----------------------------------------------------------------------
  int Nrem = BLOCKSIZE(piNum, NTHREADS_INIT_LINK);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    initTreeLink_kernel<<<Nrem, NTHREADS_INIT_LINK>>>(piNum, node.jtag);
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
      int Nsub = Nblck * NTHREADS_INIT_LINK;
      initTreeLink_kernel<<<Nblck, NTHREADS_INIT_LINK>>>(Nsub, &node.jtag[hidx]);
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  /* initTreeLink_kernel<<<BLOCKSIZE(   piNum, 1024), 1024>>>(   piNum, node.jtag); */
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->initTreeLink_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------
  Nrem = BLOCKSIZE(*cellNum, NTHREADS_INIT_CELL);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    initTreeCell_kernel<<<Nrem, NTHREADS_INIT_CELL>>>(*cellNum, cell.cell, cell.hkey, cell.parent, cell.children, cell.leaf, cell.ptag, piNum);
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
      int Nsub = Nblck * NTHREADS_INIT_CELL;
      initTreeCellOffset_kernel<<<Nblck, NTHREADS_INIT_CELL>>>(hidx, *cellNum, cell.cell, cell.hkey, cell.parent, cell.children, cell.leaf, cell.ptag, piNum);
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  /* initTreeCell_kernel<<<BLOCKSIZE(*cellNum, 1024), 1024>>>(*cellNum, cell.cell, cell.hkey, cell.parent, cell.children, cell.leaf, cell.ptag, piNum); */
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->initTreeCell_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------
  Nrem = BLOCKSIZE(*nodeNum, NTHREADS_INIT_NODE);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    initTreeNode_kernel<<<Nrem, NTHREADS_INIT_NODE>>>
      (*nodeNum, node.more, node.node2cell
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
       , node.niSub
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
       );
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
      int Nsub = Nblck * NTHREADS_INIT_NODE;
      initTreeNode_kernel<<<Nblck, NTHREADS_INIT_NODE>>>
	(Nsub, &node.more[hidx], &node.node2cell[hidx]
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	 , &node.niSub[hidx]
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	 );
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
/*   initTreeNode_kernel<<<BLOCKSIZE(*nodeNum, 1024), 1024>>>(*nodeNum, node.more, node.node2cell */
/* #   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH) */
/* 							   , node.niSub */
/* #endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH) */
/* 							   ); */
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->initTreeNode_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------
  initPHhierarchy_kernel<<<1, NUM_PHKEY_LEVEL>>>(cell.level);
  //-----------------------------------------------------------------------
  /* *cellNum = 1; */
  /* *cellRem = NUM_ALLOC_TREE_CELL - 1; */
  /* static int phead;/\* the head index of arrays to store pseudo j-particles (pj, mj) *\/ */
  /* phead = 0; */
  /* *nodeNum = 0; */
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->makeTree));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* make tree structure */
  //-----------------------------------------------------------------------
  makeTree_kernel<<<devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE, NTHREADS_MAKE_TREE>>>
    (leafLev_dev, cell.level, cellNum_dev, cell.cell, cell.hkey, cell.leaf, cell.children, cell.parent, peano,
     buf.gmem_make_tree, scanNum_dev, buf.gsync0_make_tree, buf.gsync1_make_tree, buf.gsync2_make_tree, buf.gsync3_make_tree);
  getLastCudaError("makeTree_kernel");
  checkCudaErrors(cudaMemcpy(leafLev, leafLev_dev, sizeof(int), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cellNum, cellNum_dev, sizeof(int), cudaMemcpyDeviceToHost));
  if( *cellNum > NUM_ALLOC_TREE_CELL ){
    __KILL__(stderr, "ERROR: allocated vector length for tree-cell arrays are too short, increase TREE_SAFETY_VAL(%f) defined in src/tree/make.h or use more MPI processes!\n", TREE_SAFETY_VAL);
  }/* if( *cellNum > NUM_ALLOC_TREE_CELL ){ */
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->makeTree_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* link tree structure */
  //-----------------------------------------------------------------------
  linkTree_kernel<<<devProp.numSM * NBLOCKS_PER_SM_LINK_TREE, NTHREADS_LINK_TREE>>>
    (*cellNum, cell.cell, cell.leaf, cell.ptag,
     nodeNum_dev, node.more, node.jtag,
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
     node.niSub,
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
     buf.gmem_link_tree, buf.gsync0_link_tree, buf.gsync1_link_tree
     );
  getLastCudaError("linkTree_kernel");
  checkCudaErrors(cudaMemcpy(nodeNum, nodeNum_dev, sizeof(int), cudaMemcpyDeviceToHost));
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->linkTree_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* trim tree structure */
  //-----------------------------------------------------------------------
  Nrem = BLOCKSIZE(*cellNum, NTHREADS_TRIM_TREE);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    trimTree_kernel<<<Nrem, NTHREADS_TRIM_TREE>>>(0, *cellNum, cell.cell, cell.leaf, cell.children, cell.ptag, node.node2cell, node.more);
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
      int Nsub = Nblck * NTHREADS_TRIM_TREE;
      /* trimTree_kernel<<<Nblck, NTHREADS_TRIM_TREE>>>(hidx, Nsub, cell.cell, cell.leaf, cell.children, cell.ptag, node.node2cell, node.more); */
      trimTree_kernel<<<Nblck, NTHREADS_TRIM_TREE>>>(hidx, *cellNum, cell.cell, cell.leaf, cell.children, cell.ptag, node.node2cell, node.more);
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  /* trimTree_kernel<<<BLOCKSIZE(*cellNum, NTHREADS_TRIM_TREE), NTHREADS_TRIM_TREE>>> */
  /*   (*cellNum, cell.cell, cell.leaf, cell.children, cell.ptag, node.node2cell, node.more); */
  getLastCudaError("trimTree_kernel");
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->trimTree_kernel));
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* stop stop watch */
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
#ifndef HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->makeTree));
#else///HUNT_MAKE_PARAMETER
  elapsed->makeTree += elapsed->initTreeLink_kernel + elapsed->initTreeCell_kernel + elapsed->initTreeNode_kernel + elapsed->makeTree_kernel + elapsed->linkTree_kernel + elapsed->trimTree_kernel;
#endif//HUNT_MAKE_PARAMETER
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* arrays to store properties of tree nodes (tree cells and N-body particles) */
//-------------------------------------------------------------------------
extern "C"
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
 int **more0Buf, int **more1Buf, real **rjmaxBuf, int **fail_dev, soaTreeNode *dev, soaTreeNode *hst, soaMakeTreeBuf *buf)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* const size_t num = NUM_ALLOC_TREE_NODE; */
  const size_t num = (size_t)ceilf(EXTEND_NUM_TREE_NODE * (float)NUM_ALLOC_TREE_NODE);
  //-----------------------------------------------------------------------
  mycudaMalloc    ((void **) more_dev, num * sizeof(uint     ));  alloc.device += num * sizeof(uint     );  dev->more      = * more_dev;
  mycudaMalloc    ((void **)   pj_dev, num * sizeof(jparticle));  alloc.device += num * sizeof(jparticle);  dev->jpos      = *   pj_dev;
  mycudaMalloc    ((void **)   mj_dev, num * sizeof(jmass    ));  alloc.device += num * sizeof(jmass    );  dev->mj        = *   mj_dev;
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  mycudaMalloc    ((void **)  n2c_dev, num * sizeof(int      ));  alloc.device += num * sizeof(int      );  dev->node2cell = *  n2c_dev;
#endif//CALC_MULTIPOLE_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  mycudaMalloc    ((void **)niSub_dev, num * sizeof(int      ));  alloc.device += num * sizeof(int      );  dev->niSub     = *niSub_dev;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  //-----------------------------------------------------------------------
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
  mycudaMallocHost((void **) more_hst, num * sizeof(uint     ));  alloc.host += num * sizeof(uint     );  hst->more      = * more_hst;
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
  mycudaMallocHost((void **)   pj_hst, num * sizeof(jparticle));  alloc.host += num * sizeof(jparticle);  hst->jpos      = *   pj_hst;
  mycudaMallocHost((void **)   mj_hst, num * sizeof(jmass    ));  alloc.host += num * sizeof(jmass    );  hst->mj        = *   mj_hst;
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
#ifndef MAKE_TREE_ON_DEVICE
  mycudaMallocHost((void **)  n2c_hst, num * sizeof(int      ));  alloc.host += num * sizeof(int      );  hst->node2cell = *  n2c_hst;
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  mycudaMallocHost((void **)niSub_hst, num * sizeof(int      ));  alloc.host += num * sizeof(int      );  hst->niSub     = *niSub_hst;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifdef  GADGET_MAC
  mycudaMalloc((void **)mac_dev, num * sizeof(real));  alloc.device += num * sizeof(real);  dev->mac = *mac_dev;
#endif//GADGET_MAC
  //-----------------------------------------------------------------------
#ifdef  CALC_MULTIPOLE_ON_DEVICE
#       ifdef  WS93_MAC
  mycudaMalloc    ((void **) mr2_dev, num * sizeof(real));  alloc.device += num * sizeof(real);  dev->mr2  = * mr2_dev;
#       endif//WS93_MAC
  mycudaMalloc    ((void **)bmax_dev, num * sizeof(real));  alloc.device += num * sizeof(real);  dev->bmax = *bmax_dev;
#else///CALC_MULTIPOLE_ON_DEVICE
#       ifdef  WS93_MAC
  mycudaMallocHost((void **) mr2_hst, num * sizeof(real));  alloc.host   += num * sizeof(real);  hst->mr2  = * mr2_hst;
#       endif//WS93_MAC
  mycudaMallocHost((void **)bmax_hst, num * sizeof(real));  alloc.host   += num * sizeof(real);  hst->bmax = *bmax_hst;
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  mycudaMalloc    ((void **)more0Buf, devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(int ));
  alloc.device +=                     devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(int ) ;
  mycudaMalloc    ((void **)more1Buf, devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(int ));
  alloc.device +=                     devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(int ) ;
  mycudaMalloc    ((void **)rjmaxBuf, devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(real));
  alloc.device +=                     devProp.numSM * NBLOCKS_PER_SM_MAC * NGROUPS_MAC * NUM_ALLOC_MACBUF * sizeof(real) ;
#else///CALC_MULTIPOLE_ON_DEVICE
  mycudaMallocHost((void **)more0Buf,                                                    NUM_ALLOC_MACBUF * sizeof(int ));
  alloc.host +=                                                                          NUM_ALLOC_MACBUF * sizeof(int ) ;
  mycudaMallocHost((void **)more1Buf,                                                    NUM_ALLOC_MACBUF * sizeof(int ));
  alloc.host +=                                                                          NUM_ALLOC_MACBUF * sizeof(int ) ;
  mycudaMallocHost((void **)rjmaxBuf,                                                    NUM_ALLOC_MACBUF * sizeof(real));
  alloc.host +=                                                                          NUM_ALLOC_MACBUF * sizeof(real) ;
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  mycudaMalloc((void **)gsync0, devProp.numSM * NBLOCKS_PER_SM_MAC * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_MAC * sizeof(int);
  mycudaMalloc((void **)gsync1, devProp.numSM * NBLOCKS_PER_SM_MAC * sizeof(int));  alloc.device += devProp.numSM * NBLOCKS_PER_SM_MAC * sizeof(int);
  initGsync_kernel<<<1, devProp.numSM * NBLOCKS_PER_SM_MAC>>>(devProp.numSM * NBLOCKS_PER_SM_MAC, *gsync0, *gsync1);
  getLastCudaError("initGsync_kernel");
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
  mycudaMalloc((void **)fail_dev, 1 * sizeof(int));
  alloc.device +=                 1 * sizeof(int);
  const int fail_hst = 0;
  checkCudaErrors(cudaMemcpy(*fail_dev, &fail_hst, sizeof(int), cudaMemcpyHostToDevice));
  //-----------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
  mycudaMalloc((void **)  gmem_make_tree, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int));
  alloc.device                         += devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int);
  mycudaMalloc((void **)gsync0_make_tree, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int));
  alloc.device                         += devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int);
  mycudaMalloc((void **)gsync1_make_tree, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int));
  alloc.device                         += devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int);
  mycudaMalloc((void **)gsync2_make_tree, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int));
  alloc.device                         += devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int);
  mycudaMalloc((void **)gsync3_make_tree, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int));
  alloc.device                         += devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE * sizeof(int);
  mycudaMalloc((void **)  gmem_link_tree, devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int));
  alloc.device                         += devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int);
  mycudaMalloc((void **)gsync0_link_tree, devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int));
  alloc.device                         += devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int);
  mycudaMalloc((void **)gsync1_link_tree, devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int));
  alloc.device                         += devProp.numSM * NBLOCKS_PER_SM_LINK_TREE * sizeof(int);
  initGsync_kernel<<<1, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE>>>(devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE, *gsync0_make_tree, *gsync1_make_tree);
  getLastCudaError("initGsync_kernel");
  initGsync_kernel<<<1, devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE>>>(devProp.numSM * NBLOCKS_PER_SM_MAKE_TREE, *gsync2_make_tree, *gsync3_make_tree);
  getLastCudaError("initGsync_kernel");
  initGsync_kernel<<<1, devProp.numSM * NBLOCKS_PER_SM_LINK_TREE>>>(devProp.numSM * NBLOCKS_PER_SM_LINK_TREE, *gsync0_link_tree, *gsync1_link_tree);
  getLastCudaError("initGsync_kernel");
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  buf->more0  = *more0Buf;
  buf->more1  = *more1Buf;
  buf->rjmax  = *rjmaxBuf;
  buf->fail   = *fail_dev;
  buf->gsync0 = *gsync0;
  buf->gsync1 = *gsync1;
  //-----------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
  buf->  gmem_make_tree = *  gmem_make_tree;
  buf->  gmem_link_tree = *  gmem_link_tree;
  buf->gsync0_make_tree = *gsync0_make_tree;
  buf->gsync1_make_tree = *gsync1_make_tree;
  buf->gsync2_make_tree = *gsync2_make_tree;
  buf->gsync3_make_tree = *gsync3_make_tree;
  buf->gsync0_link_tree = *gsync0_link_tree;
  buf->gsync1_link_tree = *gsync1_link_tree;
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
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
 int  *more0Buf, int  *more1Buf, real  *rjmaxBuf, int  *fail_dev)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  mycudaFree    ( more_dev);
  mycudaFree    (   pj_dev);
  mycudaFree    (   mj_dev);
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  mycudaFree    (  n2c_dev);
#endif//CALC_MULTIPOLE_ON_DEVICE
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  mycudaFree    (niSub_dev);
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  //-----------------------------------------------------------------------
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
  mycudaFreeHost( more_hst);
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(MAKE_TREE_ON_DEVICE)
#   if  (!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
  mycudaFreeHost(   pj_hst);
  mycudaFreeHost(   mj_hst);
#endif//(!defined(SERIALIZED_EXECUTION) && (!defined(BUILD_LET_ON_DEVICE) || defined(LET_COMMUNICATION_VIA_HOST))) || !defined(CALC_MULTIPOLE_ON_DEVICE)
#ifndef MAKE_TREE_ON_DEVICE
  mycudaFreeHost(  n2c_hst);
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  mycudaFreeHost(niSub_hst);
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifdef  GADGET_MAC
  mycudaFree(mac_dev);
#endif//GADGET_MAC
  //-----------------------------------------------------------------------
#ifdef  CALC_MULTIPOLE_ON_DEVICE
#       ifdef  WS93_MAC
  mycudaFree    ( mr2_dev);
#       endif//WS93_MAC
  mycudaFree    (bmax_dev);
#else///CALC_MULTIPOLE_ON_DEVICE
#       ifdef  WS93_MAC
  mycudaFreeHost( mr2_hst);
#       endif//WS93_MAC
  mycudaFreeHost(bmax_hst);
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  mycudaFree    (more0Buf);
  mycudaFree    (more1Buf);
  mycudaFree    (rjmaxBuf);
#else///CALC_MULTIPOLE_ON_DEVICE
  mycudaFreeHost(more0Buf);
  mycudaFreeHost(more1Buf);
  mycudaFreeHost(rjmaxBuf);
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  mycudaFree(gsync0);
  mycudaFree(gsync1);
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
  mycudaFree(fail_dev);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
  mycudaFree(gmem_make_tree);  mycudaFree(gsync0_make_tree);  mycudaFree(gsync1_make_tree);  mycudaFree(gsync2_make_tree);  mycudaFree(gsync3_make_tree);
  mycudaFree(gmem_link_tree);  mycudaFree(gsync0_link_tree);  mycudaFree(gsync1_link_tree);
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
extern "C"
void setTreeNode_dev
(const size_t Nj, const soaTreeNode dev, const soaTreeNode hst
#ifdef  CALC_MULTIPOLE_ON_DEVICE
 , const size_t Ni
#endif//CALC_MULTIPOLE_ON_DEVICE
#ifdef EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* send tree node as pseudo j-particles from the host to the device using the default CUDA stream */
  //-----------------------------------------------------------------------
#ifdef EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaMemcpy(dev.more, hst.more, sizeof(uint) * Nj, cudaMemcpyHostToDevice));
  //-----------------------------------------------------------------------
#ifdef  CALC_MULTIPOLE_ON_DEVICE
  checkCudaErrors(cudaMemcpy(dev.node2cell, hst.node2cell, sizeof(      int) * Nj, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(dev.jtag     , hst.jtag     , sizeof(      int) * Ni, cudaMemcpyHostToDevice));
#else///CALC_MULTIPOLE_ON_DEVICE
  checkCudaErrors(cudaMemcpy(dev.jpos     , hst.jpos     , sizeof(jparticle) * Nj, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(dev.mj       , hst.mj       , sizeof(jmass    ) * Nj, cudaMemcpyHostToDevice));
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  checkCudaErrors(cudaMemcpy(dev.niSub, hst.niSub, sizeof(int) * Nj, cudaMemcpyHostToDevice));
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  //-----------------------------------------------------------------------
#ifdef EXEC_BENCHMARK
  stopStopwatch(&(elapsed->setTreeNode_dev));
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
//-------------------------------------------------------------------------
/* arrays to store properties of tree cells */
//-------------------------------------------------------------------------
extern "C"
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
 soaTreeCell *dev, soaTreeCell *hst)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const size_t num = NUM_ALLOC_TREE_CELL;
  //-----------------------------------------------------------------------
  mycudaMalloc((void **)    cell_dev, num * sizeof(treecell));  alloc.device += num * sizeof(treecell);
  mycudaMalloc((void **)    leaf_dev, num * sizeof(    bool));  alloc.device += num * sizeof(    bool);
  mycudaMalloc((void **)    node_dev, num * sizeof(    uint));  alloc.device += num * sizeof(    uint);
#ifdef  MAKE_TREE_ON_DEVICE
  mycudaMalloc((void **)    hkey_dev, num * sizeof(   PHint));  alloc.device += num * sizeof(   PHint);
  mycudaMalloc((void **)  parent_dev, num * sizeof(    uint));  alloc.device += num * sizeof(    uint);
  mycudaMalloc((void **)children_dev, num * sizeof(    uint));  alloc.device += num * sizeof(    uint);
  mycudaMalloc((void **) leafLev_dev, num * sizeof(     int));  alloc.device +=       sizeof(     int);
  mycudaMalloc((void **) numCell_dev, num * sizeof(     int));  alloc.device +=       sizeof(     int);
  mycudaMalloc((void **) numNode_dev, num * sizeof(     int));  alloc.device +=       sizeof(     int);
  mycudaMalloc((void **) scanNum_dev, num * sizeof(     int));  alloc.device +=       sizeof(     int);
#else///MAKE_TREE_ON_DEVICE
  mycudaMallocHost((void **)cell_hst, num * sizeof(treecell));  alloc.host += num * sizeof(treecell);
  mycudaMallocHost((void **)leaf_hst, num * sizeof(    bool));  alloc.host += num * sizeof(    bool);
  mycudaMallocHost((void **)node_hst, num * sizeof(    uint));  alloc.host += num * sizeof(    uint);
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
  mycudaMallocHost((void **)cell_hst, num * sizeof(treecell));  alloc.host += num * sizeof(treecell);
  mycudaMallocHost((void **)leaf_hst, num * sizeof(    bool));  alloc.host += num * sizeof(    bool);
  mycudaMallocHost((void **)node_hst, num * sizeof(    uint));  alloc.host += num * sizeof(    uint);
#endif//defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
  //-----------------------------------------------------------------------
  mycudaMalloc    ((void **)info_dev, NUM_PHKEY_LEVEL * sizeof(PHinfo));  alloc.device += NUM_PHKEY_LEVEL * sizeof(PHinfo);
#ifdef  MAKE_TREE_ON_DEVICE
  initPHinfo_dev(*info_dev);
#else///MAKE_TREE_ON_DEVICE
  mycudaMallocHost((void **)info_hst, NUM_PHKEY_LEVEL * sizeof(PHinfo));  alloc.host   += NUM_PHKEY_LEVEL * sizeof(PHinfo);
  initPHinfo(*info_hst);
  checkCudaErrors(cudaMemcpy(*info_dev, *info_hst, NUM_PHKEY_LEVEL * sizeof(PHinfo), cudaMemcpyHostToDevice));
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
  mycudaMallocHost((void **)info_hst, NUM_PHKEY_LEVEL * sizeof(PHinfo));  alloc.host   += NUM_PHKEY_LEVEL * sizeof(PHinfo);
#endif//defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
  //-----------------------------------------------------------------------
  dev->cell  = *cell_dev;
  dev->leaf  = *leaf_dev;
  dev->ptag  = *node_dev;
  dev->level = *info_dev;
#ifdef  MAKE_TREE_ON_DEVICE
  dev->hkey     = *    hkey_dev;
  dev->parent   = *  parent_dev;
  dev->children = *children_dev;
#else///MAKE_TREE_ON_DEVICE
  hst->cell  = *cell_hst;
  hst->leaf  = *leaf_hst;
  hst->ptag  = *node_hst;
  hst->level = *info_hst;
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
  hst->cell  = *cell_hst;
  hst->leaf  = *leaf_hst;
  hst->ptag  = *node_hst;
  hst->level = *info_hst;
#endif//defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
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
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  mycudaFree(cell_dev);
  mycudaFree(leaf_dev);
  mycudaFree(node_dev);
  mycudaFree(info_dev);
#ifdef  MAKE_TREE_ON_DEVICE
  mycudaFree(    hkey_dev);
  mycudaFree(  parent_dev);
  mycudaFree(children_dev);
  mycudaFree( leafLev_dev);
  mycudaFree( numCell_dev);
  mycudaFree( numNode_dev);
  mycudaFree( scanNum_dev);
#else///MAKE_TREE_ON_DEVICE
  mycudaFreeHost(cell_hst);
  mycudaFreeHost(leaf_hst);
  mycudaFreeHost(node_hst);
  mycudaFreeHost(info_hst);
#endif//MAKE_TREE_ON_DEVICE
#   if  defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
  mycudaFreeHost(cell_hst);
  mycudaFreeHost(leaf_hst);
  mycudaFreeHost(node_hst);
  mycudaFreeHost(info_hst);
#endif//defined(MAKE_TREE_ON_DEVICE) && defined(COUNT_INTERACTIONS)
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//defined(CALC_MULTIPOLE_ON_DEVICE) || defined(MAKE_TREE_ON_DEVICE)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  CALC_MULTIPOLE_ON_DEVICE
//-------------------------------------------------------------------------
#ifndef MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------
extern "C"
void setTreeCell_dev(const size_t num, const soaTreeCell dev, const soaTreeCell hst
#ifdef EXEC_BENCHMARK
		     , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
		     )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  /* send tree cell from the host to the device using the default CUDA stream */
  checkCudaErrors(cudaMemcpy(dev.cell, hst.cell, sizeof(treecell) * num, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(dev.leaf, hst.leaf, sizeof(    bool) * num, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(dev.ptag, hst.ptag, sizeof(    uint) * num, cudaMemcpyHostToDevice));
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaMemcpy(dev.level, hst.level, sizeof(PHinfo) * NUM_PHKEY_LEVEL, cudaMemcpyHostToDevice));
  //-----------------------------------------------------------------------
#ifndef EXEC_BENCHMARK
#       ifdef  GENERATE_PHKEY_ON_DEVICE
  checkCudaErrors(cudaDeviceSynchronize());
#       endif//GENERATE_PHKEY_ON_DEVICE
#else///EXEC_BENCHMARK
  stopStopwatch(&(elapsed->setTreeCell_dev));
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  GADGET_MAC
//-------------------------------------------------------------------------
/* rewrite MAC for initial step (GADGET-MAC requires pre-estimated acceleration ) */
//-------------------------------------------------------------------------
/* Ni      :: input          :: number of i-particles */
/* ai      ::         output :: acceleration and potential of i-particles */
/* Nj      :: input          :: number of j-particles */
/* pj      :: input / output :: position and MAC of j-particles */
/* mac_bak ::         output :: tentative stock of GADGET-MAC of all j-particles */
/* bmax    :: input          :: size of distribution of j-particles */
//-------------------------------------------------------------------------
__global__ void enforceBarnesHutMAC
(const int Ni, acceleration * RESTRICT ai,
 const int Nj, jparticle * RESTRICT pj, real * RESTRICT mac_bak, READ_ONLY real * RESTRICT bmax)
{
  //-----------------------------------------------------------------------
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( gidx < Ni ){
    //---------------------------------------------------------------------
    /* set old acceleration as |a| is unity */
    const acceleration iacc = {UNITY, ZERO, ZERO, ZERO};
    //---------------------------------------------------------------------
    ai[gidx] = iacc;
    //---------------------------------------------------------------------
  }/* if( gidx < Ni ){ */
  //-----------------------------------------------------------------------
  if( gidx < Nj ){
    //---------------------------------------------------------------------
    /* load position and MAC of j-particle */
    jparticle jpos = pj[gidx];
    //---------------------------------------------------------------------
    /* save the GADGET-MAC */
    mac_bak[gidx] = jpos.w;
    //---------------------------------------------------------------------
    /* modify MAC for pseudo particles */
    /* squared size for the real particle is set to be negative (-UNITY) */
    if( jpos.w > -HALF ){
      //-------------------------------------------------------------------
      /* new MAC for initial step is Barnes-Hut criterion */
      /* jcnd.w = bmax^4 / theta^4 in this case */
      //-------------------------------------------------------------------
      const real bjmax = bmax[gidx];
#if 1
      /* theta is assumed to be 0.5 for simplicity */
      jpos.w = TWO * bjmax;
#else
      /* theta is assumed to be unity for simplicity */
      jpos.w =       bjmax;
#endif
      jpos.w *= jpos.w;
      jpos.w *= jpos.w;
      pj[gidx] = jpos;
      //-------------------------------------------------------------------
    }/* if( jpos.w > -HALF ){ */
    //---------------------------------------------------------------------
  }/* if( gidx < Nj ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
/* recover MAC for further steps */
//-------------------------------------------------------------------------
/* Nj      :: input          :: number of j-particles */
/* pj      :: input / output :: position and MAC of j-particles */
/* mac_bak :: input          :: tentative stock of GADGET-MAC of all j-particles */
//-------------------------------------------------------------------------
__global__ void recoverGADGET_MAC
(const int Nj, jparticle * RESTRICT pj, READ_ONLY real * RESTRICT mac_bak)
{
  //-----------------------------------------------------------------------
  const int jj = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  if( jj < Nj ){
    //---------------------------------------------------------------------
    /* load position and MAC of j-particle */
    jparticle jpos = pj[jj];
    const real mac = mac_bak[jj];
    //---------------------------------------------------------------------
    /* recover the original MAC */
    jpos.w = mac;
    //---------------------------------------------------------------------
    /* return the original position and MAC of j-particle */
    pj[jj] = jpos;
    //---------------------------------------------------------------------
  }/* if( jj < Nj ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
/* Ni       :: input          :: number of i-particles */
/*   ai_dev ::         output :: acceleration and potential of N-body particle */
/* Nj       :: input          :: number of j-particles */
/*   pj_dev :: input / output :: position and MAC of pseudo N-body particle as j-particles */
/*  mac_dev ::         output :: tentative stock of GADGET-MAC of all j-particles */
/* bmax_dev :: input          :: size of distribution of j-particles */
//-------------------------------------------------------------------------
extern "C"
void enforceBarnesHutMAC_dev(const int Ni, const iparticle pi, const int Nj, const soaTreeNode pj)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  /* NOTE: Nj is always greater than Ni */
  int Njrem = BLOCKSIZE(Nj, 1024);
  int Nirem = BLOCKSIZE(Ni, 1024);
  if( Njrem <= MAX_BLOCKS_PER_GRID )
    enforceBarnesHutMAC<<<Njrem, 1024>>>(Ni, pi.acc, Nj, pj.jpos, pj.mac, pj.bmax);
  else{
    //---------------------------------------------------------------------
    const int Niter = BLOCKSIZE(Njrem, MAX_BLOCKS_PER_GRID);
    int hjidx = 0;
    int hiidx = 0;
    //---------------------------------------------------------------------
    for(int iter = 0; iter < Niter; iter++){
      //-------------------------------------------------------------------
      int Njblck = MAX_BLOCKS_PER_GRID;      if( Njblck > Njrem )	Njblck = Njrem;
      int Niblck = MAX_BLOCKS_PER_GRID;      if( Niblck > Nirem )	Niblck = Nirem;
      //-------------------------------------------------------------------
      int Nisub = Niblck * 1024;
      int Njsub = Njblck * 1024;
      enforceBarnesHutMAC<<<Njblck, 1024>>>(Nisub, &pi.acc[hiidx], Njsub, &pj.jpos[hjidx], &pj.mac[hjidx], &pj.bmax[hjidx]);
      //-------------------------------------------------------------------
      hjidx += Njsub;      Njrem -= Njblck;
      hiidx += Nisub;      Nirem -= Niblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
  /* enforceBarnesHutMAC<<<BLOCKSIZE(Nj, 1024), 1024>>>(Ni, pi.acc, Nj, pj.jpos, pj.mac, pj.bmax); */
  getLastCudaError("enforceBarnesHutMAC");
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
/* Nj       :: input          :: number of j-particles */
/*   pj_dev :: input / output :: position and MAC of pseudo N-body particle as j-particles */
/*  mac_dev :: input          :: tentative stock of GADGET-MAC of all j-particles */
//-------------------------------------------------------------------------
extern "C"
void recoverGADGET_MAC_dev(const int Nj, const soaTreeNode pj)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  int Nrem = BLOCKSIZE(Nj, 1024);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    recoverGADGET_MAC<<<Nrem, 1024>>>(Nj, pj.jpos, pj.mac);
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
      int Nsub = Nblck * 1024;
      recoverGADGET_MAC<<<Nblck, 1024>>>(Nsub, &pj.jpos[hidx], &pj.mac[hidx]);
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
  /* recoverGADGET_MAC<<<BLOCKSIZE(Nj, 1024), 1024>>>(Nj, pj.jpos, pj.mac); */
  getLastCudaError("recoverGADGET_MAC");
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//GADGET_MAC
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* initialize mass and position of pseudo j-particles */
//-------------------------------------------------------------------------
/* acc ::         output :: acceleration and potential of N-body particles */
//-------------------------------------------------------------------------
__global__ void initTreeBody_kernel
(jparticle * RESTRICT pj, jmass * RESTRICT mj, real * RESTRICT bmax
#ifdef  WS93_MAC
 , real * RESTRICT mr2
#endif//WS93_MAC
)
{
  //-----------------------------------------------------------------------
  const int gidx = GLOBALIDX_X1D;
  //-----------------------------------------------------------------------
  pj  [gidx] = zero_pj;
  mj  [gidx] = zero_mj;
  bmax[gidx] = ZERO;
#ifdef  WS93_MAC
  mr2 [gidx] = ZERO;
#endif//WS93_MAC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* set N-body particles as real j-particles */
//-------------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int gidx = GLOBALIDX_X1D;/* ((BLOCKIDX_X1D) * (BLOCKDIM_X1D) + (THREADIDX_X1D)) */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  if( gidx < Ni ){
    //---------------------------------------------------------------------
    /* load an N-body particle */
    const position ipos = pi  [gidx];
    const int      jidx = jtag[gidx];
    //---------------------------------------------------------------------
    /* set the N-body particle as a real particle */
    jparticle jpos;
    jpos.x = ipos.x;
    jpos.y = ipos.y;
    jpos.z = ipos.z;
    jpos.w = -UNITY;/* squared size for the real particle is set to be negative */
    //---------------------------------------------------------------------
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
    const jmass mj_tmp = {ipos.m, eps2};
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
    const jmass mj_tmp =  ipos.m;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* commit the new particle to the j-particle array */
    pj[jidx] = jpos;
    mj[jidx] = mj_tmp;
#ifdef  WS93_MAC
    mr2[jidx] = ipos.m * (jpos.x * jpos.x + jpos.y * jpos.y + jpos.z * jpos.z);
#endif//WS93_MAC
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* parallel prefix sum within a group of TSUB threads (TSUB <= 32 to use implicit synchronization) */
/* type of prefix sum is inclusive */
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
__device__ __forceinline__  int prefixSumTsub(const int psum,                                           const int lane)
#else///USE_WARP_SHUFFLE_FUNC_MAC
__device__ __forceinline__ void prefixSumTsub(const int psum, volatile int_real * smem, const int tidx, const int lane)
#endif//USE_WARP_SHUFFLE_FUNC_MAC
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
  int val = psum;
  int tmp;
#   if  TSUB_MAC >=  2
  tmp = __shfl_up(val,  1, TSUB_MAC);  if( lane >=  1 )    val += tmp;
#   if  TSUB_MAC >=  4
  tmp = __shfl_up(val,  2, TSUB_MAC);  if( lane >=  2 )    val += tmp;
#   if  TSUB_MAC >=  8
  tmp = __shfl_up(val,  4, TSUB_MAC);  if( lane >=  4 )    val += tmp;
#   if  TSUB_MAC >= 16
  tmp = __shfl_up(val,  8, TSUB_MAC);  if( lane >=  8 )    val += tmp;
#   if  TSUB_MAC == 32
  tmp = __shfl_up(val, 16, TSUB_MAC);  if( lane >= 16 )    val += tmp;
#endif//TSUB_MAC == 32
#endif//TSUB_MAC >= 16
#endif//TSUB_MAC >=  8
#endif//TSUB_MAC >=  4
#endif//TSUB_MAC >=  2
  return (val);
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
  smem[tidx].i = psum;
  //-----------------------------------------------------------------------
# if  TSUB_MAC >=  2
  if( lane >=  1 )    smem[tidx].i += smem[tidx -  1].i;
# if  TSUB_MAC >=  4
  if( lane >=  2 )    smem[tidx].i += smem[tidx -  2].i;
# if  TSUB_MAC >=  8
  if( lane >=  4 )    smem[tidx].i += smem[tidx -  4].i;
# if  TSUB_MAC >= 16
  if( lane >=  8 )    smem[tidx].i += smem[tidx -  8].i;
# if  TSUB_MAC == 32
  if( lane >= 16 )    smem[tidx].i += smem[tidx - 16].i;
#endif//TSUB_MAC == 32
#endif//TSUB_MAC >= 16
#endif//TSUB_MAC >=  8
#endif//TSUB_MAC >=  4
#endif//TSUB_MAC >=  2
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
__device__ __forceinline__ real accumulateRealTsub(const real  sum)
#else///USE_WARP_SHUFFLE_FUNC_MAC
__device__ __forceinline__ void accumulateRealTsub(      real *sum, volatile int_real * smem, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC_MAC
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
  real val = sum;
#   if  TSUB_MAC >=  2
  real tmp;
  tmp = __shfl_xor(val,  1, TSUB_MAC);  val += tmp;
#   if  TSUB_MAC >=  4
  tmp = __shfl_xor(val,  2, TSUB_MAC);  val += tmp;
#   if  TSUB_MAC >=  8
  tmp = __shfl_xor(val,  4, TSUB_MAC);  val += tmp;
#   if  TSUB_MAC >= 16
  tmp = __shfl_xor(val,  8, TSUB_MAC);  val += tmp;
#   if  TSUB_MAC == 32
  tmp = __shfl_xor(val, 16, TSUB_MAC);  val += tmp;
#endif//TSUB_MAC == 32
#endif//TSUB_MAC >= 16
#endif//TSUB_MAC >=  8
#endif//TSUB_MAC >=  4
#endif//TSUB_MAC >=  2
  return (__shfl(val, 0, TSUB_MAC));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
  smem[tidx].r = *sum;
  //-----------------------------------------------------------------------
#   if  TSUB_MAC >=  2
  real tmp;
  tmp = smem[tidx ^  1].r;  *sum += tmp;  smem[tidx].r = *sum;
#   if  TSUB_MAC >=  4
  tmp = smem[tidx ^  2].r;  *sum += tmp;  smem[tidx].r = *sum;
#   if  TSUB_MAC >=  8
  tmp = smem[tidx ^  4].r;  *sum += tmp;  smem[tidx].r = *sum;
#   if  TSUB_MAC >= 16
  tmp = smem[tidx ^  8].r;  *sum += tmp;  smem[tidx].r = *sum;
#   if  TSUB_MAC == 32
  tmp = smem[tidx ^ 16].r;  *sum += tmp;  smem[tidx].r = *sum;
#endif//TSUB_MAC == 32
#endif//TSUB_MAC >= 16
#endif//TSUB_MAC >=  8
#endif//TSUB_MAC >=  4
#endif//TSUB_MAC >=  2
  //-----------------------------------------------------------------------
  *sum = smem[head].r;
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
__device__ __forceinline__  int accumulateIntTsub(const int  sum)
#else///USE_WARP_SHUFFLE_FUNC_MAC
__device__ __forceinline__ void accumulateIntTsub(      int *sum, volatile int_real * smem, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC_MAC
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
  int val = sum;
#   if  TSUB_MAC >=  2
  int tmp;
  tmp = __shfl_xor(val,  1, TSUB_MAC);  val += tmp;
#   if  TSUB_MAC >=  4
  tmp = __shfl_xor(val,  2, TSUB_MAC);  val += tmp;
#   if  TSUB_MAC >=  8
  tmp = __shfl_xor(val,  4, TSUB_MAC);  val += tmp;
#   if  TSUB_MAC >= 16
  tmp = __shfl_xor(val,  8, TSUB_MAC);  val += tmp;
#   if  TSUB_MAC == 32
  tmp = __shfl_xor(val, 16, TSUB_MAC);  val += tmp;
#endif//TSUB_MAC == 32
#endif//TSUB_MAC >= 16
#endif//TSUB_MAC >=  8
#endif//TSUB_MAC >=  4
#endif//TSUB_MAC >=  2
  return (__shfl(val, 0, TSUB_MAC));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
  smem[tidx].i = *sum;
  //-----------------------------------------------------------------------
#   if  TSUB_MAC >=  2
  int tmp;
  tmp = smem[tidx ^  1].i;  *sum += tmp;  smem[tidx].i = *sum;
#   if  TSUB_MAC >=  4
  tmp = smem[tidx ^  2].i;  *sum += tmp;  smem[tidx].i = *sum;
#   if  TSUB_MAC >=  8
  tmp = smem[tidx ^  4].i;  *sum += tmp;  smem[tidx].i = *sum;
#   if  TSUB_MAC >= 16
  tmp = smem[tidx ^  8].i;  *sum += tmp;  smem[tidx].i = *sum;
#   if  TSUB_MAC == 32
  tmp = smem[tidx ^ 16].i;  *sum += tmp;  smem[tidx].i = *sum;
#endif//TSUB_MAC == 32
#endif//TSUB_MAC >= 16
#endif//TSUB_MAC >=  8
#endif//TSUB_MAC >=  4
#endif//TSUB_MAC >=  2
  //-----------------------------------------------------------------------
  *sum = smem[head].i;
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
__device__ __forceinline__ real getMaximumRealTsub(const real  max)
#else///USE_WARP_SHUFFLE_FUNC_MAC
__device__ __forceinline__ void getMaximumRealTsub(      real *max, volatile int_real * smem, const int tidx, const int head)
#endif//USE_WARP_SHUFFLE_FUNC_MAC
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
  real val = max;
#   if  TSUB_MAC >=  2
  real tmp;
  tmp = __shfl_xor(val,  1, TSUB_MAC);  if( tmp > val )    val = tmp;
#   if  TSUB_MAC >=  4
  tmp = __shfl_xor(val,  2, TSUB_MAC);  if( tmp > val )    val = tmp;
#   if  TSUB_MAC >=  8
  tmp = __shfl_xor(val,  4, TSUB_MAC);  if( tmp > val )    val = tmp;
#   if  TSUB_MAC >= 16
  tmp = __shfl_xor(val,  8, TSUB_MAC);  if( tmp > val )    val = tmp;
#   if  TSUB_MAC == 32
  tmp = __shfl_xor(val, 16, TSUB_MAC);  if( tmp > val )    val = tmp;
#endif//TSUB_MAC == 32
#endif//TSUB_MAC >= 16
#endif//TSUB_MAC >=  8
#endif//TSUB_MAC >=  4
#endif//TSUB_MAC >=  2
  return (__shfl(val, 0, TSUB_MAC));
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
  smem[tidx].r = *max;
  //-----------------------------------------------------------------------
#   if  TSUB_MAC >=  2
  real tmp;
  tmp = smem[tidx ^  1].r;  if( tmp > *max ){    *max = tmp;  }  smem[tidx].r = *max;
#   if  TSUB_MAC >=  4
  tmp = smem[tidx ^  2].r;  if( tmp > *max ){    *max = tmp;  }  smem[tidx].r = *max;
#   if  TSUB_MAC >=  8
  tmp = smem[tidx ^  4].r;  if( tmp > *max ){    *max = tmp;  }  smem[tidx].r = *max;
#   if  TSUB_MAC >= 16
  tmp = smem[tidx ^  8].r;  if( tmp > *max ){    *max = tmp;  }  smem[tidx].r = *max;
#   if  TSUB_MAC == 32
  tmp = smem[tidx ^ 16].r;  if( tmp > *max ){    *max = tmp;  }  smem[tidx].r = *max;
#endif//TSUB_MAC == 32
#endif//TSUB_MAC >= 16
#endif//TSUB_MAC >=  8
#endif//TSUB_MAC >=  4
#endif//TSUB_MAC >=  2
  //-----------------------------------------------------------------------
  *max = smem[head].r;
  //-----------------------------------------------------------------------
#endif//USE_WARP_SHUFFLE_FUNC_MAC
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* calculate multipole moment(s) of tree cells based on the width-first search */
//-------------------------------------------------------------------------
/* level     :: input          :: head index and number of tree cells contained in the corresponding hierarchy */
/* cell      :: input          :: head index and number of N-body particles contained in the corresponding tree cell */
/* leaf      :: input          :: a flag to remember the corresponding tree cell is either leaf(true) of node(false) */
/* pi        :: input          :: position and mass of N-body particles */
/* node      :: input          :: head index and number of pseudo particles of the corresponding tree cell (a leaf cell contains up to NCRIT particles) */
/* more      :: input          :: head index and number of child particles of the corresponding j-particle (a leaf cell contains up to NCRIT particles) */
/* node2cell :: input          :: index of the tree cell corresponding a pseudo particle */
/* pj        :: input / output :: position and squared radius of pseudo N-body particle as j-particles */
/* mj        :: input / output :: mass of pseudo N-body particle as j-particles */
/* mr2       :: input / output :: mass times r squared of pseudo N-body particle as j-particles */
/* bmax      :: input / output :: size of pseudo N-body particle as j-particles */
//-------------------------------------------------------------------------
__global__ void __launch_bounds__(NTHREADS_MAC, NBLOCKS_PER_SM_MAC) calcMultipole_kernel
     (const int bottomLev, READ_ONLY PHinfo * RESTRICT level,
      READ_ONLY treecell * RESTRICT cell, READ_ONLY bool * RESTRICT leaf, READ_ONLY position * RESTRICT pi,
      READ_ONLY uint * RESTRICT node, READ_ONLY uint * RESTRICT more, READ_ONLY int * RESTRICT node2cell,
      jparticle * RESTRICT pj, jmass * RESTRICT mj, real * RESTRICT bmax,
      int * RESTRICT more0Buf, int * RESTRICT more1Buf, real * RESTRICT rjmaxBuf, int * RESTRICT overflow,
      int * RESTRICT gsync0, int * RESTRICT gsync1
#ifdef  WS93_MAC
      , real * RESTRICT mr2
#endif//WS93_MAC
#ifdef  COUNT_INTERACTIONS
      , tree_stats * RESTRICT stats
#endif//COUNT_INTERACTIONS
      )
{
  //-----------------------------------------------------------------------
  /* identify thread properties */
  //-----------------------------------------------------------------------
  const int tidx = THREADIDX_X1D;
  const int gidx = GLOBALIDX_X1D;
  const int bidx =  BLOCKIDX_X1D;
  const int bnum =   GRIDDIM_X1D;
  //-----------------------------------------------------------------------
  const int lane = tidx & (TSUB_MAC - 1);
  const int head = tidx - lane;
#ifndef USE_WARP_SHUFFLE_FUNC_MAC
  const int tail = head + TSUB_MAC - 1;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
  const int hbuf = (head / TSUB_MAC) * TSUB_MAC * NBUF_MAC;/* head index of the shared array close and queue within a thread group */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
  int smem;
#else///USE_WARP_SHUFFLE_FUNC_MAC
  __shared__  int_real  smem[NTHREADS_MAC];
#endif//USE_WARP_SHUFFLE_FUNC_MAC
  __shared__ jparticle pj_sm[NTHREADS_MAC];
  __shared__      real rjbuf[NTHREADS_MAC * NBUF_MAC];
  __shared__  int      list0[NTHREADS_MAC * NBUF_MAC];
  __shared__  int      list1[NTHREADS_MAC * NBUF_MAC];
  __shared__  int      pjidx[NTHREADS_MAC * NBUF_MAC];
  /* pjidx would be able to be removed by integrating with list; however, reducing SM usage from 24KB to 20KB has no effect @ Ttot = 512 */
  //-----------------------------------------------------------------------
  /* head index of remote buffers */
  const int bufHead = ((head / TSUB_MAC) + bidx * NGROUPS_MAC) * NUM_ALLOC_MACBUF;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate multipole moment(s) of pseudo j-particles */
  //-----------------------------------------------------------------------
  /* for(int levelIdx = (NUM_PHKEY_LEVEL - 2); levelIdx >= 0; levelIdx--){ */
  for(int levelIdx = bottomLev; levelIdx >= 0; levelIdx--){
    //---------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
      tree_stats local = {ZERO, ZERO, ZERO, ZERO, 0, 0};
#endif//COUNT_INTERACTIONS
    //---------------------------------------------------------------------
    /* load a tree cell and evaluate either the cell is a node or a leaf */
    /* if cidx is not less than tail, then the thread has a null leaf-cell */
    //---------------------------------------------------------------------
    const PHinfo clev = level[levelIdx];
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------
    /* loop to set a maximum number for # of blocks */
    for(int ii = 0; ii < bnum * BLOCKSIZE(clev.num, bnum * NGROUPS_MAC); ii += bnum){
      //-------------------------------------------------------------------
      const      int cidx =  clev.head + (gidx + ii * NTHREADS_MAC) / TSUB_MAC;
      const treecell root = (cidx < clev.head + clev.num) ? (cell[cidx]) : (null_cell_device);
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
      /* extend the pseudo particle chain */
      //-------------------------------------------------------------------
      if( root.head != NULL_CELL ){
	//-----------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
	local.cellNum++;
#endif//COUNT_INTERACTIONS
        //-----------------------------------------------------------------
	if( !leaf[cidx] ){
	  //---------------------------------------------------------------
	  /* when the tree cell is a node cell, then calculate multipole moment(s) of the cell */
	  //---------------------------------------------------------------
	  jparticle jcom = {ZERO, ZERO, ZERO, ZERO};
	  //---------------------------------------------------------------
#ifdef  WS93_MAC
	  real mjrj2 = ZERO;
#endif//WS93_MAC
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
	  /* sum up multipole moment(s) of child cells */
	  //---------------------------------------------------------------
	  /* only leaf cells can point multiple tree nodes */
	  const  int jidx  = node[cidx] & IDXMASK;
	  uint more_tmp = more[jidx];
	  int cnum  = (more_tmp >> IDXBITS) + 1;
	  int chead = (more_tmp  & IDXMASK);
	  //---------------------------------------------------------------
	  /* calculate multipole moments of the pseudo particle group */
	  for(int jj = lane; jj < cnum; jj += TSUB_MAC){
	    //-------------------------------------------------------------
	    const int       sidx = chead + jj;
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const real      mass = mj[sidx].mass;
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const real      mass = mj[sidx];
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const jparticle jpos = pj[sidx];
	    //-------------------------------------------------------------
	    /* calculate total mass */
	    jcom.w += mass;
	    /* calculate center-of-mass */
	    jcom.x += mass * jpos.x;
	    jcom.y += mass * jpos.y;
	    jcom.z += mass * jpos.z;
	    /* calculate trace of quadrupole moment */
#ifdef  WS93_MAC
	    mjrj2 += mr2[sidx];
#endif//WS93_MAC
	    //-------------------------------------------------------------
	  }/* for(int jj = lane; jj < cnum; jj += TSUB_MAC){ */
	  //---------------------------------------------------------------
	  /* sum up partial sums within tsub threads */
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
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
	  mjrj2 = accumulateRealTsub(mjrj2);
#else///USE_WARP_SHUFFLE_FUNC_MAC
	  accumulateRealTsub(&mjrj2, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_MAC
	  if( lane == 0 )
	    mr2[jidx] = mjrj2;
#endif//WS93_MAC
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* calculate multipole moments of the pseudo particle group */
	  //---------------------------------------------------------------
	  /* calculate center-of-mass */
	  real mtot = jcom.w;
	  real minv = UNITY / mtot;
	  jcom.x *= minv;
	  jcom.y *= minv;
	  jcom.z *= minv;
	  /* calculate trace of quadrupole moment */
#ifdef  WS93_MAC
	  real B2 = mjrj2 - mtot * (jcom.x * jcom.x + jcom.y * jcom.y + jcom.z * jcom.z);
#endif//WS93_MAC
	  //---------------------------------------------------------------
	  /* commit a pseudo particle */
	  if( lane == 0 ){
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const jmass mj_loc = {mtot, ZERO};
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const jmass mj_loc =  mtot;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    mj[jidx] =  mj_loc;
	  }
	  //---------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
	  local.nodeNum++;
	  local.mjMean +=        mtot;
	  local.mjSdev += mtot * mtot;
#endif//COUNT_INTERACTIONS
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
	  /* estimate size of particle distribution */
	  //---------------------------------------------------------------
	  /* initialize list of examined tree nodes and related variables */
	  int inum = root.num;
	  int Ntry = 1;
	  if( lane == 0 )	    list0[hbuf] = cidx;
	  //---------------------------------------------------------------
	  /* pick up NI_BMAX_ESTIMATE i-particles in maximum to estimate bmax */
	  while( inum > NI_BMAX_ESTIMATE ){
	    //-------------------------------------------------------------
	    real rmin = ZERO;
	    //-------------------------------------------------------------
	    int Nloc = 0;
	    int Nbuf = 0;
	    //-------------------------------------------------------------
	    int Niter = BLOCKSIZE(Ntry, TSUB_MAC * NBUF_MAC);
	    for(int iter = 0; iter < Niter; iter++){
	      //-----------------------------------------------------------
	      const int Nsweep = (Ntry > (TSUB_MAC * NBUF_MAC)) ? (TSUB_MAC * NBUF_MAC) : Ntry;
	      const int ibuf_loop = BLOCKSIZE(Nsweep, TSUB_MAC);
	      for(int ibuf = 0; ibuf < ibuf_loop; ibuf++){
		//---------------------------------------------------------
		cnum = 0;
		if( (lane + ibuf * TSUB_MAC) < Nsweep ){
		  //-------------------------------------------------------
		  /* load a tree node corresponding the tree cell */
		  more_tmp = node[list0[hbuf + lane + ibuf * TSUB_MAC]];
		  const int nodenum  = 1 + (more_tmp >> IDXBITS);
		  const int nodehead =      more_tmp  & IDXMASK;
		  //-------------------------------------------------------
		  /* load all child nodes of the tree cell */
		  more_tmp = more[nodehead];
		  cnum  = 1 + (more_tmp >> IDXBITS);
		  chead =      more_tmp  & IDXMASK;
		  for(int jj = 1; jj < nodenum; jj++)
		    cnum += (1 + (more[nodehead + jj] >> IDXBITS));
		  //-------------------------------------------------------
		}/* if( (lane + ibuf * TSUB_MAC) < Nsweep ){ */
		//---------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		smem = prefixSumTsub(cnum, lane);
		const int lend = BLOCKSIZE(__shfl(smem, TSUB_MAC - 1, TSUB_MAC), NBUF_MAC * TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		prefixSumTsub(cnum, smem, tidx, lane);
		const int lend = BLOCKSIZE(       smem[tail].i,                  NBUF_MAC * TSUB_MAC);
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		for(int ll = 0; ll < lend; ll++){
		  //-------------------------------------------------------
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
		  cnum -= unum;
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  const int Ntmp = smem         - (NBUF_MAC * TSUB_MAC);/* Ntmp is a temporal buffer */
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  const int Ntmp = smem[tidx].i - (NBUF_MAC * TSUB_MAC);/* Ntmp is a temporal buffer */
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		  //-------------------------------------------------------

		  //-------------------------------------------------------
		  /* pick up candidate tree nodes */
		  //-------------------------------------------------------
#   if  NBUF_MAC == 4
		  alignedFlt rjmax_loc = { REAL_MIN,  REAL_MIN,  REAL_MIN,  REAL_MIN};
		  alignedInt pjidx_loc = {NULL_NODE, NULL_NODE, NULL_NODE, NULL_NODE};
#endif//NBUF_MAC == 4
#   if  NBUF_MAC == 2
		  alignedFlt rjmax_loc = { REAL_MIN,  REAL_MIN};
		  alignedInt pjidx_loc = {NULL_NODE, NULL_NODE};
#endif//NBUF_MAC == 2
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  const int stail = (__shfl(smem, TSUB_MAC - 1, TSUB_MAC) < (TSUB_MAC * NBUF_MAC)) ? (__shfl(smem, TSUB_MAC - 1, TSUB_MAC)) : (TSUB_MAC * NBUF_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  const int stail = (       smem[tail].i                  < (TSUB_MAC * NBUF_MAC)) ? (       smem[tail].i                 ) : (TSUB_MAC * NBUF_MAC);
#endif//USE_WARP_SHUFFLE_FUNC_MAC
#if 1
#pragma unroll
		  for(int kk = 0; kk < NBUF_MAC; kk++){
		    //-----------------------------------------------------
		    const int jj = lane + kk * TSUB_MAC;
		    if( jj >= stail )
		      break;
		    //-----------------------------------------------------
		    const int kidx = pjidx[hbuf + jj];
		    const jparticle jpos = pj[kidx];
		    //-----------------------------------------------------
		    const real dx = jpos.x - jcom.x;
		    const real dy = jpos.y - jcom.y;
		    const real dz = jpos.z - jcom.z;
		    const real d2 = 1.0e-30f + dx * dx + dy * dy + dz * dz;
		    const real dr = d2 * RSQRT(d2);
		    //-----------------------------------------------------
		    const real rjmax = bmax[kidx] + dr;
		    const real rjmin = -rjmax + (TWO * (UNITY - EPSILON)) * dr;
		    rmin = (rjmin > rmin) ? rjmin : rmin;
		    //-----------------------------------------------------
		    if( rjmax > rmin ){
		      pjidx_loc.ia[kk] = kidx;
		      rjmax_loc.ra[kk] = rjmax;
		    }/* if( rjmax > rmin ){ */
		    //-----------------------------------------------------
		  }/* for(int kk = 0; kk < NBUF_MAC; kk++){ */
#else
		  for(int jj = lane; jj < stail; jj += TSUB_MAC){
		    //-----------------------------------------------------
		    const int kidx = pjidx[hbuf + jj];
		    const jparticle jpos = pj[kidx];
		    //-----------------------------------------------------
		    const real dx = jpos.x - jcom.x;
		    const real dy = jpos.y - jcom.y;
		    const real dz = jpos.z - jcom.z;
		    const real d2 = 1.0e-30f + dx * dx + dy * dy + dz * dz;
		    const real dr = d2 * RSQRT(d2);
		    //-----------------------------------------------------
		    const real rjmax = bmax[kidx] + dr;
		    const real rjmin = -rjmax + (TWO * (UNITY - EPSILON)) * dr;
		    rmin = (rjmin > rmin) ? rjmin : rmin;
		    //-----------------------------------------------------
		    if( rjmax > rmin ){
		      const int itmp = jj / TSUB_MAC;
		      pjidx_loc.ia[itmp] = kidx;
		      rjmax_loc.ra[itmp] = rjmax;
		    }/* if( rjmax > rmin ){ */
		    //-----------------------------------------------------
		  }/* for(int jj = lane; jj < stail; jj += TSUB_MAC){ */
#endif
		  //-------------------------------------------------------

		  //-------------------------------------------------------
		  /* share rmin within TSUB_MAC threads */
		  //-------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  rmin = getMaximumRealTsub(rmin);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  getMaximumRealTsub(&rmin, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		  //-------------------------------------------------------

		  //-------------------------------------------------------
		  /* recheck local buffer (is really rjmax greater than rmin ?) */
		  //-------------------------------------------------------
#pragma unroll
		  for(int jj = 0; jj < NBUF_MAC; jj++){
		    //-----------------------------------------------------
		    const int share = ( rjmax_loc.ra[jj] > rmin ) ? 1 : 0;
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		    smem = prefixSumTsub(share, lane);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		    prefixSumTsub(share, smem, tidx, lane);
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		    //-----------------------------------------------------
		    if( share ){
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		      const int dst = hbuf + Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_MAC
		      const int dst = hbuf + Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		      list1[dst] = pjidx_loc.ia[jj];
		      rjbuf[dst] = rjmax_loc.ra[jj];
		    }/* if( share ){ */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		    Nloc += __shfl(smem, TSUB_MAC - 1, TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		    Nloc +=        smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		    //-----------------------------------------------------
		    if( Nloc > ((NBUF_MAC - 1) * TSUB_MAC) ){
		      //---------------------------------------------------
		      for(int kk = lane; kk < Nloc; kk += TSUB_MAC){
			more1Buf[bufHead + Nbuf + kk] = list1[hbuf + kk];
			rjmaxBuf[bufHead + Nbuf + kk] = rjbuf[hbuf + kk];
		      }/* for(int kk = lane; kk < Nloc; kk += TSUB_MAC){ */
		      //---------------------------------------------------
		      Nbuf += Nloc;
		      Nloc = 0;
		      //---------------------------------------------------
		    }/* if( Nloc > ((NBUF_MAC - 1) * TSUB_MAC) ){ */
		    //-----------------------------------------------------
		  }/* for(int jj = 0; jj < NBUF_MAC; jj++){ */
		  //-------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  smem         = Ntmp;/* Ntmp is a temporal buffer */
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  smem[tidx].i = Ntmp;/* Ntmp is a temporal buffer */
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		  //-------------------------------------------------------

		  //-------------------------------------------------------
		}/* for(int ll = 0; ll < lend; ll++){ */
		//---------------------------------------------------------
	      }/* for(int ibuf = 0; ibuf < NBUF_MAC; ibuf++){ */
	      //-----------------------------------------------------------

	      //-----------------------------------------------------------
	      Ntry -= Nsweep;
	      //-----------------------------------------------------------
	      /* copy data from global memory to shared memory */
	      const int Ncopy = (Ntry < (TSUB_MAC * NBUF_MAC)) ? (Ntry) : (TSUB_MAC * NBUF_MAC);
	      for(int jj = lane; jj < Ncopy; jj += TSUB_MAC)
		list0[hbuf + jj] = more0Buf[(bufHead + (TSUB_MAC * NBUF_MAC) * (iter + 1)) + jj];
	      //-----------------------------------------------------------
	    }/* for(int iter = 0; iter < Niter; iter++){ */
	    //-------------------------------------------------------------


	    //-------------------------------------------------------------
	    if( Nbuf != 0 ){
	      //-----------------------------------------------------------
	      for(int ll = lane; ll < Nloc; ll += TSUB_MAC){
		more1Buf[bufHead + Nbuf + ll] = list1[hbuf + ll];
		rjmaxBuf[bufHead + Nbuf + ll] = rjbuf[hbuf + ll];
	      }/* for(int ll = lane; ll < Nloc; ll += TSUB_MAC){ */
	      //-----------------------------------------------------------
	      for(int ll = lane; ll < TSUB_MAC * NBUF_MAC; ll += TSUB_MAC){
		list1[hbuf + ll] = more1Buf[bufHead + ll];
		rjbuf[hbuf + ll] = rjmaxBuf[bufHead + ll];
	      }/* for(int ll = lane; ll < TSUB_MAC * NBUF_MAC; ll += TSUB_MAC){ */
	      //-----------------------------------------------------------
	    }/* if( Nbuf != 0 ){ */
	    //-------------------------------------------------------------
	    Ntry = Nbuf + Nloc;
	    if( (lane == 0) && (Ntry > NUM_ALLOC_MACBUF) )
	      atomicAdd(overflow, 1);
	    //-------------------------------------------------------------


	    //-------------------------------------------------------------
	    /* list up all child nodes that satisfy rjmax > rmin */
	    //-------------------------------------------------------------
	    inum = 0;
	    //-------------------------------------------------------------
	    Nloc = 0;
	    Nbuf = 0;
	    //-------------------------------------------------------------
	    Niter = BLOCKSIZE(Ntry, NBUF_MAC * TSUB_MAC);
	    for(int iter = 0; iter < Niter; iter++){
	      //-----------------------------------------------------------
	      const int krem = (Ntry < (NBUF_MAC * TSUB_MAC)) ? Ntry : (NBUF_MAC * TSUB_MAC);
	      const int knum = BLOCKSIZE(krem, TSUB_MAC);
	      //-----------------------------------------------------------
	      for(int ki = 0; ki < knum; ki++){
		//---------------------------------------------------------
		int cellIdx = lane + ki * TSUB_MAC;
		int  add = 0;
		int iadd = 0;
		//---------------------------------------------------------

		//---------------------------------------------------------
		/* select distant tree cells */
		//---------------------------------------------------------
		if( cellIdx < krem ){
		  //-------------------------------------------------------
		  /* when the current node must be taken into account */
		  cellIdx += hbuf;
		  if( rjbuf[cellIdx] > rmin ){
		    //-----------------------------------------------------
		    /* count up total number of contained i-particles */
		    cellIdx = node2cell[list1[cellIdx]];
		    iadd = cell[cellIdx].num;
		    //-----------------------------------------------------
		    add = 1;
		    //-----------------------------------------------------
		  }/* if( rjbuf[cellIdx] > rmin ){ */
		  //-------------------------------------------------------
		}/* if( cellIdx < krem ){ */
		//---------------------------------------------------------

		//---------------------------------------------------------
		/* remove duplicated tree cells */
		//---------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		smem = prefixSumTsub(add, lane);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		prefixSumTsub(add, smem, tidx, lane);
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		if( add ){
		  //-------------------------------------------------------
		  /* test uploading... */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  const int smidx = Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  const int smidx = Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		  list0[hbuf + smidx] = cellIdx;
		  //-------------------------------------------------------
		  /* if detect duplication, upload flag is turned off */
		  if( ((smidx > 0) && (list0[hbuf + smidx - 1] == cellIdx)) || ((Nbuf > 0) && (smidx == 0) && (more0Buf[bufHead + Nbuf - 1] == cellIdx)) ){
		    add  = 0;
		    iadd = 0;
		  }/* if( ((smidx > 0) && (list0[hbuf + smidx - 1] == cellIdx)) || ((Nbuf > 0) && (smidx == 0) && (more0Buf[bufHead + Nbuf - 1] == cellIdx)) ){ */
		  //-------------------------------------------------------
		}
		//---------------------------------------------------------

		//---------------------------------------------------------
		/* save tree cells on the local buffer */
		//---------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		smem = prefixSumTsub(add, lane);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		prefixSumTsub(add, smem, tidx, lane);
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		if( add ){
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		  const int smidx = Nloc + smem         - 1;
#else///USE_WARP_SHUFFLE_FUNC_MAC
		  const int smidx = Nloc + smem[tidx].i - 1;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		  list0[hbuf + smidx] = cellIdx;
		}
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		Nloc += __shfl(smem, TSUB_MAC - 1, TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		Nloc +=        smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		//---------------------------------------------------------

		//---------------------------------------------------------
		/* move data to the remote buffer if necessary */
		if( Nloc > ((NBUF_MAC - 1) * TSUB_MAC) ){
		  for(int ll = lane; ll < Nloc; ll += TSUB_MAC)
		    more0Buf[bufHead + Nbuf + ll] = list0[hbuf + ll];
		  Nbuf += Nloc;
		  Nloc  = 0;
		}/* if( Nloc > ((NBUF_MAC - 1) * TSUB_MAC) ){ */
		//---------------------------------------------------------
		/* sum up iadd within TSUB_MAC threads */
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
		iadd = accumulateIntTsub(iadd);
#else///USE_WARP_SHUFFLE_FUNC_MAC
		accumulateIntTsub(&iadd, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_MAC
		inum += iadd;
		//---------------------------------------------------------
	      }/* for(int ki = 0; ki < knum; ki++){ */
	      //-----------------------------------------------------------

	      //-----------------------------------------------------------
	      Ntry -= krem;
	      //-----------------------------------------------------------
	      /* copy data from remote buffer to local buffer */
	      const int Ncopy = (Ntry < (TSUB_MAC * NBUF_MAC)) ? (Ntry) : (TSUB_MAC * NBUF_MAC);
	      for(int jj = lane; jj < Ncopy; jj += TSUB_MAC){
		rjbuf[hbuf + jj] = rjmaxBuf[bufHead + NBUF_MAC * TSUB_MAC * (iter + 1) + jj];
		list1[hbuf + jj] = more1Buf[bufHead + NBUF_MAC * TSUB_MAC * (iter + 1) + jj];
	      }/* for(int jj = lane; jj < Ncopy; jj += TSUB_MAC){ */
	      //-----------------------------------------------------------
	    }/* for(int iter = 0; iter < Niter1; iter++){ */
	    //-------------------------------------------------------------


	    //-------------------------------------------------------------
	    if( Nbuf != 0 ){
	      //-----------------------------------------------------------
	      for(int ll = lane; ll < Nloc; ll += TSUB_MAC)
		more0Buf[bufHead + Nbuf + ll] = list0[hbuf + ll];
	      //-----------------------------------------------------------
	      for(int ll = lane; ll < NBUF_MAC * TSUB_MAC; ll += TSUB_MAC)
		list0[hbuf + ll] = more0Buf[bufHead + ll];
	      //-----------------------------------------------------------
	    }/* if( Nbuf != 0 ){ */
	    //-------------------------------------------------------------
	    Ntry = Nloc + Nbuf;
	    if( (lane == 0) && (Ntry > NUM_ALLOC_MACBUF) )
	      atomicAdd(overflow, 1);
	    //-------------------------------------------------------------
	  }/* while( inum > NI_BMAX_ESTIMATE ){ */
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
	  /* check positions of all the pick upped i-particles */
	  real jbmax = ZERO;
	  /* Since NI_BMAX_ESTIMATE <= TSUB_MAC * NBUF_MAC, Ntry1 is less than TSUB_MAC * NBUF_MAC */
	  /* load index of the pick upped i-particles to list1 */
	  const int Niter = BLOCKSIZE(Ntry, TSUB_MAC);
	  int Ncand = 0;
	  for(int iter = 0; iter < Niter; iter++){
	    //-------------------------------------------------------------
	    treecell cand;
	    int pnum = 0;
	    if( lane < Ntry ){
	      cand = cell[list0[hbuf + iter * TSUB_MAC + lane]];
	      pnum = cand.num;
	    }/* if( lane < Ntry ){ */
	    //-------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
	    smem = prefixSumTsub(pnum, lane);
	    for(int jj = 0; jj < pnum; jj++)
	      list1[hbuf + Ncand + smem         - pnum + jj] = cand.head + jj;
	    Ncand += __shfl(smem, TSUB_MAC - 1, TSUB_MAC);
#else///USE_WARP_SHUFFLE_FUNC_MAC
	    prefixSumTsub(pnum, smem, tidx, lane);
	    for(int jj = 0; jj < pnum; jj++)
	      list1[hbuf + Ncand + smem[tidx].i - pnum + jj] = cand.head + jj;
	    Ncand +=        smem[tail].i;
#endif//USE_WARP_SHUFFLE_FUNC_MAC
	    //-------------------------------------------------------------
	    Ntry -= TSUB_MAC;
	    //-------------------------------------------------------------
	  }/* for(int iter = 0; iter < Niter; iter++){ */
	  //---------------------------------------------------------------
	  for(int jj = lane; jj < Ncand; jj += TSUB_MAC){
	    //-------------------------------------------------------------
	    const position ipos = pi[list1[hbuf + jj]];
	    //-------------------------------------------------------------
	    const real dx = ipos.x - jcom.x;
	    const real dy = ipos.y - jcom.y;
	    const real dz = ipos.z - jcom.z;
	    const real r2 = 1.0e-30f + dx * dx + dy * dy + dz * dz;
	    //-------------------------------------------------------------
	    if( r2 > jbmax )
	      jbmax = r2;
	    //-------------------------------------------------------------
	  }/* for(int jj = lane; jj < Ncand; jj += TSUB_MAC){ */
	  //---------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
	  jbmax = getMaximumRealTsub(jbmax);
#else///USE_WARP_SHUFFLE_FUNC_MAC
	  getMaximumRealTsub(&jbmax, smem, tidx, head);
#endif//USE_WARP_SHUFFLE_FUNC_MAC
	  jbmax *= RSQRT(jbmax);
	  if( lane == 0 )
	    bmax[jidx] = jbmax;
	  //---------------------------------------------------------------
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
	  //---------------------------------------------------------------
	  if( lane == 0 )
	    pj[jidx] = jcom;
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
	  local.nodeNum++;
	  local.r2Mean +=          jcom.w;
	  local.r2Sdev += jcom.w * jcom.w;
#endif//COUNT_INTERACTIONS
	  //---------------------------------------------------------------
	}/* if( !leaf[cidx] ){ */
	//-----------------------------------------------------------------
      }/* if( root.head != NULL_CELL ){ */
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < bnum * BLOCKSIZE(clev.num, bnum * NGROUPS_MAC); ii += bnum){ */
    //---------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
    //---------------------------------------------------------------------
    if( lane == 0 ){
      atomicAdd(&(stats[levelIdx].nodeNum), local.nodeNum);
      atomicAdd(&(stats[levelIdx].cellNum), local.cellNum);
      atomicAdd(&(stats[levelIdx].mjMean), local. mjMean);
      atomicAdd(&(stats[levelIdx].mjSdev), local. mjSdev);
      atomicAdd(&(stats[levelIdx].r2Mean), local. r2Mean);
      atomicAdd(&(stats[levelIdx].r2Sdev), local. r2Sdev);
    }/* if( lane == 0 ){ */
    //---------------------------------------------------------------------
#endif//COUNT_INTERACTIONS
    //---------------------------------------------------------------------
  }/* for(int levelIdx = (NUM_PHKEY_LEVEL - 2); levelIdx >= 0; levelIdx--){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  DBG_CALC_MULTIPOLE
//-------------------------------------------------------------------------
__global__ void printNode2Cell_kernel(const int pjNum, READ_ONLY int * RESTRICT node2cell, READ_ONLY uint * RESTRICT more)
{
  //-----------------------------------------------------------------------
  for(int jj = 0; jj < pjNum; jj++)
    printf("jidx = %4d, cidx = %d, more.head = %d, more.num = %d\n", jj, node2cell[jj], more[jj] & IDXMASK, 1 + (more[jj] >> IDXBITS));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__global__ void printMultipole_kernel(const int pjNum, READ_ONLY jparticle * RESTRICT pj, READ_ONLY jmass * RESTRICT mj, READ_ONLY real * RESTRICT bmax)
{
  //-----------------------------------------------------------------------
  for(int jj = 0; jj < pjNum; jj++)
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
    printf("%4d\t%e\t%+e\t%+e\t%+e\t%+e\t%e\n", jj, mj[jj].mass, pj[jj].x, pj[jj].y, pj[jj].z, pj[jj].w, bmax[jj]);
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
    printf("%4d\t%e\t%+e\t%+e\t%+e\t%+e\t%e\n", jj, mj[jj]     , pj[jj].x, pj[jj].y, pj[jj].z, pj[jj].w, bmax[jj]);
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__global__ void printTreeCell_kernel(READ_ONLY PHinfo * RESTRICT level, READ_ONLY treecell * RESTRICT cell, READ_ONLY uint * RESTRICT children, READ_ONLY bool * RESTRICT leaf, READ_ONLY uint * RESTRICT node)
{
  //-----------------------------------------------------------------------
  for(int jj = 0; jj < NUM_PHKEY_LEVEL; jj++)
    printf("jj = %d, level = %u, head = %u, num = %u\n", jj, level[jj].level, level[jj].head, level[jj].num);
  //-----------------------------------------------------------------------
  printf("\n");
  for(int jj = 0; jj < NUM_PHKEY_LEVEL; jj++)
    for(int kk = 0; kk < level[jj].num; kk++){
      const int idx = level[jj].head + kk;
      if( cell[idx].head != NULL_CELL ){
	printf("cidx = %d, np = %d, ph = %d", idx, cell[idx].num, cell[idx].head);
	printf(", num = %d, head = %d", 1 + (node[idx] >> IDXBITS), node[idx] & IDXMASK);
	if( !leaf[idx] )
	  printf(", leaf = %u, nh = %d, ch = %d", leaf[idx], 1 + (children[idx] >> IDXBITS), children[idx] & IDXMASK);
	printf("\n");
      }
    }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//DBG_CALC_MULTIPOLE
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
//-------------------------------------------------------------------------
/* initialize count of Nj and Nbuf */
//-------------------------------------------------------------------------
__global__ void initStatistics_kernel(tree_stats *stats)
{
  //-----------------------------------------------------------------------
  const tree_stats zero_stats = {ZERO, ZERO, ZERO, ZERO, 0, 0};
  //-----------------------------------------------------------------------
  for(int ll = 0; ll < MAXIMUM_PHKEY_LEVEL; ll++)
    stats[ll] = zero_stats;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//COUNT_INTERACTIONS
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* calculate multipole moment(s) of tree cells based on the width-first search */
//-------------------------------------------------------------------------
/* level    :: input          :: head index and number of tree cells contained in the corresponding hierarchy */
/* cell     :: input          :: head index and number of N-body particles contained in the corresponding tree cell */
/* children :: input          :: head index and number of child cells of the corresponding tree cell */
/* leaf     :: input          :: a flag to remember the corresponding tree cell is either leaf(true) of node(false) */
/* pi       :: input          :: position and mass of N-body particles */
/* node     :: input          :: head index and number of pseudo particles of the corresponding tree cell (a leaf cell contains up to NCRIT particles) */
/* pjNum    :: input          :: number of pseudo N-body particles */
/* pj       :: input / output :: position and squared radius of pseudo N-body particle as j-particles */
/* mj       :: input / output :: mass of pseudo N-body particle as j-particles */
//-------------------------------------------------------------------------
extern "C"
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
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
  tree_stats *stats_dev;
  mycudaMalloc((void **)&stats_dev, MAXIMUM_PHKEY_LEVEL * sizeof(tree_stats));
  initStatistics_kernel<<<1, 1>>>(stats_dev);
#endif//COUNT_INTERACTIONS
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* initialize position and mass of pseudo j-particles */
  //-----------------------------------------------------------------------
  int Nrem = BLOCKSIZE(pjNum, NTHREADS_INIT_BODY);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    initTreeBody_kernel<<<Nrem, NTHREADS_INIT_BODY>>>
      (node.jpos, node.mj, node.bmax
#ifdef  WS93_MAC
       , node.mr2
#endif//WS93_MAC
       );
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
      int Nsub = Nblck * NTHREADS_INIT_BODY;
      initTreeBody_kernel<<<Nblck, NTHREADS_INIT_BODY>>>
	(&node.jpos[hidx], &node.mj[hidx], &node.bmax[hidx]
#ifdef  WS93_MAC
	 , &node.mr2[hidx]
#endif//WS93_MAC
	 );
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }/* else{ */
/*   initTreeBody_kernel<<<BLOCKSIZE(pjNum, 1024), 1024>>> */
/*     (node.jpos, node.mj, node.bmax */
/* #ifdef  WS93_MAC */
/*      , node.mr2 */
/* #endif//WS93_MAC */
/*      ); */
  getLastCudaError("initTreeBody_kernel");
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->initTreeBody_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set position and mass of j-particles */
  //-----------------------------------------------------------------------
  Nrem = BLOCKSIZE(piNum, NTHREADS_COPY_BODY);
  if( Nrem <= MAX_BLOCKS_PER_GRID )
    copyRealBody_kernel<<<Nrem, NTHREADS_COPY_BODY>>>
      (piNum, node.jtag,
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
    //---------------------------------------------------------------------
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;
    //---------------------------------------------------------------------
    for(int iter = 0; iter < Niter; iter++){
      //-------------------------------------------------------------------
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;
      //-------------------------------------------------------------------
      int Nsub = Nblck * NTHREADS_COPY_BODY;
      copyRealBody_kernel<<<Nblck, NTHREADS_COPY_BODY>>>
	(piNum, &node.jtag[hidx],
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
      //-------------------------------------------------------------------
      hidx += Nsub;
      Nrem -= Nblck;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < Niter; iter++){ */
    //---------------------------------------------------------------------
  }
  getLastCudaError("copyRealBody_kernel");
  //-----------------------------------------------------------------------
#ifdef  HUNT_MAKE_PARAMETER
  stopStopwatch(&(elapsed->copyRealBody_kernel));
  initStopwatch();
#endif//HUNT_MAKE_PARAMETER
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set pseudo j-particles */
  //-----------------------------------------------------------------------
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
  getLastCudaError("calcMultipole_kernel");
  int fail_hst;
  checkCudaErrors(cudaMemcpy(&fail_hst, buf.fail, sizeof(int), cudaMemcpyDeviceToHost));
  if( fail_hst != 0 ){
    __KILL__(stderr, "ERROR: buffer (%d elements per %d threads group) overflow at least %d times.\nPLEASE re-simulate after increasing NUM_ALLOC_MACBUF defined in src/tree/make.h.\n", NUM_ALLOC_MACBUF, TSUB_MAC, fail_hst);
  }
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->calcMultipole));
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
  //-----------------------------------------------------------------------
#ifdef  MAKE_TREE_ON_DEVICE
  checkCudaErrors(cudaMemcpy(cell_hst.level, cell.level, sizeof(PHinfo) * MAXIMUM_PHKEY_LEVEL, cudaMemcpyDeviceToHost));
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaMemcpy(stats_hst, stats_dev, MAXIMUM_PHKEY_LEVEL * sizeof(tree_stats), cudaMemcpyDeviceToHost));
  //-----------------------------------------------------------------------
  mycudaFree(stats_dev);
  //-----------------------------------------------------------------------
#endif//COUNT_INTERACTIONS
  //-----------------------------------------------------------------------
#   if  !defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
  checkCudaErrors(cudaMemcpy(node_hst.jpos, node.jpos, sizeof(jparticle) * pjNum, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(node_hst.mj  , node.mj  , sizeof(jmass    ) * pjNum, cudaMemcpyDeviceToHost));
#endif//!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
#   if  !defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
  checkCudaErrors(cudaMemcpy(bmax_root_hst, node.bmax, sizeof(     real) *     1, cudaMemcpyDeviceToHost));
#endif//!defined(SERIALIZED_EXECUTION) && !defined(BUILD_LET_ON_DEVICE)
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  DBG_CALC_MULTIPOLE
  printMultipole_kernel<<<1, 1>>>(pjNum, node.jpos, node.mj, node.bmax);
  getLastCudaError("printMultipole_kernel");
  checkCudaErrors(cudaDeviceSynchronize());
  exit(0);
#endif//DBG_CALC_MULTIPOLE
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#if 0
  jparticle *pj_hst;
  mycudaMallocHost((void **)&pj_hst, pjNum * sizeof(jparticle));

  checkCudaErrors(cudaMemcpy(pj_hst, node.jpos, pjNum * sizeof(jparticle), cudaMemcpyDeviceToHost));

  for(int jj = 0; jj < pjNum; jj++)
    fprintf(stdout, "%d:\t%e, %e, %e\n", jj, pj_hst[jj].x, pj_hst[jj].y, pj_hst[jj].z);
  fflush(stdout);
  mycudaFreeHost(pj_hst);
  exit(0);
#endif
  //-----------------------------------------------------------------------
#if 0
  printf("%s(%d): %s: tentative kill for performance tuning\n", __FILE__, __LINE__, __func__);
  exit(0);
#endif
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaFuncSetCacheConfig(calcMultipole_kernel, cudaFuncCachePreferShared));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* error checking before running the kernel */
  //-----------------------------------------------------------------------
  struct cudaFuncAttributes funcAttr;
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, calcMultipole_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAC ){
    //---------------------------------------------------------------------
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in calcMultipole_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_MAC);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_MAC = %d.\n", GPUGEN, GPUVER, NTHREADS_MAC);
#ifdef  USE_WARP_SHUFFLE_FUNC_MAC
    fprintf(stderr, "note: warp shuffle instruction is  enabled\n");
#else///USE_WARP_SHUFFLE_FUNC_MAC
    fprintf(stderr, "note: warp shuffle instruction is disabled\n");
#endif//USE_WARP_SHUFFLE_FUNC_MAC
    fflush (stderr);
    //---------------------------------------------------------------------
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAC ){ */
  int regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS_MAC);
  if( regLimit > (MAX_REGISTERS_PER_SM / NTHREADS_MAC) )
    regLimit = (MAX_REGISTERS_PER_SM / NTHREADS_MAC);
  int memLimit = (48 * 1024) / funcAttr.sharedSizeBytes;
  int Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / NTHREADS_MAC) )    Nblck = MAX_THREADS_PER_SM / NTHREADS_MAC;
  if( Nblck >   MAX_BLOCKS_PER_SM                       )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / NTHREADS_MAC) )    Nblck = ((MAX_WARPS_PER_SM * 32) / NTHREADS_MAC);
  if( Nblck != NBLOCKS_PER_SM_MAC ){
    //---------------------------------------------------------------------
    __KILL__(stderr, "ERROR: # of blocks per SM for calcMultipole_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_MAC defined in src/tree/make_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_MAC, funcAttr.numRegs, REGISTERS_PER_THREAD_MAC, GPUGEN, GPUVER);
    //---------------------------------------------------------------------
  }/* if( Nblck != NBLOCKS_PER_SM ){ */
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, makeTree_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAKE_TREE ){
    //---------------------------------------------------------------------
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in makeTree_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_MAKE_TREE);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_MAKE_TREE = %d.\n", GPUGEN, GPUVER, NTHREADS_MAKE_TREE);
    fflush (stderr);
    //---------------------------------------------------------------------
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAC ){ */
  regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS_MAKE_TREE);
  if( regLimit > (MAX_REGISTERS_PER_SM / NTHREADS_MAKE_TREE) )
    regLimit = (MAX_REGISTERS_PER_SM / NTHREADS_MAKE_TREE);
  memLimit = (16 * 1024) / funcAttr.sharedSizeBytes;
  Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / NTHREADS_MAKE_TREE) )    Nblck = MAX_THREADS_PER_SM / NTHREADS_MAKE_TREE;
  if( Nblck >   MAX_BLOCKS_PER_SM                             )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / NTHREADS_MAKE_TREE) )    Nblck = ((MAX_WARPS_PER_SM * 32) / NTHREADS_MAKE_TREE);
  if( Nblck != NBLOCKS_PER_SM_MAKE_TREE ){
    //---------------------------------------------------------------------
    __KILL__(stderr, "ERROR: # of blocks per SM for makeTree_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_MAKE_TREE defined in src/tree/make_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_MAKE_TREE, funcAttr.numRegs, REGISTERS_PER_THREAD_MAKE_TREE, GPUGEN, GPUVER);
    //---------------------------------------------------------------------
  }/* if( Nblck != NBLOCKS_PER_SM ){ */
  //-----------------------------------------------------------------------
  checkCudaErrors(cudaFuncGetAttributes(&funcAttr, linkTree_kernel));
  if( funcAttr.numRegs != REGISTERS_PER_THREAD_LINK_TREE ){
    //---------------------------------------------------------------------
    fprintf(stderr, "%s(%d): %s\n", __FILE__, __LINE__, __func__);
    fprintf(stderr, "warning: # of registers used (%d) in linkTree_kernel() is not match with the predicted value (%d).\n", funcAttr.numRegs, REGISTERS_PER_THREAD_LINK_TREE);
    fprintf(stderr, "note: GPUGEN = %d, GPUVER = %d, NTHREADS_LINK_TREE = %d.\n", GPUGEN, GPUVER, NTHREADS_LINK_TREE);
    fflush (stderr);
    //---------------------------------------------------------------------
  }/* if( funcAttr.numRegs != REGISTERS_PER_THREAD_MAC ){ */
  regLimit = MAX_REGISTERS_PER_SM / (funcAttr.numRegs * NTHREADS_LINK_TREE);
  if( regLimit > (MAX_REGISTERS_PER_SM / NTHREADS_LINK_TREE) )
    regLimit = (MAX_REGISTERS_PER_SM / NTHREADS_LINK_TREE);
  memLimit = (16 * 1024) / funcAttr.sharedSizeBytes;
  Nblck = (regLimit <= memLimit) ? regLimit : memLimit;
  if( Nblck > (MAX_THREADS_PER_SM       / NTHREADS_LINK_TREE) )    Nblck = MAX_THREADS_PER_SM / NTHREADS_LINK_TREE;
  if( Nblck >   MAX_BLOCKS_PER_SM                             )    Nblck = MAX_BLOCKS_PER_SM;
  if( Nblck > (( MAX_WARPS_PER_SM * 32) / NTHREADS_LINK_TREE) )    Nblck = ((MAX_WARPS_PER_SM * 32) / NTHREADS_LINK_TREE);
  if( Nblck != NBLOCKS_PER_SM_LINK_TREE ){
    //---------------------------------------------------------------------
    __KILL__(stderr, "ERROR: # of blocks per SM for linkTree_kernel() is mispredicted (%d).\n\tThe limits come from register and shared memory are %d and %d, respectively.\n\tHowever, the expected value of NBLOCKS_PER_SM_LINK_TREE defined in src/tree/make_dev.cu is %d.\n\tAdditional information: # of registers per thread is %d while predicted as %d (GPUGEN = %d, GPUVER = %d).\n", Nblck, regLimit, memLimit, NBLOCKS_PER_SM_LINK_TREE, funcAttr.numRegs, REGISTERS_PER_THREAD_LINK_TREE, GPUGEN, GPUVER);
    //---------------------------------------------------------------------
  }/* if( Nblck != NBLOCKS_PER_SM ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//CALC_MULTIPOLE_ON_DEVICE
//-------------------------------------------------------------------------
