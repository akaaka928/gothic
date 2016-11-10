/*************************************************************************\
 *                                                                       *
                  last updated on 2016/10/28(Fri) 16:33:46
 *                                                                       *
 *    Header File for constructing octree structure                      *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef MAKE_H
#define MAKE_H
//-------------------------------------------------------------------------
#include <stdbool.h>
//-------------------------------------------------------------------------
#include <macro.h>
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "../sort/peano.h"
//-------------------------------------------------------------------------
#include "../tree/macutil.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define SPLIT_CHILD_CELL
//-------------------------------------------------------------------------
/* #define INDIVIDUAL_GRAVITATIONAL_SOFTENING */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define TREE_SAFETY_VAL (2.0f) */
/* #define TREE_SAFETY_VAL (1.5f) */
#define TREE_SAFETY_VAL (1.0f)
//-------------------------------------------------------------------------
/* #define NUM_CELL_MAX (NUM_BODY_MAX *      TREE_SAFETY_VAL ) */
/* #define NUM_NODE_MAX (NUM_BODY_MAX * (1 + TREE_SAFETY_VAL)) */
#define NUM_CELL_MAX ((int)ceilf((float)NUM_BODY_MAX *         TREE_SAFETY_VAL ))
#define NUM_NODE_MAX ((int)ceilf((float)NUM_BODY_MAX * (1.0f + TREE_SAFETY_VAL)))
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef       NLEAF
/* #      define NLEAF (16) */
#      define NLEAF (64)
/* #      define NLEAF (128) */
#endif//      NLEAF
/* ymiki@augustus:~/tmp/141001idxmask$ ./a.out  */
/* NLEAF =   8: 0x1fffffff is 536870911 = 512M */
/* NLEAF =  16: 0x0fffffff is 268435455 = 256M */
/* NLEAF =  32: 0x07ffffff is 134217727 = 128M */
/* NLEAF =  64: 0x03ffffff is  67108863 =  64M */
/* NLEAF = 128: 0x01ffffff is  33554431 =  32M */
/* NLEAF = 256: 0x00ffffff is  16777215 =  16M */
//-------------------------------------------------------------------------
/* minimum of NLEAF is set to be 8 */
/* (32 - IDXBITS) contains # of child pseudo-particles of the corresponding pseudo particle (1, 2, 3, ..., NLEAF) */
/* MEMO: in case of NLEAF = 8, then 000 -> 1, 001 -> 2, 010 -> 3, 011 -> 4, 100 -> 5, 101 -> 6, 110 -> 7, 111 -> 8 */
#if     NLEAF <=  8
# undef NLEAF
#define NLEAF    (8)
#                   define IDXBITS (29)
#                   define IDXMASK (0x1fffffff)
#else///NLEAF <=  8
#if     NLEAF <= 16
#                   define IDXBITS (28)
#                   define IDXMASK (0xfffffff)
#else///NLEAF <= 16
#if     NLEAF <= 32
#                   define IDXBITS (27)
#                   define IDXMASK (0x7ffffff)
#else///NLEAF <= 32
#if     NLEAF <= 64
#                   define IDXBITS (26)
#                   define IDXMASK (0x3ffffff)
#else///NLEAF <= 64
#if     NLEAF <= 128
#                   define IDXBITS (25)
#                   define IDXMASK (0x1ffffff)
#else///NLEAF <= 128
# undef NLEAF
#define NLEAF   (256)
#                   define IDXBITS (24)
#                   define IDXMASK (0xffffff)
#endif//NLEAF <= 128
#endif//NLEAF <= 64
#endif//NLEAF <= 32
#endif//NLEAF <= 16
#endif//NLEAF <=  8
//-------------------------------------------------------------------------
#define NULL_CELL (IDXMASK - 1)
#define NULL_NODE (IDXMASK)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef       NCRIT
#      define NCRIT (8)
#endif//      NCRIT
//-------------------------------------------------------------------------
#   if  NCRIT > NLEAF
#undef  NCRIT
#define NCRIT   NLEAF
#endif//NCRIT > NLEAF
//-------------------------------------------------------------------------
#define CELL_UNIT (8)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define NUM_ALLOC_TREE_CELL ((((IDXMASK) + 1) < (NUM_CELL_MAX)) ? ((IDXMASK) + 1) : (NUM_CELL_MAX))
#define NUM_ALLOC_TREE_NODE ((((IDXMASK) + 1) < (NUM_NODE_MAX)) ? ((IDXMASK) + 1) : (NUM_NODE_MAX))
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define NI_BMAX_ESTIMATE (32)
#define NJ_BMAX_ESTIMATE (NLEAF)
//-------------------------------------------------------------------------
#   if  NJ_BMAX_ESTIMATE < NLEAF
#undef  NJ_BMAX_ESTIMATE
#define NJ_BMAX_ESTIMATE   NLEAF
#endif//NJ_BMAX_ESTIMATE < NLEAF
//-------------------------------------------------------------------------
#define NUM_ALLOC_MACBUF (16384)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* definition of structure: fundamental structures */
//-------------------------------------------------------------------------
typedef struct __align__(8)
{
  int num, head;/* the number of N-body particles and the index of the head particle */
} treecell;
//-------------------------------------------------------------------------
static const treecell null_cell = {0, NULL_CELL};
//-------------------------------------------------------------------------
typedef struct __align__(16)
{
  real y, x, w, z;
} jparticle;
//-------------------------------------------------------------------------
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
typedef struct __align__(8)
{
  real mass, eps2;
} jmass;
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
typedef real jmass;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
//-------------------------------------------------------------------------
typedef struct __align__(8)
{
  int head, num;
} laneinfo;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* definition of structure: list of data set */
//-------------------------------------------------------------------------
typedef struct
{
  PHinfo *level;
  treecell *cell;
  bool *leaf;
  uint *ptag;
  PHint *hkey;
  uint *parent, *children;
} soaTreeCell;
//-------------------------------------------------------------------------
typedef struct
{
  uint *more;
  int *jtag;
  int *node2cell;
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  int *niSub;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  jparticle *jpos;
  jmass *mj;
  real *bmax;
#ifdef  WS93_MAC
  real *mr2;
#endif//WS93_MAC
#ifdef  GADGET_MAC
  real *mac;
#endif//GADGET_MAC
/* #ifdef  HUNT_OPTIMAL_SEPARATION */
/*   real r_cutoff; */
/* #endif//HUNT_OPTIMAL_SEPARATION */
} soaTreeNode;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "make.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
#   if  !defined(CALC_MULTIPOLE_ON_DEVICE) && defined(WS93_MAC)
  void setGlobalConstants_make_c(const real invDelta_hst);
#endif//!defined(CALC_MULTIPOLE_ON_DEVICE) && defined(WS93_MAC)
  //-----------------------------------------------------------------------
#ifndef MAKE_TREE_ON_DEVICE
/*   muse allocTreeCell(PHint **hkey, uint **parent, uint **children */
/* #ifndef CALC_MULTIPOLE_ON_DEVICE */
/* 		     , treecell **cell, bool **leaf, uint **node */
/* #endif//CALC_MULTIPOLE_ON_DEVICE */
/* 		     ); */
  muse allocTreeCell(PHint **hkey, uint **parent, uint **children
#ifndef CALC_MULTIPOLE_ON_DEVICE
		     , treecell **cell, bool **leaf, uint **node
#endif//CALC_MULTIPOLE_ON_DEVICE
		     , soaTreeCell *hst);
  void  freeTreeCell(PHint  *hkey, uint  *parent, uint  *children
#ifndef CALC_MULTIPOLE_ON_DEVICE
		     , treecell  *cell, bool  *leaf, uint  *node
#endif//CALC_MULTIPOLE_ON_DEVICE
		     );
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifndef __CUDACC__
  //-----------------------------------------------------------------------
#ifndef MAKE_TREE_ON_DEVICE
  void makeTree
  (PHinfo * restrict level,
   int *num, int *rem, treecell * restrict cell, PHint * restrict hkey, uint * restrict parent, uint * restrict children, bool * restrict leaf,
   const int piNum, PHint * restrict peano,
   int *jnum, uint * restrict ptag, uint * restrict more, int * restrict jtag, int * restrict node2cell
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
   , int * restrict niSub
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#ifdef  EXEC_BENCHMARK
   , wall_clock_time * restrict elapsed
#endif//EXEC_BENCHMARK
   );
#endif//MAKE_TREE_ON_DEVICE
  //-----------------------------------------------------------------------
#ifndef CALC_MULTIPOLE_ON_DEVICE
  void setRealParticles
  (const int Ni, int * restrict jtag, position * restrict pi,
   const int Nj, jparticle * restrict pj, jmass * restrict mj, real * restrict bmax
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
   , const real eps2
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#ifdef  WS93_MAC
   , real * restrict mr2
#endif//WS93_MAC
#ifdef EXEC_BENCHMARK
   , wall_clock_time * restrict elapsed
#endif//EXEC_BENCHMARK
   );
  void calcMultipole
  (PHinfo * restrict level, treecell * restrict cell, bool * restrict leaf, position * restrict pi,
   uint * restrict node, uint * restrict more, int * restrict node2cell,
   jparticle * restrict pj, jmass * restrict mj, real * restrict bmax,
   int * restrict _more0Buf, int * restrict _more1Buf, real * restrict _rjmaxBuf, int * restrict overflow
#ifdef  WS93_MAC
   , real * restrict mr2
#endif//WS93_MAC
#ifdef  COUNT_INTERACTIONS
   , tree_stats * restrict stats
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
   , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
   );
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
#endif//__CUDACC__
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//MAKE_H
//-------------------------------------------------------------------------
