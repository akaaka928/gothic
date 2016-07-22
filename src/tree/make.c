/*************************************************************************\
 *                                                                       *
                  last updated on 2016/02/07(Sun) 19:17:40
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
#include <math.h>
//-------------------------------------------------------------------------
#include <macro.h>
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "../sort/peano.h"
//-------------------------------------------------------------------------
#include "macutil.h"
#include "make.h"
//-------------------------------------------------------------------------
#   if  !defined(CALC_MULTIPOLE_ON_DEVICE) && defined(WS93_MAC)
static real invDelta;
#endif//!defined(CALC_MULTIPOLE_ON_DEVICE) && defined(WS93_MAC)
//-------------------------------------------------------------------------
/* #define USE_BRENT_SCHEME_FOR_TREE_BUILD */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  !defined(CALC_MULTIPOLE_ON_DEVICE) && defined(WS93_MAC)
//-------------------------------------------------------------------------
void setGlobalConstants_make_c(const real invDelta_hst)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  invDelta = invDelta_hst;
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//!defined(CALC_MULTIPOLE_ON_DEVICE) && defined(WS93_MAC)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MAKE_TREE_ON_DEVICE
//-------------------------------------------------------------------------
/* arrays to store properties of tree cells */
//-------------------------------------------------------------------------
/* muse allocTreeCell(PHint **hkey, uint **parent, uint **children */
/* #ifndef CALC_MULTIPOLE_ON_DEVICE */
/* 		   , treecell **cell, bool **leaf, uint **node */
/* #endif//CALC_MULTIPOLE_ON_DEVICE */
/* 		   ) */
muse allocTreeCell(PHint **hkey, uint **parent, uint **children
#ifndef CALC_MULTIPOLE_ON_DEVICE
		   , treecell **cell, bool **leaf, uint **node
#endif//CALC_MULTIPOLE_ON_DEVICE
		   , soaTreeCell *hst)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const size_t num = NUM_ALLOC_TREE_CELL;
  //-----------------------------------------------------------------------
#ifndef MAKE_TREE_ON_DEVICE
  *    hkey = (   PHint *)malloc(num * sizeof(   PHint));  if( *    hkey == NULL ){    __KILL__(stderr, "ERROR: failure to allocate hkey\n");  }
  alloc.host +=                  num * sizeof(   PHint) ;
  *  parent = (    uint *)malloc(num * sizeof(    uint));  if( *  parent == NULL ){    __KILL__(stderr, "ERROR: failure to allocate parent\n");  }
  alloc.host +=                  num * sizeof(    uint) ;
  *children = (    uint *)malloc(num * sizeof(    uint));  if( *children == NULL ){    __KILL__(stderr, "ERROR: failure to allocate children\n");  }
  alloc.host +=                  num * sizeof(    uint) ;
#endif//MAKE_TREE_ON_DEVICE
#ifndef CALC_MULTIPOLE_ON_DEVICE
  *    cell = (treecell *)malloc(num * sizeof(treecell));  if( *    cell == NULL ){    __KILL__(stderr, "ERROR: failure to allocate cell\n");  }
  alloc.host +=                  num * sizeof(treecell) ;
  *    leaf = (    bool *)malloc(num * sizeof(    bool));  if( *    leaf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate leaf\n");  }
  alloc.host +=                  num * sizeof(    bool) ;
  *    node = (    uint *)malloc(num * sizeof(    uint));  if( *    node == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node\n");  }
  alloc.host +=                  num * sizeof(    uint) ;
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifndef MAKE_TREE_ON_DEVICE
  hst->hkey     = *hkey;
  hst->parent   = *parent;
  hst->children = *children;
#endif//MAKE_TREE_ON_DEVICE
#ifndef CALC_MULTIPOLE_ON_DEVICE
  hst->cell  = *cell;
  hst->leaf  = *leaf;
  host->ptag = *node;
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  freeTreeCell(PHint  *hkey, uint  *parent, uint  *children
#ifndef CALC_MULTIPOLE_ON_DEVICE
		   , treecell  *cell, bool  *leaf, uint  *node
#endif//CALC_MULTIPOLE_ON_DEVICE
		   )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#ifndef MAKE_TREE_ON_DEVICE
  free(hkey);
  free(parent);
  free(children);
#endif//MAKE_TREE_ON_DEVICE
#ifndef CALC_MULTIPOLE_ON_DEVICE
  free(cell);
  free(leaf);
  free(node);
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_BRENT_SCHEME_FOR_TREE_BUILD
//-------------------------------------------------------------------------
/* imin   ::  input :: minimum index of the checked region */
/* imax   ::  input :: maximum index of the checked region */
/* data   ::  input :: checked array */
/* goal   ::  input :: value of the array at the root */
/* aerr   ::  input :: absolute error allowed for the goal */
/* return :: output :: root index of the array */
//-------------------------------------------------------------------------
#define USE_BISECTION_AS_STARTER
/* #define USE_SECANT_METHOD */
//-------------------------------------------------------------------------
static inline real loadPHkey4BrentSearch(PHint *data, const int idx, const PHint goal)
{
  //-----------------------------------------------------------------------
  const PHint key = data[idx];
  //-----------------------------------------------------------------------
  return ((key > goal) ? ((real)(key - goal)) : (-(real)(goal - key)));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline int rootFindBrentIdx(const int imin, const int imax, PHint *data, const PHint goal, const real aerr)
{
  //-----------------------------------------------------------------------
  /* set region of root finder */
  //-----------------------------------------------------------------------
  /* pre :: previous (latest) iterate, closest to the root */
  /* old :: older iterate than pre */
  /* rev :: older iterate satisfies (pre <= root <= rev) or (rev <= root <= pre) */
  int iold = imin                       ;  real xold = (real)iold;  real fold = loadPHkey4BrentSearch(data, iold, goal);
  int irev = imax                       ;  real xrev = (real)irev;  real frev = loadPHkey4BrentSearch(data, irev, goal);
#ifdef  USE_BISECTION_AS_STARTER
  int ipre = imin + ((imax - imin) >> 1);  real xpre = (real)ipre;  real fpre = loadPHkey4BrentSearch(data, ipre, goal);
#else///USE_BISECTION_AS_STARTER
  int ipre = irev                       ;  real xpre =       xrev;  real fpre = frev;
#endif//USE_BISECTION_AS_STARTER
  //-----------------------------------------------------------------------
  if( ((fold > ZERO) && (frev > ZERO)) || ((fold < ZERO) && (frev < ZERO)) ){
    __KILL__(stderr, "ERROR: root must be bracketed in rootFindBrentOpt.\n");
  }/* if( ((fold > ZERO) && (frev > ZERO)) || ((fold < ZERO) && (frev < ZERO)) ) */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* iterate until finding the root */
  //-----------------------------------------------------------------------
  /* region :: possible region contains the root */
  /* width2 :: square of the value region at one or two step(s) past */
  /* real region = xrev - xold; */
  real region = ZERO;
  real width2 = ZERO;
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    /* exception handling to satisfy conditions about xpre, xold, and xrev */
    //---------------------------------------------------------------------
    if( ((fpre > ZERO) && (frev > ZERO)) || ((fpre < ZERO) && (frev < ZERO)) ){
      //-------------------------------------------------------------------
      /* the root must be sandwiched between xpre and xrev */
      //-------------------------------------------------------------------
      irev = iold;
      xrev = xold;
      frev = fold;
      region = xpre - xold;
      width2 = region * region;
      //-------------------------------------------------------------------
    }/* if( ((fpre > ZERO) && (frev > ZERO)) || ((fpre < ZERO) && (frev < ZERO)) ){ */
    //---------------------------------------------------------------------
    if( FABS(frev) < FABS(fpre) ){
      //-------------------------------------------------------------------
      /* xpre must be the closest guess to the root */
      //-------------------------------------------------------------------
      iold = ipre;      ipre = irev;      irev = iold;
      xold = xpre;      xpre = xrev;      xrev = xold;
      fold = fpre;      fpre = frev;      frev = fold;
      //-------------------------------------------------------------------
    }/* if( FABS(frev) < FABS(fpre) ){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* convergence check */
    //---------------------------------------------------------------------
    const real tol = TWO * EPSILON * xpre + HALF * aerr;
    const  int ibin = (irev > ipre) ? ((irev - ipre) >> 1) : -((ipre - irev) >> 1);
    const real xbin = (real)ibin;
    if( (((tol + xbin) * (tol - xbin)) > ZERO) || (fpre == ZERO) )
      return ((fpre <= ZERO) ? ipre : (ipre - 1));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* proceed the iteration */
    //---------------------------------------------------------------------
#ifdef  USE_SECANT_METHOD
    if( (width2 > tol * tol) && (((fold + fpre) * (fold - fpre)) > ZERO) )
#else///USE_SECANT_METHOD
    if( (width2 > tol * tol) && (frev != fpre) )
#endif//USE_SECANT_METHOD
      {
	//-----------------------------------------------------------------
	/* try higher order methods */
	real pp, qq;
	//-----------------------------------------------------------------
#ifdef  USE_SECANT_METHOD
	if( (frev != fpre) && (frev != fold) )
#else///USE_SECANT_METHOD
	if( (fold != fpre) && (frev != fold) )
#endif//USE_SECANT_METHOD
	  {
	    //-------------------------------------------------------------
	    /* use inverse quadratic interpolation */
	    //-------------------------------------------------------------
	    /* rpo :: ratio of fpre to fold */
	    /* rpr :: ratio of fpre to frev */
	    const real rpo = fpre / fold;/* R in the TeX note */
	    const real rpr = fpre / frev;/* S in the TeX note */
	    pp = (UNITY - rpr) * rpo * rpo * (xpre - xold) + (UNITY - rpo) * rpr * rpr * (xrev - xpre);
	    qq = (UNITY - rpo) * (UNITY - rpr) * (rpo - rpr);
	    //-------------------------------------------------------------
	  }
	//-----------------------------------------------------------------
	else{
	  //---------------------------------------------------------------
#ifdef  USE_SECANT_METHOD
	  /* use secant method */
	  pp = -fpre * (xold - xpre);
	  qq =          fold - fpre;
#else///USE_SECANT_METHOD
	  /* use false position method */
	  pp = -fpre * (xrev - xpre);
	  qq =          frev - fpre;
#endif//USE_SECANT_METHOD
	  //---------------------------------------------------------------
	}/* else{ */
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* judge whether the interpolation is acceptable or not */
	//-----------------------------------------------------------------
	/* 1. TWO * FABS(pp / qq) < length ==>> use higher order method */
	/* 2. FABS(pp / qq) < FABS(3 * (xpre - xrev) / 4) ==>> use higher order method */
	real crit2 = THREE * xbin;      crit2 *= crit2;
	crit2 = ((crit2 < QUARTER * width2) ? crit2 : QUARTER * width2);
	crit2 *= qq * qq;
	//-----------------------------------------------------------------
	if( pp * pp < crit2 ){
	  //---------------------------------------------------------------
	  /* accept higher order interpolation */
	  width2 = region * region;
	  region = pp / qq;
	  //---------------------------------------------------------------
	}/* if( (corr * corr) < crit2 ){ */
	//-----------------------------------------------------------------
	else{
	  //---------------------------------------------------------------
	  /* reject higher order interpolation, employ bisection method */
	  region = xbin;
	  width2 = region * region;
	  //---------------------------------------------------------------
	}/* else{ */
	//-----------------------------------------------------------------
      }
    //---------------------------------------------------------------------
    else{
      //-------------------------------------------------------------------
      /* reject higher order interpolation, employ bisection method */
      region = xbin;
      width2 = region * region;
      //-------------------------------------------------------------------
    }/* else{ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* set the new guess */
    //---------------------------------------------------------------------
    iold = ipre;
    xold = xpre;
    fold = fpre;
    //---------------------------------------------------------------------
    xpre += NEARBYINT(((region + tol) * (region - tol) > ZERO) ? region : ((xbin >= ZERO) ? tol : -tol));
    ipre  = (int)xpre;
    fpre  = loadPHkey4BrentSearch(data, ipre, goal);
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_BRENT_SCHEME_FOR_TREE_BUILD
//-------------------------------------------------------------------------
static inline void linkNode
(const treecell root, const bool leaf_cell, uint *ptag, uint *more, int *jtag, int *phead, int *pjNum
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
 , int *niSub
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
)
{
  //-----------------------------------------------------------------------
  const int head = root.head;
  if( head != NULL_CELL ){
    //---------------------------------------------------------------------
    /* commit Morton-key on the tree cell */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* if the cell is a leaf cell, then root.num pseudo particles are set. otherwise, a pseudo particle is set.  */
    const int nadd = (leaf_cell) ? (root.num) : (1);
    if( nadd > NLEAF ){
      __KILL__(stderr, "ERROR: nadd (%d) exceeds NLEAF(%d). Enlarge NLEAF and/or MAXIMUM_PHKEY_LEVEL(%d)\n", nadd, NLEAF, MAXIMUM_PHKEY_LEVEL);
    }/* if( nadd > NLEAF ){ */
    //---------------------------------------------------------------------
    *ptag = ((nadd - 1) << IDXBITS) + (*phead);
    //---------------------------------------------------------------------
    if( leaf_cell ){
      //-------------------------------------------------------------------
      /* when the tree cell is a leaf cell, then more index specifies the cell itself */
      //-------------------------------------------------------------------
      for(int jj = 0; jj < nadd; jj++){
	//-----------------------------------------------------------------
	/* commit the new particle to the j-particle array */
	const int jidx = (*phead) + jj;
	more[jidx] = jidx;
	//-----------------------------------------------------------------
	/* connect an i-particle with the corresponding real j-particle */
	jtag[head + jj] = jidx;
	//-----------------------------------------------------------------
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	niSub[jidx] = 1;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	//-----------------------------------------------------------------
      }/* for(int jj = 0; jj < nadd; jj++){ */
      //-------------------------------------------------------------------
    }/* if( leaf_cell ){ */
    //---------------------------------------------------------------------
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
    else
      niSub[*phead] = root.num;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* increment... */
    //---------------------------------------------------------------------
    *pjNum += nadd;
    *phead += nadd;
    //---------------------------------------------------------------------
  }/* if( head != NULL_CELL ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* make tree structure in a width-first manner */
//-------------------------------------------------------------------------
/* level     ::         output :: head index and number of tree cells contained in the corresponding hierarchy */
/* num       ::         output :: number of   used elements in an array ``cell'' */
/* rem       ::         output :: number of unused elements in an array ``cell'' */
/* cell      ::         output :: head index and number of N-body particles contained in the corresponding tree cell */
/* hkey      ::         output :: minimum value of Peano--Hilbert key could be contained in the corresponding tree cell */
/* parent    ::         output :: index of the parent cell of the corresponding tree cell */
/* children  ::         output :: head index and number of child cells of the corresponding tree cell */
/* leaf      ::         output :: a flag to remember the corresponding tree cell is either leaf(true) of node(false) */
/* piNum     :: input          :: number of N-body particles */
/* key       :: input          :: Peano--Hilbert key of N-body particles */
/* ilevel    ::         output :: Peano--Hilbert key level of the leaf cell of the corresponding N-body particles */
/* jnum      ::         output :: number of tree cells queued in an array ``cell'' to be examined in this stage */
/* ptag      ::         output :: head index and number of pseudo particles of the corresponding tree cell (a leaf cell contains up to NCRIT particles) */
/* more      ::         output :: head index and number of child particles of the corresponding j-particle (a leaf cell contains up to NCRIT particles) */
/* jtag      ::         output :: index of a j-particle corresponds to an i-particle */
/* node2cell ::         output :: index of a tree cell corresponds the pseudo j-particle */
/* niSub     ::         output :: # of contained i-particles within a pseudo j-particle */
//-------------------------------------------------------------------------
void makeTree
(PHinfo * restrict level,
 int *num, int *rem, treecell * restrict cell, PHint * restrict hkey, uint * restrict parent, uint * restrict children, bool * restrict leaf,
 const int piNum, PHint * restrict peano,
 int *jnum, uint * restrict ptag, uint * restrict more, int * restrict jtag, int * restrict node2cell
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
 , int *niSub
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < piNum; ii += 256)
    fprintf(stderr, "%20zu is key for %4d-th particle\n", peano[ii], ii);
  fflush(stderr);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  /* initialization */
  //-----------------------------------------------------------------------
  /* initialize information on hierarchy of Peano--Hilbert key */
  for(int ii = 0; ii < NUM_PHKEY_LEVEL; ii++){
    level[ii].head = NULL_CELL;
    level[ii].num  = 0;
  }/* for(int ii = 0; ii < NUM_PHKEY_LEVEL; ii++){ */
  /* initialize information on tree cells */
  for(int ii = (*num) - 1; ii > 0; ii--){
    cell    [ii] = null_cell;
    hkey    [ii] = -1;
    parent  [ii] = NULL_CELL;
    children[ii] = NULL_CELL;
    leaf    [ii] = true;
    ptag    [ii] = NULL_NODE;
  }/* for(int ii = (*num) - 1; ii > 0; ii--){ */
  *num = 0;
  *rem = NUM_ALLOC_TREE_CELL;
  /* initialize information on pseudo particles */
  for(int jj = (*jnum) - 1; jj >= 0; jj--){
    more     [jj] = NULL_NODE;
    node2cell[jj] = NULL_CELL;
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
    niSub    [jj] = 0;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
  }/* for(int jj = (*jnum) - 1; jj >= 0; jj--){ */
  /* initialize information on relation between i-particle and real j-particle */
  for(int ii = 0; ii < piNum; ii++)
    jtag[ii] = NULL_NODE;
  //-----------------------------------------------------------------------
  /* set a root cell */
  cell  [0].head = 0;
  cell  [0].num  = piNum;
  hkey  [0]      = 0;
  parent[0]      = NULL_CELL;
  *num += 1;
  *rem -= 1;
  //-----------------------------------------------------------------------
  level[0].head = 0;
  level[0].num  = 1;
  //-----------------------------------------------------------------------
  static int phead;/* the head index of arrays to store pseudo j-particles (pj, mj) */
  phead = 0;
  *jnum = 0;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* make tree structure in a width-first manner */
  //-----------------------------------------------------------------------
  for(int levelIdx = 0; levelIdx < MAXIMUM_PHKEY_LEVEL; levelIdx++){
    //---------------------------------------------------------------------
    /* set level of tree cells to examine in this procedure */
    //---------------------------------------------------------------------
    const int cellLev  = level[levelIdx].level;
    const int cellHead = level[levelIdx].head;
    const int cellNum  = level[levelIdx].num;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* set keylevel of Peano-Hilbert key hierarchy about child cell(s) */
    //---------------------------------------------------------------------
    /* PH key level is common within the threads */
    const   int leafLevel = cellLev - 1;
    const PHint leafScale = (PHint)1 << (leafLevel * 3);
#if 0
    fprintf(stderr, "%20zu is the leafScale\n", leafScale);
    fflush(stderr);
#endif
    //---------------------------------------------------------------------
    const   int leafNmax = level[levelIdx + 1].nmax;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* global properties of tree cells belong to the same level hierarchy */
    //---------------------------------------------------------------------
    static int cellTail;  /* equal to cellHead of tree cells belong to the lower level */
    static int     tail;  /* equal to the head index of the write buffer for child cells */
    //---------------------------------------------------------------------
    tail = cellTail = cellHead + cellNum;
    //---------------------------------------------------------------------
    /* if( (tail % CELL_UNIT) != 0 ){ */
    /*   tail += (CELL_UNIT - (cellTail % CELL_UNIT)); */
    /*   cellTail = tail; */
    /*   *num = tail; */
    /* }/\* if( (tail % CELL_UNIT) != 0 ){ *\/ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* examine tree cells */
    //---------------------------------------------------------------------
    for(int ii = 0; ii < cellNum; ii++){
      //-------------------------------------------------------------------
      /* load a tree cell and evaluate either the cell has child cell(s) or no */
      //-------------------------------------------------------------------
      const int cidx = cellHead + ii;
      treecell root = (cidx < cellTail) ? (cell[cidx]) : (null_cell);
      PHint keyHead = (cidx < cellTail) ? (hkey[cidx]) : ((PHint)(-1));
      //-------------------------------------------------------------------
      const int node = (root.num > NCRIT) ? (1) : (0);
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* confirmation of the memory capacity for the queue list */
      //-------------------------------------------------------------------
      if( *rem < node * NLEAF ){
	//-----------------------------------------------------------------
	__KILL__(stderr, "ERROR: allocated vector length for tree-cell arrays are too short, increase TREE_SAFETY_VAL(%f) defined in src/tree/make.h or use more MPI processes!\n", TREE_SAFETY_VAL);
	//-----------------------------------------------------------------
      }/* if( *rem < node * NLEAF ){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* divide the responsible tree cell if the cell is not a leaf cell */
      //-------------------------------------------------------------------
      if( node ){
	//-----------------------------------------------------------------
	/* the responsible tree cell is a node cell */
	leaf[cidx] = false;
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* properties of the responsible tree cell */
	//-----------------------------------------------------------------
	int head = root.head;
	int lone = root.num;
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* properties of the PH-key of child cell(s) */
	//-----------------------------------------------------------------
	PHint khead = keyHead;
	PHint ktail = khead + leafScale;
	//-----------------------------------------------------------------
	int tag = 0;
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* examine CELL_UNIT (= 8) candidate cells */
	//-----------------------------------------------------------------
	for(int jj = 0; jj < CELL_UNIT; jj++){
	  //---------------------------------------------------------------
	  /* when a child cell exists... */
	  if( (lone != 0) && (peano[head] < ktail) ){
	    //-------------------------------------------------------------
	    int ihead = 0;
	    int itail = lone - 1;
	    //-------------------------------------------------------------
	    if( peano[head + itail] >= ktail ){
	      //-----------------------------------------------------------
#ifdef  USE_BRENT_SCHEME_FOR_TREE_BUILD
	      //-----------------------------------------------------------
	      ihead = rootFindBrentIdx(ihead, itail, &peano[head], ktail, UNITY);
	      //-----------------------------------------------------------
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
		if( peano[head + ihalf] <= ktail )		  ihead = (int)ihalf;
		else		                                  itail = (int)ihalf;
		//---------------------------------------------------------
	      }/* while( true ){ */
	      //-----------------------------------------------------------
#endif//USE_BRENT_SCHEME_FOR_TREE_BUILD
	      //-----------------------------------------------------------
	      itail = ihead;
	      //-----------------------------------------------------------
	    }/* if( peano[head + itail] >= ktail ){ */
	    //-------------------------------------------------------------
	    itail++;
	    //-------------------------------------------------------------
#if 0
	    fprintf(stderr, "cellLev = %d, cidx = %d, head = %d, ihead = %d, itail = %d\n", cellLev, cidx, head, ihead, itail);
#endif
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
#ifdef  SPLIT_CHILD_CELL
	    //-------------------------------------------------------------
	    /* split child cell if the number of N-body particles exceeds the upper limit */
	    const int nsplit = (itail < leafNmax) ? (1) : (((itail) + (leafNmax - 1)) / (leafNmax));
	    /* __NOTE__("jj = %d, lone = %d, cellLev = %d, leafNmax = %d, cidx = %d, itail = %d, nsplit = %d\n", jj, lone, cellLev, leafNmax, cidx, itail, nsplit); */
	    /* const uint nchild = itail / nsplit; */
	    const int nchild = (itail + (nsplit - 1)) / nsplit;
	    int nrem = itail;
	    for(int sidx = 0; sidx < nsplit; sidx++){
	      //-----------------------------------------------------------
	      /* estimate the number of particles per child cell */
	      //-----------------------------------------------------------
	      const int nset = IMIN(nchild, nrem);
	      //-----------------------------------------------------------

	      //-----------------------------------------------------------
	      /* commit the new child cell */
	      //-----------------------------------------------------------
	      treecell childCell;
	      childCell.head = head;
	      childCell.num  = nset;
	      //-----------------------------------------------------------
	      cell  [tail + tag] = childCell;
	      hkey  [tail + tag] = khead;
	      parent[tail + tag] = cidx;
	      //-----------------------------------------------------------
	      tag++;
	      head += nset;
	      lone -= nset;
	      nrem -= nset;
	      //-----------------------------------------------------------
	    }/* for(int sidx = 0; sidx < nsplit; sidx++){ */
	    //-------------------------------------------------------------
#else///SPLIT_CHILD_CELL
	    //-------------------------------------------------------------
	    /* commit the new child cell */
	    //-------------------------------------------------------------
	    treecell childCell;
	    childCell.head =  head;
	    childCell.num  = itail;
	    //-------------------------------------------------------------
	    cell  [tail + tag] = childCell;
	    hkey  [tail + tag] = khead;
	    parent[tail + tag] = cidx;
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    /* preparation to go on the next candidate cell */
	    //-------------------------------------------------------------
	    tag++;
	    head += itail;
	    lone -= itail;
	    //-------------------------------------------------------------
#endif//SPLIT_CHILD_CELL
	    //-------------------------------------------------------------
#if 0
	    fprintf(stderr, "cellLev = %d, cidx = %d, id = %d, ihead = %d, inum = %d\n", cellLev, cidx, tail + tag - 1, cell[tail + tag - 1].head, cell[tail + tag - 1].num);
#endif
	    //-------------------------------------------------------------
	  }/* if( (lone != 0) && (peano[head] < ktail) ){ */
	  //---------------------------------------------------------------
	  /* move to the next candidate cell */
	  //---------------------------------------------------------------
	  khead  = ktail;
	  ktail += leafScale;
	  //---------------------------------------------------------------
	}/* for(int jj = 0; jj < CELL_UNIT; jj++){ */
	//-----------------------------------------------------------------
	/* error check */
	//-----------------------------------------------------------------
	if( lone != 0 ){
	  __NOTE__("the cell has %d child-cells\n", tag);
#   if  MAXIMUM_PHKEY_LEVEL <= 10
	  __NOTE__("head of the key of the first child cell is %u\n", keyHead);
	  __NOTE__("tail of the key of the  last child cell is %u\n", ktail);
	  __NOTE__("key scale of child cells times 8        is %u\n", leafScale << 3);
	  __NOTE__("key scale of child cells                is %u\n", leafScale);
#else///MAXIMUM_PHKEY_LEVEL <= 10
	  __NOTE__("head of the key of the first child cell is %zu\n", keyHead);
	  __NOTE__("tail of the key of the  last child cell is %zu\n", ktail);
	  __NOTE__("key scale of child cells times 8        is %zu\n", leafScale << 3);
	  __NOTE__("key scale of child cells                is %zu\n", leafScale);
#endif//MAXIMUM_PHKEY_LEVEL <= 10
	  __KILL__(stderr, "ERROR: tree construction failed. # of left particle(s) is %d(/%d).\n", lone, root.num);
	}/* if( lone != 0 ){ */
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* commit child cells */
	//-----------------------------------------------------------------
#if 1
	if( tag > NLEAF ){
	  __KILL__(stderr, "ERROR: # of child cells (%d) exceeds NLEAF(%d). Enlarge NLEAF and/or MAXIMUM_PHKEY_LEVEL(%d)\n", tag, NLEAF, MAXIMUM_PHKEY_LEVEL);
	}
#endif
	children[cidx] = ((tag - 1) << IDXBITS) + tail;
	//-----------------------------------------------------------------
      	/* const int increment = ((tag + (CELL_UNIT - 1)) / CELL_UNIT) * CELL_UNIT; */
      	const int increment = tag;
      	tail += increment;
      	*num += increment;
      	*rem -= increment;
	//-----------------------------------------------------------------
#if 0
	if( tag != 0 )
	  fprintf(stderr, "cellLev = %d, cidx = %d, tag = %d\n", cellLev, cidx, tag);
#endif
	//-----------------------------------------------------------------
      }/* if( node ){ */
      //-------------------------------------------------------------------


      //-------------------------------------------------------------------
      /* extend the pseudo particle chain */
      //-------------------------------------------------------------------
      linkNode(root, leaf[cidx], &ptag[cidx], more, jtag, &phead, jnum
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	       , niSub
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	       );
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < cellNum; ii++){ */
    //---------------------------------------------------------------------
    level[levelIdx + 1].head =        cellTail;
    level[levelIdx + 1].num  = tail - cellTail;
    //---------------------------------------------------------------------
  }/* for(int levelIdx = 0; levelIdx < MAXIMUM_PHKEY_LEVEL; levelIdx++){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* extend the pseudo particle chain in the lowest level */
  //-----------------------------------------------------------------------
  {
    //---------------------------------------------------------------------
    const int levelIdx = NUM_PHKEY_LEVEL - 1;
    //---------------------------------------------------------------------
    const int cellHead = level[levelIdx].head;
    const int cellNum  = level[levelIdx].num;
    const int cellTail = cellHead + cellNum;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* examine tree cells */
    //---------------------------------------------------------------------
    for(int ii = 0; ii < cellNum; ii++){
      //-------------------------------------------------------------------
      /* load a tree cell and evaluate either the cell has child cell(s) or no */
      //-------------------------------------------------------------------
      const int cidx = cellHead + ii;
      treecell root = (cidx < cellTail) ? (cell[cidx]) : (null_cell);
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* extend the pseudo particle chain */
      //-------------------------------------------------------------------
      linkNode(root, leaf[cidx], &ptag[cidx], more, jtag, &phead, jnum
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	       , niSub
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES) && !defined(FACILE_NEIGHBOR_SEARCH)
	       );
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < cellNum; ii++){ */
  }
  //-----------------------------------------------------------------------
#if 0
  for(int levelIdx = 0; levelIdx < NUM_PHKEY_LEVEL; levelIdx++)
    printf("levelIdx = %d, num = %d\n", levelIdx, level[levelIdx].num);
  fflush(stdout);
  exit(1);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* connect tree nodes in a different PH-level */
  //-----------------------------------------------------------------------
  for(int levelIdx = (NUM_PHKEY_LEVEL - 1); levelIdx >= 0; levelIdx--){
    //---------------------------------------------------------------------
    /* global properties of tree cells belong to the same level hierarchy */
    //---------------------------------------------------------------------
    static int head;/* the head index of the used array "node" */
    static int tail;/* the tail index of the used array "node" */
    //---------------------------------------------------------------------
    head = level[levelIdx].head;
    tail = level[levelIdx].head + level[levelIdx].num;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* connect pseudo particles */
    //---------------------------------------------------------------------
    for(int ii = head; ii < tail; ii++){
      //-------------------------------------------------------------------
      const int cidx = ii;
      //-------------------------------------------------------------------
      if( (cidx < tail) && (cell[cidx].head != NULL_CELL) ){
	//-----------------------------------------------------------------
	const int pidx = ptag[cidx] & IDXMASK;
#ifdef  SPLIT_CHILD_CELL
	const int pnum = 1 + (int)(ptag[cidx] >> IDXBITS);
	for(int jj = 0; jj < pnum; jj++)
	  node2cell[pidx + jj] = cidx;
#else///SPLIT_CHILD_CELL
	node2cell[pidx] = cidx;
#endif//SPLIT_CHILD_CELL
	//-----------------------------------------------------------------
	if( leaf[cidx] == false ){
	  //---------------------------------------------------------------
	  /* count up number of child nodes */
	  //---------------------------------------------------------------
	  const int chead =               children[cidx]  & IDXMASK;
	  const int ctail = chead + (1 + (children[cidx] >> IDXBITS));
	  int childNum  = 0;
	  for(int jj = chead; jj < ctail; jj++)
	    childNum += (1 + (ptag[jj] >> IDXBITS));
	  //---------------------------------------------------------------
#if 1
	  if( childNum > NLEAF ){
	    __KILL__(stderr, "ERROR: childNum(%d) exceeds NLEAF(%d) for cidx = %d @ levelIdx = %d. Enlarge NLEAF and/or MAXIMUM_PHKEY_LEVEL(%d)\n", childNum, NLEAF, cidx, levelIdx, MAXIMUM_PHKEY_LEVEL);
	  }
#endif
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* commit child particles as pseudo particle */
	  //---------------------------------------------------------------
	  more[pidx] = ((childNum - 1) << IDXBITS) + (ptag[chead] & IDXMASK);
	  //---------------------------------------------------------------
	}/* if( leaf[cidx] == false ){ */
	//-----------------------------------------------------------------
      }/* if( (cidx < tail) && (cell[cidx].head != NULL_CELL) ){ */
      //-------------------------------------------------------------------
    }/* for(int ii = head; ii < tail; ii++){ */
    //---------------------------------------------------------------------
  }/* for(int levelIdx = (NUM_PHKEY_LEVEL - 1); levelIdx >= 0; levelIdx--){ */
  //-----------------------------------------------------------------------
#ifdef EXEC_BENCHMARK
  stopStopwatch(&(elapsed->makeTree));
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
#ifndef CALC_MULTIPOLE_ON_DEVICE
//-------------------------------------------------------------------------
/* set N-body particles as real j-particles */
//-------------------------------------------------------------------------
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
  /* initialization */
  //-----------------------------------------------------------------------
  const jparticle zero_pj = {ZERO, ZERO, ZERO, ZERO};
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
  const jmass     zero_mj = {ZERO, ZERO};
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
  const jmass     zero_mj =  ZERO;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
  for(int jj = Nj - 1; jj >= 0; jj--){
    pj  [jj] = zero_pj;
    mj  [jj] = zero_mj;
    bmax[jj] = ZERO;
#ifdef  WS93_MAC
    mr2 [jj] = ZERO;
#endif//WS93_MAC
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ni; ii++){
    //---------------------------------------------------------------------
    /* load an N-body particle */
    const position ipos = pi  [ii];
    const int      jidx = jtag[ii];
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
    pj [jidx] = jpos;
    mj [jidx] = mj_tmp;
#ifdef  WS93_MAC
    mr2[jidx] = ipos.m * (jpos.x * jpos.x + jpos.y * jpos.y + jpos.z * jpos.z);
#endif//WS93_MAC
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#ifdef EXEC_BENCHMARK
  stopStopwatch(&(elapsed->calcMultipole));
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


#if 1
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
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
  {
    const tree_stats zero_tree = {ZERO, ZERO, ZERO, ZERO, 0, 0};
    for(int ii = 0; ii < MAXIMUM_PHKEY_LEVEL; ii++)
      stats[ii] = zero_tree;
  }
#endif//COUNT_INTERACTIONS
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  /* settings about remote buffers */
  int  *more0Buf = &_more0Buf[0];
  int  *more1Buf = &_more1Buf[0];
  real *rjmaxBuf = &_rjmaxBuf[0];
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate multipole moment(s) of pseudo j-particles */
  //-----------------------------------------------------------------------
  for(int levelIdx = (NUM_PHKEY_LEVEL - 2); levelIdx >= 0; levelIdx--){
    //---------------------------------------------------------------------
    /* commit information on tree cells belonging a same level of Peano--Hilbert key */
    //---------------------------------------------------------------------
    const int cellHead = level[levelIdx].head;
    const int cellNum  = level[levelIdx].num;
    const int cellTail = cellHead + cellNum;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* examine tree cells */
    //---------------------------------------------------------------------
    for(int ii = 0; ii < cellNum; ii++){
      //-------------------------------------------------------------------
      /* load a tree cell and evaluate either the cell is a node or a leaf */
      //-------------------------------------------------------------------
      /* if cidx is not less than tail, then the thread has a null leaf-cell */
      const      int cidx = cellHead + ii;
      const treecell root = (cidx < cellTail) ? (cell[cidx]) : (null_cell);
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* extend the pseudo particle chain */
      //-------------------------------------------------------------------
      if( root.head != NULL_CELL ){
	//-----------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
	if( levelIdx != (NUM_PHKEY_LEVEL - 1) )
	  stats[levelIdx].cellNum++;
#endif//COUNT_INTERACTIONS
	//-----------------------------------------------------------------
	if( !leaf[cidx] ){
	  //---------------------------------------------------------------
	  /* when the tree cell is a node cell, then calculate multipole moment(s) of the cell */
	  //---------------------------------------------------------------
	  jparticle jcom = {ZERO, ZERO, ZERO, ZERO};
	  real      mtot =  ZERO;
	  //---------------------------------------------------------------
#ifdef  WS93_MAC
	  real mjrj2 = ZERO;
#endif//WS93_MAC
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* sum up multipole moment(s) of child cells */
	  //---------------------------------------------------------------
	  const int jidx  = node[cidx] & IDXMASK;
	  const int nhead =      more[jidx]  & IDXMASK;
	  const int nnum  = 1 + (more[jidx] >> IDXBITS);
	  //---------------------------------------------------------------
	  for(int jj = nhead; jj < nhead + nnum; jj++){
	    //-------------------------------------------------------------
	    /* load a pseudo particle corresponding a child cell */
	    //-------------------------------------------------------------
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const real      mass = mj[jj].mass;
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const real      mass = mj[jj];
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
	    const jparticle jpos = pj[jj];
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    /* calculate multipole moments of the pseudo particle group */
	    //-------------------------------------------------------------
	    /* calculate total mass */
	    mtot   += mass;
	    /* calculate center-of-mass */
	    jcom.x += mass * jpos.x;
	    jcom.y += mass * jpos.y;
	    jcom.z += mass * jpos.z;
	    /* calculate trace of quadrupole moment */
#ifdef  WS93_MAC
	    mjrj2 += mr2[jj];
#endif//WS93_MAC
	    //-------------------------------------------------------------
	  }
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* calculate multipole moments of the pseudo particle group */
	  //---------------------------------------------------------------
	  /* calculate center-of-mass */
	  real minv = UNITY / mtot;
	  jcom.x *= minv;
	  jcom.y *= minv;
	  jcom.z *= minv;
	  /* calculate trace of quadrupole moment */
#ifdef  WS93_MAC
	  const real B2 = mjrj2 - mtot * (jcom.x * jcom.x + jcom.y * jcom.y + jcom.z * jcom.z);
#endif//WS93_MAC
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
	  /* estimate size of particle distribution */
	  //---------------------------------------------------------------
	  /* initialize list of examined tree nodes and related variables */
	  int inum = root.num;
	  real rjmax_tmp[NJ_BMAX_ESTIMATE];
	  int list0[NJ_BMAX_ESTIMATE], list1[NJ_BMAX_ESTIMATE];
	  int Ntry0 = 1;	  list0[0] = jidx;
	  //---------------------------------------------------------------
	  /* pick up NI_BMAX_ESTIMATE i-particles in maximum to estimate bmax */
	  while( inum > NI_BMAX_ESTIMATE ){
	    //-------------------------------------------------------------
	    /* find rmax and the corresponding rmin */
	    //-------------------------------------------------------------
	    real rmax = REAL_MIN;
	    real rmin = ZERO;
	    //-------------------------------------------------------------
	    int Nloc = 0;
	    int Nbuf = 0;
	    //-------------------------------------------------------------
	    const int Niter0 = (Ntry0 + NJ_BMAX_ESTIMATE - 1) / NJ_BMAX_ESTIMATE;
	    for(int iter = 0; iter < Niter0; iter++){
	      //-----------------------------------------------------------
	      const int Nloop = (Ntry0 < NJ_BMAX_ESTIMATE) ? Ntry0 : NJ_BMAX_ESTIMATE;
	      //-----------------------------------------------------------
	      for(int jj = 0; jj < Nloop; jj++){
		//---------------------------------------------------------
		const uint more_tmp = more[list0[jj]];
		const int khead =      more_tmp  & IDXMASK;
		const int knum  = 1 + (more_tmp >> IDXBITS);
		//-----------------------------------------------------------
		for(int kk = 0; kk < knum; kk++){
		  //---------------------------------------------------------
		  const int kidx = khead + kk;
		  const jparticle jpos = pj[kidx];
		  //---------------------------------------------------------
		  const real dx = jpos.x - jcom.x;
		  const real dy = jpos.y - jcom.y;
		  const real dz = jpos.z - jcom.z;
		  /* const real disp = SQRT(EPSILON + dx * dx + dy * dy + dz * dz); */
		  const real disp = SQRT(1.0e-30f + dx * dx + dy * dy + dz * dz);
		  const real rjmax = disp + bmax[kidx];
		  //---------------------------------------------------------
		  /* const real rjmin = -rjmax + TWO * disp; */
		  /* const real rjmin = (disp - rjmax) + 0.99999f * disp; */
		  const real rjmin = (disp - rjmax) + (UNITY - EPSILON) * disp;
		  rmin = (rjmin > rmin) ? rjmin : rmin;
		  //---------------------------------------------------------

		  //---------------------------------------------------------
		  /* if the minimum required condition is satisfied */
		  if( rjmax > rmin ){
		    //-------------------------------------------------------
		    /* save on the local buffer */
		    rjmax_tmp[Nloc] = rjmax;
		    list1    [Nloc] = kidx;
		    Nloc++;
		    //-------------------------------------------------------
		    /* move data to the remote buffer if necessary */
		    if( Nloc >= NJ_BMAX_ESTIMATE ){
		      for(int ll = 0; ll < NJ_BMAX_ESTIMATE; ll++){
			rjmaxBuf[Nbuf + ll] = rjmax_tmp[ll];
			more1Buf[Nbuf + ll] = list1    [ll];
		      }
		      Nloc -= NJ_BMAX_ESTIMATE;
		      Nbuf += NJ_BMAX_ESTIMATE;
		    }
		    //-------------------------------------------------------
		    /* when the current node is the most distant node */
		    if( rjmax > rmax ){
		      //-----------------------------------------------------
		      rmax = rjmax;
		      /* const real rjmin = -rjmax + TWO * disp; */
		      /* rmin = (rjmin > ZERO) ? rjmin : ZERO; */
		      //-----------------------------------------------------
		    }
		    //-------------------------------------------------------
		  }
		  //---------------------------------------------------------
		}
		//-----------------------------------------------------------
	      }
	      //-------------------------------------------------------------

	      //-----------------------------------------------------------
	      Ntry0 -= Nloop;
	      //-----------------------------------------------------------
	      /* copy data from remote buffer to local buffer */
	      if( Ntry0 > 0 )
		for(int jj = 0; jj < NJ_BMAX_ESTIMATE; jj++)
		  list0[jj] = more0Buf[NJ_BMAX_ESTIMATE * (iter + 1) + jj];
	      //-----------------------------------------------------------
	    }
	    //-------------------------------------------------------------
	    if( Nbuf != 0 ){
	      //-----------------------------------------------------------
	      for(int ll = 0; ll < Nloc; ll++){
		rjmaxBuf[Nbuf + ll] = rjmax_tmp[ll];
		more1Buf[Nbuf + ll] = list1    [ll];
	      }
	      //-----------------------------------------------------------
	      for(int ll = 0; ll < NJ_BMAX_ESTIMATE; ll++){
		rjmax_tmp[ll] = rjmaxBuf[ll];
		list1    [ll] = more1Buf[ll];
	      }
	      //-----------------------------------------------------------
	    }
	    //-------------------------------------------------------------


	    //-------------------------------------------------------------
	    /* list up all child nodes that satisfy rjmax > rmin */
	    //-------------------------------------------------------------
	    int Ntry1 = Nbuf + Nloc;
	    if( Ntry1 > NUM_ALLOC_MACBUF )	      *overflow += 1;
	    Nloc = 0;
	    Nbuf = 0;
	    inum = 0;
	    //-------------------------------------------------------------
	    const int Niter1 = (Ntry1 + NJ_BMAX_ESTIMATE - 1) / NJ_BMAX_ESTIMATE;
	    for(int iter = 0; iter < Niter1; iter++){
	      //-----------------------------------------------------------
	      const int Nloop = (Ntry1 < NJ_BMAX_ESTIMATE) ? Ntry1 : NJ_BMAX_ESTIMATE;
	      //-----------------------------------------------------------
	      for(int jj = 0; jj < Nloop; jj++){
		//---------------------------------------------------------
		const int   kidx = list1    [jj];
		const real rjmax = rjmax_tmp[jj];
		//---------------------------------------------------------
		/* when the current node must be taken into account */
		if( rjmax > rmin ){
		  //-------------------------------------------------------
		  /* count up total number of contained i-particles */
		  inum += cell[node2cell[kidx]].num;
		  //-------------------------------------------------------
		  /* save on the local buffer */
		  list0[Nloc] = kidx;
		  Nloc++;
		  //-------------------------------------------------------
		  /* move data to the remote buffer if necessary */
		  if( Nloc >= NJ_BMAX_ESTIMATE ){
		    for(int ll = 0; ll < NJ_BMAX_ESTIMATE; ll++){
		      more0Buf[Nbuf + ll] = list0[ll];
		    }
		    Nloc -= NJ_BMAX_ESTIMATE;
		    Nbuf += NJ_BMAX_ESTIMATE;
		  }
		  //-------------------------------------------------------
		}
		//---------------------------------------------------------
	      }
	      //-----------------------------------------------------------

	      //-----------------------------------------------------------
	      Ntry1 -= Nloop;
	      //-----------------------------------------------------------
	      /* copy data from remote buffer to local buffer */
	      if( Ntry1 > 0 )
		for(int jj = 0; jj < NJ_BMAX_ESTIMATE; jj++){
		  rjmax_tmp[jj] = rjmaxBuf[NJ_BMAX_ESTIMATE * (iter + 1) + jj];
		  list1    [jj] = more1Buf[NJ_BMAX_ESTIMATE * (iter + 1) + jj];
		}
	      //-----------------------------------------------------------
	    }
	    //-------------------------------------------------------------
	    if( Nbuf != 0 ){
	      //-----------------------------------------------------------
	      for(int ll = 0; ll < Nloc            ; ll++)		more0Buf[Nbuf + ll] = list0   [ll];
	      for(int ll = 0; ll < NJ_BMAX_ESTIMATE; ll++)		list0   [       ll] = more0Buf[ll];
	      //-----------------------------------------------------------
	    }
	    //-------------------------------------------------------------
	    Ntry0 = Nloc + Nbuf;
	    if( Ntry0 > NUM_ALLOC_MACBUF )	      *overflow += 1;
	    //-------------------------------------------------------------
	  }
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
	  /* check positions of all the pick upped i-particles */
	  real jbmax = ZERO;
	  //---------------------------------------------------------------
	  const int Niter = (Ntry0 + NJ_BMAX_ESTIMATE - 1) / NJ_BMAX_ESTIMATE;
	  for(int iter = 0; iter < Niter; iter++){
	    //-------------------------------------------------------------
	    const int Nloop = (Ntry0 < NJ_BMAX_ESTIMATE) ? Ntry0 : NJ_BMAX_ESTIMATE;
	    //-------------------------------------------------------------
	    for(int jj = 0; jj < Nloop; jj++){
	      //-----------------------------------------------------------
	      const uint more_tmp = more[list0[jj]];
	      const int khead =      more_tmp  & IDXMASK;
	      const int knum  = 1 + (more_tmp >> IDXBITS);
	      //-----------------------------------------------------------
	      for(int kk = 0; kk < knum; kk++){
		//---------------------------------------------------------
		const treecell candidate = cell[node2cell[khead + kk]];
		const int phead = candidate.head;
		const int ptail = candidate.num + phead;
		//---------------------------------------------------------
		for(int pp = phead; pp < ptail; pp++){
		  //-------------------------------------------------------
		  const position ipos = pi[pp];
		  //-------------------------------------------------------
		  const real dx = ipos.x - jcom.x;
		  const real dy = ipos.y - jcom.y;
		  const real dz = ipos.z - jcom.z;
		  const real r2 = dx * dx + dy * dy + dz * dz;
		  //-------------------------------------------------------
		  if( r2 > jbmax )		    jbmax = r2;
		  //-------------------------------------------------------
		}
		//---------------------------------------------------------
	      }
	      //-----------------------------------------------------------
	    }
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    Ntry0 -= Nloop;
	    //-------------------------------------------------------------
	    /* copy data from remote buffer to local buffer */
	    if( Ntry0 > 0 )
	      for(int jj = 0; jj < NJ_BMAX_ESTIMATE; jj++)
		list0[jj] = more0Buf[NJ_BMAX_ESTIMATE * (iter + 1) + jj];
	    //-------------------------------------------------------------
	  }
	  //---------------------------------------------------------------
	  jbmax = SQRT(jbmax);
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
#ifdef  WS93_MAC
	  const real bmax_2 = HALF * jbmax;
	  jcom.w  = bmax_2 + SQRT(bmax_2 * bmax_2 + SQRT(THREE * invDelta * B2));
	  jcom.w *= jcom.w;
#else///WS93_MAC
	  jcom.w  = jbmax * jbmax;
#endif//WS93_MAC
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
	  /* commit a pseudo particle */
	  //---------------------------------------------------------------
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
	  const jmass mj_loc = {mtot, ZERO};
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
	  const jmass mj_loc =  mtot;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
	  //---------------------------------------------------------------
	  pj  [jidx] = jcom;
	  mj  [jidx] = mj_loc;
#ifdef  WS93_MAC
	  mr2 [jidx] = mjrj2;
#endif//WS93_MAC
	  bmax[jidx] = jbmax;
	  //---------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
	  stats[levelIdx].nodeNum++;
	  stats[levelIdx].mjMean +=        mtot;	  stats[levelIdx].r2Mean +=          jcom.w;
	  stats[levelIdx].mjSdev += mtot * mtot;	  stats[levelIdx].r2Sdev += jcom.w * jcom.w;
#endif//COUNT_INTERACTIONS
	  //---------------------------------------------------------------
	}
	//-----------------------------------------------------------------
      }
      //-------------------------------------------------------------------
    }
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->calcMultipole));
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else
//-------------------------------------------------------------------------
#define SEPARATED_IMPLEMENTATION
void calcMultipole
(PHinfo * restrict level, treecell * restrict cell, bool * restrict leaf, position * restrict pi,
 uint * restrict node, uint * restrict more, int * restrict node2cell,
 jparticle * restrict pj, real * restrict mj, real * restrict bmax,
 int * restrict _more0Buf, int * restrict _more1Buf, real * restrict _rjmaxBuf, int * restrict overflow
 , real * restrict mr2
#ifndef SEPARATED_IMPLEMENTATION
 , const int pjNum
#endif//SEPARATED_IMPLEMENTATION
#ifdef  COUNT_INTERACTIONS
 , tree_stats * restrict stats
#endif//COUNT_INTERACTIONS
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
  {
    const tree_stats zero_tree = {ZERO, ZERO, ZERO, ZERO, 0, 0};
    for(int ii = 0; ii < MAXIMUM_PHKEY_LEVEL; ii++)
      stats[ii] = zero_tree;
  }
#endif//COUNT_INTERACTIONS
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  initStopwatch();
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
#ifndef SEPARATED_IMPLEMENTATION
  /* initialization */
  //-----------------------------------------------------------------------
  const jparticle zero_pj = {ZERO, ZERO, ZERO, ZERO};
  for(int jj = pjNum - 1; jj >= 0; jj--){
    pj  [jj] = zero_pj;
    mj  [jj] = ZERO;
    mr2 [jj] = ZERO;
    bmax[jj] = ZERO;
  }
#endif//SEPARATED_IMPLEMENTATION
  //-----------------------------------------------------------------------
  /* settings about remote buffers */
  int  *more0Buf = &_more0Buf[0];
  int  *more1Buf = &_more1Buf[0];
  real *rjmaxBuf = &_rjmaxBuf[0];
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate multipole moment(s) of pseudo j-particles */
  //-----------------------------------------------------------------------
#ifndef SEPARATED_IMPLEMENTATION
  for(int levelIdx = (NUM_PHKEY_LEVEL - 1); levelIdx >= 0; levelIdx--){
#else///SEPARATED_IMPLEMENTATION
  for(int levelIdx = (NUM_PHKEY_LEVEL - 2); levelIdx >= 0; levelIdx--){
#endif//SEPARATED_IMPLEMENTATION
    //---------------------------------------------------------------------
    /* commit information on tree cells belonging a same level of Peano--Hilbert key */
    //---------------------------------------------------------------------
    const int cellHead = level[levelIdx].head;
    const int cellNum  = level[levelIdx].num;
    const int cellTail = cellHead + cellNum;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* examine tree cells */
    //---------------------------------------------------------------------
    for(int ii = 0; ii < cellNum; ii++){
      //-------------------------------------------------------------------
      /* load a tree cell and evaluate either the cell is a node or a leaf */
      //-------------------------------------------------------------------
      /* if cidx is not less than tail, then the thread has a null leaf-cell */
      const      int cidx = cellHead + ii;
      const treecell root = (cidx < cellTail) ? (cell[cidx]) : (null_cell);
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* extend the pseudo particle chain */
      //-------------------------------------------------------------------
      if( root.head != NULL_CELL ){
	//-----------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
	if( levelIdx != (NUM_PHKEY_LEVEL - 1) )
	  stats[levelIdx].cellNum++;
#endif//COUNT_INTERACTIONS
	//-----------------------------------------------------------------
	if( !leaf[cidx] ){
	  //---------------------------------------------------------------
	  /* when the tree cell is a node cell, then calculate multipole moment(s) of the cell */
	  //---------------------------------------------------------------
	  jparticle jcom = {ZERO, ZERO, ZERO, ZERO};
	  real      mtot =  ZERO;
	  //---------------------------------------------------------------
#ifdef  WS93_MAC
	  real mjrj2 = ZERO;
#endif//WS93_MAC
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* sum up multipole moment(s) of child cells */
	  //---------------------------------------------------------------
	  const int jidx  = node[cidx] & IDXMASK;
	  const int nhead =      more[jidx]  & IDXMASK;
	  const int nnum  = 1 + (more[jidx] >> IDXBITS);
	  //---------------------------------------------------------------
	  for(int jj = nhead; jj < nhead + nnum; jj++){
	    //-------------------------------------------------------------
	    /* load a pseudo particle corresponding a child cell */
	    //-------------------------------------------------------------
	    const real      mass = mj[jj];
	    const jparticle jpos = pj[jj];
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    /* calculate multipole moments of the pseudo particle group */
	    //-------------------------------------------------------------
	    /* calculate total mass */
	    mtot   += mass;
	    /* calculate center-of-mass */
	    jcom.x += mass * jpos.x;
	    jcom.y += mass * jpos.y;
	    jcom.z += mass * jpos.z;
	    /* calculate trace of quadrupole moment */
#ifdef  WS93_MAC
	    mjrj2 += mr2[jj];
#endif//WS93_MAC
	    //-------------------------------------------------------------
	  }
	  //---------------------------------------------------------------

	  //---------------------------------------------------------------
	  /* calculate multipole moments of the pseudo particle group */
	  //---------------------------------------------------------------
	  /* calculate center-of-mass */
	  real minv = UNITY / mtot;
	  jcom.x *= minv;
	  jcom.y *= minv;
	  jcom.z *= minv;
	  /* calculate trace of quadrupole moment */
#ifdef  WS93_MAC
	  const real B2 = mjrj2 - mtot * (jcom.x * jcom.x + jcom.y * jcom.y + jcom.z * jcom.z);
#endif//WS93_MAC
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
	  /* estimate size of particle distribution */
	  //---------------------------------------------------------------
	  /* initialize list of examined tree nodes and related variables */
	  int inum = root.num;
	  real rjmax_tmp[NJ_BMAX_ESTIMATE];
	  int list0[NJ_BMAX_ESTIMATE], list1[NJ_BMAX_ESTIMATE];
	  int Ntry0 = 1;	  list0[0] = jidx;
	  //---------------------------------------------------------------
	  /* pick up NI_BMAX_ESTIMATE i-particles in maximum to estimate bmax */
	  while( inum > NI_BMAX_ESTIMATE ){
	    //-------------------------------------------------------------
	    /* find rmax and the corresponding rmin */
	    //-------------------------------------------------------------
	    real rmax = REAL_MIN;
	    real rmin = ZERO;
	    //-------------------------------------------------------------
	    int Nloc = 0;
	    int Nbuf = 0;
	    //-------------------------------------------------------------
	    const int Niter0 = (Ntry0 + NJ_BMAX_ESTIMATE - 1) / NJ_BMAX_ESTIMATE;
	    for(int iter = 0; iter < Niter0; iter++){
	      //-----------------------------------------------------------
	      const int Nloop = (Ntry0 < NJ_BMAX_ESTIMATE) ? Ntry0 : NJ_BMAX_ESTIMATE;
	      //-----------------------------------------------------------
	      for(int jj = 0; jj < Nloop; jj++){
		//---------------------------------------------------------
		const uint more_tmp = more[list0[jj]];
		const int khead =      more_tmp  & IDXMASK;
		const int knum  = 1 + (more_tmp >> IDXBITS);
		//-----------------------------------------------------------
		for(int kk = 0; kk < knum; kk++){
		  //---------------------------------------------------------
		  const int kidx = khead + kk;
		  const jparticle jpos = pj[kidx];
		  //---------------------------------------------------------
		  const real dx = jpos.x - jcom.x;
		  const real dy = jpos.y - jcom.y;
		  const real dz = jpos.z - jcom.z;
		  const real disp = SQRT(EPSILON + dx * dx + dy * dy + dz * dz);
		  const real rjmax = disp + bmax[kidx];
		  //---------------------------------------------------------
		  const real rjmin = -rjmax + TWO * disp;
		  rmin = (rjmin > rmin) ? rjmin : rmin;
		  //---------------------------------------------------------

		  //---------------------------------------------------------
		  /* if the minimum required condition is satisfied */
		  if( rjmax > rmin ){
		    //-------------------------------------------------------
		    /* save on the local buffer */
		    rjmax_tmp[Nloc] = rjmax;
		    list1    [Nloc] = kidx;
		    Nloc++;
		    //-------------------------------------------------------
		    /* move data to the remote buffer if necessary */
		    if( Nloc >= NJ_BMAX_ESTIMATE ){
		      for(int ll = 0; ll < NJ_BMAX_ESTIMATE; ll++){
			rjmaxBuf[Nbuf + ll] = rjmax_tmp[ll];
			more1Buf[Nbuf + ll] = list1    [ll];
		      }
		      Nloc -= NJ_BMAX_ESTIMATE;
		      Nbuf += NJ_BMAX_ESTIMATE;
		    }
		    //-------------------------------------------------------
		    /* when the current node is the most distant node */
		    if( rjmax > rmax ){
		      //-----------------------------------------------------
		      rmax = rjmax;
		      /* const real rjmin = -rjmax + TWO * disp; */
		      /* rmin = (rjmin > ZERO) ? rjmin : ZERO; */
		      //-----------------------------------------------------
		    }
		    //-------------------------------------------------------
		  }
		  //---------------------------------------------------------
		}
		//-----------------------------------------------------------
	      }
	      //-------------------------------------------------------------

	      //-----------------------------------------------------------
	      Ntry0 -= Nloop;
	      //-----------------------------------------------------------
	      /* copy data from remote buffer to local buffer */
	      if( Ntry0 > 0 )
		for(int jj = 0; jj < NJ_BMAX_ESTIMATE; jj++)
		  list0[jj] = more0Buf[NJ_BMAX_ESTIMATE * (iter + 1) + jj];
	      //-----------------------------------------------------------
	    }
	    //-------------------------------------------------------------
	    if( Nbuf != 0 ){
	      //-----------------------------------------------------------
	      for(int ll = 0; ll < Nloc; ll++){
		rjmaxBuf[Nbuf + ll] = rjmax_tmp[ll];
		more1Buf[Nbuf + ll] = list1    [ll];
	      }
	      //-----------------------------------------------------------
	      for(int ll = 0; ll < NJ_BMAX_ESTIMATE; ll++){
		rjmax_tmp[ll] = rjmaxBuf[ll];
		list1    [ll] = more1Buf[ll];
	      }
	      //-----------------------------------------------------------
	    }
	    //-------------------------------------------------------------


	    //-------------------------------------------------------------
	    /* list up all child nodes that satisfy rjmax > rmin */
	    //-------------------------------------------------------------
	    int Ntry1 = Nbuf + Nloc;
	    if( Ntry1 > NUM_ALLOC_MACBUF )	      *overflow += 1;
	    Nloc = 0;
	    Nbuf = 0;
	    inum = 0;
	    //-------------------------------------------------------------
	    const int Niter1 = (Ntry1 + NJ_BMAX_ESTIMATE - 1) / NJ_BMAX_ESTIMATE;
	    for(int iter = 0; iter < Niter1; iter++){
	      //-----------------------------------------------------------
	      const int Nloop = (Ntry1 < NJ_BMAX_ESTIMATE) ? Ntry1 : NJ_BMAX_ESTIMATE;
	      //-----------------------------------------------------------
	      for(int jj = 0; jj < Nloop; jj++){
		//---------------------------------------------------------
		const int   kidx = list1    [jj];
		const real rjmax = rjmax_tmp[jj];
		//---------------------------------------------------------
		/* when the current node must be taken into account */
		if( rjmax > rmin ){
		  //-------------------------------------------------------
		  /* count up total number of contained i-particles */
		  inum += cell[node2cell[kidx]].num;
		  //-------------------------------------------------------
		  /* save on the local buffer */
		  list0[Nloc] = kidx;
		  Nloc++;
		  //-------------------------------------------------------
		  /* move data to the remote buffer if necessary */
		  if( Nloc >= NJ_BMAX_ESTIMATE ){
		    for(int ll = 0; ll < NJ_BMAX_ESTIMATE; ll++){
		      more0Buf[Nbuf + ll] = list0[ll];
		    }
		    Nloc -= NJ_BMAX_ESTIMATE;
		    Nbuf += NJ_BMAX_ESTIMATE;
		  }
		  //-------------------------------------------------------
		}
		//---------------------------------------------------------
	      }
	      //-----------------------------------------------------------

	      //-----------------------------------------------------------
	      Ntry1 -= Nloop;
	      //-----------------------------------------------------------
	      /* copy data from remote buffer to local buffer */
	      if( Ntry1 > 0 )
		for(int jj = 0; jj < NJ_BMAX_ESTIMATE; jj++){
		  rjmax_tmp[jj] = rjmaxBuf[NJ_BMAX_ESTIMATE * (iter + 1) + jj];
		  list1    [jj] = more1Buf[NJ_BMAX_ESTIMATE * (iter + 1) + jj];
		}
	      //-----------------------------------------------------------
	    }
	    //-------------------------------------------------------------
	    if( Nbuf != 0 ){
	      //-----------------------------------------------------------
	      for(int ll = 0; ll < Nloc            ; ll++)		more0Buf[Nbuf + ll] = list0   [ll];
	      for(int ll = 0; ll < NJ_BMAX_ESTIMATE; ll++)		list0   [       ll] = more0Buf[ll];
	      //-----------------------------------------------------------
	    }
	    //-------------------------------------------------------------
	    Ntry0 = Nloc + Nbuf;
	    if( Ntry0 > NUM_ALLOC_MACBUF )	      *overflow += 1;
	    //-------------------------------------------------------------
	  }
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
	  /* check positions of all the pick upped i-particles */
	  real jbmax = ZERO;
	  //---------------------------------------------------------------
	  const int Niter = (Ntry0 + NJ_BMAX_ESTIMATE - 1) / NJ_BMAX_ESTIMATE;
	  for(int iter = 0; iter < Niter; iter++){
	    //-------------------------------------------------------------
	    const int Nloop = (Ntry0 < NJ_BMAX_ESTIMATE) ? Ntry0 : NJ_BMAX_ESTIMATE;
	    //-------------------------------------------------------------
	    for(int jj = 0; jj < Nloop; jj++){
	      //-----------------------------------------------------------
	      const uint more_tmp = more[list0[jj]];
	      const int khead =      more_tmp  & IDXMASK;
	      const int knum  = 1 + (more_tmp >> IDXBITS);
	      //-----------------------------------------------------------
	      for(int kk = 0; kk < knum; kk++){
		//---------------------------------------------------------
		const treecell candidate = cell[node2cell[khead + kk]];
		const int phead = candidate.head;
		const int ptail = candidate.num + phead;
		//---------------------------------------------------------
		for(int pp = phead; pp < ptail; pp++){
		  //-------------------------------------------------------
		  const position ipos = pi[pp];
		  //-------------------------------------------------------
		  const real dx = ipos.x - jcom.x;
		  const real dy = ipos.y - jcom.y;
		  const real dz = ipos.z - jcom.z;
		  const real r2 = dx * dx + dy * dy + dz * dz;
		  //-------------------------------------------------------
		  if( r2 > jbmax )		    jbmax = r2;
		  //-------------------------------------------------------
		}
		//---------------------------------------------------------
	      }
	      //-----------------------------------------------------------
	    }
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    Ntry0 -= Nloop;
	    //-------------------------------------------------------------
	    /* copy data from remote buffer to local buffer */
	    if( Ntry0 > 0 )
	      for(int jj = 0; jj < NJ_BMAX_ESTIMATE; jj++)
		list0[jj] = more0Buf[NJ_BMAX_ESTIMATE * (iter + 1) + jj];
	    //-------------------------------------------------------------
	  }
	  //---------------------------------------------------------------
	  jbmax = SQRT(jbmax);
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
#ifdef  WS93_MAC
	  const real bmax_2 = HALF * jbmax;
	  jcom.w = bmax_2 + SQRT(bmax_2 * bmax_2 + SQRT(THREE * invDelta * B2));
#else///WS93_MAC
	  jcom.w = jbmax * jbmax;
#endif//WS93_MAC
	  //---------------------------------------------------------------


	  //---------------------------------------------------------------
	  /* commit a pseudo particle */
	  //---------------------------------------------------------------
	  pj  [jidx] = jcom;
	  mj  [jidx] = mtot;
	  mr2 [jidx] = mjrj2;
	  bmax[jidx] = jbmax;
	  //---------------------------------------------------------------
#ifdef  COUNT_INTERACTIONS
	  stats[levelIdx].nodeNum++;
	  stats[levelIdx].mjMean +=        mtot;	  stats[levelIdx].r2Mean +=          jcom.w;
	  stats[levelIdx].mjSdev += mtot * mtot;	  stats[levelIdx].r2Sdev += jcom.w * jcom.w;
#endif//COUNT_INTERACTIONS
	  //---------------------------------------------------------------
	}
	//-----------------------------------------------------------------
#ifndef SEPARATED_IMPLEMENTATION
	else{
	  //---------------------------------------------------------------
	  /* when the tree cell is a leaf cell, then initilize corresponding j-particle(s) */
	  //---------------------------------------------------------------
	  const int jh = root.head;
	  const int nj = root.num;
	  const int hidx = node[cidx] & IDXMASK;
	  for(int jj = 0; jj < nj; jj++){
	    //-------------------------------------------------------------
	    /* load an N-body particle */
	    const position ipos = pi[jh + jj];
	    //-------------------------------------------------------------
	    /* set the N-body particle as a real particle */
	    jparticle jpos;
	    jpos.x = ipos.x;
	    jpos.y = ipos.y;
	    jpos.z = ipos.z;
	    jpos.w = -UNITY;/* squared size for the real particle is set to be negative */
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    /* commit the new particle to the j-particle array */
	    const int jidx = hidx + jj;
	    pj [jidx] = jpos;
	    mj [jidx] = ipos.m;
	    mr2[jidx] = ipos.m * (jpos.x * jpos.x + jpos.y * jpos.y + jpos.z * jpos.z);
	    //-------------------------------------------------------------
	  }
	  //---------------------------------------------------------------
	}
#endif//SEPARATED_IMPLEMENTATION
	//-----------------------------------------------------------------
      }
      //-------------------------------------------------------------------
    }
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->calcMultipole));
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif
//-------------------------------------------------------------------------
#endif//CALC_MULTIPOLE_ON_DEVICE
//-------------------------------------------------------------------------
