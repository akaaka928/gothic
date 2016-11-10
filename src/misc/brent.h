/*************************************************************************\
 *                                                                       *
                  last updated on 2016/10/28(Fri) 17:15:32
 *                                                                       *
 *    Header File for implementing Brent's algorithm                     *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef BRENT_H
#define BRENT_H
//-------------------------------------------------------------------------
#ifndef USE_BRENT_METHOD
#define TEST_BRENT_METHOD
#endif//USE_BRENT_METHOD
//-------------------------------------------------------------------------
#include <stdbool.h>
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define BRENT_METHOD_ALLOW (1.1)
#define BRENT_METHOD_LAUNCH (5)
#define BRENT_METHOD_MODIFY (2)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
  double pos;
  double val;
} brentFunc;
//-------------------------------------------------------------------------
typedef struct
{
  /* x: position provides a local minimum up to the current step */
  /* w: position provides the next local minimum */
  /* v: ww in the previous step */
  /* u: latest guess of the local minimum */
  brentFunc x, w, v, u;
  double a, b;/* a and b: edges of the region; a < b */
  double d, e;/* displacement in the previous step and the step prior to the previous one, respectively */
  double gold;/* the golden ratio */
  bool initialized;
} brentStatus;
//-------------------------------------------------------------------------
#ifdef  USE_BRENT_METHOD
typedef struct
{
  double previous;/* elapsed time in the previous sequence (a sequence includes multiple steps of tree traversal and a tree build) */
  int totNum;/* the sum of grpNum integrated in this sequence */
  int degraded;/* # of sequences increasing the elapsed time continuously */
  int interval;/* # of sequences from the previous initialization of the brentStatus */
} brentMemory;
#endif//USE_BRENT_METHOD
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "brent.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  void brentInit1st(brentStatus *brent, double xl, double xr);
  void brentInit2nd(brentStatus *brent);
  void brentPerturb(brentStatus *brent, double xl, double xr);
  //-----------------------------------------------------------------------
#ifdef  TEST_BRENT_METHOD
  int  brentCalc1st(brentStatus *brent, const double tol);
#else///TEST_BRENT_METHOD
  void brentCalc1st(brentStatus *brent, const double tol);
#endif//TEST_BRENT_METHOD
  //-----------------------------------------------------------------------
  void brentCalc2nd(brentStatus *brent);
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//BRENT_H
//-------------------------------------------------------------------------
