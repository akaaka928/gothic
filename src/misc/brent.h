/**
 * @file brent.h
 *
 * @brief Header file for Brent's algorithm in GOTHIC
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/02/24 (Fri)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef BRENT_H
#define BRENT_H


#include <stdbool.h>


/* parameters for purturbing results of Brent's method (not to lock in a local minimum) */
#define BRENT_METHOD_ALLOW (1.1)
#define BRENT_METHOD_LAUNCH (5)
#define BRENT_METHOD_MODIFY (2)


/**
 * @struct brentFunc
 *
 * @brief structure for Brent's method
 */
typedef struct
{
  double pos;
  double val;
} brentFunc;
/**
 * @struct brentStatus
 *
 * @brief structure for Brent's method
 */
typedef struct
{
  /** x: position provides a local minimum up to the current step */
  /** w: position provides the next local minimum */
  /** v: ww in the previous step */
  /** u: latest guess of the local minimum */
  brentFunc x, w, v, u;
  double a, b;/**< a and b: edges of the region; a < b */
  double d, e;/**< displacement in the previous step and the step prior to the previous one, respectively */
  double gold;/**< the golden ratio */
  bool initialized;
} brentStatus;


/**
 * @struct brentMemory
 *
 * @brief structure for Brent's method
 */
typedef struct
{
  double previous;/**< elapsed time in the previous sequence (a sequence includes multiple steps of tree traversal and a tree build) */
  int totNum;/**< the sum of grpNum integrated in this sequence */
  int degraded;/**< # of sequences increasing the elapsed time continuously */
  int interval;/**< # of sequences from the previous initialization of the brentStatus */
} brentMemory;


/* list of functions appeared in ``brent.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  void brentInit1st(brentStatus *brent, double xl, double xr);
  void brentInit2nd(brentStatus *brent);
  void brentPerturb(brentStatus *brent, double xl, double xr);

  void brentCalc1st(brentStatus *brent, const double tol);
  void brentCalc2nd(brentStatus *brent);
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//BRENT_H
