/**
 * @file brent.c
 *
 * @brief Source code for Brent's algorithm in GOTHIC
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/02/21 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "macro.h"

#include "../misc/brent.h"


/**
 * @fn sign
 *
 * @brief Returns sign(b) * abs(a).
 *
 * @param (aa) variable a
 * @param (bb) variable b
 * @return sign(bb) * fabs(aa)
 */
static inline double sign(const double aa, const double bb){  return ((bb >= 0.0) ? (fabs(aa)) : (-fabs(aa)));}


/**
 * @fn brentInit1st
 *
 * @brief Initialize the Brent's method (first half).
 *
 * @return (brent) structure required for the Brent's meghod
 * @param (xl) lower limit of x
 * @param (xr) upper limit of x
 */
void brentInit1st(brentStatus *brent, double xl, double xr)
{
  __NOTE__("%s\n", "start");

  /* set the golden ratio */
  brent->gold = 0.5 * (3.0 - sqrt(5.0));

  brent->d = 0.0;
  brent->e = 0.0;

  brent->a = xl;
  brent->b = xr;

  brent->x.pos = 0.5 * (xl + xr);
  brent->x.val = 0.0;

  __NOTE__("%s\n", "end");
}
/**
 * @fn brentInit2nd
 *
 * @brief Initialize the Brent's method (second half).
 *
 * @return (brent) structure required for the Brent's meghod
 */
void brentInit2nd(brentStatus *brent)
{
  __NOTE__("%s\n", "start");

  brent->w = brent->v = brent->x;

  __NOTE__("%s\n", "end");
}


/**
 * @fn brentPerturb
 *
 * @brief Perturb the Brent's method.
 *
 * @return (brent) structure required for the Brent's meghod
 * @param (xl) lower limit of x
 * @param (xr) upper limit of x
 *
 * @sa brentInit1st
 * @sa brentInit2nd
 */
void brentPerturb(brentStatus *brent, double xl, double xr)
{
  __NOTE__("%s\n", "start");

  const brentFunc tmp = brent->x;
  brentInit1st(brent, xl, xr);

  brent->x = tmp;
  if( tmp.pos < brent->a )    brent->a = 0.75 * tmp.pos;
  if( tmp.pos > brent->b )    brent->b = 1.25 * tmp.pos;
  brentInit2nd(brent);

  __NOTE__("%s\n", "end");
}


/**
 * @fn brentCalc1st
 *
 * @brief Execute the Brent's method (first half).
 *
 * @return (brent) structure required for the Brent's meghod
 * @param (tol) tolerance parameter
 */
void brentCalc1st(brentStatus *brent, const double tol)
{
  __NOTE__("%s\n", "start");

  const double xm = 0.5 * (brent->a + brent->b);
  const double tol0 = tol * fabs(brent->x.pos) + DBL_EPSILON;
  const double tol1 = 2.0 * tol0;


  /* execute interpolation */
  if( fabs(brent->e) > tol0 ){
    /* try parabolic interpolation */
    double rr = (brent->x.pos - brent->w.pos) * (brent->x.val - brent->v.val);
    double qq = (brent->x.pos - brent->v.pos) * (brent->x.val - brent->w.val);
    double pp = (brent->x.pos - brent->v.pos) * qq - (brent->x.pos - brent->w.pos) * rr;
    qq = 2.0 * (qq - rr);
    if( qq > 0.0 )      pp = -pp;
    else	        qq = -qq;
    const double tmp = brent->e;
    brent->e = brent->d;

    /* validate the result of the parabolic interpolation */
    if( (fabs(pp) >= fabs(0.5 * qq * tmp)) || (pp <= qq * (brent->a - brent->x.pos)) || (pp >= qq * (brent->b - brent->x.pos)) ){
      /* reject the parabolic interpolation and adopt the golden section search */
      brent->e = ((brent->x.pos >= xm) ? (brent->a - brent->x.pos) : (brent->b - brent->x.pos));
      brent->d = brent->gold * brent->e;
    }/* if( (fabs(pp) >= fabs(0.5 * qq * tmp)) || (pp <= qq * (brent->a - brent->x.pos)) || (pp >= qq * (brent->b - brent->x.pos)) ){ */
    else{
      /* accept the parabolic interpolation */
      brent->d = pp / qq;
      double tt = brent->x.pos + brent->d;
      if( ((tt - brent->a) < tol1) || ((brent->b - tt) < tol1) )
	brent->d = sign(tol0, xm - brent->x.pos);
    }/* else{ */
  }/* if( fabs(brent->e) > tol0 ){ */
  else{
    /* adopt the golden section search */
    brent->e = ((brent->x.pos >= xm) ? (brent->a - brent->x.pos) : (brent->b - brent->x.pos));
    brent->d = brent->gold * brent->e;
  }/* else{ */

  brent->u.pos = ((fabs(brent->d) >= tol0) ? (brent->x.pos + brent->d) : (brent->x.pos + sign(tol0, brent->d)));
  brent->u.val = 0.0;


  __NOTE__("%s\n", "end");
}


/**
 * @fn brentCalc2nd
 *
 * @brief Execute the Brent's method (second half).
 *
 * @return (brent) structure required for the Brent's meghod
 */
void brentCalc2nd(brentStatus *brent)
{
  __NOTE__("%s\n", "start");

  if( brent->u.val <= brent->x.val ){
    if( brent->u.pos >= brent->x.pos )      brent->a = brent->x.pos;
    else	                            brent->b = brent->x.pos;

    brent->v = brent->w;
    brent->w = brent->x;
    brent->x = brent->u;
  }/* if( brent->u.val <= brent->x.val ){ */
  else{
    if( brent->u.pos < brent->x.pos )      brent->a = brent->u.pos;
    else	                           brent->b = brent->u.pos;

    if( (brent->u.val <= brent->w.val) || (brent->w.pos == brent->x.pos) ){
      brent->v = brent->w;
      brent->w = brent->u;
    }/* if( (brent->u.val <= brent->w.val) || (brent->w.pos == brent->x.pos) ){ */
    else
      if( (brent->u.val <= brent->v.val) || (brent->v.pos == brent->x.pos) || (brent->v.pos == brent->w.pos) )
	brent->v = brent->u;
  }/* else{ */

  __NOTE__("%s\n", "end");
}
