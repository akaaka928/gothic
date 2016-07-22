/*************************************************************************\
 *                                                                       *
                  last updated on 2016/02/12(Fri) 10:37:44
 *                                                                       *
 *    Header File to execute Abel transform to deproject density profile *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef ABEL_H
#define ABEL_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MAGI_H
#include "../init/magi.h"
#endif//MAGI_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define NABEL (4096)
#define NDIVIDE_GAUSSQD4ABEL (32)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
  double *xx, *yy, *y2;/* pointers for column density profile in table form */
  double invRd;
  double ninv, bb;/* parameters for Sersic profile */
  double alpha, beta, gam, del, eps;/* parameters for two-power model (del := -1 - alpha; eps := (alpha - beta - gamma) / beta) */
  int num;/* parameter for column density profile in table form */
} profile_abel_cfg;
//-------------------------------------------------------------------------
typedef struct
{
  double (*getColumnDensityDerivative)(double, profile_abel_cfg);
  profile_abel_cfg cfg;
} abel_util;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "able.c"
//-------------------------------------------------------------------------
void execAbelTransform(profile *prf, const profile_cfg cfg, const double rmin, const double rmax, const profile_abel_cfg tmp);
//-------------------------------------------------------------------------
void readColumnDensityProfileTable(profile *prf, const double rs, char *file, const profile_cfg cfg);
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//ABEL_H
//-------------------------------------------------------------------------
