/*************************************************************************\
 *                                                                       *
                  last updated on 2016/08/11(Thu) 16:50:59
 *                                                                       *
 *    Header File to generate table f' and f'' from f by cubic spline    *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef TABLE_H
#define TABLE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef PROFILE_H
#include "../init/profile.h"
#endif//PROFILE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define NFIT (4)
#define NPUT (256)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "table.c"
//-------------------------------------------------------------------------
void getInterpolatedDensityProfile(const int num, profile * restrict prf, double * restrict xx, double * restrict ff);
void setDensityProfileTable(profile *prf, const double rs, char *file);
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//TABLE_H
//-------------------------------------------------------------------------
