/*************************************************************************\
 *                                                                       *
                  last updated on 2015/02/11(Wed) 16:15:36
 *                                                                       *
 *    Header File to examine tree statistics                             *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef STAT_H
#define STAT_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* /\* #define DEFAULT_TOLERANCE_VALUE (2.0f) *\/ */
/* /\* #define DEFAULT_TOLERANCE_VALUE (1.0f) *\/ */
/* #define DEFAULT_TOLERANCE_VALUE (0.0f) */
//-------------------------------------------------------------------------
/* #define DEEPER_SPLIT_MODE */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "stat.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
#ifndef  __CUDACC__
  //-----------------------------------------------------------------------
  void setSubRootNode(uint *more_hst, int *subrootNum, int subrootIdx[restrict]);
  //-----------------------------------------------------------------------
#endif//__CUDACC__
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//STAT_H
//-------------------------------------------------------------------------
