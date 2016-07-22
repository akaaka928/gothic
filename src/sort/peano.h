/*************************************************************************\
 *                                                                       *
                  last updated on 2015/11/10(Tue) 15:46:41
 *                                                                       *
 *    Header File for constructing octree structure                      *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef PEANO_H
#define PEANO_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef _SYS_TIME_H
#      include <sys/time.h>
#endif//_SYS_TIME_H
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


//-------------------------------------------------------------------------
/* #define MAXIMUM_PHKEY_LEVEL (5) */
/* #define MAXIMUM_PHKEY_LEVEL (6) */
/* #define MAXIMUM_PHKEY_LEVEL (7) */
/* #define MAXIMUM_PHKEY_LEVEL (8) */
/* #define MAXIMUM_PHKEY_LEVEL (9) */
/* #define MAXIMUM_PHKEY_LEVEL (10) */
#define MAXIMUM_PHKEY_LEVEL (21)
//-------------------------------------------------------------------------
#define NUM_PHKEY_LEVEL (1 + MAXIMUM_PHKEY_LEVEL)
/* 1 means the level for the root cell */
//-------------------------------------------------------------------------
#if     MAXIMUM_PHKEY_LEVEL <= 10
typedef uint  PHint;
#else
typedef ulong PHint;
#endif//MAXIMUM_PHKEY_LEVEL <= 10
//-------------------------------------------------------------------------
#ifndef GENERATE_PHKEY_ON_DEVICE
typedef struct
{
  PHint key;
    int tag;
} identity;
#endif//GENERATE_PHKEY_ON_DEVICE
//-------------------------------------------------------------------------
typedef struct __align__(16)
{
  int level, num, head, nmax;
  /* nmax is the maximum number of particles contained a tree cell in this level */
} PHinfo;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "peano.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
#ifndef GENERATE_PHKEY_ON_DEVICE
  //-----------------------------------------------------------------------
  muse allocPeanoHilbertKey
  (const int num, PHint **peano, identity **tag
#ifndef CALC_MULTIPOLE_ON_DEVICE
   , PHinfo **info
#endif//CALC_MULTIPOLE_ON_DEVICE
   );
  void  freePeanoHilbertKey
  (PHint  *peano, identity  *tag
#ifndef CALC_MULTIPOLE_ON_DEVICE
   , PHinfo  *info
#endif//CALC_MULTIPOLE_ON_DEVICE
   );
  //-----------------------------------------------------------------------
  void sortParticlesPHcurve(int num, iparticle * src, iparticle * dst, identity * tag, PHint * peano
			    , struct timeval *start
#ifdef  EXEC_BENCHMARK
			    , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
			    );
  //-----------------------------------------------------------------------
  void initPHinfo(PHinfo *info);
  //-----------------------------------------------------------------------
#endif//GENERATE_PHKEY_ON_DEVICE
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//PEANO_H
//-------------------------------------------------------------------------
