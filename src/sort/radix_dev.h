/*************************************************************************\
 *                                                                       *
                  last updated on 2015/11/22(Sun) 15:11:15
 *                                                                       *
 *    Header File for radix sort on GPU                                  *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef RADIX_DEV_H
#define RADIX_DEV_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef CUDALIB_H
#       include <cudalib.h>
#endif//CUDALIB_H
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
/* options */
//-------------------------------------------------------------------------
#define USE_WARP_SHUFFLE_FUNC_SORT
//-------------------------------------------------------------------------
/* this option slows down the sorting */
/* #define USE_MASK_INSTEAD_OF_IF */
//-------------------------------------------------------------------------
#define INCREASE_INSTRUCTION_LEVEL_PARALLELISM
//-------------------------------------------------------------------------
#define ATOMIC_BASED_JOB_ASSIGNMENT
//-------------------------------------------------------------------------
/* #define ENABLE_EARLY_EXIT_BLCK */
/* #define ENABLE_EARLY_EXIT_GRID */
//-------------------------------------------------------------------------
/* #define CHECK_2BITS */
//-------------------------------------------------------------------------
#define SMEM_PREF_FOR_GLOBAL_SORT
//-------------------------------------------------------------------------
/* #define FLIP_LOOP_ORDER */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  defined(USE_WARP_SHUFFLE_FUNC_SORT) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_SORT
#endif//defined(USE_WARP_SHUFFLE_FUNC_SORT) && (GPUGEN < 30)
//-------------------------------------------------------------------------
#ifdef  CHECK_2BITS
//-------------------------------------------------------------------------
#define RADIX_SORT_CHECK_BITS (2)
#define RADIX_SORT_CHECK_MASK (3)
#define RADIX_SORT_NUM_BUCKET (4)
#define RADIX_SORT_NUM_OF_I32 (1)
/* RADIX_SORT_CHECK_MASK = (1 << RADIX_SORT_CHECK_BITS) - 1 */
/* RADIX_SORT_NUM_BUCKET =  1 << RADIX_SORT_CHECK_BITS */
/* assumes warpSize >= RADIX_SORT_NUM_BUCKET; warpSize of 32 always satisfies the condition */
/* RADIX_SORT_NUM_OF_I32 = RADIX_SORT_NUM_BUCKET / 4; a 32bits integer can contain 4 buckets within a warp */
//-------------------------------------------------------------------------
#else///CHECK_2BITS
//-------------------------------------------------------------------------
#define RADIX_SORT_CHECK_BITS (4)
#define RADIX_SORT_CHECK_MASK (15)
#define RADIX_SORT_NUM_BUCKET (16)
#define RADIX_SORT_NUM_OF_I32 (4)
/* RADIX_SORT_CHECK_MASK = (1 << RADIX_SORT_CHECK_BITS) - 1 */
/* RADIX_SORT_NUM_BUCKET =  1 << RADIX_SORT_CHECK_BITS */
/* assumes warpSize >= RADIX_SORT_NUM_BUCKET; warpSize of 32 always satisfies the condition */
/* RADIX_SORT_NUM_OF_I32 = RADIX_SORT_NUM_BUCKET / 4; a 32bits integer can contain 4 buckets within a warp */
//-------------------------------------------------------------------------
#endif//CHECK_2BITS
//-------------------------------------------------------------------------
#define RADIX_SORT_BOX_NUM_ALLOCATED (1024)
/* #define RADIX_SORT_BOX_NUM_ALLOCATED (32) */
/* #define RADIX_SORT_BOX_NUM_ALLOCATED (64) */
//-------------------------------------------------------------------------
#define RADIX_SORT_REM_LOC_MAX_GUESS (4)
//-------------------------------------------------------------------------
#ifndef NTHREADS_SORT
/* #define NTHREADS_SORT (1024) */
/* #define NTHREADS_SORT (512) */
/* #define NTHREADS_SORT (256) */
#define NTHREADS_SORT (128)
/* #define NTHREADS_SORT (64) */
#endif//NTHREADS_SORT
/* NTHREADS_SORT must be equal or smaller than 1024 due to the capacity of shared memory */
#   if  NTHREADS_SORT > 1024
#undef  NTHREADS_SORT
#define NTHREADS_SORT  (1024)
#endif//NTHREADS_SORT > 1024
//-------------------------------------------------------------------------
#ifndef NTHREADS_SORT_ACCUMULATION
#define NTHREADS_SORT_ACCUMULATION (1024)
#endif//NTHREADS_SORT_ACCUMULATION
//-------------------------------------------------------------------------
#ifndef TSUB_SORT
#define TSUB_SORT (32)
#endif//TSUB_SORT
/* TSUB_SORT must be equal or smaller than NTHREADS_SORT */
#   if  TSUB_SORT > NTHREADS_SORT
#undef  TSUB_SORT
#define TSUB_SORT   NTHREADS_SORT
#endif//TSUB_SORT > NTHREADS_SORT
//-------------------------------------------------------------------------
#   if  NTHREADS_SORT >= 128
/* TSUB_SORT must be equal or greater than 4 */
#          if  TSUB_SORT < 4
#       undef  TSUB_SORT
#       define TSUB_SORT  (4)
#       endif//TSUB_SORT < 4
#   if  NTHREADS_SORT >= 256
/* TSUB_SORT must be equal or greater than 8 */
#          if  TSUB_SORT < 8
#       undef  TSUB_SORT
#       define TSUB_SORT  (8)
#       endif//TSUB_SORT < 8
#   if  NTHREADS_SORT == 512
/* TSUB_SORT must be equal to 32 */
#          if  TSUB_SORT != 32
#       undef  TSUB_SORT
#       define TSUB_SORT   (32)
#       endif//TSUB_SORT != 32
#endif//NTHREADS_SORT == 512
#endif//NTHREADS_SORT >= 256
#endif//NTHREADS_SORT >= 128
//-------------------------------------------------------------------------
#define NGROUPS_SORT (NTHREADS_SORT / TSUB_SORT)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef union
{
  uint u;
   int i;
} uint_int;
typedef union
{
   uint u;
  float f;
} uint_float;
typedef union
{
  ulong u;
  long i;
} ulong_long;
typedef union
{
  ulong u;
  double f;
} ulong_double;
//-------------------------------------------------------------------------
typedef union
{
  uint4 m;
  uint  a[4];
} uint4_array;
typedef union
{
  uint2 m;
  uint  a[2];
} uint2_array;
typedef union
{
  uint2 u;
   int2 i;
} uint2_int2;
typedef union
{
  uint4 u;
   int4 i;
} uint4_int4;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//RADIX_DEV_H
//-------------------------------------------------------------------------
