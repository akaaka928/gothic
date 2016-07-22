/*************************************************************************\
 *                                                                       *
                  last updated on 2015/08/24(Mon) 14:57:33
 *                                                                       *
 *    Building block of radix sort library on GPU                        *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <helper_cuda.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#include "../misc/gsync_dev.cu"
//-------------------------------------------------------------------------
#include "../sort/radix_dev.h"
#include "../sort/radix_inc.h"
//-------------------------------------------------------------------------
#ifndef RADIX_INC_CU_MULTI_CALL
#define RADIX_INC_CU_FIRST_CALL
#endif//RADIX_INC_CU_MULTI_CALL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  RADIX_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------
/* flip the sign bit for positive values, all bits for negative values */
__device__ __forceinline__ uint  flip32flt(const uint  src){  uint  mask = -int (src >> 31) | 0x80000000;          return (src ^ mask);}
__device__ __forceinline__ ulong flip64flt(const ulong src){  ulong mask = -long(src >> 63) | 0x8000000000000000;  return (src ^ mask);}
//-------------------------------------------------------------------------
/* recover the original value from the flipped key */
__device__ __forceinline__ uint  undo32flt(const uint  src){  uint  mask = ((src >> 31) - 1) | 0x80000000;          return (src ^ mask);}
__device__ __forceinline__ ulong undo64flt(const ulong src){  ulong mask = ((src >> 63) - 1) | 0x8000000000000000;  return (src ^ mask);}
//-------------------------------------------------------------------------
/* flip only the sign bit */
__device__ __forceinline__ uint  flip32int(const uint  src){  return (src ^ 0x80000000        );}
__device__ __forceinline__ ulong flip64int(const ulong src){  return (src ^ 0x8000000000000000);}
//-------------------------------------------------------------------------
#   if  RADIX_SORT_CHECK_BITS <= 2
__device__ __forceinline__ int RADIX_SORT_TSUB_GET_DSTIDX
(const int subkey, const uint psum)
{
  //-----------------------------------------------------------------------
  /* subkey & 3 == 0 ==> use    lowest 8 bits (0x000000ff) */
  /* subkey & 3 == 1 ==> use      next 8 bits (0x0000ff00) */
  /* subkey & 3 == 2 ==> use      next 8 bits (0x00ff0000) */
  /* subkey & 3 == 3 ==> use uppermost 8 bits (0xff000000) */
  return ((psum >> ((subkey & 3) * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK);
  //-----------------------------------------------------------------------
}
#else///RADIX_SORT_CHECK_BITS <= 2
/* __device__ __forceinline__ int RADIX_SORT_TSUB_GET_DSTIDX */
/* (const int subkey, const uint4_array psum) */
__device__ __forceinline__ int RADIX_SORT_TSUB_GET_DSTIDX
(const int subkey, const int bits)
{
  //-----------------------------------------------------------------------
  /* use psum?.a[subkey / 4] */
  /* const int bits = psum.a[subkey >> 2]; */
  //-----------------------------------------------------------------------
  /* subkey & 3 == 0 ==> use    lowest 8 bits (0x000000ff) */
  /* subkey & 3 == 1 ==> use      next 8 bits (0x0000ff00) */
  /* subkey & 3 == 2 ==> use      next 8 bits (0x00ff0000) */
  /* subkey & 3 == 3 ==> use uppermost 8 bits (0xff000000) */
  return ((bits >> ((subkey & 3) * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK);
  //-----------------------------------------------------------------------
}
#endif//RADIX_SORT_CHECK_BITS <= 2
//-------------------------------------------------------------------------
#   if  RADIX_SORT_CHECK_BITS <= 2
__device__ __forceinline__ int RADIX_SORT_BLCK_GET_DSTIDX
(const int subkey, const uint psum0, const uint psum1)
{
  //-----------------------------------------------------------------------
  /* subkey & 2 != 2 ==> use psum0 */
  /* subkey & 2 == 2 ==> use psum1 */
  uint psum;  psum = ((subkey & 2) >> 1) ? psum1 : psum0;
  //-----------------------------------------------------------------------
  /* subkey & 1 == 0 ==> use lower 16 bits */
  /* subkey & 1 == 1 ==> use upper 16 bits */
  return ((psum >> ((subkey & 1) * RADIX_BLOCK_SHIFT_BITS)) & RADIX_BLOCK_SHIFT_MASK);
  //-----------------------------------------------------------------------
}
#else///RADIX_SORT_CHECK_BITS <= 2
/* __device__ __forceinline__ int RADIX_SORT_BLCK_GET_DSTIDX */
/* (const int subkey, const uint4 psum0, const uint4 psum1) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   /\* subkey & 2 != 2 ==> use psum0 *\/ */
/*   /\* subkey & 2 == 2 ==> use psum1 *\/ */
/*   uint4_array psum;  psum.m = ((subkey & 2) >> 1) ? psum1 : psum0; */
/*   //----------------------------------------------------------------------- */
/*   /\* use psum?.a[subkey / 4] *\/ */
/*   const int bits = psum.a[subkey >> 2]; */
/*   //----------------------------------------------------------------------- */
/*   /\* subkey & 1 == 0 ==> use lower 16 bits *\/ */
/*   /\* subkey & 1 == 1 ==> use upper 16 bits *\/ */
/*   return ((bits >> ((subkey & 1) * RADIX_BLOCK_SHIFT_BITS)) & RADIX_BLOCK_SHIFT_MASK); */
/*   //----------------------------------------------------------------------- */
/* } */
__device__ __forceinline__ int RADIX_SORT_BLCK_GET_DSTIDX
(const int subkey, const uint4 psum0, const uint4 psum1, uint4_array *psum)
{
  //-----------------------------------------------------------------------
  /* subkey & 2 != 2 ==> use psum0 */
  /* subkey & 2 == 2 ==> use psum1 */
  psum->m = ((subkey & 2) >> 1) ? psum1 : psum0;
  //-----------------------------------------------------------------------
  /* use psum?.a[subkey / 4] */
  const int bits = psum->a[subkey >> 2];
  //-----------------------------------------------------------------------
  /* subkey & 1 == 0 ==> use lower 16 bits */
  /* subkey & 1 == 1 ==> use upper 16 bits */
  return ((bits >> ((subkey & 1) * RADIX_BLOCK_SHIFT_BITS)) & RADIX_BLOCK_SHIFT_MASK);
  //-----------------------------------------------------------------------
}
#endif//RADIX_SORT_CHECK_BITS <= 2
//-------------------------------------------------------------------------
#endif//RADIX_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* parallel prefix sum within TSUB_SORT threads (use implicit synchronization) */
/* type of prefix sum is inclusive */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_PREFIX_SUM_TSUB
(const int lane
#   if  RADIX_SORT_CHECK_BITS <= 2
 , const RADIX_SORT_UINT psum, RADIX_SORT_UINT * RESTRICT ret
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
 , const uint4 psum , uint4 * RESTRICT ret
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
 , const uint4 psum0, uint4 * RESTRICT ret0, const uint4 psum1, uint4 * RESTRICT ret1
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
 , const uint4 psum2, uint4 * RESTRICT ret2, const uint4 psum3, uint4 * RESTRICT ret3
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
 , const int hidx, volatile uint * RESTRICT smem
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
 , const int m0
#   if  TSUB_SORT >=  4
 , const int m1
#   if  TSUB_SORT >=  8
 , const int m2
#   if  TSUB_SORT >= 16
 , const int m3
#   if  TSUB_SORT == 32
 , const int m4
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
#   if  RADIX_SORT_CHECK_BITS <= 2
  /* load index */
  RADIX_SORT_UINT_INT val;  val.u = psum;
  RADIX_SORT_UINT_INT tmp;
  /* calculate inclusive prefix sum */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
  tmp.i = __shfl_up(val.i,  1, TSUB_SORT) & m0;  val.u += tmp.u;
#   if  TSUB_SORT >=  4
  tmp.i = __shfl_up(val.i,  2, TSUB_SORT) & m1;  val.u += tmp.u;
#   if  TSUB_SORT >=  8
  tmp.i = __shfl_up(val.i,  4, TSUB_SORT) & m2;  val.u += tmp.u;
#   if  TSUB_SORT >= 16
  tmp.i = __shfl_up(val.i,  8, TSUB_SORT) & m3;  val.u += tmp.u;
#   if  TSUB_SORT == 32
  tmp.i = __shfl_up(val.i, 16, TSUB_SORT) & m4;  val.u += tmp.u;
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
  tmp.i = __shfl_up(val.i,  1, TSUB_SORT);  if( lane >=  1 )    val.u += tmp.u;
#   if  TSUB_SORT >=  4
  tmp.i = __shfl_up(val.i,  2, TSUB_SORT);  if( lane >=  2 )    val.u += tmp.u;
#   if  TSUB_SORT >=  8
  tmp.i = __shfl_up(val.i,  4, TSUB_SORT);  if( lane >=  4 )    val.u += tmp.u;
#   if  TSUB_SORT >= 16
  tmp.i = __shfl_up(val.i,  8, TSUB_SORT);  if( lane >=  8 )    val.u += tmp.u;
#   if  TSUB_SORT == 32
  tmp.i = __shfl_up(val.i, 16, TSUB_SORT);  if( lane >= 16 )    val.u += tmp.u;
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
  tmp.i.x = __shfl_up(val.i.x,  1, TSUB_SORT) & m0;  tmp.i.y = __shfl_up(val.i.y,  1, TSUB_SORT) & m0;  val.u.x += tmp.u.x;  val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  1, TSUB_SORT) & m0;  tmp.i.w = __shfl_up(val.i.w,  1, TSUB_SORT) & m0;  val.u.z += tmp.u.z;  val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
#   if  TSUB_SORT >=  4
  tmp.i.x = __shfl_up(val.i.x,  2, TSUB_SORT) & m1;  tmp.i.y = __shfl_up(val.i.y,  2, TSUB_SORT) & m1;  val.u.x += tmp.u.x;  val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  2, TSUB_SORT) & m1;  tmp.i.w = __shfl_up(val.i.w,  2, TSUB_SORT) & m1;  val.u.z += tmp.u.z;  val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
#   if  TSUB_SORT >=  8
  tmp.i.x = __shfl_up(val.i.x,  4, TSUB_SORT) & m2;  tmp.i.y = __shfl_up(val.i.y,  4, TSUB_SORT) & m2;  val.u.x += tmp.u.x;  val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  4, TSUB_SORT) & m2;  tmp.i.w = __shfl_up(val.i.w,  4, TSUB_SORT) & m2;  val.u.z += tmp.u.z;  val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
#   if  TSUB_SORT >= 16
  tmp.i.x = __shfl_up(val.i.x,  8, TSUB_SORT) & m3;  tmp.i.y = __shfl_up(val.i.y,  8, TSUB_SORT) & m3;  val.u.x += tmp.u.x;  val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  8, TSUB_SORT) & m3;  tmp.i.w = __shfl_up(val.i.w,  8, TSUB_SORT) & m3;  val.u.z += tmp.u.z;  val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
#   if  TSUB_SORT == 32
  tmp.i.x = __shfl_up(val.i.x, 16, TSUB_SORT) & m4;  tmp.i.y = __shfl_up(val.i.y, 16, TSUB_SORT) & m4;  val.u.x += tmp.u.x;  val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z, 16, TSUB_SORT) & m4;  tmp.i.w = __shfl_up(val.i.w, 16, TSUB_SORT) & m4;  val.u.z += tmp.u.z;  val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
  tmp.i.x = __shfl_up(val.i.x,  1, TSUB_SORT);  tmp.i.y = __shfl_up(val.i.y,  1, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  1, TSUB_SORT);  tmp.i.w = __shfl_up(val.i.w,  1, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  if( lane >=  1 ){
    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  }/* if( lane >=  1 ){ */
#   if  TSUB_SORT >=  4
  tmp.i.x = __shfl_up(val.i.x,  2, TSUB_SORT);  tmp.i.y = __shfl_up(val.i.y,  2, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  2, TSUB_SORT);  tmp.i.w = __shfl_up(val.i.w,  2, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  if( lane >=  2 ){
    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  }/* if( lane >=  2 ){ */
#   if  TSUB_SORT >=  8
  tmp.i.x = __shfl_up(val.i.x,  4, TSUB_SORT);  tmp.i.y = __shfl_up(val.i.y,  4, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  4, TSUB_SORT);  tmp.i.w = __shfl_up(val.i.w,  4, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  if( lane >=  4 ){
    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  }/* if( lane >=  4 ){ */
#   if  TSUB_SORT >= 16
  tmp.i.x = __shfl_up(val.i.x,  8, TSUB_SORT);  tmp.i.y = __shfl_up(val.i.y,  8, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  8, TSUB_SORT);  tmp.i.w = __shfl_up(val.i.w,  8, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  if( lane >=  8 ){
    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  }/* if( lane >=  8 ){ */
#   if  TSUB_SORT == 32
  tmp.i.x = __shfl_up(val.i.x, 16, TSUB_SORT);  tmp.i.y = __shfl_up(val.i.y, 16, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z, 16, TSUB_SORT);  tmp.i.w = __shfl_up(val.i.w, 16, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  if( lane >= 16 ){
    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  }/* if( lane >= 16 ){ */
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
  /* return calculated inclusive prefix sum */
  *ret = val.u;
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  /* load index */
  uint4_int4 val;  val.u = psum;
  uint4_int4 tmp;
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
  tmp.i.z = __shfl_up(val.i.x,  1, TSUB_SORT) & m0;  tmp.i.w = __shfl_up(val.i.y,  1, TSUB_SORT) & m0;
  tmp.i.x = __shfl_up(val.i.z,  1, TSUB_SORT) & m0;  tmp.i.y = __shfl_up(val.i.w,  1, TSUB_SORT) & m0;
  val.u.x += tmp.u.z;  val.u.y += tmp.u.w;  val.u.z += tmp.u.x;  val.u.w += tmp.u.y;
#   if  TSUB_SORT >=  4
  tmp.i.z = __shfl_up(val.i.x,  2, TSUB_SORT) & m1;  tmp.i.w = __shfl_up(val.i.y,  2, TSUB_SORT) & m1;
  tmp.i.x = __shfl_up(val.i.z,  2, TSUB_SORT) & m1;  tmp.i.y = __shfl_up(val.i.w,  2, TSUB_SORT) & m1;
  val.u.x += tmp.u.z;  val.u.y += tmp.u.w;  val.u.z += tmp.u.x;  val.u.w += tmp.u.y;
#   if  TSUB_SORT >=  8
  tmp.i.z = __shfl_up(val.i.x,  4, TSUB_SORT) & m2;  tmp.i.w = __shfl_up(val.i.y,  4, TSUB_SORT) & m2;
  tmp.i.x = __shfl_up(val.i.z,  4, TSUB_SORT) & m2;  tmp.i.y = __shfl_up(val.i.w,  4, TSUB_SORT) & m2;
  val.u.x += tmp.u.z;  val.u.y += tmp.u.w;  val.u.z += tmp.u.x;  val.u.w += tmp.u.y;
#   if  TSUB_SORT >= 16
  tmp.i.z = __shfl_up(val.i.x,  8, TSUB_SORT) & m3;  tmp.i.w = __shfl_up(val.i.y,  8, TSUB_SORT) & m3;
  tmp.i.x = __shfl_up(val.i.z,  8, TSUB_SORT) & m3;  tmp.i.y = __shfl_up(val.i.w,  8, TSUB_SORT) & m3;
  val.u.x += tmp.u.z;  val.u.y += tmp.u.w;  val.u.z += tmp.u.x;  val.u.w += tmp.u.y;
#   if  TSUB_SORT == 32
  tmp.i.z = __shfl_up(val.i.x, 16, TSUB_SORT) & m4;  tmp.i.w = __shfl_up(val.i.y, 16, TSUB_SORT) & m4;
  tmp.i.x = __shfl_up(val.i.z, 16, TSUB_SORT) & m4;  tmp.i.y = __shfl_up(val.i.w, 16, TSUB_SORT) & m4;
  val.u.x += tmp.u.z;  val.u.y += tmp.u.w;  val.u.z += tmp.u.x;  val.u.w += tmp.u.y;
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
  tmp.i.x = __shfl_up(val.i.x,  1, TSUB_SORT);  tmp.i.y = __shfl_up(val.i.y,  1, TSUB_SORT);
  tmp.i.z = __shfl_up(val.i.z,  1, TSUB_SORT);  tmp.i.w = __shfl_up(val.i.w,  1, TSUB_SORT);
  if( lane >=  1 ){    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;  }
#   if  TSUB_SORT >=  4
  tmp.i.x = __shfl_up(val.i.x,  2, TSUB_SORT);  tmp.i.y = __shfl_up(val.i.y,  2, TSUB_SORT);
  tmp.i.z = __shfl_up(val.i.z,  2, TSUB_SORT);  tmp.i.w = __shfl_up(val.i.w,  2, TSUB_SORT);
  if( lane >=  2 ){    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;  }
#   if  TSUB_SORT >=  8
  tmp.i.x = __shfl_up(val.i.x,  4, TSUB_SORT);  tmp.i.y = __shfl_up(val.i.y,  4, TSUB_SORT);
  tmp.i.z = __shfl_up(val.i.z,  4, TSUB_SORT);  tmp.i.w = __shfl_up(val.i.w,  4, TSUB_SORT);
  if( lane >=  4 ){    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;  }
#   if  TSUB_SORT >= 16
  tmp.i.x = __shfl_up(val.i.x,  8, TSUB_SORT);  tmp.i.y = __shfl_up(val.i.y,  8, TSUB_SORT);
  tmp.i.z = __shfl_up(val.i.z,  8, TSUB_SORT);  tmp.i.w = __shfl_up(val.i.w,  8, TSUB_SORT);
  if( lane >=  8 ){    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;  }
#   if  TSUB_SORT == 32
  tmp.i.x = __shfl_up(val.i.x, 16, TSUB_SORT);  tmp.i.y = __shfl_up(val.i.y, 16, TSUB_SORT);
  tmp.i.z = __shfl_up(val.i.z, 16, TSUB_SORT);  tmp.i.w = __shfl_up(val.i.w, 16, TSUB_SORT);
  if( lane >= 16 ){    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;  }
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  /* load index */
  uint4_int4 val0, val1;  val0.u = psum0;  val1.u = psum1;
  uint4_int4 tmp0, tmp1;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  uint4_int4 val2, val3;  val2.u = psum2;  val3.u = psum3;
  uint4_int4 tmp2, tmp3;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
  tmp0.i.z = __shfl_up(val0.i.x,  1, TSUB_SORT) & m0;  tmp1.i.z = __shfl_up(val1.i.x,  1, TSUB_SORT) & m0;
  tmp0.i.w = __shfl_up(val0.i.y,  1, TSUB_SORT) & m0;  tmp1.i.w = __shfl_up(val1.i.y,  1, TSUB_SORT) & m0;
  tmp0.i.x = __shfl_up(val0.i.z,  1, TSUB_SORT) & m0;  tmp1.i.x = __shfl_up(val1.i.z,  1, TSUB_SORT) & m0;
  tmp0.i.y = __shfl_up(val0.i.w,  1, TSUB_SORT) & m0;  tmp1.i.y = __shfl_up(val1.i.w,  1, TSUB_SORT) & m0;
  val0.u.x += tmp0.u.z;  val1.u.x += tmp1.u.z;  val0.u.y += tmp0.u.w;  val1.u.y += tmp1.u.w;
  val0.u.z += tmp0.u.x;  val1.u.z += tmp1.u.x;  val0.u.w += tmp0.u.y;  val1.u.w += tmp1.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.z = __shfl_up(val2.i.x,  1, TSUB_SORT) & m0;  tmp3.i.z = __shfl_up(val3.i.x,  1, TSUB_SORT) & m0;
  tmp2.i.w = __shfl_up(val2.i.y,  1, TSUB_SORT) & m0;  tmp3.i.w = __shfl_up(val3.i.y,  1, TSUB_SORT) & m0;
  tmp2.i.x = __shfl_up(val2.i.z,  1, TSUB_SORT) & m0;  tmp3.i.x = __shfl_up(val3.i.z,  1, TSUB_SORT) & m0;
  tmp2.i.y = __shfl_up(val2.i.w,  1, TSUB_SORT) & m0;  tmp3.i.y = __shfl_up(val3.i.w,  1, TSUB_SORT) & m0;
  val2.u.x += tmp2.u.z;  val3.u.x += tmp3.u.z;  val2.u.y += tmp2.u.w;  val3.u.y += tmp3.u.w;
  val2.u.z += tmp2.u.x;  val3.u.z += tmp3.u.x;  val2.u.w += tmp2.u.y;  val3.u.w += tmp3.u.y;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#   if  TSUB_SORT >=  4
  tmp0.i.z = __shfl_up(val0.i.x,  2, TSUB_SORT) & m1;  tmp1.i.z = __shfl_up(val1.i.x,  2, TSUB_SORT) & m1;
  tmp0.i.w = __shfl_up(val0.i.y,  2, TSUB_SORT) & m1;  tmp1.i.w = __shfl_up(val1.i.y,  2, TSUB_SORT) & m1;
  tmp0.i.x = __shfl_up(val0.i.z,  2, TSUB_SORT) & m1;  tmp1.i.x = __shfl_up(val1.i.z,  2, TSUB_SORT) & m1;
  tmp0.i.y = __shfl_up(val0.i.w,  2, TSUB_SORT) & m1;  tmp1.i.y = __shfl_up(val1.i.w,  2, TSUB_SORT) & m1;
  val0.u.x += tmp0.u.z;  val1.u.x += tmp1.u.z;  val0.u.y += tmp0.u.w;  val1.u.y += tmp1.u.w;
  val0.u.z += tmp0.u.x;  val1.u.z += tmp1.u.x;  val0.u.w += tmp0.u.y;  val1.u.w += tmp1.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.z = __shfl_up(val2.i.x,  2, TSUB_SORT) & m1;  tmp3.i.z = __shfl_up(val3.i.x,  2, TSUB_SORT) & m1;
  tmp2.i.w = __shfl_up(val2.i.y,  2, TSUB_SORT) & m1;  tmp3.i.w = __shfl_up(val3.i.y,  2, TSUB_SORT) & m1;
  tmp2.i.x = __shfl_up(val2.i.z,  2, TSUB_SORT) & m1;  tmp3.i.x = __shfl_up(val3.i.z,  2, TSUB_SORT) & m1;
  tmp2.i.y = __shfl_up(val2.i.w,  2, TSUB_SORT) & m1;  tmp3.i.y = __shfl_up(val3.i.w,  2, TSUB_SORT) & m1;
  val2.u.x += tmp2.u.z;  val3.u.x += tmp3.u.z;  val2.u.y += tmp2.u.w;  val3.u.y += tmp3.u.w;
  val2.u.z += tmp2.u.x;  val3.u.z += tmp3.u.x;  val2.u.w += tmp2.u.y;  val3.u.w += tmp3.u.y;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#   if  TSUB_SORT >=  8
  tmp0.i.z = __shfl_up(val0.i.x,  4, TSUB_SORT) & m2;  tmp1.i.z = __shfl_up(val1.i.x,  4, TSUB_SORT) & m2;
  tmp0.i.w = __shfl_up(val0.i.y,  4, TSUB_SORT) & m2;  tmp1.i.w = __shfl_up(val1.i.y,  4, TSUB_SORT) & m2;
  tmp0.i.x = __shfl_up(val0.i.z,  4, TSUB_SORT) & m2;  tmp1.i.x = __shfl_up(val1.i.z,  4, TSUB_SORT) & m2;
  tmp0.i.y = __shfl_up(val0.i.w,  4, TSUB_SORT) & m2;  tmp1.i.y = __shfl_up(val1.i.w,  4, TSUB_SORT) & m2;
  val0.u.x += tmp0.u.z;  val1.u.x += tmp1.u.z;  val0.u.y += tmp0.u.w;  val1.u.y += tmp1.u.w;
  val0.u.z += tmp0.u.x;  val1.u.z += tmp1.u.x;  val0.u.w += tmp0.u.y;  val1.u.w += tmp1.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.z = __shfl_up(val2.i.x,  4, TSUB_SORT) & m2;  tmp3.i.z = __shfl_up(val3.i.x,  4, TSUB_SORT) & m2;
  tmp2.i.w = __shfl_up(val2.i.y,  4, TSUB_SORT) & m2;  tmp3.i.w = __shfl_up(val3.i.y,  4, TSUB_SORT) & m2;
  tmp2.i.x = __shfl_up(val2.i.z,  4, TSUB_SORT) & m2;  tmp3.i.x = __shfl_up(val3.i.z,  4, TSUB_SORT) & m2;
  tmp2.i.y = __shfl_up(val2.i.w,  4, TSUB_SORT) & m2;  tmp3.i.y = __shfl_up(val3.i.w,  4, TSUB_SORT) & m2;
  val2.u.x += tmp2.u.z;  val3.u.x += tmp3.u.z;  val2.u.y += tmp2.u.w;  val3.u.y += tmp3.u.w;
  val2.u.z += tmp2.u.x;  val3.u.z += tmp3.u.x;  val2.u.w += tmp2.u.y;  val3.u.w += tmp3.u.y;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#   if  TSUB_SORT >= 16
  tmp0.i.z = __shfl_up(val0.i.x,  8, TSUB_SORT) & m3;  tmp1.i.z = __shfl_up(val1.i.x,  8, TSUB_SORT) & m3;
  tmp0.i.w = __shfl_up(val0.i.y,  8, TSUB_SORT) & m3;  tmp1.i.w = __shfl_up(val1.i.y,  8, TSUB_SORT) & m3;
  tmp0.i.x = __shfl_up(val0.i.z,  8, TSUB_SORT) & m3;  tmp1.i.x = __shfl_up(val1.i.z,  8, TSUB_SORT) & m3;
  tmp0.i.y = __shfl_up(val0.i.w,  8, TSUB_SORT) & m3;  tmp1.i.y = __shfl_up(val1.i.w,  8, TSUB_SORT) & m3;
  val0.u.x += tmp0.u.z;  val1.u.x += tmp1.u.z;  val0.u.y += tmp0.u.w;  val1.u.y += tmp1.u.w;
  val0.u.z += tmp0.u.x;  val1.u.z += tmp1.u.x;  val0.u.w += tmp0.u.y;  val1.u.w += tmp1.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.z = __shfl_up(val2.i.x,  8, TSUB_SORT) & m3;  tmp3.i.z = __shfl_up(val3.i.x,  8, TSUB_SORT) & m3;
  tmp2.i.w = __shfl_up(val2.i.y,  8, TSUB_SORT) & m3;  tmp3.i.w = __shfl_up(val3.i.y,  8, TSUB_SORT) & m3;
  tmp2.i.x = __shfl_up(val2.i.z,  8, TSUB_SORT) & m3;  tmp3.i.x = __shfl_up(val3.i.z,  8, TSUB_SORT) & m3;
  tmp2.i.y = __shfl_up(val2.i.w,  8, TSUB_SORT) & m3;  tmp3.i.y = __shfl_up(val3.i.w,  8, TSUB_SORT) & m3;
  val2.u.x += tmp2.u.z;  val3.u.x += tmp3.u.z;  val2.u.y += tmp2.u.w;  val3.u.y += tmp3.u.w;
  val2.u.z += tmp2.u.x;  val3.u.z += tmp3.u.x;  val2.u.w += tmp2.u.y;  val3.u.w += tmp3.u.y;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#   if  TSUB_SORT == 32
  tmp0.i.z = __shfl_up(val0.i.x, 16, TSUB_SORT) & m4;  tmp1.i.z = __shfl_up(val1.i.x, 16, TSUB_SORT) & m4;
  tmp0.i.w = __shfl_up(val0.i.y, 16, TSUB_SORT) & m4;  tmp1.i.w = __shfl_up(val1.i.y, 16, TSUB_SORT) & m4;
  tmp0.i.x = __shfl_up(val0.i.z, 16, TSUB_SORT) & m4;  tmp1.i.x = __shfl_up(val1.i.z, 16, TSUB_SORT) & m4;
  tmp0.i.y = __shfl_up(val0.i.w, 16, TSUB_SORT) & m4;  tmp1.i.y = __shfl_up(val1.i.w, 16, TSUB_SORT) & m4;
  val0.u.x += tmp0.u.z;  val1.u.x += tmp1.u.z;  val0.u.y += tmp0.u.w;  val1.u.y += tmp1.u.w;
  val0.u.z += tmp0.u.x;  val1.u.z += tmp1.u.x;  val0.u.w += tmp0.u.y;  val1.u.w += tmp1.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.z = __shfl_up(val2.i.x, 16, TSUB_SORT) & m4;  tmp3.i.z = __shfl_up(val3.i.x, 16, TSUB_SORT) & m4;
  tmp2.i.w = __shfl_up(val2.i.y, 16, TSUB_SORT) & m4;  tmp3.i.w = __shfl_up(val3.i.y, 16, TSUB_SORT) & m4;
  tmp2.i.x = __shfl_up(val2.i.z, 16, TSUB_SORT) & m4;  tmp3.i.x = __shfl_up(val3.i.z, 16, TSUB_SORT) & m4;
  tmp2.i.y = __shfl_up(val2.i.w, 16, TSUB_SORT) & m4;  tmp3.i.y = __shfl_up(val3.i.w, 16, TSUB_SORT) & m4;
  val2.u.x += tmp2.u.z;  val3.u.x += tmp3.u.z;  val2.u.y += tmp2.u.w;  val3.u.y += tmp3.u.w;
  val2.u.z += tmp2.u.x;  val3.u.z += tmp3.u.x;  val2.u.w += tmp2.u.y;  val3.u.w += tmp3.u.y;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
  tmp0.i.x = __shfl_up(val0.i.x,  1, TSUB_SORT);  tmp0.i.y = __shfl_up(val0.i.y,  1, TSUB_SORT);
  tmp0.i.z = __shfl_up(val0.i.z,  1, TSUB_SORT);  tmp0.i.w = __shfl_up(val0.i.w,  1, TSUB_SORT);
  tmp1.i.x = __shfl_up(val1.i.x,  1, TSUB_SORT);  tmp1.i.y = __shfl_up(val1.i.y,  1, TSUB_SORT);
  tmp1.i.z = __shfl_up(val1.i.z,  1, TSUB_SORT);  tmp1.i.w = __shfl_up(val1.i.w,  1, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.x = __shfl_up(val2.i.x,  1, TSUB_SORT);  tmp2.i.y = __shfl_up(val2.i.y,  1, TSUB_SORT);
  tmp2.i.z = __shfl_up(val2.i.z,  1, TSUB_SORT);  tmp2.i.w = __shfl_up(val2.i.w,  1, TSUB_SORT);
  tmp3.i.x = __shfl_up(val3.i.x,  1, TSUB_SORT);  tmp3.i.y = __shfl_up(val3.i.y,  1, TSUB_SORT);
  tmp3.i.z = __shfl_up(val3.i.z,  1, TSUB_SORT);  tmp3.i.w = __shfl_up(val3.i.w,  1, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  if( lane >=  1 ){
    val0.u.x += tmp0.u.x;    val0.u.y += tmp0.u.y;    val0.u.z += tmp0.u.z;    val0.u.w += tmp0.u.w;
    val1.u.x += tmp1.u.x;    val1.u.y += tmp1.u.y;    val1.u.z += tmp1.u.z;    val1.u.w += tmp1.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    val2.u.x += tmp2.u.x;    val2.u.y += tmp2.u.y;    val2.u.z += tmp2.u.z;    val2.u.w += tmp2.u.w;
    val3.u.x += tmp3.u.x;    val3.u.y += tmp3.u.y;    val3.u.z += tmp3.u.z;    val3.u.w += tmp3.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  }
#   if  TSUB_SORT >=  4
  tmp0.i.x = __shfl_up(val0.i.x,  2, TSUB_SORT);  tmp0.i.y = __shfl_up(val0.i.y,  2, TSUB_SORT);
  tmp0.i.z = __shfl_up(val0.i.z,  2, TSUB_SORT);  tmp0.i.w = __shfl_up(val0.i.w,  2, TSUB_SORT);
  tmp1.i.x = __shfl_up(val1.i.x,  2, TSUB_SORT);  tmp1.i.y = __shfl_up(val1.i.y,  2, TSUB_SORT);
  tmp1.i.z = __shfl_up(val1.i.z,  2, TSUB_SORT);  tmp1.i.w = __shfl_up(val1.i.w,  2, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.x = __shfl_up(val2.i.x,  2, TSUB_SORT);  tmp2.i.y = __shfl_up(val2.i.y,  2, TSUB_SORT);
  tmp2.i.z = __shfl_up(val2.i.z,  2, TSUB_SORT);  tmp2.i.w = __shfl_up(val2.i.w,  2, TSUB_SORT);
  tmp3.i.x = __shfl_up(val3.i.x,  2, TSUB_SORT);  tmp3.i.y = __shfl_up(val3.i.y,  2, TSUB_SORT);
  tmp3.i.z = __shfl_up(val3.i.z,  2, TSUB_SORT);  tmp3.i.w = __shfl_up(val3.i.w,  2, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  if( lane >=  2 ){
    val0.u.x += tmp0.u.x;    val0.u.y += tmp0.u.y;    val0.u.z += tmp0.u.z;    val0.u.w += tmp0.u.w;
    val1.u.x += tmp1.u.x;    val1.u.y += tmp1.u.y;    val1.u.z += tmp1.u.z;    val1.u.w += tmp1.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    val2.u.x += tmp2.u.x;    val2.u.y += tmp2.u.y;    val2.u.z += tmp2.u.z;    val2.u.w += tmp2.u.w;
    val3.u.x += tmp3.u.x;    val3.u.y += tmp3.u.y;    val3.u.z += tmp3.u.z;    val3.u.w += tmp3.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  }
#   if  TSUB_SORT >=  8
  tmp0.i.x = __shfl_up(val0.i.x,  4, TSUB_SORT);  tmp0.i.y = __shfl_up(val0.i.y,  4, TSUB_SORT);
  tmp0.i.z = __shfl_up(val0.i.z,  4, TSUB_SORT);  tmp0.i.w = __shfl_up(val0.i.w,  4, TSUB_SORT);
  tmp1.i.x = __shfl_up(val1.i.x,  4, TSUB_SORT);  tmp1.i.y = __shfl_up(val1.i.y,  4, TSUB_SORT);
  tmp1.i.z = __shfl_up(val1.i.z,  4, TSUB_SORT);  tmp1.i.w = __shfl_up(val1.i.w,  4, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.x = __shfl_up(val2.i.x,  4, TSUB_SORT);  tmp2.i.y = __shfl_up(val2.i.y,  4, TSUB_SORT);
  tmp2.i.z = __shfl_up(val2.i.z,  4, TSUB_SORT);  tmp2.i.w = __shfl_up(val2.i.w,  4, TSUB_SORT);
  tmp3.i.x = __shfl_up(val3.i.x,  4, TSUB_SORT);  tmp3.i.y = __shfl_up(val3.i.y,  4, TSUB_SORT);
  tmp3.i.z = __shfl_up(val3.i.z,  4, TSUB_SORT);  tmp3.i.w = __shfl_up(val3.i.w,  4, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  if( lane >=  4 ){
    val0.u.x += tmp0.u.x;    val0.u.y += tmp0.u.y;    val0.u.z += tmp0.u.z;    val0.u.w += tmp0.u.w;
    val1.u.x += tmp1.u.x;    val1.u.y += tmp1.u.y;    val1.u.z += tmp1.u.z;    val1.u.w += tmp1.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    val2.u.x += tmp2.u.x;    val2.u.y += tmp2.u.y;    val2.u.z += tmp2.u.z;    val2.u.w += tmp2.u.w;
    val3.u.x += tmp3.u.x;    val3.u.y += tmp3.u.y;    val3.u.z += tmp3.u.z;    val3.u.w += tmp3.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  }
#   if  TSUB_SORT >= 16
  tmp0.i.x = __shfl_up(val0.i.x,  8, TSUB_SORT);  tmp0.i.y = __shfl_up(val0.i.y,  8, TSUB_SORT);
  tmp0.i.z = __shfl_up(val0.i.z,  8, TSUB_SORT);  tmp0.i.w = __shfl_up(val0.i.w,  8, TSUB_SORT);
  tmp1.i.x = __shfl_up(val1.i.x,  8, TSUB_SORT);  tmp1.i.y = __shfl_up(val1.i.y,  8, TSUB_SORT);
  tmp1.i.z = __shfl_up(val1.i.z,  8, TSUB_SORT);  tmp1.i.w = __shfl_up(val1.i.w,  8, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.x = __shfl_up(val2.i.x,  8, TSUB_SORT);  tmp2.i.y = __shfl_up(val2.i.y,  8, TSUB_SORT);
  tmp2.i.z = __shfl_up(val2.i.z,  8, TSUB_SORT);  tmp2.i.w = __shfl_up(val2.i.w,  8, TSUB_SORT);
  tmp3.i.x = __shfl_up(val3.i.x,  8, TSUB_SORT);  tmp3.i.y = __shfl_up(val3.i.y,  8, TSUB_SORT);
  tmp3.i.z = __shfl_up(val3.i.z,  8, TSUB_SORT);  tmp3.i.w = __shfl_up(val3.i.w,  8, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  if( lane >=  8 ){
    val0.u.x += tmp0.u.x;    val0.u.y += tmp0.u.y;    val0.u.z += tmp0.u.z;    val0.u.w += tmp0.u.w;
    val1.u.x += tmp1.u.x;    val1.u.y += tmp1.u.y;    val1.u.z += tmp1.u.z;    val1.u.w += tmp1.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    val2.u.x += tmp2.u.x;    val2.u.y += tmp2.u.y;    val2.u.z += tmp2.u.z;    val2.u.w += tmp2.u.w;
    val3.u.x += tmp3.u.x;    val3.u.y += tmp3.u.y;    val3.u.z += tmp3.u.z;    val3.u.w += tmp3.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  }
#   if  TSUB_SORT == 32
  tmp0.i.x = __shfl_up(val0.i.x, 16, TSUB_SORT);  tmp0.i.y = __shfl_up(val0.i.y, 16, TSUB_SORT);
  tmp0.i.z = __shfl_up(val0.i.z, 16, TSUB_SORT);  tmp0.i.w = __shfl_up(val0.i.w, 16, TSUB_SORT);
  tmp1.i.x = __shfl_up(val1.i.x, 16, TSUB_SORT);  tmp1.i.y = __shfl_up(val1.i.y, 16, TSUB_SORT);
  tmp1.i.z = __shfl_up(val1.i.z, 16, TSUB_SORT);  tmp1.i.w = __shfl_up(val1.i.w, 16, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.x = __shfl_up(val2.i.x, 16, TSUB_SORT);  tmp2.i.y = __shfl_up(val2.i.y, 16, TSUB_SORT);
  tmp2.i.z = __shfl_up(val2.i.z, 16, TSUB_SORT);  tmp2.i.w = __shfl_up(val2.i.w, 16, TSUB_SORT);
  tmp3.i.x = __shfl_up(val3.i.x, 16, TSUB_SORT);  tmp3.i.y = __shfl_up(val3.i.y, 16, TSUB_SORT);
  tmp3.i.z = __shfl_up(val3.i.z, 16, TSUB_SORT);  tmp3.i.w = __shfl_up(val3.i.w, 16, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  if( lane >= 16 ){
    val0.u.x += tmp0.u.x;    val0.u.y += tmp0.u.y;    val0.u.z += tmp0.u.z;    val0.u.w += tmp0.u.w;
    val1.u.x += tmp1.u.x;    val1.u.y += tmp1.u.y;    val1.u.z += tmp1.u.z;    val1.u.w += tmp1.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    val2.u.x += tmp2.u.x;    val2.u.y += tmp2.u.y;    val2.u.z += tmp2.u.z;    val2.u.w += tmp2.u.w;
    val3.u.x += tmp3.u.x;    val3.u.y += tmp3.u.y;    val3.u.z += tmp3.u.z;    val3.u.w += tmp3.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  }
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
  /* return calculated inclusive prefix sum */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  *ret  = val. u;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  *ret0 = val0.u;
  *ret1 = val1.u;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  *ret2 = val2.u;
  *ret3 = val3.u;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
  /* load index */
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  smem[hidx                ] = psum;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  smem[hidx                ] = psum.x;
  smem[hidx +     TSUB_SORT] = psum.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  smem[hidx + 2 * TSUB_SORT] = psum.z;
  smem[hidx + 3 * TSUB_SORT] = psum.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  smem[hidx                 ] = psum. x;  smem[hidx +      TSUB_SORT] = psum. y;  smem[hidx +  2 * TSUB_SORT] = psum. z;  smem[hidx +  3 * TSUB_SORT] = psum. w;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  smem[hidx                 ] = psum0.x;  smem[hidx +      TSUB_SORT] = psum0.y;  smem[hidx +  2 * TSUB_SORT] = psum0.z;  smem[hidx +  3 * TSUB_SORT] = psum0.w;
  smem[hidx +  4 * TSUB_SORT] = psum1.x;  smem[hidx +  5 * TSUB_SORT] = psum1.y;  smem[hidx +  6 * TSUB_SORT] = psum1.z;  smem[hidx +  7 * TSUB_SORT] = psum1.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  smem[hidx +  8 * TSUB_SORT] = psum2.x;  smem[hidx +  9 * TSUB_SORT] = psum2.y;  smem[hidx + 10 * TSUB_SORT] = psum2.z;  smem[hidx + 11 * TSUB_SORT] = psum2.w;
  smem[hidx + 12 * TSUB_SORT] = psum3.x;  smem[hidx + 13 * TSUB_SORT] = psum3.y;  smem[hidx + 14 * TSUB_SORT] = psum3.z;  smem[hidx + 15 * TSUB_SORT] = psum3.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
#pragma unroll
  for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
    const int idx = hidx + TSUB_SORT * ii;
    smem[idx] += smem[idx -  1] & m0;
#   if  TSUB_SORT >=  4
    smem[idx] += smem[idx -  2] & m1;
#   if  TSUB_SORT >=  8
    smem[idx] += smem[idx -  4] & m2;
#   if  TSUB_SORT >= 16
    smem[idx] += smem[idx -  8] & m3;
#   if  TSUB_SORT == 32
    smem[idx] += smem[idx - 16] & m4;
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
  }
#endif//TSUB_SORT >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
  if( lane >=  1 )
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      const int idx = hidx + TSUB_SORT * ii;      smem[idx] += smem[idx -  1];    }
#   if  TSUB_SORT >=  4
  if( lane >=  2 )
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      const int idx = hidx + TSUB_SORT * ii;      smem[idx] += smem[idx -  2];    }
#   if  TSUB_SORT >=  8
  if( lane >=  4 )
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      const int idx = hidx + TSUB_SORT * ii;      smem[idx] += smem[idx -  4];    }
#   if  TSUB_SORT >= 16
  if( lane >=  8 )
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      const int idx = hidx + TSUB_SORT * ii;      smem[idx] += smem[idx -  8];    }
#   if  TSUB_SORT == 32
  if( lane >= 16 )
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      const int idx = hidx + TSUB_SORT * ii;      smem[idx] += smem[idx - 16];    }
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#endif//USE_MASK_INSTEAD_OF_IF
  /* return calculated inclusive prefix sum */
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  *ret   = smem[hidx                ];
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  ret->x = smem[hidx                ];
  ret->y = smem[hidx +     TSUB_SORT];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  ret->z = smem[hidx + 2 * TSUB_SORT];
  ret->w = smem[hidx + 3 * TSUB_SORT];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  ret ->x = smem[hidx                 ];  ret ->y = smem[hidx +      TSUB_SORT];  ret ->z = smem[hidx +  2 * TSUB_SORT];  ret ->w = smem[hidx +  3 * TSUB_SORT];
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  ret0->x = smem[hidx                 ];  ret0->y = smem[hidx +      TSUB_SORT];  ret0->z = smem[hidx +  2 * TSUB_SORT];  ret0->w = smem[hidx +  3 * TSUB_SORT];
  ret1->x = smem[hidx +  4 * TSUB_SORT];  ret1->y = smem[hidx +  5 * TSUB_SORT];  ret1->z = smem[hidx +  6 * TSUB_SORT];  ret1->w = smem[hidx +  7 * TSUB_SORT];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  ret2->x = smem[hidx +  8 * TSUB_SORT];  ret2->y = smem[hidx +  9 * TSUB_SORT];  ret2->z = smem[hidx + 10 * TSUB_SORT];  ret2->w = smem[hidx + 11 * TSUB_SORT];
  ret3->x = smem[hidx + 12 * TSUB_SORT];  ret3->y = smem[hidx + 13 * TSUB_SORT];  ret3->z = smem[hidx + 14 * TSUB_SORT];  ret3->w = smem[hidx + 15 * TSUB_SORT];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* parallel prefix sum within a warp (use implicit synchronization) */
/* type of prefix sum is inclusive */
/* NOTE: implicit synchronization within 32 threads (a warp) is assumed */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_PREFIX_SUM_WARP
(const int lane
#   if  RADIX_SORT_CHECK_BITS <= 2
 , const RADIX_SORT_UINT psum, RADIX_SORT_UINT * RESTRICT ret
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
 , const uint4 psum , uint4 * RESTRICT ret
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
 , const uint4 psum0, uint4 * RESTRICT ret0, const uint4 psum1, uint4 * RESTRICT ret1
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
 , const uint4 psum2, uint4 * RESTRICT ret2, const uint4 psum3, uint4 * RESTRICT ret3
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
 , const int tidx, volatile uint_int * RESTRICT smem
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
#   if  RADIX_SORT_CHECK_BITS <= 2
  /* load index */
  RADIX_SORT_UINT_INT val;  val.u = psum;
  RADIX_SORT_UINT_INT tmp;
  /* calculate inclusive prefix sum */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
#ifdef  USE_MASK_INSTEAD_OF_IF
/* #   if  warpSize >=  2 */  tmp.i = __shfl_up(val.i,  1, warpSize) & m0;  val.u += tmp.u;
/* #   if  warpSize >=  4 */  tmp.i = __shfl_up(val.i,  2, warpSize) & m1;  val.u += tmp.u;
/* #   if  warpSize >=  8 */  tmp.i = __shfl_up(val.i,  4, warpSize) & m2;  val.u += tmp.u;
/* #   if  warpSize >= 16 */  tmp.i = __shfl_up(val.i,  8, warpSize) & m3;  val.u += tmp.u;
/* #   if  warpSize == 32 */  tmp.i = __shfl_up(val.i, 16, warpSize) & m4;  val.u += tmp.u;
#else///USE_MASK_INSTEAD_OF_IF
/* #   if  warpSize >=  2 */  tmp.i = __shfl_up(val.i,  1, warpSize);  if( lane >=  1 )    val.u += tmp.u;
/* #   if  warpSize >=  4 */  tmp.i = __shfl_up(val.i,  2, warpSize);  if( lane >=  2 )    val.u += tmp.u;
/* #   if  warpSize >=  8 */  tmp.i = __shfl_up(val.i,  4, warpSize);  if( lane >=  4 )    val.u += tmp.u;
/* #   if  warpSize >= 16 */  tmp.i = __shfl_up(val.i,  8, warpSize);  if( lane >=  8 )    val.u += tmp.u;
/* #   if  warpSize == 32 */  tmp.i = __shfl_up(val.i, 16, warpSize);  if( lane >= 16 )    val.u += tmp.u;
#endif//USE_MASK_INSTEAD_OF_IF
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
#ifdef  USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */
  tmp.i.x = __shfl_up(val.i.x,  1, warpSize) & m0;  tmp.i.y = __shfl_up(val.i.y,  1, warpSize) & m0;  val.u.x += tmp.u.x;  val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  1, warpSize) & m0;  tmp.i.w = __shfl_up(val.i.w,  1, warpSize) & m0;  val.u.z += tmp.u.z;  val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  /* #   if  warpSize >=  4 */
  tmp.i.x = __shfl_up(val.i.x,  2, warpSize) & m1;  tmp.i.y = __shfl_up(val.i.y,  2, warpSize) & m1;  val.u.x += tmp.u.x;  val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  2, warpSize) & m1;  tmp.i.w = __shfl_up(val.i.w,  2, warpSize) & m1;  val.u.z += tmp.u.z;  val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  /* #   if  warpSize >=  8 */
  tmp.i.x = __shfl_up(val.i.x,  4, warpSize) & m2;  tmp.i.y = __shfl_up(val.i.y,  4, warpSize) & m2;  val.u.x += tmp.u.x;  val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  4, warpSize) & m2;  tmp.i.w = __shfl_up(val.i.w,  4, warpSize) & m2;  val.u.z += tmp.u.z;  val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  /* #   if  warpSize >= 16 */
  tmp.i.x = __shfl_up(val.i.x,  8, warpSize) & m3;  tmp.i.y = __shfl_up(val.i.y,  8, warpSize) & m3;  val.u.x += tmp.u.x;  val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  8, warpSize) & m3;  tmp.i.w = __shfl_up(val.i.w,  8, warpSize) & m3;  val.u.z += tmp.u.z;  val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  /* #   if  warpSize == 32 */
  tmp.i.x = __shfl_up(val.i.x, 16, warpSize) & m4;  tmp.i.y = __shfl_up(val.i.y, 16, warpSize) & m4;  val.u.x += tmp.u.x;  val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z, 16, warpSize) & m4;  tmp.i.w = __shfl_up(val.i.w, 16, warpSize) & m4;  val.u.z += tmp.u.z;  val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
#else///USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */
  tmp.i.x = __shfl_up(val.i.x,  1, warpSize);  tmp.i.y = __shfl_up(val.i.y,  1, warpSize);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  1, warpSize);  tmp.i.w = __shfl_up(val.i.w,  1, warpSize);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  if( lane >=  1 ){
    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  }/* if( lane >=  1 ){ */
  /* #   if  warpSize >=  4 */
  tmp.i.x = __shfl_up(val.i.x,  2, warpSize);  tmp.i.y = __shfl_up(val.i.y,  2, warpSize);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  2, warpSize);  tmp.i.w = __shfl_up(val.i.w,  2, warpSize);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  if( lane >=  2 ){
    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  }/* if( lane >=  2 ){ */
  /* #   if  warpSize >=  8 */
  tmp.i.x = __shfl_up(val.i.x,  4, warpSize);  tmp.i.y = __shfl_up(val.i.y,  4, warpSize);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  4, warpSize);  tmp.i.w = __shfl_up(val.i.w,  4, warpSize);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  if( lane >=  4 ){
    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  }/* if( lane >=  4 ){ */
  /* #   if  warpSize >= 16 */
  tmp.i.x = __shfl_up(val.i.x,  8, warpSize);  tmp.i.y = __shfl_up(val.i.y,  8, warpSize);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z,  8, warpSize);  tmp.i.w = __shfl_up(val.i.w,  8, warpSize);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  if( lane >=  8 ){
    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  }/* if( lane >=  8 ){ */
  /* #   if  warpSize == 32 */
  tmp.i.x = __shfl_up(val.i.x, 16, warpSize);  tmp.i.y = __shfl_up(val.i.y, 16, warpSize);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  tmp.i.z = __shfl_up(val.i.z, 16, warpSize);  tmp.i.w = __shfl_up(val.i.w, 16, warpSize);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  if( lane >= 16 ){
    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
  }/* if( lane >= 16 ){ */
#endif//USE_MASK_INSTEAD_OF_IF
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
  /* return calculated inclusive prefix sum */
  *ret = val.u;
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  /* load index */
  uint4_int4 val;  val.u = psum;
  uint4_int4 tmp;
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */
  tmp.i.z = __shfl_up(val.i.x,  1, warpSize) & m0;  tmp.i.w = __shfl_up(val.i.y,  1, warpSize) & m0;
  tmp.i.x = __shfl_up(val.i.z,  1, warpSize) & m0;  tmp.i.y = __shfl_up(val.i.w,  1, warpSize) & m0;
  val.u.x += tmp.u.z;  val.u.y += tmp.u.w;  val.u.z += tmp.u.x;  val.u.w += tmp.u.y;
  /* #   if  warpSize >=  4 */
  tmp.i.z = __shfl_up(val.i.x,  2, warpSize) & m1;  tmp.i.w = __shfl_up(val.i.y,  2, warpSize) & m1;
  tmp.i.x = __shfl_up(val.i.z,  2, warpSize) & m1;  tmp.i.y = __shfl_up(val.i.w,  2, warpSize) & m1;
  val.u.x += tmp.u.z;  val.u.y += tmp.u.w;  val.u.z += tmp.u.x;  val.u.w += tmp.u.y;
  /* #   if  warpSize >=  8 */
  tmp.i.z = __shfl_up(val.i.x,  4, warpSize) & m2;  tmp.i.w = __shfl_up(val.i.y,  4, warpSize) & m2;
  tmp.i.x = __shfl_up(val.i.z,  4, warpSize) & m2;  tmp.i.y = __shfl_up(val.i.w,  4, warpSize) & m2;
  val.u.x += tmp.u.z;  val.u.y += tmp.u.w;  val.u.z += tmp.u.x;  val.u.w += tmp.u.y;
  /* #   if  warpSize >= 16 */
  tmp.i.z = __shfl_up(val.i.x,  8, warpSize) & m3;  tmp.i.w = __shfl_up(val.i.y,  8, warpSize) & m3;
  tmp.i.x = __shfl_up(val.i.z,  8, warpSize) & m3;  tmp.i.y = __shfl_up(val.i.w,  8, warpSize) & m3;
  val.u.x += tmp.u.z;  val.u.y += tmp.u.w;  val.u.z += tmp.u.x;  val.u.w += tmp.u.y;
  /* #   if  warpSize >= 32 */
  tmp.i.z = __shfl_up(val.i.x, 16, warpSize) & m4;  tmp.i.w = __shfl_up(val.i.y, 16, warpSize) & m4;
  tmp.i.x = __shfl_up(val.i.z, 16, warpSize) & m4;  tmp.i.y = __shfl_up(val.i.w, 16, warpSize) & m4;
  val.u.x += tmp.u.z;  val.u.y += tmp.u.w;  val.u.z += tmp.u.x;  val.u.w += tmp.u.y;
#else///USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */
  tmp.i.x = __shfl_up(val.i.x,  1, warpSize);  tmp.i.y = __shfl_up(val.i.y,  1, warpSize);
  tmp.i.z = __shfl_up(val.i.z,  1, warpSize);  tmp.i.w = __shfl_up(val.i.w,  1, warpSize);
  if( lane >=  1 ){    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;  }
  /* #   if  warpSize >=  4 */
  tmp.i.x = __shfl_up(val.i.x,  2, warpSize);  tmp.i.y = __shfl_up(val.i.y,  2, warpSize);
  tmp.i.z = __shfl_up(val.i.z,  2, warpSize);  tmp.i.w = __shfl_up(val.i.w,  2, warpSize);
  if( lane >=  2 ){    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;  }
  /* #   if  warpSize >=  8 */
  tmp.i.x = __shfl_up(val.i.x,  4, warpSize);  tmp.i.y = __shfl_up(val.i.y,  4, warpSize);
  tmp.i.z = __shfl_up(val.i.z,  4, warpSize);  tmp.i.w = __shfl_up(val.i.w,  4, warpSize);
  if( lane >=  4 ){    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;  }
  /* #   if  warpSize >= 16 */
  tmp.i.x = __shfl_up(val.i.x,  8, warpSize);  tmp.i.y = __shfl_up(val.i.y,  8, warpSize);
  tmp.i.z = __shfl_up(val.i.z,  8, warpSize);  tmp.i.w = __shfl_up(val.i.w,  8, warpSize);
  if( lane >=  8 ){    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;  }
  /* #   if  warpSize >= 32 */
  tmp.i.x = __shfl_up(val.i.x, 16, warpSize);  tmp.i.y = __shfl_up(val.i.y, 16, warpSize);
  tmp.i.z = __shfl_up(val.i.z, 16, warpSize);  tmp.i.w = __shfl_up(val.i.w, 16, warpSize);
  if( lane >= 16 ){    val.u.x += tmp.u.x;    val.u.y += tmp.u.y;    val.u.z += tmp.u.z;    val.u.w += tmp.u.w;  }
#endif//USE_MASK_INSTEAD_OF_IF
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  /* load index */
  uint4_int4 val0, val1;  val0.u = psum0;  val1.u = psum1;
  uint4_int4 tmp0, tmp1;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  uint4_int4 val2, val3;  val2.u = psum2;  val3.u = psum3;
  uint4_int4 tmp2, tmp3;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */
  tmp0.i.z = __shfl_up(val0.i.x,  1, warpSize) & m0;  tmp1.i.z = __shfl_up(val1.i.x,  1, warpSize) & m0;
  tmp0.i.w = __shfl_up(val0.i.y,  1, warpSize) & m0;  tmp1.i.w = __shfl_up(val1.i.y,  1, warpSize) & m0;
  tmp0.i.x = __shfl_up(val0.i.z,  1, warpSize) & m0;  tmp1.i.x = __shfl_up(val1.i.z,  1, warpSize) & m0;
  tmp0.i.y = __shfl_up(val0.i.w,  1, warpSize) & m0;  tmp1.i.y = __shfl_up(val1.i.w,  1, warpSize) & m0;
  val0.u.x += tmp0.u.z;  val1.u.x += tmp1.u.z;  val0.u.y += tmp0.u.w;  val1.u.y += tmp1.u.w;
  val0.u.z += tmp0.u.x;  val1.u.z += tmp1.u.x;  val0.u.w += tmp0.u.y;  val1.u.w += tmp1.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.z = __shfl_up(val2.i.x,  1, warpSize) & m0;  tmp3.i.z = __shfl_up(val3.i.x,  1, warpSize) & m0;
  tmp2.i.w = __shfl_up(val2.i.y,  1, warpSize) & m0;  tmp3.i.w = __shfl_up(val3.i.y,  1, warpSize) & m0;
  tmp2.i.x = __shfl_up(val2.i.z,  1, warpSize) & m0;  tmp3.i.x = __shfl_up(val3.i.z,  1, warpSize) & m0;
  tmp2.i.y = __shfl_up(val2.i.w,  1, warpSize) & m0;  tmp3.i.y = __shfl_up(val3.i.w,  1, warpSize) & m0;
  val2.u.x += tmp2.u.z;  val3.u.x += tmp3.u.z;  val2.u.y += tmp2.u.w;  val3.u.y += tmp3.u.w;
  val2.u.z += tmp2.u.x;  val3.u.z += tmp3.u.x;  val2.u.w += tmp2.u.y;  val3.u.w += tmp3.u.y;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  /* #   if  warpSize >=  4 */
  tmp0.i.z = __shfl_up(val0.i.x,  2, warpSize) & m1;  tmp1.i.z = __shfl_up(val1.i.x,  2, warpSize) & m1;
  tmp0.i.w = __shfl_up(val0.i.y,  2, warpSize) & m1;  tmp1.i.w = __shfl_up(val1.i.y,  2, warpSize) & m1;
  tmp0.i.x = __shfl_up(val0.i.z,  2, warpSize) & m1;  tmp1.i.x = __shfl_up(val1.i.z,  2, warpSize) & m1;
  tmp0.i.y = __shfl_up(val0.i.w,  2, warpSize) & m1;  tmp1.i.y = __shfl_up(val1.i.w,  2, warpSize) & m1;
  val0.u.x += tmp0.u.z;  val1.u.x += tmp1.u.z;  val0.u.y += tmp0.u.w;  val1.u.y += tmp1.u.w;
  val0.u.z += tmp0.u.x;  val1.u.z += tmp1.u.x;  val0.u.w += tmp0.u.y;  val1.u.w += tmp1.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.z = __shfl_up(val2.i.x,  2, warpSize) & m1;  tmp3.i.z = __shfl_up(val3.i.x,  2, warpSize) & m1;
  tmp2.i.w = __shfl_up(val2.i.y,  2, warpSize) & m1;  tmp3.i.w = __shfl_up(val3.i.y,  2, warpSize) & m1;
  tmp2.i.x = __shfl_up(val2.i.z,  2, warpSize) & m1;  tmp3.i.x = __shfl_up(val3.i.z,  2, warpSize) & m1;
  tmp2.i.y = __shfl_up(val2.i.w,  2, warpSize) & m1;  tmp3.i.y = __shfl_up(val3.i.w,  2, warpSize) & m1;
  val2.u.x += tmp2.u.z;  val3.u.x += tmp3.u.z;  val2.u.y += tmp2.u.w;  val3.u.y += tmp3.u.w;
  val2.u.z += tmp2.u.x;  val3.u.z += tmp3.u.x;  val2.u.w += tmp2.u.y;  val3.u.w += tmp3.u.y;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  /* #   if  warpSize >=  8 */
  tmp0.i.z = __shfl_up(val0.i.x,  4, warpSize) & m2;  tmp1.i.z = __shfl_up(val1.i.x,  4, warpSize) & m2;
  tmp0.i.w = __shfl_up(val0.i.y,  4, warpSize) & m2;  tmp1.i.w = __shfl_up(val1.i.y,  4, warpSize) & m2;
  tmp0.i.x = __shfl_up(val0.i.z,  4, warpSize) & m2;  tmp1.i.x = __shfl_up(val1.i.z,  4, warpSize) & m2;
  tmp0.i.y = __shfl_up(val0.i.w,  4, warpSize) & m2;  tmp1.i.y = __shfl_up(val1.i.w,  4, warpSize) & m2;
  val0.u.x += tmp0.u.z;  val1.u.x += tmp1.u.z;  val0.u.y += tmp0.u.w;  val1.u.y += tmp1.u.w;
  val0.u.z += tmp0.u.x;  val1.u.z += tmp1.u.x;  val0.u.w += tmp0.u.y;  val1.u.w += tmp1.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.z = __shfl_up(val2.i.x,  4, warpSize) & m2;  tmp3.i.z = __shfl_up(val3.i.x,  4, warpSize) & m2;
  tmp2.i.w = __shfl_up(val2.i.y,  4, warpSize) & m2;  tmp3.i.w = __shfl_up(val3.i.y,  4, warpSize) & m2;
  tmp2.i.x = __shfl_up(val2.i.z,  4, warpSize) & m2;  tmp3.i.x = __shfl_up(val3.i.z,  4, warpSize) & m2;
  tmp2.i.y = __shfl_up(val2.i.w,  4, warpSize) & m2;  tmp3.i.y = __shfl_up(val3.i.w,  4, warpSize) & m2;
  val2.u.x += tmp2.u.z;  val3.u.x += tmp3.u.z;  val2.u.y += tmp2.u.w;  val3.u.y += tmp3.u.w;
  val2.u.z += tmp2.u.x;  val3.u.z += tmp3.u.x;  val2.u.w += tmp2.u.y;  val3.u.w += tmp3.u.y;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  /* #   if  warpSize >= 16 */
  tmp0.i.z = __shfl_up(val0.i.x,  8, warpSize) & m3;  tmp1.i.z = __shfl_up(val1.i.x,  8, warpSize) & m3;
  tmp0.i.w = __shfl_up(val0.i.y,  8, warpSize) & m3;  tmp1.i.w = __shfl_up(val1.i.y,  8, warpSize) & m3;
  tmp0.i.x = __shfl_up(val0.i.z,  8, warpSize) & m3;  tmp1.i.x = __shfl_up(val1.i.z,  8, warpSize) & m3;
  tmp0.i.y = __shfl_up(val0.i.w,  8, warpSize) & m3;  tmp1.i.y = __shfl_up(val1.i.w,  8, warpSize) & m3;
  val0.u.x += tmp0.u.z;  val1.u.x += tmp1.u.z;  val0.u.y += tmp0.u.w;  val1.u.y += tmp1.u.w;
  val0.u.z += tmp0.u.x;  val1.u.z += tmp1.u.x;  val0.u.w += tmp0.u.y;  val1.u.w += tmp1.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.z = __shfl_up(val2.i.x,  8, warpSize) & m3;  tmp3.i.z = __shfl_up(val3.i.x,  8, warpSize) & m3;
  tmp2.i.w = __shfl_up(val2.i.y,  8, warpSize) & m3;  tmp3.i.w = __shfl_up(val3.i.y,  8, warpSize) & m3;
  tmp2.i.x = __shfl_up(val2.i.z,  8, warpSize) & m3;  tmp3.i.x = __shfl_up(val3.i.z,  8, warpSize) & m3;
  tmp2.i.y = __shfl_up(val2.i.w,  8, warpSize) & m3;  tmp3.i.y = __shfl_up(val3.i.w,  8, warpSize) & m3;
  val2.u.x += tmp2.u.z;  val3.u.x += tmp3.u.z;  val2.u.y += tmp2.u.w;  val3.u.y += tmp3.u.w;
  val2.u.z += tmp2.u.x;  val3.u.z += tmp3.u.x;  val2.u.w += tmp2.u.y;  val3.u.w += tmp3.u.y;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  /* #   if  warpSize >= 32 */
  tmp0.i.z = __shfl_up(val0.i.x, 16, warpSize) & m4;  tmp1.i.z = __shfl_up(val1.i.x, 16, warpSize) & m4;
  tmp0.i.w = __shfl_up(val0.i.y, 16, warpSize) & m4;  tmp1.i.w = __shfl_up(val1.i.y, 16, warpSize) & m4;
  tmp0.i.x = __shfl_up(val0.i.z, 16, warpSize) & m4;  tmp1.i.x = __shfl_up(val1.i.z, 16, warpSize) & m4;
  tmp0.i.y = __shfl_up(val0.i.w, 16, warpSize) & m4;  tmp1.i.y = __shfl_up(val1.i.w, 16, warpSize) & m4;
  val0.u.x += tmp0.u.z;  val1.u.x += tmp1.u.z;  val0.u.y += tmp0.u.w;  val1.u.y += tmp1.u.w;
  val0.u.z += tmp0.u.x;  val1.u.z += tmp1.u.x;  val0.u.w += tmp0.u.y;  val1.u.w += tmp1.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.z = __shfl_up(val2.i.x, 16, warpSize) & m4;  tmp3.i.z = __shfl_up(val3.i.x, 16, warpSize) & m4;
  tmp2.i.w = __shfl_up(val2.i.y, 16, warpSize) & m4;  tmp3.i.w = __shfl_up(val3.i.y, 16, warpSize) & m4;
  tmp2.i.x = __shfl_up(val2.i.z, 16, warpSize) & m4;  tmp3.i.x = __shfl_up(val3.i.z, 16, warpSize) & m4;
  tmp2.i.y = __shfl_up(val2.i.w, 16, warpSize) & m4;  tmp3.i.y = __shfl_up(val3.i.w, 16, warpSize) & m4;
  val2.u.x += tmp2.u.z;  val3.u.x += tmp3.u.z;  val2.u.y += tmp2.u.w;  val3.u.y += tmp3.u.w;
  val2.u.z += tmp2.u.x;  val3.u.z += tmp3.u.x;  val2.u.w += tmp2.u.y;  val3.u.w += tmp3.u.y;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#else///USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */
  tmp0.i.x = __shfl_up(val0.i.x,  1, warpSize);  tmp0.i.y = __shfl_up(val0.i.y,  1, warpSize);
  tmp0.i.z = __shfl_up(val0.i.z,  1, warpSize);  tmp0.i.w = __shfl_up(val0.i.w,  1, warpSize);
  tmp1.i.x = __shfl_up(val1.i.x,  1, warpSize);  tmp1.i.y = __shfl_up(val1.i.y,  1, warpSize);
  tmp1.i.z = __shfl_up(val1.i.z,  1, warpSize);  tmp1.i.w = __shfl_up(val1.i.w,  1, warpSize);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.x = __shfl_up(val2.i.x,  1, warpSize);  tmp2.i.y = __shfl_up(val2.i.y,  1, warpSize);
  tmp2.i.z = __shfl_up(val2.i.z,  1, warpSize);  tmp2.i.w = __shfl_up(val2.i.w,  1, warpSize);
  tmp3.i.x = __shfl_up(val3.i.x,  1, warpSize);  tmp3.i.y = __shfl_up(val3.i.y,  1, warpSize);
  tmp3.i.z = __shfl_up(val3.i.z,  1, warpSize);  tmp3.i.w = __shfl_up(val3.i.w,  1, warpSize);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  if( lane >=  1 ){
    val0.u.x += tmp0.u.x;    val0.u.y += tmp0.u.y;    val0.u.z += tmp0.u.z;    val0.u.w += tmp0.u.w;
    val1.u.x += tmp1.u.x;    val1.u.y += tmp1.u.y;    val1.u.z += tmp1.u.z;    val1.u.w += tmp1.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    val2.u.x += tmp2.u.x;    val2.u.y += tmp2.u.y;    val2.u.z += tmp2.u.z;    val2.u.w += tmp2.u.w;
    val3.u.x += tmp3.u.x;    val3.u.y += tmp3.u.y;    val3.u.z += tmp3.u.z;    val3.u.w += tmp3.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  }
  /* #   if  warpSize >=  4 */
  tmp0.i.x = __shfl_up(val0.i.x,  2, warpSize);  tmp0.i.y = __shfl_up(val0.i.y,  2, warpSize);
  tmp0.i.z = __shfl_up(val0.i.z,  2, warpSize);  tmp0.i.w = __shfl_up(val0.i.w,  2, warpSize);
  tmp1.i.x = __shfl_up(val1.i.x,  2, warpSize);  tmp1.i.y = __shfl_up(val1.i.y,  2, warpSize);
  tmp1.i.z = __shfl_up(val1.i.z,  2, warpSize);  tmp1.i.w = __shfl_up(val1.i.w,  2, warpSize);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.x = __shfl_up(val2.i.x,  2, warpSize);  tmp2.i.y = __shfl_up(val2.i.y,  2, warpSize);
  tmp2.i.z = __shfl_up(val2.i.z,  2, warpSize);  tmp2.i.w = __shfl_up(val2.i.w,  2, warpSize);
  tmp3.i.x = __shfl_up(val3.i.x,  2, warpSize);  tmp3.i.y = __shfl_up(val3.i.y,  2, warpSize);
  tmp3.i.z = __shfl_up(val3.i.z,  2, warpSize);  tmp3.i.w = __shfl_up(val3.i.w,  2, warpSize);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  if( lane >=  2 ){
    val0.u.x += tmp0.u.x;    val0.u.y += tmp0.u.y;    val0.u.z += tmp0.u.z;    val0.u.w += tmp0.u.w;
    val1.u.x += tmp1.u.x;    val1.u.y += tmp1.u.y;    val1.u.z += tmp1.u.z;    val1.u.w += tmp1.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    val2.u.x += tmp2.u.x;    val2.u.y += tmp2.u.y;    val2.u.z += tmp2.u.z;    val2.u.w += tmp2.u.w;
    val3.u.x += tmp3.u.x;    val3.u.y += tmp3.u.y;    val3.u.z += tmp3.u.z;    val3.u.w += tmp3.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  }
  /* #   if  warpSize >=  8 */
  tmp0.i.x = __shfl_up(val0.i.x,  4, warpSize);  tmp0.i.y = __shfl_up(val0.i.y,  4, warpSize);
  tmp0.i.z = __shfl_up(val0.i.z,  4, warpSize);  tmp0.i.w = __shfl_up(val0.i.w,  4, warpSize);
  tmp1.i.x = __shfl_up(val1.i.x,  4, warpSize);  tmp1.i.y = __shfl_up(val1.i.y,  4, warpSize);
  tmp1.i.z = __shfl_up(val1.i.z,  4, warpSize);  tmp1.i.w = __shfl_up(val1.i.w,  4, warpSize);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.x = __shfl_up(val2.i.x,  4, warpSize);  tmp2.i.y = __shfl_up(val2.i.y,  4, warpSize);
  tmp2.i.z = __shfl_up(val2.i.z,  4, warpSize);  tmp2.i.w = __shfl_up(val2.i.w,  4, warpSize);
  tmp3.i.x = __shfl_up(val3.i.x,  4, warpSize);  tmp3.i.y = __shfl_up(val3.i.y,  4, warpSize);
  tmp3.i.z = __shfl_up(val3.i.z,  4, warpSize);  tmp3.i.w = __shfl_up(val3.i.w,  4, warpSize);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  if( lane >=  4 ){
    val0.u.x += tmp0.u.x;    val0.u.y += tmp0.u.y;    val0.u.z += tmp0.u.z;    val0.u.w += tmp0.u.w;
    val1.u.x += tmp1.u.x;    val1.u.y += tmp1.u.y;    val1.u.z += tmp1.u.z;    val1.u.w += tmp1.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    val2.u.x += tmp2.u.x;    val2.u.y += tmp2.u.y;    val2.u.z += tmp2.u.z;    val2.u.w += tmp2.u.w;
    val3.u.x += tmp3.u.x;    val3.u.y += tmp3.u.y;    val3.u.z += tmp3.u.z;    val3.u.w += tmp3.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  }
  /* #   if  warpSize >= 16 */
  tmp0.i.x = __shfl_up(val0.i.x,  8, warpSize);  tmp0.i.y = __shfl_up(val0.i.y,  8, warpSize);
  tmp0.i.z = __shfl_up(val0.i.z,  8, warpSize);  tmp0.i.w = __shfl_up(val0.i.w,  8, warpSize);
  tmp1.i.x = __shfl_up(val1.i.x,  8, warpSize);  tmp1.i.y = __shfl_up(val1.i.y,  8, warpSize);
  tmp1.i.z = __shfl_up(val1.i.z,  8, warpSize);  tmp1.i.w = __shfl_up(val1.i.w,  8, warpSize);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.x = __shfl_up(val2.i.x,  8, warpSize);  tmp2.i.y = __shfl_up(val2.i.y,  8, warpSize);
  tmp2.i.z = __shfl_up(val2.i.z,  8, warpSize);  tmp2.i.w = __shfl_up(val2.i.w,  8, warpSize);
  tmp3.i.x = __shfl_up(val3.i.x,  8, warpSize);  tmp3.i.y = __shfl_up(val3.i.y,  8, warpSize);
  tmp3.i.z = __shfl_up(val3.i.z,  8, warpSize);  tmp3.i.w = __shfl_up(val3.i.w,  8, warpSize);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  if( lane >=  8 ){
    val0.u.x += tmp0.u.x;    val0.u.y += tmp0.u.y;    val0.u.z += tmp0.u.z;    val0.u.w += tmp0.u.w;
    val1.u.x += tmp1.u.x;    val1.u.y += tmp1.u.y;    val1.u.z += tmp1.u.z;    val1.u.w += tmp1.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    val2.u.x += tmp2.u.x;    val2.u.y += tmp2.u.y;    val2.u.z += tmp2.u.z;    val2.u.w += tmp2.u.w;
    val3.u.x += tmp3.u.x;    val3.u.y += tmp3.u.y;    val3.u.z += tmp3.u.z;    val3.u.w += tmp3.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  }
  /* #   if  warpSize >= 32 */
  tmp0.i.x = __shfl_up(val0.i.x, 16, warpSize);  tmp0.i.y = __shfl_up(val0.i.y, 16, warpSize);
  tmp0.i.z = __shfl_up(val0.i.z, 16, warpSize);  tmp0.i.w = __shfl_up(val0.i.w, 16, warpSize);
  tmp1.i.x = __shfl_up(val1.i.x, 16, warpSize);  tmp1.i.y = __shfl_up(val1.i.y, 16, warpSize);
  tmp1.i.z = __shfl_up(val1.i.z, 16, warpSize);  tmp1.i.w = __shfl_up(val1.i.w, 16, warpSize);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  tmp2.i.x = __shfl_up(val2.i.x, 16, warpSize);  tmp2.i.y = __shfl_up(val2.i.y, 16, warpSize);
  tmp2.i.z = __shfl_up(val2.i.z, 16, warpSize);  tmp2.i.w = __shfl_up(val2.i.w, 16, warpSize);
  tmp3.i.x = __shfl_up(val3.i.x, 16, warpSize);  tmp3.i.y = __shfl_up(val3.i.y, 16, warpSize);
  tmp3.i.z = __shfl_up(val3.i.z, 16, warpSize);  tmp3.i.w = __shfl_up(val3.i.w, 16, warpSize);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  if( lane >= 16 ){
    val0.u.x += tmp0.u.x;    val0.u.y += tmp0.u.y;    val0.u.z += tmp0.u.z;    val0.u.w += tmp0.u.w;
    val1.u.x += tmp1.u.x;    val1.u.y += tmp1.u.y;    val1.u.z += tmp1.u.z;    val1.u.w += tmp1.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    val2.u.x += tmp2.u.x;    val2.u.y += tmp2.u.y;    val2.u.z += tmp2.u.z;    val2.u.w += tmp2.u.w;
    val3.u.x += tmp3.u.x;    val3.u.y += tmp3.u.y;    val3.u.z += tmp3.u.z;    val3.u.w += tmp3.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
  }
#endif//USE_MASK_INSTEAD_OF_IF
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
  /* return calculated inclusive prefix sum */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  *ret  = val. u;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  *ret0 = val0.u;
  *ret1 = val1.u;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  *ret2 = val2.u;
  *ret3 = val3.u;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
  /* load index */
  const int hidx = (tidx - lane) * RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD + lane;
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  smem[hidx               ].u = psum;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  smem[hidx               ].u = psum.x;
  smem[hidx +     warpSize].u = psum.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  smem[hidx + 2 * warpSize].u = psum.z;
  smem[hidx + 3 * warpSize].u = psum.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  smem[hidx                ].u = psum. x;  smem[hidx +      warpSize].u = psum. y;  smem[hidx +  2 * warpSize].u = psum. z;  smem[hidx +  3 * warpSize].u = psum. w;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  smem[hidx                ].u = psum0.x;  smem[hidx +      warpSize].u = psum0.y;  smem[hidx +  2 * warpSize].u = psum0.z;  smem[hidx +  3 * warpSize].u = psum0.w;
  smem[hidx +  4 * warpSize].u = psum1.x;  smem[hidx +  5 * warpSize].u = psum1.y;  smem[hidx +  6 * warpSize].u = psum1.z;  smem[hidx +  7 * warpSize].u = psum1.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  smem[hidx +  8 * warpSize].u = psum2.x;  smem[hidx +  9 * warpSize].u = psum2.y;  smem[hidx + 10 * warpSize].u = psum2.z;  smem[hidx + 11 * warpSize].u = psum2.w;
  smem[hidx + 12 * warpSize].u = psum3.x;  smem[hidx + 13 * warpSize].u = psum3.y;  smem[hidx + 14 * warpSize].u = psum3.z;  smem[hidx + 15 * warpSize].u = psum3.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
#pragma unroll
  for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
    const int idx = hidx + warpSize * ii;
    smem[idx].u += smem[idx -  1].u & m0;
    smem[idx].u += smem[idx -  2].u & m1;
    smem[idx].u += smem[idx -  4].u & m2;
    smem[idx].u += smem[idx -  8].u & m3;
    smem[idx].u += smem[idx - 16].u & m4;
  }
#else///USE_MASK_INSTEAD_OF_IF
  if( lane >=  1 )
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      const int idx = hidx + warpSize * ii;      smem[idx].u += smem[idx -  1].u;    }
  if( lane >=  2 )
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      const int idx = hidx + warpSize * ii;      smem[idx].u += smem[idx -  2].u;    }
  if( lane >=  4 )
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      const int idx = hidx + warpSize * ii;      smem[idx].u += smem[idx -  4].u;    }
  if( lane >=  8 )
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      const int idx = hidx + warpSize * ii;      smem[idx].u += smem[idx -  8].u;    }
  if( lane >= 16 )
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32 * RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      const int idx = hidx + warpSize * ii;      smem[idx].u += smem[idx - 16].u;    }
#endif//USE_MASK_INSTEAD_OF_IF
  /* return calculated inclusive prefix sum */
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  *ret   = smem[hidx               ].u;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  ret->x = smem[hidx               ].u;
  ret->y = smem[hidx +     warpSize].u;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD > 2
  ret->z = smem[hidx + 2 * warpSize].u;
  ret->w = smem[hidx + 3 * warpSize].u;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD > 2
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
  ret ->x = smem[hidx                ].u;  ret ->y = smem[hidx +      warpSize].u;  ret ->z = smem[hidx +  2 * warpSize].u;  ret ->w = smem[hidx +  3 * warpSize].u;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
  ret0->x = smem[hidx                ].u;  ret0->y = smem[hidx +      warpSize].u;  ret0->z = smem[hidx +  2 * warpSize].u;  ret0->w = smem[hidx +  3 * warpSize].u;
  ret1->x = smem[hidx +  4 * warpSize].u;  ret1->y = smem[hidx +  5 * warpSize].u;  ret1->z = smem[hidx +  6 * warpSize].u;  ret1->w = smem[hidx +  7 * warpSize].u;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
  ret2->x = smem[hidx +  8 * warpSize].u;  ret2->y = smem[hidx +  9 * warpSize].u;  ret2->z = smem[hidx + 10 * warpSize].u;  ret2->w = smem[hidx + 11 * warpSize].u;
  ret3->x = smem[hidx + 12 * warpSize].u;  ret3->y = smem[hidx + 13 * warpSize].u;  ret3->z = smem[hidx + 14 * warpSize].u;  ret3->w = smem[hidx + 15 * warpSize].u;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_PREFIX_SUM_BLCK
(const int tidx, const int lane, volatile uint_int * RESTRICT psum
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  const int Ndat = (1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD;
  //-----------------------------------------------------------------------
  /* 1. prefix sum within warp */
  int val = 0;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  int tmp;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
  if( (tidx - lane) < Ndat ){
    //---------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    //---------------------------------------------------------------------
    /* load index */
    val = psum[tidx].i;
    /* calculate inclusive prefix sum */
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >= 32
#ifdef  USE_MASK_INSTEAD_OF_IF
    /* #   if  warpSize >=  2 */    tmp = __shfl_up(val,  1, warpSize) & m0;    val += tmp;
    /* #   if  warpSize >=  4 */    tmp = __shfl_up(val,  2, warpSize) & m1;    val += tmp;
    /* #   if  warpSize >=  8 */    tmp = __shfl_up(val,  4, warpSize) & m2;    val += tmp;
    /* #   if  warpSize >= 16 */    tmp = __shfl_up(val,  8, warpSize) & m3;    val += tmp;
    /* #   if  warpSize >= 32 */    tmp = __shfl_up(val, 16, warpSize) & m4;    val += tmp;
#else///USE_MASK_INSTEAD_OF_IF
    /* #   if  warpSize >=  2 */    tmp = __shfl_up(val,  1, warpSize);    if( lane >=  1 )      val += tmp;
    /* #   if  warpSize >=  4 */    tmp = __shfl_up(val,  2, warpSize);    if( lane >=  2 )      val += tmp;
    /* #   if  warpSize >=  8 */    tmp = __shfl_up(val,  4, warpSize);    if( lane >=  4 )      val += tmp;
    /* #   if  warpSize >= 16 */    tmp = __shfl_up(val,  8, warpSize);    if( lane >=  8 )      val += tmp;
    /* #   if  warpSize >= 32 */    tmp = __shfl_up(val, 16, warpSize);    if( lane >= 16 )      val += tmp;
#endif//USE_MASK_INSTEAD_OF_IF
#else///if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >= 32
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  2
    tmp = __shfl_up(val, 1, Ndat) & m0;    val += tmp;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  4
    tmp = __shfl_up(val, 2, Ndat) & m1;    val += tmp;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  8
    tmp = __shfl_up(val, 4, Ndat) & m2;    val += tmp;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) == 16
    tmp = __shfl_up(val, 8, Ndat) & m3;    val += tmp;
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) == 16
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  8
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  4
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  2
    tmp = __shfl_up(val, 1, Ndat);    if( lane >= 1 )      val += tmp;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  4
    tmp = __shfl_up(val, 2, Ndat);    if( lane >= 2 )      val += tmp;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  8
    tmp = __shfl_up(val, 4, Ndat);    if( lane >= 4 )      val += tmp;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) == 16
    tmp = __shfl_up(val, 8, Ndat);    if( lane >= 8 )      val += tmp;
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) == 16
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  8
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  4
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#endif//if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >= 32
    /* return calculated inclusive prefix sum */
    psum[tidx].i = val;
    //---------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_SORT
    //---------------------------------------------------------------------
    /* calculate inclusive prefix sum */
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >= 32
#ifdef  USE_MASK_INSTEAD_OF_IF
    /* #   if  warpSize >=  2 */    psum[tidx].i += psum[tidx -  1].i & m0;
    /* #   if  warpSize >=  4 */    psum[tidx].i += psum[tidx -  2].i & m1;
    /* #   if  warpSize >=  8 */    psum[tidx].i += psum[tidx -  4].i & m2;
    /* #   if  warpSize >= 16 */    psum[tidx].i += psum[tidx -  8].i & m3;
    /* #   if  warpSize >= 32 */    psum[tidx].i += psum[tidx - 16].i & m4;
#else///USE_MASK_INSTEAD_OF_IF
    /* #   if  warpSize >=  2 */    if( lane >=  1 )      psum[tidx].i += psum[tidx -  1].i;
    /* #   if  warpSize >=  4 */    if( lane >=  2 )      psum[tidx].i += psum[tidx -  2].i;
    /* #   if  warpSize >=  8 */    if( lane >=  4 )      psum[tidx].i += psum[tidx -  4].i;
    /* #   if  warpSize >= 16 */    if( lane >=  8 )      psum[tidx].i += psum[tidx -  8].i;
    /* #   if  warpSize >= 32 */    if( lane >= 16 )      psum[tidx].i += psum[tidx - 16].i;
#endif//USE_MASK_INSTEAD_OF_IF
#else///if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >= 32
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  2
    psum[tidx].i += psum[tidx - 1].i & m0;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  4
    psum[tidx].i += psum[tidx - 2].i & m1;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  8
    psum[tidx].i += psum[tidx - 4].i & m2;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) == 16
    psum[tidx].i += psum[tidx - 8].i & m3;
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) == 16
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  8
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  4
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  2
    if( lane >= 1 )      psum[tidx].i += psum[tidx - 1].i;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  4
    if( lane >= 2 )      psum[tidx].i += psum[tidx - 2].i;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  8
    if( lane >= 4 )      psum[tidx].i += psum[tidx - 4].i;
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) == 16
    if( lane >= 8 )      psum[tidx].i += psum[tidx - 8].i;
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) == 16
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  8
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  4
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#endif//if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >= 32
    //---------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC_SORT
    //---------------------------------------------------------------------
  }/* if( (tidx - lane) < Ndat ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 2. prefix sum about tail of each warp */
  //-----------------------------------------------------------------------
  int scan = 0;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  scan = val;
#else///USE_WARP_SHUFFLE_FUNC_SORT
  if( tidx < Ndat )
    scan = psum[tidx].i;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
  __syncthreads();
  /* warpSize = 32 = 2^5 */
#   if  ((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >= 32
  if( tidx < (Ndat >> 5) ){
    //---------------------------------------------------------------------
    val = psum[tidx * warpSize + warpSize - 1].i;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  2
    const int groupSize = Ndat >> 5;
    tmp = __shfl_up(val,  1, groupSize) & m0;    val += tmp;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  4
    tmp = __shfl_up(val,  2, groupSize) & m1;    val += tmp;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  8
    tmp = __shfl_up(val,  4, groupSize) & m2;    val += tmp;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >= 16
    tmp = __shfl_up(val,  8, groupSize) & m3;    val += tmp;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) == 32
    tmp = __shfl_up(val, 16, groupSize) & m4;    val += tmp;
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) == 32
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >= 16
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  8
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  4
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  2
    const int groupSize = Ndat >> 5;
    tmp = __shfl_up(val,  1, groupSize);    if( lane >=  1 )      val += tmp;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  4
    tmp = __shfl_up(val,  2, groupSize);    if( lane >=  2 )      val += tmp;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  8
    tmp = __shfl_up(val,  4, groupSize);    if( lane >=  4 )      val += tmp;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >= 16
    tmp = __shfl_up(val,  8, groupSize);    if( lane >=  8 )      val += tmp;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) == 32
    tmp = __shfl_up(val, 16, groupSize);    if( lane >= 16 )      val += tmp;
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) == 32
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >= 16
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  8
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  4
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#else///USE_WARP_SHUFFLE_FUNC_SORT
    psum[tidx].i = val;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  2
    psum[tidx].i += psum[tidx -  1].i & m0;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  4
    psum[tidx].i += psum[tidx -  2].i & m1;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  8
    psum[tidx].i += psum[tidx -  4].i & m2;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >= 16
    psum[tidx].i += psum[tidx -  8].i & m3;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) == 32
    psum[tidx].i += psum[tidx - 16].i & m4;
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) == 32
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >= 16
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  8
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  4
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  2
    if( lane >=  1 )      psum[tidx].i += psum[tidx -  1].i;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  4
    if( lane >=  2 )      psum[tidx].i += psum[tidx -  2].i;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  8
    if( lane >=  4 )      psum[tidx].i += psum[tidx -  4].i;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >= 16
    if( lane >=  8 )      psum[tidx].i += psum[tidx -  8].i;
#   if  (((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) == 32
    if( lane >= 16 )      psum[tidx].i += psum[tidx - 16].i;
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) == 32
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >= 16
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  8
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  4
#endif//(((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >> 5) >=  2
#endif//USE_MASK_INSTEAD_OF_IF
    val = psum[tidx].i;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
    //---------------------------------------------------------------------
    psum[tidx].i = val;
    //---------------------------------------------------------------------
  }/* if( tidx < (Ndat >> 5) ){ */
  __syncthreads();
#endif//((1 << (RADIX_SORT_CHECK_BITS - 1)) * (NTHREADS_SORT >> 5) * RADIX_SORT_ELEMENTS_PER_THREAD) >= 32
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 3. prefix sum within a block */
  //-----------------------------------------------------------------------
  /* warpSize = 32 = 2^5 */
  /* if( tidx >= warpSize ) */
  if( (tidx >= warpSize) && (tidx < Ndat) )
    scan += psum[(tidx >> 5) - 1].i;
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 4. upload calculate prefix sum */
  //-----------------------------------------------------------------------
  if( tidx < Ndat )
    psum[tidx].i = scan;
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  RADIX_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_PREFIX_SUM_GRID
(const int tidx, const int lane, volatile uint_int * RESTRICT psum
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
)
{
  //-----------------------------------------------------------------------
  /* 1. prefix sum within warp */
  int val;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  int tmp;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
  /* load index */
  val = psum[tidx].i;
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */  tmp = __shfl_up(val,  1, warpSize) & m0;  val += tmp;
  /* #   if  warpSize >=  4 */  tmp = __shfl_up(val,  2, warpSize) & m1;  val += tmp;
  /* #   if  warpSize >=  8 */  tmp = __shfl_up(val,  4, warpSize) & m2;  val += tmp;
  /* #   if  warpSize >= 16 */  tmp = __shfl_up(val,  8, warpSize) & m3;  val += tmp;
  /* #   if  warpSize >= 32 */  tmp = __shfl_up(val, 16, warpSize) & m4;  val += tmp;
#else///USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */  tmp = __shfl_up(val,  1, warpSize);  if( lane >=  1 )    val += tmp;
  /* #   if  warpSize >=  4 */  tmp = __shfl_up(val,  2, warpSize);  if( lane >=  2 )    val += tmp;
  /* #   if  warpSize >=  8 */  tmp = __shfl_up(val,  4, warpSize);  if( lane >=  4 )    val += tmp;
  /* #   if  warpSize >= 16 */  tmp = __shfl_up(val,  8, warpSize);  if( lane >=  8 )    val += tmp;
  /* #   if  warpSize >= 32 */  tmp = __shfl_up(val, 16, warpSize);  if( lane >= 16 )    val += tmp;
#endif//USE_MASK_INSTEAD_OF_IF
  /* return calculated inclusive prefix sum */
  psum[tidx].i = val;
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */  psum[tidx].i += psum[tidx -  1].i & m0;
  /* #   if  warpSize >=  4 */  psum[tidx].i += psum[tidx -  2].i & m1;
  /* #   if  warpSize >=  8 */  psum[tidx].i += psum[tidx -  4].i & m2;
  /* #   if  warpSize >= 16 */  psum[tidx].i += psum[tidx -  8].i & m3;
  /* #   if  warpSize >= 32 */  psum[tidx].i += psum[tidx - 16].i & m4;
#else///USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */  if( lane >=  1 )    psum[tidx].i += psum[tidx -  1].i;
  /* #   if  warpSize >=  4 */  if( lane >=  2 )    psum[tidx].i += psum[tidx -  2].i;
  /* #   if  warpSize >=  8 */  if( lane >=  4 )    psum[tidx].i += psum[tidx -  4].i;
  /* #   if  warpSize >= 16 */  if( lane >=  8 )    psum[tidx].i += psum[tidx -  8].i;
  /* #   if  warpSize >= 32 */  if( lane >= 16 )    psum[tidx].i += psum[tidx - 16].i;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 2. prefix sum about tail of each warp */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  int scan = val;
#else///USE_WARP_SHUFFLE_FUNC_SORT
  int scan = psum[tidx].i;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
  __syncthreads();
  /* warpSize = 32 = 2^5 */
  if( tidx < (NTHREADS_SORT >> 5) ){
    //---------------------------------------------------------------------
    val = psum[tidx * warpSize + warpSize - 1].i;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    const int groupSize = NTHREADS_SORT >> 5;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  (NTHREADS_SORT >> 5) >=  2
    tmp = __shfl_up(val,  1, groupSize) & m0;    val += tmp;
#   if  (NTHREADS_SORT >> 5) >=  4
    tmp = __shfl_up(val,  2, groupSize) & m1;    val += tmp;
#   if  (NTHREADS_SORT >> 5) >=  8
    tmp = __shfl_up(val,  4, groupSize) & m2;    val += tmp;
#   if  (NTHREADS_SORT >> 5) >= 16
    tmp = __shfl_up(val,  8, groupSize) & m3;    val += tmp;
#   if  (NTHREADS_SORT >> 5) == 32
    tmp = __shfl_up(val, 16, groupSize) & m4;    val += tmp;
#endif//(NTHREADS_SORT >> 5) == 32
#endif//(NTHREADS_SORT >> 5) >= 16
#endif//(NTHREADS_SORT >> 5) >=  8
#endif//(NTHREADS_SORT >> 5) >=  4
#endif//(NTHREADS_SORT >> 5) >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  (NTHREADS_SORT >> 5) >=  2
    tmp = __shfl_up(val,  1, groupSize);    if( lane >=  1 )      val += tmp;
#   if  (NTHREADS_SORT >> 5) >=  4
    tmp = __shfl_up(val,  2, groupSize);    if( lane >=  2 )      val += tmp;
#   if  (NTHREADS_SORT >> 5) >=  8
    tmp = __shfl_up(val,  4, groupSize);    if( lane >=  4 )      val += tmp;
#   if  (NTHREADS_SORT >> 5) >= 16
    tmp = __shfl_up(val,  8, groupSize);    if( lane >=  8 )      val += tmp;
#   if  (NTHREADS_SORT >> 5) == 32
    tmp = __shfl_up(val, 16, groupSize);    if( lane >= 16 )      val += tmp;
#endif//(NTHREADS_SORT >> 5) == 32
#endif//(NTHREADS_SORT >> 5) >= 16
#endif//(NTHREADS_SORT >> 5) >=  8
#endif//(NTHREADS_SORT >> 5) >=  4
#endif//(NTHREADS_SORT >> 5) >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#else///USE_WARP_SHUFFLE_FUNC_SORT
    psum[tidx].i = val;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  (NTHREADS_SORT >> 5) >=  2
    psum[tidx].i += psum[tidx -  1].i & m0;
#   if  (NTHREADS_SORT >> 5) >=  4
    psum[tidx].i += psum[tidx -  2].i & m1;
#   if  (NTHREADS_SORT >> 5) >=  8
    psum[tidx].i += psum[tidx -  4].i & m2;
#   if  (NTHREADS_SORT >> 5) >= 16
    psum[tidx].i += psum[tidx -  8].i & m3;
#   if  (NTHREADS_SORT >> 5) == 32
    psum[tidx].i += psum[tidx - 16].i & m4;
#endif//(NTHREADS_SORT >> 5) == 32
#endif//(NTHREADS_SORT >> 5) >= 16
#endif//(NTHREADS_SORT >> 5) >=  8
#endif//(NTHREADS_SORT >> 5) >=  4
#endif//(NTHREADS_SORT >> 5) >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  (NTHREADS_SORT >> 5) >=  2
    if( lane >=  1 )      psum[tidx].i += psum[tidx -  1].i;
#   if  (NTHREADS_SORT >> 5) >=  4
    if( lane >=  2 )      psum[tidx].i += psum[tidx -  2].i;
#   if  (NTHREADS_SORT >> 5) >=  8
    if( lane >=  4 )      psum[tidx].i += psum[tidx -  4].i;
#   if  (NTHREADS_SORT >> 5) >= 16
    if( lane >=  8 )      psum[tidx].i += psum[tidx -  8].i;
#   if  (NTHREADS_SORT >> 5) == 32
    if( lane >= 16 )      psum[tidx].i += psum[tidx - 16].i;
#endif//(NTHREADS_SORT >> 5) == 32
#endif//(NTHREADS_SORT >> 5) >= 16
#endif//(NTHREADS_SORT >> 5) >=  8
#endif//(NTHREADS_SORT >> 5) >=  4
#endif//(NTHREADS_SORT >> 5) >=  2
#endif//USE_MASK_INSTEAD_OF_IF
    val = psum[tidx].i;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
    //---------------------------------------------------------------------
    psum[tidx].i = val;
    //---------------------------------------------------------------------
  }/* if( tidx < (NTHREADS_SORT >> 5) ){ */
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 3. prefix sum within a block */
  //-----------------------------------------------------------------------
  /* warpSize = 32 = 2^5 */
  if( tidx >= warpSize )
    scan += psum[(tidx >> 5) - 1].i;
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 4. upload calculate prefix sum */
  //-----------------------------------------------------------------------
  psum[tidx].i = scan;
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_PREFIX_SUM_GRID_ACCUMULATION
(const int tidx, const int lane, volatile uint_int * RESTRICT psum
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* 1. prefix sum within warp */
  int val;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  int tmp;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
  /* load index */
  val = psum[tidx].i;
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */  tmp = __shfl_up(val,  1, warpSize) & m0;  val += tmp;
  /* #   if  warpSize >=  4 */  tmp = __shfl_up(val,  2, warpSize) & m1;  val += tmp;
  /* #   if  warpSize >=  8 */  tmp = __shfl_up(val,  4, warpSize) & m2;  val += tmp;
  /* #   if  warpSize >= 16 */  tmp = __shfl_up(val,  8, warpSize) & m3;  val += tmp;
  /* #   if  warpSize >= 32 */  tmp = __shfl_up(val, 16, warpSize) & m4;  val += tmp;
#else///USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */  tmp = __shfl_up(val,  1, warpSize);  if( lane >=  1 )    val += tmp;
  /* #   if  warpSize >=  4 */  tmp = __shfl_up(val,  2, warpSize);  if( lane >=  2 )    val += tmp;
  /* #   if  warpSize >=  8 */  tmp = __shfl_up(val,  4, warpSize);  if( lane >=  4 )    val += tmp;
  /* #   if  warpSize >= 16 */  tmp = __shfl_up(val,  8, warpSize);  if( lane >=  8 )    val += tmp;
  /* #   if  warpSize >= 32 */  tmp = __shfl_up(val, 16, warpSize);  if( lane >= 16 )    val += tmp;
#endif//USE_MASK_INSTEAD_OF_IF
  /* return calculated inclusive prefix sum */
  psum[tidx].i = val;
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */  psum[tidx].i += psum[tidx -  1].i & m0;
  /* #   if  warpSize >=  4 */  psum[tidx].i += psum[tidx -  2].i & m1;
  /* #   if  warpSize >=  8 */  psum[tidx].i += psum[tidx -  4].i & m2;
  /* #   if  warpSize >= 16 */  psum[tidx].i += psum[tidx -  8].i & m3;
  /* #   if  warpSize >= 32 */  psum[tidx].i += psum[tidx - 16].i & m4;
#else///USE_MASK_INSTEAD_OF_IF
  /* #   if  warpSize >=  2 */  if( lane >=  1 )    psum[tidx].i += psum[tidx -  1].i;
  /* #   if  warpSize >=  4 */  if( lane >=  2 )    psum[tidx].i += psum[tidx -  2].i;
  /* #   if  warpSize >=  8 */  if( lane >=  4 )    psum[tidx].i += psum[tidx -  4].i;
  /* #   if  warpSize >= 16 */  if( lane >=  8 )    psum[tidx].i += psum[tidx -  8].i;
  /* #   if  warpSize >= 32 */  if( lane >= 16 )    psum[tidx].i += psum[tidx - 16].i;
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 2. prefix sum about tail of each warp */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  int scan = val;
#else///USE_WARP_SHUFFLE_FUNC_SORT
  int scan = psum[tidx].i;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
  __syncthreads();
  /* warpSize = 32 = 2^5 */
  if( tidx < (NTHREADS_SORT_ACCUMULATION >> 5) ){
    //---------------------------------------------------------------------
    val = psum[tidx * warpSize + warpSize - 1].i;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    const int groupSize = NTHREADS_SORT_ACCUMULATION >> 5;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  2
    tmp = __shfl_up(val,  1, groupSize) & m0;    val += tmp;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  4
    tmp = __shfl_up(val,  2, groupSize) & m1;    val += tmp;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  8
    tmp = __shfl_up(val,  4, groupSize) & m2;    val += tmp;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >= 16
    tmp = __shfl_up(val,  8, groupSize) & m3;    val += tmp;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) == 32
    tmp = __shfl_up(val, 16, groupSize) & m4;    val += tmp;
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) == 32
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >= 16
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  8
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  4
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  2
    tmp = __shfl_up(val,  1, groupSize);    if( lane >=  1 )      val += tmp;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  4
    tmp = __shfl_up(val,  2, groupSize);    if( lane >=  2 )      val += tmp;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  8
    tmp = __shfl_up(val,  4, groupSize);    if( lane >=  4 )      val += tmp;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >= 16
    tmp = __shfl_up(val,  8, groupSize);    if( lane >=  8 )      val += tmp;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) == 32
    tmp = __shfl_up(val, 16, groupSize);    if( lane >= 16 )      val += tmp;
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) == 32
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >= 16
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  8
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  4
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#else///USE_WARP_SHUFFLE_FUNC_SORT
    psum[tidx].i = val;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  2
    psum[tidx].i += psum[tidx -  1].i & m0;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  4
    psum[tidx].i += psum[tidx -  2].i & m1;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  8
    psum[tidx].i += psum[tidx -  4].i & m2;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >= 16
    psum[tidx].i += psum[tidx -  8].i & m3;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) == 32
    psum[tidx].i += psum[tidx - 16].i & m4;
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) == 32
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >= 16
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  8
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  4
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  2
    if( lane >=  1 )      psum[tidx].i += psum[tidx -  1].i;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  4
    if( lane >=  2 )      psum[tidx].i += psum[tidx -  2].i;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >=  8
    if( lane >=  4 )      psum[tidx].i += psum[tidx -  4].i;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) >= 16
    if( lane >=  8 )      psum[tidx].i += psum[tidx -  8].i;
#   if  (NTHREADS_SORT_ACCUMULATION >> 5) == 32
    if( lane >= 16 )      psum[tidx].i += psum[tidx - 16].i;
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) == 32
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >= 16
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  8
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  4
#endif//(NTHREADS_SORT_ACCUMULATION >> 5) >=  2
#endif//USE_MASK_INSTEAD_OF_IF
    val = psum[tidx].i;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
    //---------------------------------------------------------------------
    psum[tidx].i = val;
    //---------------------------------------------------------------------
  }/* if( tidx < (NTHREADS_SORT_ACCUMULATION >> 5) ){ */
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 3. prefix sum within a block */
  //-----------------------------------------------------------------------
  /* warpSize = 32 = 2^5 */
  if( tidx >= warpSize )
    scan += psum[(tidx >> 5) - 1].i;
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 4. upload calculate prefix sum */
  //-----------------------------------------------------------------------
  psum[tidx].i = scan;
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//RADIX_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_TSUB
(const int numIter, const int lane, RADIX_DATA_TYPE_KEY * RESTRICT key
#ifndef SORT_ONLY
 , int * RESTRICT idx
#endif//SORT_ONLY
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
 , const int tidx, volatile uint * RESTRICT smem
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#   if  RADIX_SORT_CHECK_BITS > 2
 , uint4_array * RESTRICT sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
 )
{
  //-----------------------------------------------------------------------
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
  const int    hidx = (tidx - lane) * 4 * RADIX_SORT_ELEMENTS_PER_THREAD + lane;
  const int tailIdx = (tidx - lane) * 4 * RADIX_SORT_ELEMENTS_PER_THREAD + TSUB_SORT - 1;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
  const int m0 = (lane >=  1) ? (0xffffffff) : 0;
#   if  TSUB_SORT >=  4
  const int m1 = (lane >=  2) ? (0xffffffff) : 0;
#   if  TSUB_SORT >=  8
  const int m2 = (lane >=  4) ? (0xffffffff) : 0;
#   if  TSUB_SORT >= 16
  const int m3 = (lane >=  8) ? (0xffffffff) : 0;
#   if  TSUB_SORT == 32
  const int m4 = (lane >= 16) ? (0xffffffff) : 0;
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------
  /* LSD radix sort */
  for(int iter = 0; iter < numIter; iter++){
    //---------------------------------------------------------------------
    /* load the inputted key */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    const RADIX_SORT_UINT key_src = key[lane];
#ifndef SORT_ONLY
    const RADIX_SORT_INT  idx_src = idx[lane];
#endif//SORT_ONLY
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    const RADIX_SORT_UINT key_src =
      {
	key[lane                ],
	key[lane +     TSUB_SORT]
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	key[lane + 2 * TSUB_SORT],
	key[lane + 3 * TSUB_SORT]
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#ifndef SORT_ONLY
    const RADIX_SORT_INT  idx_src =
      {
	idx[lane                ],
	idx[lane +     TSUB_SORT]
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	idx[lane + 2 * TSUB_SORT],
	idx[lane + 3 * TSUB_SORT]
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
    //---------------------------------------------------------------------
    /* pick up RADIX_SORT_CHECK_BITS bits of key */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    const RADIX_SORT_INT subkey = (key_src >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    const RADIX_SORT_INT subkey =
      {
	(key_src.x >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK,
	(key_src.y >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	(key_src.z >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK,
	(key_src.w >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
    //---------------------------------------------------------------------
    /* setup for prefix sum within warp */
#   if  RADIX_SORT_CHECK_BITS <= 2
    RADIX_SORT_UINT psum;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    RADIX_SORT_UINT scan = 1 << ((subkey & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    RADIX_SORT_UINT scan =
      {
	1 << ((subkey.x & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS),
	1 << ((subkey.y & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS)
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	1 << ((subkey.z & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS),
	1 << ((subkey.w & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS)
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    uint4_array psum0;
    uint4_array scan0 = {0, 0, 0, 0};
    sbuf[lane] = scan0;
    sbuf[lane].a[subkey   >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey   & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan0 = sbuf[lane];
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    uint4_array psum0, psum1;
    uint4_array scan0 = {0, 0, 0, 0};
    sbuf[lane] = scan0;
    sbuf[lane].a[subkey.x >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.x & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan0 = sbuf[lane];
    uint4_array scan1 = {0, 0, 0, 0};
    sbuf[lane] = scan1;
    sbuf[lane].a[subkey.y >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.y & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan1 = sbuf[lane];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    uint4_array psum2, psum3;
    uint4_array scan2 = {0, 0, 0, 0};
    sbuf[lane] = scan2;
    sbuf[lane].a[subkey.z >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.z & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan2 = sbuf[lane];
    uint4_array scan3 = {0, 0, 0, 0};
    sbuf[lane] = scan3;
    sbuf[lane].a[subkey.w >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.w & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan3 = sbuf[lane];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* get inclusive prefix sum */
    RADIX_PREFIX_SUM_TSUB(lane
#   if  RADIX_SORT_CHECK_BITS <= 2
			  , scan, &psum
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
			  , scan0.m, &psum0.m
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
			  , scan0.m, &psum0.m, scan1.m, &psum1.m
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
			  , scan2.m, &psum2.m, scan3.m, &psum3.m
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
			  , hidx, smem
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  TSUB_SORT >=  2
			  , m0
#   if  TSUB_SORT >=  4
			  , m1
#   if  TSUB_SORT >=  8
			  , m2
#   if  TSUB_SORT >= 16
			  , m3
#   if  TSUB_SORT == 32
			  , m4
#endif//TSUB_SORT == 32
#endif//TSUB_SORT >= 16
#endif//TSUB_SORT >=  8
#endif//TSUB_SORT >=  4
#endif//TSUB_SORT >=  2
#endif//USE_MASK_INSTEAD_OF_IF
			  );
    //---------------------------------------------------------------------
    /* get number of elements in each bucket using inclusive prefix sum */
#   if  RADIX_SORT_CHECK_BITS <= 2
    RADIX_SORT_UINT_INT tmp;    tmp.u = psum;
#else///RADIX_SORT_CHECK_BITS <= 2
    uint4_int4 tmp;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    tmp.u = psum0.m;
    tmp.i.x = __shfl(tmp.i.x, TSUB_SORT - 1, TSUB_SORT);
    tmp.i.y = __shfl(tmp.i.y, TSUB_SORT - 1, TSUB_SORT);
    tmp.i.z = __shfl(tmp.i.z, TSUB_SORT - 1, TSUB_SORT);
    tmp.i.w = __shfl(tmp.i.w, TSUB_SORT - 1, TSUB_SORT);
#else///USE_WARP_SHUFFLE_FUNC_SORT
    if( lane == (TSUB_SORT - 1) ){
      smem[tailIdx    ] = psum0.m.x;
      smem[tailIdx + 1] = psum0.m.y;
      smem[tailIdx + 2] = psum0.m.z;
      smem[tailIdx + 3] = psum0.m.w;
    }/* if( lane == (TSUB_SORT - 1) ){ */
    tmp.u.x = smem[tailIdx    ];
    tmp.u.y = smem[tailIdx + 1];
    tmp.u.z = smem[tailIdx + 2];
    tmp.u.w = smem[tailIdx + 3];
#endif//USE_WARP_SHUFFLE_FUNC_SORT
    psum1.m.x += tmp.u.x;    psum1.m.y += tmp.u.y;    psum1.m.z += tmp.u.z;    psum1.m.w += tmp.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    tmp.u = psum1.m;
    tmp.i.x = __shfl(tmp.i.x, TSUB_SORT - 1, TSUB_SORT);
    tmp.i.y = __shfl(tmp.i.y, TSUB_SORT - 1, TSUB_SORT);
    tmp.i.z = __shfl(tmp.i.z, TSUB_SORT - 1, TSUB_SORT);
    tmp.i.w = __shfl(tmp.i.w, TSUB_SORT - 1, TSUB_SORT);
#else///USE_WARP_SHUFFLE_FUNC_SORT
    if( lane == (TSUB_SORT - 1) ){
      smem[tailIdx    ] = psum1.m.x;
      smem[tailIdx + 1] = psum1.m.y;
      smem[tailIdx + 2] = psum1.m.z;
      smem[tailIdx + 3] = psum1.m.w;
    }/* if( lane == (TSUB_SORT - 1) ){ */
    tmp.u.x = smem[tailIdx    ];
    tmp.u.y = smem[tailIdx + 1];
    tmp.u.z = smem[tailIdx + 2];
    tmp.u.w = smem[tailIdx + 3];
#endif//USE_WARP_SHUFFLE_FUNC_SORT
    psum2.m.x += tmp.u.x;    psum2.m.y += tmp.u.y;    psum2.m.z += tmp.u.z;    psum2.m.w += tmp.u.w;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    tmp.u = psum2.m;
    tmp.i.x = __shfl(tmp.i.x, TSUB_SORT - 1, TSUB_SORT);
    tmp.i.y = __shfl(tmp.i.y, TSUB_SORT - 1, TSUB_SORT);
    tmp.i.z = __shfl(tmp.i.z, TSUB_SORT - 1, TSUB_SORT);
    tmp.i.w = __shfl(tmp.i.w, TSUB_SORT - 1, TSUB_SORT);
#else///USE_WARP_SHUFFLE_FUNC_SORT
    if( lane == (TSUB_SORT - 1) ){
      smem[tailIdx    ] = psum2.m.x;
      smem[tailIdx + 1] = psum2.m.y;
      smem[tailIdx + 2] = psum2.m.z;
      smem[tailIdx + 3] = psum2.m.w;
    }/* if( lane == (TSUB_SORT - 1) ){ */
    tmp.u.x = smem[tailIdx    ];
    tmp.u.y = smem[tailIdx + 1];
    tmp.u.z = smem[tailIdx + 2];
    tmp.u.w = smem[tailIdx + 3];
#endif//USE_WARP_SHUFFLE_FUNC_SORT
    psum3.m.x += tmp.u.x;    psum3.m.y += tmp.u.y;    psum3.m.z += tmp.u.z;    psum3.m.w += tmp.u.w;
    tmp.u = psum3.m;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 4
    tmp.u = psum1.m;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#else///RADIX_SORT_ELEMENTS_PER_THREAD >= 2
    tmp.u = psum0.m;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
#endif//RADIX_SORT_CHECK_BITS <= 2
    //---------------------------------------------------------------------
    /* get number of elements contained in each backet */
#   if  RADIX_SORT_CHECK_BITS <= 2
    RADIX_SORT_UINT_INT tail;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    tail.i = __shfl(tmp.i, TSUB_SORT - 1, TSUB_SORT);
#else///USE_WARP_SHUFFLE_FUNC_SORT
    if( lane == (TSUB_SORT - 1) )      smem[tailIdx] = tmp.u;
    tail.u = smem[tailIdx];
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    tail.i.x = __shfl(tmp.i.x, TSUB_SORT - 1, TSUB_SORT);    tail.i.y = __shfl(tmp.i.y, TSUB_SORT - 1, TSUB_SORT);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    tail.i.z = __shfl(tmp.i.z, TSUB_SORT - 1, TSUB_SORT);    tail.i.w = __shfl(tmp.i.w, TSUB_SORT - 1, TSUB_SORT);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#else///USE_WARP_SHUFFLE_FUNC_SORT
    if( lane == (TSUB_SORT - 1) ){
      smem[tailIdx    ] = tmp.u.x;      smem[tailIdx + 1] = tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
      smem[tailIdx + 2] = tmp.u.z;      smem[tailIdx + 3] = tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
    }/* if( lane == (TSUB_SORT - 1) ){ */
    tail.u.x = smem[tailIdx    ];    tail.u.y = smem[tailIdx + 1];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    tail.u.z = smem[tailIdx + 2];    tail.u.w = smem[tailIdx + 3];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
    uint4_int4 tail;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    tail.i.x = __shfl(tmp.i.x, TSUB_SORT - 1, TSUB_SORT);
    tail.i.y = __shfl(tmp.i.y, TSUB_SORT - 1, TSUB_SORT);
    tail.i.z = __shfl(tmp.i.z, TSUB_SORT - 1, TSUB_SORT);
    tail.i.w = __shfl(tmp.i.w, TSUB_SORT - 1, TSUB_SORT);
#else///USE_WARP_SHUFFLE_FUNC_SORT
    if( lane == (TSUB_SORT - 1) ){
      smem[tailIdx    ] = tmp.u.x;
      smem[tailIdx + 1] = tmp.u.y;
      smem[tailIdx + 2] = tmp.u.z;
      smem[tailIdx + 3] = tmp.u.w;
    }/* if( lane == (TSUB_SORT - 1) ){ */
    tail.u.x = smem[tailIdx    ];
    tail.u.y = smem[tailIdx + 1];
    tail.u.z = smem[tailIdx + 2];
    tail.u.w = smem[tailIdx + 3];
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#endif//RADIX_SORT_CHECK_BITS <= 2
    //---------------------------------------------------------------------
    /* get prefix sum about RADIX_SORT_NUM_BUCKET buckets */
    /*  8bits shift toward left corresponds to multiplying 2^8  =      256 = 0x0000100; */
    /* 16bits shift toward left corresponds to multiplying 2^16 =    65536 = 0x0010000; */
    /* 24bits shift toward left corresponds to multiplying 2^24 = 16777216 = 0x1000000; */
    /*           65536 + 16777216 = 16842752 = 0x1010000 */
    /*     256 + 65536 + 16777216 = 16843008 = 0x1010100 */
    /* 1 + 256 + 65536 + 16777216 = 16843009 = 0x1010101 */
#   if  RADIX_SORT_CHECK_BITS <= 2
    /* get prefix sum about 4 buckets */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
#ifdef  INCREASE_INSTRUCTION_LEVEL_PARALLELISM
    const int3 cnts = { tail.u                                 & RADIX_WARP_SHIFT_MASK,
		       (tail.u >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
		       (tail.u >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    tmp.u = 0x1010100 * cnts.x + 0x1010000 * cnts.y + (cnts.z << 24);
#else///INCREASE_INSTRUCTION_LEVEL_PARALLELISM
    int cnts;
    tmp.u = 0;
    cnts  =  tail.u                                 & RADIX_WARP_SHIFT_MASK;    tmp.u += (cnts <<  8);
    cnts += (tail.u >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u += (cnts << 16);
    cnts += (tail.u >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u += (cnts << 24);
#endif//INCREASE_INSTRUCTION_LEVEL_PARALLELISM
    /* get inclusive prefix sum as a complete set */
    psum += tmp.u;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
#ifdef  INCREASE_INSTRUCTION_LEVEL_PARALLELISM
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    const int4 cnts0 = { tail.u.x                                 & RADIX_WARP_SHIFT_MASK,
			(tail.u.x >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
			(tail.u.x >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK,
			(tail.u.x >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    const int4 cnts1 = { tail.u.y                                 & RADIX_WARP_SHIFT_MASK,
			(tail.u.y >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
			(tail.u.y >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK,
			(tail.u.y >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    const int4 cnts2 = { tail.u.z                                 & RADIX_WARP_SHIFT_MASK,
			(tail.u.z >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
			(tail.u.z >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK,
			(tail.u.z >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    const int3 cnts3 = { tail.u.w                                 & RADIX_WARP_SHIFT_MASK,
			(tail.u.w >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
			(tail.u.w >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    tmp.u.x =                                             0x1010100 * (cnts0.x + cnts1.x + cnts2.x + cnts3.x) + 0x1010000 * (cnts0.y + cnts1.y + cnts2.y + cnts3.y) + ((cnts0.z + cnts1.z + cnts2.z + cnts3.z) << 24);
    tmp.u.y = 0x1010101 *  cnts0.x                      + 0x1010100 * (cnts1.x + cnts2.x + cnts3.x + cnts0.y) + 0x1010000 * (cnts1.y + cnts2.y + cnts3.y + cnts0.z) + ((cnts1.z + cnts2.z + cnts3.z + cnts0.w) << 24);
    tmp.u.z = 0x1010101 * (cnts0.x + cnts1.x          ) + 0x1010100 * (cnts2.x + cnts3.x + cnts0.y + cnts1.y) + 0x1010000 * (cnts2.y + cnts3.y + cnts0.z + cnts1.z) + ((cnts2.z + cnts3.z + cnts0.w + cnts1.w) << 24);
    tmp.u.w = 0x1010101 * (cnts0.x + cnts1.x + cnts2.x) + 0x1010100 * (cnts3.x + cnts0.y + cnts1.y + cnts2.y) + 0x1010000 * (cnts3.y + cnts0.z + cnts1.z + cnts2.z) + ((cnts3.z + cnts0.w + cnts1.w + cnts2.w) << 24);
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 4
    const int4 cnts0 = { tail.u.x                                 & RADIX_WARP_SHIFT_MASK,
			(tail.u.x >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
			(tail.u.x >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK,
			(tail.u.x >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    const int3 cnts1 = { tail.u.y                                 & RADIX_WARP_SHIFT_MASK,
			(tail.u.y >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
			(tail.u.y >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    tmp.u.x =                       0x1010100 * (cnts0.x + cnts1.x) + 0x1010000 * (cnts0.y + cnts1.y) + ((cnts0.z + cnts1.z) << 24);
    tmp.u.y = 0x1010101 * cnts0.x + 0x1010100 * (cnts1.x + cnts0.y) + 0x1010000 * (cnts1.y + cnts0.z) + ((cnts1.z + cnts0.w) << 24);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#else///INCREASE_INSTRUCTION_LEVEL_PARALLELISM
    tmp.u.x = tmp.u.y = 0;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    tmp.u.z = tmp.u.w = 0;
    int cnts;
    cnts  =  tail.u.x                                 & RADIX_WARP_SHIFT_MASK;    tmp.u.y +=  cnts;
    cnts +=  tail.u.y                                 & RADIX_WARP_SHIFT_MASK;    tmp.u.z +=  cnts;
    cnts +=  tail.u.z                                 & RADIX_WARP_SHIFT_MASK;    tmp.u.w +=  cnts;
    cnts +=  tail.u.w                                 & RADIX_WARP_SHIFT_MASK;    tmp.u.x += (cnts <<  8);
    cnts += (tail.u.x >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u.y += (cnts <<  8);
    cnts += (tail.u.y >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u.z += (cnts <<  8);
    cnts += (tail.u.z >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u.w += (cnts <<  8);
    cnts += (tail.u.w >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u.x += (cnts << 16);
    cnts += (tail.u.x >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.y += (cnts << 16);
    cnts += (tail.u.y >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.z += (cnts << 16);
    cnts += (tail.u.z >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.w += (cnts << 16);
    cnts += (tail.u.w >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.x += (cnts << 24);
    cnts += (tail.u.x >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.y += (cnts << 24);
    cnts += (tail.u.y >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.z += (cnts << 24);
    cnts += (tail.u.z >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.w += (cnts << 24);
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 4
    cnts  =  tail.u.x                                 & RADIX_WARP_SHIFT_MASK;    tmp.u.y +=  cnts;
    cnts +=  tail.u.y                                 & RADIX_WARP_SHIFT_MASK;    tmp.u.x += (cnts <<  8);
    cnts += (tail.u.x >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u.y += (cnts <<  8);
    cnts += (tail.u.y >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u.x += (cnts << 16);
    cnts += (tail.u.x >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.y += (cnts << 16);
    cnts += (tail.u.y >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.x += (cnts << 24);
    cnts += (tail.u.x >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.y += (cnts << 24);
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//INCREASE_INSTRUCTION_LEVEL_PARALLELISM
    /* get inclusive prefix sum as a complete set */
    psum.x += tmp.u.x;    psum.y += tmp.u.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    psum.z += tmp.u.z;    psum.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
    /* get prefix sum about 16 buckets */
#ifdef  INCREASE_INSTRUCTION_LEVEL_PARALLELISM
    const int4 cnts0 = { tail.u.x                                 & RADIX_WARP_SHIFT_MASK,
			(tail.u.x >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
			(tail.u.x >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK,
			(tail.u.x >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    const int4 cnts1 = { tail.u.y                                 & RADIX_WARP_SHIFT_MASK,
			(tail.u.y >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
			(tail.u.y >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK,
			(tail.u.y >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    const int4 cnts2 = { tail.u.z                                 & RADIX_WARP_SHIFT_MASK,
			(tail.u.z >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
			(tail.u.z >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK,
			(tail.u.z >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    const int3 cnts3 = { tail.u.w                                 & RADIX_WARP_SHIFT_MASK,
			(tail.u.w >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK,
			(tail.u.w >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK};
    tmp.u.x =                                                                                                                                       0x1010100 * cnts0.x + 0x1010000 * cnts0.y + (cnts0.z << 24);
    tmp.u.y = 0x1010101 * (cnts0.x + cnts0.y + cnts0.z + cnts0.w                                                                                ) + 0x1010100 * cnts1.x + 0x1010000 * cnts1.y + (cnts1.z << 24);
    tmp.u.z = 0x1010101 * (cnts0.x + cnts0.y + cnts0.z + cnts0.w + cnts1.x + cnts1.y + cnts1.z + cnts1.w                                        ) + 0x1010100 * cnts2.x + 0x1010000 * cnts2.y + (cnts2.z << 24);
    tmp.u.w = 0x1010101 * (cnts0.x + cnts0.y + cnts0.z + cnts0.w + cnts1.x + cnts1.y + cnts1.z + cnts1.w + cnts2.x + cnts2.y + cnts2.z + cnts2.w) + 0x1010100 * cnts3.x + 0x1010000 * cnts3.y + (cnts3.z << 24);
#else///INCREASE_INSTRUCTION_LEVEL_PARALLELISM
    tmp.u.x = tmp.u.y = tmp.u.z = tmp.u.w = 0;
    int cnts;
    cnts  =  tail.u.x                                 & RADIX_WARP_SHIFT_MASK;    tmp.u.x += (cnts <<  8);
    cnts += (tail.u.x >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u.x += (cnts << 16);
    cnts += (tail.u.x >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.x += (cnts << 24);
    cnts += (tail.u.x >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.y +=  cnts;
    cnts +=  tail.u.y                                 & RADIX_WARP_SHIFT_MASK;    tmp.u.y += (cnts <<  8);
    cnts += (tail.u.y >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u.y += (cnts << 16);
    cnts += (tail.u.y >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.y += (cnts << 24);
    cnts += (tail.u.y >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.z +=  cnts;
    cnts +=  tail.u.z                                 & RADIX_WARP_SHIFT_MASK;    tmp.u.z += (cnts <<  8);
    cnts += (tail.u.z >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u.z += (cnts << 16);
    cnts += (tail.u.z >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.z += (cnts << 24);
    cnts += (tail.u.z >> (3 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.w +=  cnts;
    cnts +=  tail.u.w                                 & RADIX_WARP_SHIFT_MASK;    tmp.u.w += (cnts <<  8);
    cnts += (tail.u.w >>      RADIX_WARP_SHIFT_BITS)  & RADIX_WARP_SHIFT_MASK;    tmp.u.w += (cnts << 16);
    cnts += (tail.u.w >> (2 * RADIX_WARP_SHIFT_BITS)) & RADIX_WARP_SHIFT_MASK;    tmp.u.w += (cnts << 24);
#endif//INCREASE_INSTRUCTION_LEVEL_PARALLELISM
    /* get inclusive prefix sum as a complete set */
    psum0.m.x += tmp.u.x;    psum0.m.y += tmp.u.y;    psum0.m.z += tmp.u.z;    psum0.m.w += tmp.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
    psum1.m.x += tmp.u.x;    psum1.m.y += tmp.u.y;    psum1.m.z += tmp.u.z;    psum1.m.w += tmp.u.w;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    psum2.m.x += tmp.u.x;    psum2.m.y += tmp.u.y;    psum2.m.z += tmp.u.z;    psum2.m.w += tmp.u.w;
    psum3.m.x += tmp.u.x;    psum3.m.y += tmp.u.y;    psum3.m.z += tmp.u.z;    psum3.m.w += tmp.u.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
#endif//RADIX_SORT_CHECK_BITS <= 2
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* prefix sum within a warp */
    /* set exclusive prefix sum */
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    psum   -= scan;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    psum.x -= scan.x;    psum.y -= scan.y;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    psum.z -= scan.z;    psum.w -= scan.w;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum0.a[ii] -= scan0.a[ii];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
      psum1.a[ii] -= scan1.a[ii];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
      psum2.a[ii] -= scan2.a[ii];
      psum3.a[ii] -= scan3.a[ii];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
    }/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_CHECK_BITS <= 2
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* execute bucket sort */
    //---------------------------------------------------------------------
    RADIX_SORT_INT dstIdx;
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    dstIdx = RADIX_SORT_TSUB_GET_DSTIDX(subkey  , psum);
    key[dstIdx] = key_src;
#ifndef SORT_ONLY
    idx[dstIdx] = idx_src;
#endif//SORT_ONLY
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    dstIdx.x = RADIX_SORT_TSUB_GET_DSTIDX(subkey.x, psum.x);
    dstIdx.y = RADIX_SORT_TSUB_GET_DSTIDX(subkey.y, psum.y);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    dstIdx.z = RADIX_SORT_TSUB_GET_DSTIDX(subkey.z, psum.z);
    dstIdx.w = RADIX_SORT_TSUB_GET_DSTIDX(subkey.w, psum.w);
    key[dstIdx.z] = key_src.z;
    key[dstIdx.w] = key_src.w;
#ifndef SORT_ONLY
    idx[dstIdx.z] = idx_src.z;
    idx[dstIdx.w] = idx_src.w;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
    key[dstIdx.x] = key_src.x;
    key[dstIdx.y] = key_src.y;
#ifndef SORT_ONLY
    idx[dstIdx.x] = idx_src.x;
    idx[dstIdx.y] = idx_src.y;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    sbuf[lane] = psum0;    dstIdx   = RADIX_SORT_TSUB_GET_DSTIDX(subkey  , sbuf[lane].a[subkey   >> 2]);
    key[dstIdx] = key_src;
#ifndef SORT_ONLY
    idx[dstIdx] = idx_src;
#endif//SORT_ONLY
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    sbuf[lane] = psum0;    dstIdx.x = RADIX_SORT_TSUB_GET_DSTIDX(subkey.x, sbuf[lane].a[subkey.x >> 2]);
    sbuf[lane] = psum1;    dstIdx.y = RADIX_SORT_TSUB_GET_DSTIDX(subkey.y, sbuf[lane].a[subkey.y >> 2]);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    sbuf[lane] = psum2;    dstIdx.z = RADIX_SORT_TSUB_GET_DSTIDX(subkey.z, sbuf[lane].a[subkey.z >> 2]);
    sbuf[lane] = psum3;    dstIdx.w = RADIX_SORT_TSUB_GET_DSTIDX(subkey.w, sbuf[lane].a[subkey.w >> 2]);
    key[dstIdx.z] = key_src.z;
    key[dstIdx.w] = key_src.w;
#ifndef SORT_ONLY
    idx[dstIdx.z] = idx_src.z;
    idx[dstIdx.w] = idx_src.w;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
    key[dstIdx.x] = key_src.x;
    key[dstIdx.y] = key_src.y;
#ifndef SORT_ONLY
    idx[dstIdx.x] = idx_src.x;
    idx[dstIdx.y] = idx_src.y;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
    //---------------------------------------------------------------------
  }/* for(int iter = 0; iter < RADIX_BLOCK_NUM_ITER; iter++){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_BLCK
(const int numIter, const int tidx
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
 , const int hp
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
 , volatile RADIX_DATA_TYPE_KEY * RESTRICT key
#ifndef SORT_ONLY
 , volatile int * RESTRICT idx
#endif//SORT_ONLY
 , volatile uint_int * RESTRICT smem
#   if  RADIX_SORT_CHECK_BITS > 2
 , uint4_array * RESTRICT sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* const int warpIdx = tidx          / warpSize; */
  /* const int warpNum = NTHREADS_SORT / warpSize; */
  /* warpSize = 32 = 2^5 */
  const int warpIdx = tidx          >> 5;
  const int warpNum = NTHREADS_SORT >> 5;
  const int lane = tidx & (warpSize - 1);
  //-----------------------------------------------------------------------
  const int dataNum = (1 << (RADIX_SORT_CHECK_BITS - 1));
  const int sameNum = warpNum * RADIX_SORT_ELEMENTS_PER_THREAD;
  const int dataIdx = (tidx + 1) / sameNum;
  //-----------------------------------------------------------------------
#ifdef  ENABLE_EARLY_EXIT_BLCK
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD ==  128 = 2^ 7 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD ==  256 = 2^ 8 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD ==  512 = 2^ 9 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD == 1024 = 2^10 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD == 2048 = 2^11 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD == 4096 = 2^12 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD == 8192 = 2^13 */
  /* 16bits shift toward left corresponds to multiplying 2^16 = 65536; */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD * (1 + 16bits left shift) = NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD * (1 + 65536) */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD * (65537) is the desired mask for 32 bits value */
  const int mask = NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD * 65537;
#endif//ENABLE_EARLY_EXIT_BLCK
  //-----------------------------------------------------------------------
  /* LSD radix sort */
  for(int iter = 0; iter < numIter; iter++){
    //---------------------------------------------------------------------
    /* load the inputted key */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    const RADIX_SORT_UINT key_src = key[tidx];
#ifndef SORT_ONLY
    const RADIX_SORT_INT  idx_src = idx[tidx];
#endif//SORT_ONLY
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    const RADIX_SORT_UINT key_src =
      {
	key[tidx                    ],
	key[tidx +     NTHREADS_SORT]
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	key[tidx + 2 * NTHREADS_SORT],
	key[tidx + 3 * NTHREADS_SORT]
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#ifndef SORT_ONLY
    const RADIX_SORT_INT  idx_src =
      {
	idx[tidx                    ],
	idx[tidx +     NTHREADS_SORT]
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	idx[tidx + 2 * NTHREADS_SORT],
	idx[tidx + 3 * NTHREADS_SORT]
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
    //---------------------------------------------------------------------
    /* pick up RADIX_SORT_CHECK_BITS bits of key */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    RADIX_SORT_INT subkey = (key_src >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    RADIX_SORT_INT subkey =
      {
	(key_src.x >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK,
	(key_src.y >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	(key_src.z >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK,
	(key_src.w >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
    //---------------------------------------------------------------------
    /* setup for prefix sum within warp */
#   if  RADIX_SORT_CHECK_BITS <= 2
    RADIX_SORT_UINT_ARRAY psum0;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    RADIX_SORT_UINT_ARRAY scan0 = 1 << ((subkey & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    RADIX_SORT_UINT_ARRAY scan0 =
      {
	1 << ((subkey.x & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS),
	1 << ((subkey.y & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS)
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	1 << ((subkey.z & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS),
	1 << ((subkey.w & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS)
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    uint4_array psum0;
    uint4_array scan0 = {0, 0, 0, 0};
    sbuf[tidx] = scan0;
    sbuf[tidx].a[subkey   >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey   & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan0 = sbuf[tidx];
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    uint4_array psum0, psum1;
    uint4_array scan0 = {0, 0, 0, 0};
    sbuf[tidx] = scan0;
    sbuf[tidx].a[subkey.x >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.x & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan0 = sbuf[tidx];
    uint4_array scan1 = {0, 0, 0, 0};
    sbuf[tidx] = scan1;
    sbuf[tidx].a[subkey.y >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.y & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan1 = sbuf[tidx];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    uint4_array psum2, psum3;
    uint4_array scan2 = {0, 0, 0, 0};
    sbuf[tidx] = scan2;
    sbuf[tidx].a[subkey.z >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.z & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan2 = sbuf[tidx];
    uint4_array scan3 = {0, 0, 0, 0};
    sbuf[tidx] = scan3;
    sbuf[tidx].a[subkey.w >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.w & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan3 = sbuf[tidx];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* get inclusive prefix sum */
    //---------------------------------------------------------------------
    /* prefix sum within a warp */
    RADIX_PREFIX_SUM_WARP(lane
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
			  , scan0  , &psum0
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
			  , scan0.m, &psum0.m
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
			  , scan0.m, &psum0.m
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
			  , scan0.m, &psum0.m, scan1.m, &psum1.m
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
			  , scan2.m, &psum2.m, scan3.m, &psum3.m
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
			  , tidx, smem
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#ifdef  USE_MASK_INSTEAD_OF_IF
			  , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
			  );
    //---------------------------------------------------------------------
#if 0
    if( tidx < warpSize )
      printf("%d\t%x\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", tidx, subkey,
	     psum0.m.x & 0xff, (psum0.m.x >> 8) & 0xff, (psum0.m.x >> 16) & 0xff, (psum0.m.x >> 24) & 0xff,
	     psum0.m.y & 0xff, (psum0.m.y >> 8) & 0xff, (psum0.m.y >> 16) & 0xff, (psum0.m.y >> 24) & 0xff,
	     psum0.m.z & 0xff, (psum0.m.z >> 8) & 0xff, (psum0.m.z >> 16) & 0xff, (psum0.m.z >> 24) & 0xff,
	     psum0.m.w & 0xff, (psum0.m.w >> 8) & 0xff, (psum0.m.w >> 16) & 0xff, (psum0.m.w >> 24) & 0xff);
#endif
    //---------------------------------------------------------------------
    /* setup for prefix sum within a block */
    /* NOTE: 8bits is not sufficient in this case --> extending 16bits is required */
#   if  RADIX_SORT_CHECK_BITS <= 2
    RADIX_SORT_UINT_ARRAY psum1, scan1;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    psum1     = RADIX_CONVERT_8x4to16x2UPPER(psum0);    scan1     = RADIX_CONVERT_8x4to16x2UPPER(scan0);
    psum0     = RADIX_CONVERT_8x4to16x2LOWER(psum0);    scan0     = RADIX_CONVERT_8x4to16x2LOWER(scan0);
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      psum1.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(psum0.a[ii]);      scan1.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(scan0.a[ii]);
    }/* for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){ */
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      psum0.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(psum0.a[ii]);      scan0.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(scan0.a[ii]);
    }/* for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){ */
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    uint4_array psum4, psum5, psum6, psum7;
    uint4_array scan4, scan5, scan6, scan7;
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum7.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(psum3.a[ii]);      scan7.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(scan3.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum6.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(psum3.a[ii]);      scan6.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(scan3.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum5.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(psum2.a[ii]);      scan5.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(scan2.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum4.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(psum2.a[ii]);      scan4.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(scan2.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD < 4
    uint4_array psum2, psum3;
    uint4_array scan2, scan3;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD < 4
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum3.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(psum1.a[ii]);      scan3.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(scan1.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum2.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(psum1.a[ii]);      scan2.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(scan1.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#else///RADIX_SORT_ELEMENTS_PER_THREAD >= 2
    uint4_array psum1;
    uint4_array scan1;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum1.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(psum0.a[ii]);      scan1.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(scan0.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum0.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(psum0.a[ii]);      scan0.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(scan0.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_CHECK_BITS <= 2
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* prefix sum within a block */
    //---------------------------------------------------------------------
    /* clear the idx_dst for safety */
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
    __syncthreads();
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#pragma unroll
    for(int ii = tidx; ii < dataNum * sameNum; ii += NTHREADS_SORT)
      smem[ii].u = 0;
    __syncthreads();
    /* store the partial prefix sum to the shared memory */
    if( lane == (warpSize - 1) )
#   if  RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
      for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	//-----------------------------------------------------------------
	const int headIdx = warpIdx + ii * warpNum;
	//-----------------------------------------------------------------
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
	//-----------------------------------------------------------------
	smem[headIdx                                                 ].u = psum0;
	smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD      * warpNum].u = psum1;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
	//-----------------------------------------------------------------
	smem[headIdx                                                 ].u = psum0.a[ii];
	smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD      * warpNum].u = psum1.a[ii];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < 4; ii++){ */
#else///RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
      for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32; ii++){
	//-----------------------------------------------------------------
	const int headIdx = warpIdx + ii * 2 * sameNum;
	//-----------------------------------------------------------------
	smem[headIdx                                                 ].u = psum0.a[ii];
	smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD      * warpNum].u = psum1.a[ii];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
	smem[headIdx +                                        warpNum].u = psum2.a[ii];
	smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 1) * warpNum].u = psum3.a[ii];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	smem[headIdx +                                   2  * warpNum].u = psum4.a[ii];
	smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 2) * warpNum].u = psum5.a[ii];
	smem[headIdx +                                   3  * warpNum].u = psum6.a[ii];
	smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 3) * warpNum].u = psum7.a[ii];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_CHECK_BITS <= 2
    __syncthreads();
    //---------------------------------------------------------------------
    /* prefix sum within block (16bits * 2 is merged) */
    RADIX_PREFIX_SUM_BLCK(tidx, lane, smem
#ifdef  USE_MASK_INSTEAD_OF_IF
			  , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
			  );
    //---------------------------------------------------------------------
#ifdef  ENABLE_EARLY_EXIT_BLCK
    /* get # of elements in each bucket */
    /* 16bits value * 4 or 16 elements --> 32 bits value * 2 or 8 elements --> uint2 or 2xuint4 */
    /* unfortunately, additional variables are necessary */
    /* or, check 4 or 16 threads and accumulate results <-- this would be the best solution */
    int stock;
    if( tidx < (RADIX_SORT_NUM_BUCKET >> 1) ){
      //-------------------------------------------------------------------
      /* get masked flag to detect the signal */
      //-------------------------------------------------------------------
      /* get # of elements in each bucket */
      int flag = smem[(tidx + 1) * sameNum - 1].i;
      if( tidx > 0 )
	flag -= smem[tidx * sameNum - 1].i;
      /* apply the mask */
      flag &= mask;
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* calculate the total sum of all flags */
      //-------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
      if( tidx == 0 )
	stock = smem[0].i;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) >= 2
      int temp;
      temp = __shfl_up(flag, 1, RADIX_SORT_NUM_BUCKET >> 1);      flag += temp;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) >= 4
      temp = __shfl_up(flag, 2, RADIX_SORT_NUM_BUCKET >> 1);      flag += temp;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) == 8
      temp = __shfl_up(flag, 4, RADIX_SORT_NUM_BUCKET >> 1);      flag += temp;
#endif//(RADIX_SORT_NUM_BUCKET >> 1) == 8
#endif//(RADIX_SORT_NUM_BUCKET >> 1) >= 4
#endif//(RADIX_SORT_NUM_BUCKET >> 1) >= 2
      if( tidx == 0 )
	smem[0].i = flag;
#else///USE_WARP_SHUFFLE_FUNC_SORT
      stock = smem[tidx].i;
      smem[tidx].i = flag;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) >= 2
      int temp;
      temp = smem[tidx ^ 1];      flag += temp;      smem[tidx].i = flag;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) >= 4
      temp = smem[tidx ^ 2];      flag += temp;      smem[tidx].i = flag;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) == 8
      temp = smem[tidx ^ 4];      flag += temp;      smem[tidx].i = flag;
#endif//(RADIX_SORT_NUM_BUCKET >> 1) == 8
#endif//(RADIX_SORT_NUM_BUCKET >> 1) >= 4
#endif//(RADIX_SORT_NUM_BUCKET >> 1) >= 2
#endif//USE_WARP_SHUFFLE_FUNC_SORT
      //-------------------------------------------------------------------
    }/* if( tidx < (RADIX_SORT_NUM_BUCKET >> 1) ){ */
    __syncthreads();
    const int skip = smem[0].i;
    //---------------------------------------------------------------------
    /* if all flags are equal to zero; then sorting must be processed */
    /* if one of flags is not equal to zero; then sorting can be skipped */
    //---------------------------------------------------------------------
#if 0
    if( (tidx == 0) && (skip == 1) )
      printf("early escape within a block\n");
#endif
    if( skip == 0 )
#endif//ENABLE_EARLY_EXIT_BLCK
      {
	//-----------------------------------------------------------------
#ifdef  ENABLE_EARLY_EXIT_BLCK
	__syncthreads();
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	if( tidx == 0 )
	  smem[0].i = stock;
#else///USE_WARP_SHUFFLE_FUNC_SORT
	if( tidx < (RADIX_SORT_NUM_BUCKET >> 1) )
	  smem[tidx].i = stock;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
	__syncthreads();
#endif//ENABLE_EARLY_EXIT_BLCK
	//-----------------------------------------------------------------
	/* prefix sum within block (complete set) */
	/* 16bits shift toward left corresponds to multiplying 2^16 = 65536; */
	const uint headLower = (smem[sameNum - 1].u & 0xffff) << 16;
	int val;
	if( dataIdx < dataNum ){
	  val  =                 (smem[(dataIdx + 1) * sameNum - 1].i & 0xffff) << 16;
	  val += (dataIdx > 0) ? (smem[ dataIdx      * sameNum - 1].i           >> 16) : 0;
	}/* if( dataIdx < dataNum ){ */
	__syncthreads();
	if( tidx < dataNum * sameNum )
	  smem[tidx].i += val;
	__syncthreads();
	/* get prefix sum as a complete set */
#   if  RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
	for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	  //---------------------------------------------------------------
	  const int headIdx = (warpIdx - 1) + ii * warpNum;
	  //---------------------------------------------------------------
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
	  psum0       += (headIdx >= 0) ? (smem[headIdx].u) : (headLower);	  psum1       += smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD * warpNum].u;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
	  psum0.a[ii] += (headIdx >= 0) ? (smem[headIdx].u) : (headLower);	  psum1.a[ii] += smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD * warpNum].u;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
          //---------------------------------------------------------------
	}/* for(int ii = 0; ii < 4; ii++){ */
#else///RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
	for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32; ii++){
	  //---------------------------------------------------------------
	  const int headIdx = (warpIdx - 1) + ii * 2 * sameNum;
	  //---------------------------------------------------------------
	  psum0.a[ii] += (headIdx >= 0) ? (smem[headIdx].u) : (headLower);
	  psum1.a[ii] += smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD      * warpNum].u;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
	  psum2.a[ii] += smem[headIdx +                                        warpNum].u;
	  psum3.a[ii] += smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 1) * warpNum].u;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	  psum4.a[ii] += smem[headIdx +                                   2  * warpNum].u;
	  psum5.a[ii] += smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 2) * warpNum].u;
	  psum6.a[ii] += smem[headIdx +                                   3  * warpNum].u;
	  psum7.a[ii] += smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 3) * warpNum].u;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
          //---------------------------------------------------------------
	}/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_CHECK_BITS <= 2
	//-----------------------------------------------------------------
	/* set exclusive prefix sum */
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
	psum0 -= scan0;
	psum1 -= scan1;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
#pragma unroll
	for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	  psum0.a[ii] -= scan0.a[ii];
	  psum1.a[ii] -= scan1.a[ii];
	}/* for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){ */
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
	for(int ii = 0; ii < 4; ii++){
	  psum0.a[ii] -= scan0.a[ii];
	  psum1.a[ii] -= scan1.a[ii];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
	  psum2.a[ii] -= scan2.a[ii];
	  psum3.a[ii] -= scan3.a[ii];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	  psum4.a[ii] -= scan4.a[ii];
	  psum5.a[ii] -= scan5.a[ii];
	  psum6.a[ii] -= scan6.a[ii];
	  psum7.a[ii] -= scan7.a[ii];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
	}/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_CHECK_BITS <= 2
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* execute bucket sort */
	//-----------------------------------------------------------------
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
	subkey = RADIX_SORT_BLCK_GET_DSTIDX(subkey  , psum0, psum1);
	key[subkey] = key_src;
#ifndef SORT_ONLY
	idx[subkey] = idx_src;
#endif//SORT_ONLY
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
	subkey.x = RADIX_SORT_BLCK_GET_DSTIDX(subkey.x, psum0.m.x, psum1.m.x);
	subkey.y = RADIX_SORT_BLCK_GET_DSTIDX(subkey.y, psum0.m.y, psum1.m.y);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	subkey.z = RADIX_SORT_BLCK_GET_DSTIDX(subkey.z, psum0.m.z, psum1.m.z);
	subkey.w = RADIX_SORT_BLCK_GET_DSTIDX(subkey.w, psum0.m.w, psum1.m.w);
	key[subkey.z] = key_src.z;
	key[subkey.w] = key_src.w;
#ifndef SORT_ONLY
	idx[subkey.z] = idx_src.z;
	idx[subkey.w] = idx_src.w;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
	key[subkey.x] = key_src.x;
	key[subkey.y] = key_src.y;
#ifndef SORT_ONLY
	idx[subkey.x] = idx_src.x;
	idx[subkey.y] = idx_src.y;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
#if 0
	printf("%d\t%x\t%d\n", tidx, subkey, RADIX_SORT_BLCK_GET_DSTIDX(subkey  , psum0.m, psum1.m));
#endif
	subkey = RADIX_SORT_BLCK_GET_DSTIDX(subkey  , psum0.m, psum1.m, &sbuf[tidx]);
	key[subkey] = key_src;
#ifndef SORT_ONLY
	idx[subkey] = idx_src;
#endif//SORT_ONLY
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
	subkey.x = RADIX_SORT_BLCK_GET_DSTIDX(subkey.x, psum0.m, psum1.m, &sbuf[tidx]);
	subkey.y = RADIX_SORT_BLCK_GET_DSTIDX(subkey.y, psum2.m, psum3.m, &sbuf[tidx]);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	subkey.z = RADIX_SORT_BLCK_GET_DSTIDX(subkey.z, psum4.m, psum5.m, &sbuf[tidx]);
	subkey.w = RADIX_SORT_BLCK_GET_DSTIDX(subkey.w, psum6.m, psum7.m, &sbuf[tidx]);
	key[subkey.z] = key_src.z;
	key[subkey.w] = key_src.w;
#ifndef SORT_ONLY
	idx[subkey.z] = idx_src.z;
	idx[subkey.w] = idx_src.w;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
	key[subkey.x] = key_src.x;
	key[subkey.y] = key_src.y;
#ifndef SORT_ONLY
	idx[subkey.x] = idx_src.x;
	idx[subkey.y] = idx_src.y;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
	//-----------------------------------------------------------------
      }
    __syncthreads();
    //---------------------------------------------------------------------
  }/* for(int iter = 0; iter < RADIX_BLOCK_NUM_ITER; iter++){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* numIter = 1; */
__device__ __forceinline__ void RADIX_SORT_CHIP
(const int iter, const int tidx
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
 , const int hp
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
 , RADIX_DATA_TYPE_KEY * RESTRICT key
#ifndef SORT_ONLY
 , int * RESTRICT idx
#endif//SORT_ONLY
 , volatile uint_int * RESTRICT smem
#   if  RADIX_SORT_CHECK_BITS > 2
 , uint4_array * RESTRICT sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* const int warpIdx = tidx          / warpSize; */
  /* const int warpNum = NTHREADS_SORT / warpSize; */
  /* warpSize = 32 = 2^5 */
  const int warpIdx = tidx          >> 5;
  const int warpNum = NTHREADS_SORT >> 5;
  const int lane = tidx & (warpSize - 1);
  //-----------------------------------------------------------------------
  const int dataNum = (1 << (RADIX_SORT_CHECK_BITS - 1));
  const int sameNum = warpNum * RADIX_SORT_ELEMENTS_PER_THREAD;
  const int dataIdx = (tidx + 1) / sameNum;
  //-----------------------------------------------------------------------
#ifdef  ENABLE_EARLY_EXIT_BLCK
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD ==  128 = 2^ 7 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD ==  256 = 2^ 8 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD ==  512 = 2^ 9 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD == 1024 = 2^10 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD == 2048 = 2^11 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD == 4096 = 2^12 */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD == 8192 = 2^13 */
  /* 16bits shift toward left corresponds to multiplying 2^16 = 65536; */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD * (1 + 16bits left shift) = NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD * (1 + 65536) */
  /* NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD * (65537) is the desired mask for 32 bits value */
  const int mask = NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD * 65537;
#endif//ENABLE_EARLY_EXIT_BLCK
  //-----------------------------------------------------------------------
  /* MSD radix sort */
  /* for(int iter = 0; iter < numIter; iter++){ */
  {
    //---------------------------------------------------------------------
    /* load the inputted key */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    const RADIX_SORT_UINT key_src = key[tidx];
#ifndef SORT_ONLY
    const RADIX_SORT_INT  idx_src = idx[tidx];
#endif//SORT_ONLY
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    const RADIX_SORT_UINT key_src =
      {
	key[tidx                    ],
	key[tidx +     NTHREADS_SORT]
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	key[tidx + 2 * NTHREADS_SORT],
	key[tidx + 3 * NTHREADS_SORT]
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#ifndef SORT_ONLY
    const RADIX_SORT_INT  idx_src =
      {
	idx[tidx                    ],
	idx[tidx +     NTHREADS_SORT]
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	idx[tidx + 2 * NTHREADS_SORT],
	idx[tidx + 3 * NTHREADS_SORT]
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
    //---------------------------------------------------------------------
    /* pick up RADIX_SORT_CHECK_BITS bits of key */
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    RADIX_SORT_INT subkey = (key_src >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    /* const RADIX_SORT_INT subkey = */
    RADIX_SORT_INT subkey =
      {
	(key_src.x >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK,
	(key_src.y >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	(key_src.z >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK,
	(key_src.w >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
    //---------------------------------------------------------------------
    /* setup for prefix sum within warp */
#   if  RADIX_SORT_CHECK_BITS <= 2
    RADIX_SORT_UINT_ARRAY psum0;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    RADIX_SORT_UINT_ARRAY scan0 = 1 << ((subkey & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    RADIX_SORT_UINT_ARRAY scan0 =
      {
	1 << ((subkey.x & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS),
	1 << ((subkey.y & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS)
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	,
	1 << ((subkey.z & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS),
	1 << ((subkey.w & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS)
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
      };
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    uint4_array psum0;
    uint4_array scan0 = {0, 0, 0, 0};
    sbuf[tidx] = scan0;
    sbuf[tidx].a[subkey   >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey   & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan0 = sbuf[tidx];
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
    uint4_array psum0, psum1;
    uint4_array scan0 = {0, 0, 0, 0};
    sbuf[tidx] = scan0;
    sbuf[tidx].a[subkey.x >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.x & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan0 = sbuf[tidx];
    uint4_array scan1 = {0, 0, 0, 0};
    sbuf[tidx] = scan1;
    sbuf[tidx].a[subkey.y >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.y & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan1 = sbuf[tidx];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    uint4_array psum2, psum3;
    uint4_array scan2 = {0, 0, 0, 0};
    sbuf[tidx] = scan2;
    sbuf[tidx].a[subkey.z >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.z & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan2 = sbuf[tidx];
    uint4_array scan3 = {0, 0, 0, 0};
    sbuf[tidx] = scan3;
    sbuf[tidx].a[subkey.w >> RADIX_WARP_CHECK_BITS] |= 1 << ((subkey.w & RADIX_WARP_CHECK_MASK) * RADIX_WARP_SHIFT_BITS);
    scan3 = sbuf[tidx];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* get inclusive prefix sum */
    //---------------------------------------------------------------------
    /* prefix sum within a warp */
    RADIX_PREFIX_SUM_WARP(lane
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
			  , scan0  , &psum0
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
			  , scan0.m, &psum0.m
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
			  , scan0.m, &psum0.m
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
			  , scan0.m, &psum0.m, scan1.m, &psum1.m
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
			  , scan2.m, &psum2.m, scan3.m, &psum3.m
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
			  , tidx, smem
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#ifdef  USE_MASK_INSTEAD_OF_IF
			  , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
			  );
    //---------------------------------------------------------------------
    /* setup for prefix sum within a block */
    /* NOTE: 8bits is not sufficient in this case --> extending 16bits is required */
#   if  RADIX_SORT_CHECK_BITS <= 2
    RADIX_SORT_UINT_ARRAY psum1, scan1;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
    psum1     = RADIX_CONVERT_8x4to16x2UPPER(psum0);    scan1     = RADIX_CONVERT_8x4to16x2UPPER(scan0);
    psum0     = RADIX_CONVERT_8x4to16x2LOWER(psum0);    scan0     = RADIX_CONVERT_8x4to16x2LOWER(scan0);
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      psum1.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(psum0.a[ii]);      scan1.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(scan0.a[ii]);
    }/* for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){ */
#pragma unroll
    for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
      psum0.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(psum0.a[ii]);      scan0.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(scan0.a[ii]);
    }/* for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){ */
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
    uint4_array psum4, psum5, psum6, psum7;
    uint4_array scan4, scan5, scan6, scan7;
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum7.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(psum3.a[ii]);      scan7.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(scan3.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum6.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(psum3.a[ii]);      scan6.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(scan3.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum5.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(psum2.a[ii]);      scan5.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(scan2.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum4.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(psum2.a[ii]);      scan4.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(scan2.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD < 4
    uint4_array psum2, psum3;
    uint4_array scan2, scan3;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD < 4
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum3.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(psum1.a[ii]);      scan3.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(scan1.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum2.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(psum1.a[ii]);      scan2.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(scan1.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#else///RADIX_SORT_ELEMENTS_PER_THREAD >= 2
    uint4_array psum1;
    uint4_array scan1;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum1.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(psum0.a[ii]);      scan1.a[ii] = RADIX_CONVERT_8x4to16x2UPPER(scan0.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#pragma unroll
    for(int ii = 0; ii < 4; ii++){
      psum0.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(psum0.a[ii]);      scan0.a[ii] = RADIX_CONVERT_8x4to16x2LOWER(scan0.a[ii]);
    }/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_CHECK_BITS <= 2
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* prefix sum within a block */
    //---------------------------------------------------------------------
    /* clear the idx_dst for safety */
#ifndef USE_WARP_SHUFFLE_FUNC_SORT
    __syncthreads();
#endif//USE_WARP_SHUFFLE_FUNC_SORT
#pragma unroll
    for(int ii = tidx; ii < dataNum * sameNum; ii += NTHREADS_SORT)
      smem[ii].u = 0;
    __syncthreads();
    /* store the partial prefix sum to the shared memory */
    if( lane == (warpSize - 1) )
#   if  RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
      for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	//-----------------------------------------------------------------
	const int headIdx = warpIdx + ii * warpNum;
	//-----------------------------------------------------------------
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
	//-----------------------------------------------------------------
	smem[headIdx                                                 ].u = psum0;
	smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD      * warpNum].u = psum1;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
	//-----------------------------------------------------------------
	smem[headIdx                                                 ].u = psum0.a[ii];
	smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD      * warpNum].u = psum1.a[ii];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < 4; ii++){ */
#else///RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
      for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32; ii++){
	//-----------------------------------------------------------------
	const int headIdx = warpIdx + ii * 2 * sameNum;
	//-----------------------------------------------------------------
	smem[headIdx                                                 ].u = psum0.a[ii];
	smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD      * warpNum].u = psum1.a[ii];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
	smem[headIdx +                                        warpNum].u = psum2.a[ii];
	smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 1) * warpNum].u = psum3.a[ii];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	smem[headIdx +                                   2  * warpNum].u = psum4.a[ii];
	smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 2) * warpNum].u = psum5.a[ii];
	smem[headIdx +                                   3  * warpNum].u = psum6.a[ii];
	smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 3) * warpNum].u = psum7.a[ii];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_CHECK_BITS <= 2
    __syncthreads();
    //---------------------------------------------------------------------
    /* prefix sum within block (16bits * 2 is merged) */
    RADIX_PREFIX_SUM_BLCK(tidx, lane, smem
#ifdef  USE_MASK_INSTEAD_OF_IF
			  , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
			  );
    //---------------------------------------------------------------------
    /* memorize # of elements in each bucket */
    int num0, num1;
    if( tidx < (RADIX_SORT_NUM_BUCKET >> 1) ){
      num0 = RADIX_CONVERT_16x2to32x1LOWER(smem[-1 + (1 + tidx) * sameNum].i);/* # of elements for subkey =     2 * tidx */
      num1 = RADIX_CONVERT_16x2to32x1UPPER(smem[-1 + (1 + tidx) * sameNum].i);/* # of elements for subkey = 1 + 2 * tidx */
    }/* if( tidx < (RADIX_SORT_NUM_BUCKET >> 1) ){ */
    //---------------------------------------------------------------------
#ifdef  ENABLE_EARLY_EXIT_BLCK
    /* get # of elements in each bucket */
    /* 16bits value * 4 or 16 elements --> 32 bits value * 2 or 8 elements --> uint2 or 2xuint4 */
    /* unfortunately, additional variables are necessary */
    /* or, check 4 or 16 threads and accumulate results <-- this would be the best solution */
    int stock;
    if( tidx < (RADIX_SORT_NUM_BUCKET >> 1) ){
      //-------------------------------------------------------------------
      /* get masked flag to detect the signal */
      //-------------------------------------------------------------------
      /* get # of elements in each bucket */
      int flag = smem[(tidx + 1) * sameNum - 1].i;
      if( tidx > 0 )
	flag -= smem[tidx * sameNum - 1].i;
      /* apply the mask */
      flag &= mask;
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* calculate the total sum of all flags */
      //-------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
      if( tidx == 0 )
	stock = smem[0].i;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) >= 2
      int temp;
      temp = __shfl_up(flag, 1, RADIX_SORT_NUM_BUCKET >> 1);      flag += temp;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) >= 4
      temp = __shfl_up(flag, 2, RADIX_SORT_NUM_BUCKET >> 1);      flag += temp;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) == 8
      temp = __shfl_up(flag, 4, RADIX_SORT_NUM_BUCKET >> 1);      flag += temp;
#endif//(RADIX_SORT_NUM_BUCKET >> 1) == 8
#endif//(RADIX_SORT_NUM_BUCKET >> 1) >= 4
#endif//(RADIX_SORT_NUM_BUCKET >> 1) >= 2
      if( tidx == 0 )
	smem[0].i = flag;
#else///USE_WARP_SHUFFLE_FUNC_SORT
      stock = smem[tidx].i;
      smem[tidx].i = flag;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) >= 2
      int temp;
      temp = smem[tidx ^ 1];      flag += temp;      smem[tidx].i = flag;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) >= 4
      temp = smem[tidx ^ 2];      flag += temp;      smem[tidx].i = flag;
#   if  (RADIX_SORT_NUM_BUCKET >> 1) == 8
      temp = smem[tidx ^ 4];      flag += temp;      smem[tidx].i = flag;
#endif//(RADIX_SORT_NUM_BUCKET >> 1) == 8
#endif//(RADIX_SORT_NUM_BUCKET >> 1) >= 4
#endif//(RADIX_SORT_NUM_BUCKET >> 1) >= 2
#endif//USE_WARP_SHUFFLE_FUNC_SORT
      //-------------------------------------------------------------------
    }/* if( tidx < (RADIX_SORT_NUM_BUCKET >> 1) ){ */
    __syncthreads();
    const int skip = smem[0].i;
    //---------------------------------------------------------------------
    /* if all flags are equal to zero; then sorting must be processed */
    /* if one of flags is not equal to zero; then sorting can be skipped */
    //---------------------------------------------------------------------
    if( skip == 0 )
#endif//ENABLE_EARLY_EXIT_BLCK
      {
	//-----------------------------------------------------------------
#ifdef  ENABLE_EARLY_EXIT_BLCK
	__syncthreads();
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	if( tidx == 0 )
	  smem[0].i = stock;
#else///USE_WARP_SHUFFLE_FUNC_SORT
	if( tidx < (RADIX_SORT_NUM_BUCKET >> 1) )
	  smem[tidx].i = stock;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
	__syncthreads();
#endif//ENABLE_EARLY_EXIT_BLCK
	//-----------------------------------------------------------------
	/* prefix sum within block (complete set) */
	/* 16bits shift toward left corresponds to multiplying 2^16 = 65536; */
	const uint headLower = (smem[sameNum - 1].u & 0xffff) << 16;
	int val;
	if( dataIdx < dataNum ){
	  val  =                 (smem[(dataIdx + 1) * sameNum - 1].i & 0xffff) << 16;
	  val += (dataIdx > 0) ? (smem[ dataIdx      * sameNum - 1].i           >> 16): 0;
	}/* if( dataIdx < dataNum ){ */
	__syncthreads();
	if( tidx < dataNum * sameNum )
	  smem[tidx].i += val;
	__syncthreads();
	/* get prefix sum as a complete set */
#   if  RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
	for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	  //---------------------------------------------------------------
	  /* const int headIdx = (warpIdx - 1) + ii * 2 * sameNum; */
	  const int headIdx = (warpIdx - 1) + ii * warpNum;
	  //---------------------------------------------------------------
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
	  psum0       += (headIdx >= 0) ? (smem[headIdx].u) : (headLower);	  psum1       += smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD * warpNum].u;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
	  psum0.a[ii] += (headIdx >= 0) ? (smem[headIdx].u) : (headLower);	  psum1.a[ii] += smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD * warpNum].u;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
          //---------------------------------------------------------------
	}/* for(int ii = 0; ii < 4; ii++){ */
#else///RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
	for(int ii = 0; ii < RADIX_SORT_NUM_OF_I32; ii++){
	  //---------------------------------------------------------------
	  const int headIdx = (warpIdx - 1) + ii * 2 * sameNum;
	  //---------------------------------------------------------------
	  psum0.a[ii] += (headIdx >= 0) ? (smem[headIdx].u) : (headLower);
	  psum1.a[ii] += smem[headIdx +  RADIX_SORT_ELEMENTS_PER_THREAD      * warpNum].u;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
	  psum2.a[ii] += smem[headIdx +                                        warpNum].u;
	  psum3.a[ii] += smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 1) * warpNum].u;
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	  psum4.a[ii] += smem[headIdx +                                   2  * warpNum].u;
	  psum5.a[ii] += smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 2) * warpNum].u;
	  psum6.a[ii] += smem[headIdx +                                   3  * warpNum].u;
	  psum7.a[ii] += smem[headIdx + (RADIX_SORT_ELEMENTS_PER_THREAD + 3) * warpNum].u;
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
          //---------------------------------------------------------------
	}/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_CHECK_BITS <= 2
	//-----------------------------------------------------------------
	/* set exclusive prefix sum */
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
	psum0 -= scan0;
	psum1 -= scan1;
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
#pragma unroll
	for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	  psum0.a[ii] -= scan0.a[ii];
	  psum1.a[ii] -= scan1.a[ii];
	}/* for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){ */
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#pragma unroll
	for(int ii = 0; ii < 4; ii++){
	  psum0.a[ii] -= scan0.a[ii];
	  psum1.a[ii] -= scan1.a[ii];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD >= 2
	  psum2.a[ii] -= scan2.a[ii];
	  psum3.a[ii] -= scan3.a[ii];
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	  psum4.a[ii] -= scan4.a[ii];
	  psum5.a[ii] -= scan5.a[ii];
	  psum6.a[ii] -= scan6.a[ii];
	  psum7.a[ii] -= scan7.a[ii];
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#endif//RADIX_SORT_ELEMENTS_PER_THREAD >= 2
	}/* for(int ii = 0; ii < 4; ii++){ */
#endif//RADIX_SORT_CHECK_BITS <= 2
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* execute bucket sort within a block */
	//-----------------------------------------------------------------
#   if  RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
	subkey = RADIX_SORT_BLCK_GET_DSTIDX(subkey  , psum0, psum1);
	key[subkey] = key_src;
#ifndef SORT_ONLY
	idx[subkey] = idx_src;
#endif//SORT_ONLY
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
	subkey.x = RADIX_SORT_BLCK_GET_DSTIDX(subkey.x, psum0.m.x, psum1.m.x);
	subkey.y = RADIX_SORT_BLCK_GET_DSTIDX(subkey.y, psum0.m.y, psum1.m.y);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	subkey.z = RADIX_SORT_BLCK_GET_DSTIDX(subkey.z, psum0.m.z, psum1.m.z);
	subkey.w = RADIX_SORT_BLCK_GET_DSTIDX(subkey.w, psum0.m.w, psum1.m.w);
	key[subkey.z] = key_src.z;
	key[subkey.w] = key_src.w;
#ifndef SORT_ONLY
	idx[subkey.z] = idx_src.z;
	idx[subkey.w] = idx_src.w;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
#if 0
	if( subkey.x > NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD )
	  printf("lane = %d, tidx = %d, iter = %d, subkey.x = %d, key_src.x = %x\n", lane, tidx, iter, subkey.x, key_src.x);
#endif
	key[subkey.x] = key_src.x;
	key[subkey.y] = key_src.y;
#ifndef SORT_ONLY
	idx[subkey.x] = idx_src.x;
	idx[subkey.y] = idx_src.y;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#else///RADIX_SORT_CHECK_BITS <= 2
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
	subkey = RADIX_SORT_BLCK_GET_DSTIDX(subkey  , psum0.m, psum1.m, &sbuf[tidx]);
	key[subkey] = key_src;
#ifndef SORT_ONLY
	idx[subkey] = idx_src;
#endif//SORT_ONLY
#else///RADIX_SORT_ELEMENTS_PER_THREAD == 1
	subkey.x = RADIX_SORT_BLCK_GET_DSTIDX(subkey.x, psum0.m, psum1.m, &sbuf[tidx]);
	subkey.y = RADIX_SORT_BLCK_GET_DSTIDX(subkey.y, psum2.m, psum3.m, &sbuf[tidx]);
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
	subkey.z = RADIX_SORT_BLCK_GET_DSTIDX(subkey.z, psum4.m, psum5.m, &sbuf[tidx]);
	subkey.w = RADIX_SORT_BLCK_GET_DSTIDX(subkey.w, psum6.m, psum7.m, &sbuf[tidx]);
	key[subkey.z] = key_src.z;
	key[subkey.w] = key_src.w;
#ifndef SORT_ONLY
	idx[subkey.z] = idx_src.z;
	idx[subkey.w] = idx_src.w;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
	key[subkey.x] = key_src.x;
	key[subkey.y] = key_src.y;
#ifndef SORT_ONLY
	idx[subkey.x] = idx_src.x;
	idx[subkey.y] = idx_src.y;
#endif//SORT_ONLY
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
#endif//RADIX_SORT_CHECK_BITS <= 2
	//-----------------------------------------------------------------
      }
    __syncthreads();
    //---------------------------------------------------------------------
    if( tidx < (RADIX_SORT_NUM_BUCKET >> 1) ){
      //-------------------------------------------------------------------
      smem[STOCK_ARRAY_SUM_HEAD + 2 * tidx    ].i = num0;
      smem[STOCK_ARRAY_SUM_HEAD + 2 * tidx + 1].i = num1;
      //-------------------------------------------------------------------
    }/* if( tidx < (RADIX_SORT_NUM_BUCKET >> 1) ){ */
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* MSD radix sort within device */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_GRID
(int * RESTRICT numIterFul, int2 * RESTRICT infoFul, int remFul, int * RESTRICT numIterLoc, int2 * RESTRICT infoLoc, const int remLocMax,
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
 int * RESTRICT checkedCounts,
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
 int * RESTRICT fail_dev, const int tidx, const int lane,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
 const int hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
 RADIX_DATA_TYPE_KEY * RESTRICT key_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_tmp_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_sm,
#ifndef SORT_ONLY
 int * RESTRICT idx_gm, int * RESTRICT idx_tmp_gm, int * RESTRICT idx_sm,
#endif//SORT_ONLY
#   if  RADIX_SORT_CHECK_BITS > 2
 uint4_array * RESTRICT sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
 volatile uint_int * RESTRICT smem, int * RESTRICT scan_gm, int * RESTRICT bucket_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul,
 const int gidx, const int bidx, const int bnum, volatile int * RESTRICT gsync0, volatile int * RESTRICT gsync1
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
)
{
  //-----------------------------------------------------------------------
  int remLoc = 0;
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
  *checkedCounts = 0;
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
  //-----------------------------------------------------------------------
  /* MSD radix sort within device using inter-block synchronization */
  //-----------------------------------------------------------------------
  while( remFul > 0 ){
    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------
    /* pull data from stack (numIterFul, infoFul) in LIFO (Last In First Out) manner */
    const int stackIdx = remFul - 1;
    const int numIter = numIterFul[stackIdx];
    const int2 info = infoFul[stackIdx];
    //---------------------------------------------------------------------
    /* evaluate # of iterations to check all elements contained in a current data set */
    const int head = info.y;
    const int  num = info.x;
    const int nloop = BLOCKSIZE(num, bnum * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
    //---------------------------------------------------------------------
    for(int loop = 0; loop < nloop; loop++){
      //-------------------------------------------------------------------
      const int hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx + loop * bnum);
      //-------------------------------------------------------------------
      int subNum = num - (hidx_gm - head);
      if( subNum < 0 )	subNum = 0;
      subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
      //-------------------------------------------------------------------
      __syncthreads();
      //-------------------------------------------------------------------
      /* a remedy to treat the case when the number of elements is not a multiple of NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD */
#pragma unroll
      for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	//-----------------------------------------------------------------
#   if  KEY_BITS == 32
	key_sm[tidx + ii * NTHREADS_SORT] = 0xffffffff;
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
	key_sm[tidx + ii * NTHREADS_SORT] = 0xffffffffffffffff;
#endif//KEY_BITS == 64
	//-----------------------------------------------------------------
      }/* for(int ii = subNum + tidx; ii < RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      /* load key_src, idx_src from global memory */
      /* set key_src and idx_src on shared memory for RADIX_SORT_CHIP */
#pragma unroll
      for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	//-----------------------------------------------------------------
	const int idx = hidx_gm + ii;
	//-----------------------------------------------------------------
	key_sm[ii] = key_gm[idx];
#ifndef SORT_ONLY
	idx_sm[ii] = idx_gm[idx];
#endif//SORT_ONLY
	//-----------------------------------------------------------------
      }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      __syncthreads();
      RADIX_SORT_CHIP(numIter - 1, tidx,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
		      hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
		      key_sm
#ifndef SORT_ONLY
		      , idx_sm
#endif//SORT_ONLY
		      , smem
#   if  RADIX_SORT_CHECK_BITS > 2
		      , sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
		      , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
		      );
      //-------------------------------------------------------------------
      /* remove the number of ghost values */
      /* target == RADIX_SORT_CHECK_MASK (= RADIX_SORT_NUM_BUCKET - 1) */
      /* # of ghosts is RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum */
      if( tidx == RADIX_SORT_CHECK_MASK )
	smem[STOCK_ARRAY_SUM_HEAD + tidx].i -= RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum;
      //-------------------------------------------------------------------
      /* put returned values of num_sm on global memory to calculate global prefix sum */
      if( tidx < RADIX_SORT_NUM_BUCKET )
      	scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, tidx, bidx + bnum * loop)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
      //-------------------------------------------------------------------
      /* if nloop > 1 && loop < nloop - 1; save num_sm, key_dst and idx_dst on the global memory */
      if( loop < (nloop - 1) ){
	//-----------------------------------------------------------------
	/* save num_sm on the global memory */
	if( tidx < RADIX_SORT_NUM_BUCKET )
	  scan_gm[INDEX2D(nloop * bnum, RADIX_SORT_NUM_BUCKET, bidx + bnum * loop, tidx)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i;
	//-----------------------------------------------------------------
	/* save key_dst and idx_dst on the global memory */
#pragma unroll
	for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	  key_tmp_gm[hidx_gm + ii] = key_sm[ii];
#ifndef SORT_ONLY
	  idx_tmp_gm[hidx_gm + ii] = idx_sm[ii];
#endif//SORT_ONLY
	}/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
	//-----------------------------------------------------------------
      }/* if( loop < (nloop - 1) ){ */
      //-------------------------------------------------------------------
    }/* for(int loop = 0; loop < nloop; loop++){ */
    //---------------------------------------------------------------------
    /* calculate global prefix sum within a device */
    /* if # of elements is 128M --> # of elements per block is 128M / (4 * Ttot) ~ 512k >> 1024. */
    /* Therefore, simple evaluation of prefix sum within a block is not sufficient */
    int tmp_sum_head;
    if( tidx < RADIX_SORT_NUM_BUCKET )
      tmp_sum_head = smem[STOCK_ARRAY_SUM_HEAD + tidx].i;
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    for(int target = bidx; target < RADIX_SORT_NUM_BUCKET; target += bnum){
      //-------------------------------------------------------------------
      int subTail = 0;
      const int niter = BLOCKSIZE(bnum * nloop, NTHREADS_SORT);
      //-------------------------------------------------------------------
      for(int iter = 0; iter < niter; iter++){
	//-----------------------------------------------------------------
	const int sidx = tidx + iter * NTHREADS_SORT;
	//-----------------------------------------------------------------
	/* load local prefix sum */
	smem[tidx].i = (sidx < bnum * nloop) ? scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, target, sidx)] : 0;
	/* prefix sum within a block */
	RADIX_PREFIX_SUM_GRID(tidx, lane, smem
#ifdef  USE_MASK_INSTEAD_OF_IF
			      , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
			      );
	/* write back global memory */
	if( sidx < bnum * nloop )
	  scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, target, sidx)] = subTail + smem[tidx].i;
	subTail += smem[NTHREADS_SORT - 1].i;
	__syncthreads();
	//-----------------------------------------------------------------
      }/* for(int iter = 0; iter < niter; iter++){ */
      //-------------------------------------------------------------------
      if( tidx == 0 ){
	tail_gm  [target] = subTail;
	/* 1 means that further global sort is required */
	bucket_gm[target] = (subTail < (NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD)) ? 0 : 1;
      }/* if( tidx == 0 ){ */
      //-------------------------------------------------------------------
    }/* for(int target = bidx; target < RADIX_SORT_NUM_BUCKET; target += bnum){ */
    if( tidx < RADIX_SORT_NUM_BUCKET )
      smem[STOCK_ARRAY_SUM_HEAD + tidx].i = tmp_sum_head;
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------
    remFul--;
    if( tidx < RADIX_SORT_NUM_BUCKET ){
      //-------------------------------------------------------------------
      /* get partition among sub-groups */
      const int snum = tail_gm[tidx];
      const int sidx = STOCK_ARRAY_PTR_HEAD + tidx;
      /* 1 means that further global sort is required */
      const int bucket_tmp = bucket_gm[tidx];
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
      /* load index */
      int v0 = snum;
      int v1 = bucket_tmp;
      int t0, t1;
      /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET) & m0;      t1 = __shfl_up(v1, 1, RADIX_SORT_NUM_BUCKET) & m0;      v0 += t0;      v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >=  4
      t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET) & m1;      t1 = __shfl_up(v1, 2, RADIX_SORT_NUM_BUCKET) & m1;      v0 += t0;      v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >=  8
      t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET) & m2;      t1 = __shfl_up(v1, 4, RADIX_SORT_NUM_BUCKET) & m2;      v0 += t0;      v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >= 16
      t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET) & m3;      t1 = __shfl_up(v1, 8, RADIX_SORT_NUM_BUCKET) & m3;      v0 += t0;      v1 += t1;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 1, RADIX_SORT_NUM_BUCKET);      if( tidx >= 1 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >=  4
      t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 2, RADIX_SORT_NUM_BUCKET);      if( tidx >= 2 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >=  8
      t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 4, RADIX_SORT_NUM_BUCKET);      if( tidx >= 4 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >= 16
      t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 8, RADIX_SORT_NUM_BUCKET);      if( tidx >= 8 ){	v0 += t0;	v1 += t1;      }
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
      /* return calculated inclusive prefix sum */
      smem[sidx].i = v0;
      smem[tidx].i = v1;
#else///USE_WARP_SHUFFLE_FUNC_SORT
      smem[sidx].i = snum;
      smem[tidx].i = bucket_tmp;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      smem[tidx].i += smem[tidx - 1].i & m0;      smem[sidx].i += smem[sidx - 1].i & m0;
#   if  RADIX_SORT_NUM_BUCKET >=  4
      smem[tidx].i += smem[tidx - 2].i & m1;      smem[sidx].i += smem[sidx - 2].i & m1;
#   if  RADIX_SORT_NUM_BUCKET >=  8
      smem[tidx].i += smem[tidx - 4].i & m2;      smem[sidx].i += smem[sidx - 4].i & m2;
#   if  RADIX_SORT_NUM_BUCKET >= 16
      smem[tidx].i += smem[tidx - 8].i & m3;      smem[sidx].i += smem[sidx - 8].i & m3;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      if( tidx >= 1 ){	smem[tidx].i += smem[tidx - 1].i;	smem[sidx].i += smem[sidx - 1].i;      }
#   if  RADIX_SORT_NUM_BUCKET >=  4
      if( tidx >= 2 ){	smem[tidx].i += smem[tidx - 2].i;	smem[sidx].i += smem[sidx - 2].i;      }
#   if  RADIX_SORT_NUM_BUCKET >=  8
      if( tidx >= 4 ){	smem[tidx].i += smem[tidx - 4].i;	smem[sidx].i += smem[sidx - 4].i;      }
#   if  RADIX_SORT_NUM_BUCKET >= 16
      if( tidx >= 8 ){	smem[tidx].i += smem[tidx - 8].i;	smem[sidx].i += smem[sidx - 8].i;      }
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#endif//USE_WARP_SHUFFLE_FUNC_SORT
      //-------------------------------------------------------------------
      smem[STOCK_ARRAY_DST_HEAD + tidx].i = bucket_tmp;
      smem[STOCK_ARRAY_IDX_HEAD + tidx].i = bucket_tmp ? (smem[tidx].i - bucket_tmp) : (tidx - smem[tidx].i);
      //-------------------------------------------------------------------
      smem[STOCK_ARRAY_PTR_NUM  + tidx].i  = snum;
      smem[STOCK_ARRAY_PTR_HEAD + tidx].i -= snum;
      //-------------------------------------------------------------------
      if( bidx == 0 ){
	int2 subInfo = {snum, head + smem[STOCK_ARRAY_PTR_HEAD + tidx].i};
	if( smem[STOCK_ARRAY_DST_HEAD + tidx].i ){
	  infoFul   [remFul + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = subInfo;
	  numIterFul[remFul + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = numIter - 1;
	}
	else{
	  infoLoc   [remLoc + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = subInfo;
	  numIterLoc[remLoc + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = numIter - 1;
	}
      }/* if( bidx == 0 ){ */
      //-------------------------------------------------------------------
    }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
    __syncthreads();
    remFul +=                         smem[RADIX_SORT_NUM_BUCKET - 1].i;
    remLoc += RADIX_SORT_NUM_BUCKET - smem[RADIX_SORT_NUM_BUCKET - 1].i;
    //---------------------------------------------------------------------
    int hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx + (nloop - 1) * bnum);
    //---------------------------------------------------------------------
    int subNum = num - (hidx_gm - head);
    if( subNum < 0 )      subNum = 0;
    subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
    for(int loop = nloop - 1; loop >= 0; loop--){
      //-------------------------------------------------------------------
      __syncthreads();
      if( tidx < RADIX_SORT_NUM_BUCKET ){
	//-----------------------------------------------------------------
	if( loop != (nloop - 1) )
	  smem[STOCK_ARRAY_SUM_HEAD + tidx].i = scan_gm[INDEX2D(nloop * bnum, RADIX_SORT_NUM_BUCKET, bidx + bnum * loop, tidx)];
	//-----------------------------------------------------------------
	/* load global prefix sum */
	smem[tidx].i = scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, tidx, bidx + loop * bnum)];
	//-----------------------------------------------------------------
	/* calculate local prefix sum of each element */
	smem[STOCK_ARRAY_NUM_HEAD + tidx].i = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
	const int tmpVal = (tidx != 0) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 1].i;
	smem[STOCK_ARRAY_SUM_HEAD + tidx].i += tmpVal;
	//-----------------------------------------------------------------
	smem[tidx].i += (smem[STOCK_ARRAY_PTR_HEAD + tidx].i - smem[STOCK_ARRAY_SUM_HEAD + tidx].i);
	//-----------------------------------------------------------------
      }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
      __syncthreads();
      //-------------------------------------------------------------------
      /* store sorted key, idx on global memory */
#pragma unroll
      for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	//-----------------------------------------------------------------
	const RADIX_DATA_TYPE_KEY key = key_sm[ii];
	const int target = (key >> ((numIter - 1) * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK;
	//-----------------------------------------------------------------
	const int dstIdx = head + smem[target].i + ii;
	key_gm[dstIdx] = key;
#ifndef SORT_ONLY
	idx_gm[dstIdx] = idx_sm[ii];
#endif//SORT_ONLY
	//-----------------------------------------------------------------
      }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      /* load temporary stock on the global memory (num_sm, key_dst_sm and idx_dst_sm) */
      if( loop != 0 ){
	//-----------------------------------------------------------------
	hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx + (loop - 1) * bnum);
	subNum = RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT;
	//-----------------------------------------------------------------
#pragma unroll
	for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	  key_sm[tidx + ii * NTHREADS_SORT] = key_tmp_gm[hidx_gm + tidx + ii * NTHREADS_SORT];
#ifndef SORT_ONLY
	  idx_sm[tidx + ii * NTHREADS_SORT] = idx_tmp_gm[hidx_gm + tidx + ii * NTHREADS_SORT];
#endif//SORT_ONLY
	}/* for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){ */
	//-----------------------------------------------------------------
      }/* if( loop != 0 ){ */
      //-------------------------------------------------------------------
    }/* for(int loop = nloop - 1; loop >= 0; loop--){ */
    //---------------------------------------------------------------------
  }/* while( remFul > 0 ){ */
  //-----------------------------------------------------------------------
  if( (gidx == 0) && (remLoc > remLocMax) )
    *fail_dev = remLoc;
  //-----------------------------------------------------------------------
  globalSync(tidx, bidx, bnum, gsync0, gsync1);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* LSD radix sort within a block */
  //-----------------------------------------------------------------------
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
  if( tidx == 0 )
    smem[0].i = atomicAdd(checkedCounts, 1);
  __syncthreads();
  int target = smem[0].i;
  while( target < remLoc )
#else///ATOMIC_BASED_JOB_ASSIGNMENT
  for(int target = bidx; target < remLoc; target += bnum)
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
    {
      //-------------------------------------------------------------------
      /* pull data from queue (numIterLoc, infoLoc) in FIFO (First In First Out) manner */
      if( tidx == 0 ){
	int2 info = infoLoc[target];
	smem[0].i = numIterLoc[target];
	smem[1].i = info.x;
	smem[2].i = info.y;
      }
      __syncthreads();
      const int numIter = smem[0].i;
      const int  rem = smem[1].i;
      const int head = smem[2].i;
      //-------------------------------------------------------------------
      /* a remedy to treat the case when the number of elements is not a multiple of NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD */
#pragma unroll
      for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	//-----------------------------------------------------------------
#   if  KEY_BITS == 32
	key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffff;
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
	key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffffffffffff;
#endif//KEY_BITS == 64
	//-----------------------------------------------------------------
      }/* for(int ii = subNum + tidx; ii < RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      /* load data set from global memory */
#pragma unroll
      for(int ii = tidx; ii < rem; ii += NTHREADS_SORT){
	//-----------------------------------------------------------------
	key_sm[ii] = key_gm[head + ii];
#ifndef SORT_ONLY
	idx_sm[ii] = idx_gm[head + ii];
#endif//SORT_ONLY
        //-----------------------------------------------------------------
      }/* for(int ii = tidx; ii < numIter.x; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      __syncthreads();
      //-------------------------------------------------------------------
      /* execute LSD radix sort within a block */
      RADIX_SORT_BLCK
	(numIter, tidx,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
	 hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
	 key_sm
#ifndef SORT_ONLY
	 , idx_sm
#endif//SORT_ONLY
	 , smem
#   if  RADIX_SORT_CHECK_BITS > 2
	 , sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
	 , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
	 );
      //-------------------------------------------------------------------
      /* write back data set to the global memory */
#pragma unroll
      for(int ii = tidx; ii < rem; ii += NTHREADS_SORT){
	//-----------------------------------------------------------------
	key_gm[head + ii] = key_sm[ii];
#ifndef SORT_ONLY
	idx_gm[head + ii] = idx_sm[ii];
#endif//SORT_ONLY
        //-----------------------------------------------------------------
      }/* for(int ii = tidx; ii < numIter.x; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
      if( tidx == 0 )
	smem[0].i = atomicAdd(checkedCounts, 1);
      __syncthreads();
      target = smem[0].i;
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
      //-------------------------------------------------------------------
    }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* MSD radix sort within device (1st half) */
/* MSD radix sort within a grid */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_GRID_MSD
(int * RESTRICT numIterFul, int2 * RESTRICT infoFul, int remFul, int * RESTRICT numIterLoc, int2 * RESTRICT infoLoc, const int remLocMax,
 int * RESTRICT remLoc_dev, int * RESTRICT fail_dev, const int tidx, const int lane,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
 const int hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
 RADIX_DATA_TYPE_KEY * RESTRICT key_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_tmp_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_sm,
#ifndef SORT_ONLY
 int * RESTRICT idx_gm, int * RESTRICT idx_tmp_gm, int * RESTRICT idx_sm,
#endif//SORT_ONLY
#   if  RADIX_SORT_CHECK_BITS > 2
 uint4_array * RESTRICT sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
 volatile uint_int * RESTRICT smem, int * RESTRICT scan_gm, int * RESTRICT bucket_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul,
 const int gidx, const int bidx, const int bnum, volatile int * RESTRICT gsync0, volatile int * RESTRICT gsync1
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
)
{
  //-----------------------------------------------------------------------
  int remLoc = 0;
  //-----------------------------------------------------------------------
  /* MSD radix sort within device using inter-block synchronization */
  //-----------------------------------------------------------------------
  while( remFul > 0 ){
    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------
    /* pull data from stack (numIterFul, infoFul) in LIFO (Last In First Out) manner */
    const int stackIdx = remFul - 1;
    const int numIter = numIterFul[stackIdx];
    const int2 info = infoFul[stackIdx];
    //---------------------------------------------------------------------
    /* evaluate # of iterations to check all elements contained in a current data set */
    const int head = info.y;
    const int  num = info.x;
    const int nloop = BLOCKSIZE(num, bnum * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
    //---------------------------------------------------------------------
    for(int loop = 0; loop < nloop; loop++){
      //-------------------------------------------------------------------
      const int hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx + loop * bnum);
      //-------------------------------------------------------------------
      int subNum = num - (hidx_gm - head);
      if( subNum < 0 )	subNum = 0;
      subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
      //-------------------------------------------------------------------
      __syncthreads();
      //-------------------------------------------------------------------
      /* a remedy to treat the case when the number of elements is not a multiple of NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD */
#pragma unroll
      for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	//-----------------------------------------------------------------
#   if  KEY_BITS == 32
	key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffff;
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
	key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffffffffffff;
#endif//KEY_BITS == 64
	//-----------------------------------------------------------------
      }/* for(int ii = subNum + tidx; ii < RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      /* load key_src, idx_src from global memory */
      /* set key_src and idx_src on shared memory for RADIX_SORT_CHIP */
#pragma unroll
      for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	//-----------------------------------------------------------------
	const int idx = hidx_gm + ii;
	//-----------------------------------------------------------------
	key_sm[ii] = key_gm[idx];
#ifndef SORT_ONLY
	idx_sm[ii] = idx_gm[idx];
#endif//SORT_ONLY
	//-----------------------------------------------------------------
      }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      __syncthreads();
      RADIX_SORT_CHIP(numIter - 1, tidx,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
		      hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
		      key_sm
#ifndef SORT_ONLY
		      , idx_sm
#endif//SORT_ONLY
		      , smem
#   if  RADIX_SORT_CHECK_BITS > 2
		      , sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
		      , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
		      );
      //-------------------------------------------------------------------
      /* remove the number of ghost values */
      /* target == RADIX_SORT_CHECK_MASK (= RADIX_SORT_NUM_BUCKET - 1) */
      /* # of ghosts is RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum */
      if( tidx == RADIX_SORT_CHECK_MASK )
	smem[STOCK_ARRAY_SUM_HEAD + tidx].i -= RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum;
      //-------------------------------------------------------------------
      /* put returned values of num_sm on global memory to calculate global prefix sum */
      if( tidx < RADIX_SORT_NUM_BUCKET )
      	scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, tidx, bidx + bnum * loop)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
      //-------------------------------------------------------------------
      /* if nloop > 1 && loop < nloop - 1; save num_sm, key_dst and idx_dst on the global memory */
      if( loop < (nloop - 1) ){
	//-----------------------------------------------------------------
	/* save num_sm on the global memory */
	if( tidx < RADIX_SORT_NUM_BUCKET )
	  scan_gm[INDEX2D(nloop * bnum, RADIX_SORT_NUM_BUCKET, bidx + bnum * loop, tidx)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i;
	//-----------------------------------------------------------------
	/* save key_dst and idx_dst on the global memory */
#pragma unroll
	for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	  key_tmp_gm[hidx_gm + ii] = key_sm[ii];
#ifndef SORT_ONLY
	  idx_tmp_gm[hidx_gm + ii] = idx_sm[ii];
#endif//SORT_ONLY
	}/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
	//-----------------------------------------------------------------
      }/* if( loop < (nloop - 1) ){ */
      //-------------------------------------------------------------------
      /* rem -= totNum; */
      //-------------------------------------------------------------------
    }/* for(int loop = 0; loop < nloop; loop++){ */
    //---------------------------------------------------------------------
    /* calculate global prefix sum within a device */
    /* if # of elements is 128M --> # of elements per block is 128M / (4 * Ttot) ~ 512k >> 1024. */
    /* Therefore, simple evaluation of prefix sum within a block is not sufficient */
    int tmp_sum_head;
    if( tidx < RADIX_SORT_NUM_BUCKET )
      tmp_sum_head = smem[STOCK_ARRAY_SUM_HEAD + tidx].i;
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    for(int target = bidx; target < RADIX_SORT_NUM_BUCKET; target += bnum){
      //-------------------------------------------------------------------
      int subTail = 0;
      const int niter = BLOCKSIZE(bnum * nloop, NTHREADS_SORT);
      //-------------------------------------------------------------------
      for(int iter = 0; iter < niter; iter++){
	//-----------------------------------------------------------------
	const int sidx = tidx + iter * NTHREADS_SORT;
	//-----------------------------------------------------------------
	/* load local prefix sum */
	smem[tidx].i = (sidx < bnum * nloop) ? scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, target, sidx)] : 0;
	/* prefix sum within a block */
	RADIX_PREFIX_SUM_GRID(tidx, lane, smem
#ifdef  USE_MASK_INSTEAD_OF_IF
			      , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
			      );
	/* write back global memory */
	if( sidx < bnum * nloop )
	  scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, target, sidx)] = subTail + smem[tidx].i;
	subTail += smem[NTHREADS_SORT - 1].i;
	__syncthreads();
	//-----------------------------------------------------------------
      }/* for(int iter = 0; iter < niter; iter++){ */
      //-------------------------------------------------------------------
      if( tidx == 0 ){
	tail_gm  [target] = subTail;
	/* 1 means that further global sort is required */
	bucket_gm[target] = (subTail < (NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD)) ? 0 : 1;
      }/* if( tidx == 0 ){ */
      //-------------------------------------------------------------------
    }/* for(int target = bidx; target < RADIX_SORT_NUM_BUCKET; target += bnum){ */
    if( tidx < RADIX_SORT_NUM_BUCKET )
      smem[STOCK_ARRAY_SUM_HEAD + tidx].i = tmp_sum_head;
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------
    remFul--;
    if( tidx < RADIX_SORT_NUM_BUCKET ){
      //-------------------------------------------------------------------
      /* get partition among sub-groups */
      const int snum = tail_gm[tidx];
      const int sidx = STOCK_ARRAY_PTR_HEAD + tidx;
      /* 1 means that further global sort is required */
      const int bucket_tmp = bucket_gm[tidx];
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
      /* load index */
      int v0 = snum;
      int v1 = bucket_tmp;
      int t0, t1;
      /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET) & m0;      t1 = __shfl_up(v1, 1, RADIX_SORT_NUM_BUCKET) & m0;      v0 += t0;      v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >=  4
      t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET) & m1;      t1 = __shfl_up(v1, 2, RADIX_SORT_NUM_BUCKET) & m1;      v0 += t0;      v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >=  8
      t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET) & m2;      t1 = __shfl_up(v1, 4, RADIX_SORT_NUM_BUCKET) & m2;      v0 += t0;      v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >= 16
      t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET) & m3;      t1 = __shfl_up(v1, 8, RADIX_SORT_NUM_BUCKET) & m3;      v0 += t0;      v1 += t1;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 1, RADIX_SORT_NUM_BUCKET);      if( tidx >= 1 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >=  4
      t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 2, RADIX_SORT_NUM_BUCKET);      if( tidx >= 2 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >=  8
      t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 4, RADIX_SORT_NUM_BUCKET);      if( tidx >= 4 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >= 16
      t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 8, RADIX_SORT_NUM_BUCKET);      if( tidx >= 8 ){	v0 += t0;	v1 += t1;      }
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
      /* return calculated inclusive prefix sum */
      smem[sidx].i = v0;
      smem[tidx].i = v1;
#else///USE_WARP_SHUFFLE_FUNC_SORT
      smem[sidx].i = snum;
      smem[tidx].i = bucket_tmp;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      smem[tidx].i += smem[tidx - 1].i & m0;      smem[sidx].i += smem[sidx - 1].i & m0;
#   if  RADIX_SORT_NUM_BUCKET >=  4
      smem[tidx].i += smem[tidx - 2].i & m1;      smem[sidx].i += smem[sidx - 2].i & m1;
#   if  RADIX_SORT_NUM_BUCKET >=  8
      smem[tidx].i += smem[tidx - 4].i & m2;      smem[sidx].i += smem[sidx - 4].i & m2;
#   if  RADIX_SORT_NUM_BUCKET >= 16
      smem[tidx].i += smem[tidx - 8].i & m3;      smem[sidx].i += smem[sidx - 8].i & m3;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      if( tidx >= 1 ){	smem[tidx].i += smem[tidx - 1].i;	smem[sidx].i += smem[sidx - 1].i;      }
#   if  RADIX_SORT_NUM_BUCKET >=  4
      if( tidx >= 2 ){	smem[tidx].i += smem[tidx - 2].i;	smem[sidx].i += smem[sidx - 2].i;      }
#   if  RADIX_SORT_NUM_BUCKET >=  8
      if( tidx >= 4 ){	smem[tidx].i += smem[tidx - 4].i;	smem[sidx].i += smem[sidx - 4].i;      }
#   if  RADIX_SORT_NUM_BUCKET >= 16
      if( tidx >= 8 ){	smem[tidx].i += smem[tidx - 8].i;	smem[sidx].i += smem[sidx - 8].i;      }
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#endif//USE_WARP_SHUFFLE_FUNC_SORT
      //-------------------------------------------------------------------
      smem[STOCK_ARRAY_DST_HEAD + tidx].i = bucket_tmp;
      smem[STOCK_ARRAY_IDX_HEAD + tidx].i = bucket_tmp ? (smem[tidx].i - bucket_tmp) : (tidx - smem[tidx].i);
      //-------------------------------------------------------------------
      smem[STOCK_ARRAY_PTR_NUM  + tidx].i  = snum;
      smem[STOCK_ARRAY_PTR_HEAD + tidx].i -= snum;
      //-------------------------------------------------------------------
      if( bidx == 0 ){
	int2 subInfo = {snum, head + smem[STOCK_ARRAY_PTR_HEAD + tidx].i};
	if( smem[STOCK_ARRAY_DST_HEAD + tidx].i ){
	  infoFul   [remFul + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = subInfo;
	  numIterFul[remFul + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = numIter - 1;
	}
	else{
	  infoLoc   [remLoc + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = subInfo;
	  numIterLoc[remLoc + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = numIter - 1;
	}
      }/* if( bidx == 0 ){ */
      //-------------------------------------------------------------------
    }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
    __syncthreads();
    remFul +=                         smem[RADIX_SORT_NUM_BUCKET - 1].i;
    remLoc += RADIX_SORT_NUM_BUCKET - smem[RADIX_SORT_NUM_BUCKET - 1].i;
    //---------------------------------------------------------------------
    int hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx + (nloop - 1) * bnum);
    //---------------------------------------------------------------------
    int subNum = num - (hidx_gm - head);
    if( subNum < 0 )      subNum = 0;
    subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
    for(int loop = nloop - 1; loop >= 0; loop--){
      //-------------------------------------------------------------------
      __syncthreads();
      //-------------------------------------------------------------------
      if( tidx < RADIX_SORT_NUM_BUCKET ){
	//-----------------------------------------------------------------
	if( loop != (nloop - 1) )
	  smem[STOCK_ARRAY_SUM_HEAD + tidx].i = scan_gm[INDEX2D(nloop * bnum, RADIX_SORT_NUM_BUCKET, bidx + bnum * loop, tidx)];
	//-----------------------------------------------------------------
	/* load global prefix sum */
	smem[tidx].i = scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, tidx, bidx + loop * bnum)];
	//-----------------------------------------------------------------
	/* calculate local prefix sum of each element */
	smem[STOCK_ARRAY_NUM_HEAD + tidx].i = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
	const int tmpVal = (tidx != 0) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 1].i;
	smem[STOCK_ARRAY_SUM_HEAD + tidx].i += tmpVal;
	//-----------------------------------------------------------------
	smem[tidx].i += (smem[STOCK_ARRAY_PTR_HEAD + tidx].i - smem[STOCK_ARRAY_SUM_HEAD + tidx].i);
	//-----------------------------------------------------------------
      }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
      __syncthreads();
      //-------------------------------------------------------------------
      /* store sorted key, idx on global memory */
#pragma unroll
      for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	//-----------------------------------------------------------------
	const RADIX_DATA_TYPE_KEY key = key_sm[ii];
	const int target = (key >> ((numIter - 1) * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK;
	//-----------------------------------------------------------------
	const int dstIdx = head + smem[target].i + ii;
	key_gm[dstIdx] = key;
#ifndef SORT_ONLY
	idx_gm[dstIdx] = idx_sm[ii];
#endif//SORT_ONLY
	//-----------------------------------------------------------------
      }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      /* load temporary stock on the global memory (num_sm, key_dst_sm and idx_dst_sm) */
      if( loop != 0 ){
	//-----------------------------------------------------------------
	hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx + (loop - 1) * bnum);
	subNum = RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT;
	//-----------------------------------------------------------------
#pragma unroll
	for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	  key_sm[tidx + ii * NTHREADS_SORT] = key_tmp_gm[hidx_gm + tidx + ii * NTHREADS_SORT];
#ifndef SORT_ONLY
	  idx_sm[tidx + ii * NTHREADS_SORT] = idx_tmp_gm[hidx_gm + tidx + ii * NTHREADS_SORT];
#endif//SORT_ONLY
	}/* for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){ */
	//-----------------------------------------------------------------
      }/* if( loop != 0 ){ */
      //-------------------------------------------------------------------
    }/* for(int loop = nloop - 1; loop >= 0; loop--){ */
    //---------------------------------------------------------------------
  }/* while( remFul > 0 ){ */
  //-----------------------------------------------------------------------
  if( gidx == 0 ){
    //---------------------------------------------------------------------
    *remLoc_dev = remLoc;
    //---------------------------------------------------------------------
    if( remLoc > remLocMax )
      *fail_dev = remLoc;
    //---------------------------------------------------------------------
  }/* if( gidx == 0 ){ */
  //-----------------------------------------------------------------------
#if 0
  if( gidx == 0 )
    for(int ii = 0; ii < remLoc; ii++)
      printf("ii = %3d: numIter = %d, info.x = %d, info.y = %d\n", ii, numIterLoc[ii], infoLoc[ii].x, infoLoc[ii].y);
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* MSD radix sort within device (2nd half) */
/* LSD radix sort within a block */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_GRID_LSD
(const int numIter, const int rem, const int head, const int tidx, const int lane,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
 const int hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
 RADIX_DATA_TYPE_KEY * RESTRICT key_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_sm,
#ifndef SORT_ONLY
 int * RESTRICT idx_gm, int * RESTRICT idx_sm,
#endif//SORT_ONLY
#   if  RADIX_SORT_CHECK_BITS > 2
 uint4_array * RESTRICT sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
 volatile uint_int * RESTRICT smem
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* a remedy to treat the case when the number of elements is not a multiple of NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD */
#pragma unroll
  for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
    //---------------------------------------------------------------------
#   if  KEY_BITS == 32
    key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffff;
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
    key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffffffffffff;
#endif//KEY_BITS == 64
    //---------------------------------------------------------------------
  }/* for(int ii = subNum + tidx; ii < RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT; ii += NTHREADS_SORT){ */
  //-----------------------------------------------------------------------
  /* load data set from global memory */
#pragma unroll
  for(int ii = tidx; ii < rem; ii += NTHREADS_SORT){
    //---------------------------------------------------------------------
    key_sm[ii] = key_gm[head + ii];
#ifndef SORT_ONLY
    idx_sm[ii] = idx_gm[head + ii];
#endif//SORT_ONLY
    //---------------------------------------------------------------------
  }/* for(int ii = tidx; ii < numIter.x; ii += NTHREADS_SORT){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  /* execute LSD radix sort within a block */
  RADIX_SORT_BLCK
    (numIter, tidx,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
     hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
     key_sm
#ifndef SORT_ONLY
     , idx_sm
#endif//SORT_ONLY
     , smem
#   if  RADIX_SORT_CHECK_BITS > 2
     , sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
     , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
     );
  //-----------------------------------------------------------------------
  /* write back data set to the global memory */
#pragma unroll
  for(int ii = tidx; ii < rem; ii += NTHREADS_SORT){
    //---------------------------------------------------------------------
    key_gm[head + ii] = key_sm[ii];
#ifndef SORT_ONLY
    idx_gm[head + ii] = idx_sm[ii];
#endif//SORT_ONLY
    //---------------------------------------------------------------------
  }/* for(int ii = tidx; ii < numIter.x; ii += NTHREADS_SORT){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* MSD radix sort within a device */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_GRID_MOD
(int * RESTRICT numIterFul, int2 * RESTRICT infoFul, int remFul, int * RESTRICT numIterLoc, int2 * RESTRICT infoLoc, const int remLocMax,
 int * RESTRICT remFulAdd, int * RESTRICT remLocAdd,
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
 int * RESTRICT checkedCounts,
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
 int * RESTRICT fail_dev, const int tidx, const int lane,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
 const int hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
 RADIX_DATA_TYPE_KEY * RESTRICT key_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_tmp_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_sm,
#ifndef SORT_ONLY
 int * RESTRICT idx_gm, int * RESTRICT idx_tmp_gm, int * RESTRICT idx_sm,
#endif//SORT_ONLY
#   if  RADIX_SORT_CHECK_BITS > 2
 uint4_array * RESTRICT sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
 volatile uint_int * RESTRICT smem, int * RESTRICT scan_gm, int * RESTRICT bucket_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul,
 const int gidx, const int bidx, const int bnum,
 volatile int * RESTRICT gsync0, volatile int * RESTRICT gsync1, int * RESTRICT gsync0_loc, int * RESTRICT gsync1_loc
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  int remLoc = 0;
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
  *checkedCounts = 0;
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
  //-----------------------------------------------------------------------
  /* MSD radix sort within device using inter-block synchronization */
  //-----------------------------------------------------------------------
  while( remFul > 0 ){
    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    if( gidx == 0 ){
      *remFulAdd = 0;
      *remLocAdd = 0;
    }/* if( gidx == 0 ){ */
    //---------------------------------------------------------------------
    /* pull data from stack (numIterFul, infoFul) in LIFO (Last In First Out) manner */
    int stackIdx = remFul - 1 - tidx;
#pragma unroll
    for(int ii = tidx; ii < 4 * bnum; ii += NTHREADS_SORT)
      smem[ii].i = 0;
    __syncthreads();
    if( (stackIdx >= 0) && (tidx < bnum) ){
      //-------------------------------------------------------------------
      smem[tidx].i = numIterFul[stackIdx];
      int2 info    =    infoFul[stackIdx];
      smem[tidx +     bnum].i = info.y;/* head */
      smem[tidx + 2 * bnum].i = info.x;/*  rem */
      //-------------------------------------------------------------------
      /* # of blocks required to contain info.x elements with nloop = 1 */
      /* smem[tidx + 3 * bnum].i = BLOCKSIZE(info.x, RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);/\* := bnum_sub *\/ */
      const int tmp = BLOCKSIZE(info.x, RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);/* := bnum_sub */
      smem[tidx + 3 * bnum].i = (tmp < bnum) ? tmp : bnum;
      //-------------------------------------------------------------------
      /* # of blocks may exceed # of threads in a warp */
      //-------------------------------------------------------------------
    }/* if( (stackIdx >= 0) && (tidx < bnum) ){ */
    __syncthreads();
    int head = 0;
    /* int  rem = bnum; */
    int rem = 1;
    stackIdx = smem[3 * bnum].i;
    for(int ii = 1; ii < bnum; ii++){
      /* head = stackIdx; */
      stackIdx += smem[ii + 3 * bnum].i;
      if( stackIdx > bnum ){
	/* if( tidx == 0 ) */
	/*   smem[ii + 3 * bnum].i -= (stackIdx - bnum); */
	rem = ii;
	break;
      }
    }/* for(int ii = 1; ii < bnum; ii++){ */
    /* remFul--; */
    remFul -= rem;

    /* determine bidx_sub */
    bool enabled = false;
    int  idx_sub = 0;
    int disp_sub = 0;
    int bidx_sub = bidx;
    int bnum_sub = bnum;
    int numIter = 0;
    int *gsync0_sub, *gsync1_sub;
    stackIdx = 0;
    head = 0;
    for(int ii = 0; ii < rem; ii++){
      const int tmp = smem[ii + 3 * bnum].i;/* bnum_sub */
      stackIdx += tmp;
      if( bidx < stackIdx ){
	/* configure sub-group */
	enabled = true;
	bidx_sub =     bidx - head;
	bnum_sub = stackIdx - head;
	gsync0_sub = &gsync0_loc[head];
	gsync1_sub = &gsync1_loc[head];
	/* then, target is ii */
	numIter = smem[ii           ].i;
	head    = smem[ii +     bnum].i;
	rem     = smem[ii + 2 * bnum].i;
	break;
      }
      head = stackIdx;
      disp_sub += RADIX_SORT_NUM_BUCKET * tmp;/* nloop is assumed to be unity */
      idx_sub++;
    }/* for(int ii = 0; ii < rem; ii++){ */
    /* const int gidx_sub = tidx + bidx_sub * NTHREADS_SORT; */


    __syncthreads();


    //---------------------------------------------------------------------
    if( enabled ){
      //-------------------------------------------------------------------
      /* evaluate # of iterations to check all elements contained in a current data set */
      /* const int numIter = numIterFul[stackIdx]; */
      /* const int2 info = infoFul[stackIdx]; */
      /* int head = info.y; */
      /* int  rem = info.x; */
      const int2 info = {rem, head};
      const int nloop = BLOCKSIZE(rem, bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
      int totNum, subNum;
      //-------------------------------------------------------------------
      for(int loop = 0; loop < nloop; loop++){
	//-----------------------------------------------------------------
	const int hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx_sub + loop * bnum_sub);
	//-----------------------------------------------------------------
	totNum = ((bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < rem) ? (bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : rem;
	subNum = totNum - bidx_sub * (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
	if( subNum < 0 )	subNum = 0;
	subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
	//-----------------------------------------------------------------
#if 0
	if( tidx == 0 )
	  printf("bidx = %d, bidx_sub = %d, totNum = %d, subNum = %d\n", bidx, bidx_sub, totNum, subNum);
	__syncthreads();
#endif
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------
	/* a remedy to treat the case when the number of elements is not a multiple of NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD */
#pragma unroll
	for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	  //---------------------------------------------------------------
#   if  KEY_BITS == 32
	  key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffff;
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
	  key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffffffffffff;
#endif//KEY_BITS == 64
	  //---------------------------------------------------------------
	}/* for(int ii = subNum + tidx; ii < RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT; ii += NTHREADS_SORT){ */
	//-----------------------------------------------------------------
	/* load key_src, idx_src from global memory */
	/* set key_src and idx_src on shared memory for RADIX_SORT_CHIP */
#pragma unroll
	for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	  //---------------------------------------------------------------
	  const int idx = hidx_gm + ii;
	  //---------------------------------------------------------------
	  key_sm[ii] = key_gm[idx];
#ifndef SORT_ONLY
	  idx_sm[ii] = idx_gm[idx];
#endif//SORT_ONLY
	  //---------------------------------------------------------------
	}/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
	//-----------------------------------------------------------------
	__syncthreads();
	RADIX_SORT_CHIP(numIter - 1, tidx,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
			hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
			key_sm
#ifndef SORT_ONLY
			, idx_sm
#endif//SORT_ONLY
			, smem
#   if  RADIX_SORT_CHECK_BITS > 2
			, sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
			, m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
			);
	//-----------------------------------------------------------------
	/* remove the number of ghost values */
	/* target == RADIX_SORT_CHECK_MASK (= RADIX_SORT_NUM_BUCKET - 1) */
	/* # of ghosts is RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum */
	if( tidx == RADIX_SORT_CHECK_MASK )
	  smem[STOCK_ARRAY_SUM_HEAD + tidx].i -= RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum;
	//-----------------------------------------------------------------
	/* put returned values of num_sm on global memory to calculate global prefix sum */
	if( tidx < RADIX_SORT_NUM_BUCKET )
	  scanFul[disp_sub + INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum_sub, tidx, bidx_sub + bnum_sub * loop)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
	//-----------------------------------------------------------------
	/* if nloop > 1 && loop < nloop - 1; save num_sm, key_dst and idx_dst on the global memory */
	if( loop < (nloop - 1) ){
	  //---------------------------------------------------------------
	  /* save num_sm on the global memory */
	  if( tidx < RADIX_SORT_NUM_BUCKET )
	    scan_gm[disp_sub + INDEX2D(nloop * bnum_sub, RADIX_SORT_NUM_BUCKET, bidx_sub + bnum_sub * loop, tidx)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i;
	  //---------------------------------------------------------------
	  /* save key_dst and idx_dst on the global memory */
#pragma unroll
	  for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	    key_tmp_gm[hidx_gm + ii] = key_sm[ii];
#ifndef SORT_ONLY
	    idx_tmp_gm[hidx_gm + ii] = idx_sm[ii];
#endif//SORT_ONLY
	  }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
	  //---------------------------------------------------------------
	}/* if( loop < (nloop - 1) ){ */
	//-----------------------------------------------------------------
	rem -= totNum;
	//-----------------------------------------------------------------
      }/* for(int loop = 0; loop < nloop; loop++){ */
      //-------------------------------------------------------------------
      /* calculate global prefix sum within a device */
      /* if # of elements is 128M --> # of elements per block is 128M / (4 * Ttot) ~ 512k >> 1024. */
      /* Therefore, simple evaluation of prefix sum within a block is not sufficient */
      int tmp_sum_head;
      if( tidx < RADIX_SORT_NUM_BUCKET )
	tmp_sum_head = smem[STOCK_ARRAY_SUM_HEAD + tidx].i;
      globalSync(tidx, bidx_sub, bnum_sub, gsync0_sub, gsync1_sub);
      for(int target = bidx_sub; target < RADIX_SORT_NUM_BUCKET; target += bnum_sub){
	//-----------------------------------------------------------------
	int subTail = 0;
	const int niter = BLOCKSIZE(bnum_sub * nloop, NTHREADS_SORT);
	//-----------------------------------------------------------------
	for(int iter = 0; iter < niter; iter++){
	  //---------------------------------------------------------------
	  const int sidx = tidx + iter * NTHREADS_SORT;
	  //---------------------------------------------------------------
	  /* load local prefix sum */
	  smem[tidx].i = (sidx < bnum_sub * nloop) ? scanFul[disp_sub + INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum_sub, target, sidx)] : 0;
	  /* prefix sum within a block */
	  RADIX_PREFIX_SUM_GRID(tidx, lane, smem
#ifdef  USE_MASK_INSTEAD_OF_IF
				, m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
				);
	  /* write back global memory */
	  if( sidx < bnum_sub * nloop )
	    scanFul[disp_sub + INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum_sub, target, sidx)] = subTail + smem[tidx].i;
	  subTail += smem[NTHREADS_SORT - 1].i;
	  __syncthreads();
	  //---------------------------------------------------------------
	}/* for(int iter = 0; iter < niter; iter++){ */
	//-----------------------------------------------------------------
	if( tidx == 0 ){
	  tail_gm  [idx_sub * RADIX_SORT_NUM_BUCKET + target] = subTail;
	  /* 1 means that further global sort is required */
	  bucket_gm[idx_sub * RADIX_SORT_NUM_BUCKET + target] = (subTail < (NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD)) ? 0 : 1;
	}/* if( tidx == 0 ){ */
	//-----------------------------------------------------------------
      }/* for(int target = bidx; target < RADIX_SORT_NUM_BUCKET; target += bnum){ */
      if( tidx < RADIX_SORT_NUM_BUCKET )
	smem[STOCK_ARRAY_SUM_HEAD + tidx].i = tmp_sum_head;
      globalSync(tidx, bidx_sub, bnum_sub, gsync0_sub, gsync1_sub);
      //-------------------------------------------------------------------
      if( tidx < RADIX_SORT_NUM_BUCKET ){
	//-----------------------------------------------------------------
	/* get partition among sub-groups */
	subNum = tail_gm[idx_sub * RADIX_SORT_NUM_BUCKET + tidx];
	/* 1 means that further global sort is required */
	const int bucket_tmp = bucket_gm[idx_sub * RADIX_SORT_NUM_BUCKET + tidx];
	const int sidx = STOCK_ARRAY_PTR_HEAD + tidx;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	/* load index */
	int v0 = subNum;
	int v1 = bucket_tmp;
	int t0, t1;
	/* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
	t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET) & m0;      t1 = __shfl_up(v1, 1, RADIX_SORT_NUM_BUCKET) & m0;      v0 += t0;	v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >=  4
	t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET) & m1;      t1 = __shfl_up(v1, 2, RADIX_SORT_NUM_BUCKET) & m1;      v0 += t0;	v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >=  8
	t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET) & m2;      t1 = __shfl_up(v1, 4, RADIX_SORT_NUM_BUCKET) & m2;      v0 += t0;	v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >= 16
	t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET) & m3;      t1 = __shfl_up(v1, 8, RADIX_SORT_NUM_BUCKET) & m3;      v0 += t0;	v1 += t1;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
	t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 1, RADIX_SORT_NUM_BUCKET);      if( tidx >= 1 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >=  4
	t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 2, RADIX_SORT_NUM_BUCKET);      if( tidx >= 2 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >=  8
	t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 4, RADIX_SORT_NUM_BUCKET);      if( tidx >= 4 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >= 16
	t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 8, RADIX_SORT_NUM_BUCKET);      if( tidx >= 8 ){	v0 += t0;	v1 += t1;      }
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
        /* return calculated inclusive prefix sum */
	smem[sidx].i = v0;
	smem[tidx].i = v1;
#else///USE_WARP_SHUFFLE_FUNC_SORT
	smem[sidx].i = subNum;
	smem[tidx].i = bucket_tmp;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
	smem[tidx].i += smem[tidx - 1].i & m0;	smem[sidx].i += smem[sidx - 1].i & m0;
#   if  RADIX_SORT_NUM_BUCKET >=  4
	smem[tidx].i += smem[tidx - 2].i & m1;	smem[sidx].i += smem[sidx - 2].i & m1;
#   if  RADIX_SORT_NUM_BUCKET >=  8
	smem[tidx].i += smem[tidx - 4].i & m2;	smem[sidx].i += smem[sidx - 4].i & m2;
#   if  RADIX_SORT_NUM_BUCKET >= 16
	smem[tidx].i += smem[tidx - 8].i & m3;	smem[sidx].i += smem[sidx - 8].i & m3;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
	if( tidx >= 1 ){	smem[tidx].i += smem[tidx - 1].i;	smem[sidx].i += smem[sidx - 1].i;      }
#   if  RADIX_SORT_NUM_BUCKET >=  4
	if( tidx >= 2 ){	smem[tidx].i += smem[tidx - 2].i;	smem[sidx].i += smem[sidx - 2].i;      }
#   if  RADIX_SORT_NUM_BUCKET >=  8
	if( tidx >= 4 ){	smem[tidx].i += smem[tidx - 4].i;	smem[sidx].i += smem[sidx - 4].i;      }
#   if  RADIX_SORT_NUM_BUCKET >= 16
	if( tidx >= 8 ){	smem[tidx].i += smem[tidx - 8].i;	smem[sidx].i += smem[sidx - 8].i;      }
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#endif//USE_WARP_SHUFFLE_FUNC_SORT
        //-----------------------------------------------------------------
	smem[STOCK_ARRAY_DST_HEAD + tidx].i = bucket_tmp;
	smem[STOCK_ARRAY_IDX_HEAD + tidx].i = bucket_tmp ? (smem[tidx].i - bucket_tmp) : (tidx - smem[tidx].i);
	//-----------------------------------------------------------------
	smem[STOCK_ARRAY_PTR_NUM  + tidx].i  = subNum;
	smem[STOCK_ARRAY_PTR_HEAD + tidx].i -= subNum;
	//-----------------------------------------------------------------
	if( bidx_sub == 0 ){
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	  int remFulHead, remLocHead;
#else///USE_WARP_SHUFFLE_FUNC_SORT
	  const int stock = smem[tidx].i;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
	  if( tidx == 0 ){
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	    remFulHead = remFul + atomicAdd(remFulAdd,                         smem[RADIX_SORT_NUM_BUCKET - 1].i);
	    remLocHead = remLoc + atomicAdd(remLocAdd, RADIX_SORT_NUM_BUCKET - smem[RADIX_SORT_NUM_BUCKET - 1].i);
#else///USE_WARP_SHUFFLE_FUNC_SORT
	    smem[0].i  = remFul + atomicAdd(remFulAdd,                         smem[RADIX_SORT_NUM_BUCKET - 1].i);
	    smem[1].i  = remLoc + atomicAdd(remLocAdd, RADIX_SORT_NUM_BUCKET - smem[RADIX_SORT_NUM_BUCKET - 1].i);
#endif//USE_WARP_SHUFFLE_FUNC_SORT
	  }
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	  remFulHead = __shfl(remFulHead, 0, RADIX_SORT_NUM_BUCKET);
	  remLocHead = __shfl(remLocHead, 0, RADIX_SORT_NUM_BUCKET);
#else///USE_WARP_SHUFFLE_FUNC_SORT
	  int remFulHead = smem[0].i;
	  int remLocHead = smem[1].i;
	  smem[tidx].i = stock;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
	  int2 subInfo = {subNum, head + smem[STOCK_ARRAY_PTR_HEAD + tidx].i};
	  if( smem[STOCK_ARRAY_DST_HEAD + tidx].i ){
	    infoFul   [remFulHead + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = subInfo;
	    numIterFul[remFulHead + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = numIter - 1;
	  }
	  else{
	    infoLoc   [remLocHead + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = subInfo;
	    numIterLoc[remLocHead + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = numIter - 1;
	  }
	}/* if( bidx == 0 ){ */
	//-----------------------------------------------------------------
      }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
      //-------------------------------------------------------------------
      int hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx_sub + (nloop - 1) * bnum_sub);
      //-------------------------------------------------------------------
#ifdef  ENABLE_EARLY_EXIT_GRID
      int stock;
      if( tidx < RADIX_SORT_NUM_BUCKET ){
	int flag = smem[tidx].i != info.x;/* 0 means sort is not required */
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	if( tidx == 0 )	  stock = smem[0].i;
#   if  RADIX_SORT_NUM_BUCKET >=  2
	int temp;
	temp = __shfl_up(flag, 1, RADIX_SORT_NUM_BUCKET);	flag *= temp;
#   if  RADIX_SORT_NUM_BUCKET >=  4
	temp = __shfl_up(flag, 2, RADIX_SORT_NUM_BUCKET);	flag *= temp;
#   if  RADIX_SORT_NUM_BUCKET >=  8
	temp = __shfl_up(flag, 4, RADIX_SORT_NUM_BUCKET);	flag *= temp;
#   if  RADIX_SORT_NUM_BUCKET == 16
	temp = __shfl_up(flag, 8, RADIX_SORT_NUM_BUCKET);	flag *= temp;
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET == 16
	if( tidx == 0 )
	  smem[0].i = flag;
#else///USE_WARP_SHUFFLE_FUNC_SORT
	stock = smem[tidx].i;
	smem[tidx].i = flag;
#   if  RADIX_SORT_NUM_BUCKET >=  2
	int temp;
	temp = smem[tidx ^ 1];	flag *= temp;	smem[tidx].i = flag;
#   if  RADIX_SORT_NUM_BUCKET >=  4
	temp = smem[tidx ^ 2];	flag *= temp;	smem[tidx].i = flag;
#   if  RADIX_SORT_NUM_BUCKET >=  8
	temp = smem[tidx ^ 4];	flag *= temp;	smem[tidx].i = flag;
#   if  RADIX_SORT_NUM_BUCKET == 16
	temp = smem[tidx ^ 8];	flag *= temp;	smem[tidx].i = flag;
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET == 16
#endif//USE_WARP_SHUFFLE_FUNC_SORT
      }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
      __syncthreads();
      const int exec = smem[0].i;
      if( exec )
#endif//ENABLE_EARLY_EXIT_GRID
	{
	  //---------------------------------------------------------------
#ifdef  ENABLE_EARLY_EXIT_GRID
	  __syncthreads();
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	  if( tidx == 0 )
	    smem[0].i = stock;
#else///USE_WARP_SHUFFLE_FUNC_SORT
	  if( tidx < RADIX_SORT_NUM_BUCKET )
	    smem[tidx].i = stock;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
	  __syncthreads();
#endif//ENABLE_EARLY_EXIT_GRID
	  //---------------------------------------------------------------
	  subNum = info.x - (bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) * (nloop - 1);
	  totNum = ((bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
	  subNum = totNum - bidx_sub * (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
	  if( subNum < 0 )	    subNum = 0;
	  subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
	  for(int loop = nloop - 1; loop >= 0; loop--){
	    //-------------------------------------------------------------
	    __syncthreads();
	    if( tidx < RADIX_SORT_NUM_BUCKET ){
	      //-----------------------------------------------------------
	      if( loop != (nloop - 1) )
		smem[STOCK_ARRAY_SUM_HEAD + tidx].i = scan_gm[INDEX2D(nloop * bnum, RADIX_SORT_NUM_BUCKET, bidx + bnum * loop, tidx)];
	      //-----------------------------------------------------------
	      /* load global prefix sum */
	      smem[tidx].i = scanFul[disp_sub + INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum_sub, tidx, bidx_sub + loop * bnum_sub)];
	      //-----------------------------------------------------------
	      /* calculate local prefix sum of each element */
	      smem[STOCK_ARRAY_NUM_HEAD + tidx].i = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
	      const int tmpVal = (tidx != 0) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 1].i;
	      smem[STOCK_ARRAY_SUM_HEAD + tidx].i += tmpVal;
	      //-----------------------------------------------------------
	      smem[tidx].i += (smem[STOCK_ARRAY_PTR_HEAD + tidx].i - smem[STOCK_ARRAY_SUM_HEAD + tidx].i);
	      //-----------------------------------------------------------
	    }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
	    __syncthreads();
	    //-------------------------------------------------------------
	    /* store sorted key, idx on global memory */
	    for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	      //-----------------------------------------------------------
	      const RADIX_DATA_TYPE_KEY key = key_sm[ii];
	      const int target = (key >> ((numIter - 1) * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK;
	      //-----------------------------------------------------------
	      const int dstIdx = head + smem[target].i + ii;
	      key_gm[dstIdx] = key;
#ifndef SORT_ONLY
	      idx_gm[dstIdx] = idx_sm[ii];
#endif//SORT_ONLY
	      //-----------------------------------------------------------
	    }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
	    //-------------------------------------------------------------
	    /* load temporary stock on the global memory (num_sm, key_dst_sm and idx_dst_sm) */
	    if( loop != 0 ){
	      //-----------------------------------------------------------
	      hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx_sub + (loop - 1) * bnum_sub);
	      subNum = RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT;
	      //-----------------------------------------------------------
#pragma unroll
	      for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
		key_sm[tidx + ii * NTHREADS_SORT] = key_tmp_gm[hidx_gm + tidx + ii * NTHREADS_SORT];
#ifndef SORT_ONLY
		idx_sm[tidx + ii * NTHREADS_SORT] = idx_tmp_gm[hidx_gm + tidx + ii * NTHREADS_SORT];
#endif//SORT_ONLY
	      }/* for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){ */
	      //-----------------------------------------------------------
	    }/* if( loop != 0 ){ */
	    //-------------------------------------------------------------
	  }/* for(int loop = nloop - 1; loop >= 0; loop--){ */
	  //---------------------------------------------------------------
	}
      //-------------------------------------------------------------------
    }/* if( enabled ){ */
    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    /* load remFul and remLoc from global memory */
    if( tidx == 0 ){
      smem[0].i = remFul + (*remFulAdd);
      smem[1].i = remLoc + (*remLocAdd);
    }
    __syncthreads();
    remFul = smem[0].i;
    remLoc = smem[1].i;
    //---------------------------------------------------------------------
  }/* while( remFul > 0 ){ */
  //-----------------------------------------------------------------------
  if( (gidx == 0) && (remLoc > remLocMax) )
    *fail_dev = remLoc;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* LSD radix sort within a block */
  //-----------------------------------------------------------------------
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
  if( tidx == 0 )
    smem[0].i = atomicAdd(checkedCounts, 1);
  __syncthreads();
  int target = smem[0].i;
  while( target < remLoc )
#else///ATOMIC_BASED_JOB_ASSIGNMENT
  for(int target = bidx; target < remLoc; target += bnum)
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
    {
      //-------------------------------------------------------------------
      /* pull data from queue (numIterLoc, infoLoc) in FIFO (First In First Out) manner */
      if( tidx == 0 ){
	int2 info = infoLoc[target];
	smem[0].i = numIterLoc[target];
	smem[1].i = info.x;
	smem[2].i = info.y;
      }
      __syncthreads();
      const int numIter = smem[0].i;
      const int  rem = smem[1].i;
      const int head = smem[2].i;
      //-------------------------------------------------------------------
      /* initialize the shared memory */
      /* a remedy to treat the case when the number of elements is not a multiple of NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD */
#pragma unroll
      for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	//-----------------------------------------------------------------
#   if  KEY_BITS == 32
	key_sm[tidx + ii * NTHREADS_SORT] = 0xffffffff;
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
	key_sm[tidx + ii * NTHREADS_SORT] = 0xffffffffffffffff;
#endif//KEY_BITS == 64
        //-----------------------------------------------------------------
      }/* for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){ */
      //-------------------------------------------------------------------
      /* load data set from global memory */
#pragma unroll
      for(int ii = tidx; ii < rem; ii += NTHREADS_SORT){
	//-----------------------------------------------------------------
	key_sm[ii] = key_gm[head + ii];
#ifndef SORT_ONLY
	idx_sm[ii] = idx_gm[head + ii];
#endif//SORT_ONLY
        //-----------------------------------------------------------------
      }/* for(int ii = tidx; ii < numIter.x; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      __syncthreads();
      //-------------------------------------------------------------------
      /* execute LSD radix sort within a block */
      RADIX_SORT_BLCK
	(numIter, tidx,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
	 hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
	 key_sm
#ifndef SORT_ONLY
	 , idx_sm
#endif//SORT_ONLY
	 , smem
#   if  RADIX_SORT_CHECK_BITS > 2
	 , sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
	 , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
	 );
      //-------------------------------------------------------------------
      /* write back data set to the global memory */
#pragma unroll
      for(int ii = tidx; ii < rem; ii += NTHREADS_SORT){
	//-----------------------------------------------------------------
	key_gm[head + ii] = key_sm[ii];
#ifndef SORT_ONLY
	idx_gm[head + ii] = idx_sm[ii];
#endif//SORT_ONLY
        //-----------------------------------------------------------------
      }/* for(int ii = tidx; ii < numIter.x; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
#ifdef  ATOMIC_BASED_JOB_ASSIGNMENT
      if( tidx == 0 )
	smem[0].i = atomicAdd(checkedCounts, 1);
      __syncthreads();
      target = smem[0].i;
#endif//ATOMIC_BASED_JOB_ASSIGNMENT
      //-------------------------------------------------------------------
    }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* MSD radix sort within a device */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_GRID_GET_PARTITION
(const int lane, const int tidx, const int bidx, const int bnum,
 const int remFul, int * RESTRICT scaned, int * RESTRICT numIterFul, int2 * RESTRICT infoFul,
 bool * RESTRICT enabled, int * RESTRICT bidx_sub, int * RESTRICT bnum_sub, int * RESTRICT numIter, int * RESTRICT head, int * RESTRICT rem, int * RESTRICT disp_sub, int * RESTRICT disp_ful, int * RESTRICT sidxLoc, int * RESTRICT sidxFul,
 volatile uint_int * RESTRICT smem
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* 4 * bnum must be equal or less than NTHREADS_SORT */
  //-----------------------------------------------------------------------
  /* pull data from stack (numIterFul, infoFul) in LIFO (Last In First Out) manner */
  int sidx = remFul - 1 - ((*scaned) + tidx);
  __syncthreads();
  if( tidx < bnum ){
    //---------------------------------------------------------------------
    if( sidx >= 0 ){
      //-------------------------------------------------------------------
      smem[tidx].i = numIterFul[sidx];
      int2 info    =    infoFul[sidx];
      smem[tidx +     bnum].i = info.y;/* head */
      smem[tidx + 2 * bnum].i = info.x;/*  rem */
      //-------------------------------------------------------------------
      /* # of blocks required to contain info.x elements with nloop = 1 */
      const int tmp = BLOCKSIZE(info.x, RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);/* := bnum_sub */
      smem[tidx + 3 * bnum].i = (tmp < bnum) ? tmp : bnum;
      //-------------------------------------------------------------------
    }/* if( sidx >= 0 ){ */
    //---------------------------------------------------------------------
    else{
      //-------------------------------------------------------------------
      smem[tidx           ].i = 0;
      smem[tidx +     bnum].i = 0;
      smem[tidx + 2 * bnum].i = 0;
      smem[tidx + 3 * bnum].i = 0;
      //-------------------------------------------------------------------
    }/* else{ */
    //---------------------------------------------------------------------
  }/* if( tidx < bnum ){ */
  __syncthreads();
  //-----------------------------------------------------------------------
  /* 1. prefix sum within warp */
  int val = 0;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  int tmp = 0;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
  /* load index */
  if( tidx < bnum )
    val = smem[tidx + 3 * bnum].i;
  /* calculate inclusive prefix sum */
  const int groupSize = (bnum < warpSize) ? bnum : warpSize;
#ifdef  USE_MASK_INSTEAD_OF_IF
  if(         bnum >=  2 ){     tmp = __shfl_up(val,  1, groupSize) & m0;       val += tmp;
    if(       bnum >=  4 ){     tmp = __shfl_up(val,  2, groupSize) & m1;       val += tmp;
      if(     bnum >=  8 ){	tmp = __shfl_up(val,  4, groupSize) & m2;	val += tmp;
  	if(   bnum >= 16 ){	tmp = __shfl_up(val,  8, groupSize) & m3;	val += tmp;
  	  if( bnum >= 32 ){	tmp = __shfl_up(val, 16, groupSize) & m4;	val += tmp;
  	  }}}}}
#else///USE_MASK_INSTEAD_OF_IF
  if(         bnum >=  2 ){     tmp = __shfl_up(val,  1, groupSize);    if( lane >=  1 )        val += tmp;
    if(       bnum >=  4 ){     tmp = __shfl_up(val,  2, groupSize);    if( lane >=  2 )	val += tmp;
      if(     bnum >=  8 ){	tmp = __shfl_up(val,  4, groupSize);	if( lane >=  4 )	val += tmp;
  	if(   bnum >= 16 ){	tmp = __shfl_up(val,  8, groupSize);	if( lane >=  8 )	val += tmp;
  	  if( bnum >= 32 ){	tmp = __shfl_up(val, 16, groupSize);	if( lane >= 16 )	val += tmp;
  	  }}}}}
#endif//USE_MASK_INSTEAD_OF_IF
  /* return calculated inclusive prefix sum */
  if( tidx < bnum )
    smem[tidx + 3 * bnum].i = val;
  //-----------------------------------------------------------------------
#else///USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------
  /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
  if(         bnum >=  2 ){     smem[tidx + 3 * bnum].i += smem[tidx -  1 + 3 * bnum].i & m0;
    if(       bnum >=  4 ){     smem[tidx + 3 * bnum].i += smem[tidx -  2 + 3 * bnum].i & m1;
      if(     bnum >=  8 ){	smem[tidx + 3 * bnum].i += smem[tidx -  4 + 3 * bnum].i & m2;
	if(   bnum >= 16 ){	smem[tidx + 3 * bnum].i += smem[tidx -  8 + 3 * bnum].i & m3;
	  if( bnum >= 32 ){	smem[tidx + 3 * bnum].i += smem[tidx - 16 + 3 * bnum].i & m4;
	  }}}}}
#else///USE_MASK_INSTEAD_OF_IF
  if(         bnum >=  2 ){     if( lane >=  1 )        smem[tidx + 3 * bnum].i += smem[tidx -  1 + 3 * bnum].i;
    if(       bnum >=  4 ){     if( lane >=  2 )	smem[tidx + 3 * bnum].i += smem[tidx -  2 + 3 * bnum].i;
      if(     bnum >=  8 ){	if( lane >=  4 )	smem[tidx + 3 * bnum].i += smem[tidx -  4 + 3 * bnum].i;
	if(   bnum >= 16 ){	if( lane >=  8 )	smem[tidx + 3 * bnum].i += smem[tidx -  8 + 3 * bnum].i;
	  if( bnum >= 32 ){	if( lane >= 16 )	smem[tidx + 3 * bnum].i += smem[tidx - 16 + 3 * bnum].i;
	  }}}}}
#endif//USE_MASK_INSTEAD_OF_IF
  //-----------------------------------------------------------------------
#endif///USE_WARP_SHUFFLE_FUNC_SORT
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 2. prefix sum about tail of each warp */
  //-----------------------------------------------------------------------
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
  int scan = val;
#else///USE_WARP_SHUFFLE_FUNC_SORT
  int scan = 0;
  if( tidx < bnum )
    scan = smem[tidx + 3 * bnum].i;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
  __syncthreads();
  /* warpSize = 32 = 2^5 */
  if( tidx < (bnum >> 5) ){
    //---------------------------------------------------------------------
    val = smem[tidx * warpSize + warpSize - 1 + 3 * bnum].i;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    const int groupSize = bnum >> 5;
#ifdef  USE_MASK_INSTEAD_OF_IF
    if(         bnum >=   64 ){ tmp = __shfl_up(val,  1, groupSize) & m0;       val += tmp;
      if(       bnum >=  128 ){	tmp = __shfl_up(val,  2, groupSize) & m1;	val += tmp;
	if(     bnum >=  256 ){	tmp = __shfl_up(val,  4, groupSize) & m2;	val += tmp;
	  if(   bnum >=  512 ){	tmp = __shfl_up(val,  8, groupSize) & m3;	val += tmp;
	    if( bnum == 1024 ){	tmp = __shfl_up(val, 16, groupSize) & m4;	val += tmp;
	    }}}}}
#else///USE_MASK_INSTEAD_OF_IF
    if(         bnum >=   64 ){ tmp = __shfl_up(val,  1, groupSize);    if( lane >=  1 )	val += tmp;
      if(       bnum >=  128 ){	tmp = __shfl_up(val,  2, groupSize);	if( lane >=  2 )	val += tmp;
	if(     bnum >=  256 ){	tmp = __shfl_up(val,  4, groupSize);	if( lane >=  4 )	val += tmp;
	  if(   bnum >=  512 ){	tmp = __shfl_up(val,  8, groupSize);	if( lane >=  8 )	val += tmp;
	    if( bnum == 1024 ){	tmp = __shfl_up(val, 16, groupSize);	if( lane >= 16 )	val += tmp;
	    }}}}}
#endif//USE_MASK_INSTEAD_OF_IF
#else///USE_WARP_SHUFFLE_FUNC_SORT
    smem[tidx + 3 * bnum].i = val;
#ifdef  USE_MASK_INSTEAD_OF_IF
    if(         bnum >=   64 ){ smem[tidx + 3 * bnum].i += smem[tidx -  1 + 3 * bnum].i & m0;
      if(       bnum >=  128 ){	smem[tidx + 3 * bnum].i += smem[tidx -  2 + 3 * bnum].i & m1;
	if(     bnum >=  256 ){	smem[tidx + 3 * bnum].i += smem[tidx -  4 + 3 * bnum].i & m2;
	  if(   bnum >=  512 ){	smem[tidx + 3 * bnum].i += smem[tidx -  8 + 3 * bnum].i & m3;
	    if( bnum == 1024 ){	smem[tidx + 3 * bnum].i += smem[tidx - 16 + 3 * bnum].i & m4;
	    }}}}}
#else///USE_MASK_INSTEAD_OF_IF
    if(         bnum >=   64 ){ if( lane >=  1 )	smem[tidx + 3 * bnum].i += smem[tidx -  1 + 3 * bnum].i;
      if(       bnum >=  128 ){	if( lane >=  2 )	smem[tidx + 3 * bnum].i += smem[tidx -  2 + 3 * bnum].i;
	if(     bnum >=  256 ){	if( lane >=  4 )	smem[tidx + 3 * bnum].i += smem[tidx -  4 + 3 * bnum].i;
	  if(   bnum >=  512 ){	if( lane >=  8 )	smem[tidx + 3 * bnum].i += smem[tidx -  8 + 3 * bnum].i;
	    if( bnum == 1024 ){	if( lane >= 16 )	smem[tidx + 3 * bnum].i += smem[tidx - 16 + 3 * bnum].i;
	    }}}}}
#endif//USE_MASK_INSTEAD_OF_IF
    val = smem[tidx + 3 * bnum].i;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
    //---------------------------------------------------------------------
    smem[tidx + 3 * bnum].i = val;
    //---------------------------------------------------------------------
  }/* if( tidx < (NTHREADS_SORT >> 5) ){ */
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 3. prefix sum within a block */
  //-----------------------------------------------------------------------
  /* warpSize = 32 = 2^5 */
  if( (tidx < bnum) && (tidx >= warpSize) )
    scan += smem[(tidx >> 5) - 1 + 3 * bnum].i;
  __syncthreads();
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* 4. upload calculate prefix sum */
  //-----------------------------------------------------------------------
  if( tidx < bnum )
    smem[tidx + 3 * bnum].i = scan;
  __syncthreads();
  //-----------------------------------------------------------------------
  int scanNum = 1;
  for(int ii = 1; ii < bnum; ii++){
    if( smem[ii + 3 * bnum].i > bnum ){
      scanNum = ii;
      break;
    }/* if( smem[ii + 3 * bnum].i > bnum ){ */
  }/* for(int ii = 1; ii < bnum; ii++){ */
  if( scanNum > (RADIX_SORT_BOX_NUM_ALLOCATED - (*scaned)) )
    scanNum = (RADIX_SORT_BOX_NUM_ALLOCATED - (*scaned));
  __syncthreads();
  /* determine sidx */
  int scanHead = 0;
  bool determined = false;
  for(int ii = 0; ii < scanNum; ii++){
    sidx = smem[ii + 3 * bnum].i;
    const int bnum_tmp = sidx - scanHead;
    if( !determined && (bidx < sidx) ){
      determined = true;
      /* configure sub-group */
      *enabled  = true;
      *bidx_sub = bidx - scanHead;
      *bnum_sub = bnum_tmp;
      /* then, target is ii */
      *numIter = smem[ii           ].i;
      *head    = smem[ii +     bnum].i;
      *rem     = smem[ii + 2 * bnum].i;
      *disp_sub = *disp_ful;
    }/* if( !determined && (bidx < sidx) ){ */
    scanHead = sidx;
    if( !determined )
      *sidxLoc += 1;
    *disp_ful += RADIX_SORT_NUM_BUCKET * bnum_tmp * BLOCKSIZE(smem[ii + 2 * bnum].i, bnum_tmp * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
  }/* for(int ii = 0; ii < scanNum; ii++){ */
  *scaned  += scanNum;
  *sidxFul += scanNum;/* necessary to set idx_sub in the next step */
  __syncthreads();
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_GRID_MOD_SCATTER
(int * RESTRICT numIterFul, int2 * RESTRICT infoFul, int remFul, int * RESTRICT snum,
 int2 * RESTRICT disp,
 const int gidx, const int bidx, const int bnum, const int tidx, const int lane,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
 const int hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
 RADIX_DATA_TYPE_KEY * RESTRICT key_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_tmp_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_sm,
#ifndef SORT_ONLY
 int * RESTRICT idx_gm, int * RESTRICT idx_tmp_gm, int * RESTRICT idx_sm,
#endif//SORT_ONLY
#   if  RADIX_SORT_CHECK_BITS > 2
 uint4_array * RESTRICT sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
 volatile uint_int * RESTRICT smem, int * RESTRICT scan_gm, int * RESTRICT scanFul
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* MSD radix sort within a device (1st scan phase) */
  int scaned = 0;
  int dispGlobal = 0;
  int sidxGlobal = 0;
  //-----------------------------------------------------------------------
  const int remMax = (remFul > RADIX_SORT_BOX_NUM_ALLOCATED) ? (RADIX_SORT_BOX_NUM_ALLOCATED) : (remFul);
  //-----------------------------------------------------------------------
  while( scaned < remMax ){
    //---------------------------------------------------------------------
    bool enabled = false;
    int bidx_sub = bidx;
    int bnum_sub = bnum;
    int numIter  = 0;
    int head = 0;
    int rem = 0;
    int disp_sub = dispGlobal;
    int sidx     = sidxGlobal;
    RADIX_SORT_GRID_GET_PARTITION(lane, tidx, bidx, bnum, remFul, &scaned, numIterFul, infoFul, &enabled, &bidx_sub, &bnum_sub, &numIter, &head, &rem, &disp_sub, &dispGlobal, &sidx, &sidxGlobal, smem
#ifdef  USE_MASK_INSTEAD_OF_IF
				  , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
				  );
    //---------------------------------------------------------------------
    if( enabled ){
      //-------------------------------------------------------------------
      /* evaluate # of iterations to check all elements contained in a current data set */
      /* const int2 info = {rem, head}; */
      const int nloop = BLOCKSIZE(rem, bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
      int totNum, subNum;
      //-------------------------------------------------------------------
      if( (bidx_sub + tidx) == 0 ){
	/* this condition corresponds to bidx_sub == 0 && tidx == 0 */
	const int2 tmp = {disp_sub, bnum_sub * nloop};
	disp[sidx] = tmp;
      }/* if( (bidx_sub + tidx) == 0 ){ */
      //-------------------------------------------------------------------
      for(int loop = 0; loop < nloop; loop++){
	//-----------------------------------------------------------------
	const int hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx_sub + loop * bnum_sub);
	//-----------------------------------------------------------------
	totNum = ((bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < rem) ? (bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : rem;
	subNum = totNum - bidx_sub * (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
	if( subNum < 0 )	  subNum = 0;
	subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------
	/* a remedy to treat the case when the number of elements is not a multiple of NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD */
#pragma unroll
	for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	  //---------------------------------------------------------------
#   if  KEY_BITS == 32
	  key_sm[tidx + ii * NTHREADS_SORT] = 0xffffffff;
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
	  key_sm[tidx + ii * NTHREADS_SORT] = 0xffffffffffffffff;
#endif//KEY_BITS == 64
	  //---------------------------------------------------------------
	}/* for(int ii = subNum + tidx; ii < RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT; ii += NTHREADS_SORT){ */
	//-----------------------------------------------------------------
	/* load key_src, idx_src from global memory */
	/* set key_src and idx_src on shared memory for RADIX_SORT_CHIP */
#pragma unroll
	for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	  //---------------------------------------------------------------
	  const int idx = hidx_gm + ii;
	  //---------------------------------------------------------------
	  key_sm[ii] = key_gm[idx];
#ifndef SORT_ONLY
	  idx_sm[ii] = idx_gm[idx];
#endif//SORT_ONLY
	  //---------------------------------------------------------------
	}/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
	//-----------------------------------------------------------------
	__syncthreads();
	RADIX_SORT_CHIP(numIter - 1, tidx,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
			hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
			key_sm
#ifndef SORT_ONLY
			, idx_sm
#endif//SORT_ONLY
			, smem
#   if  RADIX_SORT_CHECK_BITS > 2
			, sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
			, m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
			);
	//-----------------------------------------------------------------
	/* remove the number of ghost values */
	/* target == RADIX_SORT_CHECK_MASK (= RADIX_SORT_NUM_BUCKET - 1) */
	/* # of ghosts is RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum */
	if( tidx == RADIX_SORT_CHECK_MASK )
	  smem[STOCK_ARRAY_SUM_HEAD + tidx].i -= RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum;
	//-----------------------------------------------------------------
	/* put returned values of num_sm on global memory to calculate global prefix sum */
	if( tidx < RADIX_SORT_NUM_BUCKET ){
	  /* save num_sm on the global memory */
	  scanFul[disp_sub + INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum_sub, tidx, bidx_sub + bnum_sub * loop)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
	  scan_gm[disp_sub + INDEX2D(nloop * bnum_sub, RADIX_SORT_NUM_BUCKET, bidx_sub + bnum_sub * loop, tidx)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i;
	}/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
	//-----------------------------------------------------------------
	/* save key_dst and idx_dst on the global memory */
#pragma unroll
	for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	  key_tmp_gm[hidx_gm + ii] = key_sm[ii];
#ifndef SORT_ONLY
	  idx_tmp_gm[hidx_gm + ii] = idx_sm[ii];
#endif//SORT_ONLY
	}/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
	//-----------------------------------------------------------------
	rem -= totNum;
	//-----------------------------------------------------------------
      }/* for(int loop = 0; loop < nloop; loop++){ */
      //-------------------------------------------------------------------
    }/* if( enabled ){ */
    //---------------------------------------------------------------------
  }/* while( scaned < remMax ){ */
  //-----------------------------------------------------------------------
  /* save the final value of sidxGlobal */
#if 1
  if( gidx == 0 )
    atomicMax(snum, sidxGlobal);
#else
  if( tidx == 0 )
    printf("result = %d, myown = %d\n", atomicMax(snum, sidxGlobal), sidxGlobal);
#endif
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* MSD radix sort within a device */
//-------------------------------------------------------------------------
#ifdef  RADIX_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------
/* __device__ __forceinline__ void RADIX_SORT_GRID_MOD_COLLECT */
extern "C"
__global__ void radix_sort_grid_mod_collect
(int * RESTRICT bucket_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul, int2 * RESTRICT disp
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* MSD radix sort within device using inter-block synchronization */
  //-----------------------------------------------------------------------
  const int bidx =  BLOCKIDX_X1D;
  const int tidx = THREADIDX_X1D;
  const int lane = tidx & (warpSize - 1);
  //-----------------------------------------------------------------------
  __shared__ uint_int smem[NTHREADS_SORT];
  //-----------------------------------------------------------------------
  const int   sidx = bidx >> RADIX_SORT_CHECK_BITS;
  const int target = bidx & (RADIX_SORT_NUM_BUCKET - 1);
  //-----------------------------------------------------------------------
  /* if( sidx < snum ){ */
    //---------------------------------------------------------------------
    if( tidx == 0 ){
      //-------------------------------------------------------------------
      const int2 tmp = disp[sidx];
      //-------------------------------------------------------------------
      smem[0].i = tmp.x;
      smem[1].i = tmp.y;
      //-------------------------------------------------------------------
    }
    //---------------------------------------------------------------------
    __syncthreads();
    const int dispHead = smem[0].i;
    const int dispNum  = smem[1].i;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* calculate global prefix sum within a device */
    /* if # of elements is 128M --> # of elements per block is 128M / (4 * Ttot) ~ 512k >> 1024. */
    /* Therefore, simple evaluation of prefix sum within a block is not sufficient */
    int subTail = 0;
    const int niter = BLOCKSIZE(dispNum, NTHREADS_SORT);
    //---------------------------------------------------------------------
    for(int iter = 0; iter < niter; iter++){
      //-------------------------------------------------------------------
      const int idx = tidx + iter * NTHREADS_SORT;
      //-------------------------------------------------------------------
      __syncthreads();
      /* load local prefix sum */
      smem[tidx].i = (idx < dispNum) ? scanFul[dispHead + INDEX2D(RADIX_SORT_NUM_BUCKET, dispNum, target, idx)] : 0;
      /* prefix sum within a block */
      RADIX_PREFIX_SUM_GRID(tidx, lane, smem
#ifdef  USE_MASK_INSTEAD_OF_IF
			    , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
			    );
      /* write back global memory */
      if( idx < dispNum )
	scanFul[dispHead + INDEX2D(RADIX_SORT_NUM_BUCKET, dispNum, target, idx)] = subTail + smem[tidx].i;
      subTail += smem[NTHREADS_SORT - 1].i;
      //-------------------------------------------------------------------
    }/* for(int iter = 0; iter < niter; iter++){ */
    //---------------------------------------------------------------------
    if( tidx == 0 ){
      tail_gm  [sidx * RADIX_SORT_NUM_BUCKET + target] = subTail;
      /* 1 means that further global sort is required */
      bucket_gm[sidx * RADIX_SORT_NUM_BUCKET + target] = (subTail < (NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD)) ? 0 : 1;
    }/* if( tidx == 0 ){ */
    //---------------------------------------------------------------------
  /* }/\* if( sidx < snum ){ *\/ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//RADIX_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* MSD radix sort within a device */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_GRID_MOD_GATHER
(int * RESTRICT numIterFulSrc, int2 * RESTRICT infoFulSrc, int            remFul   ,
 int * RESTRICT numIterFulDst, int2 * RESTRICT infoFulDst, int * RESTRICT remFulAdd,
 int * RESTRICT numIterLoc   , int2 * RESTRICT infoLoc   , int * RESTRICT remLocAdd,
/* #ifdef  ATOMIC_BASED_JOB_ASSIGNMENT */
/*  int * RESTRICT checkedCounts, */
/* #endif//ATOMIC_BASED_JOB_ASSIGNMENT */
 int * RESTRICT snum, const int scanned,
 const int gidx, const int bidx, const int bnum, const int tidx, const int lane,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
 const int hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
 RADIX_DATA_TYPE_KEY * RESTRICT key_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_tmp_gm,
#ifndef SORT_ONLY
 int * RESTRICT idx_gm, int * RESTRICT idx_tmp_gm,
#endif//SORT_ONLY
#   if  RADIX_SORT_CHECK_BITS > 2
 uint4_array * RESTRICT sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
 volatile uint_int * RESTRICT smem, int * RESTRICT scan_gm, int * RESTRICT bucket_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  if( gidx == 0 )
    *snum = 0;
  //-----------------------------------------------------------------------
  /* MSD radix sort within a device (2nd scan phase) */
  int scaned = 0;
  int dispGlobal = 0;
  int sidxGlobal = 0;
  //-----------------------------------------------------------------------
  const int remMax = (remFul > RADIX_SORT_BOX_NUM_ALLOCATED) ? (RADIX_SORT_BOX_NUM_ALLOCATED) : (remFul);
  //-----------------------------------------------------------------------
  while( scaned < remMax ){
    //---------------------------------------------------------------------
    bool enabled = false;
    int bidx_sub = bidx;
    int bnum_sub = bnum;
    int numIter  = 0;
    int head = 0;
    int rem = 0;
    int disp_sub = dispGlobal;
    int sidx     = sidxGlobal;
    RADIX_SORT_GRID_GET_PARTITION(lane, tidx, bidx, bnum, remFul, &scaned, numIterFulSrc, infoFulSrc, &enabled, &bidx_sub, &bnum_sub, &numIter, &head, &rem, &disp_sub, &dispGlobal, &sidx, &sidxGlobal, smem
#ifdef  USE_MASK_INSTEAD_OF_IF
				  , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
				  );
    //---------------------------------------------------------------------
    if( enabled ){
      //-------------------------------------------------------------------
      const int2 info = {rem, head};
      const int nloop = BLOCKSIZE(rem, bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
      //-------------------------------------------------------------------
      if( tidx < RADIX_SORT_NUM_BUCKET ){
	//-----------------------------------------------------------------
	/* get partition among sub-groups */
	int subNum = tail_gm[sidx * RADIX_SORT_NUM_BUCKET + tidx];
	/* 1 means that further global sort is required */
	const int bucket_tmp = bucket_gm[sidx * RADIX_SORT_NUM_BUCKET + tidx];
	const int pidx = STOCK_ARRAY_PTR_HEAD + tidx;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	/* load index */
	int v0 = subNum;
	int v1 = bucket_tmp;
	int t0, t1;
	/* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
	t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET) & m0;      t1 = __shfl_up(v1, 1, RADIX_SORT_NUM_BUCKET) & m0;      v0 += t0;	v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >=  4
	t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET) & m1;      t1 = __shfl_up(v1, 2, RADIX_SORT_NUM_BUCKET) & m1;      v0 += t0;	v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >=  8
	t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET) & m2;      t1 = __shfl_up(v1, 4, RADIX_SORT_NUM_BUCKET) & m2;      v0 += t0;	v1 += t1;
#   if  RADIX_SORT_NUM_BUCKET >= 16
	t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET) & m3;      t1 = __shfl_up(v1, 8, RADIX_SORT_NUM_BUCKET) & m3;      v0 += t0;	v1 += t1;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
	t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 1, RADIX_SORT_NUM_BUCKET);      if( tidx >= 1 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >=  4
	t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 2, RADIX_SORT_NUM_BUCKET);      if( tidx >= 2 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >=  8
	t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 4, RADIX_SORT_NUM_BUCKET);      if( tidx >= 4 ){	v0 += t0;	v1 += t1;      }
#   if  RADIX_SORT_NUM_BUCKET >= 16
	t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET);      t1 = __shfl_up(v1, 8, RADIX_SORT_NUM_BUCKET);      if( tidx >= 8 ){	v0 += t0;	v1 += t1;      }
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
        /* return calculated inclusive prefix sum */
	smem[pidx].i = v0;
	smem[tidx].i = v1;
#else///USE_WARP_SHUFFLE_FUNC_SORT
	smem[pidx].i = subNum;
	smem[tidx].i = bucket_tmp;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
	smem[tidx].i += smem[tidx - 1].i & m0;	smem[pidx].i += smem[pidx - 1].i & m0;
#   if  RADIX_SORT_NUM_BUCKET >=  4
	smem[tidx].i += smem[tidx - 2].i & m1;	smem[pidx].i += smem[pidx - 2].i & m1;
#   if  RADIX_SORT_NUM_BUCKET >=  8
	smem[tidx].i += smem[tidx - 4].i & m2;	smem[pidx].i += smem[pidx - 4].i & m2;
#   if  RADIX_SORT_NUM_BUCKET >= 16
	smem[tidx].i += smem[tidx - 8].i & m3;	smem[pidx].i += smem[pidx - 8].i & m3;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
	if( tidx >= 1 ){	smem[tidx].i += smem[tidx - 1].i;	smem[pidx].i += smem[pidx - 1].i;      }
#   if  RADIX_SORT_NUM_BUCKET >=  4
	if( tidx >= 2 ){	smem[tidx].i += smem[tidx - 2].i;	smem[pidx].i += smem[pidx - 2].i;      }
#   if  RADIX_SORT_NUM_BUCKET >=  8
	if( tidx >= 4 ){	smem[tidx].i += smem[tidx - 4].i;	smem[pidx].i += smem[pidx - 4].i;      }
#   if  RADIX_SORT_NUM_BUCKET >= 16
	if( tidx >= 8 ){	smem[tidx].i += smem[tidx - 8].i;	smem[pidx].i += smem[pidx - 8].i;      }
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#endif//USE_WARP_SHUFFLE_FUNC_SORT
        //-----------------------------------------------------------------
	smem[STOCK_ARRAY_DST_HEAD + tidx].i = bucket_tmp;
	smem[STOCK_ARRAY_IDX_HEAD + tidx].i = bucket_tmp ? (smem[tidx].i - bucket_tmp) : (tidx - smem[tidx].i);
	//-----------------------------------------------------------------
	smem[STOCK_ARRAY_PTR_NUM  + tidx].i  = subNum;
	smem[STOCK_ARRAY_PTR_HEAD + tidx].i -= subNum;
	//-----------------------------------------------------------------
	if( bidx_sub == 0 ){
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	  int remFulHead, remLocHead;
#else///USE_WARP_SHUFFLE_FUNC_SORT
	  const int stock = smem[tidx].i;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
	  if( tidx == 0 ){
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	    remFulHead = atomicAdd(remFulAdd,                         smem[RADIX_SORT_NUM_BUCKET - 1].i);
	    remLocHead = atomicAdd(remLocAdd, RADIX_SORT_NUM_BUCKET - smem[RADIX_SORT_NUM_BUCKET - 1].i);
#else///USE_WARP_SHUFFLE_FUNC_SORT
	    smem[0].i  = atomicAdd(remFulAdd,                         smem[RADIX_SORT_NUM_BUCKET - 1].i);
	    smem[1].i  = atomicAdd(remLocAdd, RADIX_SORT_NUM_BUCKET - smem[RADIX_SORT_NUM_BUCKET - 1].i);
#endif//USE_WARP_SHUFFLE_FUNC_SORT
	  }
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	  remFulHead = __shfl(remFulHead, 0, RADIX_SORT_NUM_BUCKET);
	  remLocHead = __shfl(remLocHead, 0, RADIX_SORT_NUM_BUCKET);
#else///USE_WARP_SHUFFLE_FUNC_SORT
	  int remFulHead = smem[0].i;
	  int remLocHead = smem[1].i;
	  smem[tidx].i = stock;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
	  int2 subInfo = {subNum, head + smem[STOCK_ARRAY_PTR_HEAD + tidx].i};
	  if( smem[STOCK_ARRAY_DST_HEAD + tidx].i ){
	    infoFulDst   [remFulHead + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = subInfo;
	    numIterFulDst[remFulHead + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = numIter - 1;
	  }
	  else{
	    infoLoc   [remLocHead + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = subInfo;
	    numIterLoc[remLocHead + smem[STOCK_ARRAY_IDX_HEAD + tidx].i] = numIter - 1;
	  }
	}/* if( bidx == 0 ){ */
	//-----------------------------------------------------------------
      }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      rem = info.x;
      //-------------------------------------------------------------------
#ifdef  ENABLE_EARLY_EXIT_GRID
      int stock;
      if( tidx < RADIX_SORT_NUM_BUCKET ){
	int flag = smem[tidx].i != rem;/* 0 means sort is not required */
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	if( tidx == 0 )	  stock = smem[0].i;
#   if  RADIX_SORT_NUM_BUCKET >=  2
	int temp;
	temp = __shfl_up(flag, 1, RADIX_SORT_NUM_BUCKET);	flag *= temp;
#   if  RADIX_SORT_NUM_BUCKET >=  4
	temp = __shfl_up(flag, 2, RADIX_SORT_NUM_BUCKET);	flag *= temp;
#   if  RADIX_SORT_NUM_BUCKET >=  8
	temp = __shfl_up(flag, 4, RADIX_SORT_NUM_BUCKET);	flag *= temp;
#   if  RADIX_SORT_NUM_BUCKET == 16
	temp = __shfl_up(flag, 8, RADIX_SORT_NUM_BUCKET);	flag *= temp;
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET == 16
	if( tidx == 0 )	  smem[0].i = flag;
#else///USE_WARP_SHUFFLE_FUNC_SORT
	stock = smem[tidx].i;	smem[tidx].i = flag;
#   if  RADIX_SORT_NUM_BUCKET >=  2
	int temp;
	temp = smem[tidx ^ 1];	flag *= temp;	smem[tidx].i = flag;
#   if  RADIX_SORT_NUM_BUCKET >=  4
	temp = smem[tidx ^ 2];	flag *= temp;	smem[tidx].i = flag;
#   if  RADIX_SORT_NUM_BUCKET >=  8
	temp = smem[tidx ^ 4];	flag *= temp;	smem[tidx].i = flag;
#   if  RADIX_SORT_NUM_BUCKET == 16
	temp = smem[tidx ^ 8];	flag *= temp;	smem[tidx].i = flag;
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET == 16
#endif//USE_WARP_SHUFFLE_FUNC_SORT
      }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
      __syncthreads();
      const int exec = smem[0].i;
      if( exec )
#endif//ENABLE_EARLY_EXIT_GRID
	{
	  //---------------------------------------------------------------
#ifdef  ENABLE_EARLY_EXIT_GRID
	  __syncthreads();
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
	  if( tidx == 0 )
	    smem[0].i = stock;
#else///USE_WARP_SHUFFLE_FUNC_SORT
	  if( tidx < RADIX_SORT_NUM_BUCKET )
	    smem[tidx].i = stock;
#endif//USE_WARP_SHUFFLE_FUNC_SORT
	  __syncthreads();
#endif//ENABLE_EARLY_EXIT_GRID
	  //---------------------------------------------------------------
	  for(int loop = 0; loop < nloop; loop++){
	    //-------------------------------------------------------------
	    const int hidx_gm = head + RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * (bidx_sub + loop * bnum_sub);
	    //-------------------------------------------------------------
	    const int totNum = ((bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < rem) ? (bnum_sub * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : rem;
	    int subNum = totNum - bidx_sub * (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
	    if( subNum < 0 )	      subNum = 0;
	    subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    __syncthreads();
	    //-------------------------------------------------------------
	    if( tidx < RADIX_SORT_NUM_BUCKET ){
	      //-----------------------------------------------------------
	      /* load temporary stock on the global memory (num_sm, key_dst_sm and idx_dst_sm) */
	      smem[STOCK_ARRAY_SUM_HEAD + tidx].i = scan_gm[disp_sub + INDEX2D(nloop * bnum_sub, RADIX_SORT_NUM_BUCKET, bidx_sub + bnum_sub * loop, tidx)];
	      //-----------------------------------------------------------
	      /* load global prefix sum */
	      smem[tidx].i = scanFul[disp_sub + INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum_sub, tidx, bidx_sub + loop * bnum_sub)];
	      //-----------------------------------------------------------
	      smem[STOCK_ARRAY_NUM_HEAD + tidx].i = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
	      const int tmpVal = (tidx != 0) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 1].i;
	      smem[STOCK_ARRAY_SUM_HEAD + tidx].i += tmpVal;
	      //-----------------------------------------------------------
	      smem[tidx].i += (smem[STOCK_ARRAY_PTR_HEAD + tidx].i - smem[STOCK_ARRAY_SUM_HEAD + tidx].i);
	      //-----------------------------------------------------------
	    }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
	    //-------------------------------------------------------------
	    __syncthreads();
	    //-------------------------------------------------------------

	    //-------------------------------------------------------------
	    /* store sorted key, idx on global memory */
#pragma unroll
	    for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	      //-----------------------------------------------------------
	      const RADIX_DATA_TYPE_KEY key = key_tmp_gm[hidx_gm + ii];
	      const int target = (key >> ((numIter - 1) * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK;
	      //-----------------------------------------------------------
	      const int dstIdx = head + smem[target].i + ii;
	      key_gm[dstIdx] = key;
#ifndef SORT_ONLY
	      idx_gm[dstIdx] = idx_tmp_gm[hidx_gm + ii];
#endif//SORT_ONLY
	      //-----------------------------------------------------------
	    }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
	    //-------------------------------------------------------------
	    __syncthreads();
	    rem -= totNum;
	    //-------------------------------------------------------------
	  }/* for(int loop = 0; loop < nloop; loop++){ */
	}
      //-------------------------------------------------------------------
    }/* if( enabled ){ */
    //---------------------------------------------------------------------
  }/* while( scaned < remMax ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  RADIX_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------
extern "C"
__global__ void radixSortGrid_kickoff
(const int num, const int bits,
 int * RESTRICT numIterFul0, int2 * RESTRICT infoFul0,
 int * RESTRICT numIterFul1, int2 * RESTRICT infoFul1,
 int * RESTRICT remFulAdd, int * RESTRICT remLocAdd, int * RESTRICT snum)
{
  //-----------------------------------------------------------------------
  if( GLOBALIDX_X1D == 0 ){
    //---------------------------------------------------------------------
    int2 info0 = {num, 0};    infoFul0[0] = info0;    numIterFul0[0] = bits / RADIX_SORT_CHECK_BITS;
    int2 info1 = {  0, 0};    infoFul1[0] = info1;    numIterFul1[0] = 0;
    //---------------------------------------------------------------------
    *remFulAdd = 0;
    *remLocAdd = 0;
    //---------------------------------------------------------------------
    *snum = 0;
    //---------------------------------------------------------------------
  }/* if( GLOBALIDX_X1D == 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
__global__ void radixSortGrid_switch(int * RESTRICT remFulAdd)
{
  //-----------------------------------------------------------------------
  if( GLOBALIDX_X1D == 0 )
    *remFulAdd = 0;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
extern "C"
__global__ void radixSortGrid_finalize
(const int remLocMax, int * RESTRICT remFulAdd, int * RESTRICT remLocAdd, int * RESTRICT fail_dev)
{
  //-----------------------------------------------------------------------
  if( GLOBALIDX_X1D == 0 ){
    //---------------------------------------------------------------------
    const int remFulNew = *remFulAdd;    if( remFulNew > remLocMax )      atomicAdd(fail_dev, remFulNew);
    const int remLocNew = *remLocAdd;    if( remLocNew > remLocMax )      atomicAdd(fail_dev, remLocNew);
    //---------------------------------------------------------------------
  }/* if( GLOBALIDX_X1D == 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//RADIX_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* LSD radix sort within device */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_GRID_LSD_FULL
(const int numIter, const int numData,
 const int tidx, const int lane,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
 const int hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
 RADIX_DATA_TYPE_KEY * RESTRICT key_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_tmp_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_sm,
#ifndef SORT_ONLY
 int * RESTRICT idx_gm, int * RESTRICT idx_tmp_gm, int * RESTRICT idx_sm,
#endif//SORT_ONLY
#   if  RADIX_SORT_CHECK_BITS > 2
 uint4_array * RESTRICT sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
 volatile uint_int * RESTRICT smem, int * RESTRICT scan_gm, int * RESTRICT tail_gm, int * RESTRICT scanFul,
 const int gidx, const int bidx, const int bnum, volatile int * RESTRICT gsync0, volatile int * RESTRICT gsync1
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* LSD radix sort within device using inter-block synchronization */
  //-----------------------------------------------------------------------
  const int nloop = BLOCKSIZE(numData, bnum * RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT);
  //-----------------------------------------------------------------------
  for(int iter = 0; iter < numIter; iter++){
    //---------------------------------------------------------------------
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------
    for(int loop = 0; loop < nloop; loop++){
      //-------------------------------------------------------------------
#ifdef  FLIP_LOOP_ORDER
      const int grpIdx = loop + nloop * bidx;
#else///FLIP_LOOP_ORDER
      const int grpIdx = bidx +  loop * bnum;
#endif//FLIP_LOOP_ORDER
      const int hidx_gm = RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * grpIdx;
      //-------------------------------------------------------------------
      int subNum = numData - hidx_gm;
      if( subNum < 0 )	subNum = 0;
      subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
      //-------------------------------------------------------------------
      __syncthreads();
      //-------------------------------------------------------------------
      /* a remedy to treat the case when the number of elements is not a multiple of NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD */
#pragma unroll
      for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
	//-----------------------------------------------------------------
#   if  KEY_BITS == 32
	key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffff;
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
	key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffffffffffff;
#endif//KEY_BITS == 64
	//-----------------------------------------------------------------
      }/* for(int ii = subNum + tidx; ii < RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      /* load key_src, idx_src from global memory */
      /* set key_src and idx_src on shared memory for RADIX_SORT_CHIP */
#pragma unroll
      for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	//-----------------------------------------------------------------
	const int idx = hidx_gm + ii;
	//-----------------------------------------------------------------
	key_sm[ii] = key_gm[idx];
#ifndef SORT_ONLY
	idx_sm[ii] = idx_gm[idx];
#endif//SORT_ONLY
	//-----------------------------------------------------------------
      }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      __syncthreads();
      RADIX_SORT_CHIP(iter, tidx,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
		      hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
		      key_sm
#ifndef SORT_ONLY
		      , idx_sm
#endif//SORT_ONLY
		      , smem
#   if  RADIX_SORT_CHECK_BITS > 2
		      , sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
		      , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
		      );
      //-------------------------------------------------------------------
      /* remove the number of ghost values */
      /* target == RADIX_SORT_CHECK_MASK (= RADIX_SORT_NUM_BUCKET - 1) */
      /* # of ghosts is RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum */
      if( tidx == RADIX_SORT_CHECK_MASK )
	smem[STOCK_ARRAY_SUM_HEAD + tidx].i -= (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum);
      //-------------------------------------------------------------------
      /* put returned values of num_sm on global memory to calculate global prefix sum */
      if( tidx < RADIX_SORT_NUM_BUCKET )
      	scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, tidx, grpIdx)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
      //-------------------------------------------------------------------
      /* if nloop > 1 && loop < nloop - 1; save num_sm, key_dst and idx_dst on the global memory */
      if( loop < (nloop - 1) ){
	//-----------------------------------------------------------------
	/* save num_sm on the global memory */
	if( tidx < RADIX_SORT_NUM_BUCKET )
	  scan_gm[INDEX2D(nloop * bnum, RADIX_SORT_NUM_BUCKET, grpIdx, tidx)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i;
	//-----------------------------------------------------------------
	/* save key_dst and idx_dst on the global memory */
#pragma unroll
	for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	  key_tmp_gm[hidx_gm + ii] = key_sm[ii];
#ifndef SORT_ONLY
	  idx_tmp_gm[hidx_gm + ii] = idx_sm[ii];
#endif//SORT_ONLY
	}/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
	//-----------------------------------------------------------------
      }/* if( loop < (nloop - 1) ){ */
      //-------------------------------------------------------------------
    }/* for(int loop = 0; loop < nloop; loop++){ */
    //---------------------------------------------------------------------
    /* calculate global prefix sum within a device */
    /* if # of elements is 128M --> # of elements per block is 128M / (4 * Ttot) ~ 512k >> 1024. */
    /* Therefore, simple evaluation of prefix sum within a block is not sufficient */
    int tmp_sum_head;
    if( tidx < RADIX_SORT_NUM_BUCKET )
      tmp_sum_head = smem[STOCK_ARRAY_SUM_HEAD + tidx].i;
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    for(int target = bidx; target < RADIX_SORT_NUM_BUCKET; target += bnum){
      //-------------------------------------------------------------------
      int subTail = 0;
      const int niter = BLOCKSIZE(bnum * nloop, NTHREADS_SORT);
      //-------------------------------------------------------------------
      for(int ii = 0; ii < niter; ii++){
	//-----------------------------------------------------------------
	const int sidx = tidx + ii * NTHREADS_SORT;
	//-----------------------------------------------------------------
	/* load local prefix sum */
	smem[tidx].i = (sidx < bnum * nloop) ? scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, target, sidx)] : 0;
	/* prefix sum within a block */
	RADIX_PREFIX_SUM_GRID(tidx, lane, smem
#ifdef  USE_MASK_INSTEAD_OF_IF
			      , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
			      );
	/* write back global memory */
	if( sidx < bnum * nloop )
	  scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, target, sidx)] = subTail + smem[tidx].i;
	subTail += smem[NTHREADS_SORT - 1].i;
	__syncthreads();
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < niter; ii++){ */
      //-------------------------------------------------------------------
      if( tidx == 0 )
	tail_gm[target] = subTail;
      //-------------------------------------------------------------------
    }/* for(int target = bidx; target < RADIX_SORT_NUM_BUCKET; target += bnum){ */
    if( tidx < RADIX_SORT_NUM_BUCKET )
      smem[STOCK_ARRAY_SUM_HEAD + tidx].i = tmp_sum_head;
    globalSync(tidx, bidx, bnum, gsync0, gsync1);
    //---------------------------------------------------------------------
    if( tidx < RADIX_SORT_NUM_BUCKET ){
      //-------------------------------------------------------------------
      /* get partition among sub-groups */
      const int snum = tail_gm[tidx];
      const int sidx = STOCK_ARRAY_PTR_HEAD + tidx;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
      /* load index */
      int v0 = snum;
      int t0;
      /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET) & m0;      v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >=  4
      t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET) & m1;      v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >=  8
      t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET) & m2;      v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >= 16
      t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET) & m3;      v0 += t0;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET);      if( tidx >= 1 )	v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >=  4
      t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET);      if( tidx >= 2 )	v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >=  8
      t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET);      if( tidx >= 4 )	v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >= 16
      t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET);      if( tidx >= 8 )	v0 += t0;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
      /* return calculated inclusive prefix sum */
      smem[sidx].i = v0;
#else///USE_WARP_SHUFFLE_FUNC_SORT
      smem[sidx].i = snum;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      smem[sidx].i += smem[sidx - 1].i & m0;
#   if  RADIX_SORT_NUM_BUCKET >=  4
      smem[sidx].i += smem[sidx - 2].i & m1;
#   if  RADIX_SORT_NUM_BUCKET >=  8
      smem[sidx].i += smem[sidx - 4].i & m2;
#   if  RADIX_SORT_NUM_BUCKET >= 16
      smem[sidx].i += smem[sidx - 8].i & m3;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
      if( tidx >= 1 )	smem[sidx].i += smem[sidx - 1].i;
#   if  RADIX_SORT_NUM_BUCKET >=  4
      if( tidx >= 2 )	smem[sidx].i += smem[sidx - 2].i;
#   if  RADIX_SORT_NUM_BUCKET >=  8
      if( tidx >= 4 )	smem[sidx].i += smem[sidx - 4].i;
#   if  RADIX_SORT_NUM_BUCKET >= 16
      if( tidx >= 8 )	smem[sidx].i += smem[sidx - 8].i;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#endif//USE_WARP_SHUFFLE_FUNC_SORT
      //-------------------------------------------------------------------
      smem[STOCK_ARRAY_PTR_NUM  + tidx].i  = snum;
      smem[STOCK_ARRAY_PTR_HEAD + tidx].i -= snum;
      //-------------------------------------------------------------------
    }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
    //---------------------------------------------------------------------
    __syncthreads();
    //---------------------------------------------------------------------
#ifdef  FLIP_LOOP_ORDER
    int grpIdx = (nloop - 1) +  nloop      * bidx;
#else///FLIP_LOOP_ORDER
    int grpIdx = bidx        + (nloop - 1) * bnum;
#endif//FLIP_LOOP_ORDER
    int hidx_gm = RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * grpIdx;
    //---------------------------------------------------------------------
    int subNum = numData - hidx_gm;
    if( subNum < 0 )      subNum = 0;
    subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
    for(int loop = nloop - 1; loop >= 0; loop--){
      //-------------------------------------------------------------------
      if( tidx < RADIX_SORT_NUM_BUCKET ){
	//-----------------------------------------------------------------
	/* load global prefix sum */
	smem[tidx].i = scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, nloop * bnum, tidx, grpIdx)];
	//-----------------------------------------------------------------
	/* calculate local prefix sum of each element */
	smem[STOCK_ARRAY_NUM_HEAD + tidx].i = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
	const int tmpVal = (tidx != 0) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 1].i;
	smem[STOCK_ARRAY_SUM_HEAD + tidx].i += tmpVal;
	//-----------------------------------------------------------------
	smem[tidx].i += (smem[STOCK_ARRAY_PTR_HEAD + tidx].i - smem[STOCK_ARRAY_SUM_HEAD + tidx].i);
	//-----------------------------------------------------------------
      }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
      __syncthreads();
      //-------------------------------------------------------------------
      /* store sorted key, idx on global memory */
#pragma unroll
      for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	//-----------------------------------------------------------------
	const RADIX_DATA_TYPE_KEY key = key_sm[ii];
	const int target = (key >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK;
	//-----------------------------------------------------------------
	const int dstIdx = smem[target].i + ii;
	key_gm[dstIdx] = key;
#ifndef SORT_ONLY
	idx_gm[dstIdx] = idx_sm[ii];
#endif//SORT_ONLY
	//-----------------------------------------------------------------
      }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
      //-------------------------------------------------------------------
      /* load temporary stock on the global memory (num_sm, key_dst_sm and idx_dst_sm) */
      if( loop != 0 ){
	//-----------------------------------------------------------------
#ifdef  FLIP_LOOP_ORDER
	grpIdx = (loop - 1) + nloop      * bidx;
#else///FLIP_LOOP_ORDER
	grpIdx = bidx       + (loop - 1) * bnum;
#endif//FLIP_LOOP_ORDER
	hidx_gm = RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * grpIdx;
#ifdef  FLIP_LOOP_ORDER
	subNum = numData - hidx_gm;
	if( subNum < 0 )	  subNum = 0;
	subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
#else///FLIP_LOOP_ORDER
	subNum  = RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT;
#endif//FLIP_LOOP_ORDER
	//-----------------------------------------------------------------
#pragma unroll
	for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
	  key_sm[ii] = key_tmp_gm[hidx_gm + ii];
#ifndef SORT_ONLY
	  idx_sm[ii] = idx_tmp_gm[hidx_gm + ii];
#endif//SORT_ONLY
	}/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------
	if( tidx < RADIX_SORT_NUM_BUCKET )
	  smem[STOCK_ARRAY_SUM_HEAD + tidx].i = scan_gm[INDEX2D(nloop * bnum, RADIX_SORT_NUM_BUCKET, grpIdx, tidx)];
	//-----------------------------------------------------------------
	__syncthreads();
	//-----------------------------------------------------------------
      }/* if( loop != 0 ){ */
      //-------------------------------------------------------------------
    }/* for(int loop = nloop - 1; loop >= 0; loop--){ */
    //---------------------------------------------------------------------
  }/* for(int iter = 0; iter < numIter; iter++){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* LSD radix sort within device */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_GRID_LSD_1ST
(const int iter, const int numData,
 const int tidx, const int lane,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
 const int hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
 RADIX_DATA_TYPE_KEY * RESTRICT key_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_tmp_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_sm,
#ifndef SORT_ONLY
 int * RESTRICT idx_gm, int * RESTRICT idx_tmp_gm, int * RESTRICT idx_sm,
#endif//SORT_ONLY
#   if  RADIX_SORT_CHECK_BITS > 2
 uint4_array * RESTRICT sbuf,
#endif//RADIX_SORT_CHECK_BITS > 2
 volatile uint_int * RESTRICT smem, int * RESTRICT scan_gm, int * RESTRICT scanFul, const int bidx, const int bnum
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* LSD radix sort within device using inter-block synchronization */
  //-----------------------------------------------------------------------
  const int hidx_gm = RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * bidx;
  int subNum = numData - hidx_gm;
  if( subNum < 0 )    subNum = 0;
  subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
  //-----------------------------------------------------------------------
  /* a remedy to treat the case when the number of elements is not a multiple of NTHREADS_SORT * RADIX_SORT_ELEMENTS_PER_THREAD */
#pragma unroll
  for(int ii = 0; ii < RADIX_SORT_ELEMENTS_PER_THREAD; ii++){
    //---------------------------------------------------------------------
#   if  KEY_BITS == 32
    key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffff;
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
    key_sm[tidx +  ii * NTHREADS_SORT] = 0xffffffffffffffff;
#endif//KEY_BITS == 64
    //---------------------------------------------------------------------
  }/* for(int ii = subNum + tidx; ii < RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT; ii += NTHREADS_SORT){ */
  //-----------------------------------------------------------------------
  /* load key_src, idx_src from global memory */
  /* set key_src and idx_src on shared memory for RADIX_SORT_CHIP */
#pragma unroll
  for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
    //---------------------------------------------------------------------
    const int idx = hidx_gm + ii;
    //---------------------------------------------------------------------
    key_sm[ii] = key_gm[idx];
#ifndef SORT_ONLY
    idx_sm[ii] = idx_gm[idx];
#endif//SORT_ONLY
    //---------------------------------------------------------------------
  }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
  //-----------------------------------------------------------------------
  __syncthreads();
  RADIX_SORT_CHIP(iter, tidx,
#   if  RADIX_SORT_ELEMENTS_PER_THREAD != 1
		  hp,
#endif//RADIX_SORT_ELEMENTS_PER_THREAD != 1
		  key_sm
#ifndef SORT_ONLY
		  , idx_sm
#endif//SORT_ONLY
		  , smem
#   if  RADIX_SORT_CHECK_BITS > 2
		  , sbuf
#endif//RADIX_SORT_CHECK_BITS > 2
#ifdef  USE_MASK_INSTEAD_OF_IF
		  , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
		  );
  //-----------------------------------------------------------------------
  /* save key_dst and idx_dst on the global memory */
#pragma unroll
  for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
    key_tmp_gm[hidx_gm + ii] = key_sm[ii];
#ifndef SORT_ONLY
    idx_tmp_gm[hidx_gm + ii] = idx_sm[ii];
#endif//SORT_ONLY
  }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
  //-----------------------------------------------------------------------
  /* remove the number of ghost values */
  /* target == RADIX_SORT_CHECK_MASK (= RADIX_SORT_NUM_BUCKET - 1) */
  /* # of ghosts is RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum */
  if( tidx == RADIX_SORT_CHECK_MASK )
    smem[STOCK_ARRAY_SUM_HEAD + tidx].i -= (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT - subNum);
  //-----------------------------------------------------------------------
  /* put returned values of num_sm on global memory to calculate global prefix sum */
  /* save num_sm on the global memory */
  if( tidx < RADIX_SORT_NUM_BUCKET ){
    scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, bnum, tidx, bidx)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
    scan_gm[INDEX2D(bnum, RADIX_SORT_NUM_BUCKET, bidx, tidx)] = smem[STOCK_ARRAY_SUM_HEAD + tidx].i;
  }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* LSD radix sort within device */
//-------------------------------------------------------------------------
/* # of blocks is RADIX_SORT_NUM_BUCKET --> target == index of the block */
/* bnum means # of blocks for _1ST and _3RD */
__device__ __forceinline__ void RADIX_SORT_GRID_LSD_2ND
(const int tidx, const int lane, volatile uint_int * RESTRICT smem, int * RESTRICT tail_gm, int * RESTRICT scanFul, const int target, const int bnum
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* calculate global prefix sum within a device */
  /* if # of elements is 128M --> # of elements per block is 128M / (4 * Ttot) ~ 512k >> 1024. */
  /* Therefore, simple evaluation of prefix sum within a block is not sufficient */
  //-----------------------------------------------------------------------
  const int niter = BLOCKSIZE(bnum, NTHREADS_SORT_ACCUMULATION);
  //-----------------------------------------------------------------------
  int subTail = 0;
#pragma unroll
  for(int ii = 0; ii < niter; ii++){
    //---------------------------------------------------------------------
    const int sidx = tidx + ii * NTHREADS_SORT_ACCUMULATION;
    //---------------------------------------------------------------------
    /* load local prefix sum */
    smem[tidx].i = (sidx < bnum) ? scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, bnum, target, sidx)] : 0;
    /* prefix sum within a block */
    RADIX_PREFIX_SUM_GRID_ACCUMULATION(tidx, lane, smem
#ifdef  USE_MASK_INSTEAD_OF_IF
				       , m0, m1, m2, m3, m4
#endif//USE_MASK_INSTEAD_OF_IF
				       );
    /* write back global memory */
    if( sidx < bnum )
      scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, bnum, target, sidx)] = subTail + smem[tidx].i;
    subTail += smem[NTHREADS_SORT_ACCUMULATION - 1].i;
    __syncthreads();
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < niter; ii++){ */
  //-----------------------------------------------------------------------
  if( tidx == 0 )
    tail_gm[target] = subTail;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* LSD radix sort within device */
//-------------------------------------------------------------------------
__device__ __forceinline__ void RADIX_SORT_GRID_LSD_3RD
(const int iter, const int numData,
 const int tidx, const int lane,
 RADIX_DATA_TYPE_KEY * RESTRICT key_gm, RADIX_DATA_TYPE_KEY * RESTRICT key_tmp_gm,
#ifndef SORT_ONLY
 int * RESTRICT idx_gm, int * RESTRICT idx_tmp_gm,
#endif//SORT_ONLY
 volatile uint_int * RESTRICT smem, READ_ONLY int * RESTRICT scan_gm, READ_ONLY int * RESTRICT tail_gm, READ_ONLY int * RESTRICT scanFul,
 const int bidx, const int bnum
#ifdef  USE_MASK_INSTEAD_OF_IF
 , const int m0, const int m1, const int m2, const int m3, const int m4
#endif//USE_MASK_INSTEAD_OF_IF
 )
{
  //-----------------------------------------------------------------------
  /* LSD radix sort within device using inter-block synchronization */
  //-----------------------------------------------------------------------
  if( tidx < RADIX_SORT_NUM_BUCKET ){
    //---------------------------------------------------------------------
    /* get partition among sub-groups */
    const int snum = tail_gm[tidx];
    const int sidx = STOCK_ARRAY_PTR_HEAD + tidx;
#ifdef  USE_WARP_SHUFFLE_FUNC_SORT
    /* load index */
    int v0 = snum;
    int t0;
    /* calculate inclusive prefix sum */
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
    t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET) & m0;    v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >=  4
    t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET) & m1;    v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >=  8
    t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET) & m2;    v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >= 16
    t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET) & m3;    v0 += t0;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
    t0 = __shfl_up(v0, 1, RADIX_SORT_NUM_BUCKET);    if( tidx >= 1 )      v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >=  4
    t0 = __shfl_up(v0, 2, RADIX_SORT_NUM_BUCKET);    if( tidx >= 2 )      v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >=  8
    t0 = __shfl_up(v0, 4, RADIX_SORT_NUM_BUCKET);    if( tidx >= 4 )      v0 += t0;
#   if  RADIX_SORT_NUM_BUCKET >= 16
    t0 = __shfl_up(v0, 8, RADIX_SORT_NUM_BUCKET);    if( tidx >= 8 )      v0 += t0;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
    /* return calculated inclusive prefix sum */
    smem[sidx].i = v0;
#else///USE_WARP_SHUFFLE_FUNC_SORT
    smem[sidx].i = snum;
#ifdef  USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
    smem[sidx].i += smem[sidx - 1].i & m0;
#   if  RADIX_SORT_NUM_BUCKET >=  4
    smem[sidx].i += smem[sidx - 2].i & m1;
#   if  RADIX_SORT_NUM_BUCKET >=  8
    smem[sidx].i += smem[sidx - 4].i & m2;
#   if  RADIX_SORT_NUM_BUCKET >= 16
    smem[sidx].i += smem[sidx - 8].i & m3;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#else///USE_MASK_INSTEAD_OF_IF
#   if  RADIX_SORT_NUM_BUCKET >=  2
    if( tidx >= 1 )      smem[sidx].i += smem[sidx - 1].i;
#   if  RADIX_SORT_NUM_BUCKET >=  4
    if( tidx >= 2 )      smem[sidx].i += smem[sidx - 2].i;
#   if  RADIX_SORT_NUM_BUCKET >=  8
    if( tidx >= 4 )      smem[sidx].i += smem[sidx - 4].i;
#   if  RADIX_SORT_NUM_BUCKET >= 16
    if( tidx >= 8 )      smem[sidx].i += smem[sidx - 8].i;
#endif//RADIX_SORT_NUM_BUCKET >= 16
#endif//RADIX_SORT_NUM_BUCKET >=  8
#endif//RADIX_SORT_NUM_BUCKET >=  4
#endif//RADIX_SORT_NUM_BUCKET >=  2
#endif//USE_MASK_INSTEAD_OF_IF
#endif//USE_WARP_SHUFFLE_FUNC_SORT
    //---------------------------------------------------------------------
    smem[STOCK_ARRAY_PTR_NUM  + tidx].i  = snum;
    smem[STOCK_ARRAY_PTR_HEAD + tidx].i -= snum;
    //---------------------------------------------------------------------
    smem[                       tidx].i = scanFul[INDEX2D(RADIX_SORT_NUM_BUCKET, bnum, tidx, bidx)];
    smem[STOCK_ARRAY_SUM_HEAD + tidx].i = scan_gm[INDEX2D(bnum, RADIX_SORT_NUM_BUCKET, bidx, tidx)];
    /* calculate local prefix sum of each element */
    smem[STOCK_ARRAY_NUM_HEAD + tidx].i = smem[STOCK_ARRAY_SUM_HEAD + tidx].i - (tidx >= 2) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 2].i;
    const int tmpVal = (tidx != 0) * smem[STOCK_ARRAY_SUM_HEAD + tidx - 1].i;
    smem[STOCK_ARRAY_SUM_HEAD + tidx].i += tmpVal;
    //---------------------------------------------------------------------
    smem[tidx].i += (smem[STOCK_ARRAY_PTR_HEAD + tidx].i - smem[STOCK_ARRAY_SUM_HEAD + tidx].i);
    //---------------------------------------------------------------------
  }/* if( tidx < RADIX_SORT_NUM_BUCKET ){ */
  //-----------------------------------------------------------------------
  const int hidx_gm = RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT * bidx;
  int subNum = numData - hidx_gm;
  if( subNum < 0 )    subNum = 0;
  subNum = ((RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) < subNum) ? (RADIX_SORT_ELEMENTS_PER_THREAD * NTHREADS_SORT) : subNum;
  //-----------------------------------------------------------------------
  __syncthreads();
  //-----------------------------------------------------------------------
  /* store sorted key, idx on global memory */
  /* int target = 0; */
#pragma unroll
  for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){
    //---------------------------------------------------------------------
    const RADIX_DATA_TYPE_KEY key = key_tmp_gm[hidx_gm + ii];
    const int target = (key >> (iter * RADIX_SORT_CHECK_BITS)) & RADIX_SORT_CHECK_MASK;
    //---------------------------------------------------------------------
    const int dstIdx = smem[target].i + ii;
    key_gm[dstIdx] = key;
#ifndef SORT_ONLY
    idx_gm[dstIdx] = idx_tmp_gm[hidx_gm + ii];
#endif//SORT_ONLY
    //---------------------------------------------------------------------
  }/* for(int ii = tidx; ii < subNum; ii += NTHREADS_SORT){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  RADIX_INC_CU_FIRST_CALL
#undef  RADIX_INC_CU_FIRST_CALL
#define RADIX_INC_CU_MULTI_CALL
#endif//RADIX_INC_CU_FIRST_CALL
//-------------------------------------------------------------------------
