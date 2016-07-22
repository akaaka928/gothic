/*************************************************************************\
 *                                                                       *
                  last updated on 2015/08/19(Wed) 11:50:31
 *                                                                       *
 *    Header File for duplicating radix sort library on GPU              *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef RADIX_DEL_H
#       include "../sort/radix_del.h"
#endif//RADIX_DEL_H
//-------------------------------------------------------------------------
#ifndef RADIX_INC_H
#define RADIX_INC_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  KEY_BITS == 32
#define RADIX_DATA_TYPE_KEY uint
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
#define RADIX_DATA_TYPE_KEY ulong
#endif//KEY_BITS == 64
//-------------------------------------------------------------------------
#ifdef  SORT_ELEMENTS_PER_THREAD
#define RADIX_SORT_ELEMENTS_PER_THREAD SORT_ELEMENTS_PER_THREAD
#else///SORT_ELEMENTS_PER_THREAD
/* #define RADIX_SORT_ELEMENTS_PER_THREAD (1) */
/* #define RADIX_SORT_ELEMENTS_PER_THREAD (2) */
#define RADIX_SORT_ELEMENTS_PER_THREAD (4)
#endif//SORT_ELEMENTS_PER_THREAD
//-------------------------------------------------------------------------
#   if  (KEY_BITS == 64) && (NTHREADS_SORT == 1024) && (RADIX_SORT_ELEMENTS_PER_THREAD == 4) && !defined(SORT_ONLY)
#undef  RADIX_SORT_ELEMENTS_PER_THREAD
#define RADIX_SORT_ELEMENTS_PER_THREAD (2)
#endif//(KEY_BITS == 64) && (NTHREADS_SORT == 1024) && (RADIX_SORT_ELEMENTS_PER_THREAD == 4) && !defined(SORT_ONLY)
//-------------------------------------------------------------------------
#define STOCK_ARRAY_SUM_HEAD (    RADIX_SORT_NUM_BUCKET)
#define STOCK_ARRAY_NUM_HEAD (2 * RADIX_SORT_NUM_BUCKET)
#define STOCK_ARRAY_DST_HEAD (3 * RADIX_SORT_NUM_BUCKET)
#define STOCK_ARRAY_IDX_HEAD (4 * RADIX_SORT_NUM_BUCKET)
#define STOCK_ARRAY_PTR_NUM  (5 * RADIX_SORT_NUM_BUCKET)
#define STOCK_ARRAY_PTR_HEAD (6 * RADIX_SORT_NUM_BUCKET)
//-------------------------------------------------------------------------
#ifdef  SORT_ONLY
#   if  KEY_BITS == 32
/* #define RADIX_SORT_TSUB_GET_DSTIDX getDstIdxTsub32 */
/* #define RADIX_SORT_BLCK_GET_DSTIDX getDstIdxBlck32 */
#define RADIX_PREFIX_SUM_TSUB prefixSumTsub32
#define RADIX_PREFIX_SUM_WARP prefixSumWarp32
#define RADIX_PREFIX_SUM_BLCK prefixSumBlck32
/* #define RADIX_PREFIX_SUM_GRID prefixSumGrid32 */
#define RADIX_SORT_TSUB __radixSortTsub32
#define RADIX_SORT_BLCK __radixSortBlck32
#define RADIX_SORT_CHIP __radixSortChip32
#define RADIX_SORT_GRID __radixSortGrid32
#define RADIX_SORT_GRID_MSD __radixSortGrid32_msd
#define RADIX_SORT_GRID_LSD __radixSortGrid32_lsd
#define RADIX_SORT_GRID_MOD __radixSortGrid32_mod
#define RADIX_SORT_GRID_GET_PARTITION __radixSortGrid32_get_partition
#define RADIX_SORT_GRID_MOD_SCATTER   __radixSortGrid32_mod_scatter
#define RADIX_SORT_GRID_MOD_COLLECT   __radixSortGrid32_mod_collect
#define RADIX_SORT_GRID_MOD_GATHER    __radixSortGrid32_mod_gather
#define RADIX_SORT_GRID_LSD_FULL __radixSortLsdGrid32
#define RADIX_SORT_GRID_LSD_1ST  __radixSortLsdGrid32_1st
#define RADIX_SORT_GRID_LSD_2ND  __radixSortLsdGrid32_2nd
#define RADIX_SORT_GRID_LSD_3RD  __radixSortLsdGrid32_3rd
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
/* #define RADIX_SORT_TSUB_GET_DSTIDX getDstIdxTsub64 */
/* #define RADIX_SORT_BLCK_GET_DSTIDX getDstIdxBlck64 */
#define RADIX_PREFIX_SUM_TSUB prefixSumTsub64
#define RADIX_PREFIX_SUM_WARP prefixSumWarp64
#define RADIX_PREFIX_SUM_BLCK prefixSumBlck64
/* #define RADIX_PREFIX_SUM_GRID prefixSumGrid64 */
#define RADIX_SORT_TSUB __radixSortTsub64
#define RADIX_SORT_BLCK __radixSortBlck64
#define RADIX_SORT_CHIP __radixSortChip64
#define RADIX_SORT_GRID __radixSortGrid64
#define RADIX_SORT_GRID_MSD __radixSortGrid64_msd
#define RADIX_SORT_GRID_LSD __radixSortGrid64_lsd
#define RADIX_SORT_GRID_MOD __radixSortGrid64_mod
#define RADIX_SORT_GRID_GET_PARTITION __radixSortGrid64_get_partition
#define RADIX_SORT_GRID_MOD_SCATTER   __radixSortGrid64_mod_scatter
#define RADIX_SORT_GRID_MOD_COLLECT   __radixSortGrid64_mod_collect
#define RADIX_SORT_GRID_MOD_GATHER    __radixSortGrid64_mod_gather
#define RADIX_SORT_GRID_LSD_FULL __radixSortLsdGrid64
#define RADIX_SORT_GRID_LSD_1ST  __radixSortLsdGrid64_1st
#define RADIX_SORT_GRID_LSD_2ND  __radixSortLsdGrid64_2nd
#define RADIX_SORT_GRID_LSD_3RD  __radixSortLsdGrid64_3rd
#endif//KEY_BITS == 64
#else///SORT_ONLY
#   if  KEY_BITS == 32
/* #define RADIX_SORT_TSUB_GET_DSTIDX getDstIdxTsub32idx */
/* #define RADIX_SORT_BLCK_GET_DSTIDX getDstIdxBlck32idx */
#define RADIX_PREFIX_SUM_TSUB prefixSumTsub32idx
#define RADIX_PREFIX_SUM_WARP prefixSumWarp32idx
#define RADIX_PREFIX_SUM_BLCK prefixSumBlck32idx
/* #define RADIX_PREFIX_SUM_GRID prefixSumGrid32idx */
#define RADIX_SORT_TSUB __radixSortTsub32idx
#define RADIX_SORT_BLCK __radixSortBlck32idx
#define RADIX_SORT_CHIP __radixSortChip32idx
#define RADIX_SORT_GRID __radixSortGrid32idx
#define RADIX_SORT_GRID_MSD __radixSortGrid32idx_msd
#define RADIX_SORT_GRID_LSD __radixSortGrid32idx_lsd
#define RADIX_SORT_GRID_MOD __radixSortGrid32idx_mod
#define RADIX_SORT_GRID_GET_PARTITION __radixSortGrid32idx_get_partition
#define RADIX_SORT_GRID_MOD_SCATTER   __radixSortGrid32idx_mod_scatter
#define RADIX_SORT_GRID_MOD_COLLECT   __radixSortGrid32idx_mod_collect
#define RADIX_SORT_GRID_MOD_GATHER    __radixSortGrid32idx_mod_gather
#define RADIX_SORT_GRID_LSD_FULL __radixSortLsdGrid32idx
#define RADIX_SORT_GRID_LSD_1ST  __radixSortLsdGrid32idx_1st
#define RADIX_SORT_GRID_LSD_2ND  __radixSortLsdGrid32idx_2nd
#define RADIX_SORT_GRID_LSD_3RD  __radixSortLsdGrid32idx_3rd
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
/* #define RADIX_SORT_TSUB_GET_DSTIDX getDstIdxTsub64idx */
/* #define RADIX_SORT_BLCK_GET_DSTIDX getDstIdxBlck64idx */
#define RADIX_PREFIX_SUM_TSUB prefixSumTsub64idx
#define RADIX_PREFIX_SUM_WARP prefixSumWarp64idx
#define RADIX_PREFIX_SUM_BLCK prefixSumBlck64idx
/* #define RADIX_PREFIX_SUM_GRID prefixSumGrid64idx */
#define RADIX_SORT_TSUB __radixSortTsub64idx
#define RADIX_SORT_BLCK __radixSortBlck64idx
#define RADIX_SORT_CHIP __radixSortChip64idx
#define RADIX_SORT_GRID __radixSortGrid64idx
#define RADIX_SORT_GRID_MSD __radixSortGrid64idx_msd
#define RADIX_SORT_GRID_LSD __radixSortGrid64idx_lsd
#define RADIX_SORT_GRID_MOD __radixSortGrid64idx_mod
#define RADIX_SORT_GRID_GET_PARTITION __radixSortGrid64idx_get_partition
#define RADIX_SORT_GRID_MOD_SCATTER   __radixSortGrid64idx_mod_scatter
#define RADIX_SORT_GRID_MOD_COLLECT   __radixSortGrid64idx_mod_collect
#define RADIX_SORT_GRID_MOD_GATHER    __radixSortGrid64idx_mod_gather
#define RADIX_SORT_GRID_LSD_FULL __radixSortLsdGrid64idx
#define RADIX_SORT_GRID_LSD_1ST  __radixSortLsdGrid64idx_1st
#define RADIX_SORT_GRID_LSD_2ND  __radixSortLsdGrid64idx_2nd
#define RADIX_SORT_GRID_LSD_3RD  __radixSortLsdGrid64idx_3rd
#endif//KEY_BITS == 64
#endif//SORT_ONLY
//-------------------------------------------------------------------------
/* warp shuffle instruction only supports int or float, DOES NOT SUPPORT uint */
/* exploiting union reduces the limitation: use int for warp shuffle, use uint for addition */
//-------------------------------------------------------------------------
/* 8bits can count 0, 1, 2, ..., 2^8-1 = 255 */
/*   sufficient size for prefix sum within a  warp, 1 variables per thread (up to   32) */
/*   sufficient size for prefix sum within a  warp, 2 variables per thread (up to   64) */
/*   sufficient size for prefix sum within a  warp, 4 variables per thread (up to  128) */
/* insufficient size for prefix sum within a block, 1 variables per thread (up to 1024) */
#define RADIX_WARP_SHIFT_BITS (8)
#define RADIX_WARP_SHIFT_MASK (0xff)
#define RADIX_WARP_CHECK_BITS (2)
#define RADIX_WARP_CHECK_MASK (3)
/* a 32bits variable can contain 32 / RADIX_WARP_SHIFT_BITS keys */
#define RADIX_WARP_NUM_KEYS (4)
//-------------------------------------------------------------------------
/* loop length reauired to sweep KEY_BITS bits key is KEY_BITS / RADIX_WARP_CHECK_BITS */
#   if  KEY_BITS == 32
#define RADIX_WARP_NUM_ITER (16)
#endif//KEY_BITS == 32
#   if  KEY_BITS == 64
#define RADIX_WARP_NUM_ITER (32)
#endif//KEY_BITS == 64
//-------------------------------------------------------------------------
/* 15 bits can count 0, 1, 2, ..., 2^15-1 = 32k - 1 = 32767 */
/* 16 bits can count 0, 1, 2, ..., 2^16-1 = 64k - 1 = 65535 */
/*   sufficient size for prefix sum within  a block ,  1 variables per thread (up to  1024) */
/*   sufficient size for prefix sum within  a block ,  4 variables per thread (up to  4096) */
/*   sufficient size for prefix sum within  a block , 16 variables per thread (up to 16384) */
/* insufficient size for prefix sum within 16 blocks,  4 variables per thread (up to 65536) */
/* insufficient size for prefix sum within  8 blocks,  8 variables per thread (up to 65536) */
/* insufficient size for prefix sum within  4 blocks, 16 variables per thread (up to 65536) */
/* insufficient size for prefix sum within  2 blocks, 32 variables per thread (up to 65536) */
/* insufficient size for prefix sum within  a block , 64 variables per thread (up to 65536) */
#define RADIX_BLOCK_SHIFT_BITS (16)
#define RADIX_BLOCK_SHIFT_MASK (0xffff)
/* a 32bits variable can contain 32 / RADIX_BLOCK_SHIFT_BITS keys */
#define RADIX_CONVERT_8x4to16x2UPPER(a) ((((a) & 0xff000000) >> 8) | (((a) & 0xff0000)) >> 16)
#define RADIX_CONVERT_8x4to16x2LOWER(a) ((((a) & 0x0000ff00) << 8) |  ((a) & 0x0000ff)       )
#define RADIX_BLOCK_NUM_KEYS (2)
//-------------------------------------------------------------------------
/* 31 bits can count 0, 1, 2, ..., 2^31-1 = 4G - 1 */
/* sufficient size for up to 2G elements: greater size for the capacity of global memory; 2G elements = 8GiB for 32 bits variables */
#define RADIX_CONVERT_16x2to32x1UPPER(a) ((a) >> 16)
#define RADIX_CONVERT_16x2to32x1LOWER(a) ((a) & 0xffff)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  RADIX_SORT_INT
#undef  RADIX_SORT_INT
#endif//RADIX_SORT_INT
//-------------------------------------------------------------------------
#ifdef  RADIX_SORT_UINT
#undef  RADIX_SORT_UINT
#endif//RADIX_SORT_UINT
//-------------------------------------------------------------------------
#ifdef  RADIX_SORT_UINT_INT
#undef  RADIX_SORT_UINT_INT
#endif//RADIX_SORT_UINT_INT
//-------------------------------------------------------------------------
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 1
#define RADIX_SORT_INT         int
#define RADIX_SORT_UINT       uint
#define RADIX_SORT_UINT_INT   uint_int
#define RADIX_SORT_UINT_ARRAY uint
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 1
//-------------------------------------------------------------------------
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 2
#define RADIX_SORT_INT         int2
#define RADIX_SORT_UINT       uint2
#define RADIX_SORT_UINT_INT   uint2_int2
#define RADIX_SORT_UINT_ARRAY uint2_array
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 2
//-------------------------------------------------------------------------
#   if  RADIX_SORT_ELEMENTS_PER_THREAD == 4
#define RADIX_SORT_INT         int4
#define RADIX_SORT_UINT       uint4
#define RADIX_SORT_UINT_INT   uint4_int4
#define RADIX_SORT_UINT_ARRAY uint4_array
#endif//RADIX_SORT_ELEMENTS_PER_THREAD == 4
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  RADIX_DEL_H
#undef  RADIX_DEL_H
#endif//RADIX_DEL_H
//-------------------------------------------------------------------------
#endif//RADIX_INC_H
//-------------------------------------------------------------------------
