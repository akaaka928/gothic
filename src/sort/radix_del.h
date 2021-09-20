/*************************************************************************\
 *                                                                       *
                  last updated on 2015/08/06(Thu) 18:01:35
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
#define RADIX_DEL_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  RADIX_INC_H
//-------------------------------------------------------------------------
#undef  RADIX_DATA_TYPE_KEY
#undef  RADIX_SORT_ELEMENTS_PER_THREAD
#undef  STOCK_ARRAY_PTR_NUM
#undef  STOCK_ARRAY_PTR_HEAD
#undef  STOCK_ARRAY_SUM_HEAD
#undef  STOCK_ARRAY_NUM_HEAD
#undef  STOCK_ARRAY_DST_HEAD
#undef  STOCK_ARRAY_IDX_HEAD
//-------------------------------------------------------------------------
/* #undef  RADIX_SORT_TSUB_GET_DSTIDX */
/* #undef  RADIX_SORT_BLCK_GET_DSTIDX */
#undef  RADIX_PREFIX_SUM_TSUB
#undef  RADIX_PREFIX_SUM_WARP
#undef  RADIX_PREFIX_SUM_BLCK
/* #undef  RADIX_PREFIX_SUM_GRID */
#undef  RADIX_SORT_TSUB
#undef  RADIX_SORT_BLCK
#undef  RADIX_SORT_CHIP
#undef  RADIX_SORT_GRID
#undef  RADIX_SORT_GRID_MSD
#undef  RADIX_SORT_GRID_LSD
#undef  RADIX_SORT_GRID_MOD
#undef  RADIX_SORT_GRID_GET_PARTITION
#undef  RADIX_SORT_GRID_MOD_SCATTER
#undef  RADIX_SORT_GRID_MOD_COLLECT
#undef  RADIX_SORT_GRID_MOD_GATHER
#undef  RADIX_SORT_GRID_LSD_FULL
#undef  RADIX_SORT_GRID_LSD_1ST
#undef  RADIX_SORT_GRID_LSD_2ND
#undef  RADIX_SORT_GRID_LSD_3RD
//-------------------------------------------------------------------------
#undef  RADIX_WARP_SHIFT_BITS
#undef  RADIX_WARP_SHIFT_MASK
#undef  RADIX_WARP_CHECK_BITS
#undef  RADIX_WARP_CHECK_MASK
//-------------------------------------------------------------------------
#undef  RADIX_WARP_NUM_KEYS
#undef  RADIX_BLOCK_NUM_KEYS
//-------------------------------------------------------------------------
#undef  RADIX_WARP_NUM_ITER
/* #undef  RADIX_BLOCK_NUM_ITER */
//-------------------------------------------------------------------------
#undef  RADIX_BLOCK_SHIFT_BITS
#undef  RADIX_BLOCK_SHIFT_MASK
/* #undef  RADIX_SORT_CHECK_BITS */
/* #undef  RADIX_SORT_CHECK_MASK */
//-------------------------------------------------------------------------
/* #undef  RADIX_BLOCK_PICKUP_MASK0 */
/* #undef  RADIX_BLOCK_PICKUP_MASK1 */
/* #undef  RADIX_BLOCK_PICKUP_MASK2 */
/* #undef  RADIX_BLOCK_PICKUP_MASK3 */
/* #undef  RADIX_BLOCK_PICKUP_SHIFT1 */
/* #undef  RADIX_BLOCK_PICKUP_SHIFT2 */
/* #undef  RADIX_BLOCK_PICKUP_SHIFT3 */
#undef  RADIX_CONVERT_8x4to16x2UPPER
#undef  RADIX_CONVERT_8x4to16x2LOWER
//-------------------------------------------------------------------------
#undef  RADIX_CONVERT_16x2to32x1UPPER
#undef  RADIX_CONVERT_16x2to32x1LOWER
//-------------------------------------------------------------------------
#undef  RADIX_SORT_INT
#undef  RADIX_SORT_UINT
#undef  RADIX_SORT_UINT_INT
#undef  RADIX_SORT_UINT_ARRAY
//-------------------------------------------------------------------------
#undef  RADIX_INC_H
#endif//RADIX_INC_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//RADIX_DEL_H
//-------------------------------------------------------------------------
