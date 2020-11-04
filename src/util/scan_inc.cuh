/**
 * @file scan_inc.cuh
 *
 * @brief Header file for parallel prefix sum library on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/11/04 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */


#ifndef SCAN_DEL_CUH
#include "../util/scan_del.cuh"
#endif//SCAN_DEL_CUH


#ifndef SCAN_INC_CUH
#define SCAN_INC_CUH


#   if  defined(USE_WARP_SHUFFLE_FUNC_SCAN_INC) && (GPUVER < 30)
#undef          USE_WARP_SHUFFLE_FUNC_SCAN_INC
#endif//defined(USE_WARP_SHUFFLE_FUNC_SCAN_INC) && (GPUVER < 30)

#   if  defined(USE_WARP_REDUCE_FUNCTIONS_SCAN_INC) && !defined(ENABLE_WARP_REDUCE_FUNCTIONS)
#undef          USE_WARP_REDUCE_FUNCTIONS_SCAN_INC
#endif//defined(USE_WARP_REDUCE_FUNCTIONS_SCAN_INC) && !defined(ENABLE_WARP_REDUCE_FUNCTIONS)


#   if  NTHREADS_SCAN_INC ==   32
#define SHFL_MASK_SCAN_INC SHFL_MASK_01
#define PREFIX_SUM_BLCK prefixSumBlck0032
#define PREFIX_SUM_GRID prefixSumGrid0032
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition0032
#define  TOTAL_SUM_BLCK  totalSumBlck0032
#define  TOTAL_SUM_GRID  totalSumGrid0032
#endif//NTHREADS_SCAN_INC ==   32

#   if  NTHREADS_SCAN_INC ==   64
#define SHFL_MASK_SCAN_INC SHFL_MASK_02
#define PREFIX_SUM_BLCK prefixSumBlck0064
#define PREFIX_SUM_GRID prefixSumGrid0064
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition0064
#define  TOTAL_SUM_BLCK  totalSumBlck0064
#define  TOTAL_SUM_GRID  totalSumGrid0064
#endif//NTHREADS_SCAN_INC ==   64

#   if  NTHREADS_SCAN_INC ==  128
#define SHFL_MASK_SCAN_INC SHFL_MASK_04
#define PREFIX_SUM_BLCK prefixSumBlck0128
#define PREFIX_SUM_GRID prefixSumGrid0128
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition0128
#define  TOTAL_SUM_BLCK  totalSumBlck0128
#define  TOTAL_SUM_GRID  totalSumGrid0128
#endif//NTHREADS_SCAN_INC ==  128

#   if  NTHREADS_SCAN_INC ==  256
#define SHFL_MASK_SCAN_INC SHFL_MASK_08
#define PREFIX_SUM_BLCK prefixSumBlck0256
#define PREFIX_SUM_GRID prefixSumGrid0256
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition0256
#define  TOTAL_SUM_BLCK  totalSumBlck0256
#define  TOTAL_SUM_GRID  totalSumGrid0256
#endif//NTHREADS_SCAN_INC ==  256

#   if  NTHREADS_SCAN_INC ==  512
#define SHFL_MASK_SCAN_INC SHFL_MASK_16
#define PREFIX_SUM_BLCK prefixSumBlck0512
#define PREFIX_SUM_GRID prefixSumGrid0512
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition0512
#define  TOTAL_SUM_BLCK  totalSumBlck0512
#define  TOTAL_SUM_GRID  totalSumGrid0512
#endif//NTHREADS_SCAN_INC ==  512

#   if  NTHREADS_SCAN_INC == 1024
#define SHFL_MASK_SCAN_INC SHFL_MASK_32
#define PREFIX_SUM_BLCK prefixSumBlck1024
#define PREFIX_SUM_GRID prefixSumGrid1024
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition1024
#define  TOTAL_SUM_BLCK  totalSumBlck1024
#define  TOTAL_SUM_GRID  totalSumGrid1024
#endif//NTHREADS_SCAN_INC == 1024


#ifdef  SCAN_DEL_CUH
#undef  SCAN_DEL_CUH
#endif//SCAN_DEL_CUH


#endif//SCAN_INC_CUH
