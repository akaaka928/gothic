/**
 * @file scan_inc.cuh
 *
 * @brief Header file for parallel prefix sum library on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/22 (Wed)
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


#   if  defined(USE_WARP_SHUFFLE_FUNC_SCAN_INC) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_SCAN_INC
#endif//defined(USE_WARP_SHUFFLE_FUNC_SCAN_INC) && (GPUGEN < 30)


#   if  NTHREADS_SCAN_INC ==   32
#define PREFIX_SUM_BLCK prefixSumBlck0032
#define PREFIX_SUM_GRID prefixSumGrid0032
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition0032
#endif//NTHREADS_SCAN_INC ==   32

#   if  NTHREADS_SCAN_INC ==   64
#define PREFIX_SUM_BLCK prefixSumBlck0064
#define PREFIX_SUM_GRID prefixSumGrid0064
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition0064
#endif//NTHREADS_SCAN_INC ==   64

#   if  NTHREADS_SCAN_INC ==  128
#define PREFIX_SUM_BLCK prefixSumBlck0128
#define PREFIX_SUM_GRID prefixSumGrid0128
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition0128
#endif//NTHREADS_SCAN_INC ==  128

#   if  NTHREADS_SCAN_INC ==  256
#define PREFIX_SUM_BLCK prefixSumBlck0256
#define PREFIX_SUM_GRID prefixSumGrid0256
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition0256
#endif//NTHREADS_SCAN_INC ==  256

#   if  NTHREADS_SCAN_INC ==  512
#define PREFIX_SUM_BLCK prefixSumBlck0512
#define PREFIX_SUM_GRID prefixSumGrid0512
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition0512
#endif//NTHREADS_SCAN_INC ==  512

#   if  NTHREADS_SCAN_INC == 1024
#define PREFIX_SUM_BLCK prefixSumBlck1024
#define PREFIX_SUM_GRID prefixSumGrid1024
#define PREFIX_SUM_GRID_WITH_PARTITION prefixSumGridPartition1024
#endif//NTHREADS_SCAN_INC == 1024


#ifdef  SCAN_DEL_CUH
#undef  SCAN_DEL_CUH
#endif//SCAN_DEL_CUH


#endif//SCAN_INC_CUH
