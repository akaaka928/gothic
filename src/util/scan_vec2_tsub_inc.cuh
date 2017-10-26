/**
 * @file scan_vec2_tsub_inc.cuh
 *
 * @brief Header file for parallel prefix sum library for 2-components vector on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/10/26 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */


#ifndef SCAN_VEC2_TSUB_DEL_CUH
#include "../util/scan_vec2_tsub_del.cuh"
#endif//SCAN_VEC2_TSUB_DEL_CUH


#ifndef SCAN_VEC2_TSUB_INC_CUH
#define SCAN_VEC2_TSUB_INC_CUH


/**# implement shared memory version only */
/* #   if  defined(USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_INC) && (GPUGEN < 30) */
/* #undef          USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_INC */
/* #endif//defined(USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_INC) && (GPUGEN < 30) */
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_INC
#undef  USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_INC
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_INC


#   if  TSUB_SCAN_VEC2_INC ==  1
#define PREFIX_SUM_VEC2_TSUB prefixSumVec2Tsub01
#define  TOTAL_SUM_VEC2_TSUB  totalSumVec2Tsub01
#endif//TSUB_SCAN_VEC2_INC ==  1

#   if  TSUB_SCAN_VEC2_INC ==  2
#define PREFIX_SUM_VEC2_TSUB prefixSumVec2Tsub02
#define  TOTAL_SUM_VEC2_TSUB  totalSumVec2Tsub02
#endif//TSUB_SCAN_VEC2_INC ==  2

#   if  TSUB_SCAN_VEC2_INC ==  4
#define PREFIX_SUM_VEC2_TSUB prefixSumVec2Tsub04
#define  TOTAL_SUM_VEC2_TSUB  totalSumVec2Tsub04
#endif//TSUB_SCAN_VEC2_INC ==  4

#   if  TSUB_SCAN_VEC2_INC ==  8
#define PREFIX_SUM_VEC2_TSUB prefixSumVec2Tsub08
#define  TOTAL_SUM_VEC2_TSUB  totalSumVec2Tsub08
#endif//TSUB_SCAN_VEC2_INC ==  8

#   if  TSUB_SCAN_VEC2_INC == 16
#define PREFIX_SUM_VEC2_TSUB prefixSumVec2Tsub16
#define  TOTAL_SUM_VEC2_TSUB  totalSumVec2Tsub16
#endif//TSUB_SCAN_VEC2_INC == 16

#   if  TSUB_SCAN_VEC2_INC == 32
#define PREFIX_SUM_VEC2_TSUB prefixSumVec2Tsub32
#define  TOTAL_SUM_VEC2_TSUB  totalSumVec2Tsub32
#endif//TSUB_SCAN_VEC2_INC == 32


#ifdef  SCAN_VEC2_TSUB_DEL_CUH
#undef  SCAN_VEC2_TSUB_DEL_CUH
#endif//SCAN_VEC2_TSUB_DEL_CUH


#endif//SCAN_VEC2_TSUB_INC_CUH
