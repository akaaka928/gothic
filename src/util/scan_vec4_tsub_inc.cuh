/**
 * @file scan_vec4_tsub_inc.cuh
 *
 * @brief Header file for parallel prefix sum library for 4-components vector on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/04/05 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */


#ifndef SCAN_VEC4_TSUB_DEL_CUH
#include "../util/scan_vec4_tsub_del.cuh"
#endif//SCAN_VEC4_TSUB_DEL_CUH


#ifndef SCAN_VEC4_TSUB_INC_CUH
#define SCAN_VEC4_TSUB_INC_CUH


/**# implement shared memory version only */
/* #   if  defined(USE_WARP_SHUFFLE_FUNC_SCAN_VEC4_TSUB_INC) && (GPUGEN < 30) */
/* #undef          USE_WARP_SHUFFLE_FUNC_SCAN_VEC4_TSUB_INC */
/* #endif//defined(USE_WARP_SHUFFLE_FUNC_SCAN_VEC4_TSUB_INC) && (GPUGEN < 30) */
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_VEC4_TSUB_INC
#undef  USE_WARP_SHUFFLE_FUNC_SCAN_VEC4_TSUB_INC
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_VEC4_TSUB_INC


#   if  TSUB_SCAN_VEC4_INC ==  1
#define PREFIX_SUM_VEC4_TSUB prefixSumVec4Tsub01
#define  TOTAL_SUM_VEC4_TSUB  totalSumVec4Tsub01
#endif//TSUB_SCAN_VEC4_INC ==  1

#   if  TSUB_SCAN_VEC4_INC ==  2
#define PREFIX_SUM_VEC4_TSUB prefixSumVec4Tsub02
#define  TOTAL_SUM_VEC4_TSUB  totalSumVec4Tsub02
#endif//TSUB_SCAN_VEC4_INC ==  2

#   if  TSUB_SCAN_VEC4_INC ==  4
#define PREFIX_SUM_VEC4_TSUB prefixSumVec4Tsub04
#define  TOTAL_SUM_VEC4_TSUB  totalSumVec4Tsub04
#endif//TSUB_SCAN_VEC4_INC ==  4

#   if  TSUB_SCAN_VEC4_INC ==  8
#define PREFIX_SUM_VEC4_TSUB prefixSumVec4Tsub08
#define  TOTAL_SUM_VEC4_TSUB  totalSumVec4Tsub08
#endif//TSUB_SCAN_VEC4_INC ==  8

#   if  TSUB_SCAN_VEC4_INC == 16
#define PREFIX_SUM_VEC4_TSUB prefixSumVec4Tsub16
#define  TOTAL_SUM_VEC4_TSUB  totalSumVec4Tsub16
#endif//TSUB_SCAN_VEC4_INC == 16

#   if  TSUB_SCAN_VEC4_INC == 32
#define PREFIX_SUM_VEC4_TSUB prefixSumVec4Tsub32
#define  TOTAL_SUM_VEC4_TSUB  totalSumVec4Tsub32
#endif//TSUB_SCAN_VEC4_INC == 32


#ifdef  SCAN_VEC4_TSUB_DEL_CUH
#undef  SCAN_VEC4_TSUB_DEL_CUH
#endif//SCAN_VEC4_TSUB_DEL_CUH


#endif//SCAN_VEC4_TSUB_INC_CUH
