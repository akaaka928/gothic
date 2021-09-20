/**
 * @file compare_vec2_inc.cuh
 *
 * @brief Header file for comparing values in 2-components vector on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/05/23 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */


#ifndef COMPARE_VEC2_DEL_CUH
#include "../util/compare_vec2_del.cuh"
#endif//COMPARE_VEC2_DEL_CUH


#ifndef COMPARE_VEC2_INC_CUH
#define COMPARE_VEC2_INC_CUH


/**# implement shared memory version only */
/* #   if  defined(USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC) && (GPUVER < 30) */
/* #undef          USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC */
/* #endif//defined(USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC) && (GPUVER < 30) */
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC
#undef  USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_VEC2_INC


#   if  NTHREADS_COMPARE_VEC2_INC ==   32
#define GET_VEC2_MIN_BLCK getVec2MinBlck0032
#define GET_VEC2_MIN_GRID getVec2MinGrid0032
#define GET_VEC2_MAX_BLCK getVec2MaxBlck0032
#define GET_VEC2_MAX_GRID getVec2MaxGrid0032
#define GET_VEC2_MIN_MAX_BLCK getVec2MinMaxBlck0032
#define GET_VEC2_MIN_MAX_GRID getVec2MinMaxGrid0032
#endif//NTHREADS_COMPARE_VEC2_INC ==   32

#   if  NTHREADS_COMPARE_VEC2_INC ==   64
#define GET_VEC2_MIN_BLCK getVec2MinBlck0064
#define GET_VEC2_MIN_GRID getVec2MinGrid0064
#define GET_VEC2_MAX_BLCK getVec2MaxBlck0064
#define GET_VEC2_MAX_GRID getVec2MaxGrid0064
#define GET_VEC2_MIN_MAX_BLCK getVec2MinMaxBlck0064
#define GET_VEC2_MIN_MAX_GRID getVec2MinMaxGrid0064
#endif//NTHREADS_COMPARE_VEC2_INC ==   64

#   if  NTHREADS_COMPARE_VEC2_INC ==  128
#define GET_VEC2_MIN_BLCK getVec2MinBlck0128
#define GET_VEC2_MIN_GRID getVec2MinGrid0128
#define GET_VEC2_MAX_BLCK getVec2MaxBlck0128
#define GET_VEC2_MAX_GRID getVec2MaxGrid0128
#define GET_VEC2_MIN_MAX_BLCK getVec2MinMaxBlck0128
#define GET_VEC2_MIN_MAX_GRID getVec2MinMaxGrid0128
#endif//NTHREADS_COMPARE_VEC2_INC ==  128

#   if  NTHREADS_COMPARE_VEC2_INC ==  256
#define GET_VEC2_MIN_BLCK getVec2MinBlck0256
#define GET_VEC2_MIN_GRID getVec2MinGrid0256
#define GET_VEC2_MAX_BLCK getVec2MaxBlck0256
#define GET_VEC2_MAX_GRID getVec2MaxGrid0256
#define GET_VEC2_MIN_MAX_BLCK getVec2MinMaxBlck0256
#define GET_VEC2_MIN_MAX_GRID getVec2MinMaxGrid0256
#endif//NTHREADS_COMPARE_VEC2_INC ==  256

#   if  NTHREADS_COMPARE_VEC2_INC ==  512
#define GET_VEC2_MIN_BLCK getVec2MinBlck0512
#define GET_VEC2_MIN_GRID getVec2MinGrid0512
#define GET_VEC2_MAX_BLCK getVec2MaxBlck0512
#define GET_VEC2_MAX_GRID getVec2MaxGrid0512
#define GET_VEC2_MIN_MAX_BLCK getVec2MinMaxBlck0512
#define GET_VEC2_MIN_MAX_GRID getVec2MinMaxGrid0512
#endif//NTHREADS_COMPARE_VEC2_INC ==  512

#   if  NTHREADS_COMPARE_VEC2_INC == 1024
#define GET_VEC2_MIN_BLCK getVec2MinBlck1024
#define GET_VEC2_MIN_GRID getVec2MinGrid1024
#define GET_VEC2_MAX_BLCK getVec2MaxBlck1024
#define GET_VEC2_MAX_GRID getVec2MaxGrid1024
#define GET_VEC2_MIN_MAX_BLCK getVec2MinMaxBlck1024
#define GET_VEC2_MIN_MAX_GRID getVec2MinMaxGrid1024
#endif//NTHREADS_COMPARE_VEC2_INC == 1024


#ifdef  COMPARE_VEC2_DEL_CUH
#undef  COMPARE_VEC2_DEL_CUH
#endif//COMPARE_VEC2_DEL_CUH


#endif//COMPARE_VEC2_INC_CUH
