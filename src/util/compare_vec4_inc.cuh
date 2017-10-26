/**
 * @file compare_vec4_inc.cuh
 *
 * @brief Header file for comparing values in 4-components vector on GPU
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


#ifndef COMPARE_VEC4_DEL_CUH
#include "../util/compare_vec4_del.cuh"
#endif//COMPARE_VEC4_DEL_CUH


#ifndef COMPARE_VEC4_INC_CUH
#define COMPARE_VEC4_INC_CUH


/**# implement shared memory version only */
/* #   if  defined(USE_WARP_SHUFFLE_FUNC_COMPARE_VEC4_INC) && (GPUGEN < 30) */
/* #undef          USE_WARP_SHUFFLE_FUNC_COMPARE_VEC4_INC */
/* #endif//defined(USE_WARP_SHUFFLE_FUNC_COMPARE_VEC4_INC) && (GPUGEN < 30) */
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_VEC4_INC
#undef  USE_WARP_SHUFFLE_FUNC_COMPARE_VEC4_INC
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_VEC4_INC


#   if  NTHREADS_COMPARE_VEC4_INC ==   32
#define GET_VEC4_MIN_BLCK getVec4MinBlck0032
#define GET_VEC4_MIN_GRID getVec4MinGrid0032
#define GET_VEC4_MAX_BLCK getVec4MaxBlck0032
#define GET_VEC4_MAX_GRID getVec4MaxGrid0032
#define GET_VEC4_MIN_MAX_BLCK getVec4MinMaxBlck0032
#define GET_VEC4_MIN_MAX_GRID getVec4MinMaxGrid0032
#endif//NTHREADS_COMPARE_VEC4_INC ==   32

#   if  NTHREADS_COMPARE_VEC4_INC ==   64
#define GET_VEC4_MIN_BLCK getVec4MinBlck0064
#define GET_VEC4_MIN_GRID getVec4MinGrid0064
#define GET_VEC4_MAX_BLCK getVec4MaxBlck0064
#define GET_VEC4_MAX_GRID getVec4MaxGrid0064
#define GET_VEC4_MIN_MAX_BLCK getVec4MinMaxBlck0064
#define GET_VEC4_MIN_MAX_GRID getVec4MinMaxGrid0064
#endif//NTHREADS_COMPARE_VEC4_INC ==   64

#   if  NTHREADS_COMPARE_VEC4_INC ==  128
#define GET_VEC4_MIN_BLCK getVec4MinBlck0128
#define GET_VEC4_MIN_GRID getVec4MinGrid0128
#define GET_VEC4_MAX_BLCK getVec4MaxBlck0128
#define GET_VEC4_MAX_GRID getVec4MaxGrid0128
#define GET_VEC4_MIN_MAX_BLCK getVec4MinMaxBlck0128
#define GET_VEC4_MIN_MAX_GRID getVec4MinMaxGrid0128
#endif//NTHREADS_COMPARE_VEC4_INC ==  128

#   if  NTHREADS_COMPARE_VEC4_INC ==  256
#define GET_VEC4_MIN_BLCK getVec4MinBlck0256
#define GET_VEC4_MIN_GRID getVec4MinGrid0256
#define GET_VEC4_MAX_BLCK getVec4MaxBlck0256
#define GET_VEC4_MAX_GRID getVec4MaxGrid0256
#define GET_VEC4_MIN_MAX_BLCK getVec4MinMaxBlck0256
#define GET_VEC4_MIN_MAX_GRID getVec4MinMaxGrid0256
#endif//NTHREADS_COMPARE_VEC4_INC ==  256

#   if  NTHREADS_COMPARE_VEC4_INC ==  512
#define GET_VEC4_MIN_BLCK getVec4MinBlck0512
#define GET_VEC4_MIN_GRID getVec4MinGrid0512
#define GET_VEC4_MAX_BLCK getVec4MaxBlck0512
#define GET_VEC4_MAX_GRID getVec4MaxGrid0512
#define GET_VEC4_MIN_MAX_BLCK getVec4MinMaxBlck0512
#define GET_VEC4_MIN_MAX_GRID getVec4MinMaxGrid0512
#endif//NTHREADS_COMPARE_VEC4_INC ==  512

#   if  NTHREADS_COMPARE_VEC4_INC == 1024
#define GET_VEC4_MIN_BLCK getVec4MinBlck1024
#define GET_VEC4_MIN_GRID getVec4MinGrid1024
#define GET_VEC4_MAX_BLCK getVec4MaxBlck1024
#define GET_VEC4_MAX_GRID getVec4MaxGrid1024
#define GET_VEC4_MIN_MAX_BLCK getVec4MinMaxBlck1024
#define GET_VEC4_MIN_MAX_GRID getVec4MinMaxGrid1024
#endif//NTHREADS_COMPARE_VEC4_INC == 1024


#ifdef  COMPARE_VEC4_DEL_CUH
#undef  COMPARE_VEC4_DEL_CUH
#endif//COMPARE_VEC4_DEL_CUH


#endif//COMPARE_VEC4_INC_CUH
