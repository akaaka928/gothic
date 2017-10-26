/**
 * @file compare_vec3_inc.cuh
 *
 * @brief Header file for comparing values in 3-components vector on GPU
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


#ifndef COMPARE_VEC3_DEL_CUH
#include "../util/compare_vec3_del.cuh"
#endif//COMPARE_VEC3_DEL_CUH


#ifndef COMPARE_VEC3_INC_CUH
#define COMPARE_VEC3_INC_CUH


/**# implement shared memory version only */
/* #   if  defined(USE_WARP_SHUFFLE_FUNC_COMPARE_VEC3_INC) && (GPUGEN < 30) */
/* #undef          USE_WARP_SHUFFLE_FUNC_COMPARE_VEC3_INC */
/* #endif//defined(USE_WARP_SHUFFLE_FUNC_COMPARE_VEC3_INC) && (GPUGEN < 30) */
#ifdef  USE_WARP_SHUFFLE_FUNC_COMPARE_VEC3_INC
#undef  USE_WARP_SHUFFLE_FUNC_COMPARE_VEC3_INC
#endif//USE_WARP_SHUFFLE_FUNC_COMPARE_VEC3_INC


#   if  NTHREADS_COMPARE_VEC3_INC ==   32
#define GET_VEC3_MIN_BLCK getVec3MinBlck0032
#define GET_VEC3_MIN_GRID getVec3MinGrid0032
#define GET_VEC3_MAX_BLCK getVec3MaxBlck0032
#define GET_VEC3_MAX_GRID getVec3MaxGrid0032
#define GET_VEC3_MIN_MAX_BLCK getVec3MinMaxBlck0032
#define GET_VEC3_MIN_MAX_GRID getVec3MinMaxGrid0032
#endif//NTHREADS_COMPARE_VEC3_INC ==   32

#   if  NTHREADS_COMPARE_VEC3_INC ==   64
#define GET_VEC3_MIN_BLCK getVec3MinBlck0064
#define GET_VEC3_MIN_GRID getVec3MinGrid0064
#define GET_VEC3_MAX_BLCK getVec3MaxBlck0064
#define GET_VEC3_MAX_GRID getVec3MaxGrid0064
#define GET_VEC3_MIN_MAX_BLCK getVec3MinMaxBlck0064
#define GET_VEC3_MIN_MAX_GRID getVec3MinMaxGrid0064
#endif//NTHREADS_COMPARE_VEC3_INC ==   64

#   if  NTHREADS_COMPARE_VEC3_INC ==  128
#define GET_VEC3_MIN_BLCK getVec3MinBlck0128
#define GET_VEC3_MIN_GRID getVec3MinGrid0128
#define GET_VEC3_MAX_BLCK getVec3MaxBlck0128
#define GET_VEC3_MAX_GRID getVec3MaxGrid0128
#define GET_VEC3_MIN_MAX_BLCK getVec3MinMaxBlck0128
#define GET_VEC3_MIN_MAX_GRID getVec3MinMaxGrid0128
#endif//NTHREADS_COMPARE_VEC3_INC ==  128

#   if  NTHREADS_COMPARE_VEC3_INC ==  256
#define GET_VEC3_MIN_BLCK getVec3MinBlck0256
#define GET_VEC3_MIN_GRID getVec3MinGrid0256
#define GET_VEC3_MAX_BLCK getVec3MaxBlck0256
#define GET_VEC3_MAX_GRID getVec3MaxGrid0256
#define GET_VEC3_MIN_MAX_BLCK getVec3MinMaxBlck0256
#define GET_VEC3_MIN_MAX_GRID getVec3MinMaxGrid0256
#endif//NTHREADS_COMPARE_VEC3_INC ==  256

#   if  NTHREADS_COMPARE_VEC3_INC ==  512
#define GET_VEC3_MIN_BLCK getVec3MinBlck0512
#define GET_VEC3_MIN_GRID getVec3MinGrid0512
#define GET_VEC3_MAX_BLCK getVec3MaxBlck0512
#define GET_VEC3_MAX_GRID getVec3MaxGrid0512
#define GET_VEC3_MIN_MAX_BLCK getVec3MinMaxBlck0512
#define GET_VEC3_MIN_MAX_GRID getVec3MinMaxGrid0512
#endif//NTHREADS_COMPARE_VEC3_INC ==  512

#   if  NTHREADS_COMPARE_VEC3_INC == 1024
#define GET_VEC3_MIN_BLCK getVec3MinBlck1024
#define GET_VEC3_MIN_GRID getVec3MinGrid1024
#define GET_VEC3_MAX_BLCK getVec3MaxBlck1024
#define GET_VEC3_MAX_GRID getVec3MaxGrid1024
#define GET_VEC3_MIN_MAX_BLCK getVec3MinMaxBlck1024
#define GET_VEC3_MIN_MAX_GRID getVec3MinMaxGrid1024
#endif//NTHREADS_COMPARE_VEC3_INC == 1024


#ifdef  COMPARE_VEC3_DEL_CUH
#undef  COMPARE_VEC3_DEL_CUH
#endif//COMPARE_VEC3_DEL_CUH


#endif//COMPARE_VEC3_INC_CUH
