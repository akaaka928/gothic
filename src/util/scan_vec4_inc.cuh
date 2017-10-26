/**
 * @file scan_vec4_inc.cuh
 *
 * @brief Header file for parallel prefix sum library for 4-components vector on GPU
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


#ifndef SCAN_VEC4_DEL_CUH
#include "../util/scan_vec4_del.cuh"
#endif//SCAN_VEC4_DEL_CUH


#ifndef SCAN_VEC4_INC_CUH
#define SCAN_VEC4_INC_CUH


#   if  NTHREADS_SCAN_VEC4_INC ==   32
#define PREFIX_SUM_VEC4_BLCK prefixSumVec4Blck0032
#define PREFIX_SUM_VEC4_GRID prefixSumVec4Grid0032
#define PREFIX_SUM_VEC4_GRID_WITH_PARTITION prefixSumVec4GridPartition0032
#define  TOTAL_SUM_VEC4_BLCK  totalSumVec4Blck0032
#define  TOTAL_SUM_VEC4_GRID  totalSumVec4Grid0032
#endif//NTHREADS_SCAN_VEC4_INC ==   32

#   if  NTHREADS_SCAN_VEC4_INC ==   64
#define PREFIX_SUM_VEC4_BLCK prefixSumVec4Blck0064
#define PREFIX_SUM_VEC4_GRID prefixSumVec4Grid0064
#define PREFIX_SUM_VEC4_GRID_WITH_PARTITION prefixSumVec4GridPartition0064
#define  TOTAL_SUM_VEC4_BLCK  totalSumVec4Blck0064
#define  TOTAL_SUM_VEC4_GRID  totalSumVec4Grid0064
#endif//NTHREADS_SCAN_VEC4_INC ==   64

#   if  NTHREADS_SCAN_VEC4_INC ==  128
#define PREFIX_SUM_VEC4_BLCK prefixSumVec4Blck0128
#define PREFIX_SUM_VEC4_GRID prefixSumVec4Grid0128
#define PREFIX_SUM_VEC4_GRID_WITH_PARTITION prefixSumVec4GridPartition0128
#define  TOTAL_SUM_VEC4_BLCK  totalSumVec4Blck0128
#define  TOTAL_SUM_VEC4_GRID  totalSumVec4Grid0128
#endif//NTHREADS_SCAN_VEC4_INC ==  128

#   if  NTHREADS_SCAN_VEC4_INC ==  256
#define PREFIX_SUM_VEC4_BLCK prefixSumVec4Blck0256
#define PREFIX_SUM_VEC4_GRID prefixSumVec4Grid0256
#define PREFIX_SUM_VEC4_GRID_WITH_PARTITION prefixSumVec4GridPartition0256
#define  TOTAL_SUM_VEC4_BLCK  totalSumVec4Blck0256
#define  TOTAL_SUM_VEC4_GRID  totalSumVec4Grid0256
#endif//NTHREADS_SCAN_VEC4_INC ==  256

#   if  NTHREADS_SCAN_VEC4_INC ==  512
#define PREFIX_SUM_VEC4_BLCK prefixSumVec4Blck0512
#define PREFIX_SUM_VEC4_GRID prefixSumVec4Grid0512
#define PREFIX_SUM_VEC4_GRID_WITH_PARTITION prefixSumVec4GridPartition0512
#define  TOTAL_SUM_VEC4_BLCK  totalSumVec4Blck0512
#define  TOTAL_SUM_VEC4_GRID  totalSumVec4Grid0512
#endif//NTHREADS_SCAN_VEC4_INC ==  512

#   if  NTHREADS_SCAN_VEC4_INC == 1024
#define PREFIX_SUM_VEC4_BLCK prefixSumVec4Blck1024
#define PREFIX_SUM_VEC4_GRID prefixSumVec4Grid1024
#define PREFIX_SUM_VEC4_GRID_WITH_PARTITION prefixSumVec4GridPartition1024
#define  TOTAL_SUM_VEC4_BLCK  totalSumVec4Blck1024
#define  TOTAL_SUM_VEC4_GRID  totalSumVec4Grid1024
#endif//NTHREADS_SCAN_VEC4_INC == 1024


#ifdef  SCAN_VEC4_DEL_CUH
#undef  SCAN_VEC4_DEL_CUH
#endif//SCAN_VEC4_DEL_CUH


#endif//SCAN_VEC4_INC_CUH
