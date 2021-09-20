/**
 * @file scan_vec3_inc.cuh
 *
 * @brief Header file for parallel prefix sum library for 4-components vector on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/12/28 (Fri)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */


#ifndef SCAN_VEC3_DEL_CUH
#include "../util/scan_vec3_del.cuh"
#endif//SCAN_VEC3_DEL_CUH


#ifndef SCAN_VEC3_INC_CUH
#define SCAN_VEC3_INC_CUH


#   if  NTHREADS_SCAN_VEC3_INC ==   32
#define PREFIX_SUM_VEC3_BLCK prefixSumVec3Blck0032
#define PREFIX_SUM_VEC3_GRID prefixSumVec3Grid0032
#define PREFIX_SUM_VEC3_GRID_WITH_PARTITION prefixSumVec3GridPartition0032
#define  TOTAL_SUM_VEC3_BLCK  totalSumVec3Blck0032
#define  TOTAL_SUM_VEC3_GRID  totalSumVec3Grid0032
#endif//NTHREADS_SCAN_VEC3_INC ==   32

#   if  NTHREADS_SCAN_VEC3_INC ==   64
#define PREFIX_SUM_VEC3_BLCK prefixSumVec3Blck0064
#define PREFIX_SUM_VEC3_GRID prefixSumVec3Grid0064
#define PREFIX_SUM_VEC3_GRID_WITH_PARTITION prefixSumVec3GridPartition0064
#define  TOTAL_SUM_VEC3_BLCK  totalSumVec3Blck0064
#define  TOTAL_SUM_VEC3_GRID  totalSumVec3Grid0064
#endif//NTHREADS_SCAN_VEC3_INC ==   64

#   if  NTHREADS_SCAN_VEC3_INC ==  128
#define PREFIX_SUM_VEC3_BLCK prefixSumVec3Blck0128
#define PREFIX_SUM_VEC3_GRID prefixSumVec3Grid0128
#define PREFIX_SUM_VEC3_GRID_WITH_PARTITION prefixSumVec3GridPartition0128
#define  TOTAL_SUM_VEC3_BLCK  totalSumVec3Blck0128
#define  TOTAL_SUM_VEC3_GRID  totalSumVec3Grid0128
#endif//NTHREADS_SCAN_VEC3_INC ==  128

#   if  NTHREADS_SCAN_VEC3_INC ==  256
#define PREFIX_SUM_VEC3_BLCK prefixSumVec3Blck0256
#define PREFIX_SUM_VEC3_GRID prefixSumVec3Grid0256
#define PREFIX_SUM_VEC3_GRID_WITH_PARTITION prefixSumVec3GridPartition0256
#define  TOTAL_SUM_VEC3_BLCK  totalSumVec3Blck0256
#define  TOTAL_SUM_VEC3_GRID  totalSumVec3Grid0256
#endif//NTHREADS_SCAN_VEC3_INC ==  256

#   if  NTHREADS_SCAN_VEC3_INC ==  512
#define PREFIX_SUM_VEC3_BLCK prefixSumVec3Blck0512
#define PREFIX_SUM_VEC3_GRID prefixSumVec3Grid0512
#define PREFIX_SUM_VEC3_GRID_WITH_PARTITION prefixSumVec3GridPartition0512
#define  TOTAL_SUM_VEC3_BLCK  totalSumVec3Blck0512
#define  TOTAL_SUM_VEC3_GRID  totalSumVec3Grid0512
#endif//NTHREADS_SCAN_VEC3_INC ==  512

#   if  NTHREADS_SCAN_VEC3_INC == 1024
#define PREFIX_SUM_VEC3_BLCK prefixSumVec3Blck1024
#define PREFIX_SUM_VEC3_GRID prefixSumVec3Grid1024
#define PREFIX_SUM_VEC3_GRID_WITH_PARTITION prefixSumVec3GridPartition1024
#define  TOTAL_SUM_VEC3_BLCK  totalSumVec3Blck1024
#define  TOTAL_SUM_VEC3_GRID  totalSumVec3Grid1024
#endif//NTHREADS_SCAN_VEC3_INC == 1024


#ifdef  SCAN_VEC3_DEL_CUH
#undef  SCAN_VEC3_DEL_CUH
#endif//SCAN_VEC3_DEL_CUH


#endif//SCAN_VEC3_INC_CUH
