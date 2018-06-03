/**
 * @file scan_tsub_inc.cuh
 *
 * @brief Header file for parallel prefix sum library on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/05/31 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */


#ifndef SCAN_TSUB_DEL_CUH
#include "../util/scan_tsub_del.cuh"
#endif//SCAN_TSUB_DEL_CUH


#ifndef SCAN_TSUB_INC_CUH
#define SCAN_TSUB_INC_CUH


#   if  defined(USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC) && (GPUVER < 30)
#undef          USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC
#endif//defined(USE_WARP_SHUFFLE_FUNC_SCAN_TSUB_INC) && (GPUVER < 30)


#   if  TSUB_SCAN_INC ==  1
#define PREFIX_SUM_TSUB prefixSumTsub01
#define  TOTAL_SUM_TSUB  totalSumTsub01
#endif//TSUB_SCAN_INC ==  1

#   if  TSUB_SCAN_INC ==  2
#define PREFIX_SUM_TSUB prefixSumTsub02
#define  TOTAL_SUM_TSUB  totalSumTsub02
#endif//TSUB_SCAN_INC ==  2

#   if  TSUB_SCAN_INC ==  4
#define PREFIX_SUM_TSUB prefixSumTsub04
#define  TOTAL_SUM_TSUB  totalSumTsub04
#endif//TSUB_SCAN_INC ==  4

#   if  TSUB_SCAN_INC ==  8
#define PREFIX_SUM_TSUB prefixSumTsub08
#define  TOTAL_SUM_TSUB  totalSumTsub08
#endif//TSUB_SCAN_INC ==  8

#   if  TSUB_SCAN_INC == 16
#define PREFIX_SUM_TSUB prefixSumTsub16
#define  TOTAL_SUM_TSUB  totalSumTsub16
#endif//TSUB_SCAN_INC == 16

#   if  TSUB_SCAN_INC == 32
#define PREFIX_SUM_TSUB prefixSumTsub32
#define  TOTAL_SUM_TSUB  totalSumTsub32
#endif//TSUB_SCAN_INC == 32


#ifdef  SCAN_TSUB_DEL_CUH
#undef  SCAN_TSUB_DEL_CUH
#endif//SCAN_TSUB_DEL_CUH


#endif//SCAN_TSUB_INC_CUH
