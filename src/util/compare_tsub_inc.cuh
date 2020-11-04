/**
 * @file compare_tsub_inc.cuh
 *
 * @brief Header file for comparing values on GPU
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


#ifndef COMPARE_TSUB_DEL_CUH
#include "../util/compare_tsub_del.cuh"
#endif//COMPARE_TSUB_DEL_CUH


#ifndef COMPARE_TSUB_INC_CUH
#define COMPARE_TSUB_INC_CUH


#   if  defined(USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC) && (GPUVER < 30)
#undef          USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC
#endif//defined(USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_INC) && (GPUVER < 30)

#   if  defined(USE_WARP_REDUCE_FUNCTIONS_COMPARE_TSUB_INC) && !defined(ENABLE_WARP_REDUCE_FUNCTIONS)
#undef          USE_WARP_REDUCE_FUNCTIONS_COMPARE_TSUB_INC
#endif//defined(USE_WARP_REDUCE_FUNCTIONS_COMPARE_TSUB_INC) && !defined(ENABLE_WARP_REDUCE_FUNCTIONS)


#   if  TSUB_COMPARE_INC ==  1
#define GET_MIN_TSUB       getMinTsub01
#define GET_MAX_TSUB       getMaxTsub01
#define GET_MINLOC_TSUB getMinLocTsub01
#define GET_MAXLOC_TSUB getMaxLocTsub01
#endif//TSUB_COMPARE_INC ==  1

#   if  TSUB_COMPARE_INC ==  2
#define GET_MIN_TSUB       getMinTsub02
#define GET_MAX_TSUB       getMaxTsub02
#define GET_MINLOC_TSUB getMinLocTsub02
#define GET_MAXLOC_TSUB getMaxLocTsub02
#endif//TSUB_COMPARE_INC ==  2

#   if  TSUB_COMPARE_INC ==  4
#define GET_MIN_TSUB       getMinTsub04
#define GET_MAX_TSUB       getMaxTsub04
#define GET_MINLOC_TSUB getMinLocTsub04
#define GET_MAXLOC_TSUB getMaxLocTsub04
#endif//TSUB_COMPARE_INC ==  4

#   if  TSUB_COMPARE_INC ==  8
#define GET_MIN_TSUB       getMinTsub08
#define GET_MAX_TSUB       getMaxTsub08
#define GET_MINLOC_TSUB getMinLocTsub08
#define GET_MAXLOC_TSUB getMaxLocTsub08
#endif//TSUB_COMPARE_INC ==  8

#   if  TSUB_COMPARE_INC == 16
#define GET_MIN_TSUB       getMinTsub16
#define GET_MAX_TSUB       getMaxTsub16
#define GET_MINLOC_TSUB getMinLocTsub16
#define GET_MAXLOC_TSUB getMaxLocTsub16
#endif//TSUB_COMPARE_INC == 16

#   if  TSUB_COMPARE_INC == 32
#define GET_MIN_TSUB       getMinTsub32
#define GET_MAX_TSUB       getMaxTsub32
#define GET_MINLOC_TSUB getMinLocTsub32
#define GET_MAXLOC_TSUB getMaxLocTsub32
#endif//TSUB_COMPARE_INC == 32


#ifdef  COMPARE_TSUB_DEL_CUH
#undef  COMPARE_TSUB_DEL_CUH
#endif//COMPARE_TSUB_DEL_CUH


#endif//COMPARE_TSUB_INC_CUH
