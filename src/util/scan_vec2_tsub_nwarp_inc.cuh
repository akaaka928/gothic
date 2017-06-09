/**
 * @file scan_vec2_tsub_nwarp_inc.cuh
 *
 * @brief Header file for parallel prefix sum library for 2-components vector on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/04/06 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */


#ifndef SCAN_VEC2_TSUB_NWARP_DEL_CUH
#include "../util/scan_vec2_tsub_nwarp_del.cuh"
#endif//SCAN_VEC2_TSUB_NWARP_DEL_CUH


#ifndef SCAN_VEC2_TSUB_NWARP_INC_CUH
#define SCAN_VEC2_TSUB_NWARP_INC_CUH


/**# implement shared memory version only */
/* #   if  defined(USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_NWARP_INC) && (GPUGEN < 30) */
/* #undef          USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_NWARP_INC */
/* #endif//defined(USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_NWARP_INC) && (GPUGEN < 30) */
#ifdef  USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_NWARP_INC
#undef  USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_NWARP_INC
#endif//USE_WARP_SHUFFLE_FUNC_SCAN_VEC2_TSUB_NWARP_INC


#   if  (TSUB_TN_SCAN_VEC2_INC ==  1) && (NWARP_TN_SCAN_VEC2_INC ==  1)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub01Nwarp01
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub01Nwarp01
#endif//(TSUB_TN_SCAN_VEC2_INC ==  1) && (NWARP_TN_SCAN_VEC2_INC ==  1)


#   if  (TSUB_TN_SCAN_VEC2_INC ==  2) && (NWARP_TN_SCAN_VEC2_INC ==  1)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub02Nwarp01
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub02Nwarp01
#endif//(TSUB_TN_SCAN_VEC2_INC ==  2) && (NWARP_TN_SCAN_VEC2_INC ==  1)

#   if  (TSUB_TN_SCAN_VEC2_INC ==  2) && (NWARP_TN_SCAN_VEC2_INC ==  2)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub02Nwarp02
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub02Nwarp02
#endif//(TSUB_TN_SCAN_VEC2_INC ==  2) && (NWARP_TN_SCAN_VEC2_INC ==  2)


#   if  (TSUB_TN_SCAN_VEC2_INC ==  4) && (NWARP_TN_SCAN_VEC2_INC ==  1)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub04Nwarp01
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub04Nwarp01
#endif//(TSUB_TN_SCAN_VEC2_INC ==  4) && (NWARP_TN_SCAN_VEC2_INC ==  1)

#   if  (TSUB_TN_SCAN_VEC2_INC ==  4) && (NWARP_TN_SCAN_VEC2_INC ==  2)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub04Nwarp02
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub04Nwarp02
#endif//(TSUB_TN_SCAN_VEC2_INC ==  4) && (NWARP_TN_SCAN_VEC2_INC ==  2)

#   if  (TSUB_TN_SCAN_VEC2_INC ==  4) && (NWARP_TN_SCAN_VEC2_INC ==  4)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub04Nwarp04
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub04Nwarp04
#endif//(TSUB_TN_SCAN_VEC2_INC ==  4) && (NWARP_TN_SCAN_VEC2_INC ==  4)


#   if  (TSUB_TN_SCAN_VEC2_INC ==  8) && (NWARP_TN_SCAN_VEC2_INC ==  1)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub08Nwarp01
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub08Nwarp01
#endif//(TSUB_TN_SCAN_VEC2_INC ==  8) && (NWARP_TN_SCAN_VEC2_INC ==  1)

#   if  (TSUB_TN_SCAN_VEC2_INC ==  8) && (NWARP_TN_SCAN_VEC2_INC ==  2)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub08Nwarp02
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub08Nwarp02
#endif//(TSUB_TN_SCAN_VEC2_INC ==  8) && (NWARP_TN_SCAN_VEC2_INC ==  2)

#   if  (TSUB_TN_SCAN_VEC2_INC ==  8) && (NWARP_TN_SCAN_VEC2_INC ==  4)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub08Nwarp04
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub08Nwarp04
#endif//(TSUB_TN_SCAN_VEC2_INC ==  8) && (NWARP_TN_SCAN_VEC2_INC ==  4)

#   if  (TSUB_TN_SCAN_VEC2_INC ==  8) && (NWARP_TN_SCAN_VEC2_INC ==  8)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub08Nwarp08
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub08Nwarp08
#endif//(TSUB_TN_SCAN_VEC2_INC ==  8) && (NWARP_TN_SCAN_VEC2_INC ==  8)


#   if  (TSUB_TN_SCAN_VEC2_INC == 16) && (NWARP_TN_SCAN_VEC2_INC ==  1)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub16Nwarp01
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub16Nwarp01
#endif//(TSUB_TN_SCAN_VEC2_INC == 16) && (NWARP_TN_SCAN_VEC2_INC ==  1)

#   if  (TSUB_TN_SCAN_VEC2_INC == 16) && (NWARP_TN_SCAN_VEC2_INC ==  2)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub16Nwarp02
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub16Nwarp02
#endif//(TSUB_TN_SCAN_VEC2_INC == 16) && (NWARP_TN_SCAN_VEC2_INC ==  2)

#   if  (TSUB_TN_SCAN_VEC2_INC == 16) && (NWARP_TN_SCAN_VEC2_INC ==  4)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub16Nwarp04
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub16Nwarp04
#endif//(TSUB_TN_SCAN_VEC2_INC == 16) && (NWARP_TN_SCAN_VEC2_INC ==  4)

#   if  (TSUB_TN_SCAN_VEC2_INC == 16) && (NWARP_TN_SCAN_VEC2_INC ==  8)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub16Nwarp08
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub16Nwarp08
#endif//(TSUB_TN_SCAN_VEC2_INC == 16) && (NWARP_TN_SCAN_VEC2_INC ==  8)

#   if  (TSUB_TN_SCAN_VEC2_INC == 16) && (NWARP_TN_SCAN_VEC2_INC == 16)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub16Nwarp16
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub16Nwarp16
#endif//(TSUB_TN_SCAN_VEC2_INC == 16) && (NWARP_TN_SCAN_VEC2_INC == 16)


#   if  (TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC ==  1)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub32Nwarp01
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub32Nwarp01
#endif//(TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC ==  1)

#   if  (TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC ==  2)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub32Nwarp02
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub32Nwarp02
#endif//(TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC ==  2)

#   if  (TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC ==  4)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub32Nwarp04
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub32Nwarp04
#endif//(TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC ==  4)

#   if  (TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC ==  8)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub32Nwarp08
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub32Nwarp08
#endif//(TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC ==  8)

#   if  (TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC == 16)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub32Nwarp16
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub32Nwarp16
#endif//(TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC == 16)

#   if  (TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC == 32)
#define PREFIX_SUM_VEC2_TSUB_NWARP prefixSumVec2Tsub32Nwarp32
#define  TOTAL_SUM_VEC2_TSUB_NWARP  totalSumVec2Tsub32Nwarp32
#endif//(TSUB_TN_SCAN_VEC2_INC == 32) && (NWARP_TN_SCAN_VEC2_INC == 32)


#ifdef  SCAN_VEC2_TSUB_NWARP_DEL_CUH
#undef  SCAN_VEC2_TSUB_NWARP_DEL_CUH
#endif//SCAN_VEC2_TSUB_NWARP_DEL_CUH


#endif//SCAN_VEC2_TSUB_NWARP_INC_CUH
