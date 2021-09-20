/**
 * @file compare_tsub_nwarp_inc.cuh
 *
 * @brief Header file for comparing values on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2020/11/16 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */


#ifndef COMPARE_TSUB_NWARP_DEL_CUH
#include "../util/compare_tsub_nwarp_del.cuh"
#endif//COMPARE_TSUB_NWARP_DEL_CUH


#ifndef COMPARE_TSUB_NWARP_INC_CUH
#define COMPARE_TSUB_NWARP_INC_CUH


#   if  defined(USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC) && (GPUVER < 30)
#undef          USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC
#endif//defined(USE_WARP_SHUFFLE_FUNC_COMPARE_TSUB_NWARP_INC) && (GPUVER < 30)

#   if  defined(USE_WARP_REDUCE_FUNC_COMPARE_TSUB_NWARP_INC) && !defined(ENABLE_WARP_REDUCE_FUNC)
#undef          USE_WARP_REDUCE_FUNC_COMPARE_TSUB_NWARP_INC
#endif//defined(USE_WARP_REDUCE_FUNC_COMPARE_TSUB_NWARP_INC) && !defined(ENABLE_WARP_REDUCE_FUNC)


#   if  (TSUB_TN_COMPARE_INC ==  1) && (NWARP_TN_COMPARE_INC ==  1)
#define GET_MIN_TSUB_NWARP       getMinTsub01Nwarp01
#define GET_MAX_TSUB_NWARP       getMaxTsub01Nwarp01
#define GET_MINLOC_TSUB_NWARP getMinLocTsub01Nwarp01
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub01Nwarp01
#endif//(TSUB_TN_COMPARE_INC ==  1) && (NWARP_TN_COMPARE_INC ==  1)


#   if  (TSUB_TN_COMPARE_INC ==  2) && (NWARP_TN_COMPARE_INC ==  1)
#define GET_MIN_TSUB_NWARP       getMinTsub02Nwarp01
#define GET_MAX_TSUB_NWARP       getMaxTsub02Nwarp01
#define GET_MINLOC_TSUB_NWARP getMinLocTsub02Nwarp01
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub02Nwarp01
#endif//(TSUB_TN_COMPARE_INC ==  2) && (NWARP_TN_COMPARE_INC ==  1)

#   if  (TSUB_TN_COMPARE_INC ==  2) && (NWARP_TN_COMPARE_INC ==  2)
#define GET_MIN_TSUB_NWARP       getMinTsub02Nwarp02
#define GET_MAX_TSUB_NWARP       getMaxTsub02Nwarp02
#define GET_MINLOC_TSUB_NWARP getMinLocTsub02Nwarp02
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub02Nwarp02
#endif//(TSUB_TN_COMPARE_INC ==  2) && (NWARP_TN_COMPARE_INC ==  2)


#   if  (TSUB_TN_COMPARE_INC ==  4) && (NWARP_TN_COMPARE_INC ==  1)
#define GET_MIN_TSUB_NWARP       getMinTsub04Nwarp01
#define GET_MAX_TSUB_NWARP       getMaxTsub04Nwarp01
#define GET_MINLOC_TSUB_NWARP getMinLocTsub04Nwarp01
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub04Nwarp01
#endif//(TSUB_TN_COMPARE_INC ==  4) && (NWARP_TN_COMPARE_INC ==  1)

#   if  (TSUB_TN_COMPARE_INC ==  4) && (NWARP_TN_COMPARE_INC ==  2)
#define GET_MIN_TSUB_NWARP       getMinTsub04Nwarp02
#define GET_MAX_TSUB_NWARP       getMaxTsub04Nwarp02
#define GET_MINLOC_TSUB_NWARP getMinLocTsub04Nwarp02
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub04Nwarp02
#endif//(TSUB_TN_COMPARE_INC ==  4) && (NWARP_TN_COMPARE_INC ==  2)

#   if  (TSUB_TN_COMPARE_INC ==  4) && (NWARP_TN_COMPARE_INC ==  4)
#define GET_MIN_TSUB_NWARP       getMinTsub04Nwarp04
#define GET_MAX_TSUB_NWARP       getMaxTsub04Nwarp04
#define GET_MINLOC_TSUB_NWARP getMinLocTsub04Nwarp04
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub04Nwarp04
#endif//(TSUB_TN_COMPARE_INC ==  4) && (NWARP_TN_COMPARE_INC ==  4)


#   if  (TSUB_TN_COMPARE_INC ==  8) && (NWARP_TN_COMPARE_INC ==  1)
#define GET_MIN_TSUB_NWARP       getMinTsub08Nwarp01
#define GET_MAX_TSUB_NWARP       getMaxTsub08Nwarp01
#define GET_MINLOC_TSUB_NWARP getMinLocTsub08Nwarp01
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub08Nwarp01
#endif//(TSUB_TN_COMPARE_INC ==  8) && (NWARP_TN_COMPARE_INC ==  1)

#   if  (TSUB_TN_COMPARE_INC ==  8) && (NWARP_TN_COMPARE_INC ==  2)
#define GET_MIN_TSUB_NWARP       getMinTsub08Nwarp02
#define GET_MAX_TSUB_NWARP       getMaxTsub08Nwarp02
#define GET_MINLOC_TSUB_NWARP getMinLocTsub08Nwarp02
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub08Nwarp02
#endif//(TSUB_TN_COMPARE_INC ==  8) && (NWARP_TN_COMPARE_INC ==  2)

#   if  (TSUB_TN_COMPARE_INC ==  8) && (NWARP_TN_COMPARE_INC ==  4)
#define GET_MIN_TSUB_NWARP       getMinTsub08Nwarp04
#define GET_MAX_TSUB_NWARP       getMaxTsub08Nwarp04
#define GET_MINLOC_TSUB_NWARP getMinLocTsub08Nwarp04
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub08Nwarp04
#endif//(TSUB_TN_COMPARE_INC ==  8) && (NWARP_TN_COMPARE_INC ==  4)

#   if  (TSUB_TN_COMPARE_INC ==  8) && (NWARP_TN_COMPARE_INC ==  8)
#define GET_MIN_TSUB_NWARP       getMinTsub08Nwarp08
#define GET_MAX_TSUB_NWARP       getMaxTsub08Nwarp08
#define GET_MINLOC_TSUB_NWARP getMinLocTsub08Nwarp08
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub08Nwarp08
#endif//(TSUB_TN_COMPARE_INC ==  8) && (NWARP_TN_COMPARE_INC ==  8)


#   if  (TSUB_TN_COMPARE_INC == 16) && (NWARP_TN_COMPARE_INC ==  1)
#define GET_MIN_TSUB_NWARP       getMinTsub16Nwarp01
#define GET_MAX_TSUB_NWARP       getMaxTsub16Nwarp01
#define GET_MINLOC_TSUB_NWARP getMinLocTsub16Nwarp01
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub16Nwarp01
#endif//(TSUB_TN_COMPARE_INC == 16) && (NWARP_TN_COMPARE_INC ==  1)

#   if  (TSUB_TN_COMPARE_INC == 16) && (NWARP_TN_COMPARE_INC ==  2)
#define GET_MIN_TSUB_NWARP       getMinTsub16Nwarp02
#define GET_MAX_TSUB_NWARP       getMaxTsub16Nwarp02
#define GET_MINLOC_TSUB_NWARP getMinLocTsub16Nwarp02
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub16Nwarp02
#endif//(TSUB_TN_COMPARE_INC == 16) && (NWARP_TN_COMPARE_INC ==  2)

#   if  (TSUB_TN_COMPARE_INC == 16) && (NWARP_TN_COMPARE_INC ==  4)
#define GET_MIN_TSUB_NWARP       getMinTsub16Nwarp04
#define GET_MAX_TSUB_NWARP       getMaxTsub16Nwarp04
#define GET_MINLOC_TSUB_NWARP getMinLocTsub16Nwarp04
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub16Nwarp04
#endif//(TSUB_TN_COMPARE_INC == 16) && (NWARP_TN_COMPARE_INC ==  4)

#   if  (TSUB_TN_COMPARE_INC == 16) && (NWARP_TN_COMPARE_INC ==  8)
#define GET_MIN_TSUB_NWARP       getMinTsub16Nwarp08
#define GET_MAX_TSUB_NWARP       getMaxTsub16Nwarp08
#define GET_MINLOC_TSUB_NWARP getMinLocTsub16Nwarp08
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub16Nwarp08
#endif//(TSUB_TN_COMPARE_INC == 16) && (NWARP_TN_COMPARE_INC ==  8)

#   if  (TSUB_TN_COMPARE_INC == 16) && (NWARP_TN_COMPARE_INC == 16)
#define GET_MIN_TSUB_NWARP       getMinTsub16Nwarp16
#define GET_MAX_TSUB_NWARP       getMaxTsub16Nwarp16
#define GET_MINLOC_TSUB_NWARP getMinLocTsub16Nwarp16
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub16Nwarp16
#endif//(TSUB_TN_COMPARE_INC == 16) && (NWARP_TN_COMPARE_INC == 16)


#   if  (TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC ==  1)
#define GET_MIN_TSUB_NWARP       getMinTsub32Nwarp01
#define GET_MAX_TSUB_NWARP       getMaxTsub32Nwarp01
#define GET_MINLOC_TSUB_NWARP getMinLocTsub32Nwarp01
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub32Nwarp01
#endif//(TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC ==  1)

#   if  (TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC ==  2)
#define GET_MIN_TSUB_NWARP       getMinTsub32Nwarp02
#define GET_MAX_TSUB_NWARP       getMaxTsub32Nwarp02
#define GET_MINLOC_TSUB_NWARP getMinLocTsub32Nwarp02
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub32Nwarp02
#endif//(TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC ==  2)

#   if  (TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC ==  4)
#define GET_MIN_TSUB_NWARP       getMinTsub32Nwarp04
#define GET_MAX_TSUB_NWARP       getMaxTsub32Nwarp04
#define GET_MINLOC_TSUB_NWARP getMinLocTsub32Nwarp04
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub32Nwarp04
#endif//(TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC ==  4)

#   if  (TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC ==  8)
#define GET_MIN_TSUB_NWARP       getMinTsub32Nwarp08
#define GET_MAX_TSUB_NWARP       getMaxTsub32Nwarp08
#define GET_MINLOC_TSUB_NWARP getMinLocTsub32Nwarp08
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub32Nwarp08
#endif//(TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC ==  8)

#   if  (TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC == 16)
#define GET_MIN_TSUB_NWARP       getMinTsub32Nwarp16
#define GET_MAX_TSUB_NWARP       getMaxTsub32Nwarp16
#define GET_MINLOC_TSUB_NWARP getMinLocTsub32Nwarp16
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub32Nwarp16
#endif//(TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC == 16)

#   if  (TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC == 32)
#define GET_MIN_TSUB_NWARP       getMinTsub32Nwarp32
#define GET_MAX_TSUB_NWARP       getMaxTsub32Nwarp32
#define GET_MINLOC_TSUB_NWARP getMinLocTsub32Nwarp32
#define GET_MAXLOC_TSUB_NWARP getMaxLocTsub32Nwarp32
#endif//(TSUB_TN_COMPARE_INC == 32) && (NWARP_TN_COMPARE_INC == 32)


#ifdef  COMPARE_TSUB_NWARP_DEL_CUH
#undef  COMPARE_TSUB_NWARP_DEL_CUH
#endif//COMPARE_TSUB_NWARP_DEL_CUH


#endif//COMPARE_TSUB_NWARP_INC_CUH
