/**
 * @file compare_inc.cuh
 *
 * @brief Header file for comparing values on GPU
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


#ifndef COMPARE_DEL_CUH
#include "../util/compare_del.cuh"
#endif//COMPARE_DEL_CUH


#ifndef COMPARE_INC_CUH
#define COMPARE_INC_CUH


#   if  defined(USE_WARP_SHUFFLE_FUNC_COMPARE_INC) && (GPUGEN < 30)
#undef          USE_WARP_SHUFFLE_FUNC_COMPARE_INC
#endif//defined(USE_WARP_SHUFFLE_FUNC_COMPARE_INC) && (GPUGEN < 30)


#   if  NTHREADS_COMPARE_INC ==   32
#define GET_MIN_BLCK getMinBlck0032
#define GET_MIN_GRID getMinGrid0032
#define GET_MAX_BLCK getMaxBlck0032
#define GET_MAX_GRID getMaxGrid0032
#define GET_MINLOC_BLCK getMinlocBlck0032
#define GET_MINLOC_GRID getMinlocGrid0032
#define GET_MAXLOC_BLCK getMaxlocBlck0032
#define GET_MAXLOC_GRID getMaxlocGrid0032
#define GET_MINLOC_MAXLOC_BLCK getMinlocMaxlocBlck0032
#define GET_MINLOC_MAXLOC_GRID getMinlocMaxlocGrid0032
#define GET_MINLOC_2VALS_GRID getMinloc2ValsGrid0032
#define GET_MINLOC_3VALS_GRID getMinloc3ValsGrid0032
#define GET_MAXLOC_2VALS_GRID getMaxloc2ValsGrid0032
#define GET_MAXLOC_3VALS_GRID getMaxloc3ValsGrid0032
#define GET_MINLOC_MAXLOC_2VALS_GRID getMinlocMaxloc2ValsGrid0032
#define GET_MINLOC_MAXLOC_3VALS_GRID getMinlocMaxloc3ValsGrid0032
#endif//NTHREADS_COMPARE_INC ==   32

#   if  NTHREADS_COMPARE_INC ==   64
#define GET_MIN_BLCK getMinBlck0064
#define GET_MIN_GRID getMinGrid0064
#define GET_MAX_BLCK getMaxBlck0064
#define GET_MAX_GRID getMaxGrid0064
#define GET_MINLOC_BLCK getMinlocBlck0064
#define GET_MINLOC_GRID getMinlocGrid0064
#define GET_MAXLOC_BLCK getMaxlocBlck0064
#define GET_MAXLOC_GRID getMaxlocGrid0064
#define GET_MINLOC_MAXLOC_BLCK getMinlocMaxlocBlck0064
#define GET_MINLOC_MAXLOC_GRID getMinlocMaxlocGrid0064
#define GET_MINLOC_2VALS_GRID getMinloc2ValsGrid0064
#define GET_MINLOC_3VALS_GRID getMinloc3ValsGrid0064
#define GET_MAXLOC_2VALS_GRID getMaxloc2ValsGrid0064
#define GET_MAXLOC_3VALS_GRID getMaxloc3ValsGrid0064
#define GET_MINLOC_MAXLOC_2VALS_GRID getMinlocMaxloc2ValsGrid0064
#define GET_MINLOC_MAXLOC_3VALS_GRID getMinlocMaxloc3ValsGrid0064
#endif//NTHREADS_COMPARE_INC ==   64

#   if  NTHREADS_COMPARE_INC ==  128
#define GET_MIN_BLCK getMinBlck0128
#define GET_MIN_GRID getMinGrid0128
#define GET_MAX_BLCK getMaxBlck0128
#define GET_MAX_GRID getMaxGrid0128
#define GET_MINLOC_BLCK getMinlocBlck0128
#define GET_MINLOC_GRID getMinlocGrid0128
#define GET_MAXLOC_BLCK getMaxlocBlck0128
#define GET_MAXLOC_GRID getMaxlocGrid0128
#define GET_MINLOC_MAXLOC_BLCK getMinlocMaxlocBlck0128
#define GET_MINLOC_MAXLOC_GRID getMinlocMaxlocGrid0128
#define GET_MINLOC_2VALS_GRID getMinloc2ValsGrid0128
#define GET_MINLOC_3VALS_GRID getMinloc3ValsGrid0128
#define GET_MAXLOC_2VALS_GRID getMaxloc2ValsGrid0128
#define GET_MAXLOC_3VALS_GRID getMaxloc3ValsGrid0128
#define GET_MINLOC_MAXLOC_2VALS_GRID getMinlocMaxloc2ValsGrid0128
#define GET_MINLOC_MAXLOC_3VALS_GRID getMinlocMaxloc3ValsGrid0128
#endif//NTHREADS_COMPARE_INC ==  128

#   if  NTHREADS_COMPARE_INC ==  256
#define GET_MIN_BLCK getMinBlck0256
#define GET_MIN_GRID getMinGrid0256
#define GET_MAX_BLCK getMaxBlck0256
#define GET_MAX_GRID getMaxGrid0256
#define GET_MINLOC_BLCK getMinlocBlck0256
#define GET_MINLOC_GRID getMinlocGrid0256
#define GET_MAXLOC_BLCK getMaxlocBlck0256
#define GET_MAXLOC_GRID getMaxlocGrid0256
#define GET_MINLOC_MAXLOC_BLCK getMinlocMaxlocBlck0256
#define GET_MINLOC_MAXLOC_GRID getMinlocMaxlocGrid0256
#define GET_MINLOC_2VALS_GRID getMinloc2ValsGrid0256
#define GET_MINLOC_3VALS_GRID getMinloc3ValsGrid0256
#define GET_MAXLOC_2VALS_GRID getMaxloc2ValsGrid0256
#define GET_MAXLOC_3VALS_GRID getMaxloc3ValsGrid0256
#define GET_MINLOC_MAXLOC_2VALS_GRID getMinlocMaxloc2ValsGrid0256
#define GET_MINLOC_MAXLOC_3VALS_GRID getMinlocMaxloc3ValsGrid0256
#endif//NTHREADS_COMPARE_INC ==  256

#   if  NTHREADS_COMPARE_INC ==  512
#define GET_MIN_BLCK getMinBlck0512
#define GET_MIN_GRID getMinGrid0512
#define GET_MAX_BLCK getMaxBlck0512
#define GET_MAX_GRID getMaxGrid0512
#define GET_MINLOC_BLCK getMinlocBlck0512
#define GET_MINLOC_GRID getMinlocGrid0512
#define GET_MAXLOC_BLCK getMaxlocBlck0512
#define GET_MAXLOC_GRID getMaxlocGrid0512
#define GET_MINLOC_MAXLOC_BLCK getMinlocMaxlocBlck0512
#define GET_MINLOC_MAXLOC_GRID getMinlocMaxlocGrid0512
#define GET_MINLOC_2VALS_GRID getMinloc2ValsGrid0512
#define GET_MINLOC_3VALS_GRID getMinloc3ValsGrid0512
#define GET_MAXLOC_2VALS_GRID getMaxloc2ValsGrid0512
#define GET_MAXLOC_3VALS_GRID getMaxloc3ValsGrid0512
#define GET_MINLOC_MAXLOC_2VALS_GRID getMinlocMaxloc2ValsGrid0512
#define GET_MINLOC_MAXLOC_3VALS_GRID getMinlocMaxloc3ValsGrid0512
#endif//NTHREADS_COMPARE_INC ==  512

#   if  NTHREADS_COMPARE_INC == 1024
#define GET_MIN_BLCK getMinBlck1024
#define GET_MIN_GRID getMinGrid1024
#define GET_MAX_BLCK getMaxBlck1024
#define GET_MAX_GRID getMaxGrid1024
#define GET_MINLOC_BLCK getMinlocBlck1024
#define GET_MINLOC_GRID getMinlocGrid1024
#define GET_MAXLOC_BLCK getMaxlocBlck1024
#define GET_MAXLOC_GRID getMaxlocGrid1024
#define GET_MINLOC_MAXLOC_BLCK getMinlocMaxlocBlck1024
#define GET_MINLOC_MAXLOC_GRID getMinlocMaxlocGrid1024
#define GET_MINLOC_2VALS_GRID getMinloc2ValsGrid1024
#define GET_MINLOC_3VALS_GRID getMinloc3ValsGrid1024
#define GET_MAXLOC_2VALS_GRID getMaxloc2ValsGrid1024
#define GET_MAXLOC_3VALS_GRID getMaxloc3ValsGrid1024
#define GET_MINLOC_MAXLOC_2VALS_GRID getMinlocMaxloc2ValsGrid1024
#define GET_MINLOC_MAXLOC_3VALS_GRID getMinlocMaxloc3ValsGrid1024
#endif//NTHREADS_COMPARE_INC == 1024


#ifdef  COMPARE_DEL_CUH
#undef  COMPARE_DEL_CUH
#endif//COMPARE_DEL_CUH


#endif//COMPARE_INC_CUH
