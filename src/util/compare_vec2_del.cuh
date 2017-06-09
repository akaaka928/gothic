/**
 * @file compare_vec2_del.cuh
 *
 * @brief Header file for comparing values in 2-components vector form on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/04/05 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef COMPARE_VEC2_DEL_CUH
#define COMPARE_VEC2_DEL_CUH


#ifdef  COMPARE_VEC2_INC_CUH

#undef  GET_VEC2_MIN_BLCK
#undef  GET_VEC2_MIN_GRID
#undef  GET_VEC2_MAX_BLCK
#undef  GET_VEC2_MAX_GRID
#undef  GET_VEC2_MIN_MAX_BLCK
#undef  GET_VEC2_MIN_MAX_GRID

#undef  COMPARE_VEC2_INC_CUH
#endif//COMPARE_VEC2_INC_CUH


#endif//COMPARE_VEC2_DEL_CUH
