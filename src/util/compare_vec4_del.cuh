/**
 * @file compare_vec4_del.cuh
 *
 * @brief Header file for comparing values in 4-components vector form on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/29 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef COMPARE_VEC4_DEL_CUH
#define COMPARE_VEC4_DEL_CUH


#ifdef  COMPARE_VEC4_INC_CUH

#undef  GET_VEC4_MIN_BLCK
#undef  GET_VEC4_MIN_GRID
#undef  GET_VEC4_MAX_BLCK
#undef  GET_VEC4_MAX_GRID
#undef  GET_VEC4_MIN_MAX_BLCK
#undef  GET_VEC4_MIN_MAX_GRID

#undef  COMPARE_VEC4_INC_CUH
#endif//COMPARE_VEC4_INC_CUH


#endif//COMPARE_VEC4_DEL_CUH
