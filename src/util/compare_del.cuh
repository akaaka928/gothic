/**
 * @file compare_del.cuh
 *
 * @brief Header file for comparing values on GPU
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
#ifndef COMPARE_DEL_CUH
#define COMPARE_DEL_CUH


#ifdef  COMPARE_INC_CUH

#undef  SHFL_MASK_COMPARE_INC

#undef  GET_MIN_BLCK
#undef  GET_MIN_GRID
#undef  GET_MAX_BLCK
#undef  GET_MAX_GRID
#undef  GET_MINLOC_BLCK
#undef  GET_MINLOC_GRID
#undef  GET_MAXLOC_BLCK
#undef  GET_MAXLOC_GRID
#undef  GET_MINLOC_MAXLOC_BLCK
#undef  GET_MINLOC_MAXLOC_GRID
#undef  GET_MINLOC_2VALS_GRID
#undef  GET_MINLOC_3VALS_GRID
#undef  GET_MAXLOC_2VALS_GRID
#undef  GET_MAXLOC_3VALS_GRID
#undef  GET_MINLOC_MAXLOC_2VALS_GRID
#undef  GET_MINLOC_MAXLOC_3VALS_GRID

#undef  COMPARE_INC_CUH
#endif//COMPARE_INC_CUH


#endif//COMPARE_DEL_CUH
