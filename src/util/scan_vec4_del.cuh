/**
 * @file scan_vec4_del.cuh
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
#define SCAN_VEC4_DEL_CUH


#ifdef  SCAN_VEC4_INC_CUH

#undef  PREFIX_SUM_VEC4_BLCK
#undef  PREFIX_SUM_VEC4_GRID
#undef  PREFIX_SUM_VEC4_GRID_WITH_PARTITION
#undef   TOTAL_SUM_VEC4_BLCK
#undef   TOTAL_SUM_VEC4_GRID

#undef  SCAN_VEC4_INC_CUH
#endif//SCAN_VEC4_INC_CUH


#endif//SCAN_VEC4_DEL_CUH
