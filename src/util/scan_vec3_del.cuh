/**
 * @file scan_vec3_del.cuh
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
#define SCAN_VEC3_DEL_CUH


#ifdef  SCAN_VEC3_INC_CUH

#undef  PREFIX_SUM_VEC3_BLCK
#undef  PREFIX_SUM_VEC3_GRID
#undef  PREFIX_SUM_VEC3_GRID_WITH_PARTITION
#undef   TOTAL_SUM_VEC3_BLCK
#undef   TOTAL_SUM_VEC3_GRID

#undef  SCAN_VEC3_INC_CUH
#endif//SCAN_VEC3_INC_CUH


#endif//SCAN_VEC3_DEL_CUH
