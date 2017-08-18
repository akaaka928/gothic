/**
 * @file scan_del.cuh
 *
 * @brief Header file for parallel prefix sum library on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/08/16 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef SCAN_DEL_CUH
#define SCAN_DEL_CUH


#ifdef  SCAN_INC_CUH

#undef  PREFIX_SUM_BLCK
#undef  PREFIX_SUM_GRID
#undef  PREFIX_SUM_GRID_WITH_PARTITION
#undef   TOTAL_SUM_BLCK

#undef  SCAN_INC_CUH
#endif//SCAN_INC_CUH


#endif//SCAN_DEL_CUH
