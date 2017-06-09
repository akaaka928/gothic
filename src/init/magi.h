/**
 * @file magi.h
 *
 * @brief Header file for MAGI (MAny-component Galaxy Initializer)
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/16 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef MAGI_H
#define MAGI_H


#include <stdbool.h>

#include "macro.h"


/**
 * @def PROGRESS_REPORT_ON
 *
 * @brief activate progress report
 */
#ifndef PROGRESS_REPORT_ON
#define PROGRESS_REPORT_ON
#endif//PROGRESS_REPORT_ON

#define NMAX_GAUSS_QD (51)
#define NTBL_GAUSS_QD ((NMAX_GAUSS_QD >> 1) + (NMAX_GAUSS_QD & 1))
#define NINTBIN NMAX_GAUSS_QD

#define N_PRINT_LINES_ASCII (8192)


#endif//MAGI_H
