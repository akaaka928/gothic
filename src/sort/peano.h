/**
 * @file peano.h
 *
 * @brief Header file for calculating Peano--Hilbert space-filling curve
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/01 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef PEANO_H
#define PEANO_H


#include <sys/time.h>

#include "macro.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"


/**
 * @def MAXIMUM_PHKEY_LEVEL
 *
 * @brief Maximum level for PH-key hierarchy
 */
/* #define MAXIMUM_PHKEY_LEVEL (10) */
#define MAXIMUM_PHKEY_LEVEL (21)

#define NUM_PHKEY_LEVEL (1 + MAXIMUM_PHKEY_LEVEL)
/* 1 means the level for the root cell */


/**
 * @typedef PHint
 *
 * @brief data type for PH-key
 */
#if     MAXIMUM_PHKEY_LEVEL <= 10
typedef uint  PHint;
#else
typedef ulong PHint;
#endif//MAXIMUM_PHKEY_LEVEL <= 10


/**
 * @struct PHinfo
 *
 * @brief structure for PH-key hierarchy
 */
typedef struct __align__(16)
{
  int level, num, head, nmax;  /**< nmax is the maximum number of particles contained a tree cell in this level */
} PHinfo;


#endif//PEANO_H
