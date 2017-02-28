/**
 * @file eddington.h
 *
 * @brief Header file for Eddington's formula
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/02/21 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef EDDINGTON_H
#define EDDINGTON_H


#include "macro.h"

#include "../init/profile.h"


#define NDIVIDE_GAUSSQD (16)
#define NENEBIN (1048576)


/**
 * @struct dist_func
 *
 * @brief structure for distribution function
 */
typedef struct
{
  real ene, val;
} dist_func;


#define NKIND_MAX (8)


/* list of functions appeared in ``eddington.c'' */
void integrateEddingtonFormula(const int skind, profile **prf, dist_func **fene);


#endif//EDDINGTON_H
