/**
 * @file external.h
 *
 * @brief Header file for setting an external fixed potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/15 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef EXTERNAL_H
#define EXTERNAL_H


#include "../misc/structure.h"
#include "../init/profile.h"

#define N_EXT_POT_SPHE (16384)

#   if  N_EXT_POT_SPHE > NRADBIN
#undef  N_EXT_POT_SPHE
#define N_EXT_POT_SPHE   NRADBIN
#endif//N_EXT_POT_SPHE > NRADBIN


#define N_EXT_POT_DISK (2048)


/* list of functions appeared in ``external.c'' */
void genExtPotTbl1D(const int kind, profile **prf, potential_field *pot);
void superposePotFld1D(const int kind, const int skind, potential_field * restrict pot, potential_field * restrict sphe, potential_field * restrict disk);


#endif//EXTERNAL_H
