/**
 * @file external.h
 *
 * @brief Header file for setting an external fixed potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/31 (Wed)
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

#define N_EXT_POT_SPHE (1024)
/* #define N_EXT_POT_SPHE (4096) */

#   if  N_EXT_POT_SPHE > NRADBIN
#undef  N_EXT_POT_SPHE
#define N_EXT_POT_SPHE   NRADBIN
#endif//N_EXT_POT_SPHE > NRADBIN


#define N_EXT_POT_DISK (2048)


/* list of functions appeared in ``external.c'' */
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
void genExtPotTbl1D(const int kind, profile **prf, potential_field *pot);

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
void genSuperposedPotFld1D(const int kind, const int skind, profile **prf, potential_field * restrict sphe, potential_field * restrict disk);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
void superposePotFld1D(const int kind, const int skind, potential_field * restrict pot, potential_field * restrict sphe, potential_field * restrict disk);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
void extractDiskPotential(const int maxLev, const disk_data data, const potential_field sphe, disk_potential *disk);
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD


#endif//EXTERNAL_H
