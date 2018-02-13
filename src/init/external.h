/**
 * @file external.h
 *
 * @brief Header file for setting an external fixed potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/02/13 (Tue)
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

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
/* #define N_EXT_POT_SPHE (1024) */
#define N_EXT_POT_SPHE (2048)
/* #define N_EXT_POT_SPHE (4096) */
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#define N_EXT_POT_SPHE (16384)
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

#   if  N_EXT_POT_SPHE > NRADBIN
#undef  N_EXT_POT_SPHE
#define N_EXT_POT_SPHE   NRADBIN
#endif//N_EXT_POT_SPHE > NRADBIN

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
/* Nz = 2048 ==>> NR * Nz * 3 * 4B = 192MiB */
/* Nz = 4096 ==>> NR * Nz * 3 * 4B = 768MiB */
/* #define NZ_EXT_POT_DISK (2048) */
#define NZ_EXT_POT_DISK (4096)
#define NR_EXT_POT_DISK (NZ_EXT_POT_DISK << 2)
/* 1/64 = 1.5625e-2; h = zd / 64 */
#define EXT_POT_DISK_MIN_LENGTH (1.5625e-2)
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK

/* list of functions appeared in ``external.c'' */
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
void genExtPotTbl1D(const int kind, profile **prf, potential_field *pot
#ifndef ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
		    , const double invlogrbin
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
		    );

#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
void genSuperposedPotFld1D(const int kind, const int skind, profile **prf, potential_field * restrict sphe, potential_field * restrict disk);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
void superposePotFld1D(const int kind, const int skind, potential_field * restrict pot, potential_field * restrict sphe, potential_field * restrict disk);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
void extractDiskPotential(const int maxLev, const disk_data data, const potential_field sphe, disk_potential *disk);
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
void extractDiskPotential(const int maxLev, const int ndisk, disk_data *data, const potential_field sphe, disk_potential *disk);
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD


#endif//EXTERNAL_H
