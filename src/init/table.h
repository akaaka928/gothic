/**
 * @file table.h
 *
 * @brief Header file for density profiles in the machine-readable table format
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
#ifndef TABLE_H
#define TABLE_H


#include "../init/profile.h"


#define NFIT (4)
#define NPUT (256)


/* list of functions appeared in ``table.c'' */
void getInterpolatedDensityProfile(const int num, profile * restrict prf, double * restrict xx, double * restrict ff);
void setDensityProfileTable(profile *prf, const double rs, char *file);


#endif//TABLE_H
