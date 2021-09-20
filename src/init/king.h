/**
 * @file king.h
 *
 * @brief Header file for King model (lowered isothermal model)
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
#ifndef KING_H
#define KING_H


#include "../init/profile.h"
#include "../init/magi.h"


#ifdef  PROGRESS_REPORT_ON
#define KING_PROGRESS_REPORT_ON (100000)
#endif//PROGRESS_REPORT_ON

#define NADD_KING (1048576)


/* list of functions appeared in ``king.c'' */
void setDensityProfileKing(profile *prf, profile_cfg *cfg);


#endif//KING_H
