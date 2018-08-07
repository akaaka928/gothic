/**
 * @file m31coord.h
 *
 * @brief Header file for coordinate transformation related to M31
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2018/07/23 (Mon)
 *
 * Copyright (C) 2018 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef M31COORD_H
#define M31COORD_H


#include "macro.h"
#include "../misc/structure.h"


static const real zm31 = CAST_D2R(776.0);/**< same with Komiyama et al. (2018); median of NED */
/* static const real zm31 = CAST_D2R(773.0);/\**< M31 distance in Conn et al. (2016) *\/ */

static const real vm31x = CAST_D2R(0.0);
static const real vm31y = CAST_D2R(0.0);
static const real vm31z = CAST_D2R(-300.0);



/* list of functions appeared in ``m31coord.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  void setRotationMatrix(real rot[restrict][3], real inv[restrict][3]);
  void standard_coordinate(const int num, nbody_aos *body, real rot[restrict][3], real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos);
#ifdef  __CUDACC__
}
#endif//__CUDACC__



#endif//M31COORD_H
