/**
 * @file convert.h
 *
 * @brief Header file for converting type of arrays (SoA <--> AoS)
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
#ifndef CONVERT_H
#define CONVERT_H


#include "macro.h"

#include "../misc/structure.h"


/* list of functions appeared in ``convert.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  /** additional function for Leap-Frog integrator */
  void backVel(const int head, const int tail, iparticle body, real dt);

#ifdef  USE_HDF5_FORMAT
  void copyAry4Snapshot(const int head, const int Ni, iparticle src, nbody_hdf5 dst);
  void backVel4Snapshot(const int head, const int tail, nbody_hdf5 body, real dt);
#endif//USE_HDF5_FORMAT
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//CONVERT_H
