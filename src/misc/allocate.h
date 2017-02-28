/**
 * @file allocate.h
 *
 * @brief Header file for memory allocation in GOTHIC and MAGI
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
#ifndef ALLOCATE_H
#define ALLOCATE_H


#include "macro.h"

#include "../misc/structure.h"


/* list of functions appeared in ``allocate.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  muse allocParticleData
  (const int num, iparticle *body,
   ulong **idx, position **pos, acceleration **acc,
#ifdef  BLOCK_TIME_STEP
   velocity **vel, ibody_time **ti
#else///BLOCK_TIME_STEP
   real **vx, real **vy, real **vz
#endif//BLOCK_TIME_STEP
   );
  void  freeParticleData
  (ulong  *idx, position  *pos, acceleration  *acc,
#ifdef  BLOCK_TIME_STEP
   velocity  *vel, ibody_time  *ti
#else///BLOCK_TIME_STEP
   real  *vx, real  *vy, real  *vz
#endif//BLOCK_TIME_STEP
   );

#ifdef  USE_HDF5_FORMAT
  muse allocSnapshotArray(real **pos, real **vel, real **acc, real **m, real **pot, ulong **idx, const int num, nbody_hdf5 *data);
  void  freeSnapshotArray(real  *pos, real  *vel, real  *acc, real  *m, real  *pot, ulong  *idx);
#endif//USE_HDF5_FORMAT
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//ALLOCATE_H
