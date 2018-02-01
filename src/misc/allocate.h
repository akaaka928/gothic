/**
 * @file allocate.h
 *
 * @brief Header file for memory allocation in GOTHIC and MAGI
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
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
   acceleration **ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
   velocity **vel, ibody_time **ti
#else///BLOCK_TIME_STEP
   real **vx, real **vy, real **vz
#endif//BLOCK_TIME_STEP
   );
  void  freeParticleData
  (ulong  *idx, position  *pos, acceleration  *acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
   acceleration  *ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
   velocity  *vel, ibody_time  *ti
#else///BLOCK_TIME_STEP
   real  *vx, real  *vy, real  *vz
#endif//BLOCK_TIME_STEP
   );

#ifdef  USE_HDF5_FORMAT
  muse allocSnapshotArray
  (real **pos, real **vel, real **acc, real **m, real **pot, ulong **idx,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
   real **acc_ext, real **pot_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
   const int num, nbody_hdf5 *data);
  void  freeSnapshotArray
  (real  *pos, real  *vel, real  *acc, real  *m, real  *pot, ulong  *idx
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
   , real  *acc_ext, real  *pot_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
   );
#endif//USE_HDF5_FORMAT

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  muse allocPotentialField(real **rad, pot2 **Phi, potential_field **dat, const int num, const int kind, potential_field *sphe, const int skind, potential_field *disk);
  void  freePotentialField(real  *rad, pot2  *Phi, potential_field  *dat);

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
  muse allocDiskPotential(real **RR, real **zz, real **Phi, const int maxLev, const int NR, const int Nz, disk_potential *disk);
  void  freeDiskPotential(real  *RR, real  *zz, real  *Phi);
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD

#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//ALLOCATE_H
