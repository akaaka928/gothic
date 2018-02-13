/**
 * @file potential_dev.h
 *
 * @brief Header file for calculating gravity from external fixed-potential field
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
#ifndef POTENTIAL_DEV_H
#define POTENTIAL_DEV_H
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD


#include "macro.h"
#include "cudalib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#endif//USE_HDF5_FORMAT

#include "../file/io.h"


/* list of functions appeared in ``potential_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__

  void  freeSphericalPotentialTable_dev
  (pot2  *Phi
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
   , real  *rad
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
   );

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
  void  freeDiskPotentialTable_dev
  (real  *Phi,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
   real  *RR, real  *zz
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
   disk_grav *FRz
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
   );
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK

  void calcExternalGravity_dev
  (const iparticle pi, const potential_field sphe
#ifdef  BLOCK_TIME_STEP
   , const int thrd, const int grpNum, laneinfo * RESTRICT laneInfo
#else///BLOCK_TIME_STEP
   , const int Ni
#endif//BLOCK_TIME_STEP
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
   , const disk_potential disk
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
   );

  /* read fixed potential field */
  muse  readFixedPotentialTableSpherical
  (const int unit, char file[], potential_field *pot_tbl, pot2 **Phi
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
   , real **rad
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#ifdef  USE_HDF5_FORMAT
   , hdf5struct type
#endif//USE_HDF5_FORMAT
   );

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
  muse  readFixedPotentialTableDisk
  (const int unit, char file[], real **Phi_dev, pot2 **Phi_sphe_dev,
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
   real **RR_dev, real **zz_dev, real **rad_sphe_dev,
#else///ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
   disk_grav **FRz_dev,
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
   disk_potential *disk
#ifdef  USE_HDF5_FORMAT
   , hdf5struct type
#endif//USE_HDF5_FORMAT
   );
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK

#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//SET_EXTERNAL_POTENTIAL_FIELD
#endif//POTENTIAL_DEV_H
