/**
 * @file potential_dev.h
 *
 * @brief Header file for calculating gravity from external fixed-potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/17 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef POTENTIAL_DEV_H
#define POTENTIAL_DEV_H


#include "macro.h"
#include "cudalib.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "../file/io.h"
#endif//USE_HDF5_FORMAT


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
#define READ_SUPERPOSED_TABLE (-1)
#endif//SET_EXTERNAL_POTENTIAL_FIELD


/* list of functions appeared in ``potential_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__

  muse allocSphericalPotentialTable_dev(real **rad, pot2 **Phi, const int Nr);
  void  freeSphericalPotentialTable_dev(real  *rad, pot2  *Phi);

  muse allocSphericalPotentialTable_hst(real **rad, pot2 **Phi, const int Nr);
  void  freeSphericalPotentialTable_hst(real  *rad, pot2  *Phi);

  void setSphericalPotentialTable_dev(real *rad_hst, pot2 *Phi_hst, real *rad_dev, pot2 *Phi_dev, const int Nr);

  void calcExternalGravity_dev
  (const iparticle pi, const potential_field sphe
#ifdef  BLOCK_TIME_STEP
   , const int thrd, const int grpNum, laneinfo * RESTRICT laneInfo
#else///BLOCK_TIME_STEP
   , const int Ni
#endif//BLOCK_TIME_STEP
   );

  /* read fixed potential field */
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_SPHERICAL
  muse  readFixedPotentialTableSpherical(const int unit, char cfg[], potential_field *pot_tbl, real **rad, pot2 **Phi
#ifdef  USE_HDF5_FORMAT
					 , hdf5struct type
#endif//USE_HDF5_FORMAT
					 );
#endif//SET_EXTERNAL_POTENTIAL_FIELD_SPHERICAL
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//POTENTIAL_DEV_H
