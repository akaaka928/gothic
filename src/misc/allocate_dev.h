/**
 * @file allocate_dev.h
 *
 * @brief Header file for memory allocation on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/03/12 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef ALLOCATE_DEV_H
#define ALLOCATE_DEV_H


#include "macro.h"

#include "../misc/structure.h"


/* list of functions appeared in ``allocate_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  muse allocParticleDataSoA_hst
  (const int num, iparticle *body_hst,
   ulong **idx_hst, position **pos_hst, acceleration **acc_hst,
#ifdef  BLOCK_TIME_STEP
   velocity **vel_hst, ibody_time **ti_hst
#else///BLOCK_TIME_STEP
   real **vx_hst, real **vy_hst, real **vz_hst
#endif//BLOCK_TIME_STEP
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
   , acceleration **acc_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
   );

  void  freeParticleDataSoA_hst
  (ulong  *idx_hst, position  *pos_hst, acceleration  *acc_hst,
#ifdef  BLOCK_TIME_STEP
   velocity  *vel_hst, ibody_time  *ti_hst
#else///BLOCK_TIME_STEP
   real  *vx_hst, real  *vy_hst, real  *vz_hst
#endif//BLOCK_TIME_STEP
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
   , acceleration  *acc_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
   );

  muse allocParticleInfoSoA_dev
  (const int num, iparticle_treeinfo *info_dev,
   int **jtag_dev
#ifdef  COUNT_INTERACTIONS
   , int **Nj_dev, int **Nbuf_dev
#endif//COUNT_INTERACTIONS
   );
  muse allocParticleInfoSoA_hst
  (const int num
#ifdef  COUNT_INTERACTIONS
   , iparticle_treeinfo *info_hst, int **Nj_hst, int **Nbuf_hst
#endif//COUNT_INTERACTIONS
   );

  void  freeParticleInfoSoA_dev
  (int  *jtag_dev
#ifdef  COUNT_INTERACTIONS
   , int  *Nj_dev, int  *Nbuf_dev
#endif//COUNT_INTERACTIONS
   );
#ifdef  COUNT_INTERACTIONS
  void  freeParticleInfoSoA_hst(int  *Nj_hst, int  *Nbuf_hst);
#endif//COUNT_INTERACTIONS

  muse allocAccel_dev(const int num, acceleration **acc_dev, acceleration **acc_hst
#ifdef  GADGET_MAC
		      , acceleration **old_dev
#endif//GADGET_MAC
		      );
  void  freeAccel_dev(acceleration  *acc_dev, acceleration  *acc_hst
#ifdef  GADGET_MAC
		      , acceleration  *old_dev
#endif//GADGET_MAC
		      );
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//ALLOCATE_DEV_H
