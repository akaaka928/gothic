/**
 * @file allocate_dev.h
 *
 * @brief Header file for memory allocation on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/02/28 (Tue)
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
#include "../tree/walk_dev.h"


/* list of functions appeared in ``allocate_dev.cu'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  muse allocParticleDataSoA_dev
  (const int num
#ifdef  BLOCK_TIME_STEP
   , iparticle *body0, ulong **idx0, position **pos0, acceleration **acc0, velocity **vel0, ibody_time **ti0
   , iparticle *body1, ulong **idx1, position **pos1, acceleration **acc1, velocity **vel1, ibody_time **ti1
#else///BLOCK_TIME_STEP
   , iparticle *body0, ulong **idx0, position **pos0, acceleration **acc0, real **vx0, real **vy0, real **vz0
   , iparticle *body1, ulong **idx1, position **pos1, acceleration **acc1, real **vx1, real **vy1, real **vz1
#endif//BLOCK_TIME_STEP
   , real **neighbor
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
   , position **encBall, position **encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  DPADD_FOR_ACC
   , DPacc **tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
   , acceleration **res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
   );
  muse allocParticleDataSoA_hst
  (const int num, iparticle *body_hst,
   ulong **idx_hst, position **pos_hst, acceleration **acc_hst,
#ifdef  BLOCK_TIME_STEP
   velocity **vel_hst, ibody_time **ti_hst
#else///BLOCK_TIME_STEP
   real **vx_hst, real **vy_hst, real **vz_hst
#endif//BLOCK_TIME_STEP
   );

  void  freeParticleDataSoA_dev
  (ulong  *idx0, position  *pos0, acceleration  *acc0
#ifdef  BLOCK_TIME_STEP
   , velocity  *vel0, ibody_time  *ti0
   , ulong  *idx1, position  *pos1, acceleration  *acc1, velocity  *vel1, ibody_time  *ti1
#else///BLOCK_TIME_STEP
   , real  *vx0, real  *vy0, real  *vz0
   , ulong  *idx1, position  *pos1, acceleration  *acc1, real  *vx1, real  *vy1, real  *vz1
#endif//BLOCK_TIME_STEP
   , real  *neighbor
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
   , position  *encBall, position  *encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  DPADD_FOR_ACC
   , DPacc  *tmp
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
   , acceleration  *res
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
   );
  void  freeParticleDataSoA_hst
  (ulong  *idx_hst, position  *pos_hst, acceleration  *acc_hst,
#ifdef  BLOCK_TIME_STEP
   velocity  *vel_hst, ibody_time  *ti_hst
#else///BLOCK_TIME_STEP
   real  *vx_hst, real  *vy_hst, real  *vz_hst
#endif//BLOCK_TIME_STEP
   );

  muse allocParticleInfoSoA_dev
  (const int num, iparticle_treeinfo *info_dev,
   int **jtag_dev
#ifdef  COUNT_INTERACTIONS
   , int **Nj_dev, int **Nbuf_dev
#endif//COUNT_INTERACTIONS
   );
  muse allocParticleInfoSoA_hst
  (const int num, iparticle_treeinfo *info_hst
#ifdef  COUNT_INTERACTIONS
   , int **Nj_hst, int **Nbuf_hst
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
