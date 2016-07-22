/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/14(Thu) 10:59:03
 *                                                                       *
 *    Header File for memory allocation code of N-body calculation       *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef ALLOCATE_DEV_H
#define ALLOCATE_DEV_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef STRUCTURE_H
#       include "../misc/structure.h"
#endif//STRUCTURE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "allocate_dev.cu"
//-------------------------------------------------------------------------
#ifdef __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  muse allocParticleDataSoA_dev
  (const int num
#ifdef  BLOCK_TIME_STEP
   , iparticle *body0, ulong **idx0, position **pos0, acceleration **acc0, velocity **vel0, ibody_time **ti0
#       ifdef  GENERATE_PHKEY_ON_DEVICE
   , iparticle *body1, ulong **idx1, position **pos1, acceleration **acc1, velocity **vel1, ibody_time **ti1
#       else///GENERATE_PHKEY_ON_DEVICE
   , position **pj, velocity **vj
#       endif//GENERATE_PHKEY_ON_DEVICE
#else///BLOCK_TIME_STEP
   , iparticle *body0, ulong **idx0, position **pos0, acceleration **acc0, real **vx0, real **vy0, real **vz0
#       ifdef  GENERATE_PHKEY_ON_DEVICE
   , iparticle *body1, ulong **idx1, position **pos1, acceleration **acc1, real **vx1, real **vy1, real **vz1
#       endif//GENERATE_PHKEY_ON_DEVICE
#endif//BLOCK_TIME_STEP
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
   , real **neighbor
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
   , position **encBall, position **encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
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
  //-----------------------------------------------------------------------
  void  freeParticleDataSoA_dev
  (ulong  *idx0, position  *pos0, acceleration  *acc0
#ifdef  BLOCK_TIME_STEP
   , velocity  *vel0, ibody_time  *ti0
#       ifdef  GENERATE_PHKEY_ON_DEVICE
   , ulong  *idx1, position  *pos1, acceleration  *acc1, velocity  *vel1, ibody_time  *ti1
#       else///GENERATE_PHKEY_ON_DEVICE
   , position  *pj, velocity  *vj
#       endif//GENERATE_PHKEY_ON_DEVICE
#else///BLOCK_TIME_STEP
   , real  *vx0, real  *vy0, real  *vz0
#       ifdef  GENERATE_PHKEY_ON_DEVICE
   , ulong  *idx1, position  *pos1, acceleration  *acc1, real  *vx1, real  *vy1, real  *vz1
#       endif//GENERATE_PHKEY_ON_DEVICE
#endif//BLOCK_TIME_STEP
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
   , real  *neighbor
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
   , position  *encBall, position  *encBall_hst
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
   );
  void  freeParticleDataSoA_hst
  (ulong  *idx_hst, position  *pos_hst, acceleration  *acc_hst,
#ifdef  BLOCK_TIME_STEP
   velocity  *vel_hst, ibody_time  *ti_hst
#else///BLOCK_TIME_STEP
   real  *vx_hst, real  *vy_hst, real  *vz_hst
#endif//BLOCK_TIME_STEP
   );
  //-----------------------------------------------------------------------
  muse allocParticleInfoSoA_dev
  (const int num, iparticle_treeinfo *info_dev,
   int **jtag_dev
#ifdef  COUNT_INTERACTIONS
   , int **Nj_dev, int **Nbuf_dev
#endif//COUNT_INTERACTIONS
   );
  muse allocParticleInfoSoA_hst
  (const int num, iparticle_treeinfo *info_hst
#ifndef MAKE_TREE_ON_DEVICE
   , int **jtag_hst
#endif//MAKE_TREE_ON_DEVICE
#ifdef  COUNT_INTERACTIONS
   , int **Nj_hst, int **Nbuf_hst
#endif//COUNT_INTERACTIONS
   );
  //-----------------------------------------------------------------------
  void  freeParticleInfoSoA_dev
  (int  *jtag_dev
#ifdef  COUNT_INTERACTIONS
   , int  *Nj_dev, int  *Nbuf_dev
#endif//COUNT_INTERACTIONS
   );
  void  freeParticleInfoSoA_hst
  (
#ifndef MAKE_TREE_ON_DEVICE
   int  *jtag_hst
#ifdef  COUNT_INTERACTIONS
   , int  *Nj_hst, int  *Nbuf_hst
#endif//COUNT_INTERACTIONS
#else///MAKE_TREE_ON_DEVICE
#ifdef  COUNT_INTERACTIONS
   int  *Nj_hst, int  *Nbuf_hst
#else///COUNT_INTERACTIONS
   void
#endif//COUNT_INTERACTIONS
#endif//MAKE_TREE_ON_DEVICE
   );
  //-----------------------------------------------------------------------
  /* muse allocAccel_dev(const int num, acceleration **acc_dev, acceleration **acc_hst); */
  /* void  freeAccel_dev(               acceleration  *acc_dev, acceleration  *acc_hst); */
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
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//ALLOCATE_DEV_H
//-------------------------------------------------------------------------
