/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:32:21
 *                                                                       *
 *    Header File for memory allocation code of N-body calculation       *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef ALLOCATE_H
#define ALLOCATE_H
//-------------------------------------------------------------------------
#include "macro.h"
//-------------------------------------------------------------------------
#include "../misc/structure.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "allocate.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  /* muse allocParticleDataAoS(const int num, nbody_particle **body); */
  /* void  freeParticleDataAoS(               nbody_particle  *body); */
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  muse allocSnapshotArray(real **pos, real **vel, real **acc, real **m, real **pot, ulong **idx, const int num, nbody_hdf5 *data);
  void  freeSnapshotArray(real  *pos, real  *vel, real  *acc, real  *m, real  *pot, ulong  *idx);
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//ALLOCATE_H
//-------------------------------------------------------------------------
