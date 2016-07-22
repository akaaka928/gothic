/*************************************************************************\
 *                                                                       *
                  last updated on 2016/01/27(Wed) 16:53:42
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
//-- List of functions appeared in "allocate.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  muse allocParticleDataAoS(const int num, nbody_particle **body);
  void  freeParticleDataAoS(               nbody_particle  *body);
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
