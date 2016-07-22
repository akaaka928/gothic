/*************************************************************************\
 *                                                                       *
                  last updated on 2016/01/27(Wed) 16:53:26
 *                                                                       *
 *    Header File for converting type of arrays (SoA <--> AoS)           *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef CONVERT_H
#define CONVERT_H
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
//-- List of functions appeared in "convert.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  /* data packing and unpacking */
  void copyStr2Ary(const int head, const int Ni, nbody_particle *src, iparticle dst);
  void copyAry2Str(const int head, const int Ni, iparticle src, nbody_particle *dst);
  //-----------------------------------------------------------------------
  /* additional function for Leap-Frog integrator */
  void backVel(const int head, const int tail, nbody_particle *body, real dt);
  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  void copyAry4Snapshot(const int head, const int Ni, iparticle src, nbody_hdf5 dst);
  void backVel4Snapshot(const int head, const int tail, nbody_hdf5 body, real dt);
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//CONVERT_H
//-------------------------------------------------------------------------
