/**
 * @file allocate_dev.cu
 *
 * @brief Source code for memory allocation on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/12/17 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <helper_cuda.h>

#include "macro.h"
#include "cudalib.h"

#include "structure.h"
#include "allocate_dev.h"

#include "../tree/walk_dev.h"/**< to read NTHREADS */


/**
 * @fn allocParticleDataSoA_hst
 *
 * @brief Allocate memory for N-body particles as SoA on CPU.
 */
extern "C"
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
 )
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  /** the size of the array is set to be a multiple of NTHREADS */
  size_t size = (size_t)num;
  if( (num % NTHREADS) != 0 )
    size += (size_t)(NTHREADS - (num % NTHREADS));

  /** memory allocation and simple confirmation */
  mycudaMallocHost((void **)idx_hst, size * sizeof(       ulong));  alloc.host += size * sizeof(       ulong);
  mycudaMallocHost((void **)pos_hst, size * sizeof(    position));  alloc.host += size * sizeof(    position);
  mycudaMallocHost((void **)acc_hst, size * sizeof(acceleration));  alloc.host += size * sizeof(acceleration);
#ifdef  BLOCK_TIME_STEP
  mycudaMallocHost((void **)vel_hst, size * sizeof(    velocity));  alloc.host += size * sizeof(    velocity);
  mycudaMallocHost((void **) ti_hst, size * sizeof(  ibody_time));  alloc.host += size * sizeof(  ibody_time);
#else///BLOCK_TIME_STEP
  mycudaMallocHost((void **) vx_hst, size * sizeof(        real));  alloc.host += size * sizeof(        real);
  mycudaMallocHost((void **) vy_hst, size * sizeof(        real));  alloc.host += size * sizeof(        real);
  mycudaMallocHost((void **) vz_hst, size * sizeof(        real));  alloc.host += size * sizeof(        real);
#endif//BLOCK_TIME_STEP

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  mycudaMallocHost((void **)acc_ext, size * sizeof(acceleration));  alloc.host += size * sizeof(acceleration);
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  /** commit arrays to the utility structure */
  body_hst->pos  = *pos_hst;
  body_hst->acc  = *acc_hst;
#ifdef  BLOCK_TIME_STEP
  body_hst->vel  = *vel_hst;
  body_hst->time = * ti_hst;
#else///BLOCK_TIME_STEP
  body_hst->vx   = * vx_hst;
  body_hst->vy   = * vy_hst;
  body_hst->vz   = * vz_hst;
#endif//BLOCK_TIME_STEP
  body_hst->idx  = *idx_hst;

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  body_hst->acc_ext = *acc_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freeParticleDataSoA_hst
 *
 * @brief Deallocate memory for N-body particles as SoA on CPU.
 */
extern "C"
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
 )
{
  __NOTE__("%s\n", "start");

  mycudaFreeHost(idx_hst);
  mycudaFreeHost(pos_hst);
  mycudaFreeHost(acc_hst);
#ifdef  BLOCK_TIME_STEP
  mycudaFreeHost(vel_hst);
  mycudaFreeHost( ti_hst);
#else///BLOCK_TIME_STEP
  mycudaFreeHost( vx_hst);
  mycudaFreeHost( vy_hst);
  mycudaFreeHost( vz_hst);
#endif//BLOCK_TIME_STEP

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  mycudaFreeHost(acc_ext);
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
}


/**
 * @fn allocParticleInfoSoA_dev
 *
 * @brief Allocate memory for N-body particles as SoA on GPU.
 */
extern "C"
muse allocParticleInfoSoA_dev
(const int num, iparticle_treeinfo *info_dev,
 int **jtag_dev
#ifdef  COUNT_INTERACTIONS
 , int **Nj_dev, int **Nbuf_dev
#endif//COUNT_INTERACTIONS
 )
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  /** the size of the array is set to be a multiple of NTHREADS */
  size_t size = (size_t)num;
  if( (num % NTHREADS) != 0 )
    size += (size_t)(NTHREADS - (num % NTHREADS));

  /** memory allocation and simple confirmation */
  mycudaMalloc((void **)jtag_dev, size * sizeof(int));  alloc.device += size * sizeof(int);
#ifdef  COUNT_INTERACTIONS
  mycudaMalloc((void **)  Nj_dev, size * sizeof(int));  alloc.device += size * sizeof(int);
  mycudaMalloc((void **)Nbuf_dev, size * sizeof(int));  alloc.device += size * sizeof(int);
#endif//COUNT_INTERACTIONS

  /** commit arrays to the utility structure */
  info_dev->jtag = *jtag_dev;
#ifdef  COUNT_INTERACTIONS
  info_dev->  Nj = *  Nj_dev;
  info_dev->Nbuf = *Nbuf_dev;
#endif//COUNT_INTERACTIONS


  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn allocParticleInfoSoA_hst
 *
 * @brief Allocate memory for N-body particles as SoA on CPU.
 */
extern "C"
muse allocParticleInfoSoA_hst
(const int num
#ifdef  COUNT_INTERACTIONS
 , iparticle_treeinfo *info_hst, int **Nj_hst, int **Nbuf_hst
#endif//COUNT_INTERACTIONS
 )
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  /** the size of the array is set to be a multiple of NTHREADS */
  size_t size = (size_t)num;
  if( (num % NTHREADS) != 0 )
    size += (size_t)(NTHREADS - (num % NTHREADS));

#ifdef  COUNT_INTERACTIONS
  /** memory allocation and simple confirmation */
  mycudaMallocHost((void **)  Nj_hst, size * sizeof(int));  alloc.host += size * sizeof(int);
  mycudaMallocHost((void **)Nbuf_hst, size * sizeof(int));  alloc.host += size * sizeof(int);

  /** commit arrays to the utility structure */
  info_hst->  Nj = *  Nj_hst;
  info_hst->Nbuf = *Nbuf_hst;
#endif//COUNT_INTERACTIONS


  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freeParticleInfoSoA_dev
 *
 * @brief Deallocate memory for N-body particles as SoA on GPU.
 */
extern "C"
void  freeParticleInfoSoA_dev
(int  *jtag_dev
#ifdef  COUNT_INTERACTIONS
 , int  *Nj_dev, int  *Nbuf_dev
#endif//COUNT_INTERACTIONS
 )
{
  __NOTE__("%s\n", "start");

  mycudaFree(jtag_dev);
#ifdef  COUNT_INTERACTIONS
  mycudaFree(  Nj_dev);
  mycudaFree(Nbuf_dev);
#endif//COUNT_INTERACTIONS

  __NOTE__("%s\n", "end");
}


#ifdef  COUNT_INTERACTIONS
/**
 * @fn freeParticleInfoSoA_hst
 *
 * @brief Deallocate memory for N-body particles as SoA on CPU.
 */
extern "C"
void  freeParticleInfoSoA_hst(int  *Nj_hst, int  *Nbuf_hst)
{
  __NOTE__("%s\n", "start");

  mycudaFreeHost(  Nj_hst);
  mycudaFreeHost(Nbuf_hst);

  __NOTE__("%s\n", "end");
}
#endif//COUNT_INTERACTIONS


/**
 * @fn allocAccel_dev
 *
 * @brief Allocate memory for particle acceleration on GPU.
 */
extern "C"
muse allocAccel_dev
(const int num, acceleration **acc_dev, acceleration **acc_hst
#ifdef  GADGET_MAC
 , acceleration **old_dev
#endif//GADGET_MAC
)
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  /** the size of the array is set to be a multiple of NTHREADS */
  size_t size = (size_t)num;
  if( (num % NTHREADS) != 0 )
    size += (size_t)(NTHREADS - (num % NTHREADS));

  /** memory allocation and simple confirmation */
  mycudaMalloc    ((void **)acc_dev, size * sizeof(acceleration));  alloc.device += size * sizeof(acceleration);
  mycudaMallocHost((void **)acc_hst, size * sizeof(acceleration));  alloc.host   += size * sizeof(acceleration);
#ifdef  GADGET_MAC
  mycudaMalloc    ((void **)old_dev, size * sizeof(acceleration));  alloc.device += size * sizeof(acceleration);
#endif//GADGET_MAC


  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freeAccel_dev
 *
 * @brief Deallocate memory for particle acceleration on GPU.
 */
extern "C"
void  freeAccel_dev
(acceleration  *acc_dev, acceleration  *acc_hst
#ifdef  GADGET_MAC
 , acceleration  *old_dev
#endif//GADGET_MAC
)
{
  __NOTE__("%s\n", "start");

  mycudaFree    (acc_dev);
  mycudaFreeHost(acc_hst);
#ifdef  GADGET_MAC
  mycudaFree    (old_dev);
#endif//GADGET_MAC

  __NOTE__("%s\n", "end");
}


#ifdef  SET_SINK_PARTICLES
extern "C"
muse allocSinkParticleSoA_dev
(position **pold, velocity **vold, position **pnew, velocity **vnew, velocity **mom, ulong **tag, int **list, real **lmax, ulong **tag_hst,
   const int Nsink, sinkparticle *sink)
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  size_t size = (size_t)Nsink;

  /** memory allocation and simple confirmation */
  mycudaMalloc((void **)pold, size * sizeof(position));  alloc.device += size * sizeof(position);
  mycudaMalloc((void **)pnew, size * sizeof(position));  alloc.device += size * sizeof(position);
  mycudaMalloc((void **)vold, size * sizeof(velocity));  alloc.device += size * sizeof(velocity);
  mycudaMalloc((void **)vnew, size * sizeof(velocity));  alloc.device += size * sizeof(velocity);
  mycudaMalloc((void **) mom, size * sizeof(velocity));  alloc.device += size * sizeof(velocity);
  mycudaMalloc((void **) tag, size * sizeof(   ulong));  alloc.device += size * sizeof(   ulong);
  mycudaMalloc((void **)list, size * sizeof(     int));  alloc.device += size * sizeof(     int);
  mycudaMalloc((void **)lmax, size * sizeof(    real));  alloc.device += size * sizeof(    real);
  mycudaMallocHost((void **)tag_hst, size * sizeof(ulong));  alloc.host += size * sizeof(ulong);

  /** commit arrays to the utility structure */
sink->pold = *pold;
sink->pnew = *pnew;
sink->vold = *vold;
sink->vnew = *vnew;
sink-> mom = * mom;
sink-> tag = * tag;
sink->list = *list;
sink->lmax2 = *lmax;


  __NOTE__("%s\n", "end");
  return (alloc);
}
extern "C"
void  freeSinkParticleSoA_dev
(position  *pold, velocity  *vold, position  *pnew, velocity  *vnew, velocity  *mom, ulong  *tag, int  *list, real  *lmax, ulong  *tag_hst)
{
  __NOTE__("%s\n", "start");

  mycudaFree(pold);
  mycudaFree(pnew);
  mycudaFree(vold);
  mycudaFree(vnew);
  mycudaFree( mom);
  mycudaFree( tag);
  mycudaFree(list);
  mycudaFree(lmax);
  mycudaFreeHost(tag_hst);

  __NOTE__("%s\n", "end");
}
#endif//SET_SINK_PARTICLES
