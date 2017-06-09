/**
 * @file allocate_dev.cu
 *
 * @brief Source code for memory allocation on GPU
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
 * @fn allocParticleDataSoA_dev
 *
 * @brief Allocate memory for N-body particles as SoA on GPU.
 */
extern "C"
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
 )
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  /** the size of the array is set to be a multiple of NTHREADS */
  size_t size = (size_t)num;
  if( (num % NTHREADS) != 0 )
    size += (size_t)(NTHREADS - (num % NTHREADS));

  /** memory allocation and simple confirmation */
  mycudaMalloc((void **)idx0, size * sizeof(ulong));  mycudaMalloc((void **)pos0, size * sizeof(position));  mycudaMalloc((void **)acc0, size * sizeof(acceleration));
  alloc.device +=             size * sizeof(ulong) ;  alloc.device +=             size * sizeof(position) ;  alloc.device +=             size * sizeof(acceleration) ;
  mycudaMalloc((void **)idx1, size * sizeof(ulong));  mycudaMalloc((void **)pos1, size * sizeof(position));  mycudaMalloc((void **)acc1, size * sizeof(acceleration));
  alloc.device +=             size * sizeof(ulong) ;  alloc.device +=             size * sizeof(position) ;  alloc.device +=             size * sizeof(acceleration) ;
#ifdef  BLOCK_TIME_STEP
  mycudaMalloc((void **)vel0, size * sizeof(velocity));  mycudaMalloc((void **)ti0, size * sizeof(ibody_time));
  alloc.device +=             size * sizeof(velocity) ;  alloc.device +=            size * sizeof(ibody_time) ;
  mycudaMalloc((void **)vel1, size * sizeof(velocity));  mycudaMalloc((void **)ti1, size * sizeof(ibody_time));
  alloc.device +=             size * sizeof(velocity) ;  alloc.device +=            size * sizeof(ibody_time) ;
#else///BLOCK_TIME_STEP
  mycudaMalloc((void **)vx0, size * sizeof(real));  mycudaMalloc((void **)vy0, size * sizeof(real));  mycudaMalloc((void **)vz0, size * sizeof(real));
  alloc.device +=            size * sizeof(real) ;  alloc.device +=            size * sizeof(real) ;  alloc.device +=            size * sizeof(real) ;
  mycudaMalloc((void **)vx1, size * sizeof(real));  mycudaMalloc((void **)vy1, size * sizeof(real));  mycudaMalloc((void **)vz1, size * sizeof(real));
  alloc.device +=            size * sizeof(real) ;  alloc.device +=            size * sizeof(real) ;  alloc.device +=            size * sizeof(real) ;
#endif//BLOCK_TIME_STEP

  /** commit arrays to the utility structure */
  body0->idx = *idx0;  body0->pos = *pos0;  body0->acc = *acc0;
  body1->idx = *idx1;  body1->pos = *pos1;  body1->acc = *acc1;
#ifdef  BLOCK_TIME_STEP
  body0->vel = *vel0;  body0->time = *ti0;
  body1->vel = *vel1;  body1->time = *ti1;
  body1->jpos = *pos0;  body1->jvel = *vel0;
  body0->jpos = *pos1;  body0->jvel = *vel1;
#else///BLOCK_TIME_STEP
  body0->vx = *vx0;  body0->vy = *vy0;  body0->vz = *vz0;
  body1->vx = *vx1;  body1->vy = *vy1;  body1->vz = *vz1;
#endif//BLOCK_TIME_STEP

#ifdef  GADGET_MAC
  body0->acc_old = *acc1;
  body1->acc_old = *acc0;
#endif//GADGET_MAC

  mycudaMalloc((void **)neighbor, size * sizeof(real));
  alloc.device +=                 size * sizeof(real) ;
  body0->neighbor = *neighbor;
  body1->neighbor = *neighbor;

#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  mycudaMalloc    ((void **)encBall    , sizeof(position));  alloc.device +=                sizeof(position);
  mycudaMallocHost((void **)encBall_hst, sizeof(position));  alloc.host   +=                sizeof(position);
  body0->encBall = *encBall;  body0->encBall_hst = *encBall_hst;
  body1->encBall = *encBall;  body1->encBall_hst = *encBall_hst;
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR

#ifdef  DPADD_FOR_ACC
  mycudaMalloc((void **)tmp, size * sizeof(DPacc));
  alloc.device +=            size * sizeof(DPacc) ;
  body0->tmp = *tmp;
  body1->tmp = *tmp;
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
  mycudaMalloc((void **)res, size * sizeof(acceleration));
  alloc.device +=            size * sizeof(acceleration) ;
  body0->res = *res;
  body1->res = *res;
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))


  __NOTE__("%s\n", "end");
  return (alloc);
}


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


  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn freeParticleDataSoA_dev
 *
 * @brief Deallocate memory for N-body particles as SoA on GPU.
 */
extern "C"
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
 )
{
  __NOTE__("%s\n", "start");


  mycudaFree(idx0);  mycudaFree(pos0);  mycudaFree(acc0);
  mycudaFree(idx1);  mycudaFree(pos1);  mycudaFree(acc1);
#ifdef  BLOCK_TIME_STEP
  mycudaFree(vel0);  mycudaFree(ti0);
  mycudaFree(vel1);  mycudaFree(ti1);
#else///BLOCK_TIME_STEP
  mycudaFree(vx0);  mycudaFree(vy0);  mycudaFree(vz0);
  mycudaFree(vx1);  mycudaFree(vy1);  mycudaFree(vz1);
#endif//BLOCK_TIME_STEP
  mycudaFree(neighbor);

#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  mycudaFree    (encBall);
  mycudaFreeHost(encBall_hst);
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  DPADD_FOR_ACC
  mycudaFree(tmp);
#endif//DPADD_FOR_ACC
#   if  defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))
  mycudaFree(res);
#endif//defined(KAHAN_SUM_CORRECTION) && defined(ACCURATE_ACCUMULATION) && (!defined(SERIALIZED_EXECUTION) || (NWARP > 1))


  __NOTE__("%s\n", "end");
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
(const int num, iparticle_treeinfo *info_hst
#ifdef  COUNT_INTERACTIONS
 , int **Nj_hst, int **Nbuf_hst
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
#ifdef  COUNT_INTERACTIONS
  mycudaMallocHost((void **)  Nj_hst, size * sizeof(int));  alloc.host += size * sizeof(int);
  mycudaMallocHost((void **)Nbuf_hst, size * sizeof(int));  alloc.host += size * sizeof(int);
#endif//COUNT_INTERACTIONS

  /** commit arrays to the utility structure */
#ifdef  COUNT_INTERACTIONS
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
