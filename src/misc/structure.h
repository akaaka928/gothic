/**
 * @file structure.h
 *
 * @brief Header file for definition about structures in GOTHIC and MAGI
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/04/03 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef STRUCTURE_H
#define STRUCTURE_H


#include "macro.h"


/* INT_MAX is 2,147,483,647 = 2^31 - 1 */
/* float = 4 bytes, double = 8 bytes */
/*  int  = 4 bytes,  long  = 8 bytes */

/**
 * @def NUM_BODY_MAX
 *
 * @brief Maximum number of N-body particles per MPI process.
 */
/* /\* 2^26 *\/ */
/* #define NUM_BODY_MAX (67108864) */
/* /\* 2^25 *\/ */
/* #define NUM_BODY_MAX (33554432) */
/* 2^24 */
#define NUM_BODY_MAX (16777216)
/* /\* 2^23 *\/ */
/* #define NUM_BODY_MAX (8388608) */
/* /\* 2^22 *\/ */
/* #define NUM_BODY_MAX (4194304) */
/* /\* 2^21 *\/ */
/* #define NUM_BODY_MAX (2097152) */
/* /\* 2^20 *\/ */
/* #define NUM_BODY_MAX (1048576) */


#ifdef  RUN_ON_PC
/* 131072 (= 128K) is maximum value can ben executed on augustus (GTX 750 Ti) */
#   if  NUM_BODY_MAX > 131072
#undef  NUM_BODY_MAX
#define NUM_BODY_MAX  (131072)
#endif//NUM_BODY_MAX > 131072
#endif//RUN_ON_PC


/**
 * @struct muse
 *
 * @brief structure for memory allocation
 */
typedef struct
{
  size_t host, device;
} muse;


#ifdef  MONITOR_ENERGY_ERROR
/**
 * @struct energyError
 *
 * @brief structure for monitoring energy error
 */
typedef struct
{
  double E0inv, errMax;
} energyError;
#endif//MONITOR_ENERGY_ERROR


#ifdef  USE_HDF5_FORMAT
/**
 * @struct nbody_hdf5
 *
 * @brief structure for file I/O using HDF5
 */
typedef struct
{
  real *pos;/**< x0, y0, z0, w0, x1, y1, z1, ... */
  real *vel;/**< x0, y0, z0, w0, x1, y1, z1, ... */
  real *acc;/**< x0, y0, z0, w0, x1, y1, z1, ... */
  real *m, *pot;
  ulong *idx;
} nbody_hdf5;
#endif//USE_HDF5_FORMAT


/**
 * @struct position
 *
 * @brief structure for particle position (16 byte aligned)
 */
typedef struct __align__(16)
{
  real x, y, z, m;
} position;
/**
 * @struct acceleration
 *
 * @brief structure for particle acceleration (16 byte aligned)
 */
typedef struct __align__(16)
{
  real x, y, z, pot;
} acceleration;
#ifdef  DPADD_FOR_ACC
/**
 * @struct DPacc
 *
 * @brief structure for particle acceleration (32 byte aligned)
 */
typedef struct __align__(32)
{
  double x, y, z, pot;
} DPacc;
#endif//DPADD_FOR_ACC

#ifdef  BLOCK_TIME_STEP
/**
 * @struct velocity
 *
 * @brief structure for particle velocity (16 byte aligned)
 */
typedef struct __align__(16)
{
  real x, y, z, dt;
} velocity;
/**
 * @struct ibody_time
 *
 * @brief structure for particle time (16 byte aligned)
 */
typedef struct __align__(16)
{
  double t0, t1;
} ibody_time;
#endif//BLOCK_TIME_STEP


/**
 * @struct iparticle
 *
 * @brief structure for i-particles (SoA)
 */
typedef struct
{
  position     *pos;
  acceleration *acc;
#ifdef  GADGET_MAC
  acceleration *acc_old;
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
  velocity     *vel;
  position     *jpos;
  velocity     *jvel;
  ibody_time   *time;
#else///BLOCK_TIME_STEP
  real         *vx;
  real         *vy;
  real         *vz;
#endif//BLOCK_TIME_STEP
  real         *neighbor;
  ulong        *idx;
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  position *encBall, *encBall_hst;/**< center and squared radius of enclosing ball which contains all i-particles */
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  DPADD_FOR_ACC
  DPacc *tmp;
#endif//DPADD_FOR_ACC
#ifdef  KAHAN_SUM_CORRECTION
  acceleration *res;/**< residual for Kahan summation */
#endif//KAHAN_SUM_CORRECTION
#ifdef  GADGET_MAC
  real amin;
#endif//GADGET_MAC
} iparticle;


/**
 * @struct iparticle_treeinfo
 *
 * @brief structure for information on tree structure (SoA)
 */
typedef struct
{
  int *jtag;
#ifdef  COUNT_INTERACTIONS
  int *  Nj;
  int *Nbuf;
#endif//COUNT_INTERACTIONS
} iparticle_treeinfo;


#endif//STRUCTURE_H
