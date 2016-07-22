/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/21(Thu) 11:59:36
 *                                                                       *
 *    Header File for Definition about Structures of N-body Simulation   *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef STRUCTURE_H
#define STRUCTURE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* INT_MAX is 2,147,483,647 */
/* float = 4 bytes, double = 8 bytes */
/*  int  = 4 bytes,  long  = 8 bytes */
//-------------------------------------------------------------------------
/* /\* 2^26 *\/ */
/* #define NUM_BODY_MAX (67108864) */
/* /\* 2^25 *\/ */
/* #define NUM_BODY_MAX (33554432) */
/* /\* 2^24 *\/ */
/* #define NUM_BODY_MAX (16777216) */
/* /\* 2^23 *\/ */
/* #define NUM_BODY_MAX (8388608) */
/* 2^22 */
#define NUM_BODY_MAX (4194304)
/* /\* 2^21 *\/ */
/* #define NUM_BODY_MAX (2097152) */
/* /\* 2^20 *\/ */
/* #define NUM_BODY_MAX (1048576) */
//-------------------------------------------------------------------------
#ifdef  RUN_ON_PC
/* 131072 (= 128K) is maximum value can ben executed on augustus */
#   if  NUM_BODY_MAX > 131072
#undef  NUM_BODY_MAX
#define NUM_BODY_MAX  (131072)
#endif//NUM_BODY_MAX > 131072
#endif//RUN_ON_PC
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* structure for memory allocation */
//-------------------------------------------------------------------------
typedef struct
{
  size_t host, device;
} muse;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* structures for N-body particles */
//-------------------------------------------------------------------------
#ifdef  MONITOR_ENERGY_ERROR
typedef struct
{
  double E0inv, errMax;
} energyError;
#endif//MONITOR_ENERGY_ERROR
//-------------------------------------------------------------------------
typedef struct
{
  ulong idx;
#ifdef  BLOCK_TIME_STEP
  /* real t0, t1, dt; */
  double t0, t1;
  real dt;
#endif//BLOCK_TIME_STEP
  real  x,  y,  z;
  real vx, vy, vz;
  real ax, ay, az;
  real m, pot;
} nbody_particle;
//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
typedef struct
{
  /* arrays to store particle information */
  real *pos;/* x0, y0, z0, x1, y1, z1, ... */
  real *vel;/* x0, y0, z0, x1, y1, z1, ... */
  real *acc;/* x0, y0, z0, x1, y1, z1, ... */
  real *m, *pot;
  ulong *idx;
} nbody_hdf5;
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
typedef struct __align__(16)
{
  real x, y, z, m;
} position;
//-------------------------------------------------------------------------
typedef struct __align__(16)
{
  real x, y, z, pot;
} acceleration;
//-------------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
typedef struct __align__(16)
{
  real x, y, z, dt;
} velocity;
#endif//BLOCK_TIME_STEP
//-------------------------------------------------------------------------
#ifdef  BLOCK_TIME_STEP
/* typedef struct __align__(8) */
/* { */
/*   real t0, t1; */
/* } ibody_time; */
typedef struct __align__(16)
{
  double t0, t1;
} ibody_time;
#endif//BLOCK_TIME_STEP
//-------------------------------------------------------------------------
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
#   if  defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
  real         *neighbor;
#endif//defined(BRUTE_FORCE_LOCALIZATION) && defined(LOCALIZE_I_PARTICLES)
  ulong        *idx;
#ifdef  RETURN_CENTER_BY_PHKEY_GENERATOR
  position *encBall, *encBall_hst;/* center and squared radius of enclosing ball which contains all i-particles */
#endif//RETURN_CENTER_BY_PHKEY_GENERATOR
#ifdef  GADGET_MAC
  real amin;
#endif//GADGET_MAC
} iparticle;
//-------------------------------------------------------------------------
typedef struct
{
  int *jtag;
#ifdef  COUNT_INTERACTIONS
  int *  Nj;
  int *Nbuf;
#endif//COUNT_INTERACTIONS
} iparticle_treeinfo;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//STRUCTURE_H
//-------------------------------------------------------------------------
