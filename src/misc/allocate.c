/**
 * @file allocate.c
 *
 * @brief Source code for memory allocation in GOTHIC and MAGI
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/02/21 (Tue)
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

#include "macro.h"

#include "structure.h"
#include "allocate.h"


/**
 * @fn allocParticleData
 *
 * @brief Allocate arrays for N-body particles.
 *
 * @param (num) number of N-body particles
 * @return (body) structure contains N-body particle data (SoA)
 * @return (idx) index of N-body particles
 * @return (pos) position of N-body particles
 * @return (acc) acceleration of N-body particles
 * @return (vel) velocity of N-body particles (when BLOCK_TIME_STEP is enabled)
 * @return (ti) time of N-body particles (when BLOCK_TIME_STEP is enabled)
 * @return (vx) x-component of velocity of N-body particles (when BLOCK_TIME_STEP is disabled)
 * @return (vy) y-component of velocity of N-body particles (when BLOCK_TIME_STEP is disabled)
 * @return (vz) z-component of velocity of N-body particles (when BLOCK_TIME_STEP is disabled)
 */
muse allocParticleData
(const int num, iparticle *body,
 ulong **idx, position **pos, acceleration **acc,
#ifdef  BLOCK_TIME_STEP
 velocity **vel, ibody_time **ti
#else///BLOCK_TIME_STEP
 real **vx, real **vy, real **vz
#endif//BLOCK_TIME_STEP
 )
{
  __NOTE__("%s\n", "start");
  muse alloc = {0, 0};

  /* the size of the array is set to be a multiple of NTHREADS */
  size_t size = (size_t)num;
  if( (num % NSIMD) != 0 )
    size += (size_t)(NSIMD - (num % NSIMD));

  /* memory allocation and simple confirmation */
  *idx = (       ulong *)malloc(size * sizeof(       ulong));  if( *idx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate idx\n");  }
  alloc.host +=                 size * sizeof(       ulong);
  *pos = (    position *)malloc(size * sizeof(    position));  if( *pos == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pos\n");  }
  alloc.host +=                 size * sizeof(    position);
  *acc = (acceleration *)malloc(size * sizeof(acceleration));  if( *acc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate acc\n");  }
  alloc.host +=                 size * sizeof(acceleration);
#ifdef  BLOCK_TIME_STEP
  *vel = (  velocity *)malloc(size * sizeof(  velocity));  if( *vel == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vel\n");  }
  alloc.host +=               size * sizeof(  velocity);
  * ti = (ibody_time *)malloc(size * sizeof(ibody_time));  if( * ti == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ti\n");  }
  alloc.host +=               size * sizeof(ibody_time);
#else///BLOCK_TIME_STEP
  *vx = (real *)malloc(size * sizeof(real));  if( *vx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vx\n");  }
  alloc.host +=       size * sizeof(real);
  *vy = (real *)malloc(size * sizeof(real));  if( *vy == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vy\n");  }
  alloc.host +=       size * sizeof(real);
  *vz = (real *)malloc(size * sizeof(real));  if( *vz == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vz\n");  }
  alloc.host +=       size * sizeof(real);
#endif//BLOCK_TIME_STEP

  /* commit arrays to the utility structure */
  body->pos  = *pos;
  body->acc  = *acc;
#ifdef  BLOCK_TIME_STEP
  body->vel  = *vel;
  body->time = * ti;
#else///BLOCK_TIME_STEP
  body->vx   = * vx;
  body->vy   = * vy;
  body->vz   = * vz;
#endif//BLOCK_TIME_STEP
  body->idx  = *idx;

  __NOTE__("%s\n", "end");
  return (alloc);
}
/**
 * @fn freeParticleData
 *
 * @brief Deallocate arrays for N-body particles.
 *
 * @param (idx) index of N-body particles
 * @param (pos) position of N-body particles
 * @param (acc) acceleration of N-body particles
 * @param (vel) velocity of N-body particles (when BLOCK_TIME_STEP is enabled)
 * @param (ti) time of N-body particles (when BLOCK_TIME_STEP is enabled)
 * @param (vx) x-component of velocity of N-body particles (when BLOCK_TIME_STEP is disabled)
 * @param (vy) y-component of velocity of N-body particles (when BLOCK_TIME_STEP is disabled)
 * @param (vz) z-component of velocity of N-body particles (when BLOCK_TIME_STEP is disabled)
 */
void  freeParticleData
(ulong  *idx, position  *pos, acceleration  *acc,
#ifdef  BLOCK_TIME_STEP
 velocity  *vel, ibody_time  *ti
#else///BLOCK_TIME_STEP
 real  *vx, real  *vy, real  *vz
#endif//BLOCK_TIME_STEP
 )
{
  __NOTE__("%s\n", "start");

  free(idx);
  free(pos);
  free(acc);
#ifdef  BLOCK_TIME_STEP
  free(vel);
  free( ti);
#else///BLOCK_TIME_STEP
  free( vx);
  free( vy);
  free( vz);
#endif//BLOCK_TIME_STEP

  __NOTE__("%s\n", "end");
}


#ifdef  USE_HDF5_FORMAT
/**
 * @fn allocSnapshotArray
 *
 * @brief Allocate arrays for snapshots in HDF5 format.
 *
 * @return (pos) position of N-body particles
 * @return (vel) velocity of N-body particles
 * @return (acc) acceleration of N-body particles
 * @return (m) mass of N-body particles
 * @return (pot) potential of N-body particles
 * @return (idx) index of N-body particles
 * @param (num) number of N-body particles
 * @return (data) structure contains N-body particle data (SoA)
 */
muse allocSnapshotArray(real **pos, real **vel, real **acc, real **m, real **pot, ulong **idx, const int num, nbody_hdf5 *data)
{
  __NOTE__("%s\n", "start");
  muse alloc = {0, 0};

  /* the size of the array is set to be a multiple of NSIMD */
  size_t size = (size_t)num;
  if( (num % NSIMD) != 0 )
    size += (size_t)(NSIMD - (num % NSIMD));

  /* allocate particle array for output snapshot */
  *pos = (real  *)malloc(3 * size * sizeof(real) );  alloc.host += 3 * size * sizeof(real);
  *vel = (real  *)malloc(3 * size * sizeof(real) );  alloc.host += 3 * size * sizeof(real);
  *acc = (real  *)malloc(3 * size * sizeof(real) );  alloc.host += 3 * size * sizeof(real);
  * m  = (real  *)malloc(    size * sizeof(real) );  alloc.host +=     size * sizeof(real);
  *pot = (real  *)malloc(    size * sizeof(real) );  alloc.host +=     size * sizeof(real);
  *idx = (ulong *)malloc(    size * sizeof(ulong));  alloc.host +=     size * sizeof(ulong);
  if( *pos == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pos\n"  );  }
  if( *vel == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vel\n"  );  }
  if( *acc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate acc\n"  );  }
  if( * m  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate m\n"  );  }
  if( *pot == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pot\n");  }
  if( *idx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate idx\n");  }
  /* asign arrays */
  data->pos = *pos;
  data->vel = *vel;
  data->acc = *acc;
  data->m = *m;
  data->pot = *pot;
  data->idx = *idx;

  __NOTE__("%s\n", "end");
  return (alloc);
}
/**
 * @fn freeSnapshotArray
 *
 * @brief Deallocate arrays for snapshots in HDF5 format.
 *
 * @param (pos) position of N-body particles
 * @param (vel) velocity of N-body particles
 * @param (acc) acceleration of N-body particles
 * @param (m) mass of N-body particles
 * @param (pot) potential of N-body particles
 * @param (idx) index of N-body particles
 */
void  freeSnapshotArray(real  *pos, real  *vel, real  *acc, real  *m, real  *pot, ulong  *idx)
{
  __NOTE__("%s\n", "start");

  free(pos);
  free(vel);
  free(acc);
  free(m);
  free(pot);
  free(idx);

  __NOTE__("%s\n", "end");
}
#endif//USE_HDF5_FORMAT
