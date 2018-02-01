/**
 * @file allocate.c
 *
 * @brief Source code for memory allocation in GOTHIC and MAGI
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/31 (Wed)
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
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 acceleration **ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
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
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  *ext = (acceleration *)malloc(size * sizeof(acceleration));  if( *ext == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ext\n");  }
  alloc.host +=                 size * sizeof(acceleration);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
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
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  body->acc_ext = *ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
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
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 acceleration  *ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
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
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  free(ext);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
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
muse allocSnapshotArray
(real **pos, real **vel, real **acc, real **m, real **pot, ulong **idx,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 real **acc_ext, real **pot_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
 const int num, nbody_hdf5 *data)
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

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  *acc_ext = (real  *)malloc(3 * size * sizeof(real) );  alloc.host += 3 * size * sizeof(real);
  *pot_ext = (real  *)malloc(    size * sizeof(real) );  alloc.host +=     size * sizeof(real);
  if( *acc_ext == NULL ){    __KILL__(stderr, "ERROR: failure to allocate acc_ext\n"  );  }
  if( *pot_ext == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pot_ext\n");  }
  data->acc_ext = *acc_ext;
  data->pot_ext = *pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD

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
void  freeSnapshotArray
(real  *pos, real  *vel, real  *acc, real  *m, real  *pot, ulong  *idx
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
 , real  *acc_ext, real  *pot_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
)
{
  __NOTE__("%s\n", "start");

  free(pos);
  free(vel);
  free(acc);
  free(m);
  free(pot);
  free(idx);
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  free(acc_ext);
  free(pot_ext);
#endif//SET_EXTERNAL_POTENTIAL_FIELD

  __NOTE__("%s\n", "end");
}
#endif//USE_HDF5_FORMAT


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
/**
 * @fn allocPotentialField
 *
 * @brief Allocate arrays for external fixed potential field.
 *
 * @return (rad) radius
 * @return (Phi) potential and its 2nd-derivative
 * @return (dat) potential field of each component
 * @param (num) number of data points for potential field
 * @param (kind) number of components
 * @return (sphe) superposed potential field of spherical components
 * @param (kind) number of spherical symmetric components
 * @return (disk) superposed (spherical averaged) potential field of disk components
 */
muse allocPotentialField(real **rad, pot2 **Phi, potential_field **dat, const int num, const int kind, potential_field *sphe, const int skind, potential_field *disk)
{
  __NOTE__("%s\n", "start");
  muse alloc = {0, 0};

  size_t size = (size_t)num * (size_t)(kind + 1 + (kind > skind));

  /* allocate data array for external potential field */
  *rad = (real *)malloc(size * sizeof(real));  alloc.host += size * sizeof(real);
  *Phi = (pot2 *)malloc(size * sizeof(pot2));  alloc.host += size * sizeof(pot2);
  if( *rad == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rad\n"  );  }
  if( *Phi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate Phi\n"  );  }

  /* allocate structure for storing external potential field */
  *dat = (potential_field *)malloc(kind * sizeof(potential_field));  alloc.host += kind * sizeof(potential_field);
  if( *dat == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dat\n"  );  }

  /* asign arrays */
  for(int ii = 0; ii < kind; ii++){
    (*dat)[ii].rad = &((*rad)[ii * num]);
    (*dat)[ii].Phi = &((*Phi)[ii * num]);
  }/* for(int ii = 0; ii < kind; ii++){ */

  sphe->rad = &((*rad)[kind * num]);
  sphe->Phi = &((*Phi)[kind * num]);

  if( kind > skind ){
    disk->rad = &((*rad)[(kind + 1) * num]);
    disk->Phi = &((*Phi)[(kind + 1) * num]);
  }/* if( kind > skind ){ */

  __NOTE__("%s\n", "end");
  return (alloc);
}
/**
 * @fn freePotentialField
 *
 * @brief Deallocate arrays for external fixed potential field.
 *
 * @param (rad) radius
 * @param (Phi) potential and its 2nd-derivative
 * @param (dat) potential field of each component
 */
void  freePotentialField(real  *rad, pot2  *Phi, potential_field  *dat)
{
  __NOTE__("%s\n", "start");

  free(rad);
  free(Phi);
  free(dat);

  __NOTE__("%s\n", "end");
}

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK
muse allocDiskPotential(real **RR, real **zz, real **Phi, const int maxLev, const int NR, const int Nz, disk_potential *disk)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};
  size_t size;

  /* allocate data array for external potential field */
  size = maxLev *  NR                ;  *RR  = (real *)malloc(size * sizeof(real));  alloc.host += size * sizeof(real);
  size = maxLev            *  Nz     ;  *zz  = (real *)malloc(size * sizeof(real));  alloc.host += size * sizeof(real);
  size = maxLev * (NR + 1) * (Nz + 1);  *Phi = (real *)malloc(size * sizeof(real));  alloc.host += size * sizeof(real);
  if( *RR  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate RR\n"  );  }
  if( *zz  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate zz\n"  );  }
  if( *Phi == NULL ){    __KILL__(stderr, "ERROR: failure to allocate Phi\n"  );  }

  disk->RR  = *RR;
  disk->zz  = *zz;
  disk->Phi = *Phi;

  __NOTE__("%s\n", "end");
  return (alloc);
}

void  freeDiskPotential(real  *RR, real  *zz, real  *Phi)
{
  __NOTE__("%s\n", "start");

  free(RR);
  free(zz);
  free(Phi);

  __NOTE__("%s\n", "end");
}
#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
#endif//SET_EXTERNAL_POTENTIAL_FIELD
