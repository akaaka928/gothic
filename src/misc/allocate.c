/*************************************************************************\
 *                                                                       *
                  last updated on 2016/08/09(Tue) 17:22:29
 *                                                                       *
 *    Memory Allocation Code of N-body calculation                       *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
//-------------------------------------------------------------------------
#include <macro.h>
//-------------------------------------------------------------------------
#include "structure.h"
#include "allocate.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* the size of the array is set to be a multiple of NTHREADS */
  size_t size = (size_t)num;
  if( (num % NSIMD) != 0 )
    size += (size_t)(NSIMD - (num % NSIMD));
  //-----------------------------------------------------------------------
  /* memory allocation and simple confirmation */
  *idx = (       ulong *)malloc(size * sizeof(       ulong));  if( *idx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate idx");  }
  alloc.host +=                 size * sizeof(       ulong);
  *pos = (    position *)malloc(size * sizeof(    position));  if( *pos == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pos");  }
  alloc.host +=                 size * sizeof(    position);
  *acc = (acceleration *)malloc(size * sizeof(acceleration));  if( *acc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate acc");  }
  alloc.host +=                 size * sizeof(acceleration);
#ifdef  BLOCK_TIME_STEP
  *vel = (  velocity *)malloc(size * sizeof(  velocity));  if( *vel == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vel");  }
  alloc.host +=               size * sizeof(  velocity);
  * ti = (ibody_time *)malloc(size * sizeof(ibody_time));  if( * ti == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ti");  }
  alloc.host +=               size * sizeof(ibody_time);
#else///BLOCK_TIME_STEP
  vx = (real *)malloc(size * sizeof(real));  if( *vx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vx");  }
  alloc.host +=       size * sizeof(real);
  vy = (real *)malloc(size * sizeof(real));  if( *vy == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vy");  }
  alloc.host +=       size * sizeof(real);
  vz = (real *)malloc(size * sizeof(real));  if( *vz == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vz");  }
  alloc.host +=       size * sizeof(real);
#endif//BLOCK_TIME_STEP
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  freeParticleData
(ulong  *idx, position  *pos, acceleration  *acc,
#ifdef  BLOCK_TIME_STEP
 velocity  *vel, ibody_time  *ti
#else///BLOCK_TIME_STEP
 real  *vx, real  *vy, real  *vz
#endif//BLOCK_TIME_STEP
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* muse allocParticleDataAoS(const int num, nbody_particle **body) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "start"); */
/*   //----------------------------------------------------------------------- */
/*   muse alloc = {0, 0}; */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   /\* the size of the array is set to be a multiple of NSIMD *\/ */
/*   size_t size = (size_t)num; */
/*   if( (num % NSIMD) != 0 ) */
/*     size += (size_t)(NSIMD - (num % NSIMD)); */
/*   //----------------------------------------------------------------------- */
/*   /\* memory allocation and simple confirmation *\/ */
/*   alloc.host +=                    size * sizeof(nbody_particle); */
/*   *body = (nbody_particle *)malloc(size * sizeof(nbody_particle)); */
/*   if( *body == NULL ){    __KILL__(stderr, "ERROR: failure to allocate body");  } */
/*   //----------------------------------------------------------------------- */
/*   // */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "end"); */
/*   //----------------------------------------------------------------------- */
/*   return (alloc); */
/*   //----------------------------------------------------------------------- */
/* } */
/* //------------------------------------------------------------------------- */
/* void  freeParticleDataAoS(nbody_particle  *body) */
/* { */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "start"); */
/*   //----------------------------------------------------------------------- */
/*   free(body); */
/*   //----------------------------------------------------------------------- */
/*   __NOTE__("%s\n", "end"); */
/*   //----------------------------------------------------------------------- */
/* } */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
//-------------------------------------------------------------------------
muse allocSnapshotArray(real **pos, real **vel, real **acc, real **m, real **pot, ulong **idx, const int num, nbody_hdf5 *data)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* the size of the array is set to be a multiple of NSIMD */
  size_t size = (size_t)num;
  if( (num % NSIMD) != 0 )
    size += (size_t)(NSIMD - (num % NSIMD));
  //-----------------------------------------------------------------------
  /* allocate particle array for output snapshot */
  *pos = (real  *)malloc(3 * size * sizeof(real) );  alloc.host += 3 * size * sizeof(real);
  *vel = (real  *)malloc(3 * size * sizeof(real) );  alloc.host += 3 * size * sizeof(real);
  *acc = (real  *)malloc(3 * size * sizeof(real) );  alloc.host += 3 * size * sizeof(real);
  * m  = (real  *)malloc(    size * sizeof(real) );  alloc.host +=     size * sizeof(real);
  *pot = (real  *)malloc(    size * sizeof(real) );  alloc.host +=     size * sizeof(real);
  *idx = (ulong *)malloc(    size * sizeof(ulong));  alloc.host +=     size * sizeof(ulong);
  if( *pos == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pos"  );  }
  if( *vel == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vel"  );  }
  if( *acc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate acc"  );  }
  if( * m  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate m"  );  }
  if( *pot == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pot");  }
  if( *idx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate idx");  }
  /* asign arrays */
  data->pos = *pos;
  data->vel = *vel;
  data->acc = *acc;
  data->m = *m;
  data->pot = *pot;
  data->idx = *idx;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  freeSnapshotArray(real  *pos, real  *vel, real  *acc, real  *m, real  *pot, ulong  *idx)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  free(pos);
  free(vel);
  free(acc);
  free(m);
  free(pot);
  free(idx);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
