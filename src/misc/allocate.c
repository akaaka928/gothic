/*************************************************************************\
 *                                                                       *
                  last updated on 2016/01/27(Wed) 16:54:02
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
muse allocParticleDataAoS(const int num, nbody_particle **body)
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
  /* memory allocation and simple confirmation */
  alloc.host +=                    size * sizeof(nbody_particle);
  *body = (nbody_particle *)malloc(size * sizeof(nbody_particle));
  if( *body == NULL ){    __KILL__(stderr, "ERROR: failure to allocate body");  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  freeParticleDataAoS(nbody_particle  *body)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  free(body);
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
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
