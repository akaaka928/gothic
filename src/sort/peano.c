/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:42:00
 *                                                                       *
 *    Key generation of Peano-Hilbert space filling curve                *
 *    sort N-body particles to obey Peano-Hilbert space filling curve    *
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
#include <sys/time.h>
#include <math.h>
//-------------------------------------------------------------------------
#include "macro.h"
//-------------------------------------------------------------------------
#include "../misc/benchmark.h"
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "peano.h"
//-------------------------------------------------------------------------
#include "../tree/make.h"/* <-- required to read NLEAF */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef GENERATE_PHKEY_ON_DEVICE
//-------------------------------------------------------------------------
muse allocPeanoHilbertKey
(const int num, PHint **peano, identity **tag
#ifndef CALC_MULTIPOLE_ON_DEVICE
 , PHinfo **info
#endif//CALC_MULTIPOLE_ON_DEVICE
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* Peano--Hilbert key of N-body particles */
  //-----------------------------------------------------------------------
  /* the size of the array is set to be a multiple of NSIMD */
  size_t size = (size_t)num;
  if( (num % NSIMD) != 0 )
    size += (size_t)(NSIMD - (num % NSIMD));
  //-----------------------------------------------------------------------
  /* memory allocation and simple confirmation */
  *peano  = (   PHint *)malloc(size * sizeof(   PHint));  if( * peano == NULL ){    __KILL__(stderr, "ERROR: failure to allocate peano\n");  }
  alloc.host +=                size * sizeof(   PHint);
  *tag    = (identity *)malloc(size * sizeof(identity));  if( *   tag == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tag\n");  }
  alloc.host +=                size * sizeof(identity);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifndef CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
  /* Properties of Peano--Hilbert key of tree cells in the hierarchical structure */
  //-----------------------------------------------------------------------
  /* memory allocation and simple confirmation */
  *info = (PHinfo *)malloc(NUM_PHKEY_LEVEL * sizeof(PHinfo));  if( *info == NULL ){    __KILL__(stderr, "ERROR: failure to allocate info\n");  }
  alloc.host +=            NUM_PHKEY_LEVEL * sizeof(PHinfo);
  //-----------------------------------------------------------------------
  initPHinfo(*info);
  //-----------------------------------------------------------------------
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  freePeanoHilbertKey
(PHint  *peano, identity  *tag
#ifndef CALC_MULTIPOLE_ON_DEVICE
 , PHinfo  *info
#endif//CALC_MULTIPOLE_ON_DEVICE
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  free(peano);
  free(tag);
#ifndef CALC_MULTIPOLE_ON_DEVICE
  free(info);
#endif//CALC_MULTIPOLE_ON_DEVICE
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void calcBoxSize(const int num, iparticle body, real min[], real *dnew)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  real max[3];
#pragma unroll
  for(int ii = 0; ii < 3; ii++){
    min[ii] =  REAL_MAX;
    max[ii] = -REAL_MAX;
  }
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++){
    if( body.pos[ii].x < min[0] )      min[0] = body.pos[ii].x;
    if( body.pos[ii].y < min[1] )      min[1] = body.pos[ii].y;
    if( body.pos[ii].z < min[2] )      min[2] = body.pos[ii].z;
    if( body.pos[ii].x > max[0] )      max[0] = body.pos[ii].x;
    if( body.pos[ii].y > max[1] )      max[1] = body.pos[ii].y;
    if( body.pos[ii].z > max[2] )      max[2] = body.pos[ii].z;
  }
  //-----------------------------------------------------------------------
#ifndef NDEBUG
  printf("min: %f, %f, %f\n", min[0], min[1], min[2]);
  printf("max: %f, %f, %f\n", max[0], max[1], max[2]);
#if 0
  exit(1);
#endif
#endif//NDEBUG
  //-----------------------------------------------------------------------
  /* diameter */
  *dnew = ZERO;
  for(int ii = 0; ii < 3; ii++){
    real tmp = max[ii] - min[ii];
    if( *dnew < tmp )
      *dnew = tmp;
  }
  *dnew = LDEXP(UNITY, (int)CEIL(LOG2(*dnew)));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if 0
//-------------------------------------------------------------------------
/* works less equal than 63 (= 3 * 21) bits key */
/* based on Raman & Wise (2008), IEEE Trans. Comput., 57, 567-573 */
static inline PHint dilate3D(const PHint val)
{
  //-----------------------------------------------------------------------
  PHint ret = val;
  //-----------------------------------------------------------------------
#   if  MAXIMUM_PHKEY_LEVEL > 16
  ret = (ret * 0x100000001) & 0x7fff00000000ffff;/* ((x << 32) + x) = (2^32 + 1) * x = (16^8 + 1) * x */
#endif//MAXIMUM_PHKEY_LEVEL > 16
  ret = (ret * 0x000010001) & 0x00ff0000ff0000ff;/* ((x << 16) + x) = (2^16 + 1) * x = (16^4 + 1) * x */
  ret = (ret * 0x000000101) & 0x700f00f00f00f00f;/* ((x <<  8) + x) = (2^8  + 1) * x = (16^2 + 1) * x */
  ret = (ret * 0x000000011) & 0x30c30c30c30c30c3;/* ((x <<  4) + x) = (2^4  + 1) * x = (16^1 + 1) * x */
  ret = (ret * 0x000000005) & 0x1249249249249249;/* ((x <<  2) + x) = (2^2  + 1) * x = (       5) * x */
  //-----------------------------------------------------------------------
  return (ret);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline PHint morton3D(const PHint ix, const PHint iy, const PHint iz)
{
  //-----------------------------------------------------------------------
  return ((dilate3D(ix) << 2) | (dilate3D(iy) << 1) | (dilate3D(iz)));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void calcPHkey(const uint num, iparticle body, identity *phkey,
			     const int nlevel, const real dinv, real min[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* *keymax = (PHint)1 << nlevel; */
  /* real dscale = dinv * (real)(*keymax); */
  PHint keymax = (PHint)1 << nlevel;
  real dscale = dinv * (real)keymax;
  //-----------------------------------------------------------------------
  for(uint idx = 0; idx < num; idx++){
    PHint x = (PHint)(dscale * (body.pos[idx].x - min[0]));
    PHint y = (PHint)(dscale * (body.pos[idx].y - min[1]));
    PHint z = (PHint)(dscale * (body.pos[idx].z - min[2]));
    PHint key = 0;
    //---------------------------------------------------------------------
    for(int ii = nlevel - 1; ii >= 0; ii--){
      /* get xi, yi, zi from given position */
      PHint xi = (x >> ii) & 1;
      PHint yi = (y >> ii) & 1;
      PHint zi = (z >> ii) & 1;
      /* turn x, y, and z */
      x ^= -( xi & ((!yi) |   zi));
      y ^= -((xi & (  yi  |   zi)) | (yi & (!zi)));
      z ^= -((xi &  (!yi) & (!zi)) | (yi & (!zi)));
      /* append 3 bits to the key */
      key |= ((xi << 2) | ((xi ^ yi) << 1) | ((xi ^ zi) ^ yi)) << (3 * ii);
      /* if zi == 1, then rotate uncyclic (x->z->y->x) */
      if( zi ){
	PHint t = x;
	x = y;
	y = z;
	z = t;
      }
      else{
	/* if yi == 0, then exchange x and z */
	if( !yi ){
	  PHint t = x;
	  x = z;
	  z = t;
	}
      }
    }
    //---------------------------------------------------------------------
    phkey[idx].tag  = idx;
    /* phkey[idx].tail = (uint)( key        & 0x3fffffff);/\* pick up the  lower 30 bits of the key *\/ */
    /* phkey[idx].head = (uint)((key >> 30) & 0x3fffffff);/\* pick up the higher 30 bits of the key *\/ */
    /* /\* NOTE: 30 = 2 + 4 * 7 --> 0x3fffffff *\/ */
    phkey[idx].key = key;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
int keyAscendingOrder(const void *a, const void *b)
{
  //-----------------------------------------------------------------------
  if(          ((identity *)a)->key > ((identity *)b)->key ){    return ( 1);  }
  else{    if( ((identity *)a)->key < ((identity *)b)->key ){    return (-1);  }
    else                                                         return ( 0);  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void sortParticlesPHcurve(int num, iparticle * restrict src, iparticle * restrict dst, identity * restrict tag, PHint *peano
			  , struct timeval *start
#ifdef  EXEC_BENCHMARK
			  , wall_clock_time *elapsed
#endif//EXEC_BENCHMARK
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifndef EXEC_BENCHMARK
  gettimeofday(start, NULL);
#else///EXEC_BENCHMARK
  initStopwatch();
  *start = _benchIni;
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------
  real min[3];
  real diameter;
  calcBoxSize(num, *src, min, &diameter);
  real dinv = UNITY / diameter;
  //-----------------------------------------------------------------------
  calcPHkey(num, *src, tag, MAXIMUM_PHKEY_LEVEL, dinv, min);
  //-----------------------------------------------------------------------
  qsort(tag, num, sizeof(identity), keyAscendingOrder);
  //-----------------------------------------------------------------------
  /* sort N-body particle arrays */
  for(int ii = 0; ii < num; ii++)    (*dst).pos [ii] = (*src).pos [tag[ii].tag];
#ifdef  BLOCK_TIME_STEP
  for(int ii = 0; ii < num; ii++)    (*dst).vel [ii] = (*src).vel [tag[ii].tag];
  for(int ii = 0; ii < num; ii++)    (*dst).time[ii] = (*src).time[tag[ii].tag];
#else///BLOCK_TIME_STEP
  for(int ii = 0; ii < num; ii++)    (*dst).vx  [ii] = (*src).vx  [tag[ii].tag];
  for(int ii = 0; ii < num; ii++)    (*dst).vy  [ii] = (*src).vy  [tag[ii].tag];
  for(int ii = 0; ii < num; ii++)    (*dst).vz  [ii] = (*src).vz  [tag[ii].tag];
#endif//BLOCK_TIME_STEP
  for(int ii = 0; ii < num; ii++)    (*dst).idx [ii] = (*src).idx [tag[ii].tag];
  //-----------------------------------------------------------------------
  /* copy Peano--Hilbert key */
  for(int ii = 0; ii < num; ii++)
    peano[ii] = tag[ii].key;
#if 0
  for(int ii = 0; ii < num; ii++)
    printf("%d\t%zd\t%d\n", ii, peano[ii], tag[ii].tag);
  fflush(stdout);
  exit(1);
#endif
  //-----------------------------------------------------------------------
  /* swap the list structure */
  iparticle _tmp;
  _tmp = *src;
  *src = *dst;
  *dst = _tmp;
  //-----------------------------------------------------------------------
#ifdef  EXEC_BENCHMARK
  stopStopwatch(&(elapsed->sortParticlesPHcurve));
#endif//EXEC_BENCHMARK
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* Properties of Peano--Hilbert key of tree cells in the hierarchical structure */
//-------------------------------------------------------------------------
void initPHinfo(PHinfo *info)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* initialize fundamental properties */
  for(int ii = 0; ii < NUM_PHKEY_LEVEL; ii++)
    info[ii].level = MAXIMUM_PHKEY_LEVEL - ii;
  //-----------------------------------------------------------------------
  /* ulong ntmp = NLEAF; */
  ulong ntmp = 1;
  int jj = NUM_PHKEY_LEVEL - 1;
  while( true ){
    //---------------------------------------------------------------------
    info[jj].nmax = (int)ntmp;
    //---------------------------------------------------------------------
    ntmp *= NLEAF;    jj--;
    if( ntmp > INT_MAX )
      ntmp = INT_MAX;
    //---------------------------------------------------------------------
    if( jj < 0 )
      break;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//GENERATE_PHKEY_ON_DEVICE
//-------------------------------------------------------------------------
