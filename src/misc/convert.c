/*************************************************************************\
 *                                                                       *
                  last updated on 2016/01/27(Wed) 16:53:05
 *                                                                       *
 *    Converting type of arrays (SoA <--> AoS)                           *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
//-------------------------------------------------------------------------
#include <macro.h>
//-------------------------------------------------------------------------
#include "structure.h"
#include "convert.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void copyStr2Ary(const int head, const int Ni, nbody_particle *src, iparticle dst)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ni; ii++){
    //---------------------------------------------------------------------
    const int jj = head + ii;
    //---------------------------------------------------------------------
    /* position of N-body particles */
    dst.pos[ii].x   = src[jj].x;
    dst.pos[ii].y   = src[jj].y;
    dst.pos[ii].z   = src[jj].z;
    dst.pos[ii].m   = src[jj].m;
    //---------------------------------------------------------------------
    /* acceleration of N-body particles */
    dst.acc[ii].x   = src[jj].ax;
    dst.acc[ii].y   = src[jj].ay;
    dst.acc[ii].z   = src[jj].az;
    dst.acc[ii].pot = src[jj].pot;
    //---------------------------------------------------------------------
    /* velocity of N-body particles */
#ifdef  BLOCK_TIME_STEP
    dst.vel [ii].x  = src[jj].vx;
    dst.vel [ii].y  = src[jj].vy;
    dst.vel [ii].z  = src[jj].vz;
    dst.vel [ii].dt = src[jj].dt;
    dst.time[ii].t0 = src[jj].t0;
    dst.time[ii].t1 = src[jj].t1;
#else///BLOCK_TIME_STEP
    dst.vx  [ii]    = src[jj].vx;
    dst.vy  [ii]    = src[jj].vy;
    dst.vz  [ii]    = src[jj].vz;
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
    /* index of N-body particles */
    dst.idx[ii]     = src[jj].idx;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void copyAry2Str(const int head, const int Ni, iparticle src, nbody_particle *dst)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ni; ii++){
    //---------------------------------------------------------------------
    const int jj = head + ii;
    //---------------------------------------------------------------------
    /* position of N-body particles */
    dst[jj].x   = src.pos[ii].x;
    dst[jj].y   = src.pos[ii].y;
    dst[jj].z   = src.pos[ii].z;
    dst[jj].m   = src.pos[ii].m;
    //---------------------------------------------------------------------
    /* velocity of N-body particles */
#ifdef  BLOCK_TIME_STEP
    dst[jj].vx  = src.vel [ii].x;
    dst[jj].vy  = src.vel [ii].y;
    dst[jj].vz  = src.vel [ii].z;
    dst[jj].dt  = src.vel [ii].dt;
    dst[jj].t0  = src.time[ii].t0;
    dst[jj].t1  = src.time[ii].t1;
#else///BLOCK_TIME_STEP
    dst[jj].vx  = src.vx [ii];
    dst[jj].vy  = src.vy [ii];
    dst[jj].vz  = src.vz [ii];
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
    /* acceleration of N-body particles */
    dst[jj].ax  = src.acc[ii].x;
    dst[jj].ay  = src.acc[ii].y;
    dst[jj].az  = src.acc[ii].z;
    dst[jj].pot = src.acc[ii].pot;
    //---------------------------------------------------------------------
    /* index of N-body particles */
    dst[jj].idx = src.idx[ii];
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* additional function for Leap-Frog integrator */
void backVel(const int head, const int tail, nbody_particle *body, real dt)
{
  //-----------------------------------------------------------------------
  dt *= -HALF;
  for(int ii = head; ii < tail; ii++){
    body[ii].vx += dt * body[ii].ax;
    body[ii].vy += dt * body[ii].ay;
    body[ii].vz += dt * body[ii].az;
  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
//-------------------------------------------------------------------------
void copyAry4Snapshot(const int head, const int Ni, iparticle src, nbody_hdf5 dst)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ni; ii++){
    //---------------------------------------------------------------------
    const int jj = head + ii;
    //---------------------------------------------------------------------
    /* position of N-body particles */
    dst.pos[jj * 3    ] = src.pos[ii].x;
    dst.pos[jj * 3 + 1] = src.pos[ii].y;
    dst.pos[jj * 3 + 2] = src.pos[ii].z;
    dst.m  [jj        ] = src.pos[ii].m;
    //---------------------------------------------------------------------
    /* velocity of N-body particles */
#ifdef  BLOCK_TIME_STEP
    dst.vel[jj * 3    ] = src.vel[ii].x;
    dst.vel[jj * 3 + 1] = src.vel[ii].y;
    dst.vel[jj * 3 + 2] = src.vel[ii].z;
#else///BLOCK_TIME_STEP
    dst.vel[jj * 3    ] = src.vx[ii];
    dst.vel[jj * 3 + 1] = src.vy[ii];
    dst.vel[jj * 3 + 2] = src.vz[ii];
#endif//BLOCK_TIME_STEP
    //---------------------------------------------------------------------
    /* acceleration of N-body particles */
    dst.acc[jj * 3    ] = src.acc[ii].x;
    dst.acc[jj * 3 + 1] = src.acc[ii].y;
    dst.acc[jj * 3 + 2] = src.acc[ii].z;
    dst.pot[jj        ] = src.acc[ii].pot;
    //---------------------------------------------------------------------
    /* index of N-body particles */
    dst.idx[jj] = src.idx[ii];
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* additional function for Leap-Frog integrator */
void backVel4Snapshot(const int head, const int tail, nbody_hdf5 body, real dt)
{
  //-----------------------------------------------------------------------
  dt *= -HALF;
  for(int ii = 3 * head; ii < 3 * tail; ii++)
    body.vel[ii] += dt * body.acc[ii];
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
