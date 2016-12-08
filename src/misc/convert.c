/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:33:41
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
#include "macro.h"
//-------------------------------------------------------------------------
#include "structure.h"
#include "convert.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* additional function for Leap-Frog integrator */
void backVel(const int head, const int tail, iparticle body, real dt)
{
  //-----------------------------------------------------------------------
  dt *= -HALF;
  for(int ii = head; ii < tail; ii++){
#ifdef  BLOCK_TIME_STEP
    body.vel[ii].x += dt * body.acc[ii].x;
    body.vel[ii].y += dt * body.acc[ii].y;
    body.vel[ii].z += dt * body.acc[ii].z;
#else///BLOCK_TIME_STEP
    body.vx[ii] += dt * body.acc[ii].x;
    body.vy[ii] += dt * body.acc[ii].y;
    body.vz[ii] += dt * body.acc[ii].z;
#endif//BLOCK_TIME_STEP
  }/* for(int ii = head; ii < tail; ii++){ */
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
