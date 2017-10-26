/**
 * @file convert.c
 *
 * @brief Source code for converting type of arrays (SoA <--> AoS)
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/10/26 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "macro.h"

#include "structure.h"
#include "convert.h"


/**
 * @fn backVel
 *
 * @brief Additional function for Leap-Frog integrator.
 *
 * @param (head) head index of N-body particles
 * @param (tail) tail index of N-body particles
 * @return (body) structure contains N-body particle data (SoA)
 * @param (dt) interval of time steps
 */
void backVel(const int head, const int tail, iparticle body, real dt)
{
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
}


#ifdef  USE_HDF5_FORMAT
/**
 * @fn copyAry4Snapshot
 *
 * @brief Copy particle data for snapshot in HDF5 format.
 *
 * @param (head) head index of N-body particles
 * @param (Ni) number of N-body particles
 * @param (src) N-body particle data (for computation)
 * @return (dst) N-body particle data (for snapshot with HDF5)
 */
void copyAry4Snapshot(const int head, const int Ni, iparticle src, nbody_hdf5 dst)
{
  __NOTE__("%s\n", "start");

  for(int ii = 0; ii < Ni; ii++){
    const int jj = head + ii;

    /** position of N-body particles */
    dst.pos[jj * 3    ] = src.pos[ii].x;
    dst.pos[jj * 3 + 1] = src.pos[ii].y;
    dst.pos[jj * 3 + 2] = src.pos[ii].z;
    dst.m  [jj        ] = src.pos[ii].m;

    /** velocity of N-body particles */
#ifdef  BLOCK_TIME_STEP
    dst.vel[jj * 3    ] = src.vel[ii].x;
    dst.vel[jj * 3 + 1] = src.vel[ii].y;
    dst.vel[jj * 3 + 2] = src.vel[ii].z;
#else///BLOCK_TIME_STEP
    dst.vel[jj * 3    ] = src.vx[ii];
    dst.vel[jj * 3 + 1] = src.vy[ii];
    dst.vel[jj * 3 + 2] = src.vz[ii];
#endif//BLOCK_TIME_STEP

    /** acceleration of N-body particles */
    dst.acc[jj * 3    ] = src.acc[ii].x;
    dst.acc[jj * 3 + 1] = src.acc[ii].y;
    dst.acc[jj * 3 + 2] = src.acc[ii].z;
    dst.pot[jj        ] = src.acc[ii].pot;

    /** index of N-body particles */
    dst.idx[jj] = src.idx[ii];
  }/* for(int ii = 0; ii < Ni; ii++){ */

  __NOTE__("%s\n", "end");
}

/**
 * @fn backVel4Snapshot
 *
 * @brief Additional function for Leap-Frog integrator.
 *
 * @param (head) head index of N-body particles
 * @param (tail) tail index of N-body particles
 * @return (body) structure contains N-body particle data (SoA)
 * @param (dt) interval of time steps
 */
void backVel4Snapshot(const int head, const int tail, nbody_hdf5 body, real dt)
{
  dt *= -HALF;
  for(int ii = 3 * head; ii < 3 * tail; ii++)
    body.vel[ii] += dt * body.acc[ii];
}

#endif//USE_HDF5_FORMAT
