/**
 * @file m31coord.c
 *
 * @brief Coordinate transformation related to M31
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2018/07/23 (Mon)
 *
 * Copyright (C) 2018 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "macro.h"
#include "constants.h"
#include "rotate.h"

#include "../misc/structure.h"
#include "../anal/m31coord.h"

extern const double   length2astro;
extern const double velocity2astro;



/**
 * @fn scalarProd
 *
 * @brief Calculate scalar product.
 *
 * @param (aa) vector A
 * @param (bb) vector B
 * @return scalar product of A and B
 */
static inline real scalarProd(real aa[], real bb[])
{
  return (aa[0] * bb[0] + aa[1] * bb[1] + aa[2] * bb[2]);
}
/**
 * @fn vectorProd
 *
 * @brief Calculate vector product.
 *
 * @param (aa) vector A
 * @param (bb) vector B
 * @return (ans) vector product of A and B
 */
static inline void vectorProd(real aa[], real bb[], real ans[])
{
  ans[0] = aa[1] * bb[2] - aa[2] * bb[1];
  ans[1] = aa[2] * bb[0] - aa[0] * bb[2];
  ans[2] = aa[0] * bb[1] - aa[1] * bb[0];
}

static inline void rgemm(real aa[restrict][3], real bb[restrict][3], real ans[restrict][3])
{
  for(int ii = 0; ii < 3; ii++)
    for(int jj = 0; jj < 3; jj++)
      ans[ii][jj] = ZERO;

  for(int ii = 0; ii < 3; ii++)
    for(int kk = 0; kk < 3; kk++)
      for(int jj = 0; jj < 3; jj++)
	ans[ii][jj] += aa[ii][kk] * bb[kk][jj];
}

static inline void transpose(real ini[restrict][3], real fin[restrict][3])
{
  for(int ii = 0; ii < 3; ii++)
    for(int jj = 0; jj < 3; jj++)
      fin[jj][ii] = ini[ii][jj];
}


void setRotationMatrix(real rot[restrict][3], real inv[restrict][3])
{
  __NOTE__("%s\n", "start");

  const real theta = -37.0 * CAST_D2R(M_PI / 180.0);/**< position angle */
  const real   phi = -77.0 * CAST_D2R(M_PI / 180.0);/**< inclination */

  const real cost = COS(theta);  const real cosp = COS(phi);
  const real sint = SIN(theta);  const real sinp = SIN(phi);

  /* 1st rotation: rotation axis is z-axis, rotation angle is theta */
  static real rot1[3][3], inv1[3][3];
  static real axis1[3] = {ZERO, ZERO, UNITY};
  setRodriguesRotationMatrix(axis1, sint, cost, rot1, inv1);

  /* 2nd rotation: rotation axis is rotated y-axis (rot1 * (0, 1, 0)), rotation angle is phi */
  static real rot2[3][3], inv2[3][3];
  static real axis2[3];
  axis1[0] = ZERO;
  axis1[1] = UNITY;
  axis1[2] = ZERO;
  rotateVector(axis1, rot1, axis2);
  setRodriguesRotationMatrix(axis2, sinp, cosp, rot2, inv2);

  /* get rotation axis of M31's disk in observed frame */
  static real axis3[3], spin[3];
  spin[0] = ZERO;
  spin[1] = ZERO;
  spin[2] = -UNITY;
  rotateVector(spin, rot1, axis3);
  rotateVector(axis3, rot2, spin);

  /* rotate spin axis of M31's disk */
  axis3[0] = ZERO;
  axis3[1] = ZERO;
  axis3[2] = UNITY;
  static real rot3[3][3], inv3[3][3];
  initRotationMatrices(spin, axis3, rot3, inv3);

  /* major axis in M31's disk frame */
  static real major_disk[3] = {ZERO, UNITY, ZERO};
  static real major_tmp[3];
  rotateVector(major_disk, rot1, major_tmp);
  static real major_obs[3];
  rotateVector(major_tmp, rot2, major_obs);
  rotateVector(major_disk, inv3, major_tmp);
  static real rot4[3][3], inv4[3][3];
  static real axis[3];
  vectorProd(major_tmp, major_obs, axis);
  const real sintheta = -SQRT(scalarProd(axis, axis));
  const real costheta = scalarProd(major_tmp, major_obs);
  setRodriguesRotationMatrix(spin, sintheta, costheta, rot4, inv4);

  rgemm(rot4, inv3, inv);
  transpose(inv, rot);

  __NOTE__("%s\n", "end");
}


void standard_coordinate(const int num, nbody_aos *body, real rot[restrict][3], real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos)
{
  __NOTE__("%s\n", "start");

  const real rad2deg = CAST_D2R(180.0 * M_1_PI);
  static real ini[3], fin[3];

  for(int ii = 0; ii < num; ii++){
    /* coordinate rotation of particle position */
    ini[0] = CAST_D2R(CAST_R2D(body[ii].x) * length2astro);
    ini[1] = CAST_D2R(CAST_R2D(body[ii].y) * length2astro);
    ini[2] = CAST_D2R(CAST_R2D(body[ii].z) * length2astro);
    rotateVector(ini, rot, fin);
    const real xx = fin[0];
    const real yy = fin[1];
    const real zz = fin[2] + zm31;

    /* coordinate rotation of particle velocity */
    ini[0] = CAST_D2R(CAST_R2D(body[ii].vx) * velocity2astro);
    ini[1] = CAST_D2R(CAST_R2D(body[ii].vy) * velocity2astro);
    ini[2] = CAST_D2R(CAST_R2D(body[ii].vz) * velocity2astro);
    rotateVector(ini, rot, fin);
    const real vx = fin[0] + vm31x;
    const real vy = fin[1] + vm31y;
    const real vz = fin[2] + vm31z;

    const real tmp = UNITY / zz;
    xi [ii] = rad2deg * ATAN(xx * tmp);
    eta[ii] = rad2deg * ATAN(yy * tmp);
    const real d2 = xx * xx + yy * yy + zz * zz;
    const real dinv = RSQRT(d2);
    dist[ii] = d2 * dinv;

    vxi [ii] = (zz * vx - xx * vz) * dinv;
    veta[ii] = (zz * vy - yy * vz) * dinv;
    vlos[ii] = (xx * vx + yy * vy + zz * vz) * dinv;
  }/* for(int ii = 0; ii < num; ii++){ */

  __NOTE__("%s\n", "end");
}
