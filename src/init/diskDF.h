/**
 * @file diskDF.h
 *
 * @brief Header file for generating initial condition of disk component(s)
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
#ifndef DISKDF_H
#define DISKDF_H


#include "macro.h"
#include "rand.h"

#include "../misc/structure.h"
#include "../init/potdens.h"


/**
 * @def ENFORCE_EPICYCLIC_APPROXIMATION
 *
 * @brief Reduce velocity dispersion to ensure the consistency with epicyclic approximation
 * if switched off, then velocity dispersion is determined as GalactICS (Kuijken & Dubinski 1995; Widrow et al. 2003)
 */
/* #define ENFORCE_EPICYCLIC_APPROXIMATION */


/**
 * @def SWEEP_HIGH_ALTITUDE_COMPONENT
 *
 * @brief Remove structures above DISK_DIMMING_HEIGHT * zd
 */
#define SWEEP_HIGH_ALTITUDE_COMPONENT


#define SPEEDUP_CONVERGENCE
/* #define USE_POTENTIAL_SCALING_SCHEME */

#define SUBDIVIDE_NUM (4)


#ifdef  ENFORCE_EPICYCLIC_APPROXIMATION
#define DISK_PERP_VDISP(sigmaz, vcirc, frac) (fmin(sigmaz, (frac) * (vcirc)))
#else///ENFORCE_EPICYCLIC_APPROXIMATION
/** same method with GalactICS */
#define DISK_RADIAL_VDISP2(sz0_2, RR, invRd) ((sz0_2) *      exp(-(RR) * (invRd)))
#define DISK_RADIAL_VDISP( sz0  , RR, invRd) ((sz0  ) * sqrt(exp(-(RR) * (invRd))))
#endif//ENFORCE_EPICYCLIC_APPROXIMATION


/* list of functions appeared in ``diskDF.c'' */
void integrateSphericalDensityProfile(const int ndisk, const int maxLev, disk_data *disk);
void diffAxisymmetricPotential(const int maxLev, const disk_data disk);
void calcVerticalVdisp(const int ndisk, const int maxLev, disk_data *disk_info);
void distributeDiskParticles(ulong *Nuse, iparticle body, const real mass, const int maxLev, const disk_data disk, rand_state *rand);

void getEffectiveRadius(const int ndisk, const int maxLev, disk_data *disk);
void findIdx4nestedGrid(const double RR, const int maxLev, const disk_data disk, int * restrict lev, int * restrict idx, double * restrict alp);


#endif//DISKDF_H
