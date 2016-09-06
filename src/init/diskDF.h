/*************************************************************************\
 *                                                                       *
                  last updated on 2016/09/06(Tue) 15:22:13
 *                                                                       *
 *    Header File for Definition to generate initial condition of disk   *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef DISKDF_H
#define DISKDF_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MACRO_H
#include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef STRUCTURE_H
#include "../misc/structure.h"
#endif//STRUCTURE_H
//-------------------------------------------------------------------------
#ifndef POTDENS_H
#include "../init/potdens.h"
#endif//POTDENS_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define USE_ORIGINAL_VDISP_ESTIMATOR
#define SPEEDUP_CONVERGENCE
/* #define USE_POTENTIAL_SCALING_SCHEME */
//-------------------------------------------------------------------------
#define SUBDIVIDE_NUM (4)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_ORIGINAL_VDISP_ESTIMATOR
//-------------------------------------------------------------------------
#define DISK_PERP_VDISP(sigmaz, vcirc, frac) (fmin(sigmaz, (frac) * (vcirc)))
//-------------------------------------------------------------------------
#else///USE_ORIGINAL_VDISP_ESTIMATOR
//-------------------------------------------------------------------------
/* same method with GalactICS */
#define DISK_RADIAL_VDISP2(sz0_2, RR, invRd) ((sz0_2) *      exp(-(RR) * (invRd)))
#define DISK_RADIAL_VDISP( sz0  , RR, invRd) ((sz0  ) * sqrt(exp(-(RR) * (invRd))))
//-------------------------------------------------------------------------
#endif//USE_ORIGINAL_VDISP_ESTIMATOR
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "diskDF.c"
//-------------------------------------------------------------------------
void integrateSphericalDensityProfile(const int ndisk, const int maxLev, disk_data *disk);
//-------------------------------------------------------------------------
void diffAxisymmetricPotential(const int maxLev, const disk_data disk);
//-------------------------------------------------------------------------
void calcVerticalVdisp(const int ndisk, const int maxLev, disk_data *disk_info);
//-------------------------------------------------------------------------
void distributeDiskParticles(ulong *Nuse, iparticle body, const real mass, const int maxLev, const disk_data disk);
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//DISKDF_H
//-------------------------------------------------------------------------
