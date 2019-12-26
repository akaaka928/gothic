/**
 * @file observation.h
 *
 * @brief Header file for on-the-fly analysis of NW stream
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2019/12/26 (Thu)
 *
 * Copyright (C) 2019 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef OBSERVATION_H
#define OBSERVATION_H


#include "macro.h"
#include "../misc/structure.h"


/* symmetry over the disk-plane (Northern or Southern hemisphere) */
#define NFLIP (2)

/* 112 = 7 * 16 */
/* #define NANGLE (56) */
/* #define NANGLE (112) */
/* #define NANGLE (224) */
#define NANGLE (336)

#define MAP_NX (32)
#define MAP_NY (64)


#define FILE_TAG_SCORE "score"

#define GROUP_TAG_INFO "info"
#define GROUP_TAG_DATA "data"
#define GROUP_TAG_VIEW "phi"
#define GROUP_TAG_NWS  "NWstream"

#define ATTR_TAG_SCORE "score"
#define ATTR_TAG_ANGLE "best_angle"
#define ATTR_TAG_FACE  "best_face"
#define ATTR_TAG_INDEX "best_index"
#define ATTR_TAG_GRAD "score_grad"
#define ATTR_TAG_MASS "score_mass"
#define ATTR_TAG_WIDE "score_wide"
#define ATTR_TAG_DIST "score_dist"
#define ATTR_TAG_VLOS "score_vlos"
#define ATTR_TAG_TIME "time"
#define ATTR_TAG_STEP "step"

#define DATA_TAG_MAP "map"
#define DATA_TAG_BOX "box"

#define ATTR_TAG_DM_MASS "DM_mass"
#define ATTR_TAG_DM_SIZE "DM_rs"
#define ATTR_TAG_STAR_MASS "stellar_mass"
#define ATTR_TAG_STAR_CORE "stellar_r0"
#define ATTR_TAG_STAR_EDGE "stellar_rt"
#define ATTR_TAG_ORBIT_THETA "orbit_theta"
#define ATTR_TAG_ORBIT_VRAD "orbit_vrad"
#define ATTR_TAG_ORBIT_VTAN "orbit_vtan"
#define ATTR_TAG_ORBIT_VROT "orbit_vrot"

#define ATTR_TAG_ANAL_PHI "angle_phi"
#define ATTR_TAG_ANAL_INI "initial_hemisphere"

#define ATTR_TAG_NWS_ANGLE "angle"
#define ATTR_TAG_NWS_GRAD  "gradient_normalized"
#define ATTR_TAG_NWS_SIGMA "gradient_sigma"
#define ATTR_TAG_NWS_YMIN  "gradient_ymin"
#define ATTR_TAG_NWS_YMAX  "gradient_ymax"
#define ATTR_TAG_NWS_AXIS  "axis"
#define ATTR_TAG_NWS_WIDTH "width"
#define ATTR_TAG_NWS_MASS  "minimum_mass"

#define ATTR_TAG_MAP_XMIN  "Xmin"
#define ATTR_TAG_MAP_XMAX  "Xmax"
#define ATTR_TAG_MAP_YMIN  "Ymin"
#define ATTR_TAG_MAP_YMAX  "Ymax"



/* list of functions appeared in ``observation.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__

  void set_splitter
  (char *file
#ifndef ONLINE_ANALYSIS
   , ulong *Ntot, int *unit
#endif//ONLINE_ANALYSIS
   );

  void allocate_particle_arrays_for_analysis(const ulong Ntot, nbody_aos **body, real **xi, real **eta, real **dist, real **vxi, real **veta, real **vlos);
  void  release_particle_arrays_for_analysis(                  nbody_aos  *body, real  *xi, real  *eta, real  *dist, real  *vxi, real  *veta, real  *vlos);

  void allocate_observation_arrays_for_analysis(real **map, real **box, real **score);
  void  release_observation_arrays_for_analysis(real  *map, real  *box, real  *score);

  void setNWstreamProperties(void);

  void initialize_score
  (real *score_best, const int modelID, char *file,
   const real logM_dm, const real logrs_dm,
   const real logM_star, const real logr0_star, const real logrt_star,
   const real theta, const real vr, const real vt, const real vangle);

  void finalize_score(real *score_final, const int modelID, char *file);

  void prepare_for_observation(const int num, nbody_aos *body, real disk2obs[restrict][3], real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos, const real phi, const bool flip);

  void mock_observation
  (const ulong Ntot, nbody_aos *body_anal,
#ifdef  USE_HDF5_FORMAT
   nbody_hdf5 hdf5,
#else///USE_HDF5_FORMAT
   iparticle ibody,
#endif//USE_HDF5_FORMAT
    real * restrict xi_all, real * restrict eta_all, real * restrict dist_all, real * restrict vxi_all, real * restrict veta_all, real * restrict vlos_all,
    real * restrict map_all, real * restrict box_all,
    real disk2obs[restrict][3], const real dphi,
    real * restrict score_best, const int modelID, char *file, const double time, const ulong step);

#ifdef  __CUDACC__
}
#endif//__CUDACC__



#endif//OBSERVATION_H
