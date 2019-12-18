/**
 * @file observation.c
 *
 * @brief Source code for on-the-fly analysis of NW stream simulations instead of calling writeSnapshot in src/file/io.c
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2019/12/18 (Wed)
 *
 * Copyright (C) 2019 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

#include "macro.h"
#include "name.h"
#include "constants.h"
#include "rotate.h"

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#ifndef ONLINE_ANALYSIS
#include "myutil.h"
#include "../file/io.h"
#endif//ONLINE_ANALYSIS

#include "../misc/structure.h"
#include "../misc/allocate.h"

#include "../anal/m31coord.h"
#include "../anal/observation.h"


/** gradient analysis (based on gaussian fitting results for >= 50% filled data) */
static const real NWstream_angle = 69.0;/**< in units of degree */
/* static const real hh = 0.113172;/\**< in units of degree, mesh width for the gradient analysis; h_10 *\/ */
static const real NWstream_grad = -441.193 / 715.971;/**< NOTE: gradient must be normalized by the averaged count of particles */
/* static const real NWstream_sigma = 35.1089 / 715.971;/\**< +-3 sigma is successful, +-15 sigma is rejected *\/ */
static const real NWstream_sigma = 0.2;
static const real SCORE_UPPER_BOUND = 25.0;/**< +-3 sigma is acceptable, +-5 sigma is rejected */
static const real grad_Ymin = 3.05779;
static const real grad_Ymax = 5.8871;
/* NX = 32, NY = 64 grids -->> peak_pos = -3.62859, width_mean = 0.320655; width_mean is average of FWHM */
static const real NWstream_axis = -3.62859;
static const real NWstream_wide = 0.320655;
/* jj = 0, ..., 11; 52, ..., 63 is ZERO (Ndat <-- number of non-zero entries in the histogram) */
/* ver: min = 1.64314, max = 8.88616 */
/* hor: -5.72077, max = -2.02292 */
static const real Xmin = -5.72077;
static const real Xmax = -2.02292;
static const real Ymin = 1.64314;
static const real Ymax = 8.88616;

/* Gaussian fitting in h = 0.1 is necessary to estimating number of stream particles */
/* 3 parameters (axis, count, sigma) */


/** width analysis (based on particle counting on each grids) */
/* mesh width is 0.226344 degree; h_11 */
/* analysis for 3 degree <~ Y <~ 6 degree (13 or 14 points; but, the edge grid should be removed, perhaps) */
/* the most prominent grid is X ~ -3.5 degree or -3.7 degree */
/* the count of the above two grids are greater than neighbouring two grids in the both sides */
/* NX = 16, NY = 32 grids -->> peak_pos = -3.6356, width_mean = 0.311596 */
/* jj = 0, ..., 5; 26, ..., 31 is ZERO (Ndat <-- number of non-zero entries in the histogram) */
/* const real grad_Ymin = 3.11438; */
/* const real grad_Ymax = 6.50954; */
/* ver: min = 1.64314, max = 8.88616 */
/* hor: -5.72077, max = -2.02292 */

/** distance to the stream */
/* 4 fields x 3 boxes */
/* the middle grid must be the most abundant grid */
    /* nws_xi[0] = -5.891761379732909;    nws_eta[0] = 5.934934874532188 ;    nws_mod[0] = 24.64;    nws_ran[0] = 0.196;    nws_sys[0] = 0.057# Stream 1 */
    /* nws_xi[1] = -5.388266383737322;    nws_eta[1] = 4.426258888399156 ;    nws_mod[1] = 24.62;    nws_ran[1] = 0.186;    nws_sys[1] = 0.057# Stream 2 */
    /* nws_xi[2] = -4.942619286035076;    nws_eta[2] = 3.18150303805077  ;    nws_mod[2] = 24.59;    nws_ran[2] = 0.183;    nws_sys[2] = 0.057# Stream 3 */
    /* nws_xi[3] = -4.49814447496522 ;    nws_eta[3] = 1.9926453169853198;    nws_mod[3] = 24.58;    nws_ran[3] = 0.185;    nws_sys[3] = 0.057# Stream 4 */

static const real edgeNE_xi [2] = {-6.185539790195583, -3.9488436954717683};
static const real edgeNE_eta[2] = { 7.882569759192392,  1.6546637399150277};
static const real edgeSW_xi [2] = {-6.359861226297341, -4.066319246109101};
static const real edgeSW_eta[2] = { 6.40909029144571 , -0.22322771974637737};


#define DISTANCE_NFIELD (4)
#define DISTANCE_NDEPTH (3)
static const real distance_eta[DISTANCE_NFIELD - 1] = {5.101772155015601, 3.7516521430870844, 2.612617574895511};
static real nwsDmin[DISTANCE_NFIELD], nwsDmax[DISTANCE_NFIELD];

static const real nwsD_mod[DISTANCE_NFIELD] = {24.64, 24.62, 24.59, 24.58};
static const real nwsD_ran[DISTANCE_NFIELD] = {0.196, 0.186, 0.183, 0.185};
static const real nwsD_sys[DISTANCE_NFIELD] = {0.057, 0.057, 0.057, 0.057};

/** NE-edge: (RA, Dec) = (1.3, 48.8) -- (5.3, 42.8) */
/** SW-edge: (RA, Dec) = (1.3, 47.32) -- (5.3, 40.92) */
/* NE0 : -6.185539790195583 7.882569759192392 */
/* NE1 : -3.9488436954717683 1.6546637399150277 */
/* SW0 : -6.359861226297341 6.40909029144571 */
/* SW1 : -4.066319246109101 -0.22322771974637737 */
/* s12 : -5.616700610790968 5.101772155015598 */
/* s23 : -5.151110044528016 3.751652143087084 */
/* s34 : -4.734698074268049 2.6126175748955105 */
/* NE edge:  slope =  -2.784422092017102 , intercept =  -9.34058388367902 */
/* SW edge:  slope =  -2.8917360434134034 , intercept =  -11.981949647745681 */
/* edges of boundary between Stream fields 1 & 2: -5.186841492207897 5.101772155015601 , and -5.907773581780884 5.101772155015601 */
/* edges of boundary between Stream fields 2 & 3: -4.701958106244509 3.7516521430870844 , and -5.4408844910550105 3.7516521430870844 */
/* edges of boundary between Stream fields 3 & 4: -4.292884147430157 2.612617574895511 , and -5.046991496988008 2.612617574895511 */


    /* for ii in range(Nnws): */
    /*     nws_D[ii] = 10.0 ** (nws_mod[ii] / 5.0 - 2.0) */
    /*     nws_Dep[ii] = 10.0 ** ((nws_mod[ii] + nws_ran[ii] + nws_sys[ii]) / 5.0 - 2.0) - nws_D[ii] */
    /*     nws_Dem[ii] = nws_D[ii] - 10.0 ** ((nws_mod[ii] - nws_ran[ii] - nws_sys[ii]) / 5.0 - 2.0) */

        /* ax[              jj].plot(nws_s1_xi, nws_s1_eta, '-', color = col_nws, linewidth = lw_nws) */
        /* ax[              jj].plot(nws_s2_xi, nws_s2_eta, '-', color = col_nws, linewidth = lw_nws) */
        /* ax[              jj].plot(nws_s3_xi, nws_s3_eta, '-', color = col_nws, linewidth = lw_nws) */
        /* ax[              jj].plot(nws_s4_xi, nws_s4_eta, '-', color = col_nws, linewidth = lw_nws) */

    /* # Stream fields in Komiyama et al. (2018): Table 4 in Komiyama et al. (2018) and private communications with Komiyama-san */
    /* # Stream 1 */
    /* Nnws_s1 = 5 + 64 */
    /* nws_s1_xi  = [0] * Nnws_s1 */
    /* nws_s1_eta = [0] * Nnws_s1 */
    /* nws_s1_xi[0], nws_s1_eta[0] = -5.7055741988006705, 6.546142963104346 */
    /* nws_s1_xi[1], nws_s1_eta[1] = -5.186841492207897 , 5.101772155015601 */
    /* nws_s1_xi[2], nws_s1_eta[2] = -5.907773581780884 , 5.101772155015601 */
    /* nws_s1_xi[3], nws_s1_eta[3] = -6.239719954692477 , 6.061673446044403 */
    /* nws_s1_min = np.pi - np.arcsin((nws_s1_eta[3] - hsc_eta[2]) / hsc_fov) */
    /* nws_s1_max =         np.arccos((nws_s1_xi[0]  - hsc_xi[2] ) / hsc_fov) */
    /* dtheta = (nws_s1_max - nws_s1_min) / 64 */
    /* for ii in range(64): */
    /*     nws_s1_xi[4 + ii]  = hsc_xi[2]  + hsc_fov * np.cos(nws_s1_min + float(ii) * dtheta) */
    /*     nws_s1_eta[4 + ii] = hsc_eta[2] + hsc_fov * np.sin(nws_s1_min + float(ii) * dtheta) */
    /* # Stream 2 */
    /* Nnws_s2 = 5 */
    /* nws_s2_xi  = [0] * Nnws_s2 */
    /* nws_s2_eta = [0] * Nnws_s2 */
    /* nws_s2_xi[0], nws_s2_eta[0] = -5.186841492207897 , 5.101772155015601 */
    /* nws_s2_xi[1], nws_s2_eta[1] = -4.701958106244509 , 3.7516521430870844 */
    /* nws_s2_xi[2], nws_s2_eta[2] = -5.4408844910550105, 3.7516521430870844 */
    /* nws_s2_xi[3], nws_s2_eta[3] = -5.907773581780884 , 5.101772155015601 */
    /* # Stream 3 */
    /* Nnws_s3 = 5 */
    /* nws_s3_xi  = [0] * Nnws_s3 */
    /* nws_s3_eta = [0] * Nnws_s3 */
    /* nws_s3_xi[0], nws_s3_eta[0] = -4.701958106244509 , 3.7516521430870844 */
    /* nws_s3_xi[1], nws_s3_eta[1] = -4.292884147430157 , 2.612617574895511 */
    /* nws_s3_xi[2], nws_s3_eta[2] = -5.046991496988008 , 2.612617574895511 */
    /* nws_s3_xi[3], nws_s3_eta[3] = -5.4408844910550105, 3.7516521430870844 */
    /* # Stream 4 */
    /* Nnws_s4 = 5 + 64 */
    /* nws_s4_xi  = [0] * Nnws_s4 */
    /* nws_s4_eta = [0] * Nnws_s4 */
    /* nws_s4_xi[0], nws_s4_eta[0] = -4.6822224216441475, 1.5578016922010893 */
    /* nws_s4_xi[1], nws_s4_eta[1] = -5.046991496988008 , 2.612617574895511 */
    /* nws_s4_xi[2], nws_s4_eta[2] = -4.292884147430157 , 2.612617574895511 */
    /* nws_s4_xi[3], nws_s4_eta[3] = -4.113795100125332 , 2.1139580751416602 */
    /* nws_s4_min = 2.0 * np.pi + np.arcsin((nws_s4_eta[3] - hsc_eta[3]) / hsc_fov) */
    /* nws_s4_max = 2.0 * np.pi - np.arccos((nws_s4_xi[0]  - hsc_xi[3] ) / hsc_fov) */
    /* dtheta = (nws_s4_max - nws_s4_min) / 64 */
    /* for ii in range(64): */
    /*     nws_s4_xi[4 + ii]  = hsc_xi[3]  + hsc_fov * np.cos(nws_s4_min + float(ii) * dtheta) */
    /*     nws_s4_eta[4 + ii] = hsc_eta[3] + hsc_fov * np.sin(nws_s4_min + float(ii) * dtheta) */

    /* nws_s1_xi[Nnws_s1 - 1] = nws_s1_xi[0];    nws_s1_eta[Nnws_s1 - 1] = nws_s1_eta[0] */
    /* nws_s2_xi[Nnws_s2 - 1] = nws_s2_xi[0];    nws_s2_eta[Nnws_s2 - 1] = nws_s2_eta[0] */
    /* nws_s3_xi[Nnws_s3 - 1] = nws_s3_xi[0];    nws_s3_eta[Nnws_s3 - 1] = nws_s3_eta[0] */
    /* nws_s4_xi[Nnws_s4 - 1] = nws_s4_xi[0];    nws_s4_eta[Nnws_s4 - 1] = nws_s4_eta[0] */


/** minimum mass of the NW stream */
/* at least 10^5.5 Msun in the stream regions */
/* estimation by Carlberg et al. (2011): 2.2e+6 Msun (M/L = 3, PAndAS data) */
/* M/L = 1 corresponds to the lower limit: 7.3e+5 Msun */
/* mismatch of region -->> 7.3e+5 / 2 = 3.6e+5 Msun (half is very pessimistic estimation) */
/* Carlberg assumes distance to NW stream is equal to that to M31 -->> underestimate the absolute magnitude of NW stream (but, maybe slightly) */
static real NWstream_mass;


/* error bar for xi and eta?? <-- within 0.113 degree in radius? */
/* at least one particle matches with the line-of-sight velocity of the observed GC */
/** PAndAS-04 is outside of the Subaru/HSC fieles */
#define NGC_INTERNAL (4)
/** PAndAS-09, PAndAS-10, PAndAS-11, PAndAS-12 */
#define NGC_EXTERNAL (1)
/** PAndAS-04 */
static const real gc_xi  [NGC_INTERNAL + NGC_EXTERNAL] = {-5.260959997771581, -5.124223955350093, -4.946180555295868 , -4.557799679219338 ,           -6.436733608263167 };
static const real gc_eta [NGC_INTERNAL + NGC_EXTERNAL] = { 4.061434922115992,  4.137743013677189,  3.5546551505544537,  2.2086054244275655,            6.4598886185030535};
static const real gc_vmin[NGC_INTERNAL + NGC_EXTERNAL] = {-444.0 - 21.0, -435.0 - 10.0, -447.0 - 13.0, -472.0 - 5.0,		-397.0 - 7.0};
static const real gc_vmax[NGC_INTERNAL + NGC_EXTERNAL] = {-444.0 + 21.0, -435.0 + 10.0, -447.0 + 13.0, -472.0 + 5.0,		-397.0 + 7.0};





#define NFIELD_HSC (5)
/** f003, f004, f009, f022, f023 */
const real hsc_xi [NFIELD_HSC] = {-4.653744156858361, -5.7227710294896745, -5.530826768687439, -4.842739589150247, -5.912437203950868};
const real hsc_eta[NFIELD_HSC] = { 3.635683362752261,  4.47091006650441  ,  5.816784796173433,  2.290423176281251,  3.125050968157858};
const real hsc_fov = (45.0 - 2.0) / 60.0;



      /* const double inv_sqrt2sig = M_SQRT1_2 / sigma; */
      /* const double coeff_gauss = 0.5 * amp; */
      /* const double coeff_slope = 0.5 * slope; */

	/* const double val = hist[kk] * dx; */
	/* if( val > 0.0 ){ */
	/*   const double diff = coeff_gauss * (erf(inv_sqrt2sig * (xmax - mu)) - erf(inv_sqrt2sig * (xmin - mu))) - val; */
	/*   /\* const double diff = dx * (noise + coeff_slope * (xmin + xmax)) - val; *\/ */

	/*   /\* score += diff * diff; *\/ */
	/*   score += diff * diff * (isig2[kk] / (dx * dx)); */
	/* }/\* if( val > 0.0 ){ *\/ */














#ifdef  __ICC
/* Disable ICC's remark #161: unrecognized #pragma */
#     pragma warning (disable:161)
#endif//__ICC
int idxAscendingOrder(const void *a, const void *b);
int idxAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((const nbody_aos *)a)->idx > ((const nbody_aos *)b)->idx ){    return ( 1);  }
  else{    if( ((const nbody_aos *)a)->idx < ((const nbody_aos *)b)->idx ){    return (-1);  }
    else{                                                                      return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}
#ifdef  __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC


void set_splitter
(int *kind, int **bodyHead, int **bodyNum, char *file
#ifndef ONLINE_ANALYSIS
 , ulong *Ntot, int *unit
#endif//ONLINE_ANALYSIS
 )
{
  __NOTE__("%s\n", "start");

  FILE *fp;
  char filename[128];
  sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "r");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }

  int unit_tmp;
  bool checker = true;
  checker &= (1 == fscanf(fp, "%d", &unit_tmp));
  checker &= (1 == fscanf(fp, "%d\t%*d", kind));
  *bodyHead = (int *)malloc(sizeof(int) * (*kind));  if( *bodyHead == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate bodyHead");  }
  *bodyNum  = (int *)malloc(sizeof(int) * (*kind));  if( *bodyNum  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate bodyNum");  }
  for(int ii = 0; ii < (*kind); ii++)
    checker &= (1 == fscanf(fp, "%d", &((*bodyNum)[ii])));
  fclose(fp);
  if( !checker ){
    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
}
  (*bodyHead)[0] = 0;
  for(int ii = 1; ii < (*kind); ii++)
    (*bodyHead)[ii] = (*bodyHead)[ii - 1] + (*bodyNum)[ii - 1];

#ifndef ONLINE_ANALYSIS
  *unit = unit_tmp;
  *Ntot = 0;
  for(int ii = 0; ii < (*kind); ii++)
    *Ntot += (*bodyNum)[ii];
#endif//ONLINE_ANALYSIS

  __NOTE__("%s\n", "end");
}


void allocate_particle_arrays_for_analysis(const ulong Ntot, nbody_aos **body, real **xi, real **eta, real **dist, real **vxi, real **veta, real **vlos)
{
  __NOTE__("%s\n", "start");

  static int Nthreads = 1;
#ifdef  CPUS_PER_PROCESS
  Nthreads = CPUS_PER_PROCESS;
#endif//CPUS_PER_PROCESS
  __NOTE__("Nthreads = %d\n", Nthreads);

  *body = (nbody_aos *)malloc(sizeof(nbody_aos) * Ntot);
  if( *body == NULL ){    __KILL__(stderr, "ERROR: failure to allocate body\n");  }

  *xi   = (real *)malloc(sizeof(real) * Nthreads * Ntot);  if( *  xi == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xi.");  }
  *eta  = (real *)malloc(sizeof(real) * Nthreads * Ntot);  if( * eta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate eta.");  }
  *dist = (real *)malloc(sizeof(real) * Nthreads * Ntot);  if( *dist == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate dist.");  }
  *vxi  = (real *)malloc(sizeof(real) * Nthreads * Ntot);  if( *vxi  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vxi.");  }
  *veta = (real *)malloc(sizeof(real) * Nthreads * Ntot);  if( *veta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate veta.");  }
  *vlos = (real *)malloc(sizeof(real) * Nthreads * Ntot);  if( *vlos == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vlos.");  }

  __NOTE__("%s\n", "end");
}


void release_particle_arrays_for_analysis(nbody_aos *body, real *xi, real *eta, real *dist, real *vxi, real *veta, real *vlos)
{
  __NOTE__("%s\n", "start");

  free(body);

  free(xi);
  free(eta);
  free(dist);
  free(vxi);
  free(veta);
  free(vlos);

  __NOTE__("%s\n", "end");
}


void allocate_observation_arrays_for_analysis(real **map, real **box, real **score)
{
  __NOTE__("%s\n", "start");

  static int Nthreads = 1;
#ifdef  CPUS_PER_PROCESS
  Nthreads = CPUS_PER_PROCESS;
#endif//CPUS_PER_PROCESS
  __NOTE__("Nthreads = %d\n", Nthreads);

  *map = (real *)malloc(sizeof(real) * Nthreads * MAP_NY * MAP_NX);  if( *map == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate map.");  }

  *box = (real *)malloc(sizeof(real) * Nthreads * DISTANCE_NFIELD * DISTANCE_NDEPTH);  if( *box == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate box.");  }

  *score = (real *)malloc(sizeof(real) * NANGLE);  if( *score == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate score.");  }


  /** calculate distance to the NW stream */
  for(int ff = 0; ff < DISTANCE_NFIELD; ff++){
    /* const real distance = POW10(0.2 * nwsD_mod[ff] - 2.0); */
    nwsDmin[ff] = POW(TEN, 0.2 * (nwsD_mod[ff] - nwsD_ran[ff] - nwsD_sys[ff]) - TWO);
    nwsDmax[ff] = POW(TEN, 0.2 * (nwsD_mod[ff] + nwsD_ran[ff] + nwsD_sys[ff]) - TWO);
  }

  /** set the minimum mass of the NW stream */
  extern const double mass_astro2com;
  NWstream_mass = POW(TEN, 5.5) * mass_astro2com;


  __NOTE__("%s\n", "end");
}


void release_observation_arrays_for_analysis(real *map, real *box, real *score)
{
  __NOTE__("%s\n", "start");

  free(map);
  free(box);
  free(score);

  __NOTE__("%s\n", "end");
}


void initialize_score
(real *score_best, const int modelID, char *file,
 const real logM_dm, const real logrs_dm,
 const real logM_star, const real logr0_star, const real logrt_star,
 const real theta, const real vr, const real vt, const real vangle)
{
  __NOTE__("%s\n", "start");


  const real worst_score = SCORE_UPPER_BOUND * 5 * 2;/**< 5 is # of criterion, 2 is just a safety parameter indicating that analysis is not yet performed */
  for(int ii = 0; ii < NANGLE; ii++)
    score_best[ii] = worst_score;

  char filename[256];
  sprintf(filename, "%s/%s_%s%d.h5", DATAFOLDER, file, FILE_TAG_SCORE, modelID);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#ifdef  USE_DOUBLE_PRECISION
  const hid_t H5T_GOTHIC_REAL = H5T_NATIVE_DOUBLE;
#else///USE_DOUBLE_PRECISION
  const hid_t H5T_GOTHIC_REAL = H5T_NATIVE_FLOAT;
#endif//USE_DOUBLE_PRECISION


  hid_t group;
  hsize_t attr_dims;
  hid_t dspc, attr;
  /** write various information of the simulation */
  group = H5Gcreate(target, GROUP_TAG_INFO, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  attr_dims = 1;
  dspc = H5Screate_simple(1, &attr_dims, NULL);

  attr = H5Acreate(group, ATTR_TAG_DM_MASS, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &logM_dm));
  chkHDF5err(H5Aclose(attr));
  attr = H5Acreate(group, ATTR_TAG_DM_SIZE, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &logrs_dm));
  chkHDF5err(H5Aclose(attr));

  attr = H5Acreate(group, ATTR_TAG_STAR_MASS, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &logM_star));
  chkHDF5err(H5Aclose(attr));
  attr = H5Acreate(group, ATTR_TAG_STAR_CORE, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &logr0_star));
  chkHDF5err(H5Aclose(attr));
  attr = H5Acreate(group, ATTR_TAG_STAR_EDGE, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &logrt_star));
  chkHDF5err(H5Aclose(attr));

  attr = H5Acreate(group, ATTR_TAG_ORBIT_THETA, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &theta));
  chkHDF5err(H5Aclose(attr));
  attr = H5Acreate(group, ATTR_TAG_ORBIT_VRAD, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &vr));
  chkHDF5err(H5Aclose(attr));
  attr = H5Acreate(group, ATTR_TAG_ORBIT_VTAN, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &vt));
  chkHDF5err(H5Aclose(attr));
  attr = H5Acreate(group, ATTR_TAG_ORBIT_VROT, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &vangle));
  chkHDF5err(H5Aclose(attr));

  attr = H5Acreate(group, ATTR_TAG_SCORE, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &worst_score));
  chkHDF5err(H5Aclose(attr));
  const real tmp_phi = TWO * M_PI;
  attr = H5Acreate(group, ATTR_TAG_ANGLE, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &tmp_phi));
  chkHDF5err(H5Aclose(attr));
  const int tmp_idx = NANGLE;
  attr = H5Acreate(group, ATTR_TAG_INDEX, H5T_NATIVE_INT, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_NATIVE_INT, &tmp_idx));
  chkHDF5err(H5Aclose(attr));

  chkHDF5err(H5Sclose(dspc));
  chkHDF5err(H5Gclose(group));


  /** preparation to save analyzed results */
  group = H5Gcreate(target, GROUP_TAG_DATA, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  attr_dims = 1;
  dspc = H5Screate_simple(1, &attr_dims, NULL);

  const int Nangle = NANGLE;
  attr = H5Acreate(group, "NANGLE", H5T_NATIVE_INT, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &Nangle));
  chkHDF5err(H5Aclose(attr));

  const int Nx = MAP_NX;
  attr = H5Acreate(group, "MAP_NX", H5T_NATIVE_INT, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &Nx));
  chkHDF5err(H5Aclose(attr));
  const int Ny = MAP_NY;
  attr = H5Acreate(group, "MAP_NY", H5T_NATIVE_INT, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &Ny));
  chkHDF5err(H5Aclose(attr));

  const int Nfield = DISTANCE_NFIELD;
  attr = H5Acreate(group, "DISTANCE_NFIELD", H5T_NATIVE_INT, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &Nfield));
  chkHDF5err(H5Aclose(attr));
  const int Ndepth = DISTANCE_NDEPTH;
  attr = H5Acreate(group, "DISTANCE_NDEPTH", H5T_NATIVE_INT, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &Ndepth));
  chkHDF5err(H5Aclose(attr));

  for(int ii = 0; ii < NANGLE; ii++){
      char grp[16];
      sprintf(grp, "%s_%d", GROUP_TAG_VIEW, ii);
      hid_t subgrp = H5Gcreate(group, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      attr = H5Acreate(subgrp, ATTR_TAG_SCORE, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &worst_score));
      chkHDF5err(H5Aclose(attr));
      attr = H5Acreate(subgrp, ATTR_TAG_GRAD, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &worst_score));
      chkHDF5err(H5Aclose(attr));
      attr = H5Acreate(subgrp, ATTR_TAG_MASS, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &worst_score));
      chkHDF5err(H5Aclose(attr));
      attr = H5Acreate(subgrp, ATTR_TAG_WIDE, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &worst_score));
      chkHDF5err(H5Aclose(attr));
      attr = H5Acreate(subgrp, ATTR_TAG_DIST, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &worst_score));
      chkHDF5err(H5Aclose(attr));
      attr = H5Acreate(subgrp, ATTR_TAG_VLOS, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &worst_score));
      chkHDF5err(H5Aclose(attr));
      const double time = 0.0;
      attr = H5Acreate(subgrp, ATTR_TAG_TIME, H5T_NATIVE_DOUBLE, dspc, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_NATIVE_DOUBLE, &time));
      chkHDF5err(H5Aclose(attr));
      const ulong step = 0;
      attr = H5Acreate(subgrp, ATTR_TAG_STEP, H5T_NATIVE_ULONG, dspc, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_NATIVE_ULONG, &step));
      chkHDF5err(H5Aclose(attr));

      const real phi = (real)ii * TWO * M_PI / (real)NANGLE;
      attr = H5Acreate(subgrp, ATTR_TAG_ANAL_PHI, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &phi));
      chkHDF5err(H5Aclose(attr));

      chkHDF5err(H5Gclose(subgrp));
  }

  chkHDF5err(H5Sclose(dspc));
  chkHDF5err(H5Gclose(group));


  /** note properties of the NW stream */
  group = H5Gcreate(target, GROUP_TAG_NWS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  attr_dims = 1;
  dspc = H5Screate_simple(1, &attr_dims, NULL);

  attr = H5Acreate(group, ATTR_TAG_NWS_ANGLE, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &NWstream_angle));
  chkHDF5err(H5Aclose(attr));

  attr = H5Acreate(group, ATTR_TAG_NWS_GRAD, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &NWstream_grad));
  chkHDF5err(H5Aclose(attr));
  attr = H5Acreate(group, ATTR_TAG_NWS_SIGMA, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &NWstream_sigma));
  chkHDF5err(H5Aclose(attr));
  attr = H5Acreate(group, ATTR_TAG_NWS_YMIN, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &grad_Ymin));
  chkHDF5err(H5Aclose(attr));
  attr = H5Acreate(group, ATTR_TAG_NWS_YMAX, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &grad_Ymax));
  chkHDF5err(H5Aclose(attr));

  attr = H5Acreate(group, ATTR_TAG_NWS_AXIS, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &NWstream_axis));
  chkHDF5err(H5Aclose(attr));
  attr = H5Acreate(group, ATTR_TAG_NWS_WIDTH, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &NWstream_wide));
  chkHDF5err(H5Aclose(attr));

  chkHDF5err(H5Sclose(dspc));
  chkHDF5err(H5Gclose(group));


  chkHDF5err(H5Fclose(target));


  __NOTE__("%s\n", "end");
}


void finalize_score(real *score_final, const int modelID, char *file)
{
  __NOTE__("%s\n", "start");

  real best_score = SCORE_UPPER_BOUND * 5 * 2;/**< 5 is # of criterion, 2 is just a safety parameter indicating that analysis is not yet performed */
  int best_index = NANGLE;

  for(int ii = 0; ii < NANGLE; ii++){
    const real score = score_final[ii];
    if( score < best_score ){
      best_score = score;
      best_index = ii;
    }
  }


  char filename[256];
  sprintf(filename, "%s/%s_%s%d.h5", DATAFOLDER, file, FILE_TAG_SCORE, modelID);
  hid_t target = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
#ifdef  USE_DOUBLE_PRECISION
  const hid_t H5T_GOTHIC_REAL = H5T_NATIVE_DOUBLE;
#else///USE_DOUBLE_PRECISION
  const hid_t H5T_GOTHIC_REAL = H5T_NATIVE_FLOAT;
#endif//USE_DOUBLE_PRECISION

  hid_t group = H5Gopen(target, GROUP_TAG_INFO, H5P_DEFAULT);
  hid_t attr;
  attr = H5Aopen(group, ATTR_TAG_SCORE, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &best_score));
  chkHDF5err(H5Aclose(attr));
  const real best_angle = (real)best_index * TWO * M_PI / (real)NANGLE;
  attr = H5Aopen(group, ATTR_TAG_ANGLE, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &best_angle));
  chkHDF5err(H5Aclose(attr));
  attr = H5Aopen(group, ATTR_TAG_INDEX, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_NATIVE_INT, &best_index));
  chkHDF5err(H5Aclose(attr));
  chkHDF5err(H5Gclose(group));

  group = H5Gopen(target, GROUP_TAG_NWS, H5P_DEFAULT);
  hsize_t attr_dims = 1;
  hid_t dspc = H5Screate_simple(1, &attr_dims, NULL);
  attr = H5Acreate(group, ATTR_TAG_NWS_MASS, H5T_GOTHIC_REAL, dspc, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &NWstream_mass));
  chkHDF5err(H5Aclose(attr));
  chkHDF5err(H5Sclose(dspc));
  chkHDF5err(H5Gclose(group));

  chkHDF5err(H5Fclose(target));


  __NOTE__("%s\n", "end");
}


void prepare_for_observation(const int num, nbody_aos *body, real disk2obs[restrict][3], real * restrict xi, real * restrict eta, real * restrict dist, real * restrict vxi, real * restrict veta, real * restrict vlos, const real phi)
{
  __NOTE__("%s\n", "start");

  const real rad2deg = CAST_D2R(180.0 * M_1_PI);
  real ini[3], fin[3];

  extern const double length2astro, velocity2astro;

  real view[3][3];
  const real cosp = COS(phi);
  const real sinp = SIN(phi);
  view[0][0] = cosp;  view[0][1] = -sinp;  view[0][2] =	 ZERO;
  view[1][0] = sinp;  view[1][1] =  cosp;  view[1][2] =	 ZERO;
  view[2][0] = ZERO;  view[2][1] =  ZERO;  view[2][2] = UNITY;

  real rot[3][3];
  for(int ii = 0; ii < 3; ii++)
    for(int jj = 0; jj < 3; jj++)
      rot[ii][jj] = ZERO;
  for(int ii = 0; ii < 3; ii++)
    for(int kk = 0; kk < 3; kk++)
      for(int jj = 0; jj < 3; jj++)
	rot[ii][jj] += disk2obs[ii][kk] * view[kk][jj];


  for(int ii = 0; ii < num; ii++){
    /** coordinate rotation of particle position */
    ini[0] = CAST_D2R(CAST_R2D(body[ii].x) * length2astro);
    ini[1] = CAST_D2R(CAST_R2D(body[ii].y) * length2astro);
    ini[2] = CAST_D2R(CAST_R2D(body[ii].z) * length2astro);
    rotateVector(ini, rot, fin);
    const real xx = fin[0];
    const real yy = fin[1];
    const real zz = fin[2] + zm31;

    /** coordinate rotation of particle velocity */
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


static inline real normalize_score(const real score, const int score_max, const real max_value)
{
  return ((score / (real)score_max) * max_value);
}

void mock_observation
  (const ulong Ntot, nbody_aos *body_anal, int * restrict bodyHead, int * restrict bodyNum,
   real * restrict xi_all, real * restrict eta_all, real * restrict dist_all, real * restrict vxi_all, real * restrict veta_all, real * restrict vlos_all,
   real * restrict map_all, real * restrict box_all,
   real disk2obs[restrict][3], const real dphi,
   real * restrict score_best, const int modelID, char *file, const double time, const ulong step)
{
  __NOTE__("%s\n", "start");

  char filename[256];
  sprintf(filename, "%s/%s_%s%d.h5", DATAFOLDER, file, FILE_TAG_SCORE, modelID);
  hid_t target = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
#ifdef  USE_DOUBLE_PRECISION
  const hid_t H5T_GOTHIC_REAL = H5T_NATIVE_DOUBLE;
#else///USE_DOUBLE_PRECISION
  const hid_t H5T_GOTHIC_REAL = H5T_NATIVE_FLOAT;
#endif//USE_DOUBLE_PRECISION




  const real hsc_fov2 = hsc_fov * hsc_fov;

  const real NWtheta = NWstream_angle * M_PI / 180.0;
  const real cost = COS(HALF * M_PI - NWtheta);
  const real sint = SIN(HALF * M_PI - NWtheta);

  const real dXinv = (real)MAP_NX / (Xmax - Xmin);
  const int nws_imid = (int)FLOOR(dXinv * (NWstream_axis - Xmin));
  const int nws_nwid = (int)NEARBYINT(NWstream_wide * dXinv);
  int nws_imin, nws_imax;
  if( (nws_nwid & 1) ){
    nws_imin = nws_imid - ((nws_nwid - 1) >> 1);
    nws_imax = nws_imid + ((nws_nwid - 1) >> 1);
  }
  else{
    const int nws_next = (int)NEARBYINT(dXinv * (NWstream_axis - Xmin));
    nws_imin = nws_imid - (nws_nwid >> 1) + (nws_next == nws_imid);
    nws_imax = nws_imin + nws_nwid - 1;
  }

  const real dYinv = (real)MAP_NY / (Ymax - Ymin);
  const real dY = (Ymax - Ymin) / (real)MAP_NY;
  const int nws_jmin = (int) CEIL(dYinv * (grad_Ymin - Ymin));
  const int nws_jmax = (int)FLOOR(dYinv * (grad_Ymax - Ymin));
  const int nws_jnum = nws_jmax - nws_jmin + 1;

  const real gc_r2 = ((Xmax - Xmin) / (real)MAP_NX) * ((Ymax - Ymin) / (real)MAP_NY);


  const int dm_head = bodyHead[0];
  const int dm_num  = bodyNum [0];
  const int star_head = bodyHead[1];
  const int star_num  = bodyNum [1];


  /* xi = xi0 + slope * eta */
  const real slopeNE = (edgeNE_xi[1] - edgeNE_xi[0]) / (edgeNE_eta[1] - edgeNE_eta[0]);
  const real iceptNE = edgeNE_xi[0] - slopeNE * edgeNE_eta[0];
  const real slopeSW = (edgeSW_xi[1] - edgeSW_xi[0]) / (edgeSW_eta[1] - edgeSW_eta[0]);
  const real iceptSW = edgeSW_xi[0] - slopeSW * edgeSW_eta[0];




  /** index sort of the particle array */
  qsort(body_anal, Ntot, sizeof(nbody_aos), idxAscendingOrder);


  /** data analysis using OpenMP: parallelization about viewing angle */
#pragma omp parallel for num_threads(CPUS_PER_PROCESS)
  for(int pp = 0; pp < NANGLE; pp++){
    const int threadIdx = omp_get_thread_num();
    const real phi = dphi * (real)pp;

    real *pi_xi  ;    pi_xi   = &  xi_all[threadIdx * Ntot];
    real *pi_eta ;    pi_eta  = & eta_all[threadIdx * Ntot];
    real *pi_dist;    pi_dist = &dist_all[threadIdx * Ntot];
    real *pi_vxi ;    pi_vxi  = & vxi_all[threadIdx * Ntot];
    real *pi_veta;    pi_veta = &veta_all[threadIdx * Ntot];
    real *pi_vlos;    pi_vlos = &vlos_all[threadIdx * Ntot];

    /** coordinate transformation to observed frame */
    prepare_for_observation(Ntot, body_anal, disk2obs, pi_xi, pi_eta, pi_dist, pi_vxi, pi_veta, pi_vlos, phi);



    /** initialize score of the model */
    int Nmatch_vlos_gc_internal = 0;
    int Nmatch_vlos_gc_external = 0;
    int mismatch_vlos_gc[NGC_INTERNAL + NGC_EXTERNAL];
    for(int ii = 0; ii < NGC_INTERNAL + NGC_EXTERNAL; ii++)
      mismatch_vlos_gc[ii] = 1;/**< 1 means that no particle in the specified field satisfies the criterion, 0 means that there is at least one particle that satisfies the criterion */

    real *map;
    map = &map_all[threadIdx * MAP_NY * MAP_NX];
    for(int ii = 0; ii < MAP_NY * MAP_NX; ii++)
      map[ii] = ZERO;

    real *box;
    box = &box_all[threadIdx * DISTANCE_NFIELD * DISTANCE_NDEPTH];
    for(int ii = 0; ii < DISTANCE_NFIELD * DISTANCE_NDEPTH; ii++)
      box[ii] = ZERO;


    /* pick up particles in Subaru HSC fields, follow analysis for Komiyama-san's data (for stellar particles) */
    for(int ii = star_head; ii < star_head + star_num; ii++){
      const real  xi = pi_xi [ii];
      const real eta = pi_eta[ii];
      const real XX =  cost * xi + sint * eta;
      const real YY = -sint * xi + cost * eta;

      /** pick up particles locate in Subaru/HSC fields and/or its surroundings */
      if( (XX > Xmin) && (XX < Xmax) && (YY > Ymin) && (YY < Ymax) ){
	/** generate mass distribution map */
	const real mi = body_anal[ii].m;
	map[INDEX2D(MAP_NY, MAP_NX, (int)FLOOR(dYinv * (YY - Ymin)), (int)FLOOR(dXinv * (XX - Xmin)))] += mi;

	/** pick up particles locate in Subaru/HSC fields */
	bool include = false;
	for(int ff = 0; ff < NFIELD_HSC; ff++){
	  const real dx =  xi - hsc_xi [ff];
	  const real dy = eta - hsc_eta[ff];

	  include |= ((dx * dx + dy * dy) <= hsc_fov2);

	  if( include )
	    break;
	}

	/** the candidate particle locate in Subaru/HSC fields */
	if( include ){
	  /** distance analysis */
	  if( (xi >= (iceptSW + slopeSW * eta)) && (xi <= (iceptNE + slopeNE * eta)) ){
	    int ff = DISTANCE_NFIELD - 1;
	    if      ( eta >= distance_eta[0] )	      ff = 0;
	    else if ( eta >= distance_eta[1] )	      ff = 1;
	    else if ( eta >= distance_eta[2] )	      ff = 2;

	    const real dist = pi_dist[ii];
	    if      ( dist < nwsDmin[ff] )	      box[INDEX2D(DISTANCE_NFIELD, DISTANCE_NDEPTH, ff, 0)] += mi;
	    else if ( dist > nwsDmax[ff] )	      box[INDEX2D(DISTANCE_NFIELD, DISTANCE_NDEPTH, ff, 2)] += mi;
	    else	                              box[INDEX2D(DISTANCE_NFIELD, DISTANCE_NDEPTH, ff, 1)] += mi;

	  }

	  /** line-of-sight velocity analysis for stellar particles (gc_internal) */
	  if( Nmatch_vlos_gc_internal < NGC_INTERNAL )
	    for(int ff = 0; ff < NGC_INTERNAL; ff++)
	      if( mismatch_vlos_gc[ff] ){
		const real dx =  xi - gc_xi [ff];
		const real dy = eta - gc_eta[ff];

		if( (dx * dx + dy * dy) <= gc_r2 ){
		  const real vlos = pi_vlos[ii];
		  if( (vlos >= gc_vmin[ff]) && (vlos <= gc_vmax[ff]) ){
		    mismatch_vlos_gc[ff] = 0;
		    Nmatch_vlos_gc_internal++;
		  }
		}
	      }

	}
      }
      else{
	/** line-of-sight velocity analysis for stellar particles (gc_external) */
	if( Nmatch_vlos_gc_external < NGC_EXTERNAL )
	  for(int ff = NGC_INTERNAL; ff < NGC_INTERNAL + NGC_EXTERNAL; ff++)
	    if( mismatch_vlos_gc[ff] ){
	      const real dx =  xi - gc_xi [ff];
	      const real dy = eta - gc_eta[ff];

	      if( (dx * dx + dy * dy) <= gc_r2 ){
		const real vlos = pi_vlos[ii];
		if( (vlos >= gc_vmin[ff]) && (vlos <= gc_vmax[ff]) ){
		  mismatch_vlos_gc[ff] = 0;
		  Nmatch_vlos_gc_external++;
		}
	      }
	    }
      }


    }


    /** line-of-sight velocity analysis for DM particles (when no stellar particle satisfies the condition) */
    Nmatch_vlos_gc_internal += Nmatch_vlos_gc_external;
    if( Nmatch_vlos_gc_internal < (NGC_INTERNAL + NGC_EXTERNAL) )
      for(int ii = dm_head; ii < dm_head + dm_num; ii++){
	const real  xi = pi_xi [ii];
	const real eta = pi_eta[ii];

	for(int ff = 0; ff < NGC_INTERNAL + NGC_EXTERNAL; ff++)
	  if( mismatch_vlos_gc[ff] ){
	    const real dx =  xi - gc_xi [ff];
	    const real dy = eta - gc_eta[ff];

	    if( (dx * dx + dy * dy) <= gc_r2 ){
	      const real vlos = pi_vlos[ii];
	      if( (vlos >= gc_vmin[ff]) && (vlos <= gc_vmax[ff]) ){
		mismatch_vlos_gc[ff] = 0;
		Nmatch_vlos_gc_internal++;
	      }
	    }
	  }

	if( Nmatch_vlos_gc_internal == (NGC_INTERNAL + NGC_EXTERNAL) )
	  break;
      }




    /** check the stream like structure (using XY map) */
    /** check the width of the stream [0 : nws_jnum] */
    /** check the mass of the stream [0 or 1] */
    /** check the gradient of the surface-density profile */
    real mass_stream = ZERO;
    int wrong_width = 0;
    real Sx = ZERO;
    real Sy = ZERO;
    real Sxx = ZERO;
    real Sxy = ZERO;
    for(int jj = nws_jmin; jj < nws_jmax + 1; jj++){
      real west = ZERO;
      for(int ii = nws_imin - nws_nwid; ii < nws_imin; ii++)
	west += map[INDEX2D(MAP_NY, MAP_NX, ii, jj)];

      real core = ZERO;
      real peak = ZERO;
      int exist = 0;
      for(int ii = nws_imin; ii < nws_imax + 1; ii++){
	const real entry = map[INDEX2D(MAP_NY, MAP_NX, ii, jj)];
	core += entry;

	peak = FMAX(peak, entry);
	if( entry >= HALF * peak )
	  exist++;
      }
      /** check the mass of the stream */
      mass_stream += core;

      real east = ZERO;
      for(int ii = nws_imax + 1; ii < nws_imax + 1 + nws_nwid; ii++)
	east += map[INDEX2D(MAP_NY, MAP_NX, ii, jj)];


      /** check the gradient of the surface-density profile */
      /** assume sigma_i is a constant for all points */
      const real ypos = Ymin + dY * (HALF + (real)jj);
      Sx += ypos;
      Sy += core;
      Sxx += ypos * ypos;
      Sxy += ypos * core;

      /** check the width of the stream */
      core *= HALF;
      wrong_width += ((exist <= (nws_nwid >> 1)) || (west >= core) || (east >= core));
    }
    /** check the width of the stream */
    const real score_wide = (real)wrong_width;

    /** check the mass of the stream */
    const real score_mass = (real)(mass_stream < NWstream_mass);

    /** check the gradient of the surface-density profile */
    const real grad = (Sxy - Sx * Sy) / (Sxx - Sx * Sx);
    const real mean = Sy / (real)nws_jnum;
    const real grad_diff = ((grad / mean) - NWstream_grad) / NWstream_sigma;
    const real score_grad = FMIN(grad_diff * grad_diff, SCORE_UPPER_BOUND);


    /** check the distance analysis [0 : DISTANCE_NFIELD] */
    int iscore_dist = 0;
    for(int ff = 0; ff < DISTANCE_NFIELD; ff++){
      const real count = HALF * box[INDEX2D(DISTANCE_NFIELD, DISTANCE_NDEPTH, ff, 1)];
      iscore_dist += ((box[INDEX2D(DISTANCE_NFIELD, DISTANCE_NDEPTH, ff, 0)] >= count) || (box[INDEX2D(DISTANCE_NFIELD, DISTANCE_NDEPTH, ff, 2)] >= count));
    }
    const real score_dist = (real)iscore_dist;

    /** check the line-of-sight velocity [0 : NGC_INTERNAL + NGC_EXTERNAL] */
    const real score_vlos = (real)(NGC_INTERNAL + NGC_EXTERNAL - Nmatch_vlos_gc_internal);


    /** sum-up total score */
    const real score =  score_grad
      + normalize_score(score_mass,                           1, SCORE_UPPER_BOUND)
      + normalize_score(score_wide,                    nws_jnum, SCORE_UPPER_BOUND)
      + normalize_score(score_dist,             DISTANCE_NFIELD, SCORE_UPPER_BOUND)
      + normalize_score(score_vlos, NGC_INTERNAL + NGC_EXTERNAL, SCORE_UPPER_BOUND)
      ;


    /** if the score is good, update the status and dump to the HDF5 file */
    const bool update_score = score < score_best[pp];
#pragma omp critical (update)
    if( update_score ){
      char grp[16];
      sprintf(grp, "%s/%s_%d", GROUP_TAG_DATA, GROUP_TAG_VIEW, pp);
      hid_t group = H5Gopen(target, grp, H5P_DEFAULT);

      /** update score as attribute (tentative results or initial value are written in the file) */
      hid_t attr;
      attr = H5Aopen(group, ATTR_TAG_SCORE, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &score));
      chkHDF5err(H5Aclose(attr));
      attr = H5Aopen(group, ATTR_TAG_GRAD, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &score_grad));
      chkHDF5err(H5Aclose(attr));
      attr = H5Aopen(group, ATTR_TAG_MASS, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &score_mass));
      chkHDF5err(H5Aclose(attr));
      attr = H5Aopen(group, ATTR_TAG_WIDE, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &score_wide));
      chkHDF5err(H5Aclose(attr));
      attr = H5Aopen(group, ATTR_TAG_DIST, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &score_dist));
      chkHDF5err(H5Aclose(attr));
      attr = H5Aopen(group, ATTR_TAG_VLOS, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_GOTHIC_REAL, &score_vlos));
      chkHDF5err(H5Aclose(attr));
      attr = H5Aopen(group, ATTR_TAG_TIME, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_NATIVE_DOUBLE, &time));
      chkHDF5err(H5Aclose(attr));
      attr = H5Aopen(group, ATTR_TAG_STEP, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attr, H5T_NATIVE_ULONG, &step));
      chkHDF5err(H5Aclose(attr));


      /** update data (tentative results or NULL) */
      if( H5Lexists(group, DATA_TAG_MAP, H5P_DEFAULT) ){
	hid_t dset = H5Dopen(group, DATA_TAG_MAP, H5P_DEFAULT);
	chkHDF5err(H5Dwrite(dset, H5T_GOTHIC_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, map));
	chkHDF5err(H5Dclose(dset));
      }
      else{
	hsize_t dims[2] = {MAP_NY, MAP_NX};
	hid_t dspace = H5Screate_simple(2, dims, NULL);
	hid_t dset = H5Dcreate(group, DATA_TAG_MAP, H5T_GOTHIC_REAL, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	chkHDF5err(H5Dwrite(dset, H5T_GOTHIC_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, map));
	chkHDF5err(H5Dclose(dset));
	chkHDF5err(H5Sclose(dspace));
      }

      if( H5Lexists(group, DATA_TAG_BOX, H5P_DEFAULT) ){
	hid_t dset = H5Dopen(group, DATA_TAG_BOX, H5P_DEFAULT);
	chkHDF5err(H5Dwrite(dset, H5T_GOTHIC_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, box));
	chkHDF5err(H5Dclose(dset));
      }
      else{
	hsize_t dims[2] = {DISTANCE_NFIELD, DISTANCE_NDEPTH};
	hid_t dspace = H5Screate_simple(2, dims, NULL);
	hid_t dset = H5Dcreate(group, DATA_TAG_BOX, H5T_GOTHIC_REAL, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	chkHDF5err(H5Dwrite(dset, H5T_GOTHIC_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, box));
	chkHDF5err(H5Dclose(dset));
	chkHDF5err(H5Sclose(dspace));
      }

      chkHDF5err(H5Gclose(group));
    }
    score_best[pp] = score;
  }


  chkHDF5err(H5Fclose(target));


  __NOTE__("%s\n", "end");
}



/** test module for debugging analyzer */
#ifndef ONLINE_ANALYSIS
int main(int argc, char **argv)
{
#ifndef ONLINE_ANALYSIS
  if( argc < 15 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 15);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -modelID=<int>\n");
    __FPRINTF__(stderr, "          -logM_dm=<real> -logrs_dm=<real>\n");
    __FPRINTF__(stderr, "          -logM_star=<real> -logr0_star=<real> -logrt_star=<real>\n");
    __FPRINTF__(stderr, "          -theta=<real> -vrad=<real> -vtan=<real> -vangle=<real>\n");
    __FPRINTF__(stderr, "          -eps=<real>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }
  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));

  int modelID;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "modelID", &modelID));
  real logM_dm;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "logM_dm", &logM_dm));
  real logrs_dm;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "logrs_dm", &logrs_dm));
  real logM_star;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "logM_star", &logM_star));
  real logr0_star;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "logr0_star", &logr0_star));
  real logrt_star;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "logrt_star", &logrt_star));
  real theta;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "theta", &theta));
  real vrad;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vrad", &vrad));
  real vtan;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vtan", &vtan));
  real vangle;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vangle", &vangle));

  real   eps;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "eps", &eps));


  ulong Ntot;
  int   unit;

#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
  nbody_hdf5 hdf5;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  real *hdf5_acc_ext, *hdf5_pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  ulong *hdf5_idx;
  allocSnapshotArray
    (&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     &hdf5_acc_ext, &hdf5_pot_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     (int)Ntot, &hdf5);
#else///USE_HDF5_FORMAT
  iparticle ibody;
  ulong *idx;
  position *pos;
  acceleration *acc;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  acceleration *ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  velocity *vel;
  ibody_time *ti;
#else///BLOCK_TIME_STEP
  real *vx, *vy, *vz;
#endif//BLOCK_TIME_STEP
  allocParticleData
    ((int)Ntot, &ibody, &idx, &pos, &acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     &ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
     &vel, &ti
#else///BLOCK_TIME_STEP
     &vx, &vy, &vz
#endif//BLOCK_TIME_STEP
     );
#endif//USE_HDF5_FORMAT
#endif//ONLINE_ANALYSIS

  /** all values above are stored in N-body simulations or required only reading snapshot files */




  /** variables must be defined in the main code */
  static real disk2obs[3][3], obs2disk[3][3];/**< rotation matrix (disk frame <--> observed frame) */
  nbody_aos *body_anal;/**< to split DM and stellar particles by sorting particle array */

  real *xi_all, *eta_all, *dist_all, *vxi_all, *veta_all, *vlos_all;/**< for OpenMP */
  real *map_all, *box_all;/**< for OpenMP */
  real *score_all;

  int kind;
  int *bodyHead, *bodyNum;



  /** initialization of the analysis */
  /** read number of components */
  set_splitter(&kind, &bodyHead, &bodyNum, file
#ifndef ONLINE_ANALYSIS
	       , &Ntot, &unit
#endif//ONLINE_ANALYSIS
	       );
#ifndef ONLINE_ANALYSIS
  setPhysicalConstantsAndUnitSystem(unit, 1);
#endif//ONLINE_ANALYSIS

  /** memory allocation by function call, which supports OpenMP */
  allocate_particle_arrays_for_analysis(Ntot, &body_anal, &xi_all, &eta_all, &dist_all, &vxi_all, &veta_all, &vlos_all);
  allocate_observation_arrays_for_analysis(&map_all, &box_all, &score_all);


  /** set coordinate transformation from M31's disk coordinate to observerd frame */
  setRotationMatrix(obs2disk, disk2obs);


  /** initialization for the analysis */
  initialize_score(score_all, modelID, file,
		   logM_dm, logrs_dm,
		   logM_star, logr0_star, logrt_star,
		   theta, vrad, vtan, vangle);

  const real dphi = TWO * M_PI / (real)NANGLE;

  /** read snapshot instead of on-the-fly calling of writeSnapshot() */
#ifndef ONLINE_ANALYSIS
  for(int ifile = start; ifile < end + 1; ifile += interval){
    /** read snapshot */
    double time;
    ulong steps;
    int unit_read;
#ifdef  USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file, (uint)ifile, &hdf5, hdf5type);
    for(int ii = 0; ii < (int)Ntot; ii++){
      body_anal[ii]. x  = hdf5.pos[ii * 3];      body_anal[ii]. y = hdf5.pos[ii * 3 + 1];      body_anal[ii].z   = hdf5.pos[ii * 3 + 2];
      body_anal[ii].vx  = hdf5.vel[ii * 3];      body_anal[ii].vy = hdf5.vel[ii * 3 + 1];      body_anal[ii].vz  = hdf5.vel[ii * 3 + 2];
      body_anal[ii].ax  = hdf5.acc[ii * 3];      body_anal[ii].ay = hdf5.acc[ii * 3 + 1];      body_anal[ii].az  = hdf5.acc[ii * 3 + 2];
      body_anal[ii].idx = hdf5.idx[ii    ];      body_anal[ii]. m = hdf5.  m[ii        ];      body_anal[ii].pot = hdf5.pot[ii        ];
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      body_anal[ii].ax_ext = hdf5.acc_ext[ii * 3];      body_anal[ii].ay_ext = hdf5.acc_ext[ii * 3 + 1];      body_anal[ii].az_ext = hdf5.acc_ext[ii * 3 + 2];
      body_anal[ii].pot_ext = hdf5.pot_ext[ii];
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#else///USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file, ibody, (uint)ifile);
    for(int ii = 0; ii < (int)Ntot; ii++){
      body_anal[ii]. x  = ibody.pos[ii].x;      body_anal[ii]. y = ibody.pos[ii].y;      body_anal[ii]. z  = ibody.pos[ii].z;
      body_anal[ii].vx  = ibody.vel[ii].x;      body_anal[ii].vy = ibody.vel[ii].y;      body_anal[ii].vz  = ibody.vel[ii].z;
      body_anal[ii].ax  = ibody.acc[ii].x;      body_anal[ii].ay = ibody.acc[ii].y;      body_anal[ii].az  = ibody.acc[ii].z;
      body_anal[ii].idx = ibody.idx[ii]  ;      body_anal[ii]. m = ibody.pos[ii].m;      body_anal[ii].pot = ibody.acc[ii].pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      body_anal[ii].ax_ext = ibody.acc_ext[ii].x;      body_anal[ii].ay_ext = ibody.acc_ext[ii].y;      body_anal[ii].az_ext = ibody.acc_ext[ii].z;
      body_anal[ii].pot_ext = ibody.acc_ext[ii].pot;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#endif//USE_HDF5_FORMAT
    if( unit_read != unit ){
      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
    }/* if( unit_read != unit ){ */
#endif//ONLINE_ANALYSIS


    mock_observation(Ntot, body_anal, bodyHead, bodyNum,
		     xi_all, eta_all, dist_all, vxi_all, veta_all, vlos_all,
		     map_all, box_all,
		     disk2obs, dphi,
		     score_all, modelID, file, time, steps);


#ifndef ONLINE_ANALYSIS
  }
#endif//ONLINE_ANALYSIS


  finalize_score(score_all, modelID, file);


  free(bodyHead);
  free(bodyNum);

  release_particle_arrays_for_analysis(body_anal, xi_all, eta_all, dist_all, vxi_all, veta_all, vlos_all);
  release_observation_arrays_for_analysis(map_all, box_all, score_all);



#ifndef ONLINE_ANALYSIS
#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
  freeSnapshotArray
    (hdf5_pos, hdf5_vel, hdf5_acc, hdf5_m, hdf5_pot, hdf5_idx
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , hdf5_acc_ext, hdf5_pot_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
#else///USE_HDF5_FORMAT
  freeParticleData
    (idx, pos, acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
     vel, ti
#else///BLOCK_TIME_STEP
     vx, vy, vz
#endif//BLOCK_TIME_STEP
     );
#endif//USE_HDF5_FORMAT
#endif//ONLINE_ANALYSIS



  return (0);
}
#endif//ONLINE_ANALYSIS
