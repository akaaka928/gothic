/**
 * @file analyze.c
 *
 * @brief Analyzer of Subaru/HSC data using CMA-ES (Evolution Strategy with Covariance Matrix Adaptation)
 *
 * @author Yohei MIKI (University of Tokyo)
 *
 * @date 2019/12/26 (Thu)
 *
 * Copyright (C) 2019 Yohei MIKI
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */


#define LONG_GRAD_FIT


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#ifndef DISABLE_MPI
#include <mpi.h>
#include "mpilib.h"
#endif//DISABLE_MPI

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#define RESULTS "profile"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "name.h"
#include "myutil.h"
#include "rand.h"

#ifdef  USE_SFMTJUMP
#include "SFMT-jump.h"
#include "sfmtjump_polynomial.h"
#include <omp.h>
#endif//USE_SFMTJUMP

#include "../find/cmaes.h"
#include "../find/cmaes_io.h"

/**
 * @def USE_GZIP_COMPRESSION
 *
 * @brief On to enable gzip compression for HDF5 files (default is ON).
 *
 * @detail currently, h5py does not accept Szip compression in default.
 */
/* #define USE_GZIP_COMPRESSION */
#ifdef  USE_GZIP_COMPRESSION
/* The maximum number of elements in a chunk is 2^32-1 which is equal to 4,294,967,295 */
/* The maximum size for any chunk is 4GB */
#define MAXIMUM_CHUNK_SIZE      ((hsize_t)1 << 31)
#define MAXIMUM_CHUNK_SIZE_4BIT ((hsize_t)1 << 30)
#define MAXIMUM_CHUNK_SIZE_8BIT ((hsize_t)1 << 29)
#endif//USE_GZIP_COMPRESSION


/* number of parameters: (0) center of gaussian, (1) width of gaussian, (2) height of gaussian, (3) slope of background/foreground, and (4) amount of background/foreground; minize the chi-square of the profile fitting */
#define NDIM (5)
/* static inline void evaluate(struct cmaes_status *cma, const struct cmaes_config cfg, const int numx, const int numy, const int yj, double * restrict dat, double * restrict hor) */
static inline void evaluate(struct cmaes_status *cma, const struct cmaes_config cfg, const int numx, double * restrict hist, double * restrict isig2, double * restrict hor)
{
  if( cma->state == CMAES_STATE_SAMPLED ){
    const int ndim = cfg.ndim;
    const int lambda = cfg.lambda;

    for(int ii = 0; ii < lambda; ii++){
      const double    mu = cma->xx[INDEX2D(lambda, ndim, ii, 0)];/**< location of the stream axis */
      const double sigma = cma->xx[INDEX2D(lambda, ndim, ii, 1)];/**< width of the stream */
      const double   amp = cma->xx[INDEX2D(lambda, ndim, ii, 2)];/**< number count of stars */
      const double slope = cma->xx[INDEX2D(lambda, ndim, ii, 3)];/**< slope of foreground/background stars */
      const double noise = cma->xx[INDEX2D(lambda, ndim, ii, 4)];/**< floor of foreground/background stars */

      const double inv_sqrt2sig = M_SQRT1_2 / sigma;
      const double coeff_gauss = 0.5 * amp;
      const double coeff_slope = 0.5 * slope;

      double score = 0.0;
      for(int kk = 0; kk < numx; kk++){
	const double xmin = hor[kk];
	const double xmax = hor[kk + 1];
	const double dx = xmax - xmin;
	/* const double val = dat[INDEX2D(numx, numy, kk, yj)] * dx; */
	const double val = hist[kk] * dx;
	if( val > 0.0 ){
	  const double diff = dx * (noise + coeff_slope * (xmin + xmax)) + coeff_gauss * (erf(inv_sqrt2sig * (xmax - mu)) - erf(inv_sqrt2sig * (xmin - mu))) - val;
	  /* const double diff = dx * (noise + coeff_slope * (xmin + xmax)) - val; */

	  /* score += diff * diff; */
	  score += diff * diff * (isig2[kk] / (dx * dx));
	}/* if( val > 0.0 ){ */
      }/* for(int kk = 0; kk < numx; kk++){ */

      cma->score[ii] = score;
    }/* for(int ii = 0; ii < lambda; ii++){ */

    cma->state = CMAES_STATE_SCORED;
  }/* if( cma->state == CMAES_STATE_SAMPLED ){ */
}


/**
 * @fn optimize
 *
 * @brief optimize parameter space
 * @detail see cmaes_Optimize in dm:~/work/181025cmaes/c-cmaes/src/cmaes.c: line 758--792
 *
n * @return (cmaes) status of CMA-ES
 */
static inline bool optimize(struct cmaes_status *cma, struct cmaes_config cfg, struct cmaes_rand *gauss, rand_state *rand, const int numx, double * restrict hist, double * restrict isig2, double * restrict hor)
{
  __NOTE__("%s\n", "start");

  const int ndim = cfg.ndim;
  const int lambda = cfg.lambda;


  int ll = 1;
  while( true ){
    if( hist[ll] > 0.0 )
      break;
    ll++;
  }/* while( true ){ */
  int rr = numx - 2;
  while( true ){
    if( hist[rr] > 0.0 )
      break;
    rr--;
  }/* while( true ){ */
  const double xmin = 0.5 * (hor[ll] + hor[ll - 1]);
  const double xmax = 0.5 * (hor[rr] + hor[rr + 1]);

  const double dx = hor[1] - hor[0];

  /** find the most luminous grid */
  double hmax = 0.0;
  int hidx = 0;
  for(int ii = ll; ii < rr + 1; ii++)
    if( hist[ii] > hmax ){
      hmax = hist[ii];
      hidx = ii;
    }

  while( true ){
    /** generate new search points */
    bool success = sampling(cma, cfg, gauss, rand);
    if( !success ){
	goto detect_error;
	break;
      }

    /** enforce positivity on some parameters */
    for(int ii = 0; ii < lambda; ii++){
      while( true ){
	const double    mu = cma->xx[INDEX2D(lambda, ndim, ii, 0)];/**< location of the stream axis */
	const double sigma = cma->xx[INDEX2D(lambda, ndim, ii, 1)];/**< width of the stream */
	const double  star = cma->xx[INDEX2D(lambda, ndim, ii, 2)];/**< number count of stars */
	const double slope = cma->xx[INDEX2D(lambda, ndim, ii, 3)];/**< slope of foreground/background stars */
	const double noise = cma->xx[INDEX2D(lambda, ndim, ii, 4)];/**< floor of foreground/background stars */

	/** peak of the stream should be located between xmin and xmax (mu > xmin, mu < xmax) */
	/** width of the stream must be a non-zero positive value (sigma > 0.0) */
	/** contribution from the stream must not be a negative value (star >= 0.0) */
	/** contribution from the stream must not exceed total number of star counts (star <= 1.0) */
	/** noise level must not be a negative value in the fields (noise + slope * xmin >= 0.0, noise + slope * xmax >= 0.0) */
	/** width of the stream must be a narrow value [smaller than 0.3 degree in Kirihara et al. (2017)] (sigma < fmax(0.5, 0.25 * (xmax - xmin))) */

	/** additional constraints: */
	/** width of the stream should be greater than or equal to 1/4 of the grid width */
	/** edge meshes should not include the stream axis except for the case that they are the most luminous grid */

	if( (mu > xmin) && (mu < xmax) &&
	    (sigma > 0.0) &&
	    (star >= 0.0) &&
	    (star <= 1.0) &&
	    (noise + slope * xmin >= 0.0) && (noise + slope * xmax >= 0.0) &&
	    /* (sigma < 0.5) && (sigma >= 0.25 * dx) && */
	    (sigma >= 0.25 * dx) &&
	    (((mu > xmin + dx) && (mu < xmax - dx)) || ((mu < xmin + dx) && (ll == hidx)) || ((mu > xmax - dx) && (rr == hidx)))
	    )
	  break;
	resampling(cma, cfg, gauss, rand, ii);
      }/* while( true ){ */
    }/* for(int ii = 0; ii < lambda; ii++){ */


    /** evaluate the new search points */
    /* evaluate(cma, cfg, numx, numy, yj, dat, hor); */
    evaluate(cma, cfg, numx, hist, isig2, hor);

    /** update the search distribution */
    update(cma, cfg);

    /** check whether the solution is found */
    if( terminate(cma, cfg) )
      break;
  }/* while( true ){ */

  __NOTE__("%s\n", "end");
  return (true);

      detect_error:
  fprintf(stderr, "thread ID %d aborted CMA-ES precedure\n", omp_get_thread_num());
  return (false);
}


static inline void rotate(const real xx, const real yy, real * restrict XX, real * restrict YY, const real cost, const real sint)
{
  *XX =  cost * xx + sint * yy;
  *YY = -sint * xx + cost * yy;
}


static inline void reportFittingResults
(const double theta, const int mem, const real width, const real grad_ymin, const real grad_ymax, const real grad_avg, const real aa, const real sa, const real bb, const real sb, const real chisq, const real rchisq, const real aa_peak, const real sa_peak, const real chisq_peak, const real rchisq_peak,
 const int numx, real * restrict hor, const int numy, real * restrict ver, real * restrict dat, real * restrict peak_x, real * restrict fwhm, real * restrict sdens, real * restrict sigma, real * restrict slope, real * restrict noise, real * restrict fchisq, int * restrict Ndat,
#ifdef  USE_HDF5_FORMAT
 const hid_t H5T_NATIVE_REAL, const hid_t sub
#else///USE_HDF5_FORMAT
 char *file, const int skip, real * restrict peak_y, real * restrict peak_y, real * restrict grad, real * restrict num
#endif//USE_HDF5_FORMAT
 )
{
  /** dump results */
  fprintf(stdout, "%e\t%e\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", theta * M_1_PI * 180.0, ver[1] - ver[0], mem, width, aa, sa, bb, sb, rchisq, aa_peak, sa_peak, rchisq_peak, FABS(aa_peak) + FABS(sa_peak));


  if( mem > 2 ){
#ifdef  USE_HDF5_FORMAT
    /** create sub-directory */
    hid_t fit = H5Gcreate(sub, "fit", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t attr_dims = 1;
    hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    attribute = H5Acreate(fit, "grad_chisq", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &chisq));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "peak_chisq", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &chisq_peak));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "grad_rchisq", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &rchisq));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "peak_rchisq", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &rchisq_peak));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "grad_entry_average", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &grad_avg));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "grad_slope", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &aa));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "grad_slope_sigma", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &sa));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "grad_intercept", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &bb));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "grad_intercept_sigma", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &sb));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "grad_ymin", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &grad_ymin));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "grad_ymax", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &grad_ymax));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "peak_pos", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &aa_peak));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "peak_pos_sigma", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &sa_peak));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "width_mean", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &width));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(fit, "sampling_points", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &mem));
    chkHDF5err(H5Aclose(attribute));
    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(fit));
#else///USE_HDF5_FORMAT
    sprintf(filename, "%s/%s_prof_theta%.3f_hh%.2u.dat", DATAFOLDER, file, deg, ilog2(skip));
    fp = fopen(filename, "w");
    if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
    fprintf(fp, "# theta = %e, hh = %e\n", deg, ver[1] - ver[0]);
    fprintf(fp, "# gradient = %e \\pm %e\n", aa, sa);
    fprintf(fp, "# intercept = %e \\pm %e\n", bb, sb);
    fprintf(fp, "# chi square = %e\n", chisq);
    fprintf(fp, "# reduced chi square = %e\n", rchisq);
    fprintf(fp, "#peak_x(deg)\tpeak_y(deg)\tdens(deg^-2)\tfwhm(deg)\tslope(deg^-1)\tnoise\tgrad(deg^-3)\n");
    for(int jj = 0; jj < numy; jj++)
      fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", peak_x[jj], peak_y[jj], sdens[jj], fwhm[jj], slope[jj], noise[jj], grad[jj]);
    fclose(fp);

    sprintf(filename, "%s/%s_map_theta%.3f_hh%.2u.dat", DATAFOLDER, file, deg, ilog2(skip));
    fp = fopen(filename, "w");
    if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
    fprintf(fp, "#X(deg.)\tY(deg.)\tSigma(/deg.^2)\tnum\n");
    for(int ii = 0; ii < numx; ii++){
      for(int jj = 0; jj < numy; jj++)
	fprintf(fp, "%e\t%e\t%e\t%e\n", hor[ii], ver[jj], dat[INDEX2D(numx, numy, ii, jj)], num[INDEX2D(numx, numy, ii, jj)]);
      fprintf(fp, "\n");
    }/* for(int ii = 0; ii < numx; ii++){ */
    fclose(fp);
#endif//USE_HDF5_FORMAT
  }/* if( mem > 2 ){ */


  /** dump results */
#ifdef  USE_HDF5_FORMAT
  /** prepare data compression */
  hsize_t dims[2] = {numx, numy};
  hid_t dataspace = H5Screate_simple(1, &dims[1], NULL);
  hid_t property;
#ifdef  USE_GZIP_COMPRESSION
  const hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  const uint gzip_compress_level = 9;/**< 9 is the maximum compression ratio */
  hsize_t cdims_loc[2];
  const hsize_t gzip_cdims[2] = {1, 1024};
  cdims_loc[1] = gzip_cdims[1];
  property = H5Pcreate(H5P_DATASET_CREATE);
  if( dims[1] < cdims_loc[1] )
    cdims_loc[1] = dims[1];
  if( cdims_loc[1] > cdims_max )
    cdims_loc[1] = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc[1]));
  chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#else///USE_GZIP_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_GZIP_COMPRESSION

  hid_t dataset;
  dataset = H5Dcreate(sub, "peak", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, peak_x));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(sub, "fwhm", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, fwhm));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(sub, "sigma", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, sigma));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(sub, "sdens", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, sdens));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(sub, "slope", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, slope));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(sub, "noise", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, noise));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(sub, "fchisq", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, fchisq));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(sub, "Ndat", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Ndat));
  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));

  dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
  cdims_loc[0] = gzip_cdims[0];
  cdims_loc[1] = gzip_cdims[1];
  property = H5Pcreate(H5P_DATASET_CREATE);
  if( dims[1] < cdims_loc[1] ){
    cdims_loc[1] = dims[1];
    cdims_loc[0] = dims[1] / cdims_loc[1];
  }/* if( dims[1] < cdims_loc[1] ){ */
  if( dims[0] < cdims_loc[0] )
    cdims_loc[0] = dims[0];
  if( cdims_loc[0] * cdims_loc[1] > cdims_max )
    cdims_loc[0] = cdims_max / cdims_loc[1];
  chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
  chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
  dataset = H5Dcreate(sub, "dat", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat));
  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));

  dims[0] = numx + 1;
  dataspace = H5Screate_simple(1, &dims[0], NULL);
#ifdef  USE_GZIP_COMPRESSION
  cdims_loc[0] = gzip_cdims[1];
  property = H5Pcreate(H5P_DATASET_CREATE);
  if( dims[0] < cdims_loc[0] )
    cdims_loc[0] = dims[0];
  if( cdims_loc[0] > cdims_max )
    cdims_loc[0] = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc[0]));
  chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
  dataset = H5Dcreate(sub, "hor", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, hor));
  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));

  dims[0] = numy + 1;
  dataspace = H5Screate_simple(1, &dims[0], NULL);
#ifdef  USE_GZIP_COMPRESSION
  cdims_loc[0] = gzip_cdims[1];
  property = H5Pcreate(H5P_DATASET_CREATE);
  if( dims[0] < cdims_loc[0] )
    cdims_loc[0] = dims[0];
  if( cdims_loc[0] > cdims_max )
    cdims_loc[0] = cdims_max;
  chkHDF5err(H5Pset_chunk(property, 1, &cdims_loc[0]));
  chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
  dataset = H5Dcreate(sub, "ver", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, ver));
  chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
  chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
  chkHDF5err(H5Sclose(dataspace));
#endif//USE_HDF5_FORMAT
}


static inline int least_square_method
(const int ini, const int fin, bool cont[restrict], real hor[restrict], real peak_x[restrict], real peak_y[restrict], real sdens[restrict], real sigma[restrict], real fwhm[restrict], real error[restrict],
 real * restrict slope, real * restrict slope_sigma, real * restrict icept, real * restrict icept_sigma, real * restrict chisq, real * restrict rchisq, real * restrict mean_count,
 real * restrict axis, real * restrict axis_sigma, real * restrict chisq_axis, real * restrict rchisq_axis, real * restrict mean_width)
{
  real S = ZERO;
  real Sx = ZERO;
  real Sy = ZERO;
  real xsum = ZERO;
  real S_peak = ZERO;
  real Sy_peak = ZERO;
  real Ss_peak = ZERO;
  const real sinv_peak = UNITY / (HALF * (hor[1] - hor[0]));
  real width = ZERO;
  real lpos_old = peak_x[ini] - sigma[ini];
  real rpos_old = peak_x[ini] + sigma[ini];

  int mem = 0;
  real rho_min = ZERO;
  for(int jj = ini; jj < fin; jj++)
    if( sdens[jj] > ZERO ){
      rho_min += sdens[jj];
      mem++;
    }
  rho_min = QUARTER * QUARTER * (rho_min / (real)mem);

  mem = 0;
  for(int jj = ini; jj < fin; jj++){
    real lpos_new = peak_x[jj] - sigma[jj];
    real rpos_new = peak_x[jj] + sigma[jj];
    cont[jj] = false;
    if( (sdens[jj] > rho_min) && (rpos_old >= lpos_new) && (rpos_new >= lpos_old) ){
      lpos_old = lpos_new;
      rpos_old = rpos_new;
      cont[jj] = true;
      mem++;

      const real tmp = UNITY / (error[jj] * error[jj]);
      S  += tmp;
      Sx += tmp * peak_y[jj];
      Sy += tmp * sdens[jj];
      xsum += peak_y[jj];
      const real tmp_peak = sinv_peak * sinv_peak;
      S_peak  += tmp_peak;
      Sy_peak += tmp_peak * peak_x[jj];
      Ss_peak += sinv_peak;
      width += fwhm[jj];
    }
  }

  real Stt = ZERO;
  real sum = ZERO;
  real St  = ZERO;
  real avg = ZERO;
  const real invS = UNITY / S;
  for(int jj = ini; jj < fin; jj++)
    if( cont[jj] ){
      const real sinv = UNITY / error[jj];
      const real ti = (peak_y[jj] - Sx * invS) * sinv;
      St  += ti * sinv;
      Stt += ti * ti;
      sum += ti * sdens[jj] * sinv;
      avg += sdens[jj];
    }

  const real invStt = UNITY / Stt;
  const real aa = sum * invStt;
  const real bb = (Sy - aa * Sx) * invS;
  const real sa = SQRT(invStt);
  const real sb = SQRT(invS * (UNITY + Sx * invS * invStt * (Sx - TWO * St)));
  const real aa_peak = Sy_peak / S_peak;
  const real sa_peak = SQRT(Ss_peak) / S_peak;
  width /= (real)mem;
  avg /= (real)mem;

  real chisq_grad = ZERO;
  real chisq_peak = ZERO;
  for(int jj = ini; jj < fin; jj++)
    if( cont[jj] ){
      real diff = (aa * peak_y[jj] + bb - sdens[jj]) / error[jj];
      chisq_grad += diff * diff;
      diff = (aa_peak - peak_x[jj]) * sinv_peak;
      chisq_peak  += diff * diff;
    }
  const real rchisq_grad = (mem > 2) ? (chisq_grad / (real)(mem - 2)) : REAL_MAX;
  const real rchisq_peak = (mem > 1) ? (chisq_peak / (real)(mem - 1)) : REAL_MAX;

  *slope = aa;  *slope_sigma = sa;
  *icept = bb;  *icept_sigma = sb;
  *chisq = chisq_grad;
  *rchisq = rchisq_grad;
  *mean_count = avg;

  *axis = aa_peak;
  *axis_sigma = sa_peak;
  *chisq_axis = chisq_peak;
  *rchisq_axis = rchisq_peak;
  *mean_width = width;

  return (mem);
}



int main(int argc, char **argv)
{
  /** configure the details of the numerical simulation */
  /** read command line arguments */
  if( argc < 3 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 3);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -restart=<int> (optional)\n");
    __FPRINTF__(stderr, "          -ndim=<int>\n");
    __FPRINTF__(stderr, "          -lambda=<int> (optional) -mu=<int> (optional)\n");
    __FPRINTF__(stderr, "          -cm=<double> (optional)\n");
    __FPRINTF__(stderr, "          -sigma=<double> (optional)\n");
    __FPRINTF__(stderr, "          -init=<char *>\n");
    __FPRINTF__(stderr, "          -ftol=<double> (optional) -htol=<double> (optional) -xtol=<double> (optional)\n");
    __FPRINTF__(stderr, "          -maxEval=<int> (optional) -maxIter=<int> (optional)\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 3 ){ */
  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv, "file", &file));


  /** input locations of RC stars */
  const int Nstar = 10892;
  static real * xi;  xi  = (real *)malloc(sizeof(real) * Nstar);  if(  xi == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xi.");  }
  static real *eta;  eta = (real *)malloc(sizeof(real) * Nstar);  if( eta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate eta.");  }
  FILE *fp;
  char filename[128];
  sprintf(filename, "obs/RC.dat");
  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  bool checker = true;
  for(int ii = 0; ii < Nstar; ii++)
    checker &= (2 == fscanf(fp, "%le\t%le\t%*e\t%*e", &xi[ii], &eta[ii]));
  fclose(fp);
  if( !checker ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }


  /** input locations of mask regions */
  const int Nmask = 35;
  static real *  xi_mask;    xi_mask = (real *)malloc(sizeof(real) * Nmask);  if(   xi_mask == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xi_mask.");  }
  static real * eta_mask;   eta_mask = (real *)malloc(sizeof(real) * Nmask);  if(  eta_mask == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate eta_mask.");  }
  static real * rad_mask;   rad_mask = (real *)malloc(sizeof(real) * Nmask);  if(  rad_mask == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rad_mask.");  }
  sprintf(filename, "obs/mask.dat");
  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  checker = true;
  for(int ii = 0; ii < Nmask; ii++)
    checker &= (3 == fscanf(fp, "%le\t%le\t%le", &xi_mask[ii], &eta_mask[ii], &rad_mask[ii]));
  fclose(fp);
  if( !checker ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }
  static real *x_mask;  x_mask = (real *)malloc(sizeof(real) * Nmask);  if( x_mask == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate x_mask.");  }
  static real *y_mask;  y_mask = (real *)malloc(sizeof(real) * Nmask);  if( y_mask == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate y_mask.");  }


  /** set Subaru/HSC fields */
  const int Nfield = 5;
  static real * xi_field;   xi_field = (real *)malloc(sizeof(real) * Nfield);  if(  xi_field == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xi_field.");  }
  static real *eta_field;  eta_field = (real *)malloc(sizeof(real) * Nfield);  if( eta_field == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate eta_field.");  }
  static real *rad_field;  rad_field = (real *)malloc(sizeof(real) * Nfield);  if( rad_field == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate rad_field.");  }
  xi_field[0] = -4.653744156858361 ;  eta_field[0] = 3.635683362752261;  rad_field[0] = 0.75;/**< f003 */
  xi_field[1] = -5.7227710294896745;  eta_field[1] = 4.47091006650441 ;  rad_field[1] = 0.75;/**< f004 */
  xi_field[2] = -5.530826768687439 ;  eta_field[2] = 5.816784796173433;  rad_field[2] = 0.75;/**< f009 */
  xi_field[3] = -4.842739589150247 ;  eta_field[3] = 2.290423176281251;  rad_field[3] = 0.75;/**< f022 */
  xi_field[4] = -5.912437203950868 ;  eta_field[4] = 3.125050968157858;  rad_field[4] = 0.75;/**< f023 */
  for(int ii = 0; ii < Nfield; ii++)
    rad_field[ii] = (45.0 - 2.0) / 60.0;
  static real *x_field;  x_field = (real *)malloc(sizeof(real) * Nfield);  if( x_field == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate x_field.");  }
  static real *y_field;  y_field = (real *)malloc(sizeof(real) * Nfield);  if( y_field == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate y_field.");  }


  /** set bin width in the base grid (must be smaller than 3.651741e-2 in order-of-magnitude) */
  const real hh = 1.220703125e-4;/**< by using this resolution, all masked mesh do not contain RC star */
  /* const real hh = 6.103515625e-5; */
  const real hinv = UNITY / hh;


#ifdef  USE_HDF5_FORMAT
  /** create HDF5 file to dump results */
  sprintf(filename, "%s/%s_%s.h5", DATAFOLDER, file, RESULTS);
  static hid_t target;
  target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#ifdef  DOUBLE_PRECISION
  const hid_t H5T_NATIVE_REAL = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
  const hid_t H5T_NATIVE_REAL = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION
  hid_t H5T_NATIVE_STR = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(H5T_NATIVE_STR, 32));/**< memo: sizeof(char) is unity */

  hid_t input = H5Gcreate(target, "input", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t idim = 1;
  hid_t ispace = H5Screate_simple(1, &idim, NULL);
  hid_t attr = H5Acreate(input, "argc", H5T_NATIVE_INT, ispace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attr, H5T_NATIVE_INT, &argc));
  chkHDF5err(H5Aclose(attr));
  for(int ii = 0; ii < argc; ii++){
    sprintf(filename, "argv[%d]", ii);
    attr = H5Acreate(input, filename, H5T_NATIVE_STR, ispace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attr, H5T_NATIVE_STR, argv[ii]));
    chkHDF5err(H5Aclose(attr));
  }/* for(int ii = 0; ii < argc; ii++){ */
  chkHDF5err(H5Sclose(ispace));
  chkHDF5err(H5Gclose(input));
  chkHDF5err(H5Tclose(H5T_NATIVE_STR));
#endif//USE_HDF5_FORMAT


#pragma omp parallel
  {
    /** initialize pseudo random-number */
    struct cmaes_rand gauss;
    rand_state *rand;
    initRandNum(&rand);
#ifdef  USE_SFMTJUMP
    for(int ii = 0; ii < omp_get_thread_num(); ii++)
      SFMT_jump(rand, SFMTJUMP_10_100);
#endif//USE_SFMTJUMP

    /** configure the CMA-ES procedure */
    struct cmaes_config cfg;
    double *weight;
    config(argc, argv, &cfg, &weight);
#ifdef  NDIM
    if( cfg.ndim != NDIM ){
      __KILL__(stderr, "ERROR: dimension of prepared parameter-space is %d instead of %d which you specified.\n", NDIM, cfg.ndim);
    }/* if( ndim != NDIM ){ */
#endif//NDIM

    /** prepare the CMA-ES status */
    struct cmaes_status cma;
    double *mat, *offd, *tau, *diag, *orth;
    double *xx, *xmean, *xmold, *xbest, *xevol, *xtemp;
    double *ps, *pc;
    double *score, *history;
    struct dbl_int *sorted;
    allocate_cmaes(cfg, &cma, &mat, &offd, &tau, &diag, &orth, &xx, &xmean, &xmold, &xbest, &xevol, &xtemp, &ps, &pc, &score, &history, &sorted);

    double *xmean_init;
    xmean_init = (double *)malloc(sizeof(double) * cfg.ndim);
    if( xmean_init == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xmean_init\n");  }


    /** assume stream axis */
    static real degbin = 0.125;
    static real degmin = 60.0;
    static real degmax = 75.0;
    /* static real degmin = 50.0; */
    /* static real degmax = 90.0; */
    /* static real degmin = 55.0; */
    /* static real degmax = 80.0; */
    /* static int skipmax = 4096; */
    /* static int skipmin = 512; */
    /* static int skipmax = 8192; */
    /* static int skipmax = 4096; */
    /* static int skipmin = 128; */
    /* static int skipmax = 2048; */
    /* static int skipmin = 512; */
    static int skipmin = 128;
    static int skipmax = 256;
#pragma omp master
    {
      hsize_t attr_dims = 1;
      hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
      hid_t attribute = H5Acreate(target, "deg_min", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &degmin));
      chkHDF5err(H5Aclose(attribute));
      attribute = H5Acreate(target, "deg_max", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &degmax));
      chkHDF5err(H5Aclose(attribute));
      attribute = H5Acreate(target, "deg_bin", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &degbin));
      chkHDF5err(H5Aclose(attribute));
      attribute = H5Acreate(target, "skip_max", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skipmax));
      chkHDF5err(H5Aclose(attribute));
      attribute = H5Acreate(target, "skip_min", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skipmin));
      chkHDF5err(H5Aclose(attribute));
      chkHDF5err(H5Sclose(dataspace));
    }
#pragma omp barrier
    for(real deg = degmin; deg <= degmax; deg += degbin){
      const real theta = deg * M_PI / 180.0;
      const real cost = COS(0.5 * M_PI - theta);
      const real sint = SIN(0.5 * M_PI - theta);


#ifdef  USE_HDF5_FORMAT
      static hid_t group;
#pragma omp single
      {
	sprintf(filename, "deg_%.3f", deg);
	group = H5Gcreate(target, filename, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      }
#pragma omp master
      {
	hsize_t attr_dims = 1;
	hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
	hid_t attribute = H5Acreate(group, "angle", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
	chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &deg));
	chkHDF5err(H5Aclose(attribute));
	chkHDF5err(H5Sclose(dataspace));
      }
#endif//USE_HDF5_FORMAT


      /** set grid points */
      real xmin_loc = REAL_MAX;  real xmax_loc = -REAL_MAX;
      real ymin_loc = REAL_MAX;  real ymax_loc = -REAL_MAX;
#pragma omp for nowait
      for(int ii = 0; ii < Nfield; ii++){
	rotate(xi_field[ii], eta_field[ii], &x_field[ii], &y_field[ii], cost, sint);

	xmin_loc = FMIN(xmin_loc, x_field[ii] - rad_field[ii]);    xmax_loc = FMAX(xmax_loc, x_field[ii] + rad_field[ii]);
	ymin_loc = FMIN(ymin_loc, y_field[ii] - rad_field[ii]);    ymax_loc = FMAX(ymax_loc, y_field[ii] + rad_field[ii]);
      }/* for(int ii = 0; ii < Nfield; ii++){ */
      static real xmin, xmax, ymin, ymax;
#pragma omp single
      {
	xmin = REAL_MAX;	xmax = -REAL_MAX;
	ymin = REAL_MAX;	ymax = -REAL_MAX;
      }
#pragma omp critical
      xmin = FMIN(xmin, xmin_loc);
#pragma omp critical
      xmax = FMAX(xmax, xmax_loc);
#pragma omp critical
      ymin = FMIN(ymin, ymin_loc);
#pragma omp critical
      ymax = FMAX(ymax, ymax_loc);
#pragma omp barrier
      /* fprintf(stderr, "%e:%e; %e:%e\n", xmin, xmax, ymin, ymax); */
#pragma omp for
      for(int ii = 0; ii < Nmask; ii++)
	rotate(xi_mask[ii], eta_mask[ii], &x_mask[ii], &y_mask[ii], cost, sint);
      const ulong nx = (ulong)NEARBYINT(LDEXP(UNITY, (int)CEIL(LOG2((xmax - xmin) * hinv))));
      const ulong ny = (ulong)NEARBYINT(LDEXP(UNITY, (int)CEIL(LOG2((ymax - ymin) * hinv))));
#pragma omp barrier
#pragma omp single
      {
	xmin = HALF * (xmin + xmax) - hh * (real)(nx >> 1);  xmax = HALF * (xmin + xmax) + hh * (real)(nx >> 1);
	ymin = HALF * (ymin + ymax) - hh * (real)(ny >> 1);  ymax = HALF * (ymin + ymax) + hh * (real)(ny >> 1);
      }
      static real *xx;
#pragma omp single
      xx = (real *)malloc(sizeof(real) * nx);
      if( xx == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xx.");  }
      static real *yy;
#pragma omp single
      yy = (real *)malloc(sizeof(real) * ny);
      if( yy == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate yy.");  }
#pragma omp for nowait
      for(ulong ii = 0; ii < nx; ii++)    xx[ii] = xmin + hh * (HALF + (real)ii);/**< zone center */
#pragma omp for
      for(ulong jj = 0; jj < ny; jj++)    yy[jj] = ymin + hh * (HALF + (real)jj);/**< zone center */


      /** generate star count map */
      static real *star;
#pragma omp single
      star = (real *)malloc(sizeof(real) * nx * ny);
      if( star == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate star.");  }
      static real *area;
#pragma omp single
      area = (real *)malloc(sizeof(real) * nx * ny);
      if( area == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate area.");  }


      /** generate surface density map in the finnest grid */
#pragma omp for
      for(ulong ii = 0; ii < nx * ny; ii++){
	star[ii] = ZERO;
	area[ii] = hh * hh;
      }/* for(ulong ii = 0; ii < nx * ny; ii++){ */
#pragma omp for
      for(int ii = 0; ii < Nstar; ii++){
	real xx, yy;
	rotate(xi[ii], eta[ii], &xx, &yy, cost, sint);

	const int jj = (int)NEARBYINT((xx - xmin) * hinv);
	const int kk = (int)NEARBYINT((yy - ymin) * hinv);

	if( (jj >= 0) && (jj < (int)nx) && (kk >= 0) && (kk < (int)ny) )
#pragma omp critical
	  star[INDEX2D(nx, ny, (ulong)jj, (ulong)kk)] += UNITY;
      }/* for(int ii = 0; ii < Nstar; ii++){ */

      /** exclude external regions of FoV by Subaru/HSC or masked regions */
#pragma omp for
      for(ulong ii = 0; ii < nx; ii++)
	for(ulong jj = 0; jj < ny; jj++){
	  const ulong idx = INDEX2D(nx, ny, ii, jj);

	  bool interior = false;
	  /** detect grid points within observation fields */
	  for(int kk = 0; kk < Nfield; kk++){
	    const real dx = xx[ii] - x_field[kk];
	    const real dy = yy[jj] - y_field[kk];
	    if( (dx * dx + dy * dy) <= (rad_field[kk] * rad_field[kk]) ){
	      interior = true;
	      break;
	    }/* if( dx * dx + dy * dy <= rad_field[kk] * rad_field[kk] ){ */
	  }/* for(int kk = 0; kk < Nfield; kk++){ */

	  if( interior ){
	    /** exclude grid points within mask regions */
	    for(int kk = 0; kk < Nmask; kk++){
	      const real dx = xx[ii] - x_mask[kk];
	      const real dy = yy[jj] - y_mask[kk];
	      const real r2 = dx * dx + dy * dy;

	      if( r2 <= (rad_mask[kk] * rad_mask[kk]) ){
		interior = false;
		break;
	      }/* if( r2 < (rad_mask[kk] * rad_mask[kk]) ){ */
	    }/* for(int kk = 0; kk < Nmas; kk++){ */
	  }/* if( interior ){ */

	  if( !interior ){
	    star[idx] = ZERO;
	    area[idx] = ZERO;
	  }/* if( !interior ){ */
	}/* for(ulong jj = 0; jj < ny; jj++){ */


      for(int skip = skipmax; skip >= skipmin; skip >>= 1){
	/** generate subtracted surface density map */
	const int numx = nx / skip;
	const int numy = ny / skip;
	static real *dat;
#pragma omp single
	dat = (real *)malloc(sizeof(real) * numx * numy);
	if( dat == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate dat.");  }
	static real *eff;
#pragma omp single
	eff = (real *)malloc(sizeof(real) * numx * numy);
	if( eff == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate eff.");  }
	static real *hor;
#pragma omp single
	hor = (real *)malloc(sizeof(real) * (numx + 1));
	if( hor == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate hor.");  }
	static real *ver;
#pragma omp single
	ver = (real *)malloc(sizeof(real) * (numy + 1));
	if( ver == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate ver.");  }

	const real dx = (xmax - xmin) / (real)numx;
	const real dy = (ymax - ymin) / (real)numy;

#pragma omp for nowait
	for(int ii = 0; ii < numx + 1; ii++)      hor[ii] = xmin + dx * (real)ii;/**< node position */
#pragma omp for nowait
	for(int jj = 0; jj < numy + 1; jj++)      ver[jj] = ymin + dy * (real)jj;/**< node position */


#pragma omp for
	for(int ii = 0; ii < numx; ii++)
	  for(int jj = 0; jj < numy; jj++){
	    const int ll = INDEX2D(numx, numy, ii, jj);
	    dat[ll] = ZERO;
	    real val = ZERO;
	    real sub = ZERO;
	    for(int i1 = 0; i1 < skip; i1++){
	      const ulong iii = i1 + ii * skip;
	      for(int j1 = 0; j1 < skip; j1++){
		const ulong idx = INDEX2D(nx, ny, iii, (ulong)(j1 + jj * skip));
		val += star[idx];
		sub += area[idx];
	      }/* for(int j1 = 0; j1 < skip; j1++){ */
	    }/* for(int i1 = 0; i1 < skip; i1++){ */

	    /* mask dominated or edge regions should be removed from the analysis */
	    /* if( sub >= (QUARTER * dx * dy) ){ */
	    /* if( sub >= (HALF * dx * dy) ){ */
	    if( sub >= (THREE * QUARTER * dx * dy) ){
	      dat[ll] = val / (DBL_MIN + sub);
	      eff[ll] = sub;
	    }
	  }/* for(int jj = 0; jj < numy; jj++){ */

      	static real *peak_x;
#pragma omp single
	peak_x = (real *)malloc(sizeof(real) * numy);    if( peak_x == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate peak_x.");  }
      	static real *peak_y;
#pragma omp single
	peak_y = (real *)malloc(sizeof(real) * numy);    if( peak_y == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate peak_y.");  }
	static real *sdens;/**< surface density */
#pragma omp single
	sdens  = (real *)malloc(sizeof(real) * numy);    if( sdens  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate sdens.");  }
	static real *grad  ;
#pragma omp single
	grad   = (real *)malloc(sizeof(real) * numy);    if( grad   == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate grad.");  }
	static real *fwhm  ;
#pragma omp single
	fwhm   = (real *)malloc(sizeof(real) * numy);    if( fwhm   == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate fwhm.");  }
	static real *sigma ;
#pragma omp single
	sigma  = (real *)malloc(sizeof(real) * numy);    if( sigma  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate sigma.");  }
	static real *error ;
#pragma omp single
	error  = (real *)malloc(sizeof(real) * numy);    if( error  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate error.");  }
	static real *noise ;
#pragma omp single
	noise  = (real *)malloc(sizeof(real) * numy);    if( noise  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate noise.");  }
	static real *slope ;
#pragma omp single
	slope  = (real *)malloc(sizeof(real) * numy);    if( slope  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate slope.");  }
	static real *fchisq;
#pragma omp single
	fchisq = (real *)malloc(sizeof(real) * numy);    if( fchisq == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate fchisq.");  }
	static int *Ndat;
#pragma omp single
	Ndat = (int *)malloc(sizeof(int) * numy);    if( Ndat == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Ndat.");  }
	real *hist;
	hist = (real *)malloc(sizeof(real) * numx);	if( hist == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate hist.");	}
	real *isig2;
	isig2 = (real *)malloc(sizeof(real) * numx);	if( isig2 == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate isig2.");	}
	bool *cont;
	cont  = (bool *)malloc(sizeof(bool) * numy);    if( cont  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate cont.");  }


#ifdef  USE_HDF5_FORMAT
	static hid_t sub;
#pragma omp master
	{
	  sprintf(filename, "h_%.2u", ilog2(skip));
	  sub = H5Gcreate(group, filename, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  hsize_t attr_dims = 1;
	  hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
	  hid_t attribute = H5Acreate(sub, "skip", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
	  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &skip));
	  chkHDF5err(H5Aclose(attribute));
	  attribute = H5Acreate(sub, "bin", H5T_NATIVE_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
	  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_REAL, &dy));
	  chkHDF5err(H5Aclose(attribute));
	  attribute = H5Acreate(sub, "nx", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
	  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &numx));
	  chkHDF5err(H5Aclose(attribute));
	  attribute = H5Acreate(sub, "ny", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
	  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &numy));
	  chkHDF5err(H5Aclose(attribute));
	  chkHDF5err(H5Sclose(dataspace));
	}
#endif//USE_HDF5_FORMAT


#pragma omp for
	for(int jj = 0; jj < numy; jj++){
	  /** remove spikes due to mask regions */
	  hist[0] = ZERO;
	  real prev = dat[INDEX2D(numx, numy, 0, jj)];
	  real here = dat[INDEX2D(numx, numy, 1, jj)];
	  for(int ii = 1; ii < numx - 1; ii++){
	    real next = dat[INDEX2D(numx, numy, ii + 1, jj)];
	    hist[ii] = ((eff[INDEX2D(numx, numy, ii, jj)] > HALF * dx * dy) || (here < TWO * prev) || (here < TWO * next)) ? here : ZERO;
	    prev = here;
	    here = next;
	  }/* for(int ii = 1; ii < numx - 1; ii++){ */
	  hist[numx - 1] = ZERO;

	  real xl = xmin;
	  for(int ii = 0; ii < numx; ii++)
	    if( hist[ii] > ZERO ){
	      xl = hor[ii];
	      break;
	    }
	  real xr = xmax;
	  for(int ii = numx - 1; ii >= 0; ii--)
	    if( hist[ii] > ZERO ){
	      xr = hor[ii];
	      break;
	    }

	  /** initial guess of the solution */
	  /* double *xmean_init; */
	  /* number of parameters: (0) center of gaussian, (1) width of gaussian, (2) height of gaussian, (3) slope of background/foreground, and (4) amount of background/foreground; minize the chi-square of the profile fitting */
	  xmean_init[0] = HALF * (xl + xr);
	  /* xmean_init[1] = FMIN(0.2, 0.125 * (xr - xl));/\**< width is smaller than 0.3 degree in Kirihara et al. (2017) *\/ */
	  xmean_init[1] = 0.2;/**< width is smaller than 0.3 degree in Kirihara et al. (2017) */
	  xmean_init[3] = ZERO;
	  double stars = 0.0;
	  /* for(int ii = 0; ii < numx; ii++) */
	  /*   stars += dat[INDEX2D(numx, numy, ii, jj)]; */
	  for(int ii = 0; ii < numx; ii++)
	    stars += hist[ii];
	  xmean_init[2] = 1.0;
	  /* fprintf(stdout, "%e, %e\n", stars, stars * dx / (xr - xl)); */
	  xmean_init[4] = fmin(1.0, dx / (xr - xl));
	  /* for(int ii = 0; ii < numx; ii++) */
	  /*   hist[ii] = dat[INDEX2D(numx, numy, ii, jj)] / (DBL_MIN + stars); */
	  const double inv = 1.0 / (DBL_MIN + stars);
	  for(int ii = 0; ii < numx; ii++){
	    hist [ii] *= inv;
#if 0
	    isig2[ii]  = sqrt(dat[INDEX2D(numx, numy, ii, jj)] / (DBL_MIN + eff[INDEX2D(numx, numy, ii, jj)])) * inv;
#else
	    isig2[ii] = UNITY;
#endif
	  }

	  /** initialize the CMA-ES status */
	  initialize_cmaes(argc, argv, cfg, &cma, xmean_init);
	  /* if( jj == 0 ){ */
	  /*   fprintf(stdout, "# Configuration of CMA-ES (Covariance Matrix Adaptation Evolution Strategies)\n"); */
	  /*   fprintf(stdout, "#####\n"); */
	  /*   fprintf(stdout, "%d # cfg.ndim\n", cfg.ndim); */
	  /*   fprintf(stdout, "%d # cfg.lambda\n", cfg.lambda); */
	  /*   fprintf(stdout, "%d # cfg.mu\n", cfg.mu); */
	  /*   fprintf(stdout, "%e # cfg.mueff; %e indicates a reasonable setting of weight\n", cfg.mueff, 0.25 * cfg.lambda); */
	  /*   fprintf(stdout, "#####\n"); */
	  /*   fprintf(stdout, "%d # cfg.maxHistory\n", cfg.maxHistory); */
	  /*   fprintf(stdout, "#####\n"); */
	  /*   fprintf(stdout, "%e # cfg.cs; %e might be a better choice\n", cfg.cs, 4.0 / (double)cfg.ndim); */
	  /*   fprintf(stdout, "%e # cfg.ds\n", cfg.ds); */
	  /*   fprintf(stdout, "#####\n"); */
	  /*   fprintf(stdout, "%e # cfg.cm; c_m < 1 can be advantageous on noisy functions\n", cfg.cm); */
	  /*   fprintf(stdout, "%e # cfg.cc\n", cfg.cc); */
	  /*   fprintf(stdout, "%e # cfg.c1; %e might be a better choice\n", cfg.c1, 2.0 / (double)(cfg.ndim * cfg.ndim)); */
	  /*   fprintf(stdout, "%e # cfg.cmu; %e might be a better choice\n", cfg.cmu, fmin(cfg.mueff / (double)(cfg.ndim * cfg.ndim), 1.0 - cfg.c1)); */
	  /*   fprintf(stdout, "%e # 1 - c_1 - c_mu * sum w_j; close or equal to 0 is expected\n", cfg.coeff); */
	  /*   fprintf(stdout, "#####\n"); */
	  /*   fprintf(stdout, "%e # cfg.ftol\n", cfg.ftol); */
	  /*   fprintf(stdout, "%e # cfg.htol\n", cfg.htol); */
	  /*   fprintf(stdout, "%e # cfg.xtol\n", cfg.xtol); */
	  /*   fprintf(stdout, "%d # cfg.maxEval\n", cfg.maxEval); */
	  /*   fprintf(stdout, "%d # cfg.maxIter\n", cfg.maxIter); */
	  /*   fprintf(stdout, "#####\n"); */
	  /* }/\* if( jj == 0 ){ *\/ */

	  /** iteration until finding the solution (gaussian + a x + b) */
	  bool success = (stars > 0.0) ? optimize(&cma, cfg, &gauss, rand, numx, hist, isig2, hor) : false;

	  /* if( stars > 0.0 ){ */
	  /*   /\* optimize(&cma, cfg, &gauss, rand, numx, numy, jj, dat, hor); *\/ */
	  /*   optimize(&cma, cfg, &gauss, rand, numx, hist, hor); */
	  if( success ){
	    cma.xbest[2] *= stars;
	    cma.xbest[3] *= stars;
	    cma.xbest[4] *= stars;
	  }/* if( stars > 0.0 ){ */
	  else{
	    cma.xbest[0] = HALF * (xmin + xmax);
	    cma.xbest[1] = dx;
	    cma.xbest[2] = ZERO;
	    cma.xbest[3] = ZERO;
	    cma.xbest[4] = ZERO;
	  }/* else{ */


	  /** find the peak of the stream */
	  peak_x[jj] = cma.xbest[0];
	  peak_y[jj] = HALF * (ver[jj] + ver[jj + 1]);

	  /** evaluate FWHM of the stream */
	  sigma[jj] = cma.xbest[1];
	  fwhm[jj] = 2.0 * sqrt(2.0 * M_LN2) * cma.xbest[1];

	  /** evaluate surface density of the stream */
	  sdens[jj] =      cma.xbest[2];
	  /* error[jj] = sqrt(cma.xbest[2]); */
	  error[jj] = sqrt(sdens[jj] / (fwhm[jj] * dy));

	  /** record estimation of foreground/background */
	  slope[jj] = cma.xbest[3];
	  noise[jj] = cma.xbest[4];

	  /** record final score of CMA-ES procedure */
	  fchisq[jj] = cma.bestScore;
	  Ndat[jj] = 0;
	  for(int ii = 0; ii < numx; ii++)
	    Ndat[jj] += (hist[ii] > ZERO);
	}/* for(int jj = 0; jj < numy; jj++){ */


	/** additional analysis to estimate gradient of the stream */
	/** 1. estimate global gradient using the whole data */
	/** 2. remove data points which are off-axis */
	/** 3. remove Northern data points (up to 2 degrees from the Northern edge) */
	/** 4. estimate gradient of the stream with the launch point and the termination point */


	/** find (global) launch or termination point */
	static int gini, gfin;
#pragma omp single nowait
	for(int jj = numy - 1; jj >= 0; jj--)
	  if( sdens[jj] > ZERO ){
	    gfin = jj + 1;
	    break;
	  }/* if( sdens[jj] > ZERO ){ */
#pragma omp single
	for(int jj = 0; jj < numy; jj++)
	  if( sdens[jj] > ZERO ){
	    gini = jj;
	    break;
	  }/* if( sdens[jj] > ZERO ){ */

	static int gmem;
	static real gslope, gslope_sigma, gicept, gicept_sigma, gchisq, grchisq, gmean_count;
	static real gaxis, gaxis_sigma, gchisq_axis, grchisq_axis, gmean_width;
	static real gymin, gymax;
#pragma omp single
	{
	  gmem = 0;
	  gchisq = grchisq = gchisq_axis = grchisq_axis = REAL_MAX;
	  gslope = gslope_sigma = gicept = gicept_sigma = REAL_MAX;
	  gaxis = gaxis_sigma = REAL_MAX;
	  gmean_count = gmean_width = REAL_MAX;
	  gymin = gymax = REAL_MAX;
	}

	/** estimate number of trial fits */
	const int Nskip = (int)CEIL(TWO / (ver[1] - ver[0]));/**< up to 2 degrees from the Northern edge may be skipped */


#pragma omp for
	for(int remove = 0; remove < Nskip; remove++){
	  const int ini = gini;
	  const int fin = gfin - remove;
	  const real grad_ymin = HALF * (ver[ini] + ver[ini + 1]);
	  const real grad_ymax = HALF * (ver[fin] + ver[fin - 1]);

	  real slope, slope_sigma, icept, icept_sigma, chisq, rchisq, mean_count;
	  real axis, axis_sigma, chisq_axis, rchisq_axis, mean_width;

	  const int mem = least_square_method(ini, fin, cont, hor, peak_x, peak_y, sdens, sigma, fwhm, error,
					      &slope, &slope_sigma, &icept, &icept_sigma, &chisq, &rchisq, &mean_count,
					      &axis, &axis_sigma, &chisq_axis, &rchisq_axis, &mean_width);

#pragma omp critical
	  {
	    if( rchisq < grchisq ){
	      gslope       = slope;
	      gslope_sigma = slope_sigma;
	      gicept       = icept;
	      gicept_sigma = icept_sigma;
	      gchisq       = chisq;
	      grchisq      = rchisq;
	      gmean_count  = mean_count;
	      gaxis        = axis;
	      gaxis_sigma  = axis_sigma;
	      gchisq_axis  = chisq_axis;
	      grchisq_axis = rchisq_axis;
	      gmean_width  = mean_width;
	      gmem         = mem;

	      gymin = grad_ymin;
	      gymax = grad_ymax;
	    }
	  }

	}


#pragma omp master
	reportFittingResults(theta, gmem, gmean_width, gymin, gymax, gmean_count, gslope, gslope_sigma, gicept, gicept_sigma, gchisq, grchisq, gaxis, gaxis_sigma, gchisq_axis, grchisq_axis,
			     numx, hor, numy, ver, dat, peak_x, fwhm, sdens, sigma, slope, noise, fchisq, Ndat,
#ifdef  USE_HDF5_FORMAT
			     H5T_NATIVE_REAL, sub
#else///USE_HDF5_FORMAT
			     file, skip, peak_y, peak_y, grad, num
#endif//USE_HDF5_FORMAT
			     );
#pragma omp barrier


	free(hist);
	free(isig2);
	free(cont);

#pragma omp single
	{
	  free(peak_x);
	  free(peak_y);
	  free(sdens);
	  free(grad);
	  free(fwhm);
	  free(sigma);
	  free(error);
	  free(noise);
	  free(slope);
	  free(fchisq);
	  free(Ndat);

	  free(dat);
	  free(hor);
	  free(ver);
	}

#ifdef  USE_HDF5_FORMAT
#pragma omp master
	chkHDF5err(H5Gclose(sub));
#pragma omp barrier
#endif//USE_HDF5_FORMAT
      }/* for(int skip = skipmax; skip >= skipmin; skip >>= 1){ */


#pragma omp single
      {
	free(xx);
	free(yy);
	free(star);
	free(area);
      }

#ifdef  USE_HDF5_FORMAT
#pragma omp master
      chkHDF5err(H5Gclose(group));
#pragma omp barrier
#endif//USE_HDF5_FORMAT
    }/* for(real deg = degmin; deg <= degmax; deg += degbin){ */

    free(xmean_init);

    free(cma.score);  free(cma.sorted);  free(cma.history);
    free(cma.pc);  free(cma.ps);
    free(cma.xmean);  free(cma.xmold);  free(cma.xbest);  free(cma.xevol);  free(cma.xtemp);
    free(cma.xx);
    free(cma.diag);  free(cma.orth);
    free(cma.mat);  free(cma.offd);  free(cma.tau);
    free(cfg.weight);
  }

#ifdef  USE_HDF5_FORMAT
  chkHDF5err(H5Fclose(target));
#endif//USE_HDF5_FORMAT

  free(x_field);  free(y_field);
  free(xi_field);  free(eta_field);  free(rad_field);
  free(xi_mask);  free(eta_mask);  free(rad_mask);
  free(x_mask);  free(y_mask);
  free(xi);  free(eta);


  return (0);
}
