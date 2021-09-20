/**
 * @file cmaes_io.c
 *
 * @brief Source code for Input/Output functions in CMA-ES
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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifndef DISABLE_MPI
#include <mpi.h>
#include "mpilib.h"
#endif//DISABLE_MPI

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "name.h"
#include "rand.h"

#include "../find/cmaes.h"
#include "../find/cmaes_io.h"


#ifdef  USE_HDF5_FORMAT
/**
 * @fn createHDF5type
 *
 * @brief Create data type for HDF5.
 *
 * @return (sort) data type in HDF5 for struct dbl_int
 * @return (sfmt) data type in HDF5 for w128_t
 */
void createHDF5type
(hid_t *sort
#ifdef  USE_SFMT
 , hid_t *sfmt
#endif//USE_SFMT
 )
{
  __NOTE__("%s\n", "start");

  *sort = H5Tcreate(H5T_COMPOUND, sizeof(struct dbl_int));
  chkHDF5err(H5Tinsert(*sort, "val", HOFFSET(struct dbl_int, val), H5T_NATIVE_DOUBLE));
  chkHDF5err(H5Tinsert(*sort, "int", HOFFSET(struct dbl_int, idx), H5T_NATIVE_INT));

#ifdef  USE_SFMT
/*   #include <emmintrin.h> */
/* /\** 128-bit data structure *\/ */
/* union W128_T { */
/*     uint32_t u[4]; */
/*     uint64_t u64[2]; */
/*     __m128i si; */
/* }; */
/* typedef union W128_T w128_t; */

  *sfmt = H5Tcreate(H5T_COMPOUND, sizeof(w128_t));
  chkHDF5err(H5Tinsert(*sfmt, "u64[0]", HOFFSET(w128_t, u64[0]), H5T_NATIVE_ULONG));
  chkHDF5err(H5Tinsert(*sfmt, "u64[1]", HOFFSET(w128_t, u64[1]), H5T_NATIVE_ULONG));
#endif//USE_SFMT

  __NOTE__("%s\n", "end");
}
/**
 * @fn removeHDF5type
 *
 * @brief Remove data type for HDF5.
 *
 * @param (sort) data type in HDF5 for struct dbl_int
 * @param (sfmt) data type in HDF5 for w128_t
 */
void removeHDF5type
(hid_t  sort
#ifdef  USE_SFMT
 , hid_t  sfmt
#endif//USE_SFMT
 )
{
  __NOTE__("%s\n", "start");

  chkHDF5err(H5Tclose(sort));
#ifdef  USE_SFMT
  chkHDF5err(H5Tclose(sfmt));
#endif//USE_SFMT

  __NOTE__("%s\n", "end");
}
#endif//USE_HDF5_FORMAT


/**
 * @fn saveConfig
 *
 * @brief Write the configuration of the CMA-ES procedure.
 *
 * @param (file) name of the simulation
 * @param (cfg) configuration of the CMA-ES
 */
void saveConfig(char *file, const struct cmaes_config cfg)
{
  __NOTE__("%s\n", "start");

  char filename[128];
  sprintf(filename, "%s/%s_%s.%s", DATAFOLDER, file, CONFIG,
#ifdef  USE_HDF5_FORMAT
	  "h5"
#else///USE_HDF5_FORMAT
	  "dat"
#endif//USE_HDF5_FORMAT
	  );

#ifdef  USE_HDF5_FORMAT
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t group = H5Gcreate(target, CONFIG, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /** write attribute data */
  /** create the data space for the attribute */
  hsize_t attr_dims = 1;
  hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  attribute = H5Acreate(group, "mueff", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.mueff));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "cm", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.cm));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "cs", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.cs));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ds", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.ds));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "cc", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.cc));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "c1", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.c1));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "cmu", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.cmu));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "coeff", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.coeff));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "chiN", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.chiN));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ftol", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.ftol));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "htol", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.htol));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "xtol", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cfg.xtol));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "maxEval", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cfg.maxEval));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "maxIter", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cfg.maxIter));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "ndim", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cfg.ndim));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "lambda", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cfg.lambda));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "mu", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cfg.mu));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "maxHistory", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cfg.maxHistory));
  chkHDF5err(H5Aclose(attribute));
  /** close the dataspace */
  chkHDF5err(H5Sclose(dataspace));

  hsize_t dims = cfg.lambda;
  dataspace = H5Screate_simple(1, &dims, NULL);
  hid_t dataset;
  /** write weight */
  dataset = H5Dcreate(group, "weight", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cfg.weight));
  chkHDF5err(H5Dclose(dataset));
  /** close the dataspace */
  chkHDF5err(H5Sclose(dataspace));

  chkHDF5err(H5Gclose(group));
  chkHDF5err(H5Fclose(target));
#else///USE_HDF5_FORMAT
  FILE *fp;
  fp = fopen(filename, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  bool success = true;
  size_t tmp;
  tmp = 1;  if( tmp != fwrite(&cfg.mueff, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.cm   , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.cs   , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.ds   , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.cc   , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.c1   , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.cmu  , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.coeff, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.chiN , sizeof(double), tmp, fp) )    success = false;

  tmp = 1;  if( tmp != fwrite(&cfg.ftol, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.htol, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.xtol, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.maxEval, sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.maxIter, sizeof(int), tmp, fp) )    success = false;

  tmp = 1;  if( tmp != fwrite(&cfg.ndim  , sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.lambda, sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fwrite(&cfg.mu    , sizeof(int), tmp, fp) )    success = false;

  tmp = 1;  if( tmp != fwrite(&cfg.maxHistory, sizeof(int), tmp, fp) )    success = false;

  tmp = cfg.lambda;  if( tmp != fwrite(cfg.weight, sizeof(double), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }
  fclose(fp);
#endif//USE_HDF5_FORMAT

  __NOTE__("%s\n", "end");
}


/**
 * @fn loadConfig
 *
 * @brief Read the configuration of the CMA-ES procedure.
 *
 * @param (file) name of the simulation
 * @return (cfg) configuration of the CMA-ES
 * @return (array) data array to store weight
 */
void loadConfig(char *file, struct cmaes_config *cfg, double **array)
{
  __NOTE__("%s\n", "start");

  char filename[128];
  sprintf(filename, "%s/%s_%s.%s", DATAFOLDER, file, CONFIG,
#ifdef  USE_HDF5_FORMAT
	  "h5"
#else///USE_HDF5_FORMAT
	  "dat"
#endif//USE_HDF5_FORMAT
	  );

#ifdef  USE_HDF5_FORMAT
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t group = H5Gopen(target, CONFIG, H5P_DEFAULT);

  /** read attribute data */
  hid_t attribute;
  attribute = H5Aopen(group, "mueff", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->mueff));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "cm", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->cm));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "cs", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->cs));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "ds", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->ds));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "cc", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->cc));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "c1", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->c1));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "cmu", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->cmu));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "coeff", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->coeff));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "chiN", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->chiN));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "ftol", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->ftol));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "htol", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->htol));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "xtol", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cfg->xtol));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "maxEval", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &cfg->maxEval));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "maxIter", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &cfg->maxIter));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "ndim", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &cfg->ndim));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "lambda", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &cfg->lambda));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "mu", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &cfg->mu));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "maxHistory", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &cfg->maxHistory));
  chkHDF5err(H5Aclose(attribute));
#else///USE_HDF5_FORMAT
  FILE *fp;
  fp = fopen(filename, "rb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  bool success = true;
  size_t tmp;
  tmp = 1;  if( tmp != fread(&cfg->mueff, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->cm   , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->cs   , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->ds   , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->cc   , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->c1   , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->cmu  , sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->coeff, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->chiN , sizeof(double), tmp, fp) )    success = false;

  tmp = 1;  if( tmp != fread(&cfg->ftol, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->htol, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->xtol, sizeof(double), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->maxEval, sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->maxIter, sizeof(int), tmp, fp) )    success = false;

  tmp = 1;  if( tmp != fread(&cfg->ndim  , sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->lambda, sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&cfg->mu    , sizeof(int), tmp, fp) )    success = false;

  tmp = 1;  if( tmp != fread(&cfg->maxHistory, sizeof(int), tmp, fp) )    success = false;
#endif//USE_HDF5_FORMAT

  *array = (double *)malloc(sizeof(double) * cfg->lambda);
  if( *array == NULL ){    __KILL__(stderr, "ERROR: failure to allocate array\n");  }
  cfg->weight = *array;

#ifdef  USE_HDF5_FORMAT
  /** read weight */
  hid_t dataset;
  dataset = H5Dopen(group, "weight", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cfg->weight));
  chkHDF5err(H5Dclose(dataset));

  chkHDF5err(H5Gclose(group));
  chkHDF5err(H5Fclose(target));
#else///USE_HDF5_FORMAT
  tmp = cfg->lambda;  if( tmp != fread(cfg->weight, sizeof(double), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }
  fclose(fp);
#endif//USE_HDF5_FORMAT

  __NOTE__("%s\n", "end");
}


/**
 * @fn saveStatus
 *
 * @brief Write the current status of the CMA-ES procedure.
 *
 * @param (file) name of the simulation
 * @param (cfg) configuration of the CMA-ES
 * @param (cma) status of the CMA-ES
 * @param (gauss) status of the pseudo random-number in Gaussian distribution
 * @param (rand) status of the pseudo random-number
 * @param (H5T_NATIVE_DBLINT) data type of struct dbl_int for HDF5
 * @param (H5T_NATIVE_SFMT) data type of w128_t for HDF5
 */
void saveStatus
(char *file, const struct cmaes_config cfg, const struct cmaes_status cma, const struct cmaes_rand gauss, const rand_state rand
#ifdef  USE_HDF5_FORMAT
 , const hid_t H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
 , const hid_t H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
 )
{
  __NOTE__("%s\n", "start");

  char filename[128];
  sprintf(filename, "%s/%s_%s_gen%.3d.%s", DATAFOLDER, file, STATUS, cma.gen,
#ifdef  USE_HDF5_FORMAT
	  "h5"
#else///USE_HDF5_FORMAT
	  "dat"
#endif//USE_HDF5_FORMAT
	  );

  const int ndim = cfg.ndim;
  const int lambda = cfg.lambda;
  const int nhist = cfg.maxHistory;

#ifdef  USE_HDF5_FORMAT
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t group = H5Gcreate(target, STATUS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /** write attribute data */
  /** create the data space for the attribute */
  hsize_t attr_dims = 1;
  hid_t dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  attribute = H5Acreate(group, "bestScore", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cma.bestScore));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "funcRange", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cma.funcRange));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "histRange", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cma.histRange));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "sigma", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cma.sigma));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "gen", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cma.gen));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "state", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cma.state));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "count", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cma.count));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "eigen_uptodate", H5T_NATIVE_HBOOL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_HBOOL, &cma.eigen_uptodate));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(group, "initial_phase", H5T_NATIVE_HBOOL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_HBOOL, &cma.initial_phase));
  chkHDF5err(H5Aclose(attribute));
  /** close the dataspace */
  chkHDF5err(H5Sclose(dataspace));

  hid_t dataset;
  hsize_t dims = ndim * ndim;
  dataspace = H5Screate_simple(1, &dims, NULL);
  dataset = H5Dcreate(group, "mat", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.mat));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "orth", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.orth));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(dataspace));

  dims = lambda * ndim;
  dataspace = H5Screate_simple(1, &dims, NULL);
  dataset = H5Dcreate(group, "xx", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.xx));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(dataspace));

  dims = lambda;
  dataspace = H5Screate_simple(1, &dims, NULL);
  dataset = H5Dcreate(group, "sorted", H5T_NATIVE_DBLINT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DBLINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.sorted));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "score", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.score));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(dataspace));

  dims = nhist;
  dataspace = H5Screate_simple(1, &dims, NULL);
  dataset = H5Dcreate(group, "history", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.history));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(dataspace));

  dims = ndim;
  dataspace = H5Screate_simple(1, &dims, NULL);
  dataset = H5Dcreate(group, "diag", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.diag));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "xmean", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.xmean));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "xmold", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.xmold));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "xbest", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.xbest));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "xevol", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.xevol));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "xtemp", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.xtemp));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "ps", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.ps));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "pc", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.pc));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(dataspace));

  dims = ndim - 1;
  dataspace = H5Screate_simple(1, &dims, NULL);
  dataset = H5Dcreate(group, "offd", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.offd));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "tau", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma.tau));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(dataspace));

  chkHDF5err(H5Gclose(group));

  group = H5Gcreate(target, RANDOM, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dims = SFMT_N;
  dataspace = H5Screate_simple(1, &dims, NULL);
  dataset = H5Dcreate(group, "sfmt_state", H5T_NATIVE_SFMT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_SFMT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rand.state));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(dataspace));
  dims = 1;
  dataspace = H5Screate_simple(1, &dims, NULL);
  dataset = H5Dcreate(group, "sfmt_idx", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rand.idx));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "rand_hold", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &gauss.hold));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dcreate(group, "rand_save", H5T_NATIVE_HBOOL, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_HBOOL, H5S_ALL, H5S_ALL, H5P_DEFAULT, &gauss.stored));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Gclose(group));

  chkHDF5err(H5Fclose(target));
#else///USE_HDF5_FORMAT
  FILE *fp;
  fp = fopen(filename, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  bool success = true;
  size_t tmp;
  tmp = ndim * ndim;  if( tmp != fwrite(cma. mat, sizeof(double), tmp, fp) )    success = false;
  tmp = ndim * ndim;  if( tmp != fwrite(cma.orth, sizeof(double), tmp, fp) )    success = false;
  tmp = lambda * ndim;  if( tmp != fwrite(cma.xx, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fwrite(cma. diag, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fwrite(cma.xmean, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fwrite(cma.xmold, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fwrite(cma.xbest, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fwrite(cma.xevol, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fwrite(cma.xtemp, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fwrite(cma.   ps, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fwrite(cma.   pc, sizeof(double), tmp, fp) )    success = false;
  tmp = lambda;  if( tmp != fwrite(cma.score, sizeof(double), tmp, fp) )    success = false;
  tmp = lambda;  if( tmp != fwrite(cma.sorted, sizeof(struct dbl_int), tmp, fp) )    success = false;
  tmp = ndim-1;  if( tmp != fwrite(cma. offd, sizeof(double), tmp, fp) )    success = false;
  tmp = ndim-1;  if( tmp != fwrite(cma.  tau, sizeof(double), tmp, fp) )    success = false;
  tmp =  nhist;  if( tmp != fwrite(cma.history, sizeof(double), tmp, fp) )    success = false;
  tmp = SFMT_N;  if( tmp != fwrite( rand.state, sizeof(w128_t), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&rand.idx  , sizeof(int)   , tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&gauss.hold, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&gauss.stored, sizeof(bool), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.bestScore, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.funcRange, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.histRange, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.    sigma, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.      gen, sizeof(int), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.    state, sizeof(int), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.    count, sizeof(int), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.eigen_uptodate, sizeof(bool), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.initial_phase , sizeof(bool), tmp, fp) )    success = false;
  /* tmp = SFMT_N;  if( tmp != fread( rand->state, sizeof(w128_t), tmp, fp) )    success = false; */
  /* tmp =      1;  if( tmp != fread(&rand->idx  , sizeof(int)   , tmp, fp) )    success = false; */
  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }
  fclose(fp);
#endif//USE_HDF5_FORMAT


  sprintf(filename, "%s/%s_%s.%s", DATAFOLDER, file, STATUS,
#ifdef  USE_HDF5_FORMAT
	  "h5"
#else///USE_HDF5_FORMAT
	  "dat"
#endif//USE_HDF5_FORMAT
	  );
#ifdef  USE_HDF5_FORMAT
  target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /** write attribute data */
  /** create the data space for the attribute */
  attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  attribute = H5Acreate(target, "bestScore", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cma.bestScore));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "funcRange", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cma.funcRange));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "histRange", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cma.histRange));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "sigma", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &cma.sigma));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "gen", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cma.gen));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "state", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cma.state));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "count", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &cma.count));
  chkHDF5err(H5Aclose(attribute));
  /** close the dataspace */
  chkHDF5err(H5Sclose(dataspace));
  chkHDF5err(H5Fclose(target));
#else///USE_HDF5_FORMAT
  fp = fopen(filename, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  success = true;
  tmp =      1;  if( tmp != fwrite(&cma.bestScore, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.funcRange, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.histRange, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.    sigma, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.      gen, sizeof(int), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.    state, sizeof(int), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fwrite(&cma.    count, sizeof(int), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);  }
  fclose(fp);
#endif//USE_HDF5_FORMAT


  __NOTE__("%s\n", "end");
}


/**
 * @fn loadStatus
 *
 * @brief Read the current status of the CMA-ES procedure.
 *
 * @param (file) name of the simulation
 * @param (cfg) configuration of the CMA-ES
 * @param (gen) generation of CMA-ES status
 * @return (cma) status of the CMA-ES
 * @return (gauss) status of the pseudo random-number in Gaussian distribution
 * @return (rand) status of the pseudo random-number
 * @param (H5T_NATIVE_DBLINT) data type of struct dbl_int for HDF5
 * @param (H5T_NATIVE_SFMT) data type of w128_t for HDF5
 */
void loadStatus
(char *file, struct cmaes_status *cma, struct cmaes_rand *gauss, rand_state *rand, const int gen
#ifdef  USE_HDF5_FORMAT
 , const hid_t H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
 , const hid_t H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
 )
{
  __NOTE__("%s\n", "start");

  char filename[128];
  sprintf(filename, "%s/%s_%s_gen%.3d.%s", DATAFOLDER, file, STATUS, gen,
#ifdef  USE_HDF5_FORMAT
	  "h5"
#else///USE_HDF5_FORMAT
	  "dat"
#endif//USE_HDF5_FORMAT
	  );

#ifdef  USE_HDF5_FORMAT
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t group = H5Gopen(target, STATUS, H5P_DEFAULT);

  /** read attribute data */
  hid_t attribute;
  attribute = H5Aopen(group, "bestScore", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cma->bestScore));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "funcRange", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cma->funcRange));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "histRange", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cma->histRange));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "sigma", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &cma->sigma));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "gen", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &cma->gen));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "state", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &cma->state));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "count", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &cma->count));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "eigen_uptodate", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_HBOOL, &cma->eigen_uptodate));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Aopen(group, "initial_phase", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_HBOOL, &cma->initial_phase));
  chkHDF5err(H5Aclose(attribute));

  hid_t dataset;
  dataset = H5Dopen(group, "mat", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->mat));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "orth", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->orth));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "xx", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->xx));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "diag", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->diag));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "xmean", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->xmean));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "xmold", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->xmold));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "xbest", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->xbest));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "xevol", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->xevol));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "xtemp", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->xtemp));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "ps", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->ps));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "pc", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->pc));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "score", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->score));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "sorted", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DBLINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->sorted));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "offd", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->offd));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "tau", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->tau));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "history", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cma->history));
  chkHDF5err(H5Dclose(dataset));

  chkHDF5err(H5Gclose(group));
  group = H5Gopen(target, RANDOM, H5P_DEFAULT);
  dataset = H5Dopen(group, "sfmt_state", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_SFMT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rand->state));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "sfmt_idx", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rand->idx));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "rand_hold", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &gauss->hold));
  chkHDF5err(H5Dclose(dataset));
  dataset = H5Dopen(group, "rand_save", H5P_DEFAULT);
  chkHDF5err(H5Dread(dataset, H5T_NATIVE_HBOOL, H5S_ALL, H5S_ALL, H5P_DEFAULT, &gauss->stored));
  chkHDF5err(H5Dclose(dataset));
  chkHDF5err(H5Gclose(group));
  chkHDF5err(H5Fclose(target));
#else///USE_HDF5_FORMAT
  FILE *fp;
  fp = fopen(filename, "wb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }

  bool success = true;
  size_t tmp;
  tmp = ndim * ndim;  if( tmp != fread(cma-> mat, sizeof(double), tmp, fp) )    success = false;
  tmp = ndim * ndim;  if( tmp != fread(cma->orth, sizeof(double), tmp, fp) )    success = false;
  tmp = lambda * ndim;  if( tmp != fread(cma->xx, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fread(cma-> diag, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fread(cma->xmean, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fread(cma->xmold, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fread(cma->xbest, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fread(cma->xevol, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fread(cma->xtemp, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fread(cma->   ps, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fread(cma->   pc, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fread(cma->score, sizeof(double), tmp, fp) )    success = false;
  tmp =   ndim;  if( tmp != fread(cma->sorted, sizeof(struct dbl_int), tmp, fp) )    success = false;
  tmp = ndim-1;  if( tmp != fread(cma-> offd, sizeof(double), tmp, fp) )    success = false;
  tmp = ndim-1;  if( tmp != fread(cma->  tau, sizeof(double), tmp, fp) )    success = false;
  tmp =  nhist;  if( tmp != fread(cma->history, sizeof(double), tmp, fp) )    success = false;
  tmp = SFMT_N;  if( tmp != fread( rand->state, sizeof(w128_t), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&rand->idx  , sizeof(int)   , tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&gauss->hold, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&gauss->stored, sizeof(bool), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&cma->bestScore, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&cma->funcRange, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&cma->histRange, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&cma->    sigma, sizeof(double), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&cma->      gen, sizeof(int), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&cma->    state, sizeof(int), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&cma->    count, sizeof(int), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&cma->eigen_uptodate, sizeof(bool), tmp, fp) )    success = false;
  tmp =      1;  if( tmp != fread(&cma->initial_phase , sizeof(bool), tmp, fp) )    success = false;
  if( success != true ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);  }
  fclose(fp);
#endif//USE_HDF5_FORMAT

  __NOTE__("%s\n", "end");
}
