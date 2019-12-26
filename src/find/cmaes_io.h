/**
 * @file cmaes_io.h
 *
 * @brief Header file for Input/Output functions in CMA-ES
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
#ifndef CMAES_IO_H
#define CMAES_IO_H


#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#endif//USE_HDF5_FORMAT

#include "../find/cmaes.h"


/** macros to specify the file name */
#define CONFIG "cfg"
#define STATUS "state"
#define RANDOM "rand"



/** list of functions appeared in ``cmaes_io.c'' */
#ifdef  USE_HDF5_FORMAT
void createHDF5type
(hid_t *sort
#ifdef  USE_SFMT
 , hid_t *sfmt
#endif//USE_SFMT
 );
void removeHDF5type
(hid_t  sort
#ifdef  USE_SFMT
 , hid_t  sfmt
#endif//USE_SFMT
 );
#endif//USE_HDF5_FORMAT

void saveConfig(char *file, const struct cmaes_config  cfg);
void loadConfig(char *file,       struct cmaes_config *cfg, double **array);

void saveStatus
(char *file, const struct cmaes_config cfg, const struct cmaes_status cma, const struct cmaes_rand gauss, const rand_state rand
#ifdef  USE_HDF5_FORMAT
 , const hid_t H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
 , const hid_t H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
 );
void loadStatus
(char *file, struct cmaes_status *cma, struct cmaes_rand *gauss, rand_state *rand, const int gen
#ifdef  USE_HDF5_FORMAT
 , const hid_t H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
 , const hid_t H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
 );


#endif//CMAES_IO_H
