/**
 * @file mpicfg.h
 *
 * @brief Header file for MPI parallelization
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/01 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef MPICFG_H
#define MPICFG_H


#include <mpi.h>

#include "macro.h"
#include "mpilib.h"


/**
 * @def GPUS_PER_PROCESS
 *
 * @brief Number of GPUs per MPI process
 */
#ifndef GPUS_PER_PROCESS
#define GPUS_PER_PROCESS (1)
#endif//GPUS_PER_PROCESS


/**
 * @def MAX_FACTOR_INCREASE
 *
 * @brief A parameter to control maximum deviation from equidistribution
 */
#define MAX_FACTOR_INCREASE (1.3f)
/* #define MAX_FACTOR_INCREASE (1.5f) */

/**
 * @def MAX_FACTOR_SAFETY
 *
 * @brief A safety parameter to control maximum deviation from equidistribution
 */
#define MAX_FACTOR_SAFETY (1.05f)
/* #define MAX_FACTOR_SAFETY (1.1f) */

/**
 * @def MAX_FACTOR_FROM_EQUIPARTITION
 *
 * @brief Maximum deviation from equidistribution (Ntot / Ngpus)
 */
#define MAX_FACTOR_FROM_EQUIPARTITION (MAX_FACTOR_INCREASE * MAX_FACTOR_SAFETY)


/**
 * @struct MPIcfg_tree
 *
 * @brief structure for exchanging tree structures
 */
typedef struct
{
  MPI_Datatype ipos, jpos, more, mass;
#ifdef  GADGET_MAC
  MPI_Datatype iacc;
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
  MPI_Datatype ivel, time;
#endif//BLOCK_TIME_STEP
  MPI_Comm comm, cart;
  int rank, size;
  int dim[3], prd[3], pos[3];
} MPIcfg_tree;


/* list of functions appeared in ``mpicfg.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  void setNodeConfig(ulong Ntot, int *Nnode, int *Ni, MPIinfo mpi, MPIcfg_tree *let, const int devID);
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//MPICFG_H
