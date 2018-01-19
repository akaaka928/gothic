/**
 * @file exchange.h
 *
 * @brief Header file for domain decomposition using MPI
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/19 (Fri)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef EXCHANGE_H
#define EXCHANGE_H


#include <mpi.h>

#include "macro.h"
#include "mpilib.h"

#include "../misc/structure.h"
#include "../misc/tune.h"
#include "../misc/brent.h"

#include "../para/mpicfg.h"


/**
 * @def DEFAULT_SAMPLING_RATE
 *
 * @brief Constant to set particle sampling rate.
 */
#define DEFAULT_SAMPLING_RATE   (1.0e-4f)

/**
 * @def DEFAULT_SAMPLING_NUMBER
 *
 * @brief Constant to set particle sampling rate.
 */
/* #define DEFAULT_SAMPLING_NUMBER (64.0f) */
#define DEFAULT_SAMPLING_NUMBER (128.0f)
/* #define DEFAULT_SAMPLING_NUMBER (256.0f) */
/* #define DEFAULT_SAMPLING_NUMBER (512.0f) */


/**
 * @struct sampling
 *
 * @brief structure for particle sampling
 */
typedef struct
{
  int *rnum, *disp;/**< former recvNum and recvDsp */
  float *xmin, *xmax, *ymin, *ymax, *zmin, *zmax;
  float rate;/**< former samplingRate */
  int Nmax;/**< former sampleNumMax */
} sampling;

/**
 * @struct samplePos
 *
 * @brief structure for particle sampling
 */
typedef struct
{
  float *x_hst, *y_hst, *z_hst;
  float *x_dev, *y_dev, *z_dev;
  int *i_hst, *i_dev;
#ifdef  MPI_ONE_SIDED_FOR_EXCG
  MPI_Win win_x, win_y, win_z;
#endif//MPI_ONE_SIDED_FOR_EXCG
} samplePos;

/**
 * @struct particlePos
 *
 * @brief structure for particle position
 */
typedef struct
{
  float *x, *y, *z;
} particlePos;


/**
 * @struct domainCfg
 *
 * @brief structure for domain decomposition
 */
typedef struct
{
  float *xmin, *xmax, *ymin, *ymax, *zmin, *zmax;
  MPI_Request *req;
} domainCfg;

/**
 * @struct sendCfg
 *
 * @brief structure for domain decomposition (send)
 */
typedef struct
{
#ifdef  MPI_ONE_SIDED_FOR_EXCG
  sendBody body;
  MPI_Request req;
#else///MPI_ONE_SIDED_FOR_EXCG
  int num, head;
  MPI_Request pos, idx;
#ifdef  GADGET_MAC
  MPI_Request acc;
#endif//GADGET_MAC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  MPI_Request ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  MPI_Request vel, time;
#else///BLOCK_TIME_STEP
  MPI_Request vx, vy, vz;
#endif//BLOCK_TIME_STEP
#endif//MPI_ONE_SIDED_FOR_EXCG
  /* move them to sendDom domBoundary */
  /* real xmin, xmax, ymin, ymax, zmin, zmax; */
  int rank;
} sendCfg;

/**
 * @struct recvCfg
 *
 * @brief structure for domain decomposition (receive)
 */
typedef struct
{
#ifdef  MPI_ONE_SIDED_FOR_EXCG
  sendBody body;
  MPI_Request dsp;
#else///MPI_ONE_SIDED_FOR_EXCG
  int num, head;
#endif//MPI_ONE_SIDED_FOR_EXCG
  MPI_Request req;
#ifndef MPI_ONE_SIDED_FOR_EXCG
  MPI_Request pos, idx;
#ifdef  GADGET_MAC
  MPI_Request acc;
#endif//GADGET_MAC
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  MPI_Request ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  MPI_Request vel, time;
#else///BLOCK_TIME_STEP
  MPI_Request vx, vy, vz;
#endif//BLOCK_TIME_STEP
#endif//MPI_ONE_SIDED_FOR_EXCG
} recvCfg;

/**
 * @struct domainDecomposeKey
 *
 * @brief structure for domain decomposition (sorting)
 */
typedef struct
{
  int *dstRank_hst, *dstRank_dev, *bodyIdx_dev;
} domainDecomposeKey;


/* static const double loadImbalanceCrit = 0.9;/\**< constant to detect load imbalance *\/ */
static const double loadImbalanceCrit = 0.95;/**< constant to detect load imbalance */
#define SLOW_DOWN_PROCS_CRIT(tot) ((tot) >> 4)

/**
 * @struct loadImbalanceDetector
 *
 * @brief structure for detecting load imbalance
 */
typedef struct
{
  double tmin, tmax;
  bool enable, execute;
} loadImbalanceDetector;


/* list of functions appeared in ``exchange.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  muse allocateORMtopology(float **dxmin, float **dxmax, float **dymin, float **dymax, float **dzmin, float **dzmax, MPI_Request **dmreq,
			   float **sxmin, float **sxmax, float **symin, float **symax, float **szmin, float **szmax,
			   sendCfg **sendBuf, recvCfg **recvBuf, int **rnum, int **disp,
			   MPIinfo orm[], MPIinfo rep[],
			   const int Nx, const int Ny, const int Nz, MPIcfg_tree *mpi,
			   domainCfg *dom, sampling *sample, const ulong Ntot);
  void  releaseORMtopology(float  *dxmin, float  *dxmax, float  *dymin, float  *dymax, float  *dzmin, float  *dzmax, MPI_Request  *dmreq,
			   float  *sxmin, float  *sxmax, float  *symin, float  *symax, float  *szmin, float  *szmax,
			   sendCfg  *sendBuf, recvCfg  *recvBuf, int  *rnum, int  *disp,
			   MPIinfo orm[], MPIinfo rep[], const int rank);
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//EXCHANGE_H
