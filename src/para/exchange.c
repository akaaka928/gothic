/**
 * @file exchange.c
 *
 * @brief Source code for domain decomposition using MPI
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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

#include "macro.h"
#include "mpilib.h"

#include "../misc/structure.h"

#include "mpicfg.h"
#include "exchange.h"


/**
 * @fn allocateORMtopology
 *
 * @brief Memory allocation for Orthogonal Recursive Multi-section with setting MPI topology.
 */
muse allocateORMtopology(float **dxmin, float **dxmax, float **dymin, float **dymax, float **dzmin, float **dzmax, MPI_Request **dmreq,
			 float **sxmin, float **sxmax, float **symin, float **symax, float **szmin, float **szmax,
			 sendCfg **sendBuf, recvCfg **recvBuf, int **rnum, int **disp,
			 MPIinfo orm[], MPIinfo rep[],
			 const int Nx, const int Ny, const int Nz, MPIcfg_tree *mpi,
			 domainCfg *dom, sampling *sample, const ulong Ntot)
{
  __NOTE__("%s\n", "start");


  muse alloc = {0, 0};

  mpi->dim[0] = (Nx > 0) ? Nx : 1;  mpi->prd[0] = false;
  mpi->dim[1] = (Ny > 0) ? Ny : 1;  mpi->prd[1] = false;
  mpi->dim[2] = (Nz > 0) ? Nz : 1;  mpi->prd[2] = false;

  chkMPIerr(MPI_Cart_create(mpi->comm, 3, mpi->dim, mpi->prd, false, &(mpi->cart)));
  chkMPIerr(MPI_Cart_coords(mpi->cart, mpi->rank, 3, mpi->pos));


  /** MPI process decomposition */
  *dxmin = (float *)malloc(mpi->size * sizeof(float));  alloc.host += mpi->size * sizeof(float);  if( *dxmin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dxmin\n");  }
  *dxmax = (float *)malloc(mpi->size * sizeof(float));  alloc.host += mpi->size * sizeof(float);  if( *dxmax == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dxmax\n");  }
  *dymin = (float *)malloc(mpi->size * sizeof(float));  alloc.host += mpi->size * sizeof(float);  if( *dymin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dymin\n");  }
  *dymax = (float *)malloc(mpi->size * sizeof(float));  alloc.host += mpi->size * sizeof(float);  if( *dymax == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dymax\n");  }
  *dzmin = (float *)malloc(mpi->size * sizeof(float));  alloc.host += mpi->size * sizeof(float);  if( *dzmin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dzmin\n");  }
  *dzmax = (float *)malloc(mpi->size * sizeof(float));  alloc.host += mpi->size * sizeof(float);  if( *dzmax == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dzmax\n");  }
  dom->xmin = *dxmin;  dom->ymin = *dymin;  dom->zmin = *dzmin;
  dom->xmax = *dxmax;  dom->ymax = *dymax;  dom->zmax = *dzmax;
  *dmreq = (MPI_Request *)malloc(mpi->size * sizeof(MPI_Request));  alloc.host += mpi->size * sizeof(MPI_Request);
  if( *dmreq == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dmreq\n");  }
  dom->req = *dmreq;

  MPIinfo org;
  org.comm = mpi->cart;
  org.size = mpi->size;
  org.rank = mpi->rank;
  for(int ii = 0; ii < 3; ii++){
    splitMPI(org.comm, mpi->pos[ii], org.rank, &(orm[ii]));
    splitMPI(org.comm, orm[ii].rank, org.rank, &(rep[ii]));
    org = orm[ii];
  }/* for(int ii = 0; ii < 3; ii++){ */


  /** memory allocation for ORM related arrays */
  *sendBuf = (sendCfg *)malloc(mpi->size * sizeof(sendCfg));  alloc.host += mpi->size * sizeof(sendCfg);
  *recvBuf = (recvCfg *)malloc(mpi->size * sizeof(recvCfg));  alloc.host += mpi->size * sizeof(recvCfg);
  if( *sendBuf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sendBuf\n");  }
  if( *recvBuf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate recvBuf\n");  }


  /** memory allocation for arrays to set domain boundaries */
  *sxmin = (float *)malloc(mpi->dim[0] * sizeof(float));  alloc.host += mpi->dim[0] * sizeof(float);
  *sxmax = (float *)malloc(mpi->dim[0] * sizeof(float));  alloc.host += mpi->dim[0] * sizeof(float);
  *symin = (float *)malloc(mpi->dim[1] * sizeof(float));  alloc.host += mpi->dim[1] * sizeof(float);
  *symax = (float *)malloc(mpi->dim[1] * sizeof(float));  alloc.host += mpi->dim[1] * sizeof(float);
  *szmin = (float *)malloc(mpi->dim[2] * sizeof(float));  alloc.host += mpi->dim[2] * sizeof(float);
  *szmax = (float *)malloc(mpi->dim[2] * sizeof(float));  alloc.host += mpi->dim[2] * sizeof(float);
  if( *sxmin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sxmin\n");  }
  if( *sxmax == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sxmax\n");  }
  if( *symin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate symin\n");  }
  if( *symax == NULL ){    __KILL__(stderr, "ERROR: failure to allocate symax\n");  }
  if( *szmin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate szmin\n");  }
  if( *szmax == NULL ){    __KILL__(stderr, "ERROR: failure to allocate szmax\n");  }

  sample->xmin = *sxmin;  sample->ymin = *symin;  sample->zmin = *szmin;
  sample->xmax = *sxmax;  sample->ymax = *symax;  sample->zmax = *szmax;


  /** set sampling rate and allocate corresponding size of temporary array */
  sample->rate = (((float)Ntot * DEFAULT_SAMPLING_RATE) >= ((float)mpi->size * DEFAULT_SAMPLING_NUMBER)) ?
    (DEFAULT_SAMPLING_RATE) : (DEFAULT_SAMPLING_NUMBER * (float)mpi->size / (float)Ntot);


  bool root = (mpi->rank != 0) ? false : true;
  for(int ii = 0; ii < 3; ii++)
    if( (orm[ii].rank == 0) && (orm[ii].size > 1) )
      root = true;

  if( root ){
    /** recvcnts, displs */
    *rnum = (int *)malloc(mpi->size * sizeof(int));    alloc.host += mpi->size * sizeof(int);    if( *rnum == NULL ){      __KILL__(stderr, "ERROR: failure to allocate rnum\n");    }
    *disp = (int *)malloc(mpi->size * sizeof(int));    alloc.host += mpi->size * sizeof(int);    if( *disp == NULL ){      __KILL__(stderr, "ERROR: failure to allocate disp\n");    }
    sample->rnum = *rnum;
    sample->disp = *disp;
  }/* if( root ){ */

  sample->Nmax = (int)ceilf((float)Ntot * MAX_FACTOR_FROM_EQUIPARTITION * sample->rate);


  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn releaseORMtopology
 *
 * @brief Memory deallocation for Orthogonal Recursive Multi-section.
 */
void  releaseORMtopology(float  *dxmin, float  *dxmax, float  *dymin, float  *dymax, float  *dzmin, float  *dzmax, MPI_Request  *dmreq,
			 float  *sxmin, float  *sxmax, float  *symin, float  *symax, float  *szmin, float  *szmax,
			 sendCfg  *sendBuf, recvCfg  *recvBuf, int  *rnum, int  *disp,
			 MPIinfo orm[], MPIinfo rep[], const int rank)
{
  __NOTE__("%s\n", "start");


  for(int ii = 2; ii >= 0; ii--){
    freeMPIgroup(&(rep[ii]));
    freeMPIgroup(&(orm[ii]));
  }/* for(int ii = 2; ii >= 0; ii--){ */

  free(dxmin);  free(dymin);  free(dzmin);
  free(dxmax);  free(dymax);  free(dzmax);  free(dmreq);

  free(sxmin);	free(symin);  free(szmin);
  free(sxmax);	free(symax);  free(szmax);

  free(sendBuf);
  free(recvBuf);

  if( rank == 0 ){
    free(rnum);
    free(disp);
  }/* if( rank == 0 ){ */


  __NOTE__("%s\n", "end");
}
