/*************************************************************************\
 *                                                                       *
                  last updated on 2016/08/12(Fri) 11:37:41
 *                                                                       *
 *    Implementations related to OpenMP/MPI hybrid parallelization       *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <name.h>
#include <mpilib.h>
//-------------------------------------------------------------------------
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "../tree/make.h"
//-------------------------------------------------------------------------
#include "../time/adv_dev.h"
//-------------------------------------------------------------------------
#include "mpicfg.h"
#include "exchange.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* function to set MPI topology for Orthogonal Recursive Bisection */
//-------------------------------------------------------------------------
void configORBtopology
(int *ndim, MPIinfo **orb, domainCfg **box, sendCfg **sendBuf, recvCfg **recvBuf, MPIcfg_tree *mpi,
 const ulong Ntot, real *samplingRate, int *sampleNumMax, real **sampleLoc, real **sampleFul,
 int **recvNum, int **recvDsp, real **boxMin, real **boxMax)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create an MPI topology */
  /* # of MPI processes is assumed to be power of two */
  //-----------------------------------------------------------------------
  const int log2procs = (int)ilog2((uint)mpi->size);
  for(int ii = 0; ii < 3; ii++)
    mpi->dim[ii] = 0;
  *ndim = 3;
  if( *ndim > log2procs )
    *ndim = log2procs;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < *ndim; ii++){
    mpi->dim[ii] = 1 << ((log2procs + ((*ndim) - (ii + 1))) / (*ndim));
    mpi->prd[ii] = false;
  }
  //-----------------------------------------------------------------------
  chkMPIerr(MPI_Cart_create(mpi->comm, *ndim, mpi->dim, mpi->prd, false, &(mpi->cart)));
  chkMPIerr(MPI_Cart_coords(mpi->cart, mpi->rank, *ndim, mpi->pos));
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* MPI process decomposition */
  //-----------------------------------------------------------------------
  *orb = (MPIinfo   *)malloc((*    ndim) * sizeof(MPIinfo));
  *box = (domainCfg *)malloc((mpi->size) * sizeof(domainCfg));
  //-----------------------------------------------------------------------
  if( *orb == NULL ){    __KILL__(stderr, "ERROR: failure to allocate orb");  }
  if( *box == NULL ){    __KILL__(stderr, "ERROR: failure to allocate box");  }
  //-----------------------------------------------------------------------
  MPIinfo old;
  old.comm = mpi->cart;
  old.size = mpi->size;
  old.rank = mpi->rank;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < *ndim; ii++){
    //---------------------------------------------------------------------
    splitMPI(old.comm, mpi->pos[ii], old.rank, &((*orb)[ii]));
    old = (*orb)[ii];
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation for ORB related arrays */
  //-----------------------------------------------------------------------
  *sendBuf = (sendCfg *)malloc((mpi->size) * sizeof(sendCfg));
  *recvBuf = (recvCfg *)malloc((mpi->size) * sizeof(recvCfg));
  //-----------------------------------------------------------------------
  if( *sendBuf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sendBuf");  }
  if( *recvBuf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate recvBuf");  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set sampling rate and allocate corresponding size of temporary array */
  //-----------------------------------------------------------------------
#if 1
  *samplingRate = DEFAULT_SAMPLING_RATE;
  if( (float)Ntot * (*samplingRate) / (float)mpi->size < DEFAULT_SAMPLING_NUMBER )
    *samplingRate = DEFAULT_SAMPLING_NUMBER * (float)mpi->size / (float)Ntot;
#else
  *samplingRate = 1.0e-5f;
  if( Ntot / (ulong)mpi->size <=    2048 )  *samplingRate = 1.0e-1f;
  if( Ntot / (ulong)mpi->size <=   16384 )  *samplingRate = 1.0e-2f;
  if( Ntot / (ulong)mpi->size <=  131072 )  *samplingRate = 1.0e-3f;
  if( Ntot / (ulong)mpi->size <= 1048576 )  *samplingRate = 1.0e-4f;
#endif
  //-----------------------------------------------------------------------
  bool root = false;
  if( mpi->rank == 0 )
    root = true;
  for(int ii = 0; ii < *ndim; ii++)
    if( (*orb)[ii].rank == 0 )
      root = true;
  //-----------------------------------------------------------------------
  if( root ){
    //---------------------------------------------------------------------
    *sampleNumMax = (int)ceilf((float)Ntot * (*samplingRate));
    *sampleFul = (real *)malloc((*sampleNumMax) * sizeof(real));
    if( *sampleFul == NULL ){      __KILL__(stderr, "ERROR: failure to allocate sampleFul");    }
    //---------------------------------------------------------------------
    /* recvcnts, displs */
    *recvNum = (int *)malloc(mpi->size * sizeof(int));
    *recvDsp = (int *)malloc(mpi->size * sizeof(int));
    if( *recvNum == NULL ){      __KILL__(stderr, "ERROR: failure to allocate recvNum");    }
    if( *recvDsp == NULL ){      __KILL__(stderr, "ERROR: failure to allocate recvDsp");    }
    //---------------------------------------------------------------------
  }/* if( root ){ */
  //-----------------------------------------------------------------------
  *sampleNumMax = (int)ceilf((float)NUM_BODY_MAX * MAX_FACTOR_FROM_EQUIPARTITION * (*samplingRate));
  *sampleLoc = (real *)malloc((*sampleNumMax) * 3 * sizeof(real));
  if( *sampleLoc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sampleLoc");  }
  //-----------------------------------------------------------------------
  /* minimum/maximun position of the decomposed domain */
  int maxDim = (mpi->dim[0] > mpi->dim[1]) ? mpi->dim[0] : mpi->dim[1];
  if( maxDim < mpi->dim[2] )      maxDim = mpi->dim[2];
  *boxMin = (real *)malloc(maxDim * sizeof(real));
  *boxMax = (real *)malloc(maxDim * sizeof(real));
  if( *boxMin == NULL ){      __KILL__(stderr, "ERROR: failure to allocate boxMin");    }
  if( *boxMax == NULL ){      __KILL__(stderr, "ERROR: failure to allocate boxMax");    }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void removeORBtopology
(const int ndim, MPIinfo  *orb, domainCfg  *box, sendCfg  *sendBuf, recvCfg  *recvBuf,
 real  *sampleLoc, real  *sampleFul, int  *recvNum, int  *recvDsp, real  *boxMin, real  *boxMax)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndim; ii++)
    freeMPIgroup(&(orb[ii]));
  //-----------------------------------------------------------------------
  free(orb);
  free(box);
  //-----------------------------------------------------------------------
  free(sendBuf);
  free(recvBuf);
  //-----------------------------------------------------------------------
  if( sampleFul != NULL ){
    free(sampleFul);
    free(recvNum);
    free(recvDsp);
  }
  //-----------------------------------------------------------------------
  free(sampleLoc);
  free(boxMin);
  free(boxMax);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void calcLocalBoxSize(const int num, iparticle body, real min[], real max[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#pragma unroll
  for(int ii = 0; ii < 3; ii++){
    min[ii] =  REAL_MAX;
    max[ii] = -REAL_MAX;
  }/* for(int ii = 0; ii < 3; ii++){ */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++){
    if( body.pos[ii].x < min[0] )      min[0] = body.pos[ii].x;
    if( body.pos[ii].y < min[1] )      min[1] = body.pos[ii].y;
    if( body.pos[ii].z < min[2] )      min[2] = body.pos[ii].z;
    if( body.pos[ii].x > max[0] )      max[0] = body.pos[ii].x;
    if( body.pos[ii].y > max[1] )      max[1] = body.pos[ii].y;
    if( body.pos[ii].z > max[2] )      max[2] = body.pos[ii].z;
  }/* for(int ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef __ICC
/* Disable ICC's remark #161: unrecognized #pragma */
#     pragma warning (disable:161)
#endif//__ICC
//-------------------------------------------------------------------------
int posAscendingOrder(const void *a, const void *b);
int posAscendingOrder(const void *a, const void *b)
{
  //-----------------------------------------------------------------------
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
  if(          (*(real *)a) > (*(real *)b) ){    return ( 1);  }
  else{    if( (*(real *)a) < (*(real *)b) ){    return (-1);  }
    else                                         return ( 0);  }
#pragma GCC diagnostic pop
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
int ddkeyAscendingOrder(const void *a, const void *b);
int ddkeyAscendingOrder(const void *a, const void *b)
{
  //-----------------------------------------------------------------------
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
  if(          ((domainDecomposeKey *)a)->dstRank > ((domainDecomposeKey *)b)->dstRank ){    return ( 1);  }
  else{    if( ((domainDecomposeKey *)a)->dstRank < ((domainDecomposeKey *)b)->dstRank ){    return (-1);  }
    else                                                         return ( 0);  }
#pragma GCC diagnostic pop
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifdef __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void exchangeParticles
(const int numOld,              iparticle src,
#   if  defined(GENERATE_PHKEY_ON_DEVICE) || !defined(CALC_MULTIPOLE_ON_DEVICE)
 iparticle ibody_dev,
#endif//defined(GENERATE_PHKEY_ON_DEVICE) || !defined(CALC_MULTIPOLE_ON_DEVICE)
 const int numMax, int *numNew, iparticle dst,
 domainDecomposeKey * restrict key,
 sendCfg *sendBuf, recvCfg *recvBuf,
 const float samplingRate, const int sampleNumMax, const double tloc, real *sampleLoc, real *sampleFul, int *recvNum, int *recvDsp, real *boxMin, real *boxMax,
 const int ndim, MPIinfo *orb, domainCfg *domain, MPIcfg_tree mpi
#ifdef  EXEC_BENCHMARK
 , wall_clock_time *execTime
#endif//EXEC_BENCHMARK
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* copy N-body particles from device to host */
  //-----------------------------------------------------------------------
#   if  defined(GENERATE_PHKEY_ON_DEVICE) || !defined(CALC_MULTIPOLE_ON_DEVICE)
  copyParticle_dev2hst(numOld, ibody_dev, src
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );
#endif//defined(GENERATE_PHKEY_ON_DEVICE) || !defined(CALC_MULTIPOLE_ON_DEVICE)
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* get (current) box size for the local distribution of N-body particles */
  //-----------------------------------------------------------------------
  real min[3], max[3];
  calcLocalBoxSize(numOld, src, min, max);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* pick up sample particles */
  //-----------------------------------------------------------------------
  /* weight is determined using elapsed time by each process */
  __NOTE__("rank %d: tloc = %e, numOld = %d, xmin = %e, xmax = %e\n", mpi.rank, tloc, numOld, min[0], max[0]);
#if 1
  if( tloc == 0.0 ){
    __KILL__(stderr, "ERROR: tloc is %e @ rank %d\n", tloc, mpi.rank);
  }/* if( tloc == 0.0 ){ */
#endif
  double ttot;
  chkMPIerr(MPI_Allreduce(&tloc, &ttot, 1, MPI_DOUBLE, MPI_SUM, mpi.comm));
  const int Nsub = (int)ceilf((float)numOld * samplingRate * (float)(tloc / ttot));
  __NOTE__("rank %d: tloc = %e, ttot = %e, frac = %e, numOld = %d, Nsub = %d\n", mpi.rank, tloc, ttot, tloc / ttot, numOld, Nsub);
  const int iskip = numOld / Nsub;
  int sampleNum = 0;
  for(int ii = 0; ii < numOld; ii += iskip){
    sampleLoc[sampleNum                   ] = src.pos[ii].x;
    sampleLoc[sampleNum +     sampleNumMax] = src.pos[ii].y;
    sampleLoc[sampleNum + 2 * sampleNumMax] = src.pos[ii].z;
    sampleNum++;
  }
#if 1
  sampleLoc[(sampleNum - 1)                   ] = src.pos[numOld - 1].x;
  sampleLoc[(sampleNum - 1) +     sampleNumMax] = src.pos[numOld - 1].y;
  sampleLoc[(sampleNum - 1) + 2 * sampleNumMax] = src.pos[numOld - 1].z;
#endif
  __NOTE__("rank %d: iskip = %d, sampleNum = %d\n", mpi.rank, iskip, sampleNum);
  //-----------------------------------------------------------------------
  /* sort sample particles in each direction */
#if 1
  qsort(  &sampleLoc[0               ], sampleNum, sizeof(real), posAscendingOrder);
  if( ndim >= 2 )
    qsort(&sampleLoc[    sampleNumMax], sampleNum, sizeof(real), posAscendingOrder);
  if( ndim == 3 )
    qsort(&sampleLoc[2 * sampleNumMax], sampleNum, sizeof(real), posAscendingOrder);
#endif
  //-----------------------------------------------------------------------
  __NOTE__("rank %d\n", mpi.rank);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* determine local domain */
  domainCfg local;
#pragma unroll
  for(int ii = 0; ii < 3; ii++){
    local.min[ii] = HALF * REAL_MIN;
    local.max[ii] = HALF * REAL_MAX;
  }/* for(int ii = 0; ii < 3; ii++){ */
  //-----------------------------------------------------------------------
  MPIinfo root;
  root.rank = mpi.rank;
  root.size = mpi.size;
  root.comm = mpi.comm;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndim; ii++){
    //---------------------------------------------------------------------
    /* gather sampling points to the root process */
    //---------------------------------------------------------------------
    /* set receive counts */
    const int nprocs = root.size;
    chkMPIerr(MPI_Gather(&sampleNum, 1, MPI_INT, recvNum, 1, MPI_INT, 0, root.comm));
    //---------------------------------------------------------------------
    /* set receive displacements */
    recvDsp[0] = 0;
    for(int jj = 1; jj < nprocs; jj++)
      recvDsp[jj] = recvDsp[jj - 1] + recvNum[jj - 1];
    const int recvNumTot = recvDsp[nprocs - 1] + recvNum[nprocs - 1];
    //---------------------------------------------------------------------
    /* gather particle data to the root process */
    chkMPIerr(MPI_Gatherv(&sampleLoc[ii * sampleNumMax], sampleNum, MPI_REALDAT, sampleFul, recvNum, recvDsp, MPI_REALDAT, 0, root.comm));
    //---------------------------------------------------------------------
    /* sort sample particles in each direction */
    if( root.rank == 0 )
      qsort(sampleFul, recvNumTot, sizeof(real), posAscendingOrder);
#if 0
    if( root.rank == 0 ){
      fprintf(stdout, "recvNumTot = %d\n", recvNumTot);
      for(int jj = 0; jj < recvNumTot; jj++)
	fprintf(stdout, "pos[%d] = %e\n", jj, sampleFul[jj]);
      fflush(stdout);
    }
#endif
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* set decomposed domain */
    //---------------------------------------------------------------------
    /* the root process determine the partition */
    if( root.rank == 0 ){
      //-------------------------------------------------------------------
      boxMin[0] = HALF * REAL_MIN;
      for(int jj = 0; jj < mpi.dim[ii] - 1; jj++){
	//-----------------------------------------------------------------
	const int idx = (recvNumTot * (1 + jj)) / mpi.dim[ii];
	const real middle = HALF * (sampleFul[idx] + sampleFul[idx + 1]);
	//-----------------------------------------------------------------
	boxMax[jj    ] = middle;
	boxMin[jj + 1] = middle;
      }/* for(int jj = 0; jj < mpi.dim[ii] - 1; jj++){ */
      boxMax[mpi.dim[ii] - 1] = HALF * REAL_MAX;
      //-------------------------------------------------------------------
    }/* if( root.rank == 0 ){ */
    //---------------------------------------------------------------------
    /* broadcast the partition */
    chkMPIerr(MPI_Bcast(boxMin, mpi.dim[ii], MPI_REALDAT, 0, root.comm));
    chkMPIerr(MPI_Bcast(boxMax, mpi.dim[ii], MPI_REALDAT, 0, root.comm));
    //---------------------------------------------------------------------
    /* set the local domain */
    local.min[ii] = boxMin[mpi.pos[ii]];
    local.max[ii] = boxMax[mpi.pos[ii]];
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* preparation to the next dimension */
    //---------------------------------------------------------------------
    root = orb[ii];
    //---------------------------------------------------------------------
    __NOTE__("rank %d\n", mpi.rank);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndim; ii++){ */
  //-----------------------------------------------------------------------
  /* share the decomposed domain */
  chkMPIerr(MPI_Allgather(&local, 1 * sizeof(domainCfg), MPI_BYTE, domain, 1 * sizeof(domainCfg), MPI_BYTE, mpi.comm));
  //-----------------------------------------------------------------------
  __NOTE__("rank %d: [%e, %e]x[%e, %e]x[%e, %e]\n", mpi.rank,
	   domain[mpi.rank].min[0], domain[mpi.rank].max[0],
	   domain[mpi.rank].min[1], domain[mpi.rank].max[1],
	   domain[mpi.rank].min[2], domain[mpi.rank].max[2]);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* exchange N-body particles */
  const int numProcs = mpi.size;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < numProcs; ii++)
    chkMPIerr(MPI_Irecv(&(recvBuf[ii].num), 1, MPI_INT, ii, ii, mpi.comm, &(recvBuf[ii].req)));
  //-----------------------------------------------------------------------
  const int nullSend = 0;
  int overlapNum = 0;
  for(int ii = 0; ii < numProcs; ii++){
    /* if( (min[0] < EPSILON + domain[ii].max[0]) && (EPSILON + max[0] > domain[ii].min[0]) && */
    /* 	(min[1] < EPSILON + domain[ii].max[1]) && (EPSILON + max[1] > domain[ii].min[1]) && */
    /* 	(min[2] < EPSILON + domain[ii].max[2]) && (EPSILON + max[2] > domain[ii].min[2]) ){ */
    if( (min[0] <= domain[ii].max[0]) && (max[0] >= domain[ii].min[0]) &&
	(min[1] <= domain[ii].max[1]) && (max[1] >= domain[ii].min[1]) &&
	(min[2] <= domain[ii].max[2]) && (max[2] >= domain[ii].min[2]) ){
      //-------------------------------------------------------------------
      /* spatial overlap is detected */
      //-------------------------------------------------------------------
      sendBuf[overlapNum].rank = ii;
      sendBuf[overlapNum].num  =  0;
      //-------------------------------------------------------------------
      /* sendBuf[overlapNum].xmin = (min[0] > domain[ii].min[0]) ? (min[0]) : (domain[ii].min[0]); */
      /* sendBuf[overlapNum].ymin = (min[1] > domain[ii].min[1]) ? (min[1]) : (domain[ii].min[1]); */
      /* sendBuf[overlapNum].zmin = (min[2] > domain[ii].min[2]) ? (min[2]) : (domain[ii].min[2]); */
      sendBuf[overlapNum].xmin =
	(domain[ii].min[0] < QUARTER * REAL_MIN) ? (domain[ii].min[0]) : ((min[0] > domain[ii].min[0]) ? (min[0]) : (domain[ii].min[0]));
      sendBuf[overlapNum].ymin =
	(domain[ii].min[1] < QUARTER * REAL_MIN) ? (domain[ii].min[1]) : ((min[1] > domain[ii].min[1]) ? (min[1]) : (domain[ii].min[1]));
      sendBuf[overlapNum].zmin =
	(domain[ii].min[2] < QUARTER * REAL_MIN) ? (domain[ii].min[2]) : ((min[2] > domain[ii].min[2]) ? (min[2]) : (domain[ii].min[2]));
      //-------------------------------------------------------------------
      /* sendBuf[overlapNum].xmax = (max[0] < domain[ii].max[0]) ? (max[0]) : (domain[ii].max[0]); */
      /* sendBuf[overlapNum].ymax = (max[1] < domain[ii].max[1]) ? (max[1]) : (domain[ii].max[1]); */
      /* sendBuf[overlapNum].zmax = (max[2] < domain[ii].max[2]) ? (max[2]) : (domain[ii].max[2]); */
      sendBuf[overlapNum].xmax =
	(domain[ii].max[0] > QUARTER * REAL_MAX) ? (domain[ii].max[0]) : ((max[0] < domain[ii].max[0]) ? (max[0]) : (domain[ii].max[0]));
      sendBuf[overlapNum].ymax =
	(domain[ii].max[1] > QUARTER * REAL_MAX) ? (domain[ii].max[1]) : ((max[1] < domain[ii].max[1]) ? (max[1]) : (domain[ii].max[1]));
      sendBuf[overlapNum].zmax =
	(domain[ii].max[2] > QUARTER * REAL_MAX) ? (domain[ii].max[2]) : ((max[2] < domain[ii].max[2]) ? (max[2]) : (domain[ii].max[2]));
      //-------------------------------------------------------------------
      overlapNum++;
      //-------------------------------------------------------------------
    }
    else{
      //-------------------------------------------------------------------
      /* covered areas do not overlap */
      //-------------------------------------------------------------------
      chkMPIerr(MPI_Isend(&nullSend, 1, MPI_INT, ii, mpi.rank, mpi.comm, &(domain[ii].req)));
      //-------------------------------------------------------------------
    }/* else{ */
    //---------------------------------------------------------------------
    __NOTE__("rank %d, ii = %d\n", mpi.rank, ii);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < numProcs; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* determine process rank for each particle to belong */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < numOld; ii++){
    //---------------------------------------------------------------------
    key[ii].bodyIdx = ii;
#ifndef NDEBUG
    bool find = false;
#endif//NDEBUG
    //---------------------------------------------------------------------
    for(int jj = 0; jj < overlapNum; jj++){
      //-------------------------------------------------------------------
      /* if( ((EPSILON + src.pos[ii].x) > sendBuf[jj].xmin) && (src.pos[ii].x < (EPSILON + sendBuf[jj].xmax)) && */
      /* 	  ((EPSILON + src.pos[ii].y) > sendBuf[jj].ymin) && (src.pos[ii].y < (EPSILON + sendBuf[jj].ymax)) && */
      /* 	  ((EPSILON + src.pos[ii].z) > sendBuf[jj].zmin) && (src.pos[ii].z < (EPSILON + sendBuf[jj].zmax)) ){ */
      if( (src.pos[ii].x >= sendBuf[jj].xmin) && (src.pos[ii].x <= sendBuf[jj].xmax) &&
	  (src.pos[ii].y >= sendBuf[jj].ymin) && (src.pos[ii].y <= sendBuf[jj].ymax) &&
	  (src.pos[ii].z >= sendBuf[jj].zmin) && (src.pos[ii].z <= sendBuf[jj].zmax) ){
	//-----------------------------------------------------------------
	key[ii].dstRank = sendBuf[jj].rank;
	sendBuf[jj].num++;
#ifndef NDEBUG
	find = true;
#endif//NDEBUG
	break;
	//-----------------------------------------------------------------
      }
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < overlapNum; jj++){ */
    //---------------------------------------------------------------------
#ifndef NDEBUG
    if( !find ){
#if 1
      fprintf(stderr, "numOld = %d, numProcs = %d, overlapNum = %d @ rank %d\n", numOld, numProcs, overlapNum, mpi.rank);
      fprintf(stderr, "local: xmin = %f, xmax = %f, ymin = %f, ymax = %f, zmin = %f, zmax = %f @ rank %d\n",
	      min[0], max[0], min[1], max[1], min[2], max[2], mpi.rank);
      for(int jj = 0; jj < numProcs; jj++){
	fprintf(stderr, "domain[%d]: xmin = %f, xmax = %f, ymin = %f, ymax = %f, zmin = %f, zmax = %f @ rank %d\n",
		jj, domain[jj].min[0], domain[jj].max[0], domain[jj].min[1], domain[jj].max[1], domain[jj].min[2], domain[jj].max[2], mpi.rank);
      }
      for(int jj = 0; jj < overlapNum; jj++){
	fprintf(stderr, "sendBuf[%d]: xmin = %f, xmax = %f, ymin = %f, ymax = %f, zmin = %f, zmax = %f @ rank %d\n",
		jj, sendBuf[jj].xmin, sendBuf[jj].xmax, sendBuf[jj].ymin, sendBuf[jj].ymax, sendBuf[jj].zmin, sendBuf[jj].zmax, mpi.rank);
      }
      fprintf(stderr, "ii = %d: x = %f, y = %f, z = %f @ rank %d\n", ii, src.pos[ii].x, src.pos[ii].y, src.pos[ii].z, mpi.rank);
#endif
      __KILL__(stderr, "ERROR: target MPI rank is missing\n");
    }
#endif//NDEBUG
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < numOld; ii++){ */
  //-----------------------------------------------------------------------
  sendBuf[0].head = 0;
  for(int ii = 0; ii < overlapNum; ii++){
    chkMPIerr(MPI_Isend(&(sendBuf[ii].num), 1, MPI_INT, sendBuf[ii].rank, mpi.rank, mpi.comm, &(domain[sendBuf[ii].rank].req)));
    if( ii > 0 )
      sendBuf[ii].head = sendBuf[ii - 1].head + sendBuf[ii - 1].num;
  }/* for(int ii = 0; ii < overlapNum; ii++){ */
  //-----------------------------------------------------------------------
  if( (sendBuf[overlapNum - 1].head + sendBuf[overlapNum - 1].num) != numOld ){
    __KILL__(stderr, "ERROR: total number of scattered particles (%d) is differ from that of local particles (%d)\n",
	     sendBuf[overlapNum - 1].head + sendBuf[overlapNum - 1].num, numOld);
  }
  //-----------------------------------------------------------------------
  qsort(key, numOld, sizeof(domainDecomposeKey), ddkeyAscendingOrder);
  //-----------------------------------------------------------------------
  /* sort N-body particle arrays: dst is send buffer */
  for(int ii = 0; ii < numOld; ii++)    dst.pos [ii] = src.pos [key[ii].bodyIdx];
#ifdef  GADGET_MAC
  for(int ii = 0; ii < numOld; ii++)    dst.acc [ii] = src.acc [key[ii].bodyIdx];
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
  for(int ii = 0; ii < numOld; ii++)    dst.vel [ii] = src.vel [key[ii].bodyIdx];
  for(int ii = 0; ii < numOld; ii++)    dst.time[ii] = src.time[key[ii].bodyIdx];
#else///BLOCK_TIME_STEP
  for(int ii = 0; ii < numOld; ii++)    dst.vx  [ii] = src.vx  [key[ii].bodyIdx];
  for(int ii = 0; ii < numOld; ii++)    dst.vy  [ii] = src.vy  [key[ii].bodyIdx];
  for(int ii = 0; ii < numOld; ii++)    dst.vz  [ii] = src.vz  [key[ii].bodyIdx];
#endif//BLOCK_TIME_STEP
  for(int ii = 0; ii < numOld; ii++)    dst.idx [ii] = src.idx [key[ii].bodyIdx];
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < numProcs; ii++){
    MPI_Status status;
    chkMPIerr(MPI_Wait(&(domain[ii].req), &status));
  }/* for(int ii = 0; ii < numProcs; ii++){ */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < overlapNum; ii++){
    //---------------------------------------------------------------------
    chkMPIerr(MPI_Isend(&(dst.pos [sendBuf[ii].head]), sendBuf[ii].num, mpi.ipos         , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].pos)));
#ifdef  GADGET_MAC
    chkMPIerr(MPI_Isend(&(dst.acc [sendBuf[ii].head]), sendBuf[ii].num, mpi.iacc         , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].acc)));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
    chkMPIerr(MPI_Isend(&(dst.vel [sendBuf[ii].head]), sendBuf[ii].num, mpi.ivel         , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].vel)));
    chkMPIerr(MPI_Isend(&(dst.time[sendBuf[ii].head]), sendBuf[ii].num, mpi.time         , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].time)));
#else///BLOCK_TIME_STEP
    chkMPIerr(MPI_Isend(&(dst.vx  [sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT      , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].vx )));
    chkMPIerr(MPI_Isend(&(dst.vy  [sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT      , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].vy )));
    chkMPIerr(MPI_Isend(&(dst.vz  [sendBuf[ii].head]), sendBuf[ii].num, MPI_REALDAT      , sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].vz )));
#endif//BLOCK_TIME_STEP
    chkMPIerr(MPI_Isend(&(dst.idx [sendBuf[ii].head]), sendBuf[ii].num, MPI_UNSIGNED_LONG, sendBuf[ii].rank, mpi.rank, mpi.comm, &(sendBuf[ii].idx)));
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < overlapNum; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* src is receive buffer */
  //-----------------------------------------------------------------------
  *numNew = 0;
  for(int ii = 0; ii < numProcs; ii++)
    recvBuf[ii].head = 0;
  for(int ii = 0; ii < numProcs; ii++){
    //---------------------------------------------------------------------
    /* receive recvNum */
    MPI_Status status;
    chkMPIerr(MPI_Wait(&(recvBuf[ii].req), &status));
    //---------------------------------------------------------------------
    /* if recvNum != 0, then set receive buffer */
    if( recvBuf[ii].num != 0 ){
      //-------------------------------------------------------------------
      chkMPIerr(MPI_Irecv(&(src.pos [recvBuf[ii].head]), recvBuf[ii].num, mpi.ipos         , ii, ii, mpi.comm, &(recvBuf[ii].pos)));
#ifdef  GADGET_MAC
      chkMPIerr(MPI_Irecv(&(src.acc [recvBuf[ii].head]), recvBuf[ii].num, mpi.iacc         , ii, ii, mpi.comm, &(recvBuf[ii].acc)));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
      chkMPIerr(MPI_Irecv(&(src.vel [recvBuf[ii].head]), recvBuf[ii].num, mpi.ipos         , ii, ii, mpi.comm, &(recvBuf[ii].vel)));
      chkMPIerr(MPI_Irecv(&(src.time[recvBuf[ii].head]), recvBuf[ii].num, mpi.ipos         , ii, ii, mpi.comm, &(recvBuf[ii].time)));
#else///BLOCK_TIME_STEP
      chkMPIerr(MPI_Irecv(&(src.vx  [recvBuf[ii].head]), recvBuf[ii].num, MPI_REALDAT      , ii, ii, mpi.comm, &(recvBuf[ii].vx )));
      chkMPIerr(MPI_Irecv(&(src.vy  [recvBuf[ii].head]), recvBuf[ii].num, MPI_REALDAT      , ii, ii, mpi.comm, &(recvBuf[ii].vy )));
      chkMPIerr(MPI_Irecv(&(src.vz  [recvBuf[ii].head]), recvBuf[ii].num, MPI_REALDAT      , ii, ii, mpi.comm, &(recvBuf[ii].vz )));
#endif//BLOCK_TIME_STEP
      chkMPIerr(MPI_Irecv(&(src.idx [recvBuf[ii].head]), recvBuf[ii].num, MPI_UNSIGNED_LONG, ii, ii, mpi.comm, &(recvBuf[ii].idx)));
      //-------------------------------------------------------------------
      *numNew += recvBuf[ii].num;
      //-------------------------------------------------------------------
      if( ii + 1 < numProcs )
	recvBuf[ii + 1].head = recvBuf[ii].head + recvBuf[ii].num;
      //-------------------------------------------------------------------
    }/* if( recvBuf[ii].num != 0 ){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < numProcs; ii++){ */
  //-----------------------------------------------------------------------
  if( *numNew > numMax ){
    __KILL__(stderr, "ERROR: # of required receive buffer (%d) exceeds the maximum number of particles per process (%d).\n\tsuggestion: consider increasing \"MAX_FACTOR_FROM_EQUIPARTITION\" defined in src/para/mpicfg.h (current value is %f) to at least %f.\n", *numNew, numMax, MAX_FACTOR_FROM_EQUIPARTITION, MAX_FACTOR_FROM_EQUIPARTITION * (float)(*numNew) / (float)numMax);
  }/* if( *numNew > numMax ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* complete MPI communications */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < overlapNum; ii++){
    //---------------------------------------------------------------------
    MPI_Status  pos;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. pos), &pos));
#ifdef  GADGET_MAC
    MPI_Status  acc;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. acc), &acc));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
    MPI_Status  vel;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. vel), &vel));
    MPI_Status time;    chkMPIerr(MPI_Wait(&(sendBuf[ii].time), &time));
#else///BLOCK_TIME_STEP
    MPI_Status   vx;    chkMPIerr(MPI_Wait(&(sendBuf[ii].  vx), & vx));
    MPI_Status   vy;    chkMPIerr(MPI_Wait(&(sendBuf[ii].  vy), & vy));
    MPI_Status   vz;    chkMPIerr(MPI_Wait(&(sendBuf[ii].  vz), & vz));
#endif//BLOCK_TIME_STEP
    MPI_Status  idx;    chkMPIerr(MPI_Wait(&(sendBuf[ii]. idx), &idx));
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < overlapNum; ii++){ */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < numProcs; ii++)
    if( recvBuf[ii].num != 0 ){
      //-------------------------------------------------------------------
      MPI_Status  pos;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. pos), &pos));
#ifdef  GADGET_MAC
      MPI_Status  acc;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. acc), &acc));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
      MPI_Status  vel;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. vel), &vel));
      MPI_Status time;      chkMPIerr(MPI_Wait(&(recvBuf[ii].time), &time));
#else///BLOCK_TIME_STEP
      MPI_Status   vx;      chkMPIerr(MPI_Wait(&(recvBuf[ii].  vx), & vx));
      MPI_Status   vy;      chkMPIerr(MPI_Wait(&(recvBuf[ii].  vy), & vy));
      MPI_Status   vz;      chkMPIerr(MPI_Wait(&(recvBuf[ii].  vz), & vz));
#endif//BLOCK_TIME_STEP
      MPI_Status  idx;      chkMPIerr(MPI_Wait(&(recvBuf[ii]. idx), &idx));
      //-------------------------------------------------------------------
    }/* if( recvBuf[ii].num != 0 ){ */
  //-----------------------------------------------------------------------
  /* confirmation */
  const int diff = (*numNew) - numOld;
  int diff_sum;
  chkMPIerr(MPI_Reduce(&diff, &diff_sum, 1, MPI_INT, MPI_SUM, 0, mpi.comm));
  if( mpi.rank == 0 )
    if( diff_sum != 0 ){
      __KILL__(stderr, "ERROR: domain decomposition cause some error (duplication of %d particles)\n", diff_sum);
    }/* if( diff_sum != 0 ){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* copy N-body particles from device to host */
  //-----------------------------------------------------------------------
#ifdef  GENERATE_PHKEY_ON_DEVICE
  copyParticle_hst2dev(*numNew, src, ibody_dev
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );
#endif//GENERATE_PHKEY_ON_DEVICE
  //-----------------------------------------------------------------------
  __NOTE__("numOld = %d, numNew = %d @ rank %d\n", numOld, *numNew, mpi.rank);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
