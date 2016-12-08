/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:35:32
 *                                                                       *
 *    Implementations related to domain decomposition (MPI)              *
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
#include "macro.h"
#include "mpilib.h"
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#ifndef EXCHANGE_USING_GPUS
#include "../misc/tune.h"
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#include "../misc/brent.h"
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
#include "../time/adv_dev.h"
#endif//EXCHANGE_USING_GPUS
//-------------------------------------------------------------------------
#include "mpicfg.h"
#include "exchange.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* function to set MPI topology for Orthogonal Recursive Multi-section */
//-------------------------------------------------------------------------
#ifdef  EXCHANGE_USING_GPUS
//-------------------------------------------------------------------------
muse allocateORMtopology(float **dxmin, float **dxmax, float **dymin, float **dymax, float **dzmin, float **dzmax, MPI_Request **dmreq,
			 float **sxmin, float **sxmax, float **symin, float **symax, float **szmin, float **szmax,
			 sendCfg **sendBuf, recvCfg **recvBuf, int **rnum, int **disp,
			 MPIinfo orm[], MPIinfo rep[],
			 const int Nx, const int Ny, const int Nz, MPIcfg_tree *mpi,
			 domainCfg *dom, sampling *sample, const ulong Ntot)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  muse alloc = {0, 0};
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  mpi->dim[0] = (Nx > 0) ? Nx : 1;  mpi->prd[0] = false;
  mpi->dim[1] = (Ny > 0) ? Ny : 1;  mpi->prd[1] = false;
  mpi->dim[2] = (Nz > 0) ? Nz : 1;  mpi->prd[2] = false;
  //-----------------------------------------------------------------------
  chkMPIerr(MPI_Cart_create(mpi->comm, 3, mpi->dim, mpi->prd, false, &(mpi->cart)));
  chkMPIerr(MPI_Cart_coords(mpi->cart, mpi->rank, 3, mpi->pos));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* MPI process decomposition */
  //-----------------------------------------------------------------------
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
  //-----------------------------------------------------------------------
  MPIinfo org;
  org.comm = mpi->cart;
  org.size = mpi->size;
  org.rank = mpi->rank;
  for(int ii = 0; ii < 3; ii++){
    splitMPI(org.comm, mpi->pos[ii], org.rank, &(orm[ii]));
    splitMPI(org.comm, orm[ii].rank, org.rank, &(rep[ii]));
    org = orm[ii];
  }/* for(int ii = 0; ii < 3; ii++){ */
#if 0
  fprintf(stdout, "rank %d: dim[0] = %d, dim[1] = %d, dim[2] = %d, pos[0] = %d, pos[1] = %d, pos[2] = %d\n", mpi->rank, mpi->dim[0], mpi->dim[1], mpi->dim[2], mpi->pos[0], mpi->pos[1], mpi->pos[2]);
  fprintf(stdout, "rank: mpi(%d/%d), orm[0](%d/%d), orm[1](%d/%d), orm[2](%d/%d), rep[0](%d/%d), rep[1](%d/%d), rep[2](%d/%d)\n", mpi->rank, mpi->size,
	  orm[0].rank, orm[0].size, orm[1].rank, orm[1].size, orm[2].rank, orm[2].size,
	  rep[0].rank, rep[0].size, rep[1].rank, rep[1].size, rep[2].rank, rep[2].size);
  fflush(stdout);
  MPI_Finalize();
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* memory allocation for ORM related arrays */
  //-----------------------------------------------------------------------
  *sendBuf = (sendCfg *)malloc(mpi->size * sizeof(sendCfg));  alloc.host += mpi->size * sizeof(sendCfg);
  *recvBuf = (recvCfg *)malloc(mpi->size * sizeof(recvCfg));  alloc.host += mpi->size * sizeof(recvCfg);
  //-----------------------------------------------------------------------
  if( *sendBuf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sendBuf\n");  }
  if( *recvBuf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate recvBuf\n");  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* memory allocation for arrays to set domain boundaries */
  //-----------------------------------------------------------------------
  *sxmin = (float *)malloc(mpi->dim[0] * sizeof(float));  alloc.host += mpi->dim[0] * sizeof(float);
  *sxmax = (float *)malloc(mpi->dim[0] * sizeof(float));  alloc.host += mpi->dim[0] * sizeof(float);
  *symin = (float *)malloc(mpi->dim[1] * sizeof(float));  alloc.host += mpi->dim[1] * sizeof(float);
  *symax = (float *)malloc(mpi->dim[1] * sizeof(float));  alloc.host += mpi->dim[1] * sizeof(float);
  *szmin = (float *)malloc(mpi->dim[2] * sizeof(float));  alloc.host += mpi->dim[2] * sizeof(float);
  *szmax = (float *)malloc(mpi->dim[2] * sizeof(float));  alloc.host += mpi->dim[2] * sizeof(float);
  //-----------------------------------------------------------------------
  if( *sxmin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sxmin\n");  }
  if( *sxmax == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sxmax\n");  }
  if( *symin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate symin\n");  }
  if( *symax == NULL ){    __KILL__(stderr, "ERROR: failure to allocate symax\n");  }
  if( *szmin == NULL ){    __KILL__(stderr, "ERROR: failure to allocate szmin\n");  }
  if( *szmax == NULL ){    __KILL__(stderr, "ERROR: failure to allocate szmax\n");  }
  //-----------------------------------------------------------------------
  sample->xmin = *sxmin;  sample->ymin = *symin;  sample->zmin = *szmin;
  sample->xmax = *sxmax;  sample->ymax = *symax;  sample->zmax = *szmax;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set sampling rate and allocate corresponding size of temporary array */
  //-----------------------------------------------------------------------

  /* sample->rate = DEFAULT_SAMPLING_RATE; */
  /* if( Ntot * DEFAULT_SAMPLING_RATE < mpi->size * DEFAULT_SAMPLING_NUMBER ) */
  /*   sample->rate = DEFAULT_SAMPLING_NUMBER * (float)mpi->size / (float)Ntot; */

  sample->rate = (((float)Ntot * DEFAULT_SAMPLING_RATE) >= ((float)mpi->size * DEFAULT_SAMPLING_NUMBER)) ?
    (DEFAULT_SAMPLING_RATE) : (DEFAULT_SAMPLING_NUMBER * (float)mpi->size / (float)Ntot);

  //-----------------------------------------------------------------------
  bool root = (mpi->rank != 0) ? false : true;
  for(int ii = 0; ii < 3; ii++)
    if( (orm[ii].rank == 0) && (orm[ii].size > 1) )
      root = true;
  //-----------------------------------------------------------------------
  if( root ){
    //---------------------------------------------------------------------
    /* recvcnts, displs */
    *rnum = (int *)malloc(mpi->size * sizeof(int));    alloc.host += mpi->size * sizeof(int);    if( *rnum == NULL ){      __KILL__(stderr, "ERROR: failure to allocate rnum\n");    }
    *disp = (int *)malloc(mpi->size * sizeof(int));    alloc.host += mpi->size * sizeof(int);    if( *disp == NULL ){      __KILL__(stderr, "ERROR: failure to allocate disp\n");    }
    sample->rnum = *rnum;
    sample->disp = *disp;
    //---------------------------------------------------------------------
  }/* if( root ){ */
  //-----------------------------------------------------------------------
  /* sample->Nmax = (int)ceilf((float)NUM_BODY_MAX * MAX_FACTOR_FROM_EQUIPARTITION * sample->rate); */
  sample->Nmax = (int)ceilf((float)Ntot * MAX_FACTOR_FROM_EQUIPARTITION * sample->rate);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
  return (alloc);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  releaseORMtopology(float  *dxmin, float  *dxmax, float  *dymin, float  *dymax, float  *dzmin, float  *dzmax, MPI_Request  *dmreq,
			 float  *sxmin, float  *sxmax, float  *symin, float  *symax, float  *szmin, float  *szmax,
			 sendCfg  *sendBuf, recvCfg  *recvBuf, int  *rnum, int  *disp,
			 MPIinfo orm[], MPIinfo rep[], const int rank)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int ii = 2; ii >= 0; ii--){
    freeMPIgroup(&(rep[ii]));
    freeMPIgroup(&(orm[ii]));
  }/* for(int ii = 2; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  free(dxmin);  free(dymin);  free(dzmin);
  free(dxmax);  free(dymax);  free(dzmax);  free(dmreq);
  //-----------------------------------------------------------------------
  free(sxmin);	free(symin);  free(szmin);
  free(sxmax);	free(symax);  free(szmax);
  //-----------------------------------------------------------------------
  free(sendBuf);
  free(recvBuf);
  //-----------------------------------------------------------------------
  if( rank == 0 ){
    free(rnum);
    free(disp);
  }/* if( rank == 0 ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//EXCHANGE_USING_GPUS
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef EXCHANGE_USING_GPUS
//-------------------------------------------------------------------------
void configORBtopology
(domainCfg **box, sendCfg **sendBuf, recvCfg **recvBuf, MPIcfg_tree *mpi, int *ndim, MPIinfo **orb,
 const ulong Ntot, real *samplingRate, int *sampleNumMax, real **sampleLoc, real **sampleFul,
 int **recvNum, int **recvDsp, real **boxMin, real **boxMax)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

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
  }/* for(int ii = 0; ii < *ndim; ii++){ */
  //-----------------------------------------------------------------------
  const int Ndim = *ndim;
  //-----------------------------------------------------------------------
  chkMPIerr(MPI_Cart_create(mpi->comm, Ndim, mpi->dim, mpi->prd, false, &(mpi->cart)));
  chkMPIerr(MPI_Cart_coords(mpi->cart, mpi->rank, Ndim, mpi->pos));
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* MPI process decomposition */
  //-----------------------------------------------------------------------
  *orb = (MPIinfo *)malloc((Ndim) * sizeof(MPIinfo));
  if( *orb == NULL ){    __KILL__(stderr, "ERROR: failure to allocate orb\n");  }
  //-----------------------------------------------------------------------
  *box = (domainCfg *)malloc((mpi->size) * sizeof(domainCfg));
  if( *box == NULL ){    __KILL__(stderr, "ERROR: failure to allocate box\n");  }
  //-----------------------------------------------------------------------
  MPIinfo old;
  old.comm = mpi->cart;
  old.size = mpi->size;
  old.rank = mpi->rank;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < Ndim; ii++){
    //---------------------------------------------------------------------
    splitMPI(old.comm, mpi->pos[ii], old.rank, &((*orb)[ii]));    old = (*orb)[ii];
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < Ndim; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* memory allocation for ORB related arrays */
  //-----------------------------------------------------------------------
  *sendBuf = (sendCfg *)malloc((mpi->size) * sizeof(sendCfg));
  *recvBuf = (recvCfg *)malloc((mpi->size) * sizeof(recvCfg));
  //-----------------------------------------------------------------------
  if( *sendBuf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sendBuf\n");  }
  if( *recvBuf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate recvBuf\n");  }
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
  for(int ii = 0; ii < Ndim; ii++)
    if( (*orb)[ii].rank == 0 )
      root = true;
  //-----------------------------------------------------------------------
  if( root ){
    //---------------------------------------------------------------------
    *sampleNumMax = (int)ceilf((float)Ntot * (*samplingRate));
    *sampleFul = (real *)malloc((*sampleNumMax) * sizeof(real));
    if( *sampleFul == NULL ){      __KILL__(stderr, "ERROR: failure to allocate sampleFul\n");    }
    //---------------------------------------------------------------------
    /* recvcnts, displs */
    *recvNum = (int *)malloc(mpi->size * sizeof(int));
    *recvDsp = (int *)malloc(mpi->size * sizeof(int));
    if( *recvNum == NULL ){      __KILL__(stderr, "ERROR: failure to allocate recvNum\n");    }
    if( *recvDsp == NULL ){      __KILL__(stderr, "ERROR: failure to allocate recvDsp\n");    }
    //---------------------------------------------------------------------
  }/* if( root ){ */
  //-----------------------------------------------------------------------
  /* *sampleNumMax = (int)ceilf((float)NUM_BODY_MAX * MAX_FACTOR_FROM_EQUIPARTITION * (*samplingRate)); */
  *sampleNumMax = (int)ceilf((float)Ntot * MAX_FACTOR_FROM_EQUIPARTITION * (*samplingRate));
  *sampleLoc = (real *)malloc((*sampleNumMax) * 3 * sizeof(real));
  if( *sampleLoc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sampleLoc\n");  }
  //-----------------------------------------------------------------------
  /* minimum/maximun position of the decomposed domain */
  int maxDim = (mpi->dim[0] > mpi->dim[1]) ? mpi->dim[0] : mpi->dim[1];
  if( maxDim < mpi->dim[2] )      maxDim = mpi->dim[2];
  *boxMin = (real *)malloc(maxDim * sizeof(real));
  *boxMax = (real *)malloc(maxDim * sizeof(real));
  if( *boxMin == NULL ){      __KILL__(stderr, "ERROR: failure to allocate boxMin\n");    }
  if( *boxMax == NULL ){      __KILL__(stderr, "ERROR: failure to allocate boxMax\n");    }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void removeORBtopology
(domainCfg  *box, sendCfg  *sendBuf, recvCfg  *recvBuf, const int ndim, MPIinfo  *orb,
 real  *sampleLoc, real  *sampleFul, int  *recvNum, int  *recvDsp, real  *boxMin, real  *boxMax)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndim; ii++)
    freeMPIgroup(&(orb[ii]));
  free(orb);
  //-----------------------------------------------------------------------
  free(box);
  //-----------------------------------------------------------------------
  free(sendBuf);
  free(recvBuf);
  //-----------------------------------------------------------------------
  if( sampleFul != NULL ){
    free(sampleFul);
    free(recvNum);
    free(recvDsp);
  }/* if( sampleFul != NULL ){ */
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
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          (*(real *)a) > (*(real *)b) ){    return ( 1);  }
  else{    if( (*(real *)a) < (*(real *)b) ){    return (-1);  }
    else                                         return ( 0);  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
int ddkeyAscendingOrder(const void *a, const void *b);
int ddkeyAscendingOrder(const void *a, const void *b)
{
  //-----------------------------------------------------------------------
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((domainDecomposeKey *)a)->dstRank > ((domainDecomposeKey *)b)->dstRank ){    return ( 1);  }
  else{    if( ((domainDecomposeKey *)a)->dstRank < ((domainDecomposeKey *)b)->dstRank ){    return (-1);  }
    else                                                         return ( 0);  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
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
(const int numOld, const int numMax, int *numNew,  iparticle ibody_dev, iparticle src, iparticle dst,
 domainDecomposeKey *key,
 sendCfg *sendBuf, recvCfg *recvBuf,
 const float samplingRate, const int sampleNumMax, const double tloc, real *sampleLoc, real *sampleFul, int *recvNum, int *recvDsp, real *boxMin, real *boxMax,
 const int ndim, MPIinfo *orb, domainCfg *domain, MPIcfg_tree mpi
   , measuredTime *measured
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
   , brentStatus *status, brentMemory *memory
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
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
  copyParticle_dev2hst(numOld, ibody_dev, src
#ifdef  EXEC_BENCHMARK
		       , execTime
#endif//EXEC_BENCHMARK
		       );
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
#if 0
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
  /* modify parameters related to auto-tuning */
  //-----------------------------------------------------------------------
  const double scale = (double)(*numNew) / (double)numOld;
  //-----------------------------------------------------------------------
  measured->walkTree[0] *= scale;
  measured->walkTree[1] *= scale;
  measured->makeTree    *= scale;
#ifdef  WALK_TREE_TOTAL_SUM_MODEL
  measured->incSum  *= scale;
#endif//WALK_TREE_TOTAL_SUM_MODEL
  measured->genTree = 0.0;
  measured->calcAcc = 0.0;
  measured->calcMAC = 0.0;
#ifdef  MONITOR_LETGEN_TIME
  measured->makeLET = 0.0;
#endif//MONITOR_LETGEN_TIME
  //-----------------------------------------------------------------------
#   if  defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
  status->x.val *= scale;
  status->w.val *= scale;
  status->v.val *= scale;
  status->u.val *= scale;
  memory->previous *= scale;
#endif//defined(LOCALIZE_I_PARTICLES) && defined(USE_BRENT_METHOD)
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//EXCHANGE_USING_GPUS
//-------------------------------------------------------------------------
