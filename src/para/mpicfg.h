/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:37:07
 *                                                                       *
 *    Header File for N-body calculation with MPI parallelization        *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef MPICFG_H
#define MPICFG_H
//-------------------------------------------------------------------------
#include <mpi.h>
//-------------------------------------------------------------------------
#include "macro.h"
#include "mpilib.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef GPUS_PER_PROCESS
#define GPUS_PER_PROCESS (1)
#endif//GPUS_PER_PROCESS
//-------------------------------------------------------------------------
#define MAX_FACTOR_INCREASE (1.3f)
/* #define MAX_FACTOR_INCREASE (1.5f) */
//-------------------------------------------------------------------------
#define MAX_FACTOR_SAFETY (1.05f)
/* #define MAX_FACTOR_SAFETY (1.1f) */
#define MAX_FACTOR_FROM_EQUIPARTITION (MAX_FACTOR_INCREASE * MAX_FACTOR_SAFETY)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
typedef struct
{
  MPI_Datatype cell, link;
  MPI_Comm comm;
  int rank, size;
} MPIcfg_octree;
//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "mpicfg.c"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  void setNodeConfig(ulong Ntot, int *Nnode, int *Ni, MPIinfo mpi, MPIcfg_tree *let, const int devID);
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//MPICFG_H
//-------------------------------------------------------------------------
