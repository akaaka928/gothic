/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/26(Tue) 15:37:04
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


//-------------------------------------------------------------------------
#   if  !defined(MPI_INCLUDED) && !defined(OMPI_MPI_H)
#       include <mpi.h>
#endif//!defined(MPI_INCLUDED) && !defined(OMPI_MPI_H)
//-------------------------------------------------------------------------
#ifndef MACRO_H
#       include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef MPILIB_H
#       include <mpilib.h>
#endif//MPILIB_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef GPUS_PER_PROCESS
#define GPUS_PER_PROCESS (1)
#endif//GPUS_PER_PROCESS
//-------------------------------------------------------------------------
#define MAX_FACTOR_FROM_EQUIPARTITION (2.0f)
/* #define MAX_FACTOR_FROM_EQUIPARTITION (1.9f) */
/* #define MAX_FACTOR_FROM_EQUIPARTITION (1.8f) */
/* #define MAX_FACTOR_FROM_EQUIPARTITION (1.5f) */
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
  void setNodeConfig(ulong Ntot, int *Nnode, int *Ni, MPIinfo mpi, MPIcfg_tree *let);
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//MPICFG_H
//-------------------------------------------------------------------------
