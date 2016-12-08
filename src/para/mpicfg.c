/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:36:57
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
#include <math.h>
#include <mpi.h>
//-------------------------------------------------------------------------
#include "macro.h"
#include "name.h"
#include "mpilib.h"
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#include "../tree/make.h"
//-------------------------------------------------------------------------
#include "mpicfg.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void createMPIcfg_tree(MPIcfg_tree *let, MPIinfo mpi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  let->comm = mpi.comm;
  let->size = mpi.size;
  let->rank = mpi.rank;
  //-----------------------------------------------------------------------
#if 1
  int num;
  static int blck[4];
  static MPI_Aint disp[4];
  static MPI_Datatype type[4];
  /* commit position */
  /* typedef struct __align__(16){  real x, y, z, m;} position; */
  num = 1;
  type[0] = MPI_REALDAT;
  blck[0] = 4;
  disp[0] = 0;
  commitMPIstruct(num, blck, disp, type, &(let->ipos));
#ifdef  GADGET_MAC
  /* commit acceleration */
  /* typedef struct __align__(16){  real x, y, z, pot;} acceleration; */
  num = 1;
  type[0] = MPI_REALDAT;
  blck[0] = 4;
  disp[0] = 0;
  commitMPIstruct(num, blck, disp, type, &(let->iacc));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
  /* commit velocity */
  /* typedef struct __align__(16){  real x, y, z, dt;} velocity; */
  num = 1;
  type[0] = MPI_REALDAT;
  blck[0] = 4;
  disp[0] = 0;
  commitMPIstruct(num, blck, disp, type, &(let->ivel));
  /* commit ibody_time */
  /* typedef struct __align__(16){  double t0, t1;} ibody_time; */
  num = 1;
  type[0] = MPI_DOUBLE;
  blck[0] = 2;
  disp[0] = 0;
  commitMPIstruct(num, blck, disp, type, &(let->time));
#endif//BLOCK_TIME_STEP
  /* commit jparticle */
  /* typedef struct __align__(16){  real y, x, w, z;} jparticle; */
  num = 1;
  type[0] = MPI_REALDAT;
  blck[0] = 4;
  disp[0] = 0;
  commitMPIstruct(num, blck, disp, type, &(let->jpos));
  /* commit uint */
  let->more = MPI_UNSIGNED;
  /* commit jmass */
#ifdef  INDIVIDUAL_GRAVITATIONAL_SOFTENING
  /* typedef struct __align__(8){  real mass, eps2;} jmass; */
  num = 1;
  type[0] = MPI_REALDAT;
  blck[0] = 2;
  disp[0] = 0;
  commitMPIstruct(num, blck, disp, type, &(let->mass));
#else///INDIVIDUAL_GRAVITATIONAL_SOFTENING
  /* typedef real jmass; */
  let->mass = MPI_REALDAT;
#endif//INDIVIDUAL_GRAVITATIONAL_SOFTENING
#else
  commitMPIbyte(&(let->ipos), (int)sizeof(position));
#ifdef  GADGET_MAC
  commitMPIbyte(&(let->iacc), (int)sizeof(acceleration));
#endif//GADGET_MAC
#ifdef  BLOCK_TIME_STEP
  commitMPIbyte(&(let->ivel), (int)sizeof(velocity));
  commitMPIbyte(&(let->time), (int)sizeof(ibody_time));
#endif//BLOCK_TIME_STEP
  commitMPIbyte(&(let->jpos), (int)sizeof(jparticle));
  commitMPIbyte(&(let->more), (int)sizeof(uint));
  commitMPIbyte(&(let->mass), (int)sizeof(jmass));
#endif
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void setNodeConfig(ulong Ntot, int *Nnode, int *Ni, MPIinfo mpi, MPIcfg_tree *let, const int devID)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* determine # of N-body particles within an MPI process */
  //-----------------------------------------------------------------------
  static ulong Nini, Nfin;
  Nini = (Ntot * (    (ulong)mpi.rank)) / (ulong)mpi.size;
  Nfin = (Ntot * (1 + (ulong)mpi.rank)) / (ulong)mpi.size;
  Nfin -= Nini;
  *Nnode = (int)Nfin;
  *Ni = *Nnode;
  //-----------------------------------------------------------------------
  if( *Nnode > (int)ceilf((float)NUM_BODY_MAX * (float)GPUS_PER_PROCESS * MAX_FACTOR_FROM_EQUIPARTITION) ){
    __KILL__(stderr, "ERROR: # of particles per process (%d) exceeds the limit (%d)\n", *Nnode, (int)ceilf((float)NUM_BODY_MAX * (float)GPUS_PER_PROCESS * MAX_FACTOR_FROM_EQUIPARTITION));
  }/* if( *Nnode > (int)ceilf((float)NUM_BODY_MAX * (float)GPUS_PER_PROCESS * MAX_FACTOR_FROM_EQUIPARTITION) ){ */
  if( *Nnode > (int)ceilf((float)NUM_BODY_MAX * (float)GPUS_PER_PROCESS * MAX_FACTOR_FROM_EQUIPARTITION * 0.8f) ){
    __FPRINTF__(stdout, "WARNING: # of particles per process (%d) exceeds the 80%% of the limit (%d)\n", *Nnode, (int)ceilf((float)NUM_BODY_MAX * (float)GPUS_PER_PROCESS * MAX_FACTOR_FROM_EQUIPARTITION * 0.8f));
  }/* if( *Nnode > (int)ceilf((float)NUM_BODY_MAX * (float)GPUS_PER_PROCESS * MAX_FACTOR_FROM_EQUIPARTITION * 0.8f) ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* configure structures for MPI communications */
  //-----------------------------------------------------------------------
  createMPIcfg_tree(let, mpi);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* make a simple log of MPI */
  //-----------------------------------------------------------------------
  static char filename[128];
  static FILE *fp;
  sprintf(filename, "%s/mpi%.4d.txt", DOCUMENTFOLDER, mpi.rank);
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  //-----------------------------------------------------------------------
  /* write fundamental information on MPI */
  static char name[MPI_MAX_PROCESSOR_NAME];
  int tmp;
  chkMPIerr(MPI_Get_processor_name(name, &tmp));
  fprintf(fp, "host name: %s\n", name);
  fprintf(fp, "\n");
  fprintf(fp, "MPI info: rank %d out of %d processes\n", mpi.rank, mpi.size);
  fprintf(fp, "GPU info: device ID is %d\n", devID);
  fprintf(fp, "\n");
  //-----------------------------------------------------------------------
  /* write information on domain decomposition */
  fprintf(fp, "# of N-body particles is %d out of %zu\n", *Nnode, Ntot);
  fprintf(fp, "\n");
  //-----------------------------------------------------------------------
  /* write fundamental information on MPI topology */
  /* fprintf(fp, "# of dimensions for MPI topology: %d\n", jcfg->ndim); */
  /* for(int i = 0; i < jcfg->ndim; i++) */
  /*   fprintf(fp, "dim[%d] = %d\n", i, jcfg->dim[i]); */
  /* fprintf(fp, "\n"); */
  //-----------------------------------------------------------------------
  /* /\* test use of the newly created MPI topology *\/ */
  /* static int Njold, Njnew; */
  /* Njold = *Nnode; */
  /* if( jcfg->dim[1] > 1 ){ */
  /*   for(int i = 0; i < jcfg->ilog2dim1; i++){ */
  /*     int pair; */
  /*     setAccumulationPair(jcfg, 1 << i, &pair, Njold, &Njnew); */
  /*     fprintf(fp, "%d-th accumulation phase:\n", i); */
  /*     fprintf(fp, "\tprocess rank of the pair: %d\n", pair); */
  /*     fprintf(fp, "\t# of new j-particles: %d\n", Njnew); */
  /*     Njold += Njnew; */
  /*   } */
  /* } */
  /* if( jcfg->dim[0] > 1 ) */
  /*   for(int i = 0; i < jcfg->dim[0]; i++){ */
  /*     int src, dst; */
  /*     setTransferPair(*jcfg, &dst, Njold, &src, &Njnew); */
  /*     fprintf(fp, "%d-th transfer phase:\n", i); */
  /*     fprintf(fp, "\tdata of %d particles go toward %d-th process\n", Njold, dst); */
  /*     fprintf(fp, "\tdata of %d particles come from %d-th process\n", Njnew, src); */
  /*     Njold = Njnew; */
  /*   } */
  //-----------------------------------------------------------------------
  fclose(fp);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
