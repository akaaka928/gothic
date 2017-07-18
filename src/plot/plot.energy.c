/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/06(Tue) 12:40:19
 *                                                                       *
 *    Plot Code of N-body Simulations (using PLplot)                     *
 *      Time Evolution of total energy, kinetic energy, potential energy *
 *      Conservation of total energy                                     *
 *      Time Evolution of Virial Ratio (-K/W)                            *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#define NORMALIZED_MOMENTUM_ERROR
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
#include "macro.h"
#include "myutil.h"
#include "constants.h"
#include "mpilib.h"
#include "plplotlib.h"
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../file/io.h"
//-------------------------------------------------------------------------
extern const double     time2astro;extern const char     time_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double   energy2astro;extern const char   energy_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double     mass2astro;extern const char     mass_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double velocity2astro;extern const char velocity_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotTimeEvolution
(PLINT num, PLFLT *time, PLFLT *Ekin, PLFLT *Epot, PLFLT *Etot,
 PLplotPltRange evolve, PLplotPltRange virial, PLplotPltRange relerr,
 char *file, int argc, char **argv,
 real *worst, real *final);
//-------------------------------------------------------------------------
void plotMomentumError
(PLINT num, PLFLT *time, PLFLT *momx, PLFLT *momy, PLFLT *momz, PLplotPltRange momerr,
 char *file, int argc, char **argv,
 real worst[], real final[]);
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //-----------------------------------------------------------------------
  /* parallelized region employing MPI start */
  //-----------------------------------------------------------------------
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* initialization */
  //-----------------------------------------------------------------------
  if( argc < 6 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 6);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -problem=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }
  //-----------------------------------------------------------------------
  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));
  int  problem;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,  "problem", &problem));
  //-----------------------------------------------------------------------
  modifyArgcArgv4PLplot(&argc, argv, 6);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  int last;
  readConfigFile(&last, file);
  int unit;
  ulong Ntot;
  real eta, eps;
  double ft, snapshotInterval, saveInterval;
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file, mpi);
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
  nbody_hdf5 body;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
  ulong *hdf5_idx;
  allocSnapshotArray(&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx, (int)Ntot, &body);
#else///USE_HDF5_FORMAT
  iparticle body;
  ulong *idx;
  position *pos;
  acceleration *acc;
#ifdef  BLOCK_TIME_STEP
  velocity *vel;
  ibody_time *ti;
#else///BLOCK_TIME_STEP
  real *vx, *vy, *vz;
#endif//BLOCK_TIME_STEP
  allocParticleData((int)Ntot, &body, &idx, &pos, &acc,
#ifdef  BLOCK_TIME_STEP
		    &vel, &ti
#else///BLOCK_TIME_STEP
		    &vx, &vy, &vz
#endif//BLOCK_TIME_STEP
		    );
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const int nfile = (end - start + 1) / interval;
  //-----------------------------------------------------------------------
  PLFLT *time;  allocPLFLT(&time, nfile);
  PLFLT *step;  allocPLFLT(&step, nfile);
  //-----------------------------------------------------------------------
  PLFLT *Ekin;  allocPLFLT(&Ekin, nfile);
  PLFLT *Epot;  allocPLFLT(&Epot, nfile);
  PLFLT *Etot;  allocPLFLT(&Etot, nfile);
  //-----------------------------------------------------------------------
  PLFLT *momx;  allocPLFLT(&momx, nfile);
  PLFLT *momy;  allocPLFLT(&momy, nfile);
  PLFLT *momz;  allocPLFLT(&momz, nfile);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < nfile; ii++){
    //---------------------------------------------------------------------
    time[ii] = 0.0;
    step[ii] = 0.0;
    //---------------------------------------------------------------------
    Ekin[ii] = 0.0;
    Epot[ii] = 0.0;
    Etot[ii] = 0.0;
    //---------------------------------------------------------------------
    momx[ii] = 0.0;
    momy[ii] = 0.0;
    momz[ii] = 0.0;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < nfile; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  NORMALIZED_MOMENTUM_ERROR
  double MV = 0.0;
#endif//NORMALIZED_MOMENTUM_ERROR
  //-----------------------------------------------------------------------
  int ifile = mpi.rank;
  //-----------------------------------------------------------------------
  for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){
    //---------------------------------------------------------------------
    double tmp;
    ulong steps;
    int unit_read;
    readSnapshot(&unit_read, &tmp, &steps, Ntot, file, (uint)filenum
#ifdef  USE_HDF5_FORMAT
		 , &body, hdf5type
#else///USE_HDF5_FORMAT
		 , body
#endif//USE_HDF5_FORMAT
		 );
    if( unit_read != unit ){      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);    }
    //---------------------------------------------------------------------
    time[ifile] = (PLFLT)(tmp * time2astro);
    step[ifile] = (PLFLT)steps;
    //---------------------------------------------------------------------
#ifdef  NORMALIZED_MOMENTUM_ERROR
    double Mtot = 0.0;
#endif//NORMALIZED_MOMENTUM_ERROR
    //---------------------------------------------------------------------
    for(int ii = 0; ii < (int)Ntot; ii++){
      //-------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
      const double mass = (double)body.m  [ii        ];
      const double velx = (double)body.vel[ii * 3    ];
      const double vely = (double)body.vel[ii * 3 + 1];
      const double velz = (double)body.vel[ii * 3 + 2];
#else///USE_HDF5_FORMAT
      const double mass = (double)body.pos[ii].m;
#ifdef  BLOCK_TIME_STEP
      const double velx = (double)body.vel[ii].x;
      const double vely = (double)body.vel[ii].y;
      const double velz = (double)body.vel[ii].z;
#else///BLOCK_TIME_STEP
      const double velx = (double)body.vx[ii];
      const double vely = (double)body.vy[ii];
      const double velz = (double)body.vz[ii];
#endif//BLOCK_TIME_STEP
#endif//USE_HDF5_FORMAT
      //-------------------------------------------------------------------
#ifdef  NORMALIZED_MOMENTUM_ERROR
      Mtot += mass;
#endif//NORMALIZED_MOMENTUM_ERROR
      //-------------------------------------------------------------------
      Ekin[ifile] += (PLFLT)(mass * (velx * velx + vely * vely + velz * velz));
#ifdef  USE_HDF5_FORMAT
      Epot[ifile] += (PLFLT)(mass * (double)body.pot[ii]);
#else///USE_HDF5_FORMAT
      Epot[ifile] += (PLFLT)(mass * (double)body.acc[ii].pot);
#endif//USE_HDF5_FORMAT
      //-------------------------------------------------------------------
      momx[ifile] += mass * velx;
      momy[ifile] += mass * vely;
      momz[ifile] += mass * velz;
      //-------------------------------------------------------------------
/* #ifdef  NORMALIZED_MOMENTUM_ERROR */
/*       if( (mpi.rank == 0) && (ifile == 0) ) */
/* 	MV += mass * mass * (velx * velx + vely * vely + velz * velz); */
/* #endif//NORMALIZED_MOMENTUM_ERROR */
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
    //---------------------------------------------------------------------
    Ekin[ifile] *= 0.5 * energy2astro;
    Epot[ifile] *= 0.5 * energy2astro;
    Etot[ifile] = Ekin[ifile] + Epot[ifile];
    //---------------------------------------------------------------------
#ifdef  NORMALIZED_MOMENTUM_ERROR
    if( (mpi.rank == 0) && (ifile == 0) )
      MV = sqrt(-2.0 * Mtot * Etot[0]);
#endif//NORMALIZED_MOMENTUM_ERROR
    //---------------------------------------------------------------------
    momx[ifile] *= mass2astro * velocity2astro;
    momy[ifile] *= mass2astro * velocity2astro;
    momz[ifile] *= mass2astro * velocity2astro;
    //---------------------------------------------------------------------
    ifile += mpi.size;
    //---------------------------------------------------------------------
  }/* for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){ */
  //-----------------------------------------------------------------------
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, time, nfile, MPI_DOUBLE, MPI_SUM, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, step, nfile, MPI_DOUBLE, MPI_SUM, mpi.comm));
  //-----------------------------------------------------------------------
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, Ekin, nfile, MPI_DOUBLE, MPI_SUM, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, Epot, nfile, MPI_DOUBLE, MPI_SUM, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, Etot, nfile, MPI_DOUBLE, MPI_SUM, mpi.comm));
  //-----------------------------------------------------------------------
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, momx, nfile, MPI_DOUBLE, MPI_SUM, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, momy, nfile, MPI_DOUBLE, MPI_SUM, mpi.comm));
  chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, momz, nfile, MPI_DOUBLE, MPI_SUM, mpi.comm));
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  if( mpi.rank == 0 ){
    //---------------------------------------------------------------------
#ifdef  NORMALIZED_MOMENTUM_ERROR
    /* const double inv = 1.0 / (sqrt(MV) * mass2astro * velocity2astro); */
    const double inv = 1.0 / (DBL_MIN + MV * mass2astro * velocity2astro);
    for(int ii = 0; ii < nfile; ii++){
      momx[ii] *= inv;
      momy[ii] *= inv;
      momz[ii] *= inv;
    }/* for(int ii = 0; ii < nfile; ii++){ */
#endif//NORMALIZED_MOMENTUM_ERROR
    //---------------------------------------------------------------------
    printf("# time\tsteps\tEtot\tEkin\tEpot\tVirial ratio\n");
    for(int ii = 0; ii < nfile; ii++)
      printf("%e\t%e\t%e\t%e\t%e\t%e\n", time[ii], step[ii], Etot[ii], Ekin[ii], Epot[ii], Ekin[ii] / -Epot[ii]);
    //---------------------------------------------------------------------
    PLplotPltRange evolve, virial, relerr, momerr;
    //---------------------------------------------------------------------
    /* horizontal axis */
    evolve.xmin = time[0];  evolve.xmax = time[nfile - 1];  evolve.xlog = LINEAR_PLOT;
    virial.xmin = time[0];  virial.xmax = time[nfile - 1];  virial.xlog = LINEAR_PLOT;
    relerr.xmin = time[0];  relerr.xmax = time[nfile - 1];  relerr.xlog = LINEAR_PLOT;
    momerr.xmin = time[0];  momerr.xmax = time[nfile - 1];  momerr.xlog = LINEAR_PLOT;
    //---------------------------------------------------------------------
    /* vertical axis */
    switch( problem ){
    case  0:      /* cold collapse of a uniform sphere */
      evolve.ymin = -2.50  ;      evolve.ymax = 2.00  ;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.01  ;      virial.ymax = 0.99  ;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-2;      momerr.ymax = 1.0e-2;      momerr.ylog = LINEAR_PLOT;
      break;
    case  1:      /* king sphere with W0 = 3 */
      evolve.ymin = -0.6000;      evolve.ymax = 0.4000;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case  2:      /* Hernquist sphere with C = 10 */
      evolve.ymin = -0.4000;      evolve.ymax = 0.3000;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case  3:      /* NFW sphere with C = 10 */
      evolve.ymin = -0.2000;      evolve.ymax = 0.1500;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -3.0e-3;      momerr.ymax = 3.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case  4:      /* Einasto sphere with C = 10 */
      evolve.ymin = -0.4000;      evolve.ymax = 0.3000;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case  5:      /* Plummer sphere with C = 10 */
      evolve.ymin = -0.6000;      evolve.ymax = 0.4000;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case  6:      /* Burkert sphere with C = 10 */
      evolve.ymin = -0.3000;      evolve.ymax = 0.2000;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case  7:      /* Moore sphere with C = 10 */
      evolve.ymin = -0.4000;      evolve.ymax = 0.3000;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case  8:      /* Two-power sphere with C = 10 */
      evolve.ymin = -0.6000;      evolve.ymax = 0.4500;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case 10:      /* king sphere with W0 = 3 within an Einasto sphere with C = 10 */
      evolve.ymin = -1.2000;      evolve.ymax = 0.8000;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case 11:      /* A galaxy with multiple components */
      evolve.ymin = -0.0500;      evolve.ymax = 0.0400;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 12:      /* A galaxy with multiple components */
      evolve.ymin = -0.0500;      evolve.ymax = 0.0400;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 13:      /* A galaxy with multiple components */
      evolve.ymin = -3.9999;      evolve.ymax = 2.9999;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 20:      /* M31 model determined by Fardal et al. (2007) */
      evolve.ymin = -6.3e+6;      evolve.ymax = 4.2e+6;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4971;      virial.ymax = 0.5029;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 22:      /* A trial multi components galaxy model */
      evolve.ymin = -1.0e+7;      evolve.ymax = 6.0e+6;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4951;      virial.ymax = 0.5049;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 23:      /* MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008) */
      evolve.ymin = -8.0e+6;      evolve.ymax = 5.0e+6;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4951;      virial.ymax = 0.5049;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 24:      /* M31 model (Sofue 2015; Gilbert et al. 2012) */
      evolve.ymin = -1.5e+7;      evolve.ymax = 1.0e+7;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4971;      virial.ymax = 0.5029;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 25:      /* A trial multi components galaxy model */
      evolve.ymin = -1.0e+7;      evolve.ymax = 6.0e+6;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4951;      virial.ymax = 0.5049;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 26:      /* A trial multi components galaxy model (spherical model) */
      evolve.ymin = -1.0e+7;      evolve.ymax = 6.0e+6;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4951;      virial.ymax = 0.5049;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 27:      /* M31 model (NFW halo, de Vaucouleurs bulge, and exponential disk)*/
      evolve.ymin = -6.3e+6;      evolve.ymax = 4.2e+6;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4971;      virial.ymax = 0.5029;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 28:      /* A trial multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk) */
      evolve.ymin = -1.0e+7;      evolve.ymax = 6.0e+6;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4951;      virial.ymax = 0.5049;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 29:      /* Multi components galaxy model by Vasiliev & Athanassoula (2015) */
      evolve.ymin = -2.4e+1;      evolve.ymax = 1.4e+1;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4951;      virial.ymax = 0.5049;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case 30:      /* time evolution of MW/A defined in Kuijken & Dubinski (1995) */
      evolve.ymin = -1.6e+2;      evolve.ymax = 1.2e+2;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4871;      virial.ymax = 0.5049;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 31:      /* time evolution of MW/B defined in Kuijken & Dubinski (1995) */
      evolve.ymin = -3.2e+2;      evolve.ymax = 2.4e+2;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5039;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 32:      /* time evolution of MW/C defined in Kuijken & Dubinski (1995) */
      evolve.ymin = -7.5e+2;      evolve.ymax = 5.0e+2;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4941;      virial.ymax = 0.5019;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 33:      /* time evolution of MW/D defined in Kuijken & Dubinski (1995) */
      evolve.ymin = -1.6e+3;      evolve.ymax = 1.2e+3;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4951;      virial.ymax = 0.5029;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 34:      /* time evolution of M31/A defined in Widrow et al. (2003) */
      evolve.ymin = -4.0e+2;      evolve.ymax = 2.0e+2;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4801;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 35:      /* time evolution of M31/D defined in Widrow et al. (2003) */
      evolve.ymin = -4.0e+2;      evolve.ymax = 2.0e+2;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4551;      virial.ymax = 0.5349;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 36:      /* time evolution of MWa defined in Widrow & Dubinski (2005) */
      evolve.ymin = -2.0e+2;      evolve.ymax = 1.5e+2;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 37:      /* time evolution of MWb defined in Widrow & Dubinski (2005) */
      evolve.ymin = -2.0e+2;      evolve.ymax = 1.5e+2;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 38:      /* MW like galaxy (based on Widrow et al. 2008; Bedorf et al. 2014) */
      evolve.ymin = -2.0e+6;      evolve.ymax = 1.5e+6;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-3;      momerr.ymax = 1.0e-3;      momerr.ylog = LINEAR_PLOT;
      break;
    case 40:      /* Plummer sphere in table form with C = 20 */
      evolve.ymin = -9.0e+1;      evolve.ymax = 6.0e+1;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case 41:      /* Two-power sphere in table form with C = 20 */
      evolve.ymin = -1.6e+2;      evolve.ymax = 9.0e+1;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case 42:      /* de Vaucouleurs sphere in table form */
      evolve.ymin = -3.0e+4;      evolve.ymax = 2.0e+4;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case 50:      /* de Vaucouleurs sphere */
      evolve.ymin = -3.0e+4;      evolve.ymax = 2.0e+4;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case 51:      /* projected Two-power model */
      evolve.ymin = -4.0e+2;      evolve.ymax = 3.0e+2;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4901;      virial.ymax = 0.5099;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -5.0e-4;      momerr.ymax = 5.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    case 80:      /* A trial multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk) */
      evolve.ymin = -5.0e+6;      evolve.ymax = 3.0e+6;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.4951;      virial.ymax = 0.5049;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-4;      momerr.ymax = 1.0e-4;      momerr.ylog = LINEAR_PLOT;
      break;
    default:
      evolve.ymin = -2.50  ;      evolve.ymax = 2.00  ;      evolve.ylog = LINEAR_PLOT;
      virial.ymin =  0.01  ;      virial.ymax = 0.99  ;      virial.ylog = LINEAR_PLOT;
      momerr.ymin = -1.0e-2;      momerr.ymax = 1.0e-2;      momerr.ylog = LINEAR_PLOT;
      break;
    }
    relerr.ymin = (PLFLT)log10(1.0e-5);
    relerr.ymax = (PLFLT)log10(1.0e-1);
    relerr.ylog = LOGARITHMIC_PLOT;
    evolve.xgrd = virial.xgrd = relerr.xgrd = true;
    evolve.ygrd = virial.ygrd = relerr.ygrd = true;
    //---------------------------------------------------------------------
    evolve.ymin *= energy2astro;
    evolve.ymax *= energy2astro;
#ifndef NORMALIZED_MOMENTUM_ERROR
    momerr.ymin *= mass2astro * velocity2astro;
    momerr.ymax *= mass2astro * velocity2astro;
#endif//NORMALIZED_MOMENTUM_ERROR
    //---------------------------------------------------------------------
    real worst, error;
    plotTimeEvolution((PLINT)nfile, time, Ekin, Epot, Etot,
		      evolve, virial, relerr, file, argc, argv,
		      &worst, &error);
    //---------------------------------------------------------------------
    fprintf(stdout, "#eta\terror(t_fin)\terror(worst)\n");
    fprintf(stdout, "%f\t%e\t%e\n", eta, POW10(error), POW10(worst));
    //---------------------------------------------------------------------
    real momWorst[3], momError[3];
    plotMomentumError((PLINT)nfile, time, momx, momy, momz, momerr,
		      file, argc, argv, momWorst, momError);
    //---------------------------------------------------------------------
    fprintf(stdout, "#mom_dir\terror(t_fin)\terror(worst)\n");
#pragma unroll
    for(int ii = 0; ii < 3; ii++)
      fprintf(stdout, "%d\t%e\t%e\n", ii, momError[ii], momWorst[ii]);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
  freeSnapshotArray(hdf5_pos, hdf5_vel, hdf5_acc, hdf5_m, hdf5_pot, hdf5_idx);
#else///USE_HDF5_FORMAT
  freeParticleData(idx, pos, acc,
#ifdef  BLOCK_TIME_STEP
		    vel, ti
#else///BLOCK_TIME_STEP
		    vx, vy, vz
#endif//BLOCK_TIME_STEP
		    );
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  free(time);  free(step);
  free(Ekin);  free(Epot);  free(Etot);
  free(momx);  free(momy);  free(momz);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  exitMPI();
  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotTimeEvolution
(PLINT num, PLFLT *time, PLFLT *Ekin, PLFLT *Epot, PLFLT *Etot,
 PLplotPltRange evolve, PLplotPltRange virial, PLplotPltRange relerr,
 char *file, int argc, char **argv,
 real *worst, real *final)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* global setting(s) */
  //-----------------------------------------------------------------------
  /* maximum number of data kinds */
  const PLINT pkind = 3;
  const PLINT lkind = 3;
  //-----------------------------------------------------------------------
  const PLBOOL Euni = true;
  const PLBOOL Vuni = true;
  const PLBOOL Runi = true;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* preparation for dataset */
  //-----------------------------------------------------------------------
  /* memory allocation for index */
  PLINT Ekind = (Euni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  PLINT *Enum;  allocPLINT(&Enum, Ekind);
  PLINT Vkind = 1;  PLINT *Vnum;  allocPLINT(&Vnum, Vkind);
  PLINT Rkind = 1;  PLINT *Rnum;  allocPLINT(&Rnum, Rkind);
  //-----------------------------------------------------------------------
  /* set number of data points */
  for(PLINT ii = 0; ii < Ekind; ii++)    Enum[ii] = num;
  for(PLINT ii = 0; ii < Vkind; ii++)    Vnum[ii] = num;
  for(PLINT ii = 0; ii < Rkind; ii++)    Rnum[ii] = num;
  //-----------------------------------------------------------------------
  /* memory allocation for data */
  PLFLT **t;  allocPointer4PLFLT(&t, Ekind);
  PLFLT **E;  allocPointer4PLFLT(&E, Ekind);
  PLFLT **V;  allocPointer4PLFLT(&V, Vkind);
  PLFLT **R;  allocPointer4PLFLT(&R, Rkind);
  PLFLT *_V;
  PLFLT *_R;
  {
    PLINT tot = 0;
    for(PLINT ii = 0; ii < Vkind; ii++)      tot += Vnum[ii];
    allocPLFLT(&_V, tot);
    tot = 0;
    for(PLINT ii = 0; ii < Rkind; ii++)      tot += Rnum[ii];
    allocPLFLT(&_R, tot);
  }
  V[0] = _V;
  R[0] = _R;
  for(PLINT ii = 1; ii < Vkind; ii++)    V[ii] = V[ii - 1] + Vnum[ii - 1];
  for(PLINT ii = 1; ii < Rkind; ii++)    R[ii] = R[ii - 1] + Rnum[ii - 1];
  /* data preparation */
  t[0] = time;  E[0] = Etot;
  t[1] = time;  E[1] = Ekin;
  t[2] = time;  E[2] = Epot;
  for(PLINT ii = 0; ii < Vkind; ii++)
    for(PLINT jj = 0; jj < Vnum[ii]; jj++)
      V[ii][jj] = -Ekin[jj] / Epot[jj];
  for(PLINT ii = 0; ii < Rkind; ii++){
    PLFLT invEini = 1.0 / Etot[0];
    for(PLINT jj = 0; jj < Rnum[ii]; jj++)
      R[ii][jj] = log10((PLFLT)EPSILON + fabs(Etot[jj] * invEini - 1.0));
  }
  //-----------------------------------------------------------------------
  /* data analysis */
  *worst = LOG10(EPSILON);
  for(PLINT ii = 0; ii < Rkind; ii++)
    for(PLINT jj = 0; jj < Rnum[ii]; jj++)
      if( *worst < (real)R[ii][jj] )
	*worst = (real)R[ii][jj];
  *final = (real)R[0][num - 1];
  //-----------------------------------------------------------------------
  /* set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  //-----------------------------------------------------------------------
  /* set labels */
  char xlab[PLplotCharWords];  sprintf(xlab, "#fit #fr(%s)", time_astro_unit_name4plot);
  char Elab[PLplotCharWords];  sprintf(Elab, "Energy (%s)", energy_astro_unit_name4plot);
  char Vlab[PLplotCharWords];  sprintf(Vlab, "Virial ratio |#fiK#fr/#fiW#fr|");
  char Rlab[PLplotCharWords];  sprintf(Rlab, "relative error of #fiE#fr#dtot#u");
  //-----------------------------------------------------------------------
  /* set caption(s) */
  PLplotCaption Ecap;  setDefaultCaption(&Ecap);  Ecap.write = false;
  PLplotCaption Vcap;  setDefaultCaption(&Vcap);  Vcap.write = false;
  PLplotCaption Rcap;  setDefaultCaption(&Rcap);  Rcap.write = false;
  sprintf(Ecap.side, "%s", "t");
  sprintf(Vcap.side, "%s", "t");
  sprintf(Rcap.side, "%s", "t");
  sprintf(Ecap.text, "%s", "Time evolution of energy");
  sprintf(Vcap.text, "%s", "Time evolution of virial ratio");
  sprintf(Rcap.text, "%s", "Time evolution of relative error");
  //-----------------------------------------------------------------------
  /* set legends */
  PLplotLegend Eleg;  setDefaultLegend(&Eleg, false);  Eleg.write = false;
  PLplotLegend Vleg;  setDefaultLegend(&Vleg, false);  Vleg.write = false;
  PLplotLegend Rleg;  setDefaultLegend(&Rleg, false);  Rleg.write = false;
  char *ElegTxt;
  {
    allocChar4PLplot(&ElegTxt, Ekind);
    allocPointer4Char4PLplot(&(Eleg.text), Ekind);
    assignChar4PLplot(Ekind, Eleg.text, ElegTxt);
  }
  sprintf(Eleg.text[0], "%s",     "total energy");
  sprintf(Eleg.text[1], "%s",   "kinetic energy");
  sprintf(Eleg.text[2], "%s", "potential energy");
  char *VlegTxt;
  {
    allocChar4PLplot(&VlegTxt, Vkind);
    allocPointer4Char4PLplot(&(Vleg.text), Vkind);
    assignChar4PLplot(Vkind, Vleg.text, VlegTxt);
  }
  sprintf(Vleg.text[0], "%s", "virial ratio");
  char *RlegTxt;
  {
    allocChar4PLplot(&RlegTxt, Rkind);
    allocPointer4Char4PLplot(&(Rleg.text), Rkind);
    assignChar4PLplot(Rkind, Rleg.text, RlegTxt);
  }
  sprintf(Rleg.text[0], "%s", "relative error");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation to exploit multiple-plot function */
  //-----------------------------------------------------------------------
  /* maximum number of panels in each direction */
  const PLINT nxpanel = 1;
  const PLINT nypanel = 3;
  //-----------------------------------------------------------------------
  /* specify plot or skip panel */
  PLBOOL *plot;  allocPLBOOL(&plot, nxpanel * nypanel);
  /* arrays related to draw line(s) and point(s) */
  PLINT *nlkind;  allocPLINT(&nlkind, nxpanel * nypanel);
  PLINT *npkind;  allocPLINT(&npkind, nxpanel * nypanel);
  PLINT **lnum;  allocPointer4PLINT(&lnum, nxpanel * nypanel);
  PLINT **pnum;  allocPointer4PLINT(&pnum, nxpanel * nypanel);
  PLFLT ***lx;  allocDoublePointer4PLFLT(&lx, nxpanel * nypanel);
  PLFLT ***ly;  allocDoublePointer4PLFLT(&ly, nxpanel * nypanel);
  PLFLT ***px;  allocDoublePointer4PLFLT(&px, nxpanel * nypanel);
  PLFLT ***py;  allocDoublePointer4PLFLT(&py, nxpanel * nypanel);
  PLplotLineStyle ** line;  allocPointer4PLplotLineStyle(& line, nxpanel * nypanel);
  PLplotPointType **point;  allocPointer4PLplotPointType(&point, nxpanel * nypanel);
  /* arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /* arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /* arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /* arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /* arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /* array to set figure name */
  char figfile[PLplotCharWords];
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel * nypanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel * nypanel);
  assignChar4PLplot(nxpanel * nypanel, figname, _figname);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* configure to time evolution of energy */
  sprintf(figfile, "%s_%s", file, "energy");
  //-----------------------------------------------------------------------
  /* common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    for(PLINT jj = 0; jj < nypanel; jj++){
      //-------------------------------------------------------------------
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      //-------------------------------------------------------------------
      /* global setting(s) */
      plot[idx] = true;
      //-------------------------------------------------------------------
      /* line setting(s) */
      line[idx] = ls;
      lx  [idx] = t;
      //-------------------------------------------------------------------
      /* point setting(s) */
      point[idx] = pt;
      px   [idx] = t;
      //-------------------------------------------------------------------
      /* errorbar setting(s) */
      errbar[idx] = 0;
      //-------------------------------------------------------------------
      /* label setting(s) */
      xlabel[idx] = xlab;
      //-------------------------------------------------------------------
    }
  }
  //-----------------------------------------------------------------------
  /* setting(s) for time evolution of energy */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 0);
    //---------------------------------------------------------------------
    /* line setting(s) */
    nlkind[idx] = lkind;
    lnum  [idx] = Enum;
    ly    [idx] = E;
    //---------------------------------------------------------------------
    /* point setting(s) */
    npkind[idx] = pkind;
    pnum  [idx] = Enum;
    py    [idx] = E;
    //---------------------------------------------------------------------
    /* plot area */
    range[idx] = evolve;
    //---------------------------------------------------------------------
    /* label setting(s) */
    ylabel[idx] = Elab;
    //---------------------------------------------------------------------
    /* miscs */
    cap[idx] = Ecap;
    leg[idx] = Eleg;
    uni[idx] = Euni;
    //---------------------------------------------------------------------
    /* figure name */
    sprintf(figname[idx], "%s_%s", figfile, "evolve");
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  /* setting(s) for time evolution of virial ratio */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 1);
    //---------------------------------------------------------------------
    /* line setting(s) */
    nlkind[idx] = Vkind;
    lnum  [idx] = Vnum;
    ly    [idx] = V;
    //---------------------------------------------------------------------
    /* point setting(s) */
    npkind[idx] = Vkind;
    pnum  [idx] = Vnum;
    py    [idx] = V;
    //---------------------------------------------------------------------
    /* plot area */
    range[idx] = virial;
    //---------------------------------------------------------------------
    /* label setting(s) */
    ylabel[idx] = Vlab;
    //---------------------------------------------------------------------
    /* miscs */
    cap[idx] = Vcap;
    leg[idx] = Vleg;
    uni[idx] = Vuni;
    //---------------------------------------------------------------------
    /* figure name */
    sprintf(figname[idx], "%s_%s", figfile, "virial");
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  /* setting(s) for time evolution of relative error of the total energy */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 2);
    //---------------------------------------------------------------------
    /* line setting(s) */
    nlkind[idx] = Rkind;
    lnum  [idx] = Rnum;
    ly    [idx] = R;
    //---------------------------------------------------------------------
    /* point setting(s) */
    npkind[idx] = Rkind;
    pnum  [idx] = Rnum;
    py    [idx] = R;
    //---------------------------------------------------------------------
    /* plot area */
    range[idx] = relerr;
    //---------------------------------------------------------------------
    /* label setting(s) */
    ylabel[idx] = Rlab;
    //---------------------------------------------------------------------
    /* miscs */
    cap[idx] = Rcap;
    leg[idx] = Rleg;
    uni[idx] = Runi;
    //---------------------------------------------------------------------
    /* figure name */
    sprintf(figname[idx], "%s_%s", figfile, "relerr");
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create figure(s) */
  //-----------------------------------------------------------------------
  plotData(nxpanel, nypanel, plot, true, true,
  	   nlkind,  line, lnum, lx, ly,
  	   npkind, point, pnum, px, py,
  	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
  	   cap, leg, uni, range,
  	   xlabel, ylabel, "", figfile, argc, argv);
  //-----------------------------------------------------------------------
  for(PLINT idx = 0; idx < nxpanel * nypanel; idx++)
    plotData(1, 1, &plot[idx], false, false,
  	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
  	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
  	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
  	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
  	     &xlabel[idx], &ylabel[idx], "", figname[idx], argc, argv);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  free(figname);  free(_figname);
  //-----------------------------------------------------------------------
  free(Enum);  free(t);  free( E);
  free(Vnum);  free(V);  free(_V);
  free(Rnum);  free(R);  free(_R);
  free(ls);  free(pt);
  //-----------------------------------------------------------------------
  free(ElegTxt);
  free(VlegTxt);
  free(RlegTxt);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotMomentumError
(PLINT num, PLFLT *time, PLFLT *momx, PLFLT *momy, PLFLT *momz, PLplotPltRange momerr,
 char *file, int argc, char **argv,
 real worst[], real final[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* global setting(s) */
  //-----------------------------------------------------------------------
  /* maximum number of data kinds */
  const PLINT pkind = 3;
  const PLINT lkind = 3;
  //-----------------------------------------------------------------------
  const PLBOOL Euni = true;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* preparation for dataset */
  //-----------------------------------------------------------------------
  /* memory allocation for index */
  PLINT Ekind = (Euni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  PLINT *Enum;  allocPLINT(&Enum, Ekind);
  //-----------------------------------------------------------------------
  /* set number of data points */
  for(PLINT ii = 0; ii < Ekind; ii++)    Enum[ii] = num;
  //-----------------------------------------------------------------------
  /* memory allocation for data */
  PLFLT **t;  allocPointer4PLFLT(&t, Ekind);
  PLFLT **E;  allocPointer4PLFLT(&E, Ekind);
  /* data preparation and analysis */
  t[0] = time;  E[0] = momx;
  t[1] = time;  E[1] = momy;
  t[2] = time;  E[2] = momz;
  worst[0] = worst[1] = worst[2] = ZERO;
  for(PLINT ii = 0; ii < Ekind; ii++){
    //---------------------------------------------------------------------
    for(PLINT jj = 0; jj < Enum[ii]; jj++){
      //-------------------------------------------------------------------
      if( fabs(worst[ii]) < (real)fabs(E[ii][jj] - E[ii][0]))
  	worst[ii] = (real)(E[ii][jj] - E[ii][0]);
      //-------------------------------------------------------------------
    }
    //---------------------------------------------------------------------
    final[ii] = (real)(E[ii][Enum[ii] - 1] - E[ii][0]);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  /* set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  //-----------------------------------------------------------------------
  /* set labels */
  char xlab[PLplotCharWords];  sprintf(xlab, "#fit #fr(%s)", time_astro_unit_name4plot);
#ifdef  NORMALIZED_MOMENTUM_ERROR
  /* char Elab[PLplotCharWords];  sprintf(Elab, "#gS#di#u #fim#fr#di#u#fi#<0x12>v#fr#di#u / #gS#di#u |#fim#fr#di#u#fi#<0x12>v#fr#di#u|(#fit#fr=0)"); */
  char Elab[PLplotCharWords];  sprintf(Elab, "#gS#di#u #fim#fr#di#u#fi#<0x12>v#fr#di#u / (-2 #fiM#fr#dtot#u #fiE#fr#dtot#u)#u1/2#d (#fit#fr=0)");
#else///NORMALIZED_MOMENTUM_ERROR
  char Elab[PLplotCharWords];  sprintf(Elab, "Momentum (%s %s)", mass_astro_unit_name4plot, velocity_astro_unit_name4plot);
#endif//NORMALIZED_MOMENTUM_ERROR
  //-----------------------------------------------------------------------
  /* set caption(s) */
  PLplotCaption Ecap;  setDefaultCaption(&Ecap);  Ecap.write = false;
  sprintf(Ecap.side, "%s", "t");
  sprintf(Ecap.text, "%s", "Time evolution of momentum");
  //-----------------------------------------------------------------------
  /* set legends */
  PLplotLegend Eleg;  setDefaultLegend(&Eleg, false);  Eleg.write = true;
  char *ElegTxt;
  {
    allocChar4PLplot(&ElegTxt, Ekind);
    allocPointer4Char4PLplot(&(Eleg.text), Ekind);
    assignChar4PLplot(Ekind, Eleg.text, ElegTxt);
  }
  sprintf(Eleg.text[0], "%s", "#fip#fr#dx#u");
  sprintf(Eleg.text[1], "%s", "#fip#fr#dy#u");
  sprintf(Eleg.text[2], "%s", "#fip#fr#dz#u");
  Eleg.pos = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation to exploit multiple-plot function */
  //-----------------------------------------------------------------------
  /* maximum number of panels in each direction */
  const PLINT nxpanel = 1;
  const PLINT nypanel = 1;
  //-----------------------------------------------------------------------
  /* specify plot or skip panel */
  PLBOOL *plot;  allocPLBOOL(&plot, nxpanel * nypanel);
  /* arrays related to draw line(s) and point(s) */
  PLINT *nlkind;  allocPLINT(&nlkind, nxpanel * nypanel);
  PLINT *npkind;  allocPLINT(&npkind, nxpanel * nypanel);
  PLINT **lnum;  allocPointer4PLINT(&lnum, nxpanel * nypanel);
  PLINT **pnum;  allocPointer4PLINT(&pnum, nxpanel * nypanel);
  PLFLT ***lx;  allocDoublePointer4PLFLT(&lx, nxpanel * nypanel);
  PLFLT ***ly;  allocDoublePointer4PLFLT(&ly, nxpanel * nypanel);
  PLFLT ***px;  allocDoublePointer4PLFLT(&px, nxpanel * nypanel);
  PLFLT ***py;  allocDoublePointer4PLFLT(&py, nxpanel * nypanel);
  PLplotLineStyle ** line;  allocPointer4PLplotLineStyle(& line, nxpanel * nypanel);
  PLplotPointType **point;  allocPointer4PLplotPointType(&point, nxpanel * nypanel);
  /* arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /* arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /* arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /* arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /* arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /* array to set figure name */
  char figfile[PLplotCharWords];
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel * nypanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel * nypanel);
  assignChar4PLplot(nxpanel * nypanel, figname, _figname);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* configure to time evolution of momentum */
  sprintf(figfile, "%s_%s", file, "momentum");
  //-----------------------------------------------------------------------
  /* common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    for(PLINT jj = 0; jj < nypanel; jj++){
      //-------------------------------------------------------------------
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      //-------------------------------------------------------------------
      /* global setting(s) */
      plot[idx] = true;
      //-------------------------------------------------------------------
      /* line setting(s) */
      nlkind[idx] = lkind;
      line  [idx] = ls;
      lnum  [idx] = Enum;
      lx    [idx] = t;
      ly    [idx] = E;
      //-------------------------------------------------------------------
      /* point setting(s) */
      npkind[idx] = pkind;
      point [idx] = pt;
      pnum  [idx] = Enum;
      py    [idx] = E;
      px    [idx] = t;
      //-------------------------------------------------------------------
      /* errorbar setting(s) */
      errbar[idx] = 0;
      //-------------------------------------------------------------------
      /* plot area */
      range[idx] = momerr;
      //-------------------------------------------------------------------
      /* label setting(s) */
      xlabel[idx] = xlab;
      ylabel[idx] = Elab;
      //-------------------------------------------------------------------
      /* miscs */
      cap[idx] = Ecap;
      leg[idx] = Eleg;
      uni[idx] = Euni;
      //-------------------------------------------------------------------
      /* figure name */
      sprintf(figname[idx], "%s_%s", figfile, "err");
      //-------------------------------------------------------------------
    }
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* create figure(s) */
  //-----------------------------------------------------------------------
  for(PLINT idx = 0; idx < nxpanel * nypanel; idx++)
    plotData(1, 1, &plot[idx], false, false,
  	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
  	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
  	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
  	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
  	     &xlabel[idx], &ylabel[idx], "", figname[idx], argc, argv);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  free(figname);  free(_figname);
  //-----------------------------------------------------------------------
  free(Enum);  free(t);  free( E);
  free(ls);  free(pt);
  //-----------------------------------------------------------------------
  free(ElegTxt);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
