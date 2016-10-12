/*************************************************************************\
 *                                                                       *
                  last updated on 2016/10/11(Tue) 17:06:04
 *                                                                       *
 *    Plot Code of N-body Simulations (using PLplot)                     *
 *      Time Evolution of Spatial Distribution Maps                      *
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
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <myutil.h>
#include <name.h>
#include <constants.h>
#include <mpilib.h>
#include <plplotlib.h>
//-------------------------------------------------------------------------
#include "../misc/structure.h"
#include "../misc/allocate.h"
//-------------------------------------------------------------------------
#include "../file/io.h"
//-------------------------------------------------------------------------
#define NRADBIN (128)
//-------------------------------------------------------------------------
#define NMAX_GAUSS_QD (51)
#define NTBL_GAUSS_QD ((NMAX_GAUSS_QD >> 1) + (NMAX_GAUSS_QD & 1))
#define NINTBIN NMAX_GAUSS_QD
real gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
//-------------------------------------------------------------------------
/* #define PLOT_SPLITTED_DISTRIBUTION_MAPS */
//-------------------------------------------------------------------------
#define REDUCE_PARTICLE_DISTRIBUTION_MAP (8192)
//-------------------------------------------------------------------------
extern const double length2astro, time2astro, density2astro, mass2astro;
//-------------------------------------------------------------------------
static int ncrit;
//-------------------------------------------------------------------------
static const int NAllocUnit = 32;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* (l + 1)^2 components are required: */
/* l <= 1 -->  4 components */
/* l <= 2 -->  9 components */
/* l <= 3 --> 16 components */
/* l <= 4 --> 25 components */
/* #define NUM_SPH_HARMONICS (4) */
/* #define NUM_SPH_HARMONICS (9) */
/* #define NUM_SPH_HARMONICS (16) */
#define NUM_SPH_HARMONICS (25)
//-------------------------------------------------------------------------
typedef struct
{
  real Plm[NINTBIN];/* normalized Plm(x) */
  gsl_complex zz[NINTBIN];
  int l, m;
} Ylm;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* /\* this value must be the same with NDISKBIN in ../init/disk.h *\/ */
/* #define NUM_ANALYTIC (256) */
/* #define NUM_HORIZONTAL_BIN (10) */
//-------------------------------------------------------------------------
typedef struct
{
  int prf_head[NUM_SPH_HARMONICS], prf_num[NUM_SPH_HARMONICS];
  ulong head, num;/* data for N-body particles */
} model;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotRadialProfile
(const int ntot, const int ngroup, model *group, Ylm sph_harm[], PLFLT time,
 real *rad, real *rho, PLplotPltRange rhobox,
 char *file, const int filenum, int argc, char **argv);
//-------------------------------------------------------------------------
void analyzeRadialProfile
(int Np, nbody_particle *body_tot, const int kind, model *group, int *num, int *rem, const real r2max, real **rad, real **rho, real weight[], Ylm sph_harm[],
 real rad_tbl[], real costhetamin[], real costhetabin[], real phimin[], real phibin[], real rho_tbl[][NINTBIN], real rholm_tbl[][NRADBIN]);
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
int idxAscendingOrder(const void *a, const void *b)
{
  //-----------------------------------------------------------------------
  if(          ((nbody_particle *)a)->idx > ((nbody_particle *)b)->idx ){    return ( 1);  }
  else{    if( ((nbody_particle *)a)->idx < ((nbody_particle *)b)->idx ){    return (-1);  }
    else{                                                                    return ( 0);  }  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
int radAscendingOrder(const void *a, const void *b)
{
  //-----------------------------------------------------------------------
  if(          ((nbody_particle *)a)->ax > ((nbody_particle *)b)->ax ){    return ( 1);  }
  else{    if( ((nbody_particle *)a)->ax < ((nbody_particle *)b)->ax ){    return (-1);  }
    else{                                                                  return ( 0);  }  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void allocProfileArray(int num, real **rad, real **rho)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  *rad = (real *)malloc(sizeof(real) * num);
  *rho = (real *)malloc(sizeof(real) * num);
  if( (*rad == NULL) || (*rho == NULL) ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void enlargeProfileArray(int num, real **rad, real **rho)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  *rad = realloc(*rad, sizeof(real) * num);
  *rho = realloc(*rho, sizeof(real) * num);
  if( (*rad == NULL) || (*rho == NULL) ){
    __KILL__(stderr, "%s\n", "ERROR: memory allocation failed.");
  }
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void initGridSpacing(real *weight, real *xval, real *xm, real *xp, real *dx, real *phival, real *phim, real *phip, real *dphi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const real xmin = -UNITY;  /* const real phimin = -(real)M_PI; */
  const real xmax =  UNITY;  /* const real phimax =  (real)M_PI; */
  //-----------------------------------------------------------------------
  const real xmns = HALF * (xmax - xmin);  /* const real phimns = HALF * (phimax - phimin); */
  const real xpls = HALF * (xmax + xmin);  /* const real phipls = HALF * (phimax + phimin); */
  //-----------------------------------------------------------------------
  const int imid = (NINTBIN >> 1);
  if( NINTBIN & 1 ){
    weight[imid] =               gsl_gaussQD_weight[imid];
    xval  [imid] = xpls + xmns * gsl_gaussQD_pos   [imid];
  }
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < imid; ii++){
    const real tmp         =               gsl_gaussQD_weight[ii];
    xval[              ii] = xpls - xmns * gsl_gaussQD_pos   [ii];    weight[              ii] = tmp;
    xval[NINTBIN - 1 - ii] = xpls + xmns * gsl_gaussQD_pos   [ii];    weight[NINTBIN - 1 - ii] = tmp;
  }
  //-----------------------------------------------------------------------
  xm[          0] = xmin;
  for(int ii = 1; ii < NINTBIN; ii++)
    xm[ii] = HALF * (xval[ii - 1] + xval[ii]);
  xp[NINTBIN - 1] = xmax;
  for(int ii = NINTBIN - 2; ii >= 0; ii--)
    xp[ii] = HALF * (xval[ii + 1] + xval[ii]);
  dx[0] = xp[0] - xm[0];
  for(int ii = 1; ii < NINTBIN - 1; ii++)
    dx[ii] = HALF * (xval[ii + 1] - xval[ii - 1]);
  dx[NINTBIN - 1] = xp[NINTBIN - 1] - xm[NINTBIN - 1];
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < NINTBIN; ii++){
    phival[ii] = (real)M_PI * xval[ii];
    phim[ii] = (real)M_PI * xm[ii];
    phip[ii] = (real)M_PI * xp[ii];
    dphi[ii] = (real)M_PI * dx[ii];
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
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
  /* setPhysicalConstantsAndUnitSystem(UNITSYSTEM, 0); */
  //-----------------------------------------------------------------------
  if( argc < 7 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 7);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -problem=<int>\n");
    __FPRINTF__(stderr, "          -ncrit=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }
  //-----------------------------------------------------------------------
  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));
  int  problem;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,  "problem", &problem));
  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "ncrit", &ncrit));
  //-----------------------------------------------------------------------
  modifyArgcArgv4PLplot(&argc, argv, 7);
  //-----------------------------------------------------------------------
  /* initialize the table for Gaussian Quadrature provided by GSL */
  for(int ii = 0; ii < NTBL_GAUSS_QD; ii++){
    gsl_gaussQD_pos   [ii] = ZERO;
    gsl_gaussQD_weight[ii] = ZERO;
  }
  gsl_integration_glfixed_table *tab;
  tab = gsl_integration_glfixed_table_alloc(NINTBIN);
  int max = (NINTBIN >> 1) + (NINTBIN & 1);
  for(int ii = 0; ii < max; ii++){
    gsl_gaussQD_pos   [ii] = (real)(*tab).x[(max - 1) - ii];
    gsl_gaussQD_weight[ii] = (real)(*tab).w[(max - 1) - ii];
  }
  gsl_integration_glfixed_table_free(tab);
  //-----------------------------------------------------------------------
  static real   weight[NINTBIN];
  static real costheta[NINTBIN], costhetamin[NINTBIN], costhetamax[NINTBIN], costhetabin[NINTBIN];
  static real      phi[NINTBIN],      phimin[NINTBIN],      phimax[NINTBIN],      phibin[NINTBIN];
  initGridSpacing(weight, costheta, costhetamin, costhetamax, costhetabin, phi, phimin, phimax, phibin);
#if 0
  for(int ii = 0; ii < NINTBIN; ii++)
    printf("%d\t%e\t%e\t%e\t%e\t%e\n", ii, weight[ii], costheta[ii], costhetamin[ii], phi[ii], phimin[ii]);
  exit(0);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* prepare spherical harmonics */
  //-----------------------------------------------------------------------
  static Ylm sph_harm[NUM_SPH_HARMONICS];
  for(int ii = 0; ii < NUM_SPH_HARMONICS; ii++){
    //---------------------------------------------------------------------
    /* set l and m */
    //---------------------------------------------------------------------
    switch( ii ){
      /* monopole */
    case  0:      sph_harm[ii].l = 0;      sph_harm[ii].m =  0;      break;
      /* dipole */
    case  1:      sph_harm[ii].l = 1;      sph_harm[ii].m = -1;      break;
    case  2:      sph_harm[ii].l = 1;      sph_harm[ii].m =  0;      break;
    case  3:      sph_harm[ii].l = 1;      sph_harm[ii].m =  1;      break;
      /* quadrupole */
    case  4:      sph_harm[ii].l = 2;      sph_harm[ii].m = -2;      break;
    case  5:      sph_harm[ii].l = 2;      sph_harm[ii].m = -1;      break;
    case  6:      sph_harm[ii].l = 2;      sph_harm[ii].m =  0;      break;
    case  7:      sph_harm[ii].l = 2;      sph_harm[ii].m =  1;      break;
    case  8:      sph_harm[ii].l = 2;      sph_harm[ii].m =  2;      break;
      /* octupole */
    case  9:      sph_harm[ii].l = 3;      sph_harm[ii].m = -3;      break;
    case 10:      sph_harm[ii].l = 3;      sph_harm[ii].m = -2;      break;
    case 11:      sph_harm[ii].l = 3;      sph_harm[ii].m = -1;      break;
    case 12:      sph_harm[ii].l = 3;      sph_harm[ii].m =  0;      break;
    case 13:      sph_harm[ii].l = 3;      sph_harm[ii].m =  1;      break;
    case 14:      sph_harm[ii].l = 3;      sph_harm[ii].m =  2;      break;
    case 15:      sph_harm[ii].l = 3;      sph_harm[ii].m =  3;      break;
      /* hexadecapole */
    case 16:      sph_harm[ii].l = 4;      sph_harm[ii].m = -4;      break;
    case 17:      sph_harm[ii].l = 4;      sph_harm[ii].m = -3;      break;
    case 18:      sph_harm[ii].l = 4;      sph_harm[ii].m = -2;      break;
    case 19:      sph_harm[ii].l = 4;      sph_harm[ii].m = -1;      break;
    case 20:      sph_harm[ii].l = 4;      sph_harm[ii].m =  0;      break;
    case 21:      sph_harm[ii].l = 4;      sph_harm[ii].m =  1;      break;
    case 22:      sph_harm[ii].l = 4;      sph_harm[ii].m =  2;      break;
    case 23:      sph_harm[ii].l = 4;      sph_harm[ii].m =  3;      break;
    case 24:      sph_harm[ii].l = 4;      sph_harm[ii].m =  4;      break;
    default:
      __KILL__(stderr, "ERROR: currently, l > 4 is not supported\n");
      break;
    }
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* calculate sqrt((2l+1)/4pi) * sqrt((l-m)! / (l+m)!) Plm(x) */
    //---------------------------------------------------------------------
    for(int jj = 0; jj < NINTBIN; jj++)
      sph_harm[ii].Plm[jj] = (real)gsl_sf_legendre_sphPlm(sph_harm[ii].l, abs(sph_harm[ii].m), (double)costheta[jj]);
    for(int jj = 0; jj < NINTBIN; jj++)
      sph_harm[ii].zz[jj] = gsl_complex_polar(1.0, (double)(sph_harm[ii].m * phi[jj]));
    if( sph_harm[ii].m < 0 ){
      int mm = abs(sph_harm[ii].m);
      real sign = (mm & 1) ? (-UNITY) : (UNITY);
      for(int jj = 0; jj < NINTBIN; jj++)
	sph_harm[ii].Plm[jj] *= sign;
    }
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NUM_SPH_HARMONICS; ii++){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* load global settings of particle distribution */
  //-----------------------------------------------------------------------
  int last;
  readConfigFile(&last, file);
  int unit;
  ulong Ntot;
  real eta, eps;
  double saveInterval, ft, snapshotInterval;
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file, mpi);
  nbody_particle *body;
  allocParticleDataAoS((int)Ntot, &body);
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
  nbody_hdf5 hdf5;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
  ulong *hdf5_idx;
  allocSnapshotArray(&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx, (int)Ntot, &hdf5);
#endif//USE_HDF5_FORMAT
  //-----------------------------------------------------------------------
  int kind = 1;
  int skind = 1;
  model *group;
  bool multi_group = false;
  if( problem >= 4 ){
    //---------------------------------------------------------------------
    multi_group = true;
    //---------------------------------------------------------------------
    FILE *fp;
    char filename[256];
    sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
    fp = fopen(filename, "r");
    if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
    //---------------------------------------------------------------------
    fscanf(fp, "%*d");/* skip reading unit */
    fscanf(fp, "%d\t%d", &kind, &skind);
    group = (model *)malloc(kind * sizeof(model));
    if( group == NULL ){      __KILL__(stderr, "ERROR: failure to allocate group\n");    }
    for(int ii = 0; ii < kind; ii++)
      fscanf(fp, "%zu", &group[ii].num);
    //---------------------------------------------------------------------
    fclose(fp);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  if( !multi_group ){
    //---------------------------------------------------------------------
    group = (model *)malloc(kind * sizeof(model));
    if( group == NULL ){      __KILL__(stderr, "ERROR: failure to allocate group\n");    }
    //---------------------------------------------------------------------
    group[0].num = Ntot;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set plot range */
  //-----------------------------------------------------------------------
  PLFLT radius;
  double radmin, radmax;
  switch( problem ){
  case  0:    /* cold collapse of a uniform sphere */
    radius = 3.0e+0;    radmin = 1.0e-2;    radmax = 2.0e+1;    break;
  case  1:    /* king sphere with W0 = 3 */
    radius = 5.0e+0;    radmin = 1.0e-2;    radmax = 1.0e+1;    break;
  case  2:    /* Hernquist sphere with C = 10 */
    radius = 1.0e+1;    radmin = 1.0e-2;    radmax = 2.0e+1;    break;
  case  3:    /* NFW sphere with C = 5 */
    radius = 5.0e+0;    radmin = 1.0e-2;    radmax = 1.0e+1;    break;
  case 11:    /* A galaxy with multiple components */
    radius = 5.0e+1;    radmin = 1.0e-2;    radmax = 4.0e+2;    break;
  default:
    radius = 3.0e+0;    radmin = 1.0e-2;    radmax = 2.0e+1;    break;
  }
  //-----------------------------------------------------------------------
  PLplotPltRange rhorange;
  rhorange.xmin = (PLFLT)log10(radmin * (double) length2astro);  rhorange.ymin = 0.0;
  rhorange.xmax = (PLFLT)log10(radmax * (double) length2astro);  rhorange.ymax = 1.0;
  rhorange.xlog = LOGARITHMIC_PLOT;  rhorange.xgrd = true;
  rhorange.ylog =      LINEAR_PLOT;  rhorange.ygrd = true;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* temporary arrays */
  static real   rad_tbl[NRADBIN];
  static real   rho_tbl[NINTBIN][NINTBIN];
  static real rholm_tbl[NUM_SPH_HARMONICS][NRADBIN];
  //-----------------------------------------------------------------------
  /* arrays to be plotted */
  real *rad, *rholm;
  int num = 0;
  int rem = NAllocUnit;
  allocProfileArray(rem, &rad, &rholm);
  //-----------------------------------------------------------------------
  int ifile = (int)(start + mpi.rank);
  for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){
    //---------------------------------------------------------------------
    double time;
    ulong steps;
    int unit_read;
#ifdef  USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, &hdf5, file, filenum, hdf5type);
    for(int ii = 0; ii < (int)Ntot; ii++){
      //-------------------------------------------------------------------
      body[ii]. x  = hdf5.pos[ii * 3];      body[ii]. y = hdf5.pos[ii * 3 + 1];      body[ii].z   = hdf5.pos[ii * 3 + 2];
      body[ii].vx  = hdf5.vel[ii * 3];      body[ii].vy = hdf5.vel[ii * 3 + 1];      body[ii].vz  = hdf5.vel[ii * 3 + 2];
      body[ii].ax  = hdf5.acc[ii * 3];      body[ii].ay = hdf5.acc[ii * 3 + 1];      body[ii].az  = hdf5.acc[ii * 3 + 2];
      body[ii].idx = hdf5.idx[ii    ];      body[ii]. m = hdf5.  m[ii        ];      body[ii].pot = hdf5.pot[ii        ];
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#else///USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot,  body, file, filenum);
#endif//USE_HDF5_FORMAT
    if( unit_read != unit ){
      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
    }/* if( unit_read != unit ){ */
    /* sort by particle index */
    qsort(body, Ntot, sizeof(nbody_particle), idxAscendingOrder);
    //---------------------------------------------------------------------
    /* plot spherical averaged profile */
#if 0
    for(int ii = 0; ii < NUM_SPH_HARMONICS; ii++)
      for(int jj = 0; jj < NINTBIN; jj++)
	printf("l = %d, m = %+d; Plm = %e, Re(exp) = %e, Im(exp) = %e\n", sph_harm[ii].l, sph_harm[ii].m, sph_harm[ii].Plm[jj], GSL_REAL(sph_harm[ii].zz[jj]), GSL_IMAG(sph_harm[ii].zz[jj]));
    exitMPI();
    exit(0);
#endif
    analyzeRadialProfile(Ntot, body, kind, group, &num, &rem, (const real)(radius * radius), &rad, &rholm, weight, sph_harm,
			 rad_tbl, costhetamin, costhetabin, phimin, phibin, rho_tbl, rholm_tbl);
#if 0
    for(int ii = 0; ii < kind; ii++){
      fprintf(stderr, "#%d-th component\n", ii);
      for(int jj = 0; jj < NRADBIN; jj++){
	fprintf(stderr, "%e", rad[group[ii].prf_head[0] + jj]);
	for(int kk = 0; kk < NUM_SPH_HARMONICS; kk++)
	  fprintf(stderr, "\t%e", rholm[group[ii].prf_head[kk] + jj]);
	fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
    }
    exitMPI();
    exit(0);
#endif
    plotRadialProfile(num,  kind, group, sph_harm, (PLFLT)time, rad, rholm, rhorange, file, ifile, argc, argv);
    rem += num;
    num = 0;
    //---------------------------------------------------------------------
    ifile += (int)(interval * mpi.size);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  freeParticleDataAoS(body);
#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
  freeSnapshotArray(hdf5_pos, hdf5_vel, hdf5_acc, hdf5_m, hdf5_pot, hdf5_idx);
#endif//USE_HDF5_FORMAT
  free(rad);  free(rholm);
  //-----------------------------------------------------------------------
  free(group);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  exitMPI();
  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotRadialProfile
(const int ntot, const int ngroup, model *group, Ylm sph_harm[], PLFLT time,
 real *rad, real *rho, PLplotPltRange rhobox,
 char *file, const int filenum, int argc, char **argv)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* global setting(s) */
  const PLINT num_sph_harmonics = NUM_SPH_HARMONICS - 1;
  //-----------------------------------------------------------------------
  /* maximum number of data kinds */
  const PLINT  pkind = ngroup * num_sph_harmonics;
#if 1
  const PLINT  lkind = 0;
#else
  const PLINT  lkind = ngroup * num_sph_harmonics;
#endif
  const PLBOOL unify = true;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* preparation for dataset */
  //-----------------------------------------------------------------------
  /* memory allocation for index */
  PLINT *npt;  allocPLINT(&npt, pkind);
  PLINT *nln;  allocPLINT(&nln, lkind);
  //-----------------------------------------------------------------------
  /* set number of data points */
  for(PLINT ii = 0; ii < pkind; ii++)    npt[ii] = (PLINT)NRADBIN;
  for(PLINT ii = 0; ii < lkind; ii++)    nln[ii] = (PLINT)NRADBIN;
  //-----------------------------------------------------------------------
  /* memory allocation for data */
  PLFLT *_rad_p;  allocPLFLT(&_rad_p, ntot);  PLFLT **rad_p;  allocPointer4PLFLT(&rad_p, pkind);
  PLFLT *_rho_p;  allocPLFLT(&_rho_p, ntot);  PLFLT **rho_p;  allocPointer4PLFLT(&rho_p, pkind);
  PLFLT *_rad_l;  allocPLFLT(&_rad_l, ntot);  PLFLT **rad_l;  allocPointer4PLFLT(&rad_l, lkind);
  PLFLT *_rho_l;  allocPLFLT(&_rho_l, ntot);  PLFLT **rho_l;  allocPointer4PLFLT(&rho_l, lkind);
  /* set pointers */
  for(PLINT ii = 0; ii < pkind; ii++){
    rad_p[ii] = _rad_p + ii * NRADBIN;
    rho_p[ii] = _rho_p + ii * NRADBIN;
  }
  for(PLINT ii = 0; ii < lkind; ii++){
    rad_l[ii] = _rad_l + ii * NRADBIN;
    rho_l[ii] = _rho_l + ii * NRADBIN;
  }
  /* data preparation for point(s) */
  if( pkind != 0 )
    for(int ii = 0; ii < ngroup; ii++)
      for(int jj = 0; jj < num_sph_harmonics; jj++){
	//-----------------------------------------------------------------
	const int kk = num_sph_harmonics * ii + num_sph_harmonics - 1 - jj;
	//-----------------------------------------------------------------
	for(int ll = 0; ll < group[ii].prf_num[1 + jj]; ll++){
	  //---------------------------------------------------------------
	  const int idx = group[ii].prf_head[1 + jj] + ll;
	  //---------------------------------------------------------------
	  rad_p[kk][ll] = (PLFLT)log10((double)rad[idx] *  length2astro);
	  rho_p[kk][ll] = (PLFLT)rho[idx] / (PLFLT)rho[group[ii].prf_head[0] + ll];
	}
	//-----------------------------------------------------------------
      }
  /* data preparation for line(s) */
  if( lkind != 0 )
    for(int ii = 0; ii < ngroup; ii++)
      for(int jj = 0; jj < num_sph_harmonics; jj++){
	//-----------------------------------------------------------------
	const int kk = num_sph_harmonics * ii + num_sph_harmonics - 1 - jj;
	//-----------------------------------------------------------------
	for(int ll = 0; ll < group[ii].prf_num[1 + jj]; ll++){
	  //---------------------------------------------------------------
	  const int idx = group[ii].prf_head[1 + jj] + ll;
	  //---------------------------------------------------------------
	  rad_l[kk][ll] = (PLFLT)log10((double)rad[idx] *  length2astro);
	  rho_l[kk][ll] = (PLFLT)rho[idx] / (PLFLT)rho[group[ii].prf_head[0] + ll];
	}
	//-----------------------------------------------------------------
      }
  //-----------------------------------------------------------------------
  /* set symbol style */
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  if( pkind != 0 )
    for(int ii = 0; ii < ngroup; ii++)
      for(int jj = 0; jj < num_sph_harmonics; jj++){
	//-----------------------------------------------------------------
	const int kk = ii * num_sph_harmonics + jj;
	//-----------------------------------------------------------------
	pt[kk].color =       PLplotColor     [(num_sph_harmonics - 1 - jj) % 15];
	pt[kk].scale =       PLplotSymbolSize[(num_sph_harmonics - 1 - jj) % 14] * 0.5;
	sprintf(pt[kk].type, PLplotSymbolType[(num_sph_harmonics - 1 - jj) % 14]);
	/* pt[kk].type  = (char *)PLplotSymbolType[(num_sph_harmonics - 1 - jj) % 14]; */
	//-----------------------------------------------------------------
      }
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  if( lkind != 0 )
    for(int ii = 0; ii < ngroup; ii++)
      for(int jj = 0; jj < num_sph_harmonics; jj++){
	//-----------------------------------------------------------------
	const int kk = ii * num_sph_harmonics + jj;
	//-----------------------------------------------------------------
	ls[kk].color = PLplotColor[(num_sph_harmonics - 1 - jj) % 15];
	ls[kk].style =             (num_sph_harmonics - 1 - jj) %  5 ;
	ls[kk].width = THIN_LINE;
	//-----------------------------------------------------------------
      }
  //-----------------------------------------------------------------------
  /* set labels */
  char radlab[PLplotCharWords];  sprintf(radlab, "#fir");
  char rholab[PLplotCharWords];  sprintf(rholab, "|#fi#gr#dlm#u#fr / #fi#gr#fr#d00#u|");
  //-----------------------------------------------------------------------
  /* set caption(s) */
  PLplotCaption rhocap;  setDefaultCaption(&rhocap);  rhocap.write = false;
  sprintf(rhocap.side, "%s", "t");
  sprintf(rhocap.text, "%s = %e", "#fit#fr", (double)time * time2astro);
  //-----------------------------------------------------------------------
  /* set legends */
  PLplotLegend rholeg;  setDefaultLegend(&rholeg, false);  rholeg.write = true;
  rholeg.pos = PL_POSITION_RIGHT | PL_POSITION_BOTTOM | PL_POSITION_OUTSIDE;
  char *rholegTxt;
  {
    PLINT dkind = (unify) ? (IMAX(pkind, lkind)) : (pkind + lkind);
    allocChar4PLplot(&rholegTxt, dkind);    allocPointer4Char4PLplot(&(rholeg.text), dkind);    assignChar4PLplot(dkind, rholeg.text, rholegTxt);
  }
  for(int ii = 0; ii < ngroup; ii++)
    for(int jj = 0; jj < num_sph_harmonics; jj++)
      sprintf(rholeg.text[num_sph_harmonics * ii + num_sph_harmonics - 1 - jj], "(%d, %d)", sph_harm[1 + jj].l, sph_harm[1 + jj].m);
  //-----------------------------------------------------------------------
  char title[PLplotCharWords];
  sprintf(title, "#fit#fr = %6.2f", (PLFLT)((double)time * time2astro));
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation to exploit multiple-plot function */
  //-----------------------------------------------------------------------
  /* maximum number of panels in each direction */
  const PLINT nxpanel = 1;
  const PLINT nypanel = ngroup;
  //-----------------------------------------------------------------------
  /* specify plot or skip panel */
  PLBOOL *plot;
  allocPLBOOL(&plot, nxpanel * nypanel);
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
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel * nypanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel * nypanel);
  assignChar4PLplot(nxpanel * nypanel, figname, _figname);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* configuration about radial profile of particle distribution  */
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
      nlkind[idx] = (lkind != 0) ? num_sph_harmonics : 0;
      lnum  [idx] = &nln  [(lkind != 0) ? (idx * num_sph_harmonics) : 0];
      line  [idx] = &ls   [(lkind != 0) ? (idx * num_sph_harmonics) : 0];
      lx    [idx] = &rad_l[(lkind != 0) ? (idx * num_sph_harmonics) : 0];
      ly    [idx] = &rho_l[(lkind != 0) ? (idx * num_sph_harmonics) : 0];
      //-------------------------------------------------------------------
      /* point setting(s) */
      npkind[idx] = (pkind != 0) ? num_sph_harmonics : 0;
      pnum  [idx] = &npt  [(pkind != 0) ? (idx * num_sph_harmonics) : 0];
      point [idx] = &pt   [(pkind != 0) ? (idx * num_sph_harmonics) : 0];
      px    [idx] = &rad_p[(pkind != 0) ? (idx * num_sph_harmonics) : 0];
      py    [idx] = &rho_p[(pkind != 0) ? (idx * num_sph_harmonics) : 0];
      //-------------------------------------------------------------------
      /* errorbar setting(s) */
      errbar[idx] = 0;
      //-------------------------------------------------------------------
      /* plot area */
      range[idx] = rhobox;
      //-------------------------------------------------------------------
      /* legend(s) */
      uni[idx] = unify;
      //-------------------------------------------------------------------
      /* label setting(s) */
      xlabel[idx] = radlab;
      ylabel[idx] = rholab;
      //-------------------------------------------------------------------
      /* miscs */
      cap[idx] = rhocap;
      leg[idx] = rholeg;
      //-------------------------------------------------------------------
      /* file name */
      sprintf(figname[idx], "%s_%s_%s%d_%.3d", file, "multipole", "series", idx, filenum);
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
  //-----------------------------------------------------------------------
  free(_rad_p);  free(rad_p);  free(_rad_l);  free(rad_l);
  free(_rho_p);  free(rho_p);  free(_rho_l);  free(rho_l);
  free(pt);  free(ls);
  //-----------------------------------------------------------------------
  free(rholegTxt);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline gsl_complex _gaussQuad2d(real rho[][NINTBIN], real weight[], const Ylm sph_harm, const int ll)
{
  //-----------------------------------------------------------------------
  /* const double mns = 0.5 * (ymax - ymin);/\* -pi -> pi := pi *\/ */
  //-----------------------------------------------------------------------
  gsl_complex sum = gsl_complex_rect(0.0, 0.0);
  for(int mm = 0; mm < NINTBIN; mm++)
    sum = gsl_complex_add(sum, gsl_complex_mul_real(gsl_complex_conjugate(sph_harm.zz[mm]), weight[mm] * rho[ll][mm] * sph_harm.Plm[ll]));
  sum = gsl_complex_mul_real(sum, M_PI);
  //-----------------------------------------------------------------------
#if 0
  for(int mm = 0; mm < NINTBIN; mm++)
  if( fabs(weight[mm] * rho[ll][mm] * sph_harm.Plm[ll]) > 0.0 )
    printf("%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", ll, mm, weight[mm], rho[ll][mm], sph_harm.Plm[ll], GSL_REAL(gsl_complex_conjugate(sph_harm.zz[mm])), GSL_IMAG(gsl_complex_conjugate(sph_harm.zz[mm])), GSL_REAL(gsl_complex_mul_real(gsl_complex_conjugate(sph_harm.zz[mm]), weight[mm] * rho[ll][mm] * sph_harm.Plm[ll])), GSL_IMAG(gsl_complex_mul_real(gsl_complex_conjugate(sph_harm.zz[mm]), weight[mm] * rho[ll][mm] * sph_harm.Plm[ll])));
    /* printf("%d\t%e\n", mm, weight[mm] * rho[ll][mm] * sph_harm.Plm[ll]); */
  /* exit(1); */
#endif
#if 0
  if( gsl_complex_abs(sum) > 0.0 )
    printf("%e, %e\n", GSL_REAL(sum), GSL_IMAG(sum));
#endif
  return (sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
gsl_complex gaussQuad2d(real rho[][NINTBIN], real weight[], const Ylm sph_harm)
{
  //-----------------------------------------------------------------------
  /* const double mns = 0.5 * (xmax - xmin);/\* xmax = 1, xmin = -1 --> mns = 1.0 *\/ */
  //-----------------------------------------------------------------------
  gsl_complex sum = gsl_complex_rect(0.0, 0.0);
  for(int ll = 0; ll < NINTBIN; ll++)
    sum = gsl_complex_add(sum, gsl_complex_mul_real(_gaussQuad2d(rho, weight, sph_harm, ll), weight[ll]));
  //-----------------------------------------------------------------------
#if 0
  printf("%e, %e\n", GSL_REAL(sum), GSL_IMAG(sum));
#endif
  return (sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline int bisection(const real val, real min[])
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = NINTBIN - 1;
  //-----------------------------------------------------------------------
  if( val < min[ll] ){      return (ll);    }
  if( val > min[rr] ){      return (rr);    }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)(ll + rr)) >> 1;
    //---------------------------------------------------------------------
    if( (min[cc] - val) * (min[ll] - val) <= ZERO )      rr = (int)cc;
    else                                                 ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr )      return (ll);
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline real acot(real xx){
  //-----------------------------------------------------------------------
  /* based on cot(x) = -tan(x - pi/2) */
  //-----------------------------------------------------------------------
  const real sign = (xx < ZERO) ? (-UNITY) : (UNITY);
  return (sign * (real)M_PI_2 + ATAN(-xx));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void analyzeRadialProfile
(int Np, nbody_particle *body_tot, const int kind, model *group, int *num, int *rem, const real r2max, real **rad, real **rho, real weight[], Ylm sph_harm[],
 real rad_tbl[], real costhetamin[], real costhetabin[], real phimin[], real phibin[], real rho_tbl[][NINTBIN], real rholm_tbl[][NRADBIN])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* shift the center to the center-of-mass of the system */
  //-----------------------------------------------------------------------
  /* calculate the center-of-mass of the system */
  real com[3] = {ZERO, ZERO, ZERO};
  real Mtot = ZERO;
  for(int ii = 0; ii < Np; ii++){
    const real r2 = body_tot[ii].x * body_tot[ii].x + body_tot[ii].y * body_tot[ii].y + body_tot[ii].z * body_tot[ii].z;
    if( r2 < FOUR * r2max ){
      com[0] += body_tot[ii].m * body_tot[ii].x;
      com[1] += body_tot[ii].m * body_tot[ii].y;
      com[2] += body_tot[ii].m * body_tot[ii].z;
      Mtot   += body_tot[ii].m;
    }
  }
  Mtot = UNITY / Mtot;  com[0] *= Mtot;  com[1] *= Mtot;  com[2] *= Mtot;
  /* calculate the distance from the center-of-mass of the system */
  for(int ii = 0; ii < Np; ii++){
    //---------------------------------------------------------------------
    body_tot[ii].x -= com[0];
    body_tot[ii].y -= com[1];
    body_tot[ii].z -= com[2];
    //---------------------------------------------------------------------
    const real R2 = body_tot[ii].x * body_tot[ii].x + body_tot[ii].y * body_tot[ii].y;
    const real z2 = body_tot[ii].z * body_tot[ii].z;
    body_tot[ii].ax = R2 + z2;
    body_tot[ii].ay = R2;
    body_tot[ii].az = z2;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* make radial profile */
  //-----------------------------------------------------------------------
  *num = 0;
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < kind; kk++){
    //---------------------------------------------------------------------
    /* sort the array in ascending distance from the center */
    nbody_particle *body;
    body = &body_tot[group[kk].head];
    qsort(body, group[kk].num, sizeof(nbody_particle), radAscendingOrder);
    //---------------------------------------------------------------------
    const int nunit = (int)group[kk].num / NRADBIN;
    int head = 0;
    real inner = ZERO;
    for(int ii = 0; ii < NRADBIN; ii++){
      //-------------------------------------------------------------------
      /* create table of volume density */
      //-------------------------------------------------------------------
      const int tail = nunit * (ii + 1);
      if( (ulong)tail - 1 >= group[kk].num )	break;
      //-------------------------------------------------------------------
      real outer = SQRT(body[tail - 1].ax);
      rad_tbl[ii] = HALF * (inner + outer);
      //-------------------------------------------------------------------
      /* commit particles to the volume density table */
      for(int jj = 0; jj < NINTBIN; jj++)
	for(int kk = 0; kk < NINTBIN; kk++)
	  rho_tbl[jj][kk] = ZERO;
      for(int jj = head; jj < tail; jj++){
	//-----------------------------------------------------------------
	const real costheta = body[jj].z * RSQRT(body[jj].x * body[jj].x + body[jj].y * body[jj].y + body[jj].z * body[jj].z);
	const real      phi = (FABS(body[jj].y) < FABS(body[jj].x)) ? (ATAN(body[jj].y / body[jj].x)) : (acot(body[jj].x / body[jj].y));
	//-----------------------------------------------------------------
	const int ll = bisection(costheta, costhetamin);
	const int mm = bisection(     phi,      phimin);
	rho_tbl[ll][mm] += body[jj].m;
	//-----------------------------------------------------------------
      }/* for(int jj = head; jj < tail; jj++){ */
      //-------------------------------------------------------------------
      /* devide the mass by the volume of the corresponding grid */
      const real r2dr = QUARTER * (outer * outer - inner * inner) * (inner + outer);
      for(int jj = 0; jj < NINTBIN; jj++)
	for(int kk = 0; kk < NINTBIN; kk++)
	  rho_tbl[jj][kk] /= r2dr * costhetabin[jj] * phibin[kk];
      //-------------------------------------------------------------------
      inner = outer;
      head = tail;
      //-------------------------------------------------------------------
#if 0
      for(int jj = 0; jj < NINTBIN; jj++)
	for(int kk = 0; kk < NINTBIN; kk++)
	  printf("%d\t%d\t%e\n", jj, kk, rho_tbl[jj][kk]);
      exit(0);
#endif
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* calculate multipole moments */
      //-------------------------------------------------------------------
      for(int jj = 0; jj < NUM_SPH_HARMONICS; jj++)
	rholm_tbl[jj][ii] = (real)gsl_complex_abs(gaussQuad2d(rho_tbl, weight, sph_harm[jj]));
#if 0
      for(int jj = 0; jj < NUM_SPH_HARMONICS; jj++)
	printf("%d\t%e\n", jj, rholm_tbl[jj][ii]);
      exit(1);
#endif
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* commit multipole moments to the table */
    //---------------------------------------------------------------------
    for(int ii = 0; ii < NUM_SPH_HARMONICS; ii++){
      //-------------------------------------------------------------------
      real *prad, *prho;
      prad = *rad;
      prho = *rho;
      group[kk].prf_head[ii] = *num;
      //-------------------------------------------------------------------
      for(int jj = 0; jj < NRADBIN; jj++){
	//-----------------------------------------------------------------
	/* check # of unused elements */
	if( *rem == 0 ){
	  enlargeProfileArray(*num + NAllocUnit, rad, rho);
	  *rem += NAllocUnit;
	  prad = *rad;
	  prho = *rho;
	}
	//-----------------------------------------------------------------
	prad[*num] =   rad_tbl    [jj];
	prho[*num] = rholm_tbl[ii][jj];
	//-----------------------------------------------------------------
	*num += 1;
	*rem -= 1;
	//-----------------------------------------------------------------
      }/* for(int jj = 0; jj < NRADBIN; jj++){ */
      //-------------------------------------------------------------------
      group[kk].prf_num[ii] = *num - (group[kk].prf_head[ii]);
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NUM_SPH_HARMONICS; ii++){ */
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < kind; kk++){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
