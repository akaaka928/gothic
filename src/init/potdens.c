/*************************************************************************\
 *                                                                       *
                  last updated on 2016/12/10(Sat) 15:25:05
 *                                                                       *
 *    Making Initial Condition Code of N-body Simulation                 *
 *       Poisson solver to yield potential--density pair                 *
 *           adopting Full Multigrid method and BiCGSTAB method          *
 *           preconditioner for BiCGSTAB is ILU(0)                       *
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
#include <omp.h>
//-------------------------------------------------------------------------
#include <gsl/gsl_integration.h>
//-------------------------------------------------------------------------
#include "macro.h"
#include "constants.h"
//-------------------------------------------------------------------------
#include "profile.h"
#include "abel.h"
#include "blas.h"
#include "spline.h"
#include "magi.h"
#include "potdens.h"
//-------------------------------------------------------------------------
#   if  defined(USE_ELLIPTIC_INTEGRAL) && !defined(USE_GD_FORM_POTENTIAL)
#       include <gsl/gsl_sf_ellint.h>
#endif//defined(USE_ELLIPTIC_INTEGRAL) && !defined(USE_GD_FORM_POTENTIAL)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
extern const real newton;
//-------------------------------------------------------------------------
extern double gsl_gaussQD_pos[NTBL_GAUSS_QD], gsl_gaussQD_weight[NTBL_GAUSS_QD];
//-------------------------------------------------------------------------
#define NMAX_GAUSS_QD_LOW (11)
#define NTBL_GAUSS_QD_LOW ((NMAX_GAUSS_QD_LOW >> 1) + (NMAX_GAUSS_QD_LOW & 1))
#define NINTBIN_LOW NMAX_GAUSS_QD_LOW
static double gsl_gaussQD_low_pos[NTBL_GAUSS_QD_LOW], gsl_gaussQD_low_weight[NTBL_GAUSS_QD_LOW];
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#ifdef __ICC
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#     pragma warning (disable:869)
#endif//__ICC
//-------------------------------------------------------------------------
static inline double getSmoothCutoff (const double RR, const double Rt, const double invDelta){  return (0.5 * erfc(2.0 * (RR - Rt) * invDelta));}
double getColumnDensityExp   (double RR, double invRd, const disk_util disk);
double getColumnDensityExp   (double RR, double invRd, const disk_util disk){
  return (exp(-                    RR * invRd                   )   * getSmoothCutoff(RR, disk.Rcutoff, disk.invRsmooth));}
double getColumnDensitySersic(double RR, double invRd, const disk_util disk);
double getColumnDensitySersic(double RR, double invRd, const disk_util disk){
  return (exp(-disk.sersic_b * pow(RR * invRd, disk.sersic_ninv))   * getSmoothCutoff(RR, disk.Rcutoff, disk.invRsmooth));}
double getColumnDensitySpline(double RR, double invRd, const disk_util disk);
double getColumnDensitySpline(double RR, double invRd, const disk_util disk){
  return (getCubicSpline1D(RR, disk.num, disk.xx, disk.ff, disk.f2) * getSmoothCutoff(RR, disk.Rcutoff, disk.invRsmooth));}
//-------------------------------------------------------------------------
static double invzd2invz0;
static inline void setz0inv4SechProfile(void){  invzd2invz0 = 2.0 * acosh(sqrt(M_E));}
double getVerticalDensity(const double  zz, const double invzd, const disk_util disk);
double getVerticalDensity(const double  zz, const double invzd, const disk_util disk){
  const double tmp = 1.0 / cosh(0.5 * zz * invzd * invzd2invz0);/* := sech() */
  return (tmp * tmp);
}
//-------------------------------------------------------------------------
#ifdef __ICC
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#     pragma warning ( enable:869)
#endif//__ICC
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline double func4encSigma(const double RR, const disk_data disk){  return (RR * disk.getColumnDensity(RR, disk.invRd, disk.util));}
//-------------------------------------------------------------------------
static inline double gaussQD4encSigma(const double min, const double max, const disk_data disk)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  //-----------------------------------------------------------------------
  if( NINTBIN & 1 )
    sum = gsl_gaussQD_weight[(NINTBIN >> 1)] * func4encSigma(pls + mns * gsl_gaussQD_pos[(NINTBIN >> 1)], disk);
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (func4encSigma(pls + mns * gsl_gaussQD_pos[ii], disk) + func4encSigma(pls - mns * gsl_gaussQD_pos[ii], disk));
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void   freeDiskProfile
(const int ndisk, disk_data  *disk,
 double  *hor, double  *ver, double  *node_hor, double  *node_ver,
 double  *pot, double  *rho0, double  *rho1, double  *rhoTot, double  *dPhidR, double  *d2PhidR2,
 double  *Sigma, double  *vsigz, double  *enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 double  *zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 double  *radSph, double  *rhoSph, double  *encSph,
 double  *spline_xx, double  *spline_ff, double  *spline_f2, double  *spline_bp)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  free(disk);
  free(hor);  free(ver);  free(node_hor);  free(node_ver);
  free(pot);  free(rho0);  free(rho1);  free(dPhidR);  free(d2PhidR2);
  if( ndisk > 1 )
    free(rhoTot);
  free(Sigma);  free(vsigz);  free(enc);
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  free(zd);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  free(radSph);  free(rhoSph);  free(encSph);
  free(spline_xx);  free(spline_ff);  free(spline_f2);  free(spline_bp);
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void allocDiskProfile
(const int ndisk, disk_data **disk, profile_cfg *disk_cfg, int *maxLev, profile **disk_prf, const int skind, const double logrbin, const double invlogrbin,
 double **hor, double **ver, double **node_hor, double **node_ver,
 double **pot, double **rho0, double **rho1, double **rhoTot, double **dPhidR, double **d2PhidR2,
 double **Sigma, double **vsigz, double **enc,
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
 double **zd,
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
 double **radSph, double **rhoSph, double **encSph,
 double **spline_xx, double **spline_ff, double **spline_f2, double **spline_bp)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* determine maximum level of the nested grid */
  //-----------------------------------------------------------------------
  /* configuration of coarsest grid */
  double Rmax = 0.0;
  double zmax = 0.0;
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    if( disk_cfg[ii].cutoff ){      if( Rmax <  disk_cfg[ii].rc                    )	Rmax = disk_cfg[ii].rc                  ;    }
    else{                           if( Rmax < (disk_cfg[ii].rs * DISK_MAX_LENGTH) )	Rmax = disk_cfg[ii].rs * DISK_MAX_LENGTH;    }
    if(                                 zmax < (disk_cfg[ii].zd * DISK_MAX_LENGTH) )    zmax = disk_cfg[ii].zd * DISK_MAX_LENGTH;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  int log2hmax;
  double maxLz, maxLR;
  if( Rmax > ldexp(zmax, NHOR_OVER_NVER) ){    log2hmax = (int)ceil(log2(DISK_MAX_SAFETY * Rmax));    maxLR = ldexp(1.0, log2hmax);    maxLz = ldexp(maxLR, -NHOR_OVER_NVER);  }
  else{                                        log2hmax = (int)ceil(log2(DISK_MAX_SAFETY * zmax));    maxLz = ldexp(1.0, log2hmax);    maxLR = ldexp(maxLz,  NHOR_OVER_NVER);  }
  const double hh = maxLz / (double)NDISKBIN_VER;
  log2hmax = (int)nearbyint(log2(hh));
  //-----------------------------------------------------------------------
  /* configuration of finest grid */
  double Rmin = DBL_MAX;
  double zmin = DBL_MAX;
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    if( Rmin > disk_cfg[ii].rs )      Rmin = disk_cfg[ii].rs;
    if( zmin > disk_cfg[ii].zd )      zmin = disk_cfg[ii].zd;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  const int log2hmin = (int)floor(log2(DISK_MIN_LENGTH * fmin(Rmin, zmin)));
  /* const double hmin = ldexp(1.0, log2hmin); */
  //-----------------------------------------------------------------------
  /* hmin corresponds to 2^(1-Lmax) h */
  *maxLev = 1 - log2hmin + log2hmax;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* allocate arrays to store physical quantities */
  //-----------------------------------------------------------------------
  /* horizontal and vertical axes */
  *hor      = (double *)malloc((*maxLev) *  NDISKBIN_HOR      * sizeof(double));  if( *     hor == NULL ){    __KILL__(stderr, "ERROR: failure to allocate hor\n");  }
  *ver      = (double *)malloc((*maxLev) *  NDISKBIN_VER      * sizeof(double));  if( *     ver == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ver\n");  }
  *node_hor = (double *)malloc((*maxLev) * (NDISKBIN_HOR + 1) * sizeof(double));  if( *node_hor == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node_hor\n");  }
  *node_ver = (double *)malloc((*maxLev) * (NDISKBIN_VER + 1) * sizeof(double));  if( *node_ver == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node_ver\n");  }
  //-----------------------------------------------------------------------
  /* potential-density pair */
  *pot  = (double *)malloc(        (*maxLev) * NDISKBIN_HOR *  NDISKBIN_VER      * sizeof(double));  if( *pot  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pot\n" );  }
  *rho0 = (double *)malloc(ndisk * (*maxLev) * NDISKBIN_HOR * (NDISKBIN_VER + 1) * sizeof(double));  if( *rho0 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rho0\n");  }
  *rho1 = (double *)malloc(ndisk * (*maxLev) * NDISKBIN_HOR * (NDISKBIN_VER + 1) * sizeof(double));  if( *rho1 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rho1\n");  }
  /* multiple component disk */
  if( ndisk > 1 ){
    *rhoTot = (double *)malloc((*maxLev) * NDISKBIN_HOR * NDISKBIN_VER * sizeof(double));    if( *rhoTot == NULL ){      __KILL__(stderr, "ERROR: failure to allocate rhoTot\n");    }
  }/* if( ndisk > 1 ){ */
  //-----------------------------------------------------------------------
#ifndef USE_POTENTIAL_SCALING_SCHEME
  * dPhidR  = (double *)malloc((*maxLev) * NDISKBIN_HOR * NDISKBIN_VER * sizeof(double));
  *d2PhidR2 = (double *)malloc((*maxLev) * NDISKBIN_HOR * NDISKBIN_VER * sizeof(double));
#else///USE_POTENTIAL_SCALING_SCHEME
  * dPhidR  = (double *)malloc((*maxLev) * NDISKBIN_HOR                * sizeof(double));
  *d2PhidR2 = (double *)malloc((*maxLev) * NDISKBIN_HOR                * sizeof(double));
#endif//USE_POTENTIAL_SCALING_SCHEME
  if( * dPhidR  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate dPhidR\n");  }
  if( *d2PhidR2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate d2PhidR2\n");  }
  //-----------------------------------------------------------------------
  /* horizontal profile */
  *Sigma  = (double *)malloc(ndisk * (*maxLev) *  NDISKBIN_HOR      * sizeof(double));  if( *Sigma == NULL ){    __KILL__(stderr, "ERROR: failure to allocate Sigma\n");  }
  *vsigz  = (double *)malloc(ndisk * (*maxLev) *  NDISKBIN_HOR      * sizeof(double));  if( *vsigz == NULL ){    __KILL__(stderr, "ERROR: failure to allocate vsigz\n");  }
  *enc    = (double *)malloc(ndisk * (*maxLev) * (NDISKBIN_HOR + 1) * sizeof(double));  if( * enc  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate enc\n");  }
  //-----------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  *zd = (double *)malloc(ndisk * (*maxLev) * NDISKBIN_HOR * sizeof(double));
  if( *zd == NULL ){    __KILL__(stderr, "ERROR: failure to allocate zd\n");  }
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  //-----------------------------------------------------------------------
  /* for spherical averaged profile */
  *radSph = (double *)malloc(NDISKBIN_RAD * sizeof(double));  if( *radSph == NULL ){    __KILL__(stderr, "ERROR: failure to allocate radSph\n");  }
  *rhoSph = (double *)malloc(NDISKBIN_RAD * sizeof(double));  if( *rhoSph == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rhoSph\n");  }
  *encSph = (double *)malloc(NDISKBIN_RAD * sizeof(double));  if( *encSph == NULL ){    __KILL__(stderr, "ERROR: failure to allocate encSph\n");  }
  //-----------------------------------------------------------------------
  /* for spline fitting */
  /* note: NDISKBIN_HOR > NDISKBIN_VER */
  int Nspline = 2 + ((*maxLev) + 1) * (NDISKBIN_HOR >> 1);
  *spline_xx  = (double *)malloc(Nspline * CPUS_PER_PROCESS * sizeof(double));
  *spline_ff  = (double *)malloc(Nspline * CPUS_PER_PROCESS * sizeof(double));
  *spline_f2  = (double *)malloc(Nspline * CPUS_PER_PROCESS * sizeof(double));
  *spline_bp  = (double *)malloc(Nspline * CPUS_PER_PROCESS * sizeof(double));
  if( *spline_xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate spline_xx\n");  }
  if( *spline_ff == NULL ){    __KILL__(stderr, "ERROR: failure to allocate spline_ff\n");  }
  if( *spline_f2 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate spline_f2\n");  }
  if( *spline_bp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate spline_bp\n");  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set horizontal and vertical axes */
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < *maxLev; ii++){
    //---------------------------------------------------------------------
    const double bin = ldexp(hh, -ii);
    //---------------------------------------------------------------------
    (*node_hor)[INDEX2D(*maxLev, NDISKBIN_HOR + 1, ii, 0)] = 0.0;
#pragma omp parallel for
    for(int jj = 1; jj < 1 + NDISKBIN_HOR; jj++)      (*node_hor)[INDEX2D(*maxLev, 1 + NDISKBIN_HOR, ii, jj)] = bin *        (double)jj;
#pragma omp parallel for
    for(int jj = 0; jj <     NDISKBIN_HOR; jj++)      (*     hor)[INDEX2D(*maxLev,     NDISKBIN_HOR, ii, jj)] = bin * (0.5 + (double)jj);
    (*node_ver)[INDEX2D(*maxLev, NDISKBIN_VER + 1, ii, 0)] = 0.0;
#pragma omp parallel for
    for(int jj = 1; jj < 1 + NDISKBIN_VER; jj++)      (*node_ver)[INDEX2D(*maxLev, 1 + NDISKBIN_VER, ii, jj)] = bin *        (double)jj;
#pragma omp parallel for
    for(int jj = 0; jj <     NDISKBIN_VER; jj++)      (*     ver)[INDEX2D(*maxLev,     NDISKBIN_VER, ii, jj)] = bin * (0.5 + (double)jj);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < *maxLev; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* allocate utility structure and commit arrays */
  //-----------------------------------------------------------------------
  *disk = (disk_data *)malloc(ndisk * sizeof(disk_data));  if( *disk == NULL ){    __KILL__(stderr, "ERROR: failure to allocate disk\n");  }
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    /* commit disk properties */
    //---------------------------------------------------------------------
    (*disk)[ii].cfg   = &disk_cfg[ii];
    (*disk)[ii].Rmax  = maxLR;
    (*disk)[ii].zmax  = maxLz;
    (*disk)[ii].hh    = hh;
    (*disk)[ii].invRd = 1.0 / disk_cfg[ii].rs;
    /* (*disk)[ii].prf   = &((*disk_prf)[ii]); */
    (*disk)[ii].prf   = disk_prf[skind + ii];
    (*disk)[ii].logrbin = logrbin;
    (*disk)[ii].invlogrbin = invlogrbin;
    //---------------------------------------------------------------------
    /* settings for Sersic profile */
    if( disk_cfg[ii].kind == SERSIC ){
      (*disk)[ii].util.sersic_ninv = 1.0 / disk_cfg[ii].n_sersic;
      (*disk)[ii].util.sersic_b    =       disk_cfg[ii].b_sersic;
    }/* if( disk_cfg[ii].kind == SERSIC ){ */
    //---------------------------------------------------------------------
    /* additional setting about density cutoff in the horizontal direction */
    (*disk)[ii].util.Rcutoff    = Rmax * 10.0;/* i.e., the point at infinity */
    (*disk)[ii].util.invRsmooth = (*disk)[ii].invRd;
    if( disk_cfg[ii].cutoff ){
      //-------------------------------------------------------------------
      (*disk)[ii].util.Rcutoff    =       disk_cfg[ii].rc;
      (*disk)[ii].util.invRsmooth = 1.0 / disk_cfg[ii].rc_width;
      //-------------------------------------------------------------------
    }/* if( disk_cfg[ii].cutoff ){ */
    //---------------------------------------------------------------------
    /* initialize Toomre's Q-value analysis */
    (*disk)[ii].cfg->vcirc_max   = -1.0;
    (*disk)[ii].cfg->vcirc_max_R = -1.0;
    (*disk)[ii].cfg->Qmin = DBL_MAX;
    (*disk)[ii].cfg->passed = false;
    //---------------------------------------------------------------------
    /* set function pointer */
    switch( disk_cfg[ii].kind ){
    case EXP_DISK:      (*disk)[ii].getColumnDensity = getColumnDensityExp   ;      break;
    case   SERSIC:      (*disk)[ii].getColumnDensity = getColumnDensitySersic;      break;
    case TBL_DISK:      (*disk)[ii].getColumnDensity = getColumnDensitySpline;
      readColumnDensityTable4Disk((*disk)[ii].prf, disk_cfg[ii].rs, disk_cfg[ii].file, &((*disk)[ii].util.num), &((*disk)[ii].util.xx), &((*disk)[ii].util.ff), &((*disk)[ii].util.f2), &((*disk)[ii].util.bp));
      break;
    default:      __KILL__(stderr, "ERROR: undefined model(%d) is specified as disk profile\n", disk_cfg[ii].kind);      break;
    }/* switch( disk_cfg[ii].kind ){ */
    //---------------------------------------------------------------------
    /* set Sigma0 */
#ifdef  NDIVIDE_GAUSSQD4DISK
    const double Rbin = Rmax / (double)NDIVIDE_GAUSSQD4DISK;
    double sum = 0.0;
    double Rin = 0.0;
    for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){
      sum += gaussQD4encSigma(Rin, Rin + Rbin, (*disk)[ii]);
      Rin += Rbin;
    }/* for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){ */
    disk_cfg[ii].Sigma0 = disk_cfg[ii].Mtot / (2.0 * M_PI * sum);
#else///NDIVIDE_GAUSSQD4DISK
    disk_cfg[ii].Sigma0 = disk_cfg[ii].Mtot / (2.0 * M_PI * gaussQD4encSigma(0.0, Rmax, (*disk)[ii]));
#endif//NDIVIDE_GAUSSQD4DISK
    //---------------------------------------------------------------------
#if 0
    fprintf(stdout, "Sigma0 = %e, Mtot = %e\n", disk_cfg[ii].Sigma0, 2.0 * M_PI * disk_cfg[ii].Sigma0 * sum);
    fflush(stdout);
    exit(0);
#endif
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* common arrays for all components */
    //---------------------------------------------------------------------
    (*disk)[ii].hor = *hor;
    (*disk)[ii].ver = *ver;
    (*disk)[ii].node_hor = *node_hor;
    (*disk)[ii].node_ver = *node_ver;
    (*disk)[ii].pot = *pot;
    //---------------------------------------------------------------------
    (*disk)[ii]. dPhidR  = * dPhidR ;
    (*disk)[ii].d2PhidR2 = *d2PhidR2;
    //---------------------------------------------------------------------
    /* spherical average profile */
    (*disk)[ii].radSph = *radSph;
    (*disk)[ii].rhoSph = *rhoSph;
    (*disk)[ii].encSph = *encSph;
    //---------------------------------------------------------------------
    /* spline fitting */
    (*disk)[ii].spline_xx = *spline_xx;
    (*disk)[ii].spline_ff = *spline_ff;
    (*disk)[ii].spline_f2 = *spline_f2;
    (*disk)[ii].spline_bp = *spline_bp;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* individual arrays for each component */
    //---------------------------------------------------------------------
    (*disk)[ii].rho0 = &((*rho0)[INDEX(ndisk, *maxLev, NDISKBIN_HOR * (NDISKBIN_VER + 1), ii, 0, 0)]);    (*disk)[ii].rho    = &(*disk)[ii].rho0;
    (*disk)[ii].rho1 = &((*rho1)[INDEX(ndisk, *maxLev, NDISKBIN_HOR * (NDISKBIN_VER + 1), ii, 0, 0)]);    (*disk)[ii].rhoSum = &(*disk)[ii].rho1;
    //---------------------------------------------------------------------
    (*disk)[ii].Sigma  = &((*Sigma)[INDEX(ndisk, *maxLev, NDISKBIN_HOR    , ii, 0, 0)]);
    (*disk)[ii].sigmaz = &((*vsigz)[INDEX(ndisk, *maxLev, NDISKBIN_HOR    , ii, 0, 0)]);
    (*disk)[ii].enc    = &((* enc )[INDEX(ndisk, *maxLev, NDISKBIN_HOR + 1, ii, 0, 0)]);
    //---------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
    (*disk)[ii].zd = &((*zd)[INDEX(ndisk, *maxLev, NDISKBIN_HOR, ii, 0, 0)]);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* additional setting in case of multiple component disk */
    if( ndisk > 1 )
      (*disk)[ii].rhoTot = *rhoTot;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* initialization as a first touch */
    //---------------------------------------------------------------------
    for(int ll = 0; ll < *maxLev; ll++)
#pragma omp parallel for
      for(int jj = 0; jj < NDISKBIN_HOR; jj++){
	(*disk)[ii].Sigma [INDEX2D(maxLev, NDISKBIN_HOR    , ll, jj)] = 0.0;
	(*disk)[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR    , ll, jj)] = 0.0;
	(*disk)[ii].enc   [INDEX2D(maxLev, NDISKBIN_HOR + 1, ll, jj)] = 0.0;
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
	(*disk)[ii].zd    [INDEX2D(maxLev, NDISKBIN_HOR, ll, jj)] = 0.0;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
	for(int kk = 0; kk < NDISKBIN_VER; kk++){
	  (*disk)[ii].rho0[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
	  (*disk)[ii].rho1[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
	}/* for(int kk = 0; kk < NDISKBIN_VER; kk++){ */
      }/* for(int jj = 0; jj < NDISKBIN_HOR; jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* initialization as a first touch */
  //-----------------------------------------------------------------------
  for(int ll = 0; ll < *maxLev; ll++)
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_HOR; jj++)
      for(int kk = 0; kk < NDISKBIN_VER; kk++){
	(*pot)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
#ifndef USE_POTENTIAL_SCALING_SCHEME
	(* dPhidR )[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
	(*d2PhidR2)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
#endif//USE_POTENTIAL_SCALING_SCHEME
      }/* for(int kk = 0; kk < NDISKBIN_VER; kk++){ */
  //-----------------------------------------------------------------------
#pragma omp parallel num_threads(CPUS_PER_PROCESS)
  {
    const int tidx = omp_get_thread_num();
    const int num = 2 + ((*maxLev) + 1) * (NDISKBIN_HOR >> 1);
    for(int ii = 0; ii < num; ii++){
      (*spline_ff)[ii + tidx * num] = 0.0;
      (*spline_f2)[ii + tidx * num] = 0.0;
      (*spline_bp)[ii + tidx * num] = 0.0;
    }/* for(int ii = 0; ii < num; ii++){ */
  }
  //-----------------------------------------------------------------------
  if( ndisk > 1 )
    for(int ll = 0; ll < *maxLev; ll++)
#pragma omp parallel for
      for(int jj = 0; jj < NDISKBIN_HOR; jj++)
	for(int kk = 0; kk < NDISKBIN_VER; kk++)
	  (*rhoTot)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, jj, kk)] = 0.0;
  //-----------------------------------------------------------------------
#if 0
  fprintf(stdout, "# NR = %d, Nz = %d\n", NDISKBIN_HOR, NDISKBIN_VER);
  fprintf(stdout, "# Rmax = %e, zmax = %e, log2hmax = %d\n", Rmax, zmax, log2hmax);
  fprintf(stdout, "# maxLR = %e, maxLz = %e, hh = %e\n", maxLR, maxLz, hh);
  fprintf(stdout, "# Rmin = %e, zmin = %e, log2hmin = %d\n", Rmin, zmin, log2hmin);
  fprintf(stdout, "# maxLev = %d, hmin = %e\n", *maxLev, ldexp(1.0, log2hmin));
  fflush(NULL);
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline int bisection(const double val, const int num, double * restrict tab, const bool logtbl, const double invbin, double * restrict ratio)
{
  //-----------------------------------------------------------------------
  int ll =       0;
  int rr = num - 1;
  //-----------------------------------------------------------------------
  /* prohibit extraporation */
  if( val < tab[ll] + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( val > tab[rr] - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)(ll + rr)) >> 1;
    //---------------------------------------------------------------------
    if( (tab[cc] - val) * (tab[ll] - val) <= 0.0)      rr = (int)cc;
    else                                               ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr ){
      *ratio = (logtbl ? (log10(val / tab[ll])) : (val - tab[ll])) * invbin;
      return (ll);
    }/* if( (1 + ll) == rr ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline int findIdxSpherical(const double rad, profile *prf, const double invlogbin, double *ratio)
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = 3 + NRADBIN;
  //-----------------------------------------------------------------------
  if( rad < prf[ll].rad + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( rad > prf[rr].rad - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)ll + (uint)rr) >> 1;
    //---------------------------------------------------------------------
    if( (prf[cc].rad - rad) * (prf[ll].rad - rad) <= 0.0 )      rr = (int)cc;
    else                                                        ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr ){
      *ratio = log10(rad / prf[ll].rad) * invlogbin;
      return (ll);
    }/* if( (1 + ll) == rr ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double Phi_spherical(const double rad, profile *sph, const double invlogrbin_sph)
{
  //-----------------------------------------------------------------------
  double ratio;
  const int idx = findIdxSpherical(rad, sph, invlogrbin_sph, &ratio);
  //-----------------------------------------------------------------------
  return ((1.0 - ratio) * sph[idx].psi_tot + ratio * sph[1 + idx].psi_tot);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
//-------------------------------------------------------------------------
static inline int findIdxSphericalPsi(const double psi, profile *prf, double *ratio)
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = 3 + NRADBIN;
  //-----------------------------------------------------------------------
  if( psi > prf[ll].psi_tot - DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( psi < prf[rr].psi_tot + DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)ll + (uint)rr) >> 1;
    //---------------------------------------------------------------------
    if( (prf[cc].psi_tot - psi) * (prf[ll].psi_tot - psi) <= 0.0 )      rr = (int)cc;
    else                                                                ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr ){
      *ratio = (psi - prf[ll].psi_tot) / (prf[rr].psi_tot - prf[ll].psi_tot);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static const double diskDimmingHeight    = 16.0;
static const double diskDimmingHeightInv = 0.0625;
//-------------------------------------------------------------------------
static inline void getVariableDiskScaleHeight(const int ndisk, const int maxLev, disk_data * restrict disk, profile * restrict sph, const double invlogrbin_sph)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  for(int kk = 0; kk < ndisk; kk++){
    //---------------------------------------------------------------------
    for(int lev = 0; lev < maxLev; lev++){
      //-------------------------------------------------------------------
      double *RR;      RR = &(disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)]);
      double *zd;      zd = &(disk[kk]. zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)]);
      //-------------------------------------------------------------------
      const double zd0 = disk[kk].cfg->zd;
      const double zd2 = diskDimmingHeight * zd0 * diskDimmingHeight * zd0;
      //-------------------------------------------------------------------
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	//-----------------------------------------------------------------
	const double R2 = RR[ii] * RR[ii];
	//-----------------------------------------------------------------
	/* get potential @ the reference points (mid plane and dimming scale) */
	const double PsiMid = Phi_spherical(RR[ii], sph, invlogrbin_sph);
	const double PsiDim = Phi_spherical(sqrt(R2 + zd2), sph, invlogrbin_sph);
	//-----------------------------------------------------------------
	const double Psi = PsiMid + diskDimmingHeightInv * (PsiDim - PsiMid);
	//-----------------------------------------------------------------
	double ratio;
	const int irr = findIdxSphericalPsi(Psi, sph, &ratio);
	const double rr = (1.0 - ratio) * sph[irr].rad + ratio * sph[1 + irr].rad;
	//-----------------------------------------------------------------
	const double zd1 = sqrt(rr * rr - R2);
	zd[ii] = (zd0 < zd1) ? zd0 : zd1;
	//-----------------------------------------------------------------
#if 0
	fprintf(stderr, "%e\t%e\n", RR[ii], zd[ii]);
#endif
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
      //-------------------------------------------------------------------
    }/* for(int lev = 0; lev < maxLev; lev++){ */
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < ndisk; kk++){ */
  //-----------------------------------------------------------------------
#if 0
  fflush(stderr);
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline double gaussQuad1d4Rho(double (*func)(double, double, disk_util), const double min, const double max, const double xinv, const disk_util disk)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  if( NINTBIN_LOW & 1 )
    sum = gsl_gaussQD_low_weight[(NINTBIN_LOW >> 1)] * func(pls + mns * gsl_gaussQD_low_pos[(NINTBIN_LOW >> 1)], xinv, disk);
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_low_weight[ii] *
      (func(pls + mns * gsl_gaussQD_low_pos[ii], xinv, disk) +
       func(pls - mns * gsl_gaussQD_low_pos[ii], xinv, disk));
  //-----------------------------------------------------------------------
  /* finalization */
  return (sum * mns);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* assume potential generated by a Miyamoto & Nagai disk for initial guess to the solution of the Poisson equation */
//-------------------------------------------------------------------------
static inline void initPotentialField(double * restrict RR, double * restrict zz, double * restrict Phi, const int ndisk, disk_data * restrict disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* Miyamoto & Nagai disk */
#pragma omp parallel for
  for(int ii = 0; ii < NDISKBIN_HOR; ii++){
    //---------------------------------------------------------------------
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
      Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] = 0.0;
    //---------------------------------------------------------------------
    const double R2 = RR[ii] * RR[ii];
    //---------------------------------------------------------------------
    for(int jj = 0; jj < NDISKBIN_VER; jj++){
      //-------------------------------------------------------------------
      const double z2 = zz[jj] * zz[jj];
      //-------------------------------------------------------------------
      double sum = 0.0;
      for(int kk = 0; kk < ndisk; kk++){
	//-----------------------------------------------------------------
	const double Rs = disk[kk].cfg->rs + sqrt(z2 + disk[kk].cfg->zd * disk[kk].cfg->zd);
	sum += disk[kk].cfg->Mtot / sqrt(R2 + Rs * Rs);
	//-----------------------------------------------------------------
      }/* for(int kk = 0; kk < ndisk; kk++){ */
      Phi[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] -= (double)newton * sum;
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void setColumnDensityProfile(const int ndisk, const int maxLev, disk_data * restrict disk, profile * restrict sph, const double invlogrbin_sph)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  getVariableDiskScaleHeight(ndisk, maxLev, disk, sph, invlogrbin_sph);
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  //-----------------------------------------------------------------------
  /* set column density profile on the midplane */
  for(int kk = 0; kk < ndisk; kk++)
    for(int ll = 0; ll < maxLev; ll++){
      //-------------------------------------------------------------------
      const double hh = ldexp(disk[kk].hh, -ll);
      //-------------------------------------------------------------------
#if 0
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	const double R0 = disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)];
	disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)] = disk[kk].cfg->Sigma0 * gaussQD4encSigma(R0 - 0.5 * hh, R0 + 0.5 * hh, disk[kk]) / (R0 * hh);
      }
#else
      const double hinv = 1.0 / hh;
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++)
	disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)] = disk[kk].cfg->Sigma0 * gaussQuad1d4Rho(disk[kk].getColumnDensity, disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)] - 0.5 * hh, disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)] + 0.5 * hh, disk[kk].invRd, disk[kk].util) * hinv;
#endif
      //-------------------------------------------------------------------
    }/* for(int ll = 0; ll < maxLev; ll++){ */
  //-----------------------------------------------------------------------
#if 0
  for(int kk = 0; kk < ndisk; kk++)
    disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, 0, NDISKBIN_HOR - 1)] = 0.0;
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set a tentative density distribution as an initial value in the coarsest grid */
  /* initial vertical profile is assumed to be hyperbolic secant */
  //-----------------------------------------------------------------------
  setz0inv4SechProfile();
  for(int kk = 0; kk < ndisk; kk++){
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      //-------------------------------------------------------------------
      /* set scale height @ given R */
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
      const double invzd = 1.0 / disk[kk].zd[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)];
#else///ENABLE_VARIABLE_SCALE_HEIGHT
      const double invzd = 1.0 / disk[kk].cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
      const double rho0 = 0.5 * disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)] * invzd;
      //-------------------------------------------------------------------
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	(*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, jj)] = rho0 * gaussQuad1d4Rho(getVerticalDensity, disk[kk].ver[INDEX2D(maxLev, NDISKBIN_VER, 0, jj)] - 0.5 * disk[kk].hh, disk[kk].ver[INDEX2D(maxLev, NDISKBIN_VER, 0, jj)] + 0.5 * disk[kk].hh, invzd, disk[kk].util);
      //-------------------------------------------------------------------
#if 0
      (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, NDISKBIN_VER - 1)] = 0.0;
#endif
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* calculate column density using Simpson's rule */
#if 1
      double Sigma = 0.0;
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	Sigma += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, jj)];
      Sigma *= 2.0 * disk[kk].hh;/* 2.0 reflects plane symmetry about the equatorial plane (z = 0) */
#else
      double Sigma = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, 0)];
      for(int jj = 1; jj < NDISKBIN_VER - 1; jj++)
	Sigma += (double)(1 << (1 + (jj & 1))) * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, jj)];
      Sigma += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, NDISKBIN_VER - 1)];
      Sigma *= 2.0 * disk[kk].hh / 3.0;/* 2.0 reflects plane symmetry about the equatorial plane (z = 0) */
#endif
      /* calibrate column density */
      const double Mscale = disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)] / (DBL_MIN + Sigma);
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	(*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, jj)] *= Mscale;
      //-------------------------------------------------------------------
      initPotentialField(disk[0].hor, disk[0].ver, disk[0].pot, ndisk, disk);
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  }/* for(int kk = 0; kk < ndisk; kk++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void coarsePotential4boundary(const int lev_lrs, double * restrict pot)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const int lev_hrs = lev_lrs + 1;
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
    //---------------------------------------------------------------------
    const int il = ii << 1;
    const int ih = il + 1;
    //---------------------------------------------------------------------
    const int jj = (NDISKBIN_VER >> 1) - 1;
    //---------------------------------------------------------------------
    const int jl = jj << 1;
    const int jh = jl + 1;
    //---------------------------------------------------------------------
    pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] =
      (pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jh)] +
       pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jh)]) * 0.25;
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */
  //-----------------------------------------------------------------------
  {
    //---------------------------------------------------------------------
    const int ii = (NDISKBIN_HOR >> 1) - 1;
    //---------------------------------------------------------------------
    const int il = ii << 1;
    const int ih = il + 1;
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){
      //-------------------------------------------------------------------
      const int jl = jj << 1;
      const int jh = jl + 1;
      //-------------------------------------------------------------------
      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] =
	(pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jh)] +
	 pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jh)]) * 0.25;
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void coarsePotential(const int lev_lrs, double * restrict pot)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const int lev_hrs = lev_lrs + 1;
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
    //---------------------------------------------------------------------
    const int il = ii << 1;
    const int ih = il + 1;
    //---------------------------------------------------------------------
    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){
      //-------------------------------------------------------------------
      const int jl = jj << 1;
      const int jh = jl + 1;
      //-------------------------------------------------------------------
      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] =
	(pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jh)] +
	 pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jl)] + pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jh)]) * 0.25;
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void   finePotential(const int lev_hrs, double *pot)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const int lev_lrs = lev_hrs - 1;
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
    //---------------------------------------------------------------------
    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){
      //-------------------------------------------------------------------
      double dfdR, dfdz;
      dfdR = 0.125 * (pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii + 1, jj)] - pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, (ii > 0) ? (ii - 1) : (0), jj)]);
      dfdz = 0.125 * (pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj + 1)] - pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, (jj > 0) ? (jj - 1) : (0))]);
      //-------------------------------------------------------------------
      const double pot00 = pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)];
      //-------------------------------------------------------------------
      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs,      ii << 1 ,      jj << 1 )] = pot00 - dfdR - dfdz;
      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs,      ii << 1 , 1 + (jj << 1))] = pot00 - dfdR + dfdz;
      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, 1 + (ii << 1),      jj << 1 )] = pot00 + dfdR - dfdz;
      pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, 1 + (ii << 1), 1 + (jj << 1))] = pot00 + dfdR + dfdz;
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void coarseDensity(const int lev_lrs, const disk_data disk, double * restrict rho)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const int lev_hrs = lev_lrs + 1;
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
    //---------------------------------------------------------------------
    const int il = ii << 1;
    const int ih = il + 1;
    //---------------------------------------------------------------------
    const double R0 = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_lrs, ii)];
    const double Rl = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_hrs, il)];
    const double Rh = disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_hrs, ih)];
    const double Rinv = 0.25 / R0;
    //---------------------------------------------------------------------
    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){
      //-------------------------------------------------------------------
      const int jl = jj << 1;
      const int jh = jl + 1;
      //-------------------------------------------------------------------
      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] =
	((rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jl)] + rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, il, jh)]) * Rl +
	 (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jl)] + rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ih, jh)]) * Rh) * Rinv;
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */
  //-----------------------------------------------------------------------
#if 0
  const double hh_lrs = ldexp(disk.hh, -lev_lrs);
  const double hh_hrs = ldexp(disk.hh, -lev_hrs);
  double Mtot_lrs = 0.0;
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++)
    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++)
      Mtot_lrs += rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] * disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_lrs, ii)] * hh_lrs * hh_lrs;
  double Mtot_hrs = 0.0;
  for(int ii = 0; ii < NDISKBIN_HOR; ii++)
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
      Mtot_hrs += rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ii, jj)] * disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_hrs, ii)] * hh_hrs * hh_hrs;
  fprintf(stderr, "Mtot_hrs = %e\n", Mtot_hrs);
  fprintf(stderr, "Mtot_lrs = %e\n", Mtot_lrs);
  fflush(NULL);
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void   fineDensity(const int lev_hrs, const disk_data disk, double *rho)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const int lev_lrs = lev_hrs - 1;
  //-----------------------------------------------------------------------
#if 1
  //-----------------------------------------------------------------------
  const double hh_lrs = ldexp(disk.hh, -lev_lrs);
#pragma omp parallel for
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){
    //---------------------------------------------------------------------
    const double tmp = -hh_lrs / (3.0 * disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_lrs, ii)]);
    const double inner = (double)( 2 + 12 * ii) / (double)(3 + 12 * ii);
    const double outer = (double)(10 + 12 * ii) / (double)(9 + 12 * ii);
    //---------------------------------------------------------------------
    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){
      //-------------------------------------------------------------------
      double dfdR, dfdz;
      dfdR = 0.125 * (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii + 1, jj)] - rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, (ii > 0) ? (ii - 1) : (0), jj)]);
      dfdz = 0.125 * (rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj + 1)] - rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, (jj > 0) ? (jj - 1) : (0))]);
      const double rhoAvg = rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] + tmp * dfdR;
      //-------------------------------------------------------------------
      double rho00 = rhoAvg - inner * dfdR - dfdz;
      double rho01 = rhoAvg - inner * dfdR + dfdz;
      double rho10 = rhoAvg + outer * dfdR - dfdz;
      double rho11 = rhoAvg + outer * dfdR + dfdz;
      if( (rho00 < 0.0) || (rho01 < 0.0) || (rho10 < 0.0) || (rho11 < 0.0) )
      	rho00 = rho01 = rho10 = rho11 = rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)];
      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs,      ii << 1 ,      jj << 1 )] = rho00;
      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs,      ii << 1 , 1 + (jj << 1))] = rho01;
      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, 1 + (ii << 1),      jj << 1 )] = rho10;
      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, 1 + (ii << 1), 1 + (jj << 1))] = rho11;
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++){ */
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++){ */
  //-----------------------------------------------------------------------
#else
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < NDISKBIN_HOR; ii++)
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
      rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ii, jj)] = rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii >> 1, jj >> 1)];
  //-----------------------------------------------------------------------
#endif
  //-----------------------------------------------------------------------
#if 0
  const double hh_hrs = ldexp(disk.hh, -lev_hrs);
  double Mtot_lrs = 0.0;
  for(int ii = 0; ii < (NDISKBIN_HOR >> 1); ii++)
    for(int jj = 0; jj < (NDISKBIN_VER >> 1); jj++)
      Mtot_lrs += data[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_lrs, ii, jj)] * disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_lrs, ii)] * hh_lrs * hh_lrs;
  double Mtot_hrs = 0.0;
  for(int ii = 0; ii < NDISKBIN_HOR; ii++)
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
      Mtot_hrs += data[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev_hrs, ii, jj)] * disk.hor[INDEX2D(maxLev, NDISKBIN_HOR, lev_hrs, ii)] * hh_hrs * hh_hrs;
  fprintf(stderr, "Mtot_hrs = %e\n", Mtot_hrs);
  fprintf(stderr, "Mtot_lrs = %e\n", Mtot_lrs);

  fflush(NULL);
  exit(0);
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_GD_FORM_POTENTIAL
//-------------------------------------------------------------------------
static inline int findIdxSpherical4GDpot(const double rad, profile *prf, double *ratio)
{
  //-----------------------------------------------------------------------
  int ll = 0;
  int rr = NRADBIN - 1;
  //-----------------------------------------------------------------------
  if( rad < prf[ll].rho + DBL_EPSILON ){    *ratio = 0.0;    return (ll    );  }
  if( rad > prf[rr].rho - DBL_EPSILON ){    *ratio = 1.0;    return (rr - 1);  }
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    const uint cc = ((uint)ll + (uint)rr) >> 1;
    //---------------------------------------------------------------------
    if( (prf[cc].rho - rad) * (prf[ll].rho - rad) <= 0.0 )      rr = (int)cc;
    else                                                        ll = (int)cc;
    //---------------------------------------------------------------------
    if( (1 + ll) == rr ){
      *ratio = (rad - prf[ll].rho) / (prf[rr].rho - prf[ll].rho);
      return (ll);
    }/* if( (1 + ll) == rr ){ */
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double func4GDpot(const double RR, const double a2, const disk_data disk)
{
  //-----------------------------------------------------------------------
  return (RR * disk.getColumnDensity(RR, disk.invRd, disk.util) / sqrt(RR * RR - a2));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
double gaussQD4GDpot(const double min, const double max, const double a2, const disk_data disk);
double gaussQD4GDpot(const double min, const double max, const double a2, const disk_data disk)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  //-----------------------------------------------------------------------
  if( NINTBIN & 1 )
    sum = gsl_gaussQD_weight[(NINTBIN >> 1)] * func4GDpot(pls + mns * gsl_gaussQD_pos[(NINTBIN >> 1)], a2, disk);
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (func4GDpot(pls + mns * gsl_gaussQD_pos[ii], a2, disk) + func4GDpot(pls - mns * gsl_gaussQD_pos[ii], a2, disk));
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void calcGDpot
(const double zz, const int Nz, double * restrict ver, const double invdz, const double zp,
 const int ia, const double fa, const int iR, const double fR, const double aa, const double apR2, const double amR2,
 const int ndisk, disk_data * restrict disk, double * restrict invSigma, double * restrict ret)
{
  //-----------------------------------------------------------------------
  const double zpls2 = (zz - zp) * (zz - zp);  const double asinp = asin(2.0 * aa / (sqrt(zpls2 + apR2) + sqrt(zpls2 + amR2)));
  const double zmns2 = (zz + zp) * (zz + zp);  const double asinm = asin(2.0 * aa / (sqrt(zmns2 + apR2) + sqrt(zmns2 + amR2)));
  //-----------------------------------------------------------------------
  double fz;
  const int jz = bisection(zp, Nz, ver, false, invdz, &fz);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++){
    //---------------------------------------------------------------------
    const double zeta = invSigma[ii] *
      (((1.0 - fz) * (*disk[ii].rho)[INDEX2D(NR, Nz,     iR, jz)] + fz * (*disk[ii].rho)[INDEX2D(NR, Nz,     iR, 1 + jz)]) * (1.0 - fR) +
       ((1.0 - fz) * (*disk[ii].rho)[INDEX2D(NR, Nz, 1 + iR, jz)] + fz * (*disk[ii].rho)[INDEX2D(NR, Nz, 1 + iR, 1 + jz)]) *        fR   );
    //---------------------------------------------------------------------
    const double diff = (1.0 - fa) * disk[ii].prf[ia].psi + fa * disk[ii].prf[1 + ia].psi;
    //---------------------------------------------------------------------
    ret[ii] = zeta * diff * (asinp + asinm);
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void _gaussQuad2d4GDpot
(const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 const int ia, const double fa, const int iR, const double fR, const double aa, const double apR2, const double amR2,
 const int ndisk, disk_data * restrict disk, double * restrict sum, double * invSigma, double * restrict tmp)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (zmax - zmin);
  const double pls = 0.5 * (zmax + zmin);
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] = 0.0;
  if( NINTBIN & 1 ){
    //---------------------------------------------------------------------
    const double ww =             gsl_gaussQD_weight[(NINTBIN >> 1)];
    const double zp = pls + mns * gsl_gaussQD_pos   [(NINTBIN >> 1)];
    calcGDpot(zz, Nz, ver, invdz, zp, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, invSigma, tmp);
    for(int ii = 0; ii < ndisk; ii++)      sum[ii] = ww * tmp[ii];
    //---------------------------------------------------------------------
  }/* if( NINTBIN & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double ww = gsl_gaussQD_weight[ii];
    //---------------------------------------------------------------------
    double zp = pls + mns * gsl_gaussQD_pos[ii];
    calcGDpot(zz, Nz, ver, invdz, zp, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, invSigma, tmp);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
    zp = pls - mns * gsl_gaussQD_pos[ii];
    calcGDpot(zz, Nz, ver, invdz, zp, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, invSigma, tmp);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
    //---------------------------------------------------------------------
  }/* for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  /* finalization */
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] *= mns;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void _setRdep4GDpot
(const double aa, int * restrict ia, double * restrict fa, const int ndisk, disk_data * restrict disk, double * restrict invSigma,
 const double RR, int * restrict iR, double * restrict fR, const int NR, double * restrict hor, const double invdR, double * restrict apR2, double * restrict amR2)
{
  //-----------------------------------------------------------------------
  *apR2 = (aa + RR) * (aa + RR);
  *amR2 = (aa - RR) * (aa - RR);
  //-----------------------------------------------------------------------
  *iR = bisection(aa, NR, hor, false, invdR, fR);
  *ia = findIdxSpherical4GDpot(aa, disk[0].prf, fa);
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < ndisk; ii++)
    invSigma[ii] = 1.0 / (DBL_MIN + (1.0 - (*fR)) * disk[ii].Sigma[*iR] + (*fR) * disk[ii].Sigma[1 + (*iR)]);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void gaussQuad2d4GDpot
(const double RR, const int NR, double * restrict hor, const double invdR, const double Rmin, const double Rmax,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 const int ndisk, disk_data * restrict disk, double * restrict invSigma, double * restrict sub, double * restrict tmp, double * restrict sum)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (Rmax - Rmin);
  const double pls = 0.5 * (Rmax + Rmin);
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] = 0.0;
  if( NINTBIN & 1 ){
    //---------------------------------------------------------------------
    const double ww =             gsl_gaussQD_weight[(NINTBIN >> 1)];
    const double aa = pls + mns * gsl_gaussQD_pos   [(NINTBIN >> 1)];
    int ia, iR;
    double fa, fR, apR2, amR2;
    _setRdep4GDpot(aa, &ia, &fa, ndisk, disk, invSigma, RR, &iR, &fR, NR, hor, invdR, &apR2, &amR2);
    _gaussQuad2d4GDpot(zz, Nz, ver, invdz, zmin, zmax, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, tmp, invSigma, sub);
    for(int ii = 0; ii < ndisk; ii++)      sum[ii] = ww * tmp[ii];
    //---------------------------------------------------------------------
  }/* if( NINTBIN & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double ww = gsl_gaussQD_weight[ii];
    //---------------------------------------------------------------------
    double aa = pls + mns * gsl_gaussQD_pos[ii];
    int ia, iR;
    double fa, fR, apR2, amR2;
    _setRdep4GDpot(aa, &ia, &fa, ndisk, disk, invSigma, RR, &iR, &fR, NR, hor, invdR, &apR2, &amR2);
    _gaussQuad2d4GDpot(zz, Nz, ver, invdz, zmin, zmax, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, tmp, invSigma, sub);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
    //---------------------------------------------------------------------
    aa = pls - mns * gsl_gaussQD_pos[ii];
    _setRdep4GDpot(aa, &ia, &fa, ndisk, disk, invSigma, RR, &iR, &fR, NR, hor, invdR, &apR2, &amR2);
    _gaussQuad2d4GDpot(zz, Nz, ver, invdz, zmin, zmax, ia, fa, iR, fR, aa, apR2, amR2, ndisk, disk, tmp, invSigma, sub);
    for(int jj = 0; jj < ndisk; jj++)      sum[jj] += ww * tmp[jj];
    //---------------------------------------------------------------------
  }/* for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  /* finalization */
  for(int ii = 0; ii < ndisk; ii++)
    sum[ii] *= mns;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double integrateGDpot
(const double RR, const int NR, double * restrict hor, const double invdR, const double Rs, const double Rmax,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zs, const double zmax,
 const int ndisk, disk_data * restrict disk, double * restrict invSigma, double * restrict sub, double * restrict tmp, double * restrict sum)
{
  //-----------------------------------------------------------------------
  double Phi = 0.0;
  //-----------------------------------------------------------------------
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, 0.0     , 2.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, 2.0 * zs, 6.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, 6.0 * zs,       zmax, ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  //-----------------------------------------------------------------------
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 2.0 * Rs, 8.0 * Rs, zz, Nz, ver, invdz, 0.0     , 2.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 2.0 * Rs, 8.0 * Rs, zz, Nz, ver, invdz, 2.0 * zs, 6.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 2.0 * Rs, 8.0 * Rs, zz, Nz, ver, invdz, 6.0 * zs,       zmax, ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  //-----------------------------------------------------------------------
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 8.0 * Rs, Rmax, zz, Nz, ver, invdz, 0.0     , 2.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 8.0 * Rs, Rmax, zz, Nz, ver, invdz, 2.0 * zs, 6.0 * zs  , ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  gaussQuad2d4GDpot(RR, NR, hor, invdR, 8.0 * Rs, Rmax, zz, Nz, ver, invdz, 6.0 * zs,       zmax, ndisk, disk, invSigma, sub, tmp, sum);
  for(int ii = 0; ii < ndisk; ii++)    Phi += sum[ii];
  //-----------------------------------------------------------------------
  /* exit(0); */
  return (Phi);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#else///USE_GD_FORM_POTENTIAL
//-------------------------------------------------------------------------
#ifndef USE_ELLIPTIC_INTEGRAL
static inline double rinv4pot(double d2, double dd, double phi){  return (1.0 / sqrt(DBL_MIN + d2 + dd * (1.0 - cos(phi))));}
double gaussQD4rinv4pot(const double min, const double max, const double r2, const double Rmean);
double gaussQD4rinv4pot(const double min, const double max, const double r2, const double Rmean)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (max - min);
  const double pls = 0.5 * (max + min);
  double sum = 0.0;
  //-----------------------------------------------------------------------
  if( NINTBIN & 1 )
    sum = gsl_gaussQD_weight[(NINTBIN >> 1)] * rinv4pot(r2, Rmean, pls + mns * gsl_gaussQD_pos[(NINTBIN >> 1)]);
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_weight[ii] * (rinv4pot(r2, Rmean, pls + mns * gsl_gaussQD_pos[ii]) + rinv4pot(r2, Rmean, pls - mns * gsl_gaussQD_pos[ii]));
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
#endif//USE_ELLIPTIC_INTEGRAL
//-------------------------------------------------------------------------
static inline double calcPot
(const int iR, const double fR, const double Rp, const double R2, const double Rmean,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zp,
 double * restrict rho)
{
  //-----------------------------------------------------------------------
  const double zdiff = zz - zp;
#ifdef  USE_ELLIPTIC_INTEGRAL
  const double rinv = 1.0 / sqrt(R2 + zdiff * zdiff);
#endif//USE_ELLIPTIC_INTEGRAL
  //-----------------------------------------------------------------------
  double fz;
  const int jz = bisection(fabs(zp), Nz, ver, false, invdz, &fz);
  //-----------------------------------------------------------------------
  const double ret = Rp *
    (((1.0 - fz) * rho[INDEX2D(NR, Nz,     iR, jz)] + fz * rho[INDEX2D(NR, Nz,     iR, 1 + jz)]) * (1.0 - fR) +
     ((1.0 - fz) * rho[INDEX2D(NR, Nz, 1 + iR, jz)] + fz * rho[INDEX2D(NR, Nz, 1 + iR, 1 + jz)]) *        fR   )
#ifdef  USE_ELLIPTIC_INTEGRAL
    * 2.0 * rinv * gsl_sf_ellint_Kcomp(Rmean * rinv, GSL_PREC_DOUBLE)
#else///USE_ELLIPTIC_INTEGRAL
    * gaussQD4rinv4pot(0.0, M_PI, R2 + zdiff * zdiff, Rmean)
#endif//USE_ELLIPTIC_INTEGRAL
    ;
  //-----------------------------------------------------------------------
  return (ret);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
double _gaussQuad2d4calcPot
(const int iR, const double fR, const double Rp, const double R2, const double Rmean,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 double * restrict rho);
double _gaussQuad2d4calcPot
(const int iR, const double fR, const double Rp, const double R2, const double Rmean,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 double * restrict rho)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (zmax - zmin);
  const double pls = 0.5 * (zmax + zmin);
  double sum = 0.0;
  if( NINTBIN & 1 ){
    //---------------------------------------------------------------------
    const double weight =             gsl_gaussQD_weight[(NINTBIN >> 1)];
    const double     zp = pls + mns * gsl_gaussQD_pos   [(NINTBIN >> 1)];
    sum = weight * calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, zp, rho);
    //---------------------------------------------------------------------
  }/* if( NINTBIN & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double weight = gsl_gaussQD_weight[ii];
    //---------------------------------------------------------------------
    sum += weight * (calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, pls + mns * gsl_gaussQD_pos[ii], rho) +
		     calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, pls - mns * gsl_gaussQD_pos[ii], rho));
    //---------------------------------------------------------------------
  }/* for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void _setRdep4calcPot
(const double Rp, double * restrict R2, int * restrict iR, double * restrict fR, double * restrict Rmean,
 const double RR, const int NR, double * restrict hor, const double invdR)
{
  //-----------------------------------------------------------------------
#ifdef  USE_ELLIPTIC_INTEGRAL
  *Rmean = 2.0 * sqrt(RR * Rp);
  const double Rsum = RR + Rp;  *R2 = Rsum * Rsum;
#else///USE_ELLIPTIC_INTEGRAL
  *Rmean = 2.0 * RR * Rp;
  const double Rsub = RR - Rp;  *R2 = Rsub * Rsub;
#endif//USE_ELLIPTIC_INTEGRAL
  *iR = bisection(Rp, NR, hor, false, invdR, fR);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double gaussQuad2d4calcPot
(const double RR, const int NR, double * restrict hor, const double invdR, const double Rmin, const double Rmax,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zmin, const double zmax,
 double * restrict rho)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (Rmax - Rmin);
  const double pls = 0.5 * (Rmax + Rmin);
  double sum = 0.0;
  if( NINTBIN & 1 ){
    //---------------------------------------------------------------------
    const double weight =             gsl_gaussQD_weight[(NINTBIN >> 1)];
    const double     Rp = pls + mns * gsl_gaussQD_pos   [(NINTBIN >> 1)];
    double R2, fR, Rmean;
    int iR;
    _setRdep4calcPot(Rp, &R2, &iR, &fR, &Rmean, RR, NR, hor, invdR);
    sum = weight * _gaussQuad2d4calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, zmin, zmax, rho);
    //---------------------------------------------------------------------
  }/* if( NINTBIN & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    const double weight = gsl_gaussQD_weight[ii];
    //---------------------------------------------------------------------
    double Rp = pls + mns * gsl_gaussQD_pos[ii];
    double R2, fR, Rmean;
    int iR;
    _setRdep4calcPot(Rp, &R2, &iR, &fR, &Rmean, RR, NR, hor, invdR);
    sum += weight * _gaussQuad2d4calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, zmin, zmax, rho);
    //---------------------------------------------------------------------
    Rp = pls - mns * gsl_gaussQD_pos[ii];
    _setRdep4calcPot(Rp, &R2, &iR, &fR, &Rmean, RR, NR, hor, invdR);
    sum += weight * _gaussQuad2d4calcPot(iR, fR, Rp, R2, Rmean, zz, Nz, ver, invdz, zmin, zmax, rho);
    //---------------------------------------------------------------------
  }/* for(int ii = (NINTBIN >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  /* finalization */
  return (mns * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double integratePotential
(const double RR, const int NR, double * restrict hor, const double invdR, const double Rs, const double Rmax,
 const double zz, const int Nz, double * restrict ver, const double invdz, const double zs, const double zmax, double * restrict rho)
{
  //-----------------------------------------------------------------------
  double Phi = 0.0;
  //-----------------------------------------------------------------------
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, -zmax, -zs  , rho);
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz, -zs  ,  zs  , rho);
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 0.0, 2.0 * Rs, zz, Nz, ver, invdz,  zs  ,  zmax, rho);
  //-----------------------------------------------------------------------
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 2.0 * Rs, Rmax, zz, Nz, ver, invdz, -zmax, -zs  , rho);
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 2.0 * Rs, Rmax, zz, Nz, ver, invdz, -zs  ,  zs  , rho);
  Phi += gaussQuad2d4calcPot(RR, NR, hor, invdR, 2.0 * Rs, Rmax, zz, Nz, ver, invdz,  zs  ,  zmax, rho);
  //-----------------------------------------------------------------------
  return (Phi);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#endif//USE_GD_FORM_POTENTIAL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void calcOuterPotential
(double * restrict hor, const double Rmax, double * restrict ver, const double zmax,
 const double hh, const double invhh,
 const double Rs, const double zs, double * restrict Phi_NR, double * restrict Phi_Nz,
#ifdef  USE_GD_FORM_POTENTIAL
 const int ndisk, disk_data * restrict disk, double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum
#else///USE_GD_FORM_POTENTIAL
 double * restrict rho
#endif//USE_GD_FORM_POTENTIAL
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* \Phi_{  i, N_z} = potential @ R = (  i + 1/2) * R_max / N_R, z = (N_z + 1/2) * z_max / N_z; for i = 0, ..., N_R - 1 */
  //-----------------------------------------------------------------------
#pragma omp parallel for num_threads(CPUS_PER_PROCESS)
  for(int jj = 0; jj < NDISKBIN_VER; jj++){
    //---------------------------------------------------------------------
    const double RR = hor[NDISKBIN_HOR - 1] + hh;
    const double zz = ver[jj];
    //---------------------------------------------------------------------
#ifdef  USE_GD_FORM_POTENTIAL
    const int tidx = omp_get_thread_num();
    double *inv = &stock_inv[ndisk * tidx];
    double *sub = &stock_sub[ndisk * tidx];
    double *tmp = &stock_tmp[ndisk * tidx];
    double *sum = &stock_sum[ndisk * tidx];
    Phi_NR[jj] =  4.0 * (double)newton * integrateGDpot    (RR, NDISKBIN_HOR, hor, invhh, Rs, Rmax, zz, NDISKBIN_VER, ver, invhh, zs, zmax, ndisk, disk, inv, sub, tmp, sum);
#else///USE_GD_FORM_POTENTIAL
    Phi_NR[jj] = -2.0 * (double)newton * integratePotential(RR, NDISKBIN_HOR, hor, invhh, Rs, Rmax, zz, NDISKBIN_VER, ver, invhh, zs, zmax, rho);
#endif//USE_GD_FORM_POTENTIAL
    //---------------------------------------------------------------------
  }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* \Phi_{N_R,   j} = potential @ R = (N_R + 1/2) * R_max / N_R, z = (  j + 1/2) * z_max / N_z; for j = 0, ..., N_z - 1 */
  //-----------------------------------------------------------------------
#pragma omp parallel for num_threads(CPUS_PER_PROCESS)
  for(int ii = 0; ii < NDISKBIN_HOR; ii++){
    //---------------------------------------------------------------------
    const double zz = ver[NDISKBIN_VER - 1] + hh;
    const double RR = hor[ii];
    //---------------------------------------------------------------------
#ifdef  USE_GD_FORM_POTENTIAL
    const int tidx = omp_get_thread_num();
    double *inv = &stock_inv[ndisk * tidx];
    double *sub = &stock_sub[ndisk * tidx];
    double *tmp = &stock_tmp[ndisk * tidx];
    double *sum = &stock_sum[ndisk * tidx];
    Phi_Nz[ii] =  4.0 * (double)newton * integrateGDpot    (RR, NDISKBIN_HOR, hor, invhh, Rs, Rmax, zz, NDISKBIN_VER, ver, invhh, zs, zmax, ndisk, disk, inv, sub, tmp, sum);
#else///USE_GD_FORM_POTENTIAL
    Phi_Nz[ii] = -2.0 * (double)newton * integratePotential(RR, NDISKBIN_HOR, hor, invhh, Rs, Rmax, zz, NDISKBIN_VER, ver, invhh, zs, zmax, rho);
#endif//USE_GD_FORM_POTENTIAL
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void setResultantVector(const int lev, const double hh, double * restrict RR, double * restrict rho, double * restrict Phi_NR, double * restrict Phi_Nz, double * restrict vec)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  const double Phi0 = 8.0 * M_PI * (double)newton * hh;
#pragma omp parallel for
  for(int ii = 0; ii < NDISKBIN_HOR; ii++){
    //---------------------------------------------------------------------
    const double Ri = RR[ii];
    //---------------------------------------------------------------------
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
      vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] = Phi0 * Ri * rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)];
    //---------------------------------------------------------------------
    /* boundary condition @ z-edge */
    vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)] -= (double)(1 + (ii << 1)) * Phi_Nz[ii];
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
  //-----------------------------------------------------------------------
  /* boundary condition @ R-edge */
  for(int jj = 0; jj < NDISKBIN_VER; jj++)
    vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] -= (double)(NDISKBIN_HOR << 1) * Phi_NR[jj];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void setSparseMatrix(const crs mat)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* const int Nrow = NDISKBIN_HOR * NDISKBIN_VER; */
  int valIdx = 0;
  int rowIdx = 0;
  mat.row[rowIdx] = valIdx;
  //-----------------------------------------------------------------------
/*   const double Phi0 = 8.0 * M_PI * (double)newton * hh; */
/* #pragma omp parallel for */
/*   for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
/*     //--------------------------------------------------------------------- */
/*     const double Ri = RR[ii]; */
/*     //--------------------------------------------------------------------- */
/*     for(int jj = 0; jj < NDISKBIN_VER; jj++) */
/*       vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] = Phi0 * Ri * rho[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)]; */
/*     //--------------------------------------------------------------------- */
/*   }/\* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ *\/ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* inner boundary (ii = 0) */
  //-----------------------------------------------------------------------
  /* ii = 0; jj = 0; */
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, 0);  mat.val[valIdx] = -3.0;  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, 1);  mat.val[valIdx] =  1.0;  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, 0);  mat.val[valIdx] =  2.0;  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
  //-----------------------------------------------------------------------
  /* ii = 0; 1 <= jj <= NDISKBIN_VER - 2 */
  for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){
    //---------------------------------------------------------------------
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj - 1);    mat.val[valIdx] =  1.0;    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj    );    mat.val[valIdx] = -4.0;    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, jj + 1);    mat.val[valIdx] =  1.0;    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, jj    );    mat.val[valIdx] =  2.0;    valIdx++;
    //---------------------------------------------------------------------
    rowIdx++;    mat.row[rowIdx] = valIdx;
    //---------------------------------------------------------------------
  }/* for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){ */
  //-----------------------------------------------------------------------
  /* ii = 0; jj = NDISKBIN_VER - 1 */
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, NDISKBIN_VER - 2);  mat.val[valIdx] =  1.0;  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, NDISKBIN_VER - 1);  mat.val[valIdx] = -4.0;  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 1, NDISKBIN_VER - 1);  mat.val[valIdx] =  2.0;  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
  /* vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, 0, NDISKBIN_VER - 1)] -= Phi_Nz[0]; */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* main domain (1 <= ii <= NDISKBIN_HOR - 2) */
  //-----------------------------------------------------------------------
  for(int ii = 1; ii < NDISKBIN_HOR - 1; ii++){
    //---------------------------------------------------------------------
    /* inner boundary (jj = 0) */
    //---------------------------------------------------------------------
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, 0);    mat.val[valIdx] =  (double)(      ii << 1      );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , 0);    mat.val[valIdx] = -(double)((1 + (ii << 1)) * 3);    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , 1);    mat.val[valIdx] =  (double)( 1 + (ii << 1)     );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, 0);    mat.val[valIdx] =  (double)((1 +  ii     ) << 1);    valIdx++;
    //---------------------------------------------------------------------
    rowIdx++;    mat.row[rowIdx] = valIdx;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* main domain (1 <= jj <= Nz - 2) */
    //---------------------------------------------------------------------
    for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){
      //-------------------------------------------------------------------
      mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, jj    );      mat.val[valIdx] =  (double)(      ii << 1       );      valIdx++;
      mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , jj - 1);      mat.val[valIdx] =  (double)( 1 + (ii << 1)      );      valIdx++;
      mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , jj    );      mat.val[valIdx] = -(double)((1 + (ii << 1)) << 2);      valIdx++;
      mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , jj + 1);      mat.val[valIdx] =  (double)( 1 + (ii << 1)      );      valIdx++;
      mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, jj    );      mat.val[valIdx] =  (double)((1 +  ii      ) << 1);      valIdx++;
      //-------------------------------------------------------------------
      rowIdx++;      mat.row[rowIdx] = valIdx;
      //-------------------------------------------------------------------
    }/* for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){ */
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* outer boundary (jj = Nz - 1) */
    //---------------------------------------------------------------------
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii - 1, NDISKBIN_VER - 1);    mat.val[valIdx] =  (double)(      ii << 1       );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , NDISKBIN_VER - 2);    mat.val[valIdx] =  (double)( 1 + (ii << 1)      );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii    , NDISKBIN_VER - 1);    mat.val[valIdx] = -(double)((1 + (ii << 1)) << 2);    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii + 1, NDISKBIN_VER - 1);    mat.val[valIdx] =  (double)((1 +  ii      ) << 1);    valIdx++;
    //---------------------------------------------------------------------
    rowIdx++;    mat.row[rowIdx] = valIdx;
    //---------------------------------------------------------------------
    /* vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)] -= (double)(1 + (ii << 1)) * Phi_Nz[ii]; */
    //---------------------------------------------------------------------
  }/* for(int ii = 1; ii < NDISKBIN_HOR - 1; ii++){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* outer boundary (ii = NDISKBIN_HOR - 1) */
  //-----------------------------------------------------------------------
  /* ii = NDISKBIN_HOR - 1; jj = 0; */
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, 0);  mat.val[valIdx] =  (double)(      (NDISKBIN_HOR - 1) << 1      );  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, 0);  mat.val[valIdx] = -(double)((1 + ((NDISKBIN_HOR - 1) << 1)) * 3);  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, 1);  mat.val[valIdx] =  (double)( 1 + ((NDISKBIN_HOR - 1) << 1)     );  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
  /* vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, 0)] -= (double)(NDISKBIN_HOR << 1) * Phi_NR[0]; */
  //-----------------------------------------------------------------------
  /* ii = NDISKBIN_HOR - 1; 1 <= jj <= NDISKBIN_VER - 2 */
  for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){
    //---------------------------------------------------------------------
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, jj    );    mat.val[valIdx] =  (double)(      (NDISKBIN_HOR - 1) << 1       );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj - 1);    mat.val[valIdx] =  (double)( 1 + ((NDISKBIN_HOR - 1) << 1)      );    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj    );    mat.val[valIdx] = -(double)((1 + ((NDISKBIN_HOR - 1) << 1)) << 2);    valIdx++;
    mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj + 1);    mat.val[valIdx] =  (double)( 1 + ((NDISKBIN_HOR - 1) << 1)      );    valIdx++;
    //---------------------------------------------------------------------
    rowIdx++;    mat.row[rowIdx] = valIdx;
    //---------------------------------------------------------------------
    /* vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, jj)] -= (double)(NDISKBIN_HOR << 1) * Phi_NR[jj]; */
    //---------------------------------------------------------------------
  }/* for(int jj = 1; jj < NDISKBIN_VER - 1; jj++){ */
  //-----------------------------------------------------------------------
  /* ii = NR - 1; jj = NDISKBIN_VER - 1 */
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 2, NDISKBIN_VER - 1);  mat.val[valIdx] =  (double)(      (NDISKBIN_HOR - 1 ) << 1       );  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, NDISKBIN_VER - 2);  mat.val[valIdx] =  (double)( 1 + ((NDISKBIN_HOR - 1 ) << 1)      );  valIdx++;
  mat.col[valIdx] = INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, NDISKBIN_VER - 1);  mat.val[valIdx] = -(double)((1 + ((NDISKBIN_HOR - 1 ) << 1)) << 2);  valIdx++;
  rowIdx++;  mat.row[rowIdx] = valIdx;
  /* vec[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, NDISKBIN_HOR - 1, NDISKBIN_VER - 1)] -= */
  /*   (double)(NDISKBIN_HOR << 1) * Phi_NR[NDISKBIN_VER - 1] + (double)(1 + ((NDISKBIN_HOR - 1) << 1)) * Phi_Nz[NDISKBIN_HOR - 1]; */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  if( mat.row[rowIdx] != (5 * NDISKBIN_HOR * NDISKBIN_VER - 2 * (NDISKBIN_HOR + NDISKBIN_VER)) ){
    //---------------------------------------------------------------------
    __KILL__(stderr, "ERROR: Number of non-zero elements is %d, while the expected is %d\n", mat.row[rowIdx], 5 * NDISKBIN_HOR * NDISKBIN_VER - 2 * (NDISKBIN_HOR + NDISKBIN_VER));
    //---------------------------------------------------------------------
  }/* if( mat.row[rowIdx] != (5 * NDISKBIN_HOR * NDISKBIN_VER - 2 * (NDISKBIN_HOR + NDISKBIN_VER)) ){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define MONOTONIC_DIFF(ll, cc, rr) ((((cc) - (ll)) * ((rr) - (cc)) > 0.0) ? ((rr) - (ll)) : (0.0)) */
//-------------------------------------------------------------------------
static inline void getPotentialField
(const int ndisk, disk_data * restrict disk, const int lev,
 const double hh, const double invhh, double * restrict RR, double * restrict zz,
 double * restrict rho, double * restrict Phi, double * restrict Phi_NR, double * restrict Phi_Nz, soaBiCGSTAB mat, soaPreConditioning pre
#ifdef  USE_GD_FORM_POTENTIAL
 , double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum
#endif//USE_GD_FORM_POTENTIAL
 )
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* prepare Poisson equation in matrix form (vector part) */
  //-----------------------------------------------------------------------
  if( lev > 0 ){
    //---------------------------------------------------------------------
    coarsePotential4boundary(lev - 1, disk[0].pot);
    //---------------------------------------------------------------------
    /* set z-edge[NR] */
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      //-------------------------------------------------------------------
      const double mm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii > 1) ? ((ii >> 1) - 1) : (0), (NDISKBIN_VER >> 1) - 1)];
      const double mc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii > 1) ? ((ii >> 1) - 1) : (0),  NDISKBIN_VER >> 1)	 ];
      const double mp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii > 1) ? ((ii >> 1) - 1) : (0), (NDISKBIN_VER >> 1) + 1)];
      //-------------------------------------------------------------------
      const double cm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		    ii >> 1	       , (NDISKBIN_VER >> 1) - 1)];
      const double cc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		    ii >> 1	       ,  NDISKBIN_VER >> 1)	 ];
      const double cp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		    ii >> 1	       , (NDISKBIN_VER >> 1) + 1)];
      //-------------------------------------------------------------------
      const double pm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		   (ii >> 1) + 1       , (NDISKBIN_VER >> 1) - 1)];
      const double pc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		   (ii >> 1) + 1       ,  NDISKBIN_VER >> 1)	 ];
      const double pp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,		   (ii >> 1) + 1       , (NDISKBIN_VER >> 1) + 1)];
      //-------------------------------------------------------------------
      const double zmm = mc - 0.125 * (mp - mm);
      const double zcm = cc - 0.125 * (cp - cm);
      const double zpm = pc - 0.125 * (pp - pm);
      //-------------------------------------------------------------------
      Phi_Nz[ii] = zcm + (double)(((ii & 1) << 1) - 1) * 0.125 * (zpm - zmm);
      //-------------------------------------------------------------------
#if 0
      fprintf(stderr, "%d\t%e\t%e\t%e\t%e\n", ii, Phi_Nz[ii], zmm, zcm, zpm);
#endif
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
    /* Phi_Nz[0] = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, 0, NDISKBIN_VER >> 1)]; */
    /* for(int ii = 1; ii < NDISKBIN_HOR; ii++) */
    /*   Phi_Nz[ii] = */
    /* 	0.75 * disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii    ) >> 1, NDISKBIN_VER >> 1)] + */
    /* 	0.25 * disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (ii + 1) >> 1, NDISKBIN_VER >> 1)]; */
    //---------------------------------------------------------------------
    /* set R-edge[Nz] */
#pragma omp parallel for
    for(int jj = 0; jj < NDISKBIN_VER; jj++){
      //-------------------------------------------------------------------
      const double mm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) - 1, (jj > 1 ) ? ((jj >> 1) - 1) : (0))];
      const double cm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,	NDISKBIN_HOR >> 1     , (jj > 1 ) ? ((jj >> 1) - 1) : (0))];
      const double pm = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) + 1, (jj > 1 ) ? ((jj >> 1) - 1) : (0))];
      //-------------------------------------------------------------------
      const double mc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) - 1,		      jj >> 1		 )];
      const double cc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,	NDISKBIN_HOR >> 1     ,		      jj >> 1		 )];
      const double pc = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) + 1,		      jj >> 1		 )];
      //-------------------------------------------------------------------
      const double mp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) - 1,		     (jj >> 1) + 1	 )];
      const double cp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,	NDISKBIN_HOR >> 1     ,		     (jj >> 1) + 1	 )];
      const double pp = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) + 1,		     (jj >> 1) + 1	 )];
      //-------------------------------------------------------------------
      const double Rmm = cm - 0.125 * (pm - mm);
      const double Rmc = cc - 0.125 * (pc - mc);
      const double Rmp = cp - 0.125 * (pp - mp);
      //-------------------------------------------------------------------
      Phi_NR[jj] = Rmc + (double)(((jj & 1) << 1) - 1) * 0.125 * (Rmp - Rmm);
      //-------------------------------------------------------------------
#if 0
      fprintf(stderr, "%d\t%e\t%e\t%e\t%e\n", jj, Phi_NR[jj], Rmm, Rmc, Rmp);
#endif
      //-------------------------------------------------------------------
      /* const double zmm = mc - 0.125 * (mp - mm);      const double zcm = cc - 0.125 * (cp - cm);      const double zpm = pc - 0.125 * (pp - pm); */
      /* const double zmp = mc + 0.125 * (mp - mm);      const double zcp = cc + 0.125 * (cp - cm);      const double zpp = pc + 0.125 * (pp - pm); */
      //-------------------------------------------------------------------
      /* Phi_NR[jj    ] = zcm - 0.125 * (zpm - zmm); */
      /* Phi_NR[jj + 1] = zcp - 0.125 * (zpp - zmp); */
      //-------------------------------------------------------------------
    }/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
    /* Phi_NR[0] = disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, NDISKBIN_HOR >> 1, 0)]; */
    /* for(int jj = 1; jj < NDISKBIN_VER; jj++) */
    /*   Phi_NR[jj] = */
    /* 	0.75 * disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, NDISKBIN_HOR >> 1, (jj	) >> 1)] + */
    /* 	0.25 * disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, NDISKBIN_HOR >> 1, (jj + 1) >> 1)]; */
    //---------------------------------------------------------------------
#if 0
    fprintf(stderr, "# Phi_Nz[%d]\n", NDISKBIN_HOR);
    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
      fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", RR[ii], Phi_Nz[ii],
    	      disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, ii >> 1, (NDISKBIN_VER >> 1) - 1)],
    	      disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, ii >> 1,  NDISKBIN_VER >> 1     )],
    	      disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, ii >> 1, (NDISKBIN_VER >> 1) + 1)]);
    fprintf(stderr, "\n");

    fprintf(stderr, "# Phi_NR[%d]\n", NDISKBIN_VER);
    for(int jj = 0; jj < NDISKBIN_VER; jj++)
      fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", zz[jj], Phi_NR[jj],
    	      disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) - 1, jj >> 1)],
    	      disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1,  NDISKBIN_HOR >> 1     , jj >> 1)],
    	      disk[0].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev - 1, (NDISKBIN_HOR >> 1) + 1, jj >> 1)]);
    fflush(NULL);
    exit(0);
#endif
    //---------------------------------------------------------------------
  }/* if( lev > 0 ){ */
  else{
    //---------------------------------------------------------------------
    double Mtot = 0.0;
    double Rs   = 0.0;
    double zs   = 0.0;
    for(int ii = 0; ii < ndisk; ii++){
      //---------------------------------------------------------------------
      Mtot += disk[ii].cfg->Mtot;
      Rs   += disk[ii].cfg->Mtot * disk[ii].cfg->rs;
      zs   += disk[ii].cfg->Mtot * disk[ii].cfg->zd;
      //---------------------------------------------------------------------
    }/* for(int ii = 0; ii < ndisk; ii++){ */
    const double Minv = 1.0 / Mtot;
    Rs *= Minv;
    zs *= Minv;
    //-----------------------------------------------------------------------
    calcOuterPotential(RR, disk[0].Rmax, zz, disk[0].zmax, hh, invhh, Rs, zs, Phi_NR, Phi_Nz,
#ifdef  USE_GD_FORM_POTENTIAL
		       ndisk, disk, stock_inv, stock_sub, stock_tmp, stock_sum
#else///USE_GD_FORM_POTENTIAL
		       rho
#endif//USE_GD_FORM_POTENTIAL
#ifdef  CONFIRM_BUILDING_BLOCK
		       , Mtot
#endif//CONFIRM_BUILDING_BLOCK
		       );
    //---------------------------------------------------------------------
  }/* else{ */
  //-----------------------------------------------------------------------
  setResultantVector(lev, hh, RR, rho, Phi_NR, Phi_Nz, mat.vec);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* solve Poisson equation using iterative method */
  //-----------------------------------------------------------------------
  pbicgstab(mat.mat, NDISKBIN_HOR * NDISKBIN_VER, mat.vec, Phi, mat.res, mat.sdw, mat.mid, mat.tmp, mat.Api, mat.Ati, pre.ilu, pre.Kri, pre.Kpi, pre.Kti, CONVERGENCE_BICGSTAB);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline double _setzdep4calcRho
(const double zz, const int Nz, double * restrict ver, const double invdz,
 double * restrict Phi, profile * restrict sph, const double invlogrbin_sph,
 const int iR, const double fR, const double R2, const double PsiR0, const double invPsi)
{
  //-----------------------------------------------------------------------
  double fz;
  const int jz = bisection(zz, Nz, ver, false, invdz, &fz);
  //-----------------------------------------------------------------------
  const double PsiRz =
    ((1.0 - fz) * (-Phi[INDEX2D(NR, Nz,     iR, jz)]) + fz * (-Phi[INDEX2D(NR, Nz,     iR, 1 + jz)])) * (1.0 - fR) +
    ((1.0 - fz) * (-Phi[INDEX2D(NR, Nz, 1 + iR, jz)]) + fz * (-Phi[INDEX2D(NR, Nz, 1 + iR, 1 + jz)])) *        fR  +
    Phi_spherical(sqrt(R2 + zz * zz), sph, invlogrbin_sph);
  //-----------------------------------------------------------------------
  return (exp(-(PsiRz - PsiR0) * invPsi));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void _setRdep4calcRho
(const double RR, double * restrict R2, double * restrict horDep, int * restrict iR, double * restrict fR, double * restrict PsiR0, double * restrict invPsi,
 const int NR, double * restrict horRad, const double invdR, const int Nz, double * restrict ver, const double invdz, double * restrict Phi,
 const int lev, const disk_data disk, profile * restrict sph, const double invlogrbin_sph)
{
  //-----------------------------------------------------------------------
  *R2 = RR * RR;
  //-----------------------------------------------------------------------
  *iR = bisection(RR, NR, horRad, false, invdR, fR);
  *PsiR0 = (1.0 - (*fR)) * (-Phi[INDEX2D(NR, Nz, *iR, 0)]) + (*fR) * (-Phi[INDEX2D(NR, Nz, 1 + (*iR), 0)]) + Phi_spherical(RR, sph, invlogrbin_sph);
  //-----------------------------------------------------------------------
#ifdef  ENABLE_VARIABLE_SCALE_HEIGHT
  const double    zd = (1.0 - (*fR)) * disk.zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, *iR)] + (*fR) * disk.zd[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + (*iR))];
#else///ENABLE_VARIABLE_SCALE_HEIGHT
  const double    zd = disk.cfg->zd;
#endif//ENABLE_VARIABLE_SCALE_HEIGHT
  const double invzd = 1.0 / zd;
  //-----------------------------------------------------------------------
  /* *horDep = invzd * ((1.0 - (*fR)) * disk.Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, *iR)] + (*fR) * disk.Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, 1 + (*iR))]); */
  *horDep = invzd * disk.getColumnDensity(RR, disk.invRd, disk.util);
  //-----------------------------------------------------------------------
  double fz;
  const int jz = bisection(zd, Nz, ver, false, invdz, &fz);
  const double Psi_Rzd =
    ((1.0 - fz) * (-Phi[INDEX2D(NR, Nz,      *iR , jz)]) + fz * (-Phi[INDEX2D(NR, Nz,      *iR , 1 + jz)])) * (1.0 - (*fR)) +
    ((1.0 - fz) * (-Phi[INDEX2D(NR, Nz, 1 + (*iR), jz)]) + fz * (-Phi[INDEX2D(NR, Nz, 1 + (*iR), 1 + jz)])) *        (*fR)  +
    Phi_spherical(sqrt((*R2) + zd * zd), sph, invlogrbin_sph);
  *invPsi = 1.0 / (Psi_Rzd - (*PsiR0));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double _gaussQuad2d4calcRho
(const double zmin, const double zmax, const int Nz, double * restrict ver, const double invdz, const int iR, const double fR, const double R2,
 const double PsiR0, const double invPsi, double * restrict Phi, profile * restrict sph, const double invlogrbin_sph)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (zmax - zmin);
  const double pls = 0.5 * (zmax + zmin);
  double sum = 0.0;
  if( NINTBIN_LOW & 1 ){
    //---------------------------------------------------------------------
    sum = gsl_gaussQD_low_weight[(NINTBIN_LOW >> 1)] *
      _setzdep4calcRho(pls + mns * gsl_gaussQD_low_pos[(NINTBIN_LOW >> 1)], Nz, ver, invdz, Phi, sph, invlogrbin_sph, iR, fR, R2, PsiR0, invPsi);
    //---------------------------------------------------------------------
  }/* if( NINTBIN_LOW & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--)
    sum += gsl_gaussQD_low_weight[ii] *
      (_setzdep4calcRho(pls + mns * gsl_gaussQD_low_pos[ii], Nz, ver, invdz, Phi, sph, invlogrbin_sph, iR, fR, R2, PsiR0, invPsi) +
       _setzdep4calcRho(pls - mns * gsl_gaussQD_low_pos[ii], Nz, ver, invdz, Phi, sph, invlogrbin_sph, iR, fR, R2, PsiR0, invPsi));
  //-----------------------------------------------------------------------
  /* finalization */
  return (sum * mns);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline double gaussQuad2d4calcRho
(const double Rmin, const double Rmax, const int NR, double * restrict hor, const double invdR,
 const double zmin, const double zmax, const int Nz, double * restrict ver, const double invdz,
 profile * restrict sph, const double invlogrbin_sph, double * restrict Phi, const int lev, const disk_data disk)
{
  //-----------------------------------------------------------------------
  /* initialization */
  const double mns = 0.5 * (Rmax - Rmin);
  const double pls = 0.5 * (Rmax + Rmin);
  double sum = 0.0;
  if( NINTBIN_LOW & 1 ){
    //---------------------------------------------------------------------
    double R2, radVal, fR, PsiR0, invPsi;
    int iR;
    _setRdep4calcRho(pls + mns * gsl_gaussQD_low_pos[(NINTBIN_LOW >> 1)],
		     &R2, &radVal, &iR, &fR, &PsiR0, &invPsi, NR, hor, invdR, Nz, ver, invdz, Phi, lev, disk, sph, invlogrbin_sph);
    //---------------------------------------------------------------------
    sum = gsl_gaussQD_low_weight[(NINTBIN_LOW >> 1)] * radVal * _gaussQuad2d4calcRho(zmin, zmax, Nz, ver, invdz, iR, fR, R2, PsiR0, invPsi, Phi, sph, invlogrbin_sph);
    //---------------------------------------------------------------------
  }/* if( NINTBIN_LOW & 1 ){ */
  //-----------------------------------------------------------------------
  /* numerical integration */
  for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--){
    //---------------------------------------------------------------------
    double R2, radVal, fR, PsiR0, invPsi;
    int iR;
    //---------------------------------------------------------------------
    _setRdep4calcRho(pls + mns * gsl_gaussQD_low_pos[ii], &R2, &radVal, &iR, &fR, &PsiR0, &invPsi, NR, hor, invdR, Nz, ver, invdz, Phi, lev, disk, sph, invlogrbin_sph);
    sum += gsl_gaussQD_low_weight[ii] * radVal * _gaussQuad2d4calcRho(zmin, zmax, Nz, ver, invdz, iR, fR, R2, PsiR0, invPsi, Phi, sph, invlogrbin_sph);
    //---------------------------------------------------------------------
    _setRdep4calcRho(pls - mns * gsl_gaussQD_low_pos[ii], &R2, &radVal, &iR, &fR, &PsiR0, &invPsi, NR, hor, invdR, Nz, ver, invdz, Phi, lev, disk, sph, invlogrbin_sph);
    sum += gsl_gaussQD_low_weight[ii] * radVal * _gaussQuad2d4calcRho(zmin, zmax, Nz, ver, invdz, iR, fR, R2, PsiR0, invPsi, Phi, sph, invlogrbin_sph);
    //---------------------------------------------------------------------
  }/* for(int ii = (NINTBIN_LOW >> 1) - 1; ii >= 0; ii--){ */
  //-----------------------------------------------------------------------
  /* finalization */
  return (sum * mns);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void swapDblArrays(double **p0, double **p1)
{
  //-----------------------------------------------------------------------
  double *tmp;
  //-----------------------------------------------------------------------
  tmp = *p0;
  *p0 = *p1;
  *p1 = tmp;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void getPotDensPair
(const int ndisk, const int lev, const int levOld, disk_data * restrict disk,
 double * restrict Phi_NR, double * restrict Phi_Nz,
 soaBiCGSTAB mat, soaPreConditioning pre
#ifdef  USE_GD_FORM_POTENTIAL
 , double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum
#endif//USE_GD_FORM_POTENTIAL
 );
void getPotDensPair
(const int ndisk, const int lev, const int levOld, disk_data * restrict disk,
 double * restrict Phi_NR, double * restrict Phi_Nz,
 soaBiCGSTAB mat, soaPreConditioning pre
#ifdef  USE_GD_FORM_POTENTIAL
 , double * restrict stock_inv, double * restrict stock_sub, double * restrict stock_tmp, double * restrict stock_sum
#endif//USE_GD_FORM_POTENTIAL
)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* load disk data (common settings for all components) */
  //-----------------------------------------------------------------------
  double * RR;  RR  = &(disk[0].hor[INDEX2D(maxLev, NDISKBIN_HOR               , lev, 0)]);
  double * zz;  zz  = &(disk[0].ver[INDEX2D(maxLev,                NDISKBIN_VER, lev, 0)]);
  double *Phi;  Phi = &(disk[0].pot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, 0)]);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* load spherical components(s) data (common settings for all components) */
  //-----------------------------------------------------------------------
  profile *sph;  sph = disk[0].prf;
  const double invlogrbin_sph = disk[0].invlogrbin;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set grid points */
  //-----------------------------------------------------------------------
  /* for(int ii = 0; ii < ndisk; ii++){ */
  /*   //--------------------------------------------------------------------- */
  /*   disk[ii].Rbin = Rbin;    disk[ii].invRbin = invRbin; */
  /*   disk[ii].zbin = zbin;    disk[ii].invzbin = invzbin; */
  /*   //--------------------------------------------------------------------- */
  /* }/\* for(int ii = 0; ii < ndisk; ii++){ *\/ */
  const double    hh = ldexp(disk[0].hh, -lev);
  const double invhh = 1.0 / hh;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set density field and guess potential field */
  //-----------------------------------------------------------------------
  if( lev != levOld ){
    //---------------------------------------------------------------------
    if( lev > levOld )      for(int ii = 0; ii < ndisk; ii++){	  fineDensity(lev, disk[ii], *(disk[ii].rho));	  finePotential(lev, disk[ii].pot);      }
    if( lev < levOld )      for(int ii = 0; ii < ndisk; ii++){	coarseDensity(lev, disk[ii], *(disk[ii].rho));	coarsePotential(lev, disk[ii].pot);      }
    //---------------------------------------------------------------------
  }/* if( levOld != KICKOFF_POISSON_SOLVER ){ */
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* iterative process to get the potential-density pair */
  //-----------------------------------------------------------------------
  /* modify column density profile if necessary */
  if( lev > 0 ){
    //---------------------------------------------------------------------
    /* save the column density profile */
    //---------------------------------------------------------------------
    for(int ii = 0; ii < ndisk; ii++)
#pragma omp parallel for
      for(int jj = 0; jj < NDISKBIN_HOR; jj++)
	disk[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] = disk[ii].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)];
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* subtract mass above the domain in this level */
    //---------------------------------------------------------------------
    for(int ll = lev - 1; ll >= 0; ll--){
      const int ldiff = lev - ll;
      const double hh_tmp = ldexp(disk[0].hh, -ll);
#pragma omp parallel for
      for(int ii = 0; ii < (NDISKBIN_HOR >> ldiff); ii++)
	for(int kk = 0; kk < ndisk; kk++){
	  /* evaluate mass to be subtracted using Simpson's rule */
	  double mass = 0.0;
	  for(int jj = (NDISKBIN_VER >> 1); jj < NDISKBIN_VER; jj++)
	    mass += (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)];
	  mass *= 2.0 * hh_tmp;/* 2.0 reflects plane symmetry about the equatorial plane (z = 0) */
	  for(int mm = (ii << ldiff); mm < (ii << ldiff) + (1 << ldiff); mm++){
	    disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, mm)] -= mass;
	  }/* for(int mm = 0; mm < (1 << ldiff); mm++){ */
	}/* for(int kk = 0; kk < ndisk; kk++){ */
    }/* for(int ll = lev - 1; ll >= 0; ll--){ */
    //---------------------------------------------------------------------
#if 0
    if( lev == 9 ){
      for(int ii = 0; ii < NDISKBIN_HOR; ii++)
	fprintf(stderr, "%e\t%e\t%e\n", disk[0].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], disk[0].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)], disk[0].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);
      fflush(NULL);
      exit(0);
    }
#endif
    //---------------------------------------------------------------------
  }/* if( lev > 0 ){ */
  //-----------------------------------------------------------------------
#if 0
  for(int ii = 0; ii < NDISKBIN_HOR; ii++){
    fprintf(stderr, "%e", disk[0].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);
    for(int kk = 0; kk < ndisk; kk++)
      fprintf(stderr, "\t%e", disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  fflush(NULL);
  /* exit(0); */
#endif
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# procedure start: NR = %d, Nz = %d\n", NDISKBIN_HOR, NDISKBIN_VER);
  fflush(stdout);
  int steps = 0;
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    /* calculate potential field from the given density field */
    //---------------------------------------------------------------------
    /* set density field */
    double *rhoTot;
    if( ndisk > 1 ){
      //-------------------------------------------------------------------
      rhoTot = disk[0].rhoTot;
      //-------------------------------------------------------------------
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR * NDISKBIN_VER; ii++)
	rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)] = 0.0;
      for(int kk = 0; kk < ndisk; kk++)
#pragma omp parallel for
	for(int ii = 0; ii < NDISKBIN_HOR * NDISKBIN_VER; ii++)
	  rhoTot[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)] += (*disk[kk].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, ii)];
      //-------------------------------------------------------------------
    }/* if( ndisk > 1 ){ */
    else      rhoTot = *disk[0].rho;
    //---------------------------------------------------------------------
    /* solve Poisson equation using ILU(0) preconditioned BiCGSTAB method */
    getPotentialField(ndisk, disk, lev, hh, invhh, RR, zz, rhoTot, Phi, Phi_NR, Phi_Nz, mat, pre
#ifdef  USE_GD_FORM_POTENTIAL
		      , stock_inv, stock_sub, stock_tmp, stock_sum
#endif//USE_GD_FORM_POTENTIAL
		      );
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* update density field from the derived potential field */
    //---------------------------------------------------------------------
    static double errMax;
    errMax = 0.0;
    //---------------------------------------------------------------------
    for(int kk = 0; kk < ndisk; kk++){
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	//-----------------------------------------------------------------
	const double Rmin = RR[ii] - 0.5 * hh;
	const double Rmax = RR[ii] + 0.5 * hh;
	const double rho0 = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] /
	  (DBL_MIN + gaussQuad2d4calcRho(Rmin, Rmax, NDISKBIN_HOR, RR, invhh, 0.0, hh, NDISKBIN_VER, zz, invhh, sph, invlogrbin_sph, Phi, lev, disk[kk]));
	//-----------------------------------------------------------------
	/* calculate vertical density profile */
	(*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] = (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)];
	for(int jj = 1; jj < NDISKBIN_VER; jj++)
	  (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] =
	    rho0 * gaussQuad2d4calcRho(Rmin, Rmax, NDISKBIN_HOR, RR, invhh, zz[jj] - 0.5 * hh, zz[jj] + 0.5 * hh, NDISKBIN_VER, zz, invhh, sph, invlogrbin_sph, Phi, lev, disk[kk]);
	//-----------------------------------------------------------------
#if 0
#pragma omp critical
	for(int jj = 0; jj < NDISKBIN_VER; jj++)
	  if( (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] < 0.0 ){
	    fprintf(stderr, "# lev = %d, ii = %d, jj = %d, rho = %e\n", lev, ii, jj, (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)]);

	    for(int mm = 0; mm < NDISKBIN_HOR; mm++){
	      for(int nn = 0; nn < NDISKBIN_VER; nn++)
		fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", RR[mm], zz[nn], (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, mm, nn)], disk[kk].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, mm, nn)], (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, mm, nn)]);
	      fprintf(stderr, "\n");
	    }

	    fflush(NULL);
	    exit(0);
	  }
#endif
	//-----------------------------------------------------------------
	/* calculate column density */
	double Sigma = 0.0;
	for(int jj = 0; jj < NDISKBIN_VER; jj++)
	  Sigma += (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)];
	Sigma *= 2.0 * hh;/* 2.0 reflects plane symmetry about the equatorial plane (z = 0) */
	//-----------------------------------------------------------------
	/* calibrate column density */
	const double Mscale = disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] / (DBL_MIN + Sigma);
	for(int jj = 0; jj < NDISKBIN_VER; jj++)
	  (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] *= Mscale;
	//-----------------------------------------------------------------
	const double errVal = (disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, ii)] > NEGLECT_DENSITY_MINIMUM * disk[kk].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, 0)]) ? fabs(Mscale - 1.0) : (0.0);
#pragma omp flush(errMax)
	if( errVal > errMax )
#pragma omp critical
	  if( errVal > errMax )
	    errMax = errVal;
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
      //-------------------------------------------------------------------
    }/* for(int kk = 0; kk < ndisk; kk++){ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    /* convergence tests */
    //---------------------------------------------------------------------
    /* convergence test for the column density */
    bool converge = true;
    if( errMax > CONVERGENCE_POTDENSPAIR )      converge = false;
#ifdef  PROGRESS_REPORT_ON
    fprintf(stdout, "# %d-th iteration: column density error is %e, procedure is %s\n", steps, errMax, (converge ? ("    converged") : ("not converged")));
    fflush(stdout);
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------
    /* convergence test for the density field */
    if( converge ){
      //-------------------------------------------------------------------
      static double worst;
      worst = DBL_EPSILON;
      //-------------------------------------------------------------------
      for(int kk = 0; kk < ndisk; kk++){
	//-----------------------------------------------------------------
#pragma omp parallel for
	for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	  //---------------------------------------------------------------
	  double errVal = 0.0;
	  //---------------------------------------------------------------
	  if( (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)] > NEGLECT_DENSITY_MINIMUM * (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, 0, 0)] )
	    for(int jj = 0; jj < NDISKBIN_VER; jj++){
	      //-----------------------------------------------------------
	      const double err =
		((*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] > NEGLECT_DENSITY_MINIMUM * (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, 0)]) ?
		(fabs(1.0 - (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)] / ((*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, lev, ii, jj)]))) : (0.0);
	      errVal = (errVal > err) ? errVal : err;
	      //-----------------------------------------------------------
	    }/* for(int jj = 0; jj < Nz; jj++){ */
	  //---------------------------------------------------------------
#pragma omp flush(worst)
	  if( errVal > worst )
#pragma omp critical
	    if( errVal > worst )
	      worst = errVal;
	}/* for(int ii = 0; ii < NR; ii++){ */
	//-----------------------------------------------------------------
      }/* for(int kk = 0; kk < ndisk; kk++){ */
      //-------------------------------------------------------------------
      if( worst > CONVERGENCE_POTDENSPAIR )
	converge = false;
      //-------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
      fprintf(stdout, "# %d-th iteration:        density error is %e, procedure is %s\n", steps, worst, (converge ? ("    converged") : ("not converged")));
      fflush(stdout);
#endif//PROGRESS_REPORT_ON
      //-------------------------------------------------------------------
    }/* if( converge ){ */
    //---------------------------------------------------------------------
    for(int kk = 0; kk < ndisk; kk++)
      swapDblArrays(disk[kk].rho, disk[kk].rhoSum);
    //---------------------------------------------------------------------
    /* final confirmation */
    if( converge )
      break;
    //---------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
    steps++;
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
  /* recover column density profile if necessary */
  if( lev > 0 )
    for(int ii = 0; ii < ndisk; ii++)
#pragma omp parallel for
      for(int jj = 0; jj < NDISKBIN_HOR; jj++)
	disk[ii].Sigma[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] = disk[ii].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)];
  //-----------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
  fprintf(stdout, "# procedure finish after %d iteration(s): NR = %d, Nz = %d\n#\n#\n", 1 + steps, NDISKBIN_HOR, NDISKBIN_VER);
  fflush(stdout);
#endif//PROGRESS_REPORT_ON
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void makeDiskPotentialTable(const int ndisk, const int maxLev, disk_data * restrict disk)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  /* initialize the table for Gaussian Quadrature provided by GSL */
  for(int ii = 0; ii < NTBL_GAUSS_QD_LOW; ii++){
    gsl_gaussQD_low_pos   [ii] = 0.0;
    gsl_gaussQD_low_weight[ii] = 0.0;
  }/* for(int ii = 0; ii < NTBL_GAUSS_QD_LOW; ii++){ */
  gsl_integration_glfixed_table *tab;
  tab = gsl_integration_glfixed_table_alloc(NINTBIN_LOW);
  int max = (NINTBIN_LOW >> 1) + (NINTBIN_LOW & 1);
  for(int ii = 0; ii < max; ii++){
    gsl_gaussQD_low_pos   [ii] = (*tab).x[(max - 1) - ii];
    gsl_gaussQD_low_weight[ii] = (*tab).w[(max - 1) - ii];
  }/* for(int ii = 0; ii < max; ii++){ */
  gsl_integration_glfixed_table_free(tab);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* memory for the iterative solver of the Poisson equation */
  //-----------------------------------------------------------------------
  /* sparse matrix in CRS format */
  crs mat, ilu;
  static double mat_val[NNZ_CG], ilu_val[NNZ_CG];
  static    int mat_col[NNZ_CG], ilu_col[NNZ_CG], mat_row[NROW_CG + 1], ilu_row[NROW_CG + 1];
  mat.val = mat_val;  mat.col = mat_col;  mat.row = mat_row;
  ilu.val = ilu_val;  ilu.col = ilu_col;  ilu.row = ilu_row;
  //-----------------------------------------------------------------------
  /* vector used in BiCGSTAB method */
  static double vec[NCOL_CG], res[NCOL_CG], sdw[NCOL_CG], mid[NCOL_CG], tmp[NCOL_CG];
  static double Api[NCOL_CG], Ati[NCOL_CG], Kri[NCOL_CG], Kpi[NCOL_CG], Kti[NCOL_CG];
  //-----------------------------------------------------------------------
  /* outer boundary condition of the potential field */
  static double Phi_NR[NDISKBIN_VER], Phi_Nz[NDISKBIN_HOR];
  //-----------------------------------------------------------------------
  /* execute first touch for memories related to BiCGSTAB */
#pragma omp parallel
  {
#pragma omp for nowait
    for(int ii = 0; ii < NCOL_CG; ii++){
      vec[ii] = res[ii] = sdw[ii] = mid[ii] = tmp[ii] = 0.0;
      Api[ii] = Ati[ii] = Kri[ii] = Kpi[ii] = Kti[ii] = 0.0;
    }/* for(int ii = 0; ii < NCOL_CG; ii++){ */
#pragma omp for nowait
    for(int ii = 0; ii < NNZ_CG; ii++){
      mat_val[ii] = 0.0;
      mat_col[ii] = 0;
    }/* for(int ii = 0; ii < NNZ_CG; ii++){ */
#pragma omp for nowait
    for(int ii = 0; ii < NROW_CG; ii++)
      mat_row[ii] = 0;
  }
  //-----------------------------------------------------------------------
  soaBiCGSTAB smat;
  smat.mat = mat;  smat.vec = vec;  smat.res = res;  smat.sdw = sdw;
  smat.mid = mid;  smat.tmp = tmp;  smat.Api = Api;  smat.Ati = Ati;
  soaPreConditioning pre;
  pre.ilu = ilu;  pre.Kri = Kri;  pre.Kpi = Kpi;  pre.Kti = Kti;
  //-----------------------------------------------------------------------
  /* prepare Poisson equation in matrix form (matrix part) */
  setSparseMatrix(smat.mat);
  getILU0(NDISKBIN_HOR * NDISKBIN_VER, smat.mat, pre.ilu);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* set column density profile */
  //-----------------------------------------------------------------------
  /* static inline void setColumnDensityProfile(const int ndisk, const int maxLev, disk_data * restrict disk, profile * restrict sph, const double invlogrbin_sph) */
  setColumnDensityProfile(ndisk, maxLev, disk, disk[0].prf, disk[0].invlogrbin);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* calculate GD formalism */
  //-----------------------------------------------------------------------
#ifdef  USE_GD_FORM_POTENTIAL
  //-----------------------------------------------------------------------
  double *stock_inv;  stock_inv = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  double *stock_sub;  stock_sub = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  double *stock_tmp;  stock_tmp = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  double *stock_sum;  stock_sum = malloc(sizeof(double) * ndisk * (size_t)CPUS_PER_PROCESS);
  if( stock_inv == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_inv\n");  }
  if( stock_sub == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_sub\n");  }
  if( stock_tmp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_tmp\n");  }
  if( stock_sum == NULL ){    __KILL__(stderr, "ERROR: failure to allocate stock_sum\n");  }
  //-----------------------------------------------------------------------
  const double Rmax = disk[0].Rmax;
  const double    abin = Rmax / (double)NRADBIN;
  const double invabin = 1.0 / abin;
  for(int kk = 0; kk < ndisk; kk++){
    //---------------------------------------------------------------------
#pragma omp parallel for
    for(int ii = 0; ii < NRADBIN; ii++){
      //-------------------------------------------------------------------
      disk[kk].prf[ii].rho = (double)ii * abin;/* store the position ``a'' */
#ifdef  NDIVIDE_GAUSSQD4DISK
      double sum = 0.0;
      double Rin = disk[kk].prf[ii].rho;
      const double a2 = Rin * Rin;
      const double Rbin = (Rmax - Rin) / (double)NDIVIDE_GAUSSQD4DISK;
      for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){
	sum += gaussQD4GDpot(Rin, Rin + Rbin, a2, disk[kk]);
	Rin += Rbin;
      }/* for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){ */
      disk[kk].prf[ii].enc = disk[kk].cfg->Sigma0 * sum;/* store the integral */
#else///NDIVIDE_GAUSSQD4DISK
      disk[kk].prf[ii].enc = disk[kk].cfg->Sigma0 * gaussQD4GDpot(disk[kk].prf[ii].rho, Rmax, disk[kk].prf[ii].rho * disk[kk].prf[ii].rho, disk[kk]);/* store the integral */
#endif//NDIVIDE_GAUSSQD4DISK
      //-------------------------------------------------------------------
    }/* for(int ii = 0; ii < NRADBIN; ii++){ */
    //---------------------------------------------------------------------
    disk[kk].prf[0].psi = (disk[kk].prf[1].enc - disk[kk].prf[0].enc) * invabin;
#pragma omp parallel for
    for(int ii = 1; ii < NRADBIN - 1; ii++)
      disk[kk].prf[ii].psi = (disk[kk].prf[ii + 1].enc - disk[kk].prf[ii - 1].enc) * 0.5 * invabin;/* store the derivative */
    disk[kk].prf[NRADBIN - 1].psi = (disk[kk].prf[NRADBIN - 1].enc - disk[kk].prf[NRADBIN - 2].enc) * invabin;
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < ndisk; kk++){ */
  //-----------------------------------------------------------------------
#endif//USE_GD_FORM_POTENTIAL
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* iterative procedure adopting Full Multigrid Algorithm (gamma = 2) */
  //-----------------------------------------------------------------------
  int lev = 0;
  int old = lev;
  int inc = 1;
  int top = 1;
  int num = 0;
  while( true ){
    //---------------------------------------------------------------------
#ifdef  PROGRESS_REPORT_ON
    fprintf(stdout, "# grid level: lev = %d, oldLev = %d\n", lev, old);
    fflush(stdout);
#endif//PROGRESS_REPORT_ON
    //---------------------------------------------------------------------
    getPotDensPair(ndisk, lev, old, disk, Phi_NR, Phi_Nz, smat, pre
#ifdef  USE_GD_FORM_POTENTIAL
		   , stock_inv, stock_sub, stock_tmp, stock_sum
#endif//USE_GD_FORM_POTENTIAL
		   );
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* procedure proposed by Press & Teukolsky (1991), Computers in Physics 5, 514, Multigrid Methods for Boundary Value Problems. I */
    //---------------------------------------------------------------------
    if( lev == top ){   num++;  inc = -1;
      if( num == 2 ){	top++;	num =  0;      }
    }/* if( lev == top ){ */
    //---------------------------------------------------------------------
    /* in case of the coarsest grid */
    if( (lev == 0) && (inc == -1) )      inc = 1;
    //---------------------------------------------------------------------
    /* exit the loop if we found the solution in the finest grid */
    if( lev == (maxLev - 1) )      break;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* go to the next level */
    old  = lev;
    lev += inc;
    //---------------------------------------------------------------------
  }/* while( true ){ */
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < ndisk; kk++)
    for(int ll = maxLev - 2; ll >= 0; ll--){
      //-------------------------------------------------------------------
      coarsePotential(ll,             disk[kk].pot );
      coarseDensity  (ll, disk[kk], *(disk[kk].rho));
      //-------------------------------------------------------------------
    }/* for(int ll = maxLev - 2; ll >= 0; ll--){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#ifdef  USE_GD_FORM_POTENTIAL
  //-----------------------------------------------------------------------
  free(stock_inv);
  free(stock_sub);
  free(stock_tmp);
  free(stock_sum);
  //-----------------------------------------------------------------------
#endif//USE_GD_FORM_POTENTIAL
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate 2 \pi \int_0^R dR' R' Sigma(R') */
  //-----------------------------------------------------------------------
#ifdef  NDIVIDE_GAUSSQD4DISK
  for(int kk = 0; kk < ndisk; kk++){
    //---------------------------------------------------------------------
    double Min = 0.0;
    double Rin = 0.0;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    for(int lev = maxLev - 1; lev >= 0; lev--){
      //-------------------------------------------------------------------
      const int nsub = (lev != (maxLev - 1)) ? ((NDISKBIN_HOR >> 1) / NDIVIDE_GAUSSQD4DISK) : (NDISKBIN_HOR / NDIVIDE_GAUSSQD4DISK);
      int head = (lev != (maxLev - 1)) ?  (NDISKBIN_HOR >> 1)                         : (0);
      //-------------------------------------------------------------------
      disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head)] = Min;
      //-------------------------------------------------------------------
      for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){
	//-----------------------------------------------------------------
	const int tail = head + nsub;
	//-----------------------------------------------------------------
#pragma omp parallel for
	for(int ii = head; ii < tail; ii++)
	  disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii + 1)] = Min + gaussQD4encSigma(Rin, disk[kk].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii + 1)], disk[kk]);
	//-----------------------------------------------------------------
	head = tail;
	Rin = disk[kk].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head)];
	Min = disk[kk].	    enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head)];
	//-----------------------------------------------------------------
      }/* for(int iter = 0; iter < NDIVIDE_GAUSSQD4DISK; iter++){ */
      //-------------------------------------------------------------------
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR + 1; ii++)
	disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii)] *= 2.0 * M_PI * disk[kk].cfg->Sigma0;
      //-------------------------------------------------------------------
    }/* for(int lev = maxLev - 1; lev >= 0; lev--){ */
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < ndisk; kk++){ */
#else///NDIVIDE_GAUSSQD4DISK
  for(int kk = 0; kk < ndisk; kk++){
    //---------------------------------------------------------------------
    double Min = 0.0;
    double Rin = 0.0;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    for(int lev = maxLev - 1; lev >= 0; lev--){
      //-------------------------------------------------------------------
      const int nsub = (lev != (maxLev - 1)) ? (NDISKBIN_HOR >> 1) : (NDISKBIN_HOR);
      const int head = (lev != (maxLev - 1)) ? (NDISKBIN_HOR >> 1) : (0);
      //-------------------------------------------------------------------
      disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head)] = Min;
      //-------------------------------------------------------------------
#pragma omp parallel for
      for(int ii = head; ii < head + nsub; ii++)
	disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii + 1)] = Min + gaussQD4encSigma(Rin, disk[kk].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii + 1)], disk[kk]);
      //-------------------------------------------------------------------
      Rin = disk[kk].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head + nsub)];
      Min = disk[kk].	  enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, head + nsub)];
      //-------------------------------------------------------------------
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR + 1; ii++)
	disk[kk].enc[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, ii)] *= 2.0 * M_PI * disk[kk].cfg->Sigma0;
      //-------------------------------------------------------------------
    }/* for(int lev = maxLev - 1; lev >= 0; lev--){ */
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < ndisk; kk++){ */
#endif//NDIVIDE_GAUSSQD4DISK
  //-----------------------------------------------------------------------
#if 0
  extern const double length2astro, mass2astro;
  for(int ii = 0; ii < NDISKBIN_HOR + 1; ii++)
    fprintf(stderr, "%e\t%e\n", disk[0].node_hor[ii] * length2astro, disk[0].enc[ii] * mass2astro);
  fflush(NULL);
  exit(0);
#endif
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* calculate int_0^z dz' \rho(R, z') */
  //-----------------------------------------------------------------------
  for(int kk = 0; kk < ndisk; kk++){
    for(int ll = 0; ll < maxLev; ll++){
      //-------------------------------------------------------------------
      const double hh = ldexp(disk[kk].hh, -ll);
      //-------------------------------------------------------------------
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	//-----------------------------------------------------------------
	const double dV = hh * hh * disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)];
	//-----------------------------------------------------------------
	double sum = 0.0;
	(*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii, 0)] = sum;
	//-----------------------------------------------------------------
	for(int jj = 0; jj < NDISKBIN_VER; jj++){
	  //---------------------------------------------------------------
	  sum += dV * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)];
	  (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii, jj + 1)] = sum;
	  //---------------------------------------------------------------
	}/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
	//-----------------------------------------------------------------
	/* for(int jj = 0; jj < NDISKBIN_VER; jj++){ */
	/*   //--------------------------------------------------------------- */
	/*   const double mass = 0.5 * dV * (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)]; */
	/*   //--------------------------------------------------------------- */
	/*   sum += mass; */
	/*   (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)] = sum; */
	/*   sum += mass; */
	/*   //--------------------------------------------------------------- */
	/* }/\* for(int jj = 0; jj < NDISKBIN_VER; jj++){ *\/ */
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
      //-------------------------------------------------------------------
    }/* for(int ll = 0; ll < maxLev; ll++){ */
    //---------------------------------------------------------------------
#if 0
    for(int ll = 0; ll < maxLev; ll++)
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	for(int jj = 0; jj < NDISKBIN_VER; jj++)
	  fprintf(stderr, "%e\t%e\t%e\t%e\t%e\n", disk[kk].hor[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)], disk[kk].ver[INDEX2D(maxLev, NDISKBIN_VER, ll, jj)],
		  (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)],
		  (*disk[kk].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)], disk[kk].pot[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, ll, ii, jj)]);
	fprintf(stderr, "\n");
      }
    fflush(NULL);
    exit(0);
#endif
    //---------------------------------------------------------------------
#if 0
    double Mtot0 = 0.0;
    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
      Mtot0 += (*disk[0].rhoSum)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, NDISKBIN_VER - 1)];
    Mtot0 *= 4.0 * M_PI;
    double Mtot1 = 0.0;
    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	Mtot1 += (*disk[0].rho)[INDEX2D(NDISKBIN_HOR, NDISKBIN_VER, ii, jj)] * disk[0].hor[ii] * disk[0].hh * disk[0].hh;
    Mtot1 *= 4.0 * M_PI;
    double Mtot2 = 0.0;
    Mtot2 = disk[0].Sigma[0] * disk[0].hor[0];
    for(int ii = 1; ii < NDISKBIN_HOR - 1; ii++)
      Mtot2 += (double)(1 << (1 + (ii & 1))) * disk[0].Sigma[ii] * disk[0].hor[ii];
    Mtot2 += disk[0].Sigma[NDISKBIN_HOR - 1] * disk[0].hor[NDISKBIN_HOR - 1];
    Mtot2 *= disk[0].hh / 3.0;
    Mtot2 *= 2.0 * M_PI;
    extern const double mass2astro;
    fprintf(stderr, "# Mtot0 = %e, Mtot1 = %e, Mtot2 = %e\n", Mtot0 * mass2astro, Mtot1 * mass2astro, Mtot2 * mass2astro);
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      double sum = 0.0;
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	sum += (*disk[0].rho)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER, 0, ii, jj)];
      sum *= 2.0 * disk[0].hh;
      fprintf(stderr, "%e\t%e\t%e\n", disk[0].hor[ii], disk[0].Sigma[ii], sum);
    }
    fflush(NULL);
    exit(0);
#endif
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* normalization */
    //---------------------------------------------------------------------
#if 0
#pragma omp parallel for
    for(int ii = 0; ii < NDISKBIN_HOR; ii++)
      disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)] = 1.0 / (DBL_MIN + (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, 0, ii, NDISKBIN_VER)]);
    for(int ll = 1; ll < maxLev; ll++)
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++){
	const double self = (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii    , NDISKBIN_VER)];
	const double pair = (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii ^ 1, NDISKBIN_VER)];
	disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)] = disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, ll - 1, ii >> 1)] * (self + pair) / (DBL_MIN + self);
      }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
#else
    for(int ii = 0; ii < NDISKBIN_HOR; ii++){
      disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, 0, ii)] = 1.0 / (DBL_MIN + (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, 0, ii, NDISKBIN_VER)]);
      for(int ll = 1; ll < maxLev; ll++){
	disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)] = disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, ll - 1, ii >> 1)] * ((*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii, NDISKBIN_VER)] + (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii ^ 1, NDISKBIN_VER)]) / (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii, NDISKBIN_VER)];
      }/* for(int ll = 1; ll < maxLev; ll++){ */
    }/* for(int ii = 0; ii < NDISKBIN_HOR; ii++){ */
#endif
    for(int ll = 0; ll < maxLev; ll++)
#pragma omp parallel for
      for(int ii = 0; ii < NDISKBIN_HOR; ii++)
	for(int jj = 0; jj < NDISKBIN_VER + 1; jj++)
	  (*disk[kk].rhoSum)[INDEX(maxLev, NDISKBIN_HOR, NDISKBIN_VER + 1, ll, ii, jj)] *= disk[kk].sigmaz[INDEX2D(maxLev, NDISKBIN_HOR, ll, ii)];
    //---------------------------------------------------------------------
  }/* for(int kk = 0; kk < ndisk; kk++){ */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
