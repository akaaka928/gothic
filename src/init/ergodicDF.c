/*************************************************************************\
 *                                                                       *
                  last updated on 2015/01/15(Thu) 16:47:11
 *                                                                       *
 *    Making Initial Condition Code of N-body Simulation                 *
 *       King sphere, Plummer sphere                                     *
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
//-------------------------------------------------------------------------
#include <gsl/gsl_rng.h>
//-------------------------------------------------------------------------
#include <macro.h>
#include <constants.h>
#include <timer.h>
#include <name.h>
//-------------------------------------------------------------------------
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "ergodicDF.h"
//-------------------------------------------------------------------------
extern const real newton;
static const real convergence = (real)1.0e-4;
static const real extreme     = (real)0.25;
//-------------------------------------------------------------------------
static const real inv3 = ONE_THIRD;
static const real inv6 = ONE_THIRD * HALF;
//-------------------------------------------------------------------------
gsl_rng *GSLRand;
#define UNIRAND ((real)gsl_rng_uniform(GSLRand))
#define RANDVAL (TWO * (UNIRAND) - UNITY)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void allocateArray(const int num, real **array)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  *array = (real *)malloc(sizeof(real) * num);
  if( *array == NULL ){
    __KILL__(stderr, "ERROR: failure to allocate array");
  }
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#define NADD (16384)
static inline void enlargeArray(const int num, real **array)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  *array = realloc(*array, sizeof(real) * num);
  if( *array == NULL ){
    __KILL__(stderr, "ERROR: failure to rescale array");
  }
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline real simpsonIntegral(int nbin, real vmin, real vmax, real W, real (*function)(real v2, real W))
{
  //-----------------------------------------------------------------------
  /* nbin must be an even number */
  if( nbin & 1 )    nbin++;
  real wbin = (vmax - vmin) / (real)nbin;
  //-----------------------------------------------------------------------
  real sum = (*function)(vmin * vmin, W);
  for(int ii = 1; ii < nbin; ii++){
    real vel = vmin + wbin * (real)ii;
    sum += (real)((1 + (ii & 1)) << 1) * (*function)(vel * vel, W);
  }
  sum += (*function)(vmax * vmax, W);
  sum *= wbin * inv3;
  //-----------------------------------------------------------------------
  return (sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* v2min = vz2, v2max = vesc2 */
static inline real singleSimpson(int nbin, real v2min, real v2max, real W, real (*function)(real v2, real W))
{
  //-----------------------------------------------------------------------
  /* nbin must be an even number */
  if( nbin & 1 )    nbin++;
  real wbin = (v2max - v2min) / (real)nbin;
  //-----------------------------------------------------------------------
  real sum = (*function)(v2min, W);
  for(int ii = 1; ii < nbin; ii++){
    real v2 = v2min + wbin * (real)ii;
    sum += (real)((1 + (ii & 1)) << 1) * (*function)(v2, W);
  }
  sum += (*function)(v2max, W);
  sum *= wbin * inv3;
  //-----------------------------------------------------------------------
  return (sum * v2min);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline real doubleSimpson(int nbin, real vmin, real vmax, real W, real (*function)(real v2, real W))
{
  //-----------------------------------------------------------------------
  /* nbin must be an even number */
  if( nbin & 1 )    nbin++;
  real wbin = (vmax - vmin) / (real)nbin;
  real v2max = vmax * vmax;
  //-----------------------------------------------------------------------
  real sum = singleSimpson(128, vmin * vmin, v2max, W, function);
  for(int ii = 1; ii < nbin; ii++){
    real vel = vmin + wbin * (real)ii;
    sum += (real)((1 + (ii & 1)) << 1) * singleSimpson(128, vel * vel, v2max, W, function);
  }
  sum *= wbin * inv3;
  //-----------------------------------------------------------------------
  return ((real)M_PI * sum);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline void isotropicDistribution(const real rad, real *vecx, real *vecy, real *vecz)
{
  //-----------------------------------------------------------------------
  const real proj = RANDVAL;
  *vecz = rad * proj;
  real Rproj  = rad * SQRT(UNITY - proj * proj);
  //-----------------------------------------------------------------------
  real theta = TWO * (real)M_PI * UNIRAND;
  *vecx = Rproj * COS(theta);
  *vecy = Rproj * SIN(theta);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* v2 := (v / vesc)^2 */
static inline real plummerDF(real v2)
{
  return (v2 * POW(UNITY - v2, (real)3.5));
}
//-------------------------------------------------------------------------
/* Disable ICC's remark #869: parameter "hoge" was never referenced */
#pragma warning (disable:869)
//-------------------------------------------------------------------------
static inline real v4fEplummer(real v2, real W){  return (v2 * v2 * plummerDF(v2));}
static inline real v2fEplummer(real v2, real W){  return (v2      * plummerDF(v2));}
/* static inline real v0fEplummer(real v2, real W){  return (          plummerDF(v2));} */
//-------------------------------------------------------------------------
/* Enable ICC's remark #869: parameter "hoge" was never referenced */
#pragma warning (enable:869)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void rescalePlummerSphere(real Mtot, real h, int *num, real **rad, real **rho)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  *num = 131072;
  allocateArray(*num, rad);
  allocateArray(*num, rho);
  //-----------------------------------------------------------------------
  real logrmin = LOG10(h / (real)1024);
  real logrmax = LOG10(h * (real)128);
  real logdr = (logrmax - logrmin) / (real)((*num) - 1);
  //-----------------------------------------------------------------------
  real h2inv = UNITY / (h * h);
  real Mtmp = (real)0.75 * Mtot * (real)M_1_PI * h2inv * h2inv * h;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < *num; ii++){
    (*rad)[ii] = POW((real)10.0, logrmin + logdr * (real)ii);
    real tmp = UNITY / (UNITY + (*rad)[ii] * (*rad)[ii] * h2inv);
    (*rho)[ii] = Mtmp * tmp * tmp * SQRT(tmp);
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void calcPlummerDynamicalProperties(real Mtot, real h, int num, real *rad, real **vdisp, real **phi)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  allocateArray(num, phi);
  allocateArray(num, vdisp);
  //-----------------------------------------------------------------------
  real h2inv = UNITY / (h * h);
  real sum = SQRTRATIO(simpsonIntegral( 16384, ZERO, UNITY, ZERO, v4fEplummer),
		       (simpsonIntegral(16384, ZERO, UNITY, ZERO, v2fEplummer) * THREE));
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++){
    real tmp = UNITY / (UNITY + rad[ii] * rad[ii] * h2inv);
    real vesc = SQRT(TWO * newton * Mtot * h2inv * h) * POW(tmp, QUARTER);
    (*vdisp)[ii] = vesc * sum;
    (*phi)[ii]   = -newton * Mtot * SQRT(tmp * h2inv);
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* line of sight is along z-axis */
/* assume mean line of velocity is zero (satisfied for the isotropic velocity models) */
void calcPlummerObservableProperties(int num, real sigma, real *rad, real *rho, real **Sigma, real **sigmalos)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  real *nu, *s2r;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  allocateArray(num, Sigma);
  allocateArray(num, sigmalos);
  allocateArray(num, &nu);
  allocateArray(num, &s2r);
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < num; ii++){
    (*Sigma)   [ii] = ZERO;
    (*sigmalos)[ii] = ZERO;
  }
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < num - 1; ii++){
    real mass = ZERO;
    real fvel = ZERO;
    real frac = ZERO;
    for(int jj = ii; jj < num - 1; jj++){
#if 1
      real tmp = (rad[jj + 1] - rad[jj]) * rad[jj + 1] / SQRT((rad[jj + 1] + rad[ii]) * (rad[jj + 1] - rad[ii]));
#else
      real tmp =  rad[jj + 1] - rad[jj];
#endif
      mass += rho[jj] * tmp;
      fvel += s2r[jj] * tmp;
      frac +=  nu[jj] * tmp;
    }
    (*Sigma)   [ii] = TWO * mass;
    (*sigmalos)[ii] = sigma * SQRTRATIO(fvel, frac);
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  free(nu);
  free(s2r);
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void outputFundamentalInformationOfPlummerSphere(const real Mtot, const real h, const int num, const ulong Ntot, const real eps, const real snapshotInterval, const real ft, real *rad, real *rho, real *vdisp, real *phi, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  FILE *fp;
  char filename[256], date[64];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output fundamental information of the Plummer model */
  sprintf(filename, "%s/%s.info.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "%s\n", "ERROR: fundamental information file of the Plummer model couldn't open.");
  }
  //-----------------------------------------------------------------------
  getPresentDateInStrings(date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Fundamental information about the Plummer sphere\n");
  fprintf(fp, "\tgenerated on %s", date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Total mass of the Plummer sphere Mtot is   %f\n", Mtot);
  fprintf(fp, "Scale radius of the Plummer sphere h is    %f\n", h);
  fprintf(fp, "Velocity dispersion (1D) at the center is  %f\n", vdisp[0]);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Number of arrays is %d\n", num);
  fprintf(fp, "#############################################################################\n");
  const real Ms = HALF * (real)M_SQRT1_2 * Mtot;
  const real Ns = (real)Ntot * (Ms / Mtot);
  const real tff = (real)M_PI_2 * h * SQRTRATIO(h, TWO * newton * Ms);
  const real t2r = tff * Ns / ((real)32.0 * LOG(h / eps));
  fprintf(fp, "Number of N-body particles to represent Plummer sphere is   %zu (= 2^%u)\n", Ntot, ilog2((uint)Ntot));
  fprintf(fp, "Length of Plummer softening is                              %e\n", eps);
  fprintf(fp, "Number of particles within the scale radius is              %e\n", Ns);
  fprintf(fp, "Enclosed mass within the scale radius is                    %e\n", Ms);
  fprintf(fp, "Free-fall time at the scale radius in computational unit is %e\n", tff);
  fprintf(fp, "Two-body relaxation time at the scale radius is             %e\n", t2r);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Snapshot interval in the computational unit               is %e\n", snapshotInterval);
  fprintf(fp, "Snapshot interval in the unit of free-fall time           is %e\n", snapshotInterval / tff);
  fprintf(fp, "Snapshot interval in the unit of two-body relaxation time is %e\n", snapshotInterval / t2r);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Final time of the simulation in the computational unit               is %e\n", ft);
  fprintf(fp, "Final time of the simulation in the unit of free-fall time           is %e\n", ft / tff);
  fprintf(fp, "Final time of the simulation in the unit of two-body relaxation time is %e\n", ft / t2r);
  fprintf(fp, "#############################################################################\n");
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output fundamental profile of the Plummer model */
  sprintf(filename, "%s/%s.profile.dat", DATAFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "%s\n", "ERROR: fundamental profile file of the Plummer model couldn't open.");
  }
  //-----------------------------------------------------------------------
  fprintf(fp, "#r\trho(r)\tvdisp(r)\tPhi(r)\n");
  fprintf(fp, "#\tgenerated on %s", date);
  for(int i = 0; i < num; i++){
    fprintf(fp, "%e\t%e\t%e\t%e\n", rad[i], rho[i], vdisp[i], phi[i]);
  }
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
static inline real kingDF(real W)
{
  return (EXP(W) - UNITY);
}
//-------------------------------------------------------------------------
/* v2 := (v / sigma)^2 */
/* W  := Psi / sigma^2 */
static inline real v4fEking(real v2, real W){  return (v2 * v2 * kingDF(W - HALF * v2));}
static inline real v2fEking(real v2, real W){  return (v2      * kingDF(W - HALF * v2));}
static inline real v0fEking(real v2, real W){  return (          kingDF(W - HALF * v2));}
//-------------------------------------------------------------------------
static inline real calcKingDensity(real W, real rho1)
{
  //-----------------------------------------------------------------------
  return (rho1 * (EXP(W) * ERF(SQRT(W)) - SQRT((real)4.0 * W * (real)M_1_PI) * (UNITY + TWO * W * inv3)));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* static inline real kingFunc1(real rad, real uu, real dWdr) */
#define kingFunc1(rad, uu, dWdr) __kingFunc1(rad, dWdr)
static inline real __kingFunc1(real rad, real dWdr)
{
  //-----------------------------------------------------------------------
  real func;
  //-----------------------------------------------------------------------
  if( rad > EPSILON ){
    func = dWdr / (rad * rad);
  }
  else{
    func = ZERO;
  }
  //-----------------------------------------------------------------------
  return (func);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* static inline real kingFunc2(real rad, real W, real dWdr, real rho0inv, real rho1) */
#define kingFunc2(rad, W, dWdr, rho0inv, rho1) __kingFunc2(rad, W, rho0inv, rho1)
static inline real __kingFunc2(real rad, real W, real rho0inv, real rho1)
{
  //-----------------------------------------------------------------------
  return  (-(real)9.0 * rad * rad * rho0inv * calcKingDensity(W, rho1));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static void rungeKutta4thForKing(real r, real *u, real *y, real h, real rho0inv, real rho1)
{
  //-----------------------------------------------------------------------
  real k1u = kingFunc1(r, *u, *y);
  real k1y = kingFunc2(r, *u, *y, rho0inv, rho1);
  real k2u = kingFunc1(r + HALF * h, *u + HALF * h * k1u, *y + HALF * h * k1y);
  real k2y = kingFunc2(r + HALF * h, *u + HALF * h * k1u, *y + HALF * h * k1y, rho0inv, rho1);
  real k3u = kingFunc1(r + HALF * h, *u + HALF * h * k2u, *y + HALF * h * k2y);
  real k3y = kingFunc2(r + HALF * h, *u + HALF * h * k2u, *y + HALF * h * k2y, rho0inv, rho1);
  real k4u = kingFunc1(r +        h, *u +        h * k3u, *y +        h * k3y);
  real k4y = kingFunc2(r +        h, *u +        h * k3u, *y +        h * k3y, rho0inv, rho1);
  //-----------------------------------------------------------------------
  *u += (h * inv6) * (k1u + TWO * k2u + TWO * k3u + k4u);
  *y += (h * inv6) * (k1y + TWO * k2y + TWO * k3y + k4y);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* W is the dimensionless King parameter Psi / sigma^2 */
static void solvePoissonEqOfKingDF(real W, real **Rad, real **Psi, real **Rho, int *num, int *rem)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  real *rad, *psi, *rho;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* set initial condition */
  //-----------------------------------------------------------------------
  real h = UNITY / (real)1024.0;
  real hnew = h;
  //-----------------------------------------------------------------------
  /* assumption is rcal := r / r0, r0 = 1.0, sigma = 1.0 */
  //-----------------------------------------------------------------------
  /* see Eq (4.106) in Galactic Dynamics, second edtion, page 305 */
  /* rho0 = (real)2.25 * sigma * sigma / (pi * newton * r0 * r0); */
  real rho0    = (real)(2.25 * M_1_PI);
  real rho0inv = UNITY / rho0;
  //-----------------------------------------------------------------------
  /* u stores the dimensionless King parameter W, y stores dW/dr */
  real u = W;
  /* boundary condition is dW/dr(r = 0) = 0: assume core profile at the center */
  real y = ZERO;
  //-----------------------------------------------------------------------
  int ii = 0;
  rad = *Rad;  rad[ii] = ZERO;
  psi = *Psi;  psi[ii] = u;
  rho = *Rho;  rho[ii] = calcKingDensity(psi[ii], UNITY);
  //-----------------------------------------------------------------------
  real rho1 = rho0 / rho[ii];
  rho[ii] = calcKingDensity(psi[ii], rho1);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* solve Poisson Equation using 4th-order Runge-Kutta method */
  //-----------------------------------------------------------------------
  while( true ){
    //---------------------------------------------------------------------
    h = hnew;
    real uold = u;
    real yold = y;
    real uoldinv = UNITY / uold;
    while( true ){
      u = uold;
      y = yold;
      rungeKutta4thForKing(rad[ii], &u, &y, h, rho0inv, rho1);
      real diff = FABS(u - uold) * uoldinv;
      if( diff < convergence ){
	if( diff < extreme * convergence )
	  hnew = h * TWO;
	break;
      }
      else{
	h *= HALF;
	hnew = h;
      }
    }
    //---------------------------------------------------------------------
    if( *rem == 0 ){
      enlargeArray((*num) + NADD, Rad);      rad = *Rad;
      enlargeArray((*num) + NADD, Psi);      psi = *Psi;
      enlargeArray((*num) + NADD, Rho);      rho = *Rho;
      *rem += NADD;
    }
    ii++;    *num += 1;    *rem -= 1;
    rad[ii] = rad[ii - 1] + h;
    real tmp = calcKingDensity(u, rho1);
    //---------------------------------------------------------------------
    if( (u > EPSILON) && (tmp > EPSILON) ){
      psi[ii] = u;
      rho[ii] = tmp;
    }
    //---------------------------------------------------------------------
    else{
      psi[ii] = ZERO;
      rho[ii] = ZERO;
      break;
    }
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void makeKingDFTable(real W0, int *num, real **rad, real **W, real **rho)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  *num = 0;
  int rem  = 1024;
  //-----------------------------------------------------------------------
  allocateArray(rem, rad);
  allocateArray(rem, W);
  allocateArray(rem, rho);
  //-----------------------------------------------------------------------
  solvePoissonEqOfKingDF(W0, rad, W, rho, num, &rem);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void rescaleKingSphere(real Mtot, real rt, real *r0, real *sigma, int num, real *rad, real *rho, real **Menc)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  allocateArray(num, Menc);
  //-----------------------------------------------------------------------
  real MencCalc = TWO * (real)M_PI;/* 4pi / 2 */
  (*Menc)[0] = MencCalc * rad[1] * rad[1] * rad[1] * (rho[0] + rho[1]) * inv3;
  for(int ii = 1; ii < num - 1; ii++){
    (*Menc)[ii] = MencCalc * rad[ii + 1] * rad[ii + 1] * (rad[ii + 1] - rad[ii]) * (rho[ii] + rho[ii + 1]) + (*Menc)[ii - 1];
  }
  MencCalc = (*Menc)[num - 1] = (*Menc)[num - 2];
  real rtCalc = rad[num - 1];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* evaluate scaling factors */
  //-----------------------------------------------------------------------
  real lengthUnit = rt   /   rtCalc;
  real   massUnit = Mtot / MencCalc;
  //-----------------------------------------------------------------------
  real      rhoUnit = massUnit / (lengthUnit * lengthUnit * lengthUnit);
  real     timeUnit = RSQRT(newton * rhoUnit);
  real velocityUnit = lengthUnit / timeUnit;
  //-----------------------------------------------------------------------
  *r0 = lengthUnit;
  *sigma = velocityUnit;
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++){
    rad    [ii] *= lengthUnit;
    rho    [ii] *=    rhoUnit;
    (*Menc)[ii] *=   massUnit;
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void calcKingDynamicalProperties(real Mtot, real rt, real sigma, int num, real *W, real **phi, real **psi, real **vdisp)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  allocateArray(num, phi);
  allocateArray(num, psi);
  allocateArray(num, vdisp);
  //-----------------------------------------------------------------------
  real tmp = newton * Mtot / rt;
#pragma omp parallel for
  for(int ii = 0; ii < num; ii++){
    (*psi)  [ii] = W[ii] * sigma * sigma;
    (*phi)  [ii] = -(*psi)[ii] - tmp;
    (*vdisp)[ii] = sigma * SQRTRATIO(simpsonIntegral (1024, ZERO, SQRT(TWO * W[ii]), W[ii], v4fEking),
				     (simpsonIntegral(1024, ZERO, SQRT(TWO * W[ii]), W[ii], v2fEking) * THREE));
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
/* line of sight is along z-axis */
/* assume mean line of velocity is zero (satisfied for the isotropic velocity models) */
void calcKingObservableProperties(int num, real sigma, real *rad, real *rho, real *W, real **Sigma, real **sigmalos)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  real *nu, *s2r;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  allocateArray(num, Sigma);
  allocateArray(num, sigmalos);
  allocateArray(num, &nu);
  allocateArray(num, &s2r);
  //-----------------------------------------------------------------------
#pragma omp parallel for
  for(int ii = 0; ii < num; ii++){
    nu         [ii] = simpsonIntegral(1024, ZERO, SQRT(TWO * W[ii]), W[ii], v2fEking) * (real)(4.0 * M_PI);
    (*Sigma)   [ii] = ZERO;
    (*sigmalos)[ii] = ZERO;
    s2r        [ii] = TWO * doubleSimpson(128, ZERO, SQRT(TWO * W[ii]), W[ii], v0fEking);
  }
  //-----------------------------------------------------------------------
#pragma omp parallel
  for(int ii = 0; ii < num - 1; ii++){
    real mass = ZERO;
    real fvel = ZERO;
    real frac = ZERO;
    for(int jj = ii; jj < num - 1; jj++){
#if 1
      real tmp = (rad[jj + 1] - rad[jj]) * rad[jj + 1] / SQRT((rad[jj + 1] + rad[ii]) * (rad[jj + 1] - rad[ii]));
#else
      real tmp =  rad[jj + 1] - rad[jj];
#endif
      mass += rho[jj] * tmp;
      fvel += s2r[jj] * tmp;
      frac +=  nu[jj] * tmp;
    }
    (*Sigma)   [ii] = TWO * mass;
    (*sigmalos)[ii] = sigma * SQRTRATIO(fvel, frac);
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  free(nu);
  free(s2r);
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void outputFundamentalInformationOfKingSphere(const real W0, const real Mtot, const real r0, const real rt, const real sigma, const int num, const ulong Ntot, const real eps, const real snapshotInterval, const real ft, real *rad, real *rho, real *vdisp, real *Sigma, real *vdlos, real *Menc, real *phi, real *psi, real *W, char file[])
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  FILE *fp;
  char filename[256], date[64];
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output fundamental information of the King model */
  sprintf(filename, "%s/%s.info.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "%s\n", "ERROR: fundamental information file of the King model couldn't open.");
  }
  //-----------------------------------------------------------------------
  getPresentDateInStrings(date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Fundamental information about the King sphere\n");
  fprintf(fp, "\tgenerated on %s", date);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Dimensionless King parameter W0 is                  %f\n", W0);
  fprintf(fp, "Total mass of the King sphere Mtot is               %f\n", Mtot);
  fprintf(fp, "Core radius of the King sphere r0 is                %f\n", r0);
  fprintf(fp, "Tidal radius of the King sphere rt is               %f\n", rt);
  fprintf(fp, "Concentration parameter c is                        %f\n", LOG10(rt / r0));
  fprintf(fp, "Value of the parameter sigma is                     %f\n", sigma);
  fprintf(fp, "Velocity dispersion (1D) at the center is           %f\n", vdisp[0]);
  fprintf(fp, "Line-of-sight velocity dispersion at the center is  %f\n", vdlos[0]);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Required number of arrays is %d to satisfy convergence of %e\n", num, convergence);
  fprintf(fp, "#############################################################################\n");
  real Ms;
  {
    int ii = 0;
    while( true ){
      if( rad[ii] > r0 ){
	ii--;
	break;
      }
      ii++;
      if( ii == (num - 1) )
	break;
    }
    Ms = Menc[ii];
  }
  const real Ns = (real)Ntot * (Ms / Mtot);
  const real tff = (real)M_PI_2 * r0 * SQRTRATIO(r0, TWO * newton * Ms);
  const real t2r = tff * Ns / ((real)32.0 * LOG(r0 / eps));
  fprintf(fp, "Number of N-body particles to represent King sphere is     %zu (= 2^%u)\n", Ntot, ilog2((uint)Ntot));
  fprintf(fp, "Length of Plummer softening is                             %e\n", eps);
  fprintf(fp, "Number of particles within the King radius is              %e\n", Ns);
  fprintf(fp, "Enclosed mass within the King radius is                    %e\n", Ms);
  fprintf(fp, "Free-fall time at the King radius in computational unit is %e\n", tff);
  fprintf(fp, "Two-body relaxation time at the King radius is             %e\n", t2r);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Snapshot interval in the computational unit               is %e\n", snapshotInterval);
  fprintf(fp, "Snapshot interval in the unit of free-fall time           is %e\n", snapshotInterval / tff);
  fprintf(fp, "Snapshot interval in the unit of two-body relaxation time is %e\n", snapshotInterval / t2r);
  fprintf(fp, "#############################################################################\n");
  fprintf(fp, "Final time of the simulation in the computational unit               is %e\n", ft);
  fprintf(fp, "Final time of the simulation in the unit of free-fall time           is %e\n", ft / tff);
  fprintf(fp, "Final time of the simulation in the unit of two-body relaxation time is %e\n", ft / t2r);
  fprintf(fp, "#############################################################################\n");
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* output fundamental profile of the King model */
  sprintf(filename, "%s/%s.profile.dat", DATAFOLDER, file);
  fp = fopen(filename, "w");
  if( fp == NULL ){
    __KILL__(stderr, "%s\n", "ERROR: fundamental profile file of the King model couldn't open.");
  }
  //-----------------------------------------------------------------------
  fprintf(fp, "#r\trho(r)\tvdisp(r)\tSigma(R)\tvdlos(R)\tM(r)\tPhi(r)\tPsi(r)\tW(r)\n");
  fprintf(fp, "#\tgenerated on %s", date);
  for(int ii = 0; ii < num; ii++){
    fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", rad[ii], rho[ii], vdisp[ii], Sigma[ii], vdlos[ii], Menc[ii], phi[ii], psi[ii], W[ii]);
  }
  fclose(fp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void makeKingSphereByNbody(ulong num, nbody_particle *body, real Mtot, real sigma, int Nking, real *rad, real *W, real *Menc)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  real velocity;
  //-----------------------------------------------------------------------
  const gsl_rng_type *RandType;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  gsl_rng_env_setup();
  RandType = gsl_rng_mt19937;
  GSLRand  = gsl_rng_alloc(RandType);
  gsl_rng_set(GSLRand, 5489);
  //-----------------------------------------------------------------------
  real mass = Mtot / (real)num;
  for(ulong ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    /* determine spatial distribution */
    //---------------------------------------------------------------------
    real tmp = UNIRAND * Mtot;
    uint left  = 0;
    uint right = (Nking - 1);
    while( true ){
      //-------------------------------------------------------------------
      uint center = (left + right) >> 1;
      if( (Menc[center] - tmp) * (Menc[left] - tmp) < ZERO )  	right = center;
      else  	                                                 left = center;
      //-------------------------------------------------------------------
      if( (left + 1) == right )  	break;
      //-------------------------------------------------------------------
    }
    //---------------------------------------------------------------------
    real radius = ((tmp - Menc[left]) / (Menc[right] - Menc[left])) * (rad[right] - rad[left]) + rad[left];
    isotropicDistribution(radius, &(body[ii].x), &(body[ii].y), &(body[ii].z));
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    /* determine velocity distribution */
    //---------------------------------------------------------------------
    tmp = ((W[right] - W[left]) / (rad[right] - rad[left])) * (radius - rad[left]) + W[left];
    radius = SQRT(TWO * tmp);/* escape velocity */
    while( true ){
      velocity = radius * UNIRAND;
      real val = velocity * velocity * (EXP(-HALF * velocity * velocity) - EXP(-tmp));
      real try = radius * radius * UNIRAND;
      if( val > try ){
    	break;
      }
    }
    //---------------------------------------------------------------------
    isotropicDistribution(sigma * velocity, &(body[ii].vx), &(body[ii].vy), &(body[ii].vz));
    //---------------------------------------------------------------------
    body[ii].ax = body[ii].ay = body[ii].az = ZERO;
    body[ii].m   = mass;
    body[ii].pot = ZERO;
    //---------------------------------------------------------------------
    body[ii].idx = ii;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
