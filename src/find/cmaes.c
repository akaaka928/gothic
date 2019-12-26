/**
 * @file cmaes.c
 *
 * @brief Source code for CMA-ES (Evolution Strategy with Covariance Matrix Adaptation)
 *
 * @author Yohei MIKI (University of Tokyo)
 *
 * @date 2019/12/26 (Thu)
 *
 * Copyright (C) 2019 Yohei MIKI
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#ifdef  SUPPORT_CMAES_RESTART
#ifndef DISABLE_MPI
#include <mpi.h>
#include "mpilib.h"
#endif//DISABLE_MPI

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT
#endif//SUPPORT_CMAES_RESTART

#include <lapacke.h>

#include "macro.h"
#include "name.h"
#include "myutil.h"
#include "timer.h"
#include "rand.h"

#include "../find/cmaes.h"
#include "../find/cmaes_io.h"


static inline int scoreAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((const struct dbl_int *)a)->val > ((const struct dbl_int *)b)->val ){    return ( 1);  }
  else{    if( ((const struct dbl_int *)a)->val < ((const struct dbl_int *)b)->val ){    return (-1);  }
    else{                                                                                return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}


static inline bool __chkLAPACKerr(lapack_int err, const char *file, const int line)
{
  if( err != 0 ){
    fprintf(stderr, "Error is detected at %s(%d)\n", file, line);
    if( err < 0 ){
      /* __KILL__(stderr, "LAPACK error, %d-th value had an illegal value\n", -err); */
      fprintf(stderr, "LAPACK error, %d-th value had an illegal value\n", -err);
    }/* if( err < 0 ){ */
    else{
      /* __KILL__(stderr, "LAPACK error, returned value is %d\n", err); */
      fprintf(stderr, "LAPACK error, returned value is %d\n", err);
    }/* else{ */
    return (false);
  }
  return (true);
}


/**
 * @fn update_eigen
 *
 * @brief update eigensystem
 * @detail see cmaes_UpdateEigensystem: line 1885--1947
 *
 */
static inline bool update_eigen(struct cmaes_status *cma, struct cmaes_config cfg, const bool force)
{
  __NOTE__("%s\n", "start");

  if( (!force) && cma->eigen_uptodate )
    return (true);

  const int ndim = cfg.ndim;
  for(int ii = 0; ii < ndim; ii++)
    for(int jj = ii; jj < ndim; jj++){
      const int idx = INDEX2D(ndim, ndim, ii, jj);
      cma->orth[idx] = cma->mat[idx];
    }/* for(int jj = ii; jj < ndim; jj++){ */

  bool success = true;
  success &= chkLAPACKerr(LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'u', ndim, cma->orth, ndim, cma->diag, cma->offd, cma->tau));
  success &= chkLAPACKerr(LAPACKE_dorgtr(LAPACK_ROW_MAJOR, 'u', ndim, cma->orth, ndim,                       cma->tau));
  success &= chkLAPACKerr(LAPACKE_dpteqr(LAPACK_ROW_MAJOR, 'v', ndim, cma->diag, cma->offd, cma->orth, ndim));

  /** orth stores matrix with normalized eigenvectors in columns */
  /** diag stores square root of eigenvalues */
  for(int ii = 0; ii < ndim; ii++)
    cma->diag[ii] = sqrt(cma->diag[ii]);

  cma->eigen_uptodate = true;

  __NOTE__("%s\n", "end");

  return (success);
}


static inline double random_gauss(struct cmaes_rand *gauss, rand_state *rand)
{
  if( gauss->stored ){
    gauss->stored = false;
    return (gauss->hold);
  }/* if( gauss->stored ){ */

  double xx, yy, r2 = 0.0;
  while( (r2 < DBL_EPSILON) || ((r2 + EPSILON) > 1.0) ){
    xx = RANDVAL_DBL(rand);
    yy = RANDVAL_DBL(rand);

    r2 = xx * xx + yy * yy;
  }/* while( (r2 < DBL_EPSILON) || ((r2 + EPSILON) > 1.0) ){ */

  const double fac = sqrt(-2.0 * log(r2) / r2);
  gauss->stored = true;
  gauss->hold = fac * xx;
  return (fac * yy);
}


/**
 * @fn sampling
 *
 * @brief generate sampling points
 * @detail see cmaes_SamplePopulation in dm:~/work/181025cmaes/c-cmaes/src/cmaes.c: line 612--659
 *
 * @return (cmaes) status of CMA-ES
 */
bool sampling
(struct cmaes_status *cma, struct cmaes_config cfg, struct cmaes_rand *gauss, rand_state *rand
#ifdef  SUPPORT_CMAES_RESTART
 , char *file
#ifdef  USE_HDF5_FORMAT
 , const hid_t H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
 , const hid_t H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
#endif//SUPPORT_CMAES_RESTART
 )
{
  __NOTE__("%s\n", "start");

  const int state = cma->state;
  if( (state == CMAES_STATE_UPDATED) || (state == CMAES_STATE_INITIAL) ){
    const int ndim = cfg.ndim;
    /** calculate eigensystem if required */
    if( !cma->eigen_uptodate )
      if( !update_eigen(cma, cfg, false) )
	return (false);

    const int lambda = cfg.lambda;
    for(int ii = 0; ii < lambda; ii++){
      /** generate scaled random vectors */
      for(int jj = 0; jj < ndim; jj++)
	cma->xtemp[jj] = cma->diag[jj] * random_gauss(gauss, rand);

      /** add mutation */
      for(int jj = 0; jj < ndim; jj++){
	double sum = 0.0;
	for(int kk = 0; kk < ndim; kk++)
	  sum += cma->orth[INDEX2D(ndim, ndim, jj, kk)] * cma->xtemp[kk];
	cma->xx[INDEX2D(lambda, ndim, ii, jj)] = cma->xmean[jj] + cma->sigma * sum;/**< Eqs.(4) and (5) in Hansen (2016) */
      }/* for(int jj = 0; jj < ndim; jj++){ */
    }/* for(int ii = 0; ii < lambda; ii++){ */

    cma->gen++;
    cma->state = CMAES_STATE_SAMPLED;

#ifdef  SUPPORT_CMAES_RESTART
    /** dump all status for evaluating functions by using different programs */
    saveStatus(file, cfg, *cma, *gauss, *rand
#ifdef  USE_HDF5_FORMAT
	       , H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
	       , H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
	       );
#endif//SUPPORT_CMAES_RESTART
  }/* if( (state == CMAES_STATE_UPDATED) || (state == CMAES_STATE_INITIAL) ){ */

  __NOTE__("%s\n", "end");

  return (true);
}


/**
 * @fn resampling
 *
 * @brief retry sampling for a specific offspring
 * @detail see cmaes_ReSampleSingle in dm:~/work/181025cmaes/c-cmaes/src/cmaes.c: line 685--708
 *
 * @return (cmaes) status of CMA-ES
 */
void resampling(struct cmaes_status *cma, struct cmaes_config cfg, struct cmaes_rand *gauss, rand_state *rand, const int ii)
{
  const int ndim = cfg.ndim;

  /** generate scaled random vectors */
  for(int jj = 0; jj < ndim; jj++)
    cma->xtemp[jj] = cma->diag[jj] * random_gauss(gauss, rand);

  /** add mutation */
  for(int jj = 0; jj < ndim; jj++){
    double sum = 0.0;
    for(int kk = 0; kk < ndim; kk++)
      sum += cma->orth[INDEX2D(ndim, ndim, jj, kk)] * cma->xtemp[kk];
    cma->xx[INDEX2D(lambda, ndim, ii, jj)] = cma->xmean[jj] + cma->sigma * sum;/**< Eqs.(4) and (5) in Hansen (2016) */
  }/* for(int jj = 0; jj < ndim; jj++){ */
}


/**
 * @fn adapt
 *
 * @brief update covariance matrix
 * @detail see Adapt_C2 in dm:~/work/181025cmaes/c-cmaes/src/cmaes.c: line 963--1001
 */
static inline void adapt(struct cmaes_status *cma, struct cmaes_config cfg)
{
  __NOTE__("%s\n", "start");

  const int ndim = cfg.ndim;
  const int lambda = cfg.lambda;

  if( !cma->initial_phase ){
    const double c1 = cfg.c1;
    const double cmu = cfg.cmu;
    const double coeff = cfg.coeff;
    const double invsig2 = 1.0 / (DBL_MIN + cma->sigma * cma->sigma);

    cma->eigen_uptodate = false;

    /** update covariance matrix */
    for(int ii = 0; ii < ndim; ii++){
      for(int jj = ii; jj < ndim; jj++){
	double cov = coeff * cma->mat[INDEX2D(ndim, ndim, ii, jj)]
	  + c1 * cma->pc[ii] * cma->pc[jj];/**< rank-one update */
	/** additional rank mu update */
	for(int kk = 0; kk < lambda; kk++){
	  const int ll = cma->sorted[kk].idx;
	  cov += cmu * cfg.weight[kk] * (cma->xx[INDEX2D(lambda, ndim, ll, ii)] - cma->xmold[ii])  * (cma->xx[INDEX2D(lambda, ndim, ll, jj)] - cma->xmold[jj]) * invsig2;
	}/* for(int kk = 0; kk < mu; kk++){ */
	cma->mat[INDEX2D(ndim, ndim, ii, jj)] = cov;/**< Eq.(30) in Hansen (2016) */
      }/* for(int jj = ii; jj < ndim; jj++){ */
    }/* for(int ii = 0; ii < ndim; ii++){ */
  }/* if( !cma->initial_phase ){ */

  __NOTE__("%s\n", "end");
}


/**
 * @fn update
 *
 * @brief update distribution
 * @detail see cmaes_UpdateDistribution in dm:~/work/181025cmaes/c-cmaes/src/cmaes.c: line 797--958
 *
 * @return (cmaes) status of CMA-ES
 */
void update(struct cmaes_status *cma, struct cmaes_config cfg)
{
  __NOTE__("%s\n", "start");

  const int gen = cma->gen;
  const int ndim = cfg.ndim;
  const int lambda = cfg.lambda;
  if( cma->state != CMAES_STATE_SCORED ){
    __KILL__(stderr, "ERROR (generation %d): unexpected status %d while %d is assumed.\n", gen, cma->state, CMAES_STATE_SAMPLED);
  }/* if( cma->state == CMAES_STATE_SCORED ){ */

  cma->count += lambda;

  /** assign function values */
  for(int ii = 0; ii < lambda; ii++){
    cma->sorted[ii].val = cma->score[ii];
    cma->sorted[ii].idx = ii;
  }/* for(int ii = 0; ii < lambda; ii++){ */

  /** obtain index in ascending order of score (score[index[0]] <= score[index[1]] <= score[index[2]] <= ...) */
  qsort(cma->sorted, lambda, sizeof(struct dbl_int), scoreAscendingOrder);

  /** test if function values are identical, escape flat fitness */
  if( cma->sorted[0].val == cma->sorted[lambda / 2].val ){
    cma->sigma *= exp(0.2 + cfg.cs / cfg.ds);
    __FPRINTF__(stderr, "Warning (generation %d): sigma increased due to equal function values.\n", gen);
  }/* if( cma->sorted[0].val == cma->sorted[lambda >> 1].val ){ */


  /** update function value history */
  const int numHist = (gen < cfg.maxHistory) ? gen : cfg.maxHistory;
  for(int ii = numHist - 1; ii > 0; ii--)
    cma->history[ii] = cma->history[ii - 1];
  cma->history[0] = cma->sorted[0].val;

  /** update xbest */
  if( cma->bestScore > cma->sorted[0].val ){
    cma->bestScore = cma->sorted[0].val;
    const int ll = cma->sorted[0].idx;
    for(int ii = 0; ii < ndim; ii++)
      cma->xbest[ii] = cma->xx[INDEX2D(lambda, ndim, ll, ii)];
  }/* if( cma->bestScore > cma->sorted[0].val ){ */

  /** update xmean */
  const int mu = cfg.mu;
  const double tmp = sqrt(cfg.mueff) / cma->sigma;
  for(int ii = 0; ii < ndim; ii++){
    const double xmold = cma->xmold[ii] = cma->xmean[ii];
#if 1
    double xmean = xmold;
    for(int jj = 0; jj < mu; jj++)
      xmean += cfg.cm * cfg.weight[jj] * (cma->xx[INDEX2D(lambda, ndim, cma->sorted[jj].idx, ii)] - xmold);
    cma->xmean[ii] = xmean;/**< Eq.(9) in Hansen (2016) */
#else
    double xmean = 0.0;
    for(int jj = 0; jj < mu; jj++)
      xmean += cfg.weight[jj] * cma->xx[INDEX2D(lambda, ndim, cma->sorted[jj].idx, ii)];
    cma->xmean[ii] = xmean;/**< Eq.(6) in Hansen (2016) */
#endif
    cma->xevol[ii] = tmp * (xmean - xmold);/**< Eq.(24) in Hansen (2016) */
  }/* for(int ii = 0; ii < ndim; ii++){ */

  /**< calculate D^-1 * B^T * xevol into xtemp */
  for(int ii = 0; ii < ndim; ii++){
    double sum = 0.0;
    for(int jj = 0; jj < ndim; jj++)
      sum += cma->orth[INDEX2D(ndim, ndim, jj, ii)] * cma->xevol[jj];/**< B^T * xevol */
    cma->xtemp[ii] = sum / cma->diag[ii];/**< D^-1 * B^T * xevol */
  }/* for(int ii = 0; ii < ndim; ii++){ */

  /** cumulation for sigma (ps) */
  const double cs = cfg.cs;
  const double coeff = sqrt(cs * (2.0 - cs));
  double ps2 = 0.0;
  for(int ii = 0; ii < ndim; ii++){
    double sum = 0.0;
    for(int jj = 0; jj < ndim; jj++)
      sum += cma->orth[INDEX2D(ndim, ndim, ii, jj)] * cma->xtemp[jj];/**< B * D^-1 * B^T * xevol */
    const double val = (1.0 - cs) * cma->ps[ii] + coeff * sum;/**< Eq.(31) in Hansen (2016) */
    ps2 += val * val;
    cma->ps[ii] = val;
  }/* for(int ii = 0; ii < ndim; ii++){ */

  /* cumulation for covariance matrix (pc) */
  const double cc = cfg.cc;
  const double fac = sqrt(cc * (2.0 - cc));
  for(int ii = 0; ii < ndim; ii++)
    cma->pc[ii] = (1.0 - cc) * cma->pc[ii] + fac * cma->xevol[ii];/**< Eq.(24) in Hansen (2016) */

  /* stop initial phase */
  if( cma->initial_phase && ((double)gen > fmin(1.0 / cs, 1.0 + (double)ndim / cfg.mueff)) )
    if( ps2 < ((cfg.ds * (1.0 - pow(1.0 - cs, gen))) * (double)ndim * 1.05) )
      cma->initial_phase = false;

  adapt(cma, cfg);

  /** update sigma */
  cma->sigma *= exp(((sqrt(ps2) / cfg.chiN) - 1.0) * cs / cfg.ds);/**< Eq.(37) in Hansen (2016) */

  cma->state = CMAES_STATE_UPDATED;

  __NOTE__("%s\n", "end");
}


static inline void status_report(const struct cmaes_status cma, const struct cmaes_config cfg, char *file)
{
  __NOTE__("%s\n", "start");

  FILE *fp;
  char filename[128], date[64];
  sprintf(filename, "%s/%s_cmaes_gen%.3d.dat", DATAFOLDER, file, cma.gen);
  fp = fopen(filename, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
  getPresentDateInStrings(date);

  fprintf(fp, "# status report of generation %d (%d parameters, %d children):\n", cma.gen, cfg.ndim, cfg.lambda);
  fprintf(fp, "#\tbest score is %e (%d evaluations)\n", cma.bestScore, cma.count);
  fprintf(fp, "#\trange of function values is %e\n", cma.funcRange);
  fprintf(fp, "#\trange of function values in recent generations is %e\n", cma.histRange);
  fprintf(fp, "#\tcurrent guess is (%e", cma.xbest[0]);
  for(int ii = 1; ii < cfg.ndim; ii++)
    fprintf(fp, ", %e", cma.xbest[ii]);
  fprintf(fp, ")\n");
  fprintf(fp, "#\tstep size is %e\n", cma.sigma);

  fprintf(fp, "#\n#\tsampled population in this generation:\n");
  for(int ii = 0; ii < cfg.lambda; ii++){
    const int ll = cma.sorted[ii].idx;
    fprintf(fp, "%e", cma.score[ll]);
    for(int jj = 0; jj < cfg.ndim; jj++)
      fprintf(fp, " %e", cma.xx[INDEX2D(cfg.lambda, cfg.ndim, ll, jj)]);
    fprintf(fp, "\n");
  }/* for(int ii = 0; ii < cfg.lambda; ii++){ */
  fclose(fp);

  __NOTE__("%s\n", "end");
}


/**
 * @fn terminate
 *
 * @brief terminate optimization when satisfactory results are obtained
 * @detail see cmaes_TestForTermination in dm:~/work/181025cmaes/c-cmaes/src/cmaes.c: line 1535--1669
 */
bool terminate
(struct cmaes_status *cma, struct cmaes_config cfg
#ifdef  ENHANCE_STATUS_REPORT
 , char *file
#endif//ENHANCE_STATUS_REPORT
 )
{
  __NOTE__("%s\n", "start");

  const int ndim = cfg.ndim;
  const int gen = cma->gen;
  const int lambda = cfg.lambda;
  bool satisfied = false;

  /** function value converged */
  cma->funcRange = cma->score[cma->sorted[lambda - 1].idx] - cma->score[cma->sorted[0].idx];
  if( cma->funcRange < cfg.ftol ){
#ifndef SUPPRESS_STATUS_REPORT
    fprintf(stdout, "Differences of function (generation %d): range = %e while tolerance = %e\n", gen, cma->funcRange, cfg.ftol);
#endif//SUPPRESS_STATUS_REPORT
    satisfied = true;
  }/* if( cma->funcRange < cfg.ftol ){ */

  /** evolution of function value converged */
  const int numHist = (gen < cfg.maxHistory) ? gen : cfg.maxHistory;
  double maxHist = cma->bestScore;
  double minHist = DBL_MAX;
  for(int ii = 0; ii < numHist; ii++){
    const double val = cma->history[ii];
    maxHist = fmax(maxHist, val);
    minHist = fmin(minHist, val);
  }/* for(int ii = 0; ii < numHist; ii++){ */
  cma->histRange = maxHist - minHist;
  if( numHist == cfg.maxHistory ){
    if( cma->histRange < cfg.htol ){
#ifndef SUPPRESS_STATUS_REPORT
      fprintf(stdout, "History of function value (generation %d): range = %e while tolerance = %e\n", gen, cma->histRange, cfg.htol);
#endif//SUPPRESS_STATUS_REPORT
      satisfied = true;
    }/* if( cma->histRange < cfg.htol ){ */
  }/* if( numHist == cfg.maxHistory ){ */

  /** degree of variable changes */
  int unchanged = 0;
  for(int ii = 0; ii < ndim; ii++)
    unchanged += ((cma->sigma * fmax(sqrt(cma->mat[INDEX2D(ndim, ndim, ii, ii)]), cma->pc[ii])) < cfg.xtol);
  if( unchanged == ndim ){
#ifndef SUPPRESS_STATUS_REPORT
    fprintf(stdout, "Variable change is negligible in generation %d (all changes are smaller than %e)\n", gen, cfg.xtol);
#endif//SUPPRESS_STATUS_REPORT
    satisfied = true;
  }/* if( unchanged == ndim ){ */

  /** maximal function evaluations */
  if( cma->count > cfg.maxEval ){
    fprintf(stdout, "Too many function evaluations (generation %d): %d evaluations while preset maximal is %d\n", gen, cma->count, cfg.maxEval);
    satisfied = true;
  }/* if( cma->count > cfg.maxEval ){ */

  /** maximal iterations */
  if( gen > cfg.maxIter ){
    fprintf(stdout, "Too many iterations: %d generations while preset maximal is %d\n", gen, cfg.maxIter);
    satisfied = true;
  }/* if( gen > cfg.maxIter ){ */


#ifdef  ENHANCE_STATUS_REPORT
  /** report the tentative results */
  status_report(*cma, cfg, file);
#endif//ENHANCE_STATUS_REPORT

  __NOTE__("%s\n", "end");
  return (satisfied);
}


void config(const int argc, char **argv, struct cmaes_config *cfg, double **weight)
{
  __NOTE__("%s\n", "start");

  /** set problem dimension */
  requiredCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "ndim", &cfg->ndim));
  const int ndim = cfg->ndim;

  /** set number of sampling points in a generation */
  if( optionalCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "lambda", &cfg->lambda)) != myUtilAvail )
    cfg->lambda = 0;
  if( cfg->lambda < 2 )
    cfg->lambda = 4 + (int)floor(3.0 * log((double)ndim));/**< Eq.(48) in Hansen (2016) */
  const int lambda = cfg->lambda;

  /** set number of sampling points selected to calculate new mean of search distribution */
  if( optionalCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "mu", &cfg->mu)) != myUtilAvail )
    cfg->mu = -1;
  if( cfg->mu < 0 )
    cfg->mu = lambda >> 1;
  const int mu = cfg->mu;

  /** set learning rate */
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "cm", &cfg->cm)) != myUtilAvail )
    cfg->cm = 1.0;/**< Eq.(54) in Hansen (2016) */

  /** set weight to calculated a weighted average of mu (<= lambda) selected points (1st half) */
  *weight = (double *)malloc(sizeof(double) * lambda);  if( *weight == NULL ){    __KILL__(stderr, "ERROR: failure to allocate weight\n");  }
  cfg->weight = *weight;
  for(int ii = 0; ii < lambda; ii++){
    const double val = log((double)((1 + lambda) >> 1)) - log((double)(1 + ii));/**< Eq.(49) in Hansen (2016) */
    cfg->weight[ii] = val;
  }/* for(int ii = 0; ii < lambda; ii++){ */
  double s1 = 0.0;
  double s2 = 0.0;
  for(int ii = 0; ii < mu; ii++){
    const double val = cfg->weight[ii];
    s1 += val;
    s2 += val * val;
  }/* for(int ii = 0; ii < mu; ii++){ */
  cfg->mueff = s1 * s1 / s2;/**< Eq.(8) in Hansen (2016) */
  s2 = 1.0 / s1;
  for(int ii = 0; ii < mu; ii++)
    cfg->weight[ii] *= s2;/**< normalize weight (Eq. (7) in Hansen (2016) */

  /** set parameters for covariance matrix adaptation */
  cfg->cc  = (cfg->mueff + (double)(ndim << 2)) / (2.0 * cfg->mueff + (double)ndim * (double)(4 + ndim));/**< Eq.(56) in Hansen (2016) */
  cfg->c1  = 2.0 / (cfg->mueff + (1.3 + (double)ndim) * (1.3 + (double)ndim));/**< Eq.(57) in Hansen (2016) */
  cfg->cmu = fmin(1.0 - cfg->c1, 2.0 * (1.0 + cfg->mueff * (cfg->mueff - 2.0)) / (cfg->mueff * (cfg->mueff + (double)(2 + ndim) * (double)(2 + ndim))));/**< Eq.(58) in Hansen (2016) */

  /** set weight to calculated a weighted average of mu (<= lambda) selected points (2nd half) */
  s1 = 0.0;
  s2 = 0.0;
  for(int ii = mu; ii < lambda; ii++){
    const double val = fabs(cfg->weight[ii]);
    s1 += val;
    s2 += val * val;
  }/* for(int ii = mu; ii < lambda; ii++){ */
  s2 = fmin(1.0 + cfg->c1 / cfg->cmu, 1.0 + (2.0 * s1 * s1 / s2) / (2.0 + cfg->mueff));
  s2 = fmin(s2, (1.0 - cfg->c1 - cfg->cmu) / (cfg->cmu * (double)ndim)) / s1;
  for(int ii = mu; ii < lambda; ii++)
    cfg->weight[ii] *= s2;

  /** set a parameter for covariance matrix adaptation */
  s1 = 0.0;
  for(int ii = 0; ii < lambda; ii++)
    s1 += cfg->weight[ii];
  cfg->coeff = 1.0 - cfg->c1 - cfg->cmu * s1;

  if( (mu < 1) || (mu > lambda) || (mu == lambda) || (cfg->weight[0] == cfg->weight[mu - 1]) ){
    __KILL__(stderr, "ERROR: invalid setting of mu (%d), lambda (%d), or weight [%e, %e]\n", mu, lambda, cfg->weight[0], cfg->weight[mu - 1]);
  }/* if( (mu < 1) || (mu > lambda) || (mu == lambda) || (cfg->weight[0] == cfg->weight[mu - 1]) ){ */

  /** set parameters for step-size control */
  cfg->cs = (2.0 + cfg->mueff) / (5.0 + (double)ndim + cfg->mueff);/**< Eq.(55) in Hansen (2016) */
  cfg->ds = cfg->cs + 1.0 + 2.0 * fmax(0.0, sqrt((cfg->mueff - 1.0) / (double)(1 + ndim)) - 1.0);/**< Eq.(55) in Hansen (2016) */

  /** set additional parameters */
  cfg->chiN = sqrt((double)ndim) * (1.0 - 0.25 / (double)ndim + 1.0 / (21.0 * (double)ndim * (double)ndim));
  cfg->maxHistory = 10 + (int)ceil(30.0 * (double)ndim / (double)lambda);

  /** set tolerance parameters */
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "ftol", &cfg->ftol)) != myUtilAvail )
    cfg->ftol = 1.0e-12;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "htol", &cfg->htol)) != myUtilAvail )
    cfg->htol = 1.0e-12;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "xtol", &cfg->xtol)) != myUtilAvail )
    cfg->xtol = 0.0;

  if( optionalCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "maxEval", &cfg->maxEval)) != myUtilAvail )
    cfg->maxEval = 900 * (3 + ndim) * (3 + ndim);
  if( optionalCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "maxIter", &cfg->maxIter)) != myUtilAvail )
    cfg->maxIter = (int)ceil((double)cfg->maxEval / (double)lambda);

  __NOTE__("%s\n", "end");
}


void allocate_cmaes
(const struct cmaes_config cfg, struct cmaes_status *cma,
 double **mat, double **offd, double **tau, double **diag, double **orth,
 double **xx, double **xmean, double **xmold, double **xbest, double **xevol, double **xtemp,
 double **ps, double **pc, double **score, double **history, struct dbl_int **sorted)
{
  __NOTE__("%s\n", "start");

  const int ndim = cfg.ndim;

  *mat  = (double *)malloc(sizeof(double) * ndim * ndim);  if( *mat == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  mat\n");  }  cma-> mat = * mat;
  *offd = (double *)malloc(sizeof(double) * (ndim - 1));  if( *offd == NULL ){    __KILL__(stderr, "ERROR: failure to allocate offd\n");  }  cma->offd = *offd;
  *tau  = (double *)malloc(sizeof(double) * (ndim - 1));  if( * tau == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  tau\n");  }  cma-> tau = * tau;
  *diag = (double *)malloc(sizeof(double) * ndim       );  if( *diag == NULL ){    __KILL__(stderr, "ERROR: failure to allocate diag\n");  }  cma->diag = *diag;
  *orth = (double *)malloc(sizeof(double) * ndim * ndim);  if( *orth == NULL ){    __KILL__(stderr, "ERROR: failure to allocate orth\n");  }  cma->orth = *orth;

  const int lambda = cfg.lambda;
  *xx = (double *)malloc(sizeof(double) * ndim * lambda);  if( *xx == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xx\n");  }  cma->xx = *xx;
  *xmean = (double *)malloc(sizeof(double) * ndim);  if( *xmean == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xmean\n");  }  cma->xmean = *xmean;
  *xmold = (double *)malloc(sizeof(double) * ndim);  if( *xmold == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xmold\n");  }  cma->xmold = *xmold;
  *xbest = (double *)malloc(sizeof(double) * ndim);  if( *xbest == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xbest\n");  }  cma->xbest = *xbest;
  *xevol = (double *)malloc(sizeof(double) * ndim);  if( *xevol == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xevol\n");  }  cma->xevol = *xevol;
  *xtemp = (double *)malloc(sizeof(double) * ndim);  if( *xtemp == NULL ){    __KILL__(stderr, "ERROR: failure to allocate xtemp\n");  }  cma->xtemp = *xtemp;

  *pc = (double *)malloc(sizeof(double) * ndim);  if( *pc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pc\n");  }  cma->pc = *pc;
  *ps = (double *)malloc(sizeof(double) * ndim);  if( *ps == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ps\n");  }  cma->ps = *ps;
  *score = (double *)malloc(sizeof(double) * lambda);  if( *score == NULL ){    __KILL__(stderr, "ERROR: failure to allocate score\n");  }  cma->score = *score;
  *sorted = (struct dbl_int *)malloc(sizeof(struct dbl_int) * lambda);  if( *sorted == NULL ){    __KILL__(stderr, "ERROR: failure to allocate sorted\n");  }  cma->sorted = *sorted;
  *history = (double *)malloc(sizeof(double) * cfg.maxHistory);  if( *history == NULL ){    __KILL__(stderr, "ERROR: failure to allocate history\n");  }  cma->history = *history;

  __NOTE__("%s\n", "end");
}


void initialize_cmaes(const int argc, char **argv, const struct cmaes_config cfg, struct cmaes_status *cma, double *xmean_init)
{
  __NOTE__("%s\n", "start");

  /** initialize flags */
  cma->initial_phase = true;
  cma->eigen_uptodate = false;

  /** initialize status */
  cma->gen = 0;
  cma->count = 0;
  cma->state = CMAES_STATE_INITIAL;

  /** initialize score */
  cma->bestScore = DBL_MAX;
  cma->funcRange = DBL_MAX;
  cma->histRange = DBL_MAX;
  const int lambda = cfg.lambda;
  for(int ii = 0; ii < lambda; ii++){
    cma->score[ii] = DBL_MAX;
    cma->sorted[ii].val = DBL_MAX;
    cma->sorted[ii].idx = -1;
  }/* for(int ii = 0; ii < lambda; ii++){ */
  const int maxHistory = cfg.maxHistory;
  for(int ii = 0; ii < maxHistory; ii++)
    cma->history[ii] = DBL_MAX;

  /** initialize step size */
  double sigma;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "sigma", &sigma)) != myUtilAvail )
    sigma = 0.3;
  double trace = 0.0;
  const int ndim = cfg.ndim;
  for(int ii = 0; ii < ndim; ii++)
    trace += sigma * sigma;
  cma->sigma = sqrt(trace / (double)ndim);

  /** initialize matrices and vectors */
  for(int ii = 0; ii < ndim * ndim; ii++)
    cma->mat[ii] = cma->orth[ii] = 0.0;
  for(int ii = 0; ii < ndim; ii++){
    const int idx = INDEX2D(ndim, ndim, ii, ii);
    cma->orth[idx] = 1.0;
    const double val = sigma * sqrt((double)ndim / trace);
    cma-> mat[idx] = val * val;
    cma->diag[ii] = val;
    cma->pc[ii] = 0.0;
    cma->ps[ii] = 0.0;
  }/* for(int ii = 0; ii < ndim; ii++){ */

  /** initialize mean of sampling points: add option to specify a point?? */
  if( xmean_init == NULL ){
    char *init;
    if( optionalCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv, "init", &init)) == myUtilAvail ){
      FILE *fp;
      fp = fopen(init, "r");
      if( fp == NULL ){      __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", init);    }
      bool success = true;
      for(int ii = 0; ii < ndim; ii++){
	double val;
	success &= (1 == fscanf(fp, "%le", &val));
	cma->xmean[ii] = cma->xmold[ii] = val;
      }/* for(int ii = 0; ii < ndim; ii++){ */
      fclose(fp);
      if( !success ){
	__KILL__(stderr, "ERROR: failure to read \"%s\"\n", init);
      }/* if( !success ){ */
    }/* if( optionalCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv, "init", &init)) == myUtilAvail ){ */
    else
      for(int ii = 0; ii < ndim; ii++)
	cma->xmean[ii] = cma->xmold[ii] = 0.0;
  }/* if( xmean_init == NULL ){ */
  else
    for(int ii = 0; ii < ndim; ii++)
      cma->xmean[ii] = cma->xmold[ii] = xmean_init[ii];

  update_eigen(cma, cfg, true);

  __NOTE__("%s\n", "end");
}
