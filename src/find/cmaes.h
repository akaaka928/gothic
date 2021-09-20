/**
 * @file cmaes.h
 *
 * @brief Header file for CMA-ES (Evolution Strategy with Covariance Matrix Adaptation)
 *
 * @author Yohei MIKI (University of Tokyo)
 *
 * @date 2019/02/12 (Tue)
 *
 * Copyright (C) 2019 Yohei MIKI
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef CMAES_H
#define CMAES_H

#include <stdbool.h>


/**
 * @struct cmaes_config
 *
 * @brief structure for CMA-ES configuration
 */
struct cmaes_config
{
  double *weight;/**< double[lambda]: weight to calculated a weighted average of mu (<= lambda) selected points */

  double mueff;/**< mu_eff in Hansen (2016) */
  double cm;/**< learning rate, c_m in Hansen (2016) */
  double cs;/**< c_sigma (learning rate for updating conjugate evolution path) in Hansen (2016) */
  double ds;/**< d_sigma (damping parameter for updating step size) in Hansen (2016) */
  double cc;/**< c_c (learning rate for updating evolution path) in Hansen (2016) */
  double c1;/**< c_1 (learning rate for updating the covariance matrix; rank-one update) in Hansen (2016) */
  double cmu;/**< c_mu (learning rate for updating the covariance matrix; rank-mu update) in Hansen (2016) */
  double coeff;/**< 1 - c_1 - c_mu \sum w_j in Hansen (2016) */
  double chiN;/**< expectation of ||N(0, I)|| */

  double ftol;/**< tolerance parameter for checking convergence of function values in a generation */
  double htol;/**< tolerance parameter for checking convergence of multiple generations */
  double xtol;/**< tolerance parameter for checking the degree of variable changes */
  int maxEval;/**< maximal function evaluations */
  int maxIter;/**< maximal iterations */

  int ndim;/**< problem dimension (# of individual parameters) */
  int lambda;/**< number of sampling points in a generation */
  int mu;/**< number of sampling points selected to calculate new mean of search distribution */

  int maxHistory;/**< 10 * ceil(30 * ndim / lambda) */
};


/**
 * @struct cmaes_rand
 *
 * @brief structure for random numbers in CMA-ES
 */
struct cmaes_rand
{
  double hold;
  bool stored;
};


struct dbl_int
{
  double val;
  int idx;
};


/**
 * @struct cmaes_status
 *
 * @brief structure for CMA-ES status
 */
struct cmaes_status
{
  double * mat;/**< double[ndim][ndim]:  covariance matrix C[ii][jj] */
  double *offd;/**< double[ndim - 1]: off-diagonal elements of the tridiagonal matrix */
  double * tau;/**< double[ndim - 1]: scalar factors of the elementary reflectors */

  double *diag;/**< double[ndim]: square root of eigenvalues; diagonal matrix D in Hansen (2016) */
  double *orth;/**< double[ndim][ndim]: matrix with normalized eigenvectors in columns; matrix B in Hansen (2016) */

  double *xx;/**< double[lambda][ndim]: x-vectors, lambda offspring */

  double *xmean;/**< double[ndim]: mean x vector, "parent" */
  double *xmold;/**< double[ndim]: mean x vector in the previous generation */
  double *xbest;/**< double[ndim]: x-vector which generates the best fitness value */
  double *xevol;/**< double[ndim]: evolution path of x-vector */
  double *xtemp;/**< double[ndim]: temporary vector used in different places */

  double *ps;/**< double[ndim]: p_sigma (conjugate evolution path) in Hansen (2016) */
  double *pc;/**< double[ndim]: p_c (evolution path) in Hansen (2016) */

  double *score;/**< double[lambda]: fitness value for each candidate solution */
  struct dbl_int *sorted;/**< struct dbl_int [lambda]: fitness value for each candidate solution in ascending order */
  double *history;/**< double [maxHistory = 10 * ceil(30 * ndim / lambda)]: history of best score in each generation */

  double bestScore;/**< best fitness value in all generations */
  double funcRange;/**< range of function values in a generation */
  double histRange;/**< range of function values in recent generations */

  double sigma;/**< step size */

  int gen;/**< generation number */
  int state;

  int count;/**< total number of function evaluations */

  bool eigen_uptodate;/**< a flag to indicate the eigensystem is up to date */
  bool initial_phase;
};
#define CMAES_STATE_INITIAL (0)
#define CMAES_STATE_SAMPLED (1)
#define CMAES_STATE_SCORED  (2)
#define CMAES_STATE_UPDATED (3)


#define chkLAPACKerr(err) __chkLAPACKerr(err, __FILE__, __LINE__)


/** list of functions appeared in ``cmaes.c'' */
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
 );
void resampling(struct cmaes_status *cma, struct cmaes_config cfg, struct cmaes_rand *gauss, rand_state *rand, const int ii);

void update(struct cmaes_status *cma, struct cmaes_config cfg);

bool terminate
(struct cmaes_status *cma, struct cmaes_config cfg
#ifdef  ENHANCE_STATUS_REPORT
 , char *file
#endif//ENHANCE_STATUS_REPORT
 );

void config(const int argc, char **argv, struct cmaes_config *cfg, double **weight);
void allocate_cmaes
(const struct cmaes_config cfg, struct cmaes_status *cma,
 double **mat, double **offd, double **tau, double **diag, double **orth,
 double **xx, double **xmean, double **xmold, double **xbest, double **xevol, double **xtemp,
 double **ps, double **pc, double **score, double **history, struct dbl_int **sorted);
void initialize_cmaes(const int argc, char **argv, const struct cmaes_config cfg, struct cmaes_status *cma, double *xmean_init);


#endif//CMAES_H
