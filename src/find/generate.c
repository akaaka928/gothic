/**
 * @file generate.c
 *
 * @brief Generate parameter sets using CMA-ES (Evolution Strategy with Covariance Matrix Adaptation)
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

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "name.h"
#include "myutil.h"
#include "timer.h"
#include "rand.h"

#include "../find/cmaes.h"
#include "../find/cmaes_io.h"


/** list of parameters: */
/** (0) mass of DM component (log_10 M_DM) */
/** (1) scale radius of DM component (log_10 r_s) */
/** (2) mass of stellar component (log_10 M_star): [10^5 Msun : infty) */
/** (3) scale radius of stellar component (log_10 r_0) */
/** (4) dimensionless King parameter of stellar component (W_0): (0 : infty) */
/** (5) infalling velocity of the satellite (v_r): (-1.2 vesc : 0] */
/** (6) tangential velocity of the satellite (v_t): [0 : 1.2 vesc) */
/** (7) azimuthal angle of the satellite (theta): [0 : M_PI / 2] */
/** (8) direction of the velocity vector (v_theta): [0 : 2 M_PI) */


static inline void evaluate(struct cmaes_status *cma, const struct cmaes_config cfg)
{
  if( cma->state == CMAES_STATE_SAMPLED ){
    const int ndim = cfg.ndim;
    const int lambda = cfg.lambda;

    /**# generate initial-conditions */
    /**# run N-body simulations and read the final score */

    for(int ii = 0; ii < lambda; ii++){
      /** minimum value = 0.2 @ x0 = 0.45, x1 = -0.69, x2 = 0.83, and x3 = 0.23 or 2.23 */
      const double x0 = cma->xx[INDEX2D(lambda, ndim, ii, 0)] - 0.45;
      const double x1 = cma->xx[INDEX2D(lambda, ndim, ii, 1)] + 0.69;
      const double x2 = cma->xx[INDEX2D(lambda, ndim, ii, 2)] - 0.83;
      const double x3 = cma->xx[INDEX2D(lambda, ndim, ii, 3)] - 1.23;

      const double y0 = cost * x1 - sint * x2;
      const double y1 = sint * x1 + cost * x2;

      const double r2 = y0 * y0 + 25.0 * y1 * y1;
      const double xx_2 = x3 * x3 - 1.0;

      const double x4 = cma->xx[INDEX2D(lambda, ndim, ii, 3)] - 2.23;
      cma->score[ii] = (2.0 - exp(-fabs(x0))) * (3.0 - exp(-r2)) * (0.1 + xx_2 * xx_2 + fmin(0.1, x4 * x4));
    }/* for(int ii = 0; ii < lambda; ii++){ */

    cma->state = CMAES_STATE_SCORED;
  }/* if( cma->state == CMAES_STATE_SAMPLED ){ */
}


/**
 * @fn optimize
 *
 * @brief optimize parameter space
 * @detail see cmaes_Optimize in dm:~/work/181025cmaes/c-cmaes/src/cmaes.c: line 758--792
 *
 * @return (cmaes) status of CMA-ES
 */
static inline void optimize(struct cmaes_status *cma, struct cmaes_config cfg, struct cmaes_rand *gauss, rand_state *rand, char *file
#ifdef  USE_HDF5_FORMAT
			    , const hid_t H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
			    , const hid_t H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
			    )
{
  __NOTE__("%s\n", "start");

  while( true ){
    /** generate new search points */
    sampling(cma, cfg, gauss, rand
#ifdef  SUPPORT_CMAES_RESTART
	     , file
#ifdef  USE_HDF5_FORMAT
	     , H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
	     , H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
#endif//SUPPORT_CMAES_RESTART
	     );

    /** set parameter region on some parameters */
    for(int ii = 0; ii < lambda; ii++){
      while( true ){
	const double W0 = cma->xx[INDEX2D(lambda, ndim, ii, 4)];/**< dimensionless King parameter of stellar component */
	const double vr = cma->xx[INDEX2D(lambda, ndim, ii, 5)];/**< infalling velocity of the satellite */
	const double vt = cma->xx[INDEX2D(lambda, ndim, ii, 6)];/**< tangential velocity of the satellite */
	const double theta = cma->xx[INDEX2D(lambda, ndim, ii, 7)];/**< azimuthal angle of the satellite */
	const double vtheta = cma->xx[INDEX2D(lambda, ndim, ii, 8)];/**< direction of the velocity vector */




	const double    mu = cma->xx[INDEX2D(lambda, ndim, ii, 0)];/**< location of the stream axis */
	const double sigma = cma->xx[INDEX2D(lambda, ndim, ii, 1)];/**< width of the stream */
	const double  star = cma->xx[INDEX2D(lambda, ndim, ii, 2)];/**< number count of stars */
	const double slope = cma->xx[INDEX2D(lambda, ndim, ii, 3)];/**< slope of foreground/background stars */
	const double noise = cma->xx[INDEX2D(lambda, ndim, ii, 4)];/**< floor of foreground/background stars */
	if( (W0 > 0.0) &&
	    (vr >= -1.2 * vesc) && (vr <= 0.0) &&
	    (vt <=  1.2 * vesc) && (vt >= 0.0) &&
	    (theta >= 0.0) && (theta <= 0.5 * M_PI) &&
	    (vtheta >= 0.0) && (vtheta < 2.0 * M_PI)
	    )
	  break;
	resampling(cma, cfg, gauss, rand, ii);
      }/* while( true ){ */
    }/* for(int ii = 0; ii < lambda; ii++){ */


    /** if resampling works, then the intermediate results must be updated */




    /** evaluate the new search points */
    evaluate(cma, cfg);

    /** update the search distribution */
    update(cma, cfg);

    /** check whether the solution is found */
    if( terminate(cma, cfg
#ifdef  ENHANCE_STATUS_REPORT
		  , file
#endif//ENHANCE_STATUS_REPORT
		  ) )
      break;
  }/* while( true ){ */

  /** dump final results */
  saveStatus(file, cfg, *cma, *gauss, *rand
#ifdef  USE_HDF5_FORMAT
	     , H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
	     , H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
	     );

  __NOTE__("%s\n", "end");
}


int main(int argc, char **argv)
{
  /** configure the details of the numerical simulation */
  /** read command line arguments */
  if( argc < 3 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 3);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -restart=<int> (optional)\n");
    __FPRINTF__(stderr, "          -ndim=<int>\n");
    __FPRINTF__(stderr, "          -lambda=<int> (optional) -mu=<int> (optional)\n");
    __FPRINTF__(stderr, "          -cm=<double> (optional)\n");
    __FPRINTF__(stderr, "          -sigma=<double> (optional)\n");
    __FPRINTF__(stderr, "          -init=<char *>\n");
    __FPRINTF__(stderr, "          -ftol=<double> (optional) -htol=<double> (optional) -xtol=<double> (optional)\n");
    __FPRINTF__(stderr, "          -maxEval=<int> (optional) -maxIter=<int> (optional)\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 3 ){ */
  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv, "file", &file));


  /** configure the CMA-ES procedure */
  static int gen = 0;
  static struct cmaes_config cfg;
  double *weight;
  optionalCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "restart", &gen));
  if( gen > 0 )
    loadConfig(file, &cfg, &weight);
  else
    config(argc, argv, &cfg, &weight);


  /** prepare the CMA-ES status */
  static struct cmaes_status cma;
  double *mat, *offd, *tau, *diag, *orth;
  double *xx, *xmean, *xmold, *xbest, *xevol, *xtemp;
  double *ps, *pc;
  double *score, *history;
  struct dbl_int *sorted;
  allocate_cmaes(cfg, &cma, &mat, &offd, &tau, &diag, &orth, &xx, &xmean, &xmold, &xbest, &xevol, &xtemp, &ps, &pc, &score, &history, &sorted);


  /** initialize pseudo random-number */
  static struct cmaes_rand gauss;
  rand_state *rand;
  initRandNum(&rand);


  /** prepare for file I/O using HDF5 */
#ifdef  USE_HDF5_FORMAT
  hid_t H5T_NATIVE_DBLINT;
#ifdef  USE_SFMT
  hid_t H5T_NATIVE_SFMT;
#endif//USE_SFMT
  createHDF5type(&H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
		 , &H5T_NATIVE_SFMT
#endif//USE_SFMT
		 );
#endif//USE_HDF5_FORMAT


  /** configure the CMA-ES status */
  if( gen > 0 )
    loadStatus(file, &cma, &gauss, rand, gen
#ifdef  USE_HDF5_FORMAT
	       , H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
	       , H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
	       );
  else
    initialize_cmaes(argc, argv, cfg, &cma, NULL);


  /** print out configuration of CMA-ES */
  if( gen == 0 ){
    FILE *fp;
    char filename[128], date[64];
    saveConfig(file, cfg);
    sprintf(filename, "%s/%s_cmaes.cfg", DOCUMENTFOLDER, file);
    fp = fopen(filename, "w");
    if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }
    getPresentDateInStrings(date);
    fprintf(fp, "# Configuration of CMA-ES (Covariance Matrix Adaptation Evolution Strategies)\n");
    fprintf(fp, "# \tgenerated on %s", date);
    fprintf(fp, "#####\n");
    fprintf(fp, "%d # cfg.ndim\n", cfg.ndim);
    fprintf(fp, "%d # cfg.lambda\n", cfg.lambda);
    fprintf(fp, "%d # cfg.mu\n", cfg.mu);
    fprintf(fp, "%e # cfg.mueff; %e indicates a reasonable setting of weight\n", cfg.mueff, 0.25 * cfg.lambda);
    fprintf(fp, "#####\n");
    fprintf(fp, "%d # cfg.maxHistory\n", cfg.maxHistory);
    fprintf(fp, "#####\n");
    fprintf(fp, "%e # cfg.cs; %e might be a better choice\n", cfg.cs, 4.0 / (double)cfg.ndim);
    fprintf(fp, "%e # cfg.ds\n", cfg.ds);
    fprintf(fp, "#####\n");
    fprintf(fp, "%e # cfg.cm; c_m < 1 can be advantageous on noisy functions\n", cfg.cm);
    fprintf(fp, "%e # cfg.cc\n", cfg.cc);
    fprintf(fp, "%e # cfg.c1; %e might be a better choice\n", cfg.c1, 2.0 / (double)(cfg.ndim * cfg.ndim));
    fprintf(fp, "%e # cfg.cmu; %e might be a better choice\n", cfg.cmu, fmin(cfg.mueff / (double)(cfg.ndim * cfg.ndim), 1.0 - cfg.c1));
    fprintf(fp, "%e # 1 - c_1 - c_mu * sum w_j; close or equal to 0 is expected\n", cfg.coeff);
    fprintf(fp, "#####\n");
    fprintf(fp, "%e # cfg.ftol\n", cfg.ftol);
    fprintf(fp, "%e # cfg.htol\n", cfg.htol);
    fprintf(fp, "%e # cfg.xtol\n", cfg.xtol);
    fprintf(fp, "%d # cfg.maxEval\n", cfg.maxEval);
    fprintf(fp, "%d # cfg.maxIter\n", cfg.maxIter);
    fprintf(fp, "#####\n");
    fclose(fp);
  }/* if( gen == 0 ){ */


  /** iteration until finding the solution */
  optimize(&cma, cfg, &gauss, rand, file
#ifdef  USE_HDF5_FORMAT
	   , H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
	   , H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
	   );


#ifdef  USE_HDF5_FORMAT
  removeHDF5type(H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
		 , H5T_NATIVE_SFMT
#endif//USE_SFMT
		 );
#endif//USE_HDF5_FORMAT


  free(cma.score);  free(cma.sorted);  free(cma.history);
  free(cma.pc);  free(cma.ps);
  free(cma.xmean);  free(cma.xmold);  free(cma.xbest);  free(cma.xevol);  free(cma.xtemp);
  free(cma.xx);
  free(cma.diag);  free(cma.orth);
  free(cma.mat);  free(cma.offd);  free(cma.tau);
  free(cfg.weight);


  return (0);
}
