/**
 * @file generate.c
 *
 * @brief Generate parameter sets using CMA-ES (Evolution Strategy with Covariance Matrix Adaptation)
 *
 * @author Yohei MIKI (University of Tokyo)
 *
 * @date 2020/02/11 (Tue)
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

#include <unistd.h>
#include <sys/stat.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "name.h"
#include "myutil.h"
#include "timer.h"
#include "rand.h"
#include "constants.h"

#include "../find/cmaes.h"
#include "../find/cmaes_io.h"

#include "../anal/observation.h"

static const int unit_system = GalacticScale;

#define SURVEYGEN "gen"
#define SURVEYTAG "run"

#define CUSP_MODEL


/** list of parameters: */
/** (0) mass of DM component (log_10 M_DM): (-infty : infty) */
/** (1) scale radius of DM component (log_10 r_s): (-infty : infty) */
/** (2) mass of stellar component (log_10 M_star): (-infty : infty) */
/** (3) scale radius of stellar component (log_10 r_0): (-infty : infty) */
/** (4) dimensionless King parameter of stellar component (log_10 W_0): (0.4 <= W_0 <= 12.0) */
/** (5) polar angle of the satellite (theta): [0 : M_PI / 2] */
/** (6) velocity of the satellite in radial direction (v_r / 100): (infty : 0] */
/** (7) velocity of the satellite along polar angle (v_t / 100): (-infty : infty) */
/** (8) velocity of the satellite along azimuthal angle (v_p / 100): (-infty : infty) */


extern const double mass_astro2com;
extern const double rho_crit;
extern const double length2astro;


static inline void evaluate(struct cmaes_status *cma, const struct cmaes_config cfg, char *file)
{
  if( cma->state == CMAES_STATE_SAMPLED ){
    const int gen = cma->gen;
    const int ndim = cfg.ndim;
    const int lambda = cfg.lambda;

    char folder[128];
    sprintf(folder, "%s/%s/%s%d", CFGFOLDER, file, SURVEYGEN, gen);
    struct stat stat_buf;
    if( stat(folder, &stat_buf) != 0 )
      mkdir(folder, 0775);

    char summary[128];
    sprintf(summary, "%s/summary.txt", folder);
    bool prepared = (access(summary, F_OK) == 0);
    FILE *fp;
    if( prepared ){
      fp = fopen(summary, "r");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", summary);      }
      int tmp_gen = - 1;
      int tmp_ndim = 0;
      int tmp_lambda = 0;

      int checker = 1;
      checker &= (1 == fscanf(fp, "%d", &tmp_gen));
      checker &= (1 == fscanf(fp, "%d", &tmp_ndim));
      checker &= (1 == fscanf(fp, "%d", &tmp_lambda));

      if( (!checker) || (tmp_gen != gen) || (tmp_ndim != ndim) || (tmp_lambda != lambda) )
	prepared = false;

      fclose(fp);
    }

    if( !prepared ){
      /* write cfg files for MAGI and editor */
      for(int ii = 0; ii < lambda; ii++){
	const double logM_dm = cma->xx[INDEX2D(lambda, ndim, ii, 0)];
	const double logrs_dm = cma->xx[INDEX2D(lambda, ndim, ii, 1)];
	const double logM_star = cma->xx[INDEX2D(lambda, ndim, ii, 2)];
	const double logr0_star = cma->xx[INDEX2D(lambda, ndim, ii, 3)];
	const double logW0_star = cma->xx[INDEX2D(lambda, ndim, ii, 4)];

	char cfg_magi[256];
	char ini_dark[256];
	char ini_star[256];
	sprintf(cfg_magi, "%s/run%d-magi.cfg", folder, ii);
	sprintf(ini_dark, "%s/run%d-dark.ini", folder, ii);
	sprintf(ini_star, "%s/run%d-star.ini", folder, ii);
	char run_file[64];
	sprintf(run_file, "%s-gen%d-run%d", file, gen, ii);

	fp = fopen(cfg_magi, "w");
	if( fp == NULL ){	  __KILL__(stderr, "ERROR: failure to open \"%s\"\n", cfg_magi);	}
	fprintf(fp, "%d\n", unit_system);
	fprintf(fp, "%d\n", 2);/* number of components */
#ifdef  CUSP_MODEL
	fprintf(fp, "%d\t../%s\t%d\t%d\n", 4, ini_dark, 0, 0);/**< 4 means NFW model */
#else///CUSP_MODEL
	fprintf(fp, "%d\t../%s\t%d\t%d\n", 2, ini_dark, 0, 0);/**< 2 means Burkert model */
#endif//CUSP_MODEL
	fprintf(fp, "%d\t../%s\t%d\t%d\n", 1, ini_star, 0, 0);/**< 1 means King model */
	fclose(fp);

	fp = fopen(ini_dark, "w");
	if( fp == NULL ){	  __KILL__(stderr, "ERROR: failure to open \"%s\"\n", ini_dark);	}
	const double M200 = pow(10.0, logM_dm);
	const double rs = pow(10.0, logrs_dm);
	fprintf(fp, "%e\n", M200);/**< mass of DM component */
	fprintf(fp, "%e\n", rs);/**< scale radius of DM component */
	fprintf(fp, "%d\n", 1);/**< set explicit cutoff at r_200 */
	const double r200 = cbrt(3.0 * (M200 * mass_astro2com) / (800.0 * M_PI * rho_crit)) * length2astro;
	/* fprintf(stdout, "r200 = %e, M200 = %e, mass_astro2com = %e, rho_crit = %e, length2astro = %e\n", r200, M200, mass_astro2com, rho_crit, length2astro); */
	fprintf(fp, "%e %e\n", r200, fmax(0.1 * r200, rs));
	fclose(fp);

	fp = fopen(ini_star, "w");
	if( fp == NULL ){	  __KILL__(stderr, "ERROR: failure to open \"%s\"\n", ini_star);	}
	fprintf(fp, "%e\n", pow(10.0, logM_star));/**< mass of stellar component */
	fprintf(fp, "%e\n", pow(10.0, logr0_star));/**< core radius of stellar component */
	fprintf(fp, "%e\n", pow(10.0, logW0_star));/**< dimensionless King parameter of stellar component */
	fprintf(fp, "%d\n", 0);/**< does not set explicit cutoff */
	fclose(fp);

	char cfg_edit[256];
	char ini_edit[256];
	sprintf(cfg_edit, "%s/run%d-edit.cfg", folder, ii);
	sprintf(ini_edit, "%s/run%d-edit.ini", folder, ii);
	fp = fopen(cfg_edit, "w");
	if( fp == NULL ){	  __KILL__(stderr, "ERROR: failure to open \"%s\"\n", cfg_edit);	}
	fprintf(fp, "%d\n", unit_system);
	fprintf(fp, "%d\n", 1);/* number of input files */
	fprintf(fp, "%s%d-run%d-magi ../%s\n", SURVEYGEN, gen, ii, ini_edit);
	fclose(fp);

	fp = fopen(ini_edit, "w");
	if( fp == NULL ){	  __KILL__(stderr, "ERROR: failure to open \"%s\"\n", ini_edit);	}
	fprintf(fp, "0\n");/**< specify the rotation matrix (actually, no rotation will be performed) */
	fprintf(fp, "1.0 0.0 0.0\n");
	fprintf(fp, "0.0 1.0 0.0\n");
	fprintf(fp, "0.0 0.0 1.0\n");
	/** calculate initial location and velocity of the satellite */
	const double r0 = initial_separation;/**< = 10 rs */
	const double theta = cma->xx[INDEX2D(lambda, ndim, ii, 5)];
	const double phi = 0.0;
	const double vr = velocity_normalization * cma->xx[INDEX2D(lambda, ndim, ii, 6)];
	const double vt = velocity_normalization * cma->xx[INDEX2D(lambda, ndim, ii, 7)];
	const double vp = velocity_normalization * cma->xx[INDEX2D(lambda, ndim, ii, 8)];
	const double sint = sin(theta);
	const double cost = cos(theta);
	const double sinp = sin(phi);
	const double cosp = cos(phi);
	const double X0 = r0 * sint * cosp;
	const double Y0 = r0 * sint * sinp;
	const double Z0 = r0 * cost;
	fprintf(fp, "%e %e %e\n", X0, Y0, Z0);
	const double vX = sint * cosp * vr + cost * cosp * vt - sinp * vp;
	const double vY = sint * sinp * vr + cost * sinp * vt + cosp * vp;
	const double vZ = cost        * vr - sint        * vt            ;
	fprintf(fp, "%e %e %e\n", vX, vY, vZ);
	fprintf(fp, "0\n");/**< number of components to be removed */
	fclose(fp);

	register_model(run_file, logM_dm, logrs_dm, logM_star, logr0_star, logW0_star, theta, vr, vt, vp);
      }

      fp = fopen(summary, "w");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", summary);      }
      fprintf(fp, "%d\n", gen);
      fprintf(fp, "%d\n", ndim);
      fprintf(fp, "%d\n", lambda);
      fprintf(fp, "#mean value and its 1 sigma range of the input parameters:\n");
      fprintf(fp, "\tM_dm = %e Msun; (%e:%e)\n", pow(10.0, cma->xmean[0]), pow(10.0, cma->xmean[0] - cma->sigma), pow(10.0, cma->xmean[0] + cma->sigma));
      fprintf(fp, "\trs_dm = %e kpc; (%e:%e)\n", pow(10.0, cma->xmean[1]), pow(10.0, cma->xmean[1] - cma->sigma), pow(10.0, cma->xmean[1] + cma->sigma));
      fprintf(fp, "\tM_star = %e Msun; (%e:%e)\n", pow(10.0, cma->xmean[2]), pow(10.0, cma->xmean[2] - cma->sigma), pow(10.0, cma->xmean[2] + cma->sigma));
      fprintf(fp, "\tr0_star = %e kpc; (%e:%e)\n", pow(10.0, cma->xmean[3]), pow(10.0, cma->xmean[3] - cma->sigma), pow(10.0, cma->xmean[3] + cma->sigma));
      fprintf(fp, "\tW0_star = %e; (%e:%e)\n", pow(10.0, cma->xmean[4]), pow(10.0, cma->xmean[4] - cma->sigma), pow(10.0, cma->xmean[4] + cma->sigma));
      fprintf(fp, "\ttheta = %e; (%e:%e)\n", cma->xmean[5], cma->xmean[5] - cma->sigma, cma->xmean[5] + cma->sigma);
      fprintf(fp, "\tv_rad = %e km/s; (%e:%e)\n", velocity_normalization * cma->xmean[6], velocity_normalization * (cma->xmean[6] - cma->sigma), velocity_normalization * (cma->xmean[6] + cma->sigma));
      fprintf(fp, "\tv_theta = %e km/s; (%e:%e)\n", velocity_normalization * cma->xmean[7], velocity_normalization * (cma->xmean[7] - cma->sigma), velocity_normalization * (cma->xmean[7] + cma->sigma));
      fprintf(fp, "\tv_phi = %e km/s; (%e:%e)\n", velocity_normalization * cma->xmean[8], velocity_normalization * (cma->xmean[8] - cma->sigma), velocity_normalization * (cma->xmean[8] + cma->sigma));
      fclose(fp);

      __KILL__(stderr, "ERROR: generate initial conditions and run N-body simulations for generation %d (%d offsprings)\n", gen, lambda);
    }
    else{
#ifdef  USE_DOUBLE_PRECISION
      const hid_t H5T_GOTHIC_REAL = H5T_NATIVE_DOUBLE;
#else///USE_DOUBLE_PRECISION
      const hid_t H5T_GOTHIC_REAL = H5T_NATIVE_FLOAT;
#endif//USE_DOUBLE_PRECISION

      int completed = 0;
      /**# run N-body simulations and read the final score */
      /** query the score of N-body simulations of this generation (presence of the score of ALL runs) */
      for(int ii = 0; ii < lambda; ii++){
	char filename[256];
	sprintf(filename, "%s/%s-%s%d-%s%d_%s.h5", DATAFOLDER, file, SURVEYGEN, gen, SURVEYTAG, ii, FILE_TAG_SCORE);
	hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	hid_t group = H5Gopen(target, GROUP_TAG_INFO, H5P_DEFAULT);
	hid_t attr;
	attr = H5Aopen(group, ATTR_TAG_SCORE, H5P_DEFAULT);
	chkHDF5err(H5Aread(attr, H5T_GOTHIC_REAL, &(cma->score[ii])));
	chkHDF5err(H5Aclose(attr));

	double finish;
	attr = H5Aopen(group, ATTR_TAG_FINISH_TIME, H5P_DEFAULT);
	chkHDF5err(H5Aread(attr, H5T_NATIVE_DOUBLE, &finish));
	chkHDF5err(H5Aclose(attr));
	double time;
	attr = H5Aopen(group, ATTR_TAG_ANALYZED_TIME, H5P_DEFAULT);
	chkHDF5err(H5Aread(attr, H5T_NATIVE_DOUBLE, &time));
	chkHDF5err(H5Aclose(attr));

	completed += (time >= finish);

	chkHDF5err(H5Gclose(group));
	chkHDF5err(H5Fclose(target));
      }

      if( completed == lambda )
	cma->state = CMAES_STATE_SCORED;
    }


    if( cma->state != CMAES_STATE_SCORED ){
      __KILL__(stderr, "ERROR: N-body simulations for generation %d (%d offsprings) are not yet completed\n", gen, lambda);
    }

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

  const int ndim = cfg.ndim;
  const int lambda = cfg.lambda;

  /* Representative quantities at r = 10 rs: */
  /* Escape velocity            v_esc is 2.966149e+01 (= 2.900277e+02 km / s) */
  const double vmax = 360.0;/**< = 1.2 v_esc @ r = 10 r_s */
  const double vmax2 = vmax * vmax;

  /* /\** set the minimum mass of the NW stream *\/ */
  /* extern const double mass_astro2com; */
  /* NWstream_mass = POW(TEN, 5.0) * mass_astro2com; */
  const double logMstar_min = 5.0;

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
    if( cma->state == CMAES_STATE_SAMPLED ){
#ifdef  SUPPORT_CMAES_RESTART
      bool resampled = false;
#endif//SUPPORT_CMAES_RESTART
      for(int ii = 0; ii < lambda; ii++){
	while( true ){
	  const double theta = fabs(fmod(cma->xx[INDEX2D(lambda, ndim, ii, 5)], M_PI_2));/**< polar angle of the satellite */
	  cma->xx[INDEX2D(lambda, ndim, ii, 5)] = theta;

	  const double M200 = pow(10.0, cma->xx[INDEX2D(lambda, ndim, ii, 0)]) + pow(10.0, cma->xx[INDEX2D(lambda, ndim, ii, 2)]);
	  const double r200 = cbrt(3.0 * (M200 * mass_astro2com) / (800.0 * M_PI * rho_crit)) * length2astro;

	  const double rs = pow(10.0, cma->xx[INDEX2D(lambda, ndim, ii, 1)]);
	  const double r0 = pow(10.0, cma->xx[INDEX2D(lambda, ndim, ii, 3)]);

	  const double logMstar = cma->xx[INDEX2D(lambda, ndim, ii, 2)];
	  const double W0 = pow(10.0, cma->xx[INDEX2D(lambda, ndim, ii, 4)]);

	  const double vr = velocity_normalization * cma->xx[INDEX2D(lambda, ndim, ii, 6)];/**< infalling velocity of the satellite */
	  const double vt = velocity_normalization * cma->xx[INDEX2D(lambda, ndim, ii, 7)];/**< infalling velocity of the satellite */
	  const double vp = velocity_normalization * cma->xx[INDEX2D(lambda, ndim, ii, 8)];/**< infalling velocity of the satellite */

	  if( (vr <= 0.0) &&/**< consider only infalling satellite */
	      (logMstar >= logMstar_min) &&/**< consider only satellite whose mass is greater than the NW stream */
	      /* (M200 <= 1.0e+9) &&/\**< M200 is smaller than 10^9 Msun *\/ */
	      /* (rs >= 0.3) && (r0 >= 0.3) &&/\**< rs >= 300 pc and r0 >= 300 pc *\/ */
	      (rs < r200) && (r0 < r200) &&/**< C200 >= 1 for stellar and DM components */
	      ((vr * vr + vt * vt + vp * vp) <= vmax2) &&/**< v_max = 1.2 v_esc */
	      (W0 >= 0.4) && (W0 <= 12.0)/**< remove unrealistic concentration for King sphere */
	      )
	    break;
	  resampling(cma, cfg, gauss, rand, ii);
#ifdef  SUPPORT_CMAES_RESTART
	  resampled |= true;
#endif//SUPPORT_CMAES_RESTART
	}/* while( true ){ */
      }/* for(int ii = 0; ii < lambda; ii++){ */

      /** if resampling works, then the intermediate results must be updated */
#ifdef  SUPPORT_CMAES_RESTART
      if( resampled )
	saveStatus(file, cfg, *cma, *gauss, *rand
#ifdef  USE_HDF5_FORMAT
		   , H5T_NATIVE_DBLINT
#ifdef  USE_SFMT
		   , H5T_NATIVE_SFMT
#endif//USE_SFMT
#endif//USE_HDF5_FORMAT
		   );
#endif//SUPPORT_CMAES_RESTART
    }


    /** evaluate the new search points */
    evaluate(cma, cfg, file);

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
    __FPRINTF__(stderr, "          -init=<char *> (optional)\n");
    __FPRINTF__(stderr, "          -ftol=<double> (optional) -htol=<double> (optional) -xtol=<double> (optional)\n");
    __FPRINTF__(stderr, "          -maxEval=<int> (optional) -maxIter=<int> (optional)\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 3 ){ */
  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv, "file", &file));

  setPhysicalConstantsAndUnitSystem(unit_system, 1);

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
