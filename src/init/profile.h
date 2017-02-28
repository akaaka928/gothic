/**
 * @file profile.h
 *
 * @brief Header file for describing radial profile of spherical component(s)
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/02/24 (Fri)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef PROFILE_H
#define PROFILE_H


#include <stdbool.h>

#include "macro.h"



/**
 * @def ERFC_SMOOTHING
 * enable: adopt complementary error function based smoother
 * disable: adopt tangent hyperbolic based smoother
 */
#define ERFC_SMOOTHING


/** CFGFOLDER must same with CFGFOLDER defined in sample.c */
#define CFGFOLDER "cfg"

/** macros to specify the density distribution model */
/** positive value indicates spherical component(s) */
#define CENTRALBH (1000)
#define   PLUMMER  (0)
#define      KING  (1)
#define   BURKERT  (2)
#define HERNQUIST  (3)
#define       NFW  (4)
#define     MOORE  (5)
#define   EINASTO  (6)
#define TWO_POWER  (7)
#define TRI_POWER  (8)
#define  APP_KING (10)
#define APP_EVANS (11)
#define TABLE_RHO (20)
#define TABLE_SIG (21)
#define SPHSERSIC (30)
#define SIGTWOPOW (31)
/** negative value indicates disk component(s) */
#define  EXP_DISK  (-1)
#define    SERSIC  (-2)
#define  TBL_DISK (-10)


/** 2 * 2 bins are added in the both edge */
#define NRADBIN (      4194304  )
#define MINRAD  (1.0 / 1048576.0)


/**
 * @def MAKE_COLUMN_DENSITY_PROFILE
 * activate estimation for column density profile of spherical component(s)
 */
#define MAKE_COLUMN_DENSITY_PROFILE
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
#define SKIP_INTERVAL_FOR_COLUMN_DENSITY (128)
#endif//MAKE_COLUMN_DENSITY_PROFILE


/**
 * @struct profile
 *
 * @brief structure for mass distribution
 */
typedef struct
{
  double rad;
  double rho    , enc    , psi;
  double rho_tot, enc_tot, psi_tot;
  double drho_dr, d2rho_dr2;
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
  double Sigma;
#endif//MAKE_COLUMN_DENSITY_PROFILE
} profile;


/**
 * @struct profile_cfg
 *
 * @brief structure for mass distribution
 */
typedef struct
{
  char file[128];
  char table[128];/**< parameter for reading density profile in table form */
  double Mtot, rs;/**< parameters for all profiles */
  double einasto_alpha;/**< parameter for Einasto profile */
  double king_W0, king_rt, king_c;/**< parameter for King sphere */
  double twopower_alpha, twopower_beta, twopower_gamma;/**< parameters for two-power model */
  double tripower_delta, tripower_epsilon, tripower_rout;/**< additional parameters for three-power model */
  double alevans_alpha, alevans_beta, alevans_rc, alevans_rt, alevans_wt;/**< parameters for approximated lowered Evans model */
  double zd, Sigma0, vdispR0, vdispz0, vdisp_frac;/**< parameters for disk component(s) */
  double n_sersic, b_sersic;/**< parameters for Sersic profile */
  double vcirc_Rd, vcirc_max, vcirc_max_R, toomre, Qmin0, Qmin1, Qmin2, qminR0, qminR1, qminR2;/**< properties of disk component(s) */
  double retrogradeFrac;/**< fraction of retrograding disk particles */
  double rc, rc_width;
  double Ecut;
  ulong num;
  int forceNum;/**< parameter to specify number of N-body particles for the component */
  int kind;
  bool passed;/**< variable to estimate Toomre's Q-value */
  bool cutoff;
} profile_cfg;


/* list of functions appeared in ``profile.c'' */
void setDensityProfilePlummer  (profile *prf, const double rs);
void setDensityProfileBurkert  (profile *prf, const double rs);
void setDensityProfileHernquist(profile *prf, const double rs);
void setDensityProfileNFW      (profile *prf, const double rs);
void setDensityProfileMoore    (profile *prf, const double rs);

void setDensityProfileEinasto(profile *prf, const double rs, const double alpha);
void setDensityProfileAppKing(profile *prf, const double rs, const double rt);

void setDensityProfileTwoPower(profile *prf, const double rs, const double alpha, const double beta, const double gam);
void setDensityProfileTriPower(profile *prf, const double rin, const double rout, const double alp, const double bet, const double gam, const double del, const double eps);

void setDensityProfileAppLoweredEvans(profile *prf, const double rs, const double alpha, const double rc, const double beta, const double rt, const double invdelta);

void setContributionByCentralBH(profile *prf, const profile_cfg cfg);

void integrateDensityProfile(profile *prf, const double logrbin, const double Mtot, const bool cutoff, const double redge, const double width);

void readProfileCfg(char *fcfg, int *unit, int *kind, profile_cfg **cfg);

#ifdef  MAKE_COLUMN_DENSITY_PROFILE
void calcColumnDensityProfile(const int skind, profile **prf, const double logrmax, profile_cfg *cfg);
#endif//MAKE_COLUMN_DENSITY_PROFILE


#endif//PROFILE_H
