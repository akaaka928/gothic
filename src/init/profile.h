/*************************************************************************\
 *                                                                       *
                  last updated on 2016/10/24(Mon) 18:08:03
 *                                                                       *
 *    Header File to describe radial profile of spherical component(s)   *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef PROFILE_H
#define PROFILE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  !defined(_STDBOOL_H) && !defined(_STDBOOL)
#       include <stdbool.h>
#endif//!defined(_STDBOOL_H) && !defined(_STDBOOL)
//-------------------------------------------------------------------------
#ifndef MACRO_H
#      include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define ERFC_SMOOTHING
//-------------------------------------------------------------------------
/* CFGFOLDER must same with CFGFOLDER defined in sample.c */
#define CFGFOLDER "cfg"
//-------------------------------------------------------------------------
#define CENTRALBH (1000)
//-------------------------------------------------------------------------
#define   PLUMMER ( 0)
#define      KING ( 1)
#define   BURKERT ( 2)
#define HERNQUIST ( 3)
#define       NFW ( 4)
#define     MOORE ( 5)
#define   EINASTO ( 6)
#define TWO_POWER ( 7)
#define TRI_POWER ( 8)
#define  APP_KING (10)
#define APP_EVANS (11)
#define TABLE_RHO (20)
#define TABLE_SIG (21)
#define SPHSERSIC (30)
#define SIGTWOPOW (31)
//-------------------------------------------------------------------------
#define  EXP_DISK (-1)
#define    SERSIC (-2)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* 2 * 2 bins are added in the both edge */
#define NRADBIN (      4194304  )
#define MINRAD  (1.0 / 1048576.0)
//-------------------------------------------------------------------------
#define MAKE_COLUMN_DENSITY_PROFILE
//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------
typedef struct
{
  char file[128];
  char table[128];/* parameter for reading density profile in table form */
  double Mtot, rs;/* parameters for all profiles */
  double einasto_alpha;/* parameter for Einasto profile */
  double king_W0, king_rt, king_c;/* parameter for King sphere */
  double twopower_alpha, twopower_beta, twopower_gamma;/* parameters for two-power model */
  double tripower_delta, tripower_epsilon, tripower_rout;/* additional parameters for three-power model */
  double alevans_alpha, alevans_beta, alevans_rc, alevans_rt, alevans_wt;/* parameters for approximated lowered Evans model */
  double zd, Sigma0, vdispR0, vdispz0, vdisp_frac;/* parameters for disk components */
  double n_sersic, b_sersic;/* parameters for Sersic profile */
  double vcirc_Rd, vcirc_max, vcirc_max_R, toomre, Qmin;/* properties of disk components */
  bool passed;/* variable to estimate Toomre's Q-value */
  ulong num;
  int forceNum;/* parameter to specify number of N-body particles for the component */
  bool cutoff;
  double rc, rc_width;
  int kind;
  double Ecut;
} profile_cfg;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "profile.c"
//-------------------------------------------------------------------------
void setDensityProfilePlummer  (profile *prf, const double rs);
void setDensityProfileBurkert  (profile *prf, const double rs);
void setDensityProfileHernquist(profile *prf, const double rs);
void setDensityProfileNFW      (profile *prf, const double rs);
void setDensityProfileMoore    (profile *prf, const double rs);
//-------------------------------------------------------------------------
void setDensityProfileEinasto(profile *prf, const double rs, const double alpha);
void setDensityProfileAppKing(profile *prf, const double rs, const double rt);
//-------------------------------------------------------------------------
void setDensityProfileTwoPower(profile *prf, const double rs, const double alpha, const double beta, const double gam);
void setDensityProfileTriPower(profile *prf, const double rin, const double rout, const double alp, const double bet, const double gam, const double del, const double eps);
//-------------------------------------------------------------------------
void setDensityProfileAppLoweredEvans(profile *prf, const double rs, const double alpha, const double rc, const double beta, const double rt, const double invdelta);
//-------------------------------------------------------------------------
void setContributionByCentralBH(profile *prf, const profile_cfg cfg);
//-------------------------------------------------------------------------
void integrateDensityProfile(profile *prf, const double logrbin, const double Mtot, const bool cutoff, const double redge, const double width);
//-------------------------------------------------------------------------
void readProfileCfg(char *fcfg, int *unit, int *kind, profile_cfg **cfg);
//-------------------------------------------------------------------------
#ifdef  MAKE_COLUMN_DENSITY_PROFILE
void calcColumnDensityProfile(const int skind, profile **prf, const double logrmax, profile_cfg *cfg);
#endif//MAKE_COLUMN_DENSITY_PROFILE
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//PROFILE_H
//-------------------------------------------------------------------------
