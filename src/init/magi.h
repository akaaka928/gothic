/*************************************************************************\
 *                                                                       *
                  last updated on 2016/02/15(Mon) 10:45:02
 *                                                                       *
 *    Header File for Definition to generate initial condition           *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef MAGI_H
#define MAGI_H
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
#ifndef PROGRESS_REPORT_ON
#define PROGRESS_REPORT_ON
#endif//PROGRESS_REPORT_ON
//-------------------------------------------------------------------------
#define OUTPUT_ASCII_PROFILE
//-------------------------------------------------------------------------
/* 2 * 2 bins are added in the both edge */
#define NRADBIN (      4194304  )
#define MINRAD  (1.0 / 1048576.0)
typedef struct
{
  double rad;
  double rho    , enc    , psi;
  double rho_tot, enc_tot, psi_tot;
  double drho_dr, d2rho_dr2;
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
  double zd, Sigma0, vdispR0, vdispz0, vdisp_frac;/* parameters for Exponential disk */
  double n_sersic, b_sersic;/* parameters for Sersic profile */
  double vcirc_Rd, vcirc_max, vcirc_max_R, toomre, Qmin;/* properties of Exponential disk */
  ulong num;
  int forceNum;/* parameter to specify number of N-body particles for the component */
  bool cutoff;
  double rc, rc_width;
  int kind;
  double Ecut;
} profile_cfg;
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
#define NDIVIDE_GAUSSQD (16)
//-------------------------------------------------------------------------
/* #define NENEBIN (131072) */
/* #define NENEBIN (262144) */
/* #define NENEBIN (524288) */
#define NENEBIN (1048576)
typedef struct
{
  real ene, val;
} dist_func;
//-------------------------------------------------------------------------
#define NMAX_GAUSS_QD (51)
#define NTBL_GAUSS_QD ((NMAX_GAUSS_QD >> 1) + (NMAX_GAUSS_QD & 1))
#define NINTBIN NMAX_GAUSS_QD
//-------------------------------------------------------------------------
#define NKIND_MAX (8)
//-------------------------------------------------------------------------
/* #define N_SAFETY (10) */
//-------------------------------------------------------------------------
#define N_PRINT_LINES_ASCII (8192)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//MAGI_H
//-------------------------------------------------------------------------
