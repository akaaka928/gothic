/*************************************************************************\
 *                                                                       *
                  last updated on 2017/01/17(Tue) 19:46:53
 *                                                                       *
 *    Header File for applying auto-tuning                               *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef TUNE_H
#define TUNE_H
//-------------------------------------------------------------------------
#include <stdbool.h>
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
#define USE_PARABOLIC_GROWTH_MODEL
#else///WALK_TREE_COMBINED_MODEL
/* choose one of listed below: */
#define WALK_TREE_TOTAL_SUM_MODEL
/* #define WALK_TREE_ARITHMETIC_PROGRESSION_MODEL */
/* #define WALK_TREE_GEOMETRIC_PROGRESSION_MODEL */
#endif//WALK_TREE_COMBINED_MODEL
//-------------------------------------------------------------------------
#ifdef  USE_PARABOLIC_GROWTH_MODEL
/* this macro must be switched ON */
#define WALK_TREE_USE_REDUCED_CHISQ
/* this macro is optional */
#define USE_ADDITIONAL_SWITCH
#endif//USE_PARABOLIC_GROWTH_MODEL
//-------------------------------------------------------------------------
#define FORCE_ADJUSTING_PARTICLE_TIME_STEPS
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* auto-tuning related structures and functions */
//-------------------------------------------------------------------------
typedef struct
{
  double interval;
#   if  defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
  double avg, var;
#endif//defined(FORCE_ADJUSTING_PARTICLE_TIME_STEPS) && defined(BLOCK_TIME_STEP)
  int reuse;
#ifdef  BLOCK_TIME_STEP
  bool adjust;
#endif//BLOCK_TIME_STEP
} rebuildTree;
//-------------------------------------------------------------------------
typedef struct
{
  /* counters for automatic tree rebuild */
  double walkTree[2];
  double makeTree;
#ifdef  WALK_TREE_TOTAL_SUM_MODEL
  double incSum;
#endif//WALK_TREE_TOTAL_SUM_MODEL
#ifndef SERIALIZED_EXECUTION
  /* counter for automatic load balancing, reset when particle exchanging */
  double sum_excg;
/*   double genTree, calcAcc, calcMAC; */
/* #ifdef  MONITOR_LETGEN_TIME */
/*   double makeLET; */
/* #endif//MONITOR_LETGEN_TIME */
  /* counters for detecting slow-down due to particle mixing, reset when tree rebuilding */
  double sum_rebuild, excg;
#endif//SERIALIZED_EXECUTION
} measuredTime;
//-------------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
typedef struct
{
  double S, Sx, Sy, Sxx, Sxy, Syy;
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  double Sxxxx, Sxxx, Sxxy;
#endif//USE_PARABOLIC_GROWTH_MODEL
} statVal;
//-------------------------------------------------------------------------
typedef struct
{
  double slope, icept, rchisq;
  double time;/* *= scale */
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  double second;
#endif//USE_PARABOLIC_GROWTH_MODEL
} guessTime;
#endif//WALK_TREE_COMBINED_MODEL
//-------------------------------------------------------------------------
typedef struct
{
  statVal   linearStats, powerStats;
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  statVal parabolicStats;
#endif//USE_PARABOLIC_GROWTH_MODEL
  guessTime linearGuess, powerGuess;
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  guessTime parabolicGuess;
#endif//USE_PARABOLIC_GROWTH_MODEL
} autoTuningParam;
//-------------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
/* static const double tolerance4chisq = 0.125; */
static const double tolerance4chisq = 1.5625e-2;
#endif//WALK_TREE_COMBINED_MODEL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "tune.c"
//-------------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
  void initStatVal(statVal *val);
  /* void calcStatVal(statVal *val, const double xx, const double yy, const double sigma); */
  //-----------------------------------------------------------------------
  void initGuessTime(guessTime *model);
  /* void evalStatVal(guessTime *model, const statVal val); */
  //-----------------------------------------------------------------------
  void linearModel(guessTime *model, statVal *val, const double steps, const double twalk, const double reduce);
  void  powerModel(guessTime *model, statVal *val, const double steps, const double twalk, const double reduce);
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  void parabolicModel(guessTime *model, statVal *val, const double steps, const double twalk, const double reduce);
#endif//USE_PARABOLIC_GROWTH_MODEL
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------
#endif//WALK_TREE_COMBINED_MODEL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//TUNE_H
//-------------------------------------------------------------------------
