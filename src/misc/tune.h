/**
 * @file tune.h
 *
 * @brief Header file for auto-tuning of GOTHIC
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/02/28 (Tue)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef TUNE_H
#define TUNE_H


#include <stdbool.h>


/* macros to specify the auto-tuning mode */
#define USE_PARABOLIC_GROWTH_MODEL
#ifdef  USE_PARABOLIC_GROWTH_MODEL
/* this macro must be switched ON */
#define WALK_TREE_USE_REDUCED_CHISQ
/* this macro is optional */
#define USE_ADDITIONAL_SWITCH
#endif//USE_PARABOLIC_GROWTH_MODEL


#define FORCE_ADJUSTING_PARTICLE_TIME_STEPS


/**
 * @struct rebuildTree
 *
 * @brief structure for auto-tuning about tree rebuild intervals
 */
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

/**
 * @struct measuredTime
 *
 * @brief structure to record the measured execution time of various functions
 */
typedef struct
{
  /* counters for automatic tree rebuild */
  double walkTree[2];
  double makeTree;
#ifndef SERIALIZED_EXECUTION
  /* counter for automatic load balancing, reset when particle exchanging */
  double sum_excg;
  /* counters for detecting slow-down due to particle mixing, reset when tree rebuilding */
  double sum_rebuild, excg;
#endif//SERIALIZED_EXECUTION
} measuredTime;


/**
 * @struct statVal
 *
 * @brief structure for fitting with the least squares method
 */
typedef struct
{
  double S, Sx, Sy, Sxx, Sxy, Syy;
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  double Sxxxx, Sxxx, Sxxy;
#endif//USE_PARABOLIC_GROWTH_MODEL
} statVal;

/**
 * @struct guessTime
 *
 * @brief structure for guessing the execution time based on the fitting
 */
typedef struct
{
  double slope, icept, rchisq;
  double time;/* *= scale */
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  double second;
#endif//USE_PARABOLIC_GROWTH_MODEL
} guessTime;


/**
 * @struct autoTuningParam
 *
 * @brief structure for summarizing various estimations of the execution time
 */
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


/* static const double tolerance4chisq = 0.125; */
static const double tolerance4chisq = 1.5625e-2;


/* list of functions appeared in ``tune.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  void initStatVal(statVal *val);
  void initGuessTime(guessTime *model);

  void linearModel(guessTime *model, statVal *val, const double steps, const double twalk, const double reduce);
  void  powerModel(guessTime *model, statVal *val, const double steps, const double twalk, const double reduce);
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  void parabolicModel(guessTime *model, statVal *val, const double steps, const double twalk, const double reduce);
#endif//USE_PARABOLIC_GROWTH_MODEL
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//TUNE_H
