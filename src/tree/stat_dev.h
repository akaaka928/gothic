/*************************************************************************\
 *                                                                       *
                  last updated on 2016/10/28(Fri) 16:59:13
 *                                                                       *
 *    Header File to examine tree statistics                             *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef STAT_DEV_H
#define STAT_DEV_H
//-------------------------------------------------------------------------
#include <macro.h>
#include <cudalib.h>
//-------------------------------------------------------------------------
#ifdef  LOCALIZE_I_PARTICLES
#       include "../sort/peano.h"
#endif//LOCALIZE_I_PARTICLES
//-------------------------------------------------------------------------
#include "../tree/make.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* #define NUM_IGROUP_SAFETY_FACTOR (8) */
#define NUM_IGROUP_SAFETY_FACTOR (16)
//-------------------------------------------------------------------------
#ifdef  LOCALIZE_I_PARTICLES
//-------------------------------------------------------------------------
#ifndef NEIGHBOR_PHKEY_LEVEL
/* #define NEIGHBOR_PHKEY_LEVEL (0) */
#define NEIGHBOR_PHKEY_LEVEL (1)
/* #define NEIGHBOR_PHKEY_LEVEL (2) */
/* #define NEIGHBOR_PHKEY_LEVEL (3) */
/* #define NEIGHBOR_PHKEY_LEVEL (4) */
#endif//NEIGHBOR_PHKEY_LEVEL
//-------------------------------------------------------------------------
#   if  NEIGHBOR_PHKEY_LEVEL > (MAXIMUM_PHKEY_LEVEL >> 1)
#undef  NEIGHBOR_PHKEY_LEVEL
#define NEIGHBOR_PHKEY_LEVEL   (MAXIMUM_PHKEY_LEVEL >> 1)
#endif//NEIGHBOR_PHKEY_LEVEL > (MAXIMUM_PHKEY_LEVEL >> 1)
//-------------------------------------------------------------------------
#endif//LOCALIZE_I_PARTICLES
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
#   if  !defined(LOCALIZE_I_PARTICLES) || !defined(BLOCK_TIME_STEP)
#undef  USE_VARIABLE_NEIGHBOR_LEVEL
#endif//!defined(LOCALIZE_I_PARTICLES) || !defined(BLOCK_TIME_STEP)
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
//-------------------------------------------------------------------------
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
#define MINIMUM_NEIGHBOR_PHKEY_LEVEL (1)
#define MAXIMUM_NEIGHBOR_PHKEY_LEVEL (8)
typedef struct
{
  PHint jump;
  real dt;
  int num;
} histogram_dt;
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-- List of functions appeared in "stat_dev.cu"
//-------------------------------------------------------------------------
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
  //-----------------------------------------------------------------------
/*   void   freeParticleGroups(laneinfo  *laneInfo_hst, laneinfo  *laneInfo_dev, real  *laneTime_dev */
/* #ifdef  USE_VARIABLE_NEIGHBOR_LEVEL */
/* 			    , real  *dt_dev, real  *dt_hst, histogram_dt  *dtInfo */
/* #endif//USE_VARIABLE_NEIGHBOR_LEVEL */
/* 			    ); */
/*   muse  allocParticleGroups(laneinfo **laneInfo_hst, laneinfo **laneInfo_dev, real **laneTime_dev */
/* #ifdef  USE_VARIABLE_NEIGHBOR_LEVEL */
/* 			    , real **dt_dev, real **dt_hst, histogram_dt **dtInfo, int *dtInfo_num */
/* #endif//USE_VARIABLE_NEIGHBOR_LEVEL */
/* 			    , int *inumPerLane, int *maxNgrp, const int num_max, deviceProp devProp); */
  void  freeParticleGroups(laneinfo  *laneInfo_hst, laneinfo  *laneInfo_dev, double  *laneTime_dev
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
			   , real  *dt_dev, real  *dt_hst, histogram_dt  *dtInfo
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
			   );
  muse allocParticleGroups
  (laneinfo **laneInfo_hst, laneinfo **laneInfo_dev, double **laneTime_dev
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
   , real **dt_dev, real **dt_hst, histogram_dt **dtInfo, int *dtInfo_num
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
   , int *inumPerLane, int *maxNgrp, const int num_max, deviceProp devProp);
  //-----------------------------------------------------------------------
/*   void updateParticleGroups */
/*   (const int Ni, laneinfo *laneInfo_hst, const int inumPerLane, const int maxNgrp, int *Ngrp, PHint *peano, PHinfo * RESTRICT level, int *bottomLev */
/* #ifdef  USE_VARIABLE_NEIGHBOR_LEVEL */
/*    , const iparticle body, real *dt_dev, real *dt_hst, int *gnum, histogram_dt **dtInfo */
/* #endif//USE_VARIABLE_NEIGHBOR_LEVEL */
/*    ); */
  void updateParticleGroups
  (const int Ni, laneinfo *laneInfo, const int inumPerLane, const int maxNgrp, int *Ngrp
#ifdef  LOCALIZE_I_PARTICLES
   , PHint *peano
#endif//LOCALIZE_I_PARTICLES
#ifdef  USE_VARIABLE_NEIGHBOR_LEVEL
   , const iparticle body, real *dt_dev, real *dt_hst, int *gnum, histogram_dt **dtInfo
#endif//USE_VARIABLE_NEIGHBOR_LEVEL
   );
  void commitParticleGroups(const int Ngrp, laneinfo *laneInfo_hst, laneinfo *laneInfo_dev);
  //-----------------------------------------------------------------------
#ifdef  __CUDACC__
}
#endif//__CUDACC__
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//STAT_DEV_H
//-------------------------------------------------------------------------
