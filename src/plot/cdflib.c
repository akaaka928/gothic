/**
 * @file cdflib.c
 *
 * @brief Common settings for error analysis of gravitational tree code
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2017/10/26 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "macro.h"

#include "cdflib.h"



void allocPercentile(double **per)
{
  __NOTE__("%s\n", "start");

  /** allocation */
  *per = (double *)malloc(sizeof(double) * NSUMMARY);
  if( *per == NULL ){    __KILL__(stderr, "ERROR: failure to allocate per");  }

  /** initialization */
#   if  NSUMMARY == 28
  (*per)[ 0] = 1.00;
  (*per)[ 1] = 0.99999943;/**< 5 sigma */
  (*per)[ 2] = 0.99993666;/**< 4 sigma */
  (*per)[ 3] = 0.99730020;/**< 3 sigma */
  (*per)[ 4] = 0.99;
  (*per)[ 5] = 0.95449974;/**< 2 sigma */
  (*per)[ 6] = 0.95;
  (*per)[ 7] = 0.90;
  (*per)[ 8] = 0.80;
  (*per)[ 9] = 0.75;
  (*per)[10] = 0.70;
  (*per)[11] = 0.68268949;/**< 1 sigma */
  (*per)[12] = 0.67;
  (*per)[13] = 0.60;
  (*per)[14] = 0.50;
  (*per)[15] = 0.40;
  (*per)[16] = 0.33;
  (*per)[17] = 0.31731051;/**< 1 sigma */
  (*per)[18] = 0.30;
  (*per)[19] = 0.25;
  (*per)[20] = 0.20;
  (*per)[21] = 0.10;
  (*per)[22] = 0.05;
  (*per)[23] = 0.04550026;/**< 2 sigma */
  (*per)[24] = 0.01;
  (*per)[25] = 0.00269980;/**< 3 sigma */
  (*per)[26] = 0.00006334;/**< 4 sigma */
  (*per)[27] = 0.00000057;/**< 5 sigma */
#endif//NSUMMARY == 28


  __NOTE__("%s\n", "end");
}


void  freePercentile(double  *per)
{
  __NOTE__("%s\n", "start");

  free(per);

  __NOTE__("%s\n", "end");
}
