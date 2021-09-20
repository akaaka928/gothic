/**
 * @file nws.c
 *
 * @brief Data analysis for North-Western Stream
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2018/07/24 (Tue)
 *
 * Copyright (C) 2018 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "macro.h"


/**
 * @fn extract_stream_fields
 *
 * @brief Extract stellar particles locate in Stream fields by Komiyama et al. (2018)
 *
 * @param (aa) vector A
 * @param (bb) vector B
 */
void extract_stream_fields(const int head, const int num, real * restrict xi, real * restrict eta, real * restrict dist)
{
  __NOTE__("%s\n", "start");


  const real NE_slope = -2.784422092017102 ;  const real NE_intercept = -9.34058388367902;
  const real SW_slope = -2.8917360434134034;  const real SW_intercept = -11.981949647745681;

  const real hsc_fov_r2 = 0.5625;/**< = 0.75^2 */
  const real xi0_022 = -4.842739589150247;  const real eta0_022 = 2.290423176281251;
  const real xi0_009 = -5.530826768687439;  const real eta0_009 = 5.816784796173433;

  for(int ii = head; ii < head + num; ii++){
    /** stars must locate in between NE and SW edges */
    if( (eta[ii] <= (NE_intercept + NE_slope * xi[ii])) && (eta[ii] >= (SW_intercept + SW_slope * xi[ii])) ){
      int field = -1;

      if( eta[ii] <= 5.101772155015601 ){
	/* 2, 3, 4, or exterior */
	if( eta[ii] > 3.7516521430870844 )
	  field = 2;
	else{
	  /* 3, 4, or exterior */
	  if( eta[ii] > 2.612617574895511 )
	    field = 3;
	  else{
	    /* 4 or exterior: if 4, contained in 022 */
	    const real xx =  xi[ii] -  xi0_022;
	    const real yy = eta[ii] - eta0_022;
	    if( (xx * xx + yy * yy) <= hsc_fov_r2 )
	      field = 4;
	  }/* else{ */
	}/* else{ */
      }/* if( eta[ii] <= 5.101772155015601 ){ */
      else{
	/* 1 or exterior: if 1, contained in 009 or 004 */
	if( eta[ii] <= 6.061673446044403 )
	  field = 1;
	else{
	  const real xx =  xi[ii] -  xi0_009;
	  const real yy = eta[ii] - eta0_009;
	  if( (xx * xx + yy * yy) <= hsc_fov_r2 )
	    field = 1;
	}/* else{ */
      }/* else{ */


      if( field != -1 ){

	/* make histogram for each field */

	/* for visualization: small bin width */


	/* for on-the-fly analysis: very sparse (3 bins) */


      }/* if( field != -1 ){ */


    }/* if( (eta[ii] <= (NE_intercept + NE_slope * xi[ii])) && (eta[ii] >= (SW_intercept + SW_slope * xi[ii])) ){ */
  }/* for(int ii = head; ii < head + num; ii++){ */

  __NOTE__("%s\n", "end");
}




/* PFS: Field of View is ~ 1.25 deg^2 (hexagonal) */
/* for simplicity, set circles having similar area */
/* radius is 0.63 degree */



/* line-of-sight velocity of GCs */
/* center: locations of GCs */
/* radius: ??? same with Kirihara-kun?? size of GC?? width of NW Stream?? */



