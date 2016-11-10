/************************************************************************* \
 *                                                                       *
                  last updated on 2016/10/28(Fri) 17:26:33
 *                                                                       *
 *    Utility functions to apply auto-tuning                             *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
//-------------------------------------------------------------------------
#include <macro.h>
//-------------------------------------------------------------------------
#include "../misc/tune.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef  WALK_TREE_COMBINED_MODEL
//-------------------------------------------------------------------------
void initStatVal(statVal *val)
{
  //-----------------------------------------------------------------------
  val->S   = 0.0;
  val->Sx  = 0.0;
  val->Sy  = 0.0;
  val->Sxx = 0.0;
  val->Sxy = 0.0;
  val->Syy = 0.0;
  //-----------------------------------------------------------------------
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  val->Sxxy  = 0.0;
  val->Sxxx  = 0.0;
  val->Sxxxx = 0.0;
#endif//USE_PARABOLIC_GROWTH_MODEL
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void calcStatVal(statVal *val, const double xx, const double yy, const double sigma)
{
  //-----------------------------------------------------------------------
  const double s2inv = 1.0 / (sigma * sigma);
  //-----------------------------------------------------------------------
  val->S   += s2inv;
  val->Sx  += s2inv * xx;
  val->Sy  += s2inv * yy;
  val->Sxx += s2inv * xx * xx;
  val->Sxy += s2inv * xx * yy;
  val->Syy += s2inv * yy * yy;
  //-----------------------------------------------------------------------
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  val->Sxxy  += s2inv * xx * xx * yy;
  val->Sxxx  += s2inv * xx * xx * xx;
  val->Sxxxx += s2inv * xx * xx * xx * xx;
#endif//USE_PARABOLIC_GROWTH_MODEL
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void initGuessTime(guessTime *model)
{
  //-----------------------------------------------------------------------
  model->slope  = 0.0;
  model->icept  = 0.0;
#ifdef  USE_PARABOLIC_GROWTH_MODEL
  model->second = 0.0;
#endif//USE_PARABOLIC_GROWTH_MODEL
  model->rchisq = DBL_MAX;
  model->time   = -1.0;
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
static inline void evalStatVal(guessTime *model, const statVal val)
{
  //-----------------------------------------------------------------------
  model->slope = (val.S * val.Sxy - val.Sx * val.Sy) / (val.S * val.Sxx - val.Sx * val.Sx);
  model->icept = (val.Sy - val.Sx * model->slope) / val.S;
  //-----------------------------------------------------------------------
  /* model->chisq = val.Syy + model->slope * (model->slope * val.Sxx - 2.0 * val.Sxy) - val.S * model->icept * model->icept; */
  model->rchisq = val.Syy - (model->slope * (model->slope * val.Sxx + model->icept * val.Sx) +
			     model->icept * (model->icept * val.S   + model->slope * val.Sx));
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void linearModel(guessTime *model, statVal *val, const double steps, const double twalk, const double reduce)
{
  //-----------------------------------------------------------------------
  calcStatVal(val, steps, twalk, tolerance4chisq * reduce * twalk);
  //-----------------------------------------------------------------------
  if( steps > 2.5 ){
    evalStatVal(model, *val);
#ifdef  WALK_TREE_USE_REDUCED_CHISQ
    model->rchisq /= (steps - 2.0);
#endif//WALK_TREE_USE_REDUCED_CHISQ
    model->time = 0.5 * model->slope * (steps + 0.5) * (steps + 0.5);
  }/* if( steps > 2.5 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
void  powerModel(guessTime *model, statVal *val, const double steps, const double twalk, const double reduce)
{
  //-----------------------------------------------------------------------
  calcStatVal(val, steps, log(twalk), tolerance4chisq * reduce);
  //-----------------------------------------------------------------------
  if( steps > 2.5 ){
    //---------------------------------------------------------------------
    evalStatVal(model, *val);
    //---------------------------------------------------------------------
    const double np1_2 = (steps + 0.5);
    const double rr = pow(M_E, model->slope);
    const double cc = pow(M_E, model->icept);
    const double rn = pow(M_E, model->slope * np1_2);
    //---------------------------------------------------------------------
    model->time = (1.0 + (-1.0 + np1_2 * model->slope) * rn) * cc * rr / (rr - 1.0);
    //---------------------------------------------------------------------
    if( rr < 1.0 )
      model->rchisq = DBL_MAX;
#ifdef  WALK_TREE_USE_REDUCED_CHISQ
    model->rchisq /= (steps - 2.0);
#endif//WALK_TREE_USE_REDUCED_CHISQ
    //---------------------------------------------------------------------
  }/* if( steps > 2.5 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifdef  USE_PARABOLIC_GROWTH_MODEL
void parabolicModel(guessTime *model, statVal *val, const double steps, const double twalk, const double reduce)
{
  //-----------------------------------------------------------------------
  calcStatVal(val, steps, twalk, tolerance4chisq * reduce * twalk);
  //-----------------------------------------------------------------------
  if( steps > 3.5 ){
    //---------------------------------------------------------------------
    const double tmp0 = val->S * val->Sxxx - val->Sxx * val->Sx;
    const double tmp1 = val->S * val->Sxx  - val->Sx  * val->Sx;
    const double tmp2 = val->S * val->Sxy  - val->Sx  * val->Sy;
    //---------------------------------------------------------------------
    model->second = (tmp0 * tmp2 + (val->Sy * val->Sxx - val->S * val->Sxxy) * tmp1) / ((val->S * val->Sxxxx - val->Sxx * val->Sxx) * tmp1 - tmp0 * tmp0);
    model->slope  = (tmp2 - tmp0 * model->second) / tmp1;
    model->icept  = (val->Sy - model->second * val->Sxx - model->slope * val->Sx) / val->S;
    //---------------------------------------------------------------------
    model->rchisq = val->Syy - (model->second * (model->second * val->Sxxxx + 2.0 * model->slope  * val->Sxxx) +
				model->slope  * (model->slope  * val->Sxx   + 2.0 * model->icept  * val->Sx  ) +
				model->icept  * (model->icept  * val->S     + 2.0 * model->second * val->Sxx ));
    //---------------------------------------------------------------------
    const double np1_2 = (steps + 0.5);
    const double aa = model->second;
    const double bb = model->slope - 2.0 * model->second;
    model->time = (2.0 * np1_2 * aa / 3.0 + 0.5 * (bb - aa)) * np1_2 * np1_2;
    //---------------------------------------------------------------------
    if( (bb + (-1.0 + 2.0 * np1_2) * aa) < 0.0 )
      model->rchisq = DBL_MAX;
    model->rchisq /= (steps - 3.0);
    //---------------------------------------------------------------------
  }/* if( steps > 3.5 ){ */
  //-----------------------------------------------------------------------
}
#endif//USE_PARABOLIC_GROWTH_MODEL
//-------------------------------------------------------------------------
#endif//WALK_TREE_COMBINED_MODEL
//-------------------------------------------------------------------------
