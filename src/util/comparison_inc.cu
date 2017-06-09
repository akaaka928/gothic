/**
 * @file comparison_inc.cu
 *
 * @brief Source code for comparing values on GPU
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/04/05 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef COMPARISON_INC_CU
#define COMPARISON_INC_CU

#include <stdio.h>
#include <stdlib.h>
#include <helper_cuda.h>

#include "macro.h"


typedef struct __align__( 8){   float val;  int idx;}  floc;
typedef struct __align__(16){  double val;  int idx;}  dloc;
typedef struct __align__( 8){     int val;  int idx;}  iloc;
typedef struct __align__( 8){    uint val;  int idx;}  uloc;
typedef struct __align__(16){    long val;  int idx;}  lloc;
typedef struct __align__(16){   ulong val;  int idx;} ulloc;


__device__ __forceinline__  float getMinVal( float aa,  float bb){  return (fminf(aa, bb));}
__device__ __forceinline__  float getMaxVal( float aa,  float bb){  return (fmaxf(aa, bb));}
__device__ __forceinline__ double getMinVal(double aa, double bb){  return (fmin (aa, bb));}
__device__ __forceinline__ double getMaxVal(double aa, double bb){  return (fmax (aa, bb));}
__device__ __forceinline__   int getMinVal(  int aa,   int bb){  return ((aa <= bb) ? aa : bb);}
__device__ __forceinline__   int getMaxVal(  int aa,   int bb){  return ((aa >= bb) ? aa : bb);}
__device__ __forceinline__  uint getMinVal( uint aa,  uint bb){  return ((aa <= bb) ? aa : bb);}
__device__ __forceinline__  uint getMaxVal( uint aa,  uint bb){  return ((aa >= bb) ? aa : bb);}
__device__ __forceinline__  long getMinVal( long aa,  long bb){  return ((aa <= bb) ? aa : bb);}
__device__ __forceinline__  long getMaxVal( long aa,  long bb){  return ((aa >= bb) ? aa : bb);}
__device__ __forceinline__ ulong getMinVal(ulong aa, ulong bb){  return ((aa <= bb) ? aa : bb);}
__device__ __forceinline__ ulong getMaxVal(ulong aa, ulong bb){  return ((aa >= bb) ? aa : bb);}
#endif//COMPARISON_INC_CU
