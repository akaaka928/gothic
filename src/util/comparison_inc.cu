/**
 * @file comparison_inc.cu
 *
 * @brief Source code for comparing values on GPU
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/10/26 (Thu)
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


template <typename Type> __device__ __forceinline__ void stloc(Type * smem, const int tidx, const Type val){  smem[tidx] = val;}

/* #   if  defined(__CUDACC_RTC__) */
/* #define __VECTOR_FUNCTIONS_DECL__ __host__ __device__ */
/* #else///defined(__CUDACC_RTC__) */
/* #define __VECTOR_FUNCTIONS_DECL__ static __inline__ __host__ __device__ */
/* #endif//defined(__CUDACC_RTC__) */
__device__ __forceinline__  floc make_floc ( float x, int y){   floc t;  t.val = x;  t.idx = y;  return t;}
__device__ __forceinline__  dloc make_dloc (double x, int y){   dloc t;  t.val = x;  t.idx = y;  return t;}
__device__ __forceinline__  iloc make_iloc (   int x, int y){   iloc t;  t.val = x;  t.idx = y;  return t;}
__device__ __forceinline__  uloc make_uloc (  uint x, int y){   uloc t;  t.val = x;  t.idx = y;  return t;}
__device__ __forceinline__  lloc make_lloc (  long x, int y){   lloc t;  t.val = x;  t.idx = y;  return t;}
__device__ __forceinline__ ulloc make_ulloc( ulong x, int y){  ulloc t;  t.val = x;  t.idx = y;  return t;}

__device__ __forceinline__  floc ldloc(volatile  floc & smem){  return make_floc (smem.val, smem.idx);}
__device__ __forceinline__  dloc ldloc(volatile  dloc & smem){  return make_dloc (smem.val, smem.idx);}
__device__ __forceinline__  iloc ldloc(volatile  iloc & smem){  return make_iloc (smem.val, smem.idx);}
__device__ __forceinline__  uloc ldloc(volatile  uloc & smem){  return make_uloc (smem.val, smem.idx);}
__device__ __forceinline__  lloc ldloc(volatile  lloc & smem){  return make_lloc (smem.val, smem.idx);}
__device__ __forceinline__ ulloc ldloc(volatile ulloc & smem){  return make_ulloc(smem.val, smem.idx);}


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
