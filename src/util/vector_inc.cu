/**
 * @file vector_inc.cu
 *
 * @brief Source code for manipulating vectors on GPU
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
#ifndef VECTOR_INC_CU
#define VECTOR_INC_CU

#include <stdio.h>
#include <stdlib.h>
#include <helper_cuda.h>


template <typename Type> __device__ __forceinline__ void stvec(Type * smem, const int tidx, const Type val){  smem[tidx] = val;}

__device__ __forceinline__    int2 ldvec(volatile    int2 & smem){  return make_int2   (smem.x, smem.y);}
__device__ __forceinline__  float2 ldvec(volatile  float2 & smem){  return make_float2 (smem.x, smem.y);}
__device__ __forceinline__ double2 ldvec(volatile double2 & smem){  return make_double2(smem.x, smem.y);}

__device__ __forceinline__    int3 ldvec(volatile    int3 & smem){  return make_int3   (smem.x, smem.y, smem.z);}
__device__ __forceinline__  float3 ldvec(volatile  float3 & smem){  return make_float3 (smem.x, smem.y, smem.z);}
__device__ __forceinline__ double3 ldvec(volatile double3 & smem){  return make_double3(smem.x, smem.y, smem.z);}

__device__ __forceinline__    int4 ldvec(volatile    int4 & smem){  return make_int4   (smem.x, smem.y, smem.z, smem.w);}
__device__ __forceinline__  float4 ldvec(volatile  float4 & smem){  return make_float4 (smem.x, smem.y, smem.z, smem.w);}
__device__ __forceinline__ double4 ldvec(volatile double4 & smem){  return make_double4(smem.x, smem.y, smem.z, smem.w);}
#endif//COMPARISON_INC_CU
