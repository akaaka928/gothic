/**
 * @file device.h
 *
 * @brief Header file for configuration of GPUs in GOTHIC
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
#ifndef DEVICE_H
#define DEVICE_H


/**
 * @def GPUS_PER_PROCESS
 *
 * @brief Number of GPUs per MPI process
 */
#   if  !defined(GPUS_PER_PROCESS) && !defined(_OPENMP)
#        define  GPUS_PER_PROCESS (1)
#endif//!defined(GPUS_PER_PROCESS) && !defined(_OPENMP)
#   if  !defined(_OPENMP)
#undef  GPUS_PER_PROCESS
#define GPUS_PER_PROCESS (1)
#endif//!defined(_OPENMP)


/**
 * @def SMEM_SIZE
 *
 * @brief Size of allocated shared memory on GPU (must be 0)
 */
#ifndef SMEM_SIZE
#define SMEM_SIZE (0)
#endif//SMEM_SIZE

/**
 * @def SMPREF
 *
 * @brief Shared memory preferred (1) or L1 cache preferred (0); default in CUDA is 1, global default in GOTHIC is 0
 */
#ifndef SMPREF
#define SMPREF (1)
#endif//SMPREF

/**
 * @def SMPREF_LET
 *
 * @brief Shared memory preferred (1) or L1 cache preferred (0)
 */
#ifndef SMPREF_LET
#define SMPREF_LET SMPREF
#endif//SMPREF_LET

/**
 * @def WIDEBANK
 *
 * @brief 1 set shared memory bank width to be 8 bytes while 0 is 4 bytes; default in CUDA is 0, global defalut in GOTHIC is 1
 */
#ifndef WIDEBANK
#define WIDEBANK (0)
#endif//WIDEBANK


/**
 * @def GLOBAL_MEMORY_SYSBUF
 *
 * @brief Size of remained global memory for buffers
 */
/* #define GLOBAL_MEMORY_SYSBUF (128 * 1048576) */
#define GLOBAL_MEMORY_SYSBUF (256 * 1048576)
/* /\* for nvprof *\/ */
/* #define GLOBAL_MEMORY_SYSBUF (512 * 1048576) */


#endif//DEVICE_H
