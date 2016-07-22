/*************************************************************************\
 *                                                                       *
                  last updated on 2016/07/15(Fri) 11:02:12
 *                                                                       *
 *    Header File for configuration of accelerator devices (GPUs)        *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef DEVICE_H
#define DEVICE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#   if  !defined(GPUS_PER_PROCESS) && !defined(_OPENMP)
#        define  GPUS_PER_PROCESS (1)
#endif//!defined(GPUS_PER_PROCESS) && !defined(_OPENMP)
#   if  !defined(_OPENMP)
#undef  GPUS_PER_PROCESS
#define GPUS_PER_PROCESS (1)
#endif//!defined(_OPENMP)
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
/* #ifdef  __CUDACC__ */
//-------------------------------------------------------------------------
/* #ifdef  _OPENMP */
/* const int smemSize = 0; */
/* #endif//_OPENMP */
#ifndef SMEM_SIZE
#define SMEM_SIZE (0)
#endif//SMEM_SIZE
//-------------------------------------------------------------------------
#ifndef SMPREF
#define SMPREF (1)
#endif//SMPREF
//-------------------------------------------------------------------------
#ifndef SMPREF_LET
#define SMPREF_LET SMPREF
#endif//SMPREF_LET
//-------------------------------------------------------------------------
#ifndef WIDEBANK
#define WIDEBANK (1)
#endif//WIDEBANK
//-------------------------------------------------------------------------
/* #endif//__CUDACC__ */
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define GLOBAL_MEMORY_SYSBUF (128 * 1048576)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//DEVICE_H
//-------------------------------------------------------------------------
