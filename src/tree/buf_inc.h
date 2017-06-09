/**
 * @file buf_inc.h
 *
 * @brief Header file for booking buffer on global memory
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/06 (Mon)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef BUF_INC_H
#define BUF_INC_H


#define USE_SMID_TO_GET_BUFID
#define TRY_MODE_ABOUT_BUFFER


/**
 * @struct soaTreeWalkBuf
 *
 * @brief structure for buffer to walk tree
 */
typedef struct
{
  uint *freeLst, *buffer;
  int *fail;
#ifndef USE_SMID_TO_GET_BUFID
#ifndef TRY_MODE_ABOUT_BUFFER
  uint *freeNum;
  int *active;
#else///TRY_MODE_ABOUT_BUFFER
  int freeNum;
#endif//TRY_MODE_ABOUT_BUFFER
#endif//USE_SMID_TO_GET_BUFID
  int bufSize;
} soaTreeWalkBuf;


#endif//BUF_INC_H
