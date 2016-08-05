/*************************************************************************\
 *                                                                       *
                  last updated on 2016/08/05(Fri) 14:24:11
 *                                                                       *
 *    Header File for tree traversal based on octree structure           *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef BUF_INC_H
#define BUF_INC_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define USE_SMID_TO_GET_BUFID
#define TRY_MODE_ABOUT_BUFFER
//-------------------------------------------------------------------------
typedef struct
{
  uint *freeLst, *buffer;
  int *fail;
#ifndef USE_SMID_TO_GET_BUFID
#       ifndef TRY_MODE_ABOUT_BUFFER
  uint *freeNum;
  int *active;
#       else///TRY_MODE_ABOUT_BUFFER
  int freeNum;
#       endif//TRY_MODE_ABOUT_BUFFER
#endif//USE_SMID_TO_GET_BUFID
  /* size_t bufSize; */
  int bufSize;
#ifdef  ALLOCATE_LETBUFFER
  uint *bufLET;
  /* size_t bufLETsize; */
  int bufLETsize;
#endif//ALLOCATE_LETBUFFER
} soaTreeWalkBuf;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//BUF_INC_H
//-------------------------------------------------------------------------
