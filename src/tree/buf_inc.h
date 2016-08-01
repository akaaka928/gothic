/*************************************************************************\
 *                                                                       *
                  last updated on 2016/08/01(Mon) 16:58:35
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
  int bufSize;
} soaTreeWalkBuf;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//BUF_INC_H
//-------------------------------------------------------------------------
