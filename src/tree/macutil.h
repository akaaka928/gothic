/**
 * @file macutil.h
 *
 * @brief Header file for multipole acceptance criterion (MAC)
 *
 * @author Yohei Miki (University of Tsukuba)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2017/03/01 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef MACUTIL_H
#define MACUTIL_H


#include "macro.h"


#   if  defined(GADGET_MAC) && defined(WS93_MAC)
#undef WS93_MAC
#endif//defined(GADGET_MAC) && defined(WS93_MAC)


/* list of functions appeared in ``macutil.c'' */
#ifdef  __CUDACC__
extern "C"
{
#endif//__CUDACC__
#ifndef WS93_MAC
  void setGlobalConstants_macutil_c(const real preset);
#endif//WS93_MAC

  void getBufSize(
#ifdef  WS93_MAC
		  const int pp, const real delta, const int Ni,
#endif//WS93_MAC
		  int *bufSize);
#ifdef  __CUDACC__
}
#endif//__CUDACC__


#endif//MAUTIL_H
