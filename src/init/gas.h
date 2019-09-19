/**
 * @file gas.h
 *
 * @brief Header file for describing radial profile of spherical gaseous component(s)
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2019/09/19 (Thu)
 *
 * Copyright (C) 2019 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef GAS_H
#define GAS_H

#ifdef  ENABLE_GASEOUS_COMPONENT


#include "macro.h"

#include "../misc/structure.h"


/**
 * @def CONSIDER_SELF_GRAVITY_OF_GAS_COMPONENT
 *
 * @brief enable: activate iteration to include effects of self-gravity of the gas component
 *
 */
/* #define CONSIDER_SELF_GRAVITY_OF_GAS_COMPONENT */


/** macros to specify the density distribution model for gaseous components */
/** positive value indicates spherical component(s) */
#define ISOTHERMAL_GAS (10000)
/** negative value indicates disk component(s) */




#endif//ENABLE_GASEOUS_COMPONENT

#endif//GAS_H
