/**
 * @file cdflib.h
 *
 * @brief Header file for common settings for error analysis of tree code
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2017/10/26 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */
#ifndef CDFLIB_H
#define CDFLIB_H


#define NSUMMARY (28)


/* List of functions appeared in "cdflib.c" */
void allocPercentile(double **per);
void  freePercentile(double  *per);


#endif//CDFLIB_H
