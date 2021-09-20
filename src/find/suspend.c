/**
 * @file suspend.c
 *
 * @brief Suspend simulation for a specified parameter set if the corresponding simulation takes too long
 *
 * @author Yohei MIKI (University of Tokyo)
 *
 * @date 2020/01/11 (Sat)
 *
 * Copyright (C) 2020 Yohei MIKI
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "macro.h"
#include "myutil.h"

#include "../anal/observation.h"


int main(int argc, char **argv)
{
  if( argc < 2 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 2);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 2 ){ */
  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv, "file", &file));

  suspend_model(file);

  return (0);
}
