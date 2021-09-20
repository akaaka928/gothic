/**
 * @file showOptConfig.c
 *
 * @brief Source code to show optimal configuration based on results of micro-benchmark
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

#include <stdio.h>
#include <stdlib.h>

#include "macro.h"
#include "myutil.h"
#include "name.h"


#define NADD (64)
typedef struct
{
  double timeWalkTree, timeMultipole, timeMakeTree, timeIntegrate;
  int Ttot_walk, Tsub_walk, Nwarp, Nloop, NeighborLevel;/**< walk tree */
  int Ttot_cmac, Tsub_cmac;/**< calc multipole */
  int Ttot_time;/**< time integration */
} perf;


#ifdef __ICC
/** Disable ICC's remark #161: unrecognized #pragma */
#     pragma warning (disable:161)
#endif//__ICC
int walkAscendingOrder(const void *a, const void *b);
int walkAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((perf *)a)->timeWalkTree > ((perf *)b)->timeWalkTree ){    return ( 1);  }
  else{    if( ((perf *)a)->timeWalkTree < ((perf *)b)->timeWalkTree ){    return (-1);  }
    else                                                                   return ( 0);  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

int cmacAscendingOrder(const void *a, const void *b);
int cmacAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((perf *)a)->timeMultipole > ((perf *)b)->timeMultipole ){    return ( 1);  }
  else{    if( ((perf *)a)->timeMultipole < ((perf *)b)->timeMultipole ){    return (-1);  }
    else                                                                     return ( 0);  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

int makeAscendingOrder(const void *a, const void *b);
int makeAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((perf *)a)->timeMakeTree > ((perf *)b)->timeMakeTree ){    return ( 1);  }
  else{    if( ((perf *)a)->timeMakeTree < ((perf *)b)->timeMakeTree ){    return (-1);  }
    else                                                                   return ( 0);  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

int timeAscendingOrder(const void *a, const void *b);
int timeAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((perf *)a)->timeIntegrate > ((perf *)b)->timeIntegrate ){    return ( 1);  }
  else{    if( ((perf *)a)->timeIntegrate < ((perf *)b)->timeIntegrate ){    return (-1);  }
    else                                                                     return ( 0);  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

#ifdef __ICC
/** Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC


void allocPerf(perf **data, const size_t num);
void readBenchmarkResults(int *num, int *rem, perf **dat, char *tag);


int main(int argc, char **argv)
{
  /** initialization */
  if( argc < 3 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 3);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -problem=<char *>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 3 ){ */

  char *file;     requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "file", &file));
  char *problem;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "problem", &problem));


  /** load measured performance */
  /** allocate a buffer */
  perf *dat;
  static int datNum, datRem;
  datNum = 0;
  datRem = NADD;
  allocPerf(&dat, datRem);

  static char tag[128];
  sprintf(tag, "%s.%s", file, problem);

  /** read measured data */
  readBenchmarkResults(&datNum, &datRem, &dat, tag);


  /** find the optimal configuration for each function */
#   if  defined(HUNT_WALK_PARAMETER) || defined(HUNT_NODE_PARAMETER) || defined(HUNT_MAKE_PARAMETER) || defined(HUNT_TIME_PARAMETER)
  int ii;
  FILE *fp;
  static char optfile[256];
  const double tolerance = 0.1;
#endif//defined(HUNT_WALK_PARAMETER) || defined(HUNT_NODE_PARAMETER) || defined(HUNT_MAKE_PARAMETER) || defined(HUNT_TIME_PARAMETER)

#ifdef  HUNT_WALK_PARAMETER
  /** elapsed time for tree traversal */
  qsort(dat, datNum, sizeof(perf), walkAscendingOrder);
  const double bestTimeWalkTree = dat[0].timeWalkTree;
  const double goodTimeWalkTree = bestTimeWalkTree * (1.0 + tolerance);
  sprintf(optfile, "%s/%s.%s.opt", BENCH_LOG_FOLDER, tag, "walk");
  fp = fopen(optfile, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", optfile);  }
  fprintf(fp, "#time(s)\tTtot\tTsub\tNwarp\tNloop\tLevel\n");
  for(ii = 0; ii < datNum; ii++){
    if( dat[ii].timeWalkTree > goodTimeWalkTree )      break;
    fprintf(fp, "%e\t%d\t%d\t%d\t%d\t%d\n", dat[ii].timeWalkTree, dat[ii].Ttot_walk, dat[ii].Tsub_walk, dat[ii].Nwarp, dat[ii].Nloop, dat[ii].NeighborLevel);
  }
  fprintf(fp, "#%d candidates are listed within %lf%% range\n", ii, tolerance * 100.0);

  fprintf(fp, "\n#optimal configuration for NWARP = 1:\n");
  fprintf(fp, "#time(s)\tTtot\tTsub\tNwarp\tNloop\tLevel\n");
  int counter = 0;
  for(ii = 0; ii < datNum; ii++){
    if( dat[ii].Nwarp == 1 ){
      fprintf(fp, "%e\t%d\t%d\t%d\t%d\t%d\n", dat[ii].timeWalkTree, dat[ii].Ttot_walk, dat[ii].Tsub_walk, dat[ii].Nwarp, dat[ii].Nloop, dat[ii].NeighborLevel);
      counter++;
    }
    if( counter > 10 )
      break;
  }

  fprintf(fp, "\n#optimal configuration for NWARP = 2:\n");
  fprintf(fp, "#time(s)\tTtot\tTsub\tNwarp\tNloop\tLevel\n");
  counter = 0;
  for(ii = 0; ii < datNum; ii++){
    if( dat[ii].Nwarp == 2 ){
      fprintf(fp, "%e\t%d\t%d\t%d\t%d\t%d\n", dat[ii].timeWalkTree, dat[ii].Ttot_walk, dat[ii].Tsub_walk, dat[ii].Nwarp, dat[ii].Nloop, dat[ii].NeighborLevel);
      counter++;
    }
    if( counter > 10 )
      break;
  }

  fprintf(fp, "\n#optimal configuration for NWARP = 4:\n");
  fprintf(fp, "#time(s)\tTtot\tTsub\tNwarp\tNloop\tLevel\n");
  counter = 0;
  for(ii = 0; ii < datNum; ii++){
    if( dat[ii].Nwarp == 4 ){
      fprintf(fp, "%e\t%d\t%d\t%d\t%d\t%d\n", dat[ii].timeWalkTree, dat[ii].Ttot_walk, dat[ii].Tsub_walk, dat[ii].Nwarp, dat[ii].Nloop, dat[ii].NeighborLevel);
      counter++;
    }
    if( counter > 10 )
      break;
  }

  fprintf(fp, "\n#optimal configuration for NWARP = 8:\n");
  fprintf(fp, "#time(s)\tTtot\tTsub\tNwarp\tNloop\tLevel\n");
  counter = 0;
  for(ii = 0; ii < datNum; ii++){
    if( dat[ii].Nwarp == 8 ){
      fprintf(fp, "%e\t%d\t%d\t%d\t%d\t%d\n", dat[ii].timeWalkTree, dat[ii].Ttot_walk, dat[ii].Tsub_walk, dat[ii].Nwarp, dat[ii].Nloop, dat[ii].NeighborLevel);
      counter++;
    }
    if( counter > 10 )
      break;
  }

  fprintf(fp, "\n#optimal configuration for NWARP = 16:\n");
  fprintf(fp, "#time(s)\tTtot\tTsub\tNwarp\tNloop\tLevel\n");
  counter = 0;
  for(ii = 0; ii < datNum; ii++){
    if( dat[ii].Nwarp == 16 ){
      fprintf(fp, "%e\t%d\t%d\t%d\t%d\t%d\n", dat[ii].timeWalkTree, dat[ii].Ttot_walk, dat[ii].Tsub_walk, dat[ii].Nwarp, dat[ii].Nloop, dat[ii].NeighborLevel);
      counter++;
    }
    if( counter > 10 )
      break;
  }

  fprintf(fp, "\n#optimal configuration for NWARP = 32:\n");
  fprintf(fp, "#time(s)\tTtot\tTsub\tNwarp\tNloop\tLevel\n");
  counter = 0;
  for(ii = 0; ii < datNum; ii++){
    if( dat[ii].Nwarp == 32 ){
      fprintf(fp, "%e\t%d\t%d\t%d\t%d\t%d\n", dat[ii].timeWalkTree, dat[ii].Ttot_walk, dat[ii].Tsub_walk, dat[ii].Nwarp, dat[ii].Nloop, dat[ii].NeighborLevel);
      counter++;
    }
    if( counter > 10 )
      break;
  }

  fclose(fp);
#endif//HUNT_WALK_PARAMETER

#ifdef  HUNT_NODE_PARAMETER
  /* elapsed time for calculate multipole */
  qsort(dat, datNum, sizeof(perf), cmacAscendingOrder);
  const double bestTimeCmacTree = dat[0].timeMultipole;
  const double goodTimeCmacTree = bestTimeCmacTree * (1.0 + tolerance);
  sprintf(optfile, "%s/%s.%s.opt", BENCH_LOG_FOLDER, tag, "cmac");
  fp = fopen(optfile, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", optfile);  }
  fprintf(fp, "#time(s)\tTtot\tTsub\n");
  for(ii = 0; ii < datNum; ii++){
    if( dat[ii].timeMultipole > goodTimeCmacTree )      break;
    fprintf(fp, "%e\t%d\t%d\n", dat[ii].timeMultipole, dat[ii].Ttot_cmac, dat[ii].Tsub_cmac);
  }
  fprintf(fp, "#%d candidates are listed within %lf%% range\n", ii, tolerance * 100.0);
  fclose(fp);
#endif//HUNT_NODE_PARAMETER

#ifdef  HUNT_MAKE_PARAMETER
  /** elapsed time for tree construction */
  qsort(dat, datNum, sizeof(perf), makeAscendingOrder);
  const double bestTimeMakeTree = dat[0].timeMakeTree;
  const double goodTimeMakeTree = bestTimeMakeTree * (1.0 + tolerance);
  sprintf(optfile, "%s/%s.%s.opt", BENCH_LOG_FOLDER, tag, "make");
  fp = fopen(optfile, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", optfile);  }
  fprintf(fp, "#time(s)\tTtot\tTsub\tNloop\tLevel\n");
  for(ii = 0; ii < datNum; ii++){
    if( dat[ii].timeMakeTree > goodTimeMakeTree )      break;
    fprintf(fp, "%e\t%d\t%d\t%d\t%d\n", dat[ii].timeMakeTree, dat[ii].Ttot, dat[ii].Tsub, dat[ii].Nloop, dat[ii].NeighborLevel);
  }
  fprintf(fp, "#%d candidates are listed within %lf%% range\n", ii, tolerance * 100.0);
  fclose(fp);
#endif//HUNT_MAKE_PARAMETER

#ifdef  HUNT_TIME_PARAMETER
  /** elapsed time for orbit integration */
  qsort(dat, datNum, sizeof(perf), timeAscendingOrder);
  const double bestTimeOrbitInt = dat[0].timeIntegrate;
  const double goodTimeOrbitInt = bestTimeOrbitInt * (1.0 + tolerance);
  sprintf(optfile, "%s/%s.%s.opt", BENCH_LOG_FOLDER, tag, "time");
  fp = fopen(optfile, "w");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", optfile);  }
  fprintf(fp, "#time(s)\tTtot\n");
  for(ii = 0; ii < datNum; ii++){
    if( dat[ii].timeIntegrate > goodTimeOrbitInt )      break;
    fprintf(fp, "%e\t%d\n", dat[ii].timeIntegrate, dat[ii].Ttot_time);
  }
  fprintf(fp, "#%d candidates are listed within %lf%% range\n", ii, tolerance * 100.0);
  fclose(fp);
#endif//HUNT_TIME_PARAMETER


  free(dat);

  return (0);
}


void allocPerf(perf **data, const size_t num)
{
  __NOTE__("%s\n", "start");

  *data = (perf *)malloc(num * sizeof(perf));
  if( *data == NULL ){    __KILL__(stderr, "ERROR: failure to allocate data\n");  }

  __NOTE__("%s\n", "end");
}
static inline void scalePerf(perf **data, const size_t num)
{
  __NOTE__("%s\n", "start");

  *data = (perf *)realloc(*data, num * sizeof(perf));
  if( *data == NULL ){    __KILL__(stderr, "ERROR: failure to increase data\n");  }

  __NOTE__("%s\n", "end");
}

void readBenchmarkResults(int *num, int *rem, perf **dat, char *tag)
{
  __NOTE__("%s\n", "start");

  FILE *fp;
  static char log[256];
  sprintf(log, "%s/%s.log", BENCH_LOG_FOLDER, tag);
  fp = fopen(log, "r");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", log);
  }

  perf tmp;
  while( EOF != fscanf(fp, "%le\t%le\t%le\t%le\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
		       &(tmp.timeWalkTree),
		       &(tmp.timeMultipole),
		       &(tmp.timeMakeTree),
		       &(tmp.timeIntegrate),
		       &(tmp.Ttot_walk), &(tmp.Tsub_walk), &(tmp.Nwarp), &(tmp.Nloop), &(tmp.NeighborLevel),
		       &(tmp.Ttot_cmac), &(tmp.Tsub_cmac),
		       &(tmp.Ttot_time)
		       )
	 ){

    if( *rem == 0 ){
      scalePerf(dat, NADD + (*num));
      *rem += NADD;
    }

    (*dat)[*num] = tmp;
    *num += 1;
    *rem -= 1;
  }
  fclose(fp);

  __NOTE__("%s\n", "end");
}
