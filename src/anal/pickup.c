/**
 * @file pickup.c
 *
 * @brief Pickup N-body particles for data analysis
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2019/09/26 (Thu)
 *
 * Copyright (C) 2019 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#ifdef  USE_HDF5_FORMAT
#ifdef  HDF5_FOR_ZINDAIJI
#include <unistd.h>
#include <sys/stat.h>
#endif//HDF5_FOR_ZINDAIJI
#include <hdf5.h>
#include "hdf5lib.h"
/* The maximum number of elements in a chunk is 2^32-1 which is equal to 4,294,967,295 */
/* The maximum size for any chunk is 4GB */
#define MAXIMUM_CHUNK_SIZE      ((hsize_t)1 << 31)
#define MAXIMUM_CHUNK_SIZE_4BIT ((hsize_t)1 << 30)
#define MAXIMUM_CHUNK_SIZE_8BIT ((hsize_t)1 << 29)
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "myutil.h"
#include "name.h"
#include "constants.h"
#include "mpilib.h"
#include "rotate.h"

#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../file/io.h"
#include "../anal/m31coord.h"

extern const double length_astro2com;
extern const double velocity_astro2com;


#ifdef  __ICC
/* Disable ICC's remark #161: unrecognized #pragma */
#     pragma warning (disable:161)
#endif//__ICC

int idxAscendingOrder(const void *a, const void *b);
int idxAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((const nbody_aos *)a)->idx > ((const nbody_aos *)b)->idx ){    return ( 1);  }
  else{    if( ((const nbody_aos *)a)->idx < ((const nbody_aos *)b)->idx ){    return (-1);  }
    else{                                                                      return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}
#ifdef  __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC


/**
 * @struct domain
 *
 * @brief structure for each domain
 */
typedef struct
{
  real xmin, ymin, zmin;
  real xmax, ymax, zmax;
  real vxmin, vymin, vzmin;
  real vxmax, vymax, vzmax;
  int mode;
} domain;


/**
 * @fn writeCfgFormat
 *
 * @brief Print the expected format.
 *
 * @param (filename) the specified file name
 */
static inline void writeCfgFormat(char *filename)
{
  __NOTE__("%s\n", "start");

  fprintf(stderr, "ERROR: information written in \"%s\" does not match with the expected format\n", filename);
  fprintf(stderr, "Expected format is below:\n");
  fprintf(stderr, "\tmode<int>: index of selection rules: 0: (X, Y, Z), 1: (x, y, z), 2: (xi, eta, D)\n");
  fprintf(stderr, "\tNgrp<int>: number of selected groups\n");
  fprintf(stderr, "\tcfg0<char *>: filename of configuration for the 1st object\n");
  fprintf(stderr, "\tcfg1<char *>: filename of configuration for the 2nd object\n");
  fprintf(stderr, "\tcfg2<char *>: filename of configuration for the 3rd object\n");

  __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);

  __NOTE__("%s\n", "end");
}


/**
 * @fn writeBoundaryFormat
 *
 * @brief Print the expected format.
 *
 * @param (filename) the specified file name
 */
static inline void writeBoundaryFormat(char *filename)
{
  __NOTE__("%s\n", "start");

  fprintf(stderr, "ERROR: information written in \"%s\" does not match with the expected format\n", filename);
  fprintf(stderr, "Expected format is below:\n");
  fprintf(stderr, "\tmode<int>: index of selection rules: 0: (X, Y, Z), 1: (x, y, z), 2: (xi, eta, D)\n");
  fprintf(stderr, "\tuse_xmin<int: 0 or 1>\txmin<double>: specify (1) or not (0) the input value as xmin\n");
  fprintf(stderr, "\tuse_xmax<int: 0 or 1>\txmax<double>: specify (1) or not (0) the input value as xmax\n");
  fprintf(stderr, "\tuse_ymin<int: 0 or 1>\tymin<double>: specify (1) or not (0) the input value as ymin\n");
  fprintf(stderr, "\tuse_ymax<int: 0 or 1>\tymax<double>: specify (1) or not (0) the input value as ymax\n");
  fprintf(stderr, "\tuse_zmin<int: 0 or 1>\tzmin<double>: specify (1) or not (0) the input value as zmin\n");
  fprintf(stderr, "\tuse_zmax<int: 0 or 1>\tzmax<double>: specify (1) or not (0) the input value as zmax\n");
  fprintf(stderr, "\tuse_vxmin<int: 0 or 1>\tvxmin<double>: specify (1) or not (0) the input value as vxmin\n");
  fprintf(stderr, "\tuse_vxmax<int: 0 or 1>\tvxmax<double>: specify (1) or not (0) the input value as vxmax\n");
  fprintf(stderr, "\tuse_vymin<int: 0 or 1>\tvymin<double>: specify (1) or not (0) the input value as vymin\n");
  fprintf(stderr, "\tuse_vymax<int: 0 or 1>\tvymax<double>: specify (1) or not (0) the input value as vymax\n");
  fprintf(stderr, "\tuse_vzmin<int: 0 or 1>\tvzmin<double>: specify (1) or not (0) the input value as vzmin\n");
  fprintf(stderr, "\tuse_vzmax<int: 0 or 1>\tvzmax<double>: specify (1) or not (0) the input value as vzmax\n");

  __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);

  __NOTE__("%s\n", "end");
}


static inline void check_boundary(real * restrict fmin, real * restrict fmax)
{
  if( *fmin > *fmax ){
    const real temp = *fmin;
    *fmin = *fmax;
    *fmax =  temp;
  }
}


static void read_boundary(char *filename, int *mode, int *Ngrp, domain **box)
{
  /** read input */
  FILE *fp;
  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  int checker = 1;

  /** read the number of input files */
  checker &= (1 == fscanf(fp, "%d", mode));
  checker &= (1 == fscanf(fp, "%d", Ngrp));
  *box = (domain *)malloc(sizeof(domain) * (*Ngrp));
  if( *box == NULL ){
    __KILL__(stderr, "ERROR: failure to allocate box (cfg file is \"%s\", Ngrp = %d)\n", filename, *Ngrp);
  }

  for(int ii = 0; ii < *Ngrp; ii++){
    /** initialize dataset */
    (*box)[ii].mode = -2;
    (*box)[ii]. xmin = (*box)[ii]. ymin = (*box)[ii]. zmin = REAL_MIN;
    (*box)[ii]. xmax = (*box)[ii]. ymax = (*box)[ii]. zmax = REAL_MAX;
    (*box)[ii].vxmin = (*box)[ii].vymin = (*box)[ii].vzmin = REAL_MIN;
    (*box)[ii].vxmax = (*box)[ii].vymax = (*box)[ii].vzmax = REAL_MAX;

    /** read target filename */
    char boxfile[256];
    checker &= (1 == fscanf(fp, "%s", boxfile));
    FILE *fpbox;
    fpbox = fopen(boxfile, "r");
    if( fpbox == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", boxfile);    }

    int box_checker = 1;
    box_checker &= (1 == fscanf(fpbox, "%d", &(*box)[ii].mode));
    if( (*box)[ii].mode != *mode ){
      __KILL__(stderr, "ERROR: specified value of mode (%d) in \"%s\" must be identical to that in \"%s\" (%d).\n", (*box)[ii].mode, boxfile, filename, *mode);
    }

    double range;
    int specify = 0;
    /** mode: 0: (X, Y, Z), 1: (x, y, z), 2: (xi, eta, D) */
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].xmin = CAST_D2R(range * ((*mode != 2) ? length_astro2com : 1.0));    specify = 0;
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].xmax = CAST_D2R(range * ((*mode != 2) ? length_astro2com : 1.0));    specify = 0;
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].ymin = CAST_D2R(range * ((*mode != 2) ? length_astro2com : 1.0));    specify = 0;
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].ymax = CAST_D2R(range * ((*mode != 2) ? length_astro2com : 1.0));    specify = 0;
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].zmin = CAST_D2R(range * length_astro2com);    specify = 0;
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].zmax = CAST_D2R(range * length_astro2com);    specify = 0;

    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].vxmin = CAST_D2R(range * velocity_astro2com);    specify = 0;
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].vxmax = CAST_D2R(range * velocity_astro2com);    specify = 0;
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].vymin = CAST_D2R(range * velocity_astro2com);    specify = 0;
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].vymax = CAST_D2R(range * velocity_astro2com);    specify = 0;
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].vzmin = CAST_D2R(range * velocity_astro2com);    specify = 0;
    box_checker &= (2 == fscanf(fpbox, "%d\t%le", &specify, &range));
    if( specify )      (*box)[ii].vzmax = CAST_D2R(range * velocity_astro2com);    specify = 0;

    check_boundary(&(*box)[ii].xmin, &(*box)[ii].xmax);
    check_boundary(&(*box)[ii].ymin, &(*box)[ii].ymax);
    check_boundary(&(*box)[ii].zmin, &(*box)[ii].zmax);
    check_boundary(&(*box)[ii].vxmin, &(*box)[ii].vxmax);
    check_boundary(&(*box)[ii].vymin, &(*box)[ii].vymax);
    check_boundary(&(*box)[ii].vzmin, &(*box)[ii].vzmax);

    /** close the target file */
    fclose(fpbox);

    if( !box_checker )
      writeBoundaryFormat(boxfile);
  }

  fclose(fp);
  if( !checker )
    writeCfgFormat(filename);

#if 0
  fprintf(stdout, "file = %s, mode = %d, Ngrp = %d\n", filename, *mode, *Ngrp);

  for(int ii = 0; ii < *Ngrp; ii++){
    fprintf(stdout, "%d-th group:\n", ii);
    fprintf(stdout, "\txmin = %e, ymin = %e, zmin = %e\n", (*box)[ii].xmin, (*box)[ii].ymin, (*box)[ii].zmin);
    fprintf(stdout, "\txmax = %e, ymax = %e, zmax = %e\n", (*box)[ii].xmax, (*box)[ii].ymax, (*box)[ii].zmax);
    fprintf(stdout, "\tvxmin = %e, vymin = %e, vzmin = %e\n", (*box)[ii].vxmin, (*box)[ii].vymin, (*box)[ii].vzmin);
    fprintf(stdout, "\tvxmax = %e, vymax = %e, vzmax = %e\n", (*box)[ii].vxmax, (*box)[ii].vymax, (*box)[ii].vzmax);
  }

  exitMPI();
  exit(0);
#endif
}


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);


  /** initialization */
  if( argc < 7 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 7);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -target=<int>\n");/**< index of target snapshot */
    __FPRINTF__(stderr, "          -list=<char *>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 7 ){ */

  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "file", &file));
  char *fcfg;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "list", &fcfg));

  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));
  int target;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "target", &target));

  /** read settings about the simulation */
  static ulong Ntot;
  static real eps, eta;
  static double ft, snapshotInterval, saveInterval;
  static int unit;
  readSettings(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file);

  /** read selection rules */
  int mode = -1;
  int Ngrp = 0;
  domain *box;
  read_boundary(fcfg, &mode, &Ngrp, &box);

  if( (mode != 0) && (mode != 1) && (mode != 2) ){
    __KILL__(stderr, "value of mode must be 0 (simulation coordinate), 1 (observed coordinate), or 2 (M31 standard coordinate) instead of your specified value (%d)\n", mode);
  }


#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
  nbody_hdf5 hdf5;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  real *hdf5_acc_ext, *hdf5_pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  ulong *hdf5_idx;
  allocSnapshotArray
    (&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     &hdf5_acc_ext, &hdf5_pot_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     (int)Ntot, &hdf5);
#else///USE_HDF5_FORMAT
  iparticle ibody;
  ulong *idx;
  position *pos;
  acceleration *acc;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  acceleration *ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  velocity *vel;
  ibody_time *ti;
#else///BLOCK_TIME_STEP
  real *vx, *vy, *vz;
#endif//BLOCK_TIME_STEP
  allocParticleData
    ((int)Ntot, &ibody, &idx, &pos, &acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     &ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
     &vel, &ti
#else///BLOCK_TIME_STEP
     &vx, &vy, &vz
#endif//BLOCK_TIME_STEP
     );
#endif//USE_HDF5_FORMAT


  /** read number of components */
  int kind = 0;
  int *bodyHead, *bodyNum;
  FILE *fp;
  char filename[128];
  sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
  fp = fopen(filename, "r");
  if( fp == NULL ){
    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
  }
  int unit_tmp;
  bool checker = true;
  checker &= (1 == fscanf(fp, "%d", &unit_tmp));
  checker &= (1 == fscanf(fp, "%d\t%*d", &kind));
  bodyHead = (int *)malloc(sizeof(int) * kind);  if( bodyHead == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate bodyHead");  }
  bodyNum  = (int *)malloc(sizeof(int) * kind);  if( bodyNum  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate bodyNum");  }
#ifdef  HDF5_FOR_ZINDAIJI
  int *bodyType;
  bodyType = (int *)malloc(sizeof(int) * kind);  if( bodyType == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate bodyType");  }
#endif//HDF5_FOR_ZINDAIJI
  for(int ii = 0; ii < kind; ii++)
    checker &= (1 == fscanf(fp, "%d", &bodyNum[ii]));
#ifdef  HDF5_FOR_ZINDAIJI
  for(int ii = 0; ii < kind; ii++)
    checker &= (1 == fscanf(fp, "%d", &bodyType[ii]));
#endif//HDF5_FOR_ZINDAIJI
  fclose(fp);
  if( !checker ){
    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);
  }
  bodyHead[0] = 0;
  for(int ii = 1; ii < kind; ii++)
    bodyHead[ii] = bodyHead[ii - 1] + bodyNum[ii - 1];
#ifdef  HDF5_FOR_ZINDAIJI
  for(int ii = 0; ii < kind; ii++)
    bodyType[ii] &= 3;
#endif//HDF5_FOR_ZINDAIJI


  nbody_aos *body;
  body = (nbody_aos *)malloc(sizeof(nbody_aos) * Ntot);
  if( body == NULL ){    __KILL__(stderr, "ERROR: failure to allocate body\n");  }


  int numList = 0;
  ulong *list;
  list = (ulong *)malloc(sizeof(ulong) * Ntot);
  if( list == NULL ){    __KILL__(stderr, "ERROR: failure to allocate list\n");  }
  for(ulong ii = 0; ii < Ntot; ii++)
    list[ii] = Ntot;
  int *grpHead, *grpNum;
  grpHead = (int *)malloc(sizeof(int) * kind * Ngrp);  if( grpHead == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate grpHead");  }
  grpNum  = (int *)malloc(sizeof(int) * kind * Ngrp);  if( grpNum  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate grpNum");  }


  /** read target snapshot and pickup particles */
  if( mpi.rank == 0 ){
    /** read snapshot */
    double time;
    ulong steps;
    int unit_read;
#ifdef  USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file, (uint)target, &hdf5, hdf5type);
    for(int ii = 0; ii < (int)Ntot; ii++){
      body[ii]. x  = hdf5.pos[ii * 3];      body[ii]. y = hdf5.pos[ii * 3 + 1];      body[ii].z   = hdf5.pos[ii * 3 + 2];
      body[ii].vx  = hdf5.vel[ii * 3];      body[ii].vy = hdf5.vel[ii * 3 + 1];      body[ii].vz  = hdf5.vel[ii * 3 + 2];
      body[ii].ax  = hdf5.acc[ii * 3];      body[ii].ay = hdf5.acc[ii * 3 + 1];      body[ii].az  = hdf5.acc[ii * 3 + 2];
      body[ii].idx = hdf5.idx[ii    ];      body[ii]. m = hdf5.  m[ii        ];      body[ii].pot = hdf5.pot[ii        ];
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      body[ii].ax_ext = hdf5.acc_ext[ii * 3];      body[ii].ay_ext = hdf5.acc_ext[ii * 3 + 1];      body[ii].az_ext = hdf5.acc_ext[ii * 3 + 2];
      body[ii].pot_ext = hdf5.pot_ext[ii];
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#else///USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file, (uint)target, ibody);
    for(int ii = 0; ii < (int)Ntot; ii++){
      body[ii]. x  = ibody.pos[ii].x;      body[ii]. y = ibody.pos[ii].y;      body[ii]. z  = ibody.pos[ii].z;
      body[ii].vx  = ibody.vel[ii].x;      body[ii].vy = ibody.vel[ii].y;      body[ii].vz  = ibody.vel[ii].z;
      body[ii].ax  = ibody.acc[ii].x;      body[ii].ay = ibody.acc[ii].y;      body[ii].az  = ibody.acc[ii].z;
      body[ii].idx = ibody.idx[ii]  ;      body[ii]. m = ibody.pos[ii].m;      body[ii].pot = ibody.acc[ii].pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      body[ii].ax_ext = ibody.acc_ext[ii].x;      body[ii].ay_ext = ibody.acc_ext[ii].y;      body[ii].az_ext = ibody.acc_ext[ii].z;
      body[ii].pot_ext = ibody.acc_ext[ii].pot;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#endif//USE_HDF5_FORMAT
    if( unit_read != unit ){
      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
    }/* if( unit_read != unit ){ */


    /** sort particle data by index */
    qsort(body, Ntot, sizeof(nbody_aos), idxAscendingOrder);


    /** coordinate rotation when required */
    if( mode != 0 ){
      /** set coordinate transformation from M31's disk coordinate to observerd frame */
      static real rot[3][3], inv[3][3];/**< rot: M31 disk frame (X, Y, Z) ==>> observed frame (x, y, z) */
      setRotationMatrix(inv, rot);

      if( mode == 1 ){
	static real ini[3], fin[3];
	for(ulong ii = 0; ii < Ntot; ii++){
	  /** coordinate rotation of particle position */
	  ini[0] = body[ii].x;
	  ini[1] = body[ii].y;
	  ini[2] = body[ii].z;
	  rotateVector(ini, rot, fin);
	  body[ii].x = fin[0];
	  body[ii].y = fin[1];
	  body[ii].z = fin[2] + (zm31 * length_astro2com);

	  /** coordinate rotation of particle velocity */
	  ini[0] = body[ii].vx;
	  ini[1] = body[ii].vy;
	  ini[2] = body[ii].vz;
	  rotateVector(ini, rot, fin);
	  body[ii].vx = fin[0] + (vm31x * velocity_astro2com);
	  body[ii].vy = fin[1] + (vm31y * velocity_astro2com);
	  body[ii].vz = fin[2] + (vm31z * velocity_astro2com);
	}
      }
      else if( mode == 2 ){
	real *  xi;  xi   = (real *)malloc(sizeof(real) * Ntot);  if(   xi == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate xi.");  }
	real * eta;  eta  = (real *)malloc(sizeof(real) * Ntot);  if(  eta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate eta.");  }
	real *dist;  dist = (real *)malloc(sizeof(real) * Ntot);  if( dist == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate dist.");  }
	real * vxi;  vxi  = (real *)malloc(sizeof(real) * Ntot);  if( vxi  == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vxi.");  }
	real *veta;  veta = (real *)malloc(sizeof(real) * Ntot);  if( veta == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate veta.");  }
	real *vlos;  vlos = (real *)malloc(sizeof(real) * Ntot);  if( vlos == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vlos.");  }
	standard_coordinate(Ntot, body, rot, xi, eta, dist, vxi, veta, vlos);

	for(ulong ii = 0; ii < Ntot; ii++){
	  body[ii].x = xi[ii];
	  body[ii].y = eta[ii];
	  body[ii].z = dist[ii];
	  body[ii].vx = vxi[ii];
	  body[ii].vy = veta[ii];
	  body[ii].vz = vlos[ii];
	}

	free(xi);
	free(eta);
	free(dist);
	free(vxi);
	free(veta);
	free(vlos);
      }
    }


    /** pick up N-body particles contained in the specified domain */
    int grpIdx = 0;
    for(int ii = 0; ii < Ngrp; ii++){
      const real xmin = box[ii].xmin;
      const real xmax = box[ii].xmax;
      const real ymin = box[ii].ymin;
      const real ymax = box[ii].ymax;
      const real zmin = box[ii].zmin;
      const real zmax = box[ii].zmax;
      const real vxmin = box[ii].vxmin;
      const real vxmax = box[ii].vxmax;
      const real vymin = box[ii].vymin;
      const real vymax = box[ii].vymax;
      const real vzmin = box[ii].vzmin;
      const real vzmax = box[ii].vzmax;
      for(int kk = 0; kk < kind; kk++){
	grpHead[grpIdx] = numList;
	const int head = bodyHead[kk];
	const int tail = head + bodyNum[kk];
	for(int ll = head; ll < tail; ll++)
	  if( (body[ll].x >= xmin) && (body[ll].x <= xmax) )
	    if( (body[ll].y >= ymin) && (body[ll].y <= ymax) )
	      if( (body[ll].z >= zmin) && (body[ll].z <= zmax) )
		if( (body[ll].vx >= vxmin) && (body[ll].vx <= vxmax) )
		  if( (body[ll].vy >= vymin) && (body[ll].vy <= vymax) )
		    if( (body[ll].vz >= vzmin) && (body[ll].vz <= vzmax) )
		      {
			list[numList] = body[ll].idx;
			numList++;
		      }
	grpNum[grpIdx] = numList - grpHead[grpIdx];
	grpIdx++;
      }
    }

  }
  chkMPIerr(MPI_Bcast(grpHead, Ngrp * kind, MPI_INT, 0, mpi.comm));
  chkMPIerr(MPI_Bcast(grpNum , Ngrp * kind, MPI_INT, 0, mpi.comm));
  chkMPIerr(MPI_Bcast(&numList, 1, MPI_INT, 0, mpi.comm));
  chkMPIerr(MPI_Bcast(list, numList, MPI_UNSIGNED_LONG, 0, mpi.comm));

  char renamed[128];
  sprintf(renamed, "%s-chosen", file);

  if( mpi.rank == 0 ){
    const int skind = kind;
    char summary[256];
    sprintf(summary, "%s/%s.summary.txt", DOCUMENTFOLDER, renamed);
    fp = fopen(summary, "w");
    if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", summary);  }
    fprintf(fp, "%d\n", unit);
    fprintf(fp, "%d\t%d\n", Ngrp * kind, Ngrp * skind);
    for(int ii = 0; ii < Ngrp * kind; ii++)
      fprintf(fp, "%d\n", grpNum[ii]);
#ifdef  HDF5_FOR_ZINDAIJI
    for(int ii = 0; ii < Ngrp * kind; ii++)
      fprintf(fp, "%d\n", ii & 3);
#endif//HDF5_FOR_ZINDAIJI
    fclose(fp);

    ulong Nchosen = 0;
    for(int ii = 0; ii < Ngrp * kind; ii++)
      Nchosen += grpNum[ii];
    writeSettings(unit, Nchosen, eps, eta, ft, snapshotInterval, saveInterval, renamed);

    int last;
    readConfigFile(&last, file);
    updateConfigFile(last, renamed);
  }


  nbody_aos *chosen_body;
  chosen_body = (nbody_aos *)malloc(sizeof(nbody_aos) * numList);
  if( chosen_body == NULL ){    __KILL__(stderr, "ERROR: failure to allocate chosen_body\n");  }


  for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){
    /** read snapshot */
    double time;
    ulong steps;
    int unit_read;
#ifdef  USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file, (uint)filenum, &hdf5, hdf5type);
    for(int ii = 0; ii < (int)Ntot; ii++){
      body[ii]. x  = hdf5.pos[ii * 3];      body[ii]. y = hdf5.pos[ii * 3 + 1];      body[ii].z   = hdf5.pos[ii * 3 + 2];
      body[ii].vx  = hdf5.vel[ii * 3];      body[ii].vy = hdf5.vel[ii * 3 + 1];      body[ii].vz  = hdf5.vel[ii * 3 + 2];
      body[ii].ax  = hdf5.acc[ii * 3];      body[ii].ay = hdf5.acc[ii * 3 + 1];      body[ii].az  = hdf5.acc[ii * 3 + 2];
      body[ii].idx = hdf5.idx[ii    ];      body[ii]. m = hdf5.  m[ii        ];      body[ii].pot = hdf5.pot[ii        ];
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      body[ii].ax_ext = hdf5.acc_ext[ii * 3];      body[ii].ay_ext = hdf5.acc_ext[ii * 3 + 1];      body[ii].az_ext = hdf5.acc_ext[ii * 3 + 2];
      body[ii].pot_ext = hdf5.pot_ext[ii];
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#else///USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file, (uint)filenum, ibody);
    for(int ii = 0; ii < (int)Ntot; ii++){
      body[ii]. x  = ibody.pos[ii].x;      body[ii]. y = ibody.pos[ii].y;      body[ii]. z  = ibody.pos[ii].z;
      body[ii].vx  = ibody.vel[ii].x;      body[ii].vy = ibody.vel[ii].y;      body[ii].vz  = ibody.vel[ii].z;
      body[ii].ax  = ibody.acc[ii].x;      body[ii].ay = ibody.acc[ii].y;      body[ii].az  = ibody.acc[ii].z;
      body[ii].idx = ibody.idx[ii]  ;      body[ii]. m = ibody.pos[ii].m;      body[ii].pot = ibody.acc[ii].pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      body[ii].ax_ext = ibody.acc_ext[ii].x;      body[ii].ay_ext = ibody.acc_ext[ii].y;      body[ii].az_ext = ibody.acc_ext[ii].z;
      body[ii].pot_ext = ibody.acc_ext[ii].pot;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#endif//USE_HDF5_FORMAT
    if( unit_read != unit ){
      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
    }/* if( unit_read != unit ){ */

    /** sort particle data by index */
    qsort(body, Ntot, sizeof(nbody_aos), idxAscendingOrder);


    /** copy selected particles */
    for(int ii = 0; ii < numList; ii++){
      chosen_body[ii] = body[list[ii]];
      chosen_body[ii].idx = (ulong)ii;
    }

    /** output renamed snapshot */
    for(int ii = 0; ii < numList; ii++){
#ifdef  USE_HDF5_FORMAT
      hdf5.pos[ii * 3] = chosen_body[ii]. x;      hdf5.pos[ii * 3 + 1] = chosen_body[ii]. y;      hdf5.pos[ii * 3 + 2] = chosen_body[ii]. z;
      hdf5.vel[ii * 3] = chosen_body[ii].vx;      hdf5.vel[ii * 3 + 1] = chosen_body[ii].vy;      hdf5.vel[ii * 3 + 2] = chosen_body[ii].vz;
      hdf5.acc[ii * 3] = chosen_body[ii].ax;      hdf5.acc[ii * 3 + 1] = chosen_body[ii].ay;      hdf5.acc[ii * 3 + 2] = chosen_body[ii].az;
      hdf5.idx[ii] = chosen_body[ii].idx;      hdf5.m[ii] = chosen_body[ii].m;      hdf5.pot[ii] = chosen_body[ii].pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      hdf5.acc_ext[ii * 3] = chosen_body[ii].ax_ext;      hdf5.acc_ext[ii * 3 + 1] = chosen_body[ii].ay_ext;      hdf5.acc_ext[ii * 3 + 2] = chosen_body[ii].az_ext;
      hdf5.pot_ext[ii] = chosen_body[ii].pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#else///USE_HDF5_FORMAT
      ibody.pos[ii].x = chosen_body[ii]. x;      ibody.pos[ii].y = chosen_body[ii]. y;      ibody.pos[ii].z = chosen_body[ii]. z;
      ibody.vel[ii].x = chosen_body[ii].vx;      ibody.vel[ii].y = chosen_body[ii].vy;      ibody.vel[ii].z = chosen_body[ii].vz;
      ibody.acc[ii].x = chosen_body[ii].ax;      ibody.acc[ii].y = chosen_body[ii].ay;      ibody.acc[ii].z = chosen_body[ii].az;
      ibody.idx[ii] = chosen_body[ii].idx;      ibody.pos[ii].m = chosen_body[ii].m;      ibody.acc[ii].pot = chosen_body[ii].pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
      ibody.acc_ext[ii].x = chosen_body[ii].ax_ext;      ibody.acc_ext[ii].y = chosen_body[ii].ay_ext;      ibody.acc_ext[ii].z = chosen_body[ii].az_ext;
      ibody.acc_ext[ii].pot = chosen_body[ii].pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#endif//USE_HDF5_FORMAT
    }
#   if  defined(USE_HDF5_FORMAT) && defined(MONITOR_ENERGY_ERROR)
    energyError relEneErr;
#endif//defined(USE_HDF5_FORMAT) && defined(MONITOR_ENERGY_ERROR)
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
    gpu_clock deviceMonitors;
    const int monitor_step = 0;
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  REPORT_COMPUTE_RATE
    const double speed = 0.0;
    const double speed_run = 0.0;
    const double complete = 0.0;
    const double guess = 0.0;
    const double brent_avg = 0.0;
    const double rebuild_interval = 0.0;
#endif//REPORT_COMPUTE_RATE
    writeSnapshot(unit_read, time, steps, numList, renamed, (uint)filenum
#ifdef  USE_HDF5_FORMAT
		  , &hdf5, hdf5type
#ifdef  MONITOR_ENERGY_ERROR
		  , &relEneErr, false
#endif//MONITOR_ENERGY_ERROR
#else///USE_HDF5_FORMAT
		  , ibody
#endif//USE_HDF5_FORMAT
#ifdef  REPORT_GPU_CLOCK_FREQUENCY
		  , &deviceMonitors, monitor_step
#endif//REPORT_GPU_CLOCK_FREQUENCY
#ifdef  REPORT_COMPUTE_RATE
		  , speed, speed_run, complete, guess, brent_avg, rebuild_interval
#endif//REPORT_COMPUTE_RATE
		  );
  }/* for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){ */


#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
  freeSnapshotArray
    (hdf5_pos, hdf5_vel, hdf5_acc, hdf5_m, hdf5_pot, hdf5_idx
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , hdf5_acc_ext, hdf5_pot_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
#else///USE_HDF5_FORMAT
  freeParticleData
    (idx, pos, acc,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
     vel, ti
#else///BLOCK_TIME_STEP
     vx, vy, vz
#endif//BLOCK_TIME_STEP
     );
#endif//USE_HDF5_FORMAT
  free(body);
  free(bodyHead);
  free(bodyNum);


  free(list);



  free(grpHead);
  free(grpNum);
  free(box);

  free(chosen_body);



  exitMPI();

  return (0);
}
