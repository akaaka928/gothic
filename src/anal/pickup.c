/**
 * @file pickup.c
 *
 * @brief Pickup N-body particles for data analysis
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2019/09/13 (Fri)
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







void read_boundary(char *filename, int *Ngroup)
{
  /** read input */
  FILE *fp;
  fp = fopen(filename, "r");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);  }
  int checker = 1;

  /** read the number of input files */
  checker &= (1 == fscanf(fp, "%d", Ngroup));




}


static inline void check_boundary(real * restrict fmin, real * restrict fmax)
{
  if( *fmin > *fmax ){
    const real temp = *fmin;
    *fmin = *fmax;
    *fmax =  temp;
  }
}


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);


  /**# how to specify multiple groups? */


  /** initialization */
  if( argc < 7 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 7);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -target=<int>\n");/**< index of target snapshot */
    __FPRINTF__(stderr, "          -mode=<int>\n");/**< 0: (X, Y, Z), 1: (x, y, z), 2: (xi, eta, D) */
    __FPRINTF__(stderr, "          --xmin=<real, optioanl> -xmax=<real, optioanl>\n");
    __FPRINTF__(stderr, "          --ymin=<real, optioanl> -ymax=<real, optioanl>\n");
    __FPRINTF__(stderr, "          --zmin=<real, optioanl> -zmax=<real, optioanl>\n");
    __FPRINTF__(stderr, "          --vxmin=<real, optioanl> -vxmax=<real, optioanl>\n");
    __FPRINTF__(stderr, "          --vymin=<real, optioanl> -vymax=<real, optioanl>\n");
    __FPRINTF__(stderr, "          --vzmin=<real, optioanl> -vzmax=<real, optioanl>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 7 ){ */

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));

  /** read settings about the simulation */
  static ulong Ntot;
  static real eps, eta;
  static double ft, snapshotInterval, saveInterval;
  static int unit;
  readSettings(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file);

  int target;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "target", &target));
  int mode;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "mode", &mode));
  if( (mode != 0) && (mode != 1) && (mode != 2) ){
    __KILL__(stderr, "value of mode must be 0 (simulation coordinate), 1 (observed coordinate), or 2 (M31 standard coordinate) instead of your specified value (%d)\n", mode);
  }

  double tmp;
  extern const double length_astro2com;
  real xmin = REAL_MIN;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "xmin", &tmp)) == myUtilAvail )
    xmin = CAST_D2R(tmp * ((mode != 2) ? length_astro2com : 1.0));
  real xmax = REAL_MAX;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "xmax", &tmp)) == myUtilAvail )
    xmax = CAST_D2R(tmp * ((mode != 2) ? length_astro2com : 1.0));
  real ymin = REAL_MIN;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "ymin", &tmp)) == myUtilAvail )
    ymin = CAST_D2R(tmp * ((mode != 2) ? length_astro2com : 1.0));
  real ymax = REAL_MAX;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "ymax", &tmp)) == myUtilAvail )
    ymax = CAST_D2R(tmp * ((mode != 2) ? length_astro2com : 1.0));
  real zmin = REAL_MIN;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "zmin", &tmp)) == myUtilAvail )
    zmin = CAST_D2R(tmp * length_astro2com);
  real zmax = REAL_MAX;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "zmax", &tmp)) == myUtilAvail )
    zmax = CAST_D2R(tmp * length_astro2com);
  extern const double velocity_astro2com;
  real vxmin = REAL_MIN;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "vxmin", &tmp)) == myUtilAvail )
    vxmin = CAST_D2R(tmp * velocity_astro2com);
  real vxmax = REAL_MAX;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "vxmax", &tmp)) == myUtilAvail )
    vxmax = CAST_D2R(tmp * velocity_astro2com);
  real vymin = REAL_MIN;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "vymin", &tmp)) == myUtilAvail )
    vymin = CAST_D2R(tmp * velocity_astro2com);
  real vymax = REAL_MAX;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "vymax", &tmp)) == myUtilAvail )
    vymax = CAST_D2R(tmp * velocity_astro2com);
  real vzmin = REAL_MIN;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "vzmin", &tmp)) == myUtilAvail )
    vzmin = CAST_D2R(tmp * velocity_astro2com);
  real vzmax = REAL_MAX;
  if( optionalCmdArg(getCmdArgDbl(argc, (const char * const *)(void *)argv, "vzmax", &tmp)) == myUtilAvail )
    vzmax = CAST_D2R(tmp * velocity_astro2com);

  check_boundary(& xmin, & xmax);
  check_boundary(& ymin, & ymax);
  check_boundary(& zmin, & zmax);
  check_boundary(&vxmin, &vxmax);
  check_boundary(&vymin, &vymax);
  check_boundary(&vzmin, &vzmax);


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


  int numList = 0;
  ulong *list;
  list = (ulong *)malloc(sizeof(ulong) * Ntot);
  if( list == NULL ){    __KILL__(stderr, "ERROR: failure to allocate list\n");  }
  for(ulong ii = 0; ii < Ntot; ii++)
    list[ii] = Ntot;


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
    for(ulong ii = 0; ii < Ntot; ii++)
      if( (body[ii].x >= xmin) && (body[ii].x <= xmax) )
	if( (body[ii].y >= ymin) && (body[ii].y <= ymax) )
	  if( (body[ii].z >= zmin) && (body[ii].z <= zmax) )
	    if( (body[ii].vx >= vxmin) && (body[ii].vx <= vxmax) )
	      if( (body[ii].vy >= vymin) && (body[ii].vy <= vymax) )
		if( (body[ii].vz >= vzmin) && (body[ii].vz <= vzmax) )
		  {
		    list[numList] = body[ii].idx;
		    numList++;
		  }
  }


  chkMPIerr(MPI_Bcast(&numList, 1, MPI_INT, 0, mpi.comm));
  chkMPIerr(MPI_Bcast(list, numList, MPI_UNSIGNED_LONG, 0, mpi.comm));




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



  free(body);





}
