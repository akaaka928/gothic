/**
 * @file potential_dev.cu
 *
 * @brief Source code for calculating gravity from external fixed-potential field
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2018/01/18 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <helper_cuda.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
/* The maximum number of elements in a chunk is 2^32-1 which is equal to 4,294,967,295 */
/* The maximum size for any chunk is 4GB */
#define MAXIMUM_CHUNK_SIZE      ((hsize_t)1 << 31)
#define MAXIMUM_CHUNK_SIZE_4BIT ((hsize_t)1 << 30)
#define MAXIMUM_CHUNK_SIZE_8BIT ((hsize_t)1 << 29)
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "cudalib.h"
#include "name.h"

#include "../misc/benchmark.h"
#include "../misc/structure.h"

#include "make.h"
#include "walk_dev.h"
#include "potential_dev.h"


/**
 * @fn allocSphericalPotentialTable_dev
 *
 * @brief Memory allocation on the accelerator device for external fixed-potential field by spherical components.
 */
extern "C"
muse allocSphericalPotentialTable_dev(real **rad, pot2 **Phi, const int Nr)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};

  mycudaMalloc((void **)rad, Nr * sizeof(real));  alloc.device += Nr * sizeof(real);
  mycudaMalloc((void **)Phi, Nr * sizeof(pot2));  alloc.device += Nr * sizeof(pot2);

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn  freeSphericalPotentialTable_dev
 *
 * @brief Memory deallocation on the accelerator device for external fixed-potential field by spherical components.
 */
extern "C"
void  freeSphericalPotentialTable_dev(real  *rad, pot2  *Phi)
{
  __NOTE__("%s\n", "start");

  mycudaFree(rad);
  mycudaFree(Phi);

  __NOTE__("%s\n", "end");
}


/**
 * @fn allocSphericalPotentialTable_hst
 *
 * @brief Memory allocation on the host for external fixed-potential field by spherical components.
 */
extern "C"
muse allocSphericalPotentialTable_hst(real **rad, pot2 **Phi, const int Nr)
{
  __NOTE__("%s\n", "start");

  muse alloc = {0, 0};

  mycudaMallocHost((void **)rad, Nr * sizeof(real));  alloc.host += Nr * sizeof(real);
  mycudaMallocHost((void **)Phi, Nr * sizeof(pot2));  alloc.host += Nr * sizeof(pot2);

  __NOTE__("%s\n", "end");
  return (alloc);
}


/**
 * @fn  freeSphericalPotentialTable_hst
 *
 * @brief Memory deallocation on the host for external fixed-potential field by spherical components.
 */
extern "C"
void  freeSphericalPotentialTable_hst(real  *rad, pot2  *Phi)
{
  __NOTE__("%s\n", "start");

  mycudaFreeHost(rad);
  mycudaFreeHost(Phi);

  __NOTE__("%s\n", "end");
}


/**
 * @fn setSphericalPotentialTable_dev
 *
 * @brief Set external fixed-potential field by spherical components on the device.
 */
extern "C"
void setSphericalPotentialTable_dev(real *rad_hst, pot2 *Phi_hst, real *rad_dev, pot2 *Phi_dev, const int Nr)
{
  __NOTE__("%s\n", "start");

  checkCudaErrors(cudaMemcpy(rad_dev, rad_hst, sizeof(real) * Nr, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(Phi_dev, Phi_hst, sizeof(pot2) * Nr, cudaMemcpyHostToDevice));

  __NOTE__("%s\n", "end");
}














/**
 * @fn calcExternalForce_spherical
 *
 * @brief Get F_r = -\nabla \Phi(r) and \Phi(r) by an external potential field at the specified point based on cubic spline interpolation in 1D.
 *
 * @return (F_r) F_r = -\nabla \Phi(r)
 * @return (Phi) \Phi(r)
 * @param (rr) position of the specified point
 * @param (logrmin) minimum of r-table in logarithmic scale
 * @param (invlogrbin) inverse of logrbin
 * @param (num) number of data points
 * @param (xx) position of data points
 * @param (yy) value of data points and coefficients in cubic spline interpolation
 */
__device__ __forceinline__
void calcExternalForce_spherical(real * RESTRICT F_r, real * RESTRICT Phi, const real rr, const real logrmin, const real invlogrbin, const int num, real * RESTRICT xx, pot2 * RESTRICT yy)
{
  int ii = (int)FLOOR((LOG10(rr) - logrmin) * invlogrbin);
  ii = (ii >=       0 ) ? ii :       0;
  ii = (ii < (num - 1)) ? ii : num - 2;

  const real xl = xx[ii];
  const real xr = xx[ii + 1];

  const real dx = xr - xl;
  const real dx2 = dx * dx;
  const real dxinv = RSQRT(dx2);
  const real aa = (rr - xl) * dxinv;

  const pot2 yl = yy[ii];
  const pot2 yr = yy[ii + 1];

  *Phi =  ((UNITY - aa) * yl.val + aa * (yr.val + (aa - UNITY) * dx2 * ((aa + UNITY) * yr.dr2 - (aa - TWO) * yl.dr2) * ONE_SIXTH));
  *F_r = -((yr.val - yl.val) * dxinv + dx * ((THREE * aa * aa - UNITY) * yr.dr2 - (TWO + THREE * aa * (aa - TWO)) * yl.dr2) * ONE_SIXTH);
}


/**
 * @fn calcExternalGravity_kernel
 *
 * @brief Calculate gravitational force from external potential field.
 *
 * @param (logrmin) minimum of r-table in logarithmic scale
 * @param (invlogrbin) inverse of logrbin
 * @param (num) number of data points
 * @param (xx) position of data points
 * @param (yy) value of data points and coefficients in cubic spline interpolation
 */
__global__ void calcExternalGravity_kernel
(acceleration * RESTRICT acc, READ_ONLY position * RESTRICT pos,
#ifdef  BLOCK_TIME_STEP
 const int laneNum, READ_ONLY laneinfo * RESTRICT laneInfo,
#endif//BLOCK_TIME_STEP
 const real logrmin, const real invlogrbin, const int num, real * RESTRICT xx, pot2 * RESTRICT yy)
{
#ifdef  BLOCK_TIME_STEP
  const int lane    = THREADIDX_X1D & (DIV_NWARP(TSUB) - 1);
  const int laneIdx = GLOBALIDX_X1D /  DIV_NWARP(TSUB);

  laneinfo info = {NUM_BODY_MAX, 0};
  if( laneIdx < laneNum )
    info = laneInfo[laneIdx];

  if( lane < info.num ){
    const int ii = info.head + lane;
#else///BLOCK_TIME_STEP
    const int ii = GLOBALIDX_X1D;
#endif//BLOCK_TIME_STEP

    /** load particle position and acceleration */
    const position pi = pos[ii];
    acceleration ai = acc[ii];

    /** set current location */
    const real r2 = 1.0e-30f + pi.x * pi.x + pi.y * pi.y + pi.z * pi.z;
    const real rinv = RSQRT(r2);
    const real rr = r2 * rinv;

    /** evaluate gravitational field from the external potential field */
    real F_r, Phi;
    calcExternalForce_spherical(&F_r, &Phi, rr, logrmin, invlogrbin, num, xx, yy);
    F_r *= rinv;

    /** accumulate the external force */
    ai.x   += F_r * ai.x;
    ai.y   += F_r * ai.y;
    ai.z   += F_r * ai.z;
    ai.pot += Phi;

    /** store acceleration */
    acc[ii] = ai;
#ifdef  BLOCK_TIME_STEP
  }/* if( lane < info.num ){ */
#endif//BLOCK_TIME_STEP
}


/**
 * @fn calcExternalGravity_dev
 *
 * @brief Calculate gravitational force from external potential field.
 */
extern "C"
void calcExternalGravity_dev
  (const iparticle pi, const potential_field sphe
#ifdef  BLOCK_TIME_STEP
   , const int thrd, const int grpNum, laneinfo * RESTRICT laneInfo
#else///BLOCK_TIME_STEP
   , const int Ni
#endif//BLOCK_TIME_STEP
   )
{
  __NOTE__("%s\n", "start");


#ifndef BLOCK_TIME_STEP
  const int thrd = NTHREADS;
#endif//BLOCK_TIME_STEP

  /** evaluate gravitational force by external potential field */
#ifdef  BLOCK_TIME_STEP
  int Nrem = BLOCKSIZE(grpNum, NWARP * NGROUPS);
#else///BLOCK_TIME_STEP
  int Nrem = BLOCKSIZE(Ni, NTHREADS);
#endif//BLOCK_TIME_STEP

  if( Nrem <= MAX_BLOCKS_PER_GRID ){
    /** when grid splitting is not required... */
#   if  defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
    if( grpNum != 0 )
#endif//defined(BLOCK_TIME_STEP) && !defined(SERIALIZED_EXECUTION)
      calcExternalGravity_kernel<<<Nrem, thrd>>>
	(pi.acc, pi.pos,
#ifdef  BLOCK_TIME_STEP
	 BLOCKSIZE(grpNum, NGROUPS) * NGROUPS, laneInfo,
#endif//BLOCK_TIME_STEP
	 sphe.logrmin, sphe.invlogrbin, sphe.num, sphe.rad, sphe.Phi);
  }/* if( Nrem <= MAX_BLOCKS_PER_GRID ){ */
  else{
    /** when grid splitting is required... */
    const int Niter = BLOCKSIZE(Nrem, MAX_BLOCKS_PER_GRID);
    int hidx = 0;

    for(int iter = 0; iter < Niter; iter++){
      int Nblck = MAX_BLOCKS_PER_GRID;
      if( Nblck > Nrem )	Nblck = Nrem;

#ifdef  BLOCK_TIME_STEP
      int Nsub = Nblck * NWARP * NGROUPS;
#else///BLOCK_TIME_STEP
      int Nsub = Nblck * NTHREADS;
#endif//BLOCK_TIME_STEP

      calcExternalGravity_kernel<<<Nblck, thrd>>>
	(
#ifdef  BLOCK_TIME_STEP
	 pi.acc, pi.pos, BLOCKSIZE(Nsub, NGROUPS) * NGROUPS, &laneInfo[hidx],
#else///BLOCK_TIME_STEP
	 &pi.acc[hidx], &pi.pos[hidx],
#endif//BLOCK_TIME_STEP
	 sphe.logrmin, sphe.invlogrbin, sphe.num, sphe.rad, sphe.Phi);

      hidx += Nsub;
      Nrem -= Nblck;
    }/* for(int iter = 0; iter < Niter; iter++){ */
  }/* else{ */

  getLastCudaError("calcExternalGravity_kernel");


  __NOTE__("%s\n", "end");
}


#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_SPHERICAL
/**
 * @fn  readFixedPotentialTableSpherical
 *
 * @brief Read fixed potential field (of spherical symmetric components) represented in cubic spline interpolation for GOTHIC.
 *
 * @param (unit) unit system of the potential field
 * @param (cfg) name of the input file
 * @return (pot_tbl) superposed potential field for cubic spline interpolation (only for spherical symmetric components)
 * @return (rad) radius for cubic spline interpolation
 * @return (Phi) potential field for cubic spline interpolation
 */
extern "C"
muse  readFixedPotentialTableSpherical(const int unit, char file[], potential_field *pot_tbl, real **rad, pot2 **Phi
#ifdef  USE_HDF5_FORMAT
				       , hdf5struct type
#endif//USE_HDF5_FORMAT
				       )
{
  __NOTE__("%s\n", "start");


  char cfgfile[128];
  sprintf(cfgfile, "%s/%s.%s.cfg", DOCUMENTFOLDER, file, "ext_pot_sphe");
  FILE *fp_cfg;
  fp_cfg = fopen(cfgfile, "r");
  if( fp_cfg == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", cfgfile);  }

  int Nread;
  bool success_cfg = true;
  success_cfg &= (1 == fscanf(fp_cfg, "%d", &Nread));

  /* open an existing file with read only option */
  char filename[128];
#ifdef  USE_HDF5_FORMAT

  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "ext_pot");
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* read attribute data */
  hid_t attribute;
  /* read flag about DOUBLE_PRECISION */
  int useDP;
  attribute = H5Aopen(target, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));
  /* read flag about unit system */
  int unit_pot;
  hid_t group = H5Gopen(target, "unit_system", H5P_DEFAULT);
  attribute = H5Aopen(group, "unit", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &unit_pot));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));

  /* simple error checks */
#ifdef  DOUBLE_PRECISION
  if( useDP != 1 ){    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);  }
#else///DOUBLE_PRECISION
  if( useDP != 0 ){    __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);  }
#endif//DOUBLE_PRECISION
  if( unit_pot != unit ){
    __KILL__(stderr, "ERROR: unit system of the potential field (%d) differs with that in the simulation run (%d)\n", unit_pot, unit);
  }/* if( unit_pot != unit ){ */

  group = H5Gopen(target, "spherical", H5P_DEFAULT);
  /* read # of data points */
  attribute = H5Aopen(group, "num", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &pot_tbl->num));
  chkHDF5err(H5Aclose(attribute));
  /* read log_10(r_min) */
  attribute = H5Aopen(group, "logrmin", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, type.real, &pot_tbl->logrmin));
  chkHDF5err(H5Aclose(attribute));
  /* read logrbin */
  attribute = H5Aopen(group, "logrbin", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, type.real, &pot_tbl->logrbin));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));

  const int num = pot_tbl->num;

  /* memory allocation on the accelerator device */
  muse alloc_tbl = allocSphericalPotentialTable_dev(rad, Phi, num);
  real *rad_hst;
  pot2 *Phi_hst;
  allocSphericalPotentialTable_hst(&rad_hst, &Phi_hst, num);

  /* memory allocation on the host as a temporary buffer */
  pot2 *Phi_tmp;
  if( Nread == 1 )
    Phi_tmp = Phi_hst;
  else{
    Phi_tmp = (pot2 *)malloc(num * sizeof(pot2));    if( Phi_tmp == NULL ){      __KILL__(stderr, "ERROR: failure to allocate Phi_tmp\n");    }

    const pot2 zero = {ZERO, ZERO};
    for(int ii = 0; ii < num; ii++)
      Phi_hst[ii] = zero;
  }/* else{ */


  for(int ii = 0; ii < Nread; ii++){
    hid_t dataset;

    /* open an existing group */
    char list[64];
    success_cfg &= (1 == fscanf(fp_cfg, "%s", list));
    group = H5Gopen(target, list, H5P_DEFAULT);

    /* read radius */
    if( ii == 0 ){
      dataset = H5Dopen(group, "r", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, rad_hst));
      chkHDF5err(H5Dclose(dataset));
    }/* if( ii == 0 ){ */

    /* read potential */
    dataset = H5Dopen(group, "Phi(r)", H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, type.pot2, H5S_ALL, H5S_ALL, H5P_DEFAULT, Phi_tmp));
    chkHDF5err(H5Dclose(dataset));

    chkHDF5err(H5Gclose(group));

    if( Nread > 1 )
      for(int jj = 0; jj < num; jj++){
	Phi_hst[jj].val += Phi_tmp[jj].val;
	Phi_hst[jj].dr2 += Phi_tmp[jj].dr2;
      }/* for(int jj = 0; jj < num; jj++){ */
  }/* for(int ii = 0; ii < Nread; ii++){ */
  if( Nread != 1 )
    free(Phi_tmp);

  /* close the file */
  chkHDF5err(H5Fclose(target));

#else///USE_HDF5_FORMAT

  /* read numeric table for superposed spherical components */
  sprintf(filename, "%s/%s.%s.%s", DATAFOLDER, file, "pot", "sphe");
  FILE *fp;
  fp = fopen(filename, "rb");
  if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

  int unit_pot, num;
  real rmin, rbin;

  bool success = true;
  size_t tmp;
  tmp = 1;  if( tmp != fread(&unit_pot, sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&num, sizeof(int), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&rmin, sizeof(real), tmp, fp) )    success = false;
  tmp = 1;  if( tmp != fread(&rbin, sizeof(real), tmp, fp) )    success = false;

  /* simple error check */
  if( unit_pot != unit ){
    __KILL__(stderr, "ERROR: unit system of the potential field (%d) differs with that in the simulation run (%d)\n", unit_pot, unit);
  }/* if( unit_pot != unit ){ */

  pot_tbl->num = num;
  pot_tbl->logrmin = rmin;
  pot_tbl->logrbin = rbin;

  /* memory allocation on the accelerator device */
  muse alloc_tbl = allocSphericalPotentialTable_dev(rad, Phi, num);
  real *rad_hst;
  pot2 *Phi_hst;
  allocSphericalPotentialTable_hst(&rad_hst, &Phi_hst, num);

  tmp = num;  if( tmp != fread(rad_hst, sizeof(real), tmp, fp) )    success = false;

  int list;
  success_cfg &= (1 == fscanf(fp_cfg, "%d", &list));
  if( (Nread == 1) && (list == READ_SUPERPOSED_TABLE_SPHE) ){
    tmp = num;    if( tmp != fread(Phi_hst, sizeof(pot2), tmp, fp) )      success = false;

    if( success != true ){      __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);    }
    fclose(fp);
  }/* if( (Nread == 1) && (list == READ_SUPERPOSED_TABLE_SPHE) ){ */
  else{
    /* close the superposed file */
    if( success != true ){      __KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);    }
    fclose(fp);

    pot2 *Phi_tmp;
    Phi_tmp = (pot2 *)malloc(num * sizeof(pot2));    if( Phi_tmp == NULL ){      __KILL__(stderr, "ERROR: failure to allocate Phi_tmp\n");    }

    const pot2 zero = {ZERO, ZERO};
    for(int ii = 0; ii < num; ii++)
      Phi_hst[ii] = zero;

    /* open another file */
    for(int ii = 0; ii < Nread; ii++){
      sprintf(filename, "%s/%s.%s.%d", DATAFOLDER, file, "pot", list);
      fp = fopen(filename, "rb");

      tmp = num;      if( tmp != fread(Phi_tmp, sizeof(pot2), tmp, fp) )	success = false;

      for(int jj = 0; jj < num; jj++){
	Phi_hst[jj].val += Phi_tmp[jj].val;
	Phi_hst[jj].dr2 += Phi_tmp[jj].dr2;
      }/* for(int jj = 0; jj < *Nr; jj++){ */

      if( success != true ){	__KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);      }
      fclose(fp);

      if( ii < (Nread - 1) )
	success_cfg &= (1 == fscanf(fp_cfg, "%d", &list));
    }/* for(int ii = 0; ii < Nread; ii++){ */
    free(Phi_tmp);
  }/* else{ */

#endif//USE_HDF5_FORMAT

  setSphericalPotentialTable_dev(rad_hst, Phi_hst, *rad, *Phi, num);
  freeSphericalPotentialTable_hst(rad_hst, Phi_hst);

  pot_tbl->rad = *rad;
  pot_tbl->Phi = *Phi;
  pot_tbl->invlogrbin = UNITY / pot_tbl->logrbin;

  if( success_cfg != true ){    __KILL__(stderr, "ERROR: failure to read \"%s\"\n", cfgfile);  }
  fclose(fp_cfg);


  __NOTE__("%s\n", "end");

  return (alloc_tbl);
}
#endif//SET_EXTERNAL_POTENTIAL_FIELD_SPHERICAL

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD_DISK

/**
 * @fn  readFixedPotentialTableDisk
 *
 * @brief Read fixed potential field (of disk components) represented in cubic spline interpolation for GOTHIC.
 *
 * @param (file) name of the input file
 * @param (Nlead) number of components to be read
 * @param (list) list of data groups to be read
 * @return (level)
 * @return (NR) number of data points for external potential-field
 * @return (Nz) number of data points for external potential-field
 * @return (RR) R for bilinear interpolation
 * @return (zz) z for bilinear interpolation
 * @return (Phi_Rz) 2-dimensional potential field for bilinear interpolation
 * @return (Nr) number of data points for (spherically averaged) external potential-field
 * @return (rr) r for cubic spline interpolation (for outermost regions)
 * @return (Phi_rr) potential field for cubic spline interpolation (for outermost regions)
 */
extern "C"
void  readFixedPotentialTableDisk(char file[], const int Nread
#ifdef  USE_HDF5_FORMAT
				  , char *list[], hdf5struct type
#else///USE_HDF5_FORMAT
				  , int *list
#endif//USE_HDF5_FORMAT
				  , int *level, int *NR, int *Nz, real **RR, real **zz, pot2 **Phi_Rz, int *Nr, real **rr, pot2 **Phi_rr
				  )
{
  __NOTE__("%s\n", "start");
  __NOTE__("%s\n", "end");
}


# of returned components must be 1;
if the number is greater than unity, then calculate superposed potential and return it;



# of components to be read, data group of each component
# of components to be passed to GOTHIC, 

 /* use table in (R, z)-plane for short- or middle-range force and r-plane for long-range force */

/* use nested-grid table and bilinear interpolation */

 /* add R = -dR and z = -dz to describe symmetry around R = 0, z = 0 */


#endif//SET_EXTERNAL_POTENTIAL_FIELD_DISK
