/**
 * @file error.c
 *
 * @brief Source code for analyzing error of force calculation
 *
 * @author Yohei Miki (University of Tokyo)
 * @author Masayuki Umemura (University of Tsukuba)
 *
 * @date 2019/07/31 (Wed)
 *
 * Copyright (C) 2017 Yohei Miki and Masayuki Umemura
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#define USE_SZIP_COMPRESSION

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "myutil.h"
#include "constants.h"
#include "mpilib.h"
#include "name.h"

#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../file/io.h"

#ifdef  USE_HDF5_FORMAT
extern const double length2astro;extern const char length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double   time2astro;extern const char   time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
#endif//USE_HDF5_FORMAT


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);


  /** initialization */
  if( argc < 7 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 7);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file =<char *>\n");
    __FPRINTF__(stderr, "          -file0=<char *>\n");
    __FPRINTF__(stderr, "          -file1=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 7 ){ */
  char *file ;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "file" , &file ));
  char *file0;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "file0", &file0));
  char *file1;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "file1", &file1));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));


  /** prepare dataset using the reference solution (direct method or single process) */
  int last;
  readConfigFile(&last, file0);
  int unit;
  ulong Ntot;
  real eta, eps;
  double ft, snapshotInterval, saveInterval;
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file0, mpi);

#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
  nbody_hdf5 body0;  real *pos0, *vel0, *acc0, *m0, *pot0;  ulong *idx0;
  nbody_hdf5 body1;  real *pos1, *vel1, *acc1, *m1, *pot1;  ulong *idx1;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  real *acc0_ext, *pot0_ext;
  real *acc1_ext, *pot1_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  allocSnapshotArray(&pos0, &vel0, &acc0, &m0, &pot0, &idx0,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
		     &acc0_ext, &pot0_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
		     (int)Ntot, &body0);
  allocSnapshotArray(&pos1, &vel1, &acc1, &m1, &pot1, &idx1,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
		     &acc1_ext, &pot1_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
		     (int)Ntot, &body1);
#else///USE_HDF5_FORMAT
  iparticle body0, body1;
  ulong *idx0, *idx1;
  position *pos0, *pos1;
  acceleration *acc0, *acc1;
#ifdef  BLOCK_TIME_STEP
  velocity *vel0, *vel1;
  ibody_time *ti0, *ti1;
#else///BLOCK_TIME_STEP
  real *vx0, *vy0, *vz0, *vx1, *vy1, *vz1;
#endif//BLOCK_TIME_STEP
  allocParticleData((int)Ntot, &body0, &idx0, &pos0, &acc0,
#ifdef  BLOCK_TIME_STEP
		    &vel0, &ti0
#else///BLOCK_TIME_STEP
		    &vx0, &vy0, &vz0
#endif//BLOCK_TIME_STEP
		    );
  allocParticleData((int)Ntot, &body1, &idx1, &pos1, &acc1,
#ifdef  BLOCK_TIME_STEP
		    &vel1, &ti0
#else///BLOCK_TIME_STEP
		    &vx0, &vy0, &vz0
#endif//BLOCK_TIME_STEP
		    );
#endif//USE_HDF5_FORMAT

  real *axx_err;  axx_err = (real *)malloc(Ntot * sizeof(real));  if( axx_err == NULL ){    __KILL__(stderr, "ERROR: failure to allocate axx_err\n");  }
  real *ayy_err;  ayy_err = (real *)malloc(Ntot * sizeof(real));  if( ayy_err == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ayy_err\n");  }
  real *azz_err;  azz_err = (real *)malloc(Ntot * sizeof(real));  if( azz_err == NULL ){    __KILL__(stderr, "ERROR: failure to allocate azz_err\n");  }
  real *pot_err;  pot_err = (real *)malloc(Ntot * sizeof(real));  if( pot_err == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pot_err\n");  }
  real *pot_abs;  pot_abs = (real *)malloc(Ntot * sizeof(real));  if( pot_abs == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pot_abs\n");  }

  int *tag0;  tag0 = (int *)malloc(Ntot * sizeof(int));  if( tag0 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tag0\n");  }
  int *tag1;  tag1 = (int *)malloc(Ntot * sizeof(int));  if( tag1 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate tag1\n");  }

  real *pos;  pos = (real *)malloc(Ntot * 3 * sizeof(real));  if( pos == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pos\n");  }
  real *acc;  acc = (real *)malloc(Ntot * 3 * sizeof(real));  if( acc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate acc\n");  }


  for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){
    /** read simulation results */
    double time;
    ulong steps;
    int unit_read;
#ifdef  USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file1, (uint)filenum, &body1, hdf5type);
    if( unit_read != unit ){
      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
    }/* if( unit_read != unit ){ */
    readSnapshot(&unit_read, &time, &steps, Ntot, file0, (uint)filenum, &body0, hdf5type);
    if( unit_read != unit ){
      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
    }/* if( unit_read != unit ){ */
#else///USE_HDF5_FORMAT
    readSnapshot(&unit_read, &time, &steps, Ntot, file1, (uint)filenum, body1);
    if( unit_read != unit ){
      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
    }/* if( unit_read != unit ){ */
    readSnapshot(&unit_read, &time, &steps, Ntot, file0, (uint)filenum, body0);
    if( unit_read != unit ){
      __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
    }/* if( unit_read != unit ){ */
#endif//USE_HDF5_FORMAT

    for(int ii = 0; ii < (int)Ntot; ii++)      tag0[body0.idx[ii]] = ii;
    for(int ii = 0; ii < (int)Ntot; ii++)      tag1[body1.idx[ii]] = ii;


    /** error analysis */
    int potIdx = 0;    real potMax = ZERO;
    int accIdx = 0;    real accMax = ZERO;
    for(int ii = 0; ii < (int)Ntot; ii++){
      const int idx0 = tag0[ii];
      const int idx1 = tag1[ii];

#ifdef  USE_HDF5_FORMAT
      pos[    3 * ii] = body0.pos[    3 * idx0];
      pos[1 + 3 * ii] = body0.pos[1 + 3 * idx0];
      pos[2 + 3 * ii] = body0.pos[2 + 3 * idx0];
#else///USE_HDF5_FORMAT
      pos[    3 * ii] = body0.pos[idx0].x;
      pos[1 + 3 * ii] = body0.pos[idx0].y;
      pos[2 + 3 * ii] = body0.pos[idx0].z;
#endif//USE_HDF5_FORMAT

      pot_err[ii] = body1.pot[idx1] / body0.pot[idx0] - UNITY;
      pot_abs[ii] = FABS(pot_err[ii]);
#ifdef  USE_HDF5_FORMAT
      const real ax0 = body0.acc[    3 * idx0];      const real ax1 = body1.acc[    3 * idx1];
      const real ay0 = body0.acc[1 + 3 * idx0];      const real ay1 = body1.acc[1 + 3 * idx1];
      const real az0 = body0.acc[2 + 3 * idx0];      const real az1 = body1.acc[2 + 3 * idx1];
#else///USE_HDF5_FORMAT
      const real ax0 = body0.acc[idx0].x;      const real ax1 = body1.acc[idx1].x;
      const real ay0 = body0.acc[idx0].y;      const real ay1 = body1.acc[idx1].y;
      const real az0 = body0.acc[idx0].z;      const real az1 = body1.acc[idx1].z;
#endif//USE_HDF5_FORMAT
      axx_err[ii] = ax1 / ax0 - UNITY;      const real dx = ax1 - ax0;
      ayy_err[ii] = ay1 / ay0 - UNITY;      const real dy = ay1 - ay0;
      azz_err[ii] = az1 / az0 - UNITY;      const real dz = az1 - az0;
      const real a0inv = RSQRT(FLT_MIN + ax0 * ax0 + ay0 * ay0 + az0 * az0);
      acc[    3 * ii] = dx * a0inv;
      acc[1 + 3 * ii] = dy * a0inv;
      acc[2 + 3 * ii] = dz * a0inv;

      real acc_abs = (dx * dx + dy * dy + dz * dz) * a0inv * a0inv;
      if( accMax < acc_abs     ){	accMax = acc_abs    ;	accIdx = ii;      }
      if( potMax < pot_abs[ii] ){	potMax = pot_abs[ii];	potIdx = ii;      }
    }/* for(int ii = 0; ii < (int)Ntot; ii++){ */


    /** output error distribution */
    char filename[128];

#ifdef  USE_HDF5_FORMAT
    /** create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
    sprintf(filename, "%s/%s.%s%.3d.h5", DATAFOLDER, file, ACCERR, filenum);
    hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /** write particle data */
    /** create a group and data space */
    char groupname[16];
    sprintf(groupname, ACCERR);
    hid_t group = H5Gcreate(target, groupname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    /** preparation for data compression */
    hid_t dataset, dataspace, property;
#ifdef  USE_SZIP_COMPRESSION
    /** compression using szip */
    uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
    uint szip_pixels_per_block = 8;
    hsize_t cdims[2] = {32 * szip_pixels_per_block, 3};
#else///USE_SZIP_COMPRESSION
    property = H5P_DEFAULT;
#endif//USE_SZIP_COMPRESSION
    /** 2D (Ntot * 3) array */
    hsize_t dims[2] = {Ntot, 3};
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t cdims_loc[2] = {cdims[0], cdims[1]};
    if( (hsize_t)Ntot < cdims_loc[0] )
      cdims_loc[0] = (hsize_t)Ntot;
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
    /* write particle position */
    dataset = H5Dcreate(group, "position", hdf5type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos));
    chkHDF5err(H5Dclose(dataset));
    /** write acceleration error as vector */
    dataset = H5Dcreate(group, "acceleration", hdf5type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, acc));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    /** 1D (num) arrays */
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_SZIP_COMPRESSION
    property = H5Pcreate(H5P_DATASET_CREATE);
    cdims[0] = 128 * szip_pixels_per_block;
    cdims_loc[0] = cdims[0];
    if( (hsize_t)Ntot < cdims_loc[0] )
      cdims_loc[0] = (hsize_t)Ntot;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
#endif//USE_SZIP_COMPRESSION
    /** write relative error about ax */
    dataset = H5Dcreate(group, "ax", hdf5type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, axx_err));
    chkHDF5err(H5Dclose(dataset));
    /** write relative error about ay */
    dataset = H5Dcreate(group, "ay", hdf5type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, ayy_err));
    chkHDF5err(H5Dclose(dataset));
    /** write relative error about az */
    dataset = H5Dcreate(group, "az", hdf5type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, azz_err));
    chkHDF5err(H5Dclose(dataset));
    /** write relative error about pot */
    dataset = H5Dcreate(group, "pot", hdf5type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, pot_err));
    chkHDF5err(H5Dclose(dataset));
    /** write relative error about pot_abs */
    dataset = H5Dcreate(group, "pot_abs", hdf5type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, hdf5type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, pot_abs));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_SZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_SZIP_COMPRESSION
    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    /** write attribute data */
    hsize_t attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attribute;
    /** write current time */
    attribute = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time));
    chkHDF5err(H5Aclose(attribute));
    /** write # of steps */
    attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
    chkHDF5err(H5Aclose(attribute));
    /** write # of N-body particles */
    attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &Ntot));
    chkHDF5err(H5Aclose(attribute));
    /** write flag about USE_DOUBLE_PRECISION */
#ifdef  USE_DOUBLE_PRECISION
    const int useDP = 1;
#else///USE_DOUBLE_PRECISION
    const int useDP = 0;
#endif//USE_DOUBLE_PRECISION
    attribute = H5Acreate(group, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
    chkHDF5err(H5Aclose(attribute));
    /** close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    /** close the group */
    chkHDF5err(H5Gclose(group));


    /** write information on the worst particle */
    group = H5Gcreate(target, "worst", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    /** write maximum error of acceleration */
    attribute = H5Acreate(group, "error (acc)", hdf5type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, hdf5type.real, &accMax));
    chkHDF5err(H5Aclose(attribute));
    /** write particle index having maximum error of acceleration */
    attribute = H5Acreate(group, "ID (acc)", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &accIdx));
    chkHDF5err(H5Aclose(attribute));
    /** write maximum error of potential */
    attribute = H5Acreate(group, "error (pot)", hdf5type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, hdf5type.real, &potMax));
    chkHDF5err(H5Aclose(attribute));
    /** write particle index having maximum error of potential */
    attribute = H5Acreate(group, "ID (pot)", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &potIdx));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Sclose(dataspace));
    attr_dims = 3;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    /** write particle position having maximum error of acceleration */
    attribute = H5Acreate(group, "pos (acc)", hdf5type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, hdf5type.real, &pos[3 * accIdx]));
    chkHDF5err(H5Aclose(attribute));
    /** write particle position having maximum error of potential */
    attribute = H5Acreate(group, "pos (pot)", hdf5type.real, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, hdf5type.real, &pos[3 * potIdx]));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(group));

    /** write unit system */
    group = H5Gcreate(target, "unit system", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    /** create the data space for the attribute */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    /** write conversion factors */
    hid_t subgroup = H5Gcreate(group, "conversion factors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    /** create the data space for the attribute */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    attribute = H5Acreate(subgroup, "length2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &length2astro));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(subgroup, "time2astro", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time2astro));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(subgroup));
    /** write axis labels */
    subgroup = H5Gcreate(group, "axis labels", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    /** create the data space for the attribute */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    attribute = H5Acreate(subgroup, "length_astro_unit_name", hdf5type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, hdf5type.str4unit, length_astro_unit_name));
    chkHDF5err(H5Aclose(attribute));
    attribute = H5Acreate(subgroup, "time_astro_unit_name", hdf5type.str4unit, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, hdf5type.str4unit, time_astro_unit_name));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Sclose(dataspace));
    chkHDF5err(H5Gclose(subgroup));
    /** close the group for unit system */
    chkHDF5err(H5Gclose(group));
    /** close the file */
    chkHDF5err(H5Fclose(target));

#else///USE_HDF5_FORMAT

    sprintf(filename, "%s/%s.%s%.3d.dat", DATAFOLDER, file, ACCERR, filenum);
    FILE *fp;
    fp = fopen(filename, "wb");
    if( fp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);    }
    bool success = true;
    size_t tmp;
    tmp =        1;    if( tmp != fwrite(& time , sizeof(double), tmp, fp) )      success = false;
    tmp =        1;    if( tmp != fwrite(&steps , sizeof( ulong), tmp, fp) )      success = false;
    tmp =        1;    if( tmp != fwrite(&accIdx, sizeof(   int), tmp, fp) )      success = false;
    tmp =        1;    if( tmp != fwrite(&accMax, sizeof(  real), tmp, fp) )      success = false;
    tmp =        1;    if( tmp != fwrite(&potIdx, sizeof(   int), tmp, fp) )      success = false;
    tmp =        1;    if( tmp != fwrite(&potMax, sizeof(  real), tmp, fp) )      success = false;
    tmp = Ntot * 3;    if( tmp != fwrite(pos    , sizeof(  real), tmp, fp) )      success = false;
    tmp = Ntot * 3;    if( tmp != fwrite(acc    , sizeof(  real), tmp, fp) )      success = false;
    tmp = Ntot    ;    if( tmp != fwrite(pot_err, sizeof(  real), tmp, fp) )      success = false;
    tmp = Ntot    ;    if( tmp != fwrite(pot_abs, sizeof(  real), tmp, fp) )      success = false;
    tmp = Ntot    ;    if( tmp != fwrite(axx_err, sizeof(  real), tmp, fp) )      success = false;
    tmp = Ntot    ;    if( tmp != fwrite(ayy_err, sizeof(  real), tmp, fp) )      success = false;
    tmp = Ntot    ;    if( tmp != fwrite(azz_err, sizeof(  real), tmp, fp) )      success = false;
    if( success != true ){      __KILL__(stderr, "ERROR: failure to write \"%s\"\n", filename);    }
    fclose(fp);

#endif//USE_HDF5_FORMAT
  }/* for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){ */


  free(pos);  free(acc);
  free(tag0);  free(tag1);
  free(axx_err);  free(ayy_err);  free(azz_err);
  free(pot_err);  free(pot_abs);

#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
  freeSnapshotArray(pos0, vel0, acc0, m0, pot0, idx0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
		    , acc0_ext, pot0_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
		    );
  freeSnapshotArray(pos1, vel1, acc1, m1, pot1, idx1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
		    , acc1_ext, pot1_ext
#endif//SET_EXTERNAL_POTENTIAL_FIELD
		    );
#else///USE_HDF5_FORMAT
  freeParticleData(idx, pos, iacc,
#ifdef  BLOCK_TIME_STEP
		    vel, ti
#else///BLOCK_TIME_STEP
		    vx, vy, vz
#endif//BLOCK_TIME_STEP
		    );
#endif//USE_HDF5_FORMAT


  exitMPI();
  return (0);
}
