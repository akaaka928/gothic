#include <stdio.h>
#include <stdlib.h>

#include <hdf5.h>

#include "macro.h"
#include "constants.h"
#include "name.h"
#include "hdf5lib.h"

#include "profile.h"
/* #include "abel.h" */
/* #include "blas.h" */
/* #include "spline.h" */
#include "magi.h"
#include "potdens.h"


void read_disk_table(char *file, disk_data **disk, double **hor, double **ver, double **node_hor, double **node_ver, double **pot, double **rho);
void allocDiskArray(const int ndisk, disk_data **disk, const int maxLev, double **hor, double **ver, double **node_hor, double **node_ver, double **pot, double **rho);
void  freeDiskArray(disk_data *disk, double *hor, double *ver, double *node_hor, double *node_ver, double *pot, double *rho);


void read_disk_table(char *file, disk_data **disk, double **hor, double **ver, double **node_hor, double **node_ver, double **pot, double **rho)
{
  char filename[128];
  sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "disk");
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
#ifdef  DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION


  /** read attribute data */
  /** read # of components */
  int ndisk;
  hid_t attribute = H5Aopen(target, "kinds", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &ndisk));
  chkHDF5err(H5Aclose(attribute));
  /** read maximum # of nested level */
  int maxLev;
  attribute = H5Aopen(target, "maxLev", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &maxLev));
  chkHDF5err(H5Aclose(attribute));


  /** memory allocation */
  allocDiskArray(ndisk, disk, maxLev, hor, ver, node_hor, node_ver, pot, rho);


  static real tmp_array[NDISKBIN_HOR * NDISKBIN_VER];


  for(int ii = 0; ii < ndisk; ii++){
    char grp[16];    sprintf(grp, "data%d", ii);
    hid_t group = H5Gopen(target, grp, H5P_DEFAULT);

    /** write two-dimensional data */
    for(int lev = 0; lev < maxLev; lev++){
      char subgrp[16];      sprintf(subgrp, "patch%d", lev);
      hid_t patch = H5Gopen(group, subgrp, H5P_DEFAULT);

/* #pragma omp parallel for */
/*       for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++){ */
/* 	tmp_rho[jj] = CAST_D2R((*disk[ii].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * density2astro); */
/* 	tmp_Phi[jj] = CAST_D2R(	 disk[ii].pot [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] * senergy2astro); */
/*       }/\* for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++){ *\/ */

      /** read density */
      hid_t dataset;
      dataset = H5Dopen(patch, "rho", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_array));
      chkHDF5err(H5Dclose(dataset));
      for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++)
	(*(*disk)[ii].rho)[INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] = CAST_R2D(tmp_array[jj]);
      /** write potential */
      dataset = H5Dopen(patch, "Phi", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_array));
      chkHDF5err(H5Dclose(dataset));
      for(int jj = 0; jj < NDISKBIN_HOR * NDISKBIN_VER; jj++)
	(*disk)[ii].pot [INDEX2D(maxLev, NDISKBIN_HOR * NDISKBIN_VER, lev, jj)] = CAST_R2D(tmp_array[jj]);


/* #pragma omp parallel */
/*       { */
/* #pragma omp for nowait */
/* 	for(int jj = 0; jj < NDISKBIN_HOR; jj++) */
/* 	  tmp_hor[jj] = CAST_D2R(disk[ii].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] * length2astro); */
/* #pragma omp for nowait */
/* 	for(int jj = 0; jj < NDISKBIN_VER; jj++) */
/* 	  tmp_ver[jj] = CAST_D2R(disk[ii].ver[INDEX2D(maxLev, NDISKBIN_VER, lev, jj)] * length2astro); */
/* #pragma omp for nowait */
/* 	for(int jj = 0; jj < NDISKBIN_HOR + 1; jj++) */
/* 	  node_RR[jj] = CAST_D2R(disk[ii].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, jj)] * length2astro); */
/* #pragma omp for nowait */
/* 	for(int jj = 0; jj < NDISKBIN_VER + 1; jj++) */
/* 	  node_zz[jj] = CAST_D2R(disk[ii].node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev, jj)] * length2astro); */
/*       } */

      /** read horizontal position */
      dataset = H5Dopen(patch, "radius", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_array));
      chkHDF5err(H5Dclose(dataset));
      for(int jj = 0; jj < NDISKBIN_HOR; jj++)
	(*disk)[ii].hor[INDEX2D(maxLev, NDISKBIN_HOR, lev, jj)] = CAST_R2D(tmp_array[jj]);

      /** read horizontal position */
      dataset = H5Dopen(patch, "node_RR", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_array));
      chkHDF5err(H5Dclose(dataset));
      for(int jj = 0; jj < NDISKBIN_HOR + 1; jj++)
	(*disk)[ii].node_hor[INDEX2D(maxLev, NDISKBIN_HOR + 1, lev, jj)] = CAST_R2D(tmp_array[jj]);

      /** read vertical position */
      dataset = H5Dopen(patch, "height", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_array));
      chkHDF5err(H5Dclose(dataset));
      for(int jj = 0; jj < NDISKBIN_VER; jj++)
	(*disk)[ii].ver[INDEX2D(maxLev, NDISKBIN_VER, lev, jj)] = CAST_R2D(tmp_array[jj]);

      /** read vertical position */
      dataset = H5Dopen(patch, "node_zz", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_array));
      chkHDF5err(H5Dclose(dataset));
      for(int jj = 0; jj < NDISKBIN_VER + 1; jj++)
	(*disk)[ii].node_ver[INDEX2D(maxLev, NDISKBIN_VER + 1, lev, jj)] = CAST_R2D(tmp_array[jj]);


      chkHDF5err(H5Gclose(patch));
    }/* for(int lev = 0; lev < maxLev; lev++){ */


    chkHDF5err(H5Gclose(group));
  }/* for(int ii = 0; ii < ndisk; ii++){ */



  /** close the file */
  chkHDF5err(H5Fclose(target));

}


/**
 * @fn allocDiskProfile
 *
 * @brief Allocate memory for disk component(s)
 *
 * @param (ndisk) number of disk components
 * @return (disk) physical quantities of the disk component
 * @return (maxLev) maximum level of nested grid
 * @return (hor) RR, zone center value
 * @return (ver) zz, zone center value
 * @return (node_hor) RR, node center value
 * @return (node_ver) zz, node center value
 * @return (pot) potential field
 * @return (rho) density field
 */
void allocDiskArray(const int ndisk, disk_data **disk, const int maxLev, double **hor, double **ver, double **node_hor, double **node_ver, double **pot, double **rho)
{
  __NOTE__("%s\n", "start");


  /** allocate arrays to store physical quantities */
  /** horizontal and vertical axes */
  *hor      = (double *)malloc(maxLev *  NDISKBIN_HOR      * sizeof(double));  if( *     hor == NULL ){    __KILL__(stderr, "ERROR: failure to allocate hor\n");  }
  *ver      = (double *)malloc(maxLev *  NDISKBIN_VER      * sizeof(double));  if( *     ver == NULL ){    __KILL__(stderr, "ERROR: failure to allocate ver\n");  }
  *node_hor = (double *)malloc(maxLev * (NDISKBIN_HOR + 1) * sizeof(double));  if( *node_hor == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node_hor\n");  }
  *node_ver = (double *)malloc(maxLev * (NDISKBIN_VER + 1) * sizeof(double));  if( *node_ver == NULL ){    __KILL__(stderr, "ERROR: failure to allocate node_ver\n");  }

  /** potential-density pair */
  *pot = (double *)malloc(        maxLev * NDISKBIN_HOR *  NDISKBIN_VER      * sizeof(double));  if( *pot == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pot\n" );  }
  *rho = (double *)malloc(ndisk * maxLev * NDISKBIN_HOR * (NDISKBIN_VER + 1) * sizeof(double));  if( *rho == NULL ){    __KILL__(stderr, "ERROR: failure to allocate rho0\n");  }

  /** allocate utility structure and commit arrays */
  *disk = (disk_data *)malloc(ndisk * sizeof(disk_data));  if( *disk == NULL ){    __KILL__(stderr, "ERROR: failure to allocate disk\n");  }

  for(int ii = 0; ii < ndisk; ii++){
    /** common arrays for all components */
    (*disk)[ii].hor = *hor;
    (*disk)[ii].ver = *ver;
    (*disk)[ii].node_hor = *node_hor;
    (*disk)[ii].node_ver = *node_ver;
    (*disk)[ii].pot = *pot;

    /** individual arrays for each component */
    (*disk)[ii].rho0 = &((*rho)[INDEX(ndisk, maxLev, NDISKBIN_HOR * (NDISKBIN_VER + 1), ii, 0, 0)]);
    (*disk)[ii].rho = &(*disk)[ii].rho0;
  }/* for(int ii = 0; ii < ndisk; ii++){ */


  __NOTE__("%s\n", "end");
}


/**
 * @fn freeDiskProfile
 *
 * @brief Deallocate memory for disk component(s)
 *
 * @param (disk) physical quantities of the disk component
 * @param (hor) RR, zone center value
 * @param (ver) zz, zone center value
 * @param (node_hor) RR, node center value
 * @param (node_ver) zz, node center value
 * @param (pot) potential field
 * @param (rho) density field
 */
void  freeDiskArray(disk_data *disk, double *hor, double *ver, double *node_hor, double *node_ver, double *pot, double *rho)
{
  __NOTE__("%s\n", "start");


  free(hor);
  free(ver);
  free(node_hor);
  free(node_ver);

  free(pot);
  free(rho);

  free(disk);


  __NOTE__("%s\n", "end");
}


#include "myutil.h"
int main(int argc, char **argv)
{
  /** read input arguments */
  if( argc < 2 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 2);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 2 ){ */

  /** read input arguments do not depend on the unit system adopted in the numerical simulation */
  char *file;
  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv, "file", &file));


  disk_data *disk;
  double *hor, *ver, *node_hor, *node_ver;
  double *pot, *rho;
  read_disk_table(file, &disk, &hor, &ver, &node_hor, &node_ver, &pot, &rho);


  freeDiskArray(disk, hor, ver, node_hor, node_ver, pot, rho);


  return (0);
}
