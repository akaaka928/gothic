/*************************************************************************\
 *                                                                       *
                  last updated on 2016/04/20(Wed) 16:26:09
 *                                                                       *
 *    Plot Code of Cumulative distribution function for Tree code        *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
//-------------------------------------------------------------------------
#ifdef  USE_HDF5_FORMAT
#       include <hdf5.h>
#       include <hdf5lib.h>
#endif//USE_HDF5_FORMAT
//-------------------------------------------------------------------------
#include <macro.h>
#include <myutil.h>
#include <constants.h>
#include <mpilib.h>
#include <plplotlib.h>
#include <name.h>
//-------------------------------------------------------------------------
#include "../misc/structure.h"
//-------------------------------------------------------------------------
#include "../file/io.h"
//-------------------------------------------------------------------------
#include "cdflib.h"
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#define PLOT_SUMMARY_FIGURE_ONLY
//-------------------------------------------------------------------------
/* #define PRINT_VALUES */
//-------------------------------------------------------------------------
/* #define ERROR_DETECT */
//-------------------------------------------------------------------------
#define LOGPLOT_HOR
/* #define LOGPLOT_VER */
//-------------------------------------------------------------------------
#   if  defined(PRINT_VALUES) || defined(ERROR_DETECT)
#include "../misc/allocate.h"
#define ERRFILE "err"
FILE *fp;
nbody_particle *body;
#endif//defined(PRINT_VALUES) || defined(ERROR_DETECT)
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotCDF
(int num, PLFLT *cdf, PLFLT *acc, PLFLT *grv, PLFLT *pot, PLplotPltRange range,
 char file[], const int filenum, int argc, char **argv);
//-------------------------------------------------------------------------
void chkTreeError(int num, acceleration *direct, acceleration *octree,
		  double *accErr, double *grvErr, double *potErr);
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef __ICC
/* Disable ICC's remark #161: unrecognized #pragma */
#     pragma warning (disable:161)
#endif//__ICC
//-------------------------------------------------------------------------
int idxAscendingOrder(const void *a, const void *b)
{
  //-----------------------------------------------------------------------
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
  if(          ((nbody_particle *)a)->idx > ((nbody_particle *)b)->idx ){    return ( 1);  }
  else{    if( ((nbody_particle *)a)->idx < ((nbody_particle *)b)->idx ){    return (-1);  }
    else{                                                                    return ( 0);  }  }
#pragma GCC diagnostic pop
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifdef __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //-----------------------------------------------------------------------
  /* parallelized region employing MPI start */
  //-----------------------------------------------------------------------
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* initialization */
  /* setPhysicalConstantsAndUnitSystem(UNITSYSTEM, 0); */
  //-----------------------------------------------------------------------
  if( argc < 2 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 2);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }
  //-----------------------------------------------------------------------
  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  /* int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start)); */
  /* int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end)); */
  /* int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval)); */
  //-----------------------------------------------------------------------
  modifyArgcArgv4PLplot(&argc, argv, 2);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* prepare dataset */
  //-----------------------------------------------------------------------
  int last;
  readConfigFile(&last, file);
  int unit;
  ulong Ntot;
  real eta, eps;
  double ft, snapshotInterval, saveInterval;
  /* readSettings(&Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file); */
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file, mpi);
  //-----------------------------------------------------------------------
  ulong        * index;  index  = (ulong        *)malloc(sizeof(acceleration) * (size_t)Ntot);
  acceleration *accBuf;  accBuf = (acceleration *)malloc(sizeof(acceleration) * (size_t)Ntot);
  acceleration *direct;  direct = (acceleration *)malloc(sizeof(acceleration) * (size_t)Ntot);
  acceleration *octree;  octree = (acceleration *)malloc(sizeof(acceleration) * (size_t)Ntot);
  if( index  == NULL ){    __KILL__(stderr, "ERROR: failure to allocate  index");  }
  if( accBuf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate accBuf");  }
  if( direct == NULL ){    __KILL__(stderr, "ERROR: failure to allocate direct");  }
  if( octree == NULL ){    __KILL__(stderr, "ERROR: failure to allocate octree");  }
  //-----------------------------------------------------------------------
  double *potErr;  potErr = (double *)malloc(sizeof(double) * (size_t)Ntot);  if( potErr == NULL ){    __KILL__(stderr, "ERROR: failure to allocate potErr");  }
  double *grvErr;  grvErr = (double *)malloc(sizeof(double) * (size_t)Ntot);  if( grvErr == NULL ){    __KILL__(stderr, "ERROR: failure to allocate grvErr");  }
  double *accErr;  accErr = (double *)malloc(sizeof(double) * (size_t)Ntot);  if( accErr == NULL ){    __KILL__(stderr, "ERROR: failure to allocate accErr");  }
  //-----------------------------------------------------------------------
  PLFLT *cdf;  cdf = (PLFLT *)malloc(sizeof(PLFLT) * (size_t)Ntot);  if( cdf == NULL ){    __KILL__(stderr, "ERROR: failure to allocate cdf");  }
  PLFLT *pot;  pot = (PLFLT *)malloc(sizeof(PLFLT) * (size_t)Ntot);  if( pot == NULL ){    __KILL__(stderr, "ERROR: failure to allocate pot");  }
  PLFLT *grv;  grv = (PLFLT *)malloc(sizeof(PLFLT) * (size_t)Ntot);  if( grv == NULL ){    __KILL__(stderr, "ERROR: failure to allocate grv");  }
  PLFLT *acc;  acc = (PLFLT *)malloc(sizeof(PLFLT) * (size_t)Ntot);  if( acc == NULL ){    __KILL__(stderr, "ERROR: failure to allocate acc");  }
  //-----------------------------------------------------------------------
  double *percentage;
  allocPercentile(&percentage);
  //-----------------------------------------------------------------------
#   if  defined(PRINT_VALUES) || defined(ERROR_DETECT)
  allocParticleDataAoS((int)Ntot, &body);
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
  nbody_hdf5 hdf5_dat;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
  ulong *hdf5_idx;
  allocSnapshotArray(&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx, (int)Ntot, &hdf5_dat);
#endif//USE_HDF5_FORMAT
#endif//defined(PRINT_VALUES) || defined(ERROR_DETECT)
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  int checker = 1;
  int lstNum = 1;
  char lstfile[128];
  sprintf(lstfile, "%s/%s.%s.txt", DOCUMENTFOLDER, file, ERRFILE_NUM);
  FILE *lstfp;
  if( mpi.rank == 0 ){
    lstfp = fopen(lstfile, "r");
    if( lstfp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", lstfile);    }
    checker &= (1 == fscanf(lstfp, "%d", &lstNum));
    if( !checker ){      __KILL__(stderr, "ERROR: read failure from \"%s\"\n", lstfile);    }
    fclose(lstfp);
  }/* if( mpi.rank == 0 ){ */
  chkMPIerr(MPI_Bcast(&lstNum, 1, MPI_INT, 0, mpi.comm));
  //-----------------------------------------------------------------------
  sprintf(lstfile, "%s/%s.%s.txt", DOCUMENTFOLDER, file, ERRFILE_LST);
  if( mpi.rank == 0 ){
    lstfp = fopen(lstfile, "r");
    if( lstfp == NULL ){      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", lstfile);    }
  }/* if( mpi.rank == 0 ){ */
  //-----------------------------------------------------------------------
  for(int lstIdx = 0; lstIdx < lstNum; lstIdx++){
    //---------------------------------------------------------------------
    char accfile[16], macname[16];
    int nfile = 0;
    real mac = UNITY;
    if( mpi.rank == 0 ){
      checker &= (4 == fscanf(lstfp, "%s\t%d\t%s\t%e", accfile, &nfile, macname, &mac));
      if( !checker ){	__KILL__(stderr, "ERROR: read failure from \"%s\"\n", lstfile);      }
    }/* if( mpi.rank == 0 ){ */
    chkMPIerr(MPI_Bcast(& nfile,  1, MPI_INT    , 0, mpi.comm));
    chkMPIerr(MPI_Bcast(&   mac,  1, MPI_REALDAT, 0, mpi.comm));
    chkMPIerr(MPI_Bcast(accfile, 16, MPI_CHAR   , 0, mpi.comm));
    chkMPIerr(MPI_Bcast(macname, 16, MPI_CHAR   , 0, mpi.comm));
    //---------------------------------------------------------------------
    char serfile[128];
    sprintf(serfile, "%s.%s", file, accfile);
    //---------------------------------------------------------------------
    double *potErrList, *grvErrList, *accErrList;
    if( lstNum == 1 ){
      potErrList = (double *)malloc(sizeof(double) * (size_t)nfile * NSUMMARY);
      grvErrList = (double *)malloc(sizeof(double) * (size_t)nfile * NSUMMARY);
      accErrList = (double *)malloc(sizeof(double) * (size_t)nfile * NSUMMARY);
      if( potErrList == NULL ){
	__KILL__(stderr, "ERROR: failure to allocate potErrList");
      }
      if( grvErrList == NULL ){
	__KILL__(stderr, "ERROR: failure to allocate grvErrList");
      }
      if( accErrList == NULL ){
	__KILL__(stderr, "ERROR: failure to allocate accErrList");
      }
      for(int ii = 0; ii < nfile * NSUMMARY; ii++){
	potErrList[ii] = 0.0;
	grvErrList[ii] = 0.0;
	accErrList[ii] = 0.0;
      }/* for(int ii = 0; ii < nfile * NSUMMARY; ii++){ */
    }/* if( lstNum == 1 ){ */
    //---------------------------------------------------------------------
    int ifile = mpi.rank;
    for(int filenum = mpi.rank; filenum < nfile; filenum += mpi.size){
      //-------------------------------------------------------------------
      double time;
      ulong steps;
      char filename[128];
      //-------------------------------------------------------------------
      sprintf(filename, "%s/%s.%s%.3d.dat", DATAFOLDER, file, BRUTE, filenum);
      readApproxAccel(&time, &steps, (int)Ntot, index, accBuf, filename);
      for(int ii = 0; ii < (int)Ntot; ii++)
	direct[index[ii]] = accBuf[ii];
      //-------------------------------------------------------------------
      sprintf(filename, "%s/%s.%s.%.3d.dat", DATAFOLDER, file, accfile, filenum);
      readApproxAccel(&time, &steps, (int)Ntot, index, accBuf, filename);
      for(int ii = 0; ii < (int)Ntot; ii++)
	octree[index[ii]] = accBuf[ii];
      //-------------------------------------------------------------------
#   if  defined(PRINT_VALUES) || defined(ERROR_DETECT)
      //-------------------------------------------------------------------
      int unit_read;
      static char errfile[128];
      sprintf(errfile, "%s/%s.%s%.3d.dat", LOGFOLDER, serfile, ERRFILE, filenum);
      fp = fopen(errfile, "w");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", errfile);      }
#ifdef  USE_HDF5_FORMAT
      readSnapshot(&unit_read, &time, &steps, Ntot, &hdf5_dat, file, filenum, hdf5type);
      for(int ii = 0; ii < (int)Ntot; ii++){
	//-----------------------------------------------------------------
	body[ii]. x  = hdf5_dat.pos[ii * 3];      body[ii]. y = hdf5_dat.pos[ii * 3 + 1];      body[ii].z   = hdf5_dat.pos[ii * 3 + 2];
	body[ii].vx  = hdf5_dat.vel[ii * 3];      body[ii].vy = hdf5_dat.vel[ii * 3 + 1];      body[ii].vz  = hdf5_dat.vel[ii * 3 + 2];
	body[ii].ax  = hdf5_dat.acc[ii * 3];      body[ii].ay = hdf5_dat.acc[ii * 3 + 1];      body[ii].az  = hdf5_dat.acc[ii * 3 + 2];
	body[ii].idx = hdf5_dat.idx[ii    ];      body[ii]. m = hdf5_dat.  m[ii        ];      body[ii].pot = hdf5_dat.pot[ii        ];
	//-----------------------------------------------------------------
      }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#else///USE_HDF5_FORMAT
      readSnapshot(&unit_read, &time, &steps, Ntot,  body, file, filenum);
#endif//USE_HDF5_FORMAT
      //-------------------------------------------------------------------
      if( unit_read != unit ){	__KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);      }
      //-------------------------------------------------------------------
      qsort(body, (int)Ntot, sizeof(nbody_particle), idxAscendingOrder);
      //-------------------------------------------------------------------
#endif//defined(PRINT_VALUES) || defined(ERROR_DETECT)
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* data analysis */
      //-------------------------------------------------------------------
      __NOTE__("rank = %d, filenum = %d\n", mpi.rank, filenum);
      //-------------------------------------------------------------------
      chkTreeError((int)Ntot, direct, octree, accErr, grvErr, potErr);
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* analyze cumulative distribution of errors in accelerations and potentials */
      //-------------------------------------------------------------------
      PLFLT inv = 1.0 / (PLFLT)Ntot;
      for(ulong i = 0; i < Ntot; i++){
	cdf[i] = inv * (PLFLT)(Ntot - i);
	pot[i] =       (PLFLT)potErr[i];
	grv[i] =       (PLFLT)grvErr[i];
	acc[i] =       (PLFLT)accErr[i];
      }
      //-------------------------------------------------------------------
    /* output CDF of relative error */
#ifdef  USE_HDF5_FORMAT
      sprintf(filename, "%s/%s.%s.%.3d.h5", DATAFOLDER, serfile, ERRFILE_CDF, filenum);
      hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      char groupname[16];
      sprintf(groupname, "relErr");
      hid_t group = H5Gcreate(target, groupname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      /* preparation for data compression */
      hid_t dataset;
      /* compression using szip */
      uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
      uint szip_pixels_per_block = 8;
      hsize_t cdims = 128 * szip_pixels_per_block;
      hsize_t dims = Ntot;
      hid_t dataspace = H5Screate_simple(1, &dims, NULL);
      hid_t property = H5Pcreate(H5P_DATASET_CREATE);
      chkHDF5err(H5Pset_chunk(property, 1, &cdims));
      chkHDF5err(H5Pset_szip(property, szip_options_mask, szip_pixels_per_block));
      /* write CDF */
      dataset = H5Dcreate(group, "cdf", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cdf));
      chkHDF5err(H5Dclose(dataset));
      /* write potential error */
      dataset = H5Dcreate(group, "pot", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pot));
      chkHDF5err(H5Dclose(dataset));
      /* write acceleration error as a scalar */
      dataset = H5Dcreate(group, "grv", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grv));
      chkHDF5err(H5Dclose(dataset));
      /* write acceleration error as a vector */
      dataset = H5Dcreate(group, "acc", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
      chkHDF5err(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, acc));
      chkHDF5err(H5Dclose(dataset));
      /* write attribute data */
      hsize_t attr_dims = 1;
      dataspace = H5Screate_simple(1, &attr_dims, NULL);
      hid_t attribute;
      /* write current time */
      attribute = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &time));
      chkHDF5err(H5Aclose(attribute));
      /* write # of steps */
      attribute = H5Acreate(group, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
      chkHDF5err(H5Aclose(attribute));
      /* write # of N-body particles */
      attribute = H5Acreate(group, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &Ntot));
      chkHDF5err(H5Aclose(attribute));
#ifdef  DOUBLE_PRECISION
      attribute = H5Acreate(group, macname, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &mac));
#else///DOUBLE_PRECISION
      attribute = H5Acreate(group, macname, H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      chkHDF5err(H5Awrite(attribute, H5T_NATIVE_FLOAT, &mac));
#endif//DOUBLE_PRECISION
      chkHDF5err(H5Aclose(attribute));
      /* close the dataspace */
      chkHDF5err(H5Sclose(dataspace));
      chkHDF5err(H5Gclose(group));
      chkHDF5err(H5Fclose(target));
#endif//USE_HDF5_FORMAT
      //-------------------------------------------------------------------
#ifdef LOGPLOT_HOR
      for(ulong i = 0; i < Ntot; i++){
	pot[i] = log10(pot[i]);
	grv[i] = log10(grv[i]);
	acc[i] = log10(acc[i]);
      }
#endif//LOGPLOT_HOR
      //-------------------------------------------------------------------
#ifdef LOGPLOT_VER
      for(ulong i = 0; i < Ntot; i++)
	cdf[i] = log10(cdf[i]);
#endif//LOGPLOT_VER
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      /* plot the cumulative distribution function */
      //-------------------------------------------------------------------
      PLplotPltRange range;
#ifdef  LOGPLOT_HOR
      range.xmin = log10(1.0e-8);
      range.xmax = log10(1.0e-0);
      range.xlog = LOGARITHMIC_PLOT;
#else//LOGPLOT_HOR
      range.xmin = 1.0e-8;
      range.xmax = 1.0e-0;
      range.xlog = LINEAR_PLOT;
#endif//LOGPLOT_HOR
#ifdef  LOGPLOT_VER
      range.ymin = cdf[tail];
      range.ymax = log10(1.0);
      range.ylog = LOGARITHMIC_PLOT;
#else///LOGPLOT_VER
      range.ymin = 0.0;
      /* range.ymax = 9.9999e-1; */
      range.ymax = 1.0;
      range.ylog = LINEAR_PLOT;
#endif//LOGPLOT_VER
      range.xgrd = range.ygrd = true;
      //-------------------------------------------------------------------
      plotCDF((int)Ntot, cdf, acc, grv, pot, range, serfile, ifile, argc, argv);
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      FILE *logfp;
      //-------------------------------------------------------------------
      sprintf(filename, "%s/%s.%s.%s.%.3d.txt", DATAFOLDER, file, macname, "acc", filenum);
      logfp = fopen(filename, "a");
      if( logfp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
      fprintf(logfp, "%.13e", mac);
      for(int ii = 0; ii < NSUMMARY; ii++)
	fprintf(logfp, "\t%e", accErr[(ulong)floor((double)(Ntot - 1) * (1.0 - percentage[ii]))]);
      fprintf(logfp, "\n");
      fclose(logfp);
      //-------------------------------------------------------------------
      sprintf(filename, "%s/%s.%s.%s.%.3d.txt", DATAFOLDER, file, macname, "grv", filenum);
      logfp = fopen(filename, "a");
      if( logfp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
      fprintf(logfp, "%.13e", mac);
      for(int ii = 0; ii < NSUMMARY; ii++)
	fprintf(logfp, "\t%e", grvErr[(ulong)floor((double)(Ntot - 1) * (1.0 - percentage[ii]))]);
      fprintf(logfp, "\n");
      fclose(logfp);
      //-------------------------------------------------------------------
      sprintf(filename, "%s/%s.%s.%s.%.3d.txt", DATAFOLDER, file, macname, "pot", filenum);
      logfp = fopen(filename, "a");
      if( logfp == NULL ){	__KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);      }
      fprintf(logfp, "%.13e", mac);
      for(int ii = 0; ii < NSUMMARY; ii++)
	fprintf(logfp, "\t%e", potErr[(ulong)floor((double)(Ntot - 1) * (1.0 - percentage[ii]))]);
      fprintf(logfp, "\n");
      fclose(logfp);
      //-------------------------------------------------------------------
      if( lstNum == 1 )
	for(int jj = 0; jj < NSUMMARY; jj++){
	  ulong idx = (ulong)floor((double)(Ntot - 1) * (1.0 - percentage[jj]));
	  potErrList[ifile * NSUMMARY + jj] = potErr[idx];
	  grvErrList[ifile * NSUMMARY + jj] = grvErr[idx];
	  accErrList[ifile * NSUMMARY + jj] = accErr[idx];
	}/* for(int jj = 0; jj < NSUMMARY; jj++){ */
      //-------------------------------------------------------------------

      //-------------------------------------------------------------------
      ifile += mpi.size;
      //-------------------------------------------------------------------
#   if  defined(PRINT_VALUES) || defined(ERROR_DETECT)
      fclose(fp);
#endif//defined(PRINT_VALUES) || defined(ERROR_DETECT)
      //-------------------------------------------------------------------
    }/* for(int filenum = mpi.rank; filenum < nfile; filenum += mpi.size){ */
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    if( lstNum == 1 ){
      //-------------------------------------------------------------------
      chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, accErrList, nfile * NSUMMARY, MPI_DOUBLE, MPI_SUM, mpi.comm));
      chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, grvErrList, nfile * NSUMMARY, MPI_DOUBLE, MPI_SUM, mpi.comm));
      chkMPIerr(MPI_Allreduce(MPI_IN_PLACE, potErrList, nfile * NSUMMARY, MPI_DOUBLE, MPI_SUM, mpi.comm));
      //-------------------------------------------------------------------
      if( mpi.rank == 0 ){
	//-----------------------------------------------------------------
	/* print out acceleration error */
	//-----------------------------------------------------------------
	printf("# relative error of acceleration as a vector\n");
	printf("#err");
	for(int ii = 0; ii < nfile; ii++)
	  printf("\t%d", ii);
	printf("\n");
	for(int ii = 0; ii < NSUMMARY; ii++){
	  //---------------------------------------------------------------
	  printf("%.8e", percentage[ii]);
	  for(int jj = 0; jj < nfile; jj++)
	    printf("\t%e", accErrList[jj * NSUMMARY + ii]);
	  printf("\n");
	  //---------------------------------------------------------------
	}/* for(int ii = 0; ii < NSUMMARY; ii++){ */
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* print out acceleration error */
	//-----------------------------------------------------------------
	printf("# relative error of acceleration as a scalar\n");
	printf("#err");
	for(int ii = 0; ii < nfile; ii++)
	  printf("\t%d", ii);
	printf("\n");
	for(int ii = 0; ii < NSUMMARY; ii++){
	  //---------------------------------------------------------------
	  printf("%.8e", percentage[ii]);
	  for(int jj = 0; jj < nfile; jj++)
	    printf("\t%e", grvErrList[jj * NSUMMARY + ii]);
	  printf("\n");
	  //---------------------------------------------------------------
	}/* for(int ii = 0; ii < NSUMMARY; ii++){ */
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	/* print out potential error */
	//-----------------------------------------------------------------
	printf("\n");
	printf("# relative error of potential\n");
	printf("#err");
	for(int ii = 0; ii < nfile; ii++)
	  printf("\t%d", ii);
	printf("\n");
	for(int ii = 0; ii < NSUMMARY; ii++){
	  //---------------------------------------------------------------
	  printf("%.8e", percentage[ii]);
	  for(int jj = 0; jj < nfile; jj++)
	    printf("\t%e", potErrList[jj * NSUMMARY + ii]);
	  printf("\n");
	  //---------------------------------------------------------------
	}/* for(int ii = 0; ii < NSUMMARY; ii++){ */
	//-----------------------------------------------------------------
      }/* if( mpi.rank == 0 ){ */
      //-------------------------------------------------------------------
      free(accErrList);
      free(grvErrList);
      free(potErrList);
      //-------------------------------------------------------------------
    }/* if( lstNum == 1 ){ */
    //---------------------------------------------------------------------
  }/* for(int lstIdx = 0; lstIdx < lstNum; lstIdx++){ */
  //-----------------------------------------------------------------------
  if( mpi.rank == 0 )
    fclose(lstfp);
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
#   if  defined(PRINT_VALUES) || defined(ERROR_DETECT)
  freeParticleDataAoS(body);
#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
  freeSnapshotArray(hdf5_pos, hdf5_vel, hdf5_acc, hdf5_m, hdf5_pot, hdf5_idx);
#endif//USE_HDF5_FORMAT
#endif//defined(PRINT_VALUES) || defined(ERROR_DETECT)
  //-----------------------------------------------------------------------
  free( index);
  free(accBuf);
  free(direct);
  free(octree);
  //-----------------------------------------------------------------------
  free(accErr);
  free(grvErr);
  free(potErr);
  //-----------------------------------------------------------------------
  free(cdf);
  free(acc);
  free(pot);
  //-----------------------------------------------------------------------
  freePercentile(percentage);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  exitMPI();
  //-----------------------------------------------------------------------
  return (0);
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void plotCDF
(int num, PLFLT *cdf, PLFLT *acc, PLFLT *grv, PLFLT *pot, PLplotPltRange box,
 char file[], const int filenum, int argc, char **argv)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* global setting(s) */
  //-----------------------------------------------------------------------
  /* maximum number of data kinds */
  const PLINT pkind = 0;
  const PLINT lkind = 3;
  //-----------------------------------------------------------------------
  const PLBOOL Euni = true;/* E stands for error */
  const PLBOOL Auni = true;/* A stands for acceleration */
  const PLBOOL Guni = true;/* G stants for gravity */
  const PLBOOL Puni = true;/* P stands for potential */
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  /* preparation for dataset */
  //-----------------------------------------------------------------------
  /* memory allocation for index */
  PLINT Ekind = (Euni) ? (IMAX(pkind, lkind)) : (pkind + lkind);
  PLINT *Enum;  allocPLINT(&Enum, Ekind);
  PLINT Akind = 1;  PLINT *Anum;  allocPLINT(&Anum, Akind);
  PLINT Gkind = 1;  PLINT *Gnum;  allocPLINT(&Gnum, Gkind);
  PLINT Pkind = 1;  PLINT *Pnum;  allocPLINT(&Pnum, Pkind);
  //-----------------------------------------------------------------------
  /* set number of data points */
  for(PLINT ii = 0; ii < Ekind; ii++)    Enum[ii] = num;
  for(PLINT ii = 0; ii < Akind; ii++)    Anum[ii] = num;
  for(PLINT ii = 0; ii < Gkind; ii++)    Gnum[ii] = num;
  for(PLINT ii = 0; ii < Pkind; ii++)    Pnum[ii] = num;
  //-----------------------------------------------------------------------
  /* memory allocation for data */
  PLFLT **hor;  allocPointer4PLFLT(&hor, Ekind);
  PLFLT **ver;  allocPointer4PLFLT(&ver, Ekind);
  /* data preparation */
  if( Euni ){
    hor[0] = acc;    ver[0] = cdf;
    hor[1] = grv;    ver[1] = cdf;
    hor[2] = pot;    ver[2] = cdf;
  }
  else{
    hor[0] = acc;    ver[0] = cdf;
    hor[1] = grv;    ver[1] = cdf;
    hor[2] = pot;    ver[2] = cdf;
    hor[3] = acc;    ver[3] = cdf;
    hor[4] = grv;    ver[4] = cdf;
    hor[5] = pot;    ver[5] = cdf;
  }
  //-----------------------------------------------------------------------
  /* set symbol and line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  for(PLINT ii = 0; ii < lkind; ii++)
    ls[ii].width = BOLD_LINE;
  //-----------------------------------------------------------------------
  /* set labels */
  char Elab[PLplotCharWords];  sprintf(Elab, "Relative error");
  char Alab[PLplotCharWords];  sprintf(Alab, "|#fi#<0x12>a#fr#di#u#utree#d - #fi#<0x12>a#fr#di#u#udirect#d| / #fia#fr#di#u#udirect#d");
  char Glab[PLplotCharWords];  sprintf(Glab, "|#fia#fr#di#u#utree#d / #fia#fr#di#u#udirect#d - 1|");
  char Plab[PLplotCharWords];  sprintf(Plab, "|#fi#gF#fr#di#u#utree#d / #fi#gF#fr#di#u#udirect#d -1|");
  char ylab[PLplotCharWords];  sprintf(ylab, "Cumulative distribution");
  //-----------------------------------------------------------------------
  /* set caption(s) */
  PLplotCaption Ecap;  setDefaultCaption(&Ecap);  Ecap.write = false;
  PLplotCaption Acap;  setDefaultCaption(&Acap);  Acap.write = true;
  PLplotCaption Gcap;  setDefaultCaption(&Gcap);  Gcap.write = true;
  PLplotCaption Pcap;  setDefaultCaption(&Pcap);  Pcap.write = true;
  sprintf(Ecap.side, "%s", "t");
  sprintf(Acap.side, "%s", "t");
  sprintf(Gcap.side, "%s", "t");
  sprintf(Pcap.side, "%s", "t");
  sprintf(Ecap.text, "%s", "Cumulative distribution function");
  sprintf(Acap.text, "%s", "(a) acc");
  sprintf(Gcap.text, "%s", "(b) grv");
  sprintf(Pcap.text, "%s", "(c) pot");
  //-----------------------------------------------------------------------
  /* set legends */
  PLplotLegend Eleg;  setDefaultLegend(&Eleg, false);  Eleg.write = true;
  PLplotLegend Aleg;  setDefaultLegend(&Aleg, false);  Aleg.write = false;
  PLplotLegend Gleg;  setDefaultLegend(&Gleg, false);  Gleg.write = false;
  PLplotLegend Pleg;  setDefaultLegend(&Pleg, false);  Pleg.write = false;
  char *ElegTxt;
  {
    allocChar4PLplot(&ElegTxt, Ekind);
    allocPointer4Char4PLplot(&(Eleg.text), Ekind);
    assignChar4PLplot(Ekind, Eleg.text, ElegTxt);
  }
  sprintf(Eleg.text[0], "%s", "#fi#<0x12>a#fr");
  sprintf(Eleg.text[1], "%s", "#fia#fr");
  sprintf(Eleg.text[2], "%s", "#gF");
  Eleg.pos = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;
  char *AlegTxt;
  {
    allocChar4PLplot(&AlegTxt, Akind);
    allocPointer4Char4PLplot(&(Aleg.text), Akind);
    assignChar4PLplot(Akind, Aleg.text, AlegTxt);
  }
  sprintf(Aleg.text[0], "%s", "acceleration");
  Aleg.pos = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;
  char *GlegTxt;
  {
    allocChar4PLplot(&GlegTxt, Gkind);
    allocPointer4Char4PLplot(&(Gleg.text), Gkind);
    assignChar4PLplot(Gkind, Gleg.text, GlegTxt);
  }
  sprintf(Gleg.text[0], "%s", "gravity");
  Gleg.pos = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;
  char *PlegTxt;
  {
    allocChar4PLplot(&PlegTxt, Pkind);
    allocPointer4Char4PLplot(&(Pleg.text), Pkind);
    assignChar4PLplot(Pkind, Pleg.text, PlegTxt);
  }
  sprintf(Pleg.text[0], "%s", "potential");
  Pleg.pos = PL_POSITION_LEFT | PL_POSITION_TOP | PL_POSITION_INSIDE;
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory allocation to exploit multiple-plot function */
  //-----------------------------------------------------------------------
  /* maximum number of panels in each direction */
  const PLINT nxpanel = 1;
  const PLINT nypanel = 3;
  //-----------------------------------------------------------------------
  /* specify plot or skip panel */
  PLBOOL *plot;  allocPLBOOL(&plot, nxpanel * nypanel + 1);
  /* arrays related to draw line(s) and point(s) */
  PLINT *nlkind;  allocPLINT(&nlkind, nxpanel * nypanel + 1);
  PLINT *npkind;  allocPLINT(&npkind, nxpanel * nypanel + 1);
  PLINT **lnum;  allocPointer4PLINT(&lnum, nxpanel * nypanel + 1);
  PLINT **pnum;  allocPointer4PLINT(&pnum, nxpanel * nypanel + 1);
  PLFLT ***lx;  allocDoublePointer4PLFLT(&lx, nxpanel * nypanel + 1);
  PLFLT ***ly;  allocDoublePointer4PLFLT(&ly, nxpanel * nypanel + 1);
  PLFLT ***px;  allocDoublePointer4PLFLT(&px, nxpanel * nypanel + 1);
  PLFLT ***py;  allocDoublePointer4PLFLT(&py, nxpanel * nypanel + 1);
  PLplotLineStyle ** line;  allocPointer4PLplotLineStyle(& line, nxpanel * nypanel + 1);
  PLplotPointType **point;  allocPointer4PLplotPointType(&point, nxpanel * nypanel + 1);
  /* arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel + 1);
  /* arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel + 1);
  /* arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel + 1);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel + 1);
  /* arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel + 1);
  /* arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel + 1);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel + 1);
  /* array to set figure name */
  char figfile[PLplotCharWords];
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel * nypanel + 1);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel * nypanel + 1);
  assignChar4PLplot(nxpanel * nypanel + 1, figname, _figname);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* configure to cumulative distribution function of relative error */
  //-----------------------------------------------------------------------
  /* common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    for(PLINT jj = 0; jj < nypanel; jj++){
      //-------------------------------------------------------------------
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
      //-------------------------------------------------------------------
      /* global setting(s) */
      plot[idx] = true;
      //-------------------------------------------------------------------
      /* line setting(s) */
      line[idx] = ls;
      ly  [idx] = ver;
      //-------------------------------------------------------------------
      /* point setting(s) */
      point[idx] = pt;
      py   [idx] = ver;
      //-------------------------------------------------------------------
      /* errorbar setting(s) */
      errbar[idx] = false;
      //-------------------------------------------------------------------
      /* plot area */
      range[idx] = box;
      //-------------------------------------------------------------------
      /* label setting(s) */
      ylabel[idx] = ylab;
      //-------------------------------------------------------------------
    }
  }
  //-----------------------------------------------------------------------
  /* setting(s) for cumulative distribution function of relative error in potential */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 0);
    //---------------------------------------------------------------------
    /* line setting(s) */
    nlkind[idx] = IMIN(Pkind, lkind);
    lnum  [idx] = Pnum;
    lx    [idx] = &pot;
    //---------------------------------------------------------------------
    /* point setting(s) */
    npkind[idx] = IMIN(Pkind, pkind);
    pnum  [idx] = Pnum;
    px    [idx] = &pot;
    //---------------------------------------------------------------------
    /* label setting(s) */
    xlabel[idx] = Plab;
    //---------------------------------------------------------------------
    /* miscs */
    cap[idx] = Pcap;
    leg[idx] = Pleg;
    uni[idx] = Puni;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  /* setting(s) for cumulative distribution function of relative error in gravity */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 1);
    //---------------------------------------------------------------------
    /* line setting(s) */
    nlkind[idx] = IMIN(Gkind, lkind);
    lnum  [idx] = Gnum;
    lx    [idx] = &grv;
    //---------------------------------------------------------------------
    /* point setting(s) */
    npkind[idx] = IMIN(Gkind, pkind);
    pnum  [idx] = Gnum;
    px    [idx] = &grv;
    //---------------------------------------------------------------------
    /* label setting(s) */
    xlabel[idx] = Glab;
    //---------------------------------------------------------------------
    /* miscs */
    cap[idx] = Gcap;
    leg[idx] = Gleg;
    uni[idx] = Guni;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  /* setting(s) for cumulative distribution function of relative error in acceleration */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 2);
    //---------------------------------------------------------------------
    /* line setting(s) */
    nlkind[idx] = IMIN(Akind, lkind);
    lnum  [idx] = Anum;
    lx    [idx] = &acc;
    //---------------------------------------------------------------------
    /* point setting(s) */
    npkind[idx] = IMIN(Akind, pkind);
    pnum  [idx] = Anum;
    px    [idx] = &acc;
    //---------------------------------------------------------------------
    /* label setting(s) */
    xlabel[idx] = Alab;
    //---------------------------------------------------------------------
    /* miscs */
    cap[idx] = Acap;
    leg[idx] = Aleg;
    uni[idx] = Auni;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  /* setting(s) for cumulative distribution function of relative error in acceleration, gravity and acceleration */
  for(PLINT ii = 0; ii < nxpanel; ii++){
    //---------------------------------------------------------------------
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, 3);
    //---------------------------------------------------------------------
    /* global setting(s) */
    plot[idx] = true;
    //---------------------------------------------------------------------
    /* line setting(s) */
    nlkind[idx] = IMIN(Ekind, lkind);
    line  [idx] = ls;
    lnum  [idx] = Enum;
    lx    [idx] = hor;
    ly    [idx] = ver;
    //---------------------------------------------------------------------
    /* point setting(s) */
    npkind[idx] = IMIN(Ekind, pkind);
    point [idx] = pt;
    pnum  [idx] = Enum;
    px    [idx] = hor;
    py    [idx] = ver;
    //---------------------------------------------------------------------
    /* errorbar setting(s) */
    errbar[idx] = false;
    //---------------------------------------------------------------------
    /* plot area */
    range[idx] = box;
    //---------------------------------------------------------------------
    /* label setting(s) */
    xlabel[idx] = Elab;
    ylabel[idx] = ylab;
    //---------------------------------------------------------------------
    /* miscs */
    cap[idx] = Ecap;
    leg[idx] = Eleg;
    uni[idx] = Euni;
    //---------------------------------------------------------------------
  }
  //-----------------------------------------------------------------------
  /* individual file names */
  sprintf(figfile                                 , "%s_%s_%.3d", file, "cdf", filenum);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 0)], "%s_%s_%.3d", file, "pot", filenum);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 1)], "%s_%s_%.3d", file, "grv", filenum);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 2)], "%s_%s_%.3d", file, "acc", filenum);
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 3)], "%s_%s_%.3d", file, "err", filenum);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* create figure(s) */
  //-----------------------------------------------------------------------
#ifndef PLOT_SUMMARY_FIGURE_ONLY
  for(PLINT idx = 0; idx < nxpanel * nypanel; idx++)
    plotData(1, 1, &plot[idx], false, false,
  	     &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
  	     &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
  	     &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
  	     &cap[idx], &leg[idx], &uni[idx], &range[idx],
  	     &xlabel[idx], &ylabel[idx], "", figname[idx], argc, argv);
#endif//PLOT_SUMMARY_FIGURE_ONLY
  //-----------------------------------------------------------------------
  plotData(1, 1, &plot[nxpanel * nypanel], false, false,
	   &nlkind[nxpanel * nypanel], & line[nxpanel * nypanel], &lnum[nxpanel * nypanel], &lx[nxpanel * nypanel], &ly[nxpanel * nypanel],
	   &npkind[nxpanel * nypanel], &point[nxpanel * nypanel], &pnum[nxpanel * nypanel], &px[nxpanel * nypanel], &py[nxpanel * nypanel],
	   &errbar[nxpanel * nypanel], NULL, NULL, NULL, NULL, NULL, NULL,
	   &cap[nxpanel * nypanel], &leg[nxpanel * nypanel], &uni[nxpanel * nypanel], &range[nxpanel * nypanel],
	   &xlabel[nxpanel * nypanel], &ylabel[nxpanel * nypanel], "", figname[nxpanel * nypanel], argc, argv);
  //-----------------------------------------------------------------------
#ifndef PLOT_SUMMARY_FIGURE_ONLY
  /* for(PLINT idx = 0; idx < nxpanel * nypanel; idx++) */
  /*   leg[idx].write = true; */
  for(PLINT idx = 0; idx < nxpanel * nypanel; idx++)
    xlabel[idx] = Elab;
  plotData(nxpanel, nypanel, plot, true, true,
  	   nlkind,  line, lnum, lx, ly,
  	   npkind, point, pnum, px, py,
  	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
  	   cap, leg, uni, range,
  	   xlabel, ylabel, "", figfile, argc, argv);
#endif//PLOT_SUMMARY_FIGURE_ONLY
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  /* memory deallocation */
  //-----------------------------------------------------------------------
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  free(figname);  free(_figname);
  //-----------------------------------------------------------------------
  free(Enum);  free(hor);  free(ver);
  free(Anum);
  free(Gnum);
  free(Pnum);
  free(ls);  free(pt);
  //-----------------------------------------------------------------------
  free(ElegTxt);
  free(AlegTxt);
  free(GlegTxt);
  free(PlegTxt);
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifdef __ICC
/* Disable ICC's remark #161: unrecognized #pragma */
#     pragma warning (disable:161)
#endif//__ICC
//-------------------------------------------------------------------------
int dblAscendingOrder(const void *a, const void *b)
{
  //-----------------------------------------------------------------------
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
  if(          (*((double *)a)) > (*((double *)b)) ){    return ( 1);  }
  else{    if( (*((double *)a)) < (*((double *)b)) ){    return (-1);  }
    else                                                 return ( 0);  }
#pragma GCC diagnostic pop
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
int dblDescendingOrder(const void *a, const void *b)
{
  //-----------------------------------------------------------------------
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
  if(          (*((double *)a)) < (*((double *)b)) ){    return ( 1);  }
  else{    if( (*((double *)a)) > (*((double *)b)) ){    return (-1);  }
    else                                                 return ( 0);  }
#pragma GCC diagnostic pop
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
#ifdef __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void chkTreeError(int num, acceleration * restrict direct, acceleration * restrict octree,
		  double * restrict accErr, double * restrict grvErr, double * restrict potErr)
{
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "start");
  //-----------------------------------------------------------------------
  for(int ii = 0; ii < num; ii++){
    //---------------------------------------------------------------------
    const double ax = (double)direct[ii].x;
    const double ay = (double)direct[ii].y;
    const double az = (double)direct[ii].z;
    const double a2 = 1.0e-30 + ax * ax + ay * ay + az * az;
    //---------------------------------------------------------------------
    const double tx = (double)octree[ii].x;
    const double ty = (double)octree[ii].y;
    const double tz = (double)octree[ii].z;
    //---------------------------------------------------------------------
    const double ex = tx - ax;
    const double ey = ty - ay;
    const double ez = tz - az;
    //---------------------------------------------------------------------
    accErr[ii] = 1.0e-30 +       sqrt((ex * ex + ey * ey + ez * ez) /       a2);
    grvErr[ii] = 1.0e-30 + fabs((sqrt (tx * tx + ty * ty + tz * tz) / (sqrt(a2))) - 1.0);
    potErr[ii] = 1.0e-30 + fabs((double)octree[ii].pot / (double)direct[ii].pot - 1.0);
    //---------------------------------------------------------------------
#   if  defined(PRINT_VALUES) || defined(ERROR_DETECT)
#ifdef  ERROR_DETECT
    if( (grvErr[ii] > 1.0) || (potErr[ii] > 1.0) )
#endif//ERROR_DETECT
      {
	fprintf(fp, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		ii, accErr[ii], potErr[ii],
		body  [ii].x,   body[ii].y,   body[ii].z,
		octree[ii].x, octree[ii].y, octree[ii].z, octree[ii].pot,
		direct[ii].x, direct[ii].y, direct[ii].z, direct[ii].pot);
      }
#endif//defined(PRINT_VALUES) || defined(ERROR_DETECT)
    //---------------------------------------------------------------------
  }/* for(int ii = 0; ii < num; ii++){ */
  //-----------------------------------------------------------------------
  qsort(accErr, num, sizeof(double), dblDescendingOrder);
  qsort(grvErr, num, sizeof(double), dblDescendingOrder);
  qsort(potErr, num, sizeof(double), dblDescendingOrder);
  //-----------------------------------------------------------------------
  __NOTE__("%s\n", "end");
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
