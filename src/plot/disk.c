/**
 * @file disk.c
 *
 * @brief Plot code of disk components in the initial condition
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2018/03/08 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "myutil.h"
#include "name.h"
#include "constants.h"
#include "plplotlib.h"

extern const double      length2astro;extern const char      length_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double        time2astro;extern const char        time_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double        mass2astro;extern const char        mass_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double    velocity2astro;extern const char    velocity_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double     density2astro;extern const char     density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double col_density2astro;extern const char col_density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];


typedef struct
{
  real *RR, *zz, *zd;
  real *rho, *Phi;
  real *Sigma;
  real *vcirc, *sigmaR, *sigmap, *sigmaz;
  real *kappa, *Omega;
  real *QQ, *lambda;
  double Rs, zs, Mtot;
  int Nhor, Nver;
} disk_data;


void plotDistributionMaps
(const int ndisk, disk_data *disk,
 PLplotPltRange zdbox, PLplotPltRange velbox, PLplotPltRange sigbox, PLplotPltRange freqbox, PLplotPltRange Qbox,
 char *file, int argc, char **argv);


int main(int argc, char **argv)
{
  /** initialization */
  if( argc < 3 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 3);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -problem=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 3 ){ */

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,    "file", &file));
  int  problem;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "problem", &problem));

  modifyArgcArgv4PLplot(&argc, argv, 3);


  /** read # of components */
  int kind = 1;
  int skind = 1;
  if( problem >= 1 ){
    FILE *fp;
    char filename[256];
    sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file);
    fp = fopen(filename, "r");
    if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

    int unit;
    fscanf(fp, "%d", &unit);
    setPhysicalConstantsAndUnitSystem(unit, 0);
    fscanf(fp, "%d\t%d", &kind, &skind);

    fclose(fp);
  }/* if( problem >= 1 ){ */

#if 0
  fprintf(stdout, "kind = %d, skind = %d\n", kind, skind);
  exit(0);
#endif


  /** memory allocation */
  const int Ndisk = kind - skind;
  disk_data *disk;
  disk = (disk_data *)malloc(sizeof(disk_data) * Ndisk);
  if( disk == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate disk.");  }

  real *RR, *zz, *zd;
  real *rho, *Phi;
  real *Sigma;
  real *vcirc, *sigmaR, *sigmap, *sigmaz;
  real *kappa, *Omega;
  real *QQ, *lambda;


#ifdef  USE_HDF5_FORMAT
#ifdef  DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_DOUBLE;
#else///DOUBLE_PRECISION
  hid_t hdf5_real = H5T_NATIVE_FLOAT;
#endif//DOUBLE_PRECISION
#endif//USE_HDF5_FORMAT

  /** load initial profile */
  if( problem >= 1 ){

#ifdef  USE_HDF5_FORMAT
    /** read HDF5 file */
    char filename[128];
    sprintf(filename, "%s/%s.%s.h5", DATAFOLDER, file, "disk");
    hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    /** error check for precision of floating-point numbers */
    int useDP;
    hid_t attribute = H5Aopen(target, "useDP", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
    chkHDF5err(H5Aclose(attribute));
#ifdef  DOUBLE_PRECISION
    if( useDP != 1 ){
      __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, true);
    }/* if( useDP != 1 ){ */
#else///DOUBLE_PRECISION
    if( useDP != 0 ){
      __KILL__(stderr, "ERROR: useDP (%d) differs with that in the code (%d)\n", useDP, false);
    }/* if( useDP != 0 ){ */
#endif//DOUBLE_PRECISION

    for(int ii = 0; ii < Ndisk; ii++){
      char grp[16];
      sprintf(grp, "data%d", ii);
      hid_t h5group = H5Gopen(target, grp, H5P_DEFAULT);

      /** read attributes */
      int dims[2];
      /** read # of arrays */
      attribute = H5Aopen(h5group, "num", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, dims));
      chkHDF5err(H5Aclose(attribute));
      disk[ii].Nhor = dims[0];
      disk[ii].Nver = dims[1];

      /** read scale radius and scale height */
      attribute = H5Aopen(h5group, "Rs", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &disk[ii].Rs));
      chkHDF5err(H5Aclose(attribute));
      attribute = H5Aopen(h5group, "zd", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &disk[ii].zs));
      chkHDF5err(H5Aclose(attribute));
      attribute = H5Aopen(h5group, "Mtot", H5P_DEFAULT);
      chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &disk[ii].Mtot));
      chkHDF5err(H5Aclose(attribute));

      /** memory allocation if necessary */
      if( ii == 0 ){
	RR = (real *)malloc(sizeof(real) * disk[ii].Nhor);	if( RR == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate RR.");	}
	zz = (real *)malloc(sizeof(real) * disk[ii].Nver);	if( zz == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate zz.");	}

	Sigma = (real *)malloc(sizeof(real) * disk[ii].Nhor * Ndisk);	if( Sigma == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate Sigma.");	}

	vcirc  = (real *)malloc(sizeof(real) * disk[ii].Nhor * Ndisk);	if( vcirc  == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate vcirc." );	}
	sigmaR = (real *)malloc(sizeof(real) * disk[ii].Nhor * Ndisk);	if( sigmaR == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate sigmaR.");	}
	sigmap = (real *)malloc(sizeof(real) * disk[ii].Nhor * Ndisk);	if( sigmap == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate sigmap.");	}
	sigmaz = (real *)malloc(sizeof(real) * disk[ii].Nhor * Ndisk);	if( sigmaz == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate sigmaz.");	}

	kappa = (real *)malloc(sizeof(real) * disk[ii].Nhor * Ndisk);	if( kappa == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate kappa.");	}
	Omega = (real *)malloc(sizeof(real) * disk[ii].Nhor * Ndisk);	if( Omega == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate Omega.");	}

	zd     = (real *)malloc(sizeof(real) * disk[ii].Nhor * Ndisk);	if(   zd   == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate zd."    );	}
	QQ     = (real *)malloc(sizeof(real) * disk[ii].Nhor * Ndisk);	if(   QQ   == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate QQ."    );	}
	lambda = (real *)malloc(sizeof(real) * disk[ii].Nhor * Ndisk);	if( lambda == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate lambda.");	}

	rho = (real *)malloc(sizeof(real) * disk[ii].Nhor * disk[ii].Nver * Ndisk);	if( rho == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate rho.");	}
	Phi = (real *)malloc(sizeof(real) * disk[ii].Nhor * disk[ii].Nver        );	if( Phi == NULL ){	  __KILL__(stderr, "%s\n", "ERROR: failure to allocate Phi.");	}
      }/* if( ii == 0 ){ */

      disk[ii].RR = RR;      disk[ii].zz = zz;
      disk[ii].Sigma = &Sigma[disk[ii].Nhor * ii];      disk[ii].vcirc = &vcirc[disk[ii].Nhor * ii];
      disk[ii].sigmaR = &sigmaR[disk[ii].Nhor * ii];      disk[ii].sigmap = &sigmap[disk[ii].Nhor * ii];      disk[ii].sigmaz = &sigmaz[disk[ii].Nhor * ii];
      disk[ii].kappa = &kappa[disk[ii].Nhor * ii];      disk[ii].Omega = &Omega[disk[ii].Nhor * ii];
      disk[ii].zd = &zd[disk[ii].Nhor * ii];      disk[ii].QQ = &QQ[disk[ii].Nhor * ii];      disk[ii].lambda = &lambda[disk[ii].Nhor * ii];
      disk[ii].rho = &rho[disk[ii].Nhor * disk[ii].Nver * ii];
      disk[ii].Phi = Phi;

      /** read disk data */
      hid_t dataset;

      /** common data for all components */
      if( ii == 0 ){
	/** read horizontal position */
	dataset = H5Dopen(h5group, "radius", H5P_DEFAULT);
	chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].RR));
	chkHDF5err(H5Dclose(dataset));
	/** read vertical position */
	dataset = H5Dopen(h5group, "height", H5P_DEFAULT);
	chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].zz));
	chkHDF5err(H5Dclose(dataset));
	/** read potential */
	dataset = H5Dopen(h5group, "Phi", H5P_DEFAULT);
	chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].Phi));
	chkHDF5err(H5Dclose(dataset));
      }/* if( ii == 0 ){ */

      /** individual data for each component */
      /** read column density profile */
      dataset = H5Dopen(h5group, "Sigma", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].Sigma));
      chkHDF5err(H5Dclose(dataset));
      /** read circular velocity profile */
      dataset = H5Dopen(h5group, "vcirc", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].vcirc));
      chkHDF5err(H5Dclose(dataset));
      /** read sigmaR profile */
      dataset = H5Dopen(h5group, "sigmaR", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].sigmaR));
      chkHDF5err(H5Dclose(dataset));
      /** read sigmap profile */
      dataset = H5Dopen(h5group, "sigmap", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].sigmap));
      chkHDF5err(H5Dclose(dataset));
      /** read sigmaz profile */
      dataset = H5Dopen(h5group, "sigmaz", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].sigmaz));
      chkHDF5err(H5Dclose(dataset));
      /** read kappa profile */
      dataset = H5Dopen(h5group, "kappa", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].kappa));
      chkHDF5err(H5Dclose(dataset));
      /** read Omega profile */
      dataset = H5Dopen(h5group, "Omega", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].Omega));
      chkHDF5err(H5Dclose(dataset));
      /** read scale height profile */
      dataset = H5Dopen(h5group, "zd", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].zd));
      chkHDF5err(H5Dclose(dataset));
      /** read critical wavelength profile */
      dataset = H5Dopen(h5group, "lambda", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].lambda));
      chkHDF5err(H5Dclose(dataset));
      /** read Toomre's Q profile */
      dataset = H5Dopen(h5group, "Toomre's Q", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].QQ));
      chkHDF5err(H5Dclose(dataset));
      /** read density distribution */
      dataset = H5Dopen(h5group, "rho", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, hdf5_real, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk[ii].rho));

      /** close the group */
      chkHDF5err(H5Gclose(h5group));
#if 0
      fprintf(stdout, "%d-th group read\n", ii);
      fflush(stdout);
#endif
    }/* for(int ii = 0; ii < Ndisk; ii++){ */

    /** close the file */
    chkHDF5err(H5Fclose(target));

#else///USE_HDF5_FORMAT

    for(int ii = 0; ii < Ndisk; ii++){
      sprintf(filename, "%s/%s.diskdat.%d.dat", DATAFOLDER, file, ii);
      fp = fopen(filename, "rb");
      if( fp == NULL ){	__KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);      }

      bool success = true;
      success &= (fread(&group[ii].disk_nrad, sizeof(   int), 1, fp) == 1);
      success &= (fread(&group[ii].disk_nazi, sizeof(   int), 1, fp) == 1);
      success &= (fread(&group[ii].disk_nsph, sizeof(   int), 1, fp) == 1);
      success &= (fread(&group[ii].disk_Rd  , sizeof(double), 1, fp) == 1);
      success &= (fread(&group[ii].disk_zd  , sizeof(double), 1, fp) == 1);
      group[ii].disk_radius  = (real *)malloc(group[ii].disk_nrad                       * sizeof(real));
      group[ii].disk_height  = (real *)malloc(                      group[ii].disk_nazi * sizeof(real));
      group[ii].disk_rho     = (real *)malloc(group[ii].disk_nrad * group[ii].disk_nazi * sizeof(real));
      group[ii].disk_pot     = (real *)malloc(group[ii].disk_nrad * group[ii].disk_nazi * sizeof(real));
      group[ii].disk_sig     = (real *)malloc(group[ii].disk_nrad                       * sizeof(real));
      group[ii].disk_Sigma   = (real *)malloc(group[ii].disk_nrad                       * sizeof(real));
      if( group[ii].disk_radius  == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_radius\n");	}
      if( group[ii].disk_height	 == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_height\n");	}
      if( group[ii].disk_rho	 == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_rho\n"   );	}
      if( group[ii].disk_pot	 == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_pot\n"   );	}
      if( group[ii].disk_sig	 == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_sig\n"   );	}
      if( group[ii].disk_Sigma	 == NULL ){	  __KILL__(stderr, "ERROR: failure to allocate group[ii].disk_Sigma\n" );	}

      success &= (fread(group[ii].disk_radius , sizeof(real), group[ii].disk_nrad                      , fp) == group[ii].disk_nrad                      );
      success &= (fread(group[ii].disk_height , sizeof(real),                       group[ii].disk_nazi, fp) ==                       group[ii].disk_nazi);
      success &= (fread(group[ii].disk_rho    , sizeof(real), group[ii].disk_nrad * group[ii].disk_nazi, fp) == group[ii].disk_nrad * group[ii].disk_nazi);
      success &= (fread(group[ii].disk_pot    , sizeof(real), group[ii].disk_nrad * group[ii].disk_nazi, fp) == group[ii].disk_nrad * group[ii].disk_nazi);
      success &= (fread(group[ii].disk_sig    , sizeof(real), group[ii].disk_nrad                      , fp) == group[ii].disk_nrad                      );
      success &= (fread(group[ii].disk_Sigma  , sizeof(real), group[ii].disk_nrad                      , fp) == group[ii].disk_nrad                      );

      if( !success ){	__KILL__(stderr, "ERROR: failure to read \"%s\"\n", filename);      }
      fclose(fp);
    }/* for(int ii = 0; ii < Ndisk; ii++){ */

#endif//USE_HDF5_FORMAT
  }/* if( problem >= 1 ){ */

#if 0
  fprintf(stdout, "HDF5 data was read.\n");
  exit(0);
#endif


  /** set plot range */
  PLFLT Rmin, Rmax, zmin, zmax;
  PLFLT velmin, velmax, sigmin, sigmax;
  PLFLT freqmin, freqmax, Qmin, Qmax;
  PLFLT rhomin, rhomax, Sigmin, Sigmax;
  switch( problem ){
  case 11:    /**< A galaxy with multiple components */
    rhomin = 2.0e-9;    Sigmin = 1.0e-7;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 1.0e+0;    Sigmax = 1.0e-1;    Rmax = 30.0;    zmax = 1.0;
    break;
  case 12:    /**< A galaxy with multiple components */
    rhomin = 2.0e-9;    Sigmin = 1.0e-7;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 1.0e+0;    Sigmax = 1.0e-1;    Rmax = 30.0;    zmax = 1.6;
    break;
  case 13:    /**< A galaxy with multiple components */
    rhomin = 2.0e-9;    Sigmin = 1.0e-6;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 1.0e+0;    Sigmax = 1.0e-0;    Rmax = 30.0;    zmax = 2.0;
    break;
  case 20:    /**< M31 model determined by Fardal et al. (2007) */
    rhomin = 1.0e-7;    Sigmin = 4.0e-6;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 2.0e+4;    Sigmax = 5.0e+3;    Rmax = 50.0;    zmax = 1.0;
    break;
  case 22:    /**< A trial multi components galaxy model */
    rhomin = 1.0e-6;    Sigmin = 5.0e-5;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 2.0e+3;    Sigmax = 5.0e+2;    Rmax = 50.0;    zmax = 2.0;
    break;
  case 23:    /**< MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008) */
    rhomin = 1.0e-6;    Sigmin = 5.0e-5;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 2.0e+3;    Sigmax = 5.0e+2;    Rmax = 50.0;    zmax = 2.0;
    break;
  case 24:    /**< M31 model (Sofue 2015; Gilbert et al. 2012) */
    rhomin = 1.0e-7;    Sigmin = 4.0e-6;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 2.0e+3;    Sigmax = 5.0e+1;    Rmax = 50.0;    zmax = 1.0;
    break;
  case 25:    /**< A trial multi components galaxy model */
    rhomin = 1.0e-6;    Sigmin = 5.0e-5;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 2.0e+3;    Sigmax = 5.0e+2;    Rmax = 50.0;    zmax = 2.0;
    break;
  case 26:    /**< A trial multi components galaxy model (spherical model) */
    rhomin = 1.0e-6;    Sigmin = 5.0e-5;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 2.0e+3;    Sigmax = 5.0e+2;    Rmax = 50.0;    zmax = 2.0;
    break;
  case 27:    /**< M31 model (NFW halo, de Vaucouleurs bulge, and exponential disk) */
    rhomin = 1.0e-7;    Sigmin = 4.0e-6;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 2.0e+3;    Sigmax = 5.0e+2;    Rmax = 50.0;    zmax = 1.0;
    break;
  case 28:    /**< A trial multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk) */
    rhomin = 1.0e-6;    Sigmin = 5.0e-5;    Rmin =   0.0;    zmin = 0.0;    velmin = 0.0e+0;    sigmin = 0.0e+0;    freqmin = 0.0e+0;    Qmin = 1.0e-2;
    rhomax = 2.0e+3;    Sigmax = 5.0e+2;    Rmax = 100.0;    zmax = 1.2;    velmax = 3.5e+1;    sigmax = 1.5e+1;    freqmax = 1.2e+2;    Qmax = 1.0e+2;
    break;
  case 29:    /**< Multi components galaxy model by Vasiliev & Athanassoula (2015) */
    rhomin = 1.0e-7;    Sigmin = 1.0e-5;    Rmin = 0.0;    zmin = 0.0;
    rhomax = 1.0e+2;    Sigmax = 4.0e-0;    Rmax = 8.0;    zmax = 0.2;
    break;
  case 30:    /**< time evolution of MW/A defined in Kuijken & Dubinski (1995) */
    rhomin = 1.0e-5;    Sigmin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    rhomax = 1.0e+4;    Sigmax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 31:    /**< time evolution of MW/B defined in Kuijken & Dubinski (1995) */
    rhomin = 1.0e-5;    Sigmin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    rhomax = 1.0e+4;    Sigmax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 32:    /**< time evolution of MW/C defined in Kuijken & Dubinski (1995) */
    rhomin = 1.0e-5;    Sigmin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    rhomax = 1.0e+3;    Sigmax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 33:    /**< time evolution of MW/D defined in Kuijken & Dubinski (1995) */
    rhomin = 1.0e-6;    Sigmin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    rhomax = 1.0e+3;    Sigmax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 34:    /**< time evolution of M31/A defined in Widrow et al. (2003) */
    rhomin = 1.0e-5;    Sigmin = 1.0e-4;    Rmin = 0.0;    zmin = 0.00;
    rhomax = 1.0e+4;    Sigmax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 35:    /**< time evolution of M31/D defined in Widrow et al. (2003) */
    rhomin = 1.0e-5;    Sigmin = 1.0e-4;    Rmin = 0.0;    zmin = 0.0;
    rhomax = 1.0e+4;    Sigmax = 1.0e+3;    Rmax = 2.0;    zmax = 0.2;
    break;
  case 36:    /**< time evolution of MWa defined in Widrow & Dubinski (2005) */
    rhomin = 1.0e-6;    Sigmin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    rhomax = 1.0e+5;    Sigmax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 37:    /**< time evolution of MWb defined in Widrow & Dubinski (2005) */
    rhomin = 1.0e-6;    Sigmin = 1.0e-5;    Rmin = 0.0;    zmin = 0.00;
    rhomax = 1.0e+5;    Sigmax = 1.0e+3;    Rmax = 2.0;    zmax = 0.08;
    break;
  case 38:    /**< MW like galaxy (based on Widrow et al. 2008, Bedorf et al. 2014) */
    rhomin = 1.0e-6;    Sigmin = 1.0e-5;    Rmin =  0.0;    zmin = 0.0;
    rhomax = 1.0e+4;    Sigmax = 1.0e+3;    Rmax = 28.0;    zmax = 0.8;
    break;
  default:
    rhomin = 1.0e-6;    Sigmin = 1.0e-6;
    rhomax = 2.0e+2;    Sigmax = 2.0e+2;
    break;
  }

  PLplotPltRange zdbox, velbox, sigbox, freqbox, Qbox;
  zdbox.xmin = velbox.xmin = sigbox.xmin = freqbox.xmin = Qbox.xmin = log10((double)disk[0].RR[0]);
  zdbox.xmax = velbox.xmax = sigbox.xmax = freqbox.xmax = Qbox.xmax = log10(Rmax);
  zdbox.ymin = zmin * (PLFLT)length2astro;  velbox.ymin = velmin * (PLFLT)velocity2astro;  sigbox.ymin = sigmin * (PLFLT)velocity2astro;
  zdbox.ymax = zmax * (PLFLT)length2astro;  velbox.ymax = velmax * (PLFLT)velocity2astro;  sigbox.ymax = sigmax * (PLFLT)velocity2astro;
  freqbox.ymin = freqmin / (PLFLT)time2astro;  Qbox.ymin = log10(Qmin);
  freqbox.ymax = freqmax / (PLFLT)time2astro;  Qbox.ymax = log10(Qmax);
  zdbox.xgrd = velbox.xgrd = sigbox.xgrd = freqbox.xgrd = Qbox.xgrd = true;  zdbox.xlog = velbox.xlog = sigbox.xlog = freqbox.xlog = Qbox.xlog = LOGARITHMIC_PLOT;
  zdbox.ygrd = velbox.ygrd = sigbox.ygrd = freqbox.ygrd = Qbox.ygrd = true;  zdbox.ylog = velbox.ylog = sigbox.ylog = freqbox.ylog = LINEAR_PLOT;  Qbox.ylog = LOGARITHMIC_PLOT;

#if 0
  for(int ii = 0; ii < Ndisk; ii++)
    fprintf(stdout, "%d-th disk: Nhor = %d, Nver = %d\n", ii, disk[ii].Nhor, disk[ii].Nver);
  exit(0);
#endif

  plotDistributionMaps(Ndisk, disk, zdbox, velbox, sigbox, freqbox, Qbox, file, argc, argv);


#ifdef  USE_HDF5_FORMAT
  free(RR);  free(zz);  free(Sigma);
  free(vcirc);  free(sigmaR);  free(sigmap);  free(sigmaz);
  free(kappa);  free(Omega);
  free(zd);  free(QQ);  free(lambda);
  free(rho);  free(Phi);
#else///USE_HDF5_FORMAT
#endif//USE_HDF5_FORMAT
  free(disk);


  return (0);
}


void plotDistributionMaps
(const int ndisk, disk_data *disk,
 PLplotPltRange zdbox, PLplotPltRange velbox, PLplotPltRange sigbox, PLplotPltRange freqbox, PLplotPltRange Qbox,
 char *file, int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /** global setting(s) */
  /** maximum number of data kinds */
  const PLINT nkind = 9;/**< zs, zd, vcirc, sigmaR, sigmap, sigmaz, kappa, Omega, Q */
  const PLINT pkind = 0;
  const PLINT lkind = ndisk * nkind;


  /** preparation for dataset */
  /** memory allocation for index */
  PLINT *num;  allocPLINT(&num, lkind);
  for(int ii = 0; ii < ndisk; ii++)
    for(int jj = 0; jj < nkind; jj++)
      num[INDEX2D(ndisk, nkind, ii, jj)] = disk[ii].Nhor;

  /** memory allocation for data */
  PLFLT *_xx, *_yy;
  {
    PLINT tot = 0;
    for(PLINT ii = 0; ii < lkind; ii++)      tot += num[ii];
    allocPLFLT(&_xx, tot);
    allocPLFLT(&_yy, tot);
  }
  PLFLT **xx;  allocPointer4PLFLT(&xx, lkind);
  PLFLT **yy;  allocPointer4PLFLT(&yy, lkind);
  xx[0] = _xx;
  yy[0] = _yy;
  for(PLINT ii = 1; ii < lkind; ii++){
    xx[ii] = xx[ii - 1] + num[ii - 1];
    yy[ii] = yy[ii - 1] + num[ii - 1];
  }/* for(PLINT ii = 1; ii < lkind; ii++){ */

  /** data preparation */
  const double freq2astro = 1.0 / time2astro;
  for(int ii = 0; ii < ndisk; ii++){
    const int ndat = disk[ii].Nhor;
    for(int jj = 0; jj < nkind; jj++){
      const int idx = INDEX2D(nkind, ndisk, jj, ii);
      for(int kk = 0; kk < ndat; kk++)
	xx[idx][kk] = (PLFLT)log10((double)disk[ii].RR[kk] * length2astro);
    }/* for(int jj = 0; jj < nkind; jj++){ */

    for(int kk = 0; kk < ndat; kk++){
      yy[INDEX2D(nkind, ndisk, 0, ii)][kk] = (PLFLT)LOG10(disk[ii].QQ[kk]);
      yy[INDEX2D(nkind, ndisk, 1, ii)][kk] = (PLFLT)disk[ii].Omega[kk];
      yy[INDEX2D(nkind, ndisk, 2, ii)][kk] = (PLFLT)disk[ii].kappa[kk];
      yy[INDEX2D(nkind, ndisk, 3, ii)][kk] = (PLFLT)disk[ii].zd[kk];
      yy[INDEX2D(nkind, ndisk, 4, ii)][kk] = (PLFLT)((double)disk[ii].zs * length2astro);
      yy[INDEX2D(nkind, ndisk, 5, ii)][kk] = (PLFLT)disk[ii].vcirc [kk];
      yy[INDEX2D(nkind, ndisk, 6, ii)][kk] = (PLFLT)disk[ii].sigmaR[kk];
      yy[INDEX2D(nkind, ndisk, 7, ii)][kk] = (PLFLT)disk[ii].sigmap[kk];
      yy[INDEX2D(nkind, ndisk, 8, ii)][kk] = (PLFLT)disk[ii].sigmaz[kk];
    }/* for(int kk = 0; kk < ndat; kk++){ */
  }/* for(int ii = 0; ii < ndisk; ii++){ */

  /** set line style */
  PLplotLineStyle *ls;  setDefaultLineStyle(lkind, &ls);
  for(int ii = 0; ii < lkind; ii++)
    ls[ii].width = MIDDLE_LINE;
  for(int ii = 0; ii < nkind; ii++)
    for(int jj = 0; jj < ndisk; jj++)
      ls[INDEX2D(nkind, ndisk, ii, jj)].color = PLplotColor[jj];
  /** Toomre's Q */
  ls[INDEX2D(nkind, ndisk, 0, 0)].style =  SOLID_LINE;
  ls[INDEX2D(nkind, ndisk, 0, 1)].style = DASHED_LINE;
  /** epicyclic frequency */
  ls[INDEX2D(nkind, ndisk, 1, 0)].style = DOTTED_LINE;
  ls[INDEX2D(nkind, ndisk, 1, 1)].style =  SOLID_LINE;
  ls[INDEX2D(nkind, ndisk, 2, 0)].style = TRIPLE_DOT_DASHED_LINE;
  ls[INDEX2D(nkind, ndisk, 2, 1)].style = DASHED_LINE;
  /** scale height */
  ls[INDEX2D(nkind, ndisk, 3, 0)].style =  SOLID_LINE;
  ls[INDEX2D(nkind, ndisk, 4, 0)].style = DOTTED_LINE;
  ls[INDEX2D(nkind, ndisk, 3, 1)].style =  SOLID_LINE;
  ls[INDEX2D(nkind, ndisk, 4, 1)].style = DOTTED_LINE;
  /** characteristic velocity */
  for(int ii = 0; ii < ndisk; ii++){
    ls[INDEX2D(nkind, ndisk, 5, ii)].style = TRIPLE_DOT_DASHED_LINE;
    ls[INDEX2D(nkind, ndisk, 6, ii)].style = DASHED_LINE;
    ls[INDEX2D(nkind, ndisk, 7, ii)].style =  SOLID_LINE;
    ls[INDEX2D(nkind, ndisk, 8, ii)].style = DOTTED_LINE;
  }/* for(int ii = 0; ii < ndisk; ii++){ */

  /** set labels */
  char xlab[PLplotCharWords];  sprintf(xlab, "#fiR #fr(%s)", length_astro_unit_name4plot);
  char zlab[PLplotCharWords];  sprintf(zlab, "#fiz #fr(%s)", length_astro_unit_name4plot);
  char flab[PLplotCharWords];  sprintf(flab, "#fif #fr(%s#u-1#d)", time_astro_unit_name4plot);
  char vlab[PLplotCharWords];  sprintf(vlab, "#fiv #fr(%s)", velocity_astro_unit_name4plot);
  char Qlab[PLplotCharWords];  sprintf(Qlab, "#fiQ#fr");

  /** set caption(s) */
  PLplotCaption zcap;  setDefaultCaption(&zcap);  zcap.write = true;  sprintf(zcap.side, "t");
  PLplotCaption fcap;  setDefaultCaption(&fcap);  fcap.write = true;  sprintf(fcap.side, "t");
  PLplotCaption vcap;  setDefaultCaption(&vcap);  vcap.write = true;  sprintf(vcap.side, "t");
  PLplotCaption Qcap;  setDefaultCaption(&Qcap);  Qcap.write = true;  sprintf(Qcap.side, "t");
  sprintf(zcap.text, "scale height");
  sprintf(fcap.text, "epicyclic frequency");
  sprintf(vcap.text, "typical velocity");
  sprintf(Qcap.text, "Toomre's #fiQ#fr");

  /** set legends */
  PLplotLegend zleg;  setDefaultLegend(&zleg, false);  zleg.write = true;
  PLplotLegend fleg;  setDefaultLegend(&fleg, false);  fleg.write = true;
  PLplotLegend vleg;  setDefaultLegend(&vleg, false);  vleg.write = true;
  PLplotLegend sleg;  setDefaultLegend(&sleg, false);  sleg.write = true;
  PLplotLegend Qleg;  setDefaultLegend(&Qleg, false);  Qleg.write = true;
  char *zlegTxt;
  {
    allocChar4PLplot(&zlegTxt, ndisk * 2);
    allocPointer4Char4PLplot(&(zleg.text), ndisk * 2);
    assignChar4PLplot(ndisk * 2, zleg.text, zlegTxt);
  }
  for(int ii = 0; ii < ndisk; ii++){
    sprintf(zleg.text[INDEX2D(2, ndisk, 0, ii)], "#fiz#fr#ds,output#u (disk%d)", ii);
    sprintf(zleg.text[INDEX2D(2, ndisk, 1, ii)], "#fiz#fr#ds,input #u (disk%d)", ii);
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  char *flegTxt;
  {
    allocChar4PLplot(&flegTxt, ndisk * 2);
    allocPointer4Char4PLplot(&(fleg.text), ndisk * 2);
    assignChar4PLplot(ndisk * 2, fleg.text, flegTxt);
  }
  for(int ii = 0; ii < ndisk; ii++){
    sprintf(fleg.text[INDEX2D(2, ndisk, 0, ii)], "#fi#gW#fr (disk%d)", ii);
    sprintf(fleg.text[INDEX2D(2, ndisk, 1, ii)], "#fi#gk#fr (disk%d)", ii);
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  char *vlegTxt;
  {
    allocChar4PLplot(&vlegTxt, ndisk * 4);
    allocPointer4Char4PLplot(&(vleg.text), ndisk * 4);
    assignChar4PLplot(ndisk * 4, vleg.text, vlegTxt);
  }
  char *slegTxt;
  {
    allocChar4PLplot(&slegTxt, ndisk * 3);
    allocPointer4Char4PLplot(&(sleg.text), ndisk * 3);
    assignChar4PLplot(ndisk * 3, sleg.text, slegTxt);
  }
  for(int ii = 0; ii < ndisk; ii++){
    sprintf(vleg.text[INDEX2D(4, ndisk, 0, ii)], "#fiv#fr#dc#u (disk%d)", ii);
    sprintf(vleg.text[INDEX2D(4, ndisk, 1, ii)], "#fi#gs#dR#u#fr (disk%d)", ii);
    sprintf(sleg.text[INDEX2D(3, ndisk, 0, ii)], "#fi#gs#dR#u#fr (disk%d)", ii);
    sprintf(vleg.text[INDEX2D(4, ndisk, 2, ii)], "#fi#gs#dp#u#fr (disk%d)", ii);
    sprintf(sleg.text[INDEX2D(3, ndisk, 1, ii)], "#fi#gs#dp#u#fr (disk%d)", ii);
    sprintf(vleg.text[INDEX2D(4, ndisk, 3, ii)], "#fi#gs#dz#u#fr (disk%d)", ii);
    sprintf(sleg.text[INDEX2D(3, ndisk, 2, ii)], "#fi#gs#dz#u#fr (disk%d)", ii);
  }/* for(int ii = 0; ii < ndisk; ii++){ */
  char *QlegTxt;
  {
    allocChar4PLplot(&QlegTxt, ndisk);
    allocPointer4Char4PLplot(&(Qleg.text), ndisk);
    assignChar4PLplot(ndisk, Qleg.text, QlegTxt);
  }
  for(int ii = 0; ii < ndisk; ii++)
    sprintf(Qleg.text[ii], "disk%d", ii);

  PLBOOL zuni = false;
  PLBOOL funi = false;
  PLBOOL vuni = false;
  PLBOOL suni = false;
  PLBOOL Quni = false;

  char title[PLplotCharWords];
  sprintf(title, "profile on the disk plane");


  /** memory allocation to exploit multiple-plot function */
  /** maximum number of panels in each direction */
  const PLINT nxpanel = 3;
  const PLINT nypanel = 2;

  /** specify plot or skip panel */
  PLBOOL *plot;
  allocPLBOOL(&plot, nxpanel * nypanel);
  /** arrays related to draw line(s) and point(s) */
  PLINT *nlkind;  allocPLINT(&nlkind, nxpanel * nypanel);
  PLINT *npkind;  allocPLINT(&npkind, nxpanel * nypanel);
  PLINT **lnum;  allocPointer4PLINT(&lnum, nxpanel * nypanel);
  PLINT **pnum;  allocPointer4PLINT(&pnum, nxpanel * nypanel);
  PLFLT ***lx;  allocDoublePointer4PLFLT(&lx, nxpanel * nypanel);
  PLFLT ***ly;  allocDoublePointer4PLFLT(&ly, nxpanel * nypanel);
  PLFLT ***px;  allocDoublePointer4PLFLT(&px, nxpanel * nypanel);
  PLFLT ***py;  allocDoublePointer4PLFLT(&py, nxpanel * nypanel);
  PLplotLineStyle ** line;  allocPointer4PLplotLineStyle(& line, nxpanel * nypanel);
  PLplotPointType **point;  allocPointer4PLplotPointType(&point, nxpanel * nypanel);
  /** arrays related to errorbar(s) */
  PLINT *errbar;  allocPLINT(&errbar, nxpanel * nypanel);
  /** arrays for caption(s) */
  PLplotCaption *cap;  allocPLplotCaption(&cap, nxpanel * nypanel);
  /** arrays for legend(s) */
  PLplotLegend *leg;  allocPLplotLegend(&leg, nxpanel * nypanel);
  PLBOOL       *uni;  allocPLBOOL      (&uni, nxpanel * nypanel);
  /** arrays for plot range */
  PLplotPltRange *range;  allocPLplotPltRange(&range, nxpanel * nypanel);
  /** arrays related to axis labels */
  char **xlabel;  allocPointer4Char4PLplot(&xlabel, nxpanel * nypanel);
  char **ylabel;  allocPointer4Char4PLplot(&ylabel, nxpanel * nypanel);
  /** array to set figure name */
  char figfile[PLplotCharWords];
  char *_figname;  allocChar4PLplot        (&_figname, nxpanel * nypanel);
  char **figname;  allocPointer4Char4PLplot(& figname, nxpanel * nypanel);
  assignChar4PLplot(nxpanel * nypanel, figname, _figname);


  /** configure to time evolution of energy */
  /** common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

      /** global setting(s) */
      plot[idx] = false;

      /** line setting(s) */
      nlkind[idx] = 0;
      lnum  [idx] = NULL;
      line  [idx] = NULL;
      lx    [idx] = NULL;
      ly    [idx] = NULL;

      /** point setting(s) */
      npkind[idx] = pkind;
      point [idx] = NULL;
      pnum  [idx] = NULL;
      px    [idx] = NULL;
      py    [idx] = NULL;

      /** errorbar setting(s) */
      errbar[idx] = 0;

      /** label setting(s) */
      xlabel[idx] = xlab;
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */
#if 0
  fprintf(stdout, "common settings completed\n");
  fflush(stdout);
#endif

  /** individual setting(s) */
  for(PLINT jj = 0; jj < nypanel; jj++){
    const int ii = 0;
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);
    plot[idx] = true;

    /** line setting(s) */
    const int lidx = INDEX2D(nkind, ndisk, (jj == 0) ? 0 : 3, 0);
    nlkind[idx] = ndisk * (1 << jj);
    lnum  [idx] = &num[lidx];
    line  [idx] = &ls [lidx];
    lx    [idx] = &xx [lidx];
    ly    [idx] = &yy [lidx];

    /** plot area */
    range[idx] = (jj == 0) ? Qbox : zdbox;

    /** label setting(s) */
    ylabel[idx] = (jj == 0) ? Qlab : zlab;

    /** captions */
    cap[idx] = (jj == 0) ? Qcap : zcap;
    sprintf(cap[idx].text, "(%c)", 97 + INDEX2D(nxpanel, nypanel, ii, nypanel - 1 - jj));

    /** legends */
    leg[idx] = (jj == 0) ? Qleg : zleg;
    leg[idx].pos = ((jj == 0) ? PL_POSITION_RIGHT : PL_POSITION_LEFT) | PL_POSITION_BOTTOM | PL_POSITION_INSIDE;
    uni[idx] = (jj == 0) ? Quni : zuni;
  }/* for(PLINT jj = 0; jj < nypanel; jj++){ */
  for(PLINT jj = 0; jj < nypanel; jj++){
    const int ii = 1;
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

    plot[idx] = true;

    /** line setting(s) */
    const int lidx = INDEX2D(nkind, ndisk, 5 + (jj == 0), 0);
    nlkind[idx] = ndisk * (3 + jj);
    lnum  [idx] = &num[lidx];
    line  [idx] = &ls [lidx];
    lx    [idx] = &xx [lidx];
    ly    [idx] = &yy [lidx];

    /** plot area */
    range[idx] = (jj == 0) ? sigbox : velbox;

    /** label setting(s) */
    ylabel[idx] = vlab;

    /** captions */
    cap[idx] = vcap;
    sprintf(cap[idx].text, "(%c)", 97 + INDEX2D(nxpanel, nypanel, ii, nypanel - 1 - jj));

    /** legends */
    leg[idx] = (jj == 0) ? sleg : vleg;
    leg[idx].pos = PL_POSITION_RIGHT | PL_POSITION_TOP | ((jj == 0) ? PL_POSITION_INSIDE : PL_POSITION_OUTSIDE);
    uni[idx] = (jj == 0) ? suni : vuni;
  }/* for(PLINT jj = 0; jj < nypanel; jj++){ */
  {
    const int ii = 2;
    const int jj = 0;
    const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

    plot[idx] = true;

    /** line setting(s) */
    const int lidx = INDEX2D(nkind, ndisk, 1, 0);
    nlkind[idx] = ndisk * 2;
    lnum  [idx] = &num[lidx];
    line  [idx] = &ls [lidx];
    lx    [idx] = &xx [lidx];
    ly    [idx] = &yy [lidx];

    /** plot area */
    range[idx] = freqbox;

    /** label setting(s) */
    ylabel[idx] = flab;

    /** captions */
    cap[idx] = fcap;
    sprintf(cap[idx].text, "(%c)", 97 + idx);

    /** legends */
    leg[idx] = fleg;
    leg[idx].pos = PL_POSITION_RIGHT | PL_POSITION_TOP | PL_POSITION_INSIDE;
    uni[idx] = funi;
  }
#if 0
  fprintf(stdout, "individual settings completed\n");
  fflush(stdout);
#endif

  /** individual file names */
  sprintf(figfile,                                     "%s_%s", file, "disk");
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 0)], "%s_%s_%s", file, "disk", "Q");
  sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 1)], "%s_%s_%s", file, "disk", "z");
  sprintf(figname[INDEX2D(nxpanel, nypanel, 1, 0)], "%s_%s_%s", file, "disk", "sig");
  sprintf(figname[INDEX2D(nxpanel, nypanel, 1, 1)], "%s_%s_%s", file, "disk", "vel");
  sprintf(figname[INDEX2D(nxpanel, nypanel, 2, 0)], "%s_%s_%s", file, "disk", "epi");


  /** create figure(s) */
#if 0
  fprintf(stdout, "ready to plot\n");
  exit(0);
#endif
  plotData(nxpanel, nypanel, plot, true, false,
  	   nlkind,  line, lnum, lx, ly,
  	   npkind, point, pnum, px, py,
  	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
  	   cap, leg, uni, range,
  	   xlabel, ylabel, "", figfile, argc, argv);

  for(int ii = 0; ii < nxpanel; ii++)
    for(int jj = 0; jj < nypanel; jj++){
      const int idx = INDEX2D(nxpanel, nypanel, ii, jj);

      if( plot[idx] ){
	cap[idx].write = false;
	plotData(1, 1, &plot[idx], false, false,
		 &nlkind[idx], & line[idx], &lnum[idx], &lx[idx], &ly[idx],
		 &npkind[idx], &point[idx], &pnum[idx], &px[idx], &py[idx],
		 &errbar[idx], NULL, NULL, NULL, NULL, NULL, NULL,
		 &cap[idx], &leg[idx], &uni[idx], &range[idx],
		 &xlabel[idx], &ylabel[idx], "", figname[idx], argc, argv);
      }/* if( plot[idx] ){ */
    }/* for(int jj = 0; jj < nypanel; jj++){ */


  /** memory deallocation */
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  free(figname);  free(_figname);

  free(_xx);  free(xx);
  free(_yy);  free(yy);
  free(ls);

  free(zlegTxt);
  free(flegTxt);
  free(vlegTxt);
  free(slegTxt);
  free(QlegTxt);

  free(num);


  __NOTE__("%s\n", "end");
}
