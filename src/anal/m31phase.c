/**
 * @file m31phase.c
 *
 * @brief Analyzer for phase-space representation in the M31-rest frame
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2019/09/05 (Thu)
 *
 * Copyright (C) 2017 Yohei Miki
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

/**
 * @def USE_SZIP_COMPRESSION
 *
 * @brief On to enable Szip compression for HDF5 files (default is ON).
 */
/* #define USE_SZIP_COMPRESSION */

/**
 * @def USE_GZIP_COMPRESSION
 *
 * @brief On to enable gzip compression for HDF5 files (default is ON).
 *
 * @detail currently, h5py does not accept Szip compression in default.
 */
#define USE_GZIP_COMPRESSION

#   if  defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)
#define USE_FILE_COMPRESSION
#endif//defined(USE_SZIP_COMPRESSION) || defined(USE_GZIP_COMPRESSION)

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include <unistd.h>
#include <sys/stat.h>

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
#include "myutil.h"
#include "name.h"
#include "constants.h"
#include "mpilib.h"

#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../file/io.h"


extern const double      length2astro;extern const char      length_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double        time2astro;extern const char        time_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double        mass2astro;extern const char        mass_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double    velocity2astro;extern const char    velocity_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double     density2astro;extern const char     density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double col_density2astro;extern const char col_density_astro_unit_name[CONSTANTS_H_CHAR_WORDS];
extern const double     senergy2astro;extern const char     senergy_astro_unit_name[CONSTANTS_H_CHAR_WORDS];


typedef struct
{
  ulong idx;
  real  x,  y,  z;
  real vx, vy, vz;
  real ax, ay, az;
  real m, pot;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  real ax_ext, ay_ext, az_ext, pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  real rad, hor;
} nbody_particle;


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
  if(          ((const nbody_particle *)a)->idx > ((const nbody_particle *)b)->idx ){    return ( 1);  }
  else{    if( ((const nbody_particle *)a)->idx < ((const nbody_particle *)b)->idx ){    return (-1);  }
    else{                                                                                return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}

#ifdef  __ICC
/* Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC


void generateVelocityDistribution(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body, real * restrict vr, real * restrict vt, const int nvr, const real vrmin, const real dvr, const int nvt, const real vtmin, const real dvt, real * restrict vrvt);
void generatePhaseSpaceDistribution(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body, real * restrict Etot, real * restrict Jtot, const int nE, const real Emin, const real dE, const int nJ, const real dJ, real * restrict EJ);

void writeAnalyzedMaps
(const double time, const ulong steps, char file[], const uint filenum, hdf5struct type, const int kind, int * restrict bodyNum,
 int * restrict bodyHead, real * restrict vr, real * restrict vt, real * restrict Etot, real * restrict Jtot,
 const int nvr, real * restrict vvr, const int nvt, real * restrict vvt, real * restrict f_vrvt,
 const int nJ, real * restrict JJ, const int nE, real * restrict EE, real * restrict EJ);


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);

  /** initialization */
  if( argc < 17 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 17);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -start=<int> -end=<int> -interval=<int>\n");
    __FPRINTF__(stderr, "          -nvr=<int> -vrmin=<real> -vrmax=<real>\n");
    __FPRINTF__(stderr, "          -nvt=<int> -vtmin=<real> -vtmax=<real>\n");
    __FPRINTF__(stderr, "          -nE=<int> -Emin=<real> -Emax=<real>\n");
    __FPRINTF__(stderr, "          -nJ=<int> -Jmin=<real> -Jmax=<real>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 17 ){ */

  char   *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,     "file", &file));
  int    start;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,    "start", &start));
  int      end;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,      "end", &end));
  int interval;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "interval", &interval));

  int nvr;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nvr", &nvr));
  real vrmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vrmin", &vrmin));
  real vrmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vrmax", &vrmax));

  int nvt;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nvt", &nvt));
  real vtmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vtmin", &vtmin));
  real vtmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "vtmax", &vtmax));

  int nE;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nE", &nE));
  int nJ;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv, "nJ", &nJ));
  real Emin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "Emin", &Emin));
  real Emax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "Emax", &Emax));
  real Jmin;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "Jmin", &Jmin));
  real Jmax;  requiredCmdArg(getCmdArgReal(argc, (const char * const *)argv, "Jmax", &Jmax));


  /** load global settings of particle distribution */
  int last;
  readConfigFile(&last, file);
  int unit;
  ulong Ntot;
  real eta, eps;
  double ft, snapshotInterval, saveInterval;
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file, mpi);
  nbody_particle *body;
  /* allocParticleDataAoS((int)Ntot, &body); */
  body = (nbody_particle *)malloc(sizeof(nbody_particle) * Ntot);
  if( body == NULL ){    __KILL__(stderr, "ERROR: failure to allocate body\n");  }
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

  real *vr;  vr = (real *)malloc(sizeof(real) * Ntot);  if( vr == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vr");  }
  real *vt;  vt = (real *)malloc(sizeof(real) * Ntot);  if( vt == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vt");  }
  real *Etot;  Etot = (real *)malloc(sizeof(real) * Ntot);  if( Etot == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Etot");  }
  real *Jtot;  Jtot = (real *)malloc(sizeof(real) * Ntot);  if( Jtot == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate Jtot");  }


  setPhysicalConstantsAndUnitSystem(unit, 1);

  /** read number of components */
  int kind = 0;
  int skind = 0;
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
  checker &= (2 == fscanf(fp, "%d\t%d", &kind, &skind));
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
  }/* if( !checker ){ */
  bodyHead[0] = 0;
  for(int ii = 1; ii < kind; ii++)
    bodyHead[ii] = bodyHead[ii - 1] + bodyNum[ii - 1];
#ifdef  HDF5_FOR_ZINDAIJI
  for(int ii = 0; ii < kind; ii++)
    bodyType[ii] &= 3;
#endif//HDF5_FOR_ZINDAIJI


  vrmin = CAST_D2R(CAST_R2D(vrmin) / velocity2astro);
  vrmax = CAST_D2R(CAST_R2D(vrmax) / velocity2astro);
  real *vvr;  vvr = (real *)malloc(sizeof(real) * (nvr + 1));  if( vvr == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vvr");  }
  const real dvr = (vrmax - vrmin) / (real)nvr;  for(int ii = 0; ii < nvr + 1; ii++)    vvr[ii] = vrmin + dvr * (real)ii;

  vtmin = CAST_D2R(CAST_R2D(vtmin) / velocity2astro);
  vtmax = CAST_D2R(CAST_R2D(vtmax) / velocity2astro);
  real *vvt;  vvt = (real *)malloc(sizeof(real) * (nvt + 1));  if( vvt == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate vvt");  }
  const real dvt = (vtmax - vtmin) / (real)nvt;  for(int ii = 0; ii < nvt + 1; ii++)    vvt[ii] = vtmin + dvt * (real)ii;

  real *f_vrvt;  f_vrvt = (real *)malloc(sizeof(real) * kind * nvr * nvt);  if( f_vrvt == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate f_vrvt");  }

  Emin = CAST_D2R(CAST_R2D(Emin) / senergy2astro);
  Emax = CAST_D2R(CAST_R2D(Emax) / senergy2astro);
  real *EE;  EE = (real *)malloc(sizeof(real) * (nE + 1));  if( EE == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate EE");  }
  const real dE = (Emax - Emin) / (real)nE;  for(int ii = 0; ii < nE + 1; ii++)    EE[ii] = Emin + dE * (real)ii;
  Jmin = CAST_D2R(CAST_R2D(Jmin) / (length2astro * velocity2astro));
  Jmax = CAST_D2R(CAST_R2D(Jmax) / (length2astro * velocity2astro));
  real *JJ;  JJ = (real *)malloc(sizeof(real) * (nJ + 1));  if( JJ == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate JJ");  }
  const real dJ = (Jmax - Jmin) / (real)nJ;
  for(int ii = 0; ii < nJ + 1; ii++)    JJ[ii] = Jmin + dJ * (real)ii;
  real *EJ;  EJ = (real *)malloc(sizeof(real) * kind * nE * nJ);  if( EJ == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate EJ");  }


#ifdef  USE_HDF5_FORMAT
  /** obtain minimum and maximum values for matplotlib */
  double vradmin = DBL_MAX, vradmax = -DBL_MAX;
  double vtanmin = DBL_MAX, vtanmax = -DBL_MAX;
  double enemin = DBL_MAX, enemax = -DBL_MAX;
  double sammin = DBL_MAX, sammax = -DBL_MAX;
#endif//USE_HDF5_FORMAT


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
    readSnapshot(&unit_read, &time, &steps, Ntot, file, ibody, (uint)filenum);
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
    qsort(body, Ntot, sizeof(nbody_particle), idxAscendingOrder);


    /** calculate vr, vt, E, and J for all N-body particles */
    for(int kk = 0; kk < kind; kk++){
      const int ihead = bodyHead[kk];
      const int inum  = bodyNum [kk];
      for(int ii = ihead; ii < ihead + inum; ii++){
	const double xx = CAST_R2D(body[ii]. x);
	const double yy = CAST_R2D(body[ii]. y);
	const double zz = CAST_R2D(body[ii]. z);
	const double vx = CAST_R2D(body[ii].vx);
	const double vy = CAST_R2D(body[ii].vy);
	const double vz = CAST_R2D(body[ii].vz);

	const double rr = sqrt(xx * xx + yy * yy + zz * zz);

	const double cost = zz / rr;
	const double sint = sqrt(1.0 - cost * cost);
	const double cosp = xx / (rr * sint);
	const double sinp = yy / (rr * sint);

	const double vrad =  vx * sint * cosp + vy * sint * sinp + vz * cost;
	const double vtht =  vx * cost * cosp + vy * cost * sinp - vz * sint;
	const double vphi = -vx	       * sinp + vy        * cosp;
	const double vtan = sqrt(vtht * vtht + vphi * vphi);
#ifdef  USE_HDF5_FORMAT
	vradmin = fmin(vradmin, vrad);	vradmax = fmax(vradmax, vrad);
	vtanmin = fmin(vtanmin, vtan);	vtanmax = fmax(vtanmax, vtan);
#endif//USE_HDF5_FORMAT

	vr[ii] = CAST_D2R(vrad);
	vt[ii] = CAST_D2R(vtan);

#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
	const double ene = CAST_R2D(body[ii].pot + body[ii].pot_ext) + 0.5 * (vx * vx + vy * vy + vz * vz);
#else///SET_EXTERNAL_POTENTIAL_FIELD
	const double ene = CAST_R2D(body[ii].pot                   ) + 0.5 * (vx * vx + vy * vy + vz * vz);
#endif//SET_EXTERNAL_POTENTIAL_FIELD
	Etot[ii] = CAST_D2R(ene);

	const double jx = yy * vz - zz * vy;
	const double jy = zz * vx - xx * vz;
	const double jz = xx * vy - yy * vx;
	const double oam = sqrt(jx * jx + jy * jy + jz * jz);
	Jtot[ii] = CAST_D2R(oam);
#ifdef  USE_HDF5_FORMAT
	enemin = fmin(enemin, ene);
	enemax = fmax(enemax, ene);
	if( oam > 0.0 )
	  sammin = fmin(sammin, oam);
	sammax = fmax(sammax, oam);
#endif//USE_HDF5_FORMAT
      }/* for(int ii = ihead; ii < ihead + inum; ii++){ */
    }/* for(int kk = 0; kk < kind; kk++){ */
    generateVelocityDistribution  (kind, bodyHead, bodyNum, body, vr, vt, nvr, vrmin, dvr, nvt, vtmin, dvt, f_vrvt);
    generatePhaseSpaceDistribution(kind, bodyHead, bodyNum, body, Etot, Jtot, nE, Emin, dE, nJ, dJ, EJ);


#ifdef  USE_HDF5_FORMAT
    /** dump analyzed results for matplotlib and/or VisIt */
    writeAnalyzedMaps
      (time, steps, file, filenum, hdf5type, kind, bodyNum,
       bodyHead, vr, vt, Etot, Jtot,
       nvr, vvr, nvt, vvt, f_vrvt,
       nJ, JJ, nE, EE, EJ);
#endif//USE_HDF5_FORMAT
  }/* for(int filenum = start + mpi.rank * interval; filenum < end + 1; filenum += interval * mpi.size){ */


#ifdef  USE_HDF5_FORMAT
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & vradmin: MPI_IN_PLACE, &vradmin, 1, MPI_DOUBLE, MPI_MIN, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & vradmax: MPI_IN_PLACE, &vradmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & vtanmin: MPI_IN_PLACE, &vtanmin, 1, MPI_DOUBLE, MPI_MIN, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & vtanmax: MPI_IN_PLACE, &vtanmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));

  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & enemin: MPI_IN_PLACE, &enemin, 1, MPI_DOUBLE, MPI_MIN, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & enemax: MPI_IN_PLACE, &enemax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & sammin: MPI_IN_PLACE, &sammin, 1, MPI_DOUBLE, MPI_MIN, 0, mpi.comm));
  chkMPIerr(MPI_Reduce((mpi.rank != 0) ? & sammax: MPI_IN_PLACE, &sammax, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm));
#endif//USE_HDF5_FORMAT

  if( mpi.rank == 0 ){
    sprintf(filename, "%s/%s_phase_minmax.txt", DATAFOLDER, file);
    fp = fopen(filename, "w");
    if( fp == NULL ){
      __KILL__(stderr, "ERROR: failure to open \"%s\"\n", filename);
    }
    fprintf(fp, "%e\t%e\n", vradmin, vradmax);
    fprintf(fp, "%e\t%e\n", vtanmin, vtanmax);
    fprintf(fp, "%e\t%e\n", enemin, enemax);
    fprintf(fp, "%e\t%e\n", sammin, sammax);
    fclose(fp);
  }/* if( mpi.rank == 0 ){ */


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
#ifdef  HDF5_FOR_ZINDAIJI
  free(bodyType);
#endif//HDF5_FOR_ZINDAIJI

  free(vr);  free(vvr);
  free(vt);  free(vvt);
  free(f_vrvt);

  free(Etot);  free(EE);
  free(Jtot);  free(JJ);
  free(EJ);

  exitMPI();

  return (0);
}


/* 3 sigma means that neglecting component of 0.26% */
/* 5 sigma means that neglecting component of 6e-7 */
/* #define SPREAD (3) */
#define SPREAD (5)


#define SMOOTHING_FOR_VISUALIZATION TWO
/* #define SMOOTHING_FOR_VISUALIZATION THREE */


void generateVelocityDistribution(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body, real * restrict vr, real * restrict vt, const int nvr, const real vrmin, const real dvr, const int nvt, const real vtmin, const real dvt, real * restrict vrvt)
{
  __NOTE__("%s\n", "start");

  const real dvrinv = UNITY / dvr;
  const real dvtinv = UNITY / dvt;

  const real vrsig = SMOOTHING_FOR_VISUALIZATION * FABS(dvr);
  const real vtsig = SMOOTHING_FOR_VISUALIZATION * FABS(dvt);

  const real invvrsig = UNITY / vrsig;
  const real invvtsig = UNITY / vtsig;

  const int nvr_smooth = (int)CEIL(SPREAD * FABS(dvrinv) * vrsig);
  const int nvt_smooth = (int)CEIL(SPREAD * FABS(dvtinv) * vtsig);
  real *erfvr;  erfvr = (real *)malloc(sizeof(real) * (2 * nvr_smooth + 2));  if( erfvr == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfvr");  }
  real *erfvt;  erfvt = (real *)malloc(sizeof(real) * (2 * nvt_smooth + 2));  if( erfvt == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfvt");  }
  for(int ii = 0; ii < 2 * nvr_smooth + 2; ii++)    erfvr[ii] = ERF((dvr * ((real)(ii - nvr_smooth) - HALF)) * invvrsig);
  for(int ii = 0; ii < 2 * nvt_smooth + 2; ii++)    erfvt[ii] = ERF((dvt * ((real)(ii - nvt_smooth) - HALF)) * invvtsig);

  real *psfvr;  psfvr = (real *)malloc(sizeof(real) * (2 * nvr_smooth + 1));  if( psfvr == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfvr");  }
  real *psfvt;  psfvt = (real *)malloc(sizeof(real) * (2 * nvt_smooth + 1));  if( psfvt == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfvt");  }
  for(int ii = 0; ii < 2 * nvr_smooth + 1; ii++)    psfvr[ii] = HALF * (erfvr[ii + 1] - erfvr[ii]);
  for(int ii = 0; ii < 2 * nvt_smooth + 1; ii++)    psfvt[ii] = HALF * (erfvt[ii + 1] - erfvt[ii]);

  for(int ii = 0; ii < kind * nvr * nvt; ii++)
    vrvt[ii] = ZERO;

  for(int kk = 0; kk < kind; kk++){
    for(int ii = bodyHead[kk]; ii < bodyHead[kk] + bodyNum[kk]; ii++){
      const real mi = body[ii].m;
      const real ri = vr[ii];
      const real ti = vt[ii];

      const int j0 = (int)FLOOR((ri - vrmin) * dvrinv);
      const int l0 = (int)FLOOR((ti - vtmin) * dvtinv);

      for(int sr = 0; sr < 2 * nvr_smooth + 1; sr++){
	const int jj = j0 + sr - nvr_smooth;
	if( (jj >= 0) && (jj < nvr) ){
	  const real mr = mi * psfvr[sr];

	  for(int st = 0; st < 2 * nvt_smooth + 1; st++){
	    const int ll = l0 + st - nvt_smooth;
	    if( (ll >= 0) && (ll < nvt) )
	      vrvt[INDEX3D(kind, nvr, nvt, kk, jj, ll)] += mr * psfvt[st];
	  }/* for(int st = 0; st < 2 * nvt_smooth + 1; st++){ */
	}/* if( (jj >= 0) && (jj < nvr) ){ */
      }/* for(int sr = 0; sr < 2 * nvr_smooth + 1; sr++){ */

    }/* for(int ii = bodyHead[kk]; ii < bodyHead[kk] + bodyNum[kk]; ii++){ */
  }/* for(int kk = 0; kk < kind; kk++){ */

  const real dSinv = dvrinv * dvtinv;
  for(int ii = 0; ii < kind * nvr * nvt; ii++)
    vrvt[ii] *= dSinv;

  free(erfvr);  free(erfvt);
  free(psfvr);  free(psfvt);


  __NOTE__("%s\n", "end");
}


void generatePhaseSpaceDistribution(const int kind, int * restrict bodyHead, int * restrict bodyNum, nbody_particle *body, real * restrict Etot, real * restrict Jtot, const int nE, const real Emin, const real dE, const int nJ, const real dJ, real * restrict EJ)
{
  __NOTE__("%s\n", "start");

  const real dEinv = UNITY / dE;
  const real dJinv = UNITY / dJ;

  const real Esig = SMOOTHING_FOR_VISUALIZATION * FABS(dE);
  const real Jsig = SMOOTHING_FOR_VISUALIZATION * FABS(dJ);

  const real invEsig = UNITY / Esig;
  const real invJsig = UNITY / Jsig;

  const int nE_smooth = (int)CEIL(SPREAD * FABS(dEinv) * Esig);
  const int nJ_smooth = (int)CEIL(SPREAD * FABS(dJinv) * Jsig);
  real *erfE;  erfE = (real *)malloc(sizeof(real) * (2 * nE_smooth + 2));  if( erfE == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfE");  }
  real *erfJ;  erfJ = (real *)malloc(sizeof(real) * (2 * nJ_smooth + 2));  if( erfJ == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate erfJ");  }
  for(int ii = 0; ii < 2 * nE_smooth + 2; ii++)    erfE[ii] = ERF((dE * ((real)(ii - nE_smooth) - HALF)) * invEsig);
  for(int ii = 0; ii < 2 * nJ_smooth + 2; ii++)    erfJ[ii] = ERF((dJ * ((real)(ii - nJ_smooth) - HALF)) * invJsig);

  real *psfE;  psfE = (real *)malloc(sizeof(real) * (2 * nE_smooth + 1));  if( psfE == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfE");  }
  real *psfJ;  psfJ = (real *)malloc(sizeof(real) * (2 * nJ_smooth + 1));  if( psfJ == NULL ){    __KILL__(stderr, "%s\n", "ERROR: failure to allocate psfJ");  }
  for(int ii = 0; ii < 2 * nE_smooth + 1; ii++)    psfE[ii] = HALF * (erfE[ii + 1] - erfE[ii]);
  for(int ii = 0; ii < 2 * nJ_smooth + 1; ii++)    psfJ[ii] = HALF * (erfJ[ii + 1] - erfJ[ii]);

  for(int ii = 0; ii < kind * nJ * nE; ii++)
    EJ[ii] = ZERO;

  for(int kk = 0; kk < kind; kk++){
    for(int ii = bodyHead[kk]; ii < bodyHead[kk] + bodyNum[kk]; ii++){
      const real mi = body[ii].m;
      const real Ei = Etot[ii];
      const real Ji = Jtot[ii];

      const int l0 = (int)FLOOR((Ji       ) * dJinv);
      const int n0 = (int)FLOOR((Ei - Emin) * dEinv);

      for(int sE = 0; sE < 2 * nE_smooth + 1; sE++){
	const int nn = n0 + sE - nE_smooth;
	if( (nn >= 0) && (nn < nE) ){
	  const real mE = mi * psfE[sE];

	  for(int sJ = 0; sJ < 2 * nJ_smooth + 1; sJ++){
	    const int ll = l0 + sJ - nJ_smooth;
	    if( (ll >= 0) && (ll < nJ) )
	      EJ[INDEX3D(kind, nJ, nE, kk, ll, nn)] += mE * psfJ[sJ];
	  }/* for(int sJ = 0; sJ < 2 * nJ_smooth + 1; sJ++){ */
	}/* if( (nn >= 0) && (nn < nE) ){ */
      }/* for(int sE = 0; sE < 2 * nE_smooth + 1; sE++){ */

    }/* for(int ii = bodyHead[kk]; ii < bodyHead[kk] + bodyNum[kk]; ii++){ */
  }/* for(int kk = 0; kk < kind; kk++){ */

  const real dSinv = dEinv * dJinv;
  for(int ii = 0; ii < kind * nE * nJ; ii++)
    EJ[ii] *= dSinv;


  free(erfE);  free(erfJ);
  free(psfE);  free(psfJ);


  __NOTE__("%s\n", "end");
}


#ifdef  USE_HDF5_FORMAT
/**
 * @fn writeAnalyzedMaps
 *
 * @brief Write analyzed maps of the N-body simulation.
 */
void writeAnalyzedMaps
(const double time, const ulong steps, char file[], const uint filenum, hdf5struct type, const int kind, int * restrict bodyNum,
 int * restrict bodyHead, real * restrict vr, real * restrict vt, real * restrict Etot, real * restrict Jtot,
 const int nvr, real * restrict vvr, const int nvt, real * restrict vvt, real * restrict f_vrvt,
 const int nJ, real * restrict JJ, const int nE, real * restrict EE, real * restrict EJ)
{
  __NOTE__("%s\n", "start");

  static bool firstCall = true;
  if( firstCall ){
    for(int jj = 0; jj < nvr + 1; jj++)      vvr[jj] = CAST_D2R(CAST_R2D(vvr[jj]) * velocity2astro);
    for(int jj = 0; jj < nvt + 1; jj++)      vvt[jj] = CAST_D2R(CAST_R2D(vvt[jj]) * velocity2astro);

    for(int jj = 0; jj < nJ + 1; jj++)      JJ[jj] = CAST_D2R(CAST_R2D(JJ[jj]) * length2astro * velocity2astro);
    for(int jj = 0; jj < nE + 1; jj++)      EE[jj] = CAST_D2R(CAST_R2D(EE[jj]) * senergy2astro);
  }/* if( firstCall ){ */

  const double vel_density2astro = mass2astro / (velocity2astro * velocity2astro);
  for(int ii = 0; ii < kind * nvr * nvt; ii++)
    f_vrvt[ii] = CAST_D2R(CAST_R2D(f_vrvt[ii]) * vel_density2astro);
  const double cone_density2astro = mass2astro / (length2astro * velocity2astro * senergy2astro);
  for(int ii = 0; ii < kind * nJ * nE; ii++)    EJ[ii] = CAST_D2R(CAST_R2D(EJ[ii]) * cone_density2astro);


  /* create a new file (if the file already exists, the file is opened with read-write access, new data will overwrite any existing data) */
  char filename[128];
  sprintf(filename, "%s/%s.%s%.3u.h5", DATAFOLDER, file, "m31ene", filenum);
  hid_t target = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


  /* create the data space for the dataset */
  /* preparation for data compression */
  hid_t dataset, dataspace, property;
#ifdef  USE_FILE_COMPRESSION
  const hsize_t cdims_max = (MAXIMUM_CHUNK_SIZE_4BIT < MAXIMUM_CHUNK_SIZE) ? MAXIMUM_CHUNK_SIZE_4BIT : MAXIMUM_CHUNK_SIZE;
  hsize_t cdims_loc[3];
#ifdef  USE_SZIP_COMPRESSION
  /* compression using szip */
  const uint szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  const uint szip_pixels_per_block = 8;
  const hsize_t szip_cdims[2] = {1, 128 * szip_pixels_per_block};
#endif//USE_SZIP_COMPRESSION
#ifdef  USE_GZIP_COMPRESSION
  /* compression using gzip */
  const uint gzip_compress_level = 9;/**< 9 is the maximum compression ratio */
  const hsize_t gzip_cdims[2] = {1, 1024};
#endif//USE_GZIP_COMPRESSION
#else///USE_FILE_COMPRESSION
  property = H5P_DEFAULT;
#endif//USE_FILE_COMPRESSION



  /* write attribute data */
  /* create the data space for the attribute */
  hsize_t attr_dims = 1;
  dataspace = H5Screate_simple(1, &attr_dims, NULL);
  hid_t attribute;
  /* write current time */
  double wtime = time * time2astro;
  attribute = H5Acreate(target, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &wtime));
  chkHDF5err(H5Aclose(attribute));
  /* write # of steps */
  attribute = H5Acreate(target, "steps", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  /* write # of N-body particles */
  ulong num_ulong = 0;
  for(int ii = 0; ii < kind; ii++)
    num_ulong += (ulong)bodyNum[ii];
  attribute = H5Acreate(target, "number", H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_ULONG, &num_ulong));
  chkHDF5err(H5Aclose(attribute));
  /* write # of components */
  attribute = H5Acreate(target, "kinds", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &kind));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nvr", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nvr));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nvt", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nvt));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nE", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nE));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "nJ", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nJ));
  chkHDF5err(H5Aclose(attribute));
  /* write flag about USE_DOUBLE_PRECISION */
#ifdef  USE_DOUBLE_PRECISION
  const int useDP = 1;
#else///USE_DOUBLE_PRECISION
  const int useDP = 0;
#endif//USE_DOUBLE_PRECISION
  attribute = H5Acreate(target, "useDP", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));

  hid_t str4format = H5Tcopy(H5T_C_S1);
  chkHDF5err(H5Tset_size(str4format, CONSTANTS_H_CHAR_WORDS));
  attribute = H5Acreate(target, "length_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, length_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "mass_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, mass_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "time_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, time_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "velocity_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, velocity_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  attribute = H5Acreate(target, "senergy_astro_unit_name", str4format, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  chkHDF5err(H5Awrite(attribute, str4format, senergy_astro_unit_name));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Tclose(str4format));
  /* close the dataspace */
  chkHDF5err(H5Sclose(dataspace));


  for(int ii = 0; ii < kind; ii++){
    char grp[16];    sprintf(grp, "map%d", ii);
    hid_t group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t attr_dims = 1;
    hid_t attribute;

    /** 2D (nvr * nvt) array */
    hsize_t dims[2] = {nvr, nvt};
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[0];
    cdims_loc[1] = gzip_cdims[1];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[1] < cdims_loc[1] ){
      cdims_loc[1] = dims[1];
      cdims_loc[0] = dims[1] / cdims_loc[1];
    }/* if( dims[1] < cdims_loc[1] ){ */
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] > cdims_max )
      cdims_loc[0] = cdims_max / cdims_loc[1];
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    dataset = H5Dcreate(group, "vrvt", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &f_vrvt[INDEX2D(kind, nvr * nvt, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nvr) array */
    dims[0] = nvr + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    dataset = H5Dcreate(group, "vr", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, vvr));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nvt) array */
    dims[0] = nvt + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    dataset = H5Dcreate(group, "vt", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, vvt));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
    /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));


    /** 2D (nJ * nE) array */
    dims[0] = nJ;
    dims[1] = nE;
    dataspace = H5Screate_simple(2, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[0];
    cdims_loc[1] = gzip_cdims[1];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[1] < cdims_loc[1] ){
      cdims_loc[1] = dims[1];
      cdims_loc[0] = dims[1] / cdims_loc[1];
    }/* if( dims[1] < cdims_loc[1] ){ */
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] * cdims_loc[1] > cdims_max )
      cdims_loc[0] = cdims_max / cdims_loc[1];
    chkHDF5err(H5Pset_chunk(property, 2, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    dataset = H5Dcreate(group, "EJ", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &EJ[INDEX2D(kind, nE * nJ, ii, 0)]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
      /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nE + 1) arrays */
    dims[0] = nE + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    dataset = H5Dcreate(group, "E", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, EE));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
      /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /** 1D (nJ + 1) arrays */
    dims[0] = nJ + 1;
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    dataset = H5Dcreate(group, "J", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, JJ));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
      /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));

    /* attr_dims = 1; */
    /* dataspace = H5Screate_simple(1, &attr_dims, NULL); */
    /* attribute = H5Acreate(group, "Nvr", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT); */
    /* chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nvr)); */
    /* chkHDF5err(H5Aclose(attribute)); */
    /* attribute = H5Acreate(group, "Nvt", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT); */
    /* chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nvt)); */
    /* chkHDF5err(H5Aclose(attribute)); */
    /* attribute = H5Acreate(group, "NE", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT); */
    /* chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nE)); */
    /* chkHDF5err(H5Aclose(attribute)); */
    /* attribute = H5Acreate(group, "NJ", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT); */
    /* chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &nJ)); */
    /* chkHDF5err(H5Aclose(attribute)); */
    /* chkHDF5err(H5Sclose(dataspace)); */


    chkHDF5err(H5Gclose(group));
    sprintf(grp, "particle%d", ii);
    group = H5Gcreate(target, grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    /** 1D (bodyNum) arrays */
    dims[0] = bodyNum[ii];
    dataspace = H5Screate_simple(1, dims, NULL);
#ifdef  USE_GZIP_COMPRESSION
    cdims_loc[0] = gzip_cdims[1];
    property = H5Pcreate(H5P_DATASET_CREATE);
    if( dims[0] < cdims_loc[0] )
      cdims_loc[0] = dims[0];
    if( cdims_loc[0] > cdims_max )
      cdims_loc[0] = cdims_max;
    chkHDF5err(H5Pset_chunk(property, 1, cdims_loc));
    chkHDF5err(H5Pset_deflate(property, gzip_compress_level));
#endif//USE_GZIP_COMPRESSION
    for(int jj = bodyHead[ii]; jj < bodyHead[ii] + bodyNum[ii]; jj++){
      vr[jj] = CAST_D2R(CAST_R2D(vr[jj]) * velocity2astro);
      vt[jj] = CAST_D2R(CAST_R2D(vt[jj]) * velocity2astro);
      Etot[jj] = CAST_D2R(CAST_R2D(Etot[jj]) * senergy2astro);
      Jtot[jj] = CAST_D2R(CAST_R2D(Jtot[jj]) * length2astro * velocity2astro);
    }/* for(int jj = bodyHead[ii]; jj < bodyHead[ii] + bodyNum[ii]; jj++){ */
    dataset = H5Dcreate(group, "vr", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vr[bodyHead[ii]]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "vt", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vt[bodyHead[ii]]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "Etot", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Etot[bodyHead[ii]]));
    chkHDF5err(H5Dclose(dataset));
    dataset = H5Dcreate(group, "Jtot", type.real, dataspace, H5P_DEFAULT, property, H5P_DEFAULT);
    chkHDF5err(H5Dwrite(dataset, type.real, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Jtot[bodyHead[ii]]));
    chkHDF5err(H5Dclose(dataset));
#ifdef  USE_GZIP_COMPRESSION
    chkHDF5err(H5Pclose(property));
#endif//USE_GZIP_COMPRESSION
      /* close the dataspace */
    chkHDF5err(H5Sclose(dataspace));
    /* write attribute data */
    attr_dims = 1;
    dataspace = H5Screate_simple(1, &attr_dims, NULL);
    /* write # of grid points */
    attribute = H5Acreate(group, "num", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    chkHDF5err(H5Awrite(attribute, H5T_NATIVE_INT, &bodyNum[ii]));
    chkHDF5err(H5Aclose(attribute));
    chkHDF5err(H5Sclose(dataspace));


    chkHDF5err(H5Gclose(group));
  }/* for(int ii = 0; ii < kind; ii++){ */


  /* close the file */
  chkHDF5err(H5Fclose(target));

  firstCall = false;

  __NOTE__("%s\n", "end");
}
#endif//USE_HDF5_FORMAT
