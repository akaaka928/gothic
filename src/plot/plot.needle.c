/**
 * @file plot.needle.c
 *
 * @brief Plot code of needle-like structure
 *
 * @author Yohei Miki (University of Tokyo)
 *
 * @date 2017/10/26 (Thu)
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

#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../file/io.h"

/* #define REDUCE_PARTICLE_DISTRIBUTION_MAP (8192) */
/* #define REDUCE_PARTICLE_DISTRIBUTION_MAP (65536) */
#define REDUCE_PARTICLE_DISTRIBUTION_MAP (131072)
/* #define REDUCE_PARTICLE_DISTRIBUTION_MAP (262144) */
/* #define REDUCE_PARTICLE_DISTRIBUTION_MAP (524288) */

#define SKIP_IN_PLOT (1)

extern const double      length2astro;extern const char      length_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double        time2astro;extern const char        time_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double        mass2astro;extern const char        mass_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double     density2astro;extern const char     density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];
extern const double col_density2astro;extern const char col_density_astro_unit_name4plot[CONSTANTS_H_PLOT_WORDS];


typedef struct
{
  ulong head, num;/**< data for N-body particles */
} model;

typedef struct
{
  ulong idx;
  real  x,  y,  z;
  real vx, vy, vz;
  real ax, ay, az;
  real m, pot;
} nbody_particle;


#ifdef __ICC
/** Disable ICC's remark #161: unrecognized #pragma */
#     pragma warning (disable:161)
#endif//__ICC
int idxAscendingOrder(const void *a, const void *b)
{
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
  if(          ((nbody_particle *)a)->idx > ((nbody_particle *)b)->idx ){    return ( 1);  }
  else{    if( ((nbody_particle *)a)->idx < ((nbody_particle *)b)->idx ){    return (-1);  }
    else{                                                                    return ( 0);  }  }
#   if  ((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
#pragma GCC diagnostic pop
#endif//((__GNUC_MINOR__ + __GNUC__ * 10) >= 45)
}
#ifdef __ICC
/** Enable ICC's remark #161: unrecognized #pragma */
#     pragma warning (enable:161)
#endif//__ICC


void plotDistributionMaps
(const int ngroup, model *group, nbody_particle *body0, nbody_particle *body1, PLFLT time,
 PLplotPltRange xybox, PLplotPltRange xzbox, PLplotPltRange zybox,
 char *file, int argc, char **argv);


int main(int argc, char **argv)
{
  /** initialization */
  if( argc < 4 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 4);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file0=<char *>\n");
    __FPRINTF__(stderr, "          -file1=<char *>\n");
    __FPRINTF__(stderr, "          -problem=<int>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }/* if( argc < 4 ){ */

  char  *file0;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,    "file0", &file0));
  char  *file1;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)argv,    "file1", &file1));
  int  problem;  requiredCmdArg(getCmdArgInt(argc, (const char * const *)argv,  "problem", &problem));

  modifyArgcArgv4PLplot(&argc, argv, 4);


  /** load global settings of particle distribution */
  int last;
  readConfigFile(&last, file0);
  int unit;
  ulong Ntot;
  real eta, eps;
  double ft, snapshotInterval, saveInterval;
  readSettings(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file0);
  nbody_particle *body0, *body1;
  /* allocParticleDataAoS((int)Ntot, &body0); */
  /* allocParticleDataAoS((int)Ntot, &body1); */
  body0 = (nbody_particle *)malloc(sizeof(nbody_particle) * Ntot);
  if( body0 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate body0\n");  }
  body1 = (nbody_particle *)malloc(sizeof(nbody_particle) * Ntot);
  if( body1 == NULL ){    __KILL__(stderr, "ERROR: failure to allocate body1\n");  }
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
  nbody_hdf5 hdf5;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
  ulong *hdf5_idx;
  allocSnapshotArray(&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx, (int)Ntot, &hdf5);
#else///USE_HDF5_FORMAT
  iparticle ibody;
  ulong *idx;
  position *pos;
  acceleration *acc;
#ifdef  BLOCK_TIME_STEP
  velocity *vel;
  ibody_time *ti;
#else///BLOCK_TIME_STEP
  real *vx, *vy, *vz;
#endif//BLOCK_TIME_STEP
  allocParticleData((int)Ntot, &ibody, &idx, &pos, &acc,
#ifdef  BLOCK_TIME_STEP
		    &vel, &ti
#else///BLOCK_TIME_STEP
		    &vx, &vy, &vz
#endif//BLOCK_TIME_STEP
		    );
#endif//USE_HDF5_FORMAT


  /** set plot range */
  PLFLT radius;
  switch( problem ){
  case  0:    radius =  3.0;    break;    /**< cold collapse of a uniform sphere */
  case  1:    radius =  5.0;    break;    /**< king sphere with W0 = 3 */
  case  2:    radius = 10.0;    break;    /**< Hernquist sphere with C = 10 */
  case  3:    radius =  5.0;    break;    /**< NFW sphere with C = 5 */
  case  4:    radius = 10.0;    break;    /**< Einasto sphere with C = 10 */
  case  5:    radius = 10.0;    break;    /**< Plummer sphere with C = 10 */
  case  6:    radius = 10.0;    break;    /**< Burkert sphere with C = 10 */
  case  7:    radius = 10.0;    break;    /**< Moore sphere with C = 10 */
  case  8:    radius = 10.0;    break;    /**< Two-power sphere with C = 10 */
  case 10:    radius = 10.0;    break;    /**< king sphere with W0 = 3 within an Einasto sphere with C = 10 */
  case 11:    radius = 5.0e+1;    break;    /**< A galaxy with multiple components */
  case 12:    radius = 5.0e+1;    break;    /**< A galaxy with multiple components */
  case 13:    radius = 5.0e+1;    break;    /**< A galaxy with multiple components */
  case 20:    radius = 3.0e+1;    break;    /**< M31 model determined by Fardal et al. (2007) */
  case 22:    radius = 3.0e+1;    break;    /**< A trial multi components galaxy model */
  case 23:    radius = 3.0e+1;    break;    /**< MW model (Sofue 2015; Akhter et al. 2012; Gillessen et al. 2009; Juric et al. 2008) */
  case 24:    radius = 3.0e+1;    break;    /**< M31 model (Sofue 2015; Gilbert et al. 2012) */
  case 25:    radius = 3.0e+1;    break;    /**< A trial multi components galaxy model */
  case 26:    radius = 3.0e+1;    break;    /**< A trial multi components galaxy model (spherical model) */
  case 27:    radius = 3.0e+1;    break;    /**< M31 model (NFW halo, de Vaucouleurs bulge, and exponential disk) */
  case 28:    radius = 2.5e+1;    break;    /**< A trial multi components galaxy model (NFW halo, King bulge, thick Sersic disk, and thin exponential disk) */
  case 29:    radius = 1.0e+1;    break;    /**< Multi components galaxy model by Vasiliev & Athanassoula (2015) */
  case 30:    radius = 1.0e+1;    break;    /**< time evolution of MW/A defined in Kuijken & Dubinski (1995) */
  case 31:    radius = 1.0e+1;    break;    /**< time evolution of MW/B defined in Kuijken & Dubinski (1995) */
  case 32:    radius = 1.0e+1;    break;    /**< time evolution of MW/C defined in Kuijken & Dubinski (1995) */
  case 33:    radius = 1.0e+1;    break;    /**< time evolution of MW/D defined in Kuijken & Dubinski (1995) */
  case 34:    radius = 1.0e+1;    break;    /**< time evolution of M31/A defined in Widrow et al. (2003) */
  case 35:    radius = 1.0e+1;    break;    /**< time evolution of M31/D defined in Widrow et al. (2003) */
  case 36:    radius = 5.0e+1;    break;    /**< time evolution of MWa defined in Widrow & Dubinski (2005) */
  case 37:    radius = 5.0e+1;    break;    /**< time evolution of MWb defined in Widrow & Dubinski (2005) */
  case 38:    radius = 5.0e+1;    break;    /**< MW like galaxy (based on Widrow et al. 2008, Bedorf et al. 2014) */
  case 40:    radius = 20.0;    break;    /**< Plummer sphere in table form with C = 20 */
  case 41:    radius = 20.0;    break;    /**< Two-power sphere in table form with C = 20 */
  case 42:    radius = 20.0;    break;    /**< de Vaucouleurs sphere in table form */
  case 50:    radius = 20.0;    break;    /**< de Vaucouleurs sphere */
  case 51:    radius = 20.0;    break;    /**< projected Two-power model */
  default:    radius =  3.0;    break;
  }
  PLplotPltRange xybox, xzbox, zybox;
  xybox.xmin = xzbox.xmin = -radius * (PLFLT)length2astro;  xybox.xmax = xzbox.xmax =  radius * (PLFLT)length2astro;
  xybox.ymin = zybox.ymin = -radius * (PLFLT)length2astro;  xybox.ymax = zybox.ymax =  radius * (PLFLT)length2astro;
  xzbox.ymin = zybox.xmin = -radius * (PLFLT)length2astro;  xzbox.ymax = zybox.xmax =  radius * (PLFLT)length2astro;
  xybox.xlog = xzbox.xlog = zybox.xlog = LINEAR_PLOT;  xybox.xgrd = xzbox.xgrd = zybox.xgrd = false;
  xybox.ylog = xzbox.ylog = zybox.ylog = LINEAR_PLOT;  xybox.ygrd = xzbox.ygrd = zybox.ygrd = false;


  /** read analytic profile and analyze */
  int kind = 1;
  int skind = 1;
  model *group;
  bool multi_group = false;
  if( problem >= 1 ){
    multi_group = true;

    FILE *fp;
    char filename[256];
    sprintf(filename, "%s/%s.summary.txt", DOCUMENTFOLDER, file0);
    fp = fopen(filename, "r");
    if( fp == NULL ){    __KILL__(stderr, "ERROR: \"%s\" couldn't open.\n", filename);  }

    fscanf(fp, "%*d");/**< skip reading unit */
    fscanf(fp, "%d\t%d", &kind, &skind);
    group = (model *)malloc(kind * sizeof(model));
    if( group == NULL ){      __KILL__(stderr, "ERROR: failure to allocate group\n");    }
    for(int ii = 0; ii < kind; ii++)
      fscanf(fp, "%zu", &group[ii].num);

    fclose(fp);
  }/* if( problem >= 1 ){ */

  if( !multi_group ){
    group = (model *)malloc(kind * sizeof(model));
    if( group == NULL ){      __KILL__(stderr, "ERROR: failure to allocate group\n");    }

    group[0].num = Ntot;
  }/* if( !multi_group ){ */

  group[0].head = 0;
  for(int ii = 1; ii < kind; ii++)
    group[ii].head = group[ii - 1].head + group[ii - 1].num;

  /** read particle distribution and analyze */
  double time;
  ulong steps;
  int unit_read;

  /** disk galaxy w/z needle like structure */
#ifdef  USE_HDF5_FORMAT
  readSnapshot(&unit_read, &time, &steps, Ntot, file1, 0, &hdf5, hdf5type);
  for(int ii = 0; ii < (int)Ntot; ii++){
    body1[ii]. x  = hdf5.pos[ii * 3];      body1[ii]. y = hdf5.pos[ii * 3 + 1];      body1[ii].z   = hdf5.pos[ii * 3 + 2];
    body1[ii].vx  = hdf5.vel[ii * 3];      body1[ii].vy = hdf5.vel[ii * 3 + 1];      body1[ii].vz  = hdf5.vel[ii * 3 + 2];
    body1[ii].ax  = hdf5.acc[ii * 3];      body1[ii].ay = hdf5.acc[ii * 3 + 1];      body1[ii].az  = hdf5.acc[ii * 3 + 2];
    body1[ii].idx = hdf5.idx[ii    ];      body1[ii]. m = hdf5.  m[ii        ];      body1[ii].pot = hdf5.pot[ii        ];
  }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#else///USE_HDF5_FORMAT
  readSnapshot(&unit_read, &time, &steps, Ntot, file1, (uint)filenum, ibody);
  for(int ii = 0; ii < (int)Ntot; ii++){
    body1[ii]. x  = ibody.pos[ii].x;      body1[ii]. y = ibody.pos[ii].y;      body1[ii]. z  = ibody.pos[ii].z;
    body1[ii].vx  = ibody.vel[ii].x;      body1[ii].vy = ibody.vel[ii].y;      body1[ii].vz  = ibody.vel[ii].z;
    body1[ii].ax  = ibody.acc[ii].x;      body1[ii].ay = ibody.acc[ii].y;      body1[ii].az  = ibody.acc[ii].z;
    body1[ii].idx = ibody.idx[ii]  ;      body1[ii]. m = ibody.pos[ii].m;      body1[ii].pot = ibody.acc[ii].pot;
  }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#endif//USE_HDF5_FORMAT
  if( unit_read != unit ){
    __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
  }/* if( unit_read != unit ){ */
  qsort(body1, Ntot, sizeof(nbody_particle), idxAscendingOrder);


  /** disk galaxy w/o needle like structure */
#ifdef  USE_HDF5_FORMAT
  readSnapshot(&unit_read, &time, &steps, Ntot, file0, 0, &hdf5, hdf5type);
  for(int ii = 0; ii < (int)Ntot; ii++){
    body0[ii]. x  = hdf5.pos[ii * 3];      body0[ii]. y = hdf5.pos[ii * 3 + 1];      body0[ii].z   = hdf5.pos[ii * 3 + 2];
    body0[ii].vx  = hdf5.vel[ii * 3];      body0[ii].vy = hdf5.vel[ii * 3 + 1];      body0[ii].vz  = hdf5.vel[ii * 3 + 2];
    body0[ii].ax  = hdf5.acc[ii * 3];      body0[ii].ay = hdf5.acc[ii * 3 + 1];      body0[ii].az  = hdf5.acc[ii * 3 + 2];
    body0[ii].idx = hdf5.idx[ii    ];      body0[ii]. m = hdf5.  m[ii        ];      body0[ii].pot = hdf5.pot[ii        ];
  }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#else///USE_HDF5_FORMAT
  readSnapshot(&unit_read, &time, &steps, Ntot, file0, (uint)filenum, ibody);
  for(int ii = 0; ii < (int)Ntot; ii++){
    body0[ii]. x  = ibody.pos[ii].x;      body0[ii]. y = ibody.pos[ii].y;      body0[ii]. z  = ibody.pos[ii].z;
    body0[ii].vx  = ibody.vel[ii].x;      body0[ii].vy = ibody.vel[ii].y;      body0[ii].vz  = ibody.vel[ii].z;
    body0[ii].ax  = ibody.acc[ii].x;      body0[ii].ay = ibody.acc[ii].y;      body0[ii].az  = ibody.acc[ii].z;
    body0[ii].idx = ibody.idx[ii]  ;      body0[ii]. m = ibody.pos[ii].m;      body0[ii].pot = ibody.acc[ii].pot;
  }/* for(int ii = 0; ii < (int)Ntot; ii++){ */
#endif//USE_HDF5_FORMAT
  if( unit_read != unit ){
    __KILL__(stderr, "ERROR: conflict about unit system detected (unit = %d, unit_read = %d)\n", unit, unit_read);
  }/* if( unit_read != unit ){ */
  qsort(body0, Ntot, sizeof(nbody_particle), idxAscendingOrder);

  plotDistributionMaps(kind, group, body0, body1, time, xybox, xzbox, zybox, file0, argc, argv);


#ifdef  USE_HDF5_FORMAT
  removeHDF5DataType(hdf5type);
  freeSnapshotArray(hdf5_pos, hdf5_vel, hdf5_acc, hdf5_m, hdf5_pot, hdf5_idx);
#else///USE_HDF5_FORMAT
  freeParticleData(idx, pos, acc,
#ifdef  BLOCK_TIME_STEP
		    vel, ti
#else///BLOCK_TIME_STEP
		    vx, vy, vz
#endif//BLOCK_TIME_STEP
		    );
#endif//USE_HDF5_FORMAT
  free(body0);  free(body1);
  free(group);


  return (0);
}


void plotDistributionMaps
(const int ngroup, model *group, nbody_particle *body0, nbody_particle *body1, PLFLT time,
 PLplotPltRange xybox, PLplotPltRange xzbox, PLplotPltRange yzbox,
 char *file, int argc, char **argv)
{
  __NOTE__("%s\n", "start");


  /** global setting(s) */
  /** maximum number of data kinds */
  const PLINT pkind = ngroup - SKIP_IN_PLOT;
  const PLINT lkind = 0;


  /** preparation for dataset */
  /** memory allocation for index */
  PLINT * num;  allocPLINT(& num, pkind);

  ulong Np = 0;
  for(int ii = SKIP_IN_PLOT; ii < ngroup; ii++)
    Np += group[ii].num;

  /** set number of data points */
#ifdef  REDUCE_PARTICLE_DISTRIBUTION_MAP
  const PLINT Ntot = (Np > REDUCE_PARTICLE_DISTRIBUTION_MAP) ? REDUCE_PARTICLE_DISTRIBUTION_MAP : Np;
  const PLINT Ninc = (PLINT)Np / Ntot;
#else///REDUCE_PARTICLE_DISTRIBUTION_MAP
  const PLINT Ntot = (PLINT)Np;
  const PLINT Ninc = 1;
#endif//REDUCE_PARTICLE_DISTRIBUTION_MAP
  for(PLINT ii = 0; ii < pkind; ii++)    num[ii] = (PLINT)group[ngroup - 1 - ii].num / Ninc;
  PLINT Ntmp = num[0];
  for(PLINT ii = 1; ii < pkind; ii++)    Ntmp += num[ii];
  if( Ntmp != Ntot )
    num[pkind - 1] += Ntot - Ntmp;

  /** memory allocation for data */
  PLFLT *_xpt, *_ypt, *_zpt;
  {
    PLINT tot = 0;
    for(PLINT ii = 0; ii < pkind; ii++)      tot += num[ii];
    allocPLFLT(&_xpt, tot * 2);
    allocPLFLT(&_ypt, tot * 2);
    allocPLFLT(&_zpt, tot * 2);
  }
  PLFLT **xpt;  allocPointer4PLFLT(&xpt, pkind * 2);
  PLFLT **ypt;  allocPointer4PLFLT(&ypt, pkind * 2);
  PLFLT **zpt;  allocPointer4PLFLT(&zpt, pkind * 2);
  xpt[0] = _xpt;
  ypt[0] = _ypt;
  zpt[0] = _zpt;
  for(PLINT ii = 1; ii < 2 * pkind; ii++){
    xpt[ii] = xpt[ii - 1] + num[(ii - 1) % pkind];
    ypt[ii] = ypt[ii - 1] + num[(ii - 1) % pkind];
    zpt[ii] = zpt[ii - 1] + num[(ii - 1) % pkind];
  }/* for(PLINT ii = 1; ii < pkind; ii++){ */

  /** data preparation */
  for(PLINT kk = 0; kk < pkind; kk++){
    const int head = (PLINT)group[ngroup - 1 - kk].head;

    for(PLINT ii = 0; ii < num[kk]; ii++){
      xpt[        kk][ii] = (PLFLT)body0[head + ii].x * (PLFLT)length2astro;
      ypt[        kk][ii] = (PLFLT)body0[head + ii].y * (PLFLT)length2astro;
      zpt[        kk][ii] = (PLFLT)body0[head + ii].z * (PLFLT)length2astro;

      xpt[pkind + kk][ii] = (PLFLT)body1[head + ii].x * (PLFLT)length2astro;
      ypt[pkind + kk][ii] = (PLFLT)body1[head + ii].y * (PLFLT)length2astro;
      zpt[pkind + kk][ii] = (PLFLT)body1[head + ii].z * (PLFLT)length2astro;
    }/* for(PLINT ii = 0; ii < num[kk]; ii++){ */
  }/* for(PLINT kk = 0; kk < pkind; kk++){ */

  /** set symbol style */
  PLplotPointType *pt;  setDefaultPointType(pkind, &pt);
  for(PLINT ii = 0; ii < pkind; ii++){
    sprintf(pt[ii].type, PLplotSymbolType[smallDot]);
    pt[ii].scale =       PLplotSymbolSize[smallDot];
  }/* for(PLINT ii = 0; ii < pkind; ii++){ */
  if( pkind - 1 >= 0 )    pt[pkind - 1].color = BLACK;
  if( pkind - 2 >= 0 )    pt[pkind - 2].color = RED;
  if( pkind - 3 >= 0 )    pt[pkind - 3].color = BLUE;

  /** set labels */
  char xlab[PLplotCharWords];  sprintf(xlab, "#fix #fr(%s)", length_astro_unit_name4plot);
  char ylab[PLplotCharWords];  sprintf(ylab, "#fiy #fr(%s)", length_astro_unit_name4plot);
  char zlab[PLplotCharWords];  sprintf(zlab, "#fiz #fr(%s)", length_astro_unit_name4plot);

  /** set caption(s) */
  PLplotCaption xycap;  setDefaultCaption(&xycap);  xycap.write = true;  sprintf(xycap.side, "t");
  PLplotCaption xzcap;  setDefaultCaption(&xzcap);  xzcap.write = true;  sprintf(xzcap.side, "t");
  PLplotCaption yzcap;  setDefaultCaption(&yzcap);  yzcap.write = true;  sprintf(yzcap.side, "t");
  sprintf(xycap.text, "#fit#fr = %e (%s)", (double)time * time2astro, time_astro_unit_name4plot);
  sprintf(xzcap.text, "#fit#fr = %e (%s)", (double)time * time2astro, time_astro_unit_name4plot);
  sprintf(yzcap.text, "#fit#fr = %e (%s)", (double)time * time2astro, time_astro_unit_name4plot);

  /** set legends */
  PLplotLegend xyleg;  setDefaultLegend(&xyleg, false);  xyleg.write = false;
  PLplotLegend xzleg;  setDefaultLegend(&xzleg, false);  xzleg.write = false;
  PLplotLegend yzleg;  setDefaultLegend(&yzleg, false);  yzleg.write = false;
  char *xylegTxt;
  {
    allocChar4PLplot(&xylegTxt, pkind);
    allocPointer4Char4PLplot(&(xyleg.text), pkind);
    assignChar4PLplot(pkind, xyleg.text, xylegTxt);
  }
  sprintf(xyleg.text[0], "position");
  char *xzlegTxt;
  {
    allocChar4PLplot(&xzlegTxt, pkind);
    allocPointer4Char4PLplot(&(xzleg.text), pkind);
    assignChar4PLplot(pkind, xzleg.text, xzlegTxt);
  }
  sprintf(xzleg.text[0], "position");
  char *yzlegTxt;
  {
    allocChar4PLplot(&yzlegTxt, pkind);
    allocPointer4Char4PLplot(&(yzleg.text), pkind);
    assignChar4PLplot(pkind, yzleg.text, yzlegTxt);
  }
  sprintf(yzleg.text[0], "position");
  PLBOOL xyuni = false;
  PLBOOL xzuni = false;
  PLBOOL yzuni = false;

  char title[PLplotCharWords];
  sprintf(title, "#fit#fr = %6.2f (%s)", (PLFLT)((double)time * time2astro), time_astro_unit_name4plot);


  /** memory allocation to exploit multiple-plot function */
  /** maximum number of panels in each direction */
  const PLINT nxpanel = 2;
  const PLINT nypanel = 1;

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


  /** configure to time evolution of energy */
  /** common setting(s) */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

      /** global setting(s) */
      plot[idx] = true;

      /** line setting(s) */
      nlkind[idx] = lkind;
      lnum  [idx] = NULL;
      line  [idx] = NULL;
      lx    [idx] = NULL;
      ly    [idx] = NULL;

      /** point setting(s) */
      npkind[idx] = pkind;
      point [idx] = pt;
      pnum  [idx] = num;
      px    [idx] = &xpt[ii * pkind];
      py    [idx] = &zpt[ii * pkind];

      /** errorbar setting(s) */
      errbar[idx] = 0;

      /** plot area */
      range[idx] = xzbox;

      /** label setting(s) */
      xlabel[idx] = xlab;
      ylabel[idx] = zlab;

      /** captions */
      cap[idx] = xzcap;
      sprintf(cap[idx].text, "(%c)", 97 + idx);

      /** legends */
      leg[idx] = xzleg;
      uni[idx] = xzuni;
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */

  /** individual file names */
  sprintf(figfile,                                  "%s_%s_xz", file, "comparison");
  /* sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 0)], "%s_%s_%.3d", file, "xyproj", filenum); */
  /* sprintf(figname[INDEX2D(nxpanel, nypanel, 0, 1)], "%s_%s_%.3d", file, "xzproj", filenum); */
  /* sprintf(figname[INDEX2D(nxpanel, nypanel, 1, 0)], "%s_%s_%.3d", file, "zyproj", filenum); */


  /** create figure(s) */
  plotData(nxpanel, nypanel, plot, true, true,
  	   nlkind,  line, lnum, lx, ly,
  	   npkind, point, pnum, px, py,
  	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
  	   cap, leg, uni, range,
  	   xlabel, ylabel, "", figfile, argc, argv);


  /** yz-projection */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

      /** point setting(s) */
      px    [idx] = &ypt[ii * pkind];
      py    [idx] = &zpt[ii * pkind];

      /** plot area */
      range[idx] = yzbox;

      /** label setting(s) */
      xlabel[idx] = ylab;
      ylabel[idx] = zlab;

      /** captions */
      cap[idx] = yzcap;
      sprintf(cap[idx].text, "(%c)", 97 + idx);

      /** legends */
      leg[idx] = yzleg;
      uni[idx] = yzuni;
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */

  sprintf(figfile, "%s_%s_yz", file, "comparison");

  plotData(nxpanel, nypanel, plot, true, true,
  	   nlkind,  line, lnum, lx, ly,
  	   npkind, point, pnum, px, py,
  	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
  	   cap, leg, uni, range,
  	   xlabel, ylabel, "", figfile, argc, argv);


  /** yz-projection */
  for(PLINT ii = 0; ii < nxpanel; ii++)
    for(PLINT jj = 0; jj < nypanel; jj++){
      const PLINT idx = INDEX2D(nxpanel, nypanel, ii, jj);

      /** point setting(s) */
      px    [idx] = &xpt[ii * pkind];
      py    [idx] = &ypt[ii * pkind];

      /** plot area */
      range[idx] = xybox;

      /** label setting(s) */
      xlabel[idx] = xlab;
      ylabel[idx] = ylab;

      /** captions */
      cap[idx] = xycap;
      sprintf(cap[idx].text, "(%c)", 97 + idx);

      /** legends */
      leg[idx] = xyleg;
      uni[idx] = xyuni;
    }/* for(PLINT jj = 0; jj < nypanel; jj++){ */

  sprintf(figfile, "%s_%s_xy", file, "comparison");

  plotData(nxpanel, nypanel, plot, true, true,
  	   nlkind,  line, lnum, lx, ly,
  	   npkind, point, pnum, px, py,
  	   errbar, NULL, NULL, NULL, NULL, NULL, NULL,
  	   cap, leg, uni, range,
  	   xlabel, ylabel, "", figfile, argc, argv);


  /** memory deallocation */
  free(plot);
  free(nlkind);  free(lnum);  free(lx);  free(ly);  free( line);
  free(npkind);  free(pnum);  free(px);  free(py);  free(point);
  free(errbar);
  free(cap);  free(leg);  free(uni);
  free(range);
  free(xlabel);  free(ylabel);
  /* free(figname);  free(_figname); */

  free(_xpt);  free(xpt);
  free(_ypt);  free(ypt);
  free(_zpt);  free(zpt);
  free(pt);

  free(xylegTxt);
  free(xzlegTxt);
  free(yzlegTxt);

  free(num);


  __NOTE__("%s\n", "end");
}
