/**
 * @file trace.c
 *
 * @brief Source code for test-particle simulations
 *
 * @author Yohei MIKI (The University of Tokyo)
 *
 * @date 2019/06/15 (Sat)
 *
 * Copyright (C) 2019 Yohei MIKI
 * All rights reserved.
 *
 * The MIT License is applied to this software, see LICENSE.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>

#ifdef  USE_HDF5_FORMAT
#include <hdf5.h>
#include "hdf5lib.h"
#endif//USE_HDF5_FORMAT

#include "macro.h"
#include "myutil.h"
#include "name.h"
#include "constants.h"
#include "timer.h"
#include "cudalib.h"
#include "mpilib.h"

#include "../misc/device.h"
#include "../misc/structure.h"
#include "../misc/allocate.h"
#include "../misc/allocate_dev.h"

#include "../file/io.h"


int main(int argc, char **argv)
{
  /** parallelized region employing MPI start */
  static MPIinfo mpi;
  initMPI(&mpi, &argc, &argv);
  static MPIcfg_dataio iocfg;
  createMPIcfg_dataio(&iocfg, mpi);


  /** configure the details of the numerical simulation */
  /** read command line arguments */
  if( argc < 3 ){
    __FPRINTF__(stderr, "insufficient number of input parameters of %d (at least %d inputs are required).\n", argc, 3);
    __FPRINTF__(stderr, "Usage is: %s\n", argv[0]);
    __FPRINTF__(stderr, "          -file=<char *>\n");
    __FPRINTF__(stderr, "          -pot_file_sphe=<char *>\n");
    __KILL__(stderr, "%s\n", "insufficient command line arguments");
  }
  char *file;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv,  "file", &file));
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  char *pot_file_sphe;  requiredCmdArg(getCmdArgStr(argc, (const char * const *)(void *)argv,  "pot_file_sphe", &pot_file_sphe));
#endif//SET_EXTERNAL_POTENTIAL_FIELD


  /** read settings about the simulation */
  static ulong Ntot;
  static real eps, eta;
  static double ft, snapshotInterval, saveInterval;
  static int unit;
  readSettingsParallel(&unit, &Ntot, &eps, &eta, &ft, &snapshotInterval, &saveInterval, file, mpi);

  /** read setting to dump tentative results of the simulation */
  static int last;
  readConfigFileParallel(&last, file, mpi);


  /** activate the attatched accelerator device(s) */
  int devIdx;
  deviceInfo devInfo;
  deviceProp devProp;
  int headIdx;
  if( optionalCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "deviceID", &headIdx)) != myUtilAvail )
    headIdx = 0;
  openAcceleratorGPUs(&devIdx, &devInfo, &devProp, headIdx, 1, mpi.rank, mpi.size);
  __NOTE__("devIdx = %d\n", devIdx);


  /** set the number of N-body particles */
  const ulong Nini = (Ntot * (    (ulong)mpi.rank)) / (ulong)mpi.size;
  const ulong Nfin = (Ntot * (1 + (ulong)mpi.rank)) / (ulong)mpi.size;
  const int num = (int)(Nfin - Nini);


  /** set global constants */
  /** set gravitational softening and opening criterion */
  extern real newton;
  hoge;


  /** memory allocation */
  /** declaration of array to contain whole information of whole N-body particles */
#ifdef  USE_HDF5_FORMAT
  nbody_hdf5 body_snapshot;
  real *hdf5_pos, *hdf5_vel, *hdf5_acc, *hdf5_m, *hdf5_pot;
  ulong *hdf5_idx;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  real *hdf5_acc_ext, *hdf5_pot_ext;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
  const muse alloc_snap = allocSnapshotArray
    (&hdf5_pos, &hdf5_vel, &hdf5_acc, &hdf5_m, &hdf5_pot, &hdf5_idx,
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     &hdf5_acc_ext, &hdf5_pot_ext,
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     num, &body_snapshot);
#endif//USE_HDF5_FORMAT

  /** declarations of arrays to contain information of i-particles */
  /** allocParticleDataSoa_dev is defined in src/tree/walk_dev.cu --> implement a similar function in another file */
  iparticle     ibody0,  ibody1,  ibody0_dev,  ibody1_dev;
  ulong        *  idx0, *  idx1, *  idx0_dev, *  idx1_dev;
  position     *  pos0, *  pos1, *  pos0_dev, *  pos1_dev;
  acceleration *  acc0, *  acc1, *  acc0_dev, *  acc1_dev;
  real *neighbor_dev;
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
  acceleration *acc_ext0, *acc_ext1, *acc_ext_dev;
#endif//SET_EXTERNAL_POTENTIAL_FIELD
#ifdef  BLOCK_TIME_STEP
  velocity     * vel0, * vel1, * vel0_dev, * vel1_dev;
  ibody_time   *time0, *time1, *time0_dev, *time1_dev;
  const muse alloc_ibody0    = allocParticleDataSoA_hst
    (num, &ibody0, &idx0, &pos0, &acc0, &vel0, &time0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext0
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  const muse alloc_ibody1    = allocParticleDataSoA_hst
    (num, &ibody1, &idx1, &pos1, &acc1, &vel1, &time1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext1
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  const muse alloc_ibody_dev = allocParticleDataSoA_dev
    (num
     , &ibody0_dev, &idx0_dev, &pos0_dev, &acc0_dev, &vel0_dev, &time0_dev
     , &ibody1_dev, &idx1_dev, &pos1_dev, &acc1_dev, &vel1_dev, &time1_dev
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext_dev
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     , &neighbor_dev
     );
#else///BLOCK_TIME_STEP
  real *vx0, *vx1, *vx0_dev, *vx1_dev;
  real *vy0, *vy1, *vy0_dev, *vy1_dev;
  real *vz0, *vz1, *vz0_dev, *vz1_dev;
  const muse alloc_ibody0    = allocParticleDataSoA_hst
    (num, &ibody0, &idx0, &pos0, &acc0, &vx0, &vy0, &vz0
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext0
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  const muse alloc_ibody1    = allocParticleDataSoA_hst
    (num, &ibody1, &idx1, &pos1, &acc1, &vx1, &vy1, &vz1
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext1
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     );
  const muse alloc_ibody_dev = allocParticleDataSoA_dev
    (num
     , &ibody0_dev, &idx0_dev, &pos0_dev, &acc0_dev, &vx0_dev, &vy0_dev, &vz0_dev
     , &ibody1_dev, &idx1_dev, &pos1_dev, &acc1_dev, &vx1_dev, &vy1_dev, &vz1_dev
#ifdef  SET_EXTERNAL_POTENTIAL_FIELD
     , &acc_ext_dev
#endif//SET_EXTERNAL_POTENTIAL_FIELD
     , &neighbor_dev
     );
#endif//BLOCK_TIME_STEP




  /** declarations of arrays to contain information of sink particle(s) */
  /**# remove the sink particle from the initial-condition; however, read the mass of the central BH to determine the loss cone */
  /* test-particles should memory the time if the particle accretes the central BH */
  static int Nsink = 0;
  requiredCmdArg(getCmdArgInt(argc, (const char * const *)(void *)argv, "Nsink", &Nsink));

  static sinkparticle sink;
  position *pos_sink;
  velocity *vel_sink, *mom_sink;
  ulong *tag_sink, *tag_sink_hst;
  int *list_sink;
  real *lmax2_sink;
  const muse alloc_sinkparticle_dev =
    allocSinkParticleSoA_dev(&pos_sink, &vel_sink, &mom_sink, &tag_sink, &list_sink, &lmax2_sink, &tag_sink_hst, Nsink, &sink);

  read_sink_particle_indexes(Nsink, tag_sink_hst, file);
  set_sink_particle_dev(Nsink, sink, tag_sink_hst);


  /** read initial condition */
#   if  defined(MONITOR_ENERGY_ERROR) && defined(USE_HDF5_FORMAT)
  static energyError relEneErr;
#endif//defined(MONITOR_ENERGY_ERROR) && defined(USE_HDF5_FORMAT)
  static double time, dt;
  static ulong steps = 0;
  static double prevElapsed = 0.0;
#ifdef  USE_HDF5_FORMAT
  static hdf5struct hdf5type;
  createHDF5DataType(&hdf5type);
#endif//USE_HDF5_FORMAT
  readTentativeDataParallel
    (&time, &dt, &steps, &prevElapsed, &num, ibody0, file, last, &iocfg, Ntot
#ifdef  USE_HDF5_FORMAT
     , hdf5type, &dropPrevTune, &rebuild, &elapsed, &rebuildParam, &brentDistance, &brentHistory
#ifdef  MONITOR_ENERGY_ERROR
     , &relEneErr
#endif//MONITOR_ENERGY_ERROR
#endif//USE_HDF5_FORMAT
     );
#ifndef BLOCK_TIME_STEP
  real *dt_dev;
  const muse alloc_dt_dev = allocTimeStep_dev(&dt_dev);
#endif//BLOCK_TIME_STEP

  /** set up for output files */
  static double formerTime;
  if( mpi.rank == 0 )
    formerTime = getPresentTimeInMin();
  chkMPIerr(MPI_Bcast(&formerTime, 1, MPI_DOUBLE, 0, mpi.comm));
  static double invSnapshotInterval;
  invSnapshotInterval = 1.0 / snapshotInterval;
  static uint previous;
  previous = (uint)(time * invSnapshotInterval);


  /** read numeric table for fixed potential field */
  potential_field pot_tbl_sphe;
  pot2 *pot_tbl_sphe_Phi;
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  real *pot_tbl_sphe_rad;
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
  const muse alloc_ext_pot_sphe = readFixedPotentialTableSpherical
    (unit, pot_file_sphe, &pot_tbl_sphe, &pot_tbl_sphe_Phi
#ifdef  ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
     , &pot_tbl_sphe_rad
#endif//ADAPTIVE_GRIDDED_EXTERNAL_POTENTIAL_FIELD
#ifdef  USE_HDF5_FORMAT
     , hdf5type
#endif//USE_HDF5_FORMAT
     );










  exitMPI();
}
