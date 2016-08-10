/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                            avtGOTHICFileFormat.C                           //
// ************************************************************************* //

#include <avtGOTHICFileFormat.h>

#include <string>

#include <vtkCellType.h>
#include <vtkUnsignedLongArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>

#include <avtDatabaseMetaData.h>

#include <DBOptionsAttributes.h>
#include <Expression.h>

#include <InvalidVariableException.h>

// to use ostringstream
#include <iostream>

//-------------------------------------------------------------------------
// HDF5 related original routines
//-------------------------------------------------------------------------
#include <hdf5.h>
//-------------------------------------------------------------------------
#define chkHDF5err(err) __chkHDF5err(err, __FILE__, __LINE__)
//-------------------------------------------------------------------------
#define __FPRINTF__(dst, ...)				\
  {								\
    fprintf(dst, "%s(%d): %s\n", __FILE__, __LINE__, __func__);	\
    fprintf(dst, ##__VA_ARGS__);				\
    fflush(NULL);						\
  }
//-------------------------------------------------------------------------
#if defined(MPI_INCLUDED) || defined(OMPI_MPI_H)
#     define __KILL__(dst, ...)		\
  {						\
    __FPRINTF__(dst, ##__VA_ARGS__);		\
    MPI_Abort(MPI_COMM_WORLD, 1);		\
  }
#else
#     define __KILL__(dst, ...)		\
  {						\
    __FPRINTF__(dst, ##__VA_ARGS__);		\
    exit(EXIT_FAILURE);				\
  }
#endif
//-------------------------------------------------------------------------
void __chkHDF5err(herr_t err, const char *file, const int line)
{
  //-----------------------------------------------------------------------
  /* Returns a non-negative value if successful; otherwise returns a negative value. */
  if( err < 0 ){
    fprintf(stderr, "Error is detected at %s(%d)\n", file, line);
    __KILL__(stderr, "HDF5 error, error ID is %d\n", err);
  }/* if( err < 0 ){ */
  //-----------------------------------------------------------------------
}
//-------------------------------------------------------------------------
typedef unsigned long int ulong;
//-------------------------------------------------------------------------
// typedef struct
// {
//   ulong idx;
//   double t0, t1;
//   float dt;
//   float  x,  y,  z;
//   float vx, vy, vz;
//   float ax, ay, az;
//   float m, pot;
// } nbodyBlockSP;
//-------------------------------------------------------------------------
// typedef struct
// {
//   ulong idx;
//   double t0, t1;
//   double dt;
//   double  x,  y,  z;
//   double vx, vy, vz;
//   double ax, ay, az;
//   double m, pot;
// } nbodyBlockDP;
//-------------------------------------------------------------------------
// typedef struct
// {
//   ulong idx;
//   float  x,  y,  z;
//   float vx, vy, vz;
//   float ax, ay, az;
//   float m, pot;
// } nbodyShareSP;
//-------------------------------------------------------------------------
// typedef struct
// {
//   ulong idx;
//   double  x,  y,  z;
//   double vx, vy, vz;
//   double ax, ay, az;
//   double m, pot;
// } nbodyShareDP;
//-------------------------------------------------------------------------
// typedef struct
// {
//   hid_t blockSP, blockDP;
//   hid_t shareSP, shareDP;
// } hdf5struct;
//-------------------------------------------------------------------------
// void createHDF5DataType(hdf5struct *type)
// {
//   //-----------------------------------------------------------------------
//   /* commit data type of nbodyBlockSP */
//   type->blockSP = H5Tcreate(H5T_COMPOUND, sizeof(nbodyBlockSP));
//   chkHDF5err(H5Tinsert(type->blockSP,  "index", HOFFSET(nbodyBlockSP, idx), H5T_NATIVE_ULONG));
//   chkHDF5err(H5Tinsert(type->blockSP,   "time", HOFFSET(nbodyBlockSP,  t0), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockSP, "t + dt", HOFFSET(nbodyBlockSP,  t1), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockSP,     "dt", HOFFSET(nbodyBlockSP,  dt), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,      "x", HOFFSET(nbodyBlockSP,   x), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,      "y", HOFFSET(nbodyBlockSP,   y), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,      "z", HOFFSET(nbodyBlockSP,   z), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,     "vx", HOFFSET(nbodyBlockSP,  vx), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,     "vy", HOFFSET(nbodyBlockSP,  vy), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,     "vz", HOFFSET(nbodyBlockSP,  vz), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,     "ax", HOFFSET(nbodyBlockSP,  ax), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,     "ay", HOFFSET(nbodyBlockSP,  ay), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,     "az", HOFFSET(nbodyBlockSP,  az), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,      "m", HOFFSET(nbodyBlockSP,   m), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->blockSP,    "pot", HOFFSET(nbodyBlockSP, pot), H5T_NATIVE_FLOAT));
//   //-----------------------------------------------------------------------
//   /* commit data type of nbodyBlockDP */
//   type->blockDP = H5Tcreate(H5T_COMPOUND, sizeof(nbodyBlockDP));
//   chkHDF5err(H5Tinsert(type->blockDP,  "index", HOFFSET(nbodyBlockDP, idx), H5T_NATIVE_ULONG));
//   chkHDF5err(H5Tinsert(type->blockDP,   "time", HOFFSET(nbodyBlockDP,  t0), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP, "t + dt", HOFFSET(nbodyBlockDP,  t1), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,     "dt", HOFFSET(nbodyBlockDP,  dt), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,      "x", HOFFSET(nbodyBlockDP,   x), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,      "y", HOFFSET(nbodyBlockDP,   y), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,      "z", HOFFSET(nbodyBlockDP,   z), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,     "vx", HOFFSET(nbodyBlockDP,  vx), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,     "vy", HOFFSET(nbodyBlockDP,  vy), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,     "vz", HOFFSET(nbodyBlockDP,  vz), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,     "ax", HOFFSET(nbodyBlockDP,  ax), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,     "ay", HOFFSET(nbodyBlockDP,  ay), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,     "az", HOFFSET(nbodyBlockDP,  az), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,      "m", HOFFSET(nbodyBlockDP,   m), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->blockDP,    "pot", HOFFSET(nbodyBlockDP, pot), H5T_NATIVE_DOUBLE));
//   //-----------------------------------------------------------------------
//   /* commit data type of nbodyShareSP */
//   type->shareSP = H5Tcreate(H5T_COMPOUND, sizeof(nbodyShareSP));
//   chkHDF5err(H5Tinsert(type->shareSP, "index", HOFFSET(nbodyShareSP, idx), H5T_NATIVE_ULONG));
//   chkHDF5err(H5Tinsert(type->shareSP,     "x", HOFFSET(nbodyShareSP,   x), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->shareSP,     "y", HOFFSET(nbodyShareSP,   y), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->shareSP,     "z", HOFFSET(nbodyShareSP,   z), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->shareSP,    "vx", HOFFSET(nbodyShareSP,  vx), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->shareSP,    "vy", HOFFSET(nbodyShareSP,  vy), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->shareSP,    "vz", HOFFSET(nbodyShareSP,  vz), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->shareSP,    "ax", HOFFSET(nbodyShareSP,  ax), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->shareSP,    "ay", HOFFSET(nbodyShareSP,  ay), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->shareSP,    "az", HOFFSET(nbodyShareSP,  az), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->shareSP,     "m", HOFFSET(nbodyShareSP,   m), H5T_NATIVE_FLOAT));
//   chkHDF5err(H5Tinsert(type->shareSP,   "pot", HOFFSET(nbodyShareSP, pot), H5T_NATIVE_FLOAT));
//   //-----------------------------------------------------------------------
//   /* commit data type of nbodyShareDP */
//   type->shareDP = H5Tcreate(H5T_COMPOUND, sizeof(nbodyShareDP));
//   chkHDF5err(H5Tinsert(type->shareDP, "index", HOFFSET(nbodyShareDP, idx), H5T_NATIVE_ULONG));
//   chkHDF5err(H5Tinsert(type->shareDP,     "x", HOFFSET(nbodyShareDP,   x), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->shareDP,     "y", HOFFSET(nbodyShareDP,   y), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->shareDP,     "z", HOFFSET(nbodyShareDP,   z), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->shareDP,    "vx", HOFFSET(nbodyShareDP,  vx), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->shareDP,    "vy", HOFFSET(nbodyShareDP,  vy), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->shareDP,    "vz", HOFFSET(nbodyShareDP,  vz), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->shareDP,    "ax", HOFFSET(nbodyShareDP,  ax), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->shareDP,    "ay", HOFFSET(nbodyShareDP,  ay), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->shareDP,    "az", HOFFSET(nbodyShareDP,  az), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->shareDP,     "m", HOFFSET(nbodyShareDP,   m), H5T_NATIVE_DOUBLE));
//   chkHDF5err(H5Tinsert(type->shareDP,   "pot", HOFFSET(nbodyShareDP, pot), H5T_NATIVE_DOUBLE));
//   //-----------------------------------------------------------------------
// }
//-------------------------------------------------------------------------
// void removeHDF5DataType(hdf5struct  type)
// {
//   //-----------------------------------------------------------------------
//   chkHDF5err(H5Tclose(type.blockSP));
//   chkHDF5err(H5Tclose(type.blockDP));
//   chkHDF5err(H5Tclose(type.shareSP));
//   chkHDF5err(H5Tclose(type.shareDP));
//   //-----------------------------------------------------------------------
// }
//-------------------------------------------------------------------------


using     std::string;


// ****************************************************************************
//  Method: avtGOTHICFileFormat constructor
//
//  Programmer: ymiki -- generated by xml2avt
//  Creation:   Thu Feb 4 11:45:33 PDT 2016
//
// ****************************************************************************

avtGOTHICFileFormat::avtGOTHICFileFormat(const char *filename)
    : avtSTSDFileFormat(filename)
{
    // INITIALIZE DATA MEMBERS
}


// ****************************************************************************
//  Method: avtGOTHICFileFormat::FreeUpResources
//
//  Purpose:
//      When VisIt is done focusing on a particular timestep, it asks that
//      timestep to free up any resources (memory, file descriptors) that
//      it has associated with it.  This method is the mechanism for doing
//      that.
//
//  Programmer: ymiki -- generated by xml2avt
//  Creation:   Thu Feb 4 11:45:33 PDT 2016
//
// ****************************************************************************

void
avtGOTHICFileFormat::FreeUpResources(void)
{
}


// ****************************************************************************
//  Method: avtGOTHICFileFormat::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: ymiki -- generated by xml2avt
//  Creation:   Thu Feb 4 11:45:33 PDT 2016
//
// ****************************************************************************

void
avtGOTHICFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md)
{
  string meshname = "nbody";
  // add point mesh
  avtMeshMetaData *mmd = new avtMeshMetaData;
  mmd->name = meshname;
  mmd->spatialDimension = 3;
  mmd->topologicalDimension = 0;
  mmd->meshType = AVT_POINT_MESH;
  mmd->numBlocks = 1;  // <-- this must be 1 for STSD
  mmd->xLabel = "x";
  mmd->yLabel = "y";
  mmd->zLabel = "z";
  md->Add(mmd);

  // add physical variables
  std::ostringstream var;
  // add scalar variables
  avtScalarMetaData *pot = new avtScalarMetaData;  pot->meshName = meshname;  pot->centering = AVT_NODECENT;  pot->hasUnits = false;
  avtScalarMetaData *idx = new avtScalarMetaData;  idx->meshName = meshname;  idx->centering = AVT_NODECENT;  idx->hasUnits = false;
  var << meshname << "/potential";  pot->name = var.str();  md->Add(pot);  var.str("");
  var << meshname << "/index"    ;  idx->name = var.str();  md->Add(idx);  var.str("");
  // add vector variables
  avtVectorMetaData *acc = new avtVectorMetaData;  acc->meshName = meshname;  acc->centering = AVT_NODECENT;  acc->hasUnits = false;  acc->varDim = 3;
  avtVectorMetaData *vel = new avtVectorMetaData;  vel->meshName = meshname;  vel->centering = AVT_NODECENT;  vel->hasUnits = false;  vel->varDim = 3;
  var << meshname << "/acceleration";  acc->name = var.str();  md->Add(acc);  var.str("");
  var << meshname << "/velocity"    ;  vel->name = var.str();  md->Add(vel);  var.str("");
}


// ****************************************************************************
//  Method: avtGOTHICFileFormat::GetMesh
//
//  Purpose:
//      Gets the mesh associated with this file.  The mesh is returned as a
//      derived type of vtkDataSet (ie vtkRectilinearGrid, vtkStructuredGrid,
//      vtkUnstructuredGrid, etc).
//
//  Arguments:
//      meshname    The name of the mesh of interest.  This can be ignored if
//                  there is only one mesh.
//
//  Programmer: ymiki -- generated by xml2avt
//  Creation:   Thu Feb 4 11:45:33 PDT 2016
//
// ****************************************************************************

vtkDataSet *
avtGOTHICFileFormat::GetMesh(const char *meshname)
{
  if( strcmp("nbody", meshname) == 0 ){
    // open the target file
    const char* filename = GetFilename();
    hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    // read useDP
    int useDP;
    hid_t group = H5Gopen(target, meshname, H5P_DEFAULT);
    hid_t attribute = H5Aopen(group, "useDP", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
    chkHDF5err(H5Aclose(attribute));

    const int ndims = 3;
    ulong num;
    // read number of nodes from the file (written as an attribute)
    attribute = H5Aopen(group, "number", H5P_DEFAULT);
    chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, &num));
    chkHDF5err(H5Aclose(attribute));
    int nnodes = (int)num;

    hid_t dataset = H5Dopen(group, "position", H5P_DEFAULT);
    if( useDP == 0 ){
      // read the XYZ coordinates from the file
      float *position = new float[nnodes * 4];
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, position));

      // create the vtkPoints object and copy points into it
      vtkPoints *points = vtkPoints::New();
      points->SetNumberOfPoints(nnodes);
      float *pts = (float *)points->GetVoidPointer(0);
      for(int ii = 0; ii < nnodes; ii++){
	pts[    3 * ii] = position[    4 * ii];
	pts[1 + 3 * ii] = position[1 + 4 * ii];
	pts[2 + 3 * ii] = position[2 + 4 * ii];
      }// for(int ii = 0; ii < nnodes; ii++){
      delete [] position;

      // close the target file
      chkHDF5err(H5Dclose(dataset));
      chkHDF5err(H5Gclose(group));
      chkHDF5err(H5Fclose(target));

      // create a vtkUnstructuredGrid to contain the point cells
      vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
      ugrid->SetPoints(points);
      points->Delete();
      ugrid->Allocate(nnodes);
      vtkIdType onevertex;
      for(int ii = 0; ii < nnodes; ii++){
	onevertex = ii;
	ugrid->InsertNextCell(VTK_VERTEX, 1, &onevertex);
      }// for(int ii = 0; ii < nnodes; ii++){

      return ugrid;
    }// if( useDP == 0 ){
    else{
      // read the XYZ coordinates from the file
      double *position = new double[nnodes * 4];
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, position));

      // create the vtkPoints object and copy points into it
      vtkPoints *points = vtkPoints::New();
      points->SetNumberOfPoints(nnodes);
      double *pts = (double *)points->GetVoidPointer(0);
      for(int ii = 0; ii < nnodes; ii++){
	pts[    3 * ii] = position[    4 * ii];
	pts[1 + 3 * ii] = position[1 + 4 * ii];
	pts[2 + 3 * ii] = position[2 + 4 * ii];
      }// for(int ii = 0; ii < nnodes; ii++){
      delete [] position;

      // close the target file
      chkHDF5err(H5Dclose(dataset));
      chkHDF5err(H5Gclose(group));
      chkHDF5err(H5Fclose(target));

      // create a vtkUnstructuredGrid to contain the point cells
      vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
      ugrid->SetPoints(points);
      points->Delete();
      ugrid->Allocate(nnodes);
      vtkIdType onevertex;
      for(int ii = 0; ii < nnodes; ii++){
	onevertex = ii;
	ugrid->InsertNextCell(VTK_VERTEX, 1, &onevertex);
      }// for(int ii = 0; ii < nnodes; ii++){

      return ugrid;
    }// else{
  }// if( strncmp("nbody", meshname, 4) == 0 ){
  else{
    // Error exception
    EXCEPTION1(InvalidVariableException, meshname);
  }// else{
}


// ****************************************************************************
//  Method: avtGOTHICFileFormat::GetVar
//
//  Purpose:
//      Gets a scalar variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      varname    The name of the variable requested.
//
//  Programmer: ymiki -- generated by xml2avt
//  Creation:   Thu Feb 4 11:45:33 PDT 2016
//
// ****************************************************************************

vtkDataArray *
avtGOTHICFileFormat::GetVar(const char *varname)
{
  // extract meshname and variable from varname (varname = meshname/variable)
  const string str = varname;
  string::size_type pos = str.find_last_of("/");
  string meshname = str.substr(0, pos);
  string variable = str.substr(pos + 1, str.size() - pos);

  // open the target file
  const char* filename = GetFilename();
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  ulong ulnum;
  // read number of particles from the file (written as an attribute)
  hid_t group = H5Gopen(target, meshname.c_str(), H5P_DEFAULT);
  hid_t attribute = H5Aopen(group, "number", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, &ulnum));
  chkHDF5err(H5Aclose(attribute));
  int num = (int)ulnum;

  // read useDP
  int useDP;
  attribute = H5Aopen(group, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));

  if( strcmp(variable.c_str(), "index") != 0 ){
    if( useDP == 0 ){
      // read potential
      float *array = new float[num * 4];
      hid_t dataset = H5Dopen(group, "acceleration", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, array));
      chkHDF5err(H5Dclose(dataset));
      chkHDF5err(H5Gclose(group));

      // close the target file
      chkHDF5err(H5Fclose(target));

      vtkFloatArray *rv = vtkFloatArray::New();
      rv->SetNumberOfTuples(num);
      for(int ii = 0; ii < num; ii++)
	rv->SetTuple1(ii, array[3 + 4 * ii]);

      // delete temporary arrays
      delete [] array;

      return rv;
    }// if( useDP == 0 ){
    else{
      // read potential
      double *array = new double[num * 4];
      hid_t dataset = H5Dopen(group, "acceleration", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array));
      chkHDF5err(H5Dclose(dataset));
      chkHDF5err(H5Gclose(group));

      // close the target file
      chkHDF5err(H5Fclose(target));

      vtkDoubleArray *rv = vtkDoubleArray::New();
      rv->SetNumberOfTuples(num);
      for(int ii = 0; ii < num; ii++)
	rv->SetTuple1(ii, array[3 + 4 * ii]);

      // delete temporary arrays
      delete [] array;

      return rv;
    }// else{
  }// if( strcmp(variable.c_str(), "index") != 0 ){
  else{
    // read the target variable (index)
    ulong *array = new ulong[num];
    hid_t dataset = H5Dopen(group, variable.c_str(), H5P_DEFAULT);
    chkHDF5err(H5Dread(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, array));
    chkHDF5err(H5Dclose(dataset));

    // close the target file
    chkHDF5err(H5Gclose(group));
    chkHDF5err(H5Fclose(target));

    vtkUnsignedLongArray *rv = vtkUnsignedLongArray::New();
    rv->SetNumberOfTuples(num);
    for(int ii = 0; ii < num; ii++)
      rv->SetTuple1(ii, array[ii]);

    // delete temporary arrays
    delete [] array;

    return rv;
  }// else{
}


// ****************************************************************************
//  Method: avtGOTHICFileFormat::GetVectorVar
//
//  Purpose:
//      Gets a vector variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      varname    The name of the variable requested.
//
//  Programmer: ymiki -- generated by xml2avt
//  Creation:   Thu Feb 4 11:45:33 PDT 2016
//
// ****************************************************************************

vtkDataArray *
avtGOTHICFileFormat::GetVectorVar(const char *varname)
{
  // extract meshname and variable from varname (varname = meshname/variable)
  const string str = varname;
  string::size_type pos = str.find_last_of("/");
  string meshname = str.substr(0, pos);
  string variable = str.substr(pos + 1, str.size() - pos);

  // open the target file
  const char* filename = GetFilename();
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  ulong ulnum;
  // read number of particles from the file (written as an attribute)
  hid_t group = H5Gopen(target, meshname.c_str(), H5P_DEFAULT);
  hid_t attribute = H5Aopen(group, "number", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, &ulnum));
  chkHDF5err(H5Aclose(attribute));
  int num = (int)ulnum;
  const int ndims = 3;

  // read useDP
  int useDP;
  attribute = H5Aopen(group, "useDP", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &useDP));
  chkHDF5err(H5Aclose(attribute));

  // read blockTimeStep
  int blockTimeStep;
  attribute = H5Aopen(group, "blockTimeStep", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_INT, &blockTimeStep));
  chkHDF5err(H5Aclose(attribute));

  if( (blockTimeStep == 1) || (strcmp(variable.c_str(), "velocity") != 0) ){
    if( useDP == 0 ){
      // read the target variable
      float *array = new float[num * 4];
      hid_t dataset = H5Dopen(group, variable.c_str(), H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, array));
      chkHDF5err(H5Dclose(dataset));
      chkHDF5err(H5Gclose(group));
      chkHDF5err(H5Fclose(target));

      vtkFloatArray *rv = vtkFloatArray::New();
      rv->SetNumberOfComponents(ndims);
      rv->SetNumberOfTuples(num);
      float *dat = (float *)rv->GetVoidPointer(0);
      for(int ii = 0; ii < num; ii++){
	dat[    3 * ii] = array[    4 * ii];
	dat[1 + 3 * ii] = array[1 + 4 * ii];
	dat[2 + 3 * ii] = array[2 + 4 * ii];
      }// for(int ii = 0; ii < num; ii++){

      // delete temporary arrays
      delete [] array;

      return rv;
    }// if( useDP == 0 ){
    else{
      // read the target variable
      double *array = new double[num * 4];
      hid_t dataset = H5Dopen(group, variable.c_str(), H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array));
      chkHDF5err(H5Dclose(dataset));
      chkHDF5err(H5Gclose(group));
      chkHDF5err(H5Fclose(target));

      vtkDoubleArray *rv = vtkDoubleArray::New();
      rv->SetNumberOfComponents(ndims);
      rv->SetNumberOfTuples(num);
      double *dat = (double *)rv->GetVoidPointer(0);
      for(int ii = 0; ii < num; ii++){
	dat[    3 * ii] = array[    4 * ii];
	dat[1 + 3 * ii] = array[1 + 4 * ii];
	dat[2 + 3 * ii] = array[2 + 4 * ii];
      }// for(int ii = 0; ii < num; ii++){

      // delete temporary arrays
      delete [] array;

      return rv;
    }// else{
  }// if( (blockTimeStep == 1) || (strcmp(variable.c_str(), "velocity") != 0) ){
  else{
    if( useDP == 0 ){
      // read the target variable
      float *vx = new float[num];
      float *vy = new float[num];
      float *vz = new float[num];
      hid_t dataset = H5Dopen(group, "vx", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vx));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dopen(group, "vy", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vy));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dopen(group, "vz", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vz));
      chkHDF5err(H5Dclose(dataset));
      chkHDF5err(H5Gclose(group));
      chkHDF5err(H5Fclose(target));

      vtkFloatArray *rv = vtkFloatArray::New();
      rv->SetNumberOfComponents(ndims);
      rv->SetNumberOfTuples(num);
      float *dat = (float *)rv->GetVoidPointer(0);
      for(int ii = 0; ii < num; ii++){
	dat[    3 * ii] = vx[ii];
	dat[1 + 3 * ii] = vy[ii];
	dat[2 + 3 * ii] = vz[ii];
      }// for(int ii = 0; ii < num; ii++){

      // delete temporary arrays
      delete [] vx;
      delete [] vy;
      delete [] vz;

      return rv;
    }// if( useDP == 0 ){
    else{
      // read the target variable
      double *vx = new double[num];
      double *vy = new double[num];
      double *vz = new double[num];
      hid_t dataset = H5Dopen(group, "vx", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vx));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dopen(group, "vy", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vy));
      chkHDF5err(H5Dclose(dataset));
      dataset = H5Dopen(group, "vz", H5P_DEFAULT);
      chkHDF5err(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vz));
      chkHDF5err(H5Dclose(dataset));
      chkHDF5err(H5Gclose(group));
      chkHDF5err(H5Fclose(target));

      vtkDoubleArray *rv = vtkDoubleArray::New();
      rv->SetNumberOfComponents(ndims);
      rv->SetNumberOfTuples(num);
      double *dat = (double *)rv->GetVoidPointer(0);
      for(int ii = 0; ii < num; ii++){
	dat[    3 * ii] = vx[ii];
	dat[1 + 3 * ii] = vy[ii];
	dat[2 + 3 * ii] = vz[ii];
      }// for(int ii = 0; ii < num; ii++){

      // delete temporary arrays
      delete [] vx;
      delete [] vy;
      delete [] vz;

      return rv;
    }// else{
  }// else{
}


// optional function(s)
double
avtGOTHICFileFormat::GetTime(void)
{
  // open the target file
  const char* filename = GetFilename();
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  // read current time contained in the file (written as an attribute)
  double dtime;
  hid_t group = H5Gopen(target, "nbody", H5P_DEFAULT);
  hid_t attribute = H5Aopen(group, "time", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_DOUBLE, &dtime));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));
  chkHDF5err(H5Fclose(target));

  return (dtime);
}
int
avtGOTHICFileFormat::GetCycle(void)
{
  // open the target file
  const char* filename = GetFilename();
  hid_t target = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  // read current time contained in the file (written as an attribute)
  ulong steps;
  hid_t group = H5Gopen(target, "nbody", H5P_DEFAULT);
  hid_t attribute = H5Aopen(group, "steps", H5P_DEFAULT);
  chkHDF5err(H5Aread(attribute, H5T_NATIVE_ULONG, &steps));
  chkHDF5err(H5Aclose(attribute));
  chkHDF5err(H5Gclose(group));

  // close the target file
  chkHDF5err(H5Fclose(target));

  int ret = (int)steps;
  return (ret);
}
