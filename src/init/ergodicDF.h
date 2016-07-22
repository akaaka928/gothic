/*************************************************************************\
 *                                                                       *
                  last updated on 2014/08/04(Mon) 14:30:07
 *                                                                       *
 *    Header File for Definition about ergodic distribution functions    *
 *                                                                       *
 *                                                                       *
 *                                             written by Yohei MIKI     *
 *                                                                       *
\*************************************************************************/
//-------------------------------------------------------------------------
#ifndef ERGODICDF_H
#define ERGODICDF_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#ifndef MACRO_H
#      include <macro.h>
#endif//MACRO_H
//-------------------------------------------------------------------------
#ifndef STRUCTURE_H
#       include "../misc/structure.h"
#endif//STRUCTURE_H
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
/* Plummer sphere */
void rescalePlummerSphere(real Mtot, real h, int *num, real **rad, real **rho);
void calcPlummerDynamicalProperties(real Mtot, real h, int num, real *rad, real **vdisp, real **phi);
void calcPlummerObservableProperties(int num, real sigma, real *rad, real *rho, real **Sigma, real **sigmalos);
void outputFundamentalInformationOfPlummerSphere(const real Mtot, const real h, const int num, const ulong Ntot, const real eps, const real snapshotInterval, const real ft, real *rad, real *rho, real *vdisp, real *phi, char file[]);
//-------------------------------------------------------------------------
/* King sphere */
void makeKingDFTable(real W0, int *num, real **rad, real **W, real **rho);
void rescaleKingSphere(real Mtot, real rt, real *r0, real *sigma, int num, real *rad, real *rho, real **Menc);
void calcKingDynamicalProperties(real Mtot, real rt, real sigma, int num, real *W, real **phi, real **psi, real **vdisp);
void calcKingObservableProperties(int num, real sigma, real *rad, real *rho, real *W, real **Sigma, real **sigmalos);
void outputFundamentalInformationOfKingSphere(const real W0, const real Mtot, const real r0, const real rt, const real sigma, const int num, const ulong Ntot, const real eps, const real snapshotInterval, const real ft, real *rad, real *rho, real *vdisp, real *Sigma, real *vdlos, real *Menc, real *phi, real *psi, real *W, char file[]);
void makeKingSphereByNbody(ulong num, nbody_particle *body, real Mtot, real sigma, int Nking, real *rad, real *W, real *Menc);
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
#endif//ERGODICDF_H
//-------------------------------------------------------------------------
