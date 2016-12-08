#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macro.h"
#include "spline.h"
#include "spline.c"

#define NDATA (9)
#define NPLOT (256)

int main(void){

  static double xd[NDATA], yd[NDATA];
  xd[0] = 0.0;  yd[0] = 0.0;
  xd[1] = 0.1;  yd[1] = 0.0;
  xd[2] = 0.2;  yd[2] = 0.0;
  xd[3] = 0.4;  yd[3] = 0.2;
  xd[4] = 0.5;  yd[4] = 1.0;
  xd[5] = 0.6;  yd[5] = 0.2;
  xd[6] = 0.8;  yd[6] = 0.0;
  xd[7] = 0.9;  yd[7] = 0.0;
  xd[8] = 1.0;  yd[8] = 0.0;
  static double bp[NDATA], y2[NDATA];
#if 1
  genCubicSpline1D(NDATA, xd, yd, bp, 0.0, 0.0, y2);
#else
  genCubicSpline1D(NDATA, xd, yd, bp, NATURAL_CUBIC_SPLINE, NATURAL_CUBIC_SPLINE, y2);
#endif

  const double xbin = 1.0 / (double)(NPLOT - 1);
  for(int ii = 0; ii < NPLOT; ii++){
    const double xp = xbin * (double)ii;
    printf("%e\t%e\t%e\t%e\n", xp,
	   getCubicSpline1D               (xp, NDATA, xd, yd, y2),
	   getCubicSpline1stDifferential1D(xp, NDATA, xd, yd, y2),
	   getCubicSpline2ndDifferential1D(xp, NDATA, xd,     y2));
  }

  for(int ii = 0; ii < NDATA; ii++)
    fprintf(stderr, "%d\t%e\t%e\t%e\n", ii, xd[ii], yd[ii], y2[ii]);

  return (0);
}
