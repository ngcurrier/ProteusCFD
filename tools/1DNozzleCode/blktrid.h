#ifndef BLOCK_TRID_H
#define BLOCK_TRID_H
/*
Function solves a block tri-diagonal system: Tri(a,b,c) * x = r
Number of blocks = ni
Block size = nk*nk

Input: arrays a,b,c,r    (see arguments for dimensioning)
       nj,nk             (problem size)

Output: x                (the solution, see argument for dimensioning)

Note: arrays a,b,c,r are all changed in this routine and should be rebuilt
      before calling again

Note2: arrays expected to be passed in so a[i][j][k] == a[i*nj*nk + j*nk + k]
       or in 2D a[i][j] == a[i*ni + j]
*/

#include <stdlib.h>
#include <stdio.h>

void blktrid( double *a, double *b,double *c, double *x, double *r,
	      int ni, int nk);
#endif
