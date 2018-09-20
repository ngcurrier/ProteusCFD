#include "modal.h"
#include <lapacke.h>

//compute eigenvalues and vectors of problem m*x*lambda = -k*x
//  From Lapacke documentation of dggevx (conditioned general eigenvalue problem)
//  The 'x' after dggev stands for expert mode which allows for preconditioning of the problem
//  set to 'N' if you want similar answers to Matlab
//  The right eigenvector v(j) corresponding to the eigenvalue lambda(j) of (A,B) satisfies
//  A * v(j) = lambda(j) * B * v(j) .
//  The left eigenvector u(j) corresponding to the eigenvalue lambda(j) of (A,B) satisfies
//  u(j)**H * A  = lambda(j) * u(j)**H * B.
//  where u(j)**H is the conjugate-transpose of u(j).

//Documentation for LAPACK here:
//http://www.netlib.org/lapack/lapack-3.1.1/html/dggevx.f.html

//SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB,
//                   ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO,
//                   IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE,
//                   RCONDV, WORK, LWORK, IWORK, BWORK, INFO )
//-----------------------------------------------------------------
// CHARACTER         BALANC, JOBVL, JOBVR, SENSE
// INTEGER           IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N
// DOUBLE PRECISION  ABNRM, BBNRM
// All other arguments are arrays

info=LAPACKE_dggevx(LAPACK_COL_MAJOR, 'B', 'V', 'V', 'N', n, a, n, b, n,
		    alphar, alphai, beta, vl, n, vr, n, &ilo,
		    &ihi, lscale, rscale, &abnrm, &bbnrm, rconde,
		    rcondv);


//Tutorial on solving this problem with Eigen library
//http://eigen.tuxfamily.org/index.php?title=Lapack
