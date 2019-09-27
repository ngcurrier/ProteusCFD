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
       or in 2D a[i][j] == a[i*nj + j]
*/

#include "blktrid.h"

//=======================================================================
void blktrid(double *a, double *b, double *c, double *x, double *r,
	     int ni, int nk)
{
  //-----------------------------------------------------------------------
  
  int n,m,istart,i,j,k,kk,mp,mr,jb,kr ;
  double t,s ;

  //-----------------------------------------------------------------------
  // initialize index limits for c arrays
  n = ni-1;
  m = nk-1;
  int nj = nk;
  istart = 1;
  //-----------------------------------------------------------------------
  // forward elimination-blocks
  for (i = istart ; i <= n ; i++ ) {
    
    // gauss jordan for b(i,*,*)
    for (j = 0 ; j <= m-1 ; j++ ) {
      
      
      // normalize b,c,r using pivot row in b
      t = 1.0 / b[i*nj*nk + j*nk + j];
      for ( k = j+1 ; k <= m ; k++ ) {
	b[i*nj*nk + j*nk + k] = t * b[i*nj*nk + j*nk + k];
      }
      if ( i != n ) {
	for (k = 0 ; k <= m ; k++ ) {
	  c[i*nj*nk + j*nk + k] = t * c[i*nj*nk + j*nk + k];
	}
      }
      r[i*nj + j] = t * r[i*nj + j];
      
      // lower elimination in b,c,r
      for ( k = j+1 ; k <= m ; k++ ) {
	for ( kk = j+1 ; kk <= m ; kk++ ) {
	  b[i*nj*nk + k*nk + kk] = b[i*nj*nk + k*nk + kk] - b[i*nj*nk + k*nk + j] * 
	    b[i*nj*nk + j*nk + kk];
	}
	if ( i != n ) {
	  for (kk = 0 ; kk <= m ; kk++ ) {
	    c[i*nj*nk + k*nk + kk] = c[i*nj*nk + k*nk + kk] - b[i*nj*nk + k*nk + j] * 
	      c[i*nj*nk + j*nk + kk];

	  }
	}
	r[i*nj + k] = r[i*nj + k] - b[i*nj*nk + k*nk + j] * r[i*nj + j];
      } 
      // end of k loop
      
    } 
    // end of j loop
    
    // upper elimination of c,r (gauss jordan)
    t = 1.0 / b[i*nj*nk + m*nk + m];
    for ( kk = 0 ; kk <= m ; kk++ ) {
      c[i*nj*nk + m*nk + kk] = t * c[i*nj*nk + m*nk + kk];
    }
    r[i*nj + m] = t * r[i*nj + m];
    for ( j = 0 ; j <= (m-1) ; j++ ) {
      mp = m-j ;// removed a + 1 
      for ( k = 0 ; k <= mp-1 ; k++ ) {
	mr = mp-k-1 ;// added a - 1
	if( i != n ) {
	  for ( kk = 0 ; kk <= m ; kk++ ) {
	    c[i*nj*nk + mr*nk + kk] = c[i*nj*nk + mr*nk + kk] - b[i*nj*nk + mr*nk + mp] * 
	      c[i*nj*nk + mp*nk + kk];
	  }
	}
	r[i*nj + mr] = r[i*nj + mr] - b[i*nj*nk + mr*nk + mp] * r[i*nj + mp];
      }
    }
    // end of j loop
    
    // eliminate a(i+1,*,*),  b(i,*,*) is now the unit matrix
    if ( i != n ) {
      for ( j = 0 ; j <= m ; j++ ) {
	for ( k = 0 ; k <= m ; k++ ) {
	  for ( kk = 0 ; kk <= m ; kk++ ) {
	    b[(i+1)*nj*nk + k*nk + kk] = b[(i+1)*nj*nk + k*nk + kk] - 
	      a[(i+1)*nj*nk + k*nk + j] * c[i*nj*nk + j*nk + kk];
	  }
	  r[(i+1)*nj + k] = r[(i+1)*nj + k] - a[(i+1)*nj*nk + k*nk + j] * r[i*nj + j];
	}
      }
    }
    
  } 
  // end of i loop
  //-----------------------------------------------------------------------
  // back substitution
  for (k = 0 ; k <= m ; k++ ) {
    x[n*nj + k] = r[n*nj + k];
  }
  
  for (j = 0 ; j <= n-istart-1 ; j++ ) {					
    jb = n-j-1 ;								
    for (k = 0 ; k <= m ; k++ ) {
      kr = m-k ;								
      s = 0.0 ;
      for (kk = 0 ; kk <= m ; kk++ ) {
	s = c[jb*nj*nk + kr*nk + kk] * x[(jb+1)*nj + kk] + s;
      }
      x[jb*nj + kr] = r[jb*nj + kr] - s;
    }
  }

}

