//**********************************************************
// Quasi 1D Euler Equation Solver
// (sonic flow in a diverging nozzle)
// 
// Nick Currier
// October 2008
// CFD I
// University of Tennessee at Chattanooga
//
// Please, no redistribution without consent
//*********************************************************


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "blktrid.h"

using namespace std;

int main(int argc, char* argv[]){

  //general variables
  int i, j, k, iter, ni, nj, nk;
  double RMS;
  double q0,q1,q2;

  //IO variables
  string file_resid = "conv.dat";
  string file_sol = "sol.dat";
    
  //Initial conditions for the problem
  double t_step;
  const double l_bound = 0.0;
  const double extent = 10.0;
  const int points = 51;
  const double mach = 1.5;
  const double gamma = 1.4;
  const double tol = 1.0e-15;
  const double max_iter = 300000;
  //size of the blocks (i.e. 3x3 in 1d)
  nj = nk = 3;
  ni = points;

  const double pi = 3.14159265;


  //matrix variable pointers
  double *rhs = NULL;
  double *q = NULL;
  double *dq = NULL;
  double *upper = NULL;
  double *diag = NULL;
  double *lower = NULL; 
  double *A = NULL;
  double *Ah = NULL;
  double p = 0.0; //pressure
  double *fluxp = NULL; //flux from downwind (i.e. right)
  double *fluxn = NULL; // flux from upwind (i.e. left)
  double AC, AHC;
  
  if(argc != 2){
    cerr << "USAGE: ./a.out 'timestep desired' " << endl;
    return(-2);
  }
  else{
    t_step = atof(argv[1]); 
  }
  
  //setup some IO functionality
  fstream resid_out(file_resid.c_str(), ios::out);
  fstream sol_out(file_sol.c_str(), ios::out);

  if(!resid_out || !sol_out){
    cerr << "Files ~" << file_resid << " and ~" << file_sol 
	 << " cannot be opened!" << endl;
    return(-1);
  }


  double *radius = new double[points];
  //cross sectional area
  double *s = new double[points];
  //cell volume
  double *v = new double[points+1];

  double step_size = (extent - l_bound) / (points-1);

  //compute radius for each grid point
  for(i = 0; i < points; i++){
    radius[i] = 1.398 + 0.347 * tanh(0.8 * (i * step_size) - 4.0);
    s[i] = pi * radius[i] * radius[i];
  }

  //compute volume for each interior cell
  for(i = 1; i < points; i++){
    v[i] = (step_size*pi)/3.0*(radius[i-1]*radius[i-1] + 
			       radius[i-1]*radius[i]+ radius[i]*radius[i]); 
  }
  //initialize the ghost cells
  v[0] = v[1]; 
  v[points] = v[points-1];


  // setup memory for solving the equations
  // use block memory allocation for speed
  // access like a[k*nj*nk + j*nk + i] == a[k][j][i]
  // but much quicker with less indirection
  //
  // ** view as i number of j x k matrices with j 
  //    x direction and k in the y direction
  //
  //              |j->     |
  // a(i,j,k) ==  |        |
  //              |        |i
  //
  // ** note (2D) : a[i*nj + j] == a[i][j]
  

  //allocate memory for the system
  rhs = new double[(points + 1) * nj];
  q = new double[(points + 1) * nj];
  dq = new double[(points + 1) * nj];
  upper = new double[(points + 1) * nj * nj];
  diag = new double[(points + 1) * nj * nj];
  lower = new double[(points + 1) * nj * nj]; 
  A = new double[(points) * nj * nj];
  Ah = new double[(points + 1) * nj * nj];
  fluxp = new double[nj];
  fluxn = new double[nj];
  

  //initialize the diagonal blocks (i.e. set to unity)
  for(i = 0; i < points + 1; i++){
    for(j = 0; j < 3; j++){
      for(k = 0; k < 3; k++){
	//off diagonals
	diag[i*nj*nk + j*nk + k] = 0.0;
      }	
      //initialize dq
      dq[i*nj + j] = 0.0;
    }
    //initialize q
    q[i*nj + 0] = 1.0;
    q[i*nj + 1] = mach;
    q[i*nj + 2] = 1.0/(gamma*(gamma-1.0)) + 0.5*mach*mach;
  }


  iter = 0;
  do{

    //initialize bc's
    q[0*nj + 0] = 1.0; 
    q[0*nj + 1] = mach;
    q[0*nj + 2] = 1.0/(gamma*(gamma-1.0)) + 0.5*mach*mach;
    q[points*nj + 0] = q[(points-1)*nj + 0];
    q[points*nj + 1] = q[(points-1)*nj + 1];
    q[points*nj + 2] = q[(points-1)*nj + 2];;
      
 
    //compute flux jacobians A = df/dq
    for(i = 0; i < points; i++){
      //AC = s[i]/(0.5*(v[i]+v[i+1]));
      AC = s[i]/v[i];

#ifdef _SONIC
      //use upwind values
      q0 = q[i*nj + 0];
      q1 = q[i*nj + 1];
      q2 = q[i*nj + 2];
#else
      //average the q vector at the faces
      q0 = 0.5*(q[i*nj + 0] + q[(i+1)*nj + 0]);
      q1 = 0.5*(q[i*nj + 1] + q[(i+1)*nj + 1]);
      q2 = 0.5*(q[i*nj + 2] + q[(i+1)*nj + 2]);
#endif

      A[i*nj*nk + 0*nk + 0] = 0.0;
      A[i*nj*nk + 0*nk + 1] = AC * 1.0;
      A[i*nj*nk + 0*nk + 2] = 0.0;
      A[i*nj*nk + 1*nk + 0] = AC * (gamma-3.0)/2.0 * ((q1*q1)/(q0*q0));
      A[i*nj*nk + 1*nk + 1] = AC * (3.0-gamma)*(q1/q0);
      A[i*nj*nk + 1*nk + 2] = AC * (gamma-1.0);
      A[i*nj*nk + 2*nk + 0] = AC * (q1/(q0*q0))*(-gamma*q2 + (gamma-1.0)*(q1*q1/q0));
      A[i*nj*nk + 2*nk + 1] = AC * (1.0/q0)*(gamma*q2 -1.5*(gamma-1.0)*(q1*q1/q0));
      A[i*nj*nk + 2*nk + 2] = AC * gamma * q1/q0;
    
      //compute more jacobians for each cell Ah = dh/dq
      if(i > 0 && i < points-1){
	//use q vector from the cell
	q0 = q[i*nj + 0];
	q1 = q[i*nj + 1];
	q2 = q[i*nj + 2];
	
	AHC = (gamma-1.0)*(s[i]-s[i-1])/v[i];
	
	Ah[i*nj + 0*nk + 0] = 0.0;
	Ah[i*nj + 0*nk + 1] = 0.0;
	Ah[i*nj + 0*nk + 2] = 0.0;
	Ah[i*nj + 1*nk + 0] = AHC*0.5*((q1*q1)/(q0*q0));
	Ah[i*nj + 1*nk + 1] = -AHC*(q1/q0);
	Ah[i*nj + 1*nk + 2] = AHC*1.0;
	Ah[i*nj + 2*nk + 0] = 0.0;
	Ah[i*nj + 2*nk + 1] = 0.0;
	Ah[i*nj + 2*nk + 2] = 0.0;
      }
    }
    
    //update rhs (for 1D only)
    for(i = 1; i < points; i++){

#ifdef _UPWINDFLUX
      //use upwind values
      q0 = q[i*nj + 0];
      q1 = q[i*nj + 1];
      q2 = q[i*nj + 2];
#else
      //use averages at the cell faces
      q0 = 0.5*(q[i*nj + 0]+q[(i+1)*nj + 0]);
      q1 = 0.5*(q[i*nj + 1]+q[(i+1)*nj + 1]);
      q2 = 0.5*(q[i*nj + 2]+q[(i+1)*nj + 2]);
#endif

      //p = (gam-1)(r Et - 1/2 (r u^2))
      p = (gamma -1.0) * (q2 - 0.5*q1*q1/q0);
      //fluxp[0] = (r u)_{i+1/2}
      fluxp[0] = s[i] * q1;
      //fluxp[1] = (r u^2 + p)_{i+1/2}
      fluxp[1] = s[i] * ((q1*q1)/q0 + p);
      //fluxp[2] = u(r Et + P)_{i+1/2}
      fluxp[2] = s[i] * ((q1/q0)*(q2 + p));

#ifdef _UPWINDFLUX
      //use upwind values
      q0 = q[(i-1)*nj + 0];
      q1 = q[(i-1)*nj + 1];
      q2 = q[(i-1)*nj + 2]; 
#else
      //use averages at the cell faces
      q0 = 0.5*(q[i*nj + 0]+q[(i-1)*nj + 0]);
      q1 = 0.5*(q[i*nj + 1]+q[(i-1)*nj + 1]);
      q2 = 0.5*(q[i*nj + 2]+q[(i-1)*nj + 2]);
#endif
      
      p = (gamma -1.0) * (q2 - 0.5*q1*q1/q0);

      //fluxn[0] = (r u)_{i-1/2}
      fluxn[0] = s[i-1] * q1;
      //fluxn[1] = (r u^2 + p)_{i-1/2}
      fluxn[1] = s[i-1] * ((q1*q1)/q0 + p);
      //fluxn[2] = u(r Et + P)_{i-1/2}
      fluxn[2] = s[i-1] * ((q1/q0)*(q2 + p));
      
      //calculate source terms at the cell center
      q0 = q[i*nj + 0];
      q1 = q[i*nj + 1];
      q2 = q[i*nj + 2];
      p = (gamma -1.0) * (q2 - 0.5*q1*q1/q0);

      double h0 = 0.0;
      double h1 = p * (s[i] - s[i-1]);
      double h2 = 0.0;
#ifdef _SONIC
      rhs[i*nj + 0] = -1.0/v[i]*(fluxp[0] - fluxn[0] - h0);
      rhs[i*nj + 1] = -1.0/v[i]*(fluxp[1] - fluxn[1] - h1);
      rhs[i*nj + 2] = -1.0/v[i]*(fluxp[2] - fluxn[2] - h2);
#else
      rhs[i*nj + 0] = -2.0/v[i]*(fluxp[0] - fluxn[0] - h0);
      rhs[i*nj + 1] = -2.0/v[i]*(fluxp[1] - fluxn[1] - h1);
      rhs[i*nj + 2] = -2.0/v[i]*(fluxp[2] - fluxn[2] - h2);
#endif
      //      cout << " i "<< i << " fn[0] " << fluxn[0] <<" fn[1] "<< fluxn[1] 
      //   <<" fn[2] "<<  fluxn[2] << " p " << p << endl;
    }

    //build upper and lower operators
    for(i = 1; i < points; i++){
      for(j = 0; j < nj; j++){
	for(k = 0; k < nk; k++){

#ifdef _SONIC
	  //using upwind algorithm
	  lower[i*nj*nk + j*nk + k] = -A[(i-1)*nj*nk + j*nk + k];
	  //upper[i*nj*nk + j*nk + k] = A[i*nj*nk + j*nk + k];
	  upper[i*nj*nk + j*nk + k] = 0.0;
	  diag[i*nj*nk + j*nk + k] = - Ah[i*nj*nk + j*nk + k]  
	    + A[i*nj*nk + j*nk + k];
#else
	  //using standard algorithm
	  lower[i*nj*nk + j*nk + k] = -A[(i-1)*nj*nk + j*nk + k];
	  upper[i*nj*nk + j*nk + k] = A[i*nj*nk + j*nk + k];
	  diag[i*nj*nk + j*nk + k] = -2.0 * Ah[i*nj*nk + j*nk + k] + 
	    A[i*nj*nk + j*nk + k] - A[(i-1)*nj*nk + j*nk + k];
#endif	  

	}
	//add extra term on the diagonal of the diagonal blocks
#ifdef _SONIC
	diag[i*nj*nk + j*nk + j] += 1.0/t_step;
#else
	diag[i*nj*nk + j*nk + j] += 2.0/t_step;
#endif
      }
    }

    //solve for dq
    blktrid(lower, diag, upper,	dq, rhs, ni, nk);
    
    //update q vector (i.e. solution)
    //    q[i][0] = rho
    //    q[i][1] = rho * u
    //    q[i][2] = rho * E_total
    for(i = 1; i < points+1; i++){
      for(j = 0; j < nj; j++){
	q[i*nj + j] += dq[i*nj + j];
      }
    }


    //compute RMS error as a convergence criterion
    RMS = 0.0;
    for(i = 1; i < points; i++){
      for(j = 0; j < 3; j++){
	RMS += dq[i*nj + j] * dq[i*nj + j];
      }
    }
    RMS = sqrt(RMS) / (points-1);
    resid_out << RMS << endl;
    if(iter % 75 == 0){
      cout << iter << ":: RMS: " << RMS << endl;
    }
    iter++;
    
  }while(RMS > tol && iter < max_iter  && !isnan(RMS) && !isinf(RMS));

  if(isnan(RMS) || isinf(RMS)){
    cerr << endl;
    cerr << "Solution has blown up!!" << endl;
    cerr << endl;
    return (-71);
  }
  
  cout << endl;
  cout << iter << " iterations taken" << endl;
  cout << endl;

  cout.precision(13);

#ifdef _VERBOSE
  for(i = 0; i < points + 1; i++){
    for(j = 0; j < 3; j++){
      cout << q[i*nj + j] << " " ;
    }
    cout << endl;
  }
#endif

  cout << endl;
  cout << "Position   Density   Velocity   Pressure   MassFlow" << endl;
  cout << "===================================================" << endl;
  for(i = 0; i < points; i++){
    //get values averaged at grid points
    q0 = 0.5*(q[i*nj + 0] + q[(i+1)*nj + 0]);
    q1 = 0.5*(q[i*nj + 1] + q[(i+1)*nj + 1]);
    q2 = 0.5*(q[i*nj + 2] + q[(i+1)*nj + 2]);
    p = (gamma -1.0) * (q2 - 0.5*q1*q1/q0);
    cout << i*step_size << " " << q0 << " " 
	 << (q1/q0) << " " << p << " " << q1*s[i] << endl;
    sol_out << i*step_size << " " << q0 << " " 
	    << (q1/q0) << " " << p << endl;
  }
  
  resid_out.close();
  sol_out.close();

  
  delete [] radius;
  delete [] s;
  delete [] v;
  delete [] rhs;
  delete [] q;
  delete [] dq;
  delete [] upper;
  delete [] lower;
  delete [] diag;
  delete [] A;
  delete [] Ah;
  delete [] fluxp;
  delete [] fluxn;
}
