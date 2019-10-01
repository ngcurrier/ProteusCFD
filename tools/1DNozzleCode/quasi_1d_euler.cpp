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
#include <iomanip>
#include "blktrid.h"

#define R_UNIV 8.314459848      // J/mol.K

using namespace std;

template <class theType, class theType2>
inline theType2 MAX(theType x, theType2 y)
{
  return ((double(x) > double(y)) ? x : y);
};

template <class theType, class theType2>
inline theType2 MIN(theType x, theType2 y)
{
  return ((double(x) < double(y)) ? x : y);
};


int main(int argc, char* argv[]){

  //general variables
  int i, j, k, iter, ni, nj, nk;
  double RMS;
  double q0,q1,q2;

  //IO variables
  string file_resid = "conv.dat";
    
  //Initial conditions for the problem
  double cfltarget;
  double gamma;
  double ref_length, length;
  const double l_bound = 0.0;
  const int points = 51;
  const double tol = 1.0e-15;
  const double max_iter = 200000;
  //size of the blocks (i.e. 3x3 in 1d)
  nj = nk = 3;
  ni = points;

  double initialTemperature;     //K
  double initialPressure;        //Pa

  if(argc != 6){
    cerr << "USAGE: " << argv[0] << "<CFL> <nozzle Length (m)> <throat Diameter (m)> <upstream pressure (pa)> <upstream temperature (K)>" << endl;
    return(-2);
  }
  else{
    cfltarget = atof(argv[1]); 
    // we normalize the throat to be 1.0 non-dimensional lengths
    ref_length = atof(argv[3]);
    length = atof(argv[2])/ref_length;
    initialPressure = atof(argv[4]);
    initialTemperature = atof(argv[5]);
  }
  // PROBLEM SETUP REGION:
  // -------------------------------------------------------------
  // molecular weight of gas
  double MW = 28.966;                      //g/mol
  double Cp = 1006.43;                     //J/kg.K
  double initialMach = 0.5;                //Mach
  double machInit = 0.5;

  double backpressure = 325;            //Pa

  double ref_area = ref_length*ref_length; //m2
  double ref_pressure = initialPressure;   //Pa
  double backpressure_nd = backpressure / ref_pressure;
  std::cout << "PARAM: reference pressure - " << ref_pressure << " (Pa)" << std::endl;
  
  double Rs = R_UNIV/(MW/1000.0); //J/kg.K
  std::cout << "PARAM: specific gas constant - " << Rs << " (J/kg.K)" << std::endl;
  
  double Cv = Cp - Rs; //Cp = Cv + Rs - for ideal gas only
  std::cout << "PARAM: specific heat Cv - " << Cv << " (J/kg.K)" << std::endl;
  gamma = Cp/Cv;
  std::cout << "PARAM: gamma - " << gamma << " (J/kg.K)" << std::endl;

  //we use the speed of sound to dimensionalize velocity - compressible demands it
  double ref_velocity = sqrt(gamma*Rs*initialTemperature);
  std::cout << "PARAM: reference speed of sound - " << ref_velocity << " (m/s)" << std::endl;
    
  // warning: only one should ever be specified at a time (pressure or density)
  // we assume reference pressure takes precedence since it is easier to measure directly
  double ref_density = ref_pressure/(ref_velocity*ref_velocity);
  std::cout << "PARAM: reference density - " << ref_density << " (kg/m3)" << std::endl;

  double ref_time = ref_length/ref_velocity;
  std::cout << "PARAM: actual timestep reference - " << ref_time << " (s)" <<std::endl;
  double ref_enthalpy = ref_velocity*ref_velocity;
  std::cout << "PARAM: reference enthalpy - " << ref_enthalpy << std::endl;

  //we compute the reference temperature in this way b/c the compressible (Mach) based
  //solver requires this to be true
  double rho_d = initialPressure/(Rs*initialTemperature); //P = rho*Rs*T
  double P_nd = initialPressure/ref_pressure;
  double rho_nd = rho_d/ref_density;
  double c2_nd = gamma*P_nd/rho_nd;
  double ref_temperature = initialTemperature/c2_nd; //c2 = T when non-dimensionalized
  std::cout << "PARAM: reference temperature - " << ref_temperature << " (K)" << std::endl;

  //since we don't enforce it directly, check that the param object
  //contains a reasonable speed of sound c^2 = gamma*P/rho
  double T_nd = initialTemperature/ref_temperature; //from input file
  double p = initialPressure/ref_pressure; //from input file
  double rho = p*gamma/T_nd;
  double c2avg = gamma*(p/rho);
  double cavg = sqrt(c2avg);
  double cavg_dim = cavg * ref_velocity;

  const double pi = 3.14159265;

  //matrix variable pointers
  double *rhs = NULL;
  double *q = NULL;
  double *dq = NULL;
  double *tstep = NULL;
  double *upper = NULL;
  double *diag = NULL;
  double *lower = NULL; 
  double *A = NULL;
  double *Ah = NULL;
  double *fluxp = NULL; //flux from downwind (i.e. right)
  double *fluxn = NULL; // flux from upwind (i.e. left)
  double *x = NULL;
  double AC, AHC;
  
  
  //setup some IO functionality
  fstream resid_out(file_resid.c_str(), ios::out);

  if(!resid_out){
    cerr << "Files ~" << file_resid << " cannot be opened!" << endl;
    return(-1);
  }


  double *radius = new double[points];
  //cross sectional area
  double *s = new double[points];
  //cell volume
  double *v = new double[points+1];
  //coordinates
  x = new double[points];

  double step_size = (length - l_bound) / double(points-1);

  //compute x location for each grid point
  double dx = 0.0;
  double xstart = -1.0;
  for(i = 0; i < points; ++i){
    dx = i*step_size;
    x[i] = xstart + dx;
  }

  //compute radius for each grid point
  for(i = 0; i < points; i++){
  
    //equation for constantly expanding bell
    //radius[i] = 1.398 + 0.347 * tanh(0.8 * (x[i]) - 4.0);

    //equation for de laval - like nozzle
    //this is a bell-shaped nozzle
    //r = 1 − 0.868z2 + 0.432z3
    //radius[i] = 1.0 - 0.868*x[i]*x[i] + 0.432*x[i]*x[i]*x[i];

    // symmetric nozzle
    radius[i] = 1.0 + 1.5*pow((x[i]-(-0.5)),2.0);

    //compute section area
    s[i] = pi * radius[i] * radius[i];
  }

  //compute frustrum volume for each interior cell
  for(i = 1; i < points; i++){
    v[i] = ((step_size*pi)/3.0)*(radius[i-1]*radius[i-1] + 
			       radius[i-1]*radius[i]+ radius[i]*radius[i]); 
  }
  //initialize the ghost cell volumes
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
  tstep = new double[points];
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
    q[i*nj + 0] = rho;
    q[i*nj + 1] = machInit;
    q[i*nj + 2] = P_nd/(gamma-1.0) + 0.5*rho*machInit*machInit;
  }

  std::cout << "Q[0] Initial non-dimensional density: " << rho << std::endl;
  std::cout << "Q[1] Initial non-dimensional momentum: " << initialMach*rho << std::endl;
  std::cout << "Q[2] Initial non-dimensional total energy: " << P_nd/(gamma-1.0) * 0.5*rho*initialMach*initialMach << std::endl;



  iter = 0;
  do{

    // get static pressure in first inside cell
    q0 = q[1*nj + 0];
    q1 = q[1*nj + 1];
    q2 = q[1*nj + 2];
    p = (gamma -1.0) * (q2 - 0.5*q1*q1/q0);
    double gm1 = gamma-1.0;
    // use static pressure and total pressure to calculate velocity
    double u = sqrt(fabs((P_nd - p)*2.0/q0));
    if(P_nd < p){
      u = -u;
    }
    // calculate density using total temperature and velocity
    double T = T_nd/(1.0 + (gamma - 1.0)/2.0*u*u);
    rho = p*gamma/T;

    // updated bc's (inlet side)
    q[0*nj + 0] = rho;  //density
    q[0*nj + 1] = u*rho; //momentum
    q[0*nj + 2] = p/(gamma-1.0) + 0.5*rho*u*u;  //total energy
    


    //update bc's (outlet side) - this is only valid if nozzle is supersonic
    //since it is extrapolation from nozzle interior cell
    if (q[(points-1)*nj + 1]/q[(points-1)*nj + 0] > 1.0){ //supersonic outflow
      q[points*nj + 0] = q[(points-1)*nj + 0];
      q[points*nj + 1] = q[(points-1)*nj + 1];
      q[points*nj + 2] = q[(points-1)*nj + 2];
    }
    else{ //subsonic outflow - TODO: implement farfield CVBC
      rho = q[(points-1)*nj + 0];
      q[points*nj + 0] = q[(points-1)*nj + 0];
      q[points*nj + 1] = q[(points-1)*nj + 1];
      double u = q[(points-1)*nj + 1]/rho;
      q[(points)*nj + 2] = backpressure_nd/(gamma-1.0) + 0.5*rho*u*u;
    }
 
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

      //compute the eigenvalues (u, u+c, u-c)
      //p = (gam-1)(r Et - 1/2 (r u^2))
      p = (gamma -1.0) * (q2 - 0.5*q1*q1/q0);
      double c = sqrt(gamma*p/rho);
      double u = q1;
      double upc = u + c;
      double umc = u - c;
      double maxeig = MAX(fabs(upc), fabs(umc));
      double vol = v[i];
      tstep[i] = cfltarget*step_size/maxeig;

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
	diag[i*nj*nk + j*nk + j] += 1.0/tstep[i];
#else
	diag[i*nj*nk + j*nk + j] += 2.0/tstep[i];
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
	RMS += rhs[i*nj + j] * rhs[i*nj + j];
      }
    }
    RMS = sqrt(RMS) / (points-1);
    resid_out << RMS << endl;
    if(iter % 1000 == 0){
      cout << iter << ":: RMS: " << RMS << endl;
    }
    iter++;
    
  }while(RMS > tol && iter < max_iter  && !isnan(RMS) && !isinf(RMS));

  if(isnan(RMS) || isinf(RMS)){
    cerr << endl;
    cerr << "Solution has blown up!!" << endl;
    cerr << endl;
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
  cout << "Position_m\t Radius_m\t Area_m2\t Density_kgm3\t Mach\t\t  Velocity_ms\t Pressure_Pa\t  Temperature_K\t MassFlow_kgs" << endl;
  cout << "========================================================================================================================================" << endl;
  double outletMassFlow = 0.0;
  double outletVelocity = 0.0;
  double outletPressure = 0.0;
  for(i = 0; i < points; i++){
    //get values averaged at grid points
    q0 = 0.5*(q[i*nj + 0] + q[(i+1)*nj + 0]);
    q1 = 0.5*(q[i*nj + 1] + q[(i+1)*nj + 1]);
    q2 = 0.5*(q[i*nj + 2] + q[(i+1)*nj + 2]);
    double p_local = (gamma -1.0) * (q2 - 0.5*q1*q1/q0);
    double rho_local = q0;
    double T_local = gamma*p_local/rho_local;
    double u_ms = (q1/q0)*ref_velocity;
    double rho_kgm3 = rho_local*ref_density;
    double area_m2 = s[i]*ref_area;
    cout << std::fixed << std::setprecision(6) << x[i]*ref_length << "  \t" << std::setprecision(6) 
	 << radius[i]*ref_length << "   \t" <<std::setprecision(6) << s[i]*ref_area << "  \t" << std::setprecision(6) << rho_kgm3 << " \t" 
	 << std::setprecision(6) << (q1/q0) << " \t" << std::setprecision(6) << u_ms << " \t" << std::setprecision(6) << p_local*ref_pressure << " \t" 
	 << std::setprecision(6) << T_local*ref_temperature << " \t" << std::setprecision(6) << u_ms*rho_kgm3*area_m2 << endl;
    outletMassFlow = u_ms*rho_kgm3*area_m2;
    outletVelocity = u_ms;
    outletPressure = p_local*ref_pressure;
  }
  double minArea = s[0];
  for(i = 0; i < points; ++i){
    minArea = MIN(minArea, s[i]);
  }
  std::cout << "Throat area: " << minArea << std::endl;
  std::cout << "Nozzle area ratio: " << s[points-1]/minArea << std::endl;
  double pressureThrust = (outletPressure - backpressure)*s[points-1]*ref_area;
  double momentumThrust = outletMassFlow * outletVelocity;
  std::cout << "Thrust (pressure - kN): " << pressureThrust/1000.0 << std::endl;
  std::cout << "Thrust (momentum - kN): " << momentumThrust/1000.0 << std::endl;
  std::cout << "Thrust (total - kN): " << (pressureThrust + momentumThrust)/1000.0 << std::endl;

  if (backpressure > outletPressure){
    std::cout << "WARNING: Nozzle is overexpanded" << std::endl;
  }
  else{
    std::cout << "Nozzle is underexpanded" << std::endl;
  }

  std::string filename = "solution.csv";
  fstream csv_out(filename.c_str(), ios::out);
  csv_out << "Position_m,Radius_m,Density_kgm3,Mach,Velocity_ms,Pressure_pa,Temperature_K,MassFlow_kgs" << std::endl;
  for(i = 0; i < points; i++){
    //get values averaged at grid points
    q0 = 0.5*(q[i*nj + 0] + q[(i+1)*nj + 0]);
    q1 = 0.5*(q[i*nj + 1] + q[(i+1)*nj + 1]);
    q2 = 0.5*(q[i*nj + 2] + q[(i+1)*nj + 2]);
    double p_local = (gamma -1.0) * (q2 - 0.5*q1*q1/q0);
    double rho_local = q0;
    double T_local = gamma*p_local/rho_local;
    double u_ms = (q1/q0)*ref_velocity;
    double rho_kgm3 = rho_local*ref_density;
    double area_m2 = s[i]*ref_area;
    csv_out << x[i]*ref_length << "," << radius[i]*ref_length << "," << rho_kgm3 << "," 
	    << (q1/q0) << "," << u_ms << "," << p_local*ref_pressure << "," << T_local*ref_temperature << "," 
	    << u_ms*rho_kgm3*area_m2 << endl;
  }
  csv_out.close();
  

  resid_out.close();

  
  delete [] radius;
  delete [] s;
  delete [] v;
  delete [] rhs;
  delete [] q;
  delete [] dq;
  delete [] tstep;
  delete [] upper;
  delete [] lower;
  delete [] diag;
  delete [] A;
  delete [] Ah;
  delete [] fluxp;
  delete [] fluxn;
  delete [] x;
}
