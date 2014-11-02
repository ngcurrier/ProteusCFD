#include "utility.h"

namespace STRUCTDYN{

int AnalyticSolution(int dof, int nsteps, double dt, double* m, double* c, 
		     double* k, double* x, double* xd, double* xdd, SForces* forces)
{
  int i, j;
  int err = 0;
  
  //compute natural frequency (this only works in 1d)
  double wn = sqrt(k[0]/m[0]);
  //compute critical damping
  double cc = 2.0*m[0]*wn;
  //compute damping ratio
  double zeta = c[0] / cc;
  //compute damped frequency
  double wd = wn*sqrt(1.0-zeta*zeta);
  //define cyclic period
  double tau;

  //time
  double t;  
  //solution coefficients
  double A, B;
  //temporary coefficients
  double c1, c2, c3;

  //extract from forces object for readability
  double fo = forces->mag;
  double w = forces->omega;

  //if free vibration
  if(fo == 0.0){
    //underdamped
    if(zeta < 1.0){
      //solve for solution coefficients
      A = (xd[0] + zeta*wn*x[0])/wd;
      B = x[0];
      c1 = sqrt(zeta*zeta - 1.0);
      for(i = 0; i <= nsteps; i++){
	t = dt*double(i);
	c2 = exp(-zeta*wn*t);
	for(j = 0; j < dof; j++){
	  x[i*dof + j] = c2*(A*sin(wd*t) + B*cos(wd*t));
	  xd[i*dof + j] = -zeta*wn*c2*(A*sin(wd*t) + B*cos(wd*t)) + 
	    c2*(A*wd*cos(wd*t) - B*wd*sin(wd*t));
	}
      }
    }
    //overdamped
    else if(zeta > 1.0){
      //solve for solution coefficients
      c1 = sqrt(zeta*zeta - 1.0);
      c2 = (xd[0] + zeta*wn*x[0])/(c1*wn*2.0);
      A = x[0]/2.0 + c2;
      B = x[0]/2.0 - c2;
      for(i = 0; i <= nsteps; i++){
	t = dt*double(i);
	c3 = exp(-zeta*wn*t);
	for(j = 0; j < dof; j++){
	  x[i*dof + j] = c3*(A*exp(c1*wn*t) + B*exp(-c1*wn*t));
	  xd[i*dof + j] = -zeta*wn*c3*(A*exp(c1*wn*t) + B*exp(-c1*wn*t)) +
	    c3*(A*c1*wn*exp(c1*wn*t) - B*c1*wn*exp(-c1*wn*t));
	}
      }
    }
    //critically damped
    else if(zeta == 1.0){
      //solve for solution coefficients
      A = x[0];
      B = xd[0] + wn*x[0];
      for(i = 0; i <= nsteps; i++){
	t = dt*double(i);
	for(j = 0; j < dof; j++){
	  x[i*dof + j] = exp(-wn*t)*(A + B*t);
	  xd[i*dof + j] = -wn*exp(-wn*t)*(A + B*t) + exp(-wn*t)*B;
	}
      }
    }
    //no damping (we'll never get here but written for completeness)
    else{
      //solve for solution coefficients based on ICs
      A = xd[0] / wn;
      B = x[0];
      for(i = 0; i <= nsteps; i++){
	t = dt*double(i);
	for(j = 0; j < dof; j++){
	  x[i*dof + j] = A*sin(wn*t) + B*cos(wn*t);
	  xd[i*dof + j] = A*wn*cos(wn*t) - B*wn*sin(wn*t);
	}
      }
    }
  }
  //forced vibration
  else{
    double dst = fo/k[0];
    double r = w/wn;
    double cwk = 2.0*zeta*r;
    double Xp = dst/(sqrt(cwk*cwk + (1.0 - r*r)*(1.0 - r*r)));
    double phip = atan2(cwk,(1.0 - r*r));

    std::cout << "Forced vibration " << std::endl;
    std::cout << "================" << std::endl;
    std::cout << "Frequency ratio (r): " << r << std::endl;
    std::cout << "Particular solution magnitude: " << Xp << std::endl;
    std::cout << "Particular solution phase shift: " << phip << std::endl;
    std::cout << "Statically forced displacement: " << dst << std::endl;

    //underdamped
    if(zeta < 1.0){
      c1 = sqrt(zeta*zeta - 1.0);
      //solve for solution coefficients
      B = x[0] - Xp*sin(-phip);
      A = (xd[0] + zeta*wn*B - w*Xp*cos(-phip))/wd;
      for(i = 0; i <= nsteps; i++){
	t = dt*double(i);
	c2 = exp(-zeta*wn*t);
	for(j = 0; j < dof; j++){
	  x[i*dof + j] = c2*(A*sin(wd*t) + B*cos(wd*t)) + Xp*sin(w*t - phip);
	  xd[i*dof + j] = -zeta*wn*c2*(A*sin(wd*t) + B*cos(wd*t)) + 
	    c2*(A*wd*cos(wd*t) - B*wd*sin(wd*t)) + w*Xp*cos(w*t - phip);
	}
      }
    }
    //overdamped
    else if(zeta > 1.0){
      //solve for solution coefficients
      c1 = sqrt(zeta*zeta - 1.0);
      B = (-xd[0] - zeta*wn*(x[0]-Xp*sin(-phip)) )/(c1*wn*2.0) + 
	(r*Xp*cos(-phip))/(c1*2.0) + x[0]/2.0 -Xp*sin(-phip)/2.0;
      A = x[0] - B - Xp*sin(-phip);
      for(i = 0; i <= nsteps; i++){
	t = dt*double(i);
	c3 = exp(-zeta*wn*t);
	for(j = 0; j < dof; j++){
	  x[i*dof + j] = c3*(A*exp(c1*wn*t) + B*exp(-c1*wn*t)) + Xp*sin(w*t - phip);
	  xd[i*dof + j] = -zeta*wn*c3*(A*exp(c1*wn*t) + B*exp(-c1*wn*t)) +
	    c3*(A*c1*wn*exp(c1*wn*t) - B*c1*wn*exp(-c1*wn*t)) + w*Xp*cos(w*t - phip);
	}
      }
    }
    //critically damped
    else if(zeta == 1.0){
      //solve for solution coefficients
      A = x[0] - Xp*sin(-phip);
      B = wn*A + xd[0] - w*Xp*cos(-phip);
      for(i = 0; i <= nsteps; i++){
	t = dt*double(i);
	for(j = 0; j < dof; j++){
	  x[i*dof + j] = exp(-wn*t)*(A + B*t) + Xp*sin(w*t - phip);
	  xd[i*dof + j] = -wn*exp(-wn*t)*(A + B*t) + exp(-wn*t)*B + w*Xp*cos(w*t - phip);
	}
      }
    }
    //no damping (we'll never get here)
    else{
    }
  }
  return err;
}

void MSolve(double* m, double* vin, double* vout, int ndof)
{
  int i;
  for(i = 0; i < ndof; i++){
    vout[i] = vin[i]/m[i*ndof + i];
  }
  return;
}


void WriteSolution(int dof, int nsteps, double dt, double* x, double* xd, double* xdd)
{
  int i, j;
  std::ofstream fout;
  double time; 
  
  fout.open("solution.dat");

  int prec = 16;
  fout.setf(std::ios::scientific);
  //std::cout.setf(std::ios::fixed);
  fout.precision(prec);


  fout << "#time  position(1..dof) velocity(1...dof) acceleration(1..dof) .. etc." 
       << std::endl;
  for(i = 0; i <= nsteps; i++){
    time = dt*(double)i;
    fout << time << " ";
    for(j = 0; j < dof; j++){
      fout << x[i*dof + j] << " ";
    }
    for(j = 0; j < dof; j++){
      fout << xd[i*dof + j] << " ";
    }
    for(j = 0; j < dof; j++){
      fout << xdd[i*dof + j] << " ";
    }
    fout << std::endl;
  }

  fout.close();
  return;
}


void Read1DCase(int* dof_, int nsteps, double dt, double** x_, double** xd_, 
		double** xdd_,double** m_, double** c_, double** k_, 
		SForces* forces)
{
  int i, j;
  std::string filename = "config.dat";
  std::string junk;
  std::ifstream fin;
  
  //open and read config file
  fin.open(filename.c_str());
  //throw away first line, it's instructions
  getline(fin, junk);
  fin >> (*dof_);

  //create easily accessible variables
  int dof = (*dof_);
  forces->Init(dof, 1, dof);

  int tfinal = (double)nsteps*dt;
  int nptp;

  //allocate memory for the next section
  (*m_) = new double[dof*dof];
  (*c_) = new double[dof*dof];
  (*k_) = new double[dof*dof];
  
  //allocate solution memory
  (*x_) = new double[dof*(nsteps+1)];
  (*xd_) = new double[dof*(nsteps+1)];
  (*xdd_) = new double[dof*(nsteps+1)];

  //create easily accessible variables
  double* m = (*m_);
  double* c = (*c_);
  double* k = (*k_);
  double* x = (*x_);
  double* xd = (*xd_);
  double* xdd = (*xdd_);

  //read system parameters
  for(i = 0; i < dof*dof; i++){
    fin >> m[i];
  }
  for(i = 0; i < dof*dof; i++){
    fin >> c[i];
  }
  for(i = 0; i < dof*dof; i++){
    fin >> k[i];
  }

  //read force magnitude and frequency
  double fo, omega;
  fin >> fo;
  fin >> omega;
  
  //set these in the forces object
  forces->Set1DCase(omega, fo);

  //read IC's
  for(i = 0; i < dof; i++){
    fin >> x[i];
  }
  for(i = 0; i < dof; i++){
    fin >> xd[i];
  }
  fin.close();
  
  //echo config file to screen for error checking
  std::cout << "\nThe following was read from command line" << std::endl;
  std::cout << "======================================\n" << std::endl; 
  std::cout << "Timestep: " << dt << std::endl;
  std::cout << "Number of steps: " << nsteps << std::endl;
  std::cout << "Time final: " << tfinal << std::endl;
  
  std::cout << "\nThe following was read from config.dat" << std::endl;
  std::cout << "======================================\n" << std::endl;
  std::cout << "Mass:" << std::endl;
  for(i = 0; i < dof; i++){
    for(j = 0; j < dof; j++){
      std::cout << m[i*dof + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Damping:" << std::endl;
  for(i = 0; i < dof; i++){
    for(j = 0; j < dof; j++){
      std::cout << c[i*dof + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Stiffness:" << std::endl;
  for(i = 0; i < dof; i++){
    for(j = 0; j < dof; j++){
      std::cout << k[i*dof + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Initial conditions (position, velocity):" << std::endl;
  for(i = 0; i < dof; i++){
    std::cout << x[i] << " " << xd[i];
    std::cout << std::endl;
  }
  
  //compute natural frequency (this only works in 1d)
  double wn = sqrt(k[0]/m[0]);
  //compute critical damping
  double cc = 2.0*m[0]*wn;
  //compute damping ratio
  double zeta = c[0] / cc;
  //compute damped frequency
  double wd = wn*sqrt(1.0-zeta*zeta);
  //define cyclic period
  double tau;
  
  std::cout << "\nThe following values were computed based on the first DOF" << std::endl;
  std::cout << "===========================================================" << std::endl;
  std::cout << "Natural frequency (wn): " << wn << std::endl;
  std::cout << "Damped frequency (wd): " << wd << std::endl;
  std::cout << "Damping ratio (zeta): " << zeta << std::endl;
  std::cout << "Critical damping (Cc): " << cc << std::endl;
  
  //compute number of timesteps to take
  //use natural frequency if oscillatory and damped frequency if not
  if(zeta < 1.0){
    tau = 2*PI/wd;
  }
  else{
    //this is not the true period since these solutions don't oscillate
    //but works to get the timestep
    tau = 2*PI/wn;
  }
  nptp = tau/dt + 1.0;
  
  std::cout << "\nThe following was computed for timesteps:" << std::endl;
  std::cout << "==========================================\n" << std::endl;
  std::cout << "Number of points per period: " << nptp << std::endl;

  return;
}

}
