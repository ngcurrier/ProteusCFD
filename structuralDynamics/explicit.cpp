#include "explicit.h"

namespace STRUCTDYN{

int RK2ndOrder(int dof, int nsteps, double dt, double* m, double* c, double* k,
	       double* x, double* xd, double* xdd, SForces* forces, BC* bcs, int mode)
{
  int i, j;
  int err = 0;
  //time
  double t;  

  double* k1 = new double[dof*2];
  double* k2 = new double[dof*2];

  double* p = new double[dof];

  //perturbed values for RK
  double tpert;
  double* xpert = new double[dof];
  double* xdpert = new double[dof];

  double a1, a2, p1, q11;
  if(mode == 0){
    //Ralston's method
    a1 = 1.0/3.0;
    a2 = 2.0/3.0;
    p1 = 3.0/4.0;
    q11 = 3.0/4.0;
  }
  else if(mode == 1){
    //trapezoidal rule
    a1 = 1.0/2.0;
    a2 = 1.0/2.0;
    p1 = 1.0;
    q11 = 1.0;
  }
  else{
    std::cerr << "MODE: " << mode << " not recognized" << std::endl;
    return (-1);
  }

  //solve this multi-dimensional 1st orderODE
  //dw/dt = g(t,w)
  //w = [x, xdot]^T
  //g = [xdot, F]^T
  //F = [M]^-1 *{f(t) - [k]x - [c]xdot}

  //since we have to use the exact [M]^-1 this makes this
  //method equivalent to an implicit method
  //factor [M]
  Cholesky(m, p, dof);

  //stores above F
  double* F = new double[dof];

  //stores [xdot(dof x 1), F(dof x 1)]^T
  double* g = new double[dof*2];

  for(i = 0; i < nsteps; i++){
    t = dt*double(i);

    //compute k1
    GetAccelFunction(F, t, m, p, c, k, &x[i*dof], &xd[i*dof], forces, bcs, dof);
    //evaluate g at current time level
    //and F from above
    for(j = 0; j < dof; j++){
      g[0*dof + j] = xd[i*dof + j];
      g[1*dof + j] = F[j];
    }
    //k1 is g(tn, wn), copy it over (that is, it is the derivative of w)
    memcpy(k1, g, sizeof(double)*dof*2);

    //k2 is the same function as above but with the perturbations taken
    //based on k1
    tpert = t+p1*dt;
    memcpy(xpert, &x[i*dof], sizeof(double)*dof);
    memcpy(xdpert, &xd[i*dof], sizeof(double)*dof);
    for(j = 0; j < dof; j++){
      xpert[j] += q11*dt*k1[j];
      xdpert[j] += q11*dt*k1[dof + j];
    }

    //compute k2
    GetAccelFunction(F, t, m, p, c, k, xpert, xdpert, forces, bcs, dof);
    //evaluate g at perturbed time level
    //and F from above
    for(j = 0; j < dof; j++){
      g[0*dof + j] = xdpert[j];
      g[1*dof + j] = F[j];
    }
    //k2 is g(tn + p1*dt, wn + q11*k1*dt), copy it over (that is, it is the derivative of w)
    memcpy(k2, g, sizeof(double)*dof*2);
    
    //now evaluate the next time level
    //w_n+1 = w_n + (a_1*k_1 + a_2*k_2)*dt
    for(j = 0; j < dof; j++){
      x[(i+1)*dof + j] = x[i*dof + j] + (a1*k1[j] + a2*k2[j])*dt;
      xd[(i+1)*dof + j] = xd[i*dof + j] + (a1*k1[dof + j] + a2*k2[dof + j])*dt;
    }
  }

  delete [] k1; 
  delete [] k2;
  
  delete [] xpert;
  delete [] xdpert;
  
  delete [] F;
  delete [] g;
  delete [] p;
  
 return err;
}


int RK3rdOrder(int dof, int nsteps, double dt, double* m, double* c, double* k,
	       double* x, double* xd, double* xdd, SForces* forces, BC* bcs)
{
  int i, j;
  int err = 0;
  //time
  double t;  

  double* k1 = new double[dof*2];
  double* k2 = new double[dof*2];
  double* k3 = new double[dof*2];

  double* p = new double[dof];

  //perturbed values for RK
  double tpert;
  double* xpert = new double[dof];
  double* xdpert = new double[dof];

  //solve this multi-dimensional 1st orderODE
  //dw/dt = g(t,w)
  //w = [x, xdot]^T
  //g = [xdot, F]^T
  //F = [M]^-1 *{f(t) - [k]x - [c]xdot}

  //stores above F
  double* F = new double[dof];

  //stores [xdot(dof x 1), F(dof x 1)]^T
  double* g = new double[dof*2];

  //since we have to use the exact [M]^-1 this makes this
  //method equivalent to an implicit method
  //factor [M]
  Cholesky(m, p, dof);

  for(i = 0; i < nsteps; i++){
    t = dt*double(i);

    //compute k1
    GetAccelFunction(F, t, m, p, c, k, &x[i*dof], &xd[i*dof], forces, bcs, dof);
    //evaluate g at current time level
    //and F from above
    for(j = 0; j < dof; j++){
      g[0*dof + j] = xd[i*dof + j];
      g[1*dof + j] = F[j];
    }
    //k1 is g(tn, wn), copy it over (that is, it is the derivative of w)
    memcpy(k1, g, sizeof(double)*dof*2);

    //k2 is the same function as above but with the perturbations taken
    //based on k1
    tpert = t+0.5*dt;
    memcpy(xpert, &x[i*dof], sizeof(double)*dof);
    memcpy(xdpert, &xd[i*dof], sizeof(double)*dof);
    for(j = 0; j < dof; j++){
      xpert[j] += 0.5*dt*k1[j];
      xdpert[j] += 0.5*dt*k1[dof + j];
    }

    //compute k2
    GetAccelFunction(F, t, m, p, c, k, xpert, xdpert, forces, bcs, dof);
    //evaluate g at perturbed time level
    //and F from above
    for(j = 0; j < dof; j++){
      g[0*dof + j] = xdpert[j];
      g[1*dof + j] = F[j];
    }
    //k2 is g(tn + 0.5*dt, wn + 0.5*k1*dt), copy it over (that is, it is the derivative of w)
    memcpy(k2, g, sizeof(double)*dof*2);
    
    //k3 is the same function as above but with perturbations taken based on k1 and k2
    tpert = t+dt;
    memcpy(xpert, &x[i*dof], sizeof(double)*dof);
    memcpy(xdpert, &xd[i*dof], sizeof(double)*dof);
    for(j = 0; j < dof; j++){
      xpert[j] = xpert[j] - k1[j]*dt + 2.0*dt*k2[j];
      xdpert[j] = xdpert[j] - k1[dof + j]*dt + 2.0*dt*k2[dof + j];
    }

    //compute k3
    GetAccelFunction(F, t, m, p, c, k, xpert, xdpert, forces, bcs, dof);
    //evaluate g at perturbed time level
    //and F from above
    for(j = 0; j < dof; j++){
      g[0*dof + j] = xdpert[j];
      g[1*dof + j] = F[j];
    }
    //k3 is g(tn + dt, wn - k1*dt + 2*k2*dt), copy it over (that is, it is the derivative of w)
    memcpy(k3, g, sizeof(double)*dof*2);


    //now evaluate the next time level
    //w_n+1 = w_n + (1/6*k_1 + 2/3*k_2 + 1/6*k3)*dt
    double c1 = 1.0/6.0;
    double c2 = 2.0/3.0;
    for(j = 0; j < dof; j++){
      x[(i+1)*dof + j] = x[i*dof + j] + (c1*k1[j] + c2*k2[j] + c1*k3[j])*dt;
      xd[(i+1)*dof + j] = xd[i*dof + j] + (c1*k1[dof + j] + c2*k2[dof + j] + c1*k3[dof + j])*dt;
    }
  }

  delete [] k1; 
  delete [] k2;
  delete [] k3; 

  delete [] p;
  
  delete [] xpert;
  delete [] xdpert;
  
  delete [] F;
  delete [] g;
  
 return err;
}

int RK4thOrder(int dof, int nsteps, double dt, double* m, double* c, double* k,
	       double* x, double* xd, double* xdd, SForces* forces, BC* bcs)
{
  int i, j;
  int err = 0;
  //time
  double t;  

  double* k1 = new double[dof*2];
  double* k2 = new double[dof*2];
  double* k3 = new double[dof*2];
  double* k4 = new double[dof*2];

  //perturbed values for RK
  double tpert;
  double* xpert = new double[dof];
  double* xdpert = new double[dof];

  double* p = new double[dof];

  //solve this multi-dimensional 1st orderODE
  //dw/dt = g(t,w)
  //w = [x, xdot]^T
  //g = [xdot, F]^T
  //F = [M]^-1 *{f(t) - [k]x - [c]xdot}

  //stores above F
  double* F = new double[dof];

  //stores [xdot(dof x 1), F(dof x 1)]^T
  double* g = new double[dof*2];

  //since we have to use the exact [M]^-1 this makes this
  //method equivalent to an implicit method
  //factor [M]
  Cholesky(m, p, dof);

  for(i = 0; i < nsteps; i++){
    t = dt*double(i);

    //compute k1
    GetAccelFunction(F, t, m, p, c, k, &x[i*dof], &xd[i*dof], forces, bcs, dof);
    //evaluate g at current time level
    //and F from above
    for(j = 0; j < dof; j++){
      g[0*dof + j] = xd[i*dof + j];
      g[1*dof + j] = F[j];
    }
    //k1 is g(tn, wn), copy it over (that is, it is the derivative of w)
    memcpy(k1, g, sizeof(double)*dof*2);

    //k2 is the same function as above but with the perturbations taken
    //based on k1
    tpert = t+0.5*dt;
    memcpy(xpert, &x[i*dof], sizeof(double)*dof);
    memcpy(xdpert, &xd[i*dof], sizeof(double)*dof);
    for(j = 0; j < dof; j++){
      xpert[j] += 0.5*dt*k1[j];
      xdpert[j] += 0.5*dt*k1[dof + j];
    }

    //compute k2
    GetAccelFunction(F, t, m, p, c, k, xpert, xdpert, forces, bcs, dof);
    //evaluate g at perturbed time level
    //and F from above
    for(j = 0; j < dof; j++){
      g[0*dof + j] = xdpert[j];
      g[1*dof + j] = F[j];
    }
    //k2 is g(tn + 0.5*dt, wn + 0.5*k1*dt), copy it over (that is, it is the derivative of w)
    memcpy(k2, g, sizeof(double)*dof*2);
    
    //k3 is the same function as above but with perturbations taken based on k2
    tpert = t+0.5*dt;
    memcpy(xpert, &x[i*dof], sizeof(double)*dof);
    memcpy(xdpert, &xd[i*dof], sizeof(double)*dof);
    for(j = 0; j < dof; j++){
      xpert[j] += 0.5*dt*k2[j];
      xdpert[j] += 0.5*dt*k2[dof + j];
    }

    //compute k3
    GetAccelFunction(F, t, m, p, c, k, xpert, xdpert, forces, bcs, dof);
    //evaluate g at perturbed time level
    //and F from above
    for(j = 0; j < dof; j++){
      g[0*dof + j] = xdpert[j];
      g[1*dof + j] = F[j];
    }
    //k3 is g(tn + 0.5*dt, wn + 0.5*k2*dt), copy it over (that is, it is the derivative of w)
    memcpy(k3, g, sizeof(double)*dof*2);

    //k4 is the same function as above but with perturbations taken based on k2
    tpert = t+dt;
    memcpy(xpert, &x[i*dof], sizeof(double)*dof);
    memcpy(xdpert, &xd[i*dof], sizeof(double)*dof);
    for(j = 0; j < dof; j++){
      xpert[j] += dt*k3[j];
      xdpert[j] += dt*k3[dof + j];
    }

    //compute k4
    GetAccelFunction(F, t, m, p, c, k, xpert, xdpert, forces, bcs, dof);
    //evaluate g at perturbed time level
    //and F from above
    for(j = 0; j < dof; j++){
      g[0*dof + j] = xdpert[j];
      g[1*dof + j] = F[j];
    }
    //k4 is g(tn + dt, wn + k3*dt), copy it over (that is, it is the derivative of w)
    memcpy(k4, g, sizeof(double)*dof*2);

    //now evaluate the next time level
    //w_n+1 = w_n + (1/6*k1 + 1/3*k2 + 1/3*k3 + 1/6*k4)*dt
    double c1 = 1.0/6.0;
    double c2 = 1.0/3.0;
    for(j = 0; j < dof; j++){
      x[(i+1)*dof + j] = x[i*dof + j] + (c1*k1[j] + c2*k2[j] + c2*k3[j] + c1*k4[j])*dt;
      xd[(i+1)*dof + j] = xd[i*dof + j] + (c1*k1[dof + j] + c2*k2[dof + j] + 
					   c2*k3[dof + j] + c1*k4[dof + j])*dt;
    }
  }

  delete [] k1; 
  delete [] k2;
  delete [] k3; 
  delete [] k4;

  delete [] p;
  
  delete [] xpert;
  delete [] xdpert;
  
  delete [] F;
  delete [] g;
  
 return err;
}

int CentralDiff(int dof, int nsteps, double dt, double* m, double* diagm, double* c, 
		double* k, double* x, double* xd, double* xdd, SForces* forces, BC* bcs)
{

  int i, j, l;
  int err = 0;
  //time
  double t;  

  double* f = new double[dof];
  double* t1 = new double[dof*dof];
  double* t2 = new double[dof*dof];
  double* tv1 = new double[dof];
  double* tv2 = new double[dof];
  double* xnm1 = new double[dof];
  double* p = new double[dof];

  double c1 = 2.0/(dt*dt);
  double c2 = 1.0/(2.0*dt);
  double c3 = 1.0/(dt*dt);

  //method is non-self-starting... bootstrap the first timestep
  //-----------------------------------------------------------//
  t = 0;
  //get forcing function
  forces->Compute(f, t);

  //compute position at time=-1
  //technically x_(-1) = x_0 - dt*xd_0 + (dt*dt/2)*xdd_0 + O(dt^3)
  //also, technically xdd_0 = [M]^-1{f(0) - [k]*x_0 - [c]*xd_0}
  //ignore dt^2 term... BADNESS mitigated
  //this approximation is O(dt^2) in time... not awful
  for(j = 0; j < dof; j++){
    xnm1[j] = x[0*dof + j] - dt*xd[0*dof + j];
  }
  
  //build matrix for update
  for(j = 0; j < dof; j++){
    for(l = 0; l < dof; l++){
      t1[j*dof + l] = c1*m[j*dof + l] - k[j*dof + l];
    }
  }
  //perform matrix vector multiply
  MatVecMult(t1, &x[0*dof + 0], tv1, dof);
  
  //build matrix for update
  for(j = 0; j < dof; j++){
    for(l = 0; l < dof; l++){
      t1[j*dof + l] = c2*c[j*dof + l] - c3*m[j*dof + l];
    }
  }
  //perform matrix vector multiply
  MatVecMult(t1, xnm1, tv2, dof);
  
  //compute RHS.. this is not yet valid for position
  for(j = 0; j < dof; j++){
    tv1[j] = f[j] + tv1[j] + tv2[j];
  }

  //back solve assuming a diagonalized matrix
  //this routine makes that assumption... if NOT valid change this 
  //build the matrix we need
  for(j = 0; j < dof; j++){
    for(l = 0; l < dof; l++){
      t1[j*dof + l] = c2*c[j*dof + l] + c3*m[j*dof + l];
    }
  }
  //back solve for position
  Cholesky(t1, p, dof);
  CholeskySolve(t1, p, &x[(1)*dof + 0], tv1, dof);
  
  //now solve for velocity
  for(j = 0; j < dof; j++){
    xd[(1)*dof + j] = c2*(x[(1)*dof + j] - xnm1[j]);
  }

  //modify the solution matrix with unity on diagonal where IC's are specified
  bcs->ModifyMatrix(diagm, dof);

  for(i = 1; i < nsteps; i++){
    t = dt*double(i);  

    //get forcing function
    forces->Compute(f, t);
      
    //build matrix for update
    for(j = 0; j < dof; j++){
      for(l = 0; l < dof; l++){
	t1[j*dof + l] = c1*diagm[j*dof + l] - k[j*dof + l];
      }
    }
    //perform matrix vector multiply
    MatVecMult(t1, &x[i*dof + 0], tv1, dof);

    //build matrix for update
    for(j = 0; j < dof; j++){
      for(l = 0; l < dof; l++){
	t1[j*dof + l] = c2*c[j*dof + l] - c3*diagm[j*dof + l];
      }
    }
    //perform matrix vector multiply
    MatVecMult(t1, &x[(i-1)*dof + 0], tv2, dof);
    
    //perform update for position
    for(j = 0; j < dof; j++){
      //compute RHS.. this is not yet valid for position
      tv1[j] = f[j] + tv1[j] + tv2[j];
    }
    //back solve assuming a diagonalized matrix
    //this routine makes that assumption... if NOT valid change this 
    //build the matrix we need
    for(j = 0; j < dof; j++){
      for(l = 0; l < dof; l++){
	t1[j*dof + l] = c2*c[j*dof + l] + c3*diagm[j*dof + l];
      }
    }

    //apply bcs
    bcs->Apply(tv1, dof);

    //back solve for position
    MSolve(t1, tv1, &x[(i+1)*dof + 0], dof);
    //now that the position is valid update the velocity
    for(j = 0; j < dof; j++){
      xd[(i+1)*dof + j] = c2*(x[(i+1)*dof + j] - x[(i-1)*dof + j]);
     }
  }

  delete [] f;
  delete [] t1;
  delete [] t2;
  delete [] tv1;
  delete [] tv2;
  delete [] xnm1;
  delete [] p;

  return (err);
}

int HalfStepCentral(int dof, int nsteps, double dt, double* m, double* diagm, double* c,
		    double* k, double* x, double* xd, double* xdd, SForces* forces, BC* bcs)
{

  int i, j, l;
  int err = 0;
  //time
  double t;  

  double* f = new double[dof];
  double* t1 = new double[dof*dof];
  double* t2 = new double[dof*dof];
  double* tv1 = new double[dof];
  double* tv2 = new double[dof];
  double* xnm1 = new double[dof];
  double* xdnmh = new double[dof];
  double* p = new double[dof];

  double c1 = 1.0/(dt*dt);
  double c2 = 1.0/(dt);
  double c3 = 1.0/(2.0*dt);

  //method is non-self-starting... bootstrap the first timestep
  //-----------------------------------------------------------//
  t = 0;
  //get forcing function
  forces->Compute(f, t);

  //compute position at time=-1
  //technically x_(-1) = x_0 - dt*xd_0 + (dt*dt/2)*xdd_0 + O(dt^3)
  //also, technically xdd_0 = [M]^-1{f(0) - [k]*x_0 - [c]*xd_0}
  //ignore dt^2 term... BADNESS mitigated
  //this approximation is O(dt^2) in time... not awful
  for(j = 0; j < dof; j++){
    xnm1[j] = x[0*dof + j] - dt*xd[0*dof + j];
  }
  
  //build matrix for update
  for(j = 0; j < dof; j++){
    for(l = 0; l < dof; l++){
      t1[j*dof + l] = c1*m[j*dof + l] - k[j*dof + l];
    }
  }
  //perform matrix vector multiply
  MatVecMult(t1, &x[0*dof + 0], tv1, dof);
  
  //build matrix for update
  for(j = 0; j < dof; j++){
    for(l = 0; l < dof; l++){
      t1[j*dof + l] = c2*m[j*dof + l] - c[j*dof + l];
    }
  }
  //compute xdot at n-1/2 level
  for(j = 0; j < dof; j++){
    xdnmh[j] = c2*(x[0*dof + j] - xnm1[j]);
  }

  //perform matrix vector multiply
  MatVecMult(t1, xdnmh, tv2, dof);
  
  //compute RHS.. this is not yet valid for position
  for(j = 0; j < dof; j++){
    tv1[j] = f[j] + tv1[j] + tv2[j];
  }

  //back solve assuming a diagonalized matrix
  //this routine makes that assumption... if NOT valid change this 
  //build the matrix we need
  for(j = 0; j < dof; j++){
    for(l = 0; l < dof; l++){
      t1[j*dof + l] = c1*m[j*dof + l];
    }
  }
  //back solve for position
  Cholesky(t1, p, dof);
  CholeskySolve(t1, p, &x[(1)*dof + 0], tv1, dof);
  
  //now solve for velocity
  for(j = 0; j < dof; j++){
    xd[(1)*dof + j] = c3*(x[(1)*dof + j] - xnm1[j]);
  }

  //modify the solution matrix with unity on diagonal where IC's are specified
  bcs->ModifyMatrix(diagm, dof);

  for(i = 1; i < nsteps; i++){
    t = dt*double(i);  

    //get forcing function
    forces->Compute(f, t);
  
    //build matrix for update
    for(j = 0; j < dof; j++){
      for(l = 0; l < dof; l++){
	t1[j*dof + l] = c1*diagm[j*dof + l] - k[j*dof + l];
      }
    }
    //perform matrix vector multiply
    MatVecMult(t1, &x[i*dof + 0], tv1, dof);

    //build matrix for update
    for(j = 0; j < dof; j++){
      for(l = 0; l < dof; l++){
	t1[j*dof + l] = c2*diagm[j*dof + l] - c[j*dof + l];
      }
    }
    //compute xdot at n-1/2 level
    for(j = 0; j < dof; j++){
      xdnmh[j] = c2*(x[i*dof + j] - x[(i-1)*dof + j]);
    }
    //perform matrix vector multiply
    MatVecMult(t1, xdnmh, tv2, dof);
    
    //perform update for position
    for(j = 0; j < dof; j++){
      //compute RHS.. this is not yet valid for position
      tv1[j] = f[j] + tv1[j] + tv2[j];
    }
    //back solve assuming a diagonalized matrix
    //this routine makes that assumption... if NOT valid change this 
    //build the matrix we need
    for(j = 0; j < dof; j++){
      for(l = 0; l < dof; l++){
	t1[j*dof + l] = c1*diagm[j*dof + l];
      }
    }

    //apply bcs
    bcs->Apply(tv1, dof);

    //back solve for position
    MSolve(t1, tv1, &x[(i+1)*dof + 0], dof);
    //now that the position is valid update the velocity
    for(j = 0; j < dof; j++){
      xd[(i+1)*dof + j] = c3*(x[(i+1)*dof + j] - x[(i-1)*dof + j]);
     }
  }

  delete [] f;
  delete [] t1;
  delete [] t2;
  delete [] tv1;
  delete [] tv2;
  delete [] xnm1;
  delete [] xdnmh;
  delete [] p;

  return (err);
}

//this function returns 
//Takes a cholesky decomposed [M] and diagonal vector p from routine Cholesky()
//F = [M]^-1 *{f(t) - [k]x - [c]xdot}
void GetAccelFunction(double* F, double t, double* m, double* p, double* c, double* k, 
		      double* x, double* xd, SForces* forces, BC* bcs, int dof)
{
  int j;
  //stores forcing function
  double* f = new double[dof];
  //stores [k]*pos
  double* kx = new double[dof];
  //stores [c]*vel
  double* cxd = new double[dof];

  //evaluate F at current time level
  forces->Compute(f, t);

  //get [k]*pos
  MatVecMult(k, x, kx, dof);
  //get [c]*vel
  MatVecMult(c, xd, cxd, dof);

  //add the extra terms
  for(j = 0; j < dof; j++){
    f[j] += (-kx[j] - cxd[j]);
  }
  //back solve M
  CholeskySolve(m, p, F, f, dof);

  delete [] f;
  delete [] kx;
  delete [] cxd;

  return;
}

}
