#include "implicit.h"

namespace STRUCTDYN{

int NewmarkBeta(int dof, int nsteps, double dt, double gamma, double beta, double* m, 
		double* c, double* k, double* x, double* xd, double* xdd, SForces* forces, BC* bcs)
{
  int i, j, l;
  int err = 0;
  //time
  double t;  

  double* keff = new double[dof*dof];
  double* t1 = new double[dof];
  double* t2 = new double[dof];
  double* rhs = new double[dof];
  double* f = new double[dof];

  double* diag = new double[dof];

  double dt2 = dt*dt;
  double c1, c2;
  c1 = 1.0/(beta*dt2);
  c2 = gamma/(beta*dt);

  //Build Keffective matrix
  NewmarkBetaKeff(dof, dt, gamma, beta, m, c, k, keff);

  //initialize acceleration for first time step
  forces->Compute(f, 0.0);
  MatVecMult(c, &xd[0*dof], t2, dof);
  for(j = 0; j < dof; j++){
    t1[j] = f[j] - t2[j];
  }

  MatVecMult(k, &x[0*dof], t2, dof);
  for(j = 0; j < dof; j++){
    t1[j] -= t2[j];
  }

  //compute cholesky factorization
  Cholesky(m, diag, dof);
  CholeskySolve(m, diag, &xdd[0*dof], t1, dof);

  //this is kinda HACKY... watch out...
  //b/c of the way acclerations are computed above.. this method blows
  //up for higher than avg. accelerations (beta/gamma)
  //i.e. the velocity and acceleration terms for static ic's go nuclear
  //zero them if initial velocity is zero
  for(j = 0; j < dof; j++){
    if(xd[0*dof + j] == 0.0){
      xdd[0*dof + j] = 0.0;
    }
  }

  //apply bc's which will remain static to keff matrix
  bcs->ModifyMatrix(keff, dof);

  //compute cholesky factorization for keff that we intend on keeping
  Cholesky(keff, diag, dof);

  for(i = 1; i <= nsteps; i++){
    t = dt*double(i);  

    //Build RHS
    //add mass part
    for(j = 0; j < dof; j++){
      t1[j] = x[(i-1)*dof + j]/(beta*dt2) + xd[(i-1)*dof + j]/(beta*dt) + 
	(1.0/(2.0*beta) - 1.0)*xdd[(i-1)*dof +j];
    }
    //use cholesky matrix vector product since we left m decomposed above
    CholeskyMatVec(m, t1, t2, dof);
    for(j = 0; j < dof; j++){
      rhs[j] = t2[j];
    }
    //add viscous forces
    for(j = 0; j < dof; j++){
      t1[j] = x[(i-1)*dof + j]*gamma/(beta*dt) - xd[(i-1)*dof + j]*(1.0 - gamma/beta) + 
	(gamma/(2.0*beta) - 1.0)*dt*xdd[(i-1)*dof +j];
    }
    MatVecMult(c, t1, t2, dof);
    for(j = 0; j < dof; j++){
      rhs[j] += t2[j];
    }
    //add forcing function
    forces->Compute(f, t);

    for(j = 0; j < dof; j++){
      rhs[j] += f[j];
    }
    
    //apply BC's
    bcs->Apply(rhs, dof);

    //solve for updates
    CholeskySolve(keff, diag, &x[i*dof + 0], rhs, dof);

    //solve for new velocity
    for(j = 0; j < dof; j++){
      xd[i*dof + j] = c2*(x[i*dof + j] - x[(i-1)*dof + j]) + 
			 (1.0 - gamma/beta)*xd[(i-1)*dof + j] -
			 dt*(gamma/(2.0*beta) - 1.0)*xdd[(i-1)*dof + j];
    }
    
    //solve for new acceleration
    for(j = 0; j < dof; j++){
      xdd[i*dof + j] = c1*(x[i*dof + j] - x[(i-1)*dof + j] - xd[(i-1)*dof + j]*dt) -
	(1.0/(2.0*beta) - 1.0)*xdd[(i-1)*dof + j];
    }
  }


  delete [] keff;
  delete [] t1;
  delete [] t2; 
  delete [] rhs;
  delete [] f;
  delete [] diag;

  return err;
}


int WilsonTheta(int dof, int nsteps, double dt, double theta, double* m, 
		double* c, double* k, double* x, double* xd, double* xdd, SForces* forces, BC* bcs)
{
  int i, j, l;
  int err = 0;
  //time
  double t;  

  double* keff = new double[dof*dof];
  double* t1 = new double[dof];
  double* t2 = new double[dof];
  double* rhs = new double[dof];
  double* f1 = new double[dof];
  double* f2 = new double[dof];
  double* diag = new double[dof];

  double c1, c2;

  //Build Keffective matrix
  c1 = 6.0/(theta*dt*theta*dt);
  c2 = 3.0/(theta*dt);
  for(j = 0; j < dof; j++){
    for(l = 0; l < dof; l++){
      keff[j*dof + l] = c1*m[j*dof + l] + c2*c[j*dof + l] + k[j*dof + l]; 
    }
  }

  //compute acceleration for first time step
  forces->Compute(f1, 0.0);
  MatVecMult(c, &xd[0*dof], t2, dof);
  for(j = 0; j < dof; j++){
    t1[j] = f1[j] - t2[j];
  }
  MatVecMult(k, &x[0*dof], t2, dof);
  for(j = 0; j < dof; j++){
    t1[j] -= t2[j];
  }
  //compute cholesky factorization and solve
  Cholesky(m, diag, dof);
  CholeskySolve(m, diag, &xdd[0*dof], t1, dof);

  //this is kinda HACKY... watch out...
  //b/c of the way acclerations are computed above.. this method blows
  //up for higher than avg. accelerations (beta/gamma)
  //i.e. the velocity and acceleration terms for static ic's go nuclear
  //zero them if initial velocity is zero
  for(j = 0; j < dof; j++){
    if(xd[0*dof + j] == 0.0){
      xdd[0*dof + j] = 0.0;
    }
  }

  //apply bc's which will remain static to keff matrix
  bcs->ModifyMatrix(keff, dof);

  //compute cholesky factorization for keff that we intend on keeping
  Cholesky(keff, diag, dof);

  for(i = 1; i <= nsteps; i++){
    t = dt*double(i);  

    //project forcing function
    forces->Compute(f1, t-dt);
    forces->Compute(f2, t);
    for(j = 0; j < dof; j++){
      f1[j] = (1.0 - theta)*f1[j] + theta*f2[j]; 
    }

    //Build RHS
    //add mass part
    for(j = 0; j < dof; j++){
      t1[j] = x[(i-1)*dof + j]*6.0/(theta*dt*theta*dt) + xd[(i-1)*dof + j]*6.0/(theta*dt) + 
	2.0*xdd[(i-1)*dof + j];
    }
    //use cholesky matrix vector multiply since we decomped m earlier
    CholeskyMatVec(m, t1, t2, dof);
    for(j = 0; j < dof; j++){
      rhs[j] = t2[j];
    }
    //add viscous forces
    for(j = 0; j < dof; j++){
      t1[j] = x[(i-1)*dof + j]*3.0/(theta*dt) + 2.0*xd[(i-1)*dof + j] + 
	((theta*dt)/2.0)*xdd[(i-1)*dof + j];
    }
    MatVecMult(c, t1, t2, dof);
    for(j = 0; j < dof; j++){
      rhs[j] += t2[j];
    }
    //add forcing function
    for(j = 0; j < dof; j++){
      rhs[j] += f1[j];
    }

    //apply BC's
    bcs->Apply(rhs, dof);
    
    //solve for updates at n+theta
    CholeskySolve(keff, diag, &x[i*dof + 0], rhs, dof);

    //solve for new acceleration at n+theta
    for(j = 0; j < dof; j++){
      xdd[i*dof + j] = 6.0/(theta*dt*theta*dt)*
	(x[i*dof + j] - x[(i-1)*dof + j] - (theta*dt)*xd[(i-1)*dof + j]) -
	2.0*xdd[(i-1)*dof + j];
    }
    //solve for acceleration at n+1
    for(j = 0; j < dof; j++){
      xdd[i*dof + j] = xdd[(i-1)*dof + j] + (xdd[i*dof + j] - xdd[(i-1)*dof + j])/theta;
    }
    //now update velocity at n+1
    for(j = 0; j < dof; j++){
      xd[i*dof + j] = xd[(i-1)*dof + j] + dt/2.0*(xdd[(i-1)*dof + j] + xdd[i*dof + j]);
    }
    //solve for position at n+1
    for(j = 0; j < dof; j++){
      x[i*dof + j] = x[(i-1)*dof + j] + dt*xd[(i-1)*dof + j] + 
	(1.0/6.0*xdd[i*dof + j] + 1.0/3.0*xdd[(i-1)*dof + j])*dt*dt;
    }

  }

  delete [] keff;
  delete [] t1;
  delete [] t2; 
  delete [] rhs;
  delete [] f1;
  delete [] f2;
  delete [] diag;

  return err;
}


int NewmarkBetaDx(int dof, double dt, double gamma, double beta, double* m, 
		  double* c, double* k, double* dx, double* x_n, double* xd_n, double* xdd_n, 
		  double* f, BC* bcs)
{
  int i, j, l;
  int err = 0;

  double* keff = new double[dof*dof];
  double* t1 = new double[dof];
  double* t2 = new double[dof];
  double* rhs = new double[dof];
  double* p = new double[dof];
  double* diag = new double[dof];
  
  double dt2 = dt*dt;
  double c1, c2;

  //Build Keffective matrix
  NewmarkBetaKeff(dof, dt, gamma, beta, m, c, k, keff);

  //apply bc's which will remain static to keff matrix
  if(bcs != NULL){
    bcs->ModifyMatrix(keff, dof);
  }

  //compute cholesky factorization for keff that we intend on keeping
  Cholesky(keff, diag, dof);

  //Build RHS
  //add mass part
  for(j = 0; j < dof; j++){
    t1[j] = x_n[j]/(beta*dt2) + xd_n[j]/(beta*dt) + 
      (1.0/(2.0*beta) - 1.0)*xdd_n[j];
  }
  // matrix vector product since we left 
  MatVecMult(m, t1, t2, dof);
  for(j = 0; j < dof; j++){
    rhs[j] = t2[j];
  }

  //add viscous forces
  for(j = 0; j < dof; j++){
    t1[j] = x_n[j]*gamma/(beta*dt) - xd_n[j]*(1.0 - gamma/beta) + 
      (gamma/(2.0*beta) - 1.0)*dt*xdd_n[j];
  }
  MatVecMult(c, t1, t2, dof);
  for(j = 0; j < dof; j++){
    rhs[j] += t2[j];
  }

  //add external forces
  for(j = 0; j < dof; j++){
    rhs[j] += f[j];
  }
  
  //subtract unsteady piece
  CholeskyMatVec(keff, x_n, t2, dof);
  for(j = 0; j < dof; j++){
    rhs[j] -= t2[j];
  }

  //apply BC's
  if(bcs != NULL){
    bcs->Apply(rhs, dof);
  }

  //solve for updates dx
  CholeskySolve(keff, diag, dx, rhs, dof);

  delete [] keff;
  delete [] t1;
  delete [] t2; 
  delete [] rhs;
  delete [] p;
  delete [] diag;
  
  return err;
}

void NewmarkBetaUpdateVelAcc(int dof, double dt, double gamma, double beta, double* dx, 
			     double* xd_n, double* xdd_n, double* xd_np1, double* xdd_np1)
{
  int i, j;
  double c1, c2;
  double dt2 = dt*dt;
  c1 = 1.0/(beta*dt2);
  c2 = gamma/(beta*dt);
  
  //solve for new velocity
  for(j = 0; j < dof; j++){
    xd_np1[j] = c2*(dx[j]) +(1.0 - gamma/beta)*xd_n[j] -
      dt*(gamma/(2.0*beta) - 1.0)*xdd_n[j];
  }
  
  //solve for new acceleration
  for(j = 0; j < dof; j++){
    xdd_np1[j] = c1*(dx[j] - xd_n[j]*dt) - (1.0/(2.0*beta) - 1.0)*xdd_n[j];
  }
 
  return;
  
}

void NewmarkBetaKeff(int dof, double dt, double gamma, double beta, double* m, double* c, 
		     double* k, double* keff)
{
  int j, l;
  double dt2 = dt*dt;
  double c1, c2;

  //Build Keffective matrix
  c1 = 1.0/(beta*dt2);
  c2 = gamma/(beta*dt);
  for(j = 0; j < dof; j++){
    for(l = 0; l < dof; l++){
      keff[j*dof + l] = c1*m[j*dof + l] + c2*c[j*dof + l] + k[j*dof + l]; 
    }
  }

  return;
}

}
