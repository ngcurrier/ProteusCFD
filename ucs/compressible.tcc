#include "param.h"
#include "macros.h"
#include "matrix.h"
#include "solutionSpace.h"
#include "solutionField.h"
#include "pythonInterface.h"

#include <iostream>
#include <cmath>

// NATIVE VARIABLES
// rho
// rho-u
// rho-v
// rho-w
// rho-et
// NON-NATIVE AUX
// T
// P
// u
// v
// w

template <class Type>
CompressibleEqnSet<Type>::CompressibleEqnSet(SolutionSpace<Type>* space, Param<Type>* p)
{
  this->neqn = 5;
  this->nauxvars = 5;
  this->param = p;
  this->space = space;
  //variable set is conservative
  this->varsConservative = 1;

  //set types of variables for output i.e. scalar, vector, etc.
  this->idata = new DataInfo(this->neqn+this->nauxvars, std::string("variableQ"));
  this->idata->AddScalar(0, std::string("Density"));
  this->idata->AddVector(1, std::string("Momentum"));
  this->idata->AddScalar(4, std::string("TotalEnergy"));
  this->idata->AddScalar(5, std::string("Temperature"));
  this->idata->AddScalar(6, std::string("Pressure"));
  this->idata->AddVector(7, std::string("Velocity"));
  this->idata->Verify();
  
  //set gradients required
  this->gdata = new DataInfo((this->neqn+1+3)*3, "gradVariableQ");
  this->gdata->AddVector(0*3, "Grad-Density");
  this->gdata->AddVector(1*3, "Grad-rhou");
  this->gdata->AddVector(2*3, "Grad-rhov");
  this->gdata->AddVector(3*3, "Grad-rhow");
  this->gdata->AddVector(4*3, "Grad-TotalEnergy");
  this->gdata->AddVector(5*3, "Grad-Temperature");
  this->gdata->AddVector(6*3, "Grad-u");
  this->gdata->AddVector(7*3, "Grad-v");
  this->gdata->AddVector(8*3, "Grad-w");
  this->gdata->Verify();
  
  this->Qinf = new Type[this->neqn + this->nauxvars];
 }

 template <class Type>
 CompressibleEqnSet<Type>::~CompressibleEqnSet()
 {
   delete [] this->Qinf;
 }

 template <class Type>
 void CompressibleEqnSet<Type>::InitEqnSet()
 {
   //calculate Qinf values for class variables
   this->UpdateQinf();

   //allocate solution memory 
   this->space->AddField(*this->idata, FIELDS::STATE_TIME, FIELDS::VAR_EVERYWHERE);
   //set pointers in solution space for quick access
   SolutionField<Type> & field = this->space->GetField("variableQ");
   this->space->q = field.GetData(FIELDS::STATE_NP1);
   this->space->qold = field.GetData(FIELDS::STATE_N);
   this->space->qoldm1 = field.GetData(FIELDS::STATE_NM1);

   //allocate solution memory for the gradients we need
   this->space->AddField(*this->gdata, FIELDS::STATE_NONE, FIELDS::VAR_INTERIOR);
   SolutionField<Type> & gfield  = this->space->GetField("gradVariableQ");
   this->space->qgrad = gfield.GetData(FIELDS::STATE_NONE);

   //set Reynold's number - rho*V*d/mu
   Type rho = this->Qinf[0];
   Type V = this->param->GetVelocity(this->space->iter);
   this->param->Re = (rho*this->param->ref_density)*(this->param->ref_velocity*V)*
     this->param->ref_length/this->param->ref_viscosity;
 }

template <class Type>
void CompressibleEqnSet<Type>::RoeFlux(Type QL[], Type QR[], Type avec[], Type vdotn,  
				       Type gamma, Type flux[], Type beta)
{

#if 0
  OneSidedRoeFlux(QL,QR,avec, vdotn, gamma, flux, beta);
  return;
#endif
  
  Int i;
  
  Type Qroe[5];
  Type T[25];
  Type Tinv[25];
  Type eigenvalues[5];
  Type fluxL[5];
  Type fluxR[5];
  Type dQ[5];
  Type dv[5];
  Type dr[5];
  
  //magnitude stored in 4th entry in edge structures area
  Type area = avec[3];
  
  RoeVariables(QL, QR, gamma, Qroe);
  Eigensystem(Qroe, avec, vdotn, gamma, eigenvalues, T, Tinv, beta);
    
  //ENTROPY FIX
  /////////////////////////////////////////////////////////////////////////////
  // adjust the eigenvalues according to Harten and Hyman, 1983 (Hirsch,p. 469)
  // to prevent expansion shocks allowed by the Roe scheme
  // this is the second type suggested by these authors
  ////////////////////////////////////////////////////////////////////////////
  
  Type gm1 = gamma -1.0;
  Type thetaR,thetaL,eigL,eigR,eps,cR,cL,eig;
  const Type rhoL  = QL[0];
  const Type uL   = QL[1]/rhoL;
  const Type vL   = QL[2]/rhoL;
  const Type wL   = QL[3]/rhoL;
  const Type EL    = QL[4];
  const Type vmag2L= uL*uL + vL*vL + wL*wL; 
  const Type PL    = gm1*(EL - 0.5*rhoL*vmag2L);
  
  const Type rhoR  = QR[0];
  const Type uR   = QR[1]/rhoR;
  const Type vR   = QR[2]/rhoR;
  const Type wR   = QR[3]/rhoR;
  const Type ER    = QR[4];
  const Type vmag2R= uR*uR + vR*vR + wR*wR; 
  const Type PR    = gm1*(ER - 0.5*rhoR*vmag2R);
  
  thetaL = uL*avec[0] + vL*avec[1] + wL*avec[2] + vdotn;
  thetaR = uR*avec[0] + vR*avec[1] + wR*avec[2] + vdotn;
  cR = sqrt(gamma*PR/rhoR);
  cL = sqrt(gamma*PL/rhoL);
  
  // eigenvalues 1,2,3 (u)
  eigL  = thetaL;
  eigR  = thetaR;
  eig = eigenvalues[0];
  eps = MAX((eig - eigL),(eigR - eig));
  eps = MAX((Type)0.0,eps);
  if (real(CAbs(eigenvalues[0])) < real(eps)){
    eigenvalues[0] = 0.5*(eigenvalues[0]*eigenvalues[0]/eps + eps);
    eigenvalues[1] = eigenvalues[0];
    eigenvalues[2] = eigenvalues[0];
  }
  else{
    eigenvalues[0] = eigenvalues[1] = eigenvalues[2] = CAbs(eigenvalues[0]);
  }
  
  // eigenvalue 4 (u+c)
  eigL  = thetaL + cL;
  eigR  = thetaR + cR;
  eig = eigenvalues[3];
  eps = MAX((eig - eigL),(eigR - eig));
  eps = MAX((Type)0.0,eps);
  if (real(CAbs(eigenvalues[3])) < real(eps)){
    eigenvalues[3] = 0.5*(eigenvalues[3]*eigenvalues[3]/eps + eps);
  }
  else{
    eigenvalues[3] = CAbs(eigenvalues[3]);
  }
  
  // eigenvalue 5 (u-c)
  eigL  = thetaL - cL;
  eigR  = thetaR - cR;
  eig = eigenvalues[4];
  eps = MAX((eig - eigL),(eigR - eig));
  eps = MAX((Type)0.0,eps);
  if (real(CAbs(eigenvalues[4])) < real(eps)){
    eigenvalues[4] = 0.5*(eigenvalues[4]*eigenvalues[4]/eps + eps);
  }
  else{
    eigenvalues[4] = CAbs(eigenvalues[4]);
  }
  
  //compute difference
  for(i = 0; i < 5; i++){
    dQ[i] = QR[i] - QL[i];
  }
  
  //Multiply left eigenvectors by difference
  MatVecMult(Tinv, dQ, dv, 5);
  
  //scale dv by eigenvalues
  for(i = 0; i < 5; i++){
    dv[i] *= CAbs(eigenvalues[i]);
  }
  
  //Multiply right eigvenvectors by difference
  MatVecMult(T, dv, dr, 5);
  
  //Get standard fluxes on both sides of
  //control surface
  Flux(QL, avec, vdotn, gamma, fluxL, beta);
  Flux(QR, avec, vdotn, gamma, fluxR, beta);
  
  //compute roe flux
  //multiply by area as well at this point
  for(i = 0; i < 5; i++){
    flux[i] = 0.5 * area * (fluxL[i] + fluxR[i] - dr[i]);
  }

  //pseudo-2D simulations
  if(this->param->symmetry2D > 0){
    if(this->param->symmetry2D == 1){
      flux[1] = 0.0;
    }
    else if(this->param->symmetry2D == 2){
      flux[2] = 0.0;
    }
    else if(this->param->symmetry2D == 3){
      flux[3] = 0.0;
    }
  }
}  

template <class Type>
void CompressibleEqnSet<Type>::AUSMFlux(Type QL[], Type QR[], Type avec[], Type vdotn,  
					Type gamma, Type flux[], Type beta)
{
  //From paper: An Enhanced AUSM+-up Scheme for High-Speed Compressible Two-Phase Flows on Hybrid Grids by A.K. Pandare
  //http://www.iccfd.org/iccfd10/papers/ICCFD10-141-Paper.pdf
  
  Type area = avec[3];
  Type gm1 = gamma - 1.0;

  const Type rhoL   = QL[0];
  const Type uL     = QL[1]/rhoL;
  const Type vL     = QL[2]/rhoL;
  const Type wL     = QL[3]/rhoL;
  const Type EL     = QL[4];
  const Type vmag2L = uL*uL + vL*vL + wL*wL; 
  const Type PL     = gm1*(EL - 0.5*rhoL*vmag2L);
  const Type thetaL  = uL*avec[0] + vL*avec[1] + wL*avec[2];
  const Type thetabarL  = uL*avec[0] + vL*avec[1] + wL*avec[2] + vdotn;
  const Type hL = (EL + PL)/rhoL;
  
  const Type rhoR   = QR[0];
  const Type uR     = QR[1]/rhoR;
  const Type vR     = QR[2]/rhoR;
  const Type wR     = QR[3]/rhoR;
  const Type ER     = QR[4];
  const Type vmag2R = uR*uR + vR*vR + wR*wR; 
  const Type PR     = gm1*(ER - 0.5*rhoR*vmag2R);
  const Type thetaR = uR*avec[0] + vR*avec[1] + wR*avec[2];
  const Type thetabarR = uR*avec[0] + vR*avec[1] + wR*avec[2] + vdotn;
  const Type hR = (ER + PR)/rhoR;

  Type cL = sqrt(gamma*PL/rhoL);
  Type cR = sqrt(gamma*PR/rhoR);

  //compute mach numbers
  Type ML = thetabarL/cL;
  Type MR = thetabarR/cR;
  
  // eigenvalues 1,2,3 (u)
  // eigenvalue 4 (u+c)
  // eigenvalue 5 (u-c)

  //AUSM splits the flux into the pressure and convective parts

  //mass flux, mdot_i = a_i * M_i * rho_LR
  const Type V_i = 0.5*(thetabarL + thetabarR); //we use the avg. not sure if this is correct

  Type rho_LR = 0.0;
  if (real(V_i) >= 0.0){
    rho_LR = rhoL;
  }
  else{
    rho_LR = rhoR;
  }

  const Type mdot = rho_LR*V_i; //rho * V
  
  Type* psi = (Type*)alloca(sizeof(Type)*this->neqn);

  if (real(mdot) >= 0.0){
    //pick left state
    psi[0] = 1.0;
    psi[1] = uL;
    psi[2] = vL;
    psi[3] = wL;
    psi[4] = hL;
  }
  else{
    //pick right state
    psi[0] = 1.0;
    psi[1] = uR;
    psi[2] = vR;
    psi[3] = wR;
    psi[4] = hR;
  }


  //Pplus is a function of Mach L
  //Pminus is a funtions of Mach R
  Type M1plus = 0.5*(ML + CAbs(ML));
  Type M2plus = 0.25*(ML + 1.0)*(ML + 1.0);
  Type M1minus = 0.5*(MR - CAbs(MR));
  Type M2minus = 0.25*(MR - 1.0)*(MR - 1.0);
  Type M4plus = 0.0;
  if(real(CAbs(ML)) >= 1.0){
    M4plus  = M1plus;
  }
  else{
    M4plus = M2plus*(1.0 - 2.0*M2minus);
    //M4plus = M2plus*(1.0 + 2.0*M2plus); is this correct, or did I interpret the flux wrong
  }
  Type M4minus = 0.0;
  if(real(CAbs(MR)) >= 1.0){
    M4minus = M1minus;
  }
  else{
    M4minus = M2minus*(1.0 + 2.0*M2plus);
    //M4minus = M2minus*(1.0 - 2.0*M2plus); is this correct, or did I interpret the flux wrong
  }
  
  Type pdiffusion = 0.0;

  Type M = ML; //not sure if I should blend theis here?
  
  //this is all based on the left state
  Type Pplus = 0.0;
  if(real(CAbs(ML)) >= 1.0){
    Pplus = 1.0/ML*M;
  }
  else{
    Pplus = 1.0/ML*M2plus*(2.0 - ML - 3.0*ML*M2minus);
  }

  //this is all based on the right state
  Type Pminus = 0.0;
  if(real(CAbs(MR)) >= 1.0){
    Pminus = 1.0/MR*M;
  }
  else{
    Pminus = 1.0/MR*M2minus*(-2.0 - MR + 3.0*MR*M2plus);
  }

  //determine P based on the wave speeds
  Type P = Pplus*PL + Pminus*PR + pdiffusion;
  
  Type* fluxc = (Type*)alloca(sizeof(Type)*this->neqn);
  Type* fluxp = (Type*)alloca(sizeof(Type)*this->neqn);

  fluxc[0] = mdot*psi[0];
  fluxc[1] = mdot*psi[1];
  fluxc[2] = mdot*psi[2];
  fluxc[3] = mdot*psi[3];
  fluxc[4] = mdot*psi[4];
  
  fluxp[0] = 0.0;
  fluxp[1] = P*avec[0]; //P nx
  fluxp[2] = P*avec[1]; //P ny
  fluxp[3] = P*avec[2]; //P nz
  fluxp[4] = -vdotn*P;
  
  //Get standard fluxes on both sides of
  //control surface
  flux[0] = area*(fluxc[0] + fluxp[0]);
  flux[1] = area*(fluxc[1] + fluxp[1]);
  flux[2] = area*(fluxc[2] + fluxp[2]);
  flux[3] = area*(fluxc[3] + fluxp[3]);
  flux[4] = area*(fluxc[4] + fluxp[4]);
  
  //pseudo-2D simulations
  if(this->param->symmetry2D > 0){
    if(this->param->symmetry2D == 1){
      flux[1] = 0.0;
    }
    else if(this->param->symmetry2D == 2){
      flux[2] = 0.0;
    }
    else if(this->param->symmetry2D == 3){
      flux[3] = 0.0;
    }
  }
}

template <class Type>
void CompressibleEqnSet<Type>::OneSidedRoeFlux(Type* QL, Type* QR, Type* avec, 
					       Type vdotn, Type gamma, Type* flux, Type beta)
{
  Int i;
  Type fluxL[5];
  Type fluxR[5];
  Type t[5];
  Type eig[3];
  Type abseig[3];

  Type gm1 = gamma - 1.0;

  const Type rhoL   = QL[0];
  const Type uL     = QL[1]/rhoL;
  const Type vL     = QL[2]/rhoL;
  const Type wL     = QL[3]/rhoL;
  const Type EL     = QL[4];
  const Type vmag2L = uL*uL + vL*vL + wL*wL; 
  const Type PL     = gm1*(EL - 0.5*rhoL*vmag2L);
  const Type thetaL  = uL*avec[0] + vL*avec[1] + wL*avec[2];
  const Type thetabarL  = uL*avec[0] + vL*avec[1] + wL*avec[2] + vdotn;
  const Type hL = (EL + PL)/rhoL;
  
  const Type rhoR   = QR[0];
  const Type uR     = QR[1]/rhoR;
  const Type vR     = QR[2]/rhoR;
  const Type wR     = QR[3]/rhoR;
  const Type ER     = QR[4];
  const Type vmag2R = uR*uR + vR*vR + wR*wR; 
  const Type PR     = gm1*(ER - 0.5*rhoR*vmag2R);
  const Type thetaR = uR*avec[0] + vR*avec[1] + wR*avec[2];
  const Type thetabarR = uR*avec[0] + vR*avec[1] + wR*avec[2] + vdotn;
  const Type hR = (ER + PR)/rhoR;
  
  // compute the Roe averaged variables
  const Type rho   = sqrt(rhoL*rhoR);
  const Type sigma = rho/(rhoL + rho);
  const Type u     = uL + sigma*(uR - uL);
  const Type v     = vL + sigma*(vR - vL);
  const Type w     = wL + sigma*(wR - wL);
  const Type h     = hL + sigma*(hR - hL);
  const Type vmag2= 0.5*(u*u + v*v + w*w);
  const Type P     = gm1/gamma*rho*(h - vmag2);
  
  Type c2 = gamma*P/rho;
  Type c = sqrt(c2);

  const Type theta = u*avec[0] + v*avec[1] + w*avec[2];
  const Type thetabar = u*avec[0] + v*avec[1] + w*avec[2] + vdotn;

  // Now compute eigenvalues, eigenvectors, and strengths
  eig[0] = thetabar + c;
  eig[1] = thetabar - c;
  eig[2] = thetabar;

  // Limit the eigenvalues / entropy fixes here...
  abseig[0] = abs(eig[0]);
  abseig[1] = abs(eig[1]);
  abseig[2] = abs(eig[2]);

  // Compute primitive variable jumps
  Type drho   = rhoR - rhoL;
  Type dP     = PR - PL;
  Type du     = uR - uL;
  Type dv     = vR - vL;
  Type dw     = wR - wL;
  Type dtheta = thetaR - thetaL;

  // Interface speed of sound squared
  Type xc2 = 1.0/(c*c);

  // Jumps have units of density
  Type dv1 = 0.5*(dP + rho*c*dtheta)*xc2;
  Type dv2 = 0.5*(dP - rho*c*dtheta)*xc2;
  Type dv3 = rho;
  Type dv4 = (c*c*drho - dP)*xc2;

  Type r21 = u + c*avec[0];
  Type r31 = v + c*avec[1];
  Type r41 = w + c*avec[2];
  Type r51 = h + c*theta;

  Type r22 = u - c*avec[0];
  Type r32 = v - c*avec[1];
  Type r42 = w - c*avec[2];
  Type r52 = h - c*theta;

  Type r23 = du - dtheta*avec[0];
  Type r33 = dv - dtheta*avec[1];
  Type r43 = dw - dtheta*avec[2];
  Type r53 = u*du + v*dv + w*dw - theta*dtheta;

  Type r24 = u;
  Type r34 = v;
  Type r44 = w;
  Type r54 = 0.5*vmag2;

  t[0] = abseig[0]*dv1     + abseig[1]*dv2     + abseig[2]*dv4;
  t[1] = abseig[0]*r21*dv1 + abseig[1]*r22*dv2 + abseig[2]*r23*dv3 +abseig[2]*r24*dv4;
  t[2] = abseig[0]*r31*dv1 + abseig[1]*r32*dv2 + abseig[2]*r33*dv3 +abseig[2]*r34*dv4;
  t[3] = abseig[0]*r41*dv1 + abseig[1]*r42*dv2 + abseig[2]*r43*dv3 +abseig[2]*r44*dv4;
  t[4] = abseig[0]*r51*dv1 + abseig[1]*r52*dv2 + abseig[2]*r53*dv3 +abseig[2]*r54*dv4;

  // Compute flux using variables from the left side of face
  fluxL[0] = thetabarL*rhoL;
  fluxL[1] = thetabarL*rhoL*uL + avec[0]*PL;
  fluxL[2] = thetabarL*rhoL*vL + avec[1]*PL;
  fluxL[3] = thetabarL*rhoL*wL + avec[2]*PL;
  fluxL[4] = thetabarL*EL  + thetaL*PL;

  // Compute flux using variables from the right side of face
  fluxR[0] = thetabarR*rhoR;
  fluxR[1] = thetabarR*rhoR*uR + avec[0]*PR;
  fluxR[2] = thetabarR*rhoR*vR + avec[1]*PR;
  fluxR[3] = thetabarR*rhoR*wR + avec[2]*PR;
  fluxR[4] = thetabarR*ER  + thetaR*PR;

  // Compute the contribution to the flux balance
  for(i = 0; i < 5; i++){
    flux[i] = 0.5*avec[3]*(fluxL[i] + fluxR[i] - t[i]);
  }

  //pseudo-2D simulations
  if(this->param->symmetry2D > 0){
    if(this->param->symmetry2D == 1){
      flux[1] = 0.0;
    }
    else if(this->param->symmetry2D == 2){
      flux[2] = 0.0;
    }
    else if(this->param->symmetry2D == 3){
      flux[3] = 0.0;
    }
  }
  
}
 
template <class Type>
 Bool CompressibleEqnSet<Type>::RoeVariables(Type QL[], Type QR[], Type gamma, Type Qroe[])
 {
   //Q assumed to be:
   //roe
   //roe u
   //roe v
   //roe w
   //roe Et

   Type gm1 = gamma - 1.0;

   Type rhoL = QL[0];
   Type rhoR = QR[0];
   Type uL = QL[1]/QL[0];
   Type uR = QR[1]/QR[0];
   Type vL = QL[2]/QL[0];
   Type vR = QR[2]/QR[0];
   Type wL = QL[3]/QL[0];
   Type wR = QR[3]/QR[0];
   Type EL = QL[4];
   Type ER = QR[4];
   Type v2L = uL*uL + vL*vL + wL*wL;
   Type v2R = uR*uR + vR*vR + wR*wR;
   Type PL = gm1*(EL - 0.5*rhoL*v2L);
   Type PR = gm1*(ER - 0.5*rhoR*v2R);
   Type hL = (EL + PL)/rhoL;
   Type hR = (ER + PR)/rhoR;

   //Roe averaged variables
   Type rho = sqrt(rhoL*rhoR);
   Type sigma = rho/(rhoL + rho); 
   Type u = uL + sigma*(uR - uL);
   Type v = vL + sigma*(vR - vL);
   Type w = wL + sigma*(wR - wL);
   Type h = hL + sigma*(hR - hL);
   Type v2h = 0.5*(u*u + v*v + w*w);

   Qroe[0] = rho;
   Qroe[1] = rho*u;
   Qroe[2] = rho*v;
   Qroe[3] = rho*w;
   Qroe[4] = rho/gamma*(h + gm1*v2h);

   return (true);
 }

 template <class Type>
 void CompressibleEqnSet<Type>::Eigensystem(Type* Q, Type* avec, Type vdotn, Type gamma, 
					    Type* eigenvalues, Type* T, Type* Tinv, Type beta)
 {
   Type nx = avec[0];
   Type ny = avec[1];
   Type nz = avec[2];
   Type rho = Q[0];
   Type ru = Q[1];
   Type rv = Q[2];
   Type rw = Q[3];
   Type rE = Q[4];
   Type u = ru/rho;
   Type v = rv/rho;
   Type w = rw/rho;
   Type gm1 = gamma - 1.0;
   //use for eigenvectors
   Type thetaf = u*nx + v*ny + w*nz;
   //use for eigenvalues
   Type theta = thetaf + vdotn;
   Type v2h = 0.5*(u*u + v*v + w*w);
   Type P = gm1*(rE - rho*v2h);
   Type c2 = gamma*P/rho;
   Type c = sqrt(c2);
   
   //Right eigenvectors
   T[0] = nx;
   T[5] = u*nx;
   T[10] = v*nx + rho*nz;
   T[15] = w*nx - rho*ny;
   T[20] = v2h*nx + rho*(v*nz - w*ny);

   T[1] = ny;
   T[6] = u*ny - rho*nz;
   T[11] = v*ny;
   T[16] = w*ny + rho*nx;
   T[21] = v2h*ny + rho*(w*nx - u*nz);

   T[2] = nz;
   T[7] = u*nz + rho*ny;
   T[12] = v*nz - rho*nx;
   T[17] = w*nz;
   T[22] = v2h*nz + rho*(u*ny - v*nx);

   T[3] = rho/c;
   T[8] = rho*(u/c + nx);
   T[13] = rho*(v/c + ny);
   T[18] = rho*(w/c + nz);
   T[23] = rho*(v2h/c + thetaf + c/gm1);

   T[4] = rho/c;
   T[9] = rho*(u/c - nx);
   T[14] = rho*(v/c - ny);
   T[19] = rho*(w/c - nz);
   T[24] = rho*(v2h/c - thetaf + c/gm1);

   //Left eigenvectors
   Tinv[0]  = nx - nz*v/rho + ny*w/rho - nx/c2*v2h*gm1;
   Tinv[1]  = nx/c2*u*gm1;
   Tinv[2]  = nz/rho + nx/c2*v*gm1;
   Tinv[3]  = -ny/rho + nx/c2*w*gm1;
   Tinv[4]  = -nx/c2*  gm1;

   Tinv[5]  = ny + nz*u/rho - nx*w/rho - ny/c2*v2h*gm1;
   Tinv[6]  = -nz/rho + ny/c2*u*gm1;
   Tinv[7]  = ny/c2*v*gm1;
   Tinv[8]  = nx/rho + ny/c2*w*gm1;
   Tinv[9]  = -ny/c2*gm1;

   Tinv[10] = nz - ny*u/rho + nx*v/rho - nz/c2*v2h*gm1;
   Tinv[11] = ny/rho + nz/c2*u*gm1;
   Tinv[12] = -nx/rho + nz/c2*v*gm1;
   Tinv[13] = nz/c2*w*gm1;
   Tinv[14] = -nz/c2*  gm1;

   Tinv[15] = -0.5/rho*(thetaf - gm1*v2h/c);
   Tinv[16] = 0.5/rho*(nx - gm1*u/c);
   Tinv[17] = 0.5/rho*(ny - gm1*v/c);
   Tinv[18] = 0.5/rho*(nz - gm1*w/c);
   Tinv[19] = 0.5/rho*(gm1 /c);

   Tinv[20] = 0.5/rho*(thetaf + gm1*v2h/c);
   Tinv[21] = -0.5/rho*(nx + gm1*u/c);
   Tinv[22] = -0.5/rho*(ny + gm1*v/c);
   Tinv[23] = -0.5/rho*(nz + gm1*w/c);
   Tinv[24] = +0.5/rho*(gm1 /c);

 #if 0
   //check that T and Tinv really are inverses of each other
   Int yes;
   yes = MatInvCheck(T, Tinv, 5);
   if(!yes){
     Abort << "Eigensystem not valid -- inverse check failed!!" << std::endl;
   } 
 #endif

   //eigenvalues
   eigenvalues[0] = theta;
   eigenvalues[1] = theta;
   eigenvalues[2] = theta;
   eigenvalues[3] = theta + c;
   eigenvalues[4] = theta - c;

   return;
 }

 template <class Type>
 void CompressibleEqnSet<Type>::Flux(Type* Q, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta)
 {
   Type rho = Q[0];
   Type ru = Q[1];
   Type rv = Q[2];
   Type rw = Q[3];
   Type rEt = Q[4];
   Type u = ru/rho;
   Type v = rv/rho;
   Type w = rw/rho;
   Type v2h = 0.5*(u*u + v*v + w*w);
   Type gm1 = gamma - 1.0;
   Type P = gm1*(rEt - rho*v2h);
   Type ht = (rEt + P)/rho;

   //add (vdotn) if using moving grid
   Type rhotheta = rho*(avec[0]*u + avec[1]*v + avec[2]*w + vdotn);

   flux[0] = rhotheta;
   flux[1] = (u*rhotheta + P*avec[0]);
   flux[2] = (v*rhotheta + P*avec[1]);
   flux[3] = (w*rhotheta + P*avec[2]);
   flux[4] = (ht*rhotheta - vdotn*P);
}

template <class Type>
void CompressibleEqnSet<Type>::ViscousFlux(Type* Q, Type* grad, Type* avec, Type mut, Type* flux)
{
  //gradient is passed in with velocity components and then temperature
  //parse for readability
  Type ux, uy, uz, vx, vy, vz, wx, wy, wz;
  Type Tx, Ty, Tz;
  Tx = grad[15];
  Ty = grad[16];
  Tz = grad[17];
  ux = grad[18];
  uy = grad[19];
  uz = grad[20];
  vx = grad[21];
  vy = grad[22];
  vz = grad[23];
  wx = grad[24];
  wy = grad[25];
  wz = grad[26];
  
  
  Type u, v, w;
  Type rho = Q[0];
  u = Q[1]/rho;
  v = Q[2]/rho;
  w = Q[3]/rho;

  Type mu = this->ComputeViscosity(Q); 

  Type tmut = (mu + mut);
  Type tauxx, tauyy, tauzz, tauxy, tauxz, tauyz;
  Type fact = -2.0/3.0*(ux + vy + wz);
  Type tauxn, tauyn, tauzn;
  tauxx = 2.0*ux + fact;
  tauyy = 2.0*vy + fact;
  tauzz = 2.0*wz + fact;
  tauxy = uy + vx;
  tauxz = uz + wx;
  tauyz = vz + wy;

  //non-dimensionalize Reynold's number w.r.t speed of sound
  Type V = this->param->GetVelocity(this->space->iter);
  Type Mach = V;
  Type ReTilde = this->param->Re / Mach;
  Type RK = avec[3]/ReTilde;
  Type RKT = RK*tmut;

  //heat flux term
  Type cp = 1.0/(this->param->gamma - 1.0);
  //thermal conductivity
  Type k = mu/this->param->Pr*cp;
  //turbulent thermal conductivity
  Type kT = mut/this->param->PrT*cp;
  Type c1 = -(k + kT);
  Type Tn = Tx*avec[0] + Ty*avec[1] + Tz*avec[2];

  //written as on RHS of the equation
  tauxn = tauxx*avec[0] + tauxy*avec[1] + tauxz*avec[2];
  tauyn = tauxy*avec[0] + tauyy*avec[1] + tauyz*avec[2];
  tauzn = tauxz*avec[0] + tauyz*avec[1] + tauzz*avec[2];

  //this is written as if on the RHS of the equation already
  //i.e. this is the negative of the actual viscous flux
  flux[0] = 0.0;
  flux[1] = -RKT*(tauxn);
  flux[2] = -RKT*(tauyn);
  flux[3] = -RKT*(tauzn);
  flux[4] = -RKT*(tauxn*u + tauyn*v + tauzn*w) + RK*c1*Tn;

  //pseudo-2D simulations
  if(this->param->symmetry2D > 0){
    if(this->param->symmetry2D == 1){
      flux[1] = 0.0;
    }
    else if(this->param->symmetry2D == 2){
      flux[2] = 0.0;
    }
    else if(this->param->symmetry2D == 3){
      flux[3] = 0.0;
    }
  }

  return;
}


template <class Type>
Type CompressibleEqnSet<Type>::MaxEigenvalue(Type* Q, Type* avec, Type vdotn, Type gamma, Type beta)
 {
   Type maxeig, eig4, eig5;
   Type gm1 = gamma - 1.0;
   Type rho = Q[0];
   Type u = Q[1]/rho;
   Type v = Q[2]/rho;
   Type w = Q[3]/rho;
   Type rEt = Q[4];

   Type v2h = 0.5*(u*u + v*v + w*w);
   Type P = gm1*(rEt - rho*v2h);
   Type c2 = gamma*P/rho;
   Type c = sqrt(c2);
   Type theta = GetTheta(Q, avec, vdotn);

   eig4 = theta + c;
   eig5 = theta - c;

   maxeig = MAX(CAbs(eig4), CAbs(eig5));

   return (maxeig);
 }

 template <class Type>
 void CompressibleEqnSet<Type>::UpdateQinf()
 {
   //non-dimensionalization here is by rho_inf, c_inf, T_inf
   
   //use GetVelocity() function for ramping
   Type velocity = this->param->GetVelocity(this->space->iter);
   Type gamma = this->param->gamma;
   Type gm1 = gamma - 1.0;

   //since we don't enforce it directly, check that the param object
   //contains a reasonable speed of sound c^2 = gamma*P/rho
   Type T = this->param->initialTemperature/this->param->ref_temperature; //from input file
   Type p = this->param->initialPressure/this->param->ref_pressure; //from input file
   Type rho = p*gamma/T;
   Type c2avg = gamma*(p/rho);
   Type cavg = sqrt(c2avg);
   Type cavg_dim = cavg * this->param->ref_velocity;

   std::cout << "CompressibleEqnSet - computed speed of sound freestream is " << cavg_dim << std::endl;
   std::cout << "CompressibleEqnSet - freestream reference velocity is set to " << this->param->ref_velocity << " should be equivalent " << std::endl;

   if(Abs(real(cavg_dim - this->param->ref_velocity)) > 1.0e-10){
     std::stringstream ss;
     ss << "CompressibleEqnSet - non-dimensionalization not sane. Please set input variable refVelocity = ";
     ss << cavg_dim;
     Abort << ss.str();
   }
   
   Type u = this->param->flowdir[0]*velocity;
   Type v = this->param->flowdir[1]*velocity;
   Type w = this->param->flowdir[2]*velocity;

   //now set total energy
   Type E = p/gm1 + 0.5*rho*(u*u + v*v + w*w);

   //set class variables to default initialization for Qinf
   this->Qinf[0] = rho;
   this->Qinf[1] = rho*u;
   this->Qinf[2] = rho*v;
   this->Qinf[3] = rho*w;
   this->Qinf[4] = E;

   ComputeAuxiliaryVariables(this->Qinf);
 }

 template <class Type>
 void CompressibleEqnSet<Type>::SetInitialConditions()
 {
   Int i, j;
   Int neqn = this->neqn;
   Int nauxvars = this->nauxvars;
   Mesh<Type>* m = this->space->m;
   Int nnode = m->GetNumNodes();
   Int gnode = m->GetNumParallelNodes();
   Int nbnode = m->GetNumBoundaryNodes();

   //calculate Qinf values for class variables
   this->UpdateQinf();

   std::cout.setf(std::ios::scientific);
   std::cout.precision(6);

   std::cout << "Initializing flow field: " << std::endl;
   std::cout << "========================" << std::endl;
   for(i = 0; i < neqn+nauxvars; i ++){
     std::cout << "\tQinf[" << i << "] = " << this->Qinf[i] 
	       << "\t" << this->idata->GetDofName(i) << std::endl;
   }
   std::cout << "\tReynolds number: " << this->param->Re << std::endl;
   std::cout << std::endl;

   //check for setInitialConditions.py
   std::ifstream pyscript("setInitialConditions.py");
   if(pyscript){
     pyscript.close();
     std::cout << "Attempting to use python interface to set initial conditions" << std::endl;
#ifdef _HAS_PYTHON
     PythonWrapper pywrap("./", "setInitialConditions", "setInitialConditions");
     for(i = 0; i < (nnode+gnode); ++i){
       Type* iq = &this->space->q[i*(neqn+nauxvars)];
       pywrap.SetInitialConditions(this->Qinf, neqn, nauxvars, iq, &m->xyz[i*3]);
       ComputeAuxiliaryVariables(iq);
       if(real(GetPressure(iq)) < 0.0){
	 Abort << "CompressibleEqnSet::SetInitialConditions got a negative pressure - check python script";
       }
       if(real(GetTemperature(iq)) < 0.0){
	 Abort << "CompressibleEqnSet::SetInitialConditions got a negative temperature - check python script";
       }
     }
#else
     Abort << "Python not built with solver";
#endif
   }
   else{
     //set all the nodes interior and phantom
     for(i = 0; i < (nnode+gnode+nbnode); i++){
       for(j = 0; j < neqn; j++){
	 this->space->q[i*(neqn+nauxvars) + j] = this->Qinf[j];
       }
       ComputeAuxiliaryVariables(&this->space->q[i*(neqn+nauxvars)]);
     }
   }
 }

 template <class Type>
 void CompressibleEqnSet<Type>::ApplyDQ(Type* dQ, Type* Q, Type* xyz)
 {
   //specialized implementation b/c we know what variables
   //should be clipped if negative
   Type u, v, w, v2, rho, E;
   Type gm1 = this->param->gamma - 1.0;
   Type minP = 1.0e-10;
   Type minRho = 1.0e-10;
   Type minE = 1.0e-10;
   //clip density to 0.0 if it will be negative
   if(real(Q[0] + dQ[0]) < 0.0){
     std::cerr << "WARNING: density clip";
     std::cerr << "  Coordinates: " << xyz[0] 
	       << " " << xyz[1] << " " << xyz[2] << std::endl;
     Q[0] = minRho;
   }
   else{
     Q[0] += dQ[0];
   }
   //check for negative energy
   if(real(Q[4] + dQ[4]) < 0.0){
     std::cerr << "WARNING: energy clip";
     std::cerr << "  Coordinates: " << xyz[0] 
	       << " " << xyz[1] << " " << xyz[2] << std::endl;

     Q[4] = minE;
   }
   else{
     Q[4] += dQ[4];
   }
   //clip pressure to 0.0 if it will be negative
   //really we are clipping rhoEt to make pressure zero
   rho = Q[0];
   u = Q[1]/rho;
   v = Q[2]/rho;
   w = Q[3]/rho;
   E = Q[4];
   v2 = u*u + v*v + w*w;
   //pressure = gm1*(E - 0.5*rho*v2);
   if(real(E) < real(0.5*rho*v2)){
     std::cerr << "WARNING: pressure clip";
     std::cerr << "  Coordinates: " << xyz[0] 
	       << " " << xyz[1] << " " << xyz[2] << std::endl;
     std::cerr << "Dynamic pressure is : " << 0.5*rho*v2 << std::endl;
     Type v2mod = 2.0*(E - minP/gm1);
     Type frac = 0.0;
     if(real(v2mod) > 0.0){
       frac = sqrt(v2mod/v2);
     }
     //modify velocities based on fraction of minimum
     //velocities needed
     u *= frac;
     v *= frac;
     w *= frac;
     Q[1] = rho*u;
     Q[2] = rho*v;
     Q[3] = rho*w;
   }
   else{
     Q[1] += dQ[1];
     Q[2] += dQ[2];
     Q[3] += dQ[3];
   }
   ComputeAuxiliaryVariables(Q);
 }

template <class Type>
Type CompressibleEqnSet<Type>::GetTheta(Type* Q, Type* avec, Type vdotn)
{
  Type nx = avec[0];
  Type ny = avec[1];
  Type nz = avec[2];
  Type u = Q[1]/Q[0];
  Type v = Q[2]/Q[0];
  Type w = Q[3]/Q[0];
  
  return (u*nx + v*ny + w*nz + vdotn);
}

template <class Type>
Int CompressibleEqnSet<Type>::GetGradientsLocation(std::vector<Int>& gradientLoc)
{
  gradientLoc.clear();
  //density
  gradientLoc.push_back(0);
  //momentum
  gradientLoc.push_back(1);
  gradientLoc.push_back(2);
  gradientLoc.push_back(3);
  //total energy
  gradientLoc.push_back(4);
  //temperature 
  gradientLoc.push_back(5);
  //velocity for viscous terms
  gradientLoc.push_back(7);
  gradientLoc.push_back(8);
  gradientLoc.push_back(9);
  return gradientLoc.size();
}


// recall that gradients are for primitive (rho,u,v,w,rhoEt) and not conservative
// so we cannot extrpolate them directly, do some magic -- POOF!
// Qho - higher order Q 
// q - lower order Q 
// dQedge - deltaQ across the face (native, i.e. qL-qR or qR-qL)
// gradQ - gradient of lower order Q 
// dx - vector distance from CV center to face
// limiter - value 0 to 1 which limits extrapolation 
template <class Type>
void CompressibleEqnSet<Type>::ExtrapolateVariables(Type* Qho, const Type* q, const Type* dQedge, 
						    const Type* gradQ, const Type* dx, const Type* limiter)
{
  //density
  Qho[0] = q[0] + this->ExtrapolateCorrection(dQedge[0], &gradQ[0*3], dx)*limiter[0];
  //momentum
  #pragma ivdep
  for(Int i = 1; i < 4; i++){
    //get higher order ru,rv,rw
    Qho[i] = q[i] + this->ExtrapolateCorrection(dQedge[i], &gradQ[i*3], dx)*limiter[i];
  }
  //total energy
  Qho[4] = q[4] + this->ExtrapolateCorrection(dQedge[4], &gradQ[4*3], dx)*limiter[4];
}

//Q - native variables
template <class Type>
Int CompressibleEqnSet<Type>::BadExtrapolation(Type* Q)
{
  //check for negative pressure, density, and total energy
  Type r = Q[0];
  Type u = Q[1]/r;
  Type v = Q[2]/r;
  Type w = Q[3]/r;
  Type E = Q[4];
  Type gamma = this->param->gamma;
  Type v2h = 0.5*(u*u + v*v + w*w);
  Type gm1 = gamma - 1.0;
  
  Type p = gm1*(E - r*v2h);
  if(real(p) < 1.0e-10){
    return true;
  }
  if(real(r) < 0.0){
    return true;
  }
  if(real(E) < 1.0e-10){
    return true;
  }
  return false;
}

template <class Type>
Type CompressibleEqnSet<Type>::GetDensity(Type* Q)
{
  return (Q[0]);
}

//This returns absolute pressure
template <class Type>
Type CompressibleEqnSet<Type>::ComputePressure(Type* Q, Type gamma){
  Type p;
  Type r = Q[0];
  Type u = Q[1]/r;
  Type v = Q[2]/r;
  Type w = Q[3]/r;
  Type E = Q[4];
  Type v2h = 0.5*(u*u + v*v + w*w);
  Type gm1 = gamma - 1.0;
  
  p = gm1*(E - r*v2h);
  return (p);
}

template <class Type>
Type CompressibleEqnSet<Type>::ComputeTemperature(Type* Q, Type gamma)
{
  Type P = ComputePressure(Q, gamma);
  Type rho = Q[0];
  
  return (gamma*P/rho);
}

template <class Type>
void CompressibleEqnSet<Type>::ComputeStressVector(Type* vgrad, Type* avec, Type mu, Type* stress)
{
  //gradient is passed in with velocity components
  //parse for readability
  Type ux, uy, uz, vx, vy, vz, wx, wy, wz;
  ux = vgrad[3];
  uy = vgrad[4];
  uz = vgrad[5];
  vx = vgrad[6];
  vy = vgrad[7];
  vz = vgrad[8];
  wx = vgrad[9];
  wy = vgrad[10];
  wz = vgrad[11];

  Type tauxx, tauyy, tauzz, tauxy, tauxz, tauyz;
  Type div = -2.0/3.0*(ux + vy + wz);
  Type tauxn, tauyn, tauzn;
  tauxx = 2.0*ux + div;
  tauyy = 2.0*vy + div;
  tauzz = 2.0*wz + div;
  tauxy = uy + vx;
  tauxz = uz + wx;
  tauyz = vz + wy;

  //written as on RHS of the equation
  tauxn = tauxx*avec[0] + tauxy*avec[1] + tauxz*avec[2];
  tauyn = tauxy*avec[0] + tauyy*avec[1] + tauyz*avec[2];
  tauzn = tauxz*avec[0] + tauyz*avec[1] + tauzz*avec[2];

  Type V = this->param->GetVelocity(this->space->iter);
  Type Mach = V;
  Type ReTilde = this->param->Re / Mach;

  //negative sign indicates stress on the body, not the flow
  stress[0] = -(mu/ReTilde)*tauxn;
  stress[1] = -(mu/ReTilde)*tauyn;
  stress[2] = -(mu/ReTilde)*tauzn;

  return;
}

template <class Type>
Type CompressibleEqnSet<Type>::GetPressure(Type* Q)
{
  return (Q[6]);
}

template <class Type>
Type CompressibleEqnSet<Type>::GetTemperature(Type* Q)
{
  return (Q[5]);
}

template <class Type>
Type CompressibleEqnSet<Type>::GetGamma(Type* Q)
{
   return (this->param->gamma);
}

template <class Type>
Type CompressibleEqnSet<Type>::GetCp(Type* Q, Type gamma)
{
  Type P = GetPressure(Q);
  Type V = this->param->GetVelocity(this->space->iter);
  Type Mach = V;
  return ((P - 1.0/gamma)/(0.5*Mach*Mach));
}

template <class Type>
Type CompressibleEqnSet<Type>::GetCf(Type tauw, Type rho)
{
  Type V = this->param->GetVelocity(this->space->iter);
  Type Mach = V;
  return (tauw/(0.5*rho*Mach*Mach));
}

template <class Type>
Type CompressibleEqnSet<Type>::GetRe()
{
  //this function is overloaded b/c we nondimensionalize by
  //the speed of sound in the compressible eqnset, which makes the
  //mach number show up everytime Re is used
  Type V = this->param->GetVelocity(this->space->iter);
  Type Mach = V;
  return (this->param->Re / Mach);
}

template <class Type>
void CompressibleEqnSet<Type>::SourceTerm(Type* Q, Type vol, Type* source)
{
  MemBlank(source, this->neqn);
  if(this->param->gravity_on){
    Type rho = Q[0];
    Type g = this->param->gravity;
    g /= (this->param->ref_length/(this->param->ref_time*this->param->ref_time));
    source[1] += g*(rho-1.0)*vol*this->param->gravdir[0];
    source[2] += g*(rho-1.0)*vol*this->param->gravdir[1];
    source[3] += g*(rho-1.0)*vol*this->param->gravdir[2];
  }
}


template <class Type>
Int CompressibleEqnSet<Type>::GetVelocityGradLocation()
{
  return 1;
}

template <class Type>
Int CompressibleEqnSet<Type>::GetMomentumLocation()
{
  return 1;
}


template <class Type>
void CompressibleEqnSet<Type>::ComputeAuxiliaryVariables(Type* Q)
{
  Type gm1 = this->param->gamma - 1.0;
  Type u = Q[1] / Q[0];
  Type v = Q[2] / Q[0];
  Type w = Q[3] / Q[0];
  Type V2 = u*u + v*v + w*w;
  Type P = gm1 * (Q[4] - 0.5 * Q[0] * V2);
  Q[6] = P;
  Q[5] = ComputeTemperature(Q, GetGamma(Q));
  Q[7] = u;
  Q[8] = v;
  Q[9] = w;
}

template <class Type>
void CompressibleEqnSet<Type>::GetFarfieldBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta)
{
  //
  //see reference: Whitfield AIAA 1984 1552
  //for CVBC and impermeable wall BC development
  //
  Int i;
  Int neqn = this->neqn;
  Type* tempspace = (Type*)alloca(sizeof(Type)*neqn);
  Type* QLmod = (Type*)alloca(sizeof(Type)*neqn);
  Type gamma = this->param->gamma;

  //subiteration for design derivatives
  for(Int subit = 0; subit < 10; subit++){
    
    memcpy(QLmod, QL, sizeof(Type)*neqn);
    //extract grid speeds in normal direction
    //vdotn = -dot(gridspeeds, avec)
    Type u, v, w;
    Type ru, rv, rw;
    Type rhoi = QLmod[0];
    u = vdotn*avec[0];
    v = vdotn*avec[1];
    w = vdotn*avec[2];
    ru = u*rhoi;
    rv = v*rhoi;
    rw = w*rhoi;
    
    //add grid speeds to the velocity component of the internal Q
    QLmod[1] += ru;
    QLmod[2] += rv;
    QLmod[3] += rw;
    
    //average variables
    for(i = 0; i < neqn; i++){
      tempspace[i] = (QL[i] + QR[i])/2.0;
    }
    
    Type theta = GetTheta(tempspace, avec, vdotn);
    
    //exterior data
    Type rhob = QR[0];
    Type ub = QR[1]/QR[0];
    Type vb = QR[2]/QR[0];
    Type wb = QR[3]/QR[0];
    Type pb; 
    
    //interior data
    rhoi = QL[0];
    Type ui = QL[1]/QL[0] + u;
    Type vi = QL[2]/QL[0] + v;
    Type wi = QL[3]/QL[0] + w;
    Type pi = GetPressure(QL);
    
    //Freestream data
    Type rhoinf = Qinf[0];
    Type uinf = Qinf[1]/Qinf[0] + u;
    Type vinf = Qinf[2]/Qinf[0] + v;
    Type winf = Qinf[3]/Qinf[0] + w;
    Type pinf = this->GetPressure(Qinf);
    
    //use average state to calculate speed of sound
    Type pavg = ComputePressure(tempspace, gamma);
    Type rhoavg = GetDensity(tempspace);
    Type c2avg = gamma*(pavg/rhoavg);
    Type cavg = sqrt(c2avg);
    
    //Normals
    Type nx = avec[0];
    Type ny = avec[1];
    Type nz = avec[2];
    
    //no flow
    if(theta == 0.0){
      //do nothing
      return;
    }
    //supersonic outflow
    else if(real(theta) > 0.0 && real(CAbs(theta/cavg)) >= 1.0){
      for(i = 0; i < neqn; i++){
	QR[i] = QL[i];
      }
    }
    //supersonic inflow
    else if(real(theta) < 0.0 && real(CAbs(theta/cavg)) >= 1.0){
      for(i = 0; i < neqn; i++){
	QR[i] = Qinf[i];
	//TODO: investigate the effects of this hardset
	//QR[i] = QL[i] = Qinf[i];
      }
    }
    //subsonic outflow
    else if(real(theta) > 0.0 && real(CAbs(theta/cavg)) < 1.0){
      //TODO: move to eqnset class
      pb = pinf;
      rhob = rhoi + (pb - pi) / c2avg;
      Type temp = (pb - pi) / (rhoavg * cavg);
      ub = ui - nx * temp;
      vb = vi - ny * temp;
      wb = wi - nz * temp;
      QR[0] = rhob;
      QR[1] = rhob*ub;
      QR[2] = rhob*vb;
      QR[3] = rhob*wb;
      QR[4] = pb/(gamma - 1.0) + 0.5*rhob*(ub*ub + vb*vb + wb*wb);
    }
    //subsonic inflow
    else if(real(theta) < 0.0 && real(CAbs(theta/cavg)) < 1.0){
      //TODO: move to eqnset class
      pb = 0.5*(pinf + pi + rhoavg*cavg*(nx*(uinf - ui) + ny*(vinf - vi) + nz*(winf - wi)) );
      rhob = rhoinf + (pb - pinf) / c2avg;
      Type temp = (pb - pinf) / (rhoavg * cavg);
      ub = uinf + nx * temp;
      vb = vinf + ny * temp;
      wb = winf + nz * temp;
      QR[0] = rhob;
      QR[1] = rhob*ub;
      QR[2] = rhob*vb;
      QR[3] = rhob*wb;
      QR[4] = pb/(gamma - 1.0) + 0.5*rhob*(ub*ub + vb*vb + wb*wb);
    }
    else{
      //std::cerr << "CVBC error -- half edge id " << eid << "of type " <<  BCs[bcid] << std::endl;
      return;
    }
  }
}

template <class Type>
void CompressibleEqnSet<Type>::GetInviscidWallBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta)
{ 
  Int i;
  Int neqn = this->neqn;
  Type* tempspace = (Type*)alloca(sizeof(Type)*neqn);
  Type* QLmod = (Type*)alloca(sizeof(Type)*neqn);
  Type gamma = this->param->gamma;

  for(Int subit =  0; subit < 10; subit++){
    memcpy(QLmod, QL, sizeof(Type)*neqn);
    
    //extract grid speeds in normal direction
    //vdotn = -dot(gridspeeds, avec)
    Type u, v, w;
    Type ru, rv, rw;
    Type rhoi = QLmod[0];
    u = vdotn*avec[0];
    v = vdotn*avec[1];
    w = vdotn*avec[2];
    ru = u*rhoi;
    rv = v*rhoi;
    rw = w*rhoi;
    
    //add grid speeds to the velocity component of the internal Q
    QLmod[1] += ru;
    QLmod[2] += rv;
    QLmod[3] += rw;
    
    //average variables
    for(i = 0; i < neqn; i++){
      tempspace[i] = (QL[i] + QR[i])/2.0;
    }
    
    //exterior data
    Type rhob = QR[0];
    Type ub = QR[1]/QR[0];
    Type vb = QR[2]/QR[0];
    Type wb = QR[3]/QR[0];
    Type pb; 
    
    //interior data
    rhoi = QL[0];
    Type ui = QL[1]/QL[0] + ru;
    Type vi = QL[2]/QL[0] + rv;
    Type wi = QL[3]/QL[0] + rw;
    Type pi = GetPressure(QL);
    
    //use average state to calculate speed of sound
    Type rhoavg = GetDensity(tempspace);
    //add gridspeed contributions to velocity
    tempspace[1] += rhoavg*u;
    tempspace[2] += rhoavg*v;
    tempspace[3] += rhoavg*w;
    Type pavg = ComputePressure(tempspace, gamma);
    Type c2avg = gamma*(pavg/rhoavg);
    Type cavg = sqrt(c2avg);
    
    if(!this->param->no_cvbc){
      //Normals
      Type nx = avec[0];
      Type ny = avec[1];
      Type nz = avec[2];
      
      pb = pi + rhoavg*cavg*(this->GetTheta(QL, avec, vdotn));
      rhob = rhoi + (pb - pi) / c2avg;
      Type temp = (pb - pi) / (rhoavg*cavg);
      ub = ui - nx * temp;
      vb = vi - ny * temp;
      wb = wi - nz * temp;
      
      QR[0] = rhob;
      QR[1] = rhob*ub;
      QR[2] = rhob*vb;
      QR[3] = rhob*wb;
      QR[4] = pb/(gamma - 1.0) + 0.5*rhob*(ub*ub + vb*vb + wb*wb);
      
    }
    else{
      for(i = 0; i < neqn; i++) QR[i] = QLmod[i];
      //mirror the velocities
      MirrorVector(&QLmod[1], avec, &QR[1]);
    }
    
#if 0
    //check for leakage through the surface
    std::cout << std::endl;
    for(i = 0; i < neqn; i++){
      tempspace[i] = (QL[i] + QR[i])/2.0;
    }
    ui = tempspace[1]/tempspace[0];
    vi = tempspace[2]/tempspace[0];
    wi = tempspace[3]/tempspace[0];
    std::cout << "theta " << eqnset->GetTheta(tempspace, avec, vdotn) << std::endl;
    std::cout << std::endl;
#endif
    
  }
}

template <class Type>
void CompressibleEqnSet<Type>::GetInternalInflowBoundaryVariables(Type* QL, Type* QR, Type* Qinf, 
								  Type pressure, Type* densities, 
								  Type* flowDirection, Type velocity)
{
  //Here we are imposing the backpressure for the subsonic inflow case
  //hardsetting the velocity and just allow the pressure to float
  //this is not really characteristics but maybe a bit better than a hard set

  //Ttotal = Tstatic + (M*c)^2/(2*cp)
  //Ptotal = Pstatic*(1 + (gamma -1)/2*M^2)^(gamma/(gamma-1))
  //c = sqrt(Tstatic*cp*(gamma-1))

  Type gamma = this->param->gamma;

  Type pi = ComputePressure(QL, gamma);
  Type rhoi = QL[0];
  Type u, v, w;
  u = velocity*flowDirection[0];
  v = velocity*flowDirection[1];
  w = velocity*flowDirection[2];

  QR[1] = rhoi*u;
  QR[2] = rhoi*v;
  QR[3] = rhoi*w;

  //let both density and pressure float to the interior values
  QR[0] = QL[0];
  //QR[0] = Qinf[0];
  QR[4] = pi/(gamma - 1.0) + 0.5*rhoi*(u*u + v*v + w*w);
  //QR[4] = pressure;
}

template <class Type>
void CompressibleEqnSet<Type>::GetInternalOutflowBoundaryVariables(Type* QL, Type* QR, Type pressure, Type gamma)
{

  QR[0] = QL[0];
  QR[1] = QL[1];
  QR[2] = QL[2];
  QR[3] = QL[3];
  QR[4] = pressure/(gamma - 1.0) + 0.5*(QR[1]*QR[1] + QR[2]*QR[2] + QR[3]*QR[3])/QR[0];
  return;
}

template <class Type>
void CompressibleEqnSet<Type>::GetViscousWallBoundaryVariables(Type* QL, Type* QR, Type* vel, Int normalNode, 
							       Type* normalQ, Type Twall)
{
  //adiabatic
  if(real(Twall) < 0.0){
    //see notes on this for more detail
    //chain rule specifies that for adiabatic wall both density and total energy
    //are imposed from the field
    QR[0] = QL[0] = normalQ[0];
    QR[4] = QL[4] = normalQ[4];
  }
  //constant temp
  else{
    //Constant wall temp implies total energy set from the state equation
    Type v2 = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
    Type gamma = this->param->gamma;
    Type gm1 = gamma - 1.0;
    Type rhoEt = (Twall*QL[0]/(gamma*gm1) + 0.5*QR[0]*v2);
    QR[0] = QL[0];
    QR[4] = QL[4] = rhoEt;
  }

  //this is a hard set of zero on the wall since the internal nodes lie on the physical
  //surface...
  QR[1] = QL[1] = QL[0]*vel[0];
  QR[2] = QL[2] = QL[0]*vel[1];
  QR[3] = QL[3] = QL[0]*vel[2];

  return;
}

template <class Type>
void CompressibleEqnSet<Type>::ModifyViscousWallJacobian(Type* QL, Type* QR, 
							 Type* vel, Int cvid, CRS<Type>* crs,
							 Int normalNode, Type Twall)
{
  //blank all the subrows for ru, rv, rw components
  crs->A->BlankSubRow(cvid, 1);
  crs->A->BlankSubRow(cvid, 2);
  crs->A->BlankSubRow(cvid, 3);

  //modify diagonal matrix as described in Anderson's Compute Fluids Vol. 23 paper
  //"An Implicit Upwind Algorithm for Computing Turbulent Flows on Unstructured Grids"
  Type* diag = crs->A->GetPointer(cvid, cvid);
  Type gamma = this->param->gamma;
  
  //adiabatic
  if(real(Twall) < 0.0){
    crs->A->BlankSubRow(cvid, 0);
    crs->A->BlankSubRow(cvid, 4);
    Type* offJac = crs->A->GetPointer(cvid, normalNode);
    offJac[0*this->neqn + 0] = -1.0;
    offJac[4*this->neqn + 4] = -1.0;
  }
  //constant wall temperature
  else{
    Type v2 = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
    crs->A->BlankSubRow(cvid, 4);
    //this is the first value in the last row
    //we relate the change in energy to the change in density directly
    diag[4*this->neqn + 0] = -(Twall/(gamma*(gamma-1.0)) + 0.5*v2);
  }

  return;
}

template <class Type>
void CompressibleEqnSet<Type>::ModifyViscousWallResidual(Type* res, Type* vel, Int normalNode, Type Twall)
{
  //zero the rows we are setting explicitly, namely ru,rv,rw
  //this has the effect of not allowing du, dv, dw to change during the
  //linear system solve
  res[1] = 0.0;
  res[2] = 0.0;
  res[3] = 0.0;

  //this is done in Anderson's paper, d(rho*Et) ~ drho
  res[4] = 0.0;

  //if the wall is adiabatic, we zero the density term
  if(real(Twall) < 0.0){
    res[0] = 0.0;
  }

  return;
}

template <class Type>
void CompressibleEqnSet<Type>::ViscousJacobian(Type* QL, Type* QR, Type* dx, Type s2, 
					       Type* avec, Type mut, Type* aL, Type* aR)
{
  Int i;
  Int neqn = this->neqn;
  Type* Qavg = (Type*)alloca((neqn+this->nauxvars)*sizeof(Type));
  for(i = 0; i < neqn; i++){
    Qavg[i] = 0.5*(QL[i] + QR[i]);
  }
  ComputeAuxiliaryVariables(Qavg);
  Type mu = this->ComputeViscosity(Qavg);
  Type tmut = (mu + mut);

  //non-dimensionalize Reynold's number w.r.t speed of sound
  Type ReTilde = this->param->Re / this->param->GetVelocity(this->space->iter);
  Type RK = avec[3]/ReTilde;
  Type RKT = RK*tmut;
  
  //compute some partial derivatives we need to continue
  //-------------------------------------------------------------------------
  
  //first row, we need the partial of each velocity derivative Ux, Vy, Wz, ... 
  //w.r.t density

  //build some of the coefficients into the distance variables
  Type* DxL = (Type*)alloca(sizeof(Type)*3);
  Type* DxR = (Type*)alloca(sizeof(Type)*3);

  Type rhoL = QL[0];
  Type uL = QL[1]/rhoL;
  Type vL = QL[2]/rhoL;
  Type wL = QL[3]/rhoL;
  Type PL = GetPressure(QL);

  Type rhoR = QR[0];
  Type uR = QR[1]/rhoR;
  Type vR = QR[2]/rhoR;
  Type wR = QR[3]/rhoR;
  Type PR = GetPressure(QR);

  Type u = 0.5*(uL + uR);
  Type v = 0.5*(vL + vR);
  Type w = 0.5*(wL + wR);

  for(i = 0; i < 3; i++){
    //these are the directional derivative pieces
    DxL[i] = -dx[i]/(rhoL*s2);
    DxR[i] = dx[i]/(rhoR*s2);
  }

  Type dux_drL = -u*DxL[0];
  Type duy_drL = -u*DxL[1];
  Type duz_drL = -u*DxL[2];
  Type dvx_drL = -v*DxL[0];
  Type dvy_drL = -v*DxL[1];
  Type dvz_drL = -v*DxL[2];
  Type dwx_drL = -w*DxL[0];
  Type dwy_drL = -w*DxL[1];
  Type dwz_drL = -w*DxL[2];

  //compute the divergence derivative
  Type dfact_drL = -2.0/3.0*(dux_drL + dvy_drL + dwz_drL);
  Type dtauxx_drL = (2.0*dux_drL + dfact_drL);
  Type dtauyy_drL = (2.0*dvy_drL + dfact_drL);
  Type dtauzz_drL = (2.0*dwz_drL + dfact_drL);
  Type dtauxy_drL = duy_drL + dvx_drL;
  Type dtauxz_drL = duz_drL + dwx_drL;
  Type dtauyz_drL = dvz_drL + dwy_drL;

  Type dtauxn_drL = dtauxx_drL*avec[0] + dtauxy_drL*avec[1] + dtauxz_drL*avec[2];
  Type dtauyn_drL = dtauxy_drL*avec[0] + dtauyy_drL*avec[1] + dtauyz_drL*avec[2];
  Type dtauzn_drL = dtauxz_drL*avec[0] + dtauyz_drL*avec[1] + dtauzz_drL*avec[2];

  Type dux_drR = -u*DxR[0];
  Type duy_drR = -u*DxR[1];
  Type duz_drR = -u*DxR[2];
  Type dvx_drR = -v*DxR[0];
  Type dvy_drR = -v*DxR[1];
  Type dvz_drR = -v*DxR[2];
  Type dwx_drR = -w*DxR[0];
  Type dwy_drR = -w*DxR[1];
  Type dwz_drR = -w*DxR[2];

  //compute the divergence derivative
  Type dfact_drR = -2.0/3.0*(dux_drR + dvy_drR + dwz_drR);
  Type dtauxx_drR = (2.0*dux_drR + dfact_drR);
  Type dtauyy_drR = (2.0*dvy_drR + dfact_drR);
  Type dtauzz_drR = (2.0*dwz_drR + dfact_drR);
  Type dtauxy_drR = duy_drR + dvx_drR;
  Type dtauxz_drR = duz_drR + dwx_drR;
  Type dtauyz_drR = dvz_drR + dwy_drR;

  Type dtauxn_drR = dtauxx_drR*avec[0] + dtauxy_drR*avec[1] + dtauxz_drR*avec[2];
  Type dtauyn_drR = dtauxy_drR*avec[0] + dtauyy_drR*avec[1] + dtauyz_drR*avec[2];
  Type dtauzn_drR = dtauxz_drR*avec[0] + dtauyz_drR*avec[1] + dtauzz_drR*avec[2];

  //other constants
  Type c43 = 4.0/3.0;
  Type mc23 = -2.0/3.0;

  Type dR2_drhou, dR2_drhov, dR2_drhow;
  Type dR3_drhou, dR3_drhov, dR3_drhow;
  Type dR4_drhou, dR4_drhov, dR4_drhow;

  //first row   (dR1/dQr)
  aR[0] = 0.0;
  aR[1] = 0.0;
  aR[2] = 0.0;
  aR[3] = 0.0;
  aR[4] = 0.0;

  //second row  (dR2-{tauxn}/dQr)
  dR2_drhou = (c43*DxR[0]*avec[0] + DxR[1]*avec[1] + DxR[2]*avec[2]);
  dR2_drhov = (mc23*DxR[1]*avec[0] + DxR[0]*avec[1]);
  dR2_drhow = (mc23*DxR[2]*avec[0] + DxR[0]*avec[2]);
  aR[5] = -RKT*(dtauxn_drR);
  aR[6] = -RKT*dR2_drhou;
  aR[7] = -RKT*dR2_drhov;
  aR[8] = -RKT*dR2_drhow;
  aR[9] = 0.0;

  //third row   (dR3-{tauyn}/dQr)
  dR3_drhou = (mc23*DxR[0]*avec[1] + DxR[1]*avec[0]);
  dR3_drhov = (DxR[0]*avec[0] + c43*DxR[1]*avec[1] + DxR[2]*avec[2]);
  dR3_drhow = (mc23*DxR[2]*avec[1] + DxR[1]*avec[2]);
  aR[10] = -RKT*(dtauyn_drR);
  aR[11] = -RKT*dR3_drhou;
  aR[12] = -RKT*dR3_drhov;
  aR[13] = -RKT*dR3_drhow;
  aR[14] = 0.0;

  //fourth row  (dR4-{tauzn}/dQr)
  dR4_drhou = (mc23*DxR[0]*avec[2] + DxR[2]*avec[0]);
  dR4_drhov = (mc23*DxR[1]*avec[2] + DxR[2]*avec[1]);
  dR4_drhow = (DxR[0]*avec[0] + DxR[1]*avec[1] + c43*DxR[2]*avec[2]);
  aR[15] = -RKT*(dtauzn_drR);
  aR[16] = -RKT*dR4_drhou;
  aR[17] = -RKT*dR4_drhov;
  aR[18] = -RKT*dR4_drhow;
  aR[19] = 0.0;

  Type gamma = this->param->gamma;
  Type gm1 = (gamma - 1.0);

  Type v2L = (uL*uL + vL*vL + wL*wL);
  Type v2R = (uR*uR + vR*vR + wR*wR);
  
  Type dT_dPL = gamma/rhoL;
  Type dT_dPR = gamma/rhoR;

  //compute pressure derivatives
  Type dP_drL   = +gm1*0.5*v2L;
  Type dP_druL  = -gm1*uL;
  Type dP_drvL  = -gm1*vL;
  Type dP_drwL  = -gm1*wL;
  Type dP_dretL = +gm1;

  Type dP_drR   = +gm1*0.5*v2R;
  Type dP_druR  = -gm1*uR;
  Type dP_drvR  = -gm1*vR;
  Type dP_drwR  = -gm1*wR;
  Type dP_dretR = +gm1;

  //heat flux term
  Type cp = 1.0/(this->param->gamma - 1.0);
  //thermal conductivity
  Type k = mu/this->param->Pr*cp;
  //turbulent thermal conductivity
  Type kT = mut/this->param->PrT*cp;
  Type c1 = -(k + kT);

  //directional derivative part of temp grad
  Type TnL = (DxL[0]*avec[0] + DxL[1]*avec[1] + DxL[2]*avec[2])*dT_dPL*c1;
  Type TnR = (DxR[0]*avec[0] + DxR[1]*avec[1] + DxR[2]*avec[2])*dT_dPR*c1;


  //fifth row   (dR5/dQr)
  //NOTE: see Daniel's notes
  aR[20] = -RKT*(dtauxn_drR*u + dtauyn_drR*v + dtauzn_drR*w) + RK*TnR*(dP_drR - PR/rhoR);  
  aR[21] = -RKT*(dR2_drhou*u + dR3_drhou*v + dR4_drhou*w) + RK*TnR*dP_druR;
  aR[22] = -RKT*(dR2_drhov*u + dR3_drhov*v + dR4_drhov*w) + RK*TnR*dP_drvR;
  aR[23] = -RKT*(dR2_drhow*u + dR3_drhow*v + dR4_drhow*w) + RK*TnR*dP_drwR;
  aR[24] =  RK*TnR*dP_dretR;


  //Now, we do the same thing for the LHS jacobian
  //---------------------------------------------------------------------

  //first row   (dR1/dQl)
  aL[0] = 0.0;
  aL[1] = 0.0;
  aL[2] = 0.0;
  aL[3] = 0.0;
  aL[4] = 0.0;

  //second row  (dR2-{tauxn}/dQl)
  dR2_drhou = (c43*DxL[0]*avec[0] + DxL[1]*avec[1] + DxL[2]*avec[2]);
  dR2_drhov = (mc23*DxL[1]*avec[0] + DxL[0]*avec[1]);
  dR2_drhow = (mc23*DxL[2]*avec[0] + DxL[0]*avec[2]);
  aL[5] = -RKT*(dtauxn_drL);
  aL[6] = -RKT*dR2_drhou;
  aL[7] = -RKT*dR2_drhov;
  aL[8] = -RKT*dR2_drhow;
  aL[9] = 0.0;

  //third row   (dR3-{tauyn}/dQl)
  dR3_drhou = (mc23*DxL[0]*avec[1] + DxL[1]*avec[0]);
  dR3_drhov = (DxL[0]*avec[0] + c43*DxL[1]*avec[1] + DxL[2]*avec[2]);
  dR3_drhow = (mc23*DxL[2]*avec[1] + DxL[1]*avec[2]);
  aL[10] = -RKT*(dtauyn_drL);
  aL[11] = -RKT*dR3_drhou;
  aL[12] = -RKT*dR3_drhov;
  aL[13] = -RKT*dR3_drhow;
  aL[14] = 0.0;

  //fourth row  (dR4-{tauzn}/dQl)
  dR4_drhou = (mc23*DxL[0]*avec[2] + DxL[2]*avec[0]);
  dR4_drhov = (mc23*DxL[1]*avec[2] + DxL[2]*avec[1]);
  dR4_drhow = (DxL[0]*avec[0] + DxL[1]*avec[1] + c43*DxL[2]*avec[2]);
  aL[15] = -RKT*(dtauzn_drL);
  aL[16] = -RKT*dR4_drhou;
  aL[17] = -RKT*dR4_drhov;
  aL[18] = -RKT*dR4_drhow;
  aL[19] = 0.0;

  //fifth row   (dR5/dQr)
  //NOTE: see Daniel's notes
  aL[20] = -RKT*(dtauxn_drL*u + dtauyn_drL*v + dtauzn_drL*w) + RK*TnL*(dP_drL - PL/rhoL);
  aL[21] = -RKT*(dR2_drhou*u + dR3_drhou*v + dR4_drhou*w) + RK*TnL*dP_druL;
  aL[22] = -RKT*(dR2_drhov*u + dR3_drhov*v + dR4_drhov*w) + RK*TnL*dP_drvL;
  aL[23] = -RKT*(dR2_drhow*u + dR3_drhow*v + dR4_drhow*w) + RK*TnL*dP_drwL;
  aL[24] =  RK*TnL*dP_dretL;


  //okay, getting these signs correct is getting difficult.....
  //what is written above the the dR/dQ for a viscous edge with the flux
  //written on the RHS of the equations, since the flux is scattered 
  //positive for the right node and negative for the left node, we need
  //to flip the sign on the aL matrix so it returns with the correct sign
  for(i = 0; i < neqn*neqn; i++){
    aL[i] = -aL[i];
  }
  //pseudo-2D simulations
  if(this->param->symmetry2D > 0){
    if(this->param->symmetry2D == 1){
      for(i = 0; i < neqn; i++){
	aL[1*neqn + i] = aR[1*neqn + i] = 0.0;
      }
    }
    else if(this->param->symmetry2D == 2){
      for(i = 0; i < neqn; i++){
	aL[2*neqn + i] = aR[2*neqn + i] = 0.0;
      }
    }
    else if(this->param->symmetry2D == 3){
      for(i = 0; i < neqn; i++){
	aL[3*neqn + i] = aR[3*neqn + i] = 0.0;
      }
    }
  }
}
