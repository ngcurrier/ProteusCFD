#include "param.h"
#include "macros.h"
#include "matrix.h"
#include "solutionField.h"
#include <iostream>
#include <cmath>


//Variables native in Q[] 
//p
//u
//v
//w

//IMPORTANT NOTE:  BETA here is the artificial compressibility parameter, however, 
//                 all computations are performed with 1.0/beta so that a reasonable
//                 choice for beta (~Ma^2) for variable mach/preconditioned regimes
//                 is also a reasonable choice here.  Convenience only to keep a bunch
//                 of switch statements out of the code above this level.

template <class Type>
IncompressibleEqnSet<Type>::IncompressibleEqnSet(SolutionSpace<Type>* space, Param<Type>* p)
{
  this->neqn = 4;
  this->nauxvars = 0;
  this->space = space;
  this->param = p;

  //variable set is conservative
  this->varsConservative = 1;

  this->Qinf = new Type[this->neqn + this->nauxvars];

  //set types of variables for output i.e. scalar, vector, etc.
  this->idata = new DataInfo(this->neqn+this->nauxvars, std::string("variableQ"));
  this->idata->AddScalar(0, std::string("Pressure"));
  this->idata->AddVector(1, std::string("Velocity"));
  this->idata->Verify();
  
  //set gradients required
  this->gdata = new DataInfo((this->neqn)*3, "gradVariableQ");
  this->gdata->AddVector(0*3, "Grad-Pressure");
  this->gdata->AddVector(1*3, "Grad-u");
  this->gdata->AddVector(2*3, "Grad-v");
  this->gdata->AddVector(3*3, "Grad-w");
  this->gdata->Verify();

}

template <class Type>
IncompressibleEqnSet<Type>::~IncompressibleEqnSet()
{
  delete [] this->Qinf;
  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::InitEqnSet()
{
  //calculate Qinf values for class variables
  this->UpdateQinf();
  
  //allocate solution memory 
  this->space->AddField(*this->idata, FIELDS::STATE_TIME, FIELDS::VAR_EVERYWHERE);
  SolutionField<Type> & field = this->space->GetField("variableQ");
  this->space->q = field.GetData(FIELDS::STATE_NP1);
  this->space->qold = field.GetData(FIELDS::STATE_N);
  this->space->qoldm1 = field.GetData(FIELDS::STATE_NM1);

  //allocate solution memory for the gradients we need
  if(this->param->viscous || (this->param->sorder > 1)){
    this->space->AddField(*this->gdata, FIELDS::STATE_NONE, FIELDS::VAR_INTERIOR);
    SolutionField<Type> & gfield  = this->space->GetField("gradVariableQ");
    this->space->qgrad = gfield.GetData(FIELDS::STATE_NONE);
  }

  //maybe this should be different? Use ideal gas law or something to set freestream density?  Unsure.
  Type rho = 1.0;
  Type V = this->param->GetVelocity(this->space->iter);
  Type Mach = V;
  this->param->Re = (rho*this->param->ref_density)*(this->param->ref_velocity*Mach)*
    this->param->ref_length/this->param->ref_viscosity;
}

template <class Type>
void IncompressibleEqnSet<Type>::RoeFlux(Type QL[], Type QR[], Type avec[], Type vdotn, 
					 Type gamma, Type flux[], Type beta)
{
  Type ibeta = 1.0/beta;
  Type PL,uL,vL,wL;
  Type PR,uR,vR,wR;
  Type u,v,w;
  Type dp,du,dv,dw;
  Type theta,thetamat,thetaL,thetaR;
  Type eig3,eig4,eig3ob,eig4ob;
  Type alpha;
  Type tmp;
  Type c,cplus,cminus;
  Type psi,tmp2,tmp3;
  Type dudot;

  //magnitude stored in 4th entry in edge structures area
  Type area = avec[3];

  //extract the left and right states
  PL = QL[0];
  uL = QL[1];
  vL = QL[2];
  wL = QL[3];

  PR = QR[0];
  uR = QR[1];
  vR = QR[2];
  wR = QR[3];

  //compute the jumps in the dependent variables
  dp = PR - PL;
  du = uR - uL;
  dv = vR - vL;
  dw = wR - wL;

  //TODO: Call RoeVariables here instead of the mess that follows
  //

  //compute thetaL and thetaR
  thetaL = avec[0]*uL + avec[1]*vL + avec[2]*wL + vdotn;
  thetaR = avec[0]*uR + avec[1]*vR + avec[2]*wR + vdotn;

  //compute the Roe averaged variables
  u = 0.5*(uL + uR);
  v = 0.5*(vL + vR);
  w = 0.5*(wL + wR);
  theta   = 0.5*(thetaL + thetaR);
  thetamat  = theta - vdotn;

  //compute the pseudo speed of sound
  tmp = theta - 0.5*vdotn;
  c   = sqrt(tmp*tmp + ibeta);
  cplus  = 0.5*vdotn + c;
  cminus  = 0.5*vdotn - c;
 
  psi = 4.0*c*(1.0 + theta*thetamat/ibeta);

  dudot = du*avec[0] + dv*avec[1] + dw*avec[2];

  //flow from right to left
  if(real(theta) < 0.0)
    {
      eig3 = theta - cminus;

      eig3ob = eig3/ibeta;

      flux[0] = ibeta*(thetaR - vdotn);
      flux[1] = uR*thetaR + PR*avec[0];
      flux[2] = vR*thetaR + PR*avec[1];
      flux[3] = wR*thetaR + PR*avec[2];

      tmp3 = 0.25*psi/c - cplus/ibeta*thetamat;
      
      alpha = eig3*2.0/psi*(tmp3*dp + cplus*dudot);
      
      flux[0] -= -cminus*alpha;
      flux[1] -= (avec[0] + u*eig3ob)*alpha;
      flux[2] -= (avec[1] + v*eig3ob)*alpha;
      flux[3] -= (avec[2] + w*eig3ob)*alpha;
    }

  //flow from left to right
  else
    {
      eig4 = theta - cplus;
      eig4ob = eig4/ibeta;

      flux[0] = ibeta*(thetaL - vdotn);
      flux[1] = uL*thetaL + PL*avec[0];
      flux[2] = vL*thetaL + PL*avec[1];
      flux[3] = wL*thetaL + PL*avec[2];

      tmp2 = -0.25*psi/c + cminus/ibeta*thetamat;
      
      alpha = eig4*2.0/psi*(tmp2*dp - cminus*dudot);
      
      flux[0] += -cplus*alpha;
      flux[1] += (avec[0] + u*eig4ob)*alpha;
      flux[2] += (avec[1] + v*eig4ob)*alpha;
      flux[3] += (avec[2] + w*eig4ob)*alpha;
    }

  flux[0] *= area;
  flux[1] *= area;
  flux[2] *= area;
  flux[3] *= area;

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
Bool IncompressibleEqnSet<Type>::RoeVariables(Type QL[], Type QR[], Type gamma, Type Qroe[])
{
  Int i;

  Int neqn = this->neqn;
  for(i = 0; i < neqn; i++){
    Qroe[i] = 0.5*(QL[i] + QR[i]);
  }

  return (true);
}

template <class Type>
void IncompressibleEqnSet<Type>::Eigensystem(Type* Q, Type* avec, Type vdotn, Type gamma, 
					     Type* eigenvalues, Type* T, Type* Tinv, Type beta)
{
  Type ibeta = 1.0/beta;
  Type nx = avec[0];
  Type ny = avec[1];
  Type nz = avec[2];
  Type u = Q[1];
  Type v = Q[2];
  Type w = Q[3];
  Type thetamat = u*nx + v*ny + w*nz;
  Type theta = thetamat + vdotn;
  Type tmp = theta - 0.5*vdotn;
  Type c2 = tmp*tmp + ibeta;
  Type c = sqrt(c2);
  Type cplus = 0.5*vdotn + c;
  Type cminus = 0.5*vdotn - c;
  Type psi = 4.0*c*(1.0 + theta*thetamat/ibeta);
  Type dpsi = 2.0/psi;

  Type v1[3];
  Type v2[3];
  Type psi1,psi2,psi3;
  Type psi4,psi5,psi6;
  Type psi7,psi8,psi9;

  //avec = v1 x v2 ...
  PerpVectors(avec, v1, v2);

  psi1 = nx + u*theta/ibeta;
  psi2 = ny + v*theta/ibeta;
  psi3 = nz + w*theta/ibeta;
  psi4 = v2[1]*psi3 - v2[2]*psi2;
  psi5 = v2[0]*psi3 - v2[2]*psi1;
  psi6 = v2[0]*psi2 - v2[1]*psi1;
  psi7 = v1[1]*psi3 - v1[2]*psi2;
  psi8 = v1[0]*psi3 - v1[2]*psi1;
  psi9 = v1[0]*psi2 - v1[1]*psi1;

  //eigenvalues
  eigenvalues[0] = theta;
  eigenvalues[1] = theta;
  eigenvalues[2] = theta - cminus;
  eigenvalues[3] = theta - cplus;  
  
  Type eig3 = eigenvalues[2];
  Type eig4 = eigenvalues[3];

  //Right eigenvectors
  T[0] = 0.0;
  T[1] = 0.0;
  T[2] = -cminus;
  T[3] = -cplus;

  T[4] = 2.0*v1[0];
  T[5] = 2.0*v2[0];
  T[6] = nx + u*eig3/ibeta;
  T[7] = nx + u*eig4/ibeta;

  T[8] = 2.0*v1[1];
  T[9] = 2.0*v2[1];
  T[10] = ny + v*eig3/ibeta;
  T[11] = ny + v*eig4/ibeta;

  T[12] = 2.0*v1[2];
  T[13] = 2.0*v2[2];
  T[14] = nz + w*eig3/ibeta;
  T[15] = nz + w*eig4/ibeta;

  //Left eigenvectors
  Tinv[0] = dpsi*(c/ibeta*(-u*psi4 + v*psi5 - w*psi6));
  Tinv[1] = dpsi*(c*psi4);
  Tinv[2] = dpsi*(-c*psi5);
  Tinv[3] = dpsi*(c*psi6);

  Tinv[4] = dpsi*(c/ibeta*(u*psi7 - v*psi8 + w*psi9));
  Tinv[5] = dpsi*(-c*psi7);
  Tinv[6] = dpsi*(c*psi8);
  Tinv[7] = dpsi*(-c*psi9);
  
  Tinv[8] = dpsi*(0.25*psi/c - cplus/ibeta*thetamat);
  Tinv[9] = dpsi*cplus*nx;
  Tinv[10] = dpsi*cplus*ny;
  Tinv[11] = dpsi*cplus*nz;

  Tinv[12] = dpsi*(-0.25*psi/c + cminus/ibeta*thetamat);
  Tinv[13] = dpsi*(-cminus*nx);
  Tinv[14] = dpsi*(-cminus*ny);
  Tinv[15] = dpsi*(-cminus*nz);

#if 0
  //check that T and Tinv really are inverses of each other
  Int yes;
  yes = MatInvCheck(T, Tinv, this->neqn);
  if(!yes){
    std::cerr << "Eigensystem not valid -- inverse check failed!!" << std::endl;
  } 
#endif

  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::Flux(Type* Q, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta)
{
  Type ibeta = 1.0/beta;
  Type thetamat = avec[0]*Q[1] + avec[1]*Q[2] + avec[2]*Q[3];
  Type theta  = thetamat + vdotn;
  Type area = avec[3];

  flux[0] = area*ibeta*thetamat;
  flux[1] = area*(Q[1]*theta + Q[0]*avec[0]);
  flux[2] = area*(Q[2]*theta + Q[0]*avec[1]);
  flux[3] = area*(Q[3]*theta + Q[0]*avec[2]);

  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::ViscousFlux(Type* Q, Type* grad, Type* avec, Type mut, Type* flux)
{
  //gradient is passed in with velocity components and parsed for readability
  Type ux, uy, uz, vx, vy, vz, wx, wy, wz;  
  ux = grad[3];
  uy = grad[4];
  uz = grad[5];
  vx = grad[6];
  vy = grad[7];
  vz = grad[8];
  wx = grad[9];
  wy = grad[10];
  wz = grad[11];

  Type c1 = avec[3]/this->param->Re;

  Type mu = ComputeViscosity(Q);
  Type tmut = (mu + mut);

  flux[0] = 0.0;
  flux[1] = -c1*(tmut*(2.0*ux*avec[0] + (uy + vx)*avec[1] + (uz + wx)*avec[2]));
  flux[2] = -c1*(tmut*((vx + uy)*avec[0] + 2.0*vy*avec[1] + (vz + wy)*avec[2]));
  flux[3] = -c1*(tmut*((wx + uz)*avec[0] + (wy + vz)*avec[1] + 2.0*wz*avec[2]));

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
Type IncompressibleEqnSet<Type>::MaxEigenvalue(Type* Q, Type* avec, Type vdotn, Type gamma, Type beta)
{
  Type ibeta = 1.0/beta;
  Type maxeig;
  Type U;
  Type tmp;
  Type eig3,eig4;
  Type c;

  U = GetTheta(Q, avec, vdotn);
  
  tmp = U - 0.5*vdotn;

  c = sqrt(tmp*tmp + ibeta);

  eig3 = U + c - 0.5*vdotn;
  eig4 = U - c - 0.5*vdotn;
  
  maxeig = MAX(CAbs(eig3),CAbs(eig4));
  
  return (maxeig);
}

template <class Type>
void IncompressibleEqnSet<Type>::UpdateQinf()
{
  //use GetVelocity() function for ramping
  //TODO: port this to use a uinf/uref instead of Mach
  //this is the correct non-dimensionlization for all regimes involving flow
  //for compressible it is uinf/(Mach_ref/C_ref).. and so forth... oh well
  Type V = this->param->GetVelocity(this->space->iter);
  Type Mach = V;
  
  Type u = this->param->flowdir[0]*Mach;
  Type v = this->param->flowdir[1]*Mach;
  Type w = this->param->flowdir[2]*Mach;

  //set class variables to default initialization for Qinf
  this->Qinf[0] = 0.0; //p
  this->Qinf[1] = u;
  this->Qinf[2] = v;
  this->Qinf[3] = w;

  this->ComputeAuxiliaryVariables(this->Qinf);

  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::SetInitialConditions()
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
  
  //set all the nodes interior and phantom
  for(i = 0; i < (nnode+gnode+nbnode); i++){
    for(j = 0; j < neqn; j++){
      this->space->q[i*(neqn+nauxvars) + j] = this->Qinf[j];
    }
    this->ComputeAuxiliaryVariables(&this->space->q[i*(neqn+nauxvars)]);
  }
  
  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::ApplyDQ(Type* dQ, Type* Q, Type* xyz)
{
  //we shouldn't clip pressure since we assume gauge pressure and negative
  //pressures are certainly okay there
  Q[0] += dQ[0];
  //add up the velocities
  Q[1] += dQ[1];
  Q[2] += dQ[2];
  Q[3] += dQ[3];

  this->ComputeAuxiliaryVariables(Q);
}

template <class Type>
Type IncompressibleEqnSet<Type>::GetTheta(Type* Q, Type* avec, Type vdotn)
{
  Type nx = avec[0];
  Type ny = avec[1];
  Type nz = avec[2];
  Type u = Q[1];
  Type v = Q[2];
  Type w = Q[3];

  return (u*nx + v*ny + w*nz) + vdotn;
}

template <class Type>
Type IncompressibleEqnSet<Type>::GetDensity(Type* Q)
{
  return (1.0);
}

template <class Type>
Type IncompressibleEqnSet<Type>::GetPressure(Type* Q)
{
  return (Q[0]);
}

template <class Type>
Type IncompressibleEqnSet<Type>::GetTemperature(Type* Q)
{
  return (1.0);
}

template <class Type>
Type IncompressibleEqnSet<Type>::GetGamma(Type* Q)
{
  return (this->param->gamma);
}

template <class Type>
Type IncompressibleEqnSet<Type>::GetCp(Type* Q, Type gamma)
{
  return (2.0*GetPressure(Q));
}

template <class Type>
Type IncompressibleEqnSet<Type>::GetCf(Type tauw, Type rho)
{
  return (2.0*tauw);
}

template <class Type>
Type IncompressibleEqnSet<Type>::ComputeViscosity(Type* Q)
{
  return (1.0);
}

template <class Type>
Int IncompressibleEqnSet<Type>::GetVelocityGradLocation()
{
  return 1;
}

template <class Type>
Int IncompressibleEqnSet<Type>::GetMomentumLocation()
{
  return 1;
}

template <class Type>
Int IncompressibleEqnSet<Type>::GetGradientsLocation(std::vector<Int>& gradientLoc)
{
  //pressure
  gradientLoc.push_back(0);
  //velocity
  gradientLoc.push_back(1);
  gradientLoc.push_back(2);
  gradientLoc.push_back(3);
  return gradientLoc.size();
}

template <class Type>
void IncompressibleEqnSet<Type>::NativeToExtrapolated(Type* Q)
{ }

template <class Type>
void IncompressibleEqnSet<Type>::ExtrapolatedToNative(Type* Q)
{ }

template <class Type>
Int IncompressibleEqnSet<Type>::BadExtrapolation(Type* Q)
{
  //nothing to check since we use relative pressure in this eqnset
  return false;
}

template <class Type>
void IncompressibleEqnSet<Type>::ExtrapolateVariables(Type* Qho, const Type* Q, const Type* dQedge, 
						      const Type* gradQ, const Type* dx, const Type* limiter)
{
  //pressure
  Qho[0] = Q[0] + this->ExtrapolateCorrection(dQedge[0], &gradQ[0*3], dx)*limiter[0];
  //velocity
  Qho[1] = Q[1] + this->ExtrapolateCorrection(dQedge[1], &gradQ[1*3], dx)*limiter[1];
  Qho[2] = Q[2] + this->ExtrapolateCorrection(dQedge[2], &gradQ[2*3], dx)*limiter[2];
  Qho[3] = Q[3] + this->ExtrapolateCorrection(dQedge[3], &gradQ[3*3], dx)*limiter[3];
}



template <class Type>
Type IncompressibleEqnSet<Type>::ComputePressure(Type* Q, Type gamma)
{
  return (GetPressure(Q));
}

template <class Type>
void IncompressibleEqnSet<Type>::ComputeStressVector(Type* grad, Type* avec, Type mu, Type* stress)
{

  //gradient is passed in with velocity components and parsed for readability
  Type ux, uy, uz, vx, vy, vz, wx, wy, wz;  
  ux = grad[3];
  uy = grad[4];
  uz = grad[5];
  vx = grad[6];
  vy = grad[7];
  vz = grad[8];
  wx = grad[9];
  wy = grad[10];
  wz = grad[11];

  Type c1 = 1.0/this->param->Re;

  //negative sign indicates stress on the body, not the flow
  stress[0] = -c1*(mu*(2.0*ux*avec[0] + (uy + vx)*avec[1] + (uz + wx)*avec[2]));
  stress[1] = -c1*(mu*((vx + uy)*avec[0] + 2.0*vy*avec[1] + (vz + wy)*avec[2]));
  stress[2] = -c1*(mu*((wx + uz)*avec[0] + (wy + vz)*avec[1] + 2.0*wz*avec[2]));

  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::GetFarfieldBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta)
{
  Int i,j;
  Type gamma = this->param->gamma;
  Int neqn = this->neqn;
  Int nauxvars = this->nauxvars;
  Int nvars = neqn+nauxvars;
  Type* avg = (Type*)alloca(sizeof(Type)*this->neqn);
  Type* eigenvals = (Type*)alloca(sizeof(Type)*neqn);
  Type* vinv = (Type*)alloca(sizeof(Type)*neqn*neqn);
  Type* v = (Type*)alloca(sizeof(Type)*neqn*neqn);
  Type* rhs = (Type*)alloca(sizeof(Type)*neqn);
  Type* q = (Type*)alloca(sizeof(Type)*neqn);
  Type* qinf = (Type*)alloca(sizeof(Type)*neqn);
  
  Int SetFarfield_Pb = 0;
  
  //subiteration for design derivatives
  for(Int subit = 0; subit < 10; subit++){
    
    //average variables to compute eigensystem
    for(i = 0; i < neqn; i++){
      avg[i] = (QL[i] + QR[i])/2.0;
    }
    
    Eigensystem(avg, avec, vdotn, gamma, eigenvals, v, vinv, beta);
    
    if(this->param->no_cvbc){
      //outflow
      if(real(eigenvals[0]) >= 0.0){
	memcpy(QR,QL,nvars*sizeof(Type));
      }
      //inflow
      else{
	memcpy(QR,Qinf,nvars*sizeof(Type));
      }
    }
    else{      
            
      //silence compiler warnings
      for(i = 0; i < neqn; i++){
	rhs[i] = q[i] = qinf[i] = 0.0;
      }
      
      //get the primitive variables on the interior
      memcpy(q,QL,neqn*sizeof(Type));
      memcpy(qinf,Qinf,neqn*sizeof(Type));
      
      //Compute the rhs of the characteristic equations as done in Daniel's 
      //dissertation, we pick the correct Q based on the sign of the eigenvalues
      for(i = 0; i < neqn; i++){
	rhs[i] = 0.0;
	for(j = 0; j < neqn; j++){
	  rhs[i] += vinv[i*neqn + j]*(real(eigenvals[i]) > 0.0 ? q[j] : qinf[j]);
	}
      }
      
      //This can be used to explicitly set the back pressure on the domain
      if(SetFarfield_Pb && (real(eigenvals[0]) > 0.0)){
	rhs[3] = qinf[0];
	vinv[(3)*neqn + 0] = 1.0;
	vinv[(3)*neqn + 1] = 0.0;
	vinv[(3)*neqn + 2] = 0.0;
	vinv[(3)*neqn + 3] = 0.0;
      }
      
      //solve the equations
      Int p[neqn];
      LU(vinv, p, neqn);
      LuSolve(vinv, rhs, p, QR, neqn);
      //copy solution from rhs to QR
      memcpy(QR, rhs, neqn*sizeof(Type));
    }
  }
    
  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::GetInviscidWallBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta)
{
  Type gamma = this->param->gamma;
  Int neqn = this->neqn;
  Type* avg = (Type*)alloca(sizeof(Type)*this->neqn);
  Type* eigenvals = (Type*)alloca(sizeof(Type)*neqn);
  Type* Tinv = (Type*)alloca(sizeof(Type)*neqn*neqn);
  Type* T = (Type*)alloca(sizeof(Type)*neqn*neqn);
  Type* rhs = (Type*)alloca(sizeof(Type)*neqn);
  Type* q = (Type*)alloca(sizeof(Type)*neqn);
  Type* qinf = (Type*)alloca(sizeof(Type)*neqn);

  //extract grid speeds in normal direction
  //Type theta = GetTheta(QL, avec, vdotn);
  //vdotn = -dot(gridspeeds, avec)
  Type u, v, w;
  u = vdotn*avec[0];
  v = vdotn*avec[1];
  w = vdotn*avec[2];
  
  Type effVel[3];
  effVel[0] = QL[1] + u;
  effVel[1] = QL[2] + v;
  effVel[2] = QL[3] + w;

  if(!this->param->no_cvbc && false){
    //subiteration for design derivatives
    for(Int subit = 0; subit < 10; subit++){
      
      //average variables to compute eigensystem
      for(Int i = 0; i < neqn; i++){
	avg[i] = 0.5*(QL[i] + QR[i]);
      }
      
      Eigensystem(avg, avec, vdotn, gamma, eigenvals, T, Tinv, beta);
      
      //silence compiler warnings
      for(Int i = 0; i < neqn; i++){
	rhs[i] = q[i] = qinf[i] = 0.0;
      }
      
      //get the primitive variables on the interior
      memcpy(q,QL,neqn*sizeof(Type));
      
      //Compute the rhs of the characteristic equations as done in Daniel's dissertation
      for(Int i = 0; i < neqn; i++){
	rhs[i] = 0.0;
	for(Int j = 0; j < neqn; j++){
	  rhs[i] += Tinv[i*neqn + j]*q[j];
	}
      }
      
      //hard code the first equation to enforce theta = 0
      Tinv[0*neqn + 1] = avec[0];
      Tinv[0*neqn + 2] = avec[1];
      Tinv[0*neqn + 3] = avec[2];
      Tinv[0*neqn + 0] = 0.0;
      rhs[0] = vdotn;

      //solve the equations
      Int p[neqn];
      LU(Tinv, p, neqn);
      LuSolve(Tinv, rhs, p, QR, neqn);
      //copy solution from rhs to QR
      memcpy(QR, rhs, neqn*sizeof(Type));
    }
  }
  else{
    //set QR
    QR[0] = QL[0];
    MirrorVector(effVel, avec, &QR[1]);
  }
  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::GetInternalInflowBoundaryVariables(Type* QL, Type* QR, Type* Qinf, 
								    Type pressure, Type* densities,
								    Type* flowDirection, Type velocity)
{
  //pressure that is passed in is dynamic pressure
  //we must calculate the static pressure
  Type u = flowDirection[0]*velocity;
  Type v = flowDirection[1]*velocity;
  Type w = flowDirection[2]*velocity;

  Type Pi = QL[0];

  //Type Ps = pressure - 0.5*(u*u + v*v + w*w);
  //QR[0] = Ps;
  QR[0] = Pi;
  QR[1] = u;
  QR[2] = v;
  QR[3] = w;

#if 0
  //Here we are allowing the backpressure to float to the internal pressure on the boundary
  //this is a subsonic imposition of an internal inflow back pressure
  QR[0] = QL[0];
  QR[1] = Qinf[1];
  QR[2] = Qinf[2];
  QR[3] = Qinf[3];
#endif

  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::GetInternalOutflowBoundaryVariables(Type* QL, Type* QR, Type pressure, Type gamma)
{
  //pressure that is passed in is dynamic pressure
  //we must calculate the static pressure
  Type u = QL[1];
  Type v = QL[2];
  Type w = QL[3];
  Type Ps = pressure - 0.5*(u*u + v*v + w*w);

  QR[0] = Ps;
  QR[1] = u;
  QR[2] = v;
  QR[3] = w;
  return;
}

template <class Type> 
void IncompressibleEqnSet<Type>::GetViscousWallBoundaryVariables(Type* QL, Type* QR, Type* vel, Int normalNode, Type* normalQ, Type Twall)
{
  //set no slip on the surface
  QR[1] = QL[1] = vel[0];
  QR[2] = QL[2] = vel[1];
  QR[3] = QL[3] = vel[2];

  //we need to enforce the zero pressure gradient condition with the normal off wall node
  QR[0] = QL[0];

  return;
}


template <class Type>
void IncompressibleEqnSet<Type>::ModifyViscousWallJacobian(Type* QL, Type* QR, Type* vel, Int cvid, 
							   CRS<Type>* crs, Int normalNode, Type Twall)
{
  //blank all the subrows for u, v, w components
  crs->A->BlankSubRow(cvid, 1);
  crs->A->BlankSubRow(cvid, 2);
  crs->A->BlankSubRow(cvid, 3);

  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::ModifyViscousWallResidual(Type* res, Type* vel, Int normalNode, Type Twall)
{
  //zero the rows we are setting explicitly, namely u,v,w
  //this has the effect of not allowing du, dv, dw to change during the
  //linear system solve
  res[1] = 0.0;
  res[2] = 0.0;
  res[3] = 0.0;
  return;
}

template <class Type>
void IncompressibleEqnSet<Type>::ViscousJacobian(Type* QL, Type* QR, Type* dx, Type s2, 
						 Type* avec, Type mut, Type* aL, Type* aR)
{
  Type dxnx = dx[0]*avec[0];
  Type dyny = dx[1]*avec[1];
  Type dznz = dx[2]*avec[2];

  Type c1 = avec[3]/(this->param->Re*s2);

  //don't do averaging b/c this returns 1.0 anyway
  Type mu = ComputeViscosity(QL);
  Type tmut = (mu + mut);
  
  Int neqn = this->neqn;

  //first row   (p)
  aR[0*neqn + 0] = 0.0;
  aR[0*neqn + 1] = 0.0;
  aR[0*neqn + 2] = 0.0;
  aR[0*neqn + 3] = 0.0;

  //second row  (u)
  aR[1*neqn + 0] = 0.0;
  aR[1*neqn + 1] = -c1*tmut*(2.0*dxnx + dyny + dznz);
  aR[1*neqn + 2] = -c1*tmut*dx[0]*avec[1];
  aR[1*neqn + 3] = -c1*tmut*dx[0]*avec[2];

  //third row   (v)
  aR[2*neqn + 0] = 0.0;
  aR[2*neqn + 1] = -c1*tmut*dx[1]*avec[0];
  aR[2*neqn + 2] = -c1*tmut*(dxnx + 2.0*dyny + dznz);
  aR[2*neqn + 3] = -c1*tmut*dx[1]*avec[2];

  //fourth row  (w)
  aR[3*neqn + 0] = 0.0;
  aR[3*neqn + 1] = -c1*tmut*dx[2]*avec[0];
  aR[3*neqn + 2] = -c1*tmut*dx[2]*avec[1];
  aR[3*neqn + 3] = -c1*tmut*(dxnx + dyny + 2.0*dznz);

  //flip signs for the left jacobian
  for(Int i = 0; i < neqn*neqn; i++){
    aL[i] = aR[i];
  }

  return;
}
