#include "dataInfo.h"
#include "mesh.h"
#include "parallel.h"
#include "param.h"
#include "gradient.h"

template <class Type>
EqnSet<Type>::EqnSet()
{
  return;
}

template <class Type>
EqnSet<Type>::~EqnSet()
{
  delete idata;
  delete gdata;

  return;
}

template <class Type>
bool EqnSet<Type>::BoundaryFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type* flux, Type beta){
  bool error = false;
#if 0
  //use algebraic flux only
  this->Flux(QR, avec, vdotn, this->param->gamma, flux, beta);
  Type* flux2 = (Type*)alloca(sizeof(Type)*this->neqn);
  this->Flux(QL, avec, vdotn, this->param->gamma, flux2, beta);
  for(int i = 0; i < this->neqn; i++){
    flux[i] = 0.5*avec[3]*(flux[i] + flux2[i]);
    //flux[i] = avec[3]*flux[i];
  }
#else
  this->NumericalFlux(QL, QR, avec, vdotn, flux, beta);
#endif
  for(int i = 0; i < this->neqn; ++i){
    if(std::isnan(real(flux[i]))){
      std::cerr << "WARNING: EqnSet::BoundaryFlux() kneecapping a nan flux term - flux[" << i << "]" << std::endl;
      for(int j = 0; j < this->neqn; ++j){
	std::cerr << "flux[" << j << "] = " << flux[j] << std::endl;
      }      
      for(int j = 0; j < this->neqn; ++j){
	std::cerr << "QL[" << j << "] = " << QL[j] << std::endl;
      }
      for(int j = 0; j < this->neqn; ++j){
	std::cerr << "QR[" << j << "] = " << QR[j] << std::endl;
      }
      flux[i] = 0.0;
      error = true;
    }
  }
  return error;
}

template <class Type>
bool EqnSet<Type>::NumericalFlux(Type QL[], Type QR[], Type avec[], Type vdotn, Type flux[], Type beta)
{
  bool error = false;
  Int type = this->param->flux_id;
  if(type == AlgebraicFluxType){
    Type* flux2 = (Type*)alloca(sizeof(Type)*this->neqn);
    this->Flux(QL, avec, vdotn, this->param->gamma, flux, beta);
    this->Flux(QR, avec, vdotn, this->param->gamma, flux2, beta);
    for(int i = 0; i < this->neqn; i++){
      flux[i] = 0.5*avec[3]*(flux[i] + flux2[i]);
    }
  }
  else if(type == RoeFluxType){
    this->RoeFlux(QL, QR, avec, vdotn, this->param->gamma, flux, beta);
  }
  else if(type == HLLCFluxType){
    this->HLLCFlux(QL, QR, avec, vdotn, this->param->gamma, flux, beta);
  }
  for(int i = 0; i < this->neqn; ++i){
    if(std::isnan(real(flux[i]))){
      std::cerr << "WARNING: EqnSet::NumericalFlux() kneecapping a nan flux term - flux[" << i << "]" << std::endl;
      for(int j = 0; j < this->neqn; ++j){
	std::cerr << "flux[" << j << "] = " << flux[j] << std::endl;
      }      
      for(int j = 0; j < this->neqn; ++j){
	std::cerr << "QL[" << j << "] = " << QL[j] << std::endl;
      }
      for(int j = 0; j < this->neqn; ++j){
	std::cerr << "QR[" << j << "] = " << QR[j] << std::endl;
      }
      flux[i] = 0.0;
      error = true;
    }
  }
  return error;
}

template <class Type>
void EqnSet<Type>::ViscousJacobian(Type* qL, Type* qR, Type* dx, Type s2, Type* avec, 
				   Type mut, Type* aL, Type* aR)
{
  Int i, j, k;

  //perturbation
  Type h = 1.0e-8;
  
  std::vector<Int> gradLoc;
  Int nterms = GetGradientsLocation(gradLoc);

  Int nvars = neqn+nauxvars;
  Type* grad = (Type*)alloca(sizeof(Type)*nterms*3);
  Type* flux = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxL = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxR = (Type*)alloca(sizeof(Type)*neqn);
  Type* Qavg = (Type*)alloca(sizeof(Type)*nvars);
  Type* QPL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QPR = (Type*)alloca(sizeof(Type)*nvars);
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars);

  for(i = 0; i < neqn; i++){
    Qavg[i] = 0.5*(QL[i] + QR[i]);
  }
  ComputeAuxiliaryVariables(Qavg);

  memcpy(QL, qL, sizeof(Type)*nvars);
  memcpy(QR, qR, sizeof(Type)*nvars);

  DirectionalDerivativeGrad(dx, s2, avec, QL, QR, grad, gradLoc);
  ViscousFlux(Qavg, grad, avec, mut, flux);

  //this is the numerical derivative implementation of viscous jacobian
  //used for finite-rate chemistry, etc.
  for(i = 0; i < neqn; i++){
    memcpy(QPL, qL, sizeof(Type)*nvars); 
    memcpy(QPR, qR, sizeof(Type)*nvars); 

    //perturb variables
    QPL[i] += h;
    QPR[i] += h;
    ComputeAuxiliaryVariables(QPL);
    ComputeAuxiliaryVariables(QPR);

    //not sure if this is correct but since we perturb both sides by the same h
    //we don't need to compute a new qavg for each call to viscous flux
    for(k = 0; k < neqn; k++){
      Qavg[k] = 0.5*(QPL[k] + QR[k]);
    }
    ComputeAuxiliaryVariables(Qavg);

    DirectionalDerivativeGrad(dx, s2, avec, QPL, QR, grad, gradLoc);
    ViscousFlux(Qavg, grad, avec, mut, fluxL);

    DirectionalDerivativeGrad(dx, s2, avec, QL, QPR, grad, gradLoc);
    ViscousFlux(Qavg, grad, avec, mut, fluxR);
    
    for(j = 0; j < neqn; j++){
      aL[j*neqn + i] = (flux[j] - fluxL[j])/h;
      aR[j*neqn + i] = (fluxR[j] - flux[j])/h;
    }

  }

  return;
}


template <class Type>
void EqnSet<Type>::SourceTermJacobian(Type* Q, Type vol, Type* A)
{
  //this is the numerical derivative way
  
  Int i, j;
  Int neqn2 = neqn*neqn;
  Type* source = (Type*)alloca(sizeof(Type)*neqn);
  Type* sourceP = (Type*)alloca(sizeof(Type)*neqn);
  Type h = 1.0e-8;
  Type* QP = (Type*)alloca(sizeof(Type)*(neqn+nauxvars));
  
  MemBlank(A, neqn2);
  
  this->SourceTerm(Q, vol, source); 
  
  for(i = 0; i < neqn; i++){
    memcpy(QP, Q, sizeof(Type)*neqn);
    QP[i] += h;
    ComputeAuxiliaryVariables(QP);
    SourceTerm(QP, vol, sourceP);
    for(j = 0; j < neqn; j++){
      A[j*neqn + i] = (sourceP[j] - source[j])/h;      
    }
  }
}

//default is to just add v/dt terms to the diagonal
//some eqnsets will have dQ/dq terms if they are in non-conservative variables
//this is cnp1*vol/dt*M where M is dq/dQ
//dt - real time step used for residuals
//dtau - pseudotimestep used for convergence enhancement only
template <class Type>
void EqnSet<Type>::ContributeTemporalTerms(Type* Q, Type vol, Type cnp1, 
					   Type dt, Type dtau, Type* A, Type beta)
{
  for(Int i = 0; i < neqn; i++){
    if(param->useLocalTimeStepping && (real(param->dt) > 0.0)){
      //unsteady solution - add dt term
      //the second part of this term comes from pseudo timestepping
      A[i*neqn + i] += cnp1*vol/dt + vol/dtau;
    }
    else{
      A[i*neqn + i] += cnp1*vol/dtau;
    }
  }
}

template <class Type>
void EqnSet<Type>::GetPythonBoundaryVariables(Type* QL, Type* QR, Type* wallx, Type* wallAvec)
{
#ifdef _HAS_PYTHON
  PythonWrapper pywrap("./", "getBoundaryVariables", "getBoundaryVariables");
  pywrap.GetBoundaryVariables(QL, QR, this->neqn, this->nauxvars, wallx, wallAvec);
#else
  Abort << "EqnSet::GetPythonBoundaryVariables() - python not built with solver";
#endif
}
