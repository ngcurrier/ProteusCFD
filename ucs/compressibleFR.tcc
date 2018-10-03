#include "macros.h"
#include "eqnset.h"
#include "param.h"
#include "matrix.h"
#include "chem.h"
#include "solutionSpace.h"
#include "solutionField.h"
#include "pythonInterface.h"
#include <cmath>
#include <iostream>

//Variable locations
//rho_1 (kg/m^3)
//rho_2
//...
//rho_NScompute temperature from energy
//u
//v
//w
//T_static (not total)
//-----Aux Variables
//P
//rho
//cv_1
//...
//cv_NS
//mol_1 (mol/m^3)
//mol_2
//...
//mol_NS

#define N_SUBIT 10

template <class Type>
CompressibleFREqnSet<Type>::CompressibleFREqnSet(SolutionSpace<Type>* space, Param<Type>* p)
{
  this->space = space;
  this->param = p;
  this->chem = new ChemModel<Type>(this->param->casestring, this->param->chemDB);
  this->ownChem = true;
  
  nspecies = chem->nspecies;
  this->neqn = 4 + nspecies;
  this->nauxvars = 2 + 2*this->nspecies;
  //grab and store pointers

  //variable set is not conservative
  this->varsConservative = 0;

  //check to ensure number of massfractions read in (if there) 
  //is the same as number of species in model
  //default to mole fractions first
  if(this->param->molefractions.size() == this->nspecies){
    //convert to massfractions
    Type* molfrac = (Type*)alloca(sizeof(Type)*this->nspecies);
    Type* massfrac = (Type*)alloca(sizeof(Type)*this->nspecies); 
    for(Int i = 0; i < this->nspecies; ++i){
      molfrac[i] = this->param->molefractions[i];
    }
    this->chem->MoleFractionToMassFraction(molfrac, massfrac);
    this->param->massfractions.resize(this->nspecies);
    for(Int i = 0; i < this->nspecies; ++i){
      this->param->massfractions[i] = massfrac[i];
    }
  }
  else if(this->param->massfractions.size() == this->nspecies){
    //do nothing, we want to get massfractions
  }
  else{
    std::stringstream ss;
    ss << "Massfractions/molefractions were read in but do not match ";
    ss << "the number of species in the model" << std::endl;
    ss << "NSPECIES: " << this->nspecies 
       << " NMASSFRACTIONS: " << this->param->massfractions.size() << std::endl;
    ss << "Defaulting to equal distribution of density among species" << std::endl;
    Abort << ss.str();
  }


  this->Qinf = new Type[this->neqn + this->nauxvars];
  this->idata = new DataInfo(this->neqn+this->nauxvars, std::string("variableQ"));
  for(Int i = 0; i < nspecies; i++){
    //set names from chemistry model
    this->idata->AddScalar(i, this->chem->species[i].symbol + "_Density");
  }

  this->idata->AddVector(nspecies, "Velocity");
  this->idata->AddScalar(nspecies+3, "Temperature");
  this->idata->AddScalar(nspecies+4, "Pressure");
  this->idata->AddScalar(nspecies+5, "TotalDensity");
  for(Int i = 0; i < nspecies; i++){
    this->idata->AddScalar(nspecies+6+i, this->chem->species[i].symbol + "_cv");
  }
  for(Int i = 0; i < nspecies; i++){
    this->idata->AddScalar(nspecies+6+nspecies+i, this->chem->species[i].symbol + "_Concentration");
  }
  this->idata->Verify();
  
  //set gradients required
  this->gdata = new DataInfo((this->nspecies+4+this->nspecies)*3, "gradVariableQ");
  for(Int i = 0; i < nspecies; i++){
    this->gdata->AddVector(i*3, "Grad-" + this->chem->species[i].symbol + "_Density");
  }
  this->gdata->AddVector((nspecies+0)*3, "Grad-u");
  this->gdata->AddVector((nspecies+1)*3, "Grad-v");
  this->gdata->AddVector((nspecies+2)*3, "Grad-w");
  this->gdata->AddVector((nspecies+3)*3, "Grad-Temperature");
  for(Int i = 0; i < nspecies; i++){
    this->gdata->AddVector((nspecies+4+i)*3, "Grad-" + this->chem->species[i].symbol + "_Concentration");
  }
  this->gdata->Verify();

  return;
}

template <class Type>
CompressibleFREqnSet<Type>::~CompressibleFREqnSet()
{
  delete [] this->Qinf;
  if(ownChem){
    delete chem;
  }
  return;
}

template <class Type>
void CompressibleFREqnSet<Type>::InitEqnSet()
{
  this->UpdateQinf();

  //allocate solution memory for our Mesh object
  this->space->AddField(*this->idata, FIELDS::STATE_TIME, FIELDS::VAR_EVERYWHERE);
  SolutionField<Type> & field = this->space->GetField("variableQ");
  this->space->q = field.GetData(FIELDS::STATE_NP1);
  this->space->qold = field.GetData(FIELDS::STATE_N);
  this->space->qoldm1 = field.GetData(FIELDS::STATE_NM1);

  //allocate solution memory for the gradients we need
  this->space->AddField(*this->gdata, FIELDS::STATE_NONE, FIELDS::VAR_INTERIOR);
  SolutionField<Type> & gfield  = this->space->GetField("gradVariableQ");
  this->space->qgrad = gfield.GetData(FIELDS::STATE_NONE);

  //set Reynold's number - rho*v*d/mu
  Type rho = this->Qinf[nspecies+5];
  Type V = this->param->GetVelocity(this->space->iter);
  this->param->Re = (rho*this->param->ref_density)*(this->param->ref_velocity*V)*
    this->param->ref_length/this->param->ref_viscosity;
}

template <class Type>
void CompressibleFREqnSet<Type>::Eigensystem(Type* Q, Type* avec, Type vdotn, Type gamma,
					     Type* eigenvalues, Type* T, Type* Tinv, Type beta)
{
  gamma = 0.0;
  Type nx = avec[0];
  Type ny = avec[1];
  Type nz = avec[2];

  Type bm1 = beta - 1.0;
  Type bp1 = beta + 1.0;
  Type oneMBeta = 1.0 - beta;
  Type theta = GetTheta(Q, avec, vdotn);
  Type* rhoi = &Q[0];
  Type p = Q[nspecies+4];
  Type Temp = Q[nspecies+3];
  Type rho = Q[nspecies+5];
  Type* cvi = &Q[nspecies+6];
  Type* c2i = new Type[nspecies];

  Type cv, cp, R, gammaTrash, c2, RhoTrash, PTrash;
  GetFluidProperties(rhoi, Temp, cvi, cv, cp, R, gammaTrash, c2, RhoTrash, PTrash);

  //these are the modified speeds from preconditioning
  Type thetaPrime = 0.5*bp1*theta;
  Type cPrime = 0.5*sqrt(theta*theta*(oneMBeta*oneMBeta) + 4.0*beta*c2);

  //get individual species speed of sounds
  GetSpeciesSpeedOfSound(c2i, Q);
  //TODO: is it required that we reset this for stability reasons... seems c2i should be individual not bulk
  for(Int i = 0; i < nspecies; i++){
    c2i[i] = c2;
  }

  Type m[3];
  Type l[3];
  
  //Compute the remainder of a basis using the face normal
  PerpVectors(avec, l, m);

  Type lx = l[0];
  Type ly = l[1];
  Type lz = l[2];
  Type mx = m[0];
  Type my = m[1];
  Type mz = m[2];

  //Left here as a reminder, this is one by definition
  //Type ScalarTripleProd = lx*my*nz-lx*mz*ny-ly*mx*nz+ly*mz*nx+lz*mx*ny-lz*my*nx;

  Int uloc = nspecies+0;
  Int vloc = nspecies+1;
  Int wloc = nspecies+2;
  Int tloc = nspecies+3;
  for(Int i = 0; i < nspecies+2; i++){
    eigenvalues[i] = theta;
  }
  eigenvalues[wloc] = thetaPrime + cPrime;
  eigenvalues[tloc] = thetaPrime - cPrime;

  Type betam = oneMBeta*0.5;
  Type Xp = theta*betam + cPrime;
  Type Xm = theta*betam - cPrime;

  Int neqn2 = this->neqn*this->neqn;
  Int neqn = this->neqn;
  //Right eigenvectors published in my dissertation NGC
  MemBlank(T, neqn2);

  #pragma ivdep
  for(Int i = 0; i < nspecies; i++){
    T[neqn*i + i] = 1.0;
    T[neqn*i + wloc] = -(rhoi[i]*(c2i[i] + bm1*c2 - theta*bm1*Xm))/(c2i[i]*Xm);
    T[neqn*i + tloc] =  (rhoi[i]*(c2i[i] + bm1*c2 - theta*bm1*Xp))/(c2i[i]*Xp);
  }
  T[neqn*(nspecies+0) + (nspecies+0)] = lx;
  T[neqn*(nspecies+0) + (nspecies+1)] = mx;
  T[neqn*(nspecies+0) + (nspecies+2)] =  nx;
  T[neqn*(nspecies+0) + (nspecies+3)] = -nx;
  
  T[neqn*(nspecies+1) + (nspecies+0)] = ly;
  T[neqn*(nspecies+1) + (nspecies+1)] = my;
  T[neqn*(nspecies+1) + (nspecies+2)] = ny;
  T[neqn*(nspecies+1) + (nspecies+3)] = -ny;
  
  T[neqn*(nspecies+2) + (nspecies+0)] = lz;
  T[neqn*(nspecies+2) + (nspecies+1)] = mz;
  T[neqn*(nspecies+2) + (nspecies+2)] = nz;
  T[neqn*(nspecies+2) + (nspecies+3)] = -nz;
  
  T[neqn*(nspecies+3) + wloc] = -rho*Xm;
  T[neqn*(nspecies+3) + tloc] = rho*Xp;
  
  //Left eigenvectors - from maple inversion
  MemBlank(Tinv, neqn2);
  
  for(Int i = 0; i < nspecies; i++){
    Tinv[i*neqn + i] = 1.0;
    Type KK = -rhoi[i]*(c2i[i]*(Xm + Xp) - bm1*Xm*Xp*theta + bm1*c2*(Xm + Xp));
    Tinv[i*neqn + uloc] = -((ly*mz-lz*my)*KK)/(c2i[i]*Xm*Xp);
    Tinv[i*neqn + vloc] =  ((lx*mz-lz*mx)*KK)/(c2i[i]*Xm*Xp);
    Tinv[i*neqn + wloc] = -((lx*my-ly*mx)*KK)/(c2i[i]*Xm*Xp);
    Tinv[i*neqn + tloc] =  (rhoi[i]*(c2i[i] + bm1*c2))/(rho*c2i[i]*Xm*Xp);
  }
  
  Tinv[neqn*uloc + uloc] =  my*nz-mz*ny;
  Tinv[neqn*uloc + vloc] =-(mx*nz-mz*nx);
  Tinv[neqn*uloc + wloc] =  mx*ny-my*nx;
  Tinv[neqn*uloc + tloc] =  0.0;
  
  Tinv[neqn*vloc + uloc] =-(ly*nz-lz*ny);
  Tinv[neqn*vloc + vloc] =  lx*nz-lz*nx;
  Tinv[neqn*vloc + wloc] =-(lx*ny-ly*nx);
  Tinv[neqn*vloc + tloc] =  0.0;

  Tinv[neqn*wloc + uloc] =  ((Xp)*(ly*mz-lz*my))/(2.0*cPrime);
  Tinv[neqn*wloc + vloc] = -((Xp)*(lx*mz-lz*mx))/(2.0*cPrime);
  Tinv[neqn*wloc + wloc] =  ((Xp)*(lx*my-ly*mx))/(2.0*cPrime);
  Tinv[neqn*wloc + tloc] =  1.0/(2.0*rho*cPrime);

  Tinv[neqn*tloc + uloc] =  ((Xm)*(ly*mz-lz*my))/(2.0*cPrime);
  Tinv[neqn*tloc + vloc] = -((Xm)*(lx*mz-lz*mx))/(2.0*cPrime);
  Tinv[neqn*tloc + wloc] =  ((Xm)*(lx*my-ly*mx))/(2.0*cPrime);
  Tinv[neqn*tloc + tloc] =  1.0/(2.0*rho*cPrime);

#if 0
  //check that T and Tinv really are inverses of each other
  Int yes;
  yes = MatInvCheck(Tinv, T, neqn);
  if(!yes){
    std::cerr << "Eigensystem not valid -- inverse check failed!!" << std::endl;
    std::cerr << "T" << std::endl;
    MatPrintStream(T, neqn, std::cerr);
    std::cerr << "Tinv" << std::endl;
    MatPrintStream(Tinv, neqn, std::cerr);
    std::cerr << "-----------------------------" << std::endl;
  } 
#endif

  delete [] c2i;

}

template <class Type>
void CompressibleFREqnSet<Type>::RoeFlux(Type* QL, Type* QR, Type* avec, Type vdotn,
					 Type gamma, Type* flux, Type beta)
{
  //no RoeFlux here... too much work :)
  HLLCFlux(QL, QR, avec, vdotn, gamma, flux, beta);
}

template <class Type>
void CompressibleFREqnSet<Type>::HLLCFlux(Type* QL, Type* QR, Type* avec, Type vdotn, 
					  Type gamma, Type* flux, Type beta)
{
  Int i;
  gamma = 0.0; //nullify gamma to trigger errors should we not set it correctly here
  Int neqn = this->neqn;

  //get left state
  Type* rhoiL = &QL[0];
  Type uL = QL[nspecies];
  Type vL = QL[nspecies+1];
  Type wL = QL[nspecies+2];
  Type TL = QL[nspecies+3];
  Type pL = QL[nspecies+4];
  Type pgL = pL - this->Pref;
  Type rhoL = QL[nspecies+5];
  Type* cviL = &QL[nspecies+6];
  Type cvL, cpL, RL, gammaL, c2L, RhoTrashL, PTrashL;
  GetFluidProperties(rhoiL, TL, cviL, cvL, cpL, RL, gammaL, c2L, RhoTrashL, PTrashL);
  Type HTL = GetTotalEnthalpy(QL);
  Type ETL = HTL - pL;
#if 0
  Type hV2L = 0.5*(uL*uL + vL*vL + wL*wL);
  Type hL = HTL/rhoL - hv2L
  Type etL = ETL/rhoL;
  Type eL = etL - hV2L;
  Type c2L = (gammaL - 1.0)*((hL - eL) - hV2L + cvL*TL);
#endif
  Type cL = sqrt(c2L);
  Type thetaL = GetTheta(QL, avec, vdotn);

  //get right state
  Type* rhoiR = &QR[0];
  Type uR = QR[nspecies];
  Type vR = QR[nspecies+1];
  Type wR = QR[nspecies+2]; 
  Type TR = QR[nspecies+3];
  Type pR = QR[nspecies+4];
  Type pgR = pR - this->Pref;
  Type rhoR = QR[nspecies+5];
  Type* cviR = &QR[nspecies+6];
  Type cvR, cpR, RR, gammaR, c2R, RhoTrashR, PTrashR;
  GetFluidProperties(rhoiR, TR, cviR, cvR, cpR, RR, gammaR, c2R, RhoTrashR, PTrashR);
  Type HTR = GetTotalEnthalpy(QR);
  Type ETR = HTR - pR;
#if 0
  Type hV2R = 0.5*(uR*uR + vR*vR + wR*wR);
  Type hR = HTR/rhoR - hv2R
  Type etR = ETR/rhoR;
  Type eR = etR - hV2R;
  Type c2R = (gammaR - 1.0)*((hR - eR) - hV2R + cvR*TR);
#endif
  Type cR = sqrt(c2R);
  Type thetaR = GetTheta(QR, avec, vdotn);

  Type rho = sqrt(rhoL*rhoR);
  Type sigma = rho/(rhoL + rho);
  Type* roeQ = new Type[neqn + this->nauxvars];
  //compute Roe averaged variables
  for(i = 0; i < nspecies; i++){
    roeQ[i] = rhoiL[i] + sigma*(rhoiR[i] - rhoiL[i]);
  }
  Type u, v, w;
  u = roeQ[nspecies + 0] = uL + sigma*(uR - uL);
  v = roeQ[nspecies + 1] = vL + sigma*(vR - vL);
  w = roeQ[nspecies + 2] = wL + sigma*(wR - wL);
  roeQ[nspecies + 3] = TL + sigma*(TR - TL);
  ComputeAuxiliaryVariables(roeQ);

  Type theta = GetTheta(roeQ, avec, vdotn);
  
  Type* rhoi = &roeQ[0];
  Type* cvi = &roeQ[nspecies+6];
  Type T = roeQ[nspecies+3];
  Type p = roeQ[nspecies+4];
  rho = roeQ[nspecies+5];
  Type cv, cp, R, gammaTrash, c2, RhoTrash, PTrash;
  GetFluidProperties(rhoi, T, cvi, cv, cp, R, gammaTrash, c2, RhoTrash, PTrash);
  //this can be done more accurately, like:
  //Type hV2 = 0.5*(u*u + v*v + w*w);
  //c2 = (gamma -1.0)*(h - hV2 + cv*T - specific_energy);
  Type c = sqrt(c2);

  //compute max and min eigenvalues, these are modified as in Ashish Gupta dissertation p.55
  Type oneMBeta = 1.0 - beta;
  Type thetaPrime = theta*(1.0 + beta)*0.5;
  Type cPrime = 0.5*sqrt(theta*theta*(oneMBeta*oneMBeta) + 4.0*beta*c2);

  Type thetaLPrime = thetaL*(1.0 + beta)*0.5;
  Type thetaRPrime = thetaR*(1.0 + beta)*0.5;
  Type cLPrime = 0.5*sqrt(thetaL*thetaL*(oneMBeta*oneMBeta) + 4.0*beta*c2L);
  Type cRPrime = 0.5*sqrt(thetaR*thetaR*(oneMBeta*oneMBeta) + 4.0*beta*c2R);

  Type eig5L = thetaLPrime - cLPrime;
  Type eig4R = thetaRPrime + cRPrime;
  Type eig4 = thetaPrime + cPrime;
  Type eig5 = thetaPrime - cPrime;

  //compute HLLC SL, SR, and SM
  Type SL = MIN(eig5L, eig5);
  Type SR = MAX(eig4R, eig4);
  //SM is the contact wave speed, we don't modify this one
  Type SM = (pgR - pgL + rhoL*thetaL*(SL - thetaL) - rhoR*thetaR*(SR-thetaR))/
    (rhoL*(SL-thetaL) - rhoR*(SR-thetaR));

  Type* Qflux = new Type[neqn+this->nauxvars];
  Type pStar, omega, const1, const2;
  pStar = 0.0; //silence warning
  if(real(SL) >= 0.0){
    for (i = 0; i < nspecies; i++){
      Qflux[i] = rhoiL[i];
    }
    Qflux[nspecies] = rhoL*uL;
    Qflux[nspecies+1] = rhoL*vL;
    Qflux[nspecies+2] = rhoL*wL;
    Qflux[nspecies+3] = ETL;
    pStar = pgL;
    SM = thetaL;
  }
  else if(real(SR) <= 0.0){
    for (i = 0; i < nspecies; i++){
      Qflux[i] = rhoiR[i];
    }
    Qflux[nspecies] = rhoR*uR;
    Qflux[nspecies+1] = rhoR*vR;
    Qflux[nspecies+2] = rhoR*wR;
    Qflux[nspecies+3] = ETR;
    pStar = pgR;
    SM = thetaR;
  }
  else if((real(SL) <= 0.0) && (real(SM) >= 0.0)){
    pStar = pgL + rhoL*(thetaL - SL)*(thetaL - SM);
    
    omega = 1.0/(SL - SM);
    const1 = SL - thetaL;
    const2 = pStar - pgL;
    
    for (i = 0; i < nspecies; i++){
      Qflux[i] = omega*const1*rhoiL[i];
    }
    Qflux[nspecies] = omega*(const1*rhoL*uL + const2*avec[0]);
    Qflux[nspecies+1] = omega*(const1*rhoL*vL + const2*avec[1]);
    Qflux[nspecies+2] = omega*(const1*rhoL*wL + const2*avec[2]);
    Qflux[nspecies+3] = omega*(const1*ETL - pgL*thetaL + (pStar*SM)) + this->Pref*(SM - thetaL)*omega;
  }
  else if((real(SM) <= 0.0) && (real(SR) >= 0.0)){
    pStar = pgR + rhoR*(thetaR - SR)*(thetaR - SM);
    
    omega = 1.0/(SR - SM);
    const1 = SR - thetaR;
    const2 = pStar - pgR;
    
    for (i = 0; i < nspecies; i++){
      Qflux[i] = omega*const1*rhoiR[i];
    }
    Qflux[nspecies] = omega*(const1*rhoR*uR + const2*avec[0]);
    Qflux[nspecies+1] = omega*(const1*rhoR*vR + const2*avec[1]);
    Qflux[nspecies+2] = omega*(const1*rhoR*wR + const2*avec[2]);
    Qflux[nspecies+3] = omega*(const1*ETR - pgR*thetaR + (pStar*SM)) + this->Pref*(SM-thetaR)*omega;
  }
  else{
    std::stringstream ss;
    ss << "HLLC: Should never be here!!" << std::endl;
    ss << "beta = " << beta << std::endl;
    ss << "SM = " << real(SM) << std::endl;
    ss << "SL = " << real(SL) << std::endl;
    ss << "SR = " << real(SR) << std::endl;
    ss << "Theta = " << theta << std::endl;
    ss << "Thetaprime = " << thetaPrime << std::endl;
    ss << "ThetaLPrime = " << thetaLPrime << std::endl;
    ss << "ThetaRPrime = " << thetaRPrime << std::endl;
    ss << "cPrime = " << cPrime << std::endl;
    ss << "cLPrime = " << cLPrime << std::endl;
    ss << "cRPrime = " << cRPrime << std::endl;
    ss << "ThetaR = " << thetaR << std::endl;
    ss << "ThetaL = " << thetaL << std::endl;
    ss << "c2L = " << c2L << std::endl;
    ss << "c2R = " << c2R << std::endl;
    ss << "gammaL = " << gammaL << std::endl;
    ss << "gammaR = " << gammaR << std::endl;
    ss << "Left state: " << std::endl;
    for(i = 0; i < nspecies; i++){
      ss << "rhoi[" << i << "] = " << QL[i] << std::endl;
    }
    ss << "u = " << QL[nspecies] << std::endl;
    ss << "v = " << QL[nspecies+1] << std::endl;
    ss << "w = " << QL[nspecies+2] << std::endl;
    ss << "T = " << QL[nspecies+3] << std::endl;
    ss << "P = " << QL[nspecies+4] << std::endl;
    ss << "rho = " << QL[nspecies+5] << std::endl;
    ss << "R = " << RL << std::endl;
    for(i = 0; i < nspecies; i++){
      ss << "cvi[" << i << "] = " << QL[nspecies+6+i] << std::endl;
    }
    ss << "c2: " << c2L << std::endl;
    ss << "Right state: " << std::endl;
    for(i = 0; i < nspecies; i++){
      ss << "rhoi[" << i << "] = " << QR[i] << std::endl;
    }
    ss << "u = " << QR[nspecies] << std::endl;
    ss << "v = " << QR[nspecies+1] << std::endl;
    ss << "w = " << QR[nspecies+2] << std::endl;
    ss << "T = " << QR[nspecies+3] << std::endl;
    ss << "P = " << QR[nspecies+4] << std::endl;
    ss << "rho = " << QR[nspecies+5] << std::endl;
    ss << "R = " << RR << std::endl;
    for(i = 0; i < nspecies; i++){
      ss << "cvi[" << i << "] = " << QR[nspecies+6+i] << std::endl;
    }
    ss << "c2: " << c2R << std::endl;
    ss << "WARNING: HLLC died\n";
    Abort.SetSoftAbort(ss.str());
  }

  //compute the flux
  Type area = avec[3];
  Type* rhoiF = &Qflux[0];
  Type Et = Qflux[nspecies+3];
  //Type Ht = Et + (pStar+this->Pref);
  theta = SM;
  Type thetabar = theta - vdotn;
  for(i = 0; i < nspecies; i++){
    flux[i] = area*rhoiF[i]*theta;
  }
  flux[nspecies] = area*(Qflux[nspecies]*theta + pStar*avec[0]);
  flux[nspecies+1] = area*(Qflux[nspecies+1]*theta + pStar*avec[1]);
  flux[nspecies+2] = area*(Qflux[nspecies+2]*theta + pStar*avec[2]);
  //flux[nspecies+3] = area*(Ht*theta - vdotn*(pStar+this->Pref));
  flux[nspecies+3] = area*(Et*theta + pStar*thetabar) + this->Pref*thetabar*area;

  //pseudo-2D simulations
  if(this->param->symmetry2D > 0){
    if(this->param->symmetry2D == 1){
      flux[nspecies] = 0.0;
    }
    else if(this->param->symmetry2D == 2){
      flux[nspecies+1] = 0.0;
    }
    else if(this->param->symmetry2D == 3){
      flux[nspecies+2] = 0.0;
    }
  }
  
  delete [] roeQ;
  delete [] Qflux;

  return;
}

template <class Type>
void CompressibleFREqnSet<Type>::ViscousFlux(Type* Q, Type* grad, Type* avec, Type mut, Type* flux)
{
  Int i;
  //gradient is passed in with velocity components and then temperature
  //parse for readability
  Type ux, uy, uz, vx, vy, vz, wx, wy, wz;
  Type Tx, Ty, Tz;
  Int offset = nspecies*3;
  ux = grad[offset + 0];
  uy = grad[offset + 1];
  uz = grad[offset + 2];
  vx = grad[offset + 3];
  vy = grad[offset + 4];
  vz = grad[offset + 5];
  wx = grad[offset + 6];
  wy = grad[offset + 7];
  wz = grad[offset + 8];
  Tx = grad[offset + 9];
  Ty = grad[offset + 10];
  Tz = grad[offset + 11];

  Type u, v, w;
  u = Q[nspecies];
  v = Q[nspecies+1];
  w = Q[nspecies+2];

  //compute cp for the mixture
  Type* rhoi = &Q[0];
  Type rho = Q[nspecies+5];
  Type* cvi = &Q[nspecies+6];
  Type P = GetPressure(Q);
  Type T = GetTemperature(Q);
  Type cv, cp, R, gamma, c2, RhoTrash, PTrash;
  GetFluidProperties(rhoi, T, cvi, cv, cp, R, gamma, c2, RhoTrash, PTrash);

  //compute actual mixture viscosity and add turbulent viscosity
  Type mu = GetMolecularViscosity(rhoi, T);
  Type tmut = (mu + mut);
  Type tauxx, tauyy, tauzz, tauxy, tauxz, tauyz;
  Type fact = 2.0/3.0;
  Type tauxn, tauyn, tauzn;
  tauxx = 2.0*fact*ux - fact*vy - fact*wz;
  tauyy = 2.0*fact*vy - fact*ux - fact*wz;
  tauzz = 2.0*fact*wz - fact*ux - fact*vy;
  tauxy = uy + vx;
  tauxz = uz + wx;
  tauyz = vz + wy;

  Type Re = this->param->Re;
  Type RK = avec[3]/Re;
  Type RKT = RK*tmut;

  //heat flux term
  Type k = -GetThermalConductivity(rhoi, T);
  Type Tn = Tx*avec[0] + Ty*avec[1] + Tz*avec[2];

  //compute turbulent thermal conductivity
  //this term represents a turbulent contribution to thermal conductivity
  //ref Anderson p.624
  Type kT = (cp * mut)/this->param->PrT;
  k -= kT;

  tauxn = -tauxx*avec[0] - tauxy*avec[1] - tauxz*avec[2];
  tauyn = -tauxy*avec[0] - tauyy*avec[1] - tauyz*avec[2];
  tauzn = -tauxz*avec[0] - tauyz*avec[1] - tauzz*avec[2];

  for(i = 0; i < nspecies; i++){
    flux[i] = 0.0;
  }
  flux[nspecies+0] = RKT*(tauxn);
  flux[nspecies+1] = RKT*(tauyn);
  flux[nspecies+2] = RKT*(tauzn);
  flux[nspecies+3] = RKT*(tauxn*u + tauyn*v + tauzn*w) + RK*k*Tn;

  //pseudo-2D simulations
  if(this->param->symmetry2D > 0){
    if(this->param->symmetry2D == 1){
      flux[nspecies] = 0.0;
    }
    else if(this->param->symmetry2D == 2){
      flux[nspecies+1] = 0.0;
    }
    else if(this->param->symmetry2D == 3){
      flux[nspecies+2] = 0.0;
    }
  }

  return;
}

template <class Type>
Type CompressibleFREqnSet<Type>::GetDensity(Type* Q)
{
  return (Q[nspecies+5]);
}

template <class Type>
Type CompressibleFREqnSet<Type>::GetPressure(Type* Q)
{
  return (Q[nspecies+4]);
}

template <class Type>
Type CompressibleFREqnSet<Type>::GetTemperature(Type* Q)
{
  return (Q[nspecies+3]);
}

template <class Type>
Type CompressibleFREqnSet<Type>::GetGamma(Type* Q)
{
  Type* rhoi = &Q[0];
  Type* rhoiDim = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; ++i){
    rhoiDim[i] = rhoi[i]*this->param->ref_density;
  }
  Type TDim = GetTemperature(Q)*this->param->ref_temperature;
  Type cpDim = this->chem->GetCp(rhoiDim, TDim);
  Type cvDim = this->chem->GetCv(rhoiDim, TDim);
  return (cpDim/cvDim);
}

template <class Type>
Type CompressibleFREqnSet<Type>::GetCp(Type* Q, Type gamma)
{
  Type P = GetPressure(Q);
  Type Pinf = GetPressure(this->Qinf);
  Type V = this->param->GetVelocity(this->space->iter);
  return  ((P - Pinf)/(0.5*V*V));
}

template <class Type>
Int CompressibleFREqnSet<Type>::GetVelocityGradLocation()
{
  return (nspecies);
}

template <class Type>
Int CompressibleFREqnSet<Type>::GetMomentumLocation()
{
  return (nspecies);
}

template <class Type>
Int CompressibleFREqnSet<Type>::GetGradientsLocation(std::vector<Int>& gradientLoc)
{
  //densities
  for(Int i = 0; i < nspecies; i++){
    gradientLoc.push_back(i);
  }
  //velocity
  gradientLoc.push_back(nspecies+0);
  gradientLoc.push_back(nspecies+1);
  gradientLoc.push_back(nspecies+2);
  //temperature
  gradientLoc.push_back(nspecies+3);
  //moles/concentration
  for(Int i = 0; i < nspecies; i++){
    gradientLoc.push_back(nspecies+nspecies+6+i);
  }
  
  return gradientLoc.size();
}

template <class Type>
Int CompressibleFREqnSet<Type>::BadExtrapolation(Type* Q)
{
  for(Int i = 0; i < nspecies; i++){
    if(real(Q[i]) < 0.0){
      return true;
    }
  }
  if(real(GetTotalEnergy(Q)) <= 0.0){
    return true;   
  }
  if(real(GetPressure(Q)) < 1.0e-10){
    return true;
  }
  if(real(GetTemperature(Q)) < 1.0e-10){
    return true;
  }
  return false;
}

template <class Type>
void CompressibleFREqnSet<Type>::ExtrapolateVariables(Type* Qho, const Type* Q, const Type* dQedge, 
						      const Type* gradQ, const Type* dx, const Type* limiter)
{
  //densities
  #pragma ivdep
  for(Int i = 0; i < nspecies; i++){
    Qho[i] = Q[i] + this->ExtrapolateCorrection(dQedge[i], &gradQ[i*3], dx)*limiter[i];
  }
  //velocity
  Qho[nspecies+0] = Q[nspecies+0] + this->ExtrapolateCorrection(dQedge[nspecies+0], &gradQ[(nspecies+0)*3], dx)*
    limiter[nspecies+0];
  Qho[nspecies+1] = Q[nspecies+1] + this->ExtrapolateCorrection(dQedge[nspecies+1], &gradQ[(nspecies+1)*3], dx)*
    limiter[nspecies+1];
  Qho[nspecies+2] = Q[nspecies+2] + this->ExtrapolateCorrection(dQedge[nspecies+2], &gradQ[(nspecies+2)*3], dx)*
    limiter[nspecies+2];
  //temperature
  Qho[nspecies+3] = Q[nspecies+3] + this->ExtrapolateCorrection(dQedge[nspecies+3], &gradQ[(nspecies+3)*3], dx)*
    limiter[nspecies+3];
}

template <class Type>
void CompressibleFREqnSet<Type>::ComputeAuxiliaryVariables(Type* Q)
{
  Int i;
  Type R;
  Type* rho = &Q[nspecies+5];
  Type* rhoi = &Q[0];
  Type* cvi = &Q[nspecies+6];
  Type* mol = &Q[nspecies+nspecies+6];
  Type T = Q[nspecies+3];

  Type* rhoiDim = (Type*)alloca(sizeof(Type)*nspecies);
  
  //compute rho and R
  *rho = 0.0;
  R = 0.0;
  for(i = 0; i < nspecies; i++){
    *rho += rhoi[i];
    rhoiDim[i] = rhoi[i]*this->param->ref_density;
    R += this->chem->species[i].R*(rhoi[i]*this->param->ref_density);
  }
  //R is dimensional after this
  R /= (*rho)*this->param->ref_density;

  //set P via EOS
  Type TDim = T*this->param->ref_temperature;
  Type pDim = this->chem->GetP(rhoiDim, TDim);
  Q[nspecies+4] = pDim/this->param->ref_pressure;

  //compute cvi
  Type rhoDim = (*rho)*this->param->ref_density;
  Type* cpiDim = (Type*)alloca(sizeof(Type)*nspecies);
  Type s_ref = (this->param->ref_velocity * this->param->ref_velocity / 
		this->param->ref_temperature);
  for(i = 0; i < nspecies; i++){
    cpiDim[i] = chem->species[i].GetCp(TDim);
    cvi[i] = chem->eos[i]->GetCv(cpiDim[i], chem->species[i].R, rhoDim, pDim, TDim) / s_ref;
    if(real(cvi[i]) < real(0.0)){
      std::cerr << "WARNING: negative Cv computed" << std::endl;
      std::cerr << "Species: " << chem->species[i].symbol << std::endl;
      std::cerr << "Rspecific: " << chem->species[i].R << std::endl;
      std::cerr << "Cp: " << cpiDim[i] << std::endl;
      std::cerr << "Cv: " << cvi[i]*s_ref << std::endl;
      std::cerr << "T (K): " << TDim << std::endl;
    }
  }

  //compute concentrations given in mol/m^3
  chem->MassToMole(rhoi, mol);
  
#if 0
  std::cout << "In ComputeAuxiliaryVariables()" << std::endl;
  std::cout << "TDim: " << TDim << std::endl;
  std::cout << "PDim: " << pDim << std::endl;
  std::cout << "P: " << Q[nspecies+4] << std::endl;
  std::cout << "T: " << T << std::endl;
  std::cout << "R: " << R << std::endl;
  std::cout << "rhoDim: " << (*rho)*this->param->ref_density << std::endl;
  std::cout << "rho: " << (*rho) << std::endl;
#endif
}

template <class Type>
Type CompressibleFREqnSet<Type>::ComputePressure(Type* Q, Type gamma)
{
  Type* rhoi = &Q[0];
  Type TDim = Q[nspecies+3]*this->param->ref_temperature;
  Type* rhoiDim = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; i++){
    rhoiDim[i] = rhoi[i]*this->param->ref_density;
  }
  //set P via EOS
  Type pDim = this->chem->GetP(rhoiDim, TDim);
  return pDim/this->param->ref_pressure;
}


template <class Type>
void CompressibleFREqnSet<Type>::SetInitialConditions()
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
      pywrap.SetInitialConditions(this->Qinf, neqn, nauxvars, &this->space->q[i*(neqn+nauxvars)], &m->xyz[i*3]);
      ComputeAuxiliaryVariables(&this->space->q[i*(neqn+nauxvars)]);
    }
#else
    Abort << "Python not built with solver";
#endif
  }
  else{
    //set all the nodes interior and phantom to Qinf
    for(i = 0; i < (nnode+gnode+nbnode); i++){
      for(j = 0; j < neqn; j++){
	this->space->q[i*(neqn+nauxvars) + j] = this->Qinf[j];
      }
      ComputeAuxiliaryVariables(&this->space->q[i*(neqn+nauxvars)]);
    }
  }
}


template <class Type>
void CompressibleFREqnSet<Type>::ApplyDQ(Type* dQ, Type* Q, Type* xyz)
{
  //specialized implementation b/c we know what variables
  //should be clipped if negative
  Int i;
  Type minT = 1.0e-10;
  Type minRhoWarn = 1.0e-10;
  Type rho;
  //clip density to 0.0 if it will be negative
  for(i = 0; i < nspecies; i++){
    rho = Q[i] + dQ[i];
    if(real(rho) < 0.0){
      //this check silences warnings when a species density is very near to zero
      if(real(rho) < -real(minRhoWarn)){ 
	std::cerr << "WARNING: density clip species " << i;
	std::cerr << " dRho: " << dQ[i] << " newRho: " << rho; 
	std::cerr << "  Coordinates: " << xyz[0] 
		  << " " << xyz[1] << " " 
		  << xyz[2] << std::endl;
      }
      if(Q[i] == 0.0){
	Q[i] = 0.0;
      }
      else{
	//refuse to update.. does this work?
	Q[i] = Q[i];
      }
    }
    else{
      Q[i] += dQ[i];
    }
  }
  //check for negative temperature
  Type projectedT = Q[nspecies+3] + dQ[nspecies+3];
  if(real(projectedT) < 0.0){
    std::cerr << "WARNING: temperature clip - Tnew: " << projectedT;
    std::cerr << " dT: " << dQ[nspecies+3];
    std::cerr << "  Coordinates: " << xyz[0] 
	      << " " << xyz[1] << " " << xyz[2] << std::endl;
    
    Q[nspecies+3] = minT;
  }
  else{
    Q[nspecies+3] = projectedT;
  }
  //update velocities
  Q[nspecies] += dQ[nspecies];
  Q[nspecies+1] += dQ[nspecies+1];
  Q[nspecies+2] += dQ[nspecies+2];

  ComputeAuxiliaryVariables(Q);
}

template <class Type>
void CompressibleFREqnSet<Type>::GetFarfieldBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta)
{
  Int neqn = this->neqn;
  Int nvars = neqn + this->nauxvars;
  Type* qavg = new Type[nvars];
  Type* eigenvals = new Type[neqn];
  Type* Tinv = new Type[neqn*neqn];
  Type* T = new Type[neqn*neqn];
  Type* rhs = new Type[neqn];
  Type* scr = new Type[neqn];
  Type* ql = new Type[neqn];
  Type* qinf = new Type[neqn];

  //extract grid speeds in normal direction
  //vdotn = -dot(gridspeeds, avec)
  Type u, v, w;
  u = vdotn*avec[0];
  v = vdotn*avec[1];
  w = vdotn*avec[2];

  Type gamma = 0.0;

  for(Int subit = 0; subit < N_SUBIT; subit++){
    for(Int i = 0; i < neqn; i++){
      qavg[i] = 0.5*(QL[i] + QR[i]);
    }
    ComputeAuxiliaryVariables(qavg);
    
    Eigensystem(qavg, avec, vdotn, gamma, eigenvals, T, Tinv, beta);

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
      //build a vector of qinf and Ql which are in rho, u, P form which is the
      //system for which this eigenspace was developed
      memcpy(ql, QL, sizeof(Type)*neqn);
      memcpy(qinf, Qinf, sizeof(Type)*neqn);
      //store temperature as starting point for Newton solve
      Type Tguess = ql[neqn-1];
      //change last row to pressure instead of temperature
      ql[neqn-1] = GetPressure(QL);
      qinf[neqn-1] = GetPressure(Qinf);

      for(Int i = 0; i < neqn; i++){
	rhs[i] = 0.0;
	for(Int j = 0; j < neqn; j++){
	  rhs[i] += Tinv[i*neqn + j]*(real(eigenvals[i]) >= 0.0 ? ql[j] : qinf[j]);
	}
      }
      
      Int p[neqn];
      //LU(Tinv, p, neqn);
      //LuSolve(Tinv, rhs, p, scr, neqn);
      //copy solution from rhs to QR
      //memcpy(QR, rhs, neqn*sizeof(Type));
      MatVecMult(T, rhs, QR, neqn);

      //set the densities based on first eigenvalue, this is to fix some stability
      //issues with zero density inflow species
#if 0
      if(real(eigenvals[0]) >= 0.0){
	memcpy(QR, QL, nspecies*sizeof(Type));
      }
      else{
	memcpy(QR, Qinf, nspecies*sizeof(Type));
      }
#endif  
      //check for any negative densities and push them back to zero
      for(Int i = 0; i < nspecies; i++){
	if(real(QR[i]) < 0.0){
	  QR[i] = 0.0;
	}
      }

      // convert pressure solution back to temperature for storage, use Newton's method
      // to solve nonlinear mixture EOS
      Type* rhoi = &QR[0];
      Type pgoal = QR[neqn-1];
      Type T = NewtonFindTGivenP(rhoi, pgoal, Tguess);
      QR[neqn-1] = T;
    }
  }

  delete [] qavg;
  delete [] eigenvals;
  delete [] Tinv;
  delete [] T;
  delete [] rhs;
  delete [] scr;
  delete [] ql;
  delete [] qinf;

  return;
}

template <class Type>
void CompressibleFREqnSet<Type>::GetInviscidWallBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta)
{ 

  //TODO: check this for functionality.. may not be quite right.. :(
  Int i;
  Int neqn = this->neqn;
  Int nvars = neqn + this->nauxvars;
  Type gamma = 0.0;

  if(!this->param->no_cvbc){
    Type* qavg = new Type[nvars];
    Type* eigenvals = new Type[neqn];
    Type* Tinv = new Type[neqn*neqn];
    Type* T = new Type[neqn*neqn];
    Type* rhs = new Type[neqn];
    Type* scr = new Type[neqn];
    Type* ql = new Type[neqn];

    for(Int subit = 0; subit < N_SUBIT; subit++){
      for(Int i = 0; i < neqn; i++){
	qavg[i] = 0.5*(QL[i] + QR[i]);
      }
      ComputeAuxiliaryVariables(qavg);
      
      Eigensystem(qavg, avec, vdotn, gamma, eigenvals, T, Tinv, beta);
      
      //build a vector of Ql which are in [rho, u, v, w, P] form which is the
      //system for which this eigenspace was developed
      memcpy(ql, QL, sizeof(Type)*neqn);
      //store temperature as starting point for Newton solve
      Type Tguess = ql[neqn-1];
      //change last row to pressure instead of temperature
      ql[neqn-1] = GetPressure(QL);
      
      //use internal states for the characteristic equation
      for(Int i = 0; i < neqn; i++){
	rhs[i] = 0.0;
	for(Int j = 0; j < neqn; j++){
	  rhs[i] += Tinv[i*neqn + j]*ql[j];
	}
      }
      
      //hard code the last equation to enforce theta = 0, wall velocity = flow velocity
      for(Int i = 0; i < nspecies; i++){
	Tinv[(neqn-1)*neqn + i] = 0.0;
      }
      Tinv[(neqn-1)*neqn + nspecies] = avec[0];
      Tinv[(neqn-1)*neqn + nspecies+1] = avec[1];
      Tinv[(neqn-1)*neqn + nspecies+2] = avec[2];
      Tinv[(neqn-1)*neqn + nspecies+3] = 0.0;
      rhs[neqn-1] = vdotn;

      //solve the eigensystem for [rhoi, u, v, w, P] on boundary
      Int p[neqn];
      LU(Tinv, p, neqn);
      LuSolve(Tinv, rhs, p, scr, neqn);
      //copy solution from rhs to QR
      memcpy(QR, rhs, neqn*sizeof(Type));
      //MatVecMult(T, rhs, QR, neqn);

      // convert pressure solution back to temperature for storage, use Newton's method
      // to solve nonlinear mixture EOS
      Type* rhoi = &QR[0];
      Type pgoal = QR[neqn-1];
      Type T = NewtonFindTGivenP(rhoi, pgoal, Tguess);
      QR[neqn-1] = T;
    }

    delete [] qavg;
    delete [] eigenvals;
    delete [] Tinv;
    delete [] T;
    delete [] rhs;
    delete [] scr;
    delete [] ql;
  }
  else{
    Type* QLmod = (Type*)alloca(sizeof(Type)*(nvars));
    memcpy(QLmod, QL, sizeof(Type)*nvars);
    
    //extract grid speeds in normal direction
    //vdotn = -dot(gridspeeds, avec)
    Type u, v, w;
    u = vdotn*avec[0];
    v = vdotn*avec[1];
    w = vdotn*avec[2];
    
    //add grid speeds to the velocity component of the internal Q
    QLmod[nspecies] += u;
    QLmod[nspecies+1] += v;
    QLmod[nspecies+2] += w;
    ComputeAuxiliaryVariables(QLmod);
    for(i = 0; i < this->neqn; i++) QR[i] = QLmod[i];
    //mirror the velocities
    MirrorVector(&QLmod[nspecies], avec, &QR[nspecies]);
  }

#if 0
  Type* tempspace = (Type*)alloca(sizeof(Type)*(nvars));
  //check for leakage through the surface
  std::cout << std::endl;
  for(i = 0; i < neqn; i++){
    tempspace[i] = (QL[i] + QR[i])/2.0;
  }
  ComputeAuxiliaryVariables(tempspace);
  ui = tempspace[nspecies];
  vi = tempspace[nspecies+1];
  wi = tempspace[nspecies+2];
  std::cout << "theta " << GetTheta(tempspace, avec, vdotn) << std::endl;
  std::cout << std::endl;
#endif  
}

template <class Type>
void CompressibleFREqnSet<Type>::GetInternalInflowBoundaryVariables(Type* QL, Type* QR, Type* Qinf, 
								    Type pressure, Type* densities,
								    Type* flowDirection, Type velocity)
{
  Int i;
  //copy over partial densities
  if(densities == NULL){
    for(i = 0; i < nspecies; i++){
      QR[i] = Qinf[i];
    }
  }
  else{
    for(i = 0; i < nspecies; i++){
      QR[i] = QL[i] = densities[i];
    }
  }
  //copy 3 velocities from bc
  QR[nspecies] = QL[nspecies] = flowDirection[0]*velocity;
  QR[nspecies+1] = QL[nspecies+1] = flowDirection[1]*velocity;
  QR[nspecies+2] = QL[nspecies+2] = flowDirection[2]*velocity;
  //use the pressure from inside the domain, allow it to float
  Type Pgoal = GetPressure(QL);
  Type Tguess = GetTemperature(QL);
  // convert pressure solution back to temperature for storage, use Newton's method
  // to solve nonlinear mixture EOS
  Type* rhoi = &QR[0];
  Type T = NewtonFindTGivenP(rhoi, Pgoal, Tguess);
  QR[nspecies+3] = T;
}

template <class Type>
void CompressibleFREqnSet<Type>::ModifyInternalInflowResidual(Type* res, Type* densities, 
							      Type* flowDirection, Type velocity)
{
  //hardset densities and 
  if(densities != NULL){
    for(Int i = 0; i < nspecies; i++){
      res[i] = 0.0;
    }
  }
  res[nspecies+0] = 0.0;
  res[nspecies+1] = 0.0;
  res[nspecies+2] = 0.0;

  //let pressure float
}

template <class Type>
void CompressibleFREqnSet<Type>::ModifyInternalInflowJacobian(Int cvid, CRS<Type>* crs, 
							      Type* densities, Type* flowDirection, 
							      Type velocity)
{
  if(densities != NULL){
    for(Int i = 0; i < nspecies; i++){
      crs->A->BlankSubRow(cvid, i);
    }
  }
  crs->A->BlankSubRow(cvid, nspecies+0);
  crs->A->BlankSubRow(cvid, nspecies+1);
  crs->A->BlankSubRow(cvid, nspecies+2);
}

template <class Type>
Type CompressibleFREqnSet<Type>::GetTotalEnergy(Type* Q)
{
  Type Et = 0.0;
  Type h, E;
  Type* X = (Type*)alloca(sizeof(Type)*nspecies);
  Type* hi = (Type*)alloca(sizeof(Type)*nspecies);
  Type* rhoi = &Q[0];
  Type T = Q[nspecies+3];
  Type P = Q[nspecies+4];
  Type rho = Q[nspecies+5];
  Type u = Q[nspecies];
  Type v = Q[nspecies+1];
  Type w = Q[nspecies+2];
  Type v2 = u*u + v*v + w*w;
  Int i;
  for(i = 0; i < nspecies; i++){
    X[i] = rhoi[i]/rho;
  }

  h = this->chem->GetSpecificEnthalpy(X, T*this->param->ref_temperature, hi)/
    this->param->ref_enthalpy;

  E = h*rho - P;
  // add on the kinetic part for total energy
  Et = E + 0.5*rho*v2;

  return Et;
}

template <class Type>
Type CompressibleFREqnSet<Type>::GetTotalEnthalpy(Type* Q)
{
  Type Ht = 0.0;
  Type h;
  Type* X = (Type*)alloca(sizeof(Type)*nspecies);
  Type* hi = (Type*)alloca(sizeof(Type)*nspecies);
  Type* rhoi = &Q[0];
  Type T = Q[nspecies+3];
  Type rho = Q[nspecies+5];
  Type u = Q[nspecies];
  Type v = Q[nspecies+1];
  Type w = Q[nspecies+2];
  Type v2 = u*u + v*v + w*w;
  Int i;
  for(i = 0; i < nspecies; i++){
    X[i] = rhoi[i]/rho;
  }

  h = this->chem->GetSpecificEnthalpy(X, T*this->param->ref_temperature, hi)/
    this->param->ref_enthalpy;

  Ht = h*rho + 0.5*rho*v2;

  return Ht;
}

template <class Type>
void CompressibleFREqnSet<Type>::SourceTerm(Type* Q, Type vol, Type* source)
{
  Int i;
  Type* rhoi = (Type*)alloca(sizeof(Type)*nspecies);
  Type* wdot = (Type*)alloca(sizeof(Type)*nspecies);

  MemBlank(source, this->neqn);

  if(this->param->rxnOn){
    Type T = Q[nspecies+3]*this->param->ref_temperature;
    //get dimensional rhoi
    for(i = 0; i < nspecies; i++){
      rhoi[i] = Q[i]*this->param->ref_density;
    }

    //wdots are returned in (kg/m^3.s)
    chem->GetMassProductionRates(rhoi, T, wdot);
    
    //non-dimensionalize output
    for(i = 0; i < nspecies; i++){
      wdot[i] /= (this->param->ref_density/this->param->ref_time);
      source[i] = vol*wdot[i];
    }
  }

  if(this->param->gravity_on){
    Type rho = 0.0;
    for(i = 0; i < nspecies; i++){
      rho += Q[i];
    }
    Type g = this->param->gravity;
    g /= (this->param->ref_length/(this->param->ref_time*this->param->ref_time));
    source[nspecies+0] += g*(rho-1.0)*vol*this->param->gravdir[0];
    source[nspecies+1] += g*(rho-1.0)*vol*this->param->gravdir[1];
    source[nspecies+2] += g*(rho-1.0)*vol*this->param->gravdir[2];
  }

  return;
}

//some eqnsets will have dQ/dq terms if they are in non-conservative variables
//this is cnp1*vol/dt*M where M is dq/dQ - this is one of those eqnsets
//dt - real time step used for residuals
//dtau - pseudotimestep used for convergence enhancement only
template <class Type>
void CompressibleFREqnSet<Type>::ContributeTemporalTerms(Type* Q, Type vol, Type cnp1, 
							 Type dt, Type dtau, Type* A, Type beta)
{
  Int i, j;
  Int neqn = this->neqn;
  Type* dEtdRhoi = new Type[nspecies];
  Type vOverDt;
  if(this->param->useLocalTimeStepping){
    vOverDt = cnp1*vol/dt + vol/dtau;
  }
  else{
    vOverDt = cnp1*vol/dtau;
  }
  
  Int uloc = nspecies;
  Int vloc = nspecies+1;
  Int wloc = nspecies+2;
  Int tloc = nspecies+3;
  Type rho = Q[nspecies+5];
  Type T = Q[tloc];
  Type P = Q[nspecies+4];
  Type u = Q[uloc];
  Type v = Q[vloc];
  Type w = Q[wloc];
  Type rvOverDt = rho*vOverDt;
  Type v2 = u*u + v*v + w*w;
  Type* rhoiDim = new Type[nspecies];
  Type* Yi = new Type[nspecies];
  Type* c2i = new Type[nspecies];
  Type TDim, PDim;
  Type precond;

  Type* cvi = &Q[nspecies+6];
  Type* rhoi = &Q[0];
  Type cv, cp, R, gamma, c2, RhoTrash, PTrash;
  GetFluidProperties(rhoi, T, cvi, cv, cp, R, gamma, c2, RhoTrash, PTrash);
  Type c = sqrt(c2);

  Type s_ref = (this->param->ref_velocity * this->param->ref_velocity / 
		this->param->ref_temperature);

  //dimensionalize
  Type rhoDim = rho*this->param->ref_density;
  TDim = T*this->param->ref_temperature;
  PDim = P*this->param->ref_pressure;
  v2 *= this->param->ref_velocity*this->param->ref_velocity;

  #pragma ivdep
  for(i = 0; i < nspecies; i++){
    rhoiDim[i] = rhoi[i]*this->param->ref_density;
    Yi[i] = rhoi[i]/rho;
  }

  Type dEtdP = chem->dEtdP_dEtdRhoi(rhoiDim, TDim, PDim, v2, dEtdRhoi);

  //non-dimensionalize
  Type ref_detdrho = this->param->ref_enthalpy*this->param->ref_density/this->param->ref_density;
  #pragma ivdep
  for(i = 0; i < nspecies; i++){
    dEtdRhoi[i] /= ref_detdrho;
  }

  Type ref_detdP = this->param->ref_enthalpy*this->param->ref_density/this->param->ref_pressure;
  dEtdP /= ref_detdP;
  //TODO: EOSUpdate
  Type dPdT = chem->eos[0]->GetdP_dT(R*s_ref, rhoDim, PDim, TDim);
  dPdT /= (this->param->ref_pressure/this->param->ref_temperature);

  //Now, compute the value for our limiting parameter thetaprime & beta - Ashish Gupta dissertation
  //This not the same theta prime used above for defining modified eigenvalues Eriksson preconditioner
  //TODO: is it required that we reset this for stability reasons... seems c2i should be individual not bulk
  GetSpeciesSpeedOfSound(c2i, Q);
  for(Int i = 0; i < nspecies; i++){
    c2i[i] = c2;
  }

  Type oneOBeta = 1.0/beta;
  Type* thetaOBetai = new Type[nspecies];
  Type oneMbeta = 1.0 - beta;
  #pragma ivdep
  for(Int i = 0; i < nspecies; i++){
    thetaOBetai[i] = (Yi[i]*oneMbeta/c2i[i]) * oneOBeta;
  }

  //first nspecies rows
  #pragma ivdep
  for(i = 0; i < nspecies; i++){
    A[neqn*i + i] += vOverDt;
    A[neqn*i + tloc] += thetaOBetai[i]*vOverDt*dPdT;
  }
  
  //x momentum
  precond = 0.0;
  for(j = 0; j < nspecies; j++){
    A[neqn*uloc + j] += u*vOverDt;
    precond += thetaOBetai[j]*u;
  }
  A[neqn*uloc + uloc] += rvOverDt;
  A[neqn*uloc + tloc] += precond*vOverDt*dPdT;

  //y momemtum
  precond = 0.0;
  for(j = 0; j < nspecies; j++){
    A[neqn*vloc + j] += v*vOverDt;
    precond += thetaOBetai[j]*v;
  }
  A[neqn*vloc + vloc] += rvOverDt;
  A[neqn*vloc + tloc] += precond*vOverDt*dPdT;

  //z momentum
  precond = 0.0;
  for(j = 0; j < nspecies; j++){
    A[neqn*wloc + j] += w*vOverDt;
    precond += thetaOBetai[j]*w;
  }
  A[neqn*wloc + wloc] += rvOverDt;
  A[neqn*wloc + tloc] += precond*vOverDt*dPdT;

  //temperature row
  precond = 0.0;
  for(j = 0; j < nspecies; j++){
    A[neqn*tloc + j] += dEtdRhoi[j]*vOverDt;
    precond += dEtdRhoi[j]*thetaOBetai[j];
  }
  A[neqn*tloc + uloc] += u*rvOverDt;
  A[neqn*tloc + vloc] += v*rvOverDt;
  A[neqn*tloc + wloc] += w*rvOverDt;

  precond += dEtdP*(oneOBeta);

  A[neqn*tloc + tloc] += precond*vOverDt*dPdT;

#if 0
  std::cout << std::endl;
  std::cout << "v/dt: " << vOverDt << std::endl;
  MatPrint(A, neqn);
  std::cout << std::endl;
#endif

  delete [] dEtdRhoi;
  delete [] rhoiDim;
  delete [] Yi;
  delete [] c2i;
  delete [] thetaOBetai;
}

template <class Type>
void CompressibleFREqnSet<Type>::UpdateQinf()
{
  //non-dimensionalization here is by p_inf, u_in, T_inf plus massfractions

  //use GetVelocity() function for ramping
  Type V = this->param->GetVelocity(this->space->iter);
  Type u = this->param->flowdir[0]*V;
  Type v = this->param->flowdir[1]*V;
  Type w = this->param->flowdir[2]*V;

  //set class variables to default initialization for Qinf
  if(this->param->massfractions.size() == nspecies){
    for(Int i = 0; i < nspecies; i++){
      this->Qinf[i] = this->param->massfractions.at(i);
    }
  }
  else{
    //give them all equal weight
    for(Int i = 0; i < nspecies; i++){
      this->Qinf[i] = 1.0/(Type)nspecies;
    }
  }

  //we assume the pressure non-dimensionalization 
  //is set to be equal to state conditions
  Type pDim = 1.0*this->param->ref_pressure;

  //assume that the reference temperature we give is the 
  //farfield temperature, this is a reasonable assumption
  Type T = 1.0;

  //Densities are additive in a fixed volume, sum them up
  Type r = 0;
  for(Int i = 0; i < nspecies; ++i){
    Type Xi = this->Qinf[i];
    Type piDim = pDim*Xi; //for ideal gases, Dalton's law applies, TODO: check for mixtures of phases
    r += chem->eos[i]->GetRho(chem->species[i].R, piDim, T*this->param->ref_temperature);
  }
  r /= this->param->ref_density;
  
  //scale mass fractions by total density, rhoi = rho*Xi
  for(Int i = 0; i < nspecies; i++){
    this->Qinf[i] *= r;
  }

  this->Qinf[nspecies] = u;
  this->Qinf[nspecies+1] = v;
  this->Qinf[nspecies+2] = w;
  this->Qinf[nspecies+3] = T;
  
  ComputeAuxiliaryVariables(this->Qinf);
  this->Pref = GetPressure(this->Qinf);
}

// returns fluid properties given rhoi's and temperature
// all variables passed back non-dimensional
template <class Type>
void CompressibleFREqnSet<Type>::GetFluidProperties(const Type* rhoi, const Type T, const Type* cvi, Type& cv, 
						    Type& cp, Type& R, Type& gamma, Type& c2, Type& rho, Type& P)  const
{
  Int i;
  Type* rhoiDim = (Type*)alloca(sizeof(Type)*nspecies);
  Type* cviDim = (Type*)alloca(sizeof(Type)*nspecies);
  Type cvDim = 0.0;
  Type cpDim = 0.0;
  Type RDim = 0.0;
  Type PDim = 0.0;
  Type rhoDim = 0.0;
  Type c2Dim = 0.0;

  //initialize passed in outputs
  rho = 0.0;
  P = 0.0;
  gamma = 0.0;
  R = 0.0;
  cv = 0.0;
  cp = 0.0;
  c2 = 0.0;
  
  // compute reference values for dimensionalization
  Type s_ref = (this->param->ref_velocity * this->param->ref_velocity / 
		this->param->ref_temperature);
  Type rho_ref = this->param->ref_density;
  Type p_ref = this->param->ref_pressure;

  //dimensionalize rhoi and cvi and T
  #pragma ivdep
  for(i = 0; i < nspecies; ++i){
    rhoiDim[i] = rhoi[i]*this->param->ref_density;
    cviDim[i] = cvi[i]*s_ref;
  }
  Type Tdim = T*this->param->ref_temperature;
  
  //make chemistry call
  this->chem->GetFluidProperties(rhoiDim, Tdim, cviDim, cvDim, cpDim, RDim, PDim, rhoDim, c2Dim);
  gamma = (cpDim)/(cvDim);
  cv = cvDim/s_ref;
  cp = cpDim/s_ref;
  R = RDim/s_ref;
  P = PDim/p_ref;
  rho = rhoDim/rho_ref;
  //we have to nondimensionalize speed of sound squared by the same reference velocity to be internally consistent
  c2 = c2Dim/(this->param->ref_velocity*this->param->ref_velocity);
}

template <class Type>
Type CompressibleFREqnSet<Type>::ComputeViscosity(Type* Q)
{
  Type T = GetTemperature(Q);
  Type* rhoi = &Q[0];
  return GetMolecularViscosity(rhoi, T);
}

template <class Type>
Type CompressibleFREqnSet<Type>::GetMolecularViscosity(Type* rhoi, Type T)
{
  Type* rhoidim = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; i++){
    rhoidim[i] = rhoi[i] * this->param->ref_density;
  }
  return chem->GetViscosity(rhoidim, T*this->param->ref_temperature) / this->param->ref_viscosity;
}

template <class Type>
Type CompressibleFREqnSet<Type>::GetThermalConductivity(Type* rhoi, Type T)
{
  Type* rhoidim = (Type*)alloca(sizeof(Type)*nspecies);
  
  for(Int i = 0; i < nspecies; i++){
    rhoidim[i] = rhoi[i] * this->param->ref_density;
  }
  
  return chem->GetThermalConductivity(rhoidim, T*this->param->ref_temperature) / this->param->ref_k;
}


template <class Type>
Type CompressibleFREqnSet<Type>::MaxEigenvalue(Type* Q, Type* avec, Type vdotn, Type gamma, Type beta)
{
  Type eig = 1.0;
  Type rho = Q[nspecies+5];
  Type T = Q[nspecies+3];
  Type* rhoi = &Q[0];
  Type* cvi = &Q[nspecies+6];
  Type theta;
  Type eig4, eig5;

  Type P = GetPressure(Q);

  Type cv, cp, R, gammaTrash, c2, RhoTrash, PTrash;
  GetFluidProperties(rhoi, T, cvi, cv, cp, R, gammaTrash, c2, RhoTrash, PTrash);
  Type c = sqrt(c2);

  theta = GetTheta(Q, avec, vdotn);

  //we use eigenvalues modified by the preconditioning matrix, in Ashish Gupta dissertation p.66
  Type oneMBeta = 1.0 - beta;
  Type thetaPrime = theta*(1.0 + beta)*0.5;
  Type cPrime = 0.5*sqrt(theta*theta*(oneMBeta*oneMBeta) + 4.0*beta*c2);
  eig4 = thetaPrime + cPrime;
  eig5 = thetaPrime - cPrime;

  eig = MAX(CAbs(eig4), CAbs(eig5));

#if 0
  std::cout << "MAXEIG: " << eig << std::endl;
  std::cout << "theta: " << thetaPrime << std::endl;
  std::cout << "c: " << cPrime << std::endl;
#endif

  return eig;
}

template <class Type>
Type CompressibleFREqnSet<Type>::GetTheta(Type* Q, Type* avec, Type vdotn)
{
  return (Q[nspecies]*avec[0] + Q[nspecies+1]*avec[1] + Q[nspecies+2]*avec[2] + vdotn);
}


template <class Type>
void CompressibleFREqnSet<Type>::GetInternalOutflowBoundaryVariables(Type* QL, Type* QR, Type pressure, Type gamma)
{
  for(Int i = 0; i < nspecies; i++){
    QR[i] = QL[i];
  }
  QR[nspecies] = QL[nspecies];
  QR[nspecies+1] = QL[nspecies+1];
  QR[nspecies+2] = QL[nspecies+2]; 

  Type Tguess = QL[nspecies+3];
  // set backpressure
  Type pgoal = pressure;
  
  // convert pressure solution back to temperature for storage, use Newton's method
  // to solve nonlinear mixture EOS
  Type* rhoi = &QR[0];
  Type T = NewtonFindTGivenP(rhoi, pgoal, Tguess);
  QR[nspecies+3] = T;
}

template <class Type>
void CompressibleFREqnSet<Type>::Dimensionalize(Type* Q)
{
  Int i;
  for(i = 0; i < nspecies; i++){
    Q[i] *= this->param->ref_density;
  }
  Q[nspecies] *= this->param->ref_velocity;
  Q[nspecies+1] *= this->param->ref_velocity;
  Q[nspecies+2] *= this->param->ref_velocity;
  Q[nspecies+3] *= this->param->ref_temperature;
  Q[nspecies+4] *= this->param->ref_pressure;
  Q[nspecies+5] *= this->param->ref_density;
}

template <class Type>
void CompressibleFREqnSet<Type>::ViscousJacobian(Type* QL, Type* QR, Type* dx, Type s2, Type* avec, 
						 Type mut, Type* aL, Type* aR)
{
  Int i, j;
  Int neqn = this->neqn;
  Int nvars = neqn + this->nauxvars;
  
  Type TL = this->GetTemperature(QL);
  Type TR = this->GetTemperature(QR);

  Type* Q = new Type[nvars];
  for(i = 0; i < neqn; i++){
    Q[i] = 0.5*(QL[i] + QR[i]);
  }
  ComputeAuxiliaryVariables(Q);
  Type T = GetTemperature(Q);

  Type* rhoi = &Q[0];
  Type rho = Q[nspecies+5];
  Type* cvi = &Q[nspecies+6];
  Type P = GetPressure(Q);

  Type mu = GetMolecularViscosity(rhoi, T);
  Type tmut = (mu + mut);

  Type* rhoiLdim = new Type[nspecies];
  Type* rhoiRdim = new Type[nspecies];
  for(i = 0; i < nspecies; i++){
    rhoiLdim[i] = QL[i] * this->param->ref_density;
    rhoiRdim[i] = QR[i] * this->param->ref_density;
  }

  Type ReTilde = this->param->Re;
  Type RK = avec[3]/ReTilde;
  Type RKT = RK*tmut;

  //compute cp for the mixture
  Type cv, cp, R, gamma, c2, RhoTrash, PTrash;
  GetFluidProperties(rhoi, T, cvi, cv, cp, R, gamma, c2, RhoTrash, PTrash);

  //compute some partial derivatives we need to continue
  //-------------------------------------------------------------------------
  
  //first row, we need the partial of each velocity derivative Ux, Vy, Wz, ... 
  //w.r.t density

  //build some of the coefficients into the distance variables
  Type* DxL = (Type*)alloca(sizeof(Type)*3);
  Type* DxR = (Type*)alloca(sizeof(Type)*3);
  Type* DxLrho = (Type*)alloca(sizeof(Type)*3);
  Type* DxRrho = (Type*)alloca(sizeof(Type)*3);

  Type rhoL = GetDensity(QL);
  Type* rhoiL = &QL[0];
  Type uL = QL[nspecies];
  Type vL = QL[nspecies+1];
  Type wL = QL[nspecies+2];
  Type PL = GetPressure(QL);
  Type* cviL = &QL[nspecies+6];
   
  Type rhoR = GetDensity(QR); 
  Type* rhoiR = &QR[0];
  Type uR = QR[nspecies];
  Type vR = QR[nspecies+1];
  Type wR = QR[nspecies+2];
  Type PR = GetPressure(QR);
  Type* cviR = &QR[nspecies+6];

  Type cvL, cpL, RL, gammaL, c2L, RhoTrashL, PTrashL;
  GetFluidProperties(rhoiL, TL, cviL, cvL, cpL, RL, gammaL, c2L, RhoTrashL, PTrashL);
  Type cvR, cpR, RR, gammaR, c2R, RhoTrashR, PTrashR;
  GetFluidProperties(rhoiR, TR, cviR, cvR, cpR, RR, gammaR, c2R, RhoTrashR, PTrashR);

  Type u = 0.5*(uL + uR);
  Type v = 0.5*(vL + vR);
  Type w = 0.5*(wL + wR);

  for(i = 0; i < 3; i++){
    //these are the directional derivative pieces
    DxL[i] = -dx[i]/s2;
    DxR[i] = dx[i]/s2;
    DxLrho[i] = DxL[i]/rhoL;
    DxRrho[i] = DxR[i]/rhoR;
  }

  Type dux_drL = -u*DxLrho[0];
  Type duy_drL = -u*DxLrho[1];
  Type duz_drL = -u*DxLrho[2];
  Type dvx_drL = -v*DxLrho[0];
  Type dvy_drL = -v*DxLrho[1];
  Type dvz_drL = -v*DxLrho[2];
  Type dwx_drL = -w*DxLrho[0];
  Type dwy_drL = -w*DxLrho[1];
  Type dwz_drL = -w*DxLrho[2];

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

  Type dux_drR = -u*DxRrho[0];
  Type duy_drR = -u*DxRrho[1];
  Type duz_drR = -u*DxRrho[2];
  Type dvx_drR = -v*DxRrho[0];
  Type dvy_drR = -v*DxRrho[1];
  Type dvz_drR = -v*DxRrho[2];
  Type dwx_drR = -w*DxRrho[0];
  Type dwy_drR = -w*DxRrho[1];
  Type dwz_drR = -w*DxRrho[2];

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

  //here dR2 refers to the x-momentum, dR3 y-momentum and dR4 z-momentum equations
  Type dR2_drhou, dR2_drhov, dR2_drhow;
  Type dR3_drhou, dR3_drhov, dR3_drhow;
  Type dR4_drhou, dR4_drhov, dR4_drhow;

  Type* row = NULL;

  //first nspecies rows   (dR1/dQr)
  for(i = 0; i < nspecies; i++){
    row = &aR[i*neqn];
    for(j = 0; j < neqn; j++){
      row[j] = 0.0;
    }
  }

  //velocity row  (dR2-{tauxn}/dQr)
  row = &aR[nspecies*neqn];
  dR2_drhou = (c43*DxRrho[0]*avec[0] + DxRrho[1]*avec[1] + DxRrho[2]*avec[2]);
  dR2_drhov = (mc23*DxRrho[1]*avec[0] + DxRrho[0]*avec[1]);
  dR2_drhow = (mc23*DxRrho[2]*avec[0] + DxRrho[0]*avec[2]);
  for(i = 0; i < nspecies; i++){
    row[i] = -RKT*(dtauxn_drR);
  }
  row[nspecies+0] = -RKT*dR2_drhou;
  row[nspecies+1] = -RKT*dR2_drhov;
  row[nspecies+2] = -RKT*dR2_drhow;
  row[nspecies+3] = 0.0;

  //velocity+1 row   (dR3-{tauyn}/dQr)
  row = &aR[(nspecies+1)*neqn];
  dR3_drhou = (mc23*DxRrho[0]*avec[1] + DxRrho[1]*avec[0]);
  dR3_drhov = (DxRrho[0]*avec[0] + c43*DxRrho[1]*avec[1] + DxRrho[2]*avec[2]);
  dR3_drhow = (mc23*DxRrho[2]*avec[1] + DxRrho[1]*avec[2]);
  for(i = 0; i < nspecies; i++){
    row[i] = -RKT*(dtauyn_drR);
  }
  row[nspecies+0] = -RKT*dR3_drhou;
  row[nspecies+1] = -RKT*dR3_drhov;
  row[nspecies+2] = -RKT*dR3_drhow;
  row[nspecies+3] = 0.0;

  //velocity+2 row  (dR4-{tauzn}/dQr)
  row = &aR[(nspecies+2)*neqn];
  dR4_drhou = (mc23*DxRrho[0]*avec[2] + DxRrho[2]*avec[0]);
  dR4_drhov = (mc23*DxRrho[1]*avec[2] + DxRrho[2]*avec[1]);
  dR4_drhow = (DxRrho[0]*avec[0] + DxRrho[1]*avec[1] + c43*DxRrho[2]*avec[2]);
  for(i = 0; i < nspecies; i++){
    row[i] = -RKT*(dtauzn_drR);
  }
  row[nspecies+0] = -RKT*dR4_drhou;
  row[nspecies+1] = -RKT*dR4_drhov;
  row[nspecies+2] = -RKT*dR4_drhow;
  row[nspecies+3] = 0.0;

  Type gm1 = (gamma - 1.0);

  Type v2L = (uL*uL + vL*vL + wL*wL);
  Type v2R = (uR*uR + vR*vR + wR*wR);

  //TODO: EOSUpdate
  Type dT_dPL = chem->eos[0]->GetdT_dP(RL, rhoL, PL, TL);
  Type dT_dPR = chem->eos[0]->GetdT_dP(RR, rhoR, PR, TR);

  Type* dTdRhoiL = new Type[nspecies];
  Type* dTdRhoiR = new Type[nspecies];

  chem->dTdRhoi(rhoiLdim, PL*this->param->ref_pressure, dTdRhoiL);
  chem->dTdRhoi(rhoiRdim, PR*this->param->ref_pressure, dTdRhoiR);

  //compute pressure derivatives
  Type dP_druL  = RL/cvL*uL;
  Type dP_drvL  = RL/cvL*vL;
  Type dP_drwL  = RL/cvL*wL;
  Type dP_dretL = RL/cvL;

  Type dP_druR  = RR/cvR*uR;
  Type dP_drvR  = RR/cvR*vR;
  Type dP_drwR  = RR/cvR*wR;
  Type dP_dretR = RR/cvR;

  //thermal conductivity
  Type k = GetThermalConductivity(rhoi, T);
  //turbulent thermal conductivity
  Type kT = mut/this->param->PrT*cp;

  //For now, assume k & kT = const. actually k = k(T) & kT = const.
  Type c1 = -(k + kT);

  //directional derivative part of temp grad
  Type TnL = (DxL[0]*avec[0] + DxL[1]*avec[1] + DxL[2]*avec[2])*c1;
  Type TnR = (DxR[0]*avec[0] + DxR[1]*avec[1] + DxR[2]*avec[2])*c1;

  //pressure row   (dR5/dQr)
  row = &aR[(nspecies+3)*neqn];
  for(i = 0; i < nspecies; i++){
    Type dT_drhoR = dTdRhoiR[i] / (this->param->ref_temperature / this->param->ref_density);
    row[i] = -RKT*(dtauxn_drR*u + dtauyn_drR*v + dtauzn_drR*w) + RK*TnR*dT_drhoR;
  }
  row[nspecies+0] = -RKT*(dR2_drhou*u + dR3_drhou*v + dR4_drhou*w) + RK*TnR*dT_dPR*dP_druR;
  row[nspecies+1] = -RKT*(dR2_drhov*u + dR3_drhov*v + dR4_drhov*w) + RK*TnR*dT_dPR*dP_drvR;
  row[nspecies+2] = -RKT*(dR2_drhow*u + dR3_drhow*v + dR4_drhow*w) + RK*TnR*dT_dPR*dP_drwR;
  row[nspecies+3] =                                                  RK*TnR*dT_dPR*dP_dretR;


  //Now, we do the same thing for the LHS jacobian
  //---------------------------------------------------------------------

  //first nspecies rows   (dR1/dQl)
  for(i = 0; i < nspecies; i++){
    row = &aL[i*neqn];
    for(j = 0; j < neqn; j++){
      row[j] = 0.0;
    }
  }

  //velocity row  (dR2-{tauxn}/dQl)
  row = &aL[nspecies*neqn];
  dR2_drhou = (c43*DxLrho[0]*avec[0] + DxLrho[1]*avec[1] + DxLrho[2]*avec[2]);
  dR2_drhov = (mc23*DxLrho[1]*avec[0] + DxLrho[0]*avec[1]);
  dR2_drhow = (mc23*DxLrho[2]*avec[0] + DxLrho[0]*avec[2]);
  for(i = 0; i < nspecies; i++){
    row[i] = -RKT*(dtauxn_drL);
  }
  row[nspecies+0] = -RKT*dR2_drhou;
  row[nspecies+1] = -RKT*dR2_drhov;
  row[nspecies+2] = -RKT*dR2_drhow;
  row[nspecies+3] = 0.0;

  //velocity+1 row   (dR3-{tauyn}/dQl)
  row = &aL[(nspecies+1)*neqn];
  dR3_drhou = (mc23*DxLrho[0]*avec[1] + DxLrho[1]*avec[0]);
  dR3_drhov = (DxLrho[0]*avec[0] + c43*DxLrho[1]*avec[1] + DxLrho[2]*avec[2]);
  dR3_drhow = (mc23*DxLrho[2]*avec[1] + DxLrho[1]*avec[2]);
  for(i = 0; i < nspecies; i++){
    row[i] = -RKT*(dtauyn_drL);
  }
  row[nspecies+0] = -RKT*dR3_drhou;
  row[nspecies+1] = -RKT*dR3_drhov;
  row[nspecies+2] = -RKT*dR3_drhow;
  row[nspecies+3] = 0.0;

  //velocity+2 row  (dR4-{tauzn}/dQl)
  row = &aL[(nspecies+2)*neqn];
  dR4_drhou = (mc23*DxLrho[0]*avec[2] + DxLrho[2]*avec[0]);
  dR4_drhov = (mc23*DxLrho[1]*avec[2] + DxLrho[2]*avec[1]);
  dR4_drhow = (DxLrho[0]*avec[0] + DxLrho[1]*avec[1] + c43*DxLrho[2]*avec[2]);
  for(i = 0; i < nspecies; i++){
    row[i] = -RKT*(dtauzn_drL); 
  }
  row[nspecies+0] = -RKT*dR4_drhou;
  row[nspecies+1] = -RKT*dR4_drhov;
  row[nspecies+2] = -RKT*dR4_drhow;
  row[nspecies+3] = 0.0;

  //pressure row   (dR5/dQr)
  row = &aL[(nspecies+3)*neqn];
  for(i = 0; i < nspecies; i++){
    Type dT_drhoL = dTdRhoiL[i] / (this->param->ref_temperature / this->param->ref_density);
    row[i] = -RKT*(dtauxn_drL*u + dtauyn_drL*v + dtauzn_drL*w) + RK*TnL*dT_drhoL;
  }
  row[nspecies+0] = -RKT*(dR2_drhou*u + dR3_drhou*v + dR4_drhou*w) + RK*TnL*dT_dPL*dP_druL;
  row[nspecies+1] = -RKT*(dR2_drhov*u + dR3_drhov*v + dR4_drhov*w) + RK*TnL*dT_dPL*dP_drvL;
  row[nspecies+2] = -RKT*(dR2_drhow*u + dR3_drhow*v + dR4_drhow*w) + RK*TnL*dT_dPL*dP_drwL;
  row[nspecies+3] =                                                  RK*TnL*dT_dPL*dP_dretL;


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
	aL[nspecies*neqn + i] = aR[nspecies*neqn + i] = 0.0;
      }
    }
    else if(this->param->symmetry2D == 2){
      for(i = 0; i < neqn; i++){
	aL[(nspecies+1)*neqn + i] = aR[(nspecies+1)*neqn + i] = 0.0;
      }
    }
    else if(this->param->symmetry2D == 3){
      for(i = 0; i < neqn; i++){
	aL[(nspecies+2)*neqn + i] = aR[(nspecies+2)*neqn + i] = 0.0;
      }
    }
  }


  delete [] Q;
  delete [] rhoiLdim;
  delete [] rhoiRdim;
  delete [] dTdRhoiL;
  delete [] dTdRhoiR;

  return;
}

template <class Type>
void CompressibleFREqnSet<Type>::GetViscousWallBoundaryVariables(Type* QL, Type* QR, Type* vel, Int normalNode, Type* normalQ, Type Twall)
{
  //adiabatic
  if(real(Twall) < 0.0){
    for(Int i = 0; i < nspecies; i++){
      QR[i] = QL[i];
    }
    QR[nspecies+3] = QL[nspecies+3] = normalQ[nspecies+3];
  }
  //constant temperature
  else{
    //constant wall temp implies density imposed from the field
    //and total energy set from the state equation
    for(Int i = 0; i < nspecies; i++){
      QR[i] = QL[i];
    }
    QR[nspecies+3] = QL[nspecies+3] = Twall; 
  }
  QR[nspecies] = QL[nspecies] = vel[0];
  QR[nspecies+1] = QL[nspecies+1] = vel[1];
  QR[nspecies+2] = QL[nspecies+2] = vel[2];
}

template <class Type>
void CompressibleFREqnSet<Type>::ModifyViscousWallJacobian(Type* QL, Type* QR, 
							   Type* vel, Int cvid, CRS<Type>* crs,
							   Int normalNode, Type Twall)
{
  //blank all the subrows for ru, rv, rw components
  crs->A->BlankSubRow(cvid, nspecies);
  crs->A->BlankSubRow(cvid, nspecies+1);
  crs->A->BlankSubRow(cvid, nspecies+2);

  //modify diagonal matrix as described in Anderson's Compute Fluids Vol. 23 paper
  //"An Implicit Upwind Algorithm for Computing Turbulent Flows on Unstructured Grids"
  Type* diag = crs->A->GetPointer(cvid, cvid);
  Type gamma = this->param->gamma;
  
  //adiabatic
  if(real(Twall) < 0.0){
    Type* offJac = crs->A->GetPointer(cvid, normalNode);
    crs->A->BlankSubRow(cvid, nspecies+3);
    //set temperature dependence on normal node
    offJac[(nspecies+3)*this->neqn + nspecies+3] = -1.0;
  }
  //constant wall temperature
  else{
    crs->A->BlankSubRow(cvid, nspecies+3);
  }

  return;
}

template <class Type>
void CompressibleFREqnSet<Type>::ModifyViscousWallResidual(Type* res, Type* vel, Int normalNode, 
							   Type Twall)
{
  //zero the rows we are setting explicitly, namely ru,rv,rw
  //this has the effect of not allowing du, dv, dw to change during the
  //linear system solve
  res[nspecies] = 0.0;
  res[nspecies+1] = 0.0;
  res[nspecies+2] = 0.0;

  //this is done in Anderson's paper, d(rho*Et) ~ drho
  res[nspecies+3] = 0.0;

  return;
}


template <class Type>
void CompressibleFREqnSet<Type>::NativeToConservative(Type* Q)
{
  Type rho = GetDensity(Q);
  Type Et = GetTotalEnergy(Q);

  //compute ru, rv, rw
  Q[nspecies] *= rho;
  Q[nspecies+1] *= rho;
  Q[nspecies+2] *= rho;
  //replace the pressure term with Et
  Q[nspecies+3] = Et;

  return;
}

// Takes a vector of form [rhoi, rhou, rhov, rhow, rhoEt] and converts it back to
// [rhoi, u, v, w, T] for convenience as stored in the solver
template <class Type>
void CompressibleFREqnSet<Type>::ConservativeToNative(Type* Q)
{
  Int maxit = 20;
  Type tol = 1.0e-12;
  Type rho = 0.0;
  Type* rhoi = &Q[0];
  Type* hi = (Type*)alloca(sizeof(Type)*nspecies);
  Type* Y = (Type*)alloca(sizeof(Type)*nspecies);
  Type* rhoiDim = (Type*)alloca(sizeof(Type)*nspecies);
  Type R = 0.0;
  for(Int i = 0; i < nspecies; i++){
    rho += rhoi[i];
    rhoiDim[i] = rhoi[i]*this->param->ref_density;
    R += rhoiDim[i]*chem->species[i].R;
  }
  R /= this->param->ref_density*rho;
  for(Int i = 0; i < nspecies; ++i){
    Y[i] = rhoi[i]/rho;
  }
  Type u = Q[nspecies] / rho;
  Type v = Q[nspecies+1] / rho;
  Type w = Q[nspecies+2] / rho;
  Type v2 = u*u + v*v + w*w;
  // Et -  0.5*rho*v^2
  // look only at the static (internal) energy minus the kinematic part
  // we store the total energy including kinematic
  Type res = Q[nspecies+3] - 0.5*v2*rho;

  //use last pressure stored to guess at T
  Type P = GetPressure(Q);
  //TODO: EOSUpdate
  Type T = chem->eos[0]->GetT(R, rho*this->param->ref_density, P*this->param->ref_pressure)/
    this->param->ref_temperature;
  Type To = 1.0;
  Int j = 0;
  Type dT = 0.0;
  for(j = 0; j < maxit; j++){
    Type Tp = T + 1.0e-8;
    Type H = rho*(chem->GetSpecificEnthalpy(Y, T*this->param->ref_temperature, hi)/
		  this->param->ref_enthalpy);
    Type Hp = rho*(chem->GetSpecificEnthalpy(Y, Tp*this->param->ref_temperature, hi)/
		   this->param->ref_enthalpy);
    Type P = chem->GetP(rhoiDim, T*(this->param->ref_temperature))/this->param->ref_pressure;
    Type Pp = chem->GetP(rhoiDim, Tp*(this->param->ref_temperature))/this->param->ref_pressure;
    Type E = H - P;
    Type Ep = Hp - Pp;
    Type zpoint = res - E;
    Type zpointp = res - Ep;
    Type dzdT = (zpointp - zpoint)/(Tp - T);
    dT = -zpoint/dzdT;
    T += dT;
    if (real(CAbs(dT)) < real(tol)) break;
  }

  if(j == maxit){
    std::stringstream ss;
    ss << "WARNING: Newton iteration did not converge on a temperature in ConservativeToNative()" 
	      << std::endl;
    ss << "Last dT = " << dT << std::endl;
    Abort << ss.str();
  }
  
  //compute u, v, w
  Q[nspecies] = u;
  Q[nspecies+1] = v;
  Q[nspecies+2] = w;
  //replace Et with temperature
  Q[nspecies+3] = T;
}

template <class Type>
void CompressibleFREqnSet<Type>::ComputeStressVector(Type* grad, Type* avec, Type mu, Type* stress)
{
  //gradient is passed in with velocity components
  //parse for readability
  Type ux, uy, uz, vx, vy, vz, wx, wy, wz;
  Int offset = nspecies*3;
  ux = grad[offset + 0];
  uy = grad[offset + 1];
  uz = grad[offset + 2];
  vx = grad[offset + 3];
  vy = grad[offset + 4];
  vz = grad[offset + 5];
  wx = grad[offset + 6];
  wy = grad[offset + 7];
  wz = grad[offset + 8];

  Type tauxx, tauyy, tauzz, tauxy, tauxz, tauyz;
  Type div = -2.0/3.0*(ux + vy + wz);
  Type tauxn, tauyn, tauzn;
  tauxx = 2.0*ux + div;
  tauyy = 2.0*vy + div;
  tauzz = 2.0*wz + div;
  tauxy = uy + vx;
  tauxz = uz + wx;
  tauyz = vz + wy;

  Type Re = this->param->Re;

  //written as on RHS of the equation
  tauxn = tauxx*avec[0] + tauxy*avec[1] + tauxz*avec[2];
  tauyn = tauxy*avec[0] + tauyy*avec[1] + tauyz*avec[2];
  tauzn = tauxz*avec[0] + tauyz*avec[1] + tauzz*avec[2];

  //negative sign indicates stress on the body, not the flow
  stress[0] = -(mu/Re)*tauxn;
  stress[1] = -(mu/Re)*tauyn;
  stress[2] = -(mu/Re)*tauzn;

  return;

}

template <class Type>
Type CompressibleFREqnSet<Type>::GetCf(Type tauw, Type rho)
{
  //Is this correct?  Implementation from compressible perfect gas solver.
  Type V = this->param->GetVelocity(this->space->iter);
  return (tauw/(0.5*rho*V*V));
}

template <class Type>
void CompressibleFREqnSet<Type>::GetTotalTempPressureBoundaryVariables(Type* QL, Type* QR, 
								       Type Pressure, Type T, 
								       Type* flowDirection)
{
  Type TDim = T*this->param->ref_temperature;
  Type pDim = Pressure*this->param->ref_pressure;
  
  //solve for rho*R
  //TODO: EOSUpdate
  Type rhoR = chem->eos[0]->GetRhoR(pDim, TDim);

  //get massfractions
  Type* X = (Type*)alloca(sizeof(Type)*nspecies);
  if(this->param->massfractions.size() == nspecies){
    for(Int i = 0; i < nspecies; i++){
      X[i] = this->param->massfractions.at(i);
    }
  }
  else{
    //give them all equal weight
    Type rn = 1.0 / (Type)nspecies;
    for(Int i = 0; i < nspecies; i++){
      X[i] = rn;
    }
  }
  
  //compute R and rho using massfractions
  Type R = 0.0;
  for(Int i = 0; i < nspecies; i++){
    R += this->chem->species[i].R*X[i];
  }
  Type rhoDim = (rhoR/R);
  Type rho = rhoDim/this->param->ref_density;

  Type* rhoiDim = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; ++i){
    rhoiDim[i] = X[i]*rhoDim;
  }
  
  //compute cv
  Type cv = this->chem->GetCv(rhoiDim, TDim) / (this->param->ref_velocity*this->param->ref_velocity/this->param->ref_temperature);

  //Get specific energy
  Type* hi = (Type*)alloca(sizeof(Type)*nspecies);
  Type hs = chem->GetSpecificEnthalpy(X, TDim, hi);
 
  Type RoCv = R/cv;
  Type V2dim = (RoCv*hs*rhoDim)/(RoCv/2.0 * rhoDim);

  Type V2 = V2dim /(this->param->ref_velocity*this->param->ref_velocity);
  Type V = sqrt(V2);

  Type u = V*flowDirection[0];
  Type v = V*flowDirection[1];
  Type w = V*flowDirection[2];

  //load up vector
  for(Int i = 0; i < nspecies; i++){
    QR[i] = rho*X[i];
  }
  QR[nspecies+0] = u;
  QR[nspecies+1] = v;
  QR[nspecies+2] = w;
  QR[nspecies+3] = T;
}

template <class Type>
void CompressibleFREqnSet<Type>::GetSpeciesSpeedOfSound(Type* c2i, Type* Q)
{
  //by definition c^2 = dp_drho| constant entropy
  Type s_ref = (this->param->ref_velocity * this->param->ref_velocity / 
		this->param->ref_temperature);
  Type v2_ref = this->param->ref_velocity*this->param->ref_velocity;
  Type T = Q[nspecies+3];
  Type p = Q[nspecies+4];
  Type pdim = p*this->param->ref_pressure;
  Type Tdim = T*this->param->ref_temperature;
  Type* rhoi = &Q[0];
  Type* cvi = &Q[nspecies+6];
  for(Int i = 0; i < nspecies; i++){
    Type Ri = chem->species[i].R;
    //TODO: EOSUpdate
    Type cpi = chem->eos[i]->GetCp(cvi[i]*s_ref, Ri, rhoi[i]*this->param->ref_density, pdim, Tdim);
    Type gamma = cpi/(cvi[i]*s_ref);
    c2i[i] = gamma*Ri*Tdim;
    c2i[i] /= v2_ref;
  }
}

template <class Type>
Type CompressibleFREqnSet<Type>::NewtonFindTGivenP(const Type* rhoi, const Type Pgoal, const Type Tinit) const
{
  Int maxit = 28;
  Type tol = 1.0e-15;

  Type TDim = Tinit * this->param->ref_temperature;
  Type PgoalDim = Pgoal * this->param->ref_pressure;
  Type* rhoiDim = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; ++i){
    rhoiDim[i] = rhoi[i]*this->param->ref_density;
  }

  Int j = 0;
  Type dT = 0.0;
  // all dimensional values internal to the loop
  for(j = 0; j < maxit; j++){
    Type TpDim = TDim + 1.0e-8;
    Type PDim = this->chem->GetP(rhoiDim, TDim);
    Type PpDim = this->chem->GetP(rhoiDim, TpDim);
    Type zpoint = PgoalDim - PDim;
    Type zpointp = PgoalDim - PpDim;
    Type dzdT = (zpointp - zpoint)/(TpDim - TDim);
    dT = -zpoint/dzdT;
    if (real(CAbs(dT/this->param->ref_temperature)) < real(tol)) break;
    else TDim += dT;
  }

  if(j == maxit){
    std::stringstream ss;
    ss << "WARNING: Newton iteration did not converge on a temperature in NewtonFindTGivenP()" << std::endl;
    ss << "Last dT = " << dT << std::endl;
    std::cerr << ss.str() << std::endl;
    Abort << ss.str();
  }
  return TDim/this->param->ref_temperature;
}
