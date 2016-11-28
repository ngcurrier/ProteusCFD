#include <iostream>
#include <sstream>
#include "chem.h"
#include "species.h"
#include "reaction.h"
#include "h5layer.h"
#include "param.h"
#include "general.h"
#include "temporalControl.h"
#include "solutionOperations.h"

int main(int argc, char* argv[])
{
  Int i, j;
  Int ierr = 0;
  std::stringstream ss;
  Real TGiven;
  Real* rhoi;
  Real* wdot;

  if(argc != 3){
    std::cerr << "USAGE: " << argv[0] << " casename Temperature(K)" << std::endl;
    return(1);
  }
  
  //setup parameter file so we can read reference states, etc.
  std::vector<Param<double>*> paramList;
  SolutionOrdering<Real> operations;
  TemporalControl<double> temporalControl;
  std::string casestring = argv[1];
  size_t pos = casestring.rfind('/');
  std::string pathname;
  if(pos != std::string::npos){
    pathname = casestring.substr(0, pos+1);
    casestring = casestring.substr(pos);
  }
  else{
    pathname = "./";
  }
  if(ReadParamFile(paramList, casestring, pathname)){
    return (-1);
  }
  if(ReadSolutionOrdering(operations, casestring, pathname)){
    return (-1);
  }
  if(ReadTemporalControl(temporalControl, casestring, pathname)){
    return (-1);
  }

  //only use the first solution space defined in param file
  Param<double> param = *paramList.front();

  //read temperature for production rates from command line
  ss << argv[2];
  ss >> TGiven;

  //read reaction file
  Bool viscous = true;
  ChemModel<double> chem(param);

  rhoi = new Real[chem.nspecies];
  wdot = new Real[chem.nspecies];

  //using massfraction information available from param file if there
  //print out standard state conditions using chemistry
  if(param.massfractions.size() == chem.nspecies){
    for(i = 0; i < chem.nspecies; i++){
      std::cout << "rho[" << chem.species[i].symbol << "]: " 
		<< param.massfractions[i]*param.ref_density << " kg/m^3" << std::endl;
    }
  }
  else{
    std::cerr << "Number of species defined in param file does not match chem model" << std::endl;
    return(-1);
  }

  double u = param.flowdir[0]*param.GetVelocity(1)*param.ref_velocity;
  double v = param.flowdir[1]*param.GetVelocity(1)*param.ref_velocity;
  double w = param.flowdir[2]*param.GetVelocity(1)*param.ref_velocity;

  std::cout << "U: " << u  << " m/s" << std::endl;
  std::cout << "V: " << v  << " m/s" << std::endl;
  std::cout << "W: " << w  << " m/s" << std::endl;
  
  //compute rho and R
  double rho = 0.0;
  double R = 0.0;
  if(param.massfractions.size() == chem.nspecies){
    for(i = 0; i < chem.nspecies; i++){
      rhoi[i] = param.massfractions[i]*param.ref_density;
      rho += rhoi[i];
      R += chem.species[i].R*(rhoi[i]);
    }
  }
  //R is dimensional after this
  R /= rho;
  double P = chem.eos->GetP(R, rho, TGiven);
  
  std::cout << "rho: " << rho << " kg/m^3" << std::endl;
  std::cout << "Rmix: " << R << std::endl;
  std::cout << "Pressure (EOS only): " << P << " Pa" << std::endl;
  double gamma = 0.0;
  double cv = 0.0;
  double X[chem.nspecies];
  if(param.massfractions.size() == chem.nspecies){
    for(i = 0; i < chem.nspecies; i++){
      X[i] = rhoi[i]/rho;
      double cpi = chem.species[i].GetCp(TGiven);
      double cvi = chem.eos->GetCv(cpi, R, rho, P, TGiven);
      std::cout << "cv[" << chem.species[i].symbol << "]: " 
		<<  cvi << std::endl;
      cv += param.massfractions[i]*cvi;
    }
  }
  else{
    std::cerr << "Number of species defined in param file does not match chem model" << std::endl;
    return(-1);
  }
  std::cout << "cvmix: " << cv << std::endl;
  double cp = chem.eos->GetCp(cv, R, rho, P, TGiven);
  gamma = cp/cv;
  std::cout << "gammamix: " << gamma << std::endl;
  std::cout << "c (speed of sound): " << sqrt(gamma*R*TGiven) << " m/s" << std::endl;
  double hi[chem.nspecies];
  double v2 = u*u + v*v + w*w;
  double Ht = rho*chem.GetSpecificEnthalpy(X, TGiven, hi) + 0.5*rho*v2;
  double Et = rho*chem.GetSpecificEnthalpy(X, TGiven, hi) - P + 0.5*rho*v2;
  std::cout << "Total enthalpy: " << Ht << std::endl;
  std::cout << "Total energy: " << Et << std::endl;

  std::cout << "Pressure (gamma-1.0 formula): " << ((gamma - 1.0)*(Et - 0.5*rho*v2)/param.ref_enthalpy)*
    param.ref_pressure << "Pa" << std::endl;

    
  std::cout << "\nAt given temperature of " << TGiven << "K production rates are: " << std::endl;
  std::cout << "===================================================" << std::endl;
  chem.GetMassProductionRates(rhoi, TGiven, wdot);
  for(i = 0; i < chem.nspecies; i++){
    std::cout << chem.species[i].symbol << ": " << wdot[i] << " kg/(m^3 s)" << std::endl;
  }
  
  Real sum = 0.0;
  for(i = 0; i < chem.nspecies; i++){
    sum += wdot[i];
  }
  std::cout << "Mass blanance: " << sum << " kg/(m^3 s)" << std::endl;

  std::cout << "\nDerivatives at given temp: " << TGiven << std::endl;
  std::cout << "===================================================" << std::endl;
  double* dEtdRhoi = new double[chem.nspecies];
  double Pt = chem.eos->GetP(R, rho, TGiven);
  double dEtdP = chem.dEtdP_dEtdRhoi(rhoi, TGiven, Pt, v2, dEtdRhoi);
  std::cout << "dEtdP: " << dEtdP << std::endl;
  for(Int i = 0; i < chem.nspecies; i++){
    std::cout << "dEtdrho[" << i << "]: " << dEtdRhoi[i] << std::endl;
  }

  //temporal loop to compute change in makeup over time
  double dt = 0.001;
  double volume = 1.0; //m^3
  Real* source = new Real[chem.nspecies];
  Real* Y = new Real[chem.nspecies];
  Real P = Pt;
  Type tol = 1.0e-12;
  Int maxit = 20;
  for(i = 0; i < 10; ++i){
    chem.GetMassProductionRates(rhoi, T, wdot);
    for(int is = 0; is < chem.nspecies; ++is){
      source[is] = wdot[is]*vol*dt;
    }

    //update the masses/densities given the source term (production/destruction)
    rho = 0.0;
    R = 0.0;
    for(int is = 0; is < chem.nspecies; ++is){
      Real mass = rhoi[is]*vol + source[is];
      rhoi[is] = mass/vol;
      rho += rhoi[is];
      R += rhoi[is]*chem.species[i].R;
    }
    // R is dimensional after this
    R /= rho;
    for(int is = 0; is < chem.nspecies; ++is){
      Y[is] = rhoi[is]/rho;
    }
    
    //use last good temperature to guess at T
    int j = 0;
    Real T = TGiven;
    for(j = 0; j < maxit; j++){
      Type Tp = T + 1.0e-8;
      Type H = rho*(chem->GetSpecificEnthalpy(Y, T, hi));
      Type Hp = rho*(chem->GetSpecificEnthalpy(Y, Tp, hi));
      Type P = chem->eos->GetP(R, rho, T);
      Type Pp = chem->eos->GetP(R, rho,	Tp);
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
      std::cerr << "WARNING: Newton iteration did not converge on a temperature in ConservativeToNative()" 
		<< std::endl;
      std::cerr << "Last dT = " << dT << std::endl;
    }
    
  }
  delete [] source;
  delete [] Y;
  
  delete [] dEtdRhoi;
  delete [] rhoi;
  delete [] wdot;

  return(ierr);
}
