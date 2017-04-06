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
	Real* molfrac;

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
	molfrac = new Real[chem.nspecies];

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
	double gamma = 0.0;
  double cv = 0.0;
	double cp = 0.0;
  double X[chem.nspecies]; //massfraction
  if(param.massfractions.size() == chem.nspecies){
    int Tlevels = (int)3500/100.0;
    std::cout << std::endl;
    std::cout << "Temp-(K)\tCv-(J/kg.K)\tCp-(J/kg.K)\tCp/R\tmu-(Pa.s)\tk-(W/m.K)" << std::endl;
    std::cout << "----------------------------------------------------------------------------" << std::endl;
    for(j = 0; j < Tlevels; j++){
      cv = 0.0;
			cp = 0.0;
      double Ti = (double)(j*100.0 + 100.0);
      for(i = 0; i < chem.nspecies; i++){
				X[i] = rhoi[i]/rho;
				double cpi = chem.species[i].GetCp(Ti);
				double cvi = chem.eos->GetCv(cpi, R, rho, P, Ti);
				cp += param.massfractions[i]*cpi;
				cv += param.massfractions[i]*cvi;
      }
			double mu = chem.GetViscosity(rhoi, Ti);
			double k = chem.GetThermalConductivity(rhoi, Ti);
      std::cout << Ti << "\t" << cv << "\t" << cp << "\t" << cp/R << "\t" << mu << "\t" << k << std::endl;
    }
    std::cout << std::endl;
    cv = 0.0;
		cp = 0.0;
    for(i = 0; i < chem.nspecies; i++){
      X[i] = rhoi[i]/rho;
      double cpi = chem.species[i].GetCp(TGiven);
      double cvi = chem.eos->GetCv(cpi, R, rho, P, TGiven);
      std::cout << "cv[" << chem.species[i].symbol << "]: " 
								<<  cvi << " (J/kg.K)" << std::endl;
      cv += param.massfractions[i]*cvi;
			cp += param.massfractions[i]*cpi;
    }
  }
  else{
    std::cerr << "Number of species defined in param file does not match chem model" << std::endl;
    return(-1);
  }

	//compute mol fractions and MW_mix
	double MWmix = 0.0;
	std::cout << "\nMole fractions" << std::endl;
	std::cout << "========================= " << std::endl;
	double summ = 0.0;
	for(int i = 0; i < chem.nspecies; ++i){
		molfrac[i] = (rhoi[i]/rho)/chem.species[i].MW;
		summ += molfrac[i];
	}
	for(int i = 0; i < chem.nspecies; ++i){
		molfrac[i] /= summ;
		MWmix += molfrac[i] * chem.species[i].MW;
		std::cout << "xi[" << chem.species[i].symbol << "]: " << molfrac[i] << std::endl;
	}
	std::cout << "\nMass fractions" << std::endl;
	std::cout << "========================= " << std::endl;
	for(int i = 0; i < chem.nspecies; ++i){
		std::cout << "Yi[" << chem.species[i].symbol << "]: " << param.massfractions[i] << std::endl;
	}
	std::cout << std::endl;

	std::cout << "Mixture properties at " << TGiven << " (K)" << std::endl;
	std::cout << "=======================================" << std::endl;

  std::cout << "rho: " << rho << " kg/m^3" << std::endl;
  std::cout << "Rmix: " << R << " J/kg.K" << std::endl;
  std::cout << "Static pressure (EOS only): " << P << " Pa" << std::endl;
  std::cout << "cvmix: " << cv << " (J/kg.K)" << std::endl;
	std::cout << "cpmix: " << cp << " (J/kg.K)" << std::endl;
	std::cout << "mwmix: " << MWmix << " (kg/mol)" << std::endl;
  cp = chem.eos->GetCp(cv, R, rho, P, TGiven);
  gamma = cp/cv;
  std::cout << "gammamix: " << gamma << std::endl;
	std::cout << "Thermal conductivity: " << chem.GetThermalConductivity(rhoi, TGiven) << " (W/m.K)" << std::endl;
	std::cout << "Viscosity: " << chem.GetThermalConductivity(rhoi, TGiven) << " (Pa.s)" << std::endl;
  double c = sqrt(gamma*R*TGiven);
  std::cout << "c (speed of sound): " << c << " m/s" << std::endl;
  double hi[chem.nspecies];
	double u = param.flowdir[0]*param.GetVelocity(1)*param.ref_velocity;
  double v = param.flowdir[1]*param.GetVelocity(1)*param.ref_velocity;
  double w = param.flowdir[2]*param.GetVelocity(1)*param.ref_velocity;
  double v2 = u*u + v*v + w*w;
  std::cout << "U: " << u  << " m/s" << std::endl;
  std::cout << "V: " << v  << " m/s" << std::endl;
  std::cout << "W: " << w  << " m/s" << std::endl;
	std::cout << "Mach: " << sqrt(v2/(c*c)) << std::endl;
  double Ht = rho*chem.GetSpecificEnthalpy(X, TGiven, hi) + 0.5*rho*v2;
  double Et = rho*chem.GetSpecificEnthalpy(X, TGiven, hi) - P + 0.5*rho*v2;
  std::cout << "Total enthalpy: " << Ht/1000.0 << " (kJ)" << std::endl;
  std::cout << "Total energy: " << Et/1000.0 <<  " (kJ)" << std::endl;
  std::cout << "Total internal energy: " << (Et - 0.5*rho*v2)/1000.0 <<  " (kJ)" << std::endl;
  
  std::cout << "Total pressure (gamma-1.0 formula): " << ((gamma - 1.0)*(Et - 0.5*rho*v2)/param.ref_enthalpy)*
    param.ref_pressure << " Pa" << std::endl;
  std::cout << "Total temperature (gamma-1.0 formula): " << TGiven*(1.0  + (gamma-1.0)/2.0*(v2/(c*c))) << " (K) " << std::endl;
    
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

#if 0
  //temporal loop to compute change in makeup over time
  double dt = 0.0001;
  double volume = 1.0; //m^3
  Real* source = new Real[chem.nspecies];
  Real* Y = new Real[chem.nspecies];
  P = Pt;
  Real tol = 1.0e-12;
  Real T = TGiven;
  for(i = 0; i < 20; ++i){
    std::cout << dt*i << " ----------------------------------------------" << std::endl;
    std::cout << "Temp: " << T << std::endl;
    chem.GetMassProductionRates(rhoi, T, wdot);
    for(int is = 0; is < chem.nspecies; ++is){
      source[is] = wdot[is]*volume*dt;
      std::cout << "s: " << source[is] << std::endl;
    }

    //update the masses/densities given the source term (production/destruction)
    rho = 0.0;
    R = 0.0;
    for(int is = 0; is < chem.nspecies; ++is){
      Real mass = rhoi[is]*volume + source[is];
      rhoi[is] = mass/volume;
      rho += rhoi[is];
      R += rhoi[is]*chem.species[is].R;
    }
    std::cout << "Rho: " << rho << std::endl;
    // R is dimensional after this
    R /= rho;
    for(int is = 0; is < chem.nspecies; ++is){
      Y[is] = rhoi[is]/rho;
      std::cout << "y: " << Y[is] << std::endl;
    }
    
    int j = 0;
    Int maxit = 30;
    Real dT = 0.0;
    //todo: check if this is correct
    Real res = Et - 0.5*rho*v2;
    for(j = 0; j < maxit; j++){
      Real Tp = T + 1.0e-8;
      Real H = rho*(chem.GetSpecificEnthalpy(Y, T, hi));
      Real Hp = rho*(chem.GetSpecificEnthalpy(Y, Tp, hi));
      Real P = chem.eos->GetP(R, rho, T);
      Real Pp = chem.eos->GetP(R, rho, Tp);
      Real E = H - P;
      Real Ep = Hp - Pp;
      Real zpoint = res - E;
      Real zpointp = res - Ep;
      Real dzdT = (zpointp - zpoint)/(Tp - T);
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
#endif
  
  delete [] dEtdRhoi;
  delete [] rhoi;
  delete [] wdot;
	delete [] molfrac;

  return(ierr);
}
