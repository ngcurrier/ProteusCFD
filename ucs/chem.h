#ifndef CHEM_H__
#define CHEM_H__

#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include "general.h"
#include "species.h"
#include "reaction.h"
#include "chem_constants.h"
#include "param.h"
#include "EOS.h"

//NOTE ABOUT STANDARD UNITS:
//=========================
//Temperature - Kelvin
//Mass - kilogram
//Time - seconds
//Length - meters
//Concentration - mole/m^3
//Energy - Joule
//Force - Newton

//forward declarations
template <class Type> class Reaction;
template <class Type> class Species;

//Note: all inputs are expect dimensional in SI, all outputs are dimensional in SI
template <class Type>
class ChemModel
{
 public:
  Int nespecies;  //number of elemental species total
  Int nspecies;   //number of species total
  Int nreactions; //number of reactions total

  Species<Type>* species;
  Reaction<Type>* reactions;
  
  //store casestring
  std::string caseString;

  //store database filename
  std::string databaseFile;
  
  //equation of state class, one for each species
  std::vector<EOS<Type>* > eos;

  //used for testing and bootstrapping
  ChemModel();
  //used for general file defined models
  ChemModel(std::string casestring, std::string databaseFile);
  ~ChemModel();

  void GetFluidProperties(const Type* rhoi, const Type T, const Type* cvi,
			  Type& cv, Type& cp, Type& R, Type& P, Type& rho, Type& c2) const;
  void GetMassProductionRates(Type* rhoi, Type T, Type* wdot); //get production rates for chemical reaction
  //finite difference implementation
  Type GetSpecificEnthalpy(Type* massfrac, Type T, Type * hi);
  Type dEtdT_dEtdRhoi_FD(Type* rhoi, Type T, Type P, Type v2, Type* dEtdRhoi);
  Type dEtdP_dEtdRhoi_FD(Type* rhoi, Type T, Type P, Type v2, Type* dEtdRhoi);
  //analytic implementation
  Type dEtdT_dEtdRhoi(Type* rhoi, Type T, Type P, Type v2, Type* dEtdRhoi);
  Type dEtdP_dEtdRhoi(Type* rhoi, Type T, Type P, Type v2, Type* dEtdRhoi);
  void dTdRhoi(Type* rhoi, Type P, Type* dTdRhoi);
  Type dRmixdRhoi(Type* rhoi, Type rho, Int i);

  Type WilkesMixtureRule(Type* rhoi, Type* property, Type T);
  Type GetViscosity(Type* rhoi, Type T);
  Type GetThermalConductivity(Type* rhoi, Type T);
  void MassFractionToMoleFraction(Type* massfrac, Type* molefrac);
  void MoleFractionToMassFraction(Type* molefrac, Type* massfrac);
  void MassToMole(Type* rho, Type* moleConc);
  
  Type GetP(const Type* rhoi, const Type T) const;
  Type GetCp(const Type* rhoi, const Type T) const;
  Type GetCv(const Type* rhoi, const Type T) const;
  
 private:
  std::vector<std::string> GetElementsInModel();
  Int ReadReactionsFile(std::string rxnfile);
  Int CheckReactionFileUnits(std::string rxnfile);
  Int ReadNumberOfReactionsFromFile(std::string rxnfile);
  Int ReadChemkinReactionsFile(std::string rxnfile);
  Int ReadNumberOfReactionsFromFileChemkin(std::string rxnfile);
  void ConvertArrheniusCGSToSI(Type& A, Type& n, Type& Ea, const Type zp_zpp);
  void BuildGlobalSpeciesList();
};

//include implementations
#include "chem.tcc"

#endif
