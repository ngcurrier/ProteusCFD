#ifndef CHEM_SPECIES_H__
#define CHEM_SPECIES_H__

#include <string>
#include "elements.h"
#include "chem.h"
#include "chem_constants.h"
#include "h5layer.h"
#include "macros.h"

//NOTE: Most of the thermodynamics properties which are computed in this class
//      are sourced from NASA report RP-1311

template <class Type>
class Species
{
 public:
  
  Species();
  ~Species();

  Int nespecies;        //elemental species involved

  std::string name;     //human readable name, etc. Oxygen Ion
  std::string symbol;   //symbolic name of species O2
  Int charge;           //charge of the species
  Type MW;              //molecular weight kg/kmol
  Type R;               //(J /kg . K)
  Type href;            //enthalpy at 298K
  Type hf298;           //heat of formation 298K
  Type dhf298;          //Delta Hf(298K) -- used to shift curve to reference abs. zero

  Type thermo_coeff[5][7]; //coefficients for curve fits of Cp, H, S

  //these are the coefficients for Sutherlands law (low Temp.) from White p.29 and 32
  Bool hasViscousProps;
  Type k_coeff_White[3];  
  Type k_transition_White;
  Type mu_coeff_White[3];
  Type mu_transition_White;

  Int k_coeff_curves;
  Type k_coeff[3][6];  //coefficients McBride-Gordon NASA RP-1311 for thermal conductivity
                       // [*][0] - [*][1] low temp range - high temp range for fit
                       //          [*][2-5] A, B, C, D, parameters
  Int mu_coeff_curves;
  Type mu_coeff[3][6]; //coefficients McBride-Gordon NASA RP-1311 for viscosity
                       // [*][0] - [*][1] low temp range - high temp range for fit
                       //          [*][2-5] A, B, C, D, parameters

  //set name of species for lookup
  void Init(std::string name, Bool requireViscousProps, std::string database);
  
  Type GetCp(Type T); //specific heat constant pressure
  Type GetH(Type T, Bool shiftCurve = true);  //specific enthalpy
  Type GetG(Type T);  //Gibbs free energy
  Type GetS(Type T);  //specific entropy

  Type GetdHdT(Type T);

  Type GetViscosity(Type T);
  Type GetThermalConductivity(Type T);

  void GetThermoCoeff(Type T, Type* a);


  void Print(); 

 private:
  Int GetDBInfo(Bool requireViscousProps, std::string database);  //get info contained in database
};

//include implementations
#include "species.tcc"
#endif
