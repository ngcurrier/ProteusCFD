#ifndef COMPRESSIBLE_FR_H__
#define COMPRESSIBLE_FR_H__

#include "eqnset.h"
#include "dataInfo.h"

//forward declarations
template <class Type> class Param;
template <class Type> class ChemModel;

template <class Type>
class CompressibleFREqnSet : public EqnSet<Type>
{
 public:
  CompressibleFREqnSet(SolutionSpace<Type>* space, Param<Type>* p);
  ~CompressibleFREqnSet();

  Int nspecies;  //number of reacting species
  Type Pref; //this is the freestream reference pressure P = Pref + Pgauge, 
             //needed b/c we solve everything using gauge pressure for stability
  Type PrefDim; //used to interface with EOS, etc.

  ChemModel<Type>* chem;
  Bool ownChem;

  void InitEqnSet();
  void RoeFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta);
  void HLLCFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta);
  void ViscousFlux(Type* Q, Type* grad, Type* avec, Type mut, Type* flux);
  void SetInitialConditions();
  void ApplyDQ(Type* dQ, Type* Q, Type* xyz);
  Type GetDensity(Type* Q);
  Type GetPressure(Type* Q);
  Type GetTemperature(Type* Q);
  Type GetGamma(Type* Q);
  Type GetCp(Type* Q, Type gamma);
  Int GetVelocityGradLocation();
  Int GetMomentumLocation();
  Int GetGradientsLocation(std::vector<Int>& gradientLoc);
  Int BadExtrapolation(Type* Q);
  void ExtrapolateVariables(Type* Qho, const Type* Q, const Type* dQedge, const Type* gradQ, 
			    const Type* dx, const Type* limiter);
  Type ComputePressure(Type* Q, Type gamma);
  void ComputeAuxiliaryVariables(Type* Q);
  Type GetTotalEnergy(Type* Q);
  Type GetTotalEnthalpy(Type* Q);
  void SourceTerm(Type* Q, Type vol, Type* source);
  void ContributeTemporalTerms(Type* Q, Type vol, Type cnp1, Type dt, Type dtau, Type* A, Type beta);
  void UpdateQinf();
  void GetFluidProperties(const Type* rhoi, const Type T, const Type* cvi, Type& cv,
			  Type& cp, Type& R, Type& gamma, Type& c2, Type& rho, Type& P) const;
  Type ComputeViscosity(Type* Q);
  Type GetMolecularViscosity(Type* rhoi, Type T);
  Type GetThermalConductivity(Type* rhoi, Type T);
  Type MaxEigenvalue(Type* Q, Type* avec, Type vdotn, Type gamma, Type beta);
  Type GetTheta(Type* Q, Type* avec, Type vdotn);
  void Dimensionalize(Type* Q);
  void ViscousJacobian(Type* QL, Type* QR, Type* dx, Type s2, Type* avec, 
		       Type mut, Type* aL, Type* aR);
  void GetFarfieldBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta);
  void GetInviscidWallBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta);
  void GetInternalInflowBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type pressure, 
					  Type* densities, Type* flowDirection, Type velocity); 
  void ModifyInternalInflowResidual(Type* res, Type* densities, Type* flowDirection, Type velocity); 
  void ModifyInternalInflowJacobian(Int cvid, CRS<Type>* crs, Type* densities, Type* flowDirection, Type velocity); 
  void GetInternalOutflowBoundaryVariables(Type* QL, Type* QR, Type pressure, Type gamma);
  void GetViscousWallBoundaryVariables(Type* QL, Type* QR, Type* vel, Int normalNode, Type* normalQ, Type Twall);
  void GetTotalTempPressureBoundaryVariables(Type* QL, Type* QR, Type Pressure, Type T, Type* flowDirection);
  void ModifyViscousWallResidual(Type* res, Type* vel, Int normalNode, Type Twall);
  void ModifyViscousWallJacobian(Type* QL, Type* QR, Type* vel, Int cvid, CRS<Type>* crs, 
				 Int normalNode, Type Twall);
  void NativeToConservative(Type* Q);
  void ConservativeToNative(Type* Q);
  void ComputeStressVector(Type* grad, Type* avec, Type mu, Type* stress);
  Type GetCf(Type tauw, Type rho);
  void Eigensystem(Type* Q, Type* avec, Type vdotn, Type gamma,
		   Type* eigenvalues, Type* T, Type* Tinv, Type beta);
  void GetSpeciesSpeedOfSound(Type* c2i, Type* Q);
  Type NewtonFindTGivenP(const Type* rhoi, const Type Pgoal, const Type Tinit) const;
  //returns the names of the residual terms (i.e. equations) for residual printing
  virtual std::vector<std::string> GetResidualNames(){
    std::vector<std::string> names;
    names.push_back("continuity");
    names.push_back("momentumX");
    names.push_back("momentumY");
    names.push_back("momentumZ");
    names.push_back("Energy");
    return names;
  }
 private:  
  CompressibleFREqnSet();
};

//include implementations
#include "compressibleFR.tcc"

#endif
