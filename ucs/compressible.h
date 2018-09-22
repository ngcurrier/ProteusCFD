#ifndef COMPRESSIBLE_H__
#define COMPRESSIBLE_H__

#include "general.h"
#include "eqnset.h"

template <class Type>
class CompressibleEqnSet : public EqnSet<Type>
{
 public:
  CompressibleEqnSet(SolutionSpace<Type>* space, Param<Type>* p);
  ~CompressibleEqnSet();

  void InitEqnSet();
  void RoeFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta);
  void OneSidedRoeFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta);
  Bool RoeVariables(Type* QL, Type* QR, Type gamma, Type* Qroe);
  void Eigensystem(Type* Q, Type* avec, Type vdotn, Type gamma, Type* eigenvalues,
		   Type* T, Type* Tinv, Type beta);
  void Flux(Type* Q, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta);
  void ViscousFlux(Type* Q, Type* grad, Type* avec, Type mut, Type* flux);
  Type MaxEigenvalue(Type* Q, Type* avec, Type vdotn, Type gamma, Type beta);
  void UpdateQinf();
  void SetInitialConditions();
  void ApplyDQ(Type* dQ, Type* Q, Type* xyz);
  Type GetTheta(Type* Q, Type* avec, Type vdotn);
  Type GetDensity(Type* Q);
  Type GetPressure(Type* Q);
  Type GetTemperature(Type* Q);
  Type GetGamma(Type* Q);
  Type GetCp(Type* Q, Type gamma);
  Type GetCf(Type tauw, Type rho);
  Type GetRe();
  void SourceTerm(Type* Q, Type vol, Type* source);
  Int GetVelocityGradLocation();
  Int GetMomentumLocation();
  Int GetGradientsLocation(std::vector<Int>& gradientLoc);
  void ExtrapolateVariables(Type* Qho, const Type* Q, const Type* dQedge, 
			    const Type* gradQ, const Type* dx, const Type* limiter);
  Int BadExtrapolation(Type* Q);
  Type ComputePressure(Type* Q, Type gamma);
  Type ComputeTemperature(Type* Q, Type gamma);
  void ComputeStressVector(Type* vgrad, Type* avec, Type mu, Type* stress);
  void ComputeAuxiliaryVariables(Type* Q);
  void GetFarfieldBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta);
  void GetInviscidWallBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta);
  void GetInternalInflowBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type pressure, 
					  Type* densities, Type* flowDirection, Type velocity);
  void GetInternalOutflowBoundaryVariables(Type* QL, Type* QR, Type pressure, Type gamma);
  void GetViscousWallBoundaryVariables(Type* QL, Type* QR, Type* vel, Int normalNode, 
				       Type* normalQ, Type Twall);
  void ModifyViscousWallJacobian(Type* QL, Type* QR, Type* vel, Int cvid, CRS<Type>* crs, 
				 Int normalNode, Type Twall);
  void ModifyViscousWallResidual(Type* res, Type* vel, Int normalNode, Type Twall);
  void ViscousJacobian(Type* QL, Type* QR, Type* dx, Type s2, Type* avec, 
  		       Type mut, Type* aL, Type* aR);
 private:
  CompressibleEqnSet();
};

//include implementations
#include "compressible.tcc"

#endif
