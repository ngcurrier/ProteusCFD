#ifndef INCOMPRESSIBLE_H__
#define INCOMPRESSIBLE_H__

#include "general.h"
#include "eqnset.h"
#include <vector>

template <class Type>
class IncompressibleEqnSet : public EqnSet<Type>
{
protected:
public:
  IncompressibleEqnSet(SolutionSpace<Type>* space, Param<Type>* p);
  ~IncompressibleEqnSet();

  void InitEqnSet();
  void RoeFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta);
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
  Type ComputeViscosity(Type* Q);
  Int GetVelocityGradLocation();
  Int GetMomentumLocation();
  Int GetGradientsLocation(std::vector<Int>& gradientLoc);
  Int BadExtrapolation(Type* Q);
  void ExtrapolateVariables(Type* Qho, const Type* Q, const Type* dQedge, const Type* gradQ, 
			    const Type* dx, const Type* limiter);
  Type ComputePressure(Type* Q, Type gamma);
  void ComputeStressVector(Type* grad, Type* avec, Type mu, Type* stress);
  void GetFarfieldBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta);
  void GetInviscidWallBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta);
  void GetInternalInflowBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type pressure, 
					  Type* densities, Type* flowDirection, Type velocity);
  void GetInternalOutflowBoundaryVariables(Type* QL, Type* QR, Type pressure, Type gamma);
  void GetViscousWallBoundaryVariables(Type* QL, Type* QR, Type* vel, Int normalNode, 
				       Type* normalQ, Type Twall);
  void ModifyViscousWallJacobian(Type* QL, Type* QR, Type* vel, Int cvid, CRS<Type>* crs, Int normalNode, Type Twall);
  void ModifyViscousWallResidual(Type* res, Type* vel, Int normalNode, Type Twall);
  void ViscousJacobian(Type* QL, Type* QR, Type* dx, Type s2, Type* avec, 
		       Type mut, Type* aL, Type* aR);

 private:
  IncompressibleEqnSet();
};

//include implementations
#include "incompressible.tcc"

#endif
