#ifndef EQN_SET_H__
#define EQN_SET_H__

#include <iostream>
#include <fstream>
#include <string>

#include "general.h"
#include "eqnset_defines.h"
#include "param.h"
#include "exceptions.h"
#include "pythonInterface.h"

//forward declaration
template <class Type> class CRS;
template <class Type> class SolutionSpace;
class DataInfo;

template <class Type>
class EqnSet
{
 public:

  EqnSet();
  virtual ~EqnSet();

  Int neqn;
  Int nauxvars;
  Int varsConservative;
  SolutionSpace<Type>* space;
  Param<Type>* param;
  Type* Qinf;
  DataInfo* idata;
  DataInfo* gdata;

  //allocates memory and links eqnset to mesh object
  virtual void InitEqnSet()
  {
    Abort << "No default implementation for InitEqnSet()";
  };

  //Numerical Flux -- master calling function for flux type
  //see listing below for types
  virtual bool NumericalFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type* flux, Type beta);
  virtual bool BoundaryFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type* flux, Type beta);

  virtual Bool RoeVariables(Type* QL, Type* QR, Type gamma, Type* Qroe)
  {
    //returns false if this routine is not meaningful
    return (false);
  };
  virtual void Eigensystem(Type* Q, Type* avec, Type vdotn, Type gamma, Type* eigenvalues,
			   Type* T, Type* Tinv, Type beta)
  {
    Abort << "No default implementation for Eigensystem()";
  };
  virtual void Flux(Type* Q, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta)
  {
    Abort << "No default implementation for Flux()";
  };
  virtual void ViscousFlux(Type* Q, Type* grad, Type* avec, Type mut, Type* flux)
  {
    Abort << "No default implementation for ViscousFlux()";
  };
  virtual void ViscousJacobian(Type* QL, Type* QR, Type* dx, Type s2, Type* avec, 
			       Type mut, Type* aL, Type* aR);
  virtual Type MaxEigenvalue(Type* Q, Type* avec, Type vdotn, Type gamma, Type beta)
  {
    Abort << "No default implementation for MaxEigenvalue()";
    return (0.0);
  };
  virtual void RoeFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta)
  {
    std::cerr << "No default implementation for RoeFlux()";
  };
  
  virtual void AUSMFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta)
  {
    std::cerr << "No default implementation for AUSMFlux()";
  };
  virtual void OneSidedRoeFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta)
  {
    Abort << "No default implementation for OneSidedRoeFlux()";
  };
  virtual void HLLCFlux(Type* QL, Type* QR, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta)
  {
    Abort << "No default implementation for HLLCFlux()";
  };
  virtual void UpdateQinf()
  {
    Abort << "No default implementation for UpdateQinf()";
  };
  virtual void SetInitialConditions()
  {
    Abort << "No default implementation for SetInitialConditions()";
  };
  virtual void ApplyDQ(Type* dQ, Type* Q, Type* xyz)
  {
    //default action is to apply DQ
    //without any clipping action -- no idea what to clip
    for(Int j = 0; j < this->neqn; j++){
      Q[j] += dQ[j];
    }
  };
  virtual Type GetTheta(Type* Q, Type* avec, Type vdotn)
  {
    Abort << "GetTheta() not implemented";
    return (0.0);
  };
  virtual Type GetDensity(Type* Q)
  {
    return (0.0);
  };
  virtual Type GetPressure(Type* Q)
  {
    return (0.0);
  };
  virtual Type GetTemperature(Type* Q)
  {
    Abort << "GetTemperature() not implemented";
    return (0.0);
  };
  virtual Type GetGamma(Type* Q)
  {
    Abort << "GetGamma() not implemented";
    return (0.0);
  };
  virtual Type GetCp(Type* Q, Type gamma)
  {
    return (0.0);
  };
  virtual Type GetCf(Type tauw, Type rho)
  {
    Abort << "GetCf() not implemented";
    return (0.0);
  };
  virtual Type GetRe()
  {
    //the default implementation returns the reynolds number 
    //from the parameter object
    return(param->Re);
  };
  virtual void GetFarfieldBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta)
  {
    Abort << "GetFarfieldBoundaryVariables() not implemented";
  };
  virtual void GetInviscidWallBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta)
  {
    Abort << "GetInviscidWallBoundaryVariables() not implemented";
  };
  virtual void GetInternalInflowBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type pressure, 
						  Type* densities, Type* flowDirection, Type vel)
  {
    Abort << "GetInternalFlowInflowlBoundaryVariables() not implemented";
  };
  virtual void ModifyInternalInflowResidual(Type* res, Type* densities, Type* flowDirection, Type velocity)
  {
    Abort << "ModifyInternalInflowResidual() not implemented";
  }; 
  virtual void ModifyInternalInflowJacobian(Int cvid, CRS<Type>* crs, Type* densities, Type* flowDirection, Type velocity)
  {
    Abort << "ModifyInternalInflowJacobia() not implemented";
  };
  virtual void GetInternalOutflowBoundaryVariables(Type* avec, Type* QL, Type* QR, Type Pressure, Type Gamma)
  {
    Abort << "GetInternalFlowOutflowlBoundaryVariables() not implemented";
  };
  virtual void GetViscousWallBoundaryVariables(Type* QL, Type* QR, Type* vel, 
					       Int normalNode, Type* normalQ, Type Twall)
  {
    Abort << "GetViscousWallBoundaryVariables() not implemented";
  };
  virtual void GetTotalTempPressureBoundaryVariables(Type* QL, Type* QR, Type Pressure, Type T, Type* flowDirection)
  {
    Abort << "GetTotalTempPressureBoundaryVariables() not implemented";
  };
  virtual void GetIsothermalBoundaryVariables(Type* QL, Type* QR, Type Twall)
  {
    Abort << "GetIsothermalBoundaryVariables() not implemented";
  };
  virtual void GetPythonBoundaryVariables(Type* QL, Type* QR, Type* wallx, Type* wallAvec);
  virtual void GetHeatFluxBoundaryVariables(Type* QL, Type* QR, Type* normalQ, Type* normaldx, Type flux)
  {
    Abort << "GetHeatFluxBoundaryVariables() not implemented";
  };
  virtual void ModifyViscousWallJacobian(Type* QL, Type* QR, Type* vel, Int cvid, CRS<Type>* crs, 
					 Int normalNode, Type Twall)
  {
    Abort << "ModifyViscousWallJacobian() not implemented";
  };
  virtual void ModifyViscousWallResidual(Type* res, Type* vel, Int normalNode, Type Twall)
  {
  };
  virtual Int GetGradientsLocation(std::vector<Int>& gradientLoc)
  {
    //returns vector with the gradient locations expected by the eqnset being used
    Abort << "GetGradientsLocation() not implemented --- required for higher order solution or viscous sims";
    return -1;
  };
  virtual Int GetVelocityGradLocation()
  {
    //this is required for turbulence models
    Abort << "GetVelocityGradLocation() not implemented -- required by turbulence models";
    return -1;
  };
  virtual Int GetMomentumLocation()
  {
    Abort << "GetMomentumLocation() not implemented -- required by farfieldviscous BC";
    return -1;
  }

  virtual Int BadExtrapolation(Type* Q)
  {
    //takes extrapolated Q values and returns true if the limiter should
    //be clipped to zero b/c of some internal value being bad
    //for compressible eqnset an example would be negative pressure or energy
    return (false);
  };
  virtual void ExtrapolateVariables(Type* Qho, const Type* q, const Type* dQedge, 
				    const Type* gradQ, const Type* dx, const Type* limiter)
  {
    //FUNCTION always operates on extrapolated gradient, dQ and Q
    //         also must return variables in extrapolated format (non-conservative)

    //returns Qho with an extrapolated h.o. variable set, the eqnset is responsible
    //for knowing which gradient was requested and therefore doing the appropriate 
    //reconstruction -- paired with GetGradientsLocation()
    Abort << "ExtrapolateVariables() not implemented -- required for higher order solution";
  };
  //returns an unlimited extrapolation of Q using gradQ - location i
  Type ExtrapolateCorrection(const Type dQedge, const Type* gradQ, const Type* dx)
  {
    Type chi = this->param->chi;
    return (0.5*chi*dQedge + (1.0 - chi)*
	    (gradQ[0]*dx[0] +
	     gradQ[1]*dx[1] +
	     gradQ[2]*dx[2]));
  };
  virtual Type ComputePressure(Type* Q, Type gamma)
  {
    Abort << "ComputePressure() not implemented";
    return (0.0);
  };
  virtual Type ComputeTemperature(Type* Q, Type gamma)
  {
    Abort << "ComputeTemperature() not implemented";
    return (0.0);
  };
  virtual void ComputeAuxiliaryVariables(Type* Q)
  { };
  virtual Type ComputeViscosity(Type* Q)
  {
    Type T = GetTemperature(Q);
    Type Tref = param->ref_temperature;
    
    //this is sutherland's law
    Type mu;
    
    //Sutherland's temperature non-dimensional
    //Tref must be in Kelvin
    Type S = 110.4/Tref;
    
    //this is actually mu/mu_ref (non-dimensional)
    mu = (1.0 + S)*pow(T, 1.5)/(T + S);
  
    //return non-dimensional viscosity
    return mu;
  };
  virtual void SourceTerm(Type* Q, Type vol, Type* source)
  {
    //source term should be coded as if already on the RHS
    for(Int i = 0; i < neqn; i++){
      source[i] = 0.0;
    }
  };
  virtual void SourceTermJacobian(Type* Q, Type vol, Type* A);
  //this routine provides a contribution to jacobian of (dQ/dq)*(Vol/dt)I
  virtual void ContributeTemporalTerms(Type* Q, Type vol, Type cnp1, Type dt, 
				       Type dtau, Type* A, Type beta);
  //this is used to re-dimensionalize values before writing solution
  virtual void Dimensionalize(Type* Q)
  { };
  virtual void ComputeStressVector(Type* vgrad, Type* avec, Type mu, Type* stress)
  {
    Abort << "ComputeStressVector() not implemented";
  };

  virtual void NativeToConservative(Type* Q)
  {
    //by default, do nothing.... this must be implemented by each eqnset
    //only if the eqnset does not store the conservative variables by default
    //things like variable mach and compressible finite rate will need this
  };
  virtual void ConservativeToNative(Type* Q)
  {
    //by default, do nothing.... this must be implemented by each eqnset
    //only if the eqnset does not store the conservative variables by default
    //things like variable mach and compressible finite rate will need this
  };

  //returns the names of the residual terms (i.e. equations) for residual printing
  virtual std::vector<std::string> GetResidualNames(){
    std::vector<std::string> names;
    for(Int i = 0; i < neqn; ++i){
      std::stringstream ss;
      ss << "resQ_" << i;
      names.push_back(ss.str());
    }
    return names;
  }
  
 private:

};


//include implementations
#include "eqnset.tcc"

#endif
