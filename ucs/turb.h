#ifndef TURB_H__
#define TURB_H__

#include <iostream>
#include "general.h"
#include "kernel.h"
#include "dataInfo.h"
#include "crs.h"
#include "param.h"
#include "limiters.h"

//forward declaration
template <class Type> class SolutionSpace;
template <class Type> class SolutionField;

template <class Type>
class TurbulenceModel
{
public:
  Int neqn;                    //number of equations in turbulence model

  SolutionSpace<Type>* space;

  Type* tvar;              //turbulent variables
  Type* tvarold;
  Type* tvaroldm1;
  Type* tvarinf;           //farfield values for turbulent variables
  Type* tgrad;             //gradient of turbulent variables
  
  CRS<Type> crs;           //compressed row storage object.. stores matrix system
  DataInfo* idata;         //datainfo object for turbulent system
  Limiter<Type>* limiter;  //limiter for higher order turbulence model

  TurbulenceModel(SolutionSpace<Type>* space);      //initialize variables
  virtual ~TurbulenceModel();
  virtual Type Compute();                         //runs the solver for turb vars
  virtual Bool RequireWallDistance(){return false;};

  //sets infinity values for turbulence variables
  virtual void SetTinf();
  //sets initial conditions for the field
  virtual void Initialize();
  //update turbulent boundary conditions
  void UpdateBCs();             
  //kernel function with appropriate arguments for bdriver
  virtual void BC_Kernel(B_KERNEL_ARGS);
  //kernel function which modifies CRS system for solving
  virtual void BC_Jac_Kernel(B_KERNEL_ARGS);
  //compute source terms for the model 
  virtual void Source(Type nu, Type d, Type* vgrad, Type* tvars, Type vol, 
		      Type* res, Type* jac);     
  //calls the driver routines and then the appropriate member functions for diffusion
  void DiffusiveDriver();
  //compute diffusive terms for the model
  virtual void Diffusive(Type nu, Type* tgrad, Type* tvarsL, Type* tvarsR, 
			 Type* avec, Type dgrad, Type* resL, Type* resR, 
			 Type* jacL, Type* jacR);
  //compute convective terms for the model  
  virtual void Convective();                                          
  //compute mut... the point of this excercise.. 
  virtual Type ComputeEddyViscosity(Type rho, Type nu, Int node);
  virtual void TemporalResidual(Type timestep, Int iter, Int torder);
  //contributes temporal terms to jacobians
  virtual void ContributeTemporalTerms(Type vol, Type cnp1, Type dt, Type dtau, Type* A);

private:
};

template <class Type>
void Kernel_Convective(KERNEL_ARGS);

template <class Type>
void Bkernel_Convective(B_KERNEL_ARGS);

template <class Type>
void Kernel_Diffusive(KERNEL_ARGS);

template <class Type>
void Bkernel_Diffusive(B_KERNEL_ARGS);

template <class Type>
void BC_Turb_Kernel(B_KERNEL_ARGS);

template <class Type>
void BC_Turb_Jac_Kernel(B_KERNEL_ARGS);


#include "spalart.h"
#include "laminar.h"

template <class Type>
Int CreateTurbModel(TurbulenceModel<Type>** turb, SolutionSpace<Type>* space)
{
  Param<Type>* param = space->param;
  if(param->turbModel == 0){
    *turb = new Laminar<Type>(space);
  }
  else if(param->turbModel == 1){
    *turb = new Spalart<Type>(space);
  }
  else{
    std::stringstream ss;
    ss << "Turbulence model type " << space->param->turbModel << " not found.. FATAL!" << std::endl;
    Abort << ss.str();
    return (-1);
  }
  return (0);
};



//include implementations
#include "turb.tcc"

#endif
