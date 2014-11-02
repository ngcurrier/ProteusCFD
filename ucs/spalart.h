#ifndef SPALART_H_
#define SPALART_H_

#include "general.h"
#include "kernel.h"
#include "turb.h"

//forward declaration
template <class Type> class SolutionSpace;

template <class Type>
class Spalart : public TurbulenceModel<Type>
{
public:
  
  //these are model constants
  Type sigma;
  Type cb1;
  Type cb2;
  Type kappa;
  Type cw1;
  Type cw2;
  Type cw3;
  Type cv1;
  Type ct3;
  Type ct4;
  Type cw36;
  Type cv13;

  Spalart(SolutionSpace<Type>* space);
  ~Spalart();
  
  void Initialize();
  Bool RequireWallDistance(){return true;};
  void BC_Kernel(B_KERNEL_ARGS);
  void BC_Jac_Kernel(B_KERNEL_ARGS);
  void SetTinf();
  void Source(Type nu, Type d, Type* vgrad, Type* tvars, Type vol, 
	      Type* res, Type* jac);     
  void Diffusive(Type nu, Type* tgrad, Type* tvarsL, Type* tvarsR, 
		 Type* avec, Type dgrad, Type* resL, Type* resR, 
		 Type* jacL, Type* jacR);
  Type ComputeEddyViscosity(Type rho, Type nu, Int node);
  
private:
  Spalart();
};

//include implementations
#include "spalart.tcc"

#endif

