#ifndef GRADIENT_H__
#define GRADIENT_H__

#include "general.h"
#include "kernel.h"

//forward declaration
template <class Type> class SolutionSpace;

template <class Type>
class Gradient{

public:
  Gradient(Int nterms, Int stride, Int* list, Type* data, 
	   SolutionSpace<Type>* space, Int gradType=0, Type* grad=NULL, Bool weighted = true);

  ~Gradient();

  void Compute();
  Int GetNterms(){return nterms;};
  Int GetStride(){return stride;};
  Bool IsWeighted(){return weighted;};
  
  Type* data;  //pointer to things to take gradient of
  Type* grad;  //computed gradients here 
  Int* list;   //zero indexed id list of values which we need the gradient of
               //i.e. 0, 2, 4 would pick the first value, 3rd value, and 5 value
               //the gradient location however will be 0,1,2 when stored
 
 
private:
  
  Int nterms;    //number of values which should be in gradient vector
  Int nnode;    //number of nodes
  Int stride;   //distance between start of each set of Q values
  Int type;     //0 - LSQ
                //1 - Green Gauss
  Bool weighted; //boolean flag to use weighted or non-weighted LSQ-grad
  
  SolutionSpace<Type>* space;


  Bool allocated; //check to see if we should free the memory

};

//note: LSQ method comes straight from Daniel Hyams dissertation

template <class Type>
void ComputeNodeLSQCoefficients(SolutionSpace<Type>* space);

template <class Type>
void Kernel_Green_Gauss_Gradient(KERNEL_ARGS);

template <class Type>
void Bkernel_Green_Gauss_Gradient(B_KERNEL_ARGS);

template <class Type>
void Kernel_LSQ_Gradient(KERNEL_ARGS);

template <class Type>
void Bkernel_LSQ_Gradient(B_KERNEL_ARGS);

template <class Type>
void Kernel_LSQ_Coefficients(KERNEL_ARGS);

template <class Type>
void Bkernel_LSQ_CoefficientsWeighted(B_KERNEL_ARGS);

template <class Type>
void Kernel_LSQ_CoefficientsWeighted(KERNEL_ARGS);

template <class Type>
void Bkernel_LSQ_Coefficients(B_KERNEL_ARGS);

template <class Type>
void Bkernel_Symmetry_Fix(B_KERNEL_ARGS);

//this computes the gradient simply with the directional derivative method
//used for numerical derivatives, etc.
//list is an integer list with n values to define which parts of qL/qR
//we need the gradient of
//returns the approximate derivative of the computed gradient and the grad itself
template <class Type>
Type DirectionalDerivativeGrad(Type* dx, Type s2, Type* avec, Type* qL, Type* qR, Type* grad,
			       std::vector<Int> gradLoc);

//include implementations
#include "gradient.tcc"

#endif
