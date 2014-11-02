#ifndef LIMITERS_H__
#define LIMITERS_H__

#include "general.h"
#include "kernel.h"

template <class Type> class SolutionSpace;
class DataInfo;

template <class Type>
class Limiter
{
public:
  Limiter(SolutionSpace<Type>* space, Type* q, Type* grad, Int neqn, Int stride, Int gradStride, std::string name);
  ~Limiter();
  
  void Compute(SolutionSpace<Type>* space);

  //array pointers for max and min data
  Type* qmin;
  Type* qmax;

  //pointer to original data for limiter
  Type* q;

  //pointer to grads for this limiter
  Type* grad;

  //pointer to limiter array itself
  Type* l;

  //number of equations to solve for
  Int neqn;
  //size of stride between data for each node
  Int nvars;
  //number of gradient entries per control volume (stride)
  Int gradneqn;

  //number of nodes
  Int nnodes;
  //number of parallel nodes
  Int gnode;

  //data info pointer
  DataInfo* idata;
  
private:
  Limiter();

};

//Venkatakrishnan limiter
template <class Type>
void Kernel_Venkat(KERNEL_ARGS);

template <class Type>
void Bkernel_Venkat(B_KERNEL_ARGS);

//Venkat modified
template <class Type>
void Kernel_VenkatMod(KERNEL_ARGS);

template <class Type>
void Bkernel_VenkatMod(B_KERNEL_ARGS);

//Barth-Jespersen limiter
template <class Type>
void Kernel_Barth(KERNEL_ARGS);

template <class Type>
void Bkernel_Barth(B_KERNEL_ARGS);

//Find local min/max
template <class Type>
void Kernel_FindMinMax(KERNEL_ARGS);

template <class Type>
void Bkernel_FindMinMax(B_KERNEL_ARGS);

template <class Type>
void Kernel_PressureClip(KERNEL_ARGS);

//inlude implementations
#include "limiters.tcc"

#endif
