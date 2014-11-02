#ifndef KERNEL_H__
#define KERNEL_H__

#include "general.h"

template <class Type> class SolutionSpace;

#define KERNEL_ARGS SolutionSpace<Type>* space, Int cvid, Int left_cv, Int right_cv, Type* __restrict__ avec, Type vdotn, Type** __restrict__ ptrL, Type** __restrict__ ptrR,  Type* __restrict__ tempL, Type* __restrict__ tempR, Int* size, void* custom, Int eid
#define B_KERNEL_ARGS SolutionSpace<Type>* space, Int cvid, Int left_cv, Int right_cv, Type* __restrict__ avec, Type vdotn, Type* __restrict__ velw, Type** __restrict__ ptrL, Type** __restrict__ ptrR, Type* __restrict__ tempL, Type* __restrict__ tempR, Int* size, void* custom, Int eid, Int factag

//encapsulate the function pointers for the kernel/driver interface
template <class Type>
class Kernel
{
private:
  //no default constructor
  Kernel() {};

  //hold on to function pointers
  void (*interior) (KERNEL_ARGS);
  void (*boundary) (B_KERNEL_ARGS);

public:
  //**************//
  //constructors  //
  //**************//
  Kernel(void(*f)(KERNEL_ARGS)){
    interior = f;
  };
  
  Kernel(void(*f)(B_KERNEL_ARGS)){
    boundary = f;
  };

  //*************//
  //operators    //
  //*************//
  void operator() (KERNEL_ARGS) const {
    interior(space, cvid, left_cv, right_cv, avec, vdotn, ptrL, ptrR, tempL, tempR, size, custom, eid);
  }
  void operator() (B_KERNEL_ARGS) const {
    boundary(space, cvid, left_cv, right_cv, avec, vdotn, velw, ptrL, ptrR, tempL, tempR, size, custom, eid, factag);
  }
  
};

#endif
