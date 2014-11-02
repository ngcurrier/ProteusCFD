#ifndef FUNCTOR_H__
#define FUNCTOR_H__

#include "general.h"

//functor to hide away interesting member functions in a simple wrapper
//for interface to 3D integration routines
template <class Type>
class SpatialFunctor
{
public:
  SpatialFunctor(){};
  virtual Type Evaluate(Type x, Type y, Type z) = 0;
private:
};

#endif
