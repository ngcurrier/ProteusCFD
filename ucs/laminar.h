#ifndef LAMINAR_H_
#define LAMINAR_H_

#include "general.h"
#include "turb.h"

template <class Type> class SolutionSpace;

template <class Type>
class Laminar : public TurbulenceModel<Type>
{
public:
  
  Laminar(SolutionSpace<Type>* space);
  void Initialize();
  ~Laminar();
  
  void Compute();
  
private:
  Laminar();
};

//include implementations
#include "laminar.tcc"

#endif

