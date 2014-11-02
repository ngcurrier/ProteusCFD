#ifndef GAUSSIAN_H__
#define GAUSSIAN_H__

#include <cmath>
#include "general.h"
#include "solutionSpace.h"
#include "functor.h"

//provides a volume integrated source term based on a uniform 2D Gaussian
//contributes to the residual vector in solution space
template <class Type>
class GaussianSource
{
public:
  GaussianSource(Int eqnid, Int bcid, Type xloc, Type yloc, Type ampl, SolutionSpace<Type>& space);
  ~GaussianSource(){};

  void SetXY(Type xloc, Type yloc);
  void SetAmplitude(Type ampl);
  void ApplyToResidual();

  Type Evaluate(Type x, Type y, Type z);

  Type xloc;
  Type yloc;
private:
  GaussianSource();
  Type ampl;   //Amplitude
  Type sigma2; //controls diffusiveness of the source (spread)
  Int eqnid;   //controls which eqn we add the source to
  Int bcid;    //controls the bcid to which we apply the source
  SolutionSpace<Type>& space;
};

template <class Type>
Real Evaluate(Type x, Real y, Real z){return 1.0;}; 

template <class Type>
class GaussianFunctor : public SpatialFunctor<Type>
{
public:
  GaussianFunctor(GaussianSource<Type>& gaussianSource);
  Type Evaluate(Type x, Type y, Type z);
private:
  GaussianSource<Type>& gaussianSource;
};

//include implementation
#include "gaussian.tcc"

#endif
