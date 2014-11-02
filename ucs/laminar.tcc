#include "solutionSpace.h"
#include "solutionField.h"
#include "mesh.h"

template <class Type>
Laminar<Type>::Laminar()
{
  //do not call this!!!
  return;
}

template <class Type>
Laminar<Type>::~Laminar()
{
  //nothing to do here, move along
  return;
}

template <class Type>
Laminar<Type>::Laminar(SolutionSpace<Type>* space) :
  TurbulenceModel<Type>(space)
{
  this->neqn = 0;

  return;
}

template <class Type>
void Laminar<Type>::Initialize()
{
  Int i;
  Mesh<Type>* m = this->space->m;
  Type* mut = this->space->GetField("mut", FIELDS::STATE_NONE);
  //set all turbulent viscosities in the mesh to zero
  for(i = 0; i < m->nnode+m->gnode; i++){
    mut[i] = 0.0;
  }

  return;
}

template <class Type>
void Laminar<Type>::Compute()
{
  //nothing to compute, this just avoids all the gyration
  //in the default implementation
  return;
}

