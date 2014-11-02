#ifndef FORCES_H__
#define FORCES_H__

#include <iostream>
#include <fstream>

#include "kernel.h"

//forward declaration
template <class Type> class SolutionSpace;
template <class Type> class CompositeBody;

template <class Type>
class Forces
{
public:
  //this is the file pointer for lift/forces data
  std::ofstream fout;

  SolutionSpace<Type>* space;

  Forces(SolutionSpace<Type>* space);
  ~Forces();
  
  void Compute();
  void ComputeCl();
  void ComputeYpCf();

  void Report();
  void ReportCl();
  void ReportYpCf();

  //returns a list of three vectors for the forces at each node in nodelist
  //nodes in nodelist are assumed to lie on a surface, used for FSI
  void DiscreteSurfaceForces(Int* nodelist, Int npts, double* forces);

  void AllocateForcesMemory(Int num_bcs, Int num_bodies);
  void KillForcesMemory();


  //things we store only on the surface
  Type* cp;    //coefficient of pressure
  Type* yp;    //yplus
  Type* cf;    //skin friction coefficient
  Type* qdot;  //heat transfer coefficient

  Type* surfArea;

  //// Information for composite bodies, stored here so we can use it to compute
  //// Objective function effects based on CTSE
  Bool forcesAlloc;
  Int num_bcs;
  Int num_bodies;               //number of composite bodies in the bc file
  CompositeBody<Type>* bodies;     //list of the composite bodies
private:

};

template <class Type>
void ComputeSurfaceAreas(SolutionSpace<Type>* space, Int verbosity);

template <class Type>
void FORCE_Kernel(B_KERNEL_ARGS);

template <class Type>
void SURFACE_AREA_Kernel(B_KERNEL_ARGS);

template <class Type>
void YpCf_Kernel(B_KERNEL_ARGS);

//include implementations
#include "forces.tcc"

#endif
