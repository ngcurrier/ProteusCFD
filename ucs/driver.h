#ifndef DRIVER_H__
#define DRIVER_H__

#include "uns_base.h"
#include "general.h"
#include "kernel.h"
#include "threaded.h"
#include "geometry.h"
#include <iostream>
#include <cmath>
#include <vector>


//forward declarations
template <class Type> class EqnSet;
template <class Type> class Mesh;
template <class Type> class SolutionSpace;


//function which operates on a volume of a solution space
//provides scattering to the nodes which make up a volume face/edge
template <class Type>
Int Driver(SolutionSpace<Type>* space, Kernel<Type> &kernel, Int scatterSize, 
	   void* custom);

//function which operates on a surface of a solution space
//provides scattering to the nodes which make up the boundary edges
template <class Type>
Int Bdriver(SolutionSpace<Type>* space, Kernel<Type> &bkernel, Int scatterSize, 
	    void* custom);


//function which operates on a volume of a solution space
//does not provide scattering to the nodes which make up a volume face/edge
template <class Type>
Int DriverNoScatter(SolutionSpace<Type>* space, Kernel<Type> &kernel, 
		    Int scatterSize, void* custom);

//function which operates on a surface of a solution space
//does not provide scattering to the nodes which make up the boundary edges
template <class Type>
Int BdriverNoScatter(SolutionSpace<Type>* space, Kernel<Type> &bkernel, 
		     Int scatterSize, void* custom);

//function which performs scattering to the nodes of an internal/external edge
//controls mutex locking, etc. for contentious resources in memory
template <class Type>
void DriverScatter(Int left_cv, Int right_cv, Type* ptrL, Type* ptrR, 
		   Type* valL, Type* valR, Int n, std::vector<pthread_mutex_t*>& mutexes);

//include implementations
#include "driver.tcc"

#endif
