#ifndef TIMESTEP_H__
#define TIMESTEP_H__

#include <iostream>
#include <cmath>
#include "general.h"
#include "kernel.h"

//forward declarations
template <class Type> class EqnSet;

template <class Type>
Type ComputeTimesteps(EqnSet<Type>* eqnset);

template <class Type>
void Kernel_Timestep(KERNEL_ARGS);

template <class Type>
void Bkernel_Timestep(B_KERNEL_ARGS);

//include implementations
#include "timestep.tcc"

#endif
