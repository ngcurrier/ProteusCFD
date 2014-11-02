#ifndef RESIDUAL_H__
#define RESIDUAL_H__

#include "kernel.h"
#include "general.h"

template <class Type>
std::vector<Type> ComputeResiduals(SolutionSpace<Type>* space);

template <class Type>
void SpatialResidual(SolutionSpace<Type>* space);

//things which muddy the residual output to file (convergence indicator) go here
template <class Type>
void ExtraSourceResidual(SolutionSpace<Type>* space);

template <class Type>
void TemporalResidual(SolutionSpace<Type>* space, Type timestep, Int iter, Int torder);

template <class Type>
void Kernel_Inviscid_Flux(KERNEL_ARGS);

template <class Type>
void Bkernel_Inviscid_Flux(B_KERNEL_ARGS);

template <class Type>
void Kernel_Viscous_Flux(KERNEL_ARGS);

template <class Type>
void Bkernel_Viscous_Flux(B_KERNEL_ARGS);

template <class Type>
void ContributeGCL(SolutionSpace<Type>* space, Type timestep, Int iter, Int torder);

template <class Type>
void ContributeGCL2(SolutionSpace<Type>* space, Type timestep, Int iter, Int torder);

template <class Type>
void GCL_Kernel(KERNEL_ARGS);

template <class Type>
void GCL_Bkernel(B_KERNEL_ARGS);

//this is the viscous source for error transport equations
template <class Type>
void Kernel_Viscous_Src(KERNEL_ARGS);
template <class Type>
void Bkernel_Viscous_Src(B_KERNEL_ARGS);

//include implementations
#include "residual.tcc"

#endif
