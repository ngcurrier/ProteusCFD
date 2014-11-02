#ifndef JACOBIAN_H__
#define JACOBIAN_H__

#include "general.h"
#include "kernel.h"

//forward declarations
template <class Type>
void ComputeJacobians(SolutionSpace<Type>* space);

//type = 0 - upwind FD
//type = 1 - central FD
//type = 3 - complex 
//type = 4 - mixed for design
template <class Type>
void Compute_dRdQ(CRS<Type>& crs, SolutionSpace<Type>* space, Int type);

//type here is the same as Compute_dRdQ as it is simply passed through
template <class Type>
void Compute_dRdQ_Transpose(CRS<Type>& crs, SolutionSpace<Type>* space, Int type);

template <class Type>
void ComputeSpatialJacobian(SolutionSpace<Type>* space);

template <class Type>
void ContributeTemporalTerms(SolutionSpace<Type>* space);

template <class Type>
void Kernel_NumJac(KERNEL_ARGS);

template <class Type>
void Kernel_NumJac_Centered(KERNEL_ARGS);

template <class Type>
void Kernel_NumJac_Complex(KERNEL_ARGS);

template <class Type>
void Kernel_Diag_NumJac(KERNEL_ARGS);

template <class Type>
void Bkernel_NumJac(B_KERNEL_ARGS);

template <class Type>
void Bkernel_NumJac_Centered(B_KERNEL_ARGS);

template <class Type>
void Bkernel_NumJac_Complex(B_KERNEL_ARGS);

template <class Type>
void Kernel_Viscous_Jac(KERNEL_ARGS);

template <class Type>
void Bkernel_Viscous_Jac(B_KERNEL_ARGS);

#if 0
template <class Type>
void Diag_Jacobian_Check(Real* JacCheck, SolutionSpace<Type>* space, Int ptid);
#endif

//include implementations
#include "jacobian.tcc"

#endif
