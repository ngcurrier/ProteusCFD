#ifndef DERIVATIVES_H__
#define DERIVATIVES_H__

#include <fstream>
#include <string>
#include "general.h"

//not sure why but not including this causes the compiler
//to complain about the > operator all over the place... 
#include "solve.h"

//forward declarations
template <class Type> class Mesh;
template <class Type> class SolutionSpace;
template <class Type> class SolutionOperations;

template <class Type>
Type Compute_Obj_Function(SolutionSpace<Type>& space);

void Compute_dRdBeta_CTSE(Real* dRdB, SolutionSpace<Real>& space, Int beta);

void Compute_dRdBeta_FD(Real* dRdB, SolutionSpace<Real>& space, Int beta);

//WARNING: Calling these routines will completely nuke any dQ, Jacobian, etc. values stored in 
//Mesh object since it uses underlying SGS solver to iterate on dQdBeta
void Compute_dQdBeta(Real* dQdB, SolutionSpace<Real>& space, Int beta);

//compute dQ/dB using a complex perturbation
void Compute_dQdBeta_CTSE(Real* dQdB, SolutionSpace<Real>& space, Int beta);

//computes dI/dB by perturbing beta in the complex part
//for beta = 0, 1, 2, ...
Real Compute_dObjdBeta_CTSE(SolutionSpace<Real>& space, SolutionOrdering<Real>& operations, Int beta);

//computes dI/dB = dI/dQ*dQ/dB + dI/dx*dx/dB + dI/dB
Real Compute_dObjdBeta_Direct(SolutionSpace<Real>& space, Int beta);

//computes product of dI/dQ * dQ/dB - direct diff method
Real Compute_dObjdQ_dQdBeta(Real* dQdB, SolutionSpace<Real>& space);

//computes product of dI/dx * dx/dB - direct diff method
Real Compute_dObjdX_dXdBeta(SolutionSpace<Real>& space, Int beta);

//computes dI/dB where beta is a manipulation of a flow field parameter (mach, AOA, etc.)
void Compute_dObjdBeta_dParam();

//computes dQ/dB with a central difference .. grid must converged on passing
void Compute_dQdBeta_FD(Real* dQdB, SolutionSpace<Real>& space, Int beta);

//read in dX/dB for grid
void Get_dXdB(std::string path, Real* dxdb, Mesh<Real>* m, Int beta);

//compute dI/dQ for adjoint
void Compute_dObjdQ(Real* dIdQ, SolutionSpace<Real>& space);

//compute adjoint variables
void Compute_Adjoint_Variables(Real* lambda, SolutionSpace<Real>& space);

//compute adjoint variables using incremental iterative approach
void Compute_Adjoint_Variables_II(Real* lambda, SolutionSpace<Real>& space);

//compute dI/dB using adjoint variables... lambdas must be passed in
Real Compute_dObjdBeta_Adjoint(Real* lambda, SolutionSpace<Real>& space, Int beta);

//compute dX/dB (x,y,z) and write mesh sensitivities to file
void Compute_dXdB(SolutionSpace<Real>& space);

//matrix free computation of dRdQ*(value vector), vector/prod are neqn*nnodes long
void Compute_dRdQ_Product_MatrixFree(SolutionSpace<RCmplx>& cspace, Real* vector, Real* prod);

//compute dQ/dB with inremental iterative method
void Compute_dQdBeta_II(Real* dQdB, SolutionSpace<Real>& space, Int beta);

//perturb appropriate value in parameter object
void Perturb_Param(Int beta, Param<RCmplx>& param);

//perturb appropriate value in parameter object (sgn should be +1 or -1) to affect direction
void Perturb_Param(Int beta, Param<Real>& param, Int sgn);

//computes dqdbeta for bcs which are set via hard bc conditions
void ComputedQdBeta_HardsetBC(Int nnodes, Int* nodes, Real* dQdB, const SolutionSpace<Real>& space, Int beta);

//include implementations
#include "derivatives.tcc"

#endif
