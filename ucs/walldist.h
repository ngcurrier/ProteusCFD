#ifndef WALL_DIST_H__
#define WALL_DIST_H__

#include "general.h"
#include "kernel.h"

template <class Type> class PObj;

//this call will allocate and update xyzParallel with the viscous points from other processes
template <class Type>
void SyncParallelPoints(PObj<Type>* p, Type* xyzLocal, Type** xyzParallel, Int nptsLocal, Int* nptsParallel);

//this is the master call for octree based wall distance calculation
template <class Type>
void ComputeWallDistOct(Type* wallDist, SolutionSpace<Type>* space);

/*********************************************************/
//
// The code that follows is for the Laplacian based method
// of solving for the wall distance.. keep separated
//
/********************************************************/

template <class Type>
Type ComputeWallDist(Type* wallDist, SolutionSpace<Type>* space);
	
template <class Type>
void Kernel_WallDist(KERNEL_ARGS);

template <class Type>
void BKernel_WallDist(B_KERNEL_ARGS);

template <class Type>
void BKernel_WallDistBC(B_KERNEL_ARGS);

template <class Type>
void Kernel_WallDist_Jac(KERNEL_ARGS);

template <class Type>
void Kernel_WallDist_Diag(KERNEL_ARGS);

template <class Type>
void BKernel_WallDist_Diag(B_KERNEL_ARGS);

template <class Type>
class PointerStruct
{
public:
  PointerStruct(Type* A_, Type* B_, Type* C_);
  ~PointerStruct();

  Type* A;
  Type* B;
  Type* C;

private:
  PointerStruct();
};

//include implementations
#include "walldist.tcc"

#endif
