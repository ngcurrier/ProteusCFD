#ifndef MOVE_H__
#define MOVE_H__

#include <iostream>
#include "general.h"
#include "crs.h"

template <class Type> class SolutionSpace;
template <class Type> class Mesh;
template <class Type> class BoundaryConditions;

//this is master calling routine to update mesh movement during a solve
template <class Type>
void MoveMesh(SolutionSpace<Type>* space);

template <class Type, class Type2>
void UpdateMeshVelocities(Mesh<Type>* m, Type* xyz, Type* xyzold, Type* xyzoldm1, 
			  Type2 timestep, Int iter, Int torder);

template <class Type>
void UpdateMeshDisplacements(SolutionSpace<Type>* space, Type* dx, Int* reqSmoothing);

//expects nvnodes to include list of nodes which can be moved only
//no_movement_dir[0] = 0, if 1 x-dir is locked for all points
//no_movement_dir[1] = 0, if 1 y-dir is locked for all points
//no_movement_dir[2] = 0, if 1 z-dir is locked for all points
//weighting - 0 is simple, 1 is scale dependent (better)
//tol - convergence tolerance before smoothing stops
//TODO: get scaled weighting version working to tolerance
template <class Type>
void SmoothMeshLaplacian(Mesh<Type>* m, BoundaryConditions<Real>* bc, 
			 const Real tol, const Int max_iter, 
			 const bool no_movement_dir[3], const Int weighting = 0);

template <class Type>
void MoveMeshLinearElastic(Mesh<Type>* m, BoundaryConditions<Real>* bc, 
			   Type* dx, const Int iter);

template <class Type>
void SetBCLinearElastic(Mesh<Type>* m, CRS<Type>* crs, BoundaryConditions<Real>* bc);

//bumps points on a constant x, y, or z line
//dist - distance to bump points in 3 dimensions
//ptId - point is already on line
//direction - 0 (x), 1 (y), 2 (z) axis direction
template <class Type>
Int BumpNodesOnLine(const Type dist[3], const Int ptId, const Int nnodes, 
		    Type* xyz, const Int* ipsp, const Int* psp, 
		    const Int direction);

//needs bcs so we can allow symmetry planes, etc. to float
//this movement is based on updates from comp. design 
template <class Type>
void MoveMeshDesign(std::string casename, Mesh<Type>* m, BoundaryConditions<Real>* bc);

//include implementations
#include "move.tcc"

#endif
