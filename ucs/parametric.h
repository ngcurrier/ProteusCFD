#ifndef PARAMETRIC_H__
#define PARAMETRIC_H__

#include "general.h"

template <class Type> class Mesh;
template <class Type> class BoundaryConditions;

//a - maximum bump magnitude
//t1 - controls maximum point of bump
//t2 - controls width of bump (smaller number is wider)
//for a smooth function... t1 must be > 1.0
//n - number of points
//y - array holding output displacement
//x - MUST BE normalized x coordinate of points n 
//valid from 0 <= x <= 1
template <class Type>
void Hicks_Henne(Type a, Type t1, Type t2, Type* x, Type* y, Int n);

//derivative w.r.t. a
template <class Type>
void Da_Hicks_Henne(Type a, Type t1, Type t2, Type* x, Type* y, Int n);

//derivative w.r.t. t1
template <class Type>
void Dt1_Hicks_Henne(Type a, Type t1, Type t2, Type* x, Type* y, Int n);

//derivative w.r.t. t2
template <class Type>
void Dt2_Hicks_Henne(Type a, Type t1, Type t2, Type* x, Type* y, Int n);


//normalize x coordinates from 0 to 1
//if maxVal == 0 returns max point for de-normalization
template <class Type>
Type NormalizeX(Type* x, Int n, Type maxVal = 0.0);


//computes the movement of the surface for the given beta values
template <class Type>
void Compute_dX_Surface(std::string casename, Mesh<Type>* m, BoundaryConditions<Real>* bc, 
			Type* dxyz, Int beta);

//compute dX/dB (x,y,z) for surface points for the given beta values
template <class Type>
void Compute_dXdB_Surface(std::string casename, Mesh<Type>* m, BoundaryConditions<Real>* bc, 
			  Type* dxdbSurf, Int beta);

//include implementations
#include "parametric.tcc"

#endif
