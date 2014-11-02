#ifndef CRS_H__
#define CRS_H__

#include "general.h"
#include <iostream>
#include <vector>

//forward declaration
template <class Type> class PObj;
template <class Type> class CRSMatrix;

template <class Type>
class CRS
{
private:
  //number of rows/nodes in matrix
  Int nnodes;

  //number of nodes which are of a neighbor process (parallel)
  Int gnodes;

  //parallel object
  PObj<Type>* p;

public:

  //number of equations in each block
  Int neqn;

  //number of equations squared - size of block
  Int neqn2;

  //matrix system data - raw access allowed.
  CRSMatrix<Type>* A;
  Type* x;
  Type* b;
  
  CRS();
  ~CRS();
  
  //Allocate memory and build internal connectivity lists
  void Init(Int nnodes_, Int gnodes_, Int neqn_, Int* ipsp, Int* psp, PObj<Type>* pboj);

  //will return memory usage(bytes) of the deafult system A, x, pv, and b
  Int GetSystemMemUsage(Bool print);

  //Solvers
  Type SGS(Int nSgs, CRSMatrix<Type>* A, Type* x, Type* b, Bool trueResidual = false);
  Type GMRES(Int restarts, Int nSearchDir, Int precondType, CRSMatrix<Type>* A, Type* x, Type* b);
  
  //zero out all internal data storage.  Namely, x and b.
  void BlankSystem();
  void BlankMatrix();
  void BlankMatrixRow(Int row);
  void BlankVector(Type* v);
  void BlankX();
  void BlankB();

  //Allocate new vector plus parallel update space
  Int AllocateVector(Type** v);
  //Allocate new vector as above and zero entries
  Int AllocateBlankVector(Type** v);
  
  //This uses the matrix built with the crs object
  //vin/vout must be allocated such that it can be synced in parallel
  //vin should be synced before making this call
  void MatVecMultiply(CRSMatrix<Type>* A, Type* vin, Type* vout);

  //This is used in copying a real mesh to complex one...
  //and therefore uses a slow copy... if you need a better copy implement it, don't use this
  //copies to matrix Mout from matrix Min
  template <class Type2>
  void CopyMatrixSlow(CRSMatrix<Type>* Mout, CRSMatrix<Type2>* Min);

  //This computes the preconditioner N for A
  void Preconditioner(Int type, CRSMatrix<Type>** N, CRSMatrix<Type>* A);
  //This solves the system N(x) = b for preconditiong in GMRES
  //This will return a useful solution in vector x
  void PrecondBackSolve(Int type, CRSMatrix<Type>* N, Type* x, Type* b);

};

//include implementation
#include "crs.tcc"

#endif
