#ifndef CRSMATRIX_H__
#define CRSMATRIX_H__

#include "general.h"

//forward declaration
template <class Type> class PObj;

template <class Type>
class CRSMatrix
{
private:
  //number of rows/nodes in matrix
  Int nnodes;

  //number of nodes which are of a neighbor process (parallel)
  Int gnodes;

  //parallel object
  PObj<Type>* p;

  //returns the index for a block on a particular row and column 
  Int GetIndex(Int row, Int col);
  
public:
  //matrix system data - raw access allowed by friends for efficiency reasons
  Type* M;

  //permutation vector for LU decomp
  Int* pv;

  //number of equations in each block
  Int neqn;

  //number of equations squared - size of block
  Int neqn2;

  //number of total blocks in the matrix
  Int nblocks;

  //contains the beginning and ending block numbers for each row
  //accessed like ia[0] (contains starting block number for row 0)
  //and ia[1] (contains ending block number for row 0)
  Int* ia;

  //contains column in which block on the row lies
  //accessed like ja[ia[5]+2] (contains column number for the third 
  //block in row 5 
  Int* ja;

  //contains diagonal indx for each row
  //accessed like iau[k] (contains total offset into BCRS matrix
  //for the diagonal element associated with row k)
  Int* iau;

  //flag to tell if diagonal blocks have LU decomped
  Bool ludiag;

  CRSMatrix();
  ~CRSMatrix();

  Int GetNnodes();
  Int GetGnodes();
  PObj<Type>* GetPObj();

  //Allocate memory and build internal connectivity lists
  void Init(Int nnodes_, Int gnodes_, Int neqn_, Int* ipsp, Int* psp, PObj<Type>* pobj);
  //Allocate memory for second order stencil
  void Init2(Int nnodes_, Int gnodes_, Int neqn_, Int* ipsp, Int* psp, PObj<Type>* pobj);

  //Build special preconditioning structures from another matrix
  //These functions will replace standard Init() call functionality
  void BuildBlockDiagPrecond(CRSMatrix<Type>* A);
  void CopyMatrixStructure(CRSMatrix<Type>* A);
  void BuildILU0Local(CRSMatrix<Type>* A);
  void ILU0BackSub(Type* x, Type* b);

  //returns a pointer into the CRS matrix for a block on row/column
  Type* GetPointer(Int row, Int col, Bool silent = false);

  void Blank();
  void BlankRow(Int row);
  //blanks a particular row (subrow) in the blockrow and puts a one on the diagonal
  void BlankSubRow(Int blockrow, Int subrow);
  
  //prints a row of the matrix with indices to stdout
  void PrintRow(Int row);
  //prints the whole matrix to stdout... this is going to be big...
  void PrintMatrix();

  //This is a master call to transpose the system matrix
  void CRSTranspose();
  //Testing function which will fill a matrix, transpose it, and check for correct parallel result
  void TestTranspose();

  //Gets the maximum connectivity of in all rows, sets row to correct row with max
  Int GetMaxConnectivity(Int* row);

  //Gets the ghost nodes on a row in the CRS matrix
  //returns number found
  Int GetGhostsInRow(Int row, Int* ghosts);
  
  //Get matrix memory usage
  Int GetMemUsage();

  //precompute LU decomp on diagonal blocks for SGS solve
  void PrepareSGS();

  //undo the LU decomp on the diagonal to get back the original
  //matrix
  void UndoPrepareSGS();
  
  //blanks a vector of the same relative size as the matrix
  void BlankVector(Type* v);

  //writes matrix info out in binary format
  Int WriteMatrix(std::string casename);

  //reads matrix in from binary format
  //if allocate is true will also allocate memory for the matrix class
  Int ReadMatrix(std::string casename, Bool allocate);

};

//include implementation
#include "crsmatrix.tcc"

#endif
