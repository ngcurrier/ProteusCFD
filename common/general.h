#ifndef GENERAL_H__
#define GENERAL_H__

#include <complex>

typedef double Real;
//typedef float Real;
typedef int Int;
typedef unsigned int UInt;
typedef bool Bool;
typedef std::complex<Real> RCmplx;

struct IntInt
{
  int a, b;
};

//this has to be here b/c as of this date,
//local type instatiation cannot be resolved... 
template <class Type>
struct TypeInt{Type dist; Int rank; Int nodeId;};

#define NEQN 5
#define NUM_THREADS 2

#define TYPE_INT 0
#define TYPE_FLOAT 1
#define TYPE_DOUBLE 2
#define TYPE_COMPLEX_DOUBLE 3
#define TYPE_COMPLEX_FLOAT 4

//this is utilized in the MPI_MINLOC and MPI_MAXLOC operations
struct RealInt
{
  Real val;
  Int rank;
};

#endif
