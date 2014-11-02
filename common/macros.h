#ifndef MACROS_H__
#define MACROS_H__

#include "general.h"
#include <complex>
#include <cmath>
#include <iostream>

#define CHECKPT {std::cout << "Checkpoint: " << __FILE__ << ", line " << __LINE__ << std::endl;std::cerr << "Checkpoint: " << __FILE__ << ", line " << __LINE__ << std::endl;}

//used so we can overload imag() and real() for non-complex types
inline Real real(Real x)
{
  return x;
};

inline Real imag(Real x)
{
  return 0.0;
};

//this function returns true if the value
template <class Type>
inline Bool isWholeNumber(Type x)
{
  Int up = (Int)(real(x) + 0.5);
  Int flat = (Int)real(x);
  if(up == flat) return true;
  return false;
};

template <class theType, class theType2>
inline theType2 MAX(theType x, theType2 y)
{
  return ((real(x) > real(y)) ? x : y);
};

template <class theType, class theType2>
inline theType2 MIN(theType x, theType2 y)
{
  return ((real(x) < real(y)) ? x : y);
};

template <class theType>
inline theType Abs(theType x)
{
  return fabs(x);
};

template <>
inline RCmplx Abs(RCmplx x){
  return abs(x);
  //return (RCmplx) sqrt(real(x)*real(x) + imag(x)*imag(x));
};

template <class theType>
inline theType CAbs(theType x)
{
  return fabs(x);
};

template <>
inline RCmplx CAbs(RCmplx x)
{
  RCmplx b = x;
  if(real(x) < 0.0) b = -x;
  return b;
};

template <class Type>
Bool AlmostEqualRelative(Type A, Type B, Type maxRelDiff)
{
    // Calculate the difference.
    Type diff = Abs(A - B);
    A = Abs(A);
    B = Abs(B);
    // Find the largest
    Type largest = (B > A) ? B : A;
 
    if (diff <= largest * maxRelDiff){
      return true;
    }
    return false;
}

#endif
