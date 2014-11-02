#ifndef ENDIAN_UTIL_H__
#define ENDIAN_UTIL_H__

#include <iostream>
#include "general.h"

double ReverseEndian_Double8(double in);
int ReverseEndian_Int4(int i);
short int ReverseEndian_Int2(short int i);
float ReverseEndian_Float4(float in);
void ReverseEndian_Double8Array(double *in,int n);
void ReverseEndian_Int4Array(int * in,int n);
void ReverseEndian_Float4Array(float*in,int n);
void ReverseEndian_Int2Array(short int* in, int n);

int IsLittleEndian(void);

template<class T>
T swap_endian(T u)
{
  union{
    T u;
    unsigned char u8[sizeof(T)];
  } source, dest;

  source.u = u;

  for(size_t k = 0; k < sizeof(T); k++){
    dest.u8[k] = source.u8[sizeof(T) - k - 1];
  }
  return dest.u;
}

template<class theType>
theType ReverseEndianReal(theType in)
{
  Int bytes = sizeof(theType);
  if(bytes == 8){
    return (ReverseEndian_Double8((double)in));
  }
  else if(bytes == 4){
    return (ReverseEndian_Float4((float)in));
  }
  else{
    std::cout << "ReverseEndianReal: number of bytes " << bytes 
	      << " is strange.. what to do?" << std::endl;
    return(-1.0);
  }
}

template<class theType>
theType ReverseEndianInt(theType in)
{
  Int bytes = sizeof(theType);
  if(bytes == 4){
    return (ReverseEndian_Int4((int)in));
  }
  else if(bytes == 2){
    return (ReverseEndian_Int2((short int)in));
  }
  else{
    std::cout << "ReverseEndianInt: number of bytes " << bytes 
	      << " is strange.. what to do?" << std::endl;
    return(-1.0);
  }
}


template<class theType>
void ReverseEndianArrayReal(theType* in, int n)
{
  Int bytes = sizeof(theType);
  if(bytes == 8){
    ReverseEndian_Double8Array((double*)in, n);
  }
  else if(bytes == 4){
    ReverseEndian_Float4Array((float*)in, n);
  }
  else{
    std::cout << "ReverseEndianArrayReal: number of bytes " << bytes 
	      << " is strange.. what to do?" << std::endl;
  }
  return;
}

template<class theType>
void ReverseEndianArrayInt(theType* in, int n)
{
  Int bytes = sizeof(theType);
  if(bytes == 4){
    ReverseEndian_Int4Array((int*)in, n);
  }
  else if(bytes == 2){
    ReverseEndian_Int2Array((short int*)in, n);
  }
  else{
    std::cout << "ReverseEndianArrayInt: number of bytes " << bytes 
	      << " is strange.. what to do?" << std::endl;
  }
  return;
}

#endif
