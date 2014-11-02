#ifndef PARALLEL_H__
#define PARALLEL_H__

#include <iostream>				
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <typeinfo>
#include <vector>

#include <mpi.h>

#include "general.h"
#include "timer.h"

//forward declarations
template <class Type> class Mesh;
template <class Type> class CRSMatrix;

Int MPI_TypeMap(int a);
Int MPI_TypeMap(float a);
Int MPI_TypeMap(double a);
Int MPI_TypeMap(std::complex<double> a);
Int MPI_TypeMap(std::complex<float> a);


template <class Type>
MPI_Datatype MPI_GetType(Type a)
{
  Int t = MPI_TypeMap(a);
  switch(t){
  case TYPE_INT:
    return MPI_INT;
    break;
  case TYPE_FLOAT:
    return MPI_FLOAT;
    break;
  case TYPE_DOUBLE:
    return MPI_DOUBLE;
    break;
  case TYPE_COMPLEX_DOUBLE:
    return MPI_DOUBLE_COMPLEX;
    break;
  case TYPE_COMPLEX_FLOAT:
    return MPI_COMPLEX;
    break;
  default:
    std::cerr << "WARNING: MPI DATATYPE NOT FOUND!" << std::endl;
    std::cerr << "WARNING: CODE IS ABOUT TO DIE!" << std::endl;
    return MPI_BYTE;
    break;
  }
};

template <class Type>
class PObj
{
 public:

  PObj(); 
  ~PObj();

  Mesh<Type>* m;

  Int GetRank();
  Int GetNp();
  
  Int ReadPartedMesh(Mesh<Type>* m, std::string casename);
  //must be called when doing bandwidth reduction operations on mesh
  Int ReorderCommNodes(Mesh<Type>* m);

  //This function build utility mappings which are required before any of the
  //following functions are valid, mesh pointer is stored from here
  Int BuildCommMaps(Mesh<Type>* m);

  //build maps which are useful when looking at nodes in a global fashion
  //globalToLocal[i] will be -1 when the node is not local
  //localToGlobal and globalToLocal need to be unallocated, allocate inside function
  //returns number of real nodes globally
  Int BuildGlobalMaps(Int** localToGlobal, Int** globalToLocal);

  Int CheckSanityCoords(Type* xyz);

  Int UpdateGeneralVectors(Type* v, Int n);
  Int UpdateXYZ(Type* xyz);

  //This will swap the appropriate matrices for a tranpose operation using a CRS matrix
  Int TransposeCommCRS(CRSMatrix<Type>* crs);

  //timer list to measure comm times
  TimerList timers;

 private:

  Bool mapsBuilt;

  Int rank;
  Int np;

  //variables which tell of extent of parallel data passes
  Int nnode;
  Int gnode;

  //CommMaps
  Int* commCountsRecv;
  Int* commCountsSend;
  Int* commOffsetsRecv;
  Int* commOffsetsSend;
  Int** nodePackingList;
  Int* nodePackingListData;

  //requests, and other MPI stuff
  MPI_Request* sRequest;
  MPI_Request* rRequest;

};

//compute the l2 norm of a each of a series of values in a long vector
template <class Type>
std::vector<Type> StridedParallelL2Norm(PObj<Type>* p, Type* data, Int nentries, Int stride)
{
  std::vector<Type> res;
  MPI_Datatype mpit;
  Int nentries_sum = 0;
  Type* l2normGlobal = new Type[stride];
  Type* l2norm = new Type[stride];

  res.resize(stride);

  for(Int j = 0; j < stride; j++){
    l2norm[j] = l2normGlobal[j] = 0.0;
  }

  for(Int i = 0; i < nentries; i++){
    for(Int j = 0; j < stride; j++){
      l2norm[j] += data[i*stride + j]*data[i*stride + j];
    }
  }
  //get mpi datatype which is appropriate
  mpit = MPI_GetType(l2norm[0]);

  MPI_Reduce(l2norm, l2normGlobal, stride, mpit, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&nentries, &nentries_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(p->GetRank() == 0){
    for(Int j = 0; j < stride; j++){
      l2normGlobal[j] = sqrt(l2normGlobal[j])/(Type)nentries_sum;
    }
  }
  MPI_Bcast(l2normGlobal, stride, mpit, 0, MPI_COMM_WORLD); 

  for(Int j = 0; j < stride; j++){
    res[j] = l2normGlobal[j];
  }

  delete [] l2normGlobal;
  delete [] l2norm;
  return res;
};

template <class Type>
Type ParallelL2Norm(PObj<Type>* p, Type* data, Int size)
{
  Int i;
  Type l2normGlobal;
  Type l2norm = 0.0;
  Int size_sum = 0;
  MPI_Datatype mpit;
  for(i = 0; i < size; i++){
    l2norm += data[i]*data[i];
  }
  //get mpi datatype which is appropriate
  mpit = MPI_GetType(l2norm);
  MPI_Reduce(&l2norm, &l2normGlobal, 1, mpit, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&size, &size_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(p->GetRank() == 0){
    l2normGlobal = sqrt(l2normGlobal)/(Type)size_sum;
  }
  MPI_Bcast(&l2normGlobal, 1, mpit, 0, MPI_COMM_WORLD); 

  return (l2normGlobal);
};

template <class Type>
Type ParallelL2NormLittle(PObj<Type>* p, Type* data, Int size)
{
  Int i;
  Type l2normGlobal;
  Type l2norm = 0.0;
  MPI_Datatype mpit;
  for(i = 0; i < size; i++){
    l2norm += data[i]*data[i];
  }
  //get mpi datatype which is appropriate
  mpit = MPI_GetType(l2norm);
  MPI_Reduce(&l2norm, &l2normGlobal, 1, mpit, MPI_SUM, 0, MPI_COMM_WORLD);
  if(p->GetRank() == 0){
    l2normGlobal = sqrt(l2normGlobal);
  }
  MPI_Bcast(&l2normGlobal, 1, mpit, 0, MPI_COMM_WORLD); 

  return (l2normGlobal);
}

template <class Type>
Type ParallelDotProduct(PObj<Type>* p, Type* v1, Type* v2, Int size)
{
  Int i;
  Type dot = 0.0;
  Type pdot = 0.0;
  MPI_Datatype mpit;
  for(i = 0; i < size; i++){
    dot += v1[i]*v2[i];
  }
  mpit = MPI_GetType(dot);
  MPI_Allreduce(&dot, &pdot, 1, mpit, MPI_SUM, MPI_COMM_WORLD);

  return(pdot);
}


//include implementations
#include "parallel.tcc"

#endif
