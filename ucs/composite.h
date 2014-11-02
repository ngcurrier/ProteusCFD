#ifndef COMPOSITE_H__
#define COMPOSITE_H__

#include <iostream>
#include <string>
#include "strings_util.h"

template <class Type>
class CompositeBody
{
public:
  
  Int nsurfs;
  Int* list;

  Type* momentPt;
  Type* momentAxis;
  Type* surfArea;
  Type* forces;
  Type* vforces;
  Type* moments;
  Type* vmoments;

  Type cl, cm, cd;

  std::string name;

  CompositeBody();
  ~CompositeBody();
  

  void Init(Int* list, Int nsurfs);
  void Print(Int id = -1);
  void CountBodiesInFile(std::ifstream& fin);

  //return true if surfId is part of the composite body
  Bool SurfIsPart(Int surfId);
private:
};

template <class Type>
Int BodiesReadFile(CompositeBody<Type>* bodies, std::string casename);

template <class Type>
Int ParseBody(CompositeBody<Type>* bodies, std::string& line);

//include implementations
#include "composite.tcc"

#endif
